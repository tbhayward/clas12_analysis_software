#!/usr/bin/env python3
"""Particle-level additional-resolution and momentum-bias scan using cached MC.

Data and component scales come from a completed background-aware diagnostic.
All ROOT I/O happens once while building the uniform pre-cut reconstructed cache.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import copy
import csv
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import threading
import tempfile
import numpy as np
import study_exclusivity_resolution as diag

SOURCE=Path(__file__).resolve().parent
VERSION='particle-response-v1'
REGIONS=['electron_FD','proton_FD','proton_CD','photon_FD','photon_FT']
OBS=['Delta_phi','theta','theta_gamma_gamma','pTmiss','Emiss2','Mx2','Mx2_2']
TRAIN_OBS=[0,4,5]
VALID_OBS=[2,3,6]
TOPS=['FD_FD','CD_FD','CD_FT']
PARAMS=[f'{r}_{c}' for r in REGIONS for c in ('p_sigma_fraction','theta_sigma_rad','phi_sigma_rad','p_bias_fraction')]

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value):Path(path).write_text(json.dumps(diag.clean_json(value),indent=2,allow_nan=False)+'\n')
def table(path,rows):
    if not rows:return
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def parameter_grid(topologies):
    regions={0}
    for top in topologies:
        regions.add(1 if top=='FD_FD' else 2);regions.add(4 if top=='CD_FT' else 3)
    result={}
    for r in sorted(regions):
        widths=[0,.002,.006,.012] if r==0 else ([0,.005,.015,.03] if r==1 else [0,.01,.03,.06])
        for coordinate in range(4):
            result[4*r+coordinate]=widths if coordinate==0 else [0,.001,.003,.006] if coordinate==1 else [0,.005,.015,.03] if coordinate==2 else [-.006,-.002,0,.002,.006]
    return result

def compile_worker(output,reference):
    build=output/'backend';build.mkdir()
    for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):
        source=reference/'backend'/n
        if not source.is_file():raise FileNotFoundError('Missing saved global-cut source: '+str(source))
        shutil.copy2(source,build/n)
    for n in ('exclusivity_resolution_worker.cpp','exclusivity_resolution_core.h','particle_response_core.h','particle_response_worker.cpp'):shutil.copy2(SOURCE/n,build/n)
    included=build/'exclusivity_resolution_worker.cpp'
    included.write_text(included.read_text().split('\nint main(',1)[0]+'\n')
    flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
    worker=build/'particle_response_worker'
    subprocess.run([os.environ.get('CXX','c++'),*flags,'-O3','-std=c++17','-I'+str(build),str(build/'particle_response_worker.cpp'),str(build/'global_cuts.cpp'),'-o',str(worker),*libs],check=True)
    return worker

class Budget:
    def __init__(self,n):self.sem=threading.Semaphore(n);self.lock=threading.Lock();self.cancel=threading.Event();self.active=set()
    def start(self,cmd,stdout):
        self.sem.acquire()
        try:
            with self.lock:
                if self.cancel.is_set():raise RuntimeError('Cancelled after another period failed')
                p=subprocess.Popen(cmd,stdout=stdout,stderr=subprocess.STDOUT);self.active.add(p);return p
        except BaseException:self.sem.release();raise
    def finish(self,p):
        with self.lock:self.active.discard(p)
        self.sem.release()
    def cancel_all(self):
        self.cancel.set()
        with self.lock:
            for p in self.active:p.terminate()


class Pool:
    def __init__(self,worker,spec,output,workers,budget=None):self.worker=worker;self.spec=spec;self.output=output;self.workers=workers;self.budget=budget or Budget(workers)
    def run(self,mode,jobs,label):
        folder=self.output/label;folder.mkdir();active=set();lock=threading.Lock();cancel=threading.Event()
        def launch(item):
            i,job=item;jp=folder/f'{i}.json';target=folder/f'{i}.bin';log=folder/f'{i}.log';dump(jp,job)
            with log.open('w') as f:
                with lock:
                    if cancel.is_set():raise RuntimeError('Cancelled after worker failure')
                    p=self.budget.start([str(self.worker),mode,str(self.spec),str(jp),str(target)],f);active.add(p)
                code=p.wait();self.budget.finish(p)
                with lock:active.discard(p)
            if code:raise RuntimeError(f'{label}: worker failed; see {log}\n'+log.read_text()[-2000:])
            return i,job,target
        pool=ThreadPoolExecutor(max_workers=self.workers);futures=[];results=[]
        try:
            futures=[pool.submit(launch,item) for item in enumerate(jobs)]
            for future in as_completed(futures):
                record=future.result();results.append(record);job=record[1];print(f'[{self.output.name}/{label}] {len(results)}/{len(jobs)} complete: {job.get("period","")}/{job.get("sample","")}',flush=True)
        except BaseException:
            cancel.set();self.budget.cancel_all()
            with lock:
                for p in active:p.terminate()
            for f in futures:f.cancel()
            raise
        finally:pool.shutdown(wait=True,cancel_futures=True)
        return sorted(results)

def build_cache(pool,spec,manifest,reference,reuse=None):
    if reuse:
        saved=json.loads((reuse/'cache_manifest.json').read_text())
        if saved['version']!=VERSION or saved['reference_sha256']!=digest(reference/'manifest.json') or saved['seed']!=spec['seed']:raise ValueError('Cached events have incompatible provenance/seed')
        result={}
        for key,r in saved['samples'].items():
            path=reuse/r['file']
            if digest(path)!=r['sha256']:raise ValueError('Corrupt event cache: '+key)
            result[key]=dict(r,path=str(path.resolve()))
        print('[cache] Verified reusable MC cache; no ROOT reads',flush=True)
        return result
    jobs=[];inventory={}
    for p in spec['periods']:
        for sample in ('mc','background'):
            records=[r for r in manifest['inventory'] if r['period']==p and r['sample']==sample]
            if not records:raise ValueError('Missing reference inventory for '+p+'/'+sample)
            files=[]
            for record in records:
                f=Path(record['path']);st=f.stat()
                if st.st_size!=record['size'] or st.st_mtime_ns!=record['mtime_ns']:raise ValueError('ROOT input changed since reference study: '+str(f))
                files.append(str(f));inventory[str(f)]=record
            jobs.append(dict(period=p,sample=sample,files=files))
    result={}
    for i,job,path in pool.run('cache',jobs,'cache'):
        meta=json.loads(Path(str(path)+'.json').read_text())
        if not meta['records'] or path.stat().st_size!=meta['records']*meta['record_bytes']:raise ValueError('Incomplete event cache')
        if meta['record_bytes']!=200:raise ValueError('Unexpected native record layout')
        for f in job['files']:
            st=Path(f).stat();old=inventory[f]
            if (st.st_size,st.st_mtime_ns)!=(old['size'],old['mtime_ns']):raise ValueError('ROOT file changed during cache scan')
        result[job['period']+'__'+job['sample']]=dict(path=str(path.resolve()),file=str(path.resolve()),sha256=digest(path),**meta)
    dump(pool.output/'cache_manifest.json',dict(version=VERSION,reference_sha256=digest(reference/'manifest.json'),seed=spec['seed'],samples={k:{n:v for n,v in r.items() if n!='path'} for k,r in result.items()}))
    return result

class Engine:
    def __init__(self,pool,spec,cache,data,norm,fit_topologies,min_mc):
        self.pool=pool;self.spec=spec;self.cache=cache;self.data=data;self.norm=norm
        self.fit_topologies=[TOPS.index(t) for t in fit_topologies];self.min_mc=min_mc;self.calls=0;self.history=[]
    def evaluate(self,parameters,label,detailed=False):
        # Worker slot zero is always the exact nominal. Parameter rows are deduplicated.
        unique=[tuple([0.]*20)]
        for p in parameters:
            if tuple(p) not in unique:unique.append(tuple(p))
        jobs=[dict(period=p,sample=s,cache=self.cache[p+'__'+s]['path'],parameters=unique,detailed=detailed) for p in self.spec['periods'] for s in ('mc','background')]
        outputs=self.pool.run('evaluate',jobs,f'eval_{self.calls:03d}_{label}');self.calls+=1
        hist={};selected={}
        for i,job,path in outputs:
            meta=json.loads(Path(str(path)+'.json').read_text());key=job['period']+'__'+job['sample']
            hist[key]=np.memmap(path,dtype='<i8',mode='r',shape=tuple(meta['shape']))
            selected[key]=np.asarray(meta['selected'],dtype=np.int64).reshape(len(unique),2,3)
        scores,by_variable=self.score(hist,unique)
        self.history.extend(dict(call=label,scenario=i,parameters=list(p),**scores[i]) for i,p in enumerate(unique))
        return unique,scores,by_variable,hist,selected
    def weighted(self,hist,p,i,half,t,g,v):
        n=self.spec['variables'][v]['bins']+2;row=self.norm[p,t]
        totals=[];variances=[]
        for sample,scale in [('mc',row['signal_scale']),('background',row['background_scale'])]:
            count=hist[p+'__'+sample][i,half,t,g,v,:n].astype(float)
            weight=2*float(scale)/self.cache[p+'__'+sample]['probability']
            totals.append(count*weight);variances.append(count*weight*weight)
        return sum(totals),sum(variances)
    def score(self,hist,parameters):
        scores=[];records=[]
        for i,pars in enumerate(parameters):
            buckets={k:[] for k in ('train','validation_mc','heldout_observables')}
            for p in self.spec['periods']:
                for t in self.fit_topologies:
                    for v in TRAIN_OBS+VALID_OBS:
                        d=self.data[p+'__data'][t,0,v,:self.spec['variables'][v]['bins']+2].astype(float)
                        for half in (0,1):
                            m,var=self.weighted(hist,p,i,half,t,0,v)
                            raw=sum(hist[p+'__'+s][i,half,t,0,v,:len(d)].sum() for s in ('mc','background'))
                            m0,var0=self.weighted(hist,p,0,half,t,0,v)
                            raw0=sum(hist[p+'__'+s][0,half,t,0,v,:len(d)].sum() for s in ('mc','background'))
                            score,bins=shape_score(d,m,var,int(raw0),self.min_mc,m0,var0)
                            kind='heldout_observables' if v in VALID_OBS else 'validation_mc' if half else 'train'
                            if np.isfinite(score):buckets[kind].append(score)
                            records.append(dict(scenario=i,period=p,topology=TOPS[t],variable=OBS[v],half=half,kind=kind,score=score,usable_bins=bins,raw_mc=int(raw)))
            expected=len(self.spec['periods'])*len(self.fit_topologies)*len(TRAIN_OBS)
            if len(buckets['train'])!=expected:buckets['train']=[]
            scores.append({**{k:float(np.mean(values)) if values else float('inf') for k,values in buckets.items()},'training_projections':len(buckets['train'])})
        return scores,records

def shape_score(data,model,variance,raw_mc,min_mc,reference_model=None,reference_variance=None):
    # Shape-only diagnostic. The common area scale changes neither component ratio
    # nor any displayed fixed-normalization prediction; no formal p-value is used.
    if raw_mc<min_mc or data.sum()<=0 or model.sum()<=0:return float('nan'),0
    # Blocks of 4 original bins reduce MC sparsity, including finite under/overflow.
    cuts=np.arange(0,len(data),4)
    d=np.add.reduceat(data,cuts);m=np.add.reduceat(model,cuts);var=np.add.reduceat(variance,cuts)
    scale=d.sum()/m.sum();m*=scale;var*=scale*scale
    ref=m
    if reference_model is not None and reference_model.sum()>0:
        s0=data.sum()/reference_model.sum();ref=np.add.reduceat(reference_model,cuts)*s0
        var=np.add.reduceat(reference_variance,cuts)*s0*s0
    valid=(d+ref>=20)&(d+var>0)
    if valid.sum()<5:return float('nan'),int(valid.sum())
    return float(np.mean((d[valid]-m[valid])**2/(d[valid]+var[valid]))),int(valid.sum())

def best_training(parameters,scores):
    valid=[i for i,s in enumerate(scores) if np.isfinite(s['train'])]
    if not valid:raise ValueError('No usable training projections; increase cache sample or inspect inputs')
    return min(valid,key=lambda i:(scores[i]['train'],sum(abs(x) for x in parameters[i])))

def replay_check(engine,hist,selected,parameters):
    rows=[]
    for p in engine.spec['periods']:
        for sample in ('mc','background'):
            key=p+'__'+sample;prob=engine.cache[key]['probability']
            for t,top in enumerate(TOPS):
                full=engine.spec['reference_flows'][key]['selected'][t]
                estimated=float(selected[key][0,:,t].sum()/prob)
                err=np.sqrt(full*(1-prob)/prob) if prob<1 else 0.
                pull=(estimated-full)/err if err else 0. if abs(estimated-full)<1e-6 else float('inf')
                rows.append(dict(period=p,sample=sample,topology=top,reference_selected=full,sampled_estimate=estimated,sampling_error=err,pull=pull))
    return rows

def detailed_report(engine,output,parameters,scores,hist,selected):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    best=tuple(parameters[-1])
    # Caller includes nominal and final candidate only; look up final index.
    bi=len(parameters)-1;spec=engine.spec;rows=[];yield_rows=[]
    for p in spec['periods']:
        for t,top in enumerate(TOPS):
            for i,name in [(0,'nominal'),(bi,'candidate')]:
                total=0.
                for sample,scale in [('mc',engine.norm[p,t]['signal_scale']),('background',engine.norm[p,t]['background_scale'])]:
                    key=p+'__'+sample;total+=selected[key][i,:,t].sum()/engine.cache[key]['probability']*float(scale)
                yield_rows.append(dict(period=p,topology=top,scenario=name,predicted_selected=float(total),data_selected=spec['reference_flows'][p+'__data']['selected'][t]))
            fig,axes=plt.subplots(3,3,figsize=(15,11),constrained_layout=True)
            for vi,v in enumerate(spec['variables']):
                ax=axes.flat[vi];edges=diag.variable_edges(v);n=v['bins']+2
                d=engine.data[p+'__data'][t,0,vi,:n].astype(float);denom=d.sum()
                prediction=[]
                for i in (0,bi):
                    a,_=engine.weighted(hist,p,i,0,t,0,vi);b,_=engine.weighted(hist,p,i,1,t,0,vi);prediction.append((a+b)/2)
                for values,color,label in [(d,'#205cab','DATA'),(prediction[0],'#d45c18','Nominal MC'),(prediction[1],'#22946b','Candidate MC')]:
                    y=values[1:-1]/denom/np.diff(edges) if denom else values[1:-1]
                    ax.stairs(y,edges,color=color,label=label)
                q=[diag.quantile(h[1:-1],edges,x) for h in [d,*prediction] for x in (.005,.995)]
                q=[x for x in q if np.isfinite(x)]
                if q:
                    lo,hi=min(q),max(q);pad=max((hi-lo)*.1,np.min(np.diff(edges))*5)
                    ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                ax.set(xlabel=v['label'],ylabel='Density / DATA total');ax.legend(fontsize=8);ax.grid(alpha=.15)
            for ax in axes.flat[7:]:ax.axis('off')
            axes.flat[7].text(0,1,'Independent period/polarity response.\nOpening angle/pT and mass eγ\nwere not used to choose parameters.\nComponent normalizations remain fixed.\nCandidate is a detector-response test,\nnot an accepted production correction.',va='top')
            fig.suptitle(p+' · '+top+' · particle-level response test');destination=output/'plots'/p;destination.mkdir(parents=True,exist_ok=True)
            fig.savefig(destination/(top+'.png'),dpi=110);plt.close(fig)
            for g,stratum in enumerate(spec['strata']):
                for vi,v in enumerate(spec['variables']):
                    edges=diag.variable_edges(v);n=v['bins']+2;d=engine.data[p+'__data'][t,g,vi,:n].astype(float)
                    dm=diag.shape_metrics(d,edges,0,0)
                    for i,name in [(0,'nominal'),(bi,'candidate')]:
                        m0,var0=engine.weighted(hist,p,i,0,t,g,vi);m1,var1=engine.weighted(hist,p,i,1,t,g,vi);m=(m0+m1)/2
                        variance=(var0+var1)/4;neff=float(m.sum()**2/variance.sum()) if variance.sum()>0 else 0.
                        mm=diag.shape_metrics(m,edges,0,0)
                        flags=[]
                        if dm['n_finite']<300:flags.append('low DATA statistics')
                        if neff<engine.min_mc:flags.append('low effective MC statistics')
                        if min(dm['in_range_fraction'],mm['in_range_fraction'])<.99:flags.append('range truncation')
                        rows.append(dict(effective_mc=neff,status='; '.join(flags) or 'usable descriptive comparison',period=p,topology=top,axis=stratum['axis'],stratum=stratum['label'],variable=v['name'],scenario=name,
                            data_n=dm['n_finite'],predicted_n=mm['n_finite'],data_core=dm['core_spread'],mc_core=mm['core_spread'],
                            width_ratio=dm['core_spread']/mm['core_spread'] if mm['core_spread']>0 else float('nan'),
                            median_shift_MC_spreads=(dm['median']-mm['median'])/mm['core_spread'] if mm['core_spread']>0 else float('nan'),
                            data_in_range=dm['in_range_fraction'],mc_in_range=mm['in_range_fraction']))
    table(output/'final_shapes.csv',rows);table(output/'selection_yields.csv',yield_rows)
    with (output/'single_parameter_ranking.csv').open() as f:ranking=list(csv.DictReader(f))
    labels=[r['parameter'].replace('_sigma_fraction',' σ').replace('_sigma_rad',' σ').replace('_bias_fraction',' bias').replace('_',' ') for r in ranking]
    training=[1-float(r['train_score'])/float(r['nominal_train']) for r in ranking]
    validation=[1-float(r['heldout_observables'])/float(r['nominal_heldout_observables']) for r in ranking]
    fig,ax=plt.subplots(figsize=(11,7),constrained_layout=True);y=np.arange(len(labels))
    ax.barh(y-.17,training,height=.32,label='Training-observable improvement',color='#205cab')
    ax.barh(y+.17,validation,height=.32,label='Held-out-observable improvement',color='#d45c18')
    ax.set(yticks=y,yticklabels=labels,xlabel='Fractional reduction of shape-discrepancy score');ax.axvline(0,color='gray');ax.invert_yaxis();ax.legend();ax.grid(axis='x',alpha=.15)
    fig.savefig(output/'single_parameter_effects.png',dpi=120);plt.close(fig)
    links=''.join(f'<li>{p}: '+ ' · '.join(f'<a href="plots/{p}/{t}.png">{t}</a>' for t in TOPS)+'</li>' for p in spec['periods'])
    (output/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Particle response test</title><style>body{font:17px system-ui;max-width:1100px;margin:40px auto}li{margin:15px}</style><h1>Particle-response test</h1><p>Start with CD_FT. Blue: DATA; orange: nominal MC; green: candidate MC. Normalizations are fixed; no per-plot matching is applied. This is a candidate additional Gaussian response and momentum-scale model, not a production correction.</p><p><a href="single_parameter_effects.png">Which particle variations improve the distributions?</a></p><ul>'+links+'</ul><p><a href="REPORT.md">Fit definitions, limitations and parameter values</a> · <a href="single_parameter_ranking.csv">Single-parameter ranking</a> · <a href="final_scores.csv">Training/validation scores</a> · <a href="selection_yields.csv">Selection yields</a> · <a href="final_shapes.csv">All joint-cell shapes</a></p>')

def discard_histograms(hist):
    """Release mappings and unlink only this evaluation's expendable binaries."""
    files=set()
    for array in hist.values():
        if isinstance(array,np.memmap):
            files.add(Path(array.filename));array._mmap.close()
    hist.clear()
    for path in files:path.unlink()

def scratch_directory(root,output):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=Path(output).name+'-',dir=root))


def compact_results(output,periods):
    destination=output/'results';destination.mkdir(exist_ok=True)
    for period in periods:
        target=destination/period;target.mkdir(exist_ok=True)
        for path in (output/period).iterdir():
            if path.is_file() and path.suffix in ('.csv','.json','.md','.html','.png'):shutil.copy2(path,target/path.name)
        shutil.copytree(output/period/'plots',target/'plots',dirs_exist_ok=True)
    for name in ('REPORT.html','manifest.json','cache_manifest.json','reconstruction_audit.csv','scratch_location.json'):shutil.copy2(output/name,destination/name)
    return destination


def publish_results(output,spec,args,reference):
    links=''.join(f'<li><a href="{p}/REPORT.html">{p}: independent response test</a></li>' for p in spec['periods'])
    (output/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Independent period response tests</title><style>body{font:18px system-ui;max-width:1000px;margin:40px auto}li{padding:10px}</style><h1>Independent particle-response tests</h1><p>Every period and polarity has its own fitted response parameters. No Spring/Fall or inbending/outbending parameters are pooled. Each fit is checked with an independent MC subset and held-out observables from the same sample.</p><ul>'+links+'</ul>')
    dump(output/'manifest.json',dict(version=VERSION,reference_manifest_sha256=digest(reference/'manifest.json'),arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},source_hashes={n:digest(SOURCE/n) for n in ('study_particle_response.py','particle_response_worker.cpp','particle_response_core.h')},independent_periods=spec['periods']))
    compact=compact_results(output,spec['periods'])
    print('[done] Open '+str(compact/'REPORT.html')+'; transfer the results directory only, retain the MC cache on the farm',flush=True)


def recover_reports(output,reference,data,norm):
    saved=json.loads((output/'spec.json').read_text())
    cache_meta=json.loads((output/'cache_manifest.json').read_text())
    if cache_meta['reference_sha256']!=digest(reference/'manifest.json'):raise ValueError('Recovery reference differs from the completed fits')
    if cache_meta['version']!=VERSION:raise ValueError('Incompatible fit cache version')
    cache=cache_meta['samples']
    work=Path(json.loads((output/'scratch_location.json').read_text())['directory'])
    recovered=[]
    for period in saved['periods']:
        folder=output/period
        required=[folder/'REPORT.html',folder/'final_shapes.csv',folder/'selection_yields.csv',folder/'single_parameter_effects.png',*[folder/'plots'/period/(t+'.png') for t in TOPS]]
        if (folder/'report_complete.json').is_file() and all(p.is_file() for p in required):
            marker=json.loads((folder/'report_complete.json').read_text())
            if marker['candidate_sha256']!=digest(folder/'candidate_parameters.json'):raise ValueError('Candidate parameters changed after reporting: '+period)
            print('[recover] Existing completed reports retained: '+period,flush=True);continue
        candidate=json.loads((folder/'candidate_parameters.json').read_text())
        spec=json.loads((folder/'spec.json').read_text())
        if spec['periods']!=[period]:raise ValueError('Recovery period configuration mismatch')
        locations=list((work/period).glob('eval_*_final'))
        if len(locations)!=1:raise ValueError('Cannot identify completed final evaluation for '+period)
        hist={};selected={};parameters=None
        for job_path in sorted(locations[0].glob('*.json')):
            if job_path.name.endswith('.bin.json'):continue
            job=json.loads(job_path.read_text())
            if 'parameters' not in job:continue
            if job['period']!=period or not job.get('detailed'):raise ValueError('Unexpected final response job')
            path=job_path.with_suffix('.bin');meta=json.loads(Path(str(path)+'.json').read_text())
            if not path.is_file() or path.stat().st_size!=int(np.prod(meta['shape']))*8:raise ValueError('Final trial binary missing/incomplete: '+str(path))
            if parameters is None:parameters=job['parameters']
            elif parameters!=job['parameters']:raise ValueError('Signal/leakage parameter mismatch')
            key=period+'__'+job['sample'];hist[key]=np.memmap(path,dtype='<i8',mode='r',shape=tuple(meta['shape']))
            selected[key]=np.asarray(meta['selected'],dtype=np.int64).reshape(len(parameters),2,3)
        if set(hist)!={period+'__mc',period+'__background'}:raise ValueError('Incomplete final signal/leakage trials')
        expected=[candidate['parameters'].get(name,0.) for name in PARAMS]
        if list(parameters[-1])!=expected or any(parameters[0]):raise ValueError('Saved candidate differs from final trial')
        engine=Engine(None,spec,cache,data,norm,candidate['fit_topologies'],100)
        recovered.append((engine,folder,parameters,[],hist,selected))
    return saved,recovered


def render_completed(results):
    for engine,folder,pars,scores,hist,selected in results:
        detailed_report(engine,folder,pars,scores,hist,selected)
        dump(folder/'report_complete.json',dict(candidate_sha256=digest(folder/'candidate_parameters.json')))
        discard_histograms(hist)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-study',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'),help='Large MC caches and temporary evaluation files; each run gets a unique subdirectory')
    parser.add_argument('--finish-reports',action='store_true',help='Finish reports from saved final trials in the existing --output; no ROOT reads or refitting')
    parser.add_argument('--workers',type=int,default=10)
    parser.add_argument('--expected-cache-events',type=int,default=200000,help='Expected uniform pre-cut MC events per period/sample, not a selected-event cap')
    parser.add_argument('--reuse-cache',type=Path,help='Previous particle-response output containing a verified event cache')
    parser.add_argument('--seed',type=int,default=173205)
    parser.add_argument('--fit-topologies',choices=TOPS,nargs='+',default=['CD_FT'])
    parser.add_argument('--coordinate-rounds',type=int,default=1)
    parser.add_argument('--refinement-rounds',type=int,default=1)
    parser.add_argument('--min-mc',type=int,default=100)
    args=parser.parse_args(argv)
    if not 1<=args.workers<=10 or args.expected_cache_events<1000 or not 0<=args.seed<2**64 or args.coordinate_rounds<0 or args.refinement_rounds<0 or args.min_mc<10:parser.error('Invalid worker/cache/seed/fit limits')
    reference=args.reference_study.resolve();output=args.output.resolve()
    if output.exists() and not args.finish_reports:parser.error('Output already exists; choose a new directory')
    manifest=json.loads((reference/'manifest.json').read_text());spec=copy.deepcopy(manifest['spec'])
    if not spec.get('background_aware') or not spec.get('joint_grids'):parser.error('Use the completed background-aware joint diagnostic as reference')
    if digest(reference/'histograms.npz')!=manifest['histogram_sha256']:raise ValueError('Reference histogram checksum mismatch')
    spec.update(expected_cache_events=args.expected_cache_events,seed=args.seed,reference_flows=manifest['flows'])
    for k,n in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        if digest(reference/'backend'/n)!=spec['input_hashes'][k]:raise ValueError('Saved global-cut source hash mismatch')
    with np.load(reference/'histograms.npz',allow_pickle=False) as archive:
        data={p+'__data':archive[p+'__data'].copy() for p in spec['periods']}
    with (reference/'normalizations.csv').open() as f:
        norm={(r['period'],TOPS.index(r['topology'])):r for r in csv.DictReader(f)}
    if any(r['status']!='ok' for r in norm.values()):raise ValueError('Invalid reference component normalization')
    if args.finish_reports:
        if not output.is_dir():parser.error('Report recovery requires an existing output directory')
        old,finished=recover_reports(output,reference,data,norm)
        render_completed(finished);publish_results(output,old,args,reference);return
    output.mkdir(parents=True);dump(output/'spec.json',spec)
    print('[build] Compiling cached particle-response worker',flush=True)
    worker=compile_worker(output,reference)
    work=scratch_directory(args.scratch_root,output)
    dump(output/'scratch_location.json',dict(directory=str(work)))
    print('[scratch] '+str(work),flush=True)
    pool=Pool(worker,output/'spec.json',work,args.workers)
    cache=build_cache(pool,spec,manifest,reference,args.reuse_cache.resolve() if args.reuse_cache else None)
    if not args.reuse_cache:shutil.copy2(work/'cache_manifest.json',output/'cache_manifest.json')
    # Preserve cache provenance when reusing an earlier immutable cache.
    if args.reuse_cache:
        (output/'cache_manifest.json').write_text(json.dumps(dict(version=VERSION,reference_sha256=digest(reference/'manifest.json'),seed=args.seed,samples={k:dict(r,file=r['path']) for k,r in cache.items()}),indent=2))
    audit_rows=[]
    for key,record in cache.items():
        audits=record['audit'];number=sum(a['audited'] for a in audits)
        for j,name in enumerate(['x','Q2','t1','open_angle_ep2','pTmiss',*OBS]):
            audit_rows.append(dict(sample=key,variable=name,audited_events=number,
                mean_absolute_offset=sum(a['absolute_difference_sum'][j] for a in audits)/number if number else float('nan'),
                maximum_absolute_offset=max(a['absolute_difference_max'][j] for a in audits)))
    table(output/'reconstruction_audit.csv',audit_rows)
    budget=pool.budget
    def fit_period(period):
        folder=output/period;folder.mkdir()
        one=copy.deepcopy(spec);one['periods']=[period];dump(folder/'spec.json',one)
        period_work=work/period;period_work.mkdir()
        local_pool=Pool(worker,folder/'spec.json',period_work,args.workers,budget)
        engine=Engine(local_pool,one,cache,data,norm,args.fit_topologies,args.min_mc)
        return fit_candidate(engine,folder,args)
    executor=ThreadPoolExecutor(max_workers=min(args.workers,len(spec['periods'])))
    futures=[executor.submit(fit_period,p) for p in spec['periods']];results=[]
    try:
        for future in as_completed(futures):results.append(future.result())
    except BaseException:
        budget.cancel_all()
        for future in futures:future.cancel()
        raise
    finally:executor.shutdown(wait=True,cancel_futures=True)
    render_completed(results)
    publish_results(output,spec,args,reference)


def fit_candidate(engine,output,args):
    spec=engine.spec
    grid=parameter_grid(args.fit_topologies);nominal=[0.]*20;candidates=[nominal]
    for k,levels in grid.items():
        for value in levels:
            if value:
                p=list(nominal);p[k]=value;candidates.append(p)
    print(f'[fit] {len(candidates)} initial scenarios; fit {args.fit_topologies}; independent fit for {spec['periods'][0]}',flush=True)
    pars,scores,projections,hist,selected=engine.evaluate(candidates,'single_parameters')
    replay=replay_check(engine,hist,selected,pars);table(output/'nominal_replay.csv',replay)
    if any(abs(r['pull'])>7 for r in replay):raise ValueError('Nominal cache replay differs from reference beyond sampling tolerance; inspect nominal_replay.csv')
    table(output/'single_parameter_scores.csv',[dict(scenario=i,parameters=json.dumps(p),**scores[i]) for i,p in enumerate(pars)])
    table(output/'single_parameter_projection_scores.csv',projections)
    ranking=[]
    for k,levels in grid.items():
        indices=[i for i,p in enumerate(pars) if p[k]!=0 and sum(x!=0 for x in p)==1]
        b=min(indices,key=lambda i:scores[i]['train'])
        ranking.append(dict(parameter=PARAMS[k],nominal_heldout_observables=scores[0]['heldout_observables'],best_single_value=pars[b][k],train_score=scores[b]['train'],nominal_train=scores[0]['train'],improvement=scores[0]['train']-scores[b]['train'],heldout_observables=scores[b]['heldout_observables'],validation_mc=scores[b]['validation_mc']))
    table(output/'single_parameter_ranking.csv',sorted(ranking,key=lambda r:r['train_score']))
    # Finite template-response similarities flag ambiguous track attribution.
    if hist:
        vectors=[];names=[]
        for k,levels in grid.items():
            candidates=[i for i,p in enumerate(pars) if p[k]!=0 and sum(x!=0 for x in p)==1]
            idx=max(candidates,key=lambda i:pars[i][k]);features=[]
            for period in spec['periods']:
                for ti in engine.fit_topologies:
                    for vi in TRAIN_OBS:
                        n=spec['variables'][vi]['bins']+2;data=engine.data[period+'__data'][ti,0,vi,:n].astype(float)
                        m0,v0=engine.weighted(hist,period,0,0,ti,0,vi);m,_=engine.weighted(hist,period,idx,0,ti,0,vi)
                        if m0.sum()>0 and m.sum()>0:
                            scale0=data.sum()/m0.sum();scale=data.sum()/m.sum()
                            features.extend(((m*scale-m0*scale0)/np.sqrt(data+v0*scale0*scale0+1)).tolist())
            vectors.append(np.asarray(features));names.append(PARAMS[k])
        similarities=[]
        for a in range(len(names)):
            for b in range(a):
                if len(vectors[a])!=len(vectors[b]):continue
                den=np.linalg.norm(vectors[a])*np.linalg.norm(vectors[b])
                cosine=float(np.dot(vectors[a],vectors[b])/den) if den else float('nan')
                similarities.append(dict(parameter_a=names[a],parameter_b=names[b],template_response_cosine=cosine))
        table(output/'response_similarity.csv',similarities)
    best=list(pars[best_training(pars,scores)])
    discard_histograms(hist)
    del hist,selected
    steps=[]
    for phase,rounds in [('coordinate',args.coordinate_rounds),('refine',args.refinement_rounds)]:
        for iteration in range(rounds):
            for k,levels in grid.items():
                values=levels
                if phase=='refine':
                    lower=max([x for x in levels if x<best[k]],default=levels[0]);upper=min([x for x in levels if x>best[k]],default=levels[-1])
                    values=sorted(set([best[k],(lower+best[k])/2,(upper+best[k])/2]))
                trials=[list(best)]
                for value in values:
                    p=list(best);p[k]=value;trials.append(p)
                p,s,_,h,se=engine.evaluate(trials,f'{phase}_{iteration}_{k}')
                chosen=best_training(p,s);best=list(p[chosen]);steps.append(dict(phase=phase,round=iteration,parameter=PARAMS[k],chosen_value=best[k],**s[chosen]))
                discard_histograms(h)
                del h,se
    table(output/'fit_steps.csv',steps)
    p,s,projections,h,se=engine.evaluate([nominal,best],'final',detailed=True)
    table(output/'final_scores.csv',[dict(scenario='nominal',**s[0]),dict(scenario='candidate',**s[-1])])
    table(output/'final_projection_scores.csv',projections)
    active={PARAMS[k]:best[k] for k in grid};boundary=[PARAMS[k] for k,levels in grid.items() if best[k]!=0 and best[k] in (levels[0],levels[-1])]
    dump(output/'candidate_parameters.json',dict(version=VERSION,parameters=active,boundary_parameters=boundary,period=spec['periods'][0],fit_topologies=args.fit_topologies,training_observables=[OBS[i] for i in TRAIN_OBS],validation_observables=[OBS[i] for i in VALID_OBS]))
    table(output/'evaluation_history.csv',[dict(call=r['call'],scenario=r['scenario'],parameters=json.dumps(r['parameters']),**{k:v for k,v in r.items() if k not in ('call','scenario','parameters')}) for r in engine.history])
    units=[]
    for k,value in active.items():
        unit='fraction' if k.endswith('fraction') else 'radian'
        units.append(dict(parameter=k,value=value,unit=unit,display_value=value*(100 if unit=='fraction' else 1000),display_unit='percent' if unit=='fraction' else 'mrad'))
    table(output/'candidate_parameters.csv',units)
    parameter_text='\n'.join(f'- `{k}` = {v:.6g}' for k,v in active.items())
    (output/'REPORT.md').write_text('''# Particle-response test

Open REPORT.html, then CD_FT for each period. Compare blue DATA with orange
nominal MC and green candidate MC. Component scales stay fixed to the reference.
selection_yields.csv shows actual predicted selection changes, not renormalized
agreement. final_shapes.csv retains inclusive and all joint-cell comparisons.

This first search adds independent Gaussian errors to reconstructed particle
momentum magnitude, polar angle and azimuth, plus optional fractional momentum
biases. Magnitude sigmas are fractions; angular sigmas are radians. These are
ADDITIONAL widths, not multipliers of the existing resolution. No generated
residual is amplified and no physical energy-loss tail is deliberately scaled.
The model does not include correlated or non-Gaussian detector errors. response_similarity.csv compares the finite template perturbations: a cosine
near +1 flags variations producing similar changes in the training projections.
This is a sensitivity diagnostic, not a parameter covariance measurement. A good
fit does not uniquely determine which physical detector defect is responsible.

Default tuning uses CD_FT only, separately for each individual period/polarity,
with coplanarity, signed missing energy and signed missing mass squared. No
parameters are pooled between periods or polarities. Photon
opening angle, pT and missing mass eγ are held out in that same period. Every period predicts all three topologies. Within that period only, the
electron FD response is shared across topologies;
FD-proton and FD-photon parameters remain zero unless those topologies are
explicitly included with --fit-topologies. Do not infer their resolutions from
unfitted parameter values. Other topologies are background-sensitive checks.

The objective is a descriptive mean weighted shape discrepancy over equal-weight
projections, not a chi-square probability: variables share events. Four original
histogram bins form one score bin. The model sum is area-normalized ONLY in the
shape score; the signal/background ratio stays fixed. Displayed predictions and
selection yields retain their fixed reference scales. The score uses fixed
nominal sampling variance/masks to avoid favoring a trial merely because it has
less MC precision. Every training projection must have sufficient nominal MC.
No fit uncertainty or production systematic is inferred from this score.

The uniform MC cache is sampled BEFORE global/exclusivity cuts. Both signal and
leakage are varied, including events which can migrate into the selection.
Cache probability is retained in weights. Deterministic independent MC halves
provide a sampling check; both compare to the same DATA histograms, so this is
not independent DATA validation. Held-out observables and joint-cell checks provide the
additional shape checks. Common random numbers stabilize comparisons between
parameter trials. Original candidate reconstruction inefficiency cannot be
repaired by smearing already reconstructed candidates.

Derived quantities are updated by stored_value + recomputed_shift, with the
shift calculated from perturbed versus original particle four-vectors. This
preserves existing rounding/calibration offsets and makes zero response exactly
identical to the original. Cache metadata records direct recomputation offsets
for auditing. Global cuts and each omitted-cut projection are reapplied, and
joint cells are reassigned using the updated measured kinematics. A nominal
selected-yield replay must agree within sampling tolerance before fitting.

The search compares one-parameter scans, then coordinate grid updates and a
local refinement. This is a bounded candidate search, not a proof of a global
optimum. Parameters at a scan boundary need a wider test, not automatic adoption.
Reference background-transfer, radiation and within-cell physics-population
limitations persist. Inspect several observables and held-out results before
accepting any particle-level interpretation. Production files are unchanged.

## Candidate additional-response parameters

'''+parameter_text+'\n\nBoundary parameters: '+(', '.join(boundary) or 'none')+'\n')
    return engine,output,p,s,h,se

if __name__=='__main__':main()
