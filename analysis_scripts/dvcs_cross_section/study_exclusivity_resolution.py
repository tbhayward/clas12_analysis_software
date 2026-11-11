#!/usr/bin/env python3
"""Production-cut DATA/MC exclusivity diagnostics, with one cut left out at a time.

Native workers fill every projection in one streaming file pass. This measures
conditional distribution shapes, not individual track resolutions or corrections.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import html
import json
import math
import re
from pathlib import Path
import shlex
import shutil
import subprocess
import threading
import numpy as np

SOURCE=Path(__file__).resolve().parent
VERSION='exclusivity-resolution-v1'
PERIODS={'fa18_inb':('Fa18_Inb','50nA'),'fa18_out':('Fa18_Out','50nA'),
         'sp18_inb':('Sp18_Inb','50nA'),'sp18_out':('Sp18_Out','45nA'),'sp19_inb':('Sp19_Inb','50nA')}
TOPOLOGIES=['FD_FD','CD_FD','CD_FT']
# Units checked against processing_classes/src/analyzers/kinematic_variables.java:
# theta_gamma_gamma is degrees; Emiss2 is missing energy (GeV), NOT energy squared.
VARIABLES=[dict(name=n,label=l,bins=b,low=lo,high=hi) for n,l,b,lo,hi in [
    ('Delta_phi','Coplanarity Δφ (rad)',600,math.pi-.6,math.pi+.6),
    ('theta','CM angle θ (rad)',640,0,math.pi),
    ('theta_gamma_gamma','Photon opening angle (deg)',400,0,8),
    ('pTmiss','Missing transverse momentum (GeV)',400,0,.8),
    ('Emiss2','Missing energy (GeV)',600,-3,3),
    ('Mx2','Missing mass squared epγ (GeV²)',600,-.12,.12),
    ('Mx2_2','Missing mass squared eγ (GeV²)',500,-1,4)]]
AXES={
    'e_p':dict(label='Electron momentum (GeV)',edges=[0,2,3,4,5,7,12,100]),
    'e_theta':dict(label='Electron polar angle (deg)',edges=[0,8,12,16,20,25,40,180.000001]),
    'p1_p':dict(label='Proton momentum (GeV)',edges=[0,.3,.42,.5,.65,.85,1.2,2,100]),
    'p1_theta':dict(label='Proton polar angle (deg)',edges=[0,25,37,45,50,55,60,64.23,180.000001]),
    'p2_p':dict(label='Photon energy (GeV)',edges=[0,2,4,6,8,12,100]),
    'p2_theta':dict(label='Photon polar angle (deg)',edges=[0,5,8,12,20,30,40,180.000001])}

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def saved_window(stats):
    upper=stats.get('mode')=='upper_quantile'
    lo,hi=stats.get('cut_low'),stats.get('cut_high')
    if not upper and (lo is None or hi is None or not math.isfinite(lo) or not math.isfinite(hi) or hi<=lo):
        mean,sd=stats.get('mean',float('nan')),stats.get('std',float('nan'))
        if not math.isfinite(mean) or not math.isfinite(sd) or sd<=0:raise ValueError('Invalid saved cut')
        lo,hi=mean-3*sd,mean+3*sd
    if hi is None or not math.isfinite(hi):raise ValueError('Invalid saved upper cut')
    return [None if upper else lo,hi]

def configuration(path=None):
    variables=json.loads(json.dumps(VARIABLES));axes=json.loads(json.dumps(AXES))
    if path:
        custom=json.loads(Path(path).read_text())
        if set(custom)-{'variables','axes'}:raise ValueError('Configuration accepts only variables and axes')
        unknown=set(custom.get('variables',{}))-{v['name'] for v in variables}
        if unknown or set(custom.get('axes',{}))-set(axes):raise ValueError('Unknown variable or axis')
        for v in variables:
            updates=custom.get('variables',{}).get(v['name'],{})
            if set(updates)-{'bins','low','high'}:raise ValueError('Variable overrides: bins, low, high only')
            v.update(updates)
        for name,edges in custom.get('axes',{}).items():axes[name]['edges']=edges
    for v in variables:
        if not isinstance(v['bins'],int) or v['bins']<20 or v['bins']>10000 or not math.isfinite(v['low']) or not math.isfinite(v['high']) or v['high']<=v['low']:
            raise ValueError('Invalid histogram configuration')
    strata=[dict(axis='inclusive',label='Inclusive')];indices=[-1]
    for ai,(name,axis) in enumerate(axes.items()):
        e=axis['edges']
        if len(e)<2 or any(not math.isfinite(x) for x in e) or any(b<=a for a,b in zip(e,e[1:])):raise ValueError('Invalid axis edges')
        for lo,hi in zip(e,e[1:]):
            strata.append(dict(axis=name,low=lo,high=hi,label=f'{lo:g}–{hi:g}'));indices.append(ai)
    return variables,axes,strata,indices

def build_spec(run,periods,config=None):
    run=Path(run).resolve()
    paths={k:run/p for k,p in {'global_cuts':'output/jsons/global_cuts_config.json',
        'exclusivity':'output/jsons/combined_cuts.json','global_cpp':'global_cuts.cpp',
        'global_header':'global_cuts.h','domain_header':'kinematic_bin_integration.h'}.items()}
    for p in paths.values():
        if not p.is_file():raise FileNotFoundError(f'Missing saved production file: {p}')
    fields=set(re.findall(r'^\s*(?:bool|double|int|std::vector<int>)\s+(\w+)\s*(?:=|;)',paths['global_header'].read_text().split('struct GlobalCutConfig {',1)[1].split('\n};',1)[0],re.M))
    loaded=set(re.findall(r'LOAD\((\w+)\)',(SOURCE/'exclusivity_resolution_worker.cpp').read_text()))-{'name'}
    raw_global=json.loads(paths['global_cuts'].read_text())
    if fields!=loaded or fields-set(raw_global):
        raise ValueError('GlobalCutConfig schema differs from diagnostic loader; update it before proceeding')
    variables,axes,strata,indices=configuration(config)
    spec=dict(version=VERSION,production_run=str(run),periods=periods,
        period_codes={p:PERIODS[p][0] for p in periods},variables=variables,axes=axes,strata=strata,
        strata_axis_index=indices,input_hashes={k:digest(p) for k,p in paths.items()},
        global_cuts=json.loads(paths['global_cuts'].read_text()),exclusivity=json.loads(paths['exclusivity'].read_text()),
        weighting='Unit-weighted reconstructed DVCS candidates; no background subtraction or cross-section corrections')
    allowed={v['name'] for v in variables}
    for p in periods:
        for t in TOPOLOGIES:
            for sample in ('data','mc'):
                block=spec['exclusivity'][f'DVCS_{PERIODS[p][0]}_{t}'][sample]
                if 'Mx2' not in block or set(block)-allowed:raise ValueError(f'Invalid saved cuts: {p}/{t}/{sample}')
                for stats in block.values():saved_window(stats)
    return spec,paths

def inventory(root,periods,limit):
    result={}
    for p in periods:
        for sample in ('data','mc'):
            directory=(root/'data'/p/'dvcs_data/parts' if sample=='data' else root/'mc'/p/'dvcs_mc'/PERIODS[p][1]/'parts')
            pattern='*.root' if sample=='data' else '*.rec_only.root'
            files=sorted(f.resolve() for f in directory.glob(pattern) if '.partial' not in f.name)
            if limit:files=files[:limit]
            if not files:raise FileNotFoundError(f'No {sample} reconstructed parts matching {directory/pattern}')
            result[p,sample]=files
    return result

def make_jobs(files,workers):
    jobs=[]
    for (p,sample),paths in files.items():
        count=min(len(paths),max(1,math.ceil(2*workers/len(files))))
        batches=[[] for _ in range(count)];sizes=[0]*count
        for path in sorted(paths,key=lambda f:f.stat().st_size,reverse=True):
            k=int(np.argmin(sizes));batches[k].append(str(path));sizes[k]+=path.stat().st_size
        jobs.extend(dict(period=p,sample=sample,files=sorted(b)) for b in batches)
    return jobs

def compile_worker(output,paths):
    build=output/'backend';build.mkdir()
    for name in ('global_cpp','global_header','domain_header'):shutil.copy2(paths[name],build/paths[name].name)
    for name in ('exclusivity_resolution_worker.cpp','exclusivity_resolution_core.h'):shutil.copy2(SOURCE/name,build/name)
    flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True))
    libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
    worker=build/'exclusivity_resolution_worker'
    subprocess.run([os.environ.get('CXX','c++'),*flags,'-O3','-std=c++17','-I'+str(build),str(build/'exclusivity_resolution_worker.cpp'),str(build/'global_cuts.cpp'),'-o',str(worker),*libs],check=True)
    return worker

def scan(output,spec,paths,files,workers):
    print('[build] Compiling diagnostic with saved production global cuts',flush=True)
    worker=compile_worker(output,paths);jobs=make_jobs(files,workers)
    for name in ('global_cpp','global_header','domain_header'):
        if digest(output/'backend'/paths[name].name)!=spec['input_hashes'][name]:
            raise RuntimeError('Saved production sources changed while the study was starting')
    initial_inventory=[dict(period=p,sample=s,path=str(f),size=f.stat().st_size,mtime_ns=f.stat().st_mtime_ns) for (p,s),paths_ in files.items() for f in paths_]
    temp=output/'batches';temp.mkdir();write_json(output/'spec.json',spec)
    arrays={};flows={};file_records=[]
    ng=len(spec['strata']);maxbins=max(v['bins'] for v in spec['variables'])+2
    for p,s in files:
        key=f'{p}__{s}'
        arrays[key]=np.zeros((3,ng,7,maxbins),dtype=np.uint64)
        arrays[key+'__stats']=np.zeros((3,ng,7,2),dtype=np.uint64)
        flows[key]=dict(entries=0,global_pass=[0]*3,selected=[0]*3)
    active=set();lock=threading.Lock();cancel=threading.Event()
    def launch(item):
        i,job=item;jp=temp/f'{i}.json';out=temp/f'{i}.result.json';log=temp/f'{i}.log'
        write_json(jp,job)
        with log.open('w') as f:
            with lock:
                if cancel.is_set():raise RuntimeError('Study cancelled after worker failure')
                process=subprocess.Popen([str(worker),str(output/'spec.json'),str(jp),str(out)],stdout=f,stderr=subprocess.STDOUT)
                active.add(process)
            code=process.wait()
            with lock:active.discard(process)
        if code:raise RuntimeError(f'Worker failed for {job["period"]}/{job["sample"]}; see {log}\n'+log.read_text()[-3000:])
        return i,job,out
    print(f'[study] {len(jobs)} batches; at most {workers} workers TOTAL for data, MC and all periods',flush=True)
    pool=ThreadPoolExecutor(max_workers=workers)
    futures=[]
    try:
        futures=[pool.submit(launch,item) for item in enumerate(jobs)]
        for done,future in enumerate(as_completed(futures),1):
            i,job,path=future.result();payload=json.loads(path.read_text());key=f'{job["period"]}__{job["sample"]}'
            for h,record in enumerate(payload['histograms']):
                t,g,v=np.unravel_index(h,(3,ng,7));counts=record['counts']
                arrays[key][t,g,v,:len(counts)]+=np.asarray(counts,dtype=np.uint64)
                arrays[key+'__stats'][t,g,v]+=np.array([record['nonfinite'],record['pass_own_cut']],dtype=np.uint64)
            flows[key]['entries']+=payload['entries']
            for k in ('global_pass','selected'):flows[key][k]=[a+b for a,b in zip(flows[key][k],payload[k])]
            file_records.extend(dict(period=job['period'],sample=job['sample'],**f) for f in payload['files'])
            path.unlink()
            print(f'[study] {done}/{len(jobs)} batches complete: {job["period"]}/{job["sample"]}',flush=True)
    except BaseException:
        cancel.set()
        with lock:
            for p in active:p.terminate()
        for f in futures:f.cancel()
        raise
    finally:pool.shutdown(wait=True,cancel_futures=True)
    for record in initial_inventory:
        st=Path(record['path']).stat()
        if st.st_size!=record['size'] or st.st_mtime_ns!=record['mtime_ns']:
            raise RuntimeError('Input ROOT file changed while reading: '+record['path'])
    np.savez_compressed(output/'histograms.npz',**arrays)
    manifest=dict(spec=spec,flows=flows,files=sorted(file_records,key=lambda f:(f['period'],f['sample'],f['path'])),
        inventory=initial_inventory,
        histogram_sha256=digest(output/'histograms.npz'),source_hashes={n:digest(SOURCE/n) for n in ('study_exclusivity_resolution.py','exclusivity_resolution_worker.cpp','exclusivity_resolution_core.h','exclusivity_resolution_report.html')})
    write_json(output/'manifest.json',manifest)
    return arrays,manifest

# Histogram quantiles use within-bin interpolation; overflow is NEVER assigned a
# fabricated coordinate. Comparisons are flagged if appreciable mass lies outside.
def quantile(counts,edges,probability):
    total=float(np.sum(counts))
    if not total:return float('nan')
    target=probability*total;cum=np.cumsum(counts,dtype=float)
    i=min(len(counts)-1,int(np.searchsorted(cum,target,side='right')))
    prev=cum[i-1] if i else 0.
    return float(edges[i]+(edges[i+1]-edges[i])*(target-prev)/counts[i])

def shape_metrics(counts,edges,nonfinite,passed):
    n=int(np.sum(counts));inside=int(np.sum(counts[1:-1]));finite_fraction=inside/n if n else 0.
    qs=[quantile(counts[1:-1],edges,q) for q in (.025,.16,.5,.84,.975)]
    smooth=np.convolve(counts[1:-1].astype(float),[1,2,1],mode='same')
    peak=float((edges[np.argmax(smooth)]+edges[np.argmax(smooth)+1])/2) if inside else float('nan')
    return dict(n_finite=n,n_nonfinite=int(nonfinite),in_range_fraction=finite_fraction,
        underflow=int(counts[0]),overflow=int(counts[-1]),q025=qs[0],q16=qs[1],median=qs[2],q84=qs[3],q975=qs[4],
        core_spread=(qs[3]-qs[1])/2,peak=peak,pass_own_cut_fraction=passed/(n+nonfinite) if n+nonfinite else float('nan'))

def outside_fraction(counts,edges,lo,hi):
    total=float(np.sum(counts))
    if not total or not math.isfinite(lo) or not math.isfinite(hi):return float('nan')
    overlap=np.clip(np.minimum(edges[1:],hi)-np.maximum(edges[:-1],lo),0,None)/np.diff(edges)
    return float(1-np.dot(counts[1:-1],overlap)/total)

def comparisons(arrays,spec,min_entries):
    rows=[];detail={}
    for p in spec['periods']:
        for t,topology in enumerate(TOPOLOGIES):
            for g,stratum in enumerate(spec['strata']):
                for vi,v in enumerate(spec['variables']):
                    edges=np.linspace(v['low'],v['high'],v['bins']+1);metrics=[];hist=[]
                    for sample in ('data','mc'):
                        key=f'{p}__{sample}';counts=arrays[key][t,g,vi,:v['bins']+2]
                        stats=arrays[key+'__stats'][t,g,vi]
                        metrics.append(shape_metrics(counts,edges,*stats));hist.append(counts)
                    d,m=metrics;ratio=d['core_spread']/m['core_spread'] if m['core_spread']>0 else float('nan')
                    shift=(d['median']-m['median'])/m['core_spread'] if m['core_spread']>0 else float('nan')
                    dtail=outside_fraction(hist[0],edges,m['q025'],m['q975'])
                    mtail=outside_fraction(hist[1],edges,m['q025'],m['q975'])
                    flags=[]
                    if min(d['n_finite'],m['n_finite'])<min_entries:flags.append('low statistics')
                    if min(d['in_range_fraction'],m['in_range_fraction'])<.99:flags.append('range truncation')
                    if any(z['core_spread']<2*(edges[1]-edges[0]) for z in metrics):flags.append('core poorly resolved by histogram')
                    if d['n_nonfinite'] or m['n_nonfinite']:flags.append('nonfinite values')
                    row=dict(period=p,topology=topology,axis=stratum['axis'],stratum=stratum['label'],variable=v['name'],
                        width_ratio=ratio,median_shift_MC_spreads=shift,data_outside_MC95=dtail,mc_outside_MC95=mtail,
                        excess_tail_fraction=dtail-mtail,status='; '.join(flags) or 'usable shape comparison')
                    for prefix,z in zip(('data','mc'),metrics):row.update({prefix+'_'+k:value for k,value in z.items()})
                    rows.append(row)
                    detail[f'{p}|{topology}|{g}|{vi}']=dict(data=hist[0].tolist(),mc=hist[1].tolist(),summary=row)
    return rows,detail

def effective_window(spec,p,t,sample,var):
    block=spec['exclusivity'][f'DVCS_{spec["period_codes"][p]}_{t}'][sample]
    bounds=saved_window(block[var]) if var in block else [None,None]
    cfg=spec['global_cuts']
    if var=='pTmiss' and cfg['enable_pTmiss_cut']:bounds[1]=min(bounds[1] if bounds[1] is not None else float('inf'),cfg['pTmiss_max'])
    return bounds

def view_limits(spec,p,t,vi,d,m):
    v=spec['variables'][vi];edges=np.linspace(v['low'],v['high'],v['bins']+1)
    limits=[quantile(h[1:-1],edges,q) for h in (d,m) for q in (.005,.995)]
    limits.extend(x for s in ('data','mc') for x in effective_window(spec,p,t,s,v['name']) if x is not None)
    limits=[x for x in limits if math.isfinite(x)]
    if not limits:return v['low'],v['high']
    lo,hi=min(limits),max(limits);pad=max(.1*(hi-lo),5*(edges[1]-edges[0]))
    return max(v['low'],lo-pad),min(v['high'],hi+pad)

def plot_overviews(output,arrays,spec,rows,dpi):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    index={(r['period'],r['topology'],r['axis'],r['stratum'],r['variable']):r for r in rows}
    for p in spec['periods']:
        for ti,t in enumerate(TOPOLOGIES):
            destination=output/p/t;destination.mkdir(parents=True)
            fig,axs=plt.subplots(3,3,figsize=(15,11),constrained_layout=True)
            for vi,v in enumerate(spec['variables']):
                ax=axs.flat[vi];edges=np.linspace(v['low'],v['high'],v['bins']+1);hist=[]
                for sample,color in (('data','#205cab'),('mc','#d45c18')):
                    h=arrays[f'{p}__{sample}'][ti,0,vi,:v['bins']+2];hist.append(h)
                    total=float(h.sum())
                    ax.stairs(h[1:-1]/total if total else h[1:-1],edges,color=color,label=f'{sample.upper()} N={int(total):,}')
                    for bound in effective_window(spec,p,t,sample,v['name']):
                        if bound is not None:ax.axvline(bound,color=color,ls='--',alpha=.65,lw=1)
                r=index[p,t,'inclusive','Inclusive',v['name']]
                ax.set(xlim=view_limits(spec,p,t,vi,*hist),xlabel=v['label'],ylabel='Fraction of finite candidates / bin')
                ax.set_title(f"spread ratio={r['width_ratio']:.2f}; median shift={r['median_shift_MC_spreads']:+.2f} MC spreads",fontsize=10)
                ax.legend(fontsize=8);ax.grid(alpha=.15)
                if r['status']!='usable shape comparison':ax.text(.02,.96,r['status'],va='top',transform=ax.transAxes,fontsize=7,color='#a52a2a',wrap=True)
            for ax in axs.flat[7:]:ax.axis('off')
            axs.flat[7].text(0,.95,'Each panel leaves its own cut out.\nOther saved DATA/MC cuts remain.\nDashed lines: omitted cut boundaries.\nUnit-weighted candidates; no subtraction.\nWidths are conditional shape diagnostics.\nMetrics use full histogram; tails outside\nits range are recorded and flagged.',va='top',fontsize=11)
            fig.suptitle(f'{p} · {t.replace("_","–")} · inclusive exclusivity shapes')
            fig.savefig(destination/'inclusive_shapes.png',dpi=dpi);plt.close(fig)
            ng=len(spec['strata']);values=np.full((3,7,ng-1),np.nan)
            for vi,v in enumerate(spec['variables']):
                for g,s in enumerate(spec['strata'][1:]):
                    r=index[p,t,s['axis'],s['label'],v['name']]
                    if r['status']=='usable shape comparison':
                        values[:,vi,g]=[r['width_ratio'],r['median_shift_MC_spreads'],r['excess_tail_fraction']]
            fig,axs=plt.subplots(3,1,figsize=(max(16,(ng-1)*.36),10),constrained_layout=True)
            for ax,z,label,lo,hi in zip(axs,values,['DATA/MC core spread','Median shift / MC core spread','Excess fraction outside MC central 95%'],[.5,-1.5,-.2],[2,1.5,.2]):
                cm=plt.get_cmap('coolwarm').copy();cm.set_bad('#dddddd')
                im=ax.imshow(z,aspect='auto',cmap=cm,vmin=lo,vmax=hi)
                ax.set(yticks=range(7),yticklabels=[v['name'] for v in spec['variables']],xticks=range(ng-1),xticklabels=[s['label'] for s in spec['strata'][1:]],title=label)
                ax.tick_params(axis='x',rotation=90,labelsize=7)
                for ai,name in enumerate(spec['axes']):
                    positions=[g for g,s in enumerate(spec['strata'][1:]) if s['axis']==name]
                    ax.text((positions[0]+positions[-1])/2,-.7,name,ha='center',fontsize=8)
                    ax.axvline(positions[0]-.5,color='white',lw=2)
                fig.colorbar(im,ax=ax,extend='both',pad=.01)
            fig.suptitle(f'{p} · {t} · gray = insufficient statistics/range/bin resolution; color scale saturates at its limits')
            fig.savefig(destination/'shape_summary.png',dpi=dpi);plt.close(fig)

def clean_json(value):
    if isinstance(value,dict):return {k:clean_json(v) for k,v in value.items()}
    if isinstance(value,list):return [clean_json(v) for v in value]
    if isinstance(value,(float,np.floating)) and not math.isfinite(value):return None
    if isinstance(value,np.generic):return value.item()
    return value

def render(output,arrays,manifest,min_entries,dpi):
    spec=manifest['spec'];rows,detail=comparisons(arrays,spec,min_entries)
    with (output/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    plot_overviews(output,arrays,spec,rows,dpi)
    windows={f'{p}|{t}|{v["name"]}':{s:effective_window(spec,p,t,s,v['name']) for s in ('data','mc')}
        for p in spec['periods'] for t in TOPOLOGIES for v in spec['variables']}
    content=json.dumps(clean_json(dict(spec=spec,detail=detail,windows=windows)),allow_nan=False,separators=(',',':')).replace('<','\\u003c')
    template=SOURCE/'exclusivity_resolution_report.html'
    (output/'REPORT.html').write_text(template.read_text().replace('__PAYLOAD__',content))
    flowlines=[]
    for key,flow in manifest['flows'].items():
        flowlines.append(f"| {key.replace('__',' / ')} | {flow['entries']:,} | "+' | '.join(f'{n:,}' for n in flow['selected'])+' |')
    report='''# Exclusivity shape diagnostic

Open **REPORT.html** first. It lets you choose a period, topology, variable,
particle momentum/angle axis and interval, without reading ROOT again.

- `inclusive_shapes.png`: seven inclusive DATA/MC overlays per period/topology.
- `shape_summary.png`: momentum/angle dependence of width, shift and tails.
- `summary.csv`: every numerical comparison, including range/statistics flags.
- `histograms.npz` and `manifest.json`: complete reusable histograms and provenance.

Each distribution applies the saved production global cuts and all saved
sample-specific exclusivity cuts **except the variable displayed**. If that
variable has no cut, all selected cuts remain. For pTmiss its optional global
pre-cut is also omitted. Dashed boundaries are NOT applied to that panel.
These selections can sculpt other correlated variables differently in DATA/MC.
They are conditional diagnostics, not isolated detector-resolution measurements.

Core spread = (q84 − q16)/2, a robust descriptive width, not a Gaussian fit.
Median shift is expressed in MC core spreads. Peak is the histogram mode after
[1,2,1] smoothing; it is descriptive and can be unstable. Tail fraction is the
fraction outside MC's central 95% interval, including histogram under/overflow.
Quantiles use within-bin interpolation of in-range finite values. Comparisons
with <99% histogram coverage, too few entries, nonfinite values, or core spread
less than two histogram bins are flagged; do not tune smearing from those cells.
No uncertainty/significance estimate is attached to width ratios.

Candidates are unit weighted. No pi0/background subtraction, generator
reweighting, current/proton-efficiency weighting or cross-section normalization
is applied. The study uses the saved global phase space, not just populated CSV
bins. Different physics populations, backgrounds and radiation can change these
shapes. A width ratio alone does not identify a particle to smear.
The six momentum/angle subdivisions are separate 1D projections, not a joint
matching of all kinematics. A single event appears in several projections.

The producer stores theta_gamma_gamma in degrees and Emiss2 as missing energy
in GeV. Particle e/p1/p2 polar angles are converted from radians to degrees.
Default data inputs are reconstructed PhysicsEvents parts; MC uses reconstructed
rec_only parts at the period's nominal 45/50 nA, without generated-tree I/O.

Interpretation: first locate consistent shifts/broadening across observables,
then check backgrounds and matched populations. Smearing is a subsequent test;
this script neither changes events nor changes the production analysis.

## Input and cut flow

| Sample | Entries | FD_FD after all cuts | CD_FD after all cuts | CD_FT after all cuts |
|---|---:|---:|---:|---:|
'''+ '\n'.join(flowlines)+'\n'
    (output/'REPORT.md').write_text(report)
    print(f'[done] Open {output/"REPORT.html"}',flush=True)

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--production-run',type=Path,help='Preserved nominal run containing saved sources and output/jsons')
    parser.add_argument('--input-root',type=Path,default=Path('/volatile/clas12/thayward/dvcs_pass2_full_rebuild'))
    parser.add_argument('--periods',nargs='+',choices=PERIODS,default=['fa18_inb','fa18_out'])
    parser.add_argument('--workers',type=int,default=10,help='Total worker budget shared by all samples/periods')
    parser.add_argument('--max-files-per-sample',type=int,default=0,help='0: all; positive: first N sorted parts per period/sample')
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--reuse-histograms',type=Path,help='Existing completed study; no ROOT or production inputs needed')
    parser.add_argument('--config',type=Path,help='Optional JSON histogram ranges/bins and particle-axis edges')
    parser.add_argument('--min-entries',type=int,default=300)
    parser.add_argument('--dpi',type=int,default=120)
    args=parser.parse_args(argv)
    if args.workers<1 or args.workers>10 or args.max_files_per_sample<0 or args.min_entries<1 or args.dpi<50:parser.error('Require workers 1–10, nonnegative file limit, positive minimum entries, dpi >=50')
    args.periods=list(dict.fromkeys(args.periods))
    output=args.output.resolve()
    if output.exists():parser.error(f'Output already exists: {output}; choose a new directory')
    if args.reuse_histograms:
        if args.config or args.production_run:parser.error('Reuse uses its saved cuts, source provenance, ranges and strata; omit --config/--production-run')
        previous=args.reuse_histograms.resolve();manifest=json.loads((previous/'manifest.json').read_text())
        if manifest['spec']['version']!=VERSION:parser.error('Incompatible saved histogram version')
        if digest(previous/'histograms.npz')!=manifest['histogram_sha256']:parser.error('Histogram checksum mismatch')
        with np.load(previous/'histograms.npz',allow_pickle=False) as data:arrays={k:data[k] for k in data.files}
        output.mkdir(parents=True)
        for name in ('histograms.npz','manifest.json'):shutil.copy2(previous/name,output/name)
        write_json(output/'spec.json',manifest['spec'])
        print('[reuse] No ROOT inputs will be read',flush=True)
    else:
        if not args.production_run:parser.error('--production-run is required unless --reuse-histograms is supplied')
        spec,paths=build_spec(args.production_run,args.periods,args.config)
        files=inventory(args.input_root.resolve(),args.periods,args.max_files_per_sample)
        output.mkdir(parents=True)
        arrays,manifest=scan(output,spec,paths,files,args.workers)
    render(output,arrays,manifest,args.min_entries,args.dpi)

if __name__=='__main__':main()
