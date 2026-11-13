#!/usr/bin/env python3
"""Differential two-photon control closure; fixed add-pT response and period scale."""
from __future__ import annotations
import os
for _n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_n]='1'
import argparse,csv,html,json
from pathlib import Path
import numpy as np
import study_fd_background_transfer as fd
common=fd.common;base=fd.base;cuts=fd.cuts;closure=fd.closure
AXES={k:v for k,v in fd.AXES.items() if k not in ('photon_p','photon_theta_deg')}
# In a two-photon tree p2 is the composite pi0. Never equate it to the leakage photon.
EDGES={'t_abs':(9,[.11,.15,.25,.4,.6,.8,1.]),'xB':(12,[0,.2,.3,.4,.5,1.]),'Q2':(13,[1,2,3,4,6,12]),'electron_p':(14,[0,4,5,6,7,8,12]),'electron_theta_deg':(15,[0,10,12,15,20,30,180]),'proton_p':(16,[0,.42,.5,.7,1,1.5,3,10]),'proton_theta_deg':(17,[0,15,25,35,45,65,180])}


def control_ratio(d,m,scale,prob):
    if m<=0 or not np.isfinite(scale) or scale<=0:return float('nan'),float('nan')
    ratio=d/(scale*m/prob)
    # Conditional on the period normalization; bins and that anchor are correlated.
    error=prob/scale*np.sqrt(d/m**2+d*d/m**3)
    return float(ratio),float(error)


def selected_control(a,period,spec):
    cfg=next(p for p in cuts.production.PERIODS if p.key==period);keep=np.zeros(len(a),bool)
    for top,topology in enumerate(cuts.production.TOPOLOGIES):
        block=spec['exclusivity'][cuts.production.cut_json_key('eppi0',cfg,topology)]
        keep|=closure.mask(a,top,block,'data')
    return a[keep]


def analyze_controls(exports,spec,normprob,periods,out):
    import matplotlib.pyplot as plt
    rows=[];links=[];closure_rows=[]
    for period in periods:
        data=selected_control(exports[period]['pi0_data'],period,spec);mc=selected_control(exports[period]['pi0_mc'],period,spec)
        norm,prob=normprob[period];scale=norm['control_scale'];prob=prob['pi0_mc'];ratio,_=control_ratio(len(data),len(mc),scale,prob)
        closure_rows.append(dict(period=period,data=len(data),raw_mc=len(mc),predicted_mc=scale*len(mc)/prob,ratio=ratio,status='anchor reproduced' if np.isclose(ratio,1,rtol=.01) else 'ANCHOR MISMATCH: investigate before interpreting subdivisions'))
        folder=out/'plots'/period;folder.mkdir(parents=True)
        for top,label in [(None,'all_control_topologies'),*enumerate(base.TOPS)]:
            d=data if top is None else data[data[:,0]==top];m=mc if top is None else mc[mc[:,0]==top]
            fig,axes=plt.subplots(3,3,figsize=(14,10),layout='constrained')
            for ai,(axis,(col,edges)) in enumerate(EDGES.items()):
                edges=np.array(edges,dtype=float);factor=180/np.pi if axis.endswith('_deg') else 1.
                dv=d[:,col]*factor;mv=m[:,col]*factor;ax=axes.flat[ai];xx=[];yy=[];ee=[]
                for lo,hi in zip(edges[:-1],edges[1:]):
                    nd=int(((dv>=lo)&(dv<hi)).sum());nm=int(((mv>=lo)&(mv<hi)).sum());r,e=control_ratio(nd,nm,scale,prob)
                    rows.append(dict(period=period,topology=label,axis=axis,low=lo,high=hi,data=nd,raw_mc=nm,predicted_mc=scale*nm/prob,data_over_prediction=r,conditional_stat_error=e,status='empty MC' if not nm else 'low counts' if min(nd,nm)<50 else 'counting comparison; control backgrounds not subtracted'))
                    xx.append((lo+hi)/2);yy.append(r);ee.append(e)
                ax.errorbar(xx,yy,yerr=ee,fmt='o-',lw=1);ax.axhline(1,color='black',lw=.8);ax.set(xlabel=axis,ylabel='Control DATA / fixed-scale MC');ax.grid(alpha=.2)
            for ax in axes.flat[len(EDGES):]:ax.axis('off')
            fig.suptitle(period+' · '+label+' · two-photon control, common saved DATA windows')
            path=folder/(label+'.png');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
        print('[done] '+period+' differential two-photon controls',flush=True)
    base.table(out/'control_ratios.csv',rows);base.table(out/'anchor_check.csv',closure_rows)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Differential pi0 control</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>Differential two-photon control closure</h1><p>First check anchor_check.csv: the all-topology total should reproduce the old selected control anchor. Each plot then divides DATA by the MC prediction using that one fixed period scale. No per-bin renormalization occurs. Unlike the earlier cut-reference plots, each point is a disjoint kinematic bin; one is not guaranteed.</p><p>All samples use the saved full globals and common saved DATA control windows. MC use the frozen add-pT response. Error bars are conditional counting errors with the scale held fixed; bins share normalization/model uncertainty. These are reconstructed control yields, not generated cross sections or measured transfer corrections. Control DATA can contain other backgrounds.</p><p>Leading/soft photon thresholds are inherited from the input trees. This script does not impose or verify a 0.35/0.4 GeV second-photon threshold and cannot recover events already removed upstream. Composite pi0 momentum is not treated as a single detected photon. No new leakage weights or production corrections are applied. First validate the differential control, exposures and input selections.</p><p><a href="anchor_check.csv">Total-anchor check</a> · <a href="control_ratios.csv">Differential counts/ratios</a></p><ul>'+''.join('<li><a href="'+p+'">'+html.escape(p)+'</a></li>' for p in links)+'</ul>')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--common-study',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--periods',nargs='+');p.add_argument('--workers',type=int,default=10);p.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'));p.add_argument('--reuse-exports',type=Path)
    args=p.parse_args(argv)
    if not 1<=args.workers<=10:p.error('Workers must be 1–10 total')
    origin=args.common_study.resolve();out=args.output.resolve()
    if out.exists():p.error('Output exists; choose a new directory')
    cm=json.loads((origin/'manifest.json').read_text());source=Path(cm['arguments']['followup_study'])/'cut_redetermination';closed=Path(cm['arguments']['closure_study']);old=json.loads((source/'manifest.json').read_text());spec=json.loads((source/'spec.json').read_text());periods=list(dict.fromkeys(args.periods or spec['periods']))
    if not set(periods)<=set(spec['periods']):p.error('Requested period absent from source')
    if base.digest(source/'manifest.json')!=cm['source_sha256'] or base.digest(closed/'manifest.json')!=cm['closure_sha256'] or old['candidate_hashes']!=cm['candidate_hashes']:p.error('Common-study provenance differs from source/cuts')
    ref=Path(old['arguments']['reference_study']);run=Path(old['arguments']['response_run']);training=Path(old['arguments']['training_study'])
    if (run/'with_bias/cache_manifest.json').is_file():run=run/'with_bias'
    cache=json.loads((run/'cache_manifest.json').read_text())
    if base.digest(run/'cache_manifest.json')!=old['response_cache_sha256']:p.error('Response cache metadata changed')
    paths=cuts.candidate_paths({'add_pt':run},periods,training);vectors={}
    for (period,c),path in paths.items():
        if base.digest(path)!=cm['candidate_hashes'][c][period]:p.error('Candidate parameters changed')
        vectors[period]=cuts.load_candidate_parameters(path)
    norms={}
    for r in csv.DictReader((origin/'normalizations.csv').open()):
        if r['scenario']=='add_pt' and r['topology']=='FD_FD':norms[r['period']]={k:float(v) for k,v in r.items() if k not in ('period','scenario','topology','status')}
    if not set(periods)<=norms.keys():p.error('Common-study FD–FD normalization missing')
    fingerprint=dict(study='pi0-control-v1',common_sha256=base.digest(origin/'manifest.json'),normalizations_sha256=base.digest(origin/'normalizations.csv'),candidate_hashes={period:cm['candidate_hashes']['add_pt'][period] for period in periods},spec_sha256=base.digest(source/'spec.json'),worker_sources={n:base.digest(base.SOURCE/n) for n in ('particle_response_cut_arrays.h','particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h','exclusivity_resolution_core.h','exclusivity_resolution_worker.cpp')},saved_globals={n:base.digest(ref/'backend'/n) for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h')})
    exports={period:{} for period in periods};records=[]
    if args.reuse_exports:
        saved=json.loads((args.reuse_exports/'export_manifest.json').read_text())
        if saved['fingerprint']!=fingerprint:p.error('Reused exports have different inputs/worker')
        records=saved['records']
        for r in records:
            if common.export_digest(r['path'])!=r['sha256']:p.error('Reused export changed')
    else:
        # Validate every source before creating report/scratch or starting a ROOT scan.
        jobs=[];meta=json.loads((ref/'manifest.json').read_text())
        if base.digest(ref/'manifest.json')!=old['reference_sha256']:p.error('Reference inventory changed')
        for period in periods:
            entries=[r for r in meta['inventory'] if r['period']==period and r['sample']=='pi0_data']
            if not entries:p.error('Missing DATA inventory')
            for r in entries:
                st=Path(r['path']).stat()
                if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):p.error('DATA input changed')
            common_job=dict(period=period,export_kinematics=True,export_global_only=True)
            jobs.append(dict(**common_job,sample='pi0_data',files=[r['path'] for r in entries]))
            for sample in ('pi0_mc',):
                r=cache['samples'][period+'__'+sample];path=(run/r['file']).resolve()
                if common.export_digest(path)!=r['sha256']:p.error('MC cache changed')
                jobs.append(dict(**common_job,sample=sample,cache=str(path),parameters=vectors[period]))
        out.mkdir(parents=True);work=base.scratch_directory(args.scratch_root,out);base.dump(out/'scratch_location.json',dict(directory=str(work)));base.dump(out/'spec.json',spec)
        print('[export] Two-photon DATA once per period; cached morphed two-photon MC; all control topologies, no fit',flush=True)
        worker=base.compile_worker(out,ref);pool=base.Pool(worker,out/'spec.json',work,args.workers)
        for i,job,path in pool.run('cut-arrays',jobs,'control_kinematics'):
            shape=json.loads(Path(str(path)+'.json').read_text())['shape']
            records.append(dict(period=job['period'],sample=job['sample'],path=str(path),shape=shape,sha256=common.export_digest(path)))
        for r in meta['inventory']:
            if r['period'] in periods and r['sample']=='pi0_data':
                st=Path(r['path']).stat()
                if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):raise ValueError('DATA input changed during export')
    for r in records:
        shape=r['shape'];path=Path(r['path'])
        if len(shape)!=2 or shape[1]!=20 or path.stat().st_size!=int(np.prod(shape))*8:p.error('Invalid extended export layout')
        exports[r['period']][r['sample']]=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape)) if shape[0] else np.empty((0,20))
    if any(set(samples)!=set(('pi0_data','pi0_mc')) for samples in exports.values()):p.error('Incomplete export samples')
    out.mkdir(parents=True,exist_ok=True);base.dump(out/'export_manifest.json',dict(fingerprint=fingerprint,records=records));base.dump(out/'manifest.json',dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},fingerprint=fingerprint))
    normprob={period:(norms[period],{s:cache['samples'][period+'__'+s]['probability'] for s in ('pi0_mc',)}) for period in periods}
    analyze_controls(exports,spec,normprob,periods,out)
    print('[done] '+str(out/'REPORT.html')+'; transfer report directory only; extended arrays remain on scratch',flush=True)

if __name__=='__main__':main()
