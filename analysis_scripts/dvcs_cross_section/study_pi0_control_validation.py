#!/usr/bin/env python3
"""Nominal/add-pT control comparison and independent-test coarse weighting."""
from __future__ import annotations
import argparse,csv,html,json
from pathlib import Path
import numpy as np
import study_pi0_control_transfer as control
base=control.base;common=control.common;cuts=control.cuts
T=np.array([.11,.25,.4,.6,1.]);Q=np.array([1.,2.5,4.,12.])


def cell_indices(a):
    ti=np.searchsorted(T,a[:,9],side='right')-1;qi=np.searchsorted(Q,a[:,13],side='right')-1
    valid=(ti>=0)&(ti<len(T)-1)&(qi>=0)&(qi<len(Q)-1)&np.isfinite(a[:,9])&np.isfinite(a[:,13])
    return ti,qi,valid


def fit_weights(data,mc,min_counts=50):
    """Training only; fixed coarse shared coordinates; no topology/photon weights."""
    alpha=len(data)/len(mc) if len(mc) and len(data) else float('nan')
    if not np.isfinite(alpha):raise ValueError('Empty control training sample')
    td,qd,vd=cell_indices(data);tm,qm,vm=cell_indices(mc);weights=np.ones((len(T)-1,len(Q)-1));rows=[]
    for i in range(len(T)-1):
        ndt=int((vd&(td==i)).sum());nmt=int((vm&(tm==i)).sum());fallback=ndt/(alpha*nmt) if min(ndt,nmt)>=min_counts else 1.
        for j in range(len(Q)-1):
            nd=int((vd&(td==i)&(qd==j)).sum());nm=int((vm&(tm==i)&(qm==j)).sum());supported=min(nd,nm)>=min_counts
            weights[i,j]=nd/(alpha*nm) if supported else fallback
            rows.append(dict(t_low=T[i],t_high=T[i+1],Q2_low=Q[j],Q2_high=Q[j+1],train_data=nd,train_mc=nm,weight=weights[i,j],status='cell trained' if supported else 't-only fallback' if min(ndt,nmt)>=min_counts else 'unity fallback'))
    return alpha,weights,rows


def weights_for(a,weights):
    ti,qi,valid=cell_indices(a);result=np.ones(len(a));result[valid]=weights[ti[valid],qi[valid]];return result


def train_models(data,nominal,modified,min_counts):
    train=lambda a:a[a[:,11]<.5]
    alpha,weights,rows=fit_weights(train(data),train(modified),min_counts)
    if not len(train(nominal)):raise ValueError('Empty nominal training control')
    scales=dict(nominal=len(train(data))/len(train(nominal)),add_pt=alpha,add_pt_weighted=len(train(data))/weights_for(train(modified),weights).sum())
    return weights,rows,scales


def evaluate(data,cases,scales,weights,period,out):
    import matplotlib.pyplot as plt
    result=[];links=[];totals=[]
    for c,a in cases.items():
        w=weights_for(a,weights) if c=='add_pt_weighted' else np.ones(len(a));totals.append(dict(period=period,scenario=c,test_data=len(data),test_raw_mc=len(a),test_prediction=scales[c]*w.sum(),data_over_prediction=len(data)/(scales[c]*w.sum()) if w.sum() else float('nan')))
    for top,label in [(None,'all_control_topologies'),*enumerate(base.TOPS)]:
        d=data if top is None else data[data[:,0]==top];fig,axes=plt.subplots(3,3,figsize=(14,10),layout='constrained')
        for ai,(axis,(col,edges)) in enumerate(control.EDGES.items()):
            factor=180/np.pi if axis.endswith('_deg') else 1.;edges=np.array(edges);dh=np.histogram(d[:,col]*factor,edges)[0];ax=axes.flat[ai]
            for c,color in [('nominal','#888888'),('add_pt','#22946b'),('add_pt_weighted','#a343a8')]:
                a=cases[c];a=a if top is None else a[a[:,0]==top];w=weights_for(a,weights) if c=='add_pt_weighted' else np.ones(len(a));mh=np.histogram(a[:,col]*factor,edges,weights=w)[0]*scales[c];mv=np.histogram(a[:,col]*factor,edges,weights=w*w)[0]*scales[c]**2
                with np.errstate(divide='ignore',invalid='ignore'):
                    r=dh/mh;err=np.sqrt(dh/mh**2+dh**2*mv/mh**4)
                r[mh<=0]=np.nan;err[mh<=0]=np.nan;ax.errorbar((edges[:-1]+edges[1:])/2,r,yerr=err,fmt='o-',ms=3,lw=1,color=color,label=c)
                for i,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):result.append(dict(period=period,topology=label,axis=axis,low=lo,high=hi,scenario=c,test_data=int(dh[i]),test_prediction=mh[i],data_over_prediction=r[i],conditional_stat_error=err[i],status='empty prediction' if mh[i]<=0 else 'low DATA' if dh[i]<50 else 'conditional comparison'))
            ax.axhline(1,color='black',lw=.7);ax.set(xlabel=axis,ylabel='Held-out DATA / predicted MC');ax.grid(alpha=.2);ax.legend(fontsize=7)
        for ax in axes.flat[len(control.EDGES):]:ax.axis('off')
        fig.suptitle(period+' · '+label+' · independent test halves; training scales held fixed')
        folder=out/'plots'/period;folder.mkdir(parents=True,exist_ok=True);path=folder/(label+'.png');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
    return result,totals,links


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--control-study',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=10);p.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'));p.add_argument('--min-train-count',type=int,default=50)
    args=p.parse_args(argv);origin=args.control_study.resolve();out=args.output.resolve()
    if out.exists():p.error('Output exists; choose a new directory')
    if not 1<=args.workers<=10 or args.min_train_count<1:p.error('Workers must be 1–10; minimum count positive')
    manifest=json.loads((origin/'manifest.json').read_text());export=json.loads((origin/'export_manifest.json').read_text());fp=export['fingerprint'];cmroot=Path(manifest['arguments']['common_study']);cm=json.loads((cmroot/'manifest.json').read_text());source=Path(cm['arguments']['followup_study'])/'cut_redetermination';old=json.loads((source/'manifest.json').read_text());spec=json.loads((source/'spec.json').read_text());periods=list(fp['candidate_hashes'])
    if base.digest(cmroot/'manifest.json')!=fp['common_sha256'] or base.digest(source/'spec.json')!=fp['spec_sha256']:p.error('Input provenance changed')
    for n,h in fp['worker_sources'].items():
        if base.digest(base.SOURCE/n)!=h:p.error('Worker changed since control export: '+n)
    exports={period:{} for period in periods}
    for r in export['records']:
        path=Path(r['path']);shape=r['shape']
        if common.export_digest(path)!=r['sha256'] or len(shape)!=2 or shape[1]!=20 or path.stat().st_size!=int(np.prod(shape))*8:p.error('Invalid prior control export')
        exports[r['period']][r['sample']]=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape)) if shape[0] else np.empty((0,20))
    if any(set(v)!=set(('pi0_data','pi0_mc')) for v in exports.values()):p.error('Missing control exports')
    ref=Path(old['arguments']['reference_study']);run=Path(old['arguments']['response_run'])
    if (run/'with_bias/cache_manifest.json').is_file():run=run/'with_bias'
    if base.digest(run/'cache_manifest.json')!=old['response_cache_sha256']:p.error('MC cache metadata changed')
    for n,h in fp['saved_globals'].items():
        if base.digest(ref/'backend'/n)!=h:p.error('Saved globals changed')
    cache=json.loads((run/'cache_manifest.json').read_text());jobs=[]
    for period in periods:
        r=cache['samples'][period+'__pi0_mc'];path=(run/r['file']).resolve()
        if common.export_digest(path)!=r['sha256']:p.error('MC cache changed')
        jobs.append(dict(period=period,sample='pi0_mc',cache=str(path),parameters=[0.]*20,export_kinematics=True,export_global_only=True))
    out.mkdir(parents=True);work=base.scratch_directory(args.scratch_root,out);base.dump(out/'scratch_location.json',dict(directory=str(work)));base.dump(out/'spec.json',spec)
    print('[reuse] Existing DATA and corrected-control exports; only nominal MC reevaluated from cache',flush=True)
    worker=base.compile_worker(out,ref);pool=base.Pool(worker,out/'spec.json',work,args.workers);records=[]
    for _,job,path in pool.run('cut-arrays',jobs,'nominal_control'):
        shape=json.loads(Path(str(path)+'.json').read_text())['shape']
        if len(shape)!=2 or shape[1]!=20 or path.stat().st_size!=int(np.prod(shape))*8:raise ValueError('Invalid nominal export')
        exports[job['period']]['nominal']=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape)) if shape[0] else np.empty((0,20));records.append(dict(period=job['period'],path=str(path),shape=shape,sha256=common.export_digest(path)))
    rows=[];totals=[];weightrows=[];links=[]
    for period in periods:
        a=exports[period];data=control.selected_control(a['pi0_data'],period,spec);nom=control.selected_control(a['nominal'],period,spec);mod=control.selected_control(a['pi0_mc'],period,spec)
        test=lambda x:x[x[:,11]>=.5]
        weights,wr,scales=train_models(data,nom,mod,args.min_train_count)
        for r in wr:r['period']=period
        weightrows+=wr
        rr,tt,ll=evaluate(test(data),dict(nominal=test(nom),add_pt=test(mod),add_pt_weighted=test(mod)),scales,weights,period,out);rows+=rr;totals+=tt;links+=ll
        print('[done] '+period+' train-only weights and independent test comparison',flush=True)
    for name,r in [('heldout_ratios',rows),('heldout_totals',totals),('training_weights',weightrows)]:base.table(out/(name+'.csv'),r)
    base.dump(out/'manifest.json',dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},input_manifest_sha256=base.digest(origin/'manifest.json'),export_manifest_sha256=base.digest(origin/'export_manifest.json'),nominal_exports=records,training_half=0,test_half=1,t_edges=T.tolist(),Q2_edges=Q.tolist(),weighting='all control topologies together, separately per period; no topology or detected-photon weights'))
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Control response and weighting validation</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>Independent control validation</h1><p>All plots use only test-half events. Nominal and add-pT MC each receive a total scale determined on training-half DATA/MC. The purple prediction additionally uses training-only coarse |t|×Q² weights shared across control topologies. No test yields are used in scales, cell edges, fallback rules or weights. DATA/MC=1 is not imposed in the test sample.</p><p>Weights require the requested minimum counts in both training samples, otherwise fall back to the |t| marginal and finally unity. Unsupported coordinates receive unity. The weighted total scale is also determined only on training events. The grid is fixed in advance. Check heldout_totals.csv, then all_control_topologies and FD_FD plots. Other coordinates/topologies test behavior beyond the weighted projections.</p><p>Error bars are conditional on fitted training scales/weights, omit training/response/model uncertainty and are correlated across curves. Do not interpret empty or zero-DATA bins as precise exclusions. Better control closure is a diagnostic, not proof of correct missing-second-photon efficiency. Controls can contain other backgrounds. Reconstructed |t| and Q² weights can compensate response or generator differences; they are not a Born correction.</p><p>No leakage weights or production acceptance/cross-section changes are applied. A topology label in the two-photon control must not be copied to the one-photon leakage. No photon-momentum weighting is used. Branch/source provenance remains that of the saved processing inputs; this script does not inspect the development tree as processing evidence.</p><p>'+ ' · '.join('<a href="'+n+'.csv">'+n+'</a>' for n in ('heldout_totals','heldout_ratios','training_weights'))+'</p><ul>'+''.join('<li><a href="'+l+'">'+html.escape(l)+'</a></li>' for l in links)+'</ul>')
    print('[done] '+str(out/'REPORT.html')+'; no DATA ROOT read or response fit',flush=True)

if __name__=='__main__':main()
