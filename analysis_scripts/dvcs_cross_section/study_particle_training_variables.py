#!/usr/bin/env python3
"""Compare training-observable choices using a completed with-bias response cache."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,copy,csv,json
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
import study_particle_response as base

SCORED=(0,2,3,4,5,6)
CONFIGURATIONS={
    'current':(0,4,5),
    'add_pt':(0,3,4,5),
    'add_pt_angle':(0,2,3,4,5),
}
LABELS={'current':'Original training variables (refined)','add_pt':'+ missing transverse momentum','add_pt_angle':'+ transverse momentum and opening angle'}

def parse_weights(items):
    weights={base.OBS[v]:1. for v in SCORED};seen=set()
    for item in items:
        try:name,value=item.split('=',1);value=float(value)
        except ValueError:raise ValueError('Weights must be VARIABLE=positive_number')
        if name not in weights or name in seen or not np.isfinite(value) or value<=0:raise ValueError('Invalid or repeated observable weight: '+item)
        weights[name]=value;seen.add(name)
    return weights

class TrainingEngine(base.Engine):
    def __init__(self,*args,training,weights,**kwargs):
        super().__init__(*args,**kwargs)
        self.training=tuple(training);self.weights=weights
        if not set(self.training)<=set(SCORED) or len(set(self.training))!=len(self.training) or not self.training or 6 in self.training:raise ValueError('Training must be unique scored observables, reserving Mx2_2')
    def score_one(self,hist,i,multiplier,reference=False):
        _,records=super().score_one(hist,i,multiplier,reference)
        buckets={k:[] for k in ('train','validation_mc','heldout_observables')}
        for row in records:
            v=base.OBS.index(row['variable']);kind='heldout_observables' if v not in self.training else 'validation_mc' if row['half'] else 'train'
            row['kind']=kind;row['observable_weight']=self.weights[row['variable']]
            buckets[kind].append((row['score'],row['observable_weight']))
        result={}
        for kind,values in buckets.items():
            expected=len(self.fit_topologies)*(len(self.training) if kind!='heldout_observables' else 2*(len(SCORED)-len(self.training)))
            result[kind]=float(np.average([x for x,w in values],weights=[w for x,w in values])) if len(values)==expected and all(np.isfinite(x) for x,w in values) else float('inf')
        result['training_projections']=len(buckets['train'])
        return result,records

def starting_parameters(path,period):
    value=json.loads(path.read_text())
    if value.get('period')!=period or value.get('bias_mode')!='fit':raise ValueError('Use completed with-bias candidates, one per period')
    if any(name not in value['parameters'] for name in base.PARAMS):raise ValueError('Candidate does not cover all detector responses: '+str(path))
    parameters=[float(value['parameters'][name]) for name in base.PARAMS]
    if not np.all(np.isfinite(parameters)):raise ValueError('Nonfinite starting response')
    return parameters,value

def nearest_trials(grid,best,k):
    levels=sorted(set(grid[k]+[best[k]]))
    lower=max([x for x in levels if x<best[k]],default=levels[0]);upper=min([x for x in levels if x>best[k]],default=levels[-1])
    values=sorted(set([best[k],lower,(lower+best[k])/2,(best[k]+upper)/2,upper]))
    return values

def fit_local(engine,start,args,folder,reference):
    """Use training only for local refinement; optionally keep a reference fixed."""
    best=list(start);history=[];grid=base.parameter_grid(base.TOPS,'fit',args.scan_range)
    for k,v in enumerate(best):
        if not grid[k][0]<=v<=grid[k][-1]:raise ValueError('Starting fit lies outside search bounds; select --scan-range wide')
    if not reference:
        levels={k:list(v) for k,v in grid.items()}
        for iteration in range(args.rounds):
            for k in grid:
                values=nearest_trials(levels,best,k);trials=[list(best)]
                for value in values:
                    point=list(best);point[k]=value;trials.append(point)
                pars,scores,records,hist,selected=engine.evaluate(trials,f'local_{iteration}_{k}')
                # Include current point, so no local update can worsen its training objective.
                eligible=[i for i,p in enumerate(pars) if tuple(p) in {tuple(x) for x in trials}]
                chosen=min(eligible,key=lambda i:(scores[i]['train'],sum(abs(x-y) for x,y in zip(pars[i],best))))
                if not np.isfinite(scores[chosen]['train']):raise ValueError('No finite training objective; inspect MC statistics and background residuals')
                best=list(pars[chosen]);levels[k]=sorted(set(levels[k]+values))
                history.append(dict(round=iteration,parameter=base.PARAMS[k],chosen_value=best[k],**scores[chosen]))
                base.discard_histograms(hist)
    pars,scores,records,hist,selected=engine.evaluate([start,best],'final')
    start_i=pars.index(tuple(start));best_i=pars.index(tuple(best))
    try:
        p=engine.spec['periods'][0];scenario_rows=[];curves={};selection=[]
        for tag,index in [('starting',start_i),('candidate',best_i)]:
            score=scores[index];multiplier=score['background_multiplier']
            current=[dict(r,stage=tag) for r in records if r['scenario']==index]
            scenario_rows+=current
            for top in range(3):
                for v in range(len(base.OBS)):
                    curves[tag,top,v]=sum(engine.weighted(hist,p,index,half,top,0,v,multiplier)[0] for half in (0,1))/2
                raw0=selected[p+'__mc'][start_i,:,top].sum();raw=selected[p+'__mc'][index,:,top].sum()
                selection.append(dict(stage=tag,topology=base.TOPS[top],raw_signal_selection_ratio=float(raw/raw0) if raw0 else float('nan')))
        base.table(folder/'projection_scores.csv',scenario_rows);base.table(folder/'selection_changes.csv',selection)
        if history:base.table(folder/'fit_steps.csv',history)
        candidate=dict(period=p,bias_mode='fit',parameters=dict(zip(base.PARAMS,best)),starting_candidate_sha256=engine.spec['starting_candidate_sha256'],training_observables=[base.OBS[v] for v in engine.training],reserved_observables=[base.OBS[v] for v in SCORED if v not in engine.training],observable_weights=engine.weights,background_multiplier=scores[best_i]['background_multiplier'],fixed_reference=reference,scan_range=args.scan_range,rounds=0 if reference else args.rounds)
        base.dump(folder/'candidate_parameters.json',candidate)
        per_variable={}
        for tag in ['starting','candidate']:
            for v in SCORED:
                values=[r['score'] for r in scenario_rows if r['stage']==tag and r['variable']==base.OBS[v]]
                per_variable[tag,v]=float(np.mean(values)) if len(values)==6 and all(np.isfinite(x) for x in values) else float('inf')
        return dict(period=p,parameters=best,starting=start,curves=curves,scores=scores[best_i],starting_scores=scores[start_i],per_variable=per_variable,candidate=candidate,records=scenario_rows)
    finally:base.discard_histograms(hist)

def save_checkpoint(folder,result):
    arrays={'__'.join(map(str,key)):value for key,value in result['curves'].items()}
    np.savez_compressed(folder/'comparison_curves.npz',**arrays)
    payload={k:v for k,v in result.items() if k not in ('curves','records','per_variable')}
    payload['per_variable']=[[tag,v,value] for (tag,v),value in result['per_variable'].items()]
    payload['curves_sha256']=base.digest(folder/'comparison_curves.npz');payload['candidate_sha256']=base.digest(folder/'candidate_parameters.json')
    base.dump(folder/'completed_fit.json',payload)

def load_checkpoint(folder):
    meta=json.loads((folder/'completed_fit.json').read_text())
    if base.digest(folder/'comparison_curves.npz')!=meta['curves_sha256'] or base.digest(folder/'candidate_parameters.json')!=meta['candidate_sha256']:raise ValueError('Completed fit checksum mismatch: '+str(folder))
    with np.load(folder/'comparison_curves.npz',allow_pickle=False) as archive:meta['curves']={(key.split('__')[0],int(key.split('__')[1]),int(key.split('__')[2])):archive[key].copy() for key in archive.files}
    meta['per_variable']={(tag,v):float(value) if value is not None else float('inf') for tag,v,value in meta['per_variable']}
    for name in ['scores','starting_scores']:meta[name]={k:v if v is not None else float('inf') for k,v in meta[name].items()}
    return meta


def report(output,spec,data,results,weights):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summary=[];parameters=[];observables=[];plot_links=[];pnames=spec['periods']
    for (period,config),result in results.items():
        fixed=results[period,'current'];scores=result['scores']
        common=lambda z:float(np.average([z['per_variable']['candidate',v] for v in SCORED],weights=[weights[base.OBS[v]] for v in SCORED]))
        summary.append(dict(period=period,configuration=config,training_objective=scores['train'],training_shape=scores['train_shape'],validation_mc=scores['validation_mc'],reserved_Mx2_2_score=result['per_variable']['candidate',6],all_observables_common_score=common(result),starting_common_score=float(np.average([fixed['per_variable']['starting',v] for v in SCORED],weights=[weights[base.OBS[v]] for v in SCORED])),refined_original_training_common_score=common(fixed),background_multiplier=scores['background_multiplier'],parameters_at_scan_limit=sum(v!=0 and v in (base.parameter_grid(base.TOPS,'fit',result['candidate']['scan_range'])[k][0],base.parameter_grid(base.TOPS,'fit',result['candidate']['scan_range'])[k][-1]) for k,v in enumerate(result['parameters']))))
        for k,name in enumerate(base.PARAMS):parameters.append(dict(period=period,configuration=config,parameter=name,starting_value=result['starting'][k],candidate_value=result['parameters'][k],change=result['parameters'][k]-result['starting'][k],training_choice_change=result['parameters'][k]-fixed['parameters'][k],display_training_choice_change=(result['parameters'][k]-fixed['parameters'][k])*(100 if name.endswith('fraction') else 1000),display_change=(result['parameters'][k]-result['starting'][k])*(100 if name.endswith('fraction') else 1000),display_unit='percentage points' if name.endswith('fraction') else 'mrad'))
        for v in SCORED:
            before=fixed['per_variable']['starting',v];after=result['per_variable']['candidate',v]
            observables.append(dict(period=period,configuration=config,variable=base.OBS[v],used_for_training=v in CONFIGURATIONS[config],weight=weights[base.OBS[v]],starting_score=before,candidate_score=after,improvement_percent=100*(1-after/before) if before>0 else float('nan'),controlled_improvement_percent=100*(1-after/fixed['per_variable']['candidate',v]) if fixed['per_variable']['candidate',v]>0 else float('nan')))
    base.table(output/'configuration_comparison.csv',summary);base.table(output/'parameter_stability.csv',parameters);base.table(output/'observable_comparison.csv',observables)
    for period in pnames:
        for top,topology in enumerate(base.TOPS):
            fig,axes=plt.subplots(3,3,figsize=(15,11),layout='constrained')
            for v,variable in enumerate(spec['variables']):
                ax=axes.flat[v];edges=base.diag.variable_edges(variable);n=variable['bins']+2
                d=data[period+'__data'][top,0,v,:n].astype(float);curves=[('DATA',d,'#205cab'),('Saved starting response',results[period,'current']['curves']['starting',top,v],'#d45c18')]
                for config,color in [('current','#666666'),('add_pt','#22946b'),('add_pt_angle','#9467bd')]:curves.append((LABELS[config],results[period,config]['curves']['candidate',top,v],color))
                for label,values,color in curves:
                    if values.sum()>0:ax.stairs(values[1:-1]/values.sum()/np.diff(edges),edges,label=label,color=color,lw=1.2)
                if variable['name']=='Mx2':ax.set_xlim(*base.mass_plot_limits(edges,[x[1] for x in curves]))
                else:
                    q=[base.diag.quantile(values[1:-1],edges,x) for _,values,_ in curves for x in (.01,.99)];q=[x for x in q if np.isfinite(x)]
                    if q:
                        low,high=min(q),max(q);pad=max((high-low)*.08,np.min(np.diff(edges))*2);ax.set_xlim(max(edges[0],low-pad),min(edges[-1],high+pad))
                ax.set(xlabel=variable['label'],ylabel='Unit-area density');ax.grid(alpha=.15);ax.legend(fontsize=7)
            for ax in axes.flat[7:]:ax.axis('off')
            axes.flat[7].text(0,1,'Same cached MC and frozen DATA selection.\nOnly training-observable choices change.\nMissing mass squared eγ stays reserved.\nChecks use the same DATA events;\nno independent DATA validation.',va='top')
            fig.suptitle(period+' · '+topology+' · training-variable comparison')
            path=output/'plots'/(period+'_'+topology+'.png');path.parent.mkdir(exist_ok=True);fig.savefig(path,dpi=110);plt.close(fig);plot_links.append(path.relative_to(output).as_posix())
    values=np.array([[next(r['controlled_improvement_percent'] for r in observables if r['period']==p and r['configuration']==c and r['variable']==base.OBS[v]) for v in SCORED] for p in pnames for c in ['add_pt','add_pt_angle']])
    fig,ax=plt.subplots(figsize=(12,7),layout='constrained');im=ax.imshow(values,cmap='RdYlGn',vmin=-50,vmax=50,aspect='auto')
    ax.set(xticks=range(len(SCORED)),xticklabels=[base.OBS[v] for v in SCORED],yticks=range(len(values)),yticklabels=[p+' / '+c for p in pnames for c in ['add_pt','add_pt_angle']])
    for (i,j),v in np.ndenumerate(values):ax.text(j,i,f'{v:+.0f}%' if np.isfinite(v) else 'unusable',ha='center',va='center',fontsize=9)
    ax.set_title('Training-choice improvement over the equally refined original objective\nPositive = better; same observable, topology average and MC precision')
    fig.colorbar(im,ax=ax,label='Improvement (%)',shrink=.7);fig.savefig(output/'observable_improvements.png',dpi=140);plt.close(fig)
    (output/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Training-variable comparison</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}li{margin:8px}</style><h1>Training-variable comparison</h1><p>All three training choices start from the same saved broadening-plus-bias parameters and receive the same local refinement budget. The saved starting response is retained as a fixed reference. Alternatives add pTmiss, or pTmiss plus photon opening angle. Improvements relative to the refined original objective isolate training-choice effects from simply doing more optimization. Every period/polarity remains independent. The same MC caches, DATA histograms, selection windows, model bounds and background prior are used. No ROOT events are reread and production cuts are unchanged.</p><p>Default observable weights are one. Aggregate training scores with different observable sets are not directly comparable. Use the common per-observable scores and the reserved Mx2_2 check. Reserved observables share DATA events and correlated kinematics with training; this is not independent DATA validation. No configuration is automatically adopted or selected using reserved scores.</p><p><a href="observable_improvements.png">Improvements and regressions by observable</a> · <a href="configuration_comparison.csv">Common-score and reserved-score comparisons</a> · <a href="parameter_stability.csv">Parameter stability</a> · <a href="observable_comparison.csv">Per-observable comparisons</a></p><ul>'+''.join('<li><a href="'+x+'">'+Path(x).stem+'</a></li>' for x in plot_links)+'</ul>')

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-study',type=Path,required=True);parser.add_argument('--response-run',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--resume',action='store_true',help='Reuse completed period/configuration fits in this output; retry incomplete fits')
    parser.add_argument('--workers',type=int,default=10);parser.add_argument('--rounds',type=int,default=2);parser.add_argument('--min-mc',type=int,default=100)
    parser.add_argument('--scan-range',choices=['standard','wide'],default='wide');parser.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'))
    parser.add_argument('--observable-weight',action='append',default=[],metavar='VARIABLE=WEIGHT');parser.add_argument('--background-transfer-uncertainty',type=float,default=.10)
    args=parser.parse_args(argv)
    if not 1<=args.workers<=10 or args.rounds<1 or args.min_mc<10 or not np.isfinite(args.background_transfer_uncertainty) or args.background_transfer_uncertainty<0:parser.error('Invalid worker/fit settings')
    weights=parse_weights(args.observable_weight);output=args.output.resolve();run=args.response_run.resolve();reference=args.reference_study.resolve()
    if output.exists() and not args.resume:parser.error('Output exists; choose a new directory or use --resume')
    if args.resume and not (output/'manifest.json').is_file():parser.error('Resume requires a study manifest')
    if (run/'with_bias/cache_manifest.json').is_file():run=run/'with_bias'
    meta=json.loads((reference/'manifest.json').read_text());spec=copy.deepcopy(meta['spec']);cachemeta=json.loads((run/'cache_manifest.json').read_text())
    if cachemeta['version']!=base.VERSION or cachemeta['reference_sha256']!=base.digest(reference/'manifest.json'):raise ValueError('Cache/reference provenance mismatch')
    if base.digest(reference/'histograms.npz')!=meta['histogram_sha256']:raise ValueError('DATA histogram checksum mismatch')
    starts={p:starting_parameters(run/p/'candidate_parameters.json',p) for p in spec['periods']}
    cache={k:dict(r,path=str((run/r['file']).resolve())) for k,r in cachemeta['samples'].items()}
    for p in spec['periods']:
        for sample in ['mc','background','pi0_mc']:
            r=cache[p+'__'+sample]
            if base.digest(r['path'])!=r['sha256']:raise ValueError('Cache checksum mismatch: '+p+'/'+sample)
    for key,name in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        if base.digest(reference/'backend'/name)!=spec['input_hashes'][key]:raise ValueError('Saved production global-cut source changed')
    with np.load(reference/'histograms.npz',allow_pickle=False) as z:data={p+'__data':z[p+'__data'].copy() for p in spec['periods']}
    with (reference/'normalizations.csv').open() as f:norm={(r['period'],base.TOPS.index(r['topology'])):r for r in csv.DictReader(f)}
    if any(r['status']!='ok' for r in norm.values()):raise ValueError('Invalid background normalization')
    signature=dict(reference_sha256=base.digest(reference/'manifest.json'),cache_manifest_sha256=base.digest(run/'cache_manifest.json'),starting_candidates={p:base.digest(run/p/'candidate_parameters.json') for p in spec['periods']},configurations={k:[base.OBS[v] for v in vv] for k,vv in CONFIGURATIONS.items()},weights=weights,fit_settings=dict(rounds=args.rounds,min_mc=args.min_mc,scan_range=args.scan_range,background_transfer_uncertainty=args.background_transfer_uncertainty))
    if args.resume:
        previous=json.loads((output/'manifest.json').read_text())
        if any(previous.get(k)!=v for k,v in signature.items()):raise ValueError('Resume inputs or fit configuration changed')
    output.mkdir(parents=True,exist_ok=True)
    if not args.resume:base.dump(output/'manifest.json',dict(**signature,arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},source_sha256=base.digest(Path(__file__))))
    pending=any(not (output/p/c/'completed_fit.json').is_file() for p in spec['periods'] for c in CONFIGURATIONS)
    scratch=base.scratch_directory(args.scratch_root,output) if pending else None
    worker=base.compile_worker(scratch,reference) if pending else None;budget=base.Budget(args.workers)
    if pending:base.dump(output/'scratch_location.json',dict(directory=str(scratch)))
    print('[study] Existing caches only; with-bias model; fixed starting reference plus 15 equally budgeted period/configuration fits',flush=True)
    results={}
    def work(period,config):
        one=copy.deepcopy(spec);one['periods']=[period];one['reference_flows']=meta['flows'];one['starting_candidate_sha256']=base.digest(run/period/'candidate_parameters.json')
        folder=output/period/config
        if args.resume and (folder/'completed_fit.json').is_file():
            print('[reuse] '+period+'/'+config+' completed fit',flush=True);return period,config,load_checkpoint(folder)
        folder.mkdir(parents=True,exist_ok=True);base.dump(folder/'spec.json',one)
        local=scratch/period/config;local.mkdir(parents=True);pool=base.Pool(worker,folder/'spec.json',local,args.workers,budget)
        engine=TrainingEngine(pool,one,cache,data,norm,base.TOPS,args.min_mc,args.background_transfer_uncertainty,training=CONFIGURATIONS[config],weights=weights)
        result=fit_local(engine,starts[period][0],args,folder,False);save_checkpoint(folder,result);print('[fit] '+period+'/'+config+' complete',flush=True);return period,config,result
    executor=ThreadPoolExecutor(max_workers=args.workers);futures=[executor.submit(work,p,c) for p in spec['periods'] for c in CONFIGURATIONS]
    try:
        for future in as_completed(futures):p,c,result=future.result();results[p,c]=result
    except BaseException:
        budget.cancel_all()
        for future in futures:future.cancel()
        raise
    finally:executor.shutdown(wait=True,cancel_futures=True)
    report(output,spec,data,results,weights)
    print('[done] '+str(output/'REPORT.html')+'; canonical reports only, no duplicate packages',flush=True)

if __name__=='__main__':main()
