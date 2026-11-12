#!/usr/bin/env python3
"""Validate the add-pT response in kinematic subdivisions, controls and cut refits."""
from __future__ import annotations
import argparse,copy,csv,json,subprocess,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import study_particle_response as base
import study_response_cut_redetermination as cuts

SINGLE_AXES={'e_p','e_theta','p1_p','p1_theta','p2_p','p2_theta'}
AXIS_LABELS={'e_p':'Electron p (GeV)','e_theta':'Electron θ (deg)','p1_p':'Proton p (GeV)','p1_theta':'Proton θ (deg)','p2_p':'Photon E (GeV)','p2_theta':'Photon θ (deg)'}

def cell_metrics(data,models,variances,nominal,nominal_variance,raw_nominal,min_mc):
    scores=[base.shape_score(data,m,v,int(raw_nominal),min_mc,nominal,nominal_variance)[0] for m,v in zip(models,variances)]
    usable=data.sum()>=300 and all(np.isfinite(scores))
    gain=100*(1-scores[1]/scores[0]) if usable and scores[0]>0 else float('nan')
    return scores,gain,'usable descriptive comparison' if usable else 'insufficient statistics/range coverage'

def scan_period(period,spec,cache,data,norm,candidate_paths,worker,work,output,budget,workers,min_mc):
    one=copy.deepcopy(spec);one['periods']=[period];folder=output/period;folder.mkdir()
    base.dump(folder/'spec.json',one);local=work/period;local.mkdir()
    saved={name:json.loads(path.read_text()) for name,path in candidate_paths.items()}
    pars=[cuts.load_candidate_parameters(candidate_paths[name]) for name in ['current','add_pt']]
    pool=base.Pool(worker,folder/'spec.json',local,workers,budget);engine=base.Engine(pool,one,cache,data,norm,base.TOPS,min_mc)
    points,scores,detail,hist,selected=engine.evaluate(pars,'subdivisions',detailed=True)
    rows=[];examples={};indices=[points.index(tuple(p)) for p in pars]
    try:
        for top,topology in enumerate(base.TOPS):
            for g,stratum in enumerate(spec['strata']):
                for vi,variable in enumerate(spec['variables']):
                    n=variable['bins']+2;d=data[period+'__data'][top,g,vi,:n].astype(float)
                    models=[];variances=[]
                    for name,index in zip(['current','add_pt'],indices):
                        parts=[engine.weighted(hist,period,index,h,top,g,vi,saved[name]['background_multiplier']) for h in (0,1)]
                        models.append((parts[0][0]+parts[1][0])/2);variances.append((parts[0][1]+parts[1][1])/4)
                    ref=[engine.weighted(hist,period,0,h,top,g,vi,reference_scales=True) for h in (0,1)]
                    nominal=(ref[0][0]+ref[1][0])/2;variance=(ref[0][1]+ref[1][1])/4
                    raw=sum(hist[period+'__'+s][0,:,top,g,vi,:n].sum() for s in ['mc','background'])
                    ss,gain,status=cell_metrics(d,models,variances,nominal,variance,raw,min_mc)
                    row=dict(period=period,topology=topology,axis=stratum['axis'],stratum=stratum['label'],variable=variable['name'],data_entries=float(d.sum()),nominal_mc_entries=int(raw),current_score=ss[0],add_pt_score=ss[1],improvement_percent=gain,status=status)
                    rows.append(row)
                    if stratum['axis'] in SINGLE_AXES and status=='usable descriptive comparison':
                        key=(topology,vi,g);examples[key]=(row,d.copy(),[m.copy() for m in models])
        # Keep only the three largest reserved-variable regressions per topology for plots.
        retained={}
        for topology in base.TOPS:
            candidates=[(k,v) for k,v in examples.items() if k[0]==topology and v[0]['variable'] in ['theta_gamma_gamma','Mx2_2']]
            candidates.sort(key=lambda x:x[1][0]['add_pt_score']-x[1][0]['current_score'],reverse=True)
            retained.update(candidates[:3])
        base.table(folder/'subdivision_scores.csv',rows)
        return rows,retained
    finally:base.discard_histograms(hist)

def render(output,spec,rows,examples):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    links=[]
    labels={v['name']:v['label'] for v in spec['variables']}
    for period in spec['periods']:
        for topology in base.TOPS:
            z=[r for r in rows if r['period']==period and r['topology']==topology and r['axis'] in SINGLE_AXES]
            cells=list(dict.fromkeys((r['axis'],r['stratum']) for r in z));variables=base.OBS
            values=np.array([[next(r['improvement_percent'] for r in z if (r['axis'],r['stratum'])==cell and r['variable']==v) for v in variables] for cell in cells])
            fig,ax=plt.subplots(figsize=(12,max(5,len(cells)*.25)),layout='constrained');im=ax.imshow(values,cmap='RdYlGn',vmin=-50,vmax=50,aspect='auto')
            ax.set(xticks=range(7),xticklabels=variables,yticks=range(len(cells)),yticklabels=[AXIS_LABELS[a]+' '+label for a,label in cells]);ax.tick_params(axis='y',labelsize=8)
            for (i,j),v in np.ndenumerate(values):ax.text(j,i,f'{v:+.0f}' if np.isfinite(v) else 'low N',ha='center',va='center',fontsize=7)
            ax.set_title(period+' · '+topology+' · add-pT improvement (%)\nPositive = better; gray cells lack statistics');im.cmap.set_bad('#dddddd');fig.colorbar(im,ax=ax,shrink=.6,label='Score improvement (%)')
            plot=output/'plots'/(period+'_'+topology+'_subdivisions.png');plot.parent.mkdir(exist_ok=True);fig.savefig(plot,dpi=120);plt.close(fig);links.append(plot.relative_to(output).as_posix())
            cases=[v for k,v in examples[period].items() if k[0]==topology]
            if cases:
                fig,axes=plt.subplots(1,len(cases),figsize=(5*len(cases),4),squeeze=False,layout='constrained')
                for ax,(row,d,models) in zip(axes.flat,cases):
                    vi=variables.index(row['variable']);edges=base.diag.variable_edges(spec['variables'][vi])
                    for label,m,color in [('DATA',d,'#205cab'),('Original training',models[0],'#d45c18'),('+ pT training',models[1],'#22946b')]:
                        if m.sum()>0:ax.stairs(m[1:-1]/m.sum()/np.diff(edges),edges,label=label,color=color)
                    qs=[base.diag.quantile(h[1:-1],edges,q) for h in [d,*models] for q in [.01,.99]];qs=[q for q in qs if np.isfinite(q)]
                    if qs:
                        lo,hi=min(qs),max(qs);pad=max((hi-lo)*.1,.001);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                    ax.set(xlabel=labels[row['variable']],ylabel='Unit-area density',title=AXIS_LABELS[row['axis']]+' '+row['stratum']);ax.legend(fontsize=8);ax.grid(alpha=.15)
                fig.suptitle(period+' · '+topology+' · most challenging reserved-variable checks')
                plot=output/'plots'/(period+'_'+topology+'_reserved_checks.png');fig.savefig(plot,dpi=120);plt.close(fig);links.append(plot.relative_to(output).as_posix())
    summary=[]
    for p in spec['periods']:
        for t in base.TOPS:
            for variable in base.OBS:
                z=[r for r in rows if r['period']==p and r['topology']==t and r['variable']==variable and r['axis'] in SINGLE_AXES and r['status']=='usable descriptive comparison']
                summary.append(dict(period=p,topology=t,variable=variable,usable_cells=len(z),improved_cells=sum(r['improvement_percent']>0 for r in z),worsened_cells=sum(r['improvement_percent']<0 for r in z),median_improvement_percent=float(np.median([r['improvement_percent'] for r in z])) if z else float('nan')))
    base.table(output/'subdivision_summary.csv',summary)
    (output/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Add-pT follow-up</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}li{margin:8px}</style><h1>Add-pT candidate: kinematic validation and cut study</h1><p>Compares the equally refined original-training candidate with the add-pT candidate. No response fit or parameter search is repeated. DATA and MC use the same saved diagnostic windows for these subdivision checks. The checks share DATA events and are not independent DATA validation; momentum/angle cells overlap across axes. Scores use fixed nominal-MC uncertainty and shape normalization. Low-statistics cells are flagged.</p><p><a href="subdivision_summary.csv">Subdivision summary</a> · <a href="subdivision_scores.csv">All single-axis and joint-cell scores</a></p><ul>'+''.join('<li><a href="'+p+'">'+Path(p).stem+'</a></li>' for p in links)+'</ul>')

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--training-study',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--reference-study',type=Path);parser.add_argument('--response-run',type=Path);parser.add_argument('--workers',type=int,default=10);parser.add_argument('--min-mc',type=int,default=100);parser.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'))
    args=parser.parse_args(argv)
    if not 1<=args.workers<=10 or args.min_mc<10:parser.error('Invalid workers/statistics settings')
    training=args.training_study.resolve();saved=json.loads((training/'manifest.json').read_text());output=args.output.resolve()
    if output.exists():parser.error('Output exists; choose a new directory')
    reference=(args.reference_study or Path(saved['arguments']['reference_study'])).resolve();run=(args.response_run or Path(saved['arguments']['response_run'])).resolve()
    if (run/'with_bias/cache_manifest.json').is_file():run=run/'with_bias'
    meta=json.loads((reference/'manifest.json').read_text());cachemeta=json.loads((run/'cache_manifest.json').read_text());spec=copy.deepcopy(meta['spec']);spec['reference_flows']=meta['flows']
    if base.digest(reference/'manifest.json')!=saved['reference_sha256'] or base.digest(run/'cache_manifest.json')!=saved['cache_manifest_sha256']:raise ValueError('Training/reference/cache provenance differs')
    if base.digest(reference/'histograms.npz')!=meta['histogram_sha256']:raise ValueError('DATA histograms changed')
    paths=cuts.candidate_paths({c:run for c in ['current','add_pt']},spec['periods'],training)
    cache={k:dict(r,path=str((run/r['file']).resolve())) for k,r in cachemeta['samples'].items()}
    for record in cache.values():
        if base.digest(record['path'])!=record['sha256']:raise ValueError('MC cache changed')
    for key,name in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        if base.digest(reference/'backend'/name)!=spec['input_hashes'][key]:raise ValueError('Saved global-cut sources changed')
    with np.load(reference/'histograms.npz',allow_pickle=False) as archive:data={p+'__data':archive[p+'__data'].copy() for p in spec['periods']}
    with (reference/'normalizations.csv').open() as f:norm={(r['period'],base.TOPS.index(r['topology'])):r for r in csv.DictReader(f)}
    if any(r['status']!='ok' for r in norm.values()):raise ValueError('Invalid component normalization')
    output.mkdir(parents=True);work=base.scratch_directory(args.scratch_root,output);worker=base.compile_worker(work,reference);budget=base.Budget(args.workers)
    base.dump(output/'manifest.json',dict(training_manifest_sha256=base.digest(training/'manifest.json'),candidate_hashes={p:{c:base.digest(paths[p,c]) for c in ['current','add_pt']} for p in spec['periods']},arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}));base.dump(output/'scratch_location.json',dict(directory=str(work)))
    rows=[];examples={};executor=ThreadPoolExecutor(max_workers=min(args.workers,len(spec['periods'])))
    futures=[executor.submit(scan_period,p,spec,cache,data,norm,{c:paths[p,c] for c in ['current','add_pt']},worker,work,output,budget,args.workers,args.min_mc) for p in spec['periods']]
    try:
        for future in as_completed(futures):
            rr,ee=future.result();rows+=rr;examples[rr[0]['period']]=ee
    except BaseException:
        budget.cancel_all()
        for future in futures:future.cancel()
        raise
    finally:executor.shutdown(wait=True,cancel_futures=True)
    base.table(output/'subdivision_scores.csv',rows);render(output,spec,rows,examples)
    subprocess.run([sys.executable,str(Path(cuts.__file__).resolve()),'--reference-study',str(reference),'--response-run',str(run),'--training-study',str(training),'--workers',str(args.workers),'--scratch-root',str(args.scratch_root),'--output',str(output/'cut_redetermination')],check=True)
    base.add_report_link(output/'REPORT.html',output/'cut_redetermination/REPORT.html','Two-photon controls, production cut redetermination and relative acceptance')
    print('[done] '+str(output/'REPORT.html'),flush=True)

if __name__=='__main__':main()
