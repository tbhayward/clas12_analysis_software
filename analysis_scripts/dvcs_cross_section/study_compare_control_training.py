#!/usr/bin/env python3
"""Equally refined epgamma-only versus joint training; both channels validated."""
from __future__ import annotations
import argparse,copy,csv,json,subprocess,sys
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_response_cut_redetermination as cuts


def read(path):
    with path.open() as f:return list(csv.DictReader(f))


def check_joint(run):
    meta=json.loads((run/'manifest.json').read_text());args=meta['arguments']
    if not args.get('joint_control_study') or not args.get('local_refinement') or args.get('joint_control_weight',.5)!=.5:raise ValueError('Use a completed default-weight local joint fit')
    for n in ('particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h','particle_response_cut_arrays.h'):
        if base.digest(base.SOURCE/n)!=meta['source_hashes'][n]:raise ValueError('MC evaluator changed since joint fit: '+n)
    spec=json.loads((run/'spec.json').read_text())
    for period in spec['periods']:
        folder=run/period;scores=read(folder/'single_parameter_scores.csv');pars=json.loads((folder/'candidate_parameters.json').read_text())
        start=json.loads(Path(spec['joint_start_paths'][period]).read_text());values=[start['parameters'].get(k,0.) for k in base.PARAMS]
        if len(scores)!=2 or json.loads(scores[1]['parameters'])!=values or float(scores[1]['train'])>float(scores[0]['train']):raise ValueError('Cannot establish same warm start for '+period)
        if not pars.get('joint_control_fit'):raise ValueError('Not a joint candidate')
    return meta,spec


def comparison_plot(period,engine,hist,selected,candidates,out):
    import matplotlib.pyplot as plt
    records=[];summary=[];links=[]
    for name,index in [('epgamma_only',1),('joint',2)]:
        params=json.loads(candidates[name].read_text());m=float(params['background_multiplier']);score,rr=engine.score_one(hist,index,m)
        for r in rr:r.update(period=period,fit=name)
        records+=rr
        for ch in ('epgamma','epgammagamma'):
            for metric in ('train','validation_mc','heldout_observables'):summary.append(dict(period=period,fit=name,channel=ch,metric=metric,score=score[ch+'_'+metric]))
    perioddir=out/'plots'/period;perioddir.mkdir(parents=True)
    fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for ax,ch in zip(axes,('epgamma','epgammagamma')):
        for i,name in enumerate(('epgamma_only','joint')):
            vals=[next(r['score'] for r in summary if r['fit']==name and r['channel']==ch and r['metric']==metric) for metric in ('train','validation_mc','heldout_observables')];ax.bar(np.arange(3)+(i-.5)*.35,vals,width=.35,label=name,color=['#d45c18','#22946b'][i])
        ax.set(xticks=range(3),xticklabels=['Training variables','Other MC half','Reserved variables'],ylabel='Shape score (lower is better)',title=ch);ax.tick_params(axis='x',rotation=15);ax.legend();ax.grid(axis='y',alpha=.2)
    fig.suptitle(period+' · equal-start/equal-refinement comparison');path=perioddir/'scores.png';fig.savefig(path,dpi=120);plt.close(fig);links.append(path.relative_to(out).as_posix())
    for channel in ('epgamma','epgammagamma'):
        for t in engine.fit_topologies:
            fig,axes=plt.subplots(3,3,figsize=(14,10),layout='constrained')
            for vi,v in enumerate(engine.spec['variables']):
                ax=axes.flat[vi];n=v['bins']+2;edges=base.diag.variable_edges(v);d=engine.data[period+('__data' if channel=='epgamma' else '__pi0_data')][t,0,vi,:n]
                curves=[('DATA',d,'#205cab')]
                for name,index,color in [('epgamma_only',1,'#d45c18'),('joint',2,'#22946b')]:
                    if channel=='epgammagamma':mc=hist[period+'__pi0_mc'][index,:,t,0,vi,:n].sum(axis=0)
                    else:
                        multiplier=float(json.loads(candidates[name].read_text())['background_multiplier']);mc=sum(engine.weighted(hist,period,index,half,t,0,vi,multiplier)[0] for half in (0,1))
                    curves.append((name,mc,color))
                for label,h,color in curves:
                    if h.sum()>0:ax.stairs(h[1:-1]/h.sum()/np.diff(edges),edges,color=color,label=label)
                qs=[base.diag.quantile(h[1:-1],edges,q) for _,h,_ in curves if h.sum()>0 for q in (.01,.99)];qs=[x for x in qs if np.isfinite(x)]
                if qs:
                    lo,hi=min(qs),max(qs);pad=max((hi-lo)*.08,.001);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                if v['name']=='Mx2':ax.set_xlim(*base.mass_plot_limits(edges,[h for _,h,_ in curves]))
                ax.set(xlabel=v['label'],ylabel='Unit-area density');ax.legend(fontsize=7);ax.grid(alpha=.15)
            for ax in axes.flat[7:]:ax.axis('off')
            fig.suptitle(period+' · '+base.TOPS[t]+' · '+channel+' DATA versus both fitted predictions\nEpγγ uses two-photon AAOGEN MC; both curves receive identical validation definitions')
            path=perioddir/(channel+'_'+base.TOPS[t]+'.png');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
    return records,summary,links


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--joint-run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=10);args=p.parse_args(argv)
    if not 1<=args.workers<=10:p.error('Workers must be 1–10 total')
    joint=args.joint_run.resolve();out=args.output.resolve();meta,spec=check_joint(joint);old=meta['arguments']
    if out.exists():p.error('Output exists; choose a fresh directory')
    out.mkdir(parents=True);base.dump(out/'manifest.json',dict(arguments={k:str(v) for k,v in vars(args).items()},joint_sha256=base.digest(joint/'manifest.json'),comparison='Only objective control weight changes: 0 versus .5; validation includes both channels'))
    solo=out/'epgamma_only_fit'
    cmd=[sys.executable,str(base.SOURCE/'study_particle_response.py'),'--reference-study',old['reference_study'],'--joint-control-study',old['joint_control_study'],'--reuse-cache',str(joint),'--joint-control-weight','0','--local-refinement','--output',str(solo),'--workers',str(args.workers)]
    for flag in ('expected_cache_events','seed','scan_range','bias_mode','background_transfer_uncertainty','coordinate_rounds','refinement_rounds','min_mc','scratch_root'):cmd+=['--'+flag.replace('_','-'),str(old[flag])]
    cmd+=['--fit-topologies',*old['fit_topologies']]
    print('[fit] Reusing completed joint branch; fitting only the equally budgeted epgamma-only branch. Two-photon shapes remain validation.',flush=True);subprocess.run(cmd,check=True)
    solo_spec=json.loads((solo/'spec.json').read_text());solo_cache=json.loads((solo/'cache_manifest.json').read_text());joint_cache=json.loads((joint/'cache_manifest.json').read_text())
    if {k:v['sha256'] for k,v in solo_cache['samples'].items()}!={k:v['sha256'] for k,v in joint_cache['samples'].items()}:raise ValueError('Fit caches differ')
    for period in spec['periods']:
        starts=read(solo/period/'single_parameter_scores.csv');oldstarts=read(joint/period/'single_parameter_scores.csv')
        if json.loads(starts[1]['parameters'])!=json.loads(oldstarts[1]['parameters']) or float(starts[1]['train'])>float(starts[0]['train']):raise ValueError('Starting candidate differed: '+period)
    reference=Path(old['reference_study'])
    with np.load(reference/'histograms.npz',allow_pickle=False) as ar:data={period+'__data':ar[period+'__data'].copy() for period in spec['periods']}
    data.update(base.load_joint_control_data(Path(old['joint_control_study']),solo_spec,reference));norm={(r['period'],base.TOPS.index(r['topology'])):r for r in read(reference/'normalizations.csv')}
    base.TRAIN_OBS=[0,3,4,5];base.VALID_OBS=[2,6];work=base.scratch_directory(Path(old['scratch_root']),out);base.dump(out/'scratch_location.json',dict(directory=str(work)))
    worker=solo/'backend/particle_response_worker';cache={k:dict(v,path=str((solo/v['file']).resolve())) for k,v in solo_cache['samples'].items()};records=[];summaries=[];links=[]
    for period in spec['periods']:
        one=copy.deepcopy(solo_spec);one.update(periods=[period],joint_control_weight=.5);sp=work/(period+'.json');base.dump(sp,one);directory=work/period;directory.mkdir();pool=base.Pool(worker,sp,directory,args.workers)
        engine=base.Engine(pool,one,cache,data,norm,old['fit_topologies'],old['min_mc'],old['background_transfer_uncertainty']);paths={'epgamma_only':solo/period/'candidate_parameters.json','joint':joint/period/'candidate_parameters.json'}
        pars=[[0.]*20,*[cuts.load_candidate_parameters(paths[name]) for name in ('epgamma_only','joint')]]
        unique,_,_,hist,selected=engine.evaluate(pars,'paired_validation',detailed=False)
        # Identical candidates can be deduplicated by the worker; map explicitly.
        mapped={};
        for sample in ('mc','background','pi0_mc'):
            h=hist[period+'__'+sample];indices=[unique.index(tuple(x)) for x in pars];mapped[period+'__'+sample]=h[indices]
        engine.current_selected={k:v[[unique.index(tuple(x)) for x in pars]] for k,v in selected.items()};engine.control=engine.control[[unique.index(tuple(x)) for x in pars]]
        rr,ss,ll=comparison_plot(period,engine,mapped,selected,paths,out);records+=rr;summaries+=ss;links+=ll;base.discard_histograms(hist)
        print('[done] '+period+' both fitted candidates validated in both channels',flush=True)
    base.table(out/'projection_comparison.csv',records);base.table(out/'score_comparison.csv',summaries)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Paired control-training test</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>Epγ-only versus joint training</h1><p>Both branches start from the same checked pre-joint candidate, use the same event cache, train/reserved variables, bounds, period/topology setup, normalization/nuisance rules and coordinate/refinement schedule. Only the control fraction of the shape objective changes from zero to one half. The completed joint branch is reused after checking its evaluator provenance and start. Both branches still morph two-photon MC for the common background normalization.</p><p>Validation explicitly compares two-photon DATA with AAOGEN MC as well as epgamma DATA with DVCS+leakage MC. Score bars separate channels and training/other-MC-half/reserved observables. Shape overlays show DATA and both fitted models by topology. Validation MC half and reserved variables never select parameters; experimental DATA are shared, so this is not independent-DATA validation. Shapes are area-normalized and do not validate absolute efficiencies or contamination.</p><p>Combined evaluation uses frozen fitted nuisance multipliers from each branch, not a new fit on validation. The 0.5 combined score is descriptive; compare each channel and topology before accepting a tradeoff. Both branches have equal refinement schedules, not necessarily identical counts of unique proposals at parameter boundaries. No production cuts or corrections change.</p><p><a href="score_comparison.csv">Channel score comparison</a> · <a href="projection_comparison.csv">Every projection</a></p>'+''.join('<p><a href="'+x+'">'+x+'</a></p>' for x in links))
    print('[done] '+str(out/'REPORT.html'),flush=True)

if __name__=='__main__':main()
