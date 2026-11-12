#!/usr/bin/env python3
"""Isolated production-method cut refits for completed particle-response candidates."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,copy,csv,json,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import study_particle_response as response
os.environ['DVCS_EXCLUSIVITY_ARRAY_ONLY']='1'
import plot_exclusivity_data_dvcs_pi0_mc as production

T_EDGES=[0,.11,.15,.25,.40,.60,.80,1.0]

def load_candidate_parameters(path):
    """Expand intentionally fixed parameters; never hide missing fitted values."""
    candidate=json.loads(Path(path).read_text())
    mode=candidate.get('bias_mode','fit');topologies=candidate.get('fit_topologies',response.TOPS)
    if mode not in ('zero','fit') or not topologies or any(t not in response.TOPS for t in topologies):raise ValueError('Invalid candidate response configuration: '+str(path))
    active=response.parameter_grid(topologies,mode);saved=candidate['parameters'];values=[]
    if set(saved)-set(response.PARAMS):raise ValueError('Unknown candidate parameters: '+str(path))
    for index,name in enumerate(response.PARAMS):
        if name in saved:value=float(saved[name])
        elif index not in active:value=0.
        else:raise ValueError('Missing fitted parameter '+name+' in '+str(path))
        if not np.isfinite(value) or (index%4!=3 and value<0):raise ValueError('Invalid response parameter '+name+' in '+str(path))
        if mode=='zero' and index%4==3 and value!=0:raise ValueError('Broadening-only candidate has a nonzero bias: '+str(path))
        values.append(value)
    return values

def candidate_paths(models,periods,training_study=None):
    paths={(p,name):(training_study/p/name/'candidate_parameters.json' if training_study else model/p/'candidate_parameters.json') for name,model in models.items() for p in periods}
    if training_study:
        for (period,name),path in paths.items():
            if json.loads(path.read_text()).get('period')!=period:raise ValueError('Training candidate period mismatch: '+str(path))
            marker=json.loads((path.parent/'completed_fit.json').read_text())
            if response.digest(path)!=marker['candidate_sha256']:raise ValueError('Training candidate changed after completion: '+str(path))
    return paths


def publish_results(run,output):
    """Link the original cut output, including when recovering a legacy run."""
    report=run/'REPORT.html'
    if report.is_file() and 'url=results/REPORT.html' in report.read_text() and (run/'results/REPORT.html').is_file():
        report=run/'results/REPORT.html'  # Follow existing legacy navigation; copy nothing.
    if not report.is_file():return
    response.add_report_link(report,output/'REPORT.html','Redetermined DATA/MC cuts, relative acceptance and two-photon controls ('+output.name+')')


def window_mask(rows,block,kind):
    mask=np.asarray(rows[:,8]>.5)
    for branch,stats in block[kind].items():
        col=1+response.OBS.index(branch);v=rows[:,col]
        if stats.get('mode')=='upper_quantile':mask &= np.isfinite(v)&(v<=float(stats['cut_high']))
        else:mask &= production.apply_window(v,float(stats['cut_low']),float(stats['cut_high']))
    return mask

def fit_arrays(rows,top,selection='global'):
    if selection not in ('global','minimal'):raise ValueError('Unknown fit selection')
    selected=(rows[:,0]==top)&(rows[:,8 if selection=='global' else 10]>.5)
    return {name:np.asarray(rows[selected,i+1],dtype=np.float32) for i,name in enumerate(response.OBS)}

def calibration(pi0_data,pi0_mc,topology):
    """Same direct-pi0 calibrations used by the production iterative fitter.

    xF/xF2 are display-only and cannot enter any production cut or its priors.
    """
    out={}
    for variable in production.PI0_VARIABLES:
        branch='theta_gamma_gamma' if variable.branch=='theta_pi0_pi0' else variable.branch
        if branch not in response.OBS:continue
        dh=np.histogram(pi0_data[branch],bins=variable.bins,range=(variable.xmin,variable.xmax))[0]
        mh=np.histogram(pi0_mc[branch],bins=variable.bins,range=(variable.xmin,variable.xmax))[0]
        fit=production.fit_single_template(dh,mh,variable,topology,12.,20.,100,4.,8.,True,.90)
        if fit.success and np.isfinite(fit.shift) and np.isfinite(fit.sigma_add):out[variable.branch]=(fit.shift,fit.sigma_add)
    return out

def refit(period,top,arrays,folder,selection='global'):
    topology=production.TOPOLOGIES[top]
    dvcs={s:fit_arrays(arrays[key],top,selection) for s,key in [('data','data'),('signal','mc'),('background','background')]}
    pi0={s:fit_arrays(arrays[key],top,selection) for s,key in [('data','pi0_data'),('signal','pi0_mc')]}
    if any(len(a['Delta_phi'])<100 for a in [*dvcs.values(),*pi0.values()]):raise ValueError(f'{period.key}/{topology.key}: insufficient preselection events for production cut fit')
    containments=dict(nominal=.95,tight=.90,loose=.98)
    cal=calibration(pi0['data'],pi0['signal'],topology)
    steps,blocks,final=production.run_dvcs_iterative_cuts(period,topology,dvcs,12.,20.,100,['Delta_phi','theta_gamma_gamma','pTmiss'],4.,8.,True,.90,.95,.90,.95,cal,containments,.25,25.)
    order=[step.variable for step in steps if step.accepted_for_production]
    psteps,pblocks,pfinal=production.run_pi0_iterative_cuts(period,topology,pi0,order,12.,20.,100,4.,8.,True,.90,containments)
    result={name:{production.cut_json_key('dvcs',period,topology):blocks[name],production.cut_json_key('eppi0',period,topology):pblocks[name]} for name in containments}
    folder.mkdir(parents=True,exist_ok=True)
    for name,value in result.items():response.dump(folder/('combined_cuts_'+name+'.json'),value)
    flow=production.cumulative_optimization_rows(period,topology,steps)
    production.write_results_csv(folder/'cut_flow.csv',flow)
    # Invoke the production cut-development plotting helpers, under their shared lock.
    for channel,ss,hists in [('dvcs',steps,final),('pi0',psteps,pfinal)]:
        production.draw_serialized(production.draw_iterative_cut_canvas,folder/(channel+'_cut_development.png'),period,topology,channel,ss,False,110,final_histograms=hists)
    return result['nominal']

def ratio_error(old,new):
    """Paired MC counting error; fitted-response/cut uncertainty is excluded."""
    if old.shape!=new.shape:raise ValueError('Relative acceptance requires aligned event records')
    n0=int(old.sum());n1=int(new.sum())
    if not n0:return float('nan'),float('nan')
    ratio=n1/n0;error=np.sqrt(np.square(new.astype(float)-ratio*old).sum())/n0
    return ratio,float(error)

def comparison_rows(period,scenario,arrays,baseline,new_blocks,old_blocks,probabilities):
    result=[]
    for top,topology in enumerate(response.TOPS):
        code=production.period_code(next(p for p in production.PERIODS if p.key==period));old=old_blocks['DVCS_'+code+'_'+topology];new=new_blocks['DVCS_'+code+'_'+topology]
        base=baseline['mc'];mc=arrays['mc'];data=arrays['data']
        base_selected=window_mask(base,old,'mc')&(base[:,0]==top)
        selected=window_mask(mc,new,'mc')&(mc[:,0]==top)
        old_data=window_mask(data,old,'data')&(data[:,0]==top);new_data=window_mask(data,new,'data')&(data[:,0]==top)
        for lo,hi in [(0.,1.),*zip(T_EDGES[:-1],T_EDGES[1:])]:
            b=base_selected&(base[:,9]>=lo)&(base[:,9]<hi);n=selected&(mc[:,9]>=lo)&(mc[:,9]<hi)
            ratio,error=ratio_error(b,n);dm=(data[:,9]>=lo)&(data[:,9]<hi)
            result.append(dict(period=period,scenario=scenario,topology=topology,t_low=lo,t_high=hi,baseline_mc_selected=int(b.sum()),candidate_mc_selected=int(n.sum()),relative_acceptance=ratio,paired_mc_stat_error=error,inverse_acceptance_factor=1/ratio if ratio>0 else float('nan'),old_data_selected=int((old_data&dm).sum()),new_data_selected=int((new_data&dm).sum()),mc_sampling_probability=probabilities['mc']))
    return result

def plot_controls(period,scenario,arrays,baseline,blocks,output,selection='global'):
    import matplotlib.pyplot as plt
    paths=[]
    for top,topology in enumerate(response.TOPS):
        if not any(np.any((a[:,0]==top)&(a[:,8 if selection=='global' else 10]>.5)) for a in [arrays['pi0_data'],arrays['pi0_mc'],baseline['pi0_mc']]):continue
        fig,axes=plt.subplots(3,3,figsize=(15,10),constrained_layout=True)
        samples=[('DATA',arrays['pi0_data'],'#205cab'),('Nominal MC',baseline['pi0_mc'],'#d45c18'),('Candidate MC',arrays['pi0_mc'],'#22946b')]
        for i,name in enumerate(response.OBS):
            ax=axes.flat[i];values=[]
            for label,rows,color in samples:
                # Common minimal preselection, before exclusivity: avoids hiding tails behind refitted cuts.
                a=rows[(rows[:,0]==top)&(rows[:,8 if selection=='global' else 10]>.5),i+1];values.append(a[np.isfinite(a)])
            qs=[np.quantile(a,[.01,.99]) for a in values if len(a)]
            lo=min(q[0] for q in qs);hi=max(q[1] for q in qs)
            if name=='Mx2':
                cores=[np.quantile(a,[.16,.5,.84]) for a in values if len(a)];width=max(max(q[1]-q[0],q[2]-q[1],.0005) for q in cores);center=float(np.median([q[1] for q in cores]));lo=max(-.06,center-4*width);hi=min(.06,center+4*width)
            if hi<=lo:lo-=.001;hi+=.001
            edges=np.linspace(lo,hi,151)
            for (label,rows,color),a in zip(samples,values):
                counts=np.histogram(a,edges)[0];ax.stairs(counts/(max(1,len(a))*np.diff(edges)),edges,label=label,color=color)
            key='eppi0_'+production.period_code(next(p for p in production.PERIODS if p.key==period))+'_'+topology
            for kind,color in [('data','#205cab'),('mc','#22946b')]:
                stats=blocks[key][kind].get(name)
                if stats:
                    for x in [stats['cut_low'],stats['cut_high']]:
                        if lo<x<hi:ax.axvline(x,color=color,ls='--',alpha=.5,lw=.8)
            labels={'Delta_phi':'Coplanarity Δφ (rad)','theta':'CM angle θ (rad)','theta_gamma_gamma':'π⁰ opening angle (deg)','pTmiss':'Missing transverse momentum (GeV)','Emiss2':'Missing energy (GeV)','Mx2':'Missing mass squared epπ⁰ (GeV²)','Mx2_2':'Missing mass squared eπ⁰ (GeV²)'}
            ax.set(xlabel=labels[name],ylabel='Unit-area density');ax.legend(fontsize=7);ax.grid(alpha=.15)
        for ax in axes.flat[7:]:ax.axis('off')
        axes.flat[7].text(0,1,'Two-photon control: selected phase space, before exclusivity cuts.\nSeparate DATA/MC refitted boundaries are dashed.\nCurves normalized by all finite events;\nrange truncation is not renormalized.',va='top')
        fig.suptitle(period+' · '+scenario+' · '+topology+' · two-photon control')
        path=output/(period+'_'+scenario+'_'+topology+'_control.png');fig.savefig(path,dpi=110);plt.close(fig);paths.append(path.name)
    return paths

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['reference-study','response-run','output']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--fit-selection',choices=['global','minimal'],default='global',help='Selection for cut determination; global matches acceptance (default)')
    parser.add_argument('--export-only',action='store_true',help='Export cached/input events and stop before fitting cuts')
    parser.add_argument('--training-study',type=Path,help='Use checked current and add_pt candidates from the completed training-variable study')
    parser.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'));parser.add_argument('--workers',type=int,default=10)
    args=parser.parse_args(argv)
    if not 1<=args.workers<=10:parser.error('Workers must be between 1 and 10')
    output=args.output.resolve()
    if output.exists():parser.error('Output exists; choose a new directory')
    reference=args.reference_study.resolve();run=args.response_run.resolve();meta=json.loads((reference/'manifest.json').read_text());spec=copy.deepcopy(meta['spec'])
    models={name:run/name for name in ['broadening_only','with_bias'] if (run/name/'cache_manifest.json').is_file()}
    if not models:models={run.name:run}
    if args.training_study:
        if (run/'with_bias/cache_manifest.json').is_file():run=run/'with_bias'
        models={name:run for name in ['current','add_pt']}
    first=next(iter(models.values()));cachemeta=json.loads((first/'cache_manifest.json').read_text())
    for name,model in models.items():
        other=json.loads((model/'cache_manifest.json').read_text())
        if other['reference_sha256']!=cachemeta['reference_sha256'] or other['seed']!=cachemeta['seed'] or {k:r['sha256'] for k,r in other['samples'].items()}!={k:r['sha256'] for k,r in cachemeta['samples'].items()}:raise ValueError('Candidates must use the same verified event cache')
    if cachemeta['reference_sha256']!=response.digest(reference/'manifest.json'):raise ValueError('Response and reference provenance differ')
    paths=candidate_paths(models,spec['periods'],args.training_study.resolve() if args.training_study else None)
    if args.training_study:
        training_meta=json.loads((args.training_study/'manifest.json').read_text())
        if training_meta['reference_sha256']!=response.digest(reference/'manifest.json') or training_meta['cache_manifest_sha256']!=response.digest(first/'cache_manifest.json'):raise ValueError('Training/cut input provenance mismatch')
    candidate_vectors={key:load_candidate_parameters(path) for key,path in paths.items()}
    output.mkdir(parents=True);work=response.scratch_directory(args.scratch_root,output);response.dump(output/'scratch_location.json',dict(directory=str(work)))
    response.dump(output/'spec.json',spec);worker=response.compile_worker(output,reference);pool=response.Pool(worker,output/'spec.json',work,args.workers)
    cache={k:dict(r,path=str((first/r['file']).resolve())) for k,r in cachemeta['samples'].items()}
    for k,r in cache.items():
        if response.digest(r['path'])!=r['sha256']:raise ValueError('Cache checksum mismatch: '+k)
    jobs=[];descriptors=[]
    for period in spec['periods']:
        for sample in ['data','pi0_data']:
            entries=[r for r in meta['inventory'] if r['period']==period and r['sample']==sample]
            if not entries:raise ValueError('Missing DATA inventory '+period+'/'+sample)
            for r in entries:
                st=Path(r['path']).stat()
                if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):raise ValueError('DATA input changed')
            jobs.append(dict(period=period,sample=sample,files=[r['path'] for r in entries]));descriptors.append((period,'shared',sample))
        for scenario in ['nominal',*models]:
            pars=[0.]*20 if scenario=='nominal' else candidate_vectors[period,scenario]
            for sample in ['mc','background','pi0_mc']:
                jobs.append(dict(period=period,sample=sample,cache=cache[period+'__'+sample]['path'],parameters=pars));descriptors.append((period,scenario,sample))
    exported={}
    for index,job,path in pool.run('cut-arrays',jobs,'precut_arrays'):
        shape=json.loads(Path(str(path)+'.json').read_text())['shape'];exported[descriptors[index]]=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape))
    def arrays_for(period,scenario):return {sample:exported[period,'shared' if sample in ['data','pi0_data'] else scenario,sample] for sample in ['data','mc','background','pi0_data','pi0_mc']}
    if args.export_only:
        print('[done] Pre-cut event exports: '+str(work),flush=True);return
    tasks={};results={}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for period in spec['periods']:
            p=next(p for p in production.PERIODS if p.key==period)
            for scenario in ['nominal',*models]:
                a=arrays_for(period,scenario)
                for top,topology in enumerate(response.TOPS):tasks[executor.submit(refit,p,top,a,output/period/scenario/topology,args.fit_selection)]=(period,scenario,topology)
        for future in as_completed(tasks):
            p,s,t=tasks[future];results.setdefault((p,s),{}).update(future.result());print('[cuts] '+p+'/'+s+'/'+t+' complete',flush=True)
    rows=[];links=[]
    for period in spec['periods']:
        baseline=arrays_for(period,'nominal');prob={s:cache[period+'__'+s]['probability'] for s in ['mc','background','pi0_mc']}
        for scenario in ['nominal',*models]:
            blocks=results[period,scenario];folder=output/period/scenario
            response.dump(folder/'combined_cuts_95.json',blocks)
            rows+=comparison_rows(period,scenario,arrays_for(period,scenario),baseline,blocks,spec['exclusivity'],prob)
            # Also compare each candidate to a nominal refit using the same larger cache and method.
            if scenario!='nominal':
                z=comparison_rows(period,scenario,arrays_for(period,scenario),baseline,blocks,results[period,'nominal'],prob)
                for r in z:r['scenario']+='__versus_nominal_refit'
                rows+=z
            if scenario=='add_pt' and 'current' in models:
                z=comparison_rows(period,scenario,arrays_for(period,scenario),arrays_for(period,'current'),blocks,results[period,'current'],prob)
                for r in z:r['scenario']+='__versus_current_training_refit'
                rows+=z
            plots=output/'control_plots';plots.mkdir(exist_ok=True)
            links+=plot_controls(period,scenario,arrays_for(period,scenario),baseline,blocks,plots,args.fit_selection)
    response.table(output/'relative_acceptance.csv',rows)
    response.dump(output/'manifest.json',dict(reference_sha256=response.digest(reference/'manifest.json'),response_cache_sha256=response.digest(first/'cache_manifest.json'),production_fitter_sha256=response.digest(Path(production.__file__)),candidate_hashes={name:{p:response.digest(paths[p,name]) for p in spec['periods']} for name,model in models.items()},scenarios=['nominal',*models],arguments={k:str(v) for k,v in vars(args).items()}))
    (output/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Production cut redetermination</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}li{margin:8px}</style><h1>Separate DATA/MC cut redetermination</h1><p>Uses the production template fitter and iterative 90/95/98% cut algorithms. DATA are unchanged; their signal templates and fitted boundaries are reestimated for each MC response. Cut fitting uses the requested selection (global by default); relative-acceptance counts always use the saved full analysis global cuts. Production outputs are untouched.</p><p><a href="relative_acceptance.csv">Relative reconstructed acceptance by period, topology and |t|</a>. The generated denominator is unchanged and cancels in these ratios. The inverse acceptance factor alone is not a cross-section prediction: DATA cuts, background subtraction, bin migration and fitted-parameter uncertainties also matter. Statistical errors are paired MC counting errors only. The original upstream two-photon mass-window limitation remains.</p><h2>Two-photon controls before exclusivity cuts</h2><ul>'+''.join('<li><a href="control_plots/'+x+'">'+x+'</a></li>' for x in links)+'</ul><p>Per-period/scenario folders contain combined_cuts_95.json and topology-specific production cut-development figures.</p>')
    if args.training_study:
        report=output/'REPORT.html';report.write_text(report.read_text()+'<p>Training follow-up: current is the equally refined original-training candidate; add_pt adds missing transverse momentum to training. Rows ending __versus_current_training_refit isolate that training-choice change after both DATA and MC cuts are redetermined.</p>')
    else:publish_results(run,output)
    print('[done] '+str(output/'REPORT.html'),flush=True)

if __name__=='__main__':main()
