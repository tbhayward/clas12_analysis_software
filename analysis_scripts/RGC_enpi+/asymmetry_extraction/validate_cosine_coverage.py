#!/usr/bin/env python3
"""Two-dimensional u1/ll1 injections at observed kinematics; separate PNG diagnostics."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
from importlib.metadata import version
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import socket
import sys
import time

for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(variable,'1')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import validate_unpolarized_modulations as support
prod = support.prod
HERE = Path(__file__).resolve().parent
TARGETS = ('u1','ll1')
MODES = ('observed','matched_uniform','matched_forward_gap')
_WORK = None


def fit_start_values(nominal):
    # Identical starts for every grid point: injected target values never seed fits.
    result = {p:float(nominal['values'][p]) for p in prod.PHYSICS_PARAMETERS}
    result.update(u1=0.,ll1=0.)
    return result


def make_truth(nominal, u1, ll1):
    result = {p:float(nominal['values'][p]) for p in prod.PHYSICS_PARAMETERS}
    result.update(u1=float(u1),ll1=float(ll1))
    return result


def coverage_indices(events, mode, config, b, rep):
    """Match each period's counts, selecting without using run/helicity labels."""
    if mode == 'observed':
        return np.arange(len(events['phi']))
    if mode not in MODES:
        raise ValueError(mode)
    rng = np.random.default_rng(np.random.SeedSequence([config['seed'],b,rep,991]))
    ranks = rng.random(len(events['phi']))
    angle = np.abs((events['phi'] + np.pi) % (2*np.pi)-np.pi)
    selected = []
    for pi in range(len(prod.PERIODS)):
        indices = np.flatnonzero(events['period_index']==pi)
        if not len(indices):
            continue
        count = int(math.floor(config['retained_fraction']*len(indices)))
        eligible = indices if mode == 'matched_uniform' else indices[angle[indices]>=np.deg2rad(config['gap_half_width'])]
        if count < 1 or len(eligible) < count:
            raise ValueError(f'Coverage {mode}: period {pi} has {len(eligible)} eligible events, needs {count}')
        selected.extend(eligible[np.argsort(ranks[eligible],kind='stable')[:count]].tolist())
    return np.sort(np.asarray(selected,dtype=int))


def subset_events(events, indices):
    return {k:v[indices] for k,v in events.items()}


def fit_usable(fit):
    return bool(support.usable(fit) and
                all(math.isfinite(fit['errors'][p]) and fit['errors'][p]>0
                    and math.isfinite(fit['values'][p]) for p in prod.PHYSICS_PARAMETERS))


def truth_in_fit_bounds(truth, true_dilutions, records, b):
    for p in prod.PHYSICS_PARAMETERS:
        lo,hi = prod.PARAMETER_LIMITS[p]
        if not lo <= truth[p] <= hi:
            return False
    for period in prod.PERIODS:
        record = records[(period,b)]
        width = max(8*record.stat_uncertainty,.20*record.value,.02)
        if not max(1e-6,record.value-width) <= true_dilutions[period] <= record.value+width:
            return False
    return True


def evaluate_recovery(toy, states, records, b, fit, truth, true_dilutions, profile_coverage):
    quality = fit_usable(fit)
    nll,_ = prod.make_bin_nll(toy,states,records,b,'nominal')
    feasible = truth_in_fit_bounds(truth,true_dilutions,records,b)
    truth_nll = float(nll(**truth,**{f'f_{p}':true_dilutions[p] for p in prod.PERIODS}))
    difference = float(fit['minimum_nll']-truth_nll)
    # A minimum above an explicitly feasible injected parameter point is suspect.
    above = bool(feasible and math.isfinite(difference) and difference>1e-3)
    result = {'quality_usable':quality,'usable':quality and not above,
              'truth_within_fit_bounds':feasible,'fit_nll_minus_truth_nll':difference,
              'selected_minimum_above_truth':above,'profiles':{}}
    if result['usable'] and profile_coverage:
        for parameter in TARGETS:
            try:
                constrained = support.constrained_fit(toy,states,records,b,fit,
                                                      {parameter:truth[parameter]}, warm=fit_start_values(fit))
                if constrained is None:
                    result['profiles'][parameter] = {'status':'failed','reason':'No valid constrained minimum'}
                else:
                    result['profiles'][parameter] = {'status':'evaluated',
                        'two_delta_nll':max(0.,2*(constrained['minimum_nll']-fit['minimum_nll']))}
            except (ValueError,RuntimeError) as exc:
                result['profiles'][parameter] = {'status':'failed','reason':str(exc)}
    return result


def init_worker(cache, states, records, baselines, config, output, signature):
    global _WORK
    all_events = prod.load_event_cache(Path(cache))
    # Keep only requested bins in each worker, then cache those slices.
    selected = np.isin(all_events['bin_number'],config['bins'])
    all_events = subset_events(all_events,np.flatnonzero(selected))
    _WORK = (all_events,states,records,baselines,config,Path(output),signature,{})


def cell_path(output, task):
    b,mode,i,j = task
    return output/'checkpoints'/f'bin_{b:02d}'/mode/f'u1_{i:02d}_ll1_{j:02d}'


def run_cell(task):
    all_events,states,records,baselines,config,output,signature,slices = _WORK
    b,mode,i,j = task
    folder = cell_path(output,task)
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/'owner.lock').open('a+') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            return {'task':list(task),'status':'owned_elsewhere'}
        state_path = folder/'state.json'
        if state_path.exists():
            state = json.loads(state_path.read_text())
            if state['signature']!=signature or state['task']!=list(task):
                raise RuntimeError(f'Checkpoint identity mismatch: {folder}')
            if state.get('complete'):
                return {'task':list(task),'status':state['status'],'reused':True}
        if b not in slices:
            slices[b] = subset_events(all_events,np.flatnonzero(all_events['bin_number']==b))
        events = slices[b]
        nominal = baselines[b]
        if len(events['phi']) != nominal['metadata']['number_of_events']:
            raise RuntimeError(f'Bin {b}: event count differs from baseline nominal fit')
        truth = make_truth(nominal,config['u1_values'][i],config['ll1_values'][j])
        true_dilutions = {p:float(nominal['values'][f'f_{p}']) for p in prod.PERIODS}
        start = fit_start_values(nominal)
        state = {'signature':signature,'task':list(task),'bin_number':b,'mode':mode,
                 'u1_index':i,'ll1_index':j,'truth':truth,'true_dilutions':true_dilutions,
                 'start_values':start,'requested':config['replicas'],
                 'original_events':len(events['phi']),'complete':False}
        # Check physical positivity over the full observed support, even for thinning modes.
        try:
            envelopes = support.generator_envelopes(events,states,truth,true_dilutions)
        except ValueError as exc:
            state.update(status='nonphysical_injection',reason=str(exc),complete=True)
            support.atomic_json(state_path,state)
            return {'task':list(task),'status':state['status']}
        try:
            sample = subset_events(events,coverage_indices(events,mode,config,b,0))
        except ValueError as exc:
            state.update(status='coverage_unavailable',reason=str(exc),complete=True)
            support.atomic_json(state_path,state)
            return {'task':list(task),'status':state['status']}
        state['selected_events'] = len(sample['phi'])
        state['selected_events_by_period'] = {p:int(np.sum(sample['period_index']==prod.PERIOD_INDEX[p])) for p in prod.PERIODS}
        state['phi_histogram'] = np.histogram(np.mod(sample['phi'],2*np.pi),bins=np.linspace(0,2*np.pi,37))[0].tolist()
        state['phi_histogram_note'] = 'Observed coverage, or one reproducible thinning realization (replica 0); masks vary across replicas'
        state['status'] = 'running'
        support.atomic_json(state_path,state)
        for rep in range(config['replicas']):
            path = folder/f'replica_{rep:05d}.json'
            if path.exists():
                old = json.loads(path.read_text())
                if old['signature']!=signature or old['replica']!=rep or old['task']!=list(task):
                    raise RuntimeError(f'Replica identity mismatch: {path}')
                continue
            seeds = [config['seed'],b,i,j,rep,555]
            row = {'signature':signature,'task':list(task),'replica':rep,
                   'spin_seed_components':seeds,'usable':False,'profiles':{}}
            try:
                selected = coverage_indices(events,mode,config,b,rep)
                sample = subset_events(events,selected)
                # Reuse the full envelope for observed coverage; recompute on thinned samples.
                toy_envelopes = envelopes if mode=='observed' else support.generator_envelopes(sample,states,truth,true_dilutions)
                rng = np.random.default_rng(np.random.SeedSequence(seeds))
                toy = support.generate_states(sample,states,truth,true_dilutions,rng,toy_envelopes)
                aux_rng = np.random.default_rng(np.random.SeedSequence([config['seed'],b,i,j,rep,777]))
                aux,centers = support.auxiliary_records(records,b,true_dilutions,aux_rng,'fluctuate')
                row['auxiliary_centers'] = centers
                fit = prod.fit_one_variant(toy,states,aux,b,'nominal',initial_values=start)
                row.update(status='fit',fit=support.compact_fit(fit))
                row.update(evaluate_recovery(toy,states,aux,b,fit,truth,true_dilutions,config['profile_coverage']))
            except (ValueError,RuntimeError,FloatingPointError) as exc:
                row.update(status='failed',reason=str(exc))
            support.atomic_json(path,row)
            if (rep+1)%10==0 or rep+1==config['replicas']:
                print(f'[bin {b} {mode} u1={truth["u1"]:+.3f} ll1={truth["ll1"]:+.3f}] {rep+1}/{config["replicas"]}',flush=True)
        state.update(status='complete',complete=True)
        support.atomic_json(state_path,state)
        return {'task':list(task),'status':'complete'}


def summarize_cell(state, rows, nominal):
    summaries=[]
    good = [r for r in rows if r.get('usable')]
    pair_values = np.array([[r['fit']['values'][p] for p in TARGETS] for r in good])
    empirical_corr = (float(np.corrcoef(pair_values.T)[0,1]) if len(good)>1
                      and np.all(np.std(pair_values,axis=0)>0) else None)
    fitted_corr=[]
    for r in good:
        fit=r['fit'];order=fit['parameter_order']
        corr=fit['correlation']
        if corr is not None:
            fitted_corr.append(corr[order.index('u1')][order.index('ll1')])
    for p in prod.PHYSICS_PARAMETERS:
        truth=state['truth'][p]
        estimates=np.array([r['fit']['values'][p] for r in good]);errors=np.array([r['fit']['errors'][p] for r in good])
        delta=estimates-truth;z=delta/errors
        profile=[r.get('profiles',{}).get(p,{}) for r in good]
        profile=[v['two_delta_nll'] for v in profile if v.get('status')=='evaluated']
        n=len(good);attempted=len(rows)
        summary={'bin_number':state['bin_number'],'mode':state['mode'],
                 'u1_index':state['u1_index'],'ll1_index':state['ll1_index'],
                 'injected_u1':state['truth']['u1'],'injected_ll1':state['truth']['ll1'],
                 'parameter':p,'truth':truth,'cell_status':state['status'],
                 'requested':state['requested'],'attempted':attempted,'usable':n,
                 'failed_or_unusable':attempted-n,
                 'usable_fraction':n/attempted if attempted else None,
                 'selected_minimum_above_truth':sum(bool(r.get('selected_minimum_above_truth')) for r in rows),
                 'selected_events':state.get('selected_events'),
                 'mean_recovered':float(estimates.mean()) if n else None,
                 'mean_bias':float(delta.mean()) if n else None,
                 'bias_over_data_stat':float(delta.mean()/nominal['errors'][p]) if n else None,
                 'bias_mc_error':float(delta.std(ddof=1)/np.sqrt(n)) if n>1 else None,
                 'empirical_spread':float(estimates.std(ddof=1)) if n>1 else None,
                 'median_reported_sigma':float(np.median(errors)) if n else None,
                 'pull_mean':float(z.mean()) if n else None,
                 'pull_width':float(z.std(ddof=1)) if n>1 else None,
                 'pull_rms':float(np.sqrt(np.mean(z*z))) if n else None,
                 'empirical_u1_ll1_correlation':empirical_corr,
                 'median_fitted_u1_ll1_correlation':float(np.median(fitted_corr)) if fitted_corr else None,
                 'profile_evaluated':len(profile),
                 'profile_failed':sum(r.get('profiles',{}).get(p,{}).get('status')=='failed' for r in good),
                 'profile_found_lower_minimum':sum('lower minimum' in r.get('profiles',{}).get(p,{}).get('reason','') for r in good)}
        for label,threshold in [('68',1.),('95',1.96)]:
            covered=int(np.sum(np.abs(z)<=threshold))
            summary['hesse_coverage'+label+'_usable']=covered/n if n else None
            summary['hesse_coverage'+label+'_all_attempts']=covered/attempted if attempted else None
            expected=.682689 if label=='68' else .95
            summary['hesse_coverage'+label+'_mc_se_at_nominal']=math.sqrt(expected*(1-expected)/n) if n else None
        for label,level in [('68',1.),('95',3.84)]:
            covered=int(np.sum(np.asarray(profile)<=level))
            summary['profile_coverage'+label+'_evaluated']=covered/len(profile) if profile else None
            summary['profile_coverage'+label+'_all_attempts']=covered/attempted if attempted and p in TARGETS else None
        summaries.append(summary)
    return summaries


def load_results(output, manifest):
    states=[];long_rows=[];summary=[];paired={}
    for task in tasks_for(manifest['config']):
        folder=cell_path(output,task);path=folder/'state.json'
        if not path.exists():
            continue
        state=json.loads(path.read_text())
        if state['signature']!=manifest['signature']:
            raise RuntimeError(f'Stale cell {path}')
        rows=[]
        for p in sorted(folder.glob('replica_*.json')):
            row=json.loads(p.read_text())
            if row['signature']!=manifest['signature'] or row['task']!=state['task']:
                raise RuntimeError(f'Stale replica {p}')
            rows.append(row)
            for name in prod.PHYSICS_PARAMETERS:
                fit=row.get('fit',{});profile=row.get('profiles',{}).get(name,{})
                long_rows.append({'bin_number':state['bin_number'],'mode':state['mode'],
                    'u1_index':state['u1_index'],'ll1_index':state['ll1_index'],
                    'parameter':name,'replica':row['replica'],'truth':state['truth'][name],
                    'estimate':fit.get('values',{}).get(name),'sigma':fit.get('errors',{}).get(name),
                    'usable':row.get('usable'),'status':row['status'],
                    'selected_minimum_above_truth':row.get('selected_minimum_above_truth'),
                    'profile_two_delta_nll':profile.get('two_delta_nll'),'profile_status':profile.get('status')})
        states.append({**state,'attempted':len(rows),'usable':sum(r.get('usable',False) for r in rows)})
        summary.extend(summarize_cell(state,rows,manifest['baselines'][str(state['bin_number'])]))
        paired[tuple(state['task'])]=rows
    comparison=[]
    config=manifest['config']
    if 'matched_forward_gap' in config['modes']:
        for b in config['bins']:
            for i in range(len(config['u1_values'])):
                for j in range(len(config['ll1_values'])):
                    a={r['replica']:r for r in paired.get((b,'matched_uniform',i,j),[]) if r.get('usable')}
                    z={r['replica']:r for r in paired.get((b,'matched_forward_gap',i,j),[]) if r.get('usable')}
                    shared=sorted(set(a)&set(z))
                    for p in prod.PHYSICS_PARAMETERS:
                        delta=np.array([z[k]['fit']['values'][p]-a[k]['fit']['values'][p] for k in shared])
                        comparison.append({'bin_number':b,'u1_index':i,'ll1_index':j,'parameter':p,
                            'injected_u1':config['u1_values'][i],'injected_ll1':config['ll1_values'][j],
                            'paired_usable':len(shared),'uniform_usable':len(a),'gap_usable':len(z),
                            'mean_gap_minus_uniform':float(delta.mean()) if len(shared) else None,
                            'paired_mc_error':float(delta.std(ddof=1)/np.sqrt(len(delta))) if len(delta)>1 else None})
    return states,pd.DataFrame(summary),pd.DataFrame(long_rows),pd.DataFrame(comparison)


def tasks_for(config):
    return [(b,mode,i,j) for b in reversed(config['bins']) for mode in config['modes']
            for i in range(len(config['u1_values'])) for j in range(len(config['ll1_values']))]


def heatmap(ax,frame,key,config,title,low=None,high=None,cmap='viridis'):
    values=np.full((len(config['ll1_values']),len(config['u1_values'])),np.nan)
    for _,row in frame.iterrows():
        value=row.get(key)
        if value is not None and np.isfinite(value):
            values[int(row.ll1_index),int(row.u1_index)]=value
    cm=plt.get_cmap(cmap).copy();cm.set_bad('.85')
    image=ax.imshow(np.ma.masked_invalid(values),origin='lower',aspect='auto',cmap=cm,vmin=low,vmax=high)
    ax.set_xticks(range(len(config['u1_values'])),[f'{v:g}' for v in config['u1_values']])
    ax.set_yticks(range(len(config['ll1_values'])),[f'{v:g}' for v in config['ll1_values']])
    ax.set(xlabel=r'Injected $F_{UU}^{\cos\phi}/F_{UU}$',ylabel=r'Injected $F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$',title=title)
    plt.colorbar(image,ax=ax,fraction=.046,pad=.04)


def write_products(output):
    manifest=json.loads((output/'manifest.json').read_text());config=manifest['config']
    states,frame,toys,comparison=load_results(output,manifest)
    if not states:
        raise RuntimeError('No cell checkpoints available yet')
    tables=output/'tables';plots=output/'plots';tables.mkdir(exist_ok=True);plots.mkdir(exist_ok=True)
    support.atomic_json(tables/'cell_status.json',states)
    pd.DataFrame([{k:v for k,v in s.items() if k not in ('signature','truth','true_dilutions','start_values','phi_histogram')} for s in states]).to_csv(tables/'cell_status.csv',index=False)
    frame.to_csv(tables/'injection_summary.csv',index=False);toys.to_csv(tables/'toy_results.csv',index=False)
    if len(comparison):
        comparison.to_csv(tables/'matched_coverage_comparison.csv',index=False)
    labels={'u1':r'$F_{UU}^{\cos\phi}/F_{UU}$','ll1':r'$F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$'}
    for b in config['bins']:
        available=[s for s in states if s['bin_number']==b and 'phi_histogram' in s]
        if available:
            fig,ax=plt.subplots(figsize=(8,4))
            for mode in config['modes']:
                cell=next((s for s in available if s['mode']==mode),None)
                if cell:
                    ax.stairs(cell['phi_histogram'],np.linspace(0,360,37),label=f'{mode}: N={cell["selected_events"]}')
            ax.set(xlabel=r'$\phi$ (degrees)',ylabel='Events',title=f'Bin {b}: observed coverage / replica-0 thinning');ax.legend(fontsize=8)
            fig.tight_layout();fig.savefig(plots/f'phi_coverage_bin_{b:02d}.png',dpi=170);plt.close(fig)
        for mode in config['modes']:
            group=frame[(frame.bin_number==b)&(frame['mode']==mode)]
            if group.empty:
                continue
            fig,axes=plt.subplots(1,2,figsize=(12,4.6))
            for ax,p in zip(axes,TARGETS):
                part=group[(group.parameter==p)&group.mean_recovered.notna()]
                other='injected_ll1' if p=='u1' else 'injected_u1';own='injected_'+p
                for value,sub in part.groupby(other):
                    sub=sub.sort_values(own)
                    ax.errorbar(sub[own],sub.mean_recovered,yerr=sub.bias_mc_error,fmt='o-',ms=4,
                                label=f'Other cosine = {value:g}')
                limits=config[p+'_values'];ax.plot([min(limits),max(limits)],[min(limits),max(limits)],'k--',lw=1)
                ax.set(xlabel='Injected '+labels[p],ylabel='Mean recovered '+labels[p]);ax.legend(fontsize=7)
            fig.suptitle(f'Bin {b}, {mode}: recovery (bars = MC error of mean)');fig.tight_layout()
            fig.savefig(plots/f'recovery_bin_{b:02d}_{mode}.png',dpi=170);plt.close(fig)
            fig,axes=plt.subplots(1,2,figsize=(12,4.5))
            finite=group.loc[group.parameter.isin(TARGETS),'bias_over_data_stat'].dropna().abs()
            scale=max(float(finite.max()) if len(finite) else 0,.1)
            for ax,p in zip(axes,TARGETS):
                heatmap(ax,group[group.parameter==p],'bias_over_data_stat',config,labels[p]+' mean bias / data stat.',-scale,scale,'coolwarm')
            fig.suptitle(f'Bin {b}, {mode}; gray = missing/nonphysical/no usable fits');fig.tight_layout()
            fig.savefig(plots/f'bias_bin_{b:02d}_{mode}.png',dpi=170);plt.close(fig)
            fig,axes=plt.subplots(2,4,figsize=(18,8))
            metrics=[('hesse_coverage68_usable','Hesse 68%'),('hesse_coverage95_usable','Hesse 95%'),
                     ('profile_coverage68_evaluated','Profile LR 68%'),('profile_coverage95_evaluated','Profile LR 95%')]
            for row,p in zip(axes,TARGETS):
                for ax,(key,title) in zip(row,metrics):
                    heatmap(ax,group[group.parameter==p],key,config,title+' '+labels[p],0,1)
            fig.suptitle(f'Bin {b}, {mode}: usable/evaluated-fit coverage; counts in CSV');fig.tight_layout()
            fig.savefig(plots/f'coverage_bin_{b:02d}_{mode}.png',dpi=160);plt.close(fig)
            fig,axes=plt.subplots(1,3,figsize=(17,4.7));part=group[group.parameter=='u1']
            for ax,key,title,lo,hi in zip(axes,('usable_fraction','empirical_u1_ll1_correlation','median_fitted_u1_ll1_correlation'),
                ('Usable / attempted','Empirical cosine correlation','Median fitted cosine correlation'),(0,-1,-1),(1,1,1)):
                heatmap(ax,part,key,config,title,lo,hi,'viridis' if lo==0 else 'coolwarm')
            fig.suptitle(f'Bin {b}, {mode}: fit quality and cross-talk');fig.tight_layout()
            fig.savefig(plots/f'fit_quality_bin_{b:02d}_{mode}.png',dpi=170);plt.close(fig)
        if len(comparison):
            part=comparison[comparison.bin_number==b].copy()
            if len(part):
                fig,axes=plt.subplots(1,2,figsize=(12,4.6))
                for ax,p in zip(axes,TARGETS):
                    heatmap(ax,part[part.parameter==p],'mean_gap_minus_uniform',config,labels[p]+': gap minus uniform',cmap='coolwarm')
                fig.suptitle(f'Bin {b}: paired recovery difference at equal counts; MC errors in CSV');fig.tight_layout()
                fig.savefig(plots/f'matched_coverage_difference_bin_{b:02d}.png',dpi=170);plt.close(fig)
    print(f'Summary: {len(states)}/{len(tasks_for(config))} cell states; {sum(s["status"]=="complete" for s in states)} finished; {len(toys)//len(prod.PHYSICS_PARAMETERS)} attempted toys',flush=True)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bins',type=int,nargs='+',default=[4,5,6])
    p.add_argument('--u1-values',type=float,nargs='+',default=[-.5,-.25,0,.25,.5])
    p.add_argument('--ll1-values',type=float,nargs='+',default=[-.5,-.25,0,.25,.5])
    p.add_argument('--replicas',type=int,default=100)
    p.add_argument('--workers',type=int,default=10)
    p.add_argument('--seed',type=int,default=20261120)
    p.add_argument('--coverage-comparison',action='store_true',help='Also compare uniform thinning and a forward gap at matched counts')
    p.add_argument('--retained-fraction',type=float,default=.7)
    p.add_argument('--gap-half-width',type=float,default=30.,help='Remove phi within this many degrees of 0/360 in gap control')
    p.add_argument('--no-profile-coverage',action='store_false',dest='profile_coverage',help='Skip true-value profile-LR coverage for a faster Hesse-only scan')
    p.set_defaults(profile_coverage=True)
    p.add_argument('--distributed',action='store_true',help='Use a shared cell queue on two nodes; summarize separately')
    p.add_argument('--summarize-only',action='store_true')
    p.add_argument('--output-dir',type=Path,default=HERE/'output/cosine_coverage_validation')
    p.add_argument('--baseline-dir',type=Path,default=HERE/'output/unpolarized_validation')
    p.add_argument('--production-output',type=Path,default=HERE/'output/asymmetry_extraction')
    p.add_argument('--cache',type=Path)
    p.add_argument('--rebuild-events',action='store_true')
    p.add_argument('--run-info-csv',type=Path,default=HERE/prod.DEFAULT_RUN_INFO_CSV)
    p.add_argument('--cut-json',type=Path,default=HERE/prod.DEFAULT_CUT_JSON)
    p.add_argument('--dilution-json',type=Path)
    p.add_argument('--dilution-dir',type=Path,default=HERE/prod.DEFAULT_DILUTION_DIR)
    p.add_argument('--input',type=prod.parse_input_override,action='append',default=[])
    p.add_argument('--tree',default=prod.DEFAULT_TREE_NAME)
    p.add_argument('--chunk-size',default=prod.DEFAULT_CHUNK_SIZE)
    return p


def prepare(args):
    bins=sorted(set(args.bins));u1=sorted(set(args.u1_values));ll1=sorted(set(args.ll1_values))
    if not bins or any(b not in range(1,25) for b in bins) or args.replicas<2 or args.workers<1 or args.seed<0:
        raise ValueError('Valid bins, at least two replicas and positive workers required')
    for p,values in zip(TARGETS,(u1,ll1)):
        lo,hi=prod.PARAMETER_LIMITS[p]
        if not values or any(not math.isfinite(v) or not lo<=v<=hi for v in values):
            raise ValueError(f'Grid values for {p} must lie in [{lo},{hi}]')
    if not 0<args.retained_fraction<1 or not 0<args.gap_half_width<180:
        raise ValueError('Require 0<retained_fraction<1 and 0<gap_half_width<180')
    output=args.output_dir.expanduser().resolve();baseline=args.baseline_dir.expanduser().resolve()
    production=args.production_output.expanduser().resolve()
    for protected in (baseline,production):
        if output==protected or output in protected.parents or protected in output.parents:
            raise ValueError('Choose a separate output directory outside existing validation and production results')
    records=prod.parse_run_info_csv(args.run_info_csv.expanduser().resolve());states=prod.run_state_arrays(records)
    cuts_path=args.cut_json.expanduser().resolve();cuts=prod.load_channel_cuts(cuts_path,cut_label='nominal')
    dilution_path=(args.dilution_json or prod.find_default_dilution_json(args.dilution_dir.expanduser().resolve())).resolve()
    prod.validate_dilution_cut_source(dilution_path,cuts_path);dilutions=prod.load_dilution_factors(dilution_path,cut_label='nominal')
    inputs={p:Path(v).expanduser().resolve() for p,v in prod.DEFAULT_INPUTS.items()}
    inputs.update({p:v.expanduser().resolve() for p,v in args.input})
    identity=prod._cache_identity(inputs,args.tree,records)
    old=json.loads((baseline/'manifest.json').read_text())
    pair=(old['production_sha256'],support.digest(Path(prod.__file__)))
    if pair[0]!=pair[1] and pair!=support.MIGRATION_ONLY_UPGRADE:
        raise RuntimeError('Baseline production code differs from current fitter')
    if (not support.same_event_identity(old['event_identity'],identity)
        or old['cuts_sha256']!=support.digest(cuts_path)
        or old['dilution_sha256']!=support.digest(dilution_path)
        or old['run_info_sha256']!=support.digest(args.run_info_csv.expanduser().resolve())):
        raise RuntimeError('Baseline selected events, cuts, run information or dilution inputs differ')
    baselines={}
    for b in bins:
        state=json.loads((baseline/'checkpoints'/f'bin_{b:02d}.json').read_text())
        if state['bin_number']!=b or not fit_usable(state['nominal']):
            raise RuntimeError(f'Bin {b}: baseline nominal fit is unusable')
        baselines[b]=state['nominal']
    cache=args.cache.expanduser().resolve() if args.cache else production/'nominal/cache/selected_events.npz'
    actual=prod._event_cache_path(cache);actual.parent.mkdir(parents=True,exist_ok=True)
    with actual.with_suffix(actual.suffix+'.build.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if args.rebuild_events:
            if args.distributed:
                raise ValueError('Do not rebuild event caches during distributed runs')
            prod.build_event_cache(inputs,args.tree,args.chunk_size,records,cuts,cache)
        elif not actual.exists():
            prod.build_event_cache(inputs,args.tree,args.chunk_size,records,cuts,cache)
        support.validate_resume_cache(cache,cuts,identity,True)
    config={'bins':bins,'u1_values':u1,'ll1_values':ll1,'replicas':args.replicas,'seed':args.seed,
            'modes':list(MODES) if args.coverage_comparison else ['observed'],
            'retained_fraction':args.retained_fraction,'gap_half_width':args.gap_half_width,
            'profile_coverage':args.profile_coverage}
    provenance={'config':config,'runtime':{'python':sys.version,'numpy':np.__version__,'iminuit':version('iminuit')},
                'study_sha256':support.digest(Path(__file__)),'support_sha256':support.digest(Path(support.__file__)),
                'production_sha256':pair[1],'event_identity':identity,'cuts_sha256':support.digest(cuts_path),
                'dilution_sha256':support.digest(dilution_path),'run_info_sha256':support.digest(args.run_info_csv.expanduser().resolve()),
                'baseline_manifest_sha256':support.digest(baseline/'manifest.json'),
                'baselines':{str(b):fit for b,fit in baselines.items()}}
    import hashlib
    signature=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
    output.mkdir(parents=True,exist_ok=True)
    with (output/'.manifest.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);path=output/'manifest.json'
        if path.exists():
            if json.loads(path.read_text())['signature']!=signature:
                raise RuntimeError('Inputs/configuration/code changed; use a new output directory. Workers/distributed flags may change.')
        else:
            if any(p.name!='.manifest.lock' for p in output.iterdir()):
                raise RuntimeError('Nonempty output without manifest')
            support.atomic_json(path,{'signature':signature,**provenance,
                'scope':'conditional recovery under observed coverage; optional matched thinning; not detector simulation',
                'profile_coverage_definition':'profile LR at injected coefficient, thresholds 1 and 3.84; no full interval scan',
                'auxiliary_mode':'fluctuate'})
    return output,cache,states,dilutions,baselines,config,signature


def main():
    args=parser().parse_args()
    if args.summarize_only:
        write_products(args.output_dir.expanduser().resolve());return 0
    output,cache,states,records,baselines,config,signature=prepare(args)
    tasks=tasks_for(config);status=[];failures=[];started=time.monotonic()
    print(f'{len(tasks)} grid cells; {config["replicas"]} attempts per physical cell; workers={min(args.workers,10)}; profile coverage={config["profile_coverage"]}',flush=True)
    with ProcessPoolExecutor(max_workers=min(args.workers,10,len(tasks)),mp_context=mp.get_context('spawn'),
            initializer=init_worker,initargs=(str(cache),states,records,baselines,config,str(output),signature)) as pool:
        futures={pool.submit(run_cell,task):task for task in tasks}
        for future in as_completed(futures):
            task=futures[future]
            try:
                result=future.result();status.append(result)
                if result['status']!='owned_elsewhere':
                    print(f'Cell {task}: {result["status"]}; elapsed {(time.monotonic()-started)/60:.1f} min',flush=True)
            except Exception as exc:
                failures.append({'task':list(task),'reason':repr(exc)})
                print(f'Cell {task}: FAILED {exc}',flush=True)
    support.atomic_json(output/'node_status'/f'{socket.gethostname()}_{os.getpid()}.json',{'results':status,'failures':failures})
    if not args.distributed:
        write_products(output)
    print(f'Output: {output}',flush=True)
    return 2 if failures else 0


if __name__=='__main__':raise SystemExit(main())
