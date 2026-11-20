#!/usr/bin/env python3
"""Conditional u1/ll1 line injections through each bin's fitted coefficients."""
from __future__ import annotations
import os
for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(variable,'1')
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
import hashlib
from importlib.metadata import version
import json
import math
import multiprocessing as mp
from pathlib import Path
import socket
import sys
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import validate_cosine_coverage as core
support=core.support
prod=core.prod
HERE=Path(__file__).resolve().parent
DEFAULT_BINS=[4,5,6,10,11,12,16,17,18,22,23,24]
LABELS={'u1':r'$F_{UU}^{\cos\phi}/F_{UU}$','ll1':r'$F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$'}
_SETTINGS=None


def grid_for(nominal,u_values,l_values):
    # Include the actual measured pair once; share it between the two line scans.
    u=sorted(set([*u_values,float(nominal['values']['u1'])]))
    l=sorted(set([*l_values,float(nominal['values']['ll1'])]))
    return {'u1_values':u,'ll1_values':l,
            'u1_anchor':u.index(float(nominal['values']['u1'])),
            'll1_anchor':l.index(float(nominal['values']['ll1']))}


def tasks_for(config):
    # Process measured pairs first, then interleave bins, working outwards from each fit.
    by_bin={}
    for b in reversed(config['bins']):
        g=config['grids'][str(b)];a,c=g['u1_anchor'],g['ll1_anchor']
        u=sorted((i for i in range(len(g['u1_values'])) if i!=a),
                 key=lambda i:abs(g['u1_values'][i]-g['u1_values'][a]))
        l=sorted((j for j in range(len(g['ll1_values'])) if j!=c),
                 key=lambda j:abs(g['ll1_values'][j]-g['ll1_values'][c]))
        cells=[(a,c)]
        for n in range(max(len(u),len(l))):
            if n<len(u):cells.append((u[n],c))
            if n<len(l):cells.append((a,l[n]))
        by_bin[b]=[(b,mode,i,j) for i,j in cells for mode in config['modes']]
    return [by_bin[b][i] for i in range(max(map(len,by_bin.values())))
            for b in by_bin if i<len(by_bin[b])]


def scan_label(config,task):
    b,_,i,j=task;g=config['grids'][str(b)]
    if i==g['u1_anchor'] and j==g['ll1_anchor']:return 'nominal'
    return 'u1' if j==g['ll1_anchor'] else 'll1'


def init_worker(cache,states,records,baselines,config,output,signature):
    global _SETTINGS
    _SETTINGS=config
    core.init_worker(cache,states,records,baselines,config,output,signature)


def run_cell(task):
    # Reuse the exact generator/fitter/checkpoint code from the original study.
    events,states,records,baselines,_,output,signature,slices=core._WORK
    per_bin={**_SETTINGS,**_SETTINGS['grids'][str(task[0])]}
    core._WORK=(events,states,records,baselines,per_bin,output,signature,slices)
    return core.run_cell(task)


def load_results(output,manifest):
    config=manifest['config'];states=[];summaries=[];toys=[]
    for task in tasks_for(config):
        folder=core.cell_path(output,task);state_path=folder/'state.json'
        if not state_path.exists():continue
        state=json.loads(state_path.read_text())
        if state['signature']!=manifest['signature'] or state['task']!=list(task):
            raise RuntimeError(f'Checkpoint identity mismatch: {state_path}')
        rows=[]
        for path in sorted(folder.glob('replica_*.json')):
            row=json.loads(path.read_text())
            if row['signature']!=manifest['signature'] or row['task']!=list(task):
                raise RuntimeError(f'Replica identity mismatch: {path}')
            rows.append(row)
            for p in prod.PHYSICS_PARAMETERS:
                fit=row.get('fit',{});profile=row.get('profiles',{}).get(p,{})
                toys.append({'bin_number':task[0],'mode':task[1],'scan':scan_label(config,task),
                    'u1_index':task[2],'ll1_index':task[3],'injected_u1':state['truth']['u1'],
                    'injected_ll1':state['truth']['ll1'],'replica':row['replica'],'parameter':p,
                    'truth':state['truth'][p],'estimate':fit.get('values',{}).get(p),
                    'sigma':fit.get('errors',{}).get(p),'usable':row.get('usable',False),
                    'status':row['status'],'profile_status':profile.get('status'),
                    'profile_two_delta_nll':profile.get('two_delta_nll'),
                    'selected_minimum_above_truth':row.get('selected_minimum_above_truth')})
        label=scan_label(config,task)
        states.append({**state,'scan':label,'attempted':len(rows),'usable':sum(r.get('usable',False) for r in rows)})
        for summary in core.summarize_cell(state,rows,manifest['baselines'][str(task[0])]):
            summaries.append({**summary,'scan':label})
    frame=pd.DataFrame(summaries)
    # Entirely skipped/failed cells still need plottable numeric NaNs.
    for name in ('mean_recovered','bias_mc_error','bias_over_data_stat','hesse_coverage68_usable',
                 'profile_coverage68_evaluated','profile_evaluated','usable','usable_fraction'):
        if name in frame:frame[name]=pd.to_numeric(frame[name],errors='coerce')
    return states,frame,pd.DataFrame(toys)


def write_products(output):
    manifest=json.loads((output/'manifest.json').read_text());config=manifest['config']
    states,frame,toys=load_results(output,manifest)
    if not states:raise RuntimeError('No cell checkpoints available yet')
    tables=output/'tables';plots=output/'plots';tables.mkdir(exist_ok=True);plots.mkdir(exist_ok=True)
    support.atomic_json(tables/'cell_status.json',states)
    pd.DataFrame([{k:v for k,v in s.items() if k not in ('signature','truth','true_dilutions','start_values','phi_histogram')}
                  for s in states]).to_csv(tables/'cell_status.csv',index=False)
    frame.to_csv(tables/'injection_summary.csv',index=False)
    toys.to_csv(tables/'toy_results.csv',index=False)
    if not frame.empty:
        for b in config['bins']:
            for mode in config['modes']:
                group=frame[(frame.bin_number==b)&(frame['mode']==mode)]
                if group.empty:continue
                # Columns: recovered coefficients; rows: which coefficient is scanned.
                fig,axes=plt.subplots(2,2,figsize=(12,8))
                for r,axis in enumerate(core.TARGETS):
                    part=group[group['scan'].isin([axis,'nominal'])]
                    for col,p in enumerate(core.TARGETS):
                        z=part[part.parameter==p].sort_values('injected_'+axis)
                        ax=axes[r,col];valid=z[z.mean_recovered.notna()]
                        ax.plot(z['injected_'+axis],z.truth,'k--',lw=1,label='Injected value')
                        ax.errorbar(valid['injected_'+axis],valid.mean_recovered,yerr=valid.bias_mc_error,
                                    fmt='o-',color='black',label='Mean recovered (MC error of mean)')
                        ax.set(xlabel='Injected '+LABELS[axis],ylabel='Mean recovered '+LABELS[p])
                        ax.legend(fontsize=7)
                fig.suptitle(f'Bin {b}, {mode}: one coefficient scanned, other injected at fitted value; all seven refitted')
                fig.tight_layout();fig.savefig(plots/f'recovery_bin_{b:02d}_{mode}.png',dpi=170);plt.close(fig)
                for kind in ('bias','coverage','fit_quality'):
                    fig,axes=plt.subplots(1,2,figsize=(13,4.6))
                    for ax,axis in zip(axes,core.TARGETS):
                        part=group[group['scan'].isin([axis,'nominal'])]
                        for p in core.TARGETS:
                            z=part[part.parameter==p].sort_values('injected_'+axis)
                            x=z['injected_'+axis];color='tab:blue' if p=='u1' else 'tab:orange'
                            if kind=='bias':
                                ax.errorbar(x,z.bias_over_data_stat,
                                    yerr=z.bias_mc_error/manifest['baselines'][str(b)]['errors'][p],
                                    fmt='o-',color=color,label=LABELS[p])
                            elif kind=='coverage':
                                for key,style,label in [('hesse_coverage68_usable','o-','Hesse'),
                                         ('profile_coverage68_evaluated','s--','Profile LR')]:
                                    ax.plot(x,z[key],style,color=color,label=LABELS[p]+' '+label)
                            else:
                                if p=='u1':ax.plot(x,z.usable_fraction,'ko-',label='Usable / attempted')
                                fraction=np.divide(z.profile_evaluated,z.usable,
                                    out=np.full(len(z),np.nan),where=z.usable.to_numpy()>0)
                                ax.plot(x,fraction,'s--',color=color,label=LABELS[p]+' profile evaluated / usable')
                        if kind=='bias':ax.axhline(0,color='.5',ls='--');ax.set_ylabel('Mean bias / original data statistical uncertainty')
                        elif kind=='coverage':
                            ax.axhline(.682689,color='.5',ls=':');ax.set(ylim=(0,1.03),ylabel='68% coverage (usable/evaluated trials)')
                        else:ax.set(ylim=(0,1.03),ylabel='Fraction')
                        ax.set_xlabel('Injected '+LABELS[axis]);ax.legend(fontsize=7)
                    fig.suptitle(f'Bin {b}, {mode}: {kind}; counts and 95% coverage in CSV')
                    fig.tight_layout();fig.savefig(plots/f'{kind}_bin_{b:02d}_{mode}.png',dpi=170);plt.close(fig)
            cells=[s for s in states if s['bin_number']==b and 'phi_histogram' in s]
            if cells:
                fig,ax=plt.subplots(figsize=(8,4))
                for mode in config['modes']:
                    cell=next((s for s in cells if s['mode']==mode),None)
                    if cell:ax.stairs(cell['phi_histogram'],np.linspace(0,360,37),label=f'{mode}: N={cell["selected_events"]}')
                ax.set(xlabel=r'$\phi$ (degrees)',ylabel='Events',title=f'Bin {b}: accepted event coverage (replica 0 for thinning)')
                ax.legend();fig.tight_layout();fig.savefig(plots/f'phi_coverage_bin_{b:02d}.png',dpi=170);plt.close(fig)
    print(f'{len(states)}/{len(tasks_for(config))} cell states; {sum(s["status"]=="complete" for s in states)} finished; '
          f'{sum(s["attempted"] for s in states)} attempted trials',flush=True)


def parser():
    p=core.parser();p.description=__doc__
    p.set_defaults(bins=DEFAULT_BINS,output_dir=HERE/'output/cosine_axis_validation')
    return p


def prepare(args):
    bins=sorted(set(args.bins));u1=sorted(set(args.u1_values));ll1=sorted(set(args.ll1_values))
    if not bins or any(b not in range(1,25) for b in bins) or args.replicas<2 or args.workers<1 or args.seed<0:
        raise ValueError('Valid bins, at least two replicas and positive workers required')
    for p,values in zip(core.TARGETS,(u1,ll1)):
        lo,hi=prod.PARAMETER_LIMITS[p]
        if not values or any(not math.isfinite(v) or not lo<=v<=hi for v in values):
            raise ValueError(f'Grid values for {p} must lie in [{lo},{hi}]')
    if not 0<args.retained_fraction<1 or not 0<args.gap_half_width<180:
        raise ValueError('Require 0<retained_fraction<1 and 0<gap_half_width<180')
    output=args.output_dir.expanduser().resolve();baseline=args.baseline_dir.expanduser().resolve()
    production=args.production_output.expanduser().resolve()
    protected_paths=(baseline,production,HERE/'output/cosine_coverage_validation',HERE/'output/cosine_minimum_audit')
    for protected in protected_paths:
        protected=protected.resolve()
        if output==protected or output in protected.parents or protected in output.parents:
            raise ValueError('Choose a separate output directory outside previous validation and production results')
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
        if state['bin_number']!=b or not core.fit_usable(state['nominal']):
            raise RuntimeError(f'Bin {b}: baseline nominal fit is unusable')
        baselines[b]=state['nominal']
    cache=args.cache.expanduser().resolve() if args.cache else production/'nominal/cache/selected_events.npz'
    actual=prod._event_cache_path(cache);actual.parent.mkdir(parents=True,exist_ok=True)
    with actual.with_suffix(actual.suffix+'.build.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if args.rebuild_events:
            if args.distributed:raise ValueError('Do not rebuild event caches during distributed runs')
            prod.build_event_cache(inputs,args.tree,args.chunk_size,records,cuts,cache)
        elif not actual.exists():prod.build_event_cache(inputs,args.tree,args.chunk_size,records,cuts,cache)
        support.validate_resume_cache(cache,cuts,identity,True)
    config={'design':'two_lines_through_fitted_pair','bins':bins,'u1_values':u1,'ll1_values':ll1,
            'grids':{str(b):grid_for(fit,u1,ll1) for b,fit in baselines.items()},
            'replicas':args.replicas,'seed':args.seed,
            'modes':list(core.MODES) if args.coverage_comparison else ['observed'],
            'retained_fraction':args.retained_fraction,'gap_half_width':args.gap_half_width,
            'profile_coverage':args.profile_coverage}
    provenance={'config':config,'runtime':{'python':sys.version,'numpy':np.__version__,'iminuit':version('iminuit')},
                'study_sha256':support.digest(Path(__file__)),'core_sha256':support.digest(Path(core.__file__)),
                'support_sha256':support.digest(Path(support.__file__)),
                'production_sha256':pair[1],'event_identity':identity,'cuts_sha256':support.digest(cuts_path),
                'dilution_sha256':support.digest(dilution_path),'run_info_sha256':support.digest(args.run_info_csv.expanduser().resolve()),
                'baseline_manifest_sha256':support.digest(baseline/'manifest.json'),
                'baselines':{str(b):fit for b,fit in baselines.items()}}
    signature=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
    output.mkdir(parents=True,exist_ok=True)
    with (output/'.manifest.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);path=output/'manifest.json'
        if path.exists():
            if json.loads(path.read_text())['signature']!=signature:
                raise RuntimeError('Inputs/configuration/code changed; use a new output directory. Workers/distributed flags may change.')
        else:
            if any(p.name!='.manifest.lock' for p in output.iterdir()):raise RuntimeError('Nonempty output without manifest')
            support.atomic_json(path,{'signature':signature,**provenance,
                'scope':'conditional recovery on two injection lines through fitted pair; not detector simulation',
                'profile_coverage_definition':'profile LR at injected coefficient, thresholds 1 and 3.84; no full interval scan',
                'auxiliary_mode':'fluctuate'})
    return output,cache,states,dilutions,baselines,config,signature


def main():
    args=parser().parse_args()
    if args.summarize_only:write_products(args.output_dir.expanduser().resolve());return 0
    output,cache,states,records,baselines,config,signature=prepare(args)
    tasks=tasks_for(config);status=[];failures=[];started=time.monotonic()
    print(f'{len(config["bins"])} bins; {len(tasks)} injection cells including fitted pairs; '
          f'{config["replicas"]} attempts/cell; workers={min(args.workers,10)}; profiles={config["profile_coverage"]}',flush=True)
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
    if not args.distributed:write_products(output)
    print(f'Output: {output}',flush=True)
    return 2 if failures else 0

if __name__=='__main__':raise SystemExit(main())
