#!/usr/bin/env python3
"""Refit identical saved cosine toys; keep original checkpoints read-only."""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(key,'1')
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
import json
import math
import multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd
import validate_cosine_coverage as study
s=study.support
p=study.prod
_W=None


def starts_for(nominal, old, count, seed):
    # No injected truth, grid location, or neighboring toy fit enters these starts.
    base=dict(nominal['values'])
    starts=[dict(old['fit']['values'])] if 'fit' in old else []
    starts += [dict(base),{**base,**{k:0. for k in p.PHYSICS_PARAMETERS}}]
    rng=np.random.default_rng(seed)
    for _ in range(count):
        start=dict(base)
        for k in p.PHYSICS_PARAMETERS:
            lo,hi=p.PARAMETER_LIMITS[k]
            start[k]=float(rng.uniform(max(lo,-.8),min(hi,.8)))
        starts.append(start)
    return starts


def search(nll, starts, records, b, fixed=None, ncall=50000):
    """Choose by likelihood first; covariance quality must not hide a lower minimum."""
    candidates=[];attempts=[]
    for start in starts:
        values={k:float(start[k]) for k in p.PHYSICS_PARAMETERS}
        values.update({f'f_{q}':float(start.get(f'f_{q}',records[(q,b)].value)) for q in p.PERIODS})
        m=p.Minuit(nll,**values);m.errordef=p.Minuit.LIKELIHOOD;m.print_level=0;m.strategy=1
        for k in p.PHYSICS_PARAMETERS:m.limits[k]=p.PARAMETER_LIMITS[k]
        for q in p.PERIODS:
            r=records[(q,b)];width=max(8*r.stat_uncertainty,.2*r.value,.02)
            m.limits[f'f_{q}']=(max(1e-6,r.value-width),r.value+width)
            if r.stat_uncertainty==0:m.fixed[f'f_{q}']=True
        for k,v in (fixed or {}).items():m.values[k]=v;m.fixed[k]=True
        try:
            m.migrad(ncall=ncall)
            # An independent simplex step is used even for nominally valid fits.
            m.simplex(ncall=ncall//2);m.strategy=2;m.migrad(ncall=2*ncall)
            attempts.append({'nll':float(m.fval),'valid':bool(m.valid),'edm':float(m.fmin.edm)})
            if math.isfinite(m.fval) and m.fval<p.INVALID_NLL/10:candidates.append(m)
        except (ValueError,RuntimeError,FloatingPointError) as exc:
            attempts.append({'error':str(exc)})
    if not candidates:raise RuntimeError('No finite candidate minimum')
    # Retain lowest finite point even if invalid: mark unusable rather than substitute a higher valid fit.
    best=min(candidates,key=lambda m:m.fval)
    best.strategy=2;best.migrad(ncall=2*ncall);best.hesse()
    return s.compact_fit(p.minuit_result_payload(best)),attempts


def refit(toy, states, records, b, nominal, old, truth, dilutions, count, seed, ncall):
    nll,_=p.make_bin_nll(toy,states,records,b,'nominal')
    starts=starts_for(nominal,old,count,seed)
    fit,attempts=search(nll,starts,records,b,ncall=ncall)
    # Keep the original solution as a candidate, never silently accept a worse result.
    if 'fit' in old:
        old_nll=float(nll(**old['fit']['values']))
        if old_nll<fit['minimum_nll']-1e-3:
            fit=dict(old['fit']);fit['minimum_nll']=old_nll
    profiles={};consistent=False;promotions=0
    for cycle in range(3):
        profiles={};lower=[]
        for k in study.TARGETS:
            constrained,trace=search(nll,[fit['values'],*starts],records,b,{k:truth[k]},ncall)
            if constrained['minimum_nll']<fit['minimum_nll']-1e-3:lower.append(constrained['values'])
            profiles[k]={'status':'evaluated' if constrained['valid'] else 'failed',
                         'two_delta_nll':max(0.,2*(constrained['minimum_nll']-fit['minimum_nll'])),
                         'minimum_nll':constrained['minimum_nll'],'attempts':trace}
        if not lower:
            consistent=True;break
        if cycle<2:
            # Release all parameters from any better constrained point, then repeat both profiles.
            fit,trace=search(nll,[fit['values'],*lower],records,b,ncall=ncall)
            attempts.extend(trace);promotions+=len(lower)
    checks=study.evaluate_recovery(toy,states,records,b,fit,truth,dilutions,False)
    unresolved=not consistent
    checks.update(profiles=profiles,usable=checks['usable'] and not unresolved,
                  unresolved_lower_profile=unresolved,profile_promotions=promotions)
    return fit,checks,attempts


def init(cache,states,records,manifest,source,out,settings):
    global _W
    events=p.load_event_cache(Path(cache))
    _W=(events,states,records,manifest,Path(source),Path(out),settings)


def run(task):
    events,states,records,manifest,source,out,settings=_W
    b,mode,i,j,rep=task;cell=(b,mode,i,j)
    original=study.cell_path(source,cell)/f'replica_{rep:05d}.json'
    dest=study.cell_path(out,cell)/original.name;dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.with_suffix('.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return 'owned_elsewhere'
        if dest.exists():return 'reused'
        old=json.loads(original.read_text());state=json.loads((original.parent/'state.json').read_text())
        if old['signature']!=manifest['signature'] or old['task']!=list(cell):raise RuntimeError('Original checkpoint identity mismatch')
        e=study.subset_events(events,np.flatnonzero(events['bin_number']==b))
        config=manifest['config'];truth=state['truth'];d=state['true_dilutions']
        e=study.subset_events(e,study.coverage_indices(e,mode,config,b,rep))
        envelope=s.generator_envelopes(e,states,truth,d)
        toy=s.generate_states(e,states,truth,d,np.random.default_rng(np.random.SeedSequence(old['spin_seed_components'])),envelope)
        aux,centers=s.auxiliary_records(records,b,d,np.random.default_rng(np.random.SeedSequence([config['seed'],b,i,j,rep,777])),'fluctuate')
        if centers!=old.get('auxiliary_centers',centers):raise RuntimeError('Auxiliary reconstruction differs')
        if 'fit' in old:
            check_nll,_=p.make_bin_nll(toy,states,aux,b,'nominal')
            reconstructed=float(check_nll(**old['fit']['values']))
            if abs(reconstructed-old['fit']['minimum_nll'])>1e-6:
                raise RuntimeError('Reconstructed trial likelihood differs from saved trial')
        nominal=manifest['baselines'][str(b)]
        fit,checks,attempts=refit(toy,states,aux,b,nominal,old,truth,d,settings['starts'],[90210,b,rep],settings['ncall'])
        row={**old,'fit':fit,**checks,'status':'fit','audit':{'attempts':attempts,
             'original_sha256':s.digest(original),'original_usable':old.get('usable',False),
             'original_fit':old.get('fit'),'original_profiles':old.get('profiles',{}),
             'nll_improvement':old.get('fit',{}).get('minimum_nll',fit['minimum_nll'])-fit['minimum_nll']}}
        s.atomic_json(dest,row)
        return 'refitted'


def products(source,out,manifest):
    manifest=json.loads((out/'manifest.json').read_text())
    comparisons=[]
    # Build a separate view containing audited trials only, never mix in unaudited originals.
    for cell in study.tasks_for(manifest['config']):
        original=study.cell_path(source,cell)/'state.json'
        if not original.exists():continue
        dest=study.cell_path(out,cell);dest.mkdir(parents=True,exist_ok=True)
        state=json.loads(original.read_text());rows=list(dest.glob('replica_*.json'))
        state['complete']=len(rows)==state['requested'];state['status']='complete' if state['complete'] else 'running'
        if json.loads(original.read_text())['status']=='nonphysical_injection':state.update(status='nonphysical_injection',complete=True)
        s.atomic_json(dest/'state.json',state)
        for path in rows:
            r=json.loads(path.read_text());old=r['audit'].get('original_fit')
            for k in p.PHYSICS_PARAMETERS:
                comparisons.append({'bin_number':cell[0],'mode':cell[1],'u1_index':cell[2],'ll1_index':cell[3],
                    'replica':r['replica'],'parameter':k,'old_usable':r['audit']['original_usable'],'new_usable':r['usable'],
                    'nll_improvement':r['audit']['nll_improvement'],'unresolved_lower_profile':r['unresolved_lower_profile'],
                    'old_estimate':old['values'][k] if old else None,'new_estimate':r['fit']['values'][k],
                    'estimate_change':r['fit']['values'][k]-old['values'][k] if old else None})
    study.write_products(out)
    pd.DataFrame(comparisons).to_csv(out/'tables/minimum_comparison.csv',index=False)


def main():
    parser=study.parser();parser.description=__doc__
    parser.set_defaults(output_dir=study.HERE/'output/cosine_minimum_audit')
    parser.add_argument('--source-dir',type=Path,default=study.HERE/'output/cosine_coverage_validation')
    parser.add_argument('--audit-bins',type=int,nargs='+',default=[6],help='Existing bins to audit; default: bin 6')
    parser.add_argument('--extra-starts',type=int,default=8)
    parser.add_argument('--max-calls',type=int,default=50000)
    args=parser.parse_args();source=args.source_dir.resolve();out=args.output_dir.resolve()
    if source==out or source in out.parents or out in source.parents:raise ValueError('Audit and original directories must be separate')
    manifest=json.loads((source/'manifest.json').read_text())
    if args.summarize_only:products(source,out,manifest);return
    if args.extra_starts<1 or args.max_calls<100:raise ValueError('Positive start count and at least 100 calls required')
    # Use the original preparation checks/configuration, including input/cache and source-code identity.
    original_out=args.output_dir
    for key in ('bins','u1_values','ll1_values','replicas','seed','retained_fraction','gap_half_width','profile_coverage'):
        setattr(args,key,manifest['config'][key])
    args.coverage_comparison=len(manifest['config']['modes'])>1;args.output_dir=source
    _,cache,states,records,_,_,_=study.prepare(args);args.output_dir=original_out
    audit_bins=sorted(set(args.audit_bins))
    if not audit_bins or not set(audit_bins)<=set(manifest['config']['bins']):raise ValueError('Audit bins must exist in source configuration')
    view=json.loads(json.dumps(manifest));view['config']['bins']=audit_bins
    settings={'bins':audit_bins,'starts':args.extra_starts,'ncall':args.max_calls,'audit_sha256':s.digest(Path(__file__)),
              'source_signature':manifest['signature']}
    out.mkdir(parents=True,exist_ok=True)
    with (out/'.audit.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);path=out/'audit_manifest.json'
        if path.exists() and json.loads(path.read_text())!=settings:raise RuntimeError('Audit settings/code changed; use a new output directory')
        s.atomic_json(path,settings);s.atomic_json(out/'manifest.json',view)
    tasks=[]
    for cell in study.tasks_for(view['config']):
        for path in sorted(study.cell_path(source,cell).glob('replica_*.json')):
            tasks.append((*cell,int(path.stem.split('_')[1])))
    print(f'Auditing {len(tasks)} existing trials; {args.extra_starts} random starts plus saved/nominal/zero; originals retained',flush=True)
    failures=[]
    with ProcessPoolExecutor(max_workers=min(args.workers,10),mp_context=mp.get_context('spawn'),initializer=init,
            initargs=(str(cache),states,records,manifest,str(source),str(out),settings)) as pool:
        futures={pool.submit(run,t):t for t in tasks}
        for num,future in enumerate(as_completed(futures),1):
            try:status=future.result()
            except Exception as exc:failures.append({'task':futures[future],'error':repr(exc)});status='FAILED'
            if num%10==0 or status=='FAILED':print(f'{num}/{len(tasks)}: {status}',flush=True)
    s.atomic_json(out/f'node_failures_{os.getpid()}.json',failures)
    if not args.distributed:products(source,out,manifest)
    if failures:raise RuntimeError(f'{len(failures)} trials failed; rerun resumes others')

if __name__=='__main__':main()
