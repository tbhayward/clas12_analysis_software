#!/usr/bin/env python3
"""Three-tprime-bin polarized refits, with free and BM-fixed LLcosphi scenarios."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[key] = '1'
os.environ.setdefault('MPLBACKEND', 'Agg')
import argparse
import csv
import hashlib
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parents[1]
sys.path.insert(0, str(ANALYSIS))
sys.path.insert(0, str(ANALYSIS/'asymmetry_extraction'))
MODELS = ['free', 'step1_2', 'step1_2_3']
MEMBERS = [[6*ix+2*j+k+1 for k in range(2)] for ix in range(4) for j in range(3)]
PERIODS = ['RGB_spring2019_D2','RGB_fall2019_D2','RGB_spring2020_D2','RGC_summer2022_ND3','RGC_fall2022_ND3','RGC_spring2023_ND3']
GROUPS = DF = SETTINGS = None


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def initialize(groups, df, settings):
    global GROUPS, DF, SETTINGS
    GROUPS, DF, SETTINGS = groups, df, settings


def worker(task):
    from constrained_fit import fit
    model, b = task
    groups = []
    for g in GROUPS.values():
        mask = g['events']['bin'] == b
        if not mask.any():
            raise ValueError(f'Empty period in bin {b}')
        groups.append(dict(events={k:v[mask] for k,v in g['events'].items() if k != 'bin'}, states=g['states']))
    result = fit(groups, DF[b], SETTINGS[model][str(b)]['fixed_LLcos'])
    return dict(model=model, bin=b, fine_bins=MEMBERS[b-1], theory_input=SETTINGS[model][str(b)], **result)


def prepare(a):
    source = json.loads((a.source/'dilution_factor/method1.json').read_text())
    fine_df = json.loads((a.source/'dilution_factor/joint_period_method1.json').read_text())
    fine_audit = json.loads((a.source/'asymmetry_extraction/input_audit.json').read_text())
    cutset = source['cutset_id']
    if fine_df['cutset_id'] != cutset or fine_audit['identity']['cutset_id'] != cutset:
        raise ValueError('Input cut sets disagree')
    coarse = dict(source, inputs=[dict(item, selected_counts=[sum(item['selected_counts'][b-1] for b in members) for members in MEMBERS]) for item in source['inputs']])
    # Verify the RGC normalization records against current metadata. RGB additions need not change dilution.
    sys.path.insert(0, str(ANALYSIS/'asymmetry_extraction'))
    from inputs import read_rows
    rows = read_rows(a.run_info)
    for item in source['inputs']:
        for run, q in item['run_charges'].items():
            v = rows[int(run)]
            if not np.isclose(v[1]+v[2], q, rtol=1e-10, atol=1e-8):
                raise ValueError(f'Changed RGC charge for {run}')
    sums = np.zeros(12)
    counts = np.zeros(12)
    with (a.source/'kinematic_distributions/bin_kinematic_means.csv').open() as f:
        for row in csv.DictReader(f):
            if row['sample'] not in PERIODS or row['variable'] != 'tp':
                continue
            fine = int(row['analysis_bin'])-1
            b = 3*(fine//6)+(fine%6)//2
            n = int(row['combinations'])
            sums[b] += n*float(row['mean'])
            counts[b] += n
    if np.any(counts == 0):
        raise ValueError('Missing measured bin means')
    expected = np.zeros(12, dtype=int)
    for name in PERIODS:
        for b, members in enumerate(MEMBERS):
            expected[b] += sum(fine_audit['datasets'][name]['bin_counts'][i-1] for i in members)
    if not np.array_equal(counts, expected):
        raise ValueError('Kinematic mean counts do not match extraction sample')
    means = sums/counts
    settings = {m:{} for m in MODELS}
    theory_hashes = {}
    for model in MODELS:
        table = None
        if model != 'free':
            path = a.theory/(model+'.csv')
            theory_hashes[model] = digest(path)
            with path.open() as f:
                table = list(csv.DictReader(f))
        for b in range(1,13):
            xb = [.2,.3,.4,.5][(b-1)//3]
            value, q2 = None, None
            if table is not None:
                rows = sorted([r for r in table if np.isclose(float(r['xB']),xb)], key=lambda r:float(r['minus_tprime']))
                x = np.array([float(r['minus_tprime']) for r in rows])
                if len(rows) != 16 or not x[0] <= means[b-1] <= x[-1]:
                    raise ValueError('Theory grid incomplete or extrapolation required')
                value = float(np.interp(means[b-1],x,[float(r['A_LL_cosphi']) for r in rows]))
                q2 = float(rows[0]['Q2'])
            settings[model][str(b)] = dict(fixed_LLcos=value, mean_minus_tprime=float(means[b-1]), theory_xB=xb, theory_Q2=q2)
    provenance = dict(cutset_id=cutset, members=MEMBERS, minus_tprime_edges=[.05,.45,.85,1.25],
                      counts=expected.tolist(), means=means.tolist(), theory_sha256=theory_hashes,
                      source_counts_sha256=digest(a.source/'dilution_factor/method1.json'),
                      metadata_sha256=digest(a.run_info), df_replicas=a.df_replicas,
                      method='Three tprime bins per xB. UU fixed zero. Same events, charges and dilution across scenarios. BM LLcos ratios fixed at combined-sample mean tprime using supplied xB/Q2 grids. Photon-axis predictions applied directly to lab-axis fit at user request. No theory uncertainty propagated. LU six-bin nominal result remains separate.')
    for name, value in [('aggregated_target_counts.json',coarse),('LLcos_inputs.json',settings),('input_audit.json',provenance)]:
        save(a.output/name,value)
    return coarse, settings, provenance


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True,help='Completed relaxed analysis root (or review root for prepare-only)')
    p.add_argument('--output',type=Path)
    p.add_argument('--run-info',type=Path,default=ANALYSIS/'imports/clas12_run_info_ND3_source.csv')
    p.add_argument('--theory',type=Path,default=HERE/'theory')
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--df-replicas',type=int,default=2000)
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--review-output',type=Path,default=ANALYSIS/'review_three_bin_LLcos')
    a = p.parse_args()
    if a.workers < 1 or a.df_replicas < 2:
        p.error('Positive worker count and at least two dilution replicas required')
    a.source = a.source.resolve()
    a.output = (a.output or a.source/'three_bin_LLcos_study').resolve()
    a.output.mkdir(parents=True,exist_ok=True)
    coarse, settings, provenance = prepare(a)
    if a.prepare_only:
        print('Prepared 12 bins and all three LLcos scenarios: '+str(a.output))
        return
    from inputs import load
    registry = a.source/'channel_selection/selected_inputs.json'
    groups, _, audit, identity = load(registry,a.run_info,a.source/'dilution_factor/joint_period_method1.json',ANALYSIS/'channel_selection')
    if identity['cutset_id'] != provenance['cutset_id']:
        raise ValueError('Selected data mismatch')
    import uproot
    from selection import load_selected_inputs
    selected = load_selected_inputs(registry)
    # Check all 16 material count vectors against the selected trees, not just saved totals.
    for item in json.loads((a.source/'dilution_factor/method1.json').read_text())['inputs']:
        with uproot.open(selected['files'][item['dataset']]['path']) as f:
            counts = np.bincount(f[selected['tree']]['analysis_bin'].array(library='np'),minlength=25)[1:25]
        if not np.array_equal(counts,item['selected_counts']):
            raise ValueError('Control count mismatch: '+item['dataset'])
    for name,g in groups.items():
        fine = g['events']['bin'].astype(int)-1
        g['events']['bin'] = 3*(fine//6)+(fine%6)//2+1
        expected = [sum(audit[name]['bin_counts'][i-1] for i in members) for members in MEMBERS]
        if not np.array_equal(np.bincount(g['events']['bin'],minlength=13)[1:],expected):
            raise ValueError('Event rebinning mismatch')
    sys.path.insert(0,str(ANALYSIS/'dilution_factor'))
    from joint_period_fit import replica_bin
    df_signature = hashlib.sha256(json.dumps(dict(counts=coarse['inputs'],replicas=a.df_replicas,seed=20261124,code=digest(ANALYSIS/'dilution_factor/joint_period_fit.py'),model=digest(ANALYSIS/'dilution_factor/model.py')),sort_keys=True).encode()).hexdigest()
    ctx = multiprocessing.get_context('spawn')
    dfrows = []
    with ProcessPoolExecutor(max_workers=a.workers,mp_context=ctx) as pool:
        pending = {}
        for b in range(12):
            path = a.output/f'dilution_bin{b+1:02d}.json'
            if path.exists() and json.loads(path.read_text()).get('signature') == df_signature:
                dfrows.append(json.loads(path.read_text()))
            else:
                pending[pool.submit(replica_bin,(coarse,b,a.df_replicas,20261124))] = b
        for job in as_completed(pending):
            row, draws = job.result()
            row['signature'] = df_signature
            dfrows.append(row)
            save(a.output/f'dilution_bin{row["bin"]:02d}.json',row)
            np.save(a.output/f'dilution_replicas_bin{row["bin"]:02d}.npy',draws)
            print(f'Dilution {len(dfrows)}/12: bin {row["bin"]}',flush=True)
    save(a.output/'coarse_dilution.json',sorted(dfrows,key=lambda r:r['bin']))
    df = {r['bin']:(r['df'],r['stat_std']) for r in dfrows}
    signature = hashlib.sha256(json.dumps(dict(identity=identity,df=df,settings=settings,runner=digest(__file__),likelihood=digest(HERE/'constrained_fit.py'),base=digest(ANALYSIS/'asymmetry_extraction/uu_model_study/likelihood.py')),sort_keys=True).encode()).hexdigest()
    provenance.update(identity=identity,datasets=audit,signature=signature)
    save(a.output/'input_audit.json',provenance)
    (a.output/'fits').mkdir(exist_ok=True)
    results, tasks = [], []
    for model in MODELS:
        for b in range(1,13):
            path = a.output/'fits'/f'{model}_bin{b:02d}.json'
            if path.exists() and json.loads(path.read_text()).get('signature') == signature:
                results.append(json.loads(path.read_text()))
            else:
                tasks.append((model,b))
    with ProcessPoolExecutor(max_workers=a.workers,mp_context=ctx,initializer=initialize,initargs=(groups,df,settings)) as pool:
        for job in as_completed([pool.submit(worker,t) for t in tasks]):
            result = job.result()
            result['signature'] = signature
            results.append(result)
            save(a.output/'fits'/f'{result["model"]}_bin{result["bin"]:02d}.json',result)
            print(f'Fit {len(results)}/36: {result["model"]} bin {result["bin"]}, {result["status"]}',flush=True)
    from summarize import products
    products(a.output,results,provenance,dfrows,a.source/'asymmetry_extraction')
    sys.path.insert(0,str(ANALYSIS))
    from export_review import export_review
    export_review(a.output,a.review_output,a.run_info,ANALYSIS/'imports/datasets.json')
    if any(r['status'] != 'valid' for r in results):
        raise RuntimeError('Some fits invalid. Inspect summary and fit JSON files.')
    print('Study complete: '+str(a.output),flush=True)


if __name__ == '__main__':
    main()
