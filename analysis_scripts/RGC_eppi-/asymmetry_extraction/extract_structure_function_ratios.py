#!/usr/bin/env python3
"""Conditional seven-ratio extraction, combined RGB+RGC and individual periods."""
import argparse,csv,hashlib,json,os
from pathlib import Path
for variable in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[variable]='1'
import numpy as np
from inputs import load,PERIODS
from likelihood import fit,NAMES
HERE=Path(__file__).resolve().parent
GROUPS=DF=UU=None

def safe(payload):
    if isinstance(payload,dict):return {k:safe(v) for k,v in payload.items()}
    if isinstance(payload,list):return [safe(v) for v in payload]
    if isinstance(payload,float) and not np.isfinite(payload):return None
    return payload

def save(path,payload):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(safe(payload),indent=2,allow_nan=False)+'\n');temporary.replace(path)

def initialize(groups,df,uu):
    global GROUPS,DF,UU,THREAD_LIMIT
    GROUPS,DF,UU=groups,df,uu
    try:
        from threadpoolctl import threadpool_limits
        THREAD_LIMIT=threadpool_limits(limits=1)
    except ImportError:
        THREAD_LIMIT=None

def worker(task):
    b,scope,signature=task
    names=PERIODS if scope=='combined' else [scope]
    groups=[]
    for name in names:
        g=GROUPS[name];mask=g['events']['bin']==b
        if not mask.any():return dict(bin=b,scope=scope,status='empty',signature=signature)
        groups.append(dict(events={k:v[mask] for k,v in g['events'].items() if k!='bin'},states=g['states']))
    try:result=fit(groups,DF[b],fixed_uu=UU[b])
    except Exception as e:result=dict(status='failed',error=str(e))
    if UU[b] is not None:result['fixed_parameters']={NAMES[0]:UU[b][0],NAMES[1]:UU[b][1]}
    return dict(bin=b,scope=scope,signature=signature,**result)

def products(output,results,registry_path=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import chi2
    lookup={(r['scope'],r['bin']):r for r in results}
    folder=output/'figures';folder.mkdir(exist_ok=True)
    labels=['RGB Sp19','RGB Fa19','RGB Sp20','RGC Su22','RGC Fa22','RGC Sp23']
    colors=['#0072B2','#56B4E9','#009E73','#D55E00','#CC79A7','#E69F00']
    markers=['o','s','^','D','v','P']
    titles=[r'$F_{UU}^{\cos\phi}/F_{UU}$',r'$F_{UU}^{\cos2\phi}/F_{UU}$',r'$F_{LU}^{\sin\phi}/F_{UU}$',r'$F_{UL,lab}^{\sin\phi}/F_{UU}$',r'$F_{UL,lab}^{\sin2\phi}/F_{UU}$',r'$F_{LL,lab}/F_{UU}$',r'$F_{LL,lab}^{\cos\phi}/F_{UU}$']
    xedges=[.10,.25,.35,.45,.60];x=np.arange(.15,1.25,.2)
    rows=[]
    for result in results:
        if result['status']!='valid':continue
        for j,name in enumerate(result['names']):rows.append(dict(scope=result['scope'],bin=result['bin'],parameter=name,parameter_status='fixed_input' if name in result.get('fixed_parameters',{}) else 'fitted',value=result['values'][j],stat=result['errors'][j],events=result['events'],warnings='|'.join(result['warnings'])))
    if rows:
        with (output/'ratios.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    from plot_production_results import plot
    plot(output,registry_path)
    tests=[]
    for b in range(1,25):
        nominal=lookup.get(('combined',b),{})
        individual=[lookup.get((name,b),{}) for name in PERIODS]
        if nominal.get('status')=='valid' and all(r.get('status')=='valid' for r in individual):
            statistic=2*(nominal['nll']-sum(r['nll'] for r in individual))
            ndf=sum(len(r.get('free_parameters',r['names'])) for r in individual)-len(nominal.get('free_parameters',nominal['names']))
            tests.append(dict(bin=b,likelihood_ratio=statistic,ndf=ndf,p_value=float(chi2.sf(max(0,statistic),ndf)),warning='Asymptotic diagnostic; weak conditional UU identification can invalidate the chi-square approximation'))
    save(output/'period_likelihood_ratio_checks.json',tests)
    summary={scope:{status:sum(r['scope']==scope and r['status']==status for r in results) for status in ['valid','invalid','failed','empty']} for scope in ['combined']+PERIODS}
    save(output/'fit_summary.json',summary)
    from matplotlib.colors import ListedColormap
    scopes=['combined']+PERIODS
    image=np.full((7,24),3)
    for i,scope in enumerate(scopes):
        for b in range(1,25):
            r=lookup.get((scope,b),{})
            image[i,b-1]={'valid':0,'invalid':1,'failed':2}.get(r.get('status'),3)
    fig,ax=plt.subplots(figsize=(12,4))
    ax.imshow(image,aspect='auto',cmap=ListedColormap(['#009E73','#E69F00','#D55E00','#cccccc']),vmin=-.5,vmax=3.5)
    ax.set(xticks=np.arange(24),xticklabels=np.arange(1,25),yticks=np.arange(7),yticklabels=['Combined']+labels,xlabel='Analysis bin',title='Fit status: green = valid covariance, orange = invalid, red = failed, gray = missing')
    fig.tight_layout();fig.savefig(folder/'fit_status.png',dpi=180);plt.close(fig)

    print(json.dumps(summary,indent=2),flush=True)


def read_uu_inputs(path):
    with Path(path).open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    result={}
    for row in rows:
        b=int(row['bin'])
        if b in result or b not in range(1,25):raise ValueError('Duplicate or invalid UU bin: '+str(b))
        values=(float(row['UU_cos_phi']),float(row['UU_cos_2phi']))
        if not np.all(np.isfinite(values)):raise ValueError('Nonfinite imported UU ratio')
        result[b]=values
    if set(result)!=set(range(1,25)):raise ValueError('UU import must cover all 24 bins exactly once')
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selected-inputs',type=Path,default=HERE.parent/'channel_selection/output/selected_inputs.json')
    p.add_argument('--run-info',type=Path,default=HERE.parent/'imports/clas12_run_info_ND3_source.csv')
    p.add_argument('--dilution',type=Path,default=HERE.parent/'dilution_factor/output/joint_period_method1.json')
    p.add_argument('--output',type=Path)
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--uu-inputs',type=Path,help='CSV with bin,UU_cos_phi,UU_cos_2phi for all 24 bins')
    mode.add_argument('--float-uu',action='store_true',help='Reproduce the seven-free-ratio fit')
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--bins',type=int,nargs='+',default=list(range(1,25)))
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--status',action='store_true')
    a=p.parse_args()
    if a.output is None:
        folder='output_free_UU' if a.float_uu else ('output_fixed_UU_'+a.uu_inputs.stem if a.uu_inputs else 'output_fixed_UU')
        a.output=HERE/folder
    if a.status:
        summary=a.output/'fit_summary.json'
        print(summary.read_text() if summary.exists() else f'{len(list((a.output/"fits").glob("*.json")))} fit checkpoints present');return
    if a.workers<1 or any(b not in range(1,25) for b in a.bins):p.error('Need positive workers and bins 1-24')
    uu={b:None if a.float_uu else (0.,0.) for b in range(1,25)}
    if a.uu_inputs:
        uu=read_uu_inputs(a.uu_inputs)
    groups,df,audit,identity=load(a.selected_inputs,a.run_info,a.dilution,HERE.parent/'channel_selection')
    a.output.mkdir(parents=True,exist_ok=True)
    identity['UU_settings']={str(b):v for b,v in uu.items()}
    identity['UU_source']='free' if a.float_uu else (str(a.uu_inputs.resolve()) if a.uu_inputs else 'fixed_zero')
    if a.uu_inputs:identity['UU_sha256']=hashlib.sha256(a.uu_inputs.read_bytes()).hexdigest()
    save(a.output/'UU_settings.json',dict(values=identity['UU_settings'],source=identity['UU_source'],fixed=not a.float_uu))
    identity['code_sha256']=hashlib.sha256(b''.join((HERE/n).read_bytes() for n in ['likelihood.py','inputs.py','extract_structure_function_ratios.py'])).hexdigest()
    signature=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    save(a.output/'input_audit.json',dict(identity=identity,signature=signature,datasets=audit))
    for name,r in audit.items():print(f"{name}: {r['rows']:,} combinations, {r['source_runs']} source runs; {r['charge_assumption']}",flush=True)
    if a.preflight:return
    directory=a.output/'fits';directory.mkdir(exist_ok=True)
    tasks=[];results=[]
    for b in a.bins:
        for scope in ['combined']+PERIODS:
            path=directory/f'{scope}_bin{b:02d}.json'
            if path.exists():
                r=json.loads(path.read_text())
                if r.get('signature')==signature:results.append(r);continue
            tasks.append((b,scope,signature))
    from concurrent.futures import ProcessPoolExecutor,as_completed
    with ProcessPoolExecutor(max_workers=a.workers,initializer=initialize,initargs=(groups,df,uu)) as pool:
        jobs={pool.submit(worker,t):t for t in tasks}
        for job in as_completed(jobs):
            r=job.result();results.append(r)
            save(directory/f"{r['scope']}_bin{r['bin']:02d}.json",r)
            print(f"{len(results)}/{7*len(a.bins)}: {r['scope']} bin {r['bin']} {r['status']}",flush=True)
    products(a.output,results,a.selected_inputs)

if __name__=='__main__':main()
