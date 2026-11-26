#!/usr/bin/env python3
"""Repeat the identical pi-minus fit with zero, reduced-BM and full-BM UU inputs."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/eppi_uu_study_mpl')
from pathlib import Path
import argparse,json,hashlib,csv,sys
import numpy as np
HERE=Path(__file__).resolve().parent
GROUPS=DF=SETTINGS=None

def save(path,value):
 def clean(v):
  if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
  if isinstance(v,(list,tuple)):return [clean(x) for x in v]
  if isinstance(v,np.generic):return clean(v.item())
  if isinstance(v,float) and not np.isfinite(v):return None
  return v
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n');tmp.replace(path)

def initialize(groups,df,settings):
 global GROUPS,DF,SETTINGS,LIMIT
 GROUPS,DF,SETTINGS=groups,df,settings
 try:
  from threadpoolctl import threadpool_limits
  LIMIT=threadpool_limits(limits=1)
 except ImportError:LIMIT=None

def worker(task):
 from likelihood import fit
 scheme,model,b=task;groups=[]
 for name,g in GROUPS.items():
  fine=g['events']['bin'];bins=fine if scheme=='fine' else 2*((fine-1)//6)+((fine-1)%6)//3+1
  mask=bins==b
  if not mask.any():return dict(scheme=scheme,model=model,bin=b,status='empty',period=name)
  groups.append(dict(events={k:v[mask] for k,v in g['events'].items() if k!='bin'},states=g['states']))
 try:r=fit(groups,DF[scheme][b],fixed_uu=SETTINGS[scheme][model][b]['values'])
 except Exception as exc:r=dict(status='failed',error=str(exc))
 return dict(scheme=scheme,model=model,bin=b,UU_input=SETTINGS[scheme][model][b],**r)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--analysis-root',type=Path,required=True)
 p.add_argument('--selected-inputs',type=Path)
 p.add_argument('--coarse-dilution',type=Path,required=True)
 p.add_argument('--theory',type=Path,required=True,help='Directory containing photon_referenced/pim/step*.csv')
 p.add_argument('--baseline-fine',type=Path)
 p.add_argument('--baseline-coarse',type=Path)
 p.add_argument('--output',type=Path,default=HERE/'output')
 p.add_argument('--workers',type=int,default=3)
 a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=True)
 sys.path.insert(1,str(a.analysis_root/'asymmetry_extraction'))
 from inputs import load
 registry_path=a.selected_inputs or a.analysis_root/'channel_selection/output/selected_inputs.json'
 fine_df=a.analysis_root/'dilution_factor/output/joint_period_method1.json'
 groups,df,audit,identity=load(registry_path,a.analysis_root/'imports/clas12_run_info_ND3_source.csv',fine_df,registry_path.parent.parent)
 print('Validated selected samples: '+str(sum(v['rows'] for v in audit.values()))+' combinations',flush=True)
 coarse_rows=json.loads(a.coarse_dilution.read_text());coarse_df={int(r['bin']):(r['df'],r['stat_std']) for r in coarse_rows}
 if set(coarse_df)!=set(range(1,9)):raise ValueError('Need eight coarse dilution bins')
 # Existing coarse audit proves these dilution values belong to the same cut set.
 coarse_audit=a.coarse_dilution.parent/'input_audit.json'
 if not coarse_audit.exists() or json.loads(coarse_audit.read_text())['identity']['cutset_id']!=identity['cutset_id']:raise ValueError('Coarse dilution cut-set mismatch')
 for key in ['metadata_sha256','dilution_sha256']:
  if json.loads(coarse_audit.read_text())['identity'][key]!=identity[key]:raise ValueError('Coarse input mismatch: '+key)
 import uproot
 registry=json.loads(registry_path.read_text());sums={k:np.zeros(n) for k,n in [('fine',24),('coarse',8)]};counts={k:np.zeros(n,dtype=int) for k,n in [('fine',24),('coarse',8)]}
 for name,g in groups.items():
  path=Path(registry['files'][name]['path']);path=path if path.is_absolute() else registry_path.parent/path
  with uproot.open(path) as f:arr=f[registry['tree']].arrays(['analysis_bin','tprime'],library='np')
  if not np.array_equal(arr['analysis_bin'],g['events']['bin']):raise ValueError('Event ordering mismatch')
  for scheme in ['fine','coarse']:
   b=arr['analysis_bin'].astype(int)-1
   if scheme=='coarse':b=2*(b//6)+(b%6)//3
   counts[scheme]+=np.bincount(b,minlength=len(counts[scheme]));sums[scheme]+=np.bincount(b,weights=-arr['tprime'],minlength=len(counts[scheme]))
 means={k:sums[k]/counts[k] for k in sums};settings={};theory_audit={}
 for scheme,n in [('fine',24),('coarse',8)]:
  settings[scheme]={m:{} for m in ['zero','step1_2','step1_2_3']};nb=6 if scheme=='fine' else 2
  for b in range(1,n+1):
   xi=(b-1)//nb;t=float(means[scheme][b-1]);settings[scheme]['zero'][b]=dict(values=[0.,0.],minus_tprime=t,xB=[.2,.3,.4,.5][xi])
   for model in ['step1_2','step1_2_3']:
    path=a.theory/'photon_referenced/pim'/(model+'.csv')
    with path.open() as f:rows=[r for r in csv.DictReader(f) if np.isclose(float(r['xB']),[.2,.3,.4,.5][xi])]
    x=np.array([float(r['minus_tprime']) for r in rows]);q2={float(r['Q2']) for r in rows}
    if len(rows)!=16 or len(q2)!=1 or not x[0]<=t<=x[-1]:raise ValueError('UU input would require extrapolation')
    uu=[float(np.interp(t,x,[float(r[k]) for r in rows])) for k in ['A_UU_cosphi','A_UU_cos2phi']]
    settings[scheme][model][b]=dict(values=uu,minus_tprime=t,xB=[.2,.3,.4,.5][xi],Q2=q2.pop())
    theory_audit[model]=dict(path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
 for scheme in settings:
  for model,bb in settings[scheme].items():
   path=out/'UU_inputs'/f'{scheme}_{model}.csv';path.parent.mkdir(exist_ok=True)
   with path.open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['bin','UU_cos_phi','UU_cos_2phi','mean_minus_tprime','theory_xB','theory_Q2'])
    for b,v in bb.items():w.writerow([b,*v['values'],v['minus_tprime'],v['xB'],v.get('Q2','')])
 signature_payload=dict(identity=identity,coarse_dilution_sha256=hashlib.sha256(a.coarse_dilution.read_bytes()).hexdigest(),theory=theory_audit,settings=settings,likelihood_sha256=hashlib.sha256((HERE/'likelihood.py').read_bytes()).hexdigest(),runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 signature=hashlib.sha256(json.dumps(signature_payload,sort_keys=True).encode()).hexdigest()
 save(out/'input_audit.json',dict(signature=signature,**signature_payload,datasets=audit,means={k:v.tolist() for k,v in means.items()},method='All six RGB/RGC periods, combined fits. Same exposure, event selection, event depolarization ratios, and fine/coarse dilution factors across UU scenarios. UU is constant per fit bin at the measured combined-sample mean tprime and supplied xB/Q2 setting. No theory-UU uncertainty propagated.'))
 results=[];tasks=[]
 for scheme,n in [('fine',24),('coarse',8)]:
  for model in settings[scheme]:
   for b in range(1,n+1):
    file=out/'fits'/scheme/model/f'combined_bin{b:02d}.json'
    if file.exists() and (r:=json.loads(file.read_text())).get('signature')==signature:results.append(r)
    else:tasks.append((scheme,model,b))
 from concurrent.futures import ProcessPoolExecutor,as_completed
 import multiprocessing
 with ProcessPoolExecutor(max_workers=a.workers,mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(groups,dict(fine=df,coarse=coarse_df),settings)) as pool:
  jobs=[pool.submit(worker,t) for t in tasks]
  for job in as_completed(jobs):
   r=job.result();r['signature']=signature;results.append(r);save(out/'fits'/r['scheme']/r['model']/f"combined_bin{r['bin']:02d}.json",r)
   print(f"{len(results)}/96: {r['scheme']} {r['model']} bin {r['bin']} {r['status']}",flush=True)
 from summarize import summarize
 summarize(out,results,means,a.baseline_fine,a.baseline_coarse)
 print('Study complete: '+str(out),flush=True)
if __name__=='__main__':main()
