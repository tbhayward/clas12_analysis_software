"""Shared DF with independent period target rates and a shared Sp23 transfer."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import xlogy
from scipy.optimize import brentq
from model import COEFFICIENTS as C
import csv
P=['summer2022','fall2022','spring2023']

def fit_bin(inputs,b,fixed=None,start=None):
 keys=[(p,t) for p in P for t in ['ND3','C','CH2','He','ET']]+[('spring2023','CD2')]
 table={(i['dataset'].split('_')[1],i['dataset'].split('_')[2]):i for i in inputs}
 n=np.array([table[k]['selected_counts'][b] for k in keys],float)
 q=np.array([table[k]['charge'] for k in keys])
 scale=max(np.max((n+.5)/q),1e-10)**2
 def constraint(z):
  r={k:v for k,v in zip(keys,np.exp(z[:16])/q)}
  ratio=r[('spring2023','CD2')]/r[('spring2023','CH2')]
  out=[]
  for p in P:
   a,c,h,e=[r[(p,t)] for t in ['ND3','C','He','ET']]
   cd=r[('spring2023','CD2')] if p=='spring2023' else ratio*r[(p,'CH2')]
   s=cd+C['A']*c+C['B']*h+C['C']*e
   d=cd+C['D']*c+C['E']*h+C['F']*e
   out.append(((a-h)*s-z[16]*a*d)/scale)
  return np.array(out)
 def objective(z):
  mu=np.exp(z[:16]);return float(np.sum(mu-n+xlogy(n,n)-xlogy(n,mu)))
 bounds=[(-25, np.log(max(n.max()*100,1000)))]*16+([(fixed,fixed)] if fixed is not None else [(1e-5,.99999)])
 guesses=[start] if start is not None else [np.r_[np.log(n+.5),v] for v in [.25,.4,.6]]
 solutions=[]
 for z in guesses:
  z=z.copy()
  if fixed is not None:z[16]=fixed
  result=minimize(objective,z,method='SLSQP',bounds=bounds,constraints=[dict(type='eq',fun=constraint)],options=dict(maxiter=600,ftol=1e-9))
  if result.success and np.max(np.abs(constraint(result.x)))<1e-6:solutions.append(result)
 if not solutions:raise RuntimeError(f'Bin {b+1}: joint fit failed')
 return min(solutions,key=lambda x:x.fun)

def replica_bin(job):
 data,b,replicas,seed=job
 best=fit_bin(data['inputs'],b)
 rng=np.random.default_rng(np.random.SeedSequence([seed,b]))
 inputs=[dict(i,selected_counts=list(i['selected_counts'])) for i in data['inputs']]
 draws=[];failed=0
 for trial in range(replicas):
  for original,item in zip(data['inputs'],inputs):
   item['selected_counts'][b]=int(rng.poisson(original['selected_counts'][b]))
  try:
   try:result=fit_bin(inputs,b,start=best.x)
   except RuntimeError:result=fit_bin(inputs,b)
   if result.x[16]<2e-5 or result.x[16]>1-2e-5:
    alternative=fit_bin(inputs,b)
    if alternative.fun<result.fun:result=alternative
   draws.append(float(result.x[16]))
  except RuntimeError:failed+=1
 if failed/replicas>.01:raise RuntimeError(f'Bin {b+1}: more than 1% replica fits failed')
 values=np.array(draws)
 row=dict(bin=b+1,df=float(best.x[16]),stat_std=float(values.std(ddof=1)),
          replica_q16=float(np.quantile(values,.16)),replica_q84=float(np.quantile(values,.84)),
          replica_mean=float(values.mean()),replicas_requested=replicas,replicas_successful=len(values),
          replicas_failed=failed,boundary_fraction=float(np.mean((values<2e-5)|(values>1-2e-5))),
          deviance=2*float(best.fun),ndf=2)
 return row,values


def run(source,output,replicas=2000,seed=20261124,workers=6):
 from concurrent.futures import ProcessPoolExecutor,as_completed
 data=json.loads(Path(source).read_text());rows=[];samples={}
 with ProcessPoolExecutor(max_workers=workers) as pool:
  jobs=[pool.submit(replica_bin,(data,b,replicas,seed)) for b in range(24)]
  for job in as_completed(jobs):
   row,values=job.result();rows.append(row);samples['bin_'+str(row['bin'])]=values
   print(f"Bin {row['bin']}: {row['df']:.4f} +/- {row['stat_std']:.4f}, {row['replicas_failed']} failed trials, boundary {row['boundary_fraction']:.3%}",flush=True)
 rows.sort(key=lambda r:r['bin'])
 Path(output).mkdir(exist_ok=True,parents=True)
 Path(output,'joint_period_method1.json').write_text(json.dumps(dict(status='joint_period_fit',cutset_id=data['cutset_id'],
  method='Joint Poisson fit with shared DF and period-specific target rates. Symmetric statistical error is the sample standard deviation of fits to independent Poisson count replicas, matching the pi+ reporting convention.',
  replicas=replicas,seed=seed,bins=rows),indent=2)+'\n')
 with Path(output,'joint_period_method1.csv').open('w',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 np.savez_compressed(Path(output)/'joint_period_replicas.npz',**samples)
 plot_results(data,rows,Path(output))
 return rows

def plot_results(data,rows,output):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.size':11,'axes.titlesize':12,'axes.labelsize':12,'legend.fontsize':10})
 xedges=[.10,.25,.35,.45,.60];tedges=np.array([.05,.25,.45,.65,.85,1.05,1.25])
 fig,axes=plt.subplots(2,2,figsize=(11,7.5),sharex=True,sharey=True)
 for ix,ax in enumerate(axes.flat):
  values=rows[6*ix:6*ix+6];x=(tedges[:-1]+tedges[1:])/2
  y=np.array([r['df'] for r in values]);error=np.array([r['stat_std'] for r in values])
  ax.errorbar(x,y,yerr=error,fmt='o',color='#0072B2',capsize=2,linestyle='none',label='Combined RGC')
  ax.set(title=f'{xedges[ix]:.2f} ≤ xB < {xedges[ix+1]:.2f}',xlabel='−t′ (GeV²)',ylabel='Dilution factor',xlim=(.05,1.25),ylim=(0,.65))
 axes.flat[0].legend(frameon=False);fig.tight_layout();fig.savefig(output/'figures/joint_period_method1.png',dpi=200);plt.close(fig)
 # Matched display for the two period-agreement projections.
 import json
 comparison=json.loads((output/'method1_binning_comparison.json').read_text())
 fig,axes=plt.subplots(1,2,figsize=(12,4.7),sharey=True)
 for ax,scheme,xlabel,title in zip(axes,['xb_only','tprime_only'],['xB','−t′ (GeV²)'],['Integrated over −t′','Integrated over xB']):
  for ip,(period,color,marker) in enumerate(zip(['Su22','Fa22','Sp23'],['#0072B2','#D55E00','#009E73'],['o','s','^'])):
   v=[r for r in comparison['bins'] if r['scheme']==scheme and r['period']==period]
   prefix='xB' if scheme=='xb_only' else 'minus_tprime'
   offset=.004 if scheme=='xb_only' else .014
   x=np.array([(r[prefix+'_low']+r[prefix+'_high'])/2 for r in v])+(ip-1)*offset
   y=np.array([r['dilution_value'] for r in v]);error=np.array([r['dilution_stat_std'] for r in v])
   ax.errorbar(x,y,yerr=error,fmt=marker,linestyle='none',color=color,capsize=2,label=period)
  ax.set(xlabel=xlabel,ylabel='Dilution factor',title=title,ylim=(.15,.60))
 axes[0].legend(frameon=False,ncol=3);fig.tight_layout();fig.savefig(output/'figures/period_agreement.png',dpi=200);plt.close(fig)

 # Period-integrated packing cross-check from the same input target data.
 packing=data['packing_fraction_integrated']
 fig,ax=plt.subplots(figsize=(6.5,4))
 for i,(r,color,marker) in enumerate(zip(packing,['#0072B2','#D55E00','#009E73'],['o','s','^'])):
  ax.errorbar(i,r['value'],yerr=r['stat_std'],fmt=marker,color=color,capsize=4,markersize=7)
 ax.set(xticks=range(3),xticklabels=[r['period'] for r in packing],ylabel='Cell packing fraction',ylim=(.45,.68),xlim=(-.5,2.5))
 fig.tight_layout();fig.savefig(output/'figures/packing_fraction_consistency.png',dpi=200);plt.close(fig)

if __name__=='__main__':
 import sys
 run(sys.argv[1],sys.argv[2])
