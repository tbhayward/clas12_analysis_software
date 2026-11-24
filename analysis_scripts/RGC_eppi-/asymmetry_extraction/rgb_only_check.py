#!/usr/bin/env python3
"""Isolate RGB LU/UU fits and compare strict and legacy support checks."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import argparse,csv,json
from pathlib import Path
import numpy as np
from iminuit import Minuit
from inputs import load,PERIODS
from likelihood import ConditionalLikelihood,NAMES
HERE=Path(__file__).resolve().parent
G=None

class RGBObjective:
 def __init__(self,groups,legacy=False):
  self.core=ConditionalLikelihood(groups,target=False);self.legacy=legacy
  self.q=np.concatenate([g['events']['q_fraction'] for g in groups])
 def factors(self,x):
  c=self.core;n=1+c.observed@x;d=1+c.average@x
  good=np.all(n>1e-10) and np.all(d>1e-15)
  if self.legacy:good=good and np.all(self.q*n/d<=1+1e-10)
  else:good=good and np.min(1+c.allowed@x)>1e-10
  return n,d,good
 def __call__(self,*x):
  n,d,good=self.factors(np.array(x))
  return float(np.log(d).sum()-np.log(n).sum()) if good else 1e30
 def gradient(self,*x):
  n,d,good=self.factors(np.array(x));c=self.core
  return np.sum(c.average/d[:,None]-c.observed/n[:,None],axis=0) if good else np.zeros(3)


def one(task):
 b,scope,variant=task
 names=PERIODS[:3] if scope=='RGB_combined' else [scope]
 selected=[]
 for name in names:
  g=G[name];mask=g['events']['bin']==b
  selected.append(dict(events={k:v[mask] for k,v in g['events'].items() if k!='bin'},states=g['states']))
 cost=RGBObjective(selected,legacy=variant=='pi_plus_checks_and_bounds')
 fits=[]
 for seed in [[0,0,.05],[.3,-.15,.08],[-.3,.15,-.05]]:
  m=Minuit(cost,*seed,name=NAMES[:3],grad=cost.gradient);m.errordef=.5
  if variant!='strict_unbounded':m.limits=[(-1.5,1.5)]*3
  m.migrad(ncall=4000);m.hesse()
  if not(m.valid and m.accurate):m.simplex(ncall=2000);m.strategy=2;m.migrad(ncall=6000);m.hesse()
  fits.append(m)
 valid=[m for m in fits if m.valid and m.accurate]
 best=min(valid or fits,key=lambda m:m.fval);x=np.array(best.values)
 c=cost.core;errors=np.array(best.errors)
 return dict(bin=b,scope=scope,variant=variant,status='valid' if valid else 'invalid',events=c.events,
  values=x.tolist(),errors=errors.tolist() if valid else [None]*3,nll=float(best.fval),edm=float(best.fmin.edm),
  minimum_allowed_factor=float(np.min(1+c.allowed@x)),minimum_observed_factor=float(np.min(1+c.observed@x)),
  parameter_limit_proximity=bool(variant!='strict_unbounded' and np.any(1.5-np.abs(x)<2*errors)),
  correlation=np.asarray(best.covariance.correlation()).tolist() if valid and best.covariance is not None else None,
  covariance=np.asarray(best.covariance).tolist() if valid and best.covariance is not None else None)


def initialize(groups):
 global G
 G=groups


def run(registry,metadata,dilution,output,workers):
 from inputs import read_rows
 groups,df,audit,identity=load(registry,metadata,dilution,HERE.parent/'channel_selection')
 rows=read_rows(metadata)
 for name in PERIODS[:3]:
  a=groups[name]['events'];exposure=audit[name]['charge']
  a['q_fraction']=np.array([rows[int(r)][0]/2/exposure for r in a['runnum']])
 output.mkdir(parents=True,exist_ok=True)
 from concurrent.futures import ProcessPoolExecutor,as_completed
 tasks=[(b,scope,variant) for b in range(1,25) for scope in ['RGB_combined']+PERIODS[:3] for variant in (['strict_unbounded','strict_pi_plus_bounds','pi_plus_checks_and_bounds'] if scope=='RGB_combined' else ['strict_unbounded'])]
 results=[]
 with ProcessPoolExecutor(max_workers=workers,initializer=initialize,initargs=(groups,)) as pool:
  for future in as_completed([pool.submit(one,t) for t in tasks]):results.append(future.result())
 results.sort(key=lambda r:(r['scope'],r['variant'],r['bin']))
 def clean(v):
  if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
  if isinstance(v,list):return [clean(x) for x in v]
  if isinstance(v,float) and not np.isfinite(v):return None
  return v
 (output/'rgb_only_fits.json').write_text(json.dumps(clean(dict(input_identity=identity,diagnostic_only=True,results=results)),indent=2,allow_nan=False)+'\n')
 flat=[]
 for r in results:
  row={k:v for k,v in r.items() if k not in ['values','errors','correlation','covariance']}
  for j,name in enumerate(NAMES[:3]):row[name]=r['values'][j];row[name+'_stat']=r['errors'][j]
  flat.append(row)
 with (output/'rgb_only_fits.csv').open('w',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(flat[0]),lineterminator='\n');writer.writeheader();writer.writerows(flat)
 summary=[]
 for scope in ['RGB_combined']+PERIODS[:3]:
  for variant in ['strict_unbounded','strict_pi_plus_bounds','pi_plus_checks_and_bounds']:
   r=[x for x in results if x['scope']==scope and x['variant']==variant]
   if not r:continue
   summary.append(dict(scope=scope,variant=variant,valid=sum(x['status']=='valid' for x in r),invalid_bins=[x['bin'] for x in r if x['status']!='valid'],negative_unobserved_factor_bins=[x['bin'] for x in r if x['minimum_allowed_factor']<0],near_parameter_bounds=[x['bin'] for x in r if x['parameter_limit_proximity']]))
 (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
 plot(output,results)
 information_check(groups,output)
 return results


def information_check(groups,output):
 rows=[]
 for b in range(1,25):
  a={k:np.concatenate([groups[name]['events'][k][groups[name]['events']['bin']==b] for name in PERIODS[:3]]) for k in ['phi','rb','rv','rw','pb']}
  for lu in [.03,.20]:
   c1=a['rv']*np.cos(a['phi']);c2=a['rb']*np.cos(2*a['phi']);beam=a['pb']*a['rw']*np.sin(a['phi'])
   u=1+.15*c1-.15*c2;s=lu*beam;p=(1+s/u)/2
   derivative=.5*np.column_stack([-s*c1/u**2,-s*c2/u**2,beam/u])
   weighted=derivative/np.sqrt(p*(1-p))[:,None]
   covariance=np.linalg.inv(weighted.T@weighted)
   rows.append(dict(bin=b,events=len(u),illustrative_truth=True,UU_cos_phi=.15,UU_cos_2phi=-.15,LU_sin_phi=lu,
     expected_UU_cos_phi_stat=float(np.sqrt(covariance[0,0])),expected_UU_cos_2phi_stat=float(np.sqrt(covariance[1,1])),expected_LU_sin_phi_stat=float(np.sqrt(covariance[2,2]))))
 with (output/'expected_information.csv').open('w',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)


def plot(output,results):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from plot_production_results import axes_style,XB,TP
 folder=output/'figures';folder.mkdir(exist_ok=True)
 variants=['strict_unbounded','strict_pi_plus_bounds','pi_plus_checks_and_bounds']
 for ix in range(4):
  fig,axes=plt.subplots(1,3,figsize=(14,4.6))
  for j,ax in enumerate(axes):
   axes_style(ax,j)
   for k,(variant,color,marker,label) in enumerate(zip(variants,['black','#0072B2','#D55E00'],['o','s','^'],['Strict support, unbounded','Strict support, ±1.5 limits','Pi+ support checks, ±1.5 limits'])):
    r=sorted([x for x in results if x['scope']=='RGB_combined' and x['variant']==variant and ix*6<x['bin']<=(ix+1)*6 and x['status']=='valid'],key=lambda x:x['bin'])
    x=[(TP[(v['bin']-1)%6]+TP[(v['bin']-1)%6+1])/2+(k-1)*.012 for v in r]
    ax.errorbar(x,[v['values'][j] for v in r],yerr=[v['errors'][j] for v in r],fmt=marker,linestyle='none',color=color,ms=4,capsize=2,label=label)
  fig.suptitle(f'RGB only: {XB[ix]:.2f} ≤ xB < {XB[ix+1]:.2f}')
  handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.94),ncol=3,frameon=False,fontsize=9)
  fig.tight_layout(rect=(0,0,1,.85));fig.savefig(folder/f'RGB_only_implementation_comparison_xB_bin_{ix+1}.png',dpi=200);plt.close(fig)

 for ix in range(4):
  fig,axes=plt.subplots(1,3,figsize=(14,4.6))
  selected=sorted([r for r in results if r['scope']=='RGB_combined' and r['variant']=='strict_unbounded' and ix*6<r['bin']<=(ix+1)*6],key=lambda r:r['bin'])
  missing=[r['bin'] for r in selected if r['status']!='valid']
  for j,ax in enumerate(axes):
   axes_style(ax,j)
   valid=[r for r in selected if r['status']=='valid']
   x=[(TP[(r['bin']-1)%6]+TP[(r['bin']-1)%6+1])/2 for r in valid]
   ax.errorbar(x,[r['values'][j] for r in valid],yerr=[r['errors'][j] for r in valid],fmt='o',linestyle='none',color='black',ms=4,capsize=3)
   limit=.4 if j==2 else 1.
   for xx,r in zip(x,valid):
    if r['values'][j]+r['errors'][j]>limit:ax.plot(xx,.97*limit,'^',color='black',ms=4)
    if r['values'][j]-r['errors'][j]<-limit:ax.plot(xx,-.97*limit,'v',color='black',ms=4)
  fig.suptitle(f'Combined RGB: {XB[ix]:.2f} ≤ xB < {XB[ix+1]:.2f}')
  note='Statistical errors only. Points at bin centers. Edge markers indicate intervals outside the fixed axes.'
  if missing:note+=' Invalid fits omitted: '+', '.join(map(str,missing))+'.'
  fig.text(.5,.015,note,ha='center',fontsize=9)
  fig.tight_layout(rect=(0,.05,1,.93));fig.savefig(folder/f'RGB_only_xB_bin_{ix+1}.png',dpi=200);plt.close(fig)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--selected-inputs',type=Path,default=HERE.parent/'channel_selection/output/selected_inputs.json')
 p.add_argument('--run-info',type=Path,default=HERE.parent/'imports/clas12_run_info_ND3_source.csv')
 p.add_argument('--dilution',type=Path,default=HERE.parent/'dilution_factor/output/joint_period_method1.json')
 p.add_argument('--output',type=Path,default=HERE/'output/rgb_only')
 p.add_argument('--workers',type=int,default=6)
 a=p.parse_args();run(a.selected_inputs,a.run_info,a.dilution,a.output,a.workers)
