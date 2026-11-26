#!/usr/bin/env python3
"""Build the preliminary LU/UL release from six-bin LU and three-bin free-LL fits."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
os.environ.setdefault('MPLBACKEND','Agg')
from pathlib import Path
import argparse,csv,hashlib,json,multiprocessing,shutil,sys
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
HERE=Path(__file__).resolve().parent
ANALYSIS=HERE.parent
sys.path[:0]=[str(ANALYSIS/'asymmetry_extraction'),str(ANALYSIS/'asymmetry_extraction/three_bin_study'),str(ANALYSIS)]
OBS={'lu1':'LU_sin_phi','ul1':'UL_lab_sin_phi','ul2':'UL_lab_sin_2phi'}
PERIODS=['RGB_spring2019_D2','RGB_fall2019_D2','RGB_spring2020_D2','RGC_summer2022_ND3','RGC_fall2022_ND3','RGC_spring2023_ND3']
XB=[.10,.25,.35,.45,.60]
TP=[.05,.25,.45,.65,.85,1.05,1.25]
COARSE=[.05,.45,.85,1.25]
GROUPS=DF=None


def save(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');t.replace(p)


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def initialize(groups,df):
 global GROUPS,DF
 GROUPS,DF=groups,df


def period_worker(task):
 from constrained_fit import fit
 name,b=task;g=GROUPS[name];mask=g['events']['bin']==b
 result=fit([dict(events={k:v[mask] for k,v in g['events'].items() if k!='bin'},states=g['states'])],DF[b])
 return dict(scope=name,bin=b,**result)


def write_csv(p,rows):
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def build(a):
 source,study,out=a.source,a.study,a.output
 audit=json.loads((source/'asymmetry_extraction/input_audit.json').read_text())
 coarse=json.loads((study/'input_audit.json').read_text())
 if audit['identity']['cutset_id']!=coarse['identity']['cutset_id']:raise ValueError('Fine/coarse cut sets differ')
 for key in ['metadata_sha256','files']:
  if audit['identity'][key]!=coarse['identity'][key]:raise ValueError('Fine/coarse sample metadata mismatch: '+key)
 fits={}
 for obs in OBS:
  n=24 if obs=='lu1' else 12
  for b in range(1,n+1):
   p=source/'asymmetry_extraction/fits'/f'combined_bin{b:02d}.json' if obs=='lu1' else study/'fits'/f'free_bin{b:02d}.json'
   r=json.loads(p.read_text())
   if r['status']!='valid':raise ValueError('Invalid release fit: '+str(p))
   if r['fixed_parameters']!={'UU_cos_phi':0.,'UU_cos_2phi':0.}:raise ValueError('Release must use the unconstrained polarized baseline')
   fits[obs,b]=r
 raw=list(csv.DictReader((source/'kinematic_distributions/bin_kinematic_means.csv').open()))
 # Exact pooled means/variances from per-period bin summaries, for both displayed binnings.
 means=[];lookup={}
 for obs in OBS:
  n=6 if obs=='lu1' else 3
  for xi in range(4):
   for ti in range(n):
    members=[xi*6+ti+1] if n==6 else [xi*6+ti*2+1,xi*6+ti*2+2]
    b=xi*n+ti+1
    for variable in ['x','Q2','W','y','tp','pion_p']:
     rr=[r for r in raw if r['sample'] in PERIODS and int(r['analysis_bin']) in members and r['variable']==variable]
     count=sum(int(r['combinations']) for r in rr)
     if count==0:raise ValueError('Missing bin kinematics')
     mean=sum(int(r['combinations'])*float(r['mean']) for r in rr)/count
     # Input standard deviations use population moments in the histogram accumulator.
     variance=sum(int(r['combinations'])*(float(r['std'])**2+(float(r['mean'])-mean)**2) for r in rr)/count
     item=dict(observable=obs,bin=b,x_index=xi,variable=variable,combinations=count,mean=mean,std=float(np.sqrt(max(variance,0))))
     means.append(item);lookup[obs,b,variable]=item
 rows=[]
 for obs,name in OBS.items():
  n=6 if obs=='lu1' else 3;edges=TP if n==6 else COARSE
  for b in range(1,4*n+1):
   r=fits[obs,b];j=r['names'].index(name);xi=(b-1)//n;ti=(b-1)%n
   if lookup[obs,b,'tp']['combinations']!=r['events']:raise ValueError('Kinematics/event-count mismatch')
   rows.append(dict(channel='pim',observable=obs,x_index=xi,bin=b,x=lookup[obs,b,'tp']['mean'],x_low=edges[ti],x_high=edges[ti+1],value=r['values'][j],stat=r['errors'][j],syst=0.,status='measured'))
 out.mkdir(parents=True,exist_ok=True)
 write_csv(out/'measurements.csv',rows);write_csv(out/'bin_kinematics.csv',means)
 dfrows=json.loads((study/'coarse_dilution.json').read_text())
 periods={};signature=hashlib.sha256(json.dumps(dict(identity=audit['identity'],coarse=coarse['signature'],code=sha(__file__),likelihood=sha(ANALYSIS/'asymmetry_extraction/three_bin_study/constrained_fit.py')),sort_keys=True).encode()).hexdigest()
 if not a.skip_periods:
  from inputs import load
  groups,_,_,_=load(source/'channel_selection/selected_inputs.json',a.run_info,source/'dilution_factor/joint_period_method1.json',ANALYSIS/'channel_selection')
  for g in groups.values():
   fine=g['events']['bin'].astype(int)-1;g['events']['bin']=3*(fine//6)+(fine%6)//2+1
  df={r['bin']:(r['df'],r['stat_std']) for r in dfrows}
  pending=[]
  for name in PERIODS[3:]:
   for b in range(1,13):
    p=out/'period_fits'/f'{name}_bin{b:02d}.json'
    if p.exists() and json.loads(p.read_text()).get('signature')==signature:periods[name,b]=json.loads(p.read_text())
    else:pending.append((name,b))
  with ProcessPoolExecutor(max_workers=a.workers,mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(groups,df)) as pool:
   for job in as_completed([pool.submit(period_worker,t) for t in pending]):
    r=job.result();r['signature']=signature;periods[r['scope'],r['bin']]=r;save(out/'period_fits'/f'{r["scope"]}_bin{r["bin"]:02d}.json',r)
    print(f'UL period checks: {len(periods)}/36, {r["scope"]} bin {r["bin"]}: {r["status"]}',flush=True)
  if any(r['status']!='valid' for r in periods.values()):print('Some individual-period fits invalid; see period_fit_summary.json',flush=True)
 save(out/'period_fit_summary.json',dict(requested=not a.skip_periods,valid=sum(r['status']=='valid' for r in periods.values()),total=len(periods),invalid=[dict(scope=r['scope'],bin=r['bin'],status=r['status']) for r in periods.values() if r['status']!='valid']))
 # Preserve every supporting figure and metadata needed to update the source note.
 from export_review import export_review
 # Copy only review-size support, without producing another archive of the full run.
 for component in ['channel_selection','dilution_factor','kinematic_distributions']:
  for folder,dirs,files in os.walk(source/component):
   dirs[:]=[d for d in dirs if d not in ['selected','calibration_cache','__pycache__']]
   for name in files:
    p=Path(folder)/name
    if p.suffix not in ['.png','.json','.csv','.txt','.html'] or 'replica' in name:continue
    dest=out/'support'/component/p.relative_to(source/component);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
 for file in ['input_audit.json','fit_summary.json','period_likelihood_ratio_checks.json','UU_settings.json']:
  p=source/'asymmetry_extraction'/file
  if p.exists():dest=out/'support/asymmetry_extraction'/file;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
 save(out/'dilution_12_bins.json',dfrows)
 save(out/'constraint_stability.json',json.loads((study/'fit_summary.json').read_text()))
 sources={}
 for obs,b in fits:
  src=source/'asymmetry_extraction/fits'/f'combined_bin{b:02d}.json' if obs=='lu1' else study/'fits'/f'free_bin{b:02d}.json'
  sources[str(src)]=sha(src)
  dest=out/'release_fits'/('fine' if obs=='lu1' else 'three_bin')/src.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
 meta=dict(signature=signature,identity=audit['identity'],datasets=audit['datasets'],sources=sources,observables=OBS,LU_bins=24,UL_bins=12,UU='fixed_zero',polarized_constraints='none; LL and LLcosphi float as nuisance amplitudes',uncertainties='statistical, including first-order dilution. No point-to-point or polarization-scale systematics assigned.',period_checks=not a.skip_periods)
 save(out/'release_manifest.json',meta)
 from plots import make_plots
 make_plots(out,rows,periods,source,a.theory,a.dpi)
 from three_bin_kinematics import build as build_kinematics
 build_kinematics(source,out,a.dpi)
 from note_assets import make_assets
 make_assets(out,rows,means,meta,dfrows)
 export_review(out,a.review_output,a.run_info,ANALYSIS/'imports/datasets.json')
 print('Release assets: '+str(out),flush=True)


def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--source',type=Path,required=True)
 p.add_argument('--study',type=Path)
 p.add_argument('--output',type=Path)
 p.add_argument('--workers',type=int,default=6)
 p.add_argument('--skip-periods',action='store_true',help='Build from existing combined fits only; mark period checks incomplete')
 p.add_argument('--run-info',type=Path,default=ANALYSIS/'imports/clas12_run_info_ND3_source.csv')
 p.add_argument('--theory',type=Path,default=ANALYSIS.parent/'pion_theory_comparisons/inputs/photon_referenced/pim')
 p.add_argument('--review-output',type=Path,default=ANALYSIS/'review_preliminary_release')
 p.add_argument('--dpi',type=int,default=220)
 a=p.parse_args();a.source=a.source.resolve();a.study=(a.study or a.source/'three_bin_LLcos_study').resolve();a.output=(a.output or a.source/'preliminary_release').resolve()
 if a.workers<1:p.error('workers must be positive')
 build(a)


if __name__=='__main__':main()
