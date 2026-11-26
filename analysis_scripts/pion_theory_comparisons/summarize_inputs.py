"""Numerical input audits and useful model templates. No fit or sign changes."""
from pathlib import Path
import csv,json,html
import numpy as np
from theory_io import load_theory,load_measurements,OBS

def write(path,rows):
 if not rows:return
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def summarize(root,measurements,output):
 root,output=Path(root),Path(output);destination=output;output=destination/'data/diagnostics';output.mkdir(parents=True,exist_ok=True);model_inputs=destination/'data/model_UU_inputs';model_inputs.mkdir(parents=True,exist_ok=True);theory,audit=load_theory(root);data=load_measurements(measurements)
 strength=[];changes=[];residuals=[]
 for (s,o),blocks in theory.items():
  for xi,b in enumerate(blocks):
   a=b['arrays'];x=a['minus_tprime']
   for t in b['terms']:
    strength.append(dict(species=s,observable=o,xB=float(a['xB'][0]),Q2=float(a['Q2'][0]),term=t,mean=float(np.mean(a[t])),max_abs=float(np.max(np.abs(a[t]))),integrated_abs=float(np.sum(.5*(np.abs(a[t][:-1])+np.abs(a[t][1:]))*np.diff(x)))))
  legacy=root/'legacy_predictions'/(s+'_long_asymmetries_pred.csv')
  if legacy.exists():
   with legacy.open() as f:old=list(csv.DictReader(line for line in f if not line.lstrip().startswith('#')))
   for xi,b in enumerate(blocks):
    a=b['arrays'];v=[r for r in old if float(r['xB'])==[.2,.3,.4,.5][xi]];v=sorted(v,key=lambda r:float(r['minus_tprime']))
    for x,y in zip(a['minus_tprime'],a['total']):
     if float(v[0]['minus_tprime'])<=x<=float(v[-1]['minus_tprime']):
      z=float(np.interp(x,[float(r['minus_tprime']) for r in v],[float(r[OBS[o]]) for r in v]))
      changes.append(dict(species=s,observable=o,xB=[.2,.3,.4,.5][xi],minus_tprime=float(x),previous_interpolated=z,updated=float(y),difference=float(y-z)))
 # Values for a later fixed-UU fit variation, not a change to the current extraction.
 for species in ['pip','pim']:
  uu=[]
  for r in sorted([r for r in data if r['channel']==species and r['observable']=='lu1'],key=lambda r:r['bin']):
   row=dict(bin=r['bin'])
   for obs,name in [('u1','UU_cos_phi'),('u2','UU_cos_2phi')]:
    a=theory[species,obs][r['x_index']]['arrays'];x=r['x']
    if not a['minus_tprime'][0]<=x<=a['minus_tprime'][-1]:raise ValueError('UU export would extrapolate')
    row[name]=float(np.interp(x,a['minus_tprime'],a['total']))
   uu.append(row)
  write(model_inputs/(species+'_UU_inputs_BM_at_mean_tprime.csv'),uu)
 for r in data:
  if r['observable']!='lu1' or r['status']!='measured':continue
  a=theory[r['channel'],'lu1'][r['x_index']]['arrays'];x=r['x']
  if not a['minus_tprime'][0]<=x<=a['minus_tprime'][-1]:continue
  prediction=float(np.interp(x,a['minus_tprime'],a['total']))
  residuals.append(dict(channel=r['channel'],bin=r['bin'],mean_minus_tprime=x,data=r['value'],stat=r['stat'],point_to_point_syst=r['syst'],BM_at_mean_tprime=prediction,data_minus_BM=r['value']-prediction,stat_only_residual_in_sigma=(r['value']-prediction)/r['stat'],interpretation='Pointwise diagnostic at fixed model xB/Q2. Not a bin-folded goodness-of-fit test.'))
 write(output/'component_strengths.csv',strength);write(output/'updated_vs_previous_theory.csv',changes);write(output/'LU_pointwise_residuals.csv',residuals)
 summary={}
 for s in ['pip','pim','pi0']:
  summary[s]={o:max((abs(r['difference']) for r in changes if r['species']==s and r['observable']==o),default=None) for o in OBS}
 (output/'theory_update_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

def extend_gallery(output):
 import shutil
 from figure_layout import legacy_path,write_gallery
 output=Path(output);diagnostics=output/'data/diagnostics';manifest_path=diagnostics/'manifest.json'
 manifest=json.loads(manifest_path.read_text());legacy=output/'.legacy_work'
 entries=list(manifest['plots']);legacy_paths=[]
 if legacy.exists():
  for p in sorted(legacy.rglob('*.png')):
   dest=legacy_path(output,p.relative_to(legacy));dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(p),str(dest));rel=str(dest.relative_to(output));legacy_paths.append(rel)
   entries.append(dict(path=rel,title=dest.stem,category='Original pi+ reference layouts'))
  for p in legacy.glob('*.json'):shutil.copy2(p,diagnostics/('original_pi_plus_'+p.name))
  remaining=[p for p in legacy.rglob('*') if p.is_file() and p.suffix!='.json']
  if remaining:raise ValueError('Unexpected files in generated legacy staging directory')
  shutil.rmtree(legacy)
 else:
  for p in sorted((output/'figures/pi_plus/reference_original_layouts').rglob('*.png')):
   rel=str(p.relative_to(output));legacy_paths.append(rel);entries.append(dict(path=rel,title=p.stem,category='Original pi+ reference layouts'))
 manifest['legacy_plots']=legacy_paths;manifest['total_deliverable_plots']=manifest['number_of_plots']+len(legacy_paths)
 manifest_path.write_text(json.dumps(manifest,indent=2)+'\n');write_gallery(output,entries)
if __name__=='__main__':
 import sys
 summarize(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]));extend_gallery(Path(sys.argv[3]))
