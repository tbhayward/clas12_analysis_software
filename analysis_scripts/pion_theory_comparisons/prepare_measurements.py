#!/usr/bin/env python3
"""Normalize current pi+ and pi- results without changing any fit values."""
import argparse,csv,hashlib,json,shutil
from pathlib import Path
import numpy as np
OBS={'u1':'UU_cos_phi','u2':'UU_cos_2phi','lu1':'LU_sin_phi','ul1':'UL_lab_sin_phi','ul2':'UL_lab_sin_2phi','ll0':'LL_lab','ll1':'LL_lab_cos_phi'}
XB=[.10,.25,.35,.45,.60];TP=[.05,.25,.45,.65,.85,1.05,1.25]

def prepare(pip_csv,pim_fine,pim_coarse,output):
 output=Path(output);output.mkdir(parents=True,exist_ok=True);rows=[];sources={}
 def source(path):
  p=Path(path);sources[str(p.resolve())]=hashlib.sha256(p.read_bytes()).hexdigest();return p
 def row(channel,obs,xi,b,x,lo,hi,y,stat,syst,status='measured'):
  if not np.all(np.isfinite([x,y,stat,syst])) or stat<0 or syst<0:raise ValueError('Nonfinite or negative error')
  rows.append(dict(channel=channel,observable=obs,x_index=xi,bin=b,x=x,x_low=lo,x_high=hi,value=y,stat=stat,syst=syst,status=status))
 with source(pip_csv).open() as f:table=list(csv.DictReader(f))
 if sorted(int(r['bin_number']) for r in table)!=list(range(1,25)):raise ValueError('Need all 24 pi+ bins')
 for r in table:
  b=int(r['bin_number']);xi=(b-1)//6;ti=(b-1)%6
  if str(r['nominal_fit_valid']).lower() not in ['true','1']:raise ValueError('Invalid pi+ fit')
  for obs in OBS:row('pip',obs,xi,b,float(r['mean_minus_tprime_gev2']),TP[ti],TP[ti+1],float(r[obs]),float(r[obs+'_stat']),float(r[obs+'_point_to_point_systematic']))
 means=json.loads(source(Path(pim_coarse)/'mean_tprime.json').read_text())
 fine_audit=json.loads(source(Path(pim_fine)/'input_audit.json').read_text())
 coarse_audit=json.loads(source(Path(pim_coarse)/'input_audit.json').read_text())
 if fine_audit['identity']['cutset_id']!=coarse_audit['identity']['cutset_id']:raise ValueError('Mixed pi- cut sets')
 if fine_audit['identity']['metadata_sha256']!=coarse_audit['identity']['metadata_sha256']:raise ValueError('Mixed pi- polarization/charge metadata')
 if set(fine_audit['datasets'])!=set(coarse_audit['datasets']):raise ValueError('Mixed pi- sample coverage')
 for name in fine_audit['datasets']:
  if fine_audit['datasets'][name]['rows']!=coarse_audit['datasets'][name]['rows']:raise ValueError('Fine/coarse source counts differ: '+name)
  if fine_audit['identity']['files'][name]['size']!=coarse_audit['identity']['files'][name]['size']:raise ValueError('Fine/coarse selected files differ: '+name)
 for b in range(1,25):
  r=json.loads(source(Path(pim_fine)/'fits'/f'combined_bin{b:02d}.json').read_text())
  if r['status']!='valid':raise ValueError('Invalid pi- LU fit '+str(b))
  xi=(b-1)//6;ti=(b-1)%6;j=r['names'].index('LU_sin_phi')
  row('pim','lu1',xi,b,means['fine'][b-1],TP[ti],TP[ti+1],r['values'][j],r['errors'][j],0.)
 for b in range(1,9):
  r=json.loads(source(Path(pim_coarse)/'fits'/f'combined_bin{b:02d}.json').read_text())
  if r['status']!='valid':raise ValueError('Invalid pi- coarse fit '+str(b))
  xi=(b-1)//2;half=(b-1)%2
  for obs,full in OBS.items():
   if obs=='lu1':continue
   j=r['names'].index(full);fixed=full in r.get('fixed_parameters',{})
   row('pim',obs,xi,b,means['coarse'][b-1],[.05,.65][half],[.65,1.25][half],r['values'][j],r['errors'][j],0.,'fixed_input' if fixed else 'measured')
 shutil.copy2(pip_csv,output/'pip_source.csv')
 fields=list(rows[0])
 with (output/'measurements.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 provenance=dict(sources=sources,cutset_id_pim=fine_audit['identity']['cutset_id'],pip='24 bins, final result table, statistical and point-to-point uncertainties retained',pim='LU: 24 fine bins. UL/LL: 8 coarse bins. UU: fixed inputs, never plotted as measured points. Statistical errors only.',scale_systematics='not plotted',rows=len(rows))
 (output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 print('Normalized '+str(len(rows))+' measurement/input rows')
 return rows
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pip-csv',type=Path,required=True);p.add_argument('--pim-fine',type=Path,required=True);p.add_argument('--pim-coarse',type=Path,required=True);p.add_argument('--output',type=Path,default=Path(__file__).parent/'measurements');a=p.parse_args();prepare(a.pip_csv,a.pim_fine,a.pim_coarse,a.output)
