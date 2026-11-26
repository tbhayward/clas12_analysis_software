"""Strict readers for supplied CFF products and their common-denominator sums."""
import csv,json,hashlib
from pathlib import Path
import numpy as np
OBS={'u1':'A_UU_cosphi','u2':'A_UU_cos2phi','lu1':'A_LU_sinphi','ul1':'A_UL_sinphi','ul2':'A_UL_sin2phi','ll0':'A_LL','ll1':'A_LL_cosphi'}
META={'xB','Q2','minus_tprime','eps','total'}
TOKENS={'F11':r'F_{1,1}','G11':r'G_{1,1}','G14':r'G_{1,4}','F14':r'F_{1,4}','F12':r'F_{1,2}','G12':r'G_{1,2}','HT':r'H_T','EbarT':r'\bar E_T','HtL':r'\widetilde H_L','EtL_np':r'\widetilde E_L^{\rm np}','EtL_pole':r'\widetilde E_L^{\rm pole}','calE':r'\mathcal{E}','calEt':r'\widetilde{\mathcal{E}}'}

def label(term):
 if term=='other':return 'Other supplied terms'
 if term=='remainder':return 'Remaining supplied terms'
 pieces=term.split('x')
 if len(pieces)!=2 or any(t not in TOKENS for t in pieces):raise ValueError('Unknown CFF notation: '+term)
 return '$'+r'\times '.join(TOKENS[t] for t in pieces)+'$'

def load_theory(root):
 root=Path(root);theory={};audit=[]
 for species in ['pip','pim','pi0']:
  for obs,stem in OBS.items():
   path=root/species/(stem+'.csv')
   with path.open() as f:
    reader=csv.DictReader(f);names=reader.fieldnames;rows=list(reader)
   if len(names)!=len(set(names)) or not META.issubset(names):raise ValueError('Bad header: '+str(path))
   terms=[k for k in names if k not in META]
   for t in terms:label(t)
   table={k:np.array([float(r[k]) for r in rows]) for k in names}
   if not all(np.all(np.isfinite(a)) for a in table.values()):raise ValueError('Nonfinite theory '+str(path))
   if not np.array_equal(np.unique(table['xB']),[.2,.3,.4,.5]):raise ValueError('Wrong xB groups')
   data=[]
   for xb in [.2,.3,.4,.5]:
    mask=table['xB']==xb;a={k:v[mask] for k,v in table.items()}
    if np.any(np.diff(a['minus_tprime'])<=0) or len(np.unique(a['Q2']))!=1 or np.any(a['Q2']<=0):raise ValueError('Bad theory grid')
    delta=sum(a[t] for t in terms)-a['total'];bound=5e-6*(len(terms)+1)+1e-9
    if np.max(np.abs(delta))>bound:raise ValueError('CFF terms do not close to total '+str(path))
    data.append(dict(arrays=a,terms=terms))
   theory[species,obs]=data
   audit.append(dict(species=species,observable=obs,path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),rows=len(rows),terms=terms,max_abs_closure=float(np.max(np.abs(sum(table[t] for t in terms)-table['total']))),theory_uncertainty='none supplied'))
 if (root/'photon_referenced').exists():
  from step_predictions import load_steps
  stages,stage_audit=load_steps(root)
  for (species,obs),blocks in theory.items():
   for xi,block in enumerate(blocks):
    new=stages[species,'step1_2_3',xi];a=block['arrays']
    if not np.array_equal(a['minus_tprime'],new['minus_tprime']) or not np.array_equal(a['Q2'],new['Q2']):raise ValueError('New full-model grid differs from decomposition')
    if not np.array_equal(a['total'],new[OBS[obs]]):raise ValueError('Updated full model differs from old CFF decomposition. Rebuild product tables before plotting together.')
    a['total']=new[OBS[obs]].copy()
  audit.extend(dict(item,source='cumulative photon-referenced tables; full stage used for all BM totals') for item in stage_audit)
 return theory,audit

def load_measurements(path):
 with Path(path).open() as f:rows=list(csv.DictReader(f))
 for r in rows:
  for k in ['x','x_low','x_high','value','stat','syst']:r[k]=float(r[k])
  for k in ['x_index','bin']:r[k]=int(r[k])
  if r['status'] not in ['measured','fixed_input']:raise ValueError('Unknown measurement status')
  if not np.all(np.isfinite([r[k] for k in ['x','x_low','x_high','value','stat','syst']])):raise ValueError('Nonfinite data')
  if r['stat']<0 or r['syst']<0:raise ValueError('Negative uncertainty')
  if not r['x_low']<=r['x']<=r['x_high']:raise ValueError('Mean point outside its bin')
 return rows
