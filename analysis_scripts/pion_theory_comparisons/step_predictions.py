"""Cumulative photon-axis BM predictions with each stage's own denominator."""
import csv,json,hashlib,shutil
from pathlib import Path
import numpy as np
from matplotlib.lines import Line2D
from theory_io import OBS,load_measurements
STAGES=['step1','step1_2','step1_2_3']
COLORS=['#ff7f0e','#d62728','#1f77b4']
LABELS=[r'BM, $\widetilde H_L$, $\widetilde E_L$ and pion pole',r'BM, longitudinal/pole + chiral-odd',r'BM Total, + chiral-even transverse']

def load_steps(root):
 root=Path(root)/'photon_referenced';result={};audit=[]
 for sp in ['pip','pim','pi0']:
  for stage in STAGES:
   path=root/sp/(stage+'.csv')
   with path.open() as f:rows=list(csv.DictReader(f))
   expected={'xB','Q2','minus_tprime','sigL_over_sigT'}|{v+s for v in OBS.values() for s in ['', '_lo','_hi']}
   if set(rows[0])!=expected:raise ValueError('Unexpected columns '+str(path))
   a={k:np.array([float(r[k]) for r in rows]) for k in rows[0]}
   if len(rows)!=64 or not all(np.all(np.isfinite(v)) for v in a.values()):raise ValueError('Invalid table '+str(path))
   if any(np.any(v<0) for k,v in a.items() if k.endswith(('_lo','_hi'))):raise ValueError('Negative error widths')
   if not np.array_equal(np.unique(a['xB']),[.2,.3,.4,.5]):raise ValueError('Wrong xB grid')
   for xi,xb in enumerate([.2,.3,.4,.5]):
    b={k:v[a['xB']==xb] for k,v in a.items()}
    if np.any(np.diff(b['minus_tprime'])<=0) or len(np.unique(b['Q2']))!=1:raise ValueError('Unsorted grid')
    result[sp,stage,xi]=b
   audit.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),rows=len(rows)))
 for sp in ['pip','pim','pi0']:
  for xi in range(4):
   for stage in STAGES:
    for key in ['xB','Q2','minus_tprime']:
     if not np.array_equal(result[sp,stage,xi][key],result[sp,'step1',xi][key]):raise ValueError('Stage grid mismatch')
 return result,audit

def handles(gk=False):
 return ([Line2D([],[],color='#b0b0b0',lw=2.3,label='GK Total')] if gk else [])+[Line2D([],[],color=c,lw=2.4,label=l) for c,l in zip(COLORS,LABELS)]

def curves(ax,steps,sp,obs,xi,gk=None):
 key=OBS[obs]
 for stage,c in zip(STAGES,COLORS):
  b=steps[sp,stage,xi];x=b['minus_tprime'];y=b[key]
  ax.fill_between(x,y-b[key+'_lo'],y+b[key+'_hi'],color=c,alpha=.10,lw=0,zorder=2)
  ax.plot(x,y,color=c,lw=2.15 if stage!='step1_2_3' else 2.6,zorder=4)
 available=gk is not None and sp=='pip' and obs in ['lu1','ul1','ul2','ll0','ll1']
 if available:
  v=gk[np.isclose(gk[:,0],[.2,.3,.4,.5][xi])];v=v[np.argsort(v[:,1])]
  ax.plot(v[:,1],v[:,5+['lu1','ul1','ul2','ll0','ll1'].index(obs)],color='#b0b0b0',lw=2.3,zorder=3)
 return available

def focused(steps,data,sp,obs,output,dpi,alpha,gk=None,variant='bm'):
 import matplotlib.pyplot as plt
 import plot_comparisons as p
 from figure_layout import OBS_NAMES
 fig,axes=plt.subplots(2,2,figsize=(12.4,9.2));available=False
 for xi,ax in enumerate(axes.flat):
  p.axis(ax,obs,p.limits(None,None,sp,obs));lo,hi=p.XB[xi];b=steps[sp,'step1',xi]
  ax.set_title(rf'${lo:.2f}<x_B<{hi:.2f}$'+'  ·  '+rf'$Q^2={b["Q2"][0]:.1f}\;\mathrm{{GeV}}^2$')
  available=curves(ax,steps,sp,obs,xi,gk)
  p.points(ax,p.state_data(data,sp,obs,xi),'combined_errors' if 'combined_errors' in variant else 'stat_syst' if sp=='pip' else 'stat')
  if sp=='pim' and obs in ['u1','u2']:ax.text(.03,.92,'UU fixed in fit; no measured points',transform=ax.transAxes,fontsize=9)
 fig.legend(handles=handles(available)+p.error_handles('combined_errors' if 'combined_errors' in variant else 'stat_syst' if sp=='pip' else 'stat',sp),loc='lower center',bbox_to_anchor=(.5,.035),ncol=2,frameon=False,fontsize=10.2,handlelength=3)
 fig.suptitle(p.CHANNEL[sp]+'  ·  '+p.LABELS[obs],fontsize=17,y=.98)
 fig.subplots_adjust(left=.09,right=.98,bottom=.20,top=.90,hspace=.32,wspace=.24);p.watermark(fig,alpha)
 channel={'pip':'pi_plus','pim':'pi_minus'}[sp]
 folder=Path(output)/'figures'/channel/'step_predictions'/variant;folder.mkdir(parents=True,exist_ok=True)
 path=folder/f'{channel}_{OBS_NAMES[obs]}_cumulative_predictions_2x2.png';fig.savefig(path,dpi=dpi);plt.close(fig)
 return dict(path=str(path.relative_to(output)),title=channel+' '+obs+' cumulative predictions',category='cumulative BM '+variant)

def seven_panel(steps,data,sp,xi,output,dpi,alpha,gk=None,combined_errors=False):
 import matplotlib.pyplot as plt
 import plot_comparisons as p
 fig,axes=plt.subplots(3,3,figsize=(12.8,10.5));flat=axes.flat
 for panel,obs in zip(p.PANELS,OBS):
  ax=flat[panel];p.axis(ax,obs,p.limits(None,None,sp,obs));curves(ax,steps,sp,obs,xi,gk);p.points(ax,p.state_data(data,sp,obs,xi),'combined_errors' if combined_errors else 'stat_syst' if sp=='pip' else 'stat')
 flat[2].axis('off');flat[8].axis('off');lo,hi=p.XB[xi]
 flat[2].text(0,1,rf'${lo:.2f}<x_B<{hi:.2f}$'+'\n'+rf'$Q^2={steps[sp,"step1",xi]["Q2"][0]:.1f}\;\mathrm{{GeV}}^2$',va='top',fontsize=14)
 flat[8].legend(handles=handles(gk is not None and sp=='pip')+p.error_handles('combined_errors' if combined_errors else 'stat_syst' if sp=='pip' else 'stat',sp),loc='upper left',frameon=False,fontsize=8.2,handlelength=2)
 fig.subplots_adjust(left=.08,right=.98,bottom=.07,top=.92,hspace=.38,wspace=.32);fig.suptitle(p.CHANNEL[sp]+' cumulative BM predictions',fontsize=18);p.watermark(fig,alpha)
 channel={'pip':'pi_plus','pim':'pi_minus'}[sp];folder=Path(output)/'figures'/channel/'step_predictions'/('bm+gk_combined_errors' if combined_errors else 'bm+gk' if gk is not None else 'bm');folder.mkdir(parents=True,exist_ok=True)
 path=folder/f'{channel}_seven_ratios_cumulative_xB_{xi+1:02}.png';fig.savefig(path,dpi=dpi);plt.close(fig)
 return dict(path=str(path.relative_to(output)),title=channel+' cumulative seven ratios',category='cumulative BM')

def run(theory_path,measurements,output,dpi=300,alpha=.10,channels=('pip','pim')):
 import plot_comparisons as p
 from shared_plots import make_pdf
 from figure_layout import OBS_NAMES
 p.style();output=Path(output);steps,audit=load_steps(theory_path);data=load_measurements(measurements)
 gk_path=Path(theory_path)/'gk_pip.csv';gk=np.loadtxt(gk_path,delimiter=',',skiprows=1) if gk_path.exists() else None
 entries=[]
 for sp in channels:
  for obs in OBS:
   entries.append(focused(steps,data,sp,obs,output,dpi,alpha))
   if sp=='pip':entries.append(focused(steps,data,sp,obs,output,dpi,alpha,gk,'bm+gk'))
   entries.append(focused(steps,data,sp,obs,output,dpi,alpha,gk,'bm+gk_combined_errors'))
  for xi in range(4):
   entries.append(seven_panel(steps,data,sp,xi,output,dpi,alpha))
   if sp=='pip':entries.append(seven_panel(steps,data,sp,xi,output,dpi,alpha,gk))
   entries.append(seven_panel(steps,data,sp,xi,output,dpi,alpha,gk,combined_errors=True))
 from types import SimpleNamespace
 from theory_io import load_theory
 full,_=load_theory(theory_path)
 args=SimpleNamespace(output=output,watermark_alpha=alpha,dpi=dpi,ylimits={(sp,'lu1',False):(-.2,.3) for sp in ['pip','pim']})
 p.lu_charge_comparison(full,data,args)
 p.lu_charge_comparison(full,data,args,include_pi0=True)
 p.lu_charge_comparison(full,data,args,residual=True)
 p.lu_charge_comparison(full,data,args,include_pi0=True,combined_errors=True)
 entries.append(dict(path='figures/charge_comparison/pi_plus_pi_minus_LU_sin_phi_with_pi0_theory_combined_errors_2x2.png',title='Charge comparison, combined errors',category='charge comparison'))
 selected=[('pip',o) for o in OBS]+[('pim','lu1')]
 for variant in ['bm','bm+gk','bm+gk_combined_errors']:
  folder=output/'figures/shared_plots'/variant;folder.mkdir(parents=True,exist_ok=True);ordered=[]
  for sp,obs in selected:
   channel={'pip':'pi_plus','pim':'pi_minus'}[sp];filename=f'{channel}_{OBS_NAMES[obs]}_cumulative_predictions_2x2.png'
   source=output/'figures'/channel/'step_predictions'/('bm' if sp=='pim' and variant!='bm+gk_combined_errors' else variant)/filename
   if source.exists():
    destination=folder/filename;shutil.copy2(source,destination);ordered.append(destination);entries.append(dict(path=str(destination.relative_to(output)),title=channel+' '+obs+' cumulative predictions',category='shared '+variant))
  comparison=output/'figures/charge_comparison'/('pi_plus_pi_minus_LU_sin_phi_with_pi0_theory_combined_errors_2x2.png' if variant=='bm+gk_combined_errors' else 'pi_plus_pi_minus_LU_sin_phi_with_pi0_theory_2x2.png')
  if comparison.exists():
   destination=folder/comparison.name;shutil.copy2(comparison,destination);ordered.append(destination);entries.append(dict(path=str(destination.relative_to(output)),title='π+ / π− LU sinφ with π0 theory',category='shared '+variant))
  if len(ordered)==9:make_pdf(ordered,folder/'pion_chiral_sector_comparisons_9_pages.pdf')
  for stale in folder.glob('*_chiral_sector_sums_2x2.png'):stale.unlink()
 diag=output/'data/diagnostics';diag.mkdir(parents=True,exist_ok=True)
 (diag/'cumulative_predictions.json').write_text(json.dumps(dict(inputs=audit,stages=STAGES,interpretation='Complete ratios recalculated at each cumulative stage, not additive contributions. Each uses its own denominator.',frame='virtual photon',uncertainties='Supplied lo/hi treated as lower/upper error widths, consistent with previous prediction-table convention. Bands are theory fit uncertainties, not experimental systematics.',plots=entries),indent=2)+'\n')
 mf=diag/'manifest.json'
 if mf.exists():
  m=json.loads(mf.read_text());lookup={e['path']:e for e in m['plots'] if (output/e['path']).exists()};lookup.update({e['path']:e for e in entries});m['plots']=list(lookup.values());m['number_of_plots']=len(lookup);m['cumulative_theory_inputs']=audit;mf.write_text(json.dumps(m,indent=2)+'\n')
 print(f'Created {len(entries)} cumulative prediction plot entries and shared PDFs',flush=True)
 return entries
