"""Three-observable preliminary figures. No LL panels or GK overlays."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OBS=['lu1','ul1','ul2']
LABELS={'lu1':r'$F_{LU}^{\sin\phi}/F_{UU}$','ul1':r'$F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$','ul2':r'$F_{UL,\mathrm{lab}}^{\sin2\phi}/F_{UU}$'}
NAMES={'lu1':'LU_sin_phi','ul1':'UL_sin_phi','ul2':'UL_sin_2phi'}
PARAM={'lu1':'LU_sin_phi','ul1':'UL_lab_sin_phi','ul2':'UL_lab_sin_2phi'}
KEYS={'lu1':'A_LU_sinphi','ul1':'A_UL_sinphi','ul2':'A_UL_sin2phi'}
COLORS=['#0072B2','#D55E00','#009E73','#8A4B9C']
MARKERS=['o','s','^','D']
XB=[.10,.25,.35,.45,.60]
STAGES=['step1','step1_2','step1_2_3']
TC=['#ff7f0e','#d62728','#1f77b4']
TL=['BM longitudinal/pole','BM + chiral-odd','BM full']


def setup(ax,obs):
 ax.set(xlim=(0,1.3),ylim={'lu1':(-.2,.3),'ul1':(-.4,.4),'ul2':(-1,1)}[obs],xlabel=r"$-t^\prime$ (GeV$^2$)",ylabel=LABELS[obs]);ax.axhline(0,color='.6',lw=.8);ax.grid(alpha=.12)


def mark(fig):
 fig.text(.5,.51,'PRELIMINARY',ha='center',va='center',rotation=25,fontsize=fig.get_figwidth()*7.2,color='red',alpha=.10,zorder=0)


def points(ax,rr,color='black',marker='o'):
 for r in rr:
  lo,hi=ax.get_ylim();y=r['value'];visible=np.clip(y,lo+.02*(hi-lo),hi-.02*(hi-lo))
  ax.errorbar(r['x'],visible,yerr=r['stat'] if lo<=y<=hi else None,fmt=marker if lo<=y<=hi else ('^' if y>hi else 'v'),color=color,mfc='white',capsize=3,ms=5,zorder=5)
  if y+r['stat']>hi:ax.plot(r['x'],hi-.015*(hi-lo),'^',color=color,ms=5)
  if y-r['stat']<lo:ax.plot(r['x'],lo+.015*(hi-lo),'v',color=color,ms=5)


def make_plots(out,rows,periods,source,theory,dpi):
 plt.rcParams.update({'font.family':'serif','mathtext.fontset':'stix','font.size':11})
 figures=out/'figures';figures.mkdir(exist_ok=True)
 data=lambda obs,xi:[r for r in rows if r['observable']==obs and r['x_index']==xi]
 tables={}
 for stage in STAGES:
  p=theory/(stage+'.csv')
  with p.open() as f:tables[stage]=list(csv.DictReader(f))
 def curves(ax,obs,xi,stages):
  for stage in stages:
   rr=sorted([r for r in tables[stage] if np.isclose(float(r['xB']),[.2,.3,.4,.5][xi])],key=lambda r:float(r['minus_tprime']))
   x=np.array([float(r['minus_tprime']) for r in rr]);k=KEYS[obs];y=np.array([float(r[k]) for r in rr]);i=STAGES.index(stage)
   ax.plot(x,y,color=TC[i],lw=2,label=TL[i]);ax.fill_between(x,y-[float(r[k+'_lo']) for r in rr],y+[float(r[k+'_hi']) for r in rr],color=TC[i],alpha=.10)
 def finish(fig,path,legend=None):
  if legend:fig.legend(handles=legend,loc='lower center',ncol=min(4,len(legend)),frameon=False,fontsize=9)
  fig.tight_layout(rect=[0,.09 if legend else 0,1,.95]);mark(fig);fig.savefig(figures/path,dpi=dpi);plt.close(fig)
 for theory_mode in ['data','BM_total','BM_stages']:
  handles=([Line2D([],[],color=COLORS[i],marker=MARKERS[i],ls='',mfc='white',label=f'{XB[i]:.2f} < xB < {XB[i+1]:.2f}') for i in range(4)] if theory_mode=='data' else [Line2D([],[],color=TC[i],lw=2,label=TL[i]) for i in ([2] if theory_mode=='BM_total' else range(3))])
  fig,axes=plt.subplots(4,3,figsize=(13,13))
  for xi in range(4):
   for j,obs in enumerate(OBS):
    ax=axes[xi,j];setup(ax,obs);points(ax,data(obs,xi))
    if theory_mode!='data':curves(ax,obs,xi,['step1_2_3'] if theory_mode=='BM_total' else STAGES)
    ax.set_title(f'{XB[xi]:.2f} < xB < {XB[xi+1]:.2f}')
  fig.suptitle('π− preliminary spin-dependent ratios'+(' · '+theory_mode.replace('_',' ') if theory_mode!='data' else ''))
  finish(fig,f'{theory_mode}_all_xB.png',handles if theory_mode!='data' else None)
  for xi in range(4):
   fig,axes=plt.subplots(1,3,figsize=(13,4.3))
   for ax,obs in zip(axes,OBS):
    setup(ax,obs);points(ax,data(obs,xi))
    if theory_mode!='data':curves(ax,obs,xi,['step1_2_3'] if theory_mode=='BM_total' else STAGES)
   fig.suptitle(f'π− · {XB[xi]:.2f} < xB < {XB[xi+1]:.2f}')
   finish(fig,f'{theory_mode}_xB_bin_{xi+1}.png',handles if theory_mode!='data' else None)
  for obs in OBS:
   fig,axes=plt.subplots(2,2,figsize=(10.5,8))
   for xi,ax in enumerate(axes.flat):
    setup(ax,obs);points(ax,data(obs,xi))
    if theory_mode!='data':curves(ax,obs,xi,['step1_2_3'] if theory_mode=='BM_total' else STAGES)
    ax.set_title(f'{XB[xi]:.2f} < xB < {XB[xi+1]:.2f}')
   fig.suptitle('π− · '+LABELS[obs]);finish(fig,f'{theory_mode}_{NAMES[obs]}_2x2.png',handles if theory_mode!='data' else None)
 fig,axes=plt.subplots(1,3,figsize=(13,4.3))
 for ax,obs in zip(axes,OBS):
  setup(ax,obs)
  for xi in range(4):points(ax,data(obs,xi),COLORS[xi],MARKERS[xi])
 finish(fig,'data_overlay_all_xB.png',[Line2D([],[],color=COLORS[i],marker=MARKERS[i],ls='',label=f'{XB[i]:.2f} < xB < {XB[i+1]:.2f}') for i in range(4)])
 # Existing fine fits provide LU period results. New target-period fits use merged bins.
 period_rows=[]
 for xi in range(4):
  fig,axes=plt.subplots(1,3,figsize=(13,4.5))
  handles=[]
  for j,obs in enumerate(OBS):
   ax=axes[j];setup(ax,obs)
   names=['RGB_spring2019_D2','RGB_fall2019_D2','RGB_spring2020_D2','RGC_summer2022_ND3','RGC_fall2022_ND3','RGC_spring2023_ND3'] if obs=='lu1' else ['RGC_summer2022_ND3','RGC_fall2022_ND3','RGC_spring2023_ND3']
   for ni,name in enumerate(names):
    rr=[]
    for point in data(obs,xi):
     b=point['bin']
     if obs=='lu1':
      p=source/'asymmetry_extraction/fits'/f'{name}_bin{b:02d}.json'
      r=json.loads(p.read_text()) if p.exists() else None
     else:r=periods.get((name,b))
     if not r or r['status']!='valid':continue
     k=r['names'].index(PARAM[obs]);row=dict(point,value=r['values'][k],stat=r['errors'][k],x=point['x']+(ni-(len(names)-1)/2)*.008)
     rr.append(row);period_rows.append(dict(scope=name,observable=obs,bin=b,value=row['value'],stat=row['stat']))
    global_index=['RGB_spring2019_D2','RGB_fall2019_D2','RGB_spring2020_D2','RGC_summer2022_ND3','RGC_fall2022_ND3','RGC_spring2023_ND3'].index(name)
    points(ax,rr,plt.get_cmap('tab10')(global_index),'os^DvP'[global_index])
    if j==0:handles.append(Line2D([],[],color=plt.get_cmap('tab10')(ni),marker='os^DvP'[ni],ls='',label=name))
   if obs!='lu1' and not periods:ax.text(.5,.9,'Period fits pending',ha='center',transform=ax.transAxes)
  fig.suptitle(f'Period checks · {XB[xi]:.2f} < xB < {XB[xi+1]:.2f}');finish(fig,f'periods_xB_bin_{xi+1}.png',handles)
 if period_rows:
  with (out/'period_results.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(period_rows[0]));w.writeheader();w.writerows(period_rows)
 df=json.loads((out/'dilution_12_bins.json').read_text())
 fig,axes=plt.subplots(2,2,figsize=(10,7))
 for xi,ax in enumerate(axes.flat):
  rr=[r for r in df if (r['bin']-1)//3==xi];xx=[data('ul1',xi)[(r['bin']-1)%3]['x'] for r in rr]
  ax.errorbar(xx,[r['df'] for r in rr],yerr=[r['stat_std'] for r in rr],fmt='o',capsize=3)
  ax.set(xlim=(0,1.3),ylim=(0,1),xlabel=r"$-t^\prime$ (GeV$^2$)",ylabel='ND3 dilution factor',title=f'{XB[xi]:.2f} < xB < {XB[xi+1]:.2f}');ax.grid(alpha=.12)
 finish(fig,'dilution_12_bins.png')
 links=''.join(f'<p><a href="figures/{p.name}">{p.stem}</a></p>' for p in sorted(figures.glob('*.png')))
 (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><h1>π− preliminary release</h1><p>Fine LU, three-bin UL. Statistical errors. LL amplitudes float but are not displayed.</p>'+links)
