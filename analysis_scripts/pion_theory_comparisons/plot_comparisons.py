#!/usr/bin/env python3
"""Rebuild publication-style pion comparisons from measured ratios and CFF CSVs."""
import os
from pathlib import Path
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'pion_theory_mpl'))
import argparse,csv,hashlib,json,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle
from matplotlib.ticker import AutoMinorLocator,MaxNLocator
from theory_io import OBS,load_theory,load_measurements,label
from cff_style import configure_styles
from figure_layout import figure_path,write_gallery
HERE=Path(__file__).resolve().parent
ORDER=list(OBS);PANELS=[0,1,3,4,5,6,7];XB=[(.1,.25),(.25,.35),(.35,.45),(.45,.60)]
LABELS={'u1':r'$F_{UU}^{\cos\phi}/F_{UU}$','u2':r'$F_{UU}^{\cos2\phi}/F_{UU}$','lu1':r'$F_{LU}^{\sin\phi}/F_{UU}$','ul1':r'$F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$','ul2':r'$F_{UL,\mathrm{lab}}^{\sin2\phi}/F_{UU}$','ll0':r'$F_{LL,\mathrm{lab}}/F_{UU}$','ll1':r'$F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$'}
CHANNEL={'pip':r'$\pi^+$','pim':r'$\pi^-$','pi0':r'$\pi^0$'}
XCOLORS=['#0072B2','#D55E00','#009E73','#8A4B9C'];MARKERS=['o','s','^','D']
IMPORTANT=['HTxEtL_pole','F11xG14','EbarTxHtL','F11xHtL','G11xG14','HTxHT','EtL_polexcalE']
PALETTE=['#0072B2','#8C46A0','#D55E00','#009E73','#CC79A7','#996515','#56B4E9','#7A8C00','#A63150','#376795','#2A8C82','#E69F00']
FILES=[];CLIPPED=[]

def style():
 plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral','DejaVu Serif'],'mathtext.fontset':'stix','font.size':11,'axes.labelsize':14,'axes.titlesize':13,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':9,'axes.linewidth':.85,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'savefig.facecolor':'white'})

def state_data(data,species,obs,xi):return sorted([r for r in data if r['channel']==species and r['observable']==obs and r['x_index']==xi and r['status']=='measured'],key=lambda r:r['x'])

def term_style(term,allterms):
 if term in ['other','remainder']:return '#777777',(0,(2,2))
 index=allterms.index(term);return PALETTE[index%len(PALETTE)],['solid',(0,(6,2)),(0,(3,1,1,1)),(0,(1,1)),(0,(7,1,1,1,1,1))][index//len(PALETTE)%5]

def components(block,leading=None):
 a=block['arrays'];terms=block['terms'];ranked=sorted([t for t in terms if t!='other'],key=lambda t:-np.max(np.abs(a[t])))
 if leading is None:return [(t,a[t]) for t in terms]
 chosen=block.get('leading_terms',ranked)[:leading];remain=sum(a[t] for t in terms if t not in chosen)
 return [(t,a[t]) for t in chosen]+[('remainder',remain)]

def limits(theory,data,species,obs,decomp=False):
 if species=='pim' and obs=='lu1':return (-.2,.3)
 if species=='pim' and obs=='ul2':return (-1.,1.)
 return {'u1':(-.4,.8),'u2':(-.8,.6),'lu1':(-.2,.3),'ul1':(-.3,.3),'ul2':(-.4,.1),'ll0':(-.2,1.),'ll1':(-.4,.8)}[obs]

def axis(ax,obs,ylim):
 ax.set_facecolor('none');ax.set(xlim=(0,1.26),ylim=ylim,xlabel=r'$-t^\prime\;(\mathrm{GeV}^2)$',ylabel=LABELS[obs])
 ax.set_xticks(np.arange(0,1.21,.2));ax.set_yticks(np.round(np.arange(ylim[0],ylim[1]+1e-8,(.25 if obs=='ul2' and ylim[1]-ylim[0]>1.2 else {'ul1':.1,'ll0':.2,'ll1':.2,'lu1':.1,'u2':.2,'u1':.2,'ul2':.1}[obs])),10));ax.xaxis.set_minor_locator(AutoMinorLocator(2));ax.yaxis.set_minor_locator(AutoMinorLocator(2));ax.axhline(0,color='.65',lw=.75,zorder=1);ax.grid(alpha=.10,lw=.55)

def points(ax,rows,variant='stat_syst',color='black',marker='o'):
 if not rows:return
 x=[r['x'] for r in rows];y=[r['value'] for r in rows]
 e=[np.hypot(r['stat'],r['syst']) if variant=='combined_errors' else r['stat'] for r in rows]
 if variant=='stat_syst':
  for r in rows:
   if r['syst']>0:ax.add_patch(Rectangle((r['x']-.025,r['value']-r['syst']),.05,2*r['syst'],facecolor=color,edgecolor='none',alpha=.15,zorder=2))
 ax.errorbar(x,y,yerr=e,fmt=marker,color=color,mfc='white',mec=color,mew=1.05,markersize=5.4,capsize=2.1,elinewidth=1.05,linestyle='none',zorder=9)

def model(ax,block,decomp,allterms,color='#6A3D9A',linestyle='-',leading=None):
 a=block['arrays'];handles=[]
 if decomp:
  for term,y in components(block,leading):
   c,ls=block['styles'][term];line,=ax.plot(a['minus_tprime'],y,color=c,linestyle=ls,lw=2.0 if leading is not None else 1.75,zorder=3,label=label(term));handles.append(line)
 line,=ax.plot(a['minus_tprime'],a['total'],color='black' if decomp else color,lw=2.7 if decomp else 2.5,linestyle=linestyle,zorder=7,label='BM total');return [line]+handles

def watermark(fig,alpha):
 # 90% transparency. Draw behind axes and data.
 fig.text(.5,.51,'PRELIMINARY',ha='center',va='center',rotation=25,fontsize=fig.get_figwidth()*7.2,fontweight='bold',color='#D7191C',alpha=alpha,zorder=.05)

def save(fig,path,args,title,category):
 path=figure_path(args.output,path)
 path.parent.mkdir(parents=True,exist_ok=True);watermark(fig,args.watermark_alpha);fig.savefig(path,dpi=args.dpi);plt.close(fig)
 FILES.append(dict(path=str(path.relative_to(args.output)),title=title,category=category))

def error_handles(variant,species):
 result=[Line2D([],[],color='black',marker='o',mfc='white',ls='none',label='Data (stat. + ptp syst.)' if variant=='combined_errors' else 'Data (stat.)')]
 if species=='pip' and variant=='stat_syst':result.append(Patch(facecolor='.4',alpha=.2,label='Point-to-point syst.'))
 return result

def multipanel(theory,data,species,variant,args,xi=None,decomp=False,leading=None,gk=None):
 allxb=xi is None;large=decomp;fig,axs=plt.subplots(3,3,figsize=(16,15.5) if large else (12.8,10.5));flat=axs.ravel();indices=list(range(4)) if allxb else [xi]
 for panel,obs in zip(PANELS,ORDER):
  ax=flat[panel];axis(ax,obs,args.ylimits[species,obs,decomp]);handles=[]
  for idx in indices:
   block=theory[species,obs][idx];color=XCOLORS[idx] if allxb else '#6A3D9A'
   handles=model(ax,block,decomp,args.terms,color,'--' if gk is not None else '-',leading)
   if gk is not None and obs in ['lu1','ul1','ul2','ll0','ll1']:
    v=gk[np.isclose(gk[:,0],[.2,.3,.4,.5][idx])];ax.plot(v[:,1],v[:,5+['lu1','ul1','ul2','ll0','ll1'].index(obs)],color='black',lw=1.7,zorder=6)
   points(ax,state_data(data,species,obs,idx),variant,color=XCOLORS[idx] if allxb else 'black',marker=MARKERS[idx] if allxb else 'o')
  if species=='pim' and obs in ['u1','u2']:ax.text(.04,.91,'UU fixed to zero in the fit\n(no measured UU points)',transform=ax.transAxes,va='top',fontsize=10,bbox=dict(facecolor='white',edgecolor='none',alpha=.92))
  if decomp:ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.28),ncol=3,fontsize=8.1,frameon=False,columnspacing=.9,handlelength=2.7)
 flat[2].axis('off');flat[8].axis('off')
 if allxb:
  hs=[Line2D([],[],color=XCOLORS[i],marker=MARKERS[i],mfc='white',label=rf'${lo:.2f}<x_B<{hi:.2f}$') for i,(lo,hi) in enumerate(XB)]
  flat[2].legend(handles=hs,loc='upper left',frameon=False,title=r'$x_B$ ranges')
 else:
  q=theory[species,'lu1'][xi]['arrays']['Q2'][0];lo,hi=XB[xi]
  flat[2].text(.02,.95,rf'${lo:.2f}<x_B<{hi:.2f}$'+'\n\n'+rf'BM: $x_B={[.2,.3,.4,.5][xi]:.2f}$'+'\n'+rf'$Q^2={q:.1f}\;\mathrm{{GeV}}^2$',va='top',fontsize=13)
 info=error_handles(variant,species)
 if not decomp:info+=[Line2D([],[],color='#6A3D9A',lw=2.3,ls='--' if gk is not None else '-',label='BM total (updated)')]
 if gk is not None:info+=[Line2D([],[],color='black',lw=1.7,label='Goloskokov–Kroll')]
 flat[8].legend(handles=info,loc='upper left',frameon=False,fontsize=11)
 flat[8].text(.02,.51,'UL/LL data: beam-axis lab definition\nTheory supplied without frame conversion\n\n'+('LU: six −t′ bins\nUL/LL: two −t′ bins\nUU fixed to zero' if species=='pim' else 'Correlated scale errors not shown')+'\n\nCFF pieces share the total denominator.\nNo updated theory bands supplied.',va='top',fontsize=10.5,linespacing=1.45)
 title=CHANNEL[species]+' · '+('CFF decomposition' if decomp else 'Updated BM comparison')
 fig.suptitle(title,fontsize=19,y=.98)
 fig.subplots_adjust(left=.075,right=.975,bottom=.15 if large else .08,top=.91,hspace=.9 if large else .35,wspace=.32)
 kind=('cff_all' if leading is None else 'cff_leading5') if decomp else ('bm_and_gk' if gk is not None else 'totals')
 name='all_xB' if allxb else f'xB_bin_{xi+1}'
 path=args.output/species/variant/kind/(name+'.png');save(fig,path,args,title+' '+name,species+' '+kind)

def focused_limits(species,obs,decomp,base):
 return base

def focused(theory,data,species,obs,variant,args,decomp=False,leading=None):
 fig,axs=plt.subplots(2,2,figsize=(12.8,10.5) if decomp else (11.7,8.6));handles=None
 for xi,ax in enumerate(axs.flat):
  axis(ax,obs,focused_limits(species,obs,decomp,args.ylimits[species,obs,decomp]));lo,hi=XB[xi]
  ax.set_title(rf'${lo:.2f}<x_B<{hi:.2f}$'+'  ·  '+rf'$Q^2={theory[species,obs][xi]["arrays"]["Q2"][0]:.1f}\;\mathrm{{GeV}}^2$')
  handles=model(ax,theory[species,obs][xi],decomp,args.terms,leading=leading)
  points(ax,state_data(data,species,obs,xi),variant)
  if species=='pim' and obs in ['u1','u2']:ax.text(.04,.90,'Fixed UU input, no measured points',transform=ax.transAxes,fontsize=10)
 extra=error_handles(variant,species);bylabel={h.get_label():h for h in (handles or [])+extra}
 fig.legend(handles=list(bylabel.values()),loc='lower center',bbox_to_anchor=(.5,.014),ncol=3 if decomp else 3,frameon=False,fontsize=11 if decomp and leading is not None else 10.1,columnspacing=1.9,handlelength=3.0)
 fig.suptitle(CHANNEL[species]+'  ·  '+LABELS[obs]+('  ·  Every supplied CFF product' if decomp and leading is None else '  ·  Leading CFF products + remainder' if decomp else ''),fontsize=18,y=.985)
 fig.subplots_adjust(left=.09,right=.98,bottom=(.205 if leading is not None else .25) if decomp else .13,top=.91,hspace=.32,wspace=.23)
 kind='cff_all' if decomp and leading is None else 'cff_leading5' if decomp else 'totals'
 path=args.output/species/variant/'focused'/f'{obs}_{kind}.png';save(fig,path,args,CHANNEL[species]+' '+obs+' '+kind,species+' focused')

def lu_charge_comparison(theory,data,args,light=False,include_pi0=False,residual=False,combined_errors=False):
 fig,axs=plt.subplots(2,2,figsize=(11.7,8.7));lim=max(args.ylimits[s,'lu1',False][1] for s in ['pip','pim'])
 for xi,ax in enumerate(axs.flat):
  axis(ax,'lu1',(-.2,lim));lo,hi=XB[xi];ax.set_title(rf'${lo:.2f}<x_B<{hi:.2f}$')
  for species,c,m,ls in [('pip','#0072B2','o','-'),('pim','#D55E00','s','--')]:
   a=theory[species,'lu1'][xi]['arrays'];rr=state_data(data,species,'lu1',xi)
   if residual:
    rr=[dict(r,value=r['value']-float(np.interp(r['x'],a['minus_tprime'],a['total']))) for r in rr if a['minus_tprime'][0]<=r['x']<=a['minus_tprime'][-1]]
   else:ax.plot(a['minus_tprime'],a['total'],c=c,ls=ls,lw=2.2,zorder=5)
   points(ax,rr,'combined_errors' if combined_errors else 'stat_syst' if species=='pip' else 'stat',c,m)
  if include_pi0:
   a=theory['pi0','lu1'][xi]['arrays'];ax.plot(a['minus_tprime'],a['total'],color='.35',ls=':',lw=1.8,label=r'$\pi^0$ BM (no data)')
  if residual:ax.set_ylabel('Data − BM total')
 handles=[Line2D([],[],color=c,marker=m,mfc='white',ls=ls,label=CHANNEL[s]+(' data (stat. + ptp syst.) / BM total' if combined_errors and s=='pip' else ' data / BM total')) for s,c,m,ls in [('pip','#0072B2','o','-'),('pim','#D55E00','s','--')]]
 if include_pi0:handles.append(Line2D([],[],color='.35',ls=':',label=r'$\pi^0$ BM (no data)'))
 if not combined_errors:handles.append(Patch(facecolor='#0072B2',alpha=.15,label=r'$\pi^+$ point-to-point syst.'))
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.025),ncol=2,frameon=False,fontsize=11)
 fig.suptitle(LABELS['lu1']+' · '+('Data − BM at the plotted mean −t′' if residual else 'Charged-pion comparison'),fontsize=19,y=.98)
 fig.subplots_adjust(left=.095,right=.98,bottom=.16,top=.90,wspace=.23,hspace=.32)
 suffix='_residuals' if residual else '_with_pi0' if include_pi0 else ''
 if combined_errors:suffix+='_combined_errors'
 original=args.watermark_alpha
 if light:args.watermark_alpha=.10;suffix+='_light_watermark'
 save(fig,args.output/'focused'/('LU_pip_pim'+suffix+'.png'),args,'LU charged-pion comparison'+suffix,'Cross-channel LU')
 args.watermark_alpha=original

def gallery(args):
 write_gallery(args.output,FILES)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--theory',type=Path,default=HERE/'inputs');p.add_argument('--measurements',type=Path,default=HERE/'measurements/measurements.csv');p.add_argument('--output',type=Path,default=HERE/'output');p.add_argument('--watermark-alpha',type=float,default=.10);p.add_argument('--dpi',type=int,default=300);p.add_argument('--channels',nargs='+',choices=['pip','pim'],default=['pip','pim']);p.add_argument('--quick',action='store_true');a=p.parse_args()
 if not 0<=a.watermark_alpha<=1:p.error('Watermark alpha must be 0-1')
 FILES.clear()
 style();theory,audit=load_theory(a.theory);data=load_measurements(a.measurements);a.output.mkdir(parents=True,exist_ok=True)
 diagnostics=a.output/'data/diagnostics';diagnostics.mkdir(parents=True,exist_ok=True)
 for blocks in theory.values():
  ranked=sorted([t for t in blocks[0]['terms'] if t!='other'],key=lambda t:-max(np.max(np.abs(b['arrays'][t])) for b in blocks))
  for block in blocks:block['leading_terms']=ranked
 styles=configure_styles(theory)
 (diagnostics/'cff_colors.json').write_text(json.dumps(styles,indent=2)+'\n')
 terms={t for blocks in theory.values() for b in blocks for t in b['terms']};a.terms=[t for t in IMPORTANT if t in terms]+sorted(terms-set(IMPORTANT)-{'other'})
 a.ylimits={(s,o,d):limits(theory,data,s,o,d) for s in ['pip','pim'] for o in ORDER for d in [False,True]}
 # LU axes are shared across the charge-comparison variants.
 (diagnostics/'theory_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
 gk=None
 if (a.theory/'gk_pip.csv').exists():gk=np.loadtxt(a.theory/'gk_pip.csv',delimiter=',',skiprows=1)
 if a.quick:
  focused(theory,data,'pip','ul1','stat_syst',a,True);focused(theory,data,'pip','ul2','stat_syst',a,True)
  focused(theory,data,'pip','ul1','stat_syst',a,True,5);focused(theory,data,'pip','ul2','stat_syst',a,True,5)
  lu_charge_comparison(theory,data,a)
  multipanel(theory,data,'pip','stat_syst',a,xi=0,decomp=True)
 else:
  for species in a.channels:
   variants=['stat_syst','combined_errors'] if species=='pip' else ['stat']
   for variant in variants:
    for xi in range(4):
     multipanel(theory,data,species,variant,a,xi=xi)
     if species=='pip' and gk is not None:multipanel(theory,data,species,variant,a,xi=xi,gk=gk)
    multipanel(theory,data,species,variant,a)
    for obs in ORDER:focused(theory,data,species,obs,variant,a)
   variant=variants[0]
   for xi in range(4):
    multipanel(theory,data,species,variant,a,xi=xi,decomp=True)
    multipanel(theory,data,species,variant,a,xi=xi,decomp=True,leading=5)
   for obs in ORDER:
    focused(theory,data,species,obs,variant,a,decomp=True)
    focused(theory,data,species,obs,variant,a,decomp=True,leading=5)
  lu_charge_comparison(theory,data,a);lu_charge_comparison(theory,data,a,include_pi0=True);lu_charge_comparison(theory,data,a,residual=True)
 gallery(a)
 manifest=dict(plots=FILES,number_of_plots=len(FILES),watermark_opacity=a.watermark_alpha,watermark_preview_opacity=.10,measurements=str(a.measurements.resolve()),measurement_sha256=hashlib.sha256(a.measurements.read_bytes()).hexdigest(),theory_files=audit,code_sha256=hashlib.sha256(Path(__file__).read_bytes()+(HERE/'theory_io.py').read_bytes()+(HERE/'cff_style.py').read_bytes()+(HERE/'figure_layout.py').read_bytes()).hexdigest(),theory_bands='No updated uncertainties supplied. Old bands are not reused.',frame_convention='UL/LL data retain lab-axis labels. Supplied predictions plotted without frame conversion; no assertion that their target axis has changed.',curves='Supplied points connected piecewise linearly. No smoothing, extrapolation, bin averaging, sign or depolarization rescaling.',decomposition='All supplied CFF products plus supplied other term. Leading-five displays sum all omitted supplied terms into a separately labeled remainder.',scale_errors='Not included.',limits={':'.join([s,o,str(d)]):v for (s,o,d),v in a.ylimits.items()})
 (diagnostics/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'Created {len(FILES)} PNG plots in {a.output}',flush=True)
if __name__=='__main__':main()
