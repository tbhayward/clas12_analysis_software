"""Group supplied CFF-product contributions by added chiral-even content.

All subtotals retain the supplied FULL BM denominator. This is not a new fit,
not a GK curve, and not a calculation with amplitudes switched off. 'other'
remains unknown because its individual products were not supplied.
"""
import csv,json,hashlib,argparse,shutil
from pathlib import Path
from types import SimpleNamespace
import numpy as np
ADDED=frozenset(['F11','G11','F12','G12','F14','G14','calE','calEt'])
CHIRAL_ODD=frozenset(['HT','EbarT'])
LONGSTANDING=frozenset(['HtL','EtL_np','EtL_pole'])
ESTABLISHED=CHIRAL_ODD|LONGSTANDING
COLORS={'added':'#1f77b4','chiral_odd':'#d62728','longitudinal_pole':'#ff7f0e','unknown':'#777777'}
LABELS={'added':'BM, Chiral-even transverse CFFs + interf.','chiral_odd':'BM, Chiral-odd CFFs + interf.','longitudinal_pole':r'BM, $\widetilde H_L$, $\widetilde E_L$ and pion-pole products','unknown':'Unclassified supplied “other”'}

def classify(term):
 if term=='other':return 'unknown'
 tokens=term.split('x')
 if len(tokens)!=2 or not set(tokens)<=ADDED|ESTABLISHED:raise ValueError('Cannot classify '+term)
 if any(t in ADDED for t in tokens):return 'added'
 return 'chiral_odd' if any(t in CHIRAL_ODD for t in tokens) else 'longitudinal_pole'

def grouped(block):
 result={k:np.zeros_like(block['arrays']['total']) for k in COLORS}
 for term in block['terms']:result[classify(term)]+=block['arrays'][term]
 return result

def draw(theory,data,species,obs,output,dpi,alpha,individual=False,gk=None,shared_variant=None):
 import matplotlib.pyplot as plt
 from matplotlib.lines import Line2D
 import plot_comparisons as p
 from theory_io import label
 fig,axes=plt.subplots(2,2,figsize=(12.4,9.8 if individual else 9.2))
 line_handles=[]
 for xi,ax in enumerate(axes.flat):
  b=theory[species,obs][xi];a=b['arrays'];x=a['minus_tprime'];g=grouped(b)
  limit=p.focused_limits(species,obs,True,p.limits(theory,data,species,obs,True))
  p.axis(ax,obs,limit);lo,hi=p.XB[xi]
  ax.set_title(rf'${lo:.2f}<x_B<{hi:.2f}$'+'  ·  '+rf'$Q^2={a["Q2"][0]:.1f}\;\mathrm{{GeV}}^2$')
  line_handles=[]
  if individual:
   counts={k:0 for k in COLORS};patterns=['-',(0,(6,2)),(0,(3,1,1,1)),(0,(1,1)),(0,(7,1,1,1,1,1)),(0,(3,2)),(0,(8,3)),(0,(2,1,2,3)),(0,(6,1,1,1,1,1,1,1)),(0,(2,2,1,2))]
   for term in b['terms']:
    if term=='other':continue
    sector=classify(term);ls=patterns[counts[sector]%len(patterns)];counts[sector]+=1
    line,=ax.plot(x,a[term],color=COLORS[sector],ls=ls,lw=1.65,label=label(term),zorder=3);line_handles.append(line)
  else:
   for sector in ['longitudinal_pole','chiral_odd','added']:
    ax.plot(x,g[sector],color=COLORS[sector],lw=2.15,ls='--' if sector=='unknown' else '-',label=LABELS[sector],zorder=4)
  ax.plot(x,a['total'],color='black',lw=2.7,zorder=7)
  if shared_variant=='bm+gk':
   if gk is not None and species=='pip' and obs in ['lu1','ul1','ul2','ll0','ll1']:
    v=gk[np.isclose(gk[:,0],[.2,.3,.4,.5][xi])];v=v[np.argsort(v[:,1])]
    ax.plot(v[:,1],v[:,5+['lu1','ul1','ul2','ll0','ll1'].index(obs)],color='#b0b0b0',lw=2.3,zorder=5)
  p.points(ax,p.state_data(data,species,obs,xi),'stat_syst' if species=='pip' else 'stat')
  if species=='pim' and obs in ['u1','u2']:ax.text(.03,.92,'UU fixed in fit; no measured points',transform=ax.transAxes,fontsize=9)
 total=Line2D([],[],color='black',lw=2.7,label='BM Total')
 sectors=[Line2D([],[],color=COLORS[k],lw=2.15,ls='--' if k=='unknown' else '-',label=LABELS[k]) for k in ['added','chiral_odd','longitudinal_pole']]
 handles=[total]+line_handles+p.error_handles('stat_syst' if species=='pip' else 'stat',species) if individual else [total]+sectors+p.error_handles('stat_syst' if species=='pip' else 'stat',species)
 if shared_variant=='bm+gk' and gk is not None and species=='pip' and obs in ['lu1','ul1','ul2','ll0','ll1']:handles.insert(0,Line2D([],[],color='#b0b0b0',lw=2.3,label='GK Total'))
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.04),ncol=3 if individual else 2,frameon=False,fontsize=9.6 if individual else 10.4,handlelength=3.)
 fig.suptitle(p.CHANNEL[species]+'  ·  '+p.LABELS[obs]+('  ·  CFF products by sector' if individual else ''),fontsize=17,y=.98)
 fig.subplots_adjust(left=.09,right=.98,bottom=.26 if individual else .19,top=.90,hspace=.32,wspace=.24)
 p.watermark(fig,alpha)
 from figure_layout import OBS_NAMES
 channel='pi_plus' if species=='pip' else 'pi_minus'
 folder=output/'figures'/channel/('featured_stories' if species=='pip' else 'cff_breakdowns/chiral_sectors')
 if shared_variant:folder=output/'figures/shared_plots'/shared_variant
 folder.mkdir(parents=True,exist_ok=True)
 suffix='CFF_products_by_chiral_sector' if individual else 'chiral_sector_sums'
 path=folder/f'{channel}_{OBS_NAMES[obs]}_{suffix}_2x2.png'
 fig.savefig(path,dpi=dpi);plt.close(fig)
 return dict(path=str(path.relative_to(output)),title=channel+' '+obs+' '+suffix,category=channel+' chiral sectors')

def run(theory_path,measurements,output,dpi=300,alpha=.10,channels=('pip','pim')):
 if (Path(theory_path)/'photon_referenced').exists():
  from step_predictions import run as run_steps
  return run_steps(theory_path,measurements,output,dpi,alpha,channels)
 from theory_io import load_theory,load_measurements
 import plot_comparisons as p
 p.style();output=Path(output);theory,audit=load_theory(theory_path);data=load_measurements(measurements)
 entries=[];classifications=[];curves=[]
 for (species,obs),blocks in theory.items():
  for term in blocks[0]['terms']:
   classifications.append(dict(species=species,observable=obs,term=term,sector=classify(term),added_CFFs=';'.join(t for t in term.split('x') if t in ADDED)))
  for xi,b in enumerate(blocks):
   g=grouped(b);a=b['arrays'];residual=sum(g.values())-a['total']
   if np.max(np.abs(residual))>5e-6*(len(b['terms'])+1)+1e-9:raise ValueError('Sector subtotals fail closure')
   for i,x in enumerate(a['minus_tprime']):curves.append(dict(species=species,observable=obs,xB=float(a['xB'][i]),Q2=float(a['Q2'][i]),minus_tprime=float(x),full_total=float(a['total'][i]),with_added_chiral_even_CFF=float(g['added'][i]),longstanding_longitudinal_pole=float(g['longitudinal_pole'][i]),chiral_odd_and_interference=float(g['chiral_odd'][i]),unclassified_other=float(g['unknown'][i]),rounding_residual=float(residual[i])))
  if species in channels:
   entries.append(draw(theory,data,species,obs,output,dpi,alpha))
   if species=='pip' and obs in ['ul1','ul2']:entries.append(draw(theory,data,species,obs,output,dpi,alpha,individual=True))
 from shared_plots import build
 entries.extend(build(theory,data,theory_path,output,dpi,alpha,channels))
 diagnostics=output/'data/diagnostics';diagnostics.mkdir(parents=True,exist_ok=True)
 for name,rows in [('chiral_sector_classification.csv',classifications),('chiral_sector_curves.csv',curves)]:
  with (diagnostics/name).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 info=dict(added_CFFs=sorted(ADDED),chiral_odd_CFFs=sorted(CHIRAL_ODD),longstanding_longitudinal_pole_CFFs=sorted(LONGSTANDING),classification='First assign any product containing an added F/G transverse CFF or its calE/calEt kernel partner to the added sector, including interference. Of the remaining products, assign those containing HT/EbarT to the chiral-odd sector. The rest involve only longstanding HtL/EtL/pion-pole CFFs. Ordinary chiral-even longitudinal GPDs are not new.',unknown='The supplied other column is not classified and is hidden in figures. It is retained in the full total and numerical audits. It is not rounding error.',normalization='Every curve retains the full supplied BM denominator. No subtotal is an old-model or switched-off-amplitude prediction. calE/calEt are classified by the supplied transverse F11/G14 kernels, not simply because they are chiral even.',source_theory='User-supplied Brandon Manley theory section, equations defining A,B,C,D,H,E and structure functions.',source_notation='Supplied pip/notation.md defines calE as GPD E with F11 kernel, and calEt as non-pole Etilde with G14 kernel.',plots=entries,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (diagnostics/'chiral_sector_definition.json').write_text(json.dumps(info,indent=2)+'\n')
 manifest_path=diagnostics/'manifest.json'
 if manifest_path.exists():
  manifest=json.loads(manifest_path.read_text());lookup={r['path']:r for r in manifest['plots'] if (output/r['path']).exists()};lookup.update({r['path']:r for r in entries});manifest['plots']=list(lookup.values());manifest['number_of_plots']=len(lookup);manifest['chiral_sector_definition']=info;manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
 print('Created '+str(len(entries))+' chiral-sector comparison figures',flush=True)
 return entries

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--theory',type=Path,required=True);parser.add_argument('--measurements',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--dpi',type=int,default=300);a=parser.parse_args();run(a.theory,a.measurements,a.output,a.dpi)
