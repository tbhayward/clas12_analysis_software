"""Rank explicitly supplied F14/G11 products, not full parameter sensitivity."""
import csv,json,argparse,hashlib
from pathlib import Path
import numpy as np
TARGETS=['F14','G11']

def integral_abs(x,y):
 """Exact integral of |piecewise-linear y|, including zero crossings."""
 dx=np.diff(x);a=y[:-1];b=y[1:];v=.5*dx*(np.abs(a)+np.abs(b));cross=a*b<0
 v[cross]=.5*dx[cross]*(a[cross]**2+b[cross]**2)/(np.abs(a[cross])+np.abs(b[cross]))
 return float(v.sum())

def restrict(a,low=.05,high=1.2):
 original=a['minus_tprime'];low=max(low,float(original[0]));high=min(high,float(original[-1]));x=np.r_[low,original[(original>low)&(original<high)],high]
 return x,{k:np.interp(x,original,v) for k,v in a.items()}

def run(theory_path,measurements,output):
 from theory_io import load_theory,load_measurements
 from figure_layout import OBS_TEXT
 theory,_=load_theory(theory_path);data=load_measurements(measurements);output=Path(output)
 folder=output/'data/diagnostics/cff_sensitivity';folder.mkdir(parents=True,exist_ok=True)
 rows=[];overall=[]
 for (species,obs),blocks in theory.items():
  for target in TARGETS:
   binrows=[]
   for xi,block in enumerate(blocks):
    x,a=restrict(block['arrays']);terms=[t for t in block['terms'] if target in t.split('x')]
    contribution=sum((a[t] for t in terms),start=np.zeros_like(x));net=integral_abs(x,contribution);gross=sum(integral_abs(x,a[t]) for t in terms);denom=integral_abs(x,a['total']);allgross=sum(integral_abs(x,a[t]) for t in block['terms'])
    ds=[r for r in data if r['channel']==species and r['observable']==obs and r['x_index']==xi and r['status']=='measured' and x[0]<=r['x']<=x[-1]]
    leverage=sum((np.interp(r['x'],x,contribution)/np.hypot(r['stat'],r['syst']))**2 for r in ds)
    row=dict(species=species,observable=obs,target=target,xB=[.2,.3,.4,.5][xi],Q2=float(a['Q2'][0]),minus_tprime_low=float(x[0]),minus_tprime_high=float(x[-1]),listed_products=';'.join(terms),has_explicit_products=bool(terms),net_abs_integral=net,gross_abs_integral=gross,total_abs_integral=denom,all_products_abs_integral=allgross,net_fraction_of_total=net/denom if denom>0 else None,gross_fraction_of_products=gross/allgross if allgross>0 else None,mean_abs_contribution=net/(x[-1]-x[0]),pointwise_leverage_squared=leverage,measured_points=len(ds))
    rows.append(row);binrows.append(row)
   net=sum(r['net_abs_integral'] for r in binrows);gross=sum(r['gross_abs_integral'] for r in binrows);denom=sum(r['total_abs_integral'] for r in binrows);allgross=sum(r['all_products_abs_integral'] for r in binrows)
   overall.append(dict(species=species,observable=obs,target=target,listed_products=';'.join(sorted({t for r in binrows for t in r['listed_products'].split(';') if t})),has_explicit_products=any(r['has_explicit_products'] for r in binrows),net_fraction_of_total=net/denom if denom>0 else None,gross_fraction_of_products=gross/allgross if allgross>0 else None,mean_abs_contribution=net/sum(r['minus_tprime_high']-r['minus_tprime_low'] for r in binrows),pointwise_leverage=float(np.sqrt(sum(r['pointwise_leverage_squared'] for r in binrows))) if sum(r['measured_points'] for r in binrows)>0 else None,measured_points=sum(r['measured_points'] for r in binrows)))
 for name,table in [('by_xB.csv',rows),('rankings.csv',overall)]:
  with (folder/name).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
 note=dict(range_minus_tprime=[.05,1.2],range_note='The selected analysis extends to 1.25 GeV2. The supplied theory ends at 1.20, so the last 0.05 is not extrapolated.',primary_metric='Sum over four equally weighted xB settings of integral |sum of explicitly listed target-CFF products| divided by sum of integral |full model total|. Common full denominator retained. Exact piecewise-linear absolute integration avoids zero-crossing divergences.',secondary_metric='Sum of integrals of absolute individual target products divided by the corresponding sum for all listed products including other. This shows cancellation-resistant composition.',leverage='sqrt(sum [(listed net contribution at measured mean -tprime)/sqrt(stat^2+ptp_syst^2)]^2). Correlations, denominator response and bin-folding are omitted. Descriptive leverage proxy, not significance or Fisher information.',limitations=['The other column can contain F14/G11 products that are not separately supplied. No explicit product is not evidence of zero sensitivity.','Full CFF variation changes both numerator and denominator. These tables do not supply its full derivative.','Ratios can exceed 100% because signed products cancel in the full total.','Underlying F14 and G11 parameter variations are not independent of the other CFFs in this model.'],results=overall)
 note['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 (folder/'analysis.json').write_text(json.dumps(note,indent=2)+'\n')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'serif','mathtext.fontset':'stix','font.size':11})
 fig,axes=plt.subplots(2,2,figsize=(12,8.5))
 for i,species in enumerate(['pip','pim']):
  for j,target in enumerate(TARGETS):
   ax=axes[i,j];r=sorted([r for r in overall if r['species']==species and r['target']==target],key=lambda r:-(r['net_fraction_of_total'] or 0));names=[OBS_TEXT[v['observable']] for v in r];values=[100*(v['net_fraction_of_total'] or 0) for v in r]
   ax.barh(names,values,color='#1f77b4' if target=='F14' else '#d62728');ax.invert_yaxis();ax.set_title(('π+' if species=='pip' else 'π−')+' · '+('$F_{1,4}$' if target=='F14' else '$G_{1,1}$'));ax.set_xlabel('Integrated listed contribution / total (%)');ax.grid(axis='x',alpha=.15)
   for k,v in enumerate(values):ax.text(v+max(values+[1])*.015,k,f'{v:.1f}%' if r[k]['has_explicit_products'] else 'not separately listed',va='center',fontsize=9 if r[k]['has_explicit_products'] else 8)
   ax.set_xlim(0,max(values+[1])*1.24)
 fig.suptitle('Relative size of the explicitly supplied CFF products',fontsize=16);fig.tight_layout(rect=(0,0,1,.95));plots=output/'figures/cff_sensitivity';plots.mkdir(parents=True,exist_ok=True);__import__('plot_comparisons').watermark(fig,.10);fig.savefig(plots/'F14_G11_listed_contribution_ranking.png',dpi=300);plt.close(fig)
 text='CFF contribution-size ranking\n\nRange: 0.05 <= -tprime <= 1.20 GeV2, four xB settings weighted equally.\nThis ranks explicitly supplied numerator contributions at the full-model denominator, not a complete CFF sensitivity or amplitudes-off prediction.\n\n'
 for species in ['pip','pim','pi0']:
  for target in TARGETS:
   r=sorted([r for r in overall if r['species']==species and r['target']==target and r['has_explicit_products']],key=lambda r:-r['net_fraction_of_total']);text+=species+' '+target+'\n'
   for v in r:text+=f"  {v['observable']}: relative contribution {100*v['net_fraction_of_total']:.1f}%, gross-product share {100*v['gross_fraction_of_products']:.1f}%, mean absolute contribution {v['mean_abs_contribution']:.4f}, pointwise leverage {v['pointwise_leverage'] if v['pointwise_leverage'] is not None else 'not measured'}\n"
 text+='\nNo listed product does not establish zero sensitivity. The truncated other contribution and full denominator response are unavailable.\n'
 (folder/'assessment.txt').write_text(text)
 print(text,flush=True)
 manifest_path=folder.parent/'manifest.json'
 if manifest_path.exists():
  manifest=json.loads(manifest_path.read_text());entry=dict(path='figures/cff_sensitivity/F14_G11_listed_contribution_ranking.png',title='F14 and G11 contribution-size ranking',category='CFF sensitivity')
  lookup={r['path']:r for r in manifest['plots']};lookup[entry['path']]=entry;manifest['plots']=list(lookup.values());manifest['number_of_plots']=len(lookup);manifest['cff_sensitivity']=note;manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
 return overall
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--theory',type=Path,required=True);p.add_argument('--measurements',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.theory,a.measurements,a.output)
