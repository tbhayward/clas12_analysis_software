"""Bin-by-bin model-input shifts and compact PNG comparisons."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
NAMES=['UU_cos_phi','UU_cos_2phi','LU_sin_phi','UL_lab_sin_phi','UL_lab_sin_2phi','LL_lab','LL_lab_cos_phi']
MODELS=['zero','step1_2','step1_2_3'];LABELS=['UU = 0','UU = BM longitudinal/pole + chiral-odd','UU = full BM']
COLORS=['#222222','#d62728','#1f77b4'];MARKERS=['o','s','^']
TEX=[r'$F_{LU}^{\sin\phi}/F_{UU}$',r'$F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$',r'$F_{UL,\mathrm{lab}}^{\sin2\phi}/F_{UU}$',r'$F_{LL,\mathrm{lab}}/F_{UU}$',r'$F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$']
XB=[(.1,.25),(.25,.35),(.35,.45),(.45,.6)]

def summarize(out,results,means,baseline_fine,baseline_coarse):
 out=Path(out);figures=out/'figures';figures.mkdir(exist_ok=True)
 lookup={(r['scheme'],r['model'],r['bin']):r for r in results}
 rows=[];checks=[];aggregate=[]
 for scheme,n,baseline in [('fine',24,baseline_fine),('coarse',8,baseline_coarse)]:
  for b in range(1,n+1):
   zero=lookup[scheme,'zero',b]
   if baseline:
    path=Path(baseline)/'fits'/f'combined_bin{b:02d}.json'
    if path.exists():
     old=json.loads(path.read_text());check=dict(scheme=scheme,bin=b,old_status=old['status'],new_status=zero['status'])
     if old.get('status')=='valid' and zero['status']=='valid':
      check.update(events_match=old['events']==zero['events'],max_absolute_value_difference=float(np.max(np.abs(np.array(old['values'])-zero['values']))),max_absolute_error_difference=float(np.max(np.abs(np.array(old['errors'])-zero['errors']))),nll_difference=zero['nll']-old['nll'],max_relative_error_difference=float(np.max(np.abs(np.array(zero['errors'][2:])/np.array(old['errors'][2:])-1))))
     checks.append(check)
   for model in MODELS:
    r=lookup[scheme,model,b]
    for j in range(2,7):
     row=dict(scheme=scheme,model=model,bin=b,x_index=(b-1)//(6 if scheme=='fine' else 2),mean_minus_tprime=float(means[scheme][b-1]),parameter=NAMES[j],status=r['status'],events=r.get('events'),UU_cos_phi=r['UU_input']['values'][0],UU_cos_2phi=r['UU_input']['values'][1])
     if r['status']=='valid':
      row.update(value=r['values'][j],stat=r['errors'][j],nll=r['nll'],minimum_state_factor=r['minimum_state_factor'],warnings='|'.join(r['warnings']))
      if zero['status']=='valid':row.update(zero_value=zero['values'][j],zero_stat=zero['errors'][j],shift=r['values'][j]-zero['values'][j],shift_over_zero_stat=(r['values'][j]-zero['values'][j])/zero['errors'][j],error_ratio=r['errors'][j]/zero['errors'][j])
     else:row.update(error=r.get('error',''))
     rows.append(row)
  for j in range(2,7):
   for model in MODELS[1:]:
    rr=[r for r in rows if r['scheme']==scheme and r['model']==model and r['parameter']==NAMES[j] and 'shift' in r]
    if rr:
     worst=max(rr,key=lambda r:abs(r['shift_over_zero_stat']))
     aggregate.append(dict(scheme=scheme,model=model,parameter=NAMES[j],compared_bins=len(rr),median_absolute_shift=float(np.median([abs(r['shift']) for r in rr])),max_absolute_shift=max(abs(r['shift']) for r in rr),median_absolute_shift_over_zero_stat=float(np.median([abs(r['shift_over_zero_stat']) for r in rr])),max_absolute_shift_over_zero_stat=abs(worst['shift_over_zero_stat']),worst_bin=worst['bin'],median_error_ratio=float(np.median([r['error_ratio'] for r in rr]))))
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (out/'comparison.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 with (out/'shift_summary.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(aggregate[0]));w.writeheader();w.writerows(aggregate)
 (out/'baseline_reproduction.json').write_text(json.dumps(checks,indent=2)+'\n')
 status={scheme:{m:{state:sum(r['scheme']==scheme and r['model']==m and r['status']==state for r in results) for state in ['valid','invalid','failed','empty']} for m in MODELS} for scheme in ['fine','coarse']}
 summary=dict(status=status,shifts=aggregate,baseline_checks=checks,interpretation='Shifts use the same events and are correlated. Shift divided by zero-scenario statistical error is a sensitivity scale, not a significance. Fixed UU theory uncertainties are not propagated. Coarse dilution factors are reused unchanged. Do not compare absolute likelihoods between different fixed-UU scenarios as a goodness-of-fit test: UU-only terms are conditional-likelihood normalization directions.',convergence=[dict(scheme=r['scheme'],model=r['model'],bin=r['bin'],status=r['status'],warnings=r.get('warnings',[]),minimum_state_factor=r.get('minimum_state_factor'),error=r.get('error')) for r in results])
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 plt.rcParams.update({'font.family':'serif','mathtext.fontset':'stix','font.size':11})
 plot_specs=[('fine',2)]+[('coarse',j) for j in range(3,7)]
 for scheme,j in plot_specs:
  nb=6 if scheme=='fine' else 2
  for shifts in [False,True]:
   fig,axes=plt.subplots(2,2,figsize=(12,9));allvalues=[]
   for xi,ax in enumerate(axes.flat):
    for mi,model in enumerate(MODELS):
     if shifts and model=='zero':continue
     selected=[r for r in rows if r['scheme']==scheme and r['model']==model and r['parameter']==NAMES[j] and r['x_index']==xi and r['status']=='valid']
     if shifts:
      selected=[r for r in selected if 'shift_over_zero_stat' in r];x=[r['mean_minus_tprime'] for r in selected];y=[r['shift_over_zero_stat'] for r in selected]
      ax.plot(x,y,marker=MARKERS[mi],ls='none',color=COLORS[mi],ms=6,label=LABELS[mi]);allvalues+=y
     else:
      x=[r['mean_minus_tprime']+(mi-1)*.012 for r in selected];y=[r['value'] for r in selected];err=[r['stat'] for r in selected]
      ax.errorbar(x,y,yerr=err,fmt=MARKERS[mi],ls='none',color=COLORS[mi],mfc='white',ms=5,capsize=2,label=LABELS[mi]);allvalues.extend([v+e for v,e in zip(y,err)]);allvalues.extend([v-e for v,e in zip(y,err)])
    lo,hi=XB[xi];ax.set_title(rf'${lo:.2f}<x_B<{hi:.2f}$');ax.set(xlim=(0,1.26),xlabel=r'$-t^\prime\;(\mathrm{GeV}^2)$',ylabel=r'$(R_{\rm model\ UU}-R_{\rm zero\ UU})/\sigma_{\rm zero}$' if shifts else TEX[j-2]);ax.axhline(0,c='.6',lw=.8);ax.grid(alpha=.12)
   if shifts:
    lim=max(.1,max(map(abs,allvalues),default=.1)*1.15)
    for ax in axes.flat:ax.set_ylim(-lim,lim)
   else:
    lo=max(-1.,min(allvalues,default=-.1)-.05);hi=min(1.,max(allvalues,default=.1)+.05)
    if scheme=='fine' and j==2:lo,hi=-.2,.3
    if lo>=hi:lo,hi=-1,1
    for xi,ax in enumerate(axes.flat):
     ax.set_ylim(lo,hi)
     for mi,model in enumerate(MODELS):
      for row in rows:
       if row['scheme']==scheme and row['model']==model and row['parameter']==NAMES[j] and row['x_index']==xi and row['status']=='valid' and (row['value']>hi or row['value']<lo):
        yy=hi-.015*(hi-lo) if row['value']>hi else lo+.015*(hi-lo)
        ax.plot(row['mean_minus_tprime']+(mi-1)*.012,yy,marker='^' if row['value']>hi else 'v',color=COLORS[mi],ms=6)
        ax.annotate(f"{row['value']:.2f}",(row['mean_minus_tprime']+(mi-1)*.012,yy),xytext=(0,-12-11*mi if row['value']>hi else 10+11*mi),textcoords='offset points',ha='center',fontsize=8,color=COLORS[mi])
   handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=1,frameon=False,bbox_to_anchor=(.5,.01),fontsize=10)
   fig.suptitle('π− fixed-UU sensitivity · '+TEX[j-2],fontsize=17);fig.subplots_adjust(left=.10,right=.98,top=.91,bottom=.17,wspace=.25,hspace=.35)
   fig.savefig(figures/f'{scheme}_{NAMES[j]}_{"shifts" if shifts else "comparison"}.png',dpi=250);plt.close(fig)
 lines=['Fixed-UU model-input study','',json.dumps(status,indent=2),'','Shifts are correlated refits, not independent measurements. No theory-UU uncertainty is included.','']
 for r in aggregate:
  if (r['scheme']=='fine' and r['parameter']=='LU_sin_phi') or (r['scheme']=='coarse' and r['parameter']!='LU_sin_phi'):
   lines.append(f"{r['scheme']} {r['parameter']} {r['model']}: median |shift|/sigma_zero={r['median_absolute_shift_over_zero_stat']:.3f}, maximum={r['max_absolute_shift_over_zero_stat']:.3f} (bin {r['worst_bin']}); median error ratio={r['median_error_ratio']:.3f}")
 (out/'assessment.txt').write_text('\n'.join(lines)+'\n');print('\n'.join(lines),flush=True)
 from theory_check import make
 make(out)
