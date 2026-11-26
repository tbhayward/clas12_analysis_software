"""Show the UU-refit LU values alongside the corresponding BM LU curves."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make(out):
 from summarize import MODELS,COLORS,MARKERS,XB
 out=Path(out);audit=json.loads((out/'input_audit.json').read_text())
 with (out/'comparison.csv').open() as f:rows=[r for r in csv.DictReader(f) if r['scheme']=='fine' and r['parameter']=='LU_sin_phi' and r['status']=='valid']
 theory={}
 for model in MODELS[1:]:
  with Path(audit['theory'][model]['path']).open() as f:theory[model]=list(csv.DictReader(f))
 fig,axes=plt.subplots(2,2,figsize=(12,9))
 for xi,ax in enumerate(axes.flat):
  for mi,model in enumerate(MODELS):
   rr=[r for r in rows if r['model']==model and int(r['x_index'])==xi]
   ax.errorbar([float(r['mean_minus_tprime'])+(mi-1)*.012 for r in rr],[float(r['value']) for r in rr],yerr=[float(r['stat']) for r in rr],fmt=MARKERS[mi],color=COLORS[mi],mfc='white',capsize=2,ls='none',label=['Fit: UU = 0','Fit: UU = BM stage 2','Fit: UU = full BM'][mi],zorder=5)
   if model!='zero':
    t=[r for r in theory[model] if np.isclose(float(r['xB']),[.2,.3,.4,.5][xi])]
    ax.plot([float(r['minus_tprime']) for r in t],[float(r['A_LU_sinphi']) for r in t],color=COLORS[mi],ls='--',lw=1.7,label='LU prediction: '+('BM stage 2' if model=='step1_2' else 'full BM'))
  lo,hi=XB[xi];ax.set_title(rf'${lo:.2f}<x_B<{hi:.2f}$');ax.set(xlim=(0,1.26),ylim=(-.2,.3),xlabel=r'$-t^\prime\;(\mathrm{GeV}^2)$',ylabel=r'$F_{LU}^{\sin\phi}/F_{UU}$');ax.set_yticks(np.arange(-.2,.31,.1));ax.axhline(0,c='.6',lw=.8);ax.grid(alpha=.12)
 handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.01),ncol=2,frameon=False,fontsize=10)
 fig.suptitle('π− LU: fixed-UU refits and BM predictions',fontsize=17);fig.subplots_adjust(left=.10,right=.98,top=.91,bottom=.19,wspace=.25,hspace=.35)
 fig.savefig(out/'figures/fine_LU_theory_check.png',dpi=250);plt.close(fig)
