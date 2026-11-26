"""Merge saved fine-bin histograms exactly, without another ROOT scan."""
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def build(source,out,dpi):
 meta=json.loads((source/'kinematic_distributions/kinematic_summary.json').read_text())
 release=json.loads((out/'release_manifest.json').read_text())
 if meta['cutset_id']!=release['identity']['cutset_id']:raise ValueError('Histogram cut-set mismatch')
 names=list(release['datasets']);rows=[]
 with np.load(source/'kinematic_distributions/histograms.npz') as data:
  for sample,prefix in [('RGB','RGB_'),('RGC','RGC_')]:
   for variable,label in [('Q2',r'$Q^2$ (GeV$^2$)'),('W',r'$W$ (GeV)'),('phi',r'$\phi$ (degrees)')]:
    fine=sum(data[name+'__h_'+variable] for name in names if name.startswith(prefix))
    bins=fine.reshape(4,3,2,-1).sum(axis=2)
    if not np.array_equal(bins.sum(axis=(0,1)),fine.sum(axis=0)):raise ValueError('Histogram rebin lost counts')
    edges=np.array(meta['edges'][variable]);x=(edges[:-1]+edges[1:])/2
    fig,axes=plt.subplots(2,2,figsize=(10.5,7.5))
    for xi,ax in enumerate(axes.flat):
     for ti,color in enumerate(['#0072B2','#D55E00','#009E73']):
      h=bins[xi,ti];norm=h.sum();density=h/max(norm,1)
      ax.stairs(density,edges,color=color,label=f'{[.05,.45,.85][ti]:.2f}–{[.45,.85,1.25][ti]:.2f}')
      for j,n in enumerate(h):rows.append(dict(sample=sample,x_index=xi,tprime_index=ti,variable=variable,low=edges[j],high=edges[j+1],combinations=int(n),fraction_in_displayed_range=float(density[j])))
     ax.set(xlabel=label,ylabel='Fraction per display bin',title=f'{[.10,.25,.35,.45][xi]:.2f} < xB < {[.25,.35,.45,.60][xi]:.2f}');ax.grid(alpha=.12)
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,title=r'$-t^\prime$ intervals (GeV$^2$)',frameon=False)
    fig.suptitle(sample+' · three-bin kinematic distributions');fig.tight_layout(rect=[0,.12,1,.95]);fig.savefig(out/'figures'/f'{sample}_{variable}_three_bin_2x2.png',dpi=dpi);plt.close(fig)
 with (out/'three_bin_histograms.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

 # Add the new binned spectra to the existing release gallery.
 with (out/'index.html').open('a') as gallery:
  for sample in ['RGB','RGC']:
   for variable in ['Q2','W','phi']:
    name=f'{sample}_{variable}_three_bin_2x2.png'
    gallery.write(f'<p><a href="figures/{name}">{sample} {variable}, three-bin kinematics</a></p>')
