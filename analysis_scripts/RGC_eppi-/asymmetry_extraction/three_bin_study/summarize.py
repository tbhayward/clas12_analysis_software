"""Review plots and paired-refit sensitivity tables. Shift scales are not significances."""
import csv
import json
import shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES = ['LU_sin_phi','UL_lab_sin_phi','UL_lab_sin_2phi','LL_lab','LL_lab_cos_phi']
LABELS = [r'$F_{LU}^{\sin\phi}/F_{UU}$',r'$F_{UL,lab}^{\sin\phi}/F_{UU}$',r'$F_{UL,lab}^{\sin2\phi}/F_{UU}$',r'$F_{LL,lab}/F_{UU}$',r'$F_{LL,lab}^{\cos\phi}/F_{UU}$']
MODELS = ['free','step1_2','step1_2_3']
COLORS = ['#222222','#d62728','#1f77b4']
MARKERS = ['o','s','^']
LEGENDS = ['LL cosφ free','LL cosφ fixed: BM stage 2','LL cosφ fixed: BM stage 3']


def products(out, results, provenance, dfrows, fine_folder):
    out = Path(out)
    figures = out/'figures'
    figures.mkdir(exist_ok=True)
    lookup = {(r['model'],r['bin']):r for r in results}
    rows = []
    summary = {}
    for model in MODELS:
        subset = [r for r in results if r['model']==model]
        summary[model] = dict(valid=sum(r['status']=='valid' for r in subset), total=len(subset), warnings=[dict(bin=r['bin'],warnings=r.get('warnings',[])) for r in subset if r.get('warnings')])
    for name in NAMES:
        summary[name] = {}
        for model in MODELS:
            shifts, ratios = [], []
            for b in range(1,13):
                r = lookup[model,b]
                if r['status'] != 'valid':
                    continue
                j = r['names'].index(name)
                baseline = lookup['free',b]
                delta = scale = ratio = None
                if baseline['status']=='valid':
                    bj = baseline['names'].index(name)
                    delta = r['values'][j]-baseline['values'][bj]
                    scale = delta/baseline['errors'][bj]
                    ratio = r['errors'][j]/baseline['errors'][bj]
                    shifts.append(abs(scale))
                    ratios.append(ratio)
                rows.append(dict(model=model,bin=b,parameter=name,parameter_status='fixed' if name in r['fixed_parameters'] else 'fitted',value=r['values'][j],stat=r['errors'][j],shift_from_free=delta,shift_in_free_sigma=scale,error_ratio_to_free=ratio,events=r['events'],fixed_LLcos=r['theory_input']['fixed_LLcos']))
            if shifts:
                summary[name][model] = dict(median_abs_shift_in_free_sigma=float(np.median(shifts)),max_abs_shift_in_free_sigma=float(max(shifts)),median_error_ratio=float(np.median(ratios)))
        for shifted in ([False] if name=='LL_lab_cos_phi' else [False,True]):
            fig, axes = plt.subplots(2,2,figsize=(11,8),sharex=True,sharey=True)
            limits = []
            for ix,ax in enumerate(axes.flat):
                for m,model in enumerate(MODELS):
                    if shifted and model=='free':
                        continue
                    for b in range(ix*3+1,ix*3+4):
                        r = lookup[model,b]
                        if r['status']!='valid':
                            continue
                        j = r['names'].index(name)
                        x = provenance['means'][b-1]+(m-1)*.012
                        if shifted:
                            base = lookup['free',b]
                            if base['status']!='valid':continue
                            y = (r['values'][j]-base['values'][j])/base['errors'][j]
                            ax.plot(x,y,MARKERS[m],color=COLORS[m],label=LEGENDS[m] if b==ix*3+1 else None)
                            limits.append(abs(y)*1.1)
                        else:
                            y,e = r['values'][j],r['errors'][j]
                            ax.errorbar(x,np.clip(y,-.98,.98),yerr=e if -1<=y<=1 else None,fmt=MARKERS[m] if -1<=y<=1 else ('^' if y>0 else 'v'),color=COLORS[m],mfc='white',capsize=3,label=LEGENDS[m] if b==ix*3+1 else None)
                            if abs(y)>1:
                                ax.annotate(f'{y:.2f}',(x,np.sign(y)*.92),fontsize=8,ha='center')
                            limits.append(abs(y)+e)
                ax.axhline(0,color='.6',lw=1)
                ax.grid(alpha=.15)
                ax.set_xlim(0,1.3)
                ax.set_title(f'{[.10,.25,.35,.45][ix]:.2f} < xB < {[.25,.35,.45,.60][ix]:.2f}')
                ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
                ax.set_ylabel('Shift / free-fit statistical error' if shifted else LABELS[NAMES.index(name)])
            limit = max(.15, max(limits,default=.2)*1.05)
            if not shifted:limit=min(1.,limit)
            for ax in axes.flat:ax.set_ylim(-limit,limit)
            handles,labels = axes.flat[0].get_legend_handles_labels()
            fig.legend(handles,labels,loc='lower center',ncol=1,frameon=False)
            fig.suptitle('π− three-bin LL cosφ constraint study: '+LABELS[NAMES.index(name)])
            fig.tight_layout(rect=[0,.13,1,.95])
            fig.savefig(figures/(name+('_shifts' if shifted else '_comparison')+'.png'),dpi=180)
            plt.close(fig)
    with (out/'comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'fit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(9,5))
    for ix in range(4):
        rr=sorted([r for r in dfrows if (r['bin']-1)//3==ix],key=lambda r:r['bin'])
        ax.errorbar([provenance['means'][r['bin']-1] for r in rr],[r['df'] for r in rr],yerr=[r['stat_std'] for r in rr],fmt='o',capsize=3,label=f'xB interval {ix+1}')
    ax.set(xlabel=r"$-t^\prime$ (GeV$^2$)",ylabel='Shared ND3 dilution factor',ylim=(0,1))
    ax.legend();fig.tight_layout();fig.savefig(figures/'three_bin_dilution.png',dpi=180);plt.close(fig)
    for p in [fine_folder/'figures/nominal_LU_sin_phi.png',fine_folder/'ratios.csv']:
        if p.exists():shutil.copy2(p,out/('fine_LU_reference.png' if p.suffix=='.png' else 'fine_reference_ratios.csv'))
    text='''Three-bin target-sector study. Original six-bin LU remains the preferred fine result.
The three scenarios use identical selected combinations and recomputed 12-bin dilution.
Both UU ratios remain fixed zero. LLcos is free or fixed to supplied BM stage 2/3.
Predictions are used directly without photon-to-lab conversion, at user request.
Fixed inputs have no theory uncertainty. Smaller errors are conditional on the input.
Shift / free-fit error is a sensitivity scale, not a significance of correlated differences.
Dilution covariance includes the refit response of free amplitudes to a fixed physical LLcos.
Plots cap ratio axes at +/-1 and label offscale central values. See JSON/CSV for full errors.
'''
    (out/'assessment.txt').write_text(text)
    links=''.join(f'<p><a href="figures/{p.name}">{p.stem}</a></p>' for p in sorted(figures.glob('*.png')))
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Three-bin LLcos study</title><h1>Three-bin LLcos study</h1><pre>'+text+'</pre><a href="comparison.csv">Comparison table</a>'+links)
