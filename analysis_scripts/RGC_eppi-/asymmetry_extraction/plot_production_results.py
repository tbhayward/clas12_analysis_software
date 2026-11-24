#!/usr/bin/env python3
"""Statistical-only pi− plots in the pi+ final-results layout."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from inputs import PERIODS
from likelihood import NAMES
HERE=Path(__file__).resolve().parent
PANELS=[0,1,3,4,5,6,7]
LABELS=[r'$F_{UU}^{\cos\phi}/F_{UU}$',r'$F_{UU}^{\cos2\phi}/F_{UU}$',r'$F_{LU}^{\sin\phi}/F_{UU}$',r'$F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$',r'$F_{UL,\mathrm{lab}}^{\sin2\phi}/F_{UU}$',r'$F_{LL,\mathrm{lab}}/F_{UU}$',r'$F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$']
XB=[.10,.25,.35,.45,.60]
TP=np.array([.05,.25,.45,.65,.85,1.05,1.25])
COLORS=['#0072B2','#56B4E9','#009E73','#D55E00','#CC79A7','#E69F00']
MARKERS=['o','s','^','D','v','P']
PERIOD_LABELS=['RGB Sp19','RGB Fa19','RGB Sp20','RGC Su22','RGC Fa22','RGC Sp23']


def axes_style(ax,j):
    limit=.4 if j in [2,3,4] else 1.
    ax.set(xlim=(0,1.30),ylim=(-limit,limit),xlabel=r'$-t^\prime$ (GeV$^2$)',ylabel=LABELS[j])
    ax.set_yticks(np.linspace(-limit,limit,9 if limit==.4 else 11))
    ax.set_xticks(np.arange(0,1.31,.2))
    ax.axhline(0,color='0.45',lw=.8,zorder=1)
    ax.grid(alpha=.18,zorder=0)
    return limit


def load_means(registry_path):
    sys.path.insert(0,str(HERE.parent/'channel_selection'))
    from selection import load_selected_inputs
    import uproot
    registry=load_selected_inputs(registry_path)
    means={};totals=np.zeros(24);numbers=np.zeros(24)
    for scope in PERIODS:
        with uproot.open(registry['files'][scope]['path']) as root:
            a=root[registry['tree']].arrays(['analysis_bin','tprime'],library='np')
        b=a['analysis_bin']-1
        sums=np.bincount(b,weights=-a['tprime'],minlength=24)
        counts=np.bincount(b,minlength=24)
        means[scope]=np.divide(sums,counts,out=np.full(24,np.nan),where=counts>0)
        totals+=sums;numbers+=counts
    means['combined']=np.divide(totals,numbers,out=np.full(24,np.nan),where=numbers>0)
    return means,registry['cutset_id']


def plot(output,registry_path=None):
    output=Path(output);folder=output/'figures/production';folder.mkdir(parents=True,exist_ok=True)
    results=[json.loads(p.read_text()) for p in (output/'fits').glob('*.json')]
    lookup={(r['scope'],r['bin']):r for r in results}
    if registry_path:
        means,cutset=load_means(registry_path)
        audit=json.loads((output/'input_audit.json').read_text())
        if audit['identity']['cutset_id']!=cutset:raise ValueError('Plot inputs have a different cut set')
        convention='Mean selected-event minus-tprime, independently for each scope/bin'
    else:
        means={scope:np.tile((TP[:-1]+TP[1:])/2,4) for scope in ['combined']+PERIODS}
        convention='Bin centers (selected-input registry not supplied)'
    omissions=[];clipped=[];positions=[]
    def points(ax,scope,ix,j,color='black',marker='o',label=None,offset=0):
        limit=.4 if j in [2,3,4] else 1.
        xx=[];yy=[];ee=[]
        for it in range(6):
            b=ix*6+it+1;r=lookup.get((scope,b))
            if not r or r['status']!='valid' or j>=len(r.get('values',[])):
                omissions.append(dict(scope=scope,bin=b,parameter=NAMES[j],reason='missing' if r is None else r['status']))
                continue
            x=float(means[scope][b-1]);y=r['values'][j];error=r['errors'][j]
            positions.append(dict(scope=scope,bin=b,parameter=NAMES[j],x=x,horizontal_display_offset=offset))
            xx.append(x+offset);yy.append(y);ee.append(error)
            if y-error < -limit or y+error > limit:
                clipped.append(dict(scope=scope,bin=b,parameter=NAMES[j],value=y,stat=error,axis_low=-limit,axis_high=limit))
                # Edge triangles disclose any part of the error interval beyond the frame.
                if y+error>limit:ax.plot(x+offset,limit*.97,marker='^',color=color,ms=4,linestyle='none')
                if y-error < -limit:ax.plot(x+offset,-limit*.97,marker='v',color=color,ms=4,linestyle='none')
        ax.errorbar(xx,yy,yerr=ee,fmt=marker,color=color,ecolor=color,ms=4,capsize=3,lw=1.1,linestyle='none',label=label,zorder=3)
    for ix in range(4):
        fig,axes=plt.subplots(3,3,figsize=(12,10));axes=axes.ravel()
        for j,panel in enumerate(PANELS):
            axes_style(axes[panel],j);points(axes[panel],'combined',ix,j)
        axes[2].axis('off');axes[8].axis('off')
        missing=[ix*6+it+1 for it in range(6) if lookup.get(('combined',ix*6+it+1),{}).get('status')!='valid']
        if missing:axes[8].text(.05,.65,'Fits requiring validation omitted:\nbins '+', '.join(map(str,missing)),transform=axes[8].transAxes,fontsize=10)
        axes[8].text(.05,.3,'Statistical uncertainties only.\nEdge triangles: interval extends\nbeyond the fixed axis range.',transform=axes[8].transAxes,fontsize=9)
        fig.suptitle(rf'${XB[ix]:.2f}\leq x_B<{XB[ix+1]:.2f}$')
        fig.tight_layout(rect=(0,0,1,.96));fig.savefig(folder/f'final_results_xB_bin_{ix+1}.png',dpi=200);plt.close(fig)
    fig,axes=plt.subplots(3,3,figsize=(12,10));axes=axes.ravel()
    for j,panel in enumerate(PANELS):
        axes_style(axes[panel],j)
        for ix,(color,marker) in enumerate(zip(['#0072B2','#D55E00','#009E73','#CC79A7'],['o','s','^','D'])):
            points(axes[panel],'combined',ix,j,color,marker,rf'${XB[ix]:.2f}\leq x_B<{XB[ix+1]:.2f}$')
    axes[2].axis('off');axes[8].axis('off')
    handles,labels=axes[0].get_legend_handles_labels();axes[8].legend(handles,labels,loc='center',frameon=False)
    fig.tight_layout();fig.savefig(folder/'final_results_all_xB.png',dpi=200);plt.close(fig)
    for mode,indices,scopes in [('UU_LU',[0,1,2],PERIODS),('UL_LL',[3,4,5,6],PERIODS[3:])]:
        for ix in range(4):
            fig,axes=plt.subplots(1,3,figsize=(14,4.5)) if mode=='UU_LU' else plt.subplots(2,2,figsize=(10,8))
            axes=np.ravel(axes)
            for ax,j in zip(axes,indices):
                axes_style(ax,j)
                for index,scope in enumerate(scopes):
                    ip=PERIODS.index(scope)
                    points(ax,scope,ix,j,COLORS[ip],MARKERS[ip],PERIOD_LABELS[ip],(index-(len(scopes)-1)/2)*.012)
            handles,labels=axes[0].get_legend_handles_labels()
            fig.legend(handles,labels,loc='upper center',ncol=len(scopes),frameon=False,fontsize=9,bbox_to_anchor=(.5,.94))
            fig.suptitle(rf'${XB[ix]:.2f}\leq x_B<{XB[ix+1]:.2f}$')
            fig.tight_layout(rect=(0,0,1,.85 if mode=='UU_LU' else .89))
            fig.savefig(folder/f'period_{mode}_xB_bin_{ix+1}.png',dpi=200);plt.close(fig)
    for name,rows in [('plot_omissions.csv',omissions),('plot_axis_clipping.csv',clipped),('plot_positions.csv',positions)]:
        if rows:
            unique={json.dumps(r,sort_keys=True):r for r in rows};rows=list(unique.values())
            with (folder/name).open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    (folder/'plot_conventions.json').write_text(json.dumps(dict(style_reference='RGC_enpi+/asymmetry_extraction/plot_final_results.py',statistical_only=True,UU_LL_limits=[-1,1],LU_UL_limits=[-.4,.4],positions=convention),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=HERE/'output');p.add_argument('--selected-inputs',type=Path);a=p.parse_args();plot(a.output,a.selected_inputs)
