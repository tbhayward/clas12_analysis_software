#!/usr/bin/env python3
"""pi+ note's eight kinematic plot products, adapted to nominal e' p pi- data."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import csv
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'channel_selection'))
from selection import load_selected_inputs
from kinematics import XB,TP,bin_indices

PERIODS=['RGB_spring2019_D2','RGB_fall2019_D2','RGB_spring2020_D2','RGC_summer2022_ND3','RGC_fall2022_ND3','RGC_spring2023_ND3']
EDGES={'x':np.linspace(.1,.6,51),'Q2':np.linspace(1,12,55),'tp':np.linspace(.05,1.25,49),'phi':np.linspace(0,360,49),'W':np.linspace(2,4.5,50),'y':np.linspace(0,.8,49),'pion_p':np.arange(.5,10.5001,.25),'proton_p':np.linspace(0,5,101),'proton_theta':np.linspace(0,180,91)}
LABELS={'x':r'$x_B$','Q2':r'$Q^2$ (GeV$^2$)','tp':r'$-t\prime$ (GeV$^2$)','phi':r'$\phi$ (deg)','W':r'$W$ (GeV)','y':r'$y$','pion_p':r'$p_{\pi^-}$ (GeV)'}
LABELS.update(proton_p=r'$p_p$ (GeV)',proton_theta=r'$\theta_p$ (deg)')
CORE=[('x','Q2'),('x','tp'),('Q2','tp')]
PHI=[('phi','Q2'),('phi','x'),('phi','tp')]
BRANCHES=['x','Q2','tprime','phi2','W','y','p2_p','analysis_bin','p1_p','p1_theta','detector1']

def read_period(job):
    import uproot
    name,path,tree_name,cutset_id,step=job
    h={f'h_{k}':np.zeros((24,len(v)-1),dtype=np.int64) for k,v in EDGES.items()}
    h.update({f'c_{x}_{y}':np.zeros((len(EDGES[x])-1,len(EDGES[y])-1),dtype=np.int64) for x,y in CORE+PHI})
    for detector in ['FD','CD']:
        h['proton_p_'+detector]=np.zeros((len(EDGES['proton_p'])-1,),dtype=np.int64)
        h['proton_2d_'+detector]=np.zeros((len(EDGES['proton_p'])-1,len(EDGES['proton_theta'])-1),dtype=np.int64)
    h['qw']=np.zeros((4,len(EDGES['Q2'])-1,len(EDGES['W'])-1),dtype=np.int64)
    h['tails']=np.zeros((len(EDGES),2),dtype=np.int64)
    h['moments']=np.zeros((24,6,3))  # count, sum, sum2 for x,Q2,W,y,tp,pion_p
    read=0
    with uproot.open(path) as f:
        if str(f['selection_cutset_id'])!=cutset_id:raise ValueError('Mixed cut sets: '+name)
        meta=json.loads(str(f['selection_manifest']))
        if meta['dataset']!=name or meta['partial']:raise ValueError('Wrong or partial sample: '+name)
        tree=f[tree_name]
        for a in tree.iterate(BRANCHES,library='np',step_size=step):
            bins=bin_indices(a['x'],-a['tprime'])
            if np.any(bins<0) or not np.array_equal(bins+1,a['analysis_bin']):raise ValueError('Invalid analysis-bin annotations: '+name)
            values=dict(x=a['x'],Q2=a['Q2'],tp=-a['tprime'],phi=np.degrees(np.mod(a['phi2'],2*np.pi)),W=a['W'],y=a['y'],pion_p=a['p2_p'],proton_p=a['p1_p'],proton_theta=np.degrees(a['p1_theta']))
            for j,(key,edges) in enumerate(EDGES.items()):
                v=values[key]
                if not np.all(np.isfinite(v)):raise ValueError('Nonfinite '+key+' in '+name)
                index=np.searchsorted(edges,v,side='right')-1
                keep=(index>=0)&(index<len(edges)-1)
                h['h_'+key]+=np.bincount(bins[keep]*(len(edges)-1)+index[keep],minlength=24*(len(edges)-1)).reshape(24,-1)
                h['tails'][j]+=[np.count_nonzero(v<edges[0]),np.count_nonzero(v>=edges[-1])]
            if not np.all(np.isin(a['detector1'],[1,2])):raise ValueError('Invalid proton detector')
            for detector,code in [('FD',1),('CD',2)]:
                selected=a['detector1']==code
                h['proton_p_'+detector]+=np.histogram(values['proton_p'][selected],bins=EDGES['proton_p'])[0]
                h['proton_2d_'+detector]+=np.histogram2d(values['proton_p'][selected],values['proton_theta'][selected],bins=(EDGES['proton_p'],EDGES['proton_theta']))[0].astype(np.int64)
            for x,y in CORE+PHI:h[f'c_{x}_{y}']+=np.histogram2d(values[x],values[y],bins=(EDGES[x],EDGES[y]))[0].astype(np.int64)
            for xi in range(4):
                keep=bins//6==xi
                h['qw'][xi]+=np.histogram2d(values['Q2'][keep],values['W'][keep],bins=(EDGES['Q2'],EDGES['W']))[0].astype(np.int64)
            for j,key in enumerate(['x','Q2','W','y','tp','pion_p']):
                v=values[key]
                for k,weights in enumerate([None,v,v*v]):h['moments'][:,j,k]+=np.bincount(bins,weights=weights,minlength=24)
            read+=len(bins)
            print(f'{name}: {read:,}/{tree.num_entries:,} selected combinations',flush=True)
    return name,h

def shape(ax,counts,edges,label=None,color=None):
    total=counts.sum();ax.stairs(counts/max(1,total),edges,label=label,color=color,lw=1.3)
    ax.set_ylim(bottom=0)

def heatmap(fig,ax,z,x,y):
    from matplotlib.colors import LogNorm
    if z.max()>0:
        image=ax.pcolormesh(EDGES[x],EDGES[y],np.ma.masked_less(z.T,1),cmap='turbo',norm=LogNorm(1,max(2,int(z.max()))))
        fig.colorbar(image,ax=ax,label='Combinations')
    else:ax.text(.5,.5,'No selected entries',transform=ax.transAxes,ha='center')
    ax.set(xlabel=LABELS[x],ylabel=LABELS[y],xlim=(EDGES[x][0],EDGES[x][-1]),ylim=(EDGES[y][0],EDGES[y][-1]))

def plot_sample(name,h,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=True)
    def save(fig,file):
        fig.suptitle(name,fontsize=13);fig.savefig(out/file,dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(11.5,7),constrained_layout=True)
    for ax,key in zip(axes.flat,['x','Q2','tp','phi','W','y']):
        shape(ax,h['h_'+key].sum(axis=0),EDGES[key],color='black');ax.set(xlabel=LABELS[key],ylabel='Fraction of combinations / bin')
    save(fig,'01_kinematic_1d.png')
    fig,axes=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    for ax,(x,y) in zip(axes,CORE):heatmap(fig,ax,h[f'c_{x}_{y}'],x,y)
    save(fig,'02_core_correlations.png')
    fig,axes=plt.subplots(3,1,figsize=(9,10),constrained_layout=True)
    for ax,(x,y) in zip(axes,PHI):heatmap(fig,ax,h[f'c_{x}_{y}'],x,y)
    save(fig,'03_phi_correlations.png')
    fig,axes=plt.subplots(2,2,figsize=(10.5,7.5),constrained_layout=True)
    for xi,ax in enumerate(axes.flat):
        shape(ax,h['h_phi'][xi*6:(xi+1)*6].sum(axis=0),EDGES['phi'],color='black')
        ax.set(title=f'{XB[xi]:.2f} ≤ xB < {XB[xi+1]:.2f}',xlabel=LABELS['phi'],ylabel='Fraction of combinations / bin')
    save(fig,'04_phi_projections.png')
    for number,key in [('05','Q2'),('06','W')]:
        fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
        for ti,ax in enumerate(axes.flat):
            for xi in range(4):shape(ax,h['h_'+key][xi*6+ti],EDGES[key],label=f'{XB[xi]:.2f}–{XB[xi+1]:.2f}')
            ax.set(title=f'{TP[ti]:.2f} ≤ −t′ < {TP[ti+1]:.2f}',xlabel=LABELS[key],ylabel='Fraction / bin')
        axes.flat[0].legend(title='xB',fontsize=8)
        save(fig,f'{number}_{key}_by_analysis_bin.png')
    fig,axes=plt.subplots(2,2,figsize=(10.5,8),constrained_layout=True)
    for xi,ax in enumerate(axes.flat):
        heatmap(fig,ax,h['qw'][xi],'Q2','W');ax.set_title(f'{XB[xi]:.2f} ≤ xB < {XB[xi+1]:.2f}')
    save(fig,'07_Q2_W_by_xB.png')
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for xi,ax in enumerate(axes.flat):
        for ti in range(6):shape(ax,h['h_pion_p'][xi*6+ti],EDGES['pion_p'],label=f'{TP[ti]:.2f}–{TP[ti+1]:.2f}')
        ax.set(title=f'{XB[xi]:.2f} ≤ xB < {XB[xi+1]:.2f}',xlabel=LABELS['pion_p'],ylabel='Fraction / bin')
    axes.flat[0].legend(title='−t′ (GeV²)',fontsize=8)
    save(fig,'08_pion_momentum_by_analysis_bin.png')
    fig,axes=plt.subplots(1,3,figsize=(15,4),constrained_layout=True)
    total=h['proton_p_FD']+h['proton_p_CD'];normalization=max(1,total.sum())
    for label,counts,color in [('All',total,'black'),('FD',h['proton_p_FD'],'tab:blue'),('CD',h['proton_p_CD'],'tab:orange')]:
        axes[0].stairs(counts/normalization,EDGES['proton_p'],color=color,label=f'{label} (N={int(counts.sum()):,})',lw=1.3)
    axes[0].set(xlabel=LABELS['proton_p'],ylabel='Fraction of all protons / bin',ylim=(0,None));axes[0].legend(fontsize=8)
    for ax,detector in zip(axes[1:],['FD','CD']):
        heatmap(fig,ax,h['proton_2d_'+detector],'proton_p','proton_theta');ax.set_title(detector+' protons')
    save(fig,'09_proton_momentum_FD_CD.png')
    print('PLOTTED '+name,flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,default=Path(__file__).resolve().parents[1]/'channel_selection/output/selected_inputs.json')
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'output')
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--step-size',default='100 MB')
    p.add_argument('--per-period',action='store_true',help='Also produce the eight plots separately for all six periods')
    p.add_argument('--plots-only',action='store_true')
    args=p.parse_args()
    if args.workers<1:p.error('workers must be positive')
    args.output.mkdir(parents=True,exist_ok=True)
    data={}
    if args.plots_only:
        with np.load(args.output/'histograms.npz') as h:
            for name in PERIODS:data[name]={k.split('__',1)[1]:h[k] for k in h.files if k.startswith(name+'__')}
        meta=json.loads((args.output/'kinematic_summary.json').read_text())
    else:
        registry=load_selected_inputs(args.inputs,verify_files=False)
        jobs=[(n,registry['files'][n]['path'],registry['tree'],registry['cutset_id'],args.step_size) for n in PERIODS]
        with ProcessPoolExecutor(max_workers=min(args.workers,6)) as pool:
            for future in as_completed([pool.submit(read_period,j) for j in jobs]):
                name,h=future.result();data[name]=h
        np.savez_compressed(args.output/'histograms.npz',**{n+'__'+k:v for n,h in data.items() for k,v in h.items()})
        meta=dict(cutset_id=registry['cutset_id'],inputs={n:registry['files'][n] for n in PERIODS},azimuth='phi2 converted from radians to degrees',pion_momentum='p2_p',counting='Unweighted candidate combinations. Shape-normalized curves are not charge-normalized yields.',reference='RGC_enpi+/kinematic_distributions/plot_kinematic_distributions_review.py and pi+ note Sec. 3.7, Figs. 31–37',edges={k:v.tolist() for k,v in EDGES.items()},tails={n:{k:dict(underflow=int(h['tails'][j,0]),overflow=int(h['tails'][j,1])) for j,k in enumerate(EDGES)} for n,h in data.items()})
        (args.output/'kinematic_summary.json').write_text(json.dumps(meta,indent=2))
    for group in ['RGB','RGC']:
        members=[n for n in PERIODS if n.startswith(group)]
        data[group]={k:sum(data[n][k] for n in members) for k in data[members[0]]}
        plot_sample(group,data[group],args.output/group)
    if args.per_period:
        for name in PERIODS:plot_sample(name,data[name],args.output/name)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(11.5,7),constrained_layout=True)
    for ax,key in zip(axes.flat,['x','Q2','tp','phi','W','y']):
        for group,color in [('RGB','tab:blue'),('RGC','tab:red')]:shape(ax,data[group]['h_'+key].sum(axis=0),EDGES[key],label=group,color=color)
        ax.set(xlabel=LABELS[key],ylabel='Fraction of combinations / bin')
    axes.flat[0].legend();fig.suptitle('RGB D₂ and RGC ND₃: selected kinematic shapes')
    fig.savefig(args.output/'RGB_RGC_shape_comparison.png',dpi=160);plt.close(fig)
    with (args.output/'pion_momentum_histograms_by_analysis_bin.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['dataset','analysis_bin','xB_low','xB_high','minus_tprime_low','minus_tprime_high','p_low_GeV','p_high_GeV','combinations'])
        for name in PERIODS:
            for b in range(24):
                for j,n in enumerate(data[name]['h_pion_p'][b]):w.writerow([name,b+1,XB[b//6],XB[b//6+1],TP[b%6],TP[b%6+1],EDGES['pion_p'][j],EDGES['pion_p'][j+1],int(n)])
    with (args.output/'bin_kinematic_means.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['sample','analysis_bin','variable','combinations','mean','std'])
        for name,h in data.items():
            for b in range(24):
                for j,key in enumerate(['x','Q2','W','y','tp','pion_p']):
                    n,s,ss=h['moments'][b,j];mean=s/n if n else 0
                    w.writerow([name,b+1,key,int(n),mean if n else '',np.sqrt(max(0,ss/n-mean*mean)) if n else ''])
    names=['RGB','RGC']+(PERIODS if args.per_period else [])
    links=['<!doctype html><meta charset="utf-8"><h1>Selected kinematic distributions</h1><a href="RGB_RGC_shape_comparison.png">RGB/RGC shape comparison</a>']
    for name in names:
        links.append('<h2>'+name+'</h2>');links.extend(f'<a href="{name}/{f.name}">{f.stem}</a><br>' for f in sorted((args.output/name).glob('*.png')))
    (args.output/'index.html').write_text(''.join(links))
    print('Browse '+str(args.output/'index.html'),flush=True)

if __name__=='__main__':main()
