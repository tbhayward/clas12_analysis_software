#!/usr/bin/env python3
"""Pass-2 proton phase-space diagnostics. Run from dvcs_cross_section root.
Requires numpy, matplotlib, uproot. Five period jobs, max eight workers.
Generated MC is plotted at generator level; no event-level truth matching is assumed.
Reconstructed preselection approximates global preselection (not optimized exclusivity).
All yield ratios are explicitly shape-normalized, not absolute efficiencies.
"""
import argparse
import concurrent.futures as cf
import csv
import os
from pathlib import Path
import re
import numpy as np
os.environ.setdefault('MPLBACKEND','Agg')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot

PERIODS=['Sp18_inb','Sp18_out','Fa18_inb','Fa18_out','Sp19_inb']
BASE=Path('/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2')
ENERGY={'Sp18_inb':'10594','Sp18_out':'10594','Fa18_inb':'10604','Fa18_out':'10604','Sp19_inb':'10200'}
CURRENT={'Sp18_inb':'50','Sp18_out':'45','Fa18_inb':'50','Fa18_out':'50','Sp19_inb':'50'}
PBINS=np.linspace(0.15,2.5,95); TBINS=np.linspace(30,70,81)
CUTS=[(.30,64.23),(.35,62),(.40,60),(.45,60),(.50,58),(.50,50)]

def path_for(period,kind):
    if kind=='data': return BASE/'data/dvcs'/f'rga_{period.lower()}_epgamma.root'
    stem='gen' if kind=='gen' else 'rec'
    return BASE/'mc/dvcsgen'/f'{stem}_dvcsgen_rga_{period.lower()}_{CURRENT[period]}nA_{ENERGY[period]}MeV.root'

def angle(a):
    a=np.asarray(a,dtype=float)
    return np.degrees(a) if np.nanpercentile(np.abs(a),95)<3.5 else a

def collect(period,kind,step):
    file=path_for(period,kind)
    if not file.exists(): raise FileNotFoundError(str(file))
    with uproot.open(file) as f:
        tree=f['PhysicsEvents']; keys=set(tree.keys())
        required=['p1_p','p1_theta']
        if any(x not in keys for x in required): raise RuntimeError(f'{file}: missing {set(required)-keys}; available {sorted(keys)[:35]}')
        optional=['detector1','t1','open_angle_ep2','pTmiss','e_theta','p1_phi','e_phi','p2_phi']
        branches=required+[x for x in optional if x in keys]
        h2=np.zeros((len(PBINS)-1,len(TBINS)-1)); hlow=np.zeros_like(h2); hmid=np.zeros_like(h2); hhigh=np.zeros_like(h2)
        counts=np.zeros(len(CUTS)); total=0; lowtotal=0; cdtotal=0
        for a in tree.iterate(branches,library='np',step_size=step):
            p=np.asarray(a['p1_p'],float); th=angle(a['p1_theta']); mask=np.isfinite(p)&np.isfinite(th)
            if kind!='gen':
                if 'detector1' in a: mask &= a['detector1']==2
                if 't1' in a: mask &= (-a['t1']<1.)
                if 'open_angle_ep2' in a: mask &= angle(a['open_angle_ep2'])>5
            # Generator sample: no detector category exists; this is full generated CD-angle phase space.
            else: mask &= (th>=30)&(th<70)
            p=p[mask]; th=th[mask]
            t=np.asarray(a['t1'],float)[mask] if 't1' in a else None
            h2+=np.histogram2d(p,th,bins=(PBINS,TBINS))[0]
            if t is not None:
                regions=[(-t<.2),((-t>=.2)&(-t<.4)),(-t>=.4)]
                for h,sel in zip((hlow,hmid,hhigh),regions): h+=np.histogram2d(p[sel],th[sel],bins=(PBINS,TBINS))[0]
                lowtotal+=int(np.sum(regions[0]))
            total+=len(p)
            for i,(pc,tc) in enumerate(CUTS): counts[i]+=np.count_nonzero((p>pc)&(th<tc))
        return dict(kind=kind,period=period,h=h2,low=hlow,mid=hmid,high=hhigh,counts=counts,total=total,lowtotal=lowtotal)

def plot1d(results,period,out):
    fig,axs=plt.subplots(2,2,figsize=(13,9)); colors={'data':'black','rec':'tab:blue','gen':'tab:orange'}
    for kind,r in results.items():
        h=r['h']; color=colors[kind]
        for ax,v,bins in ((axs[0,0],h.sum(axis=1),PBINS),(axs[0,1],h.sum(axis=0),TBINS)):
            # Unit-area shape comparison; no claim of absolute data/MC normalization.
            y=v/max(v.sum(),1); ax.stairs(y,bins,label=kind,color=color,linewidth=1.5)
        for ax,key in ((axs[1,0],'low'),(axs[1,1],'high')):
            v=r[key].sum(axis=0); ax.stairs(v/max(v.sum(),1),TBINS,label=kind,color=color,linewidth=1.5)
    for ax in axs.flat: ax.grid(alpha=.2); ax.legend(fontsize=8)
    axs[0,0].set(xlabel='Proton momentum (GeV)',ylabel='Unit-area fraction',title='Momentum: all selected events')
    axs[0,1].set(xlabel='Proton angle (degrees)',ylabel='Unit-area fraction',title='Angle: all selected events')
    axs[1,0].set(xlabel='Proton angle (degrees)',ylabel='Unit-area fraction',title='Low |t| (<0.2 GeV²)')
    axs[1,1].set(xlabel='Proton angle (degrees)',ylabel='Unit-area fraction',title='High |t| (≥0.4 GeV²)')
    for ax in axs[:,1]:
        for x in (50,58,64.23): ax.axvline(x,color='gray',linestyle='--',alpha=.5)
    axs[0,0].axvline(.5,color='gray',linestyle='--')
    fig.suptitle(f'{period}: generated MC / reconstructed MC / reconstructed data (shape only)')
    fig.tight_layout();fig.savefig(out/'distributions.png',dpi=155);plt.close(fig)

def maps(results,period,out):
    fig,axs=plt.subplots(2,3,figsize=(16,10),layout='constrained')
    for j,kind in enumerate(('gen','rec','data')):
        h=results[kind]['h']; im=axs[0,j].pcolormesh(PBINS,TBINS,np.log10(h.T+1),shading='auto',cmap='viridis');fig.colorbar(im,ax=axs[0,j]);axs[0,j].set_title(f'{kind}: log10(count+1)')
    rec=results['rec']['h'];dat=results['data']['h'];gen=results['gen']['h']
    # Only relative shape ratios. Not detector efficiency or absolute data/MC ratio.
    for j,(num,den,title) in enumerate(((rec,gen,'rec/gen shape ratio (NOT efficiency)'),(dat,rec,'data/recMC shape ratio'),(results['data']['low'],results['rec']['low'],'low-|t| data/recMC shape ratio'))):
        norm=(den.sum()/num.sum()) if num.sum()>0 else np.nan
        val=np.divide(num*norm,den,out=np.full_like(num,np.nan),where=den>=15)
        im=axs[1,j].pcolormesh(PBINS,TBINS,val.T,shading='auto',vmin=0,vmax=2,cmap='coolwarm');fig.colorbar(im,ax=axs[1,j]);axs[1,j].set_title(title+' (unit-area normalized)')
    for ax in axs.flat:
        ax.axvline(.5,color='white',ls='--',lw=1);ax.axhline(58,color='white',ls='--',lw=1)
        ax.set(xlabel='Proton momentum (GeV)',ylabel='Proton angle (degrees)')
    fig.suptitle(period);fig.savefig(out/'maps.png',dpi=145);plt.close(fig)

def worker(args):
    period,outdir,step=args; out=Path(outdir)/period;out.mkdir(parents=True,exist_ok=True)
    results={k:collect(period,k,step) for k in ('data','rec','gen')}
    plot1d(results,period,out);maps(results,period,out)
    rows=[]
    for kind,r in results.items():
        for (pc,tc),n in zip(CUTS,r['counts']):
            rows.append(dict(period=period,sample=kind,p_min_GeV=pc,theta_max_deg=tc,n_survive=int(n),n_baseline=int(r['total']),fraction=float(n/max(r['total'],1))))
    with (out/'survival.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',default='output/proton_phase_space_root_diagnostics')
    ap.add_argument('--workers',type=int,default=5)
    ap.add_argument('--step-size',default='100 MB')
    ap.add_argument('--periods',nargs='+',choices=PERIODS,default=PERIODS)
    a=ap.parse_args();out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    jobs=[(p,str(out),a.step_size) for p in a.periods];rows=[]
    with cf.ProcessPoolExecutor(max_workers=min(8,max(1,a.workers),len(jobs))) as pool:
        for result in pool.map(worker,jobs):rows.extend(result)
    with (out/'survival_all_periods.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('Completed:',out)
if __name__=='__main__':main()
