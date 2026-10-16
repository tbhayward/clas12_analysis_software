#!/usr/bin/env python3
"""Independent fit audit: published RGA Pass-1 BSA vs Sebastian Alvarado Table II.

Both sources are fitted with identical sine-only and sine/(1+a2*cos) functions.
CAUTION: the published Pass-1 points and Sebastian Table II have DIFFERENT
3D binning. Pass-1 published points are assigned to Sebastian cells by their
published MEAN (Q2,xB,|t|). This is a diagnostic regrouping, NOT an exact
same-bin comparison or a replacement for a common event-level re-extraction.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent

def sine(phi,a):return a*np.sin(np.deg2rad(phi))
def two(phi,a,b):return a*np.sin(np.deg2rad(phi))/(1+b*np.cos(np.deg2rad(phi)))

def fit(g,kind,min_points):
    phi=g.phi.to_numpy(float); y=g.A.to_numpy(float);err=g.err.to_numpy(float)
    good=np.isfinite(phi)&np.isfinite(y)&np.isfinite(err)&(err>0)
    phi,y,err=phi[good],y[good],err[good]
    k=1 if kind=='sin' else 2
    r={'n':len(y),'a1':np.nan,'a1_err':np.nan,'a2':np.nan,'a2_err':np.nan,'chi2':np.nan,'ndf':len(y)-k,'chi2_ndf':np.nan,'boundary':False,'valid':False}
    if len(y)<max(min_points,k+2):return r
    fn=sine if kind=='sin' else two
    try:
        p,cov=curve_fit(fn,phi,y,sigma=err,absolute_sigma=True,p0=[0.2] if k==1 else [0.2,0],bounds=([-2] if k==1 else [-2,-.95],[2] if k==1 else [2,.95]),maxfev=40000)
        errs=np.sqrt(np.diag(cov)); chi=float(np.sum(((y-fn(phi,*p))/err)**2))
        r.update(a1=p[0],a1_err=errs[0],a2=p[1] if k==2 else np.nan,a2_err=errs[1] if k==2 else np.nan,chi2=chi,chi2_ndf=chi/(len(y)-k),boundary=k==2 and abs(p[1])>.94,valid=bool(np.all(np.isfinite(errs)) and np.all(errs>0) and (k==1 or abs(p[1])<.94)))
    except (ValueError,RuntimeError,TypeError,FloatingPointError):pass
    return r

def load_pass1(path):
    d=pd.read_csv(path,sep=r'\s+',comment='#',header=None,names=['phi_rad','Q2','xB','t','energy','A','err'])
    d['phi']=np.rad2deg(d.phi_rad)%360;d['t_abs']=d.t.abs()
    return d[np.isfinite(d.A)&np.isfinite(d.err)&(d.err>0)].copy()

def load_seb(path):
    d=pd.read_csv(path)
    need=['cell_id','Q2min','Q2max','xBmin','xBmax','t_abs_min','t_abs_max','phi','BSA','stat']
    if any(x not in d for x in need):raise ValueError('Sebastian CSV missing columns: '+str([x for x in need if x not in d]))
    d=d.rename(columns={'BSA':'A','stat':'err'});return d[np.isfinite(d.A)&np.isfinite(d.err)&(d.err>0)].copy()

def assign_pass1(p,s):
    # Assign by published mean kinematics, not by phi-bin boundaries.
    cols=['cell_id','Q2min','Q2max','xBmin','xBmax','t_abs_min','t_abs_max']
    cells=s[cols].drop_duplicates().sort_values('cell_id')
    if cells.cell_id.duplicated().any():raise ValueError('Ambiguous Sebastian cell boundaries')
    assigned=np.full(len(p),-1,dtype=int)
    for r in cells.itertuples(index=False):
        m=(p.Q2.to_numpy()>=r.Q2min)&(p.Q2.to_numpy()<r.Q2max)&(p.xB.to_numpy()>=r.xBmin)&(p.xB.to_numpy()<r.xBmax)&(p.t_abs.to_numpy()>=r.t_abs_min)&(p.t_abs.to_numpy()<r.t_abs_max)
        if np.any((assigned>=0)&m):raise ValueError('Overlapping Sebastian cell boundaries')
        assigned[m]=int(r.cell_id)
    p=p.copy();p['cell_id']=assigned
    return p,cells

def plot_overlays(p,s,cells,out,energy,kind,min_points):
    sub=out/f'{energy:g}_GeV'/kind;sub.mkdir(parents=True,exist_ok=True)
    fits=[]
    for r in cells.itertuples(index=False):
        cid=int(r.cell_id); a=p[(p.energy==energy)&(p.cell_id==cid)];b=s[s.cell_id==cid]
        if len(a)<min_points or len(b)<min_points:continue
        fa=fit(a,kind,min_points);fb=fit(b,kind,min_points)
        record={'energy':energy,'cell_id':cid,'fit_model':kind,'Q2min':r.Q2min,'Q2max':r.Q2max,'xBmin':r.xBmin,'xBmax':r.xBmax,'t_abs_min':r.t_abs_min,'t_abs_max':r.t_abs_max}
        for prefix,res in [('pass1',fa),('sebastian',fb)]:record.update({f'{prefix}_{k}':v for k,v in res.items()})
        record['delta_a1']=fb['a1']-fa['a1'];record['naive_delta_err']=np.hypot(fa['a1_err'],fb['a1_err']);record['naive_pull']=record['delta_a1']/record['naive_delta_err'] if record['naive_delta_err']>0 else np.nan
        record['qualified']=bool(fa['valid'] and fb['valid']);fits.append(record)
        fig,ax=plt.subplots(figsize=(7.5,4.6));ax.errorbar(a.phi,a.A,yerr=a.err,fmt='o',ms=3.5,color='#0072B2',label=f'Published Pass-1 ({len(a)} points)',alpha=.8);ax.errorbar(b.phi,b.A,yerr=b.err,fmt='s',ms=3.3,color='#D55E00',label=f'Sebastian Table II ({len(b)} points)',alpha=.8)
        grid=np.linspace(0,360,500);fn=sine if kind=='sin' else two
        for f,color,label in [(fa,'#0072B2','Pass-1'),(fb,'#D55E00','Sebastian')]:
            if f['valid']:
                params=[f['a1']] if kind=='sin' else [f['a1'],f['a2']]
                ax.plot(grid,fn(grid,*params),color=color,lw=1.7,label=f"{label} a₁={f['a1']:.3f}±{f['a1_err']:.3f}")
        ax.axhline(0,color='0.5',lw=.7);ax.set(xlim=(0,360),xlabel='φ (degrees)',ylabel='Beam-spin asymmetry',title=f'{energy:g} GeV | Sebastian cell {cid} | {kind} fit\nQ² [{r.Q2min:g},{r.Q2max:g}], xB [{r.xBmin:g},{r.xBmax:g}], |t| [{r.t_abs_min:g},{r.t_abs_max:g}]')
        ax.legend(fontsize=8,loc='best');ax.grid(alpha=.18);fig.tight_layout();fig.savefig(sub/f'cell_{cid:03d}.png',dpi=140);plt.close(fig)
    return fits

def summary(fits,out):
    d=pd.DataFrame(fits);d.to_csv(out/'all_fit_results.csv',index=False)
    qualified=d[d.qualified].copy();qualified.to_csv(out/'qualified_fit_results.csv',index=False)
    rows=[]
    for (energy,kind),g in d.groupby(['energy','fit_model']):
        q=g[g.qualified];row={'energy':energy,'fit_model':kind,'fit_pairs':len(g),'qualified_pairs':len(q),'sebastian_smaller':int((q.delta_a1<0).sum()),'mean_delta':q.delta_a1.mean(),'median_delta':q.delta_a1.median(),'median_ratio':np.median((q.sebastian_a1/q.pass1_a1).replace([np.inf,-np.inf],np.nan).dropna()) if len(q) else np.nan,'naive_chi2_per_pair':np.mean(q.naive_pull**2) if len(q) else np.nan};rows.append(row)
        if q.empty:continue
        folder=out/f'{energy:g}_GeV'/kind
        fig,ax=plt.subplots(figsize=(6,5));ax.errorbar(q.pass1_a1,q.sebastian_a1,xerr=q.pass1_a1_err,yerr=q.sebastian_a1_err,fmt='o',ms=4,alpha=.8);lo=min(q.pass1_a1.min(),q.sebastian_a1.min(),0);hi=max(q.pass1_a1.max(),q.sebastian_a1.max());ax.plot([lo,hi],[lo,hi],'k--',lw=1);ax.set(xlabel='Published Pass-1 a₁',ylabel='Sebastian Pass-2 a₁',title=f'{energy:g} GeV | {kind} | {len(q)} pairs');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(folder/'amplitude_scatter.png',dpi=160);plt.close(fig)
        fig,ax=plt.subplots(figsize=(8,4));q=q.sort_values('cell_id');ax.errorbar(q.cell_id,q.delta_a1,yerr=q.naive_delta_err,fmt='o',capsize=2);ax.axhline(0,color='k',ls='--');ax.set(xlabel='Sebastian Table II cell',ylabel='Sebastian − Pass-1 a₁',title=f'{energy:g} GeV | {kind}');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(folder/'amplitude_differences.png',dpi=160);plt.close(fig)
    sm=pd.DataFrame(rows);sm.to_csv(out/'summary.csv',index=False);print(sm.to_string(index=False));print(f'[output] {out}')

def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pass1',type=Path,default=HERE/'import'/'rga_pass1_bsa.txt')
    ap.add_argument('--sebastian',type=Path,default=HERE.parent/'imports'/'sebastian_tableII_proton_detected.csv')
    ap.add_argument('--output-dir',type=Path,default=HERE/'output'/'pass1_vs_sebastian_fit_audit')
    ap.add_argument('--min-points',type=int,default=6)
    args=ap.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    p=load_pass1(args.pass1);s=load_seb(args.sebastian);p,cells=assign_pass1(p,s)
    print(f'[inputs] Pass-1={len(p)} points, Sebastian={len(s)} points, Sebastian cells={len(cells)}')
    print(f'[mapping] Pass-1 points assigned to Sebastian cells: {int((p.cell_id>=0).sum())}/{len(p)}')
    print('[CAUTION] Published Pass-1 and Sebastian use different 3D bins. Regrouping by published point means is approximate, not a same-bin comparison.')
    allfits=[]
    for energy in [10.6,10.2]:
        for kind in ['sin','two']:allfits.extend(plot_overlays(p,s,cells,args.output_dir,energy,kind,args.min_points))
    if not allfits:raise RuntimeError('No matched cells with enough published Pass-1 points. Check input files and bin definitions.')
    summary(allfits,args.output_dir)

if __name__=='__main__':main()
