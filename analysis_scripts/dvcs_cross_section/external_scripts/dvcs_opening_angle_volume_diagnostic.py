#!/usr/bin/env python3
"""Standalone DVCS opening-angle phase-space diagnostic (v3).

Reads imports/all_bin_v3.csv directly, without any analysis output CSV.
Computes surviving geometric and optionally BH/KM15 cross-section-weighted
fractions for theta(e',gamma) > --angle degrees.

Dependencies: numpy, pandas, scipy, matplotlib; optional gepard for models.
Run from dvcs_cross_section:
  python3 external_scripts/dvcs_opening_angle_volume_diagnostic_v3.py --workers 4 --samples 4096
  python3 external_scripts/dvcs_opening_angle_volume_diagnostic_v3.py --workers 4 --model bh --samples 1024

No existing repository files are modified.
"""
from __future__ import annotations
import argparse
import math
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
from scipy.stats import qmc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

M = 0.9382720813
EDGES = ('xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax')
THREE_D = EDGES[:6]


def locate_csv(user_path):
    if user_path:
        p = Path(user_path).expanduser().resolve()
        if not p.is_file(): raise FileNotFoundError(p)
        return p
    here = Path(__file__).resolve().parent
    candidates = [here/'imports/all_bin_v3.csv', here.parent/'imports/all_bin_v3.csv',
                  here.parent.parent/'imports/all_bin_v3.csv', Path.cwd()/'imports/all_bin_v3.csv',
                  Path.cwd()/'analysis_scripts/dvcs_cross_section/imports/all_bin_v3.csv']
    for p in candidates:
        if p.is_file(): return p.resolve()
    raise FileNotFoundError('Cannot locate imports/all_bin_v3.csv; pass --csv PATH')


def load_bins(path, all_rows):
    df = pd.read_csv(path, low_memory=False)
    missing = [c for c in EDGES if c not in df.columns]
    if missing: raise ValueError(f'Missing columns in input CSV: {missing}')
    if not all_rows:
        if 'valid bin' not in df.columns:
            raise ValueError('No "valid bin" column; use --all-rows to bypass filtering')
        df = df.loc[pd.to_numeric(df['valid bin'], errors='coerce').eq(1)].copy()
    for c in EDGES: df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.dropna(subset=list(EDGES)).drop_duplicates(subset=list(EDGES))
    df = df.loc[(df.xBmax>df.xBmin)&(df.Q2max>df.Q2min)&
                (df.t_abs_max>df.t_abs_min)&(df.phimax>df.phimin)&
                ((df.phimax-df.phimin)<=360)].copy()
    return df.reset_index(drop=True)


def kinematic_mask(x, q2, t, e):
    y = q2/(2*M*x*e)
    w2 = M*M + q2*(1/x-1)
    eps2 = 4*M*M*x*x/q2
    tmin = q2*(2*(1-x)*(1-np.sqrt(1+eps2))+eps2)/(4*x*(1-x)+eps2)
    return (x>0)&(x<1)&(q2>0)&(y>0.19)&(y<0.8113919276)&(w2>4)&(t>=tmin)


def electron_photon_angle(x, q2, t, phi_deg, e):
    """Lab angle from ep->e'p'gamma, neglecting electron mass.

    Incoming electron along +z, outgoing electron in xz plane. Phi is
    the azimuth about virtual photon relative to the lepton plane.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        nu = q2/(2*M*x)
        ee = e-nu
        cos_e = 1-q2/(2*e*ee)
        valid_e = (ee>0)&(np.abs(cos_e)<=1)
        cos_e = np.clip(cos_e,-1,1)
        sin_e = np.sqrt(np.maximum(0,1-cos_e*cos_e))
        qx = -ee*sin_e
        qz = e-ee*cos_e
        qabs = np.hypot(qx,qz)
        eg = nu-t/(2*M)
        # (q-q')^2=-t implies q.q'=(t-Q2)/2.
        cos_qg = (nu*eg+(q2-t)/2)/(qabs*eg)
        valid = valid_e&(eg>0)&(np.abs(cos_qg)<=1+1e-10)
        cos_qg = np.clip(cos_qg,-1,1)
        sin_qg = np.sqrt(np.maximum(0,1-cos_qg*cos_qg))
        ux, uz = qx/qabs, qz/qabs
        cp = np.cos(np.deg2rad(phi_deg))
        gx = cos_qg*ux+sin_qg*uz*cp
        gz = cos_qg*uz-sin_qg*ux*cp
        cos_eg = sin_e*gx+cos_e*gz
        angle = np.rad2deg(np.arccos(np.clip(cos_eg,-1,1)))
    return angle, valid


def model_weights(x,q2,t,phi,e,model):
    import gepard as g
    from gepard.fits import th_KM15
    out = np.full(len(x),np.nan)
    for j in range(len(x)):
        try:
            point = g.DataPoint(xB=float(x[j]),Q2=float(q2[j]),t=-float(t[j]),
                phi=math.pi-math.radians(float(phi[j])),observable='XS',
                frame='trento',process='ep2epgamma',exptype='fixed target',
                in1energy=float(e),in1charge=-1,in1polarization=0,in2particle='p')
            point.prepare()
            pref = float(th_KM15.PreFacSigma(point))
            bh = pref*float(th_KM15.TBH2unp(point))
            sigma = bh if model=='bh' else bh+pref*(float(th_KM15.TINTunp(point))+float(th_KM15.TDVCS2unp(point)))
            if np.isfinite(sigma) and sigma>=0: out[j]=sigma
        except Exception:
            pass
    return out


def calculate(task):
    index, edge_values, energy, n, angle_cut, model = task
    edges = np.asarray(edge_values,dtype=float).reshape(4,2)
    sampler = qmc.Sobol(d=4,scramble=True,seed=19001+index+int(energy*1000))
    points = qmc.scale(sampler.random_base2(int(math.log2(n))),edges[:,0],edges[:,1])
    x,q2,t,phi = points.T
    good = kinematic_mask(x,q2,t,energy)
    angle,physical = electron_photon_angle(x,q2,t,phi,energy)
    good &= physical
    selected = good&(angle>angle_cut)
    cube = float(np.prod(edges[:,1]-edges[:,0])*math.pi/180)
    volume = cube*float(np.mean(good))
    cut_volume = cube*float(np.mean(selected))
    result = dict(row_index=index,beam_energy=energy,**dict(zip(EDGES,edge_values)),
        V_nominal=volume,V_cut=cut_volume,R_volume=cut_volume/volume if volume>0 else np.nan,
        n_total=n,n_kin=int(good.sum()),n_cut=int(selected.sum()))
    if model!='none':
        ids = np.flatnonzero(good)
        weights = model_weights(x[ids],q2[ids],t[ids],phi[ids],energy,model)
        valid_weights = np.isfinite(weights)
        result['n_model_valid'] = int(valid_weights.sum())
        result['n_model_invalid'] = int((~valid_weights).sum())
        result['R_model'] = (float(weights[angle[ids]>angle_cut].sum()/weights.sum())
                             if len(ids) and valid_weights.all() and weights.sum()>0 else np.nan)
    return result


def save_plots(df,out,model,max_plots):
    for energy, group in df.groupby('beam_energy'):
        values = group.R_volume.dropna()
        fig,ax = plt.subplots(figsize=(7,4.5))
        ax.hist(values,bins=np.linspace(0,1,51))
        ax.set(xlabel=r'$V(\theta_{e\prime\gamma}>8^\circ)/V_{\rm nominal}$',
               ylabel='Number of bins',title=f'{energy:g} GeV')
        fig.tight_layout();fig.savefig(out/f'volume_ratio_hist_{energy:g}.png',dpi=170);plt.close(fig)
        cells = group.groupby(list(THREE_D),sort=True)
        for j,(_,cell) in enumerate(cells):
            if max_plots and j>=max_plots: break
            if len(cell)<2: continue
            cell = cell.sort_values('phimin')
            ph = (cell.phimin.to_numpy()+cell.phimax.to_numpy())/2
            fig,ax = plt.subplots(figsize=(7,3.8))
            ax.plot(ph,cell.R_volume.to_numpy(),'o-',label='Geometric')
            if model!='none': ax.plot(ph,cell.R_model.to_numpy(),'s--',label=model.upper())
            ax.set(xlabel=r'$\phi$ (degrees)',ylabel='Surviving fraction',xlim=(0,360),ylim=(-0.03,1.04))
            ax.set_title(f'{energy:g} GeV; xB [{cell.xBmin.iloc[0]:.3f},{cell.xBmax.iloc[0]:.3f}], '
                f'Q² [{cell.Q2min.iloc[0]:.2f},{cell.Q2max.iloc[0]:.2f}], '
                f'|t| [{cell.t_abs_min.iloc[0]:.2f},{cell.t_abs_max.iloc[0]:.2f}]',fontsize=9)
            if model!='none':ax.legend()
            fig.tight_layout();fig.savefig(out/f'phi_E{energy:g}_cell{j:04d}.png',dpi=140);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--csv',help='Path to imports/all_bin_v3.csv; autodetected if omitted')
    parser.add_argument('--output',default='output/opening_angle_diagnostic_v3')
    parser.add_argument('--energies',type=float,nargs='+',default=[10.604,10.200])
    parser.add_argument('--angle',type=float,default=8.)
    parser.add_argument('--samples',type=int,default=4096,help='Power of two samples per bin')
    parser.add_argument('--model',choices=['none','bh','km15'],default='none')
    parser.add_argument('--workers',type=int,default=4,help='1 to 4 worker processes')
    parser.add_argument('--max-rows',type=int,default=0,help='Limit input bins per energy for tests')
    parser.add_argument('--all-rows',action='store_true',help='Include invalid bins as well')
    parser.add_argument('--max-phi-plots',type=int,default=0,help='Limit plots per energy (0 means all)')
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args()
    if args.samples<2 or args.samples&(args.samples-1):parser.error('--samples must be a power of two')
    if not 1<=args.workers<=4:parser.error('--workers must be 1, 2, 3 or 4')
    if args.model!='none':
        try:import gepard
        except ImportError:parser.error('The bh/km15 model requires the Python gepard package')
    source=locate_csv(args.csv)
    bins=load_bins(source,args.all_rows)
    if args.max_rows: bins=bins.head(args.max_rows)
    tasks=[(i,tuple(float(row[c]) for c in EDGES),e,args.samples,args.angle,args.model)
           for e in args.energies for i,row in bins.iterrows()]
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    print(f'[angle] Source: {source}\n[angle] {len(bins)} bins, {len(tasks)} tasks, '
          f'{args.workers} workers, {args.samples} samples/bin, model={args.model}',flush=True)
    rows=[]
    if args.workers==1:
        for j,task in enumerate(tasks,1):
            rows.append(calculate(task))
            if j%100==0 or j==len(tasks):print(f'[angle] {j}/{len(tasks)}',flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(calculate,task) for task in tasks]
            for j,future in enumerate(as_completed(futures),1):
                rows.append(future.result())
                if j%100==0 or j==len(tasks):print(f'[angle] {j}/{len(tasks)}',flush=True)
    if not rows:raise SystemExit('No valid input bins')
    df=pd.DataFrame(rows).sort_values(['beam_energy','row_index']).reset_index(drop=True)
    df.to_csv(out/'bin_ratios.csv',index=False)
    summaries=[]
    for energy,group in df.groupby('beam_energy'):
        v=group.R_volume.dropna()
        summaries.append(f'Energy {energy:g} GeV: {len(group)} bins, {len(v)} nonempty')
        if len(v):
            summaries.append(f'  Geometric min={v.min():.6f}, median={v.median():.6f}, max={v.max():.6f}')
            summaries.extend(f'  R_volume < {cutoff:.2f}: {int((v<cutoff).sum())}' for cutoff in (0.99,0.95,0.90))
        if args.model!='none':
            w=group.R_model.dropna()
            summaries.append(f'  {args.model.upper()} weighted: {len(w)} valid bins')
            if len(w):summaries.append(f'  Weighted min={w.min():.6f}, median={w.median():.6f}, max={w.max():.6f}')
        summaries.append('')
    report='\n'.join(summaries)
    (out/'summary.txt').write_text(report+'\n')
    print(report,flush=True)
    if not args.no_plots:save_plots(df,out,args.model,args.max_phi_plots)
    print(f'[angle] Output: {out.resolve()}',flush=True)

if __name__=='__main__':main()
