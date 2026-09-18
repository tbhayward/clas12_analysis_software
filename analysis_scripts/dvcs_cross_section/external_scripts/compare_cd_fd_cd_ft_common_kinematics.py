#!/usr/bin/env python3
"""Common-kinematics CD-FD / CD-FT topology comparison.

Reads already-produced topology CSVs. It does NOT rerun the DVCS analysis.
For each common physics bin and run period, evaluate KM15 total ep->epgamma and
pure BH at the measured mean kinematics of CD-FD and CD-FT separately.

Definitions:
  R_raw   = sigma_data(CD-FD) / sigma_data(CD-FT)
  R_kin   = sigma_model(CD-FD means) / sigma_model(CD-FT means)
  R_resid = R_raw / R_kin

Thus R_resid=1 means that the observed topology difference is fully accounted
for by the different measured mean kinematics according to that model.

Outputs CSV tables and PNG figures only. No analysis product is modified.
"""
from __future__ import annotations
import argparse, math, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DEFAULT_OUT=ROOT/'output'/'topology_common_kinematics'
PERIODS={
 'Fa18 Inb':10.604, 'Fa18 Out':10.604,
 'Sp18 Inb':10.594, 'Sp18 Out':10.594, 'Sp19 Inb':10.200,
}
T_EDGES=[0.0,0.15,0.25,0.40,0.60,np.inf]
T_LABELS=['<0.15','0.15-0.25','0.25-0.40','0.40-0.60','>=0.60']


def parse_xs(v):
    try:
        s=str(v).strip()
        if s.startswith('(') and s.endswith(')'):
            a=[float(x.strip()) for x in s[1:-1].split(',')]
            return (a+[np.nan,np.nan,np.nan])[:3]
        return [float(s),np.nan,np.nan]
    except Exception:
        return [np.nan,np.nan,np.nan]


def circdiff(a,b):
    return ((np.asarray(a,float)-np.asarray(b,float)+180.)%360.)-180.


def model_worker(task):
    i,xB,Q2,t,phi,E=task
    try:
        import gepard as g
        from gepard.fits import th_KM15
        phi_trento=math.pi-math.radians(float(phi)%360.)
        pt=g.DataPoint(xB=float(xB),t=-abs(float(t)),Q2=float(Q2),phi=phi_trento,
            observable='XS',frame='trento',process='ep2epgamma',exptype='fixed target',
            in1energy=float(E),in1charge=-1,in1polarization=0,in2particle='p')
        pt.prepare()
        th=th_KM15
        pref=float(th.PreFacSigma(pt))
        bh=pref*float(th.TBH2unp(pt))
        ep=bh+pref*float(th.TINTunp(pt))+pref*float(th.TDVCS2unp(pt))
        if not (np.isfinite(ep) and np.isfinite(bh) and ep>0 and bh>0):
            return i,np.nan,np.nan
        return i,ep,bh
    except Exception as e:
        return i,np.nan,np.nan


def load_topology(path, period, tag):
    d=pd.read_csv(path,low_memory=False)
    xscol=f'normed cross sections, ep->epg, exp, {period}, unpol'
    cols={v:f'{v}, {period}' for v in ['xBavg','Q2avg','t_abs_avg','phiavg','g_theta','p_theta']}
    req=['bin index','Bin Name','xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax',xscol,*cols.values()]
    miss=[c for c in req if c not in d.columns]
    if miss: raise RuntimeError(f'{tag} CSV missing columns for {period}: {miss}')
    q=d[req].copy()
    vals=[parse_xs(v) for v in q[xscol]]
    q[f'xs_{tag}']=[x[0] for x in vals]
    q[f'stat_{tag}']=[abs(x[1]) for x in vals]
    q[f'ptp_{tag}']=[abs(x[2]) for x in vals]
    for v,c in cols.items(): q[f'{v}_{tag}']=pd.to_numeric(q[c],errors='coerce')
    q[f'phiavg_{tag}']=np.mod(q[f'phiavg_{tag}'],360.)
    keep=['bin index','Bin Name','xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax',
          f'xs_{tag}',f'stat_{tag}',f'ptp_{tag}']+[f'{v}_{tag}' for v in cols]
    return q[keep]


def make_common(fd_path,ft_path,period):
    a=load_topology(fd_path,period,'fd'); b=load_topology(ft_path,period,'ft')
    keys=['bin index','Bin Name','xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax']
    m=a.merge(b,on=keys,how='inner',validate='one_to_one')
    finite=np.isfinite(m.xs_fd)&np.isfinite(m.xs_ft)&(m.xs_fd>0)&(m.xs_ft>0)
    m=m.loc[finite].copy()
    m['period']=period; m['ebeam']=PERIODS[period]
    m['R_raw']=m.xs_fd/m.xs_ft
    m['R_raw_stat']=m.R_raw*np.sqrt((m.stat_fd/m.xs_fd)**2+(m.stat_ft/m.xs_ft)**2)
    m['dxB_fd_minus_ft']=m.xBavg_fd-m.xBavg_ft
    m['dQ2_fd_minus_ft']=m.Q2avg_fd-m.Q2avg_ft
    m['dt_fd_minus_ft']=m.t_abs_avg_fd-m.t_abs_avg_ft
    m['dphi_fd_minus_ft']=circdiff(m.phiavg_fd,m.phiavg_ft)
    m['t_common']=0.5*(m.t_abs_avg_fd+m.t_abs_avg_ft)
    # Circular mean, robust here because the two means are in the same nominal phi bin.
    m['phi_common']=np.mod(m.phiavg_ft+0.5*m.dphi_fd_minus_ft,360.)
    return m


def add_models(d,workers):
    tasks=[]; meta=[]
    for j,r in d.iterrows():
        for tag in ['fd','ft']:
            tasks.append((len(tasks),r[f'xBavg_{tag}'],r[f'Q2avg_{tag}'],r[f't_abs_avg_{tag}'],r[f'phiavg_{tag}'],r.ebeam))
            meta.append((j,tag))
    print(f'[topology transport] evaluating {len(tasks)} mean-kinematics points with {workers} worker(s)',flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        results=list(ex.map(model_worker,tasks,chunksize=32))
    ep={}; bh={}
    for i,e,b in results:
        j,tag=meta[i]; ep[(j,tag)]=e; bh[(j,tag)]=b
    out=d.copy()
    for tag in ['fd','ft']:
        out[f'km15_{tag}']=[ep.get((j,tag),np.nan) for j in out.index]
        out[f'bh_{tag}']=[bh.get((j,tag),np.nan) for j in out.index]
    out['R_kin_km15']=out.km15_fd/out.km15_ft
    out['R_kin_bh']=out.bh_fd/out.bh_ft
    out['R_residual_km15']=out.R_raw/out.R_kin_km15
    out['R_residual_bh']=out.R_raw/out.R_kin_bh
    out['bh_fraction_fd']=out.bh_fd/out.km15_fd
    out['bh_fraction_ft']=out.bh_ft/out.km15_ft
    return out


def med(x):
    a=pd.to_numeric(x,errors='coerce').to_numpy(float); a=a[np.isfinite(a)]
    return float(np.median(a)) if len(a) else np.nan


def summarize_group(g,label):
    return dict(group=label,N=len(g),R_raw_median=med(g.R_raw),R_kin_km15_median=med(g.R_kin_km15),
        R_residual_km15_median=med(g.R_residual_km15),R_kin_bh_median=med(g.R_kin_bh),
        R_residual_bh_median=med(g.R_residual_bh),bh_fraction_fd_median=med(g.bh_fraction_fd),
        bh_fraction_ft_median=med(g.bh_fraction_ft),dxB_median=med(g.dxB_fd_minus_ft),
        dQ2_median=med(g.dQ2_fd_minus_ft),dt_median=med(g.dt_fd_minus_ft),dphi_median=med(g.dphi_fd_minus_ft))


def write_summaries(d,outdir):
    rows=[]
    for p,g in d.groupby('period',sort=False):
        rows.append(summarize_group(g,f'{p}: all'))
        cats=pd.cut(g.t_common,T_EDGES,labels=T_LABELS,right=False)
        for lab in T_LABELS:
            q=g.loc[cats==lab]
            if len(q): rows.append(summarize_group(q,f'{p}: |t| {lab}'))
    s=pd.DataFrame(rows); s.to_csv(outdir/'topology_transport_summary.csv',index=False)
    return s


def plot_phi(d,outdir):
    # Nominal phi bins are the cleanest aggregation: preserve all populated bins.
    for p,g in d.groupby('period',sort=False):
        q=g.copy(); q['phi_bin']=q.apply(lambda r:f'{r.phimin:g}-{r.phimax:g}',axis=1)
        agg=q.groupby(['phimin','phimax'],as_index=False).agg(
            phi=('phi_common','median'),raw=('R_raw','median'),km=('R_residual_km15','median'),bh=('R_residual_bh','median'),N=('R_raw','size')).sort_values('phi')
        if not len(agg): continue
        fig,ax=plt.subplots(figsize=(8.0,5.0))
        ax.plot(agg.phi,agg.raw,'o-',label='Raw CD-FD / CD-FT')
        ax.plot(agg.phi,agg.km,'o-',label='After KM15 mean-kinematics transport')
        ax.plot(agg.phi,agg.bh,'o--',label='After BH mean-kinematics transport')
        ax.axhline(1.0,linewidth=1)
        ax.set_xlabel(r'Mean $\phi$ (deg)'); ax.set_ylabel('CD-FD / CD-FT ratio')
        ax.set_title(f'{p}: topology ratio versus phi')
        ax.legend(); fig.tight_layout(); fig.savefig(outdir/f'phi_transport_{p.replace(" ","_")}.png',dpi=180); plt.close(fig)


def plot_t(d,outdir):
    for p,g in d.groupby('period',sort=False):
        cats=pd.cut(g.t_common,T_EDGES,labels=T_LABELS,right=False)
        rows=[]
        for lab in T_LABELS:
            q=g.loc[cats==lab]
            if len(q): rows.append((lab,med(q.R_raw),med(q.R_residual_km15),med(q.R_residual_bh),len(q)))
        if not rows: continue
        x=np.arange(len(rows)); labs=[r[0] for r in rows]
        fig,ax=plt.subplots(figsize=(8.0,5.0))
        ax.plot(x,[r[1] for r in rows],'o-',label='Raw')
        ax.plot(x,[r[2] for r in rows],'o-',label='KM15 transported')
        ax.plot(x,[r[3] for r in rows],'o--',label='BH transported')
        ax.axhline(1.0,linewidth=1); ax.set_xticks(x,labs); ax.set_xlabel(r'$|t|$ bin (GeV$^2$)'); ax.set_ylabel('Median CD-FD / CD-FT ratio')
        ax.set_title(f'{p}: common-kinematics topology comparison'); ax.legend(); fig.tight_layout(); fig.savefig(outdir/f't_transport_{p.replace(" ","_")}.png',dpi=180); plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cd-fd',required=True,type=Path)
    ap.add_argument('--cd-ft',required=True,type=Path)
    ap.add_argument('--outdir',type=Path,default=DEFAULT_OUT)
    ap.add_argument('--workers',type=int,default=8)
    a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    allrows=[]
    for p in PERIODS:
        q=make_common(a.cd_fd,a.cd_ft,p)
        print(f'[topology transport] {p}: {len(q)} common positive-cross-section bins',flush=True)
        allrows.append(q)
    d=pd.concat(allrows,ignore_index=True)
    d=add_models(d,max(1,a.workers))
    d.to_csv(a.outdir/'topology_common_kinematics_bin_by_bin.csv',index=False)
    s=write_summaries(d,a.outdir)
    plot_phi(d,a.outdir); plot_t(d,a.outdir)
    print('\n'+s.to_string(index=False,float_format=lambda x:f'{x:.5g}'))
    print(f'\n[topology transport] wrote outputs to {a.outdir}')

if __name__=='__main__': main()
