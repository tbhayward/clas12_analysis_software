#!/usr/bin/env python3
"""Publication-oriented plots from existing DVCS world/topology comparison CSV outputs.
No model evaluations or event-level analysis are performed.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':11,'axes.labelsize':12,'axes.titlesize':12,'legend.fontsize':9,'savefig.dpi':250,'pdf.fonttype':42,'axes.linewidth':1.0})
COLORS={'Inclusive':'#222222','CD-FD':'#0072B2','CD-FT':'#D55E00','FD-FD':'#009E73'}

def summarize(df, xcol, edges, ratio, error, label):
    rows=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        d=df.loc[(df[xcol]>=lo)&(df[xcol]<hi)].copy()
        y=pd.to_numeric(d[ratio],errors='coerce').to_numpy(float)
        s=pd.to_numeric(d[error],errors='coerce').to_numpy(float) if error else np.full(len(d),np.nan)
        x=pd.to_numeric(d[xcol],errors='coerce').to_numpy(float)
        good=np.isfinite(y)&(y>0)&np.isfinite(x)
        y,s,x=y[good],s[good],x[good]
        if len(y)==0:continue
        w=np.where(np.isfinite(s)&(s>0),1/np.maximum(s,1e-20)**2,0)
        if not np.any(w):w=np.ones(len(y))
        rows.append(dict(selection=label,variable=xcol,lower=lo,upper=hi,N=len(y),x_mean=np.mean(x),weighted_mean=np.average(y,weights=w),weighted_stat=np.sqrt(1/w.sum()),median=np.median(y),p16=np.percentile(y,16),p84=np.percentile(y,84)))
    return pd.DataFrame(rows)

def save(fig,path):
    fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight')
    fig.savefig(path.with_suffix('.png'),bbox_inches='tight',dpi=250)
    plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--world-dir',type=Path,required=True)
    p.add_argument('--topology-dir',type=Path,required=True)
    p.add_argument('--common-dir',type=Path)
    p.add_argument('--outdir',type=Path,default=Path('output/world_cross_section_note_figures'))
    p.add_argument('--bins',type=int,default=7)
    a=p.parse_args();a.outdir.mkdir(parents=True,exist_ok=True)
    f=a.world_dir/'tables'/'lee_hayward_exact_same_bin_points.csv'
    if not f.exists():
        matches=list(a.world_dir.rglob('lee_hayward_exact_same_bin_points.csv'))
        if not matches:raise FileNotFoundError(f'Cannot find exact-bin CSV under {a.world_dir}. Check --world-dir and confirm the world-comparison output is present.')
        f=matches[0]
    print('[world] reading',f)
    d=pd.read_csv(f)
    vars=[('hayward_Q2',r'$Q^2$ (GeV$^2$)'),('hayward_xB',r'$x_B$'),('hayward_t_abs',r'$|t|$ (GeV$^2$)')]
    summary=[]
    fig,axes=plt.subplots(1,3,figsize=(14.5,4.1),layout='constrained')
    for ax,(key,xlabel) in zip(axes,vars):
        v=pd.to_numeric(d[key],errors='coerce').to_numpy(float);v=v[np.isfinite(v)]
        edges=np.unique(np.quantile(v,np.linspace(0,1,a.bins+1)))
        edges[-1]=np.nextafter(edges[-1],np.inf)
        for col,lab,marker,shift in [('hayward_over_lee_direct','Direct same-bin','o',-0.008),('hayward_over_lee_after_km15_mean_transport','KM15 mean transport','s',0.008)]:
            # Statistical-only ratio errors; correlated normalization uncertainties are excluded.
            q=d.copy()
            q['_stat_ratio']=pd.to_numeric(q[col],errors='coerce')*np.sqrt((pd.to_numeric(q.hayward_stat,errors='coerce')/pd.to_numeric(q.hayward_xs,errors='coerce'))**2+(pd.to_numeric(q.lee_stat,errors='coerce')/pd.to_numeric(q.lee_xs,errors='coerce'))**2)
            t=summarize(q,key,edges,col,'_stat_ratio',lab)
            summary.append(t)
            if t.empty:continue
            width=np.median(np.diff(edges));xx=t.x_mean.to_numpy()+shift*width
            ax.errorbar(xx,t.weighted_mean,yerr=t.weighted_stat,fmt=marker,ms=5,capsize=2,label=lab)
        ax.axhline(1,color='0.4',ls='--',lw=1)
        ax.set_xlabel(xlabel);ax.set_ylabel('Pass-2 / Pass-1' if ax is axes[0] else '')
        ax.grid(alpha=.17);ax.set_ylim(bottom=0)
    axes[0].legend(frameon=False,loc='best')
    save(fig,a.outdir/'pass2_pass1_vs_kinematics')
    if summary:pd.concat(summary,ignore_index=True).to_csv(a.outdir/'pass2_pass1_binned_summary.csv',index=False)
    # Topology Pass-1 comparison: use the original matching output, not guessed CSV schema.
    f=a.topology_dir/'pass2_vs_authoritative_pass1_matched_points.csv'
    if not f.exists():
        matches=list(a.topology_dir.rglob('pass2_vs_authoritative_pass1_matched_points.csv'))
        if matches:f=matches[0]
    if f.exists():
        t=pd.read_csv(f)
        print('[topology] columns:',list(t.columns))
        # Discover expected columns conservatively; fail visibly rather than plotting incorrect quantities.
        def col(names):return next((n for n in names if n in t.columns),None)
        typ=col(['topology','selection','topology_label']);ratio=col(['pass2_over_pass1','ratio_pass2_pass1','ratio','pass2_pass1_ratio']);tabs=col(['t_abs','t_abs_avg','t','abs_t','t_pass2'])
        if typ and ratio and tabs:
            fig,ax=plt.subplots(figsize=(7.4,4.8),layout='constrained')
            rows=[]
            for tag in ['Inclusive','CD-FD','CD-FT','FD-FD']:
                q=t.loc[t[typ].astype(str).str.lower()==tag.lower()]
                if q.empty:continue
                edges=np.array([0,.2,.3,.4,1.0,np.inf]);s=summarize(q,tabs,edges,ratio,None,tag)
                rows.append(s)
                ax.plot(s.x_mean,s['median'],'o-',label=f'{tag} (N={len(q)})',color=COLORS[tag],lw=1.6)
            ax.axhline(1,color='0.4',ls='--');ax.set(xlabel=r'$|t|$ (GeV$^2$)',ylabel='Median Pass-2 / Pass-1')
            ax.grid(alpha=.17);ax.legend(frameon=False,ncol=2)
            save(fig,a.outdir/'pass2_pass1_by_topology')
            pd.concat(rows,ignore_index=True).to_csv(a.outdir/'topology_binned_summary.csv',index=False)
        else:
            print('[warning] topology CSV schema not recognized; topology plot not generated. Provide column names above.')
    if a.common_dir:
        f=a.common_dir/'topology_common_kinematics_bin_by_bin.csv'
        if not f.exists():
            matches=list(a.common_dir.rglob('topology_common_kinematics_bin_by_bin.csv'))
            if matches:f=matches[0]
        if f.exists():
            t=pd.read_csv(f)
            fig,ax=plt.subplots(figsize=(7.4,4.8),layout='constrained')
            for label,cut in [('All common bins',np.ones(len(t),bool)),('Tight mean-kinematics', (t.dxB_fd_minus_ft.abs()<.01)&(t.dQ2_fd_minus_ft.abs()<.10))]:
                q=t.loc[cut&np.isfinite(t.R_raw)&(t.R_raw>0)]
                if q.empty:continue
                edges=np.array([0,.15,.25,.4,.6,np.inf]);s=summarize(q,'t_common',edges,'R_raw','R_raw_stat',label)
                if s.empty:continue
                ax.plot(s.x_mean,s['median'],'o-',label=f'{label} (N={len(q)})')
                s.to_csv(a.outdir/('common_tight_summary.csv' if label.startswith('Tight') else 'common_all_summary.csv'),index=False)
            ax.axhline(1,color='0.4',ls='--');ax.set(xlabel=r'$|t|$ (GeV$^2$)',ylabel='Median CD-FD / CD-FT')
            ax.grid(alpha=.17);ax.legend(frameon=False)
            save(fig,a.outdir/'cd_fd_cd_ft_tight_kinematics')
            rows=[]
            for dx,dq in [(.04,.4),(.03,.3),(.02,.2),(.01,.1),(.005,.05)]:
                q=t.loc[(t.dxB_fd_minus_ft.abs()<dx)&(t.dQ2_fd_minus_ft.abs()<dq)&np.isfinite(t.R_raw)&(t.R_raw>0)]
                rows.append(dict(max_abs_dxB=dx,max_abs_dQ2=dq,N=len(q),median_raw_ratio=q.R_raw.median(),median_km15_residual=q.R_residual_km15.median()))
            pd.DataFrame(rows).to_csv(a.outdir/'kinematics_matching_threshold_scan.csv',index=False)
    print('[done]',a.outdir)

if __name__=='__main__':main()
