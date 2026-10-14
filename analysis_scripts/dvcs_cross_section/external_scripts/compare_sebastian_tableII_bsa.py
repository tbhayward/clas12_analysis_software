#!/usr/bin/env python3
"""Compare event-level Fa18 Pass-2 epgamma BSA with Alvarado Table II.

No theoretical curves or geometric phi substitution. All fit inputs use measured
mean phi. The common-bin fits use exactly the same populated azimuthal slots.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KEY = ['Q2min','Q2max','xBmin','xBmax','t_abs_min','t_abs_max','phi_slot']
CELL = KEY[:-1]
COLORS = {'pub':'#1764ad','our':'#d62728'}

def model(phi,a1,a2):
    return a1*np.sin(np.deg2rad(phi))/(1+a2*np.cos(np.deg2rad(phi)))

def parse_tuple(x):
    if pd.isna(x): return np.nan,np.nan
    try:
        a,b=str(x).strip().strip('()').split(',')[:2]
        return float(a),float(b)
    except (ValueError,IndexError): return np.nan,np.nan

def valid(df,prefix):
    return np.isfinite(df[f'{prefix}_bsa']) & np.isfinite(df[f'{prefix}_stat']) & (df[f'{prefix}_stat']>0) & np.isfinite(df[f'{prefix}_phi'])

def fit(df,prefix,mask=None):
    ok=valid(df,prefix)
    if mask is not None: ok=ok & mask
    d=df.loc[ok]
    result={'a1':np.nan,'a2':np.nan,'a1_err':np.nan,'a2_err':np.nan,'chi2_ndf':np.nan,'n':len(d),'boundary':False}
    if len(d)<5:return result
    try:
        p,cov=curve_fit(model,d[f'{prefix}_phi'].to_numpy(),d[f'{prefix}_bsa'].to_numpy(),
                        sigma=d[f'{prefix}_stat'].to_numpy(),absolute_sigma=True,
                        p0=[0.2,0],bounds=([-2,-0.95],[2,0.95]),maxfev=30000)
        errs=np.sqrt(np.diag(cov))
        chi=np.sum(((d[f'{prefix}_bsa'].to_numpy()-model(d[f'{prefix}_phi'].to_numpy(),*p))/d[f'{prefix}_stat'].to_numpy())**2)
        result.update(a1=p[0],a2=p[1],a1_err=errs[0],a2_err=errs[1],chi2_ndf=chi/(len(d)-2),boundary=abs(p[1])>0.94)
    except (ValueError,RuntimeError,TypeError):pass
    return result

def save(fig,path):
    fig.savefig(path,dpi=170,bbox_inches='tight');plt.close(fig)

def cell_label(row):
    return (f"Q² [{row['Q2min']:g},{row['Q2max']:g}], "
            f"xB [{row['xBmin']:g},{row['xBmax']:g}], "
            f"|t| [{row['t_abs_min']:g},{row['t_abs_max']:g}]")

def load(our,reference):
    raw=pd.read_csv(our,low_memory=False);ref=pd.read_csv(reference)
    phi_col='phiavg, Fa18'
    if phi_col not in raw:raise KeyError(f'Missing {phi_col!r} in {our}')
    raw['our_phi']=pd.to_numeric(raw[phi_col],errors='coerce')
    if raw['our_phi'].notna().sum()==0:raise ValueError('Measured phiavg, Fa18 is empty; C++ bin-means stage must run.')
    for name in ['phimin','phimax']:
        raw[name]=pd.to_numeric(raw[name],errors='raise')
    slot=raw['phimin']/22.5
    if not np.allclose(slot,np.rint(slot),atol=1e-7) or not np.allclose(raw['phimax']-raw['phimin'],22.5,atol=1e-7):
        raise ValueError('Unexpected azimuthal boundaries: expected 22.5-degree bins')
    raw['phi_slot']=np.rint(slot).astype(int)
    for k in CELL:
        raw[k]=pd.to_numeric(raw[k],errors='raise');ref[k]=pd.to_numeric(ref[k],errors='raise')
    for period,short in [('Fa18','our'),('Fa18 Inb','inb'),('Fa18 Out','out')]:
        col=f'BSA, counts, {period}'
        if col not in raw:raise KeyError(f'Missing {col}')
        pairs=np.array([parse_tuple(x) for x in raw[col]],dtype=float)
        raw[f'{short}_bsa']=pairs[:,0];raw[f'{short}_stat']=pairs[:,1]
    ref['pub_phi']=pd.to_numeric(ref['phi'],errors='coerce')
    ref['pub_bsa']=pd.to_numeric(ref['BSA'],errors='coerce')
    ref['pub_stat']=pd.to_numeric(ref['stat'],errors='coerce')
    if raw.duplicated(KEY).any() or ref.duplicated(KEY).any():raise ValueError('Duplicate (cell,phi_slot) keys')
    keep=KEY+['our_phi','our_bsa','our_stat','inb_bsa','inb_stat','out_bsa','out_stat']
    m=ref.merge(raw[keep],on=KEY,how='left',validate='one_to_one',indicator=True)
    m['pub_valid']=valid(m,'pub');m['our_valid']=valid(m,'our')
    m['common_valid']=m['pub_valid'] & m['our_valid']
    m['delta']=m['our_bsa']-m['pub_bsa']
    m['delta_phi']=((m['our_phi']-m['pub_phi']+180)%360)-180
    m['flagged']=False
    m['independent_error_pull']=m['delta']/np.hypot(m['our_stat'],m['pub_stat'])
    return m,raw,ref

def overview(m,out):
    common=m.loc[m.common_valid]
    if common.empty:raise ValueError('No common finite BSA measurements')
    fig,ax=plt.subplots(figsize=(7,6))
    ax.scatter(common.pub_bsa,common.our_bsa,s=13,alpha=.55,color=COLORS['our'])
    lim=max(.2,np.nanpercentile(np.abs(np.r_[common.pub_bsa,common.our_bsa]),99)*1.15)
    ax.plot([-lim,lim],[-lim,lim],'k--',lw=1,label='Equality')
    ax.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel='Alvarado Table II BSA',ylabel='Our Fa18 BSA',
           title=f'Point-by-point BSA comparison ({len(common)} common valid points)')
    ax.legend();ax.grid(alpha=.2);save(fig,out/'01_bsa_scatter.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    axes[0].hist(common.delta,bins=50,color=COLORS['our'],alpha=.85)
    axes[0].axvline(0,color='k',ls='--');axes[0].set(xlabel='Our BSA − Table II BSA',ylabel='Number of points',title='Point-by-point differences')
    axes[1].scatter(common.pub_phi,common.delta,s=12,alpha=.5,color=COLORS['our'])
    axes[1].axhline(0,color='k',ls='--');axes[1].set(xlabel='Published measured φ (deg)',ylabel='Our BSA − Table II BSA',title='Difference versus azimuth',xlim=(0,360))
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('BSA differences (shared events; not independent measurements)');save(fig,out/'02_bsa_differences.png')

def amplitudes(m,out):
    rows=[]
    for cell,d in m.groupby(CELL,sort=True):
        common=d.common_valid & ~d.flagged
        p_all=fit(d,'pub');o_all=fit(d,'our')
        p=fit(d,'pub',common);o=fit(d,'our',common)
        row=dict(zip(CELL,cell));row['cell_id']=int(d.cell_id.iloc[0]) if 'cell_id' in d and pd.notna(d.cell_id.iloc[0]) else np.nan
        row.update(n_published=len(d),n_published_valid=int(d.pub_valid.sum()),n_our_valid=int(d.our_valid.sum()),n_common=int(common.sum()),n_flagged=int(d.flagged.sum()),n_common_before_flags=int(d.common_valid.sum()))
        for tag,res in [('published_all',p_all),('our_all',o_all),('published_common',p),('our_common',o)]:
            for k,v in res.items():row[f'{tag}_{k}']=v
        for param in ['a1','a2']:
            row[f'delta_{param}_common']=o[param]-p[param]
            row[f'naive_delta_{param}_err']=np.hypot(o[f'{param}_err'],p[f'{param}_err'])
        rows.append(row)
    s=pd.DataFrame(rows);s.to_csv(out/'amplitude_comparison.csv',index=False)
    for param in ['a1','a2']:
        ok=np.isfinite(s[f'published_common_{param}']) & np.isfinite(s[f'our_common_{param}'])
        d=s.loc[ok].copy()
        if d.empty:continue
        fig,axes=plt.subplots(1,2,figsize=(13,5))
        axes[0].errorbar(d[f'published_common_{param}'],d[f'our_common_{param}'],
                         xerr=d[f'published_common_{param}_err'],yerr=d[f'our_common_{param}_err'],
                         fmt='o',ms=3.5,alpha=.65,color=COLORS['our'],ecolor='0.7',elinewidth=.7)
        vals=np.r_[d[f'published_common_{param}'],d[f'our_common_{param}']]
        lo,hi=np.nanmin(vals),np.nanmax(vals);pad=max(.05,(hi-lo)*.1)
        axes[0].plot([lo-pad,hi+pad],[lo-pad,hi+pad],'k--',lw=1)
        axes[0].set(xlim=(lo-pad,hi+pad),ylim=(lo-pad,hi+pad),xlabel=f'Table II {param}',ylabel=f'Our {param}',title='Common-φ-bin fits')
        idx=np.arange(len(d))
        axes[1].errorbar(idx,d[f'delta_{param}_common'],yerr=d[f'naive_delta_{param}_err'],fmt='o',ms=3,color=COLORS['our'],ecolor='0.65',elinewidth=.7)
        axes[1].axhline(0,color='k',ls='--');axes[1].set(xlabel='Kinematic cell (sorted by Q², xB, |t|)',ylabel=f'Our {param} − Table II {param}',title='Amplitude differences (errors ignore covariance)')
        for ax in axes:ax.grid(alpha=.2)
        fig.suptitle(f'{param} comparison: {len(d)} cells with both fits; common measured φ bins')
        save(fig,out/f'03_amplitude_{param}_comparison.png')
    return s

def missing(m,s,out):
    d=s.copy();d['missing']=d.n_published-d.n_common
    fig,axes=plt.subplots(1,2,figsize=(13,4.7))
    axes[0].hist(d.n_common,bins=np.arange(-.5,17.6,1),color=COLORS['our'])
    axes[0].set(xlabel='Number of common valid φ bins per cell',ylabel='Number of cells',title='Common-bin coverage')
    idx=np.arange(len(d));axes[1].bar(idx,d.n_common,label='Common valid',color=COLORS['our'])
    axes[1].bar(idx,d.missing,bottom=d.n_common,label='Published without our valid counterpart',color='0.7')
    axes[1].set(xlabel='Kinematic cell (sorted by Q², xB, |t|)',ylabel='Published φ bins',title='Coverage by cell')
    axes[1].legend(fontsize=8)
    fig.suptitle('Coverage is not a fit-quality statistic');save(fig,out/'04_coverage.png')
    miss=m.loc[~m.common_valid].copy();miss.to_csv(out/'missing_or_invalid_points.csv',index=False)
    return miss

def cell_plot(d,path,cell_title):
    fig,(ax,axd)=plt.subplots(2,1,figsize=(7.5,6.2),sharex=True,gridspec_kw={'height_ratios':[2.5,1]})
    for prefix,label,col,mark in [('pub','Alvarado Table II',COLORS['pub'],'o'),('our','Our Fa18 Pass-2',COLORS['our'],'s')]:
        ok=valid(d,prefix)
        ax.errorbar(d.loc[ok,f'{prefix}_phi'],d.loc[ok,f'{prefix}_bsa'],yerr=d.loc[ok,f'{prefix}_stat'],
                    fmt=mark,ms=4,capsize=2,color=col,label=f'{label} (N={ok.sum()})')
        result=fit(d,prefix,d.common_valid & ~d.flagged)
        if np.isfinite(result['a1']):
            phi=np.linspace(0,360,361)
            ax.plot(phi,model(phi,result['a1'],result['a2']),color=col,lw=1.5,
                    label=f"{label} common-bin fit: a₁={result['a1']:.3f}±{result['a1_err']:.3f}")
    c=d.loc[d.common_valid & ~d.flagged]
    if len(c):
        # Conservative independent-error bars for visual scale only; shared-event covariance unknown.
        axd.errorbar(c.pub_phi,c.delta,yerr=np.hypot(c.pub_stat,c.our_stat),fmt='o',ms=3.5,
                     color=COLORS['our'],ecolor='0.6',capsize=2)
    for slot in d.loc[d.flagged,'phi_slot']:
        ax.axvspan(slot*22.5,(slot+1)*22.5,color='0.85',alpha=.5)
    ax.axhline(0,color='0.4',lw=.8);axd.axhline(0,color='0.4',lw=.8,ls='--')
    ax.set(ylabel=r'$A_{LU}$',title=cell_title);ax.legend(fontsize=7,loc='best');ax.grid(alpha=.2)
    axd.set(xlabel=r'Measured $\phi$ (degrees)',ylabel='Our − Table II',xlim=(0,360));axd.grid(alpha=.2)
    fig.tight_layout();save(fig,path)

def diagnostics(m,s,out,all_cells=False,top_n=12):
    detail=out/'individual_cells';detail.mkdir(exist_ok=True)
    candidates=s.loc[(s.n_common>=5)&np.isfinite(s.delta_a1_common)].copy()
    candidates['abs_delta']=candidates.delta_a1_common.abs()
    selected=set()
    if len(candidates):
        selected.update(candidates.nsmallest(min(3,len(candidates)),'abs_delta').index)
        selected.update(candidates.nlargest(min(top_n-3,len(candidates)),'abs_delta').index)
    for i,row in s.iterrows():
        if not all_cells and i not in selected:continue
        key=tuple(row[k] for k in CELL)
        mask=np.ones(len(m),dtype=bool)
        for k,v in zip(CELL,key):mask &= np.isclose(m[k],v,rtol=0,atol=1e-8)
        d=m.loc[mask].sort_values('phi_slot')
        if d.empty:continue
        cell_plot(d,detail/f'cell_{i+1:03d}.png',f'Cell {i+1:03d}: {cell_label(row)}')
    # Contact sheet: the most informative selected cells, not 90 disconnected images.
    chosen=s.loc[sorted(selected)]
    if chosen.empty:return
    ncols=3;nrows=int(np.ceil(len(chosen)/ncols))
    fig,axes=plt.subplots(nrows,ncols,figsize=(15,3.2*nrows),squeeze=False)
    for ax in axes.ravel():ax.set_visible(False)
    for ax,(_,row) in zip(axes.ravel(),chosen.iterrows()):
        ax.set_visible(True);mask=np.ones(len(m),dtype=bool)
        for k in CELL:mask &= np.isclose(m[k],row[k],rtol=0,atol=1e-8)
        d=m.loc[mask]
        for prefix,col,marker in [('pub',COLORS['pub'],'o'),('our',COLORS['our'],'s')]:
            ok=valid(d,prefix)
            ax.errorbar(d.loc[ok,f'{prefix}_phi'],d.loc[ok,f'{prefix}_bsa'],yerr=d.loc[ok,f'{prefix}_stat'],
                        fmt=marker,ms=2.8,color=col,alpha=.85,elinewidth=.7)
        ax.axhline(0,color='0.5',lw=.6);ax.set(xlim=(0,360),title=f"{cell_label(row)}\nΔa₁={row.delta_a1_common:+.3f}; common N={row.n_common}")
        ax.tick_params(labelsize=8);ax.title.set_fontsize(8)
    fig.suptitle('Selected Fa18 BSA comparisons: blue = Sebastian; red = ours',fontsize=14)
    fig.tight_layout(rect=[0,0,1,.98]);save(fig,out/'05_selected_cell_overlays.png')

def apply_flags(m, exclude_cells, exclude_points, out):
    """Manual exclusions affect amplitude fits/averages, not displayed raw data."""
    m=m.copy()
    if 'cell_id' not in m:
        m['cell_id']=m.groupby(CELL,sort=True).ngroup()+1
    m['cell_id']=pd.to_numeric(m['cell_id'],errors='coerce')
    if exclude_cells:
        flags=pd.read_csv(exclude_cells)
        if 'cell_id' not in flags:raise ValueError('Exclude-cells CSV needs cell_id column')
        if 'exclude' in flags.columns:
            flags=flags.loc[pd.to_numeric(flags['exclude'],errors='coerce').fillna(0)!=0]
        ids=set(pd.to_numeric(flags.cell_id,errors='raise').astype(int))
        m.loc[m.cell_id.isin(ids),'flagged']=True
    if exclude_points:
        flags=pd.read_csv(exclude_points)
        if not {'cell_id','phi_slot'}.issubset(flags):raise ValueError('Exclude-points CSV needs cell_id,phi_slot')
        pairs=set(zip(pd.to_numeric(flags.cell_id).astype(int),pd.to_numeric(flags.phi_slot).astype(int)))
        m.loc[[((int(c),int(p)) in pairs) for c,p in zip(m.cell_id,m.phi_slot)],'flagged']=True
    template=m[['cell_id']+CELL].drop_duplicates().sort_values('cell_id').copy()
    template['exclude']=0
    template['reason']=''
    template.to_csv(out/'cell_flag_template.csv',index=False)
    m[['cell_id','phi_slot','pub_phi','our_phi','delta_phi','common_valid','flagged']].to_csv(out/'phi_point_audit.csv',index=False)
    return m

def xb_canvases(m,s,out):
    """Multi-panel xB canvases, ordered by Q2 and |t|; raw points plus common-bin fits."""
    directory=out/'xB_canvases';directory.mkdir(exist_ok=True)
    xb_groups=s.groupby(['xBmin','xBmax'],sort=True)
    for (xblo,xbhi),group in xb_groups:
        group=group.sort_values(['Q2min','t_abs_min'])
        ncols=4
        # Limit each page to 12 subplots to retain legibility.
        for page,start in enumerate(range(0,len(group),12),1):
            subset=group.iloc[start:start+12]
            nrows=int(np.ceil(len(subset)/ncols))
            fig,axes=plt.subplots(nrows,ncols,figsize=(16,3.45*nrows),squeeze=False,sharex=True,sharey=True)
            for ax in axes.flat:ax.set_visible(False)
            for ax,(_,r) in zip(axes.flat,subset.iterrows()):
                ax.set_visible(True)
                mask=np.ones(len(m),dtype=bool)
                for k in CELL:mask &= np.isclose(m[k],r[k],rtol=0,atol=1e-8)
                d=m.loc[mask].sort_values('phi_slot')
                for prefix,col,marker in [('pub',COLORS['pub'],'o'),('our',COLORS['our'],'s')]:
                    ok=valid(d,prefix)
                    ax.errorbar(d.loc[ok,f'{prefix}_phi'],d.loc[ok,f'{prefix}_bsa'],
                                yerr=d.loc[ok,f'{prefix}_stat'],fmt=marker,color=col,
                                ms=3.1,elinewidth=.8,capsize=1.5,ls='none',alpha=.9)
                    fitresult=fit(d,prefix,d.common_valid & ~d.flagged)
                    if np.isfinite(fitresult['a1']):
                        xx=np.linspace(0,360,361)
                        ax.plot(xx,model(xx,fitresult['a1'],fitresult['a2']),color=col,lw=1.3)
                if d.flagged.any():
                    ax.text(.98,.97,'FLAGGED',transform=ax.transAxes,ha='right',va='top',fontsize=7,color='darkred')
                ax.axhline(0,color='0.5',lw=.6,ls='--')
                ax.set(xlim=(0,360),ylim=(-1,1),xticks=[0,90,180,270,360])
                ax.set_title(f"Q² [{r.Q2min:g},{r.Q2max:g}]   |t| [{r.t_abs_min:g},{r.t_abs_max:g}]\n"
                             f"cell {int(r.cell_id):03d}  common {int(r.n_common)}/{int(r.n_published)}",
                             fontsize=9)
                ax.tick_params(labelsize=8);ax.grid(alpha=.16)
            for ax in axes[-1,:]:
                if ax.get_visible():ax.set_xlabel(r'$\phi$ (degrees)',fontsize=9)
            for ax in axes[:,0]:
                if ax.get_visible():ax.set_ylabel(r'$A_{LU}$',fontsize=10)
            from matplotlib.lines import Line2D
            handles=[Line2D([0],[0],color=COLORS['pub'],marker='o',lw=1.5,ms=5,label='Sebastian (Table II)'),
                     Line2D([0],[0],color=COLORS['our'],marker='s',lw=1.5,ms=5,label='Ours (Fa18 combined)')]
            fig.legend(handles=handles,loc='upper center',ncol=2,bbox_to_anchor=(.5,.972),fontsize=11)
            fig.suptitle(f'Fa18 Pass-2 BSA comparison — xB in [{xblo:g}, {xbhi:g}]  (page {page})',fontsize=15,y=.998)
            fig.tight_layout(rect=[0,0,1,.925])
            save(fig,directory/f'xB_{xblo:.3f}_{xbhi:.3f}_page{page:02d}.png')

def phi_diagnostics(m,s,out,individual=False):
    directory=out/'phi_diagnostics';directory.mkdir(exist_ok=True)
    d=m.loc[m.common_valid].copy()
    fig,axes=plt.subplots(1,2,figsize=(12.5,4.6))
    axes[0].hist(d.delta_phi.dropna(),bins=55,color='0.45')
    axes[0].set(xlabel='Our measured φ − Sebastian measured φ (degrees)',ylabel='Number of common points',title='Difference in event-weighted φ means')
    axes[1].scatter(d.delta_phi,d.delta,s=12,alpha=.5,color=COLORS['our'])
    axes[1].axhline(0,color='0.4',ls='--',lw=.8)
    axes[1].set(xlabel='Our measured φ − Sebastian measured φ (degrees)',ylabel='Our BSA − Sebastian BSA',title='BSA differences versus φ-mean differences')
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();save(fig,directory/'phi_mean_differences.png')
    rows=[]
    for cell,g in m.groupby(CELL,sort=True):
        good=g.loc[g.common_valid]
        row=dict(zip(CELL,cell));row['cell_id']=g.cell_id.iloc[0]
        row.update(n_common=len(good),n_published=len(g),n_our_valid=int(g.our_valid.sum()),
                   max_abs_delta_phi=np.nanmax(np.abs(good.delta_phi)) if len(good) else np.nan,
                   mean_abs_delta_phi=np.nanmean(np.abs(good.delta_phi)) if len(good) else np.nan,
                   rms_bsa_difference=np.sqrt(np.mean(good.delta**2)) if len(good) else np.nan)
        rows.append(row)
    pd.DataFrame(rows).to_csv(directory/'phi_coverage_by_cell.csv',index=False)
    if not individual:return
    for _,row in s.iterrows():
        mask=np.ones(len(m),dtype=bool)
        for k in CELL:mask &= np.isclose(m[k],row[k],rtol=0,atol=1e-8)
        g=m.loc[mask].sort_values('phi_slot')
        fig,(ax,ax2)=plt.subplots(2,1,figsize=(8,5.5),sharex=False,gridspec_kw={'height_ratios':[2,1]})
        for prefix,col,mark,label in [('pub',COLORS['pub'],'o','Sebastian'),('our',COLORS['our'],'s','Ours')]:
            ok=valid(g,prefix)
            ax.errorbar(g.loc[ok,f'{prefix}_phi'],g.loc[ok,f'{prefix}_bsa'],yerr=g.loc[ok,f'{prefix}_stat'],
                        fmt=mark,color=col,ms=3.5,capsize=2,label=label)
        ok=g.common_valid
        ax2.scatter(g.loc[ok,'phi_slot'],g.loc[ok,'delta_phi'],s=25,color='0.3')
        ax2.axhline(0,color='0.4',lw=.7,ls='--')
        for slot in g.loc[g.flagged,'phi_slot']:
            ax.axvspan(slot*22.5,(slot+1)*22.5,color='0.85',alpha=.5)
        ax.set(ylabel=r'$A_{LU}$',title=f"Cell {int(row.cell_id):03d}: {cell_label(row)}")
        ax.legend(fontsize=9);ax.grid(alpha=.15)
        ax2.set(xlabel='φ bin slot (0–15)',ylabel='Our − Sebastian\nmean φ (deg)',xlim=(-.5,15.5),xticks=[0,3,6,9,12,15])
        ax2.grid(alpha=.2)
        fig.tight_layout();save(fig,directory/f'phi_cell_{int(row.cell_id):03d}.png')

def filtered_amplitude_summary(s,out):
    """Simple unweighted summaries only; covariance between datasets is unknown."""
    lines=['COMMON-BIN FIT SUMMARY (descriptive, not a covariance-aware consistency test)',
           'The two analyses share events, so errors on their differences are NOT known.',
           'Cells excluded using flags are omitted from this summary.']
    for p in ['a1','a2']:
        ok=np.isfinite(s[f'delta_{p}_common']) & (s.n_flagged==0)
        vals=s.loc[ok,f'delta_{p}_common'].to_numpy()
        lines += [f'{p}: {len(vals)} unflagged jointly fitted cells; mean(ours−Sebastian)={np.mean(vals):+.6f}; median={np.median(vals):+.6f}' if len(vals) else f'{p}: no unflagged jointly fitted cells']
    (out/'amplitude_summary.txt').write_text('\n'.join(lines)+'\n')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--our',default='output/csvs/dvcs_pass2_sebastian_binning.csv')
    ap.add_argument('--reference',default='imports/sebastian_tableII_proton_detected.csv')
    ap.add_argument('--output',default='output/sebastian_bsa_cross_check/direct_comparison')
    ap.add_argument('--all-cells',action='store_true',help='Retained for compatibility: all individual cell plots are now always produced')
    ap.add_argument('--exclude-cells',default='',help='Optional CSV listing cell_id values to exclude from global amplitude averages')
    ap.add_argument('--phi-cell-plots',action='store_true',help='Generate additional per-cell phi-mean diagnostic figures')
    ap.add_argument('--exclude-points',default='',help='Optional CSV listing cell_id,phi_slot pairs to exclude from common-bin fits')
    ap.add_argument('--top-cells',type=int,default=12,help='Number of representative cell overlays (default: 12)')
    args=ap.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    m,raw,ref=load(args.our,args.reference)
    m=apply_flags(m,args.exclude_cells,args.exclude_points,out)
    m.to_csv(out/'point_by_point.csv',index=False)
    overview(m,out);s=amplitudes(m,out);miss=missing(m,s,out)
    diagnostics(m,s,out,True,args.top_cells)
    xb_canvases(m,s,out)
    phi_diagnostics(m,s,out,args.phi_cell_plots)
    filtered_amplitude_summary(s,out)
    report=[f'Published Table II rows: {len(ref)}',f'Kinematic cells: {len(s)}',
            f'Rows matched to analysis bin keys: {(m._merge=="both").sum()}',
            f'Common valid BSA points: {m.common_valid.sum()}',
            f'Published rows without common valid BSA: {len(miss)}',
            f'Cells with both common-bin a1 fits: {np.isfinite(s.delta_a1_common).sum()}',
            f'Flagged points: {int(m.flagged.sum())}',
            'xB canvases: xB_canvases/; detailed plots: individual_cells/; phi checks: phi_diagnostics/',
            'Amplitudes fitted to common valid phi slots using each dataset’s own measured phiavg.',
            'Difference/error ratios are NOT calibrated significance: datasets share experimental events.',
            'Table II includes some anomalous entries; see the reference CSV and invalid-point listing.',
            f'Output: {out}']
    (out/'README_results.txt').write_text('\n'.join(report)+'\n')
    print('\n'.join(report))

if __name__=='__main__':main()
