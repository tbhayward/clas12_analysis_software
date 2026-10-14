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
from scipy.stats import chi2 as chi2_distribution
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
    common=m.loc[m.common_valid].copy()
    if common.empty:raise ValueError('No common finite BSA measurements')
    x=common.pub_bsa.to_numpy();y=common.our_bsa.to_numpy();delta=y-x
    r=np.corrcoef(x,y)[0,1] if len(x)>2 and np.std(x)>0 and np.std(y)>0 else np.nan
    slope,intercept=np.polyfit(x,y,1) if len(x)>2 and np.std(x)>0 else (np.nan,np.nan)
    fig,ax=plt.subplots(figsize=(7,6))
    ax.scatter(x,y,s=13,alpha=.55,color=COLORS['our'],label='Matched measurements')
    lim=max(.2,np.nanpercentile(np.abs(np.r_[x,y]),99)*1.15)
    xx=np.array([-lim,lim]);ax.plot(xx,xx,'k--',lw=1,label='Equality')
    if np.isfinite(slope):ax.plot(xx,slope*xx+intercept,color='0.35',ls=':',lw=1.6,label=f'OLS: y={slope:.3f}x{intercept:+.3f}')
    ax.text(.03,.97,f'N = {len(common)}\nPearson r = {r:.3f}\nMean Δ = {np.mean(delta):+.4f}\nRMS Δ = {np.sqrt(np.mean(delta**2)):.4f}',
            transform=ax.transAxes,va='top',fontsize=10,bbox=dict(facecolor='white',alpha=.85,edgecolor='0.8'))
    ax.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel='Alvarado BSA',ylabel='Hayward BSA',title='Point-by-point BSA comparison')
    ax.legend(loc='lower right');ax.grid(alpha=.2);save(fig,out/'01_bsa_scatter.png')
    records=[]
    for category,part in [('all',common),('complete',common.loc[common.cell_id.isin(m.groupby('cell_id').filter(lambda g: int(g.common_valid.sum())==16).cell_id.unique())]),('incomplete',common.loc[~common.cell_id.isin(m.groupby('cell_id').filter(lambda g: int(g.common_valid.sum())==16).cell_id.unique())])]:
        a=part.pub_bsa.to_numpy();b=part.our_bsa.to_numpy()
        records.append(dict(coverage=category,n=len(part),pearson_r=np.corrcoef(a,b)[0,1] if len(a)>2 and np.std(a)>0 and np.std(b)>0 else np.nan,
                            mean_delta=np.mean(b-a) if len(a) else np.nan,rms_delta=np.sqrt(np.mean((b-a)**2)) if len(a) else np.nan))
    pd.DataFrame(records).to_csv(out/'bsa_correlation_by_coverage.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(12.5,4.7))
    axes[0].hist(delta,bins=50,color=COLORS['our'],alpha=.85)
    axes[0].axvline(0,color='k',ls='--')
    axes[0].set(xlabel='Hayward − Alvarado BSA',ylabel='Number of points',title=f'Raw differences (mean {np.mean(delta):+.4f}, RMS {np.sqrt(np.mean(delta**2)):.4f})')
    slot_rows=[]
    for slot,g in common.groupby('phi_slot'):
        vals=g.delta.to_numpy();n=len(vals);sem=np.std(vals,ddof=1)/np.sqrt(n) if n>1 else np.nan
        slot_rows.append(dict(phi_slot=int(slot),phi_center=(int(slot)+.5)*22.5,n=n,mean_delta=np.mean(vals),sem_delta=sem,rms_delta=np.sqrt(np.mean(vals**2))))
    slot_stats=pd.DataFrame(slot_rows).sort_values('phi_slot');slot_stats.to_csv(out/'bsa_differences_by_phi_slot.csv',index=False)
    axes[1].errorbar(slot_stats.phi_center,slot_stats.mean_delta,yerr=slot_stats.sem_delta,fmt='o-',color=COLORS['our'],capsize=3)
    axes[1].axhline(0,color='k',ls='--')
    for row in slot_stats.itertuples():axes[1].annotate(f'{row.n}',(row.phi_center,row.mean_delta),xytext=(0,8),textcoords='offset points',ha='center',fontsize=7)
    axes[1].set(xlabel='Azimuthal bin center (degrees)',ylabel='Mean Hayward − Alvarado BSA',title='Mean difference by azimuth (labels: N)',xlim=(0,360),xticks=np.arange(0,361,45))
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();save(fig,out/'02_bsa_differences.png')
    z=common.independent_error_pull.to_numpy();z=z[np.isfinite(z)]
    fig,ax=plt.subplots(figsize=(7,4.7));ax.hist(z,bins=55,color=COLORS['our'],alpha=.85)
    ax.axvline(0,color='k',ls='--');ax.set(xlabel='ΔBSA / sqrt(stat(Hayward)² + stat(Alvarado)²)',ylabel='Number of points',
      title=f'Naively normalized differences (N={len(z)}, mean={np.mean(z):+.3f}, RMS={np.sqrt(np.mean(z*z)):.3f})')
    ax.text(.02,.98,'Shared events: unknown covariance; NOT calibrated pulls',transform=ax.transAxes,va='top',fontsize=9)
    ax.grid(alpha=.2);fig.tight_layout();save(fig,out/'02b_naive_normalized_differences.png')
    fig,ax=plt.subplots(figsize=(7,4.7));ax.scatter(common.pub_phi,common.delta,s=12,alpha=.5,color=COLORS['our']);ax.axhline(0,color='k',ls='--')
    ax.set(xlabel='Alvarado measured φ (degrees)',ylabel='Hayward − Alvarado BSA',title='Individual BSA differences versus azimuth',xlim=(0,360));ax.grid(alpha=.2)
    fig.tight_layout();save(fig,out/'02c_raw_azimuth_scatter.png')

def amplitudes(m,out):
    rows=[]
    for cell,d in m.groupby(CELL,sort=True):
        common=d.common_valid & ~d.flagged
        p_all=fit(d,'pub');o_all=fit(d,'our')
        p=fit(d,'pub',common);o=fit(d,'our',common)
        row=dict(zip(CELL,cell));row['cell_id']=int(d.cell_id.iloc[0]) if 'cell_id' in d and pd.notna(d.cell_id.iloc[0]) else np.nan
        row.update(coverage_class=('complete' if int(common.sum())==16 else 'severe' if int(common.sum())<=4 else 'partial'),n_published=len(d),n_published_valid=int(d.pub_valid.sum()),n_our_valid=int(d.our_valid.sum()),n_common=int(common.sum()),n_flagged=int(d.flagged.sum()),n_common_before_flags=int(d.common_valid.sum()))
        for tag,res in [('published_all',p_all),('our_all',o_all),('published_common',p),('our_common',o)]:
            for k,v in res.items():row[f'{tag}_{k}']=v
        for param in ['a1','a2']:
            row[f'delta_{param}_common']=o[param]-p[param]
            row[f'naive_delta_{param}_err']=np.hypot(o[f'{param}_err'],p[f'{param}_err'])
        rows.append(row)
    s=pd.DataFrame(rows);s.to_csv(out/'amplitude_comparison.csv',index=False)
    stats=[]
    for cls,g in s.groupby('coverage_class'):
        for p in ['a1','a2']:
            v=g[f'delta_{p}_common'].dropna()
            stats.append(dict(coverage_class=cls,parameter=p,n_fitted=len(v),mean_delta=v.mean() if len(v) else np.nan,
                              rms_delta=np.sqrt(np.mean(v**2)) if len(v) else np.nan))
    pd.DataFrame(stats).to_csv(out/'amplitude_differences_by_coverage.csv',index=False)

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
        axes[0].set(xlim=(lo-pad,hi+pad),ylim=(lo-pad,hi+pad),xlabel=f'Alvarado {param}',ylabel=f'Hayward {param}',title='Common-φ-bin fits')
        idx=np.arange(len(d))
        axes[1].errorbar(idx,d[f'delta_{param}_common'],yerr=d[f'naive_delta_{param}_err'],fmt='o',ms=3,color=COLORS['our'],ecolor='0.65',elinewidth=.7)
        axes[1].axhline(0,color='k',ls='--');axes[1].set(xlabel='Kinematic cell (sorted by Q², xB, |t|)',ylabel=f'Hayward {param} − Alvarado {param}',title='Amplitude differences (errors ignore covariance)')
        axes[1].set_ylim(*((-0.2,0.2) if param=='a1' else (-0.3,0.3)))
        axes[1].text(.99,.98,'Points outside axis range are clipped',transform=axes[1].transAxes,ha='right',va='top',fontsize=8)
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
    # Full 0–360 degree coverage map: every published cell and all 16 azimuthal slots.
    cells=s[['cell_id']+CELL].drop_duplicates().sort_values('cell_id')
    maps=[]
    for row in cells.itertuples(index=False):
        g=m.loc[m.cell_id==row.cell_id].set_index('phi_slot')
        for slot in range(16):
            present=slot in g.index
            rec=g.loc[slot] if present else None
            maps.append(dict(cell_id=int(row.cell_id),phi_slot=slot,
                alvarado_valid=bool(rec.pub_valid) if present else False,
                hayward_valid=bool(rec.our_valid) if present else False,
                common_valid=bool(rec.common_valid) if present else False,
                excluded=bool(rec.flagged) if present else False))
    coverage=pd.DataFrame(maps);coverage.to_csv(out/'phi_coverage_map.csv',index=False)
    for field,filename,title in [('alvarado_valid','04b_alvarado_phi_coverage.png','Alvarado'),('hayward_valid','04c_hayward_phi_coverage.png','Hayward'),('common_valid','04d_common_phi_coverage.png','Common')]:
        grid=coverage.pivot(index='cell_id',columns='phi_slot',values=field).fillna(False).astype(int)
        fig,ax=plt.subplots(figsize=(11,max(6,len(grid)*.13)))
        ax.imshow(grid.to_numpy(),aspect='auto',interpolation='nearest',cmap='Greys',vmin=0,vmax=1)
        ax.set(xlabel='Azimuthal slot (0–15)',ylabel='Cell ID',title=f'{title} valid azimuthal coverage (dark = valid)',xticks=range(16))
        ax.set_yticks(np.arange(0,len(grid),max(1,len(grid)//20)))
        ax.set_yticklabels(grid.index[::max(1,len(grid)//20)])
        fig.tight_layout();save(fig,out/filename)
    miss=m.loc[~m.common_valid].copy();miss.to_csv(out/'missing_or_invalid_points.csv',index=False)
    return miss

def cell_plot(d,path,cell_title):
    fig,(ax,axd)=plt.subplots(2,1,figsize=(7.5,6.2),sharex=True,gridspec_kw={'height_ratios':[2.5,1]})
    for prefix,label,col,mark in [('pub','Alvarado',COLORS['pub'],'o'),('our','Hayward',COLORS['our'],'s')]:
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
    axd.set(xlabel=r'Measured $\phi$ (degrees)',ylabel='Hayward − Alvarado',xlim=(0,360));axd.grid(alpha=.2)
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
    fig.suptitle('Selected Fa18 BSA comparisons: blue = Alvarado; red = Hayward',fontsize=14)
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
    pd.DataFrame(columns=['cell_id','phi_slot','reason']).to_csv(out/'point_flag_template.csv',index=False)
    m[['cell_id','phi_slot','pub_phi','our_phi','delta_phi','common_valid','flagged']].to_csv(out/'phi_point_audit.csv',index=False)
    return m

def xb_canvases(m,s,out):
    """One canvas per xB interval: Q2 columns and |t| rows.

    The lower strip shows Hayward minus Alvarado with quadrature statistical
    uncertainties. Compatibility uses common measured slots, no curve fit.
    """
    from matplotlib.lines import Line2D
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    directory=out/'xB_canvases'
    directory.mkdir(exist_ok=True)
    for old in directory.glob('xB_*_page*.png'):
        old.unlink()
    handles=[Line2D([0],[0],color=COLORS['pub'],marker='o',lw=1.5,ms=5,label='Alvarado'),
             Line2D([0],[0],color=COLORS['our'],marker='s',lw=1.5,ms=5,label='Hayward')]
    compatibility=[]
    for (xblo,xbhi),group in s.groupby(['xBmin','xBmax'],sort=True):
        group=group.loc[group.n_common>0].copy()
        if group.empty:continue
        qbins=sorted(set(zip(group.Q2min,group.Q2max)))
        tbins=sorted(set(zip(group.t_abs_min,group.t_abs_max)))
        ncols,nrows=len(qbins),len(tbins)
        fig=plt.figure(figsize=(max(6,4.4*ncols),max(4.6,3.65*nrows)))
        grid=fig.add_gridspec(nrows,ncols,left=.105,right=.975,bottom=.085,top=.86,
                             wspace=.25,hspace=.52)
        slots={}
        for _,r in group.iterrows():
            i=tbins.index((r.t_abs_min,r.t_abs_max))
            j=qbins.index((r.Q2min,r.Q2max))
            sub=GridSpecFromSubplotSpec(2,1,subplot_spec=grid[i,j],height_ratios=[3,1],hspace=.05)
            ax=fig.add_subplot(sub[0]);axd=fig.add_subplot(sub[1],sharex=ax)
            slots[(i,j)]=(ax,axd)
            mask=np.ones(len(m),dtype=bool)
            for k in CELL:mask &= np.isclose(m[k],r[k],rtol=0,atol=1e-8)
            d=m.loc[mask].sort_values('phi_slot')
            common=d.loc[d.common_valid & ~d.flagged].copy()
            n=len(common)
            if not n:
                plt.delaxes(ax);plt.delaxes(axd)
                del slots[(i,j)]
                continue
            variance=common.pub_stat.to_numpy()**2+common.our_stat.to_numpy()**2
            difference=common.our_bsa.to_numpy()-common.pub_bsa.to_numpy()
            chisq=float(np.sum(difference**2/variance))
            pvalue=float(chi2_distribution.sf(chisq,n))
            compatibility.append(dict(cell_id=int(r.cell_id),xBmin=xblo,xBmax=xbhi,
                Q2min=r.Q2min,Q2max=r.Q2max,t_abs_min=r.t_abs_min,t_abs_max=r.t_abs_max,
                n_common=n,chi2=chisq,ndf=n,chi2_ndf=chisq/n,p_value_naive=pvalue))
            # Use real newlines and single backslashes in matplotlib mathtext.
            annotation=(rf'$\chi^2/\mathrm{{ndf}}={chisq:.1f}/{n}$'+'\n'+
                        rf'$p_{{\mathrm{{naive}}}}={pvalue:.2g}$')
            for prefix,col,marker in [('pub',COLORS['pub'],'o'),('our',COLORS['our'],'s')]:
                ok=valid(d,prefix)
                ax.errorbar(d.loc[ok,f'{prefix}_phi'],d.loc[ok,f'{prefix}_bsa'],
                    yerr=d.loc[ok,f'{prefix}_stat'],fmt=marker,color=col,
                    ms=3,elinewidth=.8,capsize=1.3,ls='none',alpha=.9)
                result=fit(d,prefix,d.common_valid & ~d.flagged)
                if np.isfinite(result['a1']):
                    xx=np.linspace(0,360,361)
                    ax.plot(xx,model(xx,result['a1'],result['a2']),color=col,lw=1.2)
            axd.errorbar(common.pub_phi,difference,yerr=np.sqrt(variance),
                fmt='o',color='0.3',ecolor='0.65',ms=2.6,capsize=1,elinewidth=.7)
            ax.axhline(0,color='0.5',lw=.6,ls='--')
            axd.axhline(0,color='0.5',lw=.6,ls='--')
            ax.set(xlim=(0,360),ylim=(-1,1))
            axd.set(xlim=(0,360),xticks=[0,90,180,270,360])
            # Keep the difference band compact, but avoid hiding outliers.
            spread=max(.15,float(np.nanmax(np.abs(difference)+np.sqrt(variance)))*1.1)
            axd.set_ylim(-spread,spread)
            ax.set_title(f'Cell {int(r.cell_id):03d} · common {n}/16',fontsize=9,pad=4)
            ax.text(.03,.04,annotation,transform=ax.transAxes,ha='left',va='bottom',
                fontsize=7.5,bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=1.5))
            ax.tick_params(labelbottom=False,labelsize=7)
            axd.tick_params(labelsize=7)
            ax.grid(alpha=.15);axd.grid(alpha=.15)
            axd.set_xlabel(r'$\phi$ (degrees)',fontsize=8)
            if j==0:
                ax.set_ylabel(r'$A_{LU}$',fontsize=9)
                axd.set_ylabel(r'$\Delta A$',fontsize=8)
        # Put column and row headings in figure coordinates, clear of data panels.
        for j,(lo,hi) in enumerate(qbins):
            x=.105+(j+.5)*(.87/ncols)
            fig.text(x,.89,rf'$Q^2\in[{lo:g},{hi:g}]$ GeV$^2$',ha='center',fontsize=10)
        for i,(lo,hi) in enumerate(tbins):
            y=.86-(i+.5)*(.775/nrows)
            fig.text(.025,y,rf'$|t|\in[{lo:g},{hi:g}]$ GeV$^2$',rotation=90,
                     ha='center',va='center',fontsize=9)
        fig.suptitle(rf'BSA comparison: $x_B\in[{xblo:g},{xbhi:g}]$',fontsize=14,y=.995)
        fig.legend(handles=handles,loc='upper center',ncol=2,bbox_to_anchor=(.5,.965),
                   fontsize=10,frameon=False)
        fig.text(.5,.018,'Compatibility p values assume independent statistical errors; shared-event covariance unknown.',
                 ha='center',fontsize=8,color='0.35')
        save(fig,directory/f'xB_{xblo:.3f}_{xbhi:.3f}.png')
    pd.DataFrame(compatibility,columns=['cell_id','xBmin','xBmax','Q2min','Q2max',
        't_abs_min','t_abs_max','n_common','chi2','ndf','chi2_ndf','p_value_naive']
        ).to_csv(directory/'pointwise_compatibility.csv',index=False)

def phi_diagnostics(m,s,out,individual=False):
    directory=out/'phi_diagnostics';directory.mkdir(exist_ok=True)
    d=m.loc[m.common_valid].copy()
    fig,axes=plt.subplots(1,2,figsize=(12.5,4.6))
    axes[0].hist(d.delta_phi.dropna(),bins=55,color='0.45')
    axes[0].set(xlabel='Hayward measured φ − Alvarado measured φ (degrees)',ylabel='Number of common points',title='Difference in event-weighted φ means')
    axes[1].scatter(d.delta_phi,d.delta,s=12,alpha=.5,color=COLORS['our'])
    axes[1].axhline(0,color='0.4',ls='--',lw=.8)
    axes[1].set(xlabel='Hayward measured φ − Alvarado measured φ (degrees)',ylabel='Hayward BSA − Alvarado BSA',title='BSA differences versus φ-mean differences')
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
        for prefix,col,mark,label in [('pub',COLORS['pub'],'o','Alvarado'),('our',COLORS['our'],'s','Hayward')]:
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
        ax2.set(xlabel='φ bin slot (0–15)',ylabel='Hayward − Alvarado\nmean φ (deg)',xlim=(-.5,15.5),xticks=[0,3,6,9,12,15])
        ax2.grid(alpha=.2)
        fig.tight_layout();save(fig,directory/f'phi_cell_{int(row.cell_id):03d}.png')

def filtered_amplitude_summary(s,out):
    """Simple unweighted summaries only; covariance between datasets is unknown."""
    lines=['COMMON-BIN FIT SUMMARY (descriptive, not a covariance-aware consistency test)',
           'The two analyses share events, so errors on their differences are NOT known.',
           'Flagged points are excluded from fits; fully excluded cells have no accepted fit.']
    for p in ['a1','a2']:
        ok=np.isfinite(s[f'delta_{p}_common'])
        vals=s.loc[ok,f'delta_{p}_common'].to_numpy()
        lines += [f'{p}: {len(vals)} unflagged jointly fitted cells; mean(Hayward−Alvarado)={np.mean(vals):+.6f}; median={np.median(vals):+.6f}' if len(vals) else f'{p}: no unflagged jointly fitted cells']
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
    ap.add_argument('--min-t-bin-edge',type=float,default=0.12,
                    help='Minimum lower edge of |t| bin in GeV^2 (default: 0.12; use 0 to restore all bins)')
    args=ap.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    m,raw,ref=load(args.our,args.reference)
    if not np.isfinite(args.min_t_bin_edge) or args.min_t_bin_edge < 0:
        ap.error('--min-t-bin-edge must be a finite nonnegative value')
    total_rows, total_cells = len(m), m[CELL].drop_duplicates().shape[0]
    selected = m['t_abs_min'] >= args.min_t_bin_edge - 1e-9
    m = m.loc[selected].copy()
    if m.empty:
        ap.error(f'No reference rows with |t| lower bin edge >= {args.min_t_bin_edge:g} GeV^2')
    ref = ref.loc[ref['t_abs_min'] >= args.min_t_bin_edge - 1e-9].copy()
    m=apply_flags(m,args.exclude_cells,args.exclude_points,out)
    m.to_csv(out/'point_by_point.csv',index=False)
    overview(m,out);s=amplitudes(m,out);miss=missing(m,s,out)
    uncertainty_rows=[]
    for prefix,label in [('our','Hayward'),('pub','Alvarado')]:
        d=m.loc[valid(m,prefix),f'{prefix}_stat'].to_numpy()
        uncertainty_rows.append(dict(analysis=label,n=len(d),median_stat=np.median(d),mean_stat=np.mean(d),p16_stat=np.percentile(d,16),p84_stat=np.percentile(d,84)))
    pd.DataFrame(uncertainty_rows).to_csv(out/'statistical_uncertainty_summary.csv',index=False)

    diagnostics(m,s,out,True,args.top_cells)
    xb_canvases(m,s,out)
    phi_diagnostics(m,s,out,args.phi_cell_plots)
    filtered_amplitude_summary(s,out)
    report=[f'Published Table II rows after |t| selection: {len(ref)} (of {total_rows})',
            f'Kinematic cells after |t| selection: {len(s)} (of {total_cells})',
            f'Applied |t| lower-bin-edge threshold: {args.min_t_bin_edge:g} GeV^2',
            f'Rows matched to analysis bin keys: {(m._merge=="both").sum()}',
            f'Common valid BSA points: {m.common_valid.sum()}',
            f'Published rows without common valid BSA: {len(miss)}',
            f'Cells with both common-bin a1 fits: {np.isfinite(s.delta_a1_common).sum()}',
            f'Flagged points: {int(m.flagged.sum())}',
            'xB canvases: xB_canvases/; detailed plots: individual_cells/; phi checks: phi_diagnostics/',
            'Amplitudes fitted to common valid phi slots using each dataset’s own measured phiavg.',
            'Difference/error ratios are NOT calibrated significance: datasets share experimental events.',
            'Reference data are from Table II (proton detected); statistical uncertainties only for both analyses.',
            f'Output: {out}']
    (out/'README_results.txt').write_text('\n'.join(report)+'\n')
    print('\n'.join(report))

if __name__=='__main__':main()
