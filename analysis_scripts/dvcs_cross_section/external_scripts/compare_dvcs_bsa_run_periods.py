#!/usr/bin/env python3
"""Independent-run CLAS12 DVCS BSA comparison in the regular 4D binning.

Default input: ../output/csvs/dvcs_pass2_analysis.csv (run from external_scripts).
Requires the measured phi means and 'BSA, counts, <period>' tuple columns.
"""
import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2, norm
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D

KEY=['Q2min','Q2max','xBmin','xBmax','t_abs_min','t_abs_max']
PERIODS={'Fa18':('Fa18',10.6),'Sp18':('Sp18',10.6),'Sp19':('Sp19 Inb',10.2)}
COLORS={'Fa18':'#d62728','Sp18':'#1764ad','Sp19':'#228b55'}

def modulation(phi,a1,a2):
    theta=np.deg2rad(phi)
    return a1*np.sin(theta)/(1+a2*np.cos(theta))

def find_col(columns, candidates):
    for name in candidates:
        if name in columns:return name
    return None

def parse_pair(v):
    if pd.isna(v):return (np.nan,np.nan)
    if isinstance(v,(float,int)):return (float(v),np.nan)
    nums=re.findall(r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?',str(v))
    return tuple(map(float,nums[:2])) if len(nums)>=2 else (np.nan,np.nan)

def load(path):
    df=pd.read_csv(path,low_memory=False)
    for k in KEY:
        if k not in df:raise KeyError(f'Missing bin key {k}; available: {list(df.columns)}')
        df[k]=pd.to_numeric(df[k],errors='raise')
    if not {'phimin','phimax'}.issubset(df):raise KeyError('Expected phimin and phimax')
    # The analysis uses variable-width phi bins (often 15 or 30 degrees).
    # Dividing phimin by 22.5 incorrectly aliases distinct bins.
    df['phimin']=pd.to_numeric(df['phimin'],errors='raise')
    df['phimax']=pd.to_numeric(df['phimax'],errors='raise')
    if df.duplicated(KEY+['phimin','phimax']).any():
        raise ValueError('Duplicate rows with identical 4D bin edges')
    df=df.sort_values(KEY+['phimin','phimax']).copy()
    df['phi_slot']=df.groupby(KEY,sort=False).cumcount()
    for period in PERIODS:
        source=PERIODS[period][0]
        bcol=find_col(df.columns,[f'BSA, counts, {source}',f'BSA, {source}'])
        pcol=find_col(df.columns,[f'phiavg, {source}'])
        if bcol is None:raise KeyError(f'Missing BSA column for {period}; available: {list(df.columns)}')
        if pcol is None:raise KeyError(f'Missing measured phiavg for {period}; available: {list(df.columns)}')
        pair=np.array([parse_pair(v) for v in df[bcol]],dtype=float)
        df[f'{period}_bsa']=pair[:,0];df[f'{period}_stat']=pair[:,1]
        df[f'{period}_phi']=pd.to_numeric(df[pcol],errors='coerce')
    return df

def valid(df,period):
    return (np.isfinite(df[f'{period}_bsa']) & np.isfinite(df[f'{period}_stat']) &
            (df[f'{period}_stat']>0) & np.isfinite(df[f'{period}_phi']))

def fit(d,period,mask):
    q=d.loc[mask & valid(d,period)]
    if len(q)<5:return None
    try:
        pars,cov=curve_fit(modulation,q[f'{period}_phi'],q[f'{period}_bsa'],
            sigma=q[f'{period}_stat'],absolute_sigma=True,p0=(.2,0),
            bounds=([-2,-.95],[2,.95]),maxfev=30000)
        chisq=np.sum(((q[f'{period}_bsa']-modulation(q[f'{period}_phi'],*pars))/q[f'{period}_stat'])**2)
        return dict(a1=pars[0],a2=pars[1],a1_err=np.sqrt(cov[0,0]),a2_err=np.sqrt(cov[1,1]),
                    chi2=chisq,ndf=len(q)-2,p=chi2.sf(chisq,len(q)-2))
    except (RuntimeError,ValueError):return None

def analyze(df,other,out):
    out.mkdir(parents=True,exist_ok=True)
    data=[];cells=[];fits=[]
    groups=list(df.groupby(KEY,sort=True,dropna=False))
    for idx,(key,d) in enumerate(groups):
        d=d.sort_values('phi_slot');mask=valid(d,'Fa18') & valid(d,other)
        fit_results={p:fit(d,p,mask) for p in ('Fa18',other)}
        for p,r in fit_results.items():
            if r:fits.append(dict(cell_id=idx,period=p,**dict(zip(KEY,key)),**r))
        for _,r in d.loc[mask].iterrows():
            delta=r[f'{other}_bsa']-r['Fa18_bsa']
            sigma=np.hypot(r[f'{other}_stat'],r['Fa18_stat'])
            data.append(dict(cell_id=idx,**dict(zip(KEY,key)),phi_slot=r.phi_slot,
                             phi_fa18=r['Fa18_phi'],phi_other=r[f'{other}_phi'],
                             fa18_bsa=r['Fa18_bsa'],other_bsa=r[f'{other}_bsa'],
                             fa18_stat=r['Fa18_stat'],other_stat=r[f'{other}_stat'],
                             delta=delta,sigma_delta=sigma,pull=delta/sigma))
        if mask.sum():
            z=np.array([(r[f'{other}_bsa']-r['Fa18_bsa'])/
                        np.hypot(r[f'{other}_stat'],r['Fa18_stat']) for _,r in d.loc[mask].iterrows()])
            c=float(np.dot(z,z));n=len(z)
            cells.append(dict(cell_id=idx,**dict(zip(KEY,key)),n_common=n,
                              chi2=c,ndf=n,chi2_ndf=c/n,p_value=chi2.sf(c,n)))
    pts=pd.DataFrame(data);cs=pd.DataFrame(cells);ft=pd.DataFrame(fits)
    if pts.empty:raise ValueError(f'No common valid Fa18/{other} points. Check BSA and phiavg columns.')
    pts.to_csv(out/'pointwise_comparison.csv',index=False)
    cs.to_csv(out/'cell_compatibility.csv',index=False)
    ft.to_csv(out/'modulation_fits.csv',index=False)
    cglobal=float(np.dot(pts.pull,pts.pull));n=len(pts)
    with open(out/'summary.txt','w') as f:
        f.write(f'Fa18 (10.6 GeV) versus {other} ({PERIODS[other][1]:g} GeV)\n')
        f.write(f'Common points: {n}; cells: {len(cs)}\n')
        f.write(f'Global pointwise chi2/ndf: {cglobal:.5g}/{n} = {cglobal/n:.5g}; p = {chi2.sf(cglobal,n):.5g}\n')
        f.write(f'Pull mean: {pts.pull.mean():.5g}; sample std: {pts.pull.std(ddof=1):.5g}\n')
        f.write('Only statistical uncertainties; independent periods; no energy correction applied.\n')
    # Two summary canvases, one per run comparison.
    fig,axes=plt.subplots(1,2,figsize=(12,4.7))
    z=pts.pull.to_numpy();edges=np.arange(min(-5,np.floor(z.min())-.5),max(5,np.ceil(z.max())+.5)+.25,.25)
    axes[0].hist(z,bins=edges,color=COLORS[other],alpha=.8,label='Measured pulls')
    xx=np.linspace(edges[0],edges[-1],600)
    axes[0].plot(xx,len(z)*np.diff(edges)[0]*norm.pdf(xx),color='black',ls='--',label='Standard normal')
    axes[0].set(xlabel=r'$(A_{LU}^{'+other+r'}-A_{LU}^{Fa18})/\sigma_{\Delta}$',ylabel='Number of points')
    axes[0].legend(fontsize=9)
    axes[0].text(.03,.97,f'N = {n}\nMean = {z.mean():+.3f}\nSD = {z.std(ddof=1):.3f}',
                 transform=axes[0].transAxes,va='top',fontsize=10)
    axes[1].hist(cs.chi2_ndf,bins=25,color=COLORS[other],alpha=.8)
    axes[1].axvline(1,color='black',ls='--')
    axes[1].set(xlabel=r'Cellwise $\chi^2/\mathrm{ndf}$',ylabel='Number of cells')
    fig.suptitle(f'Fa18 vs {other} BSA consistency | global χ²/ndf = {cglobal:.1f}/{n} = {cglobal/n:.2f}; p = {chi2.sf(cglobal,n):.2g}')
    fig.tight_layout(rect=(0,0,1,.92));fig.savefig(out/'summary_pulls.png',dpi=170);plt.close(fig)
    # Compact xB canvases, rows = |t|, columns = populated Q2 cells.
    canvas=out/'xB_canvases';canvas.mkdir(exist_ok=True)
    handles=[Line2D([0],[0],color=COLORS[p],marker=m,lw=1.4,label=f'{p} ({PERIODS[p][1]:g} GeV)')
             for p,m in [('Fa18','s'),(other,'o')]]
    for (xblo,xbhi),group in cs.groupby(['xBmin','xBmax'],sort=True):
        tbins=sorted(set(zip(group.t_abs_min,group.t_abs_max)))
        rows=[]
        for lo,hi in tbins:
            rows.append(group[np.isclose(group.t_abs_min,lo)&np.isclose(group.t_abs_max,hi)].sort_values(['Q2min','Q2max']))
        nr=len(rows);nc=max(map(len,rows))
        fig=plt.figure(figsize=(max(6.4,3.55*nc+1.1),max(3.8,2.85*nr+1.05)))
        grid=fig.add_gridspec(nr,nc,left=.105,right=.987,bottom=.095,top=.875,wspace=.13,hspace=.26)
        for i,(tbin,row) in enumerate(zip(tbins,rows)):
            for j,(_,r) in enumerate(row.iterrows()):
                sub=GridSpecFromSubplotSpec(2,1,subplot_spec=grid[i,j],height_ratios=[3.3,1],hspace=.035)
                ax=fig.add_subplot(sub[0]);axd=fig.add_subplot(sub[1],sharex=ax)
                d=df[np.logical_and.reduce([np.isclose(df[k],r[k],atol=1e-8,rtol=0) for k in KEY])].sort_values('phi_slot')
                mask=valid(d,'Fa18') & valid(d,other)
                q=d.loc[mask]
                for p,m in [('Fa18','s'),(other,'o')]:
                    ok=valid(d,p)
                    ax.errorbar(d.loc[ok,f'{p}_phi'],d.loc[ok,f'{p}_bsa'],yerr=d.loc[ok,f'{p}_stat'],
                                fmt=m,color=COLORS[p],ms=2.8,elinewidth=.75,capsize=1.2,ls='none',alpha=.9)
                    fr=fit(d,p,mask)
                    if fr:
                        xx=np.linspace(0,360,361)
                        ax.plot(xx,np.ma.masked_outside(modulation(xx,fr['a1'],fr['a2']),-.5,.5),color=COLORS[p],lw=1.2)
                difference=q[f'{other}_bsa'].to_numpy()-q['Fa18_bsa'].to_numpy()
                err=np.hypot(q[f'{other}_stat'],q['Fa18_stat'])
                axd.errorbar(q['Fa18_phi'],difference,yerr=err,fmt='o',color='0.3',ecolor='0.65',ms=2.5,capsize=1,elinewidth=.7)
                for a in (ax,axd):a.axhline(0,color='0.5',lw=.6,ls='--');a.grid(alpha=.15)
                ax.set(xlim=(0,360),ylim=(-.5,.5),yticks=[-.5,-.25,0,.25,.5])
                spread=max(.15,float(np.max(np.abs(difference)+err))*1.1)
                axd.set(xlim=(0,360),ylim=(-spread,spread),xticks=[0,90,180,270,360])
                ax.set_title(rf'$Q^2\in[{r.Q2min:g},{r.Q2max:g}]$ GeV$^2$ · Cell {int(r.cell_id):03d} · {len(q)}/16',fontsize=9,pad=4)
                ax.text(.025,.035,rf'$\chi^2/\mathrm{{ndf}}={r.chi2:.1f}/{int(r.ndf)}={r.chi2_ndf:.2f}$'+'\n'+rf'$p={r.p_value:.2g}$',
                        transform=ax.transAxes,ha='left',va='bottom',fontsize=7.5,bbox=dict(facecolor='white',alpha=.83,edgecolor='none',pad=1.3))
                ax.tick_params(labelbottom=False,labelsize=7);axd.tick_params(labelsize=7)
                if j>0:ax.tick_params(labelleft=False);axd.tick_params(labelleft=False)
                if i<nr-1:axd.tick_params(labelbottom=False)
            ypos=.875-(i+.5)*(.78/nr)
            fig.text(.026,ypos,rf'$|t|\in[{tbin[0]:g},{tbin[1]:g}]$ GeV$^2$',rotation=90,ha='center',va='center',fontsize=10)
        fig.text(.53,.965,rf'BSA comparison: $x_B\in[{xblo:g},{xbhi:g}]$',ha='center',va='top',fontsize=14)
        fig.legend(handles=handles,loc='upper center',ncol=2,bbox_to_anchor=(.53,.944),fontsize=10,frameon=False)
        fig.supxlabel(r'$\phi$ (degrees)',y=.052,fontsize=10)
        fig.text(.53,.016,f'Independent statistical errors; {other} beam energy {PERIODS[other][1]:g} GeV; no energy correction.',ha='center',fontsize=7.5,color='0.35')
        fig.savefig(canvas/f'xB_{xblo:.3f}_{xbhi:.3f}.png',dpi=170,bbox_inches='tight');plt.close(fig)
    return n,len(cs),cglobal,chi2.sf(cglobal,n)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-csv',default='../output/csvs/dvcs_pass2_analysis.csv',help='Override the regular DVCS analysis CSV')
    ap.add_argument('--output-dir',default='../output/bsa_run_period_consistency')
    args=ap.parse_args()
    df=load(args.input_csv)
    for other in ('Sp18','Sp19'):
        n,ncells,c,p=analyze(df,other,Path(args.output_dir)/f'Fa18_vs_{other}')
        print(f'Fa18 vs {other}: {n} common points in {ncells} cells; chi2/ndf={c:.2f}/{n}={c/n:.3f}; p={p:.4g}')

if __name__=='__main__':main()
