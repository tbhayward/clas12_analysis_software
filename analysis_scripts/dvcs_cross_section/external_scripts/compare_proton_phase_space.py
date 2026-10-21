#!/usr/bin/env python3
"""Matched-bin nominal vs restricted DVCS extraction diagnostics; PNG and CSV only."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PERIODS=['10.6 GeV','Fa18 Inb','Fa18 Out','Sp18 Inb','Sp18 Out','Sp19 Inb']

def val(s,part=0):
    a=s.astype('string').str.strip().str.strip('"').str.strip('()').str.split(',',regex=False).str[part]
    return pd.to_numeric(a,errors='coerce')

def get(d,c,part=0):
    if c not in d: return pd.Series(np.nan,index=d.index,dtype=float)
    if pd.api.types.is_numeric_dtype(d[c]): return pd.to_numeric(d[c],errors='coerce') if part==0 else pd.Series(np.nan,index=d.index,dtype=float)
    return val(d[c],part)

def ratio(a,b):
    return (a/b).where(np.isfinite(a)&np.isfinite(b)&(b>0)&(a>0))

def draw(df,x,y,out,title,ylabel,ylim=None):
    good=np.isfinite(df[x])&np.isfinite(df[y])&(df[y]>0)
    a=df.loc[good,[x,y]].sort_values(x)
    if a.empty:return
    fig,ax=plt.subplots(figsize=(9,5.4))
    ax.scatter(a[x],a[y],s=9,alpha=.24,rasterized=True)
    edges=np.linspace(a[x].min(),a[x].max(),13)
    if edges[-1]>edges[0]:
        bins=pd.cut(a[x],edges,include_lowest=True)
        med=a.groupby(bins,observed=True)[y].median()
        counts=a.groupby(bins,observed=True)[y].size()
        centers=np.array([(i.left+i.right)/2 for i in med.index])
        ax.plot(centers,med.to_numpy(),'o-',lw=2,label='Median per interval')
        for xc,ym,n in zip(centers,med,counts):
            if np.isfinite(ym):ax.annotate(str(n),(xc,ym),xytext=(0,7),textcoords='offset points',ha='center',fontsize=7)
    ax.axhline(1,color='gray',ls='--',lw=1)
    ax.set(xlabel={'t':'|t| (GeV²)','xb':'x_B','q2':'Q² (GeV²)'}[x],ylabel=ylabel,title=title)
    if ylim:ax.set_ylim(*ylim)
    ax.grid(alpha=.2);ax.legend(loc='best');fig.tight_layout();fig.savefig(out,dpi=180);plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--nominal',default='output/csvs/dvcs_pass2_INCLUSIVE.csv')
    p.add_argument('--restricted',default='output/csvs/dvcs_pass2_proton_p0p5_theta58.csv')
    p.add_argument('--outdir',default='output/proton_cut_comparison/direct_matched')
    args=p.parse_args();out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    a=pd.read_csv(args.nominal,low_memory=False).set_index('bin index',verify_integrity=True)
    b=pd.read_csv(args.restricted,low_memory=False).set_index('bin index',verify_integrity=True)
    if not a.index.equals(b.index):
        a,b=a.align(b,join='inner',axis=0)
    rows=[];summary=[]
    for period in PERIODS:
        combined=period=='10.6 GeV'
        xc=f'cross sections, ep->epg, exp, {period}, unpol'
        if xc not in a or xc not in b:continue
        d=pd.DataFrame(index=a.index)
        d['bin index']=a.index
        d['period']=period
        for short,col in [('t',f't_abs_avg, {period}'),('xb',f'xBavg, {period}'),('q2',f'Q2avg, {period}')]:
            d[short]=get(a,col).combine_first(get(b,col))
        d['sigma_nominal']=get(a,xc);d['sigma_restricted']=get(b,xc)
        d['stat_nominal']=get(a,xc,1);d['stat_restricted']=get(b,xc,1)
        d['xs_ratio']=ratio(d['sigma_restricted'],d['sigma_nominal'])
        d['stat_fraction_nominal']=ratio(d['stat_nominal'],d['sigma_nominal'])
        d['stat_fraction_restricted']=ratio(d['stat_restricted'],d['sigma_restricted'])
        d['stat_fraction_ratio']=ratio(d['stat_fraction_restricted'],d['stat_fraction_nominal'])
        # The per-period acceptance field is stored as a tuple: (value,error,...).
        if not combined:
            d['acc_nominal']=get(a,'acceptance, '+period)
            d['acc_restricted']=get(b,'acceptance, '+period)
            d['acc_ratio']=ratio(d['acc_restricted'],d['acc_nominal'])
            for kind,typ in [('data','exp'),('mc','mc')]:
                for topology in ['FD, FD','CD, FD','CD, FT']:
                    if kind=='data':
                        c=f'raw yield, ep->epg, ({topology}), exp, {period}, unpol'
                    else:
                        c=f'reconstructed current corrected yield, ep->epg, ({topology}), mc, {period}'
                    key=f'{kind}_{topology.replace(", ","_")}'
                    d[key+'_nominal']=get(a,c)
                    d[key+'_restricted']=get(b,c)
                keys=[f'{kind}_{top.replace(", ","_")}' for top in ['FD, FD','CD, FD','CD, FT']]
                for suffix in ['nominal','restricted']:
                    v=d[[k+'_'+suffix for k in keys]]
                    d[kind+'_total_'+suffix]=v.sum(axis=1,min_count=1)
                    d[kind+'_cd_'+suffix]=v[[keys[1]+'_'+suffix,keys[2]+'_'+suffix]].sum(axis=1,min_count=1)
                d[kind+'_ratio']=ratio(d[kind+'_total_restricted'],d[kind+'_total_nominal'])
                d[kind+'_cd_ratio']=ratio(d[kind+'_cd_restricted'],d[kind+'_cd_nominal'])
        rows.append(d)
        sub=d[np.isfinite(d['xs_ratio'])]
        summary.append({'period':period,'matched_positive_xs':len(sub),'median_xs_ratio':sub['xs_ratio'].median(),
                        'p10_xs_ratio':sub['xs_ratio'].quantile(.1),'p90_xs_ratio':sub['xs_ratio'].quantile(.9),
                        'median_stat_fraction_nominal':sub['stat_fraction_nominal'].median(),
                        'median_stat_fraction_restricted':sub['stat_fraction_restricted'].median(),
                        'median_acceptance_ratio':sub['acc_ratio'].median() if 'acc_ratio' in sub else np.nan,
                        'median_data_yield_ratio':sub['data_ratio'].median() if 'data_ratio' in sub else np.nan,
                        'median_mc_yield_ratio':sub['mc_ratio'].median() if 'mc_ratio' in sub else np.nan})
        folder=out/period.replace(' ','_');folder.mkdir(exist_ok=True)
        for x in ['t','xb','q2']:
            draw(d,x,'xs_ratio',folder/f'cross_section_ratio_vs_{x}.png',period+' — matched-bin cross sections','Restricted / nominal cross section')
        for y,label in [('stat_fraction_ratio','Ratio of relative statistical errors'),('acc_ratio','Restricted / nominal acceptance'),('data_ratio','Restricted / nominal reconstructed data yield'),('mc_ratio','Restricted / nominal reconstructed MC yield'),('data_cd_ratio','Restricted / nominal CD data yield'),('mc_cd_ratio','Restricted / nominal CD MC yield')]:
            if y in d:draw(d,'t',y,folder/f'{y}_vs_t.png',period+' — '+label,label)
        if not combined:
            for key in ['data_total','mc_total','data_cd','mc_cd']:
                for suffix in ['nominal','restricted']:
                    y=key+'_'+suffix
                    draw(d,'t',y,folder/f'{y}_vs_t.png',period+' — '+y.replace('_',' '),y.replace('_',' '))
    pd.concat(rows,ignore_index=True).to_csv(out/'matched_bin_diagnostics.csv',index=False)
    pd.DataFrame(summary).to_csv(out/'summary_by_period.csv',index=False)
    print(pd.DataFrame(summary).to_string(index=False))
    print('Saved PNGs and tables to',out)

if __name__=='__main__':main()
