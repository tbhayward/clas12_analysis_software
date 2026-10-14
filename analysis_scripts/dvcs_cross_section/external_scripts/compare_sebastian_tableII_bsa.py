#!/usr/bin/env python3
"""Direct Fa18 Pass-2 BSA comparison with Alvarado et al. PRL draft Table II.

No KM15/CFF model, no rebinning of already-extracted points. The C++ driver
produces event-level BSA in the exact published (Q2,xB,|t|,phi) binning.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KEY=['Q2min','Q2max','xBmin','xBmax','t_abs_min','t_abs_max','phi_slot']

def model(phi,a1,a2):
    return a1*np.sin(np.deg2rad(phi))/(1+a2*np.cos(np.deg2rad(phi)))

def fit(d):
    d=d[np.isfinite(d['BSA']) & np.isfinite(d['stat']) & (d['stat']>0)]
    if len(d)<5:return [np.nan]*5
    try:
        p,cov=curve_fit(model,d['phi'].to_numpy(),d['BSA'].to_numpy(),
                        sigma=d['stat'].to_numpy(),absolute_sigma=True,p0=[0.2,0],
                        bounds=([-2,-0.95],[2,0.95]),maxfev=20000)
        e=np.sqrt(np.diag(cov))
        chi2=np.sum(((d['BSA']-model(d['phi'],*p))/d['stat'])**2)
        return [*p,*e,chi2/(len(d)-2)]
    except (ValueError,RuntimeError):return [np.nan]*5

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--our',default='output/csvs/dvcs_pass2_sebastian_binning.csv')
    ap.add_argument('--reference',default='imports/sebastian_tableII_proton_detected.csv')
    ap.add_argument('--output',default='output/sebastian_bsa_cross_check/direct_comparison')
    args=ap.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    raw=pd.read_csv(args.our,low_memory=False)
    ref=pd.read_csv(args.reference)
    for k in KEY[:-1]:ref[k]=pd.to_numeric(ref[k])
    # Use the measured Fa18 event-weighted azimuth, never the geometric midpoint.
    phi_col = 'phiavg, Fa18'
    if phi_col not in raw.columns:
        raise KeyError(f'Missing {phi_col!r}; rerun the updated C++ analysis '
                       'to populate Sebastian-binning kinematic means.')
    raw['phi'] = pd.to_numeric(raw[phi_col], errors='coerce')
    if not raw['phi'].notna().any():
        raise ValueError(f'{phi_col!r} is entirely empty. The Sebastian-binning '
                         'CSV was not passed through update_bin_means_csv; '
                         'rerun the updated C++ analysis.')
    ref['phi']=pd.to_numeric(ref['phi'])
    def parse_tuple(x):
        if pd.isna(x):return np.nan,np.nan
        x=str(x).strip().strip('()')
        try:
            v=x.split(',');return float(v[0]),float(v[1])
        except (ValueError,IndexError):return np.nan,np.nan
    for group,label in [('Fa18','our_Fa18'),('Fa18 Inb','our_Fa18_Inb'),('Fa18 Out','our_Fa18_Out')]:
        c='BSA, counts, '+group
        if c not in raw:raise KeyError('Missing '+c)
        a=np.array([parse_tuple(x) for x in raw[c]],dtype=float)
        raw[label]=a[:,0];raw[label+'_stat']=a[:,1]
    if raw.duplicated(KEY).any():raise ValueError('Duplicate analysis bins')
    if ref.duplicated(KEY).any():raise ValueError('Duplicate published bins')
    merged=ref.merge(raw[KEY+['phi','our_Fa18','our_Fa18_stat','our_Fa18_Inb',
                              'our_Fa18_Inb_stat','our_Fa18_Out','our_Fa18_Out_stat']],
                     on=KEY,how='left',suffixes=('_published','_our'),validate='one_to_one',indicator=True)
    merged['delta']=merged['our_Fa18']-merged['BSA']
    merged['naive_pull']=merged['delta']/np.hypot(merged['our_Fa18_stat'],merged['stat'])
    merged.to_csv(out/'point_by_point.csv',index=False)
    summary=[]
    for cell,d in merged.groupby(KEY[:-1],sort=True):
        fpub=fit(d.rename(columns={'phi_published':'phi'}).assign(phi=d['phi_published']))
        four=fit(d.assign(BSA=d['our_Fa18'],stat=d['our_Fa18_stat'],phi=d['phi_our']))
        summary.append(dict(zip(KEY[:-1],cell),n_published=len(d),n_our_valid=int(d['our_Fa18'].notna().sum()),
                            published_a1=fpub[0],published_a2=fpub[1],published_a1_err=fpub[2],published_a2_err=fpub[3],published_chi2_ndf=fpub[4],
                            our_a1=four[0],our_a2=four[1],our_a1_err=four[2],our_a2_err=four[3],our_chi2_ndf=four[4],
                            delta_a1=four[0]-fpub[0],delta_a2=four[1]-fpub[1]))
        if d['our_Fa18'].notna().sum()==0:continue
        fig,ax=plt.subplots(figsize=(7.2,4.5))
        good_pub=d[np.isfinite(d['stat']) & (d['stat']>0)]
        ax.errorbar(good_pub['phi_published'],good_pub['BSA'],yerr=good_pub['stat'],fmt='o',ms=3.5,label='Alvarado et al. Table II (Fa18 Pass-2)')
        ax.errorbar(d['phi_our'],d['our_Fa18'],yerr=d['our_Fa18_stat'],fmt='s',ms=3.5,label='Our Fa18 Pass-2 combined')
        pp=np.linspace(0,360,360)
        if np.isfinite(fpub[0]):ax.plot(pp,model(pp,*fpub[:2]),alpha=.7)
        if np.isfinite(four[0]):ax.plot(pp,model(pp,*four[:2]),alpha=.7)
        ax.axhline(0,color='gray',lw=.7)
        ax.set(xlim=(0,360),xlabel=r'$\phi$ (degrees)',ylabel=r'$A_{LU}$',
               title=f'Q² [{cell[0]:g},{cell[1]:g}], xB [{cell[2]:g},{cell[3]:g}], |t| [{cell[4]:g},{cell[5]:g}]')
        ax.legend(fontsize=8);fig.tight_layout()
        fig.savefig(out/f'cell_{len(summary):03d}.png',dpi=130);plt.close(fig)
    summ=pd.DataFrame(summary);summ.to_csv(out/'amplitude_comparison.csv',index=False)
    print(f'Published points: {len(ref)}, cells: {len(summary)}')
    print(f'Matched rows: {(merged._merge=="both").sum()}; finite our BSA: {merged.our_Fa18.notna().sum()}')
    print('Table II includes anomalous entries (negative statistical error and means outside their nominal bins); see reference CSV.')
    print('Differences share experimental events: independent-error pulls are diagnostic, NOT calibrated significance.')
    print('Output:',out)

if __name__=='__main__':main()
