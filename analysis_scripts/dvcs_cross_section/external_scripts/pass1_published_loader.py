#!/usr/bin/env python3
from pathlib import Path
import math
import warnings
import numpy as np
import pandas as pd
PASS1_GLOBAL_NORM_FRAC=0.31

def _parse_e214m1(path):
    rows=[]
    for line in Path(path).open(encoding='utf-8'):
        p=line.rstrip('\n').split('\t')
        if len(p)<8: continue
        try: rows.append([int(p[0]),*map(float,p[1:8])])
        except (TypeError,ValueError): continue
    if not rows: raise ValueError(f'No E214M1 data rows parsed from {path}')
    return pd.DataFrame(rows,columns=['published_bin','published_xB','published_Q2','published_t_abs','published_phi','published_xs','published_stat','published_syst_total'])

def load_published_pass1_dataframe(published_path,legacy_boundaries_path=None,norm_frac=PASS1_GLOBAL_NORM_FRAC):
    published_path=Path(published_path)
    if legacy_boundaries_path is None: legacy_boundaries_path=published_path.with_name('all_bin_v3.csv')
    legacy=pd.read_csv(legacy_boundaries_path,low_memory=False)
    if 'valid bin' in legacy.columns: legacy=legacy.loc[pd.to_numeric(legacy['valid bin'],errors='coerce')==1].copy()
    legacy['_bin_int']=pd.to_numeric(legacy['Bin Name'],errors='coerce'); legacy['_phi_float']=pd.to_numeric(legacy['phiavg'],errors='coerce')
    pub=_parse_e214m1(published_path); merged=[]; unmatched=[]

    # The final publication re-numbered some kinematic bins relative to the
    # preliminary all_bin_v3.csv.  Therefore do NOT identify bins by integer
    # label.  Match each published (xB,Q2,|t|) group to the nearest legacy
    # kinematic group, then use the nearest legacy phi row only to recover the
    # original bin boundaries/bookkeeping.  The central kinematics and all
    # cross-section quantities below always come from E214M1.
    legacy_groups=(legacy.groupby("Bin Name",sort=False)
        .agg(_gx=("xBavg","mean"),_gq=("Q2avg","mean"),_gt=("t_abs_avg","mean"))
        .reset_index())
    pub_groups=(pub.groupby("published_bin",sort=False)
        .agg(_px=("published_xB","mean"),_pq=("published_Q2","mean"),_pt=("published_t_abs","mean"))
        .reset_index())
    bin_map={}
    for _, pg in pub_groups.iterrows():
        # Scale only makes the three coordinates commensurate; matches are
        # extremely close (publication rounding) and one-to-one.
        dist=np.sqrt(((legacy_groups["_gx"]-pg["_px"])/0.03)**2 +
                     ((legacy_groups["_gq"]-pg["_pq"])/0.30)**2 +
                     ((legacy_groups["_gt"]-pg["_pt"])/0.10)**2)
        ii=int(np.nanargmin(dist.to_numpy(float)))
        bin_map[int(pg["published_bin"])]=legacy_groups.iloc[ii]["Bin Name"]

    for p in pub.itertuples(index=False):
        legacy_bin=bin_map.get(int(p.published_bin))
        cand=legacy.loc[legacy["Bin Name"]==legacy_bin]
        if cand.empty: unmatched.append((p.published_bin,p.published_phi)); continue
        d=np.abs(cand['_phi_float'].to_numpy(float)-p.published_phi)
        if not np.any(np.isfinite(d)): unmatched.append((p.published_bin,p.published_phi)); continue
        row=legacy.loc[cand.index[int(np.nanargmin(d))]].copy()
        # Keep the legacy Bin Name because downstream pass-2 matching uses the
        # original analysis binning; retain the published bin id separately.
        row['published pass1 bin']=int(p.published_bin)
        row['xBavg']=p.published_xB; row['Q2avg']=p.published_Q2; row['t_abs_avg']=p.published_t_abs; row['phiavg']=p.published_phi
        row['cross sections, ep->epg, exp']=p.published_xs; row['cross sections, ep->epg, exp, stat. unc.']=p.published_stat
        norm=norm_frac*abs(p.published_xs); ptp=math.sqrt(max(p.published_syst_total**2-norm**2,0.0))
        row['cross sections, ep->epg, exp, syst. unc. (up)']=ptp; row['cross sections, ep->epg, exp, syst. unc. (down)']=ptp
        row['published pass1 total syst. unc.']=p.published_syst_total; row['published pass1 inferred point-to-point syst. unc.']=ptp
        row['published pass1 normalization syst. unc.']=norm; row['published pass1 source']='CLAS Physics Database E214M1'; row['valid bin']=1
        merged.append(row)
    if unmatched: warnings.warn(f'{len(unmatched)} E214M1 points could not be matched to legacy boundaries')
    out=pd.DataFrame(merged).drop(columns=['_bin_int','_phi_float'],errors='ignore').reset_index(drop=True)
    nbelow=int(np.sum(pub['published_syst_total'].to_numpy(float)<norm_frac*np.abs(pub['published_xs'].to_numpy(float))))
    print(f'[pass1-published] Loaded {len(out)}/{len(pub)} E214M1 points; legacy boundaries={legacy_boundaries_path}; normalization={100*norm_frac:.1f}%.')
    if nbelow: print(f'[pass1-published] WARNING: {nbelow} total systematic values are below {100*norm_frac:.1f}% of central value; inferred point-to-point variance clipped to zero.')
    return out
