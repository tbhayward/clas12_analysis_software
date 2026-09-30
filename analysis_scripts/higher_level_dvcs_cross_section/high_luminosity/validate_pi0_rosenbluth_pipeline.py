#!/usr/bin/env python3
"""Closure test for the pi0 native->shared bin-centering + Rosenbluth pipeline.

Builds deterministic native-coordinate pseudo-measurements from the validated GK
shared responses and the saved native->shared phi correction factors.  It then
applies those same corrections, fits U at each shared epsilon, performs the
Rosenbluth separation, and verifies recovery of the input GK L/T.
"""
from __future__ import annotations
import argparse, math
from pathlib import Path
import numpy as np
import pandas as pd

M=0.9382720813

def epsilon(Q2,xB,E):
    y=Q2/(2*M*E*xB); g2=4*M*M*xB*xB/Q2
    return (1-y-0.25*g2*y*y)/(1-y+0.5*y*y+0.25*g2*y*y)

def design(phi,eps):
    p=np.deg2rad(np.asarray(phi,float))
    return np.column_stack([np.ones(len(p)),eps*np.cos(2*p),np.sqrt(2*eps*(1+eps))*np.cos(p)])

def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser()
    ap.add_argument('--grid',type=Path,default=here/'output/pi0_gk_stage3/partons_gk/production_grid/gk_pi0_shared_physical_structure_functions.csv')
    ap.add_argument('--corrections',type=Path,default=here/'output/pi0_gk_stage3/partons_gk/06_gk_native_to_shared_phi_corrections.csv')
    ap.add_argument('--exclusions',type=Path,default=here/'output/pi0_gk_stage3/tables/00_gk_model_exclusions.csv')
    ap.add_argument('--tolerance',type=float,default=1e-8)
    a=ap.parse_args()
    g=pd.read_csv(a.grid); c=pd.read_csv(a.corrections)
    n_raw=len(g)
    excluded_ids=set()
    if a.exclusions.exists():
        ex=pd.read_csv(a.exclusions)
        if 'point_id' not in ex.columns:
            raise RuntimeError(f'Exclusion table has no point_id column: {a.exclusions}')
        excluded_ids=set(ex.point_id.dropna().astype(str))
    g=g[~g.point_id.astype(str).isin(excluded_ids)].copy()
    n_excluded=n_raw-len(g)
    cols=['dsigma_T_dt_nb_per_GeV2','dsigma_L_dt_nb_per_GeV2','dsigma_TT_dt_nb_per_GeV2','dsigma_LT_dt_nb_per_GeV2']
    if not np.isfinite(g[cols].to_numpy(float)).all():
        raise RuntimeError('Non-finite structure function in validated GK closure grid')
    if not (g.dsigma_T_dt_nb_per_GeV2.astype(float)>0).all():
        bad=g.loc[g.dsigma_T_dt_nb_per_GeV2.astype(float)<=0,'point_id'].astype(str).tolist()
        raise RuntimeError(f'Non-positive sigma_T remains in validated GK closure grid: {bad[:10]}')
    print('[Rosenbluth pipeline closure]')
    print(f'  raw GK points       : {n_raw}')
    print(f'  excluded here       : {n_excluded}')
    print(f'  validated GK points : {len(g)}')
    rows=[]
    for r in g.itertuples(index=False):
        U={}
        for camp,E,eps in [('rga',10.604,float(r.epsilon_rga)),('rgk',6.535,float(r.epsilon_rgk))]:
            cc=c[(c.point_id.astype(str)==str(r.point_id))&(c.campaign.astype(str).str.lower()==camp)].sort_values('phi_deg')
            if len(cc)<4: raise RuntimeError(f'{r.point_id}/{camp}: insufficient correction rows')
            ph=cc.phi_deg.to_numpy(float)
            shared=(float(r.dsigma_T_dt_nb_per_GeV2)+eps*float(r.dsigma_L_dt_nb_per_GeV2)
                    +eps*np.cos(2*np.deg2rad(ph))*float(r.dsigma_TT_dt_nb_per_GeV2)
                    +np.sqrt(2*eps*(1+eps))*np.cos(np.deg2rad(ph))*float(r.dsigma_LT_dt_nb_per_GeV2))
            native=shared/cc.gk_shared_over_native.to_numpy(float)
            corrected=native*cc.gk_shared_over_native.to_numpy(float)
            b=np.linalg.lstsq(design(ph,eps),corrected,rcond=None)[0]
            U[camp]=float(b[0])
        de=float(r.epsilon_rga-r.epsilon_rgk)
        L=(U['rga']-U['rgk'])/de; T=U['rga']-float(r.epsilon_rga)*L
        truth=float(r.dsigma_L_dt_nb_per_GeV2/r.dsigma_T_dt_nb_per_GeV2)
        rec=L/T
        rows.append(dict(point_id=r.point_id,input_L_over_T=truth,recovered_L_over_T=rec,difference=rec-truth))
    out=pd.DataFrame(rows)
    mx=float(np.max(np.abs(out.difference)))
    print(out.to_string(index=False))
    print(f'\nmax |recovered-input| = {mx:.3e}')
    if mx>a.tolerance: raise RuntimeError(f'closure failed: {mx:.3e} > {a.tolerance:.3e}')
    print('CLOSURE PASSED')

if __name__=='__main__': main()
