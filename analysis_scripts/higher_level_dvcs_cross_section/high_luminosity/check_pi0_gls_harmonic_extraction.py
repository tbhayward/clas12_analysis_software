#!/usr/bin/env python3
"""Audit the single-energy covariance-aware phi-harmonic extraction used by the pi0 L/T study."""
from __future__ import annotations
import importlib.util
import math
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
PROD=HERE/'pi0_LT_measurement_reach_study.py'
if not PROD.exists():
    PROD=Path('pi0_LT_measurement_reach_study.py').resolve()

def load_prod():
    spec=importlib.util.spec_from_file_location('reach',PROD)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def independent_gls(A,y,C):
    Ci=np.linalg.pinv(0.5*(C+C.T),rcond=1e-12)
    N=A.T@Ci@A
    cov=np.linalg.inv(N)
    th=cov@(A.T@Ci@y)
    return th,cov,float(np.linalg.cond(N)),np.linalg.matrix_rank(N)

def corr_from_cov(C):
    d=np.sqrt(np.clip(np.diag(C),0,np.inf)); den=np.outer(d,d)
    return np.divide(C,den,out=np.zeros_like(C),where=den>0)

def main():
    m=load_prod(); base=PROD.parent
    stage2=base/'output/pi0_gk_stage2'
    rga_file=base/'import/fa18_rosenbluth_inputs_20260924T165948Z/rga_10604/combined_reduced_cross_sections.csv'
    rgk_file=base/'import/fa18_rosenbluth_inputs_20260924T165948Z/rgk_6535/rgk6535_reduced_cross_sections.csv'
    corr_file=base/'output/pi0_gk_stage3/partons_gk/06_gk_native_to_shared_phi_corrections.csv'
    cov_file=stage2/'tables/02_within_cell_phi_correlations.npz'
    excl_file=base/'output/pi0_gk_stage3/tables/00_gk_model_exclusions.csv'
    common=pd.read_csv(stage2/'tables/03_common_rosenbluth_model_points.csv')
    rga=m._standardize_internal_cross_sections(rga_file,'RGA'); rgk=m._standardize_internal_cross_sections(rgk_file,'RGK')
    corr=pd.read_csv(corr_file); covz=np.load(cov_file,allow_pickle=False)
    excluded=set(pd.read_csv(excl_file).point_id.astype(str)) if excl_file.exists() else set()
    rows=[]; max_corr_change=0.0; max_helper_theta=0.0; max_helper_cov=0.0
    for n,r in enumerate(common.itertuples(index=False)):
        pid=f'R{n:04d}'
        if pid in excluded: continue
        q2=float(r.Q2_common_GeV2); xb=float(r.xB_common)
        for camp,df,ids,E in [('rga',rga,(int(r.iq2_rga),int(r.ixb_rga),int(r.it_rga)),10.604),('rgk',rgk,(int(r.iq2_rgk),int(r.ixb_rgk),int(r.it_rgk)),6.535)]:
            g=df[(df.iq2==ids[0])&(df.ixb==ids[1])&(df.it==ids[2])].copy()
            if len(g)<4: continue
            cc=corr[(corr.point_id.astype(str)==pid)&(corr.campaign.astype(str).str.lower()==camp)].sort_values('phi_deg')
            if len(cc)!=3: continue
            p=np.deg2rad(cc.phi_deg.to_numpy(float)); H=np.column_stack([np.ones(3),np.cos(p),np.cos(2*p)])
            nc=np.linalg.solve(H,cc.reduced_value_native.to_numpy(float)); sc=np.linalg.solve(H,cc.reduced_value_shared.to_numpy(float))
            pm=np.deg2rad(g.phi_deg.to_numpy(float)); Hm=np.column_stack([np.ones(len(g)),np.cos(pm),np.cos(2*pm)])
            ne=Hm@nc; se=Hm@sc
            scale=np.max(np.abs(ne))
            if not np.isfinite(scale) or scale==0 or np.any(np.abs(ne)<=100*np.finfo(float).eps*scale): continue
            fac=se/ne
            # Native covariance before deterministic bin-centering rescaling.
            gn,Cn=m._cell_covariance(g,covz,camp,*ids,rel_norm=0.0)
            # Apply correction in exactly the production ordering.
            fmap=dict(zip(g.iphi.to_numpy(int),fac))
            gs=gn.copy(); f=np.array([fmap[int(i)] for i in gs.iphi])
            gs['sigma']*=f; gs['delta_sigma']*=np.abs(f)
            gs['epsilon']=float(m._epsilon_from_q2_xb_E(q2,xb,E))
            _,Cs=m._cell_covariance(gs,covz,camp,*ids,rel_norm=0.0)
            D=np.diag(np.abs(f)); expected=D@Cn@D
            cov_scale_res=np.max(np.abs(Cs-expected))/max(np.max(np.abs(expected)),1e-300)
            corr_change=np.max(np.abs(corr_from_cov(Cs)-corr_from_cov(Cn))); max_corr_change=max(max_corr_change,corr_change)
            phi=np.deg2rad(gs.phi_deg.to_numpy(float)); eps=gs.epsilon.to_numpy(float)
            A=np.column_stack([np.ones(len(gs)),np.sqrt(2*eps*(1+eps))*np.cos(phi),eps*np.cos(2*phi)])
            y=gs.sigma.to_numpy(float)
            th,cov,cond,rank=independent_gls(A,y,Cs)
            thp,covp,_,_,condp=m._fit_single_energy_harmonics(gs,Cs)
            max_helper_theta=max(max_helper_theta,float(np.max(np.abs(th-thp))/max(np.max(np.abs(th)),1.0)))
            max_helper_cov=max(max_helper_cov,float(np.max(np.abs(cov-covp))/max(np.max(np.abs(cov)),1.0)))
            Cd=np.diag(np.diag(Cs)); thd,covd,_,_=independent_gls(A,y,Cd)
            sig=np.sqrt(np.clip(np.diag(cov),0,np.inf)); sigd=np.sqrt(np.clip(np.diag(covd),0,np.inf))
            rows.append(dict(point_id=pid,campaign=camp,nphi=len(gs),condition=cond,rank=rank,
                min_cov_eigenvalue=float(np.min(np.linalg.eigvalsh(0.5*(Cs+Cs.T)))),cov_scale_residual=cov_scale_res,
                corr_change=corr_change,dU=float(sig[0]),dLT=float(sig[1]),dTT=float(sig[2]),
                diag_over_full_U=float(sigd[0]/sig[0]),diag_over_full_LT=float(sigd[1]/sig[1]),diag_over_full_TT=float(sigd[2]/sig[2])))
    out=pd.DataFrame(rows)
    if out.empty: raise RuntimeError('No cells audited')
    print('\n[Step-5 GLS harmonic audit]')
    print(f'  campaign-cell fits audited             : {len(out)}')
    print(f'  max covariance DVD relative residual   : {out.cov_scale_residual.max():.3e}')
    print(f'  max correlation-matrix change after D  : {max_corr_change:.3e}')
    print(f'  max independent/production theta resid : {max_helper_theta:.3e}')
    print(f'  max independent/production cov resid   : {max_helper_cov:.3e}')
    print(f'  minimum normal-matrix rank              : {int(out["rank"].min())} / 3')
    print(f'  condition median / 95th / max           : {out.condition.median():.3g} / {out.condition.quantile(.95):.3g} / {out.condition.max():.3g}')
    print(f'  minimum covariance eigenvalue           : {out.min_cov_eigenvalue.min():.3e}')
    print('\nDiagonal-only uncertainty / full-GLS uncertainty:')
    for x in ['U','LT','TT']:
        v=out[f'diag_over_full_{x}']
        print(f'  {x:>2s}: median={v.median():.3f}, central 90%=[{v.quantile(.05):.3f},{v.quantile(.95):.3f}], min/max=[{v.min():.3f},{v.max():.3f}]')
    print('\nWorst-conditioned fits:')
    print(out.nlargest(10,'condition')[['point_id','campaign','nphi','condition','min_cov_eigenvalue']].to_string(index=False))
    effect=np.maximum.reduce([np.abs(out.diag_over_full_U-1),np.abs(out.diag_over_full_LT-1),np.abs(out.diag_over_full_TT-1)])
    tmp=out.assign(max_diag_fractional_effect=effect)
    print('\nLargest effect of dropping phi correlations:')
    print(tmp.nlargest(10,'max_diag_fractional_effect')[['point_id','campaign','nphi','diag_over_full_U','diag_over_full_LT','diag_over_full_TT','max_diag_fractional_effect']].to_string(index=False))
    dest=base/'output/pi0_gk_stage3/diagnostics/gls_harmonic_audit.csv'; dest.parent.mkdir(parents=True,exist_ok=True); out.to_csv(dest,index=False)
    print(f'\nWrote {dest}')
    return 0

if __name__=='__main__': raise SystemExit(main())
