#!/usr/bin/env python3
"""
Stage 3: GK model bridge + blinded pi0 pseudo-data + first final-observable projections.

This stage is intentionally split into two modes:

  1) No --gk-results:
     validate Stage-2 inputs and write exact model-query templates/instructions.
     No model cross sections are invented.

  2) With --gk-results:
     require a CSV with GK partial cross sections
       point_id,sigma_T,sigma_L,sigma_TT,sigma_LT
     for the Stage-3 query points. Then construct blinded phi-dependent
     pseudo-data, covariance-aware harmonic fits, Rosenbluth L/T projections,
     and workshop-oriented diagnostic figures.

The preliminary CLAS12 cross-section central values are never used.
"""

from __future__ import annotations
import argparse, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REQ_GK = ["point_id","sigma_T","sigma_L","sigma_TT","sigma_LT"]

def args():
    here=Path(__file__).resolve().parent
    p=argparse.ArgumentParser()
    p.add_argument("--stage1",type=Path,default=here/"output"/"pi0_gk_stage1")
    p.add_argument("--stage2",type=Path,default=here/"output"/"pi0_gk_stage2")
    p.add_argument("--output",type=Path,default=here/"output"/"pi0_gk_stage3")
    p.add_argument("--gk-results",type=Path,default=None)
    p.add_argument("--rga-factor",type=float,default=1.5)
    p.add_argument("--rgk-factor",type=float,default=8.0)
    p.add_argument(
        "--exposure-scan-multipliers",
        type=str,
        default="0.25,0.5,1,2,4,8,16",
        help="Comma-separated multipliers applied jointly to the nominal "
             "RGA and RGK exposure factors. A multiplier of 1 corresponds "
             "to --rga-factor and --rgk-factor.",
    )
    p.add_argument("--seed",type=int,default=20260924)
    return p.parse_args()

def design(phi_deg, eps):
    ph=np.deg2rad(np.asarray(phi_deg,float))
    # coefficients multiply [T,L,TT,LT]
    return np.column_stack([
        np.ones(len(ph)),
        np.full(len(ph),eps),
        eps*np.cos(2*ph),
        np.sqrt(2*eps*(1+eps))*np.cos(ph)
    ])

def harmonic_design(phi_deg, eps):
    ph=np.deg2rad(np.asarray(phi_deg,float))
    # coefficients multiply [U,TT,LT], U=T+eps L
    return np.column_stack([
        np.ones(len(ph)),
        eps*np.cos(2*ph),
        np.sqrt(2*eps*(1+eps))*np.cos(ph)
    ])

def pinv_cov_fit(A,y,C):
    Ci=np.linalg.pinv(C,rcond=1e-12)
    N=A.T@Ci@A
    V=np.linalg.pinv(N,rcond=1e-12)
    b=V@(A.T@Ci@y)
    return b,V

def load_corr(stage2):
    return np.load(stage2/"tables"/"02_within_cell_phi_correlations.npz",allow_pickle=False)

def campaign_rows(stage1,campaign):
    f="01_blinded_rga_bins.csv" if campaign=="rga" else "02_blinded_rgk_bins.csv"
    return pd.read_csv(stage1/f)

def build_query_points(common):
    q=common.copy()
    q.insert(0,"point_id",[f"R{n:04d}" for n in range(len(q))])
    keep=["point_id","Q2_common_GeV2","xB_common","minus_t_common_GeV2","xi_common",
          "epsilon_rga","epsilon_rgk","delta_epsilon","nphi_rga","nphi_rgk",
          "iq2_rga","ixb_rga","it_rga","iq2_rgk","ixb_rgk","it_rgk"]
    return q[keep]

def write_model_template(q,path):
    t=q[["point_id","Q2_common_GeV2","xB_common","minus_t_common_GeV2"]].copy()
    for c in ["sigma_T","sigma_L","sigma_TT","sigma_LT"]:
        t[c]=np.nan
    t.to_csv(path,index=False)

def model_merge(q,gkfile):
    g=pd.read_csv(gkfile)
    miss=[c for c in REQ_GK if c not in g.columns]
    if miss: raise RuntimeError(f"GK results missing columns: {miss}")
    if g["point_id"].duplicated().any(): raise RuntimeError("Duplicate point_id in GK results")

    # The validated PARTONS grid is authoritative for Stage-3 usability.  A
    # common Stage-2 point can be absent because it fails the shared-y cut or
    # because the requested (Q2,xB,t) point is outside physical pi0 phase space.
    # Do not invent/interpolate a model value for such points: remove them from
    # the projection and propagate that same selection through all downstream
    # pseudo-data, covariance fits and Rosenbluth results.
    supplied=set(g["point_id"].astype(str))
    requested=set(q["point_id"].astype(str))
    extra=sorted(supplied-requested)
    if extra:
        raise RuntimeError(f"GK results contain {len(extra)} unknown point_id values, e.g. {extra[:5]}")

    excluded=q.loc[~q["point_id"].isin(supplied),"point_id"].tolist()
    q_use=q[q["point_id"].isin(supplied)].copy()
    m=q_use.merge(g[REQ_GK],on="point_id",how="left",validate="one_to_one")

    if m[REQ_GK[1:]].isna().any().any():
        bad=m.loc[m[REQ_GK[1:]].isna().any(axis=1),"point_id"].tolist()
        raise RuntimeError(f"Non-finite/missing GK structure functions for {len(bad)} retained points, e.g. {bad[:5]}")

    print("[Stage-3 validated-model selection]")
    print(f"  Stage-2 candidate points : {len(q)}")
    print(f"  validated GK points      : {len(m)}")
    print(f"  excluded from projection : {len(excluded)}")
    if excluded:
        print(f"  excluded point IDs       : {', '.join(excluded)}")
    return m

def covariance_for_cell(g, corr, factor, model_y):
    # Fractional uncertainties are transferred onto GK model values.
    frac=g["relative_uncertainty"].to_numpy(float)/math.sqrt(factor)
    sig=np.abs(model_y)*frac
    return np.outer(sig,sig)*corr

def make_pseudodata(model,stage1,stage2,rga_factor,rgk_factor):
    z=load_corr(stage2)
    rows=[]; fits=[]
    rga=campaign_rows(stage1,"rga"); rgk=campaign_rows(stage1,"rgk")
    for r in model.itertuples(index=False):
        all_y=[]; all_A=[]; blocks=[]
        for camp,df,factor,eps,iq,ix,it in [
            ("rga",rga,rga_factor,r.epsilon_rga,r.iq2_rga,r.ixb_rga,r.it_rga),
            ("rgk",rgk,rgk_factor,r.epsilon_rgk,r.iq2_rgk,r.ixb_rgk,r.it_rgk)]:
            g=df[(df.iq2==iq)&(df.ixb==ix)&(df.it==it)].sort_values("iphi")
            phi=g.phi_center_deg.to_numpy(float)
            A4=design(phi,eps)
            theta=np.array([r.sigma_T,r.sigma_L,r.sigma_TT,r.sigma_LT],float)
            y=A4@theta
            key=f"{camp}_{int(iq)}_{int(ix)}_{int(it)}"
            R=z[key]
            C=covariance_for_cell(g,R,factor,y)

            A3=harmonic_design(phi,eps)
            b3,V3=pinv_cov_fit(A3,y,C)
            fits.append(dict(point_id=r.point_id,campaign=camp,epsilon=eps,
                sigma_U_fit=b3[0],sigma_U_unc=np.sqrt(max(V3[0,0],0)),
                sigma_TT_fit=b3[1],sigma_TT_unc=np.sqrt(max(V3[1,1],0)),
                sigma_LT_fit=b3[2],sigma_LT_unc=np.sqrt(max(V3[2,2],0)),
                nphi=len(g)))
            for j,rr in enumerate(g.itertuples(index=False)):
                rows.append(dict(point_id=r.point_id,campaign=camp,
                    iq2=int(iq),ixb=int(ix),it=int(it),iphi=int(rr.iphi),
                    phi_deg=float(rr.phi_center_deg),epsilon=float(eps),
                    model_cross_section=float(y[j]),
                    projected_relative_uncertainty=float(rr.relative_uncertainty/math.sqrt(factor)),
                    projected_absolute_uncertainty=float(np.sqrt(max(C[j,j],0)))))
            all_y.append(y); all_A.append(A4); blocks.append(C)

        y=np.concatenate(all_y); A=np.vstack(all_A)
        C=np.zeros((len(y),len(y)))
        n0=0
        for B in blocks:
            n=len(B); C[n0:n0+n,n0:n0+n]=B; n0+=n
        b,V=pinv_cov_fit(A,y,C)
        # overwrite later via separate table
    return pd.DataFrame(rows),pd.DataFrame(fits)

def lt_from_u(model,hfits):
    out=[]
    for r in model.itertuples(index=False):
        a=hfits[(hfits.point_id==r.point_id)&(hfits.campaign=="rga")].iloc[0]
        b=hfits[(hfits.point_id==r.point_id)&(hfits.campaign=="rgk")].iloc[0]
        de=a.epsilon-b.epsilon
        L=(a.sigma_U_fit-b.sigma_U_fit)/de
        # independent campaigns; systematics not yet included
        varL=(a.sigma_U_unc**2+b.sigma_U_unc**2)/(de**2)
        T=a.sigma_U_fit-a.epsilon*L
        # derivatives for T wrt U_a,U_b
        da=1-a.epsilon/de
        db=a.epsilon/de
        varT=da*da*a.sigma_U_unc**2+db*db*b.sigma_U_unc**2
        covTL=(da*a.sigma_U_unc**2/de + db*(-b.sigma_U_unc**2/de))
        R=L/T if T!=0 else np.nan
        if T!=0:
            dRdL=1/T; dRdT=-L/T**2
            varR=dRdL*dRdL*varL+dRdT*dRdT*varT+2*dRdL*dRdT*covTL
        else: varR=np.nan
        out.append(dict(point_id=r.point_id,Q2_GeV2=r.Q2_common_GeV2,
            xB=r.xB_common,minus_t_GeV2=r.minus_t_common_GeV2,
            delta_epsilon=de,sigma_T_truth=r.sigma_T,sigma_L_truth=r.sigma_L,
            sigma_T_proj=T,sigma_T_unc=np.sqrt(max(varT,0)),
            sigma_L_proj=L,sigma_L_unc=np.sqrt(max(varL,0)),
            R_L_over_T_truth=r.sigma_L/r.sigma_T if r.sigma_T!=0 else np.nan,
            R_L_over_T_proj=R,
            R_L_over_T_unc=np.sqrt(max(varR,0)) if np.isfinite(varR) else np.nan,
            L_significance_abs=abs(L)/np.sqrt(max(varL,1e-300))))
    return pd.DataFrame(out)

def luminosity_scan(model,stage1,stage2,nominal_rga,nominal_rgk,multipliers):
    """
    Re-run the complete pseudo-data covariance propagation and Rosenbluth
    separation while scaling both beam-energy data sets.

    A multiplier m corresponds to:
        RGA exposure = m * nominal_rga
        RGK exposure = m * nominal_rgk

    This is the relevant scan when additional beam time improves both arms
    of the Rosenbluth separation, rather than imposing an artificial fixed-RGA
    statistical floor.
    """
    rows=[]
    per_point=[]

    for multiplier in multipliers:
        rga_factor=float(nominal_rga)*float(multiplier)
        rgk_factor=float(nominal_rgk)*float(multiplier)

        pseudo,hfits=make_pseudodata(
            model,stage1,stage2,rga_factor,rgk_factor
        )
        lt=lt_from_u(model,hfits)

        relL=lt.sigma_L_unc/np.maximum(np.abs(lt.sigma_L_truth),1e-300)
        sig=lt.L_significance_abs.to_numpy(float)

        finite_rel=relL[np.isfinite(relL)]
        finite_sig=sig[np.isfinite(sig)]

        rows.append(dict(
            exposure_multiplier=float(multiplier),
            rga_factor=rga_factor,
            rgk_factor=rgk_factor,
            n_points=len(lt),
            n_ge_1sigma=int(np.sum(finite_sig >= 1.0)),
            n_ge_2sigma=int(np.sum(finite_sig >= 2.0)),
            n_ge_3sigma=int(np.sum(finite_sig >= 3.0)),
            n_ge_5sigma=int(np.sum(finite_sig >= 5.0)),
            median_L_significance=float(np.nanmedian(finite_sig)),
            max_L_significance=float(np.nanmax(finite_sig)),
            median_rel_sigma_L=float(np.nanmedian(finite_rel)),
            p25_rel_sigma_L=float(np.nanpercentile(finite_rel,25)),
            p75_rel_sigma_L=float(np.nanpercentile(finite_rel,75)),
        ))

        t=lt[[
            "point_id","Q2_GeV2","xB","minus_t_GeV2","delta_epsilon",
            "sigma_T_truth","sigma_L_truth","sigma_L_unc",
            "L_significance_abs"
        ]].copy()
        t.insert(0,"rgk_factor",rgk_factor)
        t.insert(0,"rga_factor",rga_factor)
        t.insert(0,"exposure_multiplier",float(multiplier))
        t["relative_sigma_L_uncertainty"] = (
            t.sigma_L_unc/np.maximum(np.abs(t.sigma_L_truth),1e-300)
        )
        per_point.append(t)

    return pd.DataFrame(rows),pd.concat(per_point,ignore_index=True)


def luminosity_scan_figure(scan,outfile):
    fig,ax=plt.subplots(figsize=(7.2,5.4))
    ax.plot(scan.exposure_multiplier,scan.n_ge_1sigma,marker="o",label=r"$\geq1\sigma$")
    ax.plot(scan.exposure_multiplier,scan.n_ge_2sigma,marker="o",label=r"$\geq2\sigma$")
    ax.plot(scan.exposure_multiplier,scan.n_ge_3sigma,marker="o",label=r"$\geq3\sigma$")
    ax.axvline(1.0,ls="--",label="Nominal projected exposure")
    ax.set_xscale("log",base=2)
    ax.set_xlabel("Joint exposure multiplier")
    ax.set_ylabel(r"Number of bins with resolved $\sigma_L$")
    ax.set_title(r"$\pi^0$ Rosenbluth sensitivity versus joint exposure")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile,dpi=180)
    plt.close(fig)



def figures(model,pseudo,hfits,lt,out):
    out.mkdir(parents=True,exist_ok=True)

    fig,ax=plt.subplots(figsize=(7.2,5.4))
    sc=ax.scatter(model.Q2_common_GeV2,model.sigma_L/model.sigma_T,
                  c=model.minus_t_common_GeV2,s=28)
    fig.colorbar(sc,ax=ax,label=r"$-t$ (GeV$^2$)")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)"); ax.set_ylabel(r"GK $\sigma_L/\sigma_T$")
    ax.set_title(r"GK longitudinal-to-transverse ratio")
    fig.tight_layout(); fig.savefig(out/"01_gk_L_over_T_landscape.png",dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.2,5.4))
    ax.scatter(lt.Q2_GeV2,lt.L_significance_abs,c=lt.delta_epsilon,s=28)
    ax.axhline(2,ls="--"); ax.axhline(3,ls=":")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)"); ax.set_ylabel(r"Projected $|\sigma_L|/\delta\sigma_L$")
    ax.set_title(r"Projected Rosenbluth sensitivity to $\sigma_L$")
    fig.tight_layout(); fig.savefig(out/"02_sigmaL_significance_vs_q2.png",dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.2,5.4))
    rel=lt.sigma_L_unc/np.maximum(np.abs(lt.sigma_L_truth),1e-300)
    ax.scatter(lt.Q2_GeV2,rel,c=lt.minus_t_GeV2,s=28)
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)"); ax.set_ylabel(r"Projected $\delta\sigma_L/|\sigma_L|$")
    ax.set_title(r"Projected longitudinal cross-section precision")
    ax.set_ylim(0,min(3,max(1,np.nanpercentile(rel,95))))
    fig.tight_layout(); fig.savefig(out/"03_sigmaL_relative_precision_vs_q2.png",dpi=180); plt.close(fig)

    # Representative cell: maximize a simple combination of leverage and L significance,
    # restricted to cells with useful phi coverage through successful fits.
    cand=lt.replace([np.inf,-np.inf],np.nan).dropna(subset=["L_significance_abs"])
    rid=cand.sort_values(["L_significance_abs","delta_epsilon"],ascending=False).iloc[0].point_id
    rr=model[model.point_id==rid].iloc[0]
    pp=pseudo[pseudo.point_id==rid]
    fig,ax=plt.subplots(figsize=(7.2,5.4))
    for camp,mark in [("rga","o"),("rgk","s")]:
        g=pp[pp.campaign==camp].sort_values("phi_deg")
        ax.errorbar(g.phi_deg,g.model_cross_section,yerr=g.projected_absolute_uncertainty,
                    fmt=mark,ms=4,capsize=2,label=camp.upper())
    ax.set_xlabel(r"$\phi$ (degrees)")
    ax.set_ylabel(r"Reduced cross section (nb/GeV$^2$)")
    ax.set_title(fr"Representative blinded projection: $Q^2={rr.Q2_common_GeV2:.2f}$, "
                 fr"$x_B={rr.xB_common:.3f}$, $-t={rr.minus_t_common_GeV2:.2f}$")
    ax.legend(); fig.tight_layout()
    fig.savefig(out/"04_representative_rga_rgk_phi_projection.png",dpi=180); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.2,5.4))
    ratio_tt=np.abs(model.sigma_TT)/(model.sigma_T + 0.5*(model.epsilon_rga+model.epsilon_rgk)*model.sigma_L)
    ratio_lt=np.abs(model.sigma_LT)/(model.sigma_T + 0.5*(model.epsilon_rga+model.epsilon_rgk)*model.sigma_L)
    ax.scatter(model.Q2_common_GeV2,ratio_tt,s=26,label=r"$|\sigma_{TT}|/\sigma_U$")
    ax.scatter(model.Q2_common_GeV2,ratio_lt,s=26,label=r"$|\sigma_{LT}|/\sigma_U$")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)"); ax.set_ylabel("GK transverse/interference fraction")
    ax.set_title("Model expectation for approach toward longitudinal regime")
    ax.legend(); fig.tight_layout()
    fig.savefig(out/"05_gk_transverse_interference_fractions.png",dpi=180); plt.close(fig)

def main():
    a=args()
    s1=a.stage1.resolve(); s2=a.stage2.resolve(); out=a.output.resolve()
    tabs=out/"tables"; figs=out/"figures"; queries=out/"model_queries"
    for d in (out,tabs,figs,queries): d.mkdir(parents=True,exist_ok=True)

    common=pd.read_csv(s2/"tables"/"03_common_rosenbluth_model_points.csv")
    q=build_query_points(common)
    q.to_csv(queries/"01_common_gk_query_points.csv",index=False)
    write_model_template(q,queries/"02_gk_structure_function_results_template.csv")

    instructions = """Stage-3 GK bridge
=================
Fill 02_gk_structure_function_results_template.csv with GK/PARTONS partial
cross sections for neutral-pion production at each common point.

Required columns:
  point_id, sigma_T, sigma_L, sigma_TT, sigma_LT

The validated PARTONS/GK bridge supplies all four quantities in nb/GeV^2.
Stage-3 uses only point_ids present in that validated model product; candidate
points excluded by the shared-y or physical-phase-space selections are omitted.

PARTONS documentation confirms DVMPProcessGK06 contains partial-cross-section
methods CrossSectionT, CrossSectionL, CrossSectionTT, and CrossSectionLT.
This script does not guess a PARTONS CLI/database serialization format.
Once an actual PARTONS result is available, pass the converted CSV with:
  python prepare_pi0_gk_stage3_projection.py --gk-results <file.csv>
"""
    (queries/"README_GK_RESULTS.txt").write_text(instructions)

    if a.gk_results is None:
        summary=[
            "CLAS12 pi0 blinded projection -- stage 3 model bridge",
            "="*54,"",
            f"Prepared {len(q)} common GK query points.",
            "No --gk-results file supplied, so no model cross sections or pseudo-data",
            "were invented. Fill/convert the Stage-3 GK result template and rerun.",
            "",
            "The eventual workshop observables are already encoded in the projection:",
            "  * sigma_L/sigma_T and its projected uncertainty;",
            "  * sigma_L significance versus Q2;",
            "  * covariance-aware sigma_U, sigma_TT, sigma_LT harmonic fits;",
            "  * representative RGA/RGK phi projections;",
            "  * transverse/interference fractions versus Q2.",
        ]
        (out/"summary.txt").write_text("\n".join(summary)+"\n")
        print("\n".join(summary)); print(f"\nWrote Stage-3 model bridge to {out}")
        return

    model=model_merge(q,a.gk_results.resolve())
    model.to_csv(tabs/"01_gk_common_structure_functions.csv",index=False)
    pseudo,hfits=make_pseudodata(model,s1,s2,a.rga_factor,a.rgk_factor)
    pseudo.to_csv(tabs/"02_blinded_model_pseudodata.csv",index=False)
    hfits.to_csv(tabs/"03_covariance_aware_harmonic_fits.csv",index=False)
    lt=lt_from_u(model,hfits)
    lt.to_csv(tabs/"04_rosenbluth_LT_projection.csv",index=False)
    figures(model,pseudo,hfits,lt,figs)

    try:
        exposure_scan_multipliers=[
            float(x.strip()) for x in a.exposure_scan_multipliers.split(",")
            if x.strip()
        ]
    except ValueError as exc:
        raise RuntimeError(
            "Could not parse --exposure-scan-multipliers="
            f"{a.exposure_scan_multipliers!r}"
        ) from exc
    if (not exposure_scan_multipliers
            or any(x <= 0 for x in exposure_scan_multipliers)):
        raise RuntimeError(
            "--exposure-scan-multipliers must contain positive values"
        )

    scan,scan_points=luminosity_scan(
        model,s1,s2,a.rga_factor,a.rgk_factor,exposure_scan_multipliers
    )
    scan.to_csv(tabs/"05_joint_luminosity_scan_summary.csv",index=False)
    scan_points.to_csv(tabs/"06_joint_luminosity_scan_by_point.csv",index=False)
    luminosity_scan_figure(
        scan,figs/"06_joint_luminosity_scan_sigmaL_counts.png"
    )

    print("\n[Joint RGA + RGK luminosity scan]")
    print(f"  nominal RGA factor : {a.rga_factor:g}x")
    print(f"  nominal RGK factor : {a.rgk_factor:g}x")
    print(scan.to_string(
        index=False,
        columns=[
            "exposure_multiplier","rga_factor","rgk_factor",
            "n_ge_1sigma","n_ge_2sigma","n_ge_3sigma",
            "max_L_significance","median_rel_sigma_L"
        ],
        formatters={
            "exposure_multiplier":lambda x:f"{x:g}",
            "rga_factor":lambda x:f"{x:g}",
            "rgk_factor":lambda x:f"{x:g}",
            "max_L_significance":lambda x:f"{x:.3f}",
            "median_rel_sigma_L":lambda x:f"{x:.3f}",
        }
    ))

    good2=int((lt.L_significance_abs>=2).sum()); good3=int((lt.L_significance_abs>=3).sum())
    relL=lt.sigma_L_unc/np.maximum(np.abs(lt.sigma_L_truth),1e-300)
    summary=[
        "CLAS12 pi0 blinded projection -- stage 3",
        "="*43,"",
        f"GK common points: {len(model)}",
        f"Projected cells with |sigma_L| >= 2 sigma: {good2}",
        f"Projected cells with |sigma_L| >= 3 sigma: {good3}",
        f"Median projected delta sigma_L / |sigma_L|: {np.nanmedian(relL):.4f}",
        "",
        "Projection uses model central values only. Preliminary CLAS12 central",
        "cross sections are not used. Experimental information enters through",
        "relative uncertainties and within-cell phi correlations.",
        "",
        "",
        f"Joint luminosity scan: nominal RGA={a.rga_factor:g}x, "
        f"nominal RGK={a.rgk_factor:g}x; both scaled by multipliers = "
        + ", ".join(f"{x:g}x" for x in exposure_scan_multipliers),
        "Scan outputs:",
        "  tables/05_joint_luminosity_scan_summary.csv",
        "  tables/06_joint_luminosity_scan_by_point.csv",
        "  figures/06_joint_luminosity_scan_sigmaL_counts.png",
        "",
        "Current caveat: all supplied fractional uncertainty components are scaled",
        "as 1/sqrt(exposure); finite-MC and additional systematic floors are not yet",
        "separated. RGA and RGK are treated as independent between campaigns.",
    ]
    (out/"summary.txt").write_text("\n".join(summary)+"\n")
    print("\n".join(summary)); print(f"\nWrote Stage-3 projection to {out}")

if __name__=="__main__":
    main()
