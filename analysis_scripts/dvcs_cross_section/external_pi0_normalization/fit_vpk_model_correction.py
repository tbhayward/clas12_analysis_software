#!/usr/bin/env python3
"""
Diagnose the external AAOgen/VPK pi0 model correction using the Dlamini 2021
Hall-A comparison produced by validate_aao_pi0_model.py.

This intentionally does NOT yet modify any RGA eppi0 normalization.  It asks:
  1. Is a Q2-only correction adequate at the nine Hall-A settings?
  2. Does adding xB materially improve the setting-level description?
  3. Is there reproducible t' dependence inside settings?

The fitted quantity is log(data/VPK), so the correction is positive and
multiplicative.  Fits are diagnostic; Dlamini systematic correlations are not
fully available in the published table, so do not interpret chi2 as a final
statistical test.

Run from external_pi0_normalization:
    python3 fit_vpk_model_correction.py
"""

from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
POINTS = HERE/"output"/"dlamini2021_vpk_point_comparison.csv"
OUT = HERE/"output"/"model_correction"
PNG = OUT/"png"
QREF = 4.5
XREF = 0.46


def wls(X, y, sy):
    w = 1.0/np.square(sy)
    xtw = X.T*w
    cov = np.linalg.inv(xtw@X)
    beta = cov@(xtw@y)
    resid = y-X@beta
    chi2 = float(np.sum(np.square(resid/sy)))
    return beta, cov, chi2, resid


def design(df, kind):
    lq = np.log(df.Q2_GeV2.to_numpy(float)/QREF)
    dx = df.xB_mean.to_numpy(float)-XREF
    one = np.ones(len(df))
    if kind == "constant": return np.c_[one], ["a"]
    if kind == "Q2": return np.c_[one,lq], ["a","b_lnQ2"]
    if kind == "xB": return np.c_[one,dx], ["a","c_xB"]
    if kind == "Q2_xB": return np.c_[one,lq,dx], ["a","b_lnQ2","c_xB"]
    raise ValueError(kind)


def fit_settings(points):
    # Reconstruct one ratio per setting from the point ratios. This matches the
    # first-stage summary, while keeping the implementation self-contained.
    rows=[]
    for key,g in points.groupby(["xB_mean","Q2_GeV2","Ebeam_GeV"],sort=True):
        r=g.data_over_vpk_U.to_numpy(float)
        e=g.data_over_vpk_U_err.to_numpy(float)
        w=1/e**2
        R=np.sum(w*r)/np.sum(w)
        s=math.sqrt(1/np.sum(w))
        rows.append(dict(xB_mean=key[0],Q2_GeV2=key[1],Ebeam_GeV=key[2],
                         ratio=R,ratio_err=s,n_points=len(g)))
    d=pd.DataFrame(rows)
    d["log_ratio"]=np.log(d.ratio)
    d["log_ratio_err"]=d.ratio_err/d.ratio

    fits=[]
    for kind in ["constant","Q2","xB","Q2_xB"]:
        X,names=design(d,kind)
        b,cov,chi2,res=wls(X,d.log_ratio.to_numpy(),d.log_ratio_err.to_numpy())
        k=len(b); n=len(d)
        # Small-sample corrected AIC, used only as a relative diagnostic.
        aic=chi2+2*k
        aicc=aic + (2*k*(k+1)/(n-k-1) if n>k+1 else np.inf)
        row={"model":kind,"n":n,"npar":k,"chi2":chi2,"dof":n-k,
             "chi2_per_dof":chi2/(n-k),"AICc":aicc}
        for i,name in enumerate(names):
            row[name]=b[i]; row[name+"_err"]=math.sqrt(cov[i,i])
        fits.append(row)

        if kind=="Q2_xB":
            d["pred_Q2_xB"]=np.exp(X@b)
            d["residual_Q2_xB"]=d.ratio-d.pred_Q2_xB
        if kind=="Q2":
            d["pred_Q2"]=np.exp(X@b)
            d["residual_Q2"]=d.ratio-d.pred_Q2
    return d,pd.DataFrame(fits)


def tprime_slopes(points):
    rows=[]
    for key,g in points.groupby(["xB_mean","Q2_GeV2","Ebeam_GeV"],sort=True):
        x=g.tprime_mean_GeV2.to_numpy(float)
        r=g.data_over_vpk_U.to_numpy(float)
        e=g.data_over_vpk_U_err.to_numpy(float)
        # Fit log ratio = intercept + slope*(t'-mean t').
        y=np.log(r); sy=e/r
        xc=x-np.mean(x)
        X=np.c_[np.ones(len(x)),xc]
        b,cov,chi2,res=wls(X,y,sy)
        rows.append(dict(xB_mean=key[0],Q2_GeV2=key[1],Ebeam_GeV=key[2],
                         mean_tprime_GeV2=np.mean(x),
                         log_slope_per_GeV2=b[1],
                         log_slope_err_per_GeV2=math.sqrt(cov[1,1]),
                         slope_significance=b[1]/math.sqrt(cov[1,1]),
                         chi2=chi2,dof=len(x)-2))
    return pd.DataFrame(rows)


def plots(settings, slopes):
    PNG.mkdir(parents=True,exist_ok=True)

    fig,ax=plt.subplots(figsize=(8,5.5))
    for xb,g in settings.groupby("xB_mean",sort=True):
        ax.errorbar(g.Q2_GeV2,g.ratio,yerr=g.ratio_err,fmt="o",capsize=3,
                    label=f"xB={xb:.2f}")
    q=np.linspace(settings.Q2_GeV2.min(),settings.Q2_GeV2.max(),250)
    # Show Q2-only curve from stored predictions by refitting.
    X,n=design(settings,"Q2")
    b,_,_,_=wls(X,settings.log_ratio.to_numpy(),settings.log_ratio_err.to_numpy())
    pred=np.exp(b[0]+b[1]*np.log(q/QREF))
    ax.plot(q,pred,"--",label="Q²-only diagnostic fit")
    ax.axhline(1,linewidth=1)
    ax.set(xlabel="Q² (GeV²)",ylabel="Dlamini / VPK",
           title="External VPK correction: setting-level behavior")
    ax.grid(alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(PNG/"setting_ratio_q2_diagnostic.png",dpi=220); plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,5.5))
    ax.errorbar(slopes.Q2_GeV2,slopes.log_slope_per_GeV2,
                yerr=slopes.log_slope_err_per_GeV2,fmt="o",capsize=3)
    ax.axhline(0,linewidth=1)
    for _,r in slopes.iterrows():
        ax.annotate(f"xB={r.xB_mean:.2f}",(r.Q2_GeV2,r.log_slope_per_GeV2),
                    xytext=(4,4),textcoords="offset points",fontsize=8)
    ax.set(xlabel="Q² (GeV²)",ylabel="slope of ln(data/VPK) vs t' (GeV⁻²)",
           title="Within-setting t' dependence of VPK residual")
    ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(PNG/"tprime_residual_slopes.png",dpi=220); plt.close(fig)


def main():
    if not POINTS.exists():
        raise SystemExit(f"Missing {POINTS}; run validate_aao_pi0_model.py first.")
    pts=pd.read_csv(POINTS)
    OUT.mkdir(parents=True,exist_ok=True)
    settings,fits=fit_settings(pts)
    slopes=tprime_slopes(pts)
    settings.to_csv(OUT/"setting_ratios_and_predictions.csv",index=False)
    fits.to_csv(OUT/"candidate_model_fits.csv",index=False)
    slopes.to_csv(OUT/"tprime_residual_slopes.csv",index=False)
    plots(settings,slopes)

    print("\n=== External VPK correction diagnostics ===")
    print("\nCandidate setting-level fits to ln(data/VPK):")
    print(fits[["model","npar","chi2","dof","chi2_per_dof","AICc"]].to_string(index=False))
    print("\nWithin-setting t' slopes:")
    print(slopes[["xB_mean","Q2_GeV2","log_slope_per_GeV2",
                  "log_slope_err_per_GeV2","slope_significance"]].to_string(index=False))
    print("\nInterpretation rule:")
    print("  Do not promote a correction to RGA yet. First inspect whether Q2-only")
    print("  residuals retain xB structure and whether t' slopes repeat across settings.")
    print(f"\nWrote diagnostics under {OUT}")


if __name__=="__main__":
    main()
