#!/usr/bin/env python3
"""
Workshop-facing pi0 L/T reach plots.

Uses only workshop-safe projection tables. No measured CLAS12 central values are plotted.

Outputs:
  01_phase_space_fraction_vs_Q2.png
  02_reachable_Q2_vs_precision.png
  03_assumed_L_over_T_sensitivity.png
  04_GK_assumed_sigmaL_and_sensitivity.png
  05_synthetic_rosenbluth_examples.png
  phase_space_weighted_reach.csv

"future multiplier f" means the remaining running is acquired at f times the
present effective statistical yield. With the current recorded=remaining
bookkeeping, final effective statistics are (1+f) times today's recorded data.
"""

import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def pick(df, *names):
    for n in names:
        if n in df.columns:
            return n
    return None


def edges_width(df, base):
    """Resolve a nominal-bin width from several naming conventions."""
    lo = pick(df, f"{base}_low_rga", f"{base}_min_rga", f"{base}_lo_rga",
              f"{base}_low_GeV2", f"{base}_low", f"{base}_min", f"{base}_lo")
    hi = pick(df, f"{base}_high_rga", f"{base}_max_rga", f"{base}_hi_rga",
              f"{base}_high_GeV2", f"{base}_high", f"{base}_max", f"{base}_hi")
    if lo and hi:
        return df[hi].to_numpy(float)-df[lo].to_numpy(float), lo, hi
    return None, lo, hi


def attach_phase_space(points, common):
    common = common.copy().reset_index(drop=True)
    common.insert(0, "point_id", [f"R{i:04d}" for i in range(len(common))])

    # Keep only useful metadata not already present.
    keep = ["point_id"]
    for c in common.columns:
        if c.startswith(("iq2", "ixb", "it", "Q2_", "xB_", "minus_t_")):
            if c not in points.columns:
                keep.append(c)
    d = points.merge(common[keep], on="point_id", how="left", validate="many_to_one")

    dq, qlo, qhi = edges_width(d, "Q2")
    dx, xlo, xhi = edges_width(d, "xB")
    dt, tlo, thi = edges_width(d, "minus_t")

    # Try common alternative spelling for -t.
    if dt is None:
        dt, tlo, thi = edges_width(d, "mt")

    # For Q2-binned fractions the Q2 width cancels; the important transverse
    # phase-space area is Delta xB * Delta(-t). Require real bin edges rather
    # than silently reverting to equal-bin counting.
    if dx is None or dt is None:
        print("\nAvailable common-table columns:")
        print("  " + "\n  ".join(common.columns))
        raise RuntimeError(
            "Could not identify nominal xB and -t bin boundaries. "
            "Add the actual boundary column names to edges_width(); "
            "the script intentionally refuses equal-bin weighting."
        )

    d["phase_area_xB_t"] = dx*dt
    d["phase_volume_Q2_xB_t"] = (dq*dx*dt) if dq is not None else np.nan

    iq = pick(d, "iq2_rga", "iq2")
    if iq is None:
        raise RuntimeError("Could not identify nominal Q2-bin index.")
    d["Q2_bin"] = d[iq].astype(int)

    qcol = pick(d, "Q2_GeV2", "Q2_common_GeV2", "Q2_shared_GeV2")
    xbcol = pick(d, "xB", "xB_common", "xB_shared")
    if qcol is None or xbcol is None:
        raise RuntimeError("Could not identify shared Q2/xB coordinates.")
    d["Q2_GeV2"] = d[qcol].to_numpy(float)
    d["xB"] = d[xbcol].to_numpy(float)

    # Area-weighted mean Q2 used only as the display coordinate.
    return d


def weighted_fraction(g, mask):
    w = g.phase_area_xB_t.to_numpy(float)
    m = np.asarray(mask, bool)
    good = np.isfinite(w) & (w > 0)
    return float(np.sum(w[good & m])/np.sum(w[good])) if np.any(good) else np.nan


def make_weighted_table(d, thresholds):
    rows=[]
    for (f, iq), g in d.groupby(["future_luminosity_multiplier","Q2_bin"], sort=True):
        w=g.phase_area_xB_t.to_numpy(float)
        q=g.Q2_GeV2.to_numpy(float)
        qmean=float(np.average(q,weights=w))
        row=dict(future_multiplier=float(f), Q2_bin=int(iq),
                 Q2_weighted_mean_GeV2=qmean,
                 phase_area=float(np.sum(w)))
        for th in thresholds:
            row[f"fraction_deltaR_lt_{th:g}"]=weighted_fraction(
                g, g.delta_R_L_over_T.to_numpy(float)<th)
        rows.append(row)
    return pd.DataFrame(rows)


def plot_fraction_q2(tab, thresholds, outfile):
    fs=sorted(tab.future_multiplier.unique())
    fig, axes=plt.subplots(len(thresholds),1,figsize=(7.6,2.25*len(thresholds)),
                           sharex=True, sharey=True)
    axes=np.atleast_1d(axes)
    for ax,th in zip(axes,thresholds):
        for f in fs:
            g=tab[tab.future_multiplier==f].sort_values("Q2_weighted_mean_GeV2")
            ax.plot(g.Q2_weighted_mean_GeV2,
                    100*g[f"fraction_deltaR_lt_{th:g}"],
                    marker="o", label=f"remaining ×{f:g}")
        ax.set_ylabel("Coverage (%)")
        ax.set_title(rf"$\delta(\sigma_L/\sigma_T)<{th:g}$")
        ax.set_ylim(-3,103)
        ax.grid(alpha=.25)
    axes[-1].set_xlabel(r"$Q^2$ (GeV$^2$)")
    axes[0].legend(ncol=3,fontsize=8)
    fig.suptitle(r"Phase-space-weighted $L/T$ precision coverage",y=.995)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220)
    plt.close(fig)


def plot_reachable_q2(tab, outfile):
    # Instead of "at least one bin", require an explicit fraction of the xB,t
    # phase-space area at that Q2. Show several coverage definitions.
    precisions=np.linspace(.05,.50,19)
    coverages=[.25,.50,.75]
    fs=sorted(tab.future_multiplier.unique())
    fig, axes=plt.subplots(1,len(coverages),figsize=(13.2,4.2),sharey=True)
    for ax,cov in zip(axes,coverages):
        for f in fs:
            g=tab[tab.future_multiplier==f]
            qs=[]
            for p in precisions:
                # interpolate each Q2-bin's empirical CDF from the stored point
                # table is done elsewhere; here use nearest threshold columns
                # only when available. This plot is therefore rebuilt directly
                # from a dense threshold table in main().
                col=f"fraction_deltaR_lt_{p:.3f}"
                ok=g[g[col]>=cov]
                qs.append(ok.Q2_weighted_mean_GeV2.max() if len(ok) else np.nan)
            ax.plot(precisions,qs,marker="o",ms=3,label=f"remaining ×{f:g}")
        ax.set_title(f"≥{int(100*cov)}% of $(x_B,-t)$ area")
        ax.set_xlabel(r"Target $\delta(\sigma_L/\sigma_T)$")
        ax.grid(alpha=.25)
    axes[0].set_ylabel(r"Highest reachable $Q^2$ (GeV$^2$)")
    axes[0].legend(fontsize=8)
    fig.suptitle(r"$Q^2$ reach versus desired $L/T$ precision",y=.99)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220)
    plt.close(fig)


def dense_weighted_table(d, precisions):
    rows=[]
    for (f,iq),g in d.groupby(["future_luminosity_multiplier","Q2_bin"],sort=True):
        w=g.phase_area_xB_t.to_numpy(float)
        qmean=float(np.average(g.Q2_GeV2.to_numpy(float),weights=w))
        row={"future_multiplier":float(f),"Q2_bin":int(iq),
             "Q2_weighted_mean_GeV2":qmean,"phase_area":float(w.sum())}
        for p in precisions:
            row[f"fraction_deltaR_lt_{p:.3f}"]=weighted_fraction(
                g,g.delta_R_L_over_T.to_numpy(float)<p)
        rows.append(row)
    return pd.DataFrame(rows)


def plot_assumed_ratio_sensitivity(d, outfile):
    # Significance = R_true / deltaR. Weight by actual xB,-t phase-space area.
    ratios=np.linspace(.02,.60,60)
    fs=sorted(d.future_luminosity_multiplier.unique())
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.2),sharey=True)
    for ax,nsig in zip(axes,[1,2,3]):
        for f in fs:
            g=d[d.future_luminosity_multiplier==f]
            w=g.phase_area_xB_t.to_numpy(float)
            dr=g.delta_R_L_over_T.to_numpy(float)
            frac=[]
            for r in ratios:
                good=np.isfinite(w)&(w>0)&np.isfinite(dr)
                frac.append(np.sum(w[good & (r/dr>=nsig)])/np.sum(w[good]))
            ax.plot(ratios,100*np.asarray(frac),label=f"remaining ×{f:g}")
        ax.set_title(rf"$\sigma_L>0$ sensitivity ≥ {nsig}$\sigma$")
        ax.set_xlabel(r"Assumed true $\sigma_L/\sigma_T$")
        ax.grid(alpha=.25)
    axes[0].set_ylabel("Accessible phase-space coverage (%)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Nonzero-longitudinal sensitivity versus explicit assumed strength",y=.99)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220)
    plt.close(fig)


def plot_gk_assumption(d, outfile):
    # Explicitly show the assumed sigmaL, then its significance. The model name
    # is provenance, not the logical premise of the sensitivity statement.
    fs=sorted(d.future_luminosity_multiplier.unique())
    fshow=max(fs)
    g=d[d.future_luminosity_multiplier==fshow].copy()
    fig,ax=plt.subplots(figsize=(7.5,5.0))
    sc=ax.scatter(g.Q2_GeV2,g.sigma_L_GK,
                  c=np.abs(g.sigma_L_GK/g.delta_sigma_L),s=28)
    cb=fig.colorbar(sc,ax=ax)
    cb.set_label(r"$|\sigma_L|/\delta\sigma_L$")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax.set_ylabel(r"Assumed $\sigma_L$ (nb/GeV$^2$)")
    ax.set_title(fr"Illustrative longitudinal strength and sensitivity, remaining ×{fshow:g}")
    ax.text(.02,.98,"Central values: PARTONS/GK implementation\nshown explicitly; used only as an illustrative assumption",
            transform=ax.transAxes,va="top",fontsize=9)
    ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220)
    plt.close(fig)


def plot_synthetic_rosenbluth(d, outfile):
    # Pick low/mid/high Q2 representative cells with good epsilon leverage.
    base=d[d.future_luminosity_multiplier==0].copy()
    qs=np.quantile(base.Q2_GeV2,[.15,.50,.85])
    picks=[]
    for q in qs:
        h=base.iloc[(base.Q2_GeV2-q).abs().argsort()[:12]].copy()
        picks.append(h.iloc[h.delta_epsilon.abs().argmax()])
    fs=[0,1,5,10]
    fs=[f for f in fs if f in set(d.future_luminosity_multiplier)]
    Rtrue=.20
    Ttrue=100.0

    fig,axes=plt.subplots(1,3,figsize=(13.2,4.2),sharey=True)
    for ax,r in zip(axes,picks):
        er,ek=float(r.epsilon_rga),float(r.epsilon_rgk)
        Ltrue=Rtrue*Ttrue
        xx=np.linspace(min(er,ek)-.03,max(er,ek)+.03,100)
        ax.plot(xx,Ttrue+xx*Ltrue)
        for f in fs:
            rr=d[(d.point_id==r.point_id)&(d.future_luminosity_multiplier==f)].iloc[0]
            # For a clean visual Rosenbluth example use the separated-L uncertainty
            # to show the slope band rather than fake measured U central values.
            dl=float(rr.delta_sigma_L)
            ax.fill_between(xx,Ttrue+xx*(Ltrue-dl),Ttrue+xx*(Ltrue+dl),
                            alpha=.10,label=f"remaining ×{f:g}")
        ax.scatter([ek,er],[Ttrue+ek*Ltrue,Ttrue+er*Ltrue],marker="o",zorder=5)
        ax.set_xlabel(r"$\epsilon$")
        ax.set_title(fr"$Q^2\approx{r.Q2_GeV2:.2f}$ GeV$^2$, $x_B\approx{r.xB:.2f}$")
        ax.grid(alpha=.2)
    axes[0].set_ylabel(r"Synthetic $\sigma_U=\sigma_T+\epsilon\sigma_L$ (nb/GeV$^2$)")
    axes[0].legend(fontsize=8)
    fig.suptitle(r"Illustrative Rosenbluth slope: assumed $\sigma_T=100$, $\sigma_L/\sigma_T=0.20$"
                 "\nCentral values are synthetic; CLAS12 measurements are not shown",y=1.02)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220,bbox_inches="tight")
    plt.close(fig)


def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser()
    ap.add_argument("--reach",type=Path,default=here/"output/pi0_gk_measurement_reach/tables/01_measurement_reach_by_point.csv")
    ap.add_argument("--common",type=Path,default=here/"output/pi0_gk_stage2/tables/03_common_rosenbluth_model_points.csv")
    ap.add_argument("--output",type=Path,default=here/"output/pi0_workshop_reach")
    a=ap.parse_args()

    points=pd.read_csv(a.reach)
    common=pd.read_csv(a.common)
    d=attach_phase_space(points,common)
    a.output.mkdir(parents=True,exist_ok=True)

    thresholds=[.50,.30,.20,.10]
    tab=make_weighted_table(d,thresholds)
    tab.to_csv(a.output/"phase_space_weighted_reach.csv",index=False)
    plot_fraction_q2(tab,thresholds,a.output/"01_phase_space_fraction_vs_Q2.png")

    precisions=np.linspace(.05,.50,19)
    dense=dense_weighted_table(d,precisions)
    plot_reachable_q2(dense,a.output/"02_reachable_Q2_vs_precision.png")
    plot_assumed_ratio_sensitivity(d,a.output/"03_assumed_L_over_T_sensitivity.png")
    plot_gk_assumption(d,a.output/"04_GK_assumed_sigmaL_and_sensitivity.png")
    plot_synthetic_rosenbluth(d,a.output/"05_synthetic_rosenbluth_examples.png")

    print("\nWrote workshop-safe plots (no measured CLAS12 central values):")
    for p in sorted(a.output.glob("*.png")):
        print(" ",p)
    print(" ",a.output/"phase_space_weighted_reach.csv")


if __name__=="__main__":
    main()
