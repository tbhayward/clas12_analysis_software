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


def selected_scenarios(d):
    wanted = [1.0, 5.0, 10.0]
    have = set(d.future_luminosity_multiplier.astype(float))
    missing = [f for f in wanted if f not in have]
    if missing:
        raise RuntimeError(f"Reach table is missing requested future multipliers: {missing}")
    return wanted


def decompose_variance(d, value_col, out_stat_col):
    """
    Infer the pure-statistical variance component point-by-point from the
    scenario dependence

        delta^2(f) = A/(1+f) + B,

    where A is the recorded-data statistical variance and B is the
    statistics-independent conservative floor/nuisance contribution.

    This avoids silently treating the floor-inclusive uncertainty as
    statistical reach. Negative fitted intercepts/slopes from numerical noise
    are clipped to zero.
    """
    rows = []
    for pid, g in d.groupby("point_id", sort=False):
        gg = g[np.isfinite(g[value_col])].copy()
        x = 1.0/(1.0 + gg.future_luminosity_multiplier.to_numpy(float))
        y = gg[value_col].to_numpy(float)**2
        if len(gg) < 3:
            continue
        A = np.column_stack([x, np.ones(len(x))])
        a, b = np.linalg.lstsq(A, y, rcond=None)[0]
        a = max(float(a), 0.0)
        b = max(float(b), 0.0)
        for idx, f in zip(gg.index, gg.future_luminosity_multiplier.to_numpy(float)):
            rows.append((idx, math.sqrt(a/(1.0+f)), math.sqrt(b)))
        #endfor
    #endfor
    stat = pd.Series(np.nan, index=d.index, dtype=float)
    floor = pd.Series(np.nan, index=d.index, dtype=float)
    for idx, sv, fv in rows:
        stat.loc[idx] = sv
        floor.loc[idx] = fv
    #endfor
    d[out_stat_col] = stat
    d[out_stat_col + "_floor_equiv"] = floor
    return d


def make_weighted_table(d, thresholds, uncertainty_col):
    rows=[]
    for (f, iq), g in d.groupby(["future_luminosity_multiplier","Q2_bin"], sort=True):
        w=g.phase_area_xB_t.to_numpy(float)
        q=g.Q2_GeV2.to_numpy(float)
        qmean=float(np.average(q,weights=w))
        row=dict(future_multiplier=float(f), Q2_bin=int(iq),
                 Q2_weighted_mean_GeV2=qmean,
                 phase_area=float(np.sum(w)))
        u=g[uncertainty_col].to_numpy(float)
        for th in thresholds:
            row[f"fraction_deltaR_lt_{th:g}"]=weighted_fraction(g, u<th)
        #endfor
        rows.append(row)
    #endfor
    return pd.DataFrame(rows)


def plot_fraction_q2(tab, thresholds, outfile, subtitle):
    fs=[f for f in [1.0,5.0,10.0] if f in set(tab.future_multiplier)]
    fig, axes=plt.subplots(len(thresholds),1,figsize=(7.6,2.35*len(thresholds)),
                           sharex=True,sharey=True)
    axes=np.atleast_1d(axes)
    for ax,th in zip(axes,thresholds):
        for f in fs:
            g=tab[tab.future_multiplier==f].sort_values("Q2_weighted_mean_GeV2")
            ax.plot(g.Q2_weighted_mean_GeV2,
                    100*g[f"fraction_deltaR_lt_{th:g}"],
                    marker="o",label=f"remaining ×{f:g}")
        #endfor
        ax.set_ylabel("Coverage (%)")
        ax.set_title(rf"$\delta(\sigma_L/\sigma_T)<{th:g}$")
        ax.set_ylim(-3,103)
        ax.grid(alpha=.25)
    #endfor
    axes[-1].set_xlabel(r"$Q^2$ (GeV$^2$)")
    axes[0].legend(ncol=3,fontsize=8)
    fig.suptitle(r"Phase-space coverage for projected $1\sigma$ absolute uncertainty on $\sigma_L/\sigma_T$"+ "\n"+subtitle,y=.995)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220)
    plt.close(fig)


def dense_coverage(d, uncertainty_col, precisions):
    rows=[]
    for (f,iq),g in d.groupby(["future_luminosity_multiplier","Q2_bin"],sort=True):
        w=g.phase_area_xB_t.to_numpy(float)
        qmean=float(np.average(g.Q2_GeV2.to_numpy(float),weights=w))
        u=g[uncertainty_col].to_numpy(float)
        for p in precisions:
            rows.append(dict(future_multiplier=float(f),Q2_bin=int(iq),
                             Q2_weighted_mean_GeV2=qmean,target_precision=float(p),
                             coverage=weighted_fraction(g,u<p)))
        #endfor
    #endfor
    return pd.DataFrame(rows)


def plot_coverage_heatmaps(dense, outfile, subtitle):
    fs=[1.0,5.0,10.0]
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.2),sharex=True,sharey=True)
    im=None
    for ax,f in zip(axes,fs):
        g=dense[dense.future_multiplier==f]
        piv=g.pivot(index="Q2_weighted_mean_GeV2",columns="target_precision",values="coverage")
        x=piv.columns.to_numpy(float)
        y=piv.index.to_numpy(float)
        z=100*piv.to_numpy(float)
        im=ax.pcolormesh(x,y,z,shading="nearest",vmin=0,vmax=100)
        ax.set_title(f"remaining ×{f:g}")
        ax.set_xlabel(r"Target $\delta(\sigma_L/\sigma_T)$")
    #endfor
    axes[0].set_ylabel(r"$Q^2$ (GeV$^2$)")
    # Reserve a dedicated colorbar axis well to the right of the panels.
    # Do not rely on a small `pad`, which tends to crowd the third panel.
    fig.subplots_adjust(left=.08,right=.86,bottom=.15,top=.80,wspace=.08)
    cax=fig.add_axes([.895,.17,.018,.60])
    cb=fig.colorbar(im,cax=cax)
    cb.set_label(r"$(x_B,-t)$ phase-space coverage (%)",labelpad=12)
    fig.suptitle(r"Projected $1\sigma$ absolute uncertainty on $\sigma_L/\sigma_T$"+ "\n"+subtitle,y=.99)
    fig.savefig(outfile,dpi=220)
    plt.close(fig)


def plot_assumed_ratio_sensitivity(d, outfile, subtitle):
    """
    Nonzero-L significance is sigmaL_assumed/delta_sigmaL, not
    R/delta(R).  To make the premise explicit, use the stored GK sigmaT only
    as a cross-section scale:
        sigmaL_assumed = R_true * sigmaT_GK.
    The horizontal axis remains the assumed ratio, so the model dependence is
    visible and replaceable.
    """
    if "sigma_T_GK" not in d.columns:
        raise RuntimeError("Reach table lacks sigma_T_GK needed to set the explicit cross-section scale.")
    ratios=np.linspace(.05,.60,56)
    fs=selected_scenarios(d)
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.2),sharey=True)
    for ax,nsig in zip(axes,[1,2,3]):
        for f in fs:
            g=d[d.future_luminosity_multiplier==f]
            w=g.phase_area_xB_t.to_numpy(float)
            st=g.sigma_T_GK.to_numpy(float)
            dl=g["delta_sigma_L_for_sensitivity"].to_numpy(float)
            frac=[]
            for r in ratios:
                sig=np.abs(r*st)/dl
                good=np.isfinite(w)&(w>0)&np.isfinite(sig)
                frac.append(np.sum(w[good & (sig>=nsig)])/np.sum(w[good]))
            #endfor
            ax.plot(ratios,100*np.asarray(frac),label=f"remaining ×{f:g}")
        #endfor
        ax.set_title(rf"$\sigma_L>0$ sensitivity ≥ {nsig}$\sigma$")
        ax.set_xlabel(r"Assumed true $\sigma_L/\sigma_T$")
        ax.grid(alpha=.25)
    #endfor
    axes[0].set_ylabel("Accessible phase-space coverage (%)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Nonzero-longitudinal sensitivity versus explicit assumed strength\n"
                 +subtitle+
                 "\nCross-section scale from PARTONS/GK $\\sigma_T$; ratio is scanned explicitly",y=1.04)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220,bbox_inches="tight")
    plt.close(fig)


def plot_gk_assumption(d, outfile):
    fs=selected_scenarios(d)
    fig,axes=plt.subplots(2,1,figsize=(7.6,7.2),sharex=True)
    # Explicit model premise.
    g=d[d.future_luminosity_multiplier==fs[-1]].copy()
    axes[0].scatter(g.Q2_GeV2,g.sigma_L_GK,s=22)
    axes[0].set_ylabel(r"Assumed $\sigma_L$ (nb/GeV$^2$)")
    axes[0].set_title("Illustrative PARTONS/GK longitudinal strength")
    axes[0].grid(alpha=.2)

    for f in fs:
        g=d[d.future_luminosity_multiplier==f].copy()
        sig=np.abs(g.sigma_L_GK.to_numpy(float)/g.delta_sigma_L.to_numpy(float))
        axes[1].scatter(g.Q2_GeV2,sig,s=18,label=f"remaining ×{f:g}")
    #endfor
    axes[1].axhline(1,ls="--",lw=1)
    axes[1].axhline(2,ls="--",lw=1)
    axes[1].axhline(3,ls="--",lw=1)
    axes[1].set_xlabel(r"$Q^2$ (GeV$^2$)")
    axes[1].set_ylabel(r"$|\sigma_L|/\delta\sigma_L$")
    axes[1].set_title("Sensitivity if the longitudinal response has the magnitude shown above")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220)
    plt.close(fig)


def find_cov_tl_col(d):
    return pick(d,"cov_T_L","cov_sigma_T_sigma_L","cov_sigmaT_sigmaL",
                "cov_TL","cov_T_L_GK")


def plot_synthetic_rosenbluth(d, outfile):
    """
    Synthetic central values only.  The confidence band uses the complete
    T/L covariance:
      Var[U(eps)] = V_TT + 2 eps V_TL + eps^2 V_LL.
    The covariance is rescaled to the explicit synthetic sigmaT=100 scale
    using sigma_T_GK as the reference scale.
    """
    fs=selected_scenarios(d)
    base=d[d.future_luminosity_multiplier==fs[0]].copy()
    qs=np.quantile(base.Q2_GeV2,[.15,.50,.85])
    picks=[]
    for q in qs:
        h=base.iloc[(base.Q2_GeV2-q).abs().argsort()[:12]].copy()
        picks.append(h.iloc[h.delta_epsilon.abs().argmax()])
    #endfor

    covcol=find_cov_tl_col(d)
    if covcol is None:
        print("WARNING: no T-L covariance column found; skipping synthetic Rosenbluth figure.")
        return

    dtcol=pick(d,"delta_sigma_T","delta_T")
    dlcol=pick(d,"delta_sigma_L","delta_L")
    if dtcol is None or dlcol is None:
        print("WARNING: missing delta_sigma_T/L columns; skipping synthetic Rosenbluth figure.")
        return

    Rtrue=.20
    Ttrue=100.0
    Ltrue=Rtrue*Ttrue
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.2),sharey=True)

    for ax,r in zip(axes,picks):
        er,ek=float(r.epsilon_rga),float(r.epsilon_rgk)
        xx=np.linspace(min(er,ek)-.04,max(er,ek)+.04,160)
        truth=Ttrue+xx*Ltrue
        ax.plot(xx,truth,lw=1.8,label="synthetic truth")
        ax.scatter([ek,er],[Ttrue+ek*Ltrue,Ttrue+er*Ltrue],marker="o",zorder=6)

        # Show only 1x, 5x, 10x and use the full covariance for each band.
        # Deliberately use strongly separated hues for the nested confidence
        # bands; pale default-cycle fills are too difficult to distinguish.
        band_colors={1.0:"tab:orange", 5.0:"tab:green", 10.0:"tab:purple"}
        band_alpha={1.0:.20, 5.0:.24, 10.0:.28}
        for f in fs:
            rr=d[(d.point_id==r.point_id)&
                 (d.future_luminosity_multiplier==f)].iloc[0]
            refT=abs(float(rr.sigma_T_GK))
            scale=Ttrue/refT if refT>0 else np.nan
            vtt=(float(rr[dtcol])*scale)**2
            vll=(float(rr[dlcol])*scale)**2
            vtl=float(rr[covcol])*scale**2
            vu=vtt+2*xx*vtl+(xx**2)*vll
            du=np.sqrt(np.maximum(vu,0.0))
            ax.fill_between(
                xx,truth-du,truth+du,
                color=band_colors[f],alpha=band_alpha[f],
                label=f"remaining ×{f:g}"
            )
        #endfor

        ax.set_xlabel(r"$\epsilon$")
        ax.set_title(fr"$Q^2\approx{r.Q2_GeV2:.2f}$ GeV$^2$, $x_B\approx{r.xB:.2f}$")
        ax.grid(alpha=.2)
    #endfor
    axes[0].set_ylabel(r"Synthetic $\sigma_U=\sigma_T+\epsilon\sigma_L$ (nb/GeV$^2$)")
    axes[0].legend(fontsize=8)
    fig.suptitle(r"Illustrative Rosenbluth separation: $\sigma_T=100$ nb/GeV$^2$, "
                 r"$\sigma_L/\sigma_T=0.20$"
                 "\nSynthetic central values only; bands use projected full $T$-$L$ covariance",
                 y=1.02)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220,bbox_inches="tight")
    plt.close(fig)

def plot_normalization_scan_10x(d, thresholds, outfile):
    """At fixed remaining x10, show coverage versus relative normalization uncertainty."""
    norms=sorted(d.relative_normalization_uncertainty.unique())
    fig,axes=plt.subplots(len(thresholds),1,figsize=(7.6,2.35*len(thresholds)),
                          sharex=True,sharey=True)
    axes=np.atleast_1d(axes)
    for ax,th in zip(axes,thresholds):
        for n in norms:
            g=d[d.relative_normalization_uncertainty==n]
            vals=[]
            qs=[]
            for iq,h in g.groupby("Q2_bin",sort=True):
                w=h.phase_area_xB_t.to_numpy(float)
                qs.append(float(np.average(h.Q2_GeV2.to_numpy(float),weights=w)))
                vals.append(100*weighted_fraction(
                    h,h.delta_R_L_over_T.to_numpy(float)<th))
            #endfor
            ax.plot(qs,vals,marker="o",label=f"norm. {100*n:.0f}%")
        #endfor
        ax.set_ylabel("Coverage (%)")
        ax.set_title(rf"$\delta_{{1\sigma}}(\sigma_L/\sigma_T)<{th:g}$")
        ax.set_ylim(-3,103)
        ax.grid(alpha=.25)
    #endfor
    axes[-1].set_xlabel(r"$Q^2$ (GeV$^2$)")
    axes[0].legend(ncol=3,fontsize=8)
    fig.suptitle(
        r"Remaining $\times10$: sensitivity to RGA/RGK relative normalization"
        "\nPer-campaign systematic floor held fixed",y=.995)
    fig.tight_layout()
    fig.savefig(outfile,dpi=220)
    plt.close(fig)


def plot_normalization_summary_10x(d, thresholds, outfile):
    """Global phase-space coverage versus relative-normalization uncertainty."""
    norms=sorted(d.relative_normalization_uncertainty.unique())
    fig,ax=plt.subplots(figsize=(7.2,4.8))
    for th in thresholds:
        yy=[]
        for n in norms:
            g=d[d.relative_normalization_uncertainty==n]
            w=g.phase_area_xB_t.to_numpy(float)
            good=np.isfinite(w)&(w>0)
            yy.append(100*np.sum(w[good & (g.delta_R_L_over_T.to_numpy(float)<th)])/np.sum(w[good]))
        #endfor
        ax.plot(100*np.asarray(norms),yy,marker="o",
                label=rf"$\delta_{{1\sigma}}(L/T)<{th:g}$")
    #endfor
    ax.set_xlabel("RGA/RGK relative-normalization uncertainty (%)")
    ax.set_ylabel(r"Accessible $(x_B,-t,Q^2)$ phase-space coverage (%)")
    ax.set_ylim(-3,103)
    ax.set_title(r"Remaining $\times10$: normalization-control requirement")
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile,dpi=220)
    plt.close(fig)

def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser()
    ap.add_argument("--reach",type=Path,
                    default=here/"output/pi0_gk_measurement_reach/tables/01_measurement_reach_by_point.csv")
    ap.add_argument("--common",type=Path,
                    default=here/"output/pi0_gk_stage2/tables/03_common_rosenbluth_model_points.csv")
    ap.add_argument("--stat-reach",type=Path,
                    default=here/"output/pi0_gk_measurement_reach/tables/07_measurement_reach_stat_only_by_point.csv")
    ap.add_argument("--norm-scan",type=Path,
                    default=here/"output/pi0_gk_measurement_reach/tables/08_measurement_reach_norm_scan_10x_by_point.csv")
    ap.add_argument("--output",type=Path,default=here/"output/pi0_workshop_reach")
    a=ap.parse_args()

    points=pd.read_csv(a.reach)
    stat_points=pd.read_csv(a.stat_reach)
    norm_points=pd.read_csv(a.norm_scan)
    common=pd.read_csv(a.common)

    d=attach_phase_space(points,common)
    ds=attach_phase_space(stat_points,common)
    dn=attach_phase_space(norm_points,common)
    selected_scenarios(d)
    selected_scenarios(ds)
    a.output.mkdir(parents=True,exist_ok=True)

    thresholds=[.30,.20,.10]
    conservative=make_weighted_table(d,thresholds,"delta_R_L_over_T")
    statistical=make_weighted_table(ds,thresholds,"delta_R_L_over_T")
    conservative.to_csv(a.output/"phase_space_weighted_reach_conservative.csv",index=False)
    statistical.to_csv(a.output/"phase_space_weighted_reach_stat_only.csv",index=False)

    plot_fraction_q2(
        statistical,thresholds,
        a.output/"01a_phase_space_fraction_vs_Q2_stat_only.png",
        "statistical component only")
    plot_fraction_q2(
        conservative,thresholds,
        a.output/"01b_phase_space_fraction_vs_Q2_conservative.png",
        "including current conservative floor / relative-normalization treatment")

    precisions=np.linspace(.05,.50,46)
    dense_stat=dense_coverage(ds,"delta_R_L_over_T",precisions)
    dense_cons=dense_coverage(d,"delta_R_L_over_T",precisions)
    dense_stat.to_csv(a.output/"precision_coverage_heatmap_stat_only.csv",index=False)
    dense_cons.to_csv(a.output/"precision_coverage_heatmap_conservative.csv",index=False)
    plot_coverage_heatmaps(
        dense_stat,a.output/"02a_precision_coverage_heatmap_stat_only.png",
        "statistical component only")
    plot_coverage_heatmaps(
        dense_cons,a.output/"02b_precision_coverage_heatmap_conservative.png",
        "including current conservative floor / relative-normalization treatment")

    # Make both versions of the nonzero-L sensitivity.  The significance is
    # sigmaL_assumed/delta_sigmaL; the assumed ratio is scanned explicitly.
    ds["delta_sigma_L_for_sensitivity"]=ds["delta_sigma_L"]
    plot_assumed_ratio_sensitivity(
        ds,a.output/"03a_assumed_L_over_T_sensitivity_stat_only.png",
        "statistical component only")
    d["delta_sigma_L_for_sensitivity"]=d["delta_sigma_L"]
    plot_assumed_ratio_sensitivity(
        d,a.output/"03b_assumed_L_over_T_sensitivity_conservative.png",
        "including current conservative floor / relative-normalization treatment")

    plot_gk_assumption(d,a.output/"04_GK_assumed_sigmaL_and_sensitivity.png")
    plot_synthetic_rosenbluth(d,a.output/"05_synthetic_rosenbluth_examples.png")

    # Normalization-control study at remaining x10.
    plot_normalization_scan_10x(
        dn,thresholds,a.output/"06a_normalization_scan_10x_vs_Q2.png")
    plot_normalization_summary_10x(
        dn,thresholds,a.output/"06b_normalization_scan_10x_global.png")

    print("\nWrote workshop-safe outputs (no measured CLAS12 central values).")
    print("Displayed scenarios: remaining ×1, ×5, ×10.")
    print("Stat-only curves are generated directly with systematic floor=0 and relative normalization=0.")
    print("Normalization scan holds the per-campaign systematic floor fixed and varies only RGA/RGK relative normalization.")
    for p in sorted(a.output.glob("*.png")):
        print(" ",p)
    #endfor


if __name__=="__main__":
    main()
