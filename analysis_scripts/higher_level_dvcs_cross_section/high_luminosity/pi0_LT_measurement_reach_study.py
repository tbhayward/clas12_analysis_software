#!/usr/bin/env python3
"""
Standalone pi0 Rosenbluth L/T measurement-reach study (regenerated v2).

Purpose
-------
Reorient the Stage-3 RGA+RGK projection around experimentally useful precision:
  * absolute and relative uncertainties on T, L, LT and TT;
  * projected uncertainty on R_L/T = sigma_L / sigma_T;
  * expected 95% sensitivity to an upper bound on |L/T| if L is consistent with zero;
  * counts of bins that cross useful L/T precision thresholds as luminosity increases;
  * "new capability" counts: bins that are not useful with recorded data but become
    useful after nominal/high-luminosity completion.

This intentionally retains the existing sigma_L significance information as a
complementary metric rather than replacing it.

The script imports the validated Stage-3 machinery from
prepare_pi0_gk_stage3_projection.py. It does not rerun PARTONS.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def args():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--stage3-script", type=Path,
                   default=here/"prepare_pi0_gk_stage3_projection.py")
    p.add_argument("--stage1", type=Path, default=here/"output"/"pi0_gk_stage1")
    p.add_argument("--stage2", type=Path, default=here/"output"/"pi0_gk_stage2")
    p.add_argument("--gk-results", type=Path,
                   default=here/"output"/"pi0_gk_stage3"/"model_queries"/
                           "02_gk_structure_function_results.csv")
    p.add_argument("--output", type=Path,
                   default=here/"output"/"pi0_gk_measurement_reach")
    p.add_argument("--future-lumi-scan", type=str, default="0,1,3,5",
                   help="Multiplier f applied to the remaining beam time.")
    p.add_argument("--rga-recorded", type=float, default=1.5)
    p.add_argument("--rga-remaining", type=float, default=1.5)
    p.add_argument("--rgk-recorded", type=float, default=9.0)
    p.add_argument("--rgk-remaining", type=float, default=9.0)
    return p.parse_args()


def load_stage3(path):
    spec = importlib.util.spec_from_file_location("pi0_stage3", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import Stage-3 script: {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def invvar_combine(values, uncertainties):
    values = np.asarray(values, float)
    uncertainties = np.asarray(uncertainties, float)
    good = np.isfinite(values) & np.isfinite(uncertainties) & (uncertainties > 0)
    if not np.any(good):
        return np.nan, np.nan
    w = 1.0 / uncertainties[good]**2
    return float(np.sum(w*values[good])/np.sum(w)), float(np.sqrt(1.0/np.sum(w)))


def build_reach_table(model, hfits, lt):
    """Build one row per common (Q2,xB,t) cell."""
    rows = []

    for r in model.itertuples(index=False):
        lr = lt[lt.point_id == r.point_id].iloc[0]
        hf = hfits[hfits.point_id == r.point_id]

        tt, dtt = invvar_combine(hf.sigma_TT_fit, hf.sigma_TT_unc)
        lti, dlti = invvar_combine(hf.sigma_LT_fit, hf.sigma_LT_unc)

        T = float(lr.sigma_T_proj)
        L = float(lr.sigma_L_proj)
        dT = float(lr.sigma_T_unc)
        dL = float(lr.sigma_L_unc)
        R = float(lr.R_L_over_T_proj)
        dR = float(lr.R_L_over_T_unc)

        # If the separated L result is statistically consistent with zero,
        # this is the expected Gaussian two-sided 95% sensitivity to |L/T|.
        # It is deliberately evaluated about L=0, so it does not depend on
        # GK's tiny central sigma_L.
        r95_zero = 1.96*dL/abs(T) if T != 0 else np.nan

        rows.append(dict(
            point_id=r.point_id,
            Q2_GeV2=float(r.Q2_common_GeV2),
            xB=float(r.xB_common),
            minus_t_GeV2=float(r.minus_t_common_GeV2),
            epsilon_rga=float(r.epsilon_rga),
            epsilon_rgk=float(r.epsilon_rgk),
            delta_epsilon=float(r.delta_epsilon),

            sigma_T_GK=float(r.sigma_T),
            sigma_L_GK=float(r.sigma_L),
            sigma_TT_GK=float(r.sigma_TT),
            sigma_LT_GK=float(r.sigma_LT),

            sigma_T_proj=T,
            delta_sigma_T=dT,
            rel_delta_sigma_T=dT/abs(T) if T != 0 else np.nan,

            sigma_L_proj=L,
            delta_sigma_L=dL,
            rel_delta_sigma_L=dL/abs(L) if L != 0 else np.nan,
            sigma_L_significance=abs(L)/dL if dL > 0 else np.nan,

            sigma_TT_proj=tt,
            delta_sigma_TT=dtt,
            rel_delta_sigma_TT=dtt/abs(tt) if tt != 0 else np.nan,
            sigma_TT_significance=abs(tt)/dtt if dtt > 0 else np.nan,

            sigma_LT_proj=lti,
            delta_sigma_LT=dlti,
            rel_delta_sigma_LT=dlti/abs(lti) if lti != 0 else np.nan,
            sigma_LT_significance=abs(lti)/dlti if dlti > 0 else np.nan,

            R_L_over_T_GK=float(r.sigma_L/r.sigma_T) if r.sigma_T != 0 else np.nan,
            R_L_over_T_proj=R,
            delta_R_L_over_T=dR,
            expected_95pct_abs_L_over_T_limit_if_L_zero=r95_zero,
        ))

    return pd.DataFrame(rows)


def finite_median(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else np.nan


def summarize(reach, f, rga, rgk):
    dR = reach.delta_R_L_over_T.to_numpy(float)
    r95 = reach.expected_95pct_abs_L_over_T_limit_if_L_zero.to_numpy(float)
    sigL = reach.sigma_L_significance.to_numpy(float)

    row = dict(
        future_luminosity_multiplier=f,
        rga_final_factor=rga,
        rgk_final_factor=rgk,
        n_points=len(reach),

        median_delta_sigma_T=finite_median(reach.delta_sigma_T),
        median_delta_sigma_L=finite_median(reach.delta_sigma_L),
        median_delta_sigma_LT=finite_median(reach.delta_sigma_LT),
        median_delta_sigma_TT=finite_median(reach.delta_sigma_TT),

        median_delta_R_L_over_T=finite_median(dR),
        n_delta_R_lt_0p20=int(np.sum(dR < 0.20)),
        n_delta_R_lt_0p10=int(np.sum(dR < 0.10)),
        n_delta_R_lt_0p05=int(np.sum(dR < 0.05)),

        median_expected_95pct_abs_L_over_T_limit_if_L_zero=finite_median(r95),
        n_expected_95pct_limit_lt_0p20=int(np.sum(r95 < 0.20)),
        n_expected_95pct_limit_lt_0p10=int(np.sum(r95 < 0.10)),
        n_expected_95pct_limit_lt_0p05=int(np.sum(r95 < 0.05)),

        n_sigmaL_ge_1=int(np.sum(sigL >= 1.0)),
        n_sigmaL_ge_2=int(np.sum(sigL >= 2.0)),
        n_sigmaL_ge_3=int(np.sum(sigL >= 3.0)),
    )
    return row


def add_new_capability_flags(all_points):
    """Mark bins whose L/T constraint crosses a threshold relative to f=0."""
    base = all_points[all_points.future_luminosity_multiplier == 0][
        ["point_id", "expected_95pct_abs_L_over_T_limit_if_L_zero"]
    ].rename(columns={
        "expected_95pct_abs_L_over_T_limit_if_L_zero": "r95_recorded"
    })
    out = all_points.merge(base, on="point_id", how="left", validate="many_to_one")

    for thr, tag in [(0.20, "0p20"), (0.10, "0p10"), (0.05, "0p05")]:
        out[f"newly_constrains_abs_L_over_T_below_{tag}"] = (
            (out.r95_recorded >= thr) &
            (out.expected_95pct_abs_L_over_T_limit_if_L_zero < thr)
        )
    return out


def plot_constraint_counts(summary, outfile):
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    x = summary.future_luminosity_multiplier
    ax.plot(x, summary.n_expected_95pct_limit_lt_0p20, marker="o",
            label=r"95% sensitivity: $|L/T|<0.20$")
    ax.plot(x, summary.n_expected_95pct_limit_lt_0p10, marker="o",
            label=r"95% sensitivity: $|L/T|<0.10$")
    ax.plot(x, summary.n_expected_95pct_limit_lt_0p05, marker="o",
            label=r"95% sensitivity: $|L/T|<0.05$")
    ax.axvline(1.0, ls="--", label="Nominal completion")
    ax.set_xlabel("Luminosity multiplier for remaining beam time")
    ax.set_ylabel("Number of common kinematic bins")
    ax.set_title(r"Projected ability to constrain longitudinal fraction")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=180)
    plt.close(fig)


def plot_ratio_precision_map(points, f, outfile):
    g = points[points.future_luminosity_multiplier == f].copy()
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    sc = ax.scatter(g.Q2_GeV2, g.minus_t_GeV2,
                    c=g.expected_95pct_abs_L_over_T_limit_if_L_zero, s=34)
    fig.colorbar(sc, ax=ax,
                 label=r"Expected 95% sensitivity to $|L/T|$ if $L=0$")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax.set_ylabel(r"$-t$ (GeV$^2$)")
    ax.set_title(fr"Longitudinal-fraction reach, $f={f:g}$")
    fig.tight_layout()
    fig.savefig(outfile, dpi=180)
    plt.close(fig)


def plot_absolute_uncertainties(points, f, outfile):
    g = points[points.future_luminosity_multiplier == f].copy()
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    ax.scatter(g.Q2_GeV2, g.delta_sigma_T, s=24, label=r"$\delta\sigma_T$")
    ax.scatter(g.Q2_GeV2, g.delta_sigma_L, s=24, label=r"$\delta\sigma_L$")
    ax.scatter(g.Q2_GeV2, g.delta_sigma_LT, s=24, label=r"$\delta\sigma_{LT}$")
    ax.scatter(g.Q2_GeV2, g.delta_sigma_TT, s=24, label=r"$\delta\sigma_{TT}$")
    ax.set_yscale("log")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax.set_ylabel(r"Projected uncertainty (nb/GeV$^2$)")
    ax.set_title(fr"Absolute structure-function precision, $f={f:g}$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=180)
    plt.close(fig)


def main():
    a = args()
    out = a.output.resolve()
    tabs = out/"tables"
    figs = out/"figures"
    tabs.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    try:
        fs = [float(x.strip()) for x in a.future_lumi_scan.split(",") if x.strip()]
    except ValueError as exc:
        raise RuntimeError(f"Could not parse --future-lumi-scan={a.future_lumi_scan!r}") from exc
    if not fs or any(f < 0 for f in fs):
        raise RuntimeError("--future-lumi-scan must contain non-negative values")
    if 0.0 not in fs:
        raise RuntimeError("Include f=0 so the script can define current/recorded capability.")

    st3 = load_stage3(a.stage3_script.resolve())
    common = pd.read_csv(a.stage2.resolve()/"tables"/"03_common_rosenbluth_model_points.csv")
    q = st3.build_query_points(common)
    model = st3.model_merge(q, a.gk_results.resolve())

    all_rows = []
    summaries = []

    for f in fs:
        rga = a.rga_recorded + f*a.rga_remaining
        rgk = a.rgk_recorded + f*a.rgk_remaining

        pseudo, hfits = st3.make_pseudodata(
            model, a.stage1.resolve(), a.stage2.resolve(), rga, rgk
        )
        lt = st3.lt_from_u(model, hfits)
        reach = build_reach_table(model, hfits, lt)
        reach.insert(0, "rgk_final_factor", rgk)
        reach.insert(0, "rga_final_factor", rga)
        reach.insert(0, "future_luminosity_multiplier", f)

        all_rows.append(reach)
        summaries.append(summarize(reach, f, rga, rgk))

    points = pd.concat(all_rows, ignore_index=True)
    points = add_new_capability_flags(points)
    summary = pd.DataFrame(summaries)

    points.to_csv(tabs/"01_measurement_reach_by_point.csv", index=False)
    summary.to_csv(tabs/"02_measurement_reach_summary.csv", index=False)

    unlocked = []
    for f in fs:
        if f == 0:
            continue
        g = points[points.future_luminosity_multiplier == f]
        unlocked.append(dict(
            future_luminosity_multiplier=f,
            rga_final_factor=float(g.rga_final_factor.iloc[0]),
            rgk_final_factor=float(g.rgk_final_factor.iloc[0]),
            newly_below_20pct_95CL=int(g.newly_constrains_abs_L_over_T_below_0p20.sum()),
            newly_below_10pct_95CL=int(g.newly_constrains_abs_L_over_T_below_0p10.sum()),
            newly_below_5pct_95CL=int(g.newly_constrains_abs_L_over_T_below_0p05.sum()),
        ))
    pd.DataFrame(unlocked).to_csv(tabs/"03_new_capability_vs_recorded.csv", index=False)

    plot_constraint_counts(summary, figs/"01_L_over_T_constraint_counts.png")
    for f in fs:
        tag=str(f).replace(".", "p")
        plot_ratio_precision_map(points, f, figs/f"02_L_over_T_reach_map_f{tag}.png")
    plot_absolute_uncertainties(points, 1.0 if 1.0 in fs else fs[0],
                                figs/"03_absolute_structure_function_uncertainties.png")

    cols = [
        "future_luminosity_multiplier", "rga_final_factor", "rgk_final_factor",
        "n_delta_R_lt_0p20", "n_delta_R_lt_0p10", "n_delta_R_lt_0p05",
        "n_expected_95pct_limit_lt_0p20",
        "n_expected_95pct_limit_lt_0p10",
        "n_expected_95pct_limit_lt_0p05",
        "n_sigmaL_ge_1", "n_sigmaL_ge_2", "n_sigmaL_ge_3",
    ]
    print("\n[pi0 measurement reach]")
    print("L/T limit columns are expected Gaussian 95% sensitivity if L is consistent with zero.")
    print("sigma_L significance columns retain the complementary nonzero-L question.")
    print(summary[cols].to_string(index=False))

    print("\n[capability unlocked relative to already-recorded data]")
    if unlocked:
        print(pd.DataFrame(unlocked).to_string(index=False))

    print("\nWrote:")
    print(f"  {tabs/'01_measurement_reach_by_point.csv'}")
    print(f"  {tabs/'02_measurement_reach_summary.csv'}")
    print(f"  {tabs/'03_new_capability_vs_recorded.csv'}")
    print(f"  {figs/'01_L_over_T_constraint_counts.png'}")
    print("  per-scenario L/T reach maps")
    print(f"  {figs/'03_absolute_structure_function_uncertainties.png'}")
    print("\nCaveat: this inherits the Stage-3 assumption that supplied fractional")
    print("uncertainties scale as 1/sqrt(exposure); irreducible systematic and finite-MC")
    print("floors are not yet separated.")


if __name__ == "__main__":
    main()
