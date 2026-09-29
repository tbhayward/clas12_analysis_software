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

The v3 extension explicitly targets the physics motivation for the high-luminosity
white paper: map the Q2 evolution of transverse dominance and identify kinematic
regions where additional luminosity creates a qualitatively new L/T constraint.
It also writes nearest-cell benchmark tables for the Hall-A Defurne (2016) true
Rosenbluth separation and Dlamini (2021) high-Q2 U/LT/TT/LT' measurements.

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
    p.add_argument("--future-lumi-scan", type=str, default="0,1,3,5,10",
                   help="Multiplier f applied to the remaining beam time; f=10 is the workshop high-luminosity baseline.")
    p.add_argument("--rga-recorded", type=float, default=1.5)
    p.add_argument("--rga-remaining", type=float, default=1.5)
    p.add_argument("--rgk-recorded", type=float, default=9.0)
    p.add_argument("--rgk-remaining", type=float, default=9.0)
    p.add_argument(
        "--central-values", choices=("gk", "data"), default="gk",
        help=("Central values used in the study. Default 'gk' is workshop-safe. "
              "'data' is INTERNAL ONLY and uses measured RGA/RGK cross sections."),
    )
    p.add_argument(
        "--internal-rga", type=Path,
        default=here/"import"/"fa18_rosenbluth_inputs_20260924T165948Z"/
                     "rga_10604"/"combined_reduced_cross_sections.csv",
        help="INTERNAL ONLY: measured RGA reduced-cross-section CSV.",
    )
    p.add_argument(
        "--internal-rgk", type=Path,
        default=here/"import"/"fa18_rosenbluth_inputs_20260924T165948Z"/
                     "rgk_6535"/"rgk6535_reduced_cross_sections.csv",
        help="INTERNAL ONLY: measured RGK reduced-cross-section CSV.",
    )
    p.add_argument(
        "--internal-output", type=Path,
        default=here/"output"/"pi0_LT_internal_data",
        help="Separate output directory used only by --central-values data.",
    )
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
    """Faceted (xB,-t) map of L/T reach, with one panel per Q2 setting."""
    g = points[points.future_luminosity_multiplier == f].copy()
    q2_values = sorted(g.Q2_GeV2.unique())

    ncols = 4
    nrows = int(np.ceil(len(q2_values) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(13.0, 3.8*nrows),
        sharex=True, sharey=True, squeeze=False
    )
    axes = axes.ravel()

    # Discrete physics-reach categories.  Lower values are better.
    categories = [
        (0.00, 0.05, r"$<5\%$"),
        (0.05, 0.10, r"$5$--$10\%$"),
        (0.10, 0.20, r"$10$--$20\%$"),
        (0.20, np.inf, r"$>20\%$"),
    ]
    cmap = plt.get_cmap("viridis")
    colors = [cmap(x) for x in (0.05, 0.35, 0.65, 0.92)]

    for ax, q2 in zip(axes, q2_values):
        q = g[np.isclose(g.Q2_GeV2, q2)]
        reach = q.expected_95pct_abs_L_over_T_limit_if_L_zero.to_numpy(float)

        for (lo, hi, label), color in zip(categories, colors):
            mask = (reach >= lo) & (reach < hi)
            if np.any(mask):
                ax.scatter(
                    q.loc[mask, "xB"], q.loc[mask, "minus_t_GeV2"],
                    s=55, color=color, edgecolor="black", linewidth=0.35,
                    label=label,
                )

        ax.set_title(fr"$Q^2={q2:g}$ GeV$^2$")
        ax.grid(alpha=0.18)

    for ax in axes[len(q2_values):]:
        ax.set_visible(False)

    for i, ax in enumerate(axes[:len(q2_values)]):
        if i // ncols == nrows - 1 or i + ncols >= len(q2_values):
            ax.set_xlabel(r"$x_B$")
        if i % ncols == 0:
            ax.set_ylabel(r"$-t$ (GeV$^2$)")

    # One common legend, in the desired best-to-worst order.
    handles = []
    labels = []
    for (lo, hi, label), color in zip(categories, colors):
        h = plt.Line2D([], [], linestyle="none", marker="o", markersize=7,
                       markerfacecolor=color, markeredgecolor="black",
                       markeredgewidth=0.35)
        handles.append(h)
        labels.append(label)

    fig.legend(
        handles, labels,
        title=r"Expected 95% sensitivity to $|L/T|$",
        loc="lower center", ncol=4, frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )
    fig.suptitle(
        fr"Projected $\pi^0$ longitudinal-fraction reach, $f={f:g}$",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.065, 1, 0.96))
    fig.savefig(outfile, dpi=200)
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



def q2_evolution_summary(points):
    """Summarize longitudinal-fraction reach separately at each Q2 setting."""
    rows = []
    for (f, q2), g in points.groupby(
        ["future_luminosity_multiplier", "Q2_GeV2"], sort=True
    ):
        r95 = g.expected_95pct_abs_L_over_T_limit_if_L_zero.to_numpy(float)
        dR = g.delta_R_L_over_T.to_numpy(float)
        rows.append(dict(
            future_luminosity_multiplier=float(f),
            Q2_GeV2=float(q2),
            n_bins=len(g),
            xB_min=float(g.xB.min()),
            xB_max=float(g.xB.max()),
            minus_t_min_GeV2=float(g.minus_t_GeV2.min()),
            minus_t_max_GeV2=float(g.minus_t_GeV2.max()),
            median_delta_R_L_over_T=finite_median(dR),
            median_expected_95pct_abs_L_over_T_limit_if_L_zero=finite_median(r95),
            n_expected_95pct_limit_lt_0p20=int(np.sum(r95 < 0.20)),
            n_expected_95pct_limit_lt_0p10=int(np.sum(r95 < 0.10)),
            n_expected_95pct_limit_lt_0p05=int(np.sum(r95 < 0.05)),
            n_sigmaL_ge_2=int(np.sum(g.sigma_L_significance.to_numpy(float) >= 2.0)),
            n_sigmaL_ge_3=int(np.sum(g.sigma_L_significance.to_numpy(float) >= 3.0)),
        ))
    return pd.DataFrame(rows)


def benchmark_regions(points):
    """
    Identify CLAS12 cells nearest the principal published Hall-A pi0 regions.

    These are comparison anchors, not claims of identical kinematics:
      * E07-007 / Defurne 2016: xB=0.36, Q2=1.50,1.75,2.00 GeV^2,
        true Rosenbluth T/L separation.
      * E12-06-114 / Dlamini 2021: xB=0.36,0.48,0.60 and Q2 roughly
        3.1--8.4 GeV^2, U/TT/LT/LT' but no T/L Rosenbluth separation.

    For every luminosity scenario we retain the nearest cells in (Q2,xB);
    t remains explicit so the user can judge the actual overlap.
    """
    anchors = [
        ("HallA_Defurne2016", 0.36, 1.50),
        ("HallA_Defurne2016", 0.36, 1.75),
        ("HallA_Defurne2016", 0.36, 2.00),
        ("HallA_Dlamini2021", 0.36, 3.11),
        ("HallA_Dlamini2021", 0.36, 3.57),
        ("HallA_Dlamini2021", 0.36, 4.44),
        ("HallA_Dlamini2021", 0.48, 2.67),
        ("HallA_Dlamini2021", 0.48, 4.06),
        ("HallA_Dlamini2021", 0.48, 5.16),
        ("HallA_Dlamini2021", 0.48, 6.56),
        ("HallA_Dlamini2021", 0.60, 5.49),
        ("HallA_Dlamini2021", 0.60, 8.31),
    ]

    rows = []
    for f, gf in points.groupby("future_luminosity_multiplier"):
        unique_qx = gf[["Q2_GeV2", "xB"]].drop_duplicates().copy()
        for experiment, xb0, q20 in anchors:
            # Dimensionless fractional distance prevents Q2 from trivially
            # dominating xB in the nearest-cell selection.
            d2 = ((unique_qx.Q2_GeV2-q20)/max(q20, 0.5))**2
            d2 += ((unique_qx.xB-xb0)/max(xb0, 0.1))**2
            best = unique_qx.loc[d2.idxmin()]
            sel = gf[
                np.isclose(gf.Q2_GeV2, best.Q2_GeV2) &
                np.isclose(gf.xB, best.xB)
            ].copy()
            for r in sel.itertuples(index=False):
                rows.append(dict(
                    experiment=experiment,
                    anchor_Q2_GeV2=q20,
                    anchor_xB=xb0,
                    future_luminosity_multiplier=float(f),
                    nearest_Q2_GeV2=float(r.Q2_GeV2),
                    nearest_xB=float(r.xB),
                    minus_t_GeV2=float(r.minus_t_GeV2),
                    delta_epsilon=float(r.delta_epsilon),
                    delta_sigma_T=float(r.delta_sigma_T),
                    delta_sigma_L=float(r.delta_sigma_L),
                    delta_sigma_LT=float(r.delta_sigma_LT),
                    delta_sigma_TT=float(r.delta_sigma_TT),
                    delta_R_L_over_T=float(r.delta_R_L_over_T),
                    expected_95pct_abs_L_over_T_limit_if_L_zero=float(
                        r.expected_95pct_abs_L_over_T_limit_if_L_zero
                    ),
                    sigma_L_significance=float(r.sigma_L_significance),
                ))
    return pd.DataFrame(rows)


def q2_unlock_summary(points):
    """Count genuinely new L/T constraints, relative to recorded data, versus Q2."""
    rows = []
    for (f, q2), g in points.groupby(
        ["future_luminosity_multiplier", "Q2_GeV2"], sort=True
    ):
        if f == 0:
            continue
        rows.append(dict(
            future_luminosity_multiplier=float(f),
            Q2_GeV2=float(q2),
            n_bins=len(g),
            newly_below_20pct_95CL=int(
                g.newly_constrains_abs_L_over_T_below_0p20.sum()
            ),
            newly_below_10pct_95CL=int(
                g.newly_constrains_abs_L_over_T_below_0p10.sum()
            ),
            newly_below_5pct_95CL=int(
                g.newly_constrains_abs_L_over_T_below_0p05.sum()
            ),
        ))
    return pd.DataFrame(rows)


def plot_q2_evolution(q2sum, outfile):
    """Show how far in Q2 each luminosity scenario can test transverse dominance."""
    fig, ax = plt.subplots(figsize=(8.0, 5.6))

    label_map = {
        0.0: "Recorded",
        1.0: "Nominal completion",
        3.0: r"$f=3$",
        5.0: r"$f=5$",
        10.0: r"$f=10$ (workshop baseline)",
    }

    for f, g in q2sum.groupby("future_luminosity_multiplier", sort=True):
        g = g.sort_values("Q2_GeV2")
        ax.plot(
            g.Q2_GeV2,
            100.0 * g.median_expected_95pct_abs_L_over_T_limit_if_L_zero,
            marker="o",
            linewidth=1.8,
            label=label_map.get(float(f), fr"$f={f:g}$"),
        )

    ax.axhline(20.0, ls=":", linewidth=1.5)
    ax.axhline(10.0, ls="--", linewidth=1.5)
    ax.text(4.27, 20.8, "20%", ha="right", va="bottom")
    ax.text(4.27, 10.8, "10%", ha="right", va="bottom")

    ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax.set_ylabel(r"Median expected sensitivity to $|L/T|$ (%)")
    ax.set_title(r"Projected $\pi^0$ Rosenbluth reach vs. $Q^2$")
    ax.set_xlim(1.1, 4.45)
    ax.set_yscale("log")
    ax.set_ylim(5.0, 500.0)
    ax.set_yticks([5, 10, 20, 50, 100, 200, 500])
    ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())
    ax.get_yaxis().set_minor_formatter(plt.NullFormatter())
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def plot_q2_unlocked(q2unlock, outfile):
    """Show where luminosity creates new 10%-level L/T capability."""
    fig, ax = plt.subplots(figsize=(7.6, 5.5))
    for f, g in q2unlock.groupby("future_luminosity_multiplier", sort=True):
        ax.plot(
            g.Q2_GeV2,
            g.newly_below_10pct_95CL,
            marker="o",
            label=fr"$f={f:g}$",
        )
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax.set_ylabel(r"New bins with expected $|L/T|<10\%$ sensitivity")
    ax.set_title("High-luminosity capability unlocked versus $Q^2$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=180)
    plt.close(fig)

def _pick_col(df, aliases, what, required=True):
    """Find a column by a compact list of accepted aliases."""
    lower = {c.lower(): c for c in df.columns}
    for name in aliases:
        if name.lower() in lower:
            return lower[name.lower()]
    if required:
        raise RuntimeError(
            f"Could not identify {what}. Tried {aliases}. Available columns:\n  "
            + ", ".join(df.columns)
        )
    return None


def _standardize_internal_cross_sections(path, campaign):
    """Read an INTERNAL measured reduced-cross-section table into common names."""
    df = pd.read_csv(path)

    aliases = {
        "Q2": ["Q2_flux_coordinate_GeV2", "Q2", "q2", "Q2_GeV2", "q2_mean", "mean_Q2", "Q2_mean"],
        "xB": ["xB_flux_coordinate", "xB", "xb", "x_B", "xB_mean", "mean_xB"],
        "mt": ["minus_t_center_GeV2", "minus_t", "minus_t_GeV2", "-t", "t_abs", "abs_t", "mt", "t_mean"],
        "phi": ["phi", "phi_deg", "phi_center_deg", "phi_center", "mean_phi", "phi_mean"],
        "eps": ["virtual_photon_epsilon", "epsilon", "eps", "epsilon_mean", "mean_epsilon"],
        "xs": ["reduced_cross_section_nb_per_GeV2_rad", "reduced_cross_section", "cross_section", "sigma", "xsec", "xs", "value"],
        "err": ["propagated_statistical_and_finite_MC_uncertainty_nb_per_GeV2_rad",
                "total_uncertainty", "cross_section_uncertainty", "sigma_unc", "xsec_unc",
                "uncertainty", "error", "err", "total_error"],
        "stat": ["stat_uncertainty", "stat_error", "stat_err", "sigma_stat", "xsec_stat"],
        "sys": ["sys_uncertainty", "syst_uncertainty", "sys_error", "syst_error",
                "sigma_sys", "xsec_sys"],
        "iq2": ["iq2", "iQ2", "q2_bin", "Q2_bin"],
        "ixb": ["ixb", "iXB", "xb_bin", "xB_bin"],
        "it": ["it", "iT", "t_bin", "mt_bin"],
    }

    cols = {k: _pick_col(df, v, k, required=(k in {"Q2", "xB", "mt", "phi", "eps", "xs"}))
            for k, v in aliases.items()}

    print(f"[INTERNAL {campaign} measured-data columns]")
    for key in ("Q2", "xB", "mt", "phi", "eps", "xs", "err"):
        print(f"  {key:>4s} <- {cols[key]}")

    if cols["err"] is not None:
        err = pd.to_numeric(df[cols["err"]], errors="coerce").to_numpy(float)
    elif cols["stat"] is not None:
        stat = pd.to_numeric(df[cols["stat"]], errors="coerce").to_numpy(float)
        if cols["sys"] is not None:
            sys = pd.to_numeric(df[cols["sys"]], errors="coerce").to_numpy(float)
            err = np.sqrt(stat**2 + sys**2)
        else:
            err = stat
    else:
        raise RuntimeError(
            f"Could not identify an absolute uncertainty column in {path}. "
            "For the internal extraction I need measured central values and absolute errors.\n"
            f"Available columns:\n  {', '.join(df.columns)}"
        )

    out = pd.DataFrame({
        "campaign": campaign,
        "Q2_GeV2": pd.to_numeric(df[cols["Q2"]], errors="coerce"),
        "xB": pd.to_numeric(df[cols["xB"]], errors="coerce"),
        "minus_t_GeV2": np.abs(pd.to_numeric(df[cols["mt"]], errors="coerce")),
        "phi_deg": pd.to_numeric(df[cols["phi"]], errors="coerce"),
        "epsilon": pd.to_numeric(df[cols["eps"]], errors="coerce"),
        "sigma": pd.to_numeric(df[cols["xs"]], errors="coerce"),
        "delta_sigma": err,
    })
    for key in ("iq2", "ixb", "it"):
        if cols[key] is not None:
            out[key] = df[cols[key]].to_numpy()

    good = np.all(np.isfinite(out[["Q2_GeV2", "xB", "minus_t_GeV2", "phi_deg",
                                    "epsilon", "sigma", "delta_sigma"]]), axis=1)
    good &= out.delta_sigma.to_numpy(float) > 0
    return out.loc[good].reset_index(drop=True)


def _group_internal_bins(df):
    """Return one group per native (Q2,xB,t) cell."""
    keys = [k for k in ("iq2", "ixb", "it") if k in df.columns]
    if len(keys) == 3:
        return list(df.groupby(keys, sort=True))

    # Fallback for source tables without integer bin IDs: rounded native centers.
    tmp = df.copy()
    tmp["_Q2key"] = tmp.Q2_GeV2.round(4)
    tmp["_xBkey"] = tmp.xB.round(5)
    tmp["_tkey"] = tmp.minus_t_GeV2.round(4)
    return list(tmp.groupby(["_Q2key", "_xBkey", "_tkey"], sort=True))


def _native_cells(df):
    rows = []
    for key, g in _group_internal_bins(df):
        rows.append(dict(
            native_key=str(key),
            Q2_GeV2=float(np.average(g.Q2_GeV2)),
            xB=float(np.average(g.xB)),
            minus_t_GeV2=float(np.average(g.minus_t_GeV2)),
            epsilon=float(np.average(g.epsilon)),
            n_phi=len(g),
        ))
    return pd.DataFrame(rows)


def _nearest_native_group(df, q2, xb, mt):
    """Choose the native cell nearest a requested shared cell in fractional distance."""
    cells = _native_cells(df)
    d2 = ((cells.Q2_GeV2-q2)/max(abs(q2), 0.5))**2
    d2 += ((cells.xB-xb)/max(abs(xb), 0.1))**2
    d2 += ((cells.minus_t_GeV2-mt)/max(abs(mt), 0.15))**2
    row = cells.loc[d2.idxmin()]

    # Reconstruct the selected group using its native coordinates.  This is only a
    # fallback path; integer bin IDs, when present, are preserved by the source table.
    groups = _group_internal_bins(df)
    for key, g in groups:
        qg = float(np.average(g.Q2_GeV2))
        xg = float(np.average(g.xB))
        tg = float(np.average(g.minus_t_GeV2))
        if np.isclose(qg, row.Q2_GeV2) and np.isclose(xg, row.xB) and np.isclose(tg, row.minus_t_GeV2):
            return g.copy()
    raise RuntimeError("Internal error while matching a native measured cell.")


def _joint_rosenbluth_fit(rga, rgk):
    """Fit T,L,LT,TT directly to measured phi-dependent cross sections at two epsilons."""
    frames = []
    for g in (rga, rgk):
        phi = np.deg2rad(g.phi_deg.to_numpy(float))
        eps = g.epsilon.to_numpy(float)
        A = np.column_stack([
            np.ones(len(g)),
            eps,
            np.sqrt(2.0*eps*(1.0+eps))*np.cos(phi),
            eps*np.cos(2.0*phi),
        ])
        frames.append((A, g.sigma.to_numpy(float), g.delta_sigma.to_numpy(float)))

    A = np.vstack([x[0] for x in frames])
    y = np.concatenate([x[1] for x in frames])
    dy = np.concatenate([x[2] for x in frames])
    W = 1.0/dy**2
    normal = A.T @ (W[:, None]*A)
    cov = np.linalg.pinv(normal)
    theta = cov @ (A.T @ (W*y))
    residual = y - A@theta
    chi2 = float(np.sum((residual/dy)**2))
    ndf = int(len(y)-len(theta))
    return theta, cov, chi2, ndf


def build_internal_data_extraction(common, rga_file, rgk_file):
    """INTERNAL ONLY: extract measured T,L,LT,TT and L/T at common Stage-2 cells."""
    rga = _standardize_internal_cross_sections(rga_file, "RGA")
    rgk = _standardize_internal_cross_sections(rgk_file, "RGK")
    rows = []

    for r in common.itertuples(index=False):
        q2 = float(getattr(r, "Q2_common_GeV2"))
        xb = float(getattr(r, "xB_common"))
        mt = float(getattr(r, "minus_t_common_GeV2"))
        # Common-bin identity is the six nominal bin edges.  The Stage-2
        # common table intentionally has no point_id; RGA and RGK are matched
        # by these shared edges.  Use the campaign-specific integer bin
        # indices to select the measured phi rows exactly (no nearest-neighbor
        # matching in flux-coordinate Q2/xB/t).
        iq2_rga = int(getattr(r, "iq2_rga"))
        ixb_rga = int(getattr(r, "ixb_rga"))
        it_rga = int(getattr(r, "it_rga"))
        iq2_rgk = int(getattr(r, "iq2_rgk"))
        ixb_rgk = int(getattr(r, "ixb_rgk"))
        it_rgk = int(getattr(r, "it_rgk"))

        point_id = (
            f"Q2_{float(getattr(r, 'Q2_low_GeV2')):.6g}_"
            f"{float(getattr(r, 'Q2_high_GeV2')):.6g}__"
            f"xB_{float(getattr(r, 'xB_low')):.6g}_"
            f"{float(getattr(r, 'xB_high')):.6g}__"
            f"mt_{float(getattr(r, 'minus_t_low_GeV2')):.6g}_"
            f"{float(getattr(r, 'minus_t_high_GeV2')):.6g}"
        )

        ga = rga[
            (rga["iq2"] == iq2_rga) &
            (rga["ixb"] == ixb_rga) &
            (rga["it"] == it_rga)
        ].copy()
        gk = rgk[
            (rgk["iq2"] == iq2_rgk) &
            (rgk["ixb"] == ixb_rgk) &
            (rgk["it"] == it_rgk)
        ].copy()

        if ga.empty or gk.empty:
            raise RuntimeError(
                "Missing measured rows for common cell "
                f"(Q2={q2:.6g}, xB={xb:.6g}, -t={mt:.6g}); "
                f"RGA indices=({iq2_rga},{ixb_rga},{it_rga}) n={len(ga)}, "
                f"RGK indices=({iq2_rgk},{ixb_rgk},{it_rgk}) n={len(gk)}."
            )
        theta, cov, chi2, ndf = _joint_rosenbluth_fit(ga, gk)
        T, L, LT, TT = map(float, theta)
        dT, dL, dLT, dTT = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))

        if T != 0:
            R = L/T
            grad = np.array([-L/T**2, 1.0/T, 0.0, 0.0])
            varR = float(grad @ cov @ grad)
            dR = math.sqrt(max(varR, 0.0))
        else:
            R, dR = np.nan, np.nan

        rows.append(dict(
            point_id=point_id,
            Q2_GeV2=q2,
            xB=xb,
            minus_t_GeV2=mt,
            epsilon_rga=float(np.average(ga.epsilon)),
            epsilon_rgk=float(np.average(gk.epsilon)),
            delta_epsilon=float(np.average(ga.epsilon)-np.average(gk.epsilon)),
            n_phi_rga=len(ga),
            n_phi_rgk=len(gk),
            sigma_T=T,
            delta_sigma_T=dT,
            sigma_L=L,
            delta_sigma_L=dL,
            sigma_LT=LT,
            delta_sigma_LT=dLT,
            sigma_TT=TT,
            delta_sigma_TT=dTT,
            R_L_over_T=R,
            delta_R_L_over_T=dR,
            chi2=chi2,
            ndf=ndf,
            chi2_ndf=chi2/ndf if ndf > 0 else np.nan,
        ))
    return pd.DataFrame(rows)


def _q2_offsets(values, width=0.11):
    """Small deterministic horizontal offsets so multiple cells at one Q2 remain visible."""
    values = np.asarray(values, float)
    out = np.zeros(len(values), float)
    for q2 in np.unique(values):
        idx = np.flatnonzero(np.isclose(values, q2))
        if len(idx) > 1:
            out[idx] = np.linspace(-width, width, len(idx))
    return out


def plot_internal_LT_vs_Q2(data, outfile):
    """Harut study: every measured matched-cell L/T extraction versus Q2."""
    g = data.sort_values(["Q2_GeV2", "xB", "minus_t_GeV2"]).reset_index(drop=True)
    x = g.Q2_GeV2.to_numpy(float) + _q2_offsets(g.Q2_GeV2.to_numpy(float))

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    sc = ax.scatter(
        x, g.R_L_over_T, c=g.xB, s=42 + 34*g.minus_t_GeV2,
        cmap="viridis", edgecolor="black", linewidth=0.35, zorder=3,
    )
    ax.errorbar(
        x, g.R_L_over_T, yerr=g.delta_R_L_over_T,
        fmt="none", ecolor="0.45", elinewidth=0.9, capsize=1.8, alpha=0.75, zorder=2,
    )
    ax.axhline(0.0, linewidth=1.0, color="black")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax.set_ylabel(r"Measured $\sigma_L/\sigma_T$")
    ax.set_title(r"INTERNAL: measured $\pi^0$ Rosenbluth $L/T$ vs. $Q^2$")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$x_B$")
    ax.text(
        0.99, 0.02, r"Marker size increases with $-t$",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def plot_internal_LT_by_xB(data, outfile):
    """Faceted internal view so Q2 evolution is not confused with xB/t evolution."""
    xb_values = sorted(data.xB.unique())
    ncols = min(4, max(1, len(xb_values)))
    nrows = int(np.ceil(len(xb_values)/ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0*ncols, 3.5*nrows),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()

    finite_t = data.minus_t_GeV2[np.isfinite(data.minus_t_GeV2)]
    tmin = float(finite_t.min()) if len(finite_t) else 0.0
    tmax = float(finite_t.max()) if len(finite_t) else 1.0

    for ax, xb in zip(axes, xb_values):
        g = data[np.isclose(data.xB, xb)].sort_values(["Q2_GeV2", "minus_t_GeV2"])
        sc = ax.scatter(g.Q2_GeV2, g.R_L_over_T, c=g.minus_t_GeV2,
                        cmap="viridis", vmin=tmin, vmax=tmax,
                        s=48, edgecolor="black", linewidth=0.35, zorder=3)
        ax.errorbar(g.Q2_GeV2, g.R_L_over_T, yerr=g.delta_R_L_over_T,
                    fmt="none", ecolor="0.45", elinewidth=0.9, capsize=1.8, alpha=0.75)
        ax.axhline(0.0, linewidth=0.8, color="black")
        ax.set_title(fr"$x_B={xb:g}$")
        ax.grid(alpha=0.15)

    for ax in axes[len(xb_values):]:
        ax.set_visible(False)
    for i, ax in enumerate(axes[:len(xb_values)]):
        if i//ncols == nrows-1 or i+ncols >= len(xb_values):
            ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
        if i % ncols == 0:
            ax.set_ylabel(r"$\sigma_L/\sigma_T$")

    cbar = fig.colorbar(sc, ax=list(axes[:len(xb_values)]), shrink=0.88, pad=0.02)
    cbar.set_label(r"$-t$ (GeV$^2$)")
    fig.suptitle(r"INTERNAL: measured $\pi^0$ $L/T$ evolution at fixed $x_B$", y=0.995)
    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.09, top=0.91, wspace=0.12, hspace=0.28)
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def run_internal_data_mode(a):
    """Explicitly non-default path using measured, unapproved central values."""
    out = a.internal_output.resolve()
    tabs, figs = out/"tables", out/"figures"
    tabs.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    common = pd.read_csv(a.stage2.resolve()/"tables"/"03_common_rosenbluth_model_points.csv")
    data = build_internal_data_extraction(common, a.internal_rga.resolve(), a.internal_rgk.resolve())
    data.to_csv(tabs/"INTERNAL_01_measured_LT_by_point.csv", index=False)
    plot_internal_LT_vs_Q2(data, figs/"INTERNAL_01_measured_L_over_T_vs_Q2.png")
    plot_internal_LT_by_xB(data, figs/"INTERNAL_02_measured_L_over_T_vs_Q2_by_xB.png")

    print("\n*** INTERNAL DATA MODE: measured, unapproved RGA/RGK central values ***")
    print(f"Matched Rosenbluth cells: {len(data)}")
    print("No clipping or positivity constraint is applied to sigma_L or L/T.")
    print("Negative/noisy values are retained intentionally.")
    print("\nWrote:")
    print(f"  {tabs/'INTERNAL_01_measured_LT_by_point.csv'}")
    print(f"  {figs/'INTERNAL_01_measured_L_over_T_vs_Q2.png'}")
    print(f"  {figs/'INTERNAL_02_measured_L_over_T_vs_Q2_by_xB.png'}")


def main():
    a = args()

    if a.central_values == "data":
        run_internal_data_mode(a)
        return

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

    q2sum = q2_evolution_summary(points)
    q2unlock = q2_unlock_summary(points)
    benchmarks = benchmark_regions(points)

    points.to_csv(tabs/"01_measurement_reach_by_point.csv", index=False)
    summary.to_csv(tabs/"02_measurement_reach_summary.csv", index=False)
    q2sum.to_csv(tabs/"04_Q2_evolution_summary.csv", index=False)
    q2unlock.to_csv(tabs/"05_Q2_new_capability_vs_recorded.csv", index=False)
    benchmarks.to_csv(tabs/"06_HallA_benchmark_nearest_cells.csv", index=False)

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
    plot_q2_evolution(q2sum, figs/"04_Q2_evolution_L_over_T_reach.png")
    plot_q2_unlocked(q2unlock, figs/"05_Q2_new_10pct_capability.png")

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

    print("\n[Q2-dependent 10% longitudinal-fraction reach]")
    q2print = q2sum[[
        "future_luminosity_multiplier", "Q2_GeV2", "n_bins",
        "n_expected_95pct_limit_lt_0p10",
        "median_expected_95pct_abs_L_over_T_limit_if_L_zero",
    ]]
    print(q2print.to_string(index=False))

    print("\nWrote:")
    print(f"  {tabs/'01_measurement_reach_by_point.csv'}")
    print(f"  {tabs/'02_measurement_reach_summary.csv'}")
    print(f"  {tabs/'03_new_capability_vs_recorded.csv'}")
    print(f"  {tabs/'04_Q2_evolution_summary.csv'}")
    print(f"  {tabs/'05_Q2_new_capability_vs_recorded.csv'}")
    print(f"  {tabs/'06_HallA_benchmark_nearest_cells.csv'}")
    print(f"  {figs/'01_L_over_T_constraint_counts.png'}")
    print("  per-scenario L/T reach maps")
    print(f"  {figs/'03_absolute_structure_function_uncertainties.png'}")
    print(f"  {figs/'04_Q2_evolution_L_over_T_reach.png'}")
    print(f"  {figs/'05_Q2_new_10pct_capability.png'}")
    print("\nCaveat: this inherits the Stage-3 assumption that supplied fractional")
    print("uncertainties scale as 1/sqrt(exposure); irreducible systematic and finite-MC")
    print("floors are not yet separated.")


if __name__ == "__main__":
    main()
