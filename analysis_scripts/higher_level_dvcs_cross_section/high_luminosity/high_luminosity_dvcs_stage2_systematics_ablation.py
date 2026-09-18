#!/usr/bin/env python3
"""
Stage 2C: systematic-ablation study for the pass-2 matched XS+BSA CFF projection.

This script does NOT rerun Gepard.  It reuses the exact finite-change CFF
derivatives already produced by:
    high_luminosity_dvcs_stage2_cff_sensitivity_pass2.py

It asks where luminosity stops helping because fixed systematics dominate.

Default scenarios
-----------------
baseline
    Current projection: pass-2 point-to-point systematics, provisional 10%
    fully correlated XS normalization, 4% BSA beam-polarization scale.

statistics_only
    Statistical errors only.  This is an ideal reference, not a realistic
    experimental scenario.

improved_norm
    Point-to-point systematics unchanged; XS normalization 10% -> 5%.

improved_ptp
    XS and BSA point-to-point systematics reduced by a factor 2;
    XS normalization remains 10%.

improved_both
    XS normalization 10% -> 5% AND both XS/BSA point-to-point systematics
    reduced by a factor 2.

The improvement factors are CLI options, so workshop alternatives can be
generated without editing the code.

Outputs include:
  * all-cell and high-|t| (default >=0.4 GeV^2) median CFF precision;
  * maximum |t| satisfying 25% and 50% ReH precision;
  * the highest-|t| cell precision;
  * figures showing the luminosity/systematics tradeoff.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)


def load_fit_module(path: Path):
    spec = importlib.util.spec_from_file_location("pass2_cff_fit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def scenario_inputs(
    base: pd.DataFrame,
    name: str,
    improved_norm: float,
    ptp_factor: float,
) -> pd.DataFrame:
    d = base.copy()

    if name == "baseline":
        return d

    if name == "statistics_only":
        d["xs_ptp_sys_pseudo_abs"] = 0.0
        d["bsa_ptp_sys_pseudo_abs"] = 0.0
        d["xs_scale_frac"] = 0.0
        d["bsa_scale_frac"] = 0.0
        return d

    if name == "improved_norm":
        d["xs_scale_frac"] = float(improved_norm)
        return d

    if name == "improved_ptp":
        d["xs_ptp_sys_pseudo_abs"] *= float(ptp_factor)
        d["bsa_ptp_sys_pseudo_abs"] *= float(ptp_factor)
        return d

    if name == "improved_both":
        d["xs_scale_frac"] = float(improved_norm)
        d["xs_ptp_sys_pseudo_abs"] *= float(ptp_factor)
        d["bsa_ptp_sys_pseudo_abs"] *= float(ptp_factor)
        return d

    raise ValueError(name)


def max_t_under(d: pd.DataFrame, column: str, threshold: float) -> float:
    q = d[np.isfinite(d[column]) & (d[column] < threshold)]
    return float(q["t_abs"].max()) if len(q) else np.nan


def metrics_for(d: pd.DataFrame, high_t_min: float) -> Dict[str, float]:
    hi = d[d["t_abs"] >= high_t_min]
    top = d.loc[d["t_abs"].idxmax()]
    return {
        "n_cells": len(d),
        "n_high_t_cells": len(hi),
        "median_rel_ReH_pct": 100.0 * np.nanmedian(d["relative_sigma_ReH"]),
        "median_rel_ImH_pct": 100.0 * np.nanmedian(d["relative_sigma_ImH"]),
        "median_abs_corr": np.nanmedian(np.abs(d["corr_ReH_ImH"])),
        "high_t_median_rel_ReH_pct": 100.0 * np.nanmedian(hi["relative_sigma_ReH"]),
        "high_t_median_rel_ImH_pct": 100.0 * np.nanmedian(hi["relative_sigma_ImH"]),
        "high_t_median_abs_corr": np.nanmedian(np.abs(hi["corr_ReH_ImH"])),
        "max_t_ReH_lt25pct": max_t_under(d, "relative_sigma_ReH", 0.25),
        "max_t_ReH_lt50pct": max_t_under(d, "relative_sigma_ReH", 0.50),
        "highest_t": float(top["t_abs"]),
        "highest_t_xB": float(top["xB"]),
        "highest_t_Q2": float(top["Q2"]),
        "highest_t_nphi": int(top["n_phi"]),
        "highest_t_rel_ReH_pct": 100.0 * float(top["relative_sigma_ReH"]),
        "highest_t_rel_ImH_pct": 100.0 * float(top["relative_sigma_ImH"]),
    }


def savefig(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def make_plots(summary: pd.DataFrame, figures: Path, high_t_min: float):
    order = [
        ("statistics_only", "Statistics only"),
        ("baseline", "Baseline systematics"),
        ("improved_norm", "Normalization improved"),
        ("improved_ptp", "Point-to-point improved"),
        ("improved_both", "Both improved"),
    ]

    # High-t ReH precision.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for scenario, label in order:
        d = summary[summary["scenario"] == scenario].sort_values("luminosity_factor")
        ax.plot(d["luminosity_factor"], d["high_t_median_rel_ReH_pct"],
                marker="o", label=label)
    ax.set_xscale("log")
    ax.set_xticks(LUMI_FACTORS, [f"{x}x" for x in LUMI_FACTORS])
    ax.set_xlabel("Luminosity relative to pass-2 exposure")
    ax.set_ylabel(r"Median relative uncertainty on $\mathrm{Re}\,\mathcal{H}$ (%)")
    ax.set_title(rf"High-$|t|$ cells: $|t|\geq {high_t_min:.1f}$ GeV$^2$")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "01_high_t_ReH_systematics_ablation.png")

    # High-t ImH precision.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for scenario, label in order:
        d = summary[summary["scenario"] == scenario].sort_values("luminosity_factor")
        ax.plot(d["luminosity_factor"], d["high_t_median_rel_ImH_pct"],
                marker="o", label=label)
    ax.set_xscale("log")
    ax.set_xticks(LUMI_FACTORS, [f"{x}x" for x in LUMI_FACTORS])
    ax.set_xlabel("Luminosity relative to pass-2 exposure")
    ax.set_ylabel(r"Median relative uncertainty on $\mathrm{Im}\,\mathcal{H}$ (%)")
    ax.set_title(rf"High-$|t|$ cells: $|t|\geq {high_t_min:.1f}$ GeV$^2$")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "02_high_t_ImH_systematics_ablation.png")

    # ReH 25% precision frontier.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for scenario, label in order:
        d = summary[summary["scenario"] == scenario].sort_values("luminosity_factor")
        ax.plot(d["luminosity_factor"], d["max_t_ReH_lt25pct"],
                marker="o", label=label)
    ax.set_xscale("log")
    ax.set_xticks(LUMI_FACTORS, [f"{x}x" for x in LUMI_FACTORS])
    ax.set_xlabel("Luminosity relative to pass-2 exposure")
    ax.set_ylabel(r"Maximum $|t|$ with $\delta\mathrm{Re}\mathcal{H}/|\mathrm{Re}\mathcal{H}|<25\%$ (GeV$^2$)")
    ax.set_title(r"Precision $|t|$ frontier under $|t|/Q^2<0.2$")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "03_ReH_25pct_t_frontier_systematics_ablation.png")

    # Highest-t cell.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for scenario, label in order:
        d = summary[summary["scenario"] == scenario].sort_values("luminosity_factor")
        ax.plot(d["luminosity_factor"], d["highest_t_rel_ReH_pct"],
                marker="o", label=label)
    ax.set_xscale("log")
    ax.set_xticks(LUMI_FACTORS, [f"{x}x" for x in LUMI_FACTORS])
    ax.set_xlabel("Luminosity relative to pass-2 exposure")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Re}\,\mathcal{H}$ (%)")
    ax.set_title(r"Highest-$|t|$ clean cell")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "04_highest_t_cell_ReH_systematics_ablation.png")


def main(argv: List[str] | None = None) -> int:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input-dir",
        default=str(here / "output_stage2_pass2" / "tables"),
        help="Stage-2A luminosity input tables.",
    )
    p.add_argument(
        "--cff-output-dir",
        default=str(here / "output_stage2_pass2_cff"),
        help="Existing joint CFF output containing derivative/cell tables.",
    )
    p.add_argument(
        "--cff-script",
        default=str(here / "high_luminosity_dvcs_stage2_cff_sensitivity_pass2.py"),
        help="Joint CFF script whose covariance machinery is reused.",
    )
    p.add_argument(
        "--outdir",
        default=str(here / "output_stage2_pass2_systematics_ablation"),
    )
    p.add_argument(
        "--improved-xs-normalization",
        type=float,
        default=0.05,
        help="Improved fully correlated XS normalization fraction (default 0.05).",
    )
    p.add_argument(
        "--ptp-factor",
        type=float,
        default=0.50,
        help="Multiplicative factor applied to XS and BSA point-to-point systematics.",
    )
    p.add_argument(
        "--high-t-min",
        type=float,
        default=0.40,
        help="High-|t| threshold in GeV^2.",
    )
    args = p.parse_args(argv)

    input_dir = Path(args.input_dir).resolve()
    cff_dir = Path(args.cff_output_dir).resolve()
    outdir = Path(args.outdir).resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    fitmod = load_fit_module(Path(args.cff_script).resolve())
    deriv = pd.read_csv(cff_dir / "tables" / "point_xs_bsa_cff_derivatives.csv")
    cell_meta = pd.read_csv(cff_dir / "tables" / "cell_km15_cffs.csv")
    cell_meta["bin"] = cell_meta["bin"].astype(int)

    scenarios = (
        "statistics_only",
        "baseline",
        "improved_norm",
        "improved_ptp",
        "improved_both",
    )

    all_cells = []
    summary_rows = []
    diag_rows = []

    print("=" * 94)
    print("PASS-2 XS+BSA SYSTEMATIC-ABLATION STUDY")
    print("=" * 94)
    print("baseline XS normalization       : 10.0%")
    print(f"improved XS normalization       : {100*args.improved_xs_normalization:.1f}%")
    print(f"improved point-to-point factor  : {args.ptp_factor:.3f}")
    print("BSA beam-polarization scale     : 4.0% (held fixed except statistics-only)")
    print(f"high-|t| definition             : |t| >= {args.high_t_min:.2f} GeV^2")

    for factor in LUMI_FACTORS:
        base = pd.read_csv(input_dir / f"joint_fit_input_km15_{factor}x.csv")
        base["bin"] = base["bin"].astype(int)

        for scenario in scenarios:
            d = scenario_inputs(
                base,
                scenario,
                args.improved_xs_normalization,
                args.ptp_factor,
            )
            cov, layout, diag = fitmod.covariance(
                d, deriv, cell_meta, "ReH_ImH_XS_BSA"
            )
            cells = fitmod.summarize(
                cell_meta, cov, layout, "ReH_ImH_XS_BSA", factor
            )
            cells["scenario"] = scenario
            all_cells.append(cells)

            m = metrics_for(cells, args.high_t_min)
            m["scenario"] = scenario
            m["luminosity_factor"] = factor
            summary_rows.append(m)

            diag["scenario"] = scenario
            diag["luminosity_factor"] = factor
            diag_rows.append(diag)

    cell_results = pd.concat(all_cells, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    diagnostics = pd.DataFrame(diag_rows)

    cell_results.to_csv(tables / "systematics_ablation_cell_uncertainties.csv", index=False)
    summary.to_csv(tables / "systematics_ablation_summary.csv", index=False)
    diagnostics.to_csv(tables / "systematics_ablation_global_diagnostics.csv", index=False)

    make_plots(summary, figures, args.high_t_min)

    cols = [
        "scenario", "luminosity_factor",
        "high_t_median_rel_ReH_pct", "high_t_median_rel_ImH_pct",
        "max_t_ReH_lt25pct", "highest_t_rel_ReH_pct",
        "highest_t_rel_ImH_pct",
    ]
    print("\nHigh-|t| and reach summary:")
    print(summary[cols].sort_values(["scenario", "luminosity_factor"]).to_string(
        index=False, float_format=lambda x: f"{x:.3g}"
    ))

    # Compact "what limits 10x?" comparison.
    ten = summary[summary["luminosity_factor"] == 10].copy()
    print("\n10x diagnostic:")
    print(ten[[
        "scenario", "high_t_median_rel_ReH_pct",
        "high_t_median_rel_ImH_pct", "highest_t_rel_ReH_pct",
        "highest_t_rel_ImH_pct"
    ]].to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    print(f"\n[output] {outdir}")
    print("[interpretation] compare baseline with improved_norm and improved_ptp")
    print("                 to identify which fixed systematic limits luminosity gains.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
