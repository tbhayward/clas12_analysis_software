#!/usr/bin/env python3
"""
Stage 3A: subtraction-constant sensitivity from the matched pass-2 XS+BSA study.

Purpose
-------
Move beyond independent ReH values in every (xB,Q2,t) cell and ask how well
the data constrain an xi-independent additive contribution to ReH at fixed t,
i.e. the quantity that plays the role of the dispersion-relation subtraction
constant.

This is intentionally a *sensitivity* study around KM15, not yet a claim of a
model-independent D-term extraction.  KM15 supplies the central CFFs.  We fit
a common additive shift dC(t) to ReH for all cells in a narrow |t| band:

    ReH_i -> ReH_i + dC_k,     i in t-band k.

The response dO/dC is therefore exactly the already-computed dO/dReH for each
point.  ImH remains free independently in each (xB,Q2,t) cell, so the BSA can
constrain the absorptive part rather than forcing it to KM15.

The absolute KM15 subtraction constant and the convention-dependent mapping
from C(t) to D(t) are deliberately NOT invented here.  This stage establishes
the experimental uncertainty on an additive subtraction constant versus t.
A later Stage 3B can attach the explicit dispersion-relation convention and
D(t) normalization after validating the corresponding Gepard/PARTONS model
implementation.

Systematic scenarios
--------------------
baseline:
    current pass-2 point-to-point systematics, provisional 10% fully
    correlated XS normalization, 4% BSA polarization scale.

ptp_half:
    XS and BSA point-to-point systematics reduced by factor 2.

statistics_only:
    statistical errors only.

Luminosity factors: 1, 2, 5, 10 relative to pass-2 exposure.

Inputs
------
output_stage2_pass2/tables/joint_fit_input_km15_{1,2,5,10}x.csv
output_stage2_pass2_cff/tables/point_xs_bsa_cff_derivatives.csv
output_stage2_pass2_cff/tables/cell_km15_cffs.csv

Outputs
-------
output_stage3_subtraction_constant/
    tables/subtraction_constant_uncertainties.csv
    tables/global_fit_diagnostics.csv
    tables/t_band_definition.csv
    figures/*.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)
SCENARIOS = ("baseline", "ptp_half", "statistics_only")


def assign_t_bands(cell_meta: pd.DataFrame, width: float) -> pd.DataFrame:
    """Assign cells to fixed-width |t| bands; preserve actual mean t per band."""
    c = cell_meta.copy()
    # Use [0,width), [width,2width), ... but center labels on occupied data.
    c["t_band_index"] = np.floor(c["t_abs"] / width + 1e-12).astype(int)
    return c


def scenario_input(d: pd.DataFrame, scenario: str) -> pd.DataFrame:
    q = d.copy()
    if scenario == "baseline":
        return q
    if scenario == "ptp_half":
        q["xs_ptp_sys_pseudo_abs"] *= 0.5
        q["bsa_ptp_sys_pseudo_abs"] *= 0.5
        return q
    if scenario == "statistics_only":
        q["xs_ptp_sys_pseudo_abs"] = 0.0
        q["bsa_ptp_sys_pseudo_abs"] = 0.0
        q["xs_scale_frac"] = 0.0
        q["bsa_scale_frac"] = 0.0
        return q
    raise ValueError(scenario)


def build_covariance(
    data: pd.DataFrame,
    deriv: pd.DataFrame,
    cell_meta: pd.DataFrame,
) -> Tuple[np.ndarray, List[Tuple[str, int]], Dict[str, float]]:
    """
    Parameters:
      dC_k                 one common additive ReH shift per t band
      dImH_bin             one independent ImH shift per original cell
      beta_xs_norm         common XS normalization nuisance
      beta_bsa_pol         common BSA polarization nuisance
    """
    use = data.merge(
        deriv,
        on=["point_id", "bin", "phi_deg"],
        how="inner",
        validate="one_to_one",
    )
    band_map = dict(
        zip(cell_meta["bin"].astype(int), cell_meta["t_band_index"].astype(int))
    )
    cells = sorted(cell_meta["bin"].astype(int).unique())
    bands = sorted(cell_meta["t_band_index"].astype(int).unique())

    layout: List[Tuple[str, int]] = []
    layout += [("dC", b) for b in bands]
    layout += [("dImH", b) for b in cells]
    layout += [("beta_xs_norm", -1), ("beta_bsa_pol", -1)]
    idx = {key: i for i, key in enumerate(layout)}
    npar = len(layout)

    rows, sigmas = [], []

    for r in use.itertuples(index=False):
        b = int(r.bin)
        tb = band_map[b]

        # XS row
        v = np.zeros(npar)
        v[idx[("dC", tb)]] = float(r.d_xs_d_ReH)
        v[idx[("dImH", b)]] = float(r.d_xs_d_ImH)
        v[idx[("beta_xs_norm", -1)]] = float(r.xs_scale_frac) * float(r.xs_km15)
        sx = math.hypot(float(r.xs_stat_pseudo_abs),
                        float(r.xs_ptp_sys_pseudo_abs))
        if np.isfinite(sx) and sx > 0:
            rows.append(v)
            sigmas.append(sx)

        # BSA row
        v = np.zeros(npar)
        v[idx[("dC", tb)]] = float(r.d_bsa_d_ReH)
        v[idx[("dImH", b)]] = float(r.d_bsa_d_ImH)
        v[idx[("beta_bsa_pol", -1)]] = (
            float(r.bsa_scale_frac) * float(r.bsa_km15)
        )
        sa = math.hypot(float(r.bsa_stat_pseudo_abs),
                        float(r.bsa_ptp_sys_pseudo_abs))
        if np.isfinite(sa) and sa > 0:
            rows.append(v)
            sigmas.append(sa)

    A = np.asarray(rows, dtype=float)
    sig = np.asarray(sigmas, dtype=float)
    B = A / sig[:, None]

    # Unit Gaussian priors for the two scale nuisances.
    for key in (("beta_xs_norm", -1), ("beta_bsa_pol", -1)):
        p = np.zeros(npar)
        p[idx[key]] = 1.0
        B = np.vstack((B, p))

    info = B.T @ B
    rank = int(np.linalg.matrix_rank(B))
    sv = np.linalg.svd(B, compute_uv=False)
    threshold = max(float(sv[0]), 1.0) * 1e-12
    positive = sv[sv > threshold]
    cond = float(positive[0] / positive[-1]) if len(positive) > 1 else np.inf
    cov = np.linalg.pinv(info, rcond=1e-12)

    diag = {
        "n_measurement_rows": len(A),
        "n_parameters": npar,
        "rank": rank,
        "rank_deficit": npar - rank,
        "condition_number_effective": cond,
        "xs_norm_beta_sigma": math.sqrt(
            max(cov[idx[("beta_xs_norm", -1)], idx[("beta_xs_norm", -1)]], 0)
        ),
        "bsa_pol_beta_sigma": math.sqrt(
            max(cov[idx[("beta_bsa_pol", -1)], idx[("beta_bsa_pol", -1)]], 0)
        ),
    }
    return cov, layout, diag


def band_results(
    cov: np.ndarray,
    layout: List[Tuple[str, int]],
    cell_meta: pd.DataFrame,
    scenario: str,
    lumi: int,
) -> pd.DataFrame:
    idx = {key: i for i, key in enumerate(layout)}
    rows = []
    for tb, group in cell_meta.groupby("t_band_index", sort=True):
        j = idx[("dC", int(tb))]
        sigma = math.sqrt(max(float(cov[j, j]), 0.0))
        rows.append({
            "scenario": scenario,
            "luminosity_factor": lumi,
            "t_band_index": int(tb),
            "t_min_data": float(group["t_abs"].min()),
            "t_max_data": float(group["t_abs"].max()),
            "t_mean": float(group["t_abs"].mean()),
            "n_cells": int(len(group)),
            "n_phi_points": int(group["n_phi"].sum()),
            "sigma_deltaC": sigma,
        })
    return pd.DataFrame(rows)


def savefig(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def make_plots(results: pd.DataFrame, figures: Path):
    # Baseline luminosity evolution.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for L in LUMI_FACTORS:
        d = results[
            (results["scenario"] == "baseline") &
            (results["luminosity_factor"] == L)
        ].sort_values("t_mean")
        ax.plot(d["t_mean"], d["sigma_deltaC"], marker="o", label=f"{L}x")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"Projected uncertainty on additive subtraction constant $\delta C(t)$")
    ax.set_title(r"Matched XS+BSA subtraction-constant sensitivity")
    ax.grid(alpha=0.2)
    ax.legend(title="Luminosity")
    savefig(fig, figures / "01_subtraction_constant_uncertainty_vs_t.png")

    # Systematics comparison at 10x.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for scenario, label in [
        ("baseline", "Current point-to-point systematics"),
        ("ptp_half", "Point-to-point systematics / 2"),
        ("statistics_only", "Statistics only"),
    ]:
        d = results[
            (results["scenario"] == scenario) &
            (results["luminosity_factor"] == 10)
        ].sort_values("t_mean")
        ax.plot(d["t_mean"], d["sigma_deltaC"], marker="o", label=label)
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"Projected uncertainty on additive subtraction constant $\delta C(t)$")
    ax.set_title(r"Systematic limitation on subtraction-constant sensitivity at 10x")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "02_subtraction_constant_systematics_10x.png")

    # Improvement factor from 1x to higher luminosity, avoids arbitrary
    # pass/fail precision thresholds.
    base = results[
        (results["scenario"] == "baseline") &
        (results["luminosity_factor"] == 1)
    ][["t_band_index", "sigma_deltaC"]].rename(
        columns={"sigma_deltaC": "sigma_1x"}
    )
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for L in (2, 5, 10):
        d = results[
            (results["scenario"] == "baseline") &
            (results["luminosity_factor"] == L)
        ].merge(base, on="t_band_index")
        d["improvement"] = d["sigma_1x"] / d["sigma_deltaC"]
        ax.plot(d["t_mean"], d["improvement"], marker="o", label=f"1x / {L}x")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel("Improvement factor in subtraction-constant precision")
    ax.set_title("Luminosity gain without an arbitrary precision threshold")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "03_subtraction_constant_improvement_factor.png")


def main(argv: List[str] | None = None) -> int:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input-dir",
        default=str(here / "output_stage2_pass2" / "tables"),
    )
    p.add_argument(
        "--cff-output-dir",
        default=str(here / "output_stage2_pass2_cff"),
    )
    p.add_argument(
        "--outdir",
        default=str(here / "output_stage3_subtraction_constant"),
    )
    p.add_argument(
        "--t-band-width",
        type=float,
        default=0.10,
        help="Width of |t| bands in GeV^2 (default 0.10).",
    )
    args = p.parse_args(argv)

    input_dir = Path(args.input_dir).resolve()
    cff_dir = Path(args.cff_output_dir).resolve()
    outdir = Path(args.outdir).resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    deriv = pd.read_csv(cff_dir / "tables" / "point_xs_bsa_cff_derivatives.csv")
    cell_meta = pd.read_csv(cff_dir / "tables" / "cell_km15_cffs.csv")
    cell_meta["bin"] = cell_meta["bin"].astype(int)
    cell_meta = assign_t_bands(cell_meta, args.t_band_width)

    band_def = (
        cell_meta.groupby("t_band_index", as_index=False)
        .agg(
            t_min_data=("t_abs", "min"),
            t_max_data=("t_abs", "max"),
            t_mean=("t_abs", "mean"),
            n_cells=("bin", "size"),
            n_phi_points=("n_phi", "sum"),
        )
    )
    band_def.to_csv(tables / "t_band_definition.csv", index=False)

    print("=" * 94)
    print("STAGE 3A: SUBTRACTION-CONSTANT SENSITIVITY")
    print("=" * 94)
    print(f"usable CFF cells                : {len(cell_meta)}")
    print(f"|t| band width                  : {args.t_band_width:.3f} GeV^2")
    print(f"occupied |t| bands              : {len(band_def)}")
    print("fit per band                    : one common additive dC(t) in ReH")
    print("ImH treatment                   : independent nuisance per original cell")
    print("absolute D(t) normalization     : NOT assigned in this stage")
    print("\nBand definition:")
    print(band_def.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    all_results = []
    diagnostics = []
    for L in LUMI_FACTORS:
        data = pd.read_csv(input_dir / f"joint_fit_input_km15_{L}x.csv")
        data["bin"] = data["bin"].astype(int)

        for scenario in SCENARIOS:
            q = scenario_input(data, scenario)
            cov, layout, diag = build_covariance(q, deriv, cell_meta)
            r = band_results(cov, layout, cell_meta, scenario, L)
            all_results.append(r)
            diag["scenario"] = scenario
            diag["luminosity_factor"] = L
            diagnostics.append(diag)

    results = pd.concat(all_results, ignore_index=True)
    diagnostics = pd.DataFrame(diagnostics)
    results.to_csv(tables / "subtraction_constant_uncertainties.csv", index=False)
    diagnostics.to_csv(tables / "global_fit_diagnostics.csv", index=False)
    make_plots(results, figures)

    print("\nBaseline subtraction-constant uncertainty:")
    baseline = results[results["scenario"] == "baseline"].pivot(
        index="t_mean", columns="luminosity_factor", values="sigma_deltaC"
    )
    baseline.columns = [f"{int(c)}x" for c in baseline.columns]
    print(baseline.to_string(float_format=lambda x: f"{x:.4g}"))

    print("\n10x systematic comparison:")
    ten = results[results["luminosity_factor"] == 10].pivot(
        index="t_mean", columns="scenario", values="sigma_deltaC"
    )
    print(ten.to_string(float_format=lambda x: f"{x:.4g}"))

    print(f"\n[output] {outdir}")
    print("[important] This is the uncertainty on an additive xi-independent ReH")
    print("            subtraction constant around KM15.  It is not yet labeled D(t).")
    print("[next] Validate the explicit dispersion-relation/D-term convention in")
    print("       Gepard/PARTONS, attach the KM15 central C(t), and propagate to D(t).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
