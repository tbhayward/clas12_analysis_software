#!/usr/bin/env python3
"""
Stage 2B: local CFF-H sensitivity projection for high-luminosity CLAS12 DVCS.

This is intentionally a projection, not a claim of a model-independent CFF
extraction.

Input:
  output_stage2/tables/cff_fit_input_km15_{1,2,5,10}x.csv

Pseudo-truth:
  KM15 at the exact published CLAS12 kinematics, already restricted to
  |t|/Q^2 < 0.2 by Stage 2A.

What is varied:
  A. "ReH_only": Re H is allowed to change independently in every
     (xB,Q2,t) cell; Im H and all other CFFs remain at KM15.
  B. "ReH_ImH": Re H and Im H are both allowed to change independently in
     every cell; all other CFFs remain at KM15.

A single correlated normalization nuisance is shared by ALL points.  Its prior
width is read from correlated_norm_frac in the Stage-2A input (31% for the
current pass-1 baseline).

How sensitivities are calculated:
  For each cell, Gepard/KM15 is evaluated after small finite changes to ReH
  and/or ImH.  Thus the script directly measures how much the predicted
  cross section moves when a CFF changes.  No analytic derivative formula is
  assumed.

Because the pseudo-data central values equal KM15, the best-fit CFF shifts are
zero by construction.  The object of this stage is the projected uncertainty
and parameter degeneracy as luminosity changes.

Important interpretation:
  - ReH_only is optimistic: ImH and E, Htilde, Etilde are assumed known.
  - ReH_ImH is a more conservative H-only stress test, but E/Htilde/Etilde
    are still fixed to KM15.
  - Cross sections alone cannot be expected to determine all eight leading
    twist CFF components locally.  BSA/polarized data can be added later.
  - No D-term/pressure transform is performed here.  First establish whether
    ReH is actually constrained over a useful t lever arm.

Typical use:
  python high_luminosity_dvcs_stage2_cff_sensitivity.py

The script uses the same validated Gepard point construction as
extract_emff_from_dvcs_bh.py, so the existing Python environment is expected.
"""

from __future__ import annotations

import argparse
import contextlib
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)
FIT_MODES = ("ReH_only", "ReH_ImH")
DEFAULT_MIN_NPHI = 3
DEFAULT_REL_STEP = 0.02
DEFAULT_ABS_STEP = 0.05


# =============================================================================
# Gepard interface
# =============================================================================

def make_point(g, row):
    """
    Match the validated CLAS12/Gepard convention used by the existing
    extract_emff_from_dvcs_bh.py backend.
    """
    phi_rad = math.radians(float(row.phi_deg))
    phi_trento = math.pi - phi_rad
    pt = g.DataPoint(
        xB=float(row.xB),
        t=-abs(float(row.t_abs)),
        Q2=float(row.Q2),
        phi=float(phi_trento),
        observable="XS",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(row.ebeam),
        in1charge=-1,
        in1polarization=0,
        in2particle="p",
    )
    pt.prepare()
    return pt
#enddef


def cff_owner(th, name: str):
    """
    Find the object whose CFF method is actually consumed by the DVCS formula.

    Gepard theories normally expose `th.m` as the model used internally.  In
    common KM fits th and th.m can also effectively expose the same CFF API.
    Prefer th.m, because existing validated code accesses elastic FFs there.
    """
    candidates = []
    if hasattr(th, "m"):
        candidates.append(th.m)
    #endif
    candidates.append(th)

    seen = set()
    for obj in candidates:
        if id(obj) in seen:
            continue
        #endif
        seen.add(id(obj))
        if hasattr(obj, name) and callable(getattr(obj, name)):
            return obj
        #endif
    #endfor
    raise AttributeError(f"Could not find callable Gepard CFF method {name}")
#enddef


@contextlib.contextmanager
def shifted_cff(th, name: str, shift: float):
    """Temporarily replace CFF(pt) -> CFF_nominal(pt) + shift."""
    owner = cff_owner(th, name)
    old = getattr(owner, name)

    # Instance-assigned callable receives only the explicit DataPoint argument.
    setattr(owner, name, lambda pt, _old=old, _shift=float(shift):
            float(_old(pt)) + _shift)
    try:
        yield
    finally:
        setattr(owner, name, old)
    #endtry
#enddef


def prediction(th, pt) -> float:
    return float(th.predict(pt))
#enddef


def cff_value(th, name: str, pt) -> float:
    owner = cff_owner(th, name)
    return float(getattr(owner, name)(pt))
#enddef


def finite_cff_derivative(th, pts: Sequence, name: str,
                          nominal_value: float,
                          rel_step: float,
                          abs_step: float) -> Tuple[np.ndarray, float]:
    """
    d sigma / d CFF from symmetric finite changes around KM15.
    """
    step = max(abs_step, rel_step * max(abs(nominal_value), 1.0))

    with shifted_cff(th, name, +step):
        plus = np.asarray([prediction(th, pt) for pt in pts], dtype=float)
    #endwith
    with shifted_cff(th, name, -step):
        minus = np.asarray([prediction(th, pt) for pt in pts], dtype=float)
    #endwith

    deriv = (plus - minus) / (2.0 * step)
    return deriv, step
#enddef


def validate_cff_hook(th, pt, name: str, rel_step: float,
                      abs_step: float) -> Dict[str, float]:
    nominal = prediction(th, pt)
    cff = cff_value(th, name, pt)
    step = max(abs_step, rel_step * max(abs(cff), 1.0))
    with shifted_cff(th, name, step):
        shifted = prediction(th, pt)
    #endwith
    rel_response = (shifted - nominal) / max(abs(nominal), 1e-30)
    return {
        "cff": cff,
        "step": step,
        "sigma": nominal,
        "sigma_shifted": shifted,
        "relative_response": rel_response,
    }
#enddef


# =============================================================================
# Input and cell derivative cache
# =============================================================================

def load_inputs(indir: Path) -> Dict[int, pd.DataFrame]:
    out = {}
    for factor in LUMI_FACTORS:
        path = indir / f"cff_fit_input_km15_{factor}x.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing Stage-2A input: {path}")
        #endif
        df = pd.read_csv(path)
        required = {
            "point_id", "bin", "xB", "Q2", "t_abs", "phi_deg", "ebeam",
            "xs_pseudo", "stat_pseudo_abs", "ptp_sys_pseudo_abs",
            "correlated_norm_frac", "km15_ep",
        }
        missing = required - set(df.columns)
        if missing:
            raise KeyError(f"{path} missing columns: {sorted(missing)}")
        #endif
        df["bin"] = pd.to_numeric(df["bin"], errors="raise").astype(int)
        out[factor] = df.sort_values(["bin", "phi_deg"]).reset_index(drop=True)
    #endfor

    # All luminosity scenarios must contain exactly the same kinematics.
    reference = out[1][["point_id", "bin"]].astype(str)
    for factor in LUMI_FACTORS[1:]:
        test = out[factor][["point_id", "bin"]].astype(str)
        if not reference.equals(test):
            raise RuntimeError(f"{factor}x input does not match 1x point ordering")
        #endif
    #endfor
    return out
#enddef


def build_derivative_table(df: pd.DataFrame, th, g,
                           rel_step: float, abs_step: float,
                           min_nphi: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Evaluate nominal CFFs and finite-change XS sensitivities once per cell.
    Derivatives do not depend on the luminosity scenario.
    """
    point_rows = []
    cell_rows = []

    grouped = list(df.groupby("bin", sort=True))
    print(f"[CFF] calculating finite-change sensitivities for {len(grouped)} cells")

    for icell, (bin_id, cell) in enumerate(grouped, start=1):
        cell = cell.sort_values("phi_deg")
        if len(cell) < min_nphi:
            print(f"[CFF] bin {bin_id}: skipped, only {len(cell)} phi points")
            continue
        #endif

        pts = [make_point(g, row) for row in cell.itertuples(index=False)]
        pt0 = pts[0]

        reh = cff_value(th, "ReH", pt0)
        imh = cff_value(th, "ImH", pt0)
        d_reh, step_reh = finite_cff_derivative(
            th, pts, "ReH", reh, rel_step, abs_step
        )
        d_imh, step_imh = finite_cff_derivative(
            th, pts, "ImH", imh, rel_step, abs_step
        )

        nominal_direct = np.asarray([prediction(th, pt) for pt in pts])
        nominal_input = cell["km15_ep"].to_numpy(float)
        closure = np.max(
            np.abs(nominal_direct - nominal_input)
            / np.maximum(np.abs(nominal_input), 1e-30)
        )
        if closure > 2e-6:
            raise RuntimeError(
                f"bin {bin_id}: direct Gepard prediction disagrees with "
                f"Stage-2A KM15 cache; max relative difference={closure:.3e}"
            )
        #endif

        for j, row in enumerate(cell.itertuples(index=False)):
            point_rows.append({
                "point_id": row.point_id,
                "bin": int(bin_id),
                "phi_deg": float(row.phi_deg),
                "d_sigma_d_ReH": float(d_reh[j]),
                "d_sigma_d_ImH": float(d_imh[j]),
            })
        #endfor

        # A useful scale-free measure: fractional XS response to a 1% CFF change.
        sigma = np.maximum(np.abs(nominal_input), 1e-30)
        reh_response = np.median(np.abs(d_reh * reh / sigma)) if reh != 0 else 0.0
        imh_response = np.median(np.abs(d_imh * imh / sigma)) if imh != 0 else 0.0

        cell_rows.append({
            "bin": int(bin_id),
            "xB": float(np.median(cell["xB"])),
            "Q2": float(np.median(cell["Q2"])),
            "t_abs": float(np.median(cell["t_abs"])),
            "t_over_Q2": float(np.median(cell["t_over_Q2"])) if "t_over_Q2" in cell else np.nan,
            "n_phi": int(len(cell)),
            "ReH_KM15": reh,
            "ImH_KM15": imh,
            "finite_step_ReH": step_reh,
            "finite_step_ImH": step_imh,
            "median_abs_dlnsigma_dlnReH": float(reh_response),
            "median_abs_dlnsigma_dlnImH": float(imh_response),
            "direct_cache_max_relerr": float(closure),
        })

        if icell % 10 == 0 or icell == len(grouped):
            print(f"[CFF] {icell:3d}/{len(grouped)} cells")
        #endif
    #endfor

    return pd.DataFrame(point_rows), pd.DataFrame(cell_rows)
#enddef


# =============================================================================
# Global uncertainty propagation
# =============================================================================

def parameter_layout(cell_meta: pd.DataFrame, mode: str):
    bins = cell_meta["bin"].astype(int).tolist()
    names = []
    for b in bins:
        names.append((b, "ReH"))
        if mode == "ReH_ImH":
            names.append((b, "ImH"))
        #endif
    #endfor
    names.append((-1, "norm_beta"))
    return names
#enddef


def covariance_for_scenario(df: pd.DataFrame,
                            derivatives: pd.DataFrame,
                            cell_meta: pd.DataFrame,
                            mode: str) -> Tuple[np.ndarray, List[Tuple[int, str]], Dict]:
    """
    Build the uncertainty matrix from finite observable changes.

    For each measured point:
       delta sigma ~= (d sigma/d ReH) delta ReH
                    + (d sigma/d ImH) delta ImH
                    + sigma * norm_frac * beta

    beta has a unit Gaussian prior, so beta=1 means one quoted correlated
    normalization uncertainty.

    The inverse of the weighted finite-change matrix gives the projected
    covariance. A pseudo-inverse is used only when necessary; rank diagnostics
    are exported so unconstrained cases are visible rather than hidden.
    """
    use = df.merge(derivatives, on=["point_id", "bin", "phi_deg"], how="inner")
    valid_bins = set(cell_meta["bin"].astype(int))
    use = use[use["bin"].isin(valid_bins)].copy().reset_index(drop=True)

    layout = parameter_layout(cell_meta, mode)
    index = {key: i for i, key in enumerate(layout)}
    npar = len(layout)
    npts = len(use)

    A = np.zeros((npts, npar), dtype=float)
    sigma_point = np.hypot(
        use["stat_pseudo_abs"].to_numpy(float),
        use["ptp_sys_pseudo_abs"].to_numpy(float),
    )
    y0 = use["xs_pseudo"].to_numpy(float)

    for i, row in enumerate(use.itertuples(index=False)):
        b = int(row.bin)
        A[i, index[(b, "ReH")]] = float(row.d_sigma_d_ReH)
        if mode == "ReH_ImH":
            A[i, index[(b, "ImH")]] = float(row.d_sigma_d_ImH)
        #endif
        A[i, index[(-1, "norm_beta")]] = (
            float(row.correlated_norm_frac) * float(row.xs_pseudo)
        )
    #endfor

    good = np.isfinite(sigma_point) & (sigma_point > 0)
    if not np.all(good):
        A = A[good]
        sigma_point = sigma_point[good]
    #endif

    Aw = A / sigma_point[:, None]

    # Add beta^2 prior as one extra unit-weight row.
    prior = np.zeros((1, npar))
    prior[0, index[(-1, "norm_beta")]] = 1.0
    B = np.vstack([Aw, prior])

    info = B.T @ B
    rank = int(np.linalg.matrix_rank(B))
    singular = np.linalg.svd(B, compute_uv=False)
    positive = singular[singular > max(singular[0], 1.0) * 1e-12] if len(singular) else np.array([])
    condition = (
        float(positive[0] / positive[-1])
        if len(positive) >= 2 else np.inf
    )

    cov = np.linalg.pinv(info, rcond=1e-12)

    diagnostics = {
        "n_points": int(A.shape[0]),
        "n_parameters": int(npar),
        "rank": rank,
        "rank_deficit": int(npar - rank),
        "condition_number_effective": condition,
        "norm_beta_sigma": float(math.sqrt(max(cov[index[(-1, "norm_beta")],
                                                   index[(-1, "norm_beta")]], 0.0))),
    }
    return cov, layout, diagnostics
#enddef


def summarize_scenario(df: pd.DataFrame, cell_meta: pd.DataFrame,
                       cov: np.ndarray, layout: List[Tuple[int, str]],
                       factor: int, mode: str) -> pd.DataFrame:
    idx = {key: i for i, key in enumerate(layout)}
    rows = []

    for cell in cell_meta.itertuples(index=False):
        b = int(cell.bin)
        i_re = idx[(b, "ReH")]
        sig_re = math.sqrt(max(float(cov[i_re, i_re]), 0.0))

        row = {
            "luminosity_factor": factor,
            "fit_mode": mode,
            "bin": b,
            "xB": float(cell.xB),
            "Q2": float(cell.Q2),
            "t_abs": float(cell.t_abs),
            "t_over_Q2": float(cell.t_over_Q2),
            "n_phi": int(cell.n_phi),
            "ReH_KM15": float(cell.ReH_KM15),
            "ImH_KM15": float(cell.ImH_KM15),
            "sigma_ReH": sig_re,
            "relative_sigma_ReH": (
                sig_re / abs(float(cell.ReH_KM15))
                if abs(float(cell.ReH_KM15)) > 1e-12 else np.nan
            ),
            "sigma_ImH": np.nan,
            "relative_sigma_ImH": np.nan,
            "corr_ReH_ImH": np.nan,
            "corr_ReH_norm": np.nan,
        }

        i_norm = idx[(-1, "norm_beta")]
        denom = math.sqrt(max(cov[i_re, i_re] * cov[i_norm, i_norm], 0.0))
        if denom > 0:
            row["corr_ReH_norm"] = float(cov[i_re, i_norm] / denom)
        #endif

        if mode == "ReH_ImH":
            i_im = idx[(b, "ImH")]
            sig_im = math.sqrt(max(float(cov[i_im, i_im]), 0.0))
            row["sigma_ImH"] = sig_im
            row["relative_sigma_ImH"] = (
                sig_im / abs(float(cell.ImH_KM15))
                if abs(float(cell.ImH_KM15)) > 1e-12 else np.nan
            )
            denom = math.sqrt(max(cov[i_re, i_re] * cov[i_im, i_im], 0.0))
            if denom > 0:
                row["corr_ReH_ImH"] = float(cov[i_re, i_im] / denom)
            #endif
        #endif
        rows.append(row)
    #endfor

    return pd.DataFrame(rows)
#enddef


# =============================================================================
# Diagnostics and plots
# =============================================================================

def savefig(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)
#enddef


def plot_cff_truth(cell_meta: pd.DataFrame, figures: Path):
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    ax.scatter(cell_meta["t_abs"], cell_meta["ReH_KM15"], s=32, alpha=0.75)
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"KM15 $\mathrm{Re}\,\mathcal{H}$")
    ax.set_title(r"KM15 $\mathrm{Re}\,\mathcal{H}$ at clean CLAS12 cells")
    ax.grid(alpha=0.2)
    savefig(fig, figures / "01_km15_ReH_cells.png")

    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    ax.scatter(cell_meta["t_abs"], cell_meta["ImH_KM15"], s=32, alpha=0.75)
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"KM15 $\mathrm{Im}\,\mathcal{H}$")
    ax.set_title(r"KM15 $\mathrm{Im}\,\mathcal{H}$ at clean CLAS12 cells")
    ax.grid(alpha=0.2)
    savefig(fig, figures / "02_km15_ImH_cells.png")
#enddef


def plot_response(cell_meta: pd.DataFrame, figures: Path):
    fig, ax = plt.subplots(figsize=(8.5, 5.9))
    ax.scatter(
        cell_meta["t_abs"],
        100.0 * cell_meta["median_abs_dlnsigma_dlnReH"],
        s=34, alpha=0.75, label=r"$\mathrm{Re}\,\mathcal{H}$",
    )
    ax.scatter(
        cell_meta["t_abs"],
        100.0 * cell_meta["median_abs_dlnsigma_dlnImH"],
        s=34, alpha=0.75, label=r"$\mathrm{Im}\,\mathcal{H}$",
    )
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel("Median cross-section change for a 100% CFF change (%)")
    ax.set_title("Direct finite-change CFF sensitivity of the CLAS12 cross section")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "03_cross_section_cff_response_vs_t.png")
#enddef


def plot_reh_uncertainties(results: pd.DataFrame, figures: Path):
    for mode in FIT_MODES:
        dmode = results[results["fit_mode"] == mode]
        fig, ax = plt.subplots(figsize=(8.7, 6.0))
        for factor in LUMI_FACTORS:
            d = dmode[dmode["luminosity_factor"] == factor]
            ax.scatter(
                d["t_abs"], d["sigma_ReH"],
                s=30, alpha=0.72, label=f"{factor}x",
            )
        #endfor
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(r"Projected uncertainty on $\mathrm{Re}\,\mathcal{H}$")
        title = (
            r"$\mathrm{Re}\,\mathcal{H}$ only; other CFFs fixed to KM15"
            if mode == "ReH_only"
            else r"$\mathrm{Re}\,\mathcal{H}$ and $\mathrm{Im}\,\mathcal{H}$ free per cell"
        )
        ax.set_title(title)
        ax.grid(alpha=0.2)
        ax.legend(title="Luminosity")
        savefig(fig, figures / f"04_ReH_uncertainty_{mode}.png")

        # Relative uncertainties are useful but unstable near a zero crossing.
        rel = dmode[np.isfinite(dmode["relative_sigma_ReH"])].copy()
        if len(rel):
            fig, ax = plt.subplots(figsize=(8.7, 6.0))
            for factor in LUMI_FACTORS:
                d = rel[rel["luminosity_factor"] == factor]
                ax.scatter(
                    d["t_abs"], 100.0 * d["relative_sigma_ReH"],
                    s=30, alpha=0.72, label=f"{factor}x",
                )
            #endfor
            ax.set_xlabel(r"$|t|$ (GeV$^2$)")
            ax.set_ylabel(r"Relative uncertainty on $\mathrm{Re}\,\mathcal{H}$ (%)")
            ax.set_title(title)
            ax.grid(alpha=0.2)
            ax.legend(title="Luminosity")
            savefig(fig, figures / f"05_relative_ReH_uncertainty_{mode}.png")
        #endif
    #endfor
#enddef


def plot_imh_uncertainties(results: pd.DataFrame, figures: Path):
    dmode = results[results["fit_mode"] == "ReH_ImH"]
    fig, ax = plt.subplots(figsize=(8.7, 6.0))
    for factor in LUMI_FACTORS:
        d = dmode[dmode["luminosity_factor"] == factor]
        ax.scatter(
            d["t_abs"], d["sigma_ImH"],
            s=30, alpha=0.72, label=f"{factor}x",
        )
    #endfor
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"Projected uncertainty on $\mathrm{Im}\,\mathcal{H}$")
    ax.set_title(
        r"$\mathrm{Re}\,\mathcal{H}$ and $\mathrm{Im}\,\mathcal{H}$ free per cell"
    )
    ax.grid(alpha=0.2)
    ax.legend(title="Luminosity")
    savefig(fig, figures / "06_ImH_uncertainty_ReH_ImH.png")
#enddef


def plot_correlations(results: pd.DataFrame, figures: Path):
    d = results[
        (results["fit_mode"] == "ReH_ImH")
        & (results["luminosity_factor"] == 10)
    ].copy()
    fig, ax = plt.subplots(figsize=(8.5, 5.9))
    ax.scatter(d["t_abs"], d["corr_ReH_ImH"], s=34, alpha=0.75)
    ax.axhline(0.0, linewidth=1.0, linestyle="--")
    ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"Correlation of $\mathrm{Re}\,\mathcal{H}$ and $\mathrm{Im}\,\mathcal{H}$")
    ax.set_title("CFF degeneracy at 10x luminosity")
    ax.grid(alpha=0.2)
    savefig(fig, figures / "07_ReH_ImH_correlation_10x.png")
#enddef


def make_t_summary(results: pd.DataFrame) -> pd.DataFrame:
    edges = np.arange(0.1, 1.01, 0.1)
    d = results.copy()
    d["t_interval"] = pd.cut(d["t_abs"], edges, right=False, include_lowest=True)
    rows = []
    for (mode, factor, interval), g in d.groupby(
            ["fit_mode", "luminosity_factor", "t_interval"], observed=False):
        if len(g) == 0:
            continue
        #endif
        rows.append({
            "fit_mode": mode,
            "luminosity_factor": int(factor),
            "t_low": float(interval.left),
            "t_high": float(interval.right),
            "n_cells": int(len(g)),
            "median_sigma_ReH": float(np.nanmedian(g["sigma_ReH"])),
            "median_relative_sigma_ReH_pct": float(
                100.0 * np.nanmedian(g["relative_sigma_ReH"])
            ),
            "median_sigma_ImH": float(np.nanmedian(g["sigma_ImH"]))
                if np.any(np.isfinite(g["sigma_ImH"])) else np.nan,
            "median_abs_corr_ReH_ImH": float(
                np.nanmedian(np.abs(g["corr_ReH_ImH"]))
            ) if np.any(np.isfinite(g["corr_ReH_ImH"])) else np.nan,
        })
    #endfor
    return pd.DataFrame(rows)
#enddef


def print_summary(results: pd.DataFrame, diagnostics: pd.DataFrame,
                  cell_meta: pd.DataFrame):
    print("\n" + "=" * 86)
    print("STAGE-2B LOCAL CFF-H SENSITIVITY SUMMARY")
    print("=" * 86)
    print(f"Usable clean (xB,Q2,t) cells: {len(cell_meta)}")
    print(
        f"|t| range: {cell_meta['t_abs'].min():.3f} -- "
        f"{cell_meta['t_abs'].max():.3f} GeV^2"
    )

    print("\nGlobal finite-change fit diagnostics:")
    print(
        diagnostics.to_string(
            index=False, float_format=lambda x: f"{x:.4g}"
        )
    )

    print("\nMedian projected ReH uncertainty across cells:")
    rows = []
    for mode in FIT_MODES:
        for factor in LUMI_FACTORS:
            d = results[
                (results["fit_mode"] == mode)
                & (results["luminosity_factor"] == factor)
            ]
            rows.append({
                "mode": mode,
                "L": f"{factor}x",
                "median_sigma_ReH": np.nanmedian(d["sigma_ReH"]),
                "median_relative_ReH_pct": 100*np.nanmedian(d["relative_sigma_ReH"]),
                "cells_rel_ReH_lt_50pct": int(
                    np.sum(d["relative_sigma_ReH"] < 0.50)
                ),
                "cells_rel_ReH_lt_25pct": int(
                    np.sum(d["relative_sigma_ReH"] < 0.25)
                ),
            })
        #endfor
    #endfor
    print(
        pd.DataFrame(rows).to_string(
            index=False, float_format=lambda x: f"{x:.3g}"
        )
    )

    print("\nInterpretation guardrails:")
    print("  ReH_only : optimistic; ImH and all non-H CFFs fixed to KM15.")
    print("  ReH_ImH  : ReH and ImH free independently in each cell; non-H CFFs fixed.")
    print("  The shared normalization nuisance is included in both.")
    print("  Large uncertainties/correlations are physics information, not fit failures.")
    print("  Do not proceed to D(t)/pressure until the ReH lever arm is demonstrably useful.")
#enddef


# =============================================================================
# Main
# =============================================================================

def parser():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(
        description="Finite-change CLAS12 high-luminosity CFF-H sensitivity projection"
    )
    p.add_argument(
        "--input-dir",
        default=str(here / "output_stage2" / "tables"),
    )
    p.add_argument(
        "--outdir",
        default=str(here / "output_stage2_cff"),
    )
    p.add_argument("--min-nphi", type=int, default=DEFAULT_MIN_NPHI)
    p.add_argument("--relative-step", type=float, default=DEFAULT_REL_STEP)
    p.add_argument("--absolute-step", type=float, default=DEFAULT_ABS_STEP)
    return p
#enddef


def main(argv: Optional[List[str]] = None) -> int:
    args = parser().parse_args(argv)
    indir = Path(args.input_dir).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    inputs = load_inputs(indir)
    one = inputs[1]

    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError(
            "Could not import Gepard/KM15. Run in the same Python environment "
            "used successfully by high_luminosity_dvcs_stage2_pseudodata.py."
        ) from exc
    #endtry

    th = th_KM15

    # Before doing hundreds of evaluations, prove that the CFF hook actually
    # changes the observable in this installed Gepard version.
    test_row = one.iloc[len(one)//2]
    test_pt = make_point(g, test_row)
    print("[preflight] testing direct CFF perturbations in installed Gepard")
    for name in ["ReH", "ImH"]:
        check = validate_cff_hook(
            th, test_pt, name,
            float(args.relative_step), float(args.absolute_step)
        )
        print(
            f"[preflight] {name}: value={check['cff']:.6g}, "
            f"step={check['step']:.6g}, "
            f"XS response={100*check['relative_response']:+.6g}%"
        )
        if abs(check["relative_response"]) < 1e-10:
            raise RuntimeError(
                f"Perturbing {name} did not change the cross section. "
                "The installed Gepard object wiring differs from the expected "
                "interface; stop rather than producing meaningless projections."
            )
        #endif
    #endfor

    derivatives, cell_meta = build_derivative_table(
        one, th, g,
        float(args.relative_step),
        float(args.absolute_step),
        int(args.min_nphi),
    )
    derivatives.to_csv(tables / "point_cff_derivatives.csv", index=False)
    cell_meta.to_csv(tables / "cell_km15_cffs_and_responses.csv", index=False)

    result_frames = []
    diag_rows = []
    covariance_dir = tables / "covariances"
    covariance_dir.mkdir(exist_ok=True)

    for mode in FIT_MODES:
        for factor in LUMI_FACTORS:
            print(f"[fit] mode={mode:9s} luminosity={factor}x")
            cov, layout, diag = covariance_for_scenario(
                inputs[factor], derivatives, cell_meta, mode
            )
            result_frames.append(
                summarize_scenario(
                    inputs[factor], cell_meta, cov, layout, factor, mode
                )
            )
            diag.update({
                "fit_mode": mode,
                "luminosity_factor": factor,
            })
            diag_rows.append(diag)

            # Save covariance and exact parameter ordering for reproducibility.
            np.save(
                covariance_dir / f"covariance_{mode}_{factor}x.npy",
                cov,
            )
            pd.DataFrame(
                [{"index": i, "bin": b, "parameter": p}
                 for i, (b, p) in enumerate(layout)]
            ).to_csv(
                covariance_dir / f"parameter_order_{mode}_{factor}x.csv",
                index=False,
            )
        #endfor
    #endfor

    results = pd.concat(result_frames, ignore_index=True)
    diagnostics = pd.DataFrame(diag_rows)
    t_summary = make_t_summary(results)

    results.to_csv(tables / "cff_cell_uncertainties.csv", index=False)
    diagnostics.to_csv(tables / "global_fit_diagnostics.csv", index=False)
    t_summary.to_csv(tables / "cff_t_luminosity_summary.csv", index=False)

    plot_cff_truth(cell_meta, figures)
    plot_response(cell_meta, figures)
    plot_reh_uncertainties(results, figures)
    plot_imh_uncertainties(results, figures)
    plot_correlations(results, figures)

    print_summary(results, diagnostics, cell_meta)
    print(f"\n[output] {outdir}")
    print("[next] inspect global_fit_diagnostics.csv first.")
    print("[next] then inspect ReH uncertainties and ReH-ImH correlations versus |t|.")
    print("[next] only after those pass should we build the dispersion-relation/D(t) layer.")
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
