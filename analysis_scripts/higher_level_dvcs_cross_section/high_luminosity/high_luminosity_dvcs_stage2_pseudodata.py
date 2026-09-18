#!/usr/bin/env python3
"""
Stage 2: CLAS12 DVCS high-luminosity pseudo-data preparation and KM15 closure.

Purpose
-------
This script is the bridge between the Stage-1 kinematic-reach study and the
eventual CFF / subtraction-constant extraction.

It:
  1. loads the authoritative released CLAS12 RGA DVCS cross sections;
  2. evaluates KM15 (and pure BH) at the exact published kinematics using the
     already validated Gepard backend in extract_emff_from_dvcs_bh.py;
  3. applies the preferred |t|/Q^2 < 0.2 theory-control cut;
  4. diagnoses whether KM15 is a reasonable pseudo-truth specifically in that
     clean region and in its high-|t| tail;
  5. writes deterministic Asimov pseudo-data at 1x, 2x, 5x, and 10x luminosity:
         central value = KM15 prediction
         stat(L)      = published stat / sqrt(L)
         ptp syst     = fixed at the published point-to-point value
         correlated normalization = fixed at the published fractional prior
  6. writes point-level and cell-level tables for the next CFF-fit stage.

This script deliberately DOES NOT fit CFFs yet.  The first requirement is to
establish that the chosen pseudo-truth behaves reasonably in the exact
high-Q^2/high-|t| subset that drives the mechanical-structure argument.

Expected repository layout
--------------------------
analysis_scripts/
  dvcs_cross_section/
    external_scripts/
      extract_emff_from_dvcs_bh.py
  higher_level_dvcs_cross_section/
    import/
      clasdb_E214M1.txt
    high_luminosity/
      high_luminosity_dvcs_stage2_pseudodata.py

Typical use
-----------
From higher_level_dvcs_cross_section/high_luminosity/:

  python high_luminosity_dvcs_stage2_pseudodata.py

Force a fresh KM15 evaluation:

  python high_luminosity_dvcs_stage2_pseudodata.py --force-km15

The KM15 calculation is cached.  Re-running the script after the first
successful evaluation should therefore be fast.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)
DEFAULT_CLEAN_RATIO = 0.20
DEFAULT_HIGH_T = 0.40
DEFAULT_EBEAM = 10.604
DEFAULT_NORM_FRAC = 0.31
T_EDGES = np.arange(0.10, 1.01, 0.10)


# =============================================================================
# Utilities
# =============================================================================

def resolve_existing(path: Path, alternatives: Sequence[Path]) -> Path:
    candidates = [Path(path)] + [Path(p) for p in alternatives]
    for candidate in candidates:
        candidate = candidate.expanduser()
        if candidate.exists():
            return candidate.resolve()
        #endif
    #endfor
    attempted = "\n".join(f"  {p}" for p in candidates)
    raise FileNotFoundError(f"Required input not found. Tried:\n{attempted}")
#enddef


def load_python_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load Python module from {path}")
    #endif
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
#enddef


def numeric_column(df: pd.DataFrame, candidates: Sequence[str],
                   required: bool = True, default=np.nan) -> np.ndarray:
    for name in candidates:
        if name in df.columns:
            return pd.to_numeric(df[name], errors="coerce").to_numpy(float)
        #endif
    #endfor
    if required:
        raise KeyError(
            "None of the required columns are present: "
            + ", ".join(candidates)
            + f"\nAvailable columns: {list(df.columns)}"
        )
    #endif
    return np.full(len(df), float(default), dtype=float)
#enddef


def finite_positive(x) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    return np.isfinite(a) & (a > 0.0)
#enddef


def dataset_fingerprint(df: pd.DataFrame) -> str:
    cols = ["point_id", "xB", "Q2", "t_abs", "phi_deg", "ebeam"]
    payload = df[cols].to_csv(index=False, float_format="%.17g").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:12]
#enddef


# =============================================================================
# Data loading
# =============================================================================

def canonicalize_clas12(raw: pd.DataFrame, emff, fallback_ebeam: float) -> pd.DataFrame:
    """
    Convert the validated EMFF loader output to the compact Stage-2 schema.

    We intentionally use the EMFF loader rather than reparsing E214M1 here,
    because that loader already contains the validated treatment of the
    released 31% correlated normalization and residual point-to-point
    systematic uncertainty.
    """
    out = pd.DataFrame()

    out["xB"] = numeric_column(raw, ["xB", "xb"])
    out["Q2"] = numeric_column(raw, ["Q2", "q2"])
    out["t_abs"] = np.abs(numeric_column(raw, ["t_abs", "t", "t_average"]))
    out["phi_deg"] = np.mod(
        numeric_column(raw, ["phi_deg", "phi", "phi_average"]), 360.0
    )
    out["xs_data"] = numeric_column(
        raw, ["xs", "sigma", "cross_section", "d4sigma"]
    )
    out["stat_abs"] = numeric_column(
        raw, ["xs_stat", "stat", "stat_abs", "stat_err"]
    )

    # Preferred: residual point-to-point systematic after the correlated
    # normalization component has been removed by the validated loader.
    ptp = numeric_column(
        raw,
        ["ptp_sys_abs", "xs_syst_ptp", "syst_ptp", "pointwise_syst"],
        required=False,
    )

    # Keep the released total systematic too when the loader exposes it.
    released_syst = numeric_column(
        raw,
        ["xs_syst", "syst", "syst_abs", "released_syst_abs"],
        required=False,
    )

    # If ptp_sys_abs is unavailable but the released total systematic is
    # available, recover the residual component using the same quadrature
    # convention as the validated analysis.
    norm_frac = float(getattr(emff, "PASS1_GLOBAL_SCALE_FRAC", DEFAULT_NORM_FRAC))
    missing_ptp = ~np.isfinite(ptp)
    if np.any(missing_ptp) and np.any(np.isfinite(released_syst)):
        recovered = np.sqrt(
            np.maximum(
                released_syst**2 - (norm_frac * out["xs_data"].to_numpy(float))**2,
                0.0,
            )
        )
        ptp[missing_ptp] = recovered[missing_ptp]
    #endif

    if np.any(~np.isfinite(ptp)):
        raise RuntimeError(
            "Could not determine the CLAS12 point-to-point systematic "
            "uncertainty. Use the current validated extract_emff_from_dvcs_bh.py "
            "backend, which should expose ptp_sys_abs."
        )
    #endif

    out["ptp_sys_abs"] = ptp
    out["released_syst_abs"] = released_syst

    if "ebeam" in raw.columns:
        out["ebeam"] = pd.to_numeric(raw["ebeam"], errors="coerce").to_numpy(float)
    elif "Eb" in raw.columns:
        out["ebeam"] = pd.to_numeric(raw["Eb"], errors="coerce").to_numpy(float)
    else:
        out["ebeam"] = float(fallback_ebeam)
    #endif

    if "bin" in raw.columns:
        out["bin"] = pd.to_numeric(raw["bin"], errors="coerce").to_numpy(float)
    else:
        # E214M1 itself carries the bin index. This fallback should almost
        # never be needed with the validated loader.
        out["bin"] = np.arange(len(out), dtype=int)
    #endif

    finite = (
        finite_positive(out["xB"])
        & finite_positive(out["Q2"])
        & finite_positive(out["t_abs"])
        & np.isfinite(out["phi_deg"])
        & finite_positive(out["xs_data"])
        & finite_positive(out["stat_abs"])
        & np.isfinite(out["ptp_sys_abs"])
        & finite_positive(out["ebeam"])
    )
    if not np.all(finite):
        print(f"[data] dropping {int((~finite).sum())} invalid row(s)")
    #endif
    out = out.loc[finite].copy().reset_index(drop=True)

    # Preserve the released bin number as an integer after invalid-row removal.
    out["bin"] = np.rint(out["bin"]).astype(int)
    out["source_row"] = np.arange(len(out), dtype=int)
    out["point_id"] = [f"lee2026:{i}" for i in out["source_row"]]

    out["t_over_Q2"] = out["t_abs"] / out["Q2"]
    out["stat_frac_data"] = out["stat_abs"] / out["xs_data"]
    out["ptp_sys_frac_data"] = out["ptp_sys_abs"] / out["xs_data"]
    out["point_unc_abs_1x"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    out["point_unc_frac_1x"] = out["point_unc_abs_1x"] / out["xs_data"]
    out["norm_frac"] = norm_frac

    return out
#enddef


def load_clas12(args, emff) -> pd.DataFrame:
    xs_path = resolve_existing(
        Path(args.xs),
        [
            args.here.parent / "import" / "clasdb_E214M1.txt",
            args.external_scripts / "import" / "clasdb_E214M1.txt",
            args.external_scripts.parent / "imports" / "clasdb_E214M1.txt",
        ],
    )

    print(f"[data] loading authoritative CLAS12 cross sections: {xs_path}")
    raw = emff.load_clas12_pass1_csv(xs_path)
    data = canonicalize_clas12(raw, emff, args.ebeam)

    print(f"[data] loaded {len(data)} valid 4D points")
    print(f"[data] unique released (xB,Q2,t) bins: {data['bin'].nunique()}")
    print(f"[data] correlated normalization prior: {100*data['norm_frac'].iloc[0]:.2f}%")
    return data
#enddef


# =============================================================================
# KM15 evaluation
# =============================================================================

def evaluate_km15(data: pd.DataFrame, emff, cache_path: Path,
                  force: bool = False) -> pd.DataFrame:
    """
    Evaluate KM15 total electroproduction and pure BH at exact CLAS12 kinematics.

    CLAS12 E214M1 uses the direct/identity phi convention.  This is the same
    convention validated by the PARTONS/Gepard BH closure scan.
    """
    fingerprint = dataset_fingerprint(data)
    required = ["point_id", "km15_ep", "km15_bh"]

    if cache_path.exists() and not force:
        cache = pd.read_csv(cache_path)
        cache_fp = str(cache["fingerprint"].iloc[0]) if "fingerprint" in cache.columns and len(cache) else ""
        if (
            cache_fp == fingerprint
            and len(cache) == len(data)
            and set(required).issubset(cache.columns)
            and cache["point_id"].astype(str).tolist() == data["point_id"].astype(str).tolist()
        ):
            print(f"[KM15] reusing exact-kinematics cache: {cache_path}")
            out = data.copy()
            out["km15_ep"] = pd.to_numeric(cache["km15_ep"], errors="coerce").to_numpy(float)
            out["km15_bh"] = pd.to_numeric(cache["km15_bh"], errors="coerce").to_numpy(float)
            return finalize_model_columns(out)
        #endif
        print("[KM15] existing cache does not match current kinematics; recalculating")
    #endif

    rows = []
    print(f"[KM15] evaluating {len(data)} exact CLAS12 points")
    for i, row in enumerate(data.itertuples(index=False), start=1):
        task = (
            0,
            float(row.xB),
            float(row.Q2),
            float(row.t_abs),
            float(row.phi_deg),
            float(row.ebeam),
            "identity",
        )
        result = emff.evaluate_km15_point(task)
        rows.append({
            "point_id": str(row.point_id),
            "fingerprint": fingerprint,
            "km15_ep": float(result["km15_ep"]),
            "km15_bh": float(result["km15_bh"]),
        })
        if i % 100 == 0 or i == len(data):
            print(f"[KM15] {i:4d}/{len(data)}")
        #endif
    #endfor

    cache = pd.DataFrame(rows)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache.to_csv(cache_path, index=False)
    print(f"[KM15] cache -> {cache_path}")

    out = data.copy()
    out["km15_ep"] = cache["km15_ep"].to_numpy(float)
    out["km15_bh"] = cache["km15_bh"].to_numpy(float)
    return finalize_model_columns(out)
#enddef


def finalize_model_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["data_over_km15"] = out["xs_data"] / out["km15_ep"]
    out["data_over_bh"] = out["xs_data"] / out["km15_bh"]
    out["bh_fraction_of_km15"] = out["km15_bh"] / out["km15_ep"]
    out["fractional_residual_km15"] = (out["xs_data"] - out["km15_ep"]) / out["km15_ep"]
    out["residual_over_stat"] = (out["xs_data"] - out["km15_ep"]) / out["stat_abs"]
    out["residual_over_point_unc"] = (
        (out["xs_data"] - out["km15_ep"]) / out["point_unc_abs_1x"]
    )
    return out
#enddef


# =============================================================================
# Normalization-aware closure
# =============================================================================

def fit_one_normalization_nuisance(
        data: np.ndarray,
        model: np.ndarray,
        sigma: np.ndarray,
        norm_frac: float) -> Dict[str, float]:
    """
    Fit the one quoted correlated multiplicative normalization uncertainty.

    shifted_data = data * (1 + beta * norm_frac)

    Minimize:
      sum[(shifted_data-model)^2/sigma^2] + beta^2

    beta is therefore measured in units of the quoted normalization uncertainty.
    """
    y = np.asarray(data, dtype=float)
    m = np.asarray(model, dtype=float)
    s = np.asarray(sigma, dtype=float)

    good = (
        finite_positive(y) & finite_positive(m)
        & np.isfinite(s) & (s > 0.0)
    )
    y, m, s = y[good], m[good], s[good]

    if len(y) == 0:
        return {
            "N": 0, "beta": np.nan, "scale": np.nan,
            "chi2": np.nan, "chi2_per_point": np.nan,
            "median_abs_fractional_residual": np.nan,
            "rms_pull": np.nan,
        }
    #endif

    if norm_frac <= 0.0:
        beta = 0.0
    else:
        x = norm_frac * y
        w = 1.0 / s**2
        beta = -np.sum(w * x * (y - m)) / (1.0 + np.sum(w * x**2))
    #endif

    scale = 1.0 + beta * norm_frac
    residual = scale * y - m
    pull = residual / s
    chi2 = float(np.sum(pull**2) + beta**2)

    return {
        "N": int(len(y)),
        "beta": float(beta),
        "scale": float(scale),
        "chi2": chi2,
        "chi2_per_point": float(chi2 / len(y)),
        "median_abs_fractional_residual": float(
            np.nanmedian(np.abs(residual / m))
        ),
        "rms_pull": float(np.sqrt(np.nanmean(pull**2))),
    }
#enddef


def closure_row(d: pd.DataFrame, label: str, t_low=np.nan, t_high=np.nan) -> Dict[str, float]:
    point_unc = np.hypot(
        d["stat_abs"].to_numpy(float),
        d["ptp_sys_abs"].to_numpy(float),
    )
    result = fit_one_normalization_nuisance(
        d["xs_data"].to_numpy(float),
        d["km15_ep"].to_numpy(float),
        point_unc,
        float(d["norm_frac"].iloc[0]),
    )
    result.update({
        "region": label,
        "t_low": t_low,
        "t_high": t_high,
        "median_data_over_km15_raw": float(np.nanmedian(d["data_over_km15"])),
        "median_abs_raw_fractional_residual": float(
            np.nanmedian(np.abs(d["fractional_residual_km15"]))
        ),
        "median_stat_pct": float(100.0 * np.nanmedian(d["stat_frac_data"])),
        "median_ptp_sys_pct": float(100.0 * np.nanmedian(d["ptp_sys_frac_data"])),
        "xB_min": float(d["xB"].min()),
        "xB_max": float(d["xB"].max()),
        "Q2_min": float(d["Q2"].min()),
        "Q2_max": float(d["Q2"].max()),
        "t_abs_min": float(d["t_abs"].min()),
        "t_abs_max": float(d["t_abs"].max()),
    })
    return result
#enddef


def make_closure_summary(data: pd.DataFrame, clean_ratio: float,
                         high_t: float) -> pd.DataFrame:
    clean = data[data["t_over_Q2"] < clean_ratio].copy()
    rows = [
        closure_row(data, "all_released"),
        closure_row(clean, f"clean_ratio_lt_{clean_ratio:g}"),
    ]

    high = clean[clean["t_abs"] >= high_t]
    if len(high):
        rows.append(closure_row(high, f"clean_high_t_ge_{high_t:g}"))
    #endif

    temp = clean.copy()
    temp["t_interval"] = pd.cut(
        temp["t_abs"], T_EDGES, right=False, include_lowest=True
    )
    for interval, g in temp.groupby("t_interval", observed=False):
        if len(g) == 0:
            continue
        #endif
        rows.append(
            closure_row(
                g,
                f"clean_t_{interval.left:.1f}_{interval.right:.1f}",
                interval.left,
                interval.right,
            )
        )
    #endfor

    return pd.DataFrame(rows)
#enddef


# =============================================================================
# Pseudo-data
# =============================================================================

def build_pseudodata(clean: pd.DataFrame) -> pd.DataFrame:
    """
    Deterministic (Asimov) pseudo-data.

    The central values are KM15.  No random fluctuations are introduced at
    this stage, so differences between luminosity scenarios isolate the
    information gained from smaller statistical errors.
    """
    frames = []
    for factor in LUMI_FACTORS:
        d = clean.copy()
        d["luminosity_factor"] = factor
        d["xs_pseudo"] = d["km15_ep"]
        d["stat_pseudo_abs"] = d["stat_abs"] / math.sqrt(factor)

        # Keep the *fractional* point-to-point systematic measured in the
        # current data and apply that fractional floor to the pseudo-truth.
        # This avoids encoding the data/KM15 normalization mismatch into the
        # absolute pseudo systematic.
        d["ptp_sys_pseudo_abs"] = d["ptp_sys_frac_data"] * d["xs_pseudo"]
        d["point_unc_pseudo_abs"] = np.hypot(
            d["stat_pseudo_abs"], d["ptp_sys_pseudo_abs"]
        )
        d["stat_pseudo_frac"] = d["stat_pseudo_abs"] / d["xs_pseudo"]
        d["ptp_sys_pseudo_frac"] = d["ptp_sys_pseudo_abs"] / d["xs_pseudo"]
        d["point_unc_pseudo_frac"] = d["point_unc_pseudo_abs"] / d["xs_pseudo"]

        # Correlated scale uncertainty is intentionally NOT folded into each
        # point. It remains a separate nuisance prior for the eventual fit.
        d["correlated_norm_frac"] = d["norm_frac"]
        frames.append(d)
    #endfor

    return pd.concat(frames, ignore_index=True)
#enddef


def make_cell_summary(pseudo: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (factor, bin_id), g in pseudo.groupby(
            ["luminosity_factor", "bin"], sort=True):
        rows.append({
            "luminosity_factor": int(factor),
            "bin": int(bin_id),
            "xB": float(np.nanmedian(g["xB"])),
            "Q2": float(np.nanmedian(g["Q2"])),
            "t_abs": float(np.nanmedian(g["t_abs"])),
            "t_over_Q2": float(np.nanmedian(g["t_over_Q2"])),
            "n_phi": int(len(g)),
            "phi_min": float(g["phi_deg"].min()),
            "phi_max": float(g["phi_deg"].max()),
            "median_stat_pseudo_pct": float(
                100.0 * np.nanmedian(g["stat_pseudo_frac"])
            ),
            "p90_stat_pseudo_pct": float(
                100.0 * np.nanpercentile(g["stat_pseudo_frac"], 90)
            ),
            "median_point_unc_pseudo_pct": float(
                100.0 * np.nanmedian(g["point_unc_pseudo_frac"])
            ),
            "median_bh_fraction_of_km15": float(
                np.nanmedian(g["bh_fraction_of_km15"])
            ),
        })
    #endfor
    return pd.DataFrame(rows)
#enddef


def make_t_luminosity_summary(pseudo: pd.DataFrame) -> pd.DataFrame:
    temp = pseudo.copy()
    temp["t_interval"] = pd.cut(
        temp["t_abs"], T_EDGES, right=False, include_lowest=True
    )
    rows = []
    for (factor, interval), g in temp.groupby(
            ["luminosity_factor", "t_interval"], observed=False):
        if len(g) == 0:
            continue
        #endif
        rows.append({
            "luminosity_factor": int(factor),
            "t_low": float(interval.left),
            "t_high": float(interval.right),
            "n_4d_points": int(len(g)),
            "n_cells": int(g["bin"].nunique()),
            "median_stat_pct": float(
                100.0 * np.nanmedian(g["stat_pseudo_frac"])
            ),
            "median_ptp_sys_pct": float(
                100.0 * np.nanmedian(g["ptp_sys_pseudo_frac"])
            ),
            "median_point_unc_pct": float(
                100.0 * np.nanmedian(g["point_unc_pseudo_frac"])
            ),
        })
    #endfor
    return pd.DataFrame(rows)
#enddef


# =============================================================================
# Plots
# =============================================================================

def savefig(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)
#enddef


def plot_closure_vs_kinematics(clean: pd.DataFrame, figures: Path):
    specs = [
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
        ("phi_deg", r"$\phi$ (deg)"),
    ]
    for i, (col, xlabel) in enumerate(specs, start=1):
        fig, ax = plt.subplots(figsize=(8.4, 5.8))
        ax.errorbar(
            clean[col],
            clean["data_over_km15"],
            yerr=clean["point_unc_abs_1x"] / clean["km15_ep"],
            fmt="o",
            ms=3.0,
            lw=0.65,
            capsize=0,
            alpha=0.55,
        )
        ax.axhline(1.0, linestyle="--", linewidth=1.1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\sigma_{\rm data}/\sigma_{\rm KM15}$")
        ax.set_title(r"KM15 closure for CLAS12 points with $|t|/Q^2<0.2$")
        ax.grid(alpha=0.2)
        savefig(fig, figures / f"01_{i:02d}_clean_data_over_km15_vs_{col}.png")
    #endfor
#enddef


def plot_high_t_closure(clean: pd.DataFrame, high_t: float, figures: Path):
    d = clean[clean["t_abs"] >= high_t].copy()
    if len(d) == 0:
        return
    #endif

    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    sc = ax.scatter(
        d["t_abs"], d["data_over_km15"],
        c=d["Q2"], s=28, alpha=0.75,
    )
    ax.axhline(1.0, linestyle="--", linewidth=1.1)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$Q^2$ (GeV$^2$)")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$\sigma_{\rm data}/\sigma_{\rm KM15}$")
    ax.set_title(
        rf"High-$|t|$ KM15 closure inside $|t|/Q^2<0.2$ "
        rf"($|t|\geq{high_t:.1f}$ GeV$^2$)"
    )
    ax.grid(alpha=0.2)
    savefig(fig, figures / "02_high_t_data_over_km15.png")
#enddef


def plot_closure_t_summary(summary: pd.DataFrame, figures: Path):
    d = summary[np.isfinite(summary["t_low"])].copy()
    if len(d) == 0:
        return
    #endif
    x = 0.5 * (d["t_low"] + d["t_high"])

    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    ax.plot(
        x, d["median_data_over_km15_raw"],
        marker="o", linewidth=1.8,
    )
    ax.axhline(1.0, linestyle="--", linewidth=1.1)
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"Median $\sigma_{\rm data}/\sigma_{\rm KM15}$")
    ax.set_title(r"KM15 normalization/shape check versus $|t|$")
    ax.grid(alpha=0.2)
    savefig(fig, figures / "03_km15_closure_summary_vs_t.png")

    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    ax.plot(
        x, 100.0 * d["median_abs_raw_fractional_residual"],
        marker="o", linewidth=1.8,
    )
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel("Median absolute KM15 residual (%)")
    ax.set_title(r"Size of KM15 residuals in the $|t|/Q^2<0.2$ region")
    ax.grid(alpha=0.2)
    savefig(fig, figures / "04_km15_fractional_residual_vs_t.png")
#enddef


def plot_pseudo_precision(t_summary: pd.DataFrame, figures: Path):
    fig, ax = plt.subplots(figsize=(8.6, 5.9))
    for factor in LUMI_FACTORS:
        d = t_summary[t_summary["luminosity_factor"] == factor]
        x = 0.5 * (d["t_low"] + d["t_high"])
        ax.plot(
            x, d["median_stat_pct"],
            marker="o", linewidth=1.8,
            label=f"{factor}x luminosity",
        )
    #endfor
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel("Median statistical uncertainty on pseudo-data (%)")
    ax.set_title(r"KM15 pseudo-data inside $|t|/Q^2<0.2$")
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "05_pseudodata_statistical_precision_vs_t.png")

    fig, ax = plt.subplots(figsize=(8.6, 5.9))
    for factor in LUMI_FACTORS:
        d = t_summary[t_summary["luminosity_factor"] == factor]
        x = 0.5 * (d["t_low"] + d["t_high"])
        ax.plot(
            x, d["median_point_unc_pct"],
            marker="o", linewidth=1.8,
            label=f"{factor}x luminosity",
        )
    #endfor
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel("Median stat + point-to-point systematic uncertainty (%)")
    ax.set_title(
        r"Projected pointwise precision with current systematic floor"
    )
    ax.grid(alpha=0.2)
    ax.legend()
    savefig(fig, figures / "06_pseudodata_point_uncertainty_vs_t.png")
#enddef


def plot_selected_phi_cells(clean: pd.DataFrame, figures: Path, n_cells: int = 6):
    """
    Show the highest-|t| clean cells.  This is a direct shape diagnostic:
    does KM15 reproduce the measured phi dependence where the mechanical
    lever arm is largest?
    """
    cells = (
        clean.groupby("bin", as_index=False)
        .agg(
            xB=("xB", "median"),
            Q2=("Q2", "median"),
            t_abs=("t_abs", "median"),
            n_phi=("phi_deg", "size"),
        )
        .sort_values(["t_abs", "n_phi"], ascending=[False, False])
    )
    chosen = cells.head(n_cells)

    for rank, cell in enumerate(chosen.itertuples(index=False), start=1):
        d = clean[clean["bin"] == int(cell.bin)].sort_values("phi_deg")
        fig, ax = plt.subplots(figsize=(8.2, 5.6))
        ax.errorbar(
            d["phi_deg"], d["xs_data"],
            yerr=d["point_unc_abs_1x"],
            fmt="o", ms=5, capsize=2,
            label="CLAS12 data",
        )
        ax.plot(
            d["phi_deg"], d["km15_ep"],
            marker="s", linewidth=1.6,
            label="KM15",
        )
        ax.plot(
            d["phi_deg"], d["km15_bh"],
            linestyle="--", linewidth=1.3,
            label="Pure BH",
        )
        ax.set_xlabel(r"$\phi$ (deg)")
        ax.set_ylabel(r"$d^4\sigma$ (published units)")
        ax.set_title(
            rf"bin {int(cell.bin)}: $x_B={cell.xB:.3f}$, "
            rf"$Q^2={cell.Q2:.3f}$ GeV$^2$, "
            rf"$|t|={cell.t_abs:.3f}$ GeV$^2$"
        )
        ax.grid(alpha=0.2)
        ax.legend()
        savefig(fig, figures / f"07_{rank:02d}_high_t_phi_bin_{int(cell.bin)}.png")
    #endfor
#enddef


# =============================================================================
# Reporting
# =============================================================================

def print_summary(data: pd.DataFrame, clean: pd.DataFrame,
                  closure: pd.DataFrame, clean_ratio: float):
    print("\n" + "=" * 80)
    print("STAGE-2 CLAS12 / KM15 SUMMARY")
    print("=" * 80)
    print(f"Released 4D points                  : {len(data)}")
    print(f"Released (xB,Q2,t) cells            : {data['bin'].nunique()}")
    print(
        f"|t|/Q2 < {clean_ratio:.2f} 4D points          : "
        f"{len(clean)} ({100*len(clean)/len(data):.1f}%)"
    )
    print(
        f"|t|/Q2 < {clean_ratio:.2f} cells              : "
        f"{clean['bin'].nunique()}"
    )
    print(
        "Clean-region |t| range (GeV^2)       : "
        f"{clean['t_abs'].min():.3f} -- {clean['t_abs'].max():.3f}"
    )
    print(
        "Clean-region Q2 range (GeV^2)        : "
        f"{clean['Q2'].min():.3f} -- {clean['Q2'].max():.3f}"
    )

    cols = [
        "region", "N", "beta", "scale", "chi2_per_point",
        "median_data_over_km15_raw",
        "median_abs_raw_fractional_residual",
    ]
    print("\nKM15 closure (pointwise stat+ptp systematic; one correlated norm nuisance)")
    print(
        closure[cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print("\nPseudo-data convention:")
    print("  central values             = KM15 exact-kinematics prediction")
    print("  statistical uncertainty   = current published stat / sqrt(L/L0)")
    print("  point-to-point systematic = current fractional ptp systematic, fixed")
    print(
        f"  correlated normalization  = "
        f"{100*clean['norm_frac'].iloc[0]:.2f}% nuisance, fixed"
    )
    print("  no random fluctuations are applied in this first Asimov study")
#enddef


# =============================================================================
# CLI / main
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    here = Path(__file__).resolve().parent
    external = here.parent.parent / "dvcs_cross_section" / "external_scripts"

    p = argparse.ArgumentParser(
        description="CLAS12 high-luminosity Stage-2 KM15 pseudo-data preparation"
    )
    p.add_argument(
        "--xs",
        default=str(here.parent / "import" / "clasdb_E214M1.txt"),
        help="Authoritative released CLAS12 E214M1 cross-section table",
    )
    p.add_argument(
        "--emff-script",
        default=str(external / "extract_emff_from_dvcs_bh.py"),
        help="Validated EMFF/KM15 backend",
    )
    p.add_argument(
        "--outdir",
        default=str(here / "output_stage2"),
    )
    p.add_argument("--clean-ratio", type=float, default=DEFAULT_CLEAN_RATIO)
    p.add_argument("--high-t", type=float, default=DEFAULT_HIGH_T)
    p.add_argument(
        "--ebeam", type=float, default=DEFAULT_EBEAM,
        help="Fallback beam energy only if the validated loader does not provide one",
    )
    p.add_argument("--force-km15", action="store_true")
    p.add_argument(
        "--phi-example-cells", type=int, default=6,
        help="Number of highest-|t| clean cells for data/KM15 phi-shape plots",
    )
    return p
#enddef


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    args.here = Path(__file__).resolve().parent
    args.external_scripts = args.here.parent.parent / "dvcs_cross_section" / "external_scripts"

    outdir = Path(args.outdir).expanduser()
    if not outdir.is_absolute():
        outdir = (args.here / outdir).resolve()
    #endif
    tables = outdir / "tables"
    figures = outdir / "figures"
    cache_dir = outdir / "cache"
    for directory in [tables, figures, cache_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    #endfor

    emff_path = resolve_existing(
        Path(args.emff_script),
        [args.external_scripts / "extract_emff_from_dvcs_bh.py"],
    )
    print(f"[backend] {emff_path}")
    emff = load_python_module(emff_path, "emff_high_luminosity_stage2_backend")

    data = load_clas12(args, emff)

    km15_cache = cache_dir / "clas12_km15_exact_kinematics.csv"
    data = evaluate_km15(
        data, emff, km15_cache, force=bool(args.force_km15)
    )

    data["preferred_ratio_region"] = data["t_over_Q2"] < float(args.clean_ratio)
    clean = data[data["preferred_ratio_region"]].copy().reset_index(drop=True)

    if len(clean) == 0:
        raise RuntimeError("No points survive the requested |t|/Q2 cut")
    #endif

    closure = make_closure_summary(
        data, float(args.clean_ratio), float(args.high_t)
    )
    pseudo = build_pseudodata(clean)
    cells = make_cell_summary(pseudo)
    t_summary = make_t_luminosity_summary(pseudo)

    # Core machine-readable products for the next fitting stage.
    data.to_csv(tables / "clas12_released_with_km15.csv", index=False)
    clean.to_csv(tables / "clas12_clean_with_km15.csv", index=False)
    closure.to_csv(tables / "km15_closure_summary.csv", index=False)
    pseudo.to_csv(tables / "clas12_km15_pseudodata_1x_2x_5x_10x.csv", index=False)
    cells.to_csv(tables / "pseudodata_cell_summary.csv", index=False)
    t_summary.to_csv(tables / "pseudodata_t_luminosity_summary.csv", index=False)

    # Also write one compact file per luminosity factor. These are convenient
    # direct inputs for the forthcoming CFF-fit driver.
    fit_columns = [
        "point_id", "bin", "xB", "Q2", "t_abs", "phi_deg", "ebeam",
        "t_over_Q2", "xs_pseudo", "stat_pseudo_abs",
        "ptp_sys_pseudo_abs", "point_unc_pseudo_abs",
        "correlated_norm_frac", "km15_ep", "km15_bh",
        "bh_fraction_of_km15",
    ]
    for factor in LUMI_FACTORS:
        d = pseudo[pseudo["luminosity_factor"] == factor]
        d[fit_columns].to_csv(
            tables / f"cff_fit_input_km15_{factor}x.csv", index=False
        )
    #endfor

    plot_closure_vs_kinematics(clean, figures)
    plot_high_t_closure(clean, float(args.high_t), figures)
    plot_closure_t_summary(closure, figures)
    plot_pseudo_precision(t_summary, figures)
    plot_selected_phi_cells(
        clean, figures, n_cells=max(0, int(args.phi_example_cells))
    )

    print_summary(data, clean, closure, float(args.clean_ratio))
    print(f"\n[output] {outdir}")
    print(f"[next] inspect {tables / 'km15_closure_summary.csv'}")
    print("[next] if closure is acceptable, use cff_fit_input_km15_{1,2,5,10}x.csv")
    print("       as the inputs to the Stage-2 CFF sensitivity fit.")
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
