#!/usr/bin/env python3
"""
Version 2: extract proton electromagnetic form-factor information from BH-dominated
CLAS12 ep -> e'p'gamma cross sections.

This script is designed to be run from
  /u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/external_scripts/

Default input:
  ../output/csvs/dvcs_pass2_analysis.csv

Main workflow
-------------
1. Read the preliminary CLAS12 10.6-GeV unpolarized cross section.
2. Evaluate KM15 at every measured bin center and decompose the prediction into
   BH, DVCS^2, and interference contributions using Gepard.
3. Build the five nested BH-dominated samples used in arXiv:2512.06554:
      |1 - sigma_BH/sigma_EP| <= 1%, 2%, 3%, 4%, 5%.
4. Precompute the exact quadratic BH dependence
      sigma_BH = A F1^2 + B F1 F2 + C F2^2
   at every data point. This makes fitting very fast; no KM15/CFF calculation
   is done inside the minimizer.
5. Reproduce the paper-style F1/F2 and GE/GM fits (dipole, P-pole, fixed-Kelly).
6. Perform direct Sachs-form-factor fits using the fitting families from
   Hayward et al., Phys. Rev. C 98, 045204 (2018), arXiv:1804.09150:
      P1-P4, IP1-IP4, CF2-CF4.
   CF4 is included explicitly as the preferred flexible fit.
7. Extract electric and magnetic rms radii from the slopes at Q^2=|t|=0.
8. Write fit tables, selected-point tables, and diagnostic/reproduction plots.

Physics conventions
-------------------
* q = |t| = -t > 0 is used as the form-factor argument.
* GE(0)=1 and GM(0)=mu_p for the direct Sachs fits.
* F1 = (GE + tau GM)/(1+tau), F2 = (GM-GE)/(1+tau),
  tau = q/(4 M_p^2).
* r_E^2 = -6/GE(0) dGE/dq |_{q=0};
  r_M^2 = -6/GM(0) dGM/dq |_{q=0}.
* The CSV phi is converted to the Trento angle exactly as in the existing
  km15_cli.py: phi_trento = pi - phi_csv.

Uncertainties
-------------
By default the pointwise fit uncertainty is
  sqrt(stat^2 + point_to_point_total^2 + (BH_cut * sigma)^2),
where the last term follows arXiv:2512.06554. The CSV "total scale sys" is
handled as one correlated multiplicative normalization nuisance N with a
Gaussian prior N=1 +/- median(total_scale_sys). The run-period-combination
systematic is retained as a diagnostic but is not double counted separately.

Dependencies
------------
Required: numpy, pandas, scipy, matplotlib, gepard

Typical usage
-------------
  python extract_emff_from_dvcs_bh_v1.py

  python extract_emff_from_dvcs_bh_v1.py --workers 6

  python extract_emff_from_dvcs_bh_v1.py --bh-cuts 0.01 0.02 0.03 0.04 0.05 \
      --fit-bh-cut 0.05 --sachs-families all

Outputs are written under:
  output/emff_from_bh/
"""

from __future__ import annotations

import argparse
import ast
import copy
import json
import math
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
MP = 0.9382720813          # GeV
MP2 = MP * MP
MU_P = 2.79284734463       # proton magnetic moment
KAPPA_P = MU_P - 1.0
HBARC = 0.1973269804       # GeV fm

DEFAULT_CSV = "../output/csvs/dvcs_pass2_analysis.csv"
DEFAULT_OUTDIR = "output/emff_from_bh"

COL_XB = "xBavg, 10.6 GeV"
COL_Q2 = "Q2avg, 10.6 GeV"
COL_T = "t_abs_avg, 10.6 GeV"
COL_PHI = "phiavg, 10.6 GeV"
COL_XS = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol"
COL_PTP_SYS = "Syst. err (point-to-point total)"
COL_COMB_SYS = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol, combination sys"
COL_SCALE_SYS = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol, total scale sys"


# -----------------------------------------------------------------------------
# Data containers
# -----------------------------------------------------------------------------
@dataclass
class FitResult:
    name: str
    category: str
    bh_cut: float
    npts: int
    npar: int
    chi2: float
    ndof: int
    chi2_ndof: float
    success: bool
    message: str
    params: np.ndarray
    param_names: List[str]
    covariance: np.ndarray
    rE_fm: float
    rE_err_fm: float
    rM_fm: float
    rM_err_fm: float
    norm: float
    norm_err: float
    model_kind: str
    meta: Dict[str, object]


# -----------------------------------------------------------------------------
# Parsing / I/O
# -----------------------------------------------------------------------------
def parse_tuple_cell(value) -> Tuple[float, float, float]:
    """Parse CSV cells such as '(2.01, 0.10, 0)' robustly."""
    if pd.isna(value):
        return np.nan, np.nan, np.nan
    #endif
    if isinstance(value, (tuple, list, np.ndarray)):
        vals = list(value)
    else:
        text = str(value).strip()
        try:
            vals = list(ast.literal_eval(text))
        except Exception:
            text = text.strip("()[]")
            vals = [float(x.strip()) for x in text.split(",") if x.strip()]
        #endtry
    #endif
    vals = vals + [np.nan] * (3 - len(vals))
    return float(vals[0]), float(vals[1]), float(vals[2])
#enddef


def load_clas12_csv(path: Path) -> pd.DataFrame:
    """Load only the columns needed by this analysis and unpack the XS tuple."""
    header = pd.read_csv(path, nrows=0).columns
    required = [COL_XB, COL_Q2, COL_T, COL_PHI, COL_XS]
    optional = [COL_PTP_SYS, COL_COMB_SYS, COL_SCALE_SYS, "bin index", "Bin Name",
                "xBmin", "xBmax", "Q2min", "Q2max",
                "t_abs_min", "t_abs_max", "phimin", "phimax"]
    missing = [c for c in required if c not in header]
    if missing:
        raise KeyError(f"Missing required CSV columns: {missing}")
    #endif

    usecols = required + [c for c in optional if c in header]
    df = pd.read_csv(path, usecols=usecols, low_memory=False)

    parsed = df[COL_XS].apply(parse_tuple_cell)
    df["xs"] = parsed.str[0]
    df["xs_stat"] = parsed.str[1]
    df["xs_aux"] = parsed.str[2]
    df["xB"] = pd.to_numeric(df[COL_XB], errors="coerce")
    df["Q2"] = pd.to_numeric(df[COL_Q2], errors="coerce")
    df["t_abs"] = pd.to_numeric(df[COL_T], errors="coerce")
    df["phi_deg"] = pd.to_numeric(df[COL_PHI], errors="coerce")

    if COL_PTP_SYS in df:
        # Absolute point-to-point systematic uncertainty, in the same units as
        # the measured cross section.
        df["ptp_sys_abs"] = pd.to_numeric(df[COL_PTP_SYS], errors="coerce").fillna(0.0)
    else:
        df["ptp_sys_abs"] = 0.0
    #endif

    if COL_COMB_SYS in df:
        # Fractional run-period-combination uncertainty. Retained for
        # diagnostics; it is already folded into the total-scale column.
        df["comb_sys_frac"] = pd.to_numeric(df[COL_COMB_SYS], errors="coerce").fillna(0.0)
    else:
        df["comb_sys_frac"] = 0.0
    #endif

    if COL_SCALE_SYS in df:
        # Fractional correlated scale uncertainty.
        df["scale_sys_frac"] = pd.to_numeric(df[COL_SCALE_SYS], errors="coerce")
    else:
        df["scale_sys_frac"] = np.nan
    #endif

    good = (
        np.isfinite(df["xB"]) & np.isfinite(df["Q2"]) &
        np.isfinite(df["t_abs"]) & np.isfinite(df["phi_deg"]) &
        np.isfinite(df["xs"]) & np.isfinite(df["xs_stat"]) &
        (df["xB"] > 0) & (df["xB"] < 1) &
        (df["Q2"] > 0) & (df["t_abs"] > 0) &
        (df["xs"] > 0) & (df["xs_stat"] > 0)
    )
    df = df.loc[good].copy().reset_index(drop=True)
    return df
#enddef


# -----------------------------------------------------------------------------
# Gepard / KM15 calculations
# -----------------------------------------------------------------------------
def make_gepard_point(g, xB: float, Q2: float, t_abs: float,
                      phi_deg: float, ebeam: float):
    phi_rad = math.radians(phi_deg)
    phi_trento = math.pi - phi_rad
    pt = g.DataPoint(
        xB=float(xB),
        t=-abs(float(t_abs)),
        Q2=float(Q2),
        phi=float(phi_trento),
        observable="XS",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(ebeam),
        in1charge=-1,
        in1polarization=0,
        in2particle="p",
    )
    pt.prepare()
    return pt
#enddef


def _set_eff_basis(model, f1_value: float, f2_value: float):
    """
    Temporarily replace F1/F2 on the model instance.

    Instance-assigned callables are not descriptors, so the lambda receives
    only the DataPoint argument, which matches self.m.F1(pt) usage in Gepard.
    """
    old_f1 = model.F1
    old_f2 = model.F2
    model.F1 = lambda pt: float(f1_value)
    model.F2 = lambda pt: float(f2_value)
    return old_f1, old_f2
#enddef


def _restore_eff(model, old_f1, old_f2) -> None:
    model.F1 = old_f1
    model.F2 = old_f2
#enddef


def evaluate_km15_point(args: Tuple[int, float, float, float, float, float]) -> Dict[str, float]:
    """
    Worker-safe KM15 evaluation for one point.

    Returns full EP, BH, DVCS^2, INT, R_BH, and the quadratic BH coefficients
    A/B/C such that sigma_BH=A*F1^2+B*F1*F2+C*F2^2.
    """
    idx, xB, Q2, t_abs, phi_deg, ebeam = args
    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError(
            "Could not import gepard/KM15. Activate the Python environment used "
            "by the existing km15_cli.py before running this script."
        ) from exc
    #endtry

    # Each worker process evaluates tasks sequentially. Reuse the imported KM15
    # object rather than deep-copying the full CFF/GPD model for every bin.
    # F1/F2 are restored immediately after the three BH basis evaluations.
    th = th_KM15
    pt = make_gepard_point(g, xB, Q2, t_abs, phi_deg, ebeam)

    pref = float(th.PreFacSigma(pt))
    bh_amp = float(th.TBH2unp(pt))
    int_amp = float(th.TINTunp(pt))
    dvcs_amp = float(th.TDVCS2unp(pt))

    sigma_bh = pref * bh_amp
    sigma_int = pref * int_amp
    sigma_dvcs = pref * dvcs_amp
    sigma_ep = sigma_bh + sigma_int + sigma_dvcs

    sigma_predict = float(th.predict(pt))
    ep_relerr = abs(sigma_ep - sigma_predict) / max(abs(sigma_predict), 1e-30)

    # Build exact quadratic coefficients for the BH-only fit.
    # Gepard's BH term is quadratic in F1,F2, so three basis evaluations suffice.
    model = th.m
    try:
        old_f1, old_f2 = _set_eff_basis(model, 1.0, 0.0)
        A = pref * float(th.TBH2unp(pt))
        model.F1 = lambda p: 0.0
        model.F2 = lambda p: 1.0
        C = pref * float(th.TBH2unp(pt))
        model.F1 = lambda p: 1.0
        model.F2 = lambda p: 1.0
        one_one = pref * float(th.TBH2unp(pt))
        B = one_one - A - C
    finally:
        _restore_eff(model, old_f1, old_f2)
    #endtry

    # Validate the quadratic representation using the model's nominal FFs.
    f1_nom = float(model.F1(pt))
    f2_nom = float(model.F2(pt))
    sigma_bh_reco = A * f1_nom**2 + B * f1_nom * f2_nom + C * f2_nom**2
    denom = max(abs(sigma_bh), 1e-30)
    quad_relerr = abs(sigma_bh_reco - sigma_bh) / denom

    rbh = sigma_bh / sigma_ep if sigma_ep != 0 else np.nan
    return {
        "_row": idx,
        "km15_ep": sigma_ep,
        "km15_bh": sigma_bh,
        "km15_dvcs": sigma_dvcs,
        "km15_int": sigma_int,
        "R_BH": rbh,
        "bh_delta": abs(1.0 - rbh) if np.isfinite(rbh) else np.nan,
        "bh_A": A,
        "bh_B": B,
        "bh_C": C,
        "bh_quad_relerr": quad_relerr,
        "km15_F1": f1_nom,
        "km15_F2": f2_nom,
        "km15_ep_predict": sigma_predict,
        "km15_ep_decomp_relerr": ep_relerr,
    }
#enddef


def evaluate_km15_dataframe(df: pd.DataFrame, ebeam: float, workers: int,
                            cache_path: Path, force: bool = False) -> pd.DataFrame:
    """Evaluate or load cached KM15/BH decomposition."""
    cache_cols = ["_row", "km15_ep", "km15_bh", "km15_dvcs", "km15_int",
                  "R_BH", "bh_delta", "bh_A", "bh_B", "bh_C",
                  "bh_quad_relerr", "km15_F1", "km15_F2",
                  "km15_ep_predict", "km15_ep_decomp_relerr"]
    if cache_path.exists() and not force:
        cache = pd.read_csv(cache_path)
        if len(cache) == len(df) and all(c in cache.columns for c in cache_cols):
            out = df.copy()
            for c in cache_cols:
                if c != "_row":
                    out[c] = cache[c].to_numpy()
                #endif
            #endfor
            print(f"[KM15] Loaded cache: {cache_path}")
            return out
        #endif
    #endif

    tasks = [
        (int(i), float(r.xB), float(r.Q2), float(r.t_abs), float(r.phi_deg), float(ebeam))
        for i, r in df[["xB", "Q2", "t_abs", "phi_deg"]].iterrows()
    ]

    print(f"[KM15] Evaluating {len(tasks)} CLAS12 points with {workers} worker(s)...")
    t0 = time.time()
    if workers <= 1:
        rows = []
        for j, task in enumerate(tasks, start=1):
            rows.append(evaluate_km15_point(task))
            if j % 100 == 0 or j == len(tasks):
                print(f"[KM15] {j:5d}/{len(tasks)}  elapsed={time.time()-t0:8.1f} s")
            #endif
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            rows = list(ex.map(evaluate_km15_point, tasks, chunksize=8))
        #endwith
    #endif

    km = pd.DataFrame(rows).sort_values("_row").reset_index(drop=True)
    if np.nanmax(km["bh_quad_relerr"].to_numpy()) > 1e-8:
        raise RuntimeError(
            "The quadratic BH coefficient reconstruction failed validation. "
            f"max relative error={np.nanmax(km['bh_quad_relerr']):.3e}."
        )
    #endif
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    km.to_csv(cache_path, index=False)
    print(f"[KM15] Finished in {time.time()-t0:.1f} s; cache -> {cache_path}")

    out = df.copy()
    for c in cache_cols:
        if c != "_row":
            out[c] = km[c].to_numpy()
        #endif
    #endfor
    return out
#enddef


# -----------------------------------------------------------------------------
# Elastic FF utilities
# -----------------------------------------------------------------------------
def kelly_f1_f2(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Kelly (2004) proton elastic form factors, evaluated at q=|t|.

    G_E and G_M are written in tau=q/(4 M_p^2), then converted to F1/F2.
    """
    q = np.asarray(q, dtype=float)
    tau = q / (4.0 * MP2)

    ge = (
        1.0 - 0.24 * tau
    ) / (
        1.0 + 10.98 * tau + 12.82 * tau**2 + 21.97 * tau**3
    )

    gm = MU_P * (
        1.0 + 0.12 * tau
    ) / (
        1.0 + 10.97 * tau + 18.86 * tau**2 + 6.55 * tau**3
    )

    return sachs_to_f1f2(q, ge, gm)
#enddef


def f1f2_to_sachs(q: np.ndarray, f1: np.ndarray, f2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    tau = np.asarray(q, dtype=float) / (4.0 * MP2)
    ge = f1 - tau * f2
    gm = f1 + f2
    return ge, gm
#enddef


def sachs_to_f1f2(q: np.ndarray, ge: np.ndarray, gm: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    tau = np.asarray(q, dtype=float) / (4.0 * MP2)
    f1 = (ge + tau * gm) / (1.0 + tau)
    f2 = (gm - ge) / (1.0 + tau)
    return f1, f2
#enddef


def bh_from_f1f2(A: np.ndarray, B: np.ndarray, C: np.ndarray,
                  f1: np.ndarray, f2: np.ndarray) -> np.ndarray:
    return A * f1**2 + B * f1 * f2 + C * f2**2
#enddef


# -----------------------------------------------------------------------------
# 2018 Sachs fitting families
# -----------------------------------------------------------------------------
def polynomial_shape(q: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    out = np.ones_like(q)
    for power, coeff in enumerate(coeffs, start=1):
        out = out + coeff * q**power
    #endfor
    return out
#enddef


def inverse_polynomial_shape(q: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    den = np.ones_like(q)
    for power, coeff in enumerate(coeffs, start=1):
        den = den + coeff * q**power
    #endfor
    return 1.0 / den
#enddef


def continued_fraction_shape(q: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    """
    Continued fraction used in arXiv:1804.09150, normalized to one at q=0:
      1 / [1 + c1 q/(1 + c2 q/(1 + ...))].
    """
    q = np.asarray(q, dtype=float)
    if len(coeffs) < 2:
        raise ValueError("CF family requires at least CF2.")
    #endif
    den = np.ones_like(q)
    for coeff in coeffs[:0:-1]:
        den = 1.0 + coeff * q / den
    #endfor
    return 1.0 / (1.0 + coeffs[0] * q / den)
#enddef


def sachs_shape(family: str, q: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    family = family.upper()
    if family.startswith("IP"):
        return inverse_polynomial_shape(q, coeffs)
    #endif
    if family.startswith("CF"):
        return continued_fraction_shape(q, coeffs)
    #endif
    if family.startswith("P"):
        return polynomial_shape(q, coeffs)
    #endif
    raise ValueError(f"Unknown Sachs family: {family}")
#enddef


def family_order(family: str) -> int:
    """Return the number of free shape coefficients for a Sachs fit family."""
    family = family.upper().strip()

    if family.startswith("IP"):
        suffix = family[2:]
    elif family.startswith("CF"):
        suffix = family[2:]
    elif family.startswith("P"):
        suffix = family[1:]
    else:
        raise ValueError(f"Unknown Sachs family: {family}")
    #endif

    if not suffix.isdigit():
        raise ValueError(f"Could not parse order from Sachs family: {family}")
    #endif

    n = int(suffix)
    if n < 1:
        raise ValueError(f"Sachs family order must be >= 1: {family}")
    #endif
    return n
#enddef


def family_initial(family: str, electric: bool) -> np.ndarray:
    """Reasonable low-q starting values; minimizer is not sensitive to exact choice."""
    n = family_order(family)
    family = family.upper()
    r0 = 0.84 if electric else 0.85
    slope = (r0 / HBARC)**2 / 6.0  # positive magnitude in GeV^-2
    if family.startswith("P") and not family.startswith("IP"):
        arr = np.zeros(n)
        arr[0] = -slope
    else:
        arr = np.zeros(n)
        arr[0] = slope
    #endif
    return arr
#enddef


# -----------------------------------------------------------------------------
# Generic least-squares machinery
# -----------------------------------------------------------------------------
def covariance_from_least_squares(res, ndata_residuals: int) -> np.ndarray:
    """Hessian/Jacobian covariance, with Delta-chi2=1 convention."""
    J = np.asarray(res.jac, dtype=float)
    try:
        jtj_inv = np.linalg.pinv(J.T @ J, rcond=1e-12)
    except Exception:
        jtj_inv = np.full((J.shape[1], J.shape[1]), np.nan)
    #endtry
    return jtj_inv
#enddef


def numerical_gradient(fun: Callable[[np.ndarray], float], p: np.ndarray) -> np.ndarray:
    grad = np.zeros_like(p, dtype=float)
    for i in range(len(p)):
        step = 1e-5 * max(1.0, abs(float(p[i])))
        pp = p.copy(); pm = p.copy()
        pp[i] += step; pm[i] -= step
        fp = fun(pp); fm = fun(pm)
        grad[i] = (fp - fm) / (2.0 * step)
    #endfor
    return grad
#enddef


def propagate_scalar(fun: Callable[[np.ndarray], float], p: np.ndarray,
                     cov: np.ndarray) -> Tuple[float, float]:
    value = float(fun(p))
    if not np.isfinite(value) or cov.size == 0 or not np.all(np.isfinite(cov)):
        return value, np.nan
    #endif
    grad = numerical_gradient(fun, p)
    var = float(grad @ cov @ grad)
    err = math.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan
    return value, err
#enddef


def radius_from_shape(shape_fun: Callable[[np.ndarray], np.ndarray], g0: float) -> float:
    """Numerically extract rms radius from normalized or unnormalized FF shape."""
    h = 1e-6  # GeV^2
    g_at_0 = float(np.asarray(shape_fun(np.array([0.0])))[0])
    g_h = float(np.asarray(shape_fun(np.array([h])))[0])
    deriv = (g_h - g_at_0) / h
    r2_gev = -6.0 * deriv / g_at_0
    if not np.isfinite(r2_gev) or r2_gev <= 0:
        return np.nan
    #endif
    return math.sqrt(r2_gev) * HBARC
#enddef


def fit_sigma_errors(data: pd.DataFrame, bh_cut: float,
                     include_ptp_sys: bool, include_bh_sys: bool) -> np.ndarray:
    """
    Construct the uncorrelated pointwise uncertainty.

    Default:
      sqrt(stat^2 + ptp_total^2 + (BH_cut * sigma)^2)

    The pass-2 total-scale uncertainty is handled separately as a single
    correlated normalization nuisance.
    """
    stat = data["xs_stat"].to_numpy(float)
    var = stat**2

    if include_ptp_sys:
        ptp = data["ptp_sys_abs"].to_numpy(float)
        var = var + ptp**2
    #endif

    if include_bh_sys:
        var = var + (bh_cut * data["xs"].to_numpy(float))**2
    #endif

    err = np.sqrt(var)
    bad = ~np.isfinite(err) | (err <= 0.0)
    if np.any(bad):
        err[bad] = stat[bad]
    #endif
    return err
#enddef


def scale_prior_sigma(data: pd.DataFrame, fallback: float) -> float:
    vals = data["scale_sys_frac"].to_numpy(float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if len(vals) == 0:
        return fallback
    #endif
    return float(np.median(vals))
#enddef


# -----------------------------------------------------------------------------
# Paper-style F1/F2 fits
# -----------------------------------------------------------------------------
def paper_model_f1f2(kind: str, q: np.ndarray, pars: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    kind = kind.lower()
    if kind == "dipole":
        aE, aM = pars[:2]
        f1 = (1.0 + q / aE)**-2
        f2 = KAPPA_P * (1.0 + q / aM)**-2
        return f1, f2
    #endif
    if kind == "ppole_same_a":
        a, pE, pM = pars[:3]
        f1 = (1.0 + q / a)**(-pE)
        f2 = KAPPA_P * (1.0 + q / a)**(-pM)
        return f1, f2
    #endif
    if kind == "ppole_same_p":
        aE, aM, p = pars[:3]
        f1 = (1.0 + q / aE)**(-p)
        f2 = KAPPA_P * (1.0 + q / aM)**(-p)
        return f1, f2
    #endif
    if kind == "f2_kelly_ppole_f1":
        aE, pE = pars[:2]
        f1 = (1.0 + q / aE)**(-pE)
        _, f2 = kelly_f1_f2(q)
        return f1, f2
    #endif
    raise ValueError(kind)
#enddef


def paper_model_setup(kind: str) -> Tuple[List[str], np.ndarray, np.ndarray, np.ndarray]:
    kind = kind.lower()
    if kind == "dipole":
        return ["aE", "aM"], np.array([0.95, 0.56]), np.array([0.05, 0.05]), np.array([5.0, 5.0])
    #endif
    if kind == "ppole_same_a":
        return ["a", "pE", "pM"], np.array([0.53, 1.2, 1.9]), np.array([0.03, 0.1, 0.1]), np.array([5.0, 6.0, 6.0])
    #endif
    if kind == "ppole_same_p":
        return ["aE", "aM", "p"], np.array([0.65, 0.35, 1.4]), np.array([0.03, 0.03, 0.1]), np.array([5.0, 5.0, 6.0])
    #endif
    if kind == "f2_kelly_ppole_f1":
        return ["aE", "pE"], np.array([0.56, 1.3]), np.array([0.03, 0.1]), np.array([5.0, 6.0])
    #endif
    raise ValueError(kind)
#enddef


def fit_paper_model(data: pd.DataFrame, kind: str, bh_cut: float,
                    include_ptp_sys: bool, include_bh_sys: bool,
                    use_scale_nuisance: bool, scale_fallback: float) -> FitResult:
    q = data["t_abs"].to_numpy(float)
    y = data["xs"].to_numpy(float)
    err = fit_sigma_errors(data, bh_cut, include_ptp_sys, include_bh_sys)
    A = data["bh_A"].to_numpy(float)
    B = data["bh_B"].to_numpy(float)
    C = data["bh_C"].to_numpy(float)

    names, p0, lo, hi = paper_model_setup(kind)
    scale_sig = scale_prior_sigma(data, scale_fallback)
    if use_scale_nuisance:
        names = names + ["N"]
        p0 = np.r_[p0, 1.0]
        lo = np.r_[lo, 0.5]
        hi = np.r_[hi, 1.5]
    #endif

    def residuals(p):
        core = p[:-1] if use_scale_nuisance else p
        norm = p[-1] if use_scale_nuisance else 1.0
        f1, f2 = paper_model_f1f2(kind, q, core)
        pred = norm * bh_from_f1f2(A, B, C, f1, f2)
        res = (pred - y) / err
        if use_scale_nuisance:
            res = np.r_[res, (norm - 1.0) / scale_sig]
        #endif
        return res
    #enddef

    res = least_squares(residuals, p0, bounds=(lo, hi), method="trf",
                        x_scale="jac", max_nfev=20000)
    cov = covariance_from_least_squares(res, len(y))
    chi2 = float(np.sum(res.fun**2))
    npar = len(res.x)
    ndof = max(len(res.fun) - npar, 1)

    def re_of_p(p):
        core = p[:-1] if use_scale_nuisance else p
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(qq, *paper_model_f1f2(kind, qq, core))[0], 1.0
        )
    #enddef
    def rm_of_p(p):
        core = p[:-1] if use_scale_nuisance else p
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(qq, *paper_model_f1f2(kind, qq, core))[1], MU_P
        )
    #enddef

    rE, drE = propagate_scalar(re_of_p, res.x, cov)
    rM, drM = propagate_scalar(rm_of_p, res.x, cov)
    norm = float(res.x[-1]) if use_scale_nuisance else 1.0
    norm_err = math.sqrt(max(cov[-1, -1], 0.0)) if use_scale_nuisance and np.isfinite(cov[-1, -1]) else 0.0

    return FitResult(
        name=kind,
        category="paper_F1F2",
        bh_cut=bh_cut,
        npts=len(y), npar=npar, chi2=chi2, ndof=ndof,
        chi2_ndof=chi2/ndof, success=bool(res.success), message=str(res.message),
        params=res.x.copy(), param_names=names, covariance=cov,
        rE_fm=rE, rE_err_fm=drE, rM_fm=rM, rM_err_fm=drM,
        norm=norm, norm_err=norm_err, model_kind=kind,
        meta={"scale_prior_sigma": scale_sig},
    )
#enddef


# -----------------------------------------------------------------------------
# Direct Sachs fits (YAHL18-style families)
# -----------------------------------------------------------------------------
def fit_sachs_family(data: pd.DataFrame, family: str, bh_cut: float,
                     include_ptp_sys: bool, include_bh_sys: bool,
                     use_scale_nuisance: bool, scale_fallback: float,
                     fix_gm_to_kelly: bool = False) -> FitResult:
    family = family.upper()
    n = family_order(family)
    q = data["t_abs"].to_numpy(float)
    y = data["xs"].to_numpy(float)
    err = fit_sigma_errors(data, bh_cut, include_ptp_sys, include_bh_sys)
    A = data["bh_A"].to_numpy(float)
    B = data["bh_B"].to_numpy(float)
    C = data["bh_C"].to_numpy(float)

    pE0 = family_initial(family, electric=True)
    names = [f"GE_c{i}" for i in range(1, n+1)]
    p0 = pE0.copy()

    if fix_gm_to_kelly:
        gm_n = 0
    else:
        pM0 = family_initial(family, electric=False)
        p0 = np.r_[p0, pM0]
        names += [f"GM_c{i}" for i in range(1, n+1)]
        gm_n = n
    #endif

    scale_sig = scale_prior_sigma(data, scale_fallback)
    if use_scale_nuisance:
        p0 = np.r_[p0, 1.0]
        names += ["N"]
    #endif

    # Broad finite bounds protect continued fractions from pathological poles
    # while remaining much wider than physically useful low-|t| coefficients.
    lo = np.full(len(p0), -100.0)
    hi = np.full(len(p0), +100.0)
    if use_scale_nuisance:
        lo[-1], hi[-1] = 0.5, 1.5
    #endif

    def unpack(p, qq):
        pcore = p[:-1] if use_scale_nuisance else p
        ce = pcore[:n]
        ge = sachs_shape(family, qq, ce)
        if fix_gm_to_kelly:
            kf1, kf2 = kelly_f1_f2(qq)
            _, gm = f1f2_to_sachs(qq, kf1, kf2)
        else:
            cm = pcore[n:n+gm_n]
            gm = MU_P * sachs_shape(family, qq, cm)
        #endif
        return ge, gm
    #enddef

    def residuals(p):
        norm = p[-1] if use_scale_nuisance else 1.0
        ge, gm = unpack(p, q)
        f1, f2 = sachs_to_f1f2(q, ge, gm)
        pred = norm * bh_from_f1f2(A, B, C, f1, f2)
        res = (pred - y) / err
        if use_scale_nuisance:
            res = np.r_[res, (norm - 1.0) / scale_sig]
        #endif
        if not np.all(np.isfinite(res)):
            res = np.nan_to_num(res, nan=1e12, posinf=1e12, neginf=-1e12)
        #endif
        return res
    #enddef

    res = least_squares(residuals, p0, bounds=(lo, hi), method="trf",
                        x_scale="jac", max_nfev=30000)
    cov = covariance_from_least_squares(res, len(y))
    chi2 = float(np.sum(res.fun**2))
    npar = len(res.x)
    ndof = max(len(res.fun) - npar, 1)

    def re_of_p(p):
        return radius_from_shape(lambda qq: unpack(p, qq)[0], 1.0)
    #enddef
    def rm_of_p(p):
        return radius_from_shape(lambda qq: unpack(p, qq)[1], MU_P)
    #enddef

    rE, drE = propagate_scalar(re_of_p, res.x, cov)
    rM, drM = propagate_scalar(rm_of_p, res.x, cov)
    norm = float(res.x[-1]) if use_scale_nuisance else 1.0
    norm_err = math.sqrt(max(cov[-1, -1], 0.0)) if use_scale_nuisance and np.isfinite(cov[-1, -1]) else 0.0
    suffix = "_GMKelly" if fix_gm_to_kelly else ""

    return FitResult(
        name=f"{family}{suffix}", category="Sachs_2018", bh_cut=bh_cut,
        npts=len(y), npar=npar, chi2=chi2, ndof=ndof, chi2_ndof=chi2/ndof,
        success=bool(res.success), message=str(res.message),
        params=res.x.copy(), param_names=names, covariance=cov,
        rE_fm=rE, rE_err_fm=drE, rM_fm=rM, rM_err_fm=drM,
        norm=norm, norm_err=norm_err, model_kind=family,
        meta={"fix_gm_to_kelly": fix_gm_to_kelly, "scale_prior_sigma": scale_sig},
    )
#enddef


def evaluate_fit_form_factors(fr: FitResult, q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return F1,F2,GE,GM for plotting for any FitResult."""
    q = np.asarray(q, dtype=float)
    use_n = "N" in fr.param_names
    core = fr.params[:-1] if use_n else fr.params

    if fr.category == "paper_F1F2":
        f1, f2 = paper_model_f1f2(fr.model_kind, q, core)
        ge, gm = f1f2_to_sachs(q, f1, f2)
        return f1, f2, ge, gm
    #endif

    family = fr.model_kind.upper()
    n = family_order(family)
    ge = sachs_shape(family, q, core[:n])
    if bool(fr.meta.get("fix_gm_to_kelly", False)):
        f1k, f2k = kelly_f1_f2(q)
        _, gm = f1f2_to_sachs(q, f1k, f2k)
    else:
        gm = MU_P * sachs_shape(family, q, core[n:2*n])
    #endif
    f1, f2 = sachs_to_f1f2(q, ge, gm)
    return f1, f2, ge, gm
#enddef


# -----------------------------------------------------------------------------
# Fit bookkeeping / uncertainties bands
# -----------------------------------------------------------------------------
def fitresult_to_record(fr: FitResult) -> Dict[str, object]:
    rec = {
        "name": fr.name, "category": fr.category, "bh_cut": fr.bh_cut,
        "npts": fr.npts, "npar": fr.npar, "chi2": fr.chi2,
        "ndof": fr.ndof, "chi2_ndof": fr.chi2_ndof,
        "success": fr.success, "message": fr.message,
        "rE_fm": fr.rE_fm, "rE_err_fm": fr.rE_err_fm,
        "rM_fm": fr.rM_fm, "rM_err_fm": fr.rM_err_fm,
        "norm": fr.norm, "norm_err": fr.norm_err,
    }
    for name, value in zip(fr.param_names, fr.params):
        rec[name] = value
    #endfor
    return rec
#enddef


def form_factor_band(fr: FitResult, q: np.ndarray, which: str) -> Tuple[np.ndarray, np.ndarray]:
    idx = {"F1": 0, "F2": 1, "GE": 2, "GM": 3}[which]
    central = evaluate_fit_form_factors(fr, q)[idx]
    if not np.all(np.isfinite(fr.covariance)):
        return central, np.full_like(central, np.nan)
    #endif

    sig = np.zeros_like(q)
    for j, qj in enumerate(q):
        def scalar(p):
            tmp = FitResult(**{**fr.__dict__, "params": p})
            return float(evaluate_fit_form_factors(tmp, np.array([qj]))[idx][0])
        #enddef
        grad = numerical_gradient(scalar, fr.params)
        var = float(grad @ fr.covariance @ grad)
        sig[j] = math.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan
    #endfor
    return central, sig
#enddef


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def save_bh_selection_plots(df: pd.DataFrame, bh_cuts: Sequence[float], outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    vals = df["R_BH"].to_numpy(float)
    vals = vals[np.isfinite(vals)]
    ax.hist(vals, bins=80, histtype="step", linewidth=1.5)
    for cut in bh_cuts:
        ax.axvline(1.0-cut, linestyle="--", linewidth=1.0, label=f"{100*cut:.0f}% lower bound")
        ax.axvline(1.0+cut, linestyle="--", linewidth=1.0)
    #endfor
    ax.set_xlabel(r"$R_{BH}=\sigma_{BH}/\sigma_{EP}^{KM15}$")
    ax.set_ylabel("CLAS12 bins")
    ax.set_title("KM15 BH dominance for preliminary CLAS12 bins")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "01_bh_fraction_distribution.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    sc = ax.scatter(df["phi_deg"], df["xB"], c=df["R_BH"], s=15,
                    vmin=0.8, vmax=1.05)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$R_{BH}$")
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$x_B$")
    ax.set_title("CLAS12 kinematics colored by KM15 BH fraction")
    fig.tight_layout()
    fig.savefig(outdir / "02_bh_fraction_phi_xB.pdf")
    plt.close(fig)
#enddef


def save_selected_cross_section_pages(data: pd.DataFrame, outpath: Path) -> None:
    """Data-vs-KM15/BH phi panels, grouped by the parent xB,Q2,t bin edges when available."""
    group_cols = [c for c in ["xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max"] if c in data.columns]
    if len(group_cols) < 6:
        # Robust fallback: grouping by rounded mean kinematics.
        work = data.copy()
        work["_gx"] = work["xB"].round(3)
        work["_gq"] = work["Q2"].round(2)
        work["_gt"] = work["t_abs"].round(2)
        group_cols = ["_gx", "_gq", "_gt"]
    else:
        work = data.copy()
    #endif

    groups = [(key, g.sort_values("phi_deg")) for key, g in work.groupby(group_cols, dropna=False)]
    with PdfPages(outpath) as pdf:
        for start in range(0, len(groups), 6):
            chunk = groups[start:start+6]
            fig, axes = plt.subplots(3, 2, figsize=(11, 11), squeeze=False)
            for ax, item in zip(axes.flat, chunk):
                key, g = item
                yerr = g["xs_stat"].to_numpy(float)
                ax.errorbar(g["phi_deg"], g["xs"], yerr=yerr, fmt="o", markersize=3.5,
                            capsize=2, label="CLAS12 preliminary")
                ax.plot(g["phi_deg"], g["km15_ep"], "-", linewidth=1.2, label="KM15 EP")
                ax.plot(g["phi_deg"], g["km15_bh"], "--", linewidth=1.2, label="BH only")
                ax.set_xlabel(r"$\phi$ (deg)")
                ax.set_ylabel(r"$d^4\sigma$")
                ax.set_title(
                    rf"$\langle x_B\rangle={g['xB'].mean():.3f}$, "
                    rf"$\langle Q^2\rangle={g['Q2'].mean():.2f}$, "
                    rf"$\langle|t|\rangle={g['t_abs'].mean():.3f}$"
                )
                ax.grid(alpha=0.2)
            #endfor
            for ax in axes.flat[len(chunk):]:
                ax.axis("off")
            #endfor
            handles, labels = axes.flat[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="upper center", ncol=3)
            fig.suptitle("BH-dominated CLAS12 cross sections", y=0.995)
            fig.tight_layout(rect=(0, 0, 1, 0.965))
            pdf.savefig(fig)
            plt.close(fig)
        #endfor
    #endwith
#enddef


def save_paper_ff_plots(results: List[FitResult], outdir: Path, preferred_cut: float) -> None:
    paper = [r for r in results if r.category == "paper_F1F2" and abs(r.bh_cut-preferred_cut) < 1e-12]
    if not paper:
        return
    #endif
    q = np.linspace(0.0, 1.0, 300)
    f1k, f2k = kelly_f1_f2(q)
    gek, gmk = f1f2_to_sachs(q, f1k, f2k)

    for which, ylabel, filename in [
        ("F1", r"$F_1(t)$", "04_F1_paper_style.pdf"),
        ("F2", r"$F_2(t)$", "05_F2_paper_style.pdf"),
    ]:
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        for fr in paper:
            vals = evaluate_fit_form_factors(fr, q)[0 if which == "F1" else 1]
            ax.plot(q, vals, linewidth=1.2, label=fr.name)
        #endfor
        ax.plot(q, f1k if which == "F1" else f2k, linewidth=1.7, label="Kelly")
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0, 1.0)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir / filename)
        plt.close(fig)
    #endfor

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for fr in paper:
        _, _, ge, gm = evaluate_fit_form_factors(fr, q)
        axes[0].plot(q, ge, linewidth=1.2, label=fr.name)
        axes[1].plot(q, gm, linewidth=1.2, label=fr.name)
    #endfor
    axes[0].plot(q, gek, linewidth=1.7, label="Kelly")
    axes[1].plot(q, gmk, linewidth=1.7, label="Kelly")
    axes[0].set_ylabel(r"$G_E(t)$")
    axes[1].set_ylabel(r"$G_M(t)$")
    for ax in axes:
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_xlim(0, 1.0)
        ax.grid(alpha=0.2)
    #endfor
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "06_GE_GM_paper_style.pdf")
    plt.close(fig)
#enddef


def save_sachs_family_plots(results: List[FitResult], outdir: Path, preferred_cut: float) -> None:
    fam = [r for r in results if r.category == "Sachs_2018" and abs(r.bh_cut-preferred_cut) < 1e-12
           and not bool(r.meta.get("fix_gm_to_kelly", False))]
    if not fam:
        return
    #endif
    q = np.linspace(0.0, 0.8, 300)
    f1k, f2k = kelly_f1_f2(q)
    gek, gmk = f1f2_to_sachs(q, f1k, f2k)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    for fr in fam:
        _, _, ge, gm = evaluate_fit_form_factors(fr, q)
        lw = 2.2 if fr.name == "CF4" else 1.0
        axes[0].plot(q, ge, linewidth=lw, label=fr.name)
        axes[1].plot(q, gm, linewidth=lw, label=fr.name)
    #endfor
    axes[0].plot(q, gek, "k--", linewidth=1.5, label="Kelly")
    axes[1].plot(q, gmk, "k--", linewidth=1.5, label="Kelly")
    axes[0].set_ylabel(r"$G_E$")
    axes[1].set_ylabel(r"$G_M$")
    for ax in axes:
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.grid(alpha=0.2)
    #endfor
    axes[0].legend(fontsize=7, ncol=2)
    fig.suptitle("Direct Sachs fits using the 2018 fitting families")
    fig.tight_layout()
    fig.savefig(outdir / "07_GE_GM_2018_families.pdf")
    plt.close(fig)

    cf4 = next((r for r in fam if r.name == "CF4"), None)
    if cf4 is not None:
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
        for ax, which, ylabel in zip(axes, ["GE", "GM"], [r"$G_E$", r"$G_M$"]):
            cen, sig = form_factor_band(cf4, q, which)
            ax.plot(q, cen, linewidth=1.8, label="CLAS12 CF4")
            ax.fill_between(q, cen-sig, cen+sig, alpha=0.25)
            ref = gek if which == "GE" else gmk
            ax.plot(q, ref, "k--", linewidth=1.5, label="Kelly")
            ax.set_xlabel(r"$|t|$ (GeV$^2$)")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.2)
        #endfor
        axes[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir / "08_CF4_GE_GM_with_band.pdf")
        plt.close(fig)
    #endif
#enddef


def save_radii_plot(results: List[FitResult], outdir: Path, preferred_cut: float) -> None:
    rows = [r for r in results if abs(r.bh_cut-preferred_cut) < 1e-12]
    rows = [r for r in rows if np.isfinite(r.rE_fm) or np.isfinite(r.rM_fm)]
    if not rows:
        return
    #endif
    labels = [r.name for r in rows]
    y = np.arange(len(rows))
    fig, axes = plt.subplots(1, 2, figsize=(12, max(5.0, 0.38*len(rows)+1.5)), sharey=True)
    axes[0].errorbar([r.rE_fm for r in rows], y,
                     xerr=[r.rE_err_fm for r in rows], fmt="o", capsize=2)
    axes[1].errorbar([r.rM_fm for r in rows], y,
                     xerr=[r.rM_err_fm for r in rows], fmt="o", capsize=2)
    axes[0].axvline(0.8409, linestyle="--", linewidth=1.2, label="PDG rE ~ 0.8409 fm")
    axes[1].axvline(0.851, linestyle="--", linewidth=1.2, label="reference rM ~ 0.851 fm")
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel(r"$r_E$ (fm)")
    axes[1].set_xlabel(r"$r_M$ (fm)")
    axes[0].grid(alpha=0.2); axes[1].grid(alpha=0.2)
    axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "09_extracted_radii.pdf")
    plt.close(fig)
#enddef


def save_bh_f1_f2_sensitivity(outdir: Path, ebeam: float) -> None:
    """Reproduce the qualitative Fig. 7 exercise at the paper's representative point."""
    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception:
        return
    #endtry
    q = np.linspace(0.01, 0.45, 180)
    full = np.zeros_like(q); f1zero = np.zeros_like(q); f2zero = np.zeros_like(q)
    for i, qi in enumerate(q):
        th = th_KM15
        # Paper uses phi=3.012 rad. Convert back to CSV-angle convention so
        # make_gepard_point applies phi_trento = pi - phi_csv.
        phi_csv_deg = math.degrees(math.pi - 3.012)
        pt = make_gepard_point(g, 0.1, 1.0, qi, phi_csv_deg, ebeam)
        pref = float(th.PreFacSigma(pt))
        model = th.m
        old_f1, old_f2 = model.F1, model.F2
        f1k, f2k = kelly_f1_f2(np.array([qi]))
        try:
            model.F1 = lambda p, v=float(f1k[0]): v
            model.F2 = lambda p, v=float(f2k[0]): v
            full[i] = pref * float(th.TBH2unp(pt))
            model.F1 = lambda p: 0.0
            f1zero[i] = pref * float(th.TBH2unp(pt))
            model.F1 = lambda p, v=float(f1k[0]): v
            model.F2 = lambda p: 0.0
            f2zero[i] = pref * float(th.TBH2unp(pt))
        finally:
            model.F1, model.F2 = old_f1, old_f2
        #endtry
    #endfor
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(q, full, label="Full BH")
    ax.plot(q, f1zero, "--", label=r"$F_1=0$")
    ax.plot(q, f2zero, ":", label=r"$F_2=0$")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$d^4\sigma$")
    ax.set_title(r"BH sensitivity at $x_B=0.1$, $Q^2=1$ GeV$^2$, $\phi=3.012$ rad")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "10_BH_F1_F2_sensitivity.pdf")
    plt.close(fig)
#enddef


# -----------------------------------------------------------------------------
# CLI / main
# -----------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--csv", default=DEFAULT_CSV, help="CLAS12 DVCS cross-section CSV")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR, help="Output directory")
    p.add_argument("--ebeam", type=float, default=10.604, help="Beam energy in GeV")
    p.add_argument("--workers", type=int, default=min(6, os.cpu_count() or 1), help="KM15 worker processes")
    p.add_argument("--bh-cuts", type=float, nargs="+", default=[0.01, 0.02, 0.03, 0.04, 0.05],
                   help="Nested |1-R_BH| selection thresholds")
    p.add_argument("--fit-bh-cut", type=float, default=0.05,
                   help="BH cut used for the flexible/Sachs fits")
    p.add_argument("--force-km15", action="store_true", help="Ignore the cached KM15 decomposition")
    p.add_argument("--no-point-to-point-systematics", action="store_true",
                   help="Omit the CSV absolute point-to-point systematic from the fit uncertainty")
    p.add_argument("--no-bh-selection-systematic", action="store_true",
                   help="Do not add the paper's 1-5% uncorrelated BH-selection uncertainty")
    p.add_argument("--no-scale-nuisance", action="store_true",
                   help="Do not use the correlated total-scale normalization nuisance")
    p.add_argument("--scale-fallback", type=float, default=0.077,
                   help="Fallback fractional normalization prior if total-scale column is absent")
    p.add_argument("--sachs-families", nargs="+", default=["all"],
                   help="Subset of P1-P4 IP1-IP4 CF2-CF4, or 'all'")
    p.add_argument("--fit-gm-kelly-variants", action="store_true",
                   help="Also fit each Sachs GE family with GM fixed to Kelly")
    p.add_argument("--tmin", type=float, default=None, help="Optional minimum |t| for fits")
    p.add_argument("--tmax", type=float, default=None, help="Optional maximum |t| for fits")
    return p
#enddef


def main() -> int:
    args = build_parser().parse_args()
    csv_path = Path(args.csv).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[input]  {csv_path}")
    print(f"[output] {outdir}")
    df = load_clas12_csv(csv_path)
    print(f"[data]   {len(df)} valid 10.6-GeV unpolarized cross-section points")

    cache = outdir / "km15_bh_decomposition.csv"
    df = evaluate_km15_dataframe(df, args.ebeam, max(1, args.workers), cache, args.force_km15)

    finite_model = (
        np.isfinite(df["km15_ep"]) & np.isfinite(df["km15_bh"]) &
        np.isfinite(df["R_BH"]) & (df["km15_ep"] > 0) & (df["km15_bh"] > 0)
    )
    df = df.loc[finite_model].copy().reset_index(drop=True)
    print(f"[KM15]   {len(df)} points have finite positive model predictions")

    max_ep_relerr = float(np.nanmax(df["km15_ep_decomp_relerr"]))
    max_bh_relerr = float(np.nanmax(df["bh_quad_relerr"]))
    print(f"[check]  max KM15 EP decomposition relative error = {max_ep_relerr:.3e}")
    print(f"[check]  max BH quadratic reconstruction relative error = {max_bh_relerr:.3e}")
    if max_ep_relerr > 1.0e-7:
        raise RuntimeError(
            "KM15 BH+INT+DVCS decomposition does not reproduce th_KM15.predict. "
            "Check the installed Gepard version/conventions before fitting."
        )
    #endif
    if max_bh_relerr > 1.0e-7:
        raise RuntimeError(
            "Precomputed BH quadratic coefficients do not reproduce the nominal BH cross section."
        )
    #endif

    # BH-selection summary and nested flags.
    selection_records = []
    for cut in sorted(args.bh_cuts):
        flag = df["bh_delta"] <= cut
        df[f"BH_{int(round(100*cut))}pct"] = flag
        selection_records.append({
            "bh_cut": cut,
            "bh_fraction_min": 1-cut,
            "npts": int(flag.sum()),
            "xB_min": float(df.loc[flag, "xB"].min()) if flag.any() else np.nan,
            "xB_max": float(df.loc[flag, "xB"].max()) if flag.any() else np.nan,
            "Q2_min": float(df.loc[flag, "Q2"].min()) if flag.any() else np.nan,
            "Q2_max": float(df.loc[flag, "Q2"].max()) if flag.any() else np.nan,
            "t_min": float(df.loc[flag, "t_abs"].min()) if flag.any() else np.nan,
            "t_max": float(df.loc[flag, "t_abs"].max()) if flag.any() else np.nan,
        })
        print(f"[select] |1-R_BH| <= {100*cut:.0f}% : {int(flag.sum())} points")
    #endfor
    pd.DataFrame(selection_records).to_csv(outdir / "bh_selection_summary.csv", index=False)
    df.to_csv(outdir / "clas12_with_km15_bh_decomposition.csv", index=False)

    save_bh_selection_plots(df, args.bh_cuts, outdir)

    include_ptp = not args.no_point_to_point_systematics
    include_bh_sys = not args.no_bh_selection_systematic
    use_scale_nuisance = not args.no_scale_nuisance

    fit_results: List[FitResult] = []

    # Reproduce paper Fits 1-5: dipole F1/F2 for every BH threshold.
    for cut in sorted(args.bh_cuts):
        d = df.loc[df["bh_delta"] <= cut].copy()
        if args.tmin is not None:
            d = d.loc[d["t_abs"] >= args.tmin]
        #endif
        if args.tmax is not None:
            d = d.loc[d["t_abs"] <= args.tmax]
        #endif
        if len(d) < 10:
            print(f"[fit] skip dipole at {100*cut:.0f}%: only {len(d)} points")
            continue
        #endif
        fr = fit_paper_model(d, "dipole", cut, include_ptp, include_bh_sys,
                             use_scale_nuisance, args.scale_fallback)
        fr.name = f"Dipole_{int(round(100*cut))}pct"
        fit_results.append(fr)
        print(f"[fit] {fr.name:18s} N={fr.npts:4d} chi2/dof={fr.chi2_ndof:.3f} "
              f"rE={fr.rE_fm:.4f} rM={fr.rM_fm:.4f}")
    #endfor

    # Paper's more flexible Set-5 fits.
    cut = float(args.fit_bh_cut)
    d5 = df.loc[df["bh_delta"] <= cut].copy()
    if args.tmin is not None:
        d5 = d5.loc[d5["t_abs"] >= args.tmin]
    #endif
    if args.tmax is not None:
        d5 = d5.loc[d5["t_abs"] <= args.tmax]
    #endif
    if len(d5) < 10:
        raise RuntimeError(f"Preferred BH sample has only {len(d5)} points.")
    #endif

    d5 = d5.copy()
    d5["fit_sigma_default"] = fit_sigma_errors(
        d5, cut, include_ptp, include_bh_sys
    )
    d5.to_csv(outdir / "preferred_BH_sample.csv", index=False)

    save_selected_cross_section_pages(d5, outdir / "03_selected_cross_sections_vs_KM15_BH.pdf")

    for kind, label in [
        ("ppole_same_a", "Ppole_same_a"),
        ("ppole_same_p", "Ppole_same_P"),
        ("f2_kelly_ppole_f1", "Ppole_F1_F2Kelly"),
    ]:
        fr = fit_paper_model(d5, kind, cut, include_ptp, include_bh_sys,
                             use_scale_nuisance, args.scale_fallback)
        fr.name = label
        fit_results.append(fr)
        print(f"[fit] {fr.name:18s} N={fr.npts:4d} chi2/dof={fr.chi2_ndof:.3f} "
              f"rE={fr.rE_fm:.4f} rM={fr.rM_fm:.4f}")
    #endfor

    all_families = ["P1", "P2", "P3", "P4", "IP1", "IP2", "IP3", "IP4", "CF2", "CF3", "CF4"]
    if len(args.sachs_families) == 1 and args.sachs_families[0].lower() == "all":
        families = all_families
    else:
        families = [x.upper() for x in args.sachs_families]
        bad = [x for x in families if x not in all_families]
        if bad:
            raise ValueError(f"Unknown Sachs families: {bad}")
        #endif
    #endif

    for family in families:
        fr = fit_sachs_family(d5, family, cut, include_ptp, include_bh_sys,
                              use_scale_nuisance, args.scale_fallback, fix_gm_to_kelly=False)
        fit_results.append(fr)
        print(f"[fit] {fr.name:18s} N={fr.npts:4d} chi2/dof={fr.chi2_ndof:.3f} "
              f"rE={fr.rE_fm:.4f} rM={fr.rM_fm:.4f}")
        if args.fit_gm_kelly_variants:
            frk = fit_sachs_family(d5, family, cut, include_ptp, include_bh_sys,
                                   use_scale_nuisance, args.scale_fallback, fix_gm_to_kelly=True)
            fit_results.append(frk)
            print(f"[fit] {frk.name:18s} N={frk.npts:4d} chi2/dof={frk.chi2_ndof:.3f} "
                  f"rE={frk.rE_fm:.4f} rM={frk.rM_fm:.4f}")
        #endif
    #endfor

    records = [fitresult_to_record(r) for r in fit_results]
    pd.DataFrame(records).to_csv(outdir / "fit_results.csv", index=False)

    # Full covariance matrices and metadata in a machine-readable NPZ/JSON pair.
    arrays = {}
    metadata = []
    for i, fr in enumerate(fit_results):
        key = f"fit_{i:03d}"
        arrays[key + "_params"] = fr.params
        arrays[key + "_cov"] = fr.covariance
        metadata.append({
            "key": key, "name": fr.name, "category": fr.category,
            "param_names": fr.param_names, "bh_cut": fr.bh_cut, "meta": fr.meta,
        })
    #endfor
    np.savez_compressed(outdir / "fit_covariances.npz", **arrays)
    with open(outdir / "fit_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)
    #endwith

    save_paper_ff_plots(fit_results, outdir, cut)
    save_sachs_family_plots(fit_results, outdir, cut)
    save_radii_plot(fit_results, outdir, cut)
    save_bh_f1_f2_sensitivity(outdir, args.ebeam)

    # Preferred CF4 summary, if available.
    cf4 = next((r for r in fit_results if r.name == "CF4" and abs(r.bh_cut-cut) < 1e-12), None)
    if cf4 is not None:
        print("\n[CF4 preferred direct Sachs fit]")
        print(f"  chi2/dof = {cf4.chi2:.3f}/{cf4.ndof} = {cf4.chi2_ndof:.3f}")
        print(f"  rE = {cf4.rE_fm:.5f} +/- {cf4.rE_err_fm:.5f} fm")
        print(f"  rM = {cf4.rM_fm:.5f} +/- {cf4.rM_err_fm:.5f} fm")
        print(f"  N  = {cf4.norm:.5f} +/- {cf4.norm_err:.5f}")
    #endif

    print(f"\nDone. Results are in {outdir}")
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif
