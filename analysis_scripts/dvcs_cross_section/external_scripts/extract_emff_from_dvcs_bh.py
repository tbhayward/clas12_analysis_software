#!/usr/bin/env python3
"""
CLAS12 reproduction of the analysis in

  Moradi et al., "Determination of proton electromagnetic form factors
  from DVCS measurements", arXiv:2512.06554.

This version intentionally removes the Hayward et al. (2018) P/IP/CF Sachs
fit families and follows the fitting strategy of arXiv:2512.06554 directly.

Workflow
--------
1. Read the preliminary CLAS12 10.6-GeV unpolarized ep -> e'p'gamma
   four-fold cross sections.
2. Evaluate the KM15 full EP, BH, DVCS, and interference contributions at
   each measured CLAS12 bin center.
3. Construct the five nested BH-dominated samples used in the paper:
      |1 - sigma_BH/sigma_EP| <= 1%, 2%, 3%, 4%, 5%.
4. Add the paper's corresponding 1%, 2%, 3%, 4%, or 5% uncorrelated
   BH-selection uncertainty to each retained cross-section point.
5. Reproduce Fits 1-5 using
      F1(t) = (1 - t/aE)^(-2)
      F2(t) = kappa (1 - t/aM)^(-2).
6. On Set 5 reproduce the P-pole studies:
      Fit 6: aE = aM = a, with PE and PM free.
      Fit 7: PE = PM = P, with aE and aM free.
7. Reproduce Fit 8:
      F2(t) fixed to Kelly,
      F1(t) = (1 - t/aE)^(-PE).
8. Convert F1,F2 -> GE,GM and extract rE,rM from their t=0 slopes.
9. Produce paper-style FF/radius plots plus CLAS12 observable-space plots
   showing the measured cross-section points together with the fitted BH
   predictions. No pseudo "CLAS12 GE/GM data points" are constructed.

Fit implementation
------------------
The paper uses iMinuit, the standard Hessian errors, and Delta chi^2 = 1
for 68% confidence intervals. This script does the same.

By default the CLAS12 fit uncertainty treatment is

  delta_i^2 =
      delta_stat,i^2
    + delta_point-to-point,i^2
    + (c_BH * sigma_i)^2,

where c_BH is 0.01 ... 0.05 for Fits 1 ... 5, following the additional
uncorrelated uncertainty prescription of arXiv:2512.06554.

The CSV "combination sys" is an overall fractional normalization
systematic. It is therefore treated as a correlated multiplicative
normalization nuisance N rather than being added independently to every
cross-section point:

  chi2 =
      sum_i [(N * sigma_BH,i - sigma_i)/delta_i]^2
    + [(N - 1)/delta_N]^2 .

The default delta_N is the median combination-systematic fraction of the
selected data set. Both the point-to-point and combination-normalization
systematics are included by default.

Efficiency
----------
KM15 is evaluated only once per measured point and cached. The BH cross
section is then represented exactly as

  sigma_BH = A F1^2 + B F1 F2 + C F2^2,

so the Minuit fits use only fast vectorized NumPy operations.

Default paths assume this script is run from

  /u/home/thayward/clas12_analysis_software/analysis_scripts/
      dvcs_cross_section/external_scripts/

Input:
  ../output/csvs/dvcs_pass2_analysis.csv

Output:
  output/emff_from_bh_paper_method_v5/
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
import matplotlib.pyplot as plt

try:
    from iminuit import Minuit
except ImportError:
    Minuit = None
#endtry
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
DEFAULT_OUTDIR = "output/emff_from_bh_paper_method_v5"

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
# Paper-method fit helpers
# -----------------------------------------------------------------------------
def numerical_gradient(fun: Callable[[np.ndarray], float],
                       p: np.ndarray) -> np.ndarray:
    """Central finite-difference gradient used only for propagated quantities."""
    p = np.asarray(p, dtype=float)
    grad = np.zeros_like(p)
    for i in range(len(p)):
        step = 1.0e-5 * max(1.0, abs(p[i]))
        pp = p.copy()
        pm = p.copy()
        pp[i] += step
        pm[i] -= step
        grad[i] = (fun(pp) - fun(pm)) / (2.0 * step)
    #endfor
    return grad
#enddef


def propagate_scalar(fun: Callable[[np.ndarray], float],
                     p: np.ndarray,
                     cov: np.ndarray) -> Tuple[float, float]:
    value = float(fun(p))
    if cov.size == 0 or not np.all(np.isfinite(cov)):
        return value, np.nan
    #endif
    grad = numerical_gradient(fun, p)
    var = float(grad @ cov @ grad)
    return value, math.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan
#enddef


def radius_from_shape(shape_fun: Callable[[np.ndarray], np.ndarray],
                      g0: float) -> float:
    """
    RMS radius in fm from a form factor written as a function of q=|t|=-t.

      r^2 = -6/G(0) dG/dq |_(q=0)

    which is equivalent to the paper's r^2 = +6/G(0) dG/dt at t=0.
    """
    h = 1.0e-6
    g_at_0 = float(shape_fun(np.array([0.0]))[0])
    g_at_h = float(shape_fun(np.array([h]))[0])
    slope_q = (g_at_h - g_at_0) / h
    r2_gev2 = -6.0 * slope_q / g0
    if not np.isfinite(r2_gev2) or r2_gev2 < 0.0:
        return np.nan
    #endif
    return math.sqrt(r2_gev2) * HBARC
#enddef


def fit_sigma_errors(data: pd.DataFrame,
                     bh_cut: float,
                     include_clas12_ptp_sys: bool,
                     include_bh_sys: bool) -> np.ndarray:
    """
    Pointwise uncertainty for the paper-method fit.

    Nominal CLAS12 treatment:
      sqrt(stat^2 + ptp_total^2 + (BH_cut * sigma)^2)
    """
    stat = data["xs_stat"].to_numpy(float)
    var = stat**2

    if include_clas12_ptp_sys:
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


# -----------------------------------------------------------------------------
# Moradi et al. F1/F2 models: Fits 1-8
# -----------------------------------------------------------------------------
def paper_model_f1f2(kind: str,
                     q: np.ndarray,
                     pars: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Evaluate the exact fit forms used in arXiv:2512.06554.

    q = |t| = -t, hence (1 - t/a) = (1 + q/a).
    """
    q = np.asarray(q, dtype=float)
    kind = kind.lower()

    if kind == "dipole":
        aE, aM = pars
        f1 = (1.0 + q / aE)**(-2.0)
        f2 = KAPPA_P * (1.0 + q / aM)**(-2.0)
        return f1, f2
    #endif

    if kind == "fit6_same_a":
        a, pE, pM = pars
        f1 = (1.0 + q / a)**(-pE)
        f2 = KAPPA_P * (1.0 + q / a)**(-pM)
        return f1, f2
    #endif

    if kind == "fit7_same_p":
        aE, aM, p = pars
        f1 = (1.0 + q / aE)**(-p)
        f2 = KAPPA_P * (1.0 + q / aM)**(-p)
        return f1, f2
    #endif

    if kind == "fit8_f2_kelly":
        aE, pE = pars
        f1 = (1.0 + q / aE)**(-pE)
        _, f2 = kelly_f1_f2(q)
        return f1, f2
    #endif

    raise ValueError(f"Unknown paper fit kind: {kind}")
#enddef


def paper_model_setup(kind: str) -> Tuple[List[str], np.ndarray]:
    """Initial values are seeded near the paper's reported solutions."""
    kind = kind.lower()

    if kind == "dipole":
        return ["aE", "aM"], np.array([0.95, 0.56], dtype=float)
    #endif

    if kind == "fit6_same_a":
        return ["a", "PE", "PM"], np.array([0.53, 1.2, 1.96], dtype=float)
    #endif

    if kind == "fit7_same_p":
        return ["aE", "aM", "P"], np.array([0.65, 0.35, 1.42], dtype=float)
    #endif

    if kind == "fit8_f2_kelly":
        return ["aE", "PE"], np.array([0.565, 1.318], dtype=float)
    #endif

    raise ValueError(f"Unknown paper fit kind: {kind}")
#enddef


def fit_paper_model(data: pd.DataFrame,
                    kind: str,
                    fit_name: str,
                    bh_cut: float,
                    include_clas12_ptp_sys: bool,
                    include_bh_sys: bool,
                    include_combination_norm_sys: bool) -> FitResult:
    """
    Perform one Moradi-et-al.-style chi^2 fit with iMinuit, augmented by
    the CLAS12 correlated run-period-combination normalization systematic.

    The normalization nuisance N multiplies the fitted BH prediction:

      chi2 = sum_i [(N * sigma_BH,i - sigma_i)/err_i]^2
             + [(N - 1)/sigma_N]^2

    where sigma_N is taken from the median fractional "combination sys"
    value in the selected CLAS12 sample.
    """
    q = data["t_abs"].to_numpy(float)
    y = data["xs"].to_numpy(float)
    err = fit_sigma_errors(
        data,
        bh_cut,
        include_clas12_ptp_sys,
        include_bh_sys,
    )
    A = data["bh_A"].to_numpy(float)
    B = data["bh_B"].to_numpy(float)
    C = data["bh_C"].to_numpy(float)

    names, p0 = paper_model_setup(kind)

    comb = data["comb_sys_frac"].to_numpy(float)
    good_comb = comb[np.isfinite(comb) & (comb > 0.0)]
    if len(good_comb) > 0:
        sigma_norm = float(np.median(good_comb))
    else:
        sigma_norm = 0.0
    #endif

    if include_combination_norm_sys and sigma_norm > 0.0:
        fit_names = names + ["N"]
        fit_p0 = np.concatenate([p0, np.array([1.0])])
    else:
        fit_names = names
        fit_p0 = p0.copy()
    #endif

    def split_params(p: np.ndarray) -> Tuple[np.ndarray, float]:
        if include_combination_norm_sys and sigma_norm > 0.0:
            return p[:-1], float(p[-1])
        #endif
        return p, 1.0
    #enddef

    def chi2_array(p: np.ndarray) -> float:
        ff_pars, norm = split_params(np.asarray(p, dtype=float))
        f1, f2 = paper_model_f1f2(kind, q, ff_pars)
        pred = norm * bh_from_f1f2(A, B, C, f1, f2)
        pull = (pred - y) / err
        chi2 = float(np.dot(pull, pull))

        if include_combination_norm_sys and sigma_norm > 0.0:
            chi2 += ((norm - 1.0) / sigma_norm)**2
        #endif

        return chi2
    #enddef

    def minuit_chi2(*values):
        return chi2_array(np.asarray(values, dtype=float))
    #enddef

    m = Minuit(minuit_chi2, *fit_p0, name=tuple(fit_names))
    m.errordef = Minuit.LEAST_SQUARES

    for name in names:
        m.limits[name] = (1.0e-6, None)
    #endfor

    if "N" in fit_names:
        # Keep N positive but otherwise let the Gaussian prior determine it.
        m.limits["N"] = (1.0e-6, None)
    #endif

    m.migrad()
    m.hesse()

    all_params = np.array(
        [float(m.values[name]) for name in fit_names],
        dtype=float,
    )

    full_cov = np.full((len(fit_names), len(fit_names)), np.nan, dtype=float)
    if m.covariance is not None:
        for i, ni in enumerate(fit_names):
            for j, nj in enumerate(fit_names):
                full_cov[i, j] = float(m.covariance[ni, nj])
            #endfor
        #endfor
    #endif

    if "N" in fit_names:
        params = all_params[:-1]
        cov = full_cov[:-1, :-1]
        norm = float(all_params[-1])
        norm_err = math.sqrt(max(full_cov[-1, -1], 0.0))
    else:
        params = all_params
        cov = full_cov
        norm = 1.0
        norm_err = 0.0
    #endif

    chi2 = float(m.fval)
    npar = len(all_params)
    ndof = max(len(y) - npar, 1)

    def re_of_p(p):
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(
                qq, *paper_model_f1f2(kind, qq, p)
            )[0],
            1.0,
        )
    #enddef

    def rm_of_p(p):
        return radius_from_shape(
            lambda qq: f1f2_to_sachs(
                qq, *paper_model_f1f2(kind, qq, p)
            )[1],
            MU_P,
        )
    #enddef

    rE, drE = propagate_scalar(re_of_p, params, cov)
    rM, drM = propagate_scalar(rm_of_p, params, cov)

    message = (
        f"valid={m.valid}; accurate={m.accurate}; "
        f"edm={float(m.fmin.edm):.3e}; "
        f"has_posdef_covar={bool(m.fmin.has_posdef_covar)}"
    )

    return FitResult(
        name=fit_name,
        category="paper_F1F2",
        bh_cut=bh_cut,
        npts=len(y),
        npar=npar,
        chi2=chi2,
        ndof=ndof,
        chi2_ndof=chi2 / ndof,
        success=bool(m.valid),
        message=message,
        params=params,
        param_names=names,
        covariance=cov,
        rE_fm=rE,
        rE_err_fm=drE,
        rM_fm=rM,
        rM_err_fm=drM,
        norm=norm,
        norm_err=norm_err,
        model_kind=kind,
        meta={
            "minuit_valid": bool(m.valid),
            "minuit_accurate": bool(m.accurate),
            "minuit_edm": float(m.fmin.edm),
            "has_posdef_covar": bool(m.fmin.has_posdef_covar),
            "include_clas12_ptp_sys": bool(include_clas12_ptp_sys),
            "include_bh_selection_sys": bool(include_bh_sys),
            "include_combination_norm_sys": bool(include_combination_norm_sys),
            "combination_norm_sigma": float(sigma_norm),
        },
    )
#enddef


def evaluate_fit_form_factors(fr: FitResult,
                              q: np.ndarray) -> Tuple[np.ndarray, np.ndarray,
                                                      np.ndarray, np.ndarray]:
    """Return F1, F2, GE, GM for a paper-method FitResult."""
    f1, f2 = paper_model_f1f2(fr.model_kind, q, fr.params)
    ge, gm = f1f2_to_sachs(q, f1, f2)
    return f1, f2, ge, gm
#enddef


def evaluate_fit_cross_section(fr: FitResult,
                               data: pd.DataFrame) -> np.ndarray:
    """Evaluate a fitted pure-BH prediction at each measured CLAS12 point."""
    q = data["t_abs"].to_numpy(float)
    f1, f2 = paper_model_f1f2(fr.model_kind, q, fr.params)
    return fr.norm * bh_from_f1f2(
        data["bh_A"].to_numpy(float),
        data["bh_B"].to_numpy(float),
        data["bh_C"].to_numpy(float),
        f1,
        f2,
    )
#enddef


# -----------------------------------------------------------------------------
# Fit bookkeeping and Hessian uncertainty bands
# -----------------------------------------------------------------------------
def fitresult_to_record(fr: FitResult) -> Dict[str, object]:
    rec = {
        "name": fr.name,
        "category": fr.category,
        "bh_cut": fr.bh_cut,
        "npts": fr.npts,
        "npar": fr.npar,
        "chi2": fr.chi2,
        "ndof": fr.ndof,
        "chi2_ndof": fr.chi2_ndof,
        "success": fr.success,
        "message": fr.message,
        "rE_fm": fr.rE_fm,
        "rE_err_fm": fr.rE_err_fm,
        "rM_fm": fr.rM_fm,
        "rM_err_fm": fr.rM_err_fm,
        "norm": fr.norm,
        "norm_err": fr.norm_err,
        "combination_norm_sigma": fr.meta.get("combination_norm_sigma", np.nan),
    }
    for name, value in zip(fr.param_names, fr.params):
        rec[name] = value
    #endfor
    return rec
#enddef


def form_factor_band(fr: FitResult,
                     q: np.ndarray,
                     which: str) -> Tuple[np.ndarray, np.ndarray]:
    """Central value and 1-sigma Hessian band."""
    index = {"F1": 0, "F2": 1, "GE": 2, "GM": 3}[which]
    central = evaluate_fit_form_factors(fr, q)[index]

    if not np.all(np.isfinite(fr.covariance)):
        return central, np.full_like(central, np.nan)
    #endif

    sigma = np.zeros_like(q, dtype=float)
    for iq, q_value in enumerate(q):
        def scalar(p):
            tmp = FitResult(**{**fr.__dict__, "params": np.asarray(p, dtype=float)})
            return float(
                evaluate_fit_form_factors(
                    tmp, np.array([q_value], dtype=float)
                )[index][0]
            )
        #enddef

        grad = numerical_gradient(scalar, fr.params)
        var = float(grad @ fr.covariance @ grad)
        sigma[iq] = math.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan
    #endfor

    return central, sigma
#enddef


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def shade_fit_t_range(ax, data: pd.DataFrame) -> None:
    """Mark the |t| interval actually constrained by the selected CLAS12 data."""
    if len(data) == 0:
        return
    #endif
    ax.axvspan(
        float(data["t_abs"].min()),
        float(data["t_abs"].max()),
        alpha=0.08,
        label="CLAS12 fitted |t| range",
    )
#enddef


def save_bh_selection_plots(df: pd.DataFrame,
                            bh_cuts: Sequence[float],
                            outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    vals = df["R_BH"].to_numpy(float)
    vals = vals[np.isfinite(vals)]
    ax.hist(vals, bins=80, histtype="step", linewidth=1.5)

    for cut in bh_cuts:
        ax.axvline(
            1.0 - cut,
            linestyle="--",
            linewidth=1.0,
            label=rf"$\Delta\sigma/\sigma\leq{100*cut:.0f}\%$",
        )
        ax.axvline(1.0 + cut, linestyle="--", linewidth=1.0)
    #endfor

    ax.set_xlabel(r"$R_{\rm BH}=\sigma_{\rm BH}/\sigma_{\rm EP}^{\rm KM15}$")
    ax.set_ylabel("CLAS12 cross-section points")
    ax.set_title("KM15 BH-dominance selection")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "01_bh_fraction_distribution.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    sc = ax.scatter(
        df["phi_deg"],
        df["xB"],
        c=df["R_BH"],
        s=15,
        vmin=0.8,
        vmax=1.05,
    )
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$R_{\rm BH}$")
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$x_B$")
    ax.set_title("Preliminary CLAS12 kinematics colored by KM15 BH fraction")
    fig.tight_layout()
    fig.savefig(outdir / "02_bh_fraction_phi_xB.pdf")
    plt.close(fig)
#enddef


def save_selected_cross_section_pages(data: pd.DataFrame,
                                      fit5: FitResult,
                                      fit8: FitResult,
                                      outpath: Path) -> None:
    """
    Observable-space CLAS12 plots.

    These are the appropriate place to show the CLAS12 data points: the
    measured quantity is the ep->e'p'gamma cross section, not F1/F2/GE/GM
    point-by-point. Fit 5 and Fit 8 BH predictions are overlaid.
    """
    work = data.copy()
    work["fit5_bh"] = evaluate_fit_cross_section(fit5, work)
    work["fit8_bh"] = evaluate_fit_cross_section(fit8, work)

    group_cols = [
        c for c in
        ["xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max"]
        if c in work.columns
    ]

    if len(group_cols) < 6:
        work["_gx"] = work["xB"].round(3)
        work["_gq"] = work["Q2"].round(2)
        work["_gt"] = work["t_abs"].round(2)
        group_cols = ["_gx", "_gq", "_gt"]
    #endif

    groups = [
        (key, group.sort_values("phi_deg"))
        for key, group in work.groupby(group_cols, dropna=False)
    ]

    with PdfPages(outpath) as pdf:
        for start in range(0, len(groups), 6):
            chunk = groups[start:start + 6]
            fig, axes = plt.subplots(3, 2, figsize=(11, 11), squeeze=False)

            for ax, item in zip(axes.flat, chunk):
                _, group = item
                ax.errorbar(
                    group["phi_deg"],
                    group["xs"],
                    yerr=group["fit_sigma_default"],
                    fmt="o",
                    markersize=3.5,
                    capsize=2,
                    label="CLAS12 preliminary",
                )
                ax.plot(
                    group["phi_deg"],
                    group["fit5_bh"],
                    "-",
                    linewidth=1.4,
                    label="Fit 5 BH",
                )
                ax.plot(
                    group["phi_deg"],
                    group["fit8_bh"],
                    "--",
                    linewidth=1.4,
                    label=r"Fit 8 BH ($F_2$ Kelly)",
                )
                ax.plot(
                    group["phi_deg"],
                    group["km15_bh"],
                    ":",
                    linewidth=1.1,
                    label="KM15 nominal BH",
                )
                ax.set_xlabel(r"$\phi$ (deg)")
                ax.set_ylabel(r"$d^4\sigma$")
                ax.set_title(
                    rf"$\langle x_B\rangle={group['xB'].mean():.3f}$, "
                    rf"$\langle Q^2\rangle={group['Q2'].mean():.2f}$, "
                    rf"$\langle|t|\rangle={group['t_abs'].mean():.3f}$"
                )
                ax.grid(alpha=0.2)
            #endfor

            for ax in axes.flat[len(chunk):]:
                ax.axis("off")
            #endfor

            handles, labels = axes.flat[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=8)
            fig.suptitle(
                "Set 5: preliminary CLAS12 cross sections and fitted BH predictions",
                y=0.995,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.965))
            pdf.savefig(fig)
            plt.close(fig)
        #endfor
    #endwith
#enddef


def save_fit1_to_fit5_plots(results: List[FitResult],
                            set5_data: pd.DataFrame,
                            outdir: Path) -> None:
    """Paper Fig. 4/5 analogue for CLAS12 Fits 1-5."""
    fits = [
        next(r for r in results if r.name == f"Fit {i}")
        for i in range(1, 6)
    ]
    fit5 = fits[-1]

    q = np.linspace(0.0, 1.0, 400)

    for which, ylabel, filename in [
        ("F1", r"$F_1(t)$", "04_F1_fits_1_to_5.pdf"),
        ("F2", r"$F_2(t)$", "05_F2_fits_1_to_5.pdf"),
    ]:
        fig, ax = plt.subplots(figsize=(7.2, 5.2))

        for fr in fits:
            values = evaluate_fit_form_factors(fr, q)[0 if which == "F1" else 1]
            ax.plot(q, values, linewidth=1.2, label=fr.name)
        #endfor

        central, sigma = form_factor_band(fit5, q, which)
        ax.fill_between(
            q,
            central - sigma,
            central + sigma,
            alpha=0.20,
            label="Fit 5 68% Hessian band",
        )

        shade_fit_t_range(ax, set5_data)
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.0, 1.0)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir / filename)
        plt.close(fig)
    #endfor
#enddef


def save_fit5_sachs_plot(fit5: FitResult,
                         set5_data: pd.DataFrame,
                         outdir: Path) -> None:
    """Paper Fig. 6 analogue using CLAS12 Fit 5."""
    q = np.linspace(0.0, 1.0, 400)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

    for ax, which, ylabel in zip(
        axes,
        ["GE", "GM"],
        [r"$G_E(t)$", r"$G_M(t)$"],
    ):
        central, sigma = form_factor_band(fit5, q, which)
        ax.plot(q, central, linewidth=1.5, label="Fit 5")
        ax.fill_between(
            q,
            central - sigma,
            central + sigma,
            alpha=0.20,
            label="68% Hessian band",
        )
        shade_fit_t_range(ax, set5_data)
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.0, 1.0)
        ax.grid(alpha=0.2)
    #endfor

    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "06_GE_GM_fit5.pdf")
    plt.close(fig)
#enddef


def save_bh_f1_f2_sensitivity(outdir: Path,
                              ebeam: float,
                              reference_data: pd.DataFrame) -> None:
    """
    Paper Fig. 7 analogue.

    Use a representative CLAS12 Set-5 kinematic point rather than importing
    elastic-scattering data. The point chosen is the Set-5 event closest to
    the medians in xB, Q2, and phi.
    """
    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception:
        return
    #endtry

    med_x = float(reference_data["xB"].median())
    med_q2 = float(reference_data["Q2"].median())
    med_phi = float(reference_data["phi_deg"].median())

    scale_x = max(float(reference_data["xB"].std()), 1.0e-6)
    scale_q2 = max(float(reference_data["Q2"].std()), 1.0e-6)
    scale_phi = max(float(reference_data["phi_deg"].std()), 1.0e-6)

    metric = (
        ((reference_data["xB"] - med_x) / scale_x)**2
        + ((reference_data["Q2"] - med_q2) / scale_q2)**2
        + ((reference_data["phi_deg"] - med_phi) / scale_phi)**2
    )
    row = reference_data.loc[metric.idxmin()]

    qmax = max(0.45, float(reference_data["t_abs"].max()))
    q = np.linspace(0.01, qmax, 220)
    full = np.zeros_like(q)
    f1zero = np.zeros_like(q)
    f2zero = np.zeros_like(q)

    for i, qi in enumerate(q):
        th = th_KM15
        pt = make_gepard_point(
            g,
            float(row["xB"]),
            float(row["Q2"]),
            float(qi),
            float(row["phi_deg"]),
            ebeam,
        )
        pref = float(th.PreFacSigma(pt))
        model = th.m
        old_f1, old_f2 = model.F1, model.F2
        f1k, f2k = kelly_f1_f2(np.array([qi]))

        try:
            model.F1 = lambda p, value=float(f1k[0]): value
            model.F2 = lambda p, value=float(f2k[0]): value
            full[i] = pref * float(th.TBH2unp(pt))

            model.F1 = lambda p: 0.0
            f1zero[i] = pref * float(th.TBH2unp(pt))

            model.F1 = lambda p, value=float(f1k[0]): value
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
    ax.set_title(
        rf"BH sensitivity: $x_B={row['xB']:.3f}$, "
        rf"$Q^2={row['Q2']:.2f}$ GeV$^2$, "
        rf"$\phi={row['phi_deg']:.1f}^\circ$"
    )
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "07_BH_F1_F2_sensitivity.pdf")
    plt.close(fig)
#enddef


def save_fit5_to_fit8_sachs_plot(results: List[FitResult],
                                 set5_data: pd.DataFrame,
                                 outdir: Path) -> None:
    """Paper Fig. 8 analogue: Fits 5-8 in GE and GM."""
    names = ["Fit 5", "Fit 6", "Fit 7", "Fit 8"]
    fits = [next(r for r in results if r.name == name) for name in names]
    q = np.linspace(0.0, 1.0, 400)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

    for fr in fits:
        _, _, ge, gm = evaluate_fit_form_factors(fr, q)
        axes[0].plot(q, ge, linewidth=1.3, label=fr.name)
        axes[1].plot(q, gm, linewidth=1.3, label=fr.name)
    #endfor

    # Show only the Fit-8 uncertainty band to keep the comparison legible.
    for ax, which in zip(axes, ["GE", "GM"]):
        central, sigma = form_factor_band(fits[-1], q, which)
        ax.fill_between(
            q,
            central - sigma,
            central + sigma,
            alpha=0.18,
            label="Fit 8 68% Hessian band",
        )
        shade_fit_t_range(ax, set5_data)
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_xlim(0.0, 1.0)
        ax.grid(alpha=0.2)
    #endfor

    axes[0].set_ylabel(r"$G_E(t)$")
    axes[1].set_ylabel(r"$G_M(t)$")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "08_GE_GM_fits_5_to_8.pdf")
    plt.close(fig)
#enddef


def save_radii_plot(results: List[FitResult],
                    outdir: Path) -> None:
    """Paper Table-V-style visual summary for Fits 5-8."""
    rows = [
        next(r for r in results if r.name == f"Fit {i}")
        for i in range(5, 9)
    ]
    labels = [r.name for r in rows]
    y = np.arange(len(rows))

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharey=True)

    axes[0].errorbar(
        [r.rE_fm for r in rows],
        y,
        xerr=[r.rE_err_fm for r in rows],
        fmt="o",
        capsize=2,
    )
    axes[1].errorbar(
        [r.rM_fm for r in rows],
        y,
        xerr=[r.rM_err_fm for r in rows],
        fmt="o",
        capsize=2,
    )

    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel(r"$r_E$ (fm)")
    axes[1].set_xlabel(r"$r_M$ (fm)")
    axes[0].grid(alpha=0.2)
    axes[1].grid(alpha=0.2)

    fig.tight_layout()
    fig.savefig(outdir / "09_radii_fits_5_to_8.pdf")
    plt.close(fig)
#enddef


def save_summary_tables(results: List[FitResult],
                        selection_records: List[Dict[str, float]],
                        outdir: Path) -> None:
    """Write compact paper-like Tables I-V to CSV files."""
    pd.DataFrame(selection_records).to_csv(
        outdir / "table_I_bh_datasets.csv",
        index=False,
    )

    fits15 = [
        next(r for r in results if r.name == f"Fit {i}")
        for i in range(1, 6)
    ]
    pd.DataFrame([
        {
            "Fit": fr.name,
            "BH criterion": f"<= {100*fr.bh_cut:.0f}%",
            "Npts": fr.npts,
            "chi2": fr.chi2,
            "dof": fr.ndof,
            "chi2/dof": fr.chi2_ndof,
        }
        for fr in fits15
    ]).to_csv(outdir / "table_II_fits_1_to_5_chi2.csv", index=False)

    pd.DataFrame([
        {
            "Fit": fr.name,
            "aE_GeV2": fr.params[0],
            "aE_err_GeV2": math.sqrt(fr.covariance[0, 0]),
            "aM_GeV2": fr.params[1],
            "aM_err_GeV2": math.sqrt(fr.covariance[1, 1]),
        }
        for fr in fits15
    ]).to_csv(outdir / "table_III_dipole_parameters.csv", index=False)

    fit6 = next(r for r in results if r.name == "Fit 6")
    fit7 = next(r for r in results if r.name == "Fit 7")
    fit8 = next(r for r in results if r.name == "Fit 8")

    rows = []
    for fr in [fit6, fit7, fit8]:
        row = {
            "Fit": fr.name,
            "chi2/dof": fr.chi2_ndof,
        }
        for i, name in enumerate(fr.param_names):
            row[name] = fr.params[i]
            row[name + "_err"] = math.sqrt(
                max(fr.covariance[i, i], 0.0)
            )
        #endfor
        rows.append(row)
    #endfor
    pd.DataFrame(rows).to_csv(
        outdir / "table_IV_ppole_and_fit8_parameters.csv",
        index=False,
    )

    pd.DataFrame([
        {
            "Fit": fr.name,
            "rE_fm": fr.rE_fm,
            "rE_err_fm": fr.rE_err_fm,
            "rM_fm": fr.rM_fm,
            "rM_err_fm": fr.rM_err_fm,
        }
        for fr in [fits15[-1], fit6, fit7, fit8]
    ]).to_csv(outdir / "table_V_radii.csv", index=False)
#enddef


# -----------------------------------------------------------------------------
# CLI / main
# -----------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        "--csv",
        default=DEFAULT_CSV,
        help="CLAS12 DVCS cross-section CSV",
    )
    p.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help="Output directory",
    )
    p.add_argument(
        "--ebeam",
        type=float,
        default=10.604,
        help="Beam energy in GeV",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=min(6, os.cpu_count() or 1),
        help="KM15 worker processes",
    )
    p.add_argument(
        "--bh-cuts",
        type=float,
        nargs="+",
        default=[0.01, 0.02, 0.03, 0.04, 0.05],
        help="Nested |1-R_BH| thresholds; nominal paper values are 1-5%",
    )
    p.add_argument(
        "--force-km15",
        action="store_true",
        help="Ignore the cached KM15 decomposition",
    )
    p.add_argument(
        "--no-clas12-point-to-point-systematic",
        action="store_true",
        help=(
            "Disable the preliminary CLAS12 absolute point-to-point "
            "systematic. It is included by default."
        ),
    )
    p.add_argument(
        "--no-combination-normalization-systematic",
        action="store_true",
        help=(
            "Disable the correlated normalization nuisance from the CSV "
            "combination systematic. It is included by default."
        ),
    )
    p.add_argument(
        "--no-bh-selection-systematic",
        action="store_true",
        help="Disable the paper's additional 1-5% uncorrelated uncertainty",
    )
    p.add_argument(
        "--tmin",
        type=float,
        default=None,
        help="Optional minimum |t|; normally leave unset",
    )
    p.add_argument(
        "--tmax",
        type=float,
        default=None,
        help="Optional maximum |t|; normally leave unset",
    )
    return p
#enddef


def apply_optional_t_range(data: pd.DataFrame,
                           tmin: Optional[float],
                           tmax: Optional[float]) -> pd.DataFrame:
    result = data.copy()
    if tmin is not None:
        result = result.loc[result["t_abs"] >= tmin]
    #endif
    if tmax is not None:
        result = result.loc[result["t_abs"] <= tmax]
    #endif
    return result.copy()
#enddef


def main() -> int:
    args = build_parser().parse_args()

    # Require the paper's minimizer explicitly.
    if Minuit is None:
        raise RuntimeError(
            "iminuit is required because arXiv:2512.06554 performs the "
            "chi^2 minimization with iMinuit and Hessian Delta-chi^2=1 "
            "errors. Install it in the active Python environment, e.g. "
            "`python -m pip install --user iminuit`, then rerun."
        )
    #endif

    csv_path = Path(args.csv).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[input]  {csv_path}")
    print(f"[output] {outdir}")

    df = load_clas12_csv(csv_path)
    print(
        f"[data]   {len(df)} valid 10.6-GeV unpolarized "
        "cross-section points"
    )

    # Reuse the validated cache from the earlier analysis if present.
    old_cache = Path("output/emff_from_bh/km15_bh_decomposition.csv").resolve()
    cache = outdir / "km15_bh_decomposition.csv"
    if (not cache.exists()) and old_cache.exists() and (not args.force_km15):
        cache.parent.mkdir(parents=True, exist_ok=True)
        import shutil
        shutil.copy2(old_cache, cache)
        print(f"[KM15] Copied existing validated cache: {old_cache}")
    #endif

    df = evaluate_km15_dataframe(
        df,
        args.ebeam,
        max(1, args.workers),
        cache,
        args.force_km15,
    )

    finite_model = (
        np.isfinite(df["km15_ep"])
        & np.isfinite(df["km15_bh"])
        & np.isfinite(df["R_BH"])
        & (df["km15_ep"] > 0)
        & (df["km15_bh"] > 0)
    )
    df = df.loc[finite_model].copy().reset_index(drop=True)
    print(f"[KM15]   {len(df)} points have finite positive model predictions")

    max_ep_relerr = float(np.nanmax(df["km15_ep_decomp_relerr"]))
    max_bh_relerr = float(np.nanmax(df["bh_quad_relerr"]))
    print(
        "[check]  max KM15 EP decomposition relative error "
        f"= {max_ep_relerr:.3e}"
    )
    print(
        "[check]  max BH quadratic reconstruction relative error "
        f"= {max_bh_relerr:.3e}"
    )

    if max_ep_relerr > 1.0e-7:
        raise RuntimeError(
            "KM15 BH+INT+DVCS decomposition does not reproduce "
            "th_KM15.predict."
        )
    #endif

    if max_bh_relerr > 1.0e-7:
        raise RuntimeError(
            "Precomputed BH quadratic coefficients do not reproduce "
            "the nominal BH cross section."
        )
    #endif

    bh_cuts = sorted(float(cut) for cut in args.bh_cuts)
    if bh_cuts != [0.01, 0.02, 0.03, 0.04, 0.05]:
        warnings.warn(
            "The paper defines Fits 1-5 with 1%,2%,3%,4%,5% cuts. "
            "Nonstandard --bh-cuts change that correspondence."
        )
    #endif

    selection_records = []
    selected_sets: Dict[float, pd.DataFrame] = {}

    for cut in bh_cuts:
        selected = df.loc[df["bh_delta"] <= cut].copy()
        selected = apply_optional_t_range(selected, args.tmin, args.tmax)
        selected_sets[cut] = selected

        flag = df["bh_delta"] <= cut
        df[f"BH_{int(round(100*cut))}pct"] = flag

        selection_records.append({
            "Set": int(round(100 * cut)),
            "BH_cut": cut,
            "criterion": f"|1-R_BH| <= {100*cut:.0f}%",
            "Npts": len(selected),
            "xB_min": float(selected["xB"].min()) if len(selected) else np.nan,
            "xB_max": float(selected["xB"].max()) if len(selected) else np.nan,
            "Q2_min": float(selected["Q2"].min()) if len(selected) else np.nan,
            "Q2_max": float(selected["Q2"].max()) if len(selected) else np.nan,
            "t_min": float(selected["t_abs"].min()) if len(selected) else np.nan,
            "t_max": float(selected["t_abs"].max()) if len(selected) else np.nan,
            "phi_min": float(selected["phi_deg"].min()) if len(selected) else np.nan,
            "phi_max": float(selected["phi_deg"].max()) if len(selected) else np.nan,
        })

        print(
            f"[select] |1-R_BH| <= {100*cut:.0f}% : "
            f"{len(selected)} points"
        )
    #endfor

    df.to_csv(
        outdir / "clas12_with_km15_bh_decomposition.csv",
        index=False,
    )
    save_bh_selection_plots(df, bh_cuts, outdir)

    include_ptp = not args.no_clas12_point_to_point_systematic
    include_bh_sys = not args.no_bh_selection_systematic
    include_comb_norm = not args.no_combination_normalization_systematic

    print(
        "[errors] point uncertainty = "
        "sqrt(stat^2 + point-to-point^2 + BH-selection^2)"
        if include_ptp
        else "[errors] point uncertainty = sqrt(stat^2 + BH-selection^2)"
    )
    print(
        "[errors] CLAS12 point-to-point systematic: "
        + ("INCLUDED" if include_ptp else "NOT included")
    )
    print(
        "[norm]   combination systematic as correlated normalization nuisance: "
        + ("INCLUDED" if include_comb_norm else "NOT included")
    )

    fit_results: List[FitResult] = []

    # Fits 1-5: dipole F1 and F2 on the five nested BH samples.
    for fit_index, cut in enumerate(bh_cuts, start=1):
        data = selected_sets[cut]
        if len(data) < 10:
            raise RuntimeError(
                f"Fit {fit_index} has only {len(data)} selected points."
            )
        #endif

        fr = fit_paper_model(
            data=data,
            kind="dipole",
            fit_name=f"Fit {fit_index}",
            bh_cut=cut,
            include_clas12_ptp_sys=include_ptp,
            include_bh_sys=include_bh_sys,
            include_combination_norm_sys=include_comb_norm,
        )
        fit_results.append(fr)

        aE = fr.params[0]
        aM = fr.params[1]
        daE = math.sqrt(max(fr.covariance[0, 0], 0.0))
        daM = math.sqrt(max(fr.covariance[1, 1], 0.0))

        print(
            f"[fit] Fit {fit_index}: N={fr.npts:4d} "
            f"chi2/dof={fr.chi2_ndof:.3f} "
            f"aE={aE:.4f}+/-{daE:.4f} "
            f"aM={aM:.4f}+/-{daM:.4f} "
            f"N={fr.norm:.4f}+/-{fr.norm_err:.4f}"
        )
    #endfor

    if 0.05 not in selected_sets:
        raise RuntimeError(
            "The exact paper-method Fit 6-8 stage requires the 5% Set 5."
        )
    #endif

    set5 = selected_sets[0.05].copy()
    set5["fit_sigma_default"] = fit_sigma_errors(
        set5,
        0.05,
        include_ptp,
        include_bh_sys,
    )
    set5.to_csv(outdir / "set5_selected_points.csv", index=False)

    # Fits 6-8 exactly as in the paper.
    for kind, fit_name in [
        ("fit6_same_a", "Fit 6"),
        ("fit7_same_p", "Fit 7"),
        ("fit8_f2_kelly", "Fit 8"),
    ]:
        fr = fit_paper_model(
            data=set5,
            kind=kind,
            fit_name=fit_name,
            bh_cut=0.05,
            include_clas12_ptp_sys=include_ptp,
            include_bh_sys=include_bh_sys,
            include_combination_norm_sys=include_comb_norm,
        )
        fit_results.append(fr)

        parameter_text = " ".join(
            f"{name}={value:.4f}+/-"
            f"{math.sqrt(max(fr.covariance[i, i], 0.0)):.4f}"
            for i, (name, value) in enumerate(
                zip(fr.param_names, fr.params)
            )
        )
        print(
            f"[fit] {fit_name}: N={fr.npts:4d} "
            f"chi2/dof={fr.chi2_ndof:.3f} {parameter_text} "
            f"N={fr.norm:.4f}+/-{fr.norm_err:.4f}"
        )
    #endfor

    fit5 = next(r for r in fit_results if r.name == "Fit 5")
    fit8 = next(r for r in fit_results if r.name == "Fit 8")

    # Save machine-readable outputs.
    pd.DataFrame(
        [fitresult_to_record(fr) for fr in fit_results]
    ).to_csv(outdir / "fit_results.csv", index=False)

    arrays = {}
    metadata = []
    for i, fr in enumerate(fit_results):
        key = f"fit_{i + 1}"
        arrays[key + "_params"] = fr.params
        arrays[key + "_cov"] = fr.covariance
        metadata.append({
            "key": key,
            "name": fr.name,
            "param_names": fr.param_names,
            "bh_cut": fr.bh_cut,
            "model_kind": fr.model_kind,
            "meta": fr.meta,
        })
    #endfor

    np.savez_compressed(outdir / "fit_covariances.npz", **arrays)
    with open(outdir / "fit_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
    #endwith

    save_summary_tables(fit_results, selection_records, outdir)

    # Paper-style and observable-space plots.
    save_selected_cross_section_pages(
        set5,
        fit5,
        fit8,
        outdir / "03_set5_cross_sections_with_fits.pdf",
    )
    save_fit1_to_fit5_plots(fit_results, set5, outdir)
    save_fit5_sachs_plot(fit5, set5, outdir)
    save_bh_f1_f2_sensitivity(outdir, args.ebeam, set5)
    save_fit5_to_fit8_sachs_plot(fit_results, set5, outdir)
    save_radii_plot(fit_results, outdir)

    print("\n[paper-method radii]")
    for fit_name in ["Fit 5", "Fit 6", "Fit 7", "Fit 8"]:
        fr = next(r for r in fit_results if r.name == fit_name)
        print(
            f"  {fit_name}: "
            f"rE={fr.rE_fm:.5f}+/-{fr.rE_err_fm:.5f} fm, "
            f"rM={fr.rM_fm:.5f}+/-{fr.rM_err_fm:.5f} fm, "
            f"N={fr.norm:.5f}+/-{fr.norm_err:.5f}, "
            f"chi2/dof={fr.chi2_ndof:.3f}"
        )
    #endfor

    print(f"\nDone. Results are in {outdir}")
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif
