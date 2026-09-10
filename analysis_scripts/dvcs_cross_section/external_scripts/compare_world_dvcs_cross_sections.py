"""
Compare published unpolarized proton DVCS/BH cross-section measurements.

This is the standalone analysis driver for the external-data/model-comparison
chapter of the CLAS12 RGA pass-2 DVCS analysis note.  It deliberately does NOT
load the pass-2 result yet.  The first objective is to establish the mutual
consistency of the published world datasets before the new pass-2 measurement
is introduced.

Published datasets included
---------------------------
  * CLAS6  Jo et al.       2015
  * Hall A Defurne et al.  2015
  * Hall A Defurne et al.  2017
  * CLAS6  Saylor et al.   2018
  * Hall A Georges et al.  2022
  * CLAS12 Lee et al.      2026 (pass-1)

Models used in this stage
-------------------------
  * KM15, evaluated with Gepard.
  * Pure Bethe-Heitler (BH), evaluated through the same Gepard machinery while
    retaining only the elastic BH term.

KM15 is the sole model-assisted transport prescription in this world-data
comparison stage.  PARTONS/GK16 and VGG are intentionally not run here.
GK16 can still be evaluated later for the direct pass-2 model-comparison
section, where no large world-data transport grid is required.

Analysis philosophy
-------------------
1. Keep every experimental point at its native beam energy and native mean
   kinematics for the primary model/data consistency tests.
2. Treat published correlated normalization uncertainties as one nuisance per
   experiment rather than adding them independently to every point.
3. Match different experiments only when their measured kinematics are nearby.
4. For a matched A -> B comparison, transport A to B's exact kinematics with
   the local KM15 ratio

       sigma_A_to_B = sigma_A * sigma_KM15(k_B) / sigma_KM15(k_A).

5. For common-energy presentation, transport each point to Ebeam=10.6 GeV at
   fixed (xB,Q2,t,phi):

       sigma_10p6 = sigma_data
                    * sigma_KM15(10.6) / sigma_KM15(E_native).

6. Quantitative summary tables are retained, but the primary presentation
   products are multi-panel cross-section overlays.  Raw versions show the
   measured cross sections without normalization rescaling; secondary versions
   apply one globally fitted experiment-wide normalization nuisance so the
   effect of the quoted correlated scales can be inspected transparently.

Uncertainty conventions used here
---------------------------------
CLAS6 Jo 2015:
    5% correlated elastic-normalization uncertainty.  It is removed in
    quadrature from the published total point uncertainty by the validated EMFF
    loader and then carried here as a single correlated scale nuisance.

Hall A Defurne 2015:
    Existing validated loader prescription; 2.8% correlated scale.

Hall A Defurne 2017:
    3.2% point-to-point systematic plus correlated
    sqrt(1.0%^2 + 2.0%^2 + 0.5%^2) = 2.291...% scale.

CLAS6 Saylor 2018:
    4% correlated elastic normalization, removed in quadrature from the
    published total systematic by the validated loader.

Hall A Georges 2022:
    Published symmetrized pointwise systematic is retained, and for this
    comparison study a pragmatic 5% correlated normalization prior is allowed.

CLAS12 Lee 2026:
    31% correlated normalization.  The validated authoritative E214M1 loader
    removes it in quadrature from the released total systematic to recover the
    point-to-point component.

Typical use
-----------
From dvcs_cross_section/external_scripts:

  python compare_world_dvcs_cross_sections.py --workers 8

The script is intended to live beside extract_emff_from_dvcs_bh.py in
external_scripts/.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# =============================================================================
# Configuration
# =============================================================================

TARGET_EBEAM_GEV = 10.6

class Pass2UnavailableError(RuntimeError):
    """Pass-2 input is absent or not yet finalized enough for comparison."""
    pass
#endclass


PASS2_OVERALL_NORM_FRAC = 0.021633307652784
PASS2_XS_COL = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol"
PASS2_PTP_COL = "Syst. err (point-to-point total)"
PASS2_CORR_FRAC_COL = "correlated scale sys frac, 10.6 GeV"
PASS2_NORM_FRAC_COL = "uncorrelated normalization sys frac, 10.6 GeV"

PROTON_MASS_GEV = 0.9382720813
PASS2_PERIOD_SPECS = {
    "Sp18 Inb": 10.604,
    "Sp18 Out": 10.604,
    "Fa18 Inb": 10.604,
    "Fa18 Out": 10.604,
    "Sp19 Inb": 10.200,
}


# Georges is intentionally different from the production EMFF extraction:
# for this external-data comparison we allow the ~5% normalization freedom
# indicated by the recent world-data EMFF fits.
GEORGES_COMPARISON_NORM_FRAC = 0.05

# Matching cuts are deliberately explicit and configurable.  They define only
# candidate overlap; the final score is continuous within these windows.
DEFAULT_MATCH_DXB = 0.035
DEFAULT_MATCH_DQ2 = 0.60       # GeV^2
DEFAULT_MATCH_DT = 0.10        # GeV^2
DEFAULT_MATCH_DPHI = 20.0      # degrees

DATASET_ORDER = [
    "jo2015",
    "defurne2015",
    "defurne2017",
    "saylor2018",
    "georges2022",
    "lee2026",
    "pass2",
]

DATASET_LABELS = {
    "jo2015": "CLAS6 Jo 2015",
    "defurne2015": "Hall A Defurne 2015",
    "defurne2017": "Hall A Defurne 2017",
    "saylor2018": "CLAS6 Saylor 2018",
    "georges2022": "Hall A Georges 2022",
    "lee2026": "CLAS12 Lee 2026",
    "pass2": "CLAS12 pass-2 Hayward",
}


# Universal visual identity used everywhere in this script.  These assignments
# are intentionally fixed rather than relying on matplotlib's color cycle, so a
# given experiment/model always has the same appearance from canvas to canvas.
#
# Lee/pass-1 is kept green to remain compatible with the pass-1/pass-2
# comparison style already used in this analysis; red is intentionally reserved
# for the future pass-2 measurement.
DATASET_STYLES = {
    "jo2015":       {"color": "#9467bd", "marker": "o"},  # purple
    "defurne2015":  {"color": "#8c564b", "marker": "^"},  # brown
    "defurne2017":  {"color": "#e377c2", "marker": "v"},  # pink
    "saylor2018":   {"color": "#bcbd22", "marker": "D"},  # olive
    "georges2022":  {"color": "#17becf", "marker": "P"},  # cyan
    "lee2026":      {"color": "#2ca02c", "marker": "s"},  # green
    "pass2":          {"color": "#d62728", "marker": "o"},  # red
}

MODEL_STYLES = {
    "bh":   {"color": "#1f77b4", "linestyle": "-",  "linewidth": 1.35},
    "km15": {"color": "#ff7f0e", "linestyle": "--", "linewidth": 1.55},
}

# Dense model curves for presentation.  The data-model calculations used in
# fits remain evaluated at the exact measured points; this grid is only for
# drawing smooth BH/KM15 curves.  It explicitly includes both 0 and 360 deg.
MODEL_CURVE_PHI_STEP_DEG = 15.0
PANEL_Y_SCALE_MODE = "panel"

# In the global Lee-anchor normalization study Georges is intentionally left
# unconstrained, as requested.  All other experiments receive Gaussian
# multiplicative priors based on their quoted correlated normalization.
GLOBAL_NORM_FREE_DATASETS = {"georges2022"}

GLOBAL_NORM_SCENARIO_LABELS = {
    "nominal": "nominal",
    "saylor_tmin_0p343": r"diagnostic; Saylor $|t|\geq0.343$ GeV$^2$",
    "without_saylor": "diagnostic; without Saylor 2018",
}

# These two Saylor points are treated as invalid and removed from the NOMINAL
# analysis everywhere downstream of the complete source/cache stage.
SAYLOR_NOMINAL_EXCLUDED_POINT_IDS = {
    "saylor2018:2025",  # phi ~52 deg in the current parser
    "saylor2018:2041",  # phi ~308 deg in the current parser
}

# Broader low-|t| diagnostic used in the EMFF world-data study:
# "All six; Saylor |t|>=0.343".
SAYLOR_TMIN_DIAGNOSTIC_GEV2 = 0.343

# Source convention needed for KM15 evaluation and for the existing PARTONS
# phi-mapping logic.
GEPARD_BMK_DATASETS = {"jo2015", "defurne2015", "defurne2017"}
DIRECT_PHI_DATASETS = {"saylor2018", "georges2022", "lee2026", "pass2"}

# evaluate_bh_model_selection.py expects the dataset keys used by the EMFF
# suite.  Keep this translation in exactly one place.
PARTONS_DATASET_KEYS = {
    "jo2015": "jo2015",
    "defurne2015": "halla_defurne2015",
    "defurne2017": "halla_defurne2017",
    "saylor2018": "saylor2018",
    "georges2022": "halla_georges2022",
    "lee2026": "pass1",
}


@dataclass(frozen=True)
class MatchConfig:
    dxb: float = DEFAULT_MATCH_DXB
    dq2: float = DEFAULT_MATCH_DQ2
    dt: float = DEFAULT_MATCH_DT
    dphi: float = DEFAULT_MATCH_DPHI


# =============================================================================
# Small utilities
# =============================================================================


def circular_phi_difference_deg(a, b):
    """Smallest absolute azimuthal separation in degrees, in [0,180]."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.abs((a - b + 180.0) % 360.0 - 180.0)
#enddef


def finite_positive(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return np.isfinite(arr) & (arr > 0.0)
#enddef


def safe_fraction(num, den) -> np.ndarray:
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full(np.broadcast(num, den).shape, np.nan, dtype=float)
    good = np.isfinite(num) & np.isfinite(den) & (den != 0.0)
    out[good] = num[good] / den[good]
    return out
#enddef


def resolve_existing(path: Path, alternatives: Sequence[Path]) -> Path:
    """Resolve one required input while giving useful repository-relative fallbacks."""
    candidates = [Path(path)] + [Path(p) for p in alternatives]
    for candidate in candidates:
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


# =============================================================================
# Canonicalization of the validated EMFF loaders
# =============================================================================


def _base_canonical(df: pd.DataFrame, key: str, norm_frac: float) -> pd.DataFrame:
    out = pd.DataFrame({
        "dataset": key,
        "dataset_label": DATASET_LABELS[key],
        "source_row": np.arange(len(df), dtype=int),
        "xB": pd.to_numeric(df["xB"], errors="coerce"),
        "Q2": pd.to_numeric(df["Q2"], errors="coerce"),
        "t_abs": pd.to_numeric(df["t_abs"], errors="coerce"),
        "phi_deg": np.mod(pd.to_numeric(df["phi_deg"], errors="coerce"), 360.0),
        "ebeam": pd.to_numeric(df["ebeam"], errors="coerce"),
        "xs": pd.to_numeric(df["xs"], errors="coerce"),
        "norm_frac": float(norm_frac),
    })
    return out
#enddef


def canonicalize_jo(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "jo2015", emff.JO_GLOBAL_SCALE_FRAC)

    stat = pd.to_numeric(df.get("clas6_err_stat", np.nan), errors="coerce").to_numpy(float)
    point_total = pd.to_numeric(df["clas6_err_pointwise"], errors="coerce").to_numpy(float)

    # The pointwise total from the validated loader is the published total with
    # the 5% correlated normalization removed.  Recover the residual systematic
    # only for reporting; the fits use point_unc_abs directly.
    stat_for_quad = np.where(np.isfinite(stat) & (stat >= 0.0), stat, 0.0)
    ptp_sys = np.sqrt(np.maximum(point_total**2 - stat_for_quad**2, 0.0))

    out["stat_abs"] = stat
    out["ptp_sys_abs"] = ptp_sys
    out["point_unc_abs"] = point_total
    out["phi_bmk_rad"] = pd.to_numeric(df["phi_bmk_rad"], errors="coerce")
    out["gepard_dataset_id"] = pd.to_numeric(df["gepard_dataset_id"], errors="coerce")
    return clean_canonical(out)
#enddef


def canonicalize_defurne2015(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "defurne2015", emff.HALLA_DEFURNE_GLOBAL_SCALE_FRAC)
    out["stat_abs"] = pd.to_numeric(df["halla_err_stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["halla_err_syst"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    out["phi_bmk_rad"] = pd.to_numeric(df["phi_bmk_rad"], errors="coerce")
    out["gepard_dataset_id"] = pd.to_numeric(df["gepard_dataset_id"], errors="coerce")
    return clean_canonical(out)
#enddef


def canonicalize_defurne2017(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "defurne2017", emff.HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC)
    out["stat_abs"] = pd.to_numeric(df["halla_err_stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["halla_err_syst"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    out["phi_bmk_rad"] = pd.to_numeric(df["phi_bmk_rad"], errors="coerce")
    out["gepard_dataset_id"] = pd.to_numeric(df["gepard_dataset_id"], errors="coerce")
    return clean_canonical(out)
#enddef


def canonicalize_saylor(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "saylor2018", emff.SAYLOR_GLOBAL_SCALE_FRAC)
    out["stat_abs"] = pd.to_numeric(df["stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["ptp_sys"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    return clean_canonical(out)
#enddef


def canonicalize_georges(df: pd.DataFrame) -> pd.DataFrame:
    out = _base_canonical(df, "georges2022", GEORGES_COMPARISON_NORM_FRAC)
    out["stat_abs"] = pd.to_numeric(df["stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["ptp_sys_abs"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    out["setting"] = df["setting"].astype(str).to_numpy()
    out["t_minus_tmin"] = pd.to_numeric(df["t_minus_tmin"], errors="coerce")
    out["t_min"] = pd.to_numeric(df["t_min"], errors="coerce")
    return clean_canonical(out)
#enddef


def canonicalize_lee(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "lee2026", emff.PASS1_GLOBAL_SCALE_FRAC)
    out["stat_abs"] = pd.to_numeric(df["xs_stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["ptp_sys_abs"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    if "bin" in df.columns:
        out["published_bin"] = pd.to_numeric(df["bin"], errors="coerce")
    #endif
    if "released_syst_below_nominal_scale" in df.columns:
        out["released_syst_below_nominal_scale"] = df[
            "released_syst_below_nominal_scale"
        ].astype(bool).to_numpy()
    #endif
    return clean_canonical(out)
#enddef



def enrich_lee_with_legacy_binning(
        lee: pd.DataFrame,
        legacy_path: Path) -> pd.DataFrame:
    """
    Attach the original all_bin_v3 three-dimensional bin identity and phi-bin
    bookkeeping to the authoritative E214M1 Lee points.

    The E214M1 publication renumbered some bins, so published integer bin labels
    cannot be equated directly with the legacy/pass-2 `Bin Name`.  We reproduce
    the validated pass1_published_loader mapping: first map each published
    (xB,Q2,|t|) group to the nearest legacy 3D bin, then identify the nearest
    legacy phi row inside that bin.

    These fields are used ONLY to establish exact same-analysis-bin Lee/Hayward
    comparisons.  Lee central values and uncertainties remain authoritative
    E214M1 values.
    """
    legacy = pd.read_csv(legacy_path, low_memory=False)
    if "valid bin" in legacy.columns:
        valid = pd.to_numeric(legacy["valid bin"], errors="coerce") == 1
        legacy = legacy.loc[valid].copy()
    #endif

    required = [
        "Bin Name", "xBavg", "Q2avg", "t_abs_avg", "phiavg",
        "xBmin", "xBmax", "Q2min", "Q2max",
        "t_abs_min", "t_abs_max", "phimin", "phimax",
    ]
    missing = [c for c in required if c not in legacy.columns]
    if missing:
        raise RuntimeError(
            "Legacy pass-1 binning file is missing: " + ", ".join(missing)
        )
    #endif

    work = lee.copy()
    if "published_bin" not in work.columns:
        raise RuntimeError("Lee canonical table has no published_bin column")
    #endif

    lg = (
        legacy.groupby("Bin Name", sort=False)
        .agg(
            gx=("xBavg", "mean"),
            gq=("Q2avg", "mean"),
            gt=("t_abs_avg", "mean"),
        )
        .reset_index()
    )
    pg = (
        work.groupby("published_bin", sort=False)
        .agg(
            px=("xB", "mean"),
            pq=("Q2", "mean"),
            pt=("t_abs", "mean"),
        )
        .reset_index()
    )

    bin_map = {}
    for r in pg.itertuples(index=False):
        dist = np.sqrt(
            ((lg["gx"] - float(r.px)) / 0.03)**2
            + ((lg["gq"] - float(r.pq)) / 0.30)**2
            + ((lg["gt"] - float(r.pt)) / 0.10)**2
        )
        ii = int(np.nanargmin(dist.to_numpy(float)))
        bin_map[int(round(float(r.published_bin)))] = lg.iloc[ii]["Bin Name"]
    #endfor

    attach_rows = []
    unmatched = 0
    for r in work.itertuples(index=False):
        legacy_bin = bin_map.get(int(round(float(r.published_bin))))
        cand = legacy.loc[legacy["Bin Name"] == legacy_bin].copy()
        if cand.empty:
            unmatched += 1
            attach_rows.append({})
            continue
        #endif
        dphi = np.abs(
            pd.to_numeric(cand["phiavg"], errors="coerce").to_numpy(float)
            - float(r.phi_deg)
        )
        if not np.any(np.isfinite(dphi)):
            unmatched += 1
            attach_rows.append({})
            continue
        #endif
        row = cand.iloc[int(np.nanargmin(dphi))]
        attach_rows.append({
            "analysis_bin_name": float(row["Bin Name"]),
            "legacy_phiavg": float(row["phiavg"]),
            "legacy_xBavg": float(row["xBavg"]),
            "legacy_Q2avg": float(row["Q2avg"]),
            "legacy_t_abs_avg": float(row["t_abs_avg"]),
            "xBmin": float(row["xBmin"]),
            "xBmax": float(row["xBmax"]),
            "Q2min": float(row["Q2min"]),
            "Q2max": float(row["Q2max"]),
            "t_abs_min": float(row["t_abs_min"]),
            "t_abs_max": float(row["t_abs_max"]),
            "phimin": float(row["phimin"]),
            "phimax": float(row["phimax"]),
        })
    #endfor

    attached = pd.DataFrame(attach_rows, index=work.index)
    for c in attached.columns:
        work[c] = attached[c]
    #endfor

    print(
        f"[LEE EXACT BINNING] mapped {len(work)-unmatched:,}/{len(work):,} "
        f"E214M1 points to original all_bin_v3 analysis bins using {legacy_path}",
        flush=True,
    )
    if unmatched:
        warnings.warn(f"{unmatched} Lee points could not be mapped to legacy binning")
    #endif
    return work
#enddef


def _parse_tuple3_loose(raw) -> Tuple[float, float, float]:
    """Parse a CSV tuple '(value, stat, sys)' with NaN fallback."""
    if raw is None or (isinstance(raw, float) and not np.isfinite(raw)):
        return np.nan, np.nan, np.nan
    #endif
    txt = str(raw).strip()
    if not txt:
        return np.nan, np.nan, np.nan
    #endif
    txt = txt.strip("()")
    parts = [x.strip() for x in txt.split(",")]
    try:
        vals = [float(x) for x in parts[:3]]
    except Exception:
        return np.nan, np.nan, np.nan
    #endtry
    while len(vals) < 3:
        vals.append(0.0)
    #endwhile
    return vals[0], vals[1], vals[2]
#enddef


def _parse_cross_section_tuple(raw) -> Tuple[float, float]:
    """
    Parse the pass-2 CSV tuple '(value, stat, ...)' and return (value, stat).
    """
    if raw is None or (isinstance(raw, float) and not np.isfinite(raw)):
        return np.nan, np.nan
    #endif
    s = str(raw).strip()
    if not (s.startswith("(") and s.endswith(")")):
        return np.nan, np.nan
    #endif
    fields = [x.strip() for x in s[1:-1].split(",")]
    if len(fields) < 2:
        return np.nan, np.nan
    #endif
    try:
        return float(fields[0]), float(fields[1])
    except Exception:
        return np.nan, np.nan
    #endtry
#enddef


def canonicalize_pass2_csv(path: Path) -> pd.DataFrame:
    """
    Load the finalized combined 10.6-GeV pass-2 cross section using the
    publication-level systematic decomposition:

      1. point-to-point uncertainty;
      2. one experiment-wide overall normalization nuisance;
      3. one bin-dependent correlated-scale nuisance.

    The correlated-scale nuisance is read directly from the authoritative final
    CSV column

        correlated scale sys frac, 10.6 GeV

    which the finalized C++ systematic chain constructs from the lower-level
    current-dependence and residual run-period studies.  Those lower-level
    components are intentionally NOT reopened as separate nuisance parameters
    here, because the publication treatment uses the single combined category.
    """
    raw = pd.read_csv(path, low_memory=False)

    required = [
        "bin index",
        "Bin Name",
        "xBavg, 10.6 GeV",
        "Q2avg, 10.6 GeV",
        "t_abs_avg, 10.6 GeV",
        "phiavg, 10.6 GeV",
        PASS2_XS_COL,
        PASS2_PTP_COL,
        PASS2_CORR_FRAC_COL,
        PASS2_NORM_FRAC_COL,
    ]
    missing = [c for c in required if c not in raw.columns]
    if missing:
        raise Pass2UnavailableError(
            "Pass-2 CSV exists but is not yet finalized for the world-data "
            "comparison. Missing authoritative final-systematics column(s): "
            + "; ".join(missing)
            + ". The published-world-data analysis will continue without "
              "CLAS12 pass-2 Hayward."
        )
    #endif

    xs_vals = []
    stat_vals = []
    for v in raw[PASS2_XS_COL]:
        xs, stat = _parse_cross_section_tuple(v)
        xs_vals.append(xs)
        stat_vals.append(stat)
    #endfor

    out = pd.DataFrame({
        "dataset": "pass2",
        "dataset_label": DATASET_LABELS["pass2"],
        "source_row": np.arange(len(raw), dtype=int),

        # IMPORTANT:
        #   Bin Name  = the 3D (xB,Q2,t) bin shared by the phi distribution.
        #   bin index = the unique 4D (xB,Q2,t,phi) point identifier.
        #
        # Presentation canvases must group by Bin Name, while matching and
        # nuisance fits need a unique point_id built from bin index.
        "published_bin": pd.to_numeric(raw["Bin Name"], errors="coerce"),
        "four_d_bin_index": pd.to_numeric(raw["bin index"], errors="coerce"),

        "xB": pd.to_numeric(raw["xBavg, 10.6 GeV"], errors="coerce"),
        "Q2": pd.to_numeric(raw["Q2avg, 10.6 GeV"], errors="coerce"),
        "t_abs": pd.to_numeric(raw["t_abs_avg, 10.6 GeV"], errors="coerce"),
        "phi_deg": np.mod(
            pd.to_numeric(raw["phiavg, 10.6 GeV"], errors="coerce"),
            360.0,
        ),
        "ebeam": float(TARGET_EBEAM_GEV),
        "xs": np.asarray(xs_vals, dtype=float),
        "stat_abs": np.asarray(stat_vals, dtype=float),
        "ptp_sys_abs": pd.to_numeric(raw[PASS2_PTP_COL], errors="coerce"),
        "norm_frac": pd.to_numeric(raw[PASS2_NORM_FRAC_COL], errors="coerce"),
        "corr_scale_frac": pd.to_numeric(raw[PASS2_CORR_FRAC_COL], errors="coerce"),
    })

    for out_col, csv_col in [
        ("e_theta", "e_theta, 10.6 GeV"),
        ("p_theta", "p_theta, 10.6 GeV"),
        ("g_theta", "g_theta, 10.6 GeV"),
    ]:
        if csv_col in raw.columns:
            out[out_col] = pd.to_numeric(raw[csv_col], errors="coerce")
        #endif
    #endfor

    for edge_col in [
        "xBmin", "xBmax", "Q2min", "Q2max",
        "t_abs_min", "t_abs_max", "phimin", "phimax",
    ]:
        if edge_col in raw.columns:
            out[edge_col] = pd.to_numeric(raw[edge_col], errors="coerce")
        #endif
    #endfor

    # For an exclusive recoil proton on a target at rest, |t| fixes the recoil
    # proton energy and momentum.  This is useful for detector-correlation
    # diagnostics even when the CSV does not store an explicit mean p_p.
    ep = PROTON_MASS_GEV + out["t_abs"].to_numpy(float) / (2.0 * PROTON_MASS_GEV)
    out["p_p"] = np.sqrt(np.maximum(ep*ep - PROTON_MASS_GEV**2, 0.0))

    # Ordinary plotted/fitted point uncertainty contains ONLY independent
    # statistical + point-to-point systematic terms.
    out["point_unc_abs"] = np.hypot(
        out["stat_abs"],
        out["ptp_sys_abs"],
    )

    # The overall normalization is one experiment-wide publication quantity,
    # not a per-bin measurement.  If an otherwise physical row has that CSV
    # field blank, use the finalized global value rather than dropping the row.
    n_missing_norm = int(out["norm_frac"].isna().sum())
    if n_missing_norm:
        print(
            f"[PASS2 HAYWARD] filling {n_missing_norm} blank per-row normalization field(s) "
            f"with finalized global value {100.0*PASS2_OVERALL_NORM_FRAC:.3f}%",
            flush=True,
        )
        out["norm_frac"] = out["norm_frac"].fillna(PASS2_OVERALL_NORM_FRAC)
    #endif

    # Do NOT fill missing correlated-scale responses.  Such rows remain valid
    # for raw/norm-only comparisons but cannot enter a fit using beta_corr.
    n_missing_corr_before_clean = int(out["corr_scale_frac"].isna().sum())
    if n_missing_corr_before_clean:
        bad_rows = out.loc[
            out["corr_scale_frac"].isna()
            & np.isfinite(out["xs"])
            & (out["xs"] > 0.0),
            ["published_bin", "xB", "Q2", "t_abs", "phi_deg", "xs"],
        ]
        print(
            f"[PASS2 HAYWARD] WARNING: {len(bad_rows)} physical point(s) have no finalized "
            "correlated-scale response; retained for raw/norm-only comparisons "
            "and excluded only from full correlated-scale nuisance fits.",
            flush=True,
        )
        if not bad_rows.empty:
            print(
                bad_rows.to_string(
                    index=False,
                    float_format=lambda x: f"{x:.7g}",
                ),
                flush=True,
            )
        #endif
    #endif

    # Publication-level QA: the normalization column should reproduce the
    # finalized 2.1633% common charge/target-thickness uncertainty.
    finite_norm = out["norm_frac"].to_numpy(float)
    finite_norm = finite_norm[np.isfinite(finite_norm)]
    if finite_norm.size == 0:
        raise RuntimeError(
            "Pass-2 overall-normalization column contains no finite values"
        )
    #endif
    med_norm = float(np.nanmedian(finite_norm))
    if abs(med_norm - PASS2_OVERALL_NORM_FRAC) > 5e-4:
        warnings.warn(
            f"Pass-2 median overall normalization is "
            f"{100.0 * med_norm:.3f}% instead of expected "
            f"{100.0 * PASS2_OVERALL_NORM_FRAC:.3f}%"
        )
    #endif

    out = clean_canonical(out)

    out["point_id"] = [
        f"pass2:{int(round(v))}"
        for v in out["four_d_bin_index"].to_numpy(float)
    ]

    n_corr_valid = int(np.isfinite(out["corr_scale_frac"]).sum())
    n_3d_bins = int(out["published_bin"].nunique())
    print(
        f"[PASS2 HAYWARD] loaded {len(out):,} physical 10.6-GeV points "
        f"across {n_3d_bins} three-dimensional (xB,Q2,t) bins from {path}; "
        f"{n_corr_valid:,} points have finalized correlated-scale responses",
        flush=True,
    )
    print(
        f"[PASS2 HAYWARD] median stat={100*np.nanmedian(out['stat_frac']):.2f}%, "
        f"ptp syst={100*np.nanmedian(out['ptp_sys_frac']):.2f}%, "
        f"overall norm={100*np.nanmedian(out['norm_frac']):.2f}%, "
        f"correlated scale={100*np.nanmedian(out['corr_scale_frac']):.2f}%",
        flush=True,
    )
    return out
#enddef

def clean_canonical(df: pd.DataFrame) -> pd.DataFrame:
    required = [
        "xB", "Q2", "t_abs", "phi_deg", "ebeam", "xs",
        "stat_abs", "ptp_sys_abs", "point_unc_abs", "norm_frac",
    ]
    finite = np.ones(len(df), dtype=bool)
    for col in required:
        finite &= np.isfinite(pd.to_numeric(df[col], errors="coerce").to_numpy(float))
    #endfor
    finite &= df["xs"].to_numpy(float) > 0.0
    finite &= df["point_unc_abs"].to_numpy(float) > 0.0
    finite &= df["Q2"].to_numpy(float) > 0.0
    finite &= df["t_abs"].to_numpy(float) > 0.0

    nbad = int((~finite).sum())
    if nbad:
        warnings.warn(f"Dropping {nbad} nonfinite/nonphysical canonical rows")
    #endif

    out = df.loc[finite].copy().reset_index(drop=True)
    out["source_row"] = np.arange(len(out), dtype=int)
    out["point_id"] = [
        f"{dataset}:{row}"
        for dataset, row in zip(out["dataset"], out["source_row"])
    ]
    out["stat_frac"] = out["stat_abs"] / out["xs"]
    out["ptp_sys_frac"] = out["ptp_sys_abs"] / out["xs"]
    out["point_unc_frac"] = out["point_unc_abs"] / out["xs"]

    if "corr_scale_frac" not in out.columns:
        # Published external datasets have no pass-2-style kinematic correlated
        # scale category, so their response is identically zero.
        out["corr_scale_frac"] = 0.0
    else:
        # If a dataset explicitly supplies this category (pass-2), preserve
        # missing values as missing.  They must never be silently interpreted
        # as zero response.
        out["corr_scale_frac"] = pd.to_numeric(
            out["corr_scale_frac"],
            errors="coerce",
        )
    #endif

    return out
#enddef


def load_world_data(args, emff) -> pd.DataFrame:
    """Load the six published datasets plus the finalized pass-2 measurement."""
    print("\n" + "=" * 80)
    print("LOADING PUBLISHED WORLD DATA + OPTIONAL PASS-2")
    print("=" * 80)

    jo = canonicalize_jo(emff.load_clas6_gepard_dataset(), emff)
    d15 = canonicalize_defurne2015(emff.load_halla_defurne_gepard_datasets(), emff)
    d17 = canonicalize_defurne2017(emff.load_halla_defurne2017_gepard_datasets(), emff)

    saylor_path = resolve_existing(
        Path(args.saylor_file),
        [args.script_dir / "import" / "saylor_CLAS6.txt"],
    )
    saylor = canonicalize_saylor(emff.load_saylor_supplement(str(saylor_path)), emff)

    georges_path = resolve_existing(
        Path(args.georges_file),
        [args.script_dir / "import" / "E12-06-114.xlsx"],
    )
    georges = canonicalize_georges(emff.load_georges_supplement(str(georges_path)))

    lee_path = resolve_existing(
        Path(args.lee_file),
        [
            args.script_dir / "import" / "clasdb_E214M1.txt",
            args.script_dir.parent / "imports" / "clasdb_E214M1.txt",
        ],
    )
    lee = canonicalize_lee(emff.load_clas12_pass1_csv(lee_path), emff)
    legacy_path = resolve_existing(
        Path(args.pass1_legacy_file),
        [
            args.script_dir / "import" / "all_bin_v3.csv",
            args.script_dir.parent / "imports" / "all_bin_v3.csv",
        ],
    )
    lee = enrich_lee_with_legacy_binning(lee, legacy_path)
    args.resolved_pass1_legacy_file = legacy_path

    pass2_path = args.pass2_file
    pass2 = pd.DataFrame()
    if pass2_path is not None:
        pass2_path = Path(pass2_path)
        if not pass2_path.exists():
            print(
                f"[PASS2 HAYWARD] input not found: {pass2_path}",
                flush=True,
            )
            print(
                "[PASS2 HAYWARD] skipping pass-2; continuing with "
                "published world data only",
                flush=True,
            )
        else:
            try:
                pass2 = canonicalize_pass2_csv(pass2_path)
            except Pass2UnavailableError as exc:
                print(f"[PASS2 HAYWARD] {exc}", flush=True)
                print(
                    "[PASS2 HAYWARD] skipping pass-2; continuing with "
                    "published world data only",
                    flush=True,
                )
            #endtry
        #endif
    #endif
    args.resolved_pass2_file = pass2_path

    world = pd.concat(
        [jo, d15, d17, saylor, georges, lee, pass2],
        ignore_index=True,
        sort=False,
    )
    world["dataset"] = pd.Categorical(world["dataset"], DATASET_ORDER, ordered=True)
    world = world.sort_values(["dataset", "source_row"]).reset_index(drop=True)
    world["dataset"] = world["dataset"].astype(str)

    if world["point_id"].duplicated().any():
        raise RuntimeError("Canonical point_id collision detected")
    #endif

    print(f"[WORLD] canonicalized {len(world):,} points across 7 measurements")
    return world
#enddef



def apply_nominal_data_quality_exclusions(world: pd.DataFrame) -> pd.DataFrame:
    """Remove the two known-invalid Saylor bin-87 points from the nominal sample."""
    mask = world["point_id"].astype(str).isin(SAYLOR_NOMINAL_EXCLUDED_POINT_IDS)
    removed = world.loc[mask].copy()

    print(
        f"[NOMINAL QUALITY] excluding {len(removed)} invalid Saylor bin-87 point(s)",
        flush=True,
    )
    if not removed.empty:
        cols = [
            c for c in [
                "point_id", "xB", "Q2", "t_abs", "phi_deg",
                "xs", "stat_abs", "ptp_sys_abs",
            ]
            if c in removed.columns
        ]
        print(
            removed[cols].to_string(
                index=False,
                float_format=lambda x: f"{x:.7g}",
            ),
            flush=True,
        )
    #endif

    found = set(removed["point_id"].astype(str))
    missing = set(SAYLOR_NOMINAL_EXCLUDED_POINT_IDS) - found
    if missing:
        warnings.warn(
            "Expected nominal Saylor exclusion point(s) not found: "
            + ", ".join(sorted(missing))
        )
    #endif

    return world.loc[~mask].reset_index(drop=True)
#enddef


# =============================================================================
# KM15/BH model evaluation at native and common-energy kinematics
# =============================================================================


def evaluate_one_km15(emff, row, ebeam: float) -> Dict[str, float]:
    """Use the already validated dataset-specific phi convention."""
    key = str(row.dataset)
    if key in GEPARD_BMK_DATASETS:
        task = (
            0,
            float(row.xB),
            float(row.Q2),
            float(row.t_abs),
            float(row.phi_bmk_rad),
            float(ebeam),
        )
        return emff.evaluate_km15_point_bmk(task)
    #endif

    if key in DIRECT_PHI_DATASETS:
        task = (
            0,
            float(row.xB),
            float(row.Q2),
            float(row.t_abs),
            float(row.phi_deg),
            float(ebeam),
            "identity",
        )
        return emff.evaluate_km15_point(task)
    #endif

    raise KeyError(f"No KM15 phi convention configured for dataset {key}")
#enddef



def _display_phi_to_bmk_rad(phi_deg: float) -> float:
    """Map display phi in [0,360] to the equivalent BMK angle in [-pi,pi]."""
    wrapped = (float(phi_deg) + 180.0) % 360.0 - 180.0
    return math.radians(wrapped)
#enddef


def evaluate_dense_model_curve(
        emff,
        dataset_key: str,
        ebeam: float,
        xB: float,
        Q2: float,
        t_abs: float,
        phi_step_deg: float = MODEL_CURVE_PHI_STEP_DEG) -> pd.DataFrame:
    """
    Evaluate smooth BH and KM15 curves at one fixed hadronic kinematic point.

    This deliberately does NOT connect the model values already attached to
    neighboring data points.  Those data points can have slightly different
    mean xB, Q2 and |t| values, so connecting them directly can create an
    artificial kink or apparently "wrong" model curve.  Here xB, Q2, |t| and
    Ebeam are fixed and only phi is scanned.

    The grid always contains phi=0 and phi=360 exactly.
    """
    step = float(phi_step_deg)
    if not np.isfinite(step) or step <= 0.0 or step > 90.0:
        raise ValueError(f"Invalid dense-model phi step: {phi_step_deg}")
    #endif

    nstep = int(math.ceil(360.0 / step))
    phi_grid = np.linspace(0.0, 360.0, nstep + 1)

    rows = []
    for phi in phi_grid:
        if dataset_key in GEPARD_BMK_DATASETS:
            task = (
                0,
                float(xB),
                float(Q2),
                float(t_abs),
                _display_phi_to_bmk_rad(float(phi)),
                float(ebeam),
            )
            result = emff.evaluate_km15_point_bmk(task)
        elif dataset_key in DIRECT_PHI_DATASETS:
            task = (
                0,
                float(xB),
                float(Q2),
                float(t_abs),
                float(phi),
                float(ebeam),
                "identity",
            )
            result = emff.evaluate_km15_point(task)
        else:
            raise KeyError(f"No KM15 phi convention configured for dataset {dataset_key}")
        #endif

        rows.append({
            "phi_deg": float(phi),
            "km15": float(result["km15_ep"]),
            "bh": float(result["km15_bh"]),
        })
    #endfor

    return pd.DataFrame(rows)
#enddef


def get_dense_model_curve(
        emff,
        cache: Dict[Tuple, pd.DataFrame],
        dataset_key: str,
        ebeam: float,
        xB: float,
        Q2: float,
        t_abs: float) -> pd.DataFrame:
    """
    Return a cached dense BH/KM15 phi scan.

    Rounded kinematics are used only for the in-memory cache key, not for the
    calculation itself.  This prevents repeated evaluation when the same Lee
    anchor appears in several presentation products.
    """
    key = (
        str(dataset_key),
        round(float(ebeam), 6),
        round(float(xB), 6),
        round(float(Q2), 6),
        round(float(t_abs), 6),
        round(float(MODEL_CURVE_PHI_STEP_DEG), 6),
    )
    if key not in cache:
        cache[key] = evaluate_dense_model_curve(
            emff,
            dataset_key=str(dataset_key),
            ebeam=float(ebeam),
            xB=float(xB),
            Q2=float(Q2),
            t_abs=float(t_abs),
            phi_step_deg=float(MODEL_CURVE_PHI_STEP_DEG),
        )
    #endif
    return cache[key]
#enddef


def evaluate_km15_world(
        world: pd.DataFrame,
        emff,
        cache_path: Path,
        target_ebeam: float,
        force: bool = False) -> pd.DataFrame:
    """
    Evaluate KM15 total and pure BH at native E and target E.

    The cache is incremental: existing external-world entries are reused and
    only missing point IDs (normally the new pass-2 points) are evaluated.
    """
    expected_cols = [
        "point_id",
        "km15_native", "bh_native",
        "km15_target", "bh_target",
    ]

    cache = pd.DataFrame(columns=expected_cols)
    if cache_path.exists() and not force:
        try:
            old = pd.read_csv(cache_path)
            if set(expected_cols).issubset(old.columns):
                cache = old[expected_cols].copy()
                cache["point_id"] = cache["point_id"].astype(str)
                cache = cache.drop_duplicates("point_id", keep="last")
                print(
                    f"[KM15] loaded incremental cache with {len(cache):,} points: "
                    f"{cache_path}",
                    flush=True,
                )
            #endif
        except Exception as exc:
            warnings.warn(f"Could not reuse KM15 cache: {exc}")
        #endtry
    #endif

    cache_map = set(cache["point_id"].astype(str))
    missing = world.loc[~world["point_id"].astype(str).isin(cache_map)].copy()

    if force:
        cache = pd.DataFrame(columns=expected_cols)
        missing = world.copy()
    #endif

    if not missing.empty:
        rows = []
        total = len(missing)
        print(
            f"[KM15] evaluating {total:,} missing point(s) "
            f"(native + E={target_ebeam:.3f} GeV)",
            flush=True,
        )
        for i, row in enumerate(missing.itertuples(index=False), start=1):
            native = evaluate_one_km15(emff, row, float(row.ebeam))
            target = evaluate_one_km15(emff, row, float(target_ebeam))
            rows.append({
                "point_id": str(row.point_id),
                "km15_native": float(native["km15_ep"]),
                "bh_native": float(native["km15_bh"]),
                "km15_target": float(target["km15_ep"]),
                "bh_target": float(target["km15_bh"]),
            })
            if i % 100 == 0 or i == total:
                print(
                    f"[KM15] missing-point evaluation {i}/{total} "
                    f"({100.0*i/total:5.1f}%)",
                    flush=True,
                )
            #endif
        #endfor
        cache = pd.concat([cache, pd.DataFrame(rows)], ignore_index=True)
        cache = cache.drop_duplicates("point_id", keep="last")
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache.to_csv(cache_path, index=False)
        print(f"[KM15] updated cache -> {cache_path}", flush=True)
    else:
        print("[KM15] all requested points found in cache", flush=True)
    #endif

    merged = world.merge(cache, on="point_id", how="left", validate="one_to_one")
    if merged[expected_cols[1:]].isna().any(axis=None):
        bad = merged.loc[merged[expected_cols[1:]].isna().any(axis=1), "point_id"].head(10)
        raise RuntimeError(
            "KM15 cache merge left missing predictions for: "
            + ", ".join(bad.astype(str))
        )
    #endif

    return finalize_km15_columns(merged)
#enddef

def finalize_km15_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["data_over_km15"] = out["xs"] / out["km15_native"]
    out["data_over_bh"] = out["xs"] / out["bh_native"]
    out["bh_fraction_km15"] = out["bh_native"] / out["km15_native"]
    out["km15_transport_factor"] = out["km15_target"] / out["km15_native"]
    out["bh_transport_factor"] = out["bh_target"] / out["bh_native"]
    out["xs_10p6_km15"] = out["xs"] * out["km15_transport_factor"]
    out["stat_10p6_km15"] = out["stat_abs"] * np.abs(out["km15_transport_factor"])
    out["ptp_sys_10p6_km15"] = out["ptp_sys_abs"] * np.abs(out["km15_transport_factor"])
    out["point_unc_10p6_km15"] = out["point_unc_abs"] * np.abs(out["km15_transport_factor"])
    return out
#enddef


# =============================================================================
# Dataset summaries and normalization-aware model scores
# =============================================================================


def make_dataset_summary(world: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key in DATASET_ORDER:
        d = world.loc[world["dataset"] == key]
        if d.empty:
            continue
        #endif
        rows.append({
            "dataset": key,
            "dataset_label": DATASET_LABELS[key],
            "N": len(d),
            "ebeam_min_GeV": d["ebeam"].min(),
            "ebeam_max_GeV": d["ebeam"].max(),
            "xB_min": d["xB"].min(),
            "xB_max": d["xB"].max(),
            "Q2_min_GeV2": d["Q2"].min(),
            "Q2_max_GeV2": d["Q2"].max(),
            "t_abs_min_GeV2": d["t_abs"].min(),
            "t_abs_max_GeV2": d["t_abs"].max(),
            "median_stat_pct": 100.0 * np.nanmedian(d["stat_frac"]),
            "median_ptp_sys_pct": 100.0 * np.nanmedian(d["ptp_sys_frac"]),
            "median_point_unc_pct": 100.0 * np.nanmedian(d["point_unc_frac"]),
            "correlated_norm_pct": 100.0 * float(d["norm_frac"].iloc[0]),
        })
    #endfor
    return pd.DataFrame(rows)
#enddef


def fit_one_normalization_nuisance(
        data: np.ndarray,
        model: np.ndarray,
        sigma: np.ndarray,
        norm_frac: float) -> Dict[str, float]:
    """
    Fit one correlated multiplicative scale nuisance beta.

    Convention:
        shifted_data = data * (1 + beta * norm_frac)

    and minimize
        sum[(shifted_data-model)^2/sigma^2] + beta^2.

    beta therefore has units of published normalization standard deviations.
    """
    y = np.asarray(data, dtype=float)
    m = np.asarray(model, dtype=float)
    s = np.asarray(sigma, dtype=float)
    good = np.isfinite(y) & np.isfinite(m) & np.isfinite(s) & (s > 0.0) & (y > 0.0) & (m > 0.0)
    y, m, s = y[good], m[good], s[good]

    if len(y) == 0:
        return {"N": 0, "beta": np.nan, "scale": np.nan, "chi2": np.nan, "chi2_per_point": np.nan}
    #endif

    if norm_frac <= 0.0:
        beta = 0.0
    else:
        x = norm_frac * y
        w = 1.0 / (s * s)
        beta = -np.sum(w * x * (y - m)) / (1.0 + np.sum(w * x * x))
    #endif

    scale = 1.0 + beta * norm_frac
    residual = scale * y - m
    chi2_points = np.sum((residual / s)**2)
    chi2 = chi2_points + beta**2

    return {
        "N": len(y),
        "beta": float(beta),
        "scale": float(scale),
        "chi2": float(chi2),
        "chi2_per_point": float(chi2 / len(y)),
        "median_abs_fractional_residual": float(np.nanmedian(np.abs(residual / m))),
        "rms_pull": float(np.sqrt(np.nanmean((residual / s)**2))),
    }
#enddef



def fit_pass2_model_publication_nuisances(
        data: np.ndarray,
        model: np.ndarray,
        point_unc: np.ndarray,
        norm_frac: np.ndarray,
        corr_scale_frac: np.ndarray,
        *,
        include_norm: bool,
        include_corr: bool) -> Dict[str, object]:
    """
    Compare CLAS12 pass-2 Hayward directly to one fixed model using the
    publication-level systematic decomposition.

    The nominal point residual is

        r_i = data_i - model_i

    and the two correlated response vectors are

        u_norm,i = n * data_i
        u_corr,i = f_corr,i * data_i .

    We minimize

        chi2 =
          sum_i [
            (data_i + beta_norm*u_norm,i + beta_corr*u_corr,i - model_i)
            / sigma_point,i
          ]^2
          + beta_norm^2
          + beta_corr^2,

    including only the nuisance terms enabled for the requested scenario.

    This is the standard linear correlated-systematic nuisance formulation:
    the response vectors are fixed at the nominal measured cross section, the
    point-to-point uncertainty remains in the diagonal denominator, and each
    correlated nuisance has a unit Gaussian prior.

    beta_norm and beta_corr are therefore directly in units of the assigned
    one-sigma systematic uncertainties.
    """
    y = np.asarray(data, dtype=float)
    m = np.asarray(model, dtype=float)
    s = np.asarray(point_unc, dtype=float)
    n = np.asarray(norm_frac, dtype=float)
    c = np.asarray(corr_scale_frac, dtype=float)

    good = (
        np.isfinite(y)
        & np.isfinite(m)
        & np.isfinite(s)
        & np.isfinite(n)
        & (y > 0.0)
        & (m > 0.0)
        & (s > 0.0)
    )
    if include_corr:
        good &= np.isfinite(c)
    #endif

    y = y[good]
    m = m[good]
    s = s[good]
    n = n[good]
    c = c[good]

    if len(y) == 0:
        return {
            "N": 0,
            "beta_norm": np.nan,
            "beta_corr": np.nan,
            "chi2_data": np.nan,
            "chi2_prior": np.nan,
            "chi2_total": np.nan,
            "chi2_per_point": np.nan,
            "raw_rms_pull": np.nan,
            "fitted_rms_pull": np.nan,
            "point_table": pd.DataFrame(),
        }
    #endif

    residual0 = y - m
    raw_pull = residual0 / s

    response_columns = []
    response_names = []

    if include_norm:
        response_columns.append(n * y)
        response_names.append("norm")
    #endif

    if include_corr:
        response_columns.append(c * y)
        response_names.append("corr")
    #endif

    beta_norm = 0.0
    beta_corr = 0.0

    if response_columns:
        X = np.column_stack(response_columns)
        w = 1.0 / (s * s)

        # (X^T W X + I) beta = -X^T W r0
        lhs = X.T @ (w[:, None] * X) + np.eye(X.shape[1])
        rhs = -(X.T @ (w * residual0))
        beta = np.linalg.solve(lhs, rhs)

        for name, value in zip(response_names, beta):
            if name == "norm":
                beta_norm = float(value)
            elif name == "corr":
                beta_corr = float(value)
            #endif
        #endfor
    #endif

    norm_shift_frac = beta_norm * n if include_norm else np.zeros_like(y)
    corr_shift_frac = beta_corr * c if include_corr else np.zeros_like(y)

    shifted = y * (1.0 + norm_shift_frac + corr_shift_frac)
    fitted_residual = shifted - m
    fitted_pull = fitted_residual / s

    chi2_data = float(np.sum(fitted_pull**2))
    chi2_prior = (
        (beta_norm**2 if include_norm else 0.0)
        + (beta_corr**2 if include_corr else 0.0)
    )
    chi2_total = chi2_data + chi2_prior

    # Correlation between the two response directions in the weighted data
    # space.  Values near +/-1 mean the model cannot cleanly distinguish an
    # overall normalization movement from the assigned kinematic scale shape.
    response_correlation = np.nan
    if include_norm and include_corr:
        un = n * y / s
        uc = c * y / s
        denom = math.sqrt(float(np.sum(un**2) * np.sum(uc**2)))
        if denom > 0.0:
            response_correlation = float(np.sum(un * uc) / denom)
        #endif
    #endif

    point_table = pd.DataFrame({
        "data": y,
        "model": m,
        "point_unc": s,
        "norm_frac": n,
        "corr_scale_frac": c,
        "raw_residual": residual0,
        "raw_pull": raw_pull,
        "norm_shift_frac": norm_shift_frac,
        "corr_shift_frac": corr_shift_frac,
        "total_correlated_shift_frac": norm_shift_frac + corr_shift_frac,
        "shifted_data": shifted,
        "fitted_residual": fitted_residual,
        "fitted_pull": fitted_pull,
    })

    corr_pct = 100.0 * corr_shift_frac
    total_pct = 100.0 * (norm_shift_frac + corr_shift_frac)

    return {
        "N": int(len(y)),
        "beta_norm": float(beta_norm),
        "beta_corr": float(beta_corr),
        "normalization_shift_pct": (
            100.0 * float(np.nanmedian(norm_shift_frac))
            if include_norm else 0.0
        ),
        "corr_shift_median_pct": (
            float(np.nanmedian(corr_pct))
            if include_corr else 0.0
        ),
        "corr_shift_min_pct": (
            float(np.nanmin(corr_pct))
            if include_corr else 0.0
        ),
        "corr_shift_max_pct": (
            float(np.nanmax(corr_pct))
            if include_corr else 0.0
        ),
        "total_shift_median_pct": float(np.nanmedian(total_pct)),
        "total_shift_min_pct": float(np.nanmin(total_pct)),
        "total_shift_max_pct": float(np.nanmax(total_pct)),
        "response_correlation_norm_corr": response_correlation,
        "combined_nuisance_excursion_sigma": float(
            math.sqrt(
                (beta_norm**2 if include_norm else 0.0)
                + (beta_corr**2 if include_corr else 0.0)
            )
        ),
        "chi2_data": chi2_data,
        "chi2_prior": float(chi2_prior),
        "chi2_total": float(chi2_total),
        "chi2_per_point": float(chi2_total / len(y)),
        "raw_rms_pull": float(np.sqrt(np.mean(raw_pull**2))),
        "fitted_rms_pull": float(np.sqrt(np.mean(fitted_pull**2))),
        "median_abs_fractional_residual_raw": float(
            np.nanmedian(np.abs(residual0 / m))
        ),
        "median_abs_fractional_residual_fitted": float(
            np.nanmedian(np.abs(fitted_residual / m))
        ),
        "point_table": point_table,
    }
#enddef


def make_pass2_model_publication_scores(
        world: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Produce direct Hayward-vs-KM15/BH diagnostics for three nested uncertainty
    treatments:

      raw:
        stat + point-to-point only;

      norm_only:
        raw + one 2.16% overall normalization nuisance;

      norm_plus_corr:
        raw + one 2.16% normalization nuisance + one finalized bin-dependent
        correlated-scale nuisance.

    This is the model-comparison table that should be used for pass-2 rather
    than the older one-normalization-only native-model score.
    """
    d = world.loc[world["dataset"] == "pass2"].copy()
    if d.empty:
        return pd.DataFrame(), pd.DataFrame()
    #endif

    model_specs = [
        ("KM15", "km15_native"),
        ("BH", "bh_native"),
    ]
    scenarios = [
        ("raw", False, False),
        ("norm_only", True, False),
        ("norm_plus_corr", True, True),
    ]

    summary_rows = []
    point_rows = []

    for model_name, model_col in model_specs:
        for scenario, include_norm, include_corr in scenarios:
            result = fit_pass2_model_publication_nuisances(
                d["xs"].to_numpy(float),
                d[model_col].to_numpy(float),
                d["point_unc_abs"].to_numpy(float),
                d["norm_frac"].to_numpy(float),
                d["corr_scale_frac"].to_numpy(float),
                include_norm=include_norm,
                include_corr=include_corr,
            )

            row = {
                k: v
                for k, v in result.items()
                if k != "point_table"
            }
            row.update({
                "dataset": "pass2",
                "dataset_label": DATASET_LABELS["pass2"],
                "model": model_name,
                "scenario": scenario,
                "include_norm": bool(include_norm),
                "include_corr": bool(include_corr),
            })
            summary_rows.append(row)

            pt = result["point_table"].copy()
            if not pt.empty:
                # Reproduce the exact subset selected by the fitter.
                good = (
                    np.isfinite(d["xs"].to_numpy(float))
                    & np.isfinite(d[model_col].to_numpy(float))
                    & np.isfinite(d["point_unc_abs"].to_numpy(float))
                    & np.isfinite(d["norm_frac"].to_numpy(float))
                    & (d["xs"].to_numpy(float) > 0.0)
                    & (d[model_col].to_numpy(float) > 0.0)
                    & (d["point_unc_abs"].to_numpy(float) > 0.0)
                )
                if include_corr:
                    good &= np.isfinite(d["corr_scale_frac"].to_numpy(float))
                #endif
                dd = d.loc[good].reset_index(drop=True)

                for col in [
                    "point_id", "published_bin", "xB", "Q2", "t_abs",
                    "phi_deg", "ebeam", "xs", "stat_abs", "ptp_sys_abs",
                    "point_unc_abs", "norm_frac", "corr_scale_frac",
                ]:
                    if col in dd.columns:
                        pt[col] = dd[col].to_numpy()
                    #endif
                #endfor
                pt["model_name"] = model_name
                pt["scenario"] = scenario
                pt["beta_norm"] = float(result["beta_norm"])
                pt["beta_corr"] = float(result["beta_corr"])
                point_rows.append(pt)
            #endif
        #endfor
    #endfor

    return (
        pd.DataFrame(summary_rows),
        pd.concat(point_rows, ignore_index=True)
        if point_rows else pd.DataFrame(),
    )
#enddef


def fit_pass2_pair_publication_nuisances(
        pair: pd.DataFrame) -> Dict[str, float]:
    """
    Publication-style matched comparison of one external dataset A to Hayward B.

    Work in log-ratio space:

        z_i = ln(A_i / B_i).

    One relative overall-normalization nuisance eta_rel has prior width equal to
    the quadrature combination of the two experiments' quoted normalization
    uncertainties in log space.  One independent beta_corr nuisance multiplies
    Hayward's finalized bin-dependent correlated-scale response.

    The fitted residual is

        z_i + eta_rel - beta_corr*f_corr,i.

    This preserves the fact that only the *relative* experiment-wide
    normalization is identifiable from a two-dataset comparison, while still
    allowing the publication-level Hayward correlated-scale mode to act.
    """
    A = pair["xs_a_to_b_km15"].to_numpy(float)
    B = pair["xs_b"].to_numpy(float)
    SA = pair["point_unc_a_to_b_km15"].to_numpy(float)
    SB = pair["point_unc_b"].to_numpy(float)
    C = pair["corr_scale_frac_b"].to_numpy(float)

    good = (
        np.isfinite(A) & np.isfinite(B)
        & np.isfinite(SA) & np.isfinite(SB)
        & np.isfinite(C)
        & (A > 0.0) & (B > 0.0)
        & (SA > 0.0) & (SB > 0.0)
    )

    A, B, SA, SB, C = A[good], B[good], SA[good], SB[good], C[good]
    if len(A) == 0:
        return {}
    #endif

    z = np.log(A / B)
    s = np.sqrt((SA/A)**2 + (SB/B)**2)
    w = 1.0 / (s*s)

    norm_a = float(pair["norm_frac_a"].iloc[0])
    norm_b = float(pair["norm_frac_b"].iloc[0])
    sigma_eta = math.sqrt(
        math.log1p(max(norm_a, 0.0))**2
        + math.log1p(max(norm_b, 0.0))**2
    )

    # Design columns for [eta_relative, beta_corr].
    # residual = z + eta_relative - beta_corr*C
    X = np.column_stack([
        np.ones_like(z),
        -C,
    ])

    lhs = X.T @ (w[:, None] * X)
    rhs = -(X.T @ (w*z))

    # Gaussian priors: eta/sigma_eta and beta_corr/1.
    if sigma_eta > 0.0:
        lhs[0, 0] += 1.0 / (sigma_eta*sigma_eta)
    #endif
    lhs[1, 1] += 1.0

    pars = np.linalg.solve(lhs, rhs)
    eta_rel = float(pars[0])
    beta_corr = float(pars[1])

    residual = z + eta_rel - beta_corr*C
    pulls = residual / s

    prior = (
        (eta_rel/sigma_eta)**2 if sigma_eta > 0.0 else 0.0
    ) + beta_corr**2
    chi2_data = float(np.sum(pulls**2))
    chi2_total = chi2_data + float(prior)

    corr_shift_pct = 100.0 * (
        np.exp(-beta_corr*C) - 1.0
    )

    return {
        "N": int(len(z)),
        "relative_norm_beta": (
            eta_rel/sigma_eta if sigma_eta > 0.0 else 0.0
        ),
        "relative_scale_a_to_b": float(math.exp(eta_rel)),
        "beta_pass2_corr_scale": beta_corr,
        "combined_nuisance_excursion_sigma": float(
            math.sqrt(
                (
                    (eta_rel / sigma_eta)**2
                    if sigma_eta > 0.0 else 0.0
                )
                + beta_corr**2
            )
        ),
        "chi2_data": chi2_data,
        "chi2_prior": float(prior),
        "chi2_total": chi2_total,
        "chi2_per_match": float(chi2_total / len(z)),
        "fitted_pull_rms": float(np.sqrt(np.mean(pulls**2))),
        "corr_shift_median_pct": float(np.nanmedian(corr_shift_pct)),
        "corr_shift_min_pct": float(np.nanmin(corr_shift_pct)),
        "corr_shift_max_pct": float(np.nanmax(corr_shift_pct)),
    }
#enddef


def make_pass2_pair_publication_scores(
        matches: pd.DataFrame) -> pd.DataFrame:
    """
    Compare every external dataset directly to Hayward with the publication
    correlated-scale nuisance included.

    These results complement the generic pairwise table, whose chi2 includes
    relative overall normalization but NOT the Hayward correlated-scale mode.
    """
    rows = []
    for key in DATASET_ORDER:
        if key == "pass2":
            continue
        #endif

        pair = matches.loc[
            (matches["dataset_a"] == key)
            & (matches["dataset_b"] == "pass2")
        ].copy()

        if pair.empty:
            continue
        #endif

        result = fit_pass2_pair_publication_nuisances(pair)
        if not result:
            continue
        #endif

        result.update({
            "dataset_a": key,
            "dataset_a_label": DATASET_LABELS[key],
            "dataset_b": "pass2",
            "dataset_b_label": DATASET_LABELS["pass2"],
        })
        rows.append(result)
    #endfor

    return pd.DataFrame(rows)
#enddef


def plot_pass2_model_residual_diagnostics(
        point_table: pd.DataFrame,
        outdir: Path) -> None:
    """
    Plot raw and publication-nuisance-adjusted Hayward residuals versus the four
    primary kinematic coordinates for KM15 and BH.

    These are diagnostic figures; the phi-dependent cross-section canvases
    remain the primary presentation plots.
    """
    if point_table.empty:
        return
    #endif

    outdir.mkdir(parents=True, exist_ok=True)
    variables = [
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
        ("phi_deg", r"$\phi$ (deg)"),
    ]

    for model_name in ["KM15", "BH"]:
        d = point_table.loc[
            (point_table["model_name"] == model_name)
            & (point_table["scenario"] == "norm_plus_corr")
        ].copy()
        if d.empty:
            continue
        #endif

        for variable, xlabel in variables:
            fig, ax = plt.subplots(figsize=(8.6, 5.6))
            ax.scatter(
                d[variable],
                d["raw_pull"],
                s=12,
                alpha=0.45,
                label="raw: stat ⊕ point-to-point",
            )
            ax.scatter(
                d[variable],
                d["fitted_pull"],
                s=12,
                alpha=0.45,
                label="norm + correlated-scale nuisances",
            )
            ax.axhline(0.0, lw=1.0)
            ax.axhline(+1.0, lw=0.8, ls="--")
            ax.axhline(-1.0, lw=0.8, ls="--")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("pull")
            ax.set_ylim(-6.0, 6.0)
            ax.grid(alpha=0.18)
            ax.legend(frameon=False)
            ax.set_title(
                f"CLAS12 pass-2 Hayward vs {model_name}: "
                "publication-systematic residuals"
            )
            fig.tight_layout()
            fig.savefig(
                outdir / f"hayward_vs_{model_name.lower()}_pull_vs_{variable}.png",
                dpi=220,
            )
            plt.close(fig)
        #endfor
    #endfor
#enddef


def make_native_model_scores(world: pd.DataFrame, have_gk16: bool) -> pd.DataFrame:
    """
    Legacy/simple native-model diagnostics.

    IMPORTANT: for pass-2 this table uses only the overall-normalization
    nuisance.  The authoritative publication-style pass-2 model comparison is
    written separately by make_pass2_model_publication_scores().
    """
    models = [("km15", "km15_native"), ("bh", "bh_native")]
    if have_gk16:
        models.insert(1, ("gk16", "gk16_native"))
    #endif

    rows = []
    for key in DATASET_ORDER:
        d = world.loc[world["dataset"] == key]
        if d.empty:
            continue
        #endif
        norm = float(d["norm_frac"].iloc[0])
        for model_name, col in models:
            result = fit_one_normalization_nuisance(
                d["xs"].to_numpy(float),
                d[col].to_numpy(float),
                d["point_unc_abs"].to_numpy(float),
                norm,
            )
            result.update({
                "dataset": key,
                "dataset_label": DATASET_LABELS[key],
                "model": model_name,
                "normalization_fraction": norm,
                "systematic_treatment": (
                    "overall normalization only; NOT publication-complete"
                    if key == "pass2"
                    else "overall normalization"
                ),
            })
            rows.append(result)
        #endfor
    #endfor
    return pd.DataFrame(rows)
#enddef


# =============================================================================
# Pairwise kinematic matching and direct dataset consistency
# =============================================================================


def candidate_matches(a: pd.DataFrame, b: pd.DataFrame, cfg: MatchConfig) -> pd.DataFrame:
    """Construct all candidate A/B pairs inside explicit 4D matching windows."""
    rows = []

    # The datasets are only a few thousand points; grouping by coarse xB first
    # avoids a full O(N^2) Cartesian product while retaining transparent cuts.
    for ia, ra in a.iterrows():
        xbmask = np.abs(b["xB"].to_numpy(float) - float(ra["xB"])) <= cfg.dxb
        if not np.any(xbmask):
            continue
        #endif

        sub = b.loc[xbmask]
        dx = np.abs(sub["xB"].to_numpy(float) - float(ra["xB"]))
        dq = np.abs(sub["Q2"].to_numpy(float) - float(ra["Q2"]))
        dt = np.abs(sub["t_abs"].to_numpy(float) - float(ra["t_abs"]))
        dp = circular_phi_difference_deg(sub["phi_deg"].to_numpy(float), float(ra["phi_deg"]))

        keep = (dq <= cfg.dq2) & (dt <= cfg.dt) & (dp <= cfg.dphi)
        if not np.any(keep):
            continue
        #endif

        kept = sub.loc[keep]
        dx = dx[keep]
        dq = dq[keep]
        dt = dt[keep]
        dp = dp[keep]

        score = np.sqrt(
            (dx / cfg.dxb)**2
            + (dq / cfg.dq2)**2
            + (dt / cfg.dt)**2
            + (dp / cfg.dphi)**2
        )

        for j, (ib, rb) in enumerate(kept.iterrows()):
            rows.append({
                "ia": int(ia),
                "ib": int(ib),
                "point_id_a": str(ra["point_id"]),
                "point_id_b": str(rb["point_id"]),
                "dxB": float(dx[j]),
                "dQ2": float(dq[j]),
                "dt_abs": float(dt[j]),
                "dphi_deg": float(dp[j]),
                "match_score": float(score[j]),
            })
        #endfor
    #endfor

    return pd.DataFrame(rows)
#enddef


def greedy_one_to_one_matches(candidates: pd.DataFrame) -> pd.DataFrame:
    """Choose the closest nonconflicting candidate pairs deterministically."""
    if candidates.empty:
        return candidates.copy()
    #endif

    c = candidates.sort_values(
        ["match_score", "dxB", "dQ2", "dt_abs", "dphi_deg", "point_id_a", "point_id_b"]
    )
    used_a = set()
    used_b = set()
    selected = []
    for row in c.itertuples(index=False):
        if row.point_id_a in used_a or row.point_id_b in used_b:
            continue
        #endif
        used_a.add(row.point_id_a)
        used_b.add(row.point_id_b)
        selected.append(row._asdict())
    #endfor
    return pd.DataFrame(selected)
#enddef


def fit_pair_relative_normalization(
        a_value: np.ndarray,
        b_value: np.ndarray,
        sigma_a: np.ndarray,
        sigma_b: np.ndarray,
        transport_unc: np.ndarray,
        norm_a: float,
        norm_b: float) -> Dict[str, float]:
    """
    Fit the *relative* normalization of two matched datasets in log-ratio space.

    Only the relative normalization between two independent experiments is
    identifiable from an A-vs-B comparison.  Fitting independent absolute
    scales for A and B while keeping their point uncertainties fixed creates an
    unphysical common-mode direction in which both datasets can be scaled
    toward zero.  This routine removes that degeneracy.

    For each matched point, define

        z_i = ln(A_i / B_i)

    with approximate point uncertainty

        s_i^2 = (sigma_A/A)^2 + (sigma_B/B)^2
                + (sigma_transport/A)^2.

    A relative multiplicative shift exp(eta) is applied to A, so the residual
    becomes z_i + eta.  The independent quoted normalization uncertainties are
    combined in quadrature in log space and used as the prior width on eta.

    The returned scale_a_to_b = exp(eta) is therefore the factor multiplying
    dataset A (after local kinematic transport) to align it with dataset B.
    """
    A = np.asarray(a_value, dtype=float)
    B = np.asarray(b_value, dtype=float)
    SA = np.asarray(sigma_a, dtype=float)
    SB = np.asarray(sigma_b, dtype=float)
    ST = np.asarray(transport_unc, dtype=float)

    good = (
        np.isfinite(A) & np.isfinite(B) & np.isfinite(SA) & np.isfinite(SB)
        & np.isfinite(ST) & (A > 0.0) & (B > 0.0)
        & (SA >= 0.0) & (SB >= 0.0) & (ST >= 0.0)
    )
    A, B, SA, SB, ST = A[good], B[good], SA[good], SB[good], ST[good]
    if len(A) == 0:
        return {
            "N": 0,
            "eta": np.nan,
            "beta_relative": np.nan,
            "scale_a_to_b": np.nan,
            "relative_norm_frac": np.nan,
            "chi2": np.nan,
            "chi2_per_point": np.nan,
        }
    #endif

    z = np.log(A / B)
    sfrac = np.sqrt((SA / A)**2 + (SB / B)**2 + (ST / A)**2)
    finite = np.isfinite(z) & np.isfinite(sfrac) & (sfrac > 0.0)
    z, sfrac = z[finite], sfrac[finite]
    if len(z) == 0:
        return {
            "N": 0,
            "eta": np.nan,
            "beta_relative": np.nan,
            "scale_a_to_b": np.nan,
            "relative_norm_frac": np.nan,
            "chi2": np.nan,
            "chi2_per_point": np.nan,
        }
    #endif

    # Independent multiplicative normalization uncertainties combine in
    # quadrature.  log1p maps the quoted fractional scale width onto the same
    # additive variable eta used by ln(A/B).
    sigma_eta_a = math.log1p(max(float(norm_a), 0.0))
    sigma_eta_b = math.log1p(max(float(norm_b), 0.0))
    sigma_eta = math.sqrt(sigma_eta_a**2 + sigma_eta_b**2)

    w = 1.0 / (sfrac * sfrac)
    if sigma_eta > 0.0:
        eta = -float(np.sum(w * z)) / float(np.sum(w) + 1.0 / sigma_eta**2)
        penalty = (eta / sigma_eta)**2
        beta_relative = eta / sigma_eta
    else:
        eta = 0.0
        penalty = 0.0
        beta_relative = 0.0
    #endif

    pulls = (z + eta) / sfrac
    chi2 = float(np.sum(pulls**2) + penalty)
    return {
        "N": int(len(z)),
        "eta": float(eta),
        "beta_relative": float(beta_relative),
        "scale_a_to_b": float(math.exp(eta)),
        "relative_norm_frac": float(math.sqrt(norm_a**2 + norm_b**2)),
        "chi2": chi2,
        "chi2_per_point": float(chi2 / len(z)),
    }
#enddef

def build_pairwise_comparisons(
        world: pd.DataFrame,
        cfg: MatchConfig,
        have_gk16: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Match every dataset pair and transport A to B using local model ratios.

    KM15 is the sole transport prescription in this stage.  The
    transport-model uncertainty column is therefore zero by construction and
    is retained only for backward-compatible table structure.
    """
    all_matches = []
    summaries = []

    for i, key_a in enumerate(DATASET_ORDER):
        for key_b in DATASET_ORDER[i + 1:]:
            a = world.loc[world["dataset"] == key_a].copy().reset_index(drop=True)
            b = world.loc[world["dataset"] == key_b].copy().reset_index(drop=True)
            if a.empty or b.empty:
                continue
            #endif

            candidates = candidate_matches(a, b, cfg)
            selected = greedy_one_to_one_matches(candidates)
            if selected.empty:
                summaries.append({
                    "dataset_a": key_a,
                    "dataset_b": key_b,
                    "N_candidates": len(candidates),
                    "N_matches": 0,
                })
                continue
            #endif

            rows = []
            for m in selected.itertuples(index=False):
                ra = a.iloc[int(m.ia)]
                rb = b.iloc[int(m.ib)]

                c_km15 = float(rb["km15_native"] / ra["km15_native"])
                a_to_b_km15 = float(ra["xs"] * c_km15)
                a_stat = float(ra["stat_abs"] * abs(c_km15))
                a_sys = float(ra["ptp_sys_abs"] * abs(c_km15))
                a_point = float(ra["point_unc_abs"] * abs(c_km15))

                c_gk16 = np.nan
                a_to_b_gk16 = np.nan
                # KM15 is the nominal and sole transport prescription in the
                # current world-data study, so no additional transport-model
                # systematic is assigned here.
                transport_unc = 0.0
                if have_gk16 and np.isfinite(ra["gk16_native"]) and np.isfinite(rb["gk16_native"]):
                    c_gk16 = float(rb["gk16_native"] / ra["gk16_native"])
                    a_to_b_gk16 = float(ra["xs"] * c_gk16)
                    transport_unc = abs(a_to_b_gk16 - a_to_b_km15)
                #endif

                sigma_compare = math.sqrt(
                    a_point**2 + float(rb["point_unc_abs"])**2 + transport_unc**2
                )
                raw_pull = (
                    (a_to_b_km15 - float(rb["xs"])) / sigma_compare
                    if sigma_compare > 0.0 else np.nan
                )

                rows.append({
                    "dataset_a": key_a,
                    "dataset_b": key_b,
                    "point_id_a": str(ra["point_id"]),
                    "point_id_b": str(rb["point_id"]),
                    "match_score": float(m.match_score),
                    "dxB": float(m.dxB),
                    "dQ2": float(m.dQ2),
                    "dt_abs": float(m.dt_abs),
                    "dphi_deg": float(m.dphi_deg),
                    "xB_a": float(ra["xB"]), "xB_b": float(rb["xB"]),
                    "Q2_a": float(ra["Q2"]), "Q2_b": float(rb["Q2"]),
                    "t_abs_a": float(ra["t_abs"]), "t_abs_b": float(rb["t_abs"]),
                    "phi_a": float(ra["phi_deg"]), "phi_b": float(rb["phi_deg"]),
                    "ebeam_a": float(ra["ebeam"]), "ebeam_b": float(rb["ebeam"]),
                    "xs_a": float(ra["xs"]), "xs_b": float(rb["xs"]),
                    "stat_b": float(rb["stat_abs"]),
                    "ptp_sys_b": float(rb["ptp_sys_abs"]),
                    "point_unc_a": float(ra["point_unc_abs"]),
                    "point_unc_b": float(rb["point_unc_abs"]),
                    "km15_a": float(ra["km15_native"]),
                    "km15_b": float(rb["km15_native"]),
                    "bh_b": float(rb["bh_native"]),
                    "published_bin_b": (
                        float(rb["published_bin"])
                        if "published_bin" in rb.index and np.isfinite(rb["published_bin"])
                        else np.nan
                    ),
                    "xBmin_b": float(rb.get("xBmin", np.nan)),
                    "xBmax_b": float(rb.get("xBmax", np.nan)),
                    "Q2min_b": float(rb.get("Q2min", np.nan)),
                    "Q2max_b": float(rb.get("Q2max", np.nan)),
                    "t_abs_min_b": float(rb.get("t_abs_min", np.nan)),
                    "t_abs_max_b": float(rb.get("t_abs_max", np.nan)),
                    "km15_local_transport_factor": c_km15,
                    "xs_a_to_b_km15": a_to_b_km15,
                    "stat_a_to_b_km15": a_stat,
                    "ptp_sys_a_to_b_km15": a_sys,
                    "point_unc_a_to_b_km15": a_point,
                    "gk16_local_transport_factor": c_gk16,
                    "xs_a_to_b_gk16": a_to_b_gk16,
                    "transport_model_unc_abs": transport_unc,
                    "comparison_unc_abs": sigma_compare,
                    "raw_pull": raw_pull,
                    "norm_frac_a": float(ra["norm_frac"]),
                    "norm_frac_b": float(rb["norm_frac"]),
                    "corr_scale_frac_a": float(ra.get("corr_scale_frac", 0.0)),
                    "corr_scale_frac_b": float(rb.get("corr_scale_frac", 0.0)),
                })
            #endfor

            pair = pd.DataFrame(rows)
            relfit = fit_pair_relative_normalization(
                pair["xs_a_to_b_km15"].to_numpy(float),
                pair["xs_b"].to_numpy(float),
                pair["point_unc_a_to_b_km15"].to_numpy(float),
                pair["point_unc_b"].to_numpy(float),
                pair["transport_model_unc_abs"].to_numpy(float),
                float(pair["norm_frac_a"].iloc[0]),
                float(pair["norm_frac_b"].iloc[0]),
            )
            relative_scale = float(relfit["scale_a_to_b"])
            pair["relative_norm_beta"] = float(relfit["beta_relative"])
            pair["relative_scale_a_to_b"] = relative_scale
            pair["relative_norm_frac"] = float(relfit["relative_norm_frac"])

            pair["raw_log_pull"] = np.log(pair["xs_a_to_b_km15"] / pair["xs_b"]) / np.sqrt(
                (pair["point_unc_a_to_b_km15"] / pair["xs_a_to_b_km15"])**2
                + (pair["point_unc_b"] / pair["xs_b"])**2
                + (pair["transport_model_unc_abs"] / pair["xs_a_to_b_km15"])**2
            )
            pair["profiled_pull"] = (
                np.log(pair["xs_a_to_b_km15"] / pair["xs_b"]) + math.log(relative_scale)
            ) / np.sqrt(
                (pair["point_unc_a_to_b_km15"] / pair["xs_a_to_b_km15"])**2
                + (pair["point_unc_b"] / pair["xs_b"])**2
                + (pair["transport_model_unc_abs"] / pair["xs_a_to_b_km15"])**2
            )

            summaries.append({
                "dataset_a": key_a,
                "dataset_b": key_b,
                "N_candidates": len(candidates),
                "N_matches": len(pair),
                "relative_norm_beta": float(relfit["beta_relative"]),
                "relative_scale_a_to_b": relative_scale,
                "relative_norm_pct": 100.0 * float(relfit["relative_norm_frac"]),
                "chi2_with_norm_penalty": float(relfit["chi2"]),
                "chi2_per_match": float(relfit["chi2_per_point"]),
                "raw_pull_rms": float(np.sqrt(np.nanmean(pair["raw_pull"]**2))),
                "raw_log_pull_rms": float(np.sqrt(np.nanmean(pair["raw_log_pull"]**2))),
                "profiled_pull_rms": float(np.sqrt(np.nanmean(pair["profiled_pull"]**2))),
                "median_abs_profiled_pull": float(np.nanmedian(np.abs(pair["profiled_pull"]))),
                "median_transport_model_unc_pct": (
                    100.0 * float(np.nanmedian(pair["transport_model_unc_abs"] / np.abs(pair["xs_a_to_b_km15"])))
                ),
            })
            all_matches.append(pair)
        #endfor
    #endfor

    match_table = pd.concat(all_matches, ignore_index=True) if all_matches else pd.DataFrame()
    return match_table, pd.DataFrame(summaries)
#enddef


# =============================================================================
# Plotting
# =============================================================================


# Fixed display windows keep a handful of pathological/very-low-weight points
# from determining the visual scale of summary figures.  The underlying points
# remain in all tables and quantitative calculations; matplotlib simply clips
# markers/error bars outside these display ranges.
NATIVE_RATIO_YLIMS = {
    "KM15": (0.0, 3.0),
    "GK16": (0.0, 3.0),
    "BH": (0.0, 10.0),
}
PAIRWISE_PULL_YLIM = (-5.5, 5.5)

# Cross-section presentation panels.  The 3x4 layout matches the standard
# analysis-note figure style used elsewhere in this analysis.
PANEL_NROWS = 3
PANEL_NCOLS = 4
PANEL_PER_PAGE = PANEL_NROWS * PANEL_NCOLS
PANEL_MIN_MATCHES = 4
PANEL_MAX_PAGES_PER_PAIR = 2
PANEL_MAX_WORLD_ANCHOR_PAGES = 0  # 0 = plot all qualifying Lee-anchor bins

# Nearby points in the reference (B) dataset are grouped into one phi panel
# when their hadronic kinematics are within these tighter presentation-scale
# windows.  These are deliberately much tighter than the cross-experiment
# matching windows and do not affect the quantitative matching itself.
PANEL_CELL_DXB = 0.012
PANEL_CELL_DQ2 = 0.20   # GeV^2
PANEL_CELL_DT = 0.040   # GeV^2

# A handful of published points have enormous quoted uncertainties.  They are
# retained in all quantitative calculations.  For presentation only, points
# with point_unc/xs above this threshold are drawn as open markers without an
# error bar so one pathological uncertainty cannot cover an entire panel.
PLOT_MAX_REL_POINT_UNC = 1.0


def _annotate_clipped_y(ax, values: pd.Series, ylow: float, yhigh: float) -> None:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(float)
    finite = np.isfinite(arr)
    nlow = int(np.sum(finite & (arr < ylow)))
    nhigh = int(np.sum(finite & (arr > yhigh)))
    nclip = nlow + nhigh
    if nclip <= 0:
        return
    #endif
    ax.text(
        0.012, 0.018,
        f"{nclip} central point(s) outside display range; all retained in analysis",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, alpha=0.75,
    )
#enddef


def plot_native_model_ratios(world: pd.DataFrame, outdir: Path, have_gk16: bool) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    variables = [
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
        ("phi_deg", r"$\phi$ (deg)"),
    ]

    model_specs = [("data_over_km15", "KM15"), ("data_over_bh", "BH")]
    if have_gk16:
        model_specs.insert(1, ("data_over_gk16", "GK16"))
    #endif

    for ratio_col, model_label in model_specs:
        for variable, xlabel in variables:
            fig, ax = plt.subplots(figsize=(8.4, 5.6))
            for key in DATASET_ORDER:
                d = world.loc[world["dataset"] == key]
                if d.empty:
                    continue
                #endif
                relerr = d["point_unc_abs"] / d[
                    "km15_native" if model_label == "KM15" else (
                        "gk16_native" if model_label == "GK16" else "bh_native"
                    )
                ]
                # Keep pathological published uncertainties in the analysis
                # but suppress their gigantic bars in this summary plot.
                ratio_values = d[ratio_col].to_numpy(float)
                relerr_values = np.asarray(relerr, dtype=float)
                rel_to_data = d["point_unc_abs"].to_numpy(float) / d["xs"].to_numpy(float)
                normal = np.isfinite(rel_to_data) & (rel_to_data <= PLOT_MAX_REL_POINT_UNC)
                extreme = ~normal

                style = DATASET_STYLES[key]
                if np.any(normal):
                    ax.errorbar(
                        d.loc[normal, variable], ratio_values[normal],
                        yerr=relerr_values[normal],
                        fmt=style["marker"], ms=2.8, lw=0.7, capsize=0,
                        color=style["color"],
                        markeredgecolor=style["color"],
                        alpha=0.60, label=DATASET_LABELS[key],
                    )
                elif np.any(extreme):
                    ax.plot(
                        d.loc[extreme, variable], ratio_values[extreme],
                        linestyle="none", marker=style["marker"], ms=2.8,
                        color=style["color"],
                        markeredgecolor=style["color"],
                        markerfacecolor="none", alpha=0.60,
                        label=DATASET_LABELS[key],
                    )
                #endif

                if np.any(extreme) and np.any(normal):
                    ax.plot(
                        d.loc[extreme, variable], ratio_values[extreme],
                        linestyle="none", marker=style["marker"], ms=3.0,
                        color=style["color"],
                        markeredgecolor=style["color"],
                        markerfacecolor="none", alpha=0.60,
                        label="_nolegend_",
                    )
                #endif
            #endfor
            ax.axhline(1.0, lw=1.1, linestyle="--")
            ylow, yhigh = NATIVE_RATIO_YLIMS[model_label]
            ax.set_ylim(ylow, yhigh)
            _annotate_clipped_y(ax, world[ratio_col], ylow, yhigh)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(rf"$\sigma_{{\rm data}}/\sigma_{{\rm {model_label}}}$")
            ax.set_title(f"Published world data / {model_label} at native kinematics")
            ax.grid(alpha=0.20)
            ax.legend(fontsize=8, ncol=2, frameon=False)
            fig.tight_layout()
            tag = variable.replace("_abs", "")
            fig.savefig(outdir / f"native_data_over_{model_label.lower()}_vs_{tag}.png", dpi=220)
            plt.close(fig)
        #endfor
    #endfor
#enddef


def plot_transport_uncertainty(world: pd.DataFrame, outdir: Path, have_gk16: bool) -> None:
    if not have_gk16:
        return
    #endif

    outdir.mkdir(parents=True, exist_ok=True)
    for variable, xlabel in [
        ("ebeam", r"Native $E_{\rm beam}$ (GeV)"),
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
    ]:
        fig, ax = plt.subplots(figsize=(8.4, 5.6))
        for key in DATASET_ORDER:
            d = world.loc[world["dataset"] == key]
            ax.scatter(
                d[variable], 100.0 * d["transport_model_unc_frac"],
                s=14, alpha=0.55, label=DATASET_LABELS[key],
            )
        #endfor
        ax.set_xlabel(xlabel)
        ax.set_ylabel("KM15–GK16 transport excursion (%)")
        ax.set_title(rf"Model dependence of transport to $E_{{\rm beam}}={TARGET_EBEAM_GEV:.1f}$ GeV")
        ax.grid(alpha=0.20)
        ax.legend(fontsize=8, ncol=2, frameon=False)
        fig.tight_layout()
        tag = variable.replace("_abs", "")
        fig.savefig(outdir / f"transport_model_uncertainty_vs_{tag}.png", dpi=220)
        plt.close(fig)
    #endfor
#enddef


def plot_pairwise_matrix(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    #endif
    outdir.mkdir(parents=True, exist_ok=True)

    n = len(DATASET_ORDER)
    matrix = np.full((n, n), np.nan)
    counts = np.zeros((n, n), dtype=int)
    index = {k: i for i, k in enumerate(DATASET_ORDER)}

    for row in summary.itertuples(index=False):
        i = index[row.dataset_a]
        j = index[row.dataset_b]
        if int(row.N_matches) > 0 and np.isfinite(getattr(row, "chi2_per_match", np.nan)):
            matrix[i, j] = float(row.chi2_per_match)
            matrix[j, i] = float(row.chi2_per_match)
            counts[i, j] = counts[j, i] = int(row.N_matches)
        #endif
    #endfor

    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    im = ax.imshow(matrix, origin="upper", aspect="equal")
    labels = [DATASET_LABELS[k] for k in DATASET_ORDER]
    ax.set_xticks(np.arange(n), labels=labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(n), labels=labels)
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.text(j, i, "—", ha="center", va="center")
            elif np.isfinite(matrix[i, j]):
                ax.text(j, i, f"{matrix[i,j]:.2f}\nN={counts[i,j]}", ha="center", va="center", fontsize=8)
            #endif
        #endfor
    #endfor
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\chi^2/N_{\rm match}$ after normalization shifts")
    ax.set_title("Pairwise world-data compatibility at matched kinematics")
    fig.tight_layout()
    fig.savefig(outdir / "pairwise_compatibility_matrix.png", dpi=220)
    plt.close(fig)
#enddef


def plot_pairwise_pulls(matches: pd.DataFrame, outdir: Path) -> None:
    if matches.empty:
        return
    #endif
    outdir.mkdir(parents=True, exist_ok=True)

    for variable, xlabel, col in [
        ("t", r"$|t|_B$ (GeV$^2$)", "t_abs_b"),
        ("xB", r"$x_{B,B}$", "xB_b"),
        ("Q2", r"$Q_B^2$ (GeV$^2$)", "Q2_b"),
        ("phi", r"$\phi_B$ (deg)", "phi_b"),
    ]:
        fig, ax = plt.subplots(figsize=(8.4, 5.6))
        for (ka, kb), d in matches.groupby(["dataset_a", "dataset_b"], sort=False):
            style = DATASET_STYLES.get(str(ka), {"color": None, "marker": "o"})
            ax.scatter(
                d[col], d["profiled_pull"],
                s=15, alpha=0.50,
                color=style["color"],
                marker=style["marker"],
                label=f"{DATASET_LABELS[ka]} → {DATASET_LABELS[kb]}",
            )
        #endfor
        ax.axhline(0.0, lw=1.0)
        ax.axhline(+1.0, lw=0.8, linestyle="--")
        ax.axhline(-1.0, lw=0.8, linestyle="--")
        ax.set_ylim(*PAIRWISE_PULL_YLIM)
        _annotate_clipped_y(
            ax, matches["profiled_pull"], PAIRWISE_PULL_YLIM[0], PAIRWISE_PULL_YLIM[1]
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Matched-data residual / pointwise uncertainty")
        ax.set_title("Pairwise matched world-data residuals")
        ax.grid(alpha=0.20)
        ax.legend(fontsize=6.5, ncol=2, frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / f"pairwise_profiled_pulls_vs_{variable}.png", dpi=220)
        plt.close(fig)
    #endfor
#enddef



def _draw_measurement_series(
        ax,
        x: np.ndarray,
        y: np.ndarray,
        yerr: np.ndarray,
        *,
        label: str,
        dataset_key: Optional[str] = None,
        marker: Optional[str] = None,
        color: Optional[str] = None,
        xoffset: float = 0.0,
        markersize: float = 4.2,
        alpha: float = 0.90):
    """
    Draw one measured cross-section series without allowing pathological
    published uncertainties to dominate the panel visually.

    A dataset key supplies the universal color/marker identity used throughout
    the script.  Explicit marker/color arguments remain available for special
    cases but normally should not be needed.

    Points with point_unc/xs <= PLOT_MAX_REL_POINT_UNC are drawn with their full
    pointwise (stat ⊕ point-to-point systematic) error bar.  Larger-uncertainty
    points remain visible as open markers but their error bar is suppressed for
    presentation only.  No point is removed from any fit or CSV.
    """
    if dataset_key is not None:
        style = DATASET_STYLES.get(str(dataset_key), {})
        if marker is None:
            marker = style.get("marker", "o")
        #endif
        if color is None:
            color = style.get("color", None)
        #endif
    #endif

    if marker is None:
        marker = "o"
    #endif

    x = np.asarray(x, dtype=float) + float(xoffset)
    y = np.asarray(y, dtype=float)
    yerr = np.asarray(yerr, dtype=float)

    good = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (y > 0.0) & (yerr >= 0.0)
    if not np.any(good):
        return None
    #endif

    x = x[good]
    y = y[good]
    yerr = yerr[good]
    rel = yerr / y
    normal = np.isfinite(rel) & (rel <= PLOT_MAX_REL_POINT_UNC)
    extreme = ~normal

    handle = None
    if np.any(normal):
        handle = ax.errorbar(
            x[normal], y[normal], yerr=yerr[normal],
            fmt=marker, ms=markersize, lw=0.85, capsize=2.0,
            alpha=alpha, label=label,
            color=color,
            markeredgecolor=color,
        )
    #endif

    if np.any(extreme):
        open_handle = ax.plot(
            x[extreme], y[extreme],
            linestyle="none", marker=marker, ms=markersize + 0.4,
            markerfacecolor="none", markeredgewidth=1.0,
            markeredgecolor=color,
            color=color,
            alpha=alpha,
            label=(label if handle is None else "_nolegend_"),
        )[0]
        if handle is None:
            handle = open_handle
        #endif
    #endif

    return handle
#enddef

def _robust_positive_log_limits(values: Sequence[np.ndarray]) -> Tuple[float, float]:
    """
    Choose stable log-scale limits from central values/model curves only.

    Error bars intentionally do not enter this calculation.  This prevents a
    single gigantic published uncertainty from setting the y-axis range.
    """
    arrays = []
    for value in values:
        arr = np.asarray(value, dtype=float).ravel()
        arr = arr[np.isfinite(arr) & (arr > 0.0)]
        if arr.size:
            arrays.append(arr)
        #endif
    #endfor

    if not arrays:
        return 1.0e-6, 1.0
    #endif

    v = np.concatenate(arrays)
    if len(v) >= 10:
        vlo = float(np.nanpercentile(v, 2.0))
        vhi = float(np.nanpercentile(v, 98.0))
        # Do not crop ordinary central values unless a truly wild outlier is
        # present; percentile limits are expanded generously below.
        positive_min = float(np.nanmin(v))
        positive_max = float(np.nanmax(v))
        if positive_min > 0.25 * vlo:
            vlo = positive_min
        #endif
        if positive_max < 4.0 * vhi:
            vhi = positive_max
        #endif
    else:
        vlo = float(np.nanmin(v))
        vhi = float(np.nanmax(v))
    #endif

    vlo = max(vlo * 0.45, 1.0e-12)
    vhi = max(vhi * 2.2, vlo * 10.0)
    return vlo, vhi
#enddef



def _synchronize_canvas_y_limits(
        axes,
        mode: str = "row") -> None:
    """
    Synchronize log-scale y limits within a multipanel canvas.

    mode = "panel":
        leave every subplot at its independently determined robust range.

    mode = "row":
        each horizontal row shares one common y range.  This is the default
        because adjacent panels remain directly comparable without allowing
        one extreme BH/KM15 endpoint to compress an entire 3x4 page.

    mode = "page":
        all active subplots on the page share one y range.

    Every panel has already computed a robust range that ignores pathological
    uncertainty bars.  This function combines those already-robust ranges; it
    does not inspect the raw error bars again.
    """
    mode = str(mode).strip().lower()
    if mode not in {"panel", "row", "page"}:
        raise ValueError(
            f"Unknown y-scale synchronization mode '{mode}'. "
            "Choose panel, row, or page."
        )
    #endif

    if mode == "panel":
        return
    #endif

    arr = np.asarray(axes, dtype=object)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    #endif

    def apply_group(group_axes) -> None:
        active = [
            ax for ax in group_axes
            if ax.get_visible() and ax.has_data()
        ]
        if not active:
            return
        #endif

        lows = []
        highs = []
        for ax in active:
            lo, hi = ax.get_ylim()
            if (
                np.isfinite(lo)
                and np.isfinite(hi)
                and lo > 0.0
                and hi > lo
            ):
                lows.append(float(lo))
                highs.append(float(hi))
            #endif
        #endfor

        if not lows or not highs:
            return
        #endif

        common_lo = min(lows)
        common_hi = max(highs)

        for ax in active:
            ax.set_ylim(common_lo, common_hi)
        #endfor
    #enddef

    if mode == "page":
        apply_group(list(arr.ravel()))
        return
    #endif

    # Row-wise mode.
    for irow in range(arr.shape[0]):
        apply_group(list(arr[irow, :]))
    #endfor
#enddef



def _cluster_reference_kinematic_cells(pair: pd.DataFrame) -> pd.DataFrame:
    """
    Assign matched rows to reference-dataset (B) hadronic-kinematic cells.

    If B is Lee 2026 and an authoritative published bin index is available, that
    index is used directly.  Otherwise rows are greedily clustered in
    (xB,Q2,|t|) using presentation-scale windows that are much tighter than the
    cross-experiment matching cuts.
    """
    if pair.empty:
        out = pair.copy()
        out["panel_group"] = pd.Series(dtype=str)
        return out
    #endif

    out = pair.copy().reset_index(drop=True)

    if "published_bin_b" in out.columns:
        finite_bin = np.isfinite(pd.to_numeric(out["published_bin_b"], errors="coerce"))
        if int(finite_bin.sum()) == len(out):
            out["panel_group"] = [
                f"bin_{int(round(v))}"
                for v in out["published_bin_b"].to_numpy(float)
            ]
            return out
        #endif
    #endif

    centers: List[Dict[str, object]] = []
    group_ids: List[int] = []

    order = np.lexsort((
        out["t_abs_b"].to_numpy(float),
        out["Q2_b"].to_numpy(float),
        out["xB_b"].to_numpy(float),
    ))

    assigned = np.full(len(out), -1, dtype=int)

    for idx in order:
        xb = float(out.at[idx, "xB_b"])
        q2 = float(out.at[idx, "Q2_b"])
        tt = float(out.at[idx, "t_abs_b"])

        best_group = None
        best_score = np.inf
        for ig, center in enumerate(centers):
            dx = abs(xb - float(center["xB"]))
            dq = abs(q2 - float(center["Q2"]))
            dt = abs(tt - float(center["t_abs"]))
            if dx <= PANEL_CELL_DXB and dq <= PANEL_CELL_DQ2 and dt <= PANEL_CELL_DT:
                score = math.sqrt(
                    (dx / PANEL_CELL_DXB)**2
                    + (dq / PANEL_CELL_DQ2)**2
                    + (dt / PANEL_CELL_DT)**2
                )
                if score < best_score:
                    best_group = ig
                    best_score = score
                #endif
            #endif
        #endfor

        if best_group is None:
            centers.append({
                "xB": xb,
                "Q2": q2,
                "t_abs": tt,
                "members": [int(idx)],
            })
            assigned[idx] = len(centers) - 1
        else:
            members = list(centers[best_group]["members"])
            members.append(int(idx))
            centers[best_group]["members"] = members
            centers[best_group]["xB"] = float(np.median(out.loc[members, "xB_b"]))
            centers[best_group]["Q2"] = float(np.median(out.loc[members, "Q2_b"]))
            centers[best_group]["t_abs"] = float(np.median(out.loc[members, "t_abs_b"]))
            assigned[idx] = best_group
        #endif
    #endfor

    out["panel_group"] = [f"cell_{g:04d}" for g in assigned]
    return out
#enddef


def make_pairwise_panel_summary(matches: pd.DataFrame) -> pd.DataFrame:
    """Summarize the candidate phi-dependent presentation panels."""
    if matches.empty:
        return pd.DataFrame()
    #endif

    rows = []
    for (ka, kb), pair0 in matches.groupby(["dataset_a", "dataset_b"], sort=False):
        pair = _cluster_reference_kinematic_cells(pair0)
        for group, d in pair.groupby("panel_group", sort=False):
            if len(d) < PANEL_MIN_MATCHES:
                continue
            #endif

            phi = np.sort(np.mod(d["phi_b"].to_numpy(float), 360.0))
            if len(phi) >= 2:
                gaps = np.diff(np.r_[phi, phi[0] + 360.0])
                phi_coverage = 360.0 - float(np.max(gaps))
            else:
                phi_coverage = 0.0
            #endif

            rows.append({
                "dataset_a": ka,
                "dataset_b": kb,
                "panel_group": group,
                "N_matches": int(len(d)),
                "xB_ref": float(np.median(d["xB_b"])),
                "Q2_ref": float(np.median(d["Q2_b"])),
                "t_abs_ref": float(np.median(d["t_abs_b"])),
                "ebeam_ref": float(np.median(d["ebeam_b"])),
                "phi_coverage_deg": phi_coverage,
                "raw_pull_rms": float(np.sqrt(np.nanmean(d["raw_pull"]**2))),
                "median_match_score": float(np.nanmedian(d["match_score"])),
                "panel_rank_score": float(len(d) + phi_coverage / 360.0),
            })
        #endfor
    #endfor

    return pd.DataFrame(rows)
#enddef


def _plot_one_pair_cross_section_panel(
        ax,
        d: pd.DataFrame,
        dataset_a: str,
        dataset_b: str,
        emff,
        model_curve_cache: Dict[Tuple, pd.DataFrame],
        *,
        show_legend_labels: bool = True) -> None:
    """Draw one phi-dependent A->B matched cross-section panel."""
    d = d.sort_values("phi_b").copy()

    phi = d["phi_b"].to_numpy(float)
    offset = 1.4

    norm_a = 100.0 * float(d["norm_frac_a"].iloc[0])
    norm_b = 100.0 * float(d["norm_frac_b"].iloc[0])

    label_a = (
        f"{DATASET_LABELS[dataset_a]} → reference ({norm_a:.1f}% norm)"
        if show_legend_labels else "_nolegend_"
    )
    label_b = (
        f"{DATASET_LABELS[dataset_b]} ({norm_b:.1f}% norm)"
        if show_legend_labels else "_nolegend_"
    )

    _draw_measurement_series(
        ax,
        phi,
        d["xs_a_to_b_km15"].to_numpy(float),
        d["point_unc_a_to_b_km15"].to_numpy(float),
        label=label_a,
        dataset_key=dataset_a,
        xoffset=-offset,
        markersize=4.2,
    )
    _draw_measurement_series(
        ax,
        phi,
        d["xs_b"].to_numpy(float),
        d["point_unc_b"].to_numpy(float),
        label=label_b,
        dataset_key=dataset_b,
        xoffset=+offset,
        markersize=4.0,
    )

    # IMPORTANT: evaluate the presentation curves at one fixed reference
    # kinematic point.  The previous implementation connected model values
    # evaluated at each individual matched data point, whose xB/Q2/|t| means
    # can differ slightly across phi; that can visibly distort the curve.
    xb = float(np.median(d["xB_b"]))
    q2 = float(np.median(d["Q2_b"]))
    tt = float(np.median(d["t_abs_b"]))
    ebeam = float(np.median(d["ebeam_b"]))
    model = get_dense_model_curve(
        emff,
        model_curve_cache,
        dataset_key=dataset_b,
        ebeam=ebeam,
        xB=xb,
        Q2=q2,
        t_abs=tt,
    )

    bh_style = MODEL_STYLES["bh"]
    km_style = MODEL_STYLES["km15"]
    ax.plot(
        model["phi_deg"], model["bh"],
        color=bh_style["color"],
        linestyle=bh_style["linestyle"],
        lw=bh_style["linewidth"],
        label=("BH" if show_legend_labels else "_nolegend_"),
    )
    ax.plot(
        model["phi_deg"], model["km15"],
        color=km_style["color"],
        linestyle=km_style["linestyle"],
        lw=km_style["linewidth"],
        label=("KM15" if show_legend_labels else "_nolegend_"),
    )

    ylo, yhi = _robust_positive_log_limits([
        d["xs_a_to_b_km15"].to_numpy(float),
        d["xs_b"].to_numpy(float),
        model["km15"].to_numpy(float),
        model["bh"].to_numpy(float),
    ])
    ax.set_yscale("log")
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.grid(alpha=0.18)

    ax.set_title(
        rf"$x_B={xb:.3f}$, $Q^2={q2:.2f}$, $|t|={tt:.3f}$",
        fontsize=9.0,
    )
#enddef


def _model_curve_cache_key(
        dataset_key: str,
        ebeam: float,
        xB: float,
        Q2: float,
        t_abs: float) -> Tuple:
    """Return the same rounded in-memory key used by get_dense_model_curve()."""
    return (
        str(dataset_key),
        round(float(ebeam), 6),
        round(float(xB), 6),
        round(float(Q2), 6),
        round(float(t_abs), 6),
        round(float(MODEL_CURVE_PHI_STEP_DEG), 6),
    )
#enddef


def collect_presentation_model_curve_specs(
        matches: pd.DataFrame,
        pair_panel_summary: pd.DataFrame,
        lee_anchor_summary: pd.DataFrame,
        max_pages_per_pair: int = PANEL_MAX_PAGES_PER_PAIR,
        max_lee_pages: int = PANEL_MAX_WORLD_ANCHOR_PAGES) -> List[Dict[str, float]]:
    """
    Collect and de-duplicate every fixed-kinematics BH/KM15 phi curve that will
    be needed by the pairwise and Lee-anchor presentation figures.

    Doing this before drawing any figure lets the terminal report a genuine
    overall percent-complete value for model prediction generation.
    """
    specs_by_key: Dict[Tuple, Dict[str, float]] = {}

    # Pairwise panels.
    if not matches.empty and not pair_panel_summary.empty:
        for (ka, kb), summary_pair in pair_panel_summary.groupby(
                ["dataset_a", "dataset_b"], sort=False):
            pair = matches.loc[
                (matches["dataset_a"] == ka)
                & (matches["dataset_b"] == kb)
            ].copy()
            pair = _cluster_reference_kinematic_cells(pair)

            ranked = summary_pair.sort_values(
                ["panel_rank_score", "N_matches", "phi_coverage_deg"],
                ascending=False,
            ).reset_index(drop=True)

            max_panels = min(
                len(ranked),
                int(max_pages_per_pair) * PANEL_PER_PAGE,
            )
            ranked = ranked.iloc[:max_panels].copy()

            for r in ranked.itertuples(index=False):
                group = str(r.panel_group)
                d = pair.loc[pair["panel_group"] == group]
                if d.empty:
                    continue
                #endif

                spec = {
                    "dataset_key": str(kb),
                    "ebeam": float(np.median(d["ebeam_b"])),
                    "xB": float(np.median(d["xB_b"])),
                    "Q2": float(np.median(d["Q2_b"])),
                    "t_abs": float(np.median(d["t_abs_b"])),
                    "source": f"pair:{ka}->{kb}:{group}",
                }
                key = _model_curve_cache_key(
                    spec["dataset_key"],
                    spec["ebeam"],
                    spec["xB"],
                    spec["Q2"],
                    spec["t_abs"],
                )
                specs_by_key.setdefault(key, spec)
            #endfor
        #endfor
    #endif

    # Lee-anchor panels.  The fitted-normalization variants reuse exactly the
    # same model curves, so each Lee bin is needed only once.
    if not matches.empty and not lee_anchor_summary.empty:
        lee = matches.loc[matches["dataset_b"] == "lee2026"].copy()
        if (
            not lee.empty
            and "published_bin_b" in lee.columns
            and np.isfinite(lee["published_bin_b"]).any()
        ):
            lee["anchor_group"] = [
                (
                    f"bin_{int(round(v))}"
                    if np.isfinite(v) else f"point_{pid}"
                )
                for v, pid in zip(
                    lee["published_bin_b"].to_numpy(float),
                    lee["point_id_b"].astype(str),
                )
            ]

            selected = lee_anchor_summary.copy()
            selected = selected.sort_values(
                ["published_bin", "xB_ref", "Q2_ref", "t_abs_ref"],
                ascending=True,
            ).reset_index(drop=True)
            if int(max_lee_pages) > 0:
                selected = selected.head(
                    int(max_lee_pages) * PANEL_PER_PAGE
                ).copy()
            #endif

            for r in selected.itertuples(index=False):
                group = str(r.anchor_group)
                d = lee.loc[lee["anchor_group"] == group].copy()
                if d.empty:
                    continue
                #endif
                lee_points = (
                    d.sort_values("phi_b")
                    .drop_duplicates("point_id_b")
                )
                spec = {
                    "dataset_key": "lee2026",
                    "ebeam": float(np.median(lee_points["ebeam_b"])),
                    "xB": float(np.median(lee_points["xB_b"])),
                    "Q2": float(np.median(lee_points["Q2_b"])),
                    "t_abs": float(np.median(lee_points["t_abs_b"])),
                    "source": f"lee:{group}",
                }
                key = _model_curve_cache_key(
                    spec["dataset_key"],
                    spec["ebeam"],
                    spec["xB"],
                    spec["Q2"],
                    spec["t_abs"],
                )
                specs_by_key.setdefault(key, spec)
            #endfor
        #endif
    #endif

    # Pass-2 anchor panels.  Fitted variants reuse the same model curves.
    p2_summary = make_pass2_anchor_panel_summary(matches)
    if not p2_summary.empty:
        p2 = matches.loc[matches["dataset_b"] == "pass2"].copy()
        p2["anchor_group"] = [
            f"bin_{int(round(v))}" for v in p2["published_bin_b"].to_numpy(float)
        ]
        for r in p2_summary.itertuples(index=False):
            d = p2.loc[p2["anchor_group"] == str(r.anchor_group)]
            if d.empty:
                continue
            #endif
            p2pts = d.drop_duplicates("point_id_b")
            spec = {
                "dataset_key": "pass2",
                "ebeam": float(np.median(p2pts["ebeam_b"])),
                "xB": float(np.median(p2pts["xB_b"])),
                "Q2": float(np.median(p2pts["Q2_b"])),
                "t_abs": float(np.median(p2pts["t_abs_b"])),
                "source": f"pass2:{r.anchor_group}",
            }
            key = _model_curve_cache_key(
                spec["dataset_key"], spec["ebeam"], spec["xB"], spec["Q2"], spec["t_abs"]
            )
            specs_by_key.setdefault(key, spec)
        #endfor
    #endif

    return list(specs_by_key.values())
#enddef


def precompute_presentation_model_curves(
        emff,
        cache: Dict[Tuple, pd.DataFrame],
        specs: Sequence[Dict[str, float]]) -> None:
    """
    Precompute all unique dense BH/KM15 presentation curves with a true global
    percentage-complete terminal diagnostic.

    The progress percentage counts completed fixed-kinematics phi scans.  Each
    scan itself contains the command-line-selected number of phi points and
    always spans 0--360 degrees inclusive.
    """
    total = int(len(specs))
    if total <= 0:
        print("[MODEL PREDICTIONS] no presentation curves requested", flush=True)
        return
    #endif

    nphi = int(math.ceil(360.0 / float(MODEL_CURVE_PHI_STEP_DEG))) + 1
    print(
        f"[MODEL PREDICTIONS] precomputing {total} unique BH/KM15 curves; "
        f"{nphi} phi points/curve at {MODEL_CURVE_PHI_STEP_DEG:g} deg spacing",
        flush=True,
    )

    for i, spec in enumerate(specs, start=1):
        pct_before = 100.0 * float(i - 1) / float(total)
        print(
            f"[MODEL PREDICTIONS] {i}/{total} "
            f"({pct_before:5.1f}% complete): "
            f"{DATASET_LABELS.get(str(spec['dataset_key']), str(spec['dataset_key']))}, "
            f"E={float(spec['ebeam']):.3f}, "
            f"xB={float(spec['xB']):.3f}, "
            f"Q2={float(spec['Q2']):.3f}, "
            f"|t|={float(spec['t_abs']):.3f}",
            flush=True,
        )

        get_dense_model_curve(
            emff,
            cache,
            dataset_key=str(spec["dataset_key"]),
            ebeam=float(spec["ebeam"]),
            xB=float(spec["xB"]),
            Q2=float(spec["Q2"]),
            t_abs=float(spec["t_abs"]),
        )

        pct_after = 100.0 * float(i) / float(total)
        print(
            f"[MODEL PREDICTIONS] completed {i}/{total} "
            f"({pct_after:5.1f}%)",
            flush=True,
        )
    #endfor

    print(
        f"[MODEL PREDICTIONS] complete: {total}/{total} curves (100.0%)",
        flush=True,
    )
#enddef


def plot_pairwise_cross_section_panels(
        matches: pd.DataFrame,
        panel_summary: pd.DataFrame,
        outdir: Path,
        emff,
        model_curve_cache: Dict[Tuple, pd.DataFrame],
        max_pages_per_pair: int = PANEL_MAX_PAGES_PER_PAIR) -> None:
    """
    Produce note-style 3x4 phi-dependent cross-section canvases for every
    dataset pair with sufficient matched coverage.

    Dataset A is shown after point-by-point KM15 transport to the exact
    kinematics of dataset B.  Dataset B is shown at its measured kinematics.
    Smooth BH and KM15 curves are evaluated on a dense 0--360 degree grid at
    one fixed representative B kinematic point per panel.
    """
    if matches.empty or panel_summary.empty:
        return
    #endif

    outdir.mkdir(parents=True, exist_ok=True)

    for (ka, kb), summary_pair in panel_summary.groupby(["dataset_a", "dataset_b"], sort=False):
        pair = matches.loc[
            (matches["dataset_a"] == ka) & (matches["dataset_b"] == kb)
        ].copy()
        pair = _cluster_reference_kinematic_cells(pair)

        ranked = summary_pair.sort_values(
            ["panel_rank_score", "N_matches", "phi_coverage_deg"],
            ascending=False,
        ).reset_index(drop=True)

        max_panels = min(len(ranked), int(max_pages_per_pair) * PANEL_PER_PAGE)
        ranked = ranked.iloc[:max_panels].copy()
        if ranked.empty:
            continue
        #endif

        # After selecting the most informative cells, order them physically by
        # xB, then Q2, then |t| so successive panels/canvases are predictable.
        ranked = ranked.sort_values(
            ["xB_ref", "Q2_ref", "t_abs_ref", "panel_group"],
            ascending=True,
        ).reset_index(drop=True)

        npages = int(math.ceil(len(ranked) / PANEL_PER_PAGE))
        for ipage in range(npages):
            page = ranked.iloc[
                ipage * PANEL_PER_PAGE:(ipage + 1) * PANEL_PER_PAGE
            ]
            fig, axes = plt.subplots(
                PANEL_NROWS, PANEL_NCOLS,
                figsize=(15.8, 10.8),
                squeeze=False,
            )

            for iax, ax in enumerate(axes.ravel()):
                if iax >= len(page):
                    ax.axis("off")
                    continue
                #endif

                group = str(page.iloc[iax]["panel_group"])
                d = pair.loc[pair["panel_group"] == group].copy()
                _plot_one_pair_cross_section_panel(
                    ax, d, ka, kb,
                    emff,
                    model_curve_cache,
                    show_legend_labels=(iax == 0),
                )

                row = iax // PANEL_NCOLS
                col = iax % PANEL_NCOLS
                if row == PANEL_NROWS - 1:
                    ax.set_xlabel(r"$\phi$ (deg)")
                #endif
                if col == 0:
                    ax.set_ylabel(
                        r"$d^4\sigma/(dQ^2\,dx_B\,d|t|\,d\phi)$ (nb/GeV$^4$)",
                        fontsize=8.5,
                    )
                #endif
            #endfor

            handles, labels = axes.ravel()[0].get_legend_handles_labels()
            fig.suptitle(
                f"{DATASET_LABELS[ka]} vs {DATASET_LABELS[kb]}",
                y=0.992, fontsize=14,
            )
            if handles:
                fig.legend(
                    handles, labels,
                    loc="upper center",
                    bbox_to_anchor=(0.5, 0.958),
                    ncol=4,
                    frameon=False,
                    fontsize=8.2,
                )
            #endif
            fig.text(
                0.5, 0.922,
                (
                    f"{DATASET_LABELS[ka]} transported point-by-point to "
                    f"{DATASET_LABELS[kb]} kinematics with KM15.  "
                    r"Error bars = stat $\oplus$ point-to-point syst.; "
                    "open markers have >100% point uncertainty."
                ),
                ha="center", va="top", fontsize=8.0,
            )
            # Reserve a deliberately larger top margin than the previous
            # version; this prevents the title/legend/subtitle collision.
            fig.tight_layout(rect=[0.035, 0.035, 0.995, 0.885])

            fname = (
                f"pair_{ka}_vs_{kb}_cross_sections_page{ipage + 1:02d}.png"
            )
            fig.savefig(outdir / fname, dpi=220)
            plt.close(fig)
        #endfor
    #endfor
#enddef

def make_lee_anchor_panel_summary(matches: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize Lee-centered panels containing one or more external measurements
    transported to the exact CLAS12 pass-1 point kinematics.
    """
    lee = matches.loc[matches["dataset_b"] == "lee2026"].copy()
    if lee.empty:
        return pd.DataFrame()
    #endif

    # Lee has an authoritative published bin index.  Fall back to the generic
    # reference clustering only if that index is unexpectedly unavailable.
    if "published_bin_b" in lee.columns and np.isfinite(lee["published_bin_b"]).any():
        lee["anchor_group"] = [
            (
                f"bin_{int(round(v))}"
                if np.isfinite(v) else f"point_{pid}"
            )
            for v, pid in zip(
                lee["published_bin_b"].to_numpy(float),
                lee["point_id_b"].astype(str),
            )
        ]
    else:
        pieces = []
        for (ka, kb), d in lee.groupby(["dataset_a", "dataset_b"], sort=False):
            dd = _cluster_reference_kinematic_cells(d)
            dd["anchor_group"] = dd["panel_group"]
            pieces.append(dd)
        #endfor
        lee = pd.concat(pieces, ignore_index=True)
    #endif

    rows = []
    for group, d in lee.groupby("anchor_group", sort=False):
        unique_lee_points = d.drop_duplicates("point_id_b")
        n_external = int(d["dataset_a"].nunique())
        n_external_points = int(len(d))
        n_lee_points = int(len(unique_lee_points))

        # Require enough phi information to make a meaningful cross-section
        # shape panel.
        phi = np.sort(np.mod(unique_lee_points["phi_b"].to_numpy(float), 360.0))
        if len(phi) < PANEL_MIN_MATCHES:
            continue
        #endif
        gaps = np.diff(np.r_[phi, phi[0] + 360.0])
        coverage = 360.0 - float(np.max(gaps)) if len(phi) >= 2 else 0.0

        rows.append({
            "anchor_group": group,
            "published_bin": (
                float(np.nanmedian(d["published_bin_b"]))
                if np.isfinite(d["published_bin_b"]).any() else np.nan
            ),
            "N_external_datasets": n_external,
            "N_external_matches": n_external_points,
            "N_lee_points": n_lee_points,
            "xB_ref": float(np.median(d["xB_b"])),
            "Q2_ref": float(np.median(d["Q2_b"])),
            "t_abs_ref": float(np.median(d["t_abs_b"])),
            "xBmin_ref": float(np.nanmedian(d["xBmin_b"])) if "xBmin_b" in d and np.isfinite(d["xBmin_b"]).any() else np.nan,
            "xBmax_ref": float(np.nanmedian(d["xBmax_b"])) if "xBmax_b" in d and np.isfinite(d["xBmax_b"]).any() else np.nan,
            "Q2min_ref": float(np.nanmedian(d["Q2min_b"])) if "Q2min_b" in d and np.isfinite(d["Q2min_b"]).any() else np.nan,
            "Q2max_ref": float(np.nanmedian(d["Q2max_b"])) if "Q2max_b" in d and np.isfinite(d["Q2max_b"]).any() else np.nan,
            "t_abs_min_ref": float(np.nanmedian(d["t_abs_min_b"])) if "t_abs_min_b" in d and np.isfinite(d["t_abs_min_b"]).any() else np.nan,
            "t_abs_max_ref": float(np.nanmedian(d["t_abs_max_b"])) if "t_abs_max_b" in d and np.isfinite(d["t_abs_max_b"]).any() else np.nan,
            "phi_coverage_deg": coverage,
            "datasets": ",".join(sorted(d["dataset_a"].unique())),
            "panel_rank_score": float(
                3.0 * n_external + n_lee_points + coverage / 360.0
            ),
        })
    #endfor

    return pd.DataFrame(rows)
#enddef



def _build_lee_anchor_observations(matches: pd.DataFrame) -> pd.DataFrame:
    """
    Build one observation table for the global Lee-anchor normalization fit.

    Each external point contributes its KM15-transported value at the exact Lee
    point to which it was matched.  Each Lee point is included exactly once,
    even if several external datasets matched to it.
    """
    lee = matches.loc[matches["dataset_b"] == "lee2026"].copy()
    if lee.empty:
        return pd.DataFrame()
    #endif

    rows = []

    lee_unique = lee.sort_values("point_id_b").drop_duplicates("point_id_b")
    for r in lee_unique.itertuples(index=False):
        rows.append({
            "anchor_id": str(r.point_id_b),
            "dataset": "lee2026",
            "point_id": str(r.point_id_b),
            "published_bin": float(r.published_bin_b),
            "phi_deg": float(r.phi_b),
            "t_abs": float(r.t_abs_b),
            "value": float(r.xs_b),
            "unc": float(r.point_unc_b),
            "norm_frac": float(r.norm_frac_b),
        })
    #endfor

    for r in lee.itertuples(index=False):
        rows.append({
            "anchor_id": str(r.point_id_b),
            "dataset": str(r.dataset_a),
            "point_id": str(r.point_id_a),
            "published_bin": float(r.published_bin_b),
            "phi_deg": float(r.phi_b),
            "t_abs": float(r.t_abs_a),
            "value": float(r.xs_a_to_b_km15),
            "unc": float(r.point_unc_a_to_b_km15),
            "norm_frac": float(r.norm_frac_a),
        })
    #endfor

    out = pd.DataFrame(rows)
    finite = (
        np.isfinite(out["value"].to_numpy(float))
        & np.isfinite(out["unc"].to_numpy(float))
        & (out["value"].to_numpy(float) > 0.0)
        & (out["unc"].to_numpy(float) > 0.0)
    )
    return out.loc[finite].reset_index(drop=True)
#enddef


def fit_global_lee_anchor_normalizations(
        matches: pd.DataFrame,
        *,
        exclude_datasets: Sequence[str] = (),
        exclude_point_ids: Sequence[str] = (),
        saylor_tmin: Optional[float] = None,
        scenario: str = "nominal") -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    """
    Fit experiment-wide multiplicative normalization offsets using all matched
    measurements at the Lee anchors.

    The fit is intentionally independent of KM15 *after* the external points
    have been transported to the Lee kinematics.  At each matched Lee point j,
    an unconstrained common cross section mu_j is introduced.  In log space,

        log(y_dj) = log(mu_j) + eta_d,

    where eta_d is one global normalization offset for dataset d.

    For each constrained experiment eta_d has a Gaussian prior centered on zero
    with width log(1+n_d), where n_d is its quoted correlated normalization.
    This means a +1-sigma multiplicative excursion is exactly (1+n_d).
    Georges is deliberately left free.

    Because the system is linear in log(mu_j) and eta_d, the complete global
    solution is obtained with one weighted linear least-squares solve; no
    numerical minimizer is needed.

    The reported correction applied to plotted data is exp(-eta_d).  Thus a
    correction of +2% means that dataset's displayed cross sections are
    multiplied by 1.02 in the normalization-adjusted figure.

    The total chi2 includes both measurement residuals and Gaussian nuisance
    penalties.  The reported ndf counts the Gaussian priors as independent
    constraints:

        ndf = N_data + N_priors - N_anchor_values - N_dataset_offsets.
    """
    obs = _build_lee_anchor_observations(matches)
    if obs.empty:
        return pd.DataFrame(), pd.DataFrame(), {}
    #endif

    exclude_dataset_set = {str(x) for x in exclude_datasets}
    exclude_point_set = {str(x) for x in exclude_point_ids}
    if exclude_dataset_set:
        obs = obs.loc[~obs["dataset"].isin(exclude_dataset_set)].copy()
    #endif
    if exclude_point_set:
        obs = obs.loc[~obs["point_id"].isin(exclude_point_set)].copy()
    #endif

    if saylor_tmin is not None:
        threshold = float(saylor_tmin)
        is_saylor = obs["dataset"].astype(str) == "saylor2018"
        before = int(np.sum(is_saylor))
        keep = (~is_saylor) | (obs["t_abs"].to_numpy(float) >= threshold)
        obs = obs.loc[keep].copy()
        after = int(np.sum(obs["dataset"].astype(str) == "saylor2018"))
        print(
            f"[NORM FIT] {scenario}: Saylor |t| >= {threshold:.3f} GeV^2 "
            f"retains {after}/{before} matched Saylor observations",
            flush=True,
        )
    #endif

    # A common anchor cross section is only identifiable/useful if at least two
    # different experiments remain at that point.
    counts = obs.groupby("anchor_id")["dataset"].nunique()
    valid_anchors = counts.loc[counts >= 2].index
    obs = obs.loc[obs["anchor_id"].isin(valid_anchors)].copy().reset_index(drop=True)
    if obs.empty:
        return pd.DataFrame(), pd.DataFrame(), {}
    #endif

    anchors = sorted(obs["anchor_id"].astype(str).unique())
    datasets = [key for key in DATASET_ORDER if key in set(obs["dataset"].astype(str))]

    n_anchor = len(anchors)
    n_dataset = len(datasets)
    anchor_index = {key: i for i, key in enumerate(anchors)}
    dataset_index = {key: n_anchor + i for i, key in enumerate(datasets)}
    npar = n_anchor + n_dataset

    design_rows = []
    rhs = []
    data_meta = []

    for r in obs.itertuples(index=False):
        frac_unc = float(r.unc) / float(r.value)
        if not np.isfinite(frac_unc) or frac_unc <= 0.0:
            continue
        #endif

        row = np.zeros(npar, dtype=float)
        weight = 1.0 / frac_unc
        row[anchor_index[str(r.anchor_id)]] = weight
        row[dataset_index[str(r.dataset)]] = weight
        design_rows.append(row)
        rhs.append(math.log(float(r.value)) * weight)
        data_meta.append({
            "anchor_id": str(r.anchor_id),
            "dataset": str(r.dataset),
            "point_id": str(r.point_id),
            "published_bin": float(r.published_bin),
            "phi_deg": float(r.phi_deg),
            "t_abs": float(r.t_abs),
            "value": float(r.value),
            "unc": float(r.unc),
            "frac_unc": frac_unc,
        })
    #endfor

    n_data = len(design_rows)
    n_prior = 0
    prior_datasets = []

    for dataset in datasets:
        if dataset in GLOBAL_NORM_FREE_DATASETS:
            continue
        #endif

        drows = obs.loc[obs["dataset"] == dataset]
        if drows.empty:
            continue
        #endif

        norm_frac = float(np.nanmedian(drows["norm_frac"].to_numpy(float)))
        sigma_log = math.log1p(norm_frac)
        if not np.isfinite(sigma_log) or sigma_log <= 0.0:
            continue
        #endif

        row = np.zeros(npar, dtype=float)
        row[dataset_index[dataset]] = 1.0 / sigma_log
        design_rows.append(row)
        rhs.append(0.0)
        n_prior += 1
        prior_datasets.append(dataset)
    #endfor

    A = np.asarray(design_rows, dtype=float)
    b = np.asarray(rhs, dtype=float)
    solution, _, rank, _ = np.linalg.lstsq(A, b, rcond=None)

    residual = A @ solution - b
    data_residual = residual[:n_data]
    prior_residual = residual[n_data:]

    chi2_data = float(np.sum(data_residual**2))
    chi2_prior = float(np.sum(prior_residual**2))
    chi2_total = chi2_data + chi2_prior
    ndf = int(n_data + n_prior - npar)

    # Point-level residual table.
    point_rows = []
    for meta, pull in zip(data_meta, data_residual):
        eta = float(solution[dataset_index[meta["dataset"]]])
        correction = math.exp(-eta)
        point_rows.append({
            **meta,
            "scenario": str(scenario),
            "eta_dataset": eta,
            "normalization_correction": correction,
            "normalization_correction_pct": 100.0 * (correction - 1.0),
            "adjusted_value": correction * float(meta["value"]),
            "adjusted_unc": correction * float(meta["unc"]),
            "fit_pull_log": float(pull),
        })
    #endfor

    point_table = pd.DataFrame(point_rows)

    # Dataset-level normalization summary.
    dataset_rows = []
    for dataset in datasets:
        eta = float(solution[dataset_index[dataset]])
        correction = math.exp(-eta)

        drows = obs.loc[obs["dataset"] == dataset]
        norm_frac = float(np.nanmedian(drows["norm_frac"].to_numpy(float)))
        sigma_log = math.log1p(norm_frac)
        beta = (
            eta / sigma_log
            if dataset not in GLOBAL_NORM_FREE_DATASETS and sigma_log > 0.0
            else np.nan
        )

        dpoint = point_table.loc[point_table["dataset"] == dataset]
        dataset_rows.append({
            "scenario": str(scenario),
            "dataset": dataset,
            "dataset_label": DATASET_LABELS[dataset],
            "N_points": int(len(dpoint)),
            "quoted_norm_pct": 100.0 * norm_frac,
            "normalization_constraint": (
                "free" if dataset in GLOBAL_NORM_FREE_DATASETS else "Gaussian"
            ),
            "eta_dataset": eta,
            "beta_prior_sigma": beta,
            "data_correction_scale": correction,
            "data_correction_pct": 100.0 * (correction - 1.0),
            "dataset_pull_rms": (
                float(np.sqrt(np.mean(dpoint["fit_pull_log"]**2)))
                if not dpoint.empty else np.nan
            ),
        })
    #endfor

    dataset_table = pd.DataFrame(dataset_rows)

    metrics = {
        "scenario": str(scenario),
        "N_data": int(n_data),
        "N_anchors": int(n_anchor),
        "N_datasets": int(n_dataset),
        "N_priors": int(n_prior),
        "matrix_rank": int(rank),
        "chi2_data": chi2_data,
        "chi2_prior": chi2_prior,
        "chi2_total": chi2_total,
        "ndf": int(ndf),
        "chi2_per_ndf": (chi2_total / ndf if ndf > 0 else np.nan),
    }
    return dataset_table, point_table, metrics
#enddef


def run_global_lee_anchor_normalization_scenarios(
        matches: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Nominal:
      the two invalid Saylor bin-87 points are already removed upstream.

    Diagnostics:
      * additionally require Saylor |t| >= 0.343 GeV^2;
      * remove Saylor entirely.
    """
    specs = [
        ("nominal", (), (), None),
        ("saylor_tmin_0p343", (), (), SAYLOR_TMIN_DIAGNOSTIC_GEV2),
        ("without_saylor", ("saylor2018",), (), None),
    ]

    dataset_tables = []
    point_tables = []
    metric_rows = []

    print(f"[NORM FIT] running {len(specs)} global scenarios", flush=True)

    for ispec, (scenario, excluded_datasets, excluded_points, saylor_tmin) in enumerate(specs, start=1):
        print(
            f"[NORM FIT] scenario {ispec}/{len(specs)}: "
            f"{GLOBAL_NORM_SCENARIO_LABELS.get(scenario, scenario)}",
            flush=True,
        )

        dtab, ptab, metrics = fit_global_lee_anchor_normalizations(
            matches,
            exclude_datasets=excluded_datasets,
            exclude_point_ids=excluded_points,
            saylor_tmin=saylor_tmin,
            scenario=scenario,
        )
        if not dtab.empty:
            dataset_tables.append(dtab)
        #endif
        if not ptab.empty:
            point_tables.append(ptab)
        #endif
        if metrics:
            metrics = dict(metrics)
            metrics["excluded_datasets"] = ",".join(excluded_datasets)
            metrics["excluded_point_ids"] = ",".join(excluded_points)
            metrics["saylor_tmin_GeV2"] = (
                float(saylor_tmin) if saylor_tmin is not None else np.nan
            )
            metric_rows.append(metrics)
            print(
                f"[NORM FIT] {scenario}: Ndata={metrics['N_data']}, "
                f"Nanchors={metrics['N_anchors']}, "
                f"chi2/ndf={metrics['chi2_per_ndf']:.4f}",
                flush=True,
            )
        #endif
    #endfor

    return (
        pd.concat(dataset_tables, ignore_index=True) if dataset_tables else pd.DataFrame(),
        pd.concat(point_tables, ignore_index=True) if point_tables else pd.DataFrame(),
        pd.DataFrame(metric_rows),
    )
#enddef



def make_pass2_anchor_panel_summary(matches: pd.DataFrame) -> pd.DataFrame:
    """Summarize external measurements transported to exact pass-2 kinematics."""
    p2 = matches.loc[matches["dataset_b"] == "pass2"].copy()
    if p2.empty:
        return pd.DataFrame()
    #endif

    p2["anchor_group"] = [
        (
            f"bin_{int(round(v))}"
            if np.isfinite(v) else f"point_{pid}"
        )
        for v, pid in zip(
            p2["published_bin_b"].to_numpy(float),
            p2["point_id_b"].astype(str),
        )
    ]

    rows = []
    for group, d in p2.groupby("anchor_group", sort=False):
        unique_p2 = d.drop_duplicates("point_id_b")
        phi = np.sort(np.mod(unique_p2["phi_b"].to_numpy(float), 360.0))
        if len(phi) < PANEL_MIN_MATCHES:
            continue
        #endif
        gaps = np.diff(np.r_[phi, phi[0] + 360.0])
        coverage = 360.0 - float(np.max(gaps)) if len(phi) >= 2 else 0.0

        rows.append({
            "anchor_group": group,
            "published_bin": float(np.nanmedian(d["published_bin_b"])),
            "N_external_datasets": int(d["dataset_a"].nunique()),
            "N_external_matches": int(len(d)),
            "N_pass2_points": int(len(unique_p2)),
            "xB_ref": float(np.median(d["xB_b"])),
            "Q2_ref": float(np.median(d["Q2_b"])),
            "t_abs_ref": float(np.median(d["t_abs_b"])),
            "xBmin_ref": float(np.nanmedian(d["xBmin_b"])) if "xBmin_b" in d and np.isfinite(d["xBmin_b"]).any() else np.nan,
            "xBmax_ref": float(np.nanmedian(d["xBmax_b"])) if "xBmax_b" in d and np.isfinite(d["xBmax_b"]).any() else np.nan,
            "Q2min_ref": float(np.nanmedian(d["Q2min_b"])) if "Q2min_b" in d and np.isfinite(d["Q2min_b"]).any() else np.nan,
            "Q2max_ref": float(np.nanmedian(d["Q2max_b"])) if "Q2max_b" in d and np.isfinite(d["Q2max_b"]).any() else np.nan,
            "t_abs_min_ref": float(np.nanmedian(d["t_abs_min_b"])) if "t_abs_min_b" in d and np.isfinite(d["t_abs_min_b"]).any() else np.nan,
            "t_abs_max_ref": float(np.nanmedian(d["t_abs_max_b"])) if "t_abs_max_b" in d and np.isfinite(d["t_abs_max_b"]).any() else np.nan,
            "phi_coverage_deg": coverage,
            "datasets": ",".join(sorted(d["dataset_a"].unique())),
            "panel_rank_score": float(
                3.0 * d["dataset_a"].nunique() + len(unique_p2) + coverage / 360.0
            ),
        })
    #endfor
    return pd.DataFrame(rows)
#enddef


def _build_pass2_anchor_observations(matches: pd.DataFrame) -> pd.DataFrame:
    """
    Build the observation table for pass-2-centered global consistency fits.

    External points are already transported to exact pass-2 kinematics by KM15.
    Pass-2 points are included once per anchor and retain their correlated-scale
    response magnitudes.
    """
    p2 = matches.loc[matches["dataset_b"] == "pass2"].copy()
    if p2.empty:
        return pd.DataFrame()
    #endif

    rows = []
    p2_unique = p2.sort_values("point_id_b").drop_duplicates("point_id_b")
    for r in p2_unique.itertuples(index=False):
        rows.append({
            "anchor_id": str(r.point_id_b),
            "dataset": "pass2",
            "point_id": str(r.point_id_b),
            "published_bin": float(r.published_bin_b),
            "phi_deg": float(r.phi_b),
            "t_abs": float(r.t_abs_b),
            "value": float(r.xs_b),
            "unc": float(r.point_unc_b),
            "norm_frac": float(r.norm_frac_b),
            "corr_scale_frac": float(r.corr_scale_frac_b),
        })
    #endfor

    for r in p2.itertuples(index=False):
        rows.append({
            "anchor_id": str(r.point_id_b),
            "dataset": str(r.dataset_a),
            "point_id": str(r.point_id_a),
            "published_bin": float(r.published_bin_b),
            "phi_deg": float(r.phi_b),
            "t_abs": float(r.t_abs_a),
            "value": float(r.xs_a_to_b_km15),
            "unc": float(r.point_unc_a_to_b_km15),
            "norm_frac": float(r.norm_frac_a),
            "corr_scale_frac": 0.0,
        })
    #endfor

    out = pd.DataFrame(rows)
    good = (
        np.isfinite(out["value"].to_numpy(float))
        & np.isfinite(out["unc"].to_numpy(float))
        & (out["value"].to_numpy(float) > 0.0)
        & (out["unc"].to_numpy(float) > 0.0)
    )
    return out.loc[good].reset_index(drop=True)
#enddef


def fit_pass2_anchor_nuisances(
        matches: pd.DataFrame,
        *,
        scenario: str,
        include_pass2_correlated_scale: bool,
        exclude_datasets: Sequence[str] = (),
        saylor_tmin: Optional[float] = None
        ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    """
    Global pass-2-centered consistency fit using the FINAL publication-level
    systematic categories.

    At each exact pass-2 anchor j there is a free common cross section mu_j.
    Each experiment has one experiment-wide normalization offset eta_d.

    Pass-2 therefore has exactly TWO systematic nuisance directions relevant
    here:

      1. eta_pass2:
         one overall normalization nuisance with the finalized 2.16% Gaussian
         prior;

      2. beta_corr:
         one Gaussian nuisance multiplying the finalized bin-dependent
         correlated-scale response

             beta_corr * f_corr,i ,

         where f_corr,i is read directly from
         'correlated scale sys frac, 10.6 GeV'.

    No attempt is made to decompose f_corr,i back into current-dependent and
    run-period pieces.  That lower-level decomposition was already combined in
    the finalized analysis and is not the publication covariance model.

    The fit is done in log cross section so multiplicative nuisance shifts are
    handled naturally.
    """
    obs = _build_pass2_anchor_observations(matches)
    if obs.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {}
    #endif

    exclude_set = {str(x) for x in exclude_datasets}
    if exclude_set:
        obs = obs.loc[~obs["dataset"].isin(exclude_set)].copy()
    #endif

    if saylor_tmin is not None:
        keep = (
            (obs["dataset"].astype(str) != "saylor2018")
            | (obs["t_abs"].to_numpy(float) >= float(saylor_tmin))
        )
        obs = obs.loc[keep].copy()
    #endif

    if include_pass2_correlated_scale:
        is_pass2 = obs["dataset"].astype(str) == "pass2"
        missing_corr = is_pass2 & ~np.isfinite(
            obs["corr_scale_frac"].to_numpy(float)
        )
        n_missing_corr = int(np.sum(missing_corr))
        if n_missing_corr:
            print(
                f"[PASS2 FIT] {scenario}: excluding {n_missing_corr} pass-2 "
                "observation(s) lacking a finalized correlated-scale response",
                flush=True,
            )
            obs = obs.loc[~missing_corr].copy()
        #endif
    #endif

    counts = obs.groupby("anchor_id")["dataset"].nunique()
    anchors = sorted(
        counts.loc[counts >= 2].index.astype(str)
    )
    obs = obs.loc[
        obs["anchor_id"].astype(str).isin(anchors)
    ].copy().reset_index(drop=True)

    if obs.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {}
    #endif

    datasets = [
        key for key in DATASET_ORDER
        if key in set(obs["dataset"].astype(str))
    ]

    param_names = (
        [f"mu::{anchor}" for anchor in anchors]
        + [f"eta::{dataset}" for dataset in datasets]
    )
    if include_pass2_correlated_scale:
        param_names.append("beta::pass2_corr_scale")
    #endif

    pindex = {name: i for i, name in enumerate(param_names)}

    Arows = []
    brows = []
    meta = []

    for r in obs.itertuples(index=False):
        frac_unc = float(r.unc) / float(r.value)
        if (
            not np.isfinite(frac_unc)
            or frac_unc <= 0.0
        ):
            continue
        #endif

        weight = 1.0 / frac_unc
        row = np.zeros(len(param_names), dtype=float)

        row[pindex[f"mu::{r.anchor_id}"]] = weight
        row[pindex[f"eta::{r.dataset}"]] = weight

        if (
            include_pass2_correlated_scale
            and str(r.dataset) == "pass2"
        ):
            row[pindex["beta::pass2_corr_scale"]] = (
                weight * float(r.corr_scale_frac)
            )
        #endif

        Arows.append(row)
        brows.append(math.log(float(r.value)) * weight)
        meta.append(r._asdict())
    #endfor

    n_data = len(Arows)
    n_prior = 0

    # One normalization prior per constrained experiment.  Georges remains
    # intentionally free, consistent with the external-world study.
    for dataset in datasets:
        if dataset in GLOBAL_NORM_FREE_DATASETS:
            continue
        #endif

        vals = obs.loc[
            obs["dataset"] == dataset,
            "norm_frac",
        ].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            continue
        #endif

        sigma_log = math.log1p(float(np.nanmedian(vals)))
        if sigma_log <= 0.0:
            continue
        #endif

        row = np.zeros(len(param_names), dtype=float)
        row[pindex[f"eta::{dataset}"]] = 1.0 / sigma_log
        Arows.append(row)
        brows.append(0.0)
        n_prior += 1
    #endfor

    # One publication-level pass-2 correlated-scale nuisance.
    if include_pass2_correlated_scale:
        row = np.zeros(len(param_names), dtype=float)
        row[pindex["beta::pass2_corr_scale"]] = 1.0
        Arows.append(row)
        brows.append(0.0)
        n_prior += 1
    #endif

    A = np.asarray(Arows, dtype=float)
    b = np.asarray(brows, dtype=float)

    solution, _, rank, _ = np.linalg.lstsq(A, b, rcond=None)
    residual = A @ solution - b

    data_residual = residual[:n_data]
    prior_residual = residual[n_data:]

    chi2_data = float(np.sum(data_residual**2))
    chi2_prior = float(np.sum(prior_residual**2))
    chi2_total = chi2_data + chi2_prior

    ndf = int(
        n_data
        + n_prior
        - len(param_names)
    )

    eta_by_dataset = {
        dataset: float(solution[pindex[f"eta::{dataset}"]])
        for dataset in datasets
    }

    beta_corr = (
        float(solution[pindex["beta::pass2_corr_scale"]])
        if include_pass2_correlated_scale
        else 0.0
    )

    point_rows = []
    for metadata, pull in zip(meta, data_residual):
        dataset = str(metadata["dataset"])
        eta = eta_by_dataset[dataset]

        corr_shift = 0.0
        if (
            dataset == "pass2"
            and include_pass2_correlated_scale
        ):
            corr_shift = (
                beta_corr
                * float(metadata["corr_scale_frac"])
            )
        #endif

        # The fitted model is log(y) = log(mu) + eta + beta*f.
        # To display measurements aligned to the common anchor convention,
        # apply the inverse multiplicative nuisance shift.
        correction = math.exp(-(eta + corr_shift))

        point_rows.append({
            **metadata,
            "scenario": scenario,
            "normalization_eta": eta,
            "pass2_corr_log_shift": corr_shift,
            "data_correction_scale": correction,
            "data_correction_pct": 100.0 * (correction - 1.0),
            "adjusted_value": correction * float(metadata["value"]),
            "adjusted_unc": correction * float(metadata["unc"]),
            "fit_pull_log": float(pull),
        })
    #endfor

    point_table = pd.DataFrame(point_rows)

    dataset_rows = []
    for dataset in datasets:
        eta = eta_by_dataset[dataset]

        vals = obs.loc[
            obs["dataset"] == dataset,
            "norm_frac",
        ].to_numpy(float)
        vals = vals[np.isfinite(vals)]

        norm_frac = (
            float(np.nanmedian(vals))
            if vals.size
            else np.nan
        )

        sigma_log = (
            math.log1p(norm_frac)
            if np.isfinite(norm_frac) and norm_frac >= 0.0
            else np.nan
        )

        beta_norm = (
            eta / sigma_log
            if (
                dataset not in GLOBAL_NORM_FREE_DATASETS
                and np.isfinite(sigma_log)
                and sigma_log > 0.0
            )
            else np.nan
        )

        dataset_rows.append({
            "scenario": scenario,
            "dataset": dataset,
            "dataset_label": DATASET_LABELS[dataset],
            "N_points": int(
                np.sum(point_table["dataset"] == dataset)
            ),
            "quoted_norm_pct": (
                100.0 * norm_frac
                if np.isfinite(norm_frac)
                else np.nan
            ),
            "normalization_constraint": (
                "free"
                if dataset in GLOBAL_NORM_FREE_DATASETS
                else "Gaussian"
            ),
            "normalization_eta": eta,
            "beta_norm": beta_norm,
            "global_normalization_correction": math.exp(-eta),
            "global_normalization_correction_pct": (
                100.0 * (math.exp(-eta) - 1.0)
            ),
        })
    #endfor

    dataset_table = pd.DataFrame(dataset_rows)

    nuisance_table = pd.DataFrame([{
        "scenario": scenario,
        "nuisance": "pass2_correlated_scale",
        "beta_sigma": beta_corr,
        "enabled": bool(include_pass2_correlated_scale),
    }])

    # Publication-level Hayward nuisance summary.  Keep the two nuisance
    # amplitudes separate, and also report their radial excursion in the
    # two-dimensional independent-Gaussian nuisance space.
    pass2_eta = float(eta_by_dataset.get("pass2", 0.0))
    pass2_norm_vals = obs.loc[
        obs["dataset"] == "pass2",
        "norm_frac",
    ].to_numpy(float)
    pass2_norm_vals = pass2_norm_vals[np.isfinite(pass2_norm_vals)]
    pass2_norm_frac = (
        float(np.nanmedian(pass2_norm_vals))
        if pass2_norm_vals.size else np.nan
    )
    pass2_norm_sigma_log = (
        math.log1p(pass2_norm_frac)
        if np.isfinite(pass2_norm_frac) and pass2_norm_frac > 0.0
        else np.nan
    )
    beta_pass2_norm = (
        pass2_eta / pass2_norm_sigma_log
        if np.isfinite(pass2_norm_sigma_log) and pass2_norm_sigma_log > 0.0
        else np.nan
    )

    pass2_points = point_table.loc[
        point_table["dataset"] == "pass2"
    ].copy()
    if not pass2_points.empty:
        total_shift_pct = 100.0 * (
            pass2_points["data_correction_scale"].to_numpy(float) - 1.0
        )
        corr_only_shift_pct = 100.0 * (
            np.exp(
                -pass2_points["pass2_corr_log_shift"].to_numpy(float)
            ) - 1.0
        )
        pass2_total_shift_median_pct = float(np.nanmedian(total_shift_pct))
        pass2_total_shift_min_pct = float(np.nanmin(total_shift_pct))
        pass2_total_shift_max_pct = float(np.nanmax(total_shift_pct))
        pass2_corr_shift_median_pct = float(np.nanmedian(corr_only_shift_pct))
        pass2_corr_shift_min_pct = float(np.nanmin(corr_only_shift_pct))
        pass2_corr_shift_max_pct = float(np.nanmax(corr_only_shift_pct))
    else:
        pass2_total_shift_median_pct = np.nan
        pass2_total_shift_min_pct = np.nan
        pass2_total_shift_max_pct = np.nan
        pass2_corr_shift_median_pct = np.nan
        pass2_corr_shift_min_pct = np.nan
        pass2_corr_shift_max_pct = np.nan
    #endif

    metrics = {
        "scenario": scenario,
        "include_pass2_correlated_scale": bool(
            include_pass2_correlated_scale
        ),
        "N_data": int(n_data),
        "N_anchors": int(len(anchors)),
        "N_datasets": int(len(datasets)),
        "N_priors": int(n_prior),
        "matrix_rank": int(rank),
        "chi2_data": chi2_data,
        "chi2_prior": chi2_prior,
        "chi2_total": chi2_total,
        "ndf": int(ndf),
        "chi2_per_ndf": (
            chi2_total / ndf
            if ndf > 0
            else np.nan
        ),
        "beta_pass2_norm": beta_pass2_norm,
        "beta_pass2_corr_scale": beta_corr,
        "combined_pass2_nuisance_excursion_sigma": float(
            math.sqrt(
                (beta_pass2_norm**2 if np.isfinite(beta_pass2_norm) else 0.0)
                + (beta_corr**2 if include_pass2_correlated_scale else 0.0)
            )
        ),
        "pass2_corr_shift_median_pct": pass2_corr_shift_median_pct,
        "pass2_corr_shift_min_pct": pass2_corr_shift_min_pct,
        "pass2_corr_shift_max_pct": pass2_corr_shift_max_pct,
        "pass2_total_shift_median_pct": pass2_total_shift_median_pct,
        "pass2_total_shift_min_pct": pass2_total_shift_min_pct,
        "pass2_total_shift_max_pct": pass2_total_shift_max_pct,
    }

    return (
        dataset_table,
        nuisance_table,
        point_table,
        metrics,
    )
#enddef

def run_pass2_anchor_nuisance_scenarios(matches: pd.DataFrame):
    """
    Run pass-2 consistency scenarios.

    The key comparison is:
      * norm_only_nominal:
          experiment-wide normalizations only;
      * full_corr_nominal:
          same normalization treatment plus ONE pass-2 correlated-scale
          nuisance using the finalized bin-dependent response.

    Saylor diagnostics mirror the established external-world study.
    """
    specs = [
        (
            "norm_only_nominal",
            False,
            (),
            None,
        ),
        (
            "full_corr_nominal",
            True,
            (),
            None,
        ),
        (
            "full_corr_saylor_tmin_0p343",
            True,
            (),
            SAYLOR_TMIN_DIAGNOSTIC_GEV2,
        ),
        (
            "full_corr_without_saylor",
            True,
            ("saylor2018",),
            None,
        ),
    ]

    ds_all = []
    nuisance_all = []
    point_all = []
    metrics_all = []

    for i, (
            scenario,
            include_corr,
            excluded,
            tmin) in enumerate(specs, start=1):

        print(
            f"[PASS2 FIT] scenario {i}/{len(specs)}: {scenario}",
            flush=True,
        )

        ds, nuisance, points, metrics = fit_pass2_anchor_nuisances(
            matches,
            scenario=scenario,
            include_pass2_correlated_scale=include_corr,
            exclude_datasets=excluded,
            saylor_tmin=tmin,
        )

        if not ds.empty:
            ds_all.append(ds)
        #endif
        if not nuisance.empty:
            nuisance_all.append(nuisance)
        #endif
        if not points.empty:
            point_all.append(points)
        #endif

        if metrics:
            metrics["excluded_datasets"] = ",".join(excluded)
            metrics["saylor_tmin_GeV2"] = (
                float(tmin)
                if tmin is not None
                else np.nan
            )
            metrics_all.append(metrics)

            print(
                f"[PASS2 FIT] {scenario}: "
                f"chi2/ndf={metrics['chi2_per_ndf']:.4f}, "
                f"beta_norm={metrics['beta_pass2_norm']:+.3f}, "
                f"beta_corr={metrics['beta_pass2_corr_scale']:+.3f}, "
                f"r_beta={metrics['combined_pass2_nuisance_excursion_sigma']:.3f}",
                flush=True,
            )
        #endif
    #endfor

    return (
        pd.concat(ds_all, ignore_index=True)
        if ds_all else pd.DataFrame(),
        pd.concat(nuisance_all, ignore_index=True)
        if nuisance_all else pd.DataFrame(),
        pd.concat(point_all, ignore_index=True)
        if point_all else pd.DataFrame(),
        pd.DataFrame(metrics_all),
    )
#enddef

def build_complete_pass2_anchor_legend(
        matches: pd.DataFrame,
        dataset_fit_table: Optional[pd.DataFrame] = None,
        scenario: Optional[str] = None,
        omit_datasets: Sequence[str] = ()) -> Tuple[List[Line2D], List[str]]:
    """Fixed legend for pass-2-centered canvases."""
    omit = {str(x) for x in omit_datasets}
    corr = {k: 1.0 for k in DATASET_ORDER}
    if dataset_fit_table is not None and not dataset_fit_table.empty:
        for r in dataset_fit_table.itertuples(index=False):
            corr[str(r.dataset)] = float(r.global_normalization_correction)
        #endfor
    #endif

    handles, labels = [], []
    for key in DATASET_ORDER:
        style = DATASET_STYLES[key]
        handles.append(Line2D(
            [0], [0],
            linestyle="none",
            marker=style["marker"],
            markersize=5.5,
            markerfacecolor=style["color"],
            markeredgecolor=style["color"],
            color=style["color"],
        ))
        if key in omit:
            label = f"{DATASET_LABELS[key]} (excluded)"
        elif scenario is None:
            # Quote the known overall normalization if available.
            if key == "pass2":
                label = f"{DATASET_LABELS[key]} (2.16% norm + kin. corr.)"
            else:
                # Collect ONLY the normalization values belonging to this
                # dataset itself.  The previous implementation flattened both
                # norm_frac_a and norm_frac_b from every matched row involving
                # the dataset, which mixed in the partner experiment's
                # normalization and produced nonsense labels (e.g. ~18% for
                # Lee instead of the correct 31%).
                vals_a = matches.loc[
                    matches["dataset_a"] == key,
                    "norm_frac_a",
                ].to_numpy(float)
                vals_b = matches.loc[
                    matches["dataset_b"] == key,
                    "norm_frac_b",
                ].to_numpy(float)
                vals = np.concatenate([vals_a, vals_b])
                vals = vals[np.isfinite(vals) & (vals > 0)]
                label = (
                    f"{DATASET_LABELS[key]} ({100*np.nanmedian(vals):.1f}% norm)"
                    if vals.size else DATASET_LABELS[key]
                )
            #endif
        else:
            label = f"{DATASET_LABELS[key]} ({100*(corr.get(key,1.0)-1):+.1f}% norm)"
        #endif
        labels.append(label)
    #endfor

    for mk, mlab in [("bh", "BH"), ("km15", "KM15")]:
        s = MODEL_STYLES[mk]
        handles.append(Line2D(
            [0], [0],
            color=s["color"],
            linestyle=s["linestyle"],
            linewidth=s["linewidth"],
        ))
        labels.append(mlab)
    #endfor
    return handles, labels
#enddef



def _world_canvas_specs_by_xb(selected: pd.DataFrame) -> List[Dict[str, object]]:
    """Arrange one world-data canvas per physical xB bin.

    Rows are fixed (xB,Q2) bins and columns are fixed |t| bins.  Physical bin
    boundaries are used whenever available; representative means are only a
    fallback for legacy inputs lacking boundaries.  Each panel retains its own
    y range.
    """
    if selected.empty:
        return []
    #endif

    work = selected.copy()

    def finite_key(row, lo_name, hi_name, mean_name, digits=6):
        lo = float(row.get(lo_name, np.nan))
        hi = float(row.get(hi_name, np.nan))
        if np.isfinite(lo) and np.isfinite(hi):
            return (round(lo, digits), round(hi, digits))
        #endif
        v = float(row.get(mean_name, np.nan))
        return (round(v, 3), round(v, 3))
    #enddef

    work["_xb_key"] = [
        finite_key(r, "xBmin_ref", "xBmax_ref", "xB_ref")
        for _, r in work.iterrows()
    ]
    work["_q_key"] = [
        finite_key(r, "Q2min_ref", "Q2max_ref", "Q2_ref")
        for _, r in work.iterrows()
    ]
    work["_t_key"] = [
        finite_key(r, "t_abs_min_ref", "t_abs_max_ref", "t_abs_ref")
        for _, r in work.iterrows()
    ]

    specs = []
    for xb_key, xb_group in work.groupby("_xb_key", sort=True):
        q_keys = sorted(set(xb_group["_q_key"]), key=lambda z: (z[0], z[1]))
        t_keys = sorted(set(xb_group["_t_key"]), key=lambda z: (z[0], z[1]))
        q_index = {k: i for i, k in enumerate(q_keys)}
        t_index = {k: i for i, k in enumerate(t_keys)}
        placements = {}
        for _, row in xb_group.sort_values(["Q2_ref", "t_abs_ref", "published_bin"]).iterrows():
            placements[(q_index[row["_q_key"]], t_index[row["_t_key"]])] = row
        #endfor
        specs.append({
            "xb_key": xb_key,
            "q_keys": q_keys,
            "t_keys": t_keys,
            "nrows": max(1, len(q_keys)),
            "ncols": max(1, len(t_keys)),
            "placements": placements,
        })
    #endfor
    return specs
#enddef

def plot_pass2_anchor_world_panels(
        matches: pd.DataFrame,
        panel_summary: pd.DataFrame,
        outdir: Path,
        emff,
        model_curve_cache: Dict[Tuple, pd.DataFrame],
        *,
        dataset_fit_table: Optional[pd.DataFrame] = None,
        point_fit_table: Optional[pd.DataFrame] = None,
        metrics: Optional[Dict[str, float]] = None,
        scenario: Optional[str] = None,
        omit_datasets: Sequence[str] = (),
        saylor_tmin: Optional[float] = None) -> None:
    """
    Pass-2-centered 3x4 canvases.

    Raw version: no systematic rescaling.
    Fitted version: external/global normalization shifts plus the pass-2
    bin-dependent correction implied by its single publication-level correlated-scale nuisance.
    """
    p2 = matches.loc[matches["dataset_b"] == "pass2"].copy()
    if p2.empty or panel_summary.empty:
        return
    #endif

    omit = {str(x) for x in omit_datasets}
    if omit:
        p2 = p2.loc[~p2["dataset_a"].isin(omit)].copy()
    #endif
    if saylor_tmin is not None:
        p2 = p2.loc[
            (p2["dataset_a"] != "saylor2018")
            | (p2["t_abs_a"] >= float(saylor_tmin))
        ].copy()
    #endif

    p2["anchor_group"] = [
        f"bin_{int(round(v))}"
        for v in p2["published_bin_b"].to_numpy(float)
    ]

    selected = panel_summary.loc[
        panel_summary["anchor_group"].isin(set(p2["anchor_group"]))
    ].copy()
    selected = selected.sort_values(
        ["published_bin", "xB_ref", "Q2_ref", "t_abs_ref"]
    ).reset_index(drop=True)
    if selected.empty:
        return
    #endif

    global_corr = {k: 1.0 for k in DATASET_ORDER}
    if dataset_fit_table is not None and not dataset_fit_table.empty:
        for r in dataset_fit_table.itertuples(index=False):
            global_corr[str(r.dataset)] = float(r.global_normalization_correction)
        #endfor
    #endif

    pass2_point_corr = {}
    if point_fit_table is not None and not point_fit_table.empty:
        pp = point_fit_table.loc[point_fit_table["dataset"] == "pass2"]
        pass2_point_corr = dict(zip(
            pp["point_id"].astype(str),
            pp["data_correction_scale"].to_numpy(float),
        ))
    #endif

    offsets = {
        "jo2015": -5.0,
        "defurne2015": -3.5,
        "defurne2017": -2.0,
        "saylor2018": -0.7,
        "georges2022": +0.7,
        "lee2026": +2.0,
    }

    outdir.mkdir(parents=True, exist_ok=True)
    canvas_specs = _world_canvas_specs_by_xb(selected)
    npages = len(canvas_specs)

    scenario_label = "raw" if scenario is None else str(scenario)
    print(
        f"[PLOTS PASS2] {scenario_label}: {len(selected)} panels "
        f"across {npages} xB-organized canvas(es)",
        flush=True,
    )

    for ipage, spec in enumerate(canvas_specs):
        print(
            f"[PLOTS PASS2] {scenario_label}: canvas {ipage + 1}/{npages}",
            flush=True,
        )
        nrows = int(spec["nrows"])
        ncols = int(spec["ncols"])
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(max(8.0, 3.65*ncols), max(5.8, 3.15*nrows + 1.8)),
            squeeze=False,
        )
        placements = spec["placements"]

        for iax, ax in enumerate(axes.ravel()):
            row = iax // ncols
            col = iax % ncols
            panel_row = placements.get((row, col))
            if panel_row is None:
                ax.axis("off")
                continue
            #endif

            group = str(panel_row["anchor_group"])
            d = p2.loc[p2["anchor_group"] == group].copy()
            if d.empty:
                ax.axis("off")
                continue
            #endif

            p2pts = d.sort_values("phi_b").drop_duplicates("point_id_b")
            raw_y = p2pts["xs_b"].to_numpy(float)
            raw_e = p2pts["point_unc_b"].to_numpy(float)

            if scenario is None:
                p2scale = np.ones(len(p2pts))
            else:
                p2scale = np.asarray([
                    pass2_point_corr.get(str(pid), global_corr.get("pass2", 1.0))
                    for pid in p2pts["point_id_b"].astype(str)
                ], dtype=float)
            #endif

            _draw_measurement_series(
                ax,
                p2pts["phi_b"].to_numpy(float),
                p2scale * raw_y,
                np.abs(p2scale) * raw_e,
                label="_nolegend_",
                dataset_key="pass2",
                xoffset=+4.0,
                markersize=4.4,
            )

            for ka in DATASET_ORDER:
                if ka in ("pass2",) or ka in omit:
                    continue
                #endif
                da = d.loc[d["dataset_a"] == ka].sort_values("phi_b")
                if da.empty:
                    continue
                #endif
                scale = global_corr.get(ka, 1.0) if scenario is not None else 1.0
                _draw_measurement_series(
                    ax,
                    da["phi_b"].to_numpy(float),
                    scale * da["xs_a_to_b_km15"].to_numpy(float),
                    abs(scale) * da["point_unc_a_to_b_km15"].to_numpy(float),
                    label="_nolegend_",
                    dataset_key=ka,
                    xoffset=offsets.get(ka, 0.0),
                    markersize=3.7,
                    alpha=0.84,
                )
            #endfor

            xb = float(np.median(p2pts["xB_b"]))
            q2 = float(np.median(p2pts["Q2_b"]))
            tt = float(np.median(p2pts["t_abs_b"]))
            ebeam = float(np.median(p2pts["ebeam_b"]))
            model = get_dense_model_curve(
                emff,
                model_curve_cache,
                dataset_key="pass2",
                ebeam=ebeam,
                xB=xb,
                Q2=q2,
                t_abs=tt,
            )

            for mk in ["bh", "km15"]:
                s = MODEL_STYLES[mk]
                ax.plot(
                    model["phi_deg"], model[mk],
                    color=s["color"], linestyle=s["linestyle"], lw=s["linewidth"],
                )
            #endfor

            yarrays = [
                p2scale * raw_y,
                model["bh"].to_numpy(float),
                model["km15"].to_numpy(float),
            ]
            for ka, da in d.groupby("dataset_a"):
                if ka in omit:
                    continue
                #endif
                sc = global_corr.get(str(ka), 1.0) if scenario is not None else 1.0
                yarrays.append(sc * da["xs_a_to_b_km15"].to_numpy(float))
            #endfor

            ylo, yhi = _robust_positive_log_limits(yarrays)
            ax.set_yscale("log")
            ax.set_ylim(ylo, yhi)
            ax.set_xlim(0, 360)
            ax.set_xticks([0, 90, 180, 270, 360])
            ax.grid(alpha=0.18)
            ax.set_title(
                f"bin {int(round(panel_row['published_bin']))}: "
                + rf"$x_B={xb:.3f}$, $Q^2={q2:.2f}$, $|t|={tt:.3f}$",
                fontsize=8.7,
            )

            if row == nrows - 1:
                ax.set_xlabel(r"$\phi$ (deg)")
            #endif
            if col == 0:
                ax.set_ylabel(
                    r"$d^4\sigma/(dQ^2\,dx_B\,d|t|\,d\phi)$ (nb/GeV$^4$)",
                    fontsize=8.3,
                )
            #endif
        #endfor

        handles, labels = build_complete_pass2_anchor_legend(
            matches,
            dataset_fit_table=dataset_fit_table,
            scenario=scenario,
            omit_datasets=omit_datasets,
        )

        if scenario is None:
            title = "Published world data compared at CLAS12 pass-2 Hayward kinematics"
            subtitle = (
                "External measurements are transported point-by-point to pass-2 Hayward "
                "kinematics with KM15; pass-2 is shown raw. "
                r"Error bars = stat $\oplus$ point-to-point syst."
            )
        else:
            title = f"Systematic-nuisance-adjusted world data at pass-2 Hayward kinematics — {scenario}"
            metric_text = ""
            if metrics:
                metric_text = (
                    rf"  Global $\chi^2/\mathrm{{ndf}}="
                    f"{metrics.get('chi2_per_ndf', np.nan):.2f}$; "
                    rf"$\beta_{{corr}}={metrics.get('beta_pass2_corr_scale',0):+.2f}$."
                )
            #endif
            subtitle = (
                "Experiment-wide normalization nuisances are fitted globally; "
                "full-correlation scenarios also fit one publication-level "
                "pass-2 kinematic correlated-scale nuisance."
                + metric_text
            )
        #endif

        xb_lo, xb_hi = spec["xb_key"]
        title += rf"; $x_B\in[{xb_lo:.3f},{xb_hi:.3f}]$"
        fig.suptitle(title, y=0.994, fontsize=13.5)
        fig.legend(
            handles, labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.958),
            ncol=5,
            frameon=False,
            fontsize=7.1,
        )
        fig.text(0.5, 0.900, subtitle, ha="center", va="top", fontsize=7.7)
        fig.tight_layout(rect=[0.035, 0.035, 0.995, 0.865])

        prefix = "pass2_world_raw" if scenario is None else f"pass2_world_{scenario}"
        fig.savefig(outdir / f"{prefix}_page{ipage+1:02d}.png", dpi=220)
        plt.close(fig)
    #endfor
#enddef


def build_complete_lee_anchor_legend(
        matches: pd.DataFrame,
        *,
        normalization_dataset_table: Optional[pd.DataFrame] = None,
        normalization_scenario: Optional[str] = None,
        omit_datasets: Sequence[str] = ()) -> Tuple[List[Line2D], List[str]]:
    """
    Build a canvas-independent legend containing ALL six published datasets
    plus BH and KM15.

    This intentionally does not inspect which datasets happen to occur on the
    first panel/page.  Therefore the visual key is identical from canvas to
    canvas.

    In normalization-adjusted variants the label gives the fitted global data
    correction.  A dataset deliberately excluded from that scenario is still
    shown in the legend and is marked '(excluded)'.
    """
    omit_set = {str(x) for x in omit_datasets}

    # Recover the normalization fraction from the full matched table rather
    # than from whichever particular canvas is being drawn.
    norm_frac: Dict[str, float] = {}
    for key in DATASET_ORDER:
        vals_a = matches.loc[
            matches["dataset_a"] == key,
            "norm_frac_a",
        ].to_numpy(float)
        vals_b = matches.loc[
            matches["dataset_b"] == key,
            "norm_frac_b",
        ].to_numpy(float)
        vals = np.concatenate([vals_a, vals_b])
        vals = vals[np.isfinite(vals)]
        norm_frac[key] = float(np.nanmedian(vals)) if vals.size else np.nan
    #endfor

    correction = {key: 1.0 for key in DATASET_ORDER}
    if (
        normalization_dataset_table is not None
        and not normalization_dataset_table.empty
    ):
        for r in normalization_dataset_table.itertuples(index=False):
            correction[str(r.dataset)] = float(r.data_correction_scale)
        #endfor
    #endif

    handles: List[Line2D] = []
    labels: List[str] = []

    for key in DATASET_ORDER:
        style = DATASET_STYLES[key]
        handles.append(
            Line2D(
                [0], [0],
                linestyle="none",
                marker=style["marker"],
                markersize=5.5,
                markerfacecolor=style["color"],
                markeredgecolor=style["color"],
                color=style["color"],
            )
        )

        if key in omit_set:
            label = f"{DATASET_LABELS[key]} (excluded)"
        elif normalization_scenario is None:
            n = norm_frac.get(key, np.nan)
            label = (
                f"{DATASET_LABELS[key]} ({100.0*n:.1f}% norm)"
                if np.isfinite(n)
                else DATASET_LABELS[key]
            )
        else:
            shift = 100.0 * (float(correction.get(key, 1.0)) - 1.0)
            constraint = "free" if key in GLOBAL_NORM_FREE_DATASETS else ""
            suffix = f"{shift:+.1f}%"
            if constraint:
                suffix += ", free"
            #endif
            label = f"{DATASET_LABELS[key]} ({suffix})"
        #endif

        labels.append(label)
    #endfor

    for model_key, model_label in [("bh", "BH"), ("km15", "KM15")]:
        style = MODEL_STYLES[model_key]
        handles.append(
            Line2D(
                [0], [0],
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=style["linewidth"],
            )
        )
        labels.append(model_label)
    #endfor

    return handles, labels
#enddef


def plot_lee_anchor_world_panels(
        matches: pd.DataFrame,
        panel_summary: pd.DataFrame,
        outdir: Path,
        emff,
        model_curve_cache: Dict[Tuple, pd.DataFrame],
        *,
        max_pages: int = PANEL_MAX_WORLD_ANCHOR_PAGES,
        normalization_dataset_table: Optional[pd.DataFrame] = None,
        normalization_metrics: Optional[Dict[str, float]] = None,
        normalization_scenario: Optional[str] = None,
        omit_datasets: Sequence[str] = (),
        omit_point_ids: Sequence[str] = (),
        saylor_tmin: Optional[float] = None) -> None:
    """
    Produce multi-dataset world-data overlays centered on CLAS12 Lee 2026 bins.

    With normalization_scenario=None the measured cross sections are displayed
    exactly as published/transported.

    With a normalization scenario, one globally fitted correction factor per
    experiment is applied to both its central values and pointwise error bars.
    Those factors are obtained from the model-independent common-anchor fit
    defined above; KM15 is used only to transport external measurements to Lee.
    """
    if matches.empty or panel_summary.empty:
        return
    #endif

    lee = matches.loc[matches["dataset_b"] == "lee2026"].copy()
    if lee.empty:
        return
    #endif

    omit_dataset_set = {str(x) for x in omit_datasets}
    omit_point_set = {str(x) for x in omit_point_ids}
    if omit_dataset_set:
        lee = lee.loc[~lee["dataset_a"].isin(omit_dataset_set)].copy()
    #endif
    if omit_point_set:
        lee = lee.loc[~lee["point_id_a"].astype(str).isin(omit_point_set)].copy()
    #endif

    if saylor_tmin is not None:
        threshold = float(saylor_tmin)
        keep = (
            (lee["dataset_a"].astype(str) != "saylor2018")
            | (lee["t_abs_a"].to_numpy(float) >= threshold)
        )
        lee = lee.loc[keep].copy()
    #endif

    if "published_bin_b" in lee.columns and np.isfinite(lee["published_bin_b"]).any():
        lee["anchor_group"] = [
            (
                f"bin_{int(round(v))}"
                if np.isfinite(v) else f"point_{pid}"
            )
            for v, pid in zip(
                lee["published_bin_b"].to_numpy(float),
                lee["point_id_b"].astype(str),
            )
        ]
    else:
        return
    #endif

    # Keep only panels that still contain enough information after any
    # scenario-specific omissions.
    available_groups = set(lee["anchor_group"].astype(str))
    selected = panel_summary.loc[
        panel_summary["anchor_group"].astype(str).isin(available_groups)
    ].copy()

    # For the Lee anchors, organization takes precedence over ranking:
    # every qualifying bin is shown (unless max_pages > 0 is explicitly set)
    # and panels are monotonically ordered by the published CLAS12 bin number.
    selected = selected.sort_values(
        ["published_bin", "xB_ref", "Q2_ref", "t_abs_ref"],
        ascending=True,
    ).reset_index(drop=True)

    # max_pages now limits xB-organized canvases, not an arbitrary number of
    # sequential panels; never truncate the middle of a physical xB bin.
    if selected.empty:
        return
    #endif

    # Global correction factors for the normalization-adjusted variants.
    correction = {key: 1.0 for key in DATASET_ORDER}
    if normalization_dataset_table is not None and not normalization_dataset_table.empty:
        for r in normalization_dataset_table.itertuples(index=False):
            correction[str(r.dataset)] = float(r.data_correction_scale)
        #endfor
    #endif

    outdir.mkdir(parents=True, exist_ok=True)
    canvas_specs = _world_canvas_specs_by_xb(selected)
    if int(max_pages) > 0:
        canvas_specs = canvas_specs[:int(max_pages)]
    #endif
    npages = len(canvas_specs)

    # Stable visual offsets are dataset-specific and therefore never change
    # when a canvas contains a different subset of measurements.
    external_offsets = {
        "jo2015": -4.0,
        "defurne2015": -2.5,
        "defurne2017": -1.0,
        "saylor2018": +1.0,
        "georges2022": +2.5,
    }

    for ipage, spec in enumerate(canvas_specs):
        nrows = int(spec["nrows"])
        ncols = int(spec["ncols"])
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(max(8.0, 3.65*ncols), max(5.8, 3.15*nrows + 1.8)),
            squeeze=False,
        )
        placements = spec["placements"]

        for iax, ax in enumerate(axes.ravel()):
            row = iax // ncols
            col = iax % ncols
            panel_row = placements.get((row, col))
            if panel_row is None:
                ax.axis("off")
                continue
            #endif

            group = str(panel_row["anchor_group"])
            d = lee.loc[lee["anchor_group"] == group].copy()
            if d.empty:
                ax.axis("off")
                continue
            #endif

            lee_points = d.sort_values("phi_b").drop_duplicates("point_id_b")
            first_panel = (iax == 0)

            # Lee itself is displayed once per anchor point.
            lee_scale = correction.get("lee2026", 1.0)
            norm_lee = 100.0 * float(lee_points["norm_frac_b"].iloc[0])
            lee_label = f"{DATASET_LABELS['lee2026']}"
            if normalization_scenario is None:
                lee_label += f" ({norm_lee:.1f}% norm)"
            else:
                lee_label += f" ({100.0 * (lee_scale - 1.0):+.1f}%)"
            #endif

            _draw_measurement_series(
                ax,
                lee_points["phi_b"].to_numpy(float),
                lee_scale * lee_points["xs_b"].to_numpy(float),
                abs(lee_scale) * lee_points["point_unc_b"].to_numpy(float),
                label=(lee_label if first_panel else "_nolegend_"),
                dataset_key="lee2026",
                xoffset=+4.0,
                markersize=4.2,
            )

            for ka in DATASET_ORDER:
                if ka == "lee2026" or ka in omit_dataset_set:
                    continue
                #endif

                da = d.loc[d["dataset_a"] == ka].sort_values("phi_b")
                if da.empty:
                    continue
                #endif

                scale = correction.get(ka, 1.0)
                norm_a = 100.0 * float(da["norm_frac_a"].iloc[0])
                label = DATASET_LABELS[ka]
                if normalization_scenario is None:
                    label += f" ({norm_a:.1f}% norm)"
                else:
                    label += f" ({100.0 * (scale - 1.0):+.1f}%)"
                #endif

                _draw_measurement_series(
                    ax,
                    da["phi_b"].to_numpy(float),
                    scale * da["xs_a_to_b_km15"].to_numpy(float),
                    abs(scale) * da["point_unc_a_to_b_km15"].to_numpy(float),
                    label=(label if first_panel else "_nolegend_"),
                    dataset_key=ka,
                    xoffset=external_offsets.get(ka, 0.0),
                    markersize=3.7,
                    alpha=0.84,
                )
            #endfor

            # Dense model scan at the fixed Lee-bin representative kinematics.
            xb = float(np.median(lee_points["xB_b"]))
            q2 = float(np.median(lee_points["Q2_b"]))
            tt = float(np.median(lee_points["t_abs_b"]))
            ebeam = float(np.median(lee_points["ebeam_b"]))
            model = get_dense_model_curve(
                emff,
                model_curve_cache,
                dataset_key="lee2026",
                ebeam=ebeam,
                xB=xb,
                Q2=q2,
                t_abs=tt,
            )

            bh_style = MODEL_STYLES["bh"]
            km_style = MODEL_STYLES["km15"]
            ax.plot(
                model["phi_deg"], model["bh"],
                color=bh_style["color"],
                linestyle=bh_style["linestyle"],
                lw=bh_style["linewidth"],
                label=("BH" if first_panel else "_nolegend_"),
            )
            ax.plot(
                model["phi_deg"], model["km15"],
                color=km_style["color"],
                linestyle=km_style["linestyle"],
                lw=km_style["linewidth"],
                label=("KM15" if first_panel else "_nolegend_"),
            )

            yarrays = [
                lee_scale * lee_points["xs_b"].to_numpy(float),
                model["bh"].to_numpy(float),
                model["km15"].to_numpy(float),
            ]
            for ka, da in d.groupby("dataset_a", sort=False):
                if ka in omit_dataset_set:
                    continue
                #endif
                yarrays.append(
                    correction.get(str(ka), 1.0)
                    * da["xs_a_to_b_km15"].to_numpy(float)
                )
            #endfor

            ylo, yhi = _robust_positive_log_limits(yarrays)
            ax.set_yscale("log")
            ax.set_ylim(ylo, yhi)
            ax.set_xlim(0.0, 360.0)
            ax.set_xticks([0, 90, 180, 270, 360])
            ax.grid(alpha=0.18)

            bin_text = ""
            if np.isfinite(panel_row["published_bin"]):
                bin_text = f"bin {int(round(panel_row['published_bin']))}: "
            #endif
            ax.set_title(
                bin_text + rf"$x_B={xb:.3f}$, $Q^2={q2:.2f}$, $|t|={tt:.3f}$",
                fontsize=8.8,
            )

            if row == nrows - 1:
                ax.set_xlabel(r"$\phi$ (deg)")
            #endif
            if col == 0:
                ax.set_ylabel(
                    r"$d^4\sigma/(dQ^2\,dx_B\,d|t|\,d\phi)$ (nb/GeV$^4$)",
                    fontsize=8.5,
                )
            #endif
        #endfor

        handles, labels = build_complete_lee_anchor_legend(
            matches,
            normalization_dataset_table=normalization_dataset_table,
            normalization_scenario=normalization_scenario,
            omit_datasets=omit_datasets,
        )

        if normalization_scenario is None:
            title = "Published world data mapped to CLAS12 pass-1 kinematics"
            subtitle = (
                "External measurements are transported point-by-point to the exact "
                "Lee 2026 kinematics with KM15; displayed data are not normalization-rescaled. "
                r"Error bars = stat $\oplus$ point-to-point syst."
            )
        else:
            scenario_label = GLOBAL_NORM_SCENARIO_LABELS.get(
                str(normalization_scenario),
                str(normalization_scenario),
            )
            title = (
                "Normalization-adjusted world data at CLAS12 pass-1 kinematics"
                f" — {scenario_label}"
            )
            metric_text = ""
            if normalization_metrics:
                metric_text = (
                    rf"  Global $\chi^2/\mathrm{{ndf}}="
                    f"{normalization_metrics.get('chi2_per_ndf', np.nan):.2f}$ "
                    f"({int(normalization_metrics.get('ndf', 0))} dof)."
                )
            #endif
            subtitle = (
                "One global multiplicative normalization per experiment; "
                "Gaussian priors use quoted correlated normalizations and Georges is free."
                + metric_text
            )
        #endif

        xb_lo, xb_hi = spec["xb_key"]
        title += rf"; $x_B\in[{xb_lo:.3f},{xb_hi:.3f}]$"
        fig.suptitle(title, y=0.994, fontsize=13.5)
        if handles:
            fig.legend(
                handles, labels,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.958),
                ncol=4,
                frameon=False,
                fontsize=7.4,
            )
        #endif
        fig.text(
            0.5, 0.900,
            subtitle,
            ha="center", va="top", fontsize=7.8,
        )
        _synchronize_canvas_y_limits(axes, mode=PANEL_Y_SCALE_MODE)
        fig.tight_layout(rect=[0.035, 0.035, 0.995, 0.865])

        if normalization_scenario is None:
            prefix = "world_data_at_lee_kinematics"
        else:
            prefix = f"world_data_at_lee_kinematics_normfit_{normalization_scenario}"
        #endif

        fig.savefig(
            outdir / f"{prefix}_page{ipage + 1:02d}.png",
            dpi=220,
        )
        plt.close(fig)
    #endfor
#enddef


# =============================================================================
# Output bookkeeping
# =============================================================================



def build_exact_lee_hayward_comparison(world: pd.DataFrame) -> pd.DataFrame:
    """
    Match Lee E214M1 to Hayward by the SAME original analysis bin, with NO
    model transport applied to the cross sections.

    The shared identity is (analysis Bin Name, phi bin).  The pass-1 central
    means may differ slightly from the pass-2 means because the event samples
    and reconstruction differ; those differences are recorded as diagnostics.

    A KM15 mean-kinematics transport factor is also reported, but it is NOT
    applied to the primary direct ratio.  It quantifies how much of the direct
    ratio could plausibly come from shifted bin means alone.
    """
    lee = world.loc[world["dataset"] == "lee2026"].copy()
    hay = world.loc[world["dataset"] == "pass2"].copy()
    if lee.empty or hay.empty:
        return pd.DataFrame()
    #endif
    if "analysis_bin_name" not in lee.columns:
        raise RuntimeError("Lee table lacks analysis_bin_name legacy mapping")
    #endif

    rows = []
    used_h = set()
    for lr in lee.itertuples(index=False):
        if not np.isfinite(float(lr.analysis_bin_name)):
            continue
        #endif
        cand = hay.loc[
            pd.to_numeric(hay["published_bin"], errors="coerce")
            == float(lr.analysis_bin_name)
        ].copy()
        if cand.empty:
            continue
        #endif

        # Same 3D analysis bin; choose the pass-2 phi point whose phi-bin center
        # is nearest the legacy phi row attached to the published Lee point.
        phi_ref = (
            float(lr.legacy_phiavg)
            if hasattr(lr, "legacy_phiavg") and np.isfinite(float(lr.legacy_phiavg))
            else float(lr.phi_deg)
        )
        dphi = np.abs(
            ((cand["phi_deg"].to_numpy(float) - phi_ref + 180.0) % 360.0) - 180.0
        )
        ih = int(np.nanargmin(dphi))
        hr = cand.iloc[ih]
        if float(dphi[ih]) > 12.0:
            continue
        #endif
        hid = str(hr["point_id"])
        if hid in used_h:
            continue
        #endif
        used_h.add(hid)

        ratio_direct = float(hr["xs"] / lr.xs)
        ratio_lee_over_h = float(lr.xs / hr["xs"])
        km_factor = float(hr["km15_native"] / lr.km15_native)
        lee_transported = float(lr.xs * km_factor)
        ratio_transport = float(hr["xs"] / lee_transported)

        sigma_direct = math.sqrt(
            float(hr["point_unc_abs"])**2 + float(lr.point_unc_abs)**2
        )
        pull_direct = (
            (float(hr["xs"]) - float(lr.xs)) / sigma_direct
            if sigma_direct > 0.0 else np.nan
        )

        row = {
            "analysis_bin_name": float(lr.analysis_bin_name),
            "lee_point_id": str(lr.point_id),
            "hayward_point_id": hid,
            "lee_published_bin": float(lr.published_bin),
            "lee_phi": float(lr.phi_deg),
            "hayward_phi": float(hr["phi_deg"]),
            "legacy_phi_center": phi_ref,
            "dphi_means_deg": float(
                ((float(hr["phi_deg"]) - float(lr.phi_deg) + 180.0) % 360.0) - 180.0
            ),
            "lee_xB": float(lr.xB),
            "hayward_xB": float(hr["xB"]),
            "dxB_means": float(hr["xB"] - lr.xB),
            "lee_Q2": float(lr.Q2),
            "hayward_Q2": float(hr["Q2"]),
            "dQ2_means": float(hr["Q2"] - lr.Q2),
            "lee_t_abs": float(lr.t_abs),
            "hayward_t_abs": float(hr["t_abs"]),
            "dt_abs_means": float(hr["t_abs"] - lr.t_abs),
            "lee_xs": float(lr.xs),
            "hayward_xs": float(hr["xs"]),
            "lee_stat": float(lr.stat_abs),
            "hayward_stat": float(hr["stat_abs"]),
            "lee_ptp_sys": float(lr.ptp_sys_abs),
            "hayward_ptp_sys": float(hr["ptp_sys_abs"]),
            "lee_point_unc": float(lr.point_unc_abs),
            "hayward_point_unc": float(hr["point_unc_abs"]),
            "hayward_over_lee_direct": ratio_direct,
            "lee_over_hayward_direct": ratio_lee_over_h,
            "direct_pull_pointwise": pull_direct,
            "km15_lee_to_hayward_mean_transport_factor": km_factor,
            "km15_transport_effect_pct": 100.0 * (km_factor - 1.0),
            "hayward_over_lee_after_km15_mean_transport": ratio_transport,
            "p_theta": float(hr.get("p_theta", np.nan)),
            "g_theta": float(hr.get("g_theta", np.nan)),
            "e_theta": float(hr.get("e_theta", np.nan)),
            "p_p": float(hr.get("p_p", np.nan)),
            "corr_scale_frac": float(hr.get("corr_scale_frac", np.nan)),
        }
        rows.append(row)
    #endfor

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    #endif

    # Detector-region proxies only.  These are deliberately NOT labeled as
    # event topologies because a 4D cross-section bin can contain multiple
    # topology classes.  They are useful for testing coupled p/g-angle trends.
    out["proton_angle_region"] = np.where(
        out["p_theta"] < 35.0,
        "low-p-theta (<35 deg)",
        np.where(out["p_theta"] >= 40.0, "high-p-theta (>=40 deg)", "transition (35-40 deg)"),
    )
    out["photon_angle_region"] = np.where(
        out["g_theta"] <= 5.5,
        "FT-like gamma angle",
        "FD-like gamma angle",
    )
    out["angle_region_2d"] = (
        out["proton_angle_region"].astype(str)
        + " / "
        + out["photon_angle_region"].astype(str)
    )
    return out
#enddef


def summarize_exact_lee_hayward(exact: pd.DataFrame) -> pd.DataFrame:
    if exact.empty:
        return pd.DataFrame()
    #endif
    rows = []
    groups = [("all", exact)]
    for col in ["proton_angle_region", "photon_angle_region", "angle_region_2d"]:
        for label, d in exact.groupby(col, dropna=False):
            groups.append((f"{col}: {label}", d))
        #endfor
    #endfor

    for label, d in groups:
        r = d["hayward_over_lee_direct"].to_numpy(float)
        tf = d["km15_lee_to_hayward_mean_transport_factor"].to_numpy(float)
        rows.append({
            "selection": label,
            "N": int(len(d)),
            "median_hayward_over_lee": float(np.nanmedian(r)),
            "mean_hayward_over_lee": float(np.nanmean(r)),
            "central68_low": float(np.nanpercentile(r, 16.0)),
            "central68_high": float(np.nanpercentile(r, 84.0)),
            "rms_direct_pointwise_pull": float(np.sqrt(np.nanmean(d["direct_pull_pointwise"]**2))),
            "median_abs_km15_mean_transport_pct": float(np.nanmedian(np.abs(100.0*(tf-1.0)))),
            "p95_abs_km15_mean_transport_pct": float(np.nanpercentile(np.abs(100.0*(tf-1.0)), 95.0)),
            "median_dxB": float(np.nanmedian(d["dxB_means"])),
            "median_dQ2_GeV2": float(np.nanmedian(d["dQ2_means"])),
            "median_dt_GeV2": float(np.nanmedian(d["dt_abs_means"])),
            "median_dphi_deg": float(np.nanmedian(d["dphi_means_deg"])),
        })
    #endfor
    return pd.DataFrame(rows)
#enddef


def plot_exact_lee_hayward_diagnostics(exact: pd.DataFrame, outdir: Path) -> None:
    if exact.empty:
        return
    #endif
    outdir.mkdir(parents=True, exist_ok=True)

    variables = [
        ("hayward_t_abs", r"$|t|$ (GeV$^2$)"),
        ("hayward_xB", r"$x_B$"),
        ("hayward_Q2", r"$Q^2$ (GeV$^2$)"),
        ("hayward_phi", r"$\phi$ (deg)"),
        ("p_theta", r"$\theta_p$ (deg)"),
        ("g_theta", r"$\theta_\gamma$ (deg)"),
        ("e_theta", r"$\theta_e$ (deg)"),
        ("p_p", r"$p_p$ (GeV)"),
    ]
    for col, xlabel in variables:
        if col not in exact.columns or not np.isfinite(exact[col]).any():
            continue
        #endif
        fig, ax = plt.subplots(figsize=(8.4, 5.4))
        ax.scatter(
            exact[col], exact["hayward_over_lee_direct"],
            s=15, alpha=0.45,
        )
        ax.axhline(1.0, color="black", lw=1.0)
        if col == "p_theta":
            ax.axvline(35.0, color="black", lw=0.8, ls="--")
            ax.axvline(40.0, color="black", lw=0.8, ls=":")
        #endif
        if col == "g_theta":
            ax.axvline(5.5, color="black", lw=0.8, ls="--")
        #endif
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\sigma_{\rm Hayward}/\sigma_{\rm Lee}$ (same analysis bin)")
        ax.set_ylim(0.0, min(3.0, max(2.0, float(np.nanpercentile(exact["hayward_over_lee_direct"], 99.0))*1.1)))
        ax.grid(alpha=0.18)
        ax.set_title("CLAS12 pass-2 Hayward / Lee pass-1: exact-bin central-value ratio")
        fig.tight_layout()
        fig.savefig(outdir / f"exact_bin_hayward_over_lee_vs_{col}.png", dpi=220)
        plt.close(fig)
    #endfor

    # Mean-kinematic differences versus the direct ratio: this explicitly tests
    # whether changed event-weighted means could be driving the observed shape.
    deltas = [
        ("dxB_means", r"$\Delta\langle x_B\rangle$"),
        ("dQ2_means", r"$\Delta\langle Q^2\rangle$ (GeV$^2$)"),
        ("dt_abs_means", r"$\Delta\langle |t|\rangle$ (GeV$^2$)"),
        ("dphi_means_deg", r"$\Delta\langle\phi\rangle$ (deg)"),
        ("km15_transport_effect_pct", "KM15 mean-kinematics transport effect (%)"),
    ]
    for col, xlabel in deltas:
        fig, ax = plt.subplots(figsize=(8.4, 5.4))
        ax.scatter(exact[col], exact["hayward_over_lee_direct"], s=15, alpha=0.45)
        ax.axhline(1.0, color="black", lw=1.0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\sigma_{\rm Hayward}/\sigma_{\rm Lee}$")
        ax.grid(alpha=0.18)
        ax.set_title("Same-bin pass-2/pass-1 ratio versus change in reported mean kinematics")
        fig.tight_layout()
        fig.savefig(outdir / f"exact_bin_ratio_vs_{col}.png", dpi=220)
        plt.close(fig)
    #endfor

    # 2D proton/photon-angle region summary: useful because p_theta and g_theta
    # are strongly correlated with topology and with one another.
    summary = (
        exact.groupby(["proton_angle_region", "photon_angle_region"])
        .agg(
            N=("hayward_over_lee_direct", "size"),
            median_ratio=("hayward_over_lee_direct", "median"),
            mean_ratio=("hayward_over_lee_direct", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(outdir / "exact_bin_ratio_by_proton_photon_angle_region.csv", index=False)
#enddef


def load_pass2_period_views(pass2_path: Path) -> pd.DataFrame:
    """
    Load each individual pass-2 run-period cross section from the finalized CSV.

    These are diagnostic views, not new independent published datasets.  Their
    tuple statistical uncertainties are retained exactly as stored in the CSV.
    The tuple third component is also retained for inspection but is not assumed
    to be the final publication point-to-point systematic for the individual
    period unless the CSV explicitly supplies one.
    """
    raw = pd.read_csv(pass2_path, low_memory=False)
    rows = []
    for period, ebeam in PASS2_PERIOD_SPECS.items():
        xs_col = f"normed cross sections, ep->epg, exp, {period}, unpol"
        if xs_col not in raw.columns:
            warnings.warn(f"Pass-2 period column not found: {xs_col}")
            continue
        #endif
        for i, rr in raw.iterrows():
            value, stat, tup_sys = _parse_tuple3_loose(rr.get(xs_col, np.nan))
            if not np.isfinite(value) or value <= 0.0:
                continue
            #endif
            def get_mean(name):
                cands = [
                    f"{name}, {period}",
                    f"{name}, {'10.6 GeV' if period != 'Sp19 Inb' else 'Sp19 Inb'}",
                ]
                for c in cands:
                    if c in raw.columns:
                        v = pd.to_numeric(pd.Series([rr.get(c)]), errors="coerce").iloc[0]
                        if np.isfinite(v):
                            return float(v)
                        #endif
                    #endif
                #endfor
                base = {"xBavg":"xBmin", "Q2avg":"Q2min", "t_abs_avg":"t_abs_min", "phiavg":"phimin"}
                if name in base and base[name] in raw.columns:
                    lo = float(rr.get(base[name], np.nan))
                    hi = float(rr.get(base[name].replace("min","max"), np.nan))
                    if np.isfinite(lo) and np.isfinite(hi):
                        return 0.5*(lo+hi)
                    #endif
                #endif
                return np.nan
            #enddef

            ptheta = np.nan
            gtheta = np.nan
            etheta = np.nan
            for outname, short in [("p_theta","p_theta"),("g_theta","g_theta"),("e_theta","e_theta")]:
                for c in [f"{short}, {period}", f"{short}, {'10.6 GeV' if period != 'Sp19 Inb' else 'Sp19 Inb'}"]:
                    if c in raw.columns:
                        v = pd.to_numeric(pd.Series([rr.get(c)]), errors="coerce").iloc[0]
                        if np.isfinite(v):
                            if outname == "p_theta": ptheta = float(v)
                            if outname == "g_theta": gtheta = float(v)
                            if outname == "e_theta": etheta = float(v)
                            break
                        #endif
                    #endif
                #endfor
            #endfor

            rows.append({
                "period": period,
                "ebeam": float(ebeam),
                "analysis_bin_name": float(rr.get("Bin Name", np.nan)),
                "four_d_bin_index": float(rr.get("bin index", np.nan)),
                "xB": get_mean("xBavg"),
                "Q2": get_mean("Q2avg"),
                "t_abs": get_mean("t_abs_avg"),
                "phi_deg": np.mod(get_mean("phiavg"), 360.0),
                "xs": float(value),
                "stat_abs": float(stat),
                "tuple_sys_abs": float(tup_sys),
                "p_theta": ptheta,
                "g_theta": gtheta,
                "e_theta": etheta,
            })
        #endfor
    #endfor
    return pd.DataFrame(rows)
#enddef


def compare_periods_to_lee_exact(
        periods: pd.DataFrame,
        world: pd.DataFrame,
        emff) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compare all five run periods to Lee.

    For the four 10.6-GeV periods the primary ratio is direct same-bin with NO
    energy/kinematic transport.  For Sp19, Lee is transported from 10.604 GeV
    to 10.2 GeV using KM15 evaluated at the period's mean kinematics.
    """
    lee = world.loc[world["dataset"] == "lee2026"].copy()
    if periods.empty or lee.empty:
        return pd.DataFrame(), pd.DataFrame()
    #endif
    rows = []
    for pr in periods.itertuples(index=False):
        cand = lee.loc[
            pd.to_numeric(lee["analysis_bin_name"], errors="coerce")
            == float(pr.analysis_bin_name)
        ].copy()
        if cand.empty:
            continue
        #endif
        dphi = np.abs(((cand["phi_deg"].to_numpy(float) - float(pr.phi_deg) + 180.0) % 360.0) - 180.0)
        il = int(np.nanargmin(dphi))
        if float(dphi[il]) > 12.0:
            continue
        #endif
        lr = cand.iloc[il]

        lee_ref = float(lr["xs"])

        # The four 10.6-GeV periods use direct same-bin comparison with NO
        # transport.  Sp19 alone requires transport to 10.2 GeV.  Use the
        # already-cached Lee KM15 prediction as the denominator and evaluate
        # only the 10.2-GeV target point.
        transport = 1.0
        if str(pr.period) == "Sp19 Inb":
            from types import SimpleNamespace
            period_model_row = SimpleNamespace(
                dataset="pass2",
                xB=float(pr.xB),
                Q2=float(pr.Q2),
                t_abs=float(pr.t_abs),
                phi_deg=float(pr.phi_deg),
            )
            km_lee = float(lr["km15_native"])
            km_period = evaluate_one_km15(
                emff, period_model_row, float(pr.ebeam)
            )["km15_ep"]
            transport = (
                float(km_period / km_lee)
                if np.isfinite(km_lee) and km_lee != 0.0
                else np.nan
            )
        #endif
        lee_transported = lee_ref * transport

        rows.append({
            "period": str(pr.period),
            "analysis_bin_name": float(pr.analysis_bin_name),
            "four_d_bin_index": float(pr.four_d_bin_index),
            "period_xB": float(pr.xB), "period_Q2": float(pr.Q2),
            "period_t_abs": float(pr.t_abs), "period_phi": float(pr.phi_deg),
            "period_xs": float(pr.xs), "period_stat": float(pr.stat_abs),
            "period_tuple_sys": float(pr.tuple_sys_abs),
            "lee_xB": float(lr["xB"]), "lee_Q2": float(lr["Q2"]),
            "lee_t_abs": float(lr["t_abs"]), "lee_phi": float(lr["phi_deg"]),
            "lee_xs": lee_ref, "lee_stat": float(lr["stat_abs"]),
            "lee_ptp_sys": float(lr["ptp_sys_abs"]),
            "lee_point_unc": float(lr["point_unc_abs"]),
            "direct_period_over_lee": float(pr.xs / lee_ref),
            "km15_lee_to_period_transport_factor": float(transport),
            "km15_transport_effect_pct": 100.0 * (float(transport) - 1.0),
            "period_over_lee_after_km15_transport": (
                float(pr.xs / lee_transported)
                if np.isfinite(lee_transported) and lee_transported > 0.0
                else np.nan
            ),
            "p_theta": float(pr.p_theta), "g_theta": float(pr.g_theta), "e_theta": float(pr.e_theta),
            "needs_energy_transport": str(pr.period) == "Sp19 Inb",
        })
    #endfor

    points = pd.DataFrame(rows)
    summary_rows = []
    if not points.empty:
        for period, d in points.groupby("period"):
            summary_rows.append({
                "period": period,
                "N": int(len(d)),
                "median_direct_period_over_lee": float(np.nanmedian(d["direct_period_over_lee"])),
                "mean_direct_period_over_lee": float(np.nanmean(d["direct_period_over_lee"])),
                "median_period_over_lee_after_km15_transport": float(np.nanmedian(d["period_over_lee_after_km15_transport"])),
                "median_abs_km15_transport_effect_pct": float(np.nanmedian(np.abs(d["km15_transport_effect_pct"]))),
                "p95_abs_km15_transport_effect_pct": float(np.nanpercentile(np.abs(d["km15_transport_effect_pct"]),95)),
                "central68_low": float(np.nanpercentile(d["direct_period_over_lee"], 16)),
                "central68_high": float(np.nanpercentile(d["direct_period_over_lee"], 84)),
                "median_stat_frac_period_pct": float(100*np.nanmedian(d["period_stat"]/d["period_xs"])),
                "energy_transport_required": bool(d["needs_energy_transport"].any()),
            })
        #endfor
    #endif
    return points, pd.DataFrame(summary_rows)
#enddef


def plot_period_lee_diagnostics(points: pd.DataFrame, outdir: Path) -> None:
    if points.empty:
        return
    #endif
    outdir.mkdir(parents=True, exist_ok=True)
    period_styles = {
        "Sp18 Inb": dict(marker="D"),
        "Sp18 Out": dict(marker="P"),
        "Fa18 Inb": dict(marker="^"),
        "Fa18 Out": dict(marker="v"),
        "Sp19 Inb": dict(marker="s"),
    }
    variables = [
        ("period_phi", r"$\phi$ (deg)"),
        ("period_t_abs", r"$|t|$ (GeV$^2$)"),
        ("period_xB", r"$x_B$"),
        ("period_Q2", r"$Q^2$ (GeV$^2$)"),
        ("p_theta", r"$\theta_p$ (deg)"),
        ("g_theta", r"$\theta_\gamma$ (deg)"),
    ]
    for col, xlabel in variables:
        if col not in points.columns or not np.isfinite(points[col]).any():
            continue
        #endif
        fig, ax = plt.subplots(figsize=(9.0, 5.8))
        for period, d in points.groupby("period", sort=False):
            ycol = (
                "period_over_lee_after_km15_transport"
                if period == "Sp19 Inb"
                else "direct_period_over_lee"
            )
            ax.scatter(
                d[col], d[ycol], s=17, alpha=0.45,
                marker=period_styles.get(period, {}).get("marker", "o"),
                label=period,
            )
        #endfor
        ax.axhline(1.0, color="black", lw=1.0)
        if col == "p_theta":
            ax.axvline(35.0, color="black", ls="--", lw=0.8)
            ax.axvline(40.0, color="black", ls=":", lw=0.8)
        #endif
        if col == "g_theta":
            ax.axvline(5.5, color="black", ls="--", lw=0.8)
        #endif
        ax.set_xlabel(xlabel)
        ax.set_ylabel("pass-2 period / Lee pass-1")
        ax.set_ylim(0.0, 2.5)
        ax.grid(alpha=0.18)
        ax.legend(frameon=False, ncol=3)
        ax.set_title(
            "Individual pass-2 run periods vs Lee: exact bins at 10.6 GeV; "
            "Sp19 uses KM15 energy/mean transport"
        )
        fig.tight_layout()
        fig.savefig(outdir / f"run_period_over_lee_vs_{col}.png", dpi=220)
        plt.close(fig)
    #endfor

    # A direct run-period summary plot makes coherent period offsets obvious.
    order = ["Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"]
    rows = []
    for period in order:
        d = points.loc[points["period"] == period]
        if d.empty:
            continue
        #endif
        ycol = (
            "period_over_lee_after_km15_transport"
            if period == "Sp19 Inb"
            else "direct_period_over_lee"
        )
        vals = d[ycol].to_numpy(float)
        rows.append((period, np.nanmedian(vals), np.nanpercentile(vals,16), np.nanpercentile(vals,84), len(d)))
    #endfor
    if rows:
        fig, ax = plt.subplots(figsize=(8.8, 5.4))
        x = np.arange(len(rows))
        med = np.array([r[1] for r in rows])
        lo = med - np.array([r[2] for r in rows])
        hi = np.array([r[3] for r in rows]) - med
        ax.errorbar(x, med, yerr=np.vstack([lo,hi]), fmt="o", capsize=4)
        ax.axhline(1.0, color="black", lw=1.0)
        ax.set_xticks(x, [f"{r[0]}\nN={r[4]}" for r in rows])
        ax.set_ylabel("median pass-2 period / Lee")
        ax.grid(axis="y", alpha=0.18)
        ax.set_title("Run-period dependence of the pass-2/pass-1 cross-section ratio")
        fig.tight_layout()
        fig.savefig(outdir / "run_period_over_lee_summary.png", dpi=220)
        plt.close(fig)
    #endif
#enddef


def halla_covered_region_summary(
        world: pd.DataFrame,
        matches: pd.DataFrame) -> pd.DataFrame:
    """
    Quantify CLAS behavior specifically in kinematics where Hall A has matched
    coverage.  A CLAS point is tagged Hall-A-covered if it participates in a
    selected pairwise match with either Defurne 2015 or Defurne 2017.

    This avoids inventing a rectangular Hall-A acceptance box and uses the same
    matching definition as the rest of the world-data study.
    """
    if matches.empty:
        return pd.DataFrame()
    #endif
    halla = {"defurne2015", "defurne2017"}
    clas = {"jo2015", "saylor2018", "lee2026", "pass2"}
    rows = []
    for ckey in sorted(clas):
        point_ids = set()
        partner_rows = []
        for hkey in halla:
            d1 = matches.loc[(matches["dataset_a"] == ckey) & (matches["dataset_b"] == hkey)]
            if not d1.empty:
                point_ids.update(d1["point_id_a"].astype(str))
                partner_rows.append(d1)
            #endif
            d2 = matches.loc[(matches["dataset_a"] == hkey) & (matches["dataset_b"] == ckey)]
            if not d2.empty:
                point_ids.update(d2["point_id_b"].astype(str))
                partner_rows.append(d2)
            #endif
        #endfor
        d = world.loc[(world["dataset"] == ckey) & world["point_id"].astype(str).isin(point_ids)].copy()
        if d.empty:
            continue
        #endif
        for model, col in [("KM15", "km15_native"), ("BH", "bh_native")]:
            if ckey == "pass2":
                # Hayward keeps its publication-level uncertainty model even
                # inside the Hall-A-covered phase-space restriction: one
                # overall normalization nuisance plus one finalized
                # bin-dependent correlated-scale nuisance.
                fit = fit_pass2_model_publication_nuisances(
                    d["xs"].to_numpy(float),
                    d[col].to_numpy(float),
                    d["point_unc_abs"].to_numpy(float),
                    d["norm_frac"].to_numpy(float),
                    d["corr_scale_frac"].to_numpy(float),
                    include_norm=True,
                    include_corr=True,
                )
                rows.append({
                    "dataset": ckey,
                    "dataset_label": DATASET_LABELS[ckey],
                    "model": model,
                    "N_halla_covered_points": int(fit["N"]),
                    "systematic_treatment": "overall norm + correlated scale",
                    "chi2_per_point": float(fit["chi2_per_point"]),
                    "normalization_beta": float(fit["beta_norm"]),
                    "normalization_scale": float(
                        1.0 + np.nanmedian(d["norm_frac"].to_numpy(float))
                        * float(fit["beta_norm"])
                    ),
                    "correlated_scale_beta": float(fit["beta_corr"]),
                    "combined_nuisance_excursion_sigma": float(fit["combined_nuisance_excursion_sigma"]),
                    "corr_shift_median_pct": float(fit["corr_shift_median_pct"]),
                    "corr_shift_min_pct": float(fit["corr_shift_min_pct"]),
                    "corr_shift_max_pct": float(fit["corr_shift_max_pct"]),
                    "median_abs_fractional_residual": float(fit["median_abs_fractional_residual_fitted"]),
                    "median_abs_fractional_residual_raw": float(fit["median_abs_fractional_residual_raw"]),
                })
            else:
                fit = fit_one_normalization_nuisance(
                    d["xs"].to_numpy(float),
                    d[col].to_numpy(float),
                    d["point_unc_abs"].to_numpy(float),
                    float(np.nanmedian(d["norm_frac"])),
                )
                rows.append({
                    "dataset": ckey,
                    "dataset_label": DATASET_LABELS[ckey],
                    "model": model,
                    "N_halla_covered_points": int(fit["N"]),
                    "systematic_treatment": "overall normalization",
                    "chi2_per_point": float(fit["chi2_per_point"]),
                    "normalization_beta": float(fit["beta"]),
                    "normalization_scale": float(fit["scale"]),
                    "correlated_scale_beta": np.nan,
                    "combined_nuisance_excursion_sigma": abs(float(fit["beta"])),
                    "corr_shift_median_pct": np.nan,
                    "corr_shift_min_pct": np.nan,
                    "corr_shift_max_pct": np.nan,
                    "median_abs_fractional_residual": float(fit["median_abs_fractional_residual"]),
                    "median_abs_fractional_residual_raw": float(np.nanmedian(np.abs(d["xs"].to_numpy(float)-d[col].to_numpy(float))/np.abs(d[col].to_numpy(float)))),
                })
            #endif
        #endfor
    #endfor
    return pd.DataFrame(rows)
#enddef


def clas_vs_halla_pairwise_summary(pair_summary: pd.DataFrame) -> pd.DataFrame:
    """Extract the direct CLAS-versus-Hall-A overlap metrics in one compact table."""
    if pair_summary.empty:
        return pd.DataFrame()
    #endif
    clas = {"jo2015", "saylor2018", "lee2026", "pass2"}
    halla = {"defurne2015", "defurne2017"}
    rows = []
    for r in pair_summary.itertuples(index=False):
        a = str(r.dataset_a)
        b = str(r.dataset_b)
        if not ((a in clas and b in halla) or (a in halla and b in clas)):
            continue
        #endif
        clas_key = a if a in clas else b
        hall_key = b if b in halla else a
        rows.append({
            "clas_dataset": clas_key,
            "clas_dataset_label": DATASET_LABELS[clas_key],
            "halla_dataset": hall_key,
            "halla_dataset_label": DATASET_LABELS[hall_key],
            "N_matches": int(r.N_matches),
            "chi2_per_match_with_relative_norm": float(r.chi2_per_match),
            "raw_pull_rms": float(r.raw_pull_rms),
            "profiled_pull_rms": float(r.profiled_pull_rms),
            "relative_norm_beta": float(r.relative_norm_beta),
            "relative_scale_a_to_b": float(r.relative_scale_a_to_b),
        })
    #endfor
    return pd.DataFrame(rows)
#enddef



def world_has_pass2(world: pd.DataFrame) -> bool:
    """Return True only when a valid CLAS12 pass-2 Hayward sample is loaded."""
    return (
        not world.empty
        and "dataset" in world.columns
        and bool(np.any(world["dataset"].astype(str).to_numpy() == "pass2"))
    )
#enddef


def save_outputs(
        world: pd.DataFrame,
        dataset_summary: pd.DataFrame,
        model_scores: pd.DataFrame,
        matches: pd.DataFrame,
        pair_summary: pd.DataFrame,
        outdir: Path,
        have_gk16: bool,
        emff,
        args) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not have_pass2:
        print(
            "[PASS2 HAYWARD] no finalized pass-2 sample loaded; "
            "all Hayward-specific fits/tables/plots will be skipped",
            flush=True,
        )
    #endif

    tables = outdir / "tables"
    figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    world.to_csv(tables / "canonical_world_data_with_models.csv", index=False)
    dataset_summary.to_csv(tables / "dataset_summary.csv", index=False)
    model_scores.to_csv(tables / "native_model_scores.csv", index=False)

    pass2_model_scores = pd.DataFrame()
    pass2_model_points = pd.DataFrame()
    if have_pass2:
        print(
            "[PASS2 MODEL] evaluating Hayward vs KM15/BH with publication-level "
            "normalization + correlated-scale nuisances",
            flush=True,
        )
        pass2_model_scores, pass2_model_points = (
            make_pass2_model_publication_scores(world)
        )
    #endif
    pass2_model_scores.to_csv(
        tables / "pass2_model_publication_scores.csv",
        index=False,
    )
    pass2_model_points.to_csv(
        tables / "pass2_model_publication_point_residuals.csv",
        index=False,
    )

    pass2_pair_scores = (
        make_pass2_pair_publication_scores(matches)
        if have_pass2 else pd.DataFrame()
    )
    pass2_pair_scores.to_csv(
        tables / "pass2_pairwise_publication_scores.csv",
        index=False,
    )

    if not pass2_model_scores.empty:
        print(
            "\n[PASS2 MODEL] publication-level direct model comparison:",
            flush=True,
        )
        cols = [
            "model", "scenario", "N",
            "chi2_per_point",
            "beta_norm", "beta_corr",
            "normalization_shift_pct",
            "corr_shift_median_pct",
            "corr_shift_min_pct",
            "corr_shift_max_pct",
            "total_shift_median_pct",
            "total_shift_min_pct",
            "total_shift_max_pct",
            "combined_nuisance_excursion_sigma",
            "response_correlation_norm_corr",
        ]
        print(
            pass2_model_scores[cols].to_string(
                index=False,
                float_format=lambda x: f"{x:.4f}",
            ),
            flush=True,
        )
    #endif

    if not pass2_pair_scores.empty:
        print(
            "\n[PASS2 PAIRS] publication-level matched comparisons "
            "(relative norm + Hayward correlated scale):",
            flush=True,
        )
        cols = [
            "dataset_a_label", "N",
            "chi2_per_match",
            "relative_norm_beta",
            "relative_scale_a_to_b",
            "beta_pass2_corr_scale",
            "combined_nuisance_excursion_sigma",
            "corr_shift_median_pct",
            "fitted_pull_rms",
        ]
        print(
            pass2_pair_scores[cols].to_string(
                index=False,
                float_format=lambda x: f"{x:.4f}",
            ),
            flush=True,
        )
    #endif
    pair_summary.to_csv(tables / "pairwise_overlap_summary.csv", index=False)
    if not matches.empty:
        matches.to_csv(tables / "pairwise_matched_points.csv", index=False)
    #endif

    # ------------------------------------------------------------------
    # High-priority CLAS12 pass-1/pass-2 diagnostics.
    # ------------------------------------------------------------------
    exact_lh = build_exact_lee_hayward_comparison(world)
    exact_lh.to_csv(tables / "lee_hayward_exact_same_bin_points.csv", index=False)
    exact_summary = summarize_exact_lee_hayward(exact_lh)
    exact_summary.to_csv(tables / "lee_hayward_exact_same_bin_summary.csv", index=False)
    plot_exact_lee_hayward_diagnostics(
        exact_lh, figures / "lee_hayward_exact_bin_diagnostics"
    )
    if not exact_summary.empty:
        print("\n[LEE/HAYWARD EXACT BIN] no-transport comparison:", flush=True)
        print(
            exact_summary.head(12).to_string(
                index=False, float_format=lambda x: f"{x:.4f}"
            ),
            flush=True,
        )
    #endif

    period_views = load_pass2_period_views(Path(args.resolved_pass2_file))
    period_points, period_summary = compare_periods_to_lee_exact(
        period_views, world, emff
    )
    period_points.to_csv(tables / "pass2_run_periods_vs_lee_points.csv", index=False)
    period_summary.to_csv(tables / "pass2_run_periods_vs_lee_summary.csv", index=False)
    plot_period_lee_diagnostics(
        period_points, figures / "pass2_run_periods_vs_lee"
    )
    if not period_summary.empty:
        print("\n[RUN PERIODS vs LEE]", flush=True)
        print(
            period_summary.to_string(
                index=False, float_format=lambda x: f"{x:.4f}"
            ),
            flush=True,
        )
    #endif

    hall_region = halla_covered_region_summary(world, matches)
    hall_region.to_csv(tables / "clas_in_halla_covered_region_model_scores.csv", index=False)
    clas_halla = clas_vs_halla_pairwise_summary(pair_summary)
    clas_halla.to_csv(tables / "clas_vs_halla_pairwise_summary.csv", index=False)
    if not clas_halla.empty:
        print("\n[CLAS vs HALL A DIRECT OVERLAP]", flush=True)
        print(
            clas_halla.to_string(index=False, float_format=lambda x: f"{x:.4f}"),
            flush=True,
        )
    #endif
    if not hall_region.empty:
        print("\n[HALL A COVERED REGION] CLAS model diagnostics restricted to matched Hall A phase space:", flush=True)
        print(
            hall_region.to_string(index=False, float_format=lambda x: f"{x:.4f}"),
            flush=True,
        )
    #endif

    pair_panel_summary = make_pairwise_panel_summary(matches)
    lee_anchor_summary = make_lee_anchor_panel_summary(matches)
    pass2_anchor_summary = (
        make_pass2_anchor_panel_summary(matches)
        if have_pass2 else pd.DataFrame()
    )

    n_pass2_matches = int(np.sum(matches["dataset_b"].astype(str) == "pass2")) if not matches.empty else 0
    print(
        f"[PASS2 HAYWARD PANELS] external->Hayward matched points={n_pass2_matches:,}; "
        f"qualifying 3D phi-distribution panels={len(pass2_anchor_summary):,}",
        flush=True,
    )
    if n_pass2_matches > 0 and pass2_anchor_summary.empty:
        raise RuntimeError(
            "Hayward matched points exist but no 3D panel groups were formed. "
            "Check the pass-2 Bin Name / 4D bin-index mapping."
        )
    #endif

    pair_panel_summary.to_csv(tables / "pairwise_cross_section_panel_summary.csv", index=False)
    lee_anchor_summary.to_csv(tables / "lee_anchor_panel_summary.csv", index=False)
    pass2_anchor_summary.to_csv(tables / "pass2_anchor_panel_summary.csv", index=False)

    transport_cols = [
        "point_id", "dataset", "dataset_label",
        "ebeam", "xB", "Q2", "t_abs", "phi_deg",
        "xs", "stat_abs", "ptp_sys_abs", "norm_frac",
        "km15_transport_factor", "xs_10p6_km15",
        "stat_10p6_km15", "ptp_sys_10p6_km15", "point_unc_10p6_km15",
    ]
    if have_gk16:
        transport_cols += [
            "gk16_transport_factor", "xs_10p6_gk16",
            "transport_model_unc_abs", "transport_model_unc_frac",
        ]
    #endif
    world[transport_cols].to_csv(tables / "world_data_transported_to_10p6.csv", index=False)

    # ------------------------------------------------------------------
    # Global Lee-anchor normalization fits.
    # ------------------------------------------------------------------
    norm_dataset, norm_points, norm_metrics = (
        run_global_lee_anchor_normalization_scenarios(matches)
    )
    norm_dataset.to_csv(
        tables / "lee_anchor_global_normalization_offsets.csv",
        index=False,
    )

    # Pass-2-specific nuisance treatment.
    p2_ds = pd.DataFrame()
    p2_nui = pd.DataFrame()
    p2_pts = pd.DataFrame()
    p2_metrics = pd.DataFrame()
    if have_pass2:
        p2_ds, p2_nui, p2_pts, p2_metrics = run_pass2_anchor_nuisance_scenarios(matches)
    #endif
    p2_ds.to_csv(tables / "pass2_global_normalization_offsets.csv", index=False)
    p2_nui.to_csv(tables / "pass2_correlated_scale_nuisances.csv", index=False)
    p2_pts.to_csv(tables / "pass2_global_nuisance_point_residuals.csv", index=False)
    p2_metrics.to_csv(tables / "pass2_global_nuisance_summary.csv", index=False)
    norm_points.to_csv(
        tables / "lee_anchor_global_normalization_point_residuals.csv",
        index=False,
    )
    norm_metrics.to_csv(
        tables / "lee_anchor_global_normalization_summary.csv",
        index=False,
    )

    # One in-memory cache is shared by every presentation plot so dense model
    # curves at repeated Lee anchors are evaluated only once during this run.
    model_curve_cache: Dict[Tuple, pd.DataFrame] = {}

    print("[PLOTS] collecting fixed-kinematics model curves needed by all canvases", flush=True)
    model_curve_specs = collect_presentation_model_curve_specs(
        matches,
        pair_panel_summary,
        lee_anchor_summary,
    )
    precompute_presentation_model_curves(
        emff,
        model_curve_cache,
        model_curve_specs,
    )
    print("[PLOTS] model prediction precomputation finished; drawing figures", flush=True)

    plot_native_model_ratios(world, figures / "native_model_ratios", have_gk16)
    plot_transport_uncertainty(world, figures / "transport", have_gk16)
    plot_pairwise_matrix(pair_summary, figures / "pairwise")
    plot_pairwise_pulls(matches, figures / "pairwise")

    plot_pairwise_cross_section_panels(
        matches,
        pair_panel_summary,
        figures / "cross_section_overlays" / "pairwise",
        emff,
        model_curve_cache,
    )

    # Raw, un-rescaled NOMINAL Lee-anchor figures.
    plot_lee_anchor_world_panels(
        matches,
        lee_anchor_summary,
        figures / "cross_section_overlays" / "lee_anchor" / "raw",
        emff,
        model_curve_cache,
    )

    # Requested normalization-adjusted variants.
    scenario_specs = [
        (
            "nominal",
            (),
            (),
            None,
            figures / "cross_section_overlays" / "lee_anchor" / "normfit_nominal",
        ),
        (
            "saylor_tmin_0p343",
            (),
            (),
            SAYLOR_TMIN_DIAGNOSTIC_GEV2,
            figures / "cross_section_overlays" / "lee_anchor" / "normfit_saylor_tmin_0p343",
        ),
        (
            "without_saylor",
            ("saylor2018",),
            (),
            None,
            figures / "cross_section_overlays" / "lee_anchor" / "normfit_without_saylor",
        ),
    ]

    for scenario, omit_datasets, omit_points, saylor_tmin, scenario_outdir in scenario_specs:
        dtab = norm_dataset.loc[norm_dataset["scenario"] == scenario].copy()
        mrow = norm_metrics.loc[norm_metrics["scenario"] == scenario]
        metrics = (
            mrow.iloc[0].to_dict()
            if not mrow.empty else {}
        )
        plot_lee_anchor_world_panels(
            matches,
            lee_anchor_summary,
            scenario_outdir,
            emff,
            model_curve_cache,
            normalization_dataset_table=dtab,
            normalization_metrics=metrics,
            normalization_scenario=scenario,
            omit_datasets=omit_datasets,
            omit_point_ids=omit_points,
            saylor_tmin=saylor_tmin,
        )
    #endfor

    if have_pass2:
        print(
            "[PLOTS PASS2 MODEL] residual diagnostics for KM15/BH",
            flush=True,
        )
        plot_pass2_model_residual_diagnostics(
            pass2_model_points,
            figures / "pass2_model_diagnostics",
        )
    #endif

    # ------------------------------------------------------------------
    # Dedicated CLAS12 pass-2 Hayward anchor canvases.
    #
    # Earlier versions built the pass-2 panel summary and precomputed the
    # corresponding BH/KM15 curves, but omitted these final drawing calls.
    # ------------------------------------------------------------------
    if have_pass2:
        pass2_root = (
            figures
            / "cross_section_overlays"
            / "pass2_anchor"
        )

        print(
            f"[PLOTS PASS2] generating dedicated CLAS12 pass-2 Hayward canvases "
            f"from {len(pass2_anchor_summary):,} qualifying 3D bins",
            flush=True,
        )

        # Raw comparison.
        plot_pass2_anchor_world_panels(
            matches,
            pass2_anchor_summary,
            pass2_root / "raw",
            emff,
            model_curve_cache,
        )

        # Nuisance-adjusted comparison variants.
        pass2_plot_specs = [
            (
                "norm_only_nominal",
                (),
                None,
                "norm_only_nominal",
            ),
            (
                "full_corr_nominal",
                (),
                None,
                "full_corr_nominal",
            ),
            (
                "full_corr_saylor_tmin_0p343",
                (),
                SAYLOR_TMIN_DIAGNOSTIC_GEV2,
                "full_corr_saylor_tmin_0p343",
            ),
            (
                "full_corr_without_saylor",
                ("saylor2018",),
                None,
                "full_corr_without_saylor",
            ),
        ]

        for iscenario, (
                scenario,
                omitted_datasets,
                saylor_tmin,
                dirname) in enumerate(pass2_plot_specs, start=1):

            print(
                f"[PLOTS PASS2] scenario {iscenario}/{len(pass2_plot_specs)}: "
                f"{scenario}",
                flush=True,
            )

            dataset_fit = p2_ds.loc[
                p2_ds["scenario"] == scenario
            ].copy()

            point_fit = p2_pts.loc[
                p2_pts["scenario"] == scenario
            ].copy()

            metric_rows = p2_metrics.loc[
                p2_metrics["scenario"] == scenario
            ]
            metrics = (
                metric_rows.iloc[0].to_dict()
                if not metric_rows.empty
                else {}
            )

            plot_pass2_anchor_world_panels(
                matches,
                pass2_anchor_summary,
                pass2_root / dirname,
                emff,
                model_curve_cache,
                dataset_fit_table=dataset_fit,
                point_fit_table=point_fit,
                metrics=metrics,
                scenario=scenario,
                omit_datasets=omitted_datasets,
                saylor_tmin=saylor_tmin,
            )
        #endfor

        print(
            f"[PLOTS PASS2] complete -> {pass2_root}",
            flush=True,
        )

    else:
        print(
            "[PLOTS PASS2] skipped: no finalized pass-2 sample loaded",
            flush=True,
        )
    #endif

    return norm_dataset, norm_metrics
#enddef

def print_summary(dataset_summary: pd.DataFrame, model_scores: pd.DataFrame, pair_summary: pd.DataFrame, have_gk16: bool, norm_dataset: Optional[pd.DataFrame] = None, norm_metrics: Optional[pd.DataFrame] = None) -> None:
    print("\n" + "=" * 80)
    print("DATASET SUMMARY")
    print("=" * 80)
    cols = ["dataset_label", "N", "ebeam_min_GeV", "ebeam_max_GeV", "median_point_unc_pct", "correlated_norm_pct"]
    print(dataset_summary[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\n" + "=" * 80)
    print("NATIVE-KINEMATICS MODEL SCORES")
    print("=" * 80)
    print(model_scores[[
        "dataset_label", "model", "N", "beta", "scale", "chi2_per_point", "median_abs_fractional_residual"
    ]].to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\n" + "=" * 80)
    print("PAIRWISE MATCHED-DATA SUMMARY")
    print(
        "[NOTE] chi2_per_match here includes the pair's relative overall "
        "normalization nuisance, but NOT the pass-2 kinematic correlated-scale "
        "nuisance.  See pass2_pairwise_publication_scores.csv for the "
        "publication-style Hayward comparisons."
    )
    print("=" * 80)
    if pair_summary.empty:
        print("No pairwise matches were found with the configured windows.")
    else:
        print(pair_summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    #endif

    if norm_metrics is not None and not norm_metrics.empty:
        print("\n" + "=" * 80)
        print("GLOBAL LEE-ANCHOR NORMALIZATION CONSISTENCY")
        print("=" * 80)
        metric_cols = [
            "scenario", "N_data", "N_anchors", "N_priors",
            "chi2_data", "chi2_prior", "ndf", "chi2_per_ndf",
        ]
        print(
            norm_metrics[metric_cols].to_string(
                index=False,
                float_format=lambda x: f"{x:.4f}",
            )
        )

        if norm_dataset is not None and not norm_dataset.empty:
            print("\nNormalization corrections applied to data in fitted overlays:")
            cols = [
                "scenario", "dataset_label", "N_points",
                "quoted_norm_pct", "normalization_constraint",
                "data_correction_pct", "beta_prior_sigma",
            ]
            print(
                norm_dataset[cols].to_string(
                    index=False,
                    float_format=lambda x: f"{x:.4f}",
                )
            )
        #endif
    #endif

    if not have_gk16:
        print("\n[NOTE] World-data transport uses KM15 only by design.")
        print("       No PARTONS/GK16 files are generated by this script.")
    #endif
#enddef


# =============================================================================
# CLI
# =============================================================================


def build_parser() -> argparse.ArgumentParser:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(
        description="World proton DVCS/BH cross-section consistency and model comparison"
    )
    p.add_argument("--outdir", default="../output/world_dvcs_cross_section_comparison")
    p.add_argument("--emff-script", default=str(here / "extract_emff_from_dvcs_bh.py"))
    p.add_argument("--saylor-file", default=str(here / "import" / "saylor_CLAS6.txt"))
    p.add_argument("--georges-file", default=str(here / "import" / "E12-06-114.xlsx"))
    p.add_argument("--lee-file", default=str(here / "import" / "clasdb_E214M1.txt"))
    p.add_argument(
        "--pass1-legacy-file",
        default=str(here / "import" / "all_bin_v3.csv"),
        help=(
            "Legacy/preliminary pass-1 all_bin_v3.csv used ONLY for original "
            "bin identities, boundaries, and kinematic bookkeeping. Published "
            "cross sections remain from E214M1."
        ),
    )
    p.add_argument(
        "--pass2-file",
        default=str(here.parent / "output" / "csvs" / "dvcs_pass2_analysis.csv"),
        help="Final pass-2 analysis CSV after main_systematics has materialized authoritative systematics.",
    )
    p.add_argument("--target-ebeam", type=float, default=TARGET_EBEAM_GEV)
    p.add_argument("--workers", type=int, default=1, help="Reserved for future KM15 multiprocessing; current first pass evaluates serially for model safety")
    p.add_argument(
        "--model-phi-step-deg",
        type=float,
        default=15.0,
        help=(
            "Phi spacing in degrees for BH/KM15 presentation curves "
            "(default: 15 deg). Curves always include exactly 0 and 360 deg."
        ),
    )
    p.add_argument("--force-km15", action="store_true")

    p.add_argument("--match-dxb", type=float, default=DEFAULT_MATCH_DXB)
    p.add_argument("--match-dq2", type=float, default=DEFAULT_MATCH_DQ2)
    p.add_argument("--match-dt", type=float, default=DEFAULT_MATCH_DT)
    p.add_argument("--match-dphi", type=float, default=DEFAULT_MATCH_DPHI)
    p.add_argument(
        "--panel-y-scale",
        choices=["panel", "row", "page"],
        default="row",
        help=(
            "Y-axis synchronization for multipanel cross-section canvases: "
            "'panel' gives each subplot its own range, 'row' gives each 4-panel "
            "row one common range (default), and 'page' forces all 12 panels "
            "to share one range."
        ),
    )

    return p
#enddef


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    global MODEL_CURVE_PHI_STEP_DEG
    global PANEL_Y_SCALE_MODE

    MODEL_CURVE_PHI_STEP_DEG = float(args.model_phi_step_deg)
    PANEL_Y_SCALE_MODE = str(args.panel_y_scale)
    if (
        not np.isfinite(MODEL_CURVE_PHI_STEP_DEG)
        or MODEL_CURVE_PHI_STEP_DEG <= 0.0
        or MODEL_CURVE_PHI_STEP_DEG > 90.0
    ):
        raise ValueError(
            "--model-phi-step-deg must be finite and in the interval (0, 90]."
        )
    #endif
    print(
        f"[PLOTS] BH/KM15 model-curve phi step = "
        f"{MODEL_CURVE_PHI_STEP_DEG:g} deg (0--360 deg inclusive)"
    )
    print(
        f"[PLOTS] multipanel y-scale synchronization = {PANEL_Y_SCALE_MODE}",
        flush=True,
    )
    args.script_dir = Path(__file__).resolve().parent

    outdir = Path(args.outdir).expanduser()
    if not outdir.is_absolute():
        outdir = (args.script_dir / outdir).resolve()
    #endif
    outdir.mkdir(parents=True, exist_ok=True)

    emff_path = resolve_existing(
        Path(args.emff_script),
        [args.script_dir / "extract_emff_from_dvcs_bh.py"],
    )
    emff = load_python_module(emff_path, "emff_world_compare_backend")

    # ---------------------------------------------------------------------
    # 1. Load and canonicalize the six published measurements.
    # ---------------------------------------------------------------------
    world = load_world_data(args, emff)

    # ---------------------------------------------------------------------
    # 2. Native + common-energy KM15/BH calculations.
    # ---------------------------------------------------------------------
    km15_cache = outdir / "cache" / "km15_native_and_target.csv"
    world = evaluate_km15_world(
        world,
        emff,
        cache_path=km15_cache,
        target_ebeam=float(args.target_ebeam),
        force=args.force_km15,
    )

    # Keep the complete source/cache intact, but use the quality-filtered
    # sample for every nominal score, match, fit, table and figure below.
    world = apply_nominal_data_quality_exclusions(world)
    print(
        f"[WORLD] nominal sample after quality exclusions: {len(world):,} points",
        flush=True,
    )

    # ---------------------------------------------------------------------
    # 3. KM15-only transport prescription.
    #
    # The world-data consistency study intentionally uses KM15 alone for
    # model-assisted transport.  PARTONS/GK16 is not run here: the principal
    # external cross-checks are the near-10-GeV Hall-A settings and CLAS12
    # pass-1, while lower-energy measurements are transported with KM15.
    # ---------------------------------------------------------------------
    have_gk16 = False
    print("[TRANSPORT] KM15-only prescription; no PARTONS/GK16 calculation requested")

    # ---------------------------------------------------------------------
    # 4. QA summaries and native-kinematics model comparisons.
    # ---------------------------------------------------------------------
    dataset_summary = make_dataset_summary(world)
    model_scores = make_native_model_scores(world, have_gk16)

    # ---------------------------------------------------------------------
    # 5. Pairwise experiment matching and local model-assisted transport.
    # ---------------------------------------------------------------------
    match_cfg = MatchConfig(
        dxb=float(args.match_dxb),
        dq2=float(args.match_dq2),
        dt=float(args.match_dt),
        dphi=float(args.match_dphi),
    )
    print("[STAGE] building pairwise matched-data comparisons", flush=True)
    matches, pair_summary = build_pairwise_comparisons(world, match_cfg, have_gk16)
    print(
        f"[STAGE] pairwise matching complete: {len(matches):,} matched points",
        flush=True,
    )

    # ---------------------------------------------------------------------
    # 6. Tables plus note-style multi-panel cross-section presentation.
    # ---------------------------------------------------------------------
    print("[STAGE] running normalization study and producing output products", flush=True)
    norm_dataset, norm_metrics = save_outputs(
        world,
        dataset_summary,
        model_scores,
        matches,
        pair_summary,
        outdir,
        have_gk16,
        emff,
        args,
    )

    print_summary(
        dataset_summary,
        model_scores,
        pair_summary,
        have_gk16,
        norm_dataset=norm_dataset,
        norm_metrics=norm_metrics,
    )
    print(f"\n[OUTPUT] {outdir}")
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif
