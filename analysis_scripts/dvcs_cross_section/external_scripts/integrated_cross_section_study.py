#!/usr/bin/env python3
"""
integrated_cross_section_study_v2.py

Projection and integration study for the DVCS pass-2 CSV.

Major safeguards in this version
--------------------------------

1. Differential projections use the appropriate hidden-bin widths:

     d sigma / dxB   = sum sigma_i dQ2_i dt_i dphi_i
     d sigma / dQ2   = sum sigma_i dxB_i dt_i dphi_i
     d sigma / d|t|  = sum sigma_i dxB_i dQ2_i dphi_i
     d sigma / dphi  = sum sigma_i dxB_i dQ2_i dt_i

   Phi widths are converted from degrees to radians unless --phi-degrees is
   explicitly requested.

2. Every ratio is constructed on common row support. The numerator and
   denominator are integrated over exactly the same multidimensional CSV rows.

3. Uncertainties are kept separate:

     statistical:               independent, added in quadrature
     CSV point-to-point:         independent, added in quadrature
     CSV scale:                  correlated through one nuisance parameter,
                                 delta_scale = |sum_i s_i sigma_i dV_i|
     optional extra p2p fraction independent row by row

4. Closure diagnostics compare the total cross section obtained from every
   physical projection with the direct four-dimensional integral.

5. Detector-angle diagnostics print, for every derived theta bin and period,
   the row count, integrated four-dimensional volume and summed averaging
   weight.

6. The default KM15 projection mode is row-sum. Representative mode is retained
   only as a fast diagnostic and prints a warning because evaluating a nonlinear
   model once at averaged hidden kinematics is not an integration.

7. Middle-panel ratio statistical errors are explicitly labeled approximate,
   because the CSV does not contain the covariance between an individual period
   and the combined 10.6-GeV result. Right-panel ratio errors include the exact
   covariance induced by dividing either input by their own two-point weighted
   mean, assuming the Fa18 Inb and scaled Sp19 measurements are independent.

Outputs
-------

Seven 1x3 canvases:

  integrated_xB_dependence.png
  integrated_Q2_dependence.png
  integrated_t_dependence.png
  integrated_phi_dependence.png
  integrated_e_theta_dependence.png
  integrated_p_theta_dependence.png
  integrated_g_theta_dependence.png

One 2x2 KM15 canvas unless disabled:

  km15_kinematic_dependences.png
"""

import argparse
import concurrent.futures
import math
import os
import re
import traceback
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# Constants.
# -----------------------------------------------------------------------------

FA18_INB_EBEAM_GEV = 10.6040
SP19_INB_EBEAM_GEV = 10.1998

RATIO_YMIN = 0.0
RATIO_YMAX = 1.6

DEFAULT_MAX_WORKERS = 5
DEFAULT_THETA_BINS = 7
DEFAULT_THETA_BINNING_PERIOD = "10.6 GeV"
DEFAULT_CLOSURE_TOLERANCE = 1.0e-10

SP19_TO_FA18_KM15_SCALE_COLUMN = "_sp19_to_fa18_km15_scale"
SP19_TO_FA18_BH_SCALE_COLUMN = "_sp19_to_fa18_bh_scale"

DEFAULT_KM15_PREDICTION_OUTPUT_NAME = "km15_kinematic_dependences.png"
DEFAULT_KM15_PREDICTION_PERIOD = "10.6 GeV"
DEFAULT_KM15_PREDICTION_EBEAM_GEV = FA18_INB_EBEAM_GEV
DEFAULT_KM15_PREDICTION_MODE = "row-sum"
KM15_PREDICTION_LABEL = "KM15, 10.604 GeV"

RUN_PERIODS_106 = [
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
]

LEFT_SERIES = [
    "10.6 GeV",
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
]

MIDDLE_SERIES = [
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
]

RIGHT_FA_LABEL = "Fa18 Inb, 10.604 GeV"
RIGHT_SP_LABEL = "Sp19 Inb -> 10.604 GeV, KM15"
RIGHT_SERIES_LABELS = [RIGHT_FA_LABEL, RIGHT_SP_LABEL]

SERIES_STYLE = {
    "10.6 GeV": dict(marker="o", linestyle="-", color="black"),
    "Sp18 Inb": dict(marker="D", linestyle="-", color="tab:green"),
    "Sp18 Out": dict(marker="P", linestyle="-", color="tab:purple"),
    "Fa18 Inb": dict(marker="^", linestyle="-", color="tab:blue"),
    "Fa18 Out": dict(marker="v", linestyle="-", color="tab:orange"),
    "Sp19 Inb": dict(marker="s", linestyle="-", color="tab:red"),
    RIGHT_FA_LABEL: dict(marker="^", linestyle="-", color="tab:blue"),
    RIGHT_SP_LABEL: dict(marker="s", linestyle="-", color="tab:red"),
    KM15_PREDICTION_LABEL: dict(marker="o", linestyle="-", color="black"),
}

PROJECTIONS = {
    "xB": {
        "min_col": "xBmin",
        "max_col": "xBmax",
        "avg_prefix": "xBavg",
        "xlabel": r"$x_B$",
        "ylabel": r"$d\sigma/dx_B$  (pb)",
        "title": r"$x_B$ dependence",
        "integrate_widths": ["Q2", "t", "phi"],
        "outfile_tag": "xB",
        "derived_theta": False,
    },
    "Q2": {
        "min_col": "Q2min",
        "max_col": "Q2max",
        "avg_prefix": "Q2avg",
        "xlabel": r"$Q^2$  (GeV$^2$)",
        "ylabel": r"$d\sigma/dQ^2$  (pb/GeV$^2$)",
        "title": r"$Q^2$ dependence",
        "integrate_widths": ["xB", "t", "phi"],
        "outfile_tag": "Q2",
        "derived_theta": False,
    },
    "t": {
        "min_col": "t_abs_min",
        "max_col": "t_abs_max",
        "avg_prefix": "t_abs_avg",
        "xlabel": r"$|t|$  (GeV$^2$)",
        "ylabel": r"$d\sigma/d|t|$  (pb/GeV$^2$)",
        "title": r"$|t|$ dependence",
        "integrate_widths": ["xB", "Q2", "phi"],
        "outfile_tag": "t",
        "derived_theta": False,
    },
    "phi": {
        "min_col": "phimin",
        "max_col": "phimax",
        "avg_prefix": "phiavg",
        "xlabel": r"$\phi$  (deg)",
        "ylabel": r"$d\sigma/d\phi$  (pb/rad)",
        "title": r"$\phi$ dependence",
        "integrate_widths": ["xB", "Q2", "t"],
        "outfile_tag": "phi",
        "derived_theta": False,
    },
    "e_theta": {
        "min_col": "_bin_e_theta_min",
        "max_col": "_bin_e_theta_max",
        "avg_prefix": "e_theta",
        "theta_source_prefix": "e_theta",
        "xlabel": r"$\theta_e$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": "electron polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "e_theta",
        "derived_theta": True,
    },
    "p_theta": {
        "min_col": "_bin_p_theta_min",
        "max_col": "_bin_p_theta_max",
        "avg_prefix": "p_theta",
        "theta_source_prefix": "p_theta",
        "xlabel": r"$\theta_p$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": "proton polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "p_theta",
        "derived_theta": True,
    },
    "g_theta": {
        "min_col": "_bin_g_theta_min",
        "max_col": "_bin_g_theta_max",
        "avg_prefix": "g_theta",
        "theta_source_prefix": "g_theta",
        "xlabel": r"$\theta_\gamma$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": "photon polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "g_theta",
        "derived_theta": True,
    },
}

DATA_PROJECTION_ORDER = [
    "xB",
    "Q2",
    "t",
    "phi",
    "e_theta",
    "p_theta",
    "g_theta",
]

PHYSICAL_PROJECTION_ORDER = ["xB", "Q2", "t", "phi"]
THETA_PROJECTION_ORDER = ["e_theta", "p_theta", "g_theta"]
KM15_PREDICTION_PROJECTION_ORDER = ["xB", "Q2", "t", "phi"]

FLOAT_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


# -----------------------------------------------------------------------------
# Dataclasses.
# -----------------------------------------------------------------------------

@dataclass
class Point:
    key: Tuple[float, float]
    x: float
    y: float
    stat: float
    p2p: float
    scale: float
    n_rows: int = 0
    volume: float = math.nan


@dataclass
class RepresentativeKinematics:
    xB: float
    Q2: float
    t_abs: float
    phi: float
    integration_weight_sum: float
    n_rows: int


@dataclass
class ModelDiagnostics:
    enabled: bool = False
    km15_correction_min: float = math.inf
    km15_correction_max: float = -math.inf
    km15_correction_sum: float = 0.0
    km15_correction_count: int = 0
    bh_correction_min: float = math.inf
    bh_correction_max: float = -math.inf
    bh_correction_sum: float = 0.0
    bh_correction_count: int = 0
    bh_failure_message: Optional[str] = None
    sp19_unique_correction_count: int = 0
    km15_prediction_count: int = 0
    km15_prediction_failure_count: int = 0
    km15_prediction_failure_message: Optional[str] = None


@dataclass
class ModelContext:
    enabled: bool
    g: object = None
    th_km15: object = None
    km15_prediction_cache: Dict[Tuple[float, float, float, float, float], float] = field(default_factory=dict)


# -----------------------------------------------------------------------------
# Worker globals.
# -----------------------------------------------------------------------------

_WORKER_GEPARD = None
_WORKER_TH_KM15 = None


def worker_get_model_context() -> ModelContext:
    global _WORKER_GEPARD
    global _WORKER_TH_KM15

    if _WORKER_GEPARD is None:
        import gepard as g
        _WORKER_GEPARD = g
    # endif

    if _WORKER_TH_KM15 is None:
        from gepard.fits import th_KM15
        _WORKER_TH_KM15 = th_KM15
    # endif

    return ModelContext(enabled=True, g=_WORKER_GEPARD, th_km15=_WORKER_TH_KM15)


# -----------------------------------------------------------------------------
# Gepard helpers.
# -----------------------------------------------------------------------------

def import_gepard():
    try:
        import gepard as g
    except ImportError as exc:
        raise RuntimeError(
            "Could not import Gepard. Install it in the Python environment used "
            "to run this script."
        ) from exc
    # endtry
    return g


def import_km15():
    try:
        from gepard.fits import th_KM15
    except ImportError as exc:
        raise RuntimeError("Could not import th_KM15 from gepard.fits.") from exc
    # endtry
    return th_KM15


def make_model_context(enabled: bool) -> ModelContext:
    if not enabled:
        return ModelContext(enabled=False)
    # endif
    return ModelContext(enabled=True, g=import_gepard(), th_km15=import_km15())


def make_gepard_xs_point(g, xB: float, Q2: float, t_abs: float, phi_deg_trento: float, ebeam: float):
    from gepard.kinematics import prepare

    pt = g.DataPoint(
        xB=float(xB),
        Q2=float(Q2),
        t=-abs(float(t_abs)),
        phi=math.radians(float(phi_deg_trento)),
        frame="Trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(ebeam),
        in1charge=-1,
        in1polarization=0,
        observable="XS",
    )

    if hasattr(pt, "to_conventions"):
        pt.to_conventions()
    # endif

    prepare(pt)
    return pt


def km15_xs(model_context: ModelContext, xB: float, Q2: float, t_abs: float, phi_deg_trento: float, ebeam: float) -> float:
    pt = make_gepard_xs_point(
        g=model_context.g,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi_deg_trento,
        ebeam=ebeam,
    )
    return float(model_context.th_km15.predict(pt))


def km15_xs_cached(model_context: ModelContext, xB: float, Q2: float, t_abs: float, phi_deg_trento: float, ebeam: float) -> float:
    key = (
        round(float(xB), 8),
        round(float(Q2), 8),
        round(float(t_abs), 8),
        round(float(phi_deg_trento), 8),
        round(float(ebeam), 8),
    )

    if key not in model_context.km15_prediction_cache:
        model_context.km15_prediction_cache[key] = km15_xs(
            model_context=model_context,
            xB=xB,
            Q2=Q2,
            t_abs=t_abs,
            phi_deg_trento=phi_deg_trento,
            ebeam=ebeam,
        )
    # endif

    return model_context.km15_prediction_cache[key]


def bh_xs(model_context: ModelContext, xB: float, Q2: float, t_abs: float, phi_deg_trento: float, ebeam: float) -> float:
    pt = make_gepard_xs_point(
        g=model_context.g,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi_deg_trento,
        ebeam=ebeam,
    )

    if hasattr(model_context.th_km15, "PreFacSigma") and hasattr(model_context.th_km15, "TBH2unp"):
        return float(model_context.th_km15.PreFacSigma(pt) * model_context.th_km15.TBH2unp(pt))
    # endif

    available = [
        name
        for name in dir(model_context.th_km15)
        if "BH" in name or "bh" in name or "PreFac" in name or "prefac" in name
    ]
    raise RuntimeError(
        "The loaded KM15 object does not expose PreFacSigma(pt) and TBH2unp(pt). "
        f"Candidate attributes: {available}"
    )


def correction_cache_key(xB: float, Q2: float, t_abs: float, phi: float) -> Tuple[float, float, float, float]:
    return (
        round(float(xB), 8),
        round(float(Q2), 8),
        round(float(t_abs), 8),
        round(float(phi), 8),
    )


def update_model_correction_diagnostics(diagnostics: ModelDiagnostics, correction: float, model_name: str) -> None:
    if not np.isfinite(correction):
        return
    # endif

    if model_name == "KM15":
        diagnostics.km15_correction_min = min(diagnostics.km15_correction_min, correction)
        diagnostics.km15_correction_max = max(diagnostics.km15_correction_max, correction)
        diagnostics.km15_correction_sum += correction
        diagnostics.km15_correction_count += 1
    elif model_name == "BH":
        diagnostics.bh_correction_min = min(diagnostics.bh_correction_min, correction)
        diagnostics.bh_correction_max = max(diagnostics.bh_correction_max, correction)
        diagnostics.bh_correction_sum += correction
        diagnostics.bh_correction_count += 1
    else:
        raise ValueError(f"Unknown model_name={model_name!r}")
    # endif


def km15_and_bh_sp19_to_fa18_corrections_direct(
    xB: float,
    Q2: float,
    t_abs: float,
    phi: float,
    compute_bh_diagnostic: bool,
    model_context: ModelContext,
) -> Tuple[float, float, Optional[str]]:
    km15_target = km15_xs(model_context, xB, Q2, t_abs, phi, FA18_INB_EBEAM_GEV)
    km15_source = km15_xs(model_context, xB, Q2, t_abs, phi, SP19_INB_EBEAM_GEV)

    if not np.isfinite(km15_target) or not np.isfinite(km15_source) or km15_source == 0.0:
        return math.nan, math.nan, None
    # endif

    km15_corr = km15_target / km15_source

    if not compute_bh_diagnostic:
        return km15_corr, math.nan, None
    # endif

    try:
        bh_target = bh_xs(model_context, xB, Q2, t_abs, phi, FA18_INB_EBEAM_GEV)
        bh_source = bh_xs(model_context, xB, Q2, t_abs, phi, SP19_INB_EBEAM_GEV)

        if np.isfinite(bh_target) and np.isfinite(bh_source) and bh_source != 0.0:
            return km15_corr, bh_target / bh_source, None
        # endif

        return km15_corr, math.nan, None

    except Exception as exc:
        return km15_corr, math.nan, repr(exc)
    # endtry


def worker_sp19_correction(task: Tuple[int, float, float, float, float, bool]) -> Tuple[int, float, float, Optional[str], Optional[str]]:
    task_id, xB, Q2, t_abs, phi, compute_bh_diagnostic = task
    try:
        model_context = worker_get_model_context()
        km15_corr, bh_corr, bh_failure = km15_and_bh_sp19_to_fa18_corrections_direct(
            xB=xB,
            Q2=Q2,
            t_abs=t_abs,
            phi=phi,
            compute_bh_diagnostic=compute_bh_diagnostic,
            model_context=model_context,
        )
        return task_id, km15_corr, bh_corr, bh_failure, None
    except Exception:
        return task_id, math.nan, math.nan, None, traceback.format_exc()
    # endtry


def worker_km15_prediction(task: Tuple[int, float, float, float, float, float]) -> Tuple[int, float, Optional[str]]:
    task_id, xB, Q2, t_abs, phi, ebeam = task
    try:
        model_context = worker_get_model_context()
        value = km15_xs(model_context, xB, Q2, t_abs, phi, ebeam)
        return task_id, value, None
    except Exception:
        return task_id, math.nan, traceback.format_exc()
    # endtry


# -----------------------------------------------------------------------------
# Parsing and CSV helpers.
# -----------------------------------------------------------------------------

def parse_tuple3(value) -> Tuple[float, float, float]:
    if value is None:
        return math.nan, math.nan, math.nan
    # endif

    if isinstance(value, float) and math.isnan(value):
        return math.nan, math.nan, math.nan
    # endif

    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value), 0.0, 0.0
    # endif

    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null"}:
        return math.nan, math.nan, math.nan
    # endif

    numbers = FLOAT_PATTERN.findall(text)
    if len(numbers) == 0:
        return math.nan, math.nan, math.nan
    # endif

    parsed = [float(x) for x in numbers]
    while len(parsed) < 3:
        parsed.append(0.0)
    # endwhile

    return parsed[0], parsed[1], parsed[2]


def parse_scalar_from_cell(value) -> float:
    if value is None:
        return math.nan
    # endif

    if isinstance(value, float) and math.isnan(value):
        return math.nan
    # endif

    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    # endif

    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null"}:
        return math.nan
    # endif

    numbers = FLOAT_PATTERN.findall(text)
    if len(numbers) == 0:
        return math.nan
    # endif

    return float(numbers[0])


def row_average_or_midpoint(row: pd.Series, avg_col: str, min_col: str, max_col: str) -> float:
    if avg_col in row.index:
        avg_value = parse_scalar_from_cell(row[avg_col])
        if np.isfinite(avg_value):
            return avg_value
        # endif
    # endif

    lo = parse_scalar_from_cell(row[min_col])
    hi = parse_scalar_from_cell(row[max_col])
    if np.isfinite(lo) and np.isfinite(hi):
        return 0.5 * (lo + hi)
    # endif

    return math.nan


def cross_section_column(period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {period}, unpol"


def point_to_point_total_column() -> str:
    return "Syst. err (point-to-point total)"


def total_scale_column(period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {period}, unpol, total scale sys"


def theta_average_column(theta_prefix: str, period: str) -> str:
    return f"{theta_prefix}, {period}"


def require_columns(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError("The input CSV is missing required columns:\n  " + "\n  ".join(missing))
    # endif


def projection_columns_exist(df: pd.DataFrame, projection: str) -> bool:
    info = PROJECTIONS[projection]
    return str(info["min_col"]) in df.columns and str(info["max_col"]) in df.columns


def add_width_columns(df: pd.DataFrame, phi_degrees: bool) -> pd.DataFrame:
    out = df.copy()
    out["_width_xB"] = pd.to_numeric(out["xBmax"], errors="coerce") - pd.to_numeric(out["xBmin"], errors="coerce")
    out["_width_Q2"] = pd.to_numeric(out["Q2max"], errors="coerce") - pd.to_numeric(out["Q2min"], errors="coerce")
    out["_width_t"] = pd.to_numeric(out["t_abs_max"], errors="coerce") - pd.to_numeric(out["t_abs_min"], errors="coerce")
    phi_width_deg = pd.to_numeric(out["phimax"], errors="coerce") - pd.to_numeric(out["phimin"], errors="coerce")
    out["_width_phi"] = phi_width_deg if phi_degrees else np.deg2rad(phi_width_deg)
    out["_width_4d"] = out["_width_xB"] * out["_width_Q2"] * out["_width_t"] * out["_width_phi"]
    return out


def add_derived_theta_bin_columns(df: pd.DataFrame, theta_binning_period: str, theta_bins: int) -> pd.DataFrame:
    if theta_bins <= 0:
        raise ValueError("--theta-bins must be positive.")
    # endif

    out = df.copy()

    for projection in THETA_PROJECTION_ORDER:
        info = PROJECTIONS[projection]
        source_col = theta_average_column(str(info["theta_source_prefix"]), theta_binning_period)
        if source_col not in out.columns:
            raise KeyError(
                f"Missing angle-average column needed for {projection} binning: {source_col}"
            )
        # endif

        values = out[source_col].apply(parse_scalar_from_cell).astype(float)
        finite = values[np.isfinite(values)]
        if len(finite) == 0:
            raise RuntimeError(f"No finite values found in {source_col}.")
        # endif

        vmin = float(finite.min())
        vmax = float(finite.max())
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            raise RuntimeError(f"Invalid range for {source_col}: min={vmin}, max={vmax}")
        # endif

        edges = np.linspace(vmin, vmax, theta_bins + 1)
        min_col = str(info["min_col"])
        max_col = str(info["max_col"])
        out[min_col] = math.nan
        out[max_col] = math.nan

        for idx, val in values.items():
            if not np.isfinite(val):
                continue
            # endif

            bin_index = int(np.searchsorted(edges, val, side="right") - 1)
            bin_index = max(0, min(theta_bins - 1, bin_index))
            out.at[idx, min_col] = float(edges[bin_index])
            out.at[idx, max_col] = float(edges[bin_index + 1])
        # endfor

        print(
            f"Derived {projection} bins from {source_col}: "
            f"{theta_bins} bins from {vmin:.6g} to {vmax:.6g} deg"
        )
    # endfor

    return out


def weight_for_row(row: pd.Series, integrate_widths: Sequence[str], no_width_weighting: bool) -> float:
    if no_width_weighting:
        return 1.0
    # endif

    weight = 1.0
    for axis in integrate_widths:
        width = parse_scalar_from_cell(row.get(f"_width_{axis}", math.nan))
        if not np.isfinite(width) or width <= 0.0:
            return math.nan
        # endif
        weight *= width
    # endfor
    return weight


def candidate_yield_columns_for_period(df: pd.DataFrame, period: str) -> List[str]:
    period_lower = period.lower()
    priority_phrases = [
        "pi0 corrected",
        "pi0-corrected",
        "pi0_corrected",
        "signal yield",
        "acceptance corrected yield",
        "acceptance-corrected yield",
        "acceptance_corrected_yield",
        "unfolded yield",
        "corrected yield",
        "yield",
        "counts",
        "num events",
        "numevents",
    ]

    candidates: List[str] = []
    for phrase in priority_phrases:
        for col in df.columns:
            col_lower = str(col).lower()
            if period_lower in col_lower and "unpol" in col_lower and phrase in col_lower and col not in candidates:
                candidates.append(col)
            # endif
        # endfor
    # endfor
    return candidates


def event_weight_for_average(row: pd.Series, candidate_columns: Sequence[str]) -> float:
    for col in candidate_columns:
        val = parse_scalar_from_cell(row[col])
        if np.isfinite(val) and val > 0.0:
            return float(val)
        # endif
    # endfor
    return math.nan


def average_x_for_group(group: pd.DataFrame, projection_info: Dict[str, object], period: str, yield_weight_columns: Sequence[str]) -> float:
    avg_col = f"{projection_info['avg_prefix']}, {period}"
    min_col = str(projection_info["min_col"])
    max_col = str(projection_info["max_col"])
    xlo = float(pd.to_numeric(group[min_col], errors="coerce").iloc[0])
    xhi = float(pd.to_numeric(group[max_col], errors="coerce").iloc[0])
    midpoint = 0.5 * (xlo + xhi)

    if avg_col not in group.columns:
        return midpoint
    # endif

    numerator = 0.0
    denominator = 0.0
    finite_values: List[float] = []

    for _, row in group.iterrows():
        xavg = parse_scalar_from_cell(row[avg_col])
        if not np.isfinite(xavg):
            continue
        # endif

        finite_values.append(xavg)
        weight = event_weight_for_average(row, yield_weight_columns)
        if np.isfinite(weight) and weight > 0.0:
            numerator += weight * xavg
            denominator += weight
        # endif
    # endfor

    if denominator > 0.0:
        return numerator / denominator
    # endif
    if finite_values:
        return float(np.mean(np.asarray(finite_values, dtype=float)))
    # endif
    return midpoint


def finite_cross_section_mask(df: pd.DataFrame, period: str) -> pd.Series:
    col = cross_section_column(period)
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    # endif
    return df[col].apply(lambda cell: np.isfinite(parse_tuple3(cell)[0]))


def common_support_mask(
    df: pd.DataFrame,
    periods: Sequence[str],
    require_sp19_scale: bool = False,
) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    for period in periods:
        mask &= finite_cross_section_mask(df, period)
    # endfor

    if require_sp19_scale:
        if SP19_TO_FA18_KM15_SCALE_COLUMN not in df.columns:
            return pd.Series(False, index=df.index)
        # endif
        mask &= pd.to_numeric(df[SP19_TO_FA18_KM15_SCALE_COLUMN], errors="coerce").apply(np.isfinite)
    # endif

    return mask


# -----------------------------------------------------------------------------
# Model precomputation.
# -----------------------------------------------------------------------------

def run_sp19_correction_tasks(
    tasks: List[Tuple[int, float, float, float, float, bool]],
    max_workers: int,
    no_parallel_km15: bool,
) -> Dict[int, Tuple[float, float, Optional[str], Optional[str]]]:
    results: Dict[int, Tuple[float, float, Optional[str], Optional[str]]] = {}
    if not tasks:
        return results
    # endif

    if no_parallel_km15 or max_workers <= 1:
        model_context = make_model_context(enabled=True)
        for task_id, xB, Q2, t_abs, phi, compute_bh in tasks:
            try:
                km15_corr, bh_corr, bh_failure = km15_and_bh_sp19_to_fa18_corrections_direct(
                    xB, Q2, t_abs, phi, compute_bh, model_context
                )
                results[task_id] = (km15_corr, bh_corr, bh_failure, None)
            except Exception:
                results[task_id] = (math.nan, math.nan, None, traceback.format_exc())
            # endtry
        # endfor
        return results
    # endif

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for task_id, km15_corr, bh_corr, bh_failure, fatal_error in executor.map(worker_sp19_correction, tasks):
            results[task_id] = (km15_corr, bh_corr, bh_failure, fatal_error)
        # endfor
    # endwith
    return results


def precompute_sp19_to_fa18_km15_scales(
    df: pd.DataFrame,
    enabled: bool,
    compute_bh_diagnostic: bool,
    max_workers: int,
    no_parallel_km15: bool,
) -> Tuple[pd.DataFrame, ModelDiagnostics]:
    out = df.copy()
    diagnostics = ModelDiagnostics(enabled=enabled)

    if not enabled:
        out[SP19_TO_FA18_KM15_SCALE_COLUMN] = 1.0
        out[SP19_TO_FA18_BH_SCALE_COLUMN] = math.nan
        return out, diagnostics
    # endif

    xs_col = cross_section_column("Sp19 Inb")
    if xs_col not in out.columns:
        out[SP19_TO_FA18_KM15_SCALE_COLUMN] = math.nan
        out[SP19_TO_FA18_BH_SCALE_COLUMN] = math.nan
        return out, diagnostics
    # endif

    key_by_index: Dict[int, Tuple[float, float, float, float]] = {}
    task_by_key: Dict[Tuple[float, float, float, float], Tuple[int, float, float, float, float, bool]] = {}
    task_id_by_key: Dict[Tuple[float, float, float, float], int] = {}
    next_task_id = 0

    for idx, row in out.iterrows():
        value, _, _ = parse_tuple3(row[xs_col])
        if not np.isfinite(value):
            continue
        # endif

        xB = row_average_or_midpoint(row, "xBavg, Sp19 Inb", "xBmin", "xBmax")
        Q2 = row_average_or_midpoint(row, "Q2avg, Sp19 Inb", "Q2min", "Q2max")
        t_abs = row_average_or_midpoint(row, "t_abs_avg, Sp19 Inb", "t_abs_min", "t_abs_max")
        phi = row_average_or_midpoint(row, "phiavg, Sp19 Inb", "phimin", "phimax")

        if not all(np.isfinite(v) for v in [xB, Q2, t_abs, phi]):
            continue
        # endif

        key = correction_cache_key(xB, Q2, t_abs, phi)
        key_by_index[int(idx)] = key
        if key not in task_by_key:
            task_by_key[key] = (
                next_task_id,
                float(xB),
                float(Q2),
                float(t_abs),
                float(phi),
                bool(compute_bh_diagnostic),
            )
            task_id_by_key[key] = next_task_id
            next_task_id += 1
        # endif
    # endfor

    tasks = list(task_by_key.values())
    print("\nPrecomputing Sp19 -> 10.604 GeV KM15 correction factors:")
    print(f"  unique kinematic points = {len(tasks)}")
    print(f"  parallel                = {not no_parallel_km15 and max_workers > 1}")
    print(f"  max workers             = {max_workers}")
    print(f"  BH diagnostic           = {compute_bh_diagnostic}")

    results = run_sp19_correction_tasks(tasks, max_workers, no_parallel_km15)
    scale_by_key: Dict[Tuple[float, float, float, float], float] = {}
    bh_by_key: Dict[Tuple[float, float, float, float], float] = {}

    for key, task_id in task_id_by_key.items():
        km15_corr, bh_corr, bh_failure, fatal_error = results.get(
            task_id, (math.nan, math.nan, None, "missing result")
        )
        if fatal_error is not None and diagnostics.km15_prediction_failure_message is None:
            diagnostics.km15_prediction_failure_message = fatal_error
        # endif
        if bh_failure is not None and diagnostics.bh_failure_message is None:
            diagnostics.bh_failure_message = bh_failure
        # endif

        scale_by_key[key] = km15_corr
        bh_by_key[key] = bh_corr
        update_model_correction_diagnostics(diagnostics, km15_corr, "KM15")
        update_model_correction_diagnostics(diagnostics, bh_corr, "BH")
    # endfor

    diagnostics.sp19_unique_correction_count = len(tasks)
    out[SP19_TO_FA18_KM15_SCALE_COLUMN] = math.nan
    out[SP19_TO_FA18_BH_SCALE_COLUMN] = math.nan

    for idx, key in key_by_index.items():
        out.at[idx, SP19_TO_FA18_KM15_SCALE_COLUMN] = scale_by_key.get(key, math.nan)
        out.at[idx, SP19_TO_FA18_BH_SCALE_COLUMN] = bh_by_key.get(key, math.nan)
    # endfor

    return out, diagnostics


def run_km15_prediction_tasks(
    tasks: List[Tuple[int, float, float, float, float, float]],
    max_workers: int,
    no_parallel_km15: bool,
) -> Dict[int, Tuple[float, Optional[str]]]:
    results: Dict[int, Tuple[float, Optional[str]]] = {}
    if not tasks:
        return results
    # endif

    if no_parallel_km15 or max_workers <= 1:
        model_context = make_model_context(enabled=True)
        for task_id, xB, Q2, t_abs, phi, ebeam in tasks:
            try:
                value = km15_xs_cached(model_context, xB, Q2, t_abs, phi, ebeam)
                results[task_id] = (value, None)
            except Exception:
                results[task_id] = (math.nan, traceback.format_exc())
            # endtry
        # endfor
        return results
    # endif

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for task_id, value, fatal_error in executor.map(worker_km15_prediction, tasks):
            results[task_id] = (value, fatal_error)
        # endfor
    # endwith
    return results


# -----------------------------------------------------------------------------
# Integration.
# -----------------------------------------------------------------------------

def integrate_group_for_period(
    group: pd.DataFrame,
    projection: str,
    period: str,
    no_width_weighting: bool,
    additional_p2p_frac: float,
    apply_sp19_to_fa18_scaling: bool,
) -> Optional[Point]:
    info = PROJECTIONS[projection]
    xs_col = cross_section_column(period)
    yield_columns = candidate_yield_columns_for_period(group, period)

    y_sum = 0.0
    stat2_sum = 0.0
    p2p2_sum = 0.0
    correlated_scale_sum = 0.0
    volume_sum = 0.0
    n_used = 0

    for _, row in group.iterrows():
        value, stat, tuple_sys = parse_tuple3(row[xs_col])
        if not np.isfinite(value):
            continue
        # endif

        weight = weight_for_row(row, list(info["integrate_widths"]), no_width_weighting)
        if not np.isfinite(weight):
            continue
        # endif

        if not np.isfinite(stat):
            stat = 0.0
        # endif

        csv_p2p = parse_scalar_from_cell(row.get(point_to_point_total_column(), math.nan))
        if not np.isfinite(csv_p2p):
            csv_p2p = tuple_sys if np.isfinite(tuple_sys) else 0.0
        # endif

        scale_frac = parse_scalar_from_cell(row.get(total_scale_column(period), math.nan))
        if not np.isfinite(scale_frac):
            scale_frac = 0.0
        # endif

        model_scale = 1.0
        if apply_sp19_to_fa18_scaling:
            model_scale = parse_scalar_from_cell(row.get(SP19_TO_FA18_KM15_SCALE_COLUMN, math.nan))
            if not np.isfinite(model_scale):
                continue
            # endif
        # endif

        value *= model_scale
        stat *= abs(model_scale)
        csv_p2p *= abs(model_scale)

        extra_p2p_abs = abs(additional_p2p_frac * value)
        row_p2p_abs = math.hypot(csv_p2p, extra_p2p_abs)

        y_sum += value * weight
        stat2_sum += (stat * weight) ** 2
        p2p2_sum += (row_p2p_abs * weight) ** 2
        correlated_scale_sum += scale_frac * value * weight
        volume_sum += weight
        n_used += 1
    # endfor

    if n_used == 0:
        return None
    # endif

    x = average_x_for_group(group, info, period, yield_columns)
    key = (
        float(group[str(info["min_col"])].iloc[0]),
        float(group[str(info["max_col"])].iloc[0]),
    )

    return Point(
        key=key,
        x=x,
        y=y_sum,
        stat=math.sqrt(stat2_sum),
        p2p=math.sqrt(p2p2_sum),
        scale=abs(correlated_scale_sum),
        n_rows=n_used,
        volume=volume_sum,
    )


def integrated_points_for_period(
    df: pd.DataFrame,
    projection: str,
    period: str,
    no_width_weighting: bool,
    additional_p2p_frac: float = 0.0,
    apply_sp19_to_fa18_scaling: bool = False,
    support_mask: Optional[pd.Series] = None,
) -> List[Point]:
    if not projection_columns_exist(df, projection):
        raise KeyError(f"Missing projection columns for {projection}")
    # endif

    xs_col = cross_section_column(period)
    if xs_col not in df.columns:
        return []
    # endif

    work = df
    if support_mask is not None:
        work = df.loc[support_mask].copy()
    # endif

    info = PROJECTIONS[projection]
    group_cols = [str(info["min_col"]), str(info["max_col"])]
    points: List[Point] = []

    for _, group in work.sort_values(group_cols).groupby(group_cols, sort=True, dropna=True):
        point = integrate_group_for_period(
            group=group,
            projection=projection,
            period=period,
            no_width_weighting=no_width_weighting,
            additional_p2p_frac=additional_p2p_frac,
            apply_sp19_to_fa18_scaling=apply_sp19_to_fa18_scaling,
        )
        if point is not None:
            points.append(point)
        # endif
    # endfor

    points.sort(key=lambda p: p.x)
    return points


def direct_total_for_period(
    df: pd.DataFrame,
    period: str,
    no_width_weighting: bool,
    additional_p2p_frac: float,
    apply_sp19_to_fa18_scaling: bool = False,
    support_mask: Optional[pd.Series] = None,
) -> Point:
    xs_col = cross_section_column(period)
    work = df if support_mask is None else df.loc[support_mask]

    y_sum = 0.0
    stat2_sum = 0.0
    p2p2_sum = 0.0
    scale_sum = 0.0
    volume_sum = 0.0
    n_rows = 0

    for _, row in work.iterrows():
        value, stat, tuple_sys = parse_tuple3(row.get(xs_col, math.nan))
        if not np.isfinite(value):
            continue
        # endif

        weight = 1.0 if no_width_weighting else parse_scalar_from_cell(row.get("_width_4d", math.nan))
        if not np.isfinite(weight) or weight <= 0.0:
            continue
        # endif

        if not np.isfinite(stat):
            stat = 0.0
        # endif

        p2p = parse_scalar_from_cell(row.get(point_to_point_total_column(), math.nan))
        if not np.isfinite(p2p):
            p2p = tuple_sys if np.isfinite(tuple_sys) else 0.0
        # endif

        scale_frac = parse_scalar_from_cell(row.get(total_scale_column(period), math.nan))
        if not np.isfinite(scale_frac):
            scale_frac = 0.0
        # endif

        model_scale = 1.0
        if apply_sp19_to_fa18_scaling:
            model_scale = parse_scalar_from_cell(row.get(SP19_TO_FA18_KM15_SCALE_COLUMN, math.nan))
            if not np.isfinite(model_scale):
                continue
            # endif
        # endif

        value *= model_scale
        stat *= abs(model_scale)
        p2p *= abs(model_scale)
        p2p = math.hypot(p2p, abs(additional_p2p_frac * value))

        y_sum += value * weight
        stat2_sum += (stat * weight) ** 2
        p2p2_sum += (p2p * weight) ** 2
        scale_sum += scale_frac * value * weight
        volume_sum += weight
        n_rows += 1
    # endfor

    return Point(
        key=(math.nan, math.nan),
        x=math.nan,
        y=y_sum,
        stat=math.sqrt(stat2_sum),
        p2p=math.sqrt(p2p2_sum),
        scale=abs(scale_sum),
        n_rows=n_rows,
        volume=volume_sum,
    )


def projected_total(
    points: Sequence[Point],
    projection: str,
    no_width_weighting: bool,
    phi_degrees: bool,
) -> float:
    if no_width_weighting:
        return float(sum(p.y for p in points if np.isfinite(p.y)))
    # endif

    if bool(PROJECTIONS[projection]["derived_theta"]):
        return float(sum(p.y for p in points if np.isfinite(p.y)))
    # endif

    total = 0.0
    for point in points:
        projected_width = point.key[1] - point.key[0]

        # The phi bin edges stored in the CSV are degrees.  In the default
        # analysis the differential projection is d sigma / d phi in pb/rad,
        # because the hidden phi widths entering all other integrations are
        # converted to radians.  Closure must therefore convert the displayed
        # degree bin width back to radians before reintegration.  Only the
        # explicit --phi-degrees mode should retain degree widths here.
        if projection == "phi" and not phi_degrees:
            projected_width = float(np.deg2rad(projected_width))
        # endif

        if np.isfinite(point.y) and np.isfinite(projected_width) and projected_width > 0.0:
            total += point.y * projected_width
        # endif
    # endfor
    return total


def print_projection_closure_diagnostics(
    df: pd.DataFrame,
    periods: Sequence[str],
    no_width_weighting: bool,
    phi_degrees: bool,
    additional_p2p_frac: float,
    tolerance: float,
    strict: bool,
) -> None:
    print("\nProjection-closure diagnostics")
    print("Values compare each projected reconstruction with the direct 4D integral.")

    failures: List[str] = []

    for period in periods:
        direct = direct_total_for_period(
            df=df,
            period=period,
            no_width_weighting=no_width_weighting,
            additional_p2p_frac=additional_p2p_frac,
        )

        print(f"\n{period}: direct total = {direct.y:.12g} pb from {direct.n_rows} rows")
        for projection in DATA_PROJECTION_ORDER:
            points = integrated_points_for_period(
                df=df,
                projection=projection,
                period=period,
                no_width_weighting=no_width_weighting,
                additional_p2p_frac=additional_p2p_frac,
            )
            reconstructed = projected_total(points, projection, no_width_weighting, phi_degrees)
            denom = max(abs(direct.y), 1.0e-300)
            rel = (reconstructed - direct.y) / denom
            print(
                f"  {projection:8s}: reconstructed={reconstructed:.12g}  "
                f"relative difference={rel:+.3e}"
            )
            if abs(rel) > tolerance:
                failures.append(f"{period} {projection}: relative closure difference {rel:+.6e}")
            # endif
        # endfor
    # endfor

    if failures:
        print("\nWARNING: projection closure exceeded the requested tolerance:")
        for failure in failures:
            print(f"  {failure}")
        # endfor
        if strict:
            raise RuntimeError("Projection closure failed in strict mode.")
        # endif
    else:
        print("\nAll projection closures passed.")
    # endif


def print_theta_bin_diagnostics(
    df: pd.DataFrame,
    periods: Sequence[str],
    no_width_weighting: bool,
) -> None:
    print("\nDetector-angle bin diagnostics")
    print("Each row reports finite cross-section cells, integrated 4D volume and summed averaging weight.")

    for projection in THETA_PROJECTION_ORDER:
        info = PROJECTIONS[projection]
        group_cols = [str(info["min_col"]), str(info["max_col"])]
        print(f"\n{projection}:")

        for key, group in df.sort_values(group_cols).groupby(group_cols, sort=True, dropna=True):
            lo, hi = float(key[0]), float(key[1])
            base_volume = 0.0
            for _, row in group.iterrows():
                width = 1.0 if no_width_weighting else parse_scalar_from_cell(row.get("_width_4d", math.nan))
                if np.isfinite(width) and width > 0.0:
                    base_volume += width
                # endif
            # endfor

            print(f"  bin [{lo:.6g}, {hi:.6g}] deg; all-row volume={base_volume:.6g}")
            for period in periods:
                mask = finite_cross_section_mask(group, period)
                finite_group = group.loc[mask]
                volume = 0.0
                averaging_weight_sum = 0.0
                candidate_columns = candidate_yield_columns_for_period(group, period)

                for _, row in finite_group.iterrows():
                    width = 1.0 if no_width_weighting else parse_scalar_from_cell(row.get("_width_4d", math.nan))
                    if np.isfinite(width) and width > 0.0:
                        volume += width
                    # endif
                    avg_weight = event_weight_for_average(row, candidate_columns)
                    if np.isfinite(avg_weight) and avg_weight > 0.0:
                        averaging_weight_sum += avg_weight
                    # endif
                # endfor

                print(
                    f"    {period:10s}: rows={len(finite_group):5d}  "
                    f"volume={volume:.6g}  avg-weight-sum={averaging_weight_sum:.6g}"
                )
            # endfor
        # endfor
    # endfor


# -----------------------------------------------------------------------------
# Ratio helpers.
# -----------------------------------------------------------------------------

def ratio_points_independent_approx(numerator: Sequence[Point], denominator: Sequence[Point]) -> List[Point]:
    n_by_key = {p.key: p for p in numerator}
    d_by_key = {p.key: p for p in denominator}
    ratios: List[Point] = []

    for key in sorted(set(n_by_key) & set(d_by_key)):
        n = n_by_key[key]
        d = d_by_key[key]
        if not np.isfinite(n.y) or not np.isfinite(d.y) or d.y == 0.0:
            continue
        # endif

        r = n.y / d.y
        stat2 = 0.0
        if n.y != 0.0 and np.isfinite(n.stat):
            stat2 += (n.stat / n.y) ** 2
        # endif
        if np.isfinite(d.stat):
            stat2 += (d.stat / d.y) ** 2
        # endif

        p2p2 = 0.0
        if n.y != 0.0 and np.isfinite(n.p2p):
            p2p2 += (n.p2p / n.y) ** 2
        # endif
        if np.isfinite(d.p2p):
            p2p2 += (d.p2p / d.y) ** 2
        # endif

        scale2 = 0.0
        if n.y != 0.0 and np.isfinite(n.scale):
            scale2 += (n.scale / n.y) ** 2
        # endif
        if np.isfinite(d.scale):
            scale2 += (d.scale / d.y) ** 2
        # endif

        ratios.append(
            Point(
                key=key,
                x=n.x,
                y=r,
                stat=abs(r) * math.sqrt(stat2),
                p2p=abs(r) * math.sqrt(p2p2),
                scale=abs(r) * math.sqrt(scale2),
                n_rows=min(n.n_rows, d.n_rows),
                volume=min(n.volume, d.volume),
            )
        )
    # endfor

    return ratios


def inverse_variance_weighted_mean(a: Point, b: Point) -> Optional[Tuple[float, float, float, float]]:
    if not np.isfinite(a.y) or not np.isfinite(b.y):
        return None
    # endif
    if not np.isfinite(a.stat) or not np.isfinite(b.stat) or a.stat <= 0.0 or b.stat <= 0.0:
        return None
    # endif

    wa = 1.0 / (a.stat * a.stat)
    wb = 1.0 / (b.stat * b.stat)
    wsum = wa + wb
    mean = (wa * a.y + wb * b.y) / wsum
    alpha = wa / wsum
    beta = wb / wsum
    variance = 1.0 / wsum
    return mean, variance, alpha, beta


def ratio_to_own_two_point_weighted_mean(
    fa_points: Sequence[Point],
    sp_points: Sequence[Point],
    target_label: str,
) -> List[Point]:
    fa_by_key = {p.key: p for p in fa_points}
    sp_by_key = {p.key: p for p in sp_points}
    ratios: List[Point] = []

    for key in sorted(set(fa_by_key) & set(sp_by_key)):
        fa = fa_by_key[key]
        sp = sp_by_key[key]
        mean_info = inverse_variance_weighted_mean(fa, sp)
        if mean_info is None:
            continue
        # endif

        mean, mean_var, alpha, beta = mean_info
        if mean == 0.0:
            continue
        # endif

        target = fa if target_label == RIGHT_FA_LABEL else sp
        r = target.y / mean

        # Exact statistical covariance for p / (alpha*a + beta*b), assuming
        # independent a and b and fixed inverse-variance weights.
        if target_label == RIGHT_FA_LABEL:
            cov_target_mean = alpha * fa.stat * fa.stat
        else:
            cov_target_mean = beta * sp.stat * sp.stat
        # endif

        stat_var = (
            target.stat * target.stat / (mean * mean)
            + target.y * target.y * mean_var / (mean ** 4)
            - 2.0 * target.y * cov_target_mean / (mean ** 3)
        )
        stat_var = max(0.0, stat_var)

        # P2P and scale terms are shown conservatively with independent
        # propagation because their inter-period correlations are not encoded.
        p2p_var = 0.0
        scale_var = 0.0
        if target.y != 0.0:
            p2p_var += (target.p2p / target.y) ** 2
            scale_var += (target.scale / target.y) ** 2
        # endif

        mean_p2p = math.sqrt((alpha * fa.p2p) ** 2 + (beta * sp.p2p) ** 2)
        mean_scale = math.sqrt((alpha * fa.scale) ** 2 + (beta * sp.scale) ** 2)
        p2p_var += (mean_p2p / mean) ** 2
        scale_var += (mean_scale / mean) ** 2

        ratios.append(
            Point(
                key=key,
                x=target.x,
                y=r,
                stat=math.sqrt(stat_var),
                p2p=abs(r) * math.sqrt(p2p_var),
                scale=abs(r) * math.sqrt(scale_var),
                n_rows=min(fa.n_rows, sp.n_rows),
                volume=min(fa.volume, sp.volume),
            )
        )
    # endfor

    return ratios


# -----------------------------------------------------------------------------
# KM15 prediction projections.
# -----------------------------------------------------------------------------

def representative_kinematics_for_group(
    group: pd.DataFrame,
    projection: str,
    period: str,
    no_width_weighting: bool,
    yield_weight_columns: Sequence[str],
) -> RepresentativeKinematics:
    info = PROJECTIONS[projection]
    xs_col = cross_section_column(period)

    xB_sum = 0.0
    Q2_sum = 0.0
    t_sum = 0.0
    sin_phi_sum = 0.0
    cos_phi_sum = 0.0
    rep_weight_sum = 0.0
    integration_weight_sum = 0.0
    n_rows = 0

    for _, row in group.iterrows():
        observed_value, _, _ = parse_tuple3(row.get(xs_col, math.nan))
        if not np.isfinite(observed_value):
            continue
        # endif

        xB = row_average_or_midpoint(row, f"xBavg, {period}", "xBmin", "xBmax")
        Q2 = row_average_or_midpoint(row, f"Q2avg, {period}", "Q2min", "Q2max")
        t_abs = row_average_or_midpoint(row, f"t_abs_avg, {period}", "t_abs_min", "t_abs_max")
        phi = row_average_or_midpoint(row, f"phiavg, {period}", "phimin", "phimax")
        if not all(np.isfinite(v) for v in [xB, Q2, t_abs, phi]):
            continue
        # endif

        integration_weight = weight_for_row(row, list(info["integrate_widths"]), no_width_weighting)
        if not np.isfinite(integration_weight):
            continue
        # endif

        rep_weight = event_weight_for_average(row, yield_weight_columns)
        if not np.isfinite(rep_weight) or rep_weight <= 0.0:
            rep_weight = 1.0
        # endif

        phi_rad = math.radians(phi)
        xB_sum += rep_weight * xB
        Q2_sum += rep_weight * Q2
        t_sum += rep_weight * t_abs
        sin_phi_sum += rep_weight * math.sin(phi_rad)
        cos_phi_sum += rep_weight * math.cos(phi_rad)
        rep_weight_sum += rep_weight
        integration_weight_sum += integration_weight
        n_rows += 1
    # endfor

    if rep_weight_sum <= 0.0 or n_rows == 0:
        return RepresentativeKinematics(math.nan, math.nan, math.nan, math.nan, math.nan, 0)
    # endif

    circular_phi = math.degrees(math.atan2(sin_phi_sum, cos_phi_sum)) % 360.0
    return RepresentativeKinematics(
        xB=xB_sum / rep_weight_sum,
        Q2=Q2_sum / rep_weight_sum,
        t_abs=t_sum / rep_weight_sum,
        phi=circular_phi,
        integration_weight_sum=integration_weight_sum,
        n_rows=n_rows,
    )


def build_km15_prediction_points_all_projections(
    df: pd.DataFrame,
    prediction_period: str,
    no_width_weighting: bool,
    prediction_ebeam: float,
    prediction_mode: str,
    max_workers: int,
    no_parallel_km15: bool,
    diagnostics: ModelDiagnostics,
) -> Dict[str, List[Point]]:
    if prediction_mode not in {"representative", "row-sum"}:
        raise ValueError("Invalid KM15 prediction mode.")
    # endif

    if prediction_mode == "representative":
        print(
            "\nWARNING: KM15 representative mode is diagnostic only. A nonlinear "
            "model evaluated once at averaged hidden kinematics is not a true projection."
        )
    # endif

    tasks: List[Tuple[int, float, float, float, float, float]] = []
    metadata: Dict[int, Dict[str, object]] = {}
    task_id = 0

    for projection in KM15_PREDICTION_PROJECTION_ORDER:
        info = PROJECTIONS[projection]
        group_cols = [str(info["min_col"]), str(info["max_col"])]
        yield_columns = candidate_yield_columns_for_period(df, prediction_period)

        for _, group in df.sort_values(group_cols).groupby(group_cols, sort=True, dropna=True):
            key = (
                float(group[str(info["min_col"])].iloc[0]),
                float(group[str(info["max_col"])].iloc[0]),
            )
            x_plot = average_x_for_group(group, info, prediction_period, yield_columns)

            if prediction_mode == "representative":
                rep = representative_kinematics_for_group(
                    group, projection, prediction_period, no_width_weighting, yield_columns
                )
                if rep.n_rows <= 0 or not all(
                    np.isfinite(v)
                    for v in [rep.xB, rep.Q2, rep.t_abs, rep.phi, rep.integration_weight_sum]
                ):
                    continue
                # endif

                tasks.append((task_id, rep.xB, rep.Q2, rep.t_abs, rep.phi, prediction_ebeam))
                metadata[task_id] = {
                    "projection": projection,
                    "key": key,
                    "x": x_plot,
                    "weight": rep.integration_weight_sum,
                }
                task_id += 1
            else:
                xs_col = cross_section_column(prediction_period)
                for _, row in group.iterrows():
                    observed_value, _, _ = parse_tuple3(row.get(xs_col, math.nan))
                    if not np.isfinite(observed_value):
                        continue
                    # endif

                    xB = row_average_or_midpoint(row, f"xBavg, {prediction_period}", "xBmin", "xBmax")
                    Q2 = row_average_or_midpoint(row, f"Q2avg, {prediction_period}", "Q2min", "Q2max")
                    t_abs = row_average_or_midpoint(row, f"t_abs_avg, {prediction_period}", "t_abs_min", "t_abs_max")
                    phi = row_average_or_midpoint(row, f"phiavg, {prediction_period}", "phimin", "phimax")
                    if not all(np.isfinite(v) for v in [xB, Q2, t_abs, phi]):
                        continue
                    # endif

                    weight = weight_for_row(row, list(info["integrate_widths"]), no_width_weighting)
                    if not np.isfinite(weight):
                        continue
                    # endif

                    tasks.append((task_id, xB, Q2, t_abs, phi, prediction_ebeam))
                    metadata[task_id] = {
                        "projection": projection,
                        "key": key,
                        "x": x_plot,
                        "weight": weight,
                    }
                    task_id += 1
                # endfor
            # endif
        # endfor
    # endfor

    print("\nPrecomputing KM15 kinematic-dependence prediction points:")
    print(f"  prediction mode       = {prediction_mode}")
    print(f"  prediction period     = {prediction_period}")
    print(f"  prediction Ebeam      = {prediction_ebeam:.6f} GeV")
    print(f"  KM15 tasks            = {len(tasks)}")
    print(f"  parallel              = {not no_parallel_km15 and max_workers > 1}")
    print(f"  max workers           = {max_workers}")

    results = run_km15_prediction_tasks(tasks, max_workers, no_parallel_km15)
    diagnostics.km15_prediction_count += len(tasks)

    grouped: Dict[Tuple[str, Tuple[float, float]], Dict[str, float]] = {}
    for tid, meta in metadata.items():
        value, fatal_error = results.get(tid, (math.nan, "missing result"))
        if fatal_error is not None:
            diagnostics.km15_prediction_failure_count += 1
            if diagnostics.km15_prediction_failure_message is None:
                diagnostics.km15_prediction_failure_message = fatal_error
            # endif
            continue
        # endif
        if not np.isfinite(value):
            continue
        # endif

        group_key = (str(meta["projection"]), meta["key"])
        if group_key not in grouped:
            grouped[group_key] = {"x": float(meta["x"]), "y": 0.0, "n": 0.0}
        # endif
        grouped[group_key]["y"] += value * float(meta["weight"])
        grouped[group_key]["n"] += 1.0
    # endfor

    points_by_projection: Dict[str, List[Point]] = {
        projection: [] for projection in KM15_PREDICTION_PROJECTION_ORDER
    }
    for (projection, key), values in grouped.items():
        points_by_projection[projection].append(
            Point(
                key=key,
                x=values["x"],
                y=values["y"],
                stat=0.0,
                p2p=0.0,
                scale=0.0,
                n_rows=int(values["n"]),
                volume=math.nan,
            )
        )
    # endfor

    for projection in KM15_PREDICTION_PROJECTION_ORDER:
        points_by_projection[projection].sort(key=lambda p: p.x)
    # endfor

    return points_by_projection


# -----------------------------------------------------------------------------
# Plotting.
# -----------------------------------------------------------------------------

def point_total_error(point: Point, include_systematics: bool) -> float:
    stat = point.stat if np.isfinite(point.stat) else 0.0
    if not include_systematics:
        return stat
    # endif
    p2p = point.p2p if np.isfinite(point.p2p) else 0.0
    scale = point.scale if np.isfinite(point.scale) else 0.0
    return math.sqrt(stat * stat + p2p * p2p + scale * scale)


def points_to_arrays(points: Sequence[Point], include_systematics: bool) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray([p.x for p in points], dtype=float)
    y = np.asarray([p.y for p in points], dtype=float)
    stat = np.asarray([p.stat if np.isfinite(p.stat) else 0.0 for p in points], dtype=float)
    total = np.asarray([point_total_error(p, include_systematics) for p in points], dtype=float)
    return x, y, stat, total


def plot_points(ax, points: Sequence[Point], label: str, ratio: bool, include_systematics: bool) -> None:
    if not points:
        return
    # endif

    x, y, stat_err, total_err = points_to_arrays(points, include_systematics)
    style = SERIES_STYLE.get(label, dict(marker="o", linestyle="-"))
    color = style.get("color")

    if include_systematics:
        ax.errorbar(
            x,
            y,
            yerr=total_err,
            label=None,
            linewidth=0.0,
            elinewidth=2.2,
            capsize=4.0,
            alpha=0.28,
            color=color,
            zorder=1,
        )
    # endif

    ax.errorbar(
        x,
        y,
        yerr=stat_err,
        label=label,
        markersize=5.5,
        linewidth=1.2,
        elinewidth=1.0,
        capsize=2.5,
        zorder=3,
        **style,
    )

    if ratio:
        ax.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--", zorder=0)
    # endif


def set_ratio_ylim(ax) -> None:
    ax.set_ylim(RATIO_YMIN, RATIO_YMAX)


def make_projection_plot(
    df: pd.DataFrame,
    projection: str,
    output_dir: str,
    no_width_weighting: bool,
    include_systematics: bool,
    additional_p2p_frac: float,
    apply_sp19_to_fa18_scaling: bool,
) -> bool:
    info = PROJECTIONS[projection]

    points_by_period = {
        period: integrated_points_for_period(
            df=df,
            projection=projection,
            period=period,
            no_width_weighting=no_width_weighting,
            additional_p2p_frac=additional_p2p_frac,
        )
        for period in LEFT_SERIES
    }

    # Middle ratios: each pair is integrated on its own exact common support.
    middle_ratios: Dict[str, List[Point]] = {}
    for period in MIDDLE_SERIES:
        mask = common_support_mask(df, [period, "10.6 GeV"])
        numerator = integrated_points_for_period(
            df, projection, period, no_width_weighting,
            additional_p2p_frac=additional_p2p_frac,
            support_mask=mask,
        )
        denominator = integrated_points_for_period(
            df, projection, "10.6 GeV", no_width_weighting,
            additional_p2p_frac=additional_p2p_frac,
            support_mask=mask,
        )
        middle_ratios[period] = ratio_points_independent_approx(numerator, denominator)
    # endfor

    # Right comparison: both measurements use the same finite rows and finite
    # KM15 scaling factors.
    right_mask = common_support_mask(
        df,
        ["Fa18 Inb", "Sp19 Inb"],
        require_sp19_scale=apply_sp19_to_fa18_scaling,
    )
    fa_right = integrated_points_for_period(
        df, projection, "Fa18 Inb", no_width_weighting,
        additional_p2p_frac=additional_p2p_frac,
        support_mask=right_mask,
    )
    sp_right = integrated_points_for_period(
        df, projection, "Sp19 Inb", no_width_weighting,
        additional_p2p_frac=additional_p2p_frac,
        apply_sp19_to_fa18_scaling=apply_sp19_to_fa18_scaling,
        support_mask=right_mask,
    )

    fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.5), constrained_layout=True)
    fig.suptitle(str(info["title"]), fontsize=16)

    left = axes[0]
    for period in LEFT_SERIES:
        plot_points(left, points_by_period.get(period, []), period, False, include_systematics)
    # endfor
    left.set_xlabel(str(info["xlabel"]))
    left.set_ylabel(str(info["ylabel"]))
    left.set_title("Integrated cross sections")
    left.grid(True, alpha=0.25)
    left.legend(fontsize=8, frameon=True)

    middle = axes[1]
    for period in MIDDLE_SERIES:
        plot_points(middle, middle_ratios.get(period, []), period, True, include_systematics)
    # endfor
    middle.set_xlabel(str(info["xlabel"]))
    middle.set_ylabel("run period / 10.6 GeV")
    middle.set_title("10.6-GeV period consistency\n(common support; ratio stat. errors approximate)")
    middle.grid(True, alpha=0.25)
    middle.legend(fontsize=8, frameon=True)
    set_ratio_ylim(middle)

    right = axes[2]
    for label in RIGHT_SERIES_LABELS:
        points = ratio_to_own_two_point_weighted_mean(fa_right, sp_right, label)
        plot_points(right, points, label, True, include_systematics)
    # endfor
    right.set_xlabel(str(info["xlabel"]))
    right.set_ylabel("period / weighted mean")
    right.set_title("Fa18 Inb vs Sp19 Inb scaled to 10.604 GeV\n(common support)")
    right.grid(True, alpha=0.25)
    right.legend(fontsize=7, frameon=True)
    set_ratio_ylim(right)

    outbase = os.path.join(output_dir, f"integrated_{info['outfile_tag']}_dependence")
    fig.savefig(outbase + ".png", dpi=200)
    plt.close(fig)
    return True


def make_km15_prediction_plot(
    prediction_points_by_projection: Dict[str, List[Point]],
    output_dir: str,
    output_name: str,
    prediction_period: str,
    prediction_ebeam: float,
    prediction_mode: str,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 10.0), constrained_layout=True)
    fig.suptitle(
        "KM15 projected cross-section predictions\n"
        rf"CSV kinematics from {prediction_period}; $E_{{beam}}={prediction_ebeam:.4f}$ GeV; "
        rf"mode={prediction_mode}",
        fontsize=15,
    )

    for ax, projection in zip(axes.ravel(), KM15_PREDICTION_PROJECTION_ORDER):
        info = PROJECTIONS[projection]
        plot_points(
            ax,
            prediction_points_by_projection.get(projection, []),
            KM15_PREDICTION_LABEL,
            False,
            False,
        )
        ax.set_xlabel(str(info["xlabel"]))
        ax.set_ylabel(str(info["ylabel"]))
        ax.set_title(str(info["title"]))
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8, frameon=True)
    # endfor

    output_path = os.path.join(output_dir, output_name)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Wrote KM15 prediction PNG plot to: {output_path}")


# -----------------------------------------------------------------------------
# Diagnostics summaries.
# -----------------------------------------------------------------------------

def print_single_correction_summary(title: str, count: int, correction_sum: float, correction_min: float, correction_max: float) -> None:
    if count <= 0:
        print(f"{title}: no corrections were evaluated")
        return
    # endif
    print(f"{title}:")
    print(f"  correction count   = {count}")
    print(f"  correction mean    = {correction_sum / count:.6g}")
    print(f"  correction range   = [{correction_min:.6g}, {correction_max:.6g}]")


def print_model_correction_summary(diagnostics: ModelDiagnostics) -> None:
    if not diagnostics.enabled:
        print("Sp19 -> 10.604 GeV beam-energy scaling: disabled")
        return
    # endif

    print("\nSp19 -> 10.604 GeV beam-energy scaling:")
    print(f"  source beam energy = {SP19_INB_EBEAM_GEV:.6f} GeV")
    print(f"  target beam energy = {FA18_INB_EBEAM_GEV:.6f} GeV")
    print("  applied correction = KM15")
    print("  BH correction      = diagnostic only")
    print(f"  unique correction points = {diagnostics.sp19_unique_correction_count}")
    print()

    print_single_correction_summary(
        "KM15 ratio, KM15(10.604 GeV) / KM15(10.1998 GeV)",
        diagnostics.km15_correction_count,
        diagnostics.km15_correction_sum,
        diagnostics.km15_correction_min,
        diagnostics.km15_correction_max,
    )
    print()
    print_single_correction_summary(
        "BH ratio, BH(10.604 GeV) / BH(10.1998 GeV)",
        diagnostics.bh_correction_count,
        diagnostics.bh_correction_sum,
        diagnostics.bh_correction_min,
        diagnostics.bh_correction_max,
    )

    if diagnostics.bh_correction_count <= 0 and diagnostics.bh_failure_message is not None:
        print("\nBH diagnostic failure reason:")
        print(f"  {diagnostics.bh_failure_message}")
    # endif

    if diagnostics.km15_prediction_count > 0:
        print("\nKM15 prediction evaluation:")
        print(f"  prediction tasks attempted = {diagnostics.km15_prediction_count}")
        print(f"  prediction task failures   = {diagnostics.km15_prediction_failure_count}")
        if diagnostics.km15_prediction_failure_message is not None:
            print("  first prediction failure:")
            print(diagnostics.km15_prediction_failure_message)
        # endif
    # endif


# -----------------------------------------------------------------------------
# CLI.
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make integrated DVCS cross-section projections with support, closure and uncertainty diagnostics."
    )
    parser.add_argument("csv_file", help="Input final DVCS pass-2 analysis CSV.")
    parser.add_argument(
        "--output-dir",
        default="output/integrated_study",
        help="Directory for output plots.",
    )
    parser.add_argument(
        "--no-width-weighting",
        action="store_true",
        help="Use raw sums instead of differential-bin-width integration.",
    )
    parser.add_argument(
        "--phi-degrees",
        action="store_true",
        help="Use phi widths in degrees rather than radians.",
    )
    parser.add_argument("--theta-bins", type=int, default=DEFAULT_THETA_BINS)
    parser.add_argument("--theta-binning-period", default=DEFAULT_THETA_BINNING_PERIOD)
    parser.add_argument(
        "--include-systematics",
        action="store_true",
        help="Draw outer total bars containing stat, CSV point-to-point, correlated scale and optional extra p2p terms.",
    )
    parser.add_argument(
        "--include-bin-to-bin-sys",
        action="store_true",
        help="Deprecated alias for --include-systematics.",
    )
    parser.add_argument(
        "--bin-to-bin-sys-frac",
        type=float,
        default=0.0,
        help=(
            "Additional independent row-level fractional point-to-point uncertainty. "
            "It is propagated in quadrature during integration. Default: 0.0."
        ),
    )
    parser.add_argument("--no-sp19-km15-energy-scaling", action="store_true")
    parser.add_argument("--compute-bh-diagnostic", action="store_true")
    parser.add_argument("--max-workers", type=int, default=DEFAULT_MAX_WORKERS)
    parser.add_argument("--no-parallel-km15", action="store_true")
    parser.add_argument("--no-km15-prediction-plot", action="store_true")
    parser.add_argument(
        "--km15-prediction-output-name",
        default=DEFAULT_KM15_PREDICTION_OUTPUT_NAME,
    )
    parser.add_argument(
        "--km15-prediction-period",
        default=DEFAULT_KM15_PREDICTION_PERIOD,
    )
    parser.add_argument(
        "--km15-prediction-ebeam",
        type=float,
        default=DEFAULT_KM15_PREDICTION_EBEAM_GEV,
    )
    parser.add_argument(
        "--km15-prediction-mode",
        choices=["representative", "row-sum"],
        default=DEFAULT_KM15_PREDICTION_MODE,
    )
    parser.add_argument(
        "--closure-tolerance",
        type=float,
        default=DEFAULT_CLOSURE_TOLERANCE,
        help="Maximum allowed relative projection-closure difference.",
    )
    parser.add_argument(
        "--strict-closure",
        action="store_true",
        help="Exit with an error if any closure exceeds --closure-tolerance.",
    )
    parser.add_argument(
        "--no-closure-diagnostics",
        action="store_true",
        help="Disable projection-closure diagnostics.",
    )
    parser.add_argument(
        "--no-theta-bin-diagnostics",
        action="store_true",
        help="Disable detector-angle occupancy and volume diagnostics.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.bin_to_bin_sys_frac < 0.0:
        raise ValueError("--bin-to-bin-sys-frac must be non-negative.")
    # endif
    if args.theta_bins <= 0:
        raise ValueError("--theta-bins must be positive.")
    # endif
    if args.km15_prediction_ebeam <= 0.0:
        raise ValueError("--km15-prediction-ebeam must be positive.")
    # endif
    if args.max_workers <= 0:
        raise ValueError("--max-workers must be positive.")
    # endif
    if args.closure_tolerance < 0.0:
        raise ValueError("--closure-tolerance must be non-negative.")
    # endif

    if args.max_workers > DEFAULT_MAX_WORKERS:
        print(f"Requested --max-workers {args.max_workers}; capping to {DEFAULT_MAX_WORKERS}.")
        args.max_workers = DEFAULT_MAX_WORKERS
    # endif

    include_systematics = bool(args.include_systematics or args.include_bin_to_bin_sys)

    required = [
        "xBmin", "xBmax", "Q2min", "Q2max",
        "t_abs_min", "t_abs_max", "phimin", "phimax",
    ]
    required += [
        theta_average_column("e_theta", args.theta_binning_period),
        theta_average_column("p_theta", args.theta_binning_period),
        theta_average_column("g_theta", args.theta_binning_period),
    ]
    required += [cross_section_column(period) for period in LEFT_SERIES]

    df = pd.read_csv(args.csv_file)
    require_columns(df, required)
    df = add_width_columns(df, phi_degrees=args.phi_degrees)
    df = add_derived_theta_bin_columns(df, args.theta_binning_period, args.theta_bins)
    os.makedirs(args.output_dir, exist_ok=True)

    df, diagnostics = precompute_sp19_to_fa18_km15_scales(
        df=df,
        enabled=not args.no_sp19_km15_energy_scaling,
        compute_bh_diagnostic=args.compute_bh_diagnostic,
        max_workers=args.max_workers,
        no_parallel_km15=args.no_parallel_km15,
    )

    if not args.no_closure_diagnostics:
        print_projection_closure_diagnostics(
            df=df,
            periods=LEFT_SERIES,
            no_width_weighting=args.no_width_weighting,
            phi_degrees=args.phi_degrees,
            additional_p2p_frac=args.bin_to_bin_sys_frac,
            tolerance=args.closure_tolerance,
            strict=args.strict_closure,
        )
    # endif

    if not args.no_theta_bin_diagnostics:
        print_theta_bin_diagnostics(
            df=df,
            periods=LEFT_SERIES,
            no_width_weighting=args.no_width_weighting,
        )
    # endif

    made_projection_count = 0
    for projection in DATA_PROJECTION_ORDER:
        if make_projection_plot(
            df=df,
            projection=projection,
            output_dir=args.output_dir,
            no_width_weighting=args.no_width_weighting,
            include_systematics=include_systematics,
            additional_p2p_frac=args.bin_to_bin_sys_frac,
            apply_sp19_to_fa18_scaling=not args.no_sp19_km15_energy_scaling,
        ):
            made_projection_count += 1
        # endif
    # endfor

    if not args.no_km15_prediction_plot:
        prediction_points = build_km15_prediction_points_all_projections(
            df=df,
            prediction_period=args.km15_prediction_period,
            no_width_weighting=args.no_width_weighting,
            prediction_ebeam=args.km15_prediction_ebeam,
            prediction_mode=args.km15_prediction_mode,
            max_workers=args.max_workers,
            no_parallel_km15=args.no_parallel_km15,
            diagnostics=diagnostics,
        )
        make_km15_prediction_plot(
            prediction_points_by_projection=prediction_points,
            output_dir=args.output_dir,
            output_name=args.km15_prediction_output_name,
            prediction_period=args.km15_prediction_period,
            prediction_ebeam=args.km15_prediction_ebeam,
            prediction_mode=args.km15_prediction_mode,
        )
    # endif

    if include_systematics:
        print(
            "\nPlotted inner bars: statistical. Outer bars: "
            "sqrt(stat^2 + p2p^2 + correlated_scale^2)."
        )
        print(f"Additional independent row-level p2p fraction = {args.bin_to_bin_sys_frac:.6g}")
    else:
        print("\nPlotted y-errors: statistical only")
    # endif

    print(f"All ratio y axes fixed to [{RATIO_YMIN:.2f}, {RATIO_YMAX:.2f}]")
    print(f"Derived theta projections used {args.theta_bins} bins from '{args.theta_binning_period}'.")
    print(f"KM15 prediction mode = {args.km15_prediction_mode}")
    print_model_correction_summary(diagnostics)
    print(f"\nWrote {made_projection_count} integrated dependence PNG plots to: {args.output_dir}")


if __name__ == "__main__":
    main()
# endif
