#!/usr/bin/env python3
"""
integrated_cross_section_study.py

Standalone projection/integration study for the DVCS pass-2 CSV.

The script reads the final analysis CSV and makes 1x3 canvases in
output/integrated_study/:

  integrated_xB_dependence.png
  integrated_Q2_dependence.png
  integrated_t_dependence.png
  integrated_phi_dependence.png
  integrated_e_theta_dependence.png
  integrated_p_theta_dependence.png
  integrated_g_theta_dependence.png

The detector-angle plots are made only if the corresponding CSV columns exist:

  e_theta_min, e_theta_max, e_theta_avg, <period>
  p_theta_min, p_theta_max, p_theta_avg, <period>
  g_theta_min, g_theta_max, g_theta_avg, <period>

Each 1x3 canvas contains:

  left:   integrated normed cross sections for:
            10.6 GeV
            Sp18 Inb
            Sp18 Out
            Fa18 Inb
            Fa18 Out
            Sp19 Inb

  middle: ratios of the four individual 10.6-GeV run periods to the combined
          10.6-GeV result.

  right:  comparison of:
            Fa18 Inb, 10.604 GeV
            Sp19 Inb -> 10.604 GeV, KM15

          Each is divided by the weighted mean of those two values.

The script also makes one additional 2x2 canvas:

  km15_kinematic_dependences.png

showing the KM15 model prediction for the four physical kinematic dependences:
xB, Q2, |t|, and phi.

By default the KM15 prediction plot uses representative mode:

  --km15-prediction-mode representative

which evaluates KM15 once per projected bin using count/yield-weighted average
hidden kinematics from the CSV.

A slower but more faithful option is:

  --km15-prediction-mode row-sum

which evaluates KM15 on every underlying multidimensional CSV row and sums the
model with the same bin-width integration weights used for the data projection.

By default the integration is bin-width weighted, appropriate for differential
cross sections:

  d sigma / dxB       = sum_i sigma_i * dQ2_i * dt_i * dphi_i
  d sigma / dQ2       = sum_i sigma_i * dxB_i * dt_i * dphi_i
  d sigma / d|t|      = sum_i sigma_i * dxB_i * dQ2_i * dphi_i
  d sigma / dphi      = sum_i sigma_i * dxB_i * dQ2_i * dt_i

For detector-angle projections, the script integrates over the full physics
volume:

  sigma(theta bin) = sum_i sigma_i * dxB_i * dQ2_i * dt_i * dphi_i

so those y axes are labeled as integrated cross sections rather than
d sigma / d theta.

The phi width is converted from degrees to radians by default, since the usual
DVCS differential cross section convention is per radian. Use --phi-degrees if
you intentionally want degree-weighted phi integration instead.

For the plotted x-position in each projected bin, the script tries to use the
available per-bin average kinematic columns, for example:

  xBavg, Fa18 Inb
  Q2avg, Fa18 Inb
  t_abs_avg, Fa18 Inb
  phiavg, Fa18 Inb
  e_theta_avg, Fa18 Inb
  p_theta_avg, Fa18 Inb
  g_theta_avg, Fa18 Inb

and weights those sub-bin average values by event/yield columns if available.
If no suitable yield/count column is found, it falls back to a finite-value
average of the sub-bin average kinematics. If that also fails, it falls back
to the bin midpoint.

By default, plotted y-error bars are statistical only.

If --include-bin-to-bin-sys is passed, each point shows two error bars:

  inner darker error bar:
    stat_error

  outer lighter error bar:
    total_error = sqrt(stat_error^2 + (bin_to_bin_sys_fraction * cross_section)^2)

where bin_to_bin_sys_fraction defaults to 0.10 and can be changed with:

  --bin-to-bin-sys-frac 0.10

Sp19 -> Fa18 beam-energy scaling
--------------------------------

For the right-panel Fa18/Sp19 comparison, Sp19 Inb is scaled from its beam energy
to the Fa18 beam energy before integration:

  E_source = 10.1998 GeV
  E_target = 10.6040 GeV

The row-by-row applied correction factor is:

  C_i = KM15(E_target, xB_i, Q2_i, |t|_i, phi_i)
      / KM15(E_source, xB_i, Q2_i, |t|_i, phi_i)

and the scaled Sp19 contribution is:

  sigma_i_scaled = C_i * sigma_i
  stat_i_scaled  = |C_i| * stat_i
  sys_i_scaled   = |C_i| * sys_i

This is enabled by default. Disable with:

  --no-sp19-km15-energy-scaling

The correction factors are precomputed once before plotting and stored in a
temporary dataframe column:

  _sp19_to_fa18_km15_scale

This avoids re-evaluating Gepard/KM15 inside every projected plot.

The analogous BH-only ratio can also be computed and printed as a diagnostic:

  C_i^BH = BH(E_target, xB_i, Q2_i, |t|_i, phi_i)
         / BH(E_source, xB_i, Q2_i, |t|_i, phi_i)

but the BH ratio is not applied to the data. It is disabled by default because
it is diagnostic only. Enable it with:

  --compute-bh-diagnostic

Parallelization
---------------

KM15 model evaluations are parallelized by default using process workers:

  --max-workers 5

Disable process parallelization with:

  --no-parallel-km15

All ratio-plot y axes are fixed to RATIO_YMIN--RATIO_YMAX.
"""

import argparse
import concurrent.futures
import math
import os
import re
import traceback
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Beam energies.
# -----------------------------------------------------------------------------

FA18_INB_EBEAM_GEV = 10.6040
SP19_INB_EBEAM_GEV = 10.1998


# -----------------------------------------------------------------------------
# Ratio-axis settings.
# -----------------------------------------------------------------------------

RATIO_YMIN = 0.0
RATIO_YMAX = 1.6


# -----------------------------------------------------------------------------
# Parallelism settings.
# -----------------------------------------------------------------------------

DEFAULT_MAX_WORKERS = 5


# -----------------------------------------------------------------------------
# Temporary dataframe columns.
# -----------------------------------------------------------------------------

SP19_TO_FA18_KM15_SCALE_COLUMN = "_sp19_to_fa18_km15_scale"
SP19_TO_FA18_BH_SCALE_COLUMN = "_sp19_to_fa18_bh_scale"


# -----------------------------------------------------------------------------
# KM15 prediction-plot settings.
# -----------------------------------------------------------------------------

DEFAULT_KM15_PREDICTION_OUTPUT_NAME = "km15_kinematic_dependences.png"
DEFAULT_KM15_PREDICTION_PERIOD = "10.6 GeV"
DEFAULT_KM15_PREDICTION_EBEAM_GEV = FA18_INB_EBEAM_GEV
DEFAULT_KM15_PREDICTION_MODE = "representative"
KM15_PREDICTION_LABEL = "KM15, 10.604 GeV"


# -----------------------------------------------------------------------------
# Plot order.
# -----------------------------------------------------------------------------

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

RIGHT_RAW_PERIODS = [
    "Fa18 Inb",
    "Sp19 Inb",
]

RIGHT_SERIES_LABELS = [
    "Fa18 Inb, 10.604 GeV",
    "Sp19 Inb -> 10.604 GeV, KM15",
]


# -----------------------------------------------------------------------------
# Marker and line styles.
# -----------------------------------------------------------------------------

SERIES_STYLE = {
    "10.6 GeV": dict(marker="o", linestyle="-", color="black"),
    "Sp18 Inb": dict(marker="D", linestyle="-", color="tab:green"),
    "Sp18 Out": dict(marker="P", linestyle="-", color="tab:purple"),
    "Fa18 Inb": dict(marker="^", linestyle="-", color="tab:blue"),
    "Fa18 Out": dict(marker="v", linestyle="-", color="tab:orange"),
    "Sp19 Inb": dict(marker="s", linestyle="-", color="tab:red"),
    "Fa18 Inb, 10.604 GeV": dict(marker="^", linestyle="-", color="tab:blue"),
    "Sp19 Inb -> 10.604 GeV, KM15": dict(marker="s", linestyle="-", color="tab:red"),
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
        "required": True,
        "km15_prediction": True,
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
        "required": True,
        "km15_prediction": True,
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
        "required": True,
        "km15_prediction": True,
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
        "required": True,
        "km15_prediction": True,
    },
    "e_theta": {
        "min_col": "e_theta_min",
        "max_col": "e_theta_max",
        "avg_prefix": "e_theta_avg",
        "xlabel": r"$\theta_{e}$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"electron polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "e_theta",
        "required": False,
        "km15_prediction": False,
    },
    "p_theta": {
        "min_col": "p_theta_min",
        "max_col": "p_theta_max",
        "avg_prefix": "p_theta_avg",
        "xlabel": r"$\theta_{p}$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"proton polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "p_theta",
        "required": False,
        "km15_prediction": False,
    },
    "g_theta": {
        "min_col": "g_theta_min",
        "max_col": "g_theta_max",
        "avg_prefix": "g_theta_avg",
        "xlabel": r"$\theta_{\gamma}$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"photon polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "g_theta",
        "required": False,
        "km15_prediction": False,
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

KM15_PREDICTION_PROJECTION_ORDER = [
    "xB",
    "Q2",
    "t",
    "phi",
]


FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)


# -----------------------------------------------------------------------------
# Dataclasses.
# -----------------------------------------------------------------------------

@dataclass
class Point:
    key: Tuple[float, float]
    x: float
    y: float
    stat: float
    sys: float


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

    correction_cache: Dict[
        Tuple[float, float, float, float],
        Tuple[float, float],
    ] = field(default_factory=dict)

    km15_prediction_cache: Dict[
        Tuple[float, float, float, float, float],
        float,
    ] = field(default_factory=dict)


# -----------------------------------------------------------------------------
# Worker globals for process pools.
# -----------------------------------------------------------------------------

_WORKER_GEPARD = None
_WORKER_TH_KM15 = None


def worker_get_model_context() -> ModelContext:
    """
    Lazily import Gepard/KM15 once per worker process.
    """

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

    return ModelContext(
        enabled=True,
        g=_WORKER_GEPARD,
        th_km15=_WORKER_TH_KM15,
    )


# -----------------------------------------------------------------------------
# Gepard / KM15 / BH helpers.
# -----------------------------------------------------------------------------

def import_gepard():
    try:
        import gepard as g
    except ImportError as exc:
        raise RuntimeError(
            "Could not import the Gepard package.\n\n"
            "Install it into the Python interpreter used to run this script.\n"
            "On ifarm, prefer:\n"
            "  python -m pip install --user gepard\n"
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

    return ModelContext(
        enabled=True,
        g=import_gepard(),
        th_km15=import_km15(),
    )


def make_gepard_xs_point(
    g,
    xB: float,
    Q2: float,
    t_abs: float,
    phi_deg_trento: float,
    ebeam: float,
):
    """
    Build a Gepard DataPoint at fixed-target kinematics.
    """

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


def km15_xs(
    model_context: ModelContext,
    xB: float,
    Q2: float,
    t_abs: float,
    phi_deg_trento: float,
    ebeam: float,
) -> float:
    pt = make_gepard_xs_point(
        g=model_context.g,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi_deg_trento,
        ebeam=ebeam,
    )

    return float(model_context.th_km15.predict(pt))


def km15_xs_cached(
    model_context: ModelContext,
    xB: float,
    Q2: float,
    t_abs: float,
    phi_deg_trento: float,
    ebeam: float,
) -> float:
    key = (
        round(float(xB), 8),
        round(float(Q2), 8),
        round(float(t_abs), 8),
        round(float(phi_deg_trento), 8),
        round(float(ebeam), 8),
    )

    if key in model_context.km15_prediction_cache:
        return model_context.km15_prediction_cache[key]
    # endif

    value = km15_xs(
        model_context=model_context,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi_deg_trento,
        ebeam=ebeam,
    )

    model_context.km15_prediction_cache[key] = value

    return value


def bh_xs(
    model_context: ModelContext,
    xB: float,
    Q2: float,
    t_abs: float,
    phi_deg_trento: float,
    ebeam: float,
) -> float:
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
        "The loaded KM15 theory object does not expose PreFacSigma(pt) and TBH2unp(pt). "
        f"Candidate BH/prefactor-like attributes found: {available}"
    )


def correction_cache_key(
    xB: float,
    Q2: float,
    t_abs: float,
    phi: float,
) -> Tuple[float, float, float, float]:
    return (
        round(float(xB), 8),
        round(float(Q2), 8),
        round(float(t_abs), 8),
        round(float(phi), 8),
    )


def update_model_correction_diagnostics(
    diagnostics: ModelDiagnostics,
    correction: float,
    model_name: str,
) -> None:
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
    km15_numerator = km15_xs(
        model_context=model_context,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi,
        ebeam=FA18_INB_EBEAM_GEV,
    )

    km15_denominator = km15_xs(
        model_context=model_context,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi,
        ebeam=SP19_INB_EBEAM_GEV,
    )

    if (
        not np.isfinite(km15_numerator)
        or not np.isfinite(km15_denominator)
        or km15_denominator == 0.0
    ):
        return math.nan, math.nan, None
    # endif

    km15_corr = km15_numerator / km15_denominator

    if not compute_bh_diagnostic:
        return km15_corr, math.nan, None
    # endif

    bh_corr = math.nan
    bh_failure_message = None

    try:
        bh_numerator = bh_xs(
            model_context=model_context,
            xB=xB,
            Q2=Q2,
            t_abs=t_abs,
            phi_deg_trento=phi,
            ebeam=FA18_INB_EBEAM_GEV,
        )

        bh_denominator = bh_xs(
            model_context=model_context,
            xB=xB,
            Q2=Q2,
            t_abs=t_abs,
            phi_deg_trento=phi,
            ebeam=SP19_INB_EBEAM_GEV,
        )

        if (
            np.isfinite(bh_numerator)
            and np.isfinite(bh_denominator)
            and bh_denominator != 0.0
        ):
            bh_corr = bh_numerator / bh_denominator
        # endif

    except Exception as exc:
        bh_corr = math.nan
        bh_failure_message = repr(exc)
    # endtry

    return km15_corr, bh_corr, bh_failure_message


def worker_sp19_correction(task: Tuple[int, float, float, float, float, bool]) -> Tuple[int, float, float, Optional[str], Optional[str]]:
    task_id, xB, Q2, t_abs, phi, compute_bh_diagnostic = task

    try:
        model_context = worker_get_model_context()

        km15_corr, bh_corr, bh_failure_message = km15_and_bh_sp19_to_fa18_corrections_direct(
            xB=xB,
            Q2=Q2,
            t_abs=t_abs,
            phi=phi,
            compute_bh_diagnostic=compute_bh_diagnostic,
            model_context=model_context,
        )

        return task_id, km15_corr, bh_corr, bh_failure_message, None

    except Exception:
        return task_id, math.nan, math.nan, None, traceback.format_exc()
    # endtry


def worker_km15_prediction(task: Tuple[int, float, float, float, float, float]) -> Tuple[int, float, Optional[str]]:
    task_id, xB, Q2, t_abs, phi, ebeam = task

    try:
        model_context = worker_get_model_context()

        value = km15_xs(
            model_context=model_context,
            xB=xB,
            Q2=Q2,
            t_abs=t_abs,
            phi_deg_trento=phi,
            ebeam=ebeam,
        )

        return task_id, value, None

    except Exception:
        return task_id, math.nan, traceback.format_exc()
    # endtry


# -----------------------------------------------------------------------------
# Parsing and CSV helpers.
# -----------------------------------------------------------------------------

def parse_tuple3(value) -> Tuple[float, float, float]:
    if value is None:
        return (math.nan, math.nan, math.nan)
    # endif

    if isinstance(value, float) and math.isnan(value):
        return (math.nan, math.nan, math.nan)
    # endif

    if isinstance(value, (int, float, np.integer, np.floating)):
        return (float(value), 0.0, 0.0)
    # endif

    text = str(value).strip()

    if text == "" or text.lower() in {"nan", "none", "null"}:
        return (math.nan, math.nan, math.nan)
    # endif

    numbers = FLOAT_PATTERN.findall(text)

    if len(numbers) == 0:
        return (math.nan, math.nan, math.nan)
    # endif

    parsed = [float(x) for x in numbers]

    while len(parsed) < 3:
        parsed.append(0.0)
    # endwhile

    return (parsed[0], parsed[1], parsed[2])


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


def row_average_or_midpoint(
    row: pd.Series,
    avg_col: str,
    min_col: str,
    max_col: str,
) -> float:
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


def require_columns(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [c for c in columns if c not in df.columns]

    if missing:
        raise KeyError(
            "The input CSV is missing required columns:\n  " + "\n  ".join(missing)
        )
    # endif


def projection_columns_exist(df: pd.DataFrame, projection: str) -> bool:
    info = PROJECTIONS[projection]
    needed = [str(info["min_col"]), str(info["max_col"])]

    return all(col in df.columns for col in needed)


def add_width_columns(df: pd.DataFrame, phi_degrees: bool) -> pd.DataFrame:
    out = df.copy()

    out["_width_xB"] = (
        pd.to_numeric(out["xBmax"], errors="coerce")
        - pd.to_numeric(out["xBmin"], errors="coerce")
    )

    out["_width_Q2"] = (
        pd.to_numeric(out["Q2max"], errors="coerce")
        - pd.to_numeric(out["Q2min"], errors="coerce")
    )

    out["_width_t"] = (
        pd.to_numeric(out["t_abs_max"], errors="coerce")
        - pd.to_numeric(out["t_abs_min"], errors="coerce")
    )

    phi_width_deg = (
        pd.to_numeric(out["phimax"], errors="coerce")
        - pd.to_numeric(out["phimin"], errors="coerce")
    )

    if phi_degrees:
        out["_width_phi"] = phi_width_deg
    else:
        out["_width_phi"] = np.deg2rad(phi_width_deg)
    # endif

    return out


def weight_for_row(
    row: pd.Series,
    integrate_widths: List[str],
    no_width_weighting: bool,
) -> float:
    if no_width_weighting:
        return 1.0
    # endif

    weight = 1.0

    for axis in integrate_widths:
        width = row[f"_width_{axis}"]

        if not np.isfinite(width) or width <= 0.0:
            return math.nan
        # endif

        weight *= float(width)
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
        "numEvents",
    ]

    candidates: List[str] = []

    for phrase in priority_phrases:
        phrase_lower = phrase.lower()

        for col in df.columns:
            col_lower = str(col).lower()

            if period_lower not in col_lower:
                continue
            # endif

            if "unpol" not in col_lower:
                continue
            # endif

            if phrase_lower not in col_lower:
                continue
            # endif

            if col in candidates:
                continue
            # endif

            candidates.append(col)
        # endfor
    # endfor

    return candidates


def event_weight_for_average(row: pd.Series, candidate_columns: List[str]) -> float:
    for col in candidate_columns:
        val = parse_scalar_from_cell(row[col])

        if np.isfinite(val) and val > 0.0:
            return float(val)
        # endif
    # endfor

    return math.nan


def fallback_count_weight_for_model(row: pd.Series, period: str) -> float:
    _ = row
    _ = period

    return 1.0


def average_x_for_group(
    group: pd.DataFrame,
    projection_info: Dict[str, object],
    period: str,
    yield_weight_columns: List[str],
) -> float:
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
    finite_avg_values: List[float] = []

    for _, row in group.iterrows():
        xavg = parse_scalar_from_cell(row[avg_col])

        if not np.isfinite(xavg):
            continue
        # endif

        finite_avg_values.append(float(xavg))

        w = event_weight_for_average(row, yield_weight_columns)

        if not np.isfinite(w) or w <= 0.0:
            continue
        # endif

        numerator += w * xavg
        denominator += w
    # endfor

    if denominator > 0.0:
        return numerator / denominator
    # endif

    if len(finite_avg_values) > 0:
        return float(np.mean(np.array(finite_avg_values, dtype=float)))
    # endif

    return midpoint


# -----------------------------------------------------------------------------
# Parallel model precomputation.
# -----------------------------------------------------------------------------

def run_sp19_correction_tasks(
    tasks: List[Tuple[int, float, float, float, float, bool]],
    max_workers: int,
    no_parallel_km15: bool,
) -> Dict[int, Tuple[float, float, Optional[str], Optional[str]]]:
    results: Dict[int, Tuple[float, float, Optional[str], Optional[str]]] = {}

    if len(tasks) == 0:
        return results
    # endif

    if no_parallel_km15 or max_workers <= 1:
        model_context = make_model_context(enabled=True)

        for task_id, xB, Q2, t_abs, phi, compute_bh_diagnostic in tasks:
            try:
                km15_corr, bh_corr, bh_failure_message = km15_and_bh_sp19_to_fa18_corrections_direct(
                    xB=xB,
                    Q2=Q2,
                    t_abs=t_abs,
                    phi=phi,
                    compute_bh_diagnostic=compute_bh_diagnostic,
                    model_context=model_context,
                )

                results[task_id] = (km15_corr, bh_corr, bh_failure_message, None)

            except Exception:
                results[task_id] = (math.nan, math.nan, None, traceback.format_exc())
            # endtry
        # endfor

        return results
    # endif

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for task_id, km15_corr, bh_corr, bh_failure_message, fatal_error_message in executor.map(worker_sp19_correction, tasks):
            results[task_id] = (
                km15_corr,
                bh_corr,
                bh_failure_message,
                fatal_error_message,
            )
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

        xB = row_average_or_midpoint(
            row=row,
            avg_col="xBavg, Sp19 Inb",
            min_col="xBmin",
            max_col="xBmax",
        )

        Q2 = row_average_or_midpoint(
            row=row,
            avg_col="Q2avg, Sp19 Inb",
            min_col="Q2min",
            max_col="Q2max",
        )

        t_abs = row_average_or_midpoint(
            row=row,
            avg_col="t_abs_avg, Sp19 Inb",
            min_col="t_abs_min",
            max_col="t_abs_max",
        )

        phi = row_average_or_midpoint(
            row=row,
            avg_col="phiavg, Sp19 Inb",
            min_col="phimin",
            max_col="phimax",
        )

        if not (
            np.isfinite(xB)
            and np.isfinite(Q2)
            and np.isfinite(t_abs)
            and np.isfinite(phi)
        ):
            continue
        # endif

        key = correction_cache_key(
            xB=xB,
            Q2=Q2,
            t_abs=t_abs,
            phi=phi,
        )

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

    print()
    print("Precomputing Sp19 -> 10.604 GeV KM15 correction factors:")
    print(f"  unique kinematic points = {len(tasks)}")
    print(f"  parallel                = {not no_parallel_km15 and max_workers > 1}")
    print(f"  max workers             = {max_workers}")
    print(f"  BH diagnostic           = {compute_bh_diagnostic}")

    results_by_task_id = run_sp19_correction_tasks(
        tasks=tasks,
        max_workers=max_workers,
        no_parallel_km15=no_parallel_km15,
    )

    scale_by_key: Dict[Tuple[float, float, float, float], float] = {}
    bh_by_key: Dict[Tuple[float, float, float, float], float] = {}

    for key, task_id in task_id_by_key.items():
        km15_corr, bh_corr, bh_failure_message, fatal_error_message = results_by_task_id.get(
            task_id,
            (math.nan, math.nan, None, "missing result"),
        )

        if fatal_error_message is not None:
            if diagnostics.km15_prediction_failure_message is None:
                diagnostics.km15_prediction_failure_message = fatal_error_message
            # endif

        if bh_failure_message is not None and diagnostics.bh_failure_message is None:
            diagnostics.bh_failure_message = bh_failure_message
        # endif

        scale_by_key[key] = km15_corr
        bh_by_key[key] = bh_corr

        update_model_correction_diagnostics(
            diagnostics=diagnostics,
            correction=km15_corr,
            model_name="KM15",
        )

        update_model_correction_diagnostics(
            diagnostics=diagnostics,
            correction=bh_corr,
            model_name="BH",
        )
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

    if len(tasks) == 0:
        return results
    # endif

    if no_parallel_km15 or max_workers <= 1:
        model_context = make_model_context(enabled=True)

        for task_id, xB, Q2, t_abs, phi, ebeam in tasks:
            try:
                value = km15_xs_cached(
                    model_context=model_context,
                    xB=xB,
                    Q2=Q2,
                    t_abs=t_abs,
                    phi_deg_trento=phi,
                    ebeam=ebeam,
                )

                results[task_id] = (value, None)

            except Exception:
                results[task_id] = (math.nan, traceback.format_exc())
            # endtry
        # endfor

        return results
    # endif

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for task_id, value, fatal_error_message in executor.map(worker_km15_prediction, tasks):
            results[task_id] = (value, fatal_error_message)
        # endfor
    # endwith

    return results


# -----------------------------------------------------------------------------
# Integration.
# -----------------------------------------------------------------------------

def integrated_points_for_period(
    df: pd.DataFrame,
    projection: str,
    period: str,
    no_width_weighting: bool,
    apply_sp19_to_fa18_scaling: bool = False,
) -> List[Point]:
    if not projection_columns_exist(df, projection):
        return []
    # endif

    info = PROJECTIONS[projection]
    col = cross_section_column(period)

    if col not in df.columns:
        return []
    # endif

    if apply_sp19_to_fa18_scaling and SP19_TO_FA18_KM15_SCALE_COLUMN not in df.columns:
        raise RuntimeError(
            f"Sp19 -> Fa18 scaling requested, but {SP19_TO_FA18_KM15_SCALE_COLUMN} "
            "is not present. Run precompute_sp19_to_fa18_km15_scales first."
        )
    # endif

    yield_weight_columns = candidate_yield_columns_for_period(df, period)

    group_cols = [str(info["min_col"]), str(info["max_col"])]
    points: List[Point] = []

    sorted_df = df.sort_values(group_cols)

    for _, group in sorted_df.groupby(group_cols, sort=True, dropna=True):
        y_sum = 0.0
        stat2_sum = 0.0
        sys2_sum = 0.0
        n_used = 0

        for _, row in group.iterrows():
            value, stat, sys = parse_tuple3(row[col])

            if not np.isfinite(value):
                continue
            # endif

            weight = weight_for_row(
                row=row,
                integrate_widths=list(info["integrate_widths"]),
                no_width_weighting=no_width_weighting,
            )

            if not np.isfinite(weight):
                continue
            # endif

            if not np.isfinite(stat):
                stat = 0.0
            # endif

            if not np.isfinite(sys):
                sys = 0.0
            # endif

            model_scale = 1.0

            if apply_sp19_to_fa18_scaling:
                model_scale = parse_scalar_from_cell(row[SP19_TO_FA18_KM15_SCALE_COLUMN])

                if not np.isfinite(model_scale):
                    continue
                # endif
            # endif

            value *= model_scale
            stat *= abs(model_scale)
            sys *= abs(model_scale)

            y_sum += value * weight
            stat2_sum += (stat * weight) ** 2
            sys2_sum += (sys * weight) ** 2
            n_used += 1
        # endfor

        if n_used == 0:
            continue
        # endif

        x = average_x_for_group(
            group=group,
            projection_info=info,
            period=period,
            yield_weight_columns=yield_weight_columns,
        )

        key = (
            float(group[str(info["min_col"])].iloc[0]),
            float(group[str(info["max_col"])].iloc[0]),
        )

        points.append(
            Point(
                key=key,
                x=x,
                y=y_sum,
                stat=math.sqrt(stat2_sum),
                sys=math.sqrt(sys2_sum),
            )
        )
    # endfor

    points.sort(key=lambda p: p.x)

    return points


def representative_kinematics_for_group(
    group: pd.DataFrame,
    projection: str,
    period: str,
    no_width_weighting: bool,
    yield_weight_columns: List[str],
) -> RepresentativeKinematics:
    info = PROJECTIONS[projection]

    xB_sum = 0.0
    Q2_sum = 0.0
    t_sum = 0.0
    phi_sum = 0.0
    rep_weight_sum = 0.0
    integration_weight_sum = 0.0
    n_rows = 0

    xs_col = cross_section_column(period)

    for _, row in group.iterrows():
        if xs_col in row.index:
            observed_value, _, _ = parse_tuple3(row[xs_col])

            if not np.isfinite(observed_value):
                continue
            # endif
        # endif

        xB = row_average_or_midpoint(
            row=row,
            avg_col=f"xBavg, {period}",
            min_col="xBmin",
            max_col="xBmax",
        )

        Q2 = row_average_or_midpoint(
            row=row,
            avg_col=f"Q2avg, {period}",
            min_col="Q2min",
            max_col="Q2max",
        )

        t_abs = row_average_or_midpoint(
            row=row,
            avg_col=f"t_abs_avg, {period}",
            min_col="t_abs_min",
            max_col="t_abs_max",
        )

        phi = row_average_or_midpoint(
            row=row,
            avg_col=f"phiavg, {period}",
            min_col="phimin",
            max_col="phimax",
        )

        if not (
            np.isfinite(xB)
            and np.isfinite(Q2)
            and np.isfinite(t_abs)
            and np.isfinite(phi)
        ):
            continue
        # endif

        integration_weight = weight_for_row(
            row=row,
            integrate_widths=list(info["integrate_widths"]),
            no_width_weighting=no_width_weighting,
        )

        if not np.isfinite(integration_weight):
            continue
        # endif

        rep_weight = event_weight_for_average(row, yield_weight_columns)

        if not np.isfinite(rep_weight) or rep_weight <= 0.0:
            rep_weight = fallback_count_weight_for_model(row=row, period=period)
        # endif

        if not np.isfinite(rep_weight) or rep_weight <= 0.0:
            continue
        # endif

        xB_sum += rep_weight * xB
        Q2_sum += rep_weight * Q2
        t_sum += rep_weight * t_abs
        phi_sum += rep_weight * phi
        rep_weight_sum += rep_weight
        integration_weight_sum += integration_weight
        n_rows += 1
    # endfor

    if rep_weight_sum <= 0.0 or n_rows == 0:
        return RepresentativeKinematics(
            xB=math.nan,
            Q2=math.nan,
            t_abs=math.nan,
            phi=math.nan,
            integration_weight_sum=math.nan,
            n_rows=0,
        )
    # endif

    return RepresentativeKinematics(
        xB=xB_sum / rep_weight_sum,
        Q2=Q2_sum / rep_weight_sum,
        t_abs=t_sum / rep_weight_sum,
        phi=phi_sum / rep_weight_sum,
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
        raise ValueError(
            "--km15-prediction-mode must be either 'representative' or 'row-sum'."
        )
    # endif

    projection_order = KM15_PREDICTION_PROJECTION_ORDER

    task_id = 0
    tasks: List[Tuple[int, float, float, float, float, float]] = []
    task_metadata: Dict[int, Dict[str, object]] = {}

    for projection in projection_order:
        if not projection_columns_exist(df, projection):
            continue
        # endif

        info = PROJECTIONS[projection]
        yield_weight_columns = candidate_yield_columns_for_period(df, prediction_period)
        group_cols = [str(info["min_col"]), str(info["max_col"])]
        sorted_df = df.sort_values(group_cols)

        for _, group in sorted_df.groupby(group_cols, sort=True, dropna=True):
            group_key = (
                float(group[str(info["min_col"])].iloc[0]),
                float(group[str(info["max_col"])].iloc[0]),
            )

            x_plot = average_x_for_group(
                group=group,
                projection_info=info,
                period=prediction_period,
                yield_weight_columns=yield_weight_columns,
            )

            if prediction_mode == "representative":
                rep = representative_kinematics_for_group(
                    group=group,
                    projection=projection,
                    period=prediction_period,
                    no_width_weighting=no_width_weighting,
                    yield_weight_columns=yield_weight_columns,
                )

                if rep.n_rows <= 0:
                    continue
                # endif

                if not (
                    np.isfinite(rep.xB)
                    and np.isfinite(rep.Q2)
                    and np.isfinite(rep.t_abs)
                    and np.isfinite(rep.phi)
                    and np.isfinite(rep.integration_weight_sum)
                ):
                    continue
                # endif

                tasks.append(
                    (
                        task_id,
                        float(rep.xB),
                        float(rep.Q2),
                        float(rep.t_abs),
                        float(rep.phi),
                        float(prediction_ebeam),
                    )
                )

                task_metadata[task_id] = {
                    "mode": "representative",
                    "projection": projection,
                    "key": group_key,
                    "x": float(x_plot),
                    "weight": float(rep.integration_weight_sum),
                }

                task_id += 1

            else:
                xs_col = cross_section_column(prediction_period)

                for row_index, row in group.iterrows():
                    if xs_col in row.index:
                        observed_value, _, _ = parse_tuple3(row[xs_col])

                        if not np.isfinite(observed_value):
                            continue
                        # endif
                    # endif

                    xB = row_average_or_midpoint(
                        row=row,
                        avg_col=f"xBavg, {prediction_period}",
                        min_col="xBmin",
                        max_col="xBmax",
                    )

                    Q2 = row_average_or_midpoint(
                        row=row,
                        avg_col=f"Q2avg, {prediction_period}",
                        min_col="Q2min",
                        max_col="Q2max",
                    )

                    t_abs = row_average_or_midpoint(
                        row=row,
                        avg_col=f"t_abs_avg, {prediction_period}",
                        min_col="t_abs_min",
                        max_col="t_abs_max",
                    )

                    phi = row_average_or_midpoint(
                        row=row,
                        avg_col=f"phiavg, {prediction_period}",
                        min_col="phimin",
                        max_col="phimax",
                    )

                    if not (
                        np.isfinite(xB)
                        and np.isfinite(Q2)
                        and np.isfinite(t_abs)
                        and np.isfinite(phi)
                    ):
                        continue
                    # endif

                    integration_weight = weight_for_row(
                        row=row,
                        integrate_widths=list(info["integrate_widths"]),
                        no_width_weighting=no_width_weighting,
                    )

                    if not np.isfinite(integration_weight):
                        continue
                    # endif

                    tasks.append(
                        (
                            task_id,
                            float(xB),
                            float(Q2),
                            float(t_abs),
                            float(phi),
                            float(prediction_ebeam),
                        )
                    )

                    task_metadata[task_id] = {
                        "mode": "row-sum",
                        "projection": projection,
                        "key": group_key,
                        "x": float(x_plot),
                        "weight": float(integration_weight),
                        "row_index": int(row_index),
                    }

                    task_id += 1
                # endfor
            # endif
        # endfor
    # endfor

    print()
    print("Precomputing KM15 kinematic-dependence prediction points:")
    print(f"  prediction mode       = {prediction_mode}")
    print(f"  prediction period     = {prediction_period}")
    print(f"  prediction Ebeam      = {prediction_ebeam:.6f} GeV")
    print(f"  KM15 tasks            = {len(tasks)}")
    print(f"  parallel              = {not no_parallel_km15 and max_workers > 1}")
    print(f"  max workers           = {max_workers}")

    results_by_task_id = run_km15_prediction_tasks(
        tasks=tasks,
        max_workers=max_workers,
        no_parallel_km15=no_parallel_km15,
    )

    diagnostics.km15_prediction_count += len(tasks)

    if prediction_mode == "representative":
        points_by_projection: Dict[str, List[Point]] = {
            projection: []
            for projection in projection_order
        }

        for tid, metadata in task_metadata.items():
            value, fatal_error_message = results_by_task_id.get(
                tid,
                (math.nan, "missing result"),
            )

            if fatal_error_message is not None:
                diagnostics.km15_prediction_failure_count += 1

                if diagnostics.km15_prediction_failure_message is None:
                    diagnostics.km15_prediction_failure_message = fatal_error_message
                # endif

                continue
            # endif

            if not np.isfinite(value):
                continue
            # endif

            projection = str(metadata["projection"])
            key = metadata["key"]
            x = float(metadata["x"])
            weight = float(metadata["weight"])

            y = value * weight

            points_by_projection[projection].append(
                Point(
                    key=key,
                    x=x,
                    y=y,
                    stat=0.0,
                    sys=0.0,
                )
            )
        # endfor

    else:
        grouped_sums: Dict[Tuple[str, Tuple[float, float]], Dict[str, float]] = {}

        for tid, metadata in task_metadata.items():
            value, fatal_error_message = results_by_task_id.get(
                tid,
                (math.nan, "missing result"),
            )

            if fatal_error_message is not None:
                diagnostics.km15_prediction_failure_count += 1

                if diagnostics.km15_prediction_failure_message is None:
                    diagnostics.km15_prediction_failure_message = fatal_error_message
                # endif

                continue
            # endif

            if not np.isfinite(value):
                continue
            # endif

            projection = str(metadata["projection"])
            key = metadata["key"]
            x = float(metadata["x"])
            weight = float(metadata["weight"])

            group_key = (projection, key)

            if group_key not in grouped_sums:
                grouped_sums[group_key] = {
                    "x": x,
                    "y": 0.0,
                }
            # endif

            grouped_sums[group_key]["y"] += value * weight
        # endfor

        points_by_projection = {
            projection: []
            for projection in projection_order
        }

        for (projection, key), values in grouped_sums.items():
            points_by_projection[projection].append(
                Point(
                    key=key,
                    x=float(values["x"]),
                    y=float(values["y"]),
                    stat=0.0,
                    sys=0.0,
                )
            )
        # endfor
    # endif

    for projection in projection_order:
        points_by_projection[projection].sort(key=lambda p: p.x)
    # endfor

    return points_by_projection


def point_total_error(
    point: Point,
    include_bin_to_bin_sys: bool,
    bin_to_bin_sys_frac: float,
) -> float:
    stat = point.stat if np.isfinite(point.stat) else 0.0

    if not include_bin_to_bin_sys:
        return stat
    # endif

    if not np.isfinite(point.y):
        return stat
    # endif

    sys_bin_to_bin = bin_to_bin_sys_frac * abs(point.y)

    return math.sqrt(stat * stat + sys_bin_to_bin * sys_bin_to_bin)


def points_to_arrays(
    points: List[Point],
    include_bin_to_bin_sys: bool,
    bin_to_bin_sys_frac: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.array([p.x for p in points], dtype=float)
    y = np.array([p.y for p in points], dtype=float)

    stat_err = np.array(
        [
            p.stat if np.isfinite(p.stat) else 0.0
            for p in points
        ],
        dtype=float,
    )

    total_err = np.array(
        [
            point_total_error(
                point=p,
                include_bin_to_bin_sys=include_bin_to_bin_sys,
                bin_to_bin_sys_frac=bin_to_bin_sys_frac,
            )
            for p in points
        ],
        dtype=float,
    )

    return x, y, stat_err, total_err


# -----------------------------------------------------------------------------
# Ratios.
# -----------------------------------------------------------------------------

def ratio_points(numerator: List[Point], denominator: List[Point]) -> List[Point]:
    n_by_key = {p.key: p for p in numerator}
    d_by_key = {p.key: p for p in denominator}

    ratios: List[Point] = []
    common_keys = sorted(set(n_by_key.keys()) & set(d_by_key.keys()))

    for key in common_keys:
        n = n_by_key[key]
        d = d_by_key[key]

        if not np.isfinite(n.y) or not np.isfinite(d.y) or d.y == 0.0:
            continue
        # endif

        r = n.y / d.y

        rel2 = 0.0

        if n.y != 0.0 and np.isfinite(n.stat):
            rel2 += (n.stat / n.y) ** 2
        # endif

        if np.isfinite(d.stat):
            rel2 += (d.stat / d.y) ** 2
        # endif

        ratios.append(
            Point(
                key=n.key,
                x=n.x,
                y=r,
                stat=abs(r) * math.sqrt(rel2),
                sys=0.0,
            )
        )
    # endfor

    return ratios


def weighted_mean_two_points(a: Point, b: Point) -> Optional[Point]:
    values: List[float] = []
    weights: List[float] = []

    for p in [a, b]:
        if not np.isfinite(p.y):
            continue
        # endif

        if np.isfinite(p.stat) and p.stat > 0.0:
            w = 1.0 / (p.stat * p.stat)
        else:
            w = 0.0
        # endif

        if w <= 0.0:
            continue
        # endif

        values.append(float(p.y))
        weights.append(float(w))
    # endfor

    if len(values) == 0:
        return None
    # endif

    values_arr = np.array(values, dtype=float)
    weights_arr = np.array(weights, dtype=float)

    mean = float(np.sum(weights_arr * values_arr) / np.sum(weights_arr))
    stat = float(math.sqrt(1.0 / np.sum(weights_arr)))

    return Point(
        key=a.key,
        x=a.x,
        y=mean,
        stat=stat,
        sys=0.0,
    )


def ratio_to_fa18_scaled_sp19_weighted_mean(
    right_points_by_label: Dict[str, List[Point]],
    label: str,
) -> List[Point]:
    fa_label = "Fa18 Inb, 10.604 GeV"
    sp_label = "Sp19 Inb -> 10.604 GeV, KM15"

    fa_points = {p.key: p for p in right_points_by_label.get(fa_label, [])}
    sp_points = {p.key: p for p in right_points_by_label.get(sp_label, [])}
    target_points = {p.key: p for p in right_points_by_label.get(label, [])}

    ratios: List[Point] = []
    common_keys = sorted(
        set(fa_points.keys())
        & set(sp_points.keys())
        & set(target_points.keys())
    )

    for key in common_keys:
        mean_point = weighted_mean_two_points(fa_points[key], sp_points[key])

        if mean_point is None or mean_point.y == 0.0:
            continue
        # endif

        p = target_points[key]
        r = p.y / mean_point.y

        rel2 = 0.0

        if p.y != 0.0 and np.isfinite(p.stat):
            rel2 += (p.stat / p.y) ** 2
        # endif

        if np.isfinite(mean_point.stat):
            rel2 += (mean_point.stat / mean_point.y) ** 2
        # endif

        ratios.append(
            Point(
                key=p.key,
                x=p.x,
                y=r,
                stat=abs(r) * math.sqrt(rel2),
                sys=0.0,
            )
        )
    # endfor

    return ratios


# -----------------------------------------------------------------------------
# Plotting.
# -----------------------------------------------------------------------------

def plot_points(
    ax,
    points: List[Point],
    label: str,
    ratio: bool,
    include_bin_to_bin_sys: bool,
    bin_to_bin_sys_frac: float,
) -> None:
    if len(points) == 0:
        return
    # endif

    x, y, stat_err, total_err = points_to_arrays(
        points=points,
        include_bin_to_bin_sys=include_bin_to_bin_sys,
        bin_to_bin_sys_frac=bin_to_bin_sys_frac,
    )

    style = SERIES_STYLE.get(label, dict(marker="o", linestyle="-"))
    color = style.get("color", None)

    if include_bin_to_bin_sys:
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
        ax.axhline(
            1.0,
            color="0.35",
            linewidth=1.0,
            linestyle="--",
            zorder=0,
        )
    # endif


def set_ratio_ylim(ax) -> None:
    ax.set_ylim(RATIO_YMIN, RATIO_YMAX)


def make_projection_plot(
    df: pd.DataFrame,
    projection: str,
    output_dir: str,
    no_width_weighting: bool,
    include_bin_to_bin_sys: bool,
    bin_to_bin_sys_frac: float,
    apply_sp19_to_fa18_scaling: bool,
) -> bool:
    if not projection_columns_exist(df, projection):
        print(
            f"Skipping {projection} dependence: missing "
            f"{PROJECTIONS[projection]['min_col']} and/or {PROJECTIONS[projection]['max_col']}"
        )
        return False
    # endif

    info = PROJECTIONS[projection]

    all_needed_periods = []

    for period in LEFT_SERIES + MIDDLE_SERIES + RIGHT_RAW_PERIODS:
        if period not in all_needed_periods:
            all_needed_periods.append(period)
        # endif
    # endfor

    points_by_period = {
        period: integrated_points_for_period(
            df=df,
            projection=projection,
            period=period,
            no_width_weighting=no_width_weighting,
        )
        for period in all_needed_periods
    }

    right_points_by_label: Dict[str, List[Point]] = {
        "Fa18 Inb, 10.604 GeV": integrated_points_for_period(
            df=df,
            projection=projection,
            period="Fa18 Inb",
            no_width_weighting=no_width_weighting,
        ),
        "Sp19 Inb -> 10.604 GeV, KM15": integrated_points_for_period(
            df=df,
            projection=projection,
            period="Sp19 Inb",
            no_width_weighting=no_width_weighting,
            apply_sp19_to_fa18_scaling=apply_sp19_to_fa18_scaling,
        ),
    }

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18.0, 5.5),
        constrained_layout=True,
    )

    fig.suptitle(str(info["title"]), fontsize=16)

    left = axes[0]

    for period in LEFT_SERIES:
        plot_points(
            ax=left,
            points=points_by_period.get(period, []),
            label=period,
            ratio=False,
            include_bin_to_bin_sys=include_bin_to_bin_sys,
            bin_to_bin_sys_frac=bin_to_bin_sys_frac,
        )
    # endfor

    left.set_xlabel(str(info["xlabel"]))
    left.set_ylabel(str(info["ylabel"]))
    left.set_title("Integrated cross sections")
    left.grid(True, alpha=0.25)
    left.legend(fontsize=8, frameon=True)

    middle = axes[1]
    denom_106 = points_by_period.get("10.6 GeV", [])

    for period in MIDDLE_SERIES:
        plot_points(
            ax=middle,
            points=ratio_points(points_by_period.get(period, []), denom_106),
            label=period,
            ratio=True,
            include_bin_to_bin_sys=include_bin_to_bin_sys,
            bin_to_bin_sys_frac=bin_to_bin_sys_frac,
        )
    # endfor

    middle.set_xlabel(str(info["xlabel"]))
    middle.set_ylabel(r"run period / 10.6 GeV")
    middle.set_title("10.6-GeV period consistency")
    middle.grid(True, alpha=0.25)
    middle.legend(fontsize=8, frameon=True)
    set_ratio_ylim(middle)

    right = axes[2]

    for label in RIGHT_SERIES_LABELS:
        plot_points(
            ax=right,
            points=ratio_to_fa18_scaled_sp19_weighted_mean(
                right_points_by_label=right_points_by_label,
                label=label,
            ),
            label=label,
            ratio=True,
            include_bin_to_bin_sys=include_bin_to_bin_sys,
            bin_to_bin_sys_frac=bin_to_bin_sys_frac,
        )
    # endfor

    right.set_xlabel(str(info["xlabel"]))
    right.set_ylabel(r"period / weighted mean")
    right.set_title(r"Fa18 Inb vs Sp19 Inb scaled to 10.604 GeV")
    right.grid(True, alpha=0.25)
    right.legend(fontsize=7, frameon=True)
    set_ratio_ylim(right)

    outbase = os.path.join(
        output_dir,
        f"integrated_{info['outfile_tag']}_dependence",
    )

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
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(14.0, 10.0),
        constrained_layout=True,
    )

    fig.suptitle(
        (
            "KM15 projected cross-section predictions\n"
            rf"CSV representative kinematics from {prediction_period}; "
            rf"$E_{{beam}}={prediction_ebeam:.4f}$ GeV; "
            rf"mode={prediction_mode}"
        ),
        fontsize=15,
    )

    for ax, projection in zip(axes.ravel(), KM15_PREDICTION_PROJECTION_ORDER):
        info = PROJECTIONS[projection]
        points = prediction_points_by_projection.get(projection, [])

        plot_points(
            ax=ax,
            points=points,
            label=KM15_PREDICTION_LABEL,
            ratio=False,
            include_bin_to_bin_sys=False,
            bin_to_bin_sys_frac=0.0,
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
# Diagnostics.
# -----------------------------------------------------------------------------

def print_single_correction_summary(
    title: str,
    count: int,
    correction_sum: float,
    correction_min: float,
    correction_max: float,
) -> None:
    if count <= 0:
        print(f"{title}: no corrections were evaluated")
        return
    # endif

    mean_corr = correction_sum / count

    print(f"{title}:")
    print(f"  correction count   = {count}")
    print(f"  correction mean    = {mean_corr:.6g}")
    print(f"  correction range   = [{correction_min:.6g}, {correction_max:.6g}]")


def print_model_correction_summary(diagnostics: ModelDiagnostics) -> None:
    if not diagnostics.enabled:
        print("Sp19 -> 10.604 GeV beam-energy scaling: disabled")
        return
    # endif

    print("Sp19 -> 10.604 GeV beam-energy scaling:")
    print(f"  source beam energy = {SP19_INB_EBEAM_GEV:.6f} GeV")
    print(f"  target beam energy = {FA18_INB_EBEAM_GEV:.6f} GeV")
    print("  applied correction = KM15")
    print("  BH correction      = diagnostic only")
    print(f"  unique cached correction kinematic points = {diagnostics.sp19_unique_correction_count}")
    print()

    print_single_correction_summary(
        title="KM15 ratio, KM15(10.604 GeV) / KM15(10.1998 GeV)",
        count=diagnostics.km15_correction_count,
        correction_sum=diagnostics.km15_correction_sum,
        correction_min=diagnostics.km15_correction_min,
        correction_max=diagnostics.km15_correction_max,
    )

    print()

    print_single_correction_summary(
        title="BH ratio, BH(10.604 GeV) / BH(10.1998 GeV)",
        count=diagnostics.bh_correction_count,
        correction_sum=diagnostics.bh_correction_sum,
        correction_min=diagnostics.bh_correction_min,
        correction_max=diagnostics.bh_correction_max,
    )

    if diagnostics.bh_correction_count <= 0 and diagnostics.bh_failure_message is not None:
        print()
        print("BH diagnostic failure reason:")
        print(f"  {diagnostics.bh_failure_message}")
    # endif

    if diagnostics.km15_prediction_count > 0:
        print()
        print("KM15 prediction evaluation:")
        print(f"  prediction tasks attempted = {diagnostics.km15_prediction_count}")
        print(f"  prediction task failures   = {diagnostics.km15_prediction_failure_count}")

        if diagnostics.km15_prediction_failure_message is not None:
            print("  first prediction failure:")
            print(diagnostics.km15_prediction_failure_message)
        # endif
    # endif


# -----------------------------------------------------------------------------
# Argument parsing and main.
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make integrated xB, Q2, |t|, phi, and optional detector-angle "
            "dependence plots from a final DVCS pass-2 analysis CSV."
        )
    )

    parser.add_argument(
        "csv_file",
        help="Input final DVCS pass-2 analysis CSV.",
    )

    parser.add_argument(
        "--output-dir",
        default="output/integrated_study",
        help="Directory for output plots. Default: output/integrated_study",
    )

    parser.add_argument(
        "--no-width-weighting",
        action="store_true",
        help="Use raw sums instead of bin-width weighted integrations.",
    )

    parser.add_argument(
        "--phi-degrees",
        action="store_true",
        help="Use phi bin widths in degrees instead of radians for integration weights.",
    )

    parser.add_argument(
        "--include-bin-to-bin-sys",
        action="store_true",
        help=(
            "Draw both statistical error bars and total stat+bin-to-bin-systematic "
            "error bars. The total error is sqrt(stat^2 + (f * y)^2). "
            "The default f is 0.10."
        ),
    )

    parser.add_argument(
        "--bin-to-bin-sys-frac",
        type=float,
        default=0.10,
        help=(
            "Fractional bin-to-bin systematic uncertainty used when "
            "--include-bin-to-bin-sys is enabled. Default: 0.10."
        ),
    )

    parser.add_argument(
        "--no-sp19-km15-energy-scaling",
        action="store_true",
        help=(
            "Disable the KM15 beam-energy scaling of Sp19 Inb from "
            "10.1998 GeV to 10.604 GeV in the right-panel Fa18/Sp19 comparison."
        ),
    )

    parser.add_argument(
        "--compute-bh-diagnostic",
        action="store_true",
        help=(
            "Also compute the BH-only Sp19 -> 10.604 GeV ratio as a diagnostic. "
            "This is not applied to the data and is disabled by default for speed."
        ),
    )

    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help=(
            "Maximum number of process workers for parallel KM15 evaluations. "
            f"Default: {DEFAULT_MAX_WORKERS}."
        ),
    )

    parser.add_argument(
        "--no-parallel-km15",
        action="store_true",
        help="Disable process parallelization for KM15 evaluations.",
    )

    parser.add_argument(
        "--no-km15-prediction-plot",
        action="store_true",
        help="Disable the additional 2x2 KM15 kinematic-dependence prediction plot.",
    )

    parser.add_argument(
        "--km15-prediction-output-name",
        default=DEFAULT_KM15_PREDICTION_OUTPUT_NAME,
        help=(
            "Output PNG filename for the KM15 kinematic-dependence prediction plot. "
            f"Default: {DEFAULT_KM15_PREDICTION_OUTPUT_NAME}"
        ),
    )

    parser.add_argument(
        "--km15-prediction-period",
        default=DEFAULT_KM15_PREDICTION_PERIOD,
        help=(
            "CSV period whose average kinematic columns are used for the KM15 "
            "prediction plot. Default: '10.6 GeV'."
        ),
    )

    parser.add_argument(
        "--km15-prediction-ebeam",
        type=float,
        default=DEFAULT_KM15_PREDICTION_EBEAM_GEV,
        help=(
            "Beam energy in GeV used for the KM15 prediction plot. "
            f"Default: {DEFAULT_KM15_PREDICTION_EBEAM_GEV:.4f}."
        ),
    )

    parser.add_argument(
        "--km15-prediction-mode",
        choices=["representative", "row-sum"],
        default=DEFAULT_KM15_PREDICTION_MODE,
        help=(
            "KM15 projection mode. 'representative' evaluates KM15 once per "
            "projected bin using weighted mean hidden kinematics. 'row-sum' "
            "evaluates KM15 on every underlying multidimensional CSV row and "
            "sums with the data integration weights. Default: representative."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.bin_to_bin_sys_frac < 0.0:
        raise ValueError("--bin-to-bin-sys-frac must be non-negative.")
    # endif

    if args.km15_prediction_ebeam <= 0.0:
        raise ValueError("--km15-prediction-ebeam must be positive.")
    # endif

    if args.max_workers <= 0:
        raise ValueError("--max-workers must be positive.")
    # endif

    if args.max_workers > DEFAULT_MAX_WORKERS:
        print(
            f"Requested --max-workers {args.max_workers}; capping to {DEFAULT_MAX_WORKERS}."
        )
        args.max_workers = DEFAULT_MAX_WORKERS
    # endif

    required = [
        "xBmin",
        "xBmax",
        "Q2min",
        "Q2max",
        "t_abs_min",
        "t_abs_max",
        "phimin",
        "phimax",
    ]

    required += [cross_section_column(period) for period in LEFT_SERIES]

    df = pd.read_csv(args.csv_file)
    require_columns(df, required)

    df = add_width_columns(
        df=df,
        phi_degrees=args.phi_degrees,
    )

    os.makedirs(args.output_dir, exist_ok=True)

    df, diagnostics = precompute_sp19_to_fa18_km15_scales(
        df=df,
        enabled=not args.no_sp19_km15_energy_scaling,
        compute_bh_diagnostic=args.compute_bh_diagnostic,
        max_workers=args.max_workers,
        no_parallel_km15=args.no_parallel_km15,
    )

    made_projection_count = 0

    for projection in DATA_PROJECTION_ORDER:
        made_plot = make_projection_plot(
            df=df,
            projection=projection,
            output_dir=args.output_dir,
            no_width_weighting=args.no_width_weighting,
            include_bin_to_bin_sys=args.include_bin_to_bin_sys,
            bin_to_bin_sys_frac=args.bin_to_bin_sys_frac,
            apply_sp19_to_fa18_scaling=not args.no_sp19_km15_energy_scaling,
        )

        if made_plot:
            made_projection_count += 1
        # endif
    # endfor

    if not args.no_km15_prediction_plot:
        prediction_points_by_projection = build_km15_prediction_points_all_projections(
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
            prediction_points_by_projection=prediction_points_by_projection,
            output_dir=args.output_dir,
            output_name=args.km15_prediction_output_name,
            prediction_period=args.km15_prediction_period,
            prediction_ebeam=args.km15_prediction_ebeam,
            prediction_mode=args.km15_prediction_mode,
        )
    # endif

    if args.include_bin_to_bin_sys:
        print(
            "Using plotted y-errors: statistical bars plus outer total bars "
            f"sqrt(stat^2 + ({args.bin_to_bin_sys_frac:.4f} * y)^2)"
        )
    else:
        print("Using plotted y-errors: statistical only")
    # endif

    print(f"All ratio-plot y axes fixed to [{RATIO_YMIN:.2f}, {RATIO_YMAX:.2f}]")

    if args.no_km15_prediction_plot:
        print("KM15 kinematic-dependence prediction plot: disabled")
    else:
        print(
            "KM15 kinematic-dependence prediction plot:"
            f" period = {args.km15_prediction_period},"
            f" Ebeam = {args.km15_prediction_ebeam:.6f} GeV,"
            f" mode = {args.km15_prediction_mode}"
        )
    # endif

    print_model_correction_summary(diagnostics)

    print(f"Wrote {made_projection_count} integrated dependence PNG plots to: {args.output_dir}")


if __name__ == "__main__":
    main()
# endif