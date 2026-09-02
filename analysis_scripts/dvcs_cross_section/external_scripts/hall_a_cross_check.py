#!/usr/bin/env python3
"""
hall_a_cross_check.py

Standalone Hall A / CLAS12 DVCS cross-section cross-check.

This script reads the final DVCS pass-2 analysis CSV and compares the CLAS12
cross sections in one overlapping bin against the Hall A Kin363 cross-section
table. It also overlays the previous pass-1 Hall B / CLAS12 Fa18 result from
all_bin_v3.csv by default.

The script also makes a second 1x2 canvas comparing:

  Fa18 Inb, 10.604 GeV
  Sp19 Inb -> 10.604 GeV, KM15
  Hall A

where the Sp19 Inb points are scaled row-by-row from 10.1998 GeV to 10.6040 GeV
using:

  C_i = KM15(E_target, xB_i, Q2_i, |t|_i, phi_i)
      / KM15(E_source, xB_i, Q2_i, |t|_i, phi_i)

and:

  sigma_i_scaled = C_i * sigma_i
  stat_i_scaled  = |C_i| * stat_i
  sys_i_scaled   = |C_i| * sys_i

The Hall A values are hard-coded from the right-most column of the Hall A table:

  <xB>  = 0.371
  <Q2>  = 4.568 GeV^2
  <t'>  = -0.303 GeV^2

with cross sections listed in pb/GeV^4 as:

  sigma +/- stat ^{+ sys_up}_{- sys_down}

Since the Hall A value is t' rather than t, the script computes:

  t = t_min + t'

using the exact DVCS t_min formula. For the Hall A central kinematics, this
places the comparison around |t| ~= 0.483 GeV^2, so the corresponding CLAS12
bin is selected by default using:

  Bin Name = 138

which corresponds nominally to:

  xB  in [0.357, 0.446]
  Q2  in [4.326, 5.761] GeV^2
  |t| in [0.40, 0.60] GeV^2

The first plot is a 1x2 canvas:

  left:   pass-2 integrated 10.6 GeV, pass-1 Fa18, and Hall A cross sections vs phi
  right:  pass-2 integrated 10.6 GeV / Hall A and pass-1 Fa18 / Hall A ratios vs phi

The second plot is a 1x2 canvas:

  left:   Fa18 Inb, Sp19 Inb -> 10.604 GeV with KM15, and Hall A cross sections vs phi
  right:  Fa18 Inb / Hall A and scaled Sp19 Inb / Hall A ratios vs phi

All plots use stat+sys error bars.

For Hall A:
  total upper/lower errors are:
    err_up   = sqrt(stat^2 + sys_up^2)
    err_down = sqrt(stat^2 + sys_down^2)

For pass-2 CLAS12 / Hall B:
  the CSV columns are read as:
    (value, stat, sys_csv)

  The point-to-point systematic is read directly from:
    Syst. err (point-to-point total)

  This is the CSV point-to-point total associated with the final integrated
  10.6 GeV result. The script does not add an ad hoc fractional uncertainty.

  The period-dependent scale systematic is read from:
    normed cross sections, ep->epg, exp, <period>, unpol, total scale sys

  The plotted pass-2 total uncertainty is:
    err_total = sqrt(stat^2 + sys_point_to_point^2 + (scale_frac*sigma)^2)

For pass-1 Fa18:
  the CSV columns are read as:
    cross sections, ep->epg, exp
    cross sections, ep->epg, exp, stat. unc.
    cross sections, ep->epg, exp, syst. unc. (up)
    cross sections, ep->epg, exp, syst. unc. (down)

  The script additionally includes a 31% normalization uncertainty:
    norm = 0.31 * sigma

  The pass-1 total upper/lower errors are:
    err_up   = sqrt(stat^2 + sys_up^2 + norm^2)
    err_down = sqrt(stat^2 + sys_down^2 + norm^2)

Important unit convention:
  The Hall A table is in pb/GeV^4. The current CLAS12 CSV cross-section values
  in these files appear to be in nb/GeV^4 relative to the Hall A table, so the
  script multiplies CLAS12 cross sections by 1000 by default.

  Override with:
    --clas12-scale 1.0
    --pass1-scale 1.0

Interpolation convention:
  The CLAS12 and Hall A phi points are not necessarily identical. For the ratio
  panel and chi2/ndf, Hall A is periodically linearly interpolated to each CLAS12
  phi value before forming:

    ratio = CLAS12(phi_CLAS12) / HallA_interpolated(phi_CLAS12)

Output:

  output/hall_a_cross_check/hall_a_cross_check.png
  output/hall_a_cross_check/hall_a_fa18_sp19_km15_cross_check.png
"""

import argparse
import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pass1_published_loader import load_published_pass1_dataframe


# -----------------------------------------------------------------------------
# Default CLAS12 bin corresponding to the converted Hall A t' overlap point.
# -----------------------------------------------------------------------------

DEFAULT_BIN_NAME = 138

TARGET_XB_MIN = 0.357
TARGET_XB_MAX = 0.446

TARGET_Q2_MIN = 4.326
TARGET_Q2_MAX = 5.761

TARGET_T_ABS_MIN = 0.40
TARGET_T_ABS_MAX = 0.60


# -----------------------------------------------------------------------------
# Default paths and output names.
# -----------------------------------------------------------------------------

DEFAULT_PASS1_CSV = (
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "dvcs_cross_section/imports/clasdb_E214M1.txt"
)

DEFAULT_OUTPUT_DIR = "output/hall_a_cross_check"
DEFAULT_OUTPUT_NAME = "hall_a_cross_check.png"
DEFAULT_FA18_SP19_OUTPUT_NAME = "hall_a_fa18_sp19_km15_cross_check.png"


# -----------------------------------------------------------------------------
# Hall A kinematics for the selected table column.
# -----------------------------------------------------------------------------

HALL_A_XB = 0.371
HALL_A_Q2 = 4.568
HALL_A_TPRIME = -0.303
PROTON_MASS_GEV = 0.9382720813


# -----------------------------------------------------------------------------
# Beam energies.
# -----------------------------------------------------------------------------

FA18_INB_EBEAM_GEV = 10.6040
SP19_INB_EBEAM_GEV = 10.1998


# -----------------------------------------------------------------------------
# Default uncertainty/systematic choices.
# -----------------------------------------------------------------------------

DEFAULT_CLAS12_SCALE = 1000.0
DEFAULT_PASS1_SCALE = 1000.0
DEFAULT_PASS1_NORM_SYS_FRAC = 0.31


# -----------------------------------------------------------------------------
# Plot labels and period mapping.
# -----------------------------------------------------------------------------

PASS2_COMBINED_DISPLAY_LABEL = "pass-2 10.6 GeV"
PASS2_COMBINED_CSV_PERIOD = "10.6 GeV"

PASS1_DISPLAY_LABEL = "pass-1 Fa18"

# Series shown in the primary Hall A comparison panel.  Keep this definition
# separate from the dedicated Fa18-vs-Sp19 comparison below.
TOP_COMPARISON_SERIES = [
    PASS2_COMBINED_DISPLAY_LABEL,
    PASS1_DISPLAY_LABEL,
]

FA18_INB_DISPLAY_LABEL = "Fa18 Inb"
SP19_INB_DISPLAY_LABEL = "Sp19 Inb"
SP19_KM15_DISPLAY_LABEL = "Sp19 Inb -> 10.604 GeV, KM15"

FA18_SP19_COMPARISON_SERIES = [
    FA18_INB_DISPLAY_LABEL,
    SP19_KM15_DISPLAY_LABEL,
]

ALL_PASS2_CSV_PERIODS = [
    PASS2_COMBINED_CSV_PERIOD,
    "Fa18 Inb",
    "Sp19 Inb",
]


# Hall A is always black.
# pass-1 Fa18 is drawn as an open circle with a distinct cyan color.
SERIES_STYLE = {
    "Hall A": dict(marker="s", linestyle="None", color="black"),
    PASS2_COMBINED_DISPLAY_LABEL: dict(marker="o", linestyle="None", color="tab:red"),
    PASS1_DISPLAY_LABEL: dict(marker="o", linestyle="None", color="tab:cyan"),
    "Sp18 Inb": dict(marker="D", linestyle="None", color="tab:green"),
    "Sp18 Out": dict(marker="P", linestyle="None", color="tab:purple"),
    "Fa18 Inb": dict(marker="^", linestyle="None", color="tab:blue"),
    "Fa18 Out": dict(marker="v", linestyle="None", color="tab:orange"),
    "Sp19 Inb": dict(marker="s", linestyle="None", color="tab:red"),
    SP19_KM15_DISPLAY_LABEL: dict(marker="s", linestyle="None", color="tab:red"),
}


# -----------------------------------------------------------------------------
# Hall A Kin363 right-most column:
#
# <xB>  = 0.371
# <Q2>  = 4.568 GeV^2
# <t'>  = -0.303 GeV^2
#
# sigma values in pb/GeV^4.
# -----------------------------------------------------------------------------

HALL_A_DATA = [
    # phi_deg, sigma, stat, sys_up, sys_down
    (7.5,   11.78, 3.54, 0.54, 1.05),
    (22.5,   5.90, 1.74, 0.53, 0.61),
    (37.5,   5.60, 0.91, 0.41, 0.23),
    (52.5,   4.85, 0.54, 0.19, 0.11),
    (67.5,   4.16, 0.39, 0.27, 0.09),
    (82.5,   3.39, 0.31, 0.01, 0.19),
    (97.5,   2.94, 0.28, 0.07, 0.08),
    (112.5,  2.87, 0.24, 0.07, 0.05),
    (127.5,  2.34, 0.20, 0.10, 0.00),
    (142.5,  2.67, 0.20, 0.04, 0.09),
    (157.5,  2.15, 0.19, 0.10, 0.12),
    (172.5,  2.49, 0.21, 0.08, 0.09),
    (187.5,  2.40, 0.20, 0.07, 0.06),
    (202.5,  2.11, 0.19, 0.07, 0.00),
    (217.5,  2.48, 0.20, 0.03, 0.07),
    (232.5,  2.60, 0.21, 0.13, 0.07),
    (247.5,  3.18, 0.25, 0.09, 0.07),
    (262.5,  3.48, 0.28, 0.10, 0.07),
    (277.5,  3.73, 0.32, 0.15, 0.06),
    (292.5,  4.83, 0.39, 0.11, 0.19),
    (307.5,  5.34, 0.51, 0.18, 0.18),
    (322.5,  7.01, 0.87, 0.41, 0.42),
    (337.5, 10.96, 1.91, 0.72, 0.71),
    (352.5, 12.76, 3.21, 1.39, 1.70),
]


FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)


@dataclass
class DataPoint:
    phi: float
    sigma: float
    err_low: float
    err_high: float
    stat: float = 0.0
    sys: float = 0.0
    sys_csv: float = 0.0
    sys_est: float = 0.0
    norm: float = 0.0


@dataclass
class KinematicAverages:
    xB: float
    Q2: float
    t_abs: float
    weight_period: str
    n_points: int


@dataclass
class ModelContext:
    enabled: bool
    g: object = None
    th_km15: object = None
    correction_min: float = math.inf
    correction_max: float = -math.inf
    correction_sum: float = 0.0
    correction_count: int = 0
    correction_cache: Dict[Tuple[float, float, float, float], float] = field(default_factory=dict)


# -----------------------------------------------------------------------------
# Kinematic helpers.
# -----------------------------------------------------------------------------

def compute_tmin_exact(xb: float, q2: float, mass: float = PROTON_MASS_GEV) -> float:
    """
    Compute exact DVCS t_min.

    Convention:
      t is negative.
      t_min is the least negative kinematically allowed t value.

    Formula:
      epsilon^2 = 4 M^2 xB^2 / Q^2

      t_min = -Q^2 *
        [2(1-xB)(1-sqrt(1+epsilon^2)) + epsilon^2]
        /
        [4 xB(1-xB) + epsilon^2]
    """

    eps2 = 4.0 * mass * mass * xb * xb / q2

    numerator = 2.0 * (1.0 - xb) * (1.0 - math.sqrt(1.0 + eps2)) + eps2
    denominator = 4.0 * xb * (1.0 - xb) + eps2

    return -q2 * numerator / denominator


def hall_a_converted_t() -> float:
    """
    Convert the Hall A quoted t' value into actual t.

    Convention:
      t' = t - t_min

    Therefore:
      t = t_min + t'
    """

    tmin = compute_tmin_exact(HALL_A_XB, HALL_A_Q2)
    return tmin + HALL_A_TPRIME


# -----------------------------------------------------------------------------
# Gepard / KM15 helpers.
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

    Do not pass W together with xB and Q2. Gepard treats {xB, W, Q2} as an
    overdetermined kinematic set. Instead, give xB and Q2 and let Gepard fill
    the basic missing DIS kinematics.
    """

    xB = float(xB)
    Q2 = float(Q2)
    t_abs = float(t_abs)
    phi_deg_trento = float(phi_deg_trento)
    ebeam = float(ebeam)

    pt = g.DataPoint(
        xB=xB,
        Q2=Q2,
        t=-abs(t_abs),
        phi=math.radians(phi_deg_trento),
        frame="Trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=ebeam,
        in1charge=-1,
        in1polarization=0,
        observable="XS",
    )

    if hasattr(pt, "to_conventions"):
        pt.to_conventions()
    # endif

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


def update_km15_correction_diagnostics(
    model_context: ModelContext,
    correction: float,
) -> None:
    if not np.isfinite(correction):
        return
    # endif

    model_context.correction_min = min(model_context.correction_min, correction)
    model_context.correction_max = max(model_context.correction_max, correction)
    model_context.correction_sum += correction
    model_context.correction_count += 1


def sp19_to_fa18_km15_correction(
    model_context: ModelContext,
    xB: float,
    Q2: float,
    t_abs: float,
    phi: float,
) -> float:
    """
    Compute and cache:

      KM15(10.604 GeV) / KM15(10.1998 GeV)

    at the row-level Sp19 kinematics.
    """

    if not model_context.enabled:
        return 1.0
    # endif

    if not (
        np.isfinite(xB)
        and np.isfinite(Q2)
        and np.isfinite(t_abs)
        and np.isfinite(phi)
    ):
        return math.nan
    # endif

    key = correction_cache_key(
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi=phi,
    )

    if key in model_context.correction_cache:
        correction = model_context.correction_cache[key]
        update_km15_correction_diagnostics(model_context, correction)
        return correction
    # endif

    numerator = km15_xs(
        model_context=model_context,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi,
        ebeam=FA18_INB_EBEAM_GEV,
    )

    denominator = km15_xs(
        model_context=model_context,
        xB=xB,
        Q2=Q2,
        t_abs=t_abs,
        phi_deg_trento=phi,
        ebeam=SP19_INB_EBEAM_GEV,
    )

    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator == 0.0:
        return math.nan
    # endif

    correction = numerator / denominator
    model_context.correction_cache[key] = correction
    update_km15_correction_diagnostics(model_context, correction)

    return correction


# -----------------------------------------------------------------------------
# Parsing and CSV helpers.
# -----------------------------------------------------------------------------

def parse_tuple3(value) -> Tuple[float, float, float]:
    """
    Parse a CSV cell holding either NaN, scalar, or '(value, stat, sys)'.

    Returns
    -------
    value, stat, sys : tuple of floats
    """

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
    """
    Parse the first number from a cell.

    Useful for cells that may either be scalars or tuple-like strings.
    """

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


def pass2_cross_section_column(csv_period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {csv_period}, unpol"


def pass2_point_to_point_column() -> str:
    return "Syst. err (point-to-point total)"


def pass2_total_scale_column(csv_period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {csv_period}, unpol, total scale sys"


def pass2_scale_systematic_period(csv_period: str) -> str:
    """
    Return the CSV group whose correlated scale systematic applies to a plotted
    pass-2 series.

    The production systematics schema intentionally stores correlated scale
    uncertainties for the combined groups (10.6 GeV, Fa18, Sp18, Sp19 Inb),
    not separately for Fa18 Inb/Out or Sp18 Inb/Out.  Individual-period
    cross-check plots therefore use the corresponding run-group scale nuisance.
    """
    if csv_period in {"Fa18 Inb", "Fa18 Out"}:
        return "Fa18"
    # endif
    if csv_period in {"Sp18 Inb", "Sp18 Out"}:
        return "Sp18"
    # endif
    return csv_period


def pass2_scale_systematic_column(csv_period: str) -> str:
    return pass2_total_scale_column(pass2_scale_systematic_period(csv_period))


def average_phi_column(csv_period: str) -> str:
    return f"phiavg, {csv_period}"


def avg_column(quantity: str, csv_period: str) -> str:
    return f"{quantity}, {csv_period}"


def require_columns(df: pd.DataFrame, columns: Iterable[str], context: str) -> None:
    missing = [c for c in columns if c not in df.columns]

    if missing:
        raise KeyError(
            f"The {context} CSV is missing required columns:\n  "
            + "\n  ".join(missing)
        )
    # endif


def float_close(series: pd.Series, value: float, tolerance: float = 1.0e-6) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    return np.isclose(numeric, value, rtol=0.0, atol=tolerance)


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


def select_overlap_bin(df: pd.DataFrame, bin_name: Optional[int], context: str) -> pd.DataFrame:
    """
    Select the bin overlapping the converted Hall A t' point.

    Preferred selection:
      Bin Name == 138

    Fallback selection:
      exact bin edges with a small tolerance.
    """

    selected = pd.DataFrame()

    if bin_name is not None and "Bin Name" in df.columns:
        bin_numeric = pd.to_numeric(df["Bin Name"], errors="coerce")
        selected = df.loc[bin_numeric == int(bin_name)].copy()
    # endif

    if len(selected) == 0:
        mask = (
            float_close(df["xBmin"], TARGET_XB_MIN)
            & float_close(df["xBmax"], TARGET_XB_MAX)
            & float_close(df["Q2min"], TARGET_Q2_MIN)
            & float_close(df["Q2max"], TARGET_Q2_MAX)
            & float_close(df["t_abs_min"], TARGET_T_ABS_MIN)
            & float_close(df["t_abs_max"], TARGET_T_ABS_MAX)
        )

        selected = df.loc[mask].copy()
    # endif

    if len(selected) == 0:
        raise RuntimeError(
            f"Could not find the requested Hall A overlap bin in the {context} CSV:\n"
            f"  Bin Name = {bin_name}\n"
            f"  xB  [{TARGET_XB_MIN}, {TARGET_XB_MAX}]\n"
            f"  Q2  [{TARGET_Q2_MIN}, {TARGET_Q2_MAX}]\n"
            f"  |t| [{TARGET_T_ABS_MIN}, {TARGET_T_ABS_MAX}]\n"
        )
    # endif

    selected = selected.sort_values(["phimin", "phimax"]).copy()

    return selected


# -----------------------------------------------------------------------------
# Point extraction.
# -----------------------------------------------------------------------------

def hall_a_points() -> List[DataPoint]:
    points: List[DataPoint] = []

    for phi, sigma, stat, sys_up, sys_down in HALL_A_DATA:
        err_high = math.sqrt(stat * stat + sys_up * sys_up)
        err_low = math.sqrt(stat * stat + sys_down * sys_down)

        points.append(
            DataPoint(
                phi=phi,
                sigma=sigma,
                err_low=err_low,
                err_high=err_high,
                stat=stat,
                sys=0.5 * (abs(sys_up) + abs(sys_down)),
                sys_csv=0.5 * (abs(sys_up) + abs(sys_down)),
                sys_est=0.0,
                norm=0.0,
            )
        )
    # endfor

    points.sort(key=lambda p: p.phi)

    return points


def pass2_points_for_period(
    selected: pd.DataFrame,
    csv_period: str,
    clas12_scale: float,
    model_context: Optional[ModelContext] = None,
    apply_sp19_to_fa18_km15_scaling: bool = False,
) -> List[DataPoint]:
    """
    Extract pass-2 CLAS12 phi-dependent points for one period.

    If apply_sp19_to_fa18_km15_scaling is true, the point is scaled row-by-row
    by KM15(10.604 GeV) / KM15(10.1998 GeV) before applying the external
    clas12_scale and before computing the estimated fractional systematic.
    """

    xs_col = pass2_cross_section_column(csv_period)
    phi_avg_col = average_phi_column(csv_period)

    if xs_col not in selected.columns:
        return []
    # endif

    points: List[DataPoint] = []

    sorted_df = selected.sort_values("phimin")

    for _, row in sorted_df.iterrows():
        sigma, stat, _ = parse_tuple3(row[xs_col])
        sys_csv = parse_scalar_from_cell(row.get(pass2_point_to_point_column(), math.nan))
        scale_frac = parse_scalar_from_cell(row.get(pass2_scale_systematic_column(csv_period), math.nan))

        if not np.isfinite(sigma):
            continue
        # endif

        if not np.isfinite(stat):
            stat = 0.0
        # endif

        if not np.isfinite(sys_csv):
            sys_csv = 0.0
        # endif
        if not np.isfinite(scale_frac):
            scale_frac = 0.0
        # endif

        if phi_avg_col in selected.columns:
            phi = parse_scalar_from_cell(row[phi_avg_col])
        else:
            phi = math.nan
        # endif

        if not np.isfinite(phi):
            phi_min = parse_scalar_from_cell(row["phimin"])
            phi_max = parse_scalar_from_cell(row["phimax"])
            phi = 0.5 * (phi_min + phi_max)
        # endif

        if apply_sp19_to_fa18_km15_scaling:
            if model_context is None or not model_context.enabled:
                raise RuntimeError(
                    "Sp19 -> Fa18 KM15 scaling requested, but the model context is disabled."
                )
            # endif

            xB = row_average_or_midpoint(
                row=row,
                avg_col=avg_column("xBavg", csv_period),
                min_col="xBmin",
                max_col="xBmax",
            )

            Q2 = row_average_or_midpoint(
                row=row,
                avg_col=avg_column("Q2avg", csv_period),
                min_col="Q2min",
                max_col="Q2max",
            )

            t_abs = row_average_or_midpoint(
                row=row,
                avg_col=avg_column("t_abs_avg", csv_period),
                min_col="t_abs_min",
                max_col="t_abs_max",
            )

            correction = sp19_to_fa18_km15_correction(
                model_context=model_context,
                xB=xB,
                Q2=Q2,
                t_abs=t_abs,
                phi=phi,
            )

            if not np.isfinite(correction):
                continue
            # endif

            sigma *= correction
            stat *= abs(correction)
            sys_csv *= abs(correction)
        # endif

        sigma *= clas12_scale
        stat *= clas12_scale
        sys_csv *= clas12_scale

        # The CSV point-to-point total is used directly. No ad hoc fractional
        # bin-to-bin uncertainty is added by this script.
        sys_est = 0.0
        scale_abs = abs(scale_frac * sigma)
        sys_total = math.sqrt(sys_csv * sys_csv + scale_abs * scale_abs)
        err_total = math.sqrt(stat * stat + sys_total * sys_total)

        points.append(
            DataPoint(
                phi=phi,
                sigma=sigma,
                err_low=err_total,
                err_high=err_total,
                stat=stat,
                sys=sys_total,
                sys_csv=sys_csv,
                sys_est=sys_est,
                norm=0.0,
            )
        )
    # endfor

    points.sort(key=lambda p: p.phi)

    return points


def pass1_points(
    selected: pd.DataFrame,
    pass1_scale: float,
    pass1_norm_sys_frac: float,
) -> List[DataPoint]:
    """
    Extract published pass-1 Fa18 points; legacy CSV supplies bin boundaries only.

    Required pass-1 columns:
      cross sections, ep->epg, exp
      cross sections, ep->epg, exp, stat. unc.
      cross sections, ep->epg, exp, syst. unc. (up)
      cross sections, ep->epg, exp, syst. unc. (down)

    The pass-1 values are multiplied by pass1_scale.

    The additional normalization uncertainty is:
      norm = pass1_norm_sys_frac * abs(sigma)

    Total asymmetric errors:
      err_up   = sqrt(stat^2 + sys_up^2 + norm^2)
      err_down = sqrt(stat^2 + sys_down^2 + norm^2)
    """

    required = [
        "cross sections, ep->epg, exp",
        "cross sections, ep->epg, exp, stat. unc.",
        "cross sections, ep->epg, exp, syst. unc. (up)",
        "cross sections, ep->epg, exp, syst. unc. (down)",
        "phiavg",
        "phimin",
        "phimax",
    ]

    for col in required:
        if col not in selected.columns:
            return []
        # endif
    # endfor

    points: List[DataPoint] = []

    sorted_df = selected.sort_values("phimin")

    for _, row in sorted_df.iterrows():
        sigma = parse_scalar_from_cell(row["cross sections, ep->epg, exp"])
        stat = parse_scalar_from_cell(row["cross sections, ep->epg, exp, stat. unc."])
        sys_up = parse_scalar_from_cell(row["cross sections, ep->epg, exp, syst. unc. (up)"])
        sys_down = parse_scalar_from_cell(row["cross sections, ep->epg, exp, syst. unc. (down)"])

        if not np.isfinite(sigma):
            continue
        # endif

        if not np.isfinite(stat):
            stat = 0.0
        # endif

        if not np.isfinite(sys_up):
            sys_up = 0.0
        # endif

        if not np.isfinite(sys_down):
            sys_down = 0.0
        # endif

        sigma *= pass1_scale
        stat *= pass1_scale
        sys_up *= pass1_scale
        sys_down *= pass1_scale

        norm = pass1_norm_sys_frac * abs(sigma)

        err_high = math.sqrt(stat * stat + sys_up * sys_up + norm * norm)
        err_low = math.sqrt(stat * stat + sys_down * sys_down + norm * norm)

        phi = parse_scalar_from_cell(row["phiavg"])

        if not np.isfinite(phi):
            phi_min = parse_scalar_from_cell(row["phimin"])
            phi_max = parse_scalar_from_cell(row["phimax"])
            phi = 0.5 * (phi_min + phi_max)
        # endif

        points.append(
            DataPoint(
                phi=phi,
                sigma=sigma,
                err_low=err_low,
                err_high=err_high,
                stat=stat,
                sys=0.5 * (abs(sys_up) + abs(sys_down)),
                sys_csv=0.5 * (abs(sys_up) + abs(sys_down)),
                sys_est=0.0,
                norm=norm,
            )
        )
    # endfor

    points.sort(key=lambda p: p.phi)

    return points


# -----------------------------------------------------------------------------
# Averages, interpolation, ratios, chi2.
# -----------------------------------------------------------------------------

def compute_weighted_kinematic_averages(
    selected: pd.DataFrame,
    csv_period: str = "10.6 GeV",
    fallback_periods: Optional[List[str]] = None,
) -> KinematicAverages:
    """
    Compute weighted average pass-2 CLAS12 kinematics from the selected CSV rows.
    """

    if fallback_periods is None:
        fallback_periods = [
            "10.6 GeV",
            "Fa18 Inb",
            "Fa18 Out",
            "Sp18 Inb",
            "Sp18 Out",
            "Sp19 Inb",
        ]
    # endif

    candidate_periods = [csv_period]

    for p in fallback_periods:
        if p not in candidate_periods:
            candidate_periods.append(p)
        # endif
    # endfor

    chosen_period = csv_period

    for p in candidate_periods:
        required_avg_cols = [
            avg_column("xBavg", p),
            avg_column("Q2avg", p),
            avg_column("t_abs_avg", p),
        ]

        if all(col in selected.columns for col in required_avg_cols):
            chosen_period = p
            break
        # endif
    # endfor

    xb_col = avg_column("xBavg", chosen_period)
    q2_col = avg_column("Q2avg", chosen_period)
    t_col = avg_column("t_abs_avg", chosen_period)
    xs_col = pass2_cross_section_column(chosen_period)

    def fallback_midpoint(min_col: str, max_col: str) -> float:
        lo = pd.to_numeric(selected[min_col], errors="coerce")
        hi = pd.to_numeric(selected[max_col], errors="coerce")

        if len(lo) > 0 and np.isfinite(lo.iloc[0]) and np.isfinite(hi.iloc[0]):
            return 0.5 * (float(lo.iloc[0]) + float(hi.iloc[0]))
        # endif

        return math.nan

    if xb_col not in selected.columns:
        xb_value = fallback_midpoint("xBmin", "xBmax")
        q2_value = fallback_midpoint("Q2min", "Q2max")
        t_value = fallback_midpoint("t_abs_min", "t_abs_max")

        return KinematicAverages(
            xB=xb_value,
            Q2=q2_value,
            t_abs=t_value,
            weight_period="bin midpoint fallback",
            n_points=0,
        )
    # endif

    weighted_sums = {
        "xB": 0.0,
        "Q2": 0.0,
        "t_abs": 0.0,
    }

    unweighted_values = {
        "xB": [],
        "Q2": [],
        "t_abs": [],
    }

    weight_sum = 0.0
    n_points = 0

    for _, row in selected.iterrows():
        xB = parse_scalar_from_cell(row[xb_col])
        Q2 = parse_scalar_from_cell(row[q2_col])
        t_abs = parse_scalar_from_cell(row[t_col])

        if not (np.isfinite(xB) and np.isfinite(Q2) and np.isfinite(t_abs)):
            continue
        # endif

        unweighted_values["xB"].append(xB)
        unweighted_values["Q2"].append(Q2)
        unweighted_values["t_abs"].append(t_abs)

        weight = 1.0

        if xs_col in selected.columns:
            _, stat, _ = parse_tuple3(row[xs_col])

            if np.isfinite(stat) and stat > 0.0:
                weight = 1.0 / (stat * stat)
            else:
                weight = 1.0
            # endif
        # endif

        weighted_sums["xB"] += weight * xB
        weighted_sums["Q2"] += weight * Q2
        weighted_sums["t_abs"] += weight * t_abs
        weight_sum += weight
        n_points += 1
    # endfor

    if weight_sum > 0.0:
        return KinematicAverages(
            xB=weighted_sums["xB"] / weight_sum,
            Q2=weighted_sums["Q2"] / weight_sum,
            t_abs=weighted_sums["t_abs"] / weight_sum,
            weight_period=chosen_period,
            n_points=n_points,
        )
    # endif

    if len(unweighted_values["xB"]) > 0:
        return KinematicAverages(
            xB=float(np.mean(unweighted_values["xB"])),
            Q2=float(np.mean(unweighted_values["Q2"])),
            t_abs=float(np.mean(unweighted_values["t_abs"])),
            weight_period=f"{chosen_period}, unweighted fallback",
            n_points=len(unweighted_values["xB"]),
        )
    # endif

    return KinematicAverages(
        xB=fallback_midpoint("xBmin", "xBmax"),
        Q2=fallback_midpoint("Q2min", "Q2max"),
        t_abs=fallback_midpoint("t_abs_min", "t_abs_max"),
        weight_period="bin midpoint fallback",
        n_points=0,
    )


def points_to_arrays(points: List[DataPoint]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    phi = np.array([p.phi for p in points], dtype=float)
    sigma = np.array([p.sigma for p in points], dtype=float)
    err_low = np.array([p.err_low for p in points], dtype=float)
    err_high = np.array([p.err_high for p in points], dtype=float)

    return phi, sigma, err_low, err_high


def periodic_interp(phi_query: float, phi_values: np.ndarray, y_values: np.ndarray) -> float:
    """
    Periodically linearly interpolate y(phi) at phi_query.
    """

    phi_mod = phi_query % 360.0

    phi_base = np.array(phi_values, dtype=float)
    y_base = np.array(y_values, dtype=float)

    order = np.argsort(phi_base)
    phi_base = phi_base[order]
    y_base = y_base[order]

    phi_ext = np.concatenate(
        [
            phi_base - 360.0,
            phi_base,
            phi_base + 360.0,
        ]
    )

    y_ext = np.concatenate(
        [
            y_base,
            y_base,
            y_base,
        ]
    )

    return float(np.interp(phi_mod, phi_ext, y_ext))


def interpolate_hall_a_point(phi_query: float, hall_a_points_input: List[DataPoint]) -> DataPoint:
    """
    Interpolate Hall A sigma, lower error, and upper error to the requested phi.
    """

    phi, sigma, err_low, err_high = points_to_arrays(hall_a_points_input)

    sigma_interp = periodic_interp(phi_query, phi, sigma)
    err_low_interp = periodic_interp(phi_query, phi, err_low)
    err_high_interp = periodic_interp(phi_query, phi, err_high)

    return DataPoint(
        phi=phi_query,
        sigma=sigma_interp,
        err_low=err_low_interp,
        err_high=err_high_interp,
        stat=0.0,
        sys=0.0,
        sys_csv=0.0,
        sys_est=0.0,
        norm=0.0,
    )


def symmetric_error(point: DataPoint) -> float:
    return 0.5 * (point.err_low + point.err_high)


def chi2_ndf_to_hall_a(
    comparison_points: List[DataPoint],
    hall_a_points_input: List[DataPoint],
) -> Tuple[float, int, float]:
    """
    Compute chi2/ndf comparing one CLAS12 dataset to Hall A interpolated to the
    CLAS12 phi points.
    """

    chi2 = 0.0
    ndf = 0

    for b in comparison_points:
        a = interpolate_hall_a_point(
            phi_query=b.phi,
            hall_a_points_input=hall_a_points_input,
        )

        if not np.isfinite(a.sigma):
            continue
        # endif

        err_b = symmetric_error(b)
        err_a = symmetric_error(a)

        variance = err_b * err_b + err_a * err_a

        if variance <= 0.0 or not np.isfinite(variance):
            continue
        # endif

        residual = b.sigma - a.sigma
        chi2 += residual * residual / variance
        ndf += 1
    # endfor

    if ndf <= 0:
        return (math.nan, 0, math.nan)
    # endif

    return (chi2, ndf, chi2 / ndf)


def format_label_with_chi2(label: str, chi2_info: Tuple[float, int, float]) -> str:
    chi2, ndf, chi2ndf = chi2_info

    if ndf <= 0 or not np.isfinite(chi2ndf):
        return f"{label} ($\\chi^2$/ndf=N/A)"
    # endif

    return f"{label} ($\\chi^2$/ndf={chi2ndf:.2f})"


def ratio_to_hall_a(
    comparison_points: List[DataPoint],
    hall_a_points_input: List[DataPoint],
) -> List[DataPoint]:
    """
    Compute CLAS12 / Hall A ratios.
    """

    ratios: List[DataPoint] = []

    for b in comparison_points:
        a = interpolate_hall_a_point(
            phi_query=b.phi,
            hall_a_points_input=hall_a_points_input,
        )

        if not np.isfinite(a.sigma) or a.sigma == 0.0:
            continue
        # endif

        ratio = b.sigma / a.sigma

        ratio_err_high = math.sqrt(
            (b.err_high / a.sigma) ** 2
            + (b.sigma * a.err_low / (a.sigma * a.sigma)) ** 2
        )

        ratio_err_low = math.sqrt(
            (b.err_low / a.sigma) ** 2
            + (b.sigma * a.err_high / (a.sigma * a.sigma)) ** 2
        )

        ratios.append(
            DataPoint(
                phi=b.phi,
                sigma=ratio,
                err_low=ratio_err_low,
                err_high=ratio_err_high,
                stat=0.0,
                sys=0.0,
                sys_csv=0.0,
                sys_est=0.0,
                norm=0.0,
            )
        )
    # endfor

    ratios.sort(key=lambda p: p.phi)

    return ratios


# -----------------------------------------------------------------------------
# Plotting.
# -----------------------------------------------------------------------------

def plot_dataset(
    ax,
    points: List[DataPoint],
    label: str,
    legend_label: Optional[str] = None,
) -> None:
    if len(points) == 0:
        return
    # endif

    phi, sigma, err_low, err_high = points_to_arrays(points)
    style = SERIES_STYLE.get(label, dict(marker="o", linestyle="None"))

    if legend_label is None:
        legend_label = label
    # endif

    if label == PASS1_DISPLAY_LABEL:
        ax.errorbar(
            phi,
            sigma,
            yerr=np.vstack([err_low, err_high]),
            label=legend_label,
            markersize=6.0,
            linewidth=0.0,
            elinewidth=1.0,
            capsize=2.5,
            marker=style["marker"],
            linestyle=style["linestyle"],
            color=style["color"],
            markerfacecolor="none",
            markeredgecolor=style["color"],
            markeredgewidth=1.3,
        )
    else:
        ax.errorbar(
            phi,
            sigma,
            yerr=np.vstack([err_low, err_high]),
            label=legend_label,
            markersize=5.5,
            linewidth=0.0,
            elinewidth=1.0,
            capsize=2.5,
            **style,
        )
    # endif


def auto_ratio_ylim(ax, points_by_label: dict) -> None:
    yvals: List[float] = []
    yerr_lows: List[float] = []
    yerr_highs: List[float] = []

    for points in points_by_label.values():
        for p in points:
            if np.isfinite(p.sigma):
                yvals.append(p.sigma)
                yerr_lows.append(p.err_low if np.isfinite(p.err_low) else 0.0)
                yerr_highs.append(p.err_high if np.isfinite(p.err_high) else 0.0)
            # endif
        # endfor
    # endfor

    if len(yvals) == 0:
        ax.set_ylim(0.0, 2.0)
        return
    # endif

    lows = np.array(yvals) - np.array(yerr_lows)
    highs = np.array(yvals) + np.array(yerr_highs)

    ymin = float(np.nanmin(lows))
    ymax = float(np.nanmax(highs))

    span_low = abs(1.0 - ymin)
    span_high = abs(ymax - 1.0)
    span = max(span_low, span_high, 0.35)

    ax.set_ylim(1.0 - 1.10 * span, 1.0 + 1.10 * span)


def common_cross_section_ylim(
    points_groups: List[List[DataPoint]],
) -> Tuple[float, float]:
    values: List[float] = []
    lows: List[float] = []
    highs: List[float] = []

    for points in points_groups:
        for p in points:
            if not np.isfinite(p.sigma):
                continue
            # endif

            values.append(p.sigma)
            lows.append(p.sigma - (p.err_low if np.isfinite(p.err_low) else 0.0))
            highs.append(p.sigma + (p.err_high if np.isfinite(p.err_high) else 0.0))
        # endfor
    # endfor

    if len(values) == 0:
        return (0.0, 1.0)
    # endif

    ymin = min(0.0, float(np.nanmin(lows)))
    ymax = float(np.nanmax(highs))

    if not np.isfinite(ymax) or ymax <= ymin:
        return (0.0, 1.0)
    # endif

    return (ymin, 1.10 * ymax)


# -----------------------------------------------------------------------------
# Diagnostics.
# -----------------------------------------------------------------------------

def print_dataset_diagnostics(
    points_by_label: dict,
    clas12_scale: float,
    pass1_scale: float,
    pass1_norm_sys_frac: float,
) -> None:
    print()
    print("CLAS12 dataset extraction diagnostic:")
    print("  The script uses exact CSV column names, not column ordering.")
    print(f"  pass-2 clas12_scale = {clas12_scale:g}")
    print(f"  pass-1 pass1_scale  = {pass1_scale:g}")
    print("  pass-2 point-to-point sys source = Syst. err (point-to-point total)")
    print("  additional estimated pass-2 CLAS12 sys = none")

    print(f"  pass-1 normalization sys fraction = {pass1_norm_sys_frac:.6f}")

    for label, points in points_by_label.items():
        stat_values = [
            p.stat
            for p in points
            if np.isfinite(p.stat)
        ]

        sys_csv_values = [
            p.sys_csv
            for p in points
            if np.isfinite(p.sys_csv)
        ]

        sys_est_values = [
            p.sys_est
            for p in points
            if np.isfinite(p.sys_est)
        ]

        norm_values = [
            p.norm
            for p in points
            if np.isfinite(p.norm)
        ]

        sys_total_values = [
            p.sys
            for p in points
            if np.isfinite(p.sys)
        ]

        total_values = [
            symmetric_error(p)
            for p in points
            if np.isfinite(symmetric_error(p))
        ]

        mean_stat = float(np.mean(stat_values)) if len(stat_values) > 0 else math.nan
        mean_sys_csv = float(np.mean(sys_csv_values)) if len(sys_csv_values) > 0 else math.nan
        mean_sys_est = float(np.mean(sys_est_values)) if len(sys_est_values) > 0 else math.nan
        mean_norm = float(np.mean(norm_values)) if len(norm_values) > 0 else math.nan
        mean_sys_total = float(np.mean(sys_total_values)) if len(sys_total_values) > 0 else math.nan
        mean_total = float(np.mean(total_values)) if len(total_values) > 0 else math.nan

        print(f"  {label}")
        print(f"    points              = {len(points)}")
        print(f"    mean stat error     = {mean_stat:.6g} pb/GeV^4")
        print(f"    mean CSV sys        = {mean_sys_csv:.6g} pb/GeV^4")
        print(f"    mean estimated sys  = {mean_sys_est:.6g} pb/GeV^4")
        print(f"    mean norm sys       = {mean_norm:.6g} pb/GeV^4")
        print(f"    mean stored sys     = {mean_sys_total:.6g} pb/GeV^4")
        print(f"    mean total err      = {mean_total:.6g} pb/GeV^4")

        if len(points) > 0:
            first = points[0]
            print(
                f"    first point         = phi {first.phi:.3f} deg, "
                f"sigma {first.sigma:.6g} pb/GeV^4, "
                f"stat {first.stat:.6g}, "
                f"csv_sys {first.sys_csv:.6g}, "
                f"est_sys {first.sys_est:.6g}, "
                f"norm {first.norm:.6g}, "
                f"stored_sys {first.sys:.6g}, "
                f"total_err {symmetric_error(first):.6g}"
            )
        # endif
    # endfor


def print_model_correction_summary(model_context: ModelContext) -> None:
    if not model_context.enabled:
        print()
        print("Sp19 -> 10.604 GeV KM15 scaling: disabled")
        return
    # endif

    print()
    print("Sp19 -> 10.604 GeV KM15 scaling diagnostic:")
    print(f"  source beam energy = {SP19_INB_EBEAM_GEV:.6f} GeV")
    print(f"  target beam energy = {FA18_INB_EBEAM_GEV:.6f} GeV")
    print(f"  unique cached correction kinematic points = {len(model_context.correction_cache)}")

    if model_context.correction_count <= 0:
        print("  no corrections were evaluated")
        return
    # endif

    mean_corr = model_context.correction_sum / model_context.correction_count

    print(f"  correction count   = {model_context.correction_count}")
    print(f"  correction mean    = {mean_corr:.6g}")
    print(f"  correction range   = [{model_context.correction_min:.6g}, {model_context.correction_max:.6g}]")


# -----------------------------------------------------------------------------
# Plot builders.
# -----------------------------------------------------------------------------

def make_plot(
    pass2_selected: pd.DataFrame,
    pass1_selected: Optional[pd.DataFrame],
    output_dir: str,
    output_name: str,
    clas12_scale: float,
    pass1_scale: float,
    avg_period: str,
    pass1_norm_sys_frac: float,
) -> None:
    hall_a = hall_a_points()

    clas12_avg = compute_weighted_kinematic_averages(
        selected=pass2_selected,
        csv_period=avg_period,
    )

    points_by_label = {
        PASS2_COMBINED_DISPLAY_LABEL: pass2_points_for_period(
            selected=pass2_selected,
            csv_period=PASS2_COMBINED_CSV_PERIOD,
            clas12_scale=clas12_scale,
        )
    }

    if pass1_selected is not None:
        points_by_label[PASS1_DISPLAY_LABEL] = pass1_points(
            selected=pass1_selected,
            pass1_scale=pass1_scale,
            pass1_norm_sys_frac=pass1_norm_sys_frac,
        )
    # endif

    chi2_by_label = {
        label: chi2_ndf_to_hall_a(
            comparison_points=points,
            hall_a_points_input=hall_a,
        )
        for label, points in points_by_label.items()
    }

    ratios_by_label = {
        label: ratio_to_hall_a(
            comparison_points=points,
            hall_a_points_input=hall_a,
        )
        for label, points in points_by_label.items()
    }

    hall_a_tmin = compute_tmin_exact(HALL_A_XB, HALL_A_Q2)
    hall_a_t = hall_a_converted_t()
    hall_a_t_abs = abs(hall_a_t)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15.0, 5.8),
        constrained_layout=True,
        sharex=True,
    )

    fig.suptitle(
        (
            "Hall A / CLAS12 DVCS cross-section cross-check\n"
            r"pass-2 CLAS12: "
            rf"$\langle x_B\rangle={clas12_avg.xB:.3f}$, "
            rf"$\langle Q^2\rangle={clas12_avg.Q2:.3f}~{{\rm GeV}}^2$, "
            rf"$\langle |t|\rangle={clas12_avg.t_abs:.3f}~{{\rm GeV}}^2$"
            "\n"
            r"Hall A: "
            rf"$\langle x_B\rangle={HALL_A_XB:.3f}$, "
            rf"$\langle Q^2\rangle={HALL_A_Q2:.3f}~{{\rm GeV}}^2$, "
            rf"$\langle t'\rangle={HALL_A_TPRIME:.3f}~{{\rm GeV}}^2$, "
            rf"$t_{{\min}}={hall_a_tmin:.3f}~{{\rm GeV}}^2$, "
            rf"$|t|_{{\rm equiv.}}={hall_a_t_abs:.3f}~{{\rm GeV}}^2$"
        ),
        fontsize=13,
    )

    left = axes[0]
    right = axes[1]

    plot_dataset(
        ax=left,
        points=hall_a,
        label="Hall A",
        legend_label="Hall A",
    )

    for label in TOP_COMPARISON_SERIES:
        if label not in points_by_label:
            continue
        # endif

        plot_dataset(
            ax=left,
            points=points_by_label.get(label, []),
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    left.set_xlabel(r"$\phi$ (deg)")
    left.set_ylabel(r"$d^4\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ (pb/GeV$^4$)")
    left.set_title("Cross sections: integrated pass-2 10.6 GeV and pass-1")
    left.grid(True, alpha=0.25)
    left.legend(fontsize=8, frameon=True)

    top_ratios = {}

    for label in TOP_COMPARISON_SERIES:
        if label not in ratios_by_label:
            continue
        # endif

        ratios = ratios_by_label.get(label, [])
        top_ratios[label] = ratios

        plot_dataset(
            ax=right,
            points=ratios,
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    right.axhline(
        1.0,
        color="0.35",
        linewidth=1.0,
        linestyle="--",
        zorder=0,
    )

    right.set_xlabel(r"$\phi$ (deg)")
    right.set_ylabel(r"CLAS12 / Hall A")
    right.set_title("Ratio to Hall A: integrated pass-2 10.6 GeV and pass-1")
    right.grid(True, alpha=0.25)
    right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(right, top_ratios)

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, output_name)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote: {output_path}")

    print()
    print("pass-2 CLAS12 weighted average kinematics used in title:")
    print(f"  weighting period = {clas12_avg.weight_period}")
    print(f"  points used      = {clas12_avg.n_points}")
    print(f"  <xB>             = {clas12_avg.xB:.6f}")
    print(f"  <Q2>             = {clas12_avg.Q2:.6f} GeV^2")
    print(f"  <|t|>            = {clas12_avg.t_abs:.6f} GeV^2")

    print()
    print("Hall A t' -> t diagnostic:")
    print(f"  Hall A <xB>      = {HALL_A_XB:.6f}")
    print(f"  Hall A <Q2>      = {HALL_A_Q2:.6f} GeV^2")
    print(f"  Hall A <t'>      = {HALL_A_TPRIME:.6f} GeV^2")
    print(f"  computed tmin    = {hall_a_tmin:.6f} GeV^2")
    print(f"  converted t      = tmin + t' = {hall_a_t:.6f} GeV^2")
    print(f"  converted |t|    = {hall_a_t_abs:.6f} GeV^2")
    print(f"  selected CLAS12 |t| bin = [{TARGET_T_ABS_MIN:.3f}, {TARGET_T_ABS_MAX:.3f}] GeV^2")

    print()
    print("Chi2/ndf summary against Hall A interpolation:")
    for label in TOP_COMPARISON_SERIES:
        if label not in chi2_by_label:
            continue
        # endif

        chi2, ndf, chi2ndf = chi2_by_label[label]

        if ndf > 0 and np.isfinite(chi2ndf):
            print(f"  {label:18s}: chi2 = {chi2:.4f}, ndf = {ndf:d}, chi2/ndf = {chi2ndf:.4f}")
        else:
            print(f"  {label:18s}: chi2/ndf = N/A")
        # endif
    # endfor

    print_dataset_diagnostics(
        points_by_label=points_by_label,
        clas12_scale=clas12_scale,
        pass1_scale=pass1_scale,
        pass1_norm_sys_frac=pass1_norm_sys_frac,
    )


def make_fa18_sp19_km15_plot(
    pass2_selected: pd.DataFrame,
    output_dir: str,
    output_name: str,
    clas12_scale: float,
    model_context: ModelContext,
) -> None:
    hall_a = hall_a_points()

    fa18_points = pass2_points_for_period(
        selected=pass2_selected,
        csv_period="Fa18 Inb",
        clas12_scale=clas12_scale,
    )

    sp19_scaled_points = pass2_points_for_period(
        selected=pass2_selected,
        csv_period="Sp19 Inb",
        clas12_scale=clas12_scale,
        model_context=model_context,
        apply_sp19_to_fa18_km15_scaling=model_context.enabled,
    )

    points_by_label = {
        FA18_INB_DISPLAY_LABEL: fa18_points,
        SP19_KM15_DISPLAY_LABEL: sp19_scaled_points,
    }

    chi2_by_label = {
        label: chi2_ndf_to_hall_a(
            comparison_points=points,
            hall_a_points_input=hall_a,
        )
        for label, points in points_by_label.items()
    }

    ratios_by_label = {
        label: ratio_to_hall_a(
            comparison_points=points,
            hall_a_points_input=hall_a,
        )
        for label, points in points_by_label.items()
    }

    fa18_avg = compute_weighted_kinematic_averages(
        selected=pass2_selected,
        csv_period="Fa18 Inb",
    )

    sp19_avg = compute_weighted_kinematic_averages(
        selected=pass2_selected,
        csv_period="Sp19 Inb",
    )

    hall_a_tmin = compute_tmin_exact(HALL_A_XB, HALL_A_Q2)
    hall_a_t = hall_a_converted_t()
    hall_a_t_abs = abs(hall_a_t)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15.5, 5.8),
        constrained_layout=True,
        sharex=True,
    )

    fig.suptitle(
        (
            "Hall A / CLAS12 Fa18-Sp19 energy-scaled cross-section cross-check\n"
            rf"Fa18 Inb: $\langle x_B\rangle={fa18_avg.xB:.3f}$, "
            rf"$\langle Q^2\rangle={fa18_avg.Q2:.3f}~{{\rm GeV}}^2$, "
            rf"$\langle |t|\rangle={fa18_avg.t_abs:.3f}~{{\rm GeV}}^2$"
            "\n"
            rf"Sp19 Inb scaled with KM15: $E={SP19_INB_EBEAM_GEV:.4f}\rightarrow {FA18_INB_EBEAM_GEV:.4f}~{{\rm GeV}}$, "
            rf"$\langle x_B\rangle={sp19_avg.xB:.3f}$, "
            rf"$\langle Q^2\rangle={sp19_avg.Q2:.3f}~{{\rm GeV}}^2$, "
            rf"$\langle |t|\rangle={sp19_avg.t_abs:.3f}~{{\rm GeV}}^2$"
            "\n"
            r"Hall A: "
            rf"$\langle x_B\rangle={HALL_A_XB:.3f}$, "
            rf"$\langle Q^2\rangle={HALL_A_Q2:.3f}~{{\rm GeV}}^2$, "
            rf"$|t|_{{\rm equiv.}}={hall_a_t_abs:.3f}~{{\rm GeV}}^2$"
        ),
        fontsize=12,
    )

    left = axes[0]
    right = axes[1]

    plot_dataset(
        ax=left,
        points=hall_a,
        label="Hall A",
        legend_label="Hall A",
    )

    for label in FA18_SP19_COMPARISON_SERIES:
        plot_dataset(
            ax=left,
            points=points_by_label.get(label, []),
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    left.set_xlabel(r"$\phi$ (deg)")
    left.set_ylabel(r"$d^4\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ (pb/GeV$^4$)")
    left.set_title("Cross sections: Fa18 and Sp19 scaled to 10.604 GeV")
    left.grid(True, alpha=0.25)
    left.legend(fontsize=8, frameon=True)

    ylo, yhi = common_cross_section_ylim([hall_a, fa18_points, sp19_scaled_points])
    left.set_ylim(ylo, yhi)

    for label in FA18_SP19_COMPARISON_SERIES:
        plot_dataset(
            ax=right,
            points=ratios_by_label.get(label, []),
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    right.axhline(
        1.0,
        color="0.35",
        linewidth=1.0,
        linestyle="--",
        zorder=0,
    )

    right.set_xlabel(r"$\phi$ (deg)")
    right.set_ylabel(r"CLAS12 / Hall A")
    right.set_title("Ratio to Hall A: Fa18 and scaled Sp19")
    right.grid(True, alpha=0.25)
    right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(right, ratios_by_label)

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, output_name)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote: {output_path}")

    print()
    print("Fa18 vs Sp19 KM15-scaled Hall A comparison diagnostic:")
    print(f"  Fa18 points                = {len(fa18_points)}")
    print(f"  Sp19 scaled points          = {len(sp19_scaled_points)}")
    print(f"  Hall A points               = {len(hall_a)}")
    print(f"  Hall A computed tmin        = {hall_a_tmin:.6f} GeV^2")
    print(f"  Hall A converted t          = {hall_a_t:.6f} GeV^2")
    print(f"  Hall A converted |t|        = {hall_a_t_abs:.6f} GeV^2")

    print()
    print("Fa18 vs Sp19 KM15-scaled chi2/ndf summary against Hall A interpolation:")
    for label in FA18_SP19_COMPARISON_SERIES:
        chi2, ndf, chi2ndf = chi2_by_label.get(label, (math.nan, 0, math.nan))

        if ndf > 0 and np.isfinite(chi2ndf):
            print(f"  {label:30s}: chi2 = {chi2:.4f}, ndf = {ndf:d}, chi2/ndf = {chi2ndf:.4f}")
        else:
            print(f"  {label:30s}: chi2/ndf = N/A")
        # endif
    # endfor


# -----------------------------------------------------------------------------
# Argument parsing and main.
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make a Hall A / CLAS12 DVCS cross-section cross-check plot."
    )

    parser.add_argument(
        "csv_file",
        help="Input final pass-2 DVCS analysis CSV.",
    )

    parser.add_argument(
        "--pass1-csv",
        default=DEFAULT_PASS1_CSV,
        help=(
            "Input pass-1 Hall B / CLAS12 CSV. "
            f"Default: {DEFAULT_PASS1_CSV}"
        ),
    )

    parser.add_argument(
        "--no-pass1",
        action="store_true",
        help="Disable pass-1 Fa18 overlay.",
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for output plots. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--output-name",
        default=DEFAULT_OUTPUT_NAME,
        help=f"Output PNG filename for the integrated 1x2 plot. Default: {DEFAULT_OUTPUT_NAME}",
    )

    parser.add_argument(
        "--fa18-sp19-output-name",
        default=DEFAULT_FA18_SP19_OUTPUT_NAME,
        help=(
            "Output PNG filename for the Fa18/Sp19 KM15-scaled 1x2 plot. "
            f"Default: {DEFAULT_FA18_SP19_OUTPUT_NAME}"
        ),
    )

    parser.add_argument(
        "--no-fa18-sp19-plot",
        action="store_true",
        help="Disable the additional Fa18/Sp19 KM15-scaled Hall A comparison plot.",
    )

    parser.add_argument(
        "--bin-name",
        type=int,
        default=DEFAULT_BIN_NAME,
        help=f"CLAS12 Bin Name to select. Default: {DEFAULT_BIN_NAME}",
    )

    parser.add_argument(
        "--clas12-scale",
        type=float,
        default=DEFAULT_CLAS12_SCALE,
        help=(
            "Scale factor applied to pass-2 CLAS12 cross sections and uncertainties. "
            f"Default: {DEFAULT_CLAS12_SCALE:g}, useful if the CSV is effectively "
            "in nb/GeV^4 and Hall A is in pb/GeV^4."
        ),
    )

    parser.add_argument(
        "--pass1-scale",
        type=float,
        default=DEFAULT_PASS1_SCALE,
        help=(
            "Scale factor applied to pass-1 CLAS12 cross sections and uncertainties. "
            f"Default: {DEFAULT_PASS1_SCALE:g}."
        ),
    )

    parser.add_argument(
        "--avg-period",
        default="10.6 GeV",
        help=(
            "Pass-2 period used for weighted average CLAS12 kinematics in the plot title. "
            "Default: '10.6 GeV'."
        ),
    )

    parser.add_argument(
        "--pass1-norm-sys-frac",
        type=float,
        default=DEFAULT_PASS1_NORM_SYS_FRAC,
        help=(
            "Additional fractional normalization systematic applied to pass-1 Fa18. "
            f"Default: {DEFAULT_PASS1_NORM_SYS_FRAC:.2f}."
        ),
    )

    parser.add_argument(
        "--no-sp19-km15-energy-scaling",
        action="store_true",
        help=(
            "Disable the KM15 beam-energy scaling of Sp19 Inb from "
            "10.1998 GeV to 10.604 GeV in the additional Fa18/Sp19 comparison plot."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.clas12_scale <= 0.0:
        raise ValueError("--clas12-scale must be positive.")
    # endif

    if args.pass1_scale <= 0.0:
        raise ValueError("--pass1-scale must be positive.")
    # endif

    if args.pass1_norm_sys_frac < 0.0:
        raise ValueError("--pass1-norm-sys-frac must be non-negative.")
    # endif

    pass2_required = [
        "Bin Name",
        "xBmin",
        "xBmax",
        "Q2min",
        "Q2max",
        "t_abs_min",
        "t_abs_max",
        "phimin",
        "phimax",
    ]

    pass2_required += [
        pass2_cross_section_column(csv_period)
        for csv_period in ALL_PASS2_CSV_PERIODS
    ]

    pass2_required.append(pass2_point_to_point_column())
    pass2_required += sorted({
        pass2_scale_systematic_column(csv_period)
        for csv_period in ALL_PASS2_CSV_PERIODS
    })

    for period in ["10.6 GeV", "Fa18 Inb", "Sp19 Inb"]:
        pass2_required += [
            avg_column("xBavg", period),
            avg_column("Q2avg", period),
            avg_column("t_abs_avg", period),
            average_phi_column(period),
        ]
    # endfor

    pass2_df = pd.read_csv(args.csv_file, low_memory=False)
    require_columns(pass2_df, pass2_required, context="pass-2")

    print("Pass-2 correlated scale-systematic mapping:")
    for csv_period in ALL_PASS2_CSV_PERIODS:
        scale_period = pass2_scale_systematic_period(csv_period)
        print(
            f"  {csv_period:10s} -> {scale_period:10s}: "
            f"{pass2_scale_systematic_column(csv_period)}"
        )
    # endfor

    pass2_selected = select_overlap_bin(
        df=pass2_df,
        bin_name=args.bin_name,
        context="pass-2",
    )

    pass1_selected: Optional[pd.DataFrame] = None

    if not args.no_pass1:
        if not os.path.exists(args.pass1_csv):
            print()
            print("WARNING: pass-1 CSV was requested but not found:")
            print(f"  {args.pass1_csv}")
            print("  Continuing without pass-1 Fa18.")
        else:
            pass1_required = [
                "Bin Name",
                "xBmin",
                "xBmax",
                "Q2min",
                "Q2max",
                "t_abs_min",
                "t_abs_max",
                "phimin",
                "phimax",
                "phiavg",
                "cross sections, ep->epg, exp",
                "cross sections, ep->epg, exp, stat. unc.",
                "cross sections, ep->epg, exp, syst. unc. (up)",
                "cross sections, ep->epg, exp, syst. unc. (down)",
            ]

            if str(args.pass1_csv).lower().endswith((".txt", ".dat")):
                pass1_df = load_published_pass1_dataframe(args.pass1_csv)
            else:
                pass1_df = pd.read_csv(args.pass1_csv, low_memory=False)
            require_columns(pass1_df, pass1_required, context="pass-1")

            pass1_selected = select_overlap_bin(
                df=pass1_df,
                bin_name=args.bin_name,
                context="pass-1",
            )
        # endif
    # endif

    print("Selected Hall A overlap bin from pass-2 CSV:")
    print(f"  rows: {len(pass2_selected)}")
    print(f"  Bin Name: {args.bin_name}")
    print(f"  nominal xB:  [{TARGET_XB_MIN}, {TARGET_XB_MAX}]")
    print(f"  nominal Q2:  [{TARGET_Q2_MIN}, {TARGET_Q2_MAX}] GeV^2")
    print(f"  nominal |t|: [{TARGET_T_ABS_MIN}, {TARGET_T_ABS_MAX}] GeV^2")
    print(f"  phi bins present: {len(pass2_selected)}")

    if pass1_selected is not None:
        print()
        print("Selected Hall A overlap bin from pass-1 CSV:")
        print(f"  rows: {len(pass1_selected)}")
        print(f"  Bin Name: {args.bin_name}")
        print(f"  phi bins present: {len(pass1_selected)}")
    # endif

    print("  pass-2 point-to-point systematic: read from CSV column")
    print("    Syst. err (point-to-point total)")
    print("  additional estimated pass-2 bin-to-bin systematic: none")

    if args.no_pass1:
        print("  pass-1 Fa18 overlay: disabled")
    else:
        print(f"  pass-1 CSV: {args.pass1_csv}")
        print(f"  pass-1 normalization systematic: {args.pass1_norm_sys_frac:.6f}")
    # endif

    make_plot(
        pass2_selected=pass2_selected,
        pass1_selected=pass1_selected,
        output_dir=args.output_dir,
        output_name=args.output_name,
        clas12_scale=args.clas12_scale,
        pass1_scale=args.pass1_scale,
        avg_period=args.avg_period,
        pass1_norm_sys_frac=args.pass1_norm_sys_frac,
    )

    model_context = make_model_context(
        enabled=not args.no_sp19_km15_energy_scaling,
    )

    if args.no_fa18_sp19_plot:
        print()
        print("Additional Fa18/Sp19 KM15-scaled Hall A comparison plot: disabled")
    else:
        make_fa18_sp19_km15_plot(
            pass2_selected=pass2_selected,
            output_dir=args.output_dir,
            output_name=args.fa18_sp19_output_name,
            clas12_scale=args.clas12_scale,
            model_context=model_context,
        )
    # endif

    print_model_correction_summary(model_context)


if __name__ == "__main__":
    main()
# endif