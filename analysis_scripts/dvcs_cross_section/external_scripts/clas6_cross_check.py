#!/usr/bin/env python3
"""
clas6_cross_check.py

Single-bin CLAS6 / CLAS12 DVCS cross-section cross-check.

This script compares the current pass-2 CLAS12 / Hall B result against:

  - previous pass-1 CLAS12 / Hall B Fa18 result
  - CLAS6 unpolarized cross sections from Gepard

For now this script is intentionally hard-coded to the first comparison bin
from the pass-1 analysis note:

  CLAS6:
    <xB>  = 0.185
    <Q2>  = 1.630 GeV^2
    <|t|> = 0.200 GeV^2

  CLAS12:
    Bin Name = 58
    xB  in [0.155, 0.204]
    Q2  in [1.456, 1.912] GeV^2
    |t| in [0.15, 0.25] GeV^2

The CLAS6 data are loaded from Gepard, by default dataset id 98:

  [98] CLAS 2640 XUU 1504.02009 CLAS data base E145M1

This corresponds to the unpolarized cross section measurement from:

  H.S. Jo et al., Phys. Rev. Lett. 115, 212003 (2015)

Important fixes relative to the first draft
-------------------------------------------

1. Zero-placeholder Gepard points are removed.

   Gepard dataset 98 includes entries with:

     val = 0, errstat = 0, errsyst = 0

   at phi values where no physical data point should be plotted. These are
   placeholders and badly distort the plot and interpolation if kept.

2. The CLAS12 bin is hard-coded.

   The previous pass-1 comparison did not simply choose the nearest CLAS12 bin
   independently for every CLAS6 |t| point. For the first comparison, the correct
   CLAS12 bin is Bin Name 58.

Plot layout
-----------

The script makes one 2x2 figure:

  top-left:     pass-2 10.6 GeV, pass-1 Fa18, and CLAS6 cross sections vs phi
  top-right:    pass-2 10.6 GeV / CLAS6 ratio vs phi
  bottom-left:  individual pass-2 run periods and CLAS6 cross sections vs phi
  bottom-right: individual pass-2 run periods / CLAS6 ratios vs phi

The pass-1 result is shown only in the top-left panel.

Uncertainty treatment
---------------------

CLAS6:
  total uncertainty is sqrt(errstat^2 + errsyst^2) when available. If no
  stat/sys decomposition exists, the Gepard err value is used.

pass-2 CLAS12:
  CSV tuple is read as:

    (value, stat, sys_csv)

  By default, an additional 15% bin-to-bin systematic is added:

    sys_est = 0.15 * sigma

  and the plotted uncertainty is:

    err_total = sqrt(stat^2 + sys_csv^2 + sys_est^2)

pass-1 CLAS12:
  CSV columns are read as:

    cross sections, ep->epg, exp
    cross sections, ep->epg, exp, stat. unc.
    cross sections, ep->epg, exp, syst. unc. (up)
    cross sections, ep->epg, exp, syst. unc. (down)

  A 31% normalization uncertainty is added in quadrature:

    norm = 0.31 * sigma

  so:

    err_up   = sqrt(stat^2 + sys_up^2 + norm^2)
    err_down = sqrt(stat^2 + sys_down^2 + norm^2)

Units
-----

By default, pass-2, pass-1, and CLAS6 are all left in their native units,
expected to be nb/GeV^4 for this CLAS6 comparison.

Output
------

  output/clas6_cross_check/clas6_cross_check_26_p2bin58.png
"""

import argparse
import math
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Hard-coded comparison definition.
# -----------------------------------------------------------------------------

COMPARISON_LABEL = "26"

CLAS6_TARGET_XB = 0.185
CLAS6_TARGET_Q2 = 1.630
CLAS6_TARGET_T_ABS = 0.200

PASS2_BIN_NAME = 58

PASS2_EXPECTED_XB_MIN = 0.155
PASS2_EXPECTED_XB_MAX = 0.204

PASS2_EXPECTED_Q2_MIN = 1.456
PASS2_EXPECTED_Q2_MAX = 1.912

PASS2_EXPECTED_T_ABS_MIN = 0.15
PASS2_EXPECTED_T_ABS_MAX = 0.25


# -----------------------------------------------------------------------------
# Default file locations.
# -----------------------------------------------------------------------------

DEFAULT_PASS1_CSV = (
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "dvcs_cross_section/imports/all_bin_v3.csv"
)


# -----------------------------------------------------------------------------
# Defaults.
# -----------------------------------------------------------------------------

DEFAULT_PASS2_SCALE = 1.0
DEFAULT_PASS1_SCALE = 1.0
DEFAULT_CLAS6_SCALE = 1.0

DEFAULT_PASS2_BIN_TO_BIN_SYS_FRAC = 0.15
DEFAULT_PASS1_NORM_SYS_FRAC = 0.31

DEFAULT_CLAS6_DATASET_ID = 98

DEFAULT_OUTPUT_DIR = "output/clas6_cross_check"


# -----------------------------------------------------------------------------
# Plot labels and styles.
# -----------------------------------------------------------------------------

PASS2_COMBINED_DISPLAY_LABEL = "pass-2 10.6 GeV"
PASS2_COMBINED_CSV_PERIOD = "10.6 GeV"

PASS1_DISPLAY_LABEL = "pass-1 Fa18"
CLAS6_DISPLAY_LABEL = "CLAS6"

PASS2_PERIOD_DISPLAY_TO_CSV_PERIOD = {
    "Sp18 Inb": "Sp18 Inb",
    "Sp18 Out": "Sp18 Out",
    "Fa18 Inb": "Fa18 Inb",
    "Fa18 Out": "Fa18 Out",
}

TOP_CROSS_SECTION_SERIES = [
    PASS2_COMBINED_DISPLAY_LABEL,
    PASS1_DISPLAY_LABEL,
]

TOP_RATIO_SERIES = [
    PASS2_COMBINED_DISPLAY_LABEL,
]

BOTTOM_CROSS_SECTION_SERIES = [
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
]

BOTTOM_RATIO_SERIES = [
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
]

ALL_PASS2_CSV_PERIODS = [
    PASS2_COMBINED_CSV_PERIOD,
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
]

SERIES_STYLE = {
    CLAS6_DISPLAY_LABEL: dict(marker="s", linestyle="None", color="black"),
    PASS2_COMBINED_DISPLAY_LABEL: dict(marker="o", linestyle="None", color="tab:red"),
    PASS1_DISPLAY_LABEL: dict(marker="o", linestyle="None", color="tab:cyan"),
    "Sp18 Inb": dict(marker="D", linestyle="None", color="tab:green"),
    "Sp18 Out": dict(marker="P", linestyle="None", color="tab:purple"),
    "Fa18 Inb": dict(marker="^", linestyle="None", color="tab:blue"),
    "Fa18 Out": dict(marker="v", linestyle="None", color="tab:orange"),
}

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
class KinematicPoint:
    xB: float
    Q2: float
    t_abs: float


@dataclass
class KinematicAverages:
    xB: float
    Q2: float
    t_abs: float
    weight_period: str
    n_points: int


# -----------------------------------------------------------------------------
# Generic helpers.
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


def finite_or_zero(value: float) -> float:
    if np.isfinite(value):
        return float(value)
    # endif

    return 0.0


def is_zero_placeholder(*values: float, tolerance: float = 1.0e-15) -> bool:
    finite_values = []

    for value in values:
        try:
            value_float = float(value)
        except Exception:
            continue
        # endtry

        if np.isfinite(value_float):
            finite_values.append(abs(value_float))
        # endif
    # endfor

    if len(finite_values) == 0:
        return True
    # endif

    return all(value <= tolerance for value in finite_values)


def symmetric_error(point: DataPoint) -> float:
    return 0.5 * (point.err_low + point.err_high)


def points_to_arrays(points: List[DataPoint]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    phi = np.array([p.phi for p in points], dtype=float)
    sigma = np.array([p.sigma for p in points], dtype=float)
    err_low = np.array([p.err_low for p in points], dtype=float)
    err_high = np.array([p.err_high for p in points], dtype=float)

    return phi, sigma, err_low, err_high


def require_columns(df: pd.DataFrame, columns: Iterable[str], context: str) -> None:
    missing = [c for c in columns if c not in df.columns]

    if missing:
        raise KeyError(
            f"The {context} CSV is missing required columns:\n  "
            + "\n  ".join(missing)
        )
    # endif


def value_close(a: float, b: float, tolerance: float = 1.0e-6) -> bool:
    return abs(float(a) - float(b)) <= tolerance


# -----------------------------------------------------------------------------
# CSV column helpers.
# -----------------------------------------------------------------------------

def pass2_cross_section_column(csv_period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {csv_period}, unpol"


def average_phi_column(csv_period: str) -> str:
    return f"phiavg, {csv_period}"


def avg_column(quantity: str, csv_period: str) -> str:
    return f"{quantity}, {csv_period}"


# -----------------------------------------------------------------------------
# Gepard / CLAS6.
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


def get_attr_any(obj, names: List[str], default=math.nan):
    for name in names:
        if hasattr(obj, name):
            return getattr(obj, name)
        # endif
    # endfor

    return default


def convert_gepard_phi_to_degrees(pt, phi_convention: str) -> float:
    """
    Convert Gepard point phi to degrees.

    Gepard stores available datasets internally in BMK convention. For datasets
    originally given in Trento, converting back to the published-style angle is:

      phi_Trento = pi - phi_BMK

    modulo 2pi.

    The conversion is applied only when:
      --clas6-phi-convention trento
    and the stored point reports frame='Trento'.
    """

    phi_raw = get_attr_any(pt, ["phi"])
    phi = float(phi_raw)

    is_radians = abs(phi) <= 2.0 * math.pi * 1.5
    frame = str(get_attr_any(pt, ["frame"], default="")).lower()

    if phi_convention == "trento" and is_radians and frame == "trento":
        phi = (math.pi - phi) % (2.0 * math.pi)
    # endif

    if is_radians:
        phi_deg = math.degrees(phi)
    else:
        phi_deg = phi
    # endif

    return phi_deg % 360.0


def clas6_point_from_gepard(
    pt,
    clas6_scale: float,
    phi_convention: str,
) -> Optional[Tuple[KinematicPoint, DataPoint]]:
    xB = float(get_attr_any(pt, ["xB", "xb"]))
    Q2 = float(get_attr_any(pt, ["Q2", "q2"]))
    t = float(get_attr_any(pt, ["t"]))
    val = float(get_attr_any(pt, ["val", "value"]))

    if not (np.isfinite(xB) and np.isfinite(Q2) and np.isfinite(t) and np.isfinite(val)):
        return None
    # endif

    stat = get_attr_any(pt, ["errstat", "err_stat", "stat"], default=math.nan)
    sys = get_attr_any(pt, ["errsyst", "errsys", "err_syst", "sys"], default=math.nan)
    err = get_attr_any(pt, ["err"], default=math.nan)

    stat = float(stat) if np.isfinite(stat) else math.nan
    sys = float(sys) if np.isfinite(sys) else math.nan
    err = float(err) if np.isfinite(err) else math.nan

    if is_zero_placeholder(val, stat, sys, err):
        return None
    # endif

    if np.isfinite(stat) or np.isfinite(sys):
        stat = finite_or_zero(stat)
        sys = finite_or_zero(sys)
        total = math.sqrt(stat * stat + sys * sys)
    elif np.isfinite(err):
        stat = float(err)
        sys = 0.0
        total = float(err)
    else:
        stat = 0.0
        sys = 0.0
        total = 0.0
    # endif

    val *= clas6_scale
    stat *= clas6_scale
    sys *= clas6_scale
    total *= clas6_scale

    phi_deg = convert_gepard_phi_to_degrees(pt, phi_convention=phi_convention)

    kin = KinematicPoint(
        xB=xB,
        Q2=Q2,
        t_abs=abs(t),
    )

    data = DataPoint(
        phi=phi_deg,
        sigma=val,
        err_low=total,
        err_high=total,
        stat=stat,
        sys=sys,
        sys_csv=sys,
        sys_est=0.0,
        norm=0.0,
    )

    return kin, data


def load_clas6_points_for_target(
    dataset_id: int,
    clas6_scale: float,
    phi_convention: str,
    target_xB: float,
    target_Q2: float,
    target_t_abs: float,
    tolerance_xB: float,
    tolerance_Q2: float,
    tolerance_t: float,
) -> Tuple[KinematicPoint, List[DataPoint], int, int]:
    g = import_gepard()

    if dataset_id not in g.dset:
        raise RuntimeError(
            f"Gepard dataset id {dataset_id} is not available in g.dset."
        )
    # endif

    dataset = g.dset[dataset_id]

    selected_points: List[Tuple[KinematicPoint, DataPoint]] = []
    n_raw_matching_kinematics = 0
    n_zero_placeholders_removed = 0

    for pt in dataset:
        xB = float(get_attr_any(pt, ["xB", "xb"]))
        Q2 = float(get_attr_any(pt, ["Q2", "q2"]))
        t = float(get_attr_any(pt, ["t"]))

        if not (
            abs(xB - target_xB) <= tolerance_xB
            and abs(Q2 - target_Q2) <= tolerance_Q2
            and abs(abs(t) - target_t_abs) <= tolerance_t
        ):
            continue
        # endif

        n_raw_matching_kinematics += 1

        val = float(get_attr_any(pt, ["val", "value"], default=math.nan))
        stat = get_attr_any(pt, ["errstat", "err_stat", "stat"], default=math.nan)
        sys = get_attr_any(pt, ["errsyst", "errsys", "err_syst", "sys"], default=math.nan)
        err = get_attr_any(pt, ["err"], default=math.nan)

        if is_zero_placeholder(val, stat, sys, err):
            n_zero_placeholders_removed += 1
            continue
        # endif

        converted = clas6_point_from_gepard(
            pt=pt,
            clas6_scale=clas6_scale,
            phi_convention=phi_convention,
        )

        if converted is None:
            continue
        # endif

        selected_points.append(converted)
    # endfor

    if len(selected_points) == 0:
        raise RuntimeError(
            "No nonzero CLAS6 points were selected from Gepard for target:\n"
            f"  xB={target_xB}, Q2={target_Q2}, |t|={target_t_abs}\n"
            f"  dataset id={dataset_id}"
        )
    # endif

    kin = KinematicPoint(
        xB=float(np.mean([item[0].xB for item in selected_points])),
        Q2=float(np.mean([item[0].Q2 for item in selected_points])),
        t_abs=float(np.mean([item[0].t_abs for item in selected_points])),
    )

    points = [item[1] for item in selected_points]
    points.sort(key=lambda p: p.phi)

    return kin, points, n_raw_matching_kinematics, n_zero_placeholders_removed


# -----------------------------------------------------------------------------
# CLAS12 pass-2 selection and averages.
# -----------------------------------------------------------------------------

def select_bin_by_name(df: pd.DataFrame, bin_name: int) -> pd.DataFrame:
    selected = df.loc[pd.to_numeric(df["Bin Name"], errors="coerce") == int(bin_name)].copy()
    selected = selected.sort_values(["phimin", "phimax"]).copy()

    if len(selected) == 0:
        raise RuntimeError(f"Could not find Bin Name {bin_name} in CSV.")
    # endif

    return selected


def validate_pass2_bin_edges(selected: pd.DataFrame, tolerance: float = 1.0e-6) -> None:
    first = selected.iloc[0]

    xBmin = parse_scalar_from_cell(first["xBmin"])
    xBmax = parse_scalar_from_cell(first["xBmax"])
    Q2min = parse_scalar_from_cell(first["Q2min"])
    Q2max = parse_scalar_from_cell(first["Q2max"])
    tmin = parse_scalar_from_cell(first["t_abs_min"])
    tmax = parse_scalar_from_cell(first["t_abs_max"])

    checks = [
        ("xBmin", xBmin, PASS2_EXPECTED_XB_MIN),
        ("xBmax", xBmax, PASS2_EXPECTED_XB_MAX),
        ("Q2min", Q2min, PASS2_EXPECTED_Q2_MIN),
        ("Q2max", Q2max, PASS2_EXPECTED_Q2_MAX),
        ("t_abs_min", tmin, PASS2_EXPECTED_T_ABS_MIN),
        ("t_abs_max", tmax, PASS2_EXPECTED_T_ABS_MAX),
    ]

    failed = []

    for name, actual, expected in checks:
        if not value_close(actual, expected, tolerance=tolerance):
            failed.append((name, actual, expected))
        # endif
    # endfor

    if failed:
        message = [
            f"Selected Bin Name {PASS2_BIN_NAME} does not have the expected bin edges."
        ]

        for name, actual, expected in failed:
            message.append(f"  {name}: actual={actual}, expected={expected}")
        # endfor

        raise RuntimeError("\n".join(message))
    # endif


def compute_weighted_kinematic_averages(
    selected: pd.DataFrame,
    csv_period: str = "10.6 GeV",
    fallback_periods: Optional[List[str]] = None,
) -> KinematicAverages:
    if fallback_periods is None:
        fallback_periods = [
            "10.6 GeV",
            "Fa18 Inb",
            "Fa18 Out",
            "Sp18 Inb",
            "Sp18 Out",
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
        return KinematicAverages(
            xB=fallback_midpoint("xBmin", "xBmax"),
            Q2=fallback_midpoint("Q2min", "Q2max"),
            t_abs=fallback_midpoint("t_abs_min", "t_abs_max"),
            weight_period="bin midpoint fallback",
            n_points=0,
        )
    # endif

    weighted_sums = {"xB": 0.0, "Q2": 0.0, "t_abs": 0.0}
    unweighted_values = {"xB": [], "Q2": [], "t_abs": []}

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


# -----------------------------------------------------------------------------
# CLAS12 point extraction.
# -----------------------------------------------------------------------------

def pass2_points_for_period(
    selected: pd.DataFrame,
    csv_period: str,
    pass2_scale: float,
    include_pass2_estimated_sys: bool,
    pass2_bin_to_bin_sys_frac: float,
) -> List[DataPoint]:
    xs_col = pass2_cross_section_column(csv_period)
    phi_avg_col = average_phi_column(csv_period)

    if xs_col not in selected.columns:
        return []
    # endif

    points: List[DataPoint] = []

    for _, row in selected.sort_values("phimin").iterrows():
        sigma, stat, sys_csv = parse_tuple3(row[xs_col])

        if not np.isfinite(sigma):
            continue
        # endif

        stat = finite_or_zero(stat)
        sys_csv = finite_or_zero(sys_csv)

        if is_zero_placeholder(sigma, stat, sys_csv):
            continue
        # endif

        sigma *= pass2_scale
        stat *= pass2_scale
        sys_csv *= pass2_scale

        if include_pass2_estimated_sys:
            sys_est = pass2_bin_to_bin_sys_frac * abs(sigma)
        else:
            sys_est = 0.0
        # endif

        sys_total = math.sqrt(sys_csv * sys_csv + sys_est * sys_est)

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

    for _, row in selected.sort_values("phimin").iterrows():
        sigma = parse_scalar_from_cell(row["cross sections, ep->epg, exp"])
        stat = parse_scalar_from_cell(row["cross sections, ep->epg, exp, stat. unc."])
        sys_up = parse_scalar_from_cell(row["cross sections, ep->epg, exp, syst. unc. (up)"])
        sys_down = parse_scalar_from_cell(row["cross sections, ep->epg, exp, syst. unc. (down)"])

        if not np.isfinite(sigma):
            continue
        # endif

        stat = finite_or_zero(stat)
        sys_up = finite_or_zero(sys_up)
        sys_down = finite_or_zero(sys_down)

        if is_zero_placeholder(sigma, stat, sys_up, sys_down):
            continue
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
# Interpolation, ratios, chi2.
# -----------------------------------------------------------------------------

def periodic_interp(phi_query: float, phi_values: np.ndarray, y_values: np.ndarray) -> float:
    phi_mod = phi_query % 360.0

    phi_base = np.array(phi_values, dtype=float)
    y_base = np.array(y_values, dtype=float)

    if len(phi_base) < 2:
        return math.nan
    # endif

    order = np.argsort(phi_base)
    phi_base = phi_base[order]
    y_base = y_base[order]

    phi_ext = np.concatenate([phi_base - 360.0, phi_base, phi_base + 360.0])
    y_ext = np.concatenate([y_base, y_base, y_base])

    return float(np.interp(phi_mod, phi_ext, y_ext))


def interpolate_reference_point(phi_query: float, reference_points: List[DataPoint]) -> DataPoint:
    phi, sigma, err_low, err_high = points_to_arrays(reference_points)

    return DataPoint(
        phi=phi_query,
        sigma=periodic_interp(phi_query, phi, sigma),
        err_low=periodic_interp(phi_query, phi, err_low),
        err_high=periodic_interp(phi_query, phi, err_high),
    )


def ratio_to_reference(
    numerator_points: List[DataPoint],
    reference_points: List[DataPoint],
) -> List[DataPoint]:
    ratios: List[DataPoint] = []

    for num in numerator_points:
        ref = interpolate_reference_point(
            phi_query=num.phi,
            reference_points=reference_points,
        )

        if not np.isfinite(ref.sigma) or ref.sigma == 0.0:
            continue
        # endif

        ratio = num.sigma / ref.sigma

        ratio_err_high = math.sqrt(
            (num.err_high / ref.sigma) ** 2
            + (num.sigma * ref.err_low / (ref.sigma * ref.sigma)) ** 2
        )

        ratio_err_low = math.sqrt(
            (num.err_low / ref.sigma) ** 2
            + (num.sigma * ref.err_high / (ref.sigma * ref.sigma)) ** 2
        )

        ratios.append(
            DataPoint(
                phi=num.phi,
                sigma=ratio,
                err_low=ratio_err_low,
                err_high=ratio_err_high,
            )
        )
    # endfor

    ratios.sort(key=lambda p: p.phi)
    return ratios


def chi2_ndf_to_reference(
    comparison_points: List[DataPoint],
    reference_points: List[DataPoint],
) -> Tuple[float, int, float]:
    chi2 = 0.0
    ndf = 0

    for comp in comparison_points:
        ref = interpolate_reference_point(
            phi_query=comp.phi,
            reference_points=reference_points,
        )

        if not np.isfinite(ref.sigma):
            continue
        # endif

        err_comp = symmetric_error(comp)
        err_ref = symmetric_error(ref)

        variance = err_comp * err_comp + err_ref * err_ref

        if variance <= 0.0 or not np.isfinite(variance):
            continue
        # endif

        residual = comp.sigma - ref.sigma
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


def auto_ratio_ylim(ax, points_by_label: Dict[str, List[DataPoint]]) -> None:
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


def make_plot(
    clas6_kin: KinematicPoint,
    clas6_points: List[DataPoint],
    pass2_avg: KinematicAverages,
    points_by_label: Dict[str, List[DataPoint]],
    ratios_by_label: Dict[str, List[DataPoint]],
    chi2_by_label: Dict[str, Tuple[float, int, float]],
    output_dir: str,
    y_units: str,
    use_log_cross_section: bool,
) -> None:
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(15.0, 10.0),
        constrained_layout=True,
        sharex=True,
    )

    fig.suptitle(
        (
            f"CLAS6 / CLAS12 DVCS cross-section cross-check: note label {COMPARISON_LABEL}\n"
            rf"CLAS6: $\langle x_B\rangle={clas6_kin.xB:.3f}$, "
            rf"$\langle Q^2\rangle={clas6_kin.Q2:.3f}~{{\rm GeV}}^2$, "
            rf"$\langle |t|\rangle={clas6_kin.t_abs:.3f}~{{\rm GeV}}^2$"
            "\n"
            rf"pass-2 bin {PASS2_BIN_NAME}: "
            rf"$\langle x_B\rangle={pass2_avg.xB:.3f}$, "
            rf"$\langle Q^2\rangle={pass2_avg.Q2:.3f}~{{\rm GeV}}^2$, "
            rf"$\langle |t|\rangle={pass2_avg.t_abs:.3f}~{{\rm GeV}}^2$"
        ),
        fontsize=13,
    )

    top_left = axes[0, 0]
    top_right = axes[0, 1]
    bottom_left = axes[1, 0]
    bottom_right = axes[1, 1]

    # Top-left.
    plot_dataset(
        ax=top_left,
        points=clas6_points,
        label=CLAS6_DISPLAY_LABEL,
        legend_label=CLAS6_DISPLAY_LABEL,
    )

    for label in TOP_CROSS_SECTION_SERIES:
        if label not in points_by_label:
            continue
        # endif

        plot_dataset(
            ax=top_left,
            points=points_by_label[label],
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    top_left.set_ylabel(rf"$d^4\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ [{y_units}]")
    top_left.set_title("Cross sections: pass-2 combined, pass-1, and CLAS6")
    top_left.grid(True, alpha=0.25)
    top_left.legend(fontsize=8, frameon=True)

    if use_log_cross_section:
        top_left.set_yscale("log")
    # endif

    # Top-right.
    top_ratios = {}

    for label in TOP_RATIO_SERIES:
        if label not in ratios_by_label:
            continue
        # endif

        top_ratios[label] = ratios_by_label[label]

        plot_dataset(
            ax=top_right,
            points=ratios_by_label[label],
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    top_right.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--", zorder=0)
    top_right.set_ylabel(r"pass-2 CLAS12 / CLAS6")
    top_right.set_title("Ratio to CLAS6: pass-2 combined only")
    top_right.grid(True, alpha=0.25)
    top_right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(top_right, top_ratios)

    # Bottom-left.
    plot_dataset(
        ax=bottom_left,
        points=clas6_points,
        label=CLAS6_DISPLAY_LABEL,
        legend_label=CLAS6_DISPLAY_LABEL,
    )

    for label in BOTTOM_CROSS_SECTION_SERIES:
        plot_dataset(
            ax=bottom_left,
            points=points_by_label.get(label, []),
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    bottom_left.set_xlabel(r"$\phi$ [deg]")
    bottom_left.set_ylabel(rf"$d^4\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ [{y_units}]")
    bottom_left.set_title("Cross sections: individual pass-2 periods and CLAS6")
    bottom_left.grid(True, alpha=0.25)
    bottom_left.legend(fontsize=8, frameon=True)

    if use_log_cross_section:
        bottom_left.set_yscale("log")
    # endif

    # Bottom-right.
    bottom_ratios = {}

    for label in BOTTOM_RATIO_SERIES:
        ratios = ratios_by_label.get(label, [])
        bottom_ratios[label] = ratios

        plot_dataset(
            ax=bottom_right,
            points=ratios,
            label=label,
            legend_label=format_label_with_chi2(
                label,
                chi2_by_label.get(label, (math.nan, 0, math.nan)),
            ),
        )
    # endfor

    bottom_right.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--", zorder=0)
    bottom_right.set_xlabel(r"$\phi$ [deg]")
    bottom_right.set_ylabel(r"pass-2 CLAS12 / CLAS6")
    bottom_right.set_title("Ratio to CLAS6: individual pass-2 periods")
    bottom_right.grid(True, alpha=0.25)
    bottom_right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(bottom_right, bottom_ratios)

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(
        output_dir,
        f"clas6_cross_check_{COMPARISON_LABEL}_p2bin{PASS2_BIN_NAME}.png",
    )

    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote: {output_path}")


# -----------------------------------------------------------------------------
# Diagnostics.
# -----------------------------------------------------------------------------

def print_point_summary(label: str, points: List[DataPoint]) -> None:
    total_values = [symmetric_error(p) for p in points if np.isfinite(symmetric_error(p))]
    stat_values = [p.stat for p in points if np.isfinite(p.stat)]
    sys_values = [p.sys for p in points if np.isfinite(p.sys)]

    mean_total = float(np.mean(total_values)) if len(total_values) > 0 else math.nan
    mean_stat = float(np.mean(stat_values)) if len(stat_values) > 0 else math.nan
    mean_sys = float(np.mean(sys_values)) if len(sys_values) > 0 else math.nan

    print(f"  {label}")
    print(f"    points          = {len(points)}")
    print(f"    mean stat err   = {mean_stat:.6g}")
    print(f"    mean sys err    = {mean_sys:.6g}")
    print(f"    mean total err  = {mean_total:.6g}")

    if len(points) > 0:
        first = points[0]
        print(
            f"    first point     = phi {first.phi:.3f} deg, "
            f"sigma {first.sigma:.6g}, "
            f"stat {first.stat:.6g}, "
            f"sys {first.sys:.6g}, "
            f"total {symmetric_error(first):.6g}"
        )
    # endif


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare pass-2 CLAS12 DVCS cross sections to pass-1 and CLAS6 for one hard-coded bin."
    )

    parser.add_argument(
        "pass2_csv",
        help="Input final pass-2 DVCS analysis CSV.",
    )

    parser.add_argument(
        "--pass1-csv",
        default=DEFAULT_PASS1_CSV,
        help=f"Input pass-1 CSV. Default: {DEFAULT_PASS1_CSV}",
    )

    parser.add_argument(
        "--no-pass1",
        action="store_true",
        help="Disable pass-1 Fa18 overlay.",
    )

    parser.add_argument(
        "--clas6-dataset-id",
        type=int,
        default=DEFAULT_CLAS6_DATASET_ID,
        help=f"Gepard dataset id for CLAS6 unpolarized cross sections. Default: {DEFAULT_CLAS6_DATASET_ID}.",
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--pass2-scale",
        type=float,
        default=DEFAULT_PASS2_SCALE,
        help=f"Scale factor for pass-2 CLAS12 values. Default: {DEFAULT_PASS2_SCALE:g}.",
    )

    parser.add_argument(
        "--pass1-scale",
        type=float,
        default=DEFAULT_PASS1_SCALE,
        help=f"Scale factor for pass-1 CLAS12 values. Default: {DEFAULT_PASS1_SCALE:g}.",
    )

    parser.add_argument(
        "--clas6-scale",
        type=float,
        default=DEFAULT_CLAS6_SCALE,
        help=f"Scale factor for CLAS6 Gepard values. Default: {DEFAULT_CLAS6_SCALE:g}.",
    )

    parser.add_argument(
        "--y-units",
        default=r"nb/GeV$^4$",
        help=r"Y-axis unit label. Default: nb/GeV$^4$.",
    )

    parser.add_argument(
        "--avg-period",
        default="10.6 GeV",
        help="Pass-2 period used for weighted average kinematics. Default: '10.6 GeV'.",
    )

    parser.add_argument(
        "--no-clas12-estimated-sys",
        action="store_true",
        help="Disable the extra estimated pass-2 CLAS12 bin-to-bin systematic.",
    )

    parser.add_argument(
        "--clas12-bin-to-bin-sys-frac",
        type=float,
        default=DEFAULT_PASS2_BIN_TO_BIN_SYS_FRAC,
        help=f"Estimated fractional pass-2 systematic. Default: {DEFAULT_PASS2_BIN_TO_BIN_SYS_FRAC:.2f}.",
    )

    parser.add_argument(
        "--pass1-norm-sys-frac",
        type=float,
        default=DEFAULT_PASS1_NORM_SYS_FRAC,
        help=f"Additional pass-1 normalization systematic. Default: {DEFAULT_PASS1_NORM_SYS_FRAC:.2f}.",
    )

    parser.add_argument(
        "--clas6-phi-convention",
        choices=["trento", "gepard"],
        default="trento",
        help=(
            "Use 'trento' to convert Gepard BMK phi back to published-style Trento phi "
            "for Trento-framed data; use 'gepard' to leave phi untouched. Default: trento."
        ),
    )

    parser.add_argument(
        "--clas6-tolerance-xB",
        type=float,
        default=1.0e-6,
        help="Tolerance for selecting the CLAS6 xB group. Default: 1e-6.",
    )

    parser.add_argument(
        "--clas6-tolerance-Q2",
        type=float,
        default=1.0e-6,
        help="Tolerance for selecting the CLAS6 Q2 group. Default: 1e-6.",
    )

    parser.add_argument(
        "--clas6-tolerance-t",
        type=float,
        default=1.0e-6,
        help="Tolerance for selecting the CLAS6 |t| group. Default: 1e-6.",
    )

    parser.add_argument(
        "--linear-cross-section",
        action="store_true",
        help="Use linear y-scale for cross-section panels. Default is log scale.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.pass2_scale <= 0.0:
        raise ValueError("--pass2-scale must be positive.")
    # endif

    if args.pass1_scale <= 0.0:
        raise ValueError("--pass1-scale must be positive.")
    # endif

    if args.clas6_scale <= 0.0:
        raise ValueError("--clas6-scale must be positive.")
    # endif

    if args.clas12_bin_to_bin_sys_frac < 0.0:
        raise ValueError("--clas12-bin-to-bin-sys-frac must be non-negative.")
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

    pass2_df = pd.read_csv(args.pass2_csv)
    require_columns(pass2_df, pass2_required, context="pass-2")

    pass2_selected = select_bin_by_name(
        df=pass2_df,
        bin_name=PASS2_BIN_NAME,
    )

    validate_pass2_bin_edges(pass2_selected)

    pass2_avg = compute_weighted_kinematic_averages(
        selected=pass2_selected,
        csv_period=args.avg_period,
    )

    pass1_df: Optional[pd.DataFrame] = None
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
                "phimin",
                "phimax",
                "phiavg",
                "cross sections, ep->epg, exp",
                "cross sections, ep->epg, exp, stat. unc.",
                "cross sections, ep->epg, exp, syst. unc. (up)",
                "cross sections, ep->epg, exp, syst. unc. (down)",
            ]

            pass1_df = pd.read_csv(args.pass1_csv)
            require_columns(pass1_df, pass1_required, context="pass-1")

            pass1_selected = select_bin_by_name(
                df=pass1_df,
                bin_name=PASS2_BIN_NAME,
            )
        # endif
    # endif

    print("Loading CLAS6 data from Gepard...")
    print(f"  dataset id = {args.clas6_dataset_id}")
    print("  hard-coded CLAS6 target:")
    print(f"    xB={CLAS6_TARGET_XB:.6f}")
    print(f"    Q2={CLAS6_TARGET_Q2:.6f}")
    print(f"    |t|={CLAS6_TARGET_T_ABS:.6f}")

    (
        clas6_kin,
        clas6_points,
        n_raw_matching_kinematics,
        n_zero_placeholders_removed,
    ) = load_clas6_points_for_target(
        dataset_id=args.clas6_dataset_id,
        clas6_scale=args.clas6_scale,
        phi_convention=args.clas6_phi_convention,
        target_xB=CLAS6_TARGET_XB,
        target_Q2=CLAS6_TARGET_Q2,
        target_t_abs=CLAS6_TARGET_T_ABS,
        tolerance_xB=args.clas6_tolerance_xB,
        tolerance_Q2=args.clas6_tolerance_Q2,
        tolerance_t=args.clas6_tolerance_t,
    )

    print("Selected CLAS6 group:")
    print(f"  raw matching kinematic rows = {n_raw_matching_kinematics}")
    print(f"  zero placeholders removed   = {n_zero_placeholders_removed}")
    print(f"  plotted CLAS6 points        = {len(clas6_points)}")
    print(f"  <xB>                        = {clas6_kin.xB:.6f}")
    print(f"  <Q2>                        = {clas6_kin.Q2:.6f} GeV^2")
    print(f"  <|t|>                       = {clas6_kin.t_abs:.6f} GeV^2")

    print()
    print("Selected pass-2 CLAS12 bin:")
    print(f"  Bin Name                    = {PASS2_BIN_NAME}")
    print(f"  rows                        = {len(pass2_selected)}")
    print(f"  expected xB bin             = [{PASS2_EXPECTED_XB_MIN}, {PASS2_EXPECTED_XB_MAX}]")
    print(f"  expected Q2 bin             = [{PASS2_EXPECTED_Q2_MIN}, {PASS2_EXPECTED_Q2_MAX}] GeV^2")
    print(f"  expected |t| bin            = [{PASS2_EXPECTED_T_ABS_MIN}, {PASS2_EXPECTED_T_ABS_MAX}] GeV^2")
    print(f"  weighted <xB>               = {pass2_avg.xB:.6f}")
    print(f"  weighted <Q2>               = {pass2_avg.Q2:.6f} GeV^2")
    print(f"  weighted <|t|>              = {pass2_avg.t_abs:.6f} GeV^2")
    print(f"  weighting period            = {pass2_avg.weight_period}")

    if pass1_selected is not None:
        print()
        print("Selected pass-1 CLAS12 bin:")
        print(f"  Bin Name                    = {PASS2_BIN_NAME}")
        print(f"  rows                        = {len(pass1_selected)}")
    # endif

    points_by_label: Dict[str, List[DataPoint]] = {}

    points_by_label[PASS2_COMBINED_DISPLAY_LABEL] = pass2_points_for_period(
        selected=pass2_selected,
        csv_period=PASS2_COMBINED_CSV_PERIOD,
        pass2_scale=args.pass2_scale,
        include_pass2_estimated_sys=not args.no_clas12_estimated_sys,
        pass2_bin_to_bin_sys_frac=args.clas12_bin_to_bin_sys_frac,
    )

    for display_label, csv_period in PASS2_PERIOD_DISPLAY_TO_CSV_PERIOD.items():
        points_by_label[display_label] = pass2_points_for_period(
            selected=pass2_selected,
            csv_period=csv_period,
            pass2_scale=args.pass2_scale,
            include_pass2_estimated_sys=not args.no_clas12_estimated_sys,
            pass2_bin_to_bin_sys_frac=args.clas12_bin_to_bin_sys_frac,
        )
    # endfor

    if pass1_selected is not None:
        points_by_label[PASS1_DISPLAY_LABEL] = pass1_points(
            selected=pass1_selected,
            pass1_scale=args.pass1_scale,
            pass1_norm_sys_frac=args.pass1_norm_sys_frac,
        )
    # endif

    chi2_by_label = {
        label: chi2_ndf_to_reference(
            comparison_points=points,
            reference_points=clas6_points,
        )
        for label, points in points_by_label.items()
    }

    ratios_by_label = {
        label: ratio_to_reference(
            numerator_points=points,
            reference_points=clas6_points,
        )
        for label, points in points_by_label.items()
        if label != PASS1_DISPLAY_LABEL
    }

    print()
    print("Dataset summaries:")
    print_point_summary(CLAS6_DISPLAY_LABEL, clas6_points)

    for label in TOP_CROSS_SECTION_SERIES + BOTTOM_CROSS_SECTION_SERIES:
        if label in points_by_label:
            print_point_summary(label, points_by_label[label])
        # endif
    # endfor

    print()
    print("Chi2/ndf against CLAS6:")
    for label, chi2_info in chi2_by_label.items():
        chi2, ndf, chi2ndf = chi2_info

        if ndf > 0 and np.isfinite(chi2ndf):
            print(f"  {label:18s}: chi2={chi2:.4f}, ndf={ndf:d}, chi2/ndf={chi2ndf:.4f}")
        else:
            print(f"  {label:18s}: chi2/ndf=N/A")
        # endif
    # endfor

    make_plot(
        clas6_kin=clas6_kin,
        clas6_points=clas6_points,
        pass2_avg=pass2_avg,
        points_by_label=points_by_label,
        ratios_by_label=ratios_by_label,
        chi2_by_label=chi2_by_label,
        output_dir=args.output_dir,
        y_units=args.y_units,
        use_log_cross_section=not args.linear_cross_section,
    )


if __name__ == "__main__":
    main()
# endif