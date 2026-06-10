#!/usr/bin/env python3
"""
integrated_cross_section_study.py

Standalone projection/integration study for the DVCS pass-2 CSV.

The script reads the final analysis CSV and makes four 1x3 canvases in
output/integrated_study/:

  integrated_xB_dependence.png
  integrated_Q2_dependence.png
  integrated_t_dependence.png
  integrated_phi_dependence.png

Each canvas contains:

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

By default the integration is bin-width weighted, appropriate for differential
cross sections:

  d sigma / dxB       = sum_i sigma_i * dQ2_i * dt_i * dphi_i
  d sigma / dQ2       = sum_i sigma_i * dxB_i * dt_i * dphi_i
  d sigma / d|t|      = sum_i sigma_i * dxB_i * dQ2_i * dphi_i
  d sigma / dphi      = sum_i sigma_i * dxB_i * dQ2_i * dt_i

The phi width is converted from degrees to radians by default, since the usual
DVCS differential cross section convention is per radian. Use --phi-degrees if
you intentionally want degree-weighted phi integration instead.

For the plotted x-position in each projected bin, the script tries to use the
available per-bin average kinematic columns, for example:

  xBavg, Fa18 Inb
  Q2avg, Fa18 Inb
  t_abs_avg, Fa18 Inb
  phiavg, Fa18 Inb

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

The row-by-row correction factor is:

  C_i = KM15(E_target, xB_i, Q2_i, |t|_i, phi_i)
      / KM15(E_source, xB_i, Q2_i, |t|_i, phi_i)

and the scaled Sp19 contribution is:

  sigma_i_scaled = C_i * sigma_i
  stat_i_scaled  = |C_i| * stat_i
  sys_i_scaled   = |C_i| * sys_i

This is enabled by default. Disable with:

  --no-sp19-km15-energy-scaling
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
# Beam energies.
# -----------------------------------------------------------------------------

FA18_INB_EBEAM_GEV = 10.6040
SP19_INB_EBEAM_GEV = 10.1998


# -----------------------------------------------------------------------------
# Plot order requested by user.
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

RIGHT_LABEL_TO_PERIOD = {
    "Fa18 Inb, 10.604 GeV": "Fa18 Inb",
    "Sp19 Inb -> 10.604 GeV, KM15": "Sp19 Inb",
}


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
}


PROJECTIONS = {
    "xB": {
        "min_col": "xBmin",
        "max_col": "xBmax",
        "avg_prefix": "xBavg",
        "xlabel": r"$x_B$",
        "ylabel": r"$d\sigma/dx_B$  [pb]",
        "title": r"$x_B$ dependence",
        "integrate_widths": ["Q2", "t", "phi"],
        "outfile_tag": "xB",
    },
    "Q2": {
        "min_col": "Q2min",
        "max_col": "Q2max",
        "avg_prefix": "Q2avg",
        "xlabel": r"$Q^2$  [GeV$^2$]",
        "ylabel": r"$d\sigma/dQ^2$  [pb/GeV$^2$]",
        "title": r"$Q^2$ dependence",
        "integrate_widths": ["xB", "t", "phi"],
        "outfile_tag": "Q2",
    },
    "t": {
        "min_col": "t_abs_min",
        "max_col": "t_abs_max",
        "avg_prefix": "t_abs_avg",
        "xlabel": r"$|t|$  [GeV$^2$]",
        "ylabel": r"$d\sigma/d|t|$  [pb/GeV$^2$]",
        "title": r"$|t|$ dependence",
        "integrate_widths": ["xB", "Q2", "phi"],
        "outfile_tag": "t",
    },
    "phi": {
        "min_col": "phimin",
        "max_col": "phimax",
        "avg_prefix": "phiavg",
        "xlabel": r"$\phi$  [deg]",
        "ylabel": r"$d\sigma/d\phi$  [pb/rad]",
        "title": r"$\phi$ dependence",
        "integrate_widths": ["xB", "Q2", "t"],
        "outfile_tag": "phi",
    },
}


FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)


@dataclass
class Point:
    key: Tuple[float, float]
    x: float
    y: float
    stat: float
    sys: float


@dataclass
class ModelContext:
    enabled: bool
    g: object = None
    th_km15: object = None
    correction_min: float = math.inf
    correction_max: float = -math.inf
    correction_sum: float = 0.0
    correction_count: int = 0


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


def update_model_correction_diagnostics(
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


def sp19_to_fa18_energy_correction(
    row: pd.Series,
    period: str,
    model_context: ModelContext,
) -> float:
    """
    Compute the KM15 beam-energy correction:

      C = KM15(E=10.6040 GeV) / KM15(E=10.1998 GeV)

    using the row's period-specific average kinematics.

    If a required average column is missing or invalid, the function falls back
    to the bin midpoint for that variable.
    """

    if not model_context.enabled:
        return 1.0
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
        return math.nan
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
    update_model_correction_diagnostics(model_context, correction)

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
    Parse a scalar-like CSV cell.

    If the cell is a tuple such as '(value, stat, sys)', this returns the first
    number. This is useful for yield/count columns that may be stored as tuple3.
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
    """
    Return candidate yield/count columns that can be used to weight the average
    x-position of an integrated/projection point.
    """

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
    """
    Determine the event/yield weight for placing the x-position.

    Returns NaN if no useful positive weight is found.
    """

    for col in candidate_columns:
        val = parse_scalar_from_cell(row[col])

        if np.isfinite(val) and val > 0.0:
            return float(val)
        # endif
    # endfor

    return math.nan


def average_x_for_group(
    group: pd.DataFrame,
    projection_info: Dict[str, object],
    period: str,
    yield_weight_columns: List[str],
) -> float:
    """
    Compute the plotted x-position for a projected bin.

    Preferred behavior:
      x_plot = sum_j N_j * <x>_j / sum_j N_j

    where j runs over the hidden sub-bins being integrated over. The N_j are
    taken from pi0-corrected/signal/acceptance-corrected yield columns if found.

    Fallbacks:
      1. simple finite average of the per-sub-bin average kinematic values,
      2. midpoint of the projected bin.
    """

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
# Integration.
# -----------------------------------------------------------------------------

def integrated_points_for_period(
    df: pd.DataFrame,
    projection: str,
    period: str,
    no_width_weighting: bool,
    model_context: Optional[ModelContext] = None,
    apply_sp19_to_fa18_scaling: bool = False,
) -> List[Point]:
    info = PROJECTIONS[projection]
    col = cross_section_column(period)

    if col not in df.columns:
        return []
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
                if model_context is None or not model_context.enabled:
                    raise RuntimeError(
                        "Sp19 -> Fa18 KM15 scaling requested, but the model context is disabled."
                    )
                # endif

                model_scale = sp19_to_fa18_energy_correction(
                    row=row,
                    period=period,
                    model_context=model_context,
                )

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


def point_total_error(
    point: Point,
    include_bin_to_bin_sys: bool,
    bin_to_bin_sys_frac: float,
) -> float:
    """
    Return the outer y-error to draw for a cross-section or ratio point.

    If include_bin_to_bin_sys is false:
      total_error = stat

    If include_bin_to_bin_sys is true:
      total_error = sqrt(stat^2 + (f * y)^2)

    where f defaults to 0.10.
    """

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
    """
    Weighted mean of two points using statistical uncertainties only.

    mean = sum_i y_i / stat_i^2 / sum_i 1 / stat_i^2
    stat = sqrt(1 / sum_i 1 / stat_i^2)

    If one point has invalid or zero uncertainty, it is skipped. If both are
    unusable as weighted measurements, return None.
    """

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


def auto_ratio_ylim(ax) -> None:
    ylo, yhi = ax.get_ylim()

    if not np.isfinite(yhi):
        ax.set_ylim(0.0, 2.0)
        return
    # endif

    upper = max(1.25, 1.10 * yhi)
    ax.set_ylim(0.0, upper)


def make_projection_plot(
    df: pd.DataFrame,
    projection: str,
    output_dir: str,
    no_width_weighting: bool,
    include_bin_to_bin_sys: bool,
    bin_to_bin_sys_frac: float,
    model_context: ModelContext,
) -> None:
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
            model_context=model_context,
            apply_sp19_to_fa18_scaling=model_context.enabled,
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
    auto_ratio_ylim(middle)

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
    auto_ratio_ylim(right)

    outbase = os.path.join(
        output_dir,
        f"integrated_{info['outfile_tag']}_dependence",
    )

    fig.savefig(outbase + ".png", dpi=200)
    plt.close(fig)


def print_model_correction_summary(model_context: ModelContext) -> None:
    if not model_context.enabled:
        print("Sp19 -> 10.604 GeV KM15 beam-energy scaling: disabled")
        return
    # endif

    if model_context.correction_count <= 0:
        print("Sp19 -> 10.604 GeV KM15 beam-energy scaling: enabled, but no corrections were evaluated")
        return
    # endif

    mean_corr = model_context.correction_sum / model_context.correction_count

    print("Sp19 -> 10.604 GeV KM15 beam-energy scaling:")
    print(f"  source beam energy = {SP19_INB_EBEAM_GEV:.6f} GeV")
    print(f"  target beam energy = {FA18_INB_EBEAM_GEV:.6f} GeV")
    print(f"  correction count   = {model_context.correction_count}")
    print(f"  correction mean    = {mean_corr:.6g}")
    print(f"  correction range   = [{model_context.correction_min:.6g}, {model_context.correction_max:.6g}]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make integrated xB, Q2, |t|, and phi dependence plots from a "
            "final DVCS pass-2 analysis CSV."
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

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.bin_to_bin_sys_frac < 0.0:
        raise ValueError("--bin-to-bin-sys-frac must be non-negative.")
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

    model_context = make_model_context(
        enabled=not args.no_sp19_km15_energy_scaling,
    )

    os.makedirs(args.output_dir, exist_ok=True)

    for projection in ["xB", "Q2", "t", "phi"]:
        make_projection_plot(
            df=df,
            projection=projection,
            output_dir=args.output_dir,
            no_width_weighting=args.no_width_weighting,
            include_bin_to_bin_sys=args.include_bin_to_bin_sys,
            bin_to_bin_sys_frac=args.bin_to_bin_sys_frac,
            model_context=model_context,
        )
    # endfor

    if args.include_bin_to_bin_sys:
        print(
            "Using plotted y-errors: statistical bars plus outer total bars "
            f"sqrt(stat^2 + ({args.bin_to_bin_sys_frac:.4f} * y)^2)"
        )
    else:
        print("Using plotted y-errors: statistical only")
    # endif

    print_model_correction_summary(model_context)

    print(f"Wrote integrated study PNG plots to: {args.output_dir}")


if __name__ == "__main__":
    main()
# endif