#!/usr/bin/env python3
"""
hall_a_cross_check.py

Standalone Hall A / CLAS12 DVCS cross-section cross-check.

This script reads the final DVCS pass-2 analysis CSV and compares the CLAS12
cross sections in one overlapping bin against the Hall A Kin363 cross-section
table.

The Hall A values are hard-coded from the right-most column of the Hall A table:

  <xB>  = 0.371
  <Q2>  = 4.568 GeV^2
  <t'>  = -0.303 GeV^2

with cross sections listed in pb/GeV^4 as:

  sigma +/- stat ^{+ sys_up}_{- sys_down}

The corresponding CLAS12 bin is selected by default using:

  Bin Name = 137

which corresponds to:

  xB  in [0.357, 0.446]
  Q2  in [4.326, 5.761] GeV^2
  |t| in [0.25, 0.40] GeV^2

The script makes a 1x2 canvas:

  left:  Hall B / CLAS12 and Hall A cross sections vs phi
  right: Hall B / Hall A ratios vs phi

Both panels show stat+sys error bars.

Important unit convention:
  The Hall A table is in pb/GeV^4. The current CLAS12 CSV cross-section values
  in this bin appear to be in nb/GeV^4 relative to the Hall A table, so the
  script multiplies CLAS12 cross sections by 1000 by default.

  Override with:
    --clas12-scale 1.0

The CLAS12 and Hall A phi points are not identical. For the ratio panel, Hall A
is periodically linearly interpolated to the CLAS12 phi value before forming:

  ratio = CLAS12(phi_CLAS12) / HallA_interpolated(phi_CLAS12)

Output:

  output/hall_a_cross_check/hall_a_cross_check.png
"""

import argparse
import math
import os
import re
from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Default CLAS12 bin corresponding to the Hall A overlap point.
# -----------------------------------------------------------------------------

DEFAULT_BIN_NAME = 137

TARGET_XB_MIN = 0.357
TARGET_XB_MAX = 0.446

TARGET_Q2_MIN = 4.326
TARGET_Q2_MAX = 5.761

TARGET_T_ABS_MIN = 0.25
TARGET_T_ABS_MAX = 0.40


# -----------------------------------------------------------------------------
# Plot order.
# -----------------------------------------------------------------------------

HALL_B_SERIES = [
    "10.6 GeV",
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
    "Hall A",
]


SERIES_STYLE = {
    "10.6 GeV": dict(marker="o", linestyle="-", color="black"),
    "Sp18 Inb": dict(marker="D", linestyle="-", color="tab:green"),
    "Sp18 Out": dict(marker="P", linestyle="-", color="tab:purple"),
    "Fa18 Inb": dict(marker="^", linestyle="-", color="tab:blue"),
    "Fa18 Out": dict(marker="v", linestyle="-", color="tab:orange"),
    "Hall A": dict(marker="s", linestyle="-", color="tab:red"),
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


def cross_section_column(period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {period}, unpol"


def average_phi_column(period: str) -> str:
    return f"phiavg, {period}"


def require_columns(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [c for c in columns if c not in df.columns]

    if missing:
        raise KeyError(
            "The input CSV is missing required columns:\n  " + "\n  ".join(missing)
        )
    # endif


def float_close(series: pd.Series, value: float, tolerance: float = 1.0e-6) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    return np.isclose(numeric, value, rtol=0.0, atol=tolerance)


def select_hall_b_overlap_bin(df: pd.DataFrame, bin_name: Optional[int]) -> pd.DataFrame:
    """
    Select the CLAS12 bin overlapping the Hall A Kin363 point.

    Preferred selection:
      Bin Name == 137

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
            "Could not find the requested Hall A overlap bin in the CSV:\n"
            f"  Bin Name = {bin_name}\n"
            f"  xB  [{TARGET_XB_MIN}, {TARGET_XB_MAX}]\n"
            f"  Q2  [{TARGET_Q2_MIN}, {TARGET_Q2_MAX}]\n"
            f"  |t| [{TARGET_T_ABS_MIN}, {TARGET_T_ABS_MAX}]\n"
        )
    # endif

    selected = selected.sort_values(["phimin", "phimax"]).copy()

    return selected


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
            )
        )
    # endfor

    points.sort(key=lambda p: p.phi)

    return points


def hall_b_points_for_period(
    selected: pd.DataFrame,
    period: str,
    clas12_scale: float,
) -> List[DataPoint]:
    """
    Extract Hall B / CLAS12 phi-dependent points for one period.

    Cross-section cells are assumed to be stored as:

      (value, stat, sys)

    The plotted error is:

      sqrt(stat^2 + sys^2)

    The value, stat, and sys are all multiplied by clas12_scale.

    If the period-specific phi average column exists, it is used for the x-position.
    Otherwise, the phi-bin midpoint is used.
    """

    xs_col = cross_section_column(period)
    phi_avg_col = average_phi_column(period)

    if xs_col not in selected.columns:
        return []
    # endif

    points: List[DataPoint] = []

    sorted_df = selected.sort_values("phimin")

    for _, row in sorted_df.iterrows():
        sigma, stat, sys = parse_tuple3(row[xs_col])

        if not np.isfinite(sigma):
            continue
        # endif

        if not np.isfinite(stat):
            stat = 0.0
        # endif

        if not np.isfinite(sys):
            sys = 0.0
        # endif

        sigma *= clas12_scale
        stat *= clas12_scale
        sys *= clas12_scale

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

        err = math.sqrt(stat * stat + sys * sys)

        points.append(
            DataPoint(
                phi=phi,
                sigma=sigma,
                err_low=err,
                err_high=err,
            )
        )
    # endfor

    points.sort(key=lambda p: p.phi)

    return points


def points_to_arrays(points: List[DataPoint]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    phi = np.array([p.phi for p in points], dtype=float)
    sigma = np.array([p.sigma for p in points], dtype=float)
    err_low = np.array([p.err_low for p in points], dtype=float)
    err_high = np.array([p.err_high for p in points], dtype=float)

    return phi, sigma, err_low, err_high


def periodic_interp(phi_query: float, phi_values: np.ndarray, y_values: np.ndarray) -> float:
    """
    Periodically linearly interpolate y(phi) at phi_query.

    phi is treated as periodic over 360 degrees.
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

    This is used only for the ratio panel because the CLAS12 phi averages are not
    identical to the Hall A phi centers.
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
    )


def ratio_to_hall_a(
    hall_b_points: List[DataPoint],
    hall_a_points_input: List[DataPoint],
) -> List[DataPoint]:
    """
    Compute Hall B / Hall A ratios.

    Hall A is periodically linearly interpolated to each CLAS12 phi value.

    Error propagation uses asymmetric Hall A uncertainty:

      R = B / A

      delta_R_up^2 =
          (delta_B_up / A)^2
        + (B * delta_A_down / A^2)^2

      delta_R_down^2 =
          (delta_B_down / A)^2
        + (B * delta_A_up / A^2)^2

    This convention follows the fact that a lower Hall A denominator increases
    the ratio, while an upper Hall A denominator decreases the ratio.
    """

    ratios: List[DataPoint] = []

    for b in hall_b_points:
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
            )
        )
    # endfor

    ratios.sort(key=lambda p: p.phi)

    return ratios


def plot_dataset(ax, points: List[DataPoint], label: str) -> None:
    if len(points) == 0:
        return
    # endif

    phi, sigma, err_low, err_high = points_to_arrays(points)
    style = SERIES_STYLE.get(label, dict(marker="o", linestyle="-"))

    ax.errorbar(
        phi,
        sigma,
        yerr=np.vstack([err_low, err_high]),
        label=label,
        markersize=5.5,
        linewidth=1.2,
        elinewidth=1.0,
        capsize=2.5,
        **style,
    )


def auto_ratio_ylim(ax, points_by_period: dict) -> None:
    yvals: List[float] = []

    for points in points_by_period.values():
        for p in points:
            if np.isfinite(p.sigma):
                yvals.append(p.sigma)
            # endif
        # endfor
    # endfor

    if len(yvals) == 0:
        ax.set_ylim(0.0, 2.0)
        return
    # endif

    ymin = min(yvals)
    ymax = max(yvals)

    span_low = abs(1.0 - ymin)
    span_high = abs(ymax - 1.0)
    span = max(span_low, span_high, 0.35)

    ax.set_ylim(1.0 - 1.20 * span, 1.0 + 1.20 * span)


def make_plot(
    selected: pd.DataFrame,
    output_dir: str,
    output_name: str,
    clas12_scale: float,
) -> None:
    hall_a = hall_a_points()

    hall_b_by_period = {
        period: hall_b_points_for_period(
            selected=selected,
            period=period,
            clas12_scale=clas12_scale,
        )
        for period in HALL_B_SERIES
    }

    ratios_by_period = {
        period: ratio_to_hall_a(
            hall_b_points=hall_b_by_period.get(period, []),
            hall_a_points_input=hall_a,
        )
        for period in HALL_B_SERIES
    }

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(14.0, 5.5),
        constrained_layout=True,
    )

    fig.suptitle(
        (
            "Hall A / CLAS12 DVCS cross-section cross-check\n"
            r"CLAS12 bin: "
            r"$0.357<x_B<0.446$, "
            r"$4.326<Q^2<5.761~{\rm GeV}^2$, "
            r"$0.25<|t|<0.40~{\rm GeV}^2$"
        ),
        fontsize=13,
    )

    left = axes[0]

    for label in LEFT_SERIES:
        if label == "Hall A":
            plot_dataset(left, hall_a, label)
        else:
            plot_dataset(left, hall_b_by_period.get(label, []), label)
        # endif
    # endfor

    left.set_xlabel(r"$\phi$ [deg]")
    left.set_ylabel(r"$d^4\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ [pb/GeV$^4$]")
    left.set_title("Cross sections")
    left.grid(True, alpha=0.25)
    left.legend(fontsize=8, frameon=True)

    right = axes[1]

    for label in HALL_B_SERIES:
        plot_dataset(right, ratios_by_period.get(label, []), label)
    # endfor

    right.axhline(
        1.0,
        color="0.35",
        linewidth=1.0,
        linestyle="--",
        zorder=0,
    )

    right.set_xlabel(r"$\phi$ [deg]")
    right.set_ylabel(r"CLAS12 / Hall A")
    right.set_title("Ratio to Hall A")
    right.grid(True, alpha=0.25)
    right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(right, ratios_by_period)

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, output_name)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote: {output_path}")

    print()
    print("Extracted CLAS12 points after applying scale factor:")
    print(f"  clas12_scale = {clas12_scale:g}")

    for period in HALL_B_SERIES:
        points = hall_b_by_period.get(period, [])
        print(f"  {period}: {len(points)} points")

        if len(points) > 0:
            first = points[0]
            print(
                f"    first point: phi = {first.phi:.3f} deg, "
                f"sigma = {first.sigma:.6g} pb/GeV^4"
            )
        # endif
    # endfor


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make a Hall A / CLAS12 DVCS cross-section cross-check plot."
    )

    parser.add_argument(
        "csv_file",
        help="Input final DVCS pass-2 analysis CSV.",
    )

    parser.add_argument(
        "--output-dir",
        default="output/hall_a_cross_check",
        help="Directory for output plot. Default: output/hall_a_cross_check",
    )

    parser.add_argument(
        "--output-name",
        default="hall_a_cross_check.png",
        help="Output PNG filename. Default: hall_a_cross_check.png",
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
        default=1000.0,
        help=(
            "Scale factor applied to CLAS12 cross sections and uncertainties. "
            "Default: 1000.0, useful if the CSV is effectively in nb/GeV^4 "
            "and Hall A is in pb/GeV^4."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.clas12_scale <= 0.0:
        raise ValueError("--clas12-scale must be positive.")
    # endif

    required = [
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

    required += [cross_section_column(period) for period in HALL_B_SERIES]

    df = pd.read_csv(args.csv_file)
    require_columns(df, required)

    selected = select_hall_b_overlap_bin(
        df=df,
        bin_name=args.bin_name,
    )

    print("Selected Hall A overlap bin from CSV:")
    print(f"  rows: {len(selected)}")
    print(f"  Bin Name: {args.bin_name}")
    print(f"  xB:  [{TARGET_XB_MIN}, {TARGET_XB_MAX}]")
    print(f"  Q2:  [{TARGET_Q2_MIN}, {TARGET_Q2_MAX}] GeV^2")
    print(f"  |t|: [{TARGET_T_ABS_MIN}, {TARGET_T_ABS_MAX}] GeV^2")
    print(f"  phi bins present: {len(selected)}")

    make_plot(
        selected=selected,
        output_dir=args.output_dir,
        output_name=args.output_name,
        clas12_scale=args.clas12_scale,
    )


if __name__ == "__main__":
    main()
# endif