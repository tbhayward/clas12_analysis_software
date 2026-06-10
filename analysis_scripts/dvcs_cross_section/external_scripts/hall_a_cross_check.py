#!/usr/bin/env python3
"""
hall_a_cross_check.py

Standalone Hall A / Hall B DVCS cross-section cross-check.

This script reads the final DVCS pass-2 analysis CSV and compares the CLAS12
cross sections in one overlapping bin against the Hall A Kin363 cross-section
table.

The Hall A values are hard-coded from the right-most column of Table D.5:

  <xB>  = 0.371
  <Q2>  = 4.568 GeV^2
  <t'>  = -0.303 GeV^2

with cross sections listed in pb/GeV^4 as:

  sigma +/- stat ^{+ sys_up}_{- sys_down}

The corresponding CLAS12 bin is selected as:

  xB  in [0.357, 0.446]
  Q2  in [4.326, 5.761] GeV^2
  |t| in [0.25, 0.40] GeV^2

The script makes a 1x2 canvas:

  left:  Hall B / CLAS12 and Hall A cross sections vs phi
  right: Hall B / Hall A ratios vs phi

Both panels show stat+sys error bars. Hall A errors are asymmetric because the
published systematic uncertainty is asymmetric.

Output:

  output/hall_a_cross_check/hall_a_cross_check.png
"""

import argparse
import math
import os
import re
from dataclasses import dataclass
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# CLAS12 bin corresponding to the Hall A overlap point.
# -----------------------------------------------------------------------------

TARGET_XB_MIN = 0.357
TARGET_XB_MAX = 0.446

TARGET_Q2_MIN = 4.326
TARGET_Q2_MAX = 5.761

TARGET_T_ABS_MIN = 0.25
TARGET_T_ABS_MAX = 0.40


# -----------------------------------------------------------------------------
# Plot order requested by user.
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


def select_hall_b_overlap_bin(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select the CLAS12 bin overlapping the Hall A Kin363 point.

    The selection is intentionally done using exact bin edges with a small
    tolerance. This avoids accidentally selecting nearby bins.
    """

    mask = (
        float_close(df["xBmin"], TARGET_XB_MIN)
        & float_close(df["xBmax"], TARGET_XB_MAX)
        & float_close(df["Q2min"], TARGET_Q2_MIN)
        & float_close(df["Q2max"], TARGET_Q2_MAX)
        & float_close(df["t_abs_min"], TARGET_T_ABS_MIN)
        & float_close(df["t_abs_max"], TARGET_T_ABS_MAX)
    )

    selected = df.loc[mask].copy()

    if len(selected) == 0:
        raise RuntimeError(
            "Could not find the requested Hall A overlap bin in the CSV:\n"
            f"  xB  [{TARGET_XB_MIN}, {TARGET_XB_MAX}]\n"
            f"  Q2  [{TARGET_Q2_MIN}, {TARGET_Q2_MAX}]\n"
            f"  |t| [{TARGET_T_ABS_MIN}, {TARGET_T_ABS_MAX}]\n"
        )
    # endif

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

    return points


def hall_b_points_for_period(selected: pd.DataFrame, period: str) -> List[DataPoint]:
    """
    Extract Hall B / CLAS12 phi-dependent points for one period.

    Cross-section cells are assumed to be stored as:

      (value, stat, sys)

    The plotted error is:

      sqrt(stat^2 + sys^2)

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


def hall_a_by_phi(points: List[DataPoint]) -> dict:
    return {round(p.phi, 6): p for p in points}


def ratio_to_hall_a(
    hall_b_points: List[DataPoint],
    hall_a_points_input: List[DataPoint],
    phi_tolerance: float = 1.0e-4,
) -> List[DataPoint]:
    """
    Compute Hall B / Hall A ratios.

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
        best_a = None
        best_delta = float("inf")

        for a in hall_a_points_input:
            delta = abs(b.phi - a.phi)

            if delta < best_delta:
                best_delta = delta
                best_a = a
            # endif
        # endfor

        if best_a is None or best_delta > phi_tolerance:
            continue
        # endif

        if not np.isfinite(best_a.sigma) or best_a.sigma == 0.0:
            continue
        # endif

        ratio = b.sigma / best_a.sigma

        ratio_err_high = math.sqrt(
            (b.err_high / best_a.sigma) ** 2
            + (b.sigma * best_a.err_low / (best_a.sigma * best_a.sigma)) ** 2
        )

        ratio_err_low = math.sqrt(
            (b.err_low / best_a.sigma) ** 2
            + (b.sigma * best_a.err_high / (best_a.sigma * best_a.sigma)) ** 2
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


def make_plot(
    selected: pd.DataFrame,
    output_dir: str,
    output_name: str,
) -> None:
    hall_a = hall_a_points()

    hall_b_by_period = {
        period: hall_b_points_for_period(selected, period)
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
        ratios = ratio_to_hall_a(
            hall_b_points=hall_b_by_period.get(label, []),
            hall_a_points_input=hall_a,
        )

        plot_dataset(right, ratios, label)
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

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, output_name)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote: {output_path}")


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

    return parser.parse_args()


def main() -> None:
    args = parse_args()

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

    required += [cross_section_column(period) for period in HALL_B_SERIES]

    df = pd.read_csv(args.csv_file)
    require_columns(df, required)

    selected = select_hall_b_overlap_bin(df)

    print("Selected Hall A overlap bin from CSV:")
    print(f"  rows: {len(selected)}")
    print(f"  xB:  [{TARGET_XB_MIN}, {TARGET_XB_MAX}]")
    print(f"  Q2:  [{TARGET_Q2_MIN}, {TARGET_Q2_MAX}] GeV^2")
    print(f"  |t|: [{TARGET_T_ABS_MIN}, {TARGET_T_ABS_MAX}] GeV^2")

    make_plot(
        selected=selected,
        output_dir=args.output_dir,
        output_name=args.output_name,
    )


if __name__ == "__main__":
    main()
# endif