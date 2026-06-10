#!/usr/bin/env python3
"""
integrated_cross_section_study.py

Standalone projection/integration study for the DVCS pass-2 CSV.

The script reads the final analysis CSV and makes four 1x3 canvases in
output/integrated_study/:

  integrated_xB_dependence.png/.pdf
  integrated_Q2_dependence.png/.pdf
  integrated_t_dependence.png/.pdf
  integrated_phi_dependence.png/.pdf

Each canvas contains:

  left:   integrated normed cross sections for 10.6 GeV, Sp19 Inb, and the
          four individual 10.6-GeV run periods.
  middle: ratios of the four individual 10.6-GeV run periods to the combined
          10.6-GeV result.
  right:  ratios of Fa18 Inb and Sp19 Inb to the weighted mean of Fa18 Inb
          and Sp19 Inb.

By default the integration is bin-width weighted, appropriate for differential
cross sections:

  d sigma / dxB       = sum_i sigma_i * dQ2_i * dt_i * dphi_i
  d sigma / dQ2       = sum_i sigma_i * dxB_i * dt_i * dphi_i
  d sigma / d|t|      = sum_i sigma_i * dxB_i * dQ2_i * dphi_i
  d sigma / dphi      = sum_i sigma_i * dxB_i * dQ2_i * dt_i

The phi width is converted from degrees to radians by default, since the usual
DVCS differential cross section convention is per radian. Use --phi-degrees if
you intentionally want degree-weighted phi integration instead.
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


RUN_PERIODS_106 = [
    "Fa18 Inb",
    "Fa18 Out",
    "Sp18 Inb",
    "Sp18 Out",
]

LEFT_SERIES = [
    "10.6 GeV",
    "Sp19 Inb",
    "Fa18 Inb",
    "Fa18 Out",
    "Sp18 Inb",
    "Sp18 Out",
]

MIDDLE_SERIES = RUN_PERIODS_106
RIGHT_SERIES = ["Fa18 Inb", "Sp19 Inb"]

SERIES_STYLE = {
    "10.6 GeV": dict(marker="o", linestyle="-",  color="black"),
    "Sp19 Inb": dict(marker="s", linestyle="-", color="tab:red"),
    "Fa18 Inb": dict(marker="^", linestyle="-", color="tab:blue"),
    "Fa18 Out": dict(marker="v", linestyle="-", color="tab:orange"),
    "Sp18 Inb": dict(marker="D", linestyle="-", color="tab:green"),
    "Sp18 Out": dict(marker="P", linestyle="-", color="tab:purple"),
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

WIDTH_COLUMNS = {
    "xB": ("xBmin", "xBmax"),
    "Q2": ("Q2min", "Q2max"),
    "t": ("t_abs_min", "t_abs_max"),
    "phi": ("phimin", "phimax"),
}


@dataclass
class Point:
    key: Tuple[float, float]
    x: float
    xerr: float
    y: float
    stat: float
    sys: float


FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)


def parse_tuple3(value) -> Tuple[float, float, float]:
    """Parse a CSV cell holding either NaN, scalar, or '(value, stat, sys)'."""

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

    out["_width_xB"] = pd.to_numeric(out["xBmax"], errors="coerce") - pd.to_numeric(out["xBmin"], errors="coerce")
    out["_width_Q2"] = pd.to_numeric(out["Q2max"], errors="coerce") - pd.to_numeric(out["Q2min"], errors="coerce")
    out["_width_t"] = pd.to_numeric(out["t_abs_max"], errors="coerce") - pd.to_numeric(out["t_abs_min"], errors="coerce")

    phi_width_deg = pd.to_numeric(out["phimax"], errors="coerce") - pd.to_numeric(out["phimin"], errors="coerce")
    if phi_degrees:
        out["_width_phi"] = phi_width_deg
    else:
        out["_width_phi"] = np.deg2rad(phi_width_deg)
    # endif

    return out


def weight_for_row(row: pd.Series, integrate_widths: List[str], no_width_weighting: bool) -> float:
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


def average_x_for_group(
    group: pd.DataFrame,
    projection_info: Dict[str, object],
    period: str,
) -> Tuple[float, float]:
    avg_col = f"{projection_info['avg_prefix']}, {period}"
    min_col = str(projection_info["min_col"])
    max_col = str(projection_info["max_col"])

    xlo = float(pd.to_numeric(group[min_col], errors="coerce").iloc[0])
    xhi = float(pd.to_numeric(group[max_col], errors="coerce").iloc[0])
    xcenter = 0.5 * (xlo + xhi)

    if avg_col in group.columns:
        vals = pd.to_numeric(group[avg_col], errors="coerce")
        finite = vals[np.isfinite(vals)]
        if len(finite) > 0:
            xcenter = float(np.mean(finite))
        # endif
    # endif

    return xcenter, 0.5 * (xhi - xlo)


def integrated_points_for_period(
    df: pd.DataFrame,
    projection: str,
    period: str,
    no_width_weighting: bool,
) -> List[Point]:
    info = PROJECTIONS[projection]
    col = cross_section_column(period)
    if col not in df.columns:
        return []
    # endif

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

            weight = weight_for_row(row, list(info["integrate_widths"]), no_width_weighting)
            if not np.isfinite(weight):
                continue
            # endif

            stat = 0.0 if not np.isfinite(stat) else stat
            sys = 0.0 if not np.isfinite(sys) else sys

            y_sum += value * weight
            stat2_sum += (stat * weight) ** 2
            sys2_sum += (sys * weight) ** 2
            n_used += 1
        # endfor

        if n_used == 0:
            continue
        # endif

        x, xerr = average_x_for_group(group, info, period)
        key = (float(group[str(info["min_col"])].iloc[0]), float(group[str(info["max_col"])].iloc[0]))
        points.append(Point(key=key, x=x, xerr=xerr, y=y_sum, stat=math.sqrt(stat2_sum), sys=math.sqrt(sys2_sum)))
    # endfor

    points.sort(key=lambda p: p.x)
    return points


def points_to_arrays(points: List[Point]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.array([p.x for p in points], dtype=float),
        np.array([p.xerr for p in points], dtype=float),
        np.array([p.y for p in points], dtype=float),
        np.array([p.stat for p in points], dtype=float),
        np.array([p.sys for p in points], dtype=float),
    )


def ratio_points(numerator: List[Point], denominator: List[Point]) -> List[Point]:
    n_by_x = {p.key: p for p in numerator}
    d_by_x = {p.key: p for p in denominator}

    ratios: List[Point] = []
    common_keys = sorted(set(n_by_x.keys()) & set(d_by_x.keys()))

    for key in common_keys:
        n = n_by_x[key]
        d = d_by_x[key]
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

        ratios.append(Point(key=n.key, x=n.x, xerr=n.xerr, y=r, stat=abs(r) * math.sqrt(rel2), sys=0.0))
    # endfor

    return ratios


def weighted_mean_two_points(a: Point, b: Point) -> Optional[Point]:
    values = []
    weights = []

    for p in [a, b]:
        if not np.isfinite(p.y):
            continue
        # endif
        if np.isfinite(p.stat) and p.stat > 0.0:
            w = 1.0 / (p.stat * p.stat)
        else:
            w = 0.0
        # endif
        if w > 0.0:
            values.append(p.y)
            weights.append(w)
        # endif
    # endfor

    if len(values) == 0:
        return None
    # endif

    values_arr = np.array(values, dtype=float)
    weights_arr = np.array(weights, dtype=float)
    mean = float(np.sum(weights_arr * values_arr) / np.sum(weights_arr))
    stat = float(math.sqrt(1.0 / np.sum(weights_arr)))
    return Point(key=a.key, x=a.x, xerr=a.xerr, y=mean, stat=stat, sys=0.0)


def ratio_to_fa18_sp19_weighted_mean(
    points_by_period: Dict[str, List[Point]],
    period: str,
) -> List[Point]:
    fa_points = {p.key: p for p in points_by_period.get("Fa18 Inb", [])}
    sp_points = {p.key: p for p in points_by_period.get("Sp19 Inb", [])}
    target_points = {p.key: p for p in points_by_period.get(period, [])}

    ratios: List[Point] = []
    common_keys = sorted(set(fa_points.keys()) & set(sp_points.keys()) & set(target_points.keys()))

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

        ratios.append(Point(key=p.key, x=p.x, xerr=p.xerr, y=r, stat=abs(r) * math.sqrt(rel2), sys=0.0))
    # endfor

    return ratios


def plot_points(ax, points: List[Point], label: str, ratio: bool = False) -> None:
    if len(points) == 0:
        return
    # endif

    x, xerr, y, stat, _ = points_to_arrays(points)
    style = SERIES_STYLE.get(label, dict(marker="o", linestyle="-"))

    ax.errorbar(
        x,
        y,
        xerr=xerr,
        yerr=stat,
        label=label,
        markersize=5.0,
        linewidth=1.2,
        elinewidth=1.0,
        capsize=2.5,
        **style,
    )

    if ratio:
        ax.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--", zorder=0)
    # endif


def auto_ratio_ylim(ax, center: float = 1.0, min_span: float = 0.35) -> None:
    ylo, yhi = ax.get_ylim()
    span = max(abs(yhi - center), abs(center - ylo), min_span)
    ax.set_ylim(center - 1.15 * span, center + 1.15 * span)


def make_projection_plot(
    df: pd.DataFrame,
    projection: str,
    output_dir: str,
    no_width_weighting: bool,
    show_sys_band: bool,
) -> None:
    info = PROJECTIONS[projection]

    all_needed_periods = sorted(set(LEFT_SERIES + MIDDLE_SERIES + RIGHT_SERIES))
    points_by_period = {
        period: integrated_points_for_period(df, projection, period, no_width_weighting)
        for period in all_needed_periods
    }

    fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.5), constrained_layout=True)
    fig.suptitle(str(info["title"]), fontsize=16)

    left = axes[0]
    for period in LEFT_SERIES:
        plot_points(left, points_by_period.get(period, []), period, ratio=False)
    # endfor
    left.set_xlabel(str(info["xlabel"]))
    left.set_ylabel(str(info["ylabel"]))
    left.set_title("Integrated cross sections")
    left.grid(True, alpha=0.25)
    left.legend(fontsize=8, frameon=True)

    middle = axes[1]
    denom_106 = points_by_period.get("10.6 GeV", [])
    for period in MIDDLE_SERIES:
        plot_points(middle, ratio_points(points_by_period.get(period, []), denom_106), period, ratio=True)
    # endfor
    middle.set_xlabel(str(info["xlabel"]))
    middle.set_ylabel(r"run period / 10.6 GeV")
    middle.set_title("10.6-GeV period consistency")
    middle.grid(True, alpha=0.25)
    middle.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(middle)

    right = axes[2]
    for period in RIGHT_SERIES:
        plot_points(
            right,
            ratio_to_fa18_sp19_weighted_mean(points_by_period, period),
            period,
            ratio=True,
        )
    # endfor
    right.set_xlabel(str(info["xlabel"]))
    right.set_ylabel(r"period / weighted mean")
    right.set_title("Fa18 Inb vs Sp19 Inb")
    right.grid(True, alpha=0.25)
    right.legend(fontsize=8, frameon=True)
    auto_ratio_ylim(right)

    if show_sys_band:
        # This hook is intentionally left simple: current ratio panels use stat errors
        # only, while the left panel can already be inspected with statistical error bars.
        # Add filled systematic bands here if desired later.
        pass
    # endif

    outbase = os.path.join(output_dir, f"integrated_{info['outfile_tag']}_dependence")
    fig.savefig(outbase + ".png", dpi=200)
    fig.savefig(outbase + ".pdf")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make integrated xB, Q2, |t|, and phi dependence plots from a final DVCS pass-2 CSV."
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
        "--show-sys-band",
        action="store_true",
        help="Reserved option for future systematic-band plotting; current plots use stat error bars.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    required = [
        "xBmin", "xBmax",
        "Q2min", "Q2max",
        "t_abs_min", "t_abs_max",
        "phimin", "phimax",
    ]
    required += [cross_section_column(period) for period in LEFT_SERIES]

    df = pd.read_csv(args.csv_file)
    require_columns(df, required)
    df = add_width_columns(df, phi_degrees=args.phi_degrees)

    os.makedirs(args.output_dir, exist_ok=True)

    for projection in ["xB", "Q2", "t", "phi"]:
        make_projection_plot(
            df=df,
            projection=projection,
            output_dir=args.output_dir,
            no_width_weighting=args.no_width_weighting,
            show_sys_band=args.show_sys_band,
        )
    # endfor

    print(f"Wrote integrated study plots to: {args.output_dir}")


if __name__ == "__main__":
    main()
# endif
