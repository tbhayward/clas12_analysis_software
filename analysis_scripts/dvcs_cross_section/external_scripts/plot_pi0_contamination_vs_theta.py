#!/usr/bin/env python3
"""
plot_pi0_contamination_vs_theta.py

Make a 1x3 diagnostic canvas showing pi0 contamination ratio as a function of:

  electron polar angle: theta_e
  proton polar angle:   theta_p
  photon polar angle:   theta_gamma

for the five RGA run periods:

  Fa18 Inb
  Fa18 Out
  Sp19 Inb
  Sp18 Inb
  Sp18 Out

Input CSV columns expected:

  contamination ratio, Fa18 Inb
  contamination ratio, Fa18 Out
  contamination ratio, Sp19 Inb
  contamination ratio, Sp18 Inb
  contamination ratio, Sp18 Out

and angle-average columns used for binning/plot placement:

  e_theta, <period>
  p_theta, <period>
  g_theta, <period>

By default, the theta-bin edges are derived from:

  e_theta, 10.6 GeV
  p_theta, 10.6 GeV
  g_theta, 10.6 GeV

using seven equal-width bins. Change this with:

  --theta-binning-period "10.6 GeV"
  --theta-bins 7

The plotted point for each run period is the unweighted average contamination
ratio inside each derived theta bin. The error bar is the standard error of the
mean over the CSV rows in that theta bin.
"""

import argparse
import math
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DEFAULT_OUTPUT_DIR = "output/theta_diagnostics"
DEFAULT_OUTPUT_NAME = "pi0_contamination_vs_theta.png"

DEFAULT_THETA_BINS = 7
DEFAULT_THETA_BINNING_PERIOD = "10.6 GeV"

RUN_PERIODS = [
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out",
]

PERIOD_STYLE = {
    "Fa18 Inb": dict(marker="^", linestyle="-", color="tab:blue"),
    "Fa18 Out": dict(marker="v", linestyle="-", color="tab:orange"),
    "Sp19 Inb": dict(marker="s", linestyle="-", color="tab:red"),
    "Sp18 Inb": dict(marker="D", linestyle="-", color="tab:green"),
    "Sp18 Out": dict(marker="P", linestyle="-", color="tab:purple"),
}

THETA_PROJECTIONS = {
    "e_theta": {
        "theta_prefix": "e_theta",
        "xlabel": r"$\theta_{e}$  (deg)",
        "title": r"electron polar angle",
    },
    "p_theta": {
        "theta_prefix": "p_theta",
        "xlabel": r"$\theta_{p}$  (deg)",
        "title": r"proton polar angle",
    },
    "g_theta": {
        "theta_prefix": "g_theta",
        "xlabel": r"$\theta_{\gamma}$  (deg)",
        "title": r"photon polar angle",
    },
}

FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)


@dataclass
class BinnedPoint:
    x: float
    y: float
    yerr: float
    n: int


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


def metric_column(period: str) -> str:
    return f"contamination ratio, {period}"


def theta_column(theta_prefix: str, period: str) -> str:
    return f"{theta_prefix}, {period}"


def require_columns(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [col for col in columns if col not in df.columns]

    if missing:
        raise KeyError(
            "The input CSV is missing required columns:\n  "
            + "\n  ".join(missing)
        )
    # endif


def finite_column_values(df: pd.DataFrame, column: str) -> np.ndarray:
    values = np.array(
        [
            parse_scalar_from_cell(value)
            for value in df[column].values
        ],
        dtype=float,
    )

    return values[np.isfinite(values)]


def make_theta_edges(
    df: pd.DataFrame,
    theta_prefix: str,
    theta_binning_period: str,
    theta_bins: int,
) -> np.ndarray:
    source_col = theta_column(theta_prefix, theta_binning_period)
    values = finite_column_values(df, source_col)

    if len(values) == 0:
        raise RuntimeError(f"No finite values found in {source_col}.")
    # endif

    vmin = float(np.min(values))
    vmax = float(np.max(values))

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        raise RuntimeError(
            f"Invalid theta range for {source_col}: min={vmin}, max={vmax}."
        )
    # endif

    return np.linspace(vmin, vmax, theta_bins + 1)


def bin_index_for_value(value: float, edges: np.ndarray) -> int:
    if not np.isfinite(value):
        return -1
    # endif

    idx = int(np.searchsorted(edges, value, side="right") - 1)

    if idx < 0:
        idx = 0
    # endif

    if idx >= len(edges) - 1:
        idx = len(edges) - 2
    # endif

    return idx


def binned_metric_points(
    df: pd.DataFrame,
    theta_prefix: str,
    period: str,
    theta_edges: np.ndarray,
) -> List[BinnedPoint]:
    theta_col = theta_column(theta_prefix, period)
    y_col = metric_column(period)

    points: List[BinnedPoint] = []

    values_by_bin: Dict[int, List[float]] = {
        i: []
        for i in range(len(theta_edges) - 1)
    }

    theta_values_by_bin: Dict[int, List[float]] = {
        i: []
        for i in range(len(theta_edges) - 1)
    }

    for _, row in df.iterrows():
        theta_value = parse_scalar_from_cell(row[theta_col])
        metric_value = parse_scalar_from_cell(row[y_col])

        if not np.isfinite(theta_value) or not np.isfinite(metric_value):
            continue
        # endif

        bin_index = bin_index_for_value(theta_value, theta_edges)

        if bin_index < 0:
            continue
        # endif

        values_by_bin[bin_index].append(float(metric_value))
        theta_values_by_bin[bin_index].append(float(theta_value))
    # endfor

    for bin_index in range(len(theta_edges) - 1):
        values = np.array(values_by_bin[bin_index], dtype=float)
        theta_values = np.array(theta_values_by_bin[bin_index], dtype=float)

        if len(values) == 0:
            continue
        # endif

        x = float(np.mean(theta_values))
        y = float(np.mean(values))

        if len(values) > 1:
            yerr = float(np.std(values, ddof=1) / math.sqrt(len(values)))
        else:
            yerr = 0.0
        # endif

        points.append(
            BinnedPoint(
                x=x,
                y=y,
                yerr=yerr,
                n=len(values),
            )
        )
    # endfor

    points.sort(key=lambda p: p.x)

    return points


def plot_binned_points(ax, points: List[BinnedPoint], label: str) -> None:
    if len(points) == 0:
        return
    # endif

    x = np.array([p.x for p in points], dtype=float)
    y = np.array([p.y for p in points], dtype=float)
    yerr = np.array([p.yerr for p in points], dtype=float)

    style = PERIOD_STYLE.get(label, dict(marker="o", linestyle="-"))

    ax.errorbar(
        x,
        y,
        yerr=yerr,
        label=label,
        markersize=5.5,
        linewidth=1.2,
        elinewidth=1.0,
        capsize=2.5,
        **style,
    )


def make_plot(
    df: pd.DataFrame,
    output_dir: str,
    output_name: str,
    theta_bins: int,
    theta_binning_period: str,
) -> None:
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18.0, 5.2),
        constrained_layout=True,
    )

    fig.suptitle(r"$\pi^{0}$ contamination ratio vs detector polar angle", fontsize=16)

    for ax, projection_key in zip(axes, ["e_theta", "p_theta", "g_theta"]):
        info = THETA_PROJECTIONS[projection_key]
        theta_prefix = str(info["theta_prefix"])

        theta_edges = make_theta_edges(
            df=df,
            theta_prefix=theta_prefix,
            theta_binning_period=theta_binning_period,
            theta_bins=theta_bins,
        )

        print(
            f"{projection_key}: {theta_bins} bins from "
            f"{theta_edges[0]:.6g} to {theta_edges[-1]:.6g} deg "
            f"using {theta_column(theta_prefix, theta_binning_period)}"
        )

        for period in RUN_PERIODS:
            points = binned_metric_points(
                df=df,
                theta_prefix=theta_prefix,
                period=period,
                theta_edges=theta_edges,
            )

            plot_binned_points(
                ax=ax,
                points=points,
                label=period,
            )
        # endfor

        ax.set_xlabel(str(info["xlabel"]))
        ax.set_ylabel(r"$\pi^{0}$ contamination ratio")
        ax.set_title(str(info["title"]))
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8, frameon=True)
    # endfor

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_name)

    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote: {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot DVCS pi0 contamination ratio as a function of detector polar angles."
    )

    parser.add_argument(
        "csv_file",
        help="Input final DVCS pass-2 analysis CSV.",
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--output-name",
        default=DEFAULT_OUTPUT_NAME,
        help=f"Output PNG filename. Default: {DEFAULT_OUTPUT_NAME}",
    )

    parser.add_argument(
        "--theta-bins",
        type=int,
        default=DEFAULT_THETA_BINS,
        help=f"Number of derived theta bins. Default: {DEFAULT_THETA_BINS}.",
    )

    parser.add_argument(
        "--theta-binning-period",
        default=DEFAULT_THETA_BINNING_PERIOD,
        help=(
            "Period whose theta-average columns define the common theta-bin edges. "
            f"Default: '{DEFAULT_THETA_BINNING_PERIOD}'."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.theta_bins <= 0:
        raise ValueError("--theta-bins must be positive.")
    # endif

    df = pd.read_csv(args.csv_file)

    required = []

    for period in RUN_PERIODS:
        required.append(metric_column(period))

        for theta_prefix in ["e_theta", "p_theta", "g_theta"]:
            required.append(theta_column(theta_prefix, period))
        # endfor
    # endfor

    for theta_prefix in ["e_theta", "p_theta", "g_theta"]:
        required.append(theta_column(theta_prefix, args.theta_binning_period))
    # endfor

    require_columns(df, required)

    make_plot(
        df=df,
        output_dir=args.output_dir,
        output_name=args.output_name,
        theta_bins=args.theta_bins,
        theta_binning_period=args.theta_binning_period,
    )


if __name__ == "__main__":
    main()
# endif