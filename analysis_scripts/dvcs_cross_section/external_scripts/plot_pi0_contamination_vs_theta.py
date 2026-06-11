#!/usr/bin/env python3
"""
plot_pi0_contamination_vs_theta.py

Fast diagnostic plot of pi0 contamination ratio vs detector polar angle.

Makes one 1x3 canvas:

  output/theta_diagnostics/pi0_contamination_vs_theta.png

Panels:
  electron polar angle: theta_e
  proton polar angle:   theta_p
  photon polar angle:   theta_gamma

Run periods:
  Fa18 Inb
  Fa18 Out
  Sp19 Inb
  Sp18 Inb
  Sp18 Out

Input metric columns:
  contamination ratio, Fa18 Inb
  contamination ratio, Fa18 Out
  contamination ratio, Sp19 Inb
  contamination ratio, Sp18 Inb
  contamination ratio, Sp18 Out

Input angle-average columns:
  e_theta, <period>
  p_theta, <period>
  g_theta, <period>

The theta-bin edges are derived from:
  e_theta, 10.6 GeV
  p_theta, 10.6 GeV
  g_theta, 10.6 GeV

by default, with seven equal-width bins. Change with:
  --theta-bins 7
  --theta-binning-period "10.6 GeV"

The plotted point for each run period is the mean contamination ratio inside
each derived theta bin. The error bar is the standard error of the mean over CSV
rows in that theta bin.

Speed notes:
  - Reads only required CSV columns.
  - Uses vectorized numpy/pandas binning instead of row loops.
  - Computes the three theta projections in parallel by default.
"""

import argparse
import concurrent.futures
import math
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DEFAULT_OUTPUT_DIR = "output/theta_diagnostics"
DEFAULT_OUTPUT_NAME = "pi0_contamination_vs_theta.png"

DEFAULT_THETA_BINS = 7
DEFAULT_THETA_BINNING_PERIOD = "10.6 GeV"
DEFAULT_MAX_WORKERS = 3

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

THETA_ORDER = [
    "e_theta",
    "p_theta",
    "g_theta",
]

FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)


@dataclass
class BinnedPoint:
    x: float
    y: float
    yerr: float
    n: int


@dataclass
class ProjectionResult:
    projection_key: str
    theta_edges: np.ndarray
    points_by_period: Dict[str, List[BinnedPoint]]


def metric_column(period: str) -> str:
    return f"contamination ratio, {period}"


def theta_column(theta_prefix: str, period: str) -> str:
    return f"{theta_prefix}, {period}"


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


def series_to_float(series: pd.Series) -> pd.Series:
    """
    Fast conversion for mostly-numeric columns, with regex fallback for tuple-like
    or string cells.
    """

    numeric = pd.to_numeric(series, errors="coerce")

    if numeric.notna().sum() == series.notna().sum():
        return numeric.astype(float)
    # endif

    return series.map(parse_scalar_from_cell).astype(float)


def require_columns_from_header(csv_file: str, columns: Iterable[str]) -> None:
    header = pd.read_csv(csv_file, nrows=0)
    existing = set(header.columns)
    missing = [col for col in columns if col not in existing]

    if missing:
        raise KeyError(
            "The input CSV is missing required columns:\n  "
            + "\n  ".join(missing)
        )
    # endif


def required_columns(theta_binning_period: str) -> List[str]:
    cols: List[str] = []

    for period in RUN_PERIODS:
        cols.append(metric_column(period))

        for theta_prefix in ["e_theta", "p_theta", "g_theta"]:
            cols.append(theta_column(theta_prefix, period))
        # endfor
    # endfor

    for theta_prefix in ["e_theta", "p_theta", "g_theta"]:
        cols.append(theta_column(theta_prefix, theta_binning_period))
    # endfor

    return sorted(set(cols))


def read_reduced_csv(csv_file: str, theta_binning_period: str) -> pd.DataFrame:
    cols = required_columns(theta_binning_period)
    require_columns_from_header(csv_file, cols)

    return pd.read_csv(csv_file, usecols=cols)


def make_theta_edges(
    df: pd.DataFrame,
    theta_prefix: str,
    theta_binning_period: str,
    theta_bins: int,
) -> np.ndarray:
    source_col = theta_column(theta_prefix, theta_binning_period)
    values = series_to_float(df[source_col]).to_numpy(dtype=float)
    finite = values[np.isfinite(values)]

    if len(finite) == 0:
        raise RuntimeError(f"No finite values found in {source_col}.")
    # endif

    vmin = float(np.min(finite))
    vmax = float(np.max(finite))

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        raise RuntimeError(
            f"Invalid theta range for {source_col}: min={vmin}, max={vmax}."
        )
    # endif

    return np.linspace(vmin, vmax, theta_bins + 1)


def binned_metric_points_vectorized(
    theta_values: np.ndarray,
    metric_values: np.ndarray,
    theta_edges: np.ndarray,
) -> List[BinnedPoint]:
    finite_mask = np.isfinite(theta_values) & np.isfinite(metric_values)

    if not np.any(finite_mask):
        return []
    # endif

    theta = theta_values[finite_mask]
    metric = metric_values[finite_mask]

    bin_indices = np.searchsorted(theta_edges, theta, side="right") - 1
    bin_indices = np.clip(bin_indices, 0, len(theta_edges) - 2)

    points: List[BinnedPoint] = []

    for bin_index in range(len(theta_edges) - 1):
        mask = bin_indices == bin_index

        if not np.any(mask):
            continue
        # endif

        theta_bin = theta[mask]
        metric_bin = metric[mask]

        n = int(metric_bin.size)
        x = float(np.mean(theta_bin))
        y = float(np.mean(metric_bin))

        if n > 1:
            yerr = float(np.std(metric_bin, ddof=1) / math.sqrt(n))
        else:
            yerr = 0.0
        # endif

        points.append(
            BinnedPoint(
                x=x,
                y=y,
                yerr=yerr,
                n=n,
            )
        )
    # endfor

    points.sort(key=lambda p: p.x)

    return points


def compute_projection_result(
    args_tuple: Tuple[pd.DataFrame, str, int, str],
) -> ProjectionResult:
    df, projection_key, theta_bins, theta_binning_period = args_tuple

    info = THETA_PROJECTIONS[projection_key]
    theta_prefix = str(info["theta_prefix"])

    theta_edges = make_theta_edges(
        df=df,
        theta_prefix=theta_prefix,
        theta_binning_period=theta_binning_period,
        theta_bins=theta_bins,
    )

    points_by_period: Dict[str, List[BinnedPoint]] = {}

    metric_cache: Dict[str, np.ndarray] = {}

    for period in RUN_PERIODS:
        theta_col = theta_column(theta_prefix, period)
        metric_col = metric_column(period)

        theta_values = series_to_float(df[theta_col]).to_numpy(dtype=float)

        if metric_col not in metric_cache:
            metric_cache[metric_col] = series_to_float(df[metric_col]).to_numpy(dtype=float)
        # endif

        metric_values = metric_cache[metric_col]

        points_by_period[period] = binned_metric_points_vectorized(
            theta_values=theta_values,
            metric_values=metric_values,
            theta_edges=theta_edges,
        )
    # endfor

    return ProjectionResult(
        projection_key=projection_key,
        theta_edges=theta_edges,
        points_by_period=points_by_period,
    )


def compute_all_projection_results(
    df: pd.DataFrame,
    theta_bins: int,
    theta_binning_period: str,
    max_workers: int,
    no_parallel: bool,
) -> Dict[str, ProjectionResult]:
    tasks = [
        (df, projection_key, theta_bins, theta_binning_period)
        for projection_key in THETA_ORDER
    ]

    if no_parallel or max_workers <= 1:
        results = [compute_projection_result(task) for task in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(compute_projection_result, tasks))
        # endwith
    # endif

    return {
        result.projection_key: result
        for result in results
    }


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
    results_by_projection: Dict[str, ProjectionResult],
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

    for ax, projection_key in zip(axes, THETA_ORDER):
        info = THETA_PROJECTIONS[projection_key]
        result = results_by_projection[projection_key]
        theta_prefix = str(info["theta_prefix"])

        print(
            f"{projection_key}: {theta_bins} bins from "
            f"{result.theta_edges[0]:.6g} to {result.theta_edges[-1]:.6g} deg "
            f"using {theta_column(theta_prefix, theta_binning_period)}"
        )

        for period in RUN_PERIODS:
            plot_binned_points(
                ax=ax,
                points=result.points_by_period.get(period, []),
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

    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help=f"Maximum projection workers. Default: {DEFAULT_MAX_WORKERS}.",
    )

    parser.add_argument(
        "--no-parallel",
        action="store_true",
        help="Disable parallel computation of the three theta projections.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.theta_bins <= 0:
        raise ValueError("--theta-bins must be positive.")
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

    df = read_reduced_csv(
        csv_file=args.csv_file,
        theta_binning_period=args.theta_binning_period,
    )

    print(f"Read reduced CSV with {len(df)} rows and {len(df.columns)} columns.")
    print(f"Parallel theta projections: {not args.no_parallel and args.max_workers > 1}")
    print(f"Max workers: {args.max_workers}")

    results_by_projection = compute_all_projection_results(
        df=df,
        theta_bins=args.theta_bins,
        theta_binning_period=args.theta_binning_period,
        max_workers=args.max_workers,
        no_parallel=args.no_parallel,
    )

    make_plot(
        results_by_projection=results_by_projection,
        output_dir=args.output_dir,
        output_name=args.output_name,
        theta_bins=args.theta_bins,
        theta_binning_period=args.theta_binning_period,
    )


if __name__ == "__main__":
    main()
# endif