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

The plotted point for each run period is the mean pi0 contamination ratio inside
each derived theta bin. The error bar is the standard error of the mean over CSV
rows in that theta bin.

This version intentionally avoids multiprocessing. For this task, multiprocessing
is usually slower because it copies large pandas objects between processes.
"""

import argparse
import math
import os
import re
import time
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

THETA_ORDER = [
    "e_theta",
    "p_theta",
    "g_theta",
]

FLOAT_PATTERN = re.compile(
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
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


def log(message: str) -> None:
    now = time.strftime("%H:%M:%S")
    print(f"[{now}] {message}", flush=True)


class Timer:
    def __init__(self, label: str):
        self.label = label
        self.t0 = 0.0

    def __enter__(self):
        self.t0 = time.perf_counter()
        log(f"START: {self.label}")
        return self

    def __exit__(self, exc_type, exc, tb):
        dt = time.perf_counter() - self.t0
        log(f"DONE:  {self.label}  ({dt:.3f} s)")
        return False


def metric_column(period: str) -> str:
    return f"contamination ratio, {period}"


def theta_column(theta_prefix: str, period: str) -> str:
    return f"{theta_prefix}, {period}"


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


def require_columns_from_header(csv_file: str, columns: Iterable[str]) -> None:
    with Timer("reading CSV header"):
        header = pd.read_csv(csv_file, nrows=0)
    # endwith

    existing = set(header.columns)
    missing = [col for col in columns if col not in existing]

    log(f"CSV header contains {len(header.columns)} columns")
    log(f"Script requires {len(list(columns))} unique columns")

    if missing:
        raise KeyError(
            "The input CSV is missing required columns:\n  "
            + "\n  ".join(missing)
        )
    # endif


def read_reduced_csv(csv_file: str, theta_binning_period: str) -> pd.DataFrame:
    cols = required_columns(theta_binning_period)
    require_columns_from_header(csv_file, cols)

    with Timer(f"reading reduced CSV with usecols={len(cols)}"):
        df = pd.read_csv(csv_file, usecols=cols, low_memory=False)
    # endwith

    memory_mb = df.memory_usage(deep=True).sum() / 1024.0 / 1024.0
    log(f"Reduced dataframe shape: {df.shape[0]} rows x {df.shape[1]} columns")
    log(f"Reduced dataframe memory: {memory_mb:.2f} MB")

    return df


def numeric_array_from_series(series: pd.Series, column_name: str) -> np.ndarray:
    """
    Convert a pandas Series to a float numpy array.

    First try fast pd.to_numeric. If that gives almost no finite values, use a
    vectorized regex extraction of the first number in the string. This avoids
    Python-level per-row loops.
    """

    with Timer(f"converting column to float: {column_name}"):
        numeric = pd.to_numeric(series, errors="coerce")
        arr = numeric.to_numpy(dtype=float)
        finite_count = int(np.isfinite(arr).sum())

        if finite_count > 0:
            log(
                f"  {column_name}: pd.to_numeric finite count = "
                f"{finite_count}/{len(arr)}"
            )
            return arr
        # endif

        log(
            f"  {column_name}: pd.to_numeric found zero finite values; "
            "trying vectorized regex extraction"
        )

        extracted = series.astype(str).str.extract(FLOAT_PATTERN, expand=False)
        numeric = pd.to_numeric(extracted, errors="coerce")
        arr = numeric.to_numpy(dtype=float)
        finite_count = int(np.isfinite(arr).sum())

        log(
            f"  {column_name}: regex finite count = "
            f"{finite_count}/{len(arr)}"
        )

        return arr
    # endwith


def build_numeric_cache(df: pd.DataFrame) -> Dict[str, np.ndarray]:
    cache: Dict[str, np.ndarray] = {}

    log("Building numeric array cache for all loaded columns")

    for i, col in enumerate(df.columns, start=1):
        log(f"Column {i}/{len(df.columns)}")
        cache[col] = numeric_array_from_series(df[col], col)
    # endfor

    return cache


def make_theta_edges(
    numeric_cache: Dict[str, np.ndarray],
    theta_prefix: str,
    theta_binning_period: str,
    theta_bins: int,
) -> np.ndarray:
    source_col = theta_column(theta_prefix, theta_binning_period)
    values = numeric_cache[source_col]
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

    edges = np.linspace(vmin, vmax, theta_bins + 1)

    log(
        f"  theta edges from {source_col}: "
        f"{theta_bins} bins, min={vmin:.6g}, max={vmax:.6g}"
    )

    return edges


def binned_metric_points_numpy(
    theta_values: np.ndarray,
    metric_values: np.ndarray,
    theta_edges: np.ndarray,
    debug_label: str,
) -> List[BinnedPoint]:
    finite_mask = np.isfinite(theta_values) & np.isfinite(metric_values)
    n_finite = int(np.count_nonzero(finite_mask))

    log(f"    {debug_label}: finite theta+metric rows = {n_finite}")

    if n_finite == 0:
        return []
    # endif

    theta = theta_values[finite_mask]
    metric = metric_values[finite_mask]

    bin_indices = np.searchsorted(theta_edges, theta, side="right") - 1
    bin_indices = np.clip(bin_indices, 0, len(theta_edges) - 2)

    nbins = len(theta_edges) - 1

    counts = np.bincount(bin_indices, minlength=nbins).astype(float)
    sum_metric = np.bincount(bin_indices, weights=metric, minlength=nbins)
    sum_metric2 = np.bincount(bin_indices, weights=metric * metric, minlength=nbins)
    sum_theta = np.bincount(bin_indices, weights=theta, minlength=nbins)

    points: List[BinnedPoint] = []

    for bin_index in range(nbins):
        n = int(counts[bin_index])

        if n <= 0:
            continue
        # endif

        x = float(sum_theta[bin_index] / counts[bin_index])
        y = float(sum_metric[bin_index] / counts[bin_index])

        if n > 1:
            variance_num = sum_metric2[bin_index] - (sum_metric[bin_index] * sum_metric[bin_index] / counts[bin_index])
            variance = max(0.0, variance_num / (counts[bin_index] - 1.0))
            yerr = float(math.sqrt(variance) / math.sqrt(counts[bin_index]))
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

    bin_count_text = ", ".join(str(p.n) for p in points)
    log(f"    {debug_label}: nonempty bins = {len(points)}; counts = [{bin_count_text}]")

    return points


def compute_projection_result(
    numeric_cache: Dict[str, np.ndarray],
    projection_key: str,
    theta_bins: int,
    theta_binning_period: str,
) -> ProjectionResult:
    with Timer(f"computing projection: {projection_key}"):
        info = THETA_PROJECTIONS[projection_key]
        theta_prefix = str(info["theta_prefix"])

        theta_edges = make_theta_edges(
            numeric_cache=numeric_cache,
            theta_prefix=theta_prefix,
            theta_binning_period=theta_binning_period,
            theta_bins=theta_bins,
        )

        points_by_period: Dict[str, List[BinnedPoint]] = {}

        for period in RUN_PERIODS:
            theta_col = theta_column(theta_prefix, period)
            metric_col = metric_column(period)

            log(f"  {projection_key}, {period}: theta_col={theta_col}, metric_col={metric_col}")

            points_by_period[period] = binned_metric_points_numpy(
                theta_values=numeric_cache[theta_col],
                metric_values=numeric_cache[metric_col],
                theta_edges=theta_edges,
                debug_label=f"{projection_key}, {period}",
            )
        # endfor

        return ProjectionResult(
            projection_key=projection_key,
            theta_edges=theta_edges,
            points_by_period=points_by_period,
        )
    # endwith


def compute_all_projection_results(
    numeric_cache: Dict[str, np.ndarray],
    theta_bins: int,
    theta_binning_period: str,
) -> Dict[str, ProjectionResult]:
    results: Dict[str, ProjectionResult] = {}

    for projection_key in THETA_ORDER:
        results[projection_key] = compute_projection_result(
            numeric_cache=numeric_cache,
            projection_key=projection_key,
            theta_bins=theta_bins,
            theta_binning_period=theta_binning_period,
        )
    # endfor

    return results


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
    with Timer("making matplotlib canvas"):
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

            log(
                f"{projection_key}: plotting {theta_bins} bins from "
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

        log(f"Wrote: {output_path}")
    # endwith


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
    total_t0 = time.perf_counter()
    args = parse_args()

    log("============================================================")
    log("plot_pi0_contamination_vs_theta.py")
    log("============================================================")
    log(f"Input CSV: {args.csv_file}")
    log(f"Output dir: {args.output_dir}")
    log(f"Output name: {args.output_name}")
    log(f"Theta bins: {args.theta_bins}")
    log(f"Theta binning period: {args.theta_binning_period}")

    if args.theta_bins <= 0:
        raise ValueError("--theta-bins must be positive.")
    # endif

    df = read_reduced_csv(
        csv_file=args.csv_file,
        theta_binning_period=args.theta_binning_period,
    )

    numeric_cache = build_numeric_cache(df)

    results_by_projection = compute_all_projection_results(
        numeric_cache=numeric_cache,
        theta_bins=args.theta_bins,
        theta_binning_period=args.theta_binning_period,
    )

    make_plot(
        results_by_projection=results_by_projection,
        output_dir=args.output_dir,
        output_name=args.output_name,
        theta_bins=args.theta_bins,
        theta_binning_period=args.theta_binning_period,
    )

    total_dt = time.perf_counter() - total_t0
    log(f"TOTAL RUNTIME: {total_dt:.3f} s")
    log("Done.")


if __name__ == "__main__":
    main()
# endif