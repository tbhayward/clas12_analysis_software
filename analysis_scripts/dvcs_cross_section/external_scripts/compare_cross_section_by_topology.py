#!/usr/bin/env python3
"""
compare_cross_sections_by_topology.py

Compare integrated DVCS cross sections from separate topology-resolved CSV files.

Expected usage:

  python compare_cross_sections_by_topology.py FD-FD.csv CD-FT.csv CD-FD.csv

The input CSVs should each have the usual pass-2 columns, for example:

  normed cross sections, ep->epg, exp, Fa18 Inb, unpol
  normed cross sections, ep->epg, exp, Fa18 Out, unpol
  normed cross sections, ep->epg, exp, Sp19 Inb, unpol
  normed cross sections, ep->epg, exp, Sp18 Inb, unpol
  normed cross sections, ep->epg, exp, Sp18 Out, unpol

The script makes one 2x3 canvas per projection:

  xB, Q2, |t|, phi, theta_e, theta_p, theta_gamma

Each canvas has one subplot per run period, with the three topologies overlaid.

Outputs:

  output/topology_comparison/topology_xB_comparison.png
  output/topology_comparison/topology_Q2_comparison.png
  output/topology_comparison/topology_t_comparison.png
  output/topology_comparison/topology_phi_comparison.png
  output/topology_comparison/topology_e_theta_comparison.png
  output/topology_comparison/topology_p_theta_comparison.png
  output/topology_comparison/topology_g_theta_comparison.png
"""

import argparse
import math
import os
import re
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg", force=True)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


RUN_PERIODS = ["Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"]

DEFAULT_LABELS = ["FD-FD", "CD-FT", "CD-FD"]
DEFAULT_OUTPUT_DIR = "output/topology_comparison"
DEFAULT_OUTPUT_PREFIX = "topology"
DEFAULT_THETA_BINS = 7
DEFAULT_THETA_BINNING_PERIOD = "10.6 GeV"
DEFAULT_XS_TEMPLATE = "normed cross sections, ep->epg, exp, {period}, unpol"

CATEGORY_STYLE = {
    "FD-FD": dict(marker="o", linestyle="-", color="tab:blue"),
    "CD-FT": dict(marker="s", linestyle="-", color="tab:orange"),
    "CD-FD": dict(marker="^", linestyle="-", color="tab:green"),
}

PROJECTIONS = {
    "xB": {
        "min_col": "xBmin", "max_col": "xBmax", "avg_prefix": "xBavg",
        "xlabel": r"$x_B$", "ylabel": r"$d\sigma/dx_B$  (pb)",
        "title": r"$x_B$ dependence", "integrate_widths": ["Q2", "t", "phi"],
        "outfile_tag": "xB", "derived_theta": False,
    },
    "Q2": {
        "min_col": "Q2min", "max_col": "Q2max", "avg_prefix": "Q2avg",
        "xlabel": r"$Q^2$  (GeV$^2$)", "ylabel": r"$d\sigma/dQ^2$  (pb/GeV$^2$)",
        "title": r"$Q^2$ dependence", "integrate_widths": ["xB", "t", "phi"],
        "outfile_tag": "Q2", "derived_theta": False,
    },
    "t": {
        "min_col": "t_abs_min", "max_col": "t_abs_max", "avg_prefix": "t_abs_avg",
        "xlabel": r"$|t|$  (GeV$^2$)", "ylabel": r"$d\sigma/d|t|$  (pb/GeV$^2$)",
        "title": r"$|t|$ dependence", "integrate_widths": ["xB", "Q2", "phi"],
        "outfile_tag": "t", "derived_theta": False,
    },
    "phi": {
        "min_col": "phimin", "max_col": "phimax", "avg_prefix": "phiavg",
        "xlabel": r"$\phi$  (deg)", "ylabel": r"$d\sigma/d\phi$  (pb/rad)",
        "title": r"$\phi$ dependence", "integrate_widths": ["xB", "Q2", "t"],
        "outfile_tag": "phi", "derived_theta": False,
    },
    "e_theta": {
        "min_col": "_bin_e_theta_min", "max_col": "_bin_e_theta_max",
        "avg_prefix": "e_theta", "theta_source_prefix": "e_theta",
        "xlabel": r"$\theta_{e}$  (deg)", "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"electron polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "e_theta", "derived_theta": True,
    },
    "p_theta": {
        "min_col": "_bin_p_theta_min", "max_col": "_bin_p_theta_max",
        "avg_prefix": "p_theta", "theta_source_prefix": "p_theta",
        "xlabel": r"$\theta_{p}$  (deg)", "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"proton polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "p_theta", "derived_theta": True,
    },
    "g_theta": {
        "min_col": "_bin_g_theta_min", "max_col": "_bin_g_theta_max",
        "avg_prefix": "g_theta", "theta_source_prefix": "g_theta",
        "xlabel": r"$\theta_{\gamma}$  (deg)", "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"photon polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "g_theta", "derived_theta": True,
    },
}

PROJECTION_ORDER = ["xB", "Q2", "t", "phi", "e_theta", "p_theta", "g_theta"]

FLOAT_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


@dataclass
class Point:
    key: Tuple[float, float]
    x: float
    y: float
    stat: float
    sys: float


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class Timer:
    def __init__(self, label: str):
        self.label = label
        self.t0 = 0.0

    def __enter__(self):
        self.t0 = time.perf_counter()
        log(f"START: {self.label}")
        return self

    def __exit__(self, exc_type, exc, tb):
        log(f"DONE:  {self.label}  ({time.perf_counter() - self.t0:.3f} s)")
        return False


def parse_tuple3(value) -> Tuple[float, float, float]:
    if value is None:
        return math.nan, math.nan, math.nan
    if isinstance(value, float) and math.isnan(value):
        return math.nan, math.nan, math.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value), 0.0, 0.0

    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null"}:
        return math.nan, math.nan, math.nan

    nums = FLOAT_PATTERN.findall(text)
    if not nums:
        return math.nan, math.nan, math.nan

    vals = [float(x) for x in nums]
    while len(vals) < 3:
        vals.append(0.0)
    return vals[0], vals[1], vals[2]


def parse_scalar(value) -> float:
    if value is None:
        return math.nan
    if isinstance(value, float) and math.isnan(value):
        return math.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)

    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null"}:
        return math.nan

    nums = FLOAT_PATTERN.findall(text)
    if not nums:
        return math.nan
    return float(nums[0])


def xs_col(period: str, template: str) -> str:
    return template.format(period=period)


def theta_col(prefix: str, period: str) -> str:
    return f"{prefix}, {period}"


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        log(f"Missing columns for {label}:")
        for c in missing:
            print(f"  {c}")
        relevant = [c for c in df.columns if "cross sections" in c.lower()]
        if relevant:
            log("Cross-section-like columns present:")
            for c in relevant[:80]:
                print(f"  {c}")
        raise KeyError(f"{label} is missing {len(missing)} required columns")


def read_inputs(paths: List[str], labels: List[str], template: str) -> Dict[str, pd.DataFrame]:
    required_base = [
        "xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max", "phimin", "phimax",
    ]

    for period in RUN_PERIODS:
        required_base += [
            f"xBavg, {period}", f"Q2avg, {period}", f"t_abs_avg, {period}",
            f"phiavg, {period}", f"e_theta, {period}", f"p_theta, {period}", f"g_theta, {period}",
            xs_col(period, template),
        ]

    out = {}
    for path, label in zip(paths, labels):
        with Timer(f"reading {label}: {path}"):
            df = pd.read_csv(path, low_memory=False)
        require_columns(df, required_base, label)
        out[label] = df
        log(f"{label}: {df.shape[0]} rows x {df.shape[1]} columns")
    return out


def add_width_columns(df: pd.DataFrame, phi_degrees: bool) -> pd.DataFrame:
    out = df.copy()
    out["_width_xB"] = pd.to_numeric(out["xBmax"], errors="coerce") - pd.to_numeric(out["xBmin"], errors="coerce")
    out["_width_Q2"] = pd.to_numeric(out["Q2max"], errors="coerce") - pd.to_numeric(out["Q2min"], errors="coerce")
    out["_width_t"] = pd.to_numeric(out["t_abs_max"], errors="coerce") - pd.to_numeric(out["t_abs_min"], errors="coerce")
    phi_width_deg = pd.to_numeric(out["phimax"], errors="coerce") - pd.to_numeric(out["phimin"], errors="coerce")
    out["_width_phi"] = phi_width_deg if phi_degrees else np.deg2rad(phi_width_deg)
    return out


def build_common_theta_edges(dfs: Dict[str, pd.DataFrame], theta_binning_period: str, theta_bins: int) -> Dict[str, np.ndarray]:
    edges_by_projection = {}

    for projection in ["e_theta", "p_theta", "g_theta"]:
        prefix = PROJECTIONS[projection]["theta_source_prefix"]
        col = theta_col(prefix, theta_binning_period)

        all_vals = []
        for label, df in dfs.items():
            if col not in df.columns:
                raise KeyError(f"{label} missing theta-binning column {col}")
            vals = df[col].apply(parse_scalar).astype(float).to_numpy()
            all_vals.append(vals[np.isfinite(vals)])

        merged = np.concatenate(all_vals) if all_vals else np.array([], dtype=float)
        if len(merged) == 0:
            raise RuntimeError(f"No finite theta values for {col}")

        vmin = float(np.min(merged))
        vmax = float(np.max(merged))
        edges_by_projection[projection] = np.linspace(vmin, vmax, theta_bins + 1)
        log(f"{projection}: common theta edges from {col}: {theta_bins} bins, {vmin:.6g} to {vmax:.6g} deg")

    return edges_by_projection


def add_theta_bins(df: pd.DataFrame, edges_by_projection: Dict[str, np.ndarray], theta_binning_period: str) -> pd.DataFrame:
    out = df.copy()

    for projection, edges in edges_by_projection.items():
        info = PROJECTIONS[projection]
        prefix = info["theta_source_prefix"]
        src = theta_col(prefix, theta_binning_period)

        vals = out[src].apply(parse_scalar).astype(float).to_numpy()
        finite = np.isfinite(vals)

        idx = np.searchsorted(edges, vals, side="right") - 1
        idx = np.clip(idx, 0, len(edges) - 2)

        min_col = info["min_col"]
        max_col = info["max_col"]

        out[min_col] = math.nan
        out[max_col] = math.nan
        out.loc[finite, min_col] = edges[idx[finite]]
        out.loc[finite, max_col] = edges[idx[finite] + 1]

    return out


def prepare_dataframes(dfs: Dict[str, pd.DataFrame], theta_binning_period: str, theta_bins: int, phi_degrees: bool) -> Dict[str, pd.DataFrame]:
    prepared = {label: add_width_columns(df, phi_degrees) for label, df in dfs.items()}
    edges = build_common_theta_edges(prepared, theta_binning_period, theta_bins)
    prepared = {label: add_theta_bins(df, edges, theta_binning_period) for label, df in prepared.items()}
    return prepared


def row_weight(row: pd.Series, widths: List[str], no_width_weighting: bool) -> float:
    if no_width_weighting:
        return 1.0
    w = 1.0
    for axis in widths:
        val = row[f"_width_{axis}"]
        if not np.isfinite(val) or val <= 0:
            return math.nan
        w *= float(val)
    return w


def x_position(group: pd.DataFrame, info: Dict[str, object], period: str) -> float:
    min_col = info["min_col"]
    max_col = info["max_col"]
    midpoint = 0.5 * (
        float(pd.to_numeric(group[min_col], errors="coerce").iloc[0])
        + float(pd.to_numeric(group[max_col], errors="coerce").iloc[0])
    )

    avg_col = f"{info['avg_prefix']}, {period}"
    if avg_col not in group.columns:
        return midpoint

    vals = group[avg_col].apply(parse_scalar).astype(float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return midpoint
    return float(np.mean(vals))


def integrated_points(df: pd.DataFrame, projection: str, period: str, template: str, no_width_weighting: bool) -> List[Point]:
    info = PROJECTIONS[projection]
    col = xs_col(period, template)

    points = []
    group_cols = [info["min_col"], info["max_col"]]

    for _, group in df.sort_values(group_cols).groupby(group_cols, sort=True, dropna=True):
        y = 0.0
        stat2 = 0.0
        sys2 = 0.0
        n = 0

        for _, row in group.iterrows():
            val, stat, sys = parse_tuple3(row[col])
            if not np.isfinite(val):
                continue

            w = row_weight(row, info["integrate_widths"], no_width_weighting)
            if not np.isfinite(w):
                continue

            stat = 0.0 if not np.isfinite(stat) else stat
            sys = 0.0 if not np.isfinite(sys) else sys

            y += val * w
            stat2 += (stat * w) ** 2
            sys2 += (sys * w) ** 2
            n += 1

        if n == 0:
            continue

        key = (float(group[info["min_col"]].iloc[0]), float(group[info["max_col"]].iloc[0]))
        points.append(Point(key=key, x=x_position(group, info, period), y=y, stat=math.sqrt(stat2), sys=math.sqrt(sys2)))

    points.sort(key=lambda p: p.x)
    return points


def total_error(point: Point, include_bin_to_bin_sys: bool, frac: float) -> float:
    stat = point.stat if np.isfinite(point.stat) else 0.0
    if not include_bin_to_bin_sys:
        return stat
    return math.sqrt(stat * stat + (frac * abs(point.y)) ** 2)


def plot_points(ax, points: List[Point], label: str, include_bin_to_bin_sys: bool, frac: float) -> None:
    if not points:
        return

    x = np.array([p.x for p in points])
    y = np.array([p.y for p in points])
    stat = np.array([p.stat if np.isfinite(p.stat) else 0.0 for p in points])
    outer = np.array([total_error(p, include_bin_to_bin_sys, frac) for p in points])

    style = CATEGORY_STYLE.get(label, dict(marker="o", linestyle="-"))
    color = style.get("color", None)

    if include_bin_to_bin_sys:
        ax.errorbar(x, y, yerr=outer, linewidth=0, elinewidth=2.0, capsize=4.0, alpha=0.25, color=color, zorder=1)

    ax.errorbar(x, y, yerr=stat, label=label, markersize=5.0, linewidth=1.1, elinewidth=1.0, capsize=2.0, zorder=3, **style)


def make_projection_canvas(
    dfs: Dict[str, pd.DataFrame],
    projection: str,
    output_dir: str,
    output_prefix: str,
    template: str,
    no_width_weighting: bool,
    include_bin_to_bin_sys: bool,
    frac: float,
) -> None:
    info = PROJECTIONS[projection]

    fig, axes = plt.subplots(2, 3, figsize=(18.0, 10.0), constrained_layout=False)
    axes_flat = axes.ravel()
    fig.suptitle(f"{info['title']} by topology", fontsize=16)

    for i, period in enumerate(RUN_PERIODS):
        ax = axes_flat[i]

        for label, df in dfs.items():
            pts = integrated_points(df, projection, period, template, no_width_weighting)
            plot_points(ax, pts, label, include_bin_to_bin_sys, frac)

        ax.set_title(period)
        ax.set_xlabel(info["xlabel"])
        ax.set_ylabel(info["ylabel"])
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8, frameon=True)

    axes_flat[5].axis("off")

    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.075, top=0.91, wspace=0.28, hspace=0.33)

    os.makedirs(output_dir, exist_ok=True)
    out = os.path.join(output_dir, f"{output_prefix}_{info['outfile_tag']}_comparison.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log(f"Wrote: {out}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare topology-resolved DVCS cross sections from separate CSVs.")
    p.add_argument("csv_files", nargs=3, help="CSV files in order: FD-FD CD-FT CD-FD")
    p.add_argument("--labels", nargs=3, default=DEFAULT_LABELS, help="Labels for the three CSVs.")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--output-prefix", default=DEFAULT_OUTPUT_PREFIX)
    p.add_argument("--xs-template", default=DEFAULT_XS_TEMPLATE, help="Cross-section column template with {period}.")
    p.add_argument("--no-width-weighting", action="store_true")
    p.add_argument("--phi-degrees", action="store_true")
    p.add_argument("--theta-bins", type=int, default=DEFAULT_THETA_BINS)
    p.add_argument("--theta-binning-period", default=DEFAULT_THETA_BINNING_PERIOD)
    p.add_argument("--include-bin-to-bin-sys", action="store_true")
    p.add_argument("--bin-to-bin-sys-frac", type=float, default=0.10)
    return p.parse_args()


def main() -> None:
    t0 = time.perf_counter()
    args = parse_args()

    log("============================================================")
    log("compare_cross_sections_by_topology.py")
    log("============================================================")
    log(f"CSV files: {args.csv_files}")
    log(f"Labels: {args.labels}")
    log(f"Output dir: {args.output_dir}")
    log(f"XS template: {args.xs_template}")

    dfs = read_inputs(args.csv_files, args.labels, args.xs_template)
    dfs = prepare_dataframes(dfs, args.theta_binning_period, args.theta_bins, args.phi_degrees)

    for projection in PROJECTION_ORDER:
        with Timer(f"plotting {projection}"):
            make_projection_canvas(
                dfs=dfs,
                projection=projection,
                output_dir=args.output_dir,
                output_prefix=args.output_prefix,
                template=args.xs_template,
                no_width_weighting=args.no_width_weighting,
                include_bin_to_bin_sys=args.include_bin_to_bin_sys,
                frac=args.bin_to_bin_sys_frac,
            )

    log(f"TOTAL RUNTIME: {time.perf_counter() - t0:.3f} s")
    log("Done.")


if __name__ == "__main__":
    main()
# endif