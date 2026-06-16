#!/usr/bin/env python3
"""
compare_cross_section_by_topology.py

Compare integrated DVCS cross sections from separate topology-resolved CSV files.

Expected usage:

  python compare_cross_section_by_topology.py FD-FD.csv CD-FT.csv CD-FD.csv

The input CSVs should each have the usual pass-2 columns, for example:

  normed cross sections, ep->epg, exp, Fa18 Inb, unpol
  normed cross sections, ep->epg, exp, Fa18 Out, unpol
  normed cross sections, ep->epg, exp, Sp19 Inb, unpol
  normed cross sections, ep->epg, exp, Sp18 Inb, unpol
  normed cross sections, ep->epg, exp, Sp18 Out, unpol

The script makes one 2x3 canvas per projection:

  xB, Q2, |t|, theta_e, theta_p, theta_gamma

The phi projection is intentionally omitted.

Each canvas has one subplot per run period. Each subplot has two panels:

  top:    absolute integrated cross sections
  bottom: each curve divided by the arithmetic average of the available curves
          in that same run-period panel and projected bin.

Run-period panel order:

  top row:    Sp18 Inb, Sp18 Out, empty
  bottom row: Fa18 Inb, Fa18 Out, Sp19 Inb

The ratio panel uses:

  R_i = sigma_i / <sigma>

where <sigma> is the arithmetic mean over the available CSV curves in that
period/projection bin.

A diagnostic chi2/ndf is printed in each run-period panel:

  chi2 = sum_i ((R_i - 1) / dR_i)^2
  ndf  = N_ratio_points - N_projected_bins

The uncertainty dR_i includes the statistical uncertainty on the numerator and
the propagated statistical uncertainty on the arithmetic mean. Correlations are
ignored, so this chi2 is a diagnostic consistency metric rather than a formal
hypothesis test.

For topology studies only, the script also makes an additional set of 2x3
canvases:

  output/topology_comparison/topology_period_ratio_xB_comparison.png
  output/topology_comparison/topology_period_ratio_Q2_comparison.png
  output/topology_comparison/topology_period_ratio_t_comparison.png
  output/topology_comparison/topology_period_ratio_e_theta_comparison.png
  output/topology_comparison/topology_period_ratio_p_theta_comparison.png
  output/topology_comparison/topology_period_ratio_g_theta_comparison.png

These canvases are arranged as:

  columns: FD-FD, CD-FT, CD-FD
  rows:    inbending, outbending

Each subplot shows the ratio topology / mean(topologies) for the appropriate run
periods. This is the fairer diagnostic for asking whether topology-fraction
patterns are stable between run periods with the same torus setting.

Outputs:

  output/topology_comparison/topology_xB_comparison.png
  output/topology_comparison/topology_Q2_comparison.png
  output/topology_comparison/topology_t_comparison.png
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

# Must happen before importing pyplot.
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg", force=True)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Defaults.
# -----------------------------------------------------------------------------

RUN_PERIODS = [
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out",
]

RUN_PERIOD_PANEL_ORDER = [
    "Sp18 Inb",
    "Sp18 Out",
    None,
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
]

INBENDING_PERIODS = [
    "Sp18 Inb",
    "Fa18 Inb",
    "Sp19 Inb",
]

OUTBENDING_PERIODS = [
    "Sp18 Out",
    "Fa18 Out",
]

DEFAULT_LABELS = [
    "FD-FD",
    "CD-FT",
    "CD-FD",
]

DEFAULT_OUTPUT_DIR = "output/topology_comparison"
DEFAULT_OUTPUT_PREFIX = "topology"
DEFAULT_THETA_BINS = 7
DEFAULT_THETA_BINNING_PERIOD = "10.6 GeV"
DEFAULT_XS_TEMPLATE = "normed cross sections, ep->epg, exp, {period}, unpol"

# 50/50 split between cross-section and ratio panels.
TOP_PANEL_HEIGHT = 1.0
RATIO_PANEL_HEIGHT = 1.0

# Topology stability canvases are deliberately standardized.
TOPOLOGY_STABILITY_RATIO_YMIN = 0.0
TOPOLOGY_STABILITY_RATIO_YMAX = 2.0

CATEGORY_STYLE = {
    "FD-FD": dict(marker="o", linestyle="-", color="tab:blue"),
    "CD-FT": dict(marker="s", linestyle="-", color="tab:orange"),
    "CD-FD": dict(marker="^", linestyle="-", color="tab:green"),

    "S1": dict(marker="o", linestyle="-", color="tab:blue"),
    "S2": dict(marker="s", linestyle="-", color="tab:orange"),
    "S3": dict(marker="^", linestyle="-", color="tab:green"),
    "S4": dict(marker="v", linestyle="-", color="tab:red"),
    "S5": dict(marker="D", linestyle="-", color="tab:purple"),
    "S6": dict(marker="P", linestyle="-", color="tab:brown"),

    "CD S1": dict(marker="o", linestyle="-", color="tab:blue"),
    "CD S2": dict(marker="s", linestyle="-", color="tab:orange"),
    "CD S3": dict(marker="^", linestyle="-", color="tab:green"),
}

PERIOD_STYLE = {
    "Sp18 Inb": dict(marker="D", linestyle="-", color="tab:green"),
    "Fa18 Inb": dict(marker="^", linestyle="-", color="tab:blue"),
    "Sp19 Inb": dict(marker="s", linestyle="-", color="tab:red"),
    "Sp18 Out": dict(marker="P", linestyle="-", color="tab:purple"),
    "Fa18 Out": dict(marker="v", linestyle="-", color="tab:orange"),
}

PROJECTIONS = {
    "xB": {
        "min_col": "xBmin",
        "max_col": "xBmax",
        "avg_prefix": "xBavg",
        "xlabel": r"$x_B$",
        "ylabel": r"$d\sigma/dx_B$  (pb)",
        "title": r"$x_B$ dependence",
        "integrate_widths": ["Q2", "t", "phi"],
        "outfile_tag": "xB",
        "derived_theta": False,
    },
    "Q2": {
        "min_col": "Q2min",
        "max_col": "Q2max",
        "avg_prefix": "Q2avg",
        "xlabel": r"$Q^2$  (GeV$^2$)",
        "ylabel": r"$d\sigma/dQ^2$  (pb/GeV$^2$)",
        "title": r"$Q^2$ dependence",
        "integrate_widths": ["xB", "t", "phi"],
        "outfile_tag": "Q2",
        "derived_theta": False,
    },
    "t": {
        "min_col": "t_abs_min",
        "max_col": "t_abs_max",
        "avg_prefix": "t_abs_avg",
        "xlabel": r"$|t|$  (GeV$^2$)",
        "ylabel": r"$d\sigma/d|t|$  (pb/GeV$^2$)",
        "title": r"$|t|$ dependence",
        "integrate_widths": ["xB", "Q2", "phi"],
        "outfile_tag": "t",
        "derived_theta": False,
    },
    "e_theta": {
        "min_col": "_bin_e_theta_min",
        "max_col": "_bin_e_theta_max",
        "avg_prefix": "e_theta",
        "theta_source_prefix": "e_theta",
        "xlabel": r"$\theta_{e}$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"electron polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "e_theta",
        "derived_theta": True,
    },
    "p_theta": {
        "min_col": "_bin_p_theta_min",
        "max_col": "_bin_p_theta_max",
        "avg_prefix": "p_theta",
        "theta_source_prefix": "p_theta",
        "xlabel": r"$\theta_{p}$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"proton polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "p_theta",
        "derived_theta": True,
    },
    "g_theta": {
        "min_col": "_bin_g_theta_min",
        "max_col": "_bin_g_theta_max",
        "avg_prefix": "g_theta",
        "theta_source_prefix": "g_theta",
        "xlabel": r"$\theta_{\gamma}$  (deg)",
        "ylabel": r"$\sigma_{\mathrm{int}}$  (pb)",
        "title": r"photon polar-angle dependence",
        "integrate_widths": ["xB", "Q2", "t", "phi"],
        "outfile_tag": "g_theta",
        "derived_theta": True,
    },
}

PROJECTION_ORDER = [
    "xB",
    "Q2",
    "t",
    "e_theta",
    "p_theta",
    "g_theta",
]

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


# -----------------------------------------------------------------------------
# Logging.
# -----------------------------------------------------------------------------

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
        dt = time.perf_counter() - self.t0
        log(f"DONE:  {self.label}  ({dt:.3f} s)")
        return False


# -----------------------------------------------------------------------------
# Parsing.
# -----------------------------------------------------------------------------

def parse_tuple3(value) -> Tuple[float, float, float]:
    if value is None:
        return math.nan, math.nan, math.nan
    # endif

    if isinstance(value, float) and math.isnan(value):
        return math.nan, math.nan, math.nan
    # endif

    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value), 0.0, 0.0
    # endif

    text = str(value).strip()

    if text == "" or text.lower() in {"nan", "none", "null"}:
        return math.nan, math.nan, math.nan
    # endif

    nums = FLOAT_PATTERN.findall(text)

    if not nums:
        return math.nan, math.nan, math.nan
    # endif

    vals = [float(x) for x in nums]

    while len(vals) < 3:
        vals.append(0.0)
    # endwhile

    return vals[0], vals[1], vals[2]


def parse_scalar(value) -> float:
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

    nums = FLOAT_PATTERN.findall(text)

    if not nums:
        return math.nan
    # endif

    return float(nums[0])


# -----------------------------------------------------------------------------
# Column helpers.
# -----------------------------------------------------------------------------

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
        # endfor

        relevant = [
            c for c in df.columns
            if "cross sections" in c.lower() or "cross section" in c.lower()
        ]

        if relevant:
            log("Cross-section-like columns present:")

            for c in relevant[:120]:
                print(f"  {c}")
            # endfor
        # endif

        raise KeyError(f"{label} is missing {len(missing)} required columns")
    # endif


# -----------------------------------------------------------------------------
# Input handling.
# -----------------------------------------------------------------------------

def required_input_columns(template: str) -> List[str]:
    cols = [
        "xBmin",
        "xBmax",
        "Q2min",
        "Q2max",
        "t_abs_min",
        "t_abs_max",
        "phimin",
        "phimax",
    ]

    for period in RUN_PERIODS:
        cols += [
            f"xBavg, {period}",
            f"Q2avg, {period}",
            f"t_abs_avg, {period}",
            f"phiavg, {period}",
            f"e_theta, {period}",
            f"p_theta, {period}",
            f"g_theta, {period}",
            xs_col(period, template),
        ]
    # endfor

    return sorted(set(cols))


def read_inputs(
    paths: List[str],
    labels: List[str],
    template: str,
) -> Dict[str, pd.DataFrame]:
    needed = required_input_columns(template)

    out: Dict[str, pd.DataFrame] = {}

    for path, label in zip(paths, labels):
        with Timer(f"reading {label}: {path}"):
            df = pd.read_csv(path, low_memory=False)
        # endwith

        require_columns(df, needed, label)

        out[label] = df

        log(f"{label}: {df.shape[0]} rows x {df.shape[1]} columns")
    # endfor

    return out


# -----------------------------------------------------------------------------
# Binning and preparation.
# -----------------------------------------------------------------------------

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


def build_common_theta_edges(
    dfs: Dict[str, pd.DataFrame],
    theta_binning_period: str,
    theta_bins: int,
) -> Dict[str, np.ndarray]:
    edges_by_projection: Dict[str, np.ndarray] = {}

    for projection in ["e_theta", "p_theta", "g_theta"]:
        prefix = str(PROJECTIONS[projection]["theta_source_prefix"])
        col = theta_col(prefix, theta_binning_period)

        all_vals = []

        for label, df in dfs.items():
            if col not in df.columns:
                raise KeyError(f"{label} missing theta-binning column {col}")
            # endif

            vals = df[col].apply(parse_scalar).astype(float).to_numpy()
            all_vals.append(vals[np.isfinite(vals)])
        # endfor

        merged = np.concatenate(all_vals) if all_vals else np.array([], dtype=float)

        if len(merged) == 0:
            raise RuntimeError(f"No finite theta values for {col}")
        # endif

        vmin = float(np.min(merged))
        vmax = float(np.max(merged))

        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            raise RuntimeError(f"Invalid theta range for {col}: {vmin}, {vmax}")
        # endif

        edges = np.linspace(vmin, vmax, theta_bins + 1)
        edges_by_projection[projection] = edges

        log(
            f"{projection}: common theta edges from {col}: "
            f"{theta_bins} bins, {vmin:.6g} to {vmax:.6g} deg"
        )
    # endfor

    return edges_by_projection


def add_theta_bins(
    df: pd.DataFrame,
    edges_by_projection: Dict[str, np.ndarray],
    theta_binning_period: str,
) -> pd.DataFrame:
    out = df.copy()

    for projection, edges in edges_by_projection.items():
        info = PROJECTIONS[projection]
        prefix = str(info["theta_source_prefix"])
        src = theta_col(prefix, theta_binning_period)

        vals = out[src].apply(parse_scalar).astype(float).to_numpy()
        finite = np.isfinite(vals)

        idx = np.searchsorted(edges, vals, side="right") - 1
        idx = np.clip(idx, 0, len(edges) - 2)

        min_col = str(info["min_col"])
        max_col = str(info["max_col"])

        out[min_col] = math.nan
        out[max_col] = math.nan

        out.loc[finite, min_col] = edges[idx[finite]]
        out.loc[finite, max_col] = edges[idx[finite] + 1]
    # endfor

    return out


def prepare_dataframes(
    dfs: Dict[str, pd.DataFrame],
    theta_binning_period: str,
    theta_bins: int,
    phi_degrees: bool,
) -> Dict[str, pd.DataFrame]:
    prepared = {
        label: add_width_columns(df, phi_degrees)
        for label, df in dfs.items()
    }

    edges = build_common_theta_edges(
        dfs=prepared,
        theta_binning_period=theta_binning_period,
        theta_bins=theta_bins,
    )

    prepared = {
        label: add_theta_bins(df, edges, theta_binning_period)
        for label, df in prepared.items()
    }

    return prepared


# -----------------------------------------------------------------------------
# Integration.
# -----------------------------------------------------------------------------

def row_weight(
    row: pd.Series,
    widths: List[str],
    no_width_weighting: bool,
) -> float:
    if no_width_weighting:
        return 1.0
    # endif

    w = 1.0

    for axis in widths:
        val = row[f"_width_{axis}"]

        if not np.isfinite(val) or val <= 0.0:
            return math.nan
        # endif

        w *= float(val)
    # endfor

    return w


def x_position(
    group: pd.DataFrame,
    info: Dict[str, object],
    period: str,
) -> float:
    min_col = str(info["min_col"])
    max_col = str(info["max_col"])

    midpoint = 0.5 * (
        float(pd.to_numeric(group[min_col], errors="coerce").iloc[0])
        + float(pd.to_numeric(group[max_col], errors="coerce").iloc[0])
    )

    avg_col = f"{info['avg_prefix']}, {period}"

    if avg_col not in group.columns:
        return midpoint
    # endif

    vals = group[avg_col].apply(parse_scalar).astype(float)
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return midpoint
    # endif

    return float(np.mean(vals))


def integrated_points(
    df: pd.DataFrame,
    projection: str,
    period: str,
    template: str,
    no_width_weighting: bool,
) -> List[Point]:
    info = PROJECTIONS[projection]
    col = xs_col(period, template)

    if col not in df.columns:
        return []
    # endif

    points: List[Point] = []
    group_cols = [str(info["min_col"]), str(info["max_col"])]

    sorted_df = df.sort_values(group_cols)

    for _, group in sorted_df.groupby(group_cols, sort=True, dropna=True):
        y = 0.0
        stat2 = 0.0
        sys2 = 0.0
        n = 0

        for _, row in group.iterrows():
            val, stat, sys = parse_tuple3(row[col])

            if not np.isfinite(val):
                continue
            # endif

            w = row_weight(
                row=row,
                widths=list(info["integrate_widths"]),
                no_width_weighting=no_width_weighting,
            )

            if not np.isfinite(w):
                continue
            # endif

            stat = 0.0 if not np.isfinite(stat) else stat
            sys = 0.0 if not np.isfinite(sys) else sys

            y += val * w
            stat2 += (stat * w) ** 2
            sys2 += (sys * w) ** 2
            n += 1
        # endfor

        if n == 0:
            continue
        # endif

        key = (
            float(group[str(info["min_col"])].iloc[0]),
            float(group[str(info["max_col"])].iloc[0]),
        )

        points.append(
            Point(
                key=key,
                x=x_position(group, info, period),
                y=y,
                stat=math.sqrt(stat2),
                sys=math.sqrt(sys2),
            )
        )
    # endfor

    points.sort(key=lambda p: p.x)

    return points


def precompute_all_points(
    dfs: Dict[str, pd.DataFrame],
    projection: str,
    template: str,
    no_width_weighting: bool,
) -> Dict[str, Dict[str, List[Point]]]:
    points_by_label_period: Dict[str, Dict[str, List[Point]]] = {}

    for label, df in dfs.items():
        points_by_label_period[label] = {}

        for period in RUN_PERIODS:
            points_by_label_period[label][period] = integrated_points(
                df=df,
                projection=projection,
                period=period,
                template=template,
                no_width_weighting=no_width_weighting,
            )
        # endfor
    # endfor

    return points_by_label_period


# -----------------------------------------------------------------------------
# Ratio to average and chi2.
# -----------------------------------------------------------------------------

def average_points_by_key(
    points_by_label: Dict[str, List[Point]],
) -> Dict[Tuple[float, float], Point]:
    all_keys = sorted(
        {
            p.key
            for points in points_by_label.values()
            for p in points
        }
    )

    lookup = {
        label: {p.key: p for p in points}
        for label, points in points_by_label.items()
    }

    avg_by_key: Dict[Tuple[float, float], Point] = {}

    for key in all_keys:
        used: List[Point] = []

        for label in points_by_label:
            p = lookup[label].get(key)

            if p is None:
                continue
            # endif

            if not np.isfinite(p.y):
                continue
            # endif

            used.append(p)
        # endfor

        if len(used) == 0:
            continue
        # endif

        yvals = np.array([p.y for p in used], dtype=float)
        xvals = np.array([p.x for p in used], dtype=float)
        statvals = np.array(
            [
                p.stat if np.isfinite(p.stat) else 0.0
                for p in used
            ],
            dtype=float,
        )

        mean_y = float(np.mean(yvals))
        mean_x = float(np.mean(xvals))
        mean_stat = float(math.sqrt(np.sum(statvals * statvals)) / len(used))

        avg_by_key[key] = Point(
            key=key,
            x=mean_x,
            y=mean_y,
            stat=mean_stat,
            sys=0.0,
        )
    # endfor

    return avg_by_key


def ratio_points_to_average(
    points: List[Point],
    avg_by_key: Dict[Tuple[float, float], Point],
) -> List[Point]:
    ratios: List[Point] = []

    for p in points:
        avg = avg_by_key.get(p.key)

        if avg is None:
            continue
        # endif

        if not np.isfinite(p.y) or not np.isfinite(avg.y) or avg.y == 0.0:
            continue
        # endif

        r = p.y / avg.y
        rel2 = 0.0

        if p.y != 0.0 and np.isfinite(p.stat):
            rel2 += (p.stat / p.y) ** 2
        # endif

        if avg.y != 0.0 and np.isfinite(avg.stat):
            rel2 += (avg.stat / avg.y) ** 2
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

    ratios.sort(key=lambda p: p.x)

    return ratios


def ratios_for_period(
    points_by_label_period: Dict[str, Dict[str, List[Point]]],
    period: str,
) -> Dict[str, List[Point]]:
    points_by_label = {
        label: points_by_period.get(period, [])
        for label, points_by_period in points_by_label_period.items()
    }

    avg_by_key = average_points_by_key(points_by_label)

    return {
        label: ratio_points_to_average(points, avg_by_key)
        for label, points in points_by_label.items()
    }


def chi2_ndf_for_ratio_points(
    ratio_points_by_label: Dict[str, List[Point]],
) -> Tuple[float, int, int]:
    chi2 = 0.0
    n_used = 0
    used_keys = set()

    for points in ratio_points_by_label.values():
        for p in points:
            if not np.isfinite(p.y):
                continue
            # endif

            if not np.isfinite(p.stat) or p.stat <= 0.0:
                continue
            # endif

            chi2 += ((p.y - 1.0) / p.stat) ** 2
            n_used += 1
            used_keys.add(p.key)
        # endfor
    # endfor

    n_constraints = len(used_keys)
    ndf = n_used - n_constraints

    return chi2, ndf, n_used


def chi2_label_text(
    ratio_points_by_label: Dict[str, List[Point]],
) -> str:
    chi2, ndf, n_used = chi2_ndf_for_ratio_points(ratio_points_by_label)

    if n_used <= 0:
        return r"$\chi^2/\mathrm{ndf}$: n/a"
    # endif

    if ndf <= 0:
        return rf"$\chi^2$: {chi2:.1f}"
    # endif

    return rf"$\chi^2/\mathrm{{ndf}} = {chi2:.1f}/{ndf:d} = {chi2 / ndf:.2f}$"


# -----------------------------------------------------------------------------
# Plot helpers.
# -----------------------------------------------------------------------------

def total_error(
    point: Point,
    include_bin_to_bin_sys: bool,
    frac: float,
) -> float:
    stat = point.stat if np.isfinite(point.stat) else 0.0

    if not include_bin_to_bin_sys:
        return stat
    # endif

    return math.sqrt(stat * stat + (frac * abs(point.y)) ** 2)


def collect_y_values(points_by_label: Dict[str, List[Point]]) -> List[float]:
    values: List[float] = []

    for points in points_by_label.values():
        for p in points:
            if np.isfinite(p.y):
                values.append(float(p.y))
            # endif
        # endfor
    # endfor

    return values


def collect_all_absolute_y_values(
    points_by_label_period: Dict[str, Dict[str, List[Point]]],
) -> List[float]:
    values: List[float] = []

    for points_by_period in points_by_label_period.values():
        for period in RUN_PERIODS:
            for p in points_by_period.get(period, []):
                if np.isfinite(p.y):
                    values.append(float(p.y))
                # endif
            # endfor
        # endfor
    # endfor

    return values


def collect_all_ratio_y_values(
    ratio_points_by_period_label: Dict[str, Dict[str, List[Point]]],
) -> List[float]:
    values: List[float] = []

    for ratio_points_by_label in ratio_points_by_period_label.values():
        values.extend(collect_y_values(ratio_points_by_label))
    # endfor

    return values


def dynamic_absolute_ylim_from_values(values: List[float]) -> Tuple[float, float]:
    if len(values) == 0:
        return 0.0, 0.01
    # endif

    ymax_data = max(values)

    if not np.isfinite(ymax_data) or ymax_data <= 0.0:
        return 0.0, 0.01
    # endif

    ymax = math.ceil(ymax_data / 0.01) * 0.01

    if ymax <= 0.0:
        ymax = 0.01
    # endif

    return 0.0, ymax


def dynamic_ratio_ylim_from_values(values: List[float]) -> Tuple[float, float]:
    if len(values) == 0:
        return 0.8, 1.2
    # endif

    ymin_data = min(values)
    ymax_data = max(values)

    if not np.isfinite(ymin_data) or not np.isfinite(ymax_data):
        return 0.8, 1.2
    # endif

    ymin = math.floor(ymin_data * 10.0) / 10.0
    ymax = math.ceil(ymax_data * 10.0) / 10.0

    if ymin == ymax:
        ymin -= 0.1
        ymax += 0.1
    # endif

    if ymin > 1.0:
        ymin = 1.0
    # endif

    if ymax < 1.0:
        ymax = 1.0
    # endif

    return ymin, ymax


def plot_points(
    ax,
    points: List[Point],
    label: str,
    include_bin_to_bin_sys: bool,
    frac: float,
) -> None:
    if not points:
        return
    # endif

    x = np.array([p.x for p in points], dtype=float)
    y = np.array([p.y for p in points], dtype=float)

    stat = np.array(
        [
            p.stat if np.isfinite(p.stat) else 0.0
            for p in points
        ],
        dtype=float,
    )

    outer = np.array(
        [
            total_error(p, include_bin_to_bin_sys, frac)
            for p in points
        ],
        dtype=float,
    )

    style = CATEGORY_STYLE.get(label, PERIOD_STYLE.get(label, dict(marker="o", linestyle="-")))
    color = style.get("color", None)

    if include_bin_to_bin_sys:
        ax.errorbar(
            x,
            y,
            yerr=outer,
            linewidth=0.0,
            elinewidth=2.0,
            capsize=4.0,
            alpha=0.25,
            color=color,
            zorder=1,
        )
    # endif

    ax.errorbar(
        x,
        y,
        yerr=stat,
        label=label,
        markersize=5.0,
        linewidth=1.1,
        elinewidth=1.0,
        capsize=2.0,
        zorder=3,
        **style,
    )


def plot_ratio_points(
    ax,
    points: List[Point],
    label: str,
    use_period_style: bool = False,
) -> None:
    if not points:
        return
    # endif

    x = np.array([p.x for p in points], dtype=float)
    y = np.array([p.y for p in points], dtype=float)

    stat = np.array(
        [
            p.stat if np.isfinite(p.stat) else 0.0
            for p in points
        ],
        dtype=float,
    )

    if use_period_style:
        style = PERIOD_STYLE.get(label, dict(marker="o", linestyle="-"))
    else:
        style = CATEGORY_STYLE.get(label, dict(marker="o", linestyle="-"))
    # endif

    ax.errorbar(
        x,
        y,
        yerr=stat,
        label=label,
        markersize=4.0,
        linewidth=1.0,
        elinewidth=0.9,
        capsize=2.0,
        zorder=3,
        **style,
    )


def add_chi2_box(ax, text: str) -> None:
    ax.text(
        0.025,
        0.95,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="0.65",
            alpha=0.85,
        ),
    )


# -----------------------------------------------------------------------------
# Standard 2x3 run-period canvases.
# -----------------------------------------------------------------------------

def make_projection_canvas(
    dfs: Dict[str, pd.DataFrame],
    projection: str,
    output_dir: str,
    output_prefix: str,
    template: str,
    no_width_weighting: bool,
    include_bin_to_bin_sys: bool,
    frac: float,
    comparison_name: str,
) -> None:
    info = PROJECTIONS[projection]

    with Timer(f"integrating {projection}"):
        points_by_label_period = precompute_all_points(
            dfs=dfs,
            projection=projection,
            template=template,
            no_width_weighting=no_width_weighting,
        )
    # endwith

    ratio_points_by_period_label = {
        period: ratios_for_period(
            points_by_label_period=points_by_label_period,
            period=period,
        )
        for period in RUN_PERIODS
    }

    absolute_ylim = dynamic_absolute_ylim_from_values(
        collect_all_absolute_y_values(points_by_label_period)
    )

    ratio_ylim = dynamic_ratio_ylim_from_values(
        collect_all_ratio_y_values(ratio_points_by_period_label)
    )

    fig = plt.figure(figsize=(18.5, 10.4))

    outer = fig.add_gridspec(
        2,
        3,
        left=0.055,
        right=0.992,
        bottom=0.060,
        top=0.920,
        wspace=0.205,
        hspace=0.265,
    )

    fig.suptitle(f"{info['title']} by {comparison_name}", fontsize=16)

    for panel_index, period in enumerate(RUN_PERIOD_PANEL_ORDER):
        row = panel_index // 3
        col = panel_index % 3

        inner = outer[row, col].subgridspec(
            2,
            1,
            height_ratios=[TOP_PANEL_HEIGHT, RATIO_PANEL_HEIGHT],
            hspace=0.035,
        )

        ax_top = fig.add_subplot(inner[0])
        ax_ratio = fig.add_subplot(inner[1], sharex=ax_top)

        if period is None:
            ax_top.axis("off")
            ax_ratio.axis("off")
            continue
        # endif

        points_by_label = {
            label: points_by_period.get(period, [])
            for label, points_by_period in points_by_label_period.items()
        }

        for label, points in points_by_label.items():
            plot_points(
                ax=ax_top,
                points=points,
                label=label,
                include_bin_to_bin_sys=include_bin_to_bin_sys,
                frac=frac,
            )
        # endfor

        ratio_points_by_label = ratio_points_by_period_label[period]

        for label, ratios in ratio_points_by_label.items():
            plot_ratio_points(
                ax=ax_ratio,
                points=ratios,
                label=label,
                use_period_style=False,
            )
        # endfor

        add_chi2_box(
            ax=ax_top,
            text=chi2_label_text(ratio_points_by_label),
        )

        ax_top.set_title(period, fontsize=11)
        ax_top.set_ylabel(str(info["ylabel"]))
        ax_top.set_ylim(*absolute_ylim)
        ax_top.grid(True, alpha=0.25)
        ax_top.legend(fontsize=8, frameon=True, loc="best")

        ax_ratio.axhline(
            1.0,
            color="0.35",
            linewidth=1.0,
            linestyle="--",
            zorder=0,
        )

        ax_ratio.set_ylim(*ratio_ylim)
        ax_ratio.set_xlabel(str(info["xlabel"]))
        ax_ratio.set_ylabel("ratio")
        ax_ratio.grid(True, alpha=0.25)

        plt.setp(ax_top.get_xticklabels(), visible=False)
    # endfor

    os.makedirs(output_dir, exist_ok=True)

    out = os.path.join(
        output_dir,
        f"{output_prefix}_{info['outfile_tag']}_comparison.png",
    )

    fig.savefig(out, dpi=150)
    plt.close(fig)

    log(f"Wrote: {out}")


def make_all_projection_canvases(
    dfs: Dict[str, pd.DataFrame],
    output_dir: str,
    output_prefix: str,
    template: str,
    no_width_weighting: bool,
    include_bin_to_bin_sys: bool,
    frac: float,
    comparison_name: str,
) -> None:
    for projection in PROJECTION_ORDER:
        with Timer(f"plotting {projection}"):
            make_projection_canvas(
                dfs=dfs,
                projection=projection,
                output_dir=output_dir,
                output_prefix=output_prefix,
                template=template,
                no_width_weighting=no_width_weighting,
                include_bin_to_bin_sys=include_bin_to_bin_sys,
                frac=frac,
                comparison_name=comparison_name,
            )
        # endwith
    # endfor


# -----------------------------------------------------------------------------
# Topology-specific period-comparison ratio canvases.
# -----------------------------------------------------------------------------

def make_topology_period_ratio_canvas(
    dfs: Dict[str, pd.DataFrame],
    projection: str,
    output_dir: str,
    output_prefix: str,
    template: str,
    no_width_weighting: bool,
) -> None:
    info = PROJECTIONS[projection]

    labels = list(dfs.keys())

    with Timer(f"integrating topology-period-ratio {projection}"):
        points_by_label_period = precompute_all_points(
            dfs=dfs,
            projection=projection,
            template=template,
            no_width_weighting=no_width_weighting,
        )
    # endwith

    ratios_by_period = {
        period: ratios_for_period(
            points_by_label_period=points_by_label_period,
            period=period,
        )
        for period in RUN_PERIODS
    }

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(18.5, 8.2),
        constrained_layout=False,
    )

    fig.suptitle(
        f"{info['title']}: topology ratio stability by torus setting",
        fontsize=16,
    )

    row_periods = [
        ("inbending", INBENDING_PERIODS),
        ("outbending", OUTBENDING_PERIODS),
    ]

    for row_index, (row_label, periods) in enumerate(row_periods):
        for col_index, topology_label in enumerate(labels):
            ax = axes[row_index, col_index]

            for period in periods:
                ratios = ratios_by_period.get(period, {}).get(topology_label, [])

                plot_ratio_points(
                    ax=ax,
                    points=ratios,
                    label=period,
                    use_period_style=True,
                )
            # endfor

            ax.axhline(
                1.0,
                color="0.35",
                linewidth=1.0,
                linestyle="--",
                zorder=0,
            )

            ax.set_ylim(
                TOPOLOGY_STABILITY_RATIO_YMIN,
                TOPOLOGY_STABILITY_RATIO_YMAX,
            )

            ax.set_title(f"{topology_label} / mean, {row_label}", fontsize=11)
            ax.set_xlabel(str(info["xlabel"]))
            ax.set_ylabel("topology / mean")
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=8, frameon=True, loc="best")
        # endfor
    # endfor

    fig.subplots_adjust(
        left=0.055,
        right=0.992,
        bottom=0.070,
        top=0.900,
        wspace=0.205,
        hspace=0.300,
    )

    os.makedirs(output_dir, exist_ok=True)

    out = os.path.join(
        output_dir,
        f"{output_prefix}_period_ratio_{info['outfile_tag']}_comparison.png",
    )

    fig.savefig(out, dpi=150)
    plt.close(fig)

    log(f"Wrote: {out}")


def make_all_topology_period_ratio_canvases(
    dfs: Dict[str, pd.DataFrame],
    output_dir: str,
    output_prefix: str,
    template: str,
    no_width_weighting: bool,
) -> None:
    for projection in PROJECTION_ORDER:
        with Timer(f"plotting topology-period-ratio {projection}"):
            make_topology_period_ratio_canvas(
                dfs=dfs,
                projection=projection,
                output_dir=output_dir,
                output_prefix=output_prefix,
                template=template,
                no_width_weighting=no_width_weighting,
            )
        # endwith
    # endfor


# -----------------------------------------------------------------------------
# CLI.
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare topology-resolved DVCS cross sections from separate CSVs."
    )

    parser.add_argument(
        "csv_files",
        nargs=3,
        help="CSV files in order: FD-FD CD-FT CD-FD",
    )

    parser.add_argument(
        "--labels",
        nargs=3,
        default=DEFAULT_LABELS,
        help="Labels for the three CSVs. Default: FD-FD CD-FT CD-FD",
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        help=f"Output filename prefix. Default: {DEFAULT_OUTPUT_PREFIX}",
    )

    parser.add_argument(
        "--xs-template",
        default=DEFAULT_XS_TEMPLATE,
        help=f"Cross-section column template with {{period}}. Default: {DEFAULT_XS_TEMPLATE!r}",
    )

    parser.add_argument(
        "--no-width-weighting",
        action="store_true",
        help="Use raw sums instead of bin-width weighted integrations.",
    )

    parser.add_argument(
        "--phi-degrees",
        action="store_true",
        help="Use phi bin widths in degrees instead of radians.",
    )

    parser.add_argument(
        "--theta-bins",
        type=int,
        default=DEFAULT_THETA_BINS,
        help=f"Number of derived theta bins. Default: {DEFAULT_THETA_BINS}",
    )

    parser.add_argument(
        "--theta-binning-period",
        default=DEFAULT_THETA_BINNING_PERIOD,
        help=f"Theta binning reference period. Default: {DEFAULT_THETA_BINNING_PERIOD!r}",
    )

    parser.add_argument(
        "--include-bin-to-bin-sys",
        action="store_true",
        help="Draw outer stat+fractional-bin-to-bin-systematic bars.",
    )

    parser.add_argument(
        "--bin-to-bin-sys-frac",
        type=float,
        default=0.10,
        help="Fractional bin-to-bin systematic uncertainty. Default: 0.10",
    )

    parser.add_argument(
        "--comparison-name",
        default="topology",
        help="Text used in figure titles. Default: topology",
    )

    parser.add_argument(
        "--no-topology-period-ratio-plots",
        action="store_true",
        help="Disable the additional topology ratio stability canvases.",
    )

    return parser.parse_args()


def main() -> None:
    t0 = time.perf_counter()
    args = parse_args()

    log("============================================================")
    log("compare_cross_section_by_topology.py")
    log("============================================================")
    log(f"CSV files: {args.csv_files}")
    log(f"Labels: {args.labels}")
    log(f"Output dir: {args.output_dir}")
    log(f"XS template: {args.xs_template}")
    log(f"Theta bins: {args.theta_bins}")
    log(f"Theta binning period: {args.theta_binning_period}")
    log("phi dependence plots: disabled")

    if "{period}" not in args.xs_template:
        raise ValueError("--xs-template must contain {period}")
    # endif

    if args.theta_bins <= 0:
        raise ValueError("--theta-bins must be positive")
    # endif

    if args.bin_to_bin_sys_frac < 0.0:
        raise ValueError("--bin-to-bin-sys-frac must be non-negative")
    # endif

    dfs = read_inputs(
        paths=args.csv_files,
        labels=args.labels,
        template=args.xs_template,
    )

    dfs = prepare_dataframes(
        dfs=dfs,
        theta_binning_period=args.theta_binning_period,
        theta_bins=args.theta_bins,
        phi_degrees=args.phi_degrees,
    )

    make_all_projection_canvases(
        dfs=dfs,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        template=args.xs_template,
        no_width_weighting=args.no_width_weighting,
        include_bin_to_bin_sys=args.include_bin_to_bin_sys,
        frac=args.bin_to_bin_sys_frac,
        comparison_name=args.comparison_name,
    )

    if not args.no_topology_period_ratio_plots:
        make_all_topology_period_ratio_canvases(
            dfs=dfs,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
            template=args.xs_template,
            no_width_weighting=args.no_width_weighting,
        )
    # endif

    log(f"TOTAL RUNTIME: {time.perf_counter() - t0:.3f} s")
    log("Done.")


if __name__ == "__main__":
    main()
# endif