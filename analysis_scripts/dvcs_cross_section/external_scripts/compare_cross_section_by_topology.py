#!/usr/bin/env python3
"""
compare_cross_section_by_topology.py

Compare sector/topology-resolved DVCS pass-2 CSV outputs.

Original functionality is preserved:

  python3 compare_cross_section_by_topology.py FD-FD.csv CD-FT.csv CD-FD.csv

still makes the standard normed-cross-section comparison canvases:

  topology_xB_comparison.png
  topology_Q2_comparison.png
  topology_t_comparison.png
  topology_e_theta_comparison.png
  topology_p_theta_comparison.png
  topology_g_theta_comparison.png

and, unless disabled, topology-period stability canvases.

Expanded functionality:

The same plotting engine can now also make step-by-step diagnostic comparison
plots for quantities upstream of the final cross section, including raw yields,
current factors, normalized raw yields, generated/reconstructed MC yields, pi0
contamination quantities, signal yields, acceptance, and acceptance-corrected
yields. These are useful for finding where a sector/topology dependence first
enters the analysis chain.

For projection canvases, each subplot has:

  top:    integrated/averaged value for each input CSV label
  bottom: each input CSV label divided by the arithmetic mean of available labels

Ratio panels are standardized to 0 <= ratio <= 2 by default.
"""

import argparse
import math
import os
import re
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

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

DEFAULT_RATIO_YMIN = 0.0
DEFAULT_RATIO_YMAX = 2.0

# 50/50 split between absolute and ratio panels.
TOP_PANEL_HEIGHT = 1.0
RATIO_PANEL_HEIGHT = 1.0

# Topology stability canvases are deliberately standardized.
TOPOLOGY_STABILITY_RATIO_YMIN = 0.0
TOPOLOGY_STABILITY_RATIO_YMAX = 2.0

TOPOLOGY_LABEL_TO_CSV_TOPOLOGY = {
    "FD-FD": "(FD, FD)",
    "CD-FT": "(CD, FT)",
    "CD-FD": "(CD, FD)",
    "(FD, FD)": "(FD, FD)",
    "(CD, FT)": "(CD, FT)",
    "(CD, FD)": "(CD, FD)",
}

CSV_TOPOLOGIES = [
    "(FD, FD)",
    "(CD, FD)",
    "(CD, FT)",
]

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
        "ylabel_default": r"value",
        "title": r"$x_B$ dependence",
        "integrate_widths": ["Q2", "t", "phi"],
        "outfile_tag": "xB",
        "derived_theta": False,
    },
    "Q2": {
        "min_col": "Q2min",
        "max_col": "Q2max",
        "avg_prefix": "Q2avg",
        "xlabel": r"$Q^2$ (GeV$^2$)",
        "ylabel_default": r"value",
        "title": r"$Q^2$ dependence",
        "integrate_widths": ["xB", "t", "phi"],
        "outfile_tag": "Q2",
        "derived_theta": False,
    },
    "t": {
        "min_col": "t_abs_min",
        "max_col": "t_abs_max",
        "avg_prefix": "t_abs_avg",
        "xlabel": r"$|t|$ (GeV$^2$)",
        "ylabel_default": r"value",
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
        "xlabel": r"$\theta_{e}$ (deg)",
        "ylabel_default": r"value",
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
        "xlabel": r"$\theta_{p}$ (deg)",
        "ylabel_default": r"value",
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
        "xlabel": r"$\theta_{\gamma}$ (deg)",
        "ylabel_default": r"value",
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


@dataclass
class MetricSpec:
    tag: str
    title: str
    ylabel: str
    aggregate: str  # "sum", "mean", or "scalar"
    templates: List[str]
    width_weighted: bool = False
    prefer_label_topology: bool = False
    sum_all_topologies_if_no_label_match: bool = False
    required: bool = False


# -----------------------------------------------------------------------------
# Metric definitions.
# -----------------------------------------------------------------------------

def make_template_metric(
    tag: str,
    title: str,
    ylabel: str,
    template: str,
    aggregate: str,
    width_weighted: bool = False,
    required: bool = False,
) -> MetricSpec:
    return MetricSpec(
        tag=tag,
        title=title,
        ylabel=ylabel,
        aggregate=aggregate,
        templates=[template],
        width_weighted=width_weighted,
        required=required,
    )


def topology_template(prefix: str, channel: str, sample: str, suffix: str) -> str:
    return f"{prefix}, {channel}, {{topology}}, {sample}, {{period}}{suffix}"


def all_step_metric_specs() -> List[MetricSpec]:
    """
    The expanded pass-2 chain diagnostics.

    Notes on aggregation:
      - yields are summed over integrated dimensions;
      - cross sections use bin-width-weighted integration by default;
      - efficiency/acceptance/contamination-like quantities are averaged.
    """

    specs: List[MetricSpec] = []

    # Existing/default physics output.
    specs.append(
        make_template_metric(
            tag="normed_cross_section_epg_unpol",
            title="normed DVCS cross section",
            ylabel=r"integrated normed cross section",
            template="normed cross sections, ep->epg, exp, {period}, unpol",
            aggregate="sum",
            width_weighted=True,
            required=False,
        )
    )

    specs.append(
        make_template_metric(
            tag="cross_section_epg_unpol",
            title="raw-normalization DVCS cross section",
            ylabel=r"integrated cross section",
            template="cross sections, ep->epg, exp, {period}, unpol",
            aggregate="sum",
            width_weighted=True,
            required=False,
        )
    )

    # Raw and normalized DATA yields for DVCS and eppi0.
    for channel, short in [
        ("ep->epg", "epg"),
        ("ep->eppi0", "eppi0"),
    ]:
        specs.append(
            MetricSpec(
                tag=f"raw_yield_{short}_unpol",
                title=f"raw yield, {channel}, unpol",
                ylabel="raw yield",
                aggregate="sum",
                templates=[topology_template("raw yield", channel, "exp", ", unpol")],
                width_weighted=False,
                prefer_label_topology=True,
                sum_all_topologies_if_no_label_match=True,
            )
        )

        specs.append(
            MetricSpec(
                tag=f"normalized_raw_yield_{short}_unpol",
                title=f"normalized raw yield, {channel}, unpol",
                ylabel="normalized raw yield",
                aggregate="sum",
                templates=[topology_template("normalized raw yield", channel, "exp", ", unpol")],
                width_weighted=False,
                prefer_label_topology=True,
                sum_all_topologies_if_no_label_match=True,
            )
        )

    # Current efficiency factors are period-level scalars.
    for channel, short in [
        ("ep->epg", "epg"),
        ("ep->eppi0", "eppi0"),
    ]:
        for sample, sample_tag in [
            ("exp", "data"),
            ("mc", "mc"),
        ]:
            specs.append(
                MetricSpec(
                    tag=f"current_efficiency_factor_{short}_{sample_tag}",
                    title=f"current-efficiency factor, {channel}, {sample_tag.upper()}",
                    ylabel="current-efficiency factor",
                    aggregate="scalar",
                    templates=[f"current efficiency factor, {channel}, {sample}, {{period}}"],
                    width_weighted=False,
                )
            )

    # Generated MC yields.
    for channel, short in [
        ("ep->epg", "epg"),
        ("ep->eppi0", "eppi0"),
    ]:
        specs.append(
            make_template_metric(
                tag=f"generated_yield_{short}",
                title=f"generated yield, {channel}, MC",
                ylabel="generated yield",
                template=f"generated yield, {channel}, mc, {{period}}",
                aggregate="sum",
                width_weighted=False,
            )
        )

    # Reconstructed MC yields, including pi0 background reconstructed as epg.
    for channel, short in [
        ("ep->epg", "epg"),
        ("ep->eppi0", "eppi0"),
        ("ep->eppi0->epg", "eppi0_to_epg"),
    ]:
        specs.append(
            MetricSpec(
                tag=f"reconstructed_yield_{short}",
                title=f"reconstructed yield, {channel}, MC",
                ylabel="reconstructed yield",
                aggregate="sum",
                templates=[
                    topology_template("reconstructed yield", channel, "mc", ""),
                    f"reconstructed yield, {channel}, mc, {{period}}",
                ],
                width_weighted=False,
                prefer_label_topology=True,
                sum_all_topologies_if_no_label_match=False,
            )
        )

        specs.append(
            MetricSpec(
                tag=f"reconstructed_current_corrected_yield_{short}",
                title=f"reconstructed current-corrected yield, {channel}, MC",
                ylabel="reconstructed current-corrected yield",
                aggregate="sum",
                templates=[
                    topology_template("reconstructed current corrected yield", channel, "mc", ""),
                    f"reconstructed current corrected yield, {channel}, mc, {{period}}",
                ],
                width_weighted=False,
                prefer_label_topology=True,
                sum_all_topologies_if_no_label_match=False,
            )
        )

    # Pi0 contamination and final DATA yield chain.
    specs.append(
        make_template_metric(
            tag="contamination_ratio",
            title="pi0 contamination ratio",
            ylabel="contamination ratio",
            template="contamination ratio, {period}",
            aggregate="mean",
            width_weighted=False,
        )
    )

    specs.append(
        make_template_metric(
            tag="signal_yield_epg_unpol",
            title="signal yield, ep->epg, unpol",
            ylabel="signal yield",
            template="signal yield, ep->epg, exp, {period}, unpol",
            aggregate="sum",
            width_weighted=False,
        )
    )

    specs.append(
        make_template_metric(
            tag="acceptance",
            title="acceptance",
            ylabel="acceptance",
            template="acceptance, {period}",
            aggregate="mean",
            width_weighted=False,
        )
    )

    specs.append(
        make_template_metric(
            tag="acceptance_corrected_yield_epg_unpol",
            title="acceptance-corrected yield, ep->epg, unpol",
            ylabel="acceptance-corrected yield",
            template="acceptance corrected yield, ep->epg, exp, {period}, unpol",
            aggregate="sum",
            width_weighted=False,
        )
    )

    return specs


def metric_specs_by_tag() -> Dict[str, MetricSpec]:
    return {m.tag: m for m in all_step_metric_specs()}


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

        if exc_type is None:
            log(f"DONE:  {self.label}  ({dt:.3f} s)")
        else:
            log(f"FAILED: {self.label}  ({dt:.3f} s)  {exc_type.__name__}: {exc}")

        return False


# -----------------------------------------------------------------------------
# Parsing.
# -----------------------------------------------------------------------------

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
    val, _, _ = parse_tuple3(value)
    return val


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

        relevant = [
            c for c in df.columns
            if (
                "cross sections" in c.lower()
                or "yield" in c.lower()
                or "acceptance" in c.lower()
                or "contamination" in c.lower()
                or "current efficiency" in c.lower()
            )
        ]

        if relevant:
            log("Relevant columns present:")

            for c in relevant[:180]:
                print(f"  {c}")

        raise KeyError(f"{label} is missing {len(missing)} required columns")


def base_required_columns(template: str) -> List[str]:
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
        ]

    return sorted(set(cols))


def required_input_columns(template: str) -> List[str]:
    cols = base_required_columns(template)

    for period in RUN_PERIODS:
        cols.append(xs_col(period, template))

    return sorted(set(cols))


def label_to_csv_topology(label: str) -> Optional[str]:
    return TOPOLOGY_LABEL_TO_CSV_TOPOLOGY.get(str(label).strip())


def metric_columns_for_label_period(
    df: pd.DataFrame,
    label: str,
    period: str,
    metric: MetricSpec,
) -> List[str]:
    """
    Resolve one metric to one or more actual CSV columns.

    For topology-aware columns:
      - if label is a topology label, use the matching topology column;
      - if label is not a topology label and requested, sum all topology columns;
      - for reconstructed-yield metrics, fall back to the total MC column.
    """

    resolved: List[str] = []
    csv_topology = label_to_csv_topology(label)

    for template in metric.templates:
        if "{topology}" in template:
            if metric.prefer_label_topology and csv_topology is not None:
                col = template.format(period=period, topology=csv_topology)

                if col in df.columns:
                    resolved.append(col)
                    continue

            if metric.sum_all_topologies_if_no_label_match or csv_topology is None:
                topo_cols = [
                    template.format(period=period, topology=topo)
                    for topo in CSV_TOPOLOGIES
                    if template.format(period=period, topology=topo) in df.columns
                ]

                if topo_cols:
                    resolved.extend(topo_cols)
                    continue

            # If we get here, the topology template did not resolve.
            continue

        col = template.format(period=period)

        if col in df.columns:
            resolved.append(col)
            continue

    # Remove duplicates while preserving order.
    unique = []
    seen = set()

    for c in resolved:
        if c not in seen:
            unique.append(c)
            seen.add(c)

    return unique


# -----------------------------------------------------------------------------
# Input handling.
# -----------------------------------------------------------------------------

def read_inputs(
    paths: List[str],
    labels: List[str],
    template: str,
    require_template_columns: bool = True,
) -> Dict[str, pd.DataFrame]:
    needed = required_input_columns(template) if require_template_columns else base_required_columns(template)

    out: Dict[str, pd.DataFrame] = {}

    for path, label in zip(paths, labels):
        with Timer(f"reading {label}: {path}"):
            df = pd.read_csv(path, low_memory=False)

        require_columns(df, needed, label)

        out[label] = df

        log(f"{label}: {df.shape[0]} rows x {df.shape[1]} columns")

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

            vals = df[col].apply(parse_scalar).astype(float).to_numpy()
            all_vals.append(vals[np.isfinite(vals)])

        merged = np.concatenate(all_vals) if all_vals else np.array([], dtype=float)

        if len(merged) == 0:
            raise RuntimeError(f"No finite theta values for {col}")

        vmin = float(np.min(merged))
        vmax = float(np.max(merged))

        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            raise RuntimeError(f"Invalid theta range for {col}: {vmin}, {vmax}")

        edges = np.linspace(vmin, vmax, theta_bins + 1)
        edges_by_projection[projection] = edges

        log(
            f"{projection}: common theta edges from {col}: "
            f"{theta_bins} bins, {vmin:.6g} to {vmax:.6g} deg"
        )

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
    use_width_weighting: bool,
) -> float:
    if not use_width_weighting:
        return 1.0

    w = 1.0

    for axis in widths:
        val = row[f"_width_{axis}"]

        if not np.isfinite(val) or val <= 0.0:
            return math.nan

        w *= float(val)

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

    vals = group[avg_col].apply(parse_scalar).astype(float)
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return midpoint

    return float(np.mean(vals))


def combine_columns_in_row(
    row: pd.Series,
    columns: Sequence[str],
    combine_mode: str,
) -> Tuple[float, float, float]:
    vals = []
    stat2 = 0.0
    sys2 = 0.0

    for col in columns:
        val, stat, sys = parse_tuple3(row[col])

        if not np.isfinite(val):
            continue

        vals.append(float(val))

        if np.isfinite(stat):
            stat2 += float(stat) * float(stat)

        if np.isfinite(sys):
            sys2 += float(sys) * float(sys)

    if not vals:
        return math.nan, math.nan, math.nan

    if combine_mode == "mean":
        n = float(len(vals))
        return (
            float(np.mean(vals)),
            math.sqrt(stat2) / n,
            math.sqrt(sys2) / n,
        )

    return (
        float(np.sum(vals)),
        math.sqrt(stat2),
        math.sqrt(sys2),
    )


def metric_points(
    df: pd.DataFrame,
    label: str,
    projection: str,
    period: str,
    metric: MetricSpec,
    no_width_weighting: bool,
) -> List[Point]:
    info = PROJECTIONS[projection]
    cols = metric_columns_for_label_period(df, label, period, metric)

    if not cols:
        return []

    points: List[Point] = []
    group_cols = [str(info["min_col"]), str(info["max_col"])]

    sorted_df = df.sort_values(group_cols)

    for _, group in sorted_df.groupby(group_cols, sort=True, dropna=True):
        y = 0.0
        stat2 = 0.0
        sys2 = 0.0
        vals_for_mean = []
        stat2_for_mean = 0.0
        sys2_for_mean = 0.0
        n = 0

        for _, row in group.iterrows():
            combine_mode = "mean" if metric.aggregate == "mean" else "sum"
            val, stat, sys = combine_columns_in_row(row, cols, combine_mode=combine_mode)

            if not np.isfinite(val):
                continue

            stat = 0.0 if not np.isfinite(stat) else stat
            sys = 0.0 if not np.isfinite(sys) else sys

            if metric.aggregate == "mean":
                vals_for_mean.append(val)
                stat2_for_mean += stat * stat
                sys2_for_mean += sys * sys
                n += 1
                continue

            use_width_weighting = metric.width_weighted and not no_width_weighting
            w = row_weight(
                row=row,
                widths=list(info["integrate_widths"]),
                use_width_weighting=use_width_weighting,
            )

            if not np.isfinite(w):
                continue

            y += val * w
            stat2 += (stat * w) ** 2
            sys2 += (sys * w) ** 2
            n += 1

        if n == 0:
            continue

        if metric.aggregate == "mean":
            y = float(np.mean(vals_for_mean))
            stat = math.sqrt(stat2_for_mean) / float(n)
            sys = math.sqrt(sys2_for_mean) / float(n)
        else:
            stat = math.sqrt(stat2)
            sys = math.sqrt(sys2)

        key = (
            float(group[str(info["min_col"])].iloc[0]),
            float(group[str(info["max_col"])].iloc[0]),
        )

        points.append(
            Point(
                key=key,
                x=x_position(group, info, period),
                y=y,
                stat=stat,
                sys=sys,
            )
        )

    points.sort(key=lambda p: p.x)

    return points


def integrated_points(
    df: pd.DataFrame,
    projection: str,
    period: str,
    template: str,
    no_width_weighting: bool,
) -> List[Point]:
    metric = make_template_metric(
        tag="legacy",
        title="legacy cross section",
        ylabel="cross section",
        template=template,
        aggregate="sum",
        width_weighted=True,
        required=False,
    )
    return metric_points(
        df=df,
        label="",
        projection=projection,
        period=period,
        metric=metric,
        no_width_weighting=no_width_weighting,
    )


def precompute_all_points_for_metric(
    dfs: Dict[str, pd.DataFrame],
    projection: str,
    metric: MetricSpec,
    no_width_weighting: bool,
) -> Dict[str, Dict[str, List[Point]]]:
    return {
        label: {
            period: metric_points(
                df=df,
                label=label,
                projection=projection,
                period=period,
                metric=metric,
                no_width_weighting=no_width_weighting,
            )
            for period in RUN_PERIODS
        }
        for label, df in dfs.items()
    }


def precompute_all_points(
    dfs: Dict[str, pd.DataFrame],
    projection: str,
    template: str,
    no_width_weighting: bool,
) -> Dict[str, Dict[str, List[Point]]]:
    metric = make_template_metric(
        tag="legacy",
        title="legacy cross section",
        ylabel="cross section",
        template=template,
        aggregate="sum",
        width_weighted=True,
        required=False,
    )

    return precompute_all_points_for_metric(
        dfs=dfs,
        projection=projection,
        metric=metric,
        no_width_weighting=no_width_weighting,
    )


def average_points_by_key(
    points_by_label: Dict[str, List[Point]],
) -> Dict[Tuple[float, float], Point]:
    keys = sorted(
        {
            p.key
            for points in points_by_label.values()
            for p in points
        }
    )

    avg_by_key: Dict[Tuple[float, float], Point] = {}

    for key in keys:
        used = []

        for points in points_by_label.values():
            lookup = {p.key: p for p in points}
            p = lookup.get(key)

            if p is None:
                continue

            if not np.isfinite(p.y):
                continue

            used.append(p)

        if len(used) == 0:
            continue

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

        if not np.isfinite(p.y) or not np.isfinite(avg.y) or avg.y == 0.0:
            continue

        r = p.y / avg.y
        rel2 = 0.0

        if p.y != 0.0 and np.isfinite(p.stat):
            rel2 += (p.stat / p.y) ** 2

        if avg.y != 0.0 and np.isfinite(avg.stat):
            rel2 += (avg.stat / avg.y) ** 2

        ratios.append(
            Point(
                key=p.key,
                x=p.x,
                y=r,
                stat=abs(r) * math.sqrt(rel2),
                sys=0.0,
            )
        )

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


def chi2_ndf_for_ratio_points(points_by_label: Dict[str, List[Point]]) -> Tuple[float, int]:
    chi2 = 0.0
    n = 0
    keys = set()

    for points in points_by_label.values():
        for p in points:
            keys.add(p.key)

            if not np.isfinite(p.y) or not np.isfinite(p.stat) or p.stat <= 0.0:
                continue

            chi2 += ((p.y - 1.0) / p.stat) ** 2
            n += 1

    ndf = n - len(keys)

    if ndf <= 0:
        return math.nan, 0

    return chi2, ndf


def chi2_label_text(points_by_label: Dict[str, List[Point]]) -> str:
    chi2, ndf = chi2_ndf_for_ratio_points(points_by_label)

    if not np.isfinite(chi2) or ndf <= 0:
        return r"$\chi^2/\mathrm{ndf}$ = n/a"

    return rf"$\chi^2/\mathrm{{ndf}}$ = {chi2:.1f}/{ndf:d} = {chi2 / ndf:.2f}"


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

    return math.sqrt(stat * stat + (frac * abs(point.y)) ** 2)


def collect_y_values(points_by_label: Dict[str, List[Point]]) -> List[float]:
    values: List[float] = []

    for points in points_by_label.values():
        for p in points:
            if np.isfinite(p.y):
                values.append(float(p.y))

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

    return values


def collect_all_ratio_y_values(
    ratio_points_by_period_label: Dict[str, Dict[str, List[Point]]],
) -> List[float]:
    values: List[float] = []

    for ratio_points_by_label in ratio_points_by_period_label.values():
        values.extend(collect_y_values(ratio_points_by_label))

    return values


def dynamic_absolute_ylim_from_values(values: List[float]) -> Tuple[float, float]:
    if len(values) == 0:
        return 0.0, 0.01

    finite = [float(v) for v in values if np.isfinite(v)]

    if not finite:
        return 0.0, 0.01

    ymin_data = min(finite)
    ymax_data = max(finite)

    if ymax_data <= 0.0 and ymin_data < 0.0:
        span = ymax_data - ymin_data
        pad = 0.1 * span if span > 0.0 else 0.01
        return ymin_data - pad, ymax_data + pad

    ymax = math.ceil(ymax_data / 0.01) * 0.01

    if ymax <= 0.0:
        ymax = 0.01

    ymin = 0.0 if ymin_data >= 0.0 else ymin_data - 0.1 * (ymax - ymin_data)

    return ymin, ymax


def dynamic_ratio_ylim_from_values(values: List[float]) -> Tuple[float, float]:
    return DEFAULT_RATIO_YMIN, DEFAULT_RATIO_YMAX


def plot_points(
    ax,
    points: List[Point],
    label: str,
    include_bin_to_bin_sys: bool,
    frac: float,
) -> None:
    if not points:
        return

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

def make_metric_projection_canvas(
    dfs: Dict[str, pd.DataFrame],
    projection: str,
    output_dir: str,
    output_prefix: str,
    metric: MetricSpec,
    no_width_weighting: bool,
    include_bin_to_bin_sys: bool,
    frac: float,
    comparison_name: str,
    legacy_filename: bool = False,
) -> None:
    info = PROJECTIONS[projection]

    with Timer(f"integrating {metric.tag} {projection}"):
        points_by_label_period = precompute_all_points_for_metric(
            dfs=dfs,
            projection=projection,
            metric=metric,
            no_width_weighting=no_width_weighting,
        )

    n_points = sum(
        len(points)
        for by_period in points_by_label_period.values()
        for points in by_period.values()
    )

    if n_points == 0:
        log(f"Skipping {metric.tag} {projection}: no finite points / no matching columns")
        return

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

    fig.suptitle(f"{metric.title}: {info['title']} by {comparison_name}", fontsize=16)

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

        ratio_points_by_label = ratio_points_by_period_label[period]

        for label, ratios in ratio_points_by_label.items():
            plot_ratio_points(
                ax=ax_ratio,
                points=ratios,
                label=label,
                use_period_style=False,
            )

        add_chi2_box(
            ax=ax_top,
            text=chi2_label_text(ratio_points_by_label),
        )

        ax_top.set_title(period, fontsize=11)
        ax_top.set_ylabel(str(metric.ylabel))
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

    os.makedirs(output_dir, exist_ok=True)

    if legacy_filename:
        filename = f"{output_prefix}_{info['outfile_tag']}_comparison.png"
    else:
        filename = f"{output_prefix}_{metric.tag}_{info['outfile_tag']}_comparison.png"

    out = os.path.join(output_dir, filename)

    fig.savefig(out, dpi=150)
    plt.close(fig)

    log(f"Wrote: {out}")


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
    ylabel = {
        "xB": r"$d\sigma/dx_B$ (pb)",
        "Q2": r"$d\sigma/dQ^2$ (pb/GeV$^2$)",
        "t": r"$d\sigma/d|t|$ (pb/GeV$^2$)",
        "e_theta": r"$\sigma_{\mathrm{int}}$ (pb)",
        "p_theta": r"$\sigma_{\mathrm{int}}$ (pb)",
        "g_theta": r"$\sigma_{\mathrm{int}}$ (pb)",
    }.get(projection, "cross section")

    metric = make_template_metric(
        tag="normed_cross_section_epg_unpol",
        title="normed DVCS cross section",
        ylabel=ylabel,
        template=template,
        aggregate="sum",
        width_weighted=True,
        required=False,
    )

    make_metric_projection_canvas(
        dfs=dfs,
        projection=projection,
        output_dir=output_dir,
        output_prefix=output_prefix,
        metric=metric,
        no_width_weighting=no_width_weighting,
        include_bin_to_bin_sys=include_bin_to_bin_sys,
        frac=frac,
        comparison_name=comparison_name,
        legacy_filename=True,
    )


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


# -----------------------------------------------------------------------------
# Scalar period-level canvases.
# -----------------------------------------------------------------------------

def scalar_value_for_df_label_period(
    df: pd.DataFrame,
    label: str,
    period: str,
    metric: MetricSpec,
) -> Tuple[float, float, float]:
    cols = metric_columns_for_label_period(df, label, period, metric)

    if not cols:
        return math.nan, math.nan, math.nan

    vals = []
    stat2 = 0.0
    sys2 = 0.0

    for col in cols:
        for cell in df[col].iloc[: min(20, len(df))]:
            val, stat, sys = parse_tuple3(cell)

            if np.isfinite(val):
                vals.append(val)

                if np.isfinite(stat):
                    stat2 += stat * stat

                if np.isfinite(sys):
                    sys2 += sys * sys

                break

    if not vals:
        return math.nan, math.nan, math.nan

    n = float(len(vals))

    return float(np.mean(vals)), math.sqrt(stat2) / n, math.sqrt(sys2) / n


def make_scalar_metric_canvas(
    dfs: Dict[str, pd.DataFrame],
    output_dir: str,
    output_prefix: str,
    metric: MetricSpec,
    comparison_name: str,
) -> None:
    labels = list(dfs.keys())

    fig = plt.figure(figsize=(18.5, 10.2))
    outer = fig.add_gridspec(
        2,
        3,
        left=0.055,
        right=0.992,
        bottom=0.070,
        top=0.910,
        wspace=0.205,
        hspace=0.300,
    )

    fig.suptitle(f"{metric.title} by {comparison_name}", fontsize=16)

    any_points = False

    for panel_index, period in enumerate(RUN_PERIOD_PANEL_ORDER):
        ax = fig.add_subplot(outer[panel_index // 3, panel_index % 3])

        if period is None:
            ax.axis("off")
            continue

        x = np.arange(1, len(labels) + 1, dtype=float)
        y = []
        ey = []

        for label in labels:
            val, stat, _ = scalar_value_for_df_label_period(
                df=dfs[label],
                label=label,
                period=period,
                metric=metric,
            )
            y.append(val)
            ey.append(stat if np.isfinite(stat) and stat >= 0.0 else 0.0)

        finite_y = [v for v in y if np.isfinite(v)]

        if finite_y:
            any_points = True

        ax.errorbar(
            x,
            y,
            yerr=ey,
            fmt="o",
            capsize=3,
            linewidth=1.2,
            markersize=6.5,
        )

        if finite_y:
            mean = float(np.mean(finite_y))
            ax.axhline(mean, color="0.35", linestyle="--", linewidth=1.0)

        ax.set_title(period, fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel(metric.ylabel)
        ax.grid(True, alpha=0.25)

    if not any_points:
        plt.close(fig)
        log(f"Skipping scalar {metric.tag}: no finite points / no matching columns")
        return

    os.makedirs(output_dir, exist_ok=True)

    out = os.path.join(
        output_dir,
        f"{output_prefix}_{metric.tag}_comparison.png",
    )

    fig.savefig(out, dpi=150)
    plt.close(fig)
    log(f"Wrote: {out}")


# -----------------------------------------------------------------------------
# Expanded step diagnostics.
# -----------------------------------------------------------------------------

def make_all_step_diagnostic_canvases(
    dfs: Dict[str, pd.DataFrame],
    output_dir: str,
    output_prefix: str,
    no_width_weighting: bool,
    include_bin_to_bin_sys: bool,
    frac: float,
    comparison_name: str,
    metric_tags: Optional[List[str]] = None,
) -> None:
    os.makedirs(output_dir, exist_ok=True)

    specs = all_step_metric_specs()

    if metric_tags:
        wanted = set(metric_tags)
        specs = [s for s in specs if s.tag in wanted]
        missing = wanted - {s.tag for s in specs}

        if missing:
            raise ValueError(
                "Unknown diagnostic metric tag(s): "
                + ", ".join(sorted(missing))
            )

    # Do not duplicate the legacy normed-cross-section canvases here; the
    # original make_all_projection_canvases() still produces those exact files.
    specs = [s for s in specs if s.tag != "normed_cross_section_epg_unpol"]

    diagnostic_dir = os.path.join(output_dir, "step_diagnostics")
    os.makedirs(diagnostic_dir, exist_ok=True)

    available_report_rows = []

    for metric in specs:
        if metric.aggregate == "scalar":
            with Timer(f"plotting scalar {metric.tag}"):
                make_scalar_metric_canvas(
                    dfs=dfs,
                    output_dir=diagnostic_dir,
                    output_prefix=output_prefix,
                    metric=metric,
                    comparison_name=comparison_name,
                )
            continue

        for projection in PROJECTION_ORDER:
            with Timer(f"plotting diagnostic {metric.tag} {projection}"):
                before_files = set(os.listdir(diagnostic_dir))
                make_metric_projection_canvas(
                    dfs=dfs,
                    projection=projection,
                    output_dir=diagnostic_dir,
                    output_prefix=output_prefix,
                    metric=metric,
                    no_width_weighting=no_width_weighting,
                    include_bin_to_bin_sys=include_bin_to_bin_sys,
                    frac=frac,
                    comparison_name=comparison_name,
                    legacy_filename=False,
                )
                after_files = set(os.listdir(diagnostic_dir))
                if after_files != before_files:
                    available_report_rows.append((metric.tag, projection))

    report = os.path.join(diagnostic_dir, f"{output_prefix}_step_diagnostics_index.txt")

    with open(report, "w") as fout:
        print("Step diagnostics written", file=fout)
        print("========================", file=fout)
        print("", file=fout)
        print(f"comparison_name = {comparison_name}", file=fout)
        print(f"ratio_y_range = {DEFAULT_RATIO_YMIN} to {DEFAULT_RATIO_YMAX}", file=fout)
        print("", file=fout)
        print("Available diagnostic metric tags:", file=fout)
        for metric in all_step_metric_specs():
            print(f"  {metric.tag:55s}  {metric.title}", file=fout)
        print("", file=fout)
        print("Canvases successfully attempted with finite points:", file=fout)
        for tag, projection in available_report_rows:
            print(f"  {tag:55s}  {projection}", file=fout)

    log(f"Wrote: {report}")


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


# -----------------------------------------------------------------------------
# CLI.
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare topology-resolved DVCS pass-2 quantities from separate CSVs."
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
        help="Use raw sums instead of bin-width weighted integrations for width-weighted metrics.",
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

    parser.add_argument(
        "--no-step-diagnostics",
        action="store_true",
        help="Disable expanded step-by-step diagnostic canvases.",
    )

    parser.add_argument(
        "--diagnostic-metrics",
        nargs="+",
        default=None,
        help=(
            "Optional subset of diagnostic metric tags to plot. "
            "Use the generated step_diagnostics_index.txt to see all tags."
        ),
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
    log("ratio y-axis range: 0 to 2")

    if "{period}" not in args.xs_template:
        raise ValueError("--xs-template must contain {period}")

    if args.theta_bins <= 0:
        raise ValueError("--theta-bins must be positive")

    if args.bin_to_bin_sys_frac < 0.0:
        raise ValueError("--bin-to-bin-sys-frac must be non-negative")

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

    if not args.no_step_diagnostics:
        make_all_step_diagnostic_canvases(
            dfs=dfs,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
            no_width_weighting=args.no_width_weighting,
            include_bin_to_bin_sys=args.include_bin_to_bin_sys,
            frac=args.bin_to_bin_sys_frac,
            comparison_name=args.comparison_name,
            metric_tags=args.diagnostic_metrics,
        )

    if not args.no_topology_period_ratio_plots:
        make_all_topology_period_ratio_canvases(
            dfs=dfs,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
            template=args.xs_template,
            no_width_weighting=args.no_width_weighting,
        )

    log(f"TOTAL RUNTIME: {time.perf_counter() - t0:.3f} s")
    log("Done.")


if __name__ == "__main__":
    main()
