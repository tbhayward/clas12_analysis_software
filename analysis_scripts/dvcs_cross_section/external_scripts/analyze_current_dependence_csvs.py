#!/usr/bin/env python3
"""
analyze_current_dependence_csvs.py

Compare current-dependence summary CSVs produced by current_dependence.cpp.

Expected input format is the period_summary.csv written by current_dependence.cpp,
with columns like:

  period
  data_factor
  data_factor_stat
  mc_factor
  mc_factor_stat
  data_m
  data_b
  data_sm
  data_sb
  data_cov_mb
  data_slope_percent_per_nA
  data_slope_percent_per_nA_stat
  mc_m
  mc_b
  mc_sm
  mc_sb
  mc_cov_mb

Usage for six electron/photon/proton FD sector summaries:

  python3 analyze_current_dependence_csvs.py \
    sec1/period_summary.csv \
    sec2/period_summary.csv \
    sec3/period_summary.csv \
    sec4/period_summary.csv \
    sec5/period_summary.csv \
    sec6/period_summary.csv \
    --labels S1 S2 S3 S4 S5 S6

Outputs by default:

  output/current_dependence_sector_diagnostics/current_dependence_values_long.csv
  output/current_dependence_sector_diagnostics/current_dependence_values_wide.csv
  output/current_dependence_sector_diagnostics/current_dependence_diagnostics.csv
  output/current_dependence_sector_diagnostics/current_dependence_report.txt

and PNG plots in the same directory.

Plot behavior:

  - Individual DATA-only and MC-only value canvases.
  - Individual DATA-only and MC-only ratio-to-sector-mean canvases.
  - DATA+MC overlay value canvases.
  - DATA+MC overlay ratio-to-sector-mean canvases.
  - Pull heatmaps.
  - Spread-summary canvas.

By default, the value canvases use local per-panel y-axis ranges so each period
is visually zoomed to its own sector/input spread. Use --shared-y-range to force
one common y-axis range across all panels in a given canvas.

Statistical-consistency convention:

  weighted mean:
    xbar = sum_i x_i / sigma_i^2 / sum_i 1 / sigma_i^2

  chi2:
    chi2 = sum_i (x_i - xbar)^2 / sigma_i^2
    ndf  = N - 1

This assumes the input uncertainties are independent. That is probably not
strictly true for sector splits from the same data set, so chi2/ndf should be
treated as a diagnostic tension metric, not a formal hypothesis test.
"""

import argparse
import math
import os
import re
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt


PERIOD_ORDER = [
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
]

PANEL_PERIOD_ORDER = [
    "Sp18 Inb",
    "Sp18 Out",
    None,
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
]

SAMPLE_ORDER = [
    "data",
    "mc",
]

QUANTITY_ORDER = [
    "current_efficiency_factor",
    "slope_percent_per_nA",
]

DEFAULT_OUTPUT_DIR = "output/current_dependence_sector_diagnostics"
DEFAULT_OUTPUT_PREFIX = "current_dependence"

FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)

PLOT_DPI = 200

MARKERS = [
    "o",
    "s",
    "^",
    "v",
    "D",
    "P",
    "X",
    "*",
    "<",
    ">",
]

LINESTYLES = [
    "-",
    "--",
    "-.",
    ":",
    "-",
    "--",
    "-.",
    ":",
    "-",
    "--",
]

QUANTITY_LABEL = {
    "current_efficiency_factor": "Current-efficiency factor",
    "slope_percent_per_nA": "Current-dependence slope (%/nA)",
}

QUANTITY_SHORT = {
    "current_efficiency_factor": "factor",
    "slope_percent_per_nA": "slope",
}

SAMPLE_STYLE = {
    "data": {
        "label": "DATA",
        "color": "tab:blue",
        "marker": "o",
        "offset": -0.08,
    },
    "mc": {
        "label": "MC",
        "color": "tab:orange",
        "marker": "s",
        "offset": 0.08,
    },
}


@dataclass
class ValueEntry:
    label: str
    csv_file: str
    period: str
    sample: str
    quantity: str
    value: float
    stat: float


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
# Parsing helpers.
# -----------------------------------------------------------------------------

def parse_first_number(value) -> float:
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


def require_columns(df: pd.DataFrame, columns: Iterable[str], path: str) -> None:
    missing = [c for c in columns if c not in df.columns]

    if missing:
        print()
        print(f"Missing required columns in: {path}")
        for c in missing:
            print(f"  {c}")

        print()
        print("Available columns:")
        for c in df.columns:
            print(f"  {c}")

        raise KeyError(f"{path} is missing {len(missing)} required columns.")


def safe_filename_token(text: str) -> str:
    out = text.strip().lower()
    out = out.replace("/", "_")
    out = out.replace("\\", "_")
    out = out.replace(" ", "_")
    out = out.replace("-", "_")
    out = out.replace("%", "percent")
    out = out.replace("(", "")
    out = out.replace(")", "")

    while "__" in out:
        out = out.replace("__", "_")

    return out.strip("_")


def finite_float(value, default: float = math.nan) -> float:
    try:
        v = float(value)
    except Exception:
        return default

    return v if np.isfinite(v) else default


# -----------------------------------------------------------------------------
# Current-dependence quantity helpers.
# -----------------------------------------------------------------------------

def percent_slope(m: float, b: float) -> float:
    if not np.isfinite(m) or not np.isfinite(b) or b == 0.0:
        return math.nan

    return 100.0 * m / b


def percent_slope_err(
    m: float,
    b: float,
    sm: float,
    sb: float,
    cov_mb: float,
) -> float:
    if (
        not np.isfinite(m)
        or not np.isfinite(b)
        or not np.isfinite(sm)
        or not np.isfinite(sb)
        or not np.isfinite(cov_mb)
        or b == 0.0
    ):
        return math.nan

    dm = 100.0 / b
    db = -100.0 * m / (b * b)

    var = (
        dm * dm * sm * sm
        + db * db * sb * sb
        + 2.0 * dm * db * cov_mb
    )

    if var < 0.0 and abs(var) < 1.0e-15:
        var = 0.0

    if var < 0.0:
        return math.nan

    return math.sqrt(var)


def get_scalar(row: pd.Series, column: str) -> float:
    if column not in row.index:
        return math.nan

    return parse_first_number(row[column])


def get_data_slope_value_and_error(row: pd.Series) -> Tuple[float, float]:
    value = get_scalar(row, "data_slope_percent_per_nA")
    stat = get_scalar(row, "data_slope_percent_per_nA_stat")

    if np.isfinite(value):
        return value, stat

    m = get_scalar(row, "data_m")
    b = get_scalar(row, "data_b")
    sm = get_scalar(row, "data_sm")
    sb = get_scalar(row, "data_sb")
    cov = get_scalar(row, "data_cov_mb")

    return percent_slope(m, b), percent_slope_err(m, b, sm, sb, cov)


def get_mc_slope_value_and_error(row: pd.Series) -> Tuple[float, float]:
    if "mc_slope_percent_per_nA" in row.index:
        value = get_scalar(row, "mc_slope_percent_per_nA")
        stat = get_scalar(row, "mc_slope_percent_per_nA_stat")

        if np.isfinite(value):
            return value, stat

    m = get_scalar(row, "mc_m")
    b = get_scalar(row, "mc_b")
    sm = get_scalar(row, "mc_sm")
    sb = get_scalar(row, "mc_sb")
    cov = get_scalar(row, "mc_cov_mb")

    return percent_slope(m, b), percent_slope_err(m, b, sm, sb, cov)


# -----------------------------------------------------------------------------
# Read inputs.
# -----------------------------------------------------------------------------

def read_one_summary_csv(
    path: str,
    label: str,
) -> List[ValueEntry]:
    required = [
        "period",
        "data_factor",
        "data_factor_stat",
        "mc_factor",
        "mc_factor_stat",
        "data_m",
        "data_b",
        "data_sm",
        "data_sb",
        "data_cov_mb",
        "mc_m",
        "mc_b",
        "mc_sm",
        "mc_sb",
        "mc_cov_mb",
    ]

    with Timer(f"reading {label}: {path}"):
        df = pd.read_csv(path, low_memory=False)

    require_columns(df, required, path)

    log(f"{label}: {df.shape[0]} rows x {df.shape[1]} columns")

    entries: List[ValueEntry] = []

    for _, row in df.iterrows():
        period = str(row["period"]).strip()

        if period == "" or period.lower() == "nan":
            continue

        data_factor = get_scalar(row, "data_factor")
        data_factor_stat = get_scalar(row, "data_factor_stat")

        mc_factor = get_scalar(row, "mc_factor")
        mc_factor_stat = get_scalar(row, "mc_factor_stat")

        data_slope, data_slope_stat = get_data_slope_value_and_error(row)
        mc_slope, mc_slope_stat = get_mc_slope_value_and_error(row)

        entries.append(
            ValueEntry(
                label=label,
                csv_file=path,
                period=period,
                sample="data",
                quantity="current_efficiency_factor",
                value=data_factor,
                stat=data_factor_stat,
            )
        )

        entries.append(
            ValueEntry(
                label=label,
                csv_file=path,
                period=period,
                sample="mc",
                quantity="current_efficiency_factor",
                value=mc_factor,
                stat=mc_factor_stat,
            )
        )

        entries.append(
            ValueEntry(
                label=label,
                csv_file=path,
                period=period,
                sample="data",
                quantity="slope_percent_per_nA",
                value=data_slope,
                stat=data_slope_stat,
            )
        )

        entries.append(
            ValueEntry(
                label=label,
                csv_file=path,
                period=period,
                sample="mc",
                quantity="slope_percent_per_nA",
                value=mc_slope,
                stat=mc_slope_stat,
            )
        )

    return entries


def read_all_summary_csvs(
    paths: List[str],
    labels: List[str],
) -> List[ValueEntry]:
    all_entries: List[ValueEntry] = []

    for path, label in zip(paths, labels):
        entries = read_one_summary_csv(path=path, label=label)
        all_entries.extend(entries)

    return all_entries


# -----------------------------------------------------------------------------
# Statistical diagnostics.
# -----------------------------------------------------------------------------

def chi2_pvalue(chi2: float, ndf: int) -> float:
    if not np.isfinite(chi2) or ndf <= 0:
        return math.nan

    try:
        from scipy.stats import chi2 as scipy_chi2
        return float(scipy_chi2.sf(chi2, ndf))
    except Exception:
        return math.nan


def weighted_mean_and_chi2(
    values: np.ndarray,
    stats: np.ndarray,
) -> Tuple[float, float, float, int, float, int]:
    """
    Return:

      weighted_mean, weighted_mean_stat, chi2, ndf, p_value, n_weighted

    Uses only entries with finite value and finite positive uncertainty.
    """

    mask = (
        np.isfinite(values)
        & np.isfinite(stats)
        & (stats > 0.0)
    )

    v = values[mask]
    s = stats[mask]

    if len(v) == 0:
        return math.nan, math.nan, math.nan, 0, math.nan, 0

    w = 1.0 / (s * s)
    wsum = float(np.sum(w))

    if not np.isfinite(wsum) or wsum <= 0.0:
        return math.nan, math.nan, math.nan, 0, math.nan, 0

    mean = float(np.sum(w * v) / wsum)
    mean_stat = float(math.sqrt(1.0 / wsum))

    chi2 = float(np.sum(((v - mean) / s) ** 2))
    ndf = int(len(v) - 1)
    pval = chi2_pvalue(chi2, ndf)

    return mean, mean_stat, chi2, ndf, pval, int(len(v))


def max_pairwise_pull(
    labels: List[str],
    values: np.ndarray,
    stats: np.ndarray,
) -> Tuple[float, str, str]:
    best_pull = math.nan
    best_a = ""
    best_b = ""

    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            xi = values[i]
            xj = values[j]
            si = stats[i]
            sj = stats[j]

            if not (
                np.isfinite(xi)
                and np.isfinite(xj)
                and np.isfinite(si)
                and np.isfinite(sj)
                and si >= 0.0
                and sj >= 0.0
            ):
                continue

            denom = math.sqrt(si * si + sj * sj)

            if denom <= 0.0:
                continue

            pull = abs(xi - xj) / denom

            if not np.isfinite(best_pull) or pull > best_pull:
                best_pull = float(pull)
                best_a = labels[i]
                best_b = labels[j]

    return best_pull, best_a, best_b


def empty_diagnostics() -> Dict[str, object]:
    return {
        "n_finite": 0,
        "n_weighted": 0,
        "mean": math.nan,
        "std": math.nan,
        "rms_about_mean": math.nan,
        "rms_over_mean_percent": math.nan,
        "min_label": "",
        "min_value": math.nan,
        "max_label": "",
        "max_value": math.nan,
        "max_over_min": math.nan,
        "weighted_mean": math.nan,
        "weighted_mean_stat": math.nan,
        "reference_mean": math.nan,
        "reference_mean_type": "",
        "chi2": math.nan,
        "ndf": 0,
        "chi2_ndf": math.nan,
        "p_value": math.nan,
        "max_pair_pull": math.nan,
        "max_pair_pull_a": "",
        "max_pair_pull_b": "",
    }


def diagnostics_for_group(group: pd.DataFrame) -> Dict[str, object]:
    if len(group) == 0:
        return empty_diagnostics()

    finite = group[np.isfinite(group["value"].astype(float))].copy()

    if len(finite) == 0:
        return empty_diagnostics()

    values = finite["value"].astype(float).to_numpy()
    stats = finite["stat"].astype(float).to_numpy()
    labels = finite["label"].astype(str).tolist()

    mean = float(np.mean(values))
    std = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
    rms = float(math.sqrt(np.mean((values - mean) ** 2))) if len(values) > 0 else math.nan

    if np.isfinite(mean) and mean != 0.0:
        rms_over_mean_percent = 100.0 * rms / abs(mean)
    else:
        rms_over_mean_percent = math.nan

    imin = int(np.argmin(values))
    imax = int(np.argmax(values))

    min_value = float(values[imin])
    max_value = float(values[imax])
    min_label = labels[imin]
    max_label = labels[imax]

    if np.isfinite(min_value) and min_value != 0.0:
        max_over_min = max_value / min_value
    else:
        max_over_min = math.nan

    weighted_mean, weighted_mean_stat, chi2, ndf, p_value, n_weighted = weighted_mean_and_chi2(
        values=values,
        stats=stats,
    )

    if np.isfinite(weighted_mean):
        reference_mean = weighted_mean
        reference_mean_type = "weighted"
    else:
        reference_mean = mean
        reference_mean_type = "unweighted"

    chi2_ndf = chi2 / ndf if np.isfinite(chi2) and ndf > 0 else math.nan

    max_pull, max_pull_a, max_pull_b = max_pairwise_pull(
        labels=labels,
        values=values,
        stats=stats,
    )

    return {
        "n_finite": int(len(values)),
        "n_weighted": int(n_weighted),
        "mean": mean,
        "std": std,
        "rms_about_mean": rms,
        "rms_over_mean_percent": rms_over_mean_percent,
        "min_label": min_label,
        "min_value": min_value,
        "max_label": max_label,
        "max_value": max_value,
        "max_over_min": max_over_min,
        "weighted_mean": weighted_mean,
        "weighted_mean_stat": weighted_mean_stat,
        "reference_mean": reference_mean,
        "reference_mean_type": reference_mean_type,
        "chi2": chi2,
        "ndf": int(ndf),
        "chi2_ndf": chi2_ndf,
        "p_value": p_value,
        "max_pair_pull": max_pull,
        "max_pair_pull_a": max_pull_a,
        "max_pair_pull_b": max_pull_b,
    }


def build_values_dataframe(entries: List[ValueEntry]) -> pd.DataFrame:
    rows = []

    for e in entries:
        rows.append(
            {
                "label": e.label,
                "csv_file": e.csv_file,
                "period": e.period,
                "sample": e.sample,
                "quantity": e.quantity,
                "value": e.value,
                "stat": e.stat,
            }
        )

    return pd.DataFrame(rows)


def build_diagnostics_dataframe(values_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    available_periods = set(values_df["period"].astype(str))
    known_periods = [p for p in PERIOD_ORDER if p in available_periods]
    extra_periods = sorted(available_periods - set(known_periods))
    period_order = known_periods + extra_periods

    for period in period_order:
        for sample in SAMPLE_ORDER:
            for quantity in QUANTITY_ORDER:
                group = values_df[
                    (values_df["period"] == period)
                    & (values_df["sample"] == sample)
                    & (values_df["quantity"] == quantity)
                ]

                d = diagnostics_for_group(group)

                row = {
                    "period": period,
                    "sample": sample,
                    "quantity": quantity,
                }

                row.update(d)
                rows.append(row)

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Reporting.
# -----------------------------------------------------------------------------

def format_value_err(value: float, err: float) -> str:
    if not np.isfinite(value):
        return "nan"

    if np.isfinite(err) and err >= 0.0:
        return f"{value:.8g} +/- {err:.3g}"

    return f"{value:.8g}"


def print_group_list(
    values_df: pd.DataFrame,
    period: str,
    sample: str,
    quantity: str,
    labels: List[str],
    fout,
) -> None:
    group = values_df[
        (values_df["period"] == period)
        & (values_df["sample"] == sample)
        & (values_df["quantity"] == quantity)
    ]

    lookup = {
        str(row["label"]): (float(row["value"]), float(row["stat"]))
        for _, row in group.iterrows()
    }

    print(f"  {sample}, {quantity}:", file=fout)

    for label in labels:
        value, err = lookup.get(label, (math.nan, math.nan))
        print(f"    {label:>10s}: {format_value_err(value, err)}", file=fout)


def print_diagnostic_line(
    diagnostics_df: pd.DataFrame,
    period: str,
    sample: str,
    quantity: str,
    fout,
) -> None:
    row_df = diagnostics_df[
        (diagnostics_df["period"] == period)
        & (diagnostics_df["sample"] == sample)
        & (diagnostics_df["quantity"] == quantity)
    ]

    if len(row_df) == 0:
        return

    row = row_df.iloc[0]

    print(
        "    diagnostics: "
        f"N={int(row['n_finite'])}, "
        f"N_weighted={int(row['n_weighted'])}, "
        f"mean={row['mean']:.8g}, "
        f"std={row['std']:.8g}, "
        f"RMS/mean={row['rms_over_mean_percent']:.3g}%, "
        f"weighted_mean={row['weighted_mean']:.8g} +/- {row['weighted_mean_stat']:.3g}, "
        f"chi2/ndf={row['chi2']:.3g}/{int(row['ndf'])} = {row['chi2_ndf']:.3g}, "
        f"p={row['p_value']:.3g}, "
        f"max_pair_pull={row['max_pair_pull']:.3g} "
        f"({row['max_pair_pull_a']} vs {row['max_pair_pull_b']}), "
        f"min={row['min_label']}:{row['min_value']:.8g}, "
        f"max={row['max_label']}:{row['max_value']:.8g}, "
        f"max/min={row['max_over_min']:.8g}",
        file=fout,
    )


def write_text_report(
    path: str,
    values_df: pd.DataFrame,
    diagnostics_df: pd.DataFrame,
    labels: List[str],
) -> None:
    with open(path, "w") as fout:
        print("Current-dependence sector/version comparison", file=fout)
        print("============================================", file=fout)
        print("", file=fout)

        print("Input labels:", file=fout)
        for label in labels:
            print(f"  {label}", file=fout)

        print("", file=fout)

        available_periods = set(values_df["period"].astype(str))
        known_periods = [p for p in PERIOD_ORDER if p in available_periods]
        extra_periods = sorted(available_periods - set(known_periods))
        period_order = known_periods + extra_periods

        for period in period_order:
            print("", file=fout)
            print(period, file=fout)
            print("-" * len(period), file=fout)

            for quantity in QUANTITY_ORDER:
                for sample in SAMPLE_ORDER:
                    print_group_list(
                        values_df=values_df,
                        period=period,
                        sample=sample,
                        quantity=quantity,
                        labels=labels,
                        fout=fout,
                    )

                    print_diagnostic_line(
                        diagnostics_df=diagnostics_df,
                        period=period,
                        sample=sample,
                        quantity=quantity,
                        fout=fout,
                    )

                    print("", file=fout)


def print_console_summary(diagnostics_df: pd.DataFrame) -> None:
    log("Statistical consistency summary:")

    for _, row in diagnostics_df.iterrows():
        print(
            f"  {row['period']:8s}  "
            f"{row['sample']:4s}  "
            f"{row['quantity']:26s}  "
            f"N={int(row['n_finite']):2d}  "
            f"Nw={int(row['n_weighted']):2d}  "
            f"RMS/mean={row['rms_over_mean_percent']:8.3g}%  "
            f"chi2/ndf={row['chi2_ndf']:8.3g}  "
            f"max-pull={row['max_pair_pull']:8.3g} "
            f"({row['max_pair_pull_a']} vs {row['max_pair_pull_b']})"
        )


# -----------------------------------------------------------------------------
# Wide values.
# -----------------------------------------------------------------------------

def build_wide_values_dataframe(values_df: pd.DataFrame, labels: List[str]) -> pd.DataFrame:
    rows = []

    available_periods = set(values_df["period"].astype(str))
    known_periods = [p for p in PERIOD_ORDER if p in available_periods]
    extra_periods = sorted(available_periods - set(known_periods))
    period_order = known_periods + extra_periods

    for period in period_order:
        for sample in SAMPLE_ORDER:
            for quantity in QUANTITY_ORDER:
                group = values_df[
                    (values_df["period"] == period)
                    & (values_df["sample"] == sample)
                    & (values_df["quantity"] == quantity)
                ]

                lookup = {
                    str(row["label"]): (float(row["value"]), float(row["stat"]))
                    for _, row in group.iterrows()
                }

                row = {
                    "period": period,
                    "sample": sample,
                    "quantity": quantity,
                }

                for label in labels:
                    value, stat = lookup.get(label, (math.nan, math.nan))
                    row[label] = value
                    row[f"{label}_stat"] = stat

                rows.append(row)

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Plot helpers.
# -----------------------------------------------------------------------------

def label_style(label: str, index: int) -> Dict[str, object]:
    return {
        "marker": MARKERS[index % len(MARKERS)],
        "linestyle": LINESTYLES[index % len(LINESTYLES)],
    }


def finite_min_max(values: List[float], pad_fraction: float = 0.12) -> Tuple[float, float]:
    finite = [v for v in values if np.isfinite(v)]

    if not finite:
        return 0.0, 1.0

    vmin = min(finite)
    vmax = max(finite)

    if vmin == vmax:
        pad = pad_fraction * abs(vmax) if vmax != 0.0 else 1.0
        return vmin - pad, vmax + pad

    pad = pad_fraction * (vmax - vmin)
    return vmin - pad, vmax + pad


def finite_min_max_with_errors(
    values: List[float],
    errors: List[float],
    pad_fraction: float = 0.20,
    min_abs_pad_fraction: float = 0.002,
) -> Tuple[float, float]:
    lows = []
    highs = []

    for value, err in zip(values, errors):
        if not np.isfinite(value):
            continue

        if np.isfinite(err) and err >= 0.0:
            lows.append(value - err)
            highs.append(value + err)
        else:
            lows.append(value)
            highs.append(value)

    if not lows or not highs:
        return 0.0, 1.0

    vmin = min(lows)
    vmax = max(highs)

    center = 0.5 * (vmin + vmax)
    span = vmax - vmin

    if not np.isfinite(span) or span <= 0.0:
        abs_pad = min_abs_pad_fraction * abs(center) if center != 0.0 else 1.0
        return center - abs_pad, center + abs_pad

    pad = pad_fraction * span
    min_pad = min_abs_pad_fraction * abs(center) if center != 0.0 else 0.0
    pad = max(pad, min_pad)

    return vmin - pad, vmax + pad


def period_order_from_values(values_df: pd.DataFrame) -> List[str]:
    available_periods = set(values_df["period"].astype(str))
    known_periods = [p for p in PERIOD_ORDER if p in available_periods]
    extra_periods = sorted(available_periods - set(known_periods))
    return known_periods + extra_periods


def get_group(
    values_df: pd.DataFrame,
    period: str,
    sample: str,
    quantity: str,
) -> pd.DataFrame:
    return values_df[
        (values_df["period"] == period)
        & (values_df["sample"] == sample)
        & (values_df["quantity"] == quantity)
    ].copy()


def get_diag_row(
    diagnostics_df: pd.DataFrame,
    period: str,
    sample: str,
    quantity: str,
) -> pd.Series:
    rows = diagnostics_df[
        (diagnostics_df["period"] == period)
        & (diagnostics_df["sample"] == sample)
        & (diagnostics_df["quantity"] == quantity)
    ]

    if len(rows) == 0:
        return pd.Series(dtype=object)

    return rows.iloc[0]


def value_lookup_for_group(group: pd.DataFrame) -> Dict[str, Tuple[float, float]]:
    out = {}

    for _, row in group.iterrows():
        out[str(row["label"])] = (
            float(row["value"]),
            float(row["stat"]),
        )

    return out


def draw_mean_reference(
    ax,
    mean: float,
    mean_stat: float,
    label: str = "mean",
    color: str = "black",
    alpha: float = 0.15,
) -> None:
    if not np.isfinite(mean):
        return

    ax.axhline(mean, linewidth=1.2, linestyle="--", color=color, label=label)

    if np.isfinite(mean_stat) and mean_stat > 0.0:
        ax.axhspan(mean - mean_stat, mean + mean_stat, alpha=alpha, color=color)


def collect_values_and_errors_for_period(
    values_df: pd.DataFrame,
    period: str,
    samples: List[str],
    quantity: str,
) -> Tuple[List[float], List[float]]:
    vals = []
    errs = []

    for sample in samples:
        group = get_group(values_df, period, sample, quantity)

        for _, row in group.iterrows():
            value = finite_float(row["value"])
            stat = finite_float(row["stat"], default=0.0)

            if np.isfinite(value):
                vals.append(value)
                errs.append(stat if np.isfinite(stat) and stat >= 0.0 else 0.0)

    return vals, errs


def collect_values_and_errors_for_canvas(
    values_df: pd.DataFrame,
    samples: List[str],
    quantity: str,
) -> Tuple[List[float], List[float]]:
    vals = []
    errs = []

    for period in period_order_from_values(values_df):
        period_vals, period_errs = collect_values_and_errors_for_period(
            values_df=values_df,
            period=period,
            samples=samples,
            quantity=quantity,
        )
        vals.extend(period_vals)
        errs.extend(period_errs)

    return vals, errs


def get_panel_ylim(
    values_df: pd.DataFrame,
    period: str,
    samples: List[str],
    quantity: str,
    shared_y_range: bool,
    shared_ylim: Tuple[float, float],
    y_padding_fraction: float,
) -> Tuple[float, float]:
    if shared_y_range:
        return shared_ylim

    vals, errs = collect_values_and_errors_for_period(
        values_df=values_df,
        period=period,
        samples=samples,
        quantity=quantity,
    )

    return finite_min_max_with_errors(
        values=vals,
        errors=errs,
        pad_fraction=y_padding_fraction,
    )


# -----------------------------------------------------------------------------
# Individual DATA-only or MC-only value plots.
# -----------------------------------------------------------------------------

def plot_values_by_period(
    values_df: pd.DataFrame,
    diagnostics_df: pd.DataFrame,
    labels: List[str],
    sample: str,
    quantity: str,
    output_dir: str,
    output_prefix: str,
    shared_y_range: bool,
    y_padding_fraction: float,
) -> str:
    fig = plt.figure(figsize=(18.5, 10.2))
    gs = fig.add_gridspec(
        2,
        3,
        left=0.055,
        right=0.990,
        bottom=0.075,
        top=0.900,
        wspace=0.220,
        hspace=0.300,
    )

    shared_vals, shared_errs = collect_values_and_errors_for_canvas(
        values_df=values_df,
        samples=[sample],
        quantity=quantity,
    )

    shared_ylim = finite_min_max_with_errors(
        values=shared_vals,
        errors=shared_errs,
        pad_fraction=y_padding_fraction,
    )

    for ipanel, period in enumerate(PANEL_PERIOD_ORDER):
        ax = fig.add_subplot(gs[ipanel // 3, ipanel % 3])

        if period is None:
            ax.axis("off")
            continue

        group = get_group(values_df, period, sample, quantity)
        lookup = value_lookup_for_group(group)
        diag = get_diag_row(diagnostics_df, period, sample, quantity)

        x = np.arange(1, len(labels) + 1, dtype=float)

        y = []
        ey = []

        for label in labels:
            value, stat = lookup.get(label, (math.nan, math.nan))
            y.append(value)
            ey.append(stat if np.isfinite(stat) and stat >= 0.0 else 0.0)

        ax.errorbar(
            x,
            y,
            yerr=ey,
            fmt=SAMPLE_STYLE.get(sample, {}).get("marker", "o"),
            color=SAMPLE_STYLE.get(sample, {}).get("color", None),
            capsize=3,
            linewidth=1.2,
            markersize=6.5,
        )

        ref_mean = finite_float(diag.get("reference_mean", math.nan))
        ref_type = str(diag.get("reference_mean_type", "mean"))
        ref_stat = finite_float(diag.get("weighted_mean_stat", math.nan))

        if ref_type != "weighted":
            ref_stat = math.nan

        draw_mean_reference(
            ax,
            ref_mean,
            ref_stat,
            label=ref_type,
            color="black",
            alpha=0.15,
        )

        chi2_ndf = finite_float(diag.get("chi2_ndf", math.nan))
        rms_percent = finite_float(diag.get("rms_over_mean_percent", math.nan))
        max_pull = finite_float(diag.get("max_pair_pull", math.nan))

        ax.set_title(
            f"{period}\n"
            f"RMS/mean={rms_percent:.3g}%, "
            f"chi2/ndf={chi2_ndf:.3g}, "
            f"max pull={max_pull:.3g}",
            fontsize=11,
        )

        ax.set_xticks(x)
        ax.set_xticklabels(labels)

        ymin, ymax = get_panel_ylim(
            values_df=values_df,
            period=period,
            samples=[sample],
            quantity=quantity,
            shared_y_range=shared_y_range,
            shared_ylim=shared_ylim,
            y_padding_fraction=y_padding_fraction,
        )
        ax.set_ylim(ymin, ymax)

        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=10)

        if ipanel % 3 == 0:
            ax.set_ylabel(QUANTITY_LABEL.get(quantity, quantity), fontsize=11)

        if ipanel // 3 == 1:
            ax.set_xlabel("Sector/input label", fontsize=11)

        if ipanel == 0:
            ax.legend(fontsize=9, loc="best")

    y_mode = "shared y" if shared_y_range else "local y"
    fig.suptitle(
        f"{sample.upper()} {QUANTITY_LABEL.get(quantity, quantity)} by sector/input ({y_mode})",
        fontsize=16,
        y=0.975,
    )

    suffix = "shared_y" if shared_y_range else "local_y"

    outfile = os.path.join(
        output_dir,
        f"{output_prefix}_values_by_period_{sample}_{safe_filename_token(quantity)}_{suffix}.png",
    )

    fig.savefig(outfile, dpi=PLOT_DPI)
    plt.close(fig)

    return outfile


# -----------------------------------------------------------------------------
# Combined DATA+MC value plots.
# -----------------------------------------------------------------------------

def plot_data_mc_values_by_period(
    values_df: pd.DataFrame,
    diagnostics_df: pd.DataFrame,
    labels: List[str],
    quantity: str,
    output_dir: str,
    output_prefix: str,
    shared_y_range: bool,
    y_padding_fraction: float,
) -> str:
    fig = plt.figure(figsize=(18.5, 10.2))
    gs = fig.add_gridspec(
        2,
        3,
        left=0.055,
        right=0.990,
        bottom=0.075,
        top=0.900,
        wspace=0.220,
        hspace=0.300,
    )

    shared_vals, shared_errs = collect_values_and_errors_for_canvas(
        values_df=values_df,
        samples=["data", "mc"],
        quantity=quantity,
    )

    shared_ylim = finite_min_max_with_errors(
        values=shared_vals,
        errors=shared_errs,
        pad_fraction=y_padding_fraction,
    )

    for ipanel, period in enumerate(PANEL_PERIOD_ORDER):
        ax = fig.add_subplot(gs[ipanel // 3, ipanel % 3])

        if period is None:
            ax.axis("off")
            continue

        x_base = np.arange(1, len(labels) + 1, dtype=float)

        title_lines = [period]

        for sample in SAMPLE_ORDER:
            group = get_group(values_df, period, sample, quantity)
            lookup = value_lookup_for_group(group)
            diag = get_diag_row(diagnostics_df, period, sample, quantity)

            style = SAMPLE_STYLE[sample]
            x = x_base + float(style["offset"])

            y = []
            ey = []

            for label in labels:
                value, stat = lookup.get(label, (math.nan, math.nan))
                y.append(value)
                ey.append(stat if np.isfinite(stat) and stat >= 0.0 else 0.0)

            ax.errorbar(
                x,
                y,
                yerr=ey,
                fmt=style["marker"],
                color=style["color"],
                label=style["label"],
                capsize=3,
                linewidth=1.2,
                markersize=6.2,
            )

            ref_mean = finite_float(diag.get("reference_mean", math.nan))
            ref_type = str(diag.get("reference_mean_type", "mean"))
            ref_stat = finite_float(diag.get("weighted_mean_stat", math.nan))

            if ref_type != "weighted":
                ref_stat = math.nan

            draw_mean_reference(
                ax,
                ref_mean,
                ref_stat,
                label=f"{style['label']} {ref_type}",
                color=style["color"],
                alpha=0.10,
            )

            chi2_ndf = finite_float(diag.get("chi2_ndf", math.nan))
            rms_percent = finite_float(diag.get("rms_over_mean_percent", math.nan))
            max_pull = finite_float(diag.get("max_pair_pull", math.nan))

            title_lines.append(
                f"{style['label']}: RMS/mean={rms_percent:.3g}%, "
                f"chi2/ndf={chi2_ndf:.3g}, max={max_pull:.3g}"
            )

        ax.set_title("\n".join(title_lines), fontsize=10)

        ax.set_xticks(x_base)
        ax.set_xticklabels(labels)

        ymin, ymax = get_panel_ylim(
            values_df=values_df,
            period=period,
            samples=["data", "mc"],
            quantity=quantity,
            shared_y_range=shared_y_range,
            shared_ylim=shared_ylim,
            y_padding_fraction=y_padding_fraction,
        )
        ax.set_ylim(ymin, ymax)

        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=10)

        if ipanel % 3 == 0:
            ax.set_ylabel(QUANTITY_LABEL.get(quantity, quantity), fontsize=11)

        if ipanel // 3 == 1:
            ax.set_xlabel("Sector/input label", fontsize=11)

        if ipanel == 0:
            ax.legend(fontsize=9, loc="best")

    y_mode = "shared y" if shared_y_range else "local y"
    fig.suptitle(
        f"DATA and MC {QUANTITY_LABEL.get(quantity, quantity)} by sector/input ({y_mode})",
        fontsize=16,
        y=0.975,
    )

    suffix = "shared_y" if shared_y_range else "local_y"

    outfile = os.path.join(
        output_dir,
        f"{output_prefix}_values_by_period_data_mc_{safe_filename_token(quantity)}_{suffix}.png",
    )

    fig.savefig(outfile, dpi=PLOT_DPI)
    plt.close(fig)

    return outfile


# -----------------------------------------------------------------------------
# Individual DATA-only or MC-only ratio plots.
# -----------------------------------------------------------------------------

def plot_ratio_to_mean_by_period(
    values_df: pd.DataFrame,
    diagnostics_df: pd.DataFrame,
    labels: List[str],
    sample: str,
    quantity: str,
    output_dir: str,
    output_prefix: str,
    shared_y_range: bool,
    y_padding_fraction: float,
) -> str:
    shared_ratios = []
    shared_ratio_errs = []

    for period in period_order_from_values(values_df):
        group = get_group(values_df, period, sample, quantity)
        diag = get_diag_row(diagnostics_df, period, sample, quantity)
        ref_mean = finite_float(diag.get("reference_mean", math.nan))

        if not np.isfinite(ref_mean) or ref_mean == 0.0:
            continue

        for _, row in group.iterrows():
            value = finite_float(row["value"])
            stat = finite_float(row["stat"], default=0.0)

            if np.isfinite(value):
                shared_ratios.append(value / ref_mean)
                shared_ratio_errs.append(stat / abs(ref_mean) if stat >= 0.0 else 0.0)

    shared_ylim = finite_min_max_with_errors(
        values=shared_ratios,
        errors=shared_ratio_errs,
        pad_fraction=y_padding_fraction,
    )
    shared_ylim = (
        min(shared_ylim[0], 0.995),
        max(shared_ylim[1], 1.005),
    )

    fig = plt.figure(figsize=(18.5, 10.2))
    gs = fig.add_gridspec(
        2,
        3,
        left=0.055,
        right=0.990,
        bottom=0.075,
        top=0.900,
        wspace=0.220,
        hspace=0.300,
    )

    for ipanel, period in enumerate(PANEL_PERIOD_ORDER):
        ax = fig.add_subplot(gs[ipanel // 3, ipanel % 3])

        if period is None:
            ax.axis("off")
            continue

        group = get_group(values_df, period, sample, quantity)
        lookup = value_lookup_for_group(group)
        diag = get_diag_row(diagnostics_df, period, sample, quantity)

        ref_mean = finite_float(diag.get("reference_mean", math.nan))
        ref_type = str(diag.get("reference_mean_type", "mean"))

        x = np.arange(1, len(labels) + 1, dtype=float)

        ratios = []
        ratio_errs = []

        for label in labels:
            value, stat = lookup.get(label, (math.nan, math.nan))

            if np.isfinite(value) and np.isfinite(ref_mean) and ref_mean != 0.0:
                ratios.append(value / ref_mean)

                if np.isfinite(stat) and stat >= 0.0:
                    ratio_errs.append(stat / abs(ref_mean))
                else:
                    ratio_errs.append(0.0)
            else:
                ratios.append(math.nan)
                ratio_errs.append(0.0)

        ax.errorbar(
            x,
            ratios,
            yerr=ratio_errs,
            fmt=SAMPLE_STYLE.get(sample, {}).get("marker", "o"),
            color=SAMPLE_STYLE.get(sample, {}).get("color", None),
            capsize=3,
            linewidth=1.2,
            markersize=6.5,
        )

        ax.axhline(1.0, linewidth=1.2, linestyle="--", color="black")

        chi2_ndf = finite_float(diag.get("chi2_ndf", math.nan))
        rms_percent = finite_float(diag.get("rms_over_mean_percent", math.nan))

        ax.set_title(
            f"{period}\n"
            f"ratio to {ref_type} mean; "
            f"RMS/mean={rms_percent:.3g}%, "
            f"chi2/ndf={chi2_ndf:.3g}",
            fontsize=11,
        )

        ax.set_xticks(x)
        ax.set_xticklabels(labels)

        if shared_y_range:
            ymin, ymax = shared_ylim
        else:
            ymin, ymax = finite_min_max_with_errors(
                values=ratios,
                errors=ratio_errs,
                pad_fraction=y_padding_fraction,
            )
            ymin = min(ymin, 0.995)
            ymax = max(ymax, 1.005)

        ax.set_ylim(ymin, ymax)

        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=10)

        if ipanel % 3 == 0:
            ax.set_ylabel("Value / sector mean", fontsize=11)

        if ipanel // 3 == 1:
            ax.set_xlabel("Sector/input label", fontsize=11)

    y_mode = "shared y" if shared_y_range else "local y"
    fig.suptitle(
        f"{sample.upper()} {QUANTITY_LABEL.get(quantity, quantity)} sector ratios ({y_mode})",
        fontsize=16,
        y=0.975,
    )

    suffix = "shared_y" if shared_y_range else "local_y"

    outfile = os.path.join(
        output_dir,
        f"{output_prefix}_ratio_to_mean_{sample}_{safe_filename_token(quantity)}_{suffix}.png",
    )

    fig.savefig(outfile, dpi=PLOT_DPI)
    plt.close(fig)

    return outfile


# -----------------------------------------------------------------------------
# Combined DATA+MC ratio plots.
# -----------------------------------------------------------------------------

def plot_data_mc_ratio_to_mean_by_period(
    values_df: pd.DataFrame,
    diagnostics_df: pd.DataFrame,
    labels: List[str],
    quantity: str,
    output_dir: str,
    output_prefix: str,
    shared_y_range: bool,
    y_padding_fraction: float,
) -> str:
    shared_ratios = []
    shared_ratio_errs = []

    for period in period_order_from_values(values_df):
        for sample in SAMPLE_ORDER:
            group = get_group(values_df, period, sample, quantity)
            diag = get_diag_row(diagnostics_df, period, sample, quantity)
            ref_mean = finite_float(diag.get("reference_mean", math.nan))

            if not np.isfinite(ref_mean) or ref_mean == 0.0:
                continue

            for _, row in group.iterrows():
                value = finite_float(row["value"])
                stat = finite_float(row["stat"], default=0.0)

                if np.isfinite(value):
                    shared_ratios.append(value / ref_mean)
                    shared_ratio_errs.append(stat / abs(ref_mean) if stat >= 0.0 else 0.0)

    shared_ylim = finite_min_max_with_errors(
        values=shared_ratios,
        errors=shared_ratio_errs,
        pad_fraction=y_padding_fraction,
    )
    shared_ylim = (
        min(shared_ylim[0], 0.995),
        max(shared_ylim[1], 1.005),
    )

    fig = plt.figure(figsize=(18.5, 10.2))
    gs = fig.add_gridspec(
        2,
        3,
        left=0.055,
        right=0.990,
        bottom=0.075,
        top=0.900,
        wspace=0.220,
        hspace=0.300,
    )

    for ipanel, period in enumerate(PANEL_PERIOD_ORDER):
        ax = fig.add_subplot(gs[ipanel // 3, ipanel % 3])

        if period is None:
            ax.axis("off")
            continue

        x_base = np.arange(1, len(labels) + 1, dtype=float)

        panel_ratios = []
        panel_ratio_errs = []

        title_lines = [period]

        for sample in SAMPLE_ORDER:
            style = SAMPLE_STYLE[sample]
            group = get_group(values_df, period, sample, quantity)
            lookup = value_lookup_for_group(group)
            diag = get_diag_row(diagnostics_df, period, sample, quantity)

            ref_mean = finite_float(diag.get("reference_mean", math.nan))
            ref_type = str(diag.get("reference_mean_type", "mean"))

            ratios = []
            ratio_errs = []

            for label in labels:
                value, stat = lookup.get(label, (math.nan, math.nan))

                if np.isfinite(value) and np.isfinite(ref_mean) and ref_mean != 0.0:
                    ratio = value / ref_mean
                    ratios.append(ratio)
                    panel_ratios.append(ratio)

                    if np.isfinite(stat) and stat >= 0.0:
                        ratio_err = stat / abs(ref_mean)
                    else:
                        ratio_err = 0.0

                    ratio_errs.append(ratio_err)
                    panel_ratio_errs.append(ratio_err)
                else:
                    ratios.append(math.nan)
                    ratio_errs.append(0.0)

            ax.errorbar(
                x_base + float(style["offset"]),
                ratios,
                yerr=ratio_errs,
                fmt=style["marker"],
                color=style["color"],
                label=f"{style['label']} / {ref_type} mean",
                capsize=3,
                linewidth=1.2,
                markersize=6.2,
            )

            rms_percent = finite_float(diag.get("rms_over_mean_percent", math.nan))
            chi2_ndf = finite_float(diag.get("chi2_ndf", math.nan))
            title_lines.append(
                f"{style['label']}: RMS/mean={rms_percent:.3g}%, "
                f"chi2/ndf={chi2_ndf:.3g}"
            )

        ax.axhline(1.0, linewidth=1.2, linestyle="--", color="black")

        ax.set_title("\n".join(title_lines), fontsize=10)

        ax.set_xticks(x_base)
        ax.set_xticklabels(labels)

        if shared_y_range:
            ymin, ymax = shared_ylim
        else:
            ymin, ymax = finite_min_max_with_errors(
                values=panel_ratios,
                errors=panel_ratio_errs,
                pad_fraction=y_padding_fraction,
            )
            ymin = min(ymin, 0.995)
            ymax = max(ymax, 1.005)

        ax.set_ylim(ymin, ymax)

        ax.grid(True, alpha=0.35)
        ax.tick_params(axis="both", labelsize=10)

        if ipanel % 3 == 0:
            ax.set_ylabel("Value / sector mean", fontsize=11)

        if ipanel // 3 == 1:
            ax.set_xlabel("Sector/input label", fontsize=11)

        if ipanel == 0:
            ax.legend(fontsize=9, loc="best")

    y_mode = "shared y" if shared_y_range else "local y"
    fig.suptitle(
        f"DATA and MC {QUANTITY_LABEL.get(quantity, quantity)} sector ratios ({y_mode})",
        fontsize=16,
        y=0.975,
    )

    suffix = "shared_y" if shared_y_range else "local_y"

    outfile = os.path.join(
        output_dir,
        f"{output_prefix}_ratio_to_mean_data_mc_{safe_filename_token(quantity)}_{suffix}.png",
    )

    fig.savefig(outfile, dpi=PLOT_DPI)
    plt.close(fig)

    return outfile


# -----------------------------------------------------------------------------
# Spread and heatmap plots.
# -----------------------------------------------------------------------------

def plot_spread_summary(
    diagnostics_df: pd.DataFrame,
    output_dir: str,
    output_prefix: str,
) -> str:
    period_order = [
        p for p in PERIOD_ORDER
        if p in set(diagnostics_df["period"].astype(str))
    ]

    extra_periods = sorted(set(diagnostics_df["period"].astype(str)) - set(period_order))
    period_order = period_order + extra_periods

    fig = plt.figure(figsize=(18.5, 10.0))
    gs = fig.add_gridspec(
        2,
        2,
        left=0.060,
        right=0.985,
        bottom=0.080,
        top=0.900,
        wspace=0.220,
        hspace=0.300,
    )

    combos = [
        ("data", "current_efficiency_factor"),
        ("mc", "current_efficiency_factor"),
        ("data", "slope_percent_per_nA"),
        ("mc", "slope_percent_per_nA"),
    ]

    x = np.arange(len(period_order), dtype=float)
    width = 0.36

    for i, (sample, quantity) in enumerate(combos):
        ax = fig.add_subplot(gs[i // 2, i % 2])

        subset = diagnostics_df[
            (diagnostics_df["sample"] == sample)
            & (diagnostics_df["quantity"] == quantity)
        ].copy()

        lookup = {
            str(row["period"]): row
            for _, row in subset.iterrows()
        }

        rms = []
        chi2ndf = []

        for period in period_order:
            row = lookup.get(period, None)

            if row is None:
                rms.append(math.nan)
                chi2ndf.append(math.nan)
            else:
                rms.append(float(row["rms_over_mean_percent"]))
                chi2ndf.append(float(row["chi2_ndf"]))

        ax.bar(x - 0.5 * width, rms, width=width, label="RMS/mean (%)")
        ax.bar(x + 0.5 * width, chi2ndf, width=width, label="chi2/ndf")

        ax.axhline(1.0, linewidth=1.0, linestyle="--", color="black", alpha=0.8)

        ax.set_title(
            f"{sample.upper()} {QUANTITY_LABEL.get(quantity, quantity)}",
            fontsize=12,
        )
        ax.set_xticks(x)
        ax.set_xticklabels(period_order, rotation=25, ha="right")
        ax.set_ylabel("Diagnostic value")
        ax.grid(True, axis="y", alpha=0.35)
        ax.legend(fontsize=9)

    fig.suptitle(
        "Current-dependence sector spread summary",
        fontsize=16,
        y=0.975,
    )

    outfile = os.path.join(
        output_dir,
        f"{output_prefix}_spread_summary.png",
    )

    fig.savefig(outfile, dpi=PLOT_DPI)
    plt.close(fig)

    return outfile


def plot_sector_pull_heatmap(
    values_df: pd.DataFrame,
    diagnostics_df: pd.DataFrame,
    labels: List[str],
    sample: str,
    quantity: str,
    output_dir: str,
    output_prefix: str,
) -> str:
    periods = period_order_from_values(values_df)

    matrix = np.full((len(periods), len(labels)), np.nan)

    for ip, period in enumerate(periods):
        group = get_group(values_df, period, sample, quantity)
        lookup = value_lookup_for_group(group)
        diag = get_diag_row(diagnostics_df, period, sample, quantity)

        ref_mean = finite_float(diag.get("reference_mean", math.nan))

        for il, label in enumerate(labels):
            value, stat = lookup.get(label, (math.nan, math.nan))

            if (
                np.isfinite(value)
                and np.isfinite(stat)
                and stat > 0.0
                and np.isfinite(ref_mean)
            ):
                matrix[ip, il] = (value - ref_mean) / stat

    finite = matrix[np.isfinite(matrix)]

    if finite.size > 0:
        vmax = float(np.nanmax(np.abs(finite)))
        vmax = max(vmax, 1.0)
    else:
        vmax = 1.0

    fig, ax = plt.subplots(figsize=(12.0, 6.8))

    im = ax.imshow(
        matrix,
        aspect="auto",
        vmin=-vmax,
        vmax=vmax,
        cmap="coolwarm",
    )

    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yticks(np.arange(len(periods)))
    ax.set_yticklabels(periods)

    ax.set_xlabel("Sector/input label")
    ax.set_ylabel("Run period")

    ax.set_title(
        f"{sample.upper()} {QUANTITY_LABEL.get(quantity, quantity)} signed pull to sector mean"
    )

    for ip in range(len(periods)):
        for il in range(len(labels)):
            val = matrix[ip, il]
            if np.isfinite(val):
                ax.text(
                    il,
                    ip,
                    f"{val:.1f}",
                    ha="center",
                    va="center",
                    fontsize=9,
                )
            else:
                ax.text(
                    il,
                    ip,
                    "nan",
                    ha="center",
                    va="center",
                    fontsize=8,
                )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("(value - mean) / stat")

    fig.tight_layout()

    outfile = os.path.join(
        output_dir,
        f"{output_prefix}_sector_pull_heatmap_{sample}_{safe_filename_token(quantity)}.png",
    )

    fig.savefig(outfile, dpi=PLOT_DPI)
    plt.close(fig)

    return outfile


# -----------------------------------------------------------------------------
# Plot driver.
# -----------------------------------------------------------------------------

def make_all_plots(
    values_df: pd.DataFrame,
    diagnostics_df: pd.DataFrame,
    labels: List[str],
    output_dir: str,
    output_prefix: str,
    shared_y_range: bool,
    y_padding_fraction: float,
) -> List[str]:
    os.makedirs(output_dir, exist_ok=True)

    outputs = []

    with Timer("making current-dependence plots"):
        for quantity in QUANTITY_ORDER:
            for sample in SAMPLE_ORDER:
                outputs.append(
                    plot_values_by_period(
                        values_df=values_df,
                        diagnostics_df=diagnostics_df,
                        labels=labels,
                        sample=sample,
                        quantity=quantity,
                        output_dir=output_dir,
                        output_prefix=output_prefix,
                        shared_y_range=shared_y_range,
                        y_padding_fraction=y_padding_fraction,
                    )
                )

                outputs.append(
                    plot_ratio_to_mean_by_period(
                        values_df=values_df,
                        diagnostics_df=diagnostics_df,
                        labels=labels,
                        sample=sample,
                        quantity=quantity,
                        output_dir=output_dir,
                        output_prefix=output_prefix,
                        shared_y_range=shared_y_range,
                        y_padding_fraction=y_padding_fraction,
                    )
                )

                outputs.append(
                    plot_sector_pull_heatmap(
                        values_df=values_df,
                        diagnostics_df=diagnostics_df,
                        labels=labels,
                        sample=sample,
                        quantity=quantity,
                        output_dir=output_dir,
                        output_prefix=output_prefix,
                    )
                )

            outputs.append(
                plot_data_mc_values_by_period(
                    values_df=values_df,
                    diagnostics_df=diagnostics_df,
                    labels=labels,
                    quantity=quantity,
                    output_dir=output_dir,
                    output_prefix=output_prefix,
                    shared_y_range=shared_y_range,
                    y_padding_fraction=y_padding_fraction,
                )
            )

            outputs.append(
                plot_data_mc_ratio_to_mean_by_period(
                    values_df=values_df,
                    diagnostics_df=diagnostics_df,
                    labels=labels,
                    quantity=quantity,
                    output_dir=output_dir,
                    output_prefix=output_prefix,
                    shared_y_range=shared_y_range,
                    y_padding_fraction=y_padding_fraction,
                )
            )

        outputs.append(
            plot_spread_summary(
                diagnostics_df=diagnostics_df,
                output_dir=output_dir,
                output_prefix=output_prefix,
            )
        )

    return outputs


# -----------------------------------------------------------------------------
# CLI.
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze sector/version dependence of current-dependence summary CSVs "
            "from current_dependence.cpp."
        )
    )

    parser.add_argument(
        "csv_files",
        nargs="+",
        help=(
            "Current-dependence period_summary.csv files to compare. "
            "For electron/photon/proton FD-sector studies, pass six files."
        ),
    )

    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help=(
            "Labels for the input CSVs. Must have the same number of entries "
            "as csv_files. Default: S1 S2 ... SN."
        ),
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
        "--shared-y-range",
        action="store_true",
        help=(
            "Use one common y-axis range across all panels in each value/ratio canvas. "
            "Default is local per-panel y-axis zoom."
        ),
    )

    parser.add_argument(
        "--y-padding-fraction",
        type=float,
        default=0.30,
        help=(
            "Fractional padding around the finite value/error range for zoomed y axes. "
            "Default: 0.30."
        ),
    )

    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip PNG plot creation.",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> List[str]:
    if len(args.csv_files) < 2:
        raise ValueError("Need at least two current-dependence CSV files to compare.")

    if args.labels is None:
        labels = [f"S{i}" for i in range(1, len(args.csv_files) + 1)]
    else:
        labels = list(args.labels)

    if len(labels) != len(args.csv_files):
        raise ValueError(
            f"Number of labels ({len(labels)}) must match number of CSV files "
            f"({len(args.csv_files)})."
        )

    if args.y_padding_fraction < 0.0:
        raise ValueError("--y-padding-fraction must be non-negative.")

    return labels


def main() -> None:
    t0 = time.perf_counter()
    args = parse_args()
    labels = validate_args(args)

    log("============================================================")
    log("analyze_current_dependence_csvs.py")
    log("============================================================")
    log(f"CSV files: {args.csv_files}")
    log(f"Labels: {labels}")
    log(f"Output dir: {args.output_dir}")
    log(f"Output prefix: {args.output_prefix}")
    log(f"Shared y range: {args.shared_y_range}")
    log(f"Y padding fraction: {args.y_padding_fraction}")

    os.makedirs(args.output_dir, exist_ok=True)

    entries = read_all_summary_csvs(
        paths=args.csv_files,
        labels=labels,
    )

    values_df = build_values_dataframe(entries)
    diagnostics_df = build_diagnostics_dataframe(values_df)
    wide_values_df = build_wide_values_dataframe(values_df, labels)

    values_path = os.path.join(
        args.output_dir,
        f"{args.output_prefix}_values_long.csv",
    )

    wide_values_path = os.path.join(
        args.output_dir,
        f"{args.output_prefix}_values_wide.csv",
    )

    diagnostics_path = os.path.join(
        args.output_dir,
        f"{args.output_prefix}_diagnostics.csv",
    )

    report_path = os.path.join(
        args.output_dir,
        f"{args.output_prefix}_report.txt",
    )

    values_df.to_csv(values_path, index=False)
    wide_values_df.to_csv(wide_values_path, index=False)
    diagnostics_df.to_csv(diagnostics_path, index=False)

    write_text_report(
        path=report_path,
        values_df=values_df,
        diagnostics_df=diagnostics_df,
        labels=labels,
    )

    plot_outputs = []

    if not args.no_plots:
        plot_outputs = make_all_plots(
            values_df=values_df,
            diagnostics_df=diagnostics_df,
            labels=labels,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
            shared_y_range=args.shared_y_range,
            y_padding_fraction=args.y_padding_fraction,
        )

    print_console_summary(diagnostics_df)

    log(f"Wrote long values CSV:  {values_path}")
    log(f"Wrote wide values CSV:  {wide_values_path}")
    log(f"Wrote diagnostics CSV:  {diagnostics_path}")
    log(f"Wrote text report:      {report_path}")

    if plot_outputs:
        log("Wrote plots:")
        for path in plot_outputs:
            log(f"  {path}")

    log(f"TOTAL RUNTIME: {time.perf_counter() - t0:.3f} s")
    log("Done.")


if __name__ == "__main__":
    main()