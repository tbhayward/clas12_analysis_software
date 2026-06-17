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

  python analyze_current_dependence_csvs.py \
    sec1/period_summary.csv \
    sec2/period_summary.csv \
    sec3/period_summary.csv \
    sec4/period_summary.csv \
    sec5/period_summary.csv \
    sec6/period_summary.csv \
    --labels S1 S2 S3 S4 S5 S6

Usage for three proton CD sector summaries:

  python analyze_current_dependence_csvs.py \
    cd_sec1/period_summary.csv \
    cd_sec2/period_summary.csv \
    cd_sec3/period_summary.csv \
    --labels "CD S1" "CD S2" "CD S3"

Outputs by default:

  output/current_dependence_sector_diagnostics/current_dependence_values.csv
  output/current_dependence_sector_diagnostics/current_dependence_diagnostics.csv
  output/current_dependence_sector_diagnostics/current_dependence_report.txt

The long values CSV contains one row per:

  input label, period, sample type, quantity

The diagnostics CSV contains one row per:

  period, sample type, quantity

with:

  unweighted mean
  standard deviation
  RMS spread about the mean
  RMS/mean percent
  weighted mean
  weighted-mean uncertainty
  chi2
  ndf
  chi2/ndf
  p-value if scipy is available
  largest pairwise pull
  min/max labels and values

Statistical-consistency convention:

  weighted mean:
    xbar = sum_i x_i / sigma_i^2 / sum_i 1 / sigma_i^2

  chi2:
    chi2 = sum_i (x_i - xbar)^2 / sigma_i^2
    ndf  = N - 1

This assumes the input uncertainties are independent. That is probably not
strictly true for sector splits from the same data set, so chi2/ndf should be
treated as a diagnostic tension metric, not a formal hypothesis test.

If some uncertainties are missing or non-positive, the script still reports
unweighted means/RMS spreads. Weighted chi2 diagnostics are computed only from
points with finite positive uncertainty.
"""

import argparse
import math
import os
import re
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


PERIOD_ORDER = [
    "Sp18 Inb",
    "Sp18 Out",
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


def require_columns(df: pd.DataFrame, columns: Iterable[str], path: str) -> None:
    missing = [c for c in columns if c not in df.columns]

    if missing:
        print()
        print(f"Missing required columns in: {path}")
        for c in missing:
            print(f"  {c}")
        # endfor

        print()
        print("Available columns:")
        for c in df.columns:
            print(f"  {c}")
        # endfor

        raise KeyError(f"{path} is missing {len(missing)} required columns.")
    # endif


# -----------------------------------------------------------------------------
# Current-dependence quantity helpers.
# -----------------------------------------------------------------------------

def percent_slope(m: float, b: float) -> float:
    if not np.isfinite(m) or not np.isfinite(b) or b == 0.0:
        return math.nan
    # endif

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
    # endif

    dm = 100.0 / b
    db = -100.0 * m / (b * b)

    var = (
        dm * dm * sm * sm
        + db * db * sb * sb
        + 2.0 * dm * db * cov_mb
    )

    if var < 0.0 and abs(var) < 1.0e-15:
        var = 0.0
    # endif

    if var < 0.0:
        return math.nan
    # endif

    return math.sqrt(var)


def get_scalar(row: pd.Series, column: str) -> float:
    if column not in row.index:
        return math.nan
    # endif

    return parse_first_number(row[column])


def get_data_slope_value_and_error(row: pd.Series) -> Tuple[float, float]:
    value = get_scalar(row, "data_slope_percent_per_nA")
    stat = get_scalar(row, "data_slope_percent_per_nA_stat")

    if np.isfinite(value):
        return value, stat
    # endif

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
        # endif
    # endif

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
    # endwith

    require_columns(df, required, path)

    log(f"{label}: {df.shape[0]} rows x {df.shape[1]} columns")

    entries: List[ValueEntry] = []

    for _, row in df.iterrows():
        period = str(row["period"]).strip()

        if period == "" or period.lower() == "nan":
            continue
        # endif

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
    # endfor

    return entries


def read_all_summary_csvs(
    paths: List[str],
    labels: List[str],
) -> List[ValueEntry]:
    all_entries: List[ValueEntry] = []

    for path, label in zip(paths, labels):
        entries = read_one_summary_csv(path=path, label=label)
        all_entries.extend(entries)
    # endfor

    return all_entries


# -----------------------------------------------------------------------------
# Statistical diagnostics.
# -----------------------------------------------------------------------------

def chi2_pvalue(chi2: float, ndf: int) -> float:
    if not np.isfinite(chi2) or ndf <= 0:
        return math.nan
    # endif

    try:
        from scipy.stats import chi2 as scipy_chi2
        return float(scipy_chi2.sf(chi2, ndf))
    except Exception:
        return math.nan
    # endtry


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
    # endif

    w = 1.0 / (s * s)
    wsum = float(np.sum(w))

    if not np.isfinite(wsum) or wsum <= 0.0:
        return math.nan, math.nan, math.nan, 0, math.nan, 0
    # endif

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
            # endif

            denom = math.sqrt(si * si + sj * sj)

            if denom <= 0.0:
                continue
            # endif

            pull = abs(xi - xj) / denom

            if not np.isfinite(best_pull) or pull > best_pull:
                best_pull = float(pull)
                best_a = labels[i]
                best_b = labels[j]
            # endif
        # endfor
    # endfor

    return best_pull, best_a, best_b


def diagnostics_for_group(group: pd.DataFrame) -> Dict[str, object]:
    if len(group) == 0:
        return empty_diagnostics()
    # endif

    finite = group[np.isfinite(group["value"].astype(float))].copy()

    if len(finite) == 0:
        return empty_diagnostics()
    # endif

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
    # endif

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
    # endif

    weighted_mean, weighted_mean_stat, chi2, ndf, p_value, n_weighted = weighted_mean_and_chi2(
        values=values,
        stats=stats,
    )

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
        "chi2": chi2,
        "ndf": int(ndf),
        "chi2_ndf": chi2_ndf,
        "p_value": p_value,
        "max_pair_pull": max_pull,
        "max_pair_pull_a": max_pull_a,
        "max_pair_pull_b": max_pull_b,
    }


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
        "chi2": math.nan,
        "ndf": 0,
        "chi2_ndf": math.nan,
        "p_value": math.nan,
        "max_pair_pull": math.nan,
        "max_pair_pull_a": "",
        "max_pair_pull_b": "",
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
    # endfor

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
            # endfor
        # endfor
    # endfor

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Reporting.
# -----------------------------------------------------------------------------

def format_value_err(value: float, err: float) -> str:
    if not np.isfinite(value):
        return "nan"
    # endif

    if np.isfinite(err) and err >= 0.0:
        return f"{value:.8g} +/- {err:.3g}"
    # endif

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
    # endfor


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
    # endif

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
        # endfor

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
                # endfor
            # endfor
        # endfor


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
    # endfor


# -----------------------------------------------------------------------------
# Optional wide summary table.
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
                # endfor

                rows.append(row)
            # endfor
        # endfor
    # endfor

    return pd.DataFrame(rows)


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

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> List[str]:
    if len(args.csv_files) < 2:
        raise ValueError("Need at least two current-dependence CSV files to compare.")
    # endif

    if args.labels is None:
        labels = [f"S{i}" for i in range(1, len(args.csv_files) + 1)]
    else:
        labels = list(args.labels)
    # endif

    if len(labels) != len(args.csv_files):
        raise ValueError(
            f"Number of labels ({len(labels)}) must match number of CSV files "
            f"({len(args.csv_files)})."
        )
    # endif

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

    print_console_summary(diagnostics_df)

    log(f"Wrote long values CSV:  {values_path}")
    log(f"Wrote wide values CSV:  {wide_values_path}")
    log(f"Wrote diagnostics CSV:  {diagnostics_path}")
    log(f"Wrote text report:      {report_path}")
    log(f"TOTAL RUNTIME: {time.perf_counter() - t0:.3f} s")
    log("Done.")


if __name__ == "__main__":
    main()
# endif