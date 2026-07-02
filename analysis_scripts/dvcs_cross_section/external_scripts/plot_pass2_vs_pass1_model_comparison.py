#!/usr/bin/env python3
"""
plot_pass2_vs_pass1_model_comparison.py
--------------------------------------------------------------------------
Standalone comparison plotter for pass-2 vs pass-1 DVCS cross sections
against BH and KM15 model predictions.

Usage:
    python3 plot_pass2_vs_pass1_model_comparison.py \
        pass2_dvcs.csv pass1_lee.csv \
        --output-dir output/pass2_vs_pass1_model_comparison

Optional two-panel mode:
    python3 plot_pass2_vs_pass1_model_comparison.py \
        pass2_dvcs.csv pass1_lee.csv \
        --output-dir output/pass2_vs_pass1_model_comparison \
        --two-panel

Optional removal of point-to-point systematics from plotted vertical bars:
    python3 plot_pass2_vs_pass1_model_comparison.py \
        pass2_dvcs.csv pass1_lee.csv \
        --output-dir output/pass2_vs_pass1_model_comparison \
        --two-panel \
        --no-point-to-point-systematics

Implemented uncertainty prescription:
  * Positional input order is pass-2 CSV first, pass-1 CSV second.

  * Pass-2 central values are read from the 10.6 GeV unpolarized normalized
    cross-section column:

        normed cross sections, ep->epg, exp, 10.6 GeV, unpol

    The script intentionally does NOT use:

        normed cross sections, ep->epg, exp, 10.6 GeV, unpol, combination sys

    as a cross-section input. That column is treated as an analysis-chain
    output/fill column, not as a central-value source for this comparison.

  * Pass-2 scale systematic is computed point-by-point from the four 10.6 GeV
    unpolarized period columns in the same row:

        Fa18 Inb, Fa18 Out, Sp18 Inb, Sp18 Out

    For each row:
      1. Form a stat-weighted reference from the valid period tuples.
      2. Form period/reference ratios for all valid period tuples.
      3. Compute

             s_obs(row)  = RMS_i(r_i - 1)
             s_stat(row) = RMS_i(stat_i / ref)
             s_comb(row) = sqrt(max(0, s_obs(row)^2 - s_stat(row)^2))

      4. The pass-2 scale box for that point is:

             s_comb(row) * sigma_pass2(row)

    A row must have at least two valid period entries to get a nonzero local
    scale systematic. Rows with fewer than two valid period entries get local
    s_comb = 0.

  * Cross-section panel, default:
      pass-2 vertical error bars:
          stat ⊕ 18% estimated point-to-point/systematic

      pass-2 external scale boxes:
          local point-by-point s_comb(row) * cross_section

      pass-1 vertical error bars:
          stat ⊕ provided Lee systematic

      pass-1 external scale boxes:
          31% normalization * cross_section

  * Cross-section panel with --no-point-to-point-systematics:
      pass-2 vertical error bars:
          stat only

      pass-2 external scale boxes:
          local point-by-point s_comb(row) * cross_section

      pass-1 vertical error bars:
          stat only

      pass-1 external scale boxes:
          31% normalization * cross_section

  * Ratio panel:
      pass-2/pass-1 vertical error bars:
          propagated statistical uncertainty only

      pass-2/pass-1 external scale boxes:
          propagated scale uncertainty from local pass-2 s_comb(row)
          and pass-1 31% normalization

      The pass-2 18% estimated point-to-point systematic and the pass-1 provided
      systematic columns are intentionally not included in the ratio-panel error
      bars or ratio-panel scale boxes.

  * Additional standalone stat-only pull diagnostic:
      For every pass-2/pass-1 matched point, compute

          pull(M) = (pass2 - M * pass1) /
                    sqrt(pass2_stat^2 + (M * pass1_stat)^2)

      using statistical uncertainties only. The script writes a histogram of
      the pull distribution before and after fitting a global pass-1
      normalization M. M is bounded to [0.69, 1.31], i.e. +/-31%.

      This diagnostic deliberately does not use:
          - pass-2 estimated 18% point-to-point systematic,
          - pass-1 provided systematic,
          - pass-2 local scale systematic,
          - pass-1 31% normalization uncertainty.

  * One PNG is saved for every unique (xB, Q2, |t|) bin found in either CSV.
  * The model curves are BH and KM15 only, evaluated at the bin midpoint and
    drawn as functions of phi.
  * Pass-1 and pass-2 markers are shifted slightly in phi on the cross-section
    panel so overlapping points/error bars remain visually distinguishable.
  * By default, log-y plots use a lower visible floor of 1e-5. Change with
    --log-y-min.

Parallelization:
  * Each unique (xB, Q2, |t|) panel is processed independently.
  * By default, up to 5 worker processes are used.
  * Use --workers N to change this.
  * Values above 5 are capped to 5.
  * Use --workers 1 for fully serial debugging.

External model tools:
  * BH is evaluated by running dvcsgen.
  * KM15 is evaluated by running km15_cli.py.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# ---------------------------------------------------------------------------
# Default external paths.
# ---------------------------------------------------------------------------
SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
ANALYSIS_DIR = SCRIPT_DIR.parent

DEFAULT_DVCSGEN_DIR = "/u/home/thayward/dvcsgens/dvcsgen_print"
DEFAULT_KM15_CLI = str(ANALYSIS_DIR / "km15_cli.py")


# ---------------------------------------------------------------------------
# Uncertainty prescriptions.
# ---------------------------------------------------------------------------
PASS2_ESTIMATED_SYSTEMATIC_FRACTION = 0.18
PASS1_NORMALIZATION_FRACTION = 0.31


# ---------------------------------------------------------------------------
# Plotting prescription.
# ---------------------------------------------------------------------------
PASS1_PHI_OFFSET_DEG = -1.25
PASS2_PHI_OFFSET_DEG = +1.25

PASS1_PASS2_RATIO_MATCH_TOLERANCE_DEG = 2.0

SCALE_BOX_HALF_WIDTH_DEG = 3.0
SCALE_BOX_ALPHA = 0.18

PULL_FIT_NORM_MIN = 1.0 - PASS1_NORMALIZATION_FRACTION
PULL_FIT_NORM_MAX = 1.0 + PASS1_NORMALIZATION_FRACTION


# ---------------------------------------------------------------------------
# Pass-1 / Lee CSV column names.
# ---------------------------------------------------------------------------
PASS1_XS_COL = "cross sections, ep->epg, exp"
PASS1_STAT_COL = "cross sections, ep->epg, exp, stat. unc."
PASS1_SYST_UP_COL = "cross sections, ep->epg, exp, syst. unc. (up)"
PASS1_SYST_DN_COL = "cross sections, ep->epg, exp, syst. unc. (down)"


# ---------------------------------------------------------------------------
# Pass-2 central-value column candidates.
# The combination-sys column is deliberately excluded from this list.
# ---------------------------------------------------------------------------
PASS2_CENTRAL_XS_CANDIDATES = [
    "normed cross sections, ep->epg, exp, 10.6 GeV, unpol",
    "cross sections, ep->epg, exp, 10.6 GeV, unpol",
]


PASS2_PERIODS_10P6_UNPOL = [
    "Fa18 Inb",
    "Fa18 Out",
    "Sp18 Inb",
    "Sp18 Out",
]


def pass2_period_cross_section_column(period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {period}, unpol"


# ---------------------------------------------------------------------------
# Data containers.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class BinKey:
    xb_min: float
    xb_max: float
    q2_min: float
    q2_max: float
    t_min: float
    t_max: float


@dataclass
class TupleValue:
    ok: bool = False
    value: float = 0.0
    stat: float = 0.0
    sys: float = 0.0


@dataclass
class LocalScaleResult:
    n_valid_periods: int = 0
    ref_value: float = 0.0
    ref_stat: float = 0.0
    s_obs: float = 0.0
    s_stat_exp: float = 0.0
    s_comb: float = 0.0
    min_ratio: float = 0.0
    max_ratio: float = 0.0
    half_width: float = 0.0


@dataclass
class PointScaleSummary:
    n_rows: int = 0
    n_rows_with_central_value: int = 0
    n_rows_with_local_scale: int = 0
    n_rows_with_less_than_two_periods: int = 0
    mean_s_comb: float = 0.0
    rms_s_comb: float = 0.0
    median_s_comb: float = 0.0
    max_s_comb: float = 0.0


@dataclass
class DataPoint:
    phi: float
    xs: float
    stat: float
    err_low: float
    err_high: float
    source: str
    scale_syst_abs: float = 0.0
    scale_syst_frac: float = 0.0
    local_scale_n_periods: int = 0
    local_scale_s_obs: float = 0.0
    local_scale_s_stat: float = 0.0
    local_scale_half_width: float = 0.0


@dataclass
class PanelData:
    key: BinKey
    pass2: List[DataPoint]
    pass1: List[DataPoint]


@dataclass
class ModelCurves:
    phi: List[float]
    bh: List[float]
    km15: List[float]


@dataclass
class RatioPoint:
    phi: float
    ratio: float
    stat_err_low: float
    stat_err_high: float
    scale_err: float
    pass2_scale_frac: float
    label: str


@dataclass
class StatOnlyComparisonPoint:
    key: BinKey
    phi_pass2: float
    phi_pass1: float
    pass2_xs: float
    pass2_stat: float
    pass1_xs: float
    pass1_stat: float


@dataclass
class StatOnlyPullSummary:
    n_points: int = 0
    best_pass1_norm: float = 1.0
    best_norm_min: float = PULL_FIT_NORM_MIN
    best_norm_max: float = PULL_FIT_NORM_MAX
    chi2_nominal: float = 0.0
    chi2_best: float = 0.0
    ndf_nominal: int = 0
    ndf_best: int = 0
    chi2_ndf_nominal: float = 0.0
    chi2_ndf_best: float = 0.0
    pct_within_1sigma_nominal: float = 0.0
    pct_within_3sigma_nominal: float = 0.0
    pct_within_1sigma_best: float = 0.0
    pct_within_3sigma_best: float = 0.0
    mean_pull_nominal: float = 0.0
    rms_pull_nominal: float = 0.0
    mean_pull_best: float = 0.0
    rms_pull_best: float = 0.0


@dataclass
class ModelConfig:
    e_beam: float
    dvcsgen_dir: str
    km15_cli: str
    py_km15: str
    phi_dense: int
    use_cache: bool
    allow_missing_models: bool
    verbose_worker_models: bool


@dataclass
class WorkerJob:
    index: int
    total: int
    panel: PanelData
    output_path: str
    filename: str
    pass2_label: str
    pass1_label: str
    logy: bool
    log_y_min: float
    two_panel: bool
    skip_models: bool
    model_cfg: ModelConfig
    cached_model_entry: Optional[Dict[str, List[float]]]


@dataclass
class WorkerResult:
    ok: bool
    index: int
    filename: str
    output_path: str
    manifest_entry: Dict[str, object]
    cache_key: Optional[str]
    cache_entry: Optional[Dict[str, List[float]]]
    message: str
    model_status: str
    elapsed_seconds: float
    error: str


# ---------------------------------------------------------------------------
# Status and generic helpers.
# ---------------------------------------------------------------------------
def timestamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[pass2-vs-pass1][{timestamp()}] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[pass2-vs-pass1][{timestamp()}][warn] {msg}", file=sys.stderr, flush=True)


def die(msg: str) -> None:
    raise RuntimeError(msg)


def format_seconds(seconds: float) -> str:
    if seconds < 60.0:
        return f"{seconds:.1f} s"
    #endif

    minutes = int(seconds // 60.0)
    rem = seconds - 60.0 * minutes

    if minutes < 60:
        return f"{minutes:d} min {rem:.1f} s"
    #endif

    hours = minutes // 60
    minutes = minutes % 60

    return f"{hours:d} h {minutes:d} min {rem:.1f} s"


def format_fraction_percent(frac: float) -> str:
    if not math.isfinite(frac):
        return "nan%"
    #endif

    return f"{100.0 * frac:.2f}%"


def canonical_header_map(fieldnames: Sequence[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}

    for name in fieldnames:
        out[name.strip().lower()] = name
    #endfor

    return out


def find_column(fieldnames: Sequence[str], candidates: Sequence[str], required: bool = True) -> Optional[str]:
    cmap = canonical_header_map(fieldnames)

    for cand in candidates:
        key = cand.strip().lower()

        if key in cmap:
            return cmap[key]
        #endif
    #endfor

    if required:
        formatted = "\n".join(f"  - {c}" for c in candidates)
        die(f"Could not find required column. Tried:\n{formatted}")
    #endif

    return None


def require_columns(fieldnames: Sequence[str], required: Sequence[str], context: str) -> None:
    missing: List[str] = []
    cmap = canonical_header_map(fieldnames)

    for col in required:
        if col.strip().lower() not in cmap:
            missing.append(col)
        #endif
    #endfor

    if missing:
        formatted = "\n".join(f"  - {m}" for m in missing)
        die(f"Missing required columns for {context}:\n{formatted}")
    #endif


def to_float(value: object, default: float = 0.0) -> float:
    if value is None:
        return default
    #endif

    text = str(value).strip()

    if text == "":
        return default
    #endif

    try:
        out = float(text)

        if math.isfinite(out):
            return out
        #endif

        return default
    except ValueError:
        match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)

        if match:
            try:
                out = float(match.group(0))

                if math.isfinite(out):
                    return out
                #endif

                return default
            except ValueError:
                return default
            #endtry
        #endif
    #endtry

    return default


def to_int(value: object, default: int = 0) -> int:
    if value is None:
        return default
    #endif

    text = str(value).strip()

    if text == "":
        return default
    #endif

    try:
        return int(float(text))
    except ValueError:
        return default
    #endtry


def finite_positive(value: float) -> bool:
    return math.isfinite(value) and value > 0.0


def round_key_value(value: float, ndigits: int = 6) -> float:
    return round(float(value), ndigits)


def make_bin_key(row: Dict[str, str], cols: Dict[str, str]) -> BinKey:
    return BinKey(
        xb_min=round_key_value(to_float(row[cols["xb_min"]])),
        xb_max=round_key_value(to_float(row[cols["xb_max"]])),
        q2_min=round_key_value(to_float(row[cols["q2_min"]])),
        q2_max=round_key_value(to_float(row[cols["q2_max"]])),
        t_min=round_key_value(to_float(row[cols["t_min"]])),
        t_max=round_key_value(to_float(row[cols["t_max"]])),
    )


def key_sort_tuple(key: BinKey) -> Tuple[float, float, float, float, float, float]:
    return (key.xb_min, key.xb_max, key.q2_min, key.q2_max, key.t_min, key.t_max)


def safe_filename_piece(value: float) -> str:
    text = f"{value:.6g}"
    text = text.replace("-", "m")
    text = text.replace("+", "p")
    text = text.replace(".", "p")

    return text


def panel_filename(key: BinKey, index: int) -> str:
    return (
        f"panel_{index:03d}_"
        f"xB_{safe_filename_piece(key.xb_min)}_{safe_filename_piece(key.xb_max)}_"
        f"Q2_{safe_filename_piece(key.q2_min)}_{safe_filename_piece(key.q2_max)}_"
        f"t_{safe_filename_piece(key.t_min)}_{safe_filename_piece(key.t_max)}.png"
    )


def csv_field_escape(s: str) -> str:
    text = str(s)

    need_quote = "," in text or '"' in text or "\n" in text or "\r" in text

    if not need_quote:
        return text
    #endif

    return '"' + text.replace('"', '""') + '"'


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return ""
    #endif

    return f"{value:.12g}"


def mean(values: Sequence[float]) -> float:
    vals = [v for v in values if math.isfinite(v)]

    if not vals:
        return 0.0
    #endif

    return sum(vals) / float(len(vals))


def rms(values: Sequence[float]) -> float:
    vals = [v for v in values if math.isfinite(v)]

    if not vals:
        return 0.0
    #endif

    return math.sqrt(sum(v * v for v in vals) / float(len(vals)))


def median(values: Sequence[float]) -> float:
    vals = sorted(v for v in values if math.isfinite(v))

    if not vals:
        return 0.0
    #endif

    n = len(vals)
    mid = n // 2

    if n % 2 == 1:
        return vals[mid]
    #endif

    return 0.5 * (vals[mid - 1] + vals[mid])


# ---------------------------------------------------------------------------
# CSV column detection.
# ---------------------------------------------------------------------------
def detect_common_columns(fieldnames: Sequence[str], label: str) -> Dict[str, str]:
    log(f"{label}: detecting common kinematic/valid-bin columns.")

    cols: Dict[str, str] = {}

    cols["valid"] = find_column(fieldnames, ["valid bin", "valid", "is valid"], required=False) or ""
    cols["bin_index"] = find_column(fieldnames, ["bin index", "bin", "idx"], required=False) or ""

    cols["xb_min"] = find_column(fieldnames, ["xBmin", "xbmin", "xB_min", "xb_min"])
    cols["xb_max"] = find_column(fieldnames, ["xBmax", "xbmax", "xB_max", "xb_max"])
    cols["q2_min"] = find_column(fieldnames, ["Q2min", "q2min", "Q2_min", "q2_min"])
    cols["q2_max"] = find_column(fieldnames, ["Q2max", "q2max", "Q2_max", "q2_max"])
    cols["t_min"] = find_column(fieldnames, ["t_abs_min", "tmin", "t_min", "minus_t_min", "neg_t_min"])
    cols["t_max"] = find_column(fieldnames, ["t_abs_max", "tmax", "t_max", "minus_t_max", "neg_t_max"])

    phi_avg = find_column(fieldnames, ["phiavg", "phi_avg", "phi_average", "phi"], required=False)
    phi_min = find_column(fieldnames, ["phimin", "phi_min"], required=False)
    phi_max = find_column(fieldnames, ["phimax", "phi_max"], required=False)

    if phi_avg is not None:
        cols["phi_avg"] = phi_avg
    elif phi_min is not None and phi_max is not None:
        cols["phi_min"] = phi_min
        cols["phi_max"] = phi_max
    else:
        die(f"{label}: could not find phiavg/phi_avg/phi, or phimin plus phimax columns.")
    #endif

    log(f"{label}: column map:")
    log(f"  valid    -> {cols['valid'] if cols['valid'] else '(none; treating all rows as valid)'}")
    log(f"  xB min   -> {cols['xb_min']}")
    log(f"  xB max   -> {cols['xb_max']}")
    log(f"  Q2 min   -> {cols['q2_min']}")
    log(f"  Q2 max   -> {cols['q2_max']}")
    log(f"  |t| min  -> {cols['t_min']}")
    log(f"  |t| max  -> {cols['t_max']}")

    if "phi_avg" in cols:
        log(f"  phi avg  -> {cols['phi_avg']}")
    else:
        log(f"  phi min  -> {cols['phi_min']}")
        log(f"  phi max  -> {cols['phi_max']}")
    #endif

    return cols


def row_is_valid(row: Dict[str, str], valid_col: str) -> bool:
    if valid_col == "":
        return True
    #endif

    text = str(row.get(valid_col, "")).strip()

    if text == "":
        return True
    #endif

    return to_int(text, default=0) == 1


def row_phi(row: Dict[str, str], cols: Dict[str, str]) -> float:
    if "phi_avg" in cols:
        return to_float(row[cols["phi_avg"]])
    #endif

    return 0.5 * (to_float(row[cols["phi_min"]]) + to_float(row[cols["phi_max"]]))


# ---------------------------------------------------------------------------
# Tuple parsing and stat-weighted combination.
# ---------------------------------------------------------------------------
def parse_tuple_cell(cell: object) -> TupleValue:
    text = str(cell).strip()

    if text == "":
        return TupleValue()
    #endif

    values = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)]

    if len(values) < 2:
        return TupleValue()
    #endif

    value = values[0]
    stat = abs(values[1])
    sys_value = values[2] if len(values) >= 3 else 0.0

    if not math.isfinite(value):
        return TupleValue()
    #endif

    if not math.isfinite(stat) or stat <= 0.0:
        return TupleValue()
    #endif

    if not math.isfinite(sys_value):
        sys_value = 0.0
    #endif

    return TupleValue(ok=True, value=value, stat=stat, sys=sys_value)


def valid_values_only(values: Sequence[TupleValue]) -> List[TupleValue]:
    out: List[TupleValue] = []

    for v in values:
        if v.ok and math.isfinite(v.value) and math.isfinite(v.stat) and v.stat > 0.0:
            out.append(v)
        #endif
    #endfor

    return out


def combine_stat_weighted(values: Sequence[TupleValue]) -> TupleValue:
    valid = valid_values_only(values)

    if not valid:
        return TupleValue()
    #endif

    sum_w = 0.0
    sum_wx = 0.0

    for v in valid:
        w = 1.0 / (v.stat * v.stat)
        sum_w += w
        sum_wx += w * v.value
    #endfor

    if sum_w <= 0.0:
        return TupleValue()
    #endif

    value = sum_wx / sum_w
    stat = 1.0 / math.sqrt(sum_w)

    if not math.isfinite(value) or not math.isfinite(stat) or stat <= 0.0:
        return TupleValue()
    #endif

    return TupleValue(ok=True, value=value, stat=stat, sys=0.0)


# ---------------------------------------------------------------------------
# Point-by-point pass-2 scale systematic.
# ---------------------------------------------------------------------------
def compute_local_pass2_scale(period_values: Sequence[TupleValue]) -> LocalScaleResult:
    valid = valid_values_only(period_values)

    out = LocalScaleResult(n_valid_periods=len(valid))

    if len(valid) < 2:
        return out
    #endif

    ref = combine_stat_weighted(valid)

    if not ref.ok or not math.isfinite(ref.value) or abs(ref.value) <= 0.0:
        return out
    #endif

    ratios: List[float] = []
    ratio_stats: List[float] = []

    for v in valid:
        ratio = v.value / ref.value
        ratio_stat = abs(v.stat / ref.value)

        if math.isfinite(ratio) and math.isfinite(ratio_stat) and ratio_stat > 0.0:
            ratios.append(ratio)
            ratio_stats.append(ratio_stat)
        #endif
    #endfor

    if len(ratios) < 2:
        return out
    #endif

    residuals = [r - 1.0 for r in ratios]
    s_obs = rms(residuals)
    s_stat_exp = rms(ratio_stats)
    s_comb = math.sqrt(max(0.0, s_obs * s_obs - s_stat_exp * s_stat_exp))

    min_ratio = min(ratios)
    max_ratio = max(ratios)

    out.n_valid_periods = len(ratios)
    out.ref_value = ref.value
    out.ref_stat = ref.stat
    out.s_obs = s_obs
    out.s_stat_exp = s_stat_exp
    out.s_comb = s_comb
    out.min_ratio = min_ratio
    out.max_ratio = max_ratio
    out.half_width = 0.5 * (max_ratio - min_ratio)

    return out


def summarize_local_scale_results(
    local_scales: Sequence[float],
    n_less_than_two: int,
    n_rows: int,
    n_central: int,
) -> PointScaleSummary:
    valid_nonzero = [v for v in local_scales if math.isfinite(v) and v > 0.0]

    return PointScaleSummary(
        n_rows=n_rows,
        n_rows_with_central_value=n_central,
        n_rows_with_local_scale=len(valid_nonzero),
        n_rows_with_less_than_two_periods=n_less_than_two,
        mean_s_comb=mean(valid_nonzero),
        rms_s_comb=rms(valid_nonzero),
        median_s_comb=median(valid_nonzero),
        max_s_comb=max(valid_nonzero) if valid_nonzero else 0.0,
    )


def write_pass2_local_scale_summary(output_dir: Path, summary: PointScaleSummary) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "pass2_point_by_point_scale_summary.csv"

    log(f"Pass-2 local scale systematic: writing summary CSV to {path}")

    with path.open("w", newline="") as handle:
        handle.write(
            "rows,rows with central value,rows with nonzero local scale,"
            "rows with fewer than two valid periods,mean s_comb,rms s_comb,"
            "median s_comb,max s_comb,mean percent,rms percent,median percent,max percent\n"
        )
        handle.write(
            f"{summary.n_rows},"
            f"{summary.n_rows_with_central_value},"
            f"{summary.n_rows_with_local_scale},"
            f"{summary.n_rows_with_less_than_two_periods},"
            f"{format_float(summary.mean_s_comb)},"
            f"{format_float(summary.rms_s_comb)},"
            f"{format_float(summary.median_s_comb)},"
            f"{format_float(summary.max_s_comb)},"
            f"{format_float(100.0 * summary.mean_s_comb)},"
            f"{format_float(100.0 * summary.rms_s_comb)},"
            f"{format_float(100.0 * summary.median_s_comb)},"
            f"{format_float(100.0 * summary.max_s_comb)}\n"
        )
    #endwith


def write_pass2_local_scale_points_csv(
    output_dir: Path,
    points: Sequence[Tuple[BinKey, float, LocalScaleResult, float]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "pass2_point_by_point_scale_points.csv"

    log(f"Pass-2 local scale systematic: writing point CSV to {path}")

    with path.open("w", newline="") as handle:
        handle.write(
            "xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phi,"
            "pass2_xs,n_valid_periods,ref_value,ref_stat,"
            "s_obs,s_stat,s_comb,s_comb_percent,min_ratio,max_ratio,half_width\n"
        )

        for key, phi, local, xs in points:
            handle.write(
                f"{format_float(key.xb_min)},"
                f"{format_float(key.xb_max)},"
                f"{format_float(key.q2_min)},"
                f"{format_float(key.q2_max)},"
                f"{format_float(key.t_min)},"
                f"{format_float(key.t_max)},"
                f"{format_float(phi)},"
                f"{format_float(xs)},"
                f"{local.n_valid_periods},"
                f"{format_float(local.ref_value)},"
                f"{format_float(local.ref_stat)},"
                f"{format_float(local.s_obs)},"
                f"{format_float(local.s_stat_exp)},"
                f"{format_float(local.s_comb)},"
                f"{format_float(100.0 * local.s_comb)},"
                f"{format_float(local.min_ratio)},"
                f"{format_float(local.max_ratio)},"
                f"{format_float(local.half_width)}\n"
            )
        #endfor
    #endwith


# ---------------------------------------------------------------------------
# Uncertainty prescriptions.
# ---------------------------------------------------------------------------
def pass2_cross_section_point_uncertainty(xs: float, stat: float, include_point_to_point: bool) -> float:
    if not include_point_to_point:
        return abs(stat)
    #endif

    estimated_syst = PASS2_ESTIMATED_SYSTEMATIC_FRACTION * xs

    return math.sqrt(stat * stat + estimated_syst * estimated_syst)


def pass1_cross_section_point_uncertainty(
    xs: float,
    stat: float,
    syst: float,
    include_point_to_point: bool,
) -> float:
    if not include_point_to_point:
        return abs(stat)
    #endif

    return math.sqrt(stat * stat + syst * syst)


def pass1_cross_section_scale_uncertainty(xs: float) -> float:
    return abs(PASS1_NORMALIZATION_FRACTION * xs)


def pass2_ratio_stat_uncertainty(xs: float, stat: float) -> float:
    if not finite_positive(xs):
        return 0.0
    #endif

    return abs(stat)


def pass1_ratio_stat_uncertainty(xs: float, stat: float) -> float:
    if not finite_positive(xs):
        return 0.0
    #endif

    return abs(stat)


def ratio_scale_uncertainty(ratio: float, pass2_local_scale_frac: float) -> float:
    rel = math.sqrt(
        pass2_local_scale_frac * pass2_local_scale_frac +
        PASS1_NORMALIZATION_FRACTION * PASS1_NORMALIZATION_FRACTION
    )

    return abs(ratio * rel)


# ---------------------------------------------------------------------------
# Input loaders.
# ---------------------------------------------------------------------------
def load_pass2_csv(
    path: Path,
    xs_column: Optional[str],
    output_dir: Path,
    include_point_to_point: bool,
    print_columns: bool = False,
) -> Tuple[Dict[BinKey, List[DataPoint]], PointScaleSummary]:
    t0 = time.time()

    log(f"Pass-2: opening CSV: {path}")

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            die(f"Pass-2 CSV appears empty: {path}")
        #endif

        fieldnames = reader.fieldnames

        log(f"Pass-2: detected {len(fieldnames)} columns.")

        if print_columns:
            log("Pass-2 columns:\n" + "\n".join(f"  [{i}] {name}" for i, name in enumerate(fieldnames)))
        #endif

        cols = detect_common_columns(fieldnames, "Pass-2")
        xs_col = xs_column or find_column(fieldnames, PASS2_CENTRAL_XS_CANDIDATES)

        if "combination sys" in xs_col.lower():
            die(
                "The selected pass-2 cross-section column contains 'combination sys'. "
                "This script must use a central cross-section column, not the combination-sys output column. "
                f"Selected column was: {xs_col}"
            )
        #endif

        required_period_cols = [pass2_period_cross_section_column(period) for period in PASS2_PERIODS_10P6_UNPOL]
        require_columns(fieldnames, required_period_cols, "pass-2 point-by-point local scale systematic")

        period_cols = {
            period: find_column(fieldnames, [pass2_period_cross_section_column(period)])
            for period in PASS2_PERIODS_10P6_UNPOL
        }

        log(f"Pass-2: central cross-section column -> {xs_col}")
        for period in PASS2_PERIODS_10P6_UNPOL:
            log(f"Pass-2: local scale input {period:10s} -> {period_cols[period]}")
        #endfor

        if include_point_to_point:
            log("Pass-2: cross-section vertical error bars -> sqrt(stat^2 + (0.18 * xs)^2).")
        else:
            log("Pass-2: cross-section vertical error bars -> stat only.")
        #endif

        log("Pass-2: cross-section external scale boxes -> local point-by-point s_comb(row) * xs.")
        log("Pass-2: reading rows.")

        out: Dict[BinKey, List[DataPoint]] = {}
        local_scale_points: List[Tuple[BinKey, float, LocalScaleResult, float]] = []
        local_scale_values: List[float] = []

        kept = 0
        skipped_invalid = 0
        skipped_bad_xs = 0
        total_rows = 0
        n_less_than_two = 0
        n_central = 0

        for row in reader:
            total_rows += 1

            if not row_is_valid(row, cols["valid"]):
                skipped_invalid += 1
                continue
            #endif

            key = make_bin_key(row, cols)
            phi = row_phi(row, cols)

            tuple_value = parse_tuple_cell(row.get(xs_col, ""))

            if not tuple_value.ok or not finite_positive(tuple_value.value):
                skipped_bad_xs += 1
                continue
            #endif

            n_central += 1

            period_values = [
                parse_tuple_cell(row.get(period_cols[period], ""))
                for period in PASS2_PERIODS_10P6_UNPOL
            ]

            local = compute_local_pass2_scale(period_values)

            if local.n_valid_periods < 2:
                n_less_than_two += 1
            #endif

            xs = tuple_value.value
            stat = tuple_value.stat
            local_scale_frac = local.s_comb
            scale_syst_abs = abs(local_scale_frac * xs)
            point_err = pass2_cross_section_point_uncertainty(
                xs=xs,
                stat=stat,
                include_point_to_point=include_point_to_point,
            )

            point = DataPoint(
                phi=phi,
                xs=xs,
                stat=stat,
                err_low=point_err,
                err_high=point_err,
                source="pass2",
                scale_syst_abs=scale_syst_abs,
                scale_syst_frac=local_scale_frac,
                local_scale_n_periods=local.n_valid_periods,
                local_scale_s_obs=local.s_obs,
                local_scale_s_stat=local.s_stat_exp,
                local_scale_half_width=local.half_width,
            )

            out.setdefault(key, []).append(point)
            local_scale_points.append((key, phi, local, xs))
            local_scale_values.append(local_scale_frac)
            kept += 1
        #endfor
    #endwith

    for points in out.values():
        points.sort(key=lambda p: p.phi)
    #endfor

    summary = summarize_local_scale_results(
        local_scales=local_scale_values,
        n_less_than_two=n_less_than_two,
        n_rows=total_rows,
        n_central=n_central,
    )

    write_pass2_local_scale_points_csv(output_dir, local_scale_points)
    write_pass2_local_scale_summary(output_dir, summary)

    dt = time.time() - t0

    log(f"Pass-2: finished reading in {format_seconds(dt)}.")
    log(f"Pass-2: total rows={total_rows}, kept={kept}, skipped invalid-bin={skipped_invalid}, skipped bad/nonpositive xs={skipped_bad_xs}.")
    log(f"Pass-2: rows with fewer than two valid 10.6 GeV period entries={n_less_than_two}.")
    log(
        "Pass-2: local s_comb summary over nonzero local scales -> "
        f"mean={summary.mean_s_comb:.10g} ({100.0 * summary.mean_s_comb:.4g}%), "
        f"median={summary.median_s_comb:.10g} ({100.0 * summary.median_s_comb:.4g}%), "
        f"rms={summary.rms_s_comb:.10g} ({100.0 * summary.rms_s_comb:.4g}%), "
        f"max={summary.max_s_comb:.10g} ({100.0 * summary.max_s_comb:.4g}%)."
    )
    log(f"Pass-2: unique (xB,Q2,|t|) bins with data={len(out)}.")

    return out, summary


def load_pass1_csv(
    path: Path,
    include_point_to_point: bool,
    print_columns: bool = False,
) -> Dict[BinKey, List[DataPoint]]:
    t0 = time.time()

    log(f"Pass-1: opening CSV: {path}")

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            die(f"Pass-1 CSV appears empty: {path}")
        #endif

        fieldnames = reader.fieldnames

        log(f"Pass-1: detected {len(fieldnames)} columns.")

        if print_columns:
            log("Pass-1 columns:\n" + "\n".join(f"  [{i}] {name}" for i, name in enumerate(fieldnames)))
        #endif

        cols = detect_common_columns(fieldnames, "Pass-1")

        xs_col = find_column(fieldnames, [PASS1_XS_COL])
        stat_col = find_column(fieldnames, [PASS1_STAT_COL])
        syst_up_col = find_column(fieldnames, [PASS1_SYST_UP_COL])
        syst_dn_col = find_column(fieldnames, [PASS1_SYST_DN_COL])

        log(f"Pass-1: cross-section column -> {xs_col}")
        log(f"Pass-1: stat uncertainty column -> {stat_col}")
        log(f"Pass-1: syst up uncertainty column -> {syst_up_col}")
        log(f"Pass-1: syst down uncertainty column -> {syst_dn_col}")

        if include_point_to_point:
            log("Pass-1: cross-section vertical error bars -> sqrt(stat^2 + provided_syst^2).")
        else:
            log("Pass-1: cross-section vertical error bars -> stat only.")
        #endif

        log("Pass-1: cross-section external scale boxes -> 0.31 * xs.")
        log("Pass-1: reading rows.")

        out: Dict[BinKey, List[DataPoint]] = {}
        kept = 0
        skipped_invalid = 0
        skipped_bad_xs = 0
        total_rows = 0

        for row in reader:
            total_rows += 1

            if not row_is_valid(row, cols["valid"]):
                skipped_invalid += 1
                continue
            #endif

            key = make_bin_key(row, cols)
            phi = row_phi(row, cols)

            xs = to_float(row.get(xs_col, ""))
            stat = abs(to_float(row.get(stat_col, "")))
            syst_up = abs(to_float(row.get(syst_up_col, "")))
            syst_dn = abs(to_float(row.get(syst_dn_col, "")))

            if not finite_positive(xs):
                skipped_bad_xs += 1
                continue
            #endif

            err_high = pass1_cross_section_point_uncertainty(
                xs=xs,
                stat=stat,
                syst=syst_up,
                include_point_to_point=include_point_to_point,
            )
            err_low = pass1_cross_section_point_uncertainty(
                xs=xs,
                stat=stat,
                syst=syst_dn,
                include_point_to_point=include_point_to_point,
            )

            if err_low <= 0.0:
                err_low = err_high
            #endif

            if err_high <= 0.0:
                err_high = err_low
            #endif

            scale_syst_abs = pass1_cross_section_scale_uncertainty(xs=xs)

            point = DataPoint(
                phi=phi,
                xs=xs,
                stat=stat,
                err_low=err_low,
                err_high=err_high,
                source="pass1",
                scale_syst_abs=scale_syst_abs,
                scale_syst_frac=PASS1_NORMALIZATION_FRACTION,
            )

            out.setdefault(key, []).append(point)
            kept += 1
        #endfor
    #endwith

    for points in out.values():
        points.sort(key=lambda p: p.phi)
    #endfor

    dt = time.time() - t0

    log(f"Pass-1: finished reading in {format_seconds(dt)}.")
    log(f"Pass-1: total rows={total_rows}, kept={kept}, skipped invalid-bin={skipped_invalid}, skipped bad/nonpositive xs={skipped_bad_xs}.")
    log(f"Pass-1: unique (xB,Q2,|t|) bins with data={len(out)}.")

    return out


# ---------------------------------------------------------------------------
# Model calculation helpers.
# ---------------------------------------------------------------------------
def run_command_capture(cmd: Sequence[str], env: Optional[Dict[str, str]] = None) -> str:
    completed = subprocess.run(
        list(cmd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        check=False,
    )

    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        stdout = completed.stdout.strip()

        raise RuntimeError(
            f"Command failed with code {completed.returncode}: {' '.join(shlex.quote(c) for c in cmd)}\n"
            f"stdout:\n{stdout}\n"
            f"stderr:\n{stderr}"
        )
    #endif

    return completed.stdout


def parse_last_numeric_token(text: str, which_from_end: int = 1) -> float:
    lines = [line.strip() for line in text.splitlines() if line.strip()]

    if len(lines) < which_from_end:
        return 0.0
    #endif

    target = lines[-which_from_end]
    numbers = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", target)

    if not numbers:
        return 0.0
    #endif

    return float(numbers[-1])


def bh_xs(xb: float, q2: float, t_abs: float, phi_deg: float, cfg: ModelConfig) -> float:
    dvcsgen = Path(cfg.dvcsgen_dir) / "dvcsgen"
    phi_rad = phi_deg * math.pi / 180.0

    cmd = [
        str(dvcsgen),
        "--beam", f"{cfg.e_beam:.12g}",
        "--x", f"{xb:.12g}", f"{xb:.12g}",
        "--q2", f"{q2:.12g}", f"{q2:.12g}",
        "--t", f"{t_abs:.12g}", f"{t_abs:.12g}",
        "--bh", "1",
        "--phi", f"{phi_rad:.12g}",
        "--ycol", "0.0001",
    ]

    env = os.environ.copy()
    env["PATH"] = f"{cfg.dvcsgen_dir}:{env.get('PATH', '')}"
    env["CLASDVCS_PDF"] = cfg.dvcsgen_dir

    out = run_command_capture(cmd, env=env)

    return parse_last_numeric_token(out, which_from_end=1)


def km15_xs(xb: float, q2: float, t_abs: float, phi_deg: float, cfg: ModelConfig) -> float:
    python_exe = cfg.py_km15 if cfg.py_km15 else sys.executable

    cmd = [
        python_exe,
        cfg.km15_cli,
        f"{xb:.12g}",
        f"{q2:.12g}",
        f"{t_abs:.12g}",
        f"{phi_deg:.12g}",
        f"{cfg.e_beam:.12g}",
        "0",
        "XS",
    ]

    env = os.environ.copy()
    env.pop("PYTHONPATH", None)

    out = run_command_capture(cmd, env=env).strip()

    return to_float(out, default=0.0)


def model_cache_key(key: BinKey, cfg: ModelConfig) -> str:
    return (
        f"E={cfg.e_beam:.8g}|nphi={cfg.phi_dense}|"
        f"xb={key.xb_min:.8g},{key.xb_max:.8g}|"
        f"q2={key.q2_min:.8g},{key.q2_max:.8g}|"
        f"t={key.t_min:.8g},{key.t_max:.8g}"
    )


def load_model_cache(path: Path) -> Dict[str, Dict[str, List[float]]]:
    t0 = time.time()

    log(f"Cache: checking for model cache at {path}")

    if not path.exists():
        log("Cache: no existing cache file found.")
        return {}
    #endif

    log("Cache: existing cache file found; reading JSON.")

    try:
        with path.open("r") as handle:
            raw = json.load(handle)
        #endwith

        if isinstance(raw, dict):
            good: Dict[str, Dict[str, List[float]]] = {}
            rejected = 0

            for key, entry in raw.items():
                if not isinstance(entry, dict):
                    rejected += 1
                    continue
                #endif

                if "phi" not in entry or "bh" not in entry or "km15" not in entry:
                    rejected += 1
                    continue
                #endif

                good[key] = {
                    "phi": list(entry["phi"]),
                    "bh": list(entry["bh"]),
                    "km15": list(entry["km15"]),
                }
            #endfor

            dt = time.time() - t0

            log(f"Cache: loaded {len(good)} valid entries in {format_seconds(dt)}; rejected {rejected} malformed entries.")
            return good
        #endif
    except Exception as exc:
        warn(f"Cache: could not read model cache {path}: {exc}")
    #endtry

    return {}


def save_model_cache(path: Path, cache: Dict[str, Dict[str, List[float]]]) -> None:
    t0 = time.time()

    log(f"Cache: writing {len(cache)} model entries to {path}")

    path.parent.mkdir(parents=True, exist_ok=True)

    tmp_path = path.with_suffix(path.suffix + ".tmp")

    with tmp_path.open("w") as handle:
        json.dump(cache, handle, indent=2, sort_keys=True)
    #endwith

    tmp_path.replace(path)

    dt = time.time() - t0

    log(f"Cache: finished writing cache in {format_seconds(dt)}.")


def make_model_curves_from_cache_entry(entry: Dict[str, List[float]]) -> ModelCurves:
    return ModelCurves(phi=list(entry["phi"]), bh=list(entry["bh"]), km15=list(entry["km15"]))


def make_model_cache_entry(curves: ModelCurves) -> Dict[str, List[float]]:
    return {"phi": curves.phi, "bh": curves.bh, "km15": curves.km15}


def compute_model_curves_without_cache(key: BinKey, cfg: ModelConfig) -> ModelCurves:
    xb = 0.5 * (key.xb_min + key.xb_max)
    q2 = 0.5 * (key.q2_min + key.q2_max)
    t_abs = 0.5 * (key.t_min + key.t_max)

    if cfg.phi_dense <= 1:
        phi_grid = [0.0]
    else:
        phi_grid = [360.0 * i / (cfg.phi_dense - 1) for i in range(cfg.phi_dense)]
    #endif

    bh_values: List[float] = []
    km15_values: List[float] = []

    for i_phi, phi in enumerate(phi_grid, start=1):
        if cfg.verbose_worker_models and (i_phi == 1 or i_phi == len(phi_grid) or i_phi % 25 == 0):
            print(
                f"[pass2-vs-pass1][worker {os.getpid()}] "
                f"model point {i_phi}/{len(phi_grid)} for "
                f"xB=[{key.xb_min},{key.xb_max}], Q2=[{key.q2_min},{key.q2_max}], "
                f"|t|=[{key.t_min},{key.t_max}], phi={phi:.1f}",
                flush=True,
            )
        #endif

        bh_values.append(bh_xs(xb, q2, t_abs, phi, cfg))
        km15_values.append(km15_xs(xb, q2, t_abs, phi, cfg))
    #endfor

    return ModelCurves(phi=phi_grid, bh=bh_values, km15=km15_values)


# ---------------------------------------------------------------------------
# Ratio and matching helpers.
# ---------------------------------------------------------------------------
def find_nearest_pass1_point(phi: float, pass1_points: Sequence[DataPoint]) -> Optional[DataPoint]:
    if not pass1_points:
        return None
    #endif

    best_point: Optional[DataPoint] = None
    best_delta = float("inf")

    for point in pass1_points:
        delta = abs(point.phi - phi)

        if delta < best_delta:
            best_delta = delta
            best_point = point
        #endif
    #endfor

    if best_point is None:
        return None
    #endif

    if best_delta > PASS1_PASS2_RATIO_MATCH_TOLERANCE_DEG:
        return None
    #endif

    return best_point


def pass2_over_pass1_ratio_points(panel: PanelData) -> List[RatioPoint]:
    ratio_points: List[RatioPoint] = []

    for p2 in panel.pass2:
        p1 = find_nearest_pass1_point(p2.phi, panel.pass1)

        if p1 is None:
            continue
        #endif

        if not finite_positive(p2.xs) or not finite_positive(p1.xs):
            continue
        #endif

        ratio = p2.xs / p1.xs

        p2_stat_err = pass2_ratio_stat_uncertainty(xs=p2.xs, stat=p2.stat)
        p1_stat_err = pass1_ratio_stat_uncertainty(xs=p1.xs, stat=p1.stat)

        rel_p2_stat = p2_stat_err / p2.xs
        rel_p1_stat = p1_stat_err / p1.xs
        rel_ratio_stat = math.sqrt(rel_p2_stat * rel_p2_stat + rel_p1_stat * rel_p1_stat)
        ratio_stat_err = abs(ratio * rel_ratio_stat)

        pass2_scale_frac = p2.scale_syst_frac
        scale_err = ratio_scale_uncertainty(ratio=ratio, pass2_local_scale_frac=pass2_scale_frac)

        ratio_points.append(
            RatioPoint(
                phi=p2.phi,
                ratio=ratio,
                stat_err_low=ratio_stat_err,
                stat_err_high=ratio_stat_err,
                scale_err=scale_err,
                pass2_scale_frac=pass2_scale_frac,
                label="pass2/pass1",
            )
        )
    #endfor

    return ratio_points


def ratio_axis_limits(ratio_sets: Sequence[Sequence[RatioPoint]]) -> Tuple[float, float]:
    vals: List[float] = []

    for ratio_points in ratio_sets:
        for point in ratio_points:
            point_extent = max(point.stat_err_low, point.stat_err_high, point.scale_err)
            vals.append(point.ratio - point_extent)
            vals.append(point.ratio + point_extent)
        #endfor
    #endfor

    vals = [v for v in vals if math.isfinite(v)]

    if not vals:
        return 0.0, 2.0
    #endif

    ymin = min(vals)
    ymax = max(vals)

    ymin = min(ymin, 1.0)
    ymax = max(ymax, 1.0)

    span = ymax - ymin

    if span <= 0.0:
        return ymin - 0.5, ymax + 0.5
    #endif

    pad = 0.15 * span

    return ymin - pad, ymax + pad


# ---------------------------------------------------------------------------
# Stat-only pass-1/pass-2 pull diagnostic.
# ---------------------------------------------------------------------------
def build_stat_only_comparison_points(panels: Sequence[PanelData]) -> List[StatOnlyComparisonPoint]:
    points: List[StatOnlyComparisonPoint] = []

    for panel in panels:
        for p2 in panel.pass2:
            p1 = find_nearest_pass1_point(p2.phi, panel.pass1)

            if p1 is None:
                continue
            #endif

            if not finite_positive(p2.xs) or not finite_positive(p1.xs):
                continue
            #endif

            if not finite_positive(p2.stat) or not finite_positive(p1.stat):
                continue
            #endif

            points.append(
                StatOnlyComparisonPoint(
                    key=panel.key,
                    phi_pass2=p2.phi,
                    phi_pass1=p1.phi,
                    pass2_xs=p2.xs,
                    pass2_stat=p2.stat,
                    pass1_xs=p1.xs,
                    pass1_stat=p1.stat,
                )
            )
        #endfor
    #endfor

    return points


def stat_only_pull(point: StatOnlyComparisonPoint, pass1_norm: float) -> float:
    denom2 = point.pass2_stat * point.pass2_stat
    denom2 += (pass1_norm * point.pass1_stat) * (pass1_norm * point.pass1_stat)

    if denom2 <= 0.0 or not math.isfinite(denom2):
        return float("nan")
    #endif

    numerator = point.pass2_xs - pass1_norm * point.pass1_xs

    return numerator / math.sqrt(denom2)


def stat_only_chi2(points: Sequence[StatOnlyComparisonPoint], pass1_norm: float) -> float:
    total = 0.0

    for point in points:
        pull = stat_only_pull(point, pass1_norm)

        if math.isfinite(pull):
            total += pull * pull
        #endif
    #endfor

    return total


def bounded_minimize_golden_section(
    func: Callable[[float], float],
    xmin: float,
    xmax: float,
    tolerance: float = 1.0e-8,
    max_iter: int = 200,
) -> Tuple[float, float]:
    gr = 0.5 * (math.sqrt(5.0) - 1.0)

    a = xmin
    b = xmax

    c = b - gr * (b - a)
    d = a + gr * (b - a)

    fc = func(c)
    fd = func(d)

    for _ in range(max_iter):
        if abs(b - a) <= tolerance:
            break
        #endif

        if fc < fd:
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            fc = func(c)
        else:
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            fd = func(d)
        #endif
    #endfor

    xbest = 0.5 * (a + b)
    fbest = func(xbest)

    endpoint_min = func(xmin)
    endpoint_max = func(xmax)

    if endpoint_min < fbest and endpoint_min <= endpoint_max:
        return xmin, endpoint_min
    #endif

    if endpoint_max < fbest and endpoint_max < endpoint_min:
        return xmax, endpoint_max
    #endif

    return xbest, fbest


def percent_within(pulls: Sequence[float], threshold: float) -> float:
    vals = [p for p in pulls if math.isfinite(p)]

    if not vals:
        return 0.0
    #endif

    n_within = sum(1 for p in vals if abs(p) <= threshold)

    return 100.0 * float(n_within) / float(len(vals))


def summarize_stat_only_pulls(
    points: Sequence[StatOnlyComparisonPoint],
    best_pass1_norm: float,
) -> StatOnlyPullSummary:
    nominal_pulls = [stat_only_pull(point, 1.0) for point in points]
    best_pulls = [stat_only_pull(point, best_pass1_norm) for point in points]

    nominal_pulls = [p for p in nominal_pulls if math.isfinite(p)]
    best_pulls = [p for p in best_pulls if math.isfinite(p)]

    n_points = min(len(nominal_pulls), len(best_pulls))
    chi2_nominal = sum(p * p for p in nominal_pulls)
    chi2_best = sum(p * p for p in best_pulls)

    ndf_nominal = n_points
    ndf_best = max(1, n_points - 1)

    return StatOnlyPullSummary(
        n_points=n_points,
        best_pass1_norm=best_pass1_norm,
        best_norm_min=PULL_FIT_NORM_MIN,
        best_norm_max=PULL_FIT_NORM_MAX,
        chi2_nominal=chi2_nominal,
        chi2_best=chi2_best,
        ndf_nominal=ndf_nominal,
        ndf_best=ndf_best,
        chi2_ndf_nominal=chi2_nominal / float(ndf_nominal) if ndf_nominal > 0 else 0.0,
        chi2_ndf_best=chi2_best / float(ndf_best) if ndf_best > 0 else 0.0,
        pct_within_1sigma_nominal=percent_within(nominal_pulls, 1.0),
        pct_within_3sigma_nominal=percent_within(nominal_pulls, 3.0),
        pct_within_1sigma_best=percent_within(best_pulls, 1.0),
        pct_within_3sigma_best=percent_within(best_pulls, 3.0),
        mean_pull_nominal=mean(nominal_pulls),
        rms_pull_nominal=rms(nominal_pulls),
        mean_pull_best=mean(best_pulls),
        rms_pull_best=rms(best_pulls),
    )


def write_stat_only_pull_points_csv(
    output_dir: Path,
    points: Sequence[StatOnlyComparisonPoint],
    best_pass1_norm: float,
) -> None:
    path = output_dir / "pass1_pass2_stat_only_pull_points.csv"

    log(f"Stat-only pull diagnostic: writing point CSV to {path}")

    with path.open("w", newline="") as handle:
        handle.write(
            "xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,"
            "phi_pass2,phi_pass1,pass2_xs,pass2_stat,pass1_xs,pass1_stat,"
            "pull_nominal,pull_best_norm,chi2_nominal,chi2_best_norm,best_pass1_norm\n"
        )

        for point in points:
            pull_nominal = stat_only_pull(point, 1.0)
            pull_best = stat_only_pull(point, best_pass1_norm)

            handle.write(
                f"{format_float(point.key.xb_min)},"
                f"{format_float(point.key.xb_max)},"
                f"{format_float(point.key.q2_min)},"
                f"{format_float(point.key.q2_max)},"
                f"{format_float(point.key.t_min)},"
                f"{format_float(point.key.t_max)},"
                f"{format_float(point.phi_pass2)},"
                f"{format_float(point.phi_pass1)},"
                f"{format_float(point.pass2_xs)},"
                f"{format_float(point.pass2_stat)},"
                f"{format_float(point.pass1_xs)},"
                f"{format_float(point.pass1_stat)},"
                f"{format_float(pull_nominal)},"
                f"{format_float(pull_best)},"
                f"{format_float(pull_nominal * pull_nominal)},"
                f"{format_float(pull_best * pull_best)},"
                f"{format_float(best_pass1_norm)}\n"
            )
        #endfor
    #endwith


def write_stat_only_pull_summary_csv(output_dir: Path, summary: StatOnlyPullSummary) -> None:
    path = output_dir / "pass1_pass2_stat_only_pull_summary.csv"

    log(f"Stat-only pull diagnostic: writing summary CSV to {path}")

    with path.open("w", newline="") as handle:
        handle.write(
            "n_points,best_pass1_norm,best_pass1_norm_percent,"
            "norm_min,norm_max,chi2_nominal,ndf_nominal,chi2_ndf_nominal,"
            "chi2_best,ndf_best,chi2_ndf_best,"
            "pct_within_1sigma_nominal,pct_within_3sigma_nominal,"
            "pct_within_1sigma_best,pct_within_3sigma_best,"
            "mean_pull_nominal,rms_pull_nominal,mean_pull_best,rms_pull_best\n"
        )

        handle.write(
            f"{summary.n_points},"
            f"{format_float(summary.best_pass1_norm)},"
            f"{format_float(100.0 * (summary.best_pass1_norm - 1.0))},"
            f"{format_float(summary.best_norm_min)},"
            f"{format_float(summary.best_norm_max)},"
            f"{format_float(summary.chi2_nominal)},"
            f"{summary.ndf_nominal},"
            f"{format_float(summary.chi2_ndf_nominal)},"
            f"{format_float(summary.chi2_best)},"
            f"{summary.ndf_best},"
            f"{format_float(summary.chi2_ndf_best)},"
            f"{format_float(summary.pct_within_1sigma_nominal)},"
            f"{format_float(summary.pct_within_3sigma_nominal)},"
            f"{format_float(summary.pct_within_1sigma_best)},"
            f"{format_float(summary.pct_within_3sigma_best)},"
            f"{format_float(summary.mean_pull_nominal)},"
            f"{format_float(summary.rms_pull_nominal)},"
            f"{format_float(summary.mean_pull_best)},"
            f"{format_float(summary.rms_pull_best)}\n"
        )
    #endwith


def draw_stat_only_pull_distribution(
    output_dir: Path,
    points: Sequence[StatOnlyComparisonPoint],
    summary: StatOnlyPullSummary,
) -> None:
    path = output_dir / "pass1_pass2_stat_only_pull_distribution.png"

    nominal_pulls = [stat_only_pull(point, 1.0) for point in points]
    best_pulls = [stat_only_pull(point, summary.best_pass1_norm) for point in points]

    nominal_pulls = [p for p in nominal_pulls if math.isfinite(p)]
    best_pulls = [p for p in best_pulls if math.isfinite(p)]

    if not nominal_pulls or not best_pulls:
        warn("Stat-only pull diagnostic: no valid pulls; skipping histogram.")
        return
    #endif

    finite_pulls = nominal_pulls + best_pulls

    pmin = min(finite_pulls)
    pmax = max(finite_pulls)

    xmin = min(-5.0, pmin)
    xmax = max(+5.0, pmax)

    if xmin == xmax:
        xmin -= 1.0
        xmax += 1.0
    #endif

    span = xmax - xmin
    xmin -= 0.05 * span
    xmax += 0.05 * span

    n_bins = 80
    fig, ax = plt.subplots(figsize=(9.4, 6.6))

    ax.hist(
        nominal_pulls,
        bins=n_bins,
        range=(xmin, xmax),
        histtype="step",
        linewidth=1.8,
        label="M = 1",
    )

    ax.hist(
        best_pulls,
        bins=n_bins,
        range=(xmin, xmax),
        histtype="step",
        linewidth=1.8,
        label=f"Best M = {summary.best_pass1_norm:.5f}",
    )

    ax.axvline(0.0, linewidth=1.0, linestyle="-")
    ax.axvline(-1.0, linewidth=1.0, linestyle="--")
    ax.axvline(+1.0, linewidth=1.0, linestyle="--")
    ax.axvline(-3.0, linewidth=1.0, linestyle=":")
    ax.axvline(+3.0, linewidth=1.0, linestyle=":")

    text = (
        "Stat-only pass-2/pass-1 pulls\n"
        rf"$z=(\sigma_{{p2}}-M\sigma_{{p1}})/"
        rf"\sqrt{{\delta_{{p2,stat}}^2+(M\delta_{{p1,stat}})^2}}$" "\n"
        f"Matched points: {summary.n_points}\n"
        f"Best pass-1 M in [{PULL_FIT_NORM_MIN:.2f}, {PULL_FIT_NORM_MAX:.2f}]: "
        f"{summary.best_pass1_norm:.5f} "
        f"({100.0 * (summary.best_pass1_norm - 1.0):+.2f}%)\n"
        f"χ²/ndf M=1: {summary.chi2_ndf_nominal:.3g} "
        f"({summary.chi2_nominal:.3g}/{summary.ndf_nominal})\n"
        f"χ²/ndf best: {summary.chi2_ndf_best:.3g} "
        f"({summary.chi2_best:.3g}/{summary.ndf_best})\n"
        f"|z|≤1: {summary.pct_within_1sigma_nominal:.1f}% → "
        f"{summary.pct_within_1sigma_best:.1f}%\n"
        f"|z|≤3: {summary.pct_within_3sigma_nominal:.1f}% → "
        f"{summary.pct_within_3sigma_best:.1f}%"
    )

    ax.text(
        0.97,
        0.97,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )

    ax.set_xlabel("Stat-only pull")
    ax.set_ylabel("Point count")
    ax.set_title("Pass-2 vs pass-1 point-by-point stat-only pull distribution")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="upper left", fontsize=10, frameon=True)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Stat-only pull diagnostic: wrote histogram to {path}")


def make_stat_only_pull_diagnostic(
    output_dir: Path,
    panels: Sequence[PanelData],
) -> StatOnlyPullSummary:
    t0 = time.time()

    log("Stat-only pull diagnostic: building matched pass-2/pass-1 point list.")
    points = build_stat_only_comparison_points(panels)

    if not points:
        warn("Stat-only pull diagnostic: no matched points found.")
        summary = StatOnlyPullSummary()
        write_stat_only_pull_summary_csv(output_dir, summary)
        return summary
    #endif

    log(f"Stat-only pull diagnostic: matched points={len(points)}.")
    log(
        "Stat-only pull diagnostic: minimizing global chi2 with pass-1 normalization "
        f"M constrained to [{PULL_FIT_NORM_MIN:.3f}, {PULL_FIT_NORM_MAX:.3f}]."
    )

    best_norm, best_chi2 = bounded_minimize_golden_section(
        func=lambda norm: stat_only_chi2(points, norm),
        xmin=PULL_FIT_NORM_MIN,
        xmax=PULL_FIT_NORM_MAX,
        tolerance=1.0e-10,
        max_iter=300,
    )

    summary = summarize_stat_only_pulls(points=points, best_pass1_norm=best_norm)

    log(
        "Stat-only pull diagnostic: "
        f"best pass-1 M={summary.best_pass1_norm:.10g} "
        f"({100.0 * (summary.best_pass1_norm - 1.0):+.5g}%), "
        f"chi2={best_chi2:.10g}, "
        f"chi2/ndf={summary.chi2_ndf_best:.10g}."
    )
    log(
        "Stat-only pull diagnostic: "
        f"M=1 chi2/ndf={summary.chi2_ndf_nominal:.10g}, "
        f"|z|<=1={summary.pct_within_1sigma_nominal:.4g}%, "
        f"|z|<=3={summary.pct_within_3sigma_nominal:.4g}%."
    )
    log(
        "Stat-only pull diagnostic: "
        f"best M |z|<=1={summary.pct_within_1sigma_best:.4g}%, "
        f"|z|<=3={summary.pct_within_3sigma_best:.4g}%."
    )

    write_stat_only_pull_points_csv(output_dir=output_dir, points=points, best_pass1_norm=best_norm)
    write_stat_only_pull_summary_csv(output_dir=output_dir, summary=summary)
    draw_stat_only_pull_distribution(output_dir=output_dir, points=points, summary=summary)

    dt = time.time() - t0
    log(f"Stat-only pull diagnostic: finished in {format_seconds(dt)}.")

    return summary


# ---------------------------------------------------------------------------
# Plotting.
# ---------------------------------------------------------------------------
def y_limits(
    pass2: Sequence[DataPoint],
    pass1: Sequence[DataPoint],
    curves: Optional[ModelCurves],
    logy: bool,
    log_y_min: float,
) -> Tuple[float, float]:
    vals: List[float] = []

    floor = log_y_min if logy and log_y_min > 0.0 else 1e-30

    for p in list(pass2) + list(pass1):
        point_extent = max(p.err_low, p.err_high, p.scale_syst_abs)
        vals.append(max(floor, p.xs - point_extent))
        vals.append(max(floor, p.xs + point_extent))
        vals.append(max(floor, p.xs))
    #endfor

    if curves is not None:
        vals.extend(v for v in curves.bh if finite_positive(v))
        vals.extend(v for v in curves.km15 if finite_positive(v))
    #endif

    vals = [v for v in vals if finite_positive(v)]

    if not vals:
        return floor, 1.0
    #endif

    ymin = min(vals)
    ymax = max(vals)

    if logy:
        ymin = max(log_y_min, ymin)
    #endif

    if ymin <= 0.0 or not math.isfinite(ymin):
        ymin = floor
    #endif

    if ymax <= ymin:
        ymax = 10.0 * ymin
    #endif

    return ymin, 1.75 * ymax


def draw_scale_boxes(
    ax,
    x_values: Sequence[float],
    y_values: Sequence[float],
    scale_values: Sequence[float],
    half_width: float,
    label: str,
    color,
    logy: bool,
    log_y_min: float,
) -> None:
    first = True

    for x, y, scale in zip(x_values, y_values, scale_values):
        if not math.isfinite(x) or not math.isfinite(y) or not math.isfinite(scale):
            continue
        #endif

        if scale <= 0.0:
            continue
        #endif

        bottom = y - scale
        top = y + scale

        if logy:
            bottom = max(log_y_min, bottom)
        #endif

        if top <= bottom:
            continue
        #endif

        patch = Rectangle(
            (x - half_width, bottom),
            2.0 * half_width,
            top - bottom,
            facecolor=color,
            edgecolor=color,
            alpha=SCALE_BOX_ALPHA,
            linewidth=1.0,
            label=label if first else None,
            zorder=1,
        )

        ax.add_patch(patch)
        first = False
    #endfor


def draw_cross_section_axis(
    ax,
    panel: PanelData,
    curves: Optional[ModelCurves],
    pass2_label: str,
    pass1_label: str,
    logy: bool,
    log_y_min: float,
    include_title: bool,
    include_point_to_point: bool,
) -> None:
    if curves is not None:
        ax.plot(curves.phi, curves.bh, linewidth=2.0, linestyle="-", label="BH", zorder=2)
        ax.plot(curves.phi, curves.km15, linewidth=2.0, linestyle="--", label="KM15", zorder=2)
    #endif

    pass1_error_label = "stat ⊕ syst" if include_point_to_point else "stat only"
    pass2_error_label = "stat ⊕ 18% point-to-point" if include_point_to_point else "stat only"

    if panel.pass1:
        x = [p.phi + PASS1_PHI_OFFSET_DEG for p in panel.pass1]
        y = [p.xs for p in panel.pass1]
        yerr = [[p.err_low for p in panel.pass1], [p.err_high for p in panel.pass1]]

        pass1_container = ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="s",
            markersize=5,
            capsize=2,
            linewidth=1.2,
            linestyle="None",
            label=f"{pass1_label}: {pass1_error_label}",
            zorder=4,
        )

        pass1_color = pass1_container.lines[0].get_color()

        draw_scale_boxes(
            ax=ax,
            x_values=x,
            y_values=y,
            scale_values=[p.scale_syst_abs for p in panel.pass1],
            half_width=SCALE_BOX_HALF_WIDTH_DEG,
            label=f"{pass1_label}: 31% norm box",
            color=pass1_color,
            logy=logy,
            log_y_min=log_y_min,
        )
    #endif

    if panel.pass2:
        x = [p.phi + PASS2_PHI_OFFSET_DEG for p in panel.pass2]
        y = [p.xs for p in panel.pass2]
        yerr = [p.err_high for p in panel.pass2]

        pass2_container = ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            markersize=5,
            capsize=2,
            linewidth=1.2,
            linestyle="None",
            label=f"{pass2_label}: {pass2_error_label}",
            zorder=5,
        )

        pass2_color = pass2_container.lines[0].get_color()

        draw_scale_boxes(
            ax=ax,
            x_values=x,
            y_values=y,
            scale_values=[p.scale_syst_abs for p in panel.pass2],
            half_width=SCALE_BOX_HALF_WIDTH_DEG,
            label=f"{pass2_label}: local period-spread box",
            color=pass2_color,
            logy=logy,
            log_y_min=log_y_min,
        )
    #endif

    key = panel.key

    if include_title:
        title = (
            rf"$x_B \in [{key.xb_min:.3g}, {key.xb_max:.3g}]$, "
            rf"$Q^2 \in [{key.q2_min:.3g}, {key.q2_max:.3g}]$ (GeV$^2$), "
            rf"$|t| \in [{key.t_min:.3g}, {key.t_max:.3g}]$ (GeV$^2$)"
        )

        ax.set_title(title)
    #endif

    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$d\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ (pb/GeV$^4$/rad)")
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])

    ymin, ymax = y_limits(panel.pass2, panel.pass1, curves, logy=logy, log_y_min=log_y_min)
    ax.set_ylim(ymin, ymax)

    if logy:
        ax.set_yscale("log")
    #endif

    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best", fontsize=8, frameon=True)


def draw_ratio_axis(ax, panel: PanelData) -> None:
    p2_over_p1 = pass2_over_pass1_ratio_points(panel)

    if p2_over_p1:
        container = ax.errorbar(
            [p.phi for p in p2_over_p1],
            [p.ratio for p in p2_over_p1],
            yerr=[[p.stat_err_low for p in p2_over_p1], [p.stat_err_high for p in p2_over_p1]],
            fmt="o",
            markersize=5,
            capsize=2,
            linewidth=1.2,
            linestyle="None",
            label="pass-2 / pass-1: stat only",
            zorder=5,
        )

        color = container.lines[0].get_color()

        draw_scale_boxes(
            ax=ax,
            x_values=[p.phi for p in p2_over_p1],
            y_values=[p.ratio for p in p2_over_p1],
            scale_values=[p.scale_err for p in p2_over_p1],
            half_width=SCALE_BOX_HALF_WIDTH_DEG,
            label="scale box: local pass-2 ⊕ 31% pass-1",
            color=color,
            logy=False,
            log_y_min=1.0e-30,
        )
    #endif

    ax.axhline(1.0, linewidth=1.2, linestyle="--", zorder=2)

    ymin, ymax = ratio_axis_limits([p2_over_p1])

    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel("pass-2 / pass-1")
    ax.set_title("Ratio")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best", fontsize=8, frameon=True)


def draw_panel(
    panel: PanelData,
    curves: Optional[ModelCurves],
    output_path: Path,
    pass2_label: str,
    pass1_label: str,
    logy: bool,
    log_y_min: float,
    two_panel: bool,
    include_point_to_point: bool,
) -> None:
    if two_panel:
        fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.0))

        key = panel.key
        title = (
            rf"$x_B \in [{key.xb_min:.3g}, {key.xb_max:.3g}]$, "
            rf"$Q^2 \in [{key.q2_min:.3g}, {key.q2_max:.3g}]$ (GeV$^2$), "
            rf"$|t| \in [{key.t_min:.3g}, {key.t_max:.3g}]$ (GeV$^2$)"
        )

        fig.suptitle(title, fontsize=13)

        draw_cross_section_axis(
            ax=axes[0],
            panel=panel,
            curves=curves,
            pass2_label=pass2_label,
            pass1_label=pass1_label,
            logy=logy,
            log_y_min=log_y_min,
            include_title=False,
            include_point_to_point=include_point_to_point,
        )

        axes[0].set_title("Cross sections")

        draw_ratio_axis(
            ax=axes[1],
            panel=panel,
        )

        fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    else:
        fig, ax = plt.subplots(figsize=(8.4, 6.0))

        draw_cross_section_axis(
            ax=ax,
            panel=panel,
            curves=curves,
            pass2_label=pass2_label,
            pass1_label=pass1_label,
            logy=logy,
            log_y_min=log_y_min,
            include_title=True,
            include_point_to_point=include_point_to_point,
        )

        fig.tight_layout()
    #endif

    fig.savefig(output_path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Parallel worker.
# ---------------------------------------------------------------------------
def process_one_panel(job: WorkerJob) -> WorkerResult:
    t0 = time.time()

    cache_key: Optional[str] = None
    cache_entry: Optional[Dict[str, List[float]]] = None
    curves: Optional[ModelCurves] = None
    model_status = "models skipped"

    try:
        if not job.skip_models:
            cache_key = model_cache_key(job.panel.key, job.model_cfg)

            if job.cached_model_entry is not None:
                curves = make_model_curves_from_cache_entry(job.cached_model_entry)
                model_status = "model cache hit"
            else:
                try:
                    curves = compute_model_curves_without_cache(job.panel.key, job.model_cfg)
                    cache_entry = make_model_cache_entry(curves)
                    model_status = "model computed"
                except Exception as exc:
                    message = (
                        f"Model calculation failed for "
                        f"xB=[{job.panel.key.xb_min}, {job.panel.key.xb_max}], "
                        f"Q2=[{job.panel.key.q2_min}, {job.panel.key.q2_max}], "
                        f"|t|=[{job.panel.key.t_min}, {job.panel.key.t_max}]: {exc}"
                    )

                    if job.model_cfg.allow_missing_models:
                        curves = None
                        model_status = "model failed; plot written without models"
                    else:
                        raise RuntimeError(message) from exc
                    #endif
                #endtry
            #endif
        #endif

        output_path = Path(job.output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        draw_panel(
            panel=job.panel,
            curves=curves,
            output_path=output_path,
            pass2_label=job.pass2_label,
            pass1_label=job.pass1_label,
            logy=job.logy,
            log_y_min=job.log_y_min,
            two_panel=job.two_panel,
            include_point_to_point=job.model_cfg.include_point_to_point,  # injected dynamically below
        )

        ratio_points = pass2_over_pass1_ratio_points(job.panel)
        pass2_local_scales = [p.scale_syst_frac for p in job.panel.pass2 if math.isfinite(p.scale_syst_frac)]

        manifest_entry: Dict[str, object] = {
            "file": job.filename,
            "xB": [job.panel.key.xb_min, job.panel.key.xb_max],
            "Q2": [job.panel.key.q2_min, job.panel.key.q2_max],
            "t_abs": [job.panel.key.t_min, job.panel.key.t_max],
            "n_pass2_points": len(job.panel.pass2),
            "n_pass1_points": len(job.panel.pass1),
            "has_models": curves is not None,
            "two_panel": job.two_panel,
            "n_ratio_pass2_over_pass1": len(ratio_points),
            "pass2_local_s_comb_mean_in_panel": mean(pass2_local_scales),
            "pass2_local_s_comb_rms_in_panel": rms(pass2_local_scales),
            "pass2_local_s_comb_max_in_panel": max(pass2_local_scales) if pass2_local_scales else 0.0,
            "point_to_point_systematics_in_vertical_bars": job.model_cfg.include_point_to_point,
            "left_pass2_error_bars": "stat + 18% estimated point-to-point only" if job.model_cfg.include_point_to_point else "stat only",
            "left_pass2_scale_boxes": "point-by-point local period-spread s_comb(row)",
            "left_pass1_error_bars": "stat + provided syst only" if job.model_cfg.include_point_to_point else "stat only",
            "left_pass1_scale_boxes": "31% normalization",
            "ratio_error_bars": "pass2 stat and pass1 stat only",
            "ratio_scale_boxes": "local pass2 s_comb(row) scale and pass1 31% normalization",
            "model_status": model_status,
        }

        elapsed = time.time() - t0

        message = (
            f"[{job.index}/{job.total}] wrote {output_path} "
            f"({model_status}, pass2 points={len(job.panel.pass2)}, pass1 points={len(job.panel.pass1)}, "
            f"ratios: p2/p1={len(ratio_points)}, "
            f"elapsed={format_seconds(elapsed)})"
        )

        return WorkerResult(
            ok=True,
            index=job.index,
            filename=job.filename,
            output_path=str(output_path),
            manifest_entry=manifest_entry,
            cache_key=cache_key,
            cache_entry=cache_entry,
            message=message,
            model_status=model_status,
            elapsed_seconds=elapsed,
            error="",
        )
    except Exception:
        elapsed = time.time() - t0

        return WorkerResult(
            ok=False,
            index=job.index,
            filename=job.filename,
            output_path=job.output_path,
            manifest_entry={},
            cache_key=cache_key,
            cache_entry=cache_entry,
            message="",
            model_status=model_status,
            elapsed_seconds=elapsed,
            error=traceback.format_exc(),
        )
    #endtry


# Add the dynamic field to ModelConfig without changing worker pickle behavior.
setattr(ModelConfig, "include_point_to_point", True)


# ---------------------------------------------------------------------------
# External-tool path resolution.
# ---------------------------------------------------------------------------
def resolve_km15_cli(user_value: str) -> str:
    log("Path setup: resolving km15_cli.py.")

    candidates: List[Path] = []

    if user_value:
        candidates.append(Path(user_value).expanduser())
    #endif

    candidates.append(ANALYSIS_DIR / "km15_cli.py")
    candidates.append(Path.cwd() / "km15_cli.py")

    for cand in candidates:
        log(f"Path setup: checking KM15 CLI candidate: {cand}")

        if cand.exists() and cand.is_file():
            resolved = str(cand.resolve())
            log(f"Path setup: using KM15 CLI: {resolved}")
            return resolved
        #endif
    #endfor

    tried = "\n".join(f"  - {c}" for c in candidates)

    die(f"Could not find km15_cli.py. Tried:\n{tried}\nUse --km15-cli /full/path/to/km15_cli.py.")


def resolve_python_executable(user_value: str) -> str:
    log("Path setup: resolving Python executable for KM15.")

    text = str(user_value or "").strip()

    if text:
        expanded = Path(text).expanduser()

        log(f"Path setup: requested KM15 Python: {text}")

        if expanded.exists():
            resolved = str(expanded.resolve())
            log(f"Path setup: using requested KM15 Python: {resolved}")
            return resolved
        #endif

        found = shutil.which(text)

        if found:
            log(f"Path setup: using PATH-resolved KM15 Python: {found}")
            return found
        #endif

        warn(f"Requested KM15 Python executable does not exist: {text}; using current Python instead.")
    #endif

    if sys.executable and Path(sys.executable).exists():
        log(f"Path setup: using current Python executable for KM15: {sys.executable}")
        return sys.executable
    #endif

    found = shutil.which("python3")

    if found:
        log(f"Path setup: using PATH python3 for KM15: {found}")
        return found
    #endif

    log("Path setup: falling back to literal 'python3' for KM15.")

    return "python3"


def validate_dvcsgen_dir(path_text: str) -> str:
    log("Path setup: resolving dvcsgen directory.")

    path = Path(path_text).expanduser()
    exe = path / "dvcsgen"

    log(f"Path setup: requested dvcsgen dir: {path}")

    if not exe.exists():
        warn(
            f"dvcsgen executable was not found at {exe}. "
            f"Model calculation will fail unless --dvcsgen-dir is corrected "
            f"or --allow-missing-models is used."
        )
    else:
        log(f"Path setup: found dvcsgen executable: {exe}")
    #endif

    return str(path)


# ---------------------------------------------------------------------------
# Main driver helpers.
# ---------------------------------------------------------------------------
def build_panels(pass2: Dict[BinKey, List[DataPoint]], pass1: Dict[BinKey, List[DataPoint]]) -> List[PanelData]:
    log("Panel setup: merging pass-2 and pass-1 bin keys.")

    keys = sorted(set(pass2.keys()) | set(pass1.keys()), key=key_sort_tuple)
    panels: List[PanelData] = []

    n_both = 0
    n_pass2_only = 0
    n_pass1_only = 0

    for key in keys:
        has_pass2 = key in pass2
        has_pass1 = key in pass1

        if has_pass2 and has_pass1:
            n_both += 1
        elif has_pass2:
            n_pass2_only += 1
        else:
            n_pass1_only += 1
        #endif

        panels.append(PanelData(key=key, pass2=pass2.get(key, []), pass1=pass1.get(key, [])))
    #endfor

    log(f"Panel setup: total panels={len(panels)}, both pass-2/pass-1={n_both}, pass-2 only={n_pass2_only}, pass-1 only={n_pass1_only}.")

    return panels


def build_jobs(
    panels: Sequence[PanelData],
    output_dir: Path,
    pass2_label: str,
    pass1_label: str,
    logy: bool,
    log_y_min: float,
    two_panel: bool,
    skip_models: bool,
    model_cfg: ModelConfig,
    cache: Dict[str, Dict[str, List[float]]],
) -> Tuple[List[WorkerJob], int, int]:
    log("Job setup: building per-panel worker jobs.")

    jobs: List[WorkerJob] = []
    total = len(panels)
    n_cache_hit = 0
    n_cache_miss = 0

    for i, panel in enumerate(panels, start=1):
        filename = panel_filename(panel.key, i)
        output_path = output_dir / filename

        cached_model_entry: Optional[Dict[str, List[float]]] = None

        if not skip_models and model_cfg.use_cache:
            ckey = model_cache_key(panel.key, model_cfg)
            cached_model_entry = cache.get(ckey)

            if cached_model_entry is not None:
                n_cache_hit += 1
            else:
                n_cache_miss += 1
            #endif
        elif not skip_models:
            n_cache_miss += 1
        #endif

        jobs.append(
            WorkerJob(
                index=i,
                total=total,
                panel=panel,
                output_path=str(output_path),
                filename=filename,
                pass2_label=pass2_label,
                pass1_label=pass1_label,
                logy=logy,
                log_y_min=log_y_min,
                two_panel=two_panel,
                skip_models=skip_models,
                model_cfg=model_cfg,
                cached_model_entry=cached_model_entry,
            )
        )
    #endfor

    log(f"Job setup: built {len(jobs)} jobs.")
    log(f"Job setup: model cache hits={n_cache_hit}, model cache misses/to-compute={n_cache_miss}.")

    if not skip_models:
        n_external_model_calls = n_cache_miss * model_cfg.phi_dense * 2
        log(
            f"Job setup: estimated external model subprocess calls this run="
            f"{n_external_model_calls} ({n_cache_miss} uncached panels * {model_cfg.phi_dense} phi points * BH+KM15)."
        )
    #endif

    return jobs, n_cache_hit, n_cache_miss


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make one pass-2 vs pass-1 model-comparison phi plot per unique xB,Q2,|t| bin."
    )

    parser.add_argument("pass2_csv", type=Path, help="Pass-2 CSV. This must be the first positional input.")
    parser.add_argument("pass1_csv", type=Path, help="Pass-1 / Lee CSV. This must be the second positional input.")
    parser.add_argument("--output-dir", type=Path, default=Path("output/pass2_vs_pass1_model_comparison"))
    parser.add_argument("--pass2-xs-column", default=None, help="Override the pass-2 central cross-section column name.")
    parser.add_argument("--pass2-label", default="CLAS12 pass-2")
    parser.add_argument("--pass1-label", default="CLAS12 pass-1")
    parser.add_argument("--e-beam", type=float, default=10.604, help="Beam energy used for BH/KM15 curves (GeV).")
    parser.add_argument("--dvcsgen-dir", default=os.environ.get("DVCSGEN_PATH", DEFAULT_DVCSGEN_DIR))
    parser.add_argument(
        "--km15-cli",
        default=os.environ.get("KM15_CLI", DEFAULT_KM15_CLI),
        help="Path to km15_cli.py. Default is ../km15_cli.py relative to this script.",
    )
    parser.add_argument(
        "--py-km15",
        default=os.environ.get("PY_KM15", ""),
        help="Python executable for KM15. Default is the current Python executable.",
    )
    parser.add_argument(
        "--phi-dense",
        type=int,
        default=73,
        help="Number of phi points for theory curves. Default 73 gives 5-degree spacing. Use 181 for 2-degree spacing.",
    )
    parser.add_argument("--model-cache", type=Path, default=None, help="Optional JSON cache path for model curves.")
    parser.add_argument("--no-cache", action="store_true", help="Do not read existing cached model curves.")
    parser.add_argument("--allow-missing-models", action="store_true", help="Still make data plots if BH/KM15 commands fail.")
    parser.add_argument("--skip-models", action="store_true", help="Do not compute or draw BH/KM15 curves.")
    parser.add_argument("--linear-y", action="store_true", help="Use linear y scale instead of the default log y scale.")
    parser.add_argument(
        "--log-y-min",
        type=float,
        default=1.0e-5,
        help="Visible lower y-limit floor for log-y plots. Default is 1e-5.",
    )
    parser.add_argument("--two-panel", action="store_true", help="Make a 1x2 figure: cross sections on left, pass-2/pass-1 on right.")
    parser.add_argument(
        "--no-point-to-point-systematics",
        action="store_true",
        help=(
            "Remove point-to-point systematics from plotted vertical error bars. "
            "Pass-2 bars become stat only instead of stat⊕18%; pass-1 bars become stat only instead of stat⊕Lee syst. "
            "External pass-2 local scale boxes and pass-1 31% normalization boxes are still drawn."
        ),
    )
    parser.add_argument("--print-columns", action="store_true", help="Print detected CSV columns before plotting.")
    parser.add_argument(
        "--workers",
        type=int,
        default=5,
        help="Maximum number of parallel panel workers. Default is 5. Values above 5 are capped to 5. Use 1 for serial debugging.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1,
        help="Print aggregate progress every N completed jobs. Default is 1.",
    )
    parser.add_argument(
        "--quiet-workers",
        action="store_true",
        help="Suppress per-panel completion messages, but keep setup and final summary messages.",
    )
    parser.add_argument(
        "--verbose-worker-models",
        action="store_true",
        help="Very noisy: print model-subprocess progress from inside workers every ~25 phi points.",
    )

    return parser.parse_args()


def main() -> int:
    script_t0 = time.time()

    log("Starting pass-2 vs pass-1 model-comparison plotter.")
    log(f"Script path: {SCRIPT_PATH}")
    log(f"Current working directory: {Path.cwd()}")

    args = parse_args()

    include_point_to_point = not args.no_point_to_point_systematics

    log("Command-line configuration:")
    log(f"  pass2_csv                       = {args.pass2_csv}")
    log(f"  pass1_csv                       = {args.pass1_csv}")
    log(f"  output_dir                      = {args.output_dir}")
    log(f"  e_beam                          = {args.e_beam:g} GeV")
    log(f"  requested workers               = {args.workers}")
    log(f"  phi_dense                       = {args.phi_dense}")
    log(f"  skip_models                     = {args.skip_models}")
    log(f"  no_cache                        = {args.no_cache}")
    log(f"  allow_missing_models            = {args.allow_missing_models}")
    log(f"  linear_y                        = {args.linear_y}")
    log(f"  log_y_min                       = {args.log_y_min:g}")
    log(f"  two_panel                       = {args.two_panel}")
    log(f"  quiet_workers                   = {args.quiet_workers}")
    log(f"  progress_every                  = {args.progress_every}")
    log(f"  point-to-point systematics bars = {include_point_to_point}")
    log(f"  pass2 point-to-point            = {100.0 * PASS2_ESTIMATED_SYSTEMATIC_FRACTION:.1f}%")
    log("  pass2 scale boxes               = point-by-point local period-spread s_comb(row)")
    log(f"  pass1 normalization             = {100.0 * PASS1_NORMALIZATION_FRACTION:.1f}% drawn as external boxes")
    log(
        "  stat-only pull fit              = pass-1 normalization floated in "
        f"[{PULL_FIT_NORM_MIN:.2f}, {PULL_FIT_NORM_MAX:.2f}] using stats only"
    )
    log(f"  scale box half width            = {SCALE_BOX_HALF_WIDTH_DEG:g} deg")
    log(f"  pass1 phi offset                = {PASS1_PHI_OFFSET_DEG:+.2f} deg")
    log(f"  pass2 phi offset                = {PASS2_PHI_OFFSET_DEG:+.2f} deg")

    if args.two_panel and args.skip_models:
        warn("--two-panel was requested together with --skip-models. Left-panel BH/KM15 curves will be absent.")
    #endif

    if not args.pass2_csv.exists():
        die(f"Pass-2 CSV does not exist: {args.pass2_csv}")
    #endif

    if not args.pass1_csv.exists():
        die(f"Pass-1 CSV does not exist: {args.pass1_csv}")
    #endif

    if args.log_y_min <= 0.0:
        warn(f"Requested --log-y-min {args.log_y_min:g} is not positive; using 1e-5.")
        args.log_y_min = 1.0e-5
    #endif

    log("Output setup: creating output directory if needed.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log(f"Output setup: output directory ready: {args.output_dir.resolve()}")

    cache_path = args.model_cache or (args.output_dir / "model_curve_cache.json")
    log(f"Output setup: model cache path: {cache_path}")

    args.dvcsgen_dir = validate_dvcsgen_dir(args.dvcsgen_dir)
    args.km15_cli = resolve_km15_cli(args.km15_cli)
    args.py_km15 = resolve_python_executable(args.py_km15)

    log("Input setup: loading pass-2 data and computing point-by-point local scale systematics.")
    pass2, pass2_local_scale_summary = load_pass2_csv(
        path=args.pass2_csv,
        xs_column=args.pass2_xs_column,
        output_dir=args.output_dir,
        include_point_to_point=include_point_to_point,
        print_columns=args.print_columns,
    )

    log("Input setup: loading pass-1 data.")
    pass1 = load_pass1_csv(
        path=args.pass1_csv,
        include_point_to_point=include_point_to_point,
        print_columns=args.print_columns,
    )

    panels = build_panels(pass2, pass1)

    if not panels:
        die("No panels to write after loading pass-2 and pass-1 CSVs.")
    #endif

    log(f"Panel setup: unique (xB,Q2,|t|) panels to write: {len(panels)}")

    stat_only_summary = make_stat_only_pull_diagnostic(
        output_dir=args.output_dir,
        panels=panels,
    )

    n_workers = max(1, min(int(args.workers), 5, len(panels)))

    if args.workers > 5:
        warn(f"Requested --workers {args.workers}, capped to 5 for safety.")
    #endif

    log(f"Parallel setup: using {n_workers} worker(s).")

    model_cfg = ModelConfig(
        e_beam=args.e_beam,
        dvcsgen_dir=args.dvcsgen_dir,
        km15_cli=args.km15_cli,
        py_km15=args.py_km15,
        phi_dense=max(2, args.phi_dense),
        use_cache=not args.no_cache,
        allow_missing_models=args.allow_missing_models,
        verbose_worker_models=args.verbose_worker_models,
    )
    model_cfg.include_point_to_point = include_point_to_point

    cache: Dict[str, Dict[str, List[float]]] = {}

    if not args.skip_models:
        log("Model setup:")
        log(f"  E beam       = {args.e_beam:g} GeV")
        log(f"  dvcsgen dir  = {args.dvcsgen_dir}")
        log(f"  KM15 CLI     = {args.km15_cli}")
        log(f"  KM15 Python  = {args.py_km15}")
        log(f"  phi grid     = {model_cfg.phi_dense} points")
        log(f"  cache enabled= {model_cfg.use_cache}")

        cache = load_model_cache(cache_path) if model_cfg.use_cache else {}
    else:
        log("Model setup: --skip-models was provided; no BH/KM15 curves will be computed or drawn.")
    #endif

    jobs, n_cache_hit, n_cache_miss = build_jobs(
        panels=panels,
        output_dir=args.output_dir,
        pass2_label=args.pass2_label,
        pass1_label=args.pass1_label,
        logy=not args.linear_y,
        log_y_min=args.log_y_min,
        two_panel=args.two_panel,
        skip_models=args.skip_models,
        model_cfg=model_cfg,
        cache=cache,
    )

    manifest_by_index: Dict[int, Dict[str, object]] = {}
    failures: List[WorkerResult] = []
    n_new_cache_entries = 0
    n_model_computed = 0
    n_model_cache_hit_completed = 0
    n_model_failed_but_plotted = 0

    completed = 0
    progress_every = max(1, int(args.progress_every))
    processing_t0 = time.time()

    log("Processing: starting panel jobs.")

    if n_workers == 1:
        log("Processing: running in serial mode.")

        for job in jobs:
            result = process_one_panel(job)
            completed += 1

            if result.ok:
                if not args.quiet_workers:
                    log(result.message)
                #endif

                manifest_by_index[result.index] = result.manifest_entry

                if result.cache_key is not None and result.cache_entry is not None:
                    cache[result.cache_key] = result.cache_entry
                    n_new_cache_entries += 1
                #endif

                if result.model_status == "model computed":
                    n_model_computed += 1
                elif result.model_status == "model cache hit":
                    n_model_cache_hit_completed += 1
                elif result.model_status == "model failed; plot written without models":
                    n_model_failed_but_plotted += 1
                #endif
            else:
                failures.append(result)
                print(result.error, file=sys.stderr)
            #endif

            if completed % progress_every == 0 or completed == len(jobs):
                elapsed = time.time() - processing_t0
                rate = completed / elapsed if elapsed > 0.0 else 0.0
                remaining = len(jobs) - completed
                eta = remaining / rate if rate > 0.0 else 0.0

                log(
                    f"Progress: completed {completed}/{len(jobs)} panels "
                    f"({100.0 * completed / len(jobs):.1f}%), failures={len(failures)}, "
                    f"elapsed={format_seconds(elapsed)}, ETA≈{format_seconds(eta)}."
                )
            #endif
        #endfor
    else:
        log("Processing: running in parallel mode.")

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_index = {}

            log("Processing: submitting jobs to worker pool.")

            for job in jobs:
                future = executor.submit(process_one_panel, job)
                future_to_index[future] = job.index
            #endfor

            log(f"Processing: submitted {len(future_to_index)} jobs.")

            for future in as_completed(future_to_index):
                result = future.result()
                completed += 1

                if result.ok:
                    if not args.quiet_workers:
                        log(result.message)
                    #endif

                    manifest_by_index[result.index] = result.manifest_entry

                    if result.cache_key is not None and result.cache_entry is not None:
                        cache[result.cache_key] = result.cache_entry
                        n_new_cache_entries += 1
                    #endif

                    if result.model_status == "model computed":
                        n_model_computed += 1
                    elif result.model_status == "model cache hit":
                        n_model_cache_hit_completed += 1
                    elif result.model_status == "model failed; plot written without models":
                        n_model_failed_but_plotted += 1
                    #endif
                else:
                    failures.append(result)
                    print(result.error, file=sys.stderr)
                #endif

                if completed % progress_every == 0 or completed == len(jobs):
                    elapsed = time.time() - processing_t0
                    rate = completed / elapsed if elapsed > 0.0 else 0.0
                    remaining = len(jobs) - completed
                    eta = remaining / rate if rate > 0.0 else 0.0

                    log(
                        f"Progress: completed {completed}/{len(jobs)} panels "
                        f"({100.0 * completed / len(jobs):.1f}%), failures={len(failures)}, "
                        f"elapsed={format_seconds(elapsed)}, ETA≈{format_seconds(eta)}."
                    )
                #endif
            #endfor
        #endwith
    #endif

    processing_dt = time.time() - processing_t0

    log(f"Processing: panel jobs finished in {format_seconds(processing_dt)}.")
    log(f"Processing summary: successful panels={len(manifest_by_index)}, failed panels={len(failures)}.")
    log(f"Processing summary: model cache hits completed={n_model_cache_hit_completed}.")
    log(f"Processing summary: model curves newly computed={n_model_computed}.")
    log(f"Processing summary: model failures allowed and plotted without models={n_model_failed_but_plotted}.")
    log(f"Processing summary: new cache entries returned by workers={n_new_cache_entries}.")

    if failures:
        n_fail = len(failures)
        first = failures[0]

        die(
            f"{n_fail} panel job(s) failed. First failed panel index={first.index}, "
            f"file={first.filename}. See traceback above."
        )
    #endif

    manifest: List[Dict[str, object]] = []

    log("Manifest: sorting manifest entries by panel index.")

    for index in sorted(manifest_by_index):
        entry = manifest_by_index[index]
        entry["stat_only_pull_best_pass1_norm"] = stat_only_summary.best_pass1_norm
        entry["stat_only_pull_best_pass1_norm_percent"] = 100.0 * (stat_only_summary.best_pass1_norm - 1.0)
        entry["stat_only_pull_chi2_ndf_nominal"] = stat_only_summary.chi2_ndf_nominal
        entry["stat_only_pull_chi2_ndf_best"] = stat_only_summary.chi2_ndf_best
        entry["point_to_point_systematics_in_vertical_bars"] = include_point_to_point
        manifest.append(entry)
    #endfor

    if not args.skip_models:
        save_model_cache(cache_path, cache)
        log(f"Cache summary: loaded hits before run={n_cache_hit}, misses before run={n_cache_miss}, new entries this run={n_new_cache_entries}.")
    #endif

    manifest_path = args.output_dir / "manifest.json"

    log(f"Manifest: writing manifest to {manifest_path}")

    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
    #endwith

    log(f"Manifest: wrote {len(manifest)} entries.")

    total_dt = time.time() - script_t0

    log("Final summary:")
    log(f"  output directory              = {args.output_dir.resolve()}")
    log(f"  total panels written          = {len(manifest)}")
    log(f"  workers used                  = {n_workers}")
    log(f"  two-panel mode                = {args.two_panel}")
    log(f"  model phi points              = {model_cfg.phi_dense}")
    log(f"  point-to-point syst bars      = {include_point_to_point}")
    log(f"  pass2 local scale rows        = {pass2_local_scale_summary.n_rows_with_local_scale}")
    log(f"  pass2 local s_comb mean       = {pass2_local_scale_summary.mean_s_comb:.10g} ({100.0 * pass2_local_scale_summary.mean_s_comb:.6g}%)")
    log(f"  pass2 local s_comb median     = {pass2_local_scale_summary.median_s_comb:.10g} ({100.0 * pass2_local_scale_summary.median_s_comb:.6g}%)")
    log(f"  pass2 local s_comb rms        = {pass2_local_scale_summary.rms_s_comb:.10g} ({100.0 * pass2_local_scale_summary.rms_s_comb:.6g}%)")
    log(f"  pass2 local s_comb max        = {pass2_local_scale_summary.max_s_comb:.10g} ({100.0 * pass2_local_scale_summary.max_s_comb:.6g}%)")

    if include_point_to_point:
        log("  left pass2 error bars         = stat ⊕ 18% point-to-point")
        log("  left pass1 error bars         = stat ⊕ provided syst")
    else:
        log("  left pass2 error bars         = stat only")
        log("  left pass1 error bars         = stat only")
    #endif

    log("  left pass2 scale boxes        = local point-by-point period spread")
    log("  left pass1 scale boxes        = 31% normalization")
    log("  right panel ratio             = pass-2 / pass-1 only")
    log("  right ratio error bars        = pass2 stat ⊕ pass1 stat")
    log("  right ratio scale boxes       = local pass2 scale ⊕ pass1 31% norm")
    log(f"  stat-only matched points      = {stat_only_summary.n_points}")
    log(
        "  stat-only best pass1 norm     = "
        f"{stat_only_summary.best_pass1_norm:.10g} "
        f"({100.0 * (stat_only_summary.best_pass1_norm - 1.0):+.6g}%)"
    )
    log(
        "  stat-only chi2/ndf            = "
        f"M=1 {stat_only_summary.chi2_ndf_nominal:.10g}, "
        f"best {stat_only_summary.chi2_ndf_best:.10g}"
    )
    log(
        "  stat-only |pull|<=1           = "
        f"M=1 {stat_only_summary.pct_within_1sigma_nominal:.4g}%, "
        f"best {stat_only_summary.pct_within_1sigma_best:.4g}%"
    )
    log(
        "  stat-only |pull|<=3           = "
        f"M=1 {stat_only_summary.pct_within_3sigma_nominal:.4g}%, "
        f"best {stat_only_summary.pct_within_3sigma_best:.4g}%"
    )
    log(f"  log-y visible floor           = {args.log_y_min:g}")
    log(f"  cache entries after run       = {len(cache) if not args.skip_models else 0}")
    log(f"  total elapsed                 = {format_seconds(total_dt)}")
    log("Done.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[pass2-vs-pass1][{timestamp()}][FATAL] {exc}", file=sys.stderr)
        raise SystemExit(1)
    #endtry