#!/usr/bin/env python3
"""
plot_pass2_vs_pass1_model_comparison.py
--------------------------------------------------------------------------
Focused comparison plotter for CLAS12 pass-2 vs pass-1 DVCS cross sections.

This script intentionally does only the comparison plotting:
  - pass-2 and pass-1 cross sections vs phi,
  - optional BH and KM15 curves,
  - optional pass-2/pass-1 ratio panel,
  - pass-2 point-by-point local s_comb scale boxes,
  - pass-1 31% normalization boxes,
  - one global 1-D bin-number comparison plot.

Diagnostics such as s_comb histograms, pull distributions, outlier tables,
and kinematic maps belong in diagnose_pass2_pass1_consistency.py.

Usage:
    python3 plot_pass2_vs_pass1_model_comparison.py \
        pass2_dvcs.csv pass1_lee.csv \
        --output-dir output/pass2_vs_pass1_model_comparison

Optional two-panel mode:
    python3 plot_pass2_vs_pass1_model_comparison.py \
        pass2_dvcs.csv pass1_lee.csv \
        --output-dir output/pass2_vs_pass1_model_comparison \
        --two-panel

Optional point-to-point-systematic override:
    python3 plot_pass2_vs_pass1_model_comparison.py \
        pass2_dvcs.csv pass1_lee.csv \
        --output-dir output/pass2_vs_pass1_model_comparison \
        --two-panel \
        --no-point-to-point-systematics

New global bin-number comparison:
  The script also writes:

      pass1_pass2_bin_number_ratio_comparison.png
      pass1_pass2_bin_number_ratio_comparison.csv

  For each matched pass-2/pass-1 point, ordered by (xB, Q2, |t|, phi):

      x = sequential matched-point number, starting at 1

      pass-1 y = 1

      pass-1 y error =
          sqrt( (pass1_stat / pass1_xs)^2 + 0.31^2 )

      pass-2 y =
          pass2_xs / pass1_xs

      pass-2 y error =
          (pass2_xs / pass1_xs) *
          sqrt( (pass2_stat / pass2_xs)^2 + s_comb(row)^2 )

  This global plot is intended as a compact visual check of all matched
  pass-2/pass-1 points at once. The pass-1 denominator uncertainty is not
  folded into the pass-2 ratio error here, because pass-1 is shown separately
  as the normalized reference band/points.

Uncertainty prescription:
  * Pass-2 central values are read from:

        normed cross sections, ep->epg, exp, 10.6 GeV, unpol

    or, if needed, from:

        cross sections, ep->epg, exp, 10.6 GeV, unpol

    The script deliberately refuses to use a column containing "combination sys"
    as the central cross-section source.

  * Pass-2 local scale systematic is computed row-by-row from:

        Fa18 Inb, Fa18 Out, Sp18 Inb, Sp18 Out

    using:

        ref(row)    = stat-weighted mean of valid period tuples
        r_i(row)    = sigma_i / ref
        s_obs(row)  = RMS_i(r_i - 1)
        s_stat(row) = RMS_i(delta_i / ref)
        s_comb(row) = sqrt(max(0, s_obs(row)^2 - s_stat(row)^2))

    The pass-2 scale box for that point is:

        s_comb(row) * sigma_pass2(row)

  * Cross-section panel:
      pass-2 vertical errors:
          stat ⊕ 18% estimated point-to-point systematic

      pass-2 external boxes:
          local s_comb(row) * cross section

      pass-1 vertical errors:
          stat ⊕ provided Lee systematic

      pass-1 external boxes:
          31% normalization * cross section

      If --no-point-to-point-systematics is used:
          pass-2 vertical errors = stat only
          pass-1 vertical errors = stat only

      The external boxes are still drawn.

  * Ratio panel:
      pass-2/pass-1 vertical errors:
          propagated statistical uncertainty only

      pass-2/pass-1 external boxes:
          propagated scale uncertainty from local pass-2 s_comb(row)
          and pass-1 31% normalization.

  * One PNG is written for every unique (xB, Q2, |t|) bin.
  * Default log-y visible floor is 1e-4.
  * PNG output only.
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
from typing import Dict, List, Optional, Sequence, Tuple

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
# Analysis constants.
# ---------------------------------------------------------------------------
PASS2_ESTIMATED_SYSTEMATIC_FRACTION = 0.18
PASS1_NORMALIZATION_FRACTION = 0.31

PASS1_PHI_OFFSET_DEG = -1.25
PASS2_PHI_OFFSET_DEG = +1.25
PASS1_PASS2_RATIO_MATCH_TOLERANCE_DEG = 2.0

SCALE_BOX_HALF_WIDTH_DEG = 3.0
SCALE_BOX_ALPHA = 0.18

DEFAULT_LOG_Y_MIN = 1.0e-4


# ---------------------------------------------------------------------------
# Pass-1 / Lee CSV columns.
# ---------------------------------------------------------------------------
PASS1_XS_COL = "cross sections, ep->epg, exp"
PASS1_STAT_COL = "cross sections, ep->epg, exp, stat. unc."
PASS1_SYST_UP_COL = "cross sections, ep->epg, exp, syst. unc. (up)"
PASS1_SYST_DN_COL = "cross sections, ep->epg, exp, syst. unc. (down)"


# ---------------------------------------------------------------------------
# Pass-2 columns.
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


@dataclass
class GlobalBinRatioPoint:
    bin_number: int
    key: BinKey
    phi_pass2: float
    phi_pass1: float
    pass1_y: float
    pass1_yerr: float
    pass2_y: float
    pass2_yerr: float
    pass1_xs: float
    pass1_stat: float
    pass2_xs: float
    pass2_stat: float
    pass2_s_comb: float
    pass2_n_good_periods: int


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
# Generic helpers.
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
    cmap = canonical_header_map(fieldnames)
    missing = [col for col in required if col.strip().lower() not in cmap]

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

        if not match:
            return default
        #endif

        try:
            out = float(match.group(0))

            if math.isfinite(out):
                return out
            #endif
        except ValueError:
            return default
        #endtry
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


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return ""
    #endif

    return f"{value:.12g}"


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
# Tuple parsing and local s_comb.
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

    return TupleValue(ok=True, value=value, stat=stat)


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

    out.n_valid_periods = len(ratios)
    out.ref_value = ref.value
    out.ref_stat = ref.stat
    out.s_obs = s_obs
    out.s_stat_exp = s_stat_exp
    out.s_comb = s_comb
    out.min_ratio = min(ratios)
    out.max_ratio = max(ratios)
    out.half_width = 0.5 * (out.max_ratio - out.min_ratio)

    return out


# ---------------------------------------------------------------------------
# Uncertainty prescriptions.
# ---------------------------------------------------------------------------
def pass2_vertical_uncertainty(xs: float, stat: float, no_point_to_point_systematics: bool) -> float:
    if no_point_to_point_systematics:
        return abs(stat)
    #endif

    return math.sqrt(stat * stat + (PASS2_ESTIMATED_SYSTEMATIC_FRACTION * xs) ** 2)


def pass1_vertical_uncertainty(xs: float, stat: float, syst: float, no_point_to_point_systematics: bool) -> float:
    if no_point_to_point_systematics:
        return abs(stat)
    #endif

    return math.sqrt(stat * stat + syst * syst)


def pass1_scale_uncertainty(xs: float) -> float:
    return abs(PASS1_NORMALIZATION_FRACTION * xs)


def ratio_scale_uncertainty(ratio: float, pass2_local_scale_frac: float) -> float:
    rel = math.sqrt(pass2_local_scale_frac ** 2 + PASS1_NORMALIZATION_FRACTION ** 2)

    return abs(ratio * rel)


def global_pass1_normalized_yerr(pass1_point: DataPoint) -> float:
    if not finite_positive(pass1_point.xs):
        return float("nan")
    #endif

    rel_stat = abs(pass1_point.stat / pass1_point.xs)

    return math.sqrt(rel_stat * rel_stat + PASS1_NORMALIZATION_FRACTION * PASS1_NORMALIZATION_FRACTION)


def global_pass2_ratio_yerr(pass2_point: DataPoint, pass1_point: DataPoint) -> float:
    if not finite_positive(pass2_point.xs) or not finite_positive(pass1_point.xs):
        return float("nan")
    #endif

    ratio = pass2_point.xs / pass1_point.xs
    rel_stat = abs(pass2_point.stat / pass2_point.xs)
    rel_scale = abs(pass2_point.scale_syst_frac)

    return abs(ratio * math.sqrt(rel_stat * rel_stat + rel_scale * rel_scale))


# ---------------------------------------------------------------------------
# Input loaders.
# ---------------------------------------------------------------------------
def write_pass2_local_scale_used_csv(output_dir: Path, rows: Sequence[Tuple[BinKey, float, float, LocalScaleResult]]) -> None:
    path = output_dir / "pass2_local_s_comb_used_by_comparison.csv"

    log(f"Pass-2: writing local s_comb values used by plots to {path}")

    with path.open("w", newline="") as handle:
        handle.write(
            "xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phi,pass2_xs,"
            "n_valid_periods,ref_value,ref_stat,s_obs,s_stat,s_comb,"
            "s_obs_percent,s_stat_percent,s_comb_percent,min_ratio,max_ratio,half_width\n"
        )

        for key, phi, xs, local in rows:
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
                f"{format_float(100.0 * local.s_obs)},"
                f"{format_float(100.0 * local.s_stat_exp)},"
                f"{format_float(100.0 * local.s_comb)},"
                f"{format_float(local.min_ratio)},"
                f"{format_float(local.max_ratio)},"
                f"{format_float(local.half_width)}\n"
            )
        #endfor
    #endwith


def load_pass2_csv(
    path: Path,
    output_dir: Path,
    xs_column: Optional[str],
    no_point_to_point_systematics: bool,
    print_columns: bool = False,
) -> Dict[BinKey, List[DataPoint]]:
    t0 = time.time()

    log(f"Pass-2: opening CSV: {path}")

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            die(f"Pass-2 CSV appears empty: {path}")
        #endif

        fieldnames = reader.fieldnames

        if print_columns:
            log("Pass-2 columns:\n" + "\n".join(f"  [{i}] {name}" for i, name in enumerate(fieldnames)))
        #endif

        cols = detect_common_columns(fieldnames, "Pass-2")
        xs_col = xs_column or find_column(fieldnames, PASS2_CENTRAL_XS_CANDIDATES)

        if "combination sys" in xs_col.lower():
            die(
                "Selected pass-2 central-value column contains 'combination sys'. "
                "Use the central cross-section column instead. "
                f"Selected: {xs_col}"
            )
        #endif

        required_period_cols = [pass2_period_cross_section_column(period) for period in PASS2_PERIODS_10P6_UNPOL]
        require_columns(fieldnames, required_period_cols, "pass-2 local s_comb calculation")

        period_cols = {
            period: find_column(fieldnames, [pass2_period_cross_section_column(period)])
            for period in PASS2_PERIODS_10P6_UNPOL
        }

        log(f"Pass-2: central cross-section column -> {xs_col}")

        for period in PASS2_PERIODS_10P6_UNPOL:
            log(f"Pass-2: local scale input {period:10s} -> {period_cols[period]}")
        #endfor

        if no_point_to_point_systematics:
            log("Pass-2: vertical errors -> stat only.")
        else:
            log("Pass-2: vertical errors -> stat ⊕ 18% estimated point-to-point.")
        #endif

        log("Pass-2: external boxes -> local point-by-point s_comb(row) * xs.")

        out: Dict[BinKey, List[DataPoint]] = {}
        local_rows: List[Tuple[BinKey, float, float, LocalScaleResult]] = []

        total_rows = 0
        kept = 0
        skipped_invalid = 0
        skipped_bad_xs = 0
        n_zero_local = 0
        n_nonzero_local = 0

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

            period_values = [
                parse_tuple_cell(row.get(period_cols[period], ""))
                for period in PASS2_PERIODS_10P6_UNPOL
            ]

            local = compute_local_pass2_scale(period_values)

            xs = tuple_value.value
            stat = tuple_value.stat
            point_err = pass2_vertical_uncertainty(
                xs=xs,
                stat=stat,
                no_point_to_point_systematics=no_point_to_point_systematics,
            )

            scale_abs = abs(local.s_comb * xs)

            if local.s_comb > 0.0:
                n_nonzero_local += 1
            else:
                n_zero_local += 1
            #endif

            point = DataPoint(
                phi=phi,
                xs=xs,
                stat=stat,
                err_low=point_err,
                err_high=point_err,
                source="pass2",
                scale_syst_abs=scale_abs,
                scale_syst_frac=local.s_comb,
                local_scale_n_periods=local.n_valid_periods,
                local_scale_s_obs=local.s_obs,
                local_scale_s_stat=local.s_stat_exp,
                local_scale_half_width=local.half_width,
            )

            out.setdefault(key, []).append(point)
            local_rows.append((key, phi, xs, local))
            kept += 1
        #endfor
    #endwith

    for points in out.values():
        points.sort(key=lambda p: p.phi)
    #endfor

    write_pass2_local_scale_used_csv(output_dir, local_rows)

    dt = time.time() - t0

    log(f"Pass-2: finished reading in {format_seconds(dt)}.")
    log(
        f"Pass-2: total rows={total_rows}, kept={kept}, "
        f"skipped invalid-bin={skipped_invalid}, skipped bad/nonpositive xs={skipped_bad_xs}."
    )
    log(f"Pass-2: local s_comb nonzero points={n_nonzero_local}, zero points={n_zero_local}.")
    log(f"Pass-2: unique (xB,Q2,|t|) bins with data={len(out)}.")

    return out


def load_pass1_csv(
    path: Path,
    no_point_to_point_systematics: bool,
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

        if no_point_to_point_systematics:
            log("Pass-1: vertical errors -> stat only.")
        else:
            log("Pass-1: vertical errors -> stat ⊕ provided Lee systematic.")
        #endif

        log("Pass-1: external boxes -> 31% normalization * xs.")

        out: Dict[BinKey, List[DataPoint]] = {}

        total_rows = 0
        kept = 0
        skipped_invalid = 0
        skipped_bad_xs = 0

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

            err_high = pass1_vertical_uncertainty(
                xs=xs,
                stat=stat,
                syst=syst_up,
                no_point_to_point_systematics=no_point_to_point_systematics,
            )
            err_low = pass1_vertical_uncertainty(
                xs=xs,
                stat=stat,
                syst=syst_dn,
                no_point_to_point_systematics=no_point_to_point_systematics,
            )

            if err_low <= 0.0:
                err_low = err_high
            #endif

            if err_high <= 0.0:
                err_high = err_low
            #endif

            point = DataPoint(
                phi=phi,
                xs=xs,
                stat=stat,
                err_low=err_low,
                err_high=err_high,
                source="pass1",
                scale_syst_abs=pass1_scale_uncertainty(xs),
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
    log(
        f"Pass-1: total rows={total_rows}, kept={kept}, "
        f"skipped invalid-bin={skipped_invalid}, skipped bad/nonpositive xs={skipped_bad_xs}."
    )
    log(f"Pass-1: unique (xB,Q2,|t|) bins with data={len(out)}.")

    return out


# ---------------------------------------------------------------------------
# External model helpers.
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
    if not path.exists():
        log("Cache: no existing model cache found.")
        return {}
    #endif

    try:
        with path.open("r") as handle:
            raw = json.load(handle)
        #endwith

        if not isinstance(raw, dict):
            return {}
        #endif

        good: Dict[str, Dict[str, List[float]]] = {}

        for key, entry in raw.items():
            if not isinstance(entry, dict):
                continue
            #endif

            if "phi" not in entry or "bh" not in entry or "km15" not in entry:
                continue
            #endif

            good[key] = {
                "phi": list(entry["phi"]),
                "bh": list(entry["bh"]),
                "km15": list(entry["km15"]),
            }
        #endfor

        log(f"Cache: loaded {len(good)} model entries from {path}")
        return good
    except Exception as exc:
        warn(f"Cache: could not read {path}: {exc}")
        return {}
    #endtry


def save_model_cache(path: Path, cache: Dict[str, Dict[str, List[float]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")

    with tmp_path.open("w") as handle:
        json.dump(cache, handle, indent=2, sort_keys=True)
    #endwith

    tmp_path.replace(path)
    log(f"Cache: wrote {len(cache)} model entries to {path}")


def make_model_curves_from_cache_entry(entry: Dict[str, List[float]]) -> ModelCurves:
    return ModelCurves(phi=list(entry["phi"]), bh=list(entry["bh"]), km15=list(entry["km15"]))


def make_model_cache_entry(curves: ModelCurves) -> Dict[str, List[float]]:
    return {"phi": curves.phi, "bh": curves.bh, "km15": curves.km15}


def compute_model_curves_without_cache(key: BinKey, cfg: ModelConfig) -> ModelCurves:
    xb = 0.5 * (key.xb_min + key.xb_max)
    q2 = 0.5 * (key.q2_min + key.q2_max)
    t_abs = 0.5 * (key.t_min + key.t_max)

    phi_grid = [360.0 * i / (cfg.phi_dense - 1) for i in range(cfg.phi_dense)] if cfg.phi_dense > 1 else [0.0]

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
# Ratio helpers.
# ---------------------------------------------------------------------------
def find_nearest_pass1_point(phi: float, pass1_points: Sequence[DataPoint]) -> Optional[DataPoint]:
    best_point: Optional[DataPoint] = None
    best_delta = float("inf")

    for point in pass1_points:
        delta = abs(point.phi - phi)

        if delta < best_delta:
            best_delta = delta
            best_point = point
        #endif
    #endfor

    if best_point is None or best_delta > PASS1_PASS2_RATIO_MATCH_TOLERANCE_DEG:
        return None
    #endif

    return best_point


def pass2_over_pass1_ratio_points(panel: PanelData) -> List[RatioPoint]:
    out: List[RatioPoint] = []

    for p2 in panel.pass2:
        p1 = find_nearest_pass1_point(p2.phi, panel.pass1)

        if p1 is None:
            continue
        #endif

        if not finite_positive(p2.xs) or not finite_positive(p1.xs):
            continue
        #endif

        ratio = p2.xs / p1.xs

        rel_p2_stat = abs(p2.stat / p2.xs)
        rel_p1_stat = abs(p1.stat / p1.xs)
        stat_err = abs(ratio * math.sqrt(rel_p2_stat ** 2 + rel_p1_stat ** 2))

        scale_err = ratio_scale_uncertainty(ratio=ratio, pass2_local_scale_frac=p2.scale_syst_frac)

        out.append(
            RatioPoint(
                phi=p2.phi,
                ratio=ratio,
                stat_err_low=stat_err,
                stat_err_high=stat_err,
                scale_err=scale_err,
                pass2_scale_frac=p2.scale_syst_frac,
            )
        )
    #endfor

    return out


def ratio_axis_limits(ratio_sets: Sequence[Sequence[RatioPoint]]) -> Tuple[float, float]:
    vals: List[float] = []

    for ratio_points in ratio_sets:
        for p in ratio_points:
            extent = max(p.stat_err_low, p.stat_err_high, p.scale_err)
            vals.append(p.ratio - extent)
            vals.append(p.ratio + extent)
        #endfor
    #endfor

    vals = [v for v in vals if math.isfinite(v)]

    if not vals:
        return 0.0, 2.0
    #endif

    ymin = min(min(vals), 1.0)
    ymax = max(max(vals), 1.0)
    span = ymax - ymin

    if span <= 0.0:
        return ymin - 0.5, ymax + 0.5
    #endif

    pad = 0.15 * span

    return ymin - pad, ymax + pad


# ---------------------------------------------------------------------------
# Global bin-number comparison.
# ---------------------------------------------------------------------------
def build_global_bin_ratio_points(panels: Sequence[PanelData]) -> List[GlobalBinRatioPoint]:
    out: List[GlobalBinRatioPoint] = []
    bin_number = 0

    sorted_panels = sorted(panels, key=lambda p: key_sort_tuple(p.key))

    for panel in sorted_panels:
        pass2_points = sorted(panel.pass2, key=lambda p: p.phi)

        for p2 in pass2_points:
            p1 = find_nearest_pass1_point(p2.phi, panel.pass1)

            if p1 is None:
                continue
            #endif

            if not finite_positive(p1.xs) or not finite_positive(p2.xs):
                continue
            #endif

            pass1_y = 1.0
            pass1_yerr = global_pass1_normalized_yerr(p1)
            pass2_y = p2.xs / p1.xs
            pass2_yerr = global_pass2_ratio_yerr(p2, p1)

            if not math.isfinite(pass1_yerr) or not math.isfinite(pass2_yerr):
                continue
            #endif

            bin_number += 1

            out.append(
                GlobalBinRatioPoint(
                    bin_number=bin_number,
                    key=panel.key,
                    phi_pass2=p2.phi,
                    phi_pass1=p1.phi,
                    pass1_y=pass1_y,
                    pass1_yerr=pass1_yerr,
                    pass2_y=pass2_y,
                    pass2_yerr=pass2_yerr,
                    pass1_xs=p1.xs,
                    pass1_stat=p1.stat,
                    pass2_xs=p2.xs,
                    pass2_stat=p2.stat,
                    pass2_s_comb=p2.scale_syst_frac,
                    pass2_n_good_periods=p2.local_scale_n_periods,
                )
            )
        #endfor
    #endfor

    return out


def write_global_bin_ratio_csv(output_dir: Path, points: Sequence[GlobalBinRatioPoint]) -> None:
    path = output_dir / "pass1_pass2_bin_number_ratio_comparison.csv"

    log(f"Global bin-number plot: writing CSV to {path}")

    with path.open("w", newline="") as handle:
        handle.write(
            "bin_number,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,"
            "phi_pass2,phi_pass1,"
            "pass1_normalized_y,pass1_yerr_stat_plus_31percent_norm,"
            "pass2_over_pass1,pass2_yerr_stat_plus_s_comb,"
            "pass1_xs,pass1_stat,pass1_norm_frac,"
            "pass2_xs,pass2_stat,pass2_s_comb,pass2_n_good_periods\n"
        )

        for p in points:
            key = p.key

            handle.write(
                f"{p.bin_number},"
                f"{format_float(key.xb_min)},"
                f"{format_float(key.xb_max)},"
                f"{format_float(key.q2_min)},"
                f"{format_float(key.q2_max)},"
                f"{format_float(key.t_min)},"
                f"{format_float(key.t_max)},"
                f"{format_float(p.phi_pass2)},"
                f"{format_float(p.phi_pass1)},"
                f"{format_float(p.pass1_y)},"
                f"{format_float(p.pass1_yerr)},"
                f"{format_float(p.pass2_y)},"
                f"{format_float(p.pass2_yerr)},"
                f"{format_float(p.pass1_xs)},"
                f"{format_float(p.pass1_stat)},"
                f"{format_float(PASS1_NORMALIZATION_FRACTION)},"
                f"{format_float(p.pass2_xs)},"
                f"{format_float(p.pass2_stat)},"
                f"{format_float(p.pass2_s_comb)},"
                f"{p.pass2_n_good_periods}\n"
            )
        #endfor
    #endwith


def global_ratio_y_limits(
    points: Sequence[GlobalBinRatioPoint],
    user_ymin: Optional[float],
    user_ymax: Optional[float],
) -> Tuple[float, float]:
    if user_ymin is not None and user_ymax is not None:
        if user_ymax > user_ymin:
            return user_ymin, user_ymax
        #endif
    #endif

    vals: List[float] = []

    for p in points:
        vals.append(p.pass1_y - p.pass1_yerr)
        vals.append(p.pass1_y + p.pass1_yerr)
        vals.append(p.pass2_y - p.pass2_yerr)
        vals.append(p.pass2_y + p.pass2_yerr)
        vals.append(1.0)
    #endfor

    vals = [v for v in vals if math.isfinite(v)]

    if not vals:
        ymin = 0.0
        ymax = 2.0
    else:
        ymin = min(vals)
        ymax = max(vals)

        span = ymax - ymin

        if span <= 0.0:
            ymin -= 0.5
            ymax += 0.5
        else:
            pad = 0.10 * span
            ymin -= pad
            ymax += pad
        #endif
    #endif

    if user_ymin is not None:
        ymin = user_ymin
    #endif

    if user_ymax is not None:
        ymax = user_ymax
    #endif

    if ymax <= ymin:
        ymax = ymin + 1.0
    #endif

    return ymin, ymax


def draw_global_bin_ratio_plot(
    output_dir: Path,
    panels: Sequence[PanelData],
    pass2_label: str,
    pass1_label: str,
    user_ymin: Optional[float],
    user_ymax: Optional[float],
) -> None:
    points = build_global_bin_ratio_points(panels)

    if not points:
        warn("Global bin-number plot: no matched pass-2/pass-1 points; skipping.")
        return
    #endif

    write_global_bin_ratio_csv(output_dir, points)

    path = output_dir / "pass1_pass2_bin_number_ratio_comparison.png"

    log(f"Global bin-number plot: writing PNG to {path}")

    x_pass1 = [p.bin_number - 0.12 for p in points]
    x_pass2 = [p.bin_number + 0.12 for p in points]

    y_pass1 = [p.pass1_y for p in points]
    yerr_pass1 = [p.pass1_yerr for p in points]

    y_pass2 = [p.pass2_y for p in points]
    yerr_pass2 = [p.pass2_yerr for p in points]

    fig_width = max(12.0, min(28.0, 7.0 + 0.006 * len(points)))
    fig, ax = plt.subplots(figsize=(fig_width, 6.8))

    ax.errorbar(
        x_pass1,
        y_pass1,
        yerr=yerr_pass1,
        fmt="s",
        markersize=2.8,
        capsize=1.0,
        linewidth=0.8,
        linestyle="None",
        label=f"{pass1_label}: normalized to 1, stat ⊕ 31% norm",
        zorder=4,
    )

    ax.errorbar(
        x_pass2,
        y_pass2,
        yerr=yerr_pass2,
        fmt="o",
        markersize=2.8,
        capsize=1.0,
        linewidth=0.8,
        linestyle="None",
        label=f"{pass2_label}: pass-2/pass-1, stat ⊕ local $s_{{comb}}$",
        zorder=5,
    )

    ax.axhline(1.0, linewidth=1.2, linestyle="--", zorder=2)

    ymin, ymax = global_ratio_y_limits(points, user_ymin=user_ymin, user_ymax=user_ymax)

    ax.set_xlim(0.0, float(len(points)) + 1.0)
    ax.set_ylim(ymin, ymax)

    ax.set_xlabel("Matched point number")
    ax.set_ylabel("Normalized cross section ratio")
    ax.set_title("Global pass-2 vs pass-1 comparison by matched point number")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best", fontsize=9, frameon=True)

    text = (
        f"Matched points: {len(points)}\n"
        "pass-1: y=1, error=stat ⊕ 31% norm\n"
        r"pass-2: y=$\sigma_{p2}/\sigma_{p1}$, error=pass-2 stat $\oplus$ local $s_{\mathrm{comb}}$"
    )

    ax.text(
        0.01,
        0.02,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


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
    floor = log_y_min if logy and log_y_min > 0.0 else 1.0e-30

    for p in list(pass2) + list(pass1):
        extent = max(p.err_low, p.err_high, p.scale_syst_abs)
        vals.append(max(floor, p.xs - extent))
        vals.append(max(floor, p.xs + extent))
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
        if not math.isfinite(x) or not math.isfinite(y) or not math.isfinite(scale) or scale <= 0.0:
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
) -> None:
    if curves is not None:
        ax.plot(curves.phi, curves.bh, linewidth=2.0, linestyle="-", label="BH", zorder=2)
        ax.plot(curves.phi, curves.km15, linewidth=2.0, linestyle="--", label="KM15", zorder=2)
    #endif

    if panel.pass1:
        x = [p.phi + PASS1_PHI_OFFSET_DEG for p in panel.pass1]
        y = [p.xs for p in panel.pass1]
        yerr = [[p.err_low for p in panel.pass1], [p.err_high for p in panel.pass1]]

        container = ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="s",
            markersize=5,
            capsize=2,
            linewidth=1.2,
            linestyle="None",
            label=f"{pass1_label}: vertical errors",
            zorder=4,
        )

        color = container.lines[0].get_color()

        draw_scale_boxes(
            ax=ax,
            x_values=x,
            y_values=y,
            scale_values=[p.scale_syst_abs for p in panel.pass1],
            half_width=SCALE_BOX_HALF_WIDTH_DEG,
            label=f"{pass1_label}: 31% norm box",
            color=color,
            logy=logy,
            log_y_min=log_y_min,
        )
    #endif

    if panel.pass2:
        x = [p.phi + PASS2_PHI_OFFSET_DEG for p in panel.pass2]
        y = [p.xs for p in panel.pass2]
        yerr = [[p.err_low for p in panel.pass2], [p.err_high for p in panel.pass2]]

        container = ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            markersize=5,
            capsize=2,
            linewidth=1.2,
            linestyle="None",
            label=f"{pass2_label}: vertical errors",
            zorder=5,
        )

        color = container.lines[0].get_color()

        draw_scale_boxes(
            ax=ax,
            x_values=x,
            y_values=y,
            scale_values=[p.scale_syst_abs for p in panel.pass2],
            half_width=SCALE_BOX_HALF_WIDTH_DEG,
            label=f"{pass2_label}: local $s_{{comb}}$ box",
            color=color,
            logy=logy,
            log_y_min=log_y_min,
        )
    #endif

    key = panel.key

    if include_title:
        ax.set_title(
            rf"$x_B \in [{key.xb_min:.3g}, {key.xb_max:.3g}]$, "
            rf"$Q^2 \in [{key.q2_min:.3g}, {key.q2_max:.3g}]$ (GeV$^2$), "
            rf"$|t| \in [{key.t_min:.3g}, {key.t_max:.3g}]$ (GeV$^2$)"
        )
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
    ratio_points = pass2_over_pass1_ratio_points(panel)

    if ratio_points:
        container = ax.errorbar(
            [p.phi for p in ratio_points],
            [p.ratio for p in ratio_points],
            yerr=[[p.stat_err_low for p in ratio_points], [p.stat_err_high for p in ratio_points]],
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
            x_values=[p.phi for p in ratio_points],
            y_values=[p.ratio for p in ratio_points],
            scale_values=[p.scale_err for p in ratio_points],
            half_width=SCALE_BOX_HALF_WIDTH_DEG,
            label="scale box: local pass-2 ⊕ 31% pass-1",
            color=color,
            logy=False,
            log_y_min=1.0e-30,
        )
    #endif

    ax.axhline(1.0, linewidth=1.2, linestyle="--", zorder=2)

    ymin, ymax = ratio_axis_limits([ratio_points])

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
) -> None:
    if two_panel:
        fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.0))

        key = panel.key
        fig.suptitle(
            rf"$x_B \in [{key.xb_min:.3g}, {key.xb_max:.3g}]$, "
            rf"$Q^2 \in [{key.q2_min:.3g}, {key.q2_max:.3g}]$ (GeV$^2$), "
            rf"$|t| \in [{key.t_min:.3g}, {key.t_max:.3g}]$ (GeV$^2$)",
            fontsize=13,
        )

        draw_cross_section_axis(
            ax=axes[0],
            panel=panel,
            curves=curves,
            pass2_label=pass2_label,
            pass1_label=pass1_label,
            logy=logy,
            log_y_min=log_y_min,
            include_title=False,
        )

        axes[0].set_title("Cross sections")
        draw_ratio_axis(axes[1], panel)
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
        )

        fig.tight_layout()
    #endif

    fig.savefig(output_path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Worker.
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
                    if job.model_cfg.allow_missing_models:
                        curves = None
                        model_status = "model failed; plot written without models"
                    else:
                        raise RuntimeError(
                            f"Model calculation failed for "
                            f"xB=[{job.panel.key.xb_min}, {job.panel.key.xb_max}], "
                            f"Q2=[{job.panel.key.q2_min}, {job.panel.key.q2_max}], "
                            f"|t|=[{job.panel.key.t_min}, {job.panel.key.t_max}]: {exc}"
                        ) from exc
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
            "n_ratio_pass2_over_pass1": len(ratio_points),
            "has_models": curves is not None,
            "two_panel": job.two_panel,
            "pass2_local_s_comb_mean_in_panel": mean(pass2_local_scales),
            "pass2_local_s_comb_rms_in_panel": rms(pass2_local_scales),
            "pass2_local_s_comb_max_in_panel": max(pass2_local_scales) if pass2_local_scales else 0.0,
            "model_status": model_status,
        }

        elapsed = time.time() - t0

        return WorkerResult(
            ok=True,
            index=job.index,
            filename=job.filename,
            output_path=str(output_path),
            manifest_entry=manifest_entry,
            cache_key=cache_key,
            cache_entry=cache_entry,
            message=(
                f"[{job.index}/{job.total}] wrote {output_path} "
                f"({model_status}, pass2 points={len(job.panel.pass2)}, "
                f"pass1 points={len(job.panel.pass1)}, ratios={len(ratio_points)}, "
                f"elapsed={format_seconds(elapsed)})"
            ),
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


# ---------------------------------------------------------------------------
# Setup helpers.
# ---------------------------------------------------------------------------
def resolve_km15_cli(user_value: str) -> str:
    candidates: List[Path] = []

    if user_value:
        candidates.append(Path(user_value).expanduser())
    #endif

    candidates.append(ANALYSIS_DIR / "km15_cli.py")
    candidates.append(Path.cwd() / "km15_cli.py")

    for cand in candidates:
        if cand.exists() and cand.is_file():
            resolved = str(cand.resolve())
            log(f"Path setup: using KM15 CLI: {resolved}")
            return resolved
        #endif
    #endfor

    tried = "\n".join(f"  - {c}" for c in candidates)
    die(f"Could not find km15_cli.py. Tried:\n{tried}\nUse --km15-cli /full/path/to/km15_cli.py.")


def resolve_python_executable(user_value: str) -> str:
    text = str(user_value or "").strip()

    if text:
        expanded = Path(text).expanduser()

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
        return sys.executable
    #endif

    found = shutil.which("python3")

    if found:
        return found
    #endif

    return "python3"


def validate_dvcsgen_dir(path_text: str) -> str:
    path = Path(path_text).expanduser()
    exe = path / "dvcsgen"

    if not exe.exists():
        warn(
            f"dvcsgen executable was not found at {exe}. "
            "Model calculation will fail unless --dvcsgen-dir is corrected "
            "or --allow-missing-models / --skip-models is used."
        )
    else:
        log(f"Path setup: found dvcsgen executable: {exe}")
    #endif

    return str(path)


def build_panels(pass2: Dict[BinKey, List[DataPoint]], pass1: Dict[BinKey, List[DataPoint]]) -> List[PanelData]:
    keys = sorted(set(pass2.keys()) | set(pass1.keys()), key=key_sort_tuple)
    panels = [PanelData(key=key, pass2=pass2.get(key, []), pass1=pass1.get(key, [])) for key in keys]

    n_both = sum(1 for p in panels if p.pass2 and p.pass1)
    n_pass2_only = sum(1 for p in panels if p.pass2 and not p.pass1)
    n_pass1_only = sum(1 for p in panels if p.pass1 and not p.pass2)

    log(
        f"Panel setup: total={len(panels)}, both={n_both}, "
        f"pass-2 only={n_pass2_only}, pass-1 only={n_pass1_only}."
    )

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
    jobs: List[WorkerJob] = []
    n_cache_hit = 0
    n_cache_miss = 0

    for i, panel in enumerate(panels, start=1):
        filename = panel_filename(panel.key, i)
        output_path = output_dir / filename

        cached_model_entry: Optional[Dict[str, List[float]]] = None

        if not skip_models and model_cfg.use_cache:
            ckey = model_cache_key(panel.key, model_cfg)
            cached_model_entry = cache.get(ckey)

            if cached_model_entry is None:
                n_cache_miss += 1
            else:
                n_cache_hit += 1
            #endif
        elif not skip_models:
            n_cache_miss += 1
        #endif

        jobs.append(
            WorkerJob(
                index=i,
                total=len(panels),
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

    log(f"Job setup: jobs={len(jobs)}, model cache hits={n_cache_hit}, misses={n_cache_miss}.")

    return jobs, n_cache_hit, n_cache_miss


# ---------------------------------------------------------------------------
# CLI and main.
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Focused pass-2 vs pass-1 DVCS comparison plotter."
    )

    parser.add_argument("pass2_csv", type=Path, help="Pass-2 CSV. Must be first positional input.")
    parser.add_argument("pass1_csv", type=Path, help="Pass-1 / Lee CSV. Must be second positional input.")

    parser.add_argument("--output-dir", type=Path, default=Path("output/pass2_vs_pass1_model_comparison"))
    parser.add_argument("--pass2-xs-column", default=None, help="Override pass-2 central cross-section column.")
    parser.add_argument("--pass2-label", default="CLAS12 pass-2")
    parser.add_argument("--pass1-label", default="CLAS12 pass-1")

    parser.add_argument("--e-beam", type=float, default=10.604, help="Beam energy for BH/KM15 curves (GeV).")
    parser.add_argument("--dvcsgen-dir", default=os.environ.get("DVCSGEN_PATH", DEFAULT_DVCSGEN_DIR))
    parser.add_argument("--km15-cli", default=os.environ.get("KM15_CLI", DEFAULT_KM15_CLI))
    parser.add_argument("--py-km15", default=os.environ.get("PY_KM15", ""))
    parser.add_argument("--phi-dense", type=int, default=73)

    parser.add_argument("--model-cache", type=Path, default=None)
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--skip-models", action="store_true")
    parser.add_argument("--allow-missing-models", action="store_true")

    parser.add_argument("--linear-y", action="store_true")
    parser.add_argument("--log-y-min", type=float, default=DEFAULT_LOG_Y_MIN)
    parser.add_argument("--two-panel", action="store_true")
    parser.add_argument(
        "--no-point-to-point-systematics",
        action="store_true",
        help="Use stat-only vertical error bars for pass-1 and pass-2. Scale boxes remain.",
    )

    parser.add_argument(
        "--no-global-bin-plot",
        action="store_true",
        help="Skip the global matched-point bin-number pass-1/pass-2 comparison plot.",
    )
    parser.add_argument(
        "--global-ratio-y-min",
        type=float,
        default=None,
        help="Optional fixed y-axis minimum for pass1_pass2_bin_number_ratio_comparison.png.",
    )
    parser.add_argument(
        "--global-ratio-y-max",
        type=float,
        default=None,
        help="Optional fixed y-axis maximum for pass1_pass2_bin_number_ratio_comparison.png.",
    )

    parser.add_argument("--print-columns", action="store_true")
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--progress-every", type=int, default=1)
    parser.add_argument("--quiet-workers", action="store_true")
    parser.add_argument("--verbose-worker-models", action="store_true")

    return parser.parse_args()


def main() -> int:
    t0 = time.time()
    args = parse_args()

    log("Starting focused pass-2 vs pass-1 comparison plotter.")
    log(f"Script path: {SCRIPT_PATH}")
    log(f"Current working directory: {Path.cwd()}")

    log("Configuration:")
    log(f"  pass2_csv                      = {args.pass2_csv}")
    log(f"  pass1_csv                      = {args.pass1_csv}")
    log(f"  output_dir                     = {args.output_dir}")
    log(f"  two_panel                      = {args.two_panel}")
    log(f"  no_point_to_point_systematics  = {args.no_point_to_point_systematics}")
    log(f"  no_global_bin_plot             = {args.no_global_bin_plot}")
    log(f"  global_ratio_y_min             = {args.global_ratio_y_min}")
    log(f"  global_ratio_y_max             = {args.global_ratio_y_max}")
    log(f"  skip_models                    = {args.skip_models}")
    log(f"  e_beam                         = {args.e_beam:g} GeV")
    log(f"  log_y_min                      = {args.log_y_min:g}")
    log(f"  workers requested              = {args.workers}")

    if not args.pass2_csv.exists():
        die(f"Pass-2 CSV does not exist: {args.pass2_csv}")
    #endif

    if not args.pass1_csv.exists():
        die(f"Pass-1 CSV does not exist: {args.pass1_csv}")
    #endif

    if args.log_y_min <= 0.0:
        warn(f"Requested --log-y-min {args.log_y_min:g} is invalid; using {DEFAULT_LOG_Y_MIN:g}.")
        args.log_y_min = DEFAULT_LOG_Y_MIN
    #endif

    args.output_dir.mkdir(parents=True, exist_ok=True)

    cache_path = args.model_cache or (args.output_dir / "model_curve_cache.json")

    args.dvcsgen_dir = validate_dvcsgen_dir(args.dvcsgen_dir)
    args.km15_cli = resolve_km15_cli(args.km15_cli)
    args.py_km15 = resolve_python_executable(args.py_km15)

    pass2 = load_pass2_csv(
        path=args.pass2_csv,
        output_dir=args.output_dir,
        xs_column=args.pass2_xs_column,
        no_point_to_point_systematics=args.no_point_to_point_systematics,
        print_columns=args.print_columns,
    )

    pass1 = load_pass1_csv(
        path=args.pass1_csv,
        no_point_to_point_systematics=args.no_point_to_point_systematics,
        print_columns=args.print_columns,
    )

    panels = build_panels(pass2, pass1)

    if not panels:
        die("No panels to write.")
    #endif

    if not args.no_global_bin_plot:
        draw_global_bin_ratio_plot(
            output_dir=args.output_dir,
            panels=panels,
            pass2_label=args.pass2_label,
            pass1_label=args.pass1_label,
            user_ymin=args.global_ratio_y_min,
            user_ymax=args.global_ratio_y_max,
        )
    #endif

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

    cache: Dict[str, Dict[str, List[float]]] = {}

    if not args.skip_models:
        cache = load_model_cache(cache_path) if model_cfg.use_cache else {}
    else:
        log("Model setup: --skip-models provided; no BH/KM15 curves will be computed.")
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

    n_workers = max(1, min(int(args.workers), 5, len(jobs)))

    if args.workers > 5:
        warn(f"Requested --workers {args.workers}; capped to 5.")
    #endif

    progress_every = max(1, int(args.progress_every))

    manifest_by_index: Dict[int, Dict[str, object]] = {}
    failures: List[WorkerResult] = []
    completed = 0
    n_new_cache_entries = 0
    n_model_computed = 0
    n_model_cache_hit_completed = 0
    n_model_failed_but_plotted = 0
    processing_t0 = time.time()

    log(f"Processing: starting {len(jobs)} panel jobs with {n_workers} worker(s).")

    if n_workers == 1:
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
                eta = (len(jobs) - completed) / rate if rate > 0.0 else 0.0
                log(f"Progress: {completed}/{len(jobs)} complete, failures={len(failures)}, ETA≈{format_seconds(eta)}.")
            #endif
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_index = {executor.submit(process_one_panel, job): job.index for job in jobs}

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
                    eta = (len(jobs) - completed) / rate if rate > 0.0 else 0.0
                    log(f"Progress: {completed}/{len(jobs)} complete, failures={len(failures)}, ETA≈{format_seconds(eta)}.")
                #endif
            #endfor
        #endwith
    #endif

    if failures:
        first = failures[0]
        die(f"{len(failures)} panel job(s) failed. First failed panel={first.filename}. See traceback above.")
    #endif

    if not args.skip_models:
        save_model_cache(cache_path, cache)
        log(f"Cache summary: initial hits={n_cache_hit}, initial misses={n_cache_miss}, new entries={n_new_cache_entries}.")
    #endif

    manifest = [manifest_by_index[i] for i in sorted(manifest_by_index)]
    manifest_path = args.output_dir / "manifest.json"

    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
    #endwith

    total_dt = time.time() - t0

    log("Final summary:")
    log(f"  output directory              = {args.output_dir.resolve()}")
    log(f"  panels written                = {len(manifest)}")
    log(f"  global bin-number plot        = {not args.no_global_bin_plot}")
    log(f"  no point-to-point systematics = {args.no_point_to_point_systematics}")
    log(f"  two-panel mode                = {args.two_panel}")
    log(f"  workers used                  = {n_workers}")
    log(f"  model cache hits completed    = {n_model_cache_hit_completed}")
    log(f"  model curves newly computed   = {n_model_computed}")
    log(f"  model failures allowed        = {n_model_failed_but_plotted}")
    log(f"  log-y visible floor           = {args.log_y_min:g}")
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