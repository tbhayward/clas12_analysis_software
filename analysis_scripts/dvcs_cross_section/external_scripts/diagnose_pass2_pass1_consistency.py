#!/usr/bin/env python3
"""
diagnose_pass2_pass1_consistency.py
--------------------------------------------------------------------------
Diagnostics for CLAS12 pass-2 local period-spread systematics and pass-2/pass-1
stat-only consistency.

This script intentionally does NOT make the main pass-2/pass-1 model-comparison
panels. Those are handled by plot_pass2_vs_pass1_model_comparison.py.

Diagnostics included:
  1. Point-by-point pass-2 local period-spread calculation:
       ref(row)    = stat-weighted mean of valid 10.6 GeV period tuples
       r_i(row)    = sigma_i / ref
       s_obs(row)  = RMS_i(r_i - 1)
       s_stat(row) = RMS_i(delta_i / ref)
       s_comb(row) = sqrt(max(0, s_obs(row)^2 - s_stat(row)^2))

  2. CSV outputs:
       pass2_point_by_point_scale_points.csv
       pass2_point_by_point_scale_summary.csv
       pass2_local_s_comb_kinematic_summary_by_bin.csv
       pass2_local_s_comb_top_outliers.csv

  3. Histograms:
       pass2_point_by_point_s_comb_distribution_ge2_good_periods.png
       pass2_point_by_point_s_comb_distribution_ge3_good_periods.png
       pass2_point_by_point_s_comb_distribution_4_good_periods.png
       pass2_point_by_point_s_obs_distribution_ge2_good_periods.png
       pass2_point_by_point_s_stat_distribution_ge2_good_periods.png

  4. Kinematic plots for nonzero local s_comb:
       pass2_local_s_comb_vs_phi.png
       pass2_local_s_comb_vs_xB.png
       pass2_local_s_comb_vs_Q2.png
       pass2_local_s_comb_vs_t.png
       pass2_local_s_comb_mean_by_xB_bin.png
       pass2_local_s_comb_mean_by_Q2_bin.png
       pass2_local_s_comb_mean_by_t_bin.png

  5. Stat-only pass-2/pass-1 pull diagnostic:
       z(N) = (pass2 - N * pass1) /
              sqrt(pass2_stat^2 + (N * pass1_stat)^2)

     N is a global pass-1 normalization floated within [0.69, 1.31].
     This deliberately excludes:
       - pass-2 18% point-to-point systematic,
       - pass-1 provided systematic,
       - pass-2 local s_comb scale systematic,
       - pass-1 31% normalization uncertainty.

Usage:
    python3 diagnose_pass2_pass1_consistency.py pass2.csv pass1.csv \
        --output-dir output/pass2_vs_pass1_diagnostics
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Constants.
# ---------------------------------------------------------------------------
PASS1_NORMALIZATION_FRACTION = 0.31
PULL_FIT_NORM_MIN = 1.0 - PASS1_NORMALIZATION_FRACTION
PULL_FIT_NORM_MAX = 1.0 + PASS1_NORMALIZATION_FRACTION

DEFAULT_PULL_MATCH_TOLERANCE_DEG = 2.0
DEFAULT_HIST_PERCENTILE_MAX = 99.5

PASS1_XS_COL = "cross sections, ep->epg, exp"
PASS1_STAT_COL = "cross sections, ep->epg, exp, stat. unc."
PASS1_SYST_UP_COL = "cross sections, ep->epg, exp, syst. unc. (up)"
PASS1_SYST_DN_COL = "cross sections, ep->epg, exp, syst. unc. (down)"

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
    fa18_inb_ratio: float = float("nan")
    fa18_out_ratio: float = float("nan")
    sp18_inb_ratio: float = float("nan")
    sp18_out_ratio: float = float("nan")


@dataclass
class LocalScalePoint:
    key: BinKey
    phi: float
    pass2_xs: float
    pass2_stat: float
    local: LocalScaleResult


@dataclass
class Pass1Point:
    key: BinKey
    phi: float
    xs: float
    stat: float
    syst_up: float
    syst_dn: float


@dataclass
class Pass2Point:
    key: BinKey
    phi: float
    xs: float
    stat: float


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
    best_norm: float = 1.0
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


# ---------------------------------------------------------------------------
# Logging and helpers.
# ---------------------------------------------------------------------------
def timestamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[dvcs-diagnostics][{timestamp()}] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[dvcs-diagnostics][{timestamp()}][warn] {msg}", file=sys.stderr, flush=True)


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


def canonical_header_map(fieldnames: Sequence[str]) -> Dict[str, str]:
    return {name.strip().lower(): name for name in fieldnames}


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


def percentile(values: Sequence[float], pct: float) -> float:
    vals = sorted(v for v in values if math.isfinite(v))

    if not vals:
        return 0.0
    #endif

    if pct <= 0.0:
        return vals[0]
    #endif

    if pct >= 100.0:
        return vals[-1]
    #endif

    pos = (pct / 100.0) * (len(vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))

    if lo == hi:
        return vals[lo]
    #endif

    frac = pos - lo

    return vals[lo] * (1.0 - frac) + vals[hi] * frac


def csv_escape(s: object) -> str:
    text = str(s)

    if "," not in text and '"' not in text and "\n" not in text and "\r" not in text:
        return text
    #endif

    return '"' + text.replace('"', '""') + '"'


# ---------------------------------------------------------------------------
# Common CSV kinematics.
# ---------------------------------------------------------------------------
def detect_common_columns(fieldnames: Sequence[str], label: str) -> Dict[str, str]:
    log(f"{label}: detecting common kinematic/valid-bin columns.")

    cols: Dict[str, str] = {}
    cols["valid"] = find_column(fieldnames, ["valid bin", "valid", "is valid"], required=False) or ""

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


def make_bin_key(row: Dict[str, str], cols: Dict[str, str]) -> BinKey:
    return BinKey(
        xb_min=round_key_value(to_float(row[cols["xb_min"]])),
        xb_max=round_key_value(to_float(row[cols["xb_max"]])),
        q2_min=round_key_value(to_float(row[cols["q2_min"]])),
        q2_max=round_key_value(to_float(row[cols["q2_max"]])),
        t_min=round_key_value(to_float(row[cols["t_min"]])),
        t_max=round_key_value(to_float(row[cols["t_max"]])),
    )


def key_center_xb(key: BinKey) -> float:
    return 0.5 * (key.xb_min + key.xb_max)


def key_center_q2(key: BinKey) -> float:
    return 0.5 * (key.q2_min + key.q2_max)


def key_center_t(key: BinKey) -> float:
    return 0.5 * (key.t_min + key.t_max)


def key_label(key: BinKey) -> str:
    return (
        f"xB=[{key.xb_min:.4g},{key.xb_max:.4g}], "
        f"Q2=[{key.q2_min:.4g},{key.q2_max:.4g}], "
        f"t=[{key.t_min:.4g},{key.t_max:.4g}]"
    )


# ---------------------------------------------------------------------------
# Tuple parsing and local scale.
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


def compute_local_pass2_scale(period_values_by_name: Dict[str, TupleValue]) -> LocalScaleResult:
    valid_by_name = {
        name: value
        for name, value in period_values_by_name.items()
        if value.ok and math.isfinite(value.value) and math.isfinite(value.stat) and value.stat > 0.0
    }

    out = LocalScaleResult(n_valid_periods=len(valid_by_name))

    if len(valid_by_name) < 2:
        return out
    #endif

    ref = combine_stat_weighted(list(valid_by_name.values()))

    if not ref.ok or not math.isfinite(ref.value) or abs(ref.value) <= 0.0:
        return out
    #endif

    ratios: List[float] = []
    ratio_stats: List[float] = []
    ratio_by_period: Dict[str, float] = {}

    for period, v in valid_by_name.items():
        ratio = v.value / ref.value
        ratio_stat = abs(v.stat / ref.value)

        if math.isfinite(ratio) and math.isfinite(ratio_stat) and ratio_stat > 0.0:
            ratios.append(ratio)
            ratio_stats.append(ratio_stat)
            ratio_by_period[period] = ratio
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
    out.fa18_inb_ratio = ratio_by_period.get("Fa18 Inb", float("nan"))
    out.fa18_out_ratio = ratio_by_period.get("Fa18 Out", float("nan"))
    out.sp18_inb_ratio = ratio_by_period.get("Sp18 Inb", float("nan"))
    out.sp18_out_ratio = ratio_by_period.get("Sp18 Out", float("nan"))

    return out


# ---------------------------------------------------------------------------
# Input loading.
# ---------------------------------------------------------------------------
def load_pass2_points(
    path: Path,
    xs_column: Optional[str],
    print_columns: bool,
) -> Tuple[List[Pass2Point], List[LocalScalePoint]]:
    t0 = time.time()

    log(f"Pass-2: opening {path}")

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
            die(f"Selected pass-2 central-value column contains 'combination sys': {xs_col}")
        #endif

        required_period_cols = [pass2_period_cross_section_column(period) for period in PASS2_PERIODS_10P6_UNPOL]
        require_columns(fieldnames, required_period_cols, "pass-2 local s_comb diagnostics")

        period_cols = {
            period: find_column(fieldnames, [pass2_period_cross_section_column(period)])
            for period in PASS2_PERIODS_10P6_UNPOL
        }

        pass2_points: List[Pass2Point] = []
        local_points: List[LocalScalePoint] = []

        total_rows = 0
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
            central = parse_tuple_cell(row.get(xs_col, ""))

            if not central.ok or not finite_positive(central.value):
                skipped_bad_xs += 1
                continue
            #endif

            period_values_by_name = {
                period: parse_tuple_cell(row.get(period_cols[period], ""))
                for period in PASS2_PERIODS_10P6_UNPOL
            }

            local = compute_local_pass2_scale(period_values_by_name)

            pass2_points.append(
                Pass2Point(
                    key=key,
                    phi=phi,
                    xs=central.value,
                    stat=central.stat,
                )
            )

            local_points.append(
                LocalScalePoint(
                    key=key,
                    phi=phi,
                    pass2_xs=central.value,
                    pass2_stat=central.stat,
                    local=local,
                )
            )
        #endfor
    #endwith

    dt = time.time() - t0

    log(
        f"Pass-2: read {len(pass2_points)} usable central points in {format_seconds(dt)} "
        f"(rows={total_rows}, skipped invalid={skipped_invalid}, skipped bad xs={skipped_bad_xs})."
    )

    return pass2_points, local_points


def load_pass1_points(path: Path, print_columns: bool) -> List[Pass1Point]:
    t0 = time.time()

    log(f"Pass-1: opening {path}")

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

        points: List[Pass1Point] = []

        total_rows = 0
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

            points.append(
                Pass1Point(
                    key=key,
                    phi=phi,
                    xs=xs,
                    stat=stat,
                    syst_up=syst_up,
                    syst_dn=syst_dn,
                )
            )
        #endfor
    #endwith

    dt = time.time() - t0

    log(
        f"Pass-1: read {len(points)} usable points in {format_seconds(dt)} "
        f"(rows={total_rows}, skipped invalid={skipped_invalid}, skipped bad xs={skipped_bad_xs})."
    )

    return points


# ---------------------------------------------------------------------------
# CSV outputs for local scale.
# ---------------------------------------------------------------------------
def write_local_scale_points_csv(output_dir: Path, points: Sequence[LocalScalePoint]) -> None:
    path = output_dir / "pass2_point_by_point_scale_points.csv"

    log(f"Writing local scale point CSV: {path}")

    with path.open("w", newline="") as handle:
        handle.write(
            "xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,"
            "xB_center,Q2_center,t_abs_center,phi,pass2_xs,pass2_stat,"
            "n_valid_periods,ref_value,ref_stat,"
            "s_obs,s_obs_percent,s_stat,s_stat_percent,"
            "s_comb,s_comb_percent,min_ratio,max_ratio,half_width,"
            "Fa18_Inb_ratio,Fa18_Out_ratio,Sp18_Inb_ratio,Sp18_Out_ratio\n"
        )

        for p in points:
            key = p.key
            local = p.local

            handle.write(
                f"{format_float(key.xb_min)},"
                f"{format_float(key.xb_max)},"
                f"{format_float(key.q2_min)},"
                f"{format_float(key.q2_max)},"
                f"{format_float(key.t_min)},"
                f"{format_float(key.t_max)},"
                f"{format_float(key_center_xb(key))},"
                f"{format_float(key_center_q2(key))},"
                f"{format_float(key_center_t(key))},"
                f"{format_float(p.phi)},"
                f"{format_float(p.pass2_xs)},"
                f"{format_float(p.pass2_stat)},"
                f"{local.n_valid_periods},"
                f"{format_float(local.ref_value)},"
                f"{format_float(local.ref_stat)},"
                f"{format_float(local.s_obs)},"
                f"{format_float(100.0 * local.s_obs)},"
                f"{format_float(local.s_stat_exp)},"
                f"{format_float(100.0 * local.s_stat_exp)},"
                f"{format_float(local.s_comb)},"
                f"{format_float(100.0 * local.s_comb)},"
                f"{format_float(local.min_ratio)},"
                f"{format_float(local.max_ratio)},"
                f"{format_float(local.half_width)},"
                f"{format_float(local.fa18_inb_ratio)},"
                f"{format_float(local.fa18_out_ratio)},"
                f"{format_float(local.sp18_inb_ratio)},"
                f"{format_float(local.sp18_out_ratio)}\n"
            )
        #endfor
    #endwith


def values_for(points: Sequence[LocalScalePoint], quantity: str, min_good: int = 2, exact_good: Optional[int] = None) -> List[float]:
    out: List[float] = []

    for p in points:
        n = p.local.n_valid_periods

        if exact_good is not None:
            if n != exact_good:
                continue
            #endif
        elif n < min_good:
            continue
        #endif

        if quantity == "s_comb":
            value = p.local.s_comb
        elif quantity == "s_obs":
            value = p.local.s_obs
        elif quantity == "s_stat":
            value = p.local.s_stat_exp
        elif quantity == "half_width":
            value = p.local.half_width
        else:
            die(f"Unknown quantity requested: {quantity}")
        #endif

        if math.isfinite(value) and value >= 0.0:
            out.append(value)
        #endif
    #endfor

    return out


def write_local_scale_summary_csv(output_dir: Path, points: Sequence[LocalScalePoint]) -> None:
    path = output_dir / "pass2_point_by_point_scale_summary.csv"

    log(f"Writing local scale summary CSV: {path}")

    selections = [
        ("ge2_good_periods", 2, None),
        ("ge3_good_periods", 3, None),
        ("4_good_periods", 4, 4),
        ("nonzero_s_comb_ge2", 2, None),
        ("nonzero_s_comb_ge3", 3, None),
        ("nonzero_s_comb_4", 4, 4),
    ]

    with path.open("w", newline="") as handle:
        handle.write(
            "selection,n_points,n_nonzero_s_comb,"
            "mean_s_comb,median_s_comb,rms_s_comb,max_s_comb,p90_s_comb,p95_s_comb,p99_s_comb,"
            "mean_s_obs,median_s_obs,rms_s_obs,max_s_obs,"
            "mean_s_stat,median_s_stat,rms_s_stat,max_s_stat\n"
        )

        for label, min_good, exact_good in selections:
            selected = [
                p for p in points
                if (p.local.n_valid_periods == exact_good if exact_good is not None else p.local.n_valid_periods >= min_good)
            ]

            if label.startswith("nonzero"):
                selected = [p for p in selected if p.local.s_comb > 0.0]
            #endif

            s_comb = [p.local.s_comb for p in selected if math.isfinite(p.local.s_comb)]
            s_obs = [p.local.s_obs for p in selected if math.isfinite(p.local.s_obs)]
            s_stat = [p.local.s_stat_exp for p in selected if math.isfinite(p.local.s_stat_exp)]
            nonzero = [v for v in s_comb if v > 0.0]

            handle.write(
                f"{label},"
                f"{len(selected)},"
                f"{len(nonzero)},"
                f"{format_float(mean(s_comb))},"
                f"{format_float(median(s_comb))},"
                f"{format_float(rms(s_comb))},"
                f"{format_float(max(s_comb) if s_comb else 0.0)},"
                f"{format_float(percentile(s_comb, 90.0))},"
                f"{format_float(percentile(s_comb, 95.0))},"
                f"{format_float(percentile(s_comb, 99.0))},"
                f"{format_float(mean(s_obs))},"
                f"{format_float(median(s_obs))},"
                f"{format_float(rms(s_obs))},"
                f"{format_float(max(s_obs) if s_obs else 0.0)},"
                f"{format_float(mean(s_stat))},"
                f"{format_float(median(s_stat))},"
                f"{format_float(rms(s_stat))},"
                f"{format_float(max(s_stat) if s_stat else 0.0)}\n"
            )
        #endfor
    #endwith


def write_top_outliers_csv(output_dir: Path, points: Sequence[LocalScalePoint], top_n: int) -> None:
    path = output_dir / "pass2_local_s_comb_top_outliers.csv"

    log(f"Writing top local s_comb outliers CSV: {path}")

    ranked = sorted(
        [p for p in points if math.isfinite(p.local.s_comb)],
        key=lambda p: p.local.s_comb,
        reverse=True,
    )

    with path.open("w", newline="") as handle:
        handle.write(
            "rank,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,"
            "xB_center,Q2_center,t_abs_center,phi,pass2_xs,pass2_stat,"
            "n_valid_periods,s_comb,s_comb_percent,s_obs,s_stat,"
            "min_ratio,max_ratio,half_width,"
            "Fa18_Inb_ratio,Fa18_Out_ratio,Sp18_Inb_ratio,Sp18_Out_ratio\n"
        )

        for rank, p in enumerate(ranked[:top_n], start=1):
            key = p.key
            local = p.local

            handle.write(
                f"{rank},"
                f"{format_float(key.xb_min)},"
                f"{format_float(key.xb_max)},"
                f"{format_float(key.q2_min)},"
                f"{format_float(key.q2_max)},"
                f"{format_float(key.t_min)},"
                f"{format_float(key.t_max)},"
                f"{format_float(key_center_xb(key))},"
                f"{format_float(key_center_q2(key))},"
                f"{format_float(key_center_t(key))},"
                f"{format_float(p.phi)},"
                f"{format_float(p.pass2_xs)},"
                f"{format_float(p.pass2_stat)},"
                f"{local.n_valid_periods},"
                f"{format_float(local.s_comb)},"
                f"{format_float(100.0 * local.s_comb)},"
                f"{format_float(local.s_obs)},"
                f"{format_float(local.s_stat_exp)},"
                f"{format_float(local.min_ratio)},"
                f"{format_float(local.max_ratio)},"
                f"{format_float(local.half_width)},"
                f"{format_float(local.fa18_inb_ratio)},"
                f"{format_float(local.fa18_out_ratio)},"
                f"{format_float(local.sp18_inb_ratio)},"
                f"{format_float(local.sp18_out_ratio)}\n"
            )
        #endfor
    #endwith


# ---------------------------------------------------------------------------
# Histogram plotting.
# ---------------------------------------------------------------------------
def draw_quantity_histogram(
    output_dir: Path,
    points: Sequence[LocalScalePoint],
    quantity: str,
    min_good: int,
    exact_good: Optional[int],
    filename_suffix: str,
    title_suffix: str,
    percentile_max: float,
) -> None:
    vals = values_for(points, quantity=quantity, min_good=min_good, exact_good=exact_good)

    if not vals:
        warn(f"No values for {quantity} histogram {title_suffix}; skipping.")
        return
    #endif

    path = output_dir / f"pass2_point_by_point_{quantity}_distribution_{filename_suffix}.png"

    x_shown_max = percentile(vals, percentile_max)
    x_shown_max = max(x_shown_max, 1.0e-12)

    shown = [v for v in vals if v <= x_shown_max]
    overflow = len(vals) - len(shown)

    if not shown:
        shown = vals
        x_shown_max = max(vals)
        overflow = 0
    #endif

    if x_shown_max <= 0.0:
        x_shown_max = 1.0
    #endif

    x_max = 1.10 * x_shown_max
    n_bins = 70

    label_map = {
        "s_comb": r"$s_{\mathrm{comb}}$",
        "s_obs": r"$s_{\mathrm{obs}}$",
        "s_stat": r"$s_{\mathrm{stat}}$",
        "half_width": "half width",
    }

    equation_map = {
        "s_comb": r"$s_{\mathrm{comb}}=\sqrt{\max(0,s_{\mathrm{obs}}^2-s_{\mathrm{stat}}^2)}$",
        "s_obs": r"$s_{\mathrm{obs}}=\sqrt{\langle(r_i-1)^2\rangle}$",
        "s_stat": r"$s_{\mathrm{stat}}=\sqrt{\langle(\delta_i/ref)^2\rangle}$",
        "half_width": r"$\frac{1}{2}(\max r_i-\min r_i)$",
    }

    fig, ax = plt.subplots(figsize=(10.0, 6.6))

    ax.hist(
        shown,
        bins=n_bins,
        range=(0.0, x_max),
        alpha=0.35,
        label=label_map.get(quantity, quantity),
    )

    value_mean = mean(vals)
    value_median = median(vals)
    value_rms = rms(vals)
    value_max = max(vals)

    ax.axvline(value_mean, linewidth=1.4, linestyle="--", label="Mean")
    ax.axvline(value_median, linewidth=1.4, linestyle=":", label="Median")

    text = (
        f"Point-by-point {label_map.get(quantity, quantity)}\n"
        f"{title_suffix}\n"
        f"Entries: {len(vals)}\n"
        f"Mean: {value_mean:.5g} ({100.0 * value_mean:.3g}%)\n"
        f"Median: {value_median:.5g} ({100.0 * value_median:.3g}%)\n"
        f"RMS: {value_rms:.5g} ({100.0 * value_rms:.3g}%)\n"
        f"Max: {value_max:.5g} ({100.0 * value_max:.3g}%)\n"
        f"Shown range ends at {percentile_max:g}th percentile; overflow={overflow}"
    )

    ax.text(
        0.97,
        0.97,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.88},
    )

    ax.set_xlabel(equation_map.get(quantity, quantity))
    ax.set_ylabel("Point count")
    ax.set_title(f"Pass-2 point-by-point {quantity} distribution ({title_suffix})")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="upper left", fontsize=10, frameon=True)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Wrote histogram: {path}")


def draw_all_histograms(output_dir: Path, points: Sequence[LocalScalePoint], percentile_max: float) -> None:
    draw_quantity_histogram(output_dir, points, "s_comb", 2, None, "ge2_good_periods", "at least 2 good periods", percentile_max)
    draw_quantity_histogram(output_dir, points, "s_comb", 3, None, "ge3_good_periods", "at least 3 good periods", percentile_max)
    draw_quantity_histogram(output_dir, points, "s_comb", 4, 4, "4_good_periods", "exactly 4 good periods", percentile_max)

    draw_quantity_histogram(output_dir, points, "s_obs", 2, None, "ge2_good_periods", "at least 2 good periods", percentile_max)
    draw_quantity_histogram(output_dir, points, "s_stat", 2, None, "ge2_good_periods", "at least 2 good periods", percentile_max)
    draw_quantity_histogram(output_dir, points, "half_width", 2, None, "ge2_good_periods", "at least 2 good periods", percentile_max)


# ---------------------------------------------------------------------------
# Kinematic diagnostics.
# ---------------------------------------------------------------------------
def nonzero_local_points(points: Sequence[LocalScalePoint], min_good_periods: int = 2) -> List[LocalScalePoint]:
    return [
        p for p in points
        if p.local.n_valid_periods >= min_good_periods and math.isfinite(p.local.s_comb) and p.local.s_comb > 0.0
    ]


def draw_s_comb_vs_variable(
    output_dir: Path,
    points: Sequence[LocalScalePoint],
    variable: str,
    xlabel: str,
    filename: str,
    min_good_periods: int,
) -> None:
    selected = nonzero_local_points(points, min_good_periods=min_good_periods)

    if not selected:
        warn(f"No nonzero local s_comb points for {filename}; skipping.")
        return
    #endif

    if variable == "phi":
        x = [p.phi for p in selected]
    elif variable == "xB":
        x = [key_center_xb(p.key) for p in selected]
    elif variable == "Q2":
        x = [key_center_q2(p.key) for p in selected]
    elif variable == "t":
        x = [key_center_t(p.key) for p in selected]
    else:
        die(f"Unknown kinematic variable: {variable}")
    #endif

    y = [p.local.s_comb for p in selected]

    fig, ax = plt.subplots(figsize=(8.4, 6.0))

    ax.scatter(x, y, s=14, alpha=0.55)
    ax.axhline(median(y), linewidth=1.2, linestyle=":", label=f"Median = {median(y):.4g}")
    ax.axhline(mean(y), linewidth=1.2, linestyle="--", label=f"Mean = {mean(y):.4g}")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Local $s_{\mathrm{comb}}$")
    ax.set_title(f"Nonzero local pass-2 $s_{{comb}}$ vs {xlabel}")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best", frameon=True)

    fig.tight_layout()
    path = output_dir / filename
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Wrote kinematic scatter plot: {path}")


def grouped_by_bin(points: Sequence[LocalScalePoint], variable: str, min_good_periods: int) -> Dict[str, List[float]]:
    out: Dict[str, List[float]] = {}

    for p in nonzero_local_points(points, min_good_periods=min_good_periods):
        key = p.key

        if variable == "xB":
            label = f"[{key.xb_min:.3g},{key.xb_max:.3g}]"
        elif variable == "Q2":
            label = f"[{key.q2_min:.3g},{key.q2_max:.3g}]"
        elif variable == "t":
            label = f"[{key.t_min:.3g},{key.t_max:.3g}]"
        else:
            die(f"Unknown grouped variable: {variable}")
        #endif

        out.setdefault(label, []).append(p.local.s_comb)
    #endfor

    return out


def draw_mean_s_comb_by_bin(
    output_dir: Path,
    points: Sequence[LocalScalePoint],
    variable: str,
    xlabel: str,
    filename: str,
    min_good_periods: int,
) -> None:
    groups = grouped_by_bin(points, variable=variable, min_good_periods=min_good_periods)

    if not groups:
        warn(f"No grouped values for {filename}; skipping.")
        return
    #endif

    labels = list(groups.keys())
    means = [mean(groups[label]) for label in labels]
    medians = [median(groups[label]) for label in labels]

    indices = list(range(len(labels)))

    fig, ax = plt.subplots(figsize=(10.0, 6.0))

    ax.plot(indices, means, marker="o", linestyle="-", label="Mean")
    ax.plot(indices, medians, marker="s", linestyle="--", label="Median")

    ax.set_xticks(indices)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Local $s_{\mathrm{comb}}$")
    ax.set_title(f"Nonzero local pass-2 $s_{{comb}}$ by {xlabel} bin")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best", frameon=True)

    fig.tight_layout()
    path = output_dir / filename
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Wrote grouped kinematic plot: {path}")


def write_kinematic_summary_by_bin(output_dir: Path, points: Sequence[LocalScalePoint], min_good_periods: int) -> None:
    path = output_dir / "pass2_local_s_comb_kinematic_summary_by_bin.csv"

    log(f"Writing kinematic summary CSV: {path}")

    bins: Dict[BinKey, List[LocalScalePoint]] = {}

    for p in points:
        if p.local.n_valid_periods < min_good_periods:
            continue
        #endif

        bins.setdefault(p.key, []).append(p)
    #endfor

    with path.open("w", newline="") as handle:
        handle.write(
            "xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,"
            "n_points,n_nonzero_s_comb,mean_s_comb,median_s_comb,rms_s_comb,max_s_comb,"
            "mean_s_obs,mean_s_stat\n"
        )

        for key in sorted(bins.keys(), key=lambda k: (k.xb_min, k.q2_min, k.t_min)):
            selected = bins[key]
            s_comb = [p.local.s_comb for p in selected if math.isfinite(p.local.s_comb)]
            s_comb_nonzero = [v for v in s_comb if v > 0.0]
            s_obs = [p.local.s_obs for p in selected if math.isfinite(p.local.s_obs)]
            s_stat = [p.local.s_stat_exp for p in selected if math.isfinite(p.local.s_stat_exp)]

            handle.write(
                f"{format_float(key.xb_min)},"
                f"{format_float(key.xb_max)},"
                f"{format_float(key.q2_min)},"
                f"{format_float(key.q2_max)},"
                f"{format_float(key.t_min)},"
                f"{format_float(key.t_max)},"
                f"{len(selected)},"
                f"{len(s_comb_nonzero)},"
                f"{format_float(mean(s_comb))},"
                f"{format_float(median(s_comb))},"
                f"{format_float(rms(s_comb))},"
                f"{format_float(max(s_comb) if s_comb else 0.0)},"
                f"{format_float(mean(s_obs))},"
                f"{format_float(mean(s_stat))}\n"
            )
        #endfor
    #endwith


def draw_kinematic_diagnostics(output_dir: Path, points: Sequence[LocalScalePoint], min_good_periods: int) -> None:
    draw_s_comb_vs_variable(output_dir, points, "phi", r"$\phi$ (deg)", "pass2_local_s_comb_vs_phi.png", min_good_periods)
    draw_s_comb_vs_variable(output_dir, points, "xB", r"$x_B$", "pass2_local_s_comb_vs_xB.png", min_good_periods)
    draw_s_comb_vs_variable(output_dir, points, "Q2", r"$Q^2$ (GeV$^2$)", "pass2_local_s_comb_vs_Q2.png", min_good_periods)
    draw_s_comb_vs_variable(output_dir, points, "t", r"$|t|$ (GeV$^2$)", "pass2_local_s_comb_vs_t.png", min_good_periods)

    draw_mean_s_comb_by_bin(output_dir, points, "xB", r"$x_B$", "pass2_local_s_comb_mean_by_xB_bin.png", min_good_periods)
    draw_mean_s_comb_by_bin(output_dir, points, "Q2", r"$Q^2$ (GeV$^2$)", "pass2_local_s_comb_mean_by_Q2_bin.png", min_good_periods)
    draw_mean_s_comb_by_bin(output_dir, points, "t", r"$|t|$ (GeV$^2$)", "pass2_local_s_comb_mean_by_t_bin.png", min_good_periods)

    write_kinematic_summary_by_bin(output_dir, points, min_good_periods=min_good_periods)


# ---------------------------------------------------------------------------
# Stat-only pull diagnostic.
# ---------------------------------------------------------------------------
def build_points_by_key_pass1(pass1_points: Sequence[Pass1Point]) -> Dict[BinKey, List[Pass1Point]]:
    out: Dict[BinKey, List[Pass1Point]] = {}

    for p in pass1_points:
        out.setdefault(p.key, []).append(p)
    #endfor

    for values in out.values():
        values.sort(key=lambda p: p.phi)
    #endfor

    return out


def find_nearest_pass1_point(phi: float, pass1_points: Sequence[Pass1Point], tolerance_deg: float) -> Optional[Pass1Point]:
    best: Optional[Pass1Point] = None
    best_delta = float("inf")

    for p in pass1_points:
        delta = abs(p.phi - phi)

        if delta < best_delta:
            best_delta = delta
            best = p
        #endif
    #endfor

    if best is None or best_delta > tolerance_deg:
        return None
    #endif

    return best


def build_stat_only_comparison_points(
    pass2_points: Sequence[Pass2Point],
    pass1_points: Sequence[Pass1Point],
    tolerance_deg: float,
) -> List[StatOnlyComparisonPoint]:
    pass1_by_key = build_points_by_key_pass1(pass1_points)
    out: List[StatOnlyComparisonPoint] = []

    for p2 in pass2_points:
        candidates = pass1_by_key.get(p2.key, [])

        if not candidates:
            continue
        #endif

        p1 = find_nearest_pass1_point(p2.phi, candidates, tolerance_deg=tolerance_deg)

        if p1 is None:
            continue
        #endif

        if not finite_positive(p2.xs) or not finite_positive(p1.xs):
            continue
        #endif

        if not finite_positive(p2.stat) or not finite_positive(p1.stat):
            continue
        #endif

        out.append(
            StatOnlyComparisonPoint(
                key=p2.key,
                phi_pass2=p2.phi,
                phi_pass1=p1.phi,
                pass2_xs=p2.xs,
                pass2_stat=p2.stat,
                pass1_xs=p1.xs,
                pass1_stat=p1.stat,
            )
        )
    #endfor

    return out


def stat_only_pull(point: StatOnlyComparisonPoint, pass1_norm: float) -> float:
    denom2 = point.pass2_stat * point.pass2_stat
    denom2 += (pass1_norm * point.pass1_stat) ** 2

    if denom2 <= 0.0 or not math.isfinite(denom2):
        return float("nan")
    #endif

    return (point.pass2_xs - pass1_norm * point.pass1_xs) / math.sqrt(denom2)


def stat_only_chi2(points: Sequence[StatOnlyComparisonPoint], pass1_norm: float) -> float:
    total = 0.0

    for p in points:
        z = stat_only_pull(p, pass1_norm)

        if math.isfinite(z):
            total += z * z
        #endif
    #endfor

    return total


def bounded_minimize_golden_section(func, xmin: float, xmax: float, tolerance: float = 1.0e-10, max_iter: int = 300) -> Tuple[float, float]:
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

    fmin = func(xmin)
    fmax = func(xmax)

    if fmin < fbest and fmin <= fmax:
        return xmin, fmin
    #endif

    if fmax < fbest and fmax < fmin:
        return xmax, fmax
    #endif

    return xbest, fbest


def percent_within(values: Sequence[float], threshold: float) -> float:
    vals = [v for v in values if math.isfinite(v)]

    if not vals:
        return 0.0
    #endif

    return 100.0 * sum(1 for v in vals if abs(v) <= threshold) / float(len(vals))


def summarize_stat_only_pulls(points: Sequence[StatOnlyComparisonPoint], best_norm: float) -> StatOnlyPullSummary:
    nominal_pulls = [stat_only_pull(p, 1.0) for p in points]
    best_pulls = [stat_only_pull(p, best_norm) for p in points]

    nominal_pulls = [z for z in nominal_pulls if math.isfinite(z)]
    best_pulls = [z for z in best_pulls if math.isfinite(z)]

    n_points = min(len(nominal_pulls), len(best_pulls))
    chi2_nominal = sum(z * z for z in nominal_pulls)
    chi2_best = sum(z * z for z in best_pulls)

    ndf_nominal = n_points
    ndf_best = max(1, n_points - 1)

    return StatOnlyPullSummary(
        n_points=n_points,
        best_norm=best_norm,
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


def write_pull_points_csv(output_dir: Path, points: Sequence[StatOnlyComparisonPoint], best_norm: float) -> None:
    path = output_dir / "pass1_pass2_stat_only_pull_points.csv"

    log(f"Writing stat-only pull point CSV: {path}")

    with path.open("w", newline="") as handle:
        handle.write(
            "xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,"
            "phi_pass2,phi_pass1,pass2_xs,pass2_stat,pass1_xs,pass1_stat,"
            "pull_nominal,pull_best_pass1_norm,chi2_nominal,chi2_best_pass1_norm,best_pass1_norm\n"
        )

        for p in points:
            pull_nominal = stat_only_pull(p, 1.0)
            pull_best = stat_only_pull(p, best_norm)

            handle.write(
                f"{format_float(p.key.xb_min)},"
                f"{format_float(p.key.xb_max)},"
                f"{format_float(p.key.q2_min)},"
                f"{format_float(p.key.q2_max)},"
                f"{format_float(p.key.t_min)},"
                f"{format_float(p.key.t_max)},"
                f"{format_float(p.phi_pass2)},"
                f"{format_float(p.phi_pass1)},"
                f"{format_float(p.pass2_xs)},"
                f"{format_float(p.pass2_stat)},"
                f"{format_float(p.pass1_xs)},"
                f"{format_float(p.pass1_stat)},"
                f"{format_float(pull_nominal)},"
                f"{format_float(pull_best)},"
                f"{format_float(pull_nominal * pull_nominal)},"
                f"{format_float(pull_best * pull_best)},"
                f"{format_float(best_norm)}\n"
            )
        #endfor
    #endwith


def write_pull_summary_csv(output_dir: Path, summary: StatOnlyPullSummary) -> None:
    path = output_dir / "pass1_pass2_stat_only_pull_summary.csv"

    log(f"Writing stat-only pull summary CSV: {path}")

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
            f"{format_float(summary.best_norm)},"
            f"{format_float(100.0 * (summary.best_norm - 1.0))},"
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


def draw_pull_histograms(output_dir: Path, points: Sequence[StatOnlyComparisonPoint], summary: StatOnlyPullSummary) -> None:
    nominal = [stat_only_pull(p, 1.0) for p in points]
    best = [stat_only_pull(p, summary.best_norm) for p in points]

    nominal = [z for z in nominal if math.isfinite(z)]
    best = [z for z in best if math.isfinite(z)]

    if not nominal or not best:
        warn("No valid stat-only pulls; skipping pull histograms.")
        return
    #endif

    finite = nominal + best
    xmin = min(-5.0, min(finite))
    xmax = max(5.0, max(finite))
    span = xmax - xmin

    if span <= 0.0:
        xmin -= 1.0
        xmax += 1.0
    else:
        xmin -= 0.05 * span
        xmax += 0.05 * span
    #endif

    fig, ax = plt.subplots(figsize=(9.4, 6.6))

    ax.hist(nominal, bins=80, range=(xmin, xmax), histtype="step", linewidth=1.8, label="N = 1")
    ax.hist(best, bins=80, range=(xmin, xmax), histtype="step", linewidth=1.8, label=f"Best pass-1 N = {summary.best_norm:.5f}")

    ax.axvline(0.0, linewidth=1.0, linestyle="-")
    ax.axvline(-1.0, linewidth=1.0, linestyle="--")
    ax.axvline(+1.0, linewidth=1.0, linestyle="--")
    ax.axvline(-3.0, linewidth=1.0, linestyle=":")
    ax.axvline(+3.0, linewidth=1.0, linestyle=":")

    text = (
        "Stat-only pass-2/pass-1 pulls\n"
        rf"$z=(\sigma_{{p2}}-N\sigma_{{p1}})/"
        rf"\sqrt{{\delta_{{p2,stat}}^2+(N\delta_{{p1,stat}})^2}}$" "\n"
        f"Matched points: {summary.n_points}\n"
        f"Best pass-1 N: {summary.best_norm:.5f} "
        f"({100.0 * (summary.best_norm - 1.0):+.2f}%)\n"
        f"χ²/ndf N=1: {summary.chi2_ndf_nominal:.3g} "
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
    ax.set_title("Pass-2 vs pass-1 stat-only pull distribution")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="upper left", frameon=True)

    fig.tight_layout()
    path = output_dir / "pass1_pass2_stat_only_pull_distribution.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Wrote pull histogram: {path}")

    chi2_nominal = [z * z for z in nominal if math.isfinite(z)]
    chi2_best = [z * z for z in best if math.isfinite(z)]

    xmax_chi2 = percentile(chi2_nominal + chi2_best, 99.0)
    xmax_chi2 = max(xmax_chi2, 1.0)

    fig2, ax2 = plt.subplots(figsize=(9.4, 6.6))

    ax2.hist(chi2_nominal, bins=80, range=(0.0, 1.10 * xmax_chi2), histtype="step", linewidth=1.8, label="N = 1")
    ax2.hist(chi2_best, bins=80, range=(0.0, 1.10 * xmax_chi2), histtype="step", linewidth=1.8, label=f"Best pass-1 N = {summary.best_norm:.5f}")

    ax2.set_xlabel(r"Point contribution to $\chi^2$ ($z^2$)")
    ax2.set_ylabel("Point count")
    ax2.set_title("Pass-2 vs pass-1 stat-only point-by-point χ² distribution")
    ax2.grid(True, which="major", alpha=0.25)
    ax2.legend(loc="best", frameon=True)

    fig2.tight_layout()
    path2 = output_dir / "pass1_pass2_stat_only_chi2_distribution.png"
    fig2.savefig(path2, dpi=200)
    plt.close(fig2)

    log(f"Wrote chi2 histogram: {path2}")


def run_pull_diagnostic(
    output_dir: Path,
    pass2_points: Sequence[Pass2Point],
    pass1_points: Sequence[Pass1Point],
    tolerance_deg: float,
) -> StatOnlyPullSummary:
    t0 = time.time()

    log("Pull diagnostic: building matched pass-2/pass-1 point list.")

    comparison_points = build_stat_only_comparison_points(
        pass2_points=pass2_points,
        pass1_points=pass1_points,
        tolerance_deg=tolerance_deg,
    )

    if not comparison_points:
        warn("Pull diagnostic: no matched points found.")
        summary = StatOnlyPullSummary()
        write_pull_summary_csv(output_dir, summary)
        return summary
    #endif

    log(f"Pull diagnostic: matched points={len(comparison_points)}.")
    log(
        "Pull diagnostic: minimizing stat-only global chi2 with pass-1 normalization "
        f"N in [{PULL_FIT_NORM_MIN:.3f}, {PULL_FIT_NORM_MAX:.3f}]."
    )

    best_norm, best_chi2 = bounded_minimize_golden_section(
        func=lambda norm: stat_only_chi2(comparison_points, norm),
        xmin=PULL_FIT_NORM_MIN,
        xmax=PULL_FIT_NORM_MAX,
    )

    summary = summarize_stat_only_pulls(comparison_points, best_norm)

    write_pull_points_csv(output_dir, comparison_points, best_norm)
    write_pull_summary_csv(output_dir, summary)
    draw_pull_histograms(output_dir, comparison_points, summary)

    log(
        "Pull diagnostic summary: "
        f"best pass-1 N={summary.best_norm:.10g} "
        f"({100.0 * (summary.best_norm - 1.0):+.5g}%), "
        f"best chi2={best_chi2:.10g}, "
        f"best chi2/ndf={summary.chi2_ndf_best:.10g}."
    )

    log(f"Pull diagnostic: finished in {format_seconds(time.time() - t0)}.")

    return summary


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Diagnostics for pass-2 local s_comb and pass-2/pass-1 consistency."
    )

    parser.add_argument("pass2_csv", type=Path, help="Pass-2 CSV.")
    parser.add_argument("pass1_csv", type=Path, help="Pass-1 / Lee CSV.")
    parser.add_argument("--output-dir", type=Path, default=Path("output/pass2_vs_pass1_diagnostics"))
    parser.add_argument("--pass2-xs-column", default=None, help="Override pass-2 central cross-section column.")
    parser.add_argument("--print-columns", action="store_true")

    parser.add_argument(
        "--hist-percentile-max",
        type=float,
        default=DEFAULT_HIST_PERCENTILE_MAX,
        help="Histogram shown x-range upper percentile. Default 99.5.",
    )
    parser.add_argument(
        "--top-outliers",
        type=int,
        default=100,
        help="Number of largest local s_comb points to write to outlier CSV. Default 100.",
    )
    parser.add_argument(
        "--kinematic-min-good-periods",
        type=int,
        default=2,
        help="Minimum valid periods for kinematic local s_comb diagnostics. Default 2.",
    )
    parser.add_argument(
        "--pull-match-tolerance-deg",
        type=float,
        default=DEFAULT_PULL_MATCH_TOLERANCE_DEG,
        help="Maximum phi mismatch for pass-2/pass-1 pull matching. Default 2 deg.",
    )

    parser.add_argument("--no-pull-diagnostic", action="store_true", help="Skip stat-only pass-2/pass-1 pull diagnostic.")
    parser.add_argument("--no-kinematic-plots", action="store_true", help="Skip local s_comb kinematic scatter/summary plots.")

    return parser.parse_args()


def main() -> int:
    t0 = time.time()
    args = parse_args()

    log("Starting pass-2/pass-1 diagnostic script.")
    log(f"  pass2_csv              = {args.pass2_csv}")
    log(f"  pass1_csv              = {args.pass1_csv}")
    log(f"  output_dir             = {args.output_dir}")
    log(f"  hist_percentile_max    = {args.hist_percentile_max:g}")
    log(f"  top_outliers           = {args.top_outliers}")
    log(f"  pull_match_tolerance   = {args.pull_match_tolerance_deg:g} deg")
    log(f"  no_pull_diagnostic     = {args.no_pull_diagnostic}")
    log(f"  no_kinematic_plots     = {args.no_kinematic_plots}")

    if not args.pass2_csv.exists():
        die(f"Pass-2 CSV does not exist: {args.pass2_csv}")
    #endif

    if not args.pass1_csv.exists():
        die(f"Pass-1 CSV does not exist: {args.pass1_csv}")
    #endif

    if args.hist_percentile_max <= 0.0 or args.hist_percentile_max > 100.0:
        warn(f"Invalid --hist-percentile-max {args.hist_percentile_max:g}; using {DEFAULT_HIST_PERCENTILE_MAX:g}.")
        args.hist_percentile_max = DEFAULT_HIST_PERCENTILE_MAX
    #endif

    args.output_dir.mkdir(parents=True, exist_ok=True)

    pass2_points, local_points = load_pass2_points(
        path=args.pass2_csv,
        xs_column=args.pass2_xs_column,
        print_columns=args.print_columns,
    )

    pass1_points = load_pass1_points(
        path=args.pass1_csv,
        print_columns=args.print_columns,
    )

    write_local_scale_points_csv(args.output_dir, local_points)
    write_local_scale_summary_csv(args.output_dir, local_points)
    write_top_outliers_csv(args.output_dir, local_points, top_n=max(1, int(args.top_outliers)))
    draw_all_histograms(args.output_dir, local_points, percentile_max=args.hist_percentile_max)

    if not args.no_kinematic_plots:
        draw_kinematic_diagnostics(
            output_dir=args.output_dir,
            points=local_points,
            min_good_periods=max(2, int(args.kinematic_min_good_periods)),
        )
    #endif

    pull_summary = StatOnlyPullSummary()

    if not args.no_pull_diagnostic:
        pull_summary = run_pull_diagnostic(
            output_dir=args.output_dir,
            pass2_points=pass2_points,
            pass1_points=pass1_points,
            tolerance_deg=args.pull_match_tolerance_deg,
        )
    #endif

    run_summary = {
        "pass2_csv": str(args.pass2_csv),
        "pass1_csv": str(args.pass1_csv),
        "n_pass2_points": len(pass2_points),
        "n_pass1_points": len(pass1_points),
        "n_local_scale_points": len(local_points),
        "n_local_s_comb_nonzero_ge2": len([p for p in local_points if p.local.n_valid_periods >= 2 and p.local.s_comb > 0.0]),
        "mean_s_comb_ge2": mean(values_for(local_points, "s_comb", 2, None)),
        "median_s_comb_ge2": median(values_for(local_points, "s_comb", 2, None)),
        "rms_s_comb_ge2": rms(values_for(local_points, "s_comb", 2, None)),
        "max_s_comb_ge2": max(values_for(local_points, "s_comb", 2, None)) if values_for(local_points, "s_comb", 2, None) else 0.0,
        "pull_diagnostic_enabled": not args.no_pull_diagnostic,
        "pull_matched_points": pull_summary.n_points,
        "pull_best_pass1_norm": pull_summary.best_norm,
        "pull_best_pass1_norm_percent": 100.0 * (pull_summary.best_norm - 1.0),
        "pull_chi2_ndf_nominal": pull_summary.chi2_ndf_nominal,
        "pull_chi2_ndf_best": pull_summary.chi2_ndf_best,
    }

    summary_path = args.output_dir / "diagnostic_run_summary.json"

    with summary_path.open("w") as handle:
        json.dump(run_summary, handle, indent=2)
    #endwith

    log("Final summary:")
    log(f"  output directory              = {args.output_dir.resolve()}")
    log(f"  pass2 points                  = {len(pass2_points)}")
    log(f"  pass1 points                  = {len(pass1_points)}")
    log(f"  local scale points            = {len(local_points)}")
    log(f"  nonzero local s_comb ge2      = {run_summary['n_local_s_comb_nonzero_ge2']}")
    log(f"  mean local s_comb ge2         = {run_summary['mean_s_comb_ge2']:.10g} ({100.0 * run_summary['mean_s_comb_ge2']:.6g}%)")
    log(f"  median local s_comb ge2       = {run_summary['median_s_comb_ge2']:.10g} ({100.0 * run_summary['median_s_comb_ge2']:.6g}%)")
    log(f"  rms local s_comb ge2          = {run_summary['rms_s_comb_ge2']:.10g} ({100.0 * run_summary['rms_s_comb_ge2']:.6g}%)")
    log(f"  max local s_comb ge2          = {run_summary['max_s_comb_ge2']:.10g} ({100.0 * run_summary['max_s_comb_ge2']:.6g}%)")

    if not args.no_pull_diagnostic:
        log(f"  stat-only matched points      = {pull_summary.n_points}")
        log(
            "  best pass1 normalization      = "
            f"{pull_summary.best_norm:.10g} ({100.0 * (pull_summary.best_norm - 1.0):+.6g}%)"
        )
        log(
            "  stat-only chi2/ndf            = "
            f"N=1 {pull_summary.chi2_ndf_nominal:.10g}, best {pull_summary.chi2_ndf_best:.10g}"
        )
    #endif

    log(f"  run summary JSON              = {summary_path}")
    log(f"  elapsed                       = {format_seconds(time.time() - t0)}")
    log("Done.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[dvcs-diagnostics][{timestamp()}][FATAL] {exc}", file=sys.stderr)
        raise SystemExit(1)
    #endtry