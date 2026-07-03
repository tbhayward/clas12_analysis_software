#!/usr/bin/env python3
"""
diagnose_pass2_local_scomb.py
--------------------------------------------------------------------------
Standalone pass-2 local period-spread diagnostic script.

Purpose:
  Diagnose the point-by-point local pass-2 combination scale systematic
  computed from the four 10.6 GeV unpolarized period columns:

      Fa18 Inb, Fa18 Out, Sp18 Inb, Sp18 Out

  For each pass-2 row:
    1. Read the central 10.6 GeV unpolarized cross section.
    2. Read each valid period tuple.
    3. Build a stat-weighted period reference.
    4. Form period/reference ratios.
    5. Compute:

           s_obs(row)  = RMS_i(r_i - 1)
           s_stat(row) = RMS_i(stat_i / ref)
           s_comb(row) = sqrt(max(0, s_obs(row)^2 - s_stat(row)^2))

  The diagnostic writes:
    - point-by-point CSV of local s_obs, s_stat, s_comb, and period ratios,
    - summary CSVs,
    - s_comb histograms for >=2, >=3, and exactly 4 good periods,
    - s_comb profile plots versus xB, Q2, |t|, phi, e_theta, p_theta, g_theta,
    - optional scatter plots of s_comb versus the same variables,
    - top-outlier CSVs.

Usage:
    python3 diagnose_pass2_local_scomb.py \
        pass2_dvcs.csv \
        --output-dir output/pass2_vs_pass1_diagnostics

Useful options:
    --pass2-xs-column "column name"
        Override the central pass-2 cross-section column.

    --min-good-periods-for-profiles 2
        Minimum number of valid 10.6 GeV periods required for profile plots.

    --nonzero-only-profiles
        Exclude s_comb = 0 points from profile plots and profile CSVs.

    --hist-percentile 99.5
        Histogram x-range is clipped at this percentile, with overflow counted.

    --profile-bins 16
        Number of bins for continuous variables such as theta.

    --make-scatter-plots
        Also write scatter plots of local s_comb versus each variable.

    --top-outliers 100
        Number of largest local s_comb points to write into the outlier CSV.

Notes:
  * This script intentionally does not compare to pass-1.
  * This script intentionally does not call BH/KM15/dvcsgen.
  * This script produces PNG files only.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Constants.
# ---------------------------------------------------------------------------
PASS2_PERIODS_10P6_UNPOL = [
    "Fa18 Inb",
    "Fa18 Out",
    "Sp18 Inb",
    "Sp18 Out",
]

PASS2_CENTRAL_XS_CANDIDATES = [
    "normed cross sections, ep->epg, exp, 10.6 GeV, unpol",
    "cross sections, ep->epg, exp, 10.6 GeV, unpol",
]

DEFAULT_OUTPUT_DIR = Path("output/pass2_vs_pass1_diagnostics")

DEFAULT_HIST_PERCENTILE = 99.5
DEFAULT_HIST_BINS = 70
DEFAULT_PROFILE_BINS = 16
DEFAULT_TOP_OUTLIERS = 100


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
    period_ratios: Dict[str, float] = field(default_factory=dict)
    period_ratio_stats: Dict[str, float] = field(default_factory=dict)


@dataclass
class LocalScalePoint:
    key: BinKey
    row_index: int
    bin_index: int
    phi: float
    pass2_xs: float
    pass2_stat: float
    variables: Dict[str, float]
    period_values: Dict[str, TupleValue]
    local: LocalScaleResult


@dataclass
class OverallSummary:
    n_rows: int = 0
    n_valid_rows: int = 0
    n_bad_central: int = 0
    n_invalid_bin: int = 0
    n_ge2: int = 0
    n_ge3: int = 0
    n_eq4: int = 0
    n_nonzero_ge2: int = 0
    n_zero_ge2: int = 0
    mean_s_comb_ge2: float = 0.0
    median_s_comb_ge2: float = 0.0
    rms_s_comb_ge2: float = 0.0
    p90_s_comb_ge2: float = 0.0
    p95_s_comb_ge2: float = 0.0
    p99_s_comb_ge2: float = 0.0
    max_s_comb_ge2: float = 0.0


@dataclass
class VariableSpec:
    key: str
    label: str
    unit: str
    column_candidates: List[str]
    fallback_pair: Optional[Tuple[str, str]] = None


@dataclass
class ProfileRow:
    variable: str
    bin_low: float
    bin_high: float
    x_center: float
    n_points: int
    n_nonzero: int
    nonzero_fraction: float
    mean_s_comb: float
    median_s_comb: float
    rms_s_comb: float
    p90_s_comb: float
    p95_s_comb: float
    max_s_comb: float


# ---------------------------------------------------------------------------
# Variable specifications.
# ---------------------------------------------------------------------------
VARIABLE_SPECS = [
    VariableSpec(
        key="xB",
        label=r"$x_B$",
        unit="",
        column_candidates=[
            "xBavg", "xbavg", "xB_avg", "xb_avg", "xB average", "xb average",
            "xB", "xb",
        ],
        fallback_pair=("xb_min", "xb_max"),
    ),
    VariableSpec(
        key="Q2",
        label=r"$Q^2$",
        unit=r"GeV$^2$",
        column_candidates=[
            "Q2avg", "q2avg", "Q2_avg", "q2_avg", "Q2 average", "q2 average",
            "Q2", "q2",
        ],
        fallback_pair=("q2_min", "q2_max"),
    ),
    VariableSpec(
        key="t",
        label=r"$|t|$",
        unit=r"GeV$^2$",
        column_candidates=[
            "t_abs_avg", "t_absavg", "tabs_avg", "tabsavg",
            "tavg", "t_avg", "minus_t_avg", "neg_t_avg",
            "t_abs", "tabs", "minus_t", "neg_t",
        ],
        fallback_pair=("t_min", "t_max"),
    ),
    VariableSpec(
        key="phi",
        label=r"$\phi$",
        unit="deg",
        column_candidates=[
            "phiavg", "phi_avg", "phi_average", "phi",
        ],
        fallback_pair=None,
    ),
    VariableSpec(
        key="e_theta",
        label=r"$\theta_e$",
        unit="deg",
        column_candidates=[
            "e_theta_avg", "etheta_avg", "theta_e_avg", "electron_theta_avg",
            "e_thetaavg", "ethetaavg", "theta_eavg",
            "e_theta", "etheta", "theta_e", "electron_theta",
            "e_theta_deg", "theta_e_deg", "electron theta",
        ],
        fallback_pair=None,
    ),
    VariableSpec(
        key="p_theta",
        label=r"$\theta_p$",
        unit="deg",
        column_candidates=[
            "p_theta_avg", "ptheta_avg", "theta_p_avg", "proton_theta_avg",
            "p_thetaavg", "pthetaavg", "theta_pavg",
            "p_theta", "ptheta", "theta_p", "proton_theta",
            "p_theta_deg", "theta_p_deg", "proton theta",
        ],
        fallback_pair=None,
    ),
    VariableSpec(
        key="g_theta",
        label=r"$\theta_\gamma$",
        unit="deg",
        column_candidates=[
            "g_theta_avg", "gtheta_avg", "gamma_theta_avg", "photon_theta_avg",
            "phot_theta_avg", "theta_g_avg", "theta_gamma_avg",
            "g_thetaavg", "gthetaavg", "gamma_thetaavg", "photon_thetaavg",
            "g_theta", "gtheta", "gamma_theta", "photon_theta",
            "phot_theta", "theta_g", "theta_gamma",
            "g_theta_deg", "gamma_theta_deg", "photon_theta_deg",
            "gamma theta", "photon theta",
        ],
        fallback_pair=None,
    ),
]


# ---------------------------------------------------------------------------
# Logging and formatting.
# ---------------------------------------------------------------------------
def timestamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[local-scomb][{timestamp()}] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[local-scomb][{timestamp()}][warn] {msg}", file=sys.stderr, flush=True)


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


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return ""
    #endif

    return f"{value:.12g}"


def safe_filename_piece(text: str) -> str:
    out = str(text)
    out = out.replace(" ", "_")
    out = out.replace("/", "_")
    out = out.replace("\\", "_")
    out = out.replace("|", "")
    out = out.replace("$", "")
    out = out.replace("{", "")
    out = out.replace("}", "")
    out = out.replace("^", "")
    out = out.replace("_", "_")
    out = out.replace(".", "p")
    out = out.replace("-", "m")
    return out


# ---------------------------------------------------------------------------
# Numeric helpers.
# ---------------------------------------------------------------------------
def to_float(value: object, default: float = float("nan")) -> float:
    if value is None:
        return default
    #endif

    text = str(value).strip()

    if text == "":
        return default
    #endif

    try:
        out = float(text)
        return out if math.isfinite(out) else default
    except ValueError:
        match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)

        if not match:
            return default
        #endif

        try:
            out = float(match.group(0))
            return out if math.isfinite(out) else default
        except ValueError:
            return default
        #endtry
    #endtry


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


def finite_nonnegative(value: float) -> bool:
    return math.isfinite(value) and value >= 0.0


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


def round_key_value(value: float, ndigits: int = 6) -> float:
    return round(float(value), ndigits)


# ---------------------------------------------------------------------------
# CSV column helpers.
# ---------------------------------------------------------------------------
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


def pass2_period_cross_section_column(period: str) -> str:
    return f"normed cross sections, ep->epg, exp, {period}, unpol"


def detect_common_columns(fieldnames: Sequence[str]) -> Dict[str, str]:
    cols: Dict[str, str] = {}

    cols["valid"] = find_column(fieldnames, ["valid bin", "valid", "is valid"], required=False) or ""
    cols["bin_index"] = find_column(fieldnames, ["bin index", "bin", "idx"], required=False) or ""

    cols["xb_min"] = find_column(fieldnames, ["xBmin", "xbmin", "xB_min", "xb_min"])
    cols["xb_max"] = find_column(fieldnames, ["xBmax", "xbmax", "xB_max", "xb_max"])
    cols["q2_min"] = find_column(fieldnames, ["Q2min", "q2min", "Q2_min", "q2_min"])
    cols["q2_max"] = find_column(fieldnames, ["Q2max", "q2max", "Q2_max", "q2_max"])
    cols["t_min"] = find_column(fieldnames, ["t_abs_min", "tmin", "t_min", "minus_t_min", "neg_t_min"])
    cols["t_max"] = find_column(fieldnames, ["t_abs_max", "tmax", "t_max", "minus_t_max", "neg_t_max"])

    cols["phi_avg"] = find_column(fieldnames, ["phiavg", "phi_avg", "phi_average", "phi"], required=False) or ""
    cols["phi_min"] = find_column(fieldnames, ["phimin", "phi_min"], required=False) or ""
    cols["phi_max"] = find_column(fieldnames, ["phimax", "phi_max"], required=False) or ""

    if cols["phi_avg"] == "" and (cols["phi_min"] == "" or cols["phi_max"] == ""):
        die("Could not find phiavg/phi_avg/phi, or phimin plus phimax columns.")
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


def make_bin_key(row: Dict[str, str], cols: Dict[str, str]) -> BinKey:
    return BinKey(
        xb_min=round_key_value(to_float(row.get(cols["xb_min"], ""))),
        xb_max=round_key_value(to_float(row.get(cols["xb_max"], ""))),
        q2_min=round_key_value(to_float(row.get(cols["q2_min"], ""))),
        q2_max=round_key_value(to_float(row.get(cols["q2_max"], ""))),
        t_min=round_key_value(to_float(row.get(cols["t_min"], ""))),
        t_max=round_key_value(to_float(row.get(cols["t_max"], ""))),
    )


def row_phi(row: Dict[str, str], cols: Dict[str, str]) -> float:
    if cols["phi_avg"] != "":
        return to_float(row.get(cols["phi_avg"], ""))
    #endif

    return 0.5 * (
        to_float(row.get(cols["phi_min"], "")) +
        to_float(row.get(cols["phi_max"], ""))
    )


def variable_value_from_row(
    row: Dict[str, str],
    fieldnames: Sequence[str],
    cols: Dict[str, str],
    spec: VariableSpec,
) -> float:
    if spec.key == "phi":
        return row_phi(row, cols)
    #endif

    direct_col = find_column(fieldnames, spec.column_candidates, required=False)

    if direct_col is not None:
        value = to_float(row.get(direct_col, ""))
        if math.isfinite(value):
            return value
        #endif
    #endif

    if spec.fallback_pair is not None:
        lo_key, hi_key = spec.fallback_pair
        lo_col = cols.get(lo_key, "")
        hi_col = cols.get(hi_key, "")

        if lo_col != "" and hi_col != "":
            lo = to_float(row.get(lo_col, ""))
            hi = to_float(row.get(hi_col, ""))

            if math.isfinite(lo) and math.isfinite(hi):
                return 0.5 * (lo + hi)
            #endif
        #endif
    #endif

    return float("nan")


def detect_variable_columns(fieldnames: Sequence[str], cols: Dict[str, str]) -> None:
    log("Variable-column detection:")

    for spec in VARIABLE_SPECS:
        if spec.key == "phi":
            if cols["phi_avg"] != "":
                log(f"  {spec.key:8s} -> {cols['phi_avg']}")
            else:
                log(f"  {spec.key:8s} -> midpoint({cols['phi_min']}, {cols['phi_max']})")
            #endif

            continue
        #endif

        direct_col = find_column(fieldnames, spec.column_candidates, required=False)

        if direct_col is not None:
            log(f"  {spec.key:8s} -> {direct_col}")
        elif spec.fallback_pair is not None:
            lo_key, hi_key = spec.fallback_pair
            log(f"  {spec.key:8s} -> midpoint({cols[lo_key]}, {cols[hi_key]})")
        else:
            warn(f"  {spec.key:8s} -> not found; plots for this variable will be skipped.")
        #endif
    #endfor


# ---------------------------------------------------------------------------
# Tuple parsing and local scale calculation.
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


def compute_local_pass2_scale(period_values: Dict[str, TupleValue]) -> LocalScaleResult:
    valid_items = [
        (period, value)
        for period, value in period_values.items()
        if value.ok and math.isfinite(value.value) and math.isfinite(value.stat) and value.stat > 0.0
    ]

    out = LocalScaleResult(n_valid_periods=len(valid_items))

    if len(valid_items) < 2:
        return out
    #endif

    ref = combine_stat_weighted([value for _, value in valid_items])

    if not ref.ok or not math.isfinite(ref.value) or abs(ref.value) <= 0.0:
        return out
    #endif

    ratios: List[float] = []
    ratio_stats: List[float] = []
    period_ratios: Dict[str, float] = {}
    period_ratio_stats: Dict[str, float] = {}

    for period, value in valid_items:
        ratio = value.value / ref.value
        ratio_stat = abs(value.stat / ref.value)

        if math.isfinite(ratio) and math.isfinite(ratio_stat) and ratio_stat > 0.0:
            ratios.append(ratio)
            ratio_stats.append(ratio_stat)
            period_ratios[period] = ratio
            period_ratio_stats[period] = ratio_stat
        #endif
    #endfor

    if len(ratios) < 2:
        return out
    #endif

    residuals = [ratio - 1.0 for ratio in ratios]
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
    out.period_ratios = period_ratios
    out.period_ratio_stats = period_ratio_stats

    return out


# ---------------------------------------------------------------------------
# Loading pass-2 CSV.
# ---------------------------------------------------------------------------
def load_pass2_points(
    path: Path,
    xs_column: Optional[str],
    print_columns: bool,
) -> Tuple[List[LocalScalePoint], OverallSummary]:
    t0 = time.time()

    log(f"Opening pass-2 CSV: {path}")

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            die(f"Pass-2 CSV appears empty: {path}")
        #endif

        fieldnames = reader.fieldnames

        log(f"Detected {len(fieldnames)} columns.")

        if print_columns:
            log("Pass-2 columns:")
            for i, name in enumerate(fieldnames):
                log(f"  [{i:03d}] {name}")
            #endfor
        #endif

        cols = detect_common_columns(fieldnames)
        detect_variable_columns(fieldnames, cols)

        central_col = xs_column or find_column(fieldnames, PASS2_CENTRAL_XS_CANDIDATES)

        if "combination sys" in central_col.lower():
            die(
                "The selected pass-2 central column contains 'combination sys'. "
                "Use the central cross-section column instead. "
                f"Selected column: {central_col}"
            )
        #endif

        period_cols = {
            period: find_column(fieldnames, [pass2_period_cross_section_column(period)])
            for period in PASS2_PERIODS_10P6_UNPOL
        }

        log(f"Central pass-2 cross-section column: {central_col}")

        for period in PASS2_PERIODS_10P6_UNPOL:
            log(f"Period input {period:10s}: {period_cols[period]}")
        #endfor

        points: List[LocalScalePoint] = []
        summary = OverallSummary()

        for i_row, row in enumerate(reader, start=1):
            summary.n_rows += 1

            if not row_is_valid(row, cols["valid"]):
                summary.n_invalid_bin += 1
                continue
            #endif

            central = parse_tuple_cell(row.get(central_col, ""))

            if not central.ok or not finite_positive(central.value):
                summary.n_bad_central += 1
                continue
            #endif

            key = make_bin_key(row, cols)
            phi = row_phi(row, cols)
            bin_index = to_int(row.get(cols["bin_index"], ""), default=-1) if cols["bin_index"] != "" else -1

            variables = {
                spec.key: variable_value_from_row(row, fieldnames, cols, spec)
                for spec in VARIABLE_SPECS
            }

            period_values = {
                period: parse_tuple_cell(row.get(period_cols[period], ""))
                for period in PASS2_PERIODS_10P6_UNPOL
            }

            local = compute_local_pass2_scale(period_values)

            point = LocalScalePoint(
                key=key,
                row_index=i_row,
                bin_index=bin_index,
                phi=phi,
                pass2_xs=central.value,
                pass2_stat=central.stat,
                variables=variables,
                period_values=period_values,
                local=local,
            )

            points.append(point)
            summary.n_valid_rows += 1

            if local.n_valid_periods >= 2:
                summary.n_ge2 += 1

                if local.s_comb > 0.0:
                    summary.n_nonzero_ge2 += 1
                else:
                    summary.n_zero_ge2 += 1
                #endif
            #endif

            if local.n_valid_periods >= 3:
                summary.n_ge3 += 1
            #endif

            if local.n_valid_periods == 4:
                summary.n_eq4 += 1
            #endif
        #endfor
    #endwith

    values_ge2 = [
        p.local.s_comb
        for p in points
        if p.local.n_valid_periods >= 2 and math.isfinite(p.local.s_comb)
    ]

    summary.mean_s_comb_ge2 = mean(values_ge2)
    summary.median_s_comb_ge2 = median(values_ge2)
    summary.rms_s_comb_ge2 = rms(values_ge2)
    summary.p90_s_comb_ge2 = percentile(values_ge2, 90.0)
    summary.p95_s_comb_ge2 = percentile(values_ge2, 95.0)
    summary.p99_s_comb_ge2 = percentile(values_ge2, 99.0)
    summary.max_s_comb_ge2 = max(values_ge2) if values_ge2 else 0.0

    dt = time.time() - t0

    log(f"Finished reading pass-2 CSV in {format_seconds(dt)}.")
    log(f"Rows scanned={summary.n_rows}, valid central rows kept={summary.n_valid_rows}.")
    log(f"Invalid-bin rows skipped={summary.n_invalid_bin}, bad central rows skipped={summary.n_bad_central}.")
    log(f"Good-period counts: >=2={summary.n_ge2}, >=3={summary.n_ge3}, exactly4={summary.n_eq4}.")
    log(
        "s_comb summary for >=2 good periods: "
        f"mean={summary.mean_s_comb_ge2:.6g} ({100.0 * summary.mean_s_comb_ge2:.4g}%), "
        f"median={summary.median_s_comb_ge2:.6g} ({100.0 * summary.median_s_comb_ge2:.4g}%), "
        f"p90={summary.p90_s_comb_ge2:.6g} ({100.0 * summary.p90_s_comb_ge2:.4g}%), "
        f"p95={summary.p95_s_comb_ge2:.6g} ({100.0 * summary.p95_s_comb_ge2:.4g}%), "
        f"p99={summary.p99_s_comb_ge2:.6g} ({100.0 * summary.p99_s_comb_ge2:.4g}%), "
        f"max={summary.max_s_comb_ge2:.6g} ({100.0 * summary.max_s_comb_ge2:.4g}%)."
    )

    return points, summary


# ---------------------------------------------------------------------------
# CSV outputs.
# ---------------------------------------------------------------------------
def write_points_csv(output_dir: Path, points: Sequence[LocalScalePoint]) -> None:
    path = output_dir / "pass2_point_by_point_local_scomb_points.csv"

    log(f"Writing point-by-point CSV: {path}")

    headers = [
        "row_index", "bin_index",
        "xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max",
        "pass2_xs", "pass2_stat",
    ]

    for spec in VARIABLE_SPECS:
        headers.append(spec.key)
    #endfor

    headers.extend([
        "n_valid_periods",
        "ref_value", "ref_stat",
        "s_obs", "s_obs_percent",
        "s_stat", "s_stat_percent",
        "s_comb", "s_comb_percent",
        "min_ratio", "max_ratio", "half_width",
    ])

    for period in PASS2_PERIODS_10P6_UNPOL:
        safe_period = period.replace(" ", "_")
        headers.extend([
            f"{safe_period}_value",
            f"{safe_period}_stat",
            f"{safe_period}_ratio",
            f"{safe_period}_ratio_stat",
        ])
    #endfor

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)

        for point in points:
            local = point.local
            key = point.key

            row: List[object] = [
                point.row_index,
                point.bin_index,
                format_float(key.xb_min),
                format_float(key.xb_max),
                format_float(key.q2_min),
                format_float(key.q2_max),
                format_float(key.t_min),
                format_float(key.t_max),
                format_float(point.pass2_xs),
                format_float(point.pass2_stat),
            ]

            for spec in VARIABLE_SPECS:
                row.append(format_float(point.variables.get(spec.key, float("nan"))))
            #endfor

            row.extend([
                local.n_valid_periods,
                format_float(local.ref_value),
                format_float(local.ref_stat),
                format_float(local.s_obs),
                format_float(100.0 * local.s_obs),
                format_float(local.s_stat_exp),
                format_float(100.0 * local.s_stat_exp),
                format_float(local.s_comb),
                format_float(100.0 * local.s_comb),
                format_float(local.min_ratio),
                format_float(local.max_ratio),
                format_float(local.half_width),
            ])

            for period in PASS2_PERIODS_10P6_UNPOL:
                value = point.period_values.get(period, TupleValue())
                row.extend([
                    format_float(value.value) if value.ok else "",
                    format_float(value.stat) if value.ok else "",
                    format_float(local.period_ratios.get(period, float("nan"))),
                    format_float(local.period_ratio_stats.get(period, float("nan"))),
                ])
            #endfor

            writer.writerow(row)
        #endfor
    #endwith


def write_overall_summary_csv(output_dir: Path, summary: OverallSummary) -> None:
    path = output_dir / "pass2_point_by_point_local_scomb_summary.csv"

    log(f"Writing summary CSV: {path}")

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)

        writer.writerow([
            "n_rows", "n_valid_rows", "n_bad_central", "n_invalid_bin",
            "n_ge2", "n_ge3", "n_eq4",
            "n_nonzero_ge2", "n_zero_ge2",
            "mean_s_comb_ge2", "median_s_comb_ge2", "rms_s_comb_ge2",
            "p90_s_comb_ge2", "p95_s_comb_ge2", "p99_s_comb_ge2",
            "max_s_comb_ge2",
            "mean_s_comb_ge2_percent", "median_s_comb_ge2_percent", "rms_s_comb_ge2_percent",
            "p90_s_comb_ge2_percent", "p95_s_comb_ge2_percent", "p99_s_comb_ge2_percent",
            "max_s_comb_ge2_percent",
        ])

        writer.writerow([
            summary.n_rows,
            summary.n_valid_rows,
            summary.n_bad_central,
            summary.n_invalid_bin,
            summary.n_ge2,
            summary.n_ge3,
            summary.n_eq4,
            summary.n_nonzero_ge2,
            summary.n_zero_ge2,
            format_float(summary.mean_s_comb_ge2),
            format_float(summary.median_s_comb_ge2),
            format_float(summary.rms_s_comb_ge2),
            format_float(summary.p90_s_comb_ge2),
            format_float(summary.p95_s_comb_ge2),
            format_float(summary.p99_s_comb_ge2),
            format_float(summary.max_s_comb_ge2),
            format_float(100.0 * summary.mean_s_comb_ge2),
            format_float(100.0 * summary.median_s_comb_ge2),
            format_float(100.0 * summary.rms_s_comb_ge2),
            format_float(100.0 * summary.p90_s_comb_ge2),
            format_float(100.0 * summary.p95_s_comb_ge2),
            format_float(100.0 * summary.p99_s_comb_ge2),
            format_float(100.0 * summary.max_s_comb_ge2),
        ])
    #endwith


def write_outliers_csv(output_dir: Path, points: Sequence[LocalScalePoint], top_n: int) -> None:
    path = output_dir / "pass2_point_by_point_local_scomb_top_outliers.csv"

    selected = [
        p for p in points
        if p.local.n_valid_periods >= 2 and math.isfinite(p.local.s_comb)
    ]

    selected.sort(key=lambda p: p.local.s_comb, reverse=True)
    selected = selected[:max(0, top_n)]

    log(f"Writing top-outlier CSV: {path}")

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)

        headers = [
            "rank", "row_index", "bin_index",
            "xB", "Q2", "t", "phi",
            "e_theta", "p_theta", "g_theta",
            "pass2_xs", "pass2_stat",
            "n_valid_periods",
            "s_comb", "s_comb_percent",
            "s_obs", "s_obs_percent",
            "s_stat", "s_stat_percent",
            "min_ratio", "max_ratio", "half_width",
        ]

        for period in PASS2_PERIODS_10P6_UNPOL:
            safe_period = period.replace(" ", "_")
            headers.extend([
                f"{safe_period}_value",
                f"{safe_period}_stat",
                f"{safe_period}_ratio",
            ])
        #endfor

        writer.writerow(headers)

        for rank, point in enumerate(selected, start=1):
            local = point.local

            row: List[object] = [
                rank,
                point.row_index,
                point.bin_index,
                format_float(point.variables.get("xB", float("nan"))),
                format_float(point.variables.get("Q2", float("nan"))),
                format_float(point.variables.get("t", float("nan"))),
                format_float(point.variables.get("phi", float("nan"))),
                format_float(point.variables.get("e_theta", float("nan"))),
                format_float(point.variables.get("p_theta", float("nan"))),
                format_float(point.variables.get("g_theta", float("nan"))),
                format_float(point.pass2_xs),
                format_float(point.pass2_stat),
                local.n_valid_periods,
                format_float(local.s_comb),
                format_float(100.0 * local.s_comb),
                format_float(local.s_obs),
                format_float(100.0 * local.s_obs),
                format_float(local.s_stat_exp),
                format_float(100.0 * local.s_stat_exp),
                format_float(local.min_ratio),
                format_float(local.max_ratio),
                format_float(local.half_width),
            ]

            for period in PASS2_PERIODS_10P6_UNPOL:
                value = point.period_values.get(period, TupleValue())
                row.extend([
                    format_float(value.value) if value.ok else "",
                    format_float(value.stat) if value.ok else "",
                    format_float(local.period_ratios.get(period, float("nan"))),
                ])
            #endfor

            writer.writerow(row)
        #endfor
    #endwith


# ---------------------------------------------------------------------------
# Point filtering.
# ---------------------------------------------------------------------------
def select_points(
    points: Sequence[LocalScalePoint],
    min_good_periods: int,
    exact_good_periods: Optional[int] = None,
    nonzero_only: bool = False,
) -> List[LocalScalePoint]:
    selected: List[LocalScalePoint] = []

    for point in points:
        n = point.local.n_valid_periods

        if exact_good_periods is not None:
            if n != exact_good_periods:
                continue
            #endif
        elif n < min_good_periods:
            continue
        #endif

        if not math.isfinite(point.local.s_comb):
            continue
        #endif

        if point.local.s_comb < 0.0:
            continue
        #endif

        if nonzero_only and point.local.s_comb <= 0.0:
            continue
        #endif

        selected.append(point)
    #endfor

    return selected


# ---------------------------------------------------------------------------
# Histogram plotting.
# ---------------------------------------------------------------------------
def draw_s_comb_histogram(
    output_dir: Path,
    points: Sequence[LocalScalePoint],
    min_good_periods: int,
    exact_good_periods: Optional[int],
    nonzero_only: bool,
    hist_bins: int,
    hist_percentile: float,
    filename_suffix: str,
    title_suffix: str,
) -> None:
    selected = select_points(
        points=points,
        min_good_periods=min_good_periods,
        exact_good_periods=exact_good_periods,
        nonzero_only=nonzero_only,
    )

    values = [p.local.s_comb for p in selected if finite_nonnegative(p.local.s_comb)]

    if not values:
        warn(f"No s_comb values for histogram '{title_suffix}'; skipping.")
        return
    #endif

    path = output_dir / f"pass2_point_by_point_s_comb_distribution_{filename_suffix}.png"

    shown_max = percentile(values, hist_percentile)

    if shown_max <= 0.0:
        shown_max = max(values) if values else 1.0
    #endif

    if shown_max <= 0.0:
        shown_max = 1.0
    #endif

    shown_values = [v for v in values if v <= shown_max]
    overflow = len(values) - len(shown_values)

    if not shown_values:
        shown_values = values
        shown_max = max(values)
        overflow = 0
    #endif

    x_max = 1.10 * shown_max if shown_max > 0.0 else 1.0

    value_mean = mean(values)
    value_median = median(values)
    value_rms = rms(values)
    value_p90 = percentile(values, 90.0)
    value_p95 = percentile(values, 95.0)
    value_p99 = percentile(values, 99.0)
    value_max = max(values)

    fig, ax = plt.subplots(figsize=(10.0, 6.6))

    ax.hist(
        shown_values,
        bins=hist_bins,
        range=(0.0, x_max),
        alpha=0.35,
        label=r"$s_{\mathrm{comb}}$",
    )

    ax.axvline(value_mean, linewidth=1.4, linestyle="--", label="Mean")
    ax.axvline(value_median, linewidth=1.4, linestyle=":", label="Median")

    text = (
        rf"Point-by-point $s_{{\mathrm{{comb}}}}$" "\n"
        f"{title_suffix}\n"
        f"Entries: {len(values)}\n"
        f"Mean: {value_mean:.5g} ({100.0 * value_mean:.3g}%)\n"
        f"Median: {value_median:.5g} ({100.0 * value_median:.3g}%)\n"
        f"RMS: {value_rms:.5g} ({100.0 * value_rms:.3g}%)\n"
        f"p90: {value_p90:.5g} ({100.0 * value_p90:.3g}%)\n"
        f"p95: {value_p95:.5g} ({100.0 * value_p95:.3g}%)\n"
        f"p99: {value_p99:.5g} ({100.0 * value_p99:.3g}%)\n"
        f"Max: {value_max:.5g} ({100.0 * value_max:.3g}%)\n"
        f"Shown range ends at {hist_percentile:g}th percentile; overflow={overflow}"
    )

    ax.text(
        0.97,
        0.97,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.88},
    )

    ax.set_xlabel(r"$s_{\mathrm{comb}}=\sqrt{\max(0, s_{\mathrm{obs}}^2-s_{\mathrm{stat}}^2)}$")
    ax.set_ylabel("Point count")
    ax.set_title(f"Pass-2 point-by-point combination period-spread distribution ({title_suffix})")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="upper left", fontsize=10, frameon=True)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Wrote histogram: {path}")


def draw_all_histograms(
    output_dir: Path,
    points: Sequence[LocalScalePoint],
    hist_bins: int,
    hist_percentile: float,
) -> None:
    draw_s_comb_histogram(
        output_dir=output_dir,
        points=points,
        min_good_periods=2,
        exact_good_periods=None,
        nonzero_only=False,
        hist_bins=hist_bins,
        hist_percentile=hist_percentile,
        filename_suffix="ge2_good_periods",
        title_suffix="at least 2 good periods",
    )

    draw_s_comb_histogram(
        output_dir=output_dir,
        points=points,
        min_good_periods=3,
        exact_good_periods=None,
        nonzero_only=False,
        hist_bins=hist_bins,
        hist_percentile=hist_percentile,
        filename_suffix="ge3_good_periods",
        title_suffix="at least 3 good periods",
    )

    draw_s_comb_histogram(
        output_dir=output_dir,
        points=points,
        min_good_periods=4,
        exact_good_periods=4,
        nonzero_only=False,
        hist_bins=hist_bins,
        hist_percentile=hist_percentile,
        filename_suffix="4_good_periods",
        title_suffix="exactly 4 good periods",
    )

    draw_s_comb_histogram(
        output_dir=output_dir,
        points=points,
        min_good_periods=2,
        exact_good_periods=None,
        nonzero_only=True,
        hist_bins=hist_bins,
        hist_percentile=hist_percentile,
        filename_suffix="ge2_good_periods_nonzero_only",
        title_suffix="at least 2 good periods, nonzero only",
    )


# ---------------------------------------------------------------------------
# Profile binning and plotting.
# ---------------------------------------------------------------------------
def build_profile_rows_for_variable(
    points: Sequence[LocalScalePoint],
    spec: VariableSpec,
    min_good_periods: int,
    profile_bins: int,
    nonzero_only: bool,
) -> List[ProfileRow]:
    selected: List[Tuple[float, float]] = []

    for point in points:
        if point.local.n_valid_periods < min_good_periods:
            continue
        #endif

        x = point.variables.get(spec.key, float("nan"))
        y = point.local.s_comb

        if not math.isfinite(x) or not math.isfinite(y) or y < 0.0:
            continue
        #endif

        if nonzero_only and y <= 0.0:
            continue
        #endif

        selected.append((x, y))
    #endfor

    if not selected:
        return []
    #endif

    x_values = [x for x, _ in selected]
    unique_x = sorted(set(round(x, 8) for x in x_values))

    rows: List[ProfileRow] = []

    if len(unique_x) <= max(30, profile_bins):
        for ux in unique_x:
            vals = [y for x, y in selected if round(x, 8) == ux]
            if not vals:
                continue
            #endif

            n_nonzero = sum(1 for v in vals if v > 0.0)

            rows.append(
                ProfileRow(
                    variable=spec.key,
                    bin_low=ux,
                    bin_high=ux,
                    x_center=ux,
                    n_points=len(vals),
                    n_nonzero=n_nonzero,
                    nonzero_fraction=float(n_nonzero) / float(len(vals)) if vals else 0.0,
                    mean_s_comb=mean(vals),
                    median_s_comb=median(vals),
                    rms_s_comb=rms(vals),
                    p90_s_comb=percentile(vals, 90.0),
                    p95_s_comb=percentile(vals, 95.0),
                    max_s_comb=max(vals),
                )
            )
        #endfor

        return rows
    #endif

    xmin = min(x_values)
    xmax = max(x_values)

    if xmax <= xmin:
        vals = [y for _, y in selected]
        n_nonzero = sum(1 for v in vals if v > 0.0)

        return [
            ProfileRow(
                variable=spec.key,
                bin_low=xmin,
                bin_high=xmax,
                x_center=xmin,
                n_points=len(vals),
                n_nonzero=n_nonzero,
                nonzero_fraction=float(n_nonzero) / float(len(vals)) if vals else 0.0,
                mean_s_comb=mean(vals),
                median_s_comb=median(vals),
                rms_s_comb=rms(vals),
                p90_s_comb=percentile(vals, 90.0),
                p95_s_comb=percentile(vals, 95.0),
                max_s_comb=max(vals) if vals else 0.0,
            )
        ]
    #endif

    nbins = max(1, profile_bins)
    width = (xmax - xmin) / float(nbins)

    for ibin in range(nbins):
        lo = xmin + width * ibin
        hi = xmin + width * (ibin + 1)

        if ibin == nbins - 1:
            vals = [y for x, y in selected if lo <= x <= hi]
        else:
            vals = [y for x, y in selected if lo <= x < hi]
        #endif

        if not vals:
            continue
        #endif

        n_nonzero = sum(1 for v in vals if v > 0.0)

        rows.append(
            ProfileRow(
                variable=spec.key,
                bin_low=lo,
                bin_high=hi,
                x_center=0.5 * (lo + hi),
                n_points=len(vals),
                n_nonzero=n_nonzero,
                nonzero_fraction=float(n_nonzero) / float(len(vals)) if vals else 0.0,
                mean_s_comb=mean(vals),
                median_s_comb=median(vals),
                rms_s_comb=rms(vals),
                p90_s_comb=percentile(vals, 90.0),
                p95_s_comb=percentile(vals, 95.0),
                max_s_comb=max(vals),
            )
        )
    #endfor

    return rows


def axis_label(spec: VariableSpec) -> str:
    if spec.unit:
        return f"{spec.label} ({spec.unit})"
    #endif

    return spec.label


def draw_profile_plot(
    output_dir: Path,
    spec: VariableSpec,
    rows: Sequence[ProfileRow],
    min_good_periods: int,
    nonzero_only: bool,
) -> None:
    if not rows:
        warn(f"No profile rows for {spec.key}; skipping profile plot.")
        return
    #endif

    path = output_dir / f"pass2_local_s_comb_mean_by_{safe_filename_piece(spec.key)}_bin.png"

    x = [row.x_center for row in rows]
    y_mean = [row.mean_s_comb for row in rows]
    y_median = [row.median_s_comb for row in rows]
    y_p90 = [row.p90_s_comb for row in rows]

    fig, ax = plt.subplots(figsize=(9.6, 6.2))

    ax.plot(x, y_mean, marker="o", linewidth=1.4, label="Mean")
    ax.plot(x, y_median, marker="s", linewidth=1.4, linestyle="--", label="Median")
    ax.plot(x, y_p90, marker="^", linewidth=1.4, linestyle=":", label="p90")

    title_extra = f"n good periods ≥ {min_good_periods}"
    if nonzero_only:
        title_extra += ", nonzero only"
    #endif

    ax.set_xlabel(axis_label(spec))
    ax.set_ylabel(r"$s_{\mathrm{comb}}$")
    ax.set_title(rf"Pass-2 local $s_{{\mathrm{{comb}}}}$ versus {spec.label} ({title_extra})")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best", fontsize=10, frameon=True)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Wrote profile plot: {path}")


def draw_nonzero_fraction_plot(
    output_dir: Path,
    spec: VariableSpec,
    rows: Sequence[ProfileRow],
    min_good_periods: int,
) -> None:
    if not rows:
        warn(f"No profile rows for {spec.key}; skipping nonzero-fraction plot.")
        return
    #endif

    path = output_dir / f"pass2_local_s_comb_nonzero_fraction_by_{safe_filename_piece(spec.key)}_bin.png"

    x = [row.x_center for row in rows]
    y = [100.0 * row.nonzero_fraction for row in rows]

    fig, ax = plt.subplots(figsize=(9.6, 6.2))

    ax.plot(x, y, marker="o", linewidth=1.4)

    ax.set_xlabel(axis_label(spec))
    ax.set_ylabel(r"Nonzero $s_{\mathrm{comb}}$ fraction (%)")
    ax.set_ylim(0.0, 105.0)
    ax.set_title(rf"Fraction of points with nonzero local $s_{{\mathrm{{comb}}}}$ versus {spec.label} (n good periods ≥ {min_good_periods})")
    ax.grid(True, which="major", alpha=0.25)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Wrote nonzero-fraction plot: {path}")


def draw_scatter_plot(
    output_dir: Path,
    points: Sequence[LocalScalePoint],
    spec: VariableSpec,
    min_good_periods: int,
    nonzero_only: bool,
) -> None:
    selected_x: List[float] = []
    selected_y: List[float] = []

    for point in points:
        if point.local.n_valid_periods < min_good_periods:
            continue
        #endif

        x = point.variables.get(spec.key, float("nan"))
        y = point.local.s_comb

        if not math.isfinite(x) or not math.isfinite(y) or y < 0.0:
            continue
        #endif

        if nonzero_only and y <= 0.0:
            continue
        #endif

        selected_x.append(x)
        selected_y.append(y)
    #endfor

    if not selected_x:
        warn(f"No scatter points for {spec.key}; skipping scatter plot.")
        return
    #endif

    suffix = "_nonzero_only" if nonzero_only else ""
    path = output_dir / f"pass2_local_s_comb_vs_{safe_filename_piece(spec.key)}{suffix}.png"

    fig, ax = plt.subplots(figsize=(9.6, 6.2))

    ax.scatter(selected_x, selected_y, s=14, alpha=0.45)

    ax.set_xlabel(axis_label(spec))
    ax.set_ylabel(r"$s_{\mathrm{comb}}$")
    ax.set_title(rf"Pass-2 local $s_{{\mathrm{{comb}}}}$ point cloud versus {spec.label}")
    ax.grid(True, which="major", alpha=0.25)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

    log(f"Wrote scatter plot: {path}")


def write_profile_csv(output_dir: Path, profile_rows: Sequence[ProfileRow]) -> None:
    path = output_dir / "pass2_local_s_comb_profiles.csv"

    log(f"Writing profile CSV: {path}")

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)

        writer.writerow([
            "variable",
            "bin_low", "bin_high", "x_center",
            "n_points", "n_nonzero", "nonzero_fraction",
            "mean_s_comb", "median_s_comb", "rms_s_comb",
            "p90_s_comb", "p95_s_comb", "max_s_comb",
            "mean_s_comb_percent", "median_s_comb_percent", "rms_s_comb_percent",
            "p90_s_comb_percent", "p95_s_comb_percent", "max_s_comb_percent",
        ])

        for row in profile_rows:
            writer.writerow([
                row.variable,
                format_float(row.bin_low),
                format_float(row.bin_high),
                format_float(row.x_center),
                row.n_points,
                row.n_nonzero,
                format_float(row.nonzero_fraction),
                format_float(row.mean_s_comb),
                format_float(row.median_s_comb),
                format_float(row.rms_s_comb),
                format_float(row.p90_s_comb),
                format_float(row.p95_s_comb),
                format_float(row.max_s_comb),
                format_float(100.0 * row.mean_s_comb),
                format_float(100.0 * row.median_s_comb),
                format_float(100.0 * row.rms_s_comb),
                format_float(100.0 * row.p90_s_comb),
                format_float(100.0 * row.p95_s_comb),
                format_float(100.0 * row.max_s_comb),
            ])
        #endfor
    #endwith


def draw_all_profile_plots(
    output_dir: Path,
    points: Sequence[LocalScalePoint],
    min_good_periods: int,
    profile_bins: int,
    nonzero_only: bool,
    make_scatter_plots: bool,
) -> None:
    all_profile_rows: List[ProfileRow] = []

    for spec in VARIABLE_SPECS:
        rows = build_profile_rows_for_variable(
            points=points,
            spec=spec,
            min_good_periods=min_good_periods,
            profile_bins=profile_bins,
            nonzero_only=nonzero_only,
        )

        if not rows:
            warn(f"No valid values for variable {spec.key}; skipping its profile outputs.")
            continue
        #endif

        all_profile_rows.extend(rows)

        draw_profile_plot(
            output_dir=output_dir,
            spec=spec,
            rows=rows,
            min_good_periods=min_good_periods,
            nonzero_only=nonzero_only,
        )

        draw_nonzero_fraction_plot(
            output_dir=output_dir,
            spec=spec,
            rows=rows,
            min_good_periods=min_good_periods,
        )

        if make_scatter_plots:
            draw_scatter_plot(
                output_dir=output_dir,
                points=points,
                spec=spec,
                min_good_periods=min_good_periods,
                nonzero_only=nonzero_only,
            )
        #endif
    #endfor

    write_profile_csv(output_dir=output_dir, profile_rows=all_profile_rows)


# ---------------------------------------------------------------------------
# Argument parsing.
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Diagnose point-by-point pass-2 local s_comb period-spread systematics."
    )

    parser.add_argument("pass2_csv", type=Path, help="Pass-2 CSV.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--pass2-xs-column", default=None, help="Override central pass-2 cross-section column.")
    parser.add_argument("--print-columns", action="store_true", help="Print all detected CSV columns.")
    parser.add_argument(
        "--hist-bins",
        type=int,
        default=DEFAULT_HIST_BINS,
        help=f"Number of bins in s_comb histograms. Default {DEFAULT_HIST_BINS}.",
    )
    parser.add_argument(
        "--hist-percentile",
        type=float,
        default=DEFAULT_HIST_PERCENTILE,
        help=f"Histogram x-range upper percentile. Default {DEFAULT_HIST_PERCENTILE}.",
    )
    parser.add_argument(
        "--profile-bins",
        type=int,
        default=DEFAULT_PROFILE_BINS,
        help=f"Number of bins for continuous profile variables. Default {DEFAULT_PROFILE_BINS}.",
    )
    parser.add_argument(
        "--min-good-periods-for-profiles",
        type=int,
        default=2,
        help="Minimum number of valid period entries for profile/scatter plots. Default 2.",
    )
    parser.add_argument(
        "--nonzero-only-profiles",
        action="store_true",
        help="Exclude s_comb = 0 points from profile and scatter plots.",
    )
    parser.add_argument(
        "--make-scatter-plots",
        action="store_true",
        help="Also write scatter plots of local s_comb versus each kinematic variable.",
    )
    parser.add_argument(
        "--top-outliers",
        type=int,
        default=DEFAULT_TOP_OUTLIERS,
        help=f"Number of largest local s_comb outliers to write. Default {DEFAULT_TOP_OUTLIERS}.",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------
def main() -> int:
    t0 = time.time()

    args = parse_args()

    log("Starting pass-2 local s_comb diagnostics.")
    log(f"pass2_csv                     = {args.pass2_csv}")
    log(f"output_dir                    = {args.output_dir}")
    log(f"pass2_xs_column               = {args.pass2_xs_column}")
    log(f"hist_bins                     = {args.hist_bins}")
    log(f"hist_percentile               = {args.hist_percentile}")
    log(f"profile_bins                  = {args.profile_bins}")
    log(f"min_good_periods_for_profiles = {args.min_good_periods_for_profiles}")
    log(f"nonzero_only_profiles         = {args.nonzero_only_profiles}")
    log(f"make_scatter_plots            = {args.make_scatter_plots}")
    log(f"top_outliers                  = {args.top_outliers}")

    if not args.pass2_csv.exists():
        die(f"Pass-2 CSV does not exist: {args.pass2_csv}")
    #endif

    if args.hist_bins <= 0:
        warn(f"Requested --hist-bins={args.hist_bins}; using {DEFAULT_HIST_BINS}.")
        args.hist_bins = DEFAULT_HIST_BINS
    #endif

    if args.profile_bins <= 0:
        warn(f"Requested --profile-bins={args.profile_bins}; using {DEFAULT_PROFILE_BINS}.")
        args.profile_bins = DEFAULT_PROFILE_BINS
    #endif

    if args.min_good_periods_for_profiles < 2:
        warn("Requested --min-good-periods-for-profiles < 2; using 2.")
        args.min_good_periods_for_profiles = 2
    #endif

    if args.min_good_periods_for_profiles > 4:
        warn("Requested --min-good-periods-for-profiles > 4; using 4.")
        args.min_good_periods_for_profiles = 4
    #endif

    args.output_dir.mkdir(parents=True, exist_ok=True)

    points, summary = load_pass2_points(
        path=args.pass2_csv,
        xs_column=args.pass2_xs_column,
        print_columns=args.print_columns,
    )

    write_points_csv(output_dir=args.output_dir, points=points)
    write_overall_summary_csv(output_dir=args.output_dir, summary=summary)
    write_outliers_csv(output_dir=args.output_dir, points=points, top_n=args.top_outliers)

    draw_all_histograms(
        output_dir=args.output_dir,
        points=points,
        hist_bins=args.hist_bins,
        hist_percentile=args.hist_percentile,
    )

    draw_all_profile_plots(
        output_dir=args.output_dir,
        points=points,
        min_good_periods=args.min_good_periods_for_profiles,
        profile_bins=args.profile_bins,
        nonzero_only=args.nonzero_only_profiles,
        make_scatter_plots=args.make_scatter_plots,
    )

    dt = time.time() - t0

    log("Final summary:")
    log(f"  output directory       = {args.output_dir.resolve()}")
    log(f"  valid points           = {summary.n_valid_rows}")
    log(f"  points with >=2 periods= {summary.n_ge2}")
    log(f"  points with >=3 periods= {summary.n_ge3}")
    log(f"  points with 4 periods  = {summary.n_eq4}")
    log(f"  nonzero s_comb >=2     = {summary.n_nonzero_ge2}")
    log(f"  zero s_comb >=2        = {summary.n_zero_ge2}")
    log(f"  median s_comb >=2      = {summary.median_s_comb_ge2:.10g} ({100.0 * summary.median_s_comb_ge2:.6g}%)")
    log(f"  mean s_comb >=2        = {summary.mean_s_comb_ge2:.10g} ({100.0 * summary.mean_s_comb_ge2:.6g}%)")
    log(f"  p90 s_comb >=2         = {summary.p90_s_comb_ge2:.10g} ({100.0 * summary.p90_s_comb_ge2:.6g}%)")
    log(f"  p95 s_comb >=2         = {summary.p95_s_comb_ge2:.10g} ({100.0 * summary.p95_s_comb_ge2:.6g}%)")
    log(f"  p99 s_comb >=2         = {summary.p99_s_comb_ge2:.10g} ({100.0 * summary.p99_s_comb_ge2:.6g}%)")
    log(f"  max s_comb >=2         = {summary.max_s_comb_ge2:.10g} ({100.0 * summary.max_s_comb_ge2:.6g}%)")
    log(f"  elapsed                = {format_seconds(dt)}")
    log("Done.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[local-scomb][{timestamp()}][FATAL] {exc}", file=sys.stderr)
        raise SystemExit(1)
    #endtry