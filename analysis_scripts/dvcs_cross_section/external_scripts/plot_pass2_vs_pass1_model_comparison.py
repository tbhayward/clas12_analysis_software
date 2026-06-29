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

Default assumptions implemented here:
  * Positional input order is pass-2 CSV first, pass-1 CSV second.
  * Pass-2 data use the 10.6 GeV unpolarized normalized cross-section column.
    The statistical uncertainty is read from the tuple cell, then a 18% scale
    systematic is added in quadrature.
  * Pass-1 data use the Lee-style central/stat/syst-up/syst-down columns.
    The pass-1 total uncertainty adds stat, the provided syst column, and an
    additional 31% normalization uncertainty in quadrature.
  * One PNG is saved for every unique (xB, Q2, |t|) bin found in either CSV.
  * The model curves are BH and KM15 only, evaluated at the bin midpoint and
    drawn as functions of phi.

Parallelization:
  * Each unique (xB, Q2, |t|) panel is processed independently.
  * By default, up to 5 worker processes are used.
  * Use --workers N to change this.
  * Use --workers 1 for fully serial debugging.

Status output:
  * The script prints status updates during input parsing, column detection,
    model-tool path resolution, cache loading, job construction, model-cache
    hit/miss accounting, panel processing, plot writing, cache writing, and
    manifest writing.
  * Worker-level messages are returned to the parent process and printed in
    completion order, so parallel output remains readable.
  * Use --progress-every N to print aggregate progress every N completed jobs.
    The default is 1, meaning every completed panel prints a status line.
  * Use --quiet-workers to suppress per-panel completion messages while still
    keeping setup and final summary messages.

Efficiency:
  * Model curves are cached in model_curve_cache.json.
  * Workers never write the cache file directly. They return newly computed
    cache entries to the parent process, and the parent writes one final cache.
  * The default theory curve grid is 73 phi points, i.e. 5-degree spacing.
    This is usually visually smooth while avoiding unnecessary dvcsgen/km15
    subprocess calls. Use --phi-dense 181 if you want 2-degree spacing.
  * Matplotlib uses the non-interactive Agg backend for safe ifarm/batch use.

External model tools:
  * BH is evaluated by running dvcsgen. The directory can be provided with
    --dvcsgen-dir or the DVCSGEN_PATH environment variable.
  * KM15 is evaluated by running km15_cli.py. The path can be provided with
    --km15-cli or the KM15_CLI environment variable.
  * KM15 is launched with --py-km15 / PY_KM15 if provided; otherwise the
    current Python executable is used.
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


# ---------------------------------------------------------------------------
# Default external paths.
#
# The script is normally run from:
#     dvcs_cross_section/external_scripts/
# with km15_cli.py one directory above it.
# ---------------------------------------------------------------------------
SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
ANALYSIS_DIR = SCRIPT_DIR.parent

DEFAULT_DVCSGEN_DIR = "/u/home/thayward/dvcsgens/dvcsgen_print"
DEFAULT_KM15_CLI = str(ANALYSIS_DIR / "km15_cli.py")


# Column names used by the Sangbaek Lee / pass-1 cross-check file.
PASS1_XS_COL = "cross sections, ep->epg, exp"
PASS1_STAT_COL = "cross sections, ep->epg, exp, stat. unc."
PASS1_SYST_UP_COL = "cross sections, ep->epg, exp, syst. unc. (up)"
PASS1_SYST_DN_COL = "cross sections, ep->epg, exp, syst. unc. (down)"


# Pass-2 10.6 GeV column candidates. The script picks the first existing one.
PASS2_XS_CANDIDATES = [
    "normed cross sections, ep->epg, exp, 10.6 GeV, unpol",
    "normed cross sections, ep->epg, exp, 10.6 GeV, unpol, combination sys",
    "normed cross sections, ep->epg, exp, Fa18, unpol",
    "cross sections, ep->epg, exp, 10.6 GeV, unpol",
    "cross sections, ep->epg, exp, 10.6 GeV, unpol, combination sys",
]


@dataclass(frozen=True)
class BinKey:
    xb_min: float
    xb_max: float
    q2_min: float
    q2_max: float
    t_min: float
    t_max: float


@dataclass
class DataPoint:
    phi: float
    xs: float
    stat: float
    err_low: float
    err_high: float
    source: str


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
# Status and generic helpers
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


def to_float(value: object, default: float = 0.0) -> float:
    if value is None:
        return default
    #endif

    text = str(value).strip()

    if text == "":
        return default
    #endif

    try:
        return float(text)
    except ValueError:
        match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)

        if match:
            try:
                return float(match.group(0))
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


# ---------------------------------------------------------------------------
# CSV column detection
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
# Cross-section value parsing
# ---------------------------------------------------------------------------
def parse_tuple_cell(cell: str) -> Tuple[float, float, float, bool]:
    """
    Parse tuple-like pass-2 cells such as:
        (value, stat_err, syst_err)
        value, stat_err, syst_err

    Returns:
        (value, stat, syst, ok)
    """
    text = str(cell).strip()

    if text == "":
        return 0.0, 0.0, 0.0, False
    #endif

    values = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)]

    if len(values) >= 2:
        return values[0], abs(values[1]), abs(values[2]) if len(values) >= 3 else 0.0, True
    #endif

    if len(values) == 1:
        return values[0], 0.0, 0.0, True
    #endif

    return 0.0, 0.0, 0.0, False


def load_pass2_csv(path: Path, xs_column: Optional[str], print_columns: bool = False) -> Dict[BinKey, List[DataPoint]]:
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
        xs_col = xs_column or find_column(fieldnames, PASS2_XS_CANDIDATES)

        log(f"Pass-2: cross-section column -> {xs_col}")
        log("Pass-2: uncertainty prescription -> total = sqrt(stat^2 + (0.18 * cross_section)^2).")
        log("Pass-2: reading rows.")

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
            xs, stat, _ignored_syst, ok = parse_tuple_cell(row.get(xs_col, ""))

            if not ok or not finite_positive(xs):
                skipped_bad_xs += 1
                continue
            #endif

            err = math.sqrt(stat * stat + (0.18 * xs) * (0.18 * xs))
            point = DataPoint(phi=phi, xs=xs, stat=stat, err_low=err, err_high=err, source="pass2")

            out.setdefault(key, []).append(point)
            kept += 1
        #endfor
    #endwith

    for points in out.values():
        points.sort(key=lambda p: p.phi)
    #endfor

    dt = time.time() - t0

    log(f"Pass-2: finished reading in {format_seconds(dt)}.")
    log(f"Pass-2: total rows={total_rows}, kept={kept}, skipped invalid-bin={skipped_invalid}, skipped bad/nonpositive xs={skipped_bad_xs}.")
    log(f"Pass-2: unique (xB,Q2,|t|) bins with data={len(out)}.")

    return out


def load_pass1_csv(path: Path, print_columns: bool = False) -> Dict[BinKey, List[DataPoint]]:
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
        log("Pass-1: uncertainty prescription -> total = sqrt(stat^2 + syst^2 + (0.31 * cross_section)^2).")
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

            norm = 0.31 * xs
            err_high = math.sqrt(stat * stat + syst_up * syst_up + norm * norm)
            err_low = math.sqrt(stat * stat + syst_dn * syst_dn + norm * norm)

            if err_low <= 0.0:
                err_low = err_high
            #endif

            if err_high <= 0.0:
                err_high = err_low
            #endif

            point = DataPoint(phi=phi, xs=xs, stat=stat, err_low=err_low, err_high=err_high, source="pass1")

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
# Model calculation helpers
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
                f"xB=[{key.xb_min},{key.xb_max}], Q2=[{key.q2_min},{key.q2_max}], |t|=[{key.t_min},{key.t_max}], phi={phi:.1f}",
                flush=True,
            )
        #endif

        bh_values.append(bh_xs(xb, q2, t_abs, phi, cfg))
        km15_values.append(km15_xs(xb, q2, t_abs, phi, cfg))
    #endfor

    return ModelCurves(phi=phi_grid, bh=bh_values, km15=km15_values)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def y_limits(pass2: Sequence[DataPoint], pass1: Sequence[DataPoint], curves: Optional[ModelCurves]) -> Tuple[float, float]:
    vals: List[float] = []

    for p in list(pass2) + list(pass1):
        vals.append(max(1e-30, p.xs - p.err_low))
        vals.append(max(1e-30, p.xs + p.err_high))
    #endfor

    if curves is not None:
        vals.extend(v for v in curves.bh if finite_positive(v))
        vals.extend(v for v in curves.km15 if finite_positive(v))
    #endif

    vals = [v for v in vals if finite_positive(v)]

    if not vals:
        return 1e-4, 1.0
    #endif

    ymin = min(vals)
    ymax = max(vals)

    if ymin <= 0.0 or not math.isfinite(ymin):
        positive_vals = [v for v in vals if v > 0.0]

        if positive_vals:
            ymin = min(positive_vals)
        else:
            ymin = 1e-4
        #endif
    #endif

    if ymax <= ymin:
        ymax = 10.0 * ymin
    #endif

    return 0.55 * ymin, 1.75 * ymax


def draw_panel(
    panel: PanelData,
    curves: Optional[ModelCurves],
    output_path: Path,
    pass2_label: str,
    pass1_label: str,
    logy: bool,
) -> None:
    phi_offset_deg = 1.25

    fig, ax = plt.subplots(figsize=(8.4, 6.0))

    if curves is not None:
        ax.plot(curves.phi, curves.bh, linewidth=2.0, linestyle="-", label="BH")
        ax.plot(curves.phi, curves.km15, linewidth=2.0, linestyle="--", label="KM15")
    #endif

    if panel.pass1:
        x = [p.phi - phi_offset_deg for p in panel.pass1]
        y = [p.xs for p in panel.pass1]
        yerr = [[p.err_low for p in panel.pass1], [p.err_high for p in panel.pass1]]

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="s",
            markersize=5,
            capsize=2,
            linewidth=1.2,
            linestyle="None",
            label=f"{pass1_label}: stat ⊕ syst ⊕ 31% norm",
        )
    #endif

    if panel.pass2:
        x = [p.phi + phi_offset_deg for p in panel.pass2]
        y = [p.xs for p in panel.pass2]
        yerr = [p.err_high for p in panel.pass2]

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            markersize=5,
            capsize=2,
            linewidth=1.2,
            linestyle="None",
            label=f"{pass2_label}: stat ⊕ estimated systematics",
        )
    #endif

    key = panel.key
    title = (
        rf"$x_B \in [{key.xb_min:.3g}, {key.xb_max:.3g}]$, "
        rf"$Q^2 \in [{key.q2_min:.3g}, {key.q2_max:.3g}]$ (GeV$^2$), "
        rf"$|t| \in [{key.t_min:.3g}, {key.t_max:.3g}]$ (GeV$^2$)"
    )

    ax.set_title(title)
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$d\sigma/(dx_B\,dQ^2\,d|t|\,d\phi)$ (pb/GeV$^4$/rad)")
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])

    ymin, ymax = y_limits(panel.pass2, panel.pass1, curves)
    ax.set_ylim(ymin, ymax)

    if logy:
        ax.set_yscale("log")
    #endif

    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best", fontsize=9, frameon=True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Parallel worker
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
        )

        manifest_entry: Dict[str, object] = {
            "file": job.filename,
            "xB": [job.panel.key.xb_min, job.panel.key.xb_max],
            "Q2": [job.panel.key.q2_min, job.panel.key.q2_max],
            "t_abs": [job.panel.key.t_min, job.panel.key.t_max],
            "n_pass2_points": len(job.panel.pass2),
            "n_pass1_points": len(job.panel.pass1),
            "has_models": curves is not None,
            "model_status": model_status,
        }

        elapsed = time.time() - t0

        message = (
            f"[{job.index}/{job.total}] wrote {output_path} "
            f"({model_status}, pass2 points={len(job.panel.pass2)}, pass1 points={len(job.panel.pass1)}, "
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


# ---------------------------------------------------------------------------
# External-tool path resolution
# ---------------------------------------------------------------------------
def resolve_km15_cli(user_value: str) -> str:
    """
    Resolve km15_cli.py robustly when the script is run from external_scripts/.

    Priority:
      1. explicit --km15-cli / KM15_CLI value,
      2. ../km15_cli.py relative to this script,
      3. km15_cli.py in the current working directory.
    """
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
    """
    Resolve the Python executable used for KM15.

    An explicit --py-km15/PY_KM15 is only used if it exists or can be found
    on PATH. Otherwise the current Python executable is used.
    """
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
# Main driver helpers
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
    parser.add_argument("--pass2-xs-column", default=None, help="Override the pass-2 cross-section column name.")
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

    log("Command-line configuration:")
    log(f"  pass2_csv              = {args.pass2_csv}")
    log(f"  pass1_csv              = {args.pass1_csv}")
    log(f"  output_dir             = {args.output_dir}")
    log(f"  e_beam                 = {args.e_beam:g} GeV")
    log(f"  requested workers      = {args.workers}")
    log(f"  phi_dense              = {args.phi_dense}")
    log(f"  skip_models            = {args.skip_models}")
    log(f"  no_cache               = {args.no_cache}")
    log(f"  allow_missing_models   = {args.allow_missing_models}")
    log(f"  linear_y               = {args.linear_y}")
    log(f"  quiet_workers          = {args.quiet_workers}")
    log(f"  progress_every         = {args.progress_every}")

    if not args.pass2_csv.exists():
        die(f"Pass-2 CSV does not exist: {args.pass2_csv}")
    #endif

    if not args.pass1_csv.exists():
        die(f"Pass-1 CSV does not exist: {args.pass1_csv}")
    #endif

    log("Output setup: creating output directory if needed.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log(f"Output setup: output directory ready: {args.output_dir.resolve()}")

    cache_path = args.model_cache or (args.output_dir / "model_curve_cache.json")
    log(f"Output setup: model cache path: {cache_path}")

    args.dvcsgen_dir = validate_dvcsgen_dir(args.dvcsgen_dir)
    args.km15_cli = resolve_km15_cli(args.km15_cli)
    args.py_km15 = resolve_python_executable(args.py_km15)

    log("Input setup: loading pass-2 data.")
    pass2 = load_pass2_csv(args.pass2_csv, args.pass2_xs_column, print_columns=args.print_columns)

    log("Input setup: loading pass-1 data.")
    pass1 = load_pass1_csv(args.pass1_csv, print_columns=args.print_columns)

    panels = build_panels(pass2, pass1)

    if not panels:
        die("No panels to write after loading pass-2 and pass-1 CSVs.")
    #endif

    log(f"Panel setup: unique (xB,Q2,|t|) panels to write: {len(panels)}")

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
        manifest.append(manifest_by_index[index])
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
    log(f"  output directory          = {args.output_dir.resolve()}")
    log(f"  total panels written      = {len(manifest)}")
    log(f"  workers used              = {n_workers}")
    log(f"  model phi points          = {model_cfg.phi_dense}")
    log(f"  cache entries after run   = {len(cache) if not args.skip_models else 0}")
    log(f"  total elapsed             = {format_seconds(total_dt)}")
    log("Done.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[pass2-vs-pass1][{timestamp()}][FATAL] {exc}", file=sys.stderr)
        raise SystemExit(1)
    #endtry