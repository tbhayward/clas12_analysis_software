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

External model tools:
  * BH is evaluated by running dvcsgen. The directory can be provided with
    --dvcsgen-dir or the DVCSGEN_PATH environment variable.
  * KM15 is evaluated by running km15_cli.py. The path can be provided with
    --km15-cli or the KM15_CLI environment variable.
  * On Unix-like systems, KM15 is launched with the Python executable in
    --py-km15 or the PY_KM15 environment variable, with PYTHONPATH unset.

The script is intentionally defensive about column names, because the pass-2
CSV columns have evolved during the analysis. Run with --print-columns to see
exactly what was detected.
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
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Default external paths.
#
# The script is normally run from:
#     dvcs_cross_section/external_scripts/
# with km15_cli.py one directory above it.  Do not hard-code the KM15 Python
# virtualenv here; if PY_KM15/--py-km15 is not set, the current Python
# executable is used.
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
    cache_path: Path
    use_cache: bool
    allow_missing_models: bool


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[pass2-vs-pass1] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[pass2-vs-pass1][warn] {msg}", file=sys.stderr, flush=True)


def die(msg: str) -> None:
    raise RuntimeError(msg)


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
def detect_common_columns(fieldnames: Sequence[str]) -> Dict[str, str]:
    cols: Dict[str, str] = {}

    cols["valid"] = find_column(fieldnames, ["valid bin", "valid", "is valid"], required=False) or ""
    cols["bin_index"] = find_column(fieldnames, ["bin index", "bin", "idx", ""], required=False) or ""

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
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            die(f"Pass-2 CSV appears empty: {path}")
        #endif

        fieldnames = reader.fieldnames

        if print_columns:
            log("Pass-2 columns:\n" + "\n".join(f"  [{i}] {name}" for i, name in enumerate(fieldnames)))
        #endif

        cols = detect_common_columns(fieldnames)
        xs_col = xs_column or find_column(fieldnames, PASS2_XS_CANDIDATES)

        log(f"Pass-2 cross-section column: {xs_col}")

        out: Dict[BinKey, List[DataPoint]] = {}
        kept = 0
        skipped = 0

        for row in reader:
            if not row_is_valid(row, cols["valid"]):
                continue
            #endif

            key = make_bin_key(row, cols)
            phi = row_phi(row, cols)
            xs, stat, _ignored_syst, ok = parse_tuple_cell(row.get(xs_col, ""))

            if not ok or not finite_positive(xs):
                skipped += 1
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

    log(f"Pass-2 rows kept: {kept}; skipped non-positive/unparseable xs rows: {skipped}")
    return out


def load_pass1_csv(path: Path, print_columns: bool = False) -> Dict[BinKey, List[DataPoint]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            die(f"Pass-1 CSV appears empty: {path}")
        #endif

        fieldnames = reader.fieldnames

        if print_columns:
            log("Pass-1 columns:\n" + "\n".join(f"  [{i}] {name}" for i, name in enumerate(fieldnames)))
        #endif

        cols = detect_common_columns(fieldnames)

        xs_col = find_column(fieldnames, [PASS1_XS_COL])
        stat_col = find_column(fieldnames, [PASS1_STAT_COL])
        syst_up_col = find_column(fieldnames, [PASS1_SYST_UP_COL])
        syst_dn_col = find_column(fieldnames, [PASS1_SYST_DN_COL])

        out: Dict[BinKey, List[DataPoint]] = {}
        kept = 0
        skipped = 0

        for row in reader:
            if not row_is_valid(row, cols["valid"]):
                continue
            #endif

            key = make_bin_key(row, cols)
            phi = row_phi(row, cols)

            xs = to_float(row.get(xs_col, ""))
            stat = abs(to_float(row.get(stat_col, "")))
            syst_up = abs(to_float(row.get(syst_up_col, "")))
            syst_dn = abs(to_float(row.get(syst_dn_col, "")))

            if not finite_positive(xs):
                skipped += 1
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

    log(f"Pass-1 rows kept: {kept}; skipped non-positive xs rows: {skipped}")
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
    if cfg.py_km15:
        python_exe = cfg.py_km15
    else:
        python_exe = sys.executable
    #endif

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
        return {}
    #endif

    try:
        with path.open("r") as handle:
            raw = json.load(handle)
        #endwith

        if isinstance(raw, dict):
            return raw
        #endif
    except Exception as exc:
        warn(f"Could not read model cache {path}: {exc}")
    #endtry

    return {}


def save_model_cache(path: Path, cache: Dict[str, Dict[str, List[float]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w") as handle:
        json.dump(cache, handle, indent=2, sort_keys=True)
    #endwith


def compute_model_curves(key: BinKey, cfg: ModelConfig, cache: Dict[str, Dict[str, List[float]]]) -> Optional[ModelCurves]:
    ckey = model_cache_key(key, cfg)

    if cfg.use_cache and ckey in cache:
        cached = cache[ckey]
        return ModelCurves(phi=list(cached["phi"]), bh=list(cached["bh"]), km15=list(cached["km15"]))
    #endif

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

    try:
        for phi in phi_grid:
            bh_values.append(bh_xs(xb, q2, t_abs, phi, cfg))
            km15_values.append(km15_xs(xb, q2, t_abs, phi, cfg))
        #endfor
    except Exception as exc:
        message = (
            f"Model calculation failed for xB=[{key.xb_min}, {key.xb_max}], "
            f"Q2=[{key.q2_min}, {key.q2_max}], |t|=[{key.t_min}, {key.t_max}]: {exc}"
        )

        if cfg.allow_missing_models:
            warn(message)
            return None
        #endif

        raise RuntimeError(message) from exc
    #endtry

    curves = ModelCurves(phi=phi_grid, bh=bh_values, km15=km15_values)
    cache[ckey] = {"phi": curves.phi, "bh": curves.bh, "km15": curves.km15}

    return curves


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
        ymin = min(v for v in vals if v > 0.0)
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
    fig, ax = plt.subplots(figsize=(8.4, 6.0))

    if curves is not None:
        ax.plot(curves.phi, curves.bh, linewidth=2.0, linestyle="-", label="BH")
        ax.plot(curves.phi, curves.km15, linewidth=2.0, linestyle="--", label="KM15")
    #endif

    if panel.pass1:
        x = [p.phi for p in panel.pass1]
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
        x = [p.phi for p in panel.pass2]
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
            label=f"{pass2_label}: stat ⊕ 18% syst",
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
    candidates: List[Path] = []

    if user_value:
        candidates.append(Path(user_value).expanduser())
    #endif

    candidates.append(ANALYSIS_DIR / "km15_cli.py")
    candidates.append(Path.cwd() / "km15_cli.py")

    for cand in candidates:
        if cand.exists() and cand.is_file():
            return str(cand.resolve())
        #endif
    #endfor

    tried = "\n".join(f"  - {c}" for c in candidates)
    die(f"Could not find km15_cli.py. Tried:\n{tried}\nUse --km15-cli /full/path/to/km15_cli.py.")


def resolve_python_executable(user_value: str) -> str:
    """
    Resolve the Python executable used for KM15.

    Important: the previous version defaulted to a hard-coded path under
    /u/home/thayward/venvs/km15_py312/bin/python3. That is brittle and caused
    the reported FileNotFoundError. Now an explicit --py-km15/PY_KM15 is only
    used if it exists; otherwise the current Python executable is used.
    """
    text = str(user_value or "").strip()

    if text:
        expanded = Path(text).expanduser()

        if expanded.exists():
            return str(expanded.resolve())
        #endif

        found = shutil.which(text)

        if found:
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
            f"Model calculation will fail unless --dvcsgen-dir is corrected "
            f"or --allow-missing-models is used."
        )
    #endif

    return str(path)


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def build_panels(pass2: Dict[BinKey, List[DataPoint]], pass1: Dict[BinKey, List[DataPoint]]) -> List[PanelData]:
    keys = sorted(set(pass2.keys()) | set(pass1.keys()), key=key_sort_tuple)
    panels: List[PanelData] = []

    for key in keys:
        panels.append(PanelData(key=key, pass2=pass2.get(key, []), pass1=pass1.get(key, [])))
    #endfor

    return panels


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
    parser.add_argument("--phi-dense", type=int, default=181, help="Number of phi points for theory curves.")
    parser.add_argument("--model-cache", type=Path, default=None, help="Optional JSON cache path for model curves.")
    parser.add_argument("--no-cache", action="store_true", help="Do not read existing cached model curves.")
    parser.add_argument("--allow-missing-models", action="store_true", help="Still make data plots if BH/KM15 commands fail.")
    parser.add_argument("--skip-models", action="store_true", help="Do not compute or draw BH/KM15 curves.")
    parser.add_argument("--linear-y", action="store_true", help="Use linear y scale instead of the default log y scale.")
    parser.add_argument("--print-columns", action="store_true", help="Print detected CSV columns before plotting.")

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.pass2_csv.exists():
        die(f"Pass-2 CSV does not exist: {args.pass2_csv}")
    #endif

    if not args.pass1_csv.exists():
        die(f"Pass-1 CSV does not exist: {args.pass1_csv}")
    #endif

    args.output_dir.mkdir(parents=True, exist_ok=True)

    cache_path = args.model_cache or (args.output_dir / "model_curve_cache.json")

    args.dvcsgen_dir = validate_dvcsgen_dir(args.dvcsgen_dir)
    args.km15_cli = resolve_km15_cli(args.km15_cli)
    args.py_km15 = resolve_python_executable(args.py_km15)

    pass2 = load_pass2_csv(args.pass2_csv, args.pass2_xs_column, print_columns=args.print_columns)
    pass1 = load_pass1_csv(args.pass1_csv, print_columns=args.print_columns)
    panels = build_panels(pass2, pass1)

    log(f"Unique (xB,Q2,|t|) panels to write: {len(panels)}")

    model_cfg = ModelConfig(
        e_beam=args.e_beam,
        dvcsgen_dir=args.dvcsgen_dir,
        km15_cli=args.km15_cli,
        py_km15=args.py_km15,
        phi_dense=max(2, args.phi_dense),
        cache_path=cache_path,
        use_cache=not args.no_cache,
        allow_missing_models=args.allow_missing_models,
    )

    cache: Dict[str, Dict[str, List[float]]] = {}

    if not args.skip_models:
        cache = load_model_cache(cache_path) if model_cfg.use_cache else {}
        log(f"Model settings: E={args.e_beam:g} GeV, dvcsgen_dir={args.dvcsgen_dir}, km15_cli={args.km15_cli}")
        log(f"KM15 Python: {args.py_km15}")
    #endif

    manifest: List[Dict[str, object]] = []

    for i, panel in enumerate(panels, start=1):
        curves: Optional[ModelCurves] = None

        if not args.skip_models:
            curves = compute_model_curves(panel.key, model_cfg, cache)
        #endif

        filename = panel_filename(panel.key, i)
        outpath = args.output_dir / filename

        draw_panel(
            panel=panel,
            curves=curves,
            output_path=outpath,
            pass2_label=args.pass2_label,
            pass1_label=args.pass1_label,
            logy=not args.linear_y,
        )

        manifest.append(
            {
                "file": filename,
                "xB": [panel.key.xb_min, panel.key.xb_max],
                "Q2": [panel.key.q2_min, panel.key.q2_max],
                "t_abs": [panel.key.t_min, panel.key.t_max],
                "n_pass2_points": len(panel.pass2),
                "n_pass1_points": len(panel.pass1),
                "has_models": curves is not None,
            }
        )

        log(f"[{i}/{len(panels)}] wrote {outpath}")
    #endfor

    if not args.skip_models:
        save_model_cache(cache_path, cache)
        log(f"Model cache written: {cache_path}")
    #endif

    manifest_path = args.output_dir / "manifest.json"

    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
    #endwith

    log(f"Manifest written: {manifest_path}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[pass2-vs-pass1][FATAL] {exc}", file=sys.stderr)
        raise SystemExit(1)
    #endtry