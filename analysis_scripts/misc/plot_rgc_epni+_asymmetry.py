#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ep -> en pi+ asymmetries versus -t' in several xB bins, for three run periods
and for the combined (inverse-variance weighted) file.

This version:
  - Produces a 1x3 plot ONLY (LU, UL, LL).
  - Completely omits unpolarized modulations (UUcos, UUcos2) from plots.
  - Prints LaTeX tables with \\renewcommand{\\arraystretch}{1.40} applied
    within each table environment.

Corrections / shifts (IMPORTANT):
  1) Radiative (ISR/FSR) correction convention (as before):
       A_Born = A_baseline - Delta_rad
     i.e. whenever --rad is used, central values are shifted as:
       y <- y - Delta_rad
     and statistical errors remain unchanged.

  2) Resolution (regular) migration correction (signed diffs):
     If --migration is provided and the file contains signed diffs (delta),
     we ADD that signed diff AFTER the radiative shift:
       y <- (y - Delta_rad) + diff_mig

  3) Fermi-motion migration correction (signed diffs):
     If --fermi_migration is provided (signed diffs),
     we ADD that signed diff as well:
       y <- (y - Delta_rad) + diff_mig + diff_fermi

Systematic uncertainties (IMPORTANT - UPDATED BEHAVIOR):
  - Radiative systematic (sigma_rad) comes from the --rad summary file.
  - Resolution-migration systematic (sigma_resmig) is defined as:
        sigma_resmig = |diff_mig|
  - Fermi-migration systematic (sigma_fermi) is defined as:
        sigma_fermi = |diff_fermi|

  - Total systematic per bin (drawn as solid gray rectangles) is:
        sigma_tot = sqrt( sigma_rad^2 + sigma_resmig^2 + sigma_fermi^2 )

  Notes:
    * If --fermi_migration is omitted, sigma_fermi is omitted from the quadrature sum.
    * If --migration is omitted, sigma_resmig is omitted from the quadrature sum.
    * LaTeX tables are still printed only when --rad AND --migration are provided.

LaTeX tables (stdout) for the 4 xB slices:
  enpiLowxBGE, enpiMidLowxBGE, enpiMidHighxBGE, enpiHighxBGE

Tables are printed only when --rad AND --migration are provided (same trigger as before).
If --fermi_migration is also provided, it is included in BOTH:
  - the central-value shift, and
  - the total systematic quadrature sum.

Tables use:
  central = (fit value) - (radiative Delta) + (res-mig diff) + (fermi-mig diff)
  stat    = (fit statistical uncertainty)
  syst    = sqrt( sigma_rad^2 + sigma_resmig^2 + sigma_fermi^2 )
where:
  sigma_resmig = |diff_mig|
  sigma_fermi  = |diff_fermi|

Rows are ordered in increasing -t' (i.e. increasing x axis as plotted).

Usage:
  python plot_rgc_enpi+_asymmetry.py Su22.txt Fa22.txt Sp23.txt Combined.txt --rad output/enpi+/ISR_FSR_delta_summary.txt --migration migration_diff_summary.txt --fermi_migration fermi_migration_diff_summary.txt

Input format expectations:
  1) Asymmetry fit-result files:
     contain blocks like:
       enpiLowxBGEchi2FitsALUsinphi = {{mean_tprime, value, error}, ...};
     (x is plotted as -t').

  2) Radiative delta summary file (--rad):
     Expected "block" format produced by Script 1:
       Bin: <bin_tag> ...
       Series: <KEY> | <label>
         -t'    Delta    sigma
         0.10   ...      ...
     KEY must be one of: ALUsin, AULsin, AULsin2, ALLn0, ALLcos

     NOTE: This script ALSO supports Script-1 summaries that use LaTeX series labels, e.g.:
       Series: $F_{LU}^{\\sin\\phi}/F_{UU}$
     by mapping those to the internal keys above.

  3) Migration diff summary files (--migration, --fermi_migration):
     Supported formats (explicitly detected):
       (A) Block format (preferred):
           Bin: <bin_tag> ...
           Series: <KEY> | <label>
             -t'    diff   [optional third column, ignored by this script]
             0.10   ...
           We interpret the SECOND numeric column as the signed diff to APPLY:
             y <- y + diff
           and the migration systematic uses abs(diff).

       (B) Legacy assignment format:
           NAME = {{mean_tprime, diff}, ...};
           where NAME matches:
             <bin_tag>chi2Fits<suffix>_delta
           and we interpret the second entry as the signed diff to APPLY.

Notes on systematic rectangles:
  - Rectangles span full fixed bins, not centered on data means.
  - Edges (given to you in -t, we use their absolute values on the +(-t') axis):
      * For enpiGE (integrated), 0.10-wide bins with edges:
          0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25
      * For xB sub-bins, 0.20-wide bins with edges:
          0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25
"""

import sys
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# ---------------------------------------------------------------------
# Styling / knobs
# ---------------------------------------------------------------------
COLORS = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green",
    "Combined": "black",
}

XB_COLORS = {
    "enpiLowxBGE":     "tab:blue",
    "enpiMidLowxBGE":  "tab:orange",
    "enpiMidHighxBGE": "tab:green",
    "enpiHighxBGE":    "tab:red",
}

XB_MARKERS = {
    "enpiLowxBGE":     "o",
    "enpiMidLowxBGE":  "s",
    "enpiMidHighxBGE": "^",
    "enpiHighxBGE":    "D",
}

MARKER = "o"
CAPSIZE = 3
MS = 5.0

XLIM_T = (0.0, 1.30)
X_LABEL = r"$-t'\ (\mathrm{GeV}^{2})$"

# Updated y ranges per your request:
YLIM_LU = (-0.4, 0.4)    # BSA
YLIM_UL = (-0.4, 0.4)    # TSA
YLIM_LL = (-1.0, 1.0)    # DSA

XB_BINS = {
    "enpiLowxBGE":     r"$0.10 < x_{B} < 0.25$",
    "enpiMidLowxBGE":  r"$0.25 < x_{B} < 0.35$",
    "enpiMidHighxBGE": r"$0.35 < x_{B} < 0.45$",
    "enpiHighxBGE":    r"$0.45 < x_{B} < 0.60$",
    "enpiGE":          r"$0.10 < x_{B} < 0.60$"
}

BIN_ORDER = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE", "enpiGE"]

TOP_PAD_PER_BIN = 0.92

OUT_PREFIX = "enpi"

# Systematic rectangle bin edges (x axis is +(-t')):
EDGES_SYS_GE = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25], dtype=float)
EDGES_SYS_XB = np.array([0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25], dtype=float)

# ---------------------------------------------------------------------
# Parsing helpers (fit-result .txt files)
# ---------------------------------------------------------------------
_assign_re = re.compile(r'([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};', re.DOTALL)
_triple_re = re.compile(r'\{([^{}]+)\}')

def parse_asym_file(path):
    """Parse NAME = {{mean,val,err}, ...}; blocks into dict[name] -> np.array(N,3)."""
    with open(path, "r") as f:
        text = f.read()
    #endif
    out = {}
    for m in _assign_re.finditer(text):
        name = m.group(1)
        content = m.group(2)
        triples = []
        for t in _triple_re.findall(content):
            parts = [p.strip() for p in t.split(",")]
            if len(parts) != 3:
                continue
            #endif
            try:
                triples.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
            #endif
        #endfor
        if triples:
            out[name] = np.array(triples, dtype=float)
        #endif
    #endfor
    return out

def get_series(dct, key, negate_x=True, sort_x=True):
    """
    Return dict(x,y,yerr) if key exists with finite values and positive errors; else None.
    Negates x if negate_x=True (to convert t' -> -t' for plotting). Optionally sort by x.
    """
    if key not in dct:
        return None
    #endif
    arr = np.array(dct[key], dtype=float)
    if arr.size == 0:
        return None
    #endif
    x_raw, y, e = arr[:, 0], arr[:, 1], arr[:, 2]
    mask = np.isfinite(x_raw) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    if not np.any(mask):
        return None
    #endif
    x = -x_raw[mask] if negate_x else x_raw[mask]
    y = y[mask]
    e = e[mask]
    if sort_x and x.size > 1:
        order = np.argsort(x)
        x, y, e = x[order], y[order], e[order]
    #endif
    return {"x": x, "y": y, "yerr": e}

def build_period_dict(parsed, bin_prefix):
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}"
    return {
        "ALUsin":  get_series(parsed, k("ALUsinphi")),
        "AULsin":  get_series(parsed, k("AULsinphi")),
        "AULsin2": get_series(parsed, k("AULsin2phi")),
        "ALLn0":   get_series(parsed, k("ALL")),
        "ALLcos":  get_series(parsed, k("ALLcosphi")),
        # UU intentionally omitted from plotting in this version
    }

def detect_available_bins(*dicts):
    available = []
    for bin_tag in BIN_ORDER:
        key_probe = f"{bin_tag}chi2FitsALUsinphi"
        present = any((d is not None and key_probe in d) for d in dicts)
        if present:
            available.append(bin_tag)
        #endif
    #endfor
    return available

# ---------------------------------------------------------------------
# Radiative-correction summary parser (signed Delta summary txt)
# Block format: Bin: <bin> ; Series: <KEY or label> ; numeric rows
# ---------------------------------------------------------------------
_ALLOWED_KEYS = ["ALUsin", "AULsin", "AULsin2", "ALLn0", "ALLcos"]

def _norm_series_token(s):
    """
    Normalize a Series token for robust matching:
      - strip whitespace
      - remove '$'
      - remove spaces
      - keep backslashes/braces/underscores/carets (LaTeX structure)
      - lowercase
    """
    if s is None:
        return ""
    #endif
    out = s.strip()
    out = out.replace("$", "")
    out = out.replace(" ", "")
    out = out.lower()
    return out

# Map LaTeX-style series labels to internal keys.
_RAD_LABEL_TO_KEY = {
    _norm_series_token(r"F_{LU}^{\sin\phi}/F_{UU}"):      "ALUsin",
    _norm_series_token(r"F_{UL}^{\sin\phi}/F_{UU}"):      "AULsin",
    _norm_series_token(r"F_{UL}^{\sin2\phi}/F_{UU}"):     "AULsin2",
    _norm_series_token(r"F_{LL}/F_{UU}"):                 "ALLn0",
    _norm_series_token(r"F_{LL}^{\cos\phi}/F_{UU}"):      "ALLcos",
    _norm_series_token(r"DilutionD_f"):                   "Df",
    _norm_series_token(r"DilutionD_{f}"):                 "Df",
}

def _parse_series_key_from_line(line):
    """
    Accept either:
      Series: <KEY> | <label>
    or
      Series: <label-only>
    and return a canonical internal key if possible.
    """
    s = line.strip()
    if not s.startswith("Series:"):
        return None
    #endif
    rest = s[len("Series:"):].strip()

    token = rest.split("|", 1)[0].strip() if ("|" in rest) else rest.strip()
    if not token:
        return None
    #endif

    if token in _ALLOWED_KEYS:
        return token
    #endif

    norm = _norm_series_token(token)
    mapped = _RAD_LABEL_TO_KEY.get(norm)
    if mapped in _ALLOWED_KEYS:
        return mapped
    #endif

    if mapped == "Df":
        return "Df"
    #endif

    return token

def parse_rad_summary(rad_path):
    """
    Returns:
      rad[bin_tag][series_key] = {"x": array(-t'), "delta": array, "sigma": array}
    where series_key is one of _ALLOWED_KEYS (and optionally "Df").

    IMPORTANT: This file's x is assumed already in +(-t') units (positive numbers).
    Pairing to data is by index when lengths differ.
    """
    rad = {}
    cur_bin = None

    with open(rad_path, "r") as f:
        lines = f.readlines()
    #endif

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("Bin:"):
            parts = line.split()
            if len(parts) >= 2:
                cur_bin = parts[1]
                rad.setdefault(cur_bin, {})
            else:
                raise RuntimeError(f"FATAL: malformed Bin line in rad file: '{line}'")
            #endif
            i += 1
            continue
        #endif

        if line.startswith("Series:"):
            if cur_bin is None:
                raise RuntimeError("FATAL: encountered Series before any Bin in rad file.")
            #endif

            key_raw = _parse_series_key_from_line(line)
            if key_raw is None:
                raise RuntimeError(f"FATAL: malformed Series line in rad file: '{line}'")
            #endif

            use_key = (key_raw in _ALLOWED_KEYS) or (key_raw == "Df")
            if not use_key:
                print(f"[INFO] rad: ignoring unknown series token '{key_raw}' in bin '{cur_bin}'.")
            #endif

            i += 1
            if i < len(lines):
                hdr = lines[i].strip()
                if (("-t'" in hdr) or ("Delta" in hdr) or ("delta" in hdr) or ("sigma" in hdr)):
                    i += 1
                #endif
            #endif

            xs, deltas, sigmas = [], [], []
            while i < len(lines):
                row = lines[i].strip()
                if (not row) or row.startswith("Series:") or row.startswith("Bin:") or (set(row) <= set("=-")):
                    break
                #endif
                parts = row.split()
                if len(parts) < 3:
                    raise RuntimeError(f"FATAL: rad row must have 3 numeric columns (x delta sigma). Got: '{row}'")
                #endif
                try:
                    xval = float(parts[0])
                    dval = float(parts[1])
                    sval = float(parts[2])
                except Exception:
                    raise RuntimeError(f"FATAL: non-numeric rad row: '{row}'")
                #endif
                xs.append(xval)
                deltas.append(dval)
                sigmas.append(sval)
                i += 1
            #endfor

            if use_key and xs:
                rad[cur_bin][key_raw] = {
                    "x": np.array(xs, dtype=float),
                    "delta": np.array(deltas, dtype=float),
                    "sigma": np.array(sigmas, dtype=float),
                }
            #endif
            continue
        #endif

        i += 1
    #endfor

    return rad

# ---------------------------------------------------------------------
# Migration diff parsers
#   - Block format: Bin/Series blocks, numeric rows
#   - Legacy assignment format: NAME = {{x, diff}, ...};
# We store BOTH:
#   "delta" = signed diff to APPLY (y <- y + delta)
#   "sigma" = abs(delta)  (magnitude used for sys)
# ---------------------------------------------------------------------
_pair_re = re.compile(r'\{([^{}]+)\}')

def parse_migration_file_legacy(path):
    """Parse legacy NAME = {{x, diff}, ...}; blocks into dict[name] -> np.array(N,2)."""
    with open(path, "r") as f:
        text = f.read()
    #endif
    out = {}
    for m in _assign_re.finditer(text):
        name = m.group(1)
        content = m.group(2)
        pairs = []
        for t in _pair_re.findall(content):
            parts = [p.strip() for p in t.split(",")]
            if len(parts) != 2:
                continue
            #endif
            try:
                pairs.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
            #endif
        #endfor
        if pairs:
            out[name] = np.array(pairs, dtype=float)
        #endif
    #endfor
    return out

def get_mig_delta_series_legacy(dct, key, negate_x=True, sort_x=True):
    """
    Legacy migration *_delta list: (mean_tprime, diff)
    Output:
      x     = -mean_tprime (if negate_x)
      delta = diff (signed, to APPLY)
      sigma = abs(diff) (magnitude used for sys)
    """
    if key not in dct:
        return None
    #endif
    arr = np.array(dct[key], dtype=float)
    if arr.size == 0:
        return None
    #endif
    x_raw, delta = arr[:, 0], arr[:, 1]
    mask = np.isfinite(x_raw) & np.isfinite(delta)
    if not np.any(mask):
        return None
    #endif
    x = -x_raw[mask] if negate_x else x_raw[mask]
    d = delta[mask]
    s = np.abs(d)
    if sort_x and x.size > 1:
        order = np.argsort(x)
        x, d, s = x[order], d[order], s[order]
    #endif
    return {"x": x, "delta": d, "sigma": s}

def build_migration_dict_legacy(mig_parsed, bin_prefix):
    if mig_parsed is None:
        return None
    #endif
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}_delta"
    return {
        "ALUsin":  get_mig_delta_series_legacy(mig_parsed, k("ALUsinphi")),
        "AULsin":  get_mig_delta_series_legacy(mig_parsed, k("AULsinphi")),
        "AULsin2": get_mig_delta_series_legacy(mig_parsed, k("AULsin2phi")),
        "ALLn0":   get_mig_delta_series_legacy(mig_parsed, k("ALL")),
        "ALLcos":  get_mig_delta_series_legacy(mig_parsed, k("ALLcosphi")),
    }

def parse_migration_summary_block(path):
    """
    Parse block-format migration diff summary into:
      mig[bin_tag][series_key] = {"x": array(-t'), "delta": array, "sigma": array}

    Numeric rows:
      - If 2+ columns: we interpret:
          col0 = x
          col1 = diff (signed, to APPLY)
        Any additional columns are ignored by this script.
      - sigma is defined here as abs(diff), since your definition of the migration systematic
        is based on the size of the applied shift.
    """
    mig = {}
    cur_bin = None

    with open(path, "r") as f:
        lines = f.readlines()
    #endif

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("Bin:"):
            parts = line.split()
            if len(parts) >= 2:
                cur_bin = parts[1]
                mig.setdefault(cur_bin, {})
            else:
                raise RuntimeError(f"FATAL: malformed Bin line in migration file: '{line}'")
            #endif
            i += 1
            continue
        #endif

        if line.startswith("Series:"):
            if cur_bin is None:
                raise RuntimeError("FATAL: encountered Series before any Bin in migration file.")
            #endif

            key = _parse_series_key_from_line(line)
            if key is None:
                raise RuntimeError(f"FATAL: malformed Series line in migration file: '{line}'")
            #endif

            use_key = (key in _ALLOWED_KEYS)
            if not use_key:
                print(f"[INFO] migration: ignoring unknown series key '{key}' in bin '{cur_bin}'.")
            #endif

            i += 1
            if i < len(lines):
                hdr = lines[i].strip()
                if (("-t'" in hdr) or ("delta" in hdr) or ("Delta" in hdr) or ("diff" in hdr) or ("sigma" in hdr)):
                    i += 1
                #endif
            #endif

            xs, deltas, sigmas = [], [], []
            while i < len(lines):
                row = lines[i].strip()
                if (not row) or row.startswith("Series:") or row.startswith("Bin:") or (set(row) <= set("=-")):
                    break
                #endif
                parts = row.split()
                if len(parts) < 2:
                    raise RuntimeError(f"FATAL: migration row must have at least 2 numeric columns. Got: '{row}'")
                #endif
                try:
                    xval = float(parts[0])
                    dval = float(parts[1])
                except Exception:
                    raise RuntimeError(f"FATAL: non-numeric migration row: '{row}'")
                #endif
                xs.append(xval)
                deltas.append(dval)
                sigmas.append(abs(dval))
                i += 1
            #endfor

            if use_key and xs:
                mig[cur_bin][key] = {
                    "x": np.array(xs, dtype=float),
                    "delta": np.array(deltas, dtype=float),
                    "sigma": np.array(sigmas, dtype=float),
                }
            #endif
            continue
        #endif

        i += 1
    #endfor

    return mig

def _file_has_bin_blocks(path):
    with open(path, "r") as f:
        for line in f:
            if line.strip().startswith("Bin:"):
                return True
            #endif
        #endfor
    #endif
    return False

def parse_migration_any_format(path):
    """
    Wrapper: detect block vs legacy and return:
      mig_all[bin_tag][key] = {"x","delta","sigma"}
    """
    if not os.path.isfile(path):
        raise RuntimeError(f"FATAL: migration file not found: {path}")
    #endif

    if _file_has_bin_blocks(path):
        print(f"[INFO] migration: detected block-format summary (Bin:/Series:) in: {path}")
        return parse_migration_summary_block(path)
    #endif

    print(f"[INFO] migration: detected legacy NAME={{x,diff}} format in: {path}")
    mig_parsed = parse_migration_file_legacy(path)
    mig_all = {}
    for bin_tag in BIN_ORDER:
        mig_all[bin_tag] = build_migration_dict_legacy(mig_parsed, bin_tag) or {}
    #endfor
    return mig_all

# ---------------------------------------------------------------------
# Titles / bin edges
# ---------------------------------------------------------------------
def make_title(xb_label):
    return r"$ep \rightarrow en\pi^{+}$ -- " + xb_label

def _sys_edges_for_bin(bin_tag):
    return EDGES_SYS_GE if bin_tag == "enpiGE" else EDGES_SYS_XB

# ---------------------------------------------------------------------
# Systematics drawing (solid gray only) - UPDATED to quadrature sum of:
#   sigma_rad, sigma_resmig, sigma_fermi
# ---------------------------------------------------------------------
def _quadrature_sum_arrays(arr_list, nbins):
    """
    Given a list of arrays (some may be None), truncate consistently and return:
      tot[i] = sqrt(sum_j arr_j[i]^2)
    Returns None if no usable arrays exist.
    """
    usable = []
    for a in arr_list:
        if a is None:
            continue
        #endif
        aa = np.asarray(a, dtype=float)
        if aa.size <= 0:
            continue
        #endif
        usable.append(aa)
    #endfor

    if not usable:
        return None
    #endif

    n = nbins
    for aa in usable:
        n = min(n, int(aa.size))
    #endfor

    if n <= 0:
        return None
    #endif

    acc = np.zeros(n, dtype=float)
    for aa in usable:
        acc += aa[:n] * aa[:n]
    #endfor

    return np.sqrt(acc)

def _draw_total_sys_band_solid(ax, edges, sigma_rad=None, sigma_resmig=None, sigma_fermi=None, zorder=1):
    nbins = max(0, len(edges) - 1)
    if nbins <= 0:
        return
    #endif

    tot = _quadrature_sum_arrays([sigma_rad, sigma_resmig, sigma_fermi], nbins=nbins)
    if tot is None:
        return
    #endif

    n = min(nbins, int(tot.size))
    if n <= 0:
        return
    #endif

    for i in range(n):
        sv = tot[i]
        if not np.isfinite(sv):
            continue
        #endif
        x0 = edges[i]
        x1 = edges[i + 1]
        if (not np.isfinite(x0)) or (not np.isfinite(x1)) or (x1 <= x0):
            continue
        #endif
        rect = Rectangle(
            (x0, 0.0),
            x1 - x0,
            sv,
            facecolor="0.7",
            edgecolor="none",
            linewidth=0.0,
            alpha=0.25,
            zorder=zorder
        )
        ax.add_patch(rect)
    #endfor

def _warn_len_mismatch(tag, lab, a_len, b_len, n_use):
    if a_len != b_len:
        print(f"[WARN] {tag}: length mismatch for {lab}: {a_len} vs {b_len}; pairing by index up to {n_use}.")
    #endif

def _get_delta_array(entry):
    if entry is None:
        return None
    #endif
    d = entry.get("delta", None)
    if d is None:
        return None
    #endif
    return np.asarray(d, dtype=float)

def _get_sigma_array(entry, prefer_abs_delta=True):
    """
    For migration, sigma is defined as abs(delta) (size of applied shift).
    For radiation, sigma comes from the file (prefer_abs_delta=False).
    """
    if entry is None:
        return None
    #endif

    if prefer_abs_delta and ("delta" in entry):
        d = np.asarray(entry.get("delta", []), dtype=float)
        if d.size > 0:
            return np.abs(d)
        #endif
    #endif

    s = np.asarray(entry.get("sigma", []), dtype=float)
    if s.size > 0:
        return np.abs(s)
    #endif
    return None

def _apply_total_shift(series, rad_entry=None, mig_entry=None, fermi_entry=None,
                      with_rad=False, with_mig=False, with_fermi=False, label=""):
    """
    Apply shifts in the required order:
      y <- y - Delta_rad
      y <- y + diff_mig
      y <- y + diff_fermi

    Pairing is by index if arrays differ in length.
    """
    if series is None:
        return None
    #endif

    x = series["x"]
    y = series["y"].copy()
    yerr = series["yerr"]

    if with_rad and (rad_entry is not None):
        dr = _get_delta_array(rad_entry)
        if dr is not None:
            n = min(y.size, dr.size)
            if n > 0:
                if dr.size != y.size:
                    _warn_len_mismatch("rad_shift", label, dr.size, y.size, n)
                #endif
                y[:n] = y[:n] - dr[:n]
            #endif
        #endif
    #endif

    if with_mig and (mig_entry is not None):
        dm = _get_delta_array(mig_entry)
        if dm is not None:
            n = min(y.size, dm.size)
            if n > 0:
                if dm.size != y.size:
                    _warn_len_mismatch("mig_shift", label, dm.size, y.size, n)
                #endif
                y[:n] = y[:n] + dm[:n]
            #endif
        #endif
    #endif

    if with_fermi and (fermi_entry is not None):
        df = _get_delta_array(fermi_entry)
        if df is not None:
            n = min(y.size, df.size)
            if n > 0:
                if df.size != y.size:
                    _warn_len_mismatch("fermi_shift", label, df.size, y.size, n)
                #endif
                y[:n] = y[:n] + df[:n]
            #endif
        #endif
    #endif

    return x, y, yerr

# ---------------------------------------------------------------------
# Legends
# ---------------------------------------------------------------------
def _legend_run_period(ax, labels, where="lower right"):
    handles = [Line2D([0], [0], marker=MARKER, color=COLORS[L], linestyle="", label=L) for L in labels]
    leg = ax.legend(handles=handles, title="Run Period", frameon=True, edgecolor="black",
                    loc=where, fontsize=11, title_fontsize=12)
    leg.get_frame().set_alpha(0.9)

def _legend_harmonic(ax, labels=("n=1", "n=2"), where="lower left"):
    open_marker = Line2D([0], [0], marker=MARKER, mfc="none", mec="black", linestyle="", label=labels[0])
    filled_marker = Line2D([0], [0], marker=MARKER, color="black", linestyle="", label=labels[1])
    leg = ax.legend(handles=[open_marker, filled_marker], title="Harmonic",
                    frameon=True, edgecolor="black", loc=where,
                    bbox_to_anchor=(0.02, 0.02), fontsize=11, title_fontsize=12)
    ax.add_artist(leg)

def _suffix(with_rad, with_mig, with_fermi):
    """
    Keep filenames stable as much as possible:
      - previous "_withRadMig" remains the suffix whenever rad AND (any migration) are enabled.
      - "_withMig" when migrations (any) but no rad.
      - "_withRad" when only rad.
    """
    any_mig = with_mig or with_fermi
    if with_rad and any_mig:
        return "_withRadMig"
    #endif
    if with_rad:
        return "_withRad"
    #endif
    if any_mig:
        return "_withMig"
    #endif
    return ""

# ---------------------------------------------------------------------
# Plotting (1x3 canvases: LU, UL, LL only)
# ---------------------------------------------------------------------
def _plot_panel_sets_1x3(axLU, axUL, axLL, pdata_by_label,
                        rad_for_bin=None,
                        mig_for_bin=None,
                        fermi_for_bin=None,
                        with_rad=False,
                        with_mig=False,
                        with_fermi=False,
                        bin_tag_for_width="enpiGE"):
    edges = _sys_edges_for_bin(bin_tag_for_width)

    def R(key):
        if rad_for_bin is None:
            return None
        #endif
        return rad_for_bin.get(key)

    def M(key):
        if mig_for_bin is None:
            return None
        #endif
        return mig_for_bin.get(key)

    def F(key):
        if fermi_for_bin is None:
            return None
        #endif
        return fermi_for_bin.get(key)

    def rad_sigma(key, n_target):
        e = R(key)
        if e is None:
            return None
        #endif
        s = _get_sigma_array(e, prefer_abs_delta=False)
        if s is None:
            return None
        #endif
        n = min(n_target, int(s.size))
        if n <= 0:
            return None
        #endif
        return s[:n]

    def mig_sigma_res(key, n_target):
        e = M(key)
        if e is None:
            return None
        #endif
        s = _get_sigma_array(e, prefer_abs_delta=True)
        if s is None:
            return None
        #endif
        n = min(n_target, int(s.size))
        if n <= 0:
            return None
        #endif
        return s[:n]

    def mig_sigma_fermi(key, n_target):
        e = F(key)
        if e is None:
            return None
        #endif
        s = _get_sigma_array(e, prefer_abs_delta=True)
        if s is None:
            return None
        #endif
        n = min(n_target, int(s.size))
        if n <= 0:
            return None
        #endif
        return s[:n]

    # Determine a target length for sys arrays from any available series in pdata_by_label
    def _target_len_for_key(key):
        for _, pdata in pdata_by_label.items():
            s = pdata.get(key)
            if s is not None:
                return int(s["x"].size)
            #endif
        #endfor
        return 0

    # ---------------- BSA (ALU sin phi) ----------------
    nL = _target_len_for_key("ALUsin")
    rs = rad_sigma("ALUsin", nL) if with_rad else None
    ms_res = mig_sigma_res("ALUsin", nL) if with_mig else None
    ms_fermi = mig_sigma_fermi("ALUsin", nL) if with_fermi else None
    _draw_total_sys_band_solid(axLU, edges, sigma_rad=rs, sigma_resmig=ms_res, sigma_fermi=ms_fermi, zorder=1)

    for lab, pdata in pdata_by_label.items():
        s = pdata.get("ALUsin")
        if s is None:
            continue
        #endif
        x, y, ye = _apply_total_shift(
            s,
            rad_entry=R("ALUsin"),
            mig_entry=M("ALUsin"),
            fermi_entry=F("ALUsin"),
            with_rad=with_rad,
            with_mig=with_mig,
            with_fermi=with_fermi,
            label=f"{lab}:ALUsin"
        )
        axLU.errorbar(x, y, yerr=ye, fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                      capsize=CAPSIZE, label=lab, markersize=MS, linestyle="None")
    #endfor
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)
    _legend_run_period(axLU, list(pdata_by_label.keys()), where="lower right")

    # ---------------- TSA (AUL sin phi open, sin2 phi filled) ----------------
    nU1 = _target_len_for_key("AULsin")
    rs = rad_sigma("AULsin", nU1) if with_rad else None
    ms_res = mig_sigma_res("AULsin", nU1) if with_mig else None
    ms_fermi = mig_sigma_fermi("AULsin", nU1) if with_fermi else None
    _draw_total_sys_band_solid(axUL, edges, sigma_rad=rs, sigma_resmig=ms_res, sigma_fermi=ms_fermi, zorder=1)

    for lab, pdata in pdata_by_label.items():
        s1 = pdata.get("AULsin")
        if s1 is not None:
            x1, y1, ye1 = _apply_total_shift(
                s1,
                rad_entry=R("AULsin"),
                mig_entry=M("AULsin"),
                fermi_entry=F("AULsin"),
                with_rad=with_rad,
                with_mig=with_mig,
                with_fermi=with_fermi,
                label=f"{lab}:AULsin"
            )
            axUL.errorbar(x1, y1, yerr=ye1, fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif

        s2 = pdata.get("AULsin2")
        if s2 is not None:
            x2, y2, ye2 = _apply_total_shift(
                s2,
                rad_entry=R("AULsin2"),
                mig_entry=M("AULsin2"),
                fermi_entry=F("AULsin2"),
                with_rad=with_rad,
                with_mig=with_mig,
                with_fermi=with_fermi,
                label=f"{lab}:AULsin2"
            )
            axUL.errorbar(x2, y2, yerr=ye2, fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif
    #endfor

    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axUL, labels=("n=1", "n=2"))
    _legend_run_period(axUL, list(pdata_by_label.keys()), where="lower right")

    # ---------------- DSA (ALL n=0 open, cos phi filled) ----------------
    nD = _target_len_for_key("ALLn0")
    rs = rad_sigma("ALLn0", nD) if with_rad else None
    ms_res = mig_sigma_res("ALLn0", nD) if with_mig else None
    ms_fermi = mig_sigma_fermi("ALLn0", nD) if with_fermi else None
    _draw_total_sys_band_solid(axLL, edges, sigma_rad=rs, sigma_resmig=ms_res, sigma_fermi=ms_fermi, zorder=1)

    for lab, pdata in pdata_by_label.items():
        s0 = pdata.get("ALLn0")
        if s0 is not None:
            x0, y0, ye0 = _apply_total_shift(
                s0,
                rad_entry=R("ALLn0"),
                mig_entry=M("ALLn0"),
                fermi_entry=F("ALLn0"),
                with_rad=with_rad,
                with_mig=with_mig,
                with_fermi=with_fermi,
                label=f"{lab}:ALLn0"
            )
            axLL.errorbar(x0, y0, yerr=ye0, fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif

        s1 = pdata.get("ALLcos")
        if s1 is not None:
            x1, y1, ye1 = _apply_total_shift(
                s1,
                rad_entry=R("ALLcos"),
                mig_entry=M("ALLcos"),
                fermi_entry=F("ALLcos"),
                with_rad=with_rad,
                with_mig=with_mig,
                with_fermi=with_fermi,
                label=f"{lab}:ALLcos"
            )
            axLL.errorbar(x1, y1, yerr=ye1, fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif
    #endfor

    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axLL, labels=("n=0", "n=1"))
    _legend_run_period(axLL, list(pdata_by_label.keys()), where="lower right")

def plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, out_dir):
    plt.figure(figsize=(16, 5))
    axLU = plt.subplot(1, 3, 1)
    axUL = plt.subplot(1, 3, 2)
    axLL = plt.subplot(1, 3, 3)

    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(xb_label), fontsize=16, y=0.98)

    pdata_by_label = {"Su22": p_su22, "Fa22": p_fa22, "Sp23": p_sp23}
    _plot_panel_sets_1x3(axLU, axUL, axLL, pdata_by_label,
                        rad_for_bin=None, mig_for_bin=None, fermi_for_bin=None,
                        with_rad=False, with_mig=False, with_fermi=False,
                        bin_tag_for_width=bin_tag)

    plt.tight_layout(rect=[0, 0, 1, TOP_PAD_PER_BIN])
    os.makedirs(out_dir, exist_ok=True)

    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_{bin_tag}_AllPeriods.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved all-periods figure: {out_path}")

def plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                               with_rad=False, with_mig=False, with_fermi=False,
                               rad_bin=None, mig_bin=None, fermi_bin=None):
    plt.figure(figsize=(16, 5))
    axLU = plt.subplot(1, 3, 1)
    axUL = plt.subplot(1, 3, 2)
    axLL = plt.subplot(1, 3, 3)

    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(xb_label), fontsize=16, y=0.98)

    black = COLORS["Combined"]
    edges = _sys_edges_for_bin(bin_tag)

    def R(key):
        if rad_bin is None:
            return None
        #endif
        return rad_bin.get(key)

    def M(key):
        if mig_bin is None:
            return None
        #endif
        return mig_bin.get(key)

    def F(key):
        if fermi_bin is None:
            return None
        #endif
        return fermi_bin.get(key)

    def rad_sigma(key, n_target):
        e = R(key)
        if e is None:
            return None
        #endif
        s = _get_sigma_array(e, prefer_abs_delta=False)
        if s is None:
            return None
        #endif
        n = min(n_target, int(s.size))
        if n <= 0:
            return None
        #endif
        return s[:n]

    def mig_sigma_res(key, n_target):
        e = M(key)
        if e is None:
            return None
        #endif
        s = _get_sigma_array(e, prefer_abs_delta=True)
        if s is None:
            return None
        #endif
        n = min(n_target, int(s.size))
        if n <= 0:
            return None
        #endif
        return s[:n]

    def mig_sigma_fermi(key, n_target):
        e = F(key)
        if e is None:
            return None
        #endif
        s = _get_sigma_array(e, prefer_abs_delta=True)
        if s is None:
            return None
        #endif
        n = min(n_target, int(s.size))
        if n <= 0:
            return None
        #endif
        return s[:n]

    # ---------------- BSA ----------------
    if p_comb.get("ALUsin") is not None:
        n = int(p_comb["ALUsin"]["x"].size)
        rs = rad_sigma("ALUsin", n) if with_rad else None
        ms_res = mig_sigma_res("ALUsin", n) if with_mig else None
        ms_fermi = mig_sigma_fermi("ALUsin", n) if with_fermi else None
        _draw_total_sys_band_solid(axLU, edges, sigma_rad=rs, sigma_resmig=ms_res, sigma_fermi=ms_fermi, zorder=1)

        s = p_comb["ALUsin"]
        x, y, ye = _apply_total_shift(
            s,
            rad_entry=R("ALUsin"),
            mig_entry=M("ALUsin"),
            fermi_entry=F("ALUsin"),
            with_rad=with_rad,
            with_mig=with_mig,
            with_fermi=with_fermi,
            label="Combined:ALUsin"
        )
        axLU.errorbar(x, y, yerr=ye, fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)

    # ---------------- TSA ----------------
    if p_comb.get("AULsin") is not None:
        n = int(p_comb["AULsin"]["x"].size)
        rs = rad_sigma("AULsin", n) if with_rad else None
        ms_res = mig_sigma_res("AULsin", n) if with_mig else None
        ms_fermi = mig_sigma_fermi("AULsin", n) if with_fermi else None
        _draw_total_sys_band_solid(axUL, edges, sigma_rad=rs, sigma_resmig=ms_res, sigma_fermi=ms_fermi, zorder=1)

        s1 = p_comb["AULsin"]
        x1, y1, ye1 = _apply_total_shift(
            s1,
            rad_entry=R("AULsin"),
            mig_entry=M("AULsin"),
            fermi_entry=F("AULsin"),
            with_rad=with_rad,
            with_mig=with_mig,
            with_fermi=with_fermi,
            label="Combined:AULsin"
        )
        axUL.errorbar(x1, y1, yerr=ye1, fmt=MARKER, mfc="none", mec=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif

    if p_comb.get("AULsin2") is not None:
        s2 = p_comb["AULsin2"]
        x2, y2, ye2 = _apply_total_shift(
            s2,
            rad_entry=R("AULsin2"),
            mig_entry=M("AULsin2"),
            fermi_entry=F("AULsin2"),
            with_rad=with_rad,
            with_mig=with_mig,
            with_fermi=with_fermi,
            label="Combined:AULsin2"
        )
        axUL.errorbar(x2, y2, yerr=ye2, fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif

    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)

    _leg_h = axUL.legend(
        handles=[
            Line2D([0], [0], marker=MARKER, mfc="none", mec="black", linestyle="", label="n=1"),
            Line2D([0], [0], marker=MARKER, color="black", linestyle="", label="n=2"),
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(_leg_h)

    # ---------------- DSA ----------------
    if p_comb.get("ALLn0") is not None:
        n = int(p_comb["ALLn0"]["x"].size)
        rs = rad_sigma("ALLn0", n) if with_rad else None
        ms_res = mig_sigma_res("ALLn0", n) if with_mig else None
        ms_fermi = mig_sigma_fermi("ALLn0", n) if with_fermi else None
        _draw_total_sys_band_solid(axLL, edges, sigma_rad=rs, sigma_resmig=ms_res, sigma_fermi=ms_fermi, zorder=1)

        s0 = p_comb["ALLn0"]
        x0, y0, ye0 = _apply_total_shift(
            s0,
            rad_entry=R("ALLn0"),
            mig_entry=M("ALLn0"),
            fermi_entry=F("ALLn0"),
            with_rad=with_rad,
            with_mig=with_mig,
            with_fermi=with_fermi,
            label="Combined:ALLn0"
        )
        axLL.errorbar(x0, y0, yerr=ye0, fmt=MARKER, mfc="none", mec=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif

    if p_comb.get("ALLcos") is not None:
        s1 = p_comb["ALLcos"]
        x1, y1, ye1 = _apply_total_shift(
            s1,
            rad_entry=R("ALLcos"),
            mig_entry=M("ALLcos"),
            fermi_entry=F("ALLcos"),
            with_rad=with_rad,
            with_mig=with_mig,
            with_fermi=with_fermi,
            label="Combined:ALLcos"
        )
        axLL.errorbar(x1, y1, yerr=ye1, fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif

    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)

    _leg_h2 = axLL.legend(
        handles=[
            Line2D([0], [0], marker=MARKER, mfc="none", mec="black", linestyle="", label="n=0"),
            Line2D([0], [0], marker=MARKER, color="black", linestyle="", label="n=1"),
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axLL.add_artist(_leg_h2)

    plt.tight_layout(rect=[0, 0, 1, TOP_PAD_PER_BIN])
    os.makedirs(out_dir, exist_ok=True)

    suf = _suffix(with_rad, with_mig, with_fermi)
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_{bin_tag}_CombinedOnly{suf}.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved combined-only figure: {out_path}")

# ---------------------------------------------------------------------
# LaTeX table printing (stdout)
# ---------------------------------------------------------------------
def _assert_same_x(series_list, label_list, tol=1.0e-6):
    """
    Fail-fast: all x arrays must match length and values within tol.
    series_list elements are dicts with "x" arrays.
    """
    ref = series_list[0]["x"]
    for i in range(1, len(series_list)):
        x = series_list[i]["x"]
        if x.size != ref.size:
            raise RuntimeError(f"FATAL: x-size mismatch in LaTeX table: {label_list[0]} has {ref.size}, {label_list[i]} has {x.size}")
        #endif
        if not np.all(np.abs(x - ref) <= tol):
            raise RuntimeError(f"FATAL: x-value mismatch in LaTeX table between {label_list[0]} and {label_list[i]}")
        #endif
    #endfor

def _get_rad_arrays(rad_bin, key, n_expected):
    if rad_bin is None:
        raise RuntimeError("FATAL: rad_bin is None while generating LaTeX tables.")
    #endif
    if key not in rad_bin:
        raise RuntimeError(f"FATAL: missing radiative series in rad summary for key={key}")
    #endif
    delta = np.asarray(rad_bin[key].get("delta", []), dtype=float)
    sigma = np.asarray(rad_bin[key].get("sigma", []), dtype=float)
    if delta.size < n_expected or sigma.size < n_expected:
        raise RuntimeError(f"FATAL: radiative arrays too short for key={key}: delta={delta.size}, sigma={sigma.size}, expected={n_expected}")
    #endif
    return delta[:n_expected], np.abs(sigma[:n_expected])

def _get_mig_delta_arrays(mig_bin, key, n_expected, label):
    """
    Return signed diff arrays (delta) and abs(delta) arrays for a given migration block.
    """
    if mig_bin is None:
        return None, None
    #endif
    if key not in mig_bin or mig_bin[key] is None:
        return None, None
    #endif
    d = np.asarray(mig_bin[key].get("delta", []), dtype=float)
    if d.size < n_expected:
        raise RuntimeError(f"FATAL: migration diff array too short for key={key} in {label}: delta={d.size}, expected={n_expected}")
    #endif
    d = d[:n_expected]
    return d, np.abs(d)

def _fmt_cell(val, stat, syst, ndp=3):
    """
    Wrap the full cell in $...$ to ensure math mode.
    """
    fmt = "{0:." + str(ndp) + "f}"
    core = fmt.format(val) + "^{\\pm " + fmt.format(stat) + "}_{\\pm " + fmt.format(syst) + "}"
    return "$" + core + "$"
# endif

def print_latex_table_for_bin(bin_tag, comb_parsed, rad_all, mig_all, fermi_all):
    """
    Prints one LaTeX table for a given bin_tag using Combined fit results,
    radiative deltas, and total systematics.

    Central values:
      y_corr = y - Delta_rad + diff_mig + diff_fermi

    Total systematic (UPDATED):
      sigma_tot = sqrt( sigma_rad^2 + sigma_resmig^2 + sigma_fermi^2 )
    where:
      sigma_resmig = |diff_mig|   (if --migration present)
      sigma_fermi  = |diff_fermi| (if --fermi_migration present)
    """
    p = build_period_dict(comb_parsed, bin_tag)

    required = ["ALUsin", "AULsin", "AULsin2", "ALLn0", "ALLcos"]
    missing = [k for k in required if (p.get(k) is None)]
    if missing:
        raise RuntimeError(f"FATAL: missing Combined fit series for bin={bin_tag}: {missing}")
    #endif

    series_list = [p[k] for k in required]
    _assert_same_x(series_list, required, tol=1.0e-6)

    x = p["ALUsin"]["x"]
    n = x.size
    if n <= 0:
        raise RuntimeError(f"FATAL: no rows for bin={bin_tag}")
    #endif

    if rad_all is None or mig_all is None:
        raise RuntimeError("FATAL: LaTeX tables require BOTH --rad and --migration.")
    #endif

    rad_bin = rad_all.get(bin_tag)
    mig_bin = mig_all.get(bin_tag)
    fermi_bin = fermi_all.get(bin_tag) if (fermi_all is not None) else None

    if rad_bin is None:
        raise RuntimeError(f"FATAL: missing bin in rad summary: {bin_tag}")
    #endif
    if mig_bin is None:
        raise RuntimeError(f"FATAL: missing bin in migration summary: {bin_tag}")
    #endif

    vals = {}
    stats = {}
    systs = {}

    for key in required:
        y = np.asarray(p[key]["y"], dtype=float)
        ye = np.asarray(p[key]["yerr"], dtype=float)

        d_rad, s_rad = _get_rad_arrays(rad_bin, key, n)

        d_mig,  s_res  = _get_mig_delta_arrays(mig_bin, key, n, label="--migration")
        d_fermi, s_fermi = _get_mig_delta_arrays(fermi_bin, key, n, label="--fermi_migration") if (fermi_bin is not None) else (None, None)

        # Central shifts: y - Delta_rad + diff_mig + diff_fermi
        y_shift = y - d_rad
        if d_mig is not None:
            y_shift = y_shift + d_mig
        #endif
        if d_fermi is not None:
            y_shift = y_shift + d_fermi
        #endif

        # Systematics: quadrature of rad + resmig + fermi
        s_tot_sq = s_rad * s_rad

        if s_res is not None:
            s_tot_sq = s_tot_sq + s_res * s_res
        #endif
        if s_fermi is not None:
            s_tot_sq = s_tot_sq + s_fermi * s_fermi
        #endif

        s_tot = np.sqrt(s_tot_sq)

        vals[key] = y_shift
        stats[key] = ye
        systs[key] = s_tot
    #endfor

    header = (
        "\\begin{table}[h]\n"
        "\\centering\n"
        "\\small\n"
        "\\renewcommand{\\arraystretch}{1.40}\n"
        "\\begin{tabular}{|c|c|c|c|c|c|} \\hline\n"
        "$\\langle -t' \\rangle$ & $F_{LU}^{\\sin\\phi}/F_{UU}$ & $F_{UL}^{\\sin\\phi}/F_{UU}$ & $F_{UL}^{\\sin2\\phi}/F_{UU}$ & $F_{LL}/F_{UU}$ & $F_{LL}^{\\cos\\phi}/F_{UU}$ \\\\ \\hline\n"
    )

    label = f"table:{bin_tag}_fitresults_tprime"
    caption = (
        "\\end{tabular}\n"
        "\\renewcommand{\\arraystretch}{1.0}\n"
        "\\caption{Fitted structure-function ratios per bin. Entries are $\\text{value}^{\\pm\\,\\text{stat}}_{\\pm\\,\\text{syst}}$.\\label{"
        + label +
        "}}\n"
        "\\end{table}\n"
    )

    print("")
    print("% -----------------------------------------------------------------------------")
    print(f"% LaTeX table for bin: {bin_tag}  ({XB_BINS.get(bin_tag, bin_tag)})")
    print("% Central values: y - Delta_rad + diff_mig + diff_fermi (if provided).")
    print("% Total syst: sqrt( sigma_rad^2 + sigma_resmig^2 + sigma_fermi^2 ).")
    print("%   sigma_resmig = |diff_mig|, sigma_fermi = |diff_fermi|.")
    print("% Rows ordered in increasing -t'.")
    print("% -----------------------------------------------------------------------------")
    print(header, end="")

    for i in range(n):
        tmean = float(x[i])

        c1 = _fmt_cell(float(vals["ALUsin"][i]),  float(stats["ALUsin"][i]),  float(systs["ALUsin"][i]),  ndp=3)
        c2 = _fmt_cell(float(vals["AULsin"][i]),  float(stats["AULsin"][i]),  float(systs["AULsin"][i]),  ndp=3)
        c3 = _fmt_cell(float(vals["AULsin2"][i]), float(stats["AULsin2"][i]), float(systs["AULsin2"][i]), ndp=3)
        c4 = _fmt_cell(float(vals["ALLn0"][i]),   float(stats["ALLn0"][i]),   float(systs["ALLn0"][i]),   ndp=3)
        c5 = _fmt_cell(float(vals["ALLcos"][i]),  float(stats["ALLcos"][i]),  float(systs["ALLcos"][i]),  ndp=3)

        row = f"{tmean:.3f} ~&~ {c1} ~&~ {c2} ~&~ {c3} ~&~ {c4} ~&~ {c5} \\\\ \\hline"
        print(row)
    #endfor

    print(caption, end="")

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Plot enpi+ asymmetries vs -t' with optional radiative deltas and (regular + Fermi) migration diffs.")
    ap.add_argument("su22", type=str, help="Su22 fit-results .txt")
    ap.add_argument("fa22", type=str, help="Fa22 fit-results .txt")
    ap.add_argument("sp23", type=str, help="Sp23 fit-results .txt")
    ap.add_argument("combined", type=str, help="Combined fit-results .txt")

    ap.add_argument("--rad", type=str, default=None,
                    help="Radiative Delta summary .txt (block format). Applied as A_Born = A_baseline - Delta.")
    ap.add_argument("--migration", type=str, default=None,
                    help="Regular migration DIFF summary file. Block format preferred; legacy NAME={{x,diff}} supported. Applied as y <- y + diff.")
    ap.add_argument("--fermi_migration", type=str, default=None,
                    help="Fermi-motion migration DIFF summary file. Same formats as --migration. Applied as y <- y + diff.")
    args = ap.parse_args()

    su22 = parse_asym_file(args.su22)
    fa22 = parse_asym_file(args.fa22)
    sp23 = parse_asym_file(args.sp23)
    comb = parse_asym_file(args.combined)

    available_bins = detect_available_bins(su22, fa22, sp23, comb)
    if not available_bins:
        print("[ERROR] No recognizable xB-bin sections found (e.g. enpiLowxBGEchi2FitsALUsinphi).")
        sys.exit(2)
    #endif

    rad_all = None
    if args.rad:
        if not os.path.isfile(args.rad):
            raise RuntimeError(f"FATAL: --rad file not found: {args.rad}")
        #endif
        rad_all = parse_rad_summary(args.rad)

        any_key = False
        for bt, block in rad_all.items():
            for k in _ALLOWED_KEYS:
                if k in block:
                    any_key = True
                    break
                #endif
            #endfor
        #endfor
        if not any_key:
            print("[WARN] rad: parsed file, but did not find any expected keys in any bin.")
            print("[WARN]      This usually means the Series tokens did not map correctly.")
        #endif
    #endif

    mig_all = None
    if args.migration:
        mig_all = parse_migration_any_format(args.migration)
    #endif

    fermi_all = None
    if args.fermi_migration:
        fermi_all = parse_migration_any_format(args.fermi_migration)
    #endif

    out_dir = os.path.join("output", "enpi+")
    os.makedirs(out_dir, exist_ok=True)

    # Produce plots
    for bin_tag in available_bins:
        p_su22 = build_period_dict(su22, bin_tag)
        p_fa22 = build_period_dict(fa22, bin_tag)
        p_sp23 = build_period_dict(sp23, bin_tag)
        p_comb = build_period_dict(comb, bin_tag)

        plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, out_dir)

        # Always produce the unshifted Combined-only figure (as before)
        plot_combined_only_for_bin(
            p_comb, bin_tag, out_dir,
            with_rad=False, with_mig=False, with_fermi=False,
            rad_bin=None, mig_bin=None, fermi_bin=None
        )

        # Produce shifted/systematics figure if we have at least rad+regular migration (same trigger as before),
        # with optional inclusion of fermi_migration if provided.
        if (rad_all is not None) and (mig_all is not None):
            rad_bin = rad_all.get(bin_tag, {})
            mig_bin = mig_all.get(bin_tag, {})
            fermi_bin = fermi_all.get(bin_tag, {}) if (fermi_all is not None) else None

            plot_combined_only_for_bin(
                p_comb, bin_tag, out_dir,
                with_rad=True,
                with_mig=True,
                with_fermi=(fermi_all is not None),
                rad_bin=rad_bin,
                mig_bin=mig_bin,
                fermi_bin=fermi_bin
            )
        #endif
    #endfor

    # Print LaTeX tables (four xB slices only), require at least rad+migration (same trigger as before).
    if (rad_all is not None) and (mig_all is not None):
        table_bins = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE"]
        for bt in table_bins:
            if bt in available_bins:
                print_latex_table_for_bin(bt, comb, rad_all, mig_all, fermi_all)
            else:
                raise RuntimeError(f"FATAL: requested table bin not present in input files: {bt}")
            #endif
        #endfor
    else:
        print("[INFO] Skipping LaTeX tables: requires BOTH --rad and --migration.")
    #endif

if __name__ == "__main__":
    main()
# endif