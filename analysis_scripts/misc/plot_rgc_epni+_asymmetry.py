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

Radiative correction convention (IMPORTANT):
  - We apply: A_Born = A_baseline - Delta
  - Therefore, whenever --rad is used, central values are shifted as:
      y_shift = y - Delta_rad
    and statistical errors remain unchanged.
  - Systematic rectangles represent total systematic per bin:
      sigma_tot = sqrt( sigma_rad^2 + sigma_mig^2 )
    drawn as solid gray rectangles from y=0 to y=sigma_tot spanning each full fixed x bin.

LaTeX tables (stdout) for the 4 xB slices:
  enpiLowxBGE, enpiMidLowxBGE, enpiMidHighxBGE, enpiHighxBGE

Tables are printed only when BOTH --rad and --migration are provided, and use:
  central = (fit value) - (radiative Delta)     [i.e., Born estimate]
  stat    = (fit statistical uncertainty)
  syst    = sqrt( sigma_rad^2 + sigma_mig^2 )

Rows are ordered in increasing -t' (i.e. increasing x axis as plotted).

Usage:
  python plot_rgc_enpi+_asymmetry.py Su22.txt Fa22.txt Sp23.txt Combined.txt --rad output/enpi+/ISR_FSR_delta_summary.txt --migration migration_summary.txt

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

  3) Migration summary file (--migration):
     Supported formats (explicitly detected):
       (A) Block format (preferred, Script-2 style):
           Bin: <bin_tag> ...
           Series: <KEY> | <label>
             -t'    delta_or_sigma   [optional third column sigma]
             0.10   ...
           We interpret:
             - if 2 numeric columns: x and delta; sigma_mig = abs(delta)
             - if 3 numeric columns: x, delta, sigma; sigma_mig = abs(sigma)

       (B) Legacy assignment format (previous functionality):
           NAME = {{mean_tprime, delta}, ...};
           where NAME matches:
             <bin_tag>chi2Fits<suffix>_delta
           and sigma_mig = abs(delta).

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
        # "UUcos":   get_series(parsed, k("AUUcosphi")),
        # "UUcos2":  get_series(parsed, k("AUUcos2phi")),
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

# Map LaTeX-style series labels (as seen in your ISR_FSR_delta_summary.txt) to internal keys.
# We normalize both sides with _norm_series_token before comparing.
_RAD_LABEL_TO_KEY = {
    _norm_series_token(r"F_{LU}^{\sin\phi}/F_{UU}"):      "ALUsin",
    _norm_series_token(r"F_{UL}^{\sin\phi}/F_{UU}"):      "AULsin",
    _norm_series_token(r"F_{UL}^{\sin2\phi}/F_{UU}"):     "AULsin2",
    _norm_series_token(r"F_{LL}/F_{UU}"):                 "ALLn0",
    _norm_series_token(r"F_{LL}^{\cos\phi}/F_{UU}"):      "ALLcos",
    _norm_series_token(r"DilutionD_f"):                   "Df",      # parsed but not used by this script
    _norm_series_token(r"DilutionD_{f}"):                 "Df",      # just in case of brace variant
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

    # If it is KEY | label, take the left side as "token".
    token = rest.split("|", 1)[0].strip() if ("|" in rest) else rest.strip()
    if not token:
        return None
    #endif

    # Direct key match (preferred)
    if token in _ALLOWED_KEYS:
        return token
    #endif

    # Try mapping from LaTeX label to internal key
    norm = _norm_series_token(token)
    mapped = _RAD_LABEL_TO_KEY.get(norm)
    if mapped in _ALLOWED_KEYS:
        return mapped
    #endif

    # Allow parsing dilution (ignored by asymmetry correction logic)
    if mapped == "Df":
        return "Df"
    #endif

    return token  # unknown; caller can decide to ignore

def parse_rad_summary(rad_path):
    """
    Returns:
      rad[bin_tag][series_key] = {"x": array(-t'), "delta": array, "sigma": array}
    where series_key is one of _ALLOWED_KEYS (and optionally "Df").

    This parser is strict about structure, but tolerant about "Series:" tokens:
      - If token matches allowed key, uses it.
      - Else if token matches known LaTeX label, maps to allowed key.
      - Else ignores with INFO.

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

            # skip possible header line(s)
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
# Migration systematics parsers
#   - Block format: Bin/Series blocks, numeric rows
#   - Legacy assignment format: NAME = {{x, delta}, ...};
# ---------------------------------------------------------------------
_pair_re = re.compile(r'\{([^{}]+)\}')

def parse_migration_file_legacy(path):
    """Parse legacy NAME = {{x, delta}, ...}; blocks into dict[name] -> np.array(N,2)."""
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

def get_mig_sigma_series_legacy(dct, key, negate_x=True, sort_x=True):
    """
    Legacy migration *_delta list: (mean_tprime, delta)
    Output: x = -mean_tprime (if negate_x), sigma = abs(delta)
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
    sigma = np.abs(delta[mask])
    if sort_x and x.size > 1:
        order = np.argsort(x)
        x, sigma = x[order], sigma[order]
    #endif
    return {"x": x, "sigma": sigma}

def build_migration_dict_legacy(mig_parsed, bin_prefix):
    if mig_parsed is None:
        return None
    #endif
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}_delta"
    return {
        "ALUsin":  get_mig_sigma_series_legacy(mig_parsed, k("ALUsinphi")),
        "AULsin":  get_mig_sigma_series_legacy(mig_parsed, k("AULsinphi")),
        "AULsin2": get_mig_sigma_series_legacy(mig_parsed, k("AULsin2phi")),
        "ALLn0":   get_mig_sigma_series_legacy(mig_parsed, k("ALL")),
        "ALLcos":  get_mig_sigma_series_legacy(mig_parsed, k("ALLcosphi")),
    }

def parse_migration_summary_block(path):
    """
    Parse block-format migration summary into:
      mig[bin_tag][series_key] = {"x": array(-t'), "sigma": array}

    Numeric rows:
      - If 2 columns: x delta  -> sigma = abs(delta)
      - If 3 columns: x delta sigma -> sigma = abs(sigma)
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

            # skip possible header line(s)
            i += 1
            if i < len(lines):
                hdr = lines[i].strip()
                if (("-t'" in hdr) or ("delta" in hdr) or ("Delta" in hdr) or ("sigma" in hdr)):
                    i += 1
                #endif
            #endif

            xs, sigmas = [], []
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
                    if len(parts) >= 3:
                        sval = float(parts[2])
                        sigma_val = abs(sval)
                    else:
                        sigma_val = abs(dval)
                    #endif
                except Exception:
                    raise RuntimeError(f"FATAL: non-numeric migration row: '{row}'")
                #endif
                xs.append(xval)
                sigmas.append(sigma_val)
                i += 1
            #endfor

            if use_key and xs:
                mig[cur_bin][key] = {
                    "x": np.array(xs, dtype=float),
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

# ---------------------------------------------------------------------
# Titles / bin edges
# ---------------------------------------------------------------------
def make_title(xb_label):
    return r"$ep \rightarrow en\pi^{+}$ -- " + xb_label

def _sys_edges_for_bin(bin_tag):
    return EDGES_SYS_GE if bin_tag == "enpiGE" else EDGES_SYS_XB

# ---------------------------------------------------------------------
# Systematics drawing (solid gray only)
# ---------------------------------------------------------------------
def _compute_total_sigma_per_bin(edges, rad_sigma=None, mig_sigma=None):
    nbins = max(0, len(edges) - 1)
    if nbins <= 0:
        return None
    #endif

    have_rad = (rad_sigma is not None) and (np.asarray(rad_sigma).size > 0)
    have_mig = (mig_sigma is not None) and (np.asarray(mig_sigma).size > 0)

    if not have_rad and not have_mig:
        return None
    #endif

    if have_rad:
        r = np.asarray(rad_sigma, dtype=float)
    else:
        r = None
    #endif

    if have_mig:
        m = np.asarray(mig_sigma, dtype=float)
    else:
        m = None
    #endif

    if have_rad and have_mig:
        n = min(nbins, r.size, m.size)
        tot = np.sqrt(r[:n] * r[:n] + m[:n] * m[:n])
    elif have_rad:
        n = min(nbins, r.size)
        tot = r[:n]
    else:
        n = min(nbins, m.size)
        tot = m[:n]
    #endif

    return tot

def _draw_total_sys_band_solid(ax, edges, rad_sigma=None, mig_sigma=None, zorder=1):
    tot = _compute_total_sigma_per_bin(edges, rad_sigma=rad_sigma, mig_sigma=mig_sigma)
    if tot is None:
        return
    #endif

    nbins = max(0, len(edges) - 1)
    n = min(nbins, tot.size)
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

def _apply_shift_if_available(series, rad_entry):
    """
    Apply radiative shift consistent with:
      A_Born = A_baseline - Delta
    So we do: y_shift = y - Delta.
    """
    if series is None:
        return None
    #endif
    x, y, yerr = series["x"], series["y"], series["yerr"]
    if rad_entry is None:
        return x, y, yerr
    #endif

    xr = np.asarray(rad_entry.get("x", []), dtype=float)
    dr = np.asarray(rad_entry.get("delta", []), dtype=float)

    n = min(x.size, xr.size, dr.size)
    if n == 0:
        return x, y, yerr
    #endif
    if x.size != xr.size:
        print(f"[WARN] rad Delta length ({xr.size}) != data length ({x.size}); pairing by index up to {n}.")
    #endif

    y_shift = y.copy()
    y_shift[:n] = y_shift[:n] - dr[:n]
    return x, y_shift, yerr

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

def _suffix(with_rad, with_mig):
    if with_rad and with_mig:
        return "_withRadMig"
    #endif
    if with_rad:
        return "_withRad"
    #endif
    if with_mig:
        return "_withMig"
    #endif
    return ""

# ---------------------------------------------------------------------
# Plotting (1x3 canvases: LU, UL, LL only)
# ---------------------------------------------------------------------
def _plot_panel_sets_1x3(axLU, axUL, axLL, pdata_by_label,
                        rad_for_bin=None, mig_for_bin=None,
                        with_rad=False, with_mig=False, bin_tag_for_width="enpiGE"):
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

    def rad_sigma(key):
        e = R(key)
        if e is None:
            return None
        #endif
        return np.abs(np.asarray(e.get("sigma", []), dtype=float))

    def mig_sigma(key):
        e = M(key)
        if e is None:
            return None
        #endif
        return np.abs(np.asarray(e.get("sigma", []), dtype=float))

    # ---------------- BSA (ALU sin phi) ----------------
    rs = rad_sigma("ALUsin") if with_rad else None
    ms = mig_sigma("ALUsin") if with_mig else None
    _draw_total_sys_band_solid(axLU, edges, rad_sigma=rs, mig_sigma=ms, zorder=1)

    for lab, pdata in pdata_by_label.items():
        s = pdata.get("ALUsin")
        if s is None:
            continue
        #endif
        if with_rad:
            x, y, ye = _apply_shift_if_available(s, R("ALUsin"))
        else:
            x, y, ye = s["x"], s["y"], s["yerr"]
        #endif
        axLU.errorbar(x, y, yerr=ye, fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                      capsize=CAPSIZE, label=lab, markersize=MS, linestyle="None")
    #endfor
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)
    _legend_run_period(axLU, list(pdata_by_label.keys()), where="lower right")

    # ---------------- TSA (AUL sin phi open, sin2 phi filled) ----------------
    rs = rad_sigma("AULsin") if with_rad else None
    ms = mig_sigma("AULsin") if with_mig else None
    _draw_total_sys_band_solid(axUL, edges, rad_sigma=rs, mig_sigma=ms, zorder=1)

    for lab, pdata in pdata_by_label.items():
        s1 = pdata.get("AULsin")
        if s1 is not None:
            if with_rad:
                x1, y1, ye1 = _apply_shift_if_available(s1, R("AULsin"))
            else:
                x1, y1, ye1 = s1["x"], s1["y"], s1["yerr"]
            #endif
            axUL.errorbar(x1, y1, yerr=ye1, fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif

        s2 = pdata.get("AULsin2")
        if s2 is not None:
            if with_rad and (R("AULsin2") is not None):
                x2, y2, ye2 = _apply_shift_if_available(s2, R("AULsin2"))
            else:
                x2, y2, ye2 = s2["x"], s2["y"], s2["yerr"]
            #endif
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
    rs = rad_sigma("ALLn0") if with_rad else None
    ms = mig_sigma("ALLn0") if with_mig else None
    _draw_total_sys_band_solid(axLL, edges, rad_sigma=rs, mig_sigma=ms, zorder=1)

    for lab, pdata in pdata_by_label.items():
        s0 = pdata.get("ALLn0")
        if s0 is not None:
            if with_rad:
                x0, y0, ye0 = _apply_shift_if_available(s0, R("ALLn0"))
            else:
                x0, y0, ye0 = s0["x"], s0["y"], s0["yerr"]
            #endif
            axLL.errorbar(x0, y0, yerr=ye0, fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif

        s1 = pdata.get("ALLcos")
        if s1 is not None:
            if with_rad and (R("ALLcos") is not None):
                x1, y1, ye1 = _apply_shift_if_available(s1, R("ALLcos"))
            else:
                x1, y1, ye1 = s1["x"], s1["y"], s1["yerr"]
            #endif
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
                        rad_for_bin=None, mig_for_bin=None,
                        with_rad=False, with_mig=False, bin_tag_for_width=bin_tag)

    plt.tight_layout(rect=[0, 0, 1, TOP_PAD_PER_BIN])
    os.makedirs(out_dir, exist_ok=True)

    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_{bin_tag}_AllPeriods.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved all-periods figure: {out_path}")

def plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                               with_rad=False, with_mig=False, rad_bin=None, mig_bin=None):
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

    def rad_sigma(key):
        e = R(key)
        if e is None:
            return None
        #endif
        return np.abs(np.asarray(e.get("sigma", []), dtype=float))

    def mig_sigma(key):
        e = M(key)
        if e is None:
            return None
        #endif
        return np.abs(np.asarray(e.get("sigma", []), dtype=float))

    # ---------------- BSA ----------------
    rs = rad_sigma("ALUsin") if with_rad else None
    ms = mig_sigma("ALUsin") if with_mig else None
    _draw_total_sys_band_solid(axLU, edges, rad_sigma=rs, mig_sigma=ms, zorder=1)

    if p_comb.get("ALUsin") is not None:
        s = p_comb["ALUsin"]
        if with_rad:
            x, y, ye = _apply_shift_if_available(s, R("ALUsin"))
        else:
            x, y, ye = s["x"], s["y"], s["yerr"]
        #endif
        axLU.errorbar(x, y, yerr=ye, fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)

    # ---------------- TSA ----------------
    rs = rad_sigma("AULsin") if with_rad else None
    ms = mig_sigma("AULsin") if with_mig else None
    _draw_total_sys_band_solid(axUL, edges, rad_sigma=rs, mig_sigma=ms, zorder=1)

    if p_comb.get("AULsin") is not None:
        s1 = p_comb["AULsin"]
        if with_rad:
            x1, y1, ye1 = _apply_shift_if_available(s1, R("AULsin"))
        else:
            x1, y1, ye1 = s1["x"], s1["y"], s1["yerr"]
        #endif
        axUL.errorbar(x1, y1, yerr=ye1, fmt=MARKER, mfc="none", mec=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif

    if p_comb.get("AULsin2") is not None:
        s2 = p_comb["AULsin2"]
        if with_rad and (R("AULsin2") is not None):
            x2, y2, ye2 = _apply_shift_if_available(s2, R("AULsin2"))
        else:
            x2, y2, ye2 = s2["x"], s2["y"], s2["yerr"]
        #endif
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
    rs = rad_sigma("ALLn0") if with_rad else None
    ms = mig_sigma("ALLn0") if with_mig else None
    _draw_total_sys_band_solid(axLL, edges, rad_sigma=rs, mig_sigma=ms, zorder=1)

    if p_comb.get("ALLn0") is not None:
        s0 = p_comb["ALLn0"]
        if with_rad:
            x0, y0, ye0 = _apply_shift_if_available(s0, R("ALLn0"))
        else:
            x0, y0, ye0 = s0["x"], s0["y"], s0["yerr"]
        #endif
        axLL.errorbar(x0, y0, yerr=ye0, fmt=MARKER, mfc="none", mec=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif

    if p_comb.get("ALLcos") is not None:
        s1 = p_comb["ALLcos"]
        if with_rad and (R("ALLcos") is not None):
            x1, y1, ye1 = _apply_shift_if_available(s1, R("ALLcos"))
        else:
            x1, y1, ye1 = s1["x"], s1["y"], s1["yerr"]
        #endif
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

    suf = _suffix(with_rad, with_mig)
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
    return delta[:n_expected], sigma[:n_expected]

def _get_mig_arrays(mig_bin, key, n_expected):
    if mig_bin is None:
        raise RuntimeError("FATAL: mig_bin is None while generating LaTeX tables.")
    #endif
    if key not in mig_bin or mig_bin[key] is None:
        raise RuntimeError(f"FATAL: missing migration series for key={key}")
    #endif
    sigma = np.asarray(mig_bin[key].get("sigma", []), dtype=float)
    if sigma.size < n_expected:
        raise RuntimeError(f"FATAL: migration sigma array too short for key={key}: sigma={sigma.size}, expected={n_expected}")
    #endif
    return sigma[:n_expected]

def _fmt_cell(val, stat, syst, ndp=3):
    """
    Wrap the full cell in $...$ to ensure math mode, since we use ^ and _.
    """
    fmt = "{0:." + str(ndp) + "f}"
    core = fmt.format(val) + "^{\\pm " + fmt.format(stat) + "}_{\\pm " + fmt.format(syst) + "}"
    return "$" + core + "$"
#endif

def print_latex_table_for_bin(bin_tag, comb_parsed, rad_all, mig_all):
    """
    Prints one LaTeX table for a given bin_tag using Combined fit results,
    radiative deltas, and quadrature systematics.

    Central values correspond to A_Born = A_baseline - Delta_rad.
    """
    p = build_period_dict(comb_parsed, bin_tag)

    required = ["ALUsin", "AULsin", "AULsin2", "ALLn0", "ALLcos"]
    missing = [k for k in required if (p.get(k) is None)]
    if missing:
        raise RuntimeError(f"FATAL: missing Combined fit series for bin={bin_tag}: {missing}")
    #endif

    # Ensure same x binning across all five series
    series_list = [p[k] for k in required]
    _assert_same_x(series_list, required, tol=1.0e-6)

    x = p["ALUsin"]["x"]
    n = x.size
    if n <= 0:
        raise RuntimeError(f"FATAL: no rows for bin={bin_tag}")
    #endif

    if rad_all is None or mig_all is None:
        raise RuntimeError("FATAL: LaTeX tables require BOTH --rad and --migration so systematics can be sqrt(rad^2 + mig^2).")
    #endif

    rad_bin = rad_all.get(bin_tag)
    mig_bin = mig_all.get(bin_tag)
    if rad_bin is None:
        raise RuntimeError(f"FATAL: missing bin in rad summary: {bin_tag}")
    #endif
    if mig_bin is None:
        raise RuntimeError(f"FATAL: missing bin in migration summary: {bin_tag}")
    #endif

    # Prepare arrays: central shift = y - d_rad, sys from quadrature of sigmas
    vals = {}
    stats = {}
    systs = {}

    for key in required:
        y = np.asarray(p[key]["y"], dtype=float)
        ye = np.asarray(p[key]["yerr"], dtype=float)

        d_rad, s_rad = _get_rad_arrays(rad_bin, key, n)
        s_mig = _get_mig_arrays(mig_bin, key, n)

        y_shift = y - d_rad
        s_tot = np.sqrt(s_rad * s_rad + s_mig * s_mig)

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
    print("% Central values correspond to A_Born = A_baseline - Delta_rad.")
    print("% Syst is sqrt(rad_sigma^2 + mig_sigma^2). Rows ordered in increasing -t'.")
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
    ap = argparse.ArgumentParser(description="Plot enpi+ asymmetries vs -t' with optional radiative deltas and bin-migration systematics.")
    ap.add_argument("su22", type=str, help="Su22 fit-results .txt")
    ap.add_argument("fa22", type=str, help="Fa22 fit-results .txt")
    ap.add_argument("sp23", type=str, help="Sp23 fit-results .txt")
    ap.add_argument("combined", type=str, help="Combined fit-results .txt")
    ap.add_argument("--rad", type=str, default=None, help="Radiative Delta summary .txt (block format). Applied as A_Born = A_baseline - Delta.")
    ap.add_argument("--migration", type=str, default=None, help="Migration summary file. Block format preferred; legacy NAME={{x,delta}} supported.")
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

        # Optional sanity: warn if none of the allowed keys appear for any bin.
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
            print("[WARN] rad: parsed file, but did not find any of the expected keys ALUsin/AULsin/AULsin2/ALLn0/ALLcos in any bin.")
            print("[WARN]      This usually means the Series tokens did not map correctly.")
        #endif
    #endif

    mig_all = None
    if args.migration:
        if not os.path.isfile(args.migration):
            raise RuntimeError(f"FATAL: --migration file not found: {args.migration}")
        #endif

        if _file_has_bin_blocks(args.migration):
            print("[INFO] migration: detected block-format migration summary (Bin:/Series:).")
            mig_all = parse_migration_summary_block(args.migration)
        else:
            print("[INFO] migration: detected legacy NAME={{x,delta}} format; using legacy parser.")
            mig_parsed = parse_migration_file_legacy(args.migration)
            mig_all = {}
            for bin_tag in BIN_ORDER:
                mig_all[bin_tag] = build_migration_dict_legacy(mig_parsed, bin_tag) or {}
            #endfor
        #endif
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

        plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                   with_rad=False, with_mig=False, rad_bin=None, mig_bin=None)

        if (rad_all is not None) and (mig_all is not None):
            rad_bin = rad_all.get(bin_tag, {})
            mig_bin = mig_all.get(bin_tag, {})
            plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                       with_rad=True, with_mig=True, rad_bin=rad_bin, mig_bin=mig_bin)
        #endif
    #endfor

    # Print LaTeX tables (four xB slices only), require both rad+migration
    if (rad_all is not None) and (mig_all is not None):
        table_bins = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE"]
        for bt in table_bins:
            if bt in available_bins:
                print_latex_table_for_bin(bt, comb, rad_all, mig_all)
            else:
                raise RuntimeError(f"FATAL: requested table bin not present in input files: {bt}")
            #endif
        #endfor
    else:
        print("[INFO] Skipping LaTeX tables: requires BOTH --rad and --migration.")
    #endif

if __name__ == "__main__":
    main()
#endif