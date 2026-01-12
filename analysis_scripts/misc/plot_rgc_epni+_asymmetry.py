#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ep -> en pi+ asymmetries versus -t' in several xB bins, for three run periods
and for the combined (inverse-variance weighted) file.

Usage:
  python plot_rgc_enpi_asymmetry.py Su22.txt Fa22.txt Sp23.txt Combined.txt [--rad isr_fsr_delta_summary.txt] [--migration migration_summary.txt]

Outputs (for each available xB bin):
  output/enpi+/rgc_enpi_<BinTag>_AllPeriods.pdf
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly.pdf

If --rad is provided (radiative shifts and radiative systematics):
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly_withRad.pdf

If --migration is provided (bin-migration systematics):
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly_withMig.pdf

If BOTH --rad and --migration are provided:
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly_withRadMig.pdf

Overlay (Combined-only, 1x3):
  output/enpi+/rgc_enpi_xBOverlay_CombinedOnly.pdf
  (and if --rad)        output/enpi+/rgc_enpi_xBOverlay_CombinedOnly_withRad.pdf
  (and if --migration)  output/enpi+/rgc_enpi_xBOverlay_CombinedOnly_withMig.pdf
  (and if both)         output/enpi+/rgc_enpi_xBOverlay_CombinedOnly_withRadMig.pdf

Notes:
- Fit-result text files must contain blocks like:
    enpiLowxBGEchi2FitsALUsinphi = {{mean_tprime, value, error}, ...};
  x is plotted as -t'.

Radiative option (--rad):
- Central values are shifted by the signed Delta from the summary;
  stat errors remain unchanged.
- Radiative sigma (sigmaDelta) is used as a systematic per bin.

Migration option (--migration):
- For each spin structure, uses the *_delta list:
    <name>_delta = {{mean_tprime, delta}, ...};
  and treats abs(delta) as a bin-migration systematic uncertainty per bin.
- No central-value shift is applied for migration; it is systematic only.

Systematics drawing style:
- Solid gray rectangles only.
- If both rad and migration exist, draw ONLY the quadrature total:
    sigma_total = sqrt( sigma_rad^2 + sigma_mig^2 )
- No dedicated systematic legend is drawn.

Systematics inclusion rules:
- TSA systematics use only F_{UL}^{sin phi}; sin2 phi systematics are ignored.
- DSA systematics use only F_{LL}; F_{LL}^{cos phi} systematics are ignored.
- On the 1x3 overlay, systematic bands (if present) are taken only from the
  lowest xB bin (enpiLowxBGE), and use the xB (0.20-wide) edges.

Bin edges for systematic rectangles (x axis is +(-t')):
- For enpiGE (integrated), 0.10-wide bins:
    0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25
- For xB sub-bins, 0.20-wide bins:
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

# Colors for the four xB slices on the 1x3 overlay
XB_COLORS = {
    "enpiLowxBGE":     "tab:blue",
    "enpiMidLowxBGE":  "tab:orange",
    "enpiMidHighxBGE": "tab:green",
    "enpiHighxBGE":    "tab:red",
}

# Marker shapes per xB slice (overlay)
XB_MARKERS = {
    "enpiLowxBGE":     "o",
    "enpiMidLowxBGE":  "s",
    "enpiMidHighxBGE": "^",
    "enpiHighxBGE":    "D",
}

MARKER = "o"
CAPSIZE = 3
MS = 5.0  # marker size

# Axes ranges/labels (x is -t')
XLIM_T = (0.0, 1.30)
X_LABEL = r"$-t'\ (\mathrm{GeV}^{2})$"

# Updated y ranges per your request:
# - Single-spin: [-0.4, 0.4]  (BSA + TSA)
# - Double-spin: [-1.0, 1.0]  (DSA)
YLIM_LU = (-0.4, 0.4)   # BSA
YLIM_UL = (-0.4, 0.4)   # TSA
YLIM_LL = (-1.0, 1.0)   # DSA
YLIM_UU = (-0.4, 0.4)   # UU (unchanged here; adjust if you want)

# Human-readable labels for xB bins (for legends and titles)
XB_BINS = {
    "enpiLowxBGE":     r"$0.10 < x_{B} < 0.25$",
    "enpiMidLowxBGE":  r"$0.25 < x_{B} < 0.35$",
    "enpiMidHighxBGE": r"$0.35 < x_{B} < 0.45$",
    "enpiHighxBGE":    r"$0.45 < x_{B} < 0.60$",
    "enpiGE":          r"$0.10 < x_{B} < 0.60$"
}

# Order to render canvases
BIN_ORDER = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE", "enpiGE"]

# Top padding
TOP_PAD_PER_BIN = 0.95
TOP_PAD_OVERLAY = 0.94

OUT_PREFIX = "enpi"

# Fixed bin-edge sets for systematic rectangles (on +(-t') axis)
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
    """Collect the series for one period given a bin_prefix like 'enpiLowxBGE' or 'enpiGE'."""
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}"
    return {
        "ALUsin":  get_series(parsed, k("ALUsinphi")),
        "AULsin":  get_series(parsed, k("AULsinphi")),
        "AULsin2": get_series(parsed, k("AULsin2phi")),
        "ALLn0":   get_series(parsed, k("ALL")),
        "ALLcos":  get_series(parsed, k("ALLcosphi")),
        "UUcos":   get_series(parsed, k("AUUcosphi")),
        "UUcos2":  get_series(parsed, k("AUUcos2phi")),
    }

def detect_available_bins(*dicts):
    """Detect which bin_prefix sections exist by probing for ALUsinphi keys."""
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
# Radiative-correction summary parser (from signed Delta summary txt)
# ---------------------------------------------------------------------
_SERIES_KEYWORDS = {
    "ALUsin":   r"F_{LU}^{\sin\phi}/F_{UU}",
    "AULsin":   r"F_{UL}^{\sin\phi}/F_{UU}",
    "AULsin2":  r"F_{UL}^{\sin2\phi}/F_{UU}",
    "ALLn0":    r"F_{LL}/F_{UU}",
    "ALLcos":   r"F_{LL}^{\cos\phi}/F_{UU}",
}

def parse_rad_summary(rad_path):
    """
    Parse signed-Delta summary text into:
        rad[bin_tag][series_key] = dict(x=array, delta=array, sigma=array)

    Expected pattern (whitespace-tolerant):
      Bin: enpiGE  x_B range: ...
      Series: F_{LU}^{\\sin\\phi}/F_{UU}
        -t'        Delta        sigmaDelta
        0.10     -0.0123       0.0045
        ...
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
            #endif
            i += 1
            continue
        #endif

        if line.startswith("Series:"):
            ser_line = line[len("Series:"):].strip()
            cur_series_key = None
            for key, tag in _SERIES_KEYWORDS.items():
                if tag in ser_line:
                    cur_series_key = key
                    break
                #endif
            #endfor

            i += 1
            if i < len(lines):
                hdr = lines[i]
                if ("-t'" in hdr) or ("Delta" in hdr) or ("sigma" in hdr):
                    i += 1
                #endif
            #endif

            xs, deltas, sigmas = [], [], []
            while i < len(lines):
                row = lines[i].strip()
                if (not row) or row.startswith("Series:") or row.startswith("Bin:") or set(row) <= set("=-"):
                    break
                #endif
                parts = row.split()
                try:
                    xval = float(parts[0])
                    dval = float(parts[1])
                    sval = float(parts[2])
                except Exception:
                    break
                #endif
                xs.append(xval)
                deltas.append(dval)
                sigmas.append(sval)
                i += 1
            #endfor

            if cur_bin is not None and cur_series_key is not None and xs:
                rad[cur_bin][cur_series_key] = {
                    "x": np.array(xs, dtype=float),         # already -t'
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
# Migration systematics parser (from *_delta lists)
# ---------------------------------------------------------------------
_pair_re = re.compile(r'\{([^{}]+)\}')

def parse_migration_file(path):
    """Parse NAME = {{x, val}, ...}; blocks into dict[name] -> np.array(N,2)."""
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

def get_mig_sigma_series(dct, key, negate_x=True, sort_x=True):
    """
    Read a migration *_delta list as a systematic sigma series:
      input array columns: (mean_tprime, delta)
      output: x = -mean_tprime (if negate_x), sigma = abs(delta)
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

def build_migration_dict(mig_parsed, bin_prefix):
    """Collect migration systematics for one bin_prefix from *_delta lists."""
    if mig_parsed is None:
        return None
    #endif
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}_delta"
    return {
        "ALUsin":  get_mig_sigma_series(mig_parsed, k("ALUsinphi")),
        "AULsin":  get_mig_sigma_series(mig_parsed, k("AULsinphi")),
        "AULsin2": get_mig_sigma_series(mig_parsed, k("AULsin2phi")),
        "ALLn0":   get_mig_sigma_series(mig_parsed, k("ALL")),
        "ALLcos":  get_mig_sigma_series(mig_parsed, k("ALLcosphi")),
    }

# ---------------------------------------------------------------------
# Titles (2x2 canvases only)
# ---------------------------------------------------------------------
def make_title(xb_label):
    return r"$ep \rightarrow en\pi^{+}$ -- " + xb_label

# ---------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------
def _legend_run_period(ax, labels):
    handles = [Line2D([0], [0], marker=MARKER, color=COLORS[L], linestyle="", label=L) for L in labels]
    leg = ax.legend(handles=handles, title="Run Period", frameon=True, edgecolor="black",
                    loc="lower right", fontsize=11, title_fontsize=12)
    leg.get_frame().set_alpha(0.9)

def _legend_harmonic(ax, labels=("n=1", "n=2"), where="lower left"):
    open_marker = Line2D([0], [0], marker=MARKER, mfc="none", mec="black", linestyle="", label=labels[0])
    filled_marker = Line2D([0], [0], marker=MARKER, color="black", linestyle="", label=labels[1])
    leg = ax.legend(handles=[open_marker, filled_marker], title="Harmonic",
                    frameon=True, edgecolor="black", loc=where,
                    bbox_to_anchor=(0.02, 0.02), fontsize=11, title_fontsize=12)
    ax.add_artist(leg)

def _sys_edges_for_bin(bin_tag):
    return EDGES_SYS_GE if bin_tag == "enpiGE" else EDGES_SYS_XB

def _compute_total_sigma_per_bin(edges, rad_sigma=None, mig_sigma=None):
    """
    Return total sigma array aligned to bin count:
      - If both provided: quadrature sqrt(rad^2 + mig^2)
      - Else: whichever exists
      - If neither: None
    """
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

def _draw_total_sys_band_solid(ax, edges, rad_sigma=None, mig_sigma=None, *, zorder=1):
    """
    Draw ONLY a solid gray systematic band:
      - If both rad and migration exist: quadrature total
      - Else: whichever exists
    """
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
    Align by index order (no nearest-neighbor matching).
    Returns (x, y_shift, yerr) if rad_entry exists; else (x, y, yerr).
    """
    if series is None:
        return None
    #endif
    x, y, yerr = series["x"], series["y"], series["yerr"]
    if rad_entry is None:
        return x, y, yerr
    #endif

    xr = np.asarray(rad_entry["x"])
    dr = np.asarray(rad_entry["delta"])

    n = min(x.size, xr.size, dr.size)
    if n == 0:
        return x, y, yerr
    #endif
    if x.size != xr.size:
        print(f"[WARN] rad Delta length ({xr.size}) != data length ({x.size}); pairing by index up to {n}.")
    #endif
    y_shift = y.copy()
    y_shift[:n] = y_shift[:n] + dr[:n]
    return x, y_shift, yerr

def _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label,
                     rad_for_bin=None, mig_for_bin=None,
                     with_rad=False, with_mig=False, *, bin_tag_for_width="enpiGE"):
    """
    Draw the four panels for a set of labeled period dicts.

    Systematic band style:
      - Solid gray only
      - If both rad and migration exist: quadrature total only
      - No systematic legend

    TSA systematics: only AULsin is used (ignore sin2).
    DSA systematics: only ALLn0 is used (ignore ALLcos).
    """
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
    _legend_run_period(axLU, list(pdata_by_label.keys()))

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
    _legend_run_period(axUL, list(pdata_by_label.keys()))

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
    _legend_run_period(axLL, list(pdata_by_label.keys()))

    # ---------------- UU (no sys treatment) ----------------
    for lab, pdata in pdata_by_label.items():
        s1 = pdata.get("UUcos")
        if s1 is not None:
            axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"], fmt="o", mfc="none",
                          mec=COLORS[lab], ecolor=COLORS[lab], capsize=CAPSIZE,
                          markersize=MS, linestyle="None")
        #endif
        s2 = pdata.get("UUcos2")
        if s2 is not None:
            axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"], fmt="o", color=COLORS[lab],
                          ecolor=COLORS[lab], capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif
    #endfor
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axUU, labels=("n=1", "n=2"))
    _legend_run_period(axUU, list(pdata_by_label.keys()))

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
# Per-bin canvases (2x2)
# ---------------------------------------------------------------------
def plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, out_dir):
    plt.figure(figsize=(12, 9))
    axLU = plt.subplot(2, 2, 1)
    axUL = plt.subplot(2, 2, 2)
    axLL = plt.subplot(2, 2, 3)
    axUU = plt.subplot(2, 2, 4)

    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(xb_label), fontsize=16, y=0.97)

    pdata_by_label = {"Su22": p_su22, "Fa22": p_fa22, "Sp23": p_sp23}
    _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label,
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
    plt.figure(figsize=(12, 9))
    axLU = plt.subplot(2, 2, 1)
    axUL = plt.subplot(2, 2, 2)
    axLL = plt.subplot(2, 2, 3)
    axUU = plt.subplot(2, 2, 4)

    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(xb_label), fontsize=16, y=0.97)

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

    # ---------------- TSA (bars from AULsin only) ----------------
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

    # ---------------- DSA (bars from ALLn0 only) ----------------
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

    # ---------------- UU (no sys) ----------------
    if p_comb.get("UUcos") is not None:
        s1 = p_comb["UUcos"]
        axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"], fmt=MARKER, mfc="none", mec=black,
                      ecolor=black, capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif
    if p_comb.get("UUcos2") is not None:
        s2 = p_comb["UUcos2"]
        axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"], fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    #endif
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout(rect=[0, 0, 1, TOP_PAD_PER_BIN])
    os.makedirs(out_dir, exist_ok=True)

    suf = _suffix(with_rad, with_mig)
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_{bin_tag}_CombinedOnly{suf}.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved combined-only figure: {out_path}")

# ---------------------------------------------------------------------
# xB overlay 1x3 canvas (Combined only)
# ---------------------------------------------------------------------
def plot_combined_xb_overlay_1x3(comb_parsed, out_dir, bins_to_use,
                                 with_rad=False, with_mig=False, rad_all=None, mig_all=None):
    """
    1x3 canvas (Combined only) overlaying the four xB slices:
      Left  : F_LU^{sin phi}/F_UU  (ALUsin)
      Middle: F_UL^{sin phi}/F_UU  (AULsin)
      Right : F_UL^{sin2 phi}/F_UU (AULsin2)

    Systematic band (solid gray) is taken ONLY from the lowest xB bin and:
      - For ALUsin: uses total quadrature if both provided
      - For AULsin: uses total quadrature if both provided
      - For AULsin2: no systematics drawn
    """
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharex=True)
    axL, axM, axR = axes

    low_bin_tag = "enpiLowxBGE"
    edges_overlay = EDGES_SYS_XB

    def low_rad_sigma(key):
        if (not with_rad) or (rad_all is None):
            return None
        #endif
        low_rad = rad_all.get(low_bin_tag, {})
        e = low_rad.get(key)
        if e is None:
            return None
        #endif
        return np.abs(np.asarray(e.get("sigma", []), dtype=float))

    def low_mig_sigma(key):
        if (not with_mig) or (mig_all is None):
            return None
        #endif
        low_mig = mig_all.get(low_bin_tag, {})
        e = low_mig.get(key)
        if e is None:
            return None
        #endif
        return np.abs(np.asarray(e.get("sigma", []), dtype=float))

    def _draw_component(ax, suffix_key, ylabel, ylim, draw_sys=True):
        if draw_sys:
            rs = low_rad_sigma(suffix_key)
            ms = low_mig_sigma(suffix_key)
            _draw_total_sys_band_solid(ax, edges_overlay, rad_sigma=rs, mig_sigma=ms, zorder=1)
        #endif

        for bin_tag in bins_to_use:
            pdata = build_period_dict(comb_parsed, bin_tag)
            s = pdata.get(suffix_key)
            if s is None:
                continue
            #endif

            clr = XB_COLORS.get(bin_tag, "black")
            mkr = XB_MARKERS.get(bin_tag, "o")

            if with_rad and (rad_all is not None):
                rad_bin = rad_all.get(bin_tag, {})
                rad_entry = rad_bin.get(suffix_key)
                if rad_entry is not None:
                    x, y, ye = _apply_shift_if_available(s, rad_entry)
                else:
                    x, y, ye = s["x"], s["y"], s["yerr"]
                #endif
            else:
                x, y, ye = s["x"], s["y"], s["yerr"]
            #endif

            ax.errorbar(x, y, yerr=ye, fmt=mkr, color=clr, ecolor=clr, capsize=CAPSIZE,
                        label=XB_BINS.get(bin_tag, bin_tag), markersize=MS, linestyle="None")
        #endfor

        ax.set(xlim=XLIM_T, ylim=ylim, xlabel=X_LABEL, ylabel=ylabel)
        ax.axhline(0, color="black", linestyle="--", linewidth=1.2)
        ax.grid(True, linestyle="--", alpha=0.6)
        leg = ax.legend(frameon=True, edgecolor="black", fontsize=11, loc="best")
        leg.get_frame().set_alpha(0.9)

    _draw_component(axL, "ALUsin",  r"$F_{LU}^{\sin\phi}/F_{UU}$", YLIM_LU, draw_sys=True)
    _draw_component(axM, "AULsin",  r"$F_{UL}^{\sin\phi}/F_{UU}$", YLIM_UL, draw_sys=True)
    _draw_component(axR, "AULsin2", r"$F_{UL}^{\sin2\phi}/F_{UU}$", YLIM_UL, draw_sys=False)

    fig.tight_layout(rect=[0, 0, 1, TOP_PAD_OVERLAY])
    os.makedirs(out_dir, exist_ok=True)
    suf = _suffix(with_rad, with_mig)
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_xBOverlay_CombinedOnly{suf}.pdf")
    plt.savefig(out_path)
    plt.close(fig)
    print(f"Saved xB-overlay (1x3) figure: {out_path}")

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Plot enpi+ asymmetries vs -t' with optional radiative shifts and bin-migration systematics.")
    ap.add_argument("su22", type=str, help="Su22 fit-results .txt")
    ap.add_argument("fa22", type=str, help="Fa22 fit-results .txt")
    ap.add_argument("sp23", type=str, help="Sp23 fit-results .txt")
    ap.add_argument("combined", type=str, help="Combined fit-results .txt")
    ap.add_argument("--rad", type=str, default=None, help="Signed Delta summary .txt to apply (+Delta) with sigmaDelta bars spanning full bins")
    ap.add_argument("--migration", type=str, default=None, help="Migration summary .txt; uses *_delta lists as abs(delta) systematic per bin")
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
            print(f"[WARN] --rad file not found: {args.rad}. Proceeding without radiative shifts.")
        else:
            rad_all = parse_rad_summary(args.rad)
        #endif
    #endif

    mig_all = None
    if args.migration:
        if not os.path.isfile(args.migration):
            print(f"[WARN] --migration file not found: {args.migration}. Proceeding without migration systematics.")
        else:
            mig_parsed = parse_migration_file(args.migration)
            mig_all = {}
            for bin_tag in BIN_ORDER:
                mig_all[bin_tag] = build_migration_dict(mig_parsed, bin_tag) or {}
            #endfor
        #endif
    #endif

    out_dir = os.path.join("output", "enpi+")
    os.makedirs(out_dir, exist_ok=True)

    for bin_tag in available_bins:
        p_su22 = build_period_dict(su22, bin_tag)
        p_fa22 = build_period_dict(fa22, bin_tag)
        p_sp23 = build_period_dict(sp23, bin_tag)
        p_comb = build_period_dict(comb, bin_tag)

        plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, out_dir)

        plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                   with_rad=False, with_mig=False, rad_bin=None, mig_bin=None)

        if (rad_all is not None) or (mig_all is not None):
            rad_bin = rad_all.get(bin_tag, {}) if (rad_all is not None) else None
            mig_bin = mig_all.get(bin_tag, {}) if (mig_all is not None) else None

            if (mig_all is not None) and (rad_all is None):
                plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                           with_rad=False, with_mig=True, rad_bin=None, mig_bin=mig_bin)
            #endif

            if (rad_all is not None) and (mig_all is None):
                plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                           with_rad=True, with_mig=False, rad_bin=rad_bin, mig_bin=None)
            #endif

            if (rad_all is not None) and (mig_all is not None):
                plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                           with_rad=True, with_mig=True, rad_bin=rad_bin, mig_bin=mig_bin)
            #endif
        #endif
    #endfor

    xb_overlay_candidates = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE"]
    bins_to_use = [b for b in xb_overlay_candidates if b in available_bins]
    if bins_to_use:
        plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use,
                                     with_rad=False, with_mig=False, rad_all=None, mig_all=None)

        if (mig_all is not None) and (rad_all is None):
            plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use,
                                         with_rad=False, with_mig=True, rad_all=None, mig_all=mig_all)
        #endif

        if (rad_all is not None) and (mig_all is None):
            plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use,
                                         with_rad=True, with_mig=False, rad_all=rad_all, mig_all=None)
        #endif

        if (rad_all is not None) and (mig_all is not None):
            plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use,
                                         with_rad=True, with_mig=True, rad_all=rad_all, mig_all=mig_all)
        #endif
    else:
        print("[WARN] No xB-slice data available for overlay canvas; skipping.")
    #endif

if __name__ == "__main__":
    main()
#endif