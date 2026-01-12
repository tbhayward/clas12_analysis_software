#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ep -> en pi+ asymmetries versus -t' in several xB bins, for three run periods
and for the combined (inverse-variance weighted) file.

Also prints LaTeX tables (to stdout) for the 4 xB slices:
  enpiLowxBGE, enpiMidLowxBGE, enpiMidHighxBGE, enpiHighxBGE

Tables are printed only when BOTH --rad and --migration are provided, and use:
  central = (fit value) + (radiative Delta)
  stat    = (fit statistical uncertainty)
  syst    = sqrt( sigma_rad^2 + sigma_mig^2 )

Rows are ordered in increasing -t' (i.e. increasing x axis as plotted).

Usage:
  python plot_rgc_epni+_asymmetry.py Su22.txt Fa22.txt Sp23.txt Combined.txt --rad ISR_FSR_delta_summary.txt --migration migration_summary.txt
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
YLIM_UU = (-1.0, 1.0)    # UU updated to match DSA

XB_BINS = {
    "enpiLowxBGE":     r"$0.10 < x_{B} < 0.25$",
    "enpiMidLowxBGE":  r"$0.25 < x_{B} < 0.35$",
    "enpiMidHighxBGE": r"$0.35 < x_{B} < 0.45$",
    "enpiHighxBGE":    r"$0.45 < x_{B} < 0.60$",
    "enpiGE":          r"$0.10 < x_{B} < 0.60$"
}

BIN_ORDER = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE", "enpiGE"]

TOP_PAD_PER_BIN = 0.95
TOP_PAD_OVERLAY = 0.94

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
        "UUcos":   get_series(parsed, k("AUUcosphi")),
        "UUcos2":  get_series(parsed, k("AUUcos2phi")),
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
    Returns:
      rad[bin_tag][series_key] = {"x": array(-t'), "delta": array, "sigma": array}
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
    Migration *_delta list: (mean_tprime, delta)
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

def build_migration_dict(mig_parsed, bin_prefix):
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
    if series is None:
        return None
    #endif
    x, y, yerr = series["x"], series["y"], series["yerr"]
    if rad_entry is None:
        return x, y, yerr
    #endif

    xr = np.asarray(rad_entry["x"], dtype=float)
    dr = np.asarray(rad_entry["delta"], dtype=float)

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

# ---------------------------------------------------------------------
# Legends
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
# Plotting (2x2 canvases)
# ---------------------------------------------------------------------
def _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label,
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

    # BSA (ALU sin phi)
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

    # TSA (AUL sin phi open, sin2 phi filled)
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

    # DSA (ALL n=0 open, cos phi filled)
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

    # UU
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

    # BSA
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

    # TSA (systematics from AULsin)
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

    # DSA (systematics from ALLn0)
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

    # UU
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
    IMPORTANT: wrap the full cell in $...$ to ensure math mode, since we use ^ and _.
    This prevents 'Missing $ inserted.' in LaTeX tables.
    """
    fmt = "{0:." + str(ndp) + "f}"
    core = fmt.format(val) + "^{\\pm " + fmt.format(stat) + "}_{\\pm " + fmt.format(syst) + "}"
    return "$" + core + "$"
#endif

def print_latex_table_for_bin(bin_tag, comb_parsed, rad_all, mig_all):
    """
    Prints one LaTeX table for a given bin_tag using Combined fit results,
    radiative shifts, and quadrature systematics.
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

    # Prepare arrays: central shifts from rad delta, sys from quadrature of sigmas
    vals = {}
    stats = {}
    systs = {}

    for key in required:
        y = np.asarray(p[key]["y"], dtype=float)
        ye = np.asarray(p[key]["yerr"], dtype=float)

        d_rad, s_rad = _get_rad_arrays(rad_bin, key, n)
        s_mig = _get_mig_arrays(mig_bin, key, n)

        y_shift = y + d_rad
        s_tot = np.sqrt(s_rad * s_rad + s_mig * s_mig)

        vals[key] = y_shift
        stats[key] = ye
        systs[key] = s_tot
    #endfor

    # Header strings (exactly as requested)
    header = (
        "\\begin{table}[h]\n"
        "\\centering\n"
        "\\small\n"
        "\\begin{tabular}{|c|c|c|c|c|c|} \\hline\n"
        "$\\langle -t' \\rangle$ & $F_{LU}^{\\sin\\phi}/F_{UU}$ & $F_{UL}^{\\sin\\phi}/F_{UU}$ & $F_{UL}^{\\sin2\\phi}/F_{UU}$ & $F_{LL}/F_{UU}$ & $F_{LL}^{\\cos\\phi}/F_{UU}$ \\\\ \\hline\n"
    )

    # Caption/label: deterministic
    label = f"table:{bin_tag}_fitresults_tprime"
    caption = (
        "\\end{tabular}\n"
        "\\caption{Fitted structure-function ratios per bin. Entries are $\\text{value}^{\\pm\\,\\text{stat}}_{\\pm\\,\\text{syst}}$.\\label{"
        + label +
        "}}\n"
        "\\end{table}\n"
    )

    # Print table
    print("")
    print("% -----------------------------------------------------------------------------")
    print(f"% LaTeX table for bin: {bin_tag}  ({XB_BINS.get(bin_tag, bin_tag)})")
    print("% Central values include radiative Delta shifts; syst is sqrt(rad^2 + mig^2).")
    print("% Rows are ordered in increasing -t'.")
    print("% -----------------------------------------------------------------------------")
    print(header, end="")

    # Rows (already sorted increasing -t' due to get_series sort)
    for i in range(n):
        tmean = float(x[i])

        c1 = _fmt_cell(float(vals["ALUsin"][i]),  float(stats["ALUsin"][i]),  float(systs["ALUsin"][i]),  ndp=3)
        c2 = _fmt_cell(float(vals["AULsin"][i]),  float(stats["AULsin"][i]),  float(systs["AULsin"][i]),  ndp=3)
        c3 = _fmt_cell(float(vals["AULsin2"][i]), float(stats["AULsin2"][i]), float(systs["AULsin2"][i]), ndp=3)
        c4 = _fmt_cell(float(vals["ALLn0"][i]),   float(stats["ALLn0"][i]),   float(systs["ALLn0"][i]),   ndp=3)
        c5 = _fmt_cell(float(vals["ALLcos"][i]),  float(stats["ALLcos"][i]),  float(systs["ALLcos"][i]),  ndp=3)

        row = (
            f"{tmean:.3f} ~&~ {c1} ~&~ {c2} ~&~ {c3} ~&~ {c4} ~&~ {c5} \\\\ \\hline"
        )
        print(row)
    #endfor

    print(caption, end="")

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Plot enpi+ asymmetries vs -t' with optional radiative shifts and bin-migration systematics.")
    ap.add_argument("su22", type=str, help="Su22 fit-results .txt")
    ap.add_argument("fa22", type=str, help="Fa22 fit-results .txt")
    ap.add_argument("sp23", type=str, help="Sp23 fit-results .txt")
    ap.add_argument("combined", type=str, help="Combined fit-results .txt")
    ap.add_argument("--rad", type=str, default=None, help="Signed Delta summary .txt to apply (+Delta)")
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
            raise RuntimeError(f"FATAL: --rad file not found: {args.rad}")
        #endif
        rad_all = parse_rad_summary(args.rad)
    #endif

    mig_all = None
    if args.migration:
        if not os.path.isfile(args.migration):
            raise RuntimeError(f"FATAL: --migration file not found: {args.migration}")
        #endif
        mig_parsed = parse_migration_file(args.migration)
        mig_all = {}
        for bin_tag in BIN_ORDER:
            mig_all[bin_tag] = build_migration_dict(mig_parsed, bin_tag) or {}
        #endfor
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