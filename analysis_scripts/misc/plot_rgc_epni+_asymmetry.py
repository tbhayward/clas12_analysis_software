#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ep -> en pi+ asymmetries versus -t' in several xB bins, for three run periods
and for the combined (inverse-variance weighted) file.

Usage:
  python plot_enpi_from_texts.py Su22.txt Fa22.txt Sp23.txt Combined.txt [--rad signed_delta_summary.txt]

Outputs (for each available xB bin):
  output/enpi+/rgc_enpi_<BinTag>_AllPeriods.pdf
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly.pdf

If --rad is provided, also:
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly_withRad.pdf

Overlay (Combined-only, 1×3):
  output/enpi+/rgc_enpi_xBOverlay_CombinedOnly.pdf
  (and if --rad) output/enpi+/rgc_enpi_xBOverlay_CombinedOnly_withRad.pdf

Notes:
- Fit-result text files must contain blocks like:
    enpiLowxBGEchi2FitsALUsinphi = {{mean_tprime, value, error}, ...};
  x is plotted as -t'.
- With --rad, the central values are shifted by the **signed** Δ from the summary;
  stat errors remain unchanged. σΔ is shown as a shaded "systematic" rectangle
  from y=0 to y=σΔ at each x. (We deliberately do not center the band on the point.)
- Systematic-band widths are fixed so bars **touch**:
    * enpiGE (integrated): width = 0.10 in -t'
    * xB sub-bins:        width = 0.20 in -t'
- TSA systematics use only F_{UL}^{sinφ}; sin2φ systematics are ignored.
- DSA systematics use only F_{LL}; F_{LL}^{cosφ} systematics are ignored.
- On the 1×3 overlay, systematic bands (if --rad) are taken **only** from the
  lowest xB bin (enpiLowxBGE).
"""

import sys
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# ─────────────────────────────────────────────────────────────────────
# Styling / knobs
# ─────────────────────────────────────────────────────────────────────
COLORS = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green",
    "Combined": "black",
}

# Colors for the four xB slices on the 1×3 overlay
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

YLIM_LU = (-0.2, 0.2)   # BSA
YLIM_UL = (-0.2, 0.2)   # TSA
YLIM_LL = (-0.4, 0.4)   # DSA
YLIM_UU = (-0.4, 0.4)   # UU

# Human-readable labels for xB bins (for legends & titles)
XB_BINS = {
    "enpiLowxBGE":     r"$0.10 < x_{B} < 0.25$",
    "enpiMidLowxBGE":  r"$0.25 < x_{B} < 0.35$",
    "enpiMidHighxBGE": r"$0.35 < x_{B} < 0.45$",
    "enpiHighxBGE":    r"$0.45 < x_{B} < 0.60$",
    "enpiGE":          r"$0.10 < x_{B} < 0.60$"  # inclusive slice
}

# Order to render canvases
BIN_ORDER = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE", "enpiGE"]

# Top padding
TOP_PAD_PER_BIN = 0.95   # for 2×2 canvases with a suptitle
TOP_PAD_OVERLAY = 0.94   # for 1×3 overlay (no visible title)

# Output prefix constant
OUT_PREFIX = "enpi"

# ─────────────────────────────────────────────────────────────────────
# Parsing helpers (fit-result .txt files)
# ─────────────────────────────────────────────────────────────────────
_assign_re = re.compile(r'([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};', re.DOTALL)
_triple_re = re.compile(r'\{([^{}]+)\}')

def parse_asym_file(path):
    """Parse NAME = {{mean,val,err}, ...}; blocks into dict[name] -> np.array(N,3)."""
    with open(path, "r") as f:
        text = f.read()
    out = {}
    for m in _assign_re.finditer(text):
        name = m.group(1)
        content = m.group(2)
        triples = []
        for t in _triple_re.findall(content):
            parts = [p.strip() for p in t.split(",")]
            if len(parts) != 3:
                continue
            try:
                triples.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
        if triples:
            out[name] = np.array(triples, dtype=float)
    return out

def get_series(dct, key, negate_x=True, sort_x=True):
    """
    Return dict(x,y,yerr) if key exists with finite values & positive errors; else None.
    Negates x if negate_x=True (to convert t' -> -t' for plotting). Optionally sort by x.
    """
    if key not in dct:
        return None
    arr = np.array(dct[key], dtype=float)
    if arr.size == 0:
        return None
    x_raw, y, e = arr[:,0], arr[:,1], arr[:,2]
    mask = np.isfinite(x_raw) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    if not np.any(mask):
        return None
    x = -x_raw[mask] if negate_x else x_raw[mask]
    y = y[mask]
    e = e[mask]
    if sort_x and x.size > 1:
        order = np.argsort(x)
        x, y, e = x[order], y[order], e[order]
    return {"x": x, "y": y, "yerr": e}

def build_period_dict(parsed, bin_prefix):
    """
    Collect the series for one period given a bin_prefix like 'enpiLowxBGE' or 'enpiGE'.
    We look for keys of the form f\"{bin_prefix}chi2Fits<suffix>\".
    """
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}"
    return {
        # BSA
        "ALUsin":  get_series(parsed, k("ALUsinphi")),
        # TSA (n=1 open, n=2 filled)
        "AULsin":  get_series(parsed, k("AULsinphi")),
        "AULsin2": get_series(parsed, k("AULsin2phi")),
        # DSA (n=0 open, n=1 (cosφ) filled)
        "ALLn0":   get_series(parsed, k("ALL")),
        "ALLcos":  get_series(parsed, k("ALLcosphi")),
        # UU (n=1 open, n=2 filled)
        "UUcos":   get_series(parsed, k("AUUcosphi")),
        "UUcos2":  get_series(parsed, k("AUUcos2phi")),
    }

def detect_available_bins(*dicts):
    """
    Inspect provided parsed dictionaries and return the subset of BIN_ORDER whose
    expected keys exist in ANY of the dicts. We probe ALUsinphi to decide presence.
    """
    available = []
    for bin_tag in BIN_ORDER:
        key_probe = f"{bin_tag}chi2FitsALUsinphi"
        present = any((d is not None and key_probe in d) for d in dicts)
        if present:
            available.append(bin_tag)
    return available

# ─────────────────────────────────────────────────────────────────────
# Radiative-correction summary parser (from signed Δ summary txt)
# ─────────────────────────────────────────────────────────────────────
# We map series by a distinctive LaTeX substring present in the summary file.
_SERIES_KEYWORDS = {
    "ALUsin":   r"F_{LU}^{\sin\phi}/F_{UU}",
    "AULsin":   r"F_{UL}^{\sin\phi}/F_{UU}",
    "AULsin2":  r"F_{UL}^{\sin2\phi}/F_{UU}",
    "ALLn0":    r"F_{LL}/F_{UU}",
    "ALLcos":   r"F_{LL}^{\cos\phi}/F_{UU}",
    # "Dilution" exists too, but not used here
}

def parse_rad_summary(rad_path):
    """
    Parse the signed-Δ summary text into:
        rad[bin_tag][series_key] = dict(x=array, delta=array, sigma=array)
    We expect blocks like:

      Bin: enpiGE  x_B range: ...
      Series: F_{LU}^{\sin\phi}/F_{UU}
        -t'        Δ          σΔ
        0.10     -0.0123     0.0045
        ...

    Anything outside those patterns is ignored.
    """
    rad = {}
    cur_bin = None
    cur_series_key = None
    with open(rad_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # New bin section
        if line.startswith("Bin:"):
            parts = line.split()
            if len(parts) >= 2:
                cur_bin = parts[1]
                rad.setdefault(cur_bin, {})
            cur_series_key = None
            i += 1
            continue

        # New series section
        if line.startswith("Series:"):
            ser_line = line[len("Series:"):].strip()
            cur_series_key = None
            for key, tag in _SERIES_KEYWORDS.items():
                if tag in ser_line:
                    cur_series_key = key
                    break
            # skip header line (column names), if present
            i += 1
            if i < len(lines) and any(h in lines[i] for h in ["-t'", "Δ", "sigma", "σ"]):
                i += 1
            xs, deltas, sigmas = [], [], []
            while i < len(lines):
                row = lines[i].strip()
                if (not row) or row.startswith("Series:") or row.startswith("Bin:") or set(row) <= set("=-"):
                    # end of this series block
                    break
                parts = row.split()
                # Expect 3 numeric columns: x  Δ  σ
                try:
                    xval = float(parts[0])
                    dval = float(parts[1])
                    sval = float(parts[2])
                except Exception:
                    # stop if malformed
                    break
                xs.append(xval)
                deltas.append(dval)
                sigmas.append(sval)
                i += 1

            if cur_bin is not None and cur_series_key is not None and xs:
                # Convert to arrays and store (x is already -t')
                arr = {
                    "x": np.array(xs, dtype=float),
                    "delta": np.array(deltas, dtype=float),
                    "sigma": np.array(sigmas, dtype=float),
                }
                rad[cur_bin][cur_series_key] = arr
            continue

        i += 1

    return rad

# ─────────────────────────────────────────────────────────────────────
# Titles (2×2 canvases only)
# ─────────────────────────────────────────────────────────────────────
def make_title(xb_label, with_rad=False):
    """Compose suptitle string."""
    base = r"$ep \rightarrow en\pi^{+}$ — " + xb_label
    if with_rad:
        return base + r" (with radiative shifts $\,+\Delta$ and bands $\sigma_{\Delta}$)"
    return base

# ─────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────
def _legend_run_period(ax, labels):
    handles = [Line2D([0],[0], marker=MARKER, color=COLORS[L], linestyle='', label=L) for L in labels]
    leg = ax.legend(handles=handles, title="Run Period", frameon=True, edgecolor="black",
                    loc="lower right", fontsize=11, title_fontsize=12)
    leg.get_frame().set_alpha(0.9)

def _legend_harmonic(ax, labels=("n=1","n=2"), where="lower left"):
    open_marker = Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label=labels[0])
    filled_marker = Line2D([0],[0], marker=MARKER, color='black', linestyle='', label=labels[1])
    leg = ax.legend(handles=[open_marker, filled_marker], title="Harmonic",
                    frameon=True, edgecolor="black", loc=where,
                    bbox_to_anchor=(0.02, 0.02), fontsize=11, title_fontsize=12)
    ax.add_artist(leg)

def _bin_sys_width(bin_tag: str) -> float:
    """
    Fixed systematic-band width (so bars touch):
      - enpiGE (integrated): 0.10
      - all xB sub-bins:     0.20
    """
    return 0.10 if bin_tag == "enpiGE" else 0.20

def _draw_sys_bands(ax, x, sigma, *, bin_width=None, facecolor="0.7", alpha=0.4, zorder=1):
    """
    Draw 'systematic' rectangles for each x-bin:
      rectangle spans [x - bin_width/2, x + bin_width/2] in x,
      and [0, sigma] in y (not centered at the point).

    If bin_width is None, a small fallback width is used (not recommended).
    """
    if x is None or sigma is None:
        return
    x = np.asarray(x, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    if x.size == 0 or sigma.size == 0:
        return

    w = 0.05 if (bin_width is None or bin_width <= 0) else float(bin_width)

    for xv, sv in zip(x, sigma):
        if not np.isfinite(xv) or not np.isfinite(sv):
            continue
        rect = Rectangle((xv - 0.5*w, 0.0), w, sv, facecolor=facecolor, edgecolor=None,
                         alpha=alpha, zorder=zorder)
        ax.add_patch(rect)

def _apply_shift_if_available(series, rad_entry):
    """
    Align by **index order** (no nearest-neighbor matching).
    Returns (x, y_shift, yerr, sigma_rad) or (x, y, yerr, None) if no rad_entry.
    """
    if series is None:
        return None
    x, y, yerr = series["x"], series["y"], series["yerr"]

    if rad_entry is None:
        return x, y, yerr, None

    xr = np.asarray(rad_entry["x"])
    dr = np.asarray(rad_entry["delta"])
    sr = np.asarray(rad_entry["sigma"])

    n = min(len(x), len(xr), len(dr), len(sr))
    if n == 0:
        return x, y, yerr, None
    if len(x) != len(xr):
        print(f"[WARN] Δ length ({len(xr)}) != data length ({len(x)}); pairing by index up to {n}.")
    y_shift = y.copy()
    y_shift[:n] = y_shift[:n] + dr[:n]
    sigma_out = np.full_like(y, np.nan, dtype=float)
    sigma_out[:n] = sr[:n]
    return x, y_shift, yerr, sigma_out

def _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label,
                     rad_for_bin=None, with_rad=False, *, bin_tag_for_width="enpiGE"):
    """
    Draw the four panels for a set of labeled period dicts.
    If with_rad=True and rad_for_bin is provided, shift central values by +Δ (signed)
    and draw σΔ bands from y=0 to y=σΔ (one common set per panel, independent of run period).

    Systematics:
      - TSA: use only F_{UL}^{sinφ} (ignore sin2φ systematics)
      - DSA: use only F_{LL}      (ignore F_{LL}^{cosφ} systematics)
    """
    bw = _bin_sys_width(bin_tag_for_width)

    # Helper: get the rad entry for a given series key
    def R(key):
        if rad_for_bin is None:
            return None
        return rad_for_bin.get(key)

    # ---- BSA (ALU sinφ) ----
    if with_rad and R("ALUsin") is not None:
        _draw_sys_bands(axLU, R("ALUsin")["x"], np.abs(R("ALUsin")["sigma"]), bin_width=bw)
    for lab, pdata in pdata_by_label.items():
        s = pdata["ALUsin"]
        if s is None:
            continue
        x, y, ye, _ = _apply_shift_if_available(s, R("ALUsin")) if with_rad else (s["x"], s["y"], s["yerr"], None)
        axLU.errorbar(x, y, yerr=ye, fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                      capsize=CAPSIZE, label=lab, markersize=MS, linestyle="None")
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)
    _legend_run_period(axLU, list(pdata_by_label.keys()))

    # ---- TSA (AUL sinφ open, AUL sin2φ filled) ----
    # Systematic bands: **only** from AULsin (ignore sin2 systematics)
    if with_rad and R("AULsin") is not None:
        _draw_sys_bands(axUL, R("AULsin")["x"], np.abs(R("AULsin")["sigma"]), bin_width=bw)
    for lab, pdata in pdata_by_label.items():
        s1 = pdata["AULsin"]
        if s1 is not None:
            x, y, ye, _ = _apply_shift_if_available(s1, R("AULsin")) if with_rad else (s1["x"], s1["y"], s1["yerr"], None)
            axUL.errorbar(x, y, yerr=ye, fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        s2 = pdata["AULsin2"]
        if s2 is not None:
            if with_rad and R("AULsin2") is not None:
                x2, y2, ye2, _ = _apply_shift_if_available(s2, R("AULsin2"))
            else:
                x2, y2, ye2 = s2["x"], s2["y"], s2["yerr"]
            axUL.errorbar(x2, y2, yerr=ye2, fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axUL, labels=("n=1","n=2"))
    _legend_run_period(axUL, list(pdata_by_label.keys()))

    # ---- DSA (ALL n=0 open, cosφ filled) ----
    # Systematic bands: **only** from ALL (F_{LL}), ignore F_{LL}^{cosφ} systematics
    if with_rad and R("ALLn0") is not None:
        _draw_sys_bands(axLL, R("ALLn0")["x"], np.abs(R("ALLn0")["sigma"]), bin_width=bw)
    for lab, pdata in pdata_by_label.items():
        s0 = pdata["ALLn0"]
        if s0 is not None:
            x, y, ye, _ = _apply_shift_if_available(s0, R("ALLn0")) if with_rad else (s0["x"], s0["y"], s0["yerr"], None)
            axLL.errorbar(x, y, yerr=ye, fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        s1 = pdata["ALLcos"]
        if s1 is not None:
            x1, y1, ye1, _ = _apply_shift_if_available(s1, R("ALLcos")) if with_rad else (s1["x"], s1["y"], s1["yerr"], None)
            axLL.errorbar(x1, y1, yerr=ye1, fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axLL, labels=("n=0","n=1"))
    _legend_run_period(axLL, list(pdata_by_label.keys()))

    # ---- UU (cos nφ: n=1 open, n=2 filled) ----
    for lab, pdata in pdata_by_label.items():
        s1 = pdata["UUcos"]
        if s1 is not None:
            axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"], fmt="o", mfc="none",
                          mec=COLORS[lab], ecolor=COLORS[lab], capsize=CAPSIZE,
                          markersize=MS, linestyle="None")
        s2 = pdata["UUcos2"]
        if s2 is not None:
            axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"], fmt="o", color=COLORS[lab],
                          ecolor=COLORS[lab], capsize=CAPSIZE, markersize=MS, linestyle="None")
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axUU, labels=("n=1","n=2"))
    _legend_run_period(axUU, list(pdata_by_label.keys()))

# ─────────────────────────────────────────────────────────────────────
# Per-bin canvases (2×2)
# ─────────────────────────────────────────────────────────────────────
def plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, out_dir, with_rad=False, rad_bin=None):
    plt.figure(figsize=(12, 9))  # 2×2
    axLU  = plt.subplot(2,2,1)   # BSA
    axUL  = plt.subplot(2,2,2)   # TSA
    axLL  = plt.subplot(2,2,3)   # DSA
    axUU  = plt.subplot(2,2,4)   # UU

    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(xb_label, with_rad=with_rad), fontsize=16, y=0.97)

    pdata_by_label = {"Su22": p_su22, "Fa22": p_fa22, "Sp23": p_sp23}
    _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label,
                     rad_for_bin=rad_bin, with_rad=with_rad, bin_tag_for_width=bin_tag)

    plt.tight_layout(rect=[0, 0, 1, TOP_PAD_PER_BIN])
    os.makedirs(out_dir, exist_ok=True)
    suffix = "_withRad" if with_rad else ""
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_{bin_tag}_AllPeriods{suffix}.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved all-periods figure: {out_path}")

def plot_combined_only_for_bin(p_comb, bin_tag, out_dir, with_rad=False, rad_bin=None):
    plt.figure(figsize=(12, 9))  # 2×2
    axLU  = plt.subplot(2,2,1)
    axUL  = plt.subplot(2,2,2)
    axLL  = plt.subplot(2,2,3)
    axUU  = plt.subplot(2,2,4)

    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(xb_label, with_rad=with_rad), fontsize=16, y=0.97)

    black = COLORS["Combined"]
    bw = _bin_sys_width(bin_tag)

    # Convenience accessor for rad
    def R(key):
        if rad_bin is None:
            return None
        return rad_bin.get(key)

    # BSA
    if with_rad and R("ALUsin") is not None:
        _draw_sys_bands(axLU, R("ALUsin")["x"], np.abs(R("ALUsin")["sigma"]), bin_width=bw)
    if p_comb["ALUsin"] is not None:
        s = p_comb["ALUsin"]
        if with_rad:
            x, y, ye, _ = _apply_shift_if_available(s, R("ALUsin"))
        else:
            x, y, ye = s["x"], s["y"], s["yerr"]
        axLU.errorbar(x, y, yerr=ye, fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)

    # TSA (bands from AULsin only)
    if with_rad and R("AULsin") is not None:
        _draw_sys_bands(axUL, R("AULsin")["x"], np.abs(R("AULsin")["sigma"]), bin_width=bw)
    if p_comb["AULsin"] is not None:
        s1 = p_comb["AULsin"]
        if with_rad:
            x1, y1, ye1, _ = _apply_shift_if_available(s1, R("AULsin"))
        else:
            x1, y1, ye1 = s1["x"], s1["y"], s1["yerr"]
        axUL.errorbar(x1, y1, yerr=ye1, fmt=MARKER, mfc="none", mec=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    if p_comb["AULsin2"] is not None:
        s2 = p_comb["AULsin2"]
        if with_rad and R("AULsin2") is not None:
            x2, y2, ye2, _ = _apply_shift_if_available(s2, R("AULsin2"))
        else:
            x2, y2, ye2 = s2["x"], s2["y"], s2["yerr"]
        axUL.errorbar(x2, y2, yerr=ye2, fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    # small harmonic legend
    _leg = axUL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(_leg)

    # DSA (bands from ALL only)
    if with_rad and R("ALLn0") is not None:
        _draw_sys_bands(axLL, R("ALLn0")["x"], np.abs(R("ALLn0")["sigma"]), bin_width=bw)
    if p_comb["ALLn0"] is not None:
        s0 = p_comb["ALLn0"]
        if with_rad:
            x0, y0, ye0, _ = _apply_shift_if_available(s0, R("ALLn0"))
        else:
            x0, y0, ye0 = s0["x"], s0["y"], s0["yerr"]
        axLL.errorbar(x0, y0, yerr=ye0, fmt=MARKER, mfc="none", mec=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    if p_comb["ALLcos"] is not None:
        s1 = p_comb["ALLcos"]
        if with_rad and R("ALLcos") is not None:
            x1, y1, ye1, _ = _apply_shift_if_available(s1, R("ALLcos"))
        else:
            x1, y1, ye1 = s1["x"], s1["y"], s1["yerr"]
        axLL.errorbar(x1, y1, yerr=ye1, fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)
    _leg2 = axLL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=1'),
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axLL.add_artist(_leg2)

    # UU (no rad shifts expected)
    if p_comb["UUcos"] is not None:
        s1 = p_comb["UUcos"]
        axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"], fmt=MARKER, mfc="none", mec=black,
                      ecolor=black, capsize=CAPSIZE, markersize=MS, linestyle="None")
    if p_comb["UUcos2"] is not None:
        s2 = p_comb["UUcos2"]
        axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"], fmt=MARKER, color=black, ecolor=black,
                      capsize=CAPSIZE, markersize=MS, linestyle="None")
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout(rect=[0, 0, 1, TOP_PAD_PER_BIN])
    os.makedirs(out_dir, exist_ok=True)
    suffix = "_withRad" if with_rad else ""
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_{bin_tag}_CombinedOnly{suffix}.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved combined-only figure: {out_path}")

# ─────────────────────────────────────────────────────────────────────
# xB overlay 1×3 canvas (Combined only) — NO suptitle
# ─────────────────────────────────────────────────────────────────────
def plot_combined_xb_overlay_1x3(comb_parsed, out_dir, bins_to_use,
                                 with_rad=False, rad_all=None):
    """
    1×3 canvas (Combined only) overlaying the four xB slices:
      Left  : F_LU^{sinφ}/F_UU  (ALUsinphi)
      Middle: F_UL^{sinφ}/F_UU  (AULsinphi)
      Right : F_UL^{sin2φ}/F_UU (AULsin2phi)

    If with_rad=True:
      * central values are shifted by +Δ (per series, per bin if available),
      * systematic bands are drawn **only from enpiLowxBGE** uncertainties,
        using bin width 0.20 so the bars touch,
      * TSA uses only F_{UL}^{sinφ} bands (ignore sin2φ bands).
    """
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharex=True)
    axL, axM, axR = axes

    low_bin_tag = "enpiLowxBGE"
    bw = _bin_sys_width(low_bin_tag)  # 0.20 on overlay

    def _draw_component(ax, suffix_key, ylabel, ylim, draw_sys=True, ignore_sys=False):
        handles = []

        # Bands from the **lowest xB bin only** if requested and available
        if with_rad and rad_all is not None and draw_sys and not ignore_sys:
            low_rad = rad_all.get(low_bin_tag, {})
            rad_entry = low_rad.get(suffix_key)
            if rad_entry is not None:
                _draw_sys_bands(ax, rad_entry["x"], np.abs(rad_entry["sigma"]), bin_width=bw)

        # Now the points for each xB slice
        for bin_tag in bins_to_use:
            pdata = build_period_dict(comb_parsed, bin_tag)
            s = pdata[suffix_key]
            if s is None:
                continue
            clr = XB_COLORS.get(bin_tag, "black")
            mkr = XB_MARKERS.get(bin_tag, "o")
            if with_rad and rad_all is not None:
                rad_bin = rad_all.get(bin_tag, {})
                rad_entry = rad_bin.get(suffix_key)
                if rad_entry is not None:
                    x, y, ye, _ = _apply_shift_if_available(s, rad_entry)
                else:
                    x, y, ye = s["x"], s["y"], s["yerr"]
            else:
                x, y, ye = s["x"], s["y"], s["yerr"]
            h = ax.errorbar(x, y, yerr=ye, fmt=mkr, color=clr, ecolor=clr, capsize=CAPSIZE,
                            label=XB_BINS.get(bin_tag, bin_tag), markersize=MS, linestyle="None")
            handles.append(h)

        ax.set(xlim=XLIM_T, ylim=ylim, xlabel=X_LABEL, ylabel=ylabel)
        ax.axhline(0, color="black", linestyle="--", linewidth=1.2)
        ax.grid(True, linestyle="--", alpha=0.6)
        if handles:
            leg = ax.legend(frameon=True, edgecolor="black", fontsize=11, loc="best")
            leg.get_frame().set_alpha(0.9)

    # Left: ALUsin (draw bands from low xB)
    _draw_component(axL, "ALUsin",  r"$F_{LU}^{\sin\phi}/F_{UU}$", YLIM_LU, draw_sys=True,  ignore_sys=False)
    # Middle: AULsin (draw bands from low xB; ignore sin2 here by definition)
    _draw_component(axM, "AULsin",  r"$F_{UL}^{\sin\phi}/F_{UU}$", YLIM_UL, draw_sys=True,  ignore_sys=False)
    # Right: AULsin2 (plot points; **no** systematic bands per instruction to ignore sin2φ systematics)
    _draw_component(axR, "AULsin2", r"$F_{UL}^{\sin2\phi}/F_{UU}$", YLIM_UL, draw_sys=False, ignore_sys=True)

    fig.tight_layout(rect=[0, 0, 1, TOP_PAD_OVERLAY])
    os.makedirs(out_dir, exist_ok=True)
    suffix = "_withRad" if with_rad else ""
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_xBOverlay_CombinedOnly{suffix}.pdf")
    plt.savefig(out_path)
    plt.close(fig)
    print(f"Saved xB-overlay (1×3) figure: {out_path}")

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="Plot enpi+ asymmetries vs -t' with optional radiative-correction shifts.")
    ap.add_argument("su22", type=str, help="Su22 fit-results .txt")
    ap.add_argument("fa22", type=str, help="Fa22 fit-results .txt")
    ap.add_argument("sp23", type=str, help="Sp23 fit-results .txt")
    ap.add_argument("combined", type=str, help="Combined fit-results .txt")
    ap.add_argument("--rad", type=str, default=None, help="Signed Δ summary .txt to apply (+Δ) with σΔ bands")
    args = ap.parse_args()

    # Parse fit-result files
    su22 = parse_asym_file(args.su22)
    fa22 = parse_asym_file(args.fa22)
    sp23 = parse_asym_file(args.sp23)
    comb = parse_asym_file(args.combined)

    # Which xB-bin canvases can we make?
    available_bins = detect_available_bins(su22, fa22, sp23, comb)
    if not available_bins:
        print("[ERROR] No recognizable xB-bin sections found (e.g. enpiLowxBGEchi2FitsALUsinphi).")
        sys.exit(2)

    # Optional radiative-correction summary
    rad_all = None
    if args.rad:
        if not os.path.isfile(args.rad):
            print(f"[WARN] --rad file not found: {args.rad}. Proceeding without radiative shifts.")
        else:
            rad_all = parse_rad_summary(args.rad)

    out_dir = os.path.join("output", "enpi+")
    os.makedirs(out_dir, exist_ok=True)

    # Per-bin canvases
    for bin_tag in available_bins:
        p_su22 = build_period_dict(su22, bin_tag)
        p_fa22 = build_period_dict(fa22, bin_tag)
        p_sp23 = build_period_dict(sp23, bin_tag)
        p_comb = build_period_dict(comb,  bin_tag)

        # Always make the nominal "AllPeriods" canvas (3 runs)
        plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, out_dir, with_rad=False, rad_bin=None)

        # Combined-only (nominal)
        plot_combined_only_for_bin(p_comb, bin_tag, out_dir, with_rad=False, rad_bin=None)

        # If --rad is present, ONLY make the Combined-only radiative canvas (no AllPeriods_withRad)
        if rad_all is not None:
            rad_bin = rad_all.get(bin_tag, {})
            plot_combined_only_for_bin(p_comb, bin_tag, out_dir, with_rad=True, rad_bin=rad_bin)

    # Combined-only xB overlay (1×3)
    xb_overlay_candidates = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE"]
    bins_to_use = [b for b in xb_overlay_candidates if b in available_bins]
    if bins_to_use:
        # Nominal
        plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use, with_rad=False, rad_all=None)
        # With radiative shifts (if present); sys bands from lowest xB only
        if rad_all is not None:
            plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use, with_rad=True, rad_all=rad_all)
    else:
        print("[WARN] No xB-slice data available for overlay canvas; skipping.")

if __name__ == "__main__":
    main()