#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ep -> en pi+ asymmetries versus -t in several xB bins, for three run periods
and for the combined (inverse-variance weighted) file.

Usage:
  python plot_enpi_from_texts.py Su22.txt Fa22.txt Sp23.txt Combined.txt <Prefix> "<Kinematic text>"

Examples:
  python plot_enpi_from_texts.py su22.txt fa22.txt sp23.txt combined.txt enpiGE \
    "Q^2>1, W>2, y<0.75"

Saves for each xB bin (Low, MidLow, MidHigh, High, Inclusive):
  output/enpi+/rgc_<PREFIX>_<BinTag>_AllPeriods.pdf
  output/enpi+/rgc_<PREFIX>_<BinTag>_CombinedOnly.pdf

Additional xB-overlay canvas (1×3, Combined only):
  output/enpi+/rgc_<PREFIX>_xBOverlay_CombinedOnly.pdf

Notes:
- The text files are expected to contain sections like:
    enpiLowxBGEchi2FitsALUsinphi = {{mean_t, value, error}, ...};
  where 'mean_t' is the (negative) mean of t in that bin. We plot against -t.
- We auto-detect which xB bin tags are available by inspecting keys.
- Titles no longer include a fixed Mx2 exclusivity window; those are now applied bin-by-bin upstream.
"""

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ─────────────────────────────────────────────────────────────────────
# Styling
# ─────────────────────────────────────────────────────────────────────
COLORS = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green",
    "Combined": "black",
}

# Colors used to distinguish the four xB slices on the new 1×3 overlay
XB_COLORS = {
    "enpiLowxBGE":     "tab:blue",
    "enpiMidLowxBGE":  "tab:orange",
    "enpiMidHighxBGE": "tab:green",
    "enpiHighxBGE":    "tab:red",
}

# NEW: marker shapes per xB slice for the 1×3 overlay
XB_MARKERS = {
    "enpiLowxBGE":     "o",  # circle
    "enpiMidLowxBGE":  "s",  # square
    "enpiMidHighxBGE": "^",  # triangle up
    "enpiHighxBGE":    "D",  # diamond
}

MARKER = "o"
CAPSIZE = 3
LABEL_FONTSIZE = 13
MS = 5.0  # marker size

# x-axis and labels (now -t)
XLIM_T = (0.0, 1.30)
X_LABEL = r"$-t\ (\mathrm{GeV}^{2})$"

# y-axis ranges per panel
YLIM_LU = (-0.2, 0.2)   # BSA
YLIM_UL = (-0.2, 0.2)   # TSA
YLIM_LL = (-1.0, 1.0)   # DSA
YLIM_UU = (-1.0, 1.0)   # UU

# Human-readable labels for xB bins (shown in titles / legends)
XB_BINS = {
    "enpiLowxBGE":     r"$0.10 < x_{B} < 0.25$",
    "enpiMidLowxBGE":  r"$0.25 < x_{B} < 0.35$",
    "enpiMidHighxBGE": r"$0.35 < x_{B} < 0.45$",
    "enpiHighxBGE":    r"$0.45 < x_{B} < 0.60$",
    "enpiGE":          r"$0.10 < x_{B} < 0.60$"  # inclusive slice
}

# Order in which to render canvases
BIN_ORDER = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE", "enpiGE"]

# ─────────────────────────────────────────────────────────────────────
# Parsing helpers
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
    Negates the x-values if negate_x=True (to convert t-> -t for plotting).
    Optionally sorts by x ascending for nicer plotting.
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
    We look for keys of the form f"{bin_prefix}chi2Fits{suffix}".
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
    Inspect provided parsed dictionaries and return the subset of BIN_ORDER
    whose expected keys are present in ANY of the dicts.
    We require at minimum that the ALUsin key exists for that bin in some file.
    """
    available = []
    for bin_tag in BIN_ORDER:
        key_probe = f"{bin_tag}chi2FitsALUsinphi"
        present = any((key_probe in d) for d in dicts if d is not None)
        if present:
            available.append(bin_tag)
    return available

# ─────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────
def make_title(kin_text, xb_label):
    # Titles no longer include a fixed Mx2 exclusivity window.
    common = r"$Q^{2}>1\ (\mathrm{GeV}^{2}),\ W>2\ (\mathrm{GeV}),\ y<0.75$"
    if kin_text.strip():
        return rf"$ep \rightarrow en\pi^{{+}}$ — {xb_label}; {common}; {kin_text}"
    else:
        return rf"$ep \rightarrow en\pi^{{+}}$ — {xb_label}; {common}"

def _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label):
    """Internal: draw the four panels for a set of labeled period dicts."""
    # ---- BSA (ALU sinφ) ----
    for lab, pdata in pdata_by_label.items():
        s = pdata["ALUsin"]
        if s is not None:
            axLU.errorbar(s["x"], s["y"], yerr=s["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=lab, markersize=MS, linestyle="None")
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)
    legLU = axLU.legend(title="Run Period", frameon=True, edgecolor="black",
                        fontsize=11, title_fontsize=12, loc="best")
    legLU.get_frame().set_alpha(0.9)

    # ---- TSA (AUL sinφ open, sin2φ filled) ----
    for lab, pdata in pdata_by_label.items():
        s1 = pdata["AULsin"]
        if s1 is not None:
            axUL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1", markersize=MS, linestyle="None")
        s2 = pdata["AULsin2"]
        if s2 is not None:
            axUL.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                          fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=2", markersize=MS, linestyle="None")
    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    harmUL = axUL.legend(
        handles=[
            Line2D([0],[0], marker="o", mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker="o", color='black', linestyle='', label='n=2')
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(harmUL)
    labs = list(pdata_by_label.keys())
    run_leg_handles = [Line2D([0],[0], marker="o", color=COLORS[L], linestyle='', label=L) for L in labs]
    legUL = axUL.legend(handles=run_leg_handles, title="Run Period", frameon=True, edgecolor="black",
                        loc="lower right", fontsize=11, title_fontsize=12)
    legUL.get_frame().set_alpha(0.9)

    # ---- DSA (ALL n=0 open, cosφ filled) ----
    for lab, pdata in pdata_by_label.items():
        s0 = pdata["ALLn0"]
        if s0 is not None:
            axLL.errorbar(s0["x"], s0["y"], yerr=s0["yerr"],
                          fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=0", markersize=MS, linestyle="None")
        s1 = pdata["ALLcos"]
        if s1 is not None:
            axLL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1", markersize=MS, linestyle="None")
    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)
    harmLL = axLL.legend(
        handles=[
            Line2D([0],[0], marker="o", mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker="o", color='black', linestyle='', label='n=1')
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axLL.add_artist(harmLL)
    labs = list(pdata_by_label.keys())
    run_leg_handles = [Line2D([0],[0], marker="o", color=COLORS[L], linestyle='', label=L) for L in labs]
    legLL = axLL.legend(handles=run_leg_handles, title="Run Period", frameon=True, edgecolor="black",
                        loc="lower right", fontsize=11, title_fontsize=12)
    legLL.get_frame().set_alpha(0.9)

    # ---- UU (cos nφ) - n=1 open, n=2 filled) ----
    for lab, pdata in pdata_by_label.items():
        s1 = pdata["UUcos"]
        if s1 is not None:
            axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1", markersize=MS, linestyle="None")
        s2 = pdata["UUcos2"]
        if s2 is not None:
            axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                          fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=2", markersize=MS, linestyle="None")
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)
    harmUU = axUU.legend(
        handles=[
            Line2D([0],[0], marker="o", mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker="o", color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axUU.add_artist(harmUU)
    labs = list(pdata_by_label.keys())
    run_leg_handles = [Line2D([0],[0], marker="o", color=COLORS[L], linestyle='', label=L) for L in labs]
    legUU = axUU.legend(handles=run_leg_handles, title="Run Period", frameon=True, edgecolor="black",
                        loc="lower right", fontsize=11, title_fontsize=12)
    legUU.get_frame().set_alpha(0.9)

def plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, kin_text, out_dir, out_prefix):
    plt.figure(figsize=(12, 9))
    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(kin_text, xb_label), fontsize=16, y=0.97)

    # Order: BSA TL, TSA TR, DSA BL, UU BR
    axLU  = plt.subplot(2,2,1)  # BSA
    axUL  = plt.subplot(2,2,2)  # TSA
    axLL  = plt.subplot(2,2,3)  # DSA
    axUU  = plt.subplot(2,2,4)  # UU

    pdata_by_label = {
        "Su22": p_su22,
        "Fa22": p_fa22,
        "Sp23": p_sp23,
    }
    _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{out_prefix}_{bin_tag}_AllPeriods.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved all-periods figure: {out_path}")

def plot_combined_only_for_bin(p_comb, bin_tag, kin_text, out_dir, out_prefix):
    plt.figure(figsize=(12, 9))
    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(kin_text, xb_label), fontsize=16, y=0.97)

    # Panels
    axLU  = plt.subplot(2,2,1)
    axUL  = plt.subplot(2,2,2)
    axLL  = plt.subplot(2,2,3)
    axUU  = plt.subplot(2,2,4)

    black = COLORS["Combined"]

    # BSA
    if p_comb["ALUsin"] is not None:
        s = p_comb["ALUsin"]
        axLU.errorbar(s["x"], s["y"], yerr=s["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE,
                      markersize=MS, linestyle="None")
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)

    # TSA
    if p_comb["AULsin"] is not None:
        s1 = p_comb["AULsin"]
        axUL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE,
                      markersize=MS, linestyle="None")
    if p_comb["AULsin2"] is not None:
        s2 = p_comb["AULsin2"]
        axUL.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE,
                      markersize=MS, linestyle="None")
    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    harmUL = axUL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(harmUL)

    # DSA
    if p_comb["ALLn0"] is not None:
        s0 = p_comb["ALLn0"]
        axLL.errorbar(s0["x"], s0["y"], yerr=s0["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE,
                      markersize=MS, linestyle="None")
    if p_comb["ALLcos"] is not None:
        s1 = p_comb["ALLcos"]
        axLL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE,
                      markersize=MS, linestyle="None")
    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)
    harmLL = axLL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=1'),
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axLL.add_artist(harmLL)

    # UU
    if p_comb["UUcos"] is not None:
        s1 = p_comb["UUcos"]
        axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE,
                      markersize=MS, linestyle="None")
    if p_comb["UUcos2"] is not None:
        s2 = p_comb["UUcos2"]
        axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE,
                      markersize=MS, linestyle="None")
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)
    harmUU = axUU.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        fontsize=11, title_fontsize=12
    )
    axUU.add_artist(harmUU)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{out_prefix}_{bin_tag}_CombinedOnly.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved combined-only figure: {out_path}")

# ─────────────────────────────────────────────────────────────────────
# New: xB overlay 1×3 canvas from the Combined file
# ─────────────────────────────────────────────────────────────────────
def plot_combined_xb_overlay_1x3(comb_parsed, kin_text, out_dir, out_prefix, bins_to_use):
    """
    Build a 1×3 canvas (Combined only) where each subplot overlays the four xB slices:
      Left  : F_LU^{sinφ}/F_UU  (ALUsinphi)
      Middle: F_UL^{sinφ}/F_UU  (AULsinphi)
      Right : F_UL^{sin2φ}/F_UU (AULsin2phi)

    UPDATE: each xB slice gets a unique marker shape (see XB_MARKERS) in addition to color.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharex=True)
    fig.suptitle(
        r"$ep \rightarrow en\pi^{+}$ — $x_{B}$ slices overlay; "
        r"$Q^{2}>1,\ W>2,\ y<0.75$"
        + (("; " + kin_text) if kin_text.strip() else ""),
        fontsize=15, y=0.98
    )

    axL, axM, axR = axes

    # Helper to draw one component on given axis for all xB bins
    def _draw_component(ax, suffix_key, ylabel, ylim):
        handles = []
        labels  = []
        for bin_tag in bins_to_use:
            pdata = build_period_dict(comb_parsed, bin_tag)
            s = pdata[suffix_key]
            if s is None:
                continue
            clr = XB_COLORS.get(bin_tag, "black")
            mkr = XB_MARKERS.get(bin_tag, "o")
            h = ax.errorbar(
                s["x"], s["y"], yerr=s["yerr"],
                fmt=mkr, color=clr, ecolor=clr, capsize=CAPSIZE,
                label=XB_BINS.get(bin_tag, bin_tag),
                markersize=MS, linestyle="None"
            )
            handles.append(h)
            labels.append(XB_BINS.get(bin_tag, bin_tag))

        ax.set(xlim=XLIM_T, ylim=ylim, xlabel=X_LABEL, ylabel=ylabel)
        ax.axhline(0, color="black", linestyle="--", linewidth=1.2)
        ax.grid(True, linestyle="--", alpha=0.6)
        if handles:
            leg = ax.legend(frameon=True, edgecolor="black", fontsize=11, loc="best")
            leg.get_frame().set_alpha(0.9)

    _draw_component(axL, "ALUsin",  r"$F_{LU}^{\sin\phi}/F_{UU}$", YLIM_LU)
    _draw_component(axM, "AULsin",  r"$F_{UL}^{\sin\phi}/F_{UU}$", YLIM_UL)
    _draw_component(axR, "AULsin2", r"$F_{UL}^{\sin2\phi}/F_{UU}$", YLIM_UL)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{out_prefix}_xBOverlay_CombinedOnly.pdf")
    plt.savefig(out_path)
    plt.close(fig)
    print(f"Saved xB-overlay (1×3) figure: {out_path}")

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 7:
        print("Usage: python plot_enpi_from_texts.py Su22.txt Fa22.txt Sp23.txt Combined.txt <Prefix> \"<Kinematic text>\"")
        print("Example:")
        print("  python plot_enpi_from_texts.py su22.txt fa22.txt sp23.txt combined.txt enpiGE \"Q^2>1, W>2, y<0.75\"")
        sys.exit(1)

    su22_path, fa22_path, sp23_path, comb_path, out_prefix, kin_text = sys.argv[1:7]

    # Parse
    su22 = parse_asym_file(su22_path)
    fa22 = parse_asym_file(fa22_path)
    sp23 = parse_asym_file(sp23_path)
    comb = parse_asym_file(comb_path)

    # Figure out which xB-bin canvases we can actually make
    available_bins = detect_available_bins(su22, fa22, sp23, comb)
    if not available_bins:
        print("[ERROR] No recognizable xB-bin sections found (e.g. enpiLowxBGEchi2FitsALUsinphi).")
        sys.exit(2)

    out_dir = os.path.join("output", "enpi+")
    os.makedirs(out_dir, exist_ok=True)

    # Build and plot per bin
    for bin_tag in available_bins:
        p_su22 = build_period_dict(su22, bin_tag)
        p_fa22 = build_period_dict(fa22, bin_tag)
        p_sp23 = build_period_dict(sp23, bin_tag)
        p_comb = build_period_dict(comb,  bin_tag)

        # All periods on one canvas
        plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, kin_text, out_dir, out_prefix)

        # Combined-only canvas
        plot_combined_only_for_bin(p_comb, bin_tag, kin_text, out_dir, out_prefix)

    # New: single xB overlay 1×3 canvas from the Combined file.
    xb_overlay_candidates = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE"]
    bins_to_use = [b for b in xb_overlay_candidates if b in available_bins]
    if bins_to_use:
        plot_combined_xb_overlay_1x3(comb, kin_text, out_dir, out_prefix, bins_to_use)
    else:
        print("[WARN] No xB-slice data available for overlay canvas; skipping.")

if __name__ == "__main__":
    main()