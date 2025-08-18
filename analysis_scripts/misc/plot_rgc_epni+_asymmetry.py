#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ep -> en pi+ asymmetries versus -t in several xB bins, for three run periods
and for the combined (inverse-variance weighted) file.

Usage (kept the same CLI):
  python plot_enpi_from_texts.py Su22.txt Fa22.txt Sp23.txt Combined.txt <Prefix> "<Kinematic text>"

Examples:
  python plot_enpi_from_texts.py su22.txt fa22.txt sp23.txt combined.txt enpiGE \
    "Q^2>1, W>2, y<0.75, 0.81<M_x^2<1.00 (GeV^2)"

Saves for each xB bin (Low, MidLow, MidHigh, High, Inclusive):
  output/enpi+/rgc_<PREFIX>_<BinTag>_AllPeriods.pdf
  output/enpi+/rgc_<PREFIX>_<BinTag>_CombinedOnly.pdf

Notes:
- The text files are expected to contain sections like:
    enpiLowxBGEchi2FitsALUsinphi = {{mean_t, value, error}, ...};
  where 'mean_t' is the (negative) mean of t in that bin. We plot against -t.
- We auto-detect which xB bin tags are available by inspecting keys.
- Common cuts are encoded upstream; here we only plot and annotate titles.
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
MARKER = "o"
CAPSIZE = 3
LABEL_FONTSIZE = 13

# x-axis and labels (now -t)
XLIM_T = (0.0, 1.30)
X_LABEL = r"$-t\ (\mathrm{GeV}^{2})$"

# y-axis ranges per panel (same as your previous choices)
YLIM_LU = (-0.4, 0.4)   # BSA
YLIM_UL = (-0.4, 0.4)   # TSA
YLIM_LL = (-1.0, 1.0)   # DSA
YLIM_UU = (-1.0, 1.0)   # UU

# Human-readable labels for xB bins (shown in titles)
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
        text = f.read()  # endif
    out = {}
    for m in _assign_re.finditer(text):
        name = m.group(1)
        content = m.group(2)
        triples = []
        for t in _triple_re.findall(content):
            parts = [p.strip() for p in t.split(",")]
            if len(parts) != 3:
                continue  # endif
            try:
                triples.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue  # endif
        # endfor
        if triples:
            out[name] = np.array(triples, dtype=float)
    # endfor
    return out  # endif
# endfor

def get_series(dct, key, negate_x=True, sort_x=True):
    """
    Return dict(x,y,yerr) if key exists with finite values & positive errors; else None.
    Negates the x-values if negate_x=True (to convert t-> -t for plotting).
    Optionally sorts by x ascending for nicer plotting.
    """
    if key not in dct:
        return None  # endif
    arr = np.array(dct[key], dtype=float)
    if arr.size == 0:
        return None  # endif
    x_raw, y, e = arr[:,0], arr[:,1], arr[:,2]
    mask = np.isfinite(x_raw) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    if not np.any(mask):
        return None  # endif
    x = -x_raw[mask] if negate_x else x_raw[mask]
    y = y[mask]
    e = e[mask]
    if sort_x and x.size > 1:
        order = np.argsort(x)
        x, y, e = x[order], y[order], e[order]
    # endif
    return {"x": x, "y": y, "yerr": e}  # endif
# endfor

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
    }  # endif
# endfor

def detect_available_bins(*dicts):
    """
    Inspect the provided parsed dictionaries and return the subset of BIN_ORDER
    whose expected keys are present in ANY of the dicts.
    We require at minimum that the ALUsin key exists for that bin in some file.
    """
    available = []
    for bin_tag in BIN_ORDER:
        key_probe = f"{bin_tag}chi2FitsALUsinphi"
        present = any((key_probe in d) for d in dicts if d is not None)
        if present:
            available.append(bin_tag)
        # endif
    # endfor
    return available  # endif
# endfor

# ─────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────
def make_title(kin_text, xb_label):
    # Reaction in math text, then cuts, then xB slice label
    # We also include the common cuts explicitly for clarity.
    common = r"$Q^{2}>1 (\mathrm{GeV}^{2}),\ W>2 (\mathrm{GeV}),\ y<0.75,\ 0.81<M_{x}^{2}<1.00\ (\mathrm{GeV}^{2})$"
    # Escape braces around B for format if used elsewhere; here it's raw string already.
    if kin_text.strip():
        return rf"$ep \rightarrow en\pi^{{+}}$ — {xb_label}; {common}; {kin_text}"
    else:
        return rf"$ep \rightarrow en\pi^{{+}}$ — {xb_label}; {common}"
    # endif
# endfor

def _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label):
    """Internal: draw the four panels for a set of labeled period dicts."""
    # ---- BSA (ALU sinφ) ----
    for lab, pdata in pdata_by_label.items():
        s = pdata["ALUsin"]
        if s is not None:
            axLU.errorbar(s["x"], s["y"], yerr=s["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=lab)
        # endif
    # endfor
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
                          fmt=MARKER, mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1")
        # endif
        s2 = pdata["AULsin2"]
        if s2 is not None:
            axUL.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=2")
        # endif
    # endfor
    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    harmUL = axUL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2')
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(harmUL)
    # add separate run period legend
    labs = list(pdata_by_label.keys())
    run_leg_handles = [Line2D([0],[0], marker=MARKER, color=COLORS[L], linestyle='', label=L) for L in labs]
    legUL = axUL.legend(handles=run_leg_handles, title="Run Period", frameon=True, edgecolor="black",
                        loc="lower right", fontsize=11, title_fontsize=12)
    legUL.get_frame().set_alpha(0.9)

    # ---- DSA (ALL n=0 open, cosφ filled) ----
    for lab, pdata in pdata_by_label.items():
        s0 = pdata["ALLn0"]
        if s0 is not None:
            axLL.errorbar(s0["x"], s0["y"], yerr=s0["yerr"],
                          fmt=MARKER, mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=0")
        # endif
        s1 = pdata["ALLcos"]
        if s1 is not None:
            axLL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1")
        # endif
    # endfor
    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)
    harmLL = axLL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=1')
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axLL.add_artist(harmLL)
    labs = list(pdata_by_label.keys())
    run_leg_handles = [Line2D([0],[0], marker=MARKER, color=COLORS[L], linestyle='', label=L) for L in labs]
    legLL = axLL.legend(handles=run_leg_handles, title="Run Period", frameon=True, edgecolor="black",
                        loc="lower right", fontsize=11, title_fontsize=12)
    legLL.get_frame().set_alpha(0.9)

    # ---- UU (cos nφ) - n=1 open, n=2 filled) ----
    for lab, pdata in pdata_by_label.items():
        s1 = pdata["UUcos"]
        if s1 is not None:
            axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt=MARKER, mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1")
        # endif
        s2 = pdata["UUcos2"]
        if s2 is not None:
            axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=2")
        # endif
    # endfor
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)
    harmUU = axUU.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUU.add_artist(harmUU)
    labs = list(pdata_by_label.keys())
    run_leg_handles = [Line2D([0],[0], marker=MARKER, color=COLORS[L], linestyle='', label=L) for L in labs]
    legUU = axUU.legend(handles=run_leg_handles, title="Run Period", frameon=True, edgecolor="black",
                        loc="lower right", fontsize=11, title_fontsize=12)
    legUU.get_frame().set_alpha(0.9)
# endfor

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
    print(f"Saved all-periods figure: {out_path}")  # endif
# endfor

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
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    # endif
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)

    # TSA
    if p_comb["AULsin"] is not None:
        s1 = p_comb["AULsin"]
        axUL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE)
    # endif
    if p_comb["AULsin2"] is not None:
        s2 = p_comb["AULsin2"]
        axUL.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    # endif
    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    harmUL = axUL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(harmUL)

    # DSA
    if p_comb["ALLn0"] is not None:
        s0 = p_comb["ALLn0"]
        axLL.errorbar(s0["x"], s0["y"], yerr=s0["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE)
    # endif
    if p_comb["ALLcos"] is not None:
        s1 = p_comb["ALLcos"]
        axLL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    # endif
    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)
    harmLL = axLL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=1'),
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axLL.add_artist(harmLL)

    # UU
    if p_comb["UUcos"] is not None:
        s1 = p_comb["UUcos"]
        axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE)
    # endif
    if p_comb["UUcos2"] is not None:
        s2 = p_comb["UUcos2"]
        axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    # endif
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)
    harmUU = axUU.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUU.add_artist(harmUU)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{out_prefix}_{bin_tag}_CombinedOnly.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved combined-only figure: {out_path}")  # endif
# endfor

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 7:
        print("Usage: python plot_enpi_from_texts.py Su22.txt Fa22.txt Sp23.txt Combined.txt <Prefix> \"<Kinematic text>\"")
        print("Example:")
        print("  python plot_enpi_from_texts.py su22.txt fa22.txt sp23.txt combined.txt enpiGE \"Q^2>1, W>2, y<0.75, 0.81<M_x^2<1.00 (GeV^2)\"")
        sys.exit(1)  # endif

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
        sys.exit(2)  # endif

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
    # endfor
# endfor

if __name__ == "__main__":
    main()
# endif