#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_period_study_enpi_asymmetries.py

Purpose:
  Make clean, easy-to-read comparison plots of ep -> en pi+ asymmetry fit results
  across run periods (Su22, Fa22, Sp23) plus the Combined (inverse-variance weighted)
  file, with NO radiative corrections and NO migration corrections.

Inputs (runtime args):
  1) Su22 fit-results .txt
  2) Fa22 fit-results .txt
  3) Sp23 fit-results .txt
  4) Combined fit-results .txt

Input format expectations (same as your existing script):
  Files contain assignment blocks like:
    enpiLowxBGEchi2FitsALUsinphi = {{mean_tprime, value, error}, ...};
  where mean_tprime is typically negative and we plot x = -mean_tprime (positive -t').

Plots:
  For each available xB bin tag, produce a 2x3 canvas:
    (0,0) ALUsin    : F_LU^{sin phi}/F_UU
    (0,1) AULsin    : F_UL^{sin phi}/F_UU
    (0,2) AULsin2   : F_UL^{sin 2phi}/F_UU
    (1,0) ALLn0     : F_LL/F_UU
    (1,1) ALLcos    : F_LL^{cos phi}/F_UU
    (1,2) legend pad (no data)

  Each subplot overlays:
    Su22 (color)
    Fa22 (color)
    Sp23 (color)
    Combined (black)

Output:
  Saves one PDF per bin_tag to:
    output/enpi+/run_period_study/

Example:
  python run_period_study_enpi_asymmetries.py Su22.txt Fa22.txt Sp23.txt Combined.txt
"""

import sys
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Styling / knobs
# -----------------------------------------------------------------------------
COLORS = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green",
    "Combined": "black",
}

MARKER = "o"
CAPSIZE = 3
MS = 5.0
ELINEWIDTH = 1.0

XLIM_T = (0.0, 1.30)
X_LABEL = r"$-t'\ (\mathrm{GeV}^{2})$"

# Per-series y limits (kept consistent with your prior choices, but split per panel).
YLIMS = {
    "ALUsin":  (-0.4, 0.4),
    "AULsin":  (-0.4, 0.4),
    "AULsin2": (-0.4, 0.4),
    "ALLn0":   (-1.0, 1.0),
    "ALLcos":  (-1.0, 1.0),
}

# Display labels for xB bins (ASCII-only strings inside Python; LaTeX is ASCII).
XB_BINS = {
    "enpiLowxBGE":     r"$0.10 < x_{B} < 0.25$",
    "enpiMidLowxBGE":  r"$0.25 < x_{B} < 0.35$",
    "enpiMidHighxBGE": r"$0.35 < x_{B} < 0.45$",
    "enpiHighxBGE":    r"$0.45 < x_{B} < 0.60$",
    "enpiGE":          r"$0.10 < x_{B} < 0.60$",
}

BIN_ORDER = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE", "enpiGE"]

OUT_DIR = os.path.join("output", "enpi+", "run_period_study")

# -----------------------------------------------------------------------------
# Parsing helpers (fit-result .txt files)
# -----------------------------------------------------------------------------
_assign_re = re.compile(r"([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};", re.DOTALL)
_triple_re = re.compile(r"\{([^{}]+)\}")

def parse_asym_file(path):
    """
    Parse NAME = {{mean,val,err}, ...}; blocks into dict[name] -> np.array(N,3).
    """
    if not os.path.isfile(path):
        raise RuntimeError(f"FATAL: input file not found: {path}")
    #endif

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
    """
    Build the 5 structure-function ratios for a given bin_prefix.
    Keys here match the canonical internal names used across your scripts.
    """
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}"
    return {
        "ALUsin":  get_series(parsed, k("ALUsinphi")),
        "AULsin":  get_series(parsed, k("AULsinphi")),
        "AULsin2": get_series(parsed, k("AULsin2phi")),
        "ALLn0":   get_series(parsed, k("ALL")),
        "ALLcos":  get_series(parsed, k("ALLcosphi")),
    }

def detect_available_bins(*dicts):
    """
    Determine which bin tags exist in at least one file by probing a canonical key.
    """
    available = []
    for bin_tag in BIN_ORDER:
        key_probe = f"{bin_tag}chi2FitsALUsinphi"
        present = any((d is not None and key_probe in d) for d in dicts)
        if present:
            available.append(bin_tag)
        #endif
    #endfor
    return available

# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------
def make_title(bin_tag):
    xb_label = XB_BINS.get(bin_tag, bin_tag)
    return r"$ep \rightarrow en\pi^{+}$" + "  " + xb_label

def _style_axis(ax, ylabel, ylim):
    ax.set(xlim=XLIM_T, ylim=ylim, xlabel=X_LABEL, ylabel=ylabel)
    ax.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
    ax.grid(True, linestyle="--", alpha=0.6)

def _plot_one_series(ax, series_key, ylabel, period_dicts_by_label):
    """
    Plot one series (e.g. ALUsin) overlaying all periods including Combined.
    period_dicts_by_label maps:
      label -> dict(series_key -> {x,y,yerr})
    """
    for lab, pdata in period_dicts_by_label.items():
        s = pdata.get(series_key)
        if s is None:
            continue
        #endif

        ax.errorbar(
            s["x"], s["y"], yerr=s["yerr"],
            fmt=MARKER,
            color=COLORS[lab],
            ecolor=COLORS[lab],
            capsize=CAPSIZE,
            markersize=MS,
            elinewidth=ELINEWIDTH,
            linestyle="None",
            label=lab
        )
    #endfor

    _style_axis(ax, ylabel=ylabel, ylim=YLIMS[series_key])

def plot_bin_2x3(bin_tag, p_su22, p_fa22, p_sp23, p_comb, out_dir):
    """
    Produce the 2x3 canvas for one bin_tag.
    """
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(make_title(bin_tag), fontsize=16, y=0.98)

    period_dicts_by_label = {
        "Su22": p_su22,
        "Fa22": p_fa22,
        "Sp23": p_sp23,
        "Combined": p_comb,
    }

    # Row 0
    _plot_one_series(
        axes[0, 0],
        "ALUsin",
        ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    _plot_one_series(
        axes[0, 1],
        "AULsin",
        ylabel=r"$F_{UL}^{\sin\phi}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    _plot_one_series(
        axes[0, 2],
        "AULsin2",
        ylabel=r"$F_{UL}^{\sin2\phi}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    # Row 1
    _plot_one_series(
        axes[1, 0],
        "ALLn0",
        ylabel=r"$F_{LL}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    _plot_one_series(
        axes[1, 1],
        "ALLcos",
        ylabel=r"$F_{LL}^{\cos\phi}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    # Legend pad (axes[1,2])
    ax_leg = axes[1, 2]
    ax_leg.axis("off")

    handles = []
    labels = []
    for lab in ["Su22", "Fa22", "Sp23", "Combined"]:
        handles.append(Line2D([0], [0], marker=MARKER, color=COLORS[lab], linestyle="", markersize=MS))
        labels.append(lab)
    #endfor

    ax_leg.legend(
        handles=handles,
        labels=labels,
        title="Run Period",
        loc="center",
        frameon=True,
        edgecolor="black",
        fontsize=12,
        title_fontsize=13
    )

    # Tight layout with a little room for the title
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    os.makedirs(out_dir, exist_ok=True)

    out_path = os.path.join(out_dir, f"rgc_enpi_{bin_tag}_RunPeriodStudy.pdf")
    plt.savefig(out_path)
    plt.close(fig)

    print(f"Saved: {out_path}")

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Run-period study: compare Su22/Fa22/Sp23 vs Combined for ep->en pi+ structure-function ratios (no rad/migration)."
    )
    ap.add_argument("su22", type=str, help="Su22 fit-results .txt")
    ap.add_argument("fa22", type=str, help="Fa22 fit-results .txt")
    ap.add_argument("sp23", type=str, help="Sp23 fit-results .txt")
    ap.add_argument("combined", type=str, help="Combined fit-results .txt")
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

    os.makedirs(OUT_DIR, exist_ok=True)

    for bin_tag in available_bins:
        p_su22 = build_period_dict(su22, bin_tag)
        p_fa22 = build_period_dict(fa22, bin_tag)
        p_sp23 = build_period_dict(sp23, bin_tag)
        p_comb = build_period_dict(comb, bin_tag)

        # If a bin exists but is missing everything, skip it loudly but non-fatally.
        have_any = False
        for k in ["ALUsin", "AULsin", "AULsin2", "ALLn0", "ALLcos"]:
            if (p_su22.get(k) is not None) or (p_fa22.get(k) is not None) or (p_sp23.get(k) is not None) or (p_comb.get(k) is not None):
                have_any = True
                break
            #endif
        #endfor

        if not have_any:
            print(f"[WARN] Skipping bin '{bin_tag}': no plottable series found in any input file.")
            continue
        #endif

        plot_bin_2x3(bin_tag, p_su22, p_fa22, p_sp23, p_comb, OUT_DIR)
    #endfor

if __name__ == "__main__":
    main()
# endif