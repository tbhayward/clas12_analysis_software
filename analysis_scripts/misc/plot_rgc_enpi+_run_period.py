#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_period_study_enpi_asymmetries.py

Purpose:
  Make clean comparison plots of ep -> en pi+ asymmetry fit results
  across run periods (Su22, Fa22, Sp23), with NO radiative corrections
  and NO migration corrections.

Enhancement:
  Adds a per-subplot statistical compatibility test for the THREE run periods
  (Su22, Fa22, Sp23), reporting a weighted chi-square heterogeneity statistic
  (Cochran Q) and its p-value.

  IMPORTANT CHANGE (per user request):
    We assume the series lists are already aligned in the correct bin order.
    Therefore, we pair points strictly by index:
      Su22[i] <-> Fa22[i] <-> Sp23[i]
    even if the mean -t' values differ slightly between periods.

Compatibility test definition:
  For each i (index-paired point), with three measurements y_p,i and stat errors sigma_p,i:
    w_p,i = 1/sigma_p,i^2
    mu_i  = sum_p w_p,i y_p,i / sum_p w_p,i
    Q     = sum_i sum_p w_p,i (y_p,i - mu_i)^2
    dof   = (k - 1) * N_used = 2 * N_used   (since k=3 periods)

  Under the null hypothesis that the three periods are statistically consistent
  given the stated statistical errors, Q approximately follows chi2(dof).
  The p-value is P(Chi2(dof) >= Q).

Interpretation guidance (practical):
  - Large p (e.g. p > 0.10): no evidence of incompatibility (differences look consistent with stat errors).
  - Moderate p (0.01 to 0.10): some tension; worth checking systematics, run conditions, or outliers.
  - Small p (p < 0.01): strong evidence the periods are not mutually consistent under stat-only errors.
  - Q/dof ~ 1 is "about right" for a good stat-only description. Much larger suggests extra variance.

Inputs (runtime args):
  1) Su22 fit-results .txt
  2) Fa22 fit-results .txt
  3) Sp23 fit-results .txt

Input format expectations:
  Files contain assignment blocks like:
    enpiLowxBGEchi2FitsALUsinphi = {{mean_tprime, value, error}, ...};
  We plot x = -mean_tprime (positive -t').

Plots:
  For each available xB bin tag, produce a 2x3 canvas:
    (0,0) ALUsin    : F_LU^{sin phi}/F_UU
    (0,1) AULsin    : F_UL^{sin phi}/F_UU
    (0,2) AULsin2   : F_UL^{sin 2phi}/F_UU
    (1,0) ALLn0     : F_LL/F_UU
    (1,1) ALLcos    : F_LL^{cos phi}/F_UU
    (1,2) legend pad (no data)

Output:
  Saves one PDF per bin_tag to:
    output/enpi+/run_period_study/

Example:
  python run_period_study_enpi_asymmetries.py Su22.txt Fa22.txt Sp23.txt
"""

import sys
import os
import re
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Chi-square survival function helper (p-value)
#   Prefer SciPy; fall back to mpmath; else fail fast.
# -----------------------------------------------------------------------------
def chi2_sf(x, dof):
    """
    Survival function P(Chi2(dof) >= x).
    Requires either scipy.stats or mpmath.
    """
    try:
        from scipy.stats import chi2 as _chi2
        return float(_chi2.sf(x, dof))
    except Exception:
        pass
    #endif

    try:
        import mpmath as mp
        a = 0.5 * float(dof)
        z = 0.5 * float(x)
        val = mp.gammainc(a, z, mp.inf) / mp.gamma(a)
        return float(val)
    except Exception:
        raise RuntimeError(
            "FATAL: To compute p-values you need either SciPy or mpmath installed.\n"
            "       Try one of:\n"
            "         pip install scipy\n"
            "         pip install mpmath"
        )
    #endif

# -----------------------------------------------------------------------------
# Styling / knobs
# -----------------------------------------------------------------------------
COLORS = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green",
}

MARKER = "o"
CAPSIZE = 3
MS = 5.0
ELINEWIDTH = 1.0

XLIM_T = (0.0, 1.30)
X_LABEL = r"$-t'\ (\mathrm{GeV}^{2})$"

# Per-series y limits
YLIMS = {
    "ALUsin":  (-0.4, 0.4),
    "AULsin":  (-0.4, 0.4),
    "AULsin2": (-0.4, 0.4),
    "ALLn0":   (-1.0, 1.0),
    "ALLcos":  (-1.0, 1.0),
}

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

def get_series(dct, key, negate_x=True, sort_x=False):
    """
    Return dict(x,y,yerr) if key exists with finite values and positive errors; else None.

    NOTE:
      We do NOT sort by x because the user wants index-pairing across periods.
      We preserve the file's intrinsic ordering.
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

    # Preserve order; just drop invalid rows in-place.
    x_raw = x_raw[mask]
    y = y[mask]
    e = e[mask]

    x = -x_raw if negate_x else x_raw

    if sort_x and x.size > 1:
        order = np.argsort(x)
        x, y, e = x[order], y[order], e[order]
    #endif

    return {"x": x, "y": y, "yerr": e}

def build_period_dict(parsed, bin_prefix):
    """
    Build the 5 structure-function ratios for a given bin_prefix.
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
# Statistical compatibility test across Su22/Fa22/Sp23 (index paired)
# -----------------------------------------------------------------------------
def cochran_q_three_index_paired(s_su22, s_fa22, s_sp23):
    """
    Compute Cochran Q heterogeneity statistic for exactly three periods using index pairing.

    Returns:
      (Q, dof, pval, Q_over_dof, N_used)

    We use:
      N0 = min(len(Su22), len(Fa22), len(Sp23))
      Then we keep only indices i in [0, N0) where all y and yerr are finite and yerr>0.

    If N_used == 0: returns (None, None, None, None, 0).
    """
    if (s_su22 is None) or (s_fa22 is None) or (s_sp23 is None):
        return None, None, None, None, 0
    #endif

    y1 = np.asarray(s_su22["y"], dtype=float)
    e1 = np.asarray(s_su22["yerr"], dtype=float)
    y2 = np.asarray(s_fa22["y"], dtype=float)
    e2 = np.asarray(s_fa22["yerr"], dtype=float)
    y3 = np.asarray(s_sp23["y"], dtype=float)
    e3 = np.asarray(s_sp23["yerr"], dtype=float)

    N0 = int(min(y1.size, y2.size, y3.size, e1.size, e2.size, e3.size))
    if N0 <= 0:
        return None, None, None, None, 0
    #endif

    y1 = y1[:N0]; e1 = e1[:N0]
    y2 = y2[:N0]; e2 = e2[:N0]
    y3 = y3[:N0]; e3 = e3[:N0]

    mask = (
        np.isfinite(y1) & np.isfinite(e1) & (e1 > 0.0) &
        np.isfinite(y2) & np.isfinite(e2) & (e2 > 0.0) &
        np.isfinite(y3) & np.isfinite(e3) & (e3 > 0.0)
    )

    if not np.any(mask):
        return None, None, None, None, 0
    #endif

    y1 = y1[mask]; e1 = e1[mask]
    y2 = y2[mask]; e2 = e2[mask]
    y3 = y3[mask]; e3 = e3[mask]

    N = int(y1.size)
    if N <= 0:
        return None, None, None, None, 0
    #endif

    Q = 0.0
    for i in range(N):
        w1 = 1.0 / float(e1[i] * e1[i])
        w2 = 1.0 / float(e2[i] * e2[i])
        w3 = 1.0 / float(e3[i] * e3[i])

        wsum = w1 + w2 + w3
        if (not math.isfinite(wsum)) or (wsum <= 0.0):
            continue
        #endif

        mu = (w1 * float(y1[i]) + w2 * float(y2[i]) + w3 * float(y3[i])) / wsum
        Q += w1 * (float(y1[i]) - mu) * (float(y1[i]) - mu)
        Q += w2 * (float(y2[i]) - mu) * (float(y2[i]) - mu)
        Q += w3 * (float(y3[i]) - mu) * (float(y3[i]) - mu)
    #endfor

    dof = 2 * N
    if dof <= 0:
        return None, None, None, None, N
    #endif

    pval = chi2_sf(Q, dof)
    Q_over_dof = Q / float(dof)
    return Q, dof, pval, Q_over_dof, N

def format_pval(p):
    if p is None or (not math.isfinite(p)):
        return "n/a"
    #endif
    if p < 1.0e-3:
        return "<1e-3"
    #endif
    return f"{p:.3f}"

def add_test_text(ax, Q, dof, pval, Q_over_dof, N_used):
    """
    Add a compact annotation to the axis in axes coordinates.
    """
    if (Q is None) or (dof is None) or (pval is None) or (Q_over_dof is None) or (N_used is None) or (N_used <= 0):
        txt = "Period compatibility:\n n/a"
    else:
        txt = (
            "Su22 vs Fa22 vs Sp23:\n"
            "Index paired\n"
            f"N = {int(N_used):d}\n"
            f"Q = {Q:.1f}, dof = {int(dof):d}\n"
            f"p = {format_pval(pval)}, Q/dof = {Q_over_dof:.2f}"
        )
    #endif

    ax.text(
        0.02, 0.98, txt,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=10,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.80, pad=3.0)
    )

# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------
def make_title(bin_tag):
    xb_label = XB_BINS.get(bin_tag, bin_tag)
    return r"$ep \rightarrow en\pi^{+}$" + "  " + xb_label

def style_axis(ax, ylabel, ylim):
    ax.set(xlim=XLIM_T, ylim=ylim, xlabel=X_LABEL, ylabel=ylabel)
    ax.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
    ax.grid(True, linestyle="--", alpha=0.6)

def plot_one_series(ax, series_key, ylabel, period_dicts_by_label):
    """
    Plot one series overlaying Su22/Fa22/Sp23,
    and annotate a compatibility test across those three (index paired).
    """
    # Plot points
    for lab in ["Su22", "Fa22", "Sp23"]:
        pdata = period_dicts_by_label.get(lab, {})
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

    style_axis(ax, ylabel=ylabel, ylim=YLIMS[series_key])

    # Compatibility test (index paired)
    s_su22 = period_dicts_by_label.get("Su22", {}).get(series_key)
    s_fa22 = period_dicts_by_label.get("Fa22", {}).get(series_key)
    s_sp23 = period_dicts_by_label.get("Sp23", {}).get(series_key)

    Q, dof, pval, Q_over_dof, N_used = cochran_q_three_index_paired(s_su22, s_fa22, s_sp23)
    add_test_text(ax, Q, dof, pval, Q_over_dof, N_used)

def plot_bin_2x3(bin_tag, p_su22, p_fa22, p_sp23, out_dir):
    """
    Produce the 2x3 canvas for one bin_tag.
    """
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(make_title(bin_tag), fontsize=16, y=0.98)

    period_dicts_by_label = {
        "Su22": p_su22,
        "Fa22": p_fa22,
        "Sp23": p_sp23,
    }

    # Row 0
    plot_one_series(
        axes[0, 0],
        "ALUsin",
        ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    plot_one_series(
        axes[0, 1],
        "AULsin",
        ylabel=r"$F_{UL}^{\sin\phi}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    plot_one_series(
        axes[0, 2],
        "AULsin2",
        ylabel=r"$F_{UL}^{\sin2\phi}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    # Row 1
    plot_one_series(
        axes[1, 0],
        "ALLn0",
        ylabel=r"$F_{LL}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    plot_one_series(
        axes[1, 1],
        "ALLcos",
        ylabel=r"$F_{LL}^{\cos\phi}/F_{UU}$",
        period_dicts_by_label=period_dicts_by_label
    )

    # Legend pad
    ax_leg = axes[1, 2]
    ax_leg.axis("off")

    handles = []
    labels = []
    for lab in ["Su22", "Fa22", "Sp23"]:
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
        description="Run-period study: compare Su22/Fa22/Sp23 for ep->en pi+ structure-function ratios (no rad/migration), with per-subplot compatibility p-values (index paired)."
    )
    ap.add_argument("su22", type=str, help="Su22 fit-results .txt")
    ap.add_argument("fa22", type=str, help="Fa22 fit-results .txt")
    ap.add_argument("sp23", type=str, help="Sp23 fit-results .txt")
    args = ap.parse_args()

    su22 = parse_asym_file(args.su22)
    fa22 = parse_asym_file(args.fa22)
    sp23 = parse_asym_file(args.sp23)

    available_bins = detect_available_bins(su22, fa22, sp23)
    if not available_bins:
        print("[ERROR] No recognizable xB-bin sections found (e.g. enpiLowxBGEchi2FitsALUsinphi).")
        sys.exit(2)
    #endif

    os.makedirs(OUT_DIR, exist_ok=True)

    for bin_tag in available_bins:
        p_su22 = build_period_dict(su22, bin_tag)
        p_fa22 = build_period_dict(fa22, bin_tag)
        p_sp23 = build_period_dict(sp23, bin_tag)

        have_any = False
        for k in ["ALUsin", "AULsin", "AULsin2", "ALLn0", "ALLcos"]:
            if (p_su22.get(k) is not None) or (p_fa22.get(k) is not None) or (p_sp23.get(k) is not None):
                have_any = True
                break
            #endif
        #endfor

        if not have_any:
            print(f"[WARN] Skipping bin '{bin_tag}': no plottable series found in any input file.")
            continue
        #endif

        plot_bin_2x3(bin_tag, p_su22, p_fa22, p_sp23, OUT_DIR)
    #endfor

if __name__ == "__main__":
    main()
# endif