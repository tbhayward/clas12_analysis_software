#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Misidentification vs -t for four xB slices (MC) — 2x4 canvas.

Top row (electron candidate):
  • π− mis-ID  : matching_e_pid == -211
  • K− mis-ID  : matching_e_pid == -321

Bottom row (pion candidate):
  • e+ mis-ID  : matching_p1_pid == -11
  • K+ mis-ID  : matching_p1_pid ==  321
  • p  mis-ID  : matching_p1_pid == 2212

Common cuts per event:
  0.81 < Mx2 < 1.00
  x in slice range
  -t in one of the six bins below
  (No fiducial-status cuts.)

Y-axis is a FRACTION with limits [0, 0.01].

Usage:
  python misid_vs_t_2x4.py <input.root>

Output:
  output/enpi+/misidentification.pdf
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt

try:
    import uproot
except Exception:
    print("[ERROR] Please install uproot (e.g. pip install uproot)")
    sys.exit(1)

# ─────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────
TREE_NAME = "PhysicsEvents"

# xB slices
X_SLICES = {
    "Low":     (0.10, 0.25),
    "MidLow":  (0.25, 0.35),
    "MidHigh": (0.35, 0.45),
    "High":    (0.45, 0.60),
}
SLICE_ORDER = ["Low", "MidLow", "MidHigh", "High"]

# -t bin edges in GeV^2 (ascending) → bins [edge[i], edge[i+1])
T_EDGES_POS = np.array([0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25], dtype=float)
T_CENTERS   = 0.5 * (T_EDGES_POS[:-1] + T_EDGES_POS[1:])
NTBINS      = len(T_EDGES_POS) - 1

# Plot / output
OUT_DIR  = os.path.join("output", "enpi+")
OUT_PATH = os.path.join(OUT_DIR, "misidentification.pdf")
YLIM     = (0.0, 0.03)  # fraction scale

# Styles
STYLE_E_PI  = dict(fmt="o",  color="tab:blue",   label=r"$e^- \leftarrow \pi^{-}$")
STYLE_E_K   = dict(fmt="^",  color="tab:orange", label=r"$e^- \leftarrow K^{-}$")
STYLE_P_E   = dict(fmt="D",  color="tab:purple", label=r"$\pi^- \leftarrow e^{+}$")
STYLE_P_K   = dict(fmt="s",  color="tab:green",  label=r"$\pi^- \leftarrow K^{+}$")
STYLE_P_P   = dict(fmt="v",  color="tab:red",    label=r"$\pi^- \leftarrow p$")

CAPSIZE = 3
MSIZE   = 5.5  # markersize for errorbar

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def slice_index_from_x(x):
    """Return array of indices in {0..3} for xB slice, else -1."""
    idx = np.full(x.shape, -1, dtype=np.int16)
    for k, name in enumerate(SLICE_ORDER):
        xa, xb = X_SLICES[name]
        m = (x > xa) & (x < xb)
        idx[m] = k
    return idx

def tbin_index_from_tpos(tpos):
    """Return -t bin index in {0..NTBINS-1}, else -1 if outside range."""
    bix = np.digitize(tpos, T_EDGES_POS, right=False) - 1
    bad = (bix < 0) | (bix >= NTBINS)
    bix[bad] = -1
    return bix.astype(np.int16)

def binomial_err(n, N):
    """Standard error for fraction n/N: sqrt(f*(1-f)/N). Returns np.nan if N==0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(N > 0, n / N, np.nan)
        var = np.where(N > 0, f * (1.0 - f) / N, np.nan)
    return f, np.sqrt(var)

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 2:
        print("Usage: python misid_vs_t_2x4.py <input.root>")
        sys.exit(1)

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"[ERROR] File not found: {infile}")
        sys.exit(1)

    branches = ["x", "t", "Mx2", "matching_e_pid", "matching_p1_pid"]

    # Totals per slice/bin (shared denominator for all categories)
    totals = np.zeros((4, NTBINS), dtype=np.int64)

    # Top row (electron candidate)
    n_e_from_pi = np.zeros((4, NTBINS), dtype=np.int64)   # matching_e_pid == -211
    n_e_from_k  = np.zeros((4, NTBINS), dtype=np.int64)   # matching_e_pid == -321

    # Bottom row (pion candidate)
    n_p_from_e  = np.zeros((4, NTBINS), dtype=np.int64)   # matching_p1_pid == -11
    n_p_from_k  = np.zeros((4, NTBINS), dtype=np.int64)   # matching_p1_pid == 321
    n_p_from_p  = np.zeros((4, NTBINS), dtype=np.int64)   # matching_p1_pid == 2212

    # Iterate in chunks
    step = "200 MB"
    try:
        itr = uproot.iterate({infile: TREE_NAME}, branches, step_size=step, library="np")
    except Exception as e:
        print(f"[ERROR] Failed to open '{infile}' or tree '{TREE_NAME}': {e}")
        sys.exit(1)

    for arr in itr:
        # Load
        x   = arr["x"]
        t   = arr["t"]
        mx2 = arr["Mx2"]
        ep  = arr["matching_e_pid"]
        p1  = arr["matching_p1_pid"]

        # Base cut
        base = np.isfinite(x) & np.isfinite(t) & np.isfinite(mx2) & (mx2 > 0.81) & (mx2 < 1.00)
        if not np.any(base):
            continue

        tpos = -t[base]
        xs   =  x[base]
        ep_b =  ep[base]
        p1_b =  p1[base]

        sidx = slice_index_from_x(xs)         # 0..3 or -1
        tbix = tbin_index_from_tpos(tpos)     # 0..5 or -1
        valid = (sidx >= 0) & (tbix >= 0)
        if not np.any(valid):
            continue

        si = sidx[valid]
        bi = tbix[valid]
        epv = ep_b[valid]
        p1v = p1_b[valid]

        # Denominator
        np.add.at(totals, (si, bi), 1)

        # Numerators (exact PIDs for requested channels)
        np.add.at(n_e_from_pi, (si[epv == -211], bi[epv == -211]), 1)
        np.add.at(n_e_from_k,  (si[epv == -321], bi[epv == -321]), 1)

        np.add.at(n_p_from_e,  (si[p1v ==  -11], bi[p1v ==  -11]), 1)
        np.add.at(n_p_from_k,  (si[p1v ==  321], bi[p1v ==  321]), 1)
        np.add.at(n_p_from_p,  (si[p1v == 2212], bi[p1v == 2212]), 1)

    # Fractions + binomial errors
    f_e_pi, se_e_pi = binomial_err(n_e_from_pi, totals)
    f_e_k,  se_e_k  = binomial_err(n_e_from_k,  totals)

    f_p_e,  se_p_e  = binomial_err(n_p_from_e,  totals)
    f_p_k,  se_p_k  = binomial_err(n_p_from_k,  totals)
    f_p_p,  se_p_p  = binomial_err(n_p_from_p,  totals)

    # ── Plotting (2x4) ───────────────────────────────────────────────
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, axes = plt.subplots(2, 4, figsize=(16.8, 7.2), sharey=True, sharex=True)

    fig.suptitle(
        r"Misidentification vs $-t$"
        "\n" r"$0.81<M_{x}^{2}<1.00$",
        fontsize=14, y=0.98
    )

    for col, sname in enumerate(SLICE_ORDER):
        xa, xb = X_SLICES[sname]

        # Top row: electron candidate mis-IDs
        ax_top = axes[0, col]
        # π−
        m = np.isfinite(f_e_pi[col, :])
        if np.any(m):
            ax_top.errorbar(T_CENTERS[m], f_e_pi[col, m], yerr=se_e_pi[col, m],
                            markersize=MSIZE, capsize=CAPSIZE, **STYLE_E_PI)
        # K−
        m = np.isfinite(f_e_k[col, :])
        if np.any(m):
            ax_top.errorbar(T_CENTERS[m], f_e_k[col, m], yerr=se_e_k[col, m],
                            markersize=MSIZE, capsize=CAPSIZE, **STYLE_E_K)

        ax_top.set_ylim(*YLIM)
        ax_top.grid(True, linestyle="--", alpha=0.35)
        ax_top.set_title(rf"${xa:.2f} < x_{{B}} < {xb:.2f}$", fontsize=12)
        if col == 0:
            ax_top.set_ylabel("misID fraction")

        # Put legend only on first column (top row)
        if col == 0:
            ax_top.legend(loc="upper left", frameon=True, edgecolor="black", fontsize=10)

        # Bottom row: pion candidate mis-IDs
        ax_bot = axes[1, col]
        # e+
        m = np.isfinite(f_p_e[col, :])
        if np.any(m):
            ax_bot.errorbar(T_CENTERS[m], f_p_e[col, m], yerr=se_p_e[col, m],
                            markersize=MSIZE, capsize=CAPSIZE, **STYLE_P_E)
        # K+
        m = np.isfinite(f_p_k[col, :])
        if np.any(m):
            ax_bot.errorbar(T_CENTERS[m], f_p_k[col, m], yerr=se_p_k[col, m],
                            markersize=MSIZE, capsize=CAPSIZE, **STYLE_P_K)
        # p
        m = np.isfinite(f_p_p[col, :])
        if np.any(m):
            ax_bot.errorbar(T_CENTERS[m], f_p_p[col, m], yerr=se_p_p[col, m],
                            markersize=MSIZE, capsize=CAPSIZE, **STYLE_P_P)

        ax_bot.set_ylim(*YLIM)
        ax_bot.grid(True, linestyle="--", alpha=0.35)
        ax_bot.set_xlabel(r"$-t\ (\mathrm{GeV}^{2})$")
        if col == 0:
            ax_bot.set_ylabel("misID fraction")
        # Put legend only on first column (bottom row)
        if col == 0:
            ax_bot.legend(loc="upper left", frameon=True, edgecolor="black", fontsize=10)

    # Common x-limits
    for ax in axes.ravel():
        ax.set_xlim(T_EDGES_POS[0], T_EDGES_POS[-1])

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_PATH, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")

if __name__ == "__main__":
    main()