#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the misidentification rate vs -t for four xB slices and plot as 1x4.

Definition (per xB slice, per -t bin):
  percentage = 100 * N(matching_e_pid==-211 AND matching_p1_pid==321) / N(total)

Cuts used in this study:
  0.81 < Mx2 < 1.00
  x in the slice range (see X_SLICES below)
  -t binned using edges: [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]
  (No fiducial-status cuts.)

Usage:
  python misid_vs_t.py <input.root>

Saves:
  output/enpi+/misidentification.pdf
"""

import sys
import os
import numpy as np

try:
    import uproot
except Exception:
    print("[ERROR] Please install uproot (e.g. pip install uproot)")
    sys.exit(1)

import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────
TREE_NAME = "PhysicsEvents"

# xB slices (titles show these ranges)
X_SLICES = {
    "Low":     (0.10, 0.25),
    "MidLow":  (0.25, 0.35),
    "MidHigh": (0.35, 0.45),
    "High":    (0.45, 0.60),
}
SLICE_ORDER = ["Low", "MidLow", "MidHigh", "High"]

# -t bin edges (ascending, in GeV^2); bins are [edge[i], edge[i+1])
T_EDGES_POS = np.array([0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25], dtype=float)
T_CENTERS = 0.5 * (T_EDGES_POS[:-1] + T_EDGES_POS[1:])
NTBINS = len(T_EDGES_POS) - 1

OUT_DIR = os.path.join("output", "enpi+")
OUT_PATH = os.path.join(OUT_DIR, "misidentification.pdf")

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def slice_index_from_x(x):
    """
    Return array of indices in {0..3} mapping each x to its xB slice,
    or -1 if x is outside all slice ranges (open intervals).
    """
    idx = np.full(x.shape, -1, dtype=np.int16)
    for k, name in enumerate(SLICE_ORDER):
        xa, xb = X_SLICES[name]
        m = (x > xa) & (x < xb)
        idx[m] = k
    return idx

def tbin_index_from_tpos(tpos):
    """Return -t bin index in {0..NTBINS-1}, or -1 if outside range."""
    bix = np.digitize(tpos, T_EDGES_POS, right=False) - 1  # [lo,hi)
    bad = (bix < 0) | (bix >= NTBINS)
    bix[bad] = -1
    return bix.astype(np.int16)

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 2:
        print("Usage: python misid_vs_t.py <input.root>")
        sys.exit(1)

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"[ERROR] File not found: {infile}")
        sys.exit(1)

    # Only the branches we need (no fiducial_status here)
    branches = [
        "x", "t", "Mx2",
        "matching_e_pid", "matching_p1_pid",
    ]

    totals = np.zeros((4, NTBINS), dtype=np.int64)
    bads   = np.zeros((4, NTBINS), dtype=np.int64)

    step = "200 MB"
    try:
        itr = uproot.iterate({infile: TREE_NAME}, branches, step_size=step, library="np")
    except Exception as e:
        print(f"[ERROR] Failed to open '{infile}' or tree '{TREE_NAME}': {e}")
        sys.exit(1)

    for arrays in itr:
        x   = arrays["x"]
        t   = arrays["t"]
        mx2 = arrays["Mx2"]
        e_pid  = arrays["matching_e_pid"]
        p1_pid = arrays["matching_p1_pid"]

        # Base cuts for this study (no fid cuts)
        base = (mx2 > 0.81) & (mx2 < 1.00)
        if not np.any(base):
            continue

        tpos = -t[base]   # bin in -t
        xs   = x[base]
        eok  = (e_pid[base] == -211)
        pok  = (p1_pid[base] == 321)

        sidx = slice_index_from_x(xs)          # 0..3 or -1
        tbix = tbin_index_from_tpos(tpos)      # 0..5 or -1
        valid = (sidx >= 0) & (tbix >= 0)
        if not np.any(valid):
            continue

        sidx = sidx[valid]
        tbix = tbix[valid]

        # Totals
        np.add.at(totals, (sidx, tbix), 1)

        # Mis-ID events
        mis = eok[valid] & pok[valid]
        if np.any(mis):
            np.add.at(bads, (sidx[mis], tbix[mis]), 1)

    # Percentages
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(totals > 0, bads / totals, np.nan) * 100.0

    # ── Plotting ──────────────────────────────────────────────────────
    os.makedirs(OUT_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.6), sharey=True)
    fig.suptitle(
        r"Misidentification rate vs $-t$"
        "\n" r"$matching\_e\_pid=-211,\ matching\_{p1}\_pid=321$; "
        r"$0.81<M_{x}^{2}<1.00$",
        fontsize=13, y=0.98
    )

    for k, name in enumerate(SLICE_ORDER):
        ax = axes[k]
        y = frac[k, :]
        m = np.isfinite(y)
        ax.scatter(T_CENTERS[m], y[m], s=28, marker="o")
        ax.set_xlim(T_EDGES_POS[0], T_EDGES_POS[-1])
        ax.set_xlabel(r"$-t\ (\mathrm{GeV}^{2})$")
        if k == 0:
            ax.set_ylabel("percentage of events")
        ax.grid(True, linestyle="--", alpha=0.35)
        xa, xb = X_SLICES[name]
        ax.set_title(rf"${xa:.2f} < x_{{B}} < {xb:.2f}$", fontsize=11)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_PATH, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")

if __name__ == "__main__":
    main()