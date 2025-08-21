#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Optional dependency check for uproot
try:
    import uproot
except ImportError:
    print("[ERROR] This script requires 'uproot'. Install with:  pip install uproot")
    sys.exit(1)

# ------------------------- Config -------------------------
TREE_NAME = "PhysicsEvents"

# Axis limits (matched to your analysis)
Q2_LIM  = (1.0, 7.0)
X_LIM   = (0.0, 0.70)
TP_LIM  = (0.0, 1.30)   # -t (positive)
MX2_LIM = (0.7, 1.1)

NBINS = 60  # common number of bins for all panels

LABELS = {
    "Q2":  r"$Q^{2}$ (GeV$^{2}$)",
    "x":   r"$x_{B}$",
    "tp":  r"$-t$ (GeV$^{2}$)",
    "Mx2": r"$M_{x}^{2}$ (GeV$^{2}$)",
}

COLORS = {"data": "tab:blue", "mc": "tab:orange"}

# ------------------------- Helpers -------------------------
def load_branches(path):
    """
    Load needed branches from a ROOT file (PhysicsEvents tree).
    Returns a dict with numpy arrays.
    """
    with uproot.open(path) as f:
        if TREE_NAME not in f:
            raise RuntimeError(f"Tree '{TREE_NAME}' not found in {path}")
        tree = f[TREE_NAME]
        # Read as numpy arrays
        arrs = tree.arrays(["Q2", "x", "t", "Mx2"], library="np")
    # build -t as positive quantity
    out = {
        "Q2":  np.asarray(arrs["Q2"], dtype=float),
        "x":   np.asarray(arrs["x"], dtype=float),
        "tp":  -np.asarray(arrs["t"], dtype=float),  # -t
        "Mx2": np.asarray(arrs["Mx2"], dtype=float),
    }
    # keep only finite values
    for k, v in out.items():
        m = np.isfinite(v)
        out[k] = v[m]
    return out

def edges_and_centers(lim, nb):
    edges = np.linspace(lim[0], lim[1], nb + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return edges, centers

def norm_counts(data, edges):
    """Return normalized counts (sum=1) for the given data and bin edges."""
    if data.size == 0:
        return np.zeros(len(edges) - 1, dtype=float)
    counts, _ = np.histogram(data, bins=edges)
    total = counts.sum()
    return counts / total if total > 0 else counts.astype(float)

def plot_step(ax, centers, y, *, label, color):
    ax.step(centers, y, where="mid", color=color, label=label, linewidth=1.8)

# ------------------------- Main -------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Compare basic kinematic distributions (normalized) between data and MC."
    )
    ap.add_argument("data_file", help="Input ROOT file for data")
    ap.add_argument("mc_file",   help="Input ROOT file for MC")
    args = ap.parse_args()

    if not os.path.isfile(args.data_file):
        print(f"[ERROR] Data file not found: {args.data_file}")
        sys.exit(1)
    if not os.path.isfile(args.mc_file):
        print(f"[ERROR] MC file not found: {args.mc_file}")
        sys.exit(1)

    # Load
    data = load_branches(args.data_file)
    mc   = load_branches(args.mc_file)

    # Prepare bins
    q2_edges, q2_cent = edges_and_centers(Q2_LIM,  NBINS)
    x_edges,  x_cent  = edges_and_centers(X_LIM,   NBINS)
    tp_edges, tp_cent = edges_and_centers(TP_LIM,  NBINS)
    m_edges,  m_cent  = edges_and_centers(MX2_LIM, NBINS)

    # Compute normalized counts
    q2_d = norm_counts(data["Q2"],  q2_edges);  q2_m = norm_counts(mc["Q2"],  q2_edges)
    x_d  = norm_counts(data["x"],   x_edges);   x_m  = norm_counts(mc["x"],   x_edges)
    tp_d = norm_counts(data["tp"],  tp_edges);  tp_m = norm_counts(mc["tp"],  tp_edges)
    m_d  = norm_counts(data["Mx2"], m_edges);   m_m  = norm_counts(mc["Mx2"], m_edges)

    # Plot 1x4
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.18, top=0.92, wspace=0.15)
    fig.suptitle(r"$ep \rightarrow en\pi^{+}$ — Data vs MC (normalized distributions)", fontsize=14)

    # Q2
    ax = axes[0]
    plot_step(ax, q2_cent, q2_d, label="data", color=COLORS["data"])
    plot_step(ax, q2_cent, q2_m, label="mc",   color=COLORS["mc"])
    ax.set_xlabel(LABELS["Q2"]); ax.set_xlim(Q2_LIM)

    # xB
    ax = axes[1]
    plot_step(ax, x_cent, x_d, label="data", color=COLORS["data"])
    plot_step(ax, x_cent, x_m, label="mc",   color=COLORS["mc"])
    ax.set_xlabel(LABELS["x"]); ax.set_xlim(X_LIM)

    # -t
    ax = axes[2]
    plot_step(ax, tp_cent, tp_d, label="data", color=COLORS["data"])
    plot_step(ax, tp_cent, tp_m, label="mc",   color=COLORS["mc"])
    ax.set_xlabel(LABELS["tp"]); ax.set_xlim(TP_LIM)

    # Mx2
    ax = axes[3]
    plot_step(ax, m_cent, m_d, label="data", color=COLORS["data"])
    plot_step(ax, m_cent, m_m, label="mc",   color=COLORS["mc"])
    ax.set_xlabel(LABELS["Mx2"]); ax.set_xlim(MX2_LIM)

    # Common y-label, legend, and formatting
    for ax in axes:
        ax.grid(True, linestyle="--", alpha=0.3)
    axes[0].set_ylabel("normalized counts")
    axes[0].legend(loc="best", frameon=True)

    # Save
    outdir = os.path.join("output", "enpi+")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "basic_mc_comparison.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")

if __name__ == "__main__":
    main()