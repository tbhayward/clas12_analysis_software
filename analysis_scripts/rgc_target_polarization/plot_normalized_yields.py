#!/usr/bin/env python3
"""
plot_normalized_yields.py

Module to plot high-quality normalized x_B yield histograms for three run periods
(RGC_Su22, RGC_Fa22, RGC_Sp23) and five target types.
Generates two figures:
  1) Absolute normalization (counts/nC) with dynamic y-limits
  2) Relative to Su22 yields with fixed y-axis [0,2]
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Absolute accumulated charge (in nC) per period and target
CHARGE = {
    "RGC_Su22": {
        "NH3": 3686969.636627,
        "C":   363715.413199,
        "CH2": 189723.230200,
        "He":  382585.145820,
        "ET":  470011.187943,
    },
    "RGC_Fa22": {
        "NH3": 5509572.178076,
        "C":   2006703.362768,
        "CH2": 1833777.957336,
        "He":  378165.719362,
        "ET":   57332.748981,
    },
    "RGC_Sp23": {
        # Updated Sp23 totals:
        "NH3": 1620599.347496,
        "C":    383030.166502,
        "CH2":  436816.389755,
        "He":   279407.815704,
        "ET":   171738.938831,
    },
}

# Styling parameters
PERIOD_COLORS = {
    "RGC_Su22": "C0",
    "RGC_Fa22": "C1",
    "RGC_Sp23": "C2",
}
LINE_WIDTH = 1.8
N_BINS = 100  # fine binning for smooth shapes
Y_MAX_RATIO = 2.0  # fixed y-axis upper limit for ratio plot


def plot_normalized_yields(trees, xB_bins):
    """
    Generate two multi-panel plots using precomputed histograms:
      1) Absolute normalization (counts/nC) with dynamic autoscaling
      2) Relative to Su22 (ratio) with y-axis [0, Y_MAX_RATIO]
    """
    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # Uniform bins and centers
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins = np.linspace(xmin, xmax, N_BINS + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    # Precompute normalized histograms
    norm_hist = {p: {} for p in periods}
    for p in periods:
        for t in targets:
            x = trees[p][t]["x"].array(library="np")
            counts, _ = np.histogram(x, bins=bins)
            norm_hist[p][t] = counts / CHARGE[p][t]

    # --------------------------------------------------
    # FIGURE 1: Absolute normalization (counts / nC)
    # --------------------------------------------------
    fig1, axes1 = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    axes1 = axes1.flatten()
    for idx, t in enumerate(targets):
        ax = axes1[idx]
        peak = 0.0
        for p in periods:
            y = norm_hist[p][t]
            peak = max(peak, y.max() if y.size > 0 else 0)
            ax.step(
                centers, y,
                where='mid',
                color=PERIOD_COLORS[p],
                linewidth=LINE_WIDTH,
                label=p.replace('RGC_', '')
            )
        ax.set_ylim(0, 1.2 * peak)
        ax.set_xlabel(r'$x_{B}$')
        ax.set_ylabel('counts / nC')
        ax.set_title(t)
        ax.legend(frameon=False, fontsize='small')
    axes1[-1].axis('off')
    plt.tight_layout(pad=2.0)
    os.makedirs('output', exist_ok=True)
    fig1.savefig('output/normalized_yields.pdf')
    plt.close(fig1)
    print("[Plot] Saved 'output/normalized_yields.pdf'")

    # --------------------------------------------------
    # FIGURE 2: Relative to Su22 (ratio)
    # --------------------------------------------------
    fig2, axes2 = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    axes2 = axes2.flatten()
    base = 'RGC_Su22'
    for idx, t in enumerate(targets):
        ax = axes2[idx]
        base_vals = norm_hist[base][t]
        for p in periods:
            vals = norm_hist[p][t]
            ratio = np.where(base_vals > 0, vals / base_vals, np.nan)
            ax.step(
                centers, ratio,
                where='mid',
                color=PERIOD_COLORS[p],
                linewidth=LINE_WIDTH,
                label=p.replace('RGC_', '')
            )
        ax.set_ylim(0, Y_MAX_RATIO)
        ax.set_xlabel(r'$x_{B}$')
        ax.set_ylabel('ratio / Su22')
        ax.set_title(t)
        ax.legend(frameon=False, fontsize='small')
    axes2[-1].axis('off')
    plt.tight_layout(pad=2.0)
    fig2.savefig('output/normalized_yields_ratio.pdf')
    plt.close(fig2)
    print("[Plot] Saved 'output/normalized_yields_ratio.pdf'")