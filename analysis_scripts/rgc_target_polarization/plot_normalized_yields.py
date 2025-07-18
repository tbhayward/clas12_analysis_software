#!/usr/bin/env python3
"""
plot_normalized_yields.py

Module to plot high-quality normalized x_B yield histograms for three run periods
(RGC_Su22, RGC_Fa22, RGC_Sp23) and five target types.
Generates two figures:
  1) Absolute normalization (counts/nC)
  2) Relative to Su22 yields
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Absolute accumulated charge (in nC) per period and target
CHARGE = {
    "RGC_Su22": {"NH3": 3686969.636627, "C": 363715.413199, "CH2": 189723.230200, "He": 382585.145820, "ET": 470011.187943},
    "RGC_Fa22": {"NH3": 5509572.178076, "C": 2006703.362768, "CH2": 1833777.957336, "He": 378165.719362, "ET":  57332.748981},
    "RGC_Sp23": {"NH3": 1948649.296492, "C":  552805.225420, "CH2":  436816.389755, "He": 279407.815704, "ET": 171738.938831},
}

# Styling parameters
PERIOD_COLORS = {
    "RGC_Su22": "C0",
    "RGC_Fa22": "C1",
    "RGC_Sp23": "C2",
}
LINE_WIDTH = 1.8
N_BINS = 100  # fine binning for smooth shapes


def plot_normalized_yields(trees, xB_bins):
    """
    Generate two multi-panel plots:
      - normalized yields (counts/nC)
      - relative yields (each period / Su22)
    Saves to 'output/normalized_yields.pdf' and 'output/normalized_yields_ratio.pdf'.

    Args:
        trees (dict): nested dict trees[period][target] with uproot trees
        xB_bins (list): [xmin, ..., xmax] for histogram extent
    """
    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # define uniform bins over xB
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins = np.linspace(xmin, xmax, N_BINS + 1)

    # --------------------------------------------------
    # FIGURE 1: Absolute normalization (counts / nC)
    # --------------------------------------------------
    fig1, axes1 = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=False)
    axes1 = axes1.flatten()

    for idx, target in enumerate(targets):
        ax = axes1[idx]
        max_count = 0.0
        for period in periods:
            x = trees[period][target]["x"].array(library="np")
            counts, edges = np.histogram(x, bins=bins)
            centers = 0.5 * (edges[:-1] + edges[1:])
            norm_counts = counts / CHARGE[period][target]
            ax.step(
                centers, norm_counts,
                where='mid',
                color=PERIOD_COLORS[period],
                linewidth=LINE_WIDTH,
                label=period.replace('RGC_', '')
            )
            if norm_counts.size > 0:
                max_count = max(max_count, norm_counts.max())
        ax.set_ylim(0, 1.2 * max_count)
        ax.set_xlabel(r'$x_{B}$')
        ax.set_ylabel('counts / nC')
        ax.set_title(target)
        ax.legend(frameon=False, fontsize='small')
    axes1[-1].axis('off')
    plt.tight_layout(pad=2.0)
    os.makedirs('output', exist_ok=True)
    fig1.savefig('output/normalized_yields.pdf')
    plt.close(fig1)
    print("[Plot] Saved 'output/normalized_yields.pdf'")

    # --------------------------------------------------
    # FIGURE 2: Relative to Su22 (ratio of normalized yields)
    # --------------------------------------------------
    fig2, axes2 = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=False)
    axes2 = axes2.flatten()

    for idx, target in enumerate(targets):
        ax = axes2[idx]
        # compute base (Su22) normalized counts
        x_su = trees['RGC_Su22'][target]['x'].array(library='np')
        counts_su, _ = np.histogram(x_su, bins=bins)
        norm_su = counts_su / CHARGE['RGC_Su22'][target]
        centers = 0.5 * (bins[:-1] + bins[1:])
        max_ratio = 0.0
        for period in periods:
            x = trees[period][target]['x'].array(library='np')
            counts, _ = np.histogram(x, bins=bins)
            norm = counts / CHARGE[period][target]
            # avoid zero-division
            ratio = norm / np.where(norm_su>0, norm_su, np.nan)
            ax.step(
                centers, ratio,
                where='mid',
                color=PERIOD_COLORS[period],
                linewidth=LINE_WIDTH,
                label=period.replace('RGC_', '')
            )
            # track peak ratio (exclude NaN)
            if np.isfinite(ratio).any():
                max_ratio = max(max_ratio, np.nanmax(ratio))
        ax.set_ylim(0, 1.2 * max_ratio)
        ax.set_xlabel(r'$x_{B}$')
        ax.set_ylabel('Relative to Su22')
        ax.set_title(target)
        ax.legend(frameon=False, fontsize='small')
    axes2[-1].axis('off')
    plt.tight_layout(pad=2.0)
    fig2.savefig('output/normalized_yields_ratio.pdf')
    plt.close(fig2)
    print("[Plot] Saved 'output/normalized_yields_ratio.pdf'")