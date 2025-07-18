#!/usr/bin/env python3
"""
plot_normalized_yields.py

Module to plot high-quality normalized x_B yield histograms for three run periods
(RGC_Su22, RGC_Fa22, RGC_Sp23) and five target types.
Generates two figures using precomputed normalized histograms:
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
PERIOD_COLORS = {"RGC_Su22": "C0", "RGC_Fa22": "C1", "RGC_Sp23": "C2"}
LINE_WIDTH = 1.8
N_BINS = 100  # fine binning for smooth shapes


def plot_normalized_yields(trees, xB_bins):
    """
    Generate two multi-panel plots using precomputed histograms:
      1) Absolute normalization (counts/nC)
      2) Relative to Su22 (ratio)
    Saves to 'output/normalized_yields.pdf' and 'output/normalized_yields_ratio.pdf'.

    Args:
        trees (dict): nested dict trees[period][target]
        xB_bins (list): bin edges for x_B
    """
    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # Create uniform bins and bin centers
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins = np.linspace(xmin, xmax, N_BINS + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    # Precompute normalized histograms for all period-target combinations
    norm_hist = {period: {} for period in periods}
    for period in periods:
        for target in targets:
            x = trees[period][target]["x"].array(library="np")
            counts, _ = np.histogram(x, bins=bins)
            norm_counts = counts / CHARGE[period][target]
            norm_hist[period][target] = norm_counts

    # --------------------------------------------------
    # FIGURE 1: Absolute normalization (counts / nC)
    # --------------------------------------------------
    fig1, axes1 = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=False)
    axes1 = axes1.flatten()
    for idx, target in enumerate(targets):
        ax = axes1[idx]
        max_c = 0.0
        for period in periods:
            counts = norm_hist[period][target]
            ax.step(centers, counts,
                    where='mid',
                    color=PERIOD_COLORS[period],
                    linewidth=LINE_WIDTH,
                    label=period.replace('RGC_', ''))
            max_c = max(max_c, counts.max() if counts.size>0 else 0)
        ax.set_ylim(0, 1.2 * max_c)
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
    # FIGURE 2: Relative to Su22 (ratio)
    # --------------------------------------------------
    fig2, axes2 = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=False)
    axes2 = axes2.flatten()
    base_period = 'RGC_Su22'
    for idx, target in enumerate(targets):
        ax = axes2[idx]
        base = norm_hist[base_period][target]
        max_r = 0.0
        for period in periods:
            counts = norm_hist[period][target]
            # avoid division by zero
            ratio = np.where(base>0, counts/base, np.nan)
            ax.step(centers, ratio,
                    where='mid',
                    color=PERIOD_COLORS[period],
                    linewidth=LINE_WIDTH,
                    label=period.replace('RGC_', ''))
            max_r = max(max_r, np.nanmax(ratio))
        ax.set_ylim(0, 1.2 * max_r)
        ax.set_xlabel(r'$x_{B}$')
        ax.set_ylabel('ratio / Su22')
        ax.set_title(target)
        ax.legend(frameon=False, fontsize='small')
    axes2[-1].axis('off')
    plt.tight_layout(pad=2.0)
    fig2.savefig('output/normalized_yields_ratio.pdf')
    plt.close(fig2)
    print("[Plot] Saved 'output/normalized_yields_ratio.pdf'")