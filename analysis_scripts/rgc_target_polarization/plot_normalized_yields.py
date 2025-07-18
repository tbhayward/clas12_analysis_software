#!/usr/bin/env python3
"""
plot_normalized_yields.py

Module to plot high-quality normalized x_B yield histograms for three run periods
(RGC_Su22, RGC_Fa22, RGC_Sp23) and five target types.
Improved styling: fine binning, step histograms, dynamic y-axis per subplot.
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
    For each target in [NH3, C, CH2, He, ET], plot the x_B distribution
    for RGC_Su22, RGC_Fa22, RGC_Sp23 normalized by accumulated charge.
    Arranged in a 2x3 grid; the sixth panel is blank.
    Saves to 'output/normalized_yields.pdf'.

    Args:
        trees (dict): nested dict trees[period][target] with uproot trees
        xB_bins (list): [xmin, ..., xmax] for histogram extent
    """
    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # uniform bin edges for all histograms
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins = np.linspace(xmin, xmax, N_BINS + 1)

    # Create canvas without shared y-axis so each subplot auto-scales
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=False)
    axes = axes.flatten()

    for idx, target in enumerate(targets):
        ax = axes[idx]
        max_count = 0.0
        for period in periods:
            tree = trees[period][target]
            x = tree["x"].array(library="np")
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

        # dynamic y-axis: [0, 1.2*peak]
        ax.set_ylim(0, 1.2 * max_count)
        ax.set_xlabel(r'$x_{B}$')
        ax.set_ylabel('counts / nC')
        ax.set_title(target)
        ax.legend(frameon=False, fontsize='small')

    # hide the unused 6th subplot
    axes[-1].axis('off')

    plt.tight_layout(pad=2.0)
    os.makedirs('output', exist_ok=True)
    out_path = 'output/normalized_yields.pdf'
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[Plot] Saved normalized yields to {out_path}")