#!/usr/bin/env python3
"""
plot_normalized_yields.py

Module to plot high-quality normalized x_B yield histograms for three run periods
(RGC_Su22, RGC_Fa22, RGC_Sp23) and five target types.
Generates three figures:
  1) Absolute normalization (counts/nC) with dynamic y-limits
  2) Relative to Su22 yields with fixed y-axis [0,2]
  3) Per-run He normalized yields for each period (1×3 grid)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Absolute accumulated charge (in nC) per period and target
CHARGE = {
    "RGC_Su22": {"NH3": 3686969.636627, "C": 363715.413199, "CH2": 189723.230200, "He": 382585.145820, "ET": 470011.187943},
    "RGC_Fa22": {"NH3": 5509572.178076, "C": 2006703.362768, "CH2": 1833777.957336, "He": 378165.719362, "ET":  57332.748981},
    "RGC_Sp23": {"NH3": 1620599.347496, "C":  383030.166502, "CH2":  436816.389755, "He": 279407.815704, "ET": 171738.938831},
}

# Styling parameters
PERIOD_COLORS = {"RGC_Su22": "C0", "RGC_Fa22": "C1", "RGC_Sp23": "C2"}
LINE_WIDTH = 1.8
N_BINS = 100           # fine binning for smooth shapes
Y_MAX_RATIO = 2.0      # fixed y-axis upper limit for ratio plot


def plot_normalized_yields(trees, xB_bins):
    """
    Generate three plots:
      1) Absolute normalization (counts/nC)
      2) Relative to Su22 (ratio)
      3) Per-run He yields for each period

    Args:
        trees (dict): trees[period][target] uproot trees
        xB_bins (list): edges for x_B range
    """
    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # prepare bins and centers
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins = np.linspace(xmin, xmax, N_BINS + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    # precompute normalized histograms (counts / total charge)
    norm_hist = {p: {} for p in periods}
    for p in periods:
        for t in targets:
            x = trees[p][t]["x"].array(library="np")
            counts, _ = np.histogram(x, bins=bins)
            norm_hist[p][t] = counts / CHARGE[p][t]

    # --------------------------------------------------
    # FIGURE 1: Absolute normalization
    # --------------------------------------------------
    fig1, axes1 = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    axes1 = axes1.flatten()
    for idx, t in enumerate(targets):
        ax = axes1[idx]
        peak = 0.0
        for p in periods:
            y = norm_hist[p][t]
            peak = max(peak, y.max() if y.size > 0 else 0)
            ax.step(centers, y,
                    where='mid', color=PERIOD_COLORS[p], linewidth=LINE_WIDTH,
                    label=p.replace('RGC_', ''))
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
    # FIGURE 2: Relative to Su22
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
            ax.step(centers, ratio,
                    where='mid', color=PERIOD_COLORS[p], linewidth=LINE_WIDTH,
                    label=p.replace('RGC_', ''))
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

    # --------------------------------------------------
    # FIGURE 3: Per-run He normalized yields
    # --------------------------------------------------
    # load run charges
    runinfo = '/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv'
    run_df = pd.read_csv(runinfo, header=None, comment='#', usecols=[0,1], names=['run','charge'])
    charge_map = run_df.set_index('run')['charge'].to_dict()

    fig3, axes3 = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)
    for ax, p in zip(axes3, periods):
        tree = trees[p]['He']
        runnums = tree['runnum'].array(library='np')
        xvals   = tree['x'].array(library='np')
        unique_runs = np.unique(runnums)
        # colormap for runs
        cmap = plt.get_cmap('tab20', len(unique_runs))
        for i, run in enumerate(unique_runs):
            charge = charge_map.get(run)
            if charge is None:
                print(f"Warning: missing charge for run {run}, skipping")
                continue
            mask = (runnums == run)
            counts, _ = np.histogram(xvals[mask], bins=bins)
            norm_counts = counts / charge
            ax.step(centers, norm_counts,
                    where='mid', color=cmap(i), linewidth=1.2,
                    label=str(run))
        ax.set_xlabel(r'$x_{B}$')
        ax.set_ylabel('counts / nC')
        ax.set_title(p.replace('RGC_', ''))
        ax.legend(fontsize='x-small', ncol=2, frameon=False)
    plt.tight_layout(pad=2.0)
    fig3.savefig('output/normalized_yields_runs_he.pdf')
    plt.close(fig3)
    print("[Plot] Saved 'output/normalized_yields_runs_he.pdf'")