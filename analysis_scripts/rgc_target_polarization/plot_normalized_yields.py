#!/usr/bin/env python3
"""
plot_normalized_yields.py

Module to plot high-quality normalized x_B yield histograms for three run periods
(RGC_Su22, RGC_Fa22, RGC_Sp23) and five target types.
Generates five figures:
  1) Absolute normalization (counts/nC) with dynamic y-limits
  2) Relative to Su22 ratio
  3) Per-run He yields
  4) Per-run ET yields
  5) Per-run CH2 yields

Gracefully handles basket decompression errors, cycles through distinct colors
and linestyles for per-run plots, and prints the integral, mean, and std of each
per-run histogram per target and run period.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Absolute accumulated charge (in nC) per period and target
CHARGE = {
    "RGC_Su22": {"NH3": 3686969.636627, "C": 363715.413199, "CH2": 189723.230200, "He": 266575.708348, "ET": 586020.625415},
    "RGC_Fa22": {"NH3": 5509572.178076, "C": 2006703.362768, "CH2": 1833777.957336, "He": 275693.983222, "ET":  57332.748981},
    "RGC_Sp23": {"NH3": 1620599.347496, "C":  383030.166502, "CH2":  436816.389755, "He": 266597.436226, "ET": 171738.938831},
}

# Styling parameters
PERIOD_COLORS = {"RGC_Su22": "black", "RGC_Fa22": "blue", "RGC_Sp23": "red"}
LINE_WIDTH    = 1.8
N_BINS        = 100   # fine binning for smooth shapes
Y_MAX_RATIO   = 2.0   # fixed y-axis upper limit for ratio plot

def safe_array(tree, branch):
    """
    Safely read a branch from an uproot tree, catching decompression errors.
    Returns an empty array on failure.
    """
    try:
        return tree[branch].array(library="np")
    except Exception as e:
        print(f"[Warning] Could not read '{branch}' from tree: {e}")
        return np.empty(0)

def plot_normalized_yields(trees, xB_bins):
    """
    Generate five different plots of normalized yields and per-run distributions,
    and print integrals, means, and stds per target and run period.
    """
    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # prepare bins and centers
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins    = np.linspace(xmin, xmax, N_BINS + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    # precompute normalized histograms (counts / total charge)
    norm_hist = {p: {} for p in periods}
    for p in periods:
        for t in targets:
            x = safe_array(trees[p][t], "x")
            counts = np.histogram(x, bins=bins)[0] if x.size else np.zeros(len(bins)-1)
            norm_hist[p][t] = counts / CHARGE[p][t]

    # -------- FIGURE 1: Absolute normalization --------
    fig1, axs1 = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    axs1 = axs1.flatten()
    for i, t in enumerate(targets):
        ax = axs1[i]
        peak = 0.0
        for p in periods:
            y = norm_hist[p][t]
            peak = max(peak, y.max() if y.size else 0.0)
            ax.step(
                centers, y,
                where='mid',
                color=PERIOD_COLORS[p],
                linewidth=LINE_WIDTH,
                label=p.replace("RGC_", "")
            )
        ax.set_ylim(0, 1.2 * peak)
        ax.set_title(t)
        ax.set_xlabel(r"$x_{B}$")
        ax.set_ylabel("counts / nC")
        ax.legend(frameon=False, fontsize="small")
    axs1[-1].axis("off")
    plt.tight_layout(pad=2.0)
    os.makedirs("output", exist_ok=True)
    fig1.savefig("output/normalized_yields.pdf")
    plt.close(fig1)
    print("[Plot] Saved 'output/normalized_yields.pdf'")

    # -------- FIGURE 2: Relative to Su22 --------
    fig2, axs2 = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    axs2 = axs2.flatten()
    base = "RGC_Su22"
    for i, t in enumerate(targets):
        ax = axs2[i]
        bvals = norm_hist[base][t]
        for p in periods:
            y = norm_hist[p][t]
            ratio = np.where(bvals > 0, y / bvals, np.nan)
            ax.step(
                centers, ratio,
                where='mid',
                color=PERIOD_COLORS[p],
                linewidth=LINE_WIDTH,
                label=p.replace("RGC_", "")
            )
        ax.set_ylim(0, Y_MAX_RATIO)
        ax.set_title(t)
        ax.set_xlabel(r"$x_{B}$")
        ax.set_ylabel("ratio / Su22")
        ax.legend(frameon=False, fontsize="small")
    axs2[-1].axis("off")
    plt.tight_layout(pad=2.0)
    fig2.savefig("output/normalized_yields_ratio.pdf")
    plt.close(fig2)
    print("[Plot] Saved 'output/normalized_yields_ratio.pdf'")

    # -------- load per-run accumulated charges --------
    runinfo = (
        "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
        "asymmetry_extraction/imports/clas12_run_info.csv"
    )
    run_df = pd.read_csv(
        runinfo, header=None, comment="#",
        usecols=[0,1], names=["run","charge"]
    )
    charge_map = run_df.set_index("run")["charge"].to_dict()

    # -------- helper for per-run plotting --------
    def per_run_plot(target):
        """
        Make a 1×3 grid of per-run normalized yields for the given target,
        cycle through distinct colors and linestyles, and print integrals,
        means, and stds per period.
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)
        for ax, p in zip(axes, periods):
            tree = trees[p][target]
            rn   = safe_array(tree, "runnum")
            xv   = safe_array(tree, "x")
            ur   = np.unique(rn)
            cmap = plt.get_cmap("tab20")
            linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
            integrals = []

            for i, run in enumerate(ur):
                ch = charge_map.get(run)
                if ch is None:
                    print(f"[Warning] missing charge for run {run}, skipping")
                    continue
                mask = (rn == run)
                cnt  = np.histogram(xv[mask], bins=bins)[0]
                norm = cnt / ch
                integ = norm.sum()
                integrals.append(integ)
                color = cmap(i % 20)
                style = linestyles[(i // 20) % len(linestyles)]
                ax.step(
                    centers, norm,
                    where='mid',
                    color=color,
                    linestyle=style,
                    linewidth=1.5,
                    label=str(run)
                )
                print(f"[Integral] Period={p}, Target={target}, Run={run}, Integral={integ:.4f}")

            if integrals:
                mean_int = np.mean(integrals)
                std_int  = np.std(integrals)
                print(f"[Stats]  Period={p}, Target={target}, Mean Integral={mean_int:.4f}, Std Integral={std_int:.4f}")

            ax.set_title(f"{p.replace('RGC_','')} {target}")
            ax.set_xlabel(r"$x_{B}$")
            ax.set_ylabel("counts / nC")
            ax.legend(fontsize="x-small", ncol=2, frameon=False)

        plt.tight_layout(pad=2.0)
        outname = f"output/normalized_yields_runs_{target.lower()}.pdf"
        fig.savefig(outname)
        plt.close(fig)
        print(f"[Plot] Saved '{outname}'\n")

    # -------- FIGURES 3–5: per-run for He, ET, and CH2 --------
    for T in ("He", "ET", "CH2"):
        per_run_plot(T)