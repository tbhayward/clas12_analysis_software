#!/usr/bin/env python3
"""
plot_normalized_yields.py

Module to plot high-quality normalized x_B yield histograms for three run periods
(RGC_Su22, RGC_Fa22, RGC_Sp23) and five target types.
Generates:
  1) Absolute normalization (counts/nC)
  2) Relative to Su22 ratio
  3–7) Per-run yields for each target (NH3, C, CH2, He, ET), one 1×3 canvas each

Handles basket decompression errors, parallel per-period work, distinct colors/
linestyles, and writes every run’s integral + per-period mean/σ to both console
and a single log file.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

# Absolute accumulated charge (in nC) per period and target
CHARGE = {
    "RGC_Su22": {"NH3": 3686969.636627, "C": 363715.413199, "CH2": 189723.230200, "He": 266575.708348, "ET": 586020.625415},
    "RGC_Fa22": {"NH3": 5509572.178076, "C": 2006703.362768, "CH2": 1833777.957336, "He": 275693.983222, "ET":  57332.748981},
    "RGC_Sp23": {"NH3": 1620599.347496, "C":  383030.166502, "CH2":  436816.389755, "He": 266597.436226, "ET": 171738.938831},
}

# Styling parameters
PERIOD_COLORS = {"RGC_Su22": "black", "RGC_Fa22": "blue", "RGC_Sp23": "red"}
LINE_WIDTH    = 1.8
N_BINS        = 100   # fine binning
Y_MAX_RATIO   = 2.0   # fixed y-axis upper limit for ratio plots

# Location of per-run charge info
RUNINFO = (
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "asymmetry_extraction/imports/clas12_run_info.csv"
)

def safe_array(tree, branch):
    """
    Safely read a branch from an uproot tree, catching decompression errors.
    Returns an empty numpy array on failure.
    """
    try:
        return tree[branch].array(library="np")
    except Exception as e:
        print(f"[Warning] could not read '{branch}' from tree: {e}")
        return np.empty(0)

def process_period_runs(period, target, bins, centers, charge_map, trees):
    """
    Worker for per-run plotting: for (period, target), compute each run’s
    normalized histogram, integral, then compute that period’s mean & std.
    """
    tree = trees[period][target]
    rn   = safe_array(tree, "runnum")
    xv   = safe_array(tree, "x")
    runs = np.unique(rn)
    entries = []
    for run in runs:
        ch = charge_map.get(run)
        if ch is None:
            continue
        mask = (rn == run)
        cnt  = np.histogram(xv[mask], bins=bins)[0]
        norm = cnt / ch
        integ = norm.sum()
        entries.append({"run": run, "norm": norm, "integral": integ})
    ints = [e["integral"] for e in entries]
    mean, std = (float(np.mean(ints)), float(np.std(ints))) if ints else (0.0, 0.0)
    return {"period": period, "target": target,
            "centers": centers, "entries": entries,
            "mean": mean, "std": std}

def plot_normalized_yields(trees, xB_bins):
    # create output folder & open log file
    os.makedirs("output", exist_ok=True)
    log_path = "output/per_run_integrals.txt"
    log = open(log_path, "w", buffering=1)
    log.write("Per-run integrals & stats\n")
    log.write("==========================\n\n")

    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # histogram bins and centers
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins    = np.linspace(xmin, xmax, N_BINS+1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    # load per-run charges
    run_df = pd.read_csv(RUNINFO, header=None, comment="#",
                         usecols=[0,1], names=["run","charge"])
    charge_map = run_df.set_index("run")["charge"].to_dict()

    # precompute normalized totals
    norm_hist = {p:{} for p in periods}
    for p in periods:
        for t in targets:
            x = safe_array(trees[p][t], "x")
            cnt = np.histogram(x, bins=bins)[0] if x.size else np.zeros(len(bins)-1)
            norm_hist[p][t] = cnt / CHARGE[p][t]

    # FIG 1: absolute
    fig1, axs1 = plt.subplots(2,3,figsize=(15,8), sharex=True)
    axs1 = axs1.flatten()
    for i, t in enumerate(targets):
        ax = axs1[i]
        peak = 0
        for p in periods:
            y = norm_hist[p][t]
            peak = max(peak, y.max() if y.size else 0)
            ax.step(centers, y, where='mid',
                    color=PERIOD_COLORS[p], linewidth=LINE_WIDTH,
                    label=p.replace("RGC_",""))
        ax.set_ylim(0, 1.2*peak)
        ax.set_title(t)
        ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("counts / nC")
        ax.legend(frameon=False, fontsize="small")
    axs1[-1].axis("off")
    plt.tight_layout(pad=2)
    fig1.savefig("output/normalized_yields.pdf"); plt.close(fig1)
    print("[Plot] saved output/normalized_yields.pdf")

    # FIG 2: ratio to Su22
    fig2, axs2 = plt.subplots(2,3,figsize=(15,8), sharex=True)
    axs2 = axs2.flatten()
    base = "RGC_Su22"
    for i, t in enumerate(targets):
        ax = axs2[i]
        bvals = norm_hist[base][t]
        for p in periods:
            y = norm_hist[p][t]
            r = np.where(bvals>0, y/bvals, np.nan)
            ax.step(centers, r, where='mid',
                    color=PERIOD_COLORS[p], linewidth=LINE_WIDTH,
                    label=p.replace("RGC_",""))
        ax.set_ylim(0, Y_MAX_RATIO)
        ax.set_title(t)
        ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("ratio / Su22")
        ax.legend(frameon=False, fontsize="small")
    axs2[-1].axis("off")
    plt.tight_layout(pad=2)
    fig2.savefig("output/normalized_yields_ratio.pdf"); plt.close(fig2)
    print("[Plot] saved output/normalized_yields_ratio.pdf")

    # FIGs 3–7: per-run for each target
    for target in targets:
        # build args for parallel per-period processing
        args = [(p, target, bins, centers, charge_map, trees) for p in periods]
        with ProcessPoolExecutor(max_workers=3) as exe:
            futures = {exe.submit(process_period_runs, *a): a[0] for a in args}
            results = {}
            for fut in as_completed(futures):
                period = futures[fut]
                res = fut.result()
                results[period] = res

        # plot 1x3
        fig, axes = plt.subplots(1,3,figsize=(18,5), sharex=True, sharey=True)
        for ax, p in zip(axes, periods):
            res = results[p]
            # log header
            hdr = f"=== Period {p}, Target {target} ===\n"
            print(hdr.strip()); log.write(hdr)
            for idx, e in enumerate(res["entries"]):
                run   = e["run"]
                norm  = e["norm"]
                integ = e["integral"]
                print(f" Run={run}, Integral={integ:.4f}")
                log.write(f" Run={run}, Integral={integ:.4f}\n")
                cmap  = plt.get_cmap("tab20")
                style = ['solid','dashed','dashdot','dotted'][(idx//20)%4]
                color = cmap(idx%20)
                ax.step(centers, norm, where='mid',
                        color=color, linestyle=style, linewidth=1.5,
                        label=str(run))
            stats = f" Mean={res['mean']:.4f}, Std={res['std']:.4f}\n\n"
            print(stats.strip()); log.write(stats)
            ax.set_title(f"{p.replace('RGC_','')} {target}")
            ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("counts / nC")
            ax.legend(fontsize="x-small", ncol=2, frameon=False)

        plt.tight_layout(pad=2)
        out = f"output/normalized_yields_runs_{target.lower()}.pdf"
        fig.savefig(out); plt.close(fig)
        print(f"[Plot] saved {out}\n")

    log.close()
    print(f"[Log] per-run integrals & stats written to {log_path}")