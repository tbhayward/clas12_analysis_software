#!/usr/bin/env python3
"""
plot_normalized_yields.py

Module to plot high-quality normalized x_B yield histograms for three run periods
(RGC_Su22, RGC_Fa22, RGC_Sp23) and five target types.
Generates five figures:
  1) Absolute normalization (counts/nC) with dynamic y-limits
  2) Relative to Su22 ratio
  3) Per-run NH3 yields
  4) Per-run C yields
  5) Per-run CH2 yields
  6) Per-run He yields
  7) Per-run ET yields

Gracefully handles basket decompression errors, cycles through distinct colors
and linestyles for per-run plots, parallelizes per-period work, and writes a
summary text file with every run’s integral, plus per-period mean/std and
2σ flags.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

# Absolute accumulated charge (in nC) per period and target
CHARGE = {
    "RGC_Su22": {"NH3": 3686969.636627, "C": 363715.413199, "CH2": 189723.230200, "He": 121362.943484, "ET": 586020.625415},
    "RGC_Fa22": {"NH3": 5509572.178076, "C": 2006703.362768, "CH2": 1756673.682648, "He": 275693.983222, "ET":  57332.748981},
    "RGC_Sp23": {"NH3": 1620599.347496, "C":  383030.166502, "CH2":  436816.389755, "He": 266597.436226, "ET": 171738.938831},
}

# Styling parameters
PERIOD_COLORS = {"RGC_Su22": "black", "RGC_Fa22": "blue", "RGC_Sp23": "red"}
LINE_WIDTH    = 1.8
N_BINS        = 100     # fine binning
Y_MAX_RATIO   = 2.0     # for ratio plots

def safe_array(tree, branch):
    """Read branch from uproot tree, return empty array on failure."""
    try:
        return tree[branch].array(library="np")
    except Exception as e:
        print(f"[Warning] Could not read '{branch}' from tree {getattr(tree, 'name', '')}: {e}")
        return np.empty(0)

def process_period_runs(period, target, bins, centers, charge_map, trees):
    """
    Compute per-run normalized histograms and integrals for one (period,target).
    Returns dict with: period, target, centers, list of {run,norm,integral}, mean, std.
    """
    tree = trees[period][target]
    runnums = safe_array(tree, "runnum")
    xvals   = safe_array(tree, "x")
    unique  = np.unique(runnums)
    entries = []
    for run in unique:
        charge = charge_map.get(run)
        if charge is None:
            print(f"[Warning] missing charge for {period} {target} run {run}")
            continue
        mask = (runnums == run)
        counts = np.histogram(xvals[mask], bins=bins)[0]
        norm   = counts / charge
        integ  = norm.sum()
        entries.append({"run": run, "norm": norm, "integral": integ})
    integrals = [e["integral"] for e in entries]
    mean_i = float(np.mean(integrals)) if integrals else 0.0
    std_i  = float(np.std(integrals))  if integrals else 0.0
    return {
        "period": period,
        "target": target,
        "centers": centers,
        "entries": entries,
        "mean": mean_i,
        "std": std_i
    }

def plot_normalized_yields(trees, xB_bins):
    """
    Main entrypoint: generates all figures and writes summary text file.
    """
    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # prepare bins & centers
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins    = np.linspace(xmin, xmax, N_BINS+1)
    centers = 0.5*(bins[:-1] + bins[1:])

    # 1) precompute absolute normalized histograms
    norm_hist = {p:{} for p in periods}
    for p in periods:
        for t in targets:
            x = safe_array(trees[p][t], "x")
            cnt = np.histogram(x, bins=bins)[0] if x.size else np.zeros(len(bins)-1)
            norm_hist[p][t] = cnt / CHARGE[p][t]

    # FIGURE 1: Absolute normalization
    fig1, axs1 = plt.subplots(2,3,figsize=(15,8), sharex=True)
    axs1 = axs1.flatten()
    for i,t in enumerate(targets):
        ax = axs1[i]
        peak = max((norm_hist[p][t].max() for p in periods), default=0)
        for p in periods:
            y = norm_hist[p][t]
            ax.step(centers, y, where='mid',
                    color=PERIOD_COLORS[p], linewidth=LINE_WIDTH,
                    label=p.replace("RGC_",""))
        ax.set_ylim(0, 1.2*peak)
        ax.set_title(t)
        ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("counts / nC")
        ax.legend(frameon=False, fontsize="small")
    axs1[-1].axis("off")
    plt.tight_layout(pad=2.0)
    os.makedirs("output", exist_ok=True)
    fig1.savefig("output/normalized_yields.pdf")
    plt.close(fig1)
    print("[Plot] Saved 'output/normalized_yields.pdf'")

    # FIGURE 2: Relative to Su22
    fig2, axs2 = plt.subplots(2,3,figsize=(15,8), sharex=True)
    axs2 = axs2.flatten()
    base = "RGC_Su22"
    for i,t in enumerate(targets):
        ax = axs2[i]
        bvals = norm_hist[base][t]
        for p in periods:
            y = norm_hist[p][t]
            ratio = np.where(bvals>0, y/bvals, np.nan)
            ax.step(centers, ratio, where='mid',
                    color=PERIOD_COLORS[p], linewidth=LINE_WIDTH,
                    label=p.replace("RGC_",""))
        ax.set_ylim(0, Y_MAX_RATIO)
        ax.set_title(t)
        ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("ratio / Su22")
        ax.legend(frameon=False, fontsize="small")
    axs2[-1].axis("off")
    plt.tight_layout(pad=2.0)
    fig2.savefig("output/normalized_yields_ratio.pdf")
    plt.close(fig2)
    print("[Plot] Saved 'output/normalized_yields_ratio.pdf'")

    # load run info once
    runinfo = "/u/home/thayward/.../clas12_run_info.csv"
    run_df = pd.read_csv(runinfo, header=None, comment="#", usecols=[0,1], names=["run","charge"])
    charge_map = run_df.set_index("run")["charge"].to_dict()

    # prepare summary file
    summary_path = "output/normalized_yields_integrals.txt"
    summary = []
    summary.append("Period\tTarget\tRun\tIntegral\tMean\tStd\tFlag(>2σ?)\n")

    # per-run plots for each target, in parallel per period
    for target in targets:
        # kick off one job per period
        jobs = []
        with ProcessPoolExecutor(max_workers=3) as exe:
            futures = {}
            for p in periods:
                futures[exe.submit(process_period_runs, p, target, bins, centers, charge_map, trees)] = p
            for fut in as_completed(futures):
                res = fut.result()
                period = res["period"]
                # write summary lines
                for e in res["entries"]:
                    flag = ""
                    if res["std"]>0 and abs(e["integral"]-res["mean"])>2*res["std"]:
                        flag = "YES"
                    summary.append(f"{period}\t{target}\t{e['run']}\t{e['integral']:.4f}\t"
                                   f"{res['mean']:.4f}\t{res['std']:.4f}\t{flag}\n")

                # plotting
                fig, ax = plt.subplots(1,1,figsize=(6,5), sharex=True)
                cmap = plt.get_cmap("tab20", len(res["entries"]))
                linestyles = ['solid','dashed','dashdot','dotted']
                for idx,e in enumerate(res["entries"]):
                    run = e["run"]
                    norm = e["norm"]
                    color = cmap(idx % 20)
                    ls    = linestyles[(idx//20) % len(linestyles)]
                    ax.step(centers, norm, where='mid',
                            color=color, linestyle=ls, linewidth=1.5,
                            label=str(run))
                ax.set_title(f"{period.replace('RGC_','')} {target}")
                ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("counts / nC")
                ax.legend(fontsize="x-small", ncol=2, frameon=False)
                plt.tight_layout(pad=2.0)
                outname = f"output/normalized_yields_runs_{period.lower()}_{target.lower()}.pdf"
                fig.savefig(outname)
                plt.close(fig)
                print(f"[Plot] Saved '{outname}' for {period} {target}")

    # finally write summary file
    with open(summary_path, "w") as f:
        f.writelines(summary)
    print(f"[Summary] wrote all integrals to '{summary_path}'")