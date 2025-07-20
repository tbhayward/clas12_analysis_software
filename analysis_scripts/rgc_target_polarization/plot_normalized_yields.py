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

# ---------- Configuration ----------
# Absolute accumulated charge (in nC) per period and target
CHARGE = {
    "RGC_Su22": {"NH3": 3686969.636627, "C": 363715.413199, "CH2": 189723.230200, "He": 266575.708348, "ET": 586020.625415},
    "RGC_Fa22": {"NH3": 5509572.178076, "C": 2006703.362768, "CH2": 1833777.957336, "He": 275693.983222, "ET":  57332.748981},
    "RGC_Sp23": {"NH3": 1620599.347496, "C":  383030.166502, "CH2":  436816.389755, "He": 266597.436226, "ET": 171738.938831},
}

PERIOD_COLORS = {"RGC_Su22": "black", "RGC_Fa22": "blue", "RGC_Sp23": "red"}
LINE_WIDTH    = 1.8
N_BINS        = 100   # fine binning
Y_MAX_RATIO   = 2.0   # ratio plot y‐limit

# Location of per-run charge info
RUNINFO = (
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "asymmetry_extraction/imports/clas12_run_info.csv"
)

# ---------- Helpers ----------
def safe_array(tree, branch):
    """
    Safely read a branch from an uproot tree, catching decompression errors.
    Returns an empty array on failure.
    """
    try:
        arr = tree[branch].array(library="np")
        print(f"[Read OK] {tree.name}.{branch} → {arr.size} entries")
        return arr
    except Exception as e:
        print(f"[Warning] failed to read '{branch}' from {tree.name}: {e}")
        return np.empty(0)

def process_period_runs(period, target, bins, centers, charge_map, trees):
    """
    Worker for per-run plotting. Collects each run’s normalized histogram + integral,
    then computes that period’s mean & std of integrals.
    """
    print(f"[Worker] Start processing runs for Period={period}, Target={target}")
    tree = trees[period][target]
    rn   = safe_array(tree, "runnum")
    xv   = safe_array(tree, "x")
    runs = np.unique(rn)
    entries = []
    for run in runs:
        ch = charge_map.get(run)
        if ch is None:
            print(f"[Worker]  skip run {run} (no charge info)")
            continue
        mask = (rn == run)
        cnt  = np.histogram(xv[mask], bins=bins)[0]
        norm = cnt / ch
        integ = float(norm.sum())
        print(f"[Worker]   Run={run}: integral={integ:.4f}")
        entries.append({"run": run, "norm": norm, "integral": integ})
    ints = [e["integral"] for e in entries]
    mean = float(np.mean(ints)) if ints else 0.0
    std  = float(np.std(ints))  if ints else 0.0
    print(f"[Worker] Finished {period}/{target}: mean={mean:.4f}, std={std:.4f}")
    return {"period":period, "target":target,
            "centers":centers, "entries":entries,
            "mean":mean, "std":std}

# ---------- Main Plotting Function ----------
def plot_normalized_yields(trees, xB_bins):
    print("[Start] plot_normalized_yields()")
    # ensure output dir & open log file
    os.makedirs("output", exist_ok=True)
    log_path = "output/per_run_integrals.txt"
    log = open(log_path, "w", buffering=1)
    log.write("Per-run integrals & stats\n")
    log.write("==========================\n\n")

    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # prepare bins & centers
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins    = np.linspace(xmin, xmax, N_BINS+1)
    centers = 0.5*(bins[:-1]+bins[1:])
    print(f"[Setup] xB bins: {xmin} → {xmax}, {N_BINS} bins")

    # load per-run charge info
    print(f"[Load] Reading runinfo from {RUNINFO}")
    run_df = pd.read_csv(RUNINFO, header=None, comment="#",
                         usecols=[0,1], names=["run","charge"])
    charge_map = run_df.set_index("run")["charge"].to_dict()
    print(f"[Load] {len(charge_map)} runs with charge found")

    # precompute absolute normalized histograms
    print("[Compute] Precomputing absolute normalized histograms")
    norm_hist = {p:{} for p in periods}
    for p in periods:
        for t in targets:
            print(f"[Compute]   Period={p}, Target={t}")
            x = safe_array(trees[p][t], "x")
            cnt = np.histogram(x, bins=bins)[0] if x.size else np.zeros(len(bins)-1)
            norm_hist[p][t] = cnt / CHARGE[p][t]

    # ---------- FIGURE 1: Absolute normalization ----------
    print("[Figure] Generating absolute normalization plot")
    fig1, axs1 = plt.subplots(2,3,figsize=(15,8), sharex=True)
    axs1 = axs1.flatten()
    for i, t in enumerate(targets):
        ax = axs1[i]
        peak = 0.0
        for p in periods:
            y = norm_hist[p][t]
            peak = max(peak, y.max() if y.size else 0.0)
            ax.step(centers, y, where='mid',
                    color=PERIOD_COLORS[p], linewidth=LINE_WIDTH,
                    label=p.replace("RGC_",""))
        ax.set_ylim(0, 1.2*peak)
        ax.set_title(f"Absolute: {t}")
        ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("counts / nC")
        ax.legend(frameon=False, fontsize="small")
    axs1[-1].axis("off")
    plt.tight_layout(pad=2)
    out1 = "output/normalized_yields.pdf"
    fig1.savefig(out1); plt.close(fig1)
    print(f"[Saved] {out1}")

    # ---------- FIGURE 2: Ratio to Su22 ----------
    print("[Figure] Generating ratio-to-Su22 plot")
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
        ax.set_title(f"Ratio /Su22: {t}")
        ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("ratio")
        ax.legend(frameon=False, fontsize="small")
    axs2[-1].axis("off")
    plt.tight_layout(pad=2)
    out2 = "output/normalized_yields_ratio.pdf"
    fig2.savefig(out2); plt.close(fig2)
    print(f"[Saved] {out2}")

    # ---------- FIGUREs 3–7: per-run for each target ----------
    for target in targets:
        print(f"[Per-run] Starting per-run panels for target {target}")
        # assemble args
        args = [(p, target, bins, centers, charge_map, trees) for p in periods]
        print(f"[Parallel] scheduling {len(args)} tasks for per-run {target}")
        with ProcessPoolExecutor(max_workers=3) as exe:
            futures = {exe.submit(process_period_runs, *a): a[0] for a in args}
            results = {}
            for fut in as_completed(futures):
                period = futures[fut]
                res = fut.result()
                results[period] = res
                print(f"[Parallel] Completed data for {period}/{target}")

        # now draw the 1×3 canvas
        print(f"[Figure] Plotting per-run yields for {target}")
        fig, axes = plt.subplots(1,3,figsize=(18,5), sharex=True, sharey=True)
        for ax, p in zip(axes, periods):
            res = results[p]
            log.write(f"=== Period {p}, Target {target} ===\n")
            log.write(f" Mean={res['mean']:.4f}, Std={res['std']:.4f}\n")
            for idx, e in enumerate(res["entries"]):
                run   = e["run"]
                norm  = e["norm"]
                integ = e["integral"]
                print(f"[Plot]   {p}/{target} run={run}, integral={integ:.4f}")
                log.write(f"  Run={run}, Integral={integ:.4f}\n")
                cmap  = plt.get_cmap("tab20")
                style = ['solid','dashed','dashdot','dotted'][(idx//20)%4]
                color = cmap(idx%20)
                ax.step(centers, norm, where='mid',
                        color=color, linestyle=style, linewidth=1.5,
                        label=str(run))
            log.write("\n")
            ax.set_title(f"{p.replace('RGC_','')} {target}")
            ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("counts / nC")
            ax.legend(fontsize="x-small", ncol=2, frameon=False)

        plt.tight_layout(pad=2)
        out = f"output/normalized_yields_runs_{target.lower()}.pdf"
        fig.savefig(out); plt.close(fig)
        print(f"[Saved] {out}\n")

    log.close()
    print(f"[Finished] All per-run integrals & stats logged to {log_path}")
    print("[Done] plot_normalized_yields()")