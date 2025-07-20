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
linestyles, and logs every run’s integral immediately when each period finishes.
Uses a two‐pass, chunked, vectorized histogram fill to achieve true O(N_events).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---------- Configuration ----------
CHARGE = {
    "RGC_Su22": {"NH3": 3686969.636627, "C": 363715.413199, "CH2": 189723.230200, "He": 121362.943484, "ET": 586020.625415},
    "RGC_Fa22": {"NH3": 5509572.178076, "C": 2006703.362768, "CH2": 1756673.682648, "He": 275693.983222, "ET":  57332.748981},
    "RGC_Sp23": {"NH3": 1620599.347496, "C":  383030.166502, "CH2":  436816.389755, "He": 266597.436226, "ET": 171738.938831},
}
PERIOD_COLORS = {"RGC_Su22": "black", "RGC_Fa22": "blue", "RGC_Sp23": "red"}
LINE_WIDTH    = 1.8
N_BINS        = 100
Y_MAX_RATIO   = 2.0

RUNINFO = (
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "asymmetry_extraction/imports/clas12_run_info.csv"
)

# ---------- Helpers ----------
def safe_array(tree, branch):
    """
    Safely read a branch from an uproot tree, catching decompression errors.
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
    Two-pass, chunked fill of per-run histograms:
      - Pass 1: gather unique runnums
      - Pass 2: vectorized fill via np.add.at
    Returns dict with entries, mean, std.
    """
    print(f"[Worker] START {period}/{target}")
    tree = trees[period][target]

    # Pass 1: find unique runs
    unique_runs = set()
    chunk_n = 0
    for chunk in tree.iterate("runnum", step_size=5_000_000, library="np"):
        chunk_n += 1
        rn_chunk = chunk["runnum"]
        unique_runs.update(np.unique(rn_chunk))
        print(f"[Worker:{period}/{target}] Pass1 chunk {chunk_n}: saw {rn_chunk.size} events, total runs so far {len(unique_runs)}")
    runs_sorted = sorted(unique_runs)
    run_to_idx = {r:i for i,r in enumerate(runs_sorted)}
    n_runs = len(runs_sorted)
    n_bins = len(bins) - 1
    print(f"[Worker] {period}/{target} found {n_runs} unique runs")

    # allocate histograms
    hists = np.zeros((n_runs, n_bins), dtype=np.int64)

    # Pass 2: fill
    chunk_n = 0
    for chunk in tree.iterate(["runnum","x"], step_size=5_000_000, library="np"):
        chunk_n += 1
        rn_chunk = chunk["runnum"]
        x_chunk  = chunk["x"]
        print(f"[Worker:{period}/{target}] Pass2 chunk {chunk_n}: processing {rn_chunk.size} events")
        idx_run = np.fromiter((run_to_idx.get(r, -1) for r in rn_chunk), dtype=np.int64)
        idx_bin = np.digitize(x_chunk, bins) - 1
        valid = (idx_run >= 0) & (idx_bin >= 0) & (idx_bin < n_bins)
        idx_run = idx_run[valid]
        idx_bin = idx_bin[valid]
        np.add.at(hists, (idx_run, idx_bin), 1)
    print(f"[Worker] {period}/{target} histogram fill complete")

    # normalize and compute integrals
    entries = []
    integrals = []
    for i, run in enumerate(runs_sorted):
        ch = charge_map.get(run)
        if ch is None:
            print(f"[Worker:{period}/{target}] skip run {run}: no charge info")
            continue
        norm = hists[i] / ch
        integ = float(norm.sum())
        entries.append({"run": run, "norm": norm, "integral": integ})
        integrals.append(integ)
        print(f"[Worker:{period}/{target}] Run={run} integral={integ:.4f}")

    mean = float(np.mean(integrals)) if integrals else 0.0
    std  = float(np.std(integrals))  if integrals else 0.0
    print(f"[Worker] DONE {period}/{target}: mean={mean:.4f}, std={std:.4f}")
    return {
        "period": period,
        "target":  target,
        "centers": centers,
        "entries": entries,
        "mean":    mean,
        "std":     std
    }

# ---------- Main Plotting Function ----------
def plot_normalized_yields(trees, xB_bins):
    print("[Start] plot_normalized_yields()")
    os.makedirs("output", exist_ok=True)
    log_path = "output/per_run_integrals.txt"
    log = open(log_path, "w", buffering=1)
    log.write("Per-run integrals & stats\n==========================\n\n")

    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    # bins & centers
    xmin, xmax = xB_bins[0], xB_bins[-1]
    bins    = np.linspace(xmin, xmax, N_BINS+1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    print(f"[Setup] xB range {xmin}–{xmax}, {N_BINS} bins")

    # load run charges
    print(f"[Load] reading runinfo from {RUNINFO}")
    run_df = pd.read_csv(RUNINFO, header=None, comment="#", usecols=[0,1], names=["run","charge"])
    charge_map = run_df.set_index("run")["charge"].to_dict()
    print(f"[Load] {len(charge_map)} runs with charge loaded")

    # precompute absolute histograms
    print("[Compute] absolute normalized histograms")
    norm_hist = {p:{} for p in periods}
    for p in periods:
        for t in targets:
            print(f"[Compute]   {p}/{t}")
            x = safe_array(trees[p][t], "x")
            cnt = np.histogram(x, bins=bins)[0] if x.size else np.zeros(len(bins)-1)
            norm_hist[p][t] = cnt / CHARGE[p][t]

    # FIGURE 1: Absolute normalization
    print("[Figure] generating absolute normalization")
    fig1, axs1 = plt.subplots(2,3,figsize=(15,8), sharex=True)
    axs1 = axs1.flatten()
    for i, t in enumerate(targets):
        ax = axs1[i]
        peak = max(norm_hist[p][t].max() for p in periods)
        for p in periods:
            y = norm_hist[p][t]
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

    # FIGURE 2: Ratio to Su22
    print("[Figure] generating ratio to Su22")
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

    # FIGURES 3–7: per-run for each target
    for target in targets:
        print(f"[Per-run] scheduling per-run for {target}")
        args = [(p, target, bins, centers, charge_map, trees) for p in periods]

        with ProcessPoolExecutor(max_workers=3) as executor:
            futures = {executor.submit(process_period_runs, *a): a[0] for a in args}
            results = {}
            for fut in as_completed(futures):
                period = futures[fut]
                res = fut.result()
                results[period] = res

                # immediate log
                log.write(f"=== Period {period}, Target {target} ===\n")
                log.write(f" Mean={res['mean']:.4f}, Std={res['std']:.4f}\n")
                for e in res["entries"]:
                    log.write(f"  Run={e['run']}, Integral={e['integral']:.4f}\n")
                log.write("\n")
                print(f"[Log] {period}/{target} done: mean={res['mean']:.4f}, std={res['std']:.4f}")

        # draw canvas
        print(f"[Figure] plotting per-run yields for {target}")
        fig, axes = plt.subplots(1,3,figsize=(18,5), sharex=True, sharey=True)
        for ax, p in zip(axes, periods):
            res = results[p]
            for idx, e in enumerate(res["entries"]):
                cmap  = plt.get_cmap("tab20")
                style = ['solid','dashed','dashdot','dotted'][(idx//20)%4]
                color = cmap(idx%20)
                ax.step(res["centers"], e["norm"], where='mid',
                        color=color, linestyle=style, linewidth=1.5,
                        label=str(e["run"]))
            ax.set_title(f"{p.replace('RGC_','')} {target}")
            ax.set_xlabel(r"$x_{B}$"); ax.set_ylabel("counts / nC")
            ax.legend(fontsize="x-small", ncol=2, frameon=False)
        plt.tight_layout(pad=2)
        out = f"output/normalized_yields_runs_{target.lower()}.pdf"
        fig.savefig(out); plt.close(fig)
        print(f"[Saved] {out}\n")

    log.close()
    print(f"[Finished] All per-run logs written to {log_path}")
    print("[Done] plot_normalized_yields()")