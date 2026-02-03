#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rgc_rga_mx2_comparison.py

Make a 1x3 canvas comparing Mx2 distributions for:
  - RGC Su22 NH3
  - RGC Su22 C
  - RGA Fa18 Inb H2

Histograms are charge-normalized:
  h_norm(bin) = N(bin) / Q_total

Subplot 1: NH3/charge, C/charge, H2/charge (all overlayed)
Subplot 2: Ratio (NH3/charge) / (C/charge)
Subplot 3: Difference (NH3/charge - C/charge) overlayed with H2/charge

Optional:
  --short : restrict each ROOT tree to only the first 6 run numbers (ascending),
            but keep all events belonging to those runs. Charge normalization uses
            only those runs as well.

Output:
  output/rgc_rga_mx2_comparison.png

Dependencies:
  python3, uproot, numpy, matplotlib

Example:
  python3 rgc_rga_mx2_comparison.py
  python3 rgc_rga_mx2_comparison.py --short
"""

import os
import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------
# User inputs
# -------------------------------------------------------------------------
NH3_ROOT = "/volatile/clas12/thayward/rgc_enpi+_Mx2_study/rgc_su22_inb_NH3.root"
C_ROOT   = "/volatile/clas12/thayward/rgc_enpi+_Mx2_study/rgc_su22_inb_C.root"
H2_ROOT  = "/volatile/clas12/thayward/rgc_enpi+_Mx2_study/rga_fa18_inb_H2.root"

TREE_NAME = "PhysicsEvents"

RUN_INFO_CSV = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

OUT_PNG = "output/rgc_rga_mx2_comparison.png"

# Histogram settings
MX2_MIN = 0.0
MX2_MAX = 4.0
NBINS   = 200


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
def parse_args(argv):
    short = False
    for a in argv[1:]:
        if a == "--short":
            short = True
        else:
            raise RuntimeError(f"FATAL: Unknown argument: {a}")
        #endif
    #endfor
    return short


def read_run_charge_map(csv_path):
    """
    Read a CSV file with lines like:
      runnum,charge,...

    Returns:
      dict[int, float] mapping runnum -> accumulated charge
    """
    if not os.path.isfile(csv_path):
        raise RuntimeError(f"FATAL: CSV not found: {csv_path}")

    run_charge = {}

    with open(csv_path, "r", encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            s = line.strip()
            if len(s) == 0:
                continue
            if s.startswith("#"):
                continue

            parts = [p.strip() for p in s.split(",")]
            if len(parts) < 2:
                raise RuntimeError(
                    f"FATAL: Malformed CSV line {lineno} in {csv_path}: '{line.rstrip()}'"
                )

            try:
                runnum = int(parts[0])
            except Exception:
                raise RuntimeError(
                    f"FATAL: Could not parse runnum on line {lineno} in {csv_path}: '{parts[0]}'"
                )

            try:
                charge = float(parts[1])
            except Exception:
                raise RuntimeError(
                    f"FATAL: Could not parse charge on line {lineno} in {csv_path}: '{parts[1]}'"
                )

            run_charge[runnum] = charge
        #endfor

    if len(run_charge) == 0:
        raise RuntimeError(f"FATAL: No run/charge entries parsed from {csv_path}")

    return run_charge


def load_tree_arrays(root_path, tree_name):
    """
    Load Mx2 and runnum arrays from a ROOT TTree.
    Returns:
      mx2_np (np.ndarray float64)
      run_np (np.ndarray int64)
    """
    if not os.path.isfile(root_path):
        raise RuntimeError(f"FATAL: ROOT file not found: {root_path}")

    with uproot.open(root_path) as f:
        if tree_name not in f:
            raise RuntimeError(
                f"FATAL: TTree '{tree_name}' not found in file: {root_path}"
            )

        tree = f[tree_name]

        needed = ["Mx2", "runnum"]
        for br in needed:
            if br not in tree.keys():
                raise RuntimeError(
                    f"FATAL: Branch '{br}' not found in TTree '{tree_name}' in file: {root_path}"
                )
            #endif
        #endfor

        # Load only the two needed branches
        mx2 = tree["Mx2"].array(library="np")
        run = tree["runnum"].array(library="np")

    mx2_np = np.asarray(mx2, dtype=np.float64)
    run_np = np.asarray(run, dtype=np.int64)

    if mx2_np.shape[0] != run_np.shape[0]:
        raise RuntimeError(
            f"FATAL: Array length mismatch in {root_path}: Mx2 has {mx2_np.shape[0]} entries, runnum has {run_np.shape[0]} entries"
        )

    return mx2_np, run_np


def restrict_to_first_n_runs(mx2, run, n_runs, label):
    """
    Restrict arrays to only events belonging to the first n_runs unique run numbers (ascending).
    Keeps all events for those runs.

    Returns:
      mx2_cut, run_cut, runs_kept (np.ndarray of unique runnums kept)
    """
    runs_unique = np.unique(run)
    if runs_unique.size == 0:
        raise RuntimeError(f"FATAL: No run numbers found in {label}")

    runs_kept = runs_unique[:n_runs]

    # Vectorized membership test
    mask = np.isin(run, runs_kept)
    mx2_cut = mx2[mask]
    run_cut = run[mask]

    return mx2_cut, run_cut, runs_kept


def total_charge_for_runs(runs_unique, run_charge_map, label_for_errors):
    """
    Sum charges for an array of unique run numbers.
    Fail fast if any run is missing in the CSV.
    """
    runs_unique = np.asarray(runs_unique, dtype=np.int64)

    # Find missing runs using vectorized lookup over Python dict via list comprehension
    # (unique runs is small, so this is already cheap)
    missing = [int(r) for r in runs_unique if int(r) not in run_charge_map]

    if len(missing) > 0:
        missing_sorted = sorted(set(missing))
        msg = (
            f"FATAL: Missing {len(missing_sorted)} run(s) in run-info CSV for {label_for_errors}.\n"
            f"Missing runnum values: {missing_sorted}\n"
            f"Fix: add these run numbers to {RUN_INFO_CSV} (or remove them from the ROOT file)."
        )
        raise RuntimeError(msg)

    total = 0.0
    for r in runs_unique:
        total += run_charge_map[int(r)]
    #endfor

    if not np.isfinite(total) or total <= 0.0:
        raise RuntimeError(
            f"FATAL: Computed non-positive or non-finite total charge for {label_for_errors}: {total}"
        )

    return total


def hist_counts(mx2_values, bin_edges):
    """
    Histogram Mx2 with fixed bins in [MX2_MIN, MX2_MAX].
    Drops events outside [MX2_MIN, MX2_MAX].
    """
    mask = (mx2_values >= MX2_MIN) & (mx2_values <= MX2_MAX)
    mx2_in = mx2_values[mask]
    counts, _ = np.histogram(mx2_in, bins=bin_edges)
    return counts.astype(np.float64)


def safe_ratio(num, den):
    """
    Compute num/den bin-by-bin, returning 0 where den == 0.
    """
    out = np.zeros_like(num, dtype=np.float64)
    mask = den != 0.0
    out[mask] = num[mask] / den[mask]
    return out


def step_plot(ax, bin_edges, y, label, linewidth=1.2):
    ax.step(bin_edges[:-1], y, where="post", label=label, linewidth=linewidth)


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
def main():
    short_mode = parse_args(sys.argv)

    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)

    run_charge_map = read_run_charge_map(RUN_INFO_CSV)

    # Load arrays (only 2 branches each)
    nh3_mx2, nh3_run = load_tree_arrays(NH3_ROOT, TREE_NAME)
    c_mx2,   c_run   = load_tree_arrays(C_ROOT,   TREE_NAME)
    h2_mx2,  h2_run  = load_tree_arrays(H2_ROOT,  TREE_NAME)

    # Optional short restriction
    if short_mode:
        nh3_mx2, nh3_run, nh3_runs_kept = restrict_to_first_n_runs(nh3_mx2, nh3_run, 6, "NH3")
        c_mx2,   c_run,   c_runs_kept   = restrict_to_first_n_runs(c_mx2,   c_run,   6, "C")
        h2_mx2,  h2_run,  h2_runs_kept  = restrict_to_first_n_runs(h2_mx2,  h2_run,  6, "H2")

        nh3_runs_unique = nh3_runs_kept
        c_runs_unique   = c_runs_kept
        h2_runs_unique  = h2_runs_kept
    else:
        nh3_runs_unique = np.unique(nh3_run)
        c_runs_unique   = np.unique(c_run)
        h2_runs_unique  = np.unique(h2_run)
    #endif

    # Total charge for each sample
    nh3_charge = total_charge_for_runs(nh3_runs_unique, run_charge_map, f"NH3 file '{NH3_ROOT}'")
    c_charge   = total_charge_for_runs(c_runs_unique,   run_charge_map, f"C file '{C_ROOT}'")
    h2_charge  = total_charge_for_runs(h2_runs_unique,  run_charge_map, f"H2 file '{H2_ROOT}'")

    # Fixed binning for all plots
    bin_edges = np.linspace(MX2_MIN, MX2_MAX, NBINS + 1)

    # Raw counts per bin
    nh3_counts = hist_counts(nh3_mx2, bin_edges)
    c_counts   = hist_counts(c_mx2,   bin_edges)
    h2_counts  = hist_counts(h2_mx2,  bin_edges)

    # Charge-normalized histograms (counts per unit charge)
    nh3_norm = nh3_counts / nh3_charge
    c_norm   = c_counts   / c_charge
    h2_norm  = h2_counts  / h2_charge

    # Ratio and difference
    ratio_nh3_over_c = safe_ratio(nh3_norm, c_norm)
    diff_nh3_minus_c = nh3_norm - c_norm

    # Plotting
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

    x_label = r"$M_{x}^{2}$ (GeV$^{2}$)"
    y_label = "Counts / Charge"

    # 1) NH3 vs C vs H2 (all charge-normalized)
    ax = axes[0]
    step_plot(ax, bin_edges, nh3_norm, "RGC Su22 NH3")
    step_plot(ax, bin_edges, c_norm,   "RGC Su22 C")
    step_plot(ax, bin_edges, h2_norm,  "RGA Fa18 Inb H2")
    ax.set_title("Charge-normalized Mx2")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(MX2_MIN, MX2_MAX)
    ax.legend(fontsize=10)

    # 2) Ratio (NH3/charge) / (C/charge)
    ax = axes[1]
    step_plot(ax, bin_edges, ratio_nh3_over_c, "(NH3/charge) / (C/charge)")
    ax.axhline(1.0, linewidth=1.0)
    ax.set_title("Ratio of charge-normalized spectra")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Ratio")
    ax.set_xlim(MX2_MIN, MX2_MAX)
    ax.legend(fontsize=10)

    # 3) Difference (NH3/charge - C/charge) and overlay H2/charge
    ax = axes[2]
    step_plot(ax, bin_edges, diff_nh3_minus_c, "NH3/charge - C/charge")
    step_plot(ax, bin_edges, h2_norm,          "H2/charge")
    ax.axhline(0.0, linewidth=1.0)
    ax.set_title("Difference spectrum vs H2 reference")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(MX2_MIN, MX2_MAX)
    ax.legend(fontsize=10)

    if short_mode:
        fig.suptitle("Mx2 comparison (SHORT mode: first 6 runs only)", fontsize=14)
    #endif

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"Wrote: {OUT_PNG}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
    #endif