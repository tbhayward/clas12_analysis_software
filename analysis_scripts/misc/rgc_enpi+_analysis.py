#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

# Toggle for quick debugging (process only first 5 runs per period)
QUICK_RUN = True

# Maximum run number to include
MAX_RUNNUM = 17768

# Path to the CSV of run charges
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output directory and filenames
OUTPUT_DIR         = "output/enpi+"
THREE_PANEL_OUTPUT = os.path.join(OUTPUT_DIR, "three_panel_Mx2_comparison.pdf")

# File lists: (filepath, label)
NH3_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_NH3_epi+.root", "Su22-NH3"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+.root", "Fa22-NH3"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_NH3_epi+.root", "Sp23-NH3"),
]
C_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_C_epi+.root", "Su22-C"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_C_epi+.root", "Fa22-C"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_C_epi+.root", "Sp23-C"),
]
H2_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_sp19_inb_H2_epi+.root", "Sp19-H2"),
]

# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------

def parse_run_charges(csv_path):
    """
    Read clas12_run_info.csv and return a dict { runnum: charge }, skipping any run > MAX_RUNNUM.
    Lines beginning with '#' are ignored.
    """
    run_charges = {}
    if not os.path.exists(csv_path):
        print(f"[ERROR] CSV not found: {csv_path}")
        return run_charges

    with open(csv_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            if len(parts) < 2:
                continue
            try:
                run = int(parts[0])
                charge = float(parts[1])
            except ValueError:
                continue
            if run > MAX_RUNNUM:
                continue
            run_charges[run] = charge

    print(f"[DEBUG] Parsed {len(run_charges)} runs ≤ {MAX_RUNNUM} from CSV.")
    return run_charges


def load_tree(filepath):
    """
    Open a ROOT file and return the 'PhysicsEvents' TTree.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Could not find file {filepath}")
    tree = uproot.open(filepath)["PhysicsEvents"]
    print(f"[DEBUG] Loaded TTree 'PhysicsEvents' from {filepath} with {tree.num_entries} entries.")
    return tree


def compute_charge_and_mask(tree, run_charges, quick=False):
    """
    For a given TTree, return (Q_tot, masked_Mx2_array) where:
      - Q_tot is the sum of charges for unique runs ≤ MAX_RUNNUM
        (only first 5 if quick=True),
      - masked_Mx2_array are the Mx2 values for those runs with Mx2 ≤ 2.0.
    """
    run_arr = tree["runnum"].array(library="np").astype(int)
    mx2_arr = tree["Mx2"].array(library="np")
    if run_arr.size == 0:
        return 0.0, np.array([], dtype=float)

    # Filter by runnum ≤ MAX_RUNNUM and Mx2 ≤ 2.0
    valid = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
    run_arr = run_arr[valid]
    mx2_arr = mx2_arr[valid]
    if run_arr.size == 0:
        return 0.0, np.array([], dtype=float)

    # Determine selected runs
    if quick:
        seen = []
        for r in run_arr:
            if r not in seen:
                seen.append(r)
            if len(seen) >= 5:
                break
        selected_runs = set(seen)
    else:
        selected_runs = set(np.unique(run_arr))

    # Sum total charge
    Q_tot = 0.0
    for r in selected_runs:
        q = run_charges.get(r, 0.0)
        if q <= 0.0:
            print(f"    [WARNING] run {r} has zero or missing charge.")
        Q_tot += q

    # Build mask to keep only those runs
    mask = np.isin(run_arr, list(selected_runs))
    return Q_tot, mx2_arr[mask]


def make_normalized_Mx2_plots(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 1×3 figure:
      - Panel 1: NH₃ (Su22, Fa22, Sp23) + H₂
      - Panel 2:  C   (Su22, Fa22, Sp23) + H₂
      - Panel 3: Differences [NH₃ – C] for each period (normalized to match 0–0.5), plus H₂
    Excludes runs > MAX_RUNNUM, applies Mx² ≤ 2.0 mask.
    If QUICK_RUN=True, only first 5 unique runnums per tree are used.
    """
    bins = np.linspace(0.0, 1.5, 101)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    # -------------------------------------------------------------
    # Panel 1: NH₃ + H₂
    # -------------------------------------------------------------
    for filepath, label in nh3_files:
        tree = load_tree(filepath)
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0.0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        norm_counts = counts.astype(float) / Q
        axes[0].step(bins[:-1], norm_counts, where='post', label=label)
    for filepath, label in h2_files:
        tree = load_tree(filepath)
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0.0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        norm_counts = counts.astype(float) / Q
        axes[0].step(
            bins[:-1], norm_counts, where='post',
            color='k', linestyle='-', label=label
        )
    axes[0].set(
        title="NH₃ + H₂: Normalized Mx²",
        xlabel=r"$M_x^2$ [GeV$^2$]",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    axes[0].legend(loc='upper right', fontsize='small')

    # -------------------------------------------------------------
    # Panel 2: C + H₂
    # -------------------------------------------------------------
    for filepath, label in c_files:
        tree = load_tree(filepath)
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0.0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        norm_counts = counts.astype(float) / Q
        axes[1].step(bins[:-1], norm_counts, where='post', label=label)
    for filepath, label in h2_files:
        tree = load_tree(filepath)
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0.0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        norm_counts = counts.astype(float) / Q
        axes[1].step(
            bins[:-1], norm_counts, where='post',
            color='k', linestyle='-', label=label
        )
    axes[1].set(
        title="C + H₂: Normalized Mx²",
        xlabel=r"$M_x^2$ [GeV$^2$]",
        xlim=(0.0, 1.5)
    )
    axes[1].legend(loc='upper right', fontsize='small')

    # -------------------------------------------------------------
    # Panel 3: (NH₃ – C) differences, plus H₂
    # -------------------------------------------------------------
    # Build period-level NH₃ and C histograms for each period
    period_names = ["Su22", "Fa22", "Sp23"]
    nh3_hist_dict = {}
    c_hist_dict = {}
    for period_label in period_names:
        # Find corresponding NH₃ and C filepaths
        nh3_fp = next(fp for fp, lbl in nh3_files if lbl.startswith(period_label))
        c_fp   = next(fp for fp, lbl in c_files   if lbl.startswith(period_label))

        # NH₃
        tree_n = load_tree(nh3_fp)
        Q_n, mx2_n = compute_charge_and_mask(tree_n, run_charges, quick=QUICK_RUN)
        if Q_n > 0.0:
            counts_n, _ = np.histogram(mx2_n, bins=bins)
            nh3_hist_dict[period_label] = counts_n.astype(float) / Q_n
        else:
            nh3_hist_dict[period_label] = np.zeros(len(bins)-1)

        # C
        tree_c = load_tree(c_fp)
        Q_c, mx2_c = compute_charge_and_mask(tree_c, run_charges, quick=QUICK_RUN)
        if Q_c > 0.0:
            counts_c, _ = np.histogram(mx2_c, bins=bins)
            c_hist_dict[period_label] = counts_c.astype(float) / Q_c
        else:
            c_hist_dict[period_label] = np.zeros(len(bins)-1)

    # Compute H₂ histogram once
    h2_label = h2_files[0][1]
    tree_h2  = load_tree(h2_files[0][0])
    Q_h2, mx2_h2 = compute_charge_and_mask(tree_h2, run_charges, quick=QUICK_RUN)
    if Q_h2 > 0.0:
        counts_h2, _ = np.histogram(mx2_h2, bins=bins)
        h2_norm = counts_h2.astype(float) / Q_h2
    else:
        h2_norm = np.zeros(len(bins)-1)

    # Compute differences NH₃ – scaled C for each period
    diff_dict = {}
    for period_label in period_names:
        nh3_vals = nh3_hist_dict[period_label]
        c_vals   = c_hist_dict[period_label]

        # Select bins with center < 0.5
        mask_0_5 = (bin_centers >= 0.0) & (bin_centers < 0.5)
        nh3_sum = nh3_vals[mask_0_5].sum()
        c_sum   = c_vals[mask_0_5].sum()

        if c_sum > 0.0:
            scale = nh3_sum / c_sum
        else:
            scale = 0.0
        print(f"[DEBUG] Period {period_label}: NH₃ sum(0–0.5)={nh3_sum:.3f}, C sum(0–0.5)={c_sum:.3f}, scale={scale:.3f}")

        diff = nh3_vals - (c_vals * scale)
        diff_dict[period_label] = diff

    # Plot differences and H₂ on panel 3
    linestyles = ['-', '--', '-.', ':']
    for idx, period_label in enumerate(period_names):
        diff = diff_dict[period_label]
        style = linestyles[idx % len(linestyles)]
        axes[2].step(
            bins[:-1], diff, where='post',
            linestyle=style, label=f"{period_label} NH₃–C"
        )

    # Plot H₂ on the same axes
    axes[2].step(
        bins[:-1], h2_norm, where='post',
        color='k', linestyle='-', label=h2_label
    )

    axes[2].set(
        title="(NH₃ – C) per period, and H₂",
        xlabel=r"$M_x^2$ [GeV$^2$]",
        xlim=(0.0, 1.5),
        ylim=(0.0, 0.03)
    )
    axes[2].legend(loc='upper right', fontsize='small')

    # Save the three-panel figure
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    # Load run charges, skipping runnum > MAX_RUNNUM
    run_charges = parse_run_charges(CSV_PATH)

    # Create the three‐panel comparison plot
    make_normalized_Mx2_plots(NH3_FILES, C_FILES, H2_FILES, run_charges, THREE_PANEL_OUTPUT)

if __name__ == "__main__":
    main()