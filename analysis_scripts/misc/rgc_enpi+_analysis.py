#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

# Toggle for quick debugging (process only first 5 runs per period)
QUICK_RUN = False

# Maximum run number to include
MAX_RUNNUM = 17768

# Path to the CSV of run charges
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output directory and filenames
OUTPUT_DIR             = "output/enpi+"
SP23C_COMPARE_PLOT     = os.path.join(OUTPUT_DIR, "Sp23C_period_vs_runbyrun.pdf")
NORMALIZED_OUTPUT      = os.path.join(OUTPUT_DIR, "normalized_Mx2_charges.pdf")

# Sp23-C ROOT file path
SP23C_FILE = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_C_epi+.root"

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


def compute_period_hist(tree, run_charges, bins):
    """
    Build a single “period-level” Mx² histogram for Sp23-C:
      1) Mask to runnum ≤ MAX_RUNNUM and Mx² ≤ 2.0
      2) Combine all selected events, sum all corresponding run charges
      3) Return (normalized_counts, bin_edges) and print debug info.
    """
    run_arr = tree["runnum"].array(library="np").astype(int)
    mx2_arr = tree["Mx2"].array(library="np")

    # 1) Filter out runs > MAX_RUNNUM and Mx² > 2.0
    mask = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
    run_filtered = run_arr[mask]
    mx2_filtered = mx2_arr[mask]

    if run_filtered.size == 0:
        print("[WARNING] No events remain after filtering for period-level.")
        return np.zeros(len(bins)-1), bins

    # 2) Identify unique runs, sum their charges
    unique_runs = np.unique(run_filtered)
    print("Sp23-C period includes runs:", unique_runs.tolist())

    Q_tot = 0.0
    for r in unique_runs:
        q = run_charges.get(r, 0.0)
        if q <= 0.0:
            print(f"  → WARNING: run {r} has zero or missing charge in CSV.")
        Q_tot += q
    print(f"Total combined charge for Sp23-C period: {Q_tot:.3f} nC")

    # 3) Histogram all selected Mx² values
    counts, edges = np.histogram(mx2_filtered, bins=bins)
    norm_counts = counts.astype(float) / Q_tot
    return norm_counts, edges


def compute_runbyrun_hist(tree, run_charges, bins, quick=False):
    """
    Build run-by-run Mx² histograms for Sp23-C:
      1) Mask to runnum ≤ MAX_RUNNUM and Mx² ≤ 2.0
      2) For each unique run, histogram that run's events and normalize by that run's charge
      3) Return dict { runnum: (normalized_counts, bin_edges) } and print debug info.
      If quick=True, only first 5 unique runs encountered.
    """
    run_arr = tree["runnum"].array(library="np").astype(int)
    mx2_arr = tree["Mx2"].array(library="np")

    # 1) Apply same filter
    mask = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
    run_filtered = run_arr[mask]
    mx2_filtered = mx2_arr[mask]

    if run_filtered.size == 0:
        print("[WARNING] No events remain after filtering for run-by-run.")
        return {}

    # Determine unique runs in encountered order
    if quick:
        seen = []
        for r in run_filtered:
            if r not in seen:
                seen.append(r)
            if len(seen) >= 5:
                break
        unique_runs = seen
    else:
        unique_runs = sorted(np.unique(run_filtered))

    print("Sp23-C run-by-run will process runs:", unique_runs)

    result = {}
    for runnum in unique_runs:
        mask_run = (run_filtered == runnum)
        vals = mx2_filtered[mask_run]
        Q_run = run_charges.get(runnum, 0.0)
        if Q_run <= 0.0:
            print(f"  → WARNING: run {runnum} has zero or missing charge, skipping.")
            continue
        if vals.size == 0:
            print(f"  → WARNING: run {runnum} has no events in Mx² ≤ 2.0, skipping.")
            continue

        counts, edges = np.histogram(vals, bins=bins)
        norm_counts = counts.astype(float) / Q_run
        result[runnum] = (norm_counts, edges)

    return result


def make_sp23c_comparison_plot(tree, run_charges, output_path):
    """
    Create a figure comparing:
       * Left: period-level Sp23-C Mx² (all runs combined)
       * Right: run-by-run Sp23-C Mx²
    Saves to output_path and prints the peak values.
    """
    # Common bin edges: 0 → 1.5 GeV² in 100 bins
    bins = np.linspace(0.0, 1.5, 101)

    # 1) Compute period-level histogram
    period_norm, bin_edges = compute_period_hist(tree, run_charges, bins)
    max_period = period_norm.max() if period_norm.size > 0 else 0.0

    # 2) Compute run-by-run histograms
    runbyrun_dict = compute_runbyrun_hist(tree, run_charges, bins, quick=QUICK_RUN)
    if runbyrun_dict:
        max_runbyrun = max(norm.max() for norm, _ in runbyrun_dict.values())
    else:
        max_runbyrun = 0.0

    # 3) Plot side-by-side
    fig, (ax_p, ax_r) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    # Period-level on the left (solid red)
    ax_p.step(
        bins[:-1], period_norm, where='post',
        color='red', linewidth=2, label="Sp23-C (period)"
    )
    ax_p.set(
        title="Sp23-C: Period-Level Mx²",
        xlabel=r"$M_x^2$ [GeV$^2$]",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    ax_p.legend(loc='upper right', fontsize='small')

    # Run-by-run on the right (various linestyles)
    linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 2)), (0, (3, 5, 1, 5))]
    for idx, (runnum, (norm_counts, _)) in enumerate(runbyrun_dict.items()):
        style = linestyles[idx % len(linestyles)]
        ax_r.step(
            bins[:-1], norm_counts, where='post',
            linestyle=style, label=str(runnum)
        )
    ax_r.set(
        title="Sp23-C: Run-by-Run Mx²",
        xlabel=r"$M_x^2$ [GeV$^2$]",
        xlim=(0.0, 1.5)
    )
    ax_r.legend(loc='upper left', fontsize='xx-small', ncol=3)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    # 4) Print the numerical maxima
    print(f"Max (Sp23-C period)     = {max_period:.3f} events/nC")
    print(f"Max (Sp23-C run-by-run) = {max_runbyrun:.3f} events/nC")


def make_normalized_Mx2_plots(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 1×2 figure:
      - Left: NH3 (Su22, Fa22, Sp23) + H2
      - Right: C   (Su22, Fa22, Sp23) + H2
    Exclude runs > MAX_RUNNUM, apply Mx² ≤ 2.0 mask.
    If QUICK_RUN=True, only first 5 unique runs per tree are used.
    """
    bins = np.linspace(0.0, 1.5, 101)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    def compute_charge_and_mask(tree, run_charges, quick=False):
        run_arr = tree["runnum"].array(library="np").astype(int)
        mx2_arr = tree["Mx2"].array(library="np")
        if run_arr.size == 0:
            return 0.0, np.array([], dtype=float)

        # Filter by runno ≤ MAX_RUNNUM and Mx² ≤ 2.0
        valid = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
        run_arr = run_arr[valid]
        mx2_arr = mx2_arr[valid]
        if run_arr.size == 0:
            return 0.0, np.array([], dtype=float)

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

        Q_tot = sum(run_charges.get(r, 0.0) for r in selected_runs)
        mask = np.isin(run_arr, list(selected_runs))
        return Q_tot, mx2_arr[mask]

    # --- Left subplot: NH3 + H2 ---
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

    # --- Right subplot:  C + H2 ---
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

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    # 1) Load run charges, skipping runnum > MAX_RUNNUM
    run_charges = parse_run_charges(CSV_PATH)

    # 2) Create the period vs. run-by-run comparison plot for Sp23-C
    sp23c_tree = load_tree(SP23C_FILE)
    make_sp23c_comparison_plot(sp23c_tree, run_charges, SP23C_COMPARE_PLOT)

    # 3) (Optional) Generate the two-panel normalized Mx² plot for all targets
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
    make_normalized_Mx2_plots(NH3_FILES, C_FILES, H2_FILES, run_charges, NORMALIZED_OUTPUT)


if __name__ == "__main__":
    main()