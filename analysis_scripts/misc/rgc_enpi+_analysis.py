#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# Toggle for quick debugging (process only first 5 runs per period)
QUICK_RUN = False

# Maximum run number to include
MAX_RUNNUM = 17768

# Path to the CSV of run charges
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output directory and filenames
OUTPUT_DIR               = "output/enpi+"
NORMALIZED_OUTPUT        = os.path.join(OUTPUT_DIR, "normalized_Mx2_charges.pdf")
RUNBYRUN_SP23C_OUTPUT    = os.path.join(OUTPUT_DIR, "runbyrun_Sp23C_Mx2.pdf")

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


def parse_run_charges(csv_path):
    """
    Read the CSV and return a dict { runnum: charge }.
    Skip any run > MAX_RUNNUM.
    Lines beginning with '#' are ignored.
    """
    run_charges = {}
    if not os.path.exists(csv_path):
        print(f"[ERROR] CSV not found: {csv_path}")
        return run_charges

    with open(csv_path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split(',')
            if len(parts) < 2:
                continue
            try:
                run = int(parts[0])
                charge = float(parts[1])
            except:
                continue
            if run > MAX_RUNNUM:
                continue
            run_charges[run] = charge
    return run_charges


def load_trees(file_list):
    """
    Given list of (filepath, label), return list of dicts each containing
    the opened uproot tree plus its label.
    """
    info = []
    for fp, label in file_list:
        if not os.path.exists(fp):
            print(f"[ERROR] ROOT file missing: {fp}")
            continue
        tree = uproot.open(fp)["PhysicsEvents"]
        info.append({'tree': tree, 'label': label})
    return info


def make_normalized_Mx2_plots(nh3_info, c_info, h2_info, run_charges, outpath):
    """
    Build a 1×2 figure:
    - Left: NH3 periods (Su22, Fa22, Sp23) + H2
    - Right: C periods   (Su22, Fa22, Sp23) + H2
    Exclude runs > MAX_RUNNUM. Apply mask Mx2 <= 2.0.
    If QUICK_RUN is True, only first 5 unique runnums per tree are used.
    """
    bins = np.linspace(0.0, 1.5, 101)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    def compute_charge_and_mask(tree, run_charges, quick=False):
        """
        Returns (Q_tot, masked_Mx2_array) where:
         - Q_tot = sum of charges for unique runs ≤ MAX_RUNNUM
                   (only first 5 if quick=True)
         - masked_Mx2_array = Mx2 values for those runs with Mx2 <= 2.0
        """
        run_arr = tree["runnum"].array(library="np").astype(int)
        mx2_arr = tree["Mx2"].array(library="np")
        if run_arr.size == 0:
            return 0.0, np.array([], dtype=float)

        # Filter out runs > MAX_RUNNUM
        valid_mask = (run_arr <= MAX_RUNNUM)
        run_arr = run_arr[valid_mask]
        mx2_arr = mx2_arr[valid_mask]
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
        mask = np.isin(run_arr, list(selected_runs)) & (mx2_arr <= 2.0)
        return Q_tot, mx2_arr[mask]

    # Left subplot: NH3 + H2
    for entry in nh3_info:
        tree, lbl = entry['tree'], entry['label']
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        counts = counts.astype(float) / Q
        axes[0].step(bins[:-1], counts, where='post', label=lbl)
    for entry in h2_info:
        tree, lbl = entry['tree'], entry['label']
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        counts = counts.astype(float) / Q
        axes[0].step(bins[:-1], counts, where='post', color='k', linestyle='-', label=lbl)
    axes[0].set(xlabel=r"$M_x^2$ [GeV$^2$]", ylabel="events / nC", xlim=(0.0, 1.5))
    axes[0].legend(loc='upper right', fontsize='small')

    # Right subplot: C + H2
    for entry in c_info:
        tree, lbl = entry['tree'], entry['label']
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        counts = counts.astype(float) / Q
        axes[1].step(bins[:-1], counts, where='post', label=lbl)
    for entry in h2_info:
        tree, lbl = entry['tree'], entry['label']
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        counts = counts.astype(float) / Q
        axes[1].step(bins[:-1], counts, where='post', color='k', linestyle='-', label=lbl)
    axes[1].set(xlabel=r"$M_x^2$ [GeV$^2$]", xlim=(0.0, 1.5))
    axes[1].legend(loc='upper right', fontsize='small')

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()
#enddef


def make_runbyrun_Sp23C_plot(c_info, run_charges, outpath):
    """
    Build a run-by-run histogram plot for the Spring 2023 Carbon dataset (Sp23-C):
    - One panel only, with one histogram per unique runnum.
    - Exclude runs > MAX_RUNNUM, apply Mx2 <= 2.0 mask.
    - Each run line has a unique linestyle. Legend at top-left, 3 columns, small font.
    """
    bins = np.linspace(0.0, 1.5, 101)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Predefined linestyles to cycle through
    linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 2)), (0, (3, 5, 1, 5))]

    # Find the Sp23-C entry
    sp23_entry = next((e for e in c_info if e['label'] == "Sp23-C"), None)
    if sp23_entry is None:
        print("[ERROR] Sp23-C not found in C_FILES.")
        return

    tree = sp23_entry['tree']

    # Extract runnum and Mx2 arrays
    run_arr = tree["runnum"].array(library="np").astype(int)
    mx2_arr = tree["Mx2"].array(library="np")

    # Filter by MAX_RUNNUM and Mx2 <= 2.0
    valid_mask = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
    run_arr = run_arr[valid_mask]
    mx2_arr = mx2_arr[valid_mask]
    if run_arr.size == 0:
        print("[WARNING] No valid events for Sp23-C after filtering.")
        return

    # Unique runs in ascending order
    unique_runs = sorted(np.unique(run_arr))

    # For each run, build histogram normalized by that run's charge
    for idx, runnum in enumerate(unique_runs):
        run_mask = run_arr == runnum
        vals = mx2_arr[run_mask]
        Q_run = run_charges.get(runnum, 0.0)
        if Q_run <= 0 or vals.size == 0:
            continue
        counts, _ = np.histogram(vals, bins=bins)
        counts = counts.astype(float) / Q_run
        style = linestyles[idx % len(linestyles)]
        ax.step(bins[:-1], counts, where='post', linestyle=style, label=str(runnum))

    ax.set(title="Sp23-C per run", xlabel=r"$M_x^2$ [GeV$^2$]", ylabel="events / nC", xlim=(0.0, 1.5))
    ax.legend(loc='upper left', fontsize='xx-small', ncol=3)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()
#enddef


def main():
    # Load run charges (excluding run > MAX_RUNNUM)
    run_charges = parse_run_charges(CSV_PATH)

    # Load trees
    nh3_info = load_trees(NH3_FILES)
    c_info   = load_trees(C_FILES)
    h2_info  = load_trees(H2_FILES)

    # Generate normalized Mx2 histograms
    make_normalized_Mx2_plots(nh3_info, c_info, h2_info, run_charges, NORMALIZED_OUTPUT)

    # Generate run-by-run plot for Sp23-C
    make_runbyrun_Sp23C_plot(c_info, run_charges, RUNBYRUN_SP23C_OUTPUT)
#enddef


if __name__ == '__main__':
    main()
#endif