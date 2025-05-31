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
OUTPUT_DIR        = "output/enpi+"
NORMALIZED_OUTPUT = os.path.join(OUTPUT_DIR, "normalized_Mx2_charges.pdf")

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


def print_unique_carbon_runs(c_info):
    """
    Print every unique run number found in the carbon trees (≤ MAX_RUNNUM).
    """
    for entry in c_info:
        tree, lbl = entry['tree'], entry['label']
        run_arr = tree["runnum"].array(library="np").astype(int)
        # Filter by MAX_RUNNUM
        run_arr = run_arr[run_arr <= MAX_RUNNUM]
        unique_runs = sorted(np.unique(run_arr))
        print(f"{lbl} unique runs (≤ {MAX_RUNNUM}): {unique_runs}")
#enddef


def main():
    # Load run charges (skipping runs > MAX_RUNNUM)
    run_charges = parse_run_charges(CSV_PATH)

    # Load trees
    nh3_info = load_trees(NH3_FILES)
    c_info   = load_trees(C_FILES)
    h2_info  = load_trees(H2_FILES)

    # Print every unique carbon run number (≤ MAX_RUNNUM)
    print_unique_carbon_runs(c_info)

    # Generate normalized Mx2 histograms, skipping runs > MAX_RUNNUM
    make_normalized_Mx2_plots(nh3_info, c_info, h2_info, run_charges, NORMALIZED_OUTPUT)
#enddef


if __name__ == '__main__':
    main()
#endif