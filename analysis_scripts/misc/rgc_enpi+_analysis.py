#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# Toggle for quick debugging (process only first 5 runs per period)
QUICK_RUN = False

# Path to the CSV of run charges
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output directory and filenames
OUTPUT_DIR        = "output/enpi+"
NORMALIZED_OUTPUT = os.path.join(OUTPUT_DIR, "normalized_Mx2_charges.pdf")
PER_RUN_OUTPUT    = os.path.join(OUTPUT_DIR, "per_run_Fa22_Sp23_Mx2.pdf")

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
    Lines beginning with '#' are skipped.
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
                run    = int(parts[0])
                charge = float(parts[1])
            except:
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
    Applies mask Mx2 <= 2.0 and option to limit to first 5 runs if QUICK_RUN is True.
    """
    bins = np.linspace(0.0, 1.5, 101)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    def compute_charge_and_mask(tree, run_charges, quick=False):
        """
        Returns (Q_tot, masked_Mx2_array) where:
         - Q_tot = sum of charges for unique runs (first 5 if quick=True)
         - masked_Mx2_array = Mx2 values for those runs with Mx2 <= 2.0
        """
        run_arr = tree["runnum"].array(library="np").astype(int)
        mx2_arr = tree["Mx2"].array(library="np")
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

    # -------------------
    # LEFT SUBPLOT: NH3 + H2
    # -------------------
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

    # -------------------
    # RIGHT SUBPLOT: C + H2
    # -------------------
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


def make_per_run_Fa22_Sp23_plot(nh3_info, run_charges, outpath):
    """
    Build a 1×2 figure:
    - Left: Fa22-NH3, with one histogram per unique runnum
    - Right: Sp23-NH3, with one histogram per unique runnum
    Legends are placed in 3 columns at top-left, sorted numerically by run.
    Each line has a unique linestyle. Print any runnum whose histogram
    has a bin > 0.10 events / nC.
    """
    bins = np.linspace(0.0, 1.5, 101)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    # Linestyles to cycle through
    linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 2))]
    
    def per_run_data(tree, run_charges, quick=False):
        """
        Returns a dict { runnum: (Q_run, masked_Mx2_array) }.
        - Masks Mx2 <= 2.0
        - If quick=True, only first 5 unique runnums encountered.
        """
        run_arr = tree["runnum"].array(library="np").astype(int)
        mx2_arr = tree["Mx2"].array(library="np")
        if run_arr.size == 0:
            return {}

        # Collect unique runnums in order
        unique_runs = []
        for r in run_arr:
            if r not in unique_runs:
                unique_runs.append(r)
            if quick and len(unique_runs) >= 5:
                break

        data = {}
        for run in unique_runs:
            mask = (run_arr == run) & (mx2_arr <= 2.0)
            vals = mx2_arr[mask]
            Q_run = run_charges.get(run, 0.0)
            if Q_run > 0 and vals.size > 0:
                data[run] = (Q_run, vals)
        return data

    # Extract Fa22-NH3 and Sp23-NH3 entries
    fa22_entry = next((e for e in nh3_info if e['label'] == "Fa22-NH3"), None)
    sp23_entry = next((e for e in nh3_info if e['label'] == "Sp23-NH3"), None)

    # Fa22-NH3 per-run
    if fa22_entry:
        fa22_data = per_run_data(fa22_entry['tree'], run_charges, quick=QUICK_RUN)
        flagged_fa22 = []
        for idx, runnum in enumerate(sorted(fa22_data.keys())):
            Q_run, vals = fa22_data[runnum]
            counts, _ = np.histogram(vals, bins=bins)
            counts = counts.astype(float) / Q_run
            if np.any(counts > 0.10):
                flagged_fa22.append(runnum)
            style = linestyles[idx % len(linestyles)]
            axes[0].step(bins[:-1], counts, where='post', linestyle=style, label=str(runnum))
        if flagged_fa22:
            print("Fa22-NH3 runs with bin >0.10 events/nC:", sorted(flagged_fa22))
    else:
        print("[WARNING] Fa22-NH3 entry not found!")

    axes[0].set(title="Fa22-NH3 per run",
                xlabel=r"$M_x^2$ [GeV$^2$]",
                ylabel="events / nC",
                xlim=(0.0, 1.5))
    axes[0].legend(loc='upper left', fontsize='xx-small', ncol=3)

    # Sp23-NH3 per-run
    if sp23_entry:
        sp23_data = per_run_data(sp23_entry['tree'], run_charges, quick=QUICK_RUN)
        flagged_sp23 = []
        for idx, runnum in enumerate(sorted(sp23_data.keys())):
            Q_run, vals = sp23_data[runnum]
            counts, _ = np.histogram(vals, bins=bins)
            counts = counts.astype(float) / Q_run
            if np.any(counts > 0.10):
                flagged_sp23.append(runnum)
            style = linestyles[idx % len(linestyles)]
            axes[1].step(bins[:-1], counts, where='post', linestyle=style, label=str(runnum))
        if flagged_sp23:
            print("Sp23-NH3 runs with bin >0.10 events/nC:", sorted(flagged_sp23))
    else:
        print("[WARNING] Sp23-NH3 entry not found!")

    axes[1].set(title="Sp23-NH3 per run", xlabel=r"$M_x^2$ [GeV$^2$]", xlim=(0.0, 1.5))
    axes[1].legend(loc='upper left', fontsize='xx-small', ncol=3)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()
#enddef


def main():
    # Load run charges
    run_charges = parse_run_charges(CSV_PATH)

    # Load trees
    nh3_info = load_trees(NH3_FILES)
    c_info   = load_trees(C_FILES)
    h2_info  = load_trees(H2_FILES)

    # First figure: normalized Mx2 for NH3, C, and H2
    make_normalized_Mx2_plots(nh3_info, c_info, h2_info, run_charges, NORMALIZED_OUTPUT)

    # Second figure: per-run Mx2 for Fa22-NH3 and Sp23-NH3
    make_per_run_Fa22_Sp23_plot(nh3_info, run_charges, PER_RUN_OUTPUT)
#enddef


if __name__ == '__main__':
    main()
#endif