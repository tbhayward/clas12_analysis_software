#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# Path to the CSV of run charges (update this if the file is elsewhere)
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output settings
OUTPUT_DIR  = "output/enpi+"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "normalized_Mx2_charges.pdf")

# --- File lists: (filepath, label) ---
NH3_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_NH3_epi+.root", "Su22"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+.root", "Fa22"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_NH3_epi+.root", "Sp23"),
]

C_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_C_epi+.root",  "Su22"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_C_epi+.root",  "Fa22"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_C_epi+.root",  "Sp23"),
]

H2_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_sp19_inb_H2_epi+.root", "Sp19"),
]

D2_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgb_sp19_inb_D2_epi+.root", "Sp19"),
]

def parse_run_charges(csv_path):
    """
    Read the CSV and return a dict { runnum: charge }.
    Lines beginning with '#' are skipped.
    Includes debug prints to verify file loading.
    """
    print(f"[DEBUG] parse_run_charges: opening CSV file '{csv_path}'")
    if not os.path.exists(csv_path):
        print(f"[ERROR] CSV file not found at '{csv_path}'. Please check CSV_PATH.")
        return {}
    run_charges = {}
    line_count = 0
    with open(csv_path, 'r') as f:
        for line in f:
            line_count += 1
            raw = line.strip()
            if not raw or raw.startswith('#'):
                continue
            parts = raw.split(',')
            try:
                run    = int(parts[0])
                charge = float(parts[1])
            except Exception as e:
                print(f"[WARNING] failed to parse line {line_count}: '{raw}' -> {e}")
                continue
            run_charges[run] = charge
            if len(run_charges) <= 5:
                print(f"  [DEBUG] parsed run {run} -> charge {charge}")
        #endfor
    print(f"[DEBUG] total runs parsed: {len(run_charges)}")
    return run_charges


def load_trees(file_list):
    info = []
    for fp, label in file_list:
        if not os.path.exists(fp):
            print(f"[ERROR] ROOT file not found: {fp}")
            continue
        tree = uproot.open(fp)["PhysicsEvents"]
        info.append({'tree': tree, 'label': label})
    #endfor
    return info


def make_normalized_Mx2_plots(nh3_info, c_info, run_charges, output_path):
    n_bins      = 100
    x_min, x_max = -1.0, 2.0
    bins        = np.linspace(x_min, x_max, n_bins + 1)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    def compute_total_charge(tree, label):
        runnum_arr = tree.arrays("runnum", library="np")["runnum"].astype(int)
        uniques     = np.unique(runnum_arr)
        Q_tot       = sum(run_charges.get(r, 0.0) for r in uniques)
        print(f"    [DEBUG] {label}: unique runs={len(uniques)}, Q_tot={Q_tot:.3f}")
        return Q_tot

    print("[DEBUG] NH3 total charges:")
    for entry in nh3_info:
        label = entry['label']
        tree  = entry['tree']
        Q_tot = compute_total_charge(tree, label)
        if Q_tot <= 0:
            print(f"    [WARNING] skipping {label} (Q_tot={Q_tot:.3f})")
            continue
        Mx2_arr = tree.arrays("Mx2", library="np")["Mx2"]
        counts, _ = np.histogram(Mx2_arr, bins=bins)
        counts = counts.astype(float)
        counts /= Q_tot
        axes[0].step(bins[:-1], counts, where="post", label=label)
    #endfor
    axes[0].set(title="NH₃", xlabel=r"$M_x^2$ [GeV$^2$]", ylabel="events / nC", xlim=(x_min, x_max))
    axes[0].legend(loc="upper right")

    print("[DEBUG] C total charges:")
    for entry in c_info:
        label = entry['label']
        tree  = entry['tree']
        Q_tot = compute_total_charge(tree, label)
        if Q_tot <= 0:
            print(f"    [WARNING] skipping {label} (Q_tot={Q_tot:.3f})")
            continue
        Mx2_arr = tree.arrays("Mx2", library="np")["Mx2"]
        counts, _ = np.histogram(Mx2_arr, bins=bins)
        counts = counts.astype(float)
        counts /= Q_tot
        axes[1].step(bins[:-1], counts, where="post", label=label)
    #endfor
    axes[1].set(title="C", xlabel=r"$M_x^2$ [GeV$^2$]", xlim=(x_min, x_max))
    axes[1].legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    run_charges = parse_run_charges(CSV_PATH)
    nh3_info = load_trees(NH3_FILES)
    c_info   = load_trees(C_FILES)
    # for future use: h2_info = load_trees(H2_FILES); d2_info = load_trees(D2_FILES)
    make_normalized_Mx2_plots(nh3_info, c_info, run_charges, OUTPUT_FILE)
#endif

if __name__ == "__main__":
    main()
#endif