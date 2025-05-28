#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# Path to the CSV of run charges
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
    """
    run_charges = {}
    with open(csv_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            run    = int(parts[0])
            charge = float(parts[1])
            run_charges[run] = charge
    #endfor
    return run_charges


def load_trees(file_list):
    """
    Given list of (filepath, label), return list of dicts each containing
    the opened uproot tree plus its label.
    """
    info = []
    for fp, label in file_list:
        tree = uproot.open(fp)["PhysicsEvents"]
        info.append({'tree': tree, 'label': label})
    #endfor
    return info


def make_normalized_Mx2_plots(nh3_info, c_info, run_charges, output_path):
    """
    Build and save a 1×2 figure of normalized Mx^2 histograms for NH3 and C.
    """
    # Histogram settings
    n_bins      = 100
    x_min, x_max = -1.0, 2.0
    bins        = np.linspace(x_min, x_max, n_bins + 1)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    # Left: NH3
    print("[DEBUG] NH3 total charges:")
    for entry in nh3_info:
        tree  = entry['tree']
        label = entry['label']

        runnums = tree.arrays("runnum", library="np")["runnum"]
        uniques = np.unique(runnums)
        Q_tot   = sum(run_charges.get(r, 0.0) for r in uniques)
        print(f"  {label}: Q_tot = {Q_tot:.3f} nC")

        if Q_tot <= 0:
            print(f"[WARNING] Skipping {label} (Q_tot={Q_tot:.3f})")
            continue
        #endif

        Mx2_arr = tree.arrays("Mx2", library="np")["Mx2"]
        counts, _ = np.histogram(Mx2_arr, bins=bins)
        counts = counts.astype(float)      # ensure float dtype
        counts /= Q_tot                    # normalize by total charge

        axes[0].step(bins[:-1], counts, where="post", label=label)
    #endfor
    axes[0].set_title("NH₃")
    axes[0].set_xlabel(r"$M_x^2$ [GeV$^2$")
    axes[0].set_ylabel("events / nC")
    axes[0].set_xlim(x_min, x_max)
    axes[0].legend(loc="upper right")

    # Right: C
    print("[DEBUG] C total charges:")
    for entry in c_info:
        tree  = entry['tree']
        label = entry['label']

        runnums = tree.arrays("runnum", library="np")["runnum"]
        uniques = np.unique(runnums)
        Q_tot   = sum(run_charges.get(r, 0.0) for r in uniques)
        print(f"  {label}: Q_tot = {Q_tot:.3f} nC")

        if Q_tot <= 0:
            print(f"[WARNING] Skipping {label} (Q_tot={Q_tot:.3f})")
            continue
        #endif

        Mx2_arr = tree.arrays("Mx2", library="np")["Mx2"]
        counts, _ = np.histogram(Mx2_arr, bins=bins)
        counts = counts.astype(float)
        counts /= Q_tot

        axes[1].step(bins[:-1], counts, where="post", label=label)
    #endfor
    axes[1].set_title("C")
    axes[1].set_xlabel(r"$M_x^2$ [GeV$^2$")
    axes[1].set_xlim(x_min, x_max)
    axes[1].legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


def main():
    # 1) Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 2) Load run → charge map
    run_charges = parse_run_charges(CSV_PATH)

    # 3) Open all the trees
    nh3_info = load_trees(NH3_FILES)
    c_info   = load_trees(C_FILES)
    h2_info  = load_trees(H2_FILES)  # for future hydrogen plots
    d2_info  = load_trees(D2_FILES)  # for future deuterium plots

    # 4) Generate & save the Mx2 plots
    make_normalized_Mx2_plots(nh3_info, c_info, run_charges, OUTPUT_FILE)
#endfor

if __name__ == "__main__":
    main()
#endif