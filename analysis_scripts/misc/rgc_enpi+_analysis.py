#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# Path to the CSV of run charges (update if needed)
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output settings
OUTPUT_DIR  = "output/enpi+"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "normalized_Mx2_charges.pdf")

# --- File lists: (filepath, label) ---
NH3_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_NH3_epi+.root", "Su22-NH3"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+.root", "Fa22-NH3"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_NH3_epi+.root", "Sp23-NH3"),
]
C_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_C_epi+.root",  "Su22-C"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_C_epi+.root",  "Fa22-C"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_C_epi+.root",  "Sp23-C"),
]
# Add "Sp19" to H2 and D2 labels
H2_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_sp19_inb_H2_epi+.root", "Sp19-H2"),
]
D2_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgb_sp19_inb_D2_epi+.root", "Sp19-D2"),
]

def parse_run_charges(csv_path):
    """
    Read the CSV and return a dict { runnum: charge }.
    Lines beginning with '#' are skipped.
    """
    if not os.path.exists(csv_path):
        print(f"[ERROR] CSV not found: {csv_path}")
        return {}
    run_charges = {}
    with open(csv_path, 'r') as f:
        for i, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split(',')
            if len(parts) < 2:
                print(f"[WARNING] line {i} malformed: '{s}'")
                continue
            try:
                run    = int(parts[0])
                charge = float(parts[1])
            except Exception as e:
                print(f"[WARNING] parsing error on line {i}: {e}")
                continue
            run_charges[run] = charge
    return run_charges

def load_trees(files):
    """
    Given list of (filepath, label), return list of dicts each containing
    the opened uproot tree plus its label.
    """
    info = []
    for fp, label in files:
        if not os.path.exists(fp):
            print(f"[ERROR] ROOT file missing: {fp}")
            continue
        tree = uproot.open(fp)["PhysicsEvents"]
        info.append({'tree': tree, 'label': label})
    return info

def make_normalized_Mx2_plots(nh3_info, c_info, h2_info, d2_info, run_charges, outpath):
    """
    Build and save a 1×2 figure of normalized Mx^2 histograms:
    - Left: NH3 (3 periods) + H2 + D2
    - Right: C    (3 periods) + H2 + D2
    """
    bins = np.linspace(0.0, 1.5, 101)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    def compute_charge(tree, label):
        """
        Sum total charge by unique runnum values in the tree.
        """
        run_arr = tree["runnum"].array(library="np").astype(int)
        if run_arr.size == 0:
            print(f"[WARNING] runnum array empty for {label}")
            return 0.0
        uniques = np.unique(run_arr)
        Q = sum(run_charges.get(r, 0.0) for r in uniques)
        return Q

    # Plot NH3 (left subplot)
    for entry in nh3_info:
        tree, lbl = entry['tree'], entry['label']
        Qtot = compute_charge(tree, lbl)
        if Qtot <= 0:
            continue
        arr = tree["Mx2"].array(library="np")
        counts, _ = np.histogram(arr, bins=bins)
        counts = counts.astype(float) / Qtot
        axes[0].step(bins[:-1], counts, where='post', label=lbl)
    #endfor

    # Plot H2 on left
    for entry in h2_info:
        tree, lbl = entry['tree'], entry['label']
        Qtot = compute_charge(tree, lbl)
        if Qtot <= 0:
            continue
        arr = tree["Mx2"].array(library="np")
        counts, _ = np.histogram(arr, bins=bins)
        counts = counts.astype(float) / Qtot
        axes[0].step(bins[:-1], counts, where='post', color='k', linestyle='-', label=lbl)
    #endfor

    # Plot D2 on left
    for entry in d2_info:
        tree, lbl = entry['tree'], entry['label']
        Qtot = compute_charge(tree, lbl)
        if Qtot <= 0:
            continue
        arr = tree["Mx2"].array(library="np")
        counts, _ = np.histogram(arr, bins=bins)
        counts = counts.astype(float) / Qtot
        axes[0].step(bins[:-1], counts, where='post', color='k', linestyle='--', label=lbl)
    #endfor

    axes[0].set(xlabel=r"$M_x^2$ [GeV$^2$]", ylabel="events / nC", xlim=(0.0, 1.5))
    axes[0].legend(loc='upper right')

    # Plot C (right subplot)
    for entry in c_info:
        tree, lbl = entry['tree'], entry['label']
        Qtot = compute_charge(tree, lbl)
        if Qtot <= 0:
            continue
        arr = tree["Mx2"].array(library="np")
        counts, _ = np.histogram(arr, bins=bins)
        counts = counts.astype(float) / Qtot
        axes[1].step(bins[:-1], counts, where='post', label=lbl)
    #endfor

    # Plot H2 on right
    for entry in h2_info:
        tree, lbl = entry['tree'], entry['label']
        Qtot = compute_charge(tree, lbl)
        if Qtot <= 0:
            continue
        arr = tree["Mx2"].array(library="np")
        counts, _ = np.histogram(arr, bins=bins)
        counts = counts.astype(float) / Qtot
        axes[1].step(bins[:-1], counts, where='post', color='k', linestyle='-', label=lbl)
    #endfor

    # Plot D2 on right
    for entry in d2_info:
        tree, lbl = entry['tree'], entry['label']
        Qtot = compute_charge(tree, lbl)
        if Qtot <= 0:
            continue
        arr = tree["Mx2"].array(library="np")
        counts, _ = np.histogram(arr, bins=bins)
        counts = counts.astype(float) / Qtot
        axes[1].step(bins[:-1], counts, where='post', color='k', linestyle='--', label=lbl)
    #endfor

    axes[1].set(xlabel=r"$M_x^2$ [GeV$^2$]", xlim=(0.0, 1.5))
    axes[1].legend(loc='upper right')

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()
#enddef

def main():
    run_charges = parse_run_charges(CSV_PATH)
    nh3 = load_trees(NH3_FILES)
    c   = load_trees(C_FILES)
    h2  = load_trees(H2_FILES)
    d2  = load_trees(D2_FILES)
    make_normalized_Mx2_plots(nh3, c, h2, d2, run_charges, OUTPUT_FILE)
#enddef

if __name__ == '__main__':
    main()
#endif