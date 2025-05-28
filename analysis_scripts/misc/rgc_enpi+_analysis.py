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
    print(f"[DEBUG] parse_run_charges: opening '{csv_path}'")
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
            if i <= 5:
                print(f"  parsed: run={run}, charge={charge}")
    print(f"[DEBUG] parsed total runs: {len(run_charges)}")
    return run_charges


def load_trees(files):
    info = []
    for fp, label in files:
        if not os.path.exists(fp):
            print(f"[ERROR] ROOT file missing: {fp}")
            continue
        tree = uproot.open(fp)["PhysicsEvents"]
        print(f"[DEBUG] loaded tree for {label}: {tree} ")
        info.append({'tree': tree, 'label': label})
    return info


def make_normalized_Mx2_plots(nh3_info, c_info, run_charges, outpath):
    bins = np.linspace(-1.0, 2.0, 101)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    def inspect_branches(tree, label):
        keys = tree.keys()
        print(f"[DEBUG] branches in {label}: {keys}")
        return keys

    def compute_charge(tree, label):
        run_arr = tree.arrays(tracks=["runnum"], library="np")["runnum"].astype(int)
        if run_arr.size == 0:
            print(f"[WARNING] runnum array empty for {label}")
            return 0
        uni = np.unique(run_arr)
        Q = sum(run_charges.get(r, 0.0) for r in uni)
        print(f"    {label}: unique runs={len(uni)}, Q={Q:.3f}")
        return Q

    print("[DEBUG] NH3 processing:")
    for entry in nh3_info:
        t, lbl = entry['tree'], entry['label']
        inspect_branches(t, lbl)
        Qtot = compute_charge(t, lbl)
        if Qtot <= 0:
            continue
        arr = t.arrays(["Mx2"], library="np")["Mx2"]
        cts, _ = np.histogram(arr, bins=bins)
        axes[0].step(bins[:-1], cts.astype(float)/Qtot, where='post', label=lbl)
    axes[0].set(title="NH₃", xlabel=r"$M_x^2$ [GeV$^2$]", ylabel="events/nC")
    axes[0].legend(loc='upper right')

    print("[DEBUG] C processing:")
    for entry in c_info:
        t, lbl = entry['tree'], entry['label']
        inspect_branches(t, lbl)
        Qtot = compute_charge(t, lbl)
        if Qtot <= 0:
            continue
        arr = t.arrays(["Mx2"], library="np")["Mx2"]
        cts, _ = np.histogram(arr, bins=bins)
        axes[1].step(bins[:-1], cts.astype(float)/Qtot, where='post', label=lbl)
    axes[1].set(title="C", xlabel=r"$M_x^2$ [GeV$^2$")
    axes[1].legend(loc='upper right')

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def main():
    run_charges = parse_run_charges(CSV_PATH)
    nh3 = load_trees(NH3_FILES)
    c   = load_trees(C_FILES)
    make_normalized_Mx2_plots(nh3, c, run_charges, OUTPUT_FILE)

if __name__ == '__main__':
    main()