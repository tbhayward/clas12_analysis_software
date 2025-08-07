#!/usr/bin/env python3
import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <base_file.root> <recalc_file.root>")
        sys.exit(1)

    base_file, recalc_file = sys.argv[1], sys.argv[2]

    # Open the ROOT files and extract the 't' and 'Mx2' branches
    with uproot.open(base_file) as f1:
        tree1    = f1["PhysicsEvents"]
        arr1     = tree1.arrays(["t", "Mx2"], library="np")
        t_base   = arr1["t"]
        mx2_base = arr1["Mx2"]

    with uproot.open(recalc_file) as f2:
        tree2      = f2["PhysicsEvents"]
        arr2       = tree2.arrays(["t", "Mx2"], library="np")
        t_recalc   = arr2["t"]
        mx2_recalc = arr2["Mx2"]

    # Apply Mx2 < 1.05 cut
    mask_base   = mx2_base   < 1.05
    mask_recalc = mx2_recalc < 1.05
    t_base_cut   = t_base[mask_base]
    t_recalc_cut = t_recalc[mask_recalc]

    # Define histogram bins from -10 to 1
    bins = np.linspace(-10, 1, 111)  # 0.1-wide bins

    # Plot
    plt.figure(figsize=(8, 6))
    plt.hist(t_base_cut,   bins=bins, histtype='step', label='base (Mx2 < 1.05)', color='blue')
    plt.hist(t_recalc_cut, bins=bins, histtype='step', label='recalculated (Mx2 < 1.05)', color='red')

    plt.xlabel(r'$t$')
    plt.ylabel('Counts')
    plt.legend(loc='upper right')
    plt.tight_layout()

    # Save
    outpath = "/u/home/thayward/t_example.pdf"
    plt.savefig(outpath)
    print(f"Saved histogram to {outpath}")

if __name__ == "__main__":
    main()