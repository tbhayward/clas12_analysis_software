#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Input files and tree/branch names
    file_cj11 = "/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root"
    file_cj13 = "/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root"
    tree_name = "PhysicsEvents"
    branch = "Mh"

    # Histogram settings
    xmin, xmax = 0.2, 1.4
    nbins = 100
    bins = np.linspace(xmin, xmax, nbins + 1)

    # Load arrays
    with uproot.open(file_cj11) as f1:
        mh_cj11 = f1[tree_name].arrays(branch, library="np")[branch]
    #endif
    with uproot.open(file_cj13) as f2:
        mh_cj13 = f2[tree_name].arrays(branch, library="np")[branch]
    #endif

    # Optional: filter out NaNs/Infs
    mh_cj11 = mh_cj11[np.isfinite(mh_cj11)]
    mh_cj13 = mh_cj13[np.isfinite(mh_cj13)]

    # Plot
    plt.figure(figsize=(8, 6))
    plt.hist(mh_cj11, bins=bins, histtype="step", label="cj11", color="blue", linewidth=1.5, range=(xmin, xmax))
    plt.hist(mh_cj13, bins=bins, histtype="step", label="cj13", color="red", linewidth=1.5, range=(xmin, xmax))

    plt.xlim(xmin, xmax)
    plt.xlabel(r"$M_{h}$ (GeV)")
    plt.ylabel("Counts")
    plt.legend()

    plt.tight_layout()
    outpath = "/u/home/thayward/two_pion_mass.pdf"
    plt.savefig(outpath)
    print("Saved:", outpath)

if __name__ == "__main__":
    main()
# endif