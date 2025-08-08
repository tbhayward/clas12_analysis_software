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

    # Load arrays
    with uproot.open(file_cj11) as f1:
        mh_cj11 = f1[tree_name].arrays(branch, library="np")[branch]
    with uproot.open(file_cj13) as f2:
        mh_cj13 = f2[tree_name].arrays(branch, library="np")[branch]

    # Optional: filter out NaNs/Infs
    mh_cj11 = mh_cj11[np.isfinite(mh_cj11)]
    mh_cj13 = mh_cj13[np.isfinite(mh_cj13)]

    # Create figure with 1x2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ---------- Left Panel: Full range ----------
    xmin_full, xmax_full = 0.2, 1.4
    nbins_full = 200  # doubled from previous 100
    bins_full = np.linspace(xmin_full, xmax_full, nbins_full + 1)

    axes[0].hist(mh_cj11, bins=bins_full, histtype="step", label="cj11", color="blue", linewidth=1.5)
    axes[0].hist(mh_cj13, bins=bins_full, histtype="step", label="cj13", color="red", linewidth=1.5)
    axes[0].set_xlim(xmin_full, xmax_full)
    axes[0].set_xlabel(r"$M_{h}$ (GeV)")
    axes[0].set_ylabel("Counts")
    axes[0].legend()

    # ---------- Right Panel: Zoomed range ----------
    xmin_zoom, xmax_zoom = 0.3, 0.6
    nbins_zoom = 60  # finer binning for zoom
    bins_zoom = np.linspace(xmin_zoom, xmax_zoom, nbins_zoom + 1)

    axes[1].hist(mh_cj11, bins=bins_zoom, histtype="step", label="cj11", color="blue", linewidth=1.5)
    axes[1].hist(mh_cj13, bins=bins_zoom, histtype="step", label="cj13", color="red", linewidth=1.5)
    axes[1].set_xlim(xmin_zoom, xmax_zoom)
    axes[1].set_xlabel(r"$M_{h}$ (GeV)")
    axes[1].set_ylabel("Counts")
    axes[1].legend()

    # ---------- Save ----------
    plt.tight_layout()
    outpath = "/u/home/thayward/two_pion_mass.pdf"
    plt.savefig(outpath)
    print("Saved:", outpath)

if __name__ == "__main__":
    main()
# endif