#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    # file paths and labels
    files = [
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "sidisdvcs_rgc_su22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "sidisdvcs_rgc_fa22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "sidisdvcs_rgc_sp23_inb_calibration.root",
    ]
    labels = ["Su22", "Fa22", "Sp23"]
    tree_name = "PhysicsEvents"

    # vertex cuts for negative tracks
    vz_cuts = {
        "Su22": (-7.576, 0.303),
        "Fa22": (-5.758, 1.515),
        "Sp23": (-5.758, 1.515),
    }

    # bins for nphe
    bins = np.linspace(0, 40, 100)

    # collect per-period cc_nphe arrays
    nphe_data = []

    for fname, label in zip(files, labels):
        tree = uproot.open(fname)[tree_name]
        pid       = tree["particle_pid"].array(library="np")
        vz        = tree["particle_vz"].array(library="np")
        p         = tree["p"].array(library="np")
        nphe      = tree["cc_nphe_15"].array(library="np")
        sector6   = tree["track_sector_6"].array(library="np")

        # mask negative tracks, vertex, momentum, valid sector, valid nphe
        mask = (
            ((pid == 11)  | (pid == -211) | (pid == -321)) &
            (sector6 != -9999) &
            (vz >= vz_cuts[label][0]) &
            (vz <= vz_cuts[label][1]) &
            (p > 2.0) &
            (nphe != -9999)
        )
        nphe_data.append(nphe[mask])
    # endfor

    # prepare figure
    fig, axes = plt.subplots(1, 2, figsize=(12,6))
    colors = ["C0","C1","C2"]

    # plot histograms
    for data, label, color in zip(nphe_data, labels, colors):
        axes[0].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)
        axes[1].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)
    # endfor

    # add vertical line at nphe = 2
    for ax in axes:
        ax.axvline(2, color='red', linestyle='-', linewidth=2)

    # left panel: linear
    axes[0].set_xlabel("HTCC nphe")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Negative Tracks: HTCC nphe")
    axes[0].legend()

    # right panel: log
    axes[1].set_yscale("log")
    axes[1].set_xlabel("HTCC nphe")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Negative Tracks: HTCC nphe (log scale)")
    axes[1].legend()

    fig.tight_layout()
    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(f"{outdir}/electron_htcc.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()