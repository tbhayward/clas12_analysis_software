#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    # Paths to your ROOT files
    files = [
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "sidisdvcs_rgc_su22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "sidisdvcs_rgc_fa22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "sidisdvcs_rgc_sp23_inb_calibration.root",
    ]
    labels = ["Su22", "Fa22", "Sp23"]
    colors = ["C0", "C1", "C2"]
    tree_name = "PhysicsEvents"

    # vertex‐cut thresholds you determined:
    vz_cuts = {
        "Su22":    (-7.576,  0.303),
        "Fa22":    (-5.758,  1.515),
        "Sp23":    (-5.758,  1.515),  # same as Fa22
    }

    # Collect PCal energy arrays, one per run period
    pcal_energies = []

    for fname, label in zip(files, labels):
        with uproot.open(fname)[tree_name] as tree:
            # read branches
            pid    = tree["particle_pid"].array(library="np")
            vz     = tree["particle_vz"].array(library="np")
            sector = tree["track_sector_6"].array(library="np")
            e_pcal = tree["cal_energy_1"].array(library="np")

            # apply selections
            mask = (
                (pid == 11) &
                (sector != -9999) &
                (e_pcal >= 0) &
                (vz >= vz_cuts[label][0]) &
                (vz <= vz_cuts[label][1])
            )
            pcal_energies.append(e_pcal[mask])
    #endfor

    # histogram settings: 100 bins from 0 to 1.5 GeV
    bins = np.linspace(0, 1.5, 100)

    # prepare output directory
    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    # create 1×2 figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for data, label, color in zip(pcal_energies, labels, colors):
        axes[0].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)
        axes[1].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)
    #endfor

    # left panel: linear scale
    axes[0].set_xlabel(r"$E_{\mathrm{PCal}}$ (GeV)")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Electron PCal Energy Deposition")
    axes[0].legend()

    # right panel: log scale
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$E_{\mathrm{PCal}}$ (GeV)")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Electron PCal Energy Deposition (log scale)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{outdir}/electron_pcal_deposition.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()