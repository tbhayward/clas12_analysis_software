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
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_fa22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_sp23_inb_calibration.root",
    ]
    labels = ["Su22", "Fa22", "Sp23"]
    colors = ["C0", "C1", "C2"]
    tree_name = "PhysicsEvents"

    # Collect vz arrays
    electron_vz = []
    proton_vz   = []

    for fname in files:
        with uproot.open(fname)[tree_name] as tree:
            pid    = tree["particle_pid"].array(library="np")
            vz     = tree["particle_vz"].array(library="np")
            sector = tree["track_sector_6"].array(library="np")

            mask_e = (pid == 11)   & (sector != -9999)
            mask_p = (pid == 2212) & (sector != -9999)

            electron_vz.append(vz[mask_e])
            proton_vz.append(vz[mask_p])

    # histogram bin edges from -15 to 15
    bins = np.linspace(-15, 15, 100)

    # ensure output dir exists
    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    # --- Electron vertex distribution ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for data, label, color in zip(electron_vz, labels, colors):
        axes[0].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)
        axes[1].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)

    # draw electron vertex cuts with 50% opacity
    for ax in axes:
        ax.axvline(-7, color='red', linestyle='-', alpha=0.3)
        ax.axvline(-0.5,  color='red', linestyle='-', alpha=0.3)
        ax.axvline(-6, color='red', linestyle='--', alpha=0.3)
        ax.axvline(0.5,  color='red', linestyle='--', alpha=0.3)
        ax.set_xlim(-15, 15)

    # annotate panels
    axes[0].set_xlabel(r"$v_{z}$ (cm)")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Electron Vertex Distribution")
    axes[0].legend()

    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$v_{z}$ (cm)")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Electron Vertex Distribution (log scale)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{outdir}/electron_vz.pdf")
    plt.close(fig)

    # --- Proton vertex distribution ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for data, label, color in zip(proton_vz, labels, colors):
        axes[0].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)
        axes[1].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)

    # draw proton vertex cuts with 50% opacity
    for ax in axes:
        ax.axvline(-8,   color='red', linestyle='-', alpha=0.3)
        ax.axvline(-0.5,    color='red', linestyle='-', alpha=0.3)
        ax.axvline(-7.5, color='red', linestyle='--', alpha=0.3)
        ax.axvline(0.0,  color='red', linestyle='--', alpha=0.3)
        ax.set_xlim(-15, 15)

    # annotate panels
    axes[0].set_xlabel(r"$v_{z}$ (cm)")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Proton Vertex Distribution")
    axes[0].legend()

    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$v_{z}$ (cm)")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Proton Vertex Distribution (log scale)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{outdir}/proton_vz.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()