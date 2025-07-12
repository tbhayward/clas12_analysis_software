#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main():
    # Paths to your ROOT files
    files = [
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_fa22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_sp23_inb_calibration.root",
    ]
    labels = ["Su22", "Fa22", "Sp23"]
    colors = ["C0", "C1", "C2"]

    # Name of the TTree inside each file (adjust if needed)
    tree_name = "T"

    # Containers for vertex z data
    electron_vz = []
    proton_vz = []

    # Loop over files and extract arrays
    for fname in files:
        with uproot.open(fname)[tree_name] as tree:
            pid    = tree["particle_pid"].array(library="np")
            vz     = tree["particle_vz"].array(library="np")
            sector = tree["track_sector_6"].array(library="np")

            mask_e = (pid == 11) & (sector != -9999)
            electron_vz.append(vz[mask_e])

            mask_p = (pid == 2212) & (sector != -9999)
            proton_vz.append(vz[mask_p])
    #endfor

    # Define common histogram bins
    bins = np.linspace(-15, 15, 100)

    # Ensure output directory exists
    os.makedirs("output/rgc_studies", exist_ok=True)

    # Create a multi-page PDF
    with PdfPages("output/rgc_studies/particle_vz.pdf") as pdf:
        # --- Electron vertex distributions ---
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        for data, label, color in zip(electron_vz, labels, colors):
            axes[0].hist(data, bins=bins, density=True,
                         histtype="step", color=color, label=label)
            axes[1].hist(data, bins=bins, density=True,
                         histtype="step", color=color, label=label)
        #endfor

        # Linear‐scale panel
        axes[0].set_xlabel(r"$v_{z}$ (cm)")
        axes[0].set_ylabel("Normalized Counts")
        axes[0].set_title("Electron Vertex Distribution")
        axes[0].legend()

        # Log‐scale panel
        axes[1].set_yscale("log")
        axes[1].set_xlabel(r"$v_{z}$ (cm)")
        axes[1].set_ylabel("Normalized Counts")
        axes[1].set_title("Electron Vertex Distribution (log scale)")
        axes[1].legend()

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Proton vertex distributions ---
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        for data, label, color in zip(proton_vz, labels, colors):
            axes[0].hist(data, bins=bins, density=True,
                         histtype="step", color=color, label=label)
            axes[1].hist(data, bins=bins, density=True,
                         histtype="step", color=color, label=label)
        #endfor

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
        pdf.savefig(fig)
        plt.close(fig)
    #endwith

if __name__ == "__main__":
    main()
# endif