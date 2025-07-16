#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    # Paths to the RGA ROOT files and their labels
    files = [
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_fa18_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_fa18_out_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_sp19_inb_calibration.root",
    ]
    labels = ["Fa18 Inb", "Fa18 Out", "Sp19 Inb"]
    tree_name = "PhysicsEvents"

    # Vertex‐cut thresholds for negative‐particle sample (RGA)
    vz_cuts = {
        "Fa18 Inb": (-6.364, 1.515),
        "Fa18 Out": (-7.879, 0.303),
        "Sp19 Inb": (-6.364, 1.515),
    }

    # PCal energy histogram settings: 100 bins from 0 to 1.5 GeV
    bins = np.linspace(0, 1.5, 100)

    # Prepare output directory
    outdir = "output/rga_studies"
    os.makedirs(outdir, exist_ok=True)

    # Collect PCal energies for each run period
    pcal_energies = []
    for fname, label in zip(files, labels):
        # Load needed branches into numpy arrays
        arr = uproot.open(fname)[tree_name].arrays([
            "particle_pid",
            "particle_vz",
            "track_sector_6",
            "p",
            "cc_nphe_15",
            "cal_energy_1",
            "cal_lv_1",
            "cal_lw_1",
            "traj_edge_18",
            "traj_edge_36",
            "traj_edge_6",
            "theta"
        ], library="np")

        pid       = arr["particle_pid"]
        vz        = arr["particle_vz"]
        sector6   = arr["track_sector_6"]
        p         = arr["p"]
        nphe      = arr["cc_nphe_15"]
        e1        = arr["cal_energy_1"]
        lv1       = arr["cal_lv_1"]
        lw1       = arr["cal_lw_1"]
        te18      = arr["traj_edge_18"]
        te36      = arr["traj_edge_36"]
        te6       = arr["traj_edge_6"]
        theta     = arr["theta"]

        # Fiducial cuts: out‐bending simpler, in‐bending depends on theta
        if label == "Fa18 Out":
            fid = (
                (lv1 > 9) &
                (lw1 > 9) &
                (te18 > 3) &
                (te36 > 10) &
                (te6 > 3)
            )
        else:
            fid = (
                (lv1 > 9) &
                (lw1 > 9) &
                (te18 > 3) &
                (te36 > 10) &
                (
                    ((theta > 10) & (te6 > 3)) |
                    ((theta <= 10) & (te6 > 10))
                )
            )

        valid_sector = (sector6 != -9999)

        # Build the mask for all negative species and all prior cuts
        mask = (
            ((pid == 11)   | (pid == -211) | (pid == -321) | (pid == -2212)) &
            valid_sector &
            fid &
            (vz >= vz_cuts[label][0]) &
            (vz <= vz_cuts[label][1]) &
            (p > 2.0) &
            (nphe >= 2) &
            (e1 >= 0)    # ensure valid PCal energy
        )

        # Store the PCal energies passing the mask
        pcal_energies.append(e1[mask])

    # Now plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    colors = ["C0", "C1", "C2"]

    for data, label, color in zip(pcal_energies, labels, colors):
        axes[0].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)
        axes[1].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)

    # Add vertical red line at 0.07 GeV
    for ax in axes:
        ax.axvline(0.07, color="red", linestyle="-", linewidth=2)

    # Left panel: linear scale
    axes[0].set_xlabel(r"$E_{\mathrm{PCal}}$ (GeV)")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Negative Particles: PCal Energy Deposition (RGA)")
    axes[0].legend()

    # Right panel: log scale
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$E_{\mathrm{PCal}}$ (GeV)")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Negative Particles: PCal Energy Deposition (log, RGA)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{outdir}/negative_particles_pcal_deposition_rga.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()
