#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    # Paths to your ROOT files and run labels (RGA)
    files = [
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_fa18_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_fa18_out_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_sp19_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_sp18_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_sp18_out_calibration.root",
    ]
    labels = ["Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"]
    tree_name = "PhysicsEvents"

    # vertex‐cut thresholds for negative‐particle sample (new RGA cuts)
    vz_cuts = {
        "Fa18 Inb": (-6.364, 1.515),
        "Fa18 Out": (-7.879, 0.303),
        "Sp19 Inb": (-6.364, 1.515),
        "Sp18 Inb": (-6.06060606060606, 1.8181818181818183),
        "Sp18 Out": (-7.2727272727272725, 0.9090909090909101),
    }

    # histogram settings for PCal: 100 bins 0–1.5 GeV
    bins = np.linspace(0, 1.5, 100)

    # prepare output directory
    outdir = "output/rga_studies"
    os.makedirs(outdir, exist_ok=True)

    # collect PCal energies per run
    pcal_energies = []
    for fname, label in zip(files, labels):
        arr = uproot.open(fname)[tree_name].arrays([
            "particle_pid",
            "particle_vz",
            "track_sector_6",
            "p",
            "cc_nphe_16",
            "cal_energy_1",
            "cal_lv_1",
            "cal_lw_1",
            "traj_edge_18",
            "traj_edge_36",
            "traj_edge_6",
            "theta"
        ], library="np")

        pid     = arr["particle_pid"]
        vz      = arr["particle_vz"]
        sector6 = arr["track_sector_6"]
        p       = arr["p"]
        nphe    = arr["cc_nphe_16"]
        e1      = arr["cal_energy_1"]
        lv1     = arr["cal_lv_1"]
        lw1     = arr["cal_lw_1"]
        te18    = arr["traj_edge_18"]
        te36    = arr["traj_edge_36"]
        te6     = arr["traj_edge_6"]
        theta   = arr["theta"]

        # fiducial cuts
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

        # mask: all negative particles e⁻, π⁻, K⁻, p̄
        mask = (
            ((pid == 11)   | (pid == -211) | (pid == -321) | (pid == -2212)) &
            valid_sector &
            (vz >= vz_cuts[label][0]) &
            (vz <= vz_cuts[label][1]) &  # vertex window
            (p > 2.0) &                   # momentum > 2 GeV
            (nphe >= 2) &                 # HTCC veto nphe>=2
            (e1 >= 0) &                   # valid PCal energy
            fid                            # fiducial cuts
        )

        pcal_energies.append(e1[mask])

    # now plot all five on the same canvas
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    colors = ["C0", "C1", "C2", "C3", "C4"]

    for data, label, color in zip(pcal_energies, labels, colors):
        axes[0].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)
        axes[1].hist(data, bins=bins, density=True,
                     histtype="step", color=color, label=label)

    # add vertical red line at 0.06 GeV
    for ax in axes:
        ax.axvline(0.06, color="red", linestyle="-", linewidth=2)

    # formatting panels
    axes[0].set_xlabel(r"$E_{\mathrm{PCal}}$ (GeV)")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Negative Particles: PCal Energy Deposition")
    axes[0].legend()

    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$E_{\mathrm{PCal}}$ (GeV)")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Negative Particles: PCal Energy Deposition (log scale)")
    axes[1].legend()

    fig.tight_layout()
    # keep the same save filename
    fig.savefig(f"{outdir}/negative_particles_pcal_deposition.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()