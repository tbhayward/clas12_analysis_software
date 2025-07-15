#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def main():
    # Su22 file and Su22-specific vertex cuts
    fname = (
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "sidisdvcs_rgc_su22_inb_calibration.root"
    )
    vz_cut = (-7.576, 0.303)
    tree_name = "PhysicsEvents"

    # Read the needed branches in one go
    arr = uproot.open(fname)[tree_name].arrays([
        "particle_pid", "particle_vz", "cal_sector",
        "p", "cc_nphe_15",
        "cal_energy_1", "cal_energy_4", "cal_energy_7",
        "cal_lv_1", "cal_lw_1",
        "traj_edge_18", "traj_edge_36", "traj_edge_6",
        "theta"
    ], library="np")

    pid    = arr["particle_pid"]
    vz     = arr["particle_vz"]
    sector = arr["cal_sector"]
    p      = arr["p"]
    nphe   = arr["cc_nphe_15"]
    e1     = arr["cal_energy_1"]
    e4     = arr["cal_energy_4"]
    e7     = arr["cal_energy_7"]
    lv1    = arr["cal_lv_1"]
    lw1    = arr["cal_lw_1"]
    te18   = arr["traj_edge_18"]
    te36   = arr["traj_edge_36"]
    te6    = arr["traj_edge_6"]
    theta  = arr["theta"]

    # Fiducial cuts
    fid = (
        (lv1 > 9) & (lw1 > 9) &
        (te18 > 3) & (te36 > 10) &
        (
            ((theta > 10) & (te6 > 3)) |
            ((theta <= 10) & (te6 > 10))
        )
    )
    valid_sector = (sector != -9999)

    # Base electron+cuts mask
    mask_base = (
        (pid == 11) &
        valid_sector &
        fid &
        (vz >= vz_cut[0]) & (vz <= vz_cut[1]) &
        (p > 2.0) &
        (nphe >= 2) &
        (e1 >= 0.15) &
        (e4 >= 0) &
        (e7 >= 0)
    )

    # Diagonal cut
    p_base    = p[mask_base]
    frac_pcal = e1[mask_base] / p_base
    frac_ecin = e4[mask_base] / p_base
    mask_diag = frac_ecin >= (-0.625 * frac_pcal + 0.15)

    # Combine base + diagonal
    idx0 = np.nonzero(mask_base)[0]
    keep = idx0[mask_diag]

    # Sector-1 subset
    sec1 = keep[ sector[keep] == 1 ]

    # Additional SF < 0.28 cut for sector 1
    sf_sec1 = (e1[sec1] + e4[sec1] + e7[sec1]) / p[sec1]
    sec1 = sec1[sf_sec1 < 0.50]

    # Final arrays
    p1  = p[sec1]
    sf1 = (e1[sec1] + e4[sec1] + e7[sec1]) / p1

    # Plot 2D histogram of SF vs p for sector 1
    fig, ax = plt.subplots(figsize=(8, 6))
    p_bins = np.linspace(2.0, 8.0, 60)
    sf_bins = np.linspace(0.10, 0.40, 60)

    h = ax.hist2d(
        p1, sf1,
        bins=[p_bins, sf_bins],
        cmap="jet", norm=LogNorm()
    )

    # draw the SF<0.28 cut line
    ax.axhline(0.28, color='red', linestyle='-', linewidth=2)

    ax.set_xlim(2.0, 8.0)
    ax.set_ylim(0.10, 0.40)
    ax.set_xlabel("p (GeV)")
    ax.set_ylabel("Sampling Fraction")
    ax.set_title("Su22 Sector 1: SF vs p")

    cbar = fig.colorbar(h[3], ax=ax)
    cbar.set_label("Counts (log scale)")

    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(f"{outdir}/weird_su22.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()