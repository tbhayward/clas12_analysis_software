#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def main():
    # Paths to ROOT files and run labels
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

    # Vertex cuts for negative-particle sample
    vz_cuts = {
        "Su22": (-7.576, 0.303),
        "Fa22": (-5.758, 1.515),
        "Sp23": (-5.758, 1.515),
    }

    # Binning for diagonal cut: x 0–0.25 with 150 bins, y 0–0.20 with 150 bins
    x_bins = np.linspace(0, 0.25, 150)
    y_bins = np.linspace(0, 0.20, 150)

    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    for fname, label in zip(files, labels):
        # Read all needed branches at once
        arr = uproot.open(fname)[tree_name].arrays(
            [
                "particle_pid",
                "particle_vz",
                "track_sector_6",
                "p",
                "cc_nphe_15",
                "cal_energy_1",
                "cal_energy_4",
                "cal_lv_1",
                "cal_lw_1",
                "traj_edge_18",
                "traj_edge_36",
                "traj_edge_6",
                "theta"
            ],
            library="np"
        )

        pid     = arr["particle_pid"]
        vz      = arr["particle_vz"]
        sector6 = arr["track_sector_6"]
        p       = arr["p"]
        nphe    = arr["cc_nphe_15"]
        e1      = arr["cal_energy_1"]
        e4      = arr["cal_energy_4"]
        lv1     = arr["cal_lv_1"]
        lw1     = arr["cal_lw_1"]
        te18    = arr["traj_edge_18"]
        te36    = arr["traj_edge_36"]
        te6     = arr["traj_edge_6"]
        theta   = arr["theta"]

        # Fiducial cuts
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

        # Combined selection mask, with momentum > 4.9 GeV
        mask_all = (
            ((pid == 11)   | (pid == -211) | (pid == -321) | (pid == -2212)) &
            valid_sector &
            fid &
            (vz >= vz_cuts[label][0]) &
            (vz <= vz_cuts[label][1]) &
            (p > 4.9) &               # only above 4.9 GeV
            (nphe >= 2) &
            (e1 >= 0.15) &
            (e4 >= 0)
        )

        # Compute fraction variables
        frac_pcal = e1[mask_all] / p[mask_all]
        frac_ecin = e4[mask_all] / p[mask_all]
        sectors   = sector6[mask_all]

        # Create 2x3 plot: one pad per sector
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)
        for sec in range(1, 7):
            ax = axes.flat[sec-1]
            sel = (sectors == sec)
            x = frac_pcal[sel]
            y = frac_ecin[sel]
            h = ax.hist2d(
                x, y,
                bins=[x_bins, y_bins],
                cmap="jet",
                norm=LogNorm()
            )
            # Diagonal cut line from (0,0.15) to (0.23,0)
            ax.plot([0, 0.23], [0.15, 0], color="red", linestyle="-", linewidth=2, zorder=10)

            ax.set_title(f"{label} Sector {sec}")
            ax.set_xlabel(r"$E_{\mathrm{PCal}}/p$")
            ax.set_ylabel(r"$E_{\mathrm{ECin}}/p$")
            ax.set_xlim(0, 0.25)
            ax.set_ylim(0, 0.20)

        # Shared colorbar
        cb = fig.colorbar(h[3], ax=axes.ravel().tolist(), shrink=0.9)
        cb.set_label("Counts (log scale)")

        # Save
        plt.suptitle(f"Diagonal Cut: PCal vs ECin Fractions - {label}", fontsize=16)
        outpath = os.path.join(outdir, f"diagonal_cut_{label}.pdf")
        fig.savefig(outpath)
        plt.close(fig)

if __name__ == "__main__":
    main()