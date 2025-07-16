#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def main():
    # Paths to RGA ROOT files and labels
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

    # Vertex cuts for negative tracks (RGA)
    vz_cuts = {
        "Fa18 Inb": (-6.364, 1.515),
        "Fa18 Out": (-7.879, 0.303),
        "Sp19 Inb": (-6.364, 1.515),
    }

    # Binning for diagonal cut: 150×150 over [0,0.25] × [0,0.20]
    x_bins = np.linspace(0, 0.25, 150)
    y_bins = np.linspace(0, 0.20, 150)

    outdir = "output/rga_studies"
    os.makedirs(outdir, exist_ok=True)

    for fname, label in zip(files, labels):
        # load all needed branches in one go
        arr = uproot.open(fname)[tree_name].arrays([
            "particle_pid","particle_vz","track_sector_6",
            "p","cc_nphe_15",
            "cal_energy_1","cal_energy_4",
            "cal_lv_1","cal_lw_1",
            "traj_edge_18","traj_edge_36","traj_edge_6",
            "theta",
        ], library="np")

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

        # fiducial cuts differ for Out vs Inbending
        if label == "Fa18 Out":
            fid = (
                (lv1 > 9) & (lw1 > 9) &
                (te18 > 3) & (te36 > 10) &
                (te6 > 3)
            )
        else:
            fid = (
                (lv1 > 9) & (lw1 > 9) &
                (te18 > 3) & (te36 > 10) &
                (((theta > 10) & (te6 > 3)) |
                 ((theta <= 10) & (te6 > 10)))
            )
        valid_sector = (sector6 != -9999)

        # full selection mask
        mask_all = (
            ((pid==11)|(pid==-211)|(pid==-321)|(pid==-2212)) &
            valid_sector & fid &
            (vz >= vz_cuts[label][0]) & (vz <= vz_cuts[label][1]) &
            (p > 4.9) &
            (nphe >= 2) &
            (e1 >= 0.06) &
            (e4 >= 0)
        )

        frac_pcal = e1[mask_all] / p[mask_all]
        frac_ecin = e4[mask_all] / p[mask_all]
        sectors   = sector6[mask_all]

        # 2x3 sector plot
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)
        last_quadmesh = None

        for sec in range(1, 7):
            ax = axes.flat[sec-1]
            sel = (sectors == sec)
            x   = frac_pcal[sel]
            y   = frac_ecin[sel]

            # do the 2D histogram
            counts2d, xedges, yedges, quadmesh = ax.hist2d(
                x, y,
                bins=[x_bins, y_bins],
                cmap="jet",
                norm=LogNorm()
            )

            # compute a safe vmin/vmax from nonzero bins
            nz = counts2d[counts2d > 0]
            if nz.size:
                quadmesh.set_norm(LogNorm(vmin=nz.min(), vmax=nz.max()))
            else:
                # if entirely empty, set a dummy range to avoid errors
                quadmesh.set_norm(LogNorm(vmin=1, vmax=1))

            last_quadmesh = quadmesh

            # diagonal cut line
            ax.plot([0,0.24], [0.15,0], 'r-', lw=2, zorder=10)

            ax.set_title(f"{label}, Sector {sec}")
            ax.set_xlabel(r"$E_{\mathrm{PCal}}/p$")
            ax.set_ylabel(r"$E_{\mathrm{ECin}}/p$")
            ax.set_xlim(0, 0.25)
            ax.set_ylim(0, 0.20)

        # colorbar from the last QuadMesh
        cb = fig.colorbar(last_quadmesh, ax=axes.ravel().tolist(), shrink=0.9)
        cb.set_label("Counts (log scale)")

        fig.suptitle(f"RGA Diagonal Cut: {label}", fontsize=16)
        outpath = os.path.join(outdir, f"diagonal_cut_{label.replace(' ', '_')}.pdf")
        fig.savefig(outpath)
        plt.close(fig)

if __name__ == "__main__":
    main()