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
        "Su22": (-7.576,  0.303),
        "Fa22": (-5.758,  1.515),
        "Sp23": (-5.758,  1.515),
    }

    # bins for HTCC nphe: 50 bins 0–40
    bins = np.linspace(0, 40, 51)

    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    # prepare 1×2 figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    colors = ["C0", "C1", "C2"]

    for fname, label, color in zip(files, labels, colors):
        tree = uproot.open(fname)[tree_name]
        n_entries = tree.num_entries
        counts = np.zeros(len(bins) - 1, dtype=float)
        processed = 0

        # iterate in chunks, reading only needed branches
        for arrays in tree.iterate(
            [
                "particle_pid", "particle_vz", "track_sector_6",
                "p", "cc_nphe_16",
                "cal_lv_1", "cal_lw_1",
                "traj_edge_18", "traj_edge_36", "traj_edge_6",
                "theta"
            ],
            step_size="100 MB",
            library="np"
        ):
            pid     = arrays["particle_pid"]
            vz      = arrays["particle_vz"]
            sector6 = arrays["track_sector_6"]
            p       = arrays["p"]
            nphe    = arrays["cc_nphe_16"]
            lv1     = arrays["cal_lv_1"]
            lw1     = arrays["cal_lw_1"]
            te18    = arrays["traj_edge_18"]
            te36    = arrays["traj_edge_36"]
            te6     = arrays["traj_edge_6"]
            theta   = arrays["theta"]

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

            # mask negative tracks + sector + vertex + momentum + nphe + fiducial
            mask = (
                ((pid == 11) | (pid == -211) | (pid == -321) | (pid == -2212)) &
                (sector6 != -9999) &
                (vz >= vz_cuts[label][0]) &
                (vz <= vz_cuts[label][1]) &
                (p > 2.0) &
                (nphe != -9999) &
                fid
            )

            sel = nphe[mask]
            counts += np.histogram(sel, bins=bins)[0]

            # progress
            processed += len(pid)
            pct = 100 * processed / n_entries
            print(f"{label}: processed {processed}/{n_entries} ({pct:.1f}%)", flush=True)

        # normalize and plot
        densities = counts / counts.sum()
        centers = 0.5 * (bins[:-1] + bins[1:])
        axes[0].step(centers, densities, where="mid", color=color, label=label)
        axes[1].step(centers, densities, where="mid", color=color, label=label)

    # vertical cut at nphe = 2
    for ax in axes:
        ax.axvline(2, color="red", linestyle="-", linewidth=2)
        ax.set_xlim(0, 40)

    # left: linear
    axes[0].set_xlabel("HTCC nphe")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Negative Particles: HTCC nphe")
    axes[0].legend()

    # right: log
    axes[1].set_yscale("log")
    axes[1].set_xlabel("HTCC nphe")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Negative Particles: HTCC nphe (log scale)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{outdir}/negative_particles_htcc.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()