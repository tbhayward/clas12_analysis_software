#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    # file paths and labels for RGA
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

    # RGA vertex cuts for negative particles
    vz_cuts = {
        "Fa18 Inb":    (-6.364, 1.515),
        "Fa18 Out":    (-7.879, 0.303),
        "Sp19 Inb":    (-6.364, 1.515),
    }

    # HTCC nphe histogram bins
    bins = np.linspace(0, 40, 51)

    outdir = "output/rga_studies"
    os.makedirs(outdir, exist_ok=True)

    # prepare figure with linear and log panels
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    colors = ["C0", "C1", "C2"]

    for fname, label, color in zip(files, labels, colors):
        tree = uproot.open(fname)[tree_name]
        n_entries = tree.num_entries
        counts = np.zeros(len(bins) - 1, dtype=float)
        processed = 0

        # iterate in chunks for efficiency
        for arrays in tree.iterate(
            [
                "particle_pid", "particle_vz", "track_sector_6",
                "p", "cc_nphe_15",
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
            nphe    = arrays["cc_nphe_15"]
            lv1     = arrays["cal_lv_1"]
            lw1     = arrays["cal_lw_1"]
            te18    = arrays["traj_edge_18"]
            te36    = arrays["traj_edge_36"]
            te6     = arrays["traj_edge_6"]
            theta   = arrays["theta"]

            # fiducial cuts differ for out‑bending vs in‑bending
            if label == "Fa18 Out":
                # always require te6 > 3
                fid = (
                    (lv1 > 9) &
                    (lw1 > 9) &
                    (te18 > 3) &
                    (te36 > 10) &
                    (te6 > 3)
                )
            else:
                # in‑bending: theta‑dependent
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

            # mask negative tracks + selections
            mask = (
                ((pid == 11)   |
                 (pid == -211) |
                 (pid == -321) |
                 (pid == -2212)) &
                (sector6 != -9999) &
                (vz >= vz_cuts[label][0]) &
                (vz <= vz_cuts[label][1]) &
                (p > 2.0) &
                (nphe != -9999) &
                fid
            )

            # accumulate histogram
            counts += np.histogram(nphe[mask], bins=bins)[0]

            # progress update
            processed += len(pid)
            pct = 100 * processed / n_entries
            print(f"{label}: processed {processed}/{n_entries} ({pct:.1f}%)", flush=True)

        # normalize to density
        densities = counts / counts.sum()
        centers = 0.5 * (bins[:-1] + bins[1:])
        axes[0].step(centers, densities, where="mid", color=color, label=label)
        axes[1].step(centers, densities, where="mid", color=color, label=label)

    # draw vertical nphe cut at 2
    for ax in axes:
        ax.axvline(2, color="red", linestyle="-", linewidth=2)
        ax.set_xlim(0, 40)

    # left: linear scale
    axes[0].set_xlabel("HTCC nphe")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Negative Particles: HTCC nphe (RGA)")
    axes[0].legend()

    # right: log scale
    axes[1].set_yscale("log")
    axes[1].set_xlabel("HTCC nphe")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Negative Particles: HTCC nphe (log, RGA)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{outdir}/negative_particles_htcc_rga.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()
