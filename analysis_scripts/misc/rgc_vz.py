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
        "Su22": (-7.576, 0.303),
        "Fa22": (-5.758, 1.515),
        "Sp23": (-5.758, 1.515),
    }

    # bins for nphe (50 bins from 0 to 40)
    bins = np.linspace(0, 40, 51)

    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    colors = ["C0", "C1", "C2"]

    for fname, label, color in zip(files, labels, colors):
        tree = uproot.open(fname)[tree_name]
        n_entries = tree.num_entries
        counts = np.zeros(len(bins) - 1, dtype=float)
        processed = 0

        # iterate in chunks to avoid loading entire arrays
        for arrays in tree.iterate(
            ["particle_pid", "particle_vz", "track_sector_6", "p", "cc_nphe_15"],
            step_size=1_000_000,
            library="np"
        ):
            pid     = arrays["particle_pid"]
            vz      = arrays["particle_vz"]
            sector6 = arrays["track_sector_6"]
            p       = arrays["p"]
            nphe    = arrays["cc_nphe_15"]

            # apply mask: negative tracks, vertex, momentum, valid nphe
            mask = (
                ((pid == 11)   | (pid == -211) | (pid == -321)) &
                (sector6 != -9999) &
                (vz >= vz_cuts[label][0]) &
                (vz <= vz_cuts[label][1]) &
                (p > 2.0) &
                (nphe != -9999)
            )
            sel = nphe[mask]
            counts += np.histogram(sel, bins=bins)[0]

            processed += len(pid)
            pct = 100 * processed / n_entries
            print(f"{label}: processed {processed}/{n_entries} ({pct:.1f}%)")

        # normalize
        total = counts.sum()
        densities = counts / total

        # plot as step histograms
        centers = 0.5 * (bins[:-1] + bins[1:])
        axes[0].step(centers, densities, where="mid", color=color, label=label)
        axes[1].step(centers, densities, where="mid", color=color, label=label)

    # add vertical line at nphe = 2
    for ax in axes:
        ax.axvline(2, color='red', linestyle='-', linewidth=2)
        ax.set_xlim(0, 40)

    # left panel: linear scale
    axes[0].set_xlabel("HTCC nphe")
    axes[0].set_ylabel("Normalized Counts")
    axes[0].set_title("Negative Tracks: HTCC nphe")
    axes[0].legend()

    # right panel: log scale
    axes[1].set_yscale("log")
    axes[1].set_xlabel("HTCC nphe")
    axes[1].set_ylabel("Normalized Counts")
    axes[1].set_title("Negative Tracks: HTCC nphe (log scale)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{outdir}/electron_htcc.pdf")
    plt.close(fig)

if __name__ == "__main__":
    main()