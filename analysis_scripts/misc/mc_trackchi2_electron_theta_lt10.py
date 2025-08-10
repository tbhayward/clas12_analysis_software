#!/usr/bin/env python3

import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    # -------------------------------------------------------------------------
    # Input file / tree
    # -------------------------------------------------------------------------
    mc_file = (
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "clasdis_rgc_fa22_inb_calibration.root"
    )
    tree_name = "PhysicsEvents"

    # Output path (same place as the previous script)
    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "mc_track_chi2_theta_lt10_electrons.pdf")

    # -------------------------------------------------------------------------
    # Read only the branches we need
    # -------------------------------------------------------------------------
    branches = ["particle_pid", "theta", "track_chi2_6"]

    with uproot.open(mc_file) as f:
        tree = f[tree_name]
        arr = tree.arrays(branches, library="np")
    # endif

    pid   = arr["particle_pid"]
    theta = arr["theta"]
    chi2  = arr["track_chi2_6"]

    # -------------------------------------------------------------------------
    # Selection: electrons (pid == 11) with theta < 10 degrees
    # -------------------------------------------------------------------------
    sel = (pid == 11) & (theta < 10.0)
    chi2_sel = chi2[sel]

    # Guard against NaN/inf just in case
    chi2_sel = chi2_sel[np.isfinite(chi2_sel)]

    n_total = len(pid)
    n_pass  = chi2_sel.size
    print(f"[INFO] Total entries: {n_total}")
    print(f"[INFO] Selected electrons with theta < 10°: {n_pass}")

    # -------------------------------------------------------------------------
    # Choose binning: robust range based on 99th percentile (with a floor/ceiling)
    # -------------------------------------------------------------------------
    if n_pass == 0:
        # No entries — make an empty plot with a clear message
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.text(0.5, 0.5, "No events with pid==11 and θ<10°", ha="center", va="center", fontsize=14)
        ax.set_axis_off()
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        print(f"[WARN] No matching events. Wrote empty plot to {outfile}")
        return
    # endif

    p99 = np.percentile(chi2_sel, 99)
    # Keep a sane plotting upper bound
    upper = float(np.clip(p99 * 1.1, 5.0, 500.0))
    lower = 0.0

    bins = 100
    # If the spread is tiny, fall back to an automatic bin range
    if upper <= lower + 1e-9:
        upper = lower + 1.0
    # endif

    # -------------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(chi2_sel, bins=bins, range=(lower, upper), histtype="stepfilled", alpha=0.8)
    ax.set_xlabel(r"track$\_\chi^2$ (sector 6 fit)", fontsize=12)
    ax.set_ylabel("Counts", fontsize=12)
    ax.set_title(r"MC electrons (pid==11) with $\theta<10^\circ$: track$\_\chi^2\_6$", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.4)

    # Nice annotation with counts
    ax.annotate(f"Entries: {n_pass:,}\nRange: [{lower:.2f}, {upper:.2f}]",
                xy=(0.98, 0.97), xycoords="axes fraction",
                ha="right", va="top", fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

    plt.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    print(f"[OK] Saved histogram to: {outfile}")

# endif

if __name__ == "__main__":
    main()
# endif