#!/usr/bin/env python3

import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def main():
    # -------------------------------------------------------------------------
    # Input file / tree
    # -------------------------------------------------------------------------
    mc_file = (
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "clasdis_rgc_fa22_inb_calibration.root"
    )
    tree_name = "PhysicsEvents"

    # Output paths (same directory as previous scripts)
    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)
    out_hist_1d = os.path.join(outdir, "mc_track_chi2_theta_le10_te6_le10_electrons.pdf")
    out_hist_2d = os.path.join(outdir, "mc_track_chi2_vs_theta_le10_te6_le10_electrons.pdf")

    # -------------------------------------------------------------------------
    # Read only the branches we need
    # -------------------------------------------------------------------------
    branches = ["particle_pid", "theta", "traj_edge_6", "track_chi2_6"]

    with uproot.open(mc_file) as f:
        tree = f[tree_name]
        arr = tree.arrays(branches, library="np")
    # endif

    pid   = arr["particle_pid"]
    theta = arr["theta"]
    te6   = arr["traj_edge_6"]
    chi2  = arr["track_chi2_6"]

    # -------------------------------------------------------------------------
    # Selection: electrons (pid == 11) with theta <= 10 degrees AND traj_edge_6 <= 10
    # -------------------------------------------------------------------------
    sel = (pid == 11) & (theta <= 10.0) & (te6 <= 10.0)
    chi2_sel  = chi2[sel]
    theta_sel = theta[sel]

    # Guard against NaN/inf just in case
    good = np.isfinite(chi2_sel) & np.isfinite(theta_sel)
    chi2_sel  = chi2_sel[good]
    theta_sel = theta_sel[good]

    n_total = len(pid)
    n_pass  = chi2_sel.size
    print(f"[INFO] Total entries: {n_total}")
    print(f"[INFO] Selected electrons (pid==11) with theta <= 10° and traj_edge_6 <= 10: {n_pass}")

    # -------------------------------------------------------------------------
    # If no entries, write an empty placeholder plot for both outputs and exit
    # -------------------------------------------------------------------------
    if n_pass == 0:
        for out in (out_hist_1d, out_hist_2d):
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.text(0.5, 0.5, "No events with criteria:\n"
                              "pid==11, θ≤10°, traj_edge_6≤10",
                    ha="center", va="center", fontsize=14)
            ax.set_axis_off()
            fig.savefig(out, bbox_inches="tight")
            plt.close(fig)
            print(f"[WARN] No matching events. Wrote empty plot to {out}")
        # endfor
        return
    # endif

    # -------------------------------------------------------------------------
    # Determine robust chi2 plotting range using 99th percentile
    # -------------------------------------------------------------------------
    p99 = np.percentile(chi2_sel, 99)
    chi2_upper = float(np.clip(p99 * 1.1, 5.0, 500.0))
    chi2_lower = 0.0
    if chi2_upper <= chi2_lower + 1e-9:
        chi2_upper = chi2_lower + 1.0
    # endif

    # -------------------------------------------------------------------------
    # 1D histogram of track_chi2_6
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    bins_1d = 100
    ax.hist(chi2_sel, bins=bins_1d, range=(chi2_lower, chi2_upper),
            histtype="stepfilled", alpha=0.85)
    ax.set_xlabel(r"track$\_\chi^2$ (sector 6 fit)", fontsize=12)
    ax.set_ylabel("Counts", fontsize=12)
    ax.set_title(r"MC electrons (pid==11), $\theta\le 10^\circ$, traj\_edge\_6 $\le 10$: track$\_\chi^2\_6$", fontsize=13)
    ax.grid(True, linestyle="--", alpha=0.4)

    ax.annotate(f"Entries: {n_pass:,}\n"
                f"χ² range: [{chi2_lower:.2f}, {chi2_upper:.2f}]",
                xy=(0.98, 0.97), xycoords="axes fraction",
                ha="right", va="top", fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.85))
    plt.tight_layout()
    fig.savefig(out_hist_1d)
    plt.close(fig)
    print(f"[OK] Saved 1D χ² histogram to: {out_hist_1d}")

    # -------------------------------------------------------------------------
    # 2D histogram: chi2 vs theta (theta on x, chi2 on y)
    # -------------------------------------------------------------------------
    # Binning: theta from 0→10 (since we cut), chi2 using the same robust upper
    th_lo, th_hi = 0.0, 10.0
    th_bins = 100
    chi2_bins = 100

    fig, ax = plt.subplots(figsize=(8, 6))
    h = ax.hist2d(theta_sel, chi2_sel,
                  bins=[np.linspace(th_lo, th_hi, th_bins+1),
                        np.linspace(chi2_lower, chi2_upper, chi2_bins+1)],
                  norm=LogNorm(), cmap="jet")
    cb = fig.colorbar(h[3], ax=ax)
    cb.set_label("Counts (log scale)")

    ax.set_xlabel(r"$\theta$ [deg]", fontsize=12)
    ax.set_ylabel(r"track$\_\chi^2\_6$", fontsize=12)
    ax.set_title(r"MC electrons (pid==11), $\theta\le 10^\circ$, traj\_edge\_6 $\le 10$:  $\chi^2$ vs $\theta$", fontsize=13)
    ax.grid(True, linestyle="--", alpha=0.3)

    plt.tight_layout()
    fig.savefig(out_hist_2d)
    plt.close(fig)
    print(f"[OK] Saved 2D χ² vs θ histogram to: {out_hist_2d}")

# endif

if __name__ == "__main__":
    main()
# endif