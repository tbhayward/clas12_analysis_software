#!/usr/bin/env python3

import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_1d_by_sector(chi2_sel, sector_sel, chi2_bins, chi2_xlim, outpath):
    """
    Per-sector 1D histograms of track_chi2_6 with shared x-range and bins.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for sec in range(1, 7):
        ax = axes.flat[sec-1]
        m = (sector_sel == sec)
        data = chi2_sel[m]
        ax.hist(data, bins=chi2_bins, histtype="stepfilled", alpha=0.85)
        ax.set_xlim(*chi2_xlim)
        ax.set_xlabel("track_chi2_6")
        ax.set_ylabel("Counts")
        ax.set_title(f"Sector {sec} (N={data.size:,})")
        ax.grid(True, linestyle="--", alpha=0.4)
    # endfor

    plt.suptitle("MC electrons (pid==11), $\\theta\\leq 10^\\circ$, traj_edge_6 $\\leq$ 10: track_chi2_6", fontsize=14)
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[OK] Saved per-sector 1D chi2 hist to: {outpath}")
# endif

def plot_2d_by_sector(theta_sel, chi2_sel, sector_sel, theta_bins, chi2_bins, theta_xlim, chi2_ylim, outpath):
    """
    Per-sector 2D histograms of track_chi2_6 vs theta with shared colorbar.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)
    last_h = None

    for sec in range(1, 7):
        ax = axes.flat[sec-1]
        m = (sector_sel == sec)
        x = theta_sel[m]
        y = chi2_sel[m]
        if x.size and y.size:
            h = ax.hist2d(x, y, bins=[theta_bins, chi2_bins], norm=LogNorm(), cmap="jet")
            last_h = h
        else:
            # empty container to keep axes consistent
            h = ax.hist2d([], [], bins=[theta_bins, chi2_bins], norm=LogNorm(), cmap="jet")
            last_h = h
        # endif

        ax.set_xlim(*theta_xlim)
        ax.set_ylim(*chi2_ylim)
        ax.set_xlabel(r"$\theta$ [deg]")
        ax.set_ylabel("track_chi2_6")
        ax.set_title(f"Sector {sec} (N={x.size:,})")
        ax.grid(True, linestyle="--", alpha=0.3)
    # endfor

    if last_h is not None:
        cb = fig.colorbar(last_h[3], ax=axes.ravel().tolist(), shrink=0.9)
        cb.set_label("Counts (log scale)")
    # endif

    plt.suptitle("MC electrons (pid==11), $\\theta\\leq 10^\\circ$, traj_edge_6 $\\leq$ 10:  track_chi2_6 vs $\\theta$", fontsize=14)
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[OK] Saved per-sector 2D chi2 vs theta to: {outpath}")
# endif

def main():
    # -------------------------------------------------------------------------
    # Input file / tree
    # -------------------------------------------------------------------------
    mc_file = (
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "clasdis_rgc_fa22_inb_calibration.root"
    )
    tree_name = "PhysicsEvents"

    # Output (same place as before)
    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)
    out_hist_1d = os.path.join(outdir, "mc_track_chi2_theta_le10_te6_le10_electrons_by_sector.pdf")
    out_hist_2d = os.path.join(outdir, "mc_track_chi2_vs_theta_le10_te6_le10_electrons_by_sector.pdf")

    # -------------------------------------------------------------------------
    # Read only needed branches
    # -------------------------------------------------------------------------
    branches = ["particle_pid", "theta", "traj_edge_6", "track_chi2_6", "track_sector_6"]
    with uproot.open(mc_file) as f:
        tree = f[tree_name]
        arr = tree.arrays(branches, library="np")
    # endif

    pid     = arr["particle_pid"]
    theta   = arr["theta"]
    te6     = arr["traj_edge_6"]
    chi2    = arr["track_chi2_6"]
    sector6 = arr["track_sector_6"]

    # -------------------------------------------------------------------------
    # Selection: pid==11, theta <= 10 deg, traj_edge_6 <= 10, valid sector 1..6
    # -------------------------------------------------------------------------
    valid_sector = (sector6 >= 1) & (sector6 <= 6)
    sel = (pid == 11) & (theta <= 10.0) & (te6 <= 10.0) & valid_sector

    chi2_sel   = chi2[sel]
    theta_sel  = theta[sel]
    sector_sel = sector6[sel]

    # Clean NaN/inf just in case
    good = np.isfinite(chi2_sel) & np.isfinite(theta_sel) & np.isfinite(sector_sel)
    chi2_sel   = chi2_sel[good]
    theta_sel  = theta_sel[good]
    sector_sel = sector_sel[good].astype(int)

    print(f"[INFO] Total entries: {len(pid):,}")
    print(f"[INFO] Selected electrons with theta <= 10° and traj_edge_6 <= 10 and valid sector: {chi2_sel.size:,}")

    if chi2_sel.size == 0:
        # Write empty placeholders and exit
        for out in (out_hist_1d, out_hist_2d):
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.text(0.5, 0.5, "No events with criteria:\n"
                              "pid==11, θ≤10°, traj_edge_6≤10, valid sector",
                    ha="center", va="center", fontsize=14)
            ax.set_axis_off()
            fig.savefig(out, bbox_inches="tight")
            plt.close(fig)
            print(f"[WARN] No matching events. Wrote empty plot to {out}")
        # endfor
        return
    # endif

    # -------------------------------------------------------------------------
    # Common binning/ranges
    # -------------------------------------------------------------------------
    # Robust chi2 max from 99th percentile (shared across sectors for comparability)
    p99 = np.percentile(chi2_sel, 99)
    chi2_max = float(np.clip(p99 * 1.1, 5.0, 500.0))
    chi2_min = 0.0
    if chi2_max <= chi2_min + 1e-9:
        chi2_max = chi2_min + 1.0
    # endif

    chi2_bins_1d = np.linspace(chi2_min, chi2_max, 101)
    theta_bins   = np.linspace(0.0, 10.0, 101)
    chi2_bins_2d = np.linspace(chi2_min, chi2_max, 101)

    # -------------------------------------------------------------------------
    # Plot per-sector 1D and 2D
    # -------------------------------------------------------------------------
    plot_1d_by_sector(
        chi2_sel=chi2_sel,
        sector_sel=sector_sel,
        chi2_bins=chi2_bins_1d,
        chi2_xlim=(chi2_min, chi2_max),
        outpath=out_hist_1d,
    )

    plot_2d_by_sector(
        theta_sel=theta_sel,
        chi2_sel=chi2_sel,
        sector_sel=sector_sel,
        theta_bins=theta_bins,
        chi2_bins=chi2_bins_2d,
        theta_xlim=(0.0, 10.0),
        chi2_ylim=(chi2_min, chi2_max),
        outpath=out_hist_2d,
    )

# endif

if __name__ == "__main__":
    main()
# endif