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

    plt.suptitle("MC electrons (pid==11), $\\theta\\leq 10^\\circ$, traj\\_edge\\_6 $\\leq$ 10: track\\_chi2\\_6", fontsize=14)
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
            h = ax.hist2d([], [], bins=[theta_bins, chi2_bins], norm=LogNorm(), cmap="jet")
            last_h = h
        # endif

        ax.set_xlim(*theta_xlim)   # 4° → 20°
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

    plt.suptitle("MC electrons (pid==11), $\\theta\\leq 20^\\circ$, traj\\_edge\\_6 $\\leq$ 10:  track\\_chi2\\_6 vs $\\theta$", fontsize=14)
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
    out_hist_2d = os.path.join(outdir, "mc_track_chi2_vs_theta_le20_te6_le10_electrons_by_sector.pdf")

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
    # Selections:
    #   - 1D: pid==11, theta <= 10, traj_edge_6 <= 10, valid sector 1..6
    #   - 2D: pid==11, theta <= 20, traj_edge_6 <= 10, valid sector 1..6
    # -------------------------------------------------------------------------
    valid_sector = (sector6 >= 1) & (sector6 <= 6)

    sel_1d = (pid == 11) & (theta <= 10.0) & (te6 <= 10.0) & valid_sector
    sel_2d = (pid == 11) & (theta <= 20.0) & (te6 <= 10.0) & valid_sector

    chi2_sel_1d   = chi2[sel_1d]
    theta_sel_1d  = theta[sel_1d]
    sector_sel_1d = sector6[sel_1d]

    chi2_sel_2d   = chi2[sel_2d]
    theta_sel_2d  = theta[sel_2d]
    sector_sel_2d = sector6[sel_2d]

    # Clean NaN/inf just in case
    good1 = np.isfinite(chi2_sel_1d) & np.isfinite(theta_sel_1d) & np.isfinite(sector_sel_1d)
    chi2_sel_1d   = chi2_sel_1d[good1]
    theta_sel_1d  = theta_sel_1d[good1]
    sector_sel_1d = sector_sel_1d[good1].astype(int)

    good2 = np.isfinite(chi2_sel_2d) & np.isfinite(theta_sel_2d) & np.isfinite(sector_sel_2d)
    chi2_sel_2d   = chi2_sel_2d[good2]
    theta_sel_2d  = theta_sel_2d[good2]
    sector_sel_2d = sector_sel_2d[good2].astype(int)

    print(f"[INFO] Total entries: {len(pid):,}")
    print(f"[INFO] 1D selection (θ≤10°, te6≤10): {chi2_sel_1d.size:,}")
    print(f"[INFO] 2D selection (θ≤20°, te6≤10): {chi2_sel_2d.size:,}")

    # -------------------------------------------------------------------------
    # Handle empty cases with placeholder plots
    # -------------------------------------------------------------------------
    if chi2_sel_1d.size == 0:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.text(0.5, 0.5, "No events (1D selection):\n"
                          "pid==11, θ≤10°, traj_edge_6≤10, valid sector",
                ha="center", va="center", fontsize=14)
        ax.set_axis_off()
        fig.savefig(out_hist_1d, bbox_inches="tight")
        plt.close(fig)
        print(f"[WARN] No 1D-selection events. Wrote placeholder to {out_hist_1d}")
    else:
        # Robust chi2 range for 1D (from 99th percentile)
        p99_1d = np.percentile(chi2_sel_1d, 99)
        chi2_max_1d = float(np.clip(p99_1d * 1.1, 5.0, 500.0))
        chi2_min_1d = 0.0
        if chi2_max_1d <= chi2_min_1d + 1e-9:
            chi2_max_1d = chi2_min_1d + 1.0
        # endif
        chi2_bins_1d = np.linspace(chi2_min_1d, chi2_max_1d, 101)

        plot_1d_by_sector(
            chi2_sel=chi2_sel_1d,
            sector_sel=sector_sel_1d,
            chi2_bins=chi2_bins_1d,
            chi2_xlim=(chi2_min_1d, chi2_max_1d),
            outpath=out_hist_1d,
        )
    # endif

    if chi2_sel_2d.size == 0:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.text(0.5, 0.5, "No events (2D selection):\n"
                          "pid==11, θ≤20°, traj_edge_6≤10, valid sector",
                ha="center", va="center", fontsize=14)
        ax.set_axis_off()
        fig.savefig(out_hist_2d, bbox_inches="tight")
        plt.close(fig)
        print(f"[WARN] No 2D-selection events. Wrote placeholder to {out_hist_2d}")
    else:
        # Robust chi2 range for 2D (from 99th percentile)
        p99_2d = np.percentile(chi2_sel_2d, 99)
        chi2_max_2d = float(np.clip(p99_2d * 1.1, 5.0, 500.0))
        chi2_min_2d = 0.0
        if chi2_max_2d <= chi2_min_2d + 1e-9:
            chi2_max_2d = chi2_min_2d + 1.0
        # endif

        theta_bins = np.linspace(4.0, 20.0, 101)     # 4° → 20° as requested
        chi2_bins  = np.linspace(chi2_min_2d, chi2_max_2d, 101)

        plot_2d_by_sector(
            theta_sel=theta_sel_2d,
            chi2_sel=chi2_sel_2d,
            sector_sel=sector_sel_2d,
            theta_bins=theta_bins,
            chi2_bins=chi2_bins,
            theta_xlim=(4.0, 20.0),                  # x-limits 4° → 20°
            chi2_ylim=(chi2_min_2d, chi2_max_2d),
            outpath=out_hist_2d,
        )
    # endif

# endif

if __name__ == "__main__":
    main()
# endif