#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_by_sector(x_vals, y_vals, sectors, x_bins, y_bins, title, outpath):
    """
    Helper: make a 2x3 sectorized 2D histogram panel with a shared colorbar.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)

    # Track the last hist2d handle for the colorbar
    last_h = None

    for sec in range(1, 7):
        ax = axes.flat[sec-1]
        sel = (sectors == sec)
        x = x_vals[sel]
        y = y_vals[sel]

        if x.size > 0 and y.size > 0:
            h = ax.hist2d(
                x, y,
                bins=[x_bins, y_bins],
                cmap="jet",
                norm=LogNorm()
            )
            last_h = h
        else:
            # Draw an empty image so axes still render consistently
            h = ax.hist2d([], [], bins=[x_bins, y_bins], cmap="jet", norm=LogNorm())
            last_h = h
        # endif

        # Diagonal cut line from (0, 0.15) to (0.24, 0)
        ax.plot([0, 0.24], [0.15, 0], color="red", linestyle="-", linewidth=2, zorder=10)

        ax.set_title(f"Sector {sec}")
        ax.set_xlabel(r"$E_{\mathrm{PCal}}/p$")
        ax.set_ylabel(r"$E_{\mathrm{ECin}}/p$")
        ax.set_xlim(0, 0.25)
        ax.set_ylim(0, 0.20)
    # endfor

    # Shared colorbar (use the image object returned by hist2d: last_h[3])
    if last_h is not None:
        cb = fig.colorbar(last_h[3], ax=axes.ravel().tolist(), shrink=0.9)
        cb.set_label("Counts (log scale)")
    # endif

    plt.suptitle(title, fontsize=16)
    fig.savefig(outpath)
    plt.close(fig)
    print(f"Saved: {outpath}")

def main():
    # -------------------------------------------------------------------------
    # Single MC file and tree
    # -------------------------------------------------------------------------
    mc_file = (
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "clasdis_rgc_fa22_inb_calibration.root"
    )
    tree_name = "PhysicsEvents"
    label = "MC Fa22"

    # -------------------------------------------------------------------------
    # Vertex cuts & binning (re-using your Fa22/Sp23 window and original bins)
    # -------------------------------------------------------------------------
    vz_min, vz_max = (-5.758, 1.515)

    # Binning for diagonal cut: x 0–0.25 with 150 bins, y 0–0.20 with 150 bins
    x_bins = np.linspace(0, 0.25, 150)
    y_bins = np.linspace(0, 0.20, 150)

    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    # -------------------------------------------------------------------------
    # Read branches (same variables as before + mc_particle_id)
    # -------------------------------------------------------------------------
    branches = [
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
        "theta",
        "mc_particle_id",   # <- MC truth categorization
    ]

    with uproot.open(mc_file) as f:
        arr = f[tree_name].arrays(branches, library="np")
    # endif

    # Unpack
    pid     = arr["particle_pid"]
    sector6 = arr["track_sector_6"]
    vz      = arr["particle_vz"]
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
    mc_id   = arr["mc_matching_pid"]

    # -------------------------------------------------------------------------
    # Fiducial & quality cuts (same logic as your original)
    # -------------------------------------------------------------------------
    fid = (
        (lv1 > 9) &
        (lw1 > 9) &
        (te18 > 3) &
        (te36 > 10) &
        ( ((theta > 10) & (te6 > 3)) | ((theta <= 10) & (te6 > 10)) )
    )
    valid_sector = (sector6 != -9999)

    # -------------------------------------------------------------------------
    # Your reconstructed PID subset (same as original):
    #   electrons (11), negative pions (-211), negative kaons (-321), antiprotons (-2212)
    # -------------------------------------------------------------------------
    pid_subset = ( (pid == 11) | (pid == -211) | (pid == -321) | (pid == -2212) )

    # Common selection mask with reconstructed PID subset
    mask_common = (
        pid_subset &
        valid_sector &
        fid &
        (vz >= vz_min) & (vz <= vz_max) &
        (p > 4.9) &          # only above 4.9 GeV
        (nphe >= 2) &
        (e1 >= 0.06) &
        (e4 >= 0)
    )

    # -------------------------------------------------------------------------
    # Split by MC truth:
    #  - True electrons: mc_particle_id == 11
    #  - "False electrons": mc_particle_id in {-211, -321}
    # -------------------------------------------------------------------------
    mask_true_e  = mask_common & (mc_id == 11)
    mask_false_e = mask_common & ((mc_id == -211) | (mc_id == -321))

    # Fractions & sectors for each category
    frac_pcal_true  = e1[mask_true_e]  / p[mask_true_e]
    frac_ecin_true  = e4[mask_true_e]  / p[mask_true_e]
    sectors_true    = sector6[mask_true_e]

    frac_pcal_false = e1[mask_false_e] / p[mask_false_e]
    frac_ecin_false = e4[mask_false_e] / p[mask_false_e]
    sectors_false   = sector6[mask_false_e]

    # (Optional) quick stats to sanity check counts
    print(f"[INFO] Events after common cuts: {mask_common.sum()}")
    print(f"[INFO]   True electrons (mc_id==11): {mask_true_e.sum()}")
    print(f"[INFO]   False electrons (mc_id in -211,-321): {mask_false_e.sum()}")

    # -------------------------------------------------------------------------
    # Make the two figures
    # -------------------------------------------------------------------------
    title_true  = f"Diagonal Cut: PCal vs ECin Fractions - {label} (True electrons: mc_particle_id == 11)"
    out_true    = os.path.join(outdir, "diagonal_cut_MC_true_electrons.pdf")
    plot_by_sector(frac_pcal_true, frac_ecin_true, sectors_true, x_bins, y_bins, title_true, out_true)

    title_false = f"Diagonal Cut: PCal vs ECin Fractions - {label} (False electrons: mc_particle_id ∈ {{-211, -321}})"
    out_false   = os.path.join(outdir, "diagonal_cut_MC_false_electrons.pdf")
    plot_by_sector(frac_pcal_false, frac_ecin_false, sectors_false, x_bins, y_bins, title_false, out_false)

if __name__ == "__main__":
    main()