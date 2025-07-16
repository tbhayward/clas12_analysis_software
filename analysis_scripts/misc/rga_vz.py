#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def find_thresholds(data, bins, threshold):
    """
    Find left/right bin centers where the density histogram
    crosses the given threshold around the peak.
    """
    hist, edges = np.histogram(data, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak = np.argmax(hist)
    # left crossing
    left = None
    if peak > 0:
        il = np.argmin(np.abs(hist[:peak] - threshold))
        left = centers[il]
    # right crossing
    right = None
    if peak < len(hist)-1:
        hr = hist[peak+1:]
        cr = centers[peak+1:]
        ir = np.argmin(np.abs(hr - threshold))
        right = cr[ir]
    return left, right

def plot_vz(data_list, labels, name, bins, in_thr, out_thr, outpath):
    """
    Draw vertex distributions for one track charge:
    solid red = in‑bending thresholds, dashed red = out‑bending.
    """
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    colors = ["C0", "C1", "C2"]

    # overlay the three periods
    for data, lab, c in zip(data_list, labels, colors):
        axs[0].hist(data, bins=bins, density=True,
                    histtype="step", color=c, label=lab)
        axs[1].hist(data, bins=bins, density=True,
                    histtype="step", color=c, label=lab)

    (l_in,  r_in),  (l_out,  r_out)  = in_thr,  out_thr

    for ax in axs:
        # in‑bending cuts: solid red
        ax.axvline(l_in,  color='red', linestyle='-', alpha=0.7)
        ax.axvline(r_in,  color='red', linestyle='-', alpha=0.7)
        # out‑bending cuts: dashed red
        ax.axvline(l_out, color='red', linestyle='--', alpha=0.7)
        ax.axvline(r_out, color='red', linestyle='--', alpha=0.7)
        ax.set_xlim(bins[0], bins[-1])

    # labels and titles
    axs[0].set_xlabel(r"$v_z$ (cm)")
    axs[0].set_ylabel("Normalized Counts")
    axs[0].set_title(f"{name} Vertex Distribution")
    axs[0].legend()

    axs[1].set_yscale("log")
    axs[1].set_xlabel(r"$v_z$ (cm)")
    axs[1].set_ylabel("Normalized Counts")
    axs[1].set_title(f"{name} Vertex Distribution (log scale)")
    axs[1].legend()

    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)

def main():
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
    threshold = 0.02
    bins = np.linspace(-15, 15, 100)

    negative_vz = []
    positive_vz = []

    for fname in files:
        tree = uproot.open(fname)[tree_name]
        arr = tree.arrays([
            "particle_pid", "particle_vz", "track_sector_6",
            "cal_lv_1", "cal_lw_1",
            "traj_edge_18", "traj_edge_36", "traj_edge_6",
            "theta"
        ], library="np")

        pid    = arr["particle_pid"]
        vz     = arr["particle_vz"]
        sector = arr["track_sector_6"]
        lv1    = arr["cal_lv_1"]
        lw1    = arr["cal_lw_1"]
        te18   = arr["traj_edge_18"]
        te36   = arr["traj_edge_36"]
        te6    = arr["traj_edge_6"]
        theta  = arr["theta"]

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
        valid_sec = (sector != -9999)

        # negative tracks: e⁻, π⁻, K⁻, p̄
        mask_neg = (
            ((pid == 11)   | (pid == -211) |
             (pid == -321) | (pid == -2212)) &
            valid_sec & fid
        )
        # positive tracks: e⁺, π⁺, K⁺, p
        mask_pos = (
            ((pid == -11)  | (pid == 211)  |
             (pid == 321)  | (pid == 2212)) &
            valid_sec & fid
        )

        negative_vz.append(vz[mask_neg])
        positive_vz.append(vz[mask_pos])

    # combine the two in‑bending periods (Fa18 Inb + Sp19 Inb)
    neg_inb = np.concatenate([negative_vz[0], negative_vz[2]])
    pos_inb = np.concatenate([positive_vz[0], positive_vz[2]])

    # compute thresholds
    thr_neg_in  = find_thresholds(neg_inb,   bins, threshold)
    thr_neg_out = find_thresholds(negative_vz[1], bins, threshold)
    thr_pos_in  = find_thresholds(pos_inb,   bins, threshold)
    thr_pos_out = find_thresholds(positive_vz[1], bins, threshold)

    # print out for both charges
    print("Negative particles thresholds (density~0.02):")
    print(f"  Fa18 Inb + Sp19 Inb: {thr_neg_in[0]:.3f}, {thr_neg_in[1]:.3f}")
    print(f"  Fa18 Out         : {thr_neg_out[0]:.3f}, {thr_neg_out[1]:.3f}\n")

    print("Positive particles thresholds (density~0.02):")
    print(f"  Fa18 Inb + Sp19 Inb: {thr_pos_in[0]:.3f}, {thr_pos_in[1]:.3f}")
    print(f"  Fa18 Out         : {thr_pos_out[0]:.3f}, {thr_pos_out[1]:.3f}")

    # ensure output dir
    outdir = "output/rga_studies"
    os.makedirs(outdir, exist_ok=True)

    # plot negative and positive
    plot_vz(
        negative_vz, labels, "Negative Particles",
        bins, thr_neg_in, thr_neg_out,
        os.path.join(outdir, "negative_particles_vz_rga.pdf")
    )
    plot_vz(
        positive_vz, labels, "Positive Particles",
        bins, thr_pos_in, thr_pos_out,
        os.path.join(outdir, "positive_particles_vz_rga.pdf")
    )

if __name__ == "__main__":
    main()
