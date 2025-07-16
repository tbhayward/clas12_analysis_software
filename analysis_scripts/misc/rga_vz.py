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
    left = None
    if peak > 0:
        il = np.argmin(np.abs(hist[:peak] - threshold))
        left = centers[il]
    right = None
    if peak < len(hist) - 1:
        hr = hist[peak+1:]
        cr = centers[peak+1:]
        ir = np.argmin(np.abs(hr - threshold))
        right = cr[ir]
    return left, right

def plot_vz(data_list, labels, title, bins, lefts, rights, outpath):
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    colors = ["C0", "C1", "C2"]
    for data, label, c in zip(data_list, labels, colors):
        axs[0].hist(data, bins=bins, density=True,
                    histtype="step", color=c, label=label)
        axs[1].hist(data, bins=bins, density=True,
                    histtype="step", color=c, label=label)
    for left, right, c in zip(lefts, rights, colors):
        for ax in axs:
            ax.axvline(left,  color=c, linestyle="-",  alpha=0.5)
            ax.axvline(right, color=c, linestyle="--", alpha=0.5)
            ax.set_xlim(bins[0], bins[-1])
    axs[0].set_xlabel(r"$v_z$ (cm)")
    axs[0].set_ylabel("Normalized Counts")
    axs[0].set_title(f"{title} Vertex Distribution")
    axs[0].legend()
    axs[1].set_yscale("log")
    axs[1].set_xlabel(r"$v_z$ (cm)")
    axs[1].set_ylabel("Normalized Counts")
    axs[1].set_title(f"{title} Vertex Distribution (log scale)")
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)

def main():
    # RGA run files and labels
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

    # Loop over each period, apply appropriate fiducial cuts
    for fname, label in zip(files, labels):
        arr = uproot.open(fname)[tree_name].arrays([
            "particle_pid","particle_vz","track_sector_6",
            "cal_lv_1","cal_lw_1",
            "traj_edge_18","traj_edge_36","traj_edge_6",
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

        # always require lv1>9, lw1>9, te18>3, te36>10
        fid_basic = (lv1 > 9) & (lw1 > 9) & (te18 > 3) & (te36 > 10)
        # for Fa18 Out: require te6 > 3 for all theta
        if label == "Fa18 Out":
            fid = fid_basic & (te6 > 3)
        else:
            # for inbending periods: te6 cut depends on theta
            fid = fid_basic & (
                ((theta > 10) & (te6 > 3)) |
                ((theta <= 10) & (te6 > 10))
            )

        valid_sec = (sector != -9999)

        # negative: e⁻, π⁻, K⁻, p̄
        mask_neg = ((pid == 11)|(pid == -211)|(pid == -321)|(pid == -2212)) \
                   & valid_sec & fid
        # positive: e⁺, π⁺, K⁺, p
        mask_pos = ((pid == -11)|(pid == 211)|(pid == 321)|(pid == 2212)) \
                   & valid_sec & fid

        negative_vz.append(vz[mask_neg])
        positive_vz.append(vz[mask_pos])

    # compute thresholds per period
    neg_thr = [find_thresholds(v, bins, threshold) for v in negative_vz]
    pos_thr = [find_thresholds(v, bins, threshold) for v in positive_vz]
    neg_lefts,  neg_rights  = zip(*neg_thr)
    pos_lefts,  pos_rights  = zip(*pos_thr)

    # print out
    print("Negative-particle vz thresholds (density~0.02):")
    for lab, l, r in zip(labels, neg_lefts, neg_rights):
        print(f"  {lab}: left = {l:.3f}, right = {r:.3f}")
    print("Positive-particle vz thresholds (density~0.02):")
    for lab, l, r in zip(labels, pos_lefts, pos_rights):
        print(f"  {lab}: left = {l:.3f}, right = {r:.3f}")

    # output plots
    outdir = "output/rga_studies"
    os.makedirs(outdir, exist_ok=True)

    plot_vz(negative_vz, labels, "Negative Particles", bins,
            neg_lefts, neg_rights,
            os.path.join(outdir, "negative_particles_vz.pdf"))

    plot_vz(positive_vz, labels, "Positive Particles", bins,
            pos_lefts, pos_rights,
            os.path.join(outdir, "positive_particles_vz.pdf"))

if __name__ == "__main__":
    main()
