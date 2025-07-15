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
    # left
    left = None
    if peak > 0:
        il = np.argmin(np.abs(hist[:peak] - threshold))
        left = centers[il]
    # right
    right = None
    if peak < len(hist)-1:
        hr = hist[peak+1:]
        cr = centers[peak+1:]
        ir = np.argmin(np.abs(hr - threshold))
        right = cr[ir]
    return left, right

def main():
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
    threshold = 0.02
    bins = np.linspace(-15, 15, 100)

    negative_vz = []
    positive_vz = []

    for fname, label in zip(files, labels):
        tree = uproot.open(fname)[tree_name]
        arr = tree.arrays([
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

        # fiducial
        fid = (
            (lv1 > 9) &
            (lw1 > 9) &
            (te18 > 3) &
            (te36 > 10) &
            (((theta > 10) & (te6 > 3)) |
             ((theta <= 10) & (te6 > 10)))
        )
        valid_sec = (sector != -9999)

        # negative: e⁻ (11), π⁻(-211), K⁻(-321), p̄(-2212)
        mask_neg = (
            ((pid==11)|(pid==-211)|(pid==-321)|(pid==-2212)) &
            valid_sec & fid
        )
        # positive: e⁺(-11), π⁺(211), K⁺(321), p(2212)
        mask_pos = (
            ((pid==-11)|(pid==211)|(pid==321)|(pid==2212)) &
            valid_sec & fid
        )

        negative_vz.append(vz[mask_neg])
        positive_vz.append(vz[mask_pos])

    # combine Fa22+Sp23
    neg_comb = np.concatenate([negative_vz[1], negative_vz[2]])
    pos_comb = np.concatenate([positive_vz[1], positive_vz[2]])

    # compute thresholds
    n_su_l,n_su_r = find_thresholds(negative_vz[0], bins, threshold)
    n_c_l,n_c_r   = find_thresholds(neg_comb,      bins, threshold)
    p_su_l,p_su_r = find_thresholds(positive_vz[0], bins, threshold)
    p_c_l,p_c_r   = find_thresholds(pos_comb,      bins, threshold)

    print("Thresholds (density~0.02):")
    print(f"  Negative Su22:      {n_su_l:.3f}, {n_su_r:.3f}")
    print(f"  Negative Fa22+Sp23: {n_c_l:.3f}, {n_c_r:.3f}")
    print(f"  Positive Su22:      {p_su_l:.3f}, {p_su_r:.3f}")
    print(f"  Positive Fa22+Sp23: {p_c_l:.3f}, {p_c_r:.3f}")

    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    def plot_vz(data_list, combined, name, sl, sr, cl, cr):
        fig, axs = plt.subplots(1,2,figsize=(12,6))
        colors = ["C0","C1","C2"]
        for d, lab, c in zip(data_list, labels, colors):
            axs[0].hist(d, bins=bins, density=True,
                        histtype="step", color=c, label=lab)
            axs[1].hist(d, bins=bins, density=True,
                        histtype="step", color=c, label=lab)
        for ax in axs:
            ax.axvline(sl, color="red", linestyle="-",  alpha=0.5)
            ax.axvline(sr, color="red", linestyle="-",  alpha=0.5)
            ax.axvline(cl, color="red", linestyle="--", alpha=0.5)
            ax.axvline(cr, color="red", linestyle="--", alpha=0.5)
            ax.set_xlim(-15,15)
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
        fig.savefig(f"{outdir}/{name.lower().replace(' ','_')}_vz.pdf")
        plt.close(fig)

    plot_vz(negative_vz, neg_comb,
            "Negative Particles", n_su_l,n_su_r, n_c_l,n_c_r)
    plot_vz(positive_vz, pos_comb,
            "Positive Particles", p_su_l,p_su_r, p_c_l,p_c_r)

if __name__ == "__main__":
    main()