#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def find_thresholds(data, bins, threshold):
    """
    Find the left and right bin centers where the density histogram
    is closest to the given threshold, on either side of the peak.
    """
    hist, edges = np.histogram(data, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak_idx = np.argmax(hist)
    # Left side
    left = None
    if peak_idx > 0:
        left_idx = np.argmin(np.abs(hist[:peak_idx] - threshold))
        left = centers[left_idx]
    # Right side
    right = None
    if peak_idx < len(hist) - 1:
        right_hist = hist[peak_idx+1:]
        right_centers = centers[peak_idx+1:]
        right_idx = np.argmin(np.abs(right_hist - threshold))
        right = right_centers[right_idx]
    return left, right

def main():
    # Paths to your ROOT files
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

    # Collect vz arrays
    negative_vz = []
    positive_vz = []

    for fname in files:
        with uproot.open(fname)[tree_name] as tree:
            pid     = tree["particle_pid"].array(library="np")
            vz      = tree["particle_vz"].array(library="np")
            sector  = tree["track_sector_6"].array(library="np")

            # negative particles: electron (11), pi- (-211), K- (-321)
            mask_neg = (
                ((pid == 11) | (pid == -211) | (pid == -321))
                & (sector != -9999)
            )
            # positive particles: proton (2212), positron (-11), pi+ (211), K+ (321)
            mask_pos = (
                ((pid == 2212) | (pid == -11) | (pid == 211) | (pid == 321))
                & (sector != -9999)
            )

            negative_vz.append(vz[mask_neg])
            positive_vz.append(vz[mask_pos])
    #endfor

    # Combine Fa22 and Sp23 into a "secret" fourth distribution
    neg_comb = np.concatenate([negative_vz[1], negative_vz[2]])
    pos_comb = np.concatenate([positive_vz[1], positive_vz[2]])

    # Histogram settings: -15 to 15
    bins = np.linspace(-15, 15, 100)

    # Calculate thresholds for Su22 and combined Fa22+Sp23
    n_su_l, n_su_r = find_thresholds(negative_vz[0], bins, threshold)
    n_c_l, n_c_r   = find_thresholds(neg_comb,      bins, threshold)
    p_su_l, p_su_r = find_thresholds(positive_vz[0], bins, threshold)
    p_c_l, p_c_r   = find_thresholds(pos_comb,      bins, threshold)

    print("Threshold positions (density ~ 0.02):")
    print(f"  Negative Su22:     left = {n_su_l:.3f}, right = {n_su_r:.3f}")
    print(f"  Negative Fa22+Sp23:left = {n_c_l:.3f}, right = {n_c_r:.3f}")
    print(f"  Positive Su22:     left = {p_su_l:.3f}, right = {p_su_r:.3f}")
    print(f"  Positive Fa22+Sp23:left = {p_c_l:.3f}, right = {p_c_r:.3f}")

    # Calculate percentage outside thresholds via histogram integration
    def pct_outside(data, left, right):
        counts, edges = np.histogram(data, bins=bins, density=False)
        centers = 0.5 * (edges[:-1] + edges[1:])
        total = counts.sum()
        iL = np.searchsorted(centers, left)
        iR = np.searchsorted(centers, right)
        outside = counts[:iL].sum() + counts[iR+1:].sum()
        return 100 * outside / total

    pct_neg_su = pct_outside(negative_vz[0], n_su_l, n_su_r)
    pct_neg_c  = pct_outside(neg_comb,      n_c_l,  n_c_r)
    pct_pos_su = pct_outside(positive_vz[0], p_su_l, p_su_r)
    pct_pos_c  = pct_outside(pos_comb,      p_c_l,  p_c_r)

    print("\nPercentage outside thresholds:")
    print(f"  Negative Su22:     {pct_neg_su:.2f}%")
    print(f"  Negative Fa22+Sp23:{pct_neg_c:.2f}%")
    print(f"  Positive Su22:     {pct_pos_su:.2f}%")
    print(f"  Positive Fa22+Sp23:{pct_pos_c:.2f}%")

    # Ensure output directory exists
    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    def plot_vz(data_list, combined, pname, left, right, c_left, c_right):
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        colors = ["C0", "C1", "C2"]
        for data, label, color in zip(data_list, labels, colors):
            axes[0].hist(data, bins=bins, density=True,
                         histtype="step", color=color, label=label)
            axes[1].hist(data, bins=bins, density=True,
                         histtype="step", color=color, label=label)
        #endfor

        # Plot original hard-coded lines (commented out)
        for ax in axes:
            # ax.axvline(-7,  color='red', linestyle='-',  alpha=0.5)
            # ax.axvline(0,   color='red', linestyle='-',  alpha=0.5)
            # ax.axvline(-6,  color='red', linestyle='--', alpha=0.5)
            # ax.axvline(1,   color='red', linestyle='--', alpha=0.5)

            # calculated Su22 thresholds
            ax.axvline(left,   color='red', linestyle='-',  alpha=0.5)
            ax.axvline(right,  color='red', linestyle='-',  alpha=0.5)
            # calculated Fa22+Sp23 thresholds
            ax.axvline(c_left, color='red', linestyle='--', alpha=0.5)
            ax.axvline(c_right,color='red', linestyle='--', alpha=0.5)

            ax.set_xlim(-15, 15)
        #endfor

        axes[0].set_xlabel(r"$v_z$ (cm)")
        axes[0].set_ylabel("Normalized Counts")
        axes[0].set_title(f"{pname} Vertex Distribution")
        axes[0].legend()

        axes[1].set_yscale("log")
        axes[1].set_xlabel(r"$v_z$ (cm)")
        axes[1].set_ylabel("Normalized Counts")
        axes[1].set_title(f"{pname} Vertex Distribution (log scale)")
        axes[1].legend()

        fig.tight_layout()
        fig.savefig(f"{outdir}/{pname.lower().replace(' ', '_')}_vz.pdf")
        plt.close(fig)

    # Plot negative and positive
    plot_vz(negative_vz, neg_comb, "Negative Particles",
            n_su_l, n_su_r, n_c_l, n_c_r)
    plot_vz(positive_vz, pos_comb, "Positive Particles",
            p_su_l, p_su_r, p_c_l, p_c_r)

if __name__ == "__main__":
    main()