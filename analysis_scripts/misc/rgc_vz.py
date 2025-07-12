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
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_fa22_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_sp23_inb_calibration.root",
    ]
    labels = ["Su22", "Fa22", "Sp23"]
    tree_name = "PhysicsEvents"
    threshold = 0.03

    # Collect vz arrays
    electron_vz = []
    proton_vz   = []

    for fname in files:
        with uproot.open(fname)[tree_name] as tree:
            pid    = tree["particle_pid"].array(library="np")
            vz     = tree["particle_vz"].array(library="np")
            sector = tree["track_sector_6"].array(library="np")

            mask_e = (pid == 11)   & (sector != -9999)
            mask_p = (pid == 2212) & (sector != -9999)

            electron_vz.append(vz[mask_e])
            proton_vz.append(vz[mask_p])
    #endfor

    # Combine Fa22 and Sp23 into a "secret" fourth distribution
    electron_comb = np.concatenate([electron_vz[1], electron_vz[2]])
    proton_comb   = np.concatenate([proton_vz[1],   proton_vz[2]])

    # Histogram settings: -15 to 15
    bins = np.linspace(-15, 15, 100)

    # Calculate thresholds for Su22 and combined Fa22+Sp23
    e_su_left,  e_su_right  = find_thresholds(electron_vz[0], bins, threshold)
    e_c_left,   e_c_right   = find_thresholds(electron_comb,   bins, threshold)
    p_su_left,  p_su_right  = find_thresholds(proton_vz[0],   bins, threshold)
    p_c_left,   p_c_right   = find_thresholds(proton_comb,     bins, threshold)

    print("Electron thresholds:")
    print(f"  Su22 left = {e_su_left:.3f} cm, right = {e_su_right:.3f} cm")
    print(f"  Fa22+Sp23 combined left = {e_c_left:.3f} cm, right = {e_c_right:.3f} cm")
    print("Proton thresholds:")
    print(f"  Su22 left = {p_su_left:.3f} cm, right = {p_su_right:.3f} cm")
    print(f"  Fa22+Sp23 combined left = {p_c_left:.3f} cm, right = {p_c_right:.3f} cm")

    # Ensure output directory exists
    outdir = "output/rgc_studies"
    os.makedirs(outdir, exist_ok=True)

    def plot_vz(data_list, combined, pname, su_left, su_right, c_left, c_right):
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        colors = ["C0", "C1", "C2"]
        for data, label, color in zip(data_list, labels, colors):
            axes[0].hist(data, bins=bins, density=True,
                         histtype="step", color=color, label=label)
            axes[1].hist(data, bins=bins, density=True,
                         histtype="step", color=color, label=label)
        #endfor

        for ax in axes:
            # original hard-coded lines (commented out)
            # ax.axvline(-7,  color='red', linestyle='-',  alpha=0.25)
            # ax.axvline(-0.5, color='red', linestyle='-',  alpha=0.25)
            # ax.axvline(-6,  color='red', linestyle='--', alpha=0.25)
            # ax.axvline(0.5,  color='red', linestyle='--', alpha=0.25)

            # calculated Su22 thresholds
            ax.axvline(su_left,  color='red', linestyle='-',  alpha=0.25)
            ax.axvline(su_right, color='red', linestyle='-',  alpha=0.25)
            # calculated Fa22+Sp23 thresholds
            ax.axvline(c_left,   color='red', linestyle='--', alpha=0.25)
            ax.axvline(c_right,  color='red', linestyle='--', alpha=0.25)

            ax.set_xlim(-15, 15)
        #endfor

        axes[0].set_xlabel(r"$v_{z}$ (cm)")
        axes[0].set_ylabel("Normalized Counts")
        axes[0].set_title(f"{pname} Vertex Distribution")
        axes[0].legend()

        axes[1].set_yscale("log")
        axes[1].set_xlabel(r"$v_{z}$ (cm)")
        axes[1].set_ylabel("Normalized Counts")
        axes[1].set_title(f"{pname} Vertex Distribution (log scale)")
        axes[1].legend()

        fig.tight_layout()
        fig.savefig(f"{outdir}/{pname.lower()}_vz.pdf")
        plt.close(fig)

    # Plot electrons and protons
    plot_vz(electron_vz, electron_comb, "Electron",
            e_su_left, e_su_right, e_c_left, e_c_right)
    plot_vz(proton_vz,   proton_comb,   "Proton",
            p_su_left, p_su_right, p_c_left, p_c_right)

if __name__ == "__main__":
    main()
#endif