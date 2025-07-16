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
    solid red = in-bending thresholds, dashed red = out-bending.
    """
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    colors = ["C0", "C1", "C2"]

    for data, lab, c in zip(data_list, labels, colors):
        axs[0].hist(data, bins=bins, density=True,
                    histtype="step", color=c, label=lab)
        axs[1].hist(data, bins=bins, density=True,
                    histtype="step", color=c, label=lab)

    (l_in,  r_in),  (l_out,  r_out)  = in_thr,  out_thr

    for ax in axs:
        # in-bending cuts: solid red
        ax.axvline(l_in,  color='red', linestyle='-', alpha=0.7)
        ax.axvline(r_in,  color='red', linestyle='-', alpha=0.7)
        # out-bending cuts: dashed red
        ax.axvline(l_out, color='red', linestyle='--', alpha=0.7)
        ax.axvline(r_out, color='red', linestyle='--', alpha=0.7)
        ax.set_xlim(bins[0], bins[-1])

    axs[0].set_xlabel(r"$v_z$ (cm)")
    axs[0].set_ylabel("Normalized Counts")
    axs[0].set_title(f"{name}")
    axs[0].legend()

    axs[1].set_yscale("log")
    axs[1].set_xlabel(r"$v_z$ (cm)")
    axs[1].set_ylabel("Normalized Counts")
    axs[1].set_title(f"{name} (log scale)")
    axs[1].legend()

    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)

def main():
    threshold = 0.02
    bins = np.linspace(-15, 15, 100)
    outdir = "output/rga_studies"
    os.makedirs(outdir, exist_ok=True)

    # --- Group 1: Fa18 Inb, Fa18 Out, Sp19 Inb ---
    files1 = [
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_fa18_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_fa18_out_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_sp19_inb_calibration.root",
    ]
    labels1 = ["Fa18 Inb", "Fa18 Out", "Sp19 Inb"]
    vz_cuts1 = {
        "Fa18 Inb": (-6.364, 1.515),
        "Fa18 Out": (-7.879, 0.303),
        "Sp19 Inb": (-6.364, 1.515),
    }

    neg1, pos1 = [], []
    for fname, lab in zip(files1, labels1):
        tree = uproot.open(fname)["PhysicsEvents"]
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

        fid = (
            (lv1 > 9)&(lw1 > 9)&(te18>3)&(te36>10)&
            (((theta>10)&(te6>3))|((theta<=10)&(te6>10)))
        )
        vs = (sector != -9999)

        neg1.append(vz[((pid==11)|(pid==-211)|(pid==-321)|(pid==-2212)) & vs & fid])
        pos1.append(vz[((pid==-11)|(pid==211)|(pid==321)|(pid==2212)) & vs & fid])

    # combine in-bending (Fa18 Inb + Sp19 Inb)
    neg1_in = np.concatenate([neg1[0], neg1[2]])
    pos1_in = np.concatenate([pos1[0], pos1[2]])
    neg1_out = neg1[1]
    pos1_out = pos1[1]

    thr_neg1_in  = find_thresholds(neg1_in,  bins, threshold)
    thr_neg1_out = find_thresholds(neg1_out, bins, threshold)
    thr_pos1_in  = find_thresholds(pos1_in,  bins, threshold)
    thr_pos1_out = find_thresholds(pos1_out, bins, threshold)

    print("Group 1 — Negative thresholds:")
    print(f"  In-bending (Fa18 Inb + Sp19 Inb): {thr_neg1_in}")
    print(f"  Out-bending (Fa18 Out)         : {thr_neg1_out}\n")

    print("Group 1 — Positive thresholds:")
    print(f"  In-bending (Fa18 Inb + Sp19 Inb): {thr_pos1_in}")
    print(f"  Out-bending (Fa18 Out)         : {thr_pos1_out}\n")

    # plot group 1
    plot_vz(neg1, labels1, "Negative Particles (Fa18 Inb, Fa18 Out, Sp19 Inb)",
            bins, thr_neg1_in, thr_neg1_out,
            os.path.join(outdir, "negative_particles_vz_group1.pdf"))

    plot_vz(pos1, labels1, "Positive Particles (Fa18 Inb, Fa18 Out, Sp19 Inb)",
            bins, thr_pos1_in, thr_pos1_out,
            os.path.join(outdir, "positive_particles_vz_group1.pdf"))


    # --- Group 2: Sp18 Inb, Sp18 Out ---
    files2 = [
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_sp18_inb_calibration.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/"
        "nSidis_rga_sp18_out_calibration.root",
    ]
    labels2 = ["Sp18 Inb", "Sp18 Out"]
    vz_cuts2 = {
        "Sp18 Inb": (-6.364, 1.515),
        "Sp18 Out": (-7.879, 0.303),
    }

    neg2, pos2 = [], []
    for fname, lab in zip(files2, labels2):
        tree = uproot.open(fname)["PhysicsEvents"]
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

        fid = (
            (lv1 > 9)&(lw1 > 9)&(te18>3)&(te36>10)&
            (((theta>10)&(te6>3))|((theta<=10)&(te6>10)))
        )
        vs = (sector != -9999)

        neg2.append(vz[((pid==11)|(pid==-211)|(pid==-321)|(pid==-2212)) & vs & fid])
        pos2.append(vz[((pid==-11)|(pid==211)|(pid==321)|(pid==2212)) & vs & fid])

    thr_neg2_in  = find_thresholds(neg2[0], bins, threshold)
    thr_neg2_out = find_thresholds(neg2[1], bins, threshold)
    thr_pos2_in  = find_thresholds(pos2[0], bins, threshold)
    thr_pos2_out = find_thresholds(pos2[1], bins, threshold)

    print("Group 2 — Negative thresholds:")
    print(f"  In-bending (Sp18 Inb): {thr_neg2_in}")
    print(f"  Out-bending (Sp18 Out): {thr_neg2_out}\n")

    print("Group 2 — Positive thresholds:")
    print(f"  In-bending (Sp18 Inb): {thr_pos2_in}")
    print(f"  Out-bending (Sp18 Out): {thr_pos2_out}\n")

    # plot group 2
    plot_vz(neg2, labels2, "Negative Particles (Sp18 Inb, Sp18 Out)",
            bins, thr_neg2_in, thr_neg2_out,
            os.path.join(outdir, "negative_particles_vz_group2.pdf"))

    plot_vz(pos2, labels2, "Positive Particles (Sp18 Inb, Sp18 Out)",
            bins, thr_pos2_in, thr_pos2_out,
            os.path.join(outdir, "positive_particles_vz_group2.pdf"))

if __name__ == "__main__":
    main()