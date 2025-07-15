#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def compute_means_sigmas(p_vals, sf_vals, p_bins):
    """
    Given arrays of momenta and sampling fractions,
    compute the mean and sigma of sf in each p_bin.
    Returns bin centers, means, sigmas.
    """
    bin_indices = np.digitize(p_vals, p_bins) - 1
    centers = 0.5 * (p_bins[:-1] + p_bins[1:])
    means = np.zeros_like(centers)
    sigmas = np.zeros_like(centers)
    for i in range(len(centers)):
        mask = (bin_indices == i)
        if np.count_nonzero(mask) > 2:
            data = sf_vals[mask]
            means[i] = data.mean()
            sigmas[i] = data.std()
        else:
            means[i] = np.nan
            sigmas[i] = np.nan
    return centers, means, sigmas

def make_sampling_fraction_plot(filename, label, vz_cut, outdir):
    """
    Create a 2x3 sampling fraction vs p plot,
    applying fiducial, vertex, p>2, HTCC nphe>=2,
    PCal E>=0.15, diagonal cut, then fit SF and
    print mean/sigma at specific p values.
    """
    # Load tree and branches
    tree = uproot.open(filename)["PhysicsEvents"]
    arr = tree.arrays([
        "particle_pid", "particle_vz", "cal_sector",
        "p", "cc_nphe_15",
        "cal_energy_1", "cal_energy_4", "cal_energy_7",
        "cal_lv_1", "cal_lw_1",
        "traj_edge_18", "traj_edge_36", "traj_edge_6",
        "theta"
    ], library="np")

    pid     = arr["particle_pid"]
    vz      = arr["particle_vz"]
    sector  = arr["cal_sector"]
    p       = arr["p"]
    nphe    = arr["cc_nphe_15"]
    e1      = arr["cal_energy_1"]
    e4      = arr["cal_energy_4"]
    e7      = arr["cal_energy_7"]
    lv1     = arr["cal_lv_1"]
    lw1     = arr["cal_lw_1"]
    te18    = arr["traj_edge_18"]
    te36    = arr["traj_edge_36"]
    te6     = arr["traj_edge_6"]
    theta   = arr["theta"]

    # Fiducial and sector cuts
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
    valid_sector = (sector != -9999)

    # Base mask: electron + all prior cuts
    base_mask = (
        (pid == 11) &
        valid_sector &
        fid &
        (vz >= vz_cut[0]) &
        (vz <= vz_cut[1]) &
        (p > 2.0) &
        (nphe >= 2) &
        (e1 >= 0.15) &
        (e4 >= 0)
    )

    # Fractions for diagonal cut
    p_base   = p[base_mask]
    e1_base  = e1[base_mask]
    e4_base  = e4[base_mask]
    frac_pcal = e1_base / p_base
    frac_ecin = e4_base / p_base

    # Diagonal cut: y >= -0.625*x + 0.15
    mask_diag = frac_ecin >= (-0.625 * frac_pcal + 0.15)

    # Final indices
    idx_base = np.nonzero(base_mask)[0]
    keep_idx = idx_base[mask_diag]

    # Prepare final arrays
    p_vals = p[keep_idx]
    sf_vals = (e1[keep_idx] + e4[keep_idx] + e7[keep_idx]) / p_vals
    sectors = sector[keep_idx]

    # Logging
    n_total = len(pid)
    n_after = len(keep_idx)
    print(f"{label}: total = {n_total:,}, after all cuts = {n_after:,}")

    # Prepare bins and evaluation points
    p_bins = np.linspace(2.0, 8.0, 40)
    sf_range = (0.10, 0.40)
    p_eval = np.array([2, 3, 4, 5, 6, 7, 8])

    # 2x3 canvas
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)

    for sec in range(1, 7):
        ax = axes.flat[sec-1]
        sel = (sectors == sec)
        p_sec = p_vals[sel]
        sf_sec = sf_vals[sel]

        # 2D histogram
        h = ax.hist2d(
            p_sec, sf_sec,
            bins=[p_bins, np.linspace(*sf_range, 80)],
            cmap="jet",
            norm=matplotlib.colors.LogNorm()
        )

        # Compute profile
        centers, means, sigmas = compute_means_sigmas(p_sec, sf_sec, p_bins)

        # Print mean & sigma at each p_eval
        mean_at_p = np.interp(p_eval, centers, means)
        sig_at_p  = np.interp(p_eval, centers, sigmas)
        for pv, mv, sv in zip(p_eval, mean_at_p, sig_at_p):
            print(f"{label} Sector {sec}, p={pv} GeV: mean={mv:.4f}, sigma={sv:.4f}")

        # Fit mean & sigma curves
        valid = ~np.isnan(means)
        coef_m = np.polyfit(centers[valid], means[valid], 2)
        coef_s = np.polyfit(centers[valid], sigmas[valid], 2)
        pm = np.poly1d(coef_m)
        ps = np.poly1d(coef_s)

        # Overlay mean±3σ
        p_fit = np.linspace(2.0, 8.0, 200)
        m_fit = pm(p_fit)
        s_fit = ps(p_fit)
        ax.plot(p_fit, m_fit, color='red', linestyle='-',  linewidth=2)
        ax.plot(p_fit, m_fit + 3*s_fit, color='red', linestyle='--', linewidth=2)
        ax.plot(p_fit, m_fit - 3*s_fit, color='red', linestyle='--', linewidth=2)

        ax.set_xlim(2.0, 8.0)
        ax.set_ylim(*sf_range)
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("p (GeV)")
        ax.set_ylabel("Sampling Fraction")

    # Shared colorbar
    fig.colorbar(h[3], ax=axes.ravel().tolist(), shrink=0.9).set_label("Counts (log scale)")

    # Save
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"electron_final_sampling_fraction_{label}.pdf")
    fig.suptitle(f"Final Sampling Fraction - {label}", fontsize=16)
    fig.savefig(outpath)
    plt.close(fig)

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
    vz_cuts = {
        "Su22": (-7.576, 0.303),
        "Fa22": (-5.758, 1.515),
        "Sp23": (-5.758, 1.515),
    }
    outdir = "output/rgc_studies"

    for fname, label in zip(files, labels):
        make_sampling_fraction_plot(fname, label, vz_cuts[label], outdir)

if __name__ == "__main__":
    main()