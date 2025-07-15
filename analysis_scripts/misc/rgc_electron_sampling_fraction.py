#!/usr/bin/env python3

import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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
    2×3 sampling-fraction vs p:
    apply fiducial, vertex, p>2, nphe>=2,
    e1>=0.15, e4>=0, e7>=0, diagonal cut,
    then profile & fit and print Java-style cut lines.
    """
    # load entire tree into numpy arrays
    tree = uproot.open(filename)["PhysicsEvents"]
    arr = tree.arrays([
        "particle_pid", "particle_vz", "cal_sector",
        "p", "cc_nphe_15",
        "cal_energy_1", "cal_energy_4", "cal_energy_7",
        "cal_lv_1", "cal_lw_1",
        "traj_edge_18", "traj_edge_36", "traj_edge_6",
        "theta"
    ], library="np")

    pid    = arr["particle_pid"]
    vz     = arr["particle_vz"]
    sector = arr["cal_sector"]
    p      = arr["p"]
    nphe   = arr["cc_nphe_15"]
    e1     = arr["cal_energy_1"]
    e4     = arr["cal_energy_4"]
    e7     = arr["cal_energy_7"]
    lv1    = arr["cal_lv_1"]
    lw1    = arr["cal_lw_1"]
    te18   = arr["traj_edge_18"]
    te36   = arr["traj_edge_36"]
    te6    = arr["traj_edge_6"]
    theta  = arr["theta"]

    # fiducial + valid sector mask
    fid = (
        (lv1 > 9) & (lw1 > 9) &
        (te18 > 3) & (te36 > 10) &
        (
            ((theta > 10) & (te6 > 3)) |
            ((theta <= 10) & (te6 > 10))
        )
    )
    good_sector = (sector != -9999)

    # base cuts
    base_mask = (
        (pid == 11) &
        good_sector &
        fid &
        (vz >= vz_cut[0]) & (vz <= vz_cut[1]) &
        (p  > 2.0) &
        (nphe >= 2) &
        (e1 >= 0.15) &
        (e4 >= 0)      &
        (e7 >= 0)
    )

    # diagonal cut on e4/p vs e1/p
    p_base    = p[base_mask]
    frac_pcal = e1[base_mask] / p_base
    frac_ecin = e4[base_mask] / p_base
    diag_mask = frac_ecin >= (-0.625 * frac_pcal + 0.15)

    idx0  = np.nonzero(base_mask)[0]
    keep  = idx0[diag_mask]

    p_vals  = p[keep]
    sf_vals = (e1[keep] + e4[keep] + e7[keep]) / p_vals
    secs    = sector[keep]

    # prepare binning
    p_bins   = np.linspace(2.0, 8.0, 40)
    sf_range = (0.10, 0.40)

    fig, axes = plt.subplots(2,3, figsize=(15,10), constrained_layout=True)

    for sec in range(1,7):
        ax    = axes.flat[sec-1]
        mask  = (secs == sec)
        p_sec = p_vals[mask]
        sf_sec= sf_vals[mask]

        # 2D histogram
        counts2d, xedges, yedges, im = ax.hist2d(
            p_sec, sf_sec,
            bins=[p_bins, np.linspace(*sf_range,80)],
            cmap="jet", norm=LogNorm()
        )

        # profile calculation
        centers, means, sigmas = compute_means_sigmas(p_sec, sf_sec, p_bins)

        # fit quadratic to mean and sigma
        good = ~np.isnan(means)
        cm = np.polyfit(centers[good], means[good], 2)
        cs = np.polyfit(centers[good], sigmas[good],2)

        # extract coefficients (cm[0]*p^2 + cm[1]*p + cm[2])
        c2_m, c1_m, c0_m = cm
        c2_s, c1_s, c0_s = cs

        a_minus = c0_m - 3*c0_s
        b_minus = c1_m - 3*c1_s
        c_minus = c2_m - 3*c2_s
        a_plus  = c0_m + 3*c0_s
        b_plus  = c1_m + 3*c1_s
        c_plus  = c2_m + 3*c2_s

        # Java‐compatible output
        print(
            f"Sector {sec}: sf > ({a_minus:.6f} + {b_minus:.6f}*p + {c_minus:.6f}*p*p) "
            f"&& sf < ({a_plus:.6f} + {b_plus:.6f}*p + {c_plus:.6f}*p*p);"
        )

        # overlay fits
        p_fit = np.linspace(2.0,8.0,200)
        pm   = np.poly1d(cm)
        ps   = np.poly1d(cs)
        m_fit = pm(p_fit)
        s_fit = ps(p_fit)

        ax.plot(p_fit, m_fit,      'r-',  lw=2)
        ax.plot(p_fit, m_fit+3*s_fit, 'r--', lw=2)
        ax.plot(p_fit, m_fit-3*s_fit, 'r--', lw=2)

        ax.set_xlim(2.0,8.0)
        ax.set_ylim(*sf_range)
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("p (GeV)")
        ax.set_ylabel("Sampling Fraction")

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.9)\
       .set_label("Counts (log scale)")

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"electron_final_sampling_fraction_{label}.pdf")
    fig.suptitle(f"Final Sampling Fraction – {label}", fontsize=16)
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