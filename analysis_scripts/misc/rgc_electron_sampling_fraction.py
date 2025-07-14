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
    Create a 2x3 sampling fraction plot for one run,
    enforcing PCal energy > 0.07 GeV, and print fit cuts.
    """
    # open tree
    tree = uproot.open(filename)["PhysicsEvents"]
    # read branches
    pid    = tree["particle_pid"].array(library="np")
    vz     = tree["particle_vz"].array(library="np")
    sector = tree["cal_sector"].array(library="np")
    p      = tree["p"].array(library="np")
    e1     = tree["cal_energy_1"].array(library="np")
    e4     = tree["cal_energy_4"].array(library="np")
    e7     = tree["cal_energy_7"].array(library="np")

    # apply electron, vertex, and PCal energy cuts
    mask = (
        (pid == 11) &
        (sector != -9999) &
        (vz >= vz_cut[0]) &
        (vz <= vz_cut[1]) &
        (sector >= 1) & (sector <= 6) &
        (e1 >= 0.07) &         # enforce PCal energy > 0.07 GeV
        (e4 >= 0) & (e7 >= 0) &
        (p > 2.0)
    )
    sector = sector[mask]
    p_vals = p[mask]
    sf     = (e1[mask] + e4[mask] + e7[mask]) / p_vals

    n_total = len(pid)
    n_after = np.count_nonzero(mask)
    print(f"{label}:  total entries = {n_total:,},  after cuts = {n_after:,}")

    # set up 2x3 figure with constrained_layout
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)

    # histogram & fit ranges 2 to 8 GeV
    p_bins = np.linspace(2.0, 8.0, 40)  # momentum bins from 2 to 8
    sf_range = (0.10, 0.40)

    print(f"\nSampling fraction cuts for {label}:")
    for sec in range(1, 7):
        ax = axes.flat[sec-1]
        # select sector
        mask_sec = (sector == sec)
        p_sector = p_vals[mask_sec]
        sf_sector = sf[mask_sec]

        # draw 2D histogram with ROOT-like colormap
        h = ax.hist2d(
            p_sector, sf_sector,
            bins=[p_bins, np.linspace(*sf_range, 80)],
            cmap="jet",
            norm=matplotlib.colors.LogNorm()
        )

        # compute means and sigmas
        centers, means, sigmas = compute_means_sigmas(p_sector, sf_sector, p_bins)

        # fit means and sigmas to degree-2 polynomials
        valid = ~np.isnan(means)
        coef_mean = np.polyfit(centers[valid], means[valid], 2)
        coef_sigma = np.polyfit(centers[valid], sigmas[valid], 2)
        poly_mean = np.poly1d(coef_mean)
        poly_sigma = np.poly1d(coef_sigma)

        # compute mean ± 3 sigma coefficients
        c2_m, c1_m, c0_m = coef_mean
        c2_s, c1_s, c0_s = coef_sigma
        a_minus = c0_m - 3*c0_s
        b_minus = c1_m - 3*c1_s
        c_minus = c2_m - 3*c2_s
        a_plus  = c0_m + 3*c0_s
        b_plus  = c1_m + 3*c1_s
        c_plus  = c2_m + 3*c2_s

        print(
            f"Sector {sec}: sf > ({a_minus:.6f} + {b_minus:.6f}*p + {c_minus:.6f}*p*p) "
            f"&& sf < ({a_plus:.6f} + {b_plus:.6f}*p + {c_plus:.6f}*p*p);"
        )

        # evaluate fits over 2–8 GeV
        p_fit = np.linspace(2.0, 8.0, 200)
        mean_fit = poly_mean(p_fit)
        plus3 = mean_fit + 3 * poly_sigma(p_fit)
        minus3 = mean_fit - 3 * poly_sigma(p_fit)

        # overlay fits: mean solid red, ±3σ dashed red
        ax.plot(p_fit, mean_fit,   color='red', linestyle='-',  linewidth=2, label='mean(p)')
        ax.plot(p_fit, plus3,      color='red', linestyle='--', linewidth=2, label='mean+3σ')
        ax.plot(p_fit, minus3,     color='red', linestyle='--', linewidth=2, label='mean-3σ')

        # aesthetics
        ax.set_xlim(2.0, 8.0)
        ax.set_ylim(*sf_range)
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("p (GeV)")
        ax.set_ylabel("Sampling Fraction")
        ax.legend(loc='upper right', fontsize='small')

    # colorbar
    cbar = fig.colorbar(h[3], ax=axes.ravel().tolist(), shrink=0.9)
    cbar.set_label("Counts (log scale)")

    # add suptitle
    fig.suptitle(f"Electron Sampling Fraction - {label}", fontsize=16)

    # save
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"electron_sampling_fraction_{label}.pdf")
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
    # vertex cuts from 2% threshold:
    vz_cuts = {
        "Su22": (-7.576, 0.303),
        "Fa22": (-5.758, 1.515),
        "Sp23": (-5.758, 1.515)
    }
    outdir = "output/rgc_studies"

    for fname, label in zip(files, labels):
        make_sampling_fraction_plot(fname, label, vz_cuts[label], outdir)

if __name__ == "__main__":
    main()