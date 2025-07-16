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
    means = np.full_like(centers, np.nan)
    sigmas = np.full_like(centers, np.nan)
    for i in range(len(centers)):
        mask = (bin_indices == i)
        if np.count_nonzero(mask) > 2:
            vals = sf_vals[mask]
            means[i] = vals.mean()
            sigmas[i] = vals.std()
    return centers, means, sigmas

def make_sampling_fraction_plot(filename, label, vz_cut, outdir):
    """
    2×3 sampling-fraction vs p:
      - fiducial,
      - vertex,
      - p>2,
      - nphe>=2,
      - e1>=0.15, e4>=0, e7>=0,
      - diagonal HTCC cut,
      then profile & fit and print Java‐style sf cuts.
    """
    tree = uproot.open(filename)["PhysicsEvents"]
    arr = tree.arrays([
        "particle_pid","particle_vz","cal_sector",
        "p","cc_nphe_15",
        "cal_energy_1","cal_energy_4","cal_energy_7",
        "cal_lv_1","cal_lw_1",
        "traj_edge_18","traj_edge_36","traj_edge_6",
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

    # fiducial + valid sector
    fid = (
        (lv1 > 9) & (lw1 > 9) &
        (te18 > 3) & (te36 > 10) &
        (
            ((theta > 10) & (te6 > 3)) |
            ((theta <= 10) & (te6 > 10))
        )
    )
    good_sector = (sector != -9999)

    # base mask: all negative‐track species + prior cuts
    base = (
        ((pid == 11)   | (pid == -211) |
         (pid == -321) | (pid == -2212)) &
        good_sector & fid &
        (vz >= vz_cut[0]) & (vz <= vz_cut[1]) &
        (p > 2.0) &
        (nphe >= 2) &
        (e1 >= 0.15) & (e4 >= 0) & (e7 >= 0)
    )

    # diagonal HTCC cut on e4/p vs e1/p
    p_base    = p[base]
    frac_pcal = e1[base] / p_base
    frac_ecin = e4[base] / p_base
    diag      = frac_ecin >= (-0.625 * frac_pcal + 0.15)

    idx0 = np.nonzero(base)[0]
    keep = idx0[diag]

    p_vals  = p[keep]
    sf_vals = (e1[keep] + e4[keep] + e7[keep]) / p_vals
    secs    = sector[keep]

    # binning
    p_bins   = np.linspace(2.0, 8.0, 40)
    sf_range = (0.10, 0.40)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)

    for sec in range(1, 7):
        ax    = axes.flat[sec-1]
        sel   = (secs == sec)
        p_sec = p_vals[sel]
        sf_sec= sf_vals[sel]

        # 2D histogram
        counts2d, xedges, yedges, im = ax.hist2d(
            p_sec, sf_sec,
            bins=[p_bins, np.linspace(*sf_range, 80)],
            cmap="jet", norm=LogNorm()
        )

        centers, means, sigmas = compute_means_sigmas(p_sec, sf_sec, p_bins)
        valid = ~np.isnan(means)

        # fit quadratic to mean & sigma
        cm = np.polyfit(centers[valid], means[valid], 2)
        cs = np.polyfit(centers[valid], sigmas[valid],2)

        c2_m, c1_m, c0_m = cm
        c2_s, c1_s, c0_s = cs

        low0 = c0_m - 3*c0_s
        low1 = c1_m - 3*c1_s
        low2 = c2_m - 3*c2_s
        up0  = c0_m + 3*c0_s
        up1  = c1_m + 3*c1_s
        up2  = c2_m + 3*c2_s

        # Java‐compatible output
        print(
            f"Sector {sec}: sf > ({low0:.6f} + {low1:.6f}*p + {low2:.6f}*p*p) "
            f"&& sf < ({up0:.6f} + {up1:.6f}*p + {up2:.6f}*p*p);"
        )

        # overlay fits
        pf  = np.linspace(2.0, 8.0, 200)
        pm  = np.poly1d(cm)
        ps  = np.poly1d(cs)
        mf  = pm(pf)
        sf3 = ps(pf)

        ax.plot(pf, mf,        'r-',  lw=2)
        ax.plot(pf, mf + 3*sf3, 'r--', lw=2)
        ax.plot(pf, mf - 3*sf3, 'r--', lw=2)

        ax.set_xlim(2.0, 8.0)
        ax.set_ylim(*sf_range)
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("p (GeV)")
        ax.set_ylabel("Sampling Fraction")

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.9)\
       .set_label("Counts (log scale)")

    os.makedirs(outdir, exist_ok=True)
    fig.suptitle(f"Final Sampling Fraction – {label}", fontsize=16)
    fig.savefig(os.path.join(outdir, f"sampling_fraction_{label.replace(' ','_')}.pdf"))
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

    # your actual RGA vertex cuts:
    vz_cuts = {
        "Fa18 Inb": (-6.364, 1.515),
        "Fa18 Out": (-7.879, 0.303),
        "Sp19 Inb": (-6.364, 1.515),
    }

    global outdir
    outdir = "output/rga_studies"

    for fname, label in zip(files, labels):
        make_sampling_fraction_plot(fname, label, vz_cuts[label], outdir)

if __name__ == "__main__":
    main()
