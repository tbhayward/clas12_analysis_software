#!/usr/bin/env python3

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

# List of runs with titles and file paths
RUNS = [
    {
        "name": "rga_fa18_inb",
        "title": "RGA Fa18 Inb",
        "mc_file":  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
                    "dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root",
        "data_file":"/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
                    "dvcs/rga_fa18_inb_epgamma.root"
    },
    {
        "name": "rga_fa18_out",
        "title": "RGA Fa18 Out",
        "mc_file":  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
                    "dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root",
        "data_file":"/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
                    "dvcs/rga_fa18_out_epgamma.root"
    },
    {
        "name": "rga_sp19_inb",
        "title": "RGA Sp19 Inb",
        "mc_file":  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
                    "dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root",
        "data_file":"/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
                    "dvcs/rga_sp19_inb_epgamma.root"
    },
]

# Branches with updated ranges and axis labels
BRANCH_SETTINGS = [
    ("Mx2",   (-0.04, 0.04),   r"$M_{x}^{2}$ (GeV$^{2}$)"),
    ("Mx2_1", (-0.4,  0.4),    r"$M_{x (p)}^{2}$ (GeV$^{2}$)"),
    ("Mx2_2", (0.5,  1.5),     r"$M_{x (\gamma)}^{2}$ (GeV$^{2}$)"),
]

def plot_before_smearing(runs, branch, xlim, xlabel, output_path):
    """
    For a given missing-mass branch, make a 3x2 grid of histograms:
      - rows = each run
      - left column = Data vs MC overlay (normalized, with μ & σ)
      - right column = blank (reserved for 'after' plots)
      - apply cuts to DATA: |t1|<1, theta_gamma_gamma<0.4, pTmiss<0.05
    """
    # ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 12), sharex=True)

    for i, run in enumerate(runs):
        # open trees
        tree_mc = uproot.open(run["mc_file"])["PhysicsEvents"]
        tree_dt = uproot.open(run["data_file"])["PhysicsEvents"]

        # extract arrays
        mc_vals = tree_mc[branch].array(library="np")

        # apply rudimentary cuts to data
        dt_vals_full = tree_dt[branch].array(library="np")
        t1            = tree_dt["t1"].array(library="np")
        theta_gg      = tree_dt["theta_gamma_gamma"].array(library="np")
        pt_miss       = tree_dt["pTmiss"].array(library="np")
        mask_cuts = (np.abs(t1) < 1) & (theta_gg < 0.4) & (pt_miss < 0.05)
        data_vals = dt_vals_full[mask_cuts]

        # clip for statistics within xlim
        dv = data_vals[(data_vals >= xlim[0]) & (data_vals <= xlim[1])]
        mv = mc_vals[(mc_vals       >= xlim[0]) & (mc_vals       <= xlim[1])]
        mu_dt, sigma_dt = np.mean(dv), np.std(dv)
        mu_mc, sigma_mc = np.mean(mv), np.std(mv)

        ax = axes[i, 0]
        ax.hist(data_vals, bins=100, range=xlim, density=True,
                histtype="step",
                label=f"Data (μ={mu_dt:.3f}, σ={sigma_dt:.3f})")
        ax.hist(mc_vals,   bins=100, range=xlim, density=True,
                histtype="step",
                label=f"MC   (μ={mu_mc:.3f}, σ={sigma_mc:.3f})")
        ax.set_xlim(xlim)
        ax.set_title(run["title"])
        ax.set_xlabel(xlabel)
        ax.set_ylabel("normalized counts")
        ax.legend(loc="upper right")
    #endfor

    # hide the right-hand column for now
    for ax in axes[:, 1]:
        ax.axis("off")
    #endfor

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def main():
    # loop over branches, produce one PDF each
    for branch, xlim, xlabel in BRANCH_SETTINGS:
        out_pdf = f"output/resolution_study/{branch}_before_smearing.pdf"
        plot_before_smearing(RUNS, branch, xlim, xlabel, out_pdf)
    #endfor


if __name__ == "__main__":
    main()
#endif