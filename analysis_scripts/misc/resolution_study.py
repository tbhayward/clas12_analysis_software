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

# Branches and their updated plot ranges
BRANCH_SETTINGS = [
    ("Mx2",   (-0.05, 0.05)),
    ("Mx2_1", (-0.3, 0.3)),
    ("Mx2_2", (0.6, 1.2)),
]

def plot_before_smearing(runs, branch, xlim, output_path):
    """
    For a given missing-mass branch, make a 3x2 grid of histograms:
      - rows = each run
      - left column = Data vs MC overlay (normalized)
      - right column = blank (reserved for 'after' plots)
    """
    # ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 12), sharex=True)

    for i, run in enumerate(runs):
        # open trees
        tree_mc = uproot.open(run["mc_file"])["PhysicsEvents"]
        tree_dt = uproot.open(run["data_file"])["PhysicsEvents"]

        # extract arrays
        mc_vals   = tree_mc[branch].array(library="np")
        data_vals = tree_dt[branch].array(library="np")

        # clip data for stats within xlim
        dv = data_vals[(data_vals >= xlim[0]) & (data_vals <= xlim[1])]
        mu, sigma = np.mean(dv), np.std(dv)

        ax = axes[i, 0]
        ax.hist(data_vals, bins=100, range=xlim, density=True,
                histtype="step",
                label=f"Data (μ={mu:.3f}, σ={sigma:.3f})")
        ax.hist(mc_vals,   bins=100, range=xlim, density=True,
                histtype="step",
                label="MC")
        ax.set_xlim(xlim)
        ax.set_title(run["title"])
        ax.legend()
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
    for branch, xlim in BRANCH_SETTINGS:
        out_pdf = f"output/resolution_study/{branch}_before_smearing.pdf"
        plot_before_smearing(RUNS, branch, xlim, out_pdf)
    #endfor


if __name__ == "__main__":
    main()
#endif