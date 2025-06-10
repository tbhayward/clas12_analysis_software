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

# Only Mx2_1 and Mx2_2, with updated ranges and axis labels
BRANCH_SETTINGS = [
    ("Mx2_1", (-0.4, 0.4),    r"$M_{x (p)}^{2}$ (GeV$^{2}$)"),
    ("Mx2_2", (0.5,  1.5),    r"$M_{x (\gamma)}^{2}$ (GeV$^{2}$)"),
]

# Detector topologies to split by
TOPOLOGIES = [
    {"det1": 1, "det2": 0, "label": "FD–FT"},
    {"det1": 1, "det2": 1, "label": "FD–FD"},
    {"det1": 2, "det2": 0, "label": "CD–FT"},
    {"det1": 2, "det2": 1, "label": "CD–FD"},
]

def plot_by_topology(runs, topologies, branch, xlim, xlabel, output_path):
    """
    For a given branch, make a 3x4 grid of histograms:
      - rows = each run
      - columns = each detector topology
      - Data vs MC overlay (normalized, with μ & σ)
      - apply cuts to DATA: |t1|<1, theta_gamma_gamma<0.4, pTmiss<0.05
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    n_runs = len(runs)
    n_topo = len(topologies)
    fig, axes = plt.subplots(n_runs, n_topo,
                             figsize=(4 * n_topo, 3 * n_runs),
                             sharex=True, sharey=True)

    for i, run in enumerate(runs):
        # load both trees once
        tree_mc = uproot.open(run["mc_file"])["PhysicsEvents"]
        tree_dt = uproot.open(run["data_file"])["PhysicsEvents"]

        # extract branch arrays
        mc_vals_full   = tree_mc[branch].array(library="np")
        dt_vals_full   = tree_dt[branch].array(library="np")
        t1             = tree_dt["t1"].array(library="np")
        theta_gg       = tree_dt["theta_gamma_gamma"].array(library="np")
        pt_miss        = tree_dt["pTmiss"].array(library="np")

        # apply global DATA cuts
        data_mask_cuts = (np.abs(t1) < 1) & (theta_gg < 0.4) & (pt_miss < 0.05)

        for j, topo in enumerate(topologies):
            ax = axes[i, j]

            # topology masks
            mc_mask = (
                (tree_mc["detector1"].array(library="np") == topo["det1"]) &
                (tree_mc["detector2"].array(library="np") == topo["det2"])
            )
            dt_mask = data_mask_cuts & (
                (tree_dt["detector1"].array(library="np") == topo["det1"]) &
                (tree_dt["detector2"].array(library="np") == topo["det2"])
            )

            mc_vals = mc_vals_full[mc_mask]
            data_vals = dt_vals_full[dt_mask]

            # restrict to plotting range for stats
            mv = mc_vals[(mc_vals >= xlim[0]) & (mc_vals <= xlim[1])]
            dv = data_vals[(data_vals >= xlim[0]) & (data_vals <= xlim[1])]

            mu_mc, sigma_mc = np.mean(mv), np.std(mv)
            mu_dt, sigma_dt = np.mean(dv), np.std(dv)

            # plot
            ax.hist(data_vals, bins=100, range=xlim, density=True,
                    histtype="step",
                    label=f"Data (μ={mu_dt:.3f}, σ={sigma_dt:.3f})")
            ax.hist(mc_vals,   bins=100, range=xlim, density=True,
                    histtype="step",
                    label=f"MC   (μ={mu_mc:.3f}, σ={sigma_mc:.3f})")

            ax.set_xlim(xlim)
            ax.set_title(f"{run['title']}\n{topo['label']}", fontsize=10)
            ax.legend(loc="upper right", fontsize=8)

            # only label outer axes
            if i == n_runs - 1:
                ax.set_xlabel(xlabel)
            if j == 0:
                ax.set_ylabel("normalized counts")
        #endfor
    #endfor

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def main():
    for branch, xlim, xlabel in BRANCH_SETTINGS:
        out_pdf = f"output/resolution_study/{branch}_by_topology.pdf"
        plot_by_topology(RUNS, TOPOLOGIES, branch, xlim, xlabel, out_pdf)
    #endfor


if __name__ == "__main__":
    main()
#endif