#!/usr/bin/env python3

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

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
    ("Mx2_1", (-0.5,  0.5), r"$M_{x (p)}^{2}$ (GeV$^{2}$)"),
    ("Mx2_2", (0.4,   1.6), r"$M_{x (\gamma)}^{2}$ (GeV$^{2}$)"),
]

# Detector topologies to split by
TOPOLOGIES = [
    {"det1": 1, "det2": 0, "label": "FD–FT"},
    {"det1": 1, "det2": 1, "label": "FD–FD"},
    {"det1": 2, "det2": 0, "label": "CD–FT"},
    {"det1": 2, "det2": 1, "label": "CD–FD"},
]

def process_run(run):
    """
    Read MC and data for one run, fill histograms & stats for all branches and topologies.
    Returns (run_name, results_dict).
    """
    # Prepare result container
    results = {branch: {} for branch, _, _ in BRANCH_SETTINGS}

    # Open trees once
    tree_mc = uproot.open(run['mc_file'])['PhysicsEvents']
    tree_dt = uproot.open(run['data_file'])['PhysicsEvents']

    # Pre-read data arrays for cuts and topology masks
    det1_dt      = tree_dt['detector1'].array(library='np')
    det2_dt      = tree_dt['detector2'].array(library='np')
    t1           = tree_dt['t1'].array(library='np')
    theta_gg     = tree_dt['theta_gamma_gamma'].array(library='np')
    pt_miss      = tree_dt['pTmiss'].array(library='np')
    mask_cuts_dt = (np.abs(t1) < 1) & (theta_gg < 0.4) & (pt_miss < 0.05)

    for branch, xlim, _ in BRANCH_SETTINGS:
        # Read branch and topology arrays
        mc_vals  = tree_mc[branch].array(library='np')
        dt_vals  = tree_dt[branch].array(library='np')
        det1_mc  = tree_mc['detector1'].array(library='np')
        det2_mc  = tree_mc['detector2'].array(library='np')

        # Define common bins
        bins = np.linspace(xlim[0], xlim[1], 101)

        for topo in TOPOLOGIES:
            # Masks
            mask_mc = (det1_mc == topo['det1']) & (det2_mc == topo['det2'])
            mask_dt = mask_cuts_dt & \
                      (det1_dt == topo['det1']) & (det2_dt == topo['det2'])

            mc_sel = mc_vals[mask_mc]
            dt_sel = dt_vals[mask_dt]

            # Compute stats within x-range
            mv = mc_sel[(mc_sel >= xlim[0]) & (mc_sel <= xlim[1])]
            dv = dt_sel[(dt_sel >= xlim[0]) & (dt_sel <= xlim[1])]
            mu_mc,    sigma_mc    = mv.mean(),     mv.std()
            mu_dt,    sigma_dt    = dv.mean(),     dv.std()

            # Fill normalized histograms
            counts_mc, _ = np.histogram(mc_sel, bins=bins, density=True)
            counts_dt, _ = np.histogram(dt_sel, bins=bins, density=True)

            results[branch][topo['label']] = {
                'bins':       bins,
                'mc_counts':  counts_mc,
                'dt_counts':  counts_dt,
                'mu_mc':      mu_mc,
                'sigma_mc':   sigma_mc,
                'mu_dt':      mu_dt,
                'sigma_dt':   sigma_dt
            }
        #endfor
    #endfor

    return run['name'], results


def plot_results(all_results):
    """
    Generate combined plots for all runs, branches, and topologies.
    """
    for branch, xlim, xlabel in BRANCH_SETTINGS:
        fig, axes = plt.subplots(
            len(RUNS), len(TOPOLOGIES),
            figsize=(4*len(TOPOLOGIES), 3*len(RUNS)),
            sharex=True, sharey=True
        )

        for i, run in enumerate(RUNS):
            run_name = run['name']
            res = all_results[run_name][branch]

            for j, topo in enumerate(TOPOLOGIES):
                ax = axes[i, j]
                r  = res[topo['label']]
                centers = 0.5 * (r['bins'][:-1] + r['bins'][1:])

                ax.step(centers, r['dt_counts'], where='mid',
                        label=f"Data (μ={r['mu_dt']:.3f}, σ={r['sigma_dt']:.3f})")
                ax.step(centers, r['mc_counts'], where='mid',
                        label=f"MC   (μ={r['mu_mc']:.3f}, σ={r['sigma_mc']:.3f})")

                ax.set_xlim(xlim)
                ax.set_title(f"{run['title']}\n{topo['label']}", fontsize=10)
                ax.legend(loc='upper right', fontsize=8)

                if i == len(RUNS) - 1:
                    ax.set_xlabel(xlabel)
                if j == 0:
                    ax.set_ylabel("normalized counts")
            #endfor
        #endfor

        fig.tight_layout()
        out_pdf = f"output/resolution_study/{branch}_by_topology.pdf"
        os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
        fig.savefig(out_pdf)
        plt.close(fig)
    #endfor


def main():
    # Parallel processing per run
    with ProcessPoolExecutor() as executor:
        futures     = [executor.submit(process_run, run) for run in RUNS]
        all_results = {ret[0]: ret[1] for ret in (f.result() for f in futures)}
    #endfor

    # Plot once all data is processed
    plot_results(all_results)


if __name__ == "__main__":
    main()
#endif