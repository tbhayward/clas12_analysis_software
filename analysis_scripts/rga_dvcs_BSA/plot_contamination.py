#!/usr/bin/env python3
"""
plot_contamination.py

This script reads the contamination JSON file (one per run period) and the 
binning scheme (from imports/integrated_bin_v2.csv) and produces a set of plots.
For each unique xB bin (as defined in the overall binning scheme) it produces a canvas 
whose rows are the overall unique Q² bins and columns are the overall unique t bins.
Each subplot shows the contamination (with error bars) as a function of φ (in degrees),
with the y–axis fixed from 0 to 1. The subplot title includes the average xB, Q², and t.
If no entry is found for a given (xB, Q², t, φ) bin in the JSON, c_i and its error are assumed to be 0.
"""

import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt

# Import the binning scheme loader.
from load_binning_scheme import load_binning_scheme

# Global constants: 12 phi bins from 0 to 2π.
N_PHI_BINS = 12
phi_edges = np.linspace(0, 2 * math.pi, N_PHI_BINS + 1)
phi_centers = (phi_edges[:-1] + phi_edges[1:]) / 2.0
phi_centers_deg = np.degrees(phi_centers)

def load_contamination(run_period, contamination_dir="contamination"):
    """
    Loads the contamination JSON file for a given run period.
    The keys in the JSON file are strings of tuples like "(1, 2, 3, 4)".
    This function converts them back to a tuple of ints.
    """
    filename = f"contamination_{run_period}.json"
    path = os.path.join(contamination_dir, filename)
    with open(path, "r") as f:
        data = json.load(f)
    contamination = {}
    for key_str, value in data.items():
        key_tuple = tuple(int(x.strip()) for x in key_str.strip("()").split(","))
        contamination[key_tuple] = value
    return contamination

def plot_contamination_for_run(run_period,
                                 binning_csv="imports/integrated_bin_v2.csv",
                                 contamination_dir="contamination",
                                 output_dir="contamination/plots/"):
    """
    For a given run period (e.g., 'DVCS_Fa18_inb'):
      - Loads the contamination JSON and binning scheme.
      - Uses the overall unique xB, Q², and t boundaries from the binning scheme.
      - For each xB bin, creates a canvas (figure) with rows = overall unique Q² bins and
        columns = overall unique t bins.
      - In each subplot the contamination (with error bars) is plotted versus φ (in degrees).
        If no contamination entry is found for a given bin, it is assumed to be 0.
      - Each subplot’s title shows the average xB, Q², and t for that cell.
      - The y–axis is fixed from 0 to 1.
    """
    # Load contamination data.
    contamination = load_contamination(run_period, contamination_dir)
    
    # Load the binning scheme.
    binning_scheme = load_binning_scheme(binning_csv)
    
    # Overall unique bin boundaries from the scheme.
    overall_unique_xB = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    overall_unique_Q2 = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    overall_unique_t  = sorted(set((b.tmin, b.tmax) for b in binning_scheme))
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop over every overall xB bin.
    for i_xB, (xB_min, xB_max) in enumerate(overall_unique_xB):
        xB_avg = (xB_min + xB_max) / 2.0
        # Use the overall unique Q² and t bins (even if no contamination data exists for some cells).
        nrows = len(overall_unique_Q2)
        ncols = len(overall_unique_t)
        fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols + 2, 3 * nrows + 2), squeeze=False)
        
        for i_q2 in range(nrows):
            Q2_min, Q2_max = overall_unique_Q2[i_q2]
            Q2_avg = (Q2_min + Q2_max) / 2.0
            for i_t in range(ncols):
                t_min, t_max = overall_unique_t[i_t]
                t_avg = (t_min + t_max) / 2.0
                ax = axes[i_q2, i_t]
                phi_vals = []
                c_vals = []
                c_errs = []
                # Loop over all 12 phi bins.
                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, i_q2, i_t, i_phi)
                    phi_vals.append(phi_centers_deg[i_phi])
                    if key in contamination:
                        c_vals.append(contamination[key]['c_i'])
                        c_errs.append(contamination[key]['c_i_err'])
                    else:
                        c_vals.append(0)
                        c_errs.append(0)
                ax.errorbar(phi_vals, c_vals, yerr=c_errs, fmt='o', color='black', capsize=3)
                ax.set_xlim(0, 360)
                ax.set_xticks([0, 90, 180, 270, 360])
                ax.set_ylim(0, 1)  # Fixed y-axis.
                ax.grid(True, linestyle='--', alpha=0.5)
                # Set subplot title with average xB, Q², and t.
                ax.set_title(rf"$x_B={xB_avg:.3f},\ Q^2={Q2_avg:.2f},\ -t={t_avg:.2f}$", fontsize=9)
                if i_q2 == nrows - 1:
                    ax.set_xlabel(r"$\phi\ (\deg)$", fontsize=9)
                if i_t == 0:
                    ax.set_ylabel(r"$\pi^0$ contamination", fontsize=9)
        
        fig.suptitle(rf"Run Period: {run_period} -- x_B={xB_avg:.3f}", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        outpath = os.path.join(output_dir, f"plot_contamination_{run_period}_xB_{i_xB}.png")
        plt.savefig(outpath, dpi=150)
        plt.close(fig)
        print(f"Saved plot for {run_period} xB bin {i_xB} to {outpath}")