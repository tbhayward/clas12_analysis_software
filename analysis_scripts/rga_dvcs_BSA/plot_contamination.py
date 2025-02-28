#!/usr/bin/env python3
"""
plot_contamination.py

This script reads the contamination JSON file (one per run period) and the 
binning scheme (from imports/integrated_bin_v2.csv) and produces a set of plots.
For each unique xB bin (as defined in the overall binning scheme) it produces a canvas 
whose rows are the unique Q² bins and columns are the unique t bins that appear in 
the contamination data for that xB. Each subplot shows the contamination (with error 
bars) as a function of φ (in degrees), with the y–axis fixed from 0 to 1 and a title 
showing the average xB, Q² and t for that cell.
If a given (xB, Q², t, φ) bin is missing from the JSON, it is assumed to have c_i = 0.
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
    The JSON keys are strings of tuples like "(1, 2, 3, 4)".
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
    For a given run period (e.g. 'DVCS_Fa18_inb'):
      - Loads the contamination JSON and binning scheme.
      - Uses the overall unique xB, Q² and t boundaries from the binning scheme.
      - For each overall xB bin, it filters the contamination data to find the unique
        Q² and t indices that actually appear.
      - It creates a canvas (figure) with rows = unique Q² bins and columns = unique t bins.
      - In each subplot the contamination (with error bars) is plotted versus φ (in degrees),
        using c_i = 0 and error = 0 for missing bins.
      - The subplot title shows the average xB, Q² and t for that cell.
      - The y–axis is fixed from 0 to 1.
    """
    # Load contamination data.
    contamination = load_contamination(run_period, contamination_dir)
    # Load binning scheme.
    binning_scheme = load_binning_scheme(binning_csv)
    
    # Compute overall unique boundaries.
    overall_unique_xB = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    overall_unique_Q2 = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    overall_unique_t  = sorted(set((b.tmin, b.tmax) for b in binning_scheme))
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop over each overall xB bin.
    for i_xB, (xB_min, xB_max) in enumerate(overall_unique_xB):
        # Filter contamination keys for this xB.
        keys_for_xB = [key for key in contamination.keys() if key[0] == i_xB]
        if not keys_for_xB:
            continue  # Skip if no data for this xB.
        # Determine unique Q² and t indices that appear.
        unique_q2_indices = sorted({ key[1] for key in keys_for_xB })
        unique_t_indices  = sorted({ key[2] for key in keys_for_xB })
        
        xB_avg = (xB_min + xB_max) / 2.0
        nrows = len(unique_q2_indices)
        ncols = len(unique_t_indices)
        fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols + 2, 3*nrows + 2), squeeze=False)
        
        for i, q2_idx in enumerate(unique_q2_indices):
            Q2_min, Q2_max = overall_unique_Q2[q2_idx]
            Q2_avg = (Q2_min + Q2_max) / 2.0
            for j, t_idx in enumerate(unique_t_indices):
                t_min, t_max = overall_unique_t[t_idx]
                t_avg = (t_min + t_max) / 2.0
                ax = axes[i, j]
                # Loop over all 12 phi bins and get contamination (defaulting to 0 if missing)
                phi_vals = []
                c_vals = []
                c_errs = []
                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, q2_idx, t_idx, i_phi)
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
                ax.set_ylim(0, 1)  # Fixed y–axis.
                ax.grid(True, linestyle='--', alpha=0.5)
                ax.set_title(rf"$x_B={xB_avg:.3f},\ Q^2={Q2_avg:.2f},\ -t={t_avg:.2f}$", fontsize=9)
                if i == nrows - 1:
                    ax.set_xlabel(r"$\phi\ (\deg)$", fontsize=9)
                if j == 0:
                    ax.set_ylabel(r"$\pi^0$ contamination", fontsize=9)
        
        fig.suptitle(rf"Run Period: {run_period} -- x_B={xB_avg:.3f}", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        outpath = os.path.join(output_dir, f"plot_contamination_{run_period}_xB_{i_xB}.png")
        plt.savefig(outpath, dpi=150)
        plt.close(fig)
        print(f"Saved plot for {run_period} xB bin {i_xB} to {outpath}")