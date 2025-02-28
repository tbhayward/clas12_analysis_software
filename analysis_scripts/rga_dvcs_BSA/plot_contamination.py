#!/usr/bin/env python3
"""
plot_contamination.py

This script reads the contamination JSON file (one per run period) and the 
binning scheme (from imports/integrated_bin_v2.csv) and produces a set of plots.
For each unique xB bin that is represented in the contamination data, it creates a canvas 
whose rows are the unique Q² bins and columns are the unique t bins (as they appear for that xB).
Each subplot shows the contamination (with error bars) as a function of φ (in degrees),
with the y–axis fixed from 0 to 1. In each subplot, if a contamination value is missing for a
given φ bin, 0 is assumed. The subplot’s title includes the average Q² and t for that cell,
and the canvas’s suptitle includes the run period and the average xB.
"""

import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt

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
    For a given run period (e.g. 'DVCS_Fa18_inb'):
      - Loads the contamination JSON and binning scheme.
      - Computes the overall unique bin boundaries from the scheme.
      - For each xB bin that appears in the contamination data:
          • Extracts the unique Q² and t bin indices (only those that occur for that xB).
          • Creates a canvas (figure) with rows = these unique Q² bins and columns = these unique t bins.
          • In each subplot, for each φ bin (0–11) a point (with error bar) is plotted.
            If no contamination entry is found for a given (xB, Q², t, φ) bin, c_i and its error are set to 0.
          • Each subplot’s title shows the average Q² and t (from the bin boundaries) for that cell.
          • The canvas’s suptitle shows the run period and the xB average.
      - The y–axis is fixed from 0 to 1.
    """
    # Load contamination data and binning scheme.
    contamination = load_contamination(run_period, contamination_dir)
    binning_scheme = load_binning_scheme(binning_csv)
    
    # Compute overall unique bin boundaries (sorted) from the scheme.
    overall_unique_xB = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    overall_unique_Q2 = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    overall_unique_t  = sorted(set((b.tmin, b.tmax) for b in binning_scheme))
    
    os.makedirs(output_dir, exist_ok=True)
    
    # For each overall xB bin, check if any contamination keys exist.
    for i_xB, (xB_min, xB_max) in enumerate(overall_unique_xB):
        # Get contamination keys for this xB bin.
        keys_for_xB = [key for key in contamination if key[0] == i_xB]
        if not keys_for_xB:
            # Even if no contamination is recorded, you might want a blank canvas.
            # Here we choose to create a canvas using no Q² or t bins.
            print(f"No contamination data for xB bin {i_xB} (xB: {xB_min}–{xB_max}); skipping canvas.")
            continue

        # Determine the unique Q² and t indices that occur for this xB.
        unique_q2_indices = sorted({ key[1] for key in keys_for_xB })
        unique_t_indices  = sorted({ key[2] for key in keys_for_xB })
        
        # Compute the average xB from the bin boundaries.
        xB_avg = (xB_min + xB_max) / 2.0
        
        nrows = len(unique_q2_indices)
        ncols = len(unique_t_indices)
        fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols + 2, 3 * nrows + 2), squeeze=False)
        
        # Loop over the unique Q² and t indices for this xB.
        for r, q2_idx in enumerate(unique_q2_indices):
            Q2_min, Q2_max = overall_unique_Q2[q2_idx]
            Q2_avg = (Q2_min + Q2_max) / 2.0
            for c, t_idx in enumerate(unique_t_indices):
                t_min, t_max = overall_unique_t[t_idx]
                t_avg = (t_min + t_max) / 2.0
                ax = axes[r, c]
                phi_vals = []
                c_vals = []
                c_errs = []
                # Loop over all phi bins (0 to 11).
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
                ax.set_ylim(0, 1)
                ax.grid(True, linestyle='--', alpha=0.5)
                # Set subplot title: include Q² and t averages.
                ax.set_title(rf"$Q^2={Q2_avg:.2f},\ -t={t_avg:.2f}$", fontsize=9)
                if r == nrows - 1:
                    ax.set_xlabel(r"$\phi\ (\deg)$", fontsize=9)
                if c == 0:
                    ax.set_ylabel(r"$\pi^0$ contamination", fontsize=9)
        
        fig.suptitle(rf"Run Period: {run_period} -- x_B={xB_avg:.3f}", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        outpath = os.path.join(output_dir, f"plot_contamination_{run_period}_xB_{i_xB}.png")
        plt.savefig(outpath, dpi=150)
        plt.close(fig)
        print(f"Saved plot for {run_period} xB bin {i_xB} to {outpath}")