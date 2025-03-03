#!/usr/bin/env python3
"""
plot_contamination.py

This script reads the contamination JSON file (one per run period) and the 
binning scheme (from imports/integrated_bin_v2.csv) and produces a set of plots.
For each overall xB bin (from the binning scheme) a canvas (figure) is produced.
On that canvas, only the Q² and t bins that appear in the xB‐subset of the binning scheme
are used to form the grid (with rows = Q² bins and columns = t bins).
In each subplot the 12 φ points (with error bars) are plotted;
if no contamination value is found for a given (xB,Q²,t,φ) bin, 0 is used.
Each subplot’s title shows the average Q² and t, and the canvas suptitle shows the run period and the xB average.
The y–axis is fixed from 0 to 1.
"""

import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt

from load_binning_scheme import load_binning_scheme

# Global constants: 12 phi bins from 0 to 2π.
N_PHI_BINS = 9
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
      - Loads the contamination JSON and the binning scheme.
      - From the binning scheme, computes the overall unique xB bin boundaries.
      - For each xB bin:
          • Filters the binning scheme for that xB (even if no contamination data exists, a canvas is created).
          • From the xB–subset, obtains the unique Q² and t bin boundaries.
          • Uses the overall sorted Q² and t boundaries (from the full scheme) to determine the indices
            (which are used in the contamination keys). For any missing contamination key, a value of 0 is used.
          • Creates a canvas whose rows correspond to the unique Q² bins (for that xB) and columns to the unique t bins.
          • In each subplot, the 12 φ points (with error bars) are plotted.
          • Each subplot’s title shows the average Q² and t for that cell.
          • The canvas suptitle shows the run period and the xB average.
      - The y–axis in every subplot is fixed from 0 to 1.
    """
    # Load contamination data and the binning scheme.
    contamination = load_contamination(run_period, contamination_dir)
    binning_scheme = load_binning_scheme(binning_csv)
    
    # Overall unique xB boundaries (sorted by xBmin).
    overall_unique_xB = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    
    # Also compute the overall unique Q² and t boundaries from the full scheme.
    overall_unique_Q2 = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    overall_unique_t = sorted(set((b.tmin, b.tmax) for b in binning_scheme))
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop over each overall xB bin.
    for i_xB, (xB_min, xB_max) in enumerate(overall_unique_xB):
        # Filter the binning scheme to only those bins in this xB.
        subset = [b for b in binning_scheme if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        # For this xB, obtain the unique Q² and t boundaries (in the order they appear when sorted).
        unique_Q2_bins = sorted(set((b.Q2min, b.Q2max) for b in subset))
        unique_t_bins = sorted(set((b.tmin, b.tmax) for b in subset))
        
        # For each unique Q² and t in this xB, we need their overall index.
        # Create a mapping from the Q² bin (from the subset) to its index in overall_unique_Q2.
        q2_to_overall = {q2: overall_unique_Q2.index(q2) for q2 in unique_Q2_bins}
        t_to_overall = {t: overall_unique_t.index(t) for t in unique_t_bins}
        
        xB_avg = (xB_min + xB_max) / 2.0
        nrows = len(unique_Q2_bins)
        ncols = len(unique_t_bins)
        fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols + 2, 3 * nrows + 2), squeeze=False)
        
        # Loop over each Q² and t bin for this xB.
        for r, q2 in enumerate(unique_Q2_bins):
            Q2_min, Q2_max = q2
            Q2_avg = (Q2_min + Q2_max) / 2.0
            overall_q2_idx = q2_to_overall[q2]
            for c, t in enumerate(unique_t_bins):
                t_min, t_max = t
                t_avg = (t_min + t_max) / 2.0
                overall_t_idx = t_to_overall[t]
                ax = axes[r, c]
                phi_vals = []
                c_vals = []
                c_errs = []
                # Loop over all 12 φ bins.
                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
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