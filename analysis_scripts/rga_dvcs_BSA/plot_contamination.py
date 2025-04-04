#!/usr/bin/env python3
"""
plot_contamination.py

This script reads the contamination JSON file (one per run period) and the 
binning scheme (from imports/integrated_bin_v2.csv) and produces a set of plots.
For each overall xB bin (from the binning scheme) a canvas (figure) is produced.
On that canvas, only the Q² and t bins that appear in the xB‐subset of the binning scheme
are used to form the grid (with rows = Q² bins and columns = t bins).
In each subplot the 9 φ points (with error bars) are plotted;
if no contamination value is found for a given (xB,Q²,t,φ) bin, 0 is used.
Each subplot’s title shows the actual average Q² and t (obtained from global averages),
and the canvas suptitle shows the run period and the actual xB average.
The y–axis is fixed from 0 to 1.
"""

import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt

from load_binning_scheme import load_binning_scheme

# Global constants: 12 φ bins from 0 to 2π.
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
    if not os.path.exists(path):
        raise FileNotFoundError(f"[load_contamination] File not found: {path}")
    with open(path, "r") as f:
        data = json.load(f)
    contamination = {}
    for key_str, value in data.items():
        key_tuple = tuple(int(x.strip()) for x in key_str.strip("()").split(","))
        contamination[key_tuple] = value
    return contamination
#endif

def load_global_bin_means(json_file="bin_means_global.json"):
    """
    Loads the global bin-averaged kinematics from a JSON file 
    (keys are strings like '(i_xB, i_Q2, i_t, i_phi)').
    """
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"[load_global_bin_means] File not found: {json_file}")
    with open(json_file, "r") as f:
        data = json.load(f)
    bin_means = {}
    for key_str, val in data.items():
        key_tuple = tuple(int(x.strip()) for x in key_str.strip("()").split(","))
        bin_means[key_tuple] = val  # e.g. {"xB_avg": ..., "Q2_avg": ..., "t_avg": ..., "phi_avg": ...}
    return bin_means
#endif

def plot_contamination_for_run(run_period,
                               binning_csv="imports/integrated_bin_v2.csv",
                               contamination_dir="contamination",
                               output_dir="contamination/plots/",
                               global_bin_means_json="bin_means_global.json"):
    """
    For a given run period (e.g., 'DVCS_Fa18_inb'):
      - Loads the contamination JSON and the binning scheme.
      - Loads the global bin-averaged kinematics from bin_means_global.json.
      - From the binning scheme, computes the overall unique xB bin boundaries.
      - For each xB bin:
          • We collect all relevant (Q², t, φ) bins for that xB. 
          • In each subplot, the 9 φ points (with error bars) are plotted.
          • Each subplot’s title shows the actual average Q² and t 
            (averaged over φ) from the global bin means.
          • The canvas suptitle shows the run period and the actual average xB 
            (averaged over all Q², t, φ in that xB bin).
      - The y–axis is fixed from 0 to 1.
    """
    # Load contamination data, the binning scheme, and global bin means.
    contamination = load_contamination(run_period, contamination_dir)
    binning_scheme = load_binning_scheme(binning_csv)
    bin_means_global = load_global_bin_means(global_bin_means_json)

    # Overall unique xB boundaries (sorted by xBmin).
    overall_unique_xB = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    # Also compute the overall unique Q² and t boundaries from the full scheme.
    overall_unique_Q2 = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    overall_unique_t = sorted(set((b.tmin, b.tmax) for b in binning_scheme))

    os.makedirs(output_dir, exist_ok=True)

    # Loop over each overall xB bin.
    for i_xB, (xB_min, xB_max) in enumerate(overall_unique_xB):
        # Filter the binning scheme to only those bins in this xB range.
        subset = [b for b in binning_scheme if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        
        # The unique Q² and t boundaries in this subset.
        unique_Q2_bins = sorted(set((b.Q2min, b.Q2max) for b in subset))
        unique_t_bins = sorted(set((b.tmin, b.tmax) for b in subset))

        # Mapping from (Q²min, Q²max) -> overall Q² index.
        q2_to_overall = {q2: overall_unique_Q2.index(q2) for q2 in unique_Q2_bins}
        # Mapping from (tmin, tmax) -> overall t index.
        t_to_overall = {t: overall_unique_t.index(t) for t in unique_t_bins}

        # Compute the actual xB average across all (Q², t, φ) bins for this xB slice.
        xB_values_for_slice = []
        for (Q2_min, Q2_max) in unique_Q2_bins:
            overall_q2_idx = q2_to_overall[(Q2_min, Q2_max)]
            for (t_min, t_max) in unique_t_bins:
                overall_t_idx = t_to_overall[(t_min, t_max)]
                for i_phi in range(N_PHI_BINS):
                    key_global = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                    if key_global in bin_means_global:
                        xB_values_for_slice.append(bin_means_global[key_global]["xB_avg"])
        if len(xB_values_for_slice) > 0:
            xB_avg_slice = np.mean(xB_values_for_slice)
        else:
            xB_avg_slice = (xB_min + xB_max) / 2.0

        nrows = len(unique_Q2_bins)
        ncols = len(unique_t_bins)
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(3 * ncols + 2, 3 * nrows + 2),
                                 squeeze=False)

        # Loop over each Q² and t bin for this xB.
        for r, q2 in enumerate(unique_Q2_bins):
            (Q2_min, Q2_max) = q2
            overall_q2_idx = q2_to_overall[q2]
            for c, t in enumerate(unique_t_bins):
                (t_min, t_max) = t
                overall_t_idx = t_to_overall[t]
                ax = axes[r, c]

                # Compute actual Q² and t average for this bin by averaging over φ.
                Q2_vals, t_vals = [], []
                for i_phi in range(N_PHI_BINS):
                    key_global = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                    if key_global in bin_means_global:
                        Q2_vals.append(bin_means_global[key_global]["Q2_avg"])
                        t_vals.append(bin_means_global[key_global]["t_avg"])
                if len(Q2_vals) > 0:
                    Q2_avg_bin = np.mean(Q2_vals)
                    t_avg_bin  = np.mean(t_vals)
                else:
                    Q2_avg_bin = (Q2_min + Q2_max) / 2.0
                    t_avg_bin  = (t_min + t_max) / 2.0

                # Now prepare to plot contamination vs. φ.
                phi_vals = []
                c_vals = []
                c_errs = []

                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                    # Default contamination values.
                    c_val = 0.0
                    c_err = 0.0
                    # Fallback: use midpoint if no global φ average is found.
                    this_phi_deg = phi_centers_deg[i_phi]

                    if key in contamination:
                        c_val = contamination[key]['c_i']
                        c_err = contamination[key]['c_i_err']
                    # If a global average φ exists, use it.
                    if key in bin_means_global:
                        this_phi_rad = bin_means_global[key]['phi_avg']
                        this_phi_deg = math.degrees(this_phi_rad)
                    phi_vals.append(this_phi_deg)
                    c_vals.append(c_val)
                    c_errs.append(c_err)

                ax.errorbar(phi_vals, c_vals, yerr=c_errs,
                            fmt='o', color='black', capsize=3)
                ax.set_xlim(0, 360)
                ax.set_xticks([0, 90, 180, 270, 360])
                ax.set_ylim(0, 1)
                ax.grid(True, linestyle='--', alpha=0.5)
                ax.set_title(rf"$Q^2={Q2_avg_bin:.2f},\ -t={t_avg_bin:.2f}$", fontsize=9)
                if r == nrows - 1:
                    ax.set_xlabel(r"$\phi\ (\deg)$", fontsize=9)
                if c == 0:
                    ax.set_ylabel(r"$\pi^0$ contamination", fontsize=9)
            #endfor (t)
        #endfor (Q²)

        fig.suptitle(rf"Run Period: {run_period} -- $\langle x_B \rangle = {xB_avg_slice:.3f}$",
                     fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        outpath = os.path.join(output_dir, f"plot_contamination_{run_period}_xB_{i_xB}.png")
        plt.savefig(outpath, dpi=150)
        plt.close(fig)
        print(f"[plot_contamination_for_run] Saved plot for {run_period}, xB bin {i_xB} to {outpath}")
    #endfor
#endif

if __name__ == "__main__":
    # Example call:
    run_period = "DVCS_Fa18_inb"  # update as needed
    plot_contamination_for_run(run_period)
#endif