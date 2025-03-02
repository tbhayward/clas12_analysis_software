#!/usr/bin/env python3
"""
plot_contamination.py

Updated to plot both helicity states (+1 and -1) with different colors
"""

import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt

from load_binning_scheme import load_binning_scheme

# Global constants
N_PHI_BINS = 12
phi_edges = np.linspace(0, 2 * math.pi, N_PHI_BINS + 1)
phi_centers = (phi_edges[:-1] + phi_edges[1:]) / 2.0
phi_centers_deg = np.degrees(phi_centers)

def load_contamination(run_period, contamination_dir="contamination"):
    """Load contamination data with helicity separation"""
    filename = f"contamination_{run_period}.json"
    path = os.path.join(contamination_dir, filename)
    with open(path, "r") as f:
        data = json.load(f)
    return {tuple(map(int, k.strip("()").split(","))): v for k, v in data.items()}

def plot_contamination_for_run(run_period, binning_csv="imports/integrated_bin_v2.csv",
                               contamination_dir="contamination", output_dir="contamination/plots/",
                               helicity=None):
    """Plot contamination with helicity separation"""
    # Load both helicity datasets
    contamination = load_contamination(run_period, contamination_dir)
    
    # Validate helicity parameter
    if helicity not in [None, "plus", "minus"]:
        raise ValueError("Helicity must be None, 'plus', or 'minus'")
    
    binning_scheme = load_binning_scheme(binning_csv)
    overall_unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning_scheme})
    overall_unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning_scheme})
    overall_unique_t = sorted({(b.tmin, b.tmax) for b in binning_scheme})
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Create colormap: red for plus, blue for minus
    color_map = {
        "plus": "#e41a1c",  # Red
        "minus": "#377eb8",  # Blue
        "both": ["#e41a1c", "#377eb8"]
    }
    
    for i_xB, (xB_min, xB_max) in enumerate(overall_unique_xB):
        subset = [b for b in binning_scheme if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        unique_Q2_bins = sorted({(b.Q2min, b.Q2max) for b in subset})
        unique_t_bins = sorted({(b.tmin, b.tmax) for b in subset})
        
        xB_avg = (xB_min + xB_max) / 2.0
        nrows, ncols = len(unique_Q2_bins), len(unique_t_bins)
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.5*ncols + 2, 3.5*nrows + 2), squeeze=False)
        
        for r, q2 in enumerate(unique_Q2_bins):
            Q2_min, Q2_max = q2
            Q2_avg = (Q2_min + Q2_max) / 2.0
            overall_q2_idx = overall_unique_Q2.index(q2)
            
            for c, t in enumerate(unique_t_bins):
                t_min, t_max = t
                t_avg = (t_min + t_max) / 2.0
                overall_t_idx = overall_unique_t.index(t)
                ax = axes[r, c]
                
                # Prepare data arrays
                phi_vals = phi_centers_deg.copy()
                c_plus, err_plus = np.zeros(N_PHI_BINS), np.zeros(N_PHI_BINS)
                c_minus, err_minus = np.zeros(N_PHI_BINS), np.zeros(N_PHI_BINS)
                
                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                    if key in contamination:
                        c_plus[i_phi] = contamination[key].get('c_i_plus', 0)
                        err_plus[i_phi] = contamination[key].get('c_i_plus_err', 0)
                        c_minus[i_phi] = contamination[key].get('c_i_minus', 0)
                        err_minus[i_phi] = contamination[key].get('c_i_minus_err', 0)
                
                # Plot both helicity states
                if helicity in [None, "plus"]:
                    ax.errorbar(phi_vals, c_plus, yerr=err_plus, 
                               fmt='o', ms=5, capsize=3, capthick=1,
                               color=color_map["plus"], label=r'$h^{+}$' if r==0 and c==0 else "")
                if helicity in [None, "minus"]:
                    ax.errorbar(phi_vals, c_minus, yerr=err_minus, 
                               fmt='s', ms=5, capsize=3, capthick=1,
                               color=color_map["minus"], label=r'$h^{-}$' if r==0 and c==0 else "")
                
                # Formatting
                ax.set_xlim(0, 360)
                ax.set_xticks([0, 90, 180, 270, 360])
                ax.set_ylim(0, 1)
                ax.grid(True, linestyle='--', alpha=0.3)
                ax.set_title(rf"$Q^2={Q2_avg:.2f},\ -t={t_avg:.2f}$", fontsize=10)
                
                if r == nrows - 1:
                    ax.set_xlabel(r"$\phi\ (\mathrm{deg})$", fontsize=10)
                if c == 0:
                    ax.set_ylabel(r"Contamination", fontsize=10)
                
                # Add legend to first subplot
                if r == 0 and c == 0:
                    ax.legend(loc='upper right', frameon=True, 
                             facecolor='white', framealpha=0.9)
        
        # Final figure formatting
        helicity_suffix = f"_{helicity}" if helicity else ""
        fig.suptitle(rf"{run_period} - $x_B = {xB_avg:.3f}$" + 
                    (f" ({helicity} helicity)" if helicity else ""), 
                    y=0.98, fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        out_name = f"contamination_{run_period}_xB{i_xB}{helicity_suffix}.png"
        plt.savefig(os.path.join(output_dir, out_name), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved {out_name}")