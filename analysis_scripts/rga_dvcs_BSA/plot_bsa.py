#!/usr/bin/env python3
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from load_binning_scheme import load_binning_scheme

# Constants and styling
N_PHI_BINS = 12
phi_edges = np.linspace(0, 2 * np.pi, N_PHI_BINS + 1)
phi_centers = (phi_edges[:-1] + phi_edges[1:]) / 2.0
phi_deg = np.degrees(phi_centers)
PLOT_STYLE = {
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.titlesize': 10
}

def load_bsa_data(file_path):
    """Load BSA data from JSON file with error handling"""
    try:
        with open(file_path) as f:
            data = json.load(f)
        return {tuple(map(int, k.strip("()").split(','))): v for k, v in data.items()}
    except Exception as e:
        print(f"Error loading {file_path}: {str(e)}")
        return {}

def collect_bin_data(data_dict, key_base):
    """Collect non-zero data points with errors for a bin"""
    x, y, yerr = [], [], []
    for i in range(N_PHI_BINS):
        key = key_base + (i,)
        if key in data_dict and 'bsa' in data_dict[key] and 'bsa_err' in data_dict[key]:
            if data_dict[key]['bsa'] != 0:  # Skip zero-valued data points
                x.append(phi_deg[i])
                y.append(data_dict[key]['bsa'])
                yerr.append(data_dict[key]['bsa_err'])
    return x, y, yerr

def plot_raw_bsa(binning_csv, bsa_dir="bsa_results", output_dir="bsa_plots/raw"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    for dvcs_period in ["DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb"]:
        eppi0_period = dvcs_period.replace("DVCS", "eppi0")
        dvcs_data = load_bsa_data(f"{bsa_dir}/raw_bsa_dvcs_{dvcs_period}.json")
        eppi0_data = load_bsa_data(f"{bsa_dir}/raw_bsa_eppi0_{eppi0_period}.json")

        for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
            subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
            unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
            unique_t = sorted({(b.tmin, b.tmax) for b in subset})
            
            nrows, ncols = len(unique_Q2), len(unique_t)
            fig, axs = plt.subplots(nrows, ncols, figsize=(3.5*ncols, 3.5*nrows), squeeze=False)
            
            for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
                for c, (t_min, t_max) in enumerate(unique_t):
                    ax = axs[r,c]
                    key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                    
                    # Plot DVCS data
                    dvcs_x, dvcs_y, dvcs_yerr = collect_bin_data(dvcs_data, key_base)
                    if dvcs_x:
                        ax.errorbar(dvcs_x, dvcs_y, yerr=dvcs_yerr, 
                                  fmt='o', color='black', markersize=5,
                                  capsize=3, label='DVCS')
                    
                    # Plot eppi0 data
                    eppi0_x, eppi0_y, eppi0_yerr = collect_bin_data(eppi0_data, key_base)
                    if eppi0_x:
                        ax.errorbar(eppi0_x, eppi0_y, yerr=eppi0_yerr,
                                  fmt='s', color='red', markersize=4,
                                  capsize=3, label='epπ⁰')
                    
                    # Configure axes
                    ax.set_xlim(0, 360)
                    ax.set_xticks([0, 90, 180, 270, 360])
                    ax.set_ylim(-1, 1)
                    ax.grid(True, alpha=0.3)
                    ax.set_title(f"Q²={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                    if dvcs_x or eppi0_x:
                        ax.legend(loc='upper right', frameon=False)
            
            fig.suptitle(f"{dvcs_period} - xB = {0.5*(xB_min+xB_max):.3f}")
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            plt.savefig(f"{output_dir}/{dvcs_period}_xB{i_xB}.png", dpi=150)
            plt.close()

def plot_adjusted_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/adjusted"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    colors = {
        'DVCS_Fa18_inb': ('red', 'o'),
        'DVCS_Fa18_out': ('blue', '^'),
        'DVCS_Sp19_inb': ('green', 'D')
    }
    
    os.makedirs(output_dir, exist_ok=True)
    
    for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
        subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
        unique_t = sorted({(b.tmin, b.tmax) for b in subset})
        
        nrows, ncols = len(unique_Q2), len(unique_t)
        fig, axs = plt.subplots(nrows, ncols, figsize=(3.5*ncols, 3.5*nrows), squeeze=False)
        
        for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
            for c, (t_min, t_max) in enumerate(unique_t):
                ax = axs[r,c]
                key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                
                # Plot all periods
                for period, (color, marker) in colors.items():
                    data = load_bsa_data(f"{final_dir}/adjusted_bsa_{period}.json")
                    x, y, yerr = collect_bin_data(data, key_base)
                    if x:
                        label = period.split('_')[-1].capitalize()
                        ax.errorbar(x, y, yerr=yerr, fmt=marker, color=color,
                                   markersize=5, capsize=3, label=label)
                
                # Configure axes
                ax.set_xlim(0, 360)
                ax.set_ylim(-1, 1)
                ax.grid(True, alpha=0.3)
                ax.set_title(f"Q²={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                if any(x for x in [data.keys() for data in colors.keys()]):
                    ax.legend(loc='upper right', frameon=False)
        
        fig.suptitle(f"Adjusted BSA - xB = {0.5*(xB_min+xB_max):.3f}")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(f"{output_dir}/adjusted_xB{i_xB}.png", dpi=150)
        plt.close()

def plot_combined_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/combined"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    combined_data = load_bsa_data(f"{final_dir}/combined_bsa.json")
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
        subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
        unique_t = sorted({(b.tmin, b.tmax) for b in subset})
        
        nrows, ncols = len(unique_Q2), len(unique_t)
        fig, axs = plt.subplots(nrows, ncols, figsize=(3.5*ncols, 3.5*nrows), squeeze=False)
        
        for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
            for c, (t_min, t_max) in enumerate(unique_t):
                ax = axs[r,c]
                key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                
                # Collect data
                x, y, yerr = collect_bin_data(combined_data, key_base)
                if not x:
                    continue
                
                # Perform fit with new function
                try:
                    popt, pcov = curve_fit(bsa_fit_function, 
                                         np.radians(x), y,
                                         sigma=yerr,
                                         absolute_sigma=True)
                    fit_x = np.linspace(0, 360, 100)
                    fit_y = bsa_fit_function(np.radians(fit_x), *popt)
                    ax.plot(fit_x, fit_y, 'r-', lw=1.5)
                except Exception as e:
                    popt, pcov = [np.nan]*2, None  # Now only 2 parameters
                
                # Plot data
                ax.errorbar(x, y, yerr=yerr, fmt='ko', 
                           markersize=5, capsize=3)
                
                # Add fit info for 2 parameters
                if pcov is not None:
                    text = (f"Amp = {popt[0]:.2f} ± {np.sqrt(pcov[0,0]):.2f}\n"
                            f"b1 = {popt[1]:.2f} ± {np.sqrt(pcov[1,1]):.2f}")
                    ax.text(0.95, 0.95, text, transform=ax.transAxes,
                           ha='right', va='top', fontsize=6)
                
                ax.set_xlim(0, 360)
                ax.set_ylim(-1, 1)
                ax.grid(True, alpha=0.3)
                ax.set_title(f"Q²={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
        
        fig.suptitle(f"Combined BSA with Fit - xB = {0.5*(xB_min+xB_max):.3f}")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(f"{output_dir}/combined_xB{i_xB}.png", dpi=150)
        plt.close()

def get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max):
    """Get bin indices with error checking"""
    unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})
    unique_t = sorted({(b.tmin, b.tmax) for b in binning})
    try:
        return (unique_Q2.index((Q2_min, Q2_max)), 
                unique_t.index((t_min, t_max)))
    except ValueError:
        return (-1, -1)  # Handle missing bins gracefully

def bsa_fit_function(phi, Amp, b1):
    """Simplified BSA fitting function: Amp*sin(phi)/(1 + b1*cos(phi))"""
    return Amp * np.sin(phi) / (1 + b1 * np.cos(phi))