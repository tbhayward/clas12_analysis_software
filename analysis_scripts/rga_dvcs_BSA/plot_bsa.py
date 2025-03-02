#!/usr/bin/env python3
import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Import binning scheme loading from existing code
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
    """Load BSA data from JSON file and convert keys to tuples"""
    with open(file_path) as f:
        data = json.load(f)
    return {tuple(map(int, k.strip("()").split(','))): v for k, v in data.items()}

def bsa_fit_function(phi, Amp, a1, a2, b1):
    """BSA fitting function: Amp*(a1*sin(phi) + a2*sin(2phi)) / (1 + b1*cos(phi))"""
    return Amp * (a1*np.sin(phi) + a2*np.sin(2*phi)) / (1 + b1*np.cos(phi))

def plot_raw_bsa(binning_csv, bsa_dir="bsa_results", output_dir="bsa_plots/raw"):
    """Plot raw DVCS and eppi0 BSAs on same plots"""
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    for period in ["DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb"]:
        dvcs_data = load_bsa_data(f"{bsa_dir}/raw_bsa_dvcs_{period}.json")
        eppi0_data = load_bsa_data(f"{bsa_dir}/raw_bsa_eppi0_{period}.json")
        
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
                    
                    # Plot DVCS
                    dvcs_vals = [dvcs_data.get(key_base + (i,), {}).get('bsa', 0) 
                                for i in range(N_PHI_BINS)]
                    ax.errorbar(phi_deg, dvcs_vals, fmt='o', color='black', label='DVCS')
                    
                    # Plot eppi0
                    eppi0_vals = [eppi0_data.get(key_base + (i,), {}).get('bsa', 0)
                                 for i in range(N_PHI_BINS)]
                    ax.errorbar(phi_deg, eppi0_vals, fmt='s', color='red', label='epπ⁰')
                    
                    ax.set_xlim(0, 360)
                    ax.set_xticks([0, 90, 180, 270, 360])
                    ax.set_ylim(-1, 1)
                    ax.grid(True, alpha=0.3)
                    ax.set_title(f"Q²={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                    if r == nrows-1: ax.set_xlabel("φ (deg)")
                    if c == 0: ax.set_ylabel("BSA")
            
            fig.suptitle(f"{period} - xB = {0.5*(xB_min+xB_max):.3f}")
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            plt.savefig(f"{output_dir}/{period}_xB{i_xB}.png", dpi=150)
            plt.close()

def plot_adjusted_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/adjusted"):
    """Plot adjusted BSAs for all three periods"""
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    colors = {'DVCS_Fa18_inb': 'red', 'DVCS_Fa18_out': 'blue', 'DVCS_Sp19_inb': 'green'}
    
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
                
                for period in colors.keys():
                    data = load_bsa_data(f"{final_dir}/adjusted_bsa_{period}.json")
                    vals = [data.get(key_base + (i,), {}).get('bsa', 0)
                            for i in range(N_PHI_BINS)]
                    ax.plot(phi_deg, vals, 'o-', color=colors[period], 
                           markersize=3, label=period.split('_')[-1])
                
                ax.set_xlim(0, 360)
                ax.set_ylim(-1, 1)
                ax.grid(True, alpha=0.3)
                ax.set_title(f"Q²={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                if r == nrows-1: ax.set_xlabel("φ (deg)")
                if c == 0: ax.set_ylabel("BSA")
        
        fig.suptitle(f"Adjusted BSA - xB = {0.5*(xB_min+xB_max):.3f}")
        fig.legend(['Fa18 Inb', 'Fa18 Out', 'Sp19 Inb'], loc='upper right')
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(f"{output_dir}/adjusted_xB{i_xB}.png", dpi=150)
        plt.close()

def plot_combined_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/combined"):
    """Plot combined BSA with theoretical fit"""
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
                
                # Get data points
                yvals = [combined_data.get(key_base + (i,), {}).get('bsa', 0)
                         for i in range(N_PHI_BINS)]
                yerrs = [combined_data.get(key_base + (i,), {}).get('bsa_err', 0)
                         for i in range(N_PHI_BINS)]
                
                # Perform fit
                try:
                    popt, pcov = curve_fit(bsa_fit_function, phi_edges[:-1], yvals,
                                          sigma=yerrs, absolute_sigma=True)
                    fit_phi = np.linspace(0, 2*np.pi, 100)
                    fit_deg = np.degrees(fit_phi)
                    fit_bsa = bsa_fit_function(fit_phi, *popt)
                except:
                    popt, pcov = [np.nan]*4, None
                    fit_deg, fit_bsa = [], []
                
                # Plotting
                ax.errorbar(phi_deg, yvals, yerr=yerrs, fmt='ko', capsize=3)
                if len(fit_deg) > 0:
                    ax.plot(fit_deg, fit_bsa, 'r-', linewidth=1.5)
                
                # Add fit parameters
                textstr = '\n'.join([
                    f"Amp = {popt[0]:.2f} ± {np.sqrt(pcov[0,0]):.2f}",
                    f"a1 = {popt[1]:.2f} ± {np.sqrt(pcov[1,1]):.2f}",
                    f"a2 = {popt[2]:.2f} ± {np.sqrt(pcov[2,2]):.2f}",
                    f"b1 = {popt[3]:.2f} ± {np.sqrt(pcov[3,3]):.2f}"
                ]) if pcov is not None else "Fit failed"
                ax.text(0.95, 0.95, textstr, transform=ax.transAxes,
                       ha='right', va='top', fontsize=6)
                
                ax.set_xlim(0, 360)
                ax.set_ylim(-1, 1)
                ax.grid(True, alpha=0.3)
                ax.set_title(f"Q²={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                if r == nrows-1: ax.set_xlabel("φ (deg)")
                if c == 0: ax.set_ylabel("BSA")
        
        fig.suptitle(f"Combined BSA with Fit - xB = {0.5*(xB_min+xB_max):.3f}")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(f"{output_dir}/combined_xB{i_xB}.png", dpi=150)
        plt.close()

def get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max):
    """Get overall Q² and t bin indices from bin boundaries"""
    unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})
    unique_t = sorted({(b.tmin, b.tmax) for b in binning})
    return (unique_Q2.index((Q2_min, Q2_max)), 
            unique_t.index((t_min, t_max)))