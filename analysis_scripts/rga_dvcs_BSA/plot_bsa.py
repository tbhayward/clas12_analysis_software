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
PERIOD_LABELS = {
    'DVCS_Fa18_inb': 'Fa18 Inb',
    'DVCS_Fa18_out': 'Fa18 Out',
    'DVCS_Sp19_inb': 'Sp19 Inb'
}

def load_bsa_data(file_path):
    try:
        with open(file_path) as f:
            return {tuple(map(int, k.strip("()").split(','))): v for k, v in json.load(f).items()}
    except Exception as e:
        print(f"Error loading {file_path}: {str(e)}")
        return {}

def collect_bin_data(data_dict, key_base):
    x, y, yerr = [], [], []
    for i in range(N_PHI_BINS):
        key = key_base + (i,)
        if key in data_dict and 'bsa' in data_dict[key] and 'bsa_err' in data_dict[key]:
            if data_dict[key]['bsa'] != 0:
                x.append(phi_deg[i])
                y.append(data_dict[key]['bsa'])
                yerr.append(data_dict[key]['bsa_err'])
    return x, y, yerr

def bsa_fit_function(phi, a1, b1):
    b1 = np.clip(b1, -0.99, 0.99)
    return a1 * np.sin(phi) / (1 + b1 * np.cos(phi))

def plot_raw_bsa(binning_csv, bsa_dir="bsa_results", output_dir="bsa_plots/raw"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    for dvcs_period in PERIOD_LABELS.keys():
        eppi0_period = dvcs_period.replace("DVCS", "eppi0")
        dvcs_data = load_bsa_data(f"{bsa_dir}/raw_bsa_dvcs_{dvcs_period}.json")
        eppi0_data = load_bsa_data(f"{bsa_dir}/raw_bsa_eppi0_{eppi0_period}.json")

        for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
            xB_avg = 0.5*(xB_min + xB_max)
            subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
            unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
            unique_t = sorted({(b.tmin, b.tmax) for b in subset})
            
            fig, axs = plt.subplots(len(unique_Q2), len(unique_t), 
                          figsize=(3.5*len(unique_t), 3.5*len(unique_Q2)), 
                          squeeze=False)

            for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
                for c, (t_min, t_max) in enumerate(unique_t):
                    ax = axs[r,c]
                    key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                    
                    # Plot data
                    dvcs_x, dvcs_y, dvcs_yerr = collect_bin_data(dvcs_data, key_base)
                    eppi0_x, eppi0_y, eppi0_yerr = collect_bin_data(eppi0_data, key_base)
                    
                    if dvcs_x:
                        ax.errorbar(dvcs_x, dvcs_y, dvcs_yerr, fmt='o', color='black', 
                                  markersize=5, capsize=3, label='DVCS')
                    if eppi0_x:
                        ax.errorbar(eppi0_x, eppi0_y, eppi0_yerr, fmt='s', color='red',
                                  markersize=4, capsize=3, label='epπ⁰')
                    
                    # Configure axes
                    ax.set(xlim=(0, 360), ylim=(-1, 1), 
                          xticks=[0, 90, 180, 270, 360],
                          title=f"$x_B$={xB_avg:.3f}, $Q^2$={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                    ax.grid(True, alpha=0.3)
                    if dvcs_x or eppi0_x:
                        ax.legend(loc='upper right', frameon=False)

            plt.savefig(f"{output_dir}/{dvcs_period}_xB{i_xB}.png", dpi=150)
            plt.close()

def plot_adjusted_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/adjusted"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
        xB_avg = 0.5*(xB_min + xB_max)
        subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
        unique_t = sorted({(b.tmin, b.tmax) for b in subset})
        
        fig, axs = plt.subplots(len(unique_Q2), len(unique_t),
                              figsize=(3.5*len(unique_t), 3.5*len(unique_Q2)),
                              squeeze=False)

        for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
            for c, (t_min, t_max) in enumerate(unique_t):
                ax = axs[r,c]
                key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                has_data = False
                
                for period in PERIOD_LABELS.keys():
                    data = load_bsa_data(f"{final_dir}/adjusted_bsa_{period}.json")
                    x, y, yerr = collect_bin_data(data, key_base)
                    if x:
                        has_data = True
                        ax.errorbar(x, y, yerr, fmt='o', markersize=5, capsize=3,
                                  label=PERIOD_LABELS[period])
                
                ax.set(xlim=(0, 360), ylim=(-1, 1),
                      title=f"$x_B$={xB_avg:.3f}, $Q^2$={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                ax.grid(True, alpha=0.3)
                if has_data:
                    ax.legend(loc='upper right', frameon=False)

        plt.savefig(f"{output_dir}/adjusted_xB{i_xB}.png", dpi=150)
        plt.close()

def plot_combined_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/combined"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    combined_data = load_bsa_data(f"{final_dir}/combined_bsa.json")
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
        xB_avg = 0.5*(xB_min + xB_max)
        subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
        unique_t = sorted({(b.tmin, b.tmax) for b in subset})
        
        fig, axs = plt.subplots(len(unique_Q2), len(unique_t),
                              figsize=(3.5*len(unique_t), 3.5*len(unique_Q2)),
                              squeeze=False)

        for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
            for c, (t_min, t_max) in enumerate(unique_t):
                ax = axs[r,c]
                key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                
                x, y, yerr = collect_bin_data(combined_data, key_base)
                if not x: continue
                
                try:
                    popt, pcov = curve_fit(bsa_fit_function, np.radians(x), y,
                                         sigma=yerr, bounds=([-np.inf, -1], [np.inf, 1]))
                    fit_x = np.linspace(0, 360, 100)
                    fit_y = bsa_fit_function(np.radians(fit_x), *popt)
                    
                    # Calculate chi2/ndf
                    residuals = y - bsa_fit_function(np.radians(x), *popt)
                    chi2 = np.sum((residuals/yerr)**2)
                    ndf = len(x) - 2
                    
                    ax.plot(fit_x, fit_y, 'r-', lw=1.5)
                    text = (f"$a_1$ = {popt[0]:.2f} ± {np.sqrt(pcov[0,0]):.2f}\n"
                            f"$b_1$ = {popt[1]:.2f} ± {np.sqrt(pcov[1,1]):.2f}\n"
                            f"χ²/ndf = {chi2/ndf:.2f}")
                except:
                    text = "Fit failed"
                
                ax.errorbar(x, y, yerr, fmt='ko', markersize=5, capsize=3)
                ax.text(0.95, 0.95, text, transform=ax.transAxes,
                       ha='right', va='top', fontsize=6)
                ax.set(xlim=(0, 360), ylim=(-1, 1),
                      title=f"$x_B$={xB_avg:.3f}, $Q^2$={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                ax.grid(True, alpha=0.3)

        plt.savefig(f"{output_dir}/combined_xB{i_xB}.png", dpi=150)
        plt.close()

def get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max):
    unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})
    unique_t = sorted({(b.tmin, b.tmax) for b in binning})
    try:
        return (unique_Q2.index((Q2_min, Q2_max)), 
                unique_t.index((t_min, t_max)))
    except ValueError:
        return (-1, -1)