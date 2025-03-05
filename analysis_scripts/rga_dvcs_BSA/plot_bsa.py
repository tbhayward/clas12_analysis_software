#!/usr/bin/env python3
import os
import json
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from load_binning_scheme import load_binning_scheme

# Global constants and styling
N_PHI_BINS = 9
phi_edges = np.linspace(0, 2 * np.pi, N_PHI_BINS + 1)
# Midpoints (fallback only)
phi_centers = (phi_edges[:-1] + phi_edges[1:]) / 2.0
phi_deg = np.degrees(phi_centers)
# Pre-calculate phi edges in degrees for bin assignment.
phi_edges_deg = np.degrees(phi_edges)

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
            # Convert keys from string to tuple of ints.
            return {tuple(map(int, k.strip("()").split(','))): v 
                    for k, v in json.load(f).items()}
    except Exception as e:
        print(f"Error loading {file_path}: {str(e)}")
        return {}

def load_global_bin_means(json_file="bin_means_global.json"):
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"Global bin means JSON file not found: {json_file}")
    with open(json_file, "r") as f:
        data = json.load(f)
    global_means = {}
    for key_str, val in data.items():
        key_tuple = tuple(int(x.strip()) for x in key_str.strip("()").split(","))
        global_means[key_tuple] = val
    return global_means

def collect_bin_data(data_dict, key_base, global_bin_means):
    """
    For a given bin defined by key_base, loop over the N_PHI_BINS and return lists:
      x  : actual average φ (in degrees) from global_bin_means (fallback to midpoint)
      y  : the measured BSA values
      yerr: the corresponding errors.
    """
    x, y, yerr = [], [], []
    for i in range(N_PHI_BINS):
        key = key_base + (i,)
        if key in data_dict and 'bsa' in data_dict[key] and 'bsa_err' in data_dict[key]:
            if data_dict[key]['bsa'] != 0:
                # Use the actual averaged φ (from global_bin_means) if available.
                if key in global_bin_means:
                    phi_val = math.degrees(global_bin_means[key]['phi_avg'])
                else:
                    phi_val = phi_deg[i]
                x.append(phi_val)
                y.append(data_dict[key]['bsa'])
                yerr.append(data_dict[key]['bsa_err'])
    return x, y, yerr

def bsa_fit_function(phi, c0, a1, b1):
    # Constrain parameters as needed.
    b1 = np.clip(b1, -0.99999, 0.999999)
    a1 = np.clip(a1, -0.7, 0.7)
    return c0 + (a1 * np.sin(phi)) / (1 + b1 * np.cos(phi))

def get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max):
    unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})
    unique_t = sorted({(b.tmin, b.tmax) for b in binning})
    try:
        return (unique_Q2.index((Q2_min, Q2_max)), unique_t.index((t_min, t_max)))
    except ValueError:
        return (-1, -1)

def plot_raw_bsa(binning_csv, bsa_dir="bsa_results", output_dir="bsa_plots/raw",
                 global_means_file="bin_means_global.json"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    global_bin_means = load_global_bin_means(global_means_file)
    overall_unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})
    overall_unique_t = sorted({(b.tmin, b.tmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    for dvcs_period in PERIOD_LABELS.keys():
        eppi0_period = dvcs_period.replace("DVCS", "eppi0")
        dvcs_data = load_bsa_data(f"{bsa_dir}/raw_bsa_dvcs_{dvcs_period}.json")
        eppi0_data = load_bsa_data(f"{bsa_dir}/raw_bsa_eppi0_{eppi0_period}.json")
        
        for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
            # Compute the actual xB average for this slice from global bin means.
            xB_vals = []
            subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
            unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
            unique_t = sorted({(b.tmin, b.tmax) for b in subset})
            overall_q2_dict = {q: overall_unique_Q2.index(q) for q in unique_Q2}
            overall_t_dict = {t: overall_unique_t.index(t) for t in unique_t}
            for q in unique_Q2:
                for t in unique_t:
                    overall_q2_idx = overall_q2_dict[q]
                    overall_t_idx = overall_t_dict[t]
                    for i_phi in range(N_PHI_BINS):
                        key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                        if key in global_bin_means:
                            xB_vals.append(global_bin_means[key]['xB_avg'])
            if xB_vals:
                xB_avg_slice = np.mean(xB_vals)
            else:
                xB_avg_slice = 0.5 * (xB_min + xB_max)
            
            nrows = len(unique_Q2)
            ncols = len(unique_t)
            fig, axs = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.5 * nrows), squeeze=False)
            
            for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
                for c, (t_min, t_max) in enumerate(unique_t):
                    ax = axs[r, c]
                    key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                    
                    # For the subplot title, compute the actual Q² and t averages from global bin means.
                    Q2_vals = []
                    t_vals = []
                    overall_q2_idx = overall_unique_Q2.index((Q2_min, Q2_max))
                    overall_t_idx = overall_unique_t.index((t_min, t_max))
                    for i_phi in range(N_PHI_BINS):
                        key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                        if key in global_bin_means:
                            Q2_vals.append(global_bin_means[key]['Q2_avg'])
                            t_vals.append(global_bin_means[key]['t_avg'])
                    if Q2_vals:
                        Q2_avg_cell = np.mean(Q2_vals)
                        t_avg_cell = np.mean(t_vals)
                    else:
                        Q2_avg_cell = 0.5 * (Q2_min + Q2_max)
                        t_avg_cell = 0.5 * (t_min + t_max)
                    
                    # Retrieve the BSA data using the updated collect_bin_data.
                    dvcs_x, dvcs_y, dvcs_yerr = collect_bin_data(dvcs_data, key_base, global_bin_means)
                    eppi0_x, eppi0_y, eppi0_yerr = collect_bin_data(eppi0_data, key_base, global_bin_means)
                    
                    if dvcs_x:
                        ax.errorbar(dvcs_x, dvcs_y, dvcs_yerr, fmt='o', color='black', markersize=5,
                                    capsize=3, label='DVCS')
                    if eppi0_x:
                        ax.errorbar(eppi0_x, eppi0_y, eppi0_yerr, fmt='s', color='red', markersize=4,
                                    capsize=3, label='epπ⁰')
                    
                    ax.set(xlim=(0, 360), ylim=(-1, 1),
                           title=fr"$x_B$={xB_avg_slice:.3f}, $Q^2$={Q2_avg_cell:.2f}, $-t$={t_avg_cell:.2f}")
                    ax.grid(True, alpha=0.3)
                    if r == nrows - 1:
                        ax.set_xlabel(r"$\phi$ (deg)")
                    if c == 0:
                        ax.set_ylabel(r"$A_{LU}$")
                    if dvcs_x or eppi0_x:
                        ax.legend(loc='upper right', frameon=False)
            
            plt.savefig(f"{output_dir}/{dvcs_period}_xB{i_xB}.png", dpi=150)
            plt.close()

def plot_adjusted_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/adjusted",
                      global_means_file="bin_means_global.json"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    global_bin_means = load_global_bin_means(global_means_file)
    overall_unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})
    overall_unique_t = sorted({(b.tmin, b.tmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    all_p_values = []
    
    for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
        xB_vals = []
        subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
        unique_t = sorted({(b.tmin, b.tmax) for b in subset})
        overall_q2_dict = {q: overall_unique_Q2.index(q) for q in unique_Q2}
        overall_t_dict = {t: overall_unique_t.index(t) for t in unique_t}
        for q in unique_Q2:
            for t in unique_t:
                overall_q2_idx = overall_q2_dict[q]
                overall_t_idx = overall_t_dict[t]
                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                    if key in global_bin_means:
                        xB_vals.append(global_bin_means[key]['xB_avg'])
        if xB_vals:
            xB_avg_slice = np.mean(xB_vals)
        else:
            xB_avg_slice = 0.5 * (xB_min + xB_max)
        
        nrows = len(unique_Q2)
        ncols = len(unique_t)
        xB_p_values = []
        
        fig, axs = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.5 * nrows), squeeze=False)
        
        for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
            for c, (t_min, t_max) in enumerate(unique_t):
                ax = axs[r, c]
                key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                has_data = False
                # Initialize storage for consistency test.
                phi_data = {i: {'y': [], 'yerr': []} for i in range(N_PHI_BINS)}
                
                overall_q2_idx = overall_unique_Q2.index((Q2_min, Q2_max))
                overall_t_idx = overall_unique_t.index((t_min, t_max))
                for period in PERIOD_LABELS.keys():
                    data = load_bsa_data(f"{final_dir}/adjusted_bsa_{period}.json")
                    x, y, yerr = collect_bin_data(data, key_base, global_bin_means)
                    if x:
                        has_data = True
                        ax.errorbar(x, y, yerr, fmt='o', markersize=5, capsize=3,
                                    label=PERIOD_LABELS[period])
                        # Instead of exact matching, assign phi bin indices using digitize.
                        indices = [int(np.digitize(xi, phi_edges_deg) - 1) for xi in x]
                        for idx, yi, yerri in zip(indices, y, yerr):
                            phi_data[idx]['y'].append(yi)
                            phi_data[idx]['yerr'].append(yerri)
                total_chi2 = 0.0
                total_dof = 0
                for idx in range(N_PHI_BINS):
                    y_vals = phi_data[idx]['y']
                    y_errs = phi_data[idx]['yerr']
                    if len(y_vals) < 2:
                        continue
                    weights = 1 / np.array(y_errs)**2
                    weighted_mean = np.sum(np.array(y_vals) * weights) / np.sum(weights)
                    chi2_contribution = np.sum(((np.array(y_vals) - weighted_mean)**2 / np.array(y_errs)**2))
                    dof_contribution = len(y_vals) - 1
                    total_chi2 += chi2_contribution
                    total_dof += dof_contribution
                p_value = np.nan
                if total_dof > 0:
                    from scipy.stats import chi2
                    p_value = chi2.sf(total_chi2, total_dof)
                    ax.text(0.05, 0.95, f"Consistency p={p_value:.3f}",
                            transform=ax.transAxes, ha='left', va='top', fontsize=6,
                            bbox=dict(facecolor='white', alpha=0.8))
                    xB_p_values.append(p_value)
                    all_p_values.append(p_value)
                
                # Compute actual Q² and t averages for this cell from global bin means.
                Q2_vals = []
                t_vals = []
                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                    if key in global_bin_means:
                        Q2_vals.append(global_bin_means[key]['Q2_avg'])
                        t_vals.append(global_bin_means[key]['t_avg'])
                if Q2_vals:
                    Q2_avg_cell = np.mean(Q2_vals)
                    t_avg_cell = np.mean(t_vals)
                else:
                    Q2_avg_cell = 0.5 * (Q2_min + Q2_max)
                    t_avg_cell = 0.5 * (t_min + t_max)
                ax.set(xlim=(0, 360), ylim=(-1, 1),
                       title=fr"$x_B$={xB_avg_slice:.3f}, $Q^2$={Q2_avg_cell:.2f}, $-t$={t_avg_cell:.2f}")
                ax.grid(True, alpha=0.3)
                if r == nrows - 1:
                    ax.set_xlabel(r"$\phi$ (deg)")
                if c == 0:
                    ax.set_ylabel(r"$A_{LU}$")
        plt.savefig(f"{output_dir}/adjusted_xB{i_xB}.png", dpi=150)
        plt.close()
        
        valid_pvals = [p for p in xB_p_values if not np.isnan(p)]
        if valid_pvals:
            avg_p = np.mean(valid_pvals)
            print(f"xB bin {i_xB} ({xB_avg_slice:.3f}) average consistency p-value: {avg_p:.3f} "
                  f"(Based on {len(valid_pvals)} measurements)")
        else:
            print(f"xB bin {i_xB} ({xB_avg_slice:.3f}) - no valid p-values calculated")
    
    valid_all_pvals = [p for p in all_p_values if not np.isnan(p)]
    if valid_all_pvals:
        overall_avg = np.mean(valid_all_pvals)
        print(f"\nOverall average consistency p-value: {overall_avg:.3f} "
              f"(Based on {len(valid_all_pvals)} total measurements)")
    else:
        print("\nNo valid p-values calculated in any bins")

def plot_combined_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/combined",
                      global_means_file="bin_means_global.json"):
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    combined_data = load_bsa_data(f"{final_dir}/combined_bsa.json")
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    global_bin_means = load_global_bin_means(global_means_file)
    overall_unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})
    overall_unique_t = sorted({(b.tmin, b.tmax) for b in binning})
    
    os.makedirs(output_dir, exist_ok=True)
    
    overall_chi2_ndf_list = []  # to accumulate chi2/ndf of all fits
    
    for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
        xB_vals = []
        subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
        unique_t = sorted({(b.tmin, b.tmax) for b in subset})
        overall_q2_dict = {q: overall_unique_Q2.index(q) for q in unique_Q2}
        overall_t_dict = {t: overall_unique_t.index(t) for t in unique_t}
        for q in unique_Q2:
            for t in unique_t:
                overall_q2_idx = overall_q2_dict[q]
                overall_t_idx = overall_t_dict[t]
                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                    if key in global_bin_means:
                        xB_vals.append(global_bin_means[key]['xB_avg'])
        if xB_vals:
            xB_avg_slice = np.mean(xB_vals)
        else:
            xB_avg_slice = 0.5 * (xB_min + xB_max)
        
        nrows = len(unique_Q2)
        ncols = len(unique_t)
        fig, axs = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.5 * nrows), squeeze=False)
        
        chi2_ndf_list = []  # chi2/ndf for this xB bin
        
        for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
            for c, (t_min, t_max) in enumerate(unique_t):
                ax = axs[r, c]
                key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                x, y, yerr = collect_bin_data(combined_data, key_base, global_bin_means)
                if not x:
                    continue
                try:
                    # Use initial guess: c0=0, a1=0.2, b1=-0.4.
                    popt, pcov = curve_fit(bsa_fit_function, np.radians(x), y, sigma=yerr,
                                             p0=[0, 0.2, -0.4],
                                             bounds=([-np.inf, -0.6, -0.7], [np.inf, 0.6, 0.7]))
                    fit_x = np.linspace(0, 360, 100)
                    fit_y = bsa_fit_function(np.radians(fit_x), *popt)
                    residuals = y - bsa_fit_function(np.radians(x), *popt)
                    chi2 = np.sum((residuals / yerr)**2)
                    ndf = len(x) - 3
                    chi2_ndf = chi2 / ndf if ndf > 0 else np.nan
                    chi2_ndf_list.append(chi2_ndf)
                    overall_chi2_ndf_list.append(chi2_ndf)
                    ax.plot(fit_x, fit_y, 'r-', lw=1.5)
                    text = (f"$c_0$ = {popt[0]:.2f} ± {np.sqrt(pcov[0,0]):.2f}\n"
                            f"$a_1$ = {popt[1]:.2f} ± {np.sqrt(pcov[1,1]):.2f}\n"
                            f"$b_1$ = {popt[2]:.2f} ± {np.sqrt(pcov[2,2]):.2f}\n"
                            f"χ²/ndf = {chi2_ndf:.2f}")
                except Exception as e:
                    print(f"Fit failed: {str(e)}")
                    text = "Fit failed"
                
                ax.errorbar(x, y, yerr, fmt='ko', markersize=5, capsize=3)
                ax.text(0.95, 0.95, text, transform=ax.transAxes,
                        ha='right', va='top', fontsize=6)
                
                overall_q2_idx = overall_unique_Q2.index((Q2_min, Q2_max))
                overall_t_idx = overall_unique_t.index((t_min, t_max))
                Q2_vals = []
                t_vals = []
                for i_phi in range(N_PHI_BINS):
                    key = (i_xB, overall_q2_idx, overall_t_idx, i_phi)
                    if key in global_bin_means:
                        Q2_vals.append(global_bin_means[key]['Q2_avg'])
                        t_vals.append(global_bin_means[key]['t_avg'])
                if Q2_vals:
                    Q2_avg_cell = np.mean(Q2_vals)
                    t_avg_cell = np.mean(t_vals)
                else:
                    Q2_avg_cell = 0.5 * (Q2_min + Q2_max)
                    t_avg_cell = 0.5 * (t_min + t_max)
                
                ax.set(xlim=(0, 360), ylim=(-1, 1),
                       title=fr"$x_B$={xB_avg_slice:.3f}, $Q^2$={Q2_avg_cell:.2f}, $-t$={t_avg_cell:.2f}")
                ax.grid(True, alpha=0.3)
                if r == nrows - 1:
                    ax.set_xlabel("$\phi$ (deg)")
                if c == 0:
                    ax.set_ylabel("$A_{LU}$")
        plt.savefig(f"{output_dir}/combined_xB{i_xB}.png", dpi=150)
        plt.close()
        
        # Filter out NaN values before calculating the mean for this xB bin.
        valid_chi2_ndf = [val for val in chi2_ndf_list if not np.isnan(val)]
        if valid_chi2_ndf:
            mean_chi2_ndf = np.mean(valid_chi2_ndf)
            print(f"xB bin {i_xB} (⟨x_B⟩ = {xB_avg_slice:.3f}) mean χ²/ndf: {mean_chi2_ndf:.2f} "
                  f"(from {len(valid_chi2_ndf)} fits)")
        else:
            print(f"xB bin {i_xB} (⟨x_B⟩ = {xB_avg_slice:.3f}) - no valid fits")
    
    # Filter out NaN values for overall average as well.
    valid_overall = [val for val in overall_chi2_ndf_list if not np.isnan(val)]
    if valid_overall:
        overall_mean = np.mean(valid_overall)
        print(f"\nOverall average χ²/ndf across all fits: {overall_mean:.2f} "
              f"(from {len(valid_overall)} total fits)")
    else:
        print("\nNo valid χ²/ndf values calculated in any bins")