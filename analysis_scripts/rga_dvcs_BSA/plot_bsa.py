#!/usr/bin/env python3
import os
import json
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import chi2 
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
    a1 = np.clip(a1, -0.99999, 0.999999)
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
                        # ax.set_xlabel(r"$\phi$ (deg)")
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
    
    # Define a fixed color mapping for each period.
    # For example, PERIOD_LABELS might be:
    # {"Fa18 Inb": "Fa18 Inb", "Fa18 Out": "Fa18 Out", "Sp19 Inb": "Sp19 Inb"}
    period_list = list(PERIOD_LABELS.keys())
    colors = ['C0', 'C1', 'C2']  # Extend if you have more than three periods.
    period_colors = {period: colors[i % len(colors)] for i, period in enumerate(period_list)}
    
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
                # Initialize storage for the consistency test.
                phi_data = {i: {'y': [], 'yerr': []} for i in range(N_PHI_BINS)}
                
                overall_q2_idx = overall_unique_Q2.index((Q2_min, Q2_max))
                overall_t_idx = overall_unique_t.index((t_min, t_max))
                # Loop through each period in a fixed order.
                for period in period_list:
                    data = load_bsa_data(f"{final_dir}/adjusted_bsa_{period}.json")
                    x, y, yerr = collect_bin_data(data, key_base, global_bin_means)
                    if x:
                        # Plot the data with the fixed color.
                        ax.errorbar(x, y, yerr, fmt='o', markersize=5, capsize=3,
                                    label=PERIOD_LABELS[period], color=period_colors[period])
                        # Assign phi bin indices using np.digitize.
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
                    p_value = chi2.sf(total_chi2, total_dof)
                    ax.text(0.05, 0.95, f"Consistency p={p_value:.3f}",
                            transform=ax.transAxes, ha='left', va='top', fontsize=6,
                            bbox=dict(facecolor='white', alpha=0.8))
                    xB_p_values.append(p_value)
                    all_p_values.append(p_value)
                
                # Compute the Q² and t averages for this cell.
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
                    # ax.set_xlabel(r"$\phi$ (deg)")
                    ax.set_xlabel(r"$\phi$ (deg)")
                if c == 0:
                    ax.set_ylabel(r"$A_{LU}$")
                
                # Create a consistent legend with dummy handles for all periods.
                handles = []
                for period in period_list:
                    handle = mlines.Line2D([], [], color=period_colors[period],
                                           marker='o', linestyle='None', markersize=5,
                                           label=PERIOD_LABELS[period])
                    handles.append(handle)
                ax.legend(handles=handles, loc='best', fontsize=6)
        
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
    
    overall_chi2_ndf_list = []  # Will accumulate the robust χ²/ndf from all fits
    a1_fits = {}  # Will store fitted parameters: key (i_xB, i_Q2, i_t) -> {a1, a1_err, chi2_ndf}
    
    # Loop over xB bins.
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
        
        chi2_ndf_list = []  # To store the robust χ²/ndf for each cell in this xB bin
        
        # Loop over Q² and t bins within the current xB bin.
        for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
            for c, (t_min, t_max) in enumerate(unique_t):
                ax = axs[r, c]
                key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                x, y, yerr = collect_bin_data(combined_data, key_base, global_bin_means)
                if not x:
                    continue

                # --- Pre-filtering: Remove any point with y outside [-0.6, 0.6] ---
                x_arr = np.array(x)
                y_arr = np.array(y)
                yerr_arr = np.array(yerr)
                mask = (y_arr >= -0.6) & (y_arr <= 0.6)
                x_arr = x_arr[mask]
                y_arr = y_arr[mask]
                yerr_arr = yerr_arr[mask]
                if len(x_arr) < 4:
                    continue
                
                # --- Iterative fitting with outlier removal ---
                max_iter = 10
                removed_points = 0
                for iter_idx in range(max_iter):
                    try:
                        popt, pcov = curve_fit(bsa_fit_function, np.radians(x_arr), y_arr, sigma=yerr_arr,
                                                 p0=[0, 0.2, -0.4],
                                                 bounds=([-np.inf, -0.6, -0.7], [np.inf, 0.6, 0.7]))
                    except Exception as e:
                        print(f"Fit failed for bin {key_base}: {str(e)}")
                        popt, pcov = [0, 0, 0], np.zeros((3, 3))
                        break
                    
                    fit_vals = bsa_fit_function(np.radians(x_arr), *popt)
                    residuals = y_arr - fit_vals
                    chi2_val = np.sum((residuals / yerr_arr)**2)
                    ndf = len(x_arr) - 3
                    current_chi2_ndf = chi2_val / ndf if ndf > 0 else np.nan
                    if current_chi2_ndf <= 3 or len(x_arr) < 4:
                        break
                    else:
                        contributions = (residuals / yerr_arr)**2
                        idx_to_remove = np.argmax(contributions)
                        x_arr = np.delete(x_arr, idx_to_remove)
                        y_arr = np.delete(y_arr, idx_to_remove)
                        yerr_arr = np.delete(yerr_arr, idx_to_remove)
                        removed_points += 1
                
                # --- Final fit evaluation ---
                fit_x = np.linspace(0, 360, 100)
                fit_y = bsa_fit_function(np.radians(fit_x), *popt)
                final_fit = bsa_fit_function(np.radians(x_arr), *popt)
                final_residuals = y_arr - final_fit
                final_chi2 = np.sum((final_residuals / yerr_arr)**2)
                final_ndf = len(x_arr) - 3
                final_chi2_ndf = final_chi2 / final_ndf if final_ndf > 0 else np.nan
                chi2_ndf_list.append(final_chi2_ndf)
                overall_chi2_ndf_list.append(final_chi2_ndf)
                
                # Save the fitted a1 value for this cell.
                # Here, overall_t_idx is the index for the t bin.
                key_fit = (i_xB, overall_q2_dict[(Q2_min, Q2_max)], overall_t_dict[(t_min, t_max)])
                a1_fits[key_fit] = {
                    "a1": popt[1],
                    "a1_err": np.sqrt(pcov[1,1]),
                    "chi2_ndf": final_chi2_ndf
                }
                
                ax.plot(fit_x, fit_y, 'r-', lw=1.5)
                text = (f"$c_0$ = {popt[0]:.2f} ± {np.sqrt(pcov[0,0]):.2f}\n"
                        f"$a_1$ = {popt[1]:.2f} ± {np.sqrt(pcov[1,1]):.2f}\n"
                        f"$b_1$ = {popt[2]:.2f} ± {np.sqrt(pcov[2,2]):.2f}\n"
                        f"χ²/ndf = {final_chi2_ndf:.2f}")
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
        
        valid_chi2_ndf = np.array([val for val in chi2_ndf_list if not np.isnan(val)])
        if len(valid_chi2_ndf) >= 5:
            try:
                df_fit, loc, scale = chi2.fit(valid_chi2_ndf, floc=0)
                robust_chi2_ndf = chi2.ppf(0.5, df_fit, loc=loc, scale=scale)
            except Exception as e:
                print(f"Chi2 fit failed for xB bin {i_xB}: {str(e)}")
                robust_chi2_ndf = np.mean(valid_chi2_ndf)
        elif len(valid_chi2_ndf) > 0:
            robust_chi2_ndf = np.mean(valid_chi2_ndf)
        else:
            robust_chi2_ndf = np.nan
        
        print(f"xB bin {i_xB} (⟨x_B⟩ = {xB_avg_slice:.3f}) robust χ²/ndf: {robust_chi2_ndf:.2f} "
              f"(from {len(valid_chi2_ndf)} fits)")
    
    if len(overall_chi2_ndf_list) >= 5:
        try:
            df_fit_all, loc_all, scale_all = chi2.fit(np.array(overall_chi2_ndf_list), floc=0)
            overall_robust = chi2.ppf(0.5, df_fit_all, loc=loc_all, scale=scale_all)
        except Exception as e:
            print(f"Overall chi2 fit failed: {str(e)}")
            overall_robust = np.mean(overall_chi2_ndf_list)
    elif len(overall_chi2_ndf_list) > 0:
        overall_robust = np.mean(overall_chi2_ndf_list)
    else:
        overall_robust = np.nan
    
    print(f"\nOverall robust χ²/ndf across all fits: {overall_robust:.2f} "
          f"(from {len(overall_chi2_ndf_list)} total fits)")
    
    # Save the fitted a1 parameters to a JSON file.
    a1_fits_out = os.path.join(final_dir, "a1_fits.json")
    # Convert tuple keys to strings.
    a1_fits_str = {str(k): v for k, v in a1_fits.items()}
    with open(a1_fits_out, "w") as f:
        import json
        json.dump(a1_fits_str, f, indent=2)
    print(f"Saved fitted a₁ parameters to {a1_fits_out}")


def plot_a1_vs_t_grid(binning_csv, final_dir="final_results", output_dir="bsa_plots/a1_vs_t_grid",
                      global_means_file="bin_means_global.json", a1_json_file="a1_fits.json"):
    """
    Create a grid of subplots arranged by xB (columns) and Q² (rows) such that the lowest
    xB and Q² values are in the bottom-left. Each subplot plots fitted a₁ vs. -t for that (xB, Q²) cell.
    If a cell is empty, its axis is removed.
    Finally, only the bottom-most subplot in each column displays x-axis tick labels and label,
    and only the left-most subplot in each row displays y-axis tick labels and label.
    All valid axes are forced to have y limits from -0.1 to 0.5.
    """
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    global_means = load_global_bin_means(global_means_file)
    
    a1_path = os.path.join(final_dir, a1_json_file)
    with open(a1_path) as f:
        a1_data = json.load(f)
    # Convert keys from strings like "(i_xB, i_Q2, i_t)" to tuple of ints.
    a1_data = {tuple(map(int, k.strip("()").split(","))): v for k, v in a1_data.items()}
    
    xB_indices = sorted(set(key[0] for key in a1_data.keys()))
    Q2_indices = sorted(set(key[1] for key in a1_data.keys()))
    
    # For the grid plot, we want rows (Q²) sorted in increasing order (lowest Q² first)
    Q2_indices_plot = sorted(Q2_indices)
    xB_indices_plot = sorted(xB_indices)
    
    nrows = len(Q2_indices_plot)
    ncols = len(xB_indices_plot)
    fig, axs = plt.subplots(nrows, ncols, figsize=(2.8*ncols, 3*nrows), sharex=True, sharey=True)
    
    # Reverse vertical order so that row index 0 corresponds to the bottom (lowest Q²).
    axs = axs[::-1, :]
    
    for i, i_Q2 in enumerate(Q2_indices_plot):
        for j, i_xB in enumerate(xB_indices_plot):
            ax = axs[i, j]
            # Set the y-axis limits globally.
            ax.set_ylim(-0.1, 0.5)
            cell_points = []
            for (xx, qq, tt), fit in a1_data.items():
                if xx == i_xB and qq == i_Q2:
                    global_key = (i_xB, i_Q2, tt, 0)
                    if global_key in global_means:
                        # Here we assume that the global means file already stores t as -t.
                        t_val = global_means[global_key].get("t_avg", None)
                        if t_val is not None:
                            cell_points.append((t_val, fit["a1"], fit["a1_err"]))
            if cell_points:
                cell_points.sort(key=lambda tup: tup[0])
                t_vals, a1_vals, a1_err_vals = zip(*cell_points)
                ax.errorbar(t_vals, a1_vals, yerr=a1_err_vals, fmt='o', color='black',
                            markersize=4, capsize=3, linestyle="None")
                # Set title using representative xB and Q² from global_means.
                rep_key = (i_xB, i_Q2, 0, 0)
                if rep_key in global_means:
                    xB_avg = global_means[rep_key].get("xB_avg", None)
                    Q2_avg = global_means[rep_key].get("Q2_avg", None)
                else:
                    xB_avg, Q2_avg = None, None
                title_str = ""
                if xB_avg is not None:
                    title_str += f"xB = {xB_avg:.3f}, "
                if Q2_avg is not None:
                    title_str += f"Q² ≈ {Q2_avg:.2f}"
                ax.set_title(title_str, fontsize=7)
            else:
                fig.delaxes(ax)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0, 1)
    
    # Now, set tick labels and axis labels only on the outer boundaries.
    # For each column, find the bottom-most valid axis and ensure it shows x tick labels and label.
    for j in range(ncols):
        col_axes = [ (i, axs[i,j]) for i in range(nrows) if axs[i,j] in fig.get_axes() ]
        if col_axes:
            bottom_i = min(i for i,ax in col_axes)  # Since row 0 is bottom.
            axs[bottom_i, j].set_xlabel("-t", fontsize=8)
            axs[bottom_i, j].tick_params(axis='x', labelbottom=True)
            # Turn off x tick labels for other axes in this column.
            for i, ax in col_axes:
                if i != bottom_i:
                    ax.tick_params(axis='x', labelbottom=False)
    
    # For each row, find the left-most valid axis and ensure it shows y tick labels and label.
    for i in range(nrows):
        row_axes = [ (j, axs[i,j]) for j in range(ncols) if axs[i,j] in fig.get_axes() ]
        if row_axes:
            left_j = min(j for j, ax in row_axes)
            axs[i, left_j].set_ylabel("$A_{LU}$", fontsize=8)
            axs[i, left_j].tick_params(axis='y', labelleft=True)
            # Turn off y tick labels for other axes in this row.
            for j, ax in row_axes:
                if j != left_j:
                    ax.tick_params(axis='y', labelleft=False)
    
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "a1_vs_t_grid.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved a₁ vs -t grid plot to {out_path}")

def plot_a1_vs_t_by_Q2(binning_csv, final_dir="final_results", output_dir="bsa_plots/a1_vs_t_by_Q2",
                       global_means_file="bin_means_global.json", a1_json_file="a1_fits.json"):
    """
    Create a single-row figure where each subplot corresponds to a Q² bin.
    On each subplot, curves for different xB bins are overlaid (with fixed color/marker per xB).
    The x-axis values (representing -t) are taken directly from the global means.
    The x-axis is forced to range from 0 to 1.
    The legend shows the average xB value for each curve.
    The far right subplot shows y-axis ticks on the right.
    """
    plt.style.use(PLOT_STYLE)
    binning = load_binning_scheme(binning_csv)
    global_means = load_global_bin_means(global_means_file)
    
    a1_path = os.path.join(final_dir, a1_json_file)
    with open(a1_path) as f:
        a1_data = json.load(f)
    a1_data = {tuple(map(int, k.strip("()").split(","))): v for k, v in a1_data.items()}
    # Keys are (i_xB, i_Q2, i_t)
    
    Q2_indices = sorted(set(key[1] for key in a1_data.keys()))
    xB_indices = sorted(set(key[0] for key in a1_data.keys()))
    
    bin_means_map = {}
    for i_xB in xB_indices:
        for i_Q2 in Q2_indices:
            key = (i_xB, i_Q2, 0, 0)
            if key in global_means:
                bin_means_map[(i_xB, i_Q2)] = {"xB_avg": global_means[key].get("xB_avg", None),
                                               "Q2_avg": global_means[key].get("Q2_avg", None)}
            else:
                bin_means_map[(i_xB, i_Q2)] = {"xB_avg": None, "Q2_avg": None}
    
    # Group a1 data by Q² then by xB.
    # Structure: data_by_Q2[i_Q2][i_xB] = list of (t_val, a1, a1_err)
    data_by_Q2 = {}
    for (i_xB, i_Q2, i_t), fit in a1_data.items():
        global_key = (i_xB, i_Q2, i_t, 0)
        if global_key in global_means:
            t_val = global_means[global_key].get("t_avg", None)
        else:
            t_val = None
        if t_val is None:
            continue
        # Here we assume t_val is already -t.
        if i_Q2 not in data_by_Q2:
            data_by_Q2[i_Q2] = {}
        if i_xB not in data_by_Q2[i_Q2]:
            data_by_Q2[i_Q2][i_xB] = []
        data_by_Q2[i_Q2][i_xB].append((t_val, fit["a1"], fit["a1_err"]))
    for i_Q2 in data_by_Q2:
        for i_xB in data_by_Q2[i_Q2]:
            data_by_Q2[i_Q2][i_xB].sort(key=lambda tup: tup[0])
    
    ncols = len(Q2_indices)
    fig, axs = plt.subplots(1, ncols, figsize=(3.5*ncols, 4), sharex=True, sharey=True)
    if ncols == 1:
        axs = [axs]
    
    # Assign fixed styles for xB bins.
    color_list = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
    marker_list = ['o', 's', '^', 'D', 'v', 'p']
    xB_style = {}
    for i, i_xB in enumerate(xB_indices):
        xB_style[i_xB] = (color_list[i % len(color_list)], marker_list[i % len(marker_list)])
    
    xB_avg_map = {}
    for i_xB in xB_indices:
        possible_keys = [key for key in global_means.keys() if key[0]==i_xB and key[2]==0]
        if possible_keys:
            xB_avg_map[i_xB] = global_means[possible_keys[0]].get("xB_avg", None)
        else:
            xB_avg_map[i_xB] = None
    
    for idx, i_Q2 in enumerate(sorted(Q2_indices)):
        ax = axs[idx]
        rep_keys = [key for key in global_means.keys() if key[1]==i_Q2 and key[2]==0]
        Q2_avg = global_means[rep_keys[0]].get("Q2_avg", None) if rep_keys else None
        title_str = f"Q² ≈ {Q2_avg:.2f}" if Q2_avg is not None else f"Q² index {i_Q2}"
        ax.set_title(title_str, fontsize=9)
        ax.set_xlim(0, 1)
        if i_Q2 in data_by_Q2:
            for i_xB in sorted(data_by_Q2[i_Q2].keys()):
                points = data_by_Q2[i_Q2][i_xB]
                points.sort(key=lambda tup: tup[0])
                t_vals, a1_vals, a1_err_vals = zip(*points)
                color, marker = xB_style[i_xB]
                ax.errorbar(t_vals, a1_vals, yerr=a1_err_vals, fmt=marker,
                            color=color, markersize=5, capsize=3, linestyle="None")
            handles = []
            for i_xB in sorted(data_by_Q2[i_Q2].keys()):
                color, marker = xB_style[i_xB]
                label = f"xB ≈ {xB_avg_map[i_xB]:.3f}" if xB_avg_map[i_xB] is not None else f"xB index {i_xB}"
                handle = mlines.Line2D([], [], color=color, marker=marker, linestyle="None", markersize=5, label=label)
                handles.append(handle)
            ax.legend(handles=handles, fontsize=7, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("-t")
        if idx == ncols - 1:
            ax.yaxis.tick_right()
    axs[0].set_ylabel("$A_{LU}$")
    
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "a1_vs_t_by_Q2.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved a₁ vs -t plot by Q² to {out_path}")