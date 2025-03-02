#!/usr/bin/env python3
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from load_binning_scheme import load_binning_scheme

# ... [keep constants and styling section same] ...

def load_bsa_data(file_path):
    """Load BSA data from JSON file with enhanced error handling"""
    print(f"\n[STATUS] Loading BSA data from: {file_path}")
    try:
        with open(file_path) as f:
            data = json.load(f)
        converted = {tuple(map(int, k.strip("()").split(',')): v for k, v in data.items()}
        print(f"[SUCCESS] Loaded {len(converted)} entries from {os.path.basename(file_path)}")
        return converted
    except Exception as e:
        print(f"[ERROR] Failed to load {file_path}: {str(e)}")
        return {}

def plot_raw_bsa(binning_csv, bsa_dir="bsa_results", output_dir="bsa_plots/raw"):
    print("\n" + "="*60)
    print("Starting raw BSA plotting process")
    print(f"Output directory: {output_dir}")
    
    plt.style.use(PLOT_STYLE)
    print("\n[STATUS] Loading binning scheme...")
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    print(f"Found {len(unique_xB)} xB bins in total")

    os.makedirs(output_dir, exist_ok=True)
    print(f"\n[STATUS] Created output directory: {output_dir}")

    for dvcs_period in ["DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb"]:
        print("\n" + "-"*50)
        print(f"Processing period: {dvcs_period}")
        eppi0_period = dvcs_period.replace("DVCS", "eppi0")
        
        # Load data files
        dvcs_file = f"{bsa_dir}/raw_bsa_dvcs_{dvcs_period}.json"
        eppi0_file = f"{bsa_dir}/raw_bsa_eppi0_{eppi0_period}.json"
        print(f"\n[STATUS] Loading DVCS data from: {dvcs_file}")
        dvcs_data = load_bsa_data(dvcs_file)
        print(f"[STATUS] Loading eppi0 data from: {eppi0_file}")
        eppi0_data = load_bsa_data(eppi0_file)

        for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
            xB_avg = 0.5*(xB_min + xB_max)
            print(f"\nProcessing xB bin {i_xB} ({xB_min:.3f}-{xB_max:.3f}, avg {xB_avg:.3f})")
            
            subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
            unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
            unique_t = sorted({(b.tmin, b.tmax) for b in subset})
            print(f"Found {len(unique_Q2)} Q² bins and {len(unique_t)} t bins")

            nrows, ncols = len(unique_Q2), len(unique_t)
            print(f"Creating {nrows}x{ncols} plot grid")
            fig, axs = plt.subplots(nrows, ncols, figsize=(3.5*ncols, 3.5*nrows), squeeze=False)

            for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
                for c, (t_min, t_max) in enumerate(unique_t):
                    ax = axs[r,c]
                    key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                    
                    # Plot DVCS data
                    dvcs_x, dvcs_y, dvcs_yerr = collect_bin_data(dvcs_data, key_base)
                    if dvcs_x:
                        print(f"Plotting DVCS data for Q²={0.5*(Q2_min+Q2_max):.2f}, t={0.5*(t_min+t_max):.2f}")
                        ax.errorbar(dvcs_x, dvcs_y, yerr=dvcs_yerr, 
                                  fmt='o', color='black', markersize=5,
                                  capsize=3, label='DVCS')
                    
                    # Plot eppi0 data
                    eppi0_x, eppi0_y, eppi0_yerr = collect_bin_data(eppi0_data, key_base)
                    if eppi0_x:
                        print(f"Plotting epπ⁰ data for Q²={0.5*(Q2_min+Q2_max):.2f}, t={0.5*(t_min+t_max):.2f}")
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

            plt_path = f"{output_dir}/{dvcs_period}_xB{i_xB}.png"
            fig.suptitle(f"{dvcs_period} - xB = {xB_avg:.3f}")
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            plt.savefig(plt_path, dpi=150)
            plt.close()
            print(f"[SUCCESS] Saved plot: {plt_path}")

def plot_adjusted_bsa(binning_csv, final_dir="final_results", output_dir="bsa_plots/adjusted"):
    print("\n" + "="*60)
    print("Starting adjusted BSA plotting process")
    print(f"Output directory: {output_dir}")
    
    plt.style.use(PLOT_STYLE)
    print("\n[STATUS] Loading binning scheme...")
    binning = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    colors = {
        'DVCS_Fa18_inb': ('red', 'o'),
        'DVCS_Fa18_out': ('blue', '^'),
        'DVCS_Sp19_inb': ('green', 'D')
    }
    
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n[STATUS] Created output directory: {output_dir}")

    for i_xB, (xB_min, xB_max) in enumerate(unique_xB):
        xB_avg = 0.5*(xB_min + xB_max)
        print(f"\nProcessing xB bin {i_xB} ({xB_min:.3f}-{xB_max:.3f}, avg {xB_avg:.3f})")
        
        subset = [b for b in binning if (b.xBmin, b.xBmax) == (xB_min, xB_max)]
        unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in subset})
        unique_t = sorted({(b.tmin, b.tmax) for b in subset})

        nrows, ncols = len(unique_Q2), len(unique_t)
        print(f"Creating {nrows}x{ncols} plot grid")
        fig, axs = plt.subplots(nrows, ncols, figsize=(3.5*ncols, 3.5*nrows), squeeze=False)

        for r, (Q2_min, Q2_max) in enumerate(unique_Q2):
            for c, (t_min, t_max) in enumerate(unique_t):
                ax = axs[r,c]
                key_base = (i_xB, *get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max))
                has_data = False
                
                # Plot all periods
                for period, (color, marker) in colors.items():
                    data_file = f"{final_dir}/adjusted_bsa_{period}.json"
                    print(f"[STATUS] Loading adjusted data from: {data_file}")
                    data = load_bsa_data(data_file)
                    x, y, yerr = collect_bin_data(data, key_base)
                    
                    if x:
                        has_data = True
                        label = period.split('_')[-1].capitalize()
                        print(f"Plotting {label} data for Q²={0.5*(Q2_min+Q2_max):.2f}, t={0.5*(t_min+t_max):.2f}")
                        ax.errorbar(x, y, yerr=yerr, fmt=marker, color=color,
                                   markersize=5, capsize=3, label=label)
                
                # Configure axes
                ax.set_xlim(0, 360)
                ax.set_ylim(-1, 1)
                ax.grid(True, alpha=0.3)
                ax.set_title(f"Q²={0.5*(Q2_min+Q2_max):.2f}, -t={0.5*(t_min+t_max):.2f}")
                if has_data:
                    ax.legend(loc='upper right', frameon=False)

        plt_path = f"{output_dir}/adjusted_xB{i_xB}.png"
        fig.suptitle(f"Adjusted BSA - xB = {xB_avg:.3f}")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(plt_path, dpi=150)
        plt.close()
        print(f"[SUCCESS] Saved plot: {plt_path}")

# ... [keep other functions same but ensure status messages are present] ...

def get_bin_indices(binning, Q2_min, Q2_max, t_min, t_max):
    """Get bin indices with error checking and logging"""
    unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})
    unique_t = sorted({(b.tmin, b.tmax) for b in binning})
    try:
        indices = (unique_Q2.index((Q2_min, Q2_max)), 
                  unique_t.index((t_min, t_max)))
        print(f"Found bin indices Q²: {indices[0]}, t: {indices[1]}")
        return indices
    except ValueError as e:
        print(f"[WARNING] Bin not found: Q²({Q2_min}-{Q2_max}), t({t_min}-{t_max})")
        return (-1, -1)