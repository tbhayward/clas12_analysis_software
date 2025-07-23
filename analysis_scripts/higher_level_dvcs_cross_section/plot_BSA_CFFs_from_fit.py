#!/usr/bin/env python3
import argparse
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import ROOT

ROOT.gSystem.Load("DVCS_xsec_C.so")

def parse_fit_file(filepath):
    params = {}
    with open(filepath) as f:
        for line in f:
            if '=' in line:
                key, val = line.strip().split('=')
                params[key.strip()] = float(val.strip())
    return params

def set_cff_params(params, flags):
    """Assigns parameters to ROOT global variables"""
    for cff in ['H', 'Ht', 'E', 'Et']:
        if flags[cff]:
            ROOT.gROOT.ProcessLine(f"r_{cff} = {params[f'r_{cff}']};")
            ROOT.gROOT.ProcessLine(f"alpha0_{cff} = {params[f'alpha0_{cff}']};")
            ROOT.gROOT.ProcessLine(f"alpha1_{cff} = {params[f'alpha1_{cff}']};")
            ROOT.gROOT.ProcessLine(f"n_{cff} = {params[f'n_{cff}']};")
            ROOT.gROOT.ProcessLine(f"b_{cff} = {params[f'b_{cff}']};")
            ROOT.gROOT.ProcessLine(f"Mm2_{cff} = {params[f'Mm2_{cff}']};")
            ROOT.gROOT.ProcessLine(f"P_{cff} = {params[f'P_{cff}']};")
    ROOT.renormImag = params.get('renormImag', 1.0)
    ROOT.renormReal = params.get('renormReal', 0.0)

def compute_bsa(phi_deg, Q2, xB, t, Eb):
    bsa_vals = []
    for phi in phi_deg:
        bsa = ROOT.getBSA(float(phi), Q2, xB, t, Eb)
        bsa_vals.append(bsa)
    return np.array(bsa_vals)

def set_flags(H=True, Ht=True, E=True, Et=True):
    ROOT.hasH = H
    ROOT.hasHt = Ht
    ROOT.hasE = E
    ROOT.hasEt = Et
    return {'H': H, 'Ht': Ht, 'E': E, 'Et': Et}

def plot_bsa(phi_vals, data_vals, model_curves, bin_idx, output_dir):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.errorbar(phi_vals, data_vals[:, 0], yerr=data_vals[:, 1], fmt='o', label="Data", color='black')
    for label, (yvals, style, color) in model_curves.items():
        ax.plot(phi_vals, yvals, linestyle=style, color=color, label=label, linewidth=2)

    ax.set_xlabel(r"$\phi$ (deg)", fontsize=14)
    ax.set_ylabel("BSA", fontsize=14)
    ax.set_title(f"Bin {bin_idx}", fontsize=14)
    ax.grid(True)
    ax.legend()
    ax.xaxis.set_major_locator(MultipleLocator(45))
    ax.set_xlim(0, 360)
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f"{output_dir}/bin_{bin_idx}.png")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot BSA predictions from fit file")
    parser.add_argument("fitfile", help="Path to fit results .txt file")
    parser.add_argument("--CFFs", type=int, choices=[0, 1], default=1,
                        help="CFF plot mode: 0=ImH only, 1=ImH and full fit")
    args = parser.parse_args()

    # Load and parse fit parameters
    fit_params = parse_fit_file(args.fitfile)

    # Extract bin tag from filename
    tag_match = re.search(r'_(\d{8}_\d{6})\.txt$', args.fitfile)
    tag = tag_match.group(1) if tag_match else "unknown"

    print(">> Fitted parameters:", fit_params)

    # Read bin data (from associated bin CSV file)
    csv_file = f"output/binning/bin_values_{tag}.csv"
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"Missing expected bin CSV file: {csv_file}")

    data = np.loadtxt(csv_file, delimiter=',', skiprows=1)
    n_bins = int(data[:, 0].max())
    phi_vals = np.linspace(0, 360, 100)

    print(f">> Found {n_bins} φ-bins")

    for i in range(n_bins):
        bin_data = data[data[:, 0] == i + 1]
        Q2g = bin_data[0][1]
        xBg = bin_data[0][2]
        tg = bin_data[0][3]
        Ebg = bin_data[0][4]
        phi_data = bin_data[:, 5]
        BSA_data = bin_data[:, 6]
        sig_data = bin_data[:, 7]
        data_vals = np.column_stack((BSA_data, sig_data))

        model_curves = {}

        # Original model (all CFFs on, default params)
        defaults = {
            'r_H': 0.9, 'alpha0_H': 0.43, 'alpha1_H': 0.85, 'n_H': 1.35, 'b_H': 0.4, 'Mm2_H': 0.64, 'P_H': 1.0,
            'r_Ht': 0.0, 'alpha0_Ht': 0.0, 'alpha1_Ht': 0.0, 'n_Ht': 1.0, 'b_Ht': 1.0, 'Mm2_Ht': 0.5, 'P_Ht': 1.0,
            'r_E': 0.0, 'alpha0_E': 0.0, 'alpha1_E': 0.0, 'n_E': 1.0, 'b_E': 1.0, 'Mm2_E': 0.5, 'P_E': 1.0,
            'r_Et': 0.0, 'alpha0_Et': 0.0, 'alpha1_Et': 0.0, 'n_Et': 1.0, 'b_Et': 1.0, 'Mm2_Et': 0.5, 'P_Et': 1.0,
            'renormImag': 1.0, 'renormReal': 0.0
        }
        set_flags(True, True, True, True)
        set_cff_params(defaults, {'H': 1, 'Ht': 1, 'E': 1, 'Et': 1})
        bsas_orig = compute_bsa(phi_vals, Q2g, xBg, tg, Ebg)
        model_curves["Original Model"] = (bsas_orig, '-', 'blue')

        # ImH only fit (dashed red)
        set_flags(True, False, False, False)
        set_cff_params(fit_params, {'H': 1, 'Ht': 0, 'E': 0, 'Et': 0})
        bsas_ImH = compute_bsa(phi_vals, Q2g, xBg, tg, Ebg)
        model_curves["ImH Only"] = (bsas_ImH, '--', 'red')

        if args.CFFs == 1:
            # Full fit (dot-dashed green)
            set_flags(True, True, True, True)
            set_cff_params(fit_params, {'H': 1, 'Ht': 1, 'E': 1, 'Et': 1})
            bsas_full = compute_bsa(phi_vals, Q2g, xBg, tg, Ebg)
            model_curves["Fit (All CFFs)"] = (bsas_full, '-.', 'green')

        plot_bsa(phi_data, data_vals, model_curves, i + 1, f"output/BSA_CFFs_{tag}")

if __name__ == "__main__":
    main()