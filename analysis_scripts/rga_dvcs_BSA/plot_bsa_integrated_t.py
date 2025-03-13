#!/usr/bin/env python3

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

N_PHI_BINS = 9

def load_combined_bsa_json(json_filepath):
    with open(json_filepath) as f:
        data = json.load(f)
        combined_data = {
            tuple(map(int, k.strip("()").replace(" ", "").split(","))): v
            for k, v in data.items()
        }
    return combined_data

def bsa_fit_function(phi, c0, a1, b1):
    b1 = np.clip(b1, -0.99999, 0.999999)
    a1 = np.clip(a1, -0.99999, 0.999999)
    return c0 + (a1 * np.sin(phi)) / (1 + b1 * np.cos(phi))

def integrate_t_bins(input_json, output_json):
    combined_data = load_combined_bsa_json(input_json)
    integrated_results = {}

    bin_groups = {}
    for bin_key, values in combined_data.items():
        xB_idx, Q2_idx, _, phi_idx = bin_key
        key_3d = (xB_idx, Q2_idx, phi_idx)
        bin_groups.setdefault(key_3d, []).append((values["bsa"], values["bsa_err"]))

    for key, measurements in bin_groups.items():
        bsa_vals, bsa_errs = zip(*measurements)
        weights = [1 / (err ** 2) for err in bsa_errs]
        total_weight = sum(weights)
        combined_bsa = sum(w * val for w, val in zip(weights, bsa_vals)) / total_weight
        combined_err = np.sqrt(1 / total_weight)

        integrated_results[key] = {
            "bsa": round(combined_bsa, 5),
            "bsa_err": round(combined_err, 5),
            "n_points": len(bsa_vals),
            "valid": True
        }

    with open(output_json, 'w') as f:
        json.dump({str(k): v for k, v in integrated_results.items()}, f, indent=2)
    print(f"Integrated BSA (t-integrated) results saved to: {output_json}")

def plot_integrated_bsa(json_filepath, binning_json="binning.json", output_dir="bsa_plots/integrated"):
    os.makedirs(output_dir, exist_ok=True)
    data_dict = load_combined_bsa_json(json_filepath)

    with open("bin_means/bin_means/bin_means.json") as f:
        bin_means = json.load(f)

    unique_xB = sorted({k[0] for k in data_dict})
    unique_Q2 = sorted({k[1] for k in data_dict})

    fig, axs = plt.subplots(len(unique_Q2), len(unique_xB), figsize=(20, 20), squeeze=False)

    # Track filled subplots
    has_data = np.zeros((len(unique_Q2), len(unique_xB)), dtype=bool)

    for i, xB in enumerate(unique_xB):
        for j, Q2 in enumerate(unique_Q2):
            ax = axs[len(unique_Q2)-1-j, i]  # lowest Q2 at bottom

            x, y, yerr = [], [], []
            for phi_idx in range(N_PHI_BINS):
                key = (xB, Q2, phi_idx)
                if key in data_dict:
                    phi_center = (phi_idx + 0.5) * 360.0 / N_PHI_BINS
                    bsa_val = data_dict[key]['bsa']
                    if -0.6 <= bsa_val <= 0.6:
                        x.append(phi_center)
                        y.append(bsa_val)
                        yerr.append(data_dict[key]['bsa_err'])

            if not x:
                ax.axis('off')
                continue

            ax.errorbar(x, y, yerr, fmt='ko', markersize=4, capsize=3)

            if len(x) >= 4:
                try:
                    popt, pcov = curve_fit(
                        bsa_fit_function,
                        np.radians(x),
                        y,
                        sigma=yerr,
                        p0=[0, 0.2, -0.4],
                        bounds=([-np.inf, -0.6, -0.7], [np.inf, 0.6, 0.7])
                    )
                    fit_x = np.linspace(0, 360, 100)
                    fit_y = bsa_fit_function(np.radians(fit_x), *popt)
                    ax.plot(fit_x, fit_y, 'r-', lw=1.5)

            ax.set_ylim(-1, 1)
            ax.set_xlim(0, 360)
            ax.set_xticks([0, 90, 180, 270, 360])

            # Find furthest left subplot in this row to set ylabel
            if all(not axs[len(unique_Q2)-1-j, idx].lines for idx in range(i)):
                ax.set_ylabel(r"$A_{LU}$")
            else:
                ax.set_yticklabels([])

            # Find bottom-most subplot in this column to set xlabel
            if all(not axs[len(unique_Q2)-1-idx, i].lines for idx in range(j)):
                ax.set_xlabel(r"$\phi$ (deg)")
            else:
                ax.set_xticklabels([])

            # Use bin means for titles
            mean_key = f"({xB}, {Q2}, 0, 0)"
            with open(binning_json, 'r') as f:
                bin_means = json.load(f)
            if mean_key in binning_means:
                xB_mean = binning_means[mean_key]['xB_avg']
                Q2_mean = binning_json[mean_key]['Q2_avg']
                ax.set_title(f"$x_B$={xB_mean:.2f}, $Q^2$={Q2_mean:.2f}")

            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = os.path.join(output_dir, "bsa_integrated_over_t.png")
    plt.savefig(plot_file)
    plt.close()
    print(f"Integrated BSA plots saved to {output_dir}")