#!/usr/bin/env python3

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Constants
N_PHI_BINS = 9

# Load JSON data
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

# Integrate t bins
def integrate_t_bins(input_json, output_json):
    combined_data = load_combined_bsa_json(input_json)
    integrated_results = {}

    bin_groups = {}
    for bin_key, values in combined_data.items():
        xB_idx, Q2_idx, _, phi_idx = bin_key
        key_3d = (xB_idx, Q2_idx, phi_idx)
        if key_3d not in bin_groups:
            bin_groups[key_3d] = []
        bin_groups[key_3d].append((values["bsa"], values["bsa_err"]))

    for key, measurements in bin_groups.items():
        bsa_vals, bsa_errs = zip(*measurements)
        weights = [1 / (err ** 2) for err in bsa_errs]
        total_weight = sum(weights)
        if total_weight > 0:
            combined_bsa = sum(w * val for w, val in zip(weights, bsa_vals)) / total_weight
            combined_err = np.sqrt(1 / total_weight)
        else:
            combined_bsa = np.mean(bsa_vals)
            combined_err = np.std(bsa_vals) if len(bsa_vals) > 1 else bsa_errs[0]

        integrated_results[key_3d] = {
            "bsa": round(combined_bsa, 5),
            "bsa_err": round(combined_err, 5),
            "n_points": len(bsa_vals),
            "valid": True
        }

    with open(output_json, 'w') as f:
        json.dump({str(k): v for k, v in integrated_results.items()}, f, indent=2)

    print(f"Integrated BSA (t-integrated) results saved to: {output_json}")


def plot_integrated_bsa(json_filepath, output_dir="bsa_plots/integrated"):
    os.makedirs(output_dir, exist_ok=True)
    data_dict = load_combined_bsa_json(json_filepath)

    unique_xB = sorted({k[0] for k in data_dict})
    unique_Q2 = sorted({k[1] for k in data_dict}, reverse=False)

    fig, axs = plt.subplots(len(unique_Q2), len(unique_xB), figsize=(20, 20), squeeze=False)

    for x_idx, xB in enumerate(unique_xB):
        Q2_vals = sorted({k[1] for k in data_dict if k[0] == xB})
        for q_idx, Q2 in enumerate(Q2_vals):
            ax = axs[len(Q2_vals)-1 - q_idx, x_idx]  # lowest Q2 at bottom

            key_base = (xB, Q2)
            x, y, yerr = [], [], []

            for phi_idx in range(N_PHI_BINS):
                key = key_base + (phi_idx,)
                if key in data_dict:
                    phi_center = (phi_idx + 0.5) * 360.0 / N_PHI_BINS
                    bsa_val = data_dict[key]['bsa']
                    if -0.6 <= bsa_val <= 0.6:
                        x.append(phi_center)
                        y.append(bsa_val)
                        yerr.append(data_dict[key]['bsa_err'])

            if not x:
                continue

            ax.errorbar(x, y, yerr, fmt='ko', markersize=4, capsize=3)

            if len(x) >= 4:
                try:
                    popt, _ = curve_fit(
                        bsa_fit_function,
                        np.radians(x),
                        y,
                        sigma=yerr,
                        p0=[0, 0.2, -0.4],
                        bounds=([-np.inf, -0.6, -0.7], [np.inf, 0.6, 0.7])
                    )
                    fit_x = np.linspace(0, 360, 200)
                    fit_y = bsa_fit_function(np.radians(fit_x), *popt)
                    ax.plot(fit_x, fit_y, 'r-', lw=1.5, label=f'$a_1$={popt[1]:.3f}$\pm${np.sqrt(pcov[1,1]):.3f}\n$b_1$={popt[2]:.3f}$\pm${np.sqrt(pcov[2,2]):.3f}')
                    ax.legend(fontsize='small')
                except RuntimeError:
                    print(f"Fit failed for bin {(xB, Q2)}")

            ax.set_ylim(-1, 1)
            ax.set_xlim(0, 360)
            ax.set_xticks([0, 90, 180, 270, 360])

            if x_idx == 0:
                ax.set_ylabel(r"$A_{LU}$")
            else:
                ax.set_yticklabels([])

            if q_idx == 0:
                ax.set_xlabel(r"$\phi$ (deg)")
            else:
                ax.set_xticklabels([])

            ax.set_title(f"$x_B$={xB}, $Q^2$={Q2}")
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = os.path.join(output_dir, "bsa_integrated_over_t.png")
    plt.savefig(plot_file)
    plt.close()
    print(f"Integrated BSA plots saved to {output_dir}")