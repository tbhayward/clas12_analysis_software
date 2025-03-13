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

def plot_integrated_bsa(json_filepath, output_dir="bsa_plots/integrated"):
    import os
    import numpy as np
    import json
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    os.makedirs(output_dir, exist_ok=True)
    data_dict = load_combined_bsa_json(json_filepath)

    # Load global bin means
    with open("bin_means_global.json", 'r') as f:
        bin_means = json.load(f)

    unique_xB = sorted({k[0] for k in data_dict})
    unique_Q2 = sorted({k[1] for k in data_dict})

    fig, axs = plt.subplots(len(unique_Q2), len(unique_xB), figsize=(15, 15), squeeze=False)

    populated_subplots = np.zeros_like(axs, dtype=bool)

    for x_idx, xB in enumerate(unique_xB):
        for q_idx, Q2 in enumerate(unique_Q2):
            ax = axs[len(unique_Q2) - 1 - q_idx, x_idx]

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
                    ax.plot(fit_x, fit_y, 'r-', lw=1.2)
                except RuntimeError:
                    print(f"Fit failed for bin {(xB, Q2)}")

            ax.set_ylim(-1, 1)
            ax.set_xlim(0, 360)
            ax.set_xticks([0, 90, 180, 270, 360])

            bin_mean_key = f"({xB}, {Q2}, 0, 0)"
            if bin_mean_key in bin_means:
                xB_avg = bin_means[bin_mean_key]["xB_avg"]
                Q2_avg = bin_means[bin_mean_key]["Q2_avg"]
                ax.text(0.05, 0.85, f"$x_B$={xB_avg:.2f}, $Q^2$={Q2_avg:.1f}", transform=ax.transAxes, fontsize=9)

            populated_subplots[len(unique_Q2)-1-q_idx, x_idx] = True

    # Label logic for far-left and bottom subplots
    for i in range(len(unique_Q2)):
        row_indices = np.where(populated_subplots[i])[0]
        if row_indices.size > 0:
            leftmost = row_indices[0]
            axs[i, leftmost].set_ylabel(r"$A_{LU}$")
            for idx in row_indices[1:]:
                axs[i, idx].set_yticklabels([])

    for j in range(len(unique_xB)):
        col_indices = np.where(populated_subplots[:, j])[0]
        if col_indices.size > 0:
            bottommost = col_indices[-1]
            axs[bottommost, j].set_xlabel(r"$\phi$ (deg)")
            for idx in col_indices[:-1]:
                axs[idx, j].set_xticklabels([])

    plt.subplots_adjust(hspace=0, wspace=0)
    plot_file = os.path.join(output_dir, "bsa_integrated_over_t.png")
    plt.savefig(plot_file)
    plt.close()
    print(f"Integrated BSA plots saved to {output_dir}")