import os
import json
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from load_binning_scheme import load_binning_scheme

# Define the BSA fit function.
def bsa_fit_function(phi, c0, a1, b1):
    b1 = np.clip(b1, -0.99999, 0.99999)
    a1 = np.clip(a1, -0.99999, 0.99999)
    return c0 + (a1 * np.sin(phi)) / (1 + b1 * np.cos(phi))

# Load BSA data from JSON file.
def load_bsa_data(filepath):
    with open(filepath) as f:
        return {tuple(map(int, k.strip("()").split(','))): v for k, v in json.load(f).items()}

# Load global bin averages.
def load_global_bin_means(filepath):
    with open(filepath) as f:
        return {tuple(map(int, k.strip("()").split(','))): v for k, v in json.load(f).items()}

# Combine values weighted by uncertainty.
def weighted_mean(values, errors):
    weights = 1 / np.square(errors)
    mean = np.sum(values * weights) / np.sum(weights)
    mean_err = np.sqrt(1 / np.sum(weights))
    return mean, mean_err

# Main plotting function.
def plot_bsa_integrated_t(
    binning_csv, 
    combined_bsa_file="final_results/combined_bsa.json",
    global_means_file="bin_means_global.json",
    output_dir="bsa_plots/integrated_t"
):
    os.makedirs(output_dir, exist_ok=True)
    binning = load_binning_scheme(binning_csv)
    combined_data = load_bsa_data(combined_bsa_file)
    global_bin_means = load_global_bin_means(global_means_file)

    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning})
    unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning})

    N_PHI_BINS = 9

    fig, axs = plt.subplots(len(unique_Q2), len(unique_xB), figsize=(3.5*len(unique_xB), 3.5*len(unique_Q2)), squeeze=False)

    phi_centers = np.degrees(np.linspace(0, 2*np.pi, N_PHI_BINS+1)[:-1] + np.pi/N_PHI_BINS)

    for ix, (xB_min, xB_max) in enumerate(unique_xB):
        for iq, (Q2_min, Q2_max) in enumerate(unique_Q2):

            phi_vals, bsa_vals, bsa_err_vals = [], [], []

            for iphi in range(N_PHI_BINS):
                y_vals, y_errs = [], []

                for key, data in combined_data.items():
                    bin_xB, bin_Q2, bin_t, bin_phi = key
                    binning_entry = binning[bin_xB]

                    if (binning_entry.xBmin, binning_entry.xBmax) == (xB_min, xB_max) and \
                       (binning_entry.Q2min, binning_entry.Q2max) == (Q2_min, Q2_max) and bin_phi == iphi:
                        y_vals.append(data['bsa'])
                        y_errs.append(data['bsa_err'])

                if len(y_vals) > 0:
                    combined_y, combined_err = weighted_mean(np.array(y_vals), np.array(y_errs))
                    phi_vals.append(phi_centers[iphi])
                    bsa_vals.append(combined_y)
                    bsa_err_vals.append(combined_err)

            ax = axs[len(unique_Q2)-1-iq, ix]

            if len(phi_vals) >= 3:
                phi_radians = np.radians(phi_vals)
                try:
                    popt, pcov = curve_fit(bsa_fit_function, phi_radians, bsa_vals, sigma=bsa_err_vals,
                                           p0=[0, 0.2, -0.4], bounds=([-np.inf,-0.6,-0.7],[np.inf,0.6,0.7]))
                    fit_x = np.linspace(0, 360, 100)
                    fit_y = bsa_fit_function(np.radians(fit_x), *popt)
                    ax.plot(fit_x, fit_y, 'r-', lw=1.5)

                    fit_info = f"$c_0$={popt[0]:.2f}\n$a_1$={popt[1]:.2f}\n$b_1$={popt[2]:.2f}"
                    ax.text(0.95, 0.05, fit_info, transform=ax.transAxes, fontsize=7,
                            verticalalignment='bottom', horizontalalignment='right', bbox=dict(alpha=0.7, facecolor='white'))
                except Exception as e:
                    print(f"Fit failed ({xB_min}-{xB_max}, {Q2_min}-{Q2_max}): {e}")

            ax.errorbar(phi_vals, bsa_vals, bsa_err_vals, fmt='ko', markersize=4, capsize=2)

            ax.set_xlim(0, 360)
            ax.set_ylim(-1, 1)
            ax.set_title(f"$x_B$=[{xB_min:.2f},{xB_max:.2f}] $Q^2$=[{Q2_min:.2f},{Q2_max:.2f}]")
            if iq == 0:
                ax.set_xlabel("$\phi$ (deg)")
            if ix == 0:
                ax.set_ylabel("$A_{LU}$")
            ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "integrated_t_bsa.png"), dpi=150)
    plt.close()

# Example call:
# plot_bsa_integrated_t("imports/integrated_bin_v2.csv")
