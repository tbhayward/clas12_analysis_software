#!/usr/bin/env python3

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

N_PHI_BINS = 12

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

def integrate_all_bins(input_json, output_json):
    import json
    import numpy as np

    combined_data = load_combined_bsa_json(input_json)
    integrated_results = {}

    phi_groups = {}
    for bin_key, values in combined_data.items():
        phi_idx = bin_key[3]
        bsa_val = values["bsa"]
        bsa_err = values["bsa_err"]

        if -0.5 <= bsa_val <= 0.5 and bsa_err > 0:
            phi_groups.setdefault(phi_idx, []).append((bsa_val, bsa_err))

    for phi_idx, measurements in phi_groups.items():
        bsa_vals, bsa_errs = zip(*measurements)

        weights = [1 / (err ** 2) for err in bsa_errs]
        total_weight = sum(weights)

        if total_weight == 0:
            continue

        combined_bsa = sum(w * val for w, val in zip(weights, bsa_vals)) / total_weight
        combined_err = np.sqrt(1 / total_weight)

        integrated_results[phi_idx] = {
            "bsa": round(combined_bsa, 5),
            "bsa_err": round(combined_err, 5),
            "n_points": len(bsa_vals),
            "valid": True
        }

    with open(output_json, 'w') as f:
        json.dump({str(k): v for k, v in integrated_results.items()}, f, indent=2)

    print(f"Fully integrated BSA results saved to: {output_json}")

def plot_integrated_bsa(json_filepath, output_dir="bsa_plots/integrated"):
    import os
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    os.makedirs(output_dir, exist_ok=True)

    data_dict = load_combined_bsa_json(json_filepath)

    with open("bin_means_global.json", 'r') as f:
        bin_means = json.load(f)

    unique_xB = sorted({k[0] for k in data_dict})
    unique_Q2 = sorted({k[1] for k in data_dict})

    fig, axs = plt.subplots(len(unique_Q2), len(unique_xB), figsize=(15, 15), squeeze=False)

    populated_subplots = np.zeros_like(axs, dtype=bool)

    for x_idx, xB in enumerate(unique_xB):
        for q_idx, Q2 in enumerate(unique_Q2):
            ax = axs[len(unique_Q2)-1-q_idx, x_idx]

            key_base = (xB, Q2)
            x, y, yerr = [], [], []

            for phi_idx in range(N_PHI_BINS):
                key = key_base + (phi_idx,)
                if key in data_dict:
                    phi_center = (phi_idx + 0.5) * 360.0 / N_PHI_BINS
                    bsa_val = data_dict[key]['bsa']
                    if -0.5 <= bsa_val <= 0.5:
                        x.append(phi_center)
                        y.append(bsa_val)
                        yerr.append(data_dict[key]['bsa_err'])

            if not x:
                ax.axis('off')
                continue

            populated_subplots[len(unique_Q2)-1-q_idx, x_idx] = 1

            ax.errorbar(x, y, yerr, fmt='ko', markersize=4, capsize=3)

            fitted = False
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
                    fitted = True
                except Exception as e:
                    print(f"Curve fit failed for bin ({xB}, {Q2}): {e}")

            # ax.set_ylim(-1, 1)
            ax.set_ylim(-0.6, 0.6)
            ax.set_xlim(0, 360)

            # Label x-axis ticks
            if (len(unique_Q2)-1-q_idx, x_idx) == (len(unique_Q2)-1, 0):
                ax.set_xticks([0, 90, 180, 270, 360])
                ax.set_xticklabels(["0", "90", "180", "270", "360"])
            else:
                ax.set_xticks([90, 180, 270, 360])

            # Label y-axis ticks
            if (len(unique_Q2)-1-q_idx, x_idx) == (len(unique_Q2)-1, 0):
                # ax.set_yticks([-1, -0.5, 0, 0.5, 1])
                ax.set_yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
            else:
                ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6])

            # Axis labels
            if not np.any(populated_subplots[len(unique_Q2)-1-q_idx, :x_idx]):
                ax.set_ylabel(r"$A_{LU}$")
            else:
                ax.set_yticklabels([])

            if not np.any(populated_subplots[len(unique_Q2)-q_idx:, x_idx]):
                ax.set_xlabel(r"$\phi$ (deg)")
            else:
                ax.set_xticklabels([])

            # Bin means
            bin_key = f"({xB}, {Q2}, 0, 0)"
            if bin_key in bin_means:
                xB_avg = bin_means[bin_key]["xB_avg"]
                Q2_avg = bin_means[bin_key]["Q2_avg"]
                title_label = f"$x_B$={xB_avg:.2f}, $Q^2$={Q2_avg:.1f}"
                ax.text(0.5, 0.96, title_label, ha='center', va='top', transform=ax.transAxes, fontsize='small')

            # Fit results at bottom center
            if fitted:
                a1, b1 = popt[1], popt[2]
                a1_err, b1_err = np.sqrt(pcov[1, 1]), np.sqrt(pcov[2, 2])
                # fit_label = f"$a_1$={a1:.3f}±{a1_err:.3f}\n$b_1$={b1:.3f}±{b1_err:.3f}"
                fit_label = f"$a_1$={a1:.3f}±{a1_err:.3f}"
                ax.text(0.5, 0.01, fit_label, ha='center', va='bottom', transform=ax.transAxes, fontsize='small')

            ax.grid(True, alpha=0.3)

    # Equation text in bottom-right space
    fig.text(0.90, 0.10 ,
             r"$A_{LU} = c_0 + \frac{a_1 \sin\phi}{1 + b_1 \cos\phi}$",
             ha='right', va='bottom', fontsize=20)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(os.path.join(output_dir, "bsa_integrated_over_t.png"), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "bsa_integrated_over_t.pdf"), bbox_inches='tight')
    plt.close()
    print(f"Integrated BSA plots saved to {output_dir}")

def plot_fully_integrated_bsa(json_filepath, output_dir="bsa_plots/integrated"):
    import os
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    os.makedirs(output_dir, exist_ok=True)

    integrated_data_dict = load_combined_bsa_json(json_filepath)

    N_PHI_BINS = 12
    x, y, yerr = [], [], []

    for phi_idx in range(N_PHI_BINS):
        key = (phi_idx,)
        if key in integrated_data_dict:
            phi_center = (phi_idx + 0.5) * 360.0 / N_PHI_BINS
            bsa_val = integrated_data_dict[key]['bsa']
            if -0.5 <= bsa_val <= 0.5:
                x.append(phi_center)
                y.append(bsa_val)
                yerr.append(integrated_data_dict[key]['bsa_err'])

    if not x:
        print("No valid data points to plot for fully integrated BSA.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    # Increase marker size, line width, and capsize for visibility
    ax.errorbar(x, y, yerr, fmt='ko', markersize=8, capsize=5, capthick=1.5, elinewidth=1.5)

    fitted = False
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
            fit_x = np.linspace(0, 360, 200)
            fit_y = bsa_fit_function(np.radians(fit_x), *popt)

            # Increase the line width explicitly here
            ax.plot(fit_x, fit_y, 'r-', lw=3)

            a1, b1 = popt[1], popt[2]
            a1_err, b1_err = np.sqrt(pcov[1, 1]), np.sqrt(pcov[2, 2])
            # fit_label = f"$a_1$ = {a1:.3f} ± {a1_err:.3f}\n$b_1$ = {b1:.3f} ± {b1_err:.3f}"
            fit_label = f"$a_1$ = {a1:.3f} ± {a1_err:.3f}"
            ax.text(0.5, 0.02, fit_label, ha='center', va='bottom',
                    transform=ax.transAxes, fontsize=18)

            ax.plot(fit_x, fit_y, 'r-', linewidth=3)  # Thicker red fit line
            fitted = True
        except Exception as e:
            print(f"Curve fit failed for fully integrated: {e}")

    ax.set_ylim(-0.4, 0.4)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xticklabels(["0", "90", "180", "270", "360"])
    ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4])

    # Slightly larger labels
    ax.set_ylabel(r"$A_{LU}$", fontsize=18)
    ax.set_xlabel(r"$\phi$ (deg)", fontsize=16)
    ax.set_title(r"$ep\rightarrow e'p'\gamma$", fontsize=18, pad=10)

    # Increase tick label size for readability
    ax.tick_params(axis='both', which='major', labelsize=14, width=1.5, length=7)

    # Make the grid slightly bolder
    ax.grid(True, alpha=0.3, linewidth=1.0)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "bsa_fully_integrated.pdf"))
    plt.savefig(os.path.join(output_dir, "bsa_fully_integrated.png"))
    plt.close()
    print(f"Fully integrated BSA plot saved to {output_dir}")

def plot_fully_integrated_bsa_period_comparison(
    sp18_in_json: str,
    fa18_in_json: str,
    sp18_out_json: str,
    fa18_out_json: str,
    output_dir: str = "bsa_plots/integrated"
):
    """
    1×2 panel comparing fully–integrated BSA per φ-bin:
      • Left : Sp18 vs Fa18 inbending
      • Right: Sp18 vs Fa18 outbending
    """
    import math, os, json
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    N_PHI_BINS = 12

    def load_combined_bsa_json(json_filepath):
        with open(json_filepath) as f:
            data = json.load(f)
        return {
            tuple(map(int, k.strip("()").replace(" ", "").split(","))): v
            for k, v in data.items()
        }
    #enddef

    def bsa_fit_function(phi, c0, a1, b1):
        b1 = np.clip(b1, -0.99999, 0.999999)
        a1 = np.clip(a1, -0.99999, 0.999999)
        return c0 + (a1 * np.sin(phi)) / (1 + b1 * np.cos(phi))
    #enddef

    def load_and_integrate(json_path):
        combined = load_combined_bsa_json(json_path)
        bin_groups = {}

        for (_ix, _iq, _it, iph), vals in combined.items():
            bsa, berr = vals["bsa"], vals["bsa_err"]
            if berr > 0 and abs(bsa) < 0.5:
                bin_groups.setdefault(iph, []).append((bsa, berr))
            #endif
        #endfor

        xs, ys, errs = [], [], []
        for iph in sorted(bin_groups):
            meas = bin_groups[iph]
            weights = [1.0/(e**2) for _, e in meas]
            W = sum(weights)
            avg = sum(w*v for (v,_), w in zip(meas, weights)) / W
            err = math.sqrt(1.0 / W)
            phi_center = (iph + 0.5) * 360.0 / N_PHI_BINS
            xs.append(phi_center)
            ys.append(avg)
            errs.append(err)
        #endfor

        return np.array(xs), np.array(ys), np.array(errs)
    #enddef

    # load & t-integrate each period
    xi_in, yi_in, err_in    = load_and_integrate(sp18_in_json)
    xf_in, yf_in, erf_in    = load_and_integrate(fa18_in_json)
    xi_out, yi_out, err_out = load_and_integrate(sp18_out_json)
    xf_out, yf_out, erf_out = load_and_integrate(fa18_out_json)

    os.makedirs(output_dir, exist_ok=True)
    fig, (ax_in, ax_out) = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

    def plot_and_fit(ax, x, y, yerr, color, label):
        # attempt fit first, to build legend label
        label_text = label
        if len(x) >= 4:
            popt, pcov = curve_fit(
                bsa_fit_function,
                np.radians(x),
                y,
                sigma=yerr,
                p0=[0, 0.2, -0.4],
                bounds=(
                    [-np.inf, 0.15, -np.inf],  # 0.15 ≤ a1
                    [ np.inf, 0.35,  np.inf]   # a1 ≤ 0.35
                )
            )
            a1, a1_err = popt[1], math.sqrt(pcov[1,1])
            label_text = rf"{label}, $a_1={a1:.3f}\pm{a1_err:.3f}$"
            # plot fit curve
            fx = np.linspace(0, 360, 200)
            fy = bsa_fit_function(np.radians(fx), *popt)
            ax.plot(fx, fy,
                    linestyle='--',
                    color=color,
                    linewidth=2)
        #endif

        # now plot markers with the full legend entry
        ax.errorbar(
            x, y, yerr,
            fmt='o', color=color,
            markersize=6, capsize=3,
            label=label_text
        )
    #enddef

    # Inbending panel
    plot_and_fit(ax_in, xi_in, yi_in, err_in,  'red',  'Sp18 Inb')
    plot_and_fit(ax_in, xf_in, yf_in, erf_in,  'blue', 'Fa18 Inb')
    ax_in.set(
        title="Inbending",
        xlabel=r"$\phi$ (deg)",
        ylabel=r"$A_{LU}$",
        xlim=(0,360),
        ylim=(-0.5,0.5)
    )
    ax_in.set_xticks([0,90,180,270,360])
    ax_in.grid(True, alpha=0.3)
    ax_in.legend()

    # Outbending panel
    plot_and_fit(ax_out, xi_out, yi_out, err_out, 'red',  'Sp18 Out')
    plot_and_fit(ax_out, xf_out, yf_out, erf_out, 'blue', 'Fa18 Out')
    ax_out.set(
        title="Outbending",
        xlabel=r"$\phi$ (deg)",
        xlim=(0,360)
    )
    ax_out.set_xticks([0,90,180,270,360])
    ax_out.grid(True, alpha=0.3)
    ax_out.legend()

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "bsa_period_comparison.pdf"), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, "bsa_period_comparison.png"), bbox_inches='tight')
    plt.close()
#enddef