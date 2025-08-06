#!/usr/bin/env python3
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Proton mass squared (GeV^2)
M_p2 = 0.93827**2

# Gaussian + quadratic background
def gauss_quad(x, A, mu, sigma, a, b, c):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + a*x**2 + b*x + c

# Files and title suffixes for detector1 == 1 (FD) and 2 (CD)
det_configs = {
    1: ("trihadron_fd.pdf", "#pi in FD, p in FD"),
    2: ("trihadron_cd.pdf", "#pi in FD, p in CD")
}

# e_theta bins: [5,14], [14,18], [18,22], [22,26], [26,30] degrees
edges_deg = [5, 14, 18, 22, 26, 30]
bins_deg  = list(zip(edges_deg[:-1], edges_deg[1:]))

# Histogram settings: Mx2 from -0.2 to 0.2 in 40 bins
bins_M = np.linspace(-0.2, 0.2, 41)
centers = 0.5 * (bins_M[:-1] + bins_M[1:])

# Load both ROOT files
f11 = uproot.open("/volatile/clas12/thayward/rgk_dc_study/eppi+pi-/eppi+pi-_cj11.root")
f13 = uproot.open("/volatile/clas12/thayward/rgk_dc_study/eppi+pi-/eppi+pi-_cj13.root")

arrays11 = f11["PhysicsEvents"].arrays(["Mx2", "e_theta", "detector1"], library="np")
arrays13 = f13["PhysicsEvents"].arrays(["Mx2", "e_theta", "detector1"], library="np")

Mx2_11     = arrays11["Mx2"]
e_theta_11 = arrays11["e_theta"] * (180.0 / np.pi)  # radians -> degrees
det1_11    = arrays11["detector1"]

Mx2_13     = arrays13["Mx2"]
e_theta_13 = arrays13["e_theta"] * (180.0 / np.pi)  # radians -> degrees
det1_13    = arrays13["detector1"]

# Loop over detector configurations
for det_val, (outname, title_suffix) in det_configs.items():
    # Prepare lists to store fit results
    angle_means  = []
    mu11_list    = []
    sigma11_list = []
    mu13_list    = []
    sigma13_list = []

    # Create a new figure for this detector
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes_flat = axes.flatten()

    # Loop over e_theta bins
    for i, (low_deg, high_deg) in enumerate(bins_deg):
        ax = axes_flat[i]
        # Select events in this e_theta range AND detector1 == det_val
        mask11 = ((e_theta_11 >= low_deg) & (e_theta_11 < high_deg)
                  & (det1_11 == det_val))
        mask13 = ((e_theta_13 >= low_deg) & (e_theta_13 < high_deg)
                  & (det1_13 == det_val))
        M11 = Mx2_11[mask11]
        M13 = Mx2_13[mask13]
        # Compute mean e_theta for bin
        angle_mean = np.mean(np.concatenate([
            e_theta_11[mask11],
            e_theta_13[mask13]
        ]))
        angle_means.append(angle_mean)
        # Histogram counts and errors
        counts11, _ = np.histogram(M11, bins=bins_M)
        counts13, _ = np.histogram(M13, bins=bins_M)
        errs11 = np.sqrt(counts11)
        errs13 = np.sqrt(counts13)
        # Fit over full range -0.2–0.2
        fit_mask = (centers >= -0.2) & (centers <= 0.2)
        # Initial guesses: [A, mu, sigma, a, b, c]
        p0_11 = [counts11.max(), 0.0, 0.02, 0, 0, np.median(counts11)]
        p0_13 = [counts13.max(), 0.0, 0.02, 0, 0, np.median(counts13)]
        try:
            popt11, _ = curve_fit(gauss_quad,
                                  centers[fit_mask],
                                  counts11[fit_mask],
                                  p0=p0_11)
        except:
            popt11 = [np.nan]*6
        try:
            popt13, _ = curve_fit(gauss_quad,
                                  centers[fit_mask],
                                  counts13[fit_mask],
                                  p0=p0_13)
        except:
            popt13 = [np.nan]*6
        # Extract mu and sigma
        mu11, sigma11 = popt11[1], abs(popt11[2])
        mu13, sigma13 = popt13[1], abs(popt13[2])
        mu11_list.append(mu11)
        sigma11_list.append(sigma11)
        mu13_list.append(mu13)
        sigma13_list.append(sigma13)
        # Plot histograms with error bars
        ax.errorbar(centers,
                    counts11,
                    yerr=errs11,
                    fmt='o',
                    color='blue',
                    label=f'cj 11.2.0, μ={mu11:.3f}, σ={sigma11:.3f}')
        ax.errorbar(centers + 0.001,
                    counts13,
                    yerr=errs13,
                    fmt='o',
                    color='red',
                    label=f'cj 13.0.3, μ={mu13:.3f}, σ={sigma13:.3f}')
        # Overlay fits
        ax.plot(centers, gauss_quad(centers, *popt11), '--', color='blue')
        ax.plot(centers, gauss_quad(centers, *popt13), '--', color='red')
        # Labels and title with suffix
        ax.set_xlabel(r'$M_{x}^{2}$ (GeV$^{2}$)')
        ax.set_ylabel('counts')
        ax.set_title(f'$e_θ$ ∈ [{low_deg}°, {high_deg}°] {title_suffix}')
        ax.legend()

    # Final subplot: μ vs e_theta with σ error bars, no connecting line, second offset
    ax = axes_flat[5]
    ax.errorbar(angle_means,
                mu11_list,
                yerr=sigma11_list,
                fmt='o',
                color='blue',
                label='cj 11.2.0')
    ax.errorbar(np.array(angle_means) + 0.1,
                mu13_list,
                yerr=sigma13_list,
                fmt='s',
                color='red',
                label='cj 13.0.3')
    # Horizontal line at y=0
    ax.axhline(0.0, color='grey', linestyle='--', linewidth=1)
    ax.set_xlabel(r'$e_θ$ (deg)')
    ax.set_ylabel(r'$\mu$ (GeV$^{2}$)')
    ax.set_ylim(-0.2, 0.2)
    ax.legend()

    fig.tight_layout()
    fig.savefig(f"/u/home/thayward/{outname}")