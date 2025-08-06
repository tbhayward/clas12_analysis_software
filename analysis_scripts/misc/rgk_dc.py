#!/usr/bin/env python3
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Proton mass squared (GeV^2)
M_p2 = 0.93827**2

# Gaussian + quadratic background model
def gauss_quad(x, A, mu, sigma, a, b, c):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + a*x**2 + b*x + c

# Files to use
file11 = "/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root"
file13 = "/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root"

# Branches to loop over: (branch name, LaTeX label, output filename)
angle_vars = [
    ("e_theta", r"$e_θ$", "dipion_e_theta.pdf"),
    ("p1_theta", r"$\theta_{\pi+}$", "dipion_p1_theta.pdf"),
    ("p2_theta", r"$\theta_{\pi-}$", "dipion_p2_theta.pdf"),
]

# e_theta bin edges (degrees)
edges_deg = [5, 14, 18, 22, 26, 30]
bins_deg  = list(zip(edges_deg[:-1], edges_deg[1:]))

# Histogram settings for Mx2
bins_M  = np.linspace(0.7, 1.1, 41)
centers = 0.5 * (bins_M[:-1] + bins_M[1:])

# Load trees once
f11 = uproot.open(file11)
f13 = uproot.open(file13)

# Common arrays pulled lazily for Mx2 and detector1
Mx2_11_all = f11["PhysicsEvents"]["Mx2"].array(library="np")
Mx2_13_all = f13["PhysicsEvents"]["Mx2"].array(library="np")
det1_11    = f11["PhysicsEvents"]["detector1"].array(library="np")
det1_13    = f13["PhysicsEvents"]["detector1"].array(library="np")

for branch, angle_label, outname in angle_vars:
    # Pull the angle branch and convert to degrees
    theta11 = f11["PhysicsEvents"][branch].array(library="np") * (180.0/np.pi)
    theta13 = f13["PhysicsEvents"][branch].array(library="np") * (180.0/np.pi)

    # Prepare storage for fit results
    angle_means  = []
    mu11_list    = []
    sigma11_list = []
    mu13_list    = []
    sigma13_list = []

    # Setup 2x3 canvas
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes_flat = axes.flatten()

    # Loop over e_theta bins for histograms
    for i, (low_deg, high_deg) in enumerate(bins_deg):
        ax = axes_flat[i]

        # Select events in angle bin
        mask11 = (theta11 >= low_deg) & (theta11 < high_deg)
        mask13 = (theta13 >= low_deg) & (theta13 < high_deg)

        M11 = Mx2_11_all[mask11]
        M13 = Mx2_13_all[mask13]

        # Compute mean angle in bin
        angle_mean = np.mean(np.concatenate([theta11[mask11], theta13[mask13]]))
        angle_means.append(angle_mean)

        # Build histograms and errors
        counts11, _ = np.histogram(M11, bins=bins_M)
        counts13, _ = np.histogram(M13, bins=bins_M)
        errs11 = np.sqrt(counts11)
        errs13 = np.sqrt(counts13)

        # Fit over full Mx2 range [0.7,1.1]
        fit_mask = (centers >= 0.7) & (centers <= 1.1)
        p0_11 = [counts11.max(), M_p2, 0.02, 0, 0, np.median(counts11)]
        p0_13 = [counts13.max(), M_p2, 0.02, 0, 0, np.median(counts13)]
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

        mu11, sigma11 = popt11[1], abs(popt11[2])
        mu13, sigma13 = popt13[1], abs(popt13[2])
        mu11_list.append(mu11)
        sigma11_list.append(sigma11)
        mu13_list.append(mu13)
        sigma13_list.append(sigma13)

        # Plot data + errors
        ax.errorbar(centers, counts11, yerr=errs11,
                    fmt='o', color='blue',
                    label=f'cj11.2.0, μ={mu11:.3f}, σ={sigma11:.3f}')
        ax.errorbar(centers+0.001, counts13, yerr=errs13,
                    fmt='o', color='red',
                    label=f'cj13.0.3, μ={mu13:.3f}, σ={sigma13:.3f}')

        # Overlay fits
        ax.plot(centers, gauss_quad(centers, *popt11),
                '--', color='blue')
        ax.plot(centers, gauss_quad(centers, *popt13),
                '--', color='red')

        ax.set_xlabel(r'$M_{x}^{2}$ (GeV$^{2}$)')
        ax.set_ylabel('counts')
        ax.set_title(f'{angle_label} ∈ [{low_deg}°, {high_deg}°]')
        ax.legend()

    # Final subplot: μ vs angle
    ax = axes_flat[5]
    ax.errorbar(angle_means, mu11_list, yerr=sigma11_list,
                fmt='o', color='blue', label='cj11.2.0')
    ax.errorbar(np.array(angle_means)+0.1, mu13_list, yerr=sigma13_list,
                fmt='s', color='red', label='cj13.0.3')
    ax.axhline(M_p2, color='grey', linestyle='--', linewidth=1)
    ax.set_xlabel(f'{angle_label} (deg)')
    ax.set_ylabel(r'$\mu$ (GeV$^{2}$)')
    ax.set_ylim(0.8, 1.1)
    ax.legend()

    fig.tight_layout()
    fig.savefig(f"/u/home/thayward/{outname}")