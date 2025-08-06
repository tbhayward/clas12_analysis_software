#!/usr/bin/env python3
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Proton mass squared (GeV^2)
M_p2 = 0.93827**2

# Gaussian + cubic background model
def gauss_cubic(x, A, mu, sigma, a, b, c, d):
    return (A * np.exp(-(x - mu)**2 / (2 * sigma**2))
            + a*x**3 + b*x**2 + c*x + d)

# Files to use
file11 = "/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root"
file13 = "/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root"

# Branches: (branch name, LaTeX label, output filename)
angle_vars = [
    ("e_theta",  r"$e_θ$",          "dipion_e_theta.pdf"),
    ("p1_theta", r"$\theta_{\pi+}$","dipion_p1_theta.pdf"),
    ("p2_theta", r"$\theta_{\pi-}$","dipion_p2_theta.pdf"),
]

# e_theta and pion theta bin edges (degrees)
edges_e = [5, 14, 18, 22, 26, 30]
edges_pi = [5, 15, 21, 27, 33, 40]

# Histogram settings for Mx2: from 0.80 to 1.10 in 41 bins
bins_M  = np.linspace(0.80, 1.10, 41)
centers = 0.5 * (bins_M[:-1] + bins_M[1:])

# Fit parameter bounds: A ≥ 0; μ ∈ [0.8,0.95]; σ ≥ 0; cubic coeffs free
lower_bounds = [0.0, 0.8, 0.0, -np.inf, -np.inf, -np.inf, -np.inf]
upper_bounds = [np.inf,0.95,np.inf,  np.inf,  np.inf,  np.inf,  np.inf]

# Load trees once
tree11 = uproot.open(file11)["PhysicsEvents"]
tree13 = uproot.open(file13)["PhysicsEvents"]

# Preload arrays
Mx2_11  = tree11["Mx2"].array(library="np")
Mx2_13  = tree13["Mx2"].array(library="np")
det1_11 = tree11["detector1"].array(library="np")
det2_11 = tree11["detector2"].array(library="np")
det1_13 = tree13["detector1"].array(library="np")
det2_13 = tree13["detector2"].array(library="np")

for branch, angle_label, outname in angle_vars:
    theta11 = tree11[branch].array(library="np") * (180.0/np.pi)
    theta13 = tree13[branch].array(library="np") * (180.0/np.pi)

    # choose bins per branch
    edges = edges_e if branch == "e_theta" else edges_pi
    bins_deg = list(zip(edges[:-1], edges[1:]))

    angle_means  = []
    mu11_list    = []
    sigma11_list = []
    mu13_list    = []
    sigma13_list = []

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes_flat = axes.flatten()

    for i, (low_deg, high_deg) in enumerate(bins_deg):
        ax = axes_flat[i]

        # apply angle + detector selection
        mask11 = ((theta11 >= low_deg) & (theta11 < high_deg)
                  & (det1_11 == 1) & (det2_11 == 1))
        mask13 = ((theta13 >= low_deg) & (theta13 < high_deg)
                  & (det1_13 == 1) & (det2_13 == 1))

        M11 = Mx2_11[mask11]
        M13 = Mx2_13[mask13]

        angle_mean = np.mean(np.concatenate([theta11[mask11],
                                             theta13[mask13]]))
        angle_means.append(angle_mean)

        counts11, _ = np.histogram(M11, bins=bins_M)
        counts13, _ = np.histogram(M13, bins=bins_M)
        errs11 = np.sqrt(counts11)
        errs13 = np.sqrt(counts13)

        fit_mask = (centers >= 0.80) & (centers <= 1.10)
        p0_11 = [counts11.max(), M_p2, 0.02, 0, 0, 0, np.median(counts11)]
        p0_13 = [counts13.max(), M_p2, 0.02, 0, 0, 0, np.median(counts13)]

        try:
            popt11, _ = curve_fit(
                gauss_cubic,
                centers[fit_mask],
                counts11[fit_mask],
                p0=p0_11,
                bounds=(lower_bounds, upper_bounds)
            )
        except:
            popt11 = [np.nan]*7

        try:
            popt13, _ = curve_fit(
                gauss_cubic,
                centers[fit_mask],
                counts13[fit_mask],
                p0=p0_13,
                bounds=(lower_bounds, upper_bounds)
            )
        except:
            popt13 = [np.nan]*7

        mu11, sigma11 = popt11[1], abs(popt11[2])
        mu13, sigma13 = popt13[1], abs(popt13[2])
        mu11_list.append(mu11);    sigma11_list.append(sigma11)
        mu13_list.append(mu13);    sigma13_list.append(sigma13)

        ax.errorbar(centers, counts11, yerr=errs11,
                    fmt='o', color='blue',
                    label=f'cj11.2.0, μ={mu11:.3f}, σ={sigma11:.3f}')
        ax.errorbar(centers+0.001, counts13, yerr=errs13,
                    fmt='o', color='red',
                    label=f'cj13.0.3, μ={mu13:.3f}, σ={sigma13:.3f}')

        ax.plot(centers, gauss_cubic(centers, *popt11),
                '--', color='blue')
        ax.plot(centers, gauss_cubic(centers, *popt13),
                '--', color='red')

        ax.set_xlabel(r'$M_{x}^{2}$ (GeV$^{2}$)')
        ax.set_ylabel('counts')
        ax.set_title(f'{angle_label} ∈ [{low_deg}°, {high_deg}°]')
        ax.legend()

    # final μ vs angle
    ax = axes_flat[5]
    ax.errorbar(angle_means, mu11_list, yerr=sigma11_list,
                fmt='o', color='blue', label='cj11.2.0')
    ax.errorbar(np.array(angle_means)+0.1, mu13_list, yerr=sigma13_list,
                fmt='s', color='red',  label='cj13.0.3')
    ax.axhline(M_p2, color='grey', linestyle='--', linewidth=1)
    ax.set_xlabel(f'{angle_label} (deg)')
    ax.set_ylabel(r'$\mu$ (GeV$^{2}$)')
    ax.set_ylim(0.8, 1.1)
    ax.legend()

    fig.tight_layout()
    fig.savefig(f"/u/home/thayward/{outname}")