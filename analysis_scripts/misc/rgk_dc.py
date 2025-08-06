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

# Branches to loop over
angle_vars = [
    ("e_theta",  r"$e_θ$",          "dipion_e_theta.pdf"),
    ("p1_theta", r"$\theta_{\pi+}$","dipion_p1_theta.pdf"),
    ("p2_theta", r"$\theta_{\pi-}$","dipion_p2_theta.pdf"),
]

# Bin definitions
edges_e  = [5, 12, 16, 20, 24, 30]
edges_pi = [5, 15, 21, 27, 33, 40]

# Mx² histogram: 0.80–1.10 with 40 bins → 41 edges
bins_M  = np.linspace(0.80, 1.05, 41)
centers = 0.5 * (bins_M[:-1] + bins_M[1:])

# Fit bounds: enforce 0.8 ≤ μ ≤ 0.95
lower_bounds = [0.0, 0.8,  0.0, -np.inf, -np.inf, -np.inf]
upper_bounds = [np.inf,0.95,np.inf,  np.inf,  np.inf,  np.inf]

# Open once
tree11 = uproot.open(file11)["PhysicsEvents"]
tree13 = uproot.open(file13)["PhysicsEvents"]

# Preload arrays
Mx2_11   = tree11["Mx2"].array(library="np")
Mx2_13   = tree13["Mx2"].array(library="np")
d1_11    = tree11["detector1"].array(library="np")
d2_11    = tree11["detector2"].array(library="np")
d1_13    = tree13["detector1"].array(library="np")
d2_13    = tree13["detector2"].array(library="np")

for branch, angle_label, outname in angle_vars:
    # load and convert to degrees
    th11 = tree11[branch].array(library="np") * (180/np.pi)
    th13 = tree13[branch].array(library="np") * (180/np.pi)

    # select bins depending on variable
    edges = edges_e if branch == "e_theta" else edges_pi
    bins_deg = list(zip(edges[:-1], edges[1:]))

    angle_means = []
    mu11_list, sigma11_list = [], []
    mu13_list, sigma13_list = [], []

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes_flat = axes.flatten()

    for i, (low, high) in enumerate(bins_deg):
        ax = axes_flat[i]

        # selection masks
        m11 = (th11>=low)&(th11<high)&(d1_11==1)&(d2_11==1)
        m13 = (th13>=low)&(th13<high)&(d1_13==1)&(d2_13==1)

        M11 = Mx2_11[m11]
        M13 = Mx2_13[m13]
        angle_means.append(np.mean(np.concatenate([th11[m11], th13[m13]])))

        c11, _ = np.histogram(M11, bins=bins_M)
        c13, _ = np.histogram(M13, bins=bins_M)
        e11 = np.sqrt(c11)
        e13 = np.sqrt(c13)

        fit_mask = (centers>=0.80)&(centers<=1.10)
        p0_11 = [c11.max(), M_p2, 0.02, 0, 0, np.median(c11)]
        p0_13 = [c13.max(), M_p2, 0.02, 0, 0, np.median(c13)]

        try:
            popt11, _ = curve_fit(
                gauss_quad,
                centers[fit_mask],
                c11[fit_mask],
                p0=p0_11,
                bounds=(lower_bounds, upper_bounds)
            )
        except:
            popt11 = [np.nan]*6
        try:
            popt13, _ = curve_fit(
                gauss_quad,
                centers[fit_mask],
                c13[fit_mask],
                p0=p0_13,
                bounds=(lower_bounds, upper_bounds)
            )
        except:
            popt13 = [np.nan]*6

        mu11, sigma11 = popt11[1], abs(popt11[2])
        mu13, sigma13 = popt13[1], abs(popt13[2])
        mu11_list.append(mu11); sigma11_list.append(sigma11)
        mu13_list.append(mu13); sigma13_list.append(sigma13)

        # data + errors
        ax.errorbar(centers, c11, yerr=e11, fmt='o', color='blue',
                    label=f'cj11.2.0, μ={mu11:.3f}, σ={sigma11:.3f}')
        ax.errorbar(centers+0.001, c13, yerr=e13, fmt='o', color='red',
                    label=f'cj13.0.3, μ={mu13:.3f}, σ={sigma13:.3f}')

        # fits
        ax.plot(centers, gauss_quad(centers, *popt11), '--', color='blue')
        ax.plot(centers, gauss_quad(centers, *popt13), '--', color='red')

        ax.set_xlim(0.80, 1.05)
        ax.set_ylim(0, np.max([c11.max(), c13.max()])*1.2)
        ax.set_xlabel(r'$M_{x}^{2}$ (GeV$^{2}$)')
        ax.set_ylabel('counts')
        ax.set_title(f'{angle_label} ∈ [{low}°, {high}°]')
        ax.legend()

    # final panel: μ vs angle
    ax = axes_flat[5]
    ax.errorbar(angle_means, mu11_list, yerr=sigma11_list,
                fmt='o', color='blue', label='cj11.2.0')
    ax.errorbar(np.array(angle_means)+0.1, mu13_list, yerr=sigma13_list,
                fmt='s', color='red',  label='cj13.0.3')
    ax.axhline(M_p2, linestyle='--', color='grey')
    ax.set_xlabel(f'{angle_label} (deg)')
    ax.set_ylabel(r'$\mu$ (GeV$^{2}$)')
    ax.set_ylim(0.8, 1.2)
    ax.legend()

    fig.tight_layout()
    fig.savefig(f"/u/home/thayward/{outname}")