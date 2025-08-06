#!/usr/bin/env python3
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Proton mass squared (GeV^2)
M_p2 = 0.93827**2

# Gaussian function for fitting
def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Open ROOT files and extract arrays
f11 = uproot.open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root")
f13 = uproot.open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root")

arrays11 = f11["PhysicsEvents"].arrays(["Mx2", "e_theta"], library="np")
arrays13 = f13["PhysicsEvents"].arrays(["Mx2", "e_theta"], library="np")

Mx2_11     = arrays11["Mx2"]
e_theta_11 = arrays11["e_theta"]
Mx2_13     = arrays13["Mx2"]
e_theta_13 = arrays13["e_theta"]

# Define e_theta bins
edges = [5, 10, 15, 20, 25, 40]
bins  = list(zip(edges[:-1], edges[1:]))

# Prepare lists to store fit results
angle_means  = []
mu11_list    = []
sigma11_list = []
mu13_list    = []
sigma13_list = []

# Set up a 2x3 canvas
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes_flat = axes.flatten()

for i, (low, high) in enumerate(bins):
    ax = axes_flat[i]
    # Select events in this e_theta range
    mask11 = (e_theta_11 >= low) & (e_theta_11 < high)
    mask13 = (e_theta_13 >= low) & (e_theta_13 < high)
    M11 = Mx2_11[mask11]
    M13 = Mx2_13[mask13]
    # Compute the actual mean e_theta for this bin
    angle_mean = np.mean(np.concatenate([e_theta_11[mask11], e_theta_13[mask13]]))
    angle_means.append(angle_mean)
    # Histogram Mx2 from 0.4 to 1.1 in 70 bins
    bins_M = np.linspace(0.4, 1.1, 71)
    centers = 0.5 * (bins_M[:-1] + bins_M[1:])
    counts11, _ = np.histogram(M11, bins=bins_M)
    counts13, _ = np.histogram(M13, bins=bins_M)
    errs11 = np.sqrt(counts11)
    errs13 = np.sqrt(counts13)
    # Fit only in [M_p2 - 0.1, M_p2 + 0.1]
    fit_mask = (centers >= M_p2 - 0.1) & (centers <= M_p2 + 0.1)
    # Fit for cj11
    p0 = [counts11[fit_mask].max(), M_p2, 0.02]
    try:
        popt11, _ = curve_fit(gauss, centers[fit_mask], counts11[fit_mask], p0=p0)
    except:
        popt11 = [np.nan, np.nan, np.nan]
    mu11, sigma11 = popt11[1], abs(popt11[2])
    mu11_list.append(mu11)
    sigma11_list.append(sigma11)
    # Fit for cj13
    p0 = [counts13[fit_mask].max(), M_p2, 0.02]
    try:
        popt13, _ = curve_fit(gauss, centers[fit_mask], counts13[fit_mask], p0=p0)
    except:
        popt13 = [np.nan, np.nan, np.nan]
    mu13, sigma13 = popt13[1], abs(popt13[2])
    mu13_list.append(mu13)
    sigma13_list.append(sigma13)
    # Plot data points with errors
    ax.errorbar(centers, counts11, yerr=errs11, fmt='o', color='black',
                label=f'cj 11.2.0 (N={len(M11)})')
    ax.errorbar(centers, counts13, yerr=errs13, fmt='o', color='red',
                label=f'cj 13.0.3 (N={len(M13)})')
    # Add mu and sigma to legend
    ax.legend(title=f'μ11={mu11:.3f}, σ11={sigma11:.3f}\nμ13={mu13:.3f}, σ13={sigma13:.3f}')
    ax.set_xlabel(r'$M_{x}^{2}$ (GeV$^{2}$)')
    ax.set_ylabel('counts')
#endfor

# Final subplot: mu and sigma vs e_theta
ax = axes_flat[5]
ax.plot(angle_means, mu11_list, 'o-',  label=r'$\mu$ cj 11.2.0')
ax.plot(angle_means, mu13_list, 's--', label=r'$\mu$ cj 13.0.3')
ax.plot(angle_means, sigma11_list, '^-',  label=r'$\sigma$ cj 11.2.0')
ax.plot(angle_means, sigma13_list, 'x--', label=r'$\sigma$ cj 13.0.3')
ax.set_xlabel(r'$e_\theta$ (deg)')
ax.set_ylabel('Fit parameters')
ax.legend()

fig.tight_layout()
fig.savefig("/u/home/thayward/dipion.pdf")