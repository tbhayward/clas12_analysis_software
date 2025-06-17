#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# -----------------------------------------------------------------------------
# File paths for RGA and RGC datasets
# -----------------------------------------------------------------------------
rga_epiX_atRest    = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_atRest.root"
rga_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_fermiMotion.root"
rgc_epiX_atRest    = "/volatile/clas12/thayward/fermi_motion/rgc_su22_inb_epi+_atRest.root"

rga_epiPipi_atRest = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+pi-_atRest.root"
rga_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_pi-_fermiMotion.root"
rgc_epiPipi_atRest = "/volatile/clas12/thayward/fermi_motion/rgc_su22_inb_epi+pi-_atRest.root"

# Ensure output directory exists
os.makedirs("output", exist_ok=True)

def load_array(path, branch):
    """Load a branch array from the PhysicsEvents tree."""
    with uproot.open(path) as f:
        return f["PhysicsEvents"][branch].array(library="np")

# -----------------------------------------------------------------------------
# Version lists for each channel
# -----------------------------------------------------------------------------
versions_epiX = [
    ("RGA Fa18",                   rga_epiX_atRest),
    ("RGA Fa18 sim. Fermi Motion", rga_epiX_fermi),
    ("RGC Su22",                   rgc_epiX_atRest),
]

versions_epiPipiX = [
    ("RGA Fa18",                   rga_epiPipi_atRest),
    ("RGA Fa18 sim. Fermi Motion", rga_epiPipi_fermi),
    ("RGC Su22",                   rgc_epiPipi_atRest),
]

# -----------------------------------------------------------------------------
# Fit function: Gaussian + quadratic background
# -----------------------------------------------------------------------------
def gauss_quad(x, A, mu, sigma, a0, a1, a2):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + a0 + a1*x + a2*x**2

# -----------------------------------------------------------------------------
# Histogram settings
# -----------------------------------------------------------------------------
bins      = np.linspace(-1, 3, 101)                 # 100 bins from -1 to 3
bin_width = bins[1] - bins[0]
m_p2       = 0.93827**2                              # proton mass squared for initial guess
colors     = plt.rcParams['axes.prop_cycle'].by_key()['color']

# -----------------------------------------------------------------------------
# Prepare figure with two panels
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# -----------------------------------------------------------------------------
# Loop over panels: (axes[0], versions_epiX) and (axes[1], versions_epiPipiX)
# -----------------------------------------------------------------------------
for ax, versions in zip(axes, (versions_epiX, versions_epiPipiX)):
    for idx, (lbl, path) in enumerate(versions):
        # Load and histogram data
        data   = load_array(path, "Mx2")
        counts, _ = np.histogram(data, bins=bins)
        centers   = 0.5 * (bins[:-1] + bins[1:])
        norm      = 1.0 / (data.size * bin_width)
        density   = counts * norm
        errors    = np.sqrt(counts) * norm

        # Mask for fit range [0.4, 1.2]
        mask   = (centers >= 0.4) & (centers <= 1.2)
        xfit   = centers[mask]
        yfit   = density[mask]
        errfit = errors[mask]

        # Initial guesses and bounds
        p0     = [yfit.max(), m_p2, 0.02, 0.0, 0.0, 0.0]
        bounds = ([0.0, 0.0, 0.0, -np.inf, -np.inf, -np.inf],
                  [np.inf, np.inf, np.inf,  np.inf,  np.inf,  np.inf])

        # Perform the fit
        popt, _ = curve_fit(gauss_quad, xfit, yfit,
                            p0=p0, sigma=errfit, bounds=bounds)
        mu, sigma = popt[1], popt[2]

        # Select color
        color = colors[idx % len(colors)]
        label = rf"{lbl} ($\mu={mu:.3f},\,\sigma={sigma:.3f}$)"

        # Plot data points with error bars
        ax.errorbar(centers, density, yerr=errors,
                    fmt='o', markersize=3, color=color,
                    label=label)

        # Overlay fit as solid line only over [0.4,1.2]
        xcurve = np.linspace(0.4, 1.2, 200)
        ycurve = gauss_quad(xcurve, *popt)
        ax.plot(xcurve, ycurve,
                linestyle='-', linewidth=1.5, color=color)

    # Axes labels and limits
    ax.set_xlim(-1, 3)
    ax.set_xlabel(r"$M_x^2\ \mathrm{(GeV^2)}$")

# Titles and legends
axes[0].set_title(r"$e\,\pi^{+}X:\ M_x^2$")
axes[1].set_title(r"$e\,\pi^{+}\pi^{-}X:\ M_x^2$")
axes[0].legend()
axes[1].legend()

# Final layout adjustments
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

# Save under original filename
plt.savefig("output/fermi_motion.pdf")
plt.close()