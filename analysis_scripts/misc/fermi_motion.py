#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# -----------------------------------------------------------------------------
# File paths for RGA, RGC, and RGB datasets
# -----------------------------------------------------------------------------
rga_epiX_atRest    = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_atRest.root"
# rga_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_fermiMotion.root"
rga_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_dilutionFactor.root"
rgc_epiX_atRest    = "/volatile/clas12/thayward/fermi_motion/rgc_su22_inb_epi+_atRest.root"

rga_epiPipi_atRest = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+pi-_atRest.root"
# rga_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_pi-_fermiMotion.root"
rga_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_pi-_dilutionFactor.root"
rgc_epiPipi_atRest = "/volatile/clas12/thayward/fermi_motion/rgc_su22_inb_epi+pi-_atRest.root"

rgb_epiX_atRest    = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+_atRest.root"
rgb_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+_fermiMotion.root"
rgb_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+_dilutionFactor.root"
# rgb_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+_fermiMotion_scaled30.root"
# rgb_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+_fermiMotion_nitrogenOnly.root"
rgb_epiPipi_atRest = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+pi-_atRest.root"
rgb_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+pi-_fermiMotion.root"
rgb_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+pi-_dilutionFactor.root"
# rgb_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+pi-_fermiMotion_scaled30.root"
# rgb_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+pi-_fermiMotion_nitrogenOnly.root"

# Ensure output directory exists
os.makedirs("output", exist_ok=True)

def load_array(path, branch):
    """Load a branch array from the PhysicsEvents tree."""
    with uproot.open(path) as f:
        return f["PhysicsEvents"][branch].array(library="np")

# -----------------------------------------------------------------------------
# Version lists for each panel
# -----------------------------------------------------------------------------
versions_epiX        = [
    ("RGA Fa18",                   rga_epiX_atRest),
    ("RGA Fa18 sim. Fermi Motion", rga_epiX_fermi),
    ("RGC Su22",                   rgc_epiX_atRest),
]

versions_epiPipiX    = [
    ("RGA Fa18",                   rga_epiPipi_atRest),
    ("RGA Fa18 sim. Fermi Motion", rga_epiPipi_fermi),
    ("RGC Su22",                   rgc_epiPipi_atRest),
]

versions_rgb_epiX    = [
    ("RGB Sp19",                   rgb_epiX_atRest),
    ("RGB Sp19 sim. Fermi Motion", rgb_epiX_fermi),
    ("RGC Su22",                   rgc_epiX_atRest),
]

versions_rgb_epiPipiX = [
    ("RGB Sp19",                   rgb_epiPipi_atRest),
    ("RGB Sp19 sim. Fermi Motion", rgb_epiPipi_fermi),
    ("RGC Su22",                   rgc_epiPipi_atRest),
]

# -----------------------------------------------------------------------------
# Fit function: Gaussian + quadratic background
# -----------------------------------------------------------------------------
def gauss_quad(x, A, mu, sigma, a0, a1, a2):
    return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + a0 + a1*x + a2*x**2

# -----------------------------------------------------------------------------
# Histogram settings
# -----------------------------------------------------------------------------
bins      = np.linspace(-1, 3, 101)                 # 100 bins from -1 to 3
bin_width = bins[1] - bins[0]
m_p2      = 0.93827**2                              # proton mass squared for initial guess
colors    = plt.rcParams['axes.prop_cycle'].by_key()['color']

# -----------------------------------------------------------------------------
# Prepare a 2×2 figure
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

panel_configs = [
    (axes[0,0], versions_epiX,         r"$e\,\pi^{+}X:\ M_x^2$"),
    (axes[0,1], versions_epiPipiX,     r"$e\,\pi^{+}\pi^{-}X:\ M_x^2$"),
    (axes[1,0], versions_rgb_epiX,     r"$e\,\pi^{+}X\ (\mathrm{RGB}):\ M_x^2$"),
    (axes[1,1], versions_rgb_epiPipiX, r"$e\,\pi^{+}\pi^{-}X\ (\mathrm{RGB}):\ M_x^2$"),
]

for ax, versions, title in panel_configs:
    for idx, (lbl, path) in enumerate(versions):
        # Load data and histogram
        data   = load_array(path, "Mx2")
        counts, _ = np.histogram(data, bins=bins)
        centers   = 0.5 * (bins[:-1] + bins[1:])
        norm = 1.0/(counts.sum() * bin_width)
        density = counts * norm
        errors    = np.sqrt(counts) * norm

        # Fit range [0.4, 1.2]
        mask   = (centers >= 0.4) & (centers <= 1.2)
        xfit   = centers[mask]
        yfit   = density[mask]
        errfit = errors[mask]

        # Initial guesses and bounds
        p0     = [yfit.max(), m_p2, 0.02, 0.0, 0.0, 0.0]
        bounds = ([0.0, 0.0, 0.0, -np.inf, -np.inf, -np.inf],
                  [np.inf, np.inf, np.inf,  np.inf,  np.inf,  np.inf])

        popt, _ = curve_fit(gauss_quad, xfit, yfit,
                            p0=p0, sigma=errfit, bounds=bounds)
        mu, sigma = popt[1], popt[2]

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

    ax.set_xlim(-1, 3)
    ax.set_xlabel(r"$M_x^2\ \mathrm{(GeV^2)}$")
    ax.set_title(title)
    ax.legend()

plt.tight_layout()
plt.subplots_adjust(bottom=0.15, hspace=0.3)

# Save under original filename
plt.savefig("output/fermi_motion.pdf")
plt.close()

# -----------------------------------------------------------------------------
# Plot dilution factor vs Mx2
# -----------------------------------------------------------------------------
# bin edges for Mx2
mx2_edges   = np.array([0.2, 0.4, 0.6, 0.7, 0.8, 0.9,
                        1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
                        1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0])
mx2_centers = 0.5 * (mx2_edges[:-1] + mx2_edges[1:])
dilution    = np.array([
    0.0875895, 0.0858209, 0.112782, 0.180683,
    0.322474, 0.303501, 0.197981, 0.185084,
    0.181473, 0.173468, 0.203505, 0.188717,
    0.174217, 0.177695, 0.184985, 0.171676,
    0.190922, 0.208126
])
dil_err     = np.array([
    0.0131908, 0.0117119, 0.0138437, 0.0107054,
    0.00696578,0.00666437,0.00804998,0.00804589,
    0.00781989,0.00762175,0.00689467,0.00679758,
    0.00676546,0.00652095,0.00622921,0.0061906,
    0.00240898,0.00198573
])

plt.figure()
plt.errorbar(mx2_centers, dilution, yerr=dil_err,
             fmt='o', markersize=5, color='k')
plt.xlabel(r"$M_{x}^{2}\ \mathrm{(GeV^2)}$")
plt.ylabel(r"$D_{f}$")
plt.tight_layout()
plt.savefig("output/dilution_factor.pdf")

# -----------------------------------------------------------------------------
# 9th-order polynomial fit to dilution factor
# -----------------------------------------------------------------------------
coeffs = np.polyfit(mx2_centers, dilution, 9)
poly9  = np.poly1d(coeffs)

# -----------------------------------------------------------------------------
# 9th-order polynomial fit to dilution factor
# -----------------------------------------------------------------------------
coeffs = np.polyfit(mx2_centers, dilution, 9)
poly9  = np.poly1d(coeffs)

# print the fit function in the terminal
print("9th-order polynomial fit for dilution factor D_f(Mx2):")
print(poly9)

# overlay the fit on the dilution‐factor plot
x_fit = np.linspace(mx2_centers.min(), mx2_centers.max(), 300)
y_fit = poly9(x_fit)
plt.figure()
plt.errorbar(mx2_centers, dilution, yerr=dil_err,
             fmt='o', markersize=5, color='k', label='data')
plt.plot(x_fit, y_fit, '-', linewidth=1.5, label='9th-order poly fit')
plt.xlabel(r"$M_{x}^{2}\ \mathrm{(GeV^2)}$")
plt.ylabel(r"$D_{f}$")
plt.legend()
plt.tight_layout()
plt.savefig("output/dilution_factor_with_fit.pdf")

# print the fit function in the terminal
print("9th-order polynomial fit for dilution factor D_f(Mx2):")
print(poly9)

plt.close()