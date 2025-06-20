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
rgc_epiX_atRest    = "/volatile/clas12/thayward/fermi_motion/rgc_su22_inb_epi+_atRest.root"
rgc_epiPipi_atRest = "/volatile/clas12/thayward/fermi_motion/rgc_su22_inb_epi+pi-_atRest.root"

rga_epiX_atRest    = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_atRest.root"
rga_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_dilutionFactor.root"

rga_epiPipi_atRest = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+pi-_atRest.root"
rga_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+pi-_dilutionFactor.root"

rgb_epiX_atRest    = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+_atRest.root"
rgb_epiX_fermi     = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+_dilutionFactor.root"

rgb_epiPipi_atRest = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+pi-_atRest.root"
rgb_epiPipi_fermi  = "/volatile/clas12/thayward/fermi_motion/rgb_sp19_inb_epi+pi-_dilutionFactor.root"

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
bins      = np.linspace(-2, 4, 101)   # 100 bins from -1 to 3
bin_width = bins[1] - bins[0]
m_p2      = 0.93827**2                # proton mass squared for initial guess
colors    = plt.rcParams['axes.prop_cycle'].by_key()['color']

# -----------------------------------------------------------------------------
# Prepare a 2×2 figure
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# panel_configs entries: (axis, versions, title, low, high)
panel_configs = [
    (axes[0,0], versions_epiX,        r"$e\,\pi^{+}X:\ M_x^2$",        0.653, 1.145),
    (axes[0,1], versions_epiPipiX,    r"$e\,\pi^{+}\pi^{-}X:\ M_x^2$",  0.685, 1.141),
    (axes[1,0], versions_rgb_epiX,    r"$e\,\pi^{+}X\ (\mathrm{RGB}):\ M_x^2$",        0.653, 1.145),
    (axes[1,1], versions_rgb_epiPipiX,r"$e\,\pi^{+}\pi^{-}X\ (\mathrm{RGB}):\ M_x^2$",  0.685, 1.141),
]

for ax, versions, title, low, high in panel_configs:
    # pre-load at-rest evnums and Mx2 once per panel
    at_rest_path = versions[0][1]
    mx2_at       = load_array(at_rest_path, "Mx2")
    ev_at        = load_array(at_rest_path, "evnum")

    for idx, (lbl, path) in enumerate(versions):
        # Load data and histogram
        mx2_f, ev_f = load_array(path, "Mx2"), load_array(path, "evnum")
        counts, _   = np.histogram(mx2_f, bins=bins)
        centers     = 0.5 * (bins[:-1] + bins[1:])
        norm        = 1.0/(counts.sum() * bin_width)
        density     = counts * norm
        errors      = np.sqrt(counts) * norm

        # Fit range [0.4, 1.2]
        mask   = (centers >= 0.4) & (centers <= 1.2)
        xfit   = centers[mask]
        yfit   = density[mask]
        errfit = errors[mask]

        # Initial guesses and bounds for Gaussian+quad
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

        # For the simulated-Fermi case, find events that moved into [low,high]
        if idx == 1:
            in_peak_f  = (mx2_f >= low) & (mx2_f <= high)
            ev_in_f    = ev_f[in_peak_f]
            # those evnos in fermi-peak that were outside at rest
            mask_at_out = np.isin(ev_at, ev_in_f) & ((mx2_at < low) | (mx2_at > high))
            ev_moved    = ev_at[mask_at_out]
            mx2_start   = mx2_at[mask_at_out]

            if mx2_start.size > 0:
                counts_s, _ = np.histogram(mx2_start, bins=bins)
                density_s   = counts_s * norm
                # plot original positions as dashed gray
                ax.plot(centers, density_s,
                        linestyle='--', linewidth=1.0, color='gray',
                        label=f"{lbl} orig→peak")

    ax.set_xlim(-2, 4)
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
             fmt='o', markersize=5, color='k', label='data')
plt.xlabel(r"$M_{x}^{2}\ \mathrm{(GeV^2)}$")
plt.ylabel(r"$D_{f}$")

# -----------------------------------------------------------------------------
# 4th-order polynomial + Gaussian fit to dilution factor
# -----------------------------------------------------------------------------
def poly4_gauss(x, a0, a1, a2, a3, a4, A, mu, sigma):
    return (a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4
            + A * np.exp(- (x - mu)**2 / (2 * sigma**2)))

# initial guess from 4th-order polyfit
poly4_coeffs = np.polyfit(mx2_centers, dilution, 4)
p0 = list(poly4_coeffs[::-1]) + [0.1, 0.9, 0.1]
bounds = (
    [-np.inf]*5 + [0.0, 0.7, 0.0],
    [ np.inf]*5 + [np.inf, 1.0, np.inf]
)

popt, _ = curve_fit(poly4_gauss, mx2_centers, dilution,
                    sigma=dil_err, p0=p0, bounds=bounds)

print("4th-order poly + Gaussian fit parameters for D_f(Mx2):")
print(f"a0 = {popt[0]:.6g}")
print(f"a1 = {popt[1]:.6g}")
print(f"a2 = {popt[2]:.6g}")
print(f"a3 = {popt[3]:.6g}")
print(f"a4 = {popt[4]:.6g}")
print(f"A  = {popt[5]:.6g}")
print(f"mu = {popt[6]:.6g}")
print(f"sigma = {popt[7]:.6g}")

# overlay the fit
x_fit = np.linspace(mx2_centers.min(), mx2_centers.max(), 300)
y_fit = poly4_gauss(x_fit, *popt)
plt.plot(x_fit, y_fit, '-', linewidth=1.5, label='4th-order poly + Gauss fit')
plt.legend()

plt.tight_layout()
plt.savefig("output/dilution_factor.pdf")
plt.close()



# -----------------------------------------------------------------------------
# Compute and plot, per Mx2 bin, the fraction of events that were
# outside the peak at rest but end up inside the peak after Fermi smearing
# with a uniform y-axis limit = 1.2×global max
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

panel_configs = [
    (axes[0,0], versions_epiX,         r"$e\,\pi^{+}X:\ M_x^2$",        0.653, 1.145),
    (axes[0,1], versions_epiPipiX,     r"$e\,\pi^{+}\pi^{-}X:\ M_x^2$",  0.685, 1.141),
    (axes[1,0], versions_rgb_epiX,     r"$e\,\pi^{+}X\ (\mathrm{RGB}):\ M_x^2$",        0.653, 1.145),
    (axes[1,1], versions_rgb_epiPipiX, r"$e\,\pi^{+}\pi^{-}X\ (\mathrm{RGB}):\ M_x^2$",  0.685, 1.141),
]

# bin centers for plotting
centers = 0.5 * (bins[:-1] + bins[1:])

# first pass: compute all fractions and track the global maximum
fracs_list = []
global_max = 0.0

for ax, versions, title, low, high in panel_configs:
    # at-rest values
    at_path = versions[0][1]
    mx2_at  = load_array(at_path, "Mx2")
    ev_at   = load_array(at_path, "evnum")
    # smeared values
    f_path = versions[1][1]
    mx2_f  = load_array(f_path, "Mx2")
    ev_f   = load_array(f_path, "evnum")

    # identify events that end up inside the peak after smearing
    in_peak_f = (mx2_f >= low) & (mx2_f <= high)
    ev_in_f   = set(ev_f[in_peak_f])

    # find which of those started outside the peak
    mask_out_at = [((m < low or m > high) and (e in ev_in_f))
                   for m, e in zip(mx2_at, ev_at)]
    mx2_start = mx2_at[mask_out_at]

    # histogram of their starting bins
    counts_shift, _ = np.histogram(mx2_start, bins=bins)

    # total moved into peak
    N_peak = len(ev_in_f)
    frac   = counts_shift / N_peak if N_peak > 0 else np.zeros_like(counts_shift)

    fracs_list.append(frac)
    if frac.max() > global_max:
        global_max = frac.max()

# set uniform y-limit
y_lim = global_max * 1.2

# second pass: actually draw them
for (ax, versions, title, low, high), frac in zip(panel_configs, fracs_list):
    ax.plot(centers, frac, '-o', markersize=4, color='gray')
    ax.set_xlim(-2, 4)
    ax.set_ylim(0, y_lim)
    ax.set_xlabel(r"$M_x^2\ \mathrm{(GeV^2)}$")
    ax.set_ylabel("shifted fraction")
    ax.set_title(title)

plt.tight_layout()
plt.savefig("output/shifted_percentages.pdf")
plt.close()