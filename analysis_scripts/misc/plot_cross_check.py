#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2x2 cross-check plot for F_UL/F_UU:
- Row 1: low x_B   (Avakian: xbin-1.dat) vs Hayward low-xB results
- Row 2: high x_B  (Avakian: xbin-4.dat) vs Hayward high-xB results
- Col 1: sin(phi) harmonic
- Col 2: sin(2phi) harmonic

Saves to: output/enpi+/cross_check.pdf
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# --------------------------
# Avakian (.dat-derived) values from previous calculation:
# Triples are [<t>, F_UL/F_UU, err]
# --------------------------

# Low x_B (xbin-1.dat)
avak_low_sin = [
    [0.1329,  0.0385, 0.0140],
    [0.3377, -0.1753, 0.0224],
    [0.5421, -0.2282, 0.0294],
    [0.7419, -0.0934, 0.0372],
    [0.9445, -0.0660, 0.0419],
    [1.1440,  0.1415, 0.0482],
]

avak_low_sin2 = [
    [0.1329,  0.0143, 0.0282],
    [0.3377, -0.2929, 0.0463],
    [0.5421, -0.3848, 0.0641],
    [0.7419, -0.3954, 0.0814],
    [0.9445, -0.2696, 0.0920],
    [1.1440, -0.0777, 0.1072],
]

# High x_B (xbin-4.dat)
avak_high_sin = [
    [0.3959,  0.0942, 0.0288],
    [0.5501,  0.0728, 0.0184],
    [0.7454,  0.0801, 0.0192],
    [0.9454,  0.0848, 0.0219],
    [1.1450,  0.0882, 0.0254],
]

avak_high_sin2 = [
    [0.3959,  0.1318, 0.0602],
    [0.5501, -0.1831, 0.0391],
    [0.7454, -0.1583, 0.0409],
    [0.9454, -0.0929, 0.0475],
    [1.1450, -0.2234, 0.0545],
]

# --------------------------
# Hayward personal results (triples are [-t, value, err])
# Convert -t to +t for the horizontal axis labeled "-t (GeV^{2})"
# --------------------------

hay_low_sin = [
    [-1.145094716,  0.051887670, 0.011028437],
    [-0.945014183, -0.013012838, 0.008348541],
    [-0.744714118, -0.041858839, 0.008685478],
    [-0.542435981, -0.057690416, 0.009943711],
    [-0.338121635, -0.042009072, 0.006112725],
    [-0.136911649,  0.006336992, 0.003772195],
]

hay_low_sin2 = [
    [-1.145158888, -0.096921098, 0.031273345],
    [-0.944707391, -0.076681275, 0.036980092],
    [-0.743865259, -0.144912815, 0.025745135],
    [-0.542550478, -0.065932329, 0.028436012],
    [-0.338126016, -0.069212424, 0.014637813],
    [-0.137020812, -0.013684392, 0.009972897],
]

hay_high_sin = [
    [-1.145514377,  0.028795432, 0.007978569],
    [-0.945539672,  0.045128505, 0.005936367],
    [-0.746257649,  0.034328296, 0.005664333],
    [-0.551304678,  0.022571317, 0.005742487],
    [-0.380268267,  0.028183899, 0.007550841],
]

hay_high_sin2 = [
    [-1.145548835, -0.084812590, 0.013071746],
    [-0.945571596, -0.079754103, 0.012065463],
    [-0.746230092, -0.045215948, 0.012489881],
    [-0.551421926, -0.038356521, 0.009811295],
    [-0.380486537,  0.008492162, 0.014596639],
]

# --------------------------
# Helpers
# --------------------------
def to_arrays(triples, negate_first=False):
    """
    Convert list of [x, y, err] to numpy arrays.
    If negate_first is True, x_out = -triples[i][0] (used to convert -t to +t).
    """
    triples = np.asarray(triples, dtype=float)
    x = triples[:, 0]
    if negate_first:
        x = -x
    # endif
    y = triples[:, 1]
    e = triples[:, 2]
    return x, y, e

def plot_panel(ax, avak_triples, hay_triples, ylabel):
    """
    Plot Avakian vs Hayward datasets as markers with error bars.
    No connecting lines.
    """
    xa, ya, ea = to_arrays(avak_triples, negate_first=False)
    xh, yh, eh = to_arrays(hay_triples, negate_first=True)  # convert -t to +t

    # Avakian
    ax.errorbar(
        xa, ya, yerr=ea, fmt='s', linestyle='none', capsize=3, label='Avakian'
    )
    # Hayward
    ax.errorbar(
        xh, yh, yerr=eh, fmt='o', linestyle='none', capsize=3, label='Hayward'
    )

    ax.set_xlabel(r"$-t$ (GeV^{2})")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False)
    ax.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)
    # Optional x-range guard based on available points
    xmin = 0.0
    xmax = max(1.2, np.nanmax([xa.max() if xa.size else 0.0, xh.max() if xh.size else 0.0]) * 1.05)
    ax.set_xlim(xmin, xmax)
# endfor

# --------------------------
# Figure
# --------------------------
fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)

# Titles to clarify each subplot (optional but helpful)
axes[0, 0].set_title("Low x_B: sin(phi)")
axes[0, 1].set_title("Low x_B: sin(2phi)")
axes[1, 0].set_title("High x_B: sin(phi)")
axes[1, 1].set_title("High x_B: sin(2phi)")

# Requested y-axis labels (user specified same label text for both columns)
ylabel_sin  = r"$F_{UL}^{\sin\phi}/F_{UU}$"
ylabel_sin2 = r"$F_{UL}^{\sin\phi}/F_{UU}$"  # as requested

# Row 1: Low x_B
plot_panel(axes[0, 0], avak_low_sin,  hay_low_sin,  ylabel_sin)
plot_panel(axes[0, 1], avak_low_sin2, hay_low_sin2, ylabel_sin2)

# Row 2: High x_B
plot_panel(axes[1, 0], avak_high_sin,  hay_high_sin,  ylabel_sin)
plot_panel(axes[1, 1], avak_high_sin2, hay_high_sin2, ylabel_sin2)

# --------------------------
# Save
# --------------------------
outdir = Path("output/enpi+")
outdir.mkdir(parents=True, exist_ok=True)
outfile = outdir / "cross_check.pdf"
fig.savefig(outfile, dpi=300)
print(f"Saved: {outfile}")