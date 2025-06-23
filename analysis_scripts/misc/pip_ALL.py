#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Use LaTeX for all text, serif font for a scientific look
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

# HERMES "100%" and 68% bands for Pip, 0.4<z<0.5
hermes_data = [
    # x,    Min100,    Max100,    Min68,     Max68
    [0.125, 0.237413,  0.313158,  0.255901,  0.302772],
    [0.175, 0.323555,  0.392161,  0.334571,  0.381672],
    [0.225, 0.373040,  0.469157,  0.396504,  0.453788],
    [0.275, 0.405763,  0.555878,  0.454167,  0.527462],
    [0.325, 0.441736,  0.693701,  0.515448,  0.628330],
    [0.375, 0.446821,  0.813321,  0.595255,  0.744260],
    [0.425, 0.493661,  0.926293,  0.673446,  0.865753],
    [0.475, 0.498619,  0.951409,  0.676030,  0.939469],
    [0.525, 0.504007,  0.966825,  0.644757,  0.947681],
    [0.575, 0.509595,  1.004663,  0.635088,  0.957669],
    [0.625, 0.512737,  1.012946,  0.629967,  0.964844],
    [0.675, 0.480661,  1.004594,  0.632191,  0.943644]
]
hermes = np.array(hermes_data)
h_x      = hermes[:, 0]
h_min100 = hermes[:, 1]
h_max100 = hermes[:, 2]
h_min68  = hermes[:, 3]
h_max68  = hermes[:, 4]

# CLAS12 pip and pipp datasets
meanx = [0.1175, 0.1925, 0.2675, 0.3425, 0.4175, 0.4925, 0.5675, 0.6425]
pip = [
    [meanx[0], 0.2257, 0.0157],
    [meanx[1], 0.3294, 0.0134],
    [meanx[2], 0.4268, 0.0152],
    [meanx[3], 0.5141, 0.0204],
    [meanx[4], 0.5964, 0.0277],
    [meanx[5], 0.6494, 0.0416],
    [meanx[6], 0.7270, 0.0684],
    [meanx[7], 0.8008, 0.1570]
]
pipp = [
    [0.144881754, 0.279487685, 0.011047009],
    [0.223757471, 0.407836967, 0.010338185],
    [0.319545878, 0.561530186, 0.015820949],
    [0.416155475, 0.687906307, 0.029495067],
    [0.510505127, 0.877989016, 0.071471465]
]

def to_arrays(data):
    arr = np.array(data)
    return arr[:,0], arr[:,1], arr[:,2]

p_x,  p_y,  p_yerr  = to_arrays(pip)
pp_x, pp_y, pp_yerr = to_arrays(pipp)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))

# 1) CLAS12 error‐bar data
for x, y, yerr, label, color, marker in [
    (p_x,  p_y,  p_yerr,  r"CLAS12 preliminary, $e\pi^{+}X$",   "red",   "o"),
    (pp_x, pp_y, pp_yerr, r"CLAS12 preliminary, $e\pi^{+}pX$",  "green", "^"),
]:
    ax.errorbar(
        x, y, yerr=yerr,
        fmt=marker,
        color=color,
        ecolor=color,
        capsize=4,
        elinewidth=1.5,
        markeredgewidth=1.5,
        markersize=6,
        label=label
    )
#endfor

# 2) HERMES bands
#    - lighter blue: 100% band
ax.fill_between(
    h_x, h_min100, h_max100,
    color='blue', alpha=0.2,
    label=r"HERMES $e\pi^{+}X$ 100\% CL"
)
#    - darker blue: 68% band
ax.fill_between(
    h_x, h_min68, h_max68,
    color='blue', alpha=0.4,
    label=r"HERMES $e\pi^{+}X$ 68\% CL"
)

# Axes, labels, limits
ax.set_xlabel(r"$x_{B}$",         fontsize=16)
ax.set_ylabel(r"$F_{LL}/F_{UU}$", fontsize=16)
ax.set_xlim(0.06, 0.7)
ax.set_ylim(0.0, 1.0)

# Grid and legend
ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
ax.legend(frameon=True, fontsize=14)

# Save
plt.tight_layout()
plt.savefig("/u/home/thayward/pip_ALL.pdf")