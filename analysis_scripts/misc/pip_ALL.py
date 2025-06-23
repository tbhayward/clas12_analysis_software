#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Use LaTeX for all text, serif font for a scientific look
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
pipHERMES = [
    [0.125, 0.2793365, 0.046871],
    [0.175, 0.3581215, 0.047101],
    [0.225, 0.4251460, 0.057284],
    [0.275, 0.4908145, 0.073295],
    [0.325, 0.5718890, 0.112882],
    [0.375, 0.6697575, 0.149005],
    [0.425, 0.7695995, 0.192307],
    [0.475, 0.8077495, 0.263439],
    [0.525, 0.7962190, 0.302924],
    [0.575, 0.7963785, 0.322581],
    [0.625, 0.7974055, 0.334877],
    [0.675, 0.7879175, 0.311453]
]

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

# Convert to numpy arrays for convenience
def to_arrays(data):
    arr = np.array(data)
    return arr[:,0], arr[:,1], arr[:,2]

h_x,  h_y,  h_yerr  = to_arrays(pipHERMES)
p_x,  p_y,  p_yerr  = to_arrays(pip)
pp_x, pp_y, pp_yerr = to_arrays(pipp)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))

for x, y, yerr, label, color, marker in [
    (p_x,  p_y,  p_yerr,  r"CLAS12 preliminary, $e\pi^{+}X$",  "red",   "o"),
    (pp_x, pp_y, pp_yerr, r"CLAS12 preliminary, $e\pi^{+}pX$", "green", "^"),
    (h_x,  h_y,  h_yerr,  "HERMES, 0.4 \< z \< 0.5",        "blue",  "s"),
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

ax.set_xlabel(r"$x_{B}$",         fontsize=16)
ax.set_ylabel(r"$F_{LL}/F_{UU}$", fontsize=16)

ax.set_xlim(0.06, 0.7)
ax.set_ylim(0.0, 1.0)

ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
ax.legend(frameon=True, fontsize=14)

plt.tight_layout()
plt.savefig("/u/home/thayward/pip_ALL.pdf")