#!/usr/bin/env python3

import os

import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


INPUT_FILE = (
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/"
    "rga_fa18_out_epgamma.root"
)

OUTPUT_FILE = "output/reach_plot.png"


# ------------------------------------------------------------
# Read ROOT tree
# ------------------------------------------------------------

with uproot.open(INPUT_FILE) as f:

    # Find the first TTree in the file.
    tree_name = None

    for key, obj in f.items():
        if isinstance(obj, uproot.behaviors.TTree.TTree):
            tree_name = key
            break
        #endif
    #endfor

    if tree_name is None:
        raise RuntimeError(f"No TTree found in {INPUT_FILE}")
    #endif

    print(f"Using tree: {tree_name}")

    tree = f[tree_name]

    arrays = tree.arrays(
        ["Q2", "t1"],
        library="np",
    )

#endwith


# ------------------------------------------------------------
# Extract variables
# ------------------------------------------------------------

Q2 = arrays["Q2"]
t1 = arrays["t1"]

# Mandelstam t is negative, so plot |t| = -t.
t_abs = -t1


# ------------------------------------------------------------
# Remove invalid entries
# ------------------------------------------------------------

mask = (
    np.isfinite(Q2)
    & np.isfinite(t_abs)
    & (Q2 > 0.0)
    & (t_abs > 0.0)
)

Q2 = Q2[mask]
t_abs = t_abs[mask]

print(f"Events plotted: {len(Q2):,}")


# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

os.makedirs(
    os.path.dirname(OUTPUT_FILE),
    exist_ok=True,
)


# ------------------------------------------------------------
# Make 2D histogram
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 7))

h = ax.hist2d(
    t_abs,
    Q2,
    bins=[100, 100],
    range=[
        [0.0, 1.2],
        [1.0, 8.0],
    ],
    norm=LogNorm(),
    cmap="viridis",
)


# ------------------------------------------------------------
# Color bar
# ------------------------------------------------------------

cbar = fig.colorbar(
    h[3],
    ax=ax,
)

cbar.set_label(
    "Counts",
    fontsize=16,
)

cbar.ax.tick_params(
    labelsize=13,
)


# ------------------------------------------------------------
# GPD factorization boundary
#
#       |t| / Q^2 = 0.2
#
# corresponds to
#
#       Q^2 = 5 |t|
#
# Events above this line satisfy |t|/Q^2 < 0.2.
# ------------------------------------------------------------

t_line = np.linspace(
    0.0,
    1.2,
    500,
)

Q2_line = 5.0 * t_line

ax.plot(
    t_line,
    Q2_line,
    linewidth=2.5,
    label=r"$|t|/Q^2 = 0.2$",
)


# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------

ax.set_xlim(
    0.0,
    1.2,
)

ax.set_ylim(
    1.0,
    8.0,
)

ax.set_xlabel(
    r"$|t|$ (GeV$^2$)",
    fontsize=18,
)

ax.set_ylabel(
    r"$Q^2$ (GeV$^2$)",
    fontsize=18,
)

ax.tick_params(
    axis="both",
    labelsize=14,
)

ax.legend(
    loc="upper left",
    fontsize=14,
    frameon=True,
)

ax.set_title(
    "CLAS12 RGA DVCS Kinematic Reach",
    fontsize=18,
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fig.tight_layout()

fig.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight",
)

print(f"Saved: {OUTPUT_FILE}")

plt.close(fig)