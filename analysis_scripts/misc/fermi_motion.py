#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    ("RGA Fa18 at rest",           rga_epiX_atRest),
    ("RGA Fa18 sim. Fermi Motion", rga_epiX_fermi),
    ("RGC Su22 at rest",           rgc_epiX_atRest),
]

versions_epiPipiX = [
    ("RGA Fa18 at rest",           rga_epiPipi_atRest),
    ("RGA Fa18 sim. Fermi Motion", rga_epiPipi_fermi),
    ("RGC Su22 at rest",           rgc_epiPipi_atRest),
]

# -----------------------------------------------------------------------------
# Load data arrays
# -----------------------------------------------------------------------------
# eπ⁺X → Mx²
epiX_Mx2      = [(lbl, load_array(path, "Mx2"))    for lbl, path in versions_epiX]
# eπ⁺π⁻X → Mx²
epiPipi_Mx2   = [(lbl, load_array(path, "Mx2"))    for lbl, path in versions_epiPipiX]
# eπ⁺π⁻X → Mx²_1 (target‐dependent missing mass for the π⁺)
epiPipi_Mx2_1 = [(lbl, load_array(path, "Mx2_1"))  for lbl, path in versions_epiPipiX]

# -----------------------------------------------------------------------------
# Plot in a 1×3 figure
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Subplot 1: e π⁺ X → Mx²
for lbl, data in epiX_Mx2:
    axes[0].hist(data, bins=100, range=(-1, 3), density=True,
                 histtype="step", label=lbl)
axes[0].set_xlim(-1, 3)
axes[0].set_xlabel(r"$M_x^2\ \mathrm{(GeV^2)}$")
axes[0].set_title(r"$e\,\pi^{+}X:\ M_x^2$")
axes[0].legend()

# Subplot 2: e π⁺π⁻ X → Mx²
for lbl, data in epiPipi_Mx2:
    axes[1].hist(data, bins=100, range=(-1, 3), density=True,
                 histtype="step", label=lbl)
axes[1].set_xlim(-1, 3)
axes[1].set_xlabel(r"$M_x^2\ \mathrm{(GeV^2)}$")
axes[1].set_title(r"$e\,\pi^{+}\pi^{-}X:\ M_x^2$")
axes[1].legend()

# Subplot 3: e π⁺π⁻ X → Mx²_1 (π⁺ missing mass)
for lbl, data in epiPipi_Mx2_1:
    axes[2].hist(data, bins=100, range=(-1, 3), density=True,
                 histtype="step", label=lbl)
axes[2].set_xlim(-1, 3)
axes[2].set_xlabel(r"$M_{x(\pi^{+})}^2\ \mathrm{(GeV^2)}$")
axes[2].set_title(r"$e\,\pi^{+}\pi^{-}X:\ M_{x(\pi^{+})}^2$")
axes[2].legend()

plt.tight_layout()
plt.savefig("output/fermi_motion.pdf")
plt.close()