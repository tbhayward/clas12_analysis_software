#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# File paths for RGA datasets (two files per channel)
# -----------------------------------------------------------------------------
rga_epiX_atRest   = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_atRest.root"
rga_epiX_fermi    = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_fermiMotion.root"

rga_epiPiX_atRest = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+pi-_atRest.root"
rga_epiPiX_fermi  = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_pi-_fermiMotion.root"

# Uncomment and set real RGC paths when they’re ready:
# rgc_epiX_atRest   = "/volatile/clas12/thayward/fermi_motion/rgc_fa22_inb_epi+_atRest.root"
# rgc_epiPiX_atRest = "/volatile/clas12/thayward/fermi_motion/rgc_fa22_inb_epi+pi-_atRest.root"

# Ensure output directory exists
os.makedirs("output", exist_ok=True)

def load_array(path, branch):
    """Load a branch array from the PhysicsEvents tree."""
    with uproot.open(path) as f:
        return f["PhysicsEvents"][branch].array(library="np")

# Define the versions for each channel (RGA only for now)
versions_epiX = [
    ("RGA at rest",           rga_epiX_atRest),
    ("RGA sim. Fermi Motion", rga_epiX_fermi),
    # ("RGC", rgc_epiX_atRest),
]

versions_epiPiX = [
    ("RGA at rest",           rga_epiPiX_atRest),
    ("RGA sim. Fermi Motion", rga_epiPiX_fermi),
    # ("RGC", rgc_epiPiX_atRest),
]

# Load Mx2 and Mh arrays
epiX_Mx2   = [(lbl, load_array(path, "Mx2")) for lbl, path in versions_epiX]
epiPiX_Mx2 = [(lbl, load_array(path, "Mx2")) for lbl, path in versions_epiPiX]
epiPiX_Mh  = [(lbl, load_array(path, "Mh"))   for lbl, path in versions_epiPiX]

# Create 1×3 subplot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Subplot 1: epi+X Mx²
for lbl, data in epiX_Mx2:
    axes[0].hist(data,
                 bins=100,
                 range=(-1, 3),
                 density=True,
                 histtype="step",
                 label=lbl)
axes[0].set_xlim(-1, 3)
axes[0].set_xlabel(r"$M_x^2$ (GeV$^2$)")
axes[0].set_title("epi+X: $M_x^2$")
axes[0].legend()

# Subplot 2: epi+π⁻X Mx²
for lbl, data in epiPiX_Mx2:
    axes[1].hist(data,
                 bins=100,
                 range=(-1, 3),
                 density=True,
                 histtype="step",
                 label=lbl)
axes[1].set_xlim(-1, 3)
axes[1].set_xlabel(r"$M_x^2$ (GeV$^2$)")
axes[1].set_title("epi+π⁻X: $M_x^2$")
axes[1].legend()

# Subplot 3: epi+π⁻X Mh
for lbl, data in epiPiX_Mh:
    axes[2].hist(data,
                 bins=100,
                 range=(0, 2),
                 density=True,
                 histtype="step",
                 label=lbl)
axes[2].set_xlim(0, 2)
axes[2].set_xlabel(r"$M_h$ (GeV)")
axes[2].set_title("epi+π⁻X: $M_h$")
axes[2].legend()

plt.tight_layout()
plt.savefig("output/fermi_motion.pdf")
plt.close()