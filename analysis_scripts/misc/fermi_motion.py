#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Filepaths (hard‐coded; RGC lines commented out until those ROOTs are ready)
# -----------------------------------------------------------------------------
#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Filepaths (hard‐coded; RGC lines commented out until those ROOTs are ready)
# -----------------------------------------------------------------------------
rga_epiX_atRest   = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_atRest.root"
rga_epiX_fermi    = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_fermiMotion.root"
# rgc_epiX_atRest = "/volatile/clas12/thayward/fermi_motion/rgc_fa22_inb_epi+_atRest.root"
rgc_epiX_atRest   = rga_epiX_fermi  # placeholder until RGC is ready

rga_epiPiX_atRest = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+pi-_atRest.root"
rga_epiPiX_fermi  = "/volatile/clas12/thayward/fermi_motion/rga_fa18_inb_epi+_pi-_fermiMotion.root"
# rgc_epiPiX_atRest = "/volatile/clas12/thayward/fermi_motion/rgc_fa22_inb_epi+pi-_atRest.root"
rgc_epiPiX_atRest = rga_epiPiX_fermi  # placeholder

# ensure output directory exists
os.makedirs("output", exist_ok=True)

# -----------------------------------------------------------------------------
# Helper to load a branch array from PhysicsEvents
# -----------------------------------------------------------------------------
def load_array(path, branch):
    with uproot.open(path) as f:
        return f["PhysicsEvents"][branch].array(library="np")

# -----------------------------------------------------------------------------
# 1) epi+X  → Mx2
# -----------------------------------------------------------------------------
versions = [
    ("RGA",      rga_epiX_atRest),
    ("RGA sim. Fermi Motion", rga_epiX_fermi),
    ("RGC",      rgc_epiX_atRest),
]
epiX_Mx2 = [(label, load_array(path, "Mx2")) for label, path in versions]

# -----------------------------------------------------------------------------
# 2) epi+π⁻X → Mx2
# -----------------------------------------------------------------------------
epiPiX_Mx2 = [
    (label, load_array(path, "Mx2"))
    for label, path in [
        ("RGA",      rga_epiPiX_atRest),
        ("RGA sim. Fermi Motion", rga_epiPiX_fermi),
        ("RGC",      rgc_epiPiX_atRest),
    ]
]

# -----------------------------------------------------------------------------
# 3) epi+π⁻X → Mh   (only these trees have Mh)
# -----------------------------------------------------------------------------
epiPiX_Mh = [
    (label, load_array(path, "Mh"))
    for label, path in [
        ("RGA",      rga_epiPiX_atRest),
        ("RGA sim. Fermi Motion", rga_epiPiX_fermi),
        ("RGC",      rgc_epiPiX_atRest),
    ]
]

# -----------------------------------------------------------------------------
# Plot everything in a 1×3 canvas
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# — Plot 1: epi+X Mx² —
for label, data in epiX_Mx2:
    axes[0].hist(data, bins=100, range=(-1, 3), density=True,
                 histtype="step", label=label)
axes[0].set_xlim(-1, 3)
axes[0].set_xlabel(r"$M_{x}^{2}$ (GeV$^{2}$)")
axes[0].set_title("epi+X: $M_x^2$")
axes[0].legend()

# — Plot 2: epi+π⁻X Mx² —
for label, data in epiPiX_Mx2:
    axes[1].hist(data, bins=100, range=(-1, 3), density=True,
                 histtype="step", label=label)
axes[1].set_xlim(-1, 3)
axes[1].set_xlabel(r"$M_{x}^{2}$ (GeV$^{2}$)")
axes[1].set_title("epi+π⁻X: $M_x^2$")
axes[1].legend()

# — Plot 3: epi+π⁻X Mh (0–2 GeV) —
for label, data in epiPiX_Mh:
    axes[2].hist(data, bins=100, range=(0, 2), density=True,
                 histtype="step", label=label)
axes[2].set_xlim(0, 2)
axes[2].set_xlabel(r"$M_{h}$ (GeV)")
axes[2].set_title("epi+π⁻X: $M_h$")
axes[2].legend()

plt.tight_layout()
plt.savefig("output/fermi_motion.pdf")
plt.close()

# -----------------------------------------------------------------------------
# Read arrays from "PhysicsEvents" trees
# -----------------------------------------------------------------------------
def load_array(path, branch):
    """Open ROOT file, grab PhysicsEvents, return numpy array of branch."""
    with uproot.open(path) as f:
        return f["PhysicsEvents"][branch].array(library="np")

# epi+X Mx2 for each version
versions = [
    ("RGA",      rga_epiX_atRest),
    ("RGA sim. Fermi Motion", rga_epiX_fermi),
    ("RGC",      rgc_epiX_atRest),
]
epiX_Mx2 = [(label, load_array(path, "Mx2")) for label, path in versions]

# epi+π⁻X Mx2
epiPiX_Mx2 = [
    (label, load_array(path, "Mx2"))
    for label, path in [
        ("RGA",      rga_epiPiX_atRest),
        ("RGA sim. Fermi Motion", rga_epiPiX_fermi),
        ("RGC",      rgc_epiPiX_atRest),
    ]
]

# epi+X Mh
epiX_Mh = [(label, load_array(path, "Mh")) for label, path in versions]

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 1) epi+X Mx²
for label, data in epiX_Mx2:
    axes[0].hist(data, bins=100, range=(-1, 3), density=True,
                 histtype="step", label=label)
axes[0].set_xlim(-1, 3)
axes[0].set_xlabel(r"$M_{x}^{2}$ (GeV$^{2}$)")
axes[0].set_title("epi+X: $M_x^2$")
axes[0].legend()

# 2) epi+π⁻X Mx²
for label, data in epiPiX_Mx2:
    axes[1].hist(data, bins=100, range=(-1, 3), density=True,
                 histtype="step", label=label)
axes[1].set_xlim(-1, 3)
axes[1].set_xlabel(r"$M_{x}^{2}$ (GeV$^{2}$)")
axes[1].set_title("epi+π⁻X: $M_x^2$")
axes[1].legend()

# 3) epi+X Mh
for label, data in epiX_Mh:
    axes[2].hist(data, bins=100, range=(0, 2), density=True,
                 histtype="step", label=label)
axes[2].set_xlim(0, 2)
axes[2].set_xlabel(r"$M_{h}$ (GeV)")
axes[2].set_title("epi+X: $M_h$")
axes[2].legend()

plt.tight_layout()
plt.savefig("output/fermi_motion.pdf")
plt.close()