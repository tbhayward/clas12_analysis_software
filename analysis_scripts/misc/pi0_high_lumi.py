#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt

input_file = "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root"
output_file = "/u/home/thayward/pi0_variables.png"

tree_name = "PhysicsEvents"

# ----------------------------------------------------------------------
# Read variables
# ----------------------------------------------------------------------

with uproot.open(input_file) as f:
    tree = f[tree_name]
    arrays = tree.arrays(["pT2", "t1", "z2", "Mx2"], library="np")

pT = arrays["pT2"]
t  = np.abs(arrays["t1"])
z  = arrays["z2"]
Mx2 = arrays["Mx2"]

# Remove non-finite entries.
mask = (
    np.isfinite(pT) &
    np.isfinite(t) &
    np.isfinite(z) &
    np.isfinite(Mx2)
)

pT  = pT[mask]
t   = t[mask]
z   = z[mask]
Mx2 = Mx2[mask]

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# pT vs |t|
h = axes[0, 0].hist2d(
    pT,
    t,
    bins=(100, 100),
    range=((0.0, 1.4), (0.0, 1.2)),
    cmap="coolwarm"
)

axes[0, 0].set_xlabel(r"$p_T$ (GeV)")
axes[0, 0].set_ylabel(r"$|t|$ (GeV$^2$)")
axes[0, 0].set_title(r"$p_T$ vs. $|t|$")
fig.colorbar(h[3], ax=axes[0, 0], label="Counts")

# pT vs z
h = axes[0, 1].hist2d(
    pT,
    z,
    bins=(100, 100),
    range=((0.0, 1.4), (0.4, 1.2)),
    cmap="coolwarm"
)

axes[0, 1].set_xlabel(r"$p_T$ (GeV)")
axes[0, 1].set_ylabel(r"$z$")
axes[0, 1].set_title(r"$p_T$ vs. $z$")
fig.colorbar(h[3], ax=axes[0, 1], label="Counts")

# |t| vs z
h = axes[1, 0].hist2d(
    t,
    z,
    bins=(100, 100),
    range=((0.0, 1.2), (0.4, 1.2)),
    cmap="coolwarm"
)

axes[1, 0].set_xlabel(r"$|t|$ (GeV$^2$)")
axes[1, 0].set_ylabel(r"$z$")
axes[1, 0].set_title(r"$|t|$ vs. $z$")
fig.colorbar(h[3], ax=axes[1, 0], label="Counts")

# Mx2 distribution
axes[1, 1].hist(
    Mx2,
    bins=150,
    range=(-0.25, 0.25),
    histtype="step",
    linewidth=1.5
)

axes[1, 1].set_xlim(-0.25, 0.25)
axes[1, 1].set_xlabel(r"$M_X^2$ (GeV$^2$)")
axes[1, 1].set_ylabel("Counts")
axes[1, 1].set_title(r"$M_X^2$ Distribution")

# ----------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------

plt.tight_layout()
plt.savefig(output_file, dpi=200, bbox_inches="tight")
plt.close()

print(f"Saved: {output_file}")