#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt
import os

# 1) File definitions
rec_file = "/volatile/clas12/thayward/rho0_study/rga_fa18_inb_epi+pi-X_rec.root"
gen_file = "/volatile/clas12/thayward/rho0_study/rga_fa18_inb_epi+pi-X_gen.root"

# 2) Open trees and read branches
rec_tree = uproot.open(rec_file)["PhysicsEvents"]
gen_tree = uproot.open(gen_file)["PhysicsEvents"]

rec = rec_tree.arrays(["theta", "W", "Q2", "z", "Mx2"], library="np")
gen = gen_tree.arrays(["theta", "W", "Q2", "z", "Mx2"], library="np")

# 3) Apply kinematic cuts
mask_rec = (
    (rec["W"]  > 2) &
    (rec["Q2"] > 2) &
    (rec["z"]  > 0.9) &
    (rec["Mx2"] < 1.11025)
)
mask_gen = (
    (gen["W"]  > 2) &
    (gen["Q2"] > 2) &
    (gen["z"]  > 0.9) &
    (gen["Mx2"] < 1.11025)
)

cos_rec = np.cos(rec["theta"][mask_rec])
cos_gen = np.cos(gen["theta"][mask_gen])

# 4) Binning in cos(theta)
bins    = np.linspace(-1, 1, 31)
centers = 0.5 * (bins[:-1] + bins[1:])

# 5) Histogram counts
counts_rec, _ = np.histogram(cos_rec, bins=bins)
counts_gen, _ = np.histogram(cos_gen, bins=bins)

# 6) Efficiency calculation (avoid divide-by-zero)
eff = np.zeros_like(counts_rec, dtype=float)
nonzero = counts_gen > 0
eff[nonzero] = counts_rec[nonzero] / counts_gen[nonzero]

# 7) Plotting
fig, ax = plt.subplots()
ax.plot(centers, eff, 'o-')
ax.set_yscale('log')
ax.set_ylim(10e-5, 10e-1)          # from 10e-5 to 10e-1
ax.set_xlabel(r'$\cos\theta$')
ax.set_ylabel("Rec. Eff.")
ax.grid(True, which='both', ls='--', lw=0.5)

# 8) Save output
os.makedirs("output", exist_ok=True)
plt.savefig("output/rho0_theta_acceptance.pdf")
plt.close()