#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt
import os

# 1) File definitions
rec_file = "/volatile/clas12/thayward/rho0_study/rga_fa18_inb_epi+pi-X_rec.root"
gen_file = "/volatile/clas12/thayward/rho0_study/rga_fa18_inb_epi+pi-X_gen.root"

# 2) Open Trees and read branches
rec_tree = uproot.open(rec_file)["PhysicsEvents"]
gen_tree = uproot.open(gen_file)["PhysicsEvents"]

rec = rec_tree.arrays(["theta","phi","W","Q2","z","Mx"], library="np")
gen = gen_tree.arrays(["theta","phi","W","Q2","z","Mx"], library="np")

# 3) Apply kinematic cuts
mask_rec = (
    (rec["W"]  > 2) &
    (rec["Q2"] > 2) &
    (rec["z"]  > 0.9) &
    (rec["Mx"] < 1.05)
)
mask_gen = (
    (gen["W"]  > 2) &
    (gen["Q2"] > 2) &
    (gen["z"]  > 0.9) &
    (gen["Mx"] < 1.05)
)

theta_rec = rec["theta"][mask_rec]
theta_gen = gen["theta"][mask_gen]
phi_rec   = rec["phi"][mask_rec]
phi_gen   = gen["phi"][mask_gen]

# 4) Define bins
# cos(theta) from -1 to 1, 30 bins
cos_bins    = np.linspace(-1, 1, 31)
cos_centers = 0.5 * (cos_bins[:-1] + cos_bins[1:])
# phi from 0 to 2π, 30 bins
phi_bins    = np.linspace(0, 2*np.pi, 31)
phi_centers = 0.5 * (phi_bins[:-1] + phi_bins[1:])

# 5) Histogram counts
cos_rec_vals = np.cos(theta_rec)
cos_gen_vals = np.cos(theta_gen)

counts_cos_rec, _ = np.histogram(cos_rec_vals, bins=cos_bins)
counts_cos_gen, _ = np.histogram(cos_gen_vals, bins=cos_bins)

counts_phi_rec, _ = np.histogram(phi_rec, bins=phi_bins)
counts_phi_gen, _ = np.histogram(phi_gen, bins=phi_bins)

# 6) Compute efficiencies and binomial errors
eff_cos = np.zeros_like(counts_cos_rec, dtype=float)
err_cos = np.zeros_like(counts_cos_rec, dtype=float)
mask_cos_nonzero = counts_cos_gen > 0
eff_cos[mask_cos_nonzero] = counts_cos_rec[mask_cos_nonzero] / counts_cos_gen[mask_cos_nonzero]
err_cos[mask_cos_nonzero] = np.sqrt(eff_cos[mask_cos_nonzero] * 
                                    (1 - eff_cos[mask_cos_nonzero]) / 
                                    counts_cos_gen[mask_cos_nonzero])

eff_phi = np.zeros_like(counts_phi_rec, dtype=float)
err_phi = np.zeros_like(counts_phi_rec, dtype=float)
mask_phi_nonzero = counts_phi_gen > 0
eff_phi[mask_phi_nonzero] = counts_phi_rec[mask_phi_nonzero] / counts_phi_gen[mask_phi_nonzero]
err_phi[mask_phi_nonzero] = np.sqrt(eff_phi[mask_phi_nonzero] * 
                                    (1 - eff_phi[mask_phi_nonzero]) / 
                                    counts_phi_gen[mask_phi_nonzero])

# 7) Plot side-by-side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

# -- cos(theta) efficiency
ax1.errorbar(cos_centers, eff_cos, yerr=err_cos, fmt='o', linestyle='none')
ax1.set_yscale('log')
ax1.set_ylim(1e-5, 1e-1)
ax1.set_xlabel(r'$\cos\theta$')
ax1.set_ylabel("Rec. Eff.")
ax1.grid(True, which='both', ls='--', lw=0.5)
ax1.set_title("Efficiency vs. cosθ")

# -- phi efficiency
ax2.errorbar(phi_centers, eff_phi, yerr=err_phi, fmt='o', linestyle='none')
ax2.set_yscale('log')
ax2.set_ylim(1e-5, 1e-1)
ax2.set_xlabel(r'$\phi$ (rad)')
ax2.grid(True, which='both', ls='--', lw=0.5)
ax2.set_title("Efficiency vs. φ")

plt.tight_layout()

# 8) Save output
os.makedirs("output", exist_ok=True)
plt.savefig("output/rho0_theta_acceptance.pdf")
plt.close()