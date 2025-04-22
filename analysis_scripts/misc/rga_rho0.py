#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt
import os

#--- 1) Define your files and labels
files = {
    "RGA Fa18 Inb": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+pi-X/rga_fa18_inb_epi+pi-X.root",
    "RGA Fa18 Out": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+pi-X/rga_fa18_out_epi+pi-X.root",
    "RGA Sp19 Inb": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+pi-X/rga_sp19_inb_epi+pi-X.root"
}

#--- 2) Set up bins in cos(theta)
bins = np.linspace(-0.9, 0.9, 16)
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Storage for results
results = {}

#--- 3) Loop over each dataset
for label, path in files.items():
    tree = uproot.open(path)["PhysicsEvents"]
    # Read only the branches you need
    data = tree.arrays(["helicity","beam_pol","DepA","DepW","theta","W","Q2","z","Mx"], library="np")

    # Apply kinematic cuts
    mask_cut = (
        (data["W"]  > 2) &
        (data["Q2"] > 2) &
        (data["z"]  > 0.9) &
        (data["Mx"] < 1.05)
    )

    helicity = data["helicity"][mask_cut]
    beam_pol  = data["beam_pol"][mask_cut]
    DepA      = data["DepA"][mask_cut]
    DepW      = data["DepW"][mask_cut]
    theta     = data["theta"][mask_cut]

    # Compute per‑event weight = DepA/(beam_pol*DepW)
    weights = DepA / (beam_pol * DepW)
    print("weights:", weights)

    # cos(theta)
    cos_theta = np.cos(theta)

    # Arrays for asymmetry and error
    A_vals   = np.zeros(len(bin_centers))
    err_vals = np.zeros(len(bin_centers))

    #--- 4) Loop over cos(theta) bins
    for i in range(len(bin_centers)):
        sel = (cos_theta >= bins[i]) & (cos_theta < bins[i+1])
        h = helicity[sel]
        w = weights[sel]

        # Split by helicity sign
        w_p = w[h > 0]
        w_m = w[h < 0]

        sum_p  = np.sum(w_p)
        sum_m  = np.sum(w_m)
        sum_p2 = np.sum(w_p**2)
        sum_m2 = np.sum(w_m**2)

        Ntot = sum_p + sum_m

        if Ntot > 0:
            A    = (sum_p - sum_m) / Ntot
            # effective number of events
            N_eff = Ntot**2 / (sum_p2 + sum_m2)
            err  = np.sqrt((1 - A**2) / N_eff)
        else:
            A, err = 0.0, 0.0
        # endif

        A_vals[i]   = A
        err_vals[i] = err
    # endfor

    results[label] = (A_vals, err_vals)
# endfor

#--- 5) Plotting
fig, ax = plt.subplots()
for label, (A_vals, err_vals) in results.items():
    ax.errorbar(bin_centers, A_vals, yerr=err_vals, fmt='o', label=label)
# endfor

ax.set_xlabel(r'$\cos\theta$')
ax.set_ylabel(r'$F_{LU}/F_{UU}$')
ax.legend(loc='upper right')
ax.text(
    0.02, 0.98,
    "W > 2, Q² > 2, z > 0.9, Mₓ < 1.05",
    transform=ax.transAxes, va='top', ha='left'
)

# Ensure output directory exists
os.makedirs(os.path.dirname("output/rga_rho0_ALU.pdf"), exist_ok=True)
plt.savefig("output/rga_rho0_ALU.pdf")
plt.close()