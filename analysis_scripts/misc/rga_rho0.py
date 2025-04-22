#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt
import os

#--- 1) Define your files and labels
files = {
    "RGA Fa18 Inb": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+pi-X/rga_fa18_inb_epi+pi-X.root",
    "RGA Fa18 Out": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+pi-X/rga_fa18_out_epi+pi+X.root",
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
    data = tree.arrays(
        ["helicity","beam_pol","DepA","DepW","theta","W","Q2","z","Mx"],
        library="np"
    )

    # Apply kinematic cuts
    mask = (
        (data["W"]  > 2) &
        (data["Q2"] > 2) &
        (data["z"]  > 0.9) &
        (data["Mx"] < 1.05)
    )

    helicity = data["helicity"][mask]
    beam_pol  = data["beam_pol"][mask]
    DepA      = data["DepA"][mask]
    DepW      = data["DepW"][mask]
    theta     = data["theta"][mask]

    cos_theta = np.cos(theta)

    # Prepare arrays
    A_first       = np.zeros(len(bin_centers))
    err_first     = np.zeros(len(bin_centers))
    A_weighted    = np.zeros(len(bin_centers))
    err_weighted  = np.zeros(len(bin_centers))
    mean_beam_pol = np.zeros(len(bin_centers))
    mean_DepA     = np.zeros(len(bin_centers))
    mean_DepW     = np.zeros(len(bin_centers))

    #--- 4) Loop over cos(theta) bins
    for i in range(len(bin_centers)):
        sel = (cos_theta >= bins[i]) & (cos_theta < bins[i+1])

        # Unweighted counts
        h = helicity[sel]
        n_plus  = np.sum(h > 0)
        n_minus = np.sum(h < 0)
        Ntot    = n_plus + n_minus

        if Ntot > 0:
            A_first[i]   = (n_plus - n_minus) / Ntot
            err_first[i] = np.sqrt((1 - A_first[i]**2) / Ntot)
        else:
            A_first[i], err_first[i] = 0.0, 0.0
        #endif

        # Weighted asymmetry
        bp = beam_pol[sel]
        dA = DepA[sel]
        dW = DepW[sel]
        w  = dA / (bp * dW)

        w_p   = w[h > 0]
        w_m   = w[h < 0]
        sum_p = np.sum(w_p)
        sum_m = np.sum(w_m)
        sum_p2 = np.sum(w_p**2)
        sum_m2 = np.sum(w_m**2)
        Ntot_w = sum_p + sum_m

        if Ntot_w > 0:
            A_weighted[i]   = (sum_p - sum_m) / Ntot_w
            N_eff = Ntot_w**2 / (sum_p2 + sum_m2)
            err_weighted[i] = np.sqrt((1 - A_weighted[i]**2) / N_eff)
        else:
            A_weighted[i], err_weighted[i] = 0.0, 0.0
        #endif

        # Track means
        if np.any(sel):
            mean_beam_pol[i] = np.mean(bp)
            mean_DepA[i]     = np.mean(dA)
            mean_DepW[i]     = np.mean(dW)
        else:
            mean_beam_pol[i] = np.nan
            mean_DepA[i]     = np.nan
            mean_DepW[i]     = np.nan
        #endif
    #endfor

    results[label] = {
        "A_first":       A_first,
        "err_first":     err_first,
        "A_weighted":    A_weighted,
        "err_weighted":  err_weighted,
        "mean_beam_pol": mean_beam_pol,
        "mean_DepA":     mean_DepA,
        "mean_DepW":     mean_DepW
    }
#endfor

# Ensure output directory exists
os.makedirs("output", exist_ok=True)

#--- 5) First-step (unweighted) plot
fig, ax = plt.subplots()
for label, res in results.items():
    ax.errorbar(bin_centers, res["A_first"], yerr=res["err_first"], fmt='o', label=label)
#endfor
ax.set_xlabel(r'$\cos\theta$')
ax.set_ylabel(r'$A_{LU}$')
ax.set_ylim(-1.2, 0.4)
ax.legend(loc='upper right')
ax.text(
    0.02, 0.98,
    "W > 2, Q² > 2, z > 0.9, Mₓ < 1.05",
    transform=ax.transAxes, va='top', ha='left'
)
plt.savefig("output/rga_rho0_first_step.pdf")
plt.close()

#--- 6) Weighted (final) plot
fig, ax = plt.subplots()
for label, res in results.items():
    ax.errorbar(bin_centers, res["A_weighted"], yerr=res["err_weighted"], fmt='o', label=label)
#endfor
ax.set_xlabel(r'$\cos\theta$')
ax.set_ylabel(r'$F_{LU}/F_{UU}$')
ax.set_ylim(-1.2, 0.4)
ax.legend(loc='upper right')
ax.text(
    0.02, 0.98,
    "W > 2, Q² > 2, z > 0.9, Mₓ < 1.05",
    transform=ax.transAxes, va='top', ha='left'
)
plt.savefig("output/rga_rho0_ALU.pdf")
plt.close()