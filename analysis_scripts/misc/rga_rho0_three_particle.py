#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# sinusoidal fit function
def sin_func(phi, A):
    return A * np.sin(phi)

# 1) Define your files and labels
files = {
    "RGA Fa18 Inb":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eppi+pi-X/rga_fa18_inb_eppi+pi-X.root",
    "RGA Fa18 Out":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eppi+pi-X/rga_fa18_out_eppi+pi-X.root",
    "RGA Sp19 Inb":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eppi+pi-X/rga_sp19_inb_eppi+pi-X.root"
}

# 2) Set up bins in cos(theta)
cos_bins    = np.linspace(-0.9, 0.9, 16)
cos_centers = 0.5 * (cos_bins[:-1] + cos_bins[1:])
# 3) Phi bins for fit
phi_bins    = np.linspace(0, np.pi, 13)
phi_centers = 0.5 * (phi_bins[:-1] + phi_bins[1:])

# Storage for amplitudes and uncertainties
results = {}

# 4) Loop over each dataset
for label, path in files.items():
    tree = uproot.open(path)["PhysicsEvents"]
    data = tree.arrays(
        ["helicity","beam_pol","DepA","DepW","theta","phi23","W","Q2","z23","Mx2_23","Mh23"],
        library="np"
    )

    # 5) Apply kinematic cuts (uncomment & adjust as needed)
    mask = (
        (data["W"]    > 2) &
        (data["Q2"]   > 2) &
        (data["z23"]  > 0.85) &
        (data["Mx2_23"]<1.05) &
        (data["Mh23"] > 0.65) &
        (data["Mh23"] < 0.85) 
        np.ones(len(data["helicity"]), dtype=bool)
    )

    # print("HELLO WORLD")
    helicity = data["helicity"][mask]
    beam_pol  = data["beam_pol"][mask]
    DepA      = data["DepA"][mask]
    DepW      = data["DepW"][mask]
    theta     = data["theta"][mask]
    phi       = data["phi23"][mask]

    cos_th = np.cos(theta)

    # Prepare arrays
    amps_unw  = np.zeros(len(cos_centers))
    errs_unw  = np.zeros(len(cos_centers))
    mean_bp   = np.zeros(len(cos_centers))
    mean_DepA = np.zeros(len(cos_centers))
    mean_DepW = np.zeros(len(cos_centers))

    # 6) Loop over cos(theta) bins
    for i in range(len(cos_centers)):
        sel_cos = (cos_th >= cos_bins[i]) & (cos_th < cos_bins[i+1])
        if not np.any(sel_cos):
            amps_unw[i], errs_unw[i] = 0.0, 0.0
            mean_bp[i]   = np.nan
            mean_DepA[i] = np.nan
            mean_DepW[i] = np.nan
        #endif
        else:
            # track means in this cosθ bin
            mean_bp[i]   = np.mean(beam_pol[sel_cos])
            mean_DepA[i] = np.mean(DepA[sel_cos])
            mean_DepW[i] = np.mean(DepW[sel_cos])

            # lists for phi-dependent fit
            phi_vals = []
            A_phi    = []
            err_phi  = []

            # 7) Loop over phi bins
            for j in range(len(phi_centers)):
                sel_phi = sel_cos & (phi >= phi_bins[j]) & (phi < phi_bins[j+1])
                N       = np.count_nonzero(sel_phi)
                if N < 4:
                    continue
                #endif

                hj = helicity[sel_phi]
                n_p = np.sum(hj > 0)
                n_m = np.sum(hj < 0)
                A_j   = (n_p - n_m) / N
                err_j = np.sqrt((1 - A_j**2) / N)

                phi_vals.append(phi_centers[j])
                A_phi.append(A_j)
                err_phi.append(err_j)
            #endfor

            # 8) Fit unweighted A(phi) → A·sin(phi)
            if len(A_phi) >= 4:
                try:
                    popt, pcov     = curve_fit(
                        sin_func, phi_vals, A_phi,
                        sigma=err_phi, absolute_sigma=True
                    )
                    amps_unw[i]    = popt[0]
                    errs_unw[i]    = np.sqrt(pcov[0,0])
                except:
                    amps_unw[i], errs_unw[i] = 0.0, 0.0
                #endif
            else:
                amps_unw[i], errs_unw[i] = 0.0, 0.0
            #endif
        #endif
    #endfor

    # 9) Apply DepA/DepW & beam_pol scaling
    amps_corr = amps_unw * mean_DepA / (mean_bp * mean_DepW)
    errs_corr = errs_unw * mean_DepA / (mean_bp * mean_DepW)

    results[label] = {
        "amp_unw": amps_unw,
        "err_unw": errs_unw,
        "amp_w":   amps_corr,
        "err_w":   errs_corr
    }
#endfor

# 10) Ensure output directory exists
os.makedirs("output", exist_ok=True)

# 11) First-step (unweighted amplitude) plot
fig, ax = plt.subplots()
for label, res in results.items():
    ax.errorbar(cos_centers, res["amp_unw"], yerr=res["err_unw"], fmt='o', label=label)
#endfor
ax.set_xlabel(r'$\cos\theta$')
ax.set_ylabel(r'$A_{LU}$')
ax.set_ylim(-1.2, 0.4)
ax.legend(loc='upper right')
ax.text(
    0.02, 0.98,
    "W > 2, Q² > 2, z > 0.85, Mₓ < 1.05",
    transform=ax.transAxes, va='top', ha='left'
)
plt.savefig("output/rga_rho0_first_step.pdf")
plt.close()

# 12) Corrected (final) amplitude plot
fig, ax = plt.subplots()
for label, res in results.items():
    ax.errorbar(cos_centers, res["amp_w"], yerr=res["err_w"], fmt='o', label=label)
#endfor
ax.set_xlabel(r'$\cos\theta$')
ax.set_ylabel(r'$F_{LU}/F_{UU}$')
ax.set_ylim(-1.2, 0.4)
ax.legend(loc='upper right')
ax.text(
    0.02, 0.98,
    "W > 2, Q² > 2, z > 0.85, Mₓ < 1.05",
    transform=ax.transAxes, va='top', ha='left'
)
plt.savefig("output/rga_rho0_ALU.pdf")
plt.close()