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
    # "RGA Fa18 Out":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eppi+pi-X/rga_fa18_out_eppi+pi-X.root",
    # "RGA Sp19 Inb":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eppi+pi-X/rga_sp19_inb_eppi+pi-X.root"
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
        ["helicity","beam_pol","DepA","DepW","theta","phi23","W","Q2","z23","Mx23","Mh23"],
        library="np"
    )

    # 5) Apply kinematic cuts
    mask = (
        (data["W"]  > 2) &
        (data["Q2"] > 2) &
        (data["z23"]  > 0.9) &
        (data["Mx23"] < 1.05) &
        (data["Mh23"] > 0.65) &
        (data["Mh23"] < 0.85)
    )
    helicity = data["helicity"][mask]
    beam_pol  = data["beam_pol"][mask]
    DepA      = data["DepA"][mask]
    DepW      = data["DepW"][mask]
    theta     = data["theta"][mask]
    phi       = data["phi23"][mask]

    cos_th = np.cos(theta)

    # Prepare arrays
    amps_unw = np.zeros(len(cos_centers))
    errs_unw = np.zeros(len(cos_centers))
    amps_w   = np.zeros(len(cos_centers))
    errs_w   = np.zeros(len(cos_centers))

    # 6) Loop over cos(theta) bins
    for i in range(len(cos_centers)):
        sel_cos = (cos_th >= cos_bins[i]) & (cos_th < cos_bins[i+1])
        if not np.any(sel_cos):
            amps_unw[i], errs_unw[i] = 0.0, 0.0
            amps_w[i],   errs_w[i]   = 0.0, 0.0
            continue
        #endif

        # Lists for phi-dependent asymmetry
        A_phi_unw, err_phi_unw = [], []
        A_phi_w,   err_phi_w   = [], []

        # 7) Loop over phi bins
        for j in range(len(phi_centers)):
            sel_phi = sel_cos & (phi >= phi_bins[j]) & (phi < phi_bins[j+1])

            # Unweighted asymmetry
            hj = helicity[sel_phi]
            n_p = np.sum(hj > 0)
            n_m = np.sum(hj < 0)
            N   = n_p + n_m
            if N > 0:
                A_j   = (n_p - n_m) / N
                err_j = np.sqrt((1 - A_j**2) / N)
            else:
                A_j, err_j = 0.0, 0.0
            #endif
            A_phi_unw.append(A_j)
            err_phi_unw.append(err_j)

            # Weighted asymmetry
            wj = DepA[sel_phi] / (beam_pol[sel_phi] * DepW[sel_phi])
            w_p   = wj[hj > 0]
            w_m   = wj[hj < 0]
            sp, sm   = np.sum(w_p), np.sum(w_m)
            sp2, sm2 = np.sum(w_p**2), np.sum(w_m**2)
            Nw = sp + sm
            if Nw > 0:
                A_w   = (sp - sm) / Nw
                N_eff = Nw**2 / (sp2 + sm2)
                err_w = np.sqrt((1 - A_w**2) / N_eff)
            else:
                A_w, err_w = 0.0, 0.0
            #endif
            A_phi_w.append(A_w)
            err_phi_w.append(err_w)
        #endfor

        # 8) Fit unweighted A(phi) → A·sin(phi)
        try:
            popt_u, pcov_u     = curve_fit(sin_func, phi_centers, A_phi_unw,
                                            sigma=err_phi_unw, absolute_sigma=True)
            amps_unw[i]        = popt_u[0]
            errs_unw[i]        = np.sqrt(pcov_u[0,0])
        except:
            amps_unw[i], errs_unw[i] = 0.0, 0.0
        #endif

        # 9) Fit weighted A(phi)
        try:
            popt_w, pcov_w   = curve_fit(sin_func, phi_centers, A_phi_w,
                                          sigma=err_phi_w, absolute_sigma=True)
            amps_w[i]        = popt_w[0]
            errs_w[i]        = np.sqrt(pcov_w[0,0])
        except:
            amps_w[i], errs_w[i] = 0.0, 0.0
        #endif

    #endfor

    results[label] = {
        "amp_unw": amps_unw,
        "err_unw": errs_unw,
        "amp_w":   amps_w,
        "err_w":   errs_w
    }
#endfor

# 10) Ensure output directory exists
os.makedirs("output", exist_ok=True)

# 11) First-step (unweighted amplitude) plot
fig, ax = plt.subplots()
for label, res in results.items():
    ax.errorbar(cos_centers, res["amp_unw"], yerr=res["err_unw"],
                fmt='o', label=label)
#endfor
ax.set_xlabel(r'$\cos\theta$')
ax.set_ylabel(r'$A_{LU}$')
ax.set_ylim(-1.2, 0.4)
ax.legend(loc='upper right')
ax.text(0.02, 0.98,
        "W > 2, Q² > 2, z > 0.9, Mₓ < 1.05",
        transform=ax.transAxes, va='top', ha='left')
plt.savefig("output/rga_rho0_first_step.pdf")
plt.close()

# 12) Weighted (final amplitude) plot
fig, ax = plt.subplots()
for label, res in results.items():
    ax.errorbar(cos_centers, res["amp_w"], yerr=res["err_w"],
                fmt='o', label=label)
#endfor
ax.set_xlabel(r'$\cos\theta$')
ax.set_ylabel(r'$F_{LU}/F_{UU}$')
ax.set_ylim(-1.2, 0.4)
ax.legend(loc='upper right')
ax.text(0.02, 0.98,
        "W > 2, Q² > 2, z > 0.9, Mₓ < 1.05",
        transform=ax.transAxes, va='top', ha='left')
plt.savefig("output/rga_rho0_ALU.pdf")
plt.close()