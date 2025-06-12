#!/usr/bin/env python3
"""
Script: invariant_mass_pi-_check.py
Compute and plot the invariant mass of the scattered electron and #pi^{+}
for ep -> enpi+ candidates, in parallel for RGC Su22, RGC Fa22, and RGC Sp23,
with progress updates and single-pass histogram filling.
"""
import os
import sys
import time
import numpy as np
import uproot
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

# Physical masses (GeV)
m_e  = 0.000511
m_pi = 0.13957

#===============================================================================
def beam_energy(runnum):
    """Return beam energy (GeV) based on run number."""
    if 6616  <= runnum <= 6783:
        return 10.1998
    elif 16042 <= runnum <= 17065:
        return 10.5473
    elif 17067 <= runnum <= 17724:
        return 10.5563
    elif 17725 <= runnum <= 17811:
        return 10.5593
    else:
        return 0.0

#===============================================================================
def compute_t(runnum, e_p, e_th, e_ph, p_p, p_th, p_ph):
    """Compute Mandelstam t = (q - p_pi)^2."""
    Eb = beam_energy(runnum)
    if Eb <= 0.0:
        return np.nan
    E_e  = np.sqrt(e_p**2 + m_e**2)
    sin_e = np.sin(e_th); cos_e = np.cos(e_th)
    ex = e_p * sin_e * np.cos(e_ph)
    ey = e_p * sin_e * np.sin(e_ph)
    ez = e_p * cos_e

    E_pi  = np.sqrt(p_p**2 + m_pi**2)
    sin_p = np.sin(p_th); cos_p = np.cos(p_th)
    px = p_p * sin_p * np.cos(p_ph)
    py = p_p * sin_p * np.sin(p_ph)
    pz = p_p * cos_p

    E_q = Eb - E_e
    qx, qy, qz = -ex, -ey, Eb - ez

    dE = E_q - E_pi
    dx = qx - px; dy = qy - py; dz = qz - pz
    return dE*dE - (dx*dx + dy*dy + dz*dz)

#===============================================================================
def analyze_period(period_name, file_path, branches, edges):
    """Analyze one run period and save the invariant mass canvas with progress."""
    start_time = time.time()
    print(f"Starting {period_name}...")

    # Load data
    with uproot.open(file_path) as f:
        if "PhysicsEvents" not in f:
            raise KeyError(f"Tree 'PhysicsEvents' not found in {file_path}")
        data = f["PhysicsEvents"].arrays(branches, library="np")

    # Extract arrays
    runnum = data["runnum"]
    e_p    = data["e_p"];    e_th = data["e_theta"]; e_ph = data["e_phi"]
    p_p_   = data["p_p"];    p_th = data["p_theta"]; p_ph = data["p_phi"]
    x      = data["x"]
    y      = data["y"]
    z      = data["z"]
    Q2     = data["Q2"]
    Mx2    = data["Mx2"]
    n_events = len(runnum)
    print(f"Loaded {n_events} events.")

    # Compute t-values
    print("Computing t-values...")
    compute_t_vec = np.vectorize(compute_t)
    t_vals = compute_t_vec(runnum, e_p, e_th, e_ph, p_p_, p_th, p_ph)

    # Apply kinematic cuts
    mask = (
        (np.abs(t_vals) > 0.07) & (np.abs(t_vals) < 0.7) &
        (y < 0.65) & (z > 0.55) &
        (Q2 < 8.0) & (Mx2 > 0.75) & (Mx2 < 1.05)
    )
    n_sel = np.count_nonzero(mask)
    print(f"Events after cuts: {n_sel} ({n_sel/n_events*100:.1f}%)")

    # Pre-calculate invariant mass and x for all selected events (single pass)
    e_p_m  = e_p[mask];    e_th_m = e_th[mask];    e_ph_m = e_ph[mask]
    p_p_m  = p_p_[mask];   p_th_m = p_th[mask];    p_ph_m = p_ph[mask]
    x_m    = x[mask]

    # Electron 4-vector
    E_e_m = np.sqrt(e_p_m**2 + m_e**2)
    sin_e_m = np.sin(e_th_m); cos_e_m = np.cos(e_th_m)
    ex_m = e_p_m * sin_e_m * np.cos(e_ph_m)
    ey_m = e_p_m * sin_e_m * np.sin(e_ph_m)
    ez_m = e_p_m * cos_e_m

    # Pion 4-vector
    E_pi_m = np.sqrt(p_p_m**2 + m_pi**2)
    sin_p_m = np.sin(p_th_m); cos_p_m = np.cos(p_th_m)
    px_m = p_p_m * sin_p_m * np.cos(p_ph_m)
    py_m = p_p_m * sin_p_m * np.sin(p_ph_m)
    pz_m = p_p_m * cos_p_m

    # Compute invariant mass
    E_tot   = E_e_m + E_pi_m
    px_tot  = ex_m + px_m
    py_tot  = ey_m + py_m
    pz_tot  = ez_m + pz_m
    inv_mass = np.sqrt(E_tot**2 - (px_tot**2 + py_tot**2 + pz_tot**2))

    # Assign events to xB bins
    bin_idx = np.digitize(x_m, edges) - 1  # indices 0..n_bins-1
    n_bins = len(edges) - 1

    # Prepare canvas and hist settings
    fig, axes = plt.subplots(3, 4, figsize=(16, 12), sharex=True, sharey=True)
    axes_flat = axes.flatten()

    # Fill and plot histograms per bin
    for i, ax in enumerate(axes_flat):  # loop pads #endfor
        if i < n_bins:
            # select events for this bin
            sel_mask = (bin_idx == i)
            count_i = np.count_nonzero(sel_mask)
            print(f"  Bin {i+1}/{n_bins} ({edges[i]:.2f}<x<{edges[i+1]:.2f}): {count_i} events")
            masses_i = inv_mass[sel_mask]
            ax.hist(masses_i, bins=50, range=(0, 1.5))
            ax.set_title(f"{edges[i]:.2f} < x < {edges[i+1]:.2f}")
            ax.set_xlim(0, 1.5)
        else:
            ax.axis('off')
        ax.set_xlabel(r"$M_{e\pi}$ (GeV)")
        ax.set_ylabel("Counts")

    plt.suptitle(
        fr"{period_name}: $ep \rightarrow en\pi^{{+}}$,"
        r" $0.07 < |t| < 0.7$, $z > 0.55$, $y < 0.65$,"
        r" $0.75 < M_x^2 < 1.05\ (\mathrm{{GeV}}^2)$",
        fontsize=16, y=0.96
    )
    plt.tight_layout(rect=[0,0,1,0.94])

    # Save output
    out_dir = "output/enpi+"
    os.makedirs(out_dir, exist_ok=True)
    safe = period_name.lower().replace(' ', '_')
    out_path = os.path.join(out_dir, f"invariant_mass_{safe}.pdf")
    fig.savefig(out_path)
    plt.close(fig)

    elapsed = time.time() - start_time
    print(f"Completed {period_name} in {elapsed:.1f}s, saved to {out_path}")

#===============================================================================
def main():
    periods = [
        ("RGC Su22", "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_NH3_epi+.root"),
        ("RGC Fa22", "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+.root"),
        ("RGC Sp23", "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_NH3_epi+.root"),
    ]
    branches = ["runnum","e_p","e_theta","e_phi",
                "p_p","p_theta","p_phi",
                "x","y","z","Q2","Mx2"]
    edges = [0.07,0.12,0.17,0.22,0.27,0.32,0.37,0.42,0.47,0.52,0.57]

    # Parallel execution
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(analyze_period, name, path, branches, edges): name for name, path in periods}
        for future in as_completed(futures):
            name = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"Error in {name}: {e}")
        #endfor

if __name__ == '__main__':
    main()
