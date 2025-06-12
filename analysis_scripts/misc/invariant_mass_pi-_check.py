#!/usr/bin/env python3
"""
Script: invariant_mass_pi-_check.py
Compute and plot the invariant mass of two charged pions (π⁻π⁺) for ep → enπ⁺
candidates, treating the scattered “electron” as a misidentified π⁻, in parallel for
RGC Su22, RGC Fa22, and RGC Sp23, with progress updates and optional test mode.

Usage:
  python invariant_mass_pi-_check.py [test]

If 'test' is provided, only the first 1,000,000 events per file are processed.
"""
import os
import sys
import time
import numpy as np
import uproot
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

# Check for test mode
TEST_MODE   = 'test' in sys.argv
ENTRY_STOP  = 1000000 if TEST_MODE else None

# Physical mass (GeV) for charged pions
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
    """Compute Mandelstam t = (q - p_pi)^2 (unchanged, for cuts)."""
    Eb = beam_energy(runnum)
    if Eb <= 0.0:
        return np.nan
    # Electron 4-vector (true)
    E_e  = np.sqrt(e_p**2 + 0.000511**2)
    sin_e, cos_e = np.sin(e_th), np.cos(e_th)
    ex = e_p * sin_e * np.cos(e_ph)
    ey = e_p * sin_e * np.sin(e_ph)
    ez = e_p * cos_e
    # Pion 4-vector
    E_pi = np.sqrt(p_p**2 + m_pi**2)
    sin_p, cos_p = np.sin(p_th), np.cos(p_th)
    px = p_p * sin_p * np.cos(p_ph)
    py = p_p * sin_p * np.sin(p_ph)
    pz = p_p * cos_p
    # Virtual photon q = beam - e'
    E_q = Eb - E_e
    qx, qy, qz = -ex, -ey, Eb - ez
    # Δ = q - p_pi
    dE = E_q - E_pi
    dx, dy, dz = qx - px, qy - py, qz - pz
    return dE*dE - (dx*dx + dy*dy + dz*dz)

#===============================================================================
def analyze_period(period_name, file_path, branches, edges):
    """Analyze one run period and save π⁻π⁺ invariant mass canvas with progress."""
    start_time = time.time()
    print(f"Starting {period_name} analysis{' (test mode)' if TEST_MODE else ''}...")
    # Load data
    with uproot.open(file_path) as f:
        if "PhysicsEvents" not in f:
            raise KeyError(f"Tree 'PhysicsEvents' not found in {file_path}")
        data = f["PhysicsEvents"].arrays(branches, library="np", entry_stop=ENTRY_STOP)
    runnum = data["runnum"]
    n_events = len(runnum)
    print(f"Loaded {n_events} events for {period_name}.")
    # Extract kinematics
    e_p  = data["e_p"];    e_th = data["e_theta"]; e_ph = data["e_phi"]
    p_p_ = data["p_p"];    p_th = data["p_theta"]; p_ph = data["p_phi"]
    x    = data["x"];      y    = data["y"]
    z    = data["z"];      Q2   = data["Q2"]
    Mx2  = data["Mx2"]
    # Compute t and apply cuts
    print("Computing t-values and applying cuts...")
    t_vals = np.vectorize(compute_t)(runnum, e_p, e_th, e_ph, p_p_, p_th, p_ph)
    mask = (
        (np.abs(t_vals) > 0.07) & (np.abs(t_vals) < 0.7) &
        (y < 0.65) & (z > 0.55) &
        (Q2 < 8.0) & (Mx2 > 0.75) & (Mx2 < 1.05)
    )
    n_sel = mask.sum()
    print(f"Events after cuts: {n_sel} ({n_sel/n_events*100:.1f}%)")
    # Precompute π⁻ and π⁺ four-vectors for selected events
    e_p_m  = e_p[mask];  e_th_m = e_th[mask]; e_ph_m = e_ph[mask]
    p_p_m  = p_p_[mask]; p_th_m = p_th[mask]; p_ph_m = p_ph[mask]
    x_m    = x[mask]
    # Treat electron as π⁻: use m_pi for both
    E_1 = np.sqrt(e_p_m**2 + m_pi**2)
    sin1, cos1 = np.sin(e_th_m), np.cos(e_th_m)
    px1 = e_p_m * sin1 * np.cos(e_ph_m)
    py1 = e_p_m * sin1 * np.sin(e_ph_m)
    pz1 = e_p_m * cos1
    # True π⁺
    E_2 = np.sqrt(p_p_m**2 + m_pi**2)
    sin2, cos2 = np.sin(p_th_m), np.cos(p_th_m)
    px2 = p_p_m * sin2 * np.cos(p_ph_m)
    py2 = p_p_m * sin2 * np.sin(p_ph_m)
    pz2 = p_p_m * cos2
    # Sum four-vectors and compute invariant mass
    E_tot  = E_1 + E_2
    px_tot = px1 + px2; py_tot = py1 + py2; pz_tot = pz1 + pz2
    inv_mass = np.sqrt(E_tot**2 - (px_tot**2 + py_tot**2 + pz_tot**2))
    # Bin assignment
    bin_idx = np.digitize(x_m, edges) - 1
    n_bins  = len(edges) - 1
    # Plotting
    fig, axes = plt.subplots(3, 4, figsize=(16, 12), sharex=True, sharey=True)
    axes_flat = axes.flatten()
    for i, ax in enumerate(axes_flat):
        if i < n_bins:
            sel = (bin_idx == i)
            cnt = sel.sum()
            print(f"  Bin {i+1}/{n_bins} ({edges[i]:.2f}<x<{edges[i+1]:.2f}): {cnt} events")
            ax.hist(inv_mass[sel], bins=50, range=(0, 3.5))
            ax.set_title(f"{edges[i]:.2f} < x < {edges[i+1]:.2f}")
            ax.set_xlim(0, 3.5)
        else:
            ax.axis('off')
        ax.set_xlabel(r"$M_{\pi\pi}$ (GeV)")
        ax.set_ylabel("Counts")
    plt.suptitle(
        fr"{period_name}: π⁻π⁺ invariant mass", fontsize=16, y=0.96
    )
    plt.tight_layout(rect=[0,0,1,0.94])
    out_dir = "output/enpi+"
    os.makedirs(out_dir, exist_ok=True)
    fn = period_name.lower().replace(' ', '_')
    fig.savefig(os.path.join(out_dir, f"pi-pi_mass_{fn}.pdf"))
    plt.close(fig)
    print(f"Finished {period_name} in {time.time()-start_time:.1f}s")

#===============================================================================
if __name__ == '__main__':
    branches = ["runnum","e_p","e_theta","e_phi",
                "p_p","p_theta","p_phi",
                "x","y","z","Q2","Mx2"]
    edges = [0.07,0.12,0.17,0.22,0.27,0.32,0.37,0.42,0.47,0.52,0.57]
    periods = [
        ("RGC Su22", "/work/.../rgc_su22_inb_NH3_epi+.root"),
        ("RGC Fa22", "/work/.../rgc_fa22_inb_NH3_epi+.root"),
        ("RGC Sp23", "/work/.../rgc_sp23_inb_NH3_epi+.root"),
    ]
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(analyze_period, name, path, branches, edges): name
                   for name, path in periods}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error in {futures[future]}: {e}")
