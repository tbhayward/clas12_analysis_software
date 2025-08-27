#!/usr/bin/env python3
"""
Compare baseline (no inverse-ISR) vs ISR-corrected (q -> q - R) outputs.

- Accepts text files with the 38-column layout that processing_two_particles.groovy writes,
  OR ROOT files with equivalent branches.
- Overlays Q2, x, W (baseline vs corrected).
- Shows ΔQ2 = Q2_corr - Q2_base matched on (runnum, evnum).
- Shows the sampled R-photon kinematics from the corrected set (Egamma, theta, phi).
- Plots k_p (a.k.a. e_p) with log y-scale and the label "k_p (GeV)".

Usage:
  python validate_inverse_isr.py --baseline baseline.txt --corrected corrected.txt --outdir output

If a path ends with .root, a ROOT TTree named "PhysicsEvents" is expected with branches:
  runnum, evnum, e_p, e_theta, e_phi, Q2, W, x, y, Egamma, isrTheta, isrPhi
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

try:
    import uproot  # optional, only used if you pass .root files
    HAS_UPROOT = True
except Exception:
    HAS_UPROOT = False
# endif

# -------------------------
# Column map for the 38-col text output (1-indexed in your printout, 0-indexed here)
# -------------------------
COL = {
    "fiducial_status": 0,
    "num_pos": 1,
    "num_neg": 2,
    "num_neutrals": 3,
    "runnum": 4,
    "evnum": 5,
    "helicity": 6,
    "detector": 7,
    "e_p": 8,
    "e_theta": 9,    # radians in your writer
    "e_phi": 10,     # radians in your writer
    "vz_e": 11,
    "p_p": 12,
    "p_theta": 13,   # radians
    "p_phi": 14,     # radians
    "vz_p": 15,
    "open_angle": 16,
    "Egamma": 17,    # GeV (0.0 for baseline)
    "isrTheta": 18,  # deg
    "isrPhi": 19,    # deg
    "Q2": 20,
    "W": 21,
    "Mx2": 22,
    "x": 23,
    "t": 24,
    "tmin": 25,
    "y": 26,
    "z": 27,
    "xF": 28,
    "pT": 29,
    "xi": 30,
    "eta": 31,
    "phi_trento": 32,  # radians
    "DepA": 33,
    "DepB": 34,
    "DepC": 35,
    "DepV": 36,
    "DepW": 37,
}

# -------------------------
# I/O helpers
# -------------------------
def load_text(path):
    arr = np.loadtxt(path)
    data = {k: arr[:, idx] for k, idx in COL.items()}
    # types
    data["runnum"] = data["runnum"].astype(np.int64)
    data["evnum"]  = data["evnum"].astype(np.int64)
    return data

def load_root(path):
    if not HAS_UPROOT:
        raise RuntimeError("uproot is not available, but a ROOT file was provided.")
    with uproot.open(path) as f:
        tree = f["PhysicsEvents"]
        needed = ["runnum", "evnum", "e_p", "e_theta", "e_phi", "Q2", "W", "x", "y",
                  "Egamma", "isrTheta", "isrPhi"]
        avail = [b.decode() if isinstance(b, bytes) else b for b in tree.keys()]
        # Basic presence check (be permissive)
        missing = [n for n in needed if n not in avail]
        if missing:
            raise KeyError(f"Missing branches in ROOT file {path}: {missing}")
        # Read as numpy
        arrs = tree.arrays(needed, library="np")
        data = {k: arrs[k].astype(np.float64) for k in needed}
        data["runnum"] = data["runnum"].astype(np.int64)
        data["evnum"]  = data["evnum"].astype(np.int64)
        # Fill dummies that the text path would have
        data.setdefault("Mx2", np.full_like(data["Q2"], np.nan))
        data.setdefault("t",   np.full_like(data["Q2"], np.nan))
        data.setdefault("tmin",np.full_like(data["Q2"], np.nan))
        data.setdefault("z",   np.full_like(data["Q2"], np.nan))
        data.setdefault("xF",  np.full_like(data["Q2"], np.nan))
        data.setdefault("pT",  np.full_like(data["Q2"], np.nan))
        data.setdefault("xi",  np.full_like(data["Q2"], np.nan))
        data.setdefault("eta", np.full_like(data["Q2"], np.nan))
        return data
    # endif
# enddef

def load_any(path):
    if path.lower().endswith(".root"):
        return load_root(path)
    else:
        return load_text(path)
    # endif

def degify(arr_rad):
    """Radians -> degrees, robust to NaNs."""
    return np.rad2deg(arr_rad)

# -------------------------
# Matching on (runnum, evnum) for Δ-distributions
# -------------------------
def match_by_keys(base, corr):
    """
    Return indices (i_base, i_corr) for rows with matching (runnum, evnum).
    """
    base_keys = base["runnum"].astype(np.int64) * (10**10) + base["evnum"].astype(np.int64)
    corr_keys = corr["runnum"].astype(np.int64) * (10**10) + corr["evnum"].astype(np.int64)

    # Sort & match
    idx_b = np.argsort(base_keys)
    idx_c = np.argsort(corr_keys)
    bk = base_keys[idx_b]
    ck = corr_keys[idx_c]

    # Two-pointer intersection
    i = j = 0
    sel_b = []
    sel_c = []
    while i < bk.size and j < ck.size:
        if bk[i] == ck[j]:
            sel_b.append(idx_b[i])
            sel_c.append(idx_c[j])
            i += 1
            j += 1
        elif bk[i] < ck[j]:
            i += 1
        else:
            j += 1
        # endif
    # endwhile

    return np.array(sel_b, dtype=np.int64), np.array(sel_c, dtype=np.int64)
# enddef

# -------------------------
# Plotting
# -------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True, help="Path to baseline file (.txt or .root)")
    ap.add_argument("--corrected", required=True, help="Path to ISR-corrected file (.txt or .root)")
    ap.add_argument("--outdir", default="output/validation_inverseISR", help="Output directory")
    ap.add_argument("--tag", default="", help="Optional tag for output filename")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    base = load_any(args.baseline)
    corr = load_any(args.corrected)

    # Convenience views
    base_ep   = base["e_p"]
    corr_ep   = corr["e_p"]
    base_ethd = degify(base.get("e_theta", np.nan))  # radians -> deg
    base_ephd = degify(base.get("e_phi",   np.nan))
    corr_ethd = degify(corr.get("e_theta", np.nan))
    corr_ephd = degify(corr.get("e_phi",   np.nan))

    base_Q2 = base["Q2"];    corr_Q2 = corr["Q2"]
    base_x  = base["x"];     corr_x  = corr["x"]
    base_W  = base["W"];     corr_W  = corr["W"]
    base_y  = base["y"];     corr_y  = corr["y"]

    corr_Eg = corr.get("Egamma", np.zeros_like(corr_Q2))
    corr_th = corr.get("isrTheta", np.zeros_like(corr_Q2))  # already in deg for text path
    corr_ph = corr.get("isrPhi",   np.zeros_like(corr_Q2))  # deg

    # Align for Δ distributions
    ib, ic = match_by_keys(base, corr)
    dQ2 = corr_Q2[ic] - base_Q2[ib]

    # -------------------------
    # Figure: 2x3 comparison
    # -------------------------
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    (ax_kp, ax_q2, ax_x), (ax_dq2, ax_eg, ax_th) = axes

    # (1) k_p (GeV) — baseline vs corrected (these should lie on top of each other)
    n_b, bins, _ = ax_kp.hist(base_ep, bins=100, histtype="step", label="baseline")
    n_c, _,    _ = ax_kp.hist(corr_ep, bins=bins, histtype="step", label="ISR corrected")
    ax_kp.set_xlabel(r"$k_p$ (GeV)")
    ax_kp.set_ylabel("Counts")
    ax_kp.set_yscale("log")
    # set a sensible bottom so lines don't get cut off on log scale
    nz = np.concatenate([n_b[n_b>0], n_c[n_c>0]])
    if nz.size > 0:
        ax_kp.set_ylim(bottom=max(1.0, 0.8 * nz.min()))
    ax_kp.legend(loc="best", frameon=False)

    # (2) Q2 overlay
    ax_q2.hist(base_Q2, bins=120, histtype="step", label="baseline")
    ax_q2.hist(corr_Q2, bins=120, histtype="step", label="ISR corrected")
    ax_q2.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax_q2.set_ylabel("Counts")
    ax_q2.legend(loc="best", frameon=False)

    # (3) x overlay
    ax_x.hist(base_x, bins=120, histtype="step", label="baseline")
    ax_x.hist(corr_x, bins=120, histtype="step", label="ISR corrected")
    ax_x.set_xlabel(r"$x$")
    ax_x.legend(loc="best", frameon=False)

    # (4) ΔQ2 = Q2_corr - Q2_base (matched)
    ax_dq2.hist(dQ2, bins=160, histtype="step")
    ax_dq2.set_xlabel(r"$\Delta Q^2 \equiv Q^2_{\rm corr} - Q^2_{\rm base}$ (GeV$^2$)")
    ax_dq2.set_ylabel("Counts")

    # (5) Sampled R photon energy (corrected only)
    ax_eg.hist(corr_Eg, bins=120, histtype="step")
    ax_eg.set_xlabel(r"$E_\gamma^{(R)}$ (GeV)")

    # (6) Sampled R polar angle (deg) (corrected only)
    ax_th.hist(corr_th, bins=100, histtype="step")
    ax_th.set_xlabel(r"$\theta_R$ (deg)")

    # Style
    for ax in axes.flat:
        ax.grid(True, alpha=0.25)
    # endfor

    plt.tight_layout()
    tag = f"_{args.tag}" if args.tag else ""
    outpath = os.path.join(args.outdir, f"inverseISR_validation{tag}.pdf")
    plt.savefig(outpath)
    print(f"Saved: {outpath}")

    # -------------------------
    # Optional: a small angles figure if you want φ_R as well
    # -------------------------
    fig2, (ax_phiR, ax_eAng) = plt.subplots(1, 2, figsize=(12, 4.2))

    # φ_R (deg)
    # keep it in [0,360)
    phiR = np.mod(corr_ph, 360.0)
    ax_phiR.hist(phiR, bins=np.linspace(0, 360, 121), histtype="step")
    ax_phiR.set_xlabel(r"$\phi_R$ (deg)")
    ax_phiR.grid(True, alpha=0.25)

    # Electron angular distributions (baseline vs corrected) for sanity
    # (should be identical; just here for quick QA)
    ax_eAng.hist(base_ethd, bins=120, histtype="step", label="e_θ baseline")
    ax_eAng.hist(corr_ethd, bins=120, histtype="step", label="e_θ corrected")
    ax_eAng.set_xlabel(r"$e_{\theta}$ (deg)")
    ax_eAng.legend(loc="best", frameon=False)
    ax_eAng.grid(True, alpha=0.25)

    plt.tight_layout()
    outpath2 = os.path.join(args.outdir, f"angles_supplement{tag}.pdf")
    plt.savefig(outpath2)
    print(f"Saved: {outpath2}")
# enddef

if __name__ == "__main__":
    main()
# endif