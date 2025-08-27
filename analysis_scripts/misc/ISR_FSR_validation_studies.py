#!/usr/bin/env python3
"""
Make inverseISR_validation.pdf with:
Top:   R_p (GeV), R_theta (deg), R_phi (deg)
Bottom: ΔQ^2 (GeV^2), Δx_B, Δ(-t) (GeV^2)

- Baseline vs ISR-corrected files (text or ROOT) are matched by (runnum, evnum).
- Text format is your 38-column dump; ROOT expects branches in "PhysicsEvents":
  runnum, evnum, Q2, W, x, y, t, e_p, e_theta, e_phi, Egamma, isrTheta, isrPhi.

Usage:
  python validate_inverse_isr.py --baseline baseline.txt --corrected corrected.txt --outdir output
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

try:
    import uproot
    HAS_UPROOT = True
except Exception:
    HAS_UPROOT = False
# endif

# -------- column mapping for 38-col text (0-based) --------
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
    "e_theta": 9,   # radians
    "e_phi": 10,    # radians
    "vz_e": 11,
    "p_p": 12,
    "p_theta": 13,
    "p_phi": 14,
    "vz_p": 15,
    "open_angle": 16,
    "Egamma": 17,       # GeV (0 for baseline)
    "isrTheta": 18,     # deg
    "isrPhi": 19,       # deg
    "Q2": 20,
    "W": 21,
    "Mx2": 22,
    "x": 23,
    "t": 24,            # GeV^2 (your sign convention); we will use -t
    "tmin": 25,
    "y": 26,
    "z": 27,
    "xF": 28,
    "pT": 29,
    "xi": 30,
    "eta": 31,
    "phi_trento": 32,
    "DepA": 33,
    "DepB": 34,
    "DepC": 35,
    "DepV": 36,
    "DepW": 37,
}

# --------------- loaders ---------------
def load_text(path):
    arr = np.loadtxt(path)
    data = {k: arr[:, i] for k, i in COL.items()}
    data["runnum"] = data["runnum"].astype(np.int64)
    data["evnum"]  = data["evnum"].astype(np.int64)
    return data

def load_root(path):
    if not HAS_UPROOT:
        raise RuntimeError("uproot is not available but a ROOT file was provided.")
    with uproot.open(path) as f:
        tree = f["PhysicsEvents"]
        # minimal set we actually need here
        needed = ["runnum", "evnum", "Q2", "x", "t", "Egamma", "isrTheta", "isrPhi"]
        arrs = tree.arrays(needed, library="np")
        data = {k: arrs[k].astype(np.float64) for k in needed}
        data["runnum"] = data["runnum"].astype(np.int64)
        data["evnum"]  = data["evnum"].astype(np.int64)
        return data
    # endif
# endif

def load_any(path):
    return load_root(path) if path.lower().endswith(".root") else load_text(path)

# --------------- matching ---------------
def match_by_keys(base, corr):
    """
    Return indices (i_base, i_corr) for matching (runnum, evnum).
    """
    base_keys = base["runnum"] * (10**10) + base["evnum"]
    corr_keys = corr["runnum"] * (10**10) + corr["evnum"]

    ib = np.argsort(base_keys)
    ic = np.argsort(corr_keys)
    bk = base_keys[ib]
    ck = corr_keys[ic]

    i = j = 0
    sel_b, sel_c = [], []
    while i < bk.size and j < ck.size:
        if bk[i] == ck[j]:
            sel_b.append(ib[i])
            sel_c.append(ic[j])
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

# --------------- plotting ---------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True, help="Baseline file (.txt or .root)")
    ap.add_argument("--corrected", required=True, help="ISR-corrected file (.txt or .root)")
    ap.add_argument("--outdir", default="output/validation_inverseISR", help="Output directory")
    ap.add_argument("--tag", default="", help="Optional tag for output filename")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    base = load_any(args.baseline)
    corr = load_any(args.corrected)

    # match events
    ib, ic = match_by_keys(base, corr)

    # deltas (guard NaNs)
    Q2_b, Q2_c = base["Q2"][ib], corr["Q2"][ic]
    x_b,  x_c  = base["x"][ib],  corr["x"][ic]
    t_b,  t_c  = base["t"][ib],  corr["t"][ic]

    mask_q2 = np.isfinite(Q2_b) & np.isfinite(Q2_c)
    mask_x  = np.isfinite(x_b)  & np.isfinite(x_c)
    mask_t  = np.isfinite(t_b)  & np.isfinite(t_c)

    dQ2       = (Q2_c - Q2_b)[mask_q2]
    dxB       = (x_c  - x_b )[mask_x ]
    d_minus_t = ((-t_c) - (-t_b))[mask_t]  # Δ(-t)

    # R photon (from corrected)
    Rp     = corr.get("Egamma", np.array([]))          # GeV
    Rtheta = corr.get("isrTheta", np.array([]))        # degrees
    Rphi   = corr.get("isrPhi",   np.array([]))        # degrees
    if Rphi.size:
        Rphi = np.mod(Rphi, 360.0)
    # endif

    # figure layout: 2x3
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    (ax_Rp, ax_Rth, ax_Rph), (ax_dQ2, ax_dxB, ax_dmt) = axes

    # --- top row: R_p, R_theta, R_phi ---
    if Rp.size:
        n_rp, bins, _ = ax_Rp.hist(Rp, bins=120, histtype="step")
        ax_Rp.set_xlabel(r"$R_p$ (GeV)")
        ax_Rp.set_ylabel("Counts")
        ax_Rp.set_yscale("log")
        pos = n_rp[n_rp > 0]
        if pos.size > 0:
            ax_Rp.set_ylim(bottom=max(1.0, 0.8 * pos.min()))
    else:
        ax_Rp.text(0.5, 0.5, "No $R_p$ info", ha="center", va="center", transform=ax_Rp.transAxes)
        ax_Rp.set_axis_off()
    # endif

    if Rtheta.size:
        ax_Rth.hist(Rtheta, bins=100, histtype="step")
        ax_Rth.set_xlabel(r"$R_{\theta}$ (deg)")
    else:
        ax_Rth.text(0.5, 0.5, "No $R_{\\theta}$ info", ha="center", va="center", transform=ax_Rth.transAxes)
        ax_Rth.set_axis_off()
    # endif

    if Rphi.size:
        ax_Rph.hist(Rphi, bins=np.linspace(0, 360, 121), histtype="step")
        ax_Rph.set_xlabel(r"$R_{\phi}$ (deg)")
    else:
        ax_Rph.text(0.5, 0.5, "No $R_{\\phi}$ info", ha="center", va="center", transform=ax_Rph.transAxes)
        ax_Rph.set_axis_off()
    # endif

    # --- bottom row: ΔQ2, Δx_B, Δ(-t) ---  (all log y-scale)
    if dQ2.size:
        n_dq2, _, _ = ax_dQ2.hist(dQ2, bins=160, histtype="step")
        ax_dQ2.set_xlabel(r"$\Delta Q^2$ (GeV$^2$)")
        ax_dQ2.set_ylabel("Counts")
        ax_dQ2.set_yscale("log")
        pos = n_dq2[n_dq2 > 0]
        if pos.size > 0:
            ax_dQ2.set_ylim(bottom=max(1.0, 0.8 * pos.min()))
    else:
        ax_dQ2.text(0.5, 0.5, "No matched events for $\Delta Q^2$", ha="center", va="center", transform=ax_dQ2.transAxes)
        ax_dQ2.set_axis_off()
    # endif

    if dxB.size:
        n_dxb, _, _ = ax_dxB.hist(dxB, bins=160, histtype="step")
        ax_dxB.set_xlabel(r"$\Delta x_B$")
        ax_dxB.set_yscale("log")
        pos = n_dxb[n_dxb > 0]
        if pos.size > 0:
            ax_dxB.set_ylim(bottom=max(1.0, 0.8 * pos.min()))
    else:
        ax_dxB.text(0.5, 0.5, "No matched events for $\Delta x_B$", ha="center", va="center", transform=ax_dxB.transAxes)
        ax_dxB.set_axis_off()
    # endif

    if d_minus_t.size:
        n_dmt, _, _ = ax_dmt.hist(d_minus_t, bins=160, histtype="step")
        ax_dmt.set_xlabel(r"$\Delta(-t)$ (GeV$^2$)")
        ax_dmt.set_yscale("log")
        pos = n_dmt[n_dmt > 0]
        if pos.size > 0:
            ax_dmt.set_ylim(bottom=max(1.0, 0.8 * pos.min()))
    else:
        ax_dmt.text(0.5, 0.5, "No matched events for $\Delta(-t)$", ha="center", va="center", transform=ax_dmt.transAxes)
        ax_dmt.set_axis_off()
    # endif

    # style
    for ax in axes.flat:
        if ax.get_visible():
            ax.grid(True, alpha=0.25)
    #endfor

    plt.tight_layout()
    tag = f"_{args.tag}" if args.tag else ""
    outpath = os.path.join(args.outdir, f"inverseISR_validation{tag}.pdf")
    plt.savefig(outpath)
    print(f"Saved: {outpath}")
# endif

if __name__ == "__main__":
    main()
# endif