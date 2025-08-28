#!/usr/bin/env python3
"""
Make inverseISR_validation.pdf with:
Top:   R_p (GeV), R_theta (deg), R_phi (deg)
Bottom: ΔQ^2 (GeV^2), Δx_B, ΔMx^2 (GeV^2)

- Baseline vs ISR-corrected files (text or ROOT) are matched by (runnum, evnum).
- Text format is your 38-column dump; ROOT expects branches in "PhysicsEvents":
  runnum, evnum, Q2, W, x, y, Mx2, e_p, e_theta, e_phi, Egamma, isrTheta, isrPhi.
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

def nonzero_hist_range(x, bins):
    """Return (lo_edge, hi_edge) spanning only non-empty bins; None if no entries."""
    n, edges = np.histogram(x, bins=bins)
    nz = np.where(n > 0)[0]
    if nz.size == 0:
        return None
    lo = edges[nz.min()]
    hi = edges[nz.max() + 1]
    return lo, hi

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
    "t": 24,
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
        needed = ["runnum", "evnum", "Q2", "x", "Mx2", "Egamma", "isrTheta", "isrPhi"]
        arrs = tree.arrays(needed, library="np")
        data = {k: arrs[k].astype(np.float64) for k in needed}
        data["runnum"] = data["runnum"].astype(np.int64)
        data["evnum"]  = data["evnum"].astype(np.int64)
        return data

def load_any(path):
    return load_root(path) if path.lower().endswith(".root") else load_text(path)

# --------------- matching ---------------
def match_by_keys(base, corr):
    """Return indices (i_base, i_corr) for matching (runnum, evnum)."""
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
            sel_b.append(ib[i]); sel_c.append(ic[j])
            i += 1; j += 1
        elif bk[i] < ck[j]:
            i += 1
        else:
            j += 1
    return np.array(sel_b, dtype=np.int64), np.array(sel_c, dtype=np.int64)

# --------------- plotting helper ---------------
def plot_hist_points(ax, data, bins, rng=None, xlabel="", logy=False):
    """Make a binned histogram and draw as points with √N error bars.
       For log-y, zero-count bins are dropped."""
    if data.size == 0:
        ax.text(0.5, 0.5, "No entries", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    n, edges = np.histogram(data, bins=bins, range=rng)
    centers = 0.5 * (edges[:-1] + edges[1:])
    if logy:
        mask = n > 0
    else:
        mask = n >= 0  # keep zeros on linear plots

    if not np.any(mask):
        ax.text(0.5, 0.5, "All bins are zero", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    ax.errorbar(centers[mask], n[mask], yerr=np.sqrt(n[mask]), fmt='o', ms=3, lw=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")
    if rng is not None:
        ax.set_xlim(*rng)
    if logy:
        ax.set_yscale("log")
        # put a sensible bottom just below the smallest positive count
        ymin = np.min(n[mask])
        ax.set_ylim(bottom=max(1.0, 0.8 * ymin))

    ax.grid(True, alpha=0.25)

# --------------- main ---------------
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
    Q2_b, Q2_c   = base["Q2"][ib],  corr["Q2"][ic]
    x_b,  x_c    = base["x"][ib],   corr["x"][ic]
    Mx2_b, Mx2_c = base["Mx2"][ib], corr["Mx2"][ic]

    mask_q2 = np.isfinite(Q2_b)  & np.isfinite(Q2_c)
    mask_x  = np.isfinite(x_b)   & np.isfinite(x_c)
    mask_mx = np.isfinite(Mx2_b) & np.isfinite(Mx2_c)

    dQ2  = (Q2_c  - Q2_b )[mask_q2]
    dxB  = (x_c   - x_b  )[mask_x ]
    dMx2 = (Mx2_c - Mx2_b)[mask_mx]

    # R photon (from corrected)
    Rp     = corr.get("Egamma",   np.array([]))  # GeV
    Rtheta = corr.get("isrTheta", np.array([]))  # deg
    Rphi   = corr.get("isrPhi",   np.array([]))  # deg
    if Rphi.size:
        Rphi = np.mod(Rphi, 360.0)

    # figure layout: 2x3
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    (ax_Rp, ax_Rth, ax_Rph), (ax_dQ2, ax_dxB, ax_dMx2) = axes

    # --- top row: R_p, R_theta, R_phi ---
    plot_hist_points(ax_Rp,   Rp,     bins=60, rng=None,                   xlabel=r"$R_p$ (GeV)",        logy=True)
    plot_hist_points(ax_Rth,  Rtheta, bins=50, rng=None,                   xlabel=r"$R_{\theta}$ (deg)", logy=False)

    # R_phi: force x in [0, 360], y = [0.5*min_nonzero, 1.5*max]
    if Rphi.size:
        bins_Rphi = np.linspace(0.0, 360.0, 61)  # 60 bins
        plot_hist_points(ax_Rph, Rphi, bins=bins_Rphi, rng=(0.0, 360.0),
                         xlabel=r"$R_{\phi}$ (deg)", logy=False)

        n_phi, _ = np.histogram(Rphi, bins=bins_Rphi)
        nz = n_phi[n_phi > 0]
        if nz.size:
            ymin = 0.5 * nz.min()
            ymax = 1.5 * n_phi.max()
            if ymax <= ymin:
                ymax = ymin + 1.0
            ax_Rph.set_ylim(ymin, ymax)
    else:
        ax_Rph.text(0.5, 0.5, "No $R_{\\phi}$ info", ha="center", va="center", transform=ax_Rph.transAxes)
        ax_Rph.set_axis_off()

    # --- bottom row: ΔQ2, Δx_B, ΔMx^2 ---
    plot_hist_points(ax_dQ2,  dQ2,   bins=80,  rng=(-4.0, 1.0), xlabel=r"$\Delta Q^2$ (GeV$^2$)", logy=True)
    plot_hist_points(ax_dxB,  dxB,   bins=80,  rng=(-0.3, 0.3), xlabel=r"$\Delta x_B$",           logy=True)
    plot_hist_points(ax_dMx2, dMx2,  bins=80,  rng=(-9.0, 3.0), xlabel=r"$\Delta M_x^2$ (GeV$^2$)",logy=True)

    plt.tight_layout()
    tag = f"_{args.tag}" if args.tag else ""
    outpath = os.path.join(args.outdir, f"inverseISR_validation{tag}.pdf")
    plt.savefig(outpath)
    print(f"Saved: {outpath}")

if __name__ == "__main__":
    main()