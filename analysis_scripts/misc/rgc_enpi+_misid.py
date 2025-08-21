#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Misidentification study (MC): fractions vs -t in four xB slices.

Definition (TEST MODE, per user request):
  e_is_pi_minus = (ep != 11)        # electron candidate actually π− (or anything not 11)
  p1_is_k_minus = (p1 != 211)       # hadron candidate flagged as K− (anything not 211)

Cuts:
  - Only 0.81 < Mx2 < 1.00
  - x_B in one of four bins:
      Low:     0.10 < x < 0.25
      MidLow:  0.25 < x < 0.35
      MidHigh: 0.35 < x < 0.45
      High:    0.45 < x < 0.60
  - -t binned with edges: [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]
    (branch 't' is negated for plotting/bins)

Input:
  python misid_test.py input.root

Output:
  Figure: output/enpi+/misidentification.pdf
  Console: per-bin table + lists of offending ep/p1 values.

Notes:
  - If branches 'ep'/'p1' are missing, we fall back to 'matching_e_pid'/'matching_p1_pid'.
"""

import os
import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt

# ------------------------------
# Config
# ------------------------------
X_SLICES = [
    ("Low",     (0.10, 0.25)),
    ("MidLow",  (0.25, 0.35)),
    ("MidHigh", (0.35, 0.45)),
    ("High",    (0.45, 0.60)),
]
T_EDGES_POS = np.array([0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25], dtype=float)
YLIM = (0.0, 0.01)
Y_LABEL = "fraction (normalized counts)"
X_LABEL = r"$-t\ \mathrm{(GeV^2)}$"
TITLES = {
    "Low":     r"$0.10 < x_{B} < 0.25$",
    "MidLow":  r"$0.25 < x_{B} < 0.35$",
    "MidHigh": r"$0.35 < x_{B} < 0.45$",
    "High":    r"$0.45 < x_{B} < 0.60$",
}

# Markers/colors
STYLE_E = dict(marker="o", linestyle="none", color="tab:blue",  label=r"$\pi^- \rightarrow e^-$",  ms=5)
STYLE_K = dict(marker="s", linestyle="none", color="tab:orange", label=r"$\pi^- \rightarrow K^-$", ms=5)

# ------------------------------
# Utilities
# ------------------------------
def get_branch(tree, primary, fallback=None):
    """Return the available branch name (primary if exists else fallback)."""
    names = tree.keys()
    if primary in names:
        return primary, None
    if fallback and fallback in names:
        return fallback, f"[WARN] Using fallback branch '{fallback}' (missing '{primary}')."
    return None, f"[ERROR] Neither '{primary}' nor fallback '{fallback}' found."

def hist_fraction(values_mask, base_mask, tpos, edges):
    """
    Compute per-bin fraction = (# values_mask & base in bin) / (# base in bin).
    Returns (centers, fractions, N_total, N_num) where totals are integer arrays.
    """
    # bin indices for all base events
    bidx = np.digitize(tpos, edges) - 1  # 0..nb-1
    nb = len(edges) - 1

    N_total = np.zeros(nb, dtype=int)
    N_num   = np.zeros(nb, dtype=int)
    for b in range(nb):
        in_bin = base_mask & (bidx == b)
        N_total[b] = int(np.count_nonzero(in_bin))
        if N_total[b] > 0:
            N_num[b] = int(np.count_nonzero(values_mask & (bidx == b)))
    frac = np.divide(N_num, N_total, out=np.zeros_like(N_total, dtype=float), where=N_total>0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, frac, N_total, N_num

def main():
    if len(sys.argv) != 2:
        print("Usage: python misid_test.py input.root")
        sys.exit(1)

    path = sys.argv[1]
    if not os.path.isfile(path):
        print(f"[ERROR] File not found: {path}")
        sys.exit(1)

    # Open tree
    with uproot.open(path) as f:
        if "PhysicsEvents" not in f:
            print("[ERROR] Tree 'PhysicsEvents' not found.")
            sys.exit(1)
        t = f["PhysicsEvents"]

        # Resolve PID branches per request (primary: 'ep' and 'p1')
        ep_branch, warn_ep = get_branch(t, "ep", "matching_e_pid")
        p1_branch, warn_p1 = get_branch(t, "p1", "matching_p1_pid")
        if warn_ep: print(warn_ep)
        if warn_p1: print(warn_p1)
        if ep_branch is None or p1_branch is None:
            sys.exit(2)

        # Load arrays (as numpy)
        needed = ["x", "t", "Mx2", ep_branch, p1_branch]
        arrays = t.arrays(needed, library="np")

    x   = arrays["x"].astype(float)
    tbr = arrays["t"].astype(float)
    Mx2 = arrays["Mx2"].astype(float)
    ep  = arrays[ep_branch].astype(np.int32)
    p1  = arrays[p1_branch].astype(np.int32)

    # Base masks: finite + Mx2 window + t finite
    finite = np.isfinite(x) & np.isfinite(tbr) & np.isfinite(Mx2)
    base_global = finite & (Mx2 > 0.81) & (Mx2 < 1.00)

    # -t for binning
    tpos = -tbr

    # Prepare figure
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2), sharey=True)
    fig.suptitle(r"$ep \rightarrow en\pi^{+}$ mis-ID test (MC): "
                 r"$0.81<M_{x}^{2}<1.00$; fractions vs $-t$", fontsize=14, y=0.98)

    # Header for console
    print("\n=== Mis-ID fractions by xB slice and -t bin ===")
    print("(-t bins):", ", ".join([f"[{T_EDGES_POS[i]:.2f},{T_EDGES_POS[i+1]:.2f})"
                                   for i in range(len(T_EDGES_POS)-1)]))
    print("Definition: e_is_pi_minus = (ep != 11),  p1_is_k_minus = (p1 != 211)\n")

    for ax, (tag, (xa, xb)) in zip(axes, X_SLICES):
        # xB slice+global base
        base_slice = base_global & (x > xa) & (x < xb) & (tpos > T_EDGES_POS[0]) & (tpos < T_EDGES_POS[-1])

        # TEST CONDITIONS (per-user):
        e_is_pi_minus = (ep != 11)
        p1_is_k_minus = (p1 != 211)

        # Fractions per -t bin
        xc, frac_e, Ntot, Ne = hist_fraction(e_is_pi_minus, base_slice, tpos, T_EDGES_POS)
        _,  frac_k, _,   Nk  = hist_fraction(p1_is_k_minus, base_slice, tpos, T_EDGES_POS)

        # Plot
        ax.errorbar(xc, frac_e, yerr=None, **STYLE_E)
        ax.errorbar(xc, frac_k, yerr=None, **STYLE_K)
        ax.set_xlim(T_EDGES_POS[0], T_EDGES_POS[-1])
        ax.set_ylim(*YLIM)
        ax.set_xlabel(X_LABEL)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.set_title(TITLES[tag], fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(Y_LABEL)
        if ax is axes[-1]:
            ax.legend(frameon=True, edgecolor="black", loc="upper right", fontsize=10)

        # Console printout for this slice
        print(f"\n-- xB slice {tag}: {TITLES[tag]} --")
        print("bin  [-t_min,-t_max)    N_total   N(pi->e)   frac(pi->e)   N(pi->K)   frac(pi->K)")
        for i in range(len(T_EDGES_POS)-1):
            tmin, tmax = T_EDGES_POS[i], T_EDGES_POS[i+1]
            print(f"{i+1:>3d}  [{tmin:5.2f},{tmax:5.2f})   {Ntot[i]:7d}   {Ne[i]:8d}   {frac_e[i]:11.6f}   {Nk[i]:8d}   {frac_k[i]:11.6f}")

            # Print actual offending values in this bin
            if Ntot[i] > 0:
                bmask = base_slice & (tpos >= tmin) & (tpos < tmax)
                # e: ep != 11
                hits_e = np.where(bmask & (ep != 11))[0]
                if hits_e.size:
                    vals = ep[hits_e]
                    print(f"    ep values != 11 in this bin (count {hits_e.size}): {vals.tolist()}")
                # k: p1 != 211
                hits_k = np.where(bmask & (p1 != 211))[0]
                if hits_k.size:
                    vals = p1[hits_k]
                    print(f"    p1 values != 211 in this bin (count {hits_k.size}): {vals.tolist()}")

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    outdir = os.path.join("output", "enpi+")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "misidentification.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"\nSaved: {outpath}")

if __name__ == "__main__":
    main()