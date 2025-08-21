#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make 2x3 plots of <sin(theta_gamma)> vs phi for six -t bins ( Low x_B slice ).
Cuts:
  Q2 > 1, W > 2, y < 0.75, fiducial_status >= 111,
  0.10 < x < 0.30, 0.81 < Mx2 < 1.00
Input:
  python make_sTG_phi_panels.py <input.root>
Output:
  output/enpi+/sinthetagamma_bins.pdf
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")  # ensure no GUI pops up
import matplotlib.pyplot as plt

try:
    import uproot
except ImportError:
    print("[ERROR] Please install uproot:  pip install uproot")
    sys.exit(1)

# ---------------------- config ----------------------
TREE_NAME = "PhysicsEvents"
BRANCHES  = ["Q2","W","y","Mx2","x","t","phi","fiducial_status","sinthetagamma"]
ITER_STEP = "200 MB"

# Low xB slice and 6 -t bins (converted from t: -1.25 ... -0.05)
XB_MIN, XB_MAX = 0.10, 0.30
T_EDGES_POS = np.array([0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25], dtype=float)  # -t
NB_TBINS = len(T_EDGES_POS) - 1

# phi binning for the averages
NB_PHI = 12
PHI_EDGES = np.linspace(0.0, 2.0*np.pi, NB_PHI + 1)
PHI_CEN   = 0.5 * (PHI_EDGES[:-1] + PHI_EDGES[1:])

# ---------------------- helpers ----------------------
def ensure_outdir():
    outdir = os.path.join("output", "enpi+")
    os.makedirs(outdir, exist_ok=True)
    return outdir

def wrap_phi(phi):
    return np.mod(phi, 2.0*np.pi)

# ---------------------- accumulation ----------------------
def accumulate_means(infile):
    """
    Return arrays:
      s_sum[tbin, phibin], s_cnt[tbin, phibin]
    for computing <sin(theta_gamma)>(phi) in each -t bin.
    """
    s_sum = np.zeros((NB_TBINS, NB_PHI), dtype=np.float64)
    s_cnt = np.zeros((NB_TBINS, NB_PHI), dtype=np.int64)

    it = uproot.iterate(
        {infile: TREE_NAME},
        BRANCHES,
        step_size=ITER_STEP,
        library="np",
        allow_read_errors=False
    )

    for arr in it:
        Q2  = arr["Q2"]
        W   = arr["W"]
        y   = arr["y"]
        Mx2 = arr["Mx2"]
        x   = arr["x"]
        t   = arr["t"]
        phi = wrap_phi(arr["phi"])
        fid = arr["fiducial_status"]
        sTG = arr["sinthetagamma"]

        tpos = -t  # we work in -t

        # Base + exclusivity + fid + Low xB slice
        base = (Q2 > 1.0) & (W > 2.0) & (y < 0.75)
        excl = (Mx2 > 0.81) & (Mx2 < 1.00)
        xbok = (x > XB_MIN) & (x < XB_MAX)
        fidg = (fid >= 111)

        m = base & excl & xbok & fidg
        if not np.any(m):
            continue

        tpos = tpos[m]
        phi  = phi[m]
        sTG  = sTG[m]

        # Bin in -t (six bins)
        tb = np.digitize(tpos, T_EDGES_POS) - 1
        good_t = (tb >= 0) & (tb < NB_TBINS)

        if not np.any(good_t):
            continue

        tb   = tb[good_t]
        phi  = phi[good_t]
        sTG  = sTG[good_t]

        # Bin in phi
        pb = np.digitize(phi, PHI_EDGES) - 1
        good_p = (pb >= 0) & (pb < NB_PHI)

        tb = tb[good_p]
        pb = pb[good_p]
        s  = sTG[good_p]

        # Accumulate sums and counts
        np.add.at(s_sum, (tb, pb), s)
        np.add.at(s_cnt, (tb, pb), 1)

    return s_sum, s_cnt

# ---------------------- plotting ----------------------
def make_panel_plot(s_sum, s_cnt, outdir):
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey=True)
    axes = axes.ravel()

    for i in range(NB_TBINS):
        ax = axes[i]
        counts = s_cnt[i]
        sums   = s_sum[i]
        with np.errstate(invalid="ignore", divide="ignore"):
            means = np.where(counts > 0, sums / counts, np.nan)

        # scatter (no lines)
        ax.plot(PHI_CEN, means, "o", ms=4, color="black")

        # axis labels and ticks
        ax.set_xlim(0, 2*np.pi)
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(r"$\langle \sin\theta_\gamma \rangle$")

        ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                      [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$",
                       r"$\frac{3\pi}{2}$", r"$2\pi$"])
        ax.grid(True, linestyle="--", alpha=0.25)

        tmin, tmax = T_EDGES_POS[i], T_EDGES_POS[i+1]
        ax.set_title(rf"$0.10<x_B<0.30,\quad -t\in[{tmin:.2f},{tmax:.2f})$",
                     fontsize=10)

    # Tight layout with a bit more bottom margin for tick labels
    fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])
    outpath = os.path.join(outdir, "sinthetagamma_bins.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")

# ---------------------- main ----------------------
def main():
    if len(sys.argv) != 2:
        print("Usage: python make_sTG_phi_panels.py <input.root>")
        sys.exit(1)

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"[ERROR] File not found: {infile}")
        sys.exit(1)

    outdir = ensure_outdir()
    s_sum, s_cnt = accumulate_means(infile)
    make_panel_plot(s_sum, s_cnt, outdir)

if __name__ == "__main__":
    main()