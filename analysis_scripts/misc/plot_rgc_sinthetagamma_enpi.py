#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make a 1×3 canvas of 2D histograms for sin(theta_gamma) vs Q^2, x_B, and -t.

• Input: a ROOT file path passed on the command line.
• Tree:  "PhysicsEvents"
• Branches used: Q2, x, t, sinthetagamma
    - Note: we plot -t (negate the stored 't' branch).
• Binning:
    - Q^2 in [1.0, 7.0]          (20 bins)
    - x_B in [0.0, 0.70]         (23 bins)
    - -t  in [0.0, 1.30]         (21 bins)
    - sin(theta_gamma) in [0,0.5] (40 bins)
• Style matches the earlier 2D plots:
    - pcolormesh with the "jet" colormap
    - zero-count bins masked to white (cmap.set_bad("white"))
    - no colorbar; vector-friendly (rasterized) pcolormesh
• Output: output/enpi+/sinthetagamma_2D.pdf
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Uproot for ROOT IO
try:
    import uproot
except ImportError:
    print("[ERROR] Please install uproot:  pip install uproot")
    sys.exit(1)
#endif

# ─────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────
TREE_NAME = "PhysicsEvents"
BRANCHES = ["Q2", "x", "t", "sinthetagamma"]

# Axis limits (mirroring prior scripts)
Q2_LIM  = (1.0, 7.0)
X_LIM   = (0.0, 0.70)
T_LIM   = (0.0, 1.30)     # for -t
SG_LIM  = (0.0, 0.50)     # sin(theta_gamma)

# 2D bin counts (same as earlier for X axes; choose 40 bins for sinθγ)
NB_Q2_2D = 20
NB_X_2D  = 23
NB_T_2D  = 21
NB_SG    = 40

# Labels (units on axes only)
LAB_Q2 = r"$Q^{2}$ (GeV$^{2}$)"
LAB_X  = r"$x_{B}$"
LAB_T  = r"$-t$ (GeV$^{2}$)"
LAB_SG = r"$\sin\theta_{\gamma}$"

# Uproot iteration step
ITER_STEP = "200 MB"

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def ensure_outdir():
    outdir = os.path.join("output", "enpi+")
    os.makedirs(outdir, exist_ok=True)
    return outdir
#endfor

def edges(lim, nb):
    return np.linspace(lim[0], lim[1], nb + 1)
#endfor

# Precompute bin edges
Q2_EDGES = edges(Q2_LIM, NB_Q2_2D)
X_EDGES  = edges(X_LIM,  NB_X_2D)
T_EDGES  = edges(T_LIM,  NB_T_2D)
SG_EDGES = edges(SG_LIM, NB_SG)

# ─────────────────────────────────────────────────────────────────────
# Accumulate 2D histograms in a single pass
# ─────────────────────────────────────────────────────────────────────
def accumulate_2d(infile):
    """
    Returns a dict with 3 histograms:
      H_sg_vs_q2, H_sg_vs_x, H_sg_vs_tpos
    Each is shaped (len(SG_EDGES)-1, len(X_EDGES)-1) to match pcolormesh (ny, nx).
    """
    H_sg_vs_q2   = np.zeros((len(SG_EDGES) - 1, len(Q2_EDGES) - 1), dtype=np.int64)
    H_sg_vs_x    = np.zeros((len(SG_EDGES) - 1, len(X_EDGES)  - 1), dtype=np.int64)
    H_sg_vs_tpos = np.zeros((len(SG_EDGES) - 1, len(T_EDGES)  - 1), dtype=np.int64)

    it = uproot.iterate(
        {infile: TREE_NAME},
        BRANCHES,
        step_size=ITER_STEP,
        library="np",
        allow_read_errors=False
    )

    for arrays in it:
        # Extract and prepare variables
        Q2  = arrays["Q2"]
        x   = arrays["x"]
        t   = arrays["t"]
        sg  = arrays["sinthetagamma"]

        # Derived
        tpos = -t

        # Basic finite masks and range cuts (keep everything within our plotting limits)
        m_q2  = np.isfinite(Q2) & np.isfinite(sg)  & (Q2  >= Q2_LIM[0]) & (Q2  <= Q2_LIM[1]) & (sg >= SG_LIM[0]) & (sg <= SG_LIM[1])
        m_x   = np.isfinite(x)  & np.isfinite(sg)  & (x   >= X_LIM[0])  & (x   <= X_LIM[1])  & (sg >= SG_LIM[0]) & (sg <= SG_LIM[1])
        m_t   = np.isfinite(tpos) & np.isfinite(sg) & (tpos >= T_LIM[0]) & (tpos <= T_LIM[1]) & (sg >= SG_LIM[0]) & (sg <= SG_LIM[1])

        if np.any(m_q2):
            H, xe, ye = np.histogram2d(Q2[m_q2], sg[m_q2], bins=[Q2_EDGES, SG_EDGES])
            H_sg_vs_q2 += H.T.astype(np.int64)   # transpose so shape is (ny, nx)
        #endif
        if np.any(m_x):
            H, xe, ye = np.histogram2d(x[m_x], sg[m_x], bins=[X_EDGES, SG_EDGES])
            H_sg_vs_x += H.T.astype(np.int64)
        #endif
        if np.any(m_t):
            H, xe, ye = np.histogram2d(tpos[m_t], sg[m_t], bins=[T_EDGES, SG_EDGES])
            H_sg_vs_tpos += H.T.astype(np.int64)
        #endif
    #endfor

    return {
        "q2":   H_sg_vs_q2,
        "x":    H_sg_vs_x,
        "tpos": H_sg_vs_tpos,
    }
#endfor

# ─────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────
def plot_1x3(Hs, outdir):
    """
    Create a 1×3 panel of 2D histograms:
      (left)  sinθγ vs Q^2
      (mid)   sinθγ vs x_B
      (right) sinθγ vs -t
    Styling matches the earlier 2D canvas (jet, zeros→white, pcolormesh).
    """
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.subplots_adjust(wspace=0.18, bottom=0.14)

    # Colormap like before (jet, masked zeros -> white)
    cmap = plt.get_cmap("jet").copy()
    cmap.set_bad("white")

    panels = [
        ("q2",   axes[0], Q2_EDGES, SG_EDGES, LAB_Q2, LAB_SG, Q2_LIM, SG_LIM),
        ("x",    axes[1], X_EDGES,  SG_EDGES, LAB_X,  LAB_SG, X_LIM,  SG_LIM),
        ("tpos", axes[2], T_EDGES,  SG_EDGES, LAB_T,  LAB_SG, T_LIM,  SG_LIM),
    ]

    for key, ax, xedges, yedges, xlab, ylab, xlim, ylim in panels:
        H = Hs[key].astype(float)
        Hm = np.ma.masked_where(H <= 0.0, H)  # zeros → masked → white

        qm = ax.pcolormesh(xedges, yedges, Hm, cmap=cmap, shading="flat")
        qm.set_edgecolor("face")
        qm.set_linewidth(0.0)
        qm.set_rasterized(True)  # vector backends friendly

        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    #endfor

    outpath = os.path.join(outdir, "sinthetagamma_2D.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
#endfor

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 2:
        print("Usage: python make_sinthetagamma_2D.py <input.root>")
        sys.exit(1)
    #endif

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"[ERROR] File not found: {infile}")
        sys.exit(1)
    #endif

    # Quick sanity: verify branch presence up front
    with uproot.open(infile) as f:
        if TREE_NAME not in f:
            print(f"[ERROR] TTree '{TREE_NAME}' not found in {infile}")
            sys.exit(2)
        tr = f[TREE_NAME]
        for b in BRANCHES:
            if b not in tr.keys():
                print(f"[ERROR] Branch '{b}' is missing in {TREE_NAME}.")
                sys.exit(3)
        #endif
    #endif

    outdir = ensure_outdir()
    Hs = accumulate_2d(infile)
    plot_1x3(Hs, outdir)
#endfor

if __name__ == "__main__":
    main()
#endif