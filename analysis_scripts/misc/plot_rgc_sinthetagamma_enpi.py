#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make a 2x2 canvas of 2D histograms for sin(theta_gamma) vs Q^2, x_B, -t, and phi.

• Input: a ROOT file path passed on the command line.
• Tree:  "PhysicsEvents"
• Branches used: Q2, x, t, sinthetagamma, phi
    - Note: we plot -t (negate the stored 't' branch).
    - Note: phi is in radians; we wrap it into [0, 2*pi).
• Binning:
    - Q^2 in (GeV^2): [1.0, 7.0]        (75 bins)
    - x_B in [0.0, 0.70]                 (75 bins)
    - -t in (GeV^2): [0.0, 1.30]         (75 bins)
    - sin(theta_gamma) in [0.0, 0.50]    (75 bins)
    - phi in (rad): [0.0, 2*pi)          (72 bins; ~5-degree steps)
• Style:
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
BRANCHES = ["Q2", "x", "t", "sinthetagamma", "phi"]

# Axis limits
Q2_LIM  = (1.0, 7.0)
X_LIM   = (0.0, 0.70)
T_LIM   = (0.0, 1.30)          # for -t
SG_LIM  = (0.0, 0.50)          # sin(theta_gamma)
# PHI in radians; set after numpy import to use np.pi
PHI_LIM = (0.0, 2.0 * np.pi)

# 2D bin counts
NB_Q2_2D  = 75
NB_X_2D   = 75
NB_T_2D   = 75
NB_SG     = 75
NB_PHI_2D = 72                 # ~5-degree steps in radians over [0, 2*pi)

# Labels (ASCII-only)
LAB_Q2  = r"$Q^{2}$ (GeV$^{2}$)"
LAB_X   = r"$x_{B}$"
LAB_T   = r"$-t$ (GeV$^{2}$)"
LAB_SG  = r"$\sin\theta_{\gamma}$"
LAB_PHI = r"$\phi$ (rad)"

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
Q2_EDGES  = edges(Q2_LIM,  NB_Q2_2D)
X_EDGES   = edges(X_LIM,   NB_X_2D)
T_EDGES   = edges(T_LIM,   NB_T_2D)
SG_EDGES  = edges(SG_LIM,  NB_SG)
PHI_EDGES = edges(PHI_LIM, NB_PHI_2D)

# ─────────────────────────────────────────────────────────────────────
# Accumulate 2D histograms in a single pass
# ─────────────────────────────────────────────────────────────────────
def accumulate_2d(infile):
    """
    Returns a dict with 4 histograms:
      H_sg_vs_q2, H_sg_vs_x, H_sg_vs_tpos, H_sg_vs_phi
    Each is shaped (len(SG_EDGES)-1, len(X_EDGES)-1) to match pcolormesh (ny, nx).
    """
    H_sg_vs_q2   = np.zeros((len(SG_EDGES)  - 1, len(Q2_EDGES)  - 1), dtype=np.int64)
    H_sg_vs_x    = np.zeros((len(SG_EDGES)  - 1, len(X_EDGES)   - 1), dtype=np.int64)
    H_sg_vs_tpos = np.zeros((len(SG_EDGES)  - 1, len(T_EDGES)   - 1), dtype=np.int64)
    H_sg_vs_phi  = np.zeros((len(SG_EDGES)  - 1, len(PHI_EDGES) - 1), dtype=np.int64)

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
        phi = arrays["phi"]

        # Derived
        tpos = -t
        # Wrap phi to [0, 2*pi)
        phi_mod = np.mod(phi, 2.0 * np.pi)

        # Basic finite masks and range cuts (keep everything within plotting limits)
        m_q2  = np.isfinite(Q2)       & np.isfinite(sg)       & (Q2    >= Q2_LIM[0])   & (Q2    <= Q2_LIM[1])   & (sg >= SG_LIM[0])  & (sg <= SG_LIM[1])
        m_x   = np.isfinite(x)        & np.isfinite(sg)       & (x     >= X_LIM[0])    & (x     <= X_LIM[1])    & (sg >= SG_LIM[0])  & (sg <= SG_LIM[1])
        m_t   = np.isfinite(tpos)     & np.isfinite(sg)       & (tpos  >= T_LIM[0])    & (tpos  <= T_LIM[1])    & (sg >= SG_LIM[0])  & (sg <= SG_LIM[1])
        m_phi = np.isfinite(phi_mod)  & np.isfinite(sg)       & (phi_mod >= PHI_LIM[0])& (phi_mod <  PHI_LIM[1])& (sg >= SG_LIM[0])  & (sg <= SG_LIM[1])

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
        if np.any(m_phi):
            H, xe, ye = np.histogram2d(phi_mod[m_phi], sg[m_phi], bins=[PHI_EDGES, SG_EDGES])
            H_sg_vs_phi += H.T.astype(np.int64)
        #endif
    #endfor

    return {
        "q2":   H_sg_vs_q2,
        "x":    H_sg_vs_x,
        "tpos": H_sg_vs_tpos,
        "phi":  H_sg_vs_phi,
    }
#endfor

# ─────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────
def _draw_panel(ax, H, xedges, yedges, xlim, ylim, xlab, ylab, cmap):
    Hf = H.astype(float)
    Hm = np.ma.masked_where(Hf <= 0.0, Hf)  # zeros → masked → white
    qm = ax.pcolormesh(xedges, yedges, Hm, cmap=cmap, shading="flat")
    qm.set_edgecolor("face")
    qm.set_linewidth(0.0)
    qm.set_rasterized(True)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
#endfor

def plot_2x2(Hs, outdir):
    """
    Create a 2x2 panel of 2D histograms:
      (1) sin(theta_gamma) vs Q^2
      (2) sin(theta_gamma) vs x_B
      (3) sin(theta_gamma) vs -t
      (4) sin(theta_gamma) vs phi (rad)
    """
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.subplots_adjust(wspace=0.22, hspace=0.25, bottom=0.10)

    # Colormap like before (jet, masked zeros -> white)
    cmap = plt.get_cmap("jet").copy()
    cmap.set_bad("white")

    # Arrange panels
    ax11, ax12 = axes[0, 0], axes[0, 1]
    ax21, ax22 = axes[1, 0], axes[1, 1]

    _draw_panel(ax11, Hs["q2"],  Q2_EDGES,  SG_EDGES, Q2_LIM,  SG_LIM, LAB_Q2,  LAB_SG, cmap)
    _draw_panel(ax12, Hs["x"],   X_EDGES,   SG_EDGES, X_LIM,   SG_LIM, LAB_X,   LAB_SG, cmap)
    _draw_panel(ax21, Hs["tpos"],T_EDGES,   SG_EDGES, T_LIM,   SG_LIM, LAB_T,   LAB_SG, cmap)
    _draw_panel(ax22, Hs["phi"], PHI_EDGES, SG_EDGES, PHI_LIM, SG_LIM, LAB_PHI, LAB_SG, cmap)

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
        #endif
        tr = f[TREE_NAME]
        for b in BRANCHES:
            if b not in tr.keys():
                print(f"[ERROR] Branch '{b}' is missing in {TREE_NAME}.")
                sys.exit(3)
            #endif
        #endfor
    #endif

    outdir = ensure_outdir()
    Hs = accumulate_2d(infile)
    plot_2x2(Hs, outdir)
#endfor

if __name__ == "__main__":
    main()
#endif