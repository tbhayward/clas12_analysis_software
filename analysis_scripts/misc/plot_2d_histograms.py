#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build two figures from PhysicsEvents in a ROOT file.

Global cuts on ALL plots:
  Q2 > 1, W > 2, y < 0.75, 0.81 < Mx2 < 1.00

1) 2D histograms (4x3 grid), saved to:
   output/enpi+/2d_histograms.pdf
   Grid (Y vs X):
     Q2 vs x,   Q2 vs -t,   Q2 vs phi
     x  vs Q2,  x  vs -t,   x  vs phi
     -t vs Q2,  -t vs x,    -t vs phi
     phi vs Q2, phi vs x,   phi vs -t

2) Single 4x6 canvas of 1D distributions, saved to:
   output/enpi+/binned_distributions.pdf
   Rows:  Q2, x, -t, phi
   Cols:  six (-t) bins from edges in t = [-1.25, -1.05, -0.85, -0.65, -0.45, -0.25, -0.05]
          (converted to -t: [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]; bins are consecutive pairs)
   Each subplot overlays four histograms (Low/MidLow/MidHigh/High x_B slices),
   with per-slice cuts: fiducial_status >= 111, 0.81 < Mx2 < 1.00, and x_B in that slice.
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

# ─────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────
TREE_NAME = "PhysicsEvents"

# Global axis limits
Q2_LIM  = (1.0, 8.0)
X_LIM   = (0.0, 0.70)
T_LIM   = (0.0, 1.30)  # for -t
PHI_LIM = (0.0, 2.0*np.pi)

# Binning
NB_Q2  = 60
NB_X   = 70
NB_T   = 65
NB_PHI = 72

# Labels
LAB_Q2  = r"$Q^{2}\ \mathrm{(GeV^{2})}$"
LAB_X   = r"$x_{B}$"
LAB_T   = r"$-t\ \mathrm{(GeV^{2})}$"
LAB_PHI = r"$\phi$"

# xB slices and colors
X_SLICES = {
    "Low":     (0.10, 0.25),
    "MidLow":  (0.25, 0.35),
    "MidHigh": (0.35, 0.45),
    "High":    (0.45, 0.60),
}
XB_COLORS = {
    "Low":     "tab:blue",
    "MidLow":  "tab:orange",
    "MidHigh": "tab:green",
    "High":    "tab:red",
}

# Provided t-edges (negative t), convert to positive -t edges, ascending
T_EDGES_T   = np.array([-1.25, -1.05, -0.85, -0.65, -0.45, -0.25, -0.05], dtype=float)
T_EDGES_POS = np.sort(-T_EDGES_T)  # [0.05,0.25,0.45,0.65,0.85,1.05,1.25]

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def ensure_outdir():
    outdir = os.path.join("output", "enpi+")
    os.makedirs(outdir, exist_ok=True)
    return outdir

def load_arrays(path):
    """Read branches, compute -t, wrap phi, and apply global cuts."""
    try:
        with uproot.open(path) as f:
            if TREE_NAME not in f:
                print(f"[ERROR] TTree '{TREE_NAME}' not found in {path}")
                sys.exit(2)
            tr = f[TREE_NAME]
            Q2  = tr["Q2"].array(library="np")
            W   = tr["W"].array(library="np")
            y   = tr["y"].array(library="np")
            Mx2 = tr["Mx2"].array(library="np")
            x   = tr["x"].array(library="np")
            t   = tr["t"].array(library="np")
            phi = tr["phi"].array(library="np")
            fid = tr["fiducial_status"].array(library="np")
    except Exception as e:
        print(f"[ERROR] Reading ROOT file failed: {e}")
        sys.exit(3)

    tpos = -t
    phi = np.mod(phi, 2.0*np.pi)

    mask = (Q2 > 1.0) & (W > 2.0) & (y < 0.75) & (Mx2 > 0.81) & (Mx2 < 1.00)
    return {
        "Q2":  Q2[mask],
        "W":   W[mask],
        "y":   y[mask],
        "Mx2": Mx2[mask],
        "x":   x[mask],
        "tpos": tpos[mask],
        "phi": phi[mask],
        "fid": fid[mask],
    }

def get_bins(varname):
    if varname == "Q2":  return np.linspace(Q2_LIM[0],  Q2_LIM[1],  NB_Q2+1)
    if varname == "x":   return np.linspace(X_LIM[0],   X_LIM[1],   NB_X+1)
    if varname == "t":   return np.linspace(T_LIM[0],   T_LIM[1],   NB_T+1)   # -t
    if varname == "phi": return np.linspace(PHI_LIM[0], PHI_LIM[1], NB_PHI+1)
    raise ValueError(varname)

def var_label(varname):
    return {"Q2": LAB_Q2, "x": LAB_X, "t": LAB_T, "phi": LAB_PHI}[varname]

def extract(data, varname):
    return {"Q2": data["Q2"], "x": data["x"], "t": data["tpos"], "phi": data["phi"]}[varname]

def xb_slice_mask(data, xb_range):
    xa, xb = xb_range
    return (data["x"] > xa) & (data["x"] < xb) & (data["fid"] >= 111) & (data["Mx2"] > 0.81) & (data["Mx2"] < 1.