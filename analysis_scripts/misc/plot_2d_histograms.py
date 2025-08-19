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
    return (data["x"] > xa) & (data["x"] < xb) & (data["fid"] >= 111) & (data["Mx2"] > 0.81) & (data["Mx2"] < 1.00)

def tbin_mask(data, tmin, tmax):
    return (data["tpos"] >= tmin) & (data["tpos"] < tmax)

# ─────────────────────────────────────────────────────────────────────
# Figure 1: 2D histograms (4x3)
# ─────────────────────────────────────────────────────────────────────
def make_2d_canvases(data, outdir):
    pairs = [
        ("Q2","x"),   ("Q2","t"),   ("Q2","phi"),
        ("x","Q2"),   ("x","t"),    ("x","phi"),
        ("t","Q2"),   ("t","x"),    ("t","phi"),
        ("phi","Q2"), ("phi","x"),  ("phi","t"),
    ]
    fig, axes = plt.subplots(4, 3, figsize=(14, 16), constrained_layout=True)

    for ax, (yvar, xvar) in zip(axes.flat, pairs):
        X = extract(data, xvar)
        Y = extract(data, yvar)

        xbins = get_bins(xvar)
        ybins = get_bins(yvar)

        h = ax.hist2d(X, Y, bins=[xbins, ybins], cmap="viridis")
        ax.set_xlim(*({"Q2": Q2_LIM, "x": X_LIM, "t": T_LIM, "phi": PHI_LIM}[xvar]))
        ax.set_ylim(*({"Q2": Q2_LIM, "x": X_LIM, "t": T_LIM, "phi": PHI_LIM}[yvar]))
        ax.set_xlabel(var_label(xvar))
        ax.set_ylabel(var_label(yvar))

        # Nice ticks for phi
        if xvar == "phi":
            ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                          [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
        if yvar == "phi":
            ax.set_yticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                          [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])

    fig.suptitle(r"Global cuts: $Q^{2}>1,\ W>2,\ y<0.75,\ 0.81<M_{x}^{2}<1.00$", fontsize=14, y=0.995)
    outpath = os.path.join(outdir, "2d_histograms.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"Saved: {outpath}")

# ─────────────────────────────────────────────────────────────────────
# Figure 2: Single 4x6 canvas (rows = Q2,x,-t,phi; cols = -t bins)
# ─────────────────────────────────────────────────────────────────────
def make_binned_canvas_4x6(data, outdir):
    # Precompute base masks per xB slice
    base_masks = {name: xb_slice_mask(data, rng) for name, rng in X_SLICES.items()}

    edges = T_EDGES_POS  # len 7 => 6 columns
    nbins = len(edges) - 1

    fig, axes = plt.subplots(4, nbins, figsize=(20, 12), sharey="row")
    fig.subplots_adjust(top=0.88, wspace=0.15, hspace=0.25)
    fig.suptitle(r"Overlaid distributions by $-t$ bin and $x_{B}$ slice"
                 "\nGlobal cuts: $Q^{2}>1,\ W>2,\ y<0.75,\ 0.81<M_{x}^{2}<1.00$",
                 fontsize=14)

    row_vars = ["Q2", "x", "t", "phi"]
    row_lims = {"Q2": Q2_LIM, "x": X_LIM, "t": T_LIM, "phi": PHI_LIM}
    row_bins = {"Q2": NB_Q2,  "x": NB_X,  "t": NB_T,  "phi": NB_PHI}

    # Build a single set of legend handles (one per xB slice)
    legend_handles = []
    legend_labels  = []
    for slice_name, color in XB_COLORS.items():
        h, = plt.plot([], [], color=color, lw=1.8, label=f"{slice_name} "
                       f"[{X_SLICES[slice_name][0]:.2f}, {X_SLICES[slice_name][1]:.2f}]")
        legend_handles.append(h)
        legend_labels.append(h.get_label())
    # Remove dummy line from current axes
    for h in legend_handles: h.remove()

    for col in range(nbins):
        tmin, tmax = edges[col], edges[col+1]
        col_title = rf"$-t \in [{tmin:.2f}, {tmax:.2f})\ \mathrm{{GeV}}^2$"

        for r, vname in enumerate(row_vars):
            ax = axes[r, col]
            vlim = row_lims[vname]
            nb   = row_bins[vname]

            # Axis formatting
            ax.set_xlim(*vlim)
            if col == 0:
                ax.set_ylabel("density")
            ax.set_xlabel(var_label(vname))
            if vname == "phi":
                ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                              [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
            ax.grid(True, linestyle="--", alpha=0.35)

            # Overlays: four xB slices
            for slice_name, color in XB_COLORS.items():
                m = base_masks[slice_name] & tbin_mask(data, tmin, tmax)
                vals = extract(data, vname)[m]
                if vals.size == 0:
                    continue
                ax.hist(vals,
                        bins=np.linspace(vlim[0], vlim[1], nb+1),
                        histtype="step", linewidth=1.6, density=True, color=color)

            # Column title at top row
            if r == 0:
                ax.set_title(col_title, fontsize=12)

    # One figure-level legend explaining slice colors/ranges
    fig.legend(legend_handles, legend_labels, loc="upper center",
               ncol=4, frameon=True, fontsize=11, bbox_to_anchor=(0.5, 0.93))

    outpath = os.path.join(outdir, "binned_distributions.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 2:
        print("Usage: python make_enpi_canvases.py <input.root>")
        sys.exit(1)

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"[ERROR] File not found: {infile}")
        sys.exit(1)

    outdir = ensure_outdir()
    data = load_arrays(infile)

    make_2d_canvases(data, outdir)
    make_binned_canvas_4x6(data, outdir)

if __name__ == "__main__":
    main()