#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build figures from PhysicsEvents in a ROOT file.

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
   Each subplot overlays four histograms (Low/MidLow/MidHigh/High x_B slices)
   with per-slice cuts: fiducial_status >= 111, 0.81 < Mx2 < 1.00, and x_B in that slice.

3) NEW: Single 4x6 canvas of Mx2 distributions by (x_B slice × -t bin), saved to:
   output/enpi+/Mx2_distributions.pdf
   Rows:  x_B slices = [0.10,0.25], [0.25,0.35], [0.35,0.45], [0.45,0.60]
   Cols:  same six -t bins as above
   For each subplot:
     • plot Mx2 as data points with sqrt(N) error bars over 0.7–1.1
     • fit f(Mx2) = A*exp(-(Mx2-μ)^2/(2σ^2)) + (c0 + c1*Mx2 + c2*Mx2^2)
       with initial μ = m_p^2 (GeV^2); show a thin dashed fit curve
     • legend shows only μ and σ
     • y-axis label: "counts"
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Optional SciPy for non-linear least squares (will fall back to a grid+linear solve if unavailable)
try:
    from scipy.optimize import curve_fit  # type: ignore
    _SCIPY = True
except Exception:
    _SCIPY = False
# endif

# Try uproot for ROOT IO
try:
    import uproot
except ImportError:
    print("[ERROR] Please install uproot:  pip install uproot")
    sys.exit(1)
# endif

# ─────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────
TREE_NAME = "PhysicsEvents"

# Global axis limits
Q2_LIM  = (1.0, 7.0)   # max Q2 set to 7
X_LIM   = (0.0, 0.70)
T_LIM   = (0.0, 1.30)  # for -t
PHI_LIM = (0.0, 2.0*np.pi)

# Mx2 canvas limits/bins
MX2_LIM = (0.7, 1.1)
NB_MX2  = 80

# Binning (base counts used for 1D; 2D uses 2/3 of these)
NB_Q2  = 60
NB_X   = 70
NB_T   = 65
NB_PHI = 72

# 2D binning = 0.6667 × base
NB_Q2_2D  = int(round(NB_Q2  * 2/3))   # 40
NB_X_2D   = int(round(NB_X   * 2/3))   # 47
NB_T_2D   = int(round(NB_T   * 2/3))   # 43
NB_PHI_2D = int(round(NB_PHI * 2/3))   # 48

# Labels
LAB_Q2  = r"$Q^{2}$ (GeV$^{2}$)"
LAB_X   = r"$x_{B}$"
LAB_T   = r"$-t$ (GeV$^{2}$)"
LAB_PHI = r"$\phi$"
LAB_MX2 = r"$M_{x}^{2}$ (GeV$^{2}$)"

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
T_EDGES_POS = np.sort(-T_EDGES_T)  # [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]

# Proton mass (GeV) and squared (GeV^2) for initial guess of μ
M_P  = 0.9382720813
M_P2 = M_P * M_P

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def ensure_outdir():
    outdir = os.path.join("output", "enpi+")
    os.makedirs(outdir, exist_ok=True)
    return outdir
# endfor

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
        # endwith
    except Exception as e:
        print(f"[ERROR] Reading ROOT file failed: {e}")
        sys.exit(3)
    # endif

    tpos = -t
    phi = np.mod(phi, 2.0*np.pi)

    # Global cuts (note: includes 0.81 < Mx2 < 1.00)
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
# endfor

def get_bins_1d(varname):
    if varname == "Q2":  return np.linspace(Q2_LIM[0],  Q2_LIM[1],  NB_Q2+1)
    if varname == "x":   return np.linspace(X_LIM[0],   X_LIM[1],   NB_X+1)
    if varname == "t":   return np.linspace(T_LIM[0],   T_LIM[1],   NB_T+1)   # -t
    if varname == "phi": return np.linspace(PHI_LIM[0], PHI_LIM[1], NB_PHI+1)
    raise ValueError(varname)
# endfor

def get_bins_2d(varname):
    if varname == "Q2":  return np.linspace(Q2_LIM[0],  Q2_LIM[1],  NB_Q2_2D+1)
    if varname == "x":   return np.linspace(X_LIM[0],   X_LIM[1],   NB_X_2D+1)
    if varname == "t":   return np.linspace(T_LIM[0],   T_LIM[1],   NB_T_2D+1)   # -t
    if varname == "phi": return np.linspace(PHI_LIM[0], PHI_LIM[1], NB_PHI_2D+1)
    raise ValueError(varname)
# endfor

def var_label(varname):
    return {"Q2": LAB_Q2, "x": LAB_X, "t": LAB_T, "phi": LAB_PHI}[varname]
# endfor

def extract(data, varname):
    return {"Q2": data["Q2"], "x": data["x"], "t": data["tpos"], "phi": data["phi"]}[varname]
# endfor

def xb_slice_mask(data, xb_range):
    xa, xb = xb_range
    return (data["x"] > xa) & (data["x"] < xb) & (data["fid"] >= 111) & (data["Mx2"] > 0.81) & (data["Mx2"] < 1.00)
# endfor

def xb_slice_mask_for_mx2(data, xb_range):
    """For the Mx2 canvas we keep fiducial and x_B slice; global Mx2 cut remains from load_arrays."""
    xa, xb = xb_range
    return (data["x"] > xa) & (data["x"] < xb) & (data["fid"] >= 111)
# endfor

def tbin_mask(data, tmin, tmax):
    return (data["tpos"] >= tmin) & (data["tpos"] < tmax)
# endfor

# ─────────────────────────────────────────────────────────────────────
# Figure 1: 2D histograms (4x3), with masked zeros → white
# ─────────────────────────────────────────────────────────────────────
def make_2d_canvases(data, outdir):
    pairs = [
        ("Q2","x"),   ("Q2","t"),   ("Q2","phi"),
        ("x","Q2"),   ("x","t"),    ("x","phi"),
        ("t","Q2"),   ("t","x"),    ("t","phi"),
        ("phi","Q2"), ("phi","x"),  ("phi","t"),
    ]

    fig, axes = plt.subplots(4, 3, figsize=(14, 16))
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.subplots_adjust(bottom=0.08)  # extra bottom padding

    # Build colormap: classic jet for data, pure white for empty bins
    cmap = plt.get_cmap("jet").copy()
    cmap.set_bad("white")

    for ax, (yvar, xvar) in zip(axes.flat, pairs):
        X = extract(data, xvar)
        Y = extract(data, yvar)

        xbins = get_bins_2d(xvar)
        ybins = get_bins_2d(yvar)

        # 2D histogram via numpy, then mask zeros and draw with pcolormesh
        H, xedges, yedges = np.histogram2d(X, Y, bins=[xbins, ybins])
        Hm = np.ma.masked_where(H <= 0, H)  # zeros → masked → white

        qm = ax.pcolormesh(xedges, yedges, Hm.T, cmap=cmap, shading="flat")
        qm.set_edgecolor("face")
        qm.set_linewidth(0.0)
        qm.set_rasterized(True)  # avoids hairline gaps in vector backends

        ax.set_xlim(*({"Q2": Q2_LIM, "x": X_LIM, "t": T_LIM, "phi": PHI_LIM}[xvar]))
        ax.set_ylim(*({"Q2": Q2_LIM, "x": X_LIM, "t": T_LIM, "phi": PHI_LIM}[yvar]))
        ax.set_xlabel(var_label(xvar))
        ax.set_ylabel(var_label(yvar))

        # Ticks for phi
        if xvar == "phi":
            ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                          [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
        # endif
        if yvar == "phi":
            ax.set_yticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                          [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
        # endif
    # endfor

    outpath = os.path.join(outdir, "2d_histograms.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
# endfor

# ─────────────────────────────────────────────────────────────────────
# Figure 2: Single 4x6 canvas (rows = Q2,x,-t,phi; cols = -t bins)
# ─────────────────────────────────────────────────────────────────────
def make_binned_canvas_4x6(data, outdir):
    base_masks = {name: xb_slice_mask(data, rng) for name, rng in X_SLICES.items()}

    edges = T_EDGES_POS  # len 7 => 6 columns
    nbins = len(edges) - 1

    fig, axes = plt.subplots(4, nbins, figsize=(20, 12), sharey="row")
    fig.subplots_adjust(top=0.88, wspace=0.15, hspace=0.25)
    fig.suptitle(
        r"Overlaid distributions by $-t$ bin and $x_{B}$ slice" "\n"
        r"Global cuts: $Q^{2}>1,\, W>2,\, y<0.75,\, 0.81<M_{x}^{2}<1.00$",
        fontsize=14
    )

    row_vars = ["Q2", "x", "t", "phi"]
    row_lims = {"Q2": Q2_LIM, "x": X_LIM, "t": T_LIM, "phi": PHI_LIM}
    row_bins = {"Q2": NB_Q2,  "x": NB_X,  "t": NB_T,  "phi": NB_PHI}

    # Legend handles with requested labels (x_B intervals)
    legend_handles, legend_labels = [], []
    for slice_name, color in XB_COLORS.items():
        a, b = X_SLICES[slice_name]
        lbl = rf"$x_{{B}}\in[{a:.2f},{b:.2f}]$"
        h, = plt.plot([], [], color=color, lw=1.8, label=lbl)
        legend_handles.append(h); legend_labels.append(lbl)
    # endfor
    for h in legend_handles: h.remove()

    for col in range(nbins):
        tmin, tmax = edges[col], edges[col+1]
        col_title = rf"$-t \in [{tmin:.2f}, {tmax:.2f})$ (GeV$^{{2}}$)"

        for r, vname in enumerate(row_vars):
            ax = axes[r, col]
            vlim = row_lims[vname]
            nb   = row_bins[vname]

            ax.set_xlim(*vlim)
            if col == 0:
                ax.set_ylabel("counts")
            # endif
            ax.set_xlabel(var_label(vname))
            if vname == "phi":
                ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                              [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
            # endif
            ax.grid(True, linestyle="--", alpha=0.25)

            for slice_name, color in XB_COLORS.items():
                m = base_masks[slice_name] & tbin_mask(data, tmin, tmax)
                vals = extract(data, vname)[m]
                if vals.size == 0:
                    continue
                # endif
                ax.hist(vals,
                        bins=np.linspace(vlim[0], vlim[1], nb+1),
                        histtype="step", linewidth=1.6, density=False, color=color)
            # endfor

            if r == 0:
                ax.set_title(col_title, fontsize=12)
            # endif
        # endfor
    # endfor

    fig.legend(legend_handles, legend_labels, loc="upper center",
               ncol=4, frameon=True, fontsize=11, bbox_to_anchor=(0.5, 0.93))

    outpath = os.path.join(outdir, "binned_distributions.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
# endfor

# ─────────────────────────────────────────────────────────────────────
# Figure 3: NEW — Mx2 distributions by (x_B slice × -t bin), with fits
# ─────────────────────────────────────────────────────────────────────
def _gauss_plus_quad(x, A, mu, sigma, c0, c1, c2):
    return A * np.exp(-0.5 * ((x - mu) / np.maximum(sigma, 1e-6))**2) + (c0 + c1 * x + c2 * x * x)
# endfor

def _fit_gauss_plus_quad(xc, y, yerr, mu0=M_P2, sigma0=0.03):
    """
    Fit y(x) with A*exp(-(x-mu)^2/(2*sigma^2)) + (c0 + c1*x + c2*x^2).
    Returns (A, mu, sigma, c0, c1, c2), yfit.
    Uses SciPy if available; otherwise a coarse-to-fine grid over (mu, sigma)
    with weighted linear least squares for (A,c0,c1,c2).
    """
    # Clean inputs
    sel = np.isfinite(xc) & np.isfinite(y)
    xc = xc[sel]; y = y[sel]; yerr = yerr[sel]
    yerr = np.where(yerr > 0.0, yerr, 1.0)

    # If too few points, bail gracefully
    if xc.size < 6:
        A = max(y.max(), 1.0)
        p = (A, mu0, max(sigma0, 1e-3), 0.0, 0.0, 0.0)
        return p, _gauss_plus_quad(xc, *p)
    # endif

    if _SCIPY:
        # Initial guesses
        A0 = max(y.max() - np.median(y), 1.0)
        p0 = [A0, mu0, sigma0, float(np.min(y)), 0.0, 0.0]
        # Reasonable bounds
        lower = [0.0,  0.80, 0.005, -np.inf, -np.inf, -np.inf]
        upper = [np.inf, 0.95, 0.080,  np.inf,  np.inf,  np.inf]
        try:
            popt, _ = curve_fit(_gauss_plus_quad, xc, y, p0=p0, sigma=yerr,
                                absolute_sigma=True, maxfev=20000, bounds=(lower, upper))
            yfit = _gauss_plus_quad(xc, *popt)
            return tuple(popt), yfit
        except Exception:
            # fall back to grid below
            pass
        # endif
    # endif

    # Fallback: grid search in (mu, sigma) with weighted linear solve for A,c0,c1,c2
    def solve_for(mu, sigma):
        G0 = np.exp(-0.5 * ((xc - mu) / max(sigma, 1e-6))**2)
        G = np.column_stack([G0, np.ones_like(xc), xc, xc * xc])
        w = 1.0 / yerr
        Gw = G * w[:, None]
        yw = y * w
        coeffs, *_ = np.linalg.lstsq(Gw, yw, rcond=None)
        yhat = G @ coeffs
        chi2 = np.sum(((y - yhat) / yerr)**2)
        return coeffs, yhat, chi2
    # endfor

    grid_mu = np.linspace(mu0 - 0.05, mu0 + 0.05, 21)
    grid_sg = np.linspace(0.010, 0.060, 12)

    best = (np.inf, None, None, None, None)
    for mu in grid_mu:
        for sg in grid_sg:
            coeffs, yhat, chi2 = solve_for(mu, sg)
            if chi2 < best[0]:
                best = (chi2, mu, sg, coeffs, yhat)
            # endif
        # endfor
    # endfor

    # Refine around the best
    _, mu_b, sg_b, _, _ = best
    grid_mu2 = np.linspace(mu_b - 0.02, mu_b + 0.02, 17)
    grid_sg2 = np.linspace(max(0.008, sg_b - 0.02), sg_b + 0.02, 17)

    for mu in grid_mu2:
        for sg in grid_sg2:
            coeffs, yhat, chi2 = solve_for(mu, sg)
            if chi2 < best[0]:
                best = (chi2, mu, sg, coeffs, yhat)
            # endif
        # endfor
    # endfor

    _, mu_b, sg_b, coeffs_b, yhat_b = best
    A_b, c0_b, c1_b, c2_b = coeffs_b
    params = (float(A_b), float(mu_b), float(max(sg_b, 1e-3)), float(c0_b), float(c1_b), float(c2_b))
    return params, yhat_b
# endfor

def make_mx2_canvas(data, outdir):
    """
    Build a 4x6 canvas of Mx2 distributions:
      rows = x_B slices (Low, MidLow, MidHigh, High)
      cols = -t bins based on T_EDGES_POS
    Each subplot:
      - points with sqrt(N) error bars over 0.7–1.1
      - dashed fit (Gaussian + quadratic), legend shows μ and σ only
    """
    edges = T_EDGES_POS  # length 7 => 6 columns
    nbins = len(edges) - 1
    slice_order = ["Low", "MidLow", "MidHigh", "High"]

    fig, axes = plt.subplots(len(slice_order), nbins, figsize=(20, 12), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.90, wspace=0.15, hspace=0.25)
    fig.suptitle(
        r"$M_{x}^{2}$ distributions by $x_{B}$ slice and $-t$ bin"
        "\n" + r"Global cuts: $Q^{2}>1,\, W>2,\, y<0.75,\, 0.81<M_{x}^{2}<1.00$",
        fontsize=14
    )

    mx2_bins = np.linspace(MX2_LIM[0], MX2_LIM[1], NB_MX2 + 1)
    mx2_centers = 0.5 * (mx2_bins[:-1] + mx2_bins[1:])

    # Loop over rows (x_B slices) and columns (-t bins)
    for r, sname in enumerate(slice_order):
        xb_range = X_SLICES[sname]
        base_mask = xb_slice_mask_for_mx2(data, xb_range)
        for c in range(nbins):
            tmin, tmax = edges[c], edges[c+1]
            ax = axes[r, c]

            # mask for this cell
            mask = base_mask & tbin_mask(data, tmin, tmax)
            vals = data["Mx2"][mask]

            # Histogram -> counts and sqrt(N) errors
            counts, _ = np.histogram(vals, bins=mx2_bins)
            y = counts.astype(float)
            yerr = np.sqrt(np.maximum(y, 1.0))  # avoid zero errors

            # Plot as points with error bars
            ax.errorbar(mx2_centers, y, yerr=yerr, fmt='o', ms=3, lw=1.0, capsize=2,
                        color='black', label=None)

            # Fit Gaussian + quadratic if we have enough stats
            label_txt = r"$\mu=\mathrm{N/A},\ \sigma=\mathrm{N/A}$"
            if np.count_nonzero(y) >= 6:
                params, yfit = _fit_gauss_plus_quad(mx2_centers, y, yerr, mu0=M_P2, sigma0=0.03)
                A, mu, sigma, c0, c1, c2 = params
                # Plot fit as thin dashed line
                ax.plot(mx2_centers, yfit, linestyle='--', linewidth=1.0, color='black',
                        label=rf"$\mu={mu:.4f},\ \sigma={sigma:.4f}$")
            else:
                # still add a dummy handle so legend shows N/A
                ax.plot([], [], linestyle='--', linewidth=1.0, color='black', label=label_txt)
            # endif

            # Styling
            ax.set_xlim(*MX2_LIM)
            ax.set_ylim(bottom=0)
            ax.set_xlabel(LAB_MX2)
            ax.set_ylabel("counts")
            ax.grid(True, linestyle="--", alpha=0.25)

            # Title with x_B and -t ranges (units on -t)
            xa, xb = xb_range
            ax.set_title(
                rf"$x_{{B}}\in[{xa:.2f},{xb:.2f}],\ -t\in[{tmin:.2f},{tmax:.2f})$ (GeV$^{{2}}$)",
                fontsize=10
            )

            # Legend with only mu and sigma
            ax.legend(frameon=True, fontsize=8, loc="best")
        # endfor
    # endfor

    outpath = os.path.join(outdir, "Mx2_distributions.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
# endfor

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 2:
        print("Usage: python make_enpi_canvases.py <input.root>")
        sys.exit(1)
    # endif

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"[ERROR] File not found: {infile}")
        sys.exit(1)
    # endif

    outdir = ensure_outdir()
    data = load_arrays(infile)

    make_2d_canvases(data, outdir)
    make_binned_canvas_4x6(data, outdir)
    make_mx2_canvas(data, outdir)
# endfor

if __name__ == "__main__":
    main()
# endif