#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build figures from PhysicsEvents in a ROOT file with a single-pass accumulator.

Global cuts for most plots:
  Q2 > 1, W > 2, y < 0.75, 0.81 < Mx2 < 1.00
Exception:
  The Mx2 canvas uses NO Mx2 cut (only Q2/W/y + fiducial + x_B + -t).

Outputs:
  • output/enpi+/2d_histograms.pdf
  • output/enpi+/binned_distributions.pdf
  • output/enpi+/Mx2_distributions.pdf
  • output/enpi+/Mx2_mu_vs_bin.pdf
  • output/enpi+/Mx2_fit_params.csv   <-- (μ, σ) saved for each (xB, -t) bin
"""

import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Optional SciPy for robust non-linear fits; falls back to grid search.
try:
    from scipy.optimize import curve_fit  # type: ignore
    _SCIPY = True
except Exception:
    _SCIPY = False
#endif

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
BRANCHES = ["Q2","W","y","Mx2","x","t","phi","fiducial_status"]

# Axis limits
Q2_LIM  = (1.0, 7.0)
X_LIM   = (0.0, 0.70)
T_LIM   = (0.0, 1.30)  # for -t
PHI_LIM = (0.0, 2.0*np.pi)
MX2_LIM = (0.7, 1.1)

# 1D base bin counts (HALVED vs original; 2D gets 2/3 of these)
NB_Q2  = 67
NB_X   = 67
NB_T   = 67
NB_PHI = 67
NB_MX2 = 67

# 2D binning = 0.6667 × base (rounded)
NB_Q2_2D  = int(round(NB_Q2  ))   # 20
NB_X_2D   = int(round(NB_X   ))   # 23
NB_T_2D   = int(round(NB_T   ))   # 21
NB_PHI_2D = int(round(NB_PHI ))   # 24

# Labels (units in parentheses for axes only)
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
SLICE_ORDER = ["Low","MidLow","MidHigh","High"]
XB_COLORS = {"Low":"tab:blue","MidLow":"tab:orange","MidHigh":"tab:green","High":"tab:red"}

# Provided t-edges (negative t) -> positive -t edges ascending
T_EDGES_T   = np.array([-1.25, -1.05, -0.85, -0.65, -0.45, -0.25, -0.05], dtype=float)
T_EDGES_POS = np.sort(-T_EDGES_T)  # [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]
NTCOLS = len(T_EDGES_POS) - 1      # 6

# Proton mass for initial mu (GeV^2)
M_P  = 0.9382720813
M_P2 = M_P * M_P

# Uproot iteration step size (tune for your machine)
ITER_STEP = "200 MB"

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def ensure_outdir():
    outdir = os.path.join("output", "enpi+")
    os.makedirs(outdir, exist_ok=True)
    return outdir
#endfor

def edges_and_centers(lim, nb):
    edges = np.linspace(lim[0], lim[1], nb+1)
    centers = 0.5*(edges[:-1] + edges[1:])
    return edges, centers
#endfor

# Bin edges for all axes
Q2_EDGES, Q2_CEN   = edges_and_centers(Q2_LIM,  NB_Q2)
X_EDGES,  X_CEN    = edges_and_centers(X_LIM,   NB_X)
T_EDGES,  T_CEN    = edges_and_centers(T_LIM,   NB_T)
PHI_EDGES, PHI_CEN = edges_and_centers(PHI_LIM, NB_PHI)
MX2_EDGES, MX2_CEN = edges_and_centers(MX2_LIM, NB_MX2)

# 2D edges
Q2_EDGES_2D  = np.linspace(Q2_LIM[0],  Q2_LIM[1],  NB_Q2_2D+1)
X_EDGES_2D   = np.linspace(X_LIM[0],   X_LIM[1],   NB_X_2D+1)
T_EDGES_2D   = np.linspace(T_LIM[0],   T_LIM[1],   NB_T_2D+1)
PHI_EDGES_2D = np.linspace(PHI_LIM[0], PHI_LIM[1], NB_PHI_2D+1)

# Map from name to edges (2D and 1D)
EDGES_2D = {"Q2":Q2_EDGES_2D, "x":X_EDGES_2D, "t":T_EDGES_2D, "phi":PHI_EDGES_2D}
EDGES_1D = {"Q2":Q2_EDGES, "x":X_EDGES, "t":T_EDGES, "phi":PHI_EDGES}

# 2D pairs for Figure 1
PAIRS_2D = [
    ("Q2","x"),   ("Q2","t"),   ("Q2","phi"),
    ("x","Q2"),   ("x","t"),    ("x","phi"),
    ("t","Q2"),   ("t","x"),    ("t","phi"),
    ("phi","Q2"), ("phi","x"),  ("phi","t"),
]

def var_from_arrays(arrs, name):
    if name == "t":   return -arrs["t"]  # -t
    if name == "phi": return np.mod(arrs["phi"], 2.0*np.pi)
    return arrs[name]
#endfor

def slice_index_from_x(x):
    # Returns slice index in {0..3} or -1 if out of all ranges
    idx = np.full(x.shape, -1, dtype=np.int16)
    for k,(xa,xb) in enumerate([(0.10,0.25),(0.25,0.35),(0.35,0.45),(0.45,0.60)]):
        m = (x > xa) & (x < xb)
        idx[m] = k
    #endfor
    return idx
#endfor

def tbin_index_from_tpos(tpos):
    # 6 bins from T_EDGES_POS; interval [lo, hi)
    idx = np.digitize(tpos, T_EDGES_POS) - 1
    bad = (idx < 0) | (idx >= NTCOLS)
    idx[bad] = -1
    return idx.astype(np.int16)
#endfor

def _gauss_plus_quad(x, A, mu, sigma, c0, c1, c2):
    return A * np.exp(-0.5 * ((x - mu) / np.maximum(sigma, 1e-6))**2) + (c0 + c1 * x + c2 * x * x)
#endfor

def _fit_gauss_plus_quad(xc, y, yerr, mu0=M_P2, sigma0=0.03):
    sel = np.isfinite(xc) & np.isfinite(y)
    xc = xc[sel]; y = y[sel]; yerr = yerr[sel]
    yerr = np.where(yerr > 0.0, yerr, 1.0)

    if xc.size < 6 or np.sum(y) < 1:
        A = max(y.max() if y.size else 1.0, 1.0)
        p = (A, mu0, max(sigma0, 1e-3), 0.0, 0.0, 0.0)
        return p, _gauss_plus_quad(xc, *p), False
    #endif

    if _SCIPY:
        A0 = max((y.max() - np.median(y)), 1.0)
        p0 = [A0, mu0, sigma0, float(np.min(y)), 0.0, 0.0]
        lower = [0.0,  0.80, 0.005, -np.inf, -np.inf, -np.inf]
        upper = [np.inf, 0.95, 0.080,  np.inf,  np.inf,  np.inf]
        try:
            popt, _ = curve_fit(_gauss_plus_quad, xc, y, p0=p0, sigma=yerr,
                                absolute_sigma=True, maxfev=20000, bounds=(lower, upper))
            yfit = _gauss_plus_quad(xc, *popt)
            return tuple(popt), yfit, True
        except Exception:
            pass
        #endif
    #endif

    # Fallback grid + weighted linear solve for (A,c0,c1,c2)
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
    #endfor

    grid_mu = np.linspace(mu0 - 0.05, mu0 + 0.05, 21)
    grid_sg = np.linspace(0.010, 0.060, 12)

    best = (np.inf, None, None, None, None)
    for mu in grid_mu:
        for sg in grid_sg:
            coeffs, yhat, chi2 = solve_for(mu, sg)
            if chi2 < best[0]:
                best = (chi2, mu, sg, coeffs, yhat)
            #endif
        #endfor
    #endfor

    _, mu_b, sg_b, _, _ = best
    grid_mu2 = np.linspace(mu_b - 0.02, mu_b + 0.02, 17)
    grid_sg2 = np.linspace(max(0.008, sg_b - 0.02), sg_b + 0.02, 17)

    for mu in grid_mu2:
        for sg in grid_sg2:
            coeffs, yhat, chi2 = solve_for(mu, sg)
            if chi2 < best[0]:
                best = (chi2, mu, sg, coeffs, yhat)
            #endif
        #endfor
    #endfor

    _, mu_b, sg_b, coeffs_b, yhat_b = best
    A_b, c0_b, c1_b, c2_b = coeffs_b
    params = (float(A_b), float(mu_b), float(max(sg_b, 1e-3)), float(c0_b), float(c1_b), float(c2_b))
    return params, yhat_b, True
#endfor

# ─────────────────────────────────────────────────────────────────────
# Single-pass accumulators
# ─────────────────────────────────────────────────────────────────────
def accumulate(infile):
    """
    Iterate over the tree ONCE (chunked) and fill:
      • 12× 2D histograms (with Mx2 exclusivity)
      • 4×6× {Q2,x,t,phi} 1D hist counts (with Mx2 exclusivity and fid>=111)
      • 4×6× Mx2 hist counts (NO Mx2 exclusivity, but fid>=111)
    """
    # 2D hist containers
    H2 = {}
    for (yvar, xvar) in PAIRS_2D:
        H2[(yvar, xvar)] = np.zeros((len(EDGES_2D[yvar])-1, len(EDGES_2D[xvar])-1), dtype=np.int64)
    #endfor

    # 1D overlaid distributions (with exclusivity and fid)
    H1_Q2  = np.zeros((4, NTCOLS, NB_Q2 ), dtype=np.int64)
    H1_X   = np.zeros((4, NTCOLS, NB_X  ), dtype=np.int64)
    H1_T   = np.zeros((4, NTCOLS, NB_T  ), dtype=np.int64)
    H1_PHI = np.zeros((4, NTCOLS, NB_PHI), dtype=np.int64)

    # Mx2 per cell (NO exclusivity, but fid>=111)
    H1_MX2 = np.zeros((4, NTCOLS, NB_MX2), dtype=np.int64)

    # Iterate once over the file
    it = uproot.iterate(
        {infile: TREE_NAME},
        BRANCHES,
        step_size=ITER_STEP,
        library="np",
        allow_read_errors=False
    )

    for arrays in it:
        Q2  = arrays["Q2"]
        W   = arrays["W"]
        y   = arrays["y"]
        Mx2 = arrays["Mx2"]
        x   = arrays["x"]
        t   = arrays["t"]
        phi = arrays["phi"]
        fid = arrays["fiducial_status"]

        # Precompute derived vars
        tpos = -t
        phi_wrap = np.mod(phi, 2.0*np.pi)

        # Base masks
        base = (Q2 > 1.0) & (W > 2.0) & (y < 0.75)
        mx2_win = (Mx2 > 0.81) & (Mx2 < 1.00)
        with_mx2  = base & mx2_win
        no_mx2    = base
        fid_good  = fid >= 111

        # -------------------- Figure 1: 2D hist pairs (with exclusivity) --------------------
        if np.any(with_mx2):
            arrs_w = {"Q2":Q2[with_mx2], "x":x[with_mx2], "t":t[with_mx2], "phi":phi_wrap[with_mx2]}
            for (yvar, xvar) in PAIRS_2D:
                X = var_from_arrays(arrs_w, xvar)
                Y = var_from_arrays(arrs_w, yvar)
                H, xe, ye = np.histogram2d(X, Y, bins=[EDGES_2D[xvar], EDGES_2D[yvar]])
                H2[(yvar, xvar)] += H.T.astype(np.int64)
            #endfor
        #endif

        # -------------------- Figure 2: 4x6 overlain 1D (with exclusivity + fid) --------------------
        mask_12 = with_mx2 & fid_good
        if np.any(mask_12):
            x_w   = x[mask_12]
            tpos_w = tpos[mask_12]
            phi_w  = phi_wrap[mask_12]
            Q2_w   = Q2[mask_12]

            sidx = slice_index_from_x(x_w)           # 0..3 or -1
            tbix = tbin_index_from_tpos(tpos_w)      # 0..5 or -1
            valid = (sidx >= 0) & (tbix >= 0)
            if np.any(valid):
                sidx = sidx[valid]
                tbix = tbix[valid]

                # Bin indices for each var
                q2b = np.digitize(Q2_w[valid],  Q2_EDGES)  - 1
                xb  = np.digitize(x_w[valid],   X_EDGES)   - 1
                tb  = np.digitize(tpos_w[valid],T_EDGES)   - 1
                pb  = np.digitize(phi_w[valid], PHI_EDGES) - 1

                # Clamp to valid ranges
                def clamp(b, n):
                    b[(b < 0) | (b >= n)] = -1
                    return b
                #endfor
                q2b = clamp(q2b, NB_Q2)
                xb  = clamp(xb,  NB_X)
                tb  = clamp(tb,  NB_T)
                pb  = clamp(pb,  NB_PHI)

                good_q2 = q2b >= 0
                good_x  = xb  >= 0
                good_t  = tb  >= 0
                good_p  = pb  >= 0

                np.add.at(H1_Q2,  (sidx[good_q2], tbix[good_q2], q2b[good_q2]), 1)
                np.add.at(H1_X,   (sidx[good_x ], tbix[good_x ], xb [good_x ]), 1)
                np.add.at(H1_T,   (sidx[good_t ], tbix[good_t ], tb [good_t ]), 1)
                np.add.at(H1_PHI, (sidx[good_p ], tbix[good_p ], pb [good_p ]), 1)
            #endif
        #endif

        # -------------------- Figure 3: Mx2 per cell (NO exclusivity, fid required) --------------------
        mask_mx2 = no_mx2 & fid_good
        if np.any(mask_mx2):
            x_n   = x[mask_mx2]
            tpos_n = tpos[mask_mx2]
            Mx2_n  = Mx2[mask_mx2]

            sidx = slice_index_from_x(x_n)           # 0..3 or -1
            tbix = tbin_index_from_tpos(tpos_n)      # 0..5 or -1
            valid = (sidx >= 0) & (tbix >= 0)
            if np.any(valid):
                sidx = sidx[valid]
                tbix = tbix[valid]
                mxb  = np.digitize(Mx2_n[valid], MX2_EDGES) - 1
                mxb[(mxb < 0) | (mxb >= NB_MX2)] = -1
                good = mxb >= 0
                np.add.at(H1_MX2, (sidx[good], tbix[good], mxb[good]), 1)
            #endif
        #endif
    #endfor

    return H2, H1_Q2, H1_X, H1_T, H1_PHI, H1_MX2
#endfor

# ─────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────
def make_2d_canvases(H2, outdir):
    fig, axes = plt.subplots(4, 3, figsize=(14, 16))
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.subplots_adjust(bottom=0.08)

    cmap = plt.get_cmap("jet").copy()
    cmap.set_bad("white")

    for ax, (yvar, xvar) in zip(axes.flat, PAIRS_2D):
        H = H2[(yvar, xvar)].astype(float)
        Hm = np.ma.masked_where(H <= 0.0, H)

        xedges = EDGES_2D[xvar]
        yedges = EDGES_2D[yvar]

        qm = ax.pcolormesh(xedges, yedges, Hm, cmap=cmap, shading="flat")
        qm.set_edgecolor("face"); qm.set_linewidth(0.0); qm.set_rasterized(True)

        lims = {"Q2":Q2_LIM,"x":X_LIM,"t":T_LIM,"phi":PHI_LIM}
        ax.set_xlim(*lims[xvar]); ax.set_ylim(*lims[yvar])
        ax.set_xlabel({"Q2":LAB_Q2,"x":LAB_X,"t":LAB_T,"phi":LAB_PHI}[xvar])
        ax.set_ylabel({"Q2":LAB_Q2,"x":LAB_X,"t":LAB_T,"phi":LAB_PHI}[yvar])

        if xvar == "phi":
            ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                          [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
        #endif
        if yvar == "phi":
            ax.set_yticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                          [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
        #endif
    #endfor

    outpath = os.path.join(outdir, "2d_histograms.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
#endfor

def make_binned_canvas_4x6(H1_Q2, H1_X, H1_T, H1_PHI, outdir):
    row_vars = [("Q2", H1_Q2,  Q2_EDGES,  Q2_CEN,  LAB_Q2,  Q2_LIM),
                ("x",  H1_X,   X_EDGES,   X_CEN,   LAB_X,   X_LIM),
                ("t",  H1_T,   T_EDGES,   T_CEN,   LAB_T,   T_LIM),
                ("phi",H1_PHI, PHI_EDGES, PHI_CEN, LAB_PHI, PHI_LIM)]

    fig, axes = plt.subplots(4, NTCOLS, figsize=(20, 12), sharey="row")
    # add a bit of extra room below the bottom row and slightly more row spacing
    fig.subplots_adjust(top=0.88, bottom=0.095, wspace=0.12, hspace=0.24)
    fig.suptitle(
        r"Overlaid distributions by $-t$ bin and $x_{B}$ slice" "\n"
        r"Global cuts: $Q^{2}>1,\, W>2,\, y<0.75,\, 0.81<M_{x}^{2}<1.00$",
        fontsize=14
    )

    # Prebuild legend handles (for a clean, single legend)
    legend_handles = []
    legend_labels  = []
    for s in SLICE_ORDER:
        lbl = rf"$x_{{B}}\in[{X_SLICES[s][0]:.2f},{X_SLICES[s][1]:.2f}]$"
        h, = plt.plot([], [], color=XB_COLORS[s], lw=1.8, label=lbl)
        legend_handles.append(h); legend_labels.append(lbl)
    for h in legend_handles: h.remove()

    for c in range(NTCOLS):
        tmin, tmax = T_EDGES_POS[c], T_EDGES_POS[c+1]
        col_title = rf"$-t \in [{tmin:.2f}, {tmax:.2f})$"

        for r, (name, H1, ED, CE, LAB, LIM) in enumerate(row_vars):
            ax = axes[r, c]
            ax.set_xlim(*LIM)
            ax.set_xlabel(LAB)             # ← label every panel in the row
            if c == 0:
                ax.set_ylabel("counts")     # left column only

            if name == "phi":
                ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                              [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])

            ax.grid(True, linestyle="--", alpha=0.25)

            # draw the four xB-slice curves
            for sidx, s in enumerate(SLICE_ORDER):
                counts = H1[sidx, c, :]
                if np.sum(counts) == 0:
                    continue
                ax.plot(CE, counts, drawstyle="steps-mid", linewidth=1.6,
                        color=XB_COLORS[s], label=None)

            if r == 0:
                ax.text(0.02, 0.98, col_title, transform=ax.transAxes,
                        ha='left', va='top', fontsize=11)

            if c != 0:
                ax.tick_params(labelleft=False)
            # no longer hide bottom tick labels for non-bottom rows
            # (we want ticks + labels visible on all panels now)

    fig.legend(legend_handles, legend_labels, loc="upper center",
               ncol=4, frameon=True, fontsize=11, bbox_to_anchor=(0.5, 0.93))

    outpath = os.path.join(outdir, "binned_distributions.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
#endfor

def _y_formatter_hide_zero(hide_zero):
    def _fmt(v, pos):
        if hide_zero and np.isclose(v, 0.0, atol=1e-9): return ""
        return f"{v:g}"
    return FuncFormatter(_fmt)
#endfor

def _x_formatter_hide_07(hide_07):
    def _fmt(v, pos):
        if hide_07 and np.isclose(v, 0.7, atol=1e-9): return ""
        return f"{v:g}"
    return FuncFormatter(_fmt)
#endfor

def make_mx2_canvas(H1_MX2, outdir):
    """
    Mx2 canvas:
      • own y-scale; touching subplots; selective axis labels; red dashed fit curve
      • writes CSV of (mu, sigma) and RETURNS mu/sigma grids and success mask
    """
    fig, axes = plt.subplots(4, NTCOLS, figsize=(20, 12), sharex=True, sharey=False)
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.07, top=0.92, wspace=0.0, hspace=0.0)

    fig.suptitle(
        r"$M_{x}^{2}$ distributions by $x_{B}$ slice and $-t$ bin"
        "\n" + r"Global cuts: $Q^{2}>1,\, W>2,\, y<0.75$ (no $M_{x}^{2}$ cut)",
        fontsize=14
    )

    xfit = np.linspace(MX2_LIM[0], MX2_LIM[1], 400)

    # For return and for CSV
    mu_grid    = np.full((4, NTCOLS), np.nan, dtype=float)
    sigma_grid = np.full((4, NTCOLS), np.nan, dtype=float)
    ok_grid    = np.zeros((4, NTCOLS), dtype=bool)

    csv_path = os.path.join(outdir, "Mx2_fit_params.csv")
    rows_for_csv = []
    header = ["slice_name","x_min","x_max","t_min","t_max","mu","sigma","n_entries","fit_success"]

    for r, s in enumerate(SLICE_ORDER):
        xa, xb = X_SLICES[s]
        for c in range(NTCOLS):
            tmin, tmax = T_EDGES_POS[c], T_EDGES_POS[c+1]
            ax = axes[r, c]

            y = H1_MX2[r, c, :].astype(float)
            n_entries = int(np.sum(y))
            yerr = np.sqrt(np.maximum(y, 1.0))
            ax.errorbar(MX2_CEN, y, yerr=yerr, fmt='o', ms=3, lw=1.0, capsize=2, color='black')

            fit_success = False
            mu_val = np.nan
            sigma_val = np.nan
            if np.count_nonzero(y) >= 6:
                params, _yfit_centers, fit_success = _fit_gauss_plus_quad(MX2_CEN, y, yerr, mu0=M_P2, sigma0=0.03)
                A, mu, sigma, c0, c1, c2 = params
                mu_val, sigma_val = float(mu), float(sigma)
                yfit_curve = _gauss_plus_quad(xfit, A, mu, sigma, c0, c1, c2)
                ax.plot(xfit, yfit_curve, linestyle='--', linewidth=0.9, color='red',
                        label=rf"$\mu={mu:.4f},\ \sigma={sigma:.4f}$")
            #endif

            mu_grid[r, c]    = mu_val
            sigma_grid[r, c] = sigma_val
            ok_grid[r, c]    = bool(fit_success)

            ymax = float(np.max(y)) if y.size else 1.0
            ax.set_ylim(0.0, max(1.0, 1.35*ymax))
            ax.set_xlim(*MX2_LIM)

            if c == 0: ax.set_ylabel("counts")
            if r == 3: ax.set_xlabel(LAB_MX2)

            if c != 0: ax.tick_params(labelleft=False)
            if r != 3: ax.tick_params(labelbottom=False)

            if c == 0:
                ax.yaxis.set_major_formatter(_y_formatter_hide_zero(hide_zero=(r != 3)))
            #endif
            if r == 3:
                ax.xaxis.set_major_formatter(_x_formatter_hide_07(hide_07=(c != 0)))
            #endif

            ax.grid(True, linestyle="--", alpha=0.25)
            ax.text(0.02, 0.98, rf"$x_{{B}}\in[{xa:.2f},{xb:.2f}],\ -t\in[{tmin:.2f},{tmax:.2f})$",
                    transform=ax.transAxes, ha='left', va='top', fontsize=9)
            if fit_success:
                ax.legend(frameon=True, fontsize=8, loc="best")
            #endif

            rows_for_csv.append([s, f"{xa:.4f}", f"{xb:.4f}", f"{tmin:.4f}", f"{tmax:.4f}",
                                 f"{mu_val:.6f}" if np.isfinite(mu_val) else "",
                                 f"{sigma_val:.6f}" if np.isfinite(sigma_val) else "",
                                 str(n_entries),
                                 "1" if fit_success else "0"])
        #endfor
    #endfor

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows_for_csv)
    #endif

    outpath = os.path.join(outdir, "Mx2_distributions.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
    print(f"Wrote fit parameters: {csv_path}")

    return mu_grid, sigma_grid, ok_grid
#endfor

def make_mu_index_plot(mu_grid, sigma_grid, ok_grid, outdir):
    """
    Plot μ vs indexed bin (1..24), using σ as the vertical error bar.
    Bins 1–6: Low; 7–12: MidLow; 13–18: MidHigh; 19–24: High.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    total_bins = 6 * len(SLICE_ORDER)

    # Draw faint separators between x_B blocks
    for cut in [6, 12, 18]:
        ax.axvline(cut + 0.5, linestyle=":", linewidth=0.8, color="0.8")
    #endfor

    for r, s in enumerate(SLICE_ORDER):
        xa, xb = X_SLICES[s]
        idxs = np.arange(1 + 6*r, 1 + 6*(r+1))  # 1..6, 7..12, 13..18, 19..24
        ys   = mu_grid[r, :]
        yerr = sigma_grid[r, :]
        m    = ok_grid[r, :] & np.isfinite(ys) & np.isfinite(yerr)

        if np.any(m):
            ax.errorbar(idxs[m], ys[m], yerr=yerr[m],
                        fmt='o', ms=4, lw=1.0, capsize=3,
                        color=XB_COLORS[s],
                        label=rf"$x_{{B}}\in[{xa:.2f},{xb:.2f}]$")
        else:
            # Ensure legend entry exists even if no points
            ax.plot([], [], 'o', color=XB_COLORS[s],
                    label=rf"$x_{{B}}\in[{xa:.2f},{xb:.2f}]$")
        #endif
    #endfor

    ax.set_xlim(0.5, total_bins + 0.5)
    ax.set_xlabel("bin")
    ax.set_ylabel(r"$\mu$")
    ax.grid(True, linestyle="--", alpha=0.3, axis='y')
    ax.legend(frameon=True, ncol=4, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, 1.14))
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    outpath = os.path.join(outdir, "Mx2_mu_vs_bin.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
#endfor

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) != 2:
        print("Usage: python make_enpi_canvases.py <input.root>")
        sys.exit(1)
    #endif

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"[ERROR] File not found: {infile}")
        sys.exit(1)
    #endif

    outdir = ensure_outdir()

    # Single-pass accumulation over the tree (chunked)
    H2, H1_Q2, H1_X, H1_T, H1_PHI, H1_MX2 = accumulate(infile)

    # Render just the Mx2 figure and the μ–index plot per your latest request
    # mu_grid, sigma_grid, ok_grid = make_mx2_canvas(H1_MX2, outdir)
    # make_mu_index_plot(mu_grid, sigma_grid, ok_grid, outdir)

    # If you want the others again, uncomment:
    make_2d_canvases(H2, outdir)
    # make_binned_canvas_4x6(H1_Q2, H1_X, H1_T, H1_PHI, outdir)
#endfor

if __name__ == "__main__":
    main()
#endif