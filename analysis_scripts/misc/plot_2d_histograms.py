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
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

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

# 1D base bin counts (2D gets 2/3)
NB_Q2  = 60
NB_X   = 70
NB_T   = 65
NB_PHI = 72
NB_MX2 = 80

# 2D binning = 0.6667 × base
NB_Q2_2D  = int(round(NB_Q2  * 2/3))   # 40
NB_X_2D   = int(round(NB_X   * 2/3))   # 47
NB_T_2D   = int(round(NB_T   * 2/3))   # 43
NB_PHI_2D = int(round(NB_PHI * 2/3))   # 48

# Labels (units in parentheses)
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
    edges = np.array([0.10, 0.25, 0.35, 0.45, 0.60])
    # digitize gives 1..len(edges)-1 for interior; we want 0..3
    # Handle open intervals (xa < x < xb): exclude exact edges to match prior (> and <)
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
    # valid 0..5
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
        return p, _gauss_plus_quad(xc, *p)
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
            return tuple(popt), yfit
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
    return params, yhat_b
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

                # Accumulate with np.add.at
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
    fig.subplots_adjust(top=0.88, wspace=0.15, hspace=0.25)
    fig.suptitle(
        r"Overlaid distributions by $-t$ bin and $x_{B}$ slice" "\n"
        r"Global cuts: $Q^{2}>1,\, W>2,\, y<0.75,\, 0.81<M_{x}^{2}<1.00$",
        fontsize=14
    )

    # Prebuild legend handles
    legend_handles = []
    legend_labels  = []
    for s in SLICE_ORDER:
        lbl = rf"$x_{{B}}\in[{X_SLICES[s][0]:.2f},{X_SLICES[s][1]:.2f}]$"
        h, = plt.plot([], [], color=XB_COLORS[s], lw=1.8, label=lbl)
        legend_handles.append(h); legend_labels.append(lbl)
    #endfor
    for h in legend_handles: h.remove()

    for c in range(NTCOLS):
        tmin, tmax = T_EDGES_POS[c], T_EDGES_POS[c+1]
        col_title = rf"$-t \in [{tmin:.2f}, {tmax:.2f})$ (GeV$^{{2}}$)"

        for r, (name, H1, ED, CE, LAB, LIM) in enumerate(row_vars):
            ax = axes[r, c]
            ax.set_xlim(*LIM); ax.set_xlabel(LAB)
            if c == 0: ax.set_ylabel("counts")
            if name == "phi":
                ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                              [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
            #endif
            ax.grid(True, linestyle="--", alpha=0.25)

            for sidx, s in enumerate(SLICE_ORDER):
                counts = H1[sidx, c, :]
                if np.sum(counts) == 0: 
                    continue
                #endif
                ax.plot(CE, counts, drawstyle="steps-mid", linewidth=1.6, color=XB_COLORS[s],
                        label=None)
            #endfor

            if r == 0:
                ax.set_title(col_title, fontsize=12)
            #endif
        #endfor
    #endfor

    fig.legend(legend_handles, legend_labels, loc="upper center",
               ncol=4, frameon=True, fontsize=11, bbox_to_anchor=(0.5, 0.93))

    outpath = os.path.join(outdir, "binned_distributions.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
#endfor

def make_mx2_canvas(H1_MX2, outdir):
    fig, axes = plt.subplots(4, NTCOLS, figsize=(20, 12), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.90, wspace=0.15, hspace=0.25)
    fig.suptitle(
        r"$M_{x}^{2}$ distributions by $x_{B}$ slice and $-t$ bin"
        "\n" + r"Global cuts: $Q^{2}>1,\, W>2,\, y<0.75$ (no $M_{x}^{2}$ cut)",
        fontsize=14
    )

    for r, s in enumerate(SLICE_ORDER):
        xa, xb = X_SLICES[s]
        for c in range(NTCOLS):
            tmin, tmax = T_EDGES_POS[c], T_EDGES_POS[c+1]
            ax = axes[r, c]

            y = H1_MX2[r, c, :].astype(float)
            yerr = np.sqrt(np.maximum(y, 1.0))  # simple Poisson
            ax.errorbar(MX2_CEN, y, yerr=yerr, fmt='o', ms=3, lw=1.0, capsize=2,
                        color='black', label=None)

            # Fit (Gaussian + quadratic), dashed curve; legend shows mu & sigma only
            label_txt = r"$\mu=\mathrm{N/A},\ \sigma=\mathrm{N/A}$"
            if np.count_nonzero(y) >= 6:
                params, yfit = _fit_gauss_plus_quad(MX2_CEN, y, yerr, mu0=M_P2, sigma0=0.03)
                A, mu, sigma, c0, c1, c2 = params
                ax.plot(MX2_CEN, yfit, linestyle='--', linewidth=1.0, color='black',
                        label=rf"$\mu={mu:.4f},\ \sigma={sigma:.4f}$")
            else:
                ax.plot([], [], linestyle='--', linewidth=1.0, color='black', label=label_txt)
            #endif

            ax.set_xlim(*MX2_LIM); ax.set_ylim(bottom=0)
            ax.set_xlabel(LAB_MX2); ax.set_ylabel("counts")
            ax.grid(True, linestyle="--", alpha=0.25)
            ax.set_title(
                rf"$x_{{B}}\in[{xa:.2f},{xb:.2f}],\ -t\in[{tmin:.2f},{tmax:.2f})$ (GeV$^{{2}}$)",
                fontsize=10
            )
            ax.legend(frameon=True, fontsize=8, loc="best")
        #endfor
    #endfor

    outpath = os.path.join(outdir, "Mx2_distributions.pdf")
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

    # Render figures from the accumulated histograms
    # make_2d_canvases(H2, outdir)
    # make_binned_canvas_4x6(H1_Q2, H1_X, H1_T, H1_PHI, outdir)
    make_mx2_canvas(H1_MX2, outdir)
#endfor

if __name__ == "__main__":
    main()
#endif