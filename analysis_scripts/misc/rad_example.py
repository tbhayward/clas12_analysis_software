#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DVCS radiative-correction panels (3x5) vs phi from four ROOT trees, with Poisson
error bars propagated for R_c = (gen_rad/rec_rad) / (gen_born/rec_born).

Inputs (TTree "PhysicsEvents"):
  Born:
    - gen: /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root
    - rec: /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_10594MeV.root
  Radiative:
    - gen: /volatile/clas12/thayward/temp_rad/gen_dvcsgen_sp18_inb_rad.root
    - rec: /volatile/clas12/thayward/temp_rad/rec_dvcsgen_sp18_inb_rad.root

Output:
  - output/dvcs_rad_example.pdf

Rows (top->bottom): Q^2, x_B, -t
Columns: 5 per row

Branches used: Q2, x, t1, phi2  (we define tpos = -t1)
Axis labels: x = r"$\phi$", y = r"$R_{c}$"

Poisson propagation (independent counts in a phi bin):
  Let A=gen_rad, B=rec_rad, C=gen_born, D=rec_born.
  R_c = (A*D) / (B*C)
  sigma(R_c) = R_c * sqrt(1/A + 1/D + 1/B + 1/C)
  Any bin with A==0 or B==0 or C==0 or D==0 -> NaN (not plotted).
"""

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

# -------------------- styling --------------------
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 220,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    # Use LaTeX-style mathtext without requiring a TeX install
    "mathtext.fontset": "dejavusans",
    "mathtext.default": "regular",
})

# -------------------- file paths --------------------
GEN_BORN_PATH = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root"
REC_BORN_PATH = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_10594MeV.root"
GEN_RAD_PATH  = "/volatile/clas12/thayward/temp_rad/gen_dvcsgen_sp18_inb_rad.root"
REC_RAD_PATH  = "/volatile/clas12/thayward/temp_rad/rec_dvcsgen_sp18_inb_rad.root"

TTREE_NAME = "PhysicsEvents"

# -------------------- bin definitions --------------------
Q2_BINS = [(1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.0), (5.0, 9.0)]
XB_BINS = [(0.10, 0.20), (0.20, 0.30), (0.30, 0.40), (0.40, 0.50), (0.50, 0.70)]
T_BINS  = [(0.10, 0.20), (0.20, 0.30), (0.30, 0.40), (0.40, 0.70), (0.70, 1.00)]

# -------------------- phi binning: force radians in [0, 2*pi] with 24 bins --------------------
PHI_NBINS = 24
PHI_MIN = 0.0
PHI_MAX = 2.0 * np.pi
PHI_EDGES = np.linspace(PHI_MIN, PHI_MAX, PHI_NBINS + 1)

# -------------------- helpers --------------------
def load_arrays(path):
    """
    Load required branches from ROOT file and return a dict of numpy arrays:
      'Q2', 'x', 'phi2', 'tpos'  (where tpos = -t1)
    """
    with uproot.open(path) as f:
        tree = f[TTREE_NAME]
        arr = tree.arrays(["Q2", "x", "t1", "phi2"], library="np")
    # Compute -t1 as positive
    tpos = -arr["t1"]
    return {
        "Q2": arr["Q2"],
        "x":  arr["x"],
        "phi2": arr["phi2"],
        "tpos": tpos,
    }
#endfor

def hist_counts(values, mask, bin_edges):
    """
    Histogram of 'values[mask]' into 'bin_edges'.
    """
    return np.histogram(values[mask], bins=bin_edges)[0]
#endfor

def rc_and_err(A, B, C, D):
    """
    Given arrays A,B,C,D>=0 (same shape), compute:
      R_c = (A*D) / (B*C)
      sigma(R_c) = R_c * sqrt(1/A + 1/D + 1/B + 1/C)
    If any of A,B,C,D==0 for a bin, return NaN for both R_c and sigma.
    """
    A = A.astype(float)
    B = B.astype(float)
    C = C.astype(float)
    D = D.astype(float)

    den = B * C
    num = A * D

    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(den > 0.0, num / den, np.nan)
        mask_pos = (A > 0.0) & (B > 0.0) & (C > 0.0) & (D > 0.0)
        sigma = np.full_like(R, np.nan, dtype=float)
        sigma[mask_pos] = R[mask_pos] * np.sqrt(1.0/A[mask_pos] + 1.0/D[mask_pos] + 1.0/B[mask_pos] + 1.0/C[mask_pos])
    #endif
    return R, sigma
#endfor

def compute_rc_curve_per_bin(gen_born, rec_born, gen_rad, rec_rad, varname, lo, hi, phi_edges):
    """
    For a given variable name ('Q2', 'x', 'tpos') and bin [lo, hi), compute
    the R_c(phi) curve and its Poisson error using histogrammed counts vs phi2.
    Returns phi_centers, R(phi), sigma_R(phi) (numpy arrays).
    """
    m_gb = (gen_born[varname] >= lo) & (gen_born[varname] < hi)
    m_rb = (rec_born[varname] >= lo) & (rec_born[varname] < hi)
    m_gr = (gen_rad[varname]  >= lo) & (gen_rad[varname]  < hi)
    m_rr = (rec_rad[varname]  >= lo) & (rec_rad[varname]  < hi)

    gb_phi = hist_counts(gen_born["phi2"], m_gb, phi_edges)
    rb_phi = hist_counts(rec_born["phi2"], m_rb, phi_edges)
    gr_phi = hist_counts(gen_rad["phi2"],  m_gr, phi_edges)
    rr_phi = hist_counts(rec_rad["phi2"],  m_rr, phi_edges)

    R, sigma = rc_and_err(gr_phi, rr_phi, gb_phi, rb_phi)
    centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    return centers, R, sigma
#endfor

# Nice LaTeX-like tick labels for phi in radians
_PHI_TICKS = [0.0, 0.5*np.pi, np.pi, 1.5*np.pi, 2.0*np.pi]
_PHI_TICKLABELS = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]

def format_axes(ax):
    """
    Aesthetics + LaTeX-style labels and fixed phi ticks.
    """
    ax.grid(True, alpha=0.35)
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"$\phi$", fontsize=12)
    ax.set_ylabel(r"$R_{c}$", fontsize=12)
    ax.set_xlim(PHI_MIN, PHI_MAX)
    ax.xaxis.set_major_locator(FixedLocator(_PHI_TICKS))
    ax.xaxis.set_major_formatter(FixedFormatter(_PHI_TICKLABELS))
    # Avoid y-axis scientific-offset like "1e-7 +"
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    #endfor
#endfor

# -------------------- plotting --------------------
def make_panels(gen_born, rec_born, gen_rad, rec_rad, out_pdf):
    """
    Build the 3x5 panel figure and save to 'out_pdf'.
    """
    fig, axes = plt.subplots(3, 5, figsize=(18, 10), constrained_layout=True)
    fig.suptitle(r"DVCS Radiative Correction $R_{c}$ vs $\phi$", fontsize=16)

    marker_fmt = "o-"

    # Row 0: Q^2 bins (5)
    for j, (lo, hi) in enumerate(Q2_BINS):
        ax = axes[0, j]
        phi_c, R, sR = compute_rc_curve_per_bin(gen_born, rec_born, gen_rad, rec_rad, "Q2", lo, hi, PHI_EDGES)
        ax.errorbar(phi_c, R, yerr=sR, fmt=marker_fmt, linewidth=1.6, markersize=3.0, capsize=2)
        ax.set_title(r"$Q^{2}\ \mathrm{in}\ [{:.2g}, {:.2g}]$".format(lo, hi))
        format_axes(ax)
        finite = np.isfinite(R)
        if np.any(finite):
            v = R[finite]
            lo_y, hi_y = np.nanpercentile(v, [5, 95])
            pad = 0.20 * max(1e-6, hi_y - lo_y)
            ax.set_ylim(lo_y - pad, hi_y + pad)
        #endif
    #endfor

    # Row 1: x_B bins (5)
    for j, (lo, hi) in enumerate(XB_BINS):
        ax = axes[1, j]
        phi_c, R, sR = compute_rc_curve_per_bin(gen_born, rec_born, gen_rad, rec_rad, "x", lo, hi, PHI_EDGES)
        ax.errorbar(phi_c, R, yerr=sR, fmt=marker_fmt, linewidth=1.6, markersize=3.0, capsize=2)
        ax.set_title(r"$x_{B}\ \mathrm{in}\ [{:.2g}, {:.2g}]$".format(lo, hi))
        format_axes(ax)
        finite = np.isfinite(R)
        if np.any(finite):
            v = R[finite]
            lo_y, hi_y = np.nanpercentile(v, [5, 95])
            pad = 0.20 * max(1e-6, hi_y - lo_y)
            ax.set_ylim(lo_y - pad, hi_y + pad)
        #endif
    #endfor

    # Row 2: -t bins (5) using tpos
    for j, (lo, hi) in enumerate(T_BINS):
        ax = axes[2, j]
        phi_c, R, sR = compute_rc_curve_per_bin(gen_born, rec_born, gen_rad, rec_rad, "tpos", lo, hi, PHI_EDGES)
        ax.errorbar(phi_c, R, yerr=sR, fmt=marker_fmt, linewidth=1.6, markersize=3.0, capsize=2)
        ax.set_title(r"$-t\ \mathrm{in}\ [{:.2g}, {:.2g}]$".format(lo, hi))
        format_axes(ax)
        finite = np.isfinite(R)
        if np.any(finite):
            v = R[finite]
            lo_y, hi_y = np.nanpercentile(v, [5, 95])
            pad = 0.20 * max(1e-6, hi_y - lo_y)
            ax.set_ylim(lo_y - pad, hi_y + pad)
        #endif
    #endfor

    # Save
    out_dir = os.path.dirname(out_pdf)
    if out_dir != "":
        os.makedirs(out_dir, exist_ok=True)
    #endif
    fig.savefig(out_pdf)
    print("Saved:", out_pdf)
#endfor

# -------------------- main --------------------
def main():
    # Load arrays once from each file
    gen_born = load_arrays(GEN_BORN_PATH)
    rec_born = load_arrays(REC_BORN_PATH)
    gen_rad  = load_arrays(GEN_RAD_PATH)
    rec_rad  = load_arrays(REC_RAD_PATH)

    # Build panels
    make_panels(gen_born, rec_born, gen_rad, rec_rad, out_pdf="output/dvcs_rad_example.pdf")
#endfor  # function ends; comment included per user preference

if __name__ == "__main__":
    main()
#endif