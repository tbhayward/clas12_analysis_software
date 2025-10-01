#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DVCS radiative-correction panels (3x5) vs #phi from four ROOT trees, with Poisson
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

Rows (top->bottom): Q2, xB, -t
Columns: 5 per row

Branches used: Q2, x, t1, phi2  (we define tpos = -t1)
Axis labels: x = "#phi", y = "R_{c}"

Poisson propagation (independent counts):
  Let A=gen_rad, B=rec_rad, C=gen_born, D=rec_born in a given phi bin.
  R_c = (A*D) / (B*C)
  sigma(R_c) = R_c * sqrt(1/A + 1/D + 1/B + 1/C)
  Points with any of A,B,C,D == 0 -> NaN (not plotted).
"""

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

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

# Number of bins in #phi; auto-detect degrees vs radians
PHI_NBINS_DEG = 24
PHI_NBINS_RAD = 24

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

def decide_phi_bins(phi_arrays):
    """
    Decide whether phi2 is in degrees or radians by inspecting its scale.
    Returns bin_edges (numpy array) and a label suffix string.
    """
    sample = []
    for a in phi_arrays:
        if a.size > 0:
            n = min(5000, a.size)
            sample.append(a[np.random.choice(a.size, size=n, replace=False)])
        #endif
    #endfor
    if len(sample) == 0:
        return np.linspace(-180.0, 180.0, PHI_NBINS_DEG + 1), "deg"
    #endif
    probe = np.concatenate(sample)
    max_abs = np.nanmax(np.abs(probe))
    if np.isfinite(max_abs) and max_abs < 3.5:
        return np.linspace(-np.pi, np.pi, PHI_NBINS_RAD + 1), "rad"
    else:
        return np.linspace(-180.0, 180.0, PHI_NBINS_DEG + 1), "deg"
    #endif
#endfor

def hist_counts(values, mask, bin_edges):
    """
    Fast histogram of 'values[mask]' into 'bin_edges'.
    """
    return np.histogram(values[mask], bins=bin_edges)[0]
#endfor

def rc_and_err(A, B, C, D):
    """
    Given four non-negative integer arrays A,B,C,D (same shape), compute:
      R_c = (A*D) / (B*C)
      sigma(R_c) = R_c * sqrt(1/A + 1/D + 1/B + 1/C)
    For any element where any of A,B,C,D == 0, return NaN for both R_c and sigma.
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

def format_axes(ax):
    """
    Light aesthetics for clarity.
    """
    ax.grid(True, alpha=0.35)
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel("#phi", fontsize=11)
    ax.set_ylabel("R_{c}", fontsize=11)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    #endfor
#endfor

# -------------------- plotting --------------------
def make_panels(gen_born, rec_born, gen_rad, rec_rad, out_pdf):
    """
    Build the 3x5 panel figure and save to 'out_pdf'.
    """
    # Decide phi binning from actual data
    phi_edges, _phi_units = decide_phi_bins([gen_born["phi2"], rec_born["phi2"], gen_rad["phi2"], rec_rad["phi2"]])

    # Figure setup
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 200,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
    })
    fig, axes = plt.subplots(3, 5, figsize=(18, 10), constrained_layout=True)
    fig.suptitle("DVCS Radiative Correction R_{c} vs #phi", fontsize=16)

    marker_fmt = "o-"

    # Row 0: Q2 bins (5)
    for j, (lo, hi) in enumerate(Q2_BINS):
        ax = axes[0, j]
        phi_c, R, sR = compute_rc_curve_per_bin(gen_born, rec_born, gen_rad, rec_rad, "Q2", lo, hi, phi_edges)
        ax.errorbar(phi_c, R, yerr=sR, fmt=marker_fmt, linewidth=1.6, markersize=3.0, capsize=2)
        ax.set_title("Q2 in [{:.2g}, {:.2g}]".format(lo, hi))
        format_axes(ax)
        finite = np.isfinite(R)
        if np.any(finite):
            v = R[finite]
            lo_y, hi_y = np.nanpercentile(v, [5, 95])
            pad = 0.20 * max(1e-6, hi_y - lo_y)
            ax.set_ylim(lo_y - pad, hi_y + pad)
        #endif
    #endfor

    # Row 1: xB bins (5)
    for j, (lo, hi) in enumerate(XB_BINS):
        ax = axes[1, j]
        phi_c, R, sR = compute_rc_curve_per_bin(gen_born, rec_born, gen_rad, rec_rad, "x", lo, hi, phi_edges)
        ax.errorbar(phi_c, R, yerr=sR, fmt=marker_fmt, linewidth=1.6, markersize=3.0, capsize=2)
        ax.set_title("xB in [{:.2g}, {:.2g}]".format(lo, hi))
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
        phi_c, R, sR = compute_rc_curve_per_bin(gen_born, rec_born, gen_rad, rec_rad, "tpos", lo, hi, phi_edges)
        ax.errorbar(phi_c, R, yerr=sR, fmt=marker_fmt, linewidth=1.6, markersize=3.0, capsize=2)
        ax.set_title("-t in [{:.2g}, {:.2g}]".format(lo, hi))
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