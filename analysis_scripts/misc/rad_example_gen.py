#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""
DVCS generator-level ratio plots of R_c^{(gen)}(φ):

(A) 3x5 panels of R_c^{(gen)}(φ) with Poisson error bars in fixed kinematic bins
(B) Experimental canvases: for each x_B bin, a grid with Q^2 bins as columns
    and -t bins as rows; each subplot shows R_c^{(gen)}(φ) in that (Q^2, -t) bin.

Generator-level definition (per φ bin):
  Let A = gen_rad (counts), C = gen_born (counts).
  R_c^{(gen)} = A / C
  σ(R_c^{(gen)}) = R_c^{(gen)} * sqrt(1/A + 1/C)
  If A==0 or C==0 in a φ bin, that point is undefined (NaN -> not plotted).

Inputs (TTree "PhysicsEvents"):
  Born (generator):
    - /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root
  Radiative (generator):
    - /volatile/clas12/thayward/temp_rad/gen_dvcsgen_sp18_inb_rad.root

Experimental bin CSV (two header lines, whitespace-separated columns):
  /u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_bin_v2.csv
  Fallbacks: ./integrated_bin_v2.csv or /mnt/data/integrated_bin_v2.csv

Outputs:
  - output/rad_example/genborn_phi_panels.pdf
  - output/rad_example/genborn_phi_grid_xB_<lo>_<hi>.pdf   (one per x_B bin)

Style:
  - LaTeX-like mathtext labels (no external TeX needed)
  - Closed axes boxes; marker-only points; y range fixed to [0,2]
  - φ ticks at 0, π/2, π, 3π/2, 2π
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
    # LaTeX-like mathtext without external TeX
    "mathtext.fontset": "dejavusans",
    "mathtext.default": "regular",
})

# -------------------- file paths --------------------
GEN_BORN_PATH = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root"
GEN_RAD_PATH  = "/volatile/clas12/thayward/temp_rad/gen_dvcsgen_sp18_inb_rad.root"

TTREE_NAME = "PhysicsEvents"

CSV_CANDIDATES = [
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_bin_v2.csv",
    "./integrated_bin_v2.csv",
    "/mnt/data/integrated_bin_v2.csv",
]

# -------------------- fixed panel bin definitions for part (A) --------------------
Q2_BINS_FIXED = [(1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.0), (5.0, 9.0)]
XB_BINS_FIXED = [(0.10, 0.20), (0.20, 0.30), (0.30, 0.40), (0.40, 0.50), (0.50, 0.70)]
T_BINS_FIXED  = [(0.10, 0.20), (0.20, 0.30), (0.30, 0.40), (0.40, 0.70), (0.70, 1.00)]

# -------------------- φ binning: radians in [0, 2π] with 24 bins --------------------
PHI_NBINS = 24
PHI_MIN = 0.0
PHI_MAX = 2.0 * np.pi
PHI_EDGES = np.linspace(PHI_MIN, PHI_MAX, PHI_NBINS + 1)

# -------------------- helpers: IO and math --------------------
def load_arrays(path):
    """
    Load required branches from ROOT file and return a dict of numpy arrays:
      'Q2', 'x', 'phi', 'tpos'  (phi wrapped into [0, 2π), tpos = -t1 > 0)
    """
    with uproot.open(path) as f:
        tree = f[TTREE_NAME]
        arr = tree.arrays(["Q2", "x", "t1", "phi2"], library="np")
    # Wrap φ into [0, 2π)
    phi_wrapped = np.mod(arr["phi2"], 2.0 * np.pi)
    tpos = -arr["t1"]  # positive -t
    return {"Q2": arr["Q2"], "x": arr["x"], "phi": phi_wrapped, "tpos": tpos}
#endfor

def hist_counts(values, mask, bin_edges):
    """Histogram of values[mask] into bin_edges."""
    return np.histogram(values[mask], bins=bin_edges)[0]
#endfor

def rgen_and_err(A, C):
    """
    Generator-level ratio and uncertainty:
      R_c^{(gen)} = A / C
      σ = R_c^{(gen)} * sqrt(1/A + 1/C)
    Any zero in A or C -> NaN for both.
    """
    A = A.astype(float)
    C = C.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(C > 0.0, A / C, np.nan)
        mask_pos = (A > 0.0) & (C > 0.0)
        sR = np.full_like(R, np.nan, dtype=float)
        sR[mask_pos] = R[mask_pos] * np.sqrt(1.0 / A[mask_pos] + 1.0 / C[mask_pos])
    return R, sR
#endfor

# Nice LaTeX-style φ ticks
_PHI_TICKS = [0.0, 0.5*np.pi, np.pi, 1.5*np.pi, 2.0*np.pi]
_PHI_TICKLABELS = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]

def format_axes_phi(ax):
    """Marker-only style, closed box, fixed y in [0,2], φ ticks, dashed y=1 line."""
    ax.grid(True, alpha=0.35)
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$R_{c}^{(\mathrm{gen})}$")
    ax.set_xlim(PHI_MIN, PHI_MAX)
    ax.set_ylim(0.0, 2.0)
    ax.xaxis.set_major_locator(FixedLocator(_PHI_TICKS))
    ax.xaxis.set_major_formatter(FixedFormatter(_PHI_TICKLABELS))
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(True)
    #endfor
#endfor

# -------------------- part (A): fixed 3×5 panels --------------------
def compute_rgen_curve_1D(gen_born, gen_rad, varname, lo, hi, phi_edges):
    """R_c^{(gen)}(φ) inside [lo,hi) of single variable varname."""
    m_gb = (gen_born[varname] >= lo) & (gen_born[varname] < hi)
    m_gr = (gen_rad[varname]  >= lo) & (gen_rad[varname]  < hi)

    gb_phi = hist_counts(gen_born["phi"], m_gb, phi_edges)  # C
    gr_phi = hist_counts(gen_rad["phi"],  m_gr, phi_edges)  # A

    R, sR = rgen_and_err(gr_phi, gb_phi)
    centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    return centers, R, sR
#endfor

def make_phi_panels(gen_born, gen_rad, out_pdf):
    """Create the 3x5 generator-level R_c^{(gen)}(φ) figure and save."""
    fig, axes = plt.subplots(3, 5, figsize=(18, 10), constrained_layout=True)
    fig.suptitle(r"DVCS Generator-Level Ratio $R_{c}^{(\mathrm{gen})}$ vs $\phi$")

    # Row 0: Q^2 bins
    for j, (lo, hi) in enumerate(Q2_BINS_FIXED):
        ax = axes[0, j]
        phi_c, R, sR = compute_rgen_curve_1D(gen_born, gen_rad, "Q2", lo, hi, PHI_EDGES)
        ax.errorbar(phi_c, R, yerr=sR, fmt="o", linestyle="none", linewidth=1.2, markersize=3.0, capsize=2)
        ax.set_title(rf"$Q^{{2}}\ \mathrm{{in}}\ [{lo:.2g}, {hi:.2g}]$")
        format_axes_phi(ax)
    #endfor

    # Row 1: x_B bins
    for j, (lo, hi) in enumerate(XB_BINS_FIXED):
        ax = axes[1, j]
        phi_c, R, sR = compute_rgen_curve_1D(gen_born, gen_rad, "x", lo, hi, PHI_EDGES)
        ax.errorbar(phi_c, R, yerr=sR, fmt="o", linestyle="none", linewidth=1.2, markersize=3.0, capsize=2)
        ax.set_title(rf"$x_{{B}}\ \mathrm{{in}}\ [{lo:.2g}, {hi:.2g}]$")
        format_axes_phi(ax)
    #endfor

    # Row 2: -t bins (use tpos)
    for j, (lo, hi) in enumerate(T_BINS_FIXED):
        ax = axes[2, j]
        phi_c, R, sR = compute_rgen_curve_1D(gen_born, gen_rad, "tpos", lo, hi, PHI_EDGES)
        ax.errorbar(phi_c, R, yerr=sR, fmt="o", linestyle="none", linewidth=1.2, markersize=3.0, capsize=2)
        ax.set_title(rf"$-t\ \mathrm{{in}}\ [{lo:.2g}, {hi:.2g}]$")
        format_axes_phi(ax)
    #endfor

    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf)
    print("Saved:", out_pdf)
#endfor

# -------------------- part (B): experimental canvases per x_B --------------------
class BinRow(object):
    __slots__ = ("xbmin","xbmax","q2min","q2max","tmin","tmax")
    def __init__(self, xbmin, xbmax, q2min, q2max, tmin, tmax):
        self.xbmin=xbmin; self.xbmax=xbmax
        self.q2min=q2min; self.q2max=q2max
        self.tmin=tmin;   self.tmax=tmax
#endfor

def parse_binning_csv():
    """
    Reads the experimental bin CSV (two header lines, whitespace split).
    Expected columns (0-based token index):
      4: xBmin, 5: xBmax, 7: Q2min, 8: Q2max, 10: |tmin|, 11: |tmax|
    Returns: list[BinRow]
    """
    csv_path = next((p for p in CSV_CANDIDATES if os.path.isfile(p)), None)
    if not csv_path:
        raise FileNotFoundError("Could not find integrated_bin_v2.csv in any candidate path.")
    rows = []
    with open(csv_path, "r") as f:
        lines = [ln.rstrip("\r\n") for ln in f.readlines()]
    if len(lines) <= 2:
        return rows
    for i in range(2, len(lines)):
        ln = lines[i]
        if not ln or ln.lstrip().startswith("#"):
            continue
        toks = ln.split()
        if len(toks) < 12:
            continue
        try:
            xbmin = float(toks[4]);  xbmax = float(toks[5])
            q2min = float(toks[7]);  q2max = float(toks[8])
            tmin  = abs(float(toks[10]));  tmax  = abs(float(toks[11]))
            rows.append(BinRow(xbmin, xbmax, q2min, q2max, tmin, tmax))
        except Exception:
            continue
    return rows
#endfor

def compute_rgen_curve_cell(gen_born, gen_rad, xb, q2, tt, phi_edges):
    """
    R_c^{(gen)}(φ) in the 3D window:
      x_B in xb = (xbmin, xbmax),
      Q^2 in q2 = (q2min, q2max),
      -t in tt  = (tmin,  tmax)
    """
    m_gb = ((gen_born["x"] >= xb[0]) & (gen_born["x"] < xb[1]) &
            (gen_born["Q2"] >= q2[0]) & (gen_born["Q2"] < q2[1]) &
            (gen_born["tpos"] >= tt[0]) & (gen_born["tpos"] < tt[1]))
    m_gr = ((gen_rad["x"] >= xb[0]) & (gen_rad["x"] < xb[1]) &
            (gen_rad["Q2"] >= q2[0]) & (gen_rad["Q2"] < q2[1]) &
            (gen_rad["tpos"] >= tt[0]) & (gen_rad["tpos"] < tt[1]))

    gb_phi = hist_counts(gen_born["phi"], m_gb, phi_edges)  # C
    gr_phi = hist_counts(gen_rad["phi"],  m_gr, phi_edges)  # A

    R, sR = rgen_and_err(gr_phi, gb_phi)
    centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    return centers, R, sR
#endfor

def make_phi_grid_for_xb(gen_born, gen_rad, xbmin, xbmax, binrows, out_dir):
    """
    For a fixed x_B range [xbmin,xbmax), build a grid of subplots:
      columns = sorted unique Q^2 bins in this x_B group
      rows    = sorted unique -t  bins in this x_B group
    Each subplot shows R_c^{(gen)}(φ) in that (Q^2,-t) bin.
    """
    q2_bins = sorted({(r.q2min, r.q2max) for r in binrows})
    t_bins  = sorted({(r.tmin,  r.tmax)  for r in binrows})
    nC = len(q2_bins); nR = len(t_bins)
    if nC == 0 or nR == 0:
        return

    # Figure size scales with grid density
    fig_w = max(10.0, 3.4 * nC + 2.0)
    fig_h = max(7.0, 2.8 * nR + 2.0)
    fig, axes = plt.subplots(nR, nC, figsize=(fig_w, fig_h), constrained_layout=True, squeeze=False)
    fig.suptitle(rf"Generator-Level $R_{{c}}^{{(\mathrm{{gen}})}}(\phi)$ for $x_{{B}}$ in [{xbmin:.2g}, {xbmax:.2g}]")

    for i_t, tt in enumerate(t_bins):
        for j_q, q2 in enumerate(q2_bins):
            ax = axes[i_t, j_q]
            phi_c, R, sR = compute_rgen_curve_cell(gen_born, gen_rad, (xbmin, xbmax), q2, tt, PHI_EDGES)
            ax.errorbar(phi_c, R, yerr=sR, fmt="o", linestyle="none",
                        linewidth=1.2, markersize=3.0, capsize=2)
            ax.set_title(rf"$Q^{{2}}\!\in\![{q2[0]:.2g},{q2[1]:.2g}],\ -t\!\in\![{tt[0]:.2g},{tt[1]:.2g}]$")
            format_axes_phi(ax)
        #endfor
    #endfor

    os.makedirs(out_dir, exist_ok=True)
    out_pdf = os.path.join(out_dir, f"genborn_phi_grid_xB_{xbmin:.2g}_{xbmax:.2g}.pdf")
    fig.savefig(out_pdf)
    print("Saved:", out_pdf)
#endfor

def make_all_xb_phi_grids(gen_born, gen_rad, out_dir):
    """Parse CSV, group rows by x_B bin, and write one generator-level φ-grid per x_B."""
    rows = parse_binning_csv()
    if not rows:
        print("WARNING: No rows parsed from CSV; skipping experimental canvases.")
        return
    # Group by exact x_B tuple
    groups = {}
    for r in rows:
        key = (r.xbmin, r.xbmax)
        groups.setdefault(key, []).append(r)
    #endfor
    for (xbmin, xbmax), binrows in sorted(groups.items()):
        make_phi_grid_for_xb(gen_born, gen_rad, xbmin, xbmax, binrows, out_dir)
    #endfor
#endfor

# -------------------- main --------------------
def main():
    # Load arrays once from each generator file
    gen_born = load_arrays(GEN_BORN_PATH)
    gen_rad  = load_arrays(GEN_RAD_PATH)

    # Single output root
    out_root = "output/rad_example"
    os.makedirs(out_root, exist_ok=True)

    # (A) generator-level 3×5 φ-panels
    make_phi_panels(gen_born, gen_rad, out_pdf=os.path.join(out_root, "genborn_phi_panels.pdf"))

    # (B) experimental canvases per x_B (each subplot is R_c^{(gen)}(φ) for a (Q^2, -t) cell)
    make_all_xb_phi_grids(gen_born, gen_rad, out_dir=out_root)
#endfor  # function ends; comment included per user preference

if __name__ == "__main__":
    main()
#endif