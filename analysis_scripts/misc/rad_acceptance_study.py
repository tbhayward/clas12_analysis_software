#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""
Compare Born vs Radiative acceptances A(φ) in experimental bins, with legends
showing **reconstructed-level** mean kinematics per subplot:

For each x_B bin (from CSV), make a canvas with columns = Q^2 bins and rows = -t bins.
Each subplot shows:
  - Born acceptance A_B(φ) = N_rec_B(φ) / N_gen_B(φ),
  - Rad  acceptance A_R(φ) = N_rec_R(φ) / N_gen_R(φ),
with Poisson ratio errors per φ bin, markers only, and a legend reporting:
  Born (reco):  ⟨Q²⟩, ⟨x_B⟩, ⟨−t⟩, ⟨y⟩
  Rad  (reco):  ⟨Q²⟩, ⟨x_B⟩, ⟨−t⟩, ⟨y⟩
computed from **reconstructed events inside that bin**.

Inputs (TTree "PhysicsEvents"):
  Born:
    - gen: /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_10604MeV.root
    - rec: /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_10604MeV.root
  Radiative:
    - gen: /volatile/clas12/thayward/temp_rad/gen_dvcsgen_fa18_out_rad.root
    - rec: /volatile/clas12/thayward/temp_rad/rec_dvcsgen_fa18_out_rad.root

Experimental bin CSV (two header lines, whitespace-separated columns):
  /u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_bin_v2.csv
  Fallbacks: ./integrated_bin_v2.csv or /mnt/data/integrated_bin_v2.csv

Output:
  - output/rad_example/acc_phi_grid_xB_<xb_lo>_<xb_hi>.pdf  (one canvas per x_B bin)

Style:
  - LaTeX-like mathtext (no external TeX needed)
  - Closed axes boxes; markers only; φ bins = 24 in [0, 2π]
  - y-range fixed to [0, 1.05]; φ ticks at 0, π/2, π, 3π/2, 2π
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
    "legend.fontsize": 8,
    "mathtext.fontset": "dejavusans",
    "mathtext.default": "regular",
})

# -------------------- file paths --------------------
GEN_BORN_PATH = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root"
REC_BORN_PATH = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_10594MeV.root"
GEN_RAD_PATH  = "/volatile/clas12/thayward/temp_rad/gen_sp18_inb_rad.root"
REC_RAD_PATH  = "/volatile/clas12/thayward/temp_rad/rec_sp18_inb_rad.root"

# GEN_BORN_PATH = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_10604MeV.root"
# REC_BORN_PATH = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_10604MeV.root"
# GEN_RAD_PATH  = "/volatile/clas12/thayward/temp_rad/gen_fa18_out_rad.root"
# REC_RAD_PATH  = "/volatile/clas12/thayward/temp_rad/rec_fa18_out_rad.root"

TTREE_NAME = "PhysicsEvents"

CSV_CANDIDATES = [
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_bin_v2.csv",
    "./integrated_bin_v2.csv",
    "/mnt/data/integrated_bin_v2.csv",
]

# -------------------- φ binning: radians in [0, 2π] with 24 bins --------------------
PHI_NBINS = 24
PHI_MIN = 0.0
PHI_MAX = 2.0 * np.pi
PHI_EDGES = np.linspace(PHI_MIN, PHI_MAX, PHI_NBINS + 1)

# -------------------- helpers: IO and math --------------------
def load_arrays(path):
    """
    Load branches from ROOT file and return dict of numpy arrays:
      'Q2', 'x', 'y', 'phi', 'tpos'
    φ is wrapped into [0, 2π); tpos = -t1 (>0).
    Robust to missing 'y' branch (fills with NaNs).
    """
    with uproot.open(path) as f:
        tree = f[TTREE_NAME]
        arr = tree.arrays(["Q2", "x", "t1", "phi2"], library="np")
        # Try to get 'y'; if missing, fill with NaNs of the right length
        try:
            yarr = tree["y"].array(library="np")
        except Exception:
            # Make a NaN array aligned with Q2 length
            yarr = np.full_like(arr["Q2"], np.nan, dtype=float)
        #endif
    phi_wrapped = np.mod(arr["phi2"], 2.0 * np.pi)
    tpos = -arr["t1"]
    return {"Q2": arr["Q2"], "x": arr["x"], "y": yarr, "phi": phi_wrapped, "tpos": tpos}
#endfor

def hist_counts(values, mask, bin_edges):
    """Histogram of values[mask] into bin_edges."""
    return np.histogram(values[mask], bins=bin_edges)[0]
#endfor

def ratio_with_poisson_errors(num, den):
    """
    Given two non-negative arrays num, den (same shape) of counts per bin,
    compute R = num/den with Poisson error propagation:
        σ = R * sqrt(1/num + 1/den)
    Any bin with num==0 or den==0 -> NaN for both R and σ.
    """
    num = num.astype(float)
    den = den.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(den > 0.0, num / den, np.nan)
        mask = (num > 0.0) & (den > 0.0)
        sR = np.full_like(R, np.nan, dtype=float)
        sR[mask] = R[mask] * np.sqrt(1.0 / num[mask] + 1.0 / den[mask])
    return R, sR
#endfor

def binomial_error(nrec, ngen):
    """Binomial error for p̂ = nrec/ngen; returns (p̂, σ) or (nan, nan) if ngen==0."""
    if ngen <= 0:
        return np.nan, np.nan
    p = nrec / float(ngen)
    s = np.sqrt(p * max(0.0, 1.0 - p) / float(ngen))
    return p, s
#endfor

def _safe_mean(a):
    """mean(a) with empty-guard -> np.nan."""
    return float(np.mean(a)) if a.size > 0 else np.nan
#endfor

def fmt_float(v, fmt="{:.3f}"):
    """Format float or return '--' if NaN/inf."""
    return ("--" if (v is None or not np.isfinite(v)) else fmt.format(v))
#endfor

# Pretty φ ticks
_PHI_TICKS = [0.0, 0.5*np.pi, np.pi, 1.5*np.pi, 2.0*np.pi]
_PHI_TICKLABELS = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]

def format_axes_phi(ax):
    """Closed box, grid, φ ticks, markers-only y in [0,1.05], dashed y=1 line."""
    ax.grid(True, alpha=0.35)
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$A(\phi)$")
    ax.set_xlim(PHI_MIN, PHI_MAX)
    ax.set_ylim(0.0, 1.05)
    ax.xaxis.set_major_locator(FixedLocator(_PHI_TICKS))
    ax.xaxis.set_major_formatter(FixedFormatter(_PHI_TICKLABELS))
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(True)
    #endfor
#endfor

# -------------------- CSV parsing --------------------
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

# -------------------- acceptance + reco-mean-kinematics calculations --------------------
def acceptance_phi_for_cell(gen_dict, rec_dict, xb, q2, tt, phi_edges):
    """
    Compute acceptance vs φ inside a 3D kinematic cell:
      x_B in xb = (xbmin, xbmax), Q^2 in q2 = (q2min, q2max), -t in tt = (tmin, tmax).

    Returns:
      centers, A(φ), σ_A(φ),
      Nrec_tot, Ngen_tot, A_int, σ_int,
      mean_Q2_reco, mean_x_reco, mean_tpos_reco(=−t), mean_y_reco
      (the means are computed from **reconstructed** events in the cell)
    """
    # Masks for gen & rec in the same (xB,Q2,t) window
    m_gen = ((gen_dict["x"]  >= xb[0]) & (gen_dict["x"]  < xb[1]) &
             (gen_dict["Q2"] >= q2[0]) & (gen_dict["Q2"] < q2[1]) &
             (gen_dict["tpos"]>= tt[0]) & (gen_dict["tpos"]< tt[1]))
    m_rec = ((rec_dict["x"]  >= xb[0]) & (rec_dict["x"]  < xb[1]) &
             (rec_dict["Q2"] >= q2[0]) & (rec_dict["Q2"] < q2[1]) &
             (rec_dict["tpos"]>= tt[0]) & (rec_dict["tpos"]< tt[1]))

    # Per-φ histograms
    gen_phi = hist_counts(gen_dict["phi"], m_gen, phi_edges)
    rec_phi = hist_counts(rec_dict["phi"], m_rec, phi_edges)

    # Acceptance per φ with Poisson ratio error
    Aphi, sAphi = ratio_with_poisson_errors(rec_phi, gen_phi)
    centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])

    # Integrated acceptance & binomial error from totals in this cell
    Ngen_tot = int(np.count_nonzero(m_gen))
    Nrec_tot = int(np.count_nonzero(m_rec))
    A_int, sA_int = binomial_error(Nrec_tot, Ngen_tot)

    # **Reconstructed-level** mean kinematics for legend
    mean_Q2   = _safe_mean(rec_dict["Q2"][m_rec])
    mean_x    = _safe_mean(rec_dict["x"][m_rec])
    mean_tpos = _safe_mean(rec_dict["tpos"][m_rec])   # this is ⟨−t⟩
    mean_y    = _safe_mean(rec_dict["y"][m_rec])

    return centers, Aphi, sAphi, Nrec_tot, Ngen_tot, A_int, sA_int, mean_Q2, mean_x, mean_tpos, mean_y
#endfor

# -------------------- plotting per x_B canvases --------------------
def make_acceptance_grids(gen_born, rec_born, gen_rad, rec_rad, out_dir):
    """
    Parse CSV, group rows by x_B, then for each x_B bin produce a grid:
      columns = Q^2 bins, rows = -t bins. Subplots overlay Born vs Rad A(φ),
      and the legend shows **reconstructed** mean kinematics (⟨Q²⟩, ⟨x_B⟩, ⟨−t⟩, ⟨y⟩)
      for Born and Rad within that subplot's bin.
    """
    rows = parse_binning_csv()
    if not rows:
        print("WARNING: No rows parsed from CSV; skipping canvases.")
        return
    #endif

    # Group rows by exact (xbmin, xbmax)
    groups = {}
    for r in rows:
        key = (r.xbmin, r.xbmax)
        groups.setdefault(key, []).append(r)
    #endfor

    os.makedirs(out_dir, exist_ok=True)

    for (xbmin, xbmax), binrows in sorted(groups.items()):
        q2_bins = sorted({(r.q2min, r.q2max) for r in binrows})
        t_bins  = sorted({(r.tmin,  r.tmax)  for r in binrows})
        if not q2_bins or not t_bins:
            continue
        #endif

        nC = len(q2_bins); nR = len(t_bins)
        fig_w = max(10.0, 3.4 * nC + 2.0)
        fig_h = max(7.0,  2.8 * nR + 2.0)
        fig, axes = plt.subplots(nR, nC, figsize=(fig_w, fig_h),
                                 constrained_layout=True, squeeze=False)
        fig.suptitle(rf"Acceptance $A(\phi)$, $x_{{B}}\in[{xbmin:.2g}, {xbmax:.2g}]$")

        for i_t, tt in enumerate(t_bins):
            for j_q, q2 in enumerate(q2_bins):
                ax = axes[i_t, j_q]

                # Born acceptance + **reco** means
                (c_b, A_b, sA_b, Nrec_b, Ngen_b, Aint_b, sAint_b,
                 Q2b, xb, tb, yb) = acceptance_phi_for_cell(
                    gen_born, rec_born, (xbmin, xbmax), q2, tt, PHI_EDGES
                 )
                # Radiative acceptance + **reco** means
                (c_r, A_r, sA_r, Nrec_r, Ngen_r, Aint_r, sAint_r,
                 Q2r, xr, tr, yr) = acceptance_phi_for_cell(
                    gen_rad, rec_rad, (xbmin, xbmax), q2, tt, PHI_EDGES
                 )

                # Plot: markers only, error bars
                eb1 = ax.errorbar(c_b, A_b, yerr=sA_b, fmt="o", linestyle="none",
                                  markersize=3.0, capsize=2, label="Born")
                eb2 = ax.errorbar(c_r, A_r, yerr=sA_r, fmt="s", linestyle="none",
                                  markersize=3.0, capsize=2, label="Rad")

                ax.set_title(rf"$Q^{{2}}\!\in\![{q2[0]:.2g},{q2[1]:.2g}],\ -t\!\in\![{tt[0]:.2g},{tt[1]:.2g}]$")
                format_axes_phi(ax)

                # Legend text with **reconstructed** mean kinematics
                born_text = (r"$\mathrm{Born\ (reco):}\ "
                             r"\langle Q^{2}\rangle={Q2},\ "
                             r"\langle x_{B}\rangle={xB},\ "
                             r"\langle -t\rangle={t},\ "
                             r"\langle y\rangle={y}$").format(
                                 Q2=fmt_float(Q2b, "{:.2f}"),
                                 xB=fmt_float(xb,  "{:.3f}"),
                                 t =fmt_float(tb,  "{:.3f}"),
                                 y =fmt_float(yb,  "{:.3f}"),
                             )
                rad_text  = (r"$\mathrm{Rad\ (reco):}\ "
                             r"\langle Q^{2}\rangle={Q2},\ "
                             r"\langle x_{B}\rangle={xB},\ "
                             r"\langle -t\rangle={t},\ "
                             r"\langle y\rangle={y}$").format(
                                 Q2=fmt_float(Q2r, "{:.2f}"),
                                 xB=fmt_float(xr,  "{:.3f}"),
                                 t =fmt_float(tr,  "{:.3f}"),
                                 y =fmt_float(yr,  "{:.3f}"),
                             )

                handles = [eb1[0], eb2[0]]
                labels  = [born_text, rad_text]
                ax.legend(handles, labels, frameon=False, loc="best")
            #endfor
        #endfor

        out_pdf = os.path.join(out_dir, f"acc_phi_grid_xB_{xbmin:.2g}_{xbmax:.2g}.pdf")
        fig.savefig(out_pdf)
        print("Saved:", out_pdf)
    #endfor
#endfor

# -------------------- main --------------------
def main():
    # Load arrays once from each file (full statistics)
    print("[info] Loading Born gen/reco...")
    gen_born = load_arrays(GEN_BORN_PATH)
    rec_born = load_arrays(REC_BORN_PATH)

    print("[info] Loading Radiative gen/reco...")
    gen_rad  = load_arrays(GEN_RAD_PATH)
    rec_rad  = load_arrays(REC_RAD_PATH)

    # Output root
    out_root = "output/rad_example"
    os.makedirs(out_root, exist_ok=True)

    # Per-xB canvases comparing Born vs Rad acceptances (legends from RECO means)
    make_acceptance_grids(gen_born, rec_born, gen_rad, rec_rad, out_dir=out_root)
#endfor  # function ends; comment included per your preference

if __name__ == "__main__":
    main()
#endif