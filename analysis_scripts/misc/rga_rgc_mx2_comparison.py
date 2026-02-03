#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rgc_rga_mx2_comparison.py

Mx2 comparison for three samples:
  - RGC Su22 NH3
  - RGC Su22 C
  - RGA Fa18 Inb H2

Updates requested:
  - Histogram/plot x-range: [-1, 4] (GeV^2)
  - Subplot 1: normalize each histogram by its integral over [-1, 4]
               (no charge normalization; CSV kept but commented out)
  - Subplot 2: ratio (NH3_norm / C_norm) with constant fit c in [-0.5, 0.5]
               draw c solid on [-0.5, 0.5], dashed elsewhere
  - Subplot 3: plot (NH3_norm - c * C_norm) and overlay H2_norm

Optional:
  --short : restrict each ROOT tree to only the first 6 run numbers (ascending),
            but keep all events belonging to those runs.

Output:
  output/rgc_rga_mx2_comparison.png

Dependencies:
  python3, uproot, numpy, matplotlib

Example:
  python3 rgc_rga_mx2_comparison.py
  python3 rgc_rga_mx2_comparison.py --short
"""

import os
import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------
# User inputs
# -------------------------------------------------------------------------
NH3_ROOT = "/volatile/clas12/thayward/rgc_enpi+_Mx2_study/rgc_su22_inb_NH3.root"
C_ROOT   = "/volatile/clas12/thayward/rgc_enpi+_Mx2_study/rgc_su22_inb_C.root"
H2_ROOT  = "/volatile/clas12/thayward/rgc_enpi+_Mx2_study/rga_fa18_inb_H2.root"

TREE_NAME = "PhysicsEvents"

# Left here (commented) per request, but not used for normalization anymore:
#RUN_INFO_CSV = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

OUT_PNG = "output/rgc_rga_mx2_comparison.png"

# Histogram settings
MX2_MIN = -1.0
MX2_MAX =  4.0
NBINS   = 250

# Constant fit window in subplot 2
FIT_XMIN = -0.5
FIT_XMAX =  0.5


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
def parse_args(argv):
    short = False
    for a in argv[1:]:
        if a == "--short":
            short = True
        else:
            raise RuntimeError(f"FATAL: Unknown argument: {a}")
        #endif
    #endfor
    return short


def load_tree_arrays(root_path, tree_name):
    """
    Load Mx2 and runnum arrays from a ROOT TTree.
    Returns:
      mx2_np (np.ndarray float64)
      run_np (np.ndarray int64)
    """
    if not os.path.isfile(root_path):
        raise RuntimeError(f"FATAL: ROOT file not found: {root_path}")

    with uproot.open(root_path) as f:
        if tree_name not in f:
            raise RuntimeError(f"FATAL: TTree '{tree_name}' not found in file: {root_path}")

        tree = f[tree_name]

        needed = ["Mx2", "runnum"]
        for br in needed:
            if br not in tree.keys():
                raise RuntimeError(
                    f"FATAL: Branch '{br}' not found in TTree '{tree_name}' in file: {root_path}"
                )
            #endif
        #endfor

        mx2 = tree["Mx2"].array(library="np")
        run = tree["runnum"].array(library="np")

    mx2_np = np.asarray(mx2, dtype=np.float64)
    run_np = np.asarray(run, dtype=np.int64)

    if mx2_np.shape[0] != run_np.shape[0]:
        raise RuntimeError(
            f"FATAL: Array length mismatch in {root_path}: Mx2 has {mx2_np.shape[0]} entries, runnum has {run_np.shape[0]} entries"
        )

    return mx2_np, run_np


def restrict_to_first_n_runs(mx2, run, n_runs, label):
    """
    Restrict arrays to only events belonging to the first n_runs unique run numbers (ascending).
    Keeps all events for those runs.

    Returns:
      mx2_cut, run_cut
    """
    runs_unique = np.unique(run)
    if runs_unique.size == 0:
        raise RuntimeError(f"FATAL: No run numbers found in {label}")

    runs_kept = runs_unique[:n_runs]
    mask = np.isin(run, runs_kept)

    mx2_cut = mx2[mask]
    run_cut = run[mask]

    return mx2_cut, run_cut


def hist_density(mx2_values, bin_edges):
    """
    Build a histogram over [MX2_MIN, MX2_MAX] and normalize by integral.
    We return a "density-like" array where sum(y_i * bin_width_i) = 1.

    Implementation:
      counts -> density = counts / (sum(counts) * bin_width)
    so that the area under the curve is 1 in the chosen range.
    """
    mask = (mx2_values >= MX2_MIN) & (mx2_values <= MX2_MAX)
    mx2_in = mx2_values[mask]

    counts, _ = np.histogram(mx2_in, bins=bin_edges)
    counts = counts.astype(np.float64)

    total = np.sum(counts)
    if total <= 0.0:
        raise RuntimeError("FATAL: Histogram integral is zero in the requested Mx2 range.")

    bin_width = bin_edges[1] - bin_edges[0]
    density = counts / (total * bin_width)

    return density


def safe_ratio(num, den):
    """
    Compute num/den bin-by-bin, returning NaN where den == 0.
    (NaNs won't draw in matplotlib step plots, which is often what you want.)
    """
    out = np.full_like(num, np.nan, dtype=np.float64)
    mask = den != 0.0
    out[mask] = num[mask] / den[mask]
    return out


def fit_constant_in_window(x_centers, y, xmin, xmax):
    """
    Fit y(x) to a constant in [xmin, xmax] using an unweighted mean over valid bins.
    Only bins with finite y are included.

    Returns:
      c (float), n_used (int)
    """
    win = (x_centers >= xmin) & (x_centers <= xmax) & np.isfinite(y)
    ywin = y[win]

    if ywin.size == 0:
        raise RuntimeError("FATAL: No finite ratio bins available in the fit window.")

    c = float(np.mean(ywin))
    return c, int(ywin.size)


def step_plot(ax, bin_edges, y, label, linewidth=1.2, linestyle="-"):
    ax.step(bin_edges[:-1], y, where="post", label=label, linewidth=linewidth, linestyle=linestyle)


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
def main():
    short_mode = parse_args(sys.argv)

    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)

    # Load arrays (only needed branches)
    nh3_mx2, nh3_run = load_tree_arrays(NH3_ROOT, TREE_NAME)
    c_mx2,   c_run   = load_tree_arrays(C_ROOT,   TREE_NAME)
    h2_mx2,  h2_run  = load_tree_arrays(H2_ROOT,  TREE_NAME)

    # Optional short restriction
    if short_mode:
        nh3_mx2, nh3_run = restrict_to_first_n_runs(nh3_mx2, nh3_run, 6, "NH3")
        c_mx2,   c_run   = restrict_to_first_n_runs(c_mx2,   c_run,   6, "C")
        h2_mx2,  h2_run  = restrict_to_first_n_runs(h2_mx2,  h2_run,  6, "H2")
    #endif

    # Fixed binning for all plots
    bin_edges = np.linspace(MX2_MIN, MX2_MAX, NBINS + 1)
    bin_width = bin_edges[1] - bin_edges[0]
    x_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Normalized-by-integral (area = 1) histograms
    nh3_norm = hist_density(nh3_mx2, bin_edges)
    c_norm   = hist_density(c_mx2,   bin_edges)
    h2_norm  = hist_density(h2_mx2,  bin_edges)

    # Subplot 2: ratio and constant fit
    ratio_nh3_over_c = safe_ratio(nh3_norm, c_norm)
    c_fit, n_fit_bins = fit_constant_in_window(x_centers, ratio_nh3_over_c, FIT_XMIN, FIT_XMAX)

    # Subplot 3: subtraction using fitted c
    diff_subtracted = nh3_norm - c_fit * c_norm

    # ---------------------------------------------------------------------
    # Plotting
    # ---------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

    x_label = r"$M_{x}^{2}$ (GeV$^{2}$)"
    y_label = "Arb. units (unit area on [-1, 4])"

    # 1) Overlay of normalized histograms
    ax = axes[0]
    step_plot(ax, bin_edges, nh3_norm, "RGC Su22 NH3")
    step_plot(ax, bin_edges, c_norm,   "RGC Su22 C")
    step_plot(ax, bin_edges, h2_norm,  "RGA Fa18 Inb H2")
    ax.set_title("Mx2 (each spectrum area-normalized)")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(MX2_MIN, MX2_MAX)
    ax.legend(fontsize=10)

    # 2) Ratio + constant fit line (solid in fit window, dashed elsewhere)
    ax = axes[1]
    step_plot(ax, bin_edges, ratio_nh3_over_c, "NH3_norm / C_norm")
    ax.set_title(f"Ratio and constant fit in [{FIT_XMIN}, {FIT_XMAX}]")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Ratio")
    ax.set_xlim(MX2_MIN, MX2_MAX)

    # Draw fitted constant as solid in [-0.5,0.5], dashed elsewhere
    # Left dashed segment
    ax.plot([MX2_MIN, FIT_XMIN], [c_fit, c_fit], linestyle="--", linewidth=1.8, label=None)
    # Solid fit window segment
    ax.plot([FIT_XMIN, FIT_XMAX], [c_fit, c_fit], linestyle="-", linewidth=2.2, label=f"c = {c_fit:.4f}")
    # Right dashed segment
    ax.plot([FIT_XMAX, MX2_MAX], [c_fit, c_fit], linestyle="--", linewidth=1.8, label=None)

    ax.legend(fontsize=10)

    # 3) Subtracted spectrum + H2 overlay
    ax = axes[2]
    step_plot(ax, bin_edges, diff_subtracted, "NH3_norm - c * C_norm")
    step_plot(ax, bin_edges, h2_norm,         "H2_norm")
    ax.axhline(0.0, linewidth=1.0)
    ax.set_title("Subtracted (hydrogen-like) residual vs H2")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(MX2_MIN, MX2_MAX)
    ax.legend(fontsize=10)

    if short_mode:
        fig.suptitle("Mx2 comparison (SHORT mode: first 6 runs only)", fontsize=14)
    #endif

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"Wrote: {OUT_PNG}")
    print(f"Fit constant c = {c_fit:.6f} using {n_fit_bins} bins in [{FIT_XMIN}, {FIT_XMAX}]")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
    #endif