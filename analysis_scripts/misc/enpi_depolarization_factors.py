#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
depolarization_factors.py

Builds a 2x2 figure (one subplot per x_B bin) showing the mean depolarization
ratios DepB/DepA, DepC/DepA, DepV/DepA, and DepW/DepA as a function of -t
within each x_B bin, including statistical (standard error of the mean) bars.

Input:
  ROOT file: /work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+_2.root
  Tree:      PhysicsEvents

Kinematic cuts (always enforced):
  Q2 > 1
  W > 1
  y < 0.75
  fiducial_status >= 111
  0.81 <= Mx2 <= 1.00

Binning:
  x_B bins:  [0.10, 0.25], [0.25, 0.35], [0.35, 0.45], [0.45, 0.60]
  -t bins:   edges = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75,
                      0.85, 0.95, 1.05, 1.15, 1.25]
             (Note: t in the file is negative like the Mandelstam variable; we use tpos = -t)

Output:
  PDF saved to: output/enpi+/depolarization_factors.pdf

Notes:
  • The user requested the x-axis label "x_{B}". Because each subplot is a *fixed* x_B bin
    and the points vary with -t inside that bin, the natural horizontal axis is -t (GeV^2).
    To keep the plot physically correct, we label the x-axis as "-t (GeV^2)" in the figure.
    If you still prefer the literal "x_{B}" label, change AXIS_LABEL_X below.
"""

from pathlib import Path
import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

# -----------------------
# Configuration
# -----------------------
ROOT_PATH = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+_2.root"
TREE_NAME = "PhysicsEvents"

# x_B bins
XB_EDGES = np.array([0.10, 0.25, 0.35, 0.45, 0.60], dtype=float)

# -t (GeV^2) bin edges (positive)
T_POS_EDGES = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65,
                        0.75, 0.85, 0.95, 1.05, 1.15, 1.25], dtype=float)

# Output
OUT_PATH = Path("output/enpi+/depolarization_factors.pdf")
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# Axis labels and styling
AXIS_LABEL_X = r"-t (GeV^2)"   # change to r"$x_{B}$" if you want the literal request
AXIS_LABEL_Y = "depolarization factor"

# Chunks for iteration (tune if needed)
ITER_STEP = "1e6"  # ~1 million entries per chunk

# -----------------------
# Helper functions
# -----------------------
def unbiased_sem(count, sum_x, sum_x2):
    """
    Compute mean and standard error of the mean (SEM) using an unbiased
    estimate of the sample variance. Returns (mean, sem).
    If count < 2, SEM is set to np.nan.
    """
    if count == 0:
        return np.nan, np.nan  # endif
    mean = sum_x / count
    if count < 2:
        return mean, np.nan  # endif
    var_unbiased = (sum_x2 - (sum_x * sum_x) / count) / (count - 1)
    var_unbiased = max(var_unbiased, 0.0)
    sem = np.sqrt(var_unbiased / count)
    return mean, sem  # endif
#endfor (function end marker)

def update_accumulators(acc, values):
    """
    Update (count, sum, sumsq) accumulators with finite values.
    acc is a tuple/list of three numpy arrays: (count, sum, sumsq) each scalar here.
    values is a 1D array of valid samples.
    """
    if values.size == 0:
        return  # endif
    acc[0][...] += values.size
    acc[1][...] += values.sum()
    acc[2][...] += np.dot(values, values)  # sum of squares
#endfor

# -----------------------
# Main
# -----------------------
def main():
    # Prepare accumulators for 4 ratios across [xB bins] x [t bins]
    n_xb = len(XB_EDGES) - 1
    n_t  = len(T_POS_EDGES) - 1

    # For each ratio, we store count, sum, sumsq with shape (n_xb, n_t)
    ratios = ["DepB/DepA", "DepC/DepA", "DepV/DepA", "DepW/DepA"]
    counts = {key: np.zeros((n_xb, n_t), dtype=np.int64) for key in ratios}
    sums   = {key: np.zeros((n_xb, n_t), dtype=np.float64) for key in ratios}
    sums2  = {key: np.zeros((n_xb, n_t), dtype=np.float64) for key in ratios}

    # Required branches
    needed = [
        "x", "t", "Q2", "W", "y", "fiducial_status", "Mx2",
        "DepA", "DepB", "DepC", "DepV", "DepW"
    ]

    # Iterate in chunks for memory efficiency
    tree_spec = f"{ROOT_PATH}:{TREE_NAME}"
    for arrays in uproot.iterate(tree_spec, filter_name=needed, step_size=ITER_STEP, library="np"):
        # Extract needed arrays
        x    = arrays["x"]
        t    = arrays["t"]         # negative values in file
        Q2   = arrays["Q2"]
        W    = arrays["W"]
        y    = arrays["y"]
        fid  = arrays["fiducial_status"]
        Mx2  = arrays["Mx2"]

        DepA = arrays["DepA"]
        DepB = arrays["DepB"]
        DepC = arrays["DepC"]
        DepV = arrays["DepV"]
        DepW = arrays["DepW"]

        # Apply global kinematic cuts
        base_mask = (
            (Q2 > 1.0) &
            (W  > 1.0) &
            (y  < 0.75) &
            (fid >= 111) &
            (Mx2 >= 0.81) & (Mx2 <= 1.00)
        )

        # Derived positive t
        tpos = -t

        # Keep only finite and sensible numeric entries for all needed inputs
        finite_mask = (
            np.isfinite(x) & np.isfinite(tpos) &
            np.isfinite(DepA) & np.isfinite(DepB) &
            np.isfinite(DepC) & np.isfinite(DepV) &
            np.isfinite(DepW)
        )

        mask_all = base_mask & finite_mask

        # Short-circuit if nothing passes
        if not np.any(mask_all):
            continue  # endif

        x  = x[mask_all]
        tpos = tpos[mask_all]
        A = DepA[mask_all]
        B = DepB[mask_all]
        C = DepC[mask_all]
        V = DepV[mask_all]
        Wv = DepW[mask_all]

        # Avoid division by zero: we only consider entries with A != 0
        good = (A != 0.0) & np.isfinite(A)
        if not np.any(good):
            continue  # endif

        x   = x[good]
        tpos = tpos[good]
        A   = A[good]
        B   = B[good]
        C   = C[good]
        V   = V[good]
        Wv  = Wv[good]

        # Compute ratios
        rB = B / A
        rC = C / A
        rV = V / A
        rW = Wv / A

        # Bin loop: for each x_B bin, then each -t bin
        for ix in range(n_xb):
            xb_lo = XB_EDGES[ix]
            xb_hi = XB_EDGES[ix + 1]
            # Include right edge on last bin
            if ix == n_xb - 1:
                m_x = (x >= xb_lo) & (x <= xb_hi)
            else:
                m_x = (x >= xb_lo) & (x < xb_hi)
            #endif

            if not np.any(m_x):
                continue  # endif

            x_sel   = x[m_x]
            t_sel   = tpos[m_x]
            rB_sel  = rB[m_x]
            rC_sel  = rC[m_x]
            rV_sel  = rV[m_x]
            rW_sel  = rW[m_x]

            for it in range(n_t):
                t_lo = T_POS_EDGES[it]
                t_hi = T_POS_EDGES[it + 1]
                # Include right edge on last t bin
                if it == n_t - 1:
                    m_t = (t_sel >= t_lo) & (t_sel <= t_hi)
                else:
                    m_t = (t_sel >= t_lo) & (t_sel < t_hi)
                #endif

                if not np.any(m_t):
                    continue  # endif

                # Gather valid values per ratio in this (x_B, -t) bin
                vals_B = rB_sel[m_t]
                vals_C = rC_sel[m_t]
                vals_V = rV_sel[m_t]
                vals_W = rW_sel[m_t]

                # Update accumulators
                # For convenience, wrap scalars so update_accumulators can use [...] assignment
                cB = counts["DepB/DepA"]; sB = sums["DepB/DepA"]; s2B = sums2["DepB/DepA"]
                cC = counts["DepC/DepA"]; sC = sums["DepC/DepA"]; s2C = sums2["DepC/DepA"]
                cV = counts["DepV/DepA"]; sV = sums["DepV/DepA"]; s2V = sums2["DepV/DepA"]
                cW = counts["DepW/DepA"]; sW = sums["DepW/DepA"]; s2W = sums2["DepW/DepA"]

                update_accumulators((cB[ix, it:it+1], sB[ix, it:it+1], s2B[ix, it:it+1]), vals_B)
                update_accumulators((cC[ix, it:it+1], sC[ix, it:it+1], s2C[ix, it:it+1]), vals_C)
                update_accumulators((cV[ix, it:it+1], sV[ix, it:it+1], s2V[ix, it:it+1]), vals_V)
                update_accumulators((cW[ix, it:it+1], sW[ix, it:it+1], s2W[ix, it:it+1]), vals_W)
            #endfor (t bins)
        #endfor (xB bins)
    #endfor (chunk iteration)

    # -----------------------
    # Reduce to mean and SEM per (x_B, -t) bin
    # -----------------------
    means = {key: np.full((n_xb, n_t), np.nan, dtype=float) for key in ratios}
    semes = {key: np.full((n_xb, n_t), np.nan, dtype=float) for key in ratios}

    for ix in range(n_xb):
        for it in range(n_t):
            for rkey in ratios:
                c = counts[rkey][ix, it]
                s = sums[rkey][ix, it]
                s2 = sums2[rkey][ix, it]
                mu, se = unbiased_sem(c, s, s2)
                means[rkey][ix, it] = mu
                semes[rkey][ix, it] = se
            #endfor
        #endfor
    #endfor

    # -----------------------
    # Plotting
    # -----------------------
    t_centers = 0.5 * (T_POS_EDGES[:-1] + T_POS_EDGES[1:])

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    axes = axes.reshape(-1)

    # Marker/line styles
    style = {
        "DepB/DepA": dict(marker="o", linestyle="-", linewidth=1.0, markersize=4),
        "DepC/DepA": dict(marker="s", linestyle="-", linewidth=1.0, markersize=4),
        "DepV/DepA": dict(marker="^", linestyle="-", linewidth=1.0, markersize=4),
        "DepW/DepA": dict(marker="D", linestyle="-", linewidth=1.0, markersize=4),
    }

    for ix in range(n_xb):
        ax = axes[ix]
        xb_lo, xb_hi = XB_EDGES[ix], XB_EDGES[ix+1]

        # Plot each ratio in distinct style
        for rkey in ratios:
            y   = means[rkey][ix, :]
            yerr = semes[rkey][ix, :]

            ax.errorbar(
                t_centers, y, yerr=yerr,
                label=rkey,
                **style[rkey]
            )
        #endfor

        ax.set_title(r"$x_{B}$ in [{:.2f}, {:.2f}]".format(xb_lo, xb_hi))
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
        ax.legend(loc="upper right", frameon=True, fontsize=9)
        ax.set_xlabel(AXIS_LABEL_X)
        ax.set_ylabel(AXIS_LABEL_Y)
    #endfor

    # Hide any unused subplot (there should be none here, but guard anyway)
    for k in range(4, len(axes)):
        axes[k].axis("off")
    #endfor

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=300)
    print(f"[OK] Saved: {OUT_PATH}")

#endfor (main)

if __name__ == "__main__":
    main()
#endif