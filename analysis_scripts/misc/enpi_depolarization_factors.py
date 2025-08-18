#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
enpi_olarization_factors.py

2x2 figure: each subplot is an x_B bin, showing the mean of B/A, C/A,
V/A, and W/A vs -t with statistical error bars (SEM).

Input:
  ROOT file: /work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+_2.root
  Tree:      PhysicsEvents

Kinematic cuts (always enforced):
  Q2 > 1
  W > 1
  y < 0.75
  fiducial_status >= 111
  0.81 < Mx2 < 1.00

Binning:
  x_B bins: [0.10, 0.25], [0.25, 0.35], [0.35, 0.45], [0.45, 0.60]
  -t (GeV^2) edges: [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25]
  (Note: t is stored negative in the tree; we use tpos = -t)

Output:
  output/enpi+/olarization_factors.pdf
"""

from pathlib import Path
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
T_POS_EDGES = np.array([
    0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65,
    0.75, 0.85, 0.95, 1.05, 1.15, 1.25
], dtype=float)

# Output
OUT_PATH = Path("output/enpi+/olarization_factors.pdf")
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# Labels (LaTeX-friendly; ASCII only)
AXIS_LABEL_X = r"$-t\ (\mathrm{GeV}^{2})$"
AXIS_LABEL_Y = "olarization factor"

# Fixed axes
X_LIMS = (0.0, 1.3)
Y_LIMS = (0.0, 2.0)

# Chunk size for iterate: integer number of entries per chunk
ITER_STEP = 1_000_000

# -----------------------
# Helpers
# -----------------------
def unbiased_sem(count, sum_x, sum_x2):
    """
    Compute mean and standard error of the mean (SEM) from running sums.
    Returns (mean, sem). If count == 0, returns (nan, nan). If count < 2, sem is nan.
    """
    if count == 0:
        return np.nan, np.nan  # endif
    mean = sum_x / count
    if count < 2:
        return mean, np.nan  # endif
    var_unbiased = (sum_x2 - (sum_x * sum_x) / count) / (count - 1)
    if var_unbiased < 0.0:
        var_unbiased = 0.0
    sem = np.sqrt(var_unbiased / count)
    return mean, sem  # endif
#endfor

def update_acc(acc_tuple, values):
    """
    Update (count, sum, sumsq) accumulators with finite values.
    acc_tuple is (count_ref, sum_ref, sumsq_ref) where each is a 1-element view.
    """
    if values.size == 0:
        return  # endif
    acc_tuple[0][...] += values.size
    acc_tuple[1][...] += values.sum()
    acc_tuple[2][...] += np.dot(values, values)
#endfor

# -----------------------
# Main
# -----------------------
def main():
    n_xb = len(XB_EDGES) - 1
    n_t  = len(T_POS_EDGES) - 1

    ratio_keys = ["B/A", "C/A", "V/A", "W/A"]
    counts = {k: np.zeros((n_xb, n_t), dtype=np.int64) for k in ratio_keys}
    sums   = {k: np.zeros((n_xb, n_t), dtype=np.float64) for k in ratio_keys}
    sums2  = {k: np.zeros((n_xb, n_t), dtype=np.float64) for k in ratio_keys}

    needed = [
        "x", "t", "Q2", "W", "y", "fiducial_status", "Mx2",
        "A", "B", "C", "V", "W"
    ]

    tree_spec = f"{ROOT_PATH}:{TREE_NAME}"
    for arrays in uproot.iterate(tree_spec, filter_name=needed, step_size=ITER_STEP, library="np"):
        x   = arrays["x"]
        t   = arrays["t"]
        Q2  = arrays["Q2"]
        W   = arrays["W"]
        y   = arrays["y"]
        fid = arrays["fiducial_status"]
        Mx2 = arrays["Mx2"]

        A = arrays["A"]
        B = arrays["B"]
        C = arrays["C"]
        V = arrays["V"]
        W = arrays["W"]

        # Global kinematic cuts (strict inequalities for Mx2)
        base_mask = (
            (Q2 > 1.0) &
            (W  > 1.0) &
            (y  < 0.75) &
            (fid >= 111) &
            (Mx2 > 0.81) & (Mx2 < 1.00)
        )

        tpos = -t

        finite_mask = (
            np.isfinite(x) & np.isfinite(tpos) &
            np.isfinite(A) & np.isfinite(B) &
            np.isfinite(C) & np.isfinite(V) &
            np.isfinite(W)
        )

        mask = base_mask & finite_mask
        if not np.any(mask):
            continue  # endif

        x    = x[mask]
        tpos = tpos[mask]
        A    = A[mask]
        B    = B[mask]
        C    = C[mask]
        V    = V[mask]
        Wv   = W[mask]

        # Avoid division by zero
        goodA = (A != 0.0) & np.isfinite(A)
        if not np.any(goodA):
            continue  # endif

        x    = x[goodA]
        tpos = tpos[goodA]
        A    = A[goodA]
        B    = B[goodA]
        C    = C[goodA]
        V    = V[goodA]
        Wv   = Wv[goodA]

        rB = B / A
        rC = C / A
        rV = V / A
        rW = Wv / A

        # Bin by x_B then by -t
        for ix in range(n_xb):
            xb_lo = XB_EDGES[ix]
            xb_hi = XB_EDGES[ix + 1]
            if ix == n_xb - 1:
                mx = (x >= xb_lo) & (x <= xb_hi)
            else:
                mx = (x >= xb_lo) & (x < xb_hi)
            #endif
            if not np.any(mx):
                continue  # endif

            t_sel   = tpos[mx]
            rB_sel  = rB[mx]
            rC_sel  = rC[mx]
            rV_sel  = rV[mx]
            rW_sel  = rW[mx]

            for it in range(n_t):
                t_lo = T_POS_EDGES[it]
                t_hi = T_POS_EDGES[it + 1]
                if it == n_t - 1:
                    mt = (t_sel >= t_lo) & (t_sel <= t_hi)
                else:
                    mt = (t_sel >= t_lo) & (t_sel < t_hi)
                #endif
                if not np.any(mt):
                    continue  # endif

                vals_B = rB_sel[mt]
                vals_C = rC_sel[mt]
                vals_V = rV_sel[mt]
                vals_W = rW_sel[mt]

                cB, sB, s2B = counts["B/A"], sums["B/A"], sums2["B/A"]
                cC, sC, s2C = counts["C/A"], sums["C/A"], sums2["C/A"]
                cV, sV, s2V = counts["V/A"], sums["V/A"], sums2["V/A"]
                cW, sW, s2W = counts["W/A"], sums["W/A"], sums2["W/A"]

                update_acc((cB[ix, it:it+1], sB[ix, it:it+1], s2B[ix, it:it+1]), vals_B)
                update_acc((cC[ix, it:it+1], sC[ix, it:it+1], s2C[ix, it:it+1]), vals_C)
                update_acc((cV[ix, it:it+1], sV[ix, it:it+1], s2V[ix, it:it+1]), vals_V)
                update_acc((cW[ix, it:it+1], sW[ix, it:it+1], s2W[ix, it:it+1]), vals_W)
            #endfor
        #endfor
    #endfor

    # Reduce to means and SEMs
    means = {k: np.full((n_xb, n_t), np.nan, dtype=float) for k in ratio_keys}
    semes = {k: np.full((n_xb, n_t), np.nan, dtype=float) for k in ratio_keys}

    for ix in range(n_xb):
        for it in range(n_t):
            for rkey in ratio_keys:
                c = counts[rkey][ix, it]
                s = sums[rkey][ix, it]
                s2 = sums2[rkey][ix, it]
                mu, se = unbiased_sem(c, s, s2)
                means[rkey][ix, it] = mu
                semes[rkey][ix, it] = se
            #endfor
        #endfor
    #endfor

    # Plotting
    t_centers = 0.5 * (T_POS_EDGES[:-1] + T_POS_EDGES[1:])
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    axes = axes.reshape(-1)

    # marker-only (no connecting lines)
    style = {
        "B/A": dict(fmt="o", linestyle="none", markersize=4, capsize=2),
        "C/A": dict(fmt="s", linestyle="none", markersize=4, capsize=2),
        "V/A": dict(fmt="^", linestyle="none", markersize=4, capsize=2),
        "W/A": dict(fmt="D", linestyle="none", markersize=4, capsize=2),
    }

    for ix in range(n_xb):
        ax = axes[ix]
        xb_lo, xb_hi = XB_EDGES[ix], XB_EDGES[ix + 1]
        for rkey in ratio_keys:
            y    = means[rkey][ix, :]
            yerr = semes[rkey][ix, :]
            ax.errorbar(t_centers, y, yerr=yerr, label=rkey, **style[rkey])
        #endfor

        ax.set_title(r"$x_{B}$ in [{:.2f}, {:.2f}]".format(xb_lo, xb_hi))
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
        ax.legend(loc="upper right", frameon=True, fontsize=9)
        ax.set_xlabel(AXIS_LABEL_X)
        ax.set_ylabel(AXIS_LABEL_Y)
        ax.set_xlim(*X_LIMS)
        ax.set_ylim(*Y_LIMS)
    #endfor

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=300)
    print(f"[OK] Saved: {OUT_PATH}")
#endfor

if __name__ == "__main__":
    main()
#endif