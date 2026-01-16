#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
final_test.py

Usage (5 ROOT inputs total):
  python final_test.py <root1> <root2> <root3> <root4> <root5>

This script makes a 4x6 canvas of Mx2 histograms binned in (xB, -tprime),
with the SAME binning as the normalization workflow:

  xB rows: 4 bins using XB_EDGES
  -tprime cols: 6 bins using TNEG_EDGES

It overlays, in each pad, TWO histograms of the branch "Mx2":
  - root2: labeled "MC"
  - root5: labeled "data"

Histogram settings:
  - Mx2 range: 0.4 to 2.0
  - nbins: 100 (matches the ROOT-style default used in other scripts)
  - Each histogram normalized to TRUE probability density:
        counts / (integral * bin_width)
    (equivalent to np.histogram(..., density=True) within the plotted range)

Output:
  output/Mx2_grid_mc_vs_data.png
"""

import os
import sys
import math

import numpy as np
import uproot
import matplotlib.pyplot as plt


TREE_NAME = "PhysicsEvents"

XB_EDGES = [0.10, 0.25, 0.35, 0.45, 0.60]  # 4 rows
TNEG_EDGES = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25]  # 6 cols in -tprime

MX2_MIN = 0.4
MX2_MAX = 2.0
MX2_NBINS = 100

OUT_PNG = "output/Mx2_grid_mc_vs_data.png"


# Masses (GeV)
MASS_E = 0.000511
MASS_PI = 0.139570
MASS_N = 0.9382720813

# Fixed beam energy (GeV) since there is no runnum in your newer trees.
EB_FIXED = 10.55


def fatal(msg: str, code: int = 1) -> int:
    print(f"FATAL: {msg}", file=sys.stderr)
    return code
#enddef


def require_file(path: str) -> None:
    if path is None or str(path).strip() == "":
        raise RuntimeError("Missing required input path.")
    #endif
    if not os.path.isfile(path):
        raise RuntimeError("File not found: " + str(path))
    #endif
#enddef


def find_bin(val: float, edges: list[float]) -> int:
    for i in range(len(edges) - 1):
        if val >= edges[i] and val < edges[i + 1]:
            return i
        #endif
    #endfor
    return -1
#enddef


def compute_t_scalar_fixed_Eb(Eb: float,
                             e_p: float, e_theta: float, e_phi: float,
                             p_p: float, p_theta: float, p_phi: float) -> float:
    # Reconstruct 4-vectors (assuming massless beam along +z, target at rest)
    E_e = math.sqrt(max(0.0, e_p * e_p + MASS_E * MASS_E))
    se = math.sin(e_theta)
    ce = math.cos(e_theta)
    ex = e_p * se * math.cos(e_phi)
    ey = e_p * se * math.sin(e_phi)
    ez = e_p * ce

    E_pi = math.sqrt(max(0.0, p_p * p_p + MASS_PI * MASS_PI))
    sp = math.sin(p_theta)
    cp = math.cos(p_theta)
    px = p_p * sp * math.cos(p_phi)
    py = p_p * sp * math.sin(p_phi)
    pz = p_p * cp

    dE = (Eb - E_e) - E_pi
    dx = -ex - px
    dy = -ey - py
    dz = (Eb - ez) - pz

    return dE * dE - (dx * dx + dy * dy + dz * dz)
#enddef


def compute_tmin_exact(xB: float, Q2: float) -> float:
    xb_ok = (xB > 0.0 and xB < 1.0)
    if Q2 <= 0.0 or (not xb_ok):
        if xb_ok:
            denom = (1.0 - xB)
            if denom > 0.0:
                return - (MASS_N * xB) * (MASS_N * xB) / denom
            #endif
        #endif
        return 0.0
    #endif

    eps2 = 4.0 * MASS_N * MASS_N * xB * xB / Q2
    root = math.sqrt(1.0 + eps2)
    num = Q2 * (2.0 * (1.0 - xB) * (1.0 - root) + eps2)
    den = 4.0 * xB * (1.0 - xB) + eps2
    if den == 0.0:
        return 0.0
    #endif
    return - num / den
#enddef


def load_arrays_for_grid(root_path: str) -> dict[str, np.ndarray]:
    """
    Load only the branches needed to compute (xB, -tprime) and plot Mx2:
      x, Q2, e_p,e_theta,e_phi, p_p,p_theta,p_phi, Mx2

    Returns dict of numpy arrays with the same length.
    """
    with uproot.open(root_path) as f:
        if TREE_NAME not in f:
            raise KeyError(f"TTree '{TREE_NAME}' not found in '{root_path}'.")
        #endif
        t = f[TREE_NAME]

        needed = ["x", "Q2", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "Mx2"]
        missing = [b for b in needed if b not in t.keys()]
        if len(missing) > 0:
            raise KeyError(f"Missing required branches in '{root_path}': {', '.join(missing)}")
        #endif

        arrs = t.arrays(needed, library="np")
    #endwith

    out = {}
    for k in needed:
        out[k] = np.asarray(arrs[k], dtype=np.float64)
    #endfor
    return out
#enddef


def fill_grid_mx2(arrs: dict[str, np.ndarray]) -> list[list[np.ndarray]]:
    """
    Returns a grid [r][c] where each entry is a 1D numpy array of Mx2 values
    that fell into that (xB, -tprime) bin.

    NOTE: This uses the same physics recipe as the ROOT script:
      t = compute_t_scalar(...)
      tmin = compute_tmin_exact(xB, Q2)
      tprime = t - tmin
      -tprime = -(t - tmin)
    """
    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    buckets = [[[] for _ in range(ncols)] for _ in range(nrows)]

    xB = arrs["x"]
    Q2 = arrs["Q2"]
    e_p = arrs["e_p"]
    e_theta = arrs["e_theta"]
    e_phi = arrs["e_phi"]
    p_p = arrs["p_p"]
    p_theta = arrs["p_theta"]
    p_phi = arrs["p_phi"]
    Mx2 = arrs["Mx2"]

    n = len(Mx2)

    for i in range(n):
        xb = float(xB[i])
        r = find_bin(xb, XB_EDGES)
        if r < 0:
            continue
        #endif

        q2 = float(Q2[i])

        t_val = compute_t_scalar_fixed_Eb(
            EB_FIXED,
            float(e_p[i]), float(e_theta[i]), float(e_phi[i]),
            float(p_p[i]), float(p_theta[i]), float(p_phi[i]),
        )
        tmin_val = compute_tmin_exact(xb, q2)
        tprime = t_val - tmin_val
        tneg = -tprime

        c = find_bin(tneg, TNEG_EDGES)
        if c < 0:
            continue
        #endif

        v = float(Mx2[i])
        if not math.isfinite(v):
            continue
        #endif

        buckets[r][c].append(v)
    #endfor

    # Convert to numpy arrays
    out = []
    for r in range(nrows):
        row = []
        for c in range(ncols):
            row.append(np.asarray(buckets[r][c], dtype=np.float64))
        #endfor
        out.append(row)
    #endfor
    return out
#enddef


def hist_density(arr: np.ndarray, nbins: int, xmin: float, xmax: float) -> tuple[np.ndarray, np.ndarray]:
    in_range = (arr >= xmin) & (arr <= xmax)
    x = arr[in_range]

    density, edges = np.histogram(x, bins=nbins, range=(xmin, xmax), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, density
#enddef


def make_grid_plot(mc_grid: list[list[np.ndarray]],
                   data_grid: list[list[np.ndarray]],
                   out_png: str) -> None:
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 2.8 * nrows), sharex=True, sharey=False)

    # Ensure axes is 2D
    if nrows == 1 and ncols == 1:
        axes = np.asarray([[axes]])
    elif nrows == 1:
        axes = np.asarray([axes])
    elif ncols == 1:
        axes = np.asarray([[a] for a in axes])
    #endif

    for r in range(nrows):
        xb_lo = XB_EDGES[r]
        xb_hi = XB_EDGES[r + 1]
        for c in range(ncols):
            t_lo = TNEG_EDGES[c]
            t_hi = TNEG_EDGES[c + 1]

            ax = axes[r][c]

            mc_arr = mc_grid[r][c]
            data_arr = data_grid[r][c]

            # Plot MC first, then data, matching your request.
            x_mc, y_mc = hist_density(mc_arr, MX2_NBINS, MX2_MIN, MX2_MAX)
            x_da, y_da = hist_density(data_arr, MX2_NBINS, MX2_MIN, MX2_MAX)

            ax.step(x_mc, y_mc, where="mid", linewidth=1.2, label="MC")
            ax.step(x_da, y_da, where="mid", linewidth=1.2, label="data")

            ax.grid(True)

            ax.set_xlim(MX2_MIN, MX2_MAX)

            # Title inside each pad
            ax.set_title(f"xB [{xb_lo:.2f}, {xb_hi:.2f})  -tprime [{t_lo:.2f}, {t_hi:.2f})", fontsize=9)

            # Only put axis labels on outer edges to reduce clutter
            if r == nrows - 1:
                ax.set_xlabel("Mx2 (GeV^2)")
            #endif
            if c == 0:
                ax.set_ylabel("Prob. density (1/GeV^2)")
            #endif

            # Legend: keep it small; only show on first pad to avoid clutter
            if r == 0 and c == 0:
                ax.legend(loc="upper right", frameon=True, fontsize=10)
            #endif
        #endfor
    #endfor

    fig.suptitle("Mx2 in (xB, -tprime) bins: root2 (MC) vs root5 (data)", fontsize=14)
    plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    plt.savefig(out_png, dpi=200)
    plt.close(fig)
#enddef


def main() -> int:
    if len(sys.argv) != 6:
        return fatal("expected 5 ROOT input files: <root1> <root2> <root3> <root4> <root5>", 2)
    #endif

    # root2 is MC, root5 is data per your request
    root2 = sys.argv[2]
    root5 = sys.argv[5]

    try:
        require_file(root2)
        require_file(root5)

        mc_arrs = load_arrays_for_grid(root2)
        data_arrs = load_arrays_for_grid(root5)

        mc_grid = fill_grid_mx2(mc_arrs)
        data_grid = fill_grid_mx2(data_arrs)

        make_grid_plot(mc_grid, data_grid, OUT_PNG)

    except Exception as e:
        return fatal(str(e), 3)
    #endtry

    print("Wrote: " + OUT_PNG)
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif