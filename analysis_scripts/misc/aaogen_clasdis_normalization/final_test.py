#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
final_test.py

Usage (5 ROOT inputs total):
  python final_test.py <root1> <root2> <root3> <root4> <root5>

Original behavior (PRESERVED):
Canvas 1 -> output/Mx2_comparison.png
  - root1: "At Rest"
  - root2: "Simulated Fermi Motion"

Canvas 2 -> output/Mx2_data_comparison.png
  - root3: "RGA"
  - root4: "RGA simulated Fermi Motion"
  - root5: "RGC"

Both canvases:
  - Histogram range: 0.4 to 2.0
  - TRUE probability density normalization:
        counts / (integral * bin_width)
    exactly what np.histogram(..., density=True) produces.
  - X-axis label: "M_{x}^{2} (GeV^{2})"

New behavior (ADDED, no regression):
Canvas 3 -> output/Mx2_grid_mc_vs_data.png
  - 4x6 grid in (xB, -tprime) bins (XB_EDGES x TNEG_EDGES)
  - Overlays root2 ("MC") vs root5 ("data") for the branch "Mx2"
  - Binning uses recomputed tprime from reconstructed kinematics with fixed Eb=10.55
  - Mx2 histogram range: 0.4 to 2.0
  - TRUE density normalization in each pad (shape-only)
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
MX2_NBINS_OVERLAY = 320
MX2_NBINS_GRID = 100

OUT_OVERLAY_1 = "output/Mx2_comparison.png"
OUT_OVERLAY_2 = "output/Mx2_data_comparison.png"
OUT_GRID = "output/Mx2_grid_mc_vs_data.png"

# Masses (GeV)
MASS_E = 0.000511
MASS_PI = 0.139570
MASS_N = 0.9382720813

# Fixed beam energy (GeV) for binning in (xB, -tprime) in the grid plot
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


def load_branch_array(root_path: str, tree_name: str, branch_name: str) -> np.ndarray:
    if not os.path.isfile(root_path):
        raise FileNotFoundError(f"input ROOT file not found: {root_path}")
    #endif

    with uproot.open(root_path) as f:
        if tree_name not in f:
            keys = [k for k in f.keys()]
            raise KeyError(f"TTree '{tree_name}' not found in file '{root_path}'. Keys: {keys}")
        #endif

        t = f[tree_name]
        if branch_name not in t.keys():
            brs = [k for k in t.keys()]
            raise KeyError(
                f"branch '{branch_name}' not found in TTree '{tree_name}' for file '{root_path}'. "
                f"Example branches: {brs[:50]}{' ...' if len(brs) > 50 else ''}"
            )
        #endif

        arr = t[branch_name].array(library="np")
    #endwith

    arr = np.asarray(arr, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return arr
#enddef


def hist_density(arr: np.ndarray, nbins: int, xmin: float, xmax: float) -> tuple[np.ndarray, np.ndarray]:
    in_range = (arr >= xmin) & (arr <= xmax)
    x = arr[in_range]

    density, edges = np.histogram(x, bins=nbins, range=(xmin, xmax), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, density
#enddef


def plot_overlay(out_png: str, files: list[str], labels: list[str], title: str,
                 nbins: int, xmin: float, xmax: float) -> None:
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    fig = plt.figure(figsize=(9, 6))

    tree_name = TREE_NAME
    branch_name = "Mx2"

    for path, lab in zip(files, labels):
        arr = load_branch_array(path, tree_name, branch_name)
        x, y = hist_density(arr, nbins, xmin, xmax)
        plt.step(x, y, where="mid", linewidth=1.5, label=lab)
    #endfor

    plt.xlabel("M_{x}^{2} (GeV^{2})")
    plt.ylabel("Probability density (1/GeV^{2})")
    plt.title(title)
    plt.legend(frameon=True)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)
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


def make_grid_plot(mc_grid: list[list[np.ndarray]],
                   data_grid: list[list[np.ndarray]],
                   out_png: str) -> None:
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 2.8 * nrows), sharex=True, sharey=False)

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

            x_mc, y_mc = hist_density(mc_arr, MX2_NBINS_GRID, MX2_MIN, MX2_MAX)
            x_da, y_da = hist_density(data_arr, MX2_NBINS_GRID, MX2_MIN, MX2_MAX)

            ax.step(x_mc, y_mc, where="mid", linewidth=1.2, label="MC")
            ax.step(x_da, y_da, where="mid", linewidth=1.2, label="data")

            ax.grid(True)
            ax.set_xlim(MX2_MIN, MX2_MAX)
            ax.set_title(f"xB [{xb_lo:.2f}, {xb_hi:.2f})  -tprime [{t_lo:.2f}, {t_hi:.2f})", fontsize=9)

            if r == nrows - 1:
                ax.set_xlabel("Mx2 (GeV^2)")
            #endif
            if c == 0:
                ax.set_ylabel("Prob. density (1/GeV^2)")
            #endif

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

    root1 = sys.argv[1]
    root2 = sys.argv[2]
    root3 = sys.argv[3]
    root4 = sys.argv[4]
    root5 = sys.argv[5]

    xmin = 0.4
    xmax = 2.0

    try:
        # Preserve original overlays
        plot_overlay(
            out_png=OUT_OVERLAY_1,
            files=[root1, root2],
            labels=["At Rest", "Simulated Fermi Motion"],
            title="Mx2 comparison",
            nbins=MX2_NBINS_OVERLAY,
            xmin=xmin,
            xmax=xmax,
        )

        plot_overlay(
            out_png=OUT_OVERLAY_2,
            files=[root3, root4, root5],
            labels=["RGA", "RGA simulated Fermi Motion", "RGC"],
            title="Mx2 data comparison",
            nbins=MX2_NBINS_OVERLAY,
            xmin=xmin,
            xmax=xmax,
        )

        # Add new 4x6 grid: root2 vs root5
        require_file(root2)
        require_file(root5)
        mc_arrs = load_arrays_for_grid(root2)
        data_arrs = load_arrays_for_grid(root5)
        mc_grid = fill_grid_mx2(mc_arrs)
        data_grid = fill_grid_mx2(data_arrs)
        make_grid_plot(mc_grid, data_grid, OUT_GRID)

    except Exception as e:
        return fatal(str(e), 3)
    #endtry

    print("Wrote: " + OUT_OVERLAY_1)
    print("Wrote: " + OUT_OVERLAY_2)
    print("Wrote: " + OUT_GRID)
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif