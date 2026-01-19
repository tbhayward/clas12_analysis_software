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
  - Binning uses recomputed -tprime from reconstructed kinematics with fixed Eb=10.55
  - Mx2 histogram range: 0.4 to 2.0
  - TRUE density normalization in each pad (shape-only)

Additional NEW grids (ADDED):
  output/Q2_grid_mc_vs_data.png       : per-pad Q2 shapes
  output/xB_grid_mc_vs_data.png       : per-pad xB shapes
  output/tprime_grid_mc_vs_data.png   : per-pad (-tprime) shapes (computed)
  output/phi_grid_mc_vs_data.png      : per-pad phi shapes (WRAPPED to [0, 2*pi) in radians)
All are density-normalized per pad and use the SAME pad binning (xB, -tprime).

NOTE (requested update):
  This version "fixes all the latex labels" by using Matplotlib mathtext consistently:
    - use $...$ wrappers
    - use \mathrm{} for units
    - use \phi, x_B, Q^2, t^\prime, etc.
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

# ----------------------------
# Overlay settings (unchanged)
# ----------------------------
MX2_MIN = 0.4
MX2_MAX = 2.0
MX2_NBINS_OVERLAY = 320

OUT_OVERLAY_1 = "output/Mx2_comparison.png"
OUT_OVERLAY_2 = "output/Mx2_data_comparison.png"

# ----------------------------
# Grid settings (existing Mx2 grid)
# ----------------------------
MX2_NBINS_GRID = 100
OUT_GRID_MX2 = "output/Mx2_grid_mc_vs_data.png"

# ----------------------------
# Added grid settings (new variables)
# You can tune these ranges/bins as needed.
# ----------------------------
Q2_MIN = 0.0
Q2_MAX = 10.0
Q2_NBINS_GRID = 80
OUT_GRID_Q2 = "output/Q2_grid_mc_vs_data.png"

XB_MIN = 0.10
XB_MAX = 0.60
XB_NBINS_GRID = 80
OUT_GRID_XB = "output/xB_grid_mc_vs_data.png"

TNEG_MIN = 0.0
TNEG_MAX = 2.0
TNEG_NBINS_GRID = 80
OUT_GRID_TNEG = "output/tprime_grid_mc_vs_data.png"

# phi is in radians (wrapped to [0, 2*pi))
PHI_MIN = 0.0
PHI_MAX = 2.0 * math.pi
PHI_NBINS_GRID = 90
OUT_GRID_PHI = "output/phi_grid_mc_vs_data.png"

# Masses (GeV)
MASS_E = 0.000511
MASS_PI = 0.139570
MASS_N = 0.9382720813

# Fixed beam energy (GeV) for binning in (xB, -tprime) in the grid plots
EB_FIXED = 10.55


# ----------------------------
# Label helpers (Matplotlib mathtext)
# ----------------------------
LABEL_MX2 = r"$M_{x}^{2}\;(\mathrm{GeV}^{2})$"
LABEL_Q2 = r"$Q^{2}\;(\mathrm{GeV}^{2})$"
LABEL_XB = r"$x_{B}$"
LABEL_TPRIME_NEG = r"$-t^{\prime}\;(\mathrm{GeV}^{2})$"
LABEL_PHI = r"$\phi\;(\mathrm{rad})$"

YLABEL_DENSITY_GEV2 = r"$\mathrm{Probability\ density}\;(\mathrm{GeV}^{-2})$"
YLABEL_DENSITY_RAD = r"$\mathrm{Probability\ density}\;(\mathrm{rad}^{-1})$"
YLABEL_DENSITY_DIMLESS = r"$\mathrm{Probability\ density}$"

TITLE_BIN_FMT = r"$x_{B}\in[{xb_lo:.2f},{xb_hi:.2f})\ ,\ -t^\prime\in[{t_lo:.2f},{t_hi:.2f})\ (\mathrm{{GeV}}^{{2}})$"


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

    plt.xlabel(LABEL_MX2)
    plt.ylabel(YLABEL_DENSITY_GEV2)
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


def wrap_phi_rad(phi_val: float) -> float:
    # Force into [0, 2*pi)
    p = float(phi_val)
    if not math.isfinite(p):
        return float("nan")
    #endif
    twopi = 2.0 * math.pi
    p = p % twopi
    if p < 0.0:
        p += twopi
    #endif
    return p
#enddef


def load_arrays_for_grid(root_path: str) -> dict[str, np.ndarray]:
    with uproot.open(root_path) as f:
        if TREE_NAME not in f:
            raise KeyError(f"TTree '{TREE_NAME}' not found in '{root_path}'.")
        #endif
        t = f[TREE_NAME]

        needed = [
            "x", "Q2",
            "e_p", "e_theta", "e_phi",
            "p_p", "p_theta", "p_phi",
            "Mx2",
            "phi",
        ]
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


def fill_grid_variable(arrs: dict[str, np.ndarray], var_name: str) -> list[list[np.ndarray]]:
    """
    Build 4x6 buckets in (xB, -tprime) using recomputed tprime.
    Store the requested variable's values in each pad bucket.

    var_name options:
      - "Mx2"   : store Mx2 branch
      - "Q2"    : store Q2 branch
      - "xB"    : store x branch
      - "tneg"  : store computed (-tprime)
      - "phi"   : store phi branch, wrapped to [0, 2*pi) (radians)
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
    phi = arrs["phi"]

    n = len(xB)

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

        v = float("nan")
        if var_name == "Mx2":
            v = float(Mx2[i])
        elif var_name == "Q2":
            v = float(Q2[i])
        elif var_name == "xB":
            v = float(xB[i])
        elif var_name == "tneg":
            v = float(tneg)
        elif var_name == "phi":
            v = wrap_phi_rad(float(phi[i]))
        else:
            raise RuntimeError("Unknown var_name: " + str(var_name))
        #endif

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
                   out_png: str,
                   xmin: float,
                   xmax: float,
                   nbins: int,
                   x_title: str,
                   y_title: str,
                   fig_title: str) -> None:
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    nrows = len(XB_EDGES) - 1
    ncols = len(TNEG_EDGES) - 1

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.8 * ncols, 2.8 * nrows),
        sharex=True,
        sharey=False
    )

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

            x_mc, y_mc = hist_density(mc_arr, nbins, xmin, xmax)
            x_da, y_da = hist_density(data_arr, nbins, xmin, xmax)

            ax.step(x_mc, y_mc, where="mid", linewidth=1.2, label="MC")
            ax.step(x_da, y_da, where="mid", linewidth=1.2, label="data")

            ax.grid(True)
            ax.set_xlim(xmin, xmax)

            ax.set_title(
                TITLE_BIN_FMT.format(xb_lo=xb_lo, xb_hi=xb_hi, t_lo=t_lo, t_hi=t_hi),
                fontsize=9
            )

            if r == nrows - 1:
                ax.set_xlabel(x_title)
            #endif
            if c == 0:
                ax.set_ylabel(y_title)
            #endif

            if r == 0 and c == 0:
                ax.legend(loc="upper right", frameon=True, fontsize=10)
            #endif
        #endfor
    #endfor

    fig.suptitle(fig_title, fontsize=14)

    plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.subplots_adjust(left=0.06)

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

    try:
        # Preserve original overlays (unchanged)
        plot_overlay(
            out_png=OUT_OVERLAY_1,
            files=[root1, root2],
            labels=["At Rest", "Simulated Fermi Motion"],
            title=r"$M_{x}^{2}$ comparison",
            nbins=MX2_NBINS_OVERLAY,
            xmin=MX2_MIN,
            xmax=MX2_MAX,
        )

        plot_overlay(
            out_png=OUT_OVERLAY_2,
            files=[root3, root4, root5],
            labels=["RGA", "RGA simulated Fermi Motion", "RGC"],
            title=r"$M_{x}^{2}$ data comparison",
            nbins=MX2_NBINS_OVERLAY,
            xmin=MX2_MIN,
            xmax=MX2_MAX,
        )

        # Load arrays once for the grid comparisons (root2 vs root5)
        require_file(root2)
        require_file(root5)
        mc_arrs = load_arrays_for_grid(root2)
        data_arrs = load_arrays_for_grid(root5)

        # Existing Mx2 grid (preserved)
        mc_grid_mx2 = fill_grid_variable(mc_arrs, "Mx2")
        data_grid_mx2 = fill_grid_variable(data_arrs, "Mx2")
        make_grid_plot(
            mc_grid=mc_grid_mx2,
            data_grid=data_grid_mx2,
            out_png=OUT_GRID_MX2,
            xmin=MX2_MIN,
            xmax=MX2_MAX,
            nbins=MX2_NBINS_GRID,
            x_title=LABEL_MX2,
            y_title=YLABEL_DENSITY_GEV2,
            fig_title=r"$M_{x}^{2}$ in $(x_{B},-t^\prime)$ bins: root2 (MC) vs root5 (data)"
        )

        # NEW grids: Q2, xB, -tprime, phi (radians)
        mc_grid_q2 = fill_grid_variable(mc_arrs, "Q2")
        data_grid_q2 = fill_grid_variable(data_arrs, "Q2")
        make_grid_plot(
            mc_grid=mc_grid_q2,
            data_grid=data_grid_q2,
            out_png=OUT_GRID_Q2,
            xmin=Q2_MIN,
            xmax=Q2_MAX,
            nbins=Q2_NBINS_GRID,
            x_title=LABEL_Q2,
            y_title=YLABEL_DENSITY_GEV2,
            fig_title=r"$Q^{2}$ in $(x_{B},-t^\prime)$ bins: root2 (MC) vs root5 (data)"
        )

        mc_grid_xb = fill_grid_variable(mc_arrs, "xB")
        data_grid_xb = fill_grid_variable(data_arrs, "xB")
        make_grid_plot(
            mc_grid=mc_grid_xb,
            data_grid=data_grid_xb,
            out_png=OUT_GRID_XB,
            xmin=XB_MIN,
            xmax=XB_MAX,
            nbins=XB_NBINS_GRID,
            x_title=LABEL_XB,
            y_title=YLABEL_DENSITY_DIMLESS,
            fig_title=r"$x_{B}$ in $(x_{B},-t^\prime)$ bins: root2 (MC) vs root5 (data)"
        )

        mc_grid_tneg = fill_grid_variable(mc_arrs, "tneg")
        data_grid_tneg = fill_grid_variable(data_arrs, "tneg")
        make_grid_plot(
            mc_grid=mc_grid_tneg,
            data_grid=data_grid_tneg,
            out_png=OUT_GRID_TNEG,
            xmin=TNEG_MIN,
            xmax=TNEG_MAX,
            nbins=TNEG_NBINS_GRID,
            x_title=LABEL_TPRIME_NEG,
            y_title=YLABEL_DENSITY_GEV2,
            fig_title=r"$-t^\prime$ in $(x_{B},-t^\prime)$ bins: root2 (MC) vs root5 (data)"
        )

        mc_grid_phi = fill_grid_variable(mc_arrs, "phi")
        data_grid_phi = fill_grid_variable(data_arrs, "phi")
        make_grid_plot(
            mc_grid=mc_grid_phi,
            data_grid=data_grid_phi,
            out_png=OUT_GRID_PHI,
            xmin=PHI_MIN,
            xmax=PHI_MAX,
            nbins=PHI_NBINS_GRID,
            x_title=LABEL_PHI,
            y_title=YLABEL_DENSITY_RAD,
            fig_title=r"$\phi$ in $(x_{B},-t^\prime)$ bins: root2 (MC) vs root5 (data)"
        )

    except Exception as e:
        return fatal(str(e), 3)
    #endtry

    print("Wrote: " + OUT_OVERLAY_1)
    print("Wrote: " + OUT_OVERLAY_2)
    print("Wrote: " + OUT_GRID_MX2)
    print("Wrote: " + OUT_GRID_Q2)
    print("Wrote: " + OUT_GRID_XB)
    print("Wrote: " + OUT_GRID_TNEG)
    print("Wrote: " + OUT_GRID_PHI)
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif