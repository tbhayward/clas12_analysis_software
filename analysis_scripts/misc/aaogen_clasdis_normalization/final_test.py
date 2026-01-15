#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mx2_compare.py

Usage (5 ROOT inputs total):
  python mx2_compare.py <root1> <root2> <root3> <root4> <root5>

Reads TTree "PhysicsEvents" and branch "Mx2" from each ROOT file.

Canvas 1 -> output/Mx2_comparison.png
  - root1: "At Rest"
  - root2: "Simulated Fermi Motion"

Canvas 2 -> output/Mx2_data_comparison.png
  - root3: "RGA"
  - root4: "RGA simulated Fermi Motion"
  - root5: "RGC"

Both canvases:
  - Histogram range: 0.4 to 2.0
  - Each histogram normalized to unit integral over the plotted range
  - X-axis label: "M_{x}^{2} (GeV^{2})"
"""

import os
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt


def fatal(msg: str, code: int = 1) -> int:
    print(f"FATAL: {msg}", file=sys.stderr)
    return code
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
            raise KeyError(f"branch '{branch_name}' not found in TTree '{tree_name}' for file '{root_path}'. Example branches: {brs[:50]}{' ...' if len(brs) > 50 else ''}")
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

    counts, edges = np.histogram(x, bins=nbins, range=(xmin, xmax))
    integral = float(np.sum(counts))

    if integral <= 0.0:
        density = counts.astype(np.float64)
    else:
        density = counts.astype(np.float64) / integral
    #endif

    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, density
#enddef


def plot_overlay(out_png: str, files: list[str], labels: list[str], title: str, nbins: int, xmin: float, xmax: float) -> None:
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    fig = plt.figure(figsize=(9, 6))

    tree_name = "PhysicsEvents"
    branch_name = "Mx2"

    for path, lab in zip(files, labels):
        arr = load_branch_array(path, tree_name, branch_name)
        x, y = hist_density(arr, nbins, xmin, xmax)
        plt.step(x, y, where="mid", linewidth=1.5, label=lab)
    #endfor

    plt.xlabel("M_{x}^{2} (GeV^{2})")
    plt.ylabel("Normalized counts")
    plt.title(title)
    plt.legend(frameon=True)

    plt.tight_layout()
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
    nbins = 320

    try:
        plot_overlay(
            out_png="output/Mx2_comparison.png",
            files=[root1, root2],
            labels=["At Rest", "Simulated Fermi Motion"],
            title="Mx2 comparison",
            nbins=nbins,
            xmin=xmin,
            xmax=xmax,
        )

        plot_overlay(
            out_png="output/Mx2_data_comparison.png",
            files=[root3, root4, root5],
            labels=["RGA", "RGA simulated Fermi Motion", "RGC"],
            title="Mx2 data comparison",
            nbins=nbins,
            xmin=xmin,
            xmax=xmax,
        )
    except Exception as e:
        return fatal(str(e), 3)
    #endtry

    print("Wrote: output/Mx2_comparison.png")
    print("Wrote: output/Mx2_data_comparison.png")
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif