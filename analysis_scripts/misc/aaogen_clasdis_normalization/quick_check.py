#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
quick_mx2_hist.py

Reads the TTree "PhysicsEvents" from:
  /volatile/clas12/thayward/montecarlo/fermi_motion/clasdis_rgc_sp23_inb.root

Makes a 1D histogram of the branch "Mx2" from 0.3 to 1.6 and saves:
  output/quick_test.png
"""

import os
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt


def main() -> int:
    infile = "/volatile/clas12/thayward/montecarlo/fermi_motion/clasdis_rgc_sp23_inb.root"
    tree_name = "PhysicsEvents"
    branch = "Mx2"

    xmin = 0.3
    xmax = 1.6
    nbins = 260

    out_png = "output/quick_test.png"
    os.makedirs(os.path.dirname(out_png), exist_ok=True)

    if not os.path.isfile(infile):
        print(f"FATAL: input ROOT file not found: {infile}", file=sys.stderr)
        return 1
    #endif

    try:
        with uproot.open(infile) as f:
            if tree_name not in f:
                keys = [k for k in f.keys()]
                print(f"FATAL: TTree '{tree_name}' not found in file. Keys: {keys}", file=sys.stderr)
                return 2
            #endif

            t = f[tree_name]

            if branch not in t.keys():
                brs = [k for k in t.keys()]
                print(f"FATAL: branch '{branch}' not found in TTree '{tree_name}'. Available branches include: {brs[:50]}{' ...' if len(brs) > 50 else ''}", file=sys.stderr)
                return 3
            #endif

            arr = t[branch].array(library="np")
    except Exception as e:
        print(f"FATAL: failed to read ROOT file/tree/branch: {e}", file=sys.stderr)
        return 4
    #endtry

    arr = np.asarray(arr, dtype=np.float64)
    arr = arr[np.isfinite(arr)]

    in_range = (arr >= xmin) & (arr <= xmax)
    arr_plot = arr[in_range]

    fig = plt.figure(figsize=(9, 6))
    plt.hist(arr_plot, bins=nbins, range=(xmin, xmax), histtype="step", linewidth=1.5)

    plt.xlabel("Mx2")
    plt.ylabel("Counts")
    plt.title(f"Mx2 from {tree_name} ({os.path.basename(infile)})")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"Wrote: {out_png}")
    print(f"Entries (finite): {arr.size}")
    print(f"Entries in range [{xmin}, {xmax}]: {arr_plot.size}")

    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif