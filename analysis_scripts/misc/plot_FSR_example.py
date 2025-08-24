#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
quick_fsr_hist.py

Reads a ROOT file at runtime, loads the "PhysicsEvents" TTree, selects entries
with detector == 2, and makes a histogram of the "open_angle" branch.

Saves to: output/enpi+/FSR_example.pdf

Usage:
    python quick_fsr_hist.py /path/to/file.root
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for batch environments
import matplotlib.pyplot as plt

try:
    import uproot
except ImportError:
    print("[ERROR] The 'uproot' package is required. Try: pip install uproot")
    sys.exit(1)
#endif


def main():
    # -----------------------
    # Parse arguments
    # -----------------------
    parser = argparse.ArgumentParser(description="Plot open_angle for detector==2.")
    parser.add_argument("root_file", help="Path to the input ROOT file")
    args = parser.parse_args()

    root_path = Path(args.root_file)
    if not root_path.is_file():
        print(f"[ERROR] File not found: {root_path}")
        sys.exit(1)
    #endif

    # -----------------------
    # Load tree and branches
    # -----------------------
    try:
        with uproot.open(root_path) as f:
            if "PhysicsEvents" not in f:
                print("[ERROR] TTree 'PhysicsEvents' not found in the file.")
                sys.exit(1)
            #endif
            tree = f["PhysicsEvents"]
            # Load as numpy arrays
            arrs = tree.arrays(["open_angle", "detector"], library="np")
    except Exception as e:
        print(f"[ERROR] Could not open or read the file: {e}")
        sys.exit(1)
    #endif

    if "open_angle" not in arrs or "detector" not in arrs:
        print("[ERROR] Missing required branches 'open_angle' and/or 'detector'.")
        sys.exit(1)
    #endif

    open_angle = arrs["open_angle"]
    detector = arrs["detector"]

    # Ensure flat numpy arrays
    open_angle = np.asarray(open_angle).ravel()
    detector = np.asarray(detector).ravel()

    # Some files store detector as float; cast to int safely
    # We also guard against NaNs.
    det_int = detector.astype(np.int64, copy=False)

    # -----------------------
    # Apply selection: detector == 1 (FD)
    # -----------------------
    sel = (det_int == 1)
    oa_sel = open_angle[sel]

    # Optional: keep only finite values within the plot range [0, 60]
    oa_sel = oa_sel[np.isfinite(oa_sel)]
    oa_sel = oa_sel[(oa_sel >= 0.0) & (oa_sel <= 60.0)]

    # -----------------------
    # Make histogram
    # -----------------------
    # Define bins from 0 to 60 (inclusive) with 60 bins of width 1.0
    bins = np.linspace(0.0, 60.0, 61)

    fig = plt.figure(figsize=(6.4, 4.8), dpi=150)
    ax = fig.add_subplot(111)

    ax.hist(
        oa_sel,
        bins=bins,
        histtype="stepfilled",
        alpha=0.75,
        linewidth=1.0
    )

    ax.set_xlim(0.0, 60.0)
    ax.set_xlabel(r"$\theta_{e'\gamma}$")
    ax.set_ylabel("counts")
    ax.set_title("open angle")

    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

    # -----------------------
    # Save output
    # -----------------------
    out_dir = Path("output/enpi+")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "FSR_example.pdf"
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)

    print(f"[OK] Saved histogram to: {out_path}")
# endif


if __name__ == "__main__":
    main()
# endif