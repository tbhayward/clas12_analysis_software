#!/usr/bin/env python3
"""
Plot Mx2 comparison for baseline vs ISR+FSR and show migrated events.

Usage:
  python plot_mx2_isrfsr.py BASELINE.root ISR_FSR.root

This script:
  - Opens TTree "PhysicsEvents" from both ROOT files.
  - Plots Mx2 (0 → 1.3) for baseline and ISR+FSR.
  - Adds a thin dashed gray histogram of *baseline Mx2* for events whose
    Mx2 was outside (0.81, 1.00) in baseline but inside that window in ISR+FSR.
  - Saves to: output/enpi+/ISR_FSR_Mx2.pdf
"""

import argparse
import os
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt

# Use Uproot to read ROOT TTrees
try:
    import uproot
except ImportError as e:
    raise SystemExit(
        "ERROR: This script requires the 'uproot' package.\n"
        "Install with:  pip install uproot awkward\n"
    ) from e


def load_mx2_evnnum(root_path: str, tree_name: str = "PhysicsEvents") -> Tuple[np.ndarray, np.ndarray]:
    """Load Mx2 and evnum arrays from a ROOT file's TTree."""
    with uproot.open(root_path) as f:
        if tree_name not in f:
            raise KeyError(f'Tree "{tree_name}" not found in file: {root_path}')
        tree = f[tree_name]
        # Read as numpy arrays
        arrs = tree.arrays(["Mx2", "evnum"], library="np")
        if "Mx2" not in arrs or "evnum" not in arrs:
            raise KeyError('Branches "Mx2" and/or "evnum" not found in tree.')
        mx2 = np.asarray(arrs["Mx2"]).ravel()
        evn = np.asarray(arrs["evnum"]).ravel()
        # Some branches might be floats; cast evnum to integers for matching
        # but keep original ordering/length for one-to-one matching by index.
        # We'll match using np.isin, so casting is fine.
        # If evnum is very large, use int64.
        if not np.issubdtype(evn.dtype, np.integer):
            evn = evn.astype(np.int64, copy=False)
        return mx2, evn


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Mx2 for baseline vs ISR+FSR with migrated-event overlay.")
    parser.add_argument("baseline", help="Path to baseline ROOT file")
    parser.add_argument("isrfsr", help="Path to ISR+FSR ROOT file")
    parser.add_argument(
        "-o",
        "--output",
        default="output/enpi+/ISR_FSR_Mx2.pdf",
        help='Output PDF path (default: "output/enpi+/ISR_FSR_Mx2.pdf")',
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=260,
        help="Number of bins between 0 and 2.6 (default: 130; i.e. 0.01 width)",
    )
    args = parser.parse_args()

    # Load data
    mx2_base, evn_base = load_mx2_evnnum(args.baseline, "PhysicsEvents")
    mx2_isr, evn_isr = load_mx2_evnnum(args.isrfsr, "PhysicsEvents")

    # Selection window (strict inequalities, as specified: 0.81 < Mx2 < 1.00)
    low, high = 0.81, 1.00
    base_outside_mask = (mx2_base <= low) | (mx2_base >= high)
    isr_inside_mask = (mx2_isr > low) & (mx2_isr < high)

    # Match by evnnum: find baseline events *outside* that have same evnnum as ISR+FSR events *inside*
    evn_isr_inside = evn_isr[isr_inside_mask]
    evn_base_outside = evn_base[base_outside_mask]
    # For duplicates, np.isin preserves multiplicity on the baseline side (good for histogram counts)
    migrated_on_baseline_mask = base_outside_mask.copy()
    migrated_on_baseline_mask[base_outside_mask] = np.isin(evn_base_outside, evn_isr_inside)

    # Data to plot
    data_baseline_all = mx2_base
    data_isrfsr_all = mx2_isr
    data_migrated_baseline = mx2_base[migrated_on_baseline_mask]

    # Histogram config
    x_min, x_max = 0.0, 2.6
    bins = np.linspace(x_min, x_max, args.bins + 1)

    # Plot
    plt.figure(figsize=(7.5, 5.5))
    # Baseline (solid)
    plt.hist(
        data_baseline_all,
        bins=bins,
        histtype="step",
        linewidth=1.6,
        label="baseline",
    )
    # ISR+FSR (solid)
    plt.hist(
        data_isrfsr_all,
        bins=bins,
        histtype="step",
        linewidth=1.6,
        label="sim. ISR+FSR",
    )
    # Migrated (baseline values) – thin dashed gray, no legend entry per instructions
    if data_migrated_baseline.size > 0:
        plt.hist(
            data_migrated_baseline,
            bins=bins,
            histtype="step",
            linewidth=1.0,
            linestyle="--",
            color="gray",
        )

    # Labels & legend
    plt.xlabel(r"$M_{x}^{2}$ (GeV$^{2}$)")
    plt.ylabel("counts")
    plt.xlim(x_min, x_max)
    plt.legend(loc="upper right", frameon=False)

    # Nicely mark the selection window as a visual cue (optional; remove if undesired)
    # plt.axvspan(low, high, color="k", alpha=0.06, lw=0)

    # Save
    out_path = args.output
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"Saved plot to: {out_path}")


if __name__ == "__main__":
    main()