#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import uproot
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot W distributions in bins of vz_e, normalized in 2.0 < W < 2.5."
    )
    parser.add_argument("root_file", help="Input ROOT file containing tree PhysicsEvents")
    parser.add_argument(
        "--tree",
        default="PhysicsEvents",
        help="Tree name. Default: PhysicsEvents"
    )
    parser.add_argument(
        "--out",
        default="output/W_external_radiation.png",
        help="Output image path. Default: output/W_external_radiation.png"
    )
    return parser.parse_args()


def safe_ratio(num, den):
    ratio = np.full_like(num, np.nan, dtype=float)
    mask = den > 0.0
    ratio[mask] = num[mask] / den[mask]
    return ratio


def main():
    args = parse_args()

    if not os.path.isfile(args.root_file):
        raise FileNotFoundError(f"Input ROOT file does not exist: {args.root_file}")
    #endif

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    tree_path = f"{args.root_file}:{args.tree}"

    branches = ["W", "vz_e"]
    with uproot.open(tree_path) as tree:
        missing = [b for b in branches if b not in tree.keys()]
        if missing:
            raise RuntimeError(
                "Missing required branches in tree "
                f"{args.tree}: {', '.join(missing)}"
            )
        #endif

        arrays = tree.arrays(branches, library="np")
    #endwith

    W = arrays["W"]
    vz_e = arrays["vz_e"]

    valid = np.isfinite(W) & np.isfinite(vz_e)
    W = W[valid]
    vz_e = vz_e[valid]

    # Same style of vz_e binning as the example image:
    # [-5.8, -5.2], [-5.2, -4.8], ..., [-1.2, -0.8] cm
    vz_edges = np.array([
        -5.8, -5.2, -4.8, -4.2, -3.8,
        -3.2, -2.8, -2.2, -1.8, -1.2, -0.8
    ], dtype=float)

    W_min = 2.0
    W_max = 5.0
    n_W_bins = 120
    W_edges = np.linspace(W_min, W_max, n_W_bins + 1)
    W_centers = 0.5 * (W_edges[:-1] + W_edges[1:])
    W_widths = np.diff(W_edges)

    norm_min = 2.0
    norm_max = 2.5
    norm_mask_bins = (W_centers >= norm_min) & (W_centers < norm_max)

    histograms = []
    normalized_histograms = []
    labels = []

    for i in range(len(vz_edges) - 1):
        vz_low = vz_edges[i]
        vz_high = vz_edges[i + 1]

        cut = (vz_e >= vz_low) & (vz_e < vz_high)
        counts, _ = np.histogram(W[cut], bins=W_edges)

        counts = counts.astype(float)
        norm_integral = np.sum(counts[norm_mask_bins])

        if norm_integral > 0.0:
            counts_norm = counts / norm_integral
        else:
            counts_norm = np.zeros_like(counts)
        #endif

        histograms.append(counts)
        normalized_histograms.append(counts_norm)
        labels.append(f"{vz_low:.1f} < vz_e < {vz_high:.1f} cm")
    #endfor

    normalized_histograms = np.array(normalized_histograms)

    reference_hist = normalized_histograms[0]
    ratios = np.array([
        safe_ratio(normalized_histograms[i], reference_hist)
        for i in range(len(normalized_histograms))
    ])

    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(
        nrows=3,
        ncols=5,
        width_ratios=[1.0, 1.0, 1.0, 1.0, 1.35],
        wspace=0.35,
        hspace=0.45
    )

    # Left side: individual W histograms in vz_e bins.
    for i in range(len(vz_edges) - 1):
        row = i // 4
        col = i % 4
        ax = fig.add_subplot(gs[row, col])

        ax.step(W_centers, normalized_histograms[i], where="mid", linewidth=1.5)
        ax.set_xlim(W_min, W_max)

        ymax = np.nanmax(normalized_histograms[i])
        if ymax > 0.0:
            ax.set_ylim(0.0, 1.20 * ymax)
        #endif

        ax.set_title(labels[i], fontsize=10)
        ax.set_xlabel("W (GeV)")
        ax.set_ylabel("Normalized counts")
        ax.tick_params(direction="in")
    #endfor

    # Top-right: all normalized histograms overlaid.
    ax_overlay = fig.add_subplot(gs[0, 4])

    for i in range(len(normalized_histograms)):
        ax_overlay.step(
            W_centers,
            normalized_histograms[i],
            where="mid",
            linewidth=1.3,
            label=labels[i]
        )
    #endfor

    ax_overlay.set_xlim(W_min, W_max)
    ax_overlay.set_title("W distributions normalized in 2.0 < W < 2.5", fontsize=11)
    ax_overlay.set_xlabel("W (GeV)")
    ax_overlay.set_ylabel("Normalized counts")
    ax_overlay.tick_params(direction="in")
    ax_overlay.legend(fontsize=7, frameon=False)

    # Bottom-right: ratios relative to the leftmost vz_e bin.
    ax_ratio = fig.add_subplot(gs[1:, 4])

    for i in range(len(ratios)):
        ax_ratio.step(
            W_centers,
            ratios[i],
            where="mid",
            linewidth=1.3,
            label=labels[i]
        )
    #endfor

    ax_ratio.axhline(1.0, linestyle="--", linewidth=1.0)
    ax_ratio.set_xlim(W_min, W_max)
    ax_ratio.set_ylim(0.0, 2.5)
    ax_ratio.set_title("Ratio to leftmost vz_e bin", fontsize=11)
    ax_ratio.set_xlabel("W (GeV)")
    ax_ratio.set_ylabel("Ratio")
    ax_ratio.tick_params(direction="in")
    ax_ratio.legend(fontsize=7, frameon=False)

    # Leave the unused lower-left slots empty if the number of vz_e bins is less than 12.
    for j in range(len(vz_edges) - 1, 12):
        row = j // 4
        col = j % 4
        ax_empty = fig.add_subplot(gs[row, col])
        ax_empty.axis("off")
    #endfor

    fig.suptitle("External-radiation check: W vs vz_e", fontsize=16)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
#endif