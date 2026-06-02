#!/usr/bin/env python3

import os
import argparse
import numpy as np
import uproot
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot W distributions in bins of vz_e, normalize each distribution "
            "in 2.0 < W < 2.5, plot ratios to the leftmost vz_e bin, and "
            "plot the Mx2 vs W correlation."
        )
    )
    parser.add_argument(
        "root_file",
        help="Input ROOT file containing tree PhysicsEvents"
    )
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
    parser.add_argument(
        "--n-W-bins",
        type=int,
        default=160,
        help="Number of W bins. Default: 160"
    )
    parser.add_argument(
        "--n-Mx2-bins",
        type=int,
        default=160,
        help="Number of Mx2 bins for the correlation plot. Default: 160"
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

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    #endif

    branches = ["W", "vz_e", "Mx2"]

    with uproot.open(args.root_file) as root_file:
        if args.tree not in root_file:
            raise RuntimeError(f"Tree not found in file: {args.tree}")
        #endif

        tree = root_file[args.tree]

        missing = [branch for branch in branches if branch not in tree.keys()]
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
    Mx2 = arrays["Mx2"]

    valid = np.isfinite(W) & np.isfinite(vz_e) & np.isfinite(Mx2)
    W = W[valid]
    vz_e = vz_e[valid]
    Mx2 = Mx2[valid]

    # Same vz_e binning as the example image:
    # [-5.8, -5.2], [-5.2, -4.8], ..., [-1.2, -0.8] cm
    vz_edges = np.array([
        -5.8, -5.2, -4.8, -4.2, -3.8,
        -3.2, -2.8, -2.2, -1.8, -1.2, -0.8
    ], dtype=float)

    W_min = 0.0
    W_max = 4.0
    W_edges = np.linspace(W_min, W_max, args.n_W_bins + 1)
    W_centers = 0.5 * (W_edges[:-1] + W_edges[1:])

    Mx2_min = 0.0
    Mx2_max = 4.0
    Mx2_edges = np.linspace(Mx2_min, Mx2_max, args.n_Mx2_bins + 1)

    # Normalize every vz_e-bin histogram to have unit integral in this W range.
    norm_min = 2.0
    norm_max = 2.5
    norm_mask_bins = (W_centers >= norm_min) & (W_centers < norm_max)

    normalized_histograms = []
    raw_integrals = []
    norm_integrals = []
    labels = []

    for i in range(len(vz_edges) - 1):
        vz_low = vz_edges[i]
        vz_high = vz_edges[i + 1]

        vz_cut = (vz_e >= vz_low) & (vz_e < vz_high)
        counts, _ = np.histogram(W[vz_cut], bins=W_edges)
        counts = counts.astype(float)

        raw_integral = np.sum(counts)
        norm_integral = np.sum(counts[norm_mask_bins])

        if norm_integral > 0.0:
            counts_norm = counts / norm_integral
        else:
            counts_norm = np.zeros_like(counts)
        #endif

        normalized_histograms.append(counts_norm)
        raw_integrals.append(raw_integral)
        norm_integrals.append(norm_integral)
        labels.append(f"{vz_low:.1f} < vz_e < {vz_high:.1f} cm")
    #endfor

    normalized_histograms = np.array(normalized_histograms)

    reference_histogram = normalized_histograms[0]
    ratios = np.array([
        safe_ratio(normalized_histograms[i], reference_histogram)
        for i in range(len(normalized_histograms))
    ])

    fig, (ax_overlay, ax_ratio, ax_corr) = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(12, 15),
        sharex=False,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.15], "hspace": 0.28}
    )

    for i in range(len(normalized_histograms)):
        ax_overlay.step(
            W_centers,
            normalized_histograms[i],
            where="mid",
            linewidth=1.4,
            label=labels[i]
        )
    #endfor

    ax_overlay.set_xlim(W_min, W_max)
    ax_overlay.set_ylabel("Counts normalized in 2.0 < W < 2.5")
    ax_overlay.set_title("W distributions in bins of vz_e")
    ax_overlay.tick_params(direction="in", top=True, right=True)
    ax_overlay.grid(alpha=0.25)
    ax_overlay.legend(fontsize=8, ncol=2, frameon=False)

    for i in range(len(ratios)):
        ax_ratio.step(
            W_centers,
            ratios[i],
            where="mid",
            linewidth=1.4,
            label=labels[i]
        )
    #endfor

    ax_ratio.axhline(1.0, linestyle="--", linewidth=1.0)
    ax_ratio.set_xlim(W_min, W_max)
    ax_ratio.set_ylim(0.0, 2.5)
    ax_ratio.set_xlabel("W (GeV)")
    ax_ratio.set_ylabel("Ratio to leftmost vz_e bin")
    ax_ratio.set_title("W-distribution ratios")
    ax_ratio.tick_params(direction="in", top=True, right=True)
    ax_ratio.grid(alpha=0.25)

    corr_cut = (
        (W >= W_min) & (W < W_max) &
        (Mx2 >= Mx2_min) & (Mx2 < Mx2_max)
    )

    h2, xedges, yedges = np.histogram2d(
        W[corr_cut],
        Mx2[corr_cut],
        bins=[W_edges, Mx2_edges]
    )

    h2_plot = np.where(h2.T > 0.0, h2.T, np.nan)

    mesh = ax_corr.pcolormesh(
        xedges,
        yedges,
        h2_plot,
        shading="auto"
    )

    cbar = fig.colorbar(mesh, ax=ax_corr)
    cbar.set_label("Counts")

    ax_corr.set_xlim(W_min, W_max)
    ax_corr.set_ylim(Mx2_min, Mx2_max)
    ax_corr.set_xlabel("W (GeV)")
    ax_corr.set_ylabel("Mx2 (GeV^2)")
    ax_corr.set_title("Correlation between Mx2 and W")
    ax_corr.tick_params(direction="in", top=True, right=True)
    ax_corr.grid(alpha=0.20)

    fig.suptitle("External-radiation check: W dependence vs vz_e", fontsize=16)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {args.out}")
    print("")
    print("vz_e bin summary:")
    for i in range(len(vz_edges) - 1):
        print(
            f"  {labels[i]:>26s} : "
            f"raw integral = {raw_integrals[i]:.0f}, "
            f"normalization integral 2.0 < W < 2.5 = {norm_integrals[i]:.0f}"
        )
    #endfor


if __name__ == "__main__":
    main()
#endif