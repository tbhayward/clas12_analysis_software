#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import uproot


DEFAULT_DATA = "/scratch/thayward/rga_fa18_out_epi+pi-pi0X_test.root"
DEFAULT_MC   = "/scratch/thayward/clasdis_fa18_out_epi+pi-pi0X_test.root"


def get_tree(root_file):
    """Return the first TTree-like object found in the ROOT file."""
    for key, obj in root_file.items():
        if hasattr(obj, "arrays") and hasattr(obj, "keys"):
            return obj
        #endif
    #endfor

    raise RuntimeError(f"No TTree found in {root_file.file_path}")


def load_branch(filename, branch_name):
    """Load one branch as a finite NumPy array."""
    with uproot.open(filename) as root_file:
        tree = get_tree(root_file)

        if branch_name not in tree.keys():
            available = ", ".join(tree.keys())
            raise KeyError(
                f"Branch '{branch_name}' not found in {filename}.\n"
                f"Available branches:\n{available}"
            )
        #endif

        values = tree[branch_name].array(library="np")
    #endwith

    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def normalized_hist(values, bins, hist_range):
    """Return bin centers and a unit-area histogram."""
    counts, edges = np.histogram(values, bins=bins, range=hist_range)

    total = counts.sum()
    if total > 0:
        counts = counts.astype(float) / total
    #endif

    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts


def main():
    parser = argparse.ArgumentParser(
        description="Quick data/MC comparison of Mh and Mx2."
    )
    parser.add_argument(
        "data",
        nargs="?",
        default=DEFAULT_DATA,
        help=f"Data ROOT file (default: {DEFAULT_DATA})",
    )
    parser.add_argument(
        "mc",
        nargs="?",
        default=DEFAULT_MC,
        help=f"MC ROOT file (default: {DEFAULT_MC})",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="fa18_out_eppippimpi0X_data_mc_Mh_Mx2.png",
        help="Output image filename",
    )
    args = parser.parse_args()

    print(f"Data: {args.data}")
    print(f"MC:   {args.mc}")

    data_mh = load_branch(args.data, "Mh")
    mc_mh   = load_branch(args.mc, "Mh")

    data_mx2 = load_branch(args.data, "Mx2")
    mc_mx2   = load_branch(args.mc, "Mx2")

    print(f"Loaded Mh:  data={len(data_mh):,}, MC={len(mc_mh):,}")
    print(f"Loaded Mx2: data={len(data_mx2):,}, MC={len(mc_mx2):,}")

    mh_combined = np.concatenate((data_mh, mc_mh))
    mx2_combined = np.concatenate((data_mx2, mc_mx2))

    mh_range = (
        np.percentile(mh_combined, 0.5),
        np.percentile(mh_combined, 99.5),
    )
    mx2_range = (
        np.percentile(mx2_combined, 0.5),
        np.percentile(mx2_combined, 99.5),
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    comparisons = [
        (axes[0], data_mh, mc_mh, "Mh", mh_range, r"$M_h$ (GeV)"),
        (axes[1], data_mx2, mc_mx2, "Mx2", mx2_range, r"$M_X^2$ (GeV$^2$)"),
    ]

    for ax, data_values, mc_values, branch, hist_range, xlabel in comparisons:
        centers_data, hist_data = normalized_hist(
            data_values, bins=100, hist_range=hist_range
        )
        centers_mc, hist_mc = normalized_hist(
            mc_values, bins=100, hist_range=hist_range
        )

        ax.step(
            centers_data,
            hist_data,
            where="mid",
            linewidth=1.5,
            label="Data",
        )
        ax.step(
            centers_mc,
            hist_mc,
            where="mid",
            linewidth=1.5,
            label="MC",
        )

        ax.set_xlim(hist_range)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Normalized counts")
        ax.set_title(branch)
        ax.legend(frameon=False)
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(r"Fa18 Out: $ep\pi^+\pi^-\pi^0X$ Data vs. MC")
    fig.tight_layout()

    plt.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
#endif
