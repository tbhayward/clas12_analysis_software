#!/usr/bin/env python3

import argparse
from pathlib import Path

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
    """Load one branch as a NumPy array."""
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

    return np.asarray(values, dtype=float)


def load_parent_theta(filename):
    """
    Combine p1, p2, and p3 three-momenta and return the polar angle of the
    parent momentum in degrees.

    Assumes p*_theta and p*_phi are stored in radians.
    """
    with uproot.open(filename) as root_file:
        tree = get_tree(root_file)

        required = [
            "p1_p", "p1_theta", "p1_phi",
            "p2_p", "p2_theta", "p2_phi",
            "p3_p", "p3_theta", "p3_phi",
        ]

        missing = [branch for branch in required if branch not in tree.keys()]
        if missing:
            available = ", ".join(tree.keys())
            raise KeyError(
                f"Missing branches in {filename}: {', '.join(missing)}\n"
                f"Available branches:\n{available}"
            )
        #endif

        arrays = tree.arrays(required, library="np")
    #endwith

    px_parent = np.zeros_like(arrays["p1_p"], dtype=float)
    py_parent = np.zeros_like(arrays["p1_p"], dtype=float)
    pz_parent = np.zeros_like(arrays["p1_p"], dtype=float)

    for i in (1, 2, 3):
        p = np.asarray(arrays[f"p{i}_p"], dtype=float)
        theta = np.asarray(arrays[f"p{i}_theta"], dtype=float)
        phi = np.asarray(arrays[f"p{i}_phi"], dtype=float)

        px_parent += p * np.sin(theta) * np.cos(phi)
        py_parent += p * np.sin(theta) * np.sin(phi)
        pz_parent += p * np.cos(theta)
    #endfor

    pt_parent = np.sqrt(px_parent**2 + py_parent**2)
    theta_parent = np.degrees(np.arctan2(pt_parent, pz_parent))

    return theta_parent[np.isfinite(theta_parent)]


def finite(values):
    """Return only finite values."""
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


def percentile_range(data_values, mc_values, low=0.5, high=99.5):
    """Common plotting range from data and MC."""
    combined = np.concatenate((data_values, mc_values))
    return (
        np.percentile(combined, low),
        np.percentile(combined, high),
    )


def main():
    parser = argparse.ArgumentParser(
        description="Quick data/MC comparison of Mh, Mx2, and three-pion parent theta."
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
        default="omega_study_fa18_out_data_mc.png",
        help="Output image filename; it will be written under output/.",
    )
    args = parser.parse_args()

    print(f"Data: {args.data}")
    print(f"MC:   {args.mc}")

    data_mh = finite(load_branch(args.data, "Mh"))
    mc_mh   = finite(load_branch(args.mc, "Mh"))

    data_mx2 = finite(load_branch(args.data, "Mx2"))
    mc_mx2   = finite(load_branch(args.mc, "Mx2"))

    data_parent_theta = load_parent_theta(args.data)
    mc_parent_theta   = load_parent_theta(args.mc)

    print(f"Loaded Mh:           data={len(data_mh):,}, MC={len(mc_mh):,}")
    print(f"Loaded Mx2:          data={len(data_mx2):,}, MC={len(mc_mx2):,}")
    print(
        f"Loaded parent theta: data={len(data_parent_theta):,}, "
        f"MC={len(mc_parent_theta):,}"
    )

    mh_range = percentile_range(data_mh, mc_mh)
    mx2_range = percentile_range(data_mx2, mc_mx2)
    parent_theta_range = percentile_range(
        data_parent_theta,
        mc_parent_theta,
        low=0.0,
        high=100.0,
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    comparisons = [
        (
            axes[0],
            data_mh,
            mc_mh,
            "Mh",
            mh_range,
            r"$M_h$ (GeV)",
        ),
        (
            axes[1],
            data_mx2,
            mc_mx2,
            "Mx2",
            mx2_range,
            r"$M_X^2$ (GeV$^2$)",
        ),
        (
            axes[2],
            data_parent_theta,
            mc_parent_theta,
            "Three-pion parent polar angle",
            parent_theta_range,
            r"$\theta_{\pi^+\pi^-\pi^0}$ (deg)",
        ),
    ]

    for ax, data_values, mc_values, title, hist_range, xlabel in comparisons:
        centers_data, hist_data = normalized_hist(
            data_values,
            bins=100,
            hist_range=hist_range,
        )
        centers_mc, hist_mc = normalized_hist(
            mc_values,
            bins=100,
            hist_range=hist_range,
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
        ax.set_title(title)
        ax.legend(frameon=False)
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(r"Fa18 Out: $ep\pi^+\pi^-\pi^0X$ Data vs. MC")
    fig.tight_layout()

    output_dir = Path("output")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / Path(args.output).name
    fig.savefig(output_path, dpi=200, bbox_inches="tight")

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
#endif
