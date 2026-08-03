#!/usr/bin/env python3
"""
plot_first_100k_dvcsgen_particle_kinematics.py

Read the first 100,000 entries from the specified ROOT file and plot:

    p1_p
    p1_theta
    p2_p
    p2_theta

The theta branches are assumed to be stored in radians and are converted to
degrees before plotting.

Dependencies:
    pip install uproot numpy matplotlib
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import uproot


DEFAULT_INPUT = Path(
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
    "dvcsgen/gen_dvcsgen_rga_fa18_out_50nA_10604MeV.root"
)
DEFAULT_OUTPUT = Path(
    "first_100k_dvcsgen_particle_kinematics.png"
)
DEFAULT_MAX_ENTRIES = 100_000

BRANCHES: tuple[str, ...] = (
    "p1_p",
    "p1_theta",
    "p2_p",
    "p2_theta",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot p1_p, p1_theta, p2_p, and p2_theta for the first "
            "entries of a ROOT tree."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input ROOT file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--tree",
        default=None,
        help=(
            "ROOT tree name. If omitted, the script automatically uses the "
            "first TTree found in the file."
        ),
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=DEFAULT_MAX_ENTRIES,
        help=f"Maximum number of entries to read (default: {DEFAULT_MAX_ENTRIES:,})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output PNG path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def find_first_tree(
    root_file: uproot.ReadOnlyDirectory,
) -> tuple[str, uproot.behaviors.TTree.TTree]:
    for key, class_name in root_file.classnames().items():
        if class_name.startswith("TTree"):
            clean_name = key.split(";")[0]
            return clean_name, root_file[key]
        # endif
    # endfor

    available = ", ".join(root_file.keys())
    raise RuntimeError(
        "No TTree was found in the ROOT file. "
        f"Available keys: {available}"
    )


def require_branches(
    tree: uproot.behaviors.TTree.TTree,
    branches: Sequence[str],
) -> None:
    available = set(tree.keys())
    missing = [branch for branch in branches if branch not in available]

    if missing:
        raise KeyError(
            "Required branches are missing:\n  "
            + "\n  ".join(missing)
            + "\nAvailable branches include:\n  "
            + ", ".join(sorted(available))
        )
    # endif


def finite_values(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    return values[np.isfinite(values)]


def main() -> int:
    args = parse_args()

    if args.max_entries <= 0:
        raise ValueError("--max-entries must be positive.")
    # endif

    if not args.input.exists():
        raise FileNotFoundError(f"Input ROOT file does not exist: {args.input}")
    # endif

    with uproot.open(args.input) as root_file:
        if args.tree is None:
            tree_name, tree = find_first_tree(root_file)
        else:
            if args.tree not in root_file:
                available = ", ".join(root_file.keys())
                raise KeyError(
                    f"Tree '{args.tree}' was not found. "
                    f"Available keys: {available}"
                )
            # endif
            tree_name = args.tree
            tree = root_file[args.tree]
        # endif

        require_branches(tree, BRANCHES)

        number_to_read = min(args.max_entries, tree.num_entries)
        arrays = tree.arrays(
            list(BRANCHES),
            entry_start=0,
            entry_stop=number_to_read,
            library="np",
        )
    # endwith

    p1_p = finite_values(arrays["p1_p"])
    p1_theta_deg = np.degrees(finite_values(arrays["p1_theta"]))
    p2_p = finite_values(arrays["p2_p"])
    p2_theta_deg = np.degrees(finite_values(arrays["p2_theta"]))

    plot_definitions = (
        (
            p1_p,
            r"$p_1$ momentum (GeV)",
            (0.0, 5.0),
            120,
        ),
        (
            p1_theta_deg,
            r"$\theta_1$ (deg)",
            (0.0, 70.0),
            120,
        ),
        (
            p2_p,
            r"$p_2$ momentum (GeV)",
            (0.0, 10.0),
            120,
        ),
        (
            p2_theta_deg,
            r"$\theta_2$ (deg)",
            (0.0, 40.0),
            120,
        ),
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    axes_flat = axes.ravel()

    for axis, (values, label, plot_range, bins) in zip(
        axes_flat,
        plot_definitions,
    ):
        axis.hist(
            values,
            bins=bins,
            range=plot_range,
            histtype="step",
            linewidth=1.5,
        )
        axis.set_xlabel(label)
        axis.set_ylabel("events / bin")
        axis.set_xlim(*plot_range)
        axis.set_ylim(bottom=0.0)
        axis.grid(axis="y", alpha=0.25)
        axis.set_title(
            f"{label}: {values.size:,} finite entries"
        )
    # endfor

    fig.suptitle(
        "First "
        f"{number_to_read:,} entries from {args.input.name}\n"
        f"Tree: {tree_name}",
        fontsize=15,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        args.output,
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"Read {number_to_read:,} entries from tree '{tree_name}'.")
    print(f"Wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
