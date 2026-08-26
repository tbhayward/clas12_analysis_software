#!/usr/bin/env python3
"""
quick_check_clasdis_photon_pi0_topology.py

Very small standalone sanity-check script for the CLASDIS Fa18 Inb/Out files.

It does NO BDT training and NO cross-tree event matching.

For each of:
    fa18_inb
    fa18_out

it reads:
    rec_clasdis_rga_<period>_epgammaX.root
    rec_clasdis_rga_<period>_eppi0X.root

and directly inspects:
    e'p'gammaX:
        p2_theta
        p2_p
        detector2

    e'p'pi0X:
        p2_theta          (reconstructed pi0 theta)
        p2_p              (reconstructed pi0 momentum)
        detector_gamma1
        detector_gamma2
        FT-FT / FT-FD / FD-FD topology from detector_gamma1/gamma2

Outputs are plot-focused and saved under:
    output/clasdis_photon_pi0_quick_check/

Angle branches are auto-detected as radians/degrees for plotting.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import uproot


BASE = Path(
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/clasdis"
)

PERIODS: Dict[str, Dict[str, Path]] = {
    "fa18_inb": {
        "epg": BASE / "rec_clasdis_rga_fa18_inb_epgammaX.root",
        "eppi0": BASE / "rec_clasdis_rga_fa18_inb_eppi0X.root",
    },
    "fa18_out": {
        "epg": BASE / "rec_clasdis_rga_fa18_out_epgammaX.root",
        "eppi0": BASE / "rec_clasdis_rga_fa18_out_eppi0X.root",
    },
}

OUTPUT = Path("output/clasdis_photon_pi0_quick_check")


def first_tree_name(filename: Path) -> str:
    with uproot.open(filename) as f:
        for key, obj in f.items():
            if isinstance(obj, uproot.behaviors.TTree.TTree):
                return key.split(";")[0]
            #endif
        #endfor

    raise RuntimeError(f"No TTree found in {filename}")


def load_arrays(filename: Path, branches: list[str]) -> dict[str, np.ndarray]:
    tree_name = first_tree_name(filename)

    print(f"[load] {filename}")
    print(f"[load] tree = {tree_name}")

    with uproot.open(filename) as f:
        tree = f[tree_name]

        missing = [b for b in branches if b not in tree.keys()]
        if missing:
            raise RuntimeError(
                f"{filename} is missing required branches: {missing}"
            )
        #endif

        print(f"[load] entries = {tree.num_entries:,}")
        arrays = tree.arrays(branches, library="np")

    return {key: np.asarray(value) for key, value in arrays.items()}


def infer_angle_unit(values: np.ndarray) -> str:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]

    if len(finite) == 0:
        return "rad"
    #endif

    q995 = float(np.quantile(np.abs(finite), 0.995))

    if q995 <= 2.0 * np.pi + 0.25:
        return "rad"
    #endif

    return "deg"


def theta_deg(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)

    if infer_angle_unit(values) == "rad":
        return np.rad2deg(values)
    #endif

    return values


def detector_counts(values: np.ndarray) -> Tuple[int, int, int]:
    values = np.asarray(values, dtype=int)

    n_ft = int(np.count_nonzero(values == 0))
    n_fd = int(np.count_nonzero(values == 1))
    n_other = int(len(values) - n_ft - n_fd)

    return n_ft, n_fd, n_other


def pi0_topology(
    detector_gamma1: np.ndarray,
    detector_gamma2: np.ndarray,
) -> np.ndarray:
    g1 = np.asarray(detector_gamma1, dtype=int)
    g2 = np.asarray(detector_gamma2, dtype=int)

    labels = np.full(len(g1), "other", dtype=object)

    labels[(g1 == 0) & (g2 == 0)] = "FT-FT"

    labels[
        ((g1 == 0) & (g2 == 1))
        | ((g1 == 1) & (g2 == 0))
    ] = "FT-FD"

    labels[(g1 == 1) & (g2 == 1)] = "FD-FD"

    return labels


def add_count_labels(ax, values):
    for i, value in enumerate(values):
        ax.text(
            i,
            value,
            f"{int(value):,}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    #endfor


def plot_period(
    period: str,
    epg: dict[str, np.ndarray],
    eppi0: dict[str, np.ndarray],
) -> None:
    outdir = OUTPUT / period
    outdir.mkdir(parents=True, exist_ok=True)

    epg_theta = theta_deg(epg["p2_theta"])
    pi0_theta = theta_deg(eppi0["p2_theta"])

    epg_detector = np.asarray(epg["detector2"], dtype=int)
    g1_detector = np.asarray(eppi0["detector_gamma1"], dtype=int)
    g2_detector = np.asarray(eppi0["detector_gamma2"], dtype=int)

    topology = pi0_topology(g1_detector, g2_detector)

    # ---------------------------------------------------------------
    # Canvas 1: raw angular and detector distributions
    # ---------------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    axes[0, 0].hist(
        epg_theta[np.isfinite(epg_theta)],
        bins=100,
        histtype="step",
        linewidth=1.5,
    )
    axes[0, 0].set_xlabel(r"$\theta_\gamma$ (deg)")
    axes[0, 0].set_ylabel("e'p'gammaX candidates")
    axes[0, 0].set_title("e'p'gammaX photon theta")
    axes[0, 0].grid(alpha=0.2)

    ft_mask = epg_detector == 0
    fd_mask = epg_detector == 1

    axes[0, 1].hist(
        epg_theta[ft_mask & np.isfinite(epg_theta)],
        bins=100,
        histtype="step",
        linewidth=1.5,
        label="detector2 = 0 (FT)",
    )
    axes[0, 1].hist(
        epg_theta[fd_mask & np.isfinite(epg_theta)],
        bins=100,
        histtype="step",
        linewidth=1.5,
        label="detector2 = 1 (FD)",
    )
    axes[0, 1].set_xlabel(r"$\theta_\gamma$ (deg)")
    axes[0, 1].set_ylabel("Candidates")
    axes[0, 1].set_title("Photon theta by detector2")
    axes[0, 1].legend()
    axes[0, 1].grid(alpha=0.2)

    epg_counts = detector_counts(epg_detector)
    axes[0, 2].bar(
        ["FT (0)", "FD (1)", "Other"],
        epg_counts,
    )
    add_count_labels(axes[0, 2], epg_counts)
    axes[0, 2].set_ylabel("e'p'gammaX candidates")
    axes[0, 2].set_title("Raw detector2")
    axes[0, 2].grid(axis="y", alpha=0.2)

    axes[1, 0].hist(
        pi0_theta[np.isfinite(pi0_theta)],
        bins=100,
        histtype="step",
        linewidth=1.5,
    )
    axes[1, 0].set_xlabel(r"$\theta_{\pi^0}$ (deg)")
    axes[1, 0].set_ylabel("e'p'pi0X candidates")
    axes[1, 0].set_title("Reconstructed pi0 theta")
    axes[1, 0].grid(alpha=0.2)

    daughter_counts = detector_counts(
        np.concatenate([g1_detector, g2_detector])
    )
    axes[1, 1].bar(
        ["FT (0)", "FD (1)", "Other"],
        daughter_counts,
    )
    add_count_labels(axes[1, 1], daughter_counts)
    axes[1, 1].set_ylabel("Pi0 daughter photons")
    axes[1, 1].set_title("detector_gamma1 + detector_gamma2")
    axes[1, 1].grid(axis="y", alpha=0.2)

    topology_labels = ["FT-FT", "FT-FD", "FD-FD", "other"]
    topology_counts = [
        int(np.count_nonzero(topology == label))
        for label in topology_labels
    ]

    axes[1, 2].bar(
        topology_labels,
        topology_counts,
    )
    add_count_labels(axes[1, 2], topology_counts)
    axes[1, 2].set_ylabel("e'p'pi0X candidates")
    axes[1, 2].set_title("Direct pi0 daughter topology")
    axes[1, 2].grid(axis="y", alpha=0.2)

    fig.suptitle(
        f"{period}: raw CLASDIS photon / pi0 detector sanity check\n"
        f"N(epgammaX)={len(epg_theta):,}, "
        f"N(eppi0X)={len(pi0_theta):,}"
    )
    fig.tight_layout()

    fig.savefig(
        outdir / "raw_theta_and_detector_distributions.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    # ---------------------------------------------------------------
    # Canvas 2: inspect gamma1 and gamma2 detector IDs separately
    # ---------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    g1_counts = detector_counts(g1_detector)
    g2_counts = detector_counts(g2_detector)

    axes[0].bar(["FT (0)", "FD (1)", "Other"], g1_counts)
    add_count_labels(axes[0], g1_counts)
    axes[0].set_ylabel("e'p'pi0X candidates")
    axes[0].set_title("detector_gamma1")
    axes[0].grid(axis="y", alpha=0.2)

    axes[1].bar(["FT (0)", "FD (1)", "Other"], g2_counts)
    add_count_labels(axes[1], g2_counts)
    axes[1].set_ylabel("e'p'pi0X candidates")
    axes[1].set_title("detector_gamma2")
    axes[1].grid(axis="y", alpha=0.2)

    detector_values = sorted(
        set(g1_detector.tolist())
        | set(g2_detector.tolist())
    )

    if len(detector_values) == 0:
        detector_values = [0, 1]
    #endif

    matrix = np.zeros(
        (len(detector_values), len(detector_values)),
        dtype=int,
    )

    value_to_index = {
        value: index
        for index, value in enumerate(detector_values)
    }

    for d1, d2 in zip(g1_detector, g2_detector):
        matrix[
            value_to_index[int(d1)],
            value_to_index[int(d2)],
        ] += 1
    #endfor

    image = axes[2].imshow(matrix, aspect="auto")

    axes[2].set_xticks(
        np.arange(len(detector_values)),
        [str(v) for v in detector_values],
    )
    axes[2].set_yticks(
        np.arange(len(detector_values)),
        [str(v) for v in detector_values],
    )
    axes[2].set_xlabel("detector_gamma2")
    axes[2].set_ylabel("detector_gamma1")
    axes[2].set_title("gamma1 vs gamma2 detector ID")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            axes[2].text(
                j,
                i,
                f"{matrix[i, j]:,}",
                ha="center",
                va="center",
            )
        #endfor
    #endfor

    fig.colorbar(image, ax=axes[2], label="Candidates")

    fig.suptitle(f"{period}: raw reconstructed-pi0 daughter detector IDs")
    fig.tight_layout()

    fig.savefig(
        outdir / "pi0_daughter_detector_ids.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    # ---------------------------------------------------------------
    # Canvas 3: momentum distributions by detector category
    # ---------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    epg_p = np.asarray(epg["p2_p"], dtype=float)
    pi0_p = np.asarray(eppi0["p2_p"], dtype=float)

    axes[0].hist(
        epg_p[ft_mask & np.isfinite(epg_p)],
        bins=100,
        histtype="step",
        linewidth=1.5,
        label="FT",
    )
    axes[0].hist(
        epg_p[fd_mask & np.isfinite(epg_p)],
        bins=100,
        histtype="step",
        linewidth=1.5,
        label="FD",
    )
    axes[0].set_xlabel(r"$p_\gamma$ (GeV)")
    axes[0].set_ylabel("e'p'gammaX candidates")
    axes[0].set_title("Photon momentum by detector")
    axes[0].legend()
    axes[0].grid(alpha=0.2)

    for label in topology_labels:
        mask = topology == label

        if np.count_nonzero(mask) == 0:
            continue
        #endif

        axes[1].hist(
            pi0_p[mask & np.isfinite(pi0_p)],
            bins=100,
            histtype="step",
            linewidth=1.5,
            label=label,
        )
    #endfor

    axes[1].set_xlabel(r"$p_{\pi^0}$ (GeV)")
    axes[1].set_ylabel("e'p'pi0X candidates")
    axes[1].set_title("Pi0 momentum by daughter topology")
    axes[1].legend()
    axes[1].grid(alpha=0.2)

    fig.suptitle(f"{period}: momentum distributions")
    fig.tight_layout()

    fig.savefig(
        outdir / "momentum_distributions.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(
        f"[{period}] epgammaX detector2: "
        f"FT={epg_counts[0]:,}, FD={epg_counts[1]:,}, other={epg_counts[2]:,}"
    )
    print(
        f"[{period}] pi0 daughters: "
        f"FT={daughter_counts[0]:,}, FD={daughter_counts[1]:,}, "
        f"other={daughter_counts[2]:,}"
    )
    print(
        f"[{period}] pi0 topology: "
        + ", ".join(
            f"{label}={count:,}"
            for label, count in zip(topology_labels, topology_counts)
        )
    )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    for period, files in PERIODS.items():
        print("")
        print("=" * 80)
        print(period)
        print("=" * 80)

        epg = load_arrays(
            files["epg"],
            [
                "p2_p",
                "p2_theta",
                "detector2",
            ],
        )

        eppi0 = load_arrays(
            files["eppi0"],
            [
                "p2_p",
                "p2_theta",
                "detector_gamma1",
                "detector_gamma2",
            ],
        )

        plot_period(
            period=period,
            epg=epg,
            eppi0=eppi0,
        )
    #endfor

    print("")
    print(f"Done. Plots written under {OUTPUT}/")


if __name__ == "__main__":
    main()
