#!/usr/bin/env python3
"""
quick_check_photon_truth_samples.py

Quick validation of the new Fa18 Inb photon-truth ROOT samples.

This is intentionally diagnostic and plot-focused.  It checks that:
  * all expected trees/branches are readable;
  * reconstructed kinematics look sensible across DVCSgen, AAOgen, CLASDIS,
    and data;
  * REC::Particle -> MC::RecMatch -> MC::Particle matching behaves sensibly;
  * MC::Lund parent/grandparent PID information is populated where expected;
  * CLASDIS photons can be separated into pi0 / eta / other / fake/unmatched
    truth categories;
  * epgammagamma samples show a pi0 mass peak when both photons share the same
    pi0 parent;
  * reconstructed/generated kinematics close where gen_valid == 1;
  * FT/FD detector populations look sensible.

Outputs are written under:
    output/photon_truth_quick_check/

No CSV or NPZ files are produced.

Example:
    python quick_check_photon_truth_samples.py

Optional:
    python quick_check_photon_truth_samples.py --output output/my_check
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import uproot


# -----------------------------------------------------------------------------
# Hard-coded Fa18 Inb 100k samples
# -----------------------------------------------------------------------------

BASE = Path("/work/clas12/thayward/pi0_BDT/training_sample/fa18_inb")

FILES: Dict[str, Path] = {
    "dvcsgen_epg": BASE / "dvcsgen_rga_fa18_inb_epgammaX_100k.root",
    "aaogen_epg": BASE / "aaogen_rga_fa18_inb_epgammaX_100k.root",
    "aaogen_epgg": BASE / "aaogen_rga_fa18_inb_epgammagammaX_100k.root",
    "clasdis_epg": BASE / "clasdis_rga_fa18_inb_epgammaX_100k.root",
    "clasdis_epgg": BASE / "clasdis_rga_fa18_inb_epgammagammaX_100k.root",
    "data_epg": BASE / "data_rga_fa18_inb_epgammaX_100k.root",
    "data_epgg": BASE / "data_rga_fa18_inb_epgammagammaX_100k.root",
}

DISPLAY_NAMES = {
    "dvcsgen_epg": "DVCSgen epgamma",
    "aaogen_epg": "AAOgen epgamma",
    "aaogen_epgg": "AAOgen epgammagamma",
    "clasdis_epg": "CLASDIS epgamma",
    "clasdis_epgg": "CLASDIS epgammagamma",
    "data_epg": "Data epgamma",
    "data_epgg": "Data epgammagamma",
}

EPG_KEYS = ["dvcsgen_epg", "aaogen_epg", "clasdis_epg", "data_epg"]
EPGG_KEYS = ["aaogen_epgg", "clasdis_epgg", "data_epgg"]
MC_EPG_KEYS = ["dvcsgen_epg", "aaogen_epg", "clasdis_epg"]
MC_EPGG_KEYS = ["aaogen_epgg", "clasdis_epgg"]


# -----------------------------------------------------------------------------
# Branch requirements
# -----------------------------------------------------------------------------

COMMON_EPG_BRANCHES = [
    "is_mc",
    "runnum",
    "detector2",
    "p2_p",
    "p2_theta",
    "p2_phi",
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    "Emiss2",
    "theta_gamma_gamma",
    "pTmiss",
]

MC_EPG_BRANCHES = [
    "matching_gamma_pid",
    "gamma_mcindex",
    "gamma_parent_index",
    "gamma_parent_pid",
    "gamma_grandparent_index",
    "gamma_grandparent_pid",
    "gen_valid",
]

COMMON_EPGG_BRANCHES = [
    "is_mc",
    "runnum",
    "detector_gamma1",
    "detector_gamma2",
    "p2_p",
    "p2_theta",
    "p2_phi",
    "p3_p",
    "p3_theta",
    "p3_phi",
    "Mh_gammagamma",
]

MC_EPGG_BRANCHES = [
    "matching_gamma1_pid",
    "gamma1_mcindex",
    "gamma1_parent_index",
    "gamma1_parent_pid",
    "gamma1_grandparent_pid",
    "matching_gamma2_pid",
    "gamma2_mcindex",
    "gamma2_parent_index",
    "gamma2_parent_pid",
    "gamma2_grandparent_pid",
    "gen_valid",
]


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Quick validation of new photon-truth ROOT samples."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/photon_truth_quick_check"),
        help="Output directory for plots.",
    )
    parser.add_argument(
        "--tree",
        type=str,
        default=None,
        help="Optional explicit tree name. Default: automatically find first TTree.",
    )
    return parser.parse_args()


def first_tree_name(path: Path) -> str:
    with uproot.open(path) as root_file:
        for key, obj in root_file.items():
            if hasattr(obj, "num_entries") and hasattr(obj, "arrays"):
                return key.split(";")[0]
            #endif
        #endfor
    raise RuntimeError(f"No TTree-like object found in {path}")


def resolve_tree_name(path: Path, requested: Optional[str]) -> str:
    if requested:
        return requested
    #endif
    return first_tree_name(path)


def load_arrays(
    key: str,
    path: Path,
    requested_tree: Optional[str],
    required: Sequence[str],
    optional: Sequence[str] = (),
) -> Dict[str, np.ndarray]:
    tree_name = resolve_tree_name(path, requested_tree)

    with uproot.open(path) as root_file:
        tree = root_file[tree_name]
        available = {name.split(";")[0] for name in tree.keys()}

        missing = [branch for branch in required if branch not in available]
        if missing:
            raise RuntimeError(
                f"{key}: missing required branches in {path}: {missing}"
            )
        #endif

        branches = list(required)
        branches.extend(
            branch for branch in optional
            if branch in available and branch not in branches
        )

        arrays_raw = tree.arrays(branches, library="np")

    arrays = {name: np.asarray(arrays_raw[name]) for name in branches}
    arrays["_entries"] = np.asarray([len(next(iter(arrays.values())))], dtype=np.int64)
    arrays["_tree_name"] = np.asarray([tree_name], dtype=object)
    return arrays


def finite(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def deg_if_needed(theta: np.ndarray) -> np.ndarray:
    values = np.asarray(theta, dtype=float)
    good = values[np.isfinite(values)]

    if len(good) == 0:
        return values
    #endif

    # Existing processing normally stores radians. Be robust if that changes.
    if np.nanpercentile(np.abs(good), 99) <= 3.5:
        return np.degrees(values)
    #endif

    return values


def normalized_hist(
    ax: plt.Axes,
    values: np.ndarray,
    bins: int,
    hist_range: Optional[Tuple[float, float]],
    label: str,
    *,
    linewidth: float = 1.5,
    linestyle: str = "-",
) -> None:
    vals = finite(values)
    if len(vals) == 0:
        return
    #endif

    ax.hist(
        vals,
        bins=bins,
        range=hist_range,
        density=True,
        histtype="step",
        linewidth=linewidth,
        linestyle=linestyle,
        label=label,
    )


def robust_range(
    arrays: Iterable[np.ndarray],
    qlo: float = 0.005,
    qhi: float = 0.995,
) -> Optional[Tuple[float, float]]:
    vals: List[np.ndarray] = []

    for array in arrays:
        clean = finite(array)
        if len(clean) > 0:
            vals.append(clean)
        #endif
    #endfor

    if not vals:
        return None
    #endif

    merged = np.concatenate(vals)
    lo, hi = np.quantile(merged, [qlo, qhi])

    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        return None
    #endif

    return float(lo), float(hi)


def save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] {path}")


def top_pid_values(values: np.ndarray, n: int = 12) -> Tuple[np.ndarray, np.ndarray]:
    vals = np.asarray(values, dtype=np.int64)
    unique, counts = np.unique(vals, return_counts=True)
    order = np.argsort(counts)[::-1][:n]
    return unique[order], counts[order]


def truth_category_epg(sample: Mapping[str, np.ndarray]) -> np.ndarray:
    matching = np.asarray(sample["matching_gamma_pid"], dtype=np.int64)
    parent = np.asarray(sample["gamma_parent_pid"], dtype=np.int64)
    mcindex = np.asarray(sample["gamma_mcindex"], dtype=np.int64)

    categories = np.full(len(matching), "other", dtype=object)

    categories[mcindex < 0] = "unmatched"
    categories[(mcindex >= 0) & (matching != 22)] = "non-photon truth match"
    categories[(matching == 22) & (parent == 111)] = "pi0 -> gamma"
    categories[(matching == 22) & (parent == 221)] = "eta -> gamma"
    categories[
        (matching == 22)
        & (parent != 111)
        & (parent != 221)
    ] = "other true gamma"

    return categories


def truth_pair_category(sample: Mapping[str, np.ndarray]) -> np.ndarray:
    m1 = np.asarray(sample["matching_gamma1_pid"], dtype=np.int64)
    m2 = np.asarray(sample["matching_gamma2_pid"], dtype=np.int64)
    p1 = np.asarray(sample["gamma1_parent_pid"], dtype=np.int64)
    p2 = np.asarray(sample["gamma2_parent_pid"], dtype=np.int64)
    i1 = np.asarray(sample["gamma1_parent_index"], dtype=np.int64)
    i2 = np.asarray(sample["gamma2_parent_index"], dtype=np.int64)

    categories = np.full(len(m1), "other/mixed", dtype=object)

    same_parent = (i1 > 0) & (i2 > 0) & (i1 == i2)
    both_true_gamma = (m1 == 22) & (m2 == 22)

    categories[
        both_true_gamma & same_parent & (p1 == 111) & (p2 == 111)
    ] = "same pi0 parent"

    categories[
        both_true_gamma & same_parent & (p1 == 221) & (p2 == 221)
    ] = "same eta parent"

    categories[
        both_true_gamma
        & (p1 == 111)
        & (p2 == 111)
        & ~same_parent
    ] = "different pi0 parents"

    categories[
        both_true_gamma
        & same_parent
        & (p1 == p2)
        & (p1 != 111)
        & (p1 != 221)
    ] = "same other parent"

    categories[
        (m1 != 22) | (m2 != 22)
    ] = "non-photon/unmatched daughter"

    return categories


# -----------------------------------------------------------------------------
# Terminal summaries
# -----------------------------------------------------------------------------

def print_sample_summary(samples: Mapping[str, Mapping[str, np.ndarray]]) -> None:
    print("\n" + "=" * 92)
    print("SAMPLE SUMMARY")
    print("=" * 92)

    for key, sample in samples.items():
        n = int(sample["_entries"][0])
        is_mc_values = np.asarray(sample["is_mc"], dtype=np.int64)
        mc_fraction = np.mean(is_mc_values == 1) if n > 0 else float("nan")
        print(
            f"{DISPLAY_NAMES[key]:28s}  entries={n:8,d}  "
            f"is_mc fraction={mc_fraction:7.4f}"
        )
    #endfor


def print_epg_truth_summary(
    samples: Mapping[str, Mapping[str, np.ndarray]]
) -> None:
    print("\n" + "=" * 92)
    print("SINGLE-PHOTON MC TRUTH SUMMARY")
    print("=" * 92)

    for key in MC_EPG_KEYS:
        sample = samples[key]
        n = int(sample["_entries"][0])

        matching = np.asarray(sample["matching_gamma_pid"], dtype=np.int64)
        parent = np.asarray(sample["gamma_parent_pid"], dtype=np.int64)
        gen_valid = np.asarray(sample["gen_valid"], dtype=np.int64)

        matched = sample["gamma_mcindex"] >= 0
        true_gamma = matching == 22
        pi0 = true_gamma & (parent == 111)
        eta = true_gamma & (parent == 221)

        print(f"\n{DISPLAY_NAMES[key]}")
        print(f"  truth matched          : {np.sum(matched):8,d} / {n:,}")
        print(f"  matched as PID 22      : {np.sum(true_gamma):8,d} / {n:,}")
        print(f"  parent PID 111 (pi0)   : {np.sum(pi0):8,d} / {n:,}")
        print(f"  parent PID 221 (eta)   : {np.sum(eta):8,d} / {n:,}")
        print(f"  gen_valid == 1         : {np.sum(gen_valid == 1):8,d} / {n:,}")

        pids, counts = top_pid_values(parent, n=8)
        formatted = ", ".join(
            f"{int(pid)}:{int(count)}" for pid, count in zip(pids, counts)
        )
        print(f"  most common parent PIDs: {formatted}")
    #endfor


def print_epgg_truth_summary(
    samples: Mapping[str, Mapping[str, np.ndarray]]
) -> None:
    print("\n" + "=" * 92)
    print("TWO-PHOTON MC TRUTH SUMMARY")
    print("=" * 92)

    for key in MC_EPGG_KEYS:
        sample = samples[key]
        categories = truth_pair_category(sample)
        unique, counts = np.unique(categories, return_counts=True)
        order = np.argsort(counts)[::-1]

        print(f"\n{DISPLAY_NAMES[key]}")
        for category, count in zip(unique[order], counts[order]):
            print(f"  {category:30s}: {int(count):8,d}")
        #endfor
    #endfor


# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

def plot_epg_kinematics(
    samples: Mapping[str, Mapping[str, np.ndarray]],
    output: Path,
) -> None:
    specs = [
        ("p2_p", r"$p_\gamma$ (GeV)", None),
        ("p2_theta", r"$\theta_\gamma$ (deg)", "angle"),
        ("Mx2", r"$M_X^2(ep\gamma)$ (GeV$^2$)", (-1.0, 2.0)),
        ("Mx2_1", r"$M_X^2(ep)$ (GeV$^2$)", None),
        ("Mx2_2", r"$M_X^2(e\gamma)$ (GeV$^2$)", None),
        ("Emiss2", r"$E_{\rm miss}$", None),
        ("theta_gamma_gamma", r"$\theta_{\gamma\gamma}$", None),
        ("pTmiss", r"$p_T^{\rm miss}$ (GeV)", None),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(18.0, 8.8))
    axes = axes.ravel()

    for ax, (branch, xlabel, fixed_range) in zip(axes, specs):
        converted = {}

        for key in EPG_KEYS:
            values = samples[key][branch]
            if branch == "p2_theta":
                values = deg_if_needed(values)
            #endif
            converted[key] = values
        #endfor

        hist_range = fixed_range
        if hist_range is None:
            hist_range = robust_range(converted.values())
        #endif

        for key in EPG_KEYS:
            normalized_hist(
                ax,
                converted[key],
                bins=80,
                hist_range=hist_range,
                label=DISPLAY_NAMES[key],
                linewidth=1.5 if key != "data_epg" else 1.8,
            )
        #endfor

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.2)
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=4,
        fontsize=9,
        frameon=True,
    )
    fig.suptitle("Fa18 Inb: reconstructed epgammaX kinematic sanity check", y=0.995)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    save(fig, output / "01_epgamma_reconstructed_kinematics.png")


def plot_epg_truth_overview(
    samples: Mapping[str, Mapping[str, np.ndarray]],
    output: Path,
) -> None:
    fig, axes = plt.subplots(3, 4, figsize=(18.0, 12.0))

    for row, key in enumerate(MC_EPG_KEYS):
        sample = samples[key]

        fields = [
            ("matching_gamma_pid", "Matched truth PID"),
            ("gamma_parent_pid", "Immediate parent PID"),
            ("gamma_grandparent_pid", "Grandparent PID"),
        ]

        for col, (branch, title) in enumerate(fields):
            ax = axes[row, col]
            pids, counts = top_pid_values(sample[branch], n=12)
            ax.bar(np.arange(len(pids)), counts)
            ax.set_xticks(np.arange(len(pids)))
            ax.set_xticklabels([str(int(pid)) for pid in pids], rotation=45)
            ax.set_title(title)
            ax.set_ylabel("Candidates")
            ax.grid(axis="y", alpha=0.2)
        #endfor

        ax = axes[row, 3]
        gen = np.asarray(sample["gen_valid"], dtype=np.int64)
        labels = ["0", "1"]
        counts = [np.sum(gen == 0), np.sum(gen == 1)]
        ax.bar(labels, counts)
        ax.set_title("gen_valid")
        ax.set_ylabel("Candidates")
        ax.grid(axis="y", alpha=0.2)

        axes[row, 0].text(
            -0.33,
            0.5,
            DISPLAY_NAMES[key],
            transform=axes[row, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=11,
        )
    #endfor

    fig.suptitle(
        "Single-photon MC truth plumbing: PID, ancestry, and generated-event validity",
        y=0.995,
    )
    fig.tight_layout(rect=[0.03, 0.0, 1.0, 0.97])
    save(fig, output / "02_epgamma_truth_pid_overview.png")


def plot_clasdis_truth_categories(
    samples: Mapping[str, Mapping[str, np.ndarray]],
    output: Path,
) -> None:
    sample = samples["clasdis_epg"]
    categories = truth_category_epg(sample)

    order = [
        "pi0 -> gamma",
        "eta -> gamma",
        "other true gamma",
        "non-photon truth match",
        "unmatched",
        "other",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.5))

    ax = axes[0, 0]
    counts = [np.sum(categories == category) for category in order]
    ax.bar(np.arange(len(order)), counts)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, rotation=30, ha="right")
    ax.set_ylabel("Candidates")
    ax.set_title("CLASDIS epgamma truth categories")
    ax.grid(axis="y", alpha=0.2)

    specs = [
        ("p2_p", r"$p_\gamma$ (GeV)", None),
        ("p2_theta", r"$\theta_\gamma$ (deg)", "angle"),
        ("Mx2", r"$M_X^2(ep\gamma)$ (GeV$^2$)", (-1.0, 2.0)),
    ]

    for ax, (branch, xlabel, fixed_range) in zip(axes.ravel()[1:], specs):
        values = sample[branch]
        if branch == "p2_theta":
            values = deg_if_needed(values)
        #endif

        hist_range = fixed_range
        if hist_range is None:
            hist_range = robust_range([values])
        #endif

        for category in order[:5]:
            mask = categories == category
            if np.sum(mask) == 0:
                continue
            #endif

            normalized_hist(
                ax,
                values[mask],
                bins=70,
                hist_range=hist_range,
                label=f"{category} (N={np.sum(mask):,})",
                linewidth=1.5,
            )
        #endfor

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.2)
    #endfor

    axes[0, 1].legend(fontsize=7.5)
    fig.suptitle("CLASDIS single-photon truth composition", y=0.995)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.965])
    save(fig, output / "03_clasdis_epgamma_truth_categories.png")


def plot_detector_populations(
    samples: Mapping[str, Mapping[str, np.ndarray]],
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))

    ax = axes[0]
    x = np.arange(len(EPG_KEYS))
    ft = []
    fd = []

    for key in EPG_KEYS:
        det = np.asarray(samples[key]["detector2"], dtype=np.int64)
        ft.append(np.sum(det == 0))
        fd.append(np.sum(det == 1))
    #endfor

    width = 0.38
    ax.bar(x - width / 2, ft, width=width, label="FT (detector2=0)")
    ax.bar(x + width / 2, fd, width=width, label="FD (detector2=1)")
    ax.set_xticks(x)
    ax.set_xticklabels([DISPLAY_NAMES[key] for key in EPG_KEYS], rotation=25, ha="right")
    ax.set_ylabel("Candidates")
    ax.set_title("epgamma detector population")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.2)

    for panel, gamma_branch in enumerate(["detector_gamma1", "detector_gamma2"], start=1):
        ax = axes[panel]
        x = np.arange(len(EPGG_KEYS))
        ft = []
        fd = []

        for key in EPGG_KEYS:
            det = np.asarray(samples[key][gamma_branch], dtype=np.int64)
            ft.append(np.sum(det == 0))
            fd.append(np.sum(det == 1))
        #endfor

        ax.bar(x - width / 2, ft, width=width, label="FT")
        ax.bar(x + width / 2, fd, width=width, label="FD")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [DISPLAY_NAMES[key] for key in EPGG_KEYS],
            rotation=25,
            ha="right",
        )
        ax.set_ylabel("Candidates")
        ax.set_title(gamma_branch)
        ax.grid(axis="y", alpha=0.2)
    #endfor

    fig.suptitle("Photon detector assignments", y=0.995)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save(fig, output / "04_detector_populations.png")


def plot_mgg_overview(
    samples: Mapping[str, Mapping[str, np.ndarray]],
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))

    ax = axes[0]
    for key in EPGG_KEYS:
        normalized_hist(
            ax,
            samples[key]["Mh_gammagamma"],
            bins=100,
            hist_range=(0.0, 0.8),
            label=DISPLAY_NAMES[key],
            linewidth=1.7 if key == "data_epgg" else 1.5,
        )
    #endfor

    ax.axvline(0.13498, linestyle="--", linewidth=1.0)
    ax.axvline(0.54786, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    ax.set_ylabel("Normalized density")
    ax.set_title("All reconstructed photon pairs")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)

    ax = axes[1]
    sample = samples["clasdis_epgg"]
    categories = truth_pair_category(sample)
    mgg = sample["Mh_gammagamma"]

    pair_order = [
        "same pi0 parent",
        "same eta parent",
        "different pi0 parents",
        "same other parent",
        "other/mixed",
        "non-photon/unmatched daughter",
    ]

    for category in pair_order:
        mask = categories == category
        if np.sum(mask) == 0:
            continue
        #endif

        normalized_hist(
            ax,
            mgg[mask],
            bins=100,
            hist_range=(0.0, 0.8),
            label=f"{category} (N={np.sum(mask):,})",
            linewidth=1.5,
        )
    #endfor

    ax.axvline(0.13498, linestyle="--", linewidth=1.0)
    ax.axvline(0.54786, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    ax.set_ylabel("Normalized density")
    ax.set_title("CLASDIS separated by daughter ancestry")
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.2)

    fig.suptitle("Two-photon mass and truth-parent validation", y=0.995)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save(fig, output / "05_mgammagamma_truth_validation.png")


def plot_parent_pair_matrix(
    samples: Mapping[str, Mapping[str, np.ndarray]],
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5))

    for ax, key in zip(axes, MC_EPGG_KEYS):
        sample = samples[key]
        p1 = np.asarray(sample["gamma1_parent_pid"], dtype=np.int64)
        p2 = np.asarray(sample["gamma2_parent_pid"], dtype=np.int64)

        combined = np.concatenate([p1, p2])
        top_pids, _ = top_pid_values(combined, n=10)
        pid_to_bin = {int(pid): i for i, pid in enumerate(top_pids)}

        matrix = np.zeros((len(top_pids), len(top_pids)), dtype=int)

        for pid1, pid2 in zip(p1, p2):
            if int(pid1) in pid_to_bin and int(pid2) in pid_to_bin:
                matrix[pid_to_bin[int(pid1)], pid_to_bin[int(pid2)]] += 1
            #endif
        #endfor

        image = ax.imshow(matrix, origin="lower", aspect="auto")
        ax.set_xticks(np.arange(len(top_pids)))
        ax.set_yticks(np.arange(len(top_pids)))
        ax.set_xticklabels([str(int(pid)) for pid in top_pids], rotation=45)
        ax.set_yticklabels([str(int(pid)) for pid in top_pids])
        ax.set_xlabel("gamma2 parent PID")
        ax.set_ylabel("gamma1 parent PID")
        ax.set_title(DISPLAY_NAMES[key])
        fig.colorbar(image, ax=ax, label="Pairs")
    #endfor

    fig.suptitle("Immediate-parent PID pair matrix", y=0.995)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save(fig, output / "06_epgammagamma_parent_pid_matrix.png")


def plot_epg_reco_gen_closure(
    samples: Mapping[str, Mapping[str, np.ndarray]],
    output: Path,
) -> None:
    optional_pairs = [
        ("p2_p", "gen_p2_p", r"$(p_\gamma^{rec}-p_\gamma^{gen})/p_\gamma^{gen}$", "relative"),
        ("p2_theta", "gen_p2_theta", r"$\theta_\gamma^{rec}-\theta_\gamma^{gen}$ (deg)", "angle"),
        ("Mx2", "gen_Mx2", r"$M_X^{2,rec}(ep\gamma)-M_X^{2,gen}(ep\gamma)$", "difference"),
        ("pTmiss", "gen_pTmiss", r"$p_T^{miss,rec}-p_T^{miss,gen}$ (GeV)", "difference"),
    ]

    fig, axes = plt.subplots(3, 4, figsize=(18.0, 11.5))

    for row, key in enumerate(MC_EPG_KEYS):
        sample = samples[key]
        gen_valid = np.asarray(sample["gen_valid"], dtype=np.int64) == 1

        for col, (rec_branch, gen_branch, xlabel, mode) in enumerate(optional_pairs):
            ax = axes[row, col]

            if gen_branch not in sample:
                ax.text(0.5, 0.5, f"{gen_branch}\nnot present", ha="center", va="center")
                ax.axis("off")
                continue
            #endif

            rec = np.asarray(sample[rec_branch], dtype=float)
            gen = np.asarray(sample[gen_branch], dtype=float)
            mask = gen_valid & np.isfinite(rec) & np.isfinite(gen)

            if mode == "relative":
                mask &= np.abs(gen) > 1.0e-9
                values = (rec[mask] - gen[mask]) / gen[mask]
            elif mode == "angle":
                rec_deg = deg_if_needed(rec)
                gen_deg = deg_if_needed(gen)
                values = rec_deg[mask] - gen_deg[mask]
            else:
                values = rec[mask] - gen[mask]
            #endif

            hist_range = robust_range([values], qlo=0.01, qhi=0.99)
            normalized_hist(
                ax,
                values,
                bins=80,
                hist_range=hist_range,
                label=f"N={len(values):,}",
            )
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Normalized density")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
        #endfor

        axes[row, 0].text(
            -0.33,
            0.5,
            DISPLAY_NAMES[key],
            transform=axes[row, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=11,
        )
    #endfor

    fig.suptitle("Single-photon reconstructed/generated closure (gen_valid == 1)", y=0.995)
    fig.tight_layout(rect=[0.03, 0.0, 1.0, 0.965])
    save(fig, output / "07_epgamma_reco_gen_closure.png")


def plot_epgg_reco_gen_closure(
    samples: Mapping[str, Mapping[str, np.ndarray]],
    output: Path,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.5))

    for row, key in enumerate(MC_EPGG_KEYS):
        sample = samples[key]
        valid = np.asarray(sample["gen_valid"], dtype=np.int64) == 1

        specs = [
            ("p2_p", "gen_p2_p", r"$(p_{\gamma1}^{rec}-p_{\gamma1}^{gen})/p_{\gamma1}^{gen}$", "relative"),
            ("p3_p", "gen_p3_p", r"$(p_{\gamma2}^{rec}-p_{\gamma2}^{gen})/p_{\gamma2}^{gen}$", "relative"),
            ("Mh_gammagamma", "gen_Mh_gammagamma", r"$M_{\gamma\gamma}^{rec}-M_{\gamma\gamma}^{gen}$ (GeV)", "difference"),
        ]

        for col, (rec_branch, gen_branch, xlabel, mode) in enumerate(specs):
            ax = axes[row, col]

            if gen_branch not in sample:
                ax.text(0.5, 0.5, f"{gen_branch}\nnot present", ha="center", va="center")
                ax.axis("off")
                continue
            #endif

            rec = np.asarray(sample[rec_branch], dtype=float)
            gen = np.asarray(sample[gen_branch], dtype=float)
            mask = valid & np.isfinite(rec) & np.isfinite(gen)

            if mode == "relative":
                mask &= np.abs(gen) > 1.0e-9
                values = (rec[mask] - gen[mask]) / gen[mask]
            else:
                values = rec[mask] - gen[mask]
            #endif

            hist_range = robust_range([values], qlo=0.01, qhi=0.99)
            normalized_hist(
                ax,
                values,
                bins=80,
                hist_range=hist_range,
                label=f"N={len(values):,}",
            )
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Normalized density")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
        #endfor

        axes[row, 0].text(
            -0.30,
            0.5,
            DISPLAY_NAMES[key],
            transform=axes[row, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=11,
        )
    #endfor

    fig.suptitle("Two-photon reconstructed/generated closure (gen_valid == 1)", y=0.995)
    fig.tight_layout(rect=[0.035, 0.0, 1.0, 0.965])
    save(fig, output / "08_epgammagamma_reco_gen_closure.png")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    print("=" * 92)
    print("PHOTON TRUTH QUICK CHECK")
    print("=" * 92)

    missing_files = [path for path in FILES.values() if not path.exists()]
    if missing_files:
        print("ERROR: missing input files:")
        for path in missing_files:
            print(f"  {path}")
        #endfor
        raise SystemExit(1)
    #endif

    samples: Dict[str, Dict[str, np.ndarray]] = {}

    # Load epgamma trees.
    for key in EPG_KEYS:
        required = list(COMMON_EPG_BRANCHES)

        if key in MC_EPG_KEYS:
            required += MC_EPG_BRANCHES
        #endif

        optional = [
            "gen_p2_p",
            "gen_p2_theta",
            "gen_Mx2",
            "gen_pTmiss",
        ]

        print(f"[load] {DISPLAY_NAMES[key]}: {FILES[key]}")
        samples[key] = load_arrays(
            key=key,
            path=FILES[key],
            requested_tree=args.tree,
            required=required,
            optional=optional,
        )
    #endfor

    # Load epgammagamma trees.
    for key in EPGG_KEYS:
        required = list(COMMON_EPGG_BRANCHES)

        if key in MC_EPGG_KEYS:
            required += MC_EPGG_BRANCHES
        #endif

        optional = [
            "gen_p2_p",
            "gen_p3_p",
            "gen_Mh_gammagamma",
        ]

        print(f"[load] {DISPLAY_NAMES[key]}: {FILES[key]}")
        samples[key] = load_arrays(
            key=key,
            path=FILES[key],
            requested_tree=args.tree,
            required=required,
            optional=optional,
        )
    #endfor

    print_sample_summary(samples)
    print_epg_truth_summary(samples)
    print_epgg_truth_summary(samples)

    plot_epg_kinematics(samples, args.output)
    plot_epg_truth_overview(samples, args.output)
    plot_clasdis_truth_categories(samples, args.output)
    plot_detector_populations(samples, args.output)
    plot_mgg_overview(samples, args.output)
    plot_parent_pair_matrix(samples, args.output)
    plot_epg_reco_gen_closure(samples, args.output)
    plot_epgg_reco_gen_closure(samples, args.output)

    print("\n" + "=" * 92)
    print("DONE")
    print("=" * 92)
    print(f"Plots written to: {args.output}")
    print("Most important first checks:")
    print("  02_epgamma_truth_pid_overview.png")
    print("  03_clasdis_epgamma_truth_categories.png")
    print("  05_mgammagamma_truth_validation.png")
    print("  07_epgamma_reco_gen_closure.png")
    print("  08_epgammagamma_reco_gen_closure.png")


if __name__ == "__main__":
    main()
