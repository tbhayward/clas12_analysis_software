#!/usr/bin/env python3
"""
train_pi0_bdt.py

Truth-labelled BDT study for identifying reconstructed photons that originate
from pi0 decays in e'p'gammaX events.

This version is built around the new truth-enabled processing output.  It no
longer infers CLASDIS pi0 content by matching the e'p'gammaX and e'p'pi0X
trees.  Instead, the MC labels are defined directly from reconstructed-to-MC
matching and, for CLASDIS, the validated MC::Lund ancestry.

Physics question
----------------
For a reconstructed photon candidate in an e'p'gammaX event:

    "Is this reconstructed photon actually a daughter of a pi0?"

Positive class:
    1. AAOgen, matching_gamma_pid == 22
       AAOgen is generated as exclusive ep pi0, so a genuinely truth-matched
       generated photon is a pi0 daughter by generator construction.

    2. CLASDIS, matching_gamma_pid == 22
                  && gamma_parent_pid == 111
       Here the pi0 ancestry is explicit through MC::Lund.

Negative class:
    A. DVCSgen genuine photons:
           matching_gamma_pid == 22

    B. Genuine CLASDIS photons not from pi0:
           matching_gamma_pid == 22
           && gamma_parent_pid != 111
       These are further split diagnostically into eta and other photons.

    C. Reconstructed photon candidates matched to a non-photon truth particle:
           gamma_mcindex >= 0 && matching_gamma_pid != 22
       Kept separately for DVCSgen, AAOgen, and CLASDIS.

    D. Reconstructed photon candidates with no MC::RecMatch association:
           gamma_mcindex < 0
       Also kept separately by source.

The negative categories are NOT allowed to acquire arbitrary importance merely
because one generator produces many more examples.  The training weights are
hierarchical:
    * total positive weight = total negative weight = 0.5;
    * AAOgen pi0 and CLASDIS pi0 share positive weight equally (0.25 each);
    * the nominal negative-class target weights are:
          genuine DVCS/BH gamma                         0.300,
          genuine CLASDIS non-pi0 gamma                 0.050,
          matched non-photon reconstruction artifacts   0.075,
          unmatched reconstruction artifacts            0.075;
    * within any negative macro-category containing multiple fine categories,
      the macro weight is divided in proportion to the observed retained event
      populations rather than equally.  In particular eta and other genuine
      CLASDIS photons follow their actual CLASDIS population ratio.

Separate BDTs are trained for:
    * FT: detector2 == 0
    * FD: detector2 == 1

Common candidate selections:
    * open_angle_ep2 > 7 degrees;
    * optional loose M_X^2(epgamma) window, default -0.5 < Mx2 < 0.5 GeV^2.

The BDT score is a classifier score only.  Because the training priors and
sample weights are deliberately artificial, predict_proba[:,1] must NOT be
interpreted as an absolute physical P(pi0 | event).

Pair-level validation
---------------------
The e'p'gammagammaX trees are never used for training.  After both FT and FD
models are trained, each photon in a reconstructed pair is scored using the
model corresponding to its own detector region.  CLASDIS pairs are then
separated by truth ancestry (same pi0 parent, different pi0 parents, same eta
parent, etc.) and compared against reconstructed M_gamma_gamma.  This provides
an independent closure test of the single-photon classifier.

Primary outputs
---------------
Plots are preferred over text/CSV products.  The script writes:

    output/pi0_bdt_truth/<period>/<sample_tag>/
        truth_category_counts.png

        FT/
            01_training_category_counts.png
            02_input_features_core.png
            03_input_features_reco_artifacts.png
            04_roc_curves.png
            04b_core_operating_points.png
            05_bdt_score_by_truth_category.png
            06_bdt_score_core_and_data.png
            06b_low_mx2ep_score_shift.png
            07_feature_importance.png
            08_confusion_matrix_score_0p5.png
            pi0_bdt_model.joblib

        FD/
            same set of outputs

        09_clasdis_pair_score_validation.png
        10_clasdis_mgg_vs_pair_score.png
        11_data_mgg_by_pair_score.png
        12_aaogen_mgg_by_pair_score.png
        13_data_pi0_mass_fits_vs_score.png
        14_data_pi0_peak_fraction_vs_score.png

No CSV files are produced.

Examples
--------
Current in-progress temporary samples:

    python train_pi0_bdt.py --period fa18_inb

When the uncapped/full processing products exist:

    python train_pi0_bdt.py --period fa18_inb --sample-tag full

A larger category cap:

    python train_pi0_bdt.py --period fa18_inb --max-events-per-category 50000

Disable the exploratory Mx2(epgamma) selection:

    python train_pi0_bdt.py --period fa18_inb --disable-mx2-cut
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import joblib
import matplotlib.pyplot as plt
import numpy as np
import uproot

from scipy.optimize import curve_fit

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)


# =============================================================================
# Configuration
# =============================================================================

PROGRAM_START = time.perf_counter()

FEATURES: List[str] = [
    "p2_p",
    "p2_theta",
    "open_angle_ep2",
    "open_angle_p1p2",
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    "Emiss2",
    "theta_gamma_gamma",
    "pTmiss",
]

FEATURE_LABELS: Dict[str, str] = {
    "p2_p": r"$p_\gamma$ (GeV)",
    "p2_theta": r"$\theta_\gamma$ (deg)",
    "open_angle_ep2": r"$\angle(e',\gamma)$ (deg)",
    "open_angle_p1p2": r"$\angle(p',\gamma)$ (deg)",
    "Mx2": r"$M_X^2(ep\gamma)$ (GeV$^2$)",
    "Mx2_1": r"$M_X^2(ep)$ (GeV$^2$)",
    "Mx2_2": r"$M_X^2(e\gamma)$ (GeV$^2$)",
    "Emiss2": r"$E_{\rm miss}$ (GeV)",
    "theta_gamma_gamma": r"$\theta_{\gamma\gamma}$ (deg)",
    "pTmiss": r"$p_T^{\rm miss}$ (GeV)",
}

OPEN_ANGLE_MIN_DEG = 7.0
DEFAULT_MX2_MIN = -0.5
DEFAULT_MX2_MAX = 0.5

# Diagnostic only: study classifier response in a low-Mx2(ep) region below
# the eta peak. This does NOT change the training sample.
ETA_MASS_GEV = 0.547862
DEFAULT_LOW_MX2_EP_MAX = 0.25

# Reconstructed pi0 mass-fit diagnostics in real epgammagamma data.
PI0_MASS_GEV = 0.134977
PI0_FIT_MIN = 0.070
PI0_FIT_MAX = 0.200
PI0_FIT_BINS = 65
PI0_FIT_NSIGMA = 2.0
PAIR_SCORE_THRESHOLDS = np.arange(0.0, 0.81, 0.1)

DEFAULT_BDT_PARAMS = {
    "n_estimators": 100,
    "learning_rate": 0.05,
    "max_depth": 2,
    "min_samples_leaf": 20,
    "subsample": 0.80,
    "random_state": 42,
    "verbose": 1,
}

REGIONS = {
    "FT": 0,
    "FD": 1,
}

BASE_DIRS = {
    "fa18_inb": Path(
        "/work/clas12/thayward/pi0_BDT/training_sample/fa18_inb"
    ),
}

# Fine-grained truth categories.
CATEGORY_ORDER = [
    "aaogen_pi0_gamma",
    "clasdis_pi0_gamma",
    "dvcsgen_true_gamma",
    "clasdis_eta_gamma",
    "clasdis_other_true_gamma",
    "dvcsgen_nonphoton_match",
    "aaogen_nonphoton_match",
    "clasdis_nonphoton_match",
    "dvcsgen_unmatched",
    "aaogen_unmatched",
    "clasdis_unmatched",
]

CATEGORY_LABELS = {
    "aaogen_pi0_gamma": r"AAOgen true $\pi^0\to\gamma$",
    "clasdis_pi0_gamma": r"CLASDIS true $\pi^0\to\gamma$",
    "dvcsgen_true_gamma": r"DVCSgen true $\gamma$",
    "clasdis_eta_gamma": r"CLASDIS true $\eta\to\gamma$",
    "clasdis_other_true_gamma": r"CLASDIS other true $\gamma$",
    "dvcsgen_nonphoton_match": r"DVCSgen matched non-$\gamma$",
    "aaogen_nonphoton_match": r"AAOgen matched non-$\gamma$",
    "clasdis_nonphoton_match": r"CLASDIS matched non-$\gamma$",
    "dvcsgen_unmatched": r"DVCSgen unmatched",
    "aaogen_unmatched": r"AAOgen unmatched",
    "clasdis_unmatched": r"CLASDIS unmatched",
}

POSITIVE_CATEGORIES = {
    "aaogen_pi0_gamma",
    "clasdis_pi0_gamma",
}

NEGATIVE_MACRO_GROUPS = {
    "dvcs_true_gamma": ["dvcsgen_true_gamma"],
    "clasdis_nonpi0_true_gamma": [
        "clasdis_eta_gamma",
        "clasdis_other_true_gamma",
    ],
    "nonphoton_match": [
        "dvcsgen_nonphoton_match",
        "aaogen_nonphoton_match",
        "clasdis_nonphoton_match",
    ],
    "unmatched": [
        "dvcsgen_unmatched",
        "aaogen_unmatched",
        "clasdis_unmatched",
    ],
}

# Nominal total training weight assigned to each negative macro-category.
# These sum to 0.5, matching the total positive-class weight.
NEGATIVE_MACRO_TARGETS = {
    "dvcs_true_gamma": 0.300,
    "clasdis_nonpi0_true_gamma": 0.050,
    "nonphoton_match": 0.075,
    "unmatched": 0.075,
}

CORE_PLOT_CATEGORIES = [
    "aaogen_pi0_gamma",
    "clasdis_pi0_gamma",
    "dvcsgen_true_gamma",
    "clasdis_eta_gamma",
    "clasdis_other_true_gamma",
]

ARTIFACT_PLOT_CATEGORIES = [
    "dvcsgen_nonphoton_match",
    "aaogen_nonphoton_match",
    "clasdis_nonphoton_match",
    "dvcsgen_unmatched",
    "aaogen_unmatched",
    "clasdis_unmatched",
]

EPG_REQUIRED = sorted(
    set(
        FEATURES
        + [
            "detector2",
            "matching_gamma_pid",
            "gamma_mcindex",
            "gamma_parent_pid",
        ]
    )
)

DATA_EPG_REQUIRED = sorted(
    set(
        FEATURES
        + [
            "detector2",
        ]
    )
)

# Branches needed to create the two single-photon feature vectors from the
# e'p'gammagammaX tree.
EPGG_REQUIRED_COMMON = [
    "detector_gamma1",
    "detector_gamma2",
    "p2_p",
    "p2_theta",
    "p3_p",
    "p3_theta",
    "open_angle_ep2",
    "open_angle_ep3",
    "open_angle_p1p2",
    "open_angle_p1p3",
    "gamma1_epgamma_Mx2",
    "gamma1_ep_Mx2",
    "gamma1_egamma_Mx2",
    "gamma1_Emiss2",
    "gamma1_theta_gamma_gamma",
    "gamma1_pTmiss",
    "gamma2_epgamma_Mx2",
    "gamma2_ep_Mx2",
    "gamma2_egamma_Mx2",
    "gamma2_Emiss2",
    "gamma2_theta_gamma_gamma",
    "gamma2_pTmiss",
    "Mh_gammagamma",
]

EPGG_TRUTH = [
    "matching_gamma1_pid",
    "gamma1_mcindex",
    "gamma1_parent_index",
    "gamma1_parent_pid",
    "matching_gamma2_pid",
    "gamma2_mcindex",
    "gamma2_parent_index",
    "gamma2_parent_pid",
]


# =============================================================================
# Small data containers
# =============================================================================

@dataclass
class CategorySample:
    name: str
    X: np.ndarray
    y: int
    source: str


@dataclass
class SplitData:
    X: np.ndarray
    y: np.ndarray
    categories: np.ndarray
    weights: np.ndarray


@dataclass
class RegionResult:
    name: str
    detector_value: int
    model: GradientBoostingClassifier
    test_auc: float
    core_auc: float
    balanced_accuracy: float
    category_arrays: Dict[str, np.ndarray]


# =============================================================================
# Utilities
# =============================================================================

def progress(message: str) -> None:
    elapsed = time.perf_counter() - PROGRAM_START
    print(f"[+{elapsed:8.1f} s] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Train FT/FD pi0-origin BDTs from direct truth labels in the new "
            "epgammaX processing samples."
        )
    )

    parser.add_argument(
        "--period",
        default="fa18_inb",
        choices=sorted(BASE_DIRS),
        help="Run period. Current truth-enabled defaults are available for fa18_inb.",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=None,
        help="Override the directory containing the processed ROOT samples.",
    )
    parser.add_argument(
        "--sample-tag",
        default="temp",
        help=(
            "Filename suffix between X and .root. Default: temp. "
            "Use 'full' for uncapped files with no suffix."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/pi0_bdt_truth"),
        help="Top-level output directory.",
    )

    parser.add_argument(
        "--max-events-per-category",
        type=int,
        default=20000,
        help=(
            "Maximum number of events retained from each fine truth category "
            "per FT/FD model. <=0 means no category cap."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed.",
    )

    parser.add_argument(
        "--open-angle-min",
        type=float,
        default=OPEN_ANGLE_MIN_DEG,
        help="Minimum e'-gamma opening angle in degrees.",
    )
    parser.add_argument(
        "--mx2-min",
        type=float,
        default=DEFAULT_MX2_MIN,
        help="Lower Mx2(epgamma) selection boundary.",
    )
    parser.add_argument(
        "--mx2-max",
        type=float,
        default=DEFAULT_MX2_MAX,
        help="Upper Mx2(epgamma) selection boundary.",
    )
    parser.add_argument(
        "--disable-mx2-cut",
        action="store_true",
        help="Do not apply the exploratory Mx2(epgamma) window.",
    )
    parser.add_argument(
        "--low-mx2-ep-max",
        type=float,
        default=DEFAULT_LOW_MX2_EP_MAX,
        help=(
            "Upper Mx2(ep) boundary used ONLY for the post-training low-Mx2(ep) "
            "score diagnostic. Default 0.25 GeV^2, below the eta peak at "
            f"m_eta^2={ETA_MASS_GEV**2:.3f} GeV^2."
        ),
    )

    parser.add_argument(
        "--n-estimators",
        type=int,
        default=DEFAULT_BDT_PARAMS["n_estimators"],
        help="Number of boosting trees.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=DEFAULT_BDT_PARAMS["max_depth"],
        help="Maximum depth of each shallow boosting tree.",
    )
    parser.add_argument(
        "--learning-rate",
        type=float,
        default=DEFAULT_BDT_PARAMS["learning_rate"],
        help="Gradient-boosting learning rate.",
    )
    parser.add_argument(
        "--min-samples-leaf",
        type=int,
        default=DEFAULT_BDT_PARAMS["min_samples_leaf"],
        help="Minimum samples in a tree leaf.",
    )
    parser.add_argument(
        "--subsample",
        type=float,
        default=DEFAULT_BDT_PARAMS["subsample"],
        help="Training subsample fraction for stochastic gradient boosting.",
    )

    return parser.parse_args()


def sample_suffix(sample_tag: str) -> str:
    tag = sample_tag.strip()

    if tag.lower() in {"", "full", "none"}:
        return ""
    #endif

    if tag.startswith("_"):
        return tag
    #endif

    return "_" + tag


def build_paths(args: argparse.Namespace) -> Dict[str, Path]:
    base = args.base_dir if args.base_dir is not None else BASE_DIRS[args.period]
    suffix = sample_suffix(args.sample_tag)

    return {
        "dvcsgen_epg": base / f"dvcsgen_rga_fa18_inb_epgammaX{suffix}.root",
        "aaogen_epg": base / f"aaogen_rga_fa18_inb_epgammaX{suffix}.root",
        "clasdis_epg": base / f"clasdis_rga_fa18_inb_epgammaX{suffix}.root",
        "data_epg": base / f"data_rga_fa18_inb_epgammaX{suffix}.root",
        "aaogen_epgg": base / f"aaogen_rga_fa18_inb_epgammagammaX{suffix}.root",
        "clasdis_epgg": base / f"clasdis_rga_fa18_inb_epgammagammaX{suffix}.root",
        "data_epgg": base / f"data_rga_fa18_inb_epgammagammaX{suffix}.root",
    }


def first_tree_name(path: Path) -> str:
    with uproot.open(path) as root_file:
        for key, obj in root_file.items():
            if hasattr(obj, "arrays") and hasattr(obj, "num_entries"):
                return key.split(";")[0]
            #endif
        #endfor
    raise RuntimeError(f"No TTree-like object found in {path}")


def read_tree_arrays(
    path: Path,
    branches: Sequence[str],
) -> Dict[str, np.ndarray]:
    tree_name = first_tree_name(path)

    progress(f"LOAD: {path}")
    start = time.perf_counter()

    with uproot.open(path) as root_file:
        tree = root_file[tree_name]
        available = {name.split(";")[0] for name in tree.keys()}
        missing = [name for name in branches if name not in available]

        if missing:
            raise RuntimeError(
                f"Missing required branches in {path}: {missing}"
            )
        #endif

        raw = tree.arrays(list(branches), library="np")
    #endwith

    arrays = {name: np.asarray(raw[name]) for name in branches}
    n = len(next(iter(arrays.values()))) if arrays else 0

    progress(
        f"LOAD COMPLETE: {n:,} rows, {len(branches)} branches "
        f"({time.perf_counter() - start:.1f} s)"
    )
    return arrays


def ensure_files(paths: Mapping[str, Path]) -> None:
    missing = [path for path in paths.values() if not path.exists()]

    if missing:
        print("ERROR: missing required input files:")
        for path in missing:
            print(f"  {path}")
        #endfor
        raise SystemExit(1)
    #endif


def maybe_deg(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]

    if len(finite) == 0:
        return arr
    #endif

    # The processing analyzers normally store angular variables in radians for
    # theta/phi, while the open-angle observables are already degrees.
    if np.nanpercentile(np.abs(finite), 99.0) <= 3.5:
        return np.degrees(arr)
    #endif

    return arr


def feature_matrix(
    arrays: Mapping[str, np.ndarray],
) -> np.ndarray:
    columns = []

    for feature in FEATURES:
        values = np.asarray(arrays[feature], dtype=float)

        if feature == "p2_theta":
            values = maybe_deg(values)
        elif feature == "theta_gamma_gamma":
            values = maybe_deg(values)
        #endif

        columns.append(values)
    #endfor

    return np.column_stack(columns)


def common_candidate_mask(
    arrays: Mapping[str, np.ndarray],
    args: argparse.Namespace,
    X: Optional[np.ndarray] = None,
) -> np.ndarray:
    n = len(arrays["detector2"])
    mask = np.ones(n, dtype=bool)

    open_angle = np.asarray(arrays["open_angle_ep2"], dtype=float)
    mask &= np.isfinite(open_angle)
    mask &= open_angle > args.open_angle_min

    mx2 = np.asarray(arrays["Mx2"], dtype=float)

    if not args.disable_mx2_cut:
        mask &= np.isfinite(mx2)
        mask &= mx2 > args.mx2_min
        mask &= mx2 < args.mx2_max
    #endif

    if X is None:
        X = feature_matrix(arrays)
    #endif

    mask &= np.all(np.isfinite(X), axis=1)

    return mask


def category_masks(
    source: str,
    arrays: Mapping[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    matching = np.asarray(arrays["matching_gamma_pid"], dtype=np.int64)
    mcindex = np.asarray(arrays["gamma_mcindex"], dtype=np.int64)
    parent = np.asarray(arrays["gamma_parent_pid"], dtype=np.int64)

    matched = mcindex >= 0
    true_gamma = matched & (matching == 22)
    nonphoton = matched & (matching != 22)
    unmatched = mcindex < 0

    if source == "dvcsgen":
        return {
            "dvcsgen_true_gamma": true_gamma,
            "dvcsgen_nonphoton_match": nonphoton,
            "dvcsgen_unmatched": unmatched,
        }
    #endif

    if source == "aaogen":
        return {
            "aaogen_pi0_gamma": true_gamma,
            "aaogen_nonphoton_match": nonphoton,
            "aaogen_unmatched": unmatched,
        }
    #endif

    if source == "clasdis":
        return {
            "clasdis_pi0_gamma": true_gamma & (parent == 111),
            "clasdis_eta_gamma": true_gamma & (parent == 221),
            "clasdis_other_true_gamma": (
                true_gamma
                & (parent != 111)
                & (parent != 221)
            ),
            "clasdis_nonphoton_match": nonphoton,
            "clasdis_unmatched": unmatched,
        }
    #endif

    raise ValueError(f"Unknown source: {source}")


def capped_indices(
    indices: np.ndarray,
    cap: int,
    rng: np.random.Generator,
) -> np.ndarray:
    indices = np.asarray(indices, dtype=np.int64)

    if cap <= 0 or len(indices) <= cap:
        return indices
    #endif

    return np.sort(rng.choice(indices, size=cap, replace=False))


def robust_hist_range(
    arrays: Iterable[np.ndarray],
    qlo: float = 0.005,
    qhi: float = 0.995,
) -> Optional[Tuple[float, float]]:
    clean = []

    for values in arrays:
        vals = np.asarray(values, dtype=float)
        vals = vals[np.isfinite(vals)]

        if len(vals) > 0:
            clean.append(vals)
        #endif
    #endfor

    if not clean:
        return None
    #endif

    merged = np.concatenate(clean)
    lo, hi = np.quantile(merged, [qlo, qhi])

    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        return None
    #endif

    return float(lo), float(hi)


def normalized_hist(
    ax: plt.Axes,
    values: np.ndarray,
    *,
    label: str,
    bins: int = 70,
    hist_range: Optional[Tuple[float, float]] = None,
    linewidth: float = 1.5,
    linestyle: str = "-",
    color: Optional[str] = None,
) -> None:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]

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
        color=color,
        label=label,
    )


def save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    progress(f"PLOT: {path}")


# =============================================================================
# Category construction and splitting
# =============================================================================

def build_region_categories(
    arrays_by_source: Mapping[str, Mapping[str, np.ndarray]],
    region_name: str,
    detector_value: int,
    args: argparse.Namespace,
) -> Dict[str, CategorySample]:
    rng = np.random.default_rng(args.seed + detector_value)
    result: Dict[str, CategorySample] = {}

    for source in ["dvcsgen", "aaogen", "clasdis"]:
        arrays = arrays_by_source[source]
        X_all = feature_matrix(arrays)
        base = common_candidate_mask(arrays, args, X=X_all)
        base &= np.asarray(arrays["detector2"], dtype=np.int64) == detector_value

        masks = category_masks(source, arrays)

        for category, truth_mask in masks.items():
            indices = np.flatnonzero(base & truth_mask)
            indices = capped_indices(
                indices,
                args.max_events_per_category,
                rng,
            )

            y = 1 if category in POSITIVE_CATEGORIES else 0

            result[category] = CategorySample(
                name=category,
                X=X_all[indices],
                y=y,
                source=source,
            )

            progress(
                f"{region_name}: {category:28s} "
                f"N={len(indices):,}"
            )
        #endfor
    #endfor

    return result


def split_one_category(
    sample: CategorySample,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(sample.X)

    if n == 0:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty, empty
    #endif

    perm = rng.permutation(n)

    if n < 10:
        # Tiny categories should not drive evaluation metrics. Keep them in
        # training so the model at least sees the topology.
        return perm, np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    #endif

    n_train = max(1, int(round(0.70 * n)))
    n_val = max(1, int(round(0.15 * n)))

    if n_train + n_val >= n:
        n_val = max(1, n - n_train - 1)
    #endif

    n_test = n - n_train - n_val

    if n_test < 1:
        n_test = 1
        n_train = max(1, n - n_val - n_test)
    #endif

    train = perm[:n_train]
    val = perm[n_train:n_train + n_val]
    test = perm[n_train + n_val:]

    return train, val, test


def training_weight_targets(
    categories: Mapping[str, CategorySample],
) -> Dict[str, float]:
    """
    Return desired total training weight for each fine category.

    The positive and negative classes each receive total target weight 0.5.

    Positive:
        AAOgen true pi0 photons and CLASDIS true pi0 photons share the positive
        weight equally when both are available.

    Negative:
        Macro-category priorities are fixed by NEGATIVE_MACRO_TARGETS. If a
        whole macro-category is absent, the available macro targets are
        renormalized back to total negative weight 0.5 while preserving their
        nominal relative priorities.

        Within a macro-category, its weight is divided among the available fine
        categories in proportion to their retained event populations. Thus
        rare eta/other/artifact categories are represented, but are not given
        the same total influence merely because they have distinct truth labels.
    """
    targets: Dict[str, float] = {}

    available_positive = [
        name
        for name in ["aaogen_pi0_gamma", "clasdis_pi0_gamma"]
        if name in categories and len(categories[name].X) > 0
    ]

    if available_positive:
        per_positive = 0.5 / len(available_positive)

        for name in available_positive:
            targets[name] = per_positive
        #endfor
    #endif

    available_macros = []

    for macro, fine_names in NEGATIVE_MACRO_GROUPS.items():
        available_fine = [
            name
            for name in fine_names
            if name in categories and len(categories[name].X) > 0
        ]

        if available_fine:
            available_macros.append((macro, available_fine))
        #endif
    #endfor

    if available_macros:
        available_nominal_total = sum(
            NEGATIVE_MACRO_TARGETS[macro]
            for macro, _ in available_macros
        )

        for macro, fine_names in available_macros:
            macro_target = (
                0.5
                * NEGATIVE_MACRO_TARGETS[macro]
                / available_nominal_total
            )

            counts = np.asarray(
                [len(categories[name].X) for name in fine_names],
                dtype=float,
            )
            count_total = float(np.sum(counts))

            if count_total <= 0.0:
                continue
            #endif

            for name, count in zip(fine_names, counts):
                targets[name] = macro_target * float(count) / count_total
            #endfor
        #endfor
    #endif

    return targets


def build_splits(
    categories: Mapping[str, CategorySample],
    seed: int,
) -> Tuple[SplitData, SplitData, SplitData]:
    rng = np.random.default_rng(seed)
    indices_by_category = {}

    for name, sample in categories.items():
        indices_by_category[name] = split_one_category(sample, rng)
    #endfor

    targets = training_weight_targets(categories)

    def assemble(split_index: int, weighted: bool) -> SplitData:
        X_parts = []
        y_parts = []
        category_parts = []
        weight_parts = []

        for name in CATEGORY_ORDER:
            if name not in categories:
                continue
            #endif

            sample = categories[name]
            indices = indices_by_category[name][split_index]

            if len(indices) == 0:
                continue
            #endif

            X_parts.append(sample.X[indices])
            y_parts.append(
                np.full(len(indices), sample.y, dtype=np.int8)
            )
            category_parts.append(
                np.full(len(indices), name, dtype=object)
            )

            if weighted:
                total_target = targets.get(name, 0.0)
                per_event = (
                    total_target / len(indices)
                    if len(indices) > 0
                    else 0.0
                )
                weight_parts.append(
                    np.full(len(indices), per_event, dtype=float)
                )
            else:
                weight_parts.append(
                    np.ones(len(indices), dtype=float)
                )
            #endif
        #endfor

        if not X_parts:
            raise RuntimeError("No events available for requested split.")
        #endif

        X = np.concatenate(X_parts, axis=0)
        y = np.concatenate(y_parts, axis=0)
        cats = np.concatenate(category_parts, axis=0)
        weights = np.concatenate(weight_parts, axis=0)

        # Normalize mean weight to one. This preserves all relative weights
        # while keeping sklearn's effective weight scale numerically ordinary.
        if weighted and np.sum(weights) > 0:
            weights *= len(weights) / np.sum(weights)
        #endif

        order = rng.permutation(len(X))

        return SplitData(
            X=X[order],
            y=y[order],
            categories=cats[order],
            weights=weights[order],
        )

    train = assemble(0, weighted=True)
    validation = assemble(1, weighted=False)
    test = assemble(2, weighted=False)

    return train, validation, test


# =============================================================================
# Plotting: counts and features
# =============================================================================

def plot_global_truth_counts(
    categories_by_region: Mapping[str, Mapping[str, CategorySample]],
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 7.5))

    for ax, region_name in zip(axes, ["FT", "FD"]):
        categories = categories_by_region[region_name]
        names = [
            name
            for name in CATEGORY_ORDER
            if name in categories
        ]
        counts = [len(categories[name].X) for name in names]
        labels = [CATEGORY_LABELS[name] for name in names]

        y = np.arange(len(names))
        ax.barh(y, counts)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Candidates retained after common selections/cap")
        ax.set_title(region_name)
        ax.grid(axis="x", alpha=0.2)
    #endfor

    fig.suptitle("Truth-labelled e'p'gammaX category populations", y=0.995)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.965])
    save_figure(fig, output_path)


def plot_region_training_counts(
    categories: Mapping[str, CategorySample],
    output_path: Path,
    region_name: str,
) -> None:
    targets = training_weight_targets(categories)

    names = [
        name
        for name in CATEGORY_ORDER
        if name in categories and len(categories[name].X) > 0
    ]

    fig, axes = plt.subplots(1, 2, figsize=(16.0, 7.5))

    ax = axes[0]
    counts = [len(categories[name].X) for name in names]
    y = np.arange(len(names))
    ax.barh(y, counts)
    ax.set_yticks(y)
    ax.set_yticklabels([CATEGORY_LABELS[name] for name in names], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Candidates")
    ax.set_title("Fine truth-category counts")
    ax.grid(axis="x", alpha=0.2)

    ax = axes[1]
    weights = [targets.get(name, 0.0) for name in names]
    ax.barh(y, weights)
    ax.set_yticks(y)
    ax.set_yticklabels([CATEGORY_LABELS[name] for name in names], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Target total training weight")
    ax.set_title("Controlled training composition")
    ax.grid(axis="x", alpha=0.2)

    fig.suptitle(
        f"{region_name}: training composition before train/validation/test split",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.965])
    save_figure(fig, output_path)


def plot_feature_distributions(
    category_arrays: Mapping[str, np.ndarray],
    data_X: np.ndarray,
    category_names: Sequence[str],
    output_path: Path,
    title: str,
    *,
    include_data: bool,
) -> None:
    fig, axes = plt.subplots(2, 5, figsize=(20.0, 8.2))
    axes = axes.ravel()

    for i, (ax, feature) in enumerate(zip(axes, FEATURES)):
        values_for_range = []

        for category in category_names:
            if category in category_arrays and len(category_arrays[category]) > 0:
                values_for_range.append(category_arrays[category][:, i])
            #endif
        #endfor

        if include_data and len(data_X) > 0:
            values_for_range.append(data_X[:, i])
        #endif

        hist_range = robust_hist_range(values_for_range)

        # Keep the known narrow Mx2(epgamma) selection visually stable.
        if feature == "Mx2":
            hist_range = robust_hist_range(values_for_range, qlo=0.001, qhi=0.999)
        #endif

        for category in category_names:
            if category not in category_arrays:
                continue
            #endif

            values = category_arrays[category]

            if len(values) == 0:
                continue
            #endif

            normalized_hist(
                ax,
                values[:, i],
                label=f"{CATEGORY_LABELS[category]} (N={len(values):,})",
                hist_range=hist_range,
            )
        #endfor

        if include_data and len(data_X) > 0:
            normalized_hist(
                ax,
                data_X[:, i],
                label=f"Real data (N={len(data_X):,})",
                hist_range=hist_range,
                linewidth=2.0,
                linestyle="--",
                color="black",
            )
        #endif

        ax.set_xlabel(FEATURE_LABELS[feature])
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.2)
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()

    # Reserve separate, non-overlapping vertical bands for the figure title
    # and the shared legend.  In particular, do not let fig.legend() occupy
    # the same top-of-canvas region as fig.suptitle().
    fig.suptitle(title, y=0.985)

    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.935),
            ncol=min(3, len(labels)),
            fontsize=8,
            frameon=True,
        )
    #endif

    # Keep the subplot grid entirely below the title+legend header band.
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.82])
    save_figure(fig, output_path)


# =============================================================================
# Model training and evaluation
# =============================================================================

def train_model(
    train: SplitData,
    args: argparse.Namespace,
    region_name: str,
) -> GradientBoostingClassifier:
    params = dict(DEFAULT_BDT_PARAMS)
    params.update(
        {
            "n_estimators": args.n_estimators,
            "learning_rate": args.learning_rate,
            "max_depth": args.max_depth,
            "min_samples_leaf": args.min_samples_leaf,
            "subsample": args.subsample,
            "random_state": args.seed,
        }
    )

    model = GradientBoostingClassifier(**params)

    progress(
        f"{region_name} BDT FIT: {len(train.X):,} training rows x "
        f"{train.X.shape[1]} features; "
        f"{args.n_estimators} trees, max_depth={args.max_depth}"
    )

    start = time.perf_counter()
    model.fit(
        train.X,
        train.y,
        sample_weight=train.weights,
    )
    progress(
        f"{region_name} BDT FIT COMPLETE: "
        f"{time.perf_counter() - start:.1f} s"
    )

    return model


def score(model: GradientBoostingClassifier, X: np.ndarray) -> np.ndarray:
    if len(X) == 0:
        return np.empty(0, dtype=float)
    #endif

    return model.predict_proba(X)[:, 1]


def compute_core_auc(
    model: GradientBoostingClassifier,
    test: SplitData,
) -> float:
    positive = np.isin(
        test.categories,
        ["aaogen_pi0_gamma", "clasdis_pi0_gamma"],
    )
    negative = np.isin(
        test.categories,
        ["dvcsgen_true_gamma"],
    )
    mask = positive | negative

    if np.sum(mask) < 2 or len(np.unique(test.y[mask])) < 2:
        return float("nan")
    #endif

    return float(
        roc_auc_score(
            test.y[mask],
            score(model, test.X[mask]),
        )
    )


def plot_roc(
    model: GradientBoostingClassifier,
    test: SplitData,
    output_path: Path,
    region_name: str,
) -> Tuple[float, float]:
    scores = score(model, test.X)

    overall_auc = float(roc_auc_score(test.y, scores))
    fpr, tpr, _ = roc_curve(test.y, scores)

    positive_core = np.isin(
        test.categories,
        ["aaogen_pi0_gamma", "clasdis_pi0_gamma"],
    )
    negative_core = test.categories == "dvcsgen_true_gamma"
    core_mask = positive_core | negative_core

    core_auc = float("nan")
    core_fpr = np.empty(0)
    core_tpr = np.empty(0)

    if (
        np.sum(core_mask) >= 2
        and len(np.unique(test.y[core_mask])) == 2
    ):
        core_scores = scores[core_mask]
        core_auc = float(
            roc_auc_score(test.y[core_mask], core_scores)
        )
        core_fpr, core_tpr, _ = roc_curve(
            test.y[core_mask],
            core_scores,
        )
    #endif

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.plot(
        fpr,
        tpr,
        linewidth=2.0,
        label=f"All truth categories (AUC={overall_auc:.4f})",
    )

    if len(core_fpr) > 0:
        ax.plot(
            core_fpr,
            core_tpr,
            linewidth=2.0,
            linestyle="--",
            label=(
                r"$\pi^0$ daughters vs genuine DVCS/BH $\gamma$ "
                f"(AUC={core_auc:.4f})"
            ),
        )
    #endif

    ax.plot([0, 1], [0, 1], linestyle=":", linewidth=1.2)
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.02)
    ax.set_title(f"{region_name}: held-out test ROC")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9)
    fig.tight_layout()
    save_figure(fig, output_path)

    return overall_auc, core_auc


def plot_core_operating_points(
    model: GradientBoostingClassifier,
    test: SplitData,
    output_path: Path,
    region_name: str,
) -> None:
    """
    Show the physically interpretable core operating curve.

    x = fraction of genuine DVCS/BH photons that would be falsely tagged pi0
        when requiring BDT score >= threshold.
    y = fraction of true pi0-daughter photons retained by the same threshold.

    The labelled points are evaluated DIRECTLY at the stated thresholds rather
    than snapped to the nearest internal sklearn ROC threshold.
    """
    positive = np.isin(
        test.categories,
        ["aaogen_pi0_gamma", "clasdis_pi0_gamma"],
    )
    negative = test.categories == "dvcsgen_true_gamma"
    mask = positive | negative

    if np.sum(mask) < 2 or len(np.unique(test.y[mask])) < 2:
        return
    #endif

    y_true = test.y[mask]
    scores = score(model, test.X[mask])

    fpr, tpr, _ = roc_curve(y_true, scores)

    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    ax.plot(
        100.0 * fpr,
        100.0 * tpr,
        linewidth=2.0,
        label=r"True $\pi^0$ daughter vs genuine DVCS/BH $\gamma$",
    )

    for threshold in np.arange(0.1, 1.0, 0.1):
        predicted_positive = scores >= threshold

        pos_eff = (
            np.sum(predicted_positive & (y_true == 1))
            / np.sum(y_true == 1)
        )
        dvcs_mistag = (
            np.sum(predicted_positive & (y_true == 0))
            / np.sum(y_true == 0)
        )

        ax.scatter(
            [100.0 * dvcs_mistag],
            [100.0 * pos_eff],
            s=34,
            zorder=4,
        )
        ax.annotate(
            f"{threshold:.1f}",
            (100.0 * dvcs_mistag, 100.0 * pos_eff),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )
    #endfor

    ax.set_xlabel("Genuine DVCS/BH photon mis-tag rate [%]")
    ax.set_ylabel(r"True $\pi^0$-daughter efficiency [%]")
    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(0.0, 101.0)
    ax.set_title(
        f"{region_name}: core BDT operating points\n"
        r"label = requirement BDT score $\geq$ threshold"
    )
    ax.grid(alpha=0.2)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_scores_by_category(
    model: GradientBoostingClassifier,
    category_arrays: Mapping[str, np.ndarray],
    output_path: Path,
    region_name: str,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 9.5), sharex=True)

    positive_and_true = [
        "aaogen_pi0_gamma",
        "clasdis_pi0_gamma",
        "dvcsgen_true_gamma",
        "clasdis_eta_gamma",
        "clasdis_other_true_gamma",
    ]

    artifacts = ARTIFACT_PLOT_CATEGORIES

    for ax, categories, subtitle in [
        (
            axes[0],
            positive_and_true,
            "True photons: pi0 daughters versus non-pi0 photons",
        ),
        (
            axes[1],
            artifacts,
            "Reconstruction-artifact categories",
        ),
    ]:
        for category in categories:
            X = category_arrays.get(category)

            if X is None or len(X) == 0:
                continue
            #endif

            values = score(model, X)

            normalized_hist(
                ax,
                values,
                label=f"{CATEGORY_LABELS[category]} (N={len(values):,})",
                bins=60,
                hist_range=(0.0, 1.0),
            )
        #endfor

        ax.set_ylabel("Normalized density")
        ax.set_title(subtitle)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8, ncol=2)
    #endfor

    axes[1].set_xlabel(r"BDT $\pi^0$ score")
    fig.suptitle(
        f"{region_name}: classifier response by direct MC truth category",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.965])
    save_figure(fig, output_path)


def plot_core_scores_and_data(
    model: GradientBoostingClassifier,
    category_arrays: Mapping[str, np.ndarray],
    data_X: np.ndarray,
    output_path: Path,
    region_name: str,
) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 6.7))

    categories = [
        "aaogen_pi0_gamma",
        "clasdis_pi0_gamma",
        "dvcsgen_true_gamma",
        "clasdis_eta_gamma",
        "clasdis_other_true_gamma",
    ]

    for category in categories:
        X = category_arrays.get(category)

        if X is None or len(X) == 0:
            continue
        #endif

        normalized_hist(
            ax,
            score(model, X),
            label=f"{CATEGORY_LABELS[category]} (N={len(X):,})",
            bins=60,
            hist_range=(0.0, 1.0),
        )
    #endfor

    if len(data_X) > 0:
        normalized_hist(
            ax,
            score(model, data_X),
            label=f"Ordinary real data (N={len(data_X):,})",
            bins=60,
            hist_range=(0.0, 1.0),
            linewidth=2.2,
            linestyle="--",
            color="black",
        )
    #endif

    ax.set_xlabel(r"BDT $\pi^0$ score")
    ax.set_ylabel("Normalized density")
    ax.set_xlim(0.0, 1.0)
    ax.set_title(f"{region_name}: core truth samples and real data")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, output_path)



def plot_low_mx2ep_score_shift(
    model: GradientBoostingClassifier,
    category_arrays: Mapping[str, np.ndarray],
    data_X: np.ndarray,
    output_path: Path,
    region_name: str,
    mx2_ep_max: float,
) -> None:
    """
    Compare BDT-score distributions before/after a LOW Mx2(ep) requirement.

    This is deliberately a post-training diagnostic only. Mx2(ep) is itself a
    BDT input feature, so any score shift is NOT an independent validation of
    the classifier. The useful question is whether the inclusive CLASDIS pi0
    population becomes more AAOgen-like in an exclusive-like low-missing-mass
    region.
    """
    feature_index = FEATURES.index("Mx2_1")

    panels = [
        ("aaogen_pi0_gamma", r"AAOgen true $\pi^0\to\gamma$"),
        ("clasdis_pi0_gamma", r"CLASDIS true $\pi^0\to\gamma$"),
        ("dvcsgen_true_gamma", r"DVCSgen genuine $\gamma$"),
        ("data", "Real data"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.0), sharex=True)
    bins = np.linspace(0.0, 1.0, 61)

    for ax, (category, label) in zip(axes.flat, panels):
        if category == "data":
            X = data_X
        else:
            X = category_arrays.get(category, np.empty((0, len(FEATURES))))
        #endif

        if len(X) == 0:
            ax.set_axis_off()
            continue
        #endif

        scores_all = score(model, X)
        low_mask = (
            np.isfinite(X[:, feature_index])
            & (X[:, feature_index] < mx2_ep_max)
        )
        scores_low = scores_all[low_mask]

        if len(scores_all) > 0:
            ax.hist(
                scores_all,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=1.7,
                label=(
                    f"All selected: N={len(scores_all):,}, "
                    f"mean={np.mean(scores_all):.3f}"
                ),
            )
        #endif

        if len(scores_low) > 0:
            retained = 100.0 * len(scores_low) / len(scores_all)
            ax.hist(
                scores_low,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=1.7,
                label=(
                    rf"$M_X^2(ep)<{mx2_ep_max:.2f}$: "
                    f"N={len(scores_low):,} ({retained:.1f}%), "
                    f"mean={np.mean(scores_low):.3f}"
                ),
            )
        #endif

        ax.set_title(label)
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    #endfor

    for ax in axes[-1, :]:
        ax.set_xlabel(r"BDT $\pi^0$ score")
    #endfor

    fig.suptitle(
        (
            f"{region_name}: BDT response in a low-$M_X^2(ep)$ region "
            f"(diagnostic only; $m_\\eta^2={ETA_MASS_GEV**2:.3f}$ GeV$^2$)"
        ),
        y=0.985,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    save_figure(fig, output_path)



def plot_feature_importance(
    model: GradientBoostingClassifier,
    output_path: Path,
    region_name: str,
) -> None:
    importance = np.asarray(model.feature_importances_, dtype=float)
    order = np.argsort(importance)

    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    ax.barh(
        np.arange(len(FEATURES)),
        importance[order],
    )
    ax.set_yticks(np.arange(len(FEATURES)))
    ax.set_yticklabels(
        [FEATURE_LABELS[FEATURES[i]] for i in order],
        fontsize=9,
    )
    ax.set_xlabel("Gradient-boosting feature importance")
    ax.set_title(
        f"{region_name}: native tree feature importance\n"
        "(relative split-improvement contribution)"
    )
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_confusion(
    model: GradientBoostingClassifier,
    test: SplitData,
    output_path: Path,
    region_name: str,
) -> float:
    scores = score(model, test.X)
    pred = (scores >= 0.5).astype(np.int8)

    matrix = confusion_matrix(
        test.y,
        pred,
        labels=[0, 1],
        normalize="true",
    )
    balanced = float(
        balanced_accuracy_score(test.y, pred)
    )

    fig, ax = plt.subplots(figsize=(6.7, 5.8))
    image = ax.imshow(
        matrix,
        vmin=0.0,
        vmax=1.0,
        aspect="equal",
    )

    for i in range(2):
        for j in range(2):
            ax.text(
                j,
                i,
                f"{100.0 * matrix[i, j]:.1f}%",
                ha="center",
                va="center",
                fontsize=14,
            )
        #endfor
    #endfor

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Predicted negative", r"Predicted $\pi^0$"])
    ax.set_yticklabels(["True negative", r"True $\pi^0$"])
    ax.set_xlabel("Classifier decision at score = 0.5")
    ax.set_ylabel("Direct MC truth")
    ax.set_title(
        f"{region_name}: pedagogical score=0.5 confusion matrix\n"
        f"balanced accuracy={balanced:.4f}"
    )
    fig.colorbar(image, ax=ax, label="Fraction within true class")
    fig.tight_layout()
    save_figure(fig, output_path)

    return balanced


# =============================================================================
# Data scoring
# =============================================================================

def selected_data_X(
    arrays: Mapping[str, np.ndarray],
    detector_value: int,
    args: argparse.Namespace,
) -> np.ndarray:
    X = feature_matrix(arrays)
    mask = common_candidate_mask(arrays, args, X=X)
    mask &= np.asarray(arrays["detector2"], dtype=np.int64) == detector_value
    return X[mask]


# =============================================================================
# e'p'gammagammaX validation
# =============================================================================

def epgg_photon_feature_matrix(
    arrays: Mapping[str, np.ndarray],
    which: int,
) -> np.ndarray:
    if which == 1:
        mapping = {
            "p2_p": "p2_p",
            "p2_theta": "p2_theta",
            "open_angle_ep2": "open_angle_ep2",
            "open_angle_p1p2": "open_angle_p1p2",
            "Mx2": "gamma1_epgamma_Mx2",
            "Mx2_1": "gamma1_ep_Mx2",
            "Mx2_2": "gamma1_egamma_Mx2",
            "Emiss2": "gamma1_Emiss2",
            "theta_gamma_gamma": "gamma1_theta_gamma_gamma",
            "pTmiss": "gamma1_pTmiss",
        }
    elif which == 2:
        mapping = {
            "p2_p": "p3_p",
            "p2_theta": "p3_theta",
            "open_angle_ep2": "open_angle_ep3",
            "open_angle_p1p2": "open_angle_p1p3",
            "Mx2": "gamma2_epgamma_Mx2",
            "Mx2_1": "gamma2_ep_Mx2",
            "Mx2_2": "gamma2_egamma_Mx2",
            "Emiss2": "gamma2_Emiss2",
            "theta_gamma_gamma": "gamma2_theta_gamma_gamma",
            "pTmiss": "gamma2_pTmiss",
        }
    else:
        raise ValueError("which must be 1 or 2")
    #endif

    columns = []

    for feature in FEATURES:
        values = np.asarray(arrays[mapping[feature]], dtype=float)

        if feature in {"p2_theta", "theta_gamma_gamma"}:
            values = maybe_deg(values)
        #endif

        columns.append(values)
    #endfor

    return np.column_stack(columns)


def epgg_single_mask(
    X: np.ndarray,
    detector: np.ndarray,
    args: argparse.Namespace,
) -> np.ndarray:
    mask = np.all(np.isfinite(X), axis=1)

    open_idx = FEATURES.index("open_angle_ep2")
    mx2_idx = FEATURES.index("Mx2")

    mask &= X[:, open_idx] > args.open_angle_min
    mask &= np.isin(detector, [0, 1])

    if not args.disable_mx2_cut:
        mask &= X[:, mx2_idx] > args.mx2_min
        mask &= X[:, mx2_idx] < args.mx2_max
    #endif

    return mask


def score_epgg_photons(
    arrays: Mapping[str, np.ndarray],
    models: Mapping[str, GradientBoostingClassifier],
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    X1 = epgg_photon_feature_matrix(arrays, 1)
    X2 = epgg_photon_feature_matrix(arrays, 2)

    d1 = np.asarray(arrays["detector_gamma1"], dtype=np.int64)
    d2 = np.asarray(arrays["detector_gamma2"], dtype=np.int64)

    valid1 = epgg_single_mask(X1, d1, args)
    valid2 = epgg_single_mask(X2, d2, args)

    s1 = np.full(len(X1), np.nan, dtype=float)
    s2 = np.full(len(X2), np.nan, dtype=float)

    for region_name, detector_value in REGIONS.items():
        mask1 = valid1 & (d1 == detector_value)
        mask2 = valid2 & (d2 == detector_value)

        if np.any(mask1):
            s1[mask1] = score(models[region_name], X1[mask1])
        #endif

        if np.any(mask2):
            s2[mask2] = score(models[region_name], X2[mask2])
        #endif
    #endfor

    pair_valid = np.isfinite(s1) & np.isfinite(s2)
    return s1, s2, pair_valid


def clasdis_pair_categories(
    arrays: Mapping[str, np.ndarray],
) -> np.ndarray:
    m1 = np.asarray(arrays["matching_gamma1_pid"], dtype=np.int64)
    m2 = np.asarray(arrays["matching_gamma2_pid"], dtype=np.int64)
    mc1 = np.asarray(arrays["gamma1_mcindex"], dtype=np.int64)
    mc2 = np.asarray(arrays["gamma2_mcindex"], dtype=np.int64)
    p1 = np.asarray(arrays["gamma1_parent_pid"], dtype=np.int64)
    p2 = np.asarray(arrays["gamma2_parent_pid"], dtype=np.int64)
    i1 = np.asarray(arrays["gamma1_parent_index"], dtype=np.int64)
    i2 = np.asarray(arrays["gamma2_parent_index"], dtype=np.int64)

    n = len(m1)
    out = np.full(n, "other_mixed", dtype=object)

    unmatched = (mc1 < 0) | (mc2 < 0)
    nonphoton = (
        ((mc1 >= 0) & (m1 != 22))
        | ((mc2 >= 0) & (m2 != 22))
    )

    both_gamma = (
        (mc1 >= 0)
        & (mc2 >= 0)
        & (m1 == 22)
        & (m2 == 22)
    )
    same_parent = (i1 > 0) & (i2 > 0) & (i1 == i2)

    out[unmatched] = "unmatched_daughter"
    out[~unmatched & nonphoton] = "nonphoton_daughter"

    out[
        both_gamma
        & same_parent
        & (p1 == 111)
        & (p2 == 111)
    ] = "same_pi0_parent"

    out[
        both_gamma
        & ~same_parent
        & (p1 == 111)
        & (p2 == 111)
    ] = "different_pi0_parents"

    out[
        both_gamma
        & same_parent
        & (p1 == 221)
        & (p2 == 221)
    ] = "same_eta_parent"

    out[
        both_gamma
        & same_parent
        & (p1 == p2)
        & (p1 != 111)
        & (p1 != 221)
    ] = "same_other_parent"

    return out


PAIR_LABELS = {
    "same_pi0_parent": r"Same true $\pi^0$ parent",
    "different_pi0_parents": r"Different true $\pi^0$ parents",
    "same_eta_parent": r"Same true $\eta$ parent",
    "same_other_parent": "Same other true parent",
    "other_mixed": "Other/mixed true photons",
    "nonphoton_daughter": r"Matched non-$\gamma$ daughter",
    "unmatched_daughter": "Unmatched daughter",
}


def plot_clasdis_pair_validation(
    arrays: Mapping[str, np.ndarray],
    models: Mapping[str, GradientBoostingClassifier],
    args: argparse.Namespace,
    output_dir: Path,
) -> None:
    s1, s2, valid = score_epgg_photons(arrays, models, args)
    pair_category = clasdis_pair_categories(arrays)
    mgg = np.asarray(arrays["Mh_gammagamma"], dtype=float)

    min_score = np.minimum(s1, s2)
    mean_score = 0.5 * (s1 + s2)

    categories = [
        "same_pi0_parent",
        "different_pi0_parents",
        "same_eta_parent",
        "other_mixed",
        "nonphoton_daughter",
        "unmatched_daughter",
    ]

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.8))

    for ax, values, xlabel in [
        (axes[0], min_score, r"min(BDT score$_1$, BDT score$_2$)"),
        (axes[1], mean_score, r"mean(BDT score$_1$, BDT score$_2$)"),
    ]:
        for category in categories:
            mask = valid & (pair_category == category)

            if np.sum(mask) == 0:
                continue
            #endif

            normalized_hist(
                ax,
                values[mask],
                label=f"{PAIR_LABELS[category]} (N={np.sum(mask):,})",
                bins=60,
                hist_range=(0.0, 1.0),
            )
        #endfor

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7.5)
    #endfor

    fig.suptitle(
        "CLASDIS e'p'gammagammaX: independent pair-level BDT truth validation",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(
        fig,
        output_dir / "09_clasdis_pair_score_validation.png",
    )

    # Reconstructed mass versus conservative pair score.
    mask = valid & np.isfinite(mgg) & np.isfinite(min_score)

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.8))

    ax = axes[0]
    image = ax.hist2d(
        mgg[mask],
        min_score[mask],
        bins=(100, 70),
        range=((0.0, 0.8), (0.0, 1.0)),
    )
    ax.axvline(0.13498, linestyle="--", linewidth=1.0)
    ax.axvline(0.54786, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    ax.set_ylabel(r"min(BDT score$_1$, BDT score$_2$)")
    ax.set_title("All selected CLASDIS pairs")
    fig.colorbar(image[3], ax=ax, label="Pairs")

    ax = axes[1]
    same_pi0 = (
        mask
        & (pair_category == "same_pi0_parent")
    )
    image = ax.hist2d(
        mgg[same_pi0],
        min_score[same_pi0],
        bins=(80, 70),
        range=((0.0, 0.3), (0.0, 1.0)),
    )
    ax.axvline(0.13498, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    ax.set_ylabel(r"min(BDT score$_1$, BDT score$_2$)")
    ax.set_title(
        rf"Same true $\pi^0$ parent (N={np.sum(same_pi0):,})"
    )
    fig.colorbar(image[3], ax=ax, label="Pairs")

    fig.suptitle(
        "CLASDIS reconstructed mass versus pair-level pi0 score",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(
        fig,
        output_dir / "10_clasdis_mgg_vs_pair_score.png",
    )


def plot_mgg_by_pair_score(
    arrays: Mapping[str, np.ndarray],
    models: Mapping[str, GradientBoostingClassifier],
    args: argparse.Namespace,
    output_path: Path,
    title: str,
) -> None:
    s1, s2, valid = score_epgg_photons(arrays, models, args)
    mgg = np.asarray(arrays["Mh_gammagamma"], dtype=float)
    min_score = np.minimum(s1, s2)

    finite = valid & np.isfinite(mgg) & np.isfinite(min_score)

    score_bins = [
        (0.0, 0.2),
        (0.2, 0.5),
        (0.5, 0.8),
        (0.8, 1.000001),
    ]

    fig, ax = plt.subplots(figsize=(10.5, 6.7))

    for lo, hi in score_bins:
        mask = finite & (min_score >= lo) & (min_score < hi)

        if np.sum(mask) == 0:
            continue
        #endif

        normalized_hist(
            ax,
            mgg[mask],
            label=(
                rf"{lo:.1f} $\leq$ min score $<$ {min(hi,1.0):.1f} "
                f"(N={np.sum(mask):,})"
            ),
            bins=100,
            hist_range=(0.0, 0.8),
        )
    #endfor

    ax.axvline(0.13498, linestyle="--", linewidth=1.0)
    ax.axvline(0.54786, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    ax.set_ylabel("Normalized density")
    ax.set_title(title)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, output_path)


def pi0_mass_model(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    sigma: float,
    background0: float,
    background1: float,
) -> np.ndarray:
    """Gaussian pi0 peak plus locally linear combinatorial background."""
    gaussian = amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
    background = background0 + background1 * (x - PI0_MASS_GEV)
    return gaussian + background


def fit_pi0_mass_peak(
    masses: np.ndarray,
) -> Optional[Dict[str, object]]:
    """
    Fit the reconstructed pi0 peak in one cumulative BDT-score selection.

    The reported purity is S/(S+B) from the fitted Gaussian and fitted linear
    background inside +/- PI0_FIT_NSIGMA fitted sigma around the fitted mean.
    """
    values = np.asarray(masses, dtype=float)
    values = values[
        np.isfinite(values)
        & (values >= PI0_FIT_MIN)
        & (values <= PI0_FIT_MAX)
    ]

    if len(values) < 250:
        return None
    #endif

    counts, edges = np.histogram(
        values,
        bins=PI0_FIT_BINS,
        range=(PI0_FIT_MIN, PI0_FIT_MAX),
    )
    centers = 0.5 * (edges[:-1] + edges[1:])

    sideband = (
        (centers < 0.105)
        | (centers > 0.165)
    )
    background_guess = float(
        np.median(counts[sideband])
        if np.any(sideband)
        else np.median(counts)
    )
    amplitude_guess = max(
        float(np.max(counts) - background_guess),
        1.0,
    )

    p0 = [
        amplitude_guess,
        PI0_MASS_GEV,
        0.010,
        max(background_guess, 0.1),
        0.0,
    ]

    lower = [
        0.0,
        0.120,
        0.003,
        0.0,
        -1.0e7,
    ]
    upper = [
        np.inf,
        0.150,
        0.030,
        np.inf,
        1.0e7,
    ]

    sigma_counts = np.sqrt(np.maximum(counts, 1.0))

    try:
        params, covariance = curve_fit(
            pi0_mass_model,
            centers,
            counts,
            p0=p0,
            sigma=sigma_counts,
            absolute_sigma=True,
            bounds=(lower, upper),
            maxfev=50000,
        )
    except Exception:
        return None
    #endtry

    amplitude, mean, sigma, background0, background1 = params

    if not (
        np.isfinite(mean)
        and np.isfinite(sigma)
        and sigma > 0.0
    ):
        return None
    #endif

    peak_window = (
        np.abs(centers - mean)
        <= PI0_FIT_NSIGMA * sigma
    )

    gaussian_counts = amplitude * np.exp(
        -0.5 * ((centers - mean) / sigma) ** 2
    )
    background_counts = (
        background0
        + background1 * (centers - PI0_MASS_GEV)
    )
    background_counts = np.maximum(background_counts, 0.0)

    signal = float(np.sum(gaussian_counts[peak_window]))
    background = float(np.sum(background_counts[peak_window]))
    denominator = signal + background
    purity = signal / denominator if denominator > 0.0 else np.nan

    return {
        "params": params,
        "covariance": covariance,
        "centers": centers,
        "counts": counts,
        "edges": edges,
        "mean": float(mean),
        "sigma": float(sigma),
        "signal": signal,
        "background": background,
        "purity": float(purity),
        "n_fit_range": int(len(values)),
    }


def plot_data_pi0_mass_fits_vs_score(
    arrays: Mapping[str, np.ndarray],
    models: Mapping[str, GradientBoostingClassifier],
    args: argparse.Namespace,
    output_dir: Path,
) -> None:
    """
    Fit the real-data pi0 mass peak for cumulative pair-score requirements.

    Uses the conservative pair discriminator min(score1, score2).
    """
    s1, s2, valid = score_epgg_photons(arrays, models, args)
    mgg = np.asarray(arrays["Mh_gammagamma"], dtype=float)
    min_score = np.minimum(s1, s2)

    finite = valid & np.isfinite(mgg) & np.isfinite(min_score)

    thresholds = np.asarray(PAIR_SCORE_THRESHOLDS, dtype=float)
    fit_results: List[Optional[Dict[str, object]]] = []

    fig, axes = plt.subplots(3, 3, figsize=(15.0, 12.0))

    for ax, threshold in zip(axes.flat, thresholds):
        mask = finite & (min_score >= threshold)
        result = fit_pi0_mass_peak(mgg[mask])
        fit_results.append(result)

        masses = mgg[
            mask
            & (mgg >= PI0_FIT_MIN)
            & (mgg <= PI0_FIT_MAX)
        ]

        if len(masses) > 0:
            counts, edges = np.histogram(
                masses,
                bins=PI0_FIT_BINS,
                range=(PI0_FIT_MIN, PI0_FIT_MAX),
            )
            centers = 0.5 * (edges[:-1] + edges[1:])
            ax.errorbar(
                centers,
                counts,
                yerr=np.sqrt(np.maximum(counts, 1.0)),
                fmt=".",
                markersize=3.5,
                linewidth=0.8,
                label="Data",
            )
        #endif

        if result is not None:
            params = result["params"]
            x_dense = np.linspace(PI0_FIT_MIN, PI0_FIT_MAX, 500)
            total = pi0_mass_model(x_dense, *params)
            background = (
                params[3]
                + params[4] * (x_dense - PI0_MASS_GEV)
            )

            ax.plot(x_dense, total, linewidth=1.5, label="Gaussian + linear bg")
            ax.plot(
                x_dense,
                np.maximum(background, 0.0),
                linestyle="--",
                linewidth=1.2,
                label="Fitted bg",
            )

            mean = float(result["mean"])
            sigma = float(result["sigma"])
            purity = float(result["purity"])
            ax.axvspan(
                mean - PI0_FIT_NSIGMA * sigma,
                mean + PI0_FIT_NSIGMA * sigma,
                alpha=0.10,
            )
            annotation = (
                f"N={np.sum(mask):,}\n"
                rf"$\mu={mean*1000:.1f}$ MeV" + "\n"
                rf"$\sigma={sigma*1000:.1f}$ MeV" + "\n"
                rf"$S/(S+B)={100.0*purity:.1f}\%$"
            )
        else:
            annotation = f"N={np.sum(mask):,}\nfit unavailable"
        #endif

        ax.text(
            0.03,
            0.95,
            annotation,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )
        ax.set_title(rf"min pair score $\geq {threshold:.1f}$")
        ax.set_xlim(PI0_FIT_MIN, PI0_FIT_MAX)
        ax.grid(alpha=0.2)
    #endfor

    for ax in axes[-1, :]:
        ax.set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    #endfor

    for ax in axes[:, 0]:
        ax.set_ylabel("Counts / bin")
    #endfor

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.suptitle(
        (
            r"Real data: fitted $\pi^0\to\gamma\gamma$ mass peak versus "
            r"minimum single-photon BDT score"
        ),
        y=0.992,
    )
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.955),
            ncol=3,
            fontsize=8,
        )
    #endif
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.90])
    save_figure(
        fig,
        output_dir / "13_data_pi0_mass_fits_vs_score.png",
    )

    good = [
        (threshold, result)
        for threshold, result in zip(thresholds, fit_results)
        if result is not None and np.isfinite(result["purity"])
    ]

    if not good:
        return
    #endif

    threshold_values = np.asarray(
        [item[0] for item in good],
        dtype=float,
    )
    purity_values = np.asarray(
        [100.0 * float(item[1]["purity"]) for item in good],
        dtype=float,
    )
    signal_values = np.asarray(
        [float(item[1]["signal"]) for item in good],
        dtype=float,
    )

    if signal_values[0] > 0.0:
        retained_signal = 100.0 * signal_values / signal_values[0]
    else:
        retained_signal = np.full_like(signal_values, np.nan)
    #endif

    fig, ax = plt.subplots(figsize=(8.5, 6.4))
    ax.plot(
        threshold_values,
        purity_values,
        marker="o",
        linewidth=1.8,
        label=(
            rf"Fitted $\pi^0$ fraction in $\pm{PI0_FIT_NSIGMA:.0f}\sigma$ "
            "mass window"
        ),
    )
    ax.set_xlabel(r"Minimum required pair score: min(score$_1$, score$_2$)")
    ax.set_ylabel(r"Fitted $\pi^0$ fraction $S/(S+B)$ [%]")
    ax.set_ylim(0.0, 105.0)
    ax.grid(alpha=0.2)

    ax2 = ax.twinx()
    ax2.plot(
        threshold_values,
        retained_signal,
        marker="s",
        linestyle="--",
        linewidth=1.5,
        label=r"Fitted $\pi^0$ signal retained",
    )
    ax2.set_ylabel(r"Fitted $\pi^0$ signal retained [% of score $\geq0$]")
    ax2.set_ylim(0.0, 105.0)

    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(
        handles1 + handles2,
        labels1 + labels2,
        loc="center right",
        fontsize=8.5,
    )
    ax.set_title(
        (
            r"Real data: $\pi^0$ peak purity and retained signal versus "
            "pair-score requirement"
        )
    )
    fig.tight_layout()
    save_figure(
        fig,
        output_dir / "14_data_pi0_peak_fraction_vs_score.png",
    )




# =============================================================================
# Region workflow
# =============================================================================

def run_region(
    region_name: str,
    detector_value: int,
    categories: Mapping[str, CategorySample],
    data_arrays: Mapping[str, np.ndarray],
    args: argparse.Namespace,
    output_dir: Path,
) -> Tuple[RegionResult, SplitData]:
    progress("=" * 78)
    progress(f"REGION {region_name}: starting")
    progress("=" * 78)

    train, validation, test = build_splits(
        categories,
        seed=args.seed + 100 * detector_value,
    )

    progress(
        f"{region_name}: split sizes — "
        f"train={len(train.X):,}, "
        f"validation={len(validation.X):,}, "
        f"test={len(test.X):,}"
    )

    data_X = selected_data_X(
        data_arrays,
        detector_value,
        args,
    )
    progress(
        f"{region_name}: ordinary real-data selected candidates={len(data_X):,}"
    )

    region_dir = output_dir / region_name
    region_dir.mkdir(parents=True, exist_ok=True)

    category_arrays = {
        name: sample.X
        for name, sample in categories.items()
    }

    plot_region_training_counts(
        categories,
        region_dir / "01_training_category_counts.png",
        region_name,
    )

    plot_feature_distributions(
        category_arrays,
        data_X,
        CORE_PLOT_CATEGORIES,
        region_dir / "02_input_features_core.png",
        (
            f"{region_name}: core physics samples — "
            "independent pi0 sources, genuine gamma backgrounds, and data"
        ),
        include_data=True,
    )

    plot_feature_distributions(
        category_arrays,
        np.empty((0, len(FEATURES))),
        ARTIFACT_PLOT_CATEGORIES,
        region_dir / "03_input_features_reco_artifacts.png",
        (
            f"{region_name}: reconstructed-photon artifact categories — "
            "matched non-photon versus unmatched"
        ),
        include_data=False,
    )

    model = train_model(train, args, region_name)

    overall_auc, core_auc = plot_roc(
        model,
        test,
        region_dir / "04_roc_curves.png",
        region_name,
    )

    plot_core_operating_points(
        model,
        test,
        region_dir / "04b_core_operating_points.png",
        region_name,
    )

    plot_scores_by_category(
        model,
        category_arrays,
        region_dir / "05_bdt_score_by_truth_category.png",
        region_name,
    )

    plot_core_scores_and_data(
        model,
        category_arrays,
        data_X,
        region_dir / "06_bdt_score_core_and_data.png",
        region_name,
    )

    plot_low_mx2ep_score_shift(
        model,
        category_arrays,
        data_X,
        region_dir / "06b_low_mx2ep_score_shift.png",
        region_name,
        args.low_mx2_ep_max,
    )

    plot_feature_importance(
        model,
        region_dir / "07_feature_importance.png",
        region_name,
    )

    balanced = plot_confusion(
        model,
        test,
        region_dir / "08_confusion_matrix_score_0p5.png",
        region_name,
    )

    bundle = {
        "model": model,
        "features": list(FEATURES),
        "feature_labels": dict(FEATURE_LABELS),
        "period": args.period,
        "region": region_name,
        "detector2": detector_value,
        "open_angle_min_deg": args.open_angle_min,
        "mx2_cut_enabled": not args.disable_mx2_cut,
        "mx2_min": args.mx2_min,
        "mx2_max": args.mx2_max,
        "score_interpretation": (
            "Classifier score only; not a calibrated physical probability."
        ),
        "positive_categories": sorted(POSITIVE_CATEGORIES),
        "negative_macro_groups": NEGATIVE_MACRO_GROUPS,
        "negative_macro_targets": NEGATIVE_MACRO_TARGETS,
        "negative_fine_weight_rule": (
            "Population-proportional within each negative macro-category."
        ),
    }
    joblib.dump(bundle, region_dir / "pi0_bdt_model.joblib")
    progress(f"MODEL: {region_dir / 'pi0_bdt_model.joblib'}")

    result = RegionResult(
        name=region_name,
        detector_value=detector_value,
        model=model,
        test_auc=overall_auc,
        core_auc=core_auc,
        balanced_accuracy=balanced,
        category_arrays=category_arrays,
    )

    progress(
        f"REGION {region_name} COMPLETE: "
        f"test AUC={overall_auc:.4f}, "
        f"core AUC={core_auc:.4f}, "
        f"balanced accuracy@0.5={balanced:.4f}"
    )

    return result, test


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()
    paths = build_paths(args)
    ensure_files(paths)

    output_dir = (
        args.output
        / args.period
        / args.sample_tag.replace("/", "_")
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    progress("=" * 90)
    progress("TRUTH-LABELLED PI0 BDT")
    progress(
        f"period={args.period}, sample_tag={args.sample_tag}, "
        f"category_cap={args.max_events_per_category}, seed={args.seed}"
    )
    progress(
        f"common cuts: open_angle_ep2>{args.open_angle_min:.2f} deg; "
        + (
            f"{args.mx2_min:.3f}<Mx2(epgamma)<{args.mx2_max:.3f} GeV^2"
            if not args.disable_mx2_cut
            else "Mx2(epgamma) cut DISABLED"
        )
    )
    progress("=" * 90)

    for name, path in paths.items():
        progress(f"INPUT {name}: {path}")
    #endfor

    # Read each epgamma tree only once, and only the branches needed for
    # training/diagnostics.
    dvcsgen = read_tree_arrays(paths["dvcsgen_epg"], EPG_REQUIRED)
    aaogen = read_tree_arrays(paths["aaogen_epg"], EPG_REQUIRED)
    clasdis = read_tree_arrays(paths["clasdis_epg"], EPG_REQUIRED)
    data = read_tree_arrays(paths["data_epg"], DATA_EPG_REQUIRED)

    arrays_by_source = {
        "dvcsgen": dvcsgen,
        "aaogen": aaogen,
        "clasdis": clasdis,
    }

    categories_by_region = {}

    # Build once for the global count diagnostic. The region workflow repeats
    # only the cheap in-memory masking, not the ROOT reads.
    for region_name, detector_value in REGIONS.items():
        categories_by_region[region_name] = build_region_categories(
            arrays_by_source,
            region_name,
            detector_value,
            args,
        )
    #endfor

    plot_global_truth_counts(
        categories_by_region,
        output_dir / "truth_category_counts.png",
    )

    region_results: Dict[str, RegionResult] = {}
    models: Dict[str, GradientBoostingClassifier] = {}

    for region_name, detector_value in REGIONS.items():
        result, _ = run_region(
            region_name,
            detector_value,
            categories_by_region[region_name],
            data,
            args,
            output_dir,
        )
        region_results[region_name] = result
        models[region_name] = result.model
    #endfor

    # The pair-level validation files are deliberately loaded only after both
    # models have finished training. They never enter model fitting.
    progress("PAIR VALIDATION: loading CLASDIS epgammagammaX")
    clasdis_epgg = read_tree_arrays(
        paths["clasdis_epgg"],
        EPGG_REQUIRED_COMMON + EPGG_TRUTH,
    )

    progress("PAIR VALIDATION: loading data epgammagammaX")
    data_epgg = read_tree_arrays(
        paths["data_epgg"],
        EPGG_REQUIRED_COMMON,
    )

    progress("PAIR VALIDATION: loading AAOgen epgammagammaX")
    aaogen_epgg = read_tree_arrays(
        paths["aaogen_epgg"],
        EPGG_REQUIRED_COMMON,
    )

    plot_clasdis_pair_validation(
        clasdis_epgg,
        models,
        args,
        output_dir,
    )

    plot_mgg_by_pair_score(
        data_epgg,
        models,
        args,
        output_dir / "11_data_mgg_by_pair_score.png",
        (
            "Real data: reconstructed M_gamma_gamma versus "
            "single-photon BDT pair score"
        ),
    )

    plot_mgg_by_pair_score(
        aaogen_epgg,
        models,
        args,
        output_dir / "12_aaogen_mgg_by_pair_score.png",
        (
            "AAOgen: reconstructed M_gamma_gamma versus "
            "single-photon BDT pair score"
        ),
    )

    plot_data_pi0_mass_fits_vs_score(
        data_epgg,
        models,
        args,
        output_dir,
    )

    print("\n" + "=" * 90)
    print("FINAL SUMMARY")
    print("=" * 90)

    for region_name in ["FT", "FD"]:
        result = region_results[region_name]
        print(
            f"{region_name}: "
            f"all-category test AUC={result.test_auc:.4f}, "
            f"pi0-vs-DVCS-gamma AUC={result.core_auc:.4f}, "
            f"balanced accuracy@0.5={result.balanced_accuracy:.4f}"
        )
    #endfor

    print(f"\nPlots/models written under:\n  {output_dir}")
    print("\nMost important plots to inspect first:")
    print("  truth_category_counts.png")
    print("  FT/02_input_features_core.png")
    print("  FD/02_input_features_core.png")
    print("  FT/05_bdt_score_by_truth_category.png")
    print("  FD/05_bdt_score_by_truth_category.png")
    print("  09_clasdis_pair_score_validation.png")
    print("  10_clasdis_mgg_vs_pair_score.png")
    print("  11_data_mgg_by_pair_score.png")
    print("  12_aaogen_mgg_by_pair_score.png")
    print("  13_data_pi0_mass_fits_vs_score.png")
    print("  14_data_pi0_peak_fraction_vs_score.png")
    print("=" * 90)


if __name__ == "__main__":
    main()
