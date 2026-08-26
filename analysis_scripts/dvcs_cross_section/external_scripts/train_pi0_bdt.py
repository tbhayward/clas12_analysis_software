#!/usr/bin/env python3
"""
train_pi0_bdt_intro.py

Introductory BDT study for identifying reconstructed e'p'gammaX candidates
that are pi0-like.

Version-5 introductory philosophy
--------------------
Training labels are intentionally conservative:

    positive (y = 1):
        AAOgen e'p'gammaX
        -> the generator is exclusive pi0, so the reconstructed photon sample
           is treated as the clean pi0-origin class.

    negative (y = 0):
        DVCSgen e'p'gammaX
        -> genuine single-photon DVCS/BH events.

CLASDIS is NOT used as a negative training sample because, without generator
ancestry, an unmatched e'p'gammaX event can still be a pi0 event in which the
second photon was not reconstructed.

Instead, CLASDIS e'p'gammaX events that can be matched to reconstructed
CLASDIS e'p'pi0X events are used as an independent pi0-enriched positive
control sample.

The script trains separate models for:
    * FT photons: detector2 == 0
    * FD photons: detector2 == 1

The MC runnum/evnum branches are never used for matching or event identity.
CLASDIS cross-tree matching relies only on reconstructed electron/proton
kinematics.

For each topology it trains two feature configurations:

    topology:
        Uses only reconstructed e'p'gammaX topology/exclusivity variables.

    full:
        Uses the topology variables plus broad production kinematics
        (Q2, W, x, y, t, tmin).

The comparison is useful because a large improvement from "topology" to
"full" can indicate that the classifier is exploiting generator phase-space
differences rather than genuinely learning the one-photon remnant of a pi0.

No reconstructed e'p'pi0X-only variable (for example M_gamma_gamma) is ever
passed to the BDT.

Requirements
------------
    python -m pip install uproot pandas numpy matplotlib scikit-learn joblib

Example
-------
python train_pi0_bdt.py --period fa18_inb --max-events-per-class 200000 --max-control-events 100000

Run once for fa18_inb and once for fa18_out. The standard ROOT-file paths are
hard-coded by period; explicit file arguments remain available as overrides.

Important
---------
The BDT predict_proba output is called a "BDT score" here. It is NOT
interpreted as an absolute physical P(pi0 | event), because the training
sample class priors are artificial and the score has not been calibrated.
"""

from __future__ import annotations

import argparse
import json
import math
import time
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

from sklearn.metrics import (
    auc,
    balanced_accuracy_score,
    confusion_matrix,
    roc_curve,
)
from sklearn.model_selection import train_test_split
from sklearn.neighbors import NearestNeighbors
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance


# ---------------------------------------------------------------------------
# Progress reporting
# ---------------------------------------------------------------------------

PROGRAM_START_TIME = time.perf_counter()


def progress(message: str) -> None:
    """Print a timestamped progress message immediately."""
    elapsed = time.perf_counter() - PROGRAM_START_TIME
    print(f"[+{elapsed:8.1f} s] {message}", flush=True)


def elapsed_since(start_time: float) -> str:
    return f"{time.perf_counter() - start_time:.1f} s"


def dataframe_memory_mb(df: pd.DataFrame) -> float:
    return float(df.memory_usage(index=True, deep=True).sum()) / (1024.0 ** 2)


# ---------------------------------------------------------------------------
# Feature definitions
# ---------------------------------------------------------------------------

# These are quantities available in the e'p'gammaX tree itself.
#
# p2 is the reconstructed photon in e'p'gammaX.
#
# We intentionally do not assign physics names to Mx2/Mx2_1/Mx2_2 here,
# because their precise definitions are analysis-code-specific. The BDT can
# use the branches directly; plot labels retain the branch names so there is
# no possibility of silently mislabeling them.
TOPOLOGY_FEATURES: List[str] = [
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

PRODUCTION_FEATURES: List[str] = [
    "Q2",
    "W",
    "x",
    "y",
    "t",
    "tmin",
]

FEATURE_SETS: Dict[str, List[str]] = {
    "topology": TOPOLOGY_FEATURES,
    "full": TOPOLOGY_FEATURES + PRODUCTION_FEATURES,
}

# Branches needed for detector-region selection, kinematic matching,
# and all possible model inputs. MC runnum/evnum are deliberately not read.
EPG_REQUIRED_BRANCHES: List[str] = sorted(
    set(
        [
            "detector2",
            "e_p",
            "e_theta",
            "e_phi",
            "p1_p",
            "p1_theta",
            "p1_phi",
            "p2_p",
            "p2_theta",
            "p2_phi",
        ]
        + TOPOLOGY_FEATURES
        + PRODUCTION_FEATURES
    )
)

EPPI0_MATCH_BRANCHES: List[str] = [
    "e_p",
    "e_theta",
    "e_phi",
    "p1_p",
    "p1_theta",
    "p1_phi",
    "Mh_gammagamma",
]

# Matching scales are deliberately fairly tight. Cross-tree MC matching uses
# only reconstructed electron/proton kinematics; runnum/evnum are never used.
MATCH_SCALES = {
    "e_p": 0.010,        # GeV
    "e_theta": 0.10,     # degrees
    "e_phi": 0.20,       # degrees
    "p1_p": 0.010,       # GeV
    "p1_theta": 0.10,    # degrees
    "p1_phi": 0.20,      # degrees
}

# Default scikit-learn histogram gradient-boosting settings are intentionally
# modest. HistGradientBoostingClassifier is available with scikit-learn and
# avoids requiring the external xgboost package on ifarm.
DEFAULT_BDT_PARAMS = {
    "max_iter": 300,
    "max_leaf_nodes": 15,
    "max_depth": 3,
    "learning_rate": 0.05,
    "min_samples_leaf": 20,
    "l2_regularization": 1.0,
    "early_stopping": True,
    "validation_fraction": 0.15,
    "n_iter_no_change": 20,
    "random_state": 42,
}

# Period-specific files. Command-line overrides are still accepted, but the
# normal invocation only needs --period.
PERIOD_FILES = {
    "fa18_inb": {
        "aaogen_epg": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
        "dvcsgen_epg": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root",
        "clasdis_epg": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/clasdis/rec_clasdis_rga_fa18_inb_epgammaX.root",
        "clasdis_eppi0": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/clasdis/rec_clasdis_rga_fa18_inb_eppi0X.root",
    },
    "fa18_out": {
        "aaogen_epg": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_out_epgamma_0.40GeV.root",
        "dvcsgen_epg": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_fa18_out_epgamma_0.40GeV.root",
        "clasdis_epg": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/clasdis/rec_clasdis_rga_fa18_out_epgammaX.root",
        "clasdis_eppi0": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/clasdis/rec_clasdis_rga_fa18_out_eppi0X.root",
    },
}


@dataclass
class DatasetSummary:
    period: str
    region: str
    feature_set: str
    n_positive_total: int
    n_negative_total: int
    n_train: int
    n_validation: int
    n_test: int
    n_clasdis_control: int
    auc_test: float
    balanced_accuracy_test_at_050: float


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Train introductory FT/FD BDTs to identify pi0-like e'p'gammaX "
            "candidates using AAOgen positives and DVCSgen negatives."
        )
    )

    parser.add_argument(
        "--period",
        required=True,
        choices=["fa18_inb", "fa18_out"],
        help="Run period label used for output organization.",
    )

    parser.add_argument(
        "--aaogen-epg",
        default=None,
        help="Optional AAOgen e'p'gammaX ROOT-file override.",
    )
    parser.add_argument(
        "--dvcsgen-epg",
        default=None,
        help="Optional DVCSgen e'p'gammaX ROOT-file override.",
    )

    parser.add_argument(
        "--clasdis-epg",
        default=None,
        help=(
            "Optional CLASDIS e'p'gammaX ROOT file. Used only for an independent "
            "pi0-enriched control sample."
        ),
    )
    parser.add_argument(
        "--clasdis-eppi0",
        default=None,
        help=(
            "Optional CLASDIS e'p'pi0X ROOT file. Must be supplied together with "
            "--clasdis-epg to construct the control sample."
        ),
    )

    parser.add_argument(
        "--tree",
        default=None,
        help=(
            "ROOT TTree name. If omitted, the script automatically uses the first "
            "TTree found in each file."
        ),
    )

    parser.add_argument(
        "--output-dir",
        default="output/pi0_bdt_intro",
        help="Top-level output directory.",
    )

    parser.add_argument(
        "--max-events-per-class",
        type=int,
        default=None,
        help=(
            "Optional per-region training cap. For quick tests this also limits "
            "the AAOgen/DVCSgen ROOT read to 2x this value per file, so the full "
            "trees are not needlessly materialized before capping."
        ),
    )

    parser.add_argument(
        "--max-clasdis-read",
        type=int,
        default=None,
        help=(
            "Optional maximum number of ROOT entries read from EACH CLASDIS tree "
            "before kinematic matching. Recommended for quick pipeline tests. "
            "Omit for the eventual full CLASDIS control-sample study."
        ),
    )

    parser.add_argument(
        "--max-control-events",
        type=int,
        default=None,
        help="Optional random cap on the matched CLASDIS control sample.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed.",
    )

    parser.add_argument(
        "--skip-full",
        action="store_true",
        help="Train only the topology feature set.",
    )

    parser.add_argument(
        "--skip-clasdis-control",
        action="store_true",
        help="Do not construct the optional CLASDIS matched control sample.",
    )

    return parser.parse_args()


def first_tree_name(filename: str) -> str:
    with uproot.open(filename) as root_file:
        for key, obj in root_file.items():
            if isinstance(obj, uproot.behaviors.TTree.TTree):
                return key.split(";")[0]
            #endif
        #endfor
    raise RuntimeError(f"No TTree found in {filename}")


def resolve_tree_name(filename: str, requested_tree: Optional[str]) -> str:
    if requested_tree is not None:
        return requested_tree
    #endif
    return first_tree_name(filename)


def validate_branches(
    filename: str,
    tree_name: str,
    requested_branches: Sequence[str],
) -> List[str]:
    with uproot.open(filename) as root_file:
        tree = root_file[tree_name]
        available = set(tree.keys())

    missing = [branch for branch in requested_branches if branch not in available]
    if missing:
        raise RuntimeError(
            f"\nFile: {filename}\n"
            f"Tree: {tree_name}\n"
            f"Missing required branches:\n  " + "\n  ".join(missing)
        )
    #endif

    return list(requested_branches)


def load_dataframe(
    filename: str,
    requested_tree: Optional[str],
    branches: Sequence[str],
    source: str,
    entry_limit: Optional[int] = None,
) -> pd.DataFrame:
    stage_start = time.perf_counter()
    progress(f"LOAD {source}: opening {filename}")

    tree_name = resolve_tree_name(filename, requested_tree)
    progress(f"LOAD {source}: using tree '{tree_name}'; validating {len(branches)} branches")
    validate_branches(filename, tree_name, branches)

    with uproot.open(filename) as root_file:
        tree = root_file[tree_name]
        n_total = int(tree.num_entries)
        if entry_limit is None:
            n_read = n_total
        else:
            n_read = min(int(entry_limit), n_total)
        #endif

        progress(
            f"LOAD {source}: tree contains {n_total:,} rows; "
            f"reading {n_read:,} ({100.0*n_read/max(n_total,1):.2f}%)"
        )
        progress(
            f"LOAD {source}: reading {len(branches)} branches into memory "
            f"(entry_start=0, entry_stop={n_read:,})"
        )

        df = tree.arrays(
            list(branches),
            entry_start=0,
            entry_stop=n_read,
            library="pd",
        )
    #endwith

    df = df.reset_index(drop=True)
    df["source"] = source

    progress(
        f"LOAD {source}: complete — {len(df):,} rows, "
        f"~{dataframe_memory_mb(df):.1f} MB in memory ({elapsed_since(stage_start)})"
    )
    return df


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def wrap_delta_phi_deg(delta_phi: np.ndarray) -> np.ndarray:
    return (delta_phi + 180.0) % 360.0 - 180.0


def clean_feature_rows(
    df: pd.DataFrame,
    features: Sequence[str],
) -> pd.DataFrame:
    out = df.copy()

    # Replace +/-inf with NaN. scikit-learn gradient boosting can technically handle NaN, but for this
    # first diagnostic study we require all chosen inputs to be finite so that
    # feature plots and class comparisons are straightforward.
    out[list(features)] = out[list(features)].replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=list(features))

    return out.reset_index(drop=True)


def random_cap(
    df: pd.DataFrame,
    max_rows: Optional[int],
    seed: int,
) -> pd.DataFrame:
    if max_rows is None or len(df) <= max_rows:
        return df.reset_index(drop=True)
    #endif

    return df.sample(n=max_rows, random_state=seed).reset_index(drop=True)


def detector_region_mask(
    df: pd.DataFrame,
    region: str,
) -> pd.Series:
    """Use the reconstructed photon detector assignment directly."""
    if region == "FT":
        return df["detector2"] == 0
    elif region == "FD":
        return df["detector2"] == 1
    else:
        raise ValueError(f"Unknown detector region: {region}")
    #endif


def build_training_sample(
    aaogen: pd.DataFrame,
    dvcsgen: pd.DataFrame,
    region: str,
    max_events_per_class: Optional[int],
    seed: int,
) -> pd.DataFrame:
    pos = aaogen.loc[detector_region_mask(aaogen, region)].copy()
    neg = dvcsgen.loc[detector_region_mask(dvcsgen, region)].copy()

    pos["label"] = 1
    neg["label"] = 0

    pos = random_cap(pos, max_events_per_class, seed)
    neg = random_cap(neg, max_events_per_class, seed + 1)

    return pd.concat([pos, neg], ignore_index=True)


def angle_difference_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0


def _matching_embedding(df: pd.DataFrame) -> np.ndarray:
    """
    Build a KD/nearest-neighbor-friendly embedding of the electron/proton
    kinematics. Azimuths are represented by sin/cos components so the
    +/-180-degree boundary is handled correctly.

    For small angular differences, the chord-distance scaling below is
    approximately the wrapped delta-phi divided by the requested phi scale.
    """
    e_phi = np.deg2rad(df["e_phi"].to_numpy(dtype=float))
    p_phi = np.deg2rad(df["p1_phi"].to_numpy(dtype=float))

    e_phi_factor = (180.0 / np.pi) / MATCH_SCALES["e_phi"]
    p_phi_factor = (180.0 / np.pi) / MATCH_SCALES["p1_phi"]

    return np.column_stack(
        [
            df["e_p"].to_numpy(dtype=float) / MATCH_SCALES["e_p"],
            df["e_theta"].to_numpy(dtype=float) / MATCH_SCALES["e_theta"],
            np.cos(e_phi) * e_phi_factor,
            np.sin(e_phi) * e_phi_factor,
            df["p1_p"].to_numpy(dtype=float) / MATCH_SCALES["p1_p"],
            df["p1_theta"].to_numpy(dtype=float) / MATCH_SCALES["p1_theta"],
            np.cos(p_phi) * p_phi_factor,
            np.sin(p_phi) * p_phi_factor,
        ]
    )


def match_epgamma_to_eppi0_by_kinematics(
    epg: pd.DataFrame,
    eppi0: pd.DataFrame,
    max_normalized_distance: float = 4.0,
    query_chunk_size: int = 100_000,
) -> pd.DataFrame:
    """
    Construct a pi0-enriched CLASDIS control sample using ONLY reconstructed
    electron/proton kinematics.

    MC runnum and evnum are not used anywhere in this matching.

    The nearest-neighbor lookup is performed in chunks so that long jobs emit
    useful progress messages rather than appearing to hang. The nearest-neighbor
    stage only proposes a candidate; the final acceptance uses the exact
    six-dimensional normalized distance with wrapped azimuthal differences.
    """

    stage_start = time.perf_counter()
    match_columns = [
        "e_p",
        "e_theta",
        "e_phi",
        "p1_p",
        "p1_theta",
        "p1_phi",
    ]

    progress(
        f"CLASDIS MATCH: preparing kinematic samples — "
        f"epgammaX={len(epg):,}, eppi0X={len(eppi0):,}"
    )

    left = epg.dropna(subset=match_columns).copy().reset_index(names="_epg_index")
    right = eppi0.dropna(subset=match_columns).copy().reset_index(drop=True)

    progress(
        f"CLASDIS MATCH: finite e/p kinematics — "
        f"epgammaX={len(left):,}, eppi0X={len(right):,}"
    )

    if len(left) == 0 or len(right) == 0:
        progress("CLASDIS MATCH: no usable rows; returning empty control sample")
        return epg.iloc[0:0].copy().reset_index(drop=True)
    #endif

    progress("CLASDIS MATCH: building scaled e/p kinematic embeddings")
    left_embedding = _matching_embedding(left)
    right_embedding = _matching_embedding(right)

    progress(
        f"CLASDIS MATCH: fitting nearest-neighbor index to "
        f"{len(right_embedding):,} reconstructed eppi0X candidates"
    )
    nn_start = time.perf_counter()
    neighbor_model = NearestNeighbors(n_neighbors=1, algorithm="auto")
    neighbor_model.fit(right_embedding)
    progress(f"CLASDIS MATCH: nearest-neighbor index ready ({elapsed_since(nn_start)})")

    n_left = len(left_embedding)
    n_chunks = int(math.ceil(n_left / query_chunk_size))
    nearest_indices = np.empty(n_left, dtype=np.int64)

    progress(
        f"CLASDIS MATCH: querying {n_left:,} epgammaX candidates in "
        f"{n_chunks} chunk(s) of <= {query_chunk_size:,}"
    )
    query_start = time.perf_counter()

    for chunk_index in range(n_chunks):
        lo = chunk_index * query_chunk_size
        hi = min((chunk_index + 1) * query_chunk_size, n_left)
        _, chunk_indices = neighbor_model.kneighbors(
            left_embedding[lo:hi],
            return_distance=True,
        )
        nearest_indices[lo:hi] = chunk_indices[:, 0]

        # Report every chunk for small jobs; approximately every 10% for large jobs.
        report_every = max(1, n_chunks // 10)
        if (chunk_index + 1) % report_every == 0 or chunk_index + 1 == n_chunks:
            frac = 100.0 * hi / n_left
            progress(
                f"CLASDIS MATCH: nearest-neighbor query {hi:,}/{n_left:,} "
                f"({frac:.1f}%)"
            )
        #endif
    #endfor

    progress(f"CLASDIS MATCH: neighbor queries complete ({elapsed_since(query_start)})")
    progress("CLASDIS MATCH: applying exact wrapped-angle six-dimensional distance")

    nearest = right.iloc[nearest_indices].reset_index(drop=True)

    components = np.column_stack(
        [
            (left["e_p"].to_numpy() - nearest["e_p"].to_numpy())
            / MATCH_SCALES["e_p"],
            (left["e_theta"].to_numpy() - nearest["e_theta"].to_numpy())
            / MATCH_SCALES["e_theta"],
            angle_difference_deg(
                left["e_phi"].to_numpy(),
                nearest["e_phi"].to_numpy(),
            )
            / MATCH_SCALES["e_phi"],
            (left["p1_p"].to_numpy() - nearest["p1_p"].to_numpy())
            / MATCH_SCALES["p1_p"],
            (left["p1_theta"].to_numpy() - nearest["p1_theta"].to_numpy())
            / MATCH_SCALES["p1_theta"],
            angle_difference_deg(
                left["p1_phi"].to_numpy(),
                nearest["p1_phi"].to_numpy(),
            )
            / MATCH_SCALES["p1_phi"],
        ]
    )

    normalized_distance = np.sqrt(np.sum(np.square(components), axis=1))
    accepted = normalized_distance <= max_normalized_distance

    matched_indices = left.loc[accepted, "_epg_index"].to_numpy(dtype=int)
    matched = epg.iloc[matched_indices].copy()
    matched["control_match"] = True
    matched["control_match_distance"] = normalized_distance[accepted]

    n_accepted = int(np.count_nonzero(accepted))
    accepted_fraction = 100.0 * n_accepted / len(left)
    if n_accepted > 0:
        accepted_distances = normalized_distance[accepted]
        q50, q90, q99 = np.quantile(accepted_distances, [0.50, 0.90, 0.99])
        progress(
            f"CLASDIS MATCH: accepted {n_accepted:,}/{len(left):,} "
            f"({accepted_fraction:.2f}%); distance median={q50:.3f}, "
            f"90%={q90:.3f}, 99%={q99:.3f}"
        )
    else:
        progress(
            f"CLASDIS MATCH: accepted 0/{len(left):,} candidates at "
            f"distance <= {max_normalized_distance}"
        )
    #endif

    progress(f"CLASDIS MATCH: complete ({elapsed_since(stage_start)})")
    return matched.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Train/validation/test splitting and weighting
# ---------------------------------------------------------------------------

def train_validation_test_split(
    df: pd.DataFrame,
    seed: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Approximately 70% train, 15% validation, 15% test with class stratification.

    MC runnum/evnum are deliberately not used for grouping or identity. This
    introductory version therefore performs a row-level split. If later study
    shows that a single generated event can contribute multiple highly
    correlated e'p'gammaX rows, the production version should construct a
    kinematic event-group identifier without using runnum/evnum.
    """

    train, remainder = train_test_split(
        df,
        train_size=0.70,
        random_state=seed,
        stratify=df["label"],
        shuffle=True,
    )

    validation, test = train_test_split(
        remainder,
        train_size=0.50,
        random_state=seed + 1,
        stratify=remainder["label"],
        shuffle=True,
    )

    return (
        train.reset_index(drop=True),
        validation.reset_index(drop=True),
        test.reset_index(drop=True),
    )


def class_balancing_weights(labels: np.ndarray) -> np.ndarray:
    """
    Give the positive and negative classes equal total training weight.

    This prevents the arbitrary generated event counts from becoming an
    implicit prior P(pi0).
    """

    labels = np.asarray(labels, dtype=int)

    n_pos = np.count_nonzero(labels == 1)
    n_neg = np.count_nonzero(labels == 0)

    if n_pos == 0 or n_neg == 0:
        raise RuntimeError(
            f"Cannot construct class weights: n_pos={n_pos}, n_neg={n_neg}"
        )
    #endif

    weights = np.ones(len(labels), dtype=float)
    weights[labels == 1] = 0.5 / n_pos
    weights[labels == 0] = 0.5 / n_neg

    # Rescale to mean weight ~1 for numerical convenience.
    weights *= len(labels)

    return weights


# ---------------------------------------------------------------------------
# Model training
# ---------------------------------------------------------------------------

def make_model(seed: int) -> HistGradientBoostingClassifier:
    params = dict(DEFAULT_BDT_PARAMS)
    params["random_state"] = seed
    return HistGradientBoostingClassifier(**params)


def train_model(
    train: pd.DataFrame,
    validation: pd.DataFrame,
    features: Sequence[str],
    seed: int,
) -> HistGradientBoostingClassifier:
    model = make_model(seed)

    X_train = train[list(features)].to_numpy(dtype=float)
    y_train = train["label"].to_numpy(dtype=int)
    w_train = class_balancing_weights(y_train)

    # HistGradientBoostingClassifier performs its own internal stratified
    # validation for early stopping. The separately held-out validation sample
    # remains available as an independent diagnostic dataset.
    fit_start = time.perf_counter()
    progress(
        f"BDT FIT: fitting HistGradientBoostingClassifier on "
        f"{len(train):,} rows x {len(features)} features"
    )
    model.fit(
        X_train,
        y_train,
        sample_weight=w_train,
    )
    progress(
        f"BDT FIT: complete after {getattr(model, 'n_iter_', '?')} boosting iterations "
        f"({elapsed_since(fit_start)})"
    )

    return model


def model_scores(
    model: HistGradientBoostingClassifier,
    df: pd.DataFrame,
    features: Sequence[str],
) -> np.ndarray:
    if len(df) == 0:
        return np.asarray([], dtype=float)
    #endif

    X = df[list(features)].to_numpy(dtype=float)
    return model.predict_proba(X)[:, 1]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def finite_common_range(
    values_a: np.ndarray,
    values_b: np.ndarray,
    lower_quantile: float = 0.005,
    upper_quantile: float = 0.995,
) -> Tuple[float, float]:
    values = np.concatenate([values_a, values_b])
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return 0.0, 1.0
    #endif

    low = float(np.quantile(values, lower_quantile))
    high = float(np.quantile(values, upper_quantile))

    if not np.isfinite(low) or not np.isfinite(high) or low == high:
        low = float(np.min(values))
        high = float(np.max(values))
    #endif

    if low == high:
        high = low + 1.0
    #endif

    return low, high


def plot_input_feature_distributions(
    train: pd.DataFrame,
    features: Sequence[str],
    output_path: Path,
    title: str,
) -> None:
    ncols = 3
    nrows = int(math.ceil(len(features) / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(5.0 * ncols, 3.8 * nrows),
    )
    axes = np.atleast_1d(axes).ravel()

    pos = train.loc[train["label"] == 1]
    neg = train.loc[train["label"] == 0]

    for ax, feature in zip(axes, features):
        a = neg[feature].to_numpy(dtype=float)
        b = pos[feature].to_numpy(dtype=float)

        low, high = finite_common_range(a, b)

        ax.hist(
            a,
            bins=60,
            range=(low, high),
            density=True,
            histtype="step",
            linewidth=1.5,
            label="DVCSgen: y=0",
        )
        ax.hist(
            b,
            bins=60,
            range=(low, high),
            density=True,
            histtype="step",
            linewidth=1.5,
            label="AAOgen: y=1",
        )

        ax.set_xlabel(feature)
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.20)
    #endfor

    for ax in axes[len(features):]:
        ax.axis("off")
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2)
    fig.suptitle(title, y=1.01)
    fig.tight_layout()

    save_figure(fig, output_path)


def plot_score_distribution(
    test: pd.DataFrame,
    test_scores: np.ndarray,
    control: Optional[pd.DataFrame],
    control_scores: Optional[np.ndarray],
    output_path: Path,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.5))

    y = test["label"].to_numpy(dtype=int)

    ax.hist(
        test_scores[y == 0],
        bins=50,
        range=(0.0, 1.0),
        density=True,
        histtype="step",
        linewidth=1.8,
        label="DVCSgen test: genuine gamma",
    )
    ax.hist(
        test_scores[y == 1],
        bins=50,
        range=(0.0, 1.0),
        density=True,
        histtype="step",
        linewidth=1.8,
        label="AAOgen test: pi0",
    )

    if (
        control is not None
        and control_scores is not None
        and len(control_scores) > 0
    ):
        ax.hist(
            control_scores,
            bins=50,
            range=(0.0, 1.0),
            density=True,
            histtype="step",
            linewidth=1.8,
            linestyle="--",
            label="CLASDIS matched eppi0 control",
        )
    #endif

    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("BDT pi0 score")
    ax.set_ylabel("Normalized density")
    ax.set_title(title)
    ax.grid(alpha=0.20)
    ax.legend()

    save_figure(fig, output_path)


def plot_roc_curve(
    test: pd.DataFrame,
    test_scores: np.ndarray,
    output_path: Path,
    title: str,
) -> float:
    y = test["label"].to_numpy(dtype=int)
    fpr, tpr, _ = roc_curve(y, test_scores)
    roc_auc = float(auc(fpr, tpr))

    fig, ax = plt.subplots(figsize=(6.5, 6.0))

    ax.plot(fpr, tpr, linewidth=2.0, label=f"BDT AUC = {roc_auc:.4f}")
    ax.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", label="Random classifier")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("False-positive rate: genuine gamma tagged as pi0")
    ax.set_ylabel("True-positive rate: pi0 tagged as pi0")
    ax.set_title(title)
    ax.grid(alpha=0.20)
    ax.legend(loc="lower right")

    save_figure(fig, output_path)
    return roc_auc


def plot_feature_importance(
    model: HistGradientBoostingClassifier,
    test: pd.DataFrame,
    features: Sequence[str],
    output_path: Path,
    title: str,
    seed: int,
) -> None:
    # HistGradientBoostingClassifier does not expose impurity-based feature
    # importances. Use permutation importance on a capped held-out test sample
    # instead. This is also a more directly interpretable diagnostic: how much
    # does shuffling each variable degrade the classifier's ROC AUC?
    importance_df = test.copy()
    if len(importance_df) > 20000:
        importance_df = importance_df.sample(n=20000, random_state=seed)
    #endif

    X = importance_df[list(features)].to_numpy(dtype=float)
    y = importance_df["label"].to_numpy(dtype=int)

    progress(
        f"PERMUTATION IMPORTANCE: starting on {len(importance_df):,} held-out rows, "
        f"{len(features)} features x 3 repeats"
    )
    importance_start = time.perf_counter()
    result = permutation_importance(
        model,
        X,
        y,
        scoring="roc_auc",
        n_repeats=3,
        random_state=seed,
        n_jobs=1,
    )
    progress(f"PERMUTATION IMPORTANCE: complete ({elapsed_since(importance_start)})")

    importance = np.asarray(result.importances_mean, dtype=float)
    uncertainty = np.asarray(result.importances_std, dtype=float)
    order = np.argsort(importance)

    fig, ax = plt.subplots(figsize=(8.0, max(5.0, 0.42 * len(features))))

    ax.barh(
        np.asarray(features)[order],
        importance[order],
        xerr=uncertainty[order],
    )

    ax.set_xlabel("Permutation importance: decrease in ROC AUC")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.20)

    save_figure(fig, output_path)


def plot_score_vs_selected_features(
    test: pd.DataFrame,
    test_scores: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    selected = [
        feature
        for feature in ["p2_p", "p2_theta", "Mx2", "Mx2_1", "Mx2_2", "pTmiss"]
        if feature in test.columns
    ]

    ncols = 3
    nrows = int(math.ceil(len(selected) / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(5.0 * ncols, 4.0 * nrows),
    )
    axes = np.atleast_1d(axes).ravel()

    y = test["label"].to_numpy(dtype=int)

    for ax, feature in zip(axes, selected):
        values = test[feature].to_numpy(dtype=float)

        # Plot a deterministic subset if the test sample is huge.
        if len(values) > 100_000:
            rng = np.random.default_rng(12345)
            idx = rng.choice(len(values), size=100_000, replace=False)
        else:
            idx = np.arange(len(values))
        #endif

        neg_idx = idx[y[idx] == 0]
        pos_idx = idx[y[idx] == 1]

        ax.scatter(
            values[neg_idx],
            test_scores[neg_idx],
            s=3,
            alpha=0.12,
            label="DVCSgen",
        )
        ax.scatter(
            values[pos_idx],
            test_scores[pos_idx],
            s=3,
            alpha=0.12,
            label="AAOgen",
        )

        ax.set_xlabel(feature)
        ax.set_ylabel("BDT pi0 score")
        ax.set_ylim(0.0, 1.0)
        ax.grid(alpha=0.20)
    #endfor

    for ax in axes[len(selected):]:
        ax.axis("off")
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2)
    fig.suptitle(title, y=1.01)
    fig.tight_layout()

    save_figure(fig, output_path)


def plot_confusion_matrix_at_half(
    test: pd.DataFrame,
    test_scores: np.ndarray,
    output_path: Path,
    title: str,
) -> float:
    y_true = test["label"].to_numpy(dtype=int)
    y_pred = (test_scores >= 0.5).astype(int)

    cm = confusion_matrix(y_true, y_pred, labels=[0, 1], normalize="true")
    bal_acc = float(balanced_accuracy_score(y_true, y_pred))

    fig, ax = plt.subplots(figsize=(6.2, 5.5))

    image = ax.imshow(cm, vmin=0.0, vmax=1.0)

    for i in range(2):
        for j in range(2):
            ax.text(
                j,
                i,
                f"{100.0 * cm[i, j]:.1f}%",
                ha="center",
                va="center",
                fontsize=13,
            )
        #endfor
    #endfor

    ax.set_xticks([0, 1], ["tag gamma", "tag pi0"])
    ax.set_yticks([0, 1], ["true gamma", "true pi0"])
    ax.set_xlabel("Prediction at score = 0.5")
    ax.set_ylabel("Truth label")
    ax.set_title(f"{title}\nBalanced accuracy = {bal_acc:.4f}")
    fig.colorbar(image, ax=ax, label="Row-normalized fraction")

    save_figure(fig, output_path)
    return bal_acc


def plot_control_score_vs_kinematics(
    control: pd.DataFrame,
    control_scores: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    if len(control) == 0:
        return
    #endif

    selected = ["p2_p", "p2_theta", "x", "Q2", "t", "pTmiss"]
    selected = [feature for feature in selected if feature in control.columns]

    ncols = 3
    nrows = int(math.ceil(len(selected) / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(5.0 * ncols, 4.0 * nrows),
    )
    axes = np.atleast_1d(axes).ravel()

    if len(control) > 100_000:
        rng = np.random.default_rng(67890)
        idx = rng.choice(len(control), size=100_000, replace=False)
    else:
        idx = np.arange(len(control))
    #endif

    for ax, feature in zip(axes, selected):
        ax.scatter(
            control[feature].to_numpy(dtype=float)[idx],
            control_scores[idx],
            s=3,
            alpha=0.15,
        )
        ax.set_xlabel(feature)
        ax.set_ylabel("BDT pi0 score")
        ax.set_ylim(0.0, 1.0)
        ax.grid(alpha=0.20)
    #endfor

    for ax in axes[len(selected):]:
        ax.axis("off")
    #endfor

    fig.suptitle(title)
    fig.tight_layout()

    save_figure(fig, output_path)


def plot_topology_vs_full_roc(
    roc_data: Dict[str, Tuple[np.ndarray, np.ndarray, float]],
    output_path: Path,
    title: str,
) -> None:
    if len(roc_data) < 2:
        return
    #endif

    fig, ax = plt.subplots(figsize=(6.5, 6.0))

    for feature_set_name, (fpr, tpr, roc_auc) in roc_data.items():
        ax.plot(
            fpr,
            tpr,
            linewidth=2.0,
            label=f"{feature_set_name}: AUC = {roc_auc:.4f}",
        )
    #endfor

    ax.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", label="Random classifier")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_title(title)
    ax.grid(alpha=0.20)
    ax.legend(loc="lower right")

    save_figure(fig, output_path)


# ---------------------------------------------------------------------------
# Main per-region training
# ---------------------------------------------------------------------------

def run_region(
    args: argparse.Namespace,
    aaogen: pd.DataFrame,
    dvcsgen: pd.DataFrame,
    clasdis_control: Optional[pd.DataFrame],
    region: str,
    period_output: Path,
) -> List[DatasetSummary]:

    region_start = time.perf_counter()
    progress(f"REGION {region}: starting")
    region_output = period_output / region.lower()
    region_output.mkdir(parents=True, exist_ok=True)

    raw_pos = int(np.count_nonzero(detector_region_mask(aaogen, region)))
    raw_neg = int(np.count_nonzero(detector_region_mask(dvcsgen, region)))
    progress(
        f"REGION {region}: raw detector-selected rows — "
        f"AAOgen={raw_pos:,}, DVCSgen={raw_neg:,}"
    )

    combined = build_training_sample(
        aaogen=aaogen,
        dvcsgen=dvcsgen,
        region=region,
        max_events_per_class=args.max_events_per_class,
        seed=args.seed,
    )

    n_pos = int(np.count_nonzero(combined["label"].to_numpy() == 1))
    n_neg = int(np.count_nonzero(combined["label"].to_numpy() == 0))

    progress(
        f"REGION {region}: training pool after optional cap — "
        f"{n_pos:,} AAOgen positives, {n_neg:,} DVCSgen negatives"
    )

    if n_pos < 100 or n_neg < 100:
        warnings.warn(
            f"{args.period} {region}: very small training sample "
            f"(positive={n_pos}, negative={n_neg})."
        )
    #endif

    control_region = None

    if clasdis_control is not None:
        control_region = clasdis_control.loc[
            detector_region_mask(
                clasdis_control,
                region,
            )
        ].copy()

        control_region = random_cap(
            control_region,
            args.max_control_events,
            args.seed + 10,
        )

        progress(
            f"REGION {region}: CLASDIS matched control after detector selection/cap — "
            f"{len(control_region):,} candidates"
        )
    #endif

    feature_set_names = ["topology"]
    if not args.skip_full:
        feature_set_names.append("full")
    #endif

    summaries: List[DatasetSummary] = []
    roc_comparison: Dict[str, Tuple[np.ndarray, np.ndarray, float]] = {}

    for feature_set_name in feature_set_names:
        feature_start = time.perf_counter()
        features = FEATURE_SETS[feature_set_name]
        progress(
            f"REGION {region} / {feature_set_name}: starting with "
            f"{len(features)} input features"
        )

        model_output = region_output / feature_set_name
        model_output.mkdir(parents=True, exist_ok=True)

        working = clean_feature_rows(combined, features)
        progress(
            f"REGION {region} / {feature_set_name}: finite-feature cleaning kept "
            f"{len(working):,}/{len(combined):,} training-pool rows"
        )

        # Clean control rows with the same feature requirements so the model is
        # always evaluated on a compatible input matrix.
        if control_region is not None:
            control_clean = clean_feature_rows(control_region, features)
        else:
            control_clean = None
        #endif

        progress(f"REGION {region} / {feature_set_name}: splitting train/validation/test")
        train, validation, test = train_validation_test_split(
            working,
            args.seed,
        )

        progress(
            f"REGION {region} / {feature_set_name}: split sizes — "
            f"train={len(train):,}, validation={len(validation):,}, test={len(test):,}"
        )

        progress(f"REGION {region} / {feature_set_name}: plotting input feature distributions")
        plot_input_feature_distributions(
            train=train,
            features=features,
            output_path=model_output / "input_feature_distributions.png",
            title=(
                f"{args.period} {region}: training inputs "
                f"({feature_set_name} feature set)"
            ),
        )

        progress(f"REGION {region} / {feature_set_name}: input plots complete; starting BDT fit")
        model = train_model(
            train=train,
            validation=validation,
            features=features,
            seed=args.seed,
        )

        progress(f"REGION {region} / {feature_set_name}: scoring held-out test sample")
        test_scores = model_scores(model, test, features)

        if control_clean is not None:
            control_scores = model_scores(model, control_clean, features)
        else:
            control_scores = None
        #endif

        progress(f"REGION {region} / {feature_set_name}: making ROC and score-distribution plots")
        roc_auc = plot_roc_curve(
            test=test,
            test_scores=test_scores,
            output_path=model_output / "roc_curve.png",
            title=f"{args.period} {region}: {feature_set_name} BDT",
        )

        y_test = test["label"].to_numpy(dtype=int)
        fpr, tpr, _ = roc_curve(y_test, test_scores)
        roc_comparison[feature_set_name] = (fpr, tpr, roc_auc)

        plot_score_distribution(
            test=test,
            test_scores=test_scores,
            control=control_clean,
            control_scores=control_scores,
            output_path=model_output / "bdt_score_distribution.png",
            title=f"{args.period} {region}: {feature_set_name} BDT score",
        )

        progress(f"REGION {region} / {feature_set_name}: starting permutation-importance diagnostic")
        plot_feature_importance(
            model=model,
            test=test,
            features=features,
            output_path=model_output / "feature_importance.png",
            title=f"{args.period} {region}: {feature_set_name} feature importance",
            seed=args.seed,
        )

        progress(f"REGION {region} / {feature_set_name}: permutation importance complete; making remaining diagnostics")
        plot_score_vs_selected_features(
            test=test,
            test_scores=test_scores,
            output_path=model_output / "score_vs_selected_features.png",
            title=f"{args.period} {region}: score behavior on held-out MC",
        )

        balanced_accuracy = plot_confusion_matrix_at_half(
            test=test,
            test_scores=test_scores,
            output_path=model_output / "confusion_matrix_score_0p5.png",
            title=f"{args.period} {region}: {feature_set_name} BDT",
        )

        if (
            control_clean is not None
            and control_scores is not None
            and len(control_clean) > 0
        ):
            plot_control_score_vs_kinematics(
                control=control_clean,
                control_scores=control_scores,
                output_path=model_output / "clasdis_control_score_vs_kinematics.png",
                title=(
                    f"{args.period} {region}: CLASDIS matched-pi0 control "
                    f"({feature_set_name})"
                ),
            )
        #endif

        # Persist only the compact model/metadata needed to reproduce the
        # classifier. The scientific output of this introductory script is
        # intentionally plot-focused rather than CSV-focused.
        progress(f"REGION {region} / {feature_set_name}: saving trained model")
        joblib.dump(
            {
                "model": model,
                "features": list(features),
                "period": args.period,
                "region": region,
                "feature_set": feature_set_name,
            },
            model_output / "bdt_model.joblib",
        )

        summary = DatasetSummary(
            period=args.period,
            region=region,
            feature_set=feature_set_name,
            n_positive_total=n_pos,
            n_negative_total=n_neg,
            n_train=len(train),
            n_validation=len(validation),
            n_test=len(test),
            n_clasdis_control=0 if control_clean is None else len(control_clean),
            auc_test=roc_auc,
            balanced_accuracy_test_at_050=balanced_accuracy,
        )
        summaries.append(summary)

        progress(
            f"REGION {region} / {feature_set_name}: complete — "
            f"test AUC={roc_auc:.5f}, balanced accuracy@0.5={balanced_accuracy:.5f} "
            f"({elapsed_since(feature_start)})"
        )
    #endfor

    progress(f"REGION {region}: making topology-vs-full ROC comparison")
    plot_topology_vs_full_roc(
        roc_data=roc_comparison,
        output_path=region_output / "roc_topology_vs_full.png",
        title=f"{args.period} {region}: topology vs full feature sets",
    )

    progress(f"REGION {region}: complete ({elapsed_since(region_start)})")
    return summaries


# ---------------------------------------------------------------------------
# Program entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    main_start = time.perf_counter()

    np.random.seed(args.seed)
    progress(
        f"START: period={args.period}, max-events-per-class={args.max_events_per_class}, "
        f"max-clasdis-read={args.max_clasdis_read}, "
        f"max-control-events={args.max_control_events}, seed={args.seed}"
    )

    period_output = Path(args.output_dir) / args.period
    period_output.mkdir(parents=True, exist_ok=True)

    defaults = PERIOD_FILES[args.period]
    aaogen_epg_file = args.aaogen_epg or defaults["aaogen_epg"]
    dvcsgen_epg_file = args.dvcsgen_epg or defaults["dvcsgen_epg"]
    clasdis_epg_file = args.clasdis_epg or defaults["clasdis_epg"]
    clasdis_eppi0_file = args.clasdis_eppi0 or defaults["clasdis_eppi0"]

    progress(f"INPUT AAOgen e'p'gammaX: {aaogen_epg_file}")
    progress(f"INPUT DVCSgen e'p'gammaX: {dvcsgen_epg_file}")
    if not args.skip_clasdis_control:
        progress(f"INPUT CLASDIS e'p'gammaX: {clasdis_epg_file}")
        progress(f"INPUT CLASDIS e'p'pi0X: {clasdis_eppi0_file}")
    #endif

    # IMPORTANT: for quick tests, constrain ROOT I/O itself rather than
    # reading the complete generator trees and only capping afterward.
    if args.max_events_per_class is None:
        training_read_limit = None
    else:
        training_read_limit = 2 * int(args.max_events_per_class)
    #endif

    if training_read_limit is None:
        progress("I/O MODE: full AAOgen/DVCSgen ROOT trees will be read")
    else:
        progress(
            f"I/O MODE: quick training read enabled — at most "
            f"{training_read_limit:,} rows from each AAOgen/DVCSgen file"
        )
    #endif

    if args.skip_clasdis_control:
        progress("I/O MODE: CLASDIS control disabled")
    elif args.max_clasdis_read is None:
        progress("I/O MODE: full CLASDIS trees will be read for the control sample")
    else:
        progress(
            f"I/O MODE: quick CLASDIS read enabled — at most "
            f"{args.max_clasdis_read:,} rows from each CLASDIS tree"
        )
    #endif

    aaogen = load_dataframe(
        filename=aaogen_epg_file,
        requested_tree=args.tree,
        branches=EPG_REQUIRED_BRANCHES,
        source="aaogen",
        entry_limit=training_read_limit,
    )

    dvcsgen = load_dataframe(
        filename=dvcsgen_epg_file,
        requested_tree=args.tree,
        branches=EPG_REQUIRED_BRANCHES,
        source="dvcsgen",
        entry_limit=training_read_limit,
    )

    clasdis_control: Optional[pd.DataFrame] = None

    have_clasdis_pair = not args.skip_clasdis_control

    if have_clasdis_pair:
        clasdis_epg = load_dataframe(
            filename=clasdis_epg_file,
            requested_tree=args.tree,
            branches=EPG_REQUIRED_BRANCHES,
            source="clasdis_epg",
            entry_limit=args.max_clasdis_read,
        )

        clasdis_eppi0 = load_dataframe(
            filename=clasdis_eppi0_file,
            requested_tree=args.tree,
            branches=EPPI0_MATCH_BRANCHES,
            source="clasdis_eppi0",
            entry_limit=args.max_clasdis_read,
        )

        progress(
            "CLASDIS CONTROL: starting cross-tree e/p kinematic matching "
            "(no MC runnum/evnum usage)"
        )

        clasdis_control = match_epgamma_to_eppi0_by_kinematics(
            epg=clasdis_epg,
            eppi0=clasdis_eppi0,
        )

        progress(
            f"CLASDIS CONTROL: matched sample contains "
            f"{len(clasdis_control):,} e'p'gammaX candidates"
        )
    #endif

    all_summaries: List[DatasetSummary] = []

    for region in ["FT", "FD"]:
        region_summaries = run_region(
            args=args,
            aaogen=aaogen,
            dvcsgen=dvcsgen,
            clasdis_control=clasdis_control,
            region=region,
            period_output=period_output,
        )
        all_summaries.extend(region_summaries)
    #endfor

    # A tiny JSON summary is kept only so the numerical headline results can
    # be inspected without reading plots by eye. The workflow remains
    # deliberately plot-centered.
    summary_path = period_output / "training_summary.json"

    with summary_path.open("w") as f:
        json.dump(
            [asdict(summary) for summary in all_summaries],
            f,
            indent=2,
        )
    #endwith

    progress(f"DONE: complete analysis runtime {elapsed_since(main_start)}")
    print("")
    print(f"Plots/models: {period_output}", flush=True)
    print(f"Numerical summary: {summary_path}", flush=True)


if __name__ == "__main__":
    main()
