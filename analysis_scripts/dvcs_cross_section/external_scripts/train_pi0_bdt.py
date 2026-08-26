#!/usr/bin/env python3
"""
train_pi0_bdt_intro.py

Introductory BDT study for identifying reconstructed e'p'gammaX candidates
that are pi0-like.

Version-8 diagnostic philosophy
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

Only the topology/exclusivity BDT is trained in this version. Earlier
iterations showed no meaningful performance gain from additionally supplying
Q2, W, x, y, t, and tmin, so the production-kinematics BDT has been removed
to reduce runtime and keep the diagnostic question focused.

This version adds detailed CLASDIS cross-tree matching diagnostics:
    * six individual e/p nearest-neighbor residuals;
    * real-versus-scrambled nearest-neighbor distance;
    * accepted-match multiplicity per reconstructed e'p'pi0X candidate;
    * detector composition of matched e'p'gammaX candidates;
    * matched e'p'gammaX detector versus reconstructed pi0 daughter topology;
    * BDT score versus CLASDIS match distance;
    * CLASDIS BDT score for progressively tighter match-quality selections.

The scrambled control destroys the electron-proton correlation in the e'p'pi0X
sample before repeating the nearest-neighbor search. It therefore provides a
combinatorial reference for how close unrelated events can appear in the same
six-dimensional e/p kinematic space.

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
from sklearn.ensemble import GradientBoostingClassifier
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
    "detector_gamma1",
    "detector_gamma2",
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

# Introductory boosted-decision-tree settings.
#
# We use the classic GradientBoostingClassifier rather than
# HistGradientBoostingClassifier because the latter behaved pathologically on
# the current ifarm scikit-learn build. This configuration uses only 100
# shallow trees and enables sklearn's native verbose fit output so the user can
# watch the boosting iterations proceed.
DEFAULT_BDT_PARAMS = {
    "n_estimators": 100,
    "learning_rate": 0.05,
    "max_depth": 2,
    "min_samples_leaf": 20,
    "subsample": 0.80,
    "random_state": 42,
    "verbose": 1,
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
            "before kinematic matching. Omit for the recommended full CLASDIS "
            "matching diagnostic; partial cross-tree reads can lose true matches."
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
        "--n-estimators",
        type=int,
        default=100,
        help="Number of boosted trees per model. Default: 100.",
    )

    parser.add_argument(
        "--max-depth",
        type=int,
        default=2,
        help="Maximum depth of each individual decision tree. Default: 2.",
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


def _nearest_match_diagnostics(
    left: pd.DataFrame,
    right: pd.DataFrame,
    query_chunk_size: int,
    label: str,
) -> pd.DataFrame:
    """
    For every row in left, find the nearest right-side e/p candidate and return
    the exact physical residuals plus the six-dimensional normalized distance.
    """

    left_embedding = _matching_embedding(left)
    right_embedding = _matching_embedding(right)

    progress(
        f"CLASDIS MATCH [{label}]: fitting nearest-neighbor index to "
        f"{len(right_embedding):,} right-side candidates"
    )
    neighbor_model = NearestNeighbors(n_neighbors=1, algorithm="auto")
    neighbor_model.fit(right_embedding)

    n_left = len(left_embedding)
    n_chunks = int(math.ceil(n_left / query_chunk_size))
    nearest_indices = np.empty(n_left, dtype=np.int64)

    progress(
        f"CLASDIS MATCH [{label}]: querying {n_left:,} epgammaX candidates in "
        f"{n_chunks} chunk(s)"
    )

    for chunk_index in range(n_chunks):
        lo = chunk_index * query_chunk_size
        hi = min((chunk_index + 1) * query_chunk_size, n_left)

        _, chunk_indices = neighbor_model.kneighbors(
            left_embedding[lo:hi],
            return_distance=True,
        )
        nearest_indices[lo:hi] = chunk_indices[:, 0]

        report_every = max(1, n_chunks // 10)
        if (chunk_index + 1) % report_every == 0 or chunk_index + 1 == n_chunks:
            progress(
                f"CLASDIS MATCH [{label}]: nearest-neighbor query "
                f"{hi:,}/{n_left:,} ({100.0 * hi / n_left:.1f}%)"
            )
        #endif
    #endfor

    nearest = right.iloc[nearest_indices].reset_index(drop=True)

    delta_e_p = (
        left["e_p"].to_numpy(dtype=float)
        - nearest["e_p"].to_numpy(dtype=float)
    )
    delta_e_theta = (
        left["e_theta"].to_numpy(dtype=float)
        - nearest["e_theta"].to_numpy(dtype=float)
    )
    delta_e_phi = angle_difference_deg(
        left["e_phi"].to_numpy(dtype=float),
        nearest["e_phi"].to_numpy(dtype=float),
    )
    delta_p_p = (
        left["p1_p"].to_numpy(dtype=float)
        - nearest["p1_p"].to_numpy(dtype=float)
    )
    delta_p_theta = (
        left["p1_theta"].to_numpy(dtype=float)
        - nearest["p1_theta"].to_numpy(dtype=float)
    )
    delta_p_phi = angle_difference_deg(
        left["p1_phi"].to_numpy(dtype=float),
        nearest["p1_phi"].to_numpy(dtype=float),
    )

    components = np.column_stack(
        [
            delta_e_p / MATCH_SCALES["e_p"],
            delta_e_theta / MATCH_SCALES["e_theta"],
            delta_e_phi / MATCH_SCALES["e_phi"],
            delta_p_p / MATCH_SCALES["p1_p"],
            delta_p_theta / MATCH_SCALES["p1_theta"],
            delta_p_phi / MATCH_SCALES["p1_phi"],
        ]
    )
    normalized_distance = np.sqrt(np.sum(np.square(components), axis=1))

    diagnostics = pd.DataFrame(
        {
            "_epg_index": left["_epg_index"].to_numpy(dtype=int),
            "_nearest_eppi0_index": nearest["_eppi0_index"].to_numpy(dtype=int),
            "delta_e_p": delta_e_p,
            "delta_e_theta": delta_e_theta,
            "delta_e_phi": delta_e_phi,
            "delta_p_p": delta_p_p,
            "delta_p_theta": delta_p_theta,
            "delta_p_phi": delta_p_phi,
            "match_distance": normalized_distance,
            "detector2": left["detector2"].to_numpy(dtype=int),
            "p2_p": left["p2_p"].to_numpy(dtype=float),
            "p2_theta": left["p2_theta"].to_numpy(dtype=float),
        }
    )

    if "Mh_gammagamma" in nearest.columns:
        diagnostics["matched_Mh_gammagamma"] = nearest[
            "Mh_gammagamma"
        ].to_numpy(dtype=float)
    #endif

    if "detector_gamma1" in nearest.columns:
        diagnostics["matched_detector_gamma1"] = nearest[
            "detector_gamma1"
        ].to_numpy(dtype=int)
    #endif

    if "detector_gamma2" in nearest.columns:
        diagnostics["matched_detector_gamma2"] = nearest[
            "detector_gamma2"
        ].to_numpy(dtype=int)
    #endif

    return diagnostics


def match_epgamma_to_eppi0_by_kinematics(
    epg: pd.DataFrame,
    eppi0: pd.DataFrame,
    max_normalized_distance: float = 4.0,
    query_chunk_size: int = 100_000,
    seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Construct a pi0-enriched CLASDIS control sample using ONLY reconstructed
    electron/proton kinematics.

    Returns
    -------
    matched:
        Accepted real CLASDIS e'p'gammaX candidates. These retain the original
        e'p'gammaX branches and receive detailed match-quality columns.

    real_diagnostics:
        One row per finite e'p'gammaX candidate, containing its nearest real
        reconstructed e'p'pi0X candidate and the full e/p residual information.

    scrambled_diagnostics:
        The same nearest-neighbor exercise after destroying the electron-proton
        association inside the e'p'pi0X sample. This is a combinatorial control.

    Notes
    -----
    MC runnum/evnum are never used.

    Multiple e'p'gammaX candidates are explicitly allowed to match the same
    e'p'pi0X candidate. That multiplicity is one of the diagnostics produced
    later because one physical e/p event may contribute multiple photon
    candidates.
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

    left = (
        epg.dropna(subset=match_columns)
        .copy()
        .reset_index(names="_epg_index")
    )
    right = (
        eppi0.dropna(subset=match_columns)
        .copy()
        .reset_index(names="_eppi0_index")
    )

    progress(
        f"CLASDIS MATCH: finite e/p kinematics — "
        f"epgammaX={len(left):,}, eppi0X={len(right):,}"
    )

    if len(left) == 0 or len(right) == 0:
        progress("CLASDIS MATCH: no usable rows")
        empty_match = epg.iloc[0:0].copy().reset_index(drop=True)
        empty_diag = pd.DataFrame()
        return empty_match, empty_diag, empty_diag.copy()
    #endif

    progress("CLASDIS MATCH: running REAL nearest-neighbor association")
    real_diagnostics = _nearest_match_diagnostics(
        left=left,
        right=right,
        query_chunk_size=query_chunk_size,
        label="REAL",
    )

    real_diagnostics["accepted"] = (
        real_diagnostics["match_distance"] <= max_normalized_distance
    )

    # Scrambled control:
    # Keep the reconstructed electron rows fixed, but independently permute the
    # proton quantities. This destroys the true e-p event correlation while
    # retaining the one-particle marginal distributions.
    progress(
        "CLASDIS MATCH: constructing SCRAMBLED control by permuting proton "
        "kinematics relative to the electron"
    )
    rng = np.random.default_rng(seed)
    scrambled_right = right.copy()
    proton_permutation = rng.permutation(len(scrambled_right))

    for column in ["p1_p", "p1_theta", "p1_phi"]:
        scrambled_right[column] = (
            right[column].to_numpy()[proton_permutation]
        )
    #endfor

    scrambled_diagnostics = _nearest_match_diagnostics(
        left=left,
        right=scrambled_right,
        query_chunk_size=query_chunk_size,
        label="SCRAMBLED",
    )
    scrambled_diagnostics["accepted"] = (
        scrambled_diagnostics["match_distance"] <= max_normalized_distance
    )

    accepted = real_diagnostics["accepted"].to_numpy(dtype=bool)
    matched_diag = real_diagnostics.loc[accepted].copy().reset_index(drop=True)

    matched_indices = matched_diag["_epg_index"].to_numpy(dtype=int)
    matched = epg.iloc[matched_indices].copy().reset_index(drop=True)

    matched["control_match"] = True
    matched["control_match_distance"] = matched_diag[
        "match_distance"
    ].to_numpy(dtype=float)
    matched["control_nearest_eppi0_index"] = matched_diag[
        "_nearest_eppi0_index"
    ].to_numpy(dtype=int)

    for column in [
        "delta_e_p",
        "delta_e_theta",
        "delta_e_phi",
        "delta_p_p",
        "delta_p_theta",
        "delta_p_phi",
        "matched_Mh_gammagamma",
        "matched_detector_gamma1",
        "matched_detector_gamma2",
    ]:
        if column in matched_diag.columns:
            matched[f"control_{column}"] = matched_diag[column].to_numpy()
        #endif
    #endfor

    n_accepted = len(matched)
    accepted_fraction = 100.0 * n_accepted / len(real_diagnostics)

    n_scrambled_accepted = int(
        np.count_nonzero(scrambled_diagnostics["accepted"].to_numpy())
    )
    scrambled_fraction = (
        100.0 * n_scrambled_accepted / len(scrambled_diagnostics)
    )

    if n_accepted > 0:
        accepted_distances = matched["control_match_distance"].to_numpy()
        q50, q90, q99 = np.quantile(
            accepted_distances,
            [0.50, 0.90, 0.99],
        )
        progress(
            f"CLASDIS MATCH: REAL accepted {n_accepted:,}/"
            f"{len(real_diagnostics):,} ({accepted_fraction:.2f}%); "
            f"distance median={q50:.3f}, 90%={q90:.3f}, 99%={q99:.3f}"
        )
    else:
        progress("CLASDIS MATCH: REAL accepted zero candidates")
    #endif

    progress(
        f"CLASDIS MATCH: SCRAMBLED accepted {n_scrambled_accepted:,}/"
        f"{len(scrambled_diagnostics):,} ({scrambled_fraction:.2f}%) "
        f"at the same d <= {max_normalized_distance:.2f} requirement"
    )

    if n_accepted > 0:
        unique_pi0 = matched[
            "control_nearest_eppi0_index"
        ].nunique()
        progress(
            f"CLASDIS MATCH: accepted epgammaX candidates map to "
            f"{unique_pi0:,} unique reconstructed eppi0X candidates; "
            f"mean multiplicity={n_accepted / max(unique_pi0, 1):.2f}"
        )
    #endif

    progress(f"CLASDIS MATCH: complete ({elapsed_since(stage_start)})")

    return (
        matched,
        real_diagnostics.reset_index(drop=True),
        scrambled_diagnostics.reset_index(drop=True),
    )


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

def make_model(
    seed: int,
    n_estimators: int,
    max_depth: int,
) -> GradientBoostingClassifier:
    params = dict(DEFAULT_BDT_PARAMS)
    params["random_state"] = seed
    params["n_estimators"] = int(n_estimators)
    params["max_depth"] = int(max_depth)
    return GradientBoostingClassifier(**params)


def train_model(
    train: pd.DataFrame,
    validation: pd.DataFrame,
    features: Sequence[str],
    seed: int,
    n_estimators: int,
    max_depth: int,
) -> GradientBoostingClassifier:
    model = make_model(
        seed=seed,
        n_estimators=n_estimators,
        max_depth=max_depth,
    )

    X_train = train[list(features)].to_numpy(dtype=float)
    y_train = train["label"].to_numpy(dtype=int)
    w_train = class_balancing_weights(y_train)

    fit_start = time.perf_counter()
    progress(
        f"BDT FIT: fitting GradientBoostingClassifier on "
        f"{len(train):,} rows x {len(features)} features; "
        f"{n_estimators} shallow trees (max_depth={max_depth})"
    )
    progress(
        "BDT FIT: sklearn iteration progress will print below "
        "(Iter / Train Loss / Remaining Time)"
    )
    model.fit(
        X_train,
        y_train,
        sample_weight=w_train,
    )
    progress(
        f"BDT FIT: complete after {getattr(model, 'n_estimators_', DEFAULT_BDT_PARAMS['n_estimators'])} "
        f"boosting trees ({elapsed_since(fit_start)})"
    )

    return model


def model_scores(
    model: GradientBoostingClassifier,
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

def _safe_hist_range(
    values: np.ndarray,
    central_quantile: float = 0.995,
    symmetric: bool = False,
) -> Tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return (-1.0, 1.0)
    #endif

    if symmetric:
        bound = float(np.quantile(np.abs(values), central_quantile))
        if bound <= 0.0 or not np.isfinite(bound):
            bound = 1.0
        #endif
        return (-bound, bound)
    #endif

    low = float(np.quantile(values, 1.0 - central_quantile))
    high = float(np.quantile(values, central_quantile))
    if low == high:
        high = low + 1.0
    #endif
    return (low, high)


def pi0_detector_topology(
    gamma1: np.ndarray,
    gamma2: np.ndarray,
) -> np.ndarray:
    """
    Encode reconstructed pi0 daughter detector topology using detector IDs:
        0 -> FT
        1 -> FD
    """
    gamma1 = np.asarray(gamma1, dtype=int)
    gamma2 = np.asarray(gamma2, dtype=int)

    labels = np.full(len(gamma1), "other", dtype=object)

    ftft = (gamma1 == 0) & (gamma2 == 0)
    ftfd = (
        ((gamma1 == 0) & (gamma2 == 1))
        | ((gamma1 == 1) & (gamma2 == 0))
    )
    fdfd = (gamma1 == 1) & (gamma2 == 1)

    labels[ftft] = "FT-FT"
    labels[ftfd] = "FT-FD"
    labels[fdfd] = "FD-FD"

    return labels


def plot_clasdis_match_residuals(
    real_diag: pd.DataFrame,
    scrambled_diag: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    residual_specs = [
        ("delta_e_p", r"$\Delta p_e$ (GeV)"),
        ("delta_e_theta", r"$\Delta\theta_e$ (deg)"),
        ("delta_e_phi", r"$\Delta\phi_e$ (deg)"),
        ("delta_p_p", r"$\Delta p_p$ (GeV)"),
        ("delta_p_theta", r"$\Delta\theta_p$ (deg)"),
        ("delta_p_phi", r"$\Delta\phi_p$ (deg)"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.5))
    axes = axes.ravel()

    for ax, (column, xlabel) in zip(axes, residual_specs):
        real_values = real_diag[column].to_numpy(dtype=float)
        scrambled_values = scrambled_diag[column].to_numpy(dtype=float)

        combined = np.concatenate([real_values, scrambled_values])
        low, high = _safe_hist_range(
            combined,
            central_quantile=0.995,
            symmetric=True,
        )

        ax.hist(
            scrambled_values,
            bins=100,
            range=(low, high),
            density=True,
            histtype="step",
            linewidth=1.5,
            label="Scrambled nearest neighbor",
        )
        ax.hist(
            real_values,
            bins=100,
            range=(low, high),
            density=True,
            histtype="step",
            linewidth=1.5,
            label="Real nearest neighbor",
        )

        ax.axvline(0.0, linestyle="--", linewidth=1.0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.20)
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2)
    fig.suptitle(title, y=1.01)
    fig.tight_layout()

    save_figure(fig, output_path)


def plot_clasdis_match_distance(
    real_diag: pd.DataFrame,
    scrambled_diag: pd.DataFrame,
    max_distance: float,
    output_path: Path,
    title: str,
) -> None:
    real = real_diag["match_distance"].to_numpy(dtype=float)
    scrambled = scrambled_diag["match_distance"].to_numpy(dtype=float)

    finite = np.concatenate(
        [
            real[np.isfinite(real)],
            scrambled[np.isfinite(scrambled)],
        ]
    )

    if len(finite) == 0:
        return
    #endif

    x_max = max(
        max_distance * 2.0,
        float(np.quantile(finite, 0.995)),
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))

    for ax in axes:
        ax.hist(
            scrambled,
            bins=120,
            range=(0.0, x_max),
            density=True,
            histtype="step",
            linewidth=1.6,
            label="Scrambled",
        )
        ax.hist(
            real,
            bins=120,
            range=(0.0, x_max),
            density=True,
            histtype="step",
            linewidth=1.6,
            label="Real",
        )
        ax.axvline(
            max_distance,
            linestyle="--",
            linewidth=1.4,
            label=f"Current cut d={max_distance:g}",
        )
        ax.set_xlabel("Six-dimensional normalized match distance d")
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.20)
    #endfor

    axes[1].set_yscale("log")
    axes[0].set_title("Linear scale")
    axes[1].set_title("Log scale")
    axes[0].legend()

    fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_clasdis_match_multiplicity(
    real_diag: pd.DataFrame,
    n_eppi0: int,
    max_distance: float,
    output_path: Path,
    title: str,
) -> None:
    accepted = real_diag.loc[
        real_diag["match_distance"] <= max_distance
    ]

    nearest = accepted["_nearest_eppi0_index"].to_numpy(dtype=int)

    counts = np.bincount(
        nearest,
        minlength=max(n_eppi0, 1),
    )

    max_display = max(10, min(50, int(np.quantile(counts, 0.995)) + 2))
    clipped = np.minimum(counts, max_display)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))

    bins = np.arange(-0.5, max_display + 1.5, 1.0)
    axes[0].hist(clipped, bins=bins)
    axes[0].set_xlabel(
        f"Accepted epgammaX candidates per eppi0X candidate"
        + (f" (last bin >= {max_display})" if np.any(counts > max_display) else "")
    )
    axes[0].set_ylabel("Number of reconstructed eppi0X candidates")
    axes[0].set_yscale("log")
    axes[0].grid(alpha=0.20)

    nonzero = counts[counts > 0]
    if len(nonzero) > 0:
        sorted_counts = np.sort(nonzero)[::-1]
        axes[1].plot(
            np.arange(1, len(sorted_counts) + 1),
            sorted_counts,
            marker=".",
            linestyle="none",
            markersize=3,
        )
        axes[1].set_yscale("log")
    #endif
    axes[1].set_xlabel("Matched eppi0X candidate rank")
    axes[1].set_ylabel("Accepted epgammaX multiplicity")
    axes[1].grid(alpha=0.20)

    fig.suptitle(
        f"{title}\n"
        f"{len(accepted):,} accepted epgammaX candidates map to "
        f"{np.count_nonzero(counts):,}/{n_eppi0:,} eppi0X candidates"
    )
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_clasdis_detector_composition(
    real_diag: pd.DataFrame,
    max_distance: float,
    output_path: Path,
    title: str,
) -> None:
    accepted = real_diag.loc[
        real_diag["match_distance"] <= max_distance
    ].copy()

    if len(accepted) == 0:
        return
    #endif

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.3))

    detector_categories = ["FT (0)", "FD (1)", "Other"]
    detector_counts = [
        int(np.count_nonzero(accepted["detector2"].to_numpy() == 0)),
        int(np.count_nonzero(accepted["detector2"].to_numpy() == 1)),
        int(np.count_nonzero(
            ~np.isin(accepted["detector2"].to_numpy(), [0, 1])
        )),
    ]

    axes[0].bar(detector_categories, detector_counts)
    axes[0].set_ylabel("Accepted matched epgammaX candidates")
    axes[0].set_title("Matched epgammaX photon detector")
    axes[0].grid(axis="y", alpha=0.20)

    if (
        "matched_detector_gamma1" in accepted.columns
        and "matched_detector_gamma2" in accepted.columns
    ):
        pi0_topology = pi0_detector_topology(
            accepted["matched_detector_gamma1"].to_numpy(),
            accepted["matched_detector_gamma2"].to_numpy(),
        )

        row_labels = ["FT epgamma", "FD epgamma"]
        column_labels = ["FT-FT", "FT-FD", "FD-FD", "other"]
        matrix = np.zeros((2, 4), dtype=int)

        for row_index, epg_detector in enumerate([0, 1]):
            row_mask = accepted["detector2"].to_numpy() == epg_detector
            for column_index, topology in enumerate(column_labels):
                matrix[row_index, column_index] = int(
                    np.count_nonzero(row_mask & (pi0_topology == topology))
                )
            #endfor
        #endfor

        image = axes[1].imshow(matrix, aspect="auto")
        axes[1].set_xticks(
            np.arange(len(column_labels)),
            column_labels,
        )
        axes[1].set_yticks(
            np.arange(len(row_labels)),
            row_labels,
        )
        axes[1].set_xlabel("Matched reconstructed pi0 daughter topology")
        axes[1].set_ylabel("epgammaX p2 detector")
        axes[1].set_title("Detector-topology association")

        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                axes[1].text(
                    j,
                    i,
                    f"{matrix[i, j]:,}",
                    ha="center",
                    va="center",
                )
            #endfor
        #endfor

        fig.colorbar(image, ax=axes[1], label="Candidate count")
    else:
        axes[1].axis("off")
    #endif

    fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_control_score_vs_match_distance(
    control: pd.DataFrame,
    control_scores: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    if len(control) == 0 or "control_match_distance" not in control.columns:
        return
    #endif

    distance = control["control_match_distance"].to_numpy(dtype=float)
    score = np.asarray(control_scores, dtype=float)

    fig, ax = plt.subplots(figsize=(8.0, 5.8))

    if len(control) > 5_000:
        hb = ax.hexbin(
            distance,
            score,
            gridsize=55,
            mincnt=1,
            bins="log",
        )
        fig.colorbar(hb, ax=ax, label="log10(candidate count)")
    else:
        ax.scatter(
            distance,
            score,
            s=8,
            alpha=0.30,
        )
    #endif

    ax.set_xlabel("CLASDIS e/p match distance d")
    ax.set_ylabel("BDT pi0 score")
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.20)
    ax.set_title(title)

    save_figure(fig, output_path)


def plot_control_score_by_match_quality(
    control: pd.DataFrame,
    control_scores: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    if len(control) == 0 or "control_match_distance" not in control.columns:
        return
    #endif

    distance = control["control_match_distance"].to_numpy(dtype=float)
    score = np.asarray(control_scores, dtype=float)

    fig, ax = plt.subplots(figsize=(8.0, 5.8))

    thresholds = [1.0, 2.0, 3.0, 4.0]
    for threshold in thresholds:
        mask = distance <= threshold
        if np.count_nonzero(mask) == 0:
            continue
        #endif

        ax.hist(
            score[mask],
            bins=50,
            range=(0.0, 1.0),
            density=True,
            histtype="step",
            linewidth=1.5,
            label=f"d <= {threshold:g} (N={np.count_nonzero(mask):,})",
        )
    #endfor

    ax.set_xlabel("BDT pi0 score")
    ax.set_ylabel("Normalized density")
    ax.set_xlim(0.0, 1.0)
    ax.grid(alpha=0.20)
    ax.legend()
    ax.set_title(title)

    save_figure(fig, output_path)


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
    model: GradientBoostingClassifier,
    test: pd.DataFrame,
    features: Sequence[str],
    output_path: Path,
    title: str,
    seed: int,
) -> None:
    # Use permutation importance on a capped held-out test sample. This asks
    # how much shuffling each variable degrades the classifier's ROC AUC and is
    # more directly useful here than the tree's built-in impurity importance.
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
    model_output = region_output / "topology"
    model_output.mkdir(parents=True, exist_ok=True)

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
            detector_region_mask(clasdis_control, region)
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

    feature_start = time.perf_counter()
    features = TOPOLOGY_FEATURES

    progress(
        f"REGION {region} / topology: starting with "
        f"{len(features)} input features"
    )

    working = clean_feature_rows(combined, features)

    progress(
        f"REGION {region} / topology: finite-feature cleaning kept "
        f"{len(working):,}/{len(combined):,} training-pool rows"
    )

    if control_region is not None:
        control_clean = clean_feature_rows(control_region, features)
    else:
        control_clean = None
    #endif

    progress(f"REGION {region} / topology: splitting train/validation/test")

    train, validation, test = train_validation_test_split(
        working,
        args.seed,
    )

    progress(
        f"REGION {region} / topology: split sizes — "
        f"train={len(train):,}, validation={len(validation):,}, test={len(test):,}"
    )

    progress(
        f"REGION {region} / topology: plotting input feature distributions"
    )

    plot_input_feature_distributions(
        train=train,
        features=features,
        output_path=model_output / "input_feature_distributions.png",
        title=f"{args.period} {region}: topology-BDT training inputs",
    )

    progress(
        f"REGION {region} / topology: input plots complete; starting BDT fit"
    )

    model = train_model(
        train=train,
        validation=validation,
        features=features,
        seed=args.seed,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
    )

    progress(f"REGION {region} / topology: scoring held-out test sample")
    test_scores = model_scores(model, test, features)

    if control_clean is not None:
        control_scores = model_scores(model, control_clean, features)
    else:
        control_scores = None
    #endif

    progress(
        f"REGION {region} / topology: making ROC and score-distribution plots"
    )

    roc_auc = plot_roc_curve(
        test=test,
        test_scores=test_scores,
        output_path=model_output / "roc_curve.png",
        title=f"{args.period} {region}: topology BDT",
    )

    plot_score_distribution(
        test=test,
        test_scores=test_scores,
        control=control_clean,
        control_scores=control_scores,
        output_path=model_output / "bdt_score_distribution.png",
        title=f"{args.period} {region}: topology BDT score",
    )

    progress(
        f"REGION {region} / topology: starting permutation-importance diagnostic"
    )

    plot_feature_importance(
        model=model,
        test=test,
        features=features,
        output_path=model_output / "feature_importance.png",
        title=f"{args.period} {region}: topology feature importance",
        seed=args.seed,
    )

    progress(
        f"REGION {region} / topology: permutation importance complete; "
        f"making remaining diagnostics"
    )

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
        title=f"{args.period} {region}: topology BDT",
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
            title=f"{args.period} {region}: CLASDIS matched control",
        )

        plot_control_score_vs_match_distance(
            control=control_clean,
            control_scores=control_scores,
            output_path=model_output / "clasdis_score_vs_match_distance.png",
            title=(
                f"{args.period} {region}: CLASDIS BDT score vs e/p match quality"
            ),
        )

        plot_control_score_by_match_quality(
            control=control_clean,
            control_scores=control_scores,
            output_path=model_output / "clasdis_score_by_match_quality.png",
            title=(
                f"{args.period} {region}: CLASDIS score under tighter match cuts"
            ),
        )
    #endif

    progress(f"REGION {region} / topology: saving trained model")

    joblib.dump(
        {
            "model": model,
            "features": list(features),
            "period": args.period,
            "region": region,
            "feature_set": "topology",
        },
        model_output / "bdt_model.joblib",
    )

    summary = DatasetSummary(
        period=args.period,
        region=region,
        feature_set="topology",
        n_positive_total=n_pos,
        n_negative_total=n_neg,
        n_train=len(train),
        n_validation=len(validation),
        n_test=len(test),
        n_clasdis_control=0 if control_clean is None else len(control_clean),
        auc_test=roc_auc,
        balanced_accuracy_test_at_050=balanced_accuracy,
    )

    progress(
        f"REGION {region} / topology: complete — "
        f"test AUC={roc_auc:.5f}, "
        f"balanced accuracy@0.5={balanced_accuracy:.5f} "
        f"({elapsed_since(feature_start)})"
    )

    progress(f"REGION {region}: complete ({elapsed_since(region_start)})")
    return [summary]


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
        f"max-control-events={args.max_control_events}, "
        f"n-estimators={args.n_estimators}, max-depth={args.max_depth}, seed={args.seed}"
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
        progress(
            "CLASDIS WARNING: partial cross-tree reads can remove genuine "
            "eppi0X<->epgammaX partners because the two trees need not share "
            "entry ordering. Full CLASDIS reads are preferred for diagnostics."
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
    clasdis_real_diagnostics: Optional[pd.DataFrame] = None
    clasdis_scrambled_diagnostics: Optional[pd.DataFrame] = None
    n_clasdis_eppi0 = 0

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

        n_clasdis_eppi0 = len(clasdis_eppi0)

        (
            clasdis_control,
            clasdis_real_diagnostics,
            clasdis_scrambled_diagnostics,
        ) = match_epgamma_to_eppi0_by_kinematics(
            epg=clasdis_epg,
            eppi0=clasdis_eppi0,
            seed=args.seed,
        )

        progress(
            f"CLASDIS CONTROL: matched sample contains "
            f"{len(clasdis_control):,} e'p'gammaX candidates"
        )

        match_output = period_output / "clasdis_matching"
        match_output.mkdir(parents=True, exist_ok=True)

        progress("CLASDIS DIAGNOSTICS: plotting real-vs-scrambled residuals")
        plot_clasdis_match_residuals(
            real_diag=clasdis_real_diagnostics,
            scrambled_diag=clasdis_scrambled_diagnostics,
            output_path=match_output / "match_residuals_real_vs_scrambled.png",
            title=f"{args.period}: CLASDIS e/p nearest-neighbor residuals",
        )

        progress("CLASDIS DIAGNOSTICS: plotting nearest-neighbor distance")
        plot_clasdis_match_distance(
            real_diag=clasdis_real_diagnostics,
            scrambled_diag=clasdis_scrambled_diagnostics,
            max_distance=4.0,
            output_path=match_output / "match_distance_real_vs_scrambled.png",
            title=f"{args.period}: CLASDIS real vs scrambled matching",
        )

        progress("CLASDIS DIAGNOSTICS: plotting accepted-match multiplicity")
        plot_clasdis_match_multiplicity(
            real_diag=clasdis_real_diagnostics,
            n_eppi0=n_clasdis_eppi0,
            max_distance=4.0,
            output_path=match_output / "match_multiplicity.png",
            title=f"{args.period}: CLASDIS matched-candidate multiplicity",
        )

        progress("CLASDIS DIAGNOSTICS: plotting detector composition")
        plot_clasdis_detector_composition(
            real_diag=clasdis_real_diagnostics,
            max_distance=4.0,
            output_path=match_output / "match_detector_composition.png",
            title=f"{args.period}: CLASDIS matched detector composition",
        )

        progress("CLASDIS DIAGNOSTICS: matching diagnostic plots complete")
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
