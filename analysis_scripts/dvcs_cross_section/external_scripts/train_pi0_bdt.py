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
    * DIRECT reconstructed-pi0 topology from detector_gamma1/detector_gamma2;
    * detector composition of all and matched e'p'gammaX candidates;
    * matched e'p'gammaX detector versus reconstructed pi0 daughter topology;
    * multiplicity separated by direct FT-FT / FT-FD / FD-FD pi0 topology;
    * exact-two-photon M(gamma gamma), pi0 momentum, and detector-pair closure;
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

# Diagnostic feature configurations. The nominal topology BDT remains the
# production/reference classifier. The additional models are trained only to
# understand where its discrimination comes from.
ABLATION_FEATURE_SETS: Dict[str, List[str]] = {
    "nominal": list(TOPOLOGY_FEATURES),

    # Leave-one-important-variable-out tests.
    "no_p2_p": [
        feature for feature in TOPOLOGY_FEATURES
        if feature != "p2_p"
    ],
    "no_theta_gamma_gamma": [
        feature for feature in TOPOLOGY_FEATURES
        if feature != "theta_gamma_gamma"
    ],
    "no_Emiss2": [
        feature for feature in TOPOLOGY_FEATURES
        if feature != "Emiss2"
    ],
    "no_pTmiss": [
        feature for feature in TOPOLOGY_FEATURES
        if feature != "pTmiss"
    ],

    # Remove all three stored missing-mass-like branches together.
    "no_Mx2_group": [
        feature for feature in TOPOLOGY_FEATURES
        if feature not in {"Mx2", "Mx2_1", "Mx2_2"}
    ],

    # Compact model containing the variables that have repeatedly appeared
    # among the strongest diagnostics.
    "minimal4": [
        "p2_p",
        "theta_gamma_gamma",
        "Emiss2",
        "pTmiss",
    ],
}

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
    "p2_p",
    "p2_theta",
    "p2_phi",
    "Mh_gammagamma",
    "detector_gamma1",
    "detector_gamma2",
]

# CLASDIS parent-event matching configuration.
#
# This intentionally follows the established matching used in the photon-
# efficiency study: convert reconstructed electron/proton momenta to Cartesian
# (px, py, pz), build a six-dimensional nearest-neighbor index, and require
# both a small scaled Euclidean distance and a tight per-component momentum
# agreement.
MATCH_COMPONENT_TOL_GEV = 0.002
MATCH_MAX_SCALED_DISTANCE = 2.0


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
        "--skip-ablations",
        action="store_true",
        help=(
            "Train only the nominal topology BDT and skip the diagnostic "
            "feature-ablation models."
        ),
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


def infer_angle_unit(values: np.ndarray, name: str) -> str:
    """
    Infer whether an angular branch is stored in radians or degrees.

    CLAS12 analysis trees used here are normally in radians, but this helper
    keeps the matcher robust and mirrors the explicit angle handling used in
    the photon-efficiency analysis.

    Heuristic:
        values comfortably within [-2*pi, 2*pi] -> radians
        otherwise                                -> degrees
    """
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]

    if len(finite) == 0:
        raise RuntimeError(f"Cannot infer angle unit for empty/non-finite branch {name}")
    #endif

    q995 = float(np.quantile(np.abs(finite), 0.995))

    if q995 <= (2.0 * np.pi + 0.25):
        return "rad"
    #endif

    return "deg"


def angles_to_radians(
    theta: np.ndarray,
    phi: np.ndarray,
    theta_name: str,
    phi_name: str,
) -> Tuple[np.ndarray, np.ndarray]:
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)

    theta_unit = infer_angle_unit(theta, theta_name)
    phi_unit = infer_angle_unit(phi, phi_name)

    if theta_unit == "deg":
        theta = np.deg2rad(theta)
    #endif

    if phi_unit == "deg":
        phi = np.deg2rad(phi)
    #endif

    return theta, phi


def spherical_to_cartesian(
    p: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    theta_name: str,
    phi_name: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert momentum magnitude and polar/azimuthal angles to Cartesian
    momentum components. Angles are auto-detected as radians or degrees.
    """
    p = np.asarray(p, dtype=float)
    theta_rad, phi_rad = angles_to_radians(
        theta,
        phi,
        theta_name=theta_name,
        phi_name=phi_name,
    )

    sin_theta = np.sin(theta_rad)

    px = p * sin_theta * np.cos(phi_rad)
    py = p * sin_theta * np.sin(phi_rad)
    pz = p * np.cos(theta_rad)

    return px, py, pz


def build_ep_cartesian_components(
    df: pd.DataFrame,
) -> Dict[str, np.ndarray]:
    e_px, e_py, e_pz = spherical_to_cartesian(
        p=df["e_p"].to_numpy(dtype=float),
        theta=df["e_theta"].to_numpy(dtype=float),
        phi=df["e_phi"].to_numpy(dtype=float),
        theta_name="e_theta",
        phi_name="e_phi",
    )

    p_px, p_py, p_pz = spherical_to_cartesian(
        p=df["p1_p"].to_numpy(dtype=float),
        theta=df["p1_theta"].to_numpy(dtype=float),
        phi=df["p1_phi"].to_numpy(dtype=float),
        theta_name="p1_theta",
        phi_name="p1_phi",
    )

    return {
        "e_px": e_px,
        "e_py": e_py,
        "e_pz": e_pz,
        "p_px": p_px,
        "p_py": p_py,
        "p_pz": p_pz,
    }


def _matching_embedding(df: pd.DataFrame) -> np.ndarray:
    components = build_ep_cartesian_components(df)

    matrix = np.column_stack(
        [
            components["e_px"],
            components["e_py"],
            components["e_pz"],
            components["p_px"],
            components["p_py"],
            components["p_pz"],
        ]
    )

    return matrix / MATCH_COMPONENT_TOL_GEV


def _nearest_match_diagnostics(
    left: pd.DataFrame,
    right: pd.DataFrame,
    query_chunk_size: int,
    label: str,
) -> pd.DataFrame:
    """
    For every e'p'gammaX candidate in left, find the nearest e'p'pi0X parent
    candidate in six-dimensional Cartesian electron/proton momentum space.

    Distances are measured in units of MATCH_COMPONENT_TOL_GEV, i.e. each
    Cartesian momentum component is scaled by 0.002 GeV before constructing
    the Euclidean nearest-neighbor distance.
    """

    left_embedding = _matching_embedding(left)
    right_embedding = _matching_embedding(right)

    progress(
        f"CLASDIS MATCH [{label}]: fitting nearest-neighbor index to "
        f"{len(right_embedding):,} right-side candidates in Cartesian e/p momentum space"
    )

    neighbor_model = NearestNeighbors(
        n_neighbors=1,
        algorithm="auto",
    )
    neighbor_model.fit(right_embedding)

    n_left = len(left_embedding)
    n_chunks = int(math.ceil(n_left / query_chunk_size))
    nearest_indices = np.empty(n_left, dtype=np.int64)
    nearest_scaled_distances = np.empty(n_left, dtype=float)

    progress(
        f"CLASDIS MATCH [{label}]: querying {n_left:,} epgammaX candidates in "
        f"{n_chunks} chunk(s)"
    )

    for chunk_index in range(n_chunks):
        lo = chunk_index * query_chunk_size
        hi = min((chunk_index + 1) * query_chunk_size, n_left)

        chunk_distances, chunk_indices = neighbor_model.kneighbors(
            left_embedding[lo:hi],
            return_distance=True,
        )

        nearest_indices[lo:hi] = chunk_indices[:, 0]
        nearest_scaled_distances[lo:hi] = chunk_distances[:, 0]

        report_every = max(1, n_chunks // 10)
        if (chunk_index + 1) % report_every == 0 or chunk_index + 1 == n_chunks:
            progress(
                f"CLASDIS MATCH [{label}]: nearest-neighbor query "
                f"{hi:,}/{n_left:,} ({100.0 * hi / n_left:.1f}%)"
            )
        #endif
    #endfor

    nearest = right.iloc[nearest_indices].reset_index(drop=True)

    left_components = build_ep_cartesian_components(left)
    right_components = build_ep_cartesian_components(nearest)

    residual_columns = {}
    for component in [
        "e_px",
        "e_py",
        "e_pz",
        "p_px",
        "p_py",
        "p_pz",
    ]:
        residual_columns[f"delta_{component}"] = (
            left_components[component]
            - right_components[component]
        )
    #endfor

    residual_matrix = np.column_stack(
        [
            residual_columns["delta_e_px"],
            residual_columns["delta_e_py"],
            residual_columns["delta_e_pz"],
            residual_columns["delta_p_px"],
            residual_columns["delta_p_py"],
            residual_columns["delta_p_pz"],
        ]
    )

    max_abs_component_delta = np.max(
        np.abs(residual_matrix),
        axis=1,
    )

    # Recompute the exact scaled Euclidean distance from the Cartesian
    # residuals. This should agree with the nearest-neighbor distance returned
    # above up to floating-point precision.
    exact_scaled_distance = np.sqrt(
        np.sum(
            np.square(
                residual_matrix / MATCH_COMPONENT_TOL_GEV
            ),
            axis=1,
        )
    )

    diagnostics_dict = {
        "_epg_index": left["_epg_index"].to_numpy(dtype=int),
        "_nearest_eppi0_index": nearest["_eppi0_index"].to_numpy(dtype=int),
        "match_distance": exact_scaled_distance,
        "neighbor_distance": nearest_scaled_distances,
        "max_abs_component_delta": max_abs_component_delta,
        "detector2": left["detector2"].to_numpy(dtype=int),
        "p2_p": left["p2_p"].to_numpy(dtype=float),
        "p2_theta": left["p2_theta"].to_numpy(dtype=float),
    }
    diagnostics_dict.update(residual_columns)

    diagnostics = pd.DataFrame(diagnostics_dict)

    for column in ["p2_p", "p2_theta", "p2_phi", "Mh_gammagamma"]:
        if column in nearest.columns:
            diagnostics[f"matched_{column}"] = nearest[
                column
            ].to_numpy(dtype=float)
        #endif
    #endfor

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
    max_scaled_distance: float = MATCH_MAX_SCALED_DISTANCE,
    component_tolerance_gev: float = MATCH_COMPONENT_TOL_GEV,
    query_chunk_size: int = 100_000,
    seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Construct a pi0-enriched CLASDIS control sample using the established
    photon-efficiency parent-event matching prescription.

    Matching uses ONLY reconstructed electron/proton kinematics:
        electron: (px, py, pz)
        proton:   (px, py, pz)

    The spherical branches (p, theta, phi) are converted to Cartesian momentum
    components with explicit radian/degree handling.

    A real match is accepted only if BOTH:
        1. scaled six-dimensional Euclidean distance <= max_scaled_distance
        2. every individual Cartesian component residual satisfies
           |Delta p_component| <= component_tolerance_gev

    Default values follow the established photon-efficiency matcher:
        component_tolerance_gev = 0.002 GeV
        max_scaled_distance     = 2.0

    MC runnum/evnum are never used.

    Returns
    -------
    matched:
        Accepted CLASDIS e'p'gammaX candidates with detailed matching columns.

    real_diagnostics:
        One row per finite e'p'gammaX candidate with its nearest real
        reconstructed e'p'pi0X parent candidate.

    scrambled_diagnostics:
        The same nearest-neighbor exercise after independently permuting the
        proton Cartesian kinematics relative to the electron in the e'p'pi0X
        sample. This destroys the true e-p parent association while preserving
        one-particle marginal distributions.
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
        f"CLASDIS MATCH: preparing established Cartesian parent matcher — "
        f"epgammaX={len(epg):,}, eppi0X={len(eppi0):,}"
    )
    progress(
        f"CLASDIS MATCH: criteria — max scaled distance <= "
        f"{max_scaled_distance:.3f}; every Cartesian component |Delta p| <= "
        f"{component_tolerance_gev:.4f} GeV"
    )

    left = (
        epg.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=match_columns)
        .copy()
        .reset_index(names="_epg_index")
    )
    right = (
        eppi0.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=match_columns)
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

    # Report inferred angular units once for transparency.
    for branch in ["e_theta", "e_phi", "p1_theta", "p1_phi"]:
        unit = infer_angle_unit(
            pd.concat([left[branch], right[branch]], ignore_index=True).to_numpy(),
            branch,
        )
        progress(f"CLASDIS MATCH: inferred {branch} unit = {unit}")
    #endfor

    progress(
        "CLASDIS MATCH: running REAL Cartesian nearest-neighbor association"
    )

    real_diagnostics = _nearest_match_diagnostics(
        left=left,
        right=right,
        query_chunk_size=query_chunk_size,
        label="REAL",
    )

    real_distance_ok = (
        real_diagnostics["match_distance"].to_numpy(dtype=float)
        <= max_scaled_distance
    )
    real_component_ok = (
        real_diagnostics["max_abs_component_delta"].to_numpy(dtype=float)
        <= component_tolerance_gev
    )

    real_diagnostics["distance_ok"] = real_distance_ok
    real_diagnostics["component_ok"] = real_component_ok
    real_diagnostics["accepted"] = real_distance_ok & real_component_ok

    # Scrambled control.
    #
    # We scramble proton rows relative to electron rows BEFORE Cartesian
    # matching. That removes any true e-p parent correlation but preserves the
    # marginal electron and proton distributions.
    progress(
        "CLASDIS MATCH: constructing SCRAMBLED control by permuting proton "
        "kinematics relative to electron kinematics"
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

    scrambled_distance_ok = (
        scrambled_diagnostics["match_distance"].to_numpy(dtype=float)
        <= max_scaled_distance
    )
    scrambled_component_ok = (
        scrambled_diagnostics["max_abs_component_delta"].to_numpy(dtype=float)
        <= component_tolerance_gev
    )

    scrambled_diagnostics["distance_ok"] = scrambled_distance_ok
    scrambled_diagnostics["component_ok"] = scrambled_component_ok
    scrambled_diagnostics["accepted"] = (
        scrambled_distance_ok & scrambled_component_ok
    )

    accepted = real_diagnostics["accepted"].to_numpy(dtype=bool)
    matched_diag = (
        real_diagnostics.loc[accepted]
        .copy()
        .reset_index(drop=True)
    )

    matched_indices = matched_diag["_epg_index"].to_numpy(dtype=int)
    matched = (
        epg.iloc[matched_indices]
        .copy()
        .reset_index(drop=True)
    )

    matched["control_match"] = True
    matched["control_match_distance"] = matched_diag[
        "match_distance"
    ].to_numpy(dtype=float)
    matched["control_max_abs_component_delta"] = matched_diag[
        "max_abs_component_delta"
    ].to_numpy(dtype=float)
    matched["control_nearest_eppi0_index"] = matched_diag[
        "_nearest_eppi0_index"
    ].to_numpy(dtype=int)

    for column in [
        "delta_e_px",
        "delta_e_py",
        "delta_e_pz",
        "delta_p_px",
        "delta_p_py",
        "delta_p_pz",
        "matched_p2_p",
        "matched_p2_theta",
        "matched_p2_phi",
        "matched_Mh_gammagamma",
        "matched_detector_gamma1",
        "matched_detector_gamma2",
    ]:
        if column in matched_diag.columns:
            matched[f"control_{column}"] = (
                matched_diag[column].to_numpy()
            )
        #endif
    #endfor

    n_real_total = len(real_diagnostics)
    n_real_distance = int(np.count_nonzero(real_distance_ok))
    n_real_component = int(np.count_nonzero(real_component_ok))
    n_real_accepted = int(np.count_nonzero(accepted))

    n_scrambled_total = len(scrambled_diagnostics)
    n_scrambled_accepted = int(
        np.count_nonzero(
            scrambled_diagnostics["accepted"].to_numpy(dtype=bool)
        )
    )

    progress(
        f"CLASDIS MATCH: REAL distance criterion alone accepts "
        f"{n_real_distance:,}/{n_real_total:,} "
        f"({100.0*n_real_distance/max(n_real_total,1):.3f}%)"
    )
    progress(
        f"CLASDIS MATCH: REAL component criterion alone accepts "
        f"{n_real_component:,}/{n_real_total:,} "
        f"({100.0*n_real_component/max(n_real_total,1):.3f}%)"
    )
    progress(
        f"CLASDIS MATCH: REAL final accepted "
        f"{n_real_accepted:,}/{n_real_total:,} "
        f"({100.0*n_real_accepted/max(n_real_total,1):.3f}%)"
    )
    progress(
        f"CLASDIS MATCH: SCRAMBLED final accepted "
        f"{n_scrambled_accepted:,}/{n_scrambled_total:,} "
        f"({100.0*n_scrambled_accepted/max(n_scrambled_total,1):.5f}%)"
    )

    if n_real_accepted > 0:
        accepted_distances = matched[
            "control_match_distance"
        ].to_numpy(dtype=float)
        accepted_components = matched[
            "control_max_abs_component_delta"
        ].to_numpy(dtype=float)

        qd50, qd90, qd99 = np.quantile(
            accepted_distances,
            [0.50, 0.90, 0.99],
        )
        qc50, qc90, qc99 = np.quantile(
            accepted_components,
            [0.50, 0.90, 0.99],
        )

        progress(
            f"CLASDIS MATCH: accepted scaled-distance quantiles — "
            f"median={qd50:.4f}, 90%={qd90:.4f}, 99%={qd99:.4f}"
        )
        progress(
            f"CLASDIS MATCH: accepted max-component |Delta p| quantiles (GeV) — "
            f"median={qc50:.6f}, 90%={qc90:.6f}, 99%={qc99:.6f}"
        )

        unique_pi0 = matched[
            "control_nearest_eppi0_index"
        ].nunique()

        progress(
            f"CLASDIS MATCH: accepted epgammaX candidates map to "
            f"{unique_pi0:,}/{len(right):,} unique reconstructed eppi0X candidates; "
            f"mean multiplicity={n_real_accepted/max(unique_pi0,1):.3f}"
        )

        n_ft = int(np.count_nonzero(matched["detector2"].to_numpy() == 0))
        n_fd = int(np.count_nonzero(matched["detector2"].to_numpy() == 1))
        n_other = n_real_accepted - n_ft - n_fd

        progress(
            f"CLASDIS MATCH: accepted epgammaX detector composition — "
            f"FT={n_ft:,}, FD={n_fd:,}, other={n_other:,}"
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
        ("delta_e_px", r"$\Delta p_{x,e}$ (GeV)"),
        ("delta_e_py", r"$\Delta p_{y,e}$ (GeV)"),
        ("delta_e_pz", r"$\Delta p_{z,e}$ (GeV)"),
        ("delta_p_px", r"$\Delta p_{x,p}$ (GeV)"),
        ("delta_p_py", r"$\Delta p_{y,p}$ (GeV)"),
        ("delta_p_pz", r"$\Delta p_{z,p}$ (GeV)"),
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

        # Ensure the physically relevant +/- few-MeV region is always visible.
        display_bound = max(
            5.0 * MATCH_COMPONENT_TOL_GEV,
            abs(low),
            abs(high),
        )
        low = -display_bound
        high = display_bound

        ax.hist(
            scrambled_values,
            bins=120,
            range=(low, high),
            density=True,
            histtype="step",
            linewidth=1.4,
            label="Scrambled nearest neighbor",
        )
        ax.hist(
            real_values,
            bins=120,
            range=(low, high),
            density=True,
            histtype="step",
            linewidth=1.6,
            label="Real nearest neighbor",
        )

        ax.axvline(
            -MATCH_COMPONENT_TOL_GEV,
            linestyle="--",
            linewidth=1.0,
        )
        ax.axvline(
            MATCH_COMPONENT_TOL_GEV,
            linestyle="--",
            linewidth=1.0,
            label=r"$|\Delta p_i|=0.002$ GeV",
        )
        ax.axvline(0.0, linestyle=":", linewidth=1.0)

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.20)
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3)
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


def plot_clasdis_max_component_delta(
    real_diag: pd.DataFrame,
    scrambled_diag: pd.DataFrame,
    component_tolerance_gev: float,
    output_path: Path,
    title: str,
) -> None:
    real = real_diag["max_abs_component_delta"].to_numpy(dtype=float)
    scrambled = scrambled_diag["max_abs_component_delta"].to_numpy(dtype=float)

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
        5.0 * component_tolerance_gev,
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
            linewidth=1.5,
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
            component_tolerance_gev,
            linestyle="--",
            linewidth=1.4,
            label=(
                f"Current component cut "
                f"{component_tolerance_gev:.3f} GeV"
            ),
        )
        ax.set_xlabel(
            r"Maximum Cartesian component mismatch "
            r"$\max_i |\Delta p_i|$ (GeV)"
        )
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


def accepted_match_rows(
    diagnostics: pd.DataFrame,
    max_distance: float = MATCH_MAX_SCALED_DISTANCE,
) -> pd.DataFrame:
    """
    Return rows satisfying the complete parent-match prescription.

    Prefer the explicit boolean 'accepted' column, which includes BOTH the
    scaled-distance requirement and the individual Cartesian-component cut.
    Fall back to the distance cut only for older diagnostic frames.
    """
    if "accepted" in diagnostics.columns:
        return diagnostics.loc[
            diagnostics["accepted"].to_numpy(dtype=bool)
        ].copy()
    #endif

    return diagnostics.loc[
        diagnostics["match_distance"].to_numpy(dtype=float) <= max_distance
    ].copy()


def plot_clasdis_match_multiplicity(
    real_diag: pd.DataFrame,
    n_eppi0: int,
    max_distance: float,
    output_path: Path,
    title: str,
) -> None:
    accepted = accepted_match_rows(real_diag, max_distance=max_distance)

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


def plot_direct_eppi0_detector_topology(
    eppi0: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    """
    Detector topology of the reconstructed pi0 sample itself.

    This plot uses ONLY detector_gamma1 and detector_gamma2 from the e'p'pi0X
    tree. No cross-tree matching enters this diagnostic.
    """
    g1 = eppi0["detector_gamma1"].to_numpy(dtype=int)
    g2 = eppi0["detector_gamma2"].to_numpy(dtype=int)
    topology = pi0_detector_topology(g1, g2)

    topology_labels = ["FT-FT", "FT-FD", "FD-FD", "other"]
    topology_counts = [
        int(np.count_nonzero(topology == label))
        for label in topology_labels
    ]

    all_daughters = np.concatenate([g1, g2])
    daughter_labels = ["FT (0)", "FD (1)", "Other"]
    daughter_counts = [
        int(np.count_nonzero(all_daughters == 0)),
        int(np.count_nonzero(all_daughters == 1)),
        int(np.count_nonzero(~np.isin(all_daughters, [0, 1]))),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.3))

    axes[0].bar(topology_labels, topology_counts)
    axes[0].set_ylabel("Reconstructed eppi0X candidates")
    axes[0].set_title("Direct pi0 detector topology")
    axes[0].grid(axis="y", alpha=0.20)

    for i, count in enumerate(topology_counts):
        axes[0].text(i, count, f"{count:,}", ha="center", va="bottom")
    #endfor

    axes[1].bar(daughter_labels, daughter_counts)
    axes[1].set_ylabel("Reconstructed pi0 daughter photons")
    axes[1].set_title("Direct daughter-detector composition")
    axes[1].grid(axis="y", alpha=0.20)

    for i, count in enumerate(daughter_counts):
        axes[1].text(i, count, f"{count:,}", ha="center", va="bottom")
    #endfor

    n_total = len(eppi0)
    fig.suptitle(
        f"{title}\n"
        f"N(pi0 candidates)={n_total:,}; classification comes directly from "
        f"detector_gamma1/detector_gamma2"
    )
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_all_vs_matched_epgamma_detector_composition(
    epg: pd.DataFrame,
    matched: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    """
    Compare detector2 for every CLASDIS e'p'gammaX row with detector2 for the
    subset passing the e/p parent match.
    """
    labels = ["FT (0)", "FD (1)", "Other"]

    def counts(df: pd.DataFrame) -> np.ndarray:
        detector = df["detector2"].to_numpy(dtype=int)
        return np.asarray(
            [
                np.count_nonzero(detector == 0),
                np.count_nonzero(detector == 1),
                np.count_nonzero(~np.isin(detector, [0, 1])),
            ],
            dtype=int,
        )
    #enddef

    all_counts = counts(epg)
    matched_counts = counts(matched)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.3))

    x = np.arange(len(labels))
    width = 0.38

    axes[0].bar(x - width / 2.0, all_counts, width=width, label="All CLASDIS epgammaX")
    axes[0].bar(x + width / 2.0, matched_counts, width=width, label="Matched subset")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("Candidate count")
    axes[0].set_title("Absolute counts")
    axes[0].grid(axis="y", alpha=0.20)
    axes[0].legend()

    all_frac = all_counts / max(np.sum(all_counts), 1)
    matched_frac = matched_counts / max(np.sum(matched_counts), 1)

    axes[1].bar(x - width / 2.0, all_frac, width=width, label="All CLASDIS epgammaX")
    axes[1].bar(x + width / 2.0, matched_frac, width=width, label="Matched subset")
    axes[1].set_xticks(x, labels)
    axes[1].set_ylabel("Fraction")
    axes[1].set_ylim(0.0, 1.0)
    axes[1].set_title("Normalized composition")
    axes[1].grid(axis="y", alpha=0.20)
    axes[1].legend()

    fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_match_multiplicity_by_pi0_topology(
    real_diag: pd.DataFrame,
    eppi0: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    """
    For each reconstructed eppi0X candidate, count how many accepted epgammaX
    candidates map to it, and separate that multiplicity by the pi0's DIRECT
    detector_gamma1/detector_gamma2 topology.
    """
    accepted = accepted_match_rows(real_diag)

    n_pi0 = len(eppi0)
    nearest = accepted["_nearest_eppi0_index"].to_numpy(dtype=int)
    multiplicity = np.bincount(nearest, minlength=max(n_pi0, 1))

    pi0_topology = pi0_detector_topology(
        eppi0["detector_gamma1"].to_numpy(dtype=int),
        eppi0["detector_gamma2"].to_numpy(dtype=int),
    )

    labels = ["FT-FT", "FT-FD", "FD-FD", "other"]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))

    max_mult = max(5, min(10, int(np.quantile(multiplicity, 0.995)) + 1))
    bins = np.arange(-0.5, max_mult + 1.5, 1.0)

    for label in labels:
        vals = multiplicity[pi0_topology == label]
        if len(vals) == 0:
            continue
        #endif
        vals = np.minimum(vals, max_mult)
        axes[0].hist(
            vals,
            bins=bins,
            histtype="step",
            linewidth=1.6,
            label=f"{label} (Npi0={len(vals):,})",
        )
    #endfor

    axes[0].set_xlabel(
        "Accepted epgammaX candidates per reconstructed pi0"
        + (f" (last bin >= {max_mult})" if np.any(multiplicity > max_mult) else "")
    )
    axes[0].set_ylabel("Number of eppi0X candidates")
    axes[0].set_yscale("log")
    axes[0].set_title("Multiplicity by direct pi0 topology")
    axes[0].grid(alpha=0.20)
    axes[0].legend()

    mean_mult = []
    zero_frac = []
    for label in labels:
        vals = multiplicity[pi0_topology == label]
        if len(vals) == 0:
            mean_mult.append(0.0)
            zero_frac.append(0.0)
        else:
            mean_mult.append(float(np.mean(vals)))
            zero_frac.append(float(np.mean(vals == 0)))
        #endif
    #endfor

    x = np.arange(len(labels))
    axes[1].bar(x - 0.18, mean_mult, width=0.36, label="Mean matched multiplicity")
    axes[1].bar(x + 0.18, zero_frac, width=0.36, label="Fraction with zero matches")
    axes[1].set_xticks(x, labels)
    axes[1].set_ylabel("Value")
    axes[1].set_title("Matching efficiency/multiplicity by topology")
    axes[1].grid(axis="y", alpha=0.20)
    axes[1].legend()

    fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, output_path)


def build_two_photon_closure_dataframe(
    epg: pd.DataFrame,
    eppi0: pd.DataFrame,
    real_diag: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build closure quantities for reconstructed pi0 candidates with exactly TWO
    accepted e'p'gammaX partners.

    The two matched p2 photons are treated as massless four-vectors. We compare:
      * their reconstructed M(gamma gamma) to stored Mh_gammagamma;
      * their summed three-momentum to the stored eppi0X p2 pi0 vector;
      * their detector2 pair to detector_gamma1/detector_gamma2.
    """
    accepted = accepted_match_rows(real_diag)

    grouped = accepted.groupby("_nearest_eppi0_index", sort=False)

    rows = []

    for pi0_index, group in grouped:
        if len(group) != 2:
            continue
        #endif

        epg_indices = group["_epg_index"].to_numpy(dtype=int)
        photons = epg.iloc[epg_indices].copy().reset_index(drop=True)

        p = photons["p2_p"].to_numpy(dtype=float)
        theta = photons["p2_theta"].to_numpy(dtype=float)
        phi = photons["p2_phi"].to_numpy(dtype=float)

        px, py, pz = spherical_to_cartesian(
            p=p,
            theta=theta,
            phi=phi,
            theta_name="p2_theta",
            phi_name="p2_phi",
        )

        energy = p.copy()  # massless photons

        e_sum = float(np.sum(energy))
        px_sum = float(np.sum(px))
        py_sum = float(np.sum(py))
        pz_sum = float(np.sum(pz))

        m2 = (
            e_sum * e_sum
            - px_sum * px_sum
            - py_sum * py_sum
            - pz_sum * pz_sum
        )
        mgg = math.sqrt(max(m2, 0.0))

        pi0_row = eppi0.iloc[int(pi0_index)]

        pi0_px, pi0_py, pi0_pz = spherical_to_cartesian(
            p=np.asarray([pi0_row["p2_p"]], dtype=float),
            theta=np.asarray([pi0_row["p2_theta"]], dtype=float),
            phi=np.asarray([pi0_row["p2_phi"]], dtype=float),
            theta_name="p2_theta",
            phi_name="p2_phi",
        )

        epg_detectors = np.sort(
            photons["detector2"].to_numpy(dtype=int)
        )
        stored_detectors = np.sort(
            np.asarray(
                [
                    pi0_row["detector_gamma1"],
                    pi0_row["detector_gamma2"],
                ],
                dtype=int,
            )
        )

        rows.append(
            {
                "pi0_index": int(pi0_index),
                "mgg_from_two_epgamma": mgg,
                "stored_Mh_gammagamma": float(pi0_row["Mh_gammagamma"]),
                "delta_mgg": mgg - float(pi0_row["Mh_gammagamma"]),
                "delta_pi0_px": px_sum - float(pi0_px[0]),
                "delta_pi0_py": py_sum - float(pi0_py[0]),
                "delta_pi0_pz": pz_sum - float(pi0_pz[0]),
                "epg_detector_a": int(epg_detectors[0]),
                "epg_detector_b": int(epg_detectors[1]),
                "stored_detector_a": int(stored_detectors[0]),
                "stored_detector_b": int(stored_detectors[1]),
                "detector_pair_matches": bool(
                    np.array_equal(epg_detectors, stored_detectors)
                ),
                "stored_pi0_topology": str(
                    pi0_detector_topology(
                        np.asarray([pi0_row["detector_gamma1"]], dtype=int),
                        np.asarray([pi0_row["detector_gamma2"]], dtype=int),
                    )[0]
                ),
            }
        )
    #endfor

    return pd.DataFrame(rows)


def plot_two_photon_closure(
    closure: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    if len(closure) == 0:
        return
    #endif

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 9.0))
    axes = axes.ravel()

    axes[0].hist(
        closure["mgg_from_two_epgamma"],
        bins=80,
        range=(0.05, 0.25),
        histtype="step",
        linewidth=1.6,
        label="From two matched epgammaX photons",
    )
    axes[0].hist(
        closure["stored_Mh_gammagamma"],
        bins=80,
        range=(0.05, 0.25),
        histtype="step",
        linewidth=1.6,
        label="Stored eppi0X Mh_gammagamma",
    )
    axes[0].set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    axes[0].set_ylabel("Candidates")
    axes[0].grid(alpha=0.20)
    axes[0].legend()

    axes[1].hist(
        closure["delta_mgg"],
        bins=100,
        histtype="step",
        linewidth=1.6,
    )
    axes[1].set_xlabel(
        r"$M_{\gamma\gamma}^{2\,ep\gamma}-M_{\gamma\gamma}^{eppi0}$ (GeV)"
    )
    axes[1].set_ylabel("Candidates")
    axes[1].grid(alpha=0.20)

    axes[2].scatter(
        closure["stored_Mh_gammagamma"],
        closure["mgg_from_two_epgamma"],
        s=8,
        alpha=0.35,
    )
    axes[2].plot([0.10, 0.17], [0.10, 0.17], linestyle="--")
    axes[2].set_xlim(0.10, 0.17)
    axes[2].set_ylim(0.10, 0.17)
    axes[2].set_xlabel("Stored eppi0X Mgg (GeV)")
    axes[2].set_ylabel("Two matched epgammaX Mgg (GeV)")
    axes[2].grid(alpha=0.20)

    for ax, column, label in [
        (axes[3], "delta_pi0_px", r"$\Delta p_{x,\pi^0}$ (GeV)"),
        (axes[4], "delta_pi0_py", r"$\Delta p_{y,\pi^0}$ (GeV)"),
        (axes[5], "delta_pi0_pz", r"$\Delta p_{z,\pi^0}$ (GeV)"),
    ]:
        ax.hist(
            closure[column],
            bins=100,
            histtype="step",
            linewidth=1.6,
        )
        ax.set_xlabel(label)
        ax.set_ylabel("Candidates")
        ax.grid(alpha=0.20)
    #endfor

    detector_match_fraction = float(
        np.mean(closure["detector_pair_matches"].to_numpy(dtype=bool))
    )

    fig.suptitle(
        f"{title}\n"
        f"N(exactly-two matched photons)={len(closure):,}; "
        f"detector-pair agreement={100.0*detector_match_fraction:.2f}%"
    )
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_two_photon_detector_pair_confusion(
    closure: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    if len(closure) == 0:
        return
    #endif

    labels = ["FT-FT", "FT-FD", "FD-FD", "other"]

    def pair_topology(a: int, b: int) -> str:
        return str(
            pi0_detector_topology(
                np.asarray([a], dtype=int),
                np.asarray([b], dtype=int),
            )[0]
        )
    #enddef

    reconstructed = np.asarray(
        [
            pair_topology(a, b)
            for a, b in zip(
                closure["epg_detector_a"],
                closure["epg_detector_b"],
            )
        ],
        dtype=object,
    )
    stored = closure["stored_pi0_topology"].to_numpy(dtype=object)

    matrix = np.zeros((4, 4), dtype=int)
    for i, stored_label in enumerate(labels):
        for j, reco_label in enumerate(labels):
            matrix[i, j] = int(
                np.count_nonzero(
                    (stored == stored_label)
                    & (reconstructed == reco_label)
                )
            )
        #endfor
    #endfor

    fig, ax = plt.subplots(figsize=(7.0, 6.3))
    image = ax.imshow(matrix, aspect="auto")

    ax.set_xticks(np.arange(4), labels)
    ax.set_yticks(np.arange(4), labels)
    ax.set_xlabel("Topology from the two matched epgammaX detector2 values")
    ax.set_ylabel("Stored eppi0X detector_gamma1/gamma2 topology")
    ax.set_title(title)

    for i in range(4):
        for j in range(4):
            ax.text(
                j,
                i,
                f"{matrix[i, j]:,}",
                ha="center",
                va="center",
            )
        #endfor
    #endfor

    fig.colorbar(image, ax=ax, label="Candidate count")
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_clasdis_detector_composition(
    real_diag: pd.DataFrame,
    max_distance: float,
    output_path: Path,
    title: str,
) -> None:
    accepted = accepted_match_rows(real_diag, max_distance=max_distance)

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
        axes[1].set_xlabel("Nearest eppi0X direct detector_gamma1/gamma2 topology")
        axes[1].set_ylabel("epgammaX p2 detector")
        axes[1].set_title("Candidate-level parent association")

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

    thresholds = [0.5, 1.0, 1.5, 2.0]
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




def plot_ablation_auc_summary(
    results: Dict[str, float],
    output_path: Path,
    title: str,
) -> None:
    labels = list(results.keys())
    auc_values = np.asarray([results[label] for label in labels], dtype=float)

    nominal_auc = float(results["nominal"])
    delta_auc = auc_values - nominal_auc

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.5))

    x = np.arange(len(labels))

    axes[0].bar(x, auc_values)
    axes[0].axhline(nominal_auc, linestyle="--", linewidth=1.2)
    axes[0].set_xticks(x, labels, rotation=35, ha="right")
    axes[0].set_ylabel("Held-out ROC AUC")
    axes[0].set_ylim(
        max(0.5, float(np.min(auc_values)) - 0.02),
        min(1.0, float(np.max(auc_values)) + 0.01),
    )
    axes[0].set_title("Classifier performance")
    axes[0].grid(axis="y", alpha=0.20)

    axes[1].bar(x, delta_auc)
    axes[1].axhline(0.0, linestyle="--", linewidth=1.2)
    axes[1].set_xticks(x, labels, rotation=35, ha="right")
    axes[1].set_ylabel(r"$\Delta$AUC relative to nominal")
    axes[1].set_title("Loss/gain relative to nominal")
    axes[1].grid(axis="y", alpha=0.20)

    for ax in axes:
        for tick in ax.get_xticklabels():
            tick.set_fontsize(9)
        #endfor
    #endfor

    fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_ablation_score_overlays(
    score_map: Dict[str, Tuple[np.ndarray, np.ndarray]],
    output_path: Path,
    title: str,
) -> None:
    """
    Compare AAOgen and DVCSgen test-score shapes for nominal and selected
    ablation models. To avoid an unreadable plot, show nominal, minimal4, and
    the two ablations with the largest AUC loss (supplied by caller).
    """
    names = list(score_map.keys())
    n = len(names)

    fig, axes = plt.subplots(
        nrows=n,
        ncols=1,
        figsize=(8.5, max(4.0, 3.2 * n)),
        sharex=True,
    )
    axes = np.atleast_1d(axes)

    for ax, name in zip(axes, names):
        neg_scores, pos_scores = score_map[name]

        ax.hist(
            neg_scores,
            bins=50,
            range=(0.0, 1.0),
            density=True,
            histtype="step",
            linewidth=1.5,
            label="DVCSgen",
        )
        ax.hist(
            pos_scores,
            bins=50,
            range=(0.0, 1.0),
            density=True,
            histtype="step",
            linewidth=1.5,
            label="AAOgen",
        )

        ax.set_ylabel("Density")
        ax.set_title(name)
        ax.grid(alpha=0.20)
        ax.legend()
    #endfor

    axes[-1].set_xlabel("BDT pi0 score")
    fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, output_path)


def binned_score_statistics(
    values: np.ndarray,
    scores: np.ndarray,
    edges: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    scores = np.asarray(scores, dtype=float)

    centers = 0.5 * (edges[:-1] + edges[1:])
    median = np.full(len(centers), np.nan)
    mean = np.full(len(centers), np.nan)
    q16 = np.full(len(centers), np.nan)
    q84 = np.full(len(centers), np.nan)

    for i in range(len(centers)):
        if i == len(centers) - 1:
            mask = (values >= edges[i]) & (values <= edges[i + 1])
        else:
            mask = (values >= edges[i]) & (values < edges[i + 1])
        #endif

        if np.count_nonzero(mask) < 5:
            continue
        #endif

        selected = scores[mask]
        median[i] = np.median(selected)
        mean[i] = np.mean(selected)
        q16[i], q84[i] = np.quantile(selected, [0.16, 0.84])
    #endfor

    return centers, mean, median, np.vstack([q16, q84])


def plot_clasdis_transfer_vs_kinematics(
    control: pd.DataFrame,
    control_scores: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    """
    Quantify whether the independent CLASDIS pi0-enriched control remains
    high-score across photon energy and polar angle.
    """
    if len(control) == 0:
        return
    #endif

    energy = control["p2_p"].to_numpy(dtype=float)
    theta_raw = control["p2_theta"].to_numpy(dtype=float)

    theta_unit = infer_angle_unit(theta_raw, "p2_theta")
    if theta_unit == "rad":
        theta_deg_values = np.rad2deg(theta_raw)
    else:
        theta_deg_values = theta_raw
    #endif

    score = np.asarray(control_scores, dtype=float)

    energy_finite = energy[np.isfinite(energy)]
    theta_finite = theta_deg_values[np.isfinite(theta_deg_values)]

    if len(energy_finite) == 0 or len(theta_finite) == 0:
        return
    #endif

    e_lo, e_hi = np.quantile(energy_finite, [0.01, 0.99])
    t_lo, t_hi = np.quantile(theta_finite, [0.01, 0.99])

    e_edges = np.linspace(e_lo, e_hi, 9)
    t_edges = np.linspace(t_lo, t_hi, 9)

    e_centers, e_mean, e_median, e_band = binned_score_statistics(
        energy, score, e_edges
    )
    t_centers, t_mean, t_median, t_band = binned_score_statistics(
        theta_deg_values, score, t_edges
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))

    for ax, centers, mean, median, band, xlabel in [
        (
            axes[0],
            e_centers,
            e_mean,
            e_median,
            e_band,
            r"$E_\gamma$ / $p_\gamma$ (GeV)",
        ),
        (
            axes[1],
            t_centers,
            t_mean,
            t_median,
            t_band,
            r"$\theta_\gamma$ (deg)",
        ),
    ]:
        valid = np.isfinite(median)

        ax.plot(
            centers[valid],
            median[valid],
            marker="o",
            label="Median score",
        )
        ax.plot(
            centers[valid],
            mean[valid],
            marker="s",
            linestyle="--",
            label="Mean score",
        )
        ax.fill_between(
            centers[valid],
            band[0, valid],
            band[1, valid],
            alpha=0.20,
            label="16–84% interval",
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel("CLASDIS BDT pi0 score")
        ax.set_ylim(0.0, 1.0)
        ax.grid(alpha=0.20)
        ax.legend()
    #endfor

    fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, output_path)


def plot_clasdis_vs_aaogen_transfer(
    test: pd.DataFrame,
    test_scores: np.ndarray,
    control: pd.DataFrame,
    control_scores: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    """
    Directly compare the score distribution of AAOgen pi0 positives with the
    independent CLASDIS pi0-enriched control.
    """
    if len(control) == 0:
        return
    #endif

    y = test["label"].to_numpy(dtype=int)
    aaogen_scores = test_scores[y == 1]

    fig, ax = plt.subplots(figsize=(8.0, 5.5))

    ax.hist(
        aaogen_scores,
        bins=50,
        range=(0.0, 1.0),
        density=True,
        histtype="step",
        linewidth=1.7,
        label=f"AAOgen held-out (N={len(aaogen_scores):,})",
    )
    ax.hist(
        control_scores,
        bins=50,
        range=(0.0, 1.0),
        density=True,
        histtype="step",
        linewidth=1.7,
        label=f"CLASDIS matched pi0 control (N={len(control_scores):,})",
    )

    ax.set_xlabel("BDT pi0 score")
    ax.set_ylabel("Normalized density")
    ax.set_xlim(0.0, 1.0)
    ax.grid(alpha=0.20)
    ax.legend()
    ax.set_title(title)

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
    nominal_output = region_output / "topology"
    ablation_output = region_output / "ablations"

    nominal_output.mkdir(parents=True, exist_ok=True)
    ablation_output.mkdir(parents=True, exist_ok=True)

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
        f"REGION {region}: training pool — "
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
            f"REGION {region}: CLASDIS matched control — "
            f"{len(control_region):,} candidates"
        )
    #endif

    # Use the NOMINAL feature list to define the common finite sample and one
    # common train/validation/test split for all ablation models. This makes AUC
    # differences directly interpretable.
    working = clean_feature_rows(
        combined,
        ABLATION_FEATURE_SETS["nominal"],
    )

    if control_region is not None:
        control_clean = clean_feature_rows(
            control_region,
            ABLATION_FEATURE_SETS["nominal"],
        )
    else:
        control_clean = None
    #endif

    train, validation, test = train_validation_test_split(
        working,
        args.seed,
    )

    progress(
        f"REGION {region}: common split — "
        f"train={len(train):,}, validation={len(validation):,}, "
        f"test={len(test):,}"
    )

    # ---------------------------------------------------------------------
    # Nominal topology model
    # ---------------------------------------------------------------------
    features = ABLATION_FEATURE_SETS["nominal"]

    progress(f"REGION {region}: training nominal topology BDT")

    plot_input_feature_distributions(
        train=train,
        features=features,
        output_path=nominal_output / "input_feature_distributions.png",
        title=f"{args.period} {region}: topology-BDT training inputs",
    )

    nominal_model = train_model(
        train=train,
        validation=validation,
        features=features,
        seed=args.seed,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
    )

    test_scores = model_scores(
        nominal_model,
        test,
        features,
    )

    if control_clean is not None:
        control_scores = model_scores(
            nominal_model,
            control_clean,
            features,
        )
    else:
        control_scores = None
    #endif

    roc_auc = plot_roc_curve(
        test=test,
        test_scores=test_scores,
        output_path=nominal_output / "roc_curve.png",
        title=f"{args.period} {region}: nominal topology BDT",
    )

    plot_score_distribution(
        test=test,
        test_scores=test_scores,
        control=control_clean,
        control_scores=control_scores,
        output_path=nominal_output / "bdt_score_distribution.png",
        title=f"{args.period} {region}: nominal topology BDT score",
    )

    plot_feature_importance(
        model=nominal_model,
        test=test,
        features=features,
        output_path=nominal_output / "feature_importance.png",
        title=f"{args.period} {region}: permutation feature importance",
        seed=args.seed,
    )

    balanced_accuracy = plot_confusion_matrix_at_half(
        test=test,
        test_scores=test_scores,
        output_path=nominal_output / "confusion_matrix_score_0p5.png",
        title=f"{args.period} {region}: nominal topology BDT",
    )

    if (
        control_clean is not None
        and control_scores is not None
        and len(control_clean) > 0
    ):
        plot_clasdis_vs_aaogen_transfer(
            test=test,
            test_scores=test_scores,
            control=control_clean,
            control_scores=control_scores,
            output_path=nominal_output / "clasdis_vs_aaogen_score.png",
            title=f"{args.period} {region}: AAOgen vs CLASDIS pi0-score transfer",
        )

        plot_clasdis_transfer_vs_kinematics(
            control=control_clean,
            control_scores=control_scores,
            output_path=nominal_output / "clasdis_score_vs_kinematics.png",
            title=f"{args.period} {region}: CLASDIS score stability vs photon kinematics",
        )
    #endif

    joblib.dump(
        {
            "model": nominal_model,
            "features": list(features),
            "period": args.period,
            "region": region,
            "feature_set": "topology_nominal",
        },
        nominal_output / "bdt_model.joblib",
    )

    # ---------------------------------------------------------------------
    # Ablation models
    # ---------------------------------------------------------------------
    ablation_auc: Dict[str, float] = {"nominal": roc_auc}
    ablation_score_cache: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    y_test = test["label"].to_numpy(dtype=int)
    ablation_score_cache["nominal"] = (
        test_scores[y_test == 0],
        test_scores[y_test == 1],
    )

    if not args.skip_ablations:
        for index, (name, ablation_features) in enumerate(
            ABLATION_FEATURE_SETS.items()
        ):
            if name == "nominal":
                continue
            #endif

            progress(
                f"REGION {region}: ablation '{name}' — "
                f"{len(ablation_features)} features"
            )

            model = train_model(
                train=train,
                validation=validation,
                features=ablation_features,
                seed=args.seed + index + 1,
                n_estimators=args.n_estimators,
                max_depth=args.max_depth,
            )

            scores = model_scores(
                model,
                test,
                ablation_features,
            )

            fpr, tpr, _ = roc_curve(
                y_test,
                scores,
            )
            model_auc = float(auc(fpr, tpr))

            ablation_auc[name] = model_auc
            ablation_score_cache[name] = (
                scores[y_test == 0],
                scores[y_test == 1],
            )

            progress(
                f"REGION {region}: ablation '{name}' AUC={model_auc:.5f} "
                f"(delta={model_auc-roc_auc:+.5f})"
            )
        #endfor

        plot_ablation_auc_summary(
            results=ablation_auc,
            output_path=ablation_output / "ablation_auc_summary.png",
            title=f"{args.period} {region}: feature-ablation performance",
        )

        # Show nominal, minimal4, and the two non-minimal ablations with the
        # largest degradation relative to nominal.
        degradation = sorted(
            [
                (name, value - roc_auc)
                for name, value in ablation_auc.items()
                if name not in {"nominal", "minimal4"}
            ],
            key=lambda item: item[1],
        )

        overlay_names = ["nominal"]
        if "minimal4" in ablation_score_cache:
            overlay_names.append("minimal4")
        #endif
        overlay_names.extend(
            [name for name, _ in degradation[:2]]
        )

        overlay_map = {
            name: ablation_score_cache[name]
            for name in overlay_names
            if name in ablation_score_cache
        }

        plot_ablation_score_overlays(
            score_map=overlay_map,
            output_path=ablation_output / "ablation_score_overlays.png",
            title=f"{args.period} {region}: score-shape sensitivity to feature removal",
        )
    #endif

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
        f"REGION {region}: complete — nominal AUC={roc_auc:.5f}, "
        f"balanced accuracy@0.5={balanced_accuracy:.5f} "
        f"({elapsed_since(region_start)})"
    )

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

        match_output = period_output / "clasdis_matching"
        match_output.mkdir(parents=True, exist_ok=True)

        progress(
            "CLASDIS DIAGNOSTICS: plotting DIRECT reconstructed-pi0 detector "
            "topology from detector_gamma1/detector_gamma2"
        )
        plot_direct_eppi0_detector_topology(
            eppi0=clasdis_eppi0,
            output_path=match_output / "direct_eppi0_detector_topology.png",
            title=f"{args.period}: direct reconstructed-pi0 detector topology",
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

        progress(
            "CLASDIS DIAGNOSTICS: comparing all CLASDIS epgammaX detector "
            "composition to accepted matched candidates"
        )
        plot_all_vs_matched_epgamma_detector_composition(
            epg=clasdis_epg,
            matched=clasdis_control,
            output_path=match_output / "all_vs_matched_epgamma_detector_composition.png",
            title=f"{args.period}: CLASDIS epgammaX detector composition",
        )

        progress(
            f"CLASDIS DIAGNOSTICS: exact-two-photon closure sample contains "
            f"{len(closure):,} reconstructed pi0 candidates"
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
