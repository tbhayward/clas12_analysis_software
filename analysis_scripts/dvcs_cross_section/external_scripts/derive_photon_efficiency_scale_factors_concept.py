#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_concept.py

Fast, standalone proof-of-concept script for the RGA photon-efficiency study.

Purpose
-------
Make only the first-stage epgamma shape-comparison canvases, with essentially
no analysis machinery beyond the minimum needed to compare the samples.

For each of the five RGA periods and separately for FT and FD photons, make a
3x8 canvas comparing:

    1. Mx2
    2. Mx2_1
    3. Mx2_2
    4. Emiss2
    5. E_tag (= p2_p, because p2 is always the photon)
    6. Delta_phi (tree branch, displayed as a residual about 180 degrees)
    7. pTmiss
    8. theta_gamma_gamma

Samples:
    * nSidis data
    * exclusive-pi0 AAO MC reconstructed as epgamma
    * BH/DVCS DVCSgen MC

CLASDIS is intentionally not used anywhere in this script.

The first row uses ONLY

    angle(e, gamma) > 5 degrees,

apart from assigning the reconstructed photon to FT or FD by its polar angle.

The second row cumulatively adds

    Mx2_1 < 0.15 GeV^2.

The third row cumulatively adds

    Emiss2 > 1.0 GeV.

The FT/FD angular assignment is:

    FT: 2.4 <= theta_gamma < 5.0 degrees
    FD: 6.0 <= theta_gamma < 35.0 degrees

There are deliberately no missing-mass, missing-energy, pTmiss, Delta_phi,
electron-momentum, tag-energy, or other exclusivity cuts.

Efficiency design
-----------------
This script is intentionally optimized for speed:

    * ROOT files are streamed in chunks with uproot.iterate().
    * Events are never accumulated in large in-memory feature tables.
    * Histograms are filled immediately and the chunk is discarded.
    * A given ROOT sample is read only once; FT and FD histograms are filled
      simultaneously from the same chunk.
    * Run periods are processed concurrently with ProcessPoolExecutor.
    * ROOT files within one period are still read sequentially to avoid nested
      I/O contention and excessive memory pressure.
    * Only the 12 branches actually needed here are read.
    * No eppi0 files, event association, CLASDIS, Stage-II, Stage-III, or
      grand diagnostics are touched.
    * A lightweight AAO-vs-DVCS template-fit closure test is run after every
      cumulative optimizer step to quantify whether the surviving templates
      remain distinguishable enough to fit their mixture fraction.

Examples
--------
Run all five periods concurrently with up to 4 million events per ROOT sample:

    python derive_photon_efficiency_scale_factors_concept_v2.py --workers 5

Run one period:

    python derive_photon_efficiency_scale_factors_concept.py \
        --period fa18_inb

Use all entries:

    python derive_photon_efficiency_scale_factors_concept.py \
        --max-entries 0

By default all shape-comparison panels use a logarithmic y-axis. Add

    --linear

for a linear y-axis.

The script runs the greedy rectangular-cut optimizer BEFORE making the
shape-comparison canvases. The optimizer is NEVER combined across periods:
each of the five run periods is optimized independently, and FT/FD are
optimized independently within each period.

Before optimization, a loose exclusivity preselection is imposed:

    Mx2_1 < 0.30 GeV^2

by default. The real-data interval

    0.20 <= Mx2_1 < 0.30 GeV^2

is used as an empirical nonexclusive/SIDIS-like sideband. Because Mx2_1
defines both the loose preselection and the sideband, Mx2_1 itself is excluded
from the greedy threshold scan. The optimizer therefore seeks reconstructed
cuts that simultaneously retain AAO pi0, suppress DVCSgen, and suppress this
data sideband.

After the optimizer selects its period- and detector-specific cut sequence,
the shape-comparison canvas is constructed with one cumulative row per selected
cut. The same selected cuts are then shown on data, AAO, and DVCSgen.

The optimizer is run independently in four predicted-probe-energy bins:

    0.4 <= E_probe < 1.0 GeV
    1.0 <= E_probe < 2.0 GeV
    2.0 <= E_probe < 4.0 GeV
    4.0 <= E_probe < 9.5 GeV

with E_probe identified here with the stored Emiss2 quantity used by the
missing-photon construction. Because E_probe now defines the optimization bin,
Emiss2 is removed completely from the optimizer variables: the optimizer may
not sculpt the probe-energy spectrum inside or across the requested bins.

Shape canvases are written under:

    output/photon_efficiency_concept/stage1_shape_comparison/

Optimizer CSV/TXT/progression plots are written under:

    output/photon_efficiency_concept/rectangular_optimizer/

Template-fit closure CSV/TXT/summary canvases are written under:

    output/photon_efficiency_concept/template_fit_closure/

For every cumulative cut step, the closure test automatically chooses the
single surviving reconstructed variable with the largest AAO-vs-DVCS
Jensen-Shannon divergence, builds independent training and held-out templates,
creates pseudo-data at known pi0 fractions, and refits those fractions. This
tests whether increasing purity has made the remaining AAO and DVCS shapes too
degenerate to determine their normalization from data.
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


# =============================================================================
# Configuration
# =============================================================================

TREE_DEFAULT = "PhysicsEvents"

THETA_EGAMMA_MIN_DEG = 5.0

# Keep the same explicit photon angular regions used in the current study.
FT_THETA_MIN_DEG = 2.4
FT_THETA_MAX_DEG = 5.0
FD_THETA_MIN_DEG = 6.0
FD_THETA_MAX_DEG = 35.0


# Predicted missing-photon energy bins used for the purity optimization.
# These are intentionally kept fixed for this concept iteration.
#
# tuple format:
#   (key, lower edge [GeV], upper edge [GeV])
E_PROBE_BINS: Tuple[Tuple[str, float, float], ...] = (
    ("0p4_1p0", 0.4, 1.0),
    ("1p0_2p0", 1.0, 2.0),
    ("2p0_4p0", 2.0, 4.0),
    ("4p0_9p5", 4.0, 9.5),
)

# Only these branches are read from the very large epgamma ROOT files.
BRANCHES: Tuple[str, ...] = (
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    "Emiss2",
    "Delta_phi",
    "pTmiss",
    "theta_gamma_gamma",
    "e_theta",
    "e_phi",
    "p2_p",
    "p2_theta",
    "p2_phi",
)

# Plot ranges are intentionally centralized here so that the concept study is
# trivial to iterate on.
#
# tuple format:
#   (branch/key, title, x label, x_min, x_max, n_bins)
PLOT_SPECS = (
    (
        "Mx2",
        r"$MM^2(ep\gamma X)$",
        r"$MM^2(ep\gamma X)$ (GeV$^2$)",
        -0.05,
        0.10,
        100,
    ),
    (
        "Mx2_1",
        r"$MM^2(epX)$",
        r"$MM^2(epX)$ (GeV$^2$)",
        -0.30,
        0.70,
        140,
    ),
    (
        "Mx2_2",
        r"$MM^2(e\gamma X)$",
        r"$MM^2(e\gamma X)$ (GeV$^2$)",
        0.00,
        5.00,
        140,
    ),
    (
        "Emiss2",
        r"$E_{\rm miss}(ep\gamma X)$",
        r"$E_{\rm miss}(ep\gamma X)$ (GeV)",
        -1.0,
        5.0,
        120,
    ),
    (
        "E_tag",
        r"$E_{\rm tag}$",
        r"$E_{\rm tag}$ (GeV)",
        0.0,
        9.5,
        120,
    ),
    (
        "Delta_phi_residual_deg",
        r"$\Delta\phi(p,\gamma)$",
        r"$\Delta\phi(p,\gamma)$ residual from $180^\circ$ (deg)",
        -25.0,
        25.0,
        120,
    ),
    (
        "pTmiss",
        r"$p_{T,\rm miss}$",
        r"$p_{T,\rm miss}$ (GeV)",
        0.0,
        1.60,
        160,
    ),
    (
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$",
        r"$\theta_{\gamma\gamma}$ (deg)",
        0.0,
        6.0,
        120,
    ),
)


@dataclass(frozen=True)
class Period:
    key: str
    label: str
    data: str
    pi0_mc: str
    dvcs_mc: str


_BASE = (
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
    "dvcsgen_files_greater_than_0.40GeV"
)

_DATA_BASE = (
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
    "dvcs/efficiency_study"
)

PERIODS: Tuple[Period, ...] = (
    Period(
        "fa18_inb",
        "Fa18 Inb",
        f"{_DATA_BASE}/nSidis_rga_fa18_inb_epgamma.root",
        f"{_BASE}/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root",
    ),
    Period(
        "fa18_out",
        "Fa18 Out",
        f"{_DATA_BASE}/nSidis_rga_fa18_out_epgamma.root",
        f"{_BASE}/bkg_rga_fa18_out_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_fa18_out_epgamma_0.40GeV.root",
    ),
    Period(
        "sp18_inb",
        "Sp18 Inb",
        f"{_DATA_BASE}/nSidis_rga_sp18_inb_epgamma.root",
        f"{_BASE}/bkg_rga_sp18_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp18_inb_epgamma_0.40GeV.root",
    ),
    Period(
        "sp18_out",
        "Sp18 Out",
        f"{_DATA_BASE}/nSidis_rga_sp18_out_epgamma.root",
        f"{_BASE}/bkg_rga_sp18_out_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp18_out_epgamma_0.40GeV.root",
    ),
    Period(
        "sp19_inb",
        "Sp19 Inb",
        f"{_DATA_BASE}/nSidis_rga_sp19_inb_epgamma.root",
        f"{_BASE}/bkg_rga_sp19_inb_epgamma_0.40GeV.root",
        f"{_BASE}/dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root",
    ),
)


SAMPLES = (
    ("data", "data", "black"),
    ("pi0", r"exclusive $\pi^0$ (AAO)", "tab:red"),
    ("dvcs", "BH/DVCS (DVCSgen)", "tab:blue"),
)


# Reconstructed variables available to the first rectangular-cut optimizer.
# The optimizer is run separately for every period and independently in FT/FD.
OPTIMIZER_FEATURES: Tuple[str, ...] = (
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    "Emiss2",
    "E_tag",
    "Delta_phi_residual_deg",
    "pTmiss",
    "theta_gamma_gamma",
)

# Mx2_1 defines the loose exclusivity preselection/data sideband, while
# Emiss2 defines the E_probe bin. Both are deliberately excluded from the
# greedy threshold scan.
OPTIMIZER_SCAN_FEATURES: Tuple[str, ...] = tuple(
    feature
    for feature in OPTIMIZER_FEATURES
    if feature not in ("Mx2_1", "Emiss2")
)

# Every remaining optimizer variable may use either one-sided direction.
# Emiss2/E_probe is not present here because it defines the energy bin itself.
OPTIMIZER_ALLOWED_DIRECTIONS = {
    feature: ("lt", "gt")
    for feature in OPTIMIZER_SCAN_FEATURES
}


# =============================================================================
# Small utilities
# =============================================================================

def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fast standalone Stage-I photon-efficiency shape comparison: "
            "data vs AAO pi0 vs DVCSgen, separately in FT and FD."
        )
    )
    parser.add_argument(
        "--period",
        action="append",
        choices=[p.key for p in PERIODS],
        help=(
            "Run only this period. May be supplied more than once. "
            "Default: all five RGA periods."
        ),
    )
    parser.add_argument(
        "--tree",
        default=TREE_DEFAULT,
        help=f"ROOT tree name. Default: {TREE_DEFAULT}.",
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=4_000_000,
        help=(
            "Maximum entries read from each ROOT sample. "
            "Use 0 for the entire tree. Default: 4000000."
        ),
    )
    parser.add_argument(
        "--step-size",
        type=int,
        default=500_000,
        help="Entries per uproot chunk. Default: 500000.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=5,
        help=(
            "Maximum number of run periods processed concurrently. "
            "Useful range is 1-5 because there are five RGA periods. "
            "Default: 5."
        ),
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help=(
            "Use a linear y-axis for the shape-comparison canvases. "
            "Default: logarithmic y-axis."
        ),
    )
    parser.add_argument(
        "--no-optimize-cuts",
        action="store_true",
        help=(
            "Disable the iterative rectangular-cut optimizer. "
            "Default: run it independently for every selected period and "
            "separately in FT and FD."
        ),
    )
    parser.add_argument(
        "--optimizer-max-events",
        type=int,
        default=200_000,
        help=(
            "Maximum events retained in memory per sample, period, detector region, "
            "and E_probe bin for cut optimization. A uniform "
            "priority reservoir is used while streaming. Default: 200000."
        ),
    )
    parser.add_argument(
        "--optimizer-steps",
        type=int,
        default=4,
        help="Maximum number of greedy rectangular cuts. Default: 4.",
    )
    parser.add_argument(
        "--optimizer-quantiles",
        type=int,
        default=80,
        help=(
            "Number of quantile positions scanned per variable and direction "
            "at each optimizer step. Default: 80."
        ),
    )
    parser.add_argument(
        "--optimizer-min-step-pi0-eff",
        type=float,
        default=0.75,
        help=(
            "Each newly added cut must retain at least this fraction of the "
            "AAO events surviving the previous step. Default: 0.75."
        ),
    )
    parser.add_argument(
        "--optimizer-min-total-pi0-eff",
        type=float,
        default=0.20,
        help=(
            "Never accept a cut that reduces cumulative AAO retention below "
            "this fraction of the baseline. Default: 0.20."
        ),
    )
    parser.add_argument(
        "--optimizer-min-improvement",
        type=float,
        default=1.02,
        help=(
            "Stop if the best new cut improves the current separation score "
            "by less than this multiplicative factor. Default: 1.02."
        ),
    )
    parser.add_argument(
        "--optimizer-mx2-1-preselection-max",
        type=float,
        default=0.30,
        help=(
            "Loose exclusivity preselection Mx2_1 < value in GeV^2, applied "
            "before optimization and in the post-optimizer shape canvases. "
            "Default: 0.30."
        ),
    )
    parser.add_argument(
        "--optimizer-data-sideband-min",
        type=float,
        default=0.20,
        help=(
            "Lower Mx2_1 edge in GeV^2 of the real-data sideband used as an "
            "empirical nonexclusive/SIDIS-like background proxy. Default: 0.20."
        ),
    )
    parser.add_argument(
        "--optimizer-data-sideband-max",
        type=float,
        default=0.30,
        help=(
            "Upper Mx2_1 edge in GeV^2 of the real-data sideband. Default: 0.30."
        ),
    )
    parser.add_argument(
        "--optimizer-data-sideband-weight",
        type=float,
        default=1.0,
        help=(
            "Relative weight lambda of data-sideband retention in the "
            "optimizer score denominator. Default: 1.0."
        ),
    )
    parser.add_argument(
        "--closure-toys",
        type=int,
        default=100,
        help=(
            "Pseudoexperiments per true pi0 fraction and cumulative optimizer "
            "step for template-fit closure. Default: 100."
        ),
    )
    parser.add_argument(
        "--closure-events",
        type=int,
        default=5000,
        help=(
            "Reference pseudo-data event count used in each closure toy. "
            "This fixes the statistical scale for comparing identifiability "
            "between cut steps. Default: 5000."
        ),
    )
    parser.add_argument(
        "--closure-bins",
        type=int,
        default=24,
        help=(
            "Maximum quantile bins for the automatically selected 1D closure "
            "discriminator. Default: 24."
        ),
    )
    parser.add_argument(
        "--output",
        default="output/photon_efficiency_concept/stage1_shape_comparison",
        help=(
            "Output directory. Default: "
            "output/photon_efficiency_concept/stage1_shape_comparison"
        ),
    )
    return parser.parse_args()


def get_tree(root_file: uproot.ReadOnlyDirectory, requested: str):
    if requested in root_file:
        return root_file[requested], requested
    #endif

    # Be forgiving if the requested name is absent, but do not silently choose
    # an arbitrary non-tree object.
    for key, classname in root_file.classnames().items():
        if classname.startswith("TTree") or "RNTuple" in classname:
            clean = key.split(";")[0]
            return root_file[key], clean
        #endif
    #endfor

    raise RuntimeError(
        f"No TTree/RNTuple found in ROOT file; requested '{requested}'."
    )


def preflight_file(path: str, tree_name: str) -> Tuple[str, int]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    #endif

    with uproot.open(path) as root_file:
        tree, found_name = get_tree(root_file, tree_name)
        available = set(tree.keys())
        missing = [branch for branch in BRANCHES if branch not in available]
        if missing:
            raise RuntimeError(
                f"{path}: tree '{found_name}' is missing required branches: "
                + ", ".join(missing)
            )
        #endif
        return found_name, int(tree.num_entries)
    #endwith


def infer_angle_unit(
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    g_theta: np.ndarray,
    g_phi: np.ndarray,
) -> dict:
    # The production files have historically contained angular variables in
    # either radians or degrees depending on the processing path. Infer the
    # convention locally and robustly.
    probe = np.concatenate(
        [
            np.abs(np.asarray(e_theta, dtype=float)[:10000]),
            np.abs(np.asarray(e_phi, dtype=float)[:10000]),
            np.abs(np.asarray(g_theta, dtype=float)[:10000]),
            np.abs(np.asarray(g_phi, dtype=float)[:10000]),
        ]
    )
    probe = probe[np.isfinite(probe)]
    if probe.size == 0:
        return "rad"
    #endif
    return "deg" if float(np.nanpercentile(probe, 99.0)) > 7.0 else "rad"


def opening_angle_deg(
    theta1: np.ndarray,
    phi1: np.ndarray,
    theta2: np.ndarray,
    phi2: np.ndarray,
    unit: str,
) -> np.ndarray:
    if unit == "deg":
        theta1 = np.deg2rad(theta1)
        phi1 = np.deg2rad(phi1)
        theta2 = np.deg2rad(theta2)
        phi2 = np.deg2rad(phi2)
    else:
        theta1 = np.asarray(theta1, dtype=float)
        phi1 = np.asarray(phi1, dtype=float)
        theta2 = np.asarray(theta2, dtype=float)
        phi2 = np.asarray(phi2, dtype=float)
    #endif

    cos_angle = (
        np.cos(theta1) * np.cos(theta2)
        + np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
    )
    return np.rad2deg(np.arccos(np.clip(cos_angle, -1.0, 1.0)))


def photon_theta_deg(values: np.ndarray, unit: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values if unit == "deg" else np.rad2deg(values)


def delta_phi_residual_deg(values: np.ndarray) -> np.ndarray:
    """
    Use the stored ROOT-tree Delta_phi branch; do NOT reconstruct Delta_phi
    from particle azimuths here.

    For plotting only, express that stored coplanarity angle as its signed
    residual about 180 degrees, so a perfectly back-to-back p-gamma pair is
    displayed at 0 degrees. The current processing stores Delta_phi in radians,
    centered near pi. A defensive degree fallback is included for portability.
    """
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size and float(np.nanpercentile(np.abs(finite[:10000]), 99.0)) > 7.0:
        residual = (values - 180.0 + 180.0) % 360.0 - 180.0
        return residual
    #endif

    residual_rad = (values - math.pi + math.pi) % (2.0 * math.pi) - math.pi
    return np.rad2deg(residual_rad)


def empty_histograms() -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
    out: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    for region in ("FT", "FD"):
        out[region] = {}
        for row in ("minimal", "mx2_1_cut", "emiss_cut"):
            out[region][row] = {}
            for key, _title, _xlabel, x_min, x_max, n_bins in PLOT_SPECS:
                out[region][row][key] = np.zeros(n_bins, dtype=np.int64)
            #endfor
        #endfor
    #endfor
    return out


def histogram_edges() -> Dict[str, np.ndarray]:
    return {
        key: np.linspace(x_min, x_max, n_bins + 1)
        for key, _title, _xlabel, x_min, x_max, n_bins in PLOT_SPECS
    }



# =============================================================================
# Rectangular-cut optimizer support
# =============================================================================

class PriorityReservoir:
    """
    Uniform fixed-size event sample maintained with random priorities.

    Every accepted event receives an independent U(0,1) priority. We retain
    only the K smallest priorities seen so far. This is a vectorized uniform
    reservoir across all streamed chunks and avoids keeping the full event
    sample in memory.
    """

    def __init__(self, capacity: int, n_features: int, seed: int) -> None:
        self.capacity = int(capacity)
        self.n_features = int(n_features)
        self.rng = np.random.default_rng(int(seed))
        self.priorities = np.empty(0, dtype=np.float32)
        self.values = np.empty((0, n_features), dtype=np.float32)

    def update(self, values: np.ndarray) -> None:
        if self.capacity <= 0 or values.size == 0:
            return
        #endif

        values = np.asarray(values, dtype=np.float32)
        finite = np.all(np.isfinite(values), axis=1)
        values = values[finite]
        if values.shape[0] == 0:
            return
        #endif

        priorities = self.rng.random(values.shape[0]).astype(np.float32)

        if self.values.shape[0] == 0:
            combined_values = values
            combined_priorities = priorities
        else:
            combined_values = np.concatenate((self.values, values), axis=0)
            combined_priorities = np.concatenate(
                (self.priorities, priorities),
                axis=0,
            )
        #endif

        if combined_values.shape[0] <= self.capacity:
            self.values = combined_values
            self.priorities = combined_priorities
            return
        #endif

        keep = np.argpartition(
            combined_priorities,
            self.capacity - 1,
        )[: self.capacity]
        self.values = combined_values[keep]
        self.priorities = combined_priorities[keep]

    def array(self) -> np.ndarray:
        return np.asarray(self.values, dtype=np.float32)


def optimizer_matrix(values: Dict[str, np.ndarray], mask: np.ndarray) -> np.ndarray:
    """Build the reconstructed-feature matrix used by the optimizer."""
    if not np.any(mask):
        return np.empty((0, len(OPTIMIZER_FEATURES)), dtype=np.float32)
    #endif

    return np.column_stack(
        [
            np.asarray(values[key][mask], dtype=np.float32)
            for key in OPTIMIZER_FEATURES
        ]
    )


def _candidate_thresholds(
    pi_values: np.ndarray,
    dvcs_values: np.ndarray,
    n_quantiles: int,
) -> np.ndarray:
    """
    Quantile-based thresholds avoid wasting scan points in empty phase space.
    """
    combined = np.concatenate((pi_values, dvcs_values))
    combined = combined[np.isfinite(combined)]

    if combined.size < 10:
        return np.empty(0, dtype=float)
    #endif

    q = np.linspace(0.02, 0.98, max(3, int(n_quantiles)))
    thresholds = np.quantile(combined, q)
    return np.unique(np.asarray(thresholds, dtype=float))


def _candidate_thresholds_three(
    pi_values: np.ndarray,
    dvcs_values: np.ndarray,
    side_values: np.ndarray,
    n_quantiles: int,
) -> np.ndarray:
    """Quantile thresholds from surviving AAO, DVCSgen, and data sideband."""
    combined = np.concatenate((pi_values, dvcs_values, side_values))
    combined = combined[np.isfinite(combined)]

    if combined.size < 10:
        return np.empty(0, dtype=float)
    #endif

    q = np.linspace(0.02, 0.98, max(3, int(n_quantiles)))
    return np.unique(np.asarray(np.quantile(combined, q), dtype=float))


def optimize_rectangular_cuts(
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    data_sideband_events: np.ndarray,
    max_steps: int,
    n_quantiles: int,
    min_step_pi0_eff: float,
    min_total_pi0_eff: float,
    min_improvement: float,
    data_sideband_weight: float,
) -> list:
    """
    Greedy one-sided optimizer balancing AAO retention against two backgrounds.

    At each step it scans one-sided rectangular thresholds over
    OPTIMIZER_SCAN_FEATURES and maximizes

        F = eps_pi0 / sqrt(
            eps_DVCS + lambda * eps_data_sideband + 1e-9
        )

    using cumulative efficiencies relative to the corresponding baseline
    populations.

    Mx2_1 is excluded from the scan because it defines both the loose
    exclusivity preselection and the real-data sideband.

    Emiss2 is also excluded from the scan because it defines the E_probe bin.
    Thus the optimizer cannot improve purity by discarding part of the
    requested photon-energy interval.
    """
    if (
        pi0_events.shape[0] == 0
        or dvcs_events.shape[0] == 0
        or data_sideband_events.shape[0] == 0
    ):
        return []
    #endif

    pi_mask = np.ones(pi0_events.shape[0], dtype=bool)
    dvcs_mask = np.ones(dvcs_events.shape[0], dtype=bool)
    side_mask = np.ones(data_sideband_events.shape[0], dtype=bool)

    n_pi0_0 = int(pi0_events.shape[0])
    n_dvcs_0 = int(dvcs_events.shape[0])
    n_side_0 = int(data_sideband_events.shape[0])

    current_score = 1.0 / math.sqrt(1.0 + float(data_sideband_weight))
    used_features = set()
    results = []

    for step in range(1, int(max_steps) + 1):
        n_pi_before = int(np.count_nonzero(pi_mask))
        n_dvcs_before = int(np.count_nonzero(dvcs_mask))
        n_side_before = int(np.count_nonzero(side_mask))

        if min(n_pi_before, n_dvcs_before, n_side_before) == 0:
            break
        #endif

        best = None

        for feature in OPTIMIZER_SCAN_FEATURES:
            if feature in used_features:
                continue
            #endif

            ifeature = FEATURE_INDEX[feature]
            thresholds = _candidate_thresholds_three(
                pi0_events[pi_mask, ifeature],
                dvcs_events[dvcs_mask, ifeature],
                data_sideband_events[side_mask, ifeature],
                n_quantiles,
            )

            for threshold in thresholds:
                for direction in OPTIMIZER_ALLOWED_DIRECTIONS[feature]:
                    if direction == "lt":
                        pi_local = pi0_events[:, ifeature] < threshold
                        dvcs_local = dvcs_events[:, ifeature] < threshold
                        side_local = data_sideband_events[:, ifeature] < threshold
                    else:
                        pi_local = pi0_events[:, ifeature] > threshold
                        dvcs_local = dvcs_events[:, ifeature] > threshold
                        side_local = data_sideband_events[:, ifeature] > threshold
                    #endif

                    pi_test = pi_mask & pi_local
                    dvcs_test = dvcs_mask & dvcs_local
                    side_test = side_mask & side_local

                    n_pi = int(np.count_nonzero(pi_test))
                    n_dvcs = int(np.count_nonzero(dvcs_test))
                    n_side = int(np.count_nonzero(side_test))

                    step_pi_eff = n_pi / n_pi_before
                    step_dvcs_eff = n_dvcs / n_dvcs_before
                    step_side_eff = n_side / n_side_before

                    total_pi_eff = n_pi / n_pi0_0
                    total_dvcs_eff = n_dvcs / n_dvcs_0
                    total_side_eff = n_side / n_side_0

                    if step_pi_eff < min_step_pi0_eff:
                        continue
                    #endif
                    if total_pi_eff < min_total_pi0_eff:
                        continue
                    #endif

                    score = total_pi_eff / math.sqrt(
                        total_dvcs_eff
                        + float(data_sideband_weight) * total_side_eff
                        + 1.0e-9
                    )

                    if best is None or score > best["score"]:
                        best = {
                            "step": step,
                            "feature": feature,
                            "direction": direction,
                            "threshold": float(threshold),
                            "score": float(score),
                            "step_pi0_eff": float(step_pi_eff),
                            "step_dvcs_eff": float(step_dvcs_eff),
                            "step_data_sideband_eff": float(step_side_eff),
                            "total_pi0_eff": float(total_pi_eff),
                            "total_dvcs_eff": float(total_dvcs_eff),
                            "total_data_sideband_eff": float(total_side_eff),
                            "n_pi0": n_pi,
                            "n_dvcs": n_dvcs,
                            "n_data_sideband": n_side,
                            "pi_mask": pi_test,
                            "dvcs_mask": dvcs_test,
                            "side_mask": side_test,
                        }
                    #endif
                #endfor
            #endfor
        #endfor

        if best is None:
            break
        #endif

        improvement = best["score"] / max(current_score, 1.0e-12)
        if improvement < min_improvement:
            break
        #endif

        results.append(
            {
                key: value
                for key, value in best.items()
                if key not in ("pi_mask", "dvcs_mask", "side_mask")
            }
        )

        pi_mask = best["pi_mask"]
        dvcs_mask = best["dvcs_mask"]
        side_mask = best["side_mask"]
        current_score = best["score"]
        used_features.add(best["feature"])
    #endfor

    return results


def cut_expression(result: dict) -> str:
    op = "<" if result["direction"] == "lt" else ">"
    return f"{result['feature']} {op} {result['threshold']:.6g}"



# =============================================================================
# Template-fit closure / identifiability study
# =============================================================================

CLOSURE_TRUE_FRACTIONS: Tuple[float, ...] = (0.20, 0.50, 0.80)


def jensen_shannon_divergence(
    p: np.ndarray,
    q: np.ndarray,
) -> float:
    """Jensen-Shannon divergence in natural-log units."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)

    p = p / np.sum(p)
    q = q / np.sum(q)
    m = 0.5 * (p + q)

    return float(
        0.5 * np.sum(p * np.log(p / m))
        + 0.5 * np.sum(q * np.log(q / m))
    )


def closure_histogram_probabilities(
    values: np.ndarray,
    edges: np.ndarray,
    pseudocount: float = 0.5,
) -> np.ndarray:
    """
    Histogram probability vector with a small Jeffreys-like pseudocount so
    finite MC statistics never produce impossible zero-probability bins.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    counts, _ = np.histogram(values, bins=edges)
    probs = counts.astype(float) + float(pseudocount)
    return probs / np.sum(probs)


def closure_quantile_edges(
    pi_values: np.ndarray,
    dvcs_values: np.ndarray,
    n_bins: int,
) -> np.ndarray:
    """
    Quantile bins from the combined training templates.

    Quantile binning gives approximately useful occupancy everywhere and is
    more stable than a fixed-width histogram for long-tailed variables.
    """
    combined = np.concatenate((pi_values, dvcs_values))
    combined = combined[np.isfinite(combined)]

    if combined.size < 20:
        return np.empty(0, dtype=float)
    #endif

    edges = np.quantile(
        combined,
        np.linspace(0.0, 1.0, max(4, int(n_bins)) + 1),
    )
    edges = np.unique(np.asarray(edges, dtype=float))

    if edges.size < 5:
        return np.empty(0, dtype=float)
    #endif

    # Include the extrema robustly in np.histogram.
    edges[0] = np.nextafter(edges[0], -np.inf)
    edges[-1] = np.nextafter(edges[-1], np.inf)
    return edges


def fit_two_template_fraction(
    counts: np.ndarray,
    p_pi0: np.ndarray,
    p_dvcs: np.ndarray,
) -> Tuple[float, float]:
    """
    Maximum-likelihood fit of

        p_i(f) = f p_pi0,i + (1-f) p_dvcs,i

    to multinomial-binned pseudo-data.

    The score equation is monotonic in f, so a bounded bisection gives the
    exact one-parameter maximum without scipy.  The returned uncertainty is
    the local inverse-curvature estimate.
    """
    counts = np.asarray(counts, dtype=float)
    p_pi0 = np.asarray(p_pi0, dtype=float)
    p_dvcs = np.asarray(p_dvcs, dtype=float)
    delta = p_pi0 - p_dvcs

    def gradient(f_value: float) -> float:
        mixture = p_dvcs + f_value * delta
        mixture = np.clip(mixture, 1.0e-15, None)
        return float(np.sum(counts * delta / mixture))

    g0 = gradient(0.0)
    g1 = gradient(1.0)

    if g0 <= 0.0:
        f_hat = 0.0
    elif g1 >= 0.0:
        f_hat = 1.0
    else:
        lo = 0.0
        hi = 1.0
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            if gradient(mid) > 0.0:
                lo = mid
            else:
                hi = mid
            #endif
        #endfor
        f_hat = 0.5 * (lo + hi)
    #endif

    mixture = p_dvcs + f_hat * delta
    mixture = np.clip(mixture, 1.0e-15, None)
    information = float(
        np.sum(counts * (delta * delta) / (mixture * mixture))
    )
    sigma = (
        1.0 / math.sqrt(information)
        if information > 0.0
        else float("inf")
    )

    return float(f_hat), float(sigma)


def build_cumulative_optimizer_masks(
    events: np.ndarray,
    optimizer_results: list,
) -> list:
    """Baseline plus one cumulative mask after each optimizer-selected cut."""
    masks = [np.ones(events.shape[0], dtype=bool)]
    current = masks[0]

    for result in optimizer_results:
        current = apply_optimizer_cut(events, current, result)
        masks.append(current.copy())
    #endfor

    return masks


def run_template_fit_closure(
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    optimizer_results: list,
    n_toys: int,
    n_pseudodata_events: int,
    n_bins: int,
    seed: int,
) -> list:
    """
    Test AAO/DVCS mixture identifiability at every cumulative cut step.

    Procedure at each step:
      1. Apply all cuts through that step.
      2. Split surviving AAO and DVCS events into independent template-training
         and held-out pseudo-data halves.
      3. Among the allowed reconstructed variables, find the single 1D
         variable with the largest AAO-vs-DVCS Jensen-Shannon divergence in the
         training templates.
      4. Build pseudo-data from the held-out shapes at true pi0 fractions
         0.20, 0.50, and 0.80.
      5. Fit each pseudo-data sample with the independent training templates.

    Large fit RMS/bias after later cuts means that purity optimization has made
    the residual AAO and DVCS templates too degenerate for a reliable fraction
    extraction.
    """
    if (
        pi0_events.shape[0] < 40
        or dvcs_events.shape[0] < 40
        or n_toys <= 0
        or n_pseudodata_events <= 0
    ):
        return []
    #endif

    rng = np.random.default_rng(int(seed))

    # Split once before cut application so the training/held-out samples remain
    # statistically independent at every cumulative optimizer step.
    pi_perm = rng.permutation(pi0_events.shape[0])
    dvcs_perm = rng.permutation(dvcs_events.shape[0])

    pi_train_flag = np.zeros(pi0_events.shape[0], dtype=bool)
    dvcs_train_flag = np.zeros(dvcs_events.shape[0], dtype=bool)
    pi_train_flag[pi_perm[: pi_perm.size // 2]] = True
    dvcs_train_flag[dvcs_perm[: dvcs_perm.size // 2]] = True

    pi_test_flag = ~pi_train_flag
    dvcs_test_flag = ~dvcs_train_flag

    pi_masks = build_cumulative_optimizer_masks(
        pi0_events,
        optimizer_results,
    )
    dvcs_masks = build_cumulative_optimizer_masks(
        dvcs_events,
        optimizer_results,
    )

    closure_results = []

    for step_index, (pi_mask, dvcs_mask) in enumerate(
        zip(pi_masks, dvcs_masks)
    ):
        pi_train = pi_mask & pi_train_flag
        pi_test = pi_mask & pi_test_flag
        dvcs_train = dvcs_mask & dvcs_train_flag
        dvcs_test = dvcs_mask & dvcs_test_flag

        n_pi_train = int(np.count_nonzero(pi_train))
        n_pi_test = int(np.count_nonzero(pi_test))
        n_dvcs_train = int(np.count_nonzero(dvcs_train))
        n_dvcs_test = int(np.count_nonzero(dvcs_test))

        if min(n_pi_train, n_pi_test, n_dvcs_train, n_dvcs_test) < 20:
            continue
        #endif

        best_feature = None
        best_edges = None
        best_p_pi_train = None
        best_p_dvcs_train = None
        best_js = -1.0

        # Emiss2 and Mx2_1 are still valid closure diagnostics even though they
        # are not optimizer variables.  The closure study asks only whether
        # any surviving reconstructed 1D shape retains normalization leverage.
        for feature in OPTIMIZER_FEATURES:
            ifeature = FEATURE_INDEX[feature]

            pi_values = pi0_events[pi_train, ifeature]
            dvcs_values = dvcs_events[dvcs_train, ifeature]

            edges = closure_quantile_edges(
                pi_values,
                dvcs_values,
                n_bins,
            )
            if edges.size == 0:
                continue
            #endif

            p_pi_train = closure_histogram_probabilities(
                pi_values,
                edges,
            )
            p_dvcs_train = closure_histogram_probabilities(
                dvcs_values,
                edges,
            )

            js = jensen_shannon_divergence(
                p_pi_train,
                p_dvcs_train,
            )

            if js > best_js:
                best_js = js
                best_feature = feature
                best_edges = edges
                best_p_pi_train = p_pi_train
                best_p_dvcs_train = p_dvcs_train
            #endif
        #endfor

        if best_feature is None:
            continue
        #endif

        ifeature = FEATURE_INDEX[best_feature]
        p_pi_test = closure_histogram_probabilities(
            pi0_events[pi_test, ifeature],
            best_edges,
        )
        p_dvcs_test = closure_histogram_probabilities(
            dvcs_events[dvcs_test, ifeature],
            best_edges,
        )

        for true_fraction in CLOSURE_TRUE_FRACTIONS:
            truth_probability = (
                true_fraction * p_pi_test
                + (1.0 - true_fraction) * p_dvcs_test
            )
            truth_probability /= np.sum(truth_probability)

            fitted = np.empty(n_toys, dtype=float)
            fitted_sigma = np.empty(n_toys, dtype=float)

            for itoy in range(n_toys):
                toy_counts = rng.multinomial(
                    int(n_pseudodata_events),
                    truth_probability,
                )
                fitted[itoy], fitted_sigma[itoy] = (
                    fit_two_template_fraction(
                        toy_counts,
                        best_p_pi_train,
                        best_p_dvcs_train,
                    )
                )
            #endfor

            delta = fitted - float(true_fraction)

            closure_results.append(
                {
                    "step": int(step_index),
                    "feature": str(best_feature),
                    "js_divergence": float(best_js),
                    "true_fraction": float(true_fraction),
                    "mean_fit": float(np.mean(fitted)),
                    "bias": float(np.mean(delta)),
                    "rms_error": float(np.sqrt(np.mean(delta * delta))),
                    "std_fit": float(np.std(fitted, ddof=1)),
                    "mean_sigma": float(np.mean(fitted_sigma)),
                    "boundary_fraction": float(
                        np.mean(
                            (fitted <= 1.0e-6)
                            | (fitted >= 1.0 - 1.0e-6)
                        )
                    ),
                    "n_pi0_train": n_pi_train,
                    "n_pi0_test": n_pi_test,
                    "n_dvcs_train": n_dvcs_train,
                    "n_dvcs_test": n_dvcs_test,
                    "n_toys": int(n_toys),
                    "n_pseudodata_events": int(n_pseudodata_events),
                }
            )
        #endfor
    #endfor

    return closure_results


def write_template_fit_closure_outputs(
    period: Period,
    region: str,
    energy_bin_key: str,
    energy_min: float,
    energy_max: float,
    closure_results: list,
    output_dir: Path,
) -> None:
    """Write detailed closure results for one period/region/E_probe bin."""
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = (
        f"template_fit_closure_{period.key}_{region.lower()}_"
        f"eprobe_{energy_bin_key}"
    )

    csv_path = output_dir / f"{stem}.csv"
    txt_path = output_dir / f"{stem}.txt"

    fields = (
        "step",
        "feature",
        "js_divergence",
        "true_fraction",
        "mean_fit",
        "bias",
        "rms_error",
        "std_fit",
        "mean_sigma",
        "boundary_fraction",
        "n_pi0_train",
        "n_pi0_test",
        "n_dvcs_train",
        "n_dvcs_test",
        "n_toys",
        "n_pseudodata_events",
    )

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in closure_results:
            writer.writerow({field: result[field] for field in fields})
        #endfor
    #endwith

    lines = [
        (
            f"{period.label} {region}: template-fit closure, "
            f"{energy_min:g} <= E_probe < {energy_max:g} GeV"
        ),
        "=" * 88,
        (
            "At each cumulative optimizer step, the single reconstructed "
            "variable with maximum AAO-vs-DVCS Jensen-Shannon divergence is "
            "used for an independent train/held-out two-template closure fit."
        ),
        "",
    ]

    if not closure_results:
        lines.append("Insufficient statistics for closure.")
    else:
        steps = sorted({int(result["step"]) for result in closure_results})
        for step in steps:
            step_rows = [
                result
                for result in closure_results
                if int(result["step"]) == step
            ]
            if not step_rows:
                continue
            #endif

            lines.append(
                f"Step {step}: best feature={step_rows[0]['feature']}, "
                f"JS={step_rows[0]['js_divergence']:.5f}"
            )

            for result in step_rows:
                lines.append(
                    f"  f_true={result['true_fraction']:.2f}: "
                    f"<f_fit>={result['mean_fit']:.4f}, "
                    f"bias={result['bias']:+.4f}, "
                    f"RMS={result['rms_error']:.4f}, "
                    f"<sigma_fit>={result['mean_sigma']:.4f}, "
                    f"boundary={100.0*result['boundary_fraction']:.1f}%"
                )
            #endfor
            lines.append("")
        #endfor
    #endif

    txt_path.write_text("\n".join(lines) + "\n")


def make_combined_closure_canvas(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """
    One 4x2 closure canvas per period.

    Each panel shows RMS(f_fit-f_true) versus cumulative cut step for the three
    injected pi0 fractions.  Growth of the RMS after later cuts is the direct
    warning that template identifiability is being lost.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        2,
        figsize=(13.8, 4.0 * len(E_PROBE_BINS) + 1.0),
        squeeze=False,
    )

    legend_handles = []
    legend_labels = []

    for irow, (bin_key, e_min, e_max) in enumerate(E_PROBE_BINS):
        bin_summary = optimizer_summary.get("energy_bins", {}).get(
            bin_key,
            {},
        )

        for icol, region in enumerate(("FD", "FT")):
            ax = axes[irow, icol]
            region_summary = bin_summary.get("regions", {}).get(region)
            closure_results = (
                region_summary.get("closure_results", [])
                if region_summary is not None
                else []
            )

            if not closure_results:
                ax.text(
                    0.5,
                    0.5,
                    "Insufficient closure statistics",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(
                    f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV"
                )
                continue
            #endif

            for true_fraction in CLOSURE_TRUE_FRACTIONS:
                rows = [
                    result
                    for result in closure_results
                    if abs(
                        result["true_fraction"] - true_fraction
                    ) < 1.0e-9
                ]
                rows.sort(key=lambda result: result["step"])

                line, = ax.plot(
                    [result["step"] for result in rows],
                    [result["rms_error"] for result in rows],
                    marker="o",
                    label=rf"$f_{{\pi^0}}^{{true}}={true_fraction:.1f}$",
                )

                if irow == 0 and icol == 0:
                    legend_handles.append(line)
                    legend_labels.append(
                        rf"$f_{{\pi^0}}^{{true}}={true_fraction:.1f}$"
                    )
                #endif
            #endfor

            # Annotate which 1D shape carried the most separation at each step.
            step_rows = {}
            for result in closure_results:
                step_rows.setdefault(
                    int(result["step"]),
                    result,
                )
            #endfor

            annotation = "\n".join(
                (
                    f"{step}: {row['feature']} "
                    f"(JS={row['js_divergence']:.3f})"
                )
                for step, row in sorted(step_rows.items())
            )
            ax.text(
                0.98,
                0.97,
                annotation,
                transform=ax.transAxes,
                fontsize=7.4,
                va="top",
                ha="right",
                bbox={
                    "facecolor": "white",
                    "alpha": 0.82,
                    "edgecolor": "0.6",
                },
            )

            ax.set_xlabel("Cumulative cut step")
            ax.set_ylabel(r"closure RMS of $f_{\pi^0}$")
            ax.set_ylim(bottom=0.0)
            ax.grid(alpha=0.25)
            ax.set_title(
                f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV"
            )
        #endfor
    #endfor

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.975),
            ncol=3,
            frameon=True,
        )
    #endif

    fig.suptitle(
        f"{period.label}: AAO/DVCS template-fit closure by cut step",
        fontsize=13,
        y=0.998,
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.985,
        bottom=0.055,
        top=0.93,
        wspace=0.18,
        hspace=0.42,
    )

    out = output_dir / f"template_fit_closure_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def write_optimizer_outputs(
    period: Period,
    region: str,
    energy_bin_key: str,
    energy_min: float,
    energy_max: float,
    pi0_events: np.ndarray,
    dvcs_events: np.ndarray,
    data_sideband_events: np.ndarray,
    results: list,
    output_dir: Path,
) -> None:
    """
    Write detailed per-region CSV/TXT optimizer results.

    Plotting is handled separately so each period gets one combined FD+FT
    progression canvas rather than separate FD and FT PNG files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = (
        f"rectangular_optimizer_{period.key}_{region.lower()}_"
        f"eprobe_{energy_bin_key}"
    )

    csv_path = output_dir / f"{stem}.csv"
    txt_path = output_dir / f"{stem}.txt"

    fields = (
        "step",
        "feature",
        "direction",
        "threshold",
        "step_pi0_eff",
        "step_dvcs_eff",
        "step_data_sideband_eff",
        "total_pi0_eff",
        "total_dvcs_eff",
        "total_data_sideband_eff",
        "score",
        "n_pi0",
        "n_dvcs",
        "n_data_sideband",
    )

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            writer.writerow({key: result[key] for key in fields})
        #endfor
    #endwith

    lines = [
        (
            f"{period.label} {region} rectangular-cut optimizer: "
            f"{energy_min:g} <= E_probe < {energy_max:g} GeV"
        ),
        "=" * 72,
        (
            "Baseline: "
            f"AAO={pi0_events.shape[0]:,}, "
            f"DVCSgen={dvcs_events.shape[0]:,}, "
            f"data sideband={data_sideband_events.shape[0]:,}"
        ),
        (
            "Objective: maximize cumulative eps_pi0 / "
            "sqrt(eps_DVCS + lambda*eps_data_sideband), with configured "
            "AAO-retention constraints."
        ),
        "",
    ]

    if not results:
        lines.append("No accepted optimizer step.")
    else:
        for result in results:
            lines.extend(
                [
                    f"Step {result['step']}: {cut_expression(result)}",
                    (
                        f"  step retention: "
                        f"AAO={100.0*result['step_pi0_eff']:.2f}%  "
                        f"DVCSgen={100.0*result['step_dvcs_eff']:.2f}%  "
                        f"data-SB={100.0*result['step_data_sideband_eff']:.2f}%"
                    ),
                    (
                        f"  cumulative:     "
                        f"AAO={100.0*result['total_pi0_eff']:.2f}%  "
                        f"DVCSgen={100.0*result['total_dvcs_eff']:.2f}%  "
                        f"data-SB={100.0*result['total_data_sideband_eff']:.2f}%  "
                        f"F={result['score']:.4g}"
                    ),
                    "",
                ]
            )
        #endfor
    #endif

    txt_path.write_text("\n".join(lines) + "\n")


def make_combined_optimizer_canvas(
    period: Period,
    optimizer_summary: dict,
    output_dir: Path,
) -> Path:
    """
    One optimizer progression canvas per period.

    Rows are E_probe bins; columns are FD and FT. This keeps all eight
    period-specific optimizations in one compact diagnostic.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(E_PROBE_BINS),
        2,
        figsize=(13.8, 4.0 * len(E_PROBE_BINS) + 1.0),
        sharey=True,
        squeeze=False,
    )

    legend_handles = []
    legend_labels = []

    for irow, (bin_key, e_min, e_max) in enumerate(E_PROBE_BINS):
        bin_summary = optimizer_summary.get("energy_bins", {}).get(bin_key, {})

        for icol, region in enumerate(("FD", "FT")):
            ax = axes[irow, icol]
            region_summary = bin_summary.get("regions", {}).get(region)

            if region_summary is None:
                ax.text(
                    0.5,
                    0.5,
                    "No optimizer result",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(
                    f"{region}: {e_min:g}-{e_max:g} GeV"
                )
                continue
            #endif

            results = region_summary.get("results", [])
            x = [0] + [int(result["step"]) for result in results]
            pi_eff = [1.0] + [
                result["total_pi0_eff"]
                for result in results
            ]
            dvcs_eff = [1.0] + [
                result["total_dvcs_eff"]
                for result in results
            ]
            side_eff = [1.0] + [
                result["total_data_sideband_eff"]
                for result in results
            ]

            line_pi, = ax.plot(
                x,
                pi_eff,
                marker="o",
                label=r"AAO $\pi^0$",
            )
            line_dvcs, = ax.plot(
                x,
                dvcs_eff,
                marker="o",
                label="DVCSgen",
            )
            line_side, = ax.plot(
                x,
                side_eff,
                marker="o",
                label="data sideband",
            )

            if irow == 0 and icol == 0:
                legend_handles = [line_pi, line_dvcs, line_side]
                legend_labels = [
                    r"AAO $\pi^0$",
                    "DVCSgen",
                    "data sideband",
                ]
            #endif

            ax.set_yscale("log")
            ax.set_ylim(1.0e-4, 1.1)
            ax.set_xticks(x)
            ax.set_xlabel("Greedy cut step")
            ax.set_title(
                f"{region}: {e_min:g} <= E_probe < {e_max:g} GeV\n"
                f"AAO={region_summary['n_pi0_baseline']:,}, "
                f"DVCS={region_summary['n_dvcs_baseline']:,}, "
                f"data-SB={region_summary['n_data_sideband_baseline']:,}",
                fontsize=10.0,
            )
            ax.grid(alpha=0.25)

            if results:
                annotation = "\n".join(
                    f"{r['step']}. {cut_expression(r)}"
                    for r in results
                )
                final = results[-1]
                annotation += (
                    "\n"
                    f"final: AAO {100.0*final['total_pi0_eff']:.1f}%, "
                    f"DVCS {100.0*final['total_dvcs_eff']:.2f}%, "
                    f"data-SB "
                    f"{100.0*final['total_data_sideband_eff']:.2f}%"
                )
                ax.text(
                    0.02,
                    0.03,
                    annotation,
                    transform=ax.transAxes,
                    fontsize=7.8,
                    va="bottom",
                    ha="left",
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.88,
                        "edgecolor": "0.5",
                    },
                )
            #endif
        #endfor
    #endfor

    for ax in axes[:, 0]:
        ax.set_ylabel("Cumulative retained fraction")
    #endfor

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.975),
            ncol=3,
            frameon=True,
        )
    #endif

    fig.suptitle(
        f"{period.label}: E_probe-binned rectangular-cut optimization",
        fontsize=13,
        y=0.998,
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.985,
        bottom=0.055,
        top=0.93,
        wspace=0.15,
        hspace=0.42,
    )

    out = output_dir / f"rectangular_optimizer_{period.key}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


# =============================================================================
# Streaming histogrammer
# =============================================================================

def stream_sample_histograms(
    path: str,
    sample_label: str,
    tree_name: str,
    max_entries: int,
    step_size: int,
    collect_optimizer: bool = False,
    optimizer_max_events: int = 0,
    optimizer_seed: int = 12345,
) -> Tuple[
    Dict[str, Dict[str, Dict[str, np.ndarray]]],
    Dict[str, int],
    Dict[str, Dict[str, np.ndarray]],
]:
    """
    Read one ROOT sample exactly once and fill both FT and FD histograms.

    When optimizer collection is enabled, maintain an independent bounded
    uniform reservoir for every detector region AND every E_probe bin. This
    avoids undersampling a sparse energy bin when another energy region is much
    more populated, while still requiring only one ROOT-file pass.
    """
    hist = empty_histograms()
    edges = histogram_edges()
    counts = {
        "read": 0,
        "opening_angle": 0,
        "FT": 0,
        "FD": 0,
        "FT_mx2_1_cut": 0,
        "FD_mx2_1_cut": 0,
        "FT_emiss_cut": 0,
        "FD_emiss_cut": 0,
    }

    region_seed_offset = {"FT": 11, "FD": 29}
    reservoirs = {
        region: {
            bin_key: PriorityReservoir(
                optimizer_max_events,
                len(OPTIMIZER_FEATURES),
                optimizer_seed
                + region_seed_offset[region]
                + 1000 * ibin,
            )
            for ibin, (bin_key, _emin, _emax) in enumerate(E_PROBE_BINS)
        }
        for region in ("FT", "FD")
    }

    with uproot.open(path) as root_file:
        tree, found_tree = get_tree(root_file, tree_name)
        total = int(tree.num_entries)
        entry_stop = total if max_entries <= 0 else min(total, max_entries)

        log(
            f"{sample_label}: tree '{found_tree}', streaming "
            f"{entry_stop:,}/{total:,} entries in {step_size:,}-entry chunks."
        )

        for arrays in tree.iterate(
            expressions=list(BRANCHES),
            entry_start=0,
            entry_stop=entry_stop,
            step_size=step_size,
            library="np",
        ):
            n = len(arrays["Mx2"])
            counts["read"] += n

            angle_unit = infer_angle_unit(
                arrays["e_theta"],
                arrays["e_phi"],
                arrays["p2_theta"],
                arrays["p2_phi"],
            )

            theta_egamma = opening_angle_deg(
                arrays["e_theta"],
                arrays["e_phi"],
                arrays["p2_theta"],
                arrays["p2_phi"],
                angle_unit,
            )
            photon_theta = photon_theta_deg(arrays["p2_theta"], angle_unit)

            base = np.isfinite(theta_egamma) & (
                theta_egamma > THETA_EGAMMA_MIN_DEG
            )
            counts["opening_angle"] += int(np.count_nonzero(base))

            region_masks = {
                "FT": (
                    base
                    & np.isfinite(photon_theta)
                    & (photon_theta >= FT_THETA_MIN_DEG)
                    & (photon_theta < FT_THETA_MAX_DEG)
                ),
                "FD": (
                    base
                    & np.isfinite(photon_theta)
                    & (photon_theta >= FD_THETA_MIN_DEG)
                    & (photon_theta < FD_THETA_MAX_DEG)
                ),
            }

            values = {
                "Mx2": np.asarray(arrays["Mx2"], dtype=float),
                "Mx2_1": np.asarray(arrays["Mx2_1"], dtype=float),
                "Mx2_2": np.asarray(arrays["Mx2_2"], dtype=float),
                "Emiss2": np.asarray(arrays["Emiss2"], dtype=float),
                "E_tag": np.asarray(arrays["p2_p"], dtype=float),
                "Delta_phi_residual_deg": delta_phi_residual_deg(
                    arrays["Delta_phi"]
                ),
                "pTmiss": np.asarray(arrays["pTmiss"], dtype=float),
                "theta_gamma_gamma": np.asarray(
                    arrays["theta_gamma_gamma"],
                    dtype=float,
                ),
            }

            e_probe = values["Emiss2"]

            for region, region_mask in region_masks.items():
                counts[region] += int(np.count_nonzero(region_mask))

                if collect_optimizer:
                    for bin_key, e_min, e_max in E_PROBE_BINS:
                        energy_mask = (
                            region_mask
                            & np.isfinite(e_probe)
                            & (e_probe >= e_min)
                            & (e_probe < e_max)
                        )
                        reservoirs[region][bin_key].update(
                            optimizer_matrix(values, energy_mask)
                        )
                    #endfor
                #endif

                # Retain the old fixed diagnostic histograms only for the
                # --no-optimize-cuts fallback path.
                mx2_1_mask = (
                    region_mask
                    & np.isfinite(values["Mx2_1"])
                    & (values["Mx2_1"] < 0.15)
                )

                emiss_mask = (
                    mx2_1_mask
                    & np.isfinite(values["Emiss2"])
                    & (values["Emiss2"] > 1.0)
                )

                row_masks = {
                    "minimal": region_mask,
                    "mx2_1_cut": mx2_1_mask,
                    "emiss_cut": emiss_mask,
                }

                counts[f"{region}_mx2_1_cut"] += int(
                    np.count_nonzero(mx2_1_mask)
                )
                counts[f"{region}_emiss_cut"] += int(
                    np.count_nonzero(emiss_mask)
                )

                for row_name, row_mask in row_masks.items():
                    for key in values:
                        v = values[key]
                        mask = row_mask & np.isfinite(v)
                        if np.any(mask):
                            h, _ = np.histogram(v[mask], bins=edges[key])
                            hist[region][row_name][key] += h.astype(
                                np.int64,
                                copy=False,
                            )
                        #endif
                    #endfor
                #endfor
            #endfor

            if (
                counts["read"] == entry_stop
                or counts["read"] % max(step_size, 1_000_000) == 0
            ):
                frac = (
                    100.0 * counts["read"] / entry_stop
                    if entry_stop
                    else 100.0
                )
                log(
                    f"{sample_label}: {counts['read']:,}/{entry_stop:,} "
                    f"({frac:.1f}%)"
                )
            #endif
        #endfor
    #endwith

    if collect_optimizer:
        energy_counts = []
        for region in ("FT", "FD"):
            pieces = []
            for bin_key, e_min, e_max in E_PROBE_BINS:
                n_bin = reservoirs[region][bin_key].array().shape[0]
                pieces.append(
                    f"{e_min:g}-{e_max:g} GeV={n_bin:,}"
                )
            #endfor
            energy_counts.append(f"{region}: " + ", ".join(pieces))
        #endfor
        log(
            f"{sample_label}: optimizer reservoirs after "
            f"Angle(e,gamma)>{THETA_EGAMMA_MIN_DEG:g} deg -> "
            + " | ".join(energy_counts)
        )
    else:
        log(
            f"{sample_label}: retained after Angle(e,gamma)>5 deg: "
            f"FT={counts['FT']:,}, FD={counts['FD']:,}."
        )
    #endif

    optimizer_events = {
        region: {
            bin_key: reservoirs[region][bin_key].array()
            for bin_key, _e_min, _e_max in E_PROBE_BINS
        }
        for region in ("FT", "FD")
    }
    return hist, counts, optimizer_events


# =============================================================================
# Dynamic post-optimizer shape-comparison support
# =============================================================================

FEATURE_INDEX = {
    feature: i
    for i, feature in enumerate(OPTIMIZER_FEATURES)
}


def apply_optimizer_cut(
    events: np.ndarray,
    current_mask: np.ndarray,
    result: dict,
) -> np.ndarray:
    """Apply one selected rectangular cut cumulatively to one sample."""
    feature = result["feature"]
    threshold = float(result["threshold"])
    values = events[:, FEATURE_INDEX[feature]]

    if result["direction"] == "lt":
        return current_mask & np.isfinite(values) & (values < threshold)
    #endif

    return current_mask & np.isfinite(values) & (values > threshold)


def optimizer_row_label(result: dict) -> str:
    """Compact human-readable label for one optimizer-selected cut."""
    op = "<" if result["direction"] == "lt" else ">"
    feature = result["feature"]

    pretty = {
        "Mx2": r"$Mx2$",
        "Mx2_1": r"$Mx2_1$",
        "Mx2_2": r"$Mx2_2$",
        "Emiss2": r"$E_{\rm miss}$",
        "E_tag": r"$E_{\rm tag}$",
        "Delta_phi_residual_deg": r"$\Delta\phi(p,\gamma)$",
        "pTmiss": r"$p_{T,\rm miss}$",
        "theta_gamma_gamma": r"$\theta_{\gamma\gamma}$",
    }[feature]

    return rf"$+\,{pretty}{op}{result['threshold']:.3g}$"


def make_optimizer_cut_shape_canvas(
    period: Period,
    region: str,
    energy_bin_key: str,
    energy_min: float,
    energy_max: float,
    sample_events: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    optimizer_results: list,
    output_dir: Path,
    linear_y: bool,
    mx2_1_preselection_max: float,
) -> Path:
    """
    Build one post-optimizer shape canvas for one period, detector region, and
    E_probe bin.

    Row 1 is the energy-bin baseline plus the loose Mx2_1 preselection. Each
    subsequent row cumulatively applies the optimizer-selected cut sequence.
    """
    n_rows = 1 + len(optimizer_results)
    n_cols = len(PLOT_SPECS)

    fig_height = max(4.2, 3.25 * n_rows + 1.8)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(27.5, fig_height),
        squeeze=False,
    )

    sample_meta = {key: (label, color) for key, label, color in SAMPLES}
    edges_by_key = histogram_edges()

    masks_by_sample = {}
    counts_by_sample = {}

    for sample_key, _label, _color in SAMPLES:
        events = sample_events[sample_key][region][energy_bin_key]
        mx2_1_values = events[:, FEATURE_INDEX["Mx2_1"]]
        baseline_mask = (
            np.isfinite(mx2_1_values)
            & (mx2_1_values < mx2_1_preselection_max)
        )
        masks = [baseline_mask]

        current = masks[0]
        for result in optimizer_results:
            current = apply_optimizer_cut(events, current, result)
            masks.append(current.copy())
        #endfor

        masks_by_sample[sample_key] = masks
        counts_by_sample[sample_key] = [
            int(np.count_nonzero(mask))
            for mask in masks
        ]
    #endfor

    legend_handles = []
    legend_labels = []

    row_labels = [
        (
            rf"${energy_min:g}\leq E_{{\rm probe}}<{energy_max:g}$ GeV"
            + "\n"
            + rf"$Mx2_1<{mx2_1_preselection_max:.3g}$ GeV$^2$"
        )
    ]
    row_labels.extend(
        [
            f"step {i}: {optimizer_row_label(result)}"
            for i, result in enumerate(optimizer_results, start=1)
        ]
    )

    for irow in range(n_rows):
        for icol, spec in enumerate(PLOT_SPECS):
            ax = axes[irow, icol]
            key, title, xlabel, _xmin, _xmax, _nbins = spec
            edges = edges_by_key[key]
            centers = 0.5 * (edges[:-1] + edges[1:])
            ifeature = FEATURE_INDEX[key]

            for sample_key, _label, _color in SAMPLES:
                label, color = sample_meta[sample_key]
                events = sample_events[sample_key][region][energy_bin_key]
                mask = masks_by_sample[sample_key][irow]

                values = events[:, ifeature]
                finite = mask & np.isfinite(values)
                h, _ = np.histogram(values[finite], bins=edges)

                y = normalized_hist(h)
                line = ax.step(
                    centers,
                    y,
                    where="mid",
                    linewidth=1.45,
                    color=color,
                    label=label,
                )[0]

                if irow == 0 and icol == 0:
                    legend_handles.append(line)
                    legend_labels.append(label)
                #endif
            #endfor

            if irow == 0:
                ax.set_title(title, fontsize=10.5)
            #endif

            ax.set_xlabel(xlabel, fontsize=8.8)
            if icol == 0:
                ax.set_ylabel(
                    row_labels[irow] + "\nunit-normalized",
                    fontsize=8.8,
                )
            else:
                ax.set_ylabel("unit-normalized", fontsize=8.8)
            #endif

            ax.set_xlim(edges[0], edges[-1])

            if not linear_y:
                ax.set_yscale("log")
                ax.set_ylim(1.0e-3, 1.0)
            #endif

            ax.tick_params(axis="both", labelsize=7.4)
            ax.grid(alpha=0.18)

            if key == "Mx2_1":
                ax.axvline(
                    mx2_1_preselection_max,
                    linestyle="--",
                    linewidth=0.9,
                    color="0.35",
                )
            #endif

            # Emiss2 is the E_probe binning variable, not an optimizer cut.
            if key == "Emiss2":
                ax.axvline(
                    energy_min,
                    linestyle=":",
                    linewidth=1.0,
                    color="0.35",
                )
                ax.axvline(
                    energy_max,
                    linestyle=":",
                    linewidth=1.0,
                    color="0.35",
                )
            #endif

            for applied_result in optimizer_results[:irow]:
                if applied_result["feature"] != key:
                    continue
                #endif
                ax.axvline(
                    float(applied_result["threshold"]),
                    linestyle="--",
                    linewidth=0.9,
                    color="0.35",
                )
            #endfor
        #endfor
    #endfor

    region_text = (
        rf"FT: ${FT_THETA_MIN_DEG:.1f}^\circ\leq\theta_\gamma"
        rf"<{FT_THETA_MAX_DEG:.1f}^\circ$"
        if region == "FT"
        else
        rf"FD: ${FD_THETA_MIN_DEG:.1f}^\circ\leq\theta_\gamma"
        rf"<{FD_THETA_MAX_DEG:.1f}^\circ$"
    )

    final_count_text = ", ".join(
        f"{sample_meta[key][0]}={counts_by_sample[key][-1]:,}"
        for key, _label, _color in SAMPLES
    )

    cut_sequence = (
        " -> ".join(cut_expression(result) for result in optimizer_results)
        if optimizer_results
        else "no optimizer cut accepted"
    )

    fig.suptitle(
        f"{period.label}: {region}, "
        rf"${energy_min:g}\leq E_{{\rm probe}}<{energy_max:g}$ GeV\n"
        rf"baseline: $\mathrm{{Angle}}(e,\gamma)>"
        rf"{THETA_EGAMMA_MIN_DEG:.0f}^\circ$; "
        + rf"$Mx2_1<{mx2_1_preselection_max:.3g}$ GeV$^2$; "
        + region_text
        + "\n"
        + f"optimizer sequence: {cut_sequence}"
        + "\n"
        + f"final sampled counts: {final_count_text}",
        fontsize=10.5,
        y=0.995,
    )

    fig.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.915),
        ncol=3,
        fontsize=9.2,
        frameon=True,
        fancybox=False,
        edgecolor="black",
    )

    fig.subplots_adjust(
        left=0.047,
        right=0.995,
        bottom=0.055,
        top=0.83,
        wspace=0.35,
        hspace=0.48,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out = (
        output_dir
        / (
            f"canvas_shape_comparison_{period.key}_{region.lower()}_"
            f"eprobe_{energy_bin_key}_optimized.png"
        )
    )
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out



# =============================================================================
# Plotting
# =============================================================================

def normalized_hist(counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    total = float(np.sum(counts))
    if total <= 0.0:
        return np.zeros_like(counts)
    #endif
    return counts / total


def make_canvas(
    period: Period,
    region: str,
    sample_hists: Dict[str, Dict[str, Dict[str, Dict[str, np.ndarray]]]],
    sample_counts: Dict[str, Dict[str, int]],
    output_dir: Path,
    linear_y: bool,
) -> Path:
    fig, axes = plt.subplots(3, 8, figsize=(27.5, 11.8))
    edges_by_key = histogram_edges()

    sample_meta = {key: (label, color) for key, label, color in SAMPLES}
    rows = (
        (
            "minimal",
            r"Row 1: only $\mathrm{Angle}(e,\gamma)>5^\circ$",
        ),
        (
            "mx2_1_cut",
            r"Row 2: additionally $Mx2_1<0.15$ GeV$^2$",
        ),
        (
            "emiss_cut",
            r"Row 3: additionally $E_{\rm miss}>1$ GeV",
        ),
    )

    legend_handles = []
    legend_labels = []

    for irow, (row_key, row_label) in enumerate(rows):
        for icol, spec in enumerate(PLOT_SPECS):
            ax = axes[irow, icol]
            key, title, xlabel, _x_min, _x_max, _n_bins = spec
            edges = edges_by_key[key]
            centers = 0.5 * (edges[:-1] + edges[1:])

            for sample_key, _label, _color in SAMPLES:
                label, color = sample_meta[sample_key]
                y = normalized_hist(
                    sample_hists[sample_key][region][row_key][key]
                )
                line = ax.step(
                    centers,
                    y,
                    where="mid",
                    linewidth=1.5,
                    color=color,
                    label=label,
                )[0]

                if irow == 0 and icol == 0:
                    legend_handles.append(line)
                    legend_labels.append(label)
                #endif
            #endfor

            if irow == 0:
                ax.set_title(title, fontsize=11)
            #endif

            ax.set_xlabel(xlabel, fontsize=9.2)
            if icol == 0:
                ax.set_ylabel(row_label + "\nunit-normalized", fontsize=9.2)
            else:
                ax.set_ylabel("unit-normalized", fontsize=9.2)
            #endif

            ax.set_xlim(edges[0], edges[-1])

            if not linear_y:
                # The histograms are unit-normalized. A small positive floor
                # keeps empty bins from causing log-scale warnings while not
                # distorting any nonzero bin content.
                ax.set_yscale("log")
                ax.set_ylim(1.0e-3, 1.0)
            #endif

            ax.tick_params(axis="both", labelsize=7.8)
            ax.grid(alpha=0.18)

            # Show the cumulative selection thresholds on the relevant panels.
            if key == "Mx2_1":
                ax.axvline(
                    0.15,
                    linestyle="--",
                    linewidth=1.0,
                    color="0.35",
                )
            #endif
            if key == "Emiss2":
                ax.axvline(
                    1.0,
                    linestyle="--",
                    linewidth=1.0,
                    color="0.35",
                )
            #endif
        #endfor
    #endfor

    region_text = (
        rf"FT: ${FT_THETA_MIN_DEG:.1f}^\circ\leq\theta_\gamma"
        rf"<{FT_THETA_MAX_DEG:.1f}^\circ$"
        if region == "FT"
        else
        rf"FD: ${FD_THETA_MIN_DEG:.1f}^\circ\leq\theta_\gamma"
        rf"<{FD_THETA_MAX_DEG:.1f}^\circ$"
    )

    row1_counts = ", ".join(
        f"{sample_meta[key][0]}={sample_counts[key][region]:,}"
        for key, _label, _color in SAMPLES
    )
    row2_counts = ", ".join(
        f"{sample_meta[key][0]}="
        f"{sample_counts[key][region + '_mx2_1_cut']:,}"
        for key, _label, _color in SAMPLES
    )
    row3_counts = ", ".join(
        f"{sample_meta[key][0]}="
        f"{sample_counts[key][region + '_emiss_cut']:,}"
        for key, _label, _color in SAMPLES
    )

    fig.suptitle(
        f"{period.label}: {region} epgamma shape comparison\n"
        rf"Row 1: $\mathrm{{Angle}}(e,\gamma)>"
        rf"{THETA_EGAMMA_MIN_DEG:.0f}^\circ$; "
        + region_text
        + r"; Row 2: $+\,Mx2_1<0.15$ GeV$^2$"
        + r"; Row 3: $+\,E_{\rm miss}>1$ GeV"
        + "\n"
        + ("\nLinear y scale" if linear_y else "\nLogarithmic y scale")
        + "\n"
        + f"Counts — row 1: {row1_counts}    |    "
        + f"row 2: {row2_counts}    |    row 3: {row3_counts}",
        fontsize=10.8,
        y=0.992,
    )

    # One boxed canvas-level legend centered directly beneath the title.
    fig.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.905),
        ncol=3,
        fontsize=9.5,
        frameon=True,
        fancybox=False,
        edgecolor="black",
    )

    fig.subplots_adjust(
        left=0.048,
        right=0.995,
        bottom=0.075,
        top=0.82,
        wspace=0.36,
        hspace=0.46,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"canvas_shape_comparison_{period.key}_{region.lower()}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


# =============================================================================
# Period driver
# =============================================================================

def process_period(
    period: Period,
    tree_name: str,
    max_entries: int,
    step_size: int,
    output_dir: Path,
    linear_y: bool,
    run_optimizer: bool,
    optimizer_max_events: int,
    optimizer_steps: int,
    optimizer_quantiles: int,
    optimizer_min_step_pi0_eff: float,
    optimizer_min_total_pi0_eff: float,
    optimizer_min_improvement: float,
    optimizer_mx2_1_preselection_max: float,
    optimizer_data_sideband_min: float,
    optimizer_data_sideband_max: float,
    optimizer_data_sideband_weight: float,
    closure_toys: int,
    closure_events: int,
    closure_bins: int,
) -> dict:
    t0 = time.perf_counter()
    log(f"{period.label}: starting.")

    optimizer_summary = {
        "period_key": period.key,
        "period_label": period.label,
        "regions": {},
    }

    paths = {
        "data": period.data,
        "pi0": period.pi0_mc,
        "dvcs": period.dvcs_mc,
    }

    for sample_key, path in paths.items():
        found_tree, total = preflight_file(path, tree_name)
        log(
            f"{period.label}: preflight {sample_key}: "
            f"{Path(path).name}, tree '{found_tree}', {total:,} entries."
        )
    #endfor

    sample_hists: Dict[
        str,
        Dict[str, Dict[str, Dict[str, np.ndarray]]],
    ] = {}
    sample_counts: Dict[str, Dict[str, int]] = {}
    optimizer_events: Dict[str, Dict[str, np.ndarray]] = {}

    period_index = [p.key for p in PERIODS].index(period.key)

    for sample_index, (sample_key, label, _color) in enumerate(SAMPLES):
        # In optimizer mode, retain a bounded uniform event reservoir for
        # all three samples. AAO and DVCSgen determine the cuts; the data
        # reservoir is used only afterward to visualize those same cuts.
        collect_optimizer = run_optimizer

        hists, counts, opt_events = stream_sample_histograms(
            paths[sample_key],
            f"{period.label} {label}",
            tree_name,
            max_entries,
            step_size,
            collect_optimizer=collect_optimizer,
            optimizer_max_events=optimizer_max_events,
            optimizer_seed=(
                100_000
                + 10_000 * period_index
                + 100 * sample_index
            ),
        )
        sample_hists[sample_key] = hists
        sample_counts[sample_key] = counts

        if collect_optimizer:
            optimizer_events[sample_key] = opt_events
        #endif
    #endfor

    if run_optimizer:
        optimizer_outdir = output_dir.parent / "rectangular_optimizer"

        optimizer_summary["energy_bins"] = {}

        for bin_key, energy_min, energy_max in E_PROBE_BINS:
            optimizer_summary["energy_bins"][bin_key] = {
                "energy_min": float(energy_min),
                "energy_max": float(energy_max),
                "regions": {},
            }

            for region in ("FT", "FD"):
                pi0_all = optimizer_events["pi0"][region][bin_key]
                dvcs_all = optimizer_events["dvcs"][region][bin_key]
                data_all = optimizer_events["data"][region][bin_key]

                mx2_idx = FEATURE_INDEX["Mx2_1"]

                pi0_pre = (
                    np.isfinite(pi0_all[:, mx2_idx])
                    & (
                        pi0_all[:, mx2_idx]
                        < optimizer_mx2_1_preselection_max
                    )
                )
                dvcs_pre = (
                    np.isfinite(dvcs_all[:, mx2_idx])
                    & (
                        dvcs_all[:, mx2_idx]
                        < optimizer_mx2_1_preselection_max
                    )
                )
                sideband_mask = (
                    np.isfinite(data_all[:, mx2_idx])
                    & (
                        data_all[:, mx2_idx]
                        >= optimizer_data_sideband_min
                    )
                    & (
                        data_all[:, mx2_idx]
                        < optimizer_data_sideband_max
                    )
                )

                pi0_events = pi0_all[pi0_pre]
                dvcs_events = dvcs_all[dvcs_pre]
                data_sideband_events = data_all[sideband_mask]

                log(
                    f"{period.label} {region}, "
                    f"{energy_min:g}<=E_probe<{energy_max:g} GeV: "
                    f"after Mx2_1<"
                    f"{optimizer_mx2_1_preselection_max:.3g} GeV^2 -> "
                    f"AAO={pi0_events.shape[0]:,}, "
                    f"DVCSgen={dvcs_events.shape[0]:,}; "
                    f"data sideband "
                    f"{optimizer_data_sideband_min:.3g}<=Mx2_1<"
                    f"{optimizer_data_sideband_max:.3g} GeV^2 -> "
                    f"{data_sideband_events.shape[0]:,}."
                )

                results = optimize_rectangular_cuts(
                    pi0_events,
                    dvcs_events,
                    data_sideband_events,
                    max_steps=optimizer_steps,
                    n_quantiles=optimizer_quantiles,
                    min_step_pi0_eff=optimizer_min_step_pi0_eff,
                    min_total_pi0_eff=optimizer_min_total_pi0_eff,
                    min_improvement=optimizer_min_improvement,
                    data_sideband_weight=optimizer_data_sideband_weight,
                )

                write_optimizer_outputs(
                    period,
                    region,
                    bin_key,
                    energy_min,
                    energy_max,
                    pi0_events,
                    dvcs_events,
                    data_sideband_events,
                    results,
                    optimizer_outdir,
                )

                closure_results = run_template_fit_closure(
                    pi0_events,
                    dvcs_events,
                    results,
                    n_toys=closure_toys,
                    n_pseudodata_events=closure_events,
                    n_bins=closure_bins,
                    seed=(
                        700_000
                        + 10_000 * period_index
                        + 1000 * [key for key, _lo, _hi in E_PROBE_BINS].index(
                            bin_key
                        )
                        + (0 if region == "FT" else 100)
                    ),
                )

                closure_outdir = (
                    output_dir.parent / "template_fit_closure"
                )
                write_template_fit_closure_outputs(
                    period,
                    region,
                    bin_key,
                    energy_min,
                    energy_max,
                    closure_results,
                    closure_outdir,
                )

                optimizer_summary["energy_bins"][bin_key]["regions"][region] = {
                    "n_pi0_baseline": int(pi0_events.shape[0]),
                    "n_dvcs_baseline": int(dvcs_events.shape[0]),
                    "n_data_sideband_baseline": int(
                        data_sideband_events.shape[0]
                    ),
                    "results": results,
                    "closure_results": closure_results,
                }
            #endfor
        #endfor

        combined_optimizer_plot = make_combined_optimizer_canvas(
            period,
            optimizer_summary,
            optimizer_outdir,
        )
        log(
            f"{period.label}: wrote E_probe-binned FD+FT optimizer canvas "
            f"{combined_optimizer_plot}"
        )

        combined_closure_plot = make_combined_closure_canvas(
            period,
            optimizer_summary,
            output_dir.parent / "template_fit_closure",
        )
        log(
            f"{period.label}: wrote combined template-fit closure canvas "
            f"{combined_closure_plot}"
        )

        # Make one cumulative shape canvas for every period x detector x
        # E_probe bin, using exactly the cuts selected in that bin.
        for bin_key, energy_min, energy_max in E_PROBE_BINS:
            for region in ("FT", "FD"):
                region_results = (
                    optimizer_summary["energy_bins"][bin_key]
                    ["regions"][region]["results"]
                )
                shape_out = make_optimizer_cut_shape_canvas(
                    period,
                    region,
                    bin_key,
                    energy_min,
                    energy_max,
                    optimizer_events,
                    region_results,
                    output_dir,
                    linear_y,
                    optimizer_mx2_1_preselection_max,
                )
                log(
                    f"{period.label}: wrote {region} "
                    f"{energy_min:g}-{energy_max:g} GeV "
                    f"post-optimizer shape canvas {shape_out}"
                )
            #endfor
        #endfor
    else:
        # If optimization is explicitly disabled, preserve the old fixed
        # three-row unbinned concept canvases.
        for region in ("FT", "FD"):
            out = make_canvas(
                period,
                region,
                sample_hists,
                sample_counts,
                output_dir,
                linear_y,
            )
            log(f"{period.label}: wrote {out}")
        #endfor
    #endif

    log(
        f"{period.label}: complete in "
        f"{time.perf_counter() - t0:.1f} s."
    )

    return optimizer_summary


def process_period_worker(
    period: Period,
    tree_name: str,
    max_entries: int,
    step_size: int,
    output_dir_str: str,
    linear_y: bool,
    run_optimizer: bool,
    optimizer_max_events: int,
    optimizer_steps: int,
    optimizer_quantiles: int,
    optimizer_min_step_pi0_eff: float,
    optimizer_min_total_pi0_eff: float,
    optimizer_min_improvement: float,
    optimizer_mx2_1_preselection_max: float,
    optimizer_data_sideband_min: float,
    optimizer_data_sideband_max: float,
    optimizer_data_sideband_weight: float,
    closure_toys: int,
    closure_events: int,
    closure_bins: int,
) -> str:
    """
    Picklable process-level wrapper for one run period.

    Each worker streams that period's data, AAO, and DVCSgen files
    sequentially. FT and FD are filled simultaneously within each file,
    avoiding nested ROOT I/O parallelism.
    """
    return process_period(
        period,
        tree_name,
        max_entries,
        step_size,
        Path(output_dir_str),
        linear_y,
        run_optimizer,
        optimizer_max_events,
        optimizer_steps,
        optimizer_quantiles,
        optimizer_min_step_pi0_eff,
        optimizer_min_total_pi0_eff,
        optimizer_min_improvement,
        optimizer_mx2_1_preselection_max,
        optimizer_data_sideband_min,
        optimizer_data_sideband_max,
        optimizer_data_sideband_weight,
        closure_toys,
        closure_events,
        closure_bins,
    )



def print_optimizer_summary(
    summaries_by_period: Dict[str, dict],
    selected_periods: Sequence[Period],
) -> None:
    """
    Print all optimizer results once, serially, after parallel processing ends.
    Results are grouped by period, then E_probe bin, then detector region.
    """
    print("", flush=True)
    print("=" * 100, flush=True)
    print("E_PROBE-BINNED RECTANGULAR CUT OPTIMIZER SUMMARY", flush=True)
    print("=" * 100, flush=True)

    for period in selected_periods:
        summary = summaries_by_period.get(period.key)
        if summary is None:
            continue
        #endif

        print(f"\n{period.label}", flush=True)
        print("-" * len(period.label), flush=True)

        energy_summaries = summary.get("energy_bins", {})
        if not energy_summaries:
            print(
                "  Optimizer disabled or no optimizer result.",
                flush=True,
            )
            continue
        #endif

        for bin_key, energy_min, energy_max in E_PROBE_BINS:
            bin_summary = energy_summaries.get(bin_key)
            if bin_summary is None:
                continue
            #endif

            print(
                f"\n  {energy_min:g} <= E_probe < {energy_max:g} GeV",
                flush=True,
            )

            for region in ("FT", "FD"):
                region_summary = bin_summary.get("regions", {}).get(region)
                if region_summary is None:
                    continue
                #endif

                print(
                    f"    {region}: baseline "
                    f"AAO={region_summary['n_pi0_baseline']:,}, "
                    f"DVCSgen={region_summary['n_dvcs_baseline']:,}, "
                    f"data-SB="
                    f"{region_summary['n_data_sideband_baseline']:,}",
                    flush=True,
                )

                results = region_summary.get("results", [])
                if not results:
                    print(
                        "      No cut passed the configured "
                        "retention/improvement constraints.",
                        flush=True,
                    )
                    continue
                #endif

                for result in results:
                    print(
                        f"      {result['step']}. "
                        f"{cut_expression(result)}"
                        f"  | step AAO="
                        f"{100.0*result['step_pi0_eff']:.2f}%"
                        f", DVCS="
                        f"{100.0*result['step_dvcs_eff']:.2f}%"
                        f", data-SB="
                        f"{100.0*result['step_data_sideband_eff']:.2f}%"
                        f"  | cumulative AAO="
                        f"{100.0*result['total_pi0_eff']:.2f}%"
                        f", DVCS="
                        f"{100.0*result['total_dvcs_eff']:.2f}%"
                        f", data-SB="
                        f"{100.0*result['total_data_sideband_eff']:.2f}%",
                        flush=True,
                    )
                #endfor

                final = results[-1]
                print(
                    f"      FINAL: AAO="
                    f"{100.0*final['total_pi0_eff']:.2f}%, "
                    f"DVCSgen="
                    f"{100.0*final['total_dvcs_eff']:.2f}%, "
                    f"data-SB="
                    f"{100.0*final['total_data_sideband_eff']:.2f}%, "
                    f"F={final['score']:.4g}",
                    flush=True,
                )

                closure_results = region_summary.get(
                    "closure_results",
                    [],
                )
                final_closure = [
                    result
                    for result in closure_results
                    if int(result["step"]) == len(results)
                    and abs(result["true_fraction"] - 0.50) < 1.0e-9
                ]
                if final_closure:
                    closure = final_closure[0]
                    print(
                        f"      CLOSURE @ final step, f_pi0=0.50: "
                        f"best={closure['feature']}, "
                        f"JS={closure['js_divergence']:.4f}, "
                        f"<f_fit>={closure['mean_fit']:.4f}, "
                        f"bias={closure['bias']:+.4f}, "
                        f"RMS={closure['rms_error']:.4f}",
                        flush=True,
                    )
                #endif
            #endfor
        #endfor
    #endfor

    print("\n" + "=" * 100, flush=True)


def main() -> int:
    args = parse_args()

    if args.max_entries < 0:
        raise ValueError("--max-entries must be >= 0.")
    #endif
    if args.step_size <= 0:
        raise ValueError("--step-size must be > 0.")
    #endif
    if args.workers <= 0:
        raise ValueError("--workers must be > 0.")
    #endif
    if args.optimizer_max_events <= 0:
        raise ValueError("--optimizer-max-events must be > 0.")
    #endif
    if args.optimizer_steps < 0:
        raise ValueError("--optimizer-steps must be >= 0.")
    #endif
    if args.optimizer_quantiles < 3:
        raise ValueError("--optimizer-quantiles must be >= 3.")
    #endif
    if not 0.0 < args.optimizer_min_step_pi0_eff <= 1.0:
        raise ValueError("--optimizer-min-step-pi0-eff must be in (0,1].")
    #endif
    if not 0.0 < args.optimizer_min_total_pi0_eff <= 1.0:
        raise ValueError("--optimizer-min-total-pi0-eff must be in (0,1].")
    #endif
    if args.optimizer_min_improvement < 1.0:
        raise ValueError("--optimizer-min-improvement must be >= 1.")
    #endif
    if (
        args.optimizer_data_sideband_min
        >= args.optimizer_data_sideband_max
    ):
        raise ValueError(
            "--optimizer-data-sideband-min must be < "
            "--optimizer-data-sideband-max."
        )
    #endif
    if (
        args.optimizer_data_sideband_max
        > args.optimizer_mx2_1_preselection_max
    ):
        raise ValueError(
            "--optimizer-data-sideband-max must be <= "
            "--optimizer-mx2-1-preselection-max."
        )
    #endif
    if args.optimizer_data_sideband_weight < 0.0:
        raise ValueError("--optimizer-data-sideband-weight must be >= 0.")
    #endif
    if args.closure_toys <= 0:
        raise ValueError("--closure-toys must be > 0.")
    #endif
    if args.closure_events <= 0:
        raise ValueError("--closure-events must be > 0.")
    #endif
    if args.closure_bins < 4:
        raise ValueError("--closure-bins must be >= 4.")
    #endif

    selected_keys = set(args.period or [p.key for p in PERIODS])
    selected_periods = [p for p in PERIODS if p.key in selected_keys]

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    log(
        "Standalone photon-efficiency concept study: "
        "data + AAO pi0 + DVCSgen only; CLASDIS disabled."
    )
    log(
        "Shape-comparison y-axis: "
        + ("linear (--linear override)." if args.linear else "logarithmic (default).")
    )
    if args.no_optimize_cuts:
        log("Rectangular-cut optimizer: disabled.")
    else:
        log(
            "Rectangular-cut optimizer: enabled independently for every "
            "period and separately in FT/FD; "
            f"reservoir <= {args.optimizer_max_events:,} events/class/region, "
            f"max steps={args.optimizer_steps}; "
            f"loose Mx2_1<{args.optimizer_mx2_1_preselection_max:.3g} GeV^2; "
            f"data sideband {args.optimizer_data_sideband_min:.3g}<=Mx2_1<"
            f"{args.optimizer_data_sideband_max:.3g} GeV^2; "
            f"sideband weight={args.optimizer_data_sideband_weight:.3g}; "
            "E_probe (= Emiss2) bins: "
            "0.4-1, 1-2, 2-4, 4-9.5 GeV; "
            "Emiss2 excluded from optimizer thresholds; "
            f"closure={args.closure_toys} toys x "
            f"{args.closure_events} events at f_pi0=0.2,0.5,0.8."
        )
    #endif
    log(
        "Common pre-binning selection: "
        f"Angle(e,gamma) > {THETA_EGAMMA_MIN_DEG:.1f} deg. "
        f"Photon regions: FT {FT_THETA_MIN_DEG:.1f}-{FT_THETA_MAX_DEG:.1f} deg, "
        f"FD {FD_THETA_MIN_DEG:.1f}-{FD_THETA_MAX_DEG:.1f} deg."
    )

    n_workers = min(int(args.workers), len(selected_periods))
    log(
        f"Period-level parallelism: {n_workers} process(es) for "
        f"{len(selected_periods)} selected period(s)."
    )

    summaries_by_period: Dict[str, dict] = {}

    if n_workers == 1:
        for period in selected_periods:
            summaries_by_period[period.key] = process_period(
                period,
                args.tree,
                args.max_entries,
                args.step_size,
                output_dir,
                args.linear,
                not args.no_optimize_cuts,
                args.optimizer_max_events,
                args.optimizer_steps,
                args.optimizer_quantiles,
                args.optimizer_min_step_pi0_eff,
                args.optimizer_min_total_pi0_eff,
                args.optimizer_min_improvement,
                args.optimizer_mx2_1_preselection_max,
                args.optimizer_data_sideband_min,
                args.optimizer_data_sideband_max,
                args.optimizer_data_sideband_weight,
                args.closure_toys,
                args.closure_events,
                args.closure_bins,
            )
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_period = {
                executor.submit(
                    process_period_worker,
                    period,
                    args.tree,
                    args.max_entries,
                    args.step_size,
                    str(output_dir),
                    args.linear,
                    not args.no_optimize_cuts,
                    args.optimizer_max_events,
                    args.optimizer_steps,
                    args.optimizer_quantiles,
                    args.optimizer_min_step_pi0_eff,
                    args.optimizer_min_total_pi0_eff,
                    args.optimizer_min_improvement,
                    args.optimizer_mx2_1_preselection_max,
                    args.optimizer_data_sideband_min,
                    args.optimizer_data_sideband_max,
                    args.optimizer_data_sideband_weight,
                    args.closure_toys,
                    args.closure_events,
                    args.closure_bins,
                ): period
                for period in selected_periods
            }

            for future in as_completed(future_to_period):
                period = future_to_period[future]
                try:
                    summaries_by_period[period.key] = future.result()
                except Exception as exc:
                    raise RuntimeError(
                        f"Parallel concept-study processing failed for "
                        f"{period.label}: {exc}"
                    ) from exc
                #endtry
            #endfor
        #endwith
    #endif

    if not args.no_optimize_cuts:
        print_optimizer_summary(
            summaries_by_period,
            selected_periods,
        )
    #endif

    log(f"Done. Outputs are in {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
