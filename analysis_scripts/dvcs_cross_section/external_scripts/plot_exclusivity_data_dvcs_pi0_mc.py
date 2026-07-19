#!/usr/bin/env python3
"""
plot_exclusivity_data_dvcs_pi0_mc.py

For each of the five run periods and three detector topologies, compare the
eight exclusivity/kinematic-variable shapes and perform a shared two-template fit across all eight fitted projections by
default, with optional configurable driver subsets:

    data = (1-f_pi0) * shifted-and-smeared DVCS MC
         + f_pi0     * eppi0 MC reconstructed as DVCS.

The detector topology plus the minimal global preselection are required:

    (-t1) < 1.0
    open_angle_ep2 > 5 degrees
    all reconstructed FD particles occupy distinct FD sectors

No hard exclusivity cuts are applied. The DVCS nuisance parameters are fitted
inside MC-defined signal regions by default: 90% for nuisance/core fits and
95% for the shared DVCS pi0-fraction fit. The fit uses raw binned data counts and a Poisson likelihood.
Each projection is normalized to the observed data yield inside the fit mask
used by that projection, so the extrapolation outside the mask is a genuine
shape validation:

    f_pi0        : one fraction shared by all selected fit histograms
    shift        : variable-specific DVCS-template shift or log-scale shift
    sigma_add    : variable-specific additive or log-space smearing

The eppi0 template is fixed. Signed variables use additive shift plus Gaussian
smearing. Positive-definite pTmiss and theta_gamma_gamma use log-space shift
and smearing. By default Delta_phi, theta_gamma_gamma and pTmiss determine the shared fraction using one common in-range event population. When a
restricted --fraction-variable subset is supplied, the remaining variables
are validation projections with profiled nuisance parameters. Optional Gaussian nuisance penalties discourage extreme template
shifts and broadenings.

The script also compares reconstructed eppi0 data directly with reconstructed eppi0 MC using theta_pi0_pi0 in place of theta_gamma_gamma. The nuisance fits use xF, xF2, z2 and theta in place of the missing-mass projections. The shape-comparison canvases retain those new variables and additionally show Mx2, Mx2_1, Mx2_2 and pT in a third row. The direct eppi0 nuisance fit is restricted to an MC-defined signal core so that out-of-core backgrounds and tails remain diagnostics rather than forcing excessive template morphing.

Outputs:

  * one unit-area 3x4 shape-comparison canvas for each period/topology combination
  * one 4x2 DVCS-core two-template-fit canvas for each period/topology combination
  * fit_results.csv containing all fitted parameters and diagnostics
  * compact 2x4 iterative cut-development canvases
  * compact 2x4 post-cut summary canvases
  * main-suite-compatible combined_cuts JSON files for 97%, 99% and 95%

The five run periods are processed in parallel with at most five worker processes. Dependencies: Python 3, numpy, matplotlib, scipy and either uproot or PyROOT.
"""

import argparse
import concurrent.futures
import multiprocessing
import csv
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

# Configure Matplotlib before importing pyplot. This avoids slow or locked font
# caches on ifarm-like systems.
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.optimize import minimize, minimize_scalar
    from scipy.special import ndtr
except ImportError:
    gaussian_filter1d = None
    minimize = None
    minimize_scalar = None
    ndtr = None
#endif

try:
    import uproot
except ImportError:
    uproot = None
#endif

try:
    import ROOT
except ImportError:
    ROOT = None
#endif


TREE_NAME = "PhysicsEvents"
DEFAULT_OUTPUT_DIR = "output/exclusivity_data_dvcs_pi0_mc"
DEFAULT_STEP_SIZE = "250 MB"

T1_ABS_MAX = 1.0
OPEN_ANGLE_MIN_DEG = 5.0
REQUIRE_DISTINCT_FD_SECTORS = True


@dataclass(frozen=True)
class PeriodConfig:
    key: str
    label: str
    data_file: str
    dvcs_mc_file: str
    pi0_as_dvcs_mc_file: str
    eppi0_data_file: str
    eppi0_mc_file: str


@dataclass(frozen=True)
class TopologyConfig:
    key: str
    label: str
    detector1: int
    detector2: int


@dataclass(frozen=True)
class VariableConfig:
    branch: str
    label: str
    bins: int
    xmin: float
    xmax: float
    aliases: Tuple[str, ...] = ()


@dataclass
class FitResult:
    success: bool
    message: str
    f_pi0: float = math.nan
    f_pi0_err: float = math.nan
    shift: float = math.nan
    shift_err: float = math.nan
    sigma_add: float = math.nan
    sigma_add_err: float = math.nan
    sigma_right: float = math.nan
    sigma_right_err: float = math.nan
    deviance: float = math.nan
    ndf: int = 0
    data_total: float = 0.0
    model_counts: Optional[np.ndarray] = None
    dvcs_component_counts: Optional[np.ndarray] = None
    pi0_component_counts: Optional[np.ndarray] = None
    transformed_dvcs_shape: Optional[np.ndarray] = None
    fit_mask: Optional[np.ndarray] = None
    morph_label: str = "additive"
    excluded_data_counts: float = 0.0
    excluded_model_counts: float = 0.0
    excluded_excess_counts: float = 0.0
    excluded_excess_fraction: float = 0.0


@dataclass
class SharedFitSummary:
    success: bool
    message: str
    f_pi0: float = math.nan
    f_pi0_err: float = math.nan
    deviance: float = math.nan
    ndf: int = 0
    variable_results: Optional[Dict[str, FitResult]] = None
    fraction_variables: Tuple[str, ...] = ()


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        key="fa18_inb",
        label="Fa18 Inb",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_50nA_10604MeV_epgamma.root",
        eppi0_data_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root",
        eppi0_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_50nA_10604MeV.root",
    ),
    PeriodConfig(
        key="fa18_out",
        label="Fa18 Out",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_50nA_10604MeV_epgamma.root",
        eppi0_data_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root",
        eppi0_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_50nA_10604MeV.root",
    ),
    PeriodConfig(
        key="sp19_inb",
        label="Sp19 Inb",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_50nA_10200MeV_epgamma.root",
        eppi0_data_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root",
        eppi0_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_50nA_10200MeV.root",
    ),
    PeriodConfig(
        key="sp18_inb",
        label="Sp18 Inb",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_50nA_10594MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_inb_50nA_10594MeV_epgamma.root",
        eppi0_data_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root",
        eppi0_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_50nA_10594MeV.root",
    ),
    PeriodConfig(
        key="sp18_out",
        label="Sp18 Out",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_45nA_10594MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_out_45nA_10594MeV_epgamma.root",
        eppi0_data_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root",
        eppi0_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_45nA_10594MeV.root",
    ),
)


TOPOLOGIES: Tuple[TopologyConfig, ...] = (
    TopologyConfig("FD_FD", "(FD, FD)", detector1=1, detector2=1),
    TopologyConfig("CD_FD", "(CD, FD)", detector1=2, detector2=1),
    TopologyConfig("CD_FT", "(CD, FT)", detector1=2, detector2=0),
)


VARIABLES: Tuple[VariableConfig, ...] = (
    VariableConfig("Delta_phi", r"$\Delta\phi$ (rad)", 100, 2.84159, 3.44159),
    VariableConfig("theta_gamma_gamma", r"$\theta_{\gamma\gamma}$ (rad)", 120, 0.0, 3.0),
    VariableConfig("pTmiss", r"$p_{T}^{\mathrm{miss}}$ (GeV)", 125, 0.0, 0.5),
    VariableConfig("z2", r"$z_2$", 100, 0.0, 1.1),
    VariableConfig("Emiss2", r"$E_{\mathrm{miss}}^{2}$ (GeV$^2$)", 100, -1.0, 2.0),
    VariableConfig("xF", r"$x_F(ep+\gamma)$", 100, -0.5, 0.2),
    VariableConfig(
        "xF2",
        r"$x_{F,2}(\gamma)$",
        100,
        0.0,
        1.0,
        aliases=("xF_2",),
    ),
    VariableConfig(
        "theta",
        r"$\theta_{p\gamma}^{\mathrm{CM}}$ (rad)",
        100,
        2.0,
        math.pi,
        aliases=("theta2", "theta_2"),
    ),
)


PI0_VARIABLES: Tuple[VariableConfig, ...] = (
    VariableConfig("Delta_phi", r"$\Delta\phi$ (rad)", 100, 2.84159, 3.44159),
    VariableConfig("theta_pi0_pi0", r"$\theta_{\pi^0\pi^0}$ (rad)", 120, 0.0, 3.0),
    VariableConfig("pTmiss", r"$p_{T}^{\mathrm{miss}}$ (GeV)", 125, 0.0, 0.5),
    VariableConfig("z2", r"$z_2$", 100, 0.0, 1.1),
    VariableConfig("Emiss2", r"$E_{\mathrm{miss}}^{2}$ (GeV$^2$)", 100, -1.0, 2.0),
    VariableConfig("xF", r"$x_F(ep+\pi^0)$", 100, -0.5, 0.2),
    VariableConfig(
        "xF2",
        r"$x_{F,2}(\pi^0)$",
        100,
        0.0,
        1.0,
        aliases=("xF_2",),
    ),
    VariableConfig(
        "theta",
        r"$\theta_{p\pi^0}^{\mathrm{CM}}$ (rad)",
        100,
        2.0,
        math.pi,
        aliases=("theta2", "theta_2"),
    ),
)


# All eight variables remain available for nuisance validation and automatic cut selection, while only the configured fraction drivers determine f_pi0.
# The unit-area shape-comparison canvases additionally include the three
# missing-mass projections and the event transverse momentum pT.
SHAPE_ONLY_VARIABLES: Tuple[VariableConfig, ...] = (
    VariableConfig(
        "Mx2",
        r"$M_x^2$ (GeV$^2$)",
        100,
        -0.03,
        0.03,
        aliases=("Mx2_epg", "Mx2_eg", "Mx2_epgamma", "Mx2_epi0", "Mx2_eppi0"),
    ),
    VariableConfig(
        "Mx2_1",
        r"$M_{x1}^2$ (GeV$^2$)",
        100,
        -1.5,
        1.5,
        aliases=("Mx2_ep", "Mx2_x1", "Mx2_proton", "Mx2_p"),
    ),
    VariableConfig(
        "Mx2_2",
        r"$M_{x2}^2$ (GeV$^2$)",
        125,
        -1.0,
        4.0,
        aliases=("Mx2_egamma", "Mx2_gamma", "Mx2_pi0", "Mx2_x2"),
    ),
    VariableConfig("pT", r"$p_T$ (GeV)", 100, 0.0, 0.3),
)

DVCS_SHAPE_VARIABLES: Tuple[VariableConfig, ...] = VARIABLES + SHAPE_ONLY_VARIABLES
PI0_SHAPE_VARIABLES: Tuple[VariableConfig, ...] = PI0_VARIABLES + SHAPE_ONLY_VARIABLES


SAMPLE_LABELS = {
    "data": r"$e'p'\gamma$ data",
    "dvcs_mc": "DVCS MC",
    "pi0_mc": r"$e\pi^0$ MC as DVCS",
}


SAMPLE_COLORS = {
    "data": "black",
    "dvcs_mc": "tab:blue",
    "pi0_mc": "tab:red",
    "fit": "tab:green",
    "raw_dvcs": "0.55",
}


def log(message: str) -> None:
    now = time.strftime("%H:%M:%S")
    print(f"[{now}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot unit-area exclusivity shapes and two-template fits for DVCS data, reconstructed DVCS MC and eppi0 MC reconstructed as DVCS after the minimal global preselection."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--period",
        choices=[p.key for p in PERIODS],
        action="append",
        help="Restrict to one or more periods. May be supplied repeatedly.",
    )
    parser.add_argument(
        "--topology",
        choices=[t.key for t in TOPOLOGIES],
        action="append",
        help="Restrict to one or more topologies. May be supplied repeatedly.",
    )
    parser.add_argument(
        "--step-size",
        default=DEFAULT_STEP_SIZE,
        help=f"uproot iteration chunk size (default: {DEFAULT_STEP_SIZE})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=5,
        help=(
            "Number of run-period worker processes. The value is hard-capped "
            "at 5 (default: 5)."
        ),
    )
    parser.add_argument(
        "--max-shift-bins",
        type=float,
        default=12.0,
        help="Maximum absolute DVCS-template shift in histogram bins (default: 12).",
    )
    parser.add_argument(
        "--max-smear-bins",
        type=float,
        default=20.0,
        help="Maximum added Gaussian smearing in histogram bins (default: 20).",
    )
    parser.add_argument(
        "--fit-min-counts",
        type=int,
        default=100,
        help="Minimum data entries inside a variable range required for a fit (default: 100).",
    )
    parser.add_argument(
        "--dvcs-core-containment",
        type=float,
        default=0.90,
        help=(
            "DVCS-MC containment used to determine variable-specific shift and "
            "smearing nuisances (default: 0.90)."
        ),
    )
    parser.add_argument(
        "--dvcs-fraction-containment",
        type=float,
        default=0.95,
        help=(
            "DVCS-MC containment used by the selected histograms to determine "
            "the shared pi0 fraction (default: 0.95)."
        ),
    )
    parser.add_argument(
        "--pi0-core-containment",
        type=float,
        default=0.90,
        help=(
            "MC containment used by the direct eppi0 data--MC core fits "
            "(default: 0.90)."
        ),
    )
    parser.add_argument(
        "--cut-nominal-containment",
        type=float,
        default=0.95,
        help="Nominal iterative signal containment (default: 0.95).",
    )
    parser.add_argument(
        "--cut-loose-containment",
        type=float,
        default=0.99,
        help="Loose iterative signal containment (default: 0.99).",
    )
    parser.add_argument(
        "--cut-tight-containment",
        type=float,
        default=0.90,
        help="Tight iterative signal containment (default: 0.90).",
    )
    parser.add_argument(
        "--skip-iterative-cuts",
        action="store_true",
        help="Skip iterative exclusivity-cut development and JSON generation.",
    )
    parser.add_argument(
        "--fraction-variable",
        choices=[v.branch for v in VARIABLES],
        action="append",
        help=(
            "Variable used to determine the shared pi0 fraction. May be supplied "
            "repeatedly. Default: Delta_phi, theta_gamma_gamma and pTmiss."
        ),
    )
    parser.add_argument(
        "--shift-prior-bins",
        type=float,
        default=4.0,
        help=(
            "Gaussian prior width for additive template shifts, in histogram bins "
            "(default: 4). Positive-definite variables use a log-shift prior of 0.20."
        ),
    )
    parser.add_argument(
        "--smear-prior-bins",
        type=float,
        default=8.0,
        help=(
            "Half-Gaussian prior width for additive smearing, in histogram bins "
            "(default: 8). Positive-definite variables use a log-smearing prior of 0.40."
        ),
    )
    parser.add_argument(
        "--disable-nuisance-penalties",
        action="store_true",
        help="Disable Gaussian penalties on DVCS-template shift and smearing nuisances.",
    )
    parser.add_argument(
        "--outside-overshoot-penalty",
        type=float,
        default=0.25,
        help=(
            "Weight of the one-sided penalty for model counts exceeding data "
            "outside the active fit region (default: 0.25; 0 disables it)."
        ),
    )
    parser.add_argument(
        "--emiss2-mean-order-penalty",
        type=float,
        default=25.0,
        help=(
            "Weight of the E_miss^2 physical-ordering penalty. It discourages "
            "the morphed DVCS mean from lying farther from zero than the pi0 "
            "background mean (default: 25; 0 disables it)."
        ),
    )
    parser.add_argument(
        "--log-y",
        action="store_true",
        help="Use logarithmic y axes.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="PNG resolution in dots per inch (default: 180).",
    )
    parser.add_argument(
        "--no-clean-output",
        action="store_true",
        help=(
            "Keep previous iterative-cut outputs. By default they are removed "
            "before running so aborted jobs cannot leave stale canvases."
        ),
    )
    return parser.parse_args()


def selected_periods(keys: Optional[Sequence[str]]) -> List[PeriodConfig]:
    if not keys:
        return list(PERIODS)
    # endif

    wanted = set(keys)
    return [period for period in PERIODS if period.key in wanted]


def selected_topologies(keys: Optional[Sequence[str]]) -> List[TopologyConfig]:
    if not keys:
        return list(TOPOLOGIES)
    # endif

    wanted = set(keys)
    return [topology for topology in TOPOLOGIES if topology.key in wanted]


def require_input_file(path: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input ROOT file does not exist: {path}")
    # endif


def io_backend() -> str:
    if uproot is not None:
        return "uproot"
    # endif

    if ROOT is not None:
        return "pyroot"
    # endif

    raise RuntimeError(
        "Neither uproot nor PyROOT is available. Run inside the CLAS12 ROOT "
        "environment or install uproot."
    )


def available_tree_branches(path: str) -> set[str]:
    backend = io_backend()

    if backend == "uproot":
        with uproot.open(path) as root_file:
            if TREE_NAME not in root_file:
                raise KeyError(f"Tree '{TREE_NAME}' not found in {path}")
            # endif
            return set(root_file[TREE_NAME].keys())
        # endwith
    # endif

    root_file = ROOT.TFile.Open(path, "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    # endif

    tree = root_file.Get(TREE_NAME)
    if tree is None:
        root_file.Close()
        raise KeyError(f"Tree '{TREE_NAME}' not found in {path}")
    # endif

    branches = {branch.GetName() for branch in tree.GetListOfBranches()}
    root_file.Close()
    return branches


def resolve_variable_branches(path: str, variables: Sequence[VariableConfig] = VARIABLES) -> Dict[str, str]:
    branches = available_tree_branches(path)

    required_selection = {
        "detector1", "detector2", "t1", "open_angle_ep2",
        "e_phi", "p1_phi", "p2_phi",
    }
    missing_selection = sorted(required_selection - branches)
    if missing_selection:
        raise KeyError(
            f"Missing selection branches in {path}: {', '.join(missing_selection)}"
        )
    # endif

    resolved: Dict[str, str] = {}
    missing_variables: List[str] = []

    for variable in variables:
        candidates = (variable.branch,) + variable.aliases
        match = next((name for name in candidates if name in branches), None)

        if match is None:
            missing_variables.append(
                f"{variable.branch} (tried: {', '.join(candidates)})"
            )
        else:
            resolved[variable.branch] = match
        # endif
    # endfor

    if missing_variables:
        raise KeyError(
            f"Missing exclusivity branches in {path}:\n  "
            + "\n  ".join(missing_variables)
        )
    # endif

    return resolved


def fd_sector_from_phi(phi_rad: np.ndarray) -> np.ndarray:
    """Return the CLAS12 FD sector (1--6) from an azimuth in radians."""

    phi_deg = np.mod(np.degrees(phi_rad), 360.0)
    sector = np.full(phi_deg.shape, -1, dtype=np.int16)
    sector[(phi_deg >= 330.0) | (phi_deg < 30.0)] = 1
    sector[(phi_deg >= 30.0) & (phi_deg < 90.0)] = 2
    sector[(phi_deg >= 90.0) & (phi_deg < 150.0)] = 3
    sector[(phi_deg >= 150.0) & (phi_deg < 210.0)] = 4
    sector[(phi_deg >= 210.0) & (phi_deg < 270.0)] = 5
    sector[(phi_deg >= 270.0) & (phi_deg < 330.0)] = 6
    return sector


def empty_histograms(variables: Sequence[VariableConfig] = VARIABLES) -> Dict[str, np.ndarray]:
    return {
        variable.branch: np.zeros(variable.bins, dtype=np.float64)
        for variable in variables
    }


def fill_histograms_uproot(
    path: str,
    topologies: Sequence[TopologyConfig],
    step_size: str,
    apply_mx2_1_cut: bool = False,
    variables: Sequence[VariableConfig] = VARIABLES,
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, int]]:
    """Read one ROOT file once and fill all requested topology/variable histograms."""

    require_input_file(path)
    resolved = resolve_variable_branches(path, variables)

    expressions = [
        "detector1", "detector2", "t1", "open_angle_ep2",
        "e_phi", "p1_phi", "p2_phi",
    ] + sorted(set(resolved.values()))
    histograms = {topology.key: empty_histograms(variables) for topology in topologies}
    selected_events = {topology.key: 0 for topology in topologies}

    tree_spec = f"{path}:{TREE_NAME}"
    log(f"Reading {tree_spec}")

    for arrays in uproot.iterate(
        tree_spec,
        expressions=expressions,
        step_size=step_size,
        library="np",
    ):
        detector1 = np.asarray(arrays["detector1"])
        detector2 = np.asarray(arrays["detector2"])
        t1 = np.asarray(arrays["t1"], dtype=np.float64)
        open_angle = np.asarray(arrays["open_angle_ep2"], dtype=np.float64)

        base_mask = (
            np.isfinite(t1)
            & np.isfinite(open_angle)
            & ((-t1) < T1_ABS_MAX)
            & (open_angle > OPEN_ANGLE_MIN_DEG)
        )

        if REQUIRE_DISTINCT_FD_SECTORS:
            e_sector = fd_sector_from_phi(np.asarray(arrays["e_phi"], dtype=np.float64))
            p_sector = fd_sector_from_phi(np.asarray(arrays["p1_phi"], dtype=np.float64))
            g_sector = fd_sector_from_phi(np.asarray(arrays["p2_phi"], dtype=np.float64))
        # endif

        for topology in topologies:
            mask = (
                base_mask
                & (detector1 == topology.detector1)
                & (detector2 == topology.detector2)
            )

            if REQUIRE_DISTINCT_FD_SECTORS:
                mask &= e_sector >= 1
                if topology.detector1 == 1:
                    mask &= (p_sector >= 1) & (p_sector != e_sector)
                # endif
                if topology.detector2 == 1:
                    mask &= (g_sector >= 1) & (g_sector != e_sector)
                # endif
                if topology.detector1 == 1 and topology.detector2 == 1:
                    mask &= p_sector != g_sector
                # endif
            # endif

            if apply_mx2_1_cut:
                mx2_1_values = np.asarray(arrays[resolved["Mx2_1"]], dtype=np.float64)
                mask &= np.isfinite(mx2_1_values) & (mx2_1_values < mx2_1_upper_cut(topology))
            # endif

            selected_events[topology.key] += int(np.count_nonzero(mask))

            if not np.any(mask):
                continue
            # endif

            for variable in variables:
                source_branch = resolved[variable.branch]
                values = np.asarray(arrays[source_branch], dtype=np.float64)[mask]
                values = values[np.isfinite(values)]

                counts, _ = np.histogram(
                    values,
                    bins=variable.bins,
                    range=(variable.xmin, variable.xmax),
                )
                histograms[topology.key][variable.branch] += counts
            # endfor
        # endfor
    # endfor

    return histograms, selected_events


def fill_histograms_pyroot(
    path: str,
    topologies: Sequence[TopologyConfig],
    step_size: str,
    apply_mx2_1_cut: bool = False,
    variables: Sequence[VariableConfig] = VARIABLES,
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, int]]:
    del step_size  # PyROOT performs a direct TTree event loop.

    require_input_file(path)
    resolved = resolve_variable_branches(path, variables)
    histograms = {topology.key: empty_histograms(variables) for topology in topologies}
    selected_events = {topology.key: 0 for topology in topologies}

    root_file = ROOT.TFile.Open(path, "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    # endif

    tree = root_file.Get(TREE_NAME)
    if tree is None:
        root_file.Close()
        raise KeyError(f"Tree '{TREE_NAME}' not found in {path}")
    # endif

    log(f"Reading {path}:{TREE_NAME} with PyROOT")

    try:
        for event in tree:
            detector1 = int(getattr(event, "detector1"))
            detector2 = int(getattr(event, "detector2"))

            topology = next(
                (
                    candidate
                    for candidate in topologies
                    if detector1 == candidate.detector1
                    and detector2 == candidate.detector2
                ),
                None,
            )

            if topology is None:
                continue
            # endif

            t1 = float(getattr(event, "t1"))
            open_angle = float(getattr(event, "open_angle_ep2"))
            if not math.isfinite(t1) or not math.isfinite(open_angle):
                continue
            # endif
            if (-t1) >= T1_ABS_MAX or open_angle <= OPEN_ANGLE_MIN_DEG:
                continue
            # endif

            if REQUIRE_DISTINCT_FD_SECTORS:
                sectors = [
                    int(fd_sector_from_phi(np.asarray([float(getattr(event, "e_phi"))]))[0])
                ]
                if topology.detector1 == 1:
                    sectors.append(
                        int(fd_sector_from_phi(np.asarray([float(getattr(event, "p1_phi"))]))[0])
                    )
                # endif
                if topology.detector2 == 1:
                    sectors.append(
                        int(fd_sector_from_phi(np.asarray([float(getattr(event, "p2_phi"))]))[0])
                    )
                # endif
                if any(sector < 1 for sector in sectors) or len(set(sectors)) != len(sectors):
                    continue
                # endif
            # endif

            if apply_mx2_1_cut:
                mx2_1_value = float(getattr(event, resolved["Mx2_1"]))
                if not math.isfinite(mx2_1_value) or mx2_1_value >= mx2_1_upper_cut(topology):
                    continue
                # endif
            # endif

            selected_events[topology.key] += 1

            for variable in variables:
                value = float(getattr(event, resolved[variable.branch]))
                if not math.isfinite(value):
                    continue
                # endif

                if value < variable.xmin or value >= variable.xmax:
                    continue
                # endif

                fraction = (value - variable.xmin) / (variable.xmax - variable.xmin)
                index = int(fraction * variable.bins)
                if 0 <= index < variable.bins:
                    histograms[topology.key][variable.branch][index] += 1.0
                # endif
            # endfor
        # endfor
    finally:
        root_file.Close()
    # endtry

    return histograms, selected_events


def fill_histograms_for_file(
    path: str,
    topologies: Sequence[TopologyConfig],
    step_size: str,
    apply_mx2_1_cut: bool = False,
    variables: Sequence[VariableConfig] = VARIABLES,
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, int]]:
    backend = io_backend()

    if backend == "uproot":
        return fill_histograms_uproot(path, topologies, step_size, apply_mx2_1_cut, variables)
    # endif

    return fill_histograms_pyroot(path, topologies, step_size, apply_mx2_1_cut, variables)


def normalize_density(
    counts: np.ndarray,
    variable: VariableConfig,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return unit-integral density and Poisson uncertainty inside plot range."""

    counts = np.asarray(counts, dtype=np.float64)
    total = float(np.sum(counts))
    bin_width = (variable.xmax - variable.xmin) / float(variable.bins)

    if total <= 0.0 or not math.isfinite(total):
        zeros = np.zeros_like(counts)
        return zeros, zeros
    # endif

    density = counts / (total * bin_width)
    uncertainty = np.sqrt(counts) / (total * bin_width)
    return density, uncertainty


def bin_geometry(variable: VariableConfig) -> Tuple[np.ndarray, np.ndarray]:
    edges = np.linspace(variable.xmin, variable.xmax, variable.bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return edges, centers


def positive_y_floor(*arrays: np.ndarray) -> Optional[float]:
    positive_parts = [array[array > 0.0] for array in arrays if array.size > 0]
    positive_parts = [array for array in positive_parts if array.size > 0]

    if not positive_parts:
        return None
    # endif

    minimum = min(float(np.min(array)) for array in positive_parts)
    return max(minimum * 0.5, 1.0e-8)


def normalized_shape(counts: np.ndarray) -> Optional[np.ndarray]:
    values = np.asarray(counts, dtype=np.float64)
    total = float(np.sum(values))
    if total <= 0.0 or not math.isfinite(total):
        return None
    # endif
    return values / total


def is_lower_bounded_morph_variable(variable: VariableConfig) -> bool:
    """Variables peaked near the lower boundary with a tail toward larger values."""
    return variable.branch in {
        "pTmiss",
        "theta_gamma_gamma",
        "theta_pi0_pi0",
    }


def is_upper_bounded_morph_variable(variable: VariableConfig) -> bool:
    """Variables peaked near the upper boundary with a tail toward smaller values."""
    return variable.branch in {
        "z2",
        "xF2",
    }


def is_logit_morph_variable(variable: VariableConfig) -> bool:
    """Variables better represented by a bounded logit coordinate."""
    return variable.branch == "z2"


def is_log_morph_variable(variable: VariableConfig) -> bool:
    return (
        is_lower_bounded_morph_variable(variable)
        or is_upper_bounded_morph_variable(variable)
        or is_logit_morph_variable(variable)
    )



def is_asymmetric_additive_variable(variable: VariableConfig) -> bool:
    """Interior-peaked variables needing different left/right broadening."""
    return variable.branch in {"theta", "xF"}


def mx2_1_upper_cut(topology: TopologyConfig) -> float:
    """Topology-dependent hard upper cut used for the second-stage plots and fits."""
    if topology.key == "FD_FD":
        return 0.20
    # endif
    if topology.key == "CD_FD":
        return 0.30
    # endif
    if topology.key == "CD_FT":
        return 0.40
    # endif
    raise ValueError(f"Unsupported topology for Mx2_1 cut: {topology.key}")


def fit_mask_for_variable(
    variable: VariableConfig,
    topology: TopologyConfig,
) -> np.ndarray:
    """Bins used in the shared likelihood; all bins remain visible in plots."""
    _, centers = bin_geometry(variable)
    mask = np.ones(variable.bins, dtype=bool)

    if variable.branch == "xF":
        # Avoid allowing small range-clipping effects to drive the nuisance fit.
        mask[:2] = False
        mask[-2:] = False
    # endif

    return mask


def transform_additive_shape(
    base_shape: np.ndarray,
    variable: VariableConfig,
    shift: float,
    sigma_add: float,
) -> Optional[np.ndarray]:
    """Shift and Gaussian-smear a template without truncating its visible tail.

    Each source-bin probability is transported into every target bin using the
    integrated Gaussian probability, avoiding hard interpolation cutoffs.
    """
    edges, centers = bin_geometry(variable)
    source_weights = np.asarray(base_shape, dtype=np.float64)

    if sigma_add <= 1.0e-10:
        transformed_centers = centers + shift
        target, _ = np.histogram(
            transformed_centers, bins=edges, weights=source_weights
        )
    else:
        lower = edges[:-1, None]
        upper = edges[1:, None]
        source_means = (centers + shift)[None, :]
        probabilities = ndtr((upper - source_means) / sigma_add) - ndtr(
            (lower - source_means) / sigma_add
        )
        target = probabilities @ source_weights
    # endif

    target = np.clip(target, 0.0, None)
    total = float(np.sum(target))
    if total <= 0.0 or not math.isfinite(total):
        return None
    # endif
    return target / total


def transform_asymmetric_additive_shape(
    base_shape: np.ndarray,
    variable: VariableConfig,
    shift: float,
    sigma_left: float,
    sigma_right: float,
) -> Optional[np.ndarray]:
    """Shift a template and smear it with a normalized split Gaussian.

    The left and right widths are independent. For theta this allows a broad
    low-angle shoulder while retaining the comparatively sharp high-angle edge.
    """
    edges, centers = bin_geometry(variable)
    source_weights = np.asarray(base_shape, dtype=np.float64)
    sigma_left = max(float(sigma_left), 1.0e-10)
    sigma_right = max(float(sigma_right), 1.0e-10)
    source_means = (centers + shift)[None, :]

    def split_cdf(x: np.ndarray) -> np.ndarray:
        z_left = (x - source_means) / sigma_left
        z_right = (x - source_means) / sigma_right
        left_weight = sigma_left / (sigma_left + sigma_right)
        right_weight = sigma_right / (sigma_left + sigma_right)

        left_value = 2.0 * left_weight * ndtr(z_left)
        right_value = (
            left_weight
            + 2.0 * right_weight * (ndtr(z_right) - 0.5)
        )
        return np.where(x <= source_means, left_value, right_value)

    lower = edges[:-1, None]
    upper = edges[1:, None]
    probabilities = split_cdf(upper) - split_cdf(lower)
    probabilities = np.clip(probabilities, 0.0, 1.0)
    target = probabilities @ source_weights

    target = np.clip(target, 0.0, None)
    total = float(np.sum(target))
    if total <= 0.0 or not math.isfinite(total):
        return None
    # endif
    return target / total


def transform_positive_shape(
    base_shape: np.ndarray,
    variable: VariableConfig,
    log_shift: float,
    log_sigma: float,
) -> Optional[np.ndarray]:
    """Morph a positive template through log(x+epsilon) shift and smearing."""
    edges, centers = bin_geometry(variable)
    epsilon = max(0.25 * (edges[1] - edges[0]), 1.0e-8)
    source_weights = np.asarray(base_shape, dtype=np.float64)
    target = np.zeros(variable.bins, dtype=np.float64)

    if log_sigma <= 1.0e-8:
        transformed_centers = np.exp(np.log(centers + epsilon) + log_shift) - epsilon
        target, _ = np.histogram(transformed_centers, bins=edges, weights=source_weights)
    else:
        lower_log = np.log(np.maximum(edges[:-1] + epsilon, 1.0e-15))[:, None]
        upper_log = np.log(np.maximum(edges[1:] + epsilon, 1.0e-15))[:, None]
        source_means = (np.log(centers + epsilon) + log_shift)[None, :]
        probabilities = ndtr((upper_log - source_means) / log_sigma) - ndtr(
            (lower_log - source_means) / log_sigma
        )
        target = probabilities @ source_weights
    # endif

    target = np.clip(target, 0.0, None)
    total = float(np.sum(target))
    if total <= 0.0 or not math.isfinite(total):
        return None
    # endif
    return target / total


def transform_logit_bounded_shape(
    base_shape: np.ndarray,
    variable: VariableConfig,
    logit_shift: float,
    logit_sigma: float,
) -> Optional[np.ndarray]:
    """Morph a bounded template in logit((x-xmin)/(xmax-x)) space."""

    edges, centers = bin_geometry(variable)
    width = float(variable.xmax - variable.xmin)
    epsilon = max(0.25 * (edges[1] - edges[0]), 1.0e-8)
    source_weights = np.asarray(base_shape, dtype=np.float64)

    def to_logit(values: np.ndarray) -> np.ndarray:
        scaled = (values - variable.xmin + epsilon) / (width + 2.0 * epsilon)
        scaled = np.clip(scaled, 1.0e-12, 1.0 - 1.0e-12)
        return np.log(scaled / (1.0 - scaled))

    source_means = to_logit(centers) + logit_shift
    if logit_sigma <= 1.0e-8:
        transformed_scaled = 1.0 / (1.0 + np.exp(-source_means))
        transformed_centers = variable.xmin - epsilon + transformed_scaled * (width + 2.0 * epsilon)
        target, _ = np.histogram(transformed_centers, bins=edges, weights=source_weights)
    else:
        lower_logit = to_logit(edges[:-1])[:, None]
        upper_logit = to_logit(edges[1:])[:, None]
        probabilities = ndtr((upper_logit - source_means[None, :]) / logit_sigma) - ndtr((lower_logit - source_means[None, :]) / logit_sigma)
        target = probabilities @ source_weights
    # endif

    target = np.clip(target, 0.0, None)
    total = float(np.sum(target))
    if total <= 0.0 or not math.isfinite(total):
        return None
    # endif
    return target / total


def transform_upper_bounded_shape(
    base_shape: np.ndarray,
    variable: VariableConfig,
    log_shift: float,
    log_sigma: float,
) -> Optional[np.ndarray]:
    """Morph an upper-edge-peaked template in log(xmax-x+epsilon).

    This preserves the upper plotting boundary and broadens primarily toward
    smaller values, which is appropriate for z2 and xF2.
    """
    edges, centers = bin_geometry(variable)
    upper_bound = float(variable.xmax)
    epsilon = max(0.25 * (edges[1] - edges[0]), 1.0e-8)
    source_weights = np.asarray(base_shape, dtype=np.float64)

    source_distance = np.maximum(upper_bound - centers + epsilon, 1.0e-15)
    source_log_means = np.log(source_distance) + log_shift

    if log_sigma <= 1.0e-8:
        transformed_distance = np.exp(source_log_means)
        transformed_centers = upper_bound + epsilon - transformed_distance
        target, _ = np.histogram(
            transformed_centers,
            bins=edges,
            weights=source_weights,
        )
    else:
        # For target x in [edge_i, edge_{i+1}], the complementary distance
        # u=xmax+epsilon-x lies in [u_low, u_high] with reversed x ordering.
        distance_low = np.maximum(
            upper_bound + epsilon - edges[1:],
            1.0e-15,
        )[:, None]
        distance_high = np.maximum(
            upper_bound + epsilon - edges[:-1],
            1.0e-15,
        )[:, None]

        lower_log = np.log(distance_low)
        upper_log = np.log(distance_high)
        probabilities = ndtr(
            (upper_log - source_log_means[None, :]) / log_sigma
        ) - ndtr(
            (lower_log - source_log_means[None, :]) / log_sigma
        )
        target = probabilities @ source_weights
    # endif

    target = np.clip(target, 0.0, None)
    total = float(np.sum(target))
    if total <= 0.0 or not math.isfinite(total):
        return None
    # endif
    return target / total


def transform_dvcs_shape(
    base_shape: np.ndarray,
    variable: VariableConfig,
    shift: float,
    sigma_add: float,
    sigma_right: Optional[float] = None,
) -> Optional[np.ndarray]:
    if is_asymmetric_additive_variable(variable):
        right_width = sigma_add if sigma_right is None else sigma_right
        return transform_asymmetric_additive_shape(
            base_shape,
            variable,
            shift,
            sigma_add,
            right_width,
        )
    # endif
    if is_lower_bounded_morph_variable(variable):
        return transform_positive_shape(base_shape, variable, shift, sigma_add)
    # endif
    if is_logit_morph_variable(variable):
        return transform_logit_bounded_shape(base_shape, variable, shift, sigma_add)
    # endif
    if is_upper_bounded_morph_variable(variable):
        return transform_upper_bounded_shape(
            base_shape,
            variable,
            shift,
            sigma_add,
        )
    # endif
    return transform_additive_shape(base_shape, variable, shift, sigma_add)


def poisson_deviance(observed: np.ndarray, expected: np.ndarray) -> float:
    expected = np.clip(np.asarray(expected, dtype=np.float64), 1.0e-12, None)
    observed = np.asarray(observed, dtype=np.float64)
    positive = observed > 0.0
    terms = expected - observed
    terms[positive] += observed[positive] * np.log(
        observed[positive] / expected[positive]
    )
    return 2.0 * float(np.sum(terms))


def fit_shared_two_templates(
    data_histograms: Mapping[str, np.ndarray],
    dvcs_histograms: Mapping[str, np.ndarray],
    pi0_histograms: Mapping[str, np.ndarray],
    topology: TopologyConfig,
    max_shift_bins: float,
    max_smear_bins: float,
    min_counts: int,
    fraction_variable_branches: Sequence[str],
    shift_prior_bins: float,
    smear_prior_bins: float,
    use_nuisance_penalties: bool,
    core_containment: float,
    fraction_containment: float,
    pi0_core_calibration: Optional[Mapping[str, Tuple[float, float]]] = None,
    outside_overshoot_penalty_weight: float = 0.25,
    emiss2_mean_order_penalty_weight: float = 25.0,
) -> SharedFitSummary:
    """Fit a shared f_pi0 while separating resolution and fraction regions.

    Variable-specific DVCS shifts and smearings are fitted only in the narrower
    DVCS-MC signal core. The shared fraction is then updated using the broader
    containment regions of the selected discriminator variables. Validation
    variables profile only their nuisance parameters at the fitted fraction.
    """
    requested_fraction_variables = tuple(dict.fromkeys(fraction_variable_branches))
    requested_set = set(requested_fraction_variables)
    prepared: Dict[str, Dict[str, object]] = {}
    active_variables: List[VariableConfig] = []

    for variable in VARIABLES:
        data = np.asarray(data_histograms[variable.branch], dtype=np.float64)
        dvcs_counts = np.asarray(dvcs_histograms[variable.branch], dtype=np.float64)
        dvcs_shape = normalized_shape(dvcs_counts)
        pi0_shape = normalized_shape(pi0_histograms[variable.branch])
        data_total = float(np.sum(data))
        if data_total < min_counts or dvcs_shape is None or pi0_shape is None:
            continue
        # endif
        core_mask = mc_signal_containment_mask(
            dvcs_counts, variable, topology, core_containment
        )
        fraction_mask = mc_signal_containment_mask(
            dvcs_counts, variable, topology, fraction_containment
        )
        # The fraction region must include the nuisance core.
        fraction_mask = fraction_mask | core_mask
        prepared[variable.branch] = {
            "data": data,
            "data_total": data_total,
            "dvcs_shape": dvcs_shape,
            "pi0_shape": pi0_shape,
            "core_mask": core_mask,
            "fraction_mask": fraction_mask,
        }
        active_variables.append(variable)
    # endfor

    fraction_variables = [v for v in active_variables if v.branch in requested_set]
    if not fraction_variables:
        return SharedFitSummary(False, "none of the requested fraction variables has sufficient data and MC")
    # endif

    calibration = dict(pi0_core_calibration or {})

    def calibration_key(variable: VariableConfig) -> str:
        return "theta_pi0_pi0" if variable.branch == "theta_gamma_gamma" else variable.branch

    def calibrated_shift_center(variable: VariableConfig) -> float:
        entry = calibration.get(calibration_key(variable))
        if entry is None:
            return 0.0
        # endif
        value = float(entry[0])
        return value if math.isfinite(value) else 0.0

    def nuisance_bounds(variable: VariableConfig) -> List[Tuple[float, float]]:
        if variable.branch in {"Mx2", "Mx2_1"}:
            bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
            return [(0.0, 0.0), (0.0, max_smear_bins * bin_width)]
        # endif
        if is_log_morph_variable(variable):
            center = calibrated_shift_center(variable)
            return [
                (max(-0.70, center - 0.40), min(0.70, center + 0.40)),
                (0.0, 1.00),
            ]
        # endif
        bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
        center = calibrated_shift_center(variable)
        half_range = max_shift_bins * bin_width
        if is_asymmetric_additive_variable(variable):
            return [
                (center - half_range, center + half_range),
                (0.0, max_smear_bins * bin_width),
                (0.0, max_smear_bins * bin_width),
            ]
        # endif
        return [
            (center - half_range, center + half_range),
            (0.0, max_smear_bins * bin_width),
        ]

    def nuisance_start(variable: VariableConfig) -> np.ndarray:
        center = calibrated_shift_center(variable)
        if is_log_morph_variable(variable):
            return np.asarray([center, 0.10], dtype=np.float64)
        # endif
        bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
        if is_asymmetric_additive_variable(variable):
            if variable.branch == "xF":
                return np.asarray([center, 1.5 * bin_width, 1.5 * bin_width], dtype=np.float64)
            # endif
            return np.asarray([center, 2.0 * bin_width, 0.75 * bin_width], dtype=np.float64)
        # endif
        return np.asarray([center, 2.0 * bin_width], dtype=np.float64)

    def nuisance_penalty(
        variable: VariableConfig,
        nuisance: np.ndarray,
    ) -> float:
        if not use_nuisance_penalties:
            return 0.0
        # endif
        shift = float(nuisance[0])
        sigma_left = float(nuisance[1])
        sigma_right = (
            float(nuisance[2])
            if is_asymmetric_additive_variable(variable)
            else sigma_left
        )
        shift_center = calibrated_shift_center(variable)
        if is_log_morph_variable(variable):
            shift_width = 0.20
            smear_width = 0.40
        else:
            bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
            shift_width = max(shift_prior_bins * bin_width, 1.0e-12)
            smear_width = max(smear_prior_bins * bin_width, 1.0e-12)
        # endif
        penalty = 0.5 * ((shift - shift_center) / shift_width) ** 2
        penalty += 0.5 * (sigma_left / smear_width) ** 2
        if is_asymmetric_additive_variable(variable):
            if variable.branch == "theta":
                right_prior = max(0.5 * smear_width, 1.0e-12)
            else:
                right_prior = smear_width
            # endif
            penalty += 0.5 * (sigma_right / right_prior) ** 2
        # endif
        return penalty

    def build_variable_model(
        variable: VariableConfig,
        fraction: float,
        nuisance: np.ndarray,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        info = prepared[variable.branch]
        shift = float(nuisance[0])
        sigma_left = float(nuisance[1])
        sigma_right = (
            float(nuisance[2])
            if is_asymmetric_additive_variable(variable)
            else None
        )
        transformed = transform_dvcs_shape(
            info["dvcs_shape"],
            variable,
            shift,
            sigma_left,
            sigma_right,
        )
        if transformed is None:
            return None
        # endif
        pi0_shape = info["pi0_shape"]
        dvcs_shape_component = np.clip(
            (1.0 - fraction) * transformed,
            0.0,
            None,
        )
        pi0_shape_component = np.clip(
            fraction * pi0_shape,
            0.0,
            None,
        )
        total_shape = dvcs_shape_component + pi0_shape_component
        if float(np.sum(total_shape)) <= 0.0:
            return None
        # endif
        return (
            total_shape,
            dvcs_shape_component,
            pi0_shape_component,
            transformed,
        )

    def physical_shape_penalty(
        variable: VariableConfig,
        transformed_signal: np.ndarray,
        pi0_shape: np.ndarray,
    ) -> float:
        """Apply targeted physical constraints to pathological nuisance morphs."""

        if (
            variable.branch != "Emiss2"
            or emiss2_mean_order_penalty_weight <= 0.0
        ):
            return 0.0
        # endif

        dvcs_mean, _ = distribution_moments(transformed_signal, variable)
        pi0_mean, _ = distribution_moments(pi0_shape, variable)
        bin_width = (variable.xmax - variable.xmin) / float(variable.bins)

        # The DVCS signal should not be shifted farther from zero than the
        # missing-particle pi0 background. Permit one bin of numerical tolerance.
        violation = abs(dvcs_mean) - abs(pi0_mean) - bin_width
        if violation <= 0.0:
            return 0.0
        # endif

        return (
            0.5
            * emiss2_mean_order_penalty_weight
            * (violation / max(bin_width, 1.0e-12)) ** 2
        )

    def outside_region_overshoot_penalty(
        data: np.ndarray,
        total_shape: np.ndarray,
        active_mask: np.ndarray,
    ) -> float:
        """Penalize only model overshoot outside the fitted region."""

        if outside_overshoot_penalty_weight <= 0.0 or np.all(active_mask):
            return 0.0
        # endif

        shape_in_mask = float(np.sum(total_shape[active_mask]))
        data_in_mask = float(np.sum(data[active_mask]))
        if shape_in_mask <= 0.0 or data_in_mask <= 0.0:
            return 0.0
        # endif

        scale = data_in_mask / shape_in_mask
        expected = scale * total_shape
        excluded = ~active_mask
        excess = np.clip(expected[excluded] - data[excluded], 0.0, None)
        variance = np.maximum(data[excluded] + 1.0, 1.0)
        return (
            0.5
            * outside_overshoot_penalty_weight
            * float(np.sum((excess ** 2) / variance))
        )

    def objective_for_mask(
        variable: VariableConfig,
        fraction: float,
        nuisance: np.ndarray,
        mask_name: str,
        include_penalty: bool,
    ) -> float:
        built = build_variable_model(variable, fraction, nuisance)
        if built is None:
            return 1.0e100
        # endif
        info = prepared[variable.branch]
        mask = np.asarray(info[mask_name], dtype=bool)
        total_shape = built[0]
        shape_in_mask = float(np.sum(total_shape[mask]))
        data_in_mask = float(np.sum(info["data"][mask]))
        if shape_in_mask <= 0.0 or data_in_mask <= 0.0:
            return 1.0e100
        # endif
        model_in_mask = data_in_mask * total_shape[mask] / shape_in_mask
        value = 0.5 * poisson_deviance(
            info["data"][mask],
            model_in_mask,
        )

        # The outside-region term is one-sided: unmodelled positive data tails
        # are tolerated, but an extrapolated model is not allowed to vastly
        # exceed the observed data where those bins were excluded from the fit.
        value += outside_region_overshoot_penalty(
            info["data"],
            total_shape,
            mask,
        )

        if include_penalty:
            value += nuisance_penalty(variable, nuisance)
            value += physical_shape_penalty(
                variable,
                built[3],
                info["pi0_shape"],
            )
        # endif
        return value

    best_value = math.inf
    best_fraction = math.nan
    best_driver_nuisances: Dict[str, np.ndarray] = {}

    for initial_fraction in (0.10, 0.30, 0.60, 0.85):
        fraction = initial_fraction
        nuisances = {v.branch: nuisance_start(v) for v in fraction_variables}
        previous_value = math.inf

        for iteration in range(15):
            # Resolution step: fit each DVCS nuisance only in its 90% core.
            for variable in fraction_variables:
                result = minimize(
                    lambda values, v=variable, f=fraction: objective_for_mask(
                        v, f, values, "core_mask", True
                    ),
                    nuisances[variable.branch],
                    method="L-BFGS-B",
                    bounds=nuisance_bounds(variable),
                    options={"maxiter": 400, "ftol": 1.0e-9},
                )
                if result.success and np.all(np.isfinite(result.x)):
                    nuisances[variable.branch] = np.asarray(result.x, dtype=np.float64)
                # endif
            # endfor

            # Fraction step: hold the core-derived nuisances fixed and use the
            # broader 95% discriminator regions to update the shared fraction.
            def fraction_objective(candidate_fraction: float) -> float:
                return sum(
                    objective_for_mask(
                        v, candidate_fraction, nuisances[v.branch],
                        "fraction_mask", False,
                    )
                    for v in fraction_variables
                )

            fraction_result = minimize_scalar(
                fraction_objective,
                bounds=(0.0, 1.0),
                method="bounded",
                options={"xatol": 2.0e-5, "maxiter": 200},
            )
            if fraction_result.success and math.isfinite(float(fraction_result.x)):
                fraction = float(fraction_result.x)
            # endif
            current_value = fraction_objective(fraction)
            if abs(previous_value - current_value) <= 1.0e-7 * max(1.0, current_value):
                break
            # endif
            previous_value = current_value
        # endfor

        value = sum(
            objective_for_mask(v, fraction, nuisances[v.branch], "fraction_mask", False)
            for v in fraction_variables
        )
        if value < best_value:
            best_value = value
            best_fraction = fraction
            best_driver_nuisances = {k: val.copy() for k, val in nuisances.items()}
        # endif
    # endfor

    if not math.isfinite(best_fraction):
        return SharedFitSummary(False, "selected-discriminator core/profile fit failed")
    # endif

    # Profile every validation variable in its own DVCS core at fixed f_pi0.
    all_nuisances: Dict[str, np.ndarray] = {}
    for variable in active_variables:
        start = best_driver_nuisances.get(variable.branch, nuisance_start(variable))
        result = minimize(
            lambda values, v=variable: objective_for_mask(
                v, best_fraction, values, "core_mask", True
            ),
            start,
            method="L-BFGS-B",
            bounds=nuisance_bounds(variable),
            options={"maxiter": 500, "ftol": 1.0e-10},
        )
        all_nuisances[variable.branch] = (
            np.asarray(result.x, dtype=np.float64)
            if result.success and np.all(np.isfinite(result.x))
            else start
        )
    # endfor

    fraction_error = math.nan
    try:
        step = 1.0e-3
        if step < best_fraction < 1.0 - step:
            def fixed_nuisance_fraction_objective(candidate_fraction: float) -> float:
                return sum(
                    objective_for_mask(
                        v, candidate_fraction, all_nuisances[v.branch],
                        "fraction_mask", False,
                    )
                    for v in fraction_variables
                )
            left = fixed_nuisance_fraction_objective(best_fraction - step)
            center = fixed_nuisance_fraction_objective(best_fraction)
            right = fixed_nuisance_fraction_objective(best_fraction + step)
            curvature = (left - 2.0 * center + right) / (step ** 2)
            if curvature > 0.0:
                fraction_error = 1.0 / math.sqrt(curvature)
            # endif
        # endif
    except Exception:
        pass
    # endtry

    results: Dict[str, FitResult] = {}
    driver_deviance = 0.0
    driver_bins = 0
    for variable in VARIABLES:
        if variable.branch not in prepared:
            results[variable.branch] = FitResult(False, "insufficient data or empty template")
            continue
        # endif
        info = prepared[variable.branch]
        nuisance = all_nuisances[variable.branch]
        shift = float(nuisance[0])
        sigma_add = float(nuisance[1])
        sigma_right = (
            float(nuisance[2])
            if is_asymmetric_additive_variable(variable)
            else math.nan
        )
        built = build_variable_model(variable, best_fraction, nuisance)
        if built is None:
            results[variable.branch] = FitResult(False, "invalid fitted template")
            continue
        # endif
        total_shape, dvcs_shape_component, pi0_shape_component, transformed = built
        display_mask = np.asarray(
            info["fraction_mask"] if variable.branch in requested_set else info["core_mask"],
            dtype=bool,
        )
        display_shape_sum = float(np.sum(total_shape[display_mask]))
        display_data_sum = float(np.sum(info["data"][display_mask]))
        if display_shape_sum <= 0.0 or display_data_sum <= 0.0:
            results[variable.branch] = FitResult(
                False,
                "invalid fit-region normalization",
            )
            continue
        # endif
        scale = display_data_sum / display_shape_sum
        model = scale * total_shape
        dvcs_component = scale * dvcs_shape_component
        pi0_component = scale * pi0_shape_component
        variable_deviance = poisson_deviance(
            info["data"][display_mask],
            model[display_mask],
        )
        used_bins = int(np.count_nonzero(display_mask))
        if variable.branch in requested_set:
            driver_deviance += variable_deviance
            driver_bins += used_bins
        # endif
        excluded = ~display_mask
        excluded_data = float(np.sum(info["data"][excluded]))
        excluded_model = float(np.sum(model[excluded]))
        excluded_excess = excluded_data - excluded_model
        excluded_fraction = (
            excluded_excess / float(info["data_total"])
            if float(info["data_total"]) > 0.0 else 0.0
        )
        role = (
            f"fraction driver ({100.0 * fraction_containment:.0f}% region)"
            if variable.branch in requested_set
            else f"validation core ({100.0 * core_containment:.0f}%)"
        )
        results[variable.branch] = FitResult(
            success=True,
            message=role,
            f_pi0=best_fraction,
            f_pi0_err=fraction_error,
            shift=shift,
            sigma_add=sigma_add,
            sigma_right=sigma_right,
            deviance=variable_deviance,
            ndf=max(0, used_bins - 2),
            data_total=float(info["data_total"]),
            model_counts=model,
            dvcs_component_counts=dvcs_component,
            pi0_component_counts=pi0_component,
            transformed_dvcs_shape=transformed,
            fit_mask=display_mask,
            morph_label=("asymmetric-additive" if is_asymmetric_additive_variable(variable) else "lower-log-space" if is_lower_bounded_morph_variable(variable) else "logit-space" if is_logit_morph_variable(variable) else "upper-log-space" if is_upper_bounded_morph_variable(variable) else "additive"),
            excluded_data_counts=excluded_data,
            excluded_model_counts=excluded_model,
            excluded_excess_counts=excluded_excess,
            excluded_excess_fraction=excluded_fraction,
        )
    # endfor

    n_driver_parameters = 1 + sum(
        3 if is_asymmetric_additive_variable(variable) else 2
        for variable in fraction_variables
    )
    return SharedFitSummary(
        success=True,
        message="DVCS-core nuisance / broader-region fraction fit converged",
        f_pi0=best_fraction,
        f_pi0_err=fraction_error,
        deviance=driver_deviance,
        ndf=max(0, driver_bins - n_driver_parameters),
        variable_results=results,
        fraction_variables=tuple(v.branch for v in fraction_variables),
    )



def draw_shape_canvas(
    output_path: Path,
    period: PeriodConfig,
    topology: TopologyConfig,
    data_histograms: Mapping[str, np.ndarray],
    dvcs_histograms: Mapping[str, np.ndarray],
    pi0_histograms: Mapping[str, np.ndarray],
    selected_counts: Mapping[str, int],
    log_y: bool,
    dpi: int,
    selection_label: str = "minimal preselection",
) -> None:
    """Draw the original independently unit-normalized shape comparison."""
    fig, axes = plt.subplots(3, 4, figsize=(18.0, 13.0))
    flat_axes = axes.ravel()

    for axis, variable in zip(flat_axes, DVCS_SHAPE_VARIABLES):
        edges, centers = bin_geometry(variable)
        data_shape = normalized_shape(data_histograms[variable.branch])
        dvcs_shape = normalized_shape(dvcs_histograms[variable.branch])
        pi0_shape = normalized_shape(pi0_histograms[variable.branch])

        if data_shape is not None:
            counts = np.asarray(data_histograms[variable.branch], dtype=np.float64)
            total = float(np.sum(counts))
            errors = np.sqrt(counts) / total if total > 0.0 else np.zeros_like(counts)
            axis.errorbar(
                centers, data_shape, yerr=errors, fmt="o", markersize=2.4,
                linewidth=0.8, capsize=0.0, color=SAMPLE_COLORS["data"],
                label=SAMPLE_LABELS["data"], zorder=4,
            )
        # endif
        if dvcs_shape is not None:
            axis.stairs(
                dvcs_shape, edges, color=SAMPLE_COLORS["dvcs_mc"],
                linewidth=1.7, label="DVCS MC", zorder=3,
            )
        # endif
        if pi0_shape is not None:
            axis.stairs(
                pi0_shape, edges, color=SAMPLE_COLORS["pi0_mc"],
                linewidth=1.7, label=r"$e\pi^0$ MC as DVCS", zorder=2,
            )
        # endif

        axis.set_xlim(variable.xmin, variable.xmax)
        axis.set_xlabel(variable.label)
        axis.set_ylabel("unit-normalized entries / bin")
        axis.grid(axis="y", alpha=0.25)
        if log_y:
            floor = positive_y_floor(
                *[array for array in (data_shape, dvcs_shape, pi0_shape) if array is not None]
            )
            if floor is not None:
                axis.set_yscale("log")
                axis.set_ylim(bottom=max(1.0e-7, floor))
            # endif
        else:
            axis.set_ylim(bottom=0.0)
        # endif
    # endfor

    handles, labels = flat_axes[0].get_legend_handles_labels()
    fig.suptitle(
        f"Exclusivity shapes after {selection_label}: {period.label}, {topology.label}\n"
        f"selected entries: data={selected_counts['data']:,}, "
        f"DVCS MC={selected_counts['dvcs_mc']:,}, "
        rf"$e\pi^0$ MC as DVCS={selected_counts['pi0_mc']:,}",
        fontsize=15, y=0.975,
    )
    fig.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.905),
        ncol=len(labels), frameon=False,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.89))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def draw_fit_canvas(
    output_path: Path,
    period: PeriodConfig,
    topology: TopologyConfig,
    data_histograms: Mapping[str, np.ndarray],
    dvcs_histograms: Mapping[str, np.ndarray],
    pi0_histograms: Mapping[str, np.ndarray],
    selected_counts: Mapping[str, int],
    fit_results: Mapping[str, FitResult],
    shared_summary: SharedFitSummary,
    log_y: bool,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(18.0, 9.3))
    flat_axes = axes.ravel()

    for axis, variable in zip(flat_axes, VARIABLES):
        edges, centers = bin_geometry(variable)
        data_counts = np.asarray(data_histograms[variable.branch], dtype=np.float64)
        data_error = np.sqrt(data_counts)
        dvcs_shape = normalized_shape(dvcs_histograms[variable.branch])
        result = fit_results[variable.branch]

        axis.errorbar(
            centers, data_counts, yerr=data_error, fmt="o", markersize=2.4,
            linewidth=0.8, capsize=0.0, color=SAMPLE_COLORS["data"],
            label=SAMPLE_LABELS["data"], zorder=5,
        )

        if dvcs_shape is not None:
            axis.stairs(
                result.data_total * dvcs_shape, edges, color=SAMPLE_COLORS["raw_dvcs"],
                linewidth=1.2, linestyle=":", label="raw DVCS MC shape", zorder=1,
            )
        # endif

        if result.success:
            if result.fit_mask is not None and not np.all(result.fit_mask):
                excluded = ~result.fit_mask
                for index in np.flatnonzero(excluded):
                    axis.axvspan(edges[index], edges[index + 1], color="0.85", alpha=0.35, linewidth=0)
                # endfor
            # endif
            axis.stairs(
                result.dvcs_component_counts, edges, color=SAMPLE_COLORS["dvcs_mc"],
                linewidth=1.6, label="fitted DVCS component", zorder=3,
            )
            axis.stairs(
                result.pi0_component_counts, edges, color=SAMPLE_COLORS["pi0_mc"],
                linewidth=1.6, label=r"fitted $e\pi^0$ component", zorder=2,
            )
            axis.stairs(
                result.model_counts, edges, color=SAMPLE_COLORS["fit"],
                linewidth=2.0, linestyle="--", label="total two-template fit", zorder=4,
            )
            quality = (
                rf"$f_{{\pi^0}}={result.f_pi0:.3f}$ ({result.message})" + "\n"
                rf"$\Delta={result.shift:.4g}$" + "\n"
                + (
                    rf"$\sigma_L={result.sigma_add:.4g},\ "
                    rf"\sigma_R={result.sigma_right:.4g}$ "
                    rf"({result.morph_label})"
                    if is_asymmetric_additive_variable(variable)
                    else rf"$\sigma={result.sigma_add:.4g}$ "
                    rf"({result.morph_label})"
                )
                + "\n"
                + f"$D/ndf={result.deviance:.1f}/{result.ndf}$"
            )
            if result.excluded_data_counts > 0.0 and result.fit_mask is not None and not np.all(result.fit_mask):
                quality += (
                    "\n"
                    f"masked excess={result.excluded_excess_counts:.0f} "
                    f"({100.0 * result.excluded_excess_fraction:.2f}% of data)"
                )
            # endif
        else:
            quality = f"fit failed: {result.message}"
        # endif

        axis.text(
            0.98, 0.96, quality, transform=axis.transAxes, ha="right", va="top",
            fontsize=8.5, bbox=dict(facecolor="white", alpha=0.78, edgecolor="none"),
        )
        axis.set_xlim(variable.xmin, variable.xmax)
        axis.set_xlabel(variable.label)
        axis.set_ylabel("events / bin")
        axis.grid(axis="y", alpha=0.25)

        if log_y:
            arrays = [data_counts]
            if result.success:
                arrays.extend([result.model_counts, result.dvcs_component_counts, result.pi0_component_counts])
            # endif
            floor = positive_y_floor(*arrays)
            if floor is not None:
                axis.set_yscale("log")
                axis.set_ylim(bottom=max(0.5, floor))
            # endif
        else:
            axis.set_ylim(bottom=0.0)
        # endif
    # endfor

    handles, labels = flat_axes[0].get_legend_handles_labels()
    fig.suptitle(
        f"DVCS-core two-template fits after minimal preselection: {period.label}, {topology.label}\n"
        f"fraction drivers (common population): {', '.join(shared_summary.fraction_variables)}; "
        f"topology-selected entries: data={selected_counts['data']:,}, "
        f"DVCS MC={selected_counts['dvcs_mc']:,}, "
        rf"$e\pi^0$ MC as DVCS={selected_counts['pi0_mc']:,}",
        fontsize=15, y=0.992,
    )
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.875),
               ncol=len(labels), frameon=False)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.82))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)




def mc_signal_containment_mask(
    mc_counts: np.ndarray,
    variable: VariableConfig,
    topology: TopologyConfig,
    containment: float,
) -> np.ndarray:
    """Return an MC-defined containment mask for a reconstructed signal template.

    Signed variables retain the central equal-tail containment interval.
    Lower-edge-peaked variables retain the interval from the lower plot edge
    upward. Upper-edge-peaked variables retain the interval from the upper
    plot edge downward. The CM theta variable uses a central additive core. Existing variable-specific
    edge masks are also respected.
    """
    if containment >= 0.999999:
        return np.ones(variable.bins, dtype=bool)
    # endif

    counts = np.asarray(mc_counts, dtype=np.float64)
    total = float(np.sum(counts))
    base_mask = fit_mask_for_variable(variable, topology)
    if total <= 0.0 or not math.isfinite(total):
        return base_mask
    # endif

    cumulative = np.cumsum(counts) / total
    if is_lower_bounded_morph_variable(variable):
        lower_index = 0
        upper_index = int(np.searchsorted(cumulative, containment, side="left"))
    elif is_upper_bounded_morph_variable(variable):
        lower_index = int(
            np.searchsorted(cumulative, 1.0 - containment, side="left")
        )
        upper_index = variable.bins - 1
    else:
        tail = 0.5 * (1.0 - containment)
        lower_index = int(np.searchsorted(cumulative, tail, side="left"))
        upper_index = int(np.searchsorted(cumulative, 1.0 - tail, side="left"))
    # endif

    lower_index = max(0, min(lower_index, variable.bins - 1))
    upper_index = max(lower_index, min(upper_index, variable.bins - 1))
    core_mask = np.zeros(variable.bins, dtype=bool)
    core_mask[lower_index : upper_index + 1] = True
    core_mask &= base_mask

    # Protect against very narrow or sparse MC cores.
    if np.count_nonzero(core_mask) < 5:
        return base_mask
    # endif
    return core_mask


def fit_single_template(
    data_counts: np.ndarray,
    mc_counts: np.ndarray,
    variable: VariableConfig,
    topology: TopologyConfig,
    max_shift_bins: float,
    max_smear_bins: float,
    min_counts: int,
    shift_prior_bins: float,
    smear_prior_bins: float,
    use_nuisance_penalties: bool,
    core_containment: float,
) -> FitResult:
    """Fit the MC-defined exclusive-pi0 signal core with nuisance morphing."""

    data = np.asarray(data_counts, dtype=np.float64)
    mc_shape = normalized_shape(mc_counts)
    data_total = float(np.sum(data))
    if data_total < min_counts or mc_shape is None:
        return FitResult(False, "insufficient counts", data_total=data_total)
    # endif

    mask = mc_signal_containment_mask(
        mc_counts,
        variable,
        topology,
        core_containment,
    )
    core_data_total = float(np.sum(data[mask]))
    if core_data_total < min_counts:
        return FitResult(
            False,
            "insufficient counts in MC-defined signal core",
            data_total=data_total,
            fit_mask=mask,
        )
    # endif

    bin_width = (variable.xmax - variable.xmin) / variable.bins
    log_morph = is_log_morph_variable(variable)
    asymmetric = is_asymmetric_additive_variable(variable)

    if log_morph:
        bounds = [(-0.8, 0.8), (0.0, 1.2)]
        x0 = np.asarray([0.0, 0.10], dtype=np.float64)
    elif asymmetric:
        bounds = [
            (-max_shift_bins * bin_width, max_shift_bins * bin_width),
            (0.0, max_smear_bins * bin_width),
            (0.0, max_smear_bins * bin_width),
        ]
        x0 = np.asarray(
            [0.0, 2.0 * bin_width, 0.75 * bin_width],
            dtype=np.float64,
        )
    else:
        bounds = [
            (-max_shift_bins * bin_width, max_shift_bins * bin_width),
            (0.0, max_smear_bins * bin_width),
        ]
        x0 = np.asarray([0.0, 0.5 * bin_width], dtype=np.float64)
    # endif

    def build(params: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        sigma_right = float(params[2]) if asymmetric else None
        transformed = transform_dvcs_shape(
            mc_shape,
            variable,
            float(params[0]),
            float(params[1]),
            sigma_right,
        )
        if transformed is None:
            return np.zeros_like(data), np.zeros_like(data)
        # endif
        norm = float(np.sum(transformed[mask]))
        if norm <= 0.0:
            return np.zeros_like(data), transformed
        # endif
        model = core_data_total * transformed / norm
        return model, transformed

    def objective(params: np.ndarray) -> float:
        model, _ = build(params)
        value = poisson_deviance(
            data[mask],
            np.clip(model[mask], 1.0e-12, None),
        )
        if use_nuisance_penalties:
            if log_morph:
                value += (float(params[0]) / 0.20) ** 2
                value += (float(params[1]) / 0.40) ** 2
            else:
                shift_width = max(shift_prior_bins * bin_width, 1.0e-12)
                smear_width = max(smear_prior_bins * bin_width, 1.0e-12)
                value += (float(params[0]) / shift_width) ** 2
                value += (float(params[1]) / smear_width) ** 2
                if asymmetric:
                    right_prior = max(0.5 * smear_width, 1.0e-12) if variable.branch == "theta" else smear_width
                    value += (float(params[2]) / right_prior) ** 2
                # endif
            # endif
        # endif
        return float(value)

    result = minimize(
        objective,
        x0,
        method="L-BFGS-B",
        bounds=bounds,
    )
    fitted = np.asarray(result.x, dtype=np.float64)
    model, transformed = build(fitted)
    nfit = int(np.count_nonzero(mask))
    parameter_count = 3 if asymmetric else 2

    if asymmetric:
        morph_label = "asymmetric-additive"
    elif is_lower_bounded_morph_variable(variable):
        morph_label = "lower-log-space"
    elif is_logit_morph_variable(variable):
        morph_label = "logit-space"
    elif is_upper_bounded_morph_variable(variable):
        morph_label = "upper-log-space"
    else:
        morph_label = "additive"
    # endif

    return FitResult(
        success=bool(result.success),
        message=(
            "full-range fit"
            if core_containment >= 0.999999
            else f"MC-defined {100.0 * core_containment:.1f}% signal-core fit"
        ),
        shift=float(fitted[0]),
        sigma_add=float(fitted[1]),
        sigma_right=float(fitted[2]) if asymmetric else math.nan,
        deviance=poisson_deviance(
            data[mask],
            np.clip(model[mask], 1.0e-12, None),
        ),
        ndf=max(1, nfit - parameter_count),
        data_total=data_total,
        model_counts=model,
        dvcs_component_counts=model,
        transformed_dvcs_shape=transformed,
        fit_mask=mask,
        morph_label=morph_label,
        excluded_data_counts=float(np.sum(data[~mask])),
        excluded_model_counts=float(np.sum(model[~mask])),
        excluded_excess_counts=float(
            np.sum(data[~mask]) - np.sum(model[~mask])
        ),
        excluded_excess_fraction=(
            float(np.sum(data[~mask]) - np.sum(model[~mask])) / data_total
            if data_total > 0.0
            else 0.0
        ),
    )


def draw_pi0_shape_canvas(
    output_path: Path,
    period: PeriodConfig,
    topology: TopologyConfig,
    data_histograms: Mapping[str, np.ndarray],
    mc_histograms: Mapping[str, np.ndarray],
    selected_counts: Mapping[str, int],
    log_y: bool,
    dpi: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 4, figsize=(20, 14))
    axes_flat = axes.ravel()

    for axis, variable in zip(axes_flat, PI0_SHAPE_VARIABLES):
        edges, centers = bin_geometry(variable)
        data_shape, data_err = normalize_density(data_histograms[variable.branch], variable)
        mc_shape, _ = normalize_density(mc_histograms[variable.branch], variable)
        axis.errorbar(centers, data_shape, yerr=data_err, fmt="o", markersize=3.0,
                      color="black", label=r"$e'p'\pi^0$ data")
        axis.stairs(mc_shape, edges, color="tab:blue", linewidth=2.0,
                    label=r"$e'p'\pi^0$ MC")
        axis.set_xlabel(variable.label)
        axis.set_ylabel("unit-area density")
        axis.grid(alpha=0.20)
        if log_y:
            axis.set_yscale("log")
            floor = positive_y_floor(data_shape, mc_shape)
            if floor is not None:
                axis.set_ylim(bottom=floor * 0.5)
            # endif
        # endif
    # endfor

    fig.suptitle(
        rf"$e'p'\pi^0$ exclusivity shapes after minimal preselection: {period.label}, {topology.label}" + "\n" +
        rf"selected entries: data={selected_counts['data']:,}, MC={selected_counts['mc']:,}",
        fontsize=17, y=0.975,
    )
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.915), ncol=2, frameon=False)
    fig.subplots_adjust(top=0.885, left=0.06, right=0.985, bottom=0.06, wspace=0.28, hspace=0.34)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def draw_pi0_fit_canvas(
    output_path: Path,
    period: PeriodConfig,
    topology: TopologyConfig,
    data_histograms: Mapping[str, np.ndarray],
    mc_histograms: Mapping[str, np.ndarray],
    selected_counts: Mapping[str, int],
    fit_results: Mapping[str, FitResult],
    log_y: bool,
    dpi: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes_flat = axes.ravel()

    for axis, variable in zip(axes_flat, PI0_VARIABLES):
        edges, centers = bin_geometry(variable)
        data = np.asarray(data_histograms[variable.branch], dtype=np.float64)
        raw_mc = normalized_shape(mc_histograms[variable.branch])
        result = fit_results[variable.branch]
        if result.fit_mask is not None and not np.all(result.fit_mask):
            for index in np.flatnonzero(~result.fit_mask):
                axis.axvspan(edges[index], edges[index + 1], color="0.88", alpha=0.42, linewidth=0)
            # endfor
        # endif
        axis.errorbar(centers, data, yerr=np.sqrt(np.clip(data, 0.0, None)), fmt="o",
                      markersize=3.0, color="black", label=r"$e'p'\pi^0$ data")
        if raw_mc is not None:
            axis.stairs(float(np.sum(data)) * raw_mc, edges, color="0.55", linestyle=":",
                        linewidth=1.8, label=r"raw $e'p'\pi^0$ MC")
        # endif
        if result.model_counts is not None:
            axis.stairs(result.model_counts, edges, color="tab:green", linewidth=2.2,
                        label="shifted/smeared MC fit")
        # endif
        axis.text(
            0.97, 0.95,
            rf"$\Delta={result.shift:.4g}$" + "\n" +
            (
                rf"$\sigma_L={result.sigma_add:.4g},\ "
                rf"\sigma_R={result.sigma_right:.4g}$ "
                rf"({result.morph_label})"
                if is_asymmetric_additive_variable(variable)
                else rf"$\sigma_{{add}}={result.sigma_add:.4g}$ "
                rf"({result.morph_label})"
            ) + "\n" +
            rf"$D/ndf={result.deviance:.1f}/{result.ndf}$" +
            (
                "\n" + f"outside-core excess={result.excluded_excess_counts:.0f}"
                if result.fit_mask is not None and not np.all(result.fit_mask)
                else ""
            ),
            transform=axis.transAxes, ha="right", va="top", fontsize=10,
        )
        axis.set_xlabel(variable.label)
        axis.set_ylabel("events / bin")
        axis.grid(alpha=0.20)
        if log_y:
            axis.set_yscale("log")
            floor = positive_y_floor(data, result.model_counts if result.model_counts is not None else np.asarray([]))
            if floor is not None:
                axis.set_ylim(bottom=max(0.5, floor * 0.5))
            # endif
        # endif
    # endfor

    fig.suptitle(
        rf"$e'p'\pi^0$ data--MC nuisance fits after minimal preselection: {period.label}, {topology.label}" + "\n" +
        rf"selected entries: data={selected_counts['data']:,}, MC={selected_counts['mc']:,}",
        fontsize=18, y=0.985,
    )
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.905), ncol=3, frameon=False)
    fig.subplots_adjust(top=0.82, left=0.06, right=0.985, bottom=0.08, wspace=0.28, hspace=0.30)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def process_pi0_period(
    period: PeriodConfig,
    topologies: Sequence[TopologyConfig],
    output_dir: Path,
    step_size: str,
    max_shift_bins: float,
    max_smear_bins: float,
    fit_min_counts: int,
    shift_prior_bins: float,
    smear_prior_bins: float,
    use_nuisance_penalties: bool,
    pi0_core_containment: float,
    log_y: bool,
    dpi: int,
) -> Tuple[List[Dict[str, object]], Dict[str, Dict[str, Tuple[float, float]]]]:
    data_hists, data_counts = fill_histograms_for_file(
        period.eppi0_data_file, topologies, step_size, False, PI0_SHAPE_VARIABLES
    )
    mc_hists, mc_counts = fill_histograms_for_file(
        period.eppi0_mc_file, topologies, step_size, False, PI0_SHAPE_VARIABLES
    )
    rows: List[Dict[str, object]] = []
    calibrations: Dict[str, Dict[str, Tuple[float, float]]] = {}

    for topology in topologies:
        selected = {"data": data_counts[topology.key], "mc": mc_counts[topology.key]}
        shape_path = output_dir / "pi0_channel" / "shape_comparisons" / f"eppi0_shapes_{period.key}_{topology.key.lower()}.png"
        draw_pi0_shape_canvas(shape_path, period, topology, data_hists[topology.key], mc_hists[topology.key], selected, log_y, dpi)
        log(f"Wrote {shape_path}")

        fit_results: Dict[str, FitResult] = {}
        for variable in PI0_VARIABLES:
            result = fit_single_template(
                data_hists[topology.key][variable.branch],
                mc_hists[topology.key][variable.branch],
                variable, topology, max_shift_bins, max_smear_bins, fit_min_counts,
                shift_prior_bins, smear_prior_bins, use_nuisance_penalties,
                pi0_core_containment,
            )
            fit_results[variable.branch] = result
            bin_width = (variable.xmax - variable.xmin) / variable.bins
            additive = not is_log_morph_variable(variable)
            rows.append({
                "period": period.key,
                "period_label": period.label,
                "topology": topology.key,
                "variable": variable.branch,
                "success": int(result.success),
                "message": result.message,
                "data_entries_in_range": result.data_total,
                "mc_entries_in_range": float(np.sum(mc_hists[topology.key][variable.branch])),
                "morph_type": result.morph_label,
                "nuisance_boundary_flags": ";".join(
                    nuisance_boundary_flags(
                        result,
                        variable,
                        max_shift_bins,
                        max_smear_bins,
                    )
                ),
                "core_containment": pi0_core_containment,
                "fit_bins_used": int(np.count_nonzero(result.fit_mask)) if result.fit_mask is not None else 0,
                "outside_core_data_counts": result.excluded_data_counts,
                "outside_core_model_counts": result.excluded_model_counts,
                "outside_core_excess_counts": result.excluded_excess_counts,
                "shift_or_log_shift": result.shift,
                "shift_bins": result.shift / bin_width if result.success and additive else math.nan,
                "sigma_or_log_sigma": result.sigma_add,
                "sigma_right": result.sigma_right,
                "sigma_bins": result.sigma_add / bin_width if result.success and additive else math.nan,
                "sigma_right_bins": result.sigma_right / bin_width if result.success and additive and math.isfinite(result.sigma_right) else math.nan,
                "poisson_deviance": result.deviance,
                "ndf": result.ndf,
                "deviance_per_ndf": result.deviance / result.ndf if result.success and result.ndf > 0 else math.nan,
            })
        # endfor

        calibrations[topology.key] = {
            branch: (result.shift, result.sigma_add)
            for branch, result in fit_results.items()
            if result.success and math.isfinite(result.shift)
        }

        fit_path = output_dir / "pi0_channel" / "template_fits" / f"eppi0_template_fit_{period.key}_{topology.key.lower()}.png"
        draw_pi0_fit_canvas(fit_path, period, topology, data_hists[topology.key], mc_hists[topology.key], selected, fit_results, log_y, dpi)
        log(f"Wrote {fit_path}")
    # endfor

    return rows, calibrations

def process_period(
    period: PeriodConfig,
    topologies: Sequence[TopologyConfig],
    output_dir: Path,
    step_size: str,
    max_shift_bins: float,
    max_smear_bins: float,
    fit_min_counts: int,
    fraction_variables: Sequence[str],
    shift_prior_bins: float,
    smear_prior_bins: float,
    use_nuisance_penalties: bool,
    dvcs_core_containment: float,
    dvcs_fraction_containment: float,
    outside_overshoot_penalty_weight: float,
    emiss2_mean_order_penalty_weight: float,
    log_y: bool,
    dpi: int,
    pi0_core_calibrations: Optional[Mapping[str, Mapping[str, Tuple[float, float]]]] = None,
    fraction_population_arrays: Optional[
        Mapping[str, Mapping[str, Mapping[str, np.ndarray]]]
    ] = None,
) -> List[Dict[str, object]]:
    log(f"Starting period {period.label}")

    # Minimal global preselection only. No hard exclusivity cuts are imposed.
    data_hists, data_counts = fill_histograms_for_file(
        period.data_file, topologies, step_size, apply_mx2_1_cut=False,
        variables=DVCS_SHAPE_VARIABLES,
    )
    dvcs_hists, dvcs_counts = fill_histograms_for_file(
        period.dvcs_mc_file, topologies, step_size, apply_mx2_1_cut=False,
        variables=DVCS_SHAPE_VARIABLES,
    )
    pi0_hists, pi0_counts = fill_histograms_for_file(
        period.pi0_as_dvcs_mc_file, topologies, step_size, apply_mx2_1_cut=False,
        variables=DVCS_SHAPE_VARIABLES,
    )

    rows: List[Dict[str, object]] = []

    for topology in topologies:
        selected = {
            "data": data_counts[topology.key],
            "dvcs_mc": dvcs_counts[topology.key],
            "pi0_mc": pi0_counts[topology.key],
        }
        shape_path = (
            output_dir / "dvcs_channel" / "shape_comparisons" /
            f"exclusivity_shapes_{period.key}_{topology.key.lower()}.png"
        )
        draw_shape_canvas(
            shape_path, period, topology, data_hists[topology.key],
            dvcs_hists[topology.key], pi0_hists[topology.key],
            selected, log_y, dpi, "minimal preselection",
        )
        log(f"Wrote {shape_path}")

        fit_data_hists = data_hists[topology.key]
        fit_dvcs_hists = dvcs_hists[topology.key]
        fit_pi0_hists = pi0_hists[topology.key]
        common_population_counts = {
            "data": math.nan,
            "dvcs_mc": math.nan,
            "pi0_mc": math.nan,
        }
        if fraction_population_arrays is not None:
            fit_data_hists, data_common_mask = common_population_histograms(
                fraction_population_arrays["data"][topology.key],
                fraction_variables,
                VARIABLES,
            )
            fit_dvcs_hists, dvcs_common_mask = common_population_histograms(
                fraction_population_arrays["signal"][topology.key],
                fraction_variables,
                VARIABLES,
            )
            fit_pi0_hists, pi0_common_mask = common_population_histograms(
                fraction_population_arrays["background"][topology.key],
                fraction_variables,
                VARIABLES,
            )
            common_population_counts = {
                "data": int(np.count_nonzero(data_common_mask)),
                "dvcs_mc": int(np.count_nonzero(dvcs_common_mask)),
                "pi0_mc": int(np.count_nonzero(pi0_common_mask)),
            }
        # endif

        shared_summary = fit_shared_two_templates(
            fit_data_hists,
            fit_dvcs_hists,
            fit_pi0_hists,
            topology,
            max_shift_bins, max_smear_bins, fit_min_counts,
            fraction_variables, shift_prior_bins, smear_prior_bins,
            use_nuisance_penalties,
            dvcs_core_containment,
            dvcs_fraction_containment,
            (pi0_core_calibrations or {}).get(topology.key, {}),
            outside_overshoot_penalty_weight,
            emiss2_mean_order_penalty_weight,
        )
        if not shared_summary.success or shared_summary.variable_results is None:
            raise RuntimeError(
                f"Shared fit failed for {period.label} {topology.label}: {shared_summary.message}"
            )
        # endif
        fit_results = shared_summary.variable_results

        individual_summaries: Dict[str, SharedFitSummary] = {}
        for branch in fraction_variables:
            individual_summaries[branch] = fit_shared_two_templates(
                fit_data_hists, fit_dvcs_hists,
                fit_pi0_hists, topology,
                max_shift_bins, max_smear_bins, fit_min_counts,
                [branch], shift_prior_bins, smear_prior_bins,
                use_nuisance_penalties,
                dvcs_core_containment,
                dvcs_fraction_containment,
                (pi0_core_calibrations or {}).get(topology.key, {}),
                outside_overshoot_penalty_weight,
                emiss2_mean_order_penalty_weight,
            )
        # endfor

        individual_fraction_values = np.asarray(
            [
                summary.f_pi0
                for summary in individual_summaries.values()
                if summary.success and math.isfinite(summary.f_pi0)
            ],
            dtype=np.float64,
        )
        if individual_fraction_values.size:
            individual_fraction_median = float(np.median(individual_fraction_values))
            individual_fraction_mad = float(
                np.median(
                    np.abs(
                        individual_fraction_values
                        - individual_fraction_median
                    )
                )
            )
            individual_fraction_min = float(np.min(individual_fraction_values))
            individual_fraction_max = float(np.max(individual_fraction_values))
        else:
            individual_fraction_median = math.nan
            individual_fraction_mad = math.nan
            individual_fraction_min = math.nan
            individual_fraction_max = math.nan
        # endif

        for variable in VARIABLES:
            result = fit_results[variable.branch]
            bin_width = (variable.xmax - variable.xmin) / variable.bins
            additive = not is_log_morph_variable(variable)
            rows.append({
                "period": period.key,
                "period_label": period.label,
                "topology": topology.key,
                "variable": variable.branch,
                "dvcs_core_containment": dvcs_core_containment,
                "dvcs_fraction_containment": dvcs_fraction_containment,
                "success": int(result.success),
                "message": result.message,
                "fit_role": result.message,
                "fraction_variables": ";".join(shared_summary.fraction_variables),
                "shared_f_pi0": shared_summary.f_pi0,
                "shared_f_pi0_err": shared_summary.f_pi0_err,
                "individual_f_pi0_median": individual_fraction_median,
                "individual_f_pi0_mad": individual_fraction_mad,
                "individual_f_pi0_min": individual_fraction_min,
                "individual_f_pi0_max": individual_fraction_max,
                "global_poisson_deviance": shared_summary.deviance,
                "global_ndf": shared_summary.ndf,
                "data_entries_in_range": result.data_total,
                "dvcs_mc_entries_in_range": float(np.sum(fit_dvcs_hists[variable.branch])),
                "pi0_mc_entries_in_range": float(np.sum(fit_pi0_hists[variable.branch])),
                "common_fraction_population_data": common_population_counts["data"],
                "common_fraction_population_dvcs_mc": common_population_counts["dvcs_mc"],
                "common_fraction_population_pi0_mc": common_population_counts["pi0_mc"],
                "morph_type": result.morph_label,
                "shift_or_log_shift": result.shift,
                "shift_err": result.shift_err,
                "shift_bins": result.shift / bin_width if result.success and additive else math.nan,
                "sigma_or_log_sigma": result.sigma_add,
                "sigma_right": result.sigma_right,
                "sigma_err": result.sigma_add_err,
                "sigma_right_err": result.sigma_right_err,
                "sigma_bins": result.sigma_add / bin_width if result.success and additive else math.nan,
                "fit_bins_used": int(np.count_nonzero(result.fit_mask)) if result.fit_mask is not None else 0,
                "outside_region_data_counts": result.excluded_data_counts,
                "outside_region_model_counts": result.excluded_model_counts,
                "outside_region_excess_counts": result.excluded_excess_counts,
                "outside_region_excess_fraction_of_data": result.excluded_excess_fraction,
                "variable_poisson_deviance": result.deviance,
                "variable_ndf": result.ndf,
                "variable_deviance_per_ndf": result.deviance / result.ndf if result.success and result.ndf > 0 else math.nan,
                **{
                    f"individual_f_pi0_{branch}": summary.f_pi0
                    for branch, summary in individual_summaries.items()
                },
                **{
                    f"individual_f_pi0_err_{branch}": summary.f_pi0_err
                    for branch, summary in individual_summaries.items()
                },
            })
        # endfor

        fit_path = (
            output_dir / "dvcs_channel" / "template_fits" /
            f"exclusivity_template_fit_{period.key}_{topology.key.lower()}.png"
        )
        draw_fit_canvas(
            fit_path, period, topology, data_hists[topology.key],
            dvcs_hists[topology.key], pi0_hists[topology.key],
            selected, fit_results, shared_summary, log_y, dpi,
        )
        log(f"Wrote {fit_path}")
    # endfor

    return rows


def write_results_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    # endif
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    # endwith



@dataclass
class IterativeCutStep:
    iteration: int
    variable: str
    score: float
    data_efficiency: float
    dvcs_mc_efficiency: float
    pi0_mc_efficiency: float
    f_pi0_before: float
    f_pi0_after: float
    f_pi0_raw_propagated: float
    f_pi0_diagnostic: float
    boundaries: Dict[str, Dict[str, Tuple[float, float]]]
    before_histograms: Dict[str, np.ndarray]
    after_histograms: Dict[str, np.ndarray]
    fit_result: FitResult
    calibration_source: str = "fitted nuisance"


def period_code(period: PeriodConfig) -> str:
    mapping = {
        "fa18_inb": "Fa18_Inb",
        "fa18_out": "Fa18_Out",
        "sp18_inb": "Sp18_Inb",
        "sp18_out": "Sp18_Out",
        "sp19_inb": "Sp19_Inb",
    }
    return mapping[period.key]


def cut_json_key(channel: str, period: PeriodConfig, topology: TopologyConfig) -> str:
    prefix = "DVCS" if channel == "dvcs" else "eppi0"
    return f"{prefix}_{period_code(period)}_{topology.key}"


def containment_window_from_shape(
    shape: np.ndarray,
    variable: VariableConfig,
    containment: float,
) -> Tuple[float, float, str]:
    """Return a one- or two-sided histogram-edge containment window."""
    values = np.asarray(shape, dtype=np.float64)
    total = float(np.sum(values))
    edges, _ = bin_geometry(variable)
    if total <= 0.0 or not math.isfinite(total):
        return variable.xmin, variable.xmax, "quantile_window"
    # endif

    cdf = np.cumsum(values) / total
    if is_lower_bounded_morph_variable(variable):
        upper_index = int(np.searchsorted(cdf, containment, side="left"))
        upper_index = max(0, min(upper_index, variable.bins - 1))
        return variable.xmin, float(edges[upper_index + 1]), "upper_quantile"
    # endif

    if is_upper_bounded_morph_variable(variable):
        lower_index = int(np.searchsorted(cdf, 1.0 - containment, side="left"))
        lower_index = max(0, min(lower_index, variable.bins - 1))
        return float(edges[lower_index]), variable.xmax, "quantile_window"
    # endif

    tail = 0.5 * (1.0 - containment)
    lower_index = int(np.searchsorted(cdf, tail, side="left"))
    upper_index = int(np.searchsorted(cdf, 1.0 - tail, side="left"))
    lower_index = max(0, min(lower_index, variable.bins - 1))
    upper_index = max(lower_index, min(upper_index, variable.bins - 1))
    return float(edges[lower_index]), float(edges[upper_index + 1]), "quantile_window"


def apply_window(values: np.ndarray, low: float, high: float) -> np.ndarray:
    values = np.asarray(values)
    return np.isfinite(values) & (values >= low) & (values <= high)


def common_in_range_mask(
    arrays: Mapping[str, np.ndarray],
    variable_branches: Sequence[str],
    variable_lookup: Mapping[str, VariableConfig],
    base_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Require one common event population across the selected variables."""

    if not arrays:
        return np.asarray([], dtype=bool)
    # endif

    first = np.asarray(next(iter(arrays.values())))
    mask = (
        np.ones(first.size, dtype=bool)
        if base_mask is None
        else np.asarray(base_mask, dtype=bool).copy()
    )
    for branch in variable_branches:
        variable = variable_lookup[branch]
        values = np.asarray(arrays[branch], dtype=np.float64)
        mask &= (
            np.isfinite(values)
            & (values >= variable.xmin)
            & (values < variable.xmax)
        )
    # endfor
    return mask


def common_population_histograms(
    arrays: Mapping[str, np.ndarray],
    variable_branches: Sequence[str],
    variables: Sequence[VariableConfig],
    base_mask: Optional[np.ndarray] = None,
) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    """Histogram all variables after a common driver-range requirement."""

    lookup = {variable.branch: variable for variable in variables}
    mask = common_in_range_mask(
        arrays,
        variable_branches,
        lookup,
        base_mask,
    )
    return arrays_to_histograms(arrays, mask, variables), mask


def arrays_to_histograms(
    arrays: Mapping[str, np.ndarray],
    mask: np.ndarray,
    variables: Sequence[VariableConfig],
) -> Dict[str, np.ndarray]:
    output: Dict[str, np.ndarray] = {}
    for variable in variables:
        values = np.asarray(arrays[variable.branch])[mask]
        values = values[np.isfinite(values)]
        output[variable.branch], _ = np.histogram(
            values,
            bins=variable.bins,
            range=(variable.xmin, variable.xmax),
        )
        output[variable.branch] = output[variable.branch].astype(np.float64)
    # endfor
    return output


def load_selected_arrays_uproot(
    path: str,
    topologies: Sequence[TopologyConfig],
    step_size: str,
    variables: Sequence[VariableConfig],
) -> Dict[str, Dict[str, np.ndarray]]:
    """Load minimally selected fit-variable arrays once for iterative cuts."""
    require_input_file(path)
    resolved = resolve_variable_branches(path, variables)
    expressions = [
        "detector1", "detector2", "t1", "open_angle_ep2",
        "e_phi", "p1_phi", "p2_phi",
    ] + sorted(set(resolved.values()))
    pieces: Dict[str, Dict[str, List[np.ndarray]]] = {
        topology.key: {variable.branch: [] for variable in variables}
        for topology in topologies
    }

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=step_size,
        library="np",
    ):
        detector1 = np.asarray(arrays["detector1"])
        detector2 = np.asarray(arrays["detector2"])
        t1 = np.asarray(arrays["t1"], dtype=np.float64)
        open_angle = np.asarray(arrays["open_angle_ep2"], dtype=np.float64)
        base_mask = (
            np.isfinite(t1)
            & np.isfinite(open_angle)
            & ((-t1) < T1_ABS_MAX)
            & (open_angle > OPEN_ANGLE_MIN_DEG)
        )
        if REQUIRE_DISTINCT_FD_SECTORS:
            e_sector = fd_sector_from_phi(np.asarray(arrays["e_phi"], dtype=np.float64))
            p_sector = fd_sector_from_phi(np.asarray(arrays["p1_phi"], dtype=np.float64))
            g_sector = fd_sector_from_phi(np.asarray(arrays["p2_phi"], dtype=np.float64))
        # endif

        for topology in topologies:
            mask = (
                base_mask
                & (detector1 == topology.detector1)
                & (detector2 == topology.detector2)
            )
            if REQUIRE_DISTINCT_FD_SECTORS:
                mask &= e_sector >= 1
                if topology.detector1 == 1:
                    mask &= (p_sector >= 1) & (p_sector != e_sector)
                # endif
                if topology.detector2 == 1:
                    mask &= (g_sector >= 1) & (g_sector != e_sector)
                # endif
                if topology.detector1 == 1 and topology.detector2 == 1:
                    mask &= p_sector != g_sector
                # endif
            # endif
            if not np.any(mask):
                continue
            # endif
            for variable in variables:
                values = np.asarray(
                    arrays[resolved[variable.branch]],
                    dtype=np.float32,
                )[mask]
                pieces[topology.key][variable.branch].append(values)
            # endfor
        # endfor
    # endfor

    result: Dict[str, Dict[str, np.ndarray]] = {}
    for topology in topologies:
        result[topology.key] = {}
        for variable in variables:
            chunks = pieces[topology.key][variable.branch]
            result[topology.key][variable.branch] = (
                np.concatenate(chunks)
                if chunks
                else np.asarray([], dtype=np.float32)
            )
        # endfor
    # endfor
    return result


def load_selected_arrays_pyroot(
    path: str,
    topologies: Sequence[TopologyConfig],
    step_size: str,
    variables: Sequence[VariableConfig],
) -> Dict[str, Dict[str, np.ndarray]]:
    del step_size
    require_input_file(path)
    resolved = resolve_variable_branches(path, variables)
    pieces: Dict[str, Dict[str, List[float]]] = {
        topology.key: {variable.branch: [] for variable in variables}
        for topology in topologies
    }
    root_file = ROOT.TFile.Open(path, "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    # endif
    tree = root_file.Get(TREE_NAME)
    if tree is None:
        root_file.Close()
        raise KeyError(f"Tree '{TREE_NAME}' not found in {path}")
    # endif
    try:
        for event in tree:
            detector1 = int(getattr(event, "detector1"))
            detector2 = int(getattr(event, "detector2"))
            topology = next(
                (
                    candidate for candidate in topologies
                    if candidate.detector1 == detector1
                    and candidate.detector2 == detector2
                ),
                None,
            )
            if topology is None:
                continue
            # endif
            t1 = float(getattr(event, "t1"))
            opening = float(getattr(event, "open_angle_ep2"))
            if not math.isfinite(t1) or not math.isfinite(opening):
                continue
            # endif
            if (-t1) >= T1_ABS_MAX or opening <= OPEN_ANGLE_MIN_DEG:
                continue
            # endif
            if REQUIRE_DISTINCT_FD_SECTORS:
                sectors = [int(fd_sector_from_phi(np.asarray([float(getattr(event, "e_phi"))]))[0])]
                if topology.detector1 == 1:
                    sectors.append(int(fd_sector_from_phi(np.asarray([float(getattr(event, "p1_phi"))]))[0]))
                # endif
                if topology.detector2 == 1:
                    sectors.append(int(fd_sector_from_phi(np.asarray([float(getattr(event, "p2_phi"))]))[0]))
                # endif
                if any(value < 1 for value in sectors) or len(sectors) != len(set(sectors)):
                    continue
                # endif
            # endif
            for variable in variables:
                pieces[topology.key][variable.branch].append(
                    float(getattr(event, resolved[variable.branch]))
                )
            # endfor
        # endfor
    finally:
        root_file.Close()
    # endtry
    return {
        topology.key: {
            variable.branch: np.asarray(
                pieces[topology.key][variable.branch],
                dtype=np.float32,
            )
            for variable in variables
        }
        for topology in topologies
    }


def load_selected_arrays_for_file(
    path: str,
    topologies: Sequence[TopologyConfig],
    step_size: str,
    variables: Sequence[VariableConfig],
) -> Dict[str, Dict[str, np.ndarray]]:
    if io_backend() == "uproot":
        return load_selected_arrays_uproot(path, topologies, step_size, variables)
    # endif
    return load_selected_arrays_pyroot(path, topologies, step_size, variables)


def distribution_moments(
    shape: np.ndarray,
    variable: VariableConfig,
) -> Tuple[float, float]:
    _, centers = bin_geometry(variable)
    weights = np.asarray(shape, dtype=np.float64)
    total = float(np.sum(weights))
    if total <= 0.0:
        return 0.5 * (variable.xmin + variable.xmax), 0.0
    # endif
    mean = float(np.sum(weights * centers) / total)
    variance = float(np.sum(weights * (centers - mean) ** 2) / total)
    return mean, math.sqrt(max(0.0, variance))


def draw_iterative_cut_canvas(
    output_path: Path,
    period: PeriodConfig,
    topology: TopologyConfig,
    channel: str,
    steps: Sequence[IterativeCutStep],
    summary: bool,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(20.0, 10.0))
    flat_axes = axes.ravel()
    variable_lookup = {
        variable.branch: variable
        for variable in (VARIABLES if channel == "dvcs" else PI0_VARIABLES)
    }
    ordered_steps = list(steps)
    if summary:
        template_order = [
            variable.branch
            for variable in (VARIABLES if channel == "dvcs" else PI0_VARIABLES)
        ]
        step_by_variable = {step.variable: step for step in steps}
        ordered_steps = [
            step_by_variable[branch]
            for branch in template_order
            if branch in step_by_variable
        ]
    # endif

    for axis, step in zip(flat_axes, ordered_steps):
        variable = variable_lookup[step.variable]
        edges, centers = bin_geometry(variable)
        histograms = step.after_histograms if summary else step.before_histograms
        data = histograms["data"]
        signal = histograms["signal"]
        background = histograms.get("background")
        axis.errorbar(
            centers,
            data,
            yerr=np.sqrt(np.clip(data, 0.0, None)),
            fmt="o",
            markersize=2.8,
            linewidth=0.8,
            color="black",
            label="data",
            zorder=5,
        )
        if signal is not None:
            signal_shape = normalized_shape(signal)
            if signal_shape is not None:
                axis.stairs(
                    float(np.sum(data)) * signal_shape,
                    edges,
                    color="0.55",
                    linestyle=":",
                    linewidth=1.4,
                    label="raw signal MC",
                )
            # endif
        # endif
        if not summary and step.fit_result.success:
            if channel == "dvcs":
                axis.stairs(
                    step.fit_result.dvcs_component_counts,
                    edges,
                    color="tab:blue",
                    linewidth=1.6,
                    label="fitted DVCS",
                )
                axis.stairs(
                    step.fit_result.pi0_component_counts,
                    edges,
                    color="tab:red",
                    linewidth=1.6,
                    label=r"fitted $e\pi^0$",
                )
                axis.stairs(
                    step.fit_result.model_counts,
                    edges,
                    color="tab:green",
                    linestyle="--",
                    linewidth=2.0,
                    label="total fit",
                )
            else:
                axis.stairs(
                    step.fit_result.model_counts,
                    edges,
                    color="tab:green",
                    linewidth=2.0,
                    label="morphed signal MC",
                )
            # endif
        elif summary and background is not None:
            background_shape = normalized_shape(background)
            if background_shape is not None:
                axis.stairs(
                    float(np.sum(data)) * background_shape,
                    edges,
                    color="tab:red",
                    linewidth=1.4,
                    label=r"$e\pi^0$ MC",
                )
            # endif
        # endif

        low, high = step.boundaries["nominal"]["data"]
        mc_low, mc_high = step.boundaries["nominal"]["mc"]
        if low > variable.xmin + 1.0e-12:
            axis.axvline(
                low,
                color="tab:purple",
                linewidth=1.8,
                label="data boundary",
            )
        # endif
        if high < variable.xmax - 1.0e-12:
            axis.axvline(high, color="tab:purple", linewidth=1.8)
        # endif
        if mc_low > variable.xmin + 1.0e-12:
            axis.axvline(
                mc_low,
                color="tab:purple",
                linewidth=0.9,
                linestyle="--",
                label="MC boundary",
            )
        # endif
        if mc_high < variable.xmax - 1.0e-12:
            axis.axvline(
                mc_high,
                color="tab:purple",
                linewidth=0.9,
                linestyle="--",
            )
        # endif
        annotation = (
            f"step {step.iteration + 1}: {step.variable}\n"
            f"calibration={step.calibration_source}\n"
            f"score={step.score:.3f}\n"
            f"data cut=[{low:.5g}, {high:.5g}]\n"
            f"MC cut=[{step.boundaries['nominal']['mc'][0]:.5g}, "
            f"{step.boundaries['nominal']['mc'][1]:.5g}]\n"
            f"signal eff={100.0 * step.dvcs_mc_efficiency:.1f}%"
        )
        if channel == "dvcs":
            annotation += (
                f"\npi0 eff={100.0 * step.pi0_mc_efficiency:.1f}%"
                f"\nf_pi0: {step.f_pi0_before:.3f}"
                f" -> {step.f_pi0_after:.3f}"
            )
        # endif
        axis.text(
            0.98,
            0.96,
            annotation,
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=8.5,
            bbox=dict(facecolor="white", alpha=0.80, edgecolor="none"),
        )
        axis.set_xlim(variable.xmin, variable.xmax)
        axis.set_ylim(bottom=0.0)
        axis.set_xlabel(variable.label)
        axis.set_ylabel("events / bin")
        axis.grid(axis="y", alpha=0.22)
    # endfor

    title_channel = r"$e'p'\gamma$" if channel == "dvcs" else r"$e'p'\pi^0$"
    title_kind = "final iterative-cut summary" if summary else "iterative cut development"
    fig.suptitle(
        f"{title_channel} {title_kind}: {period.label}, {topology.label}\n"
        "automatic order; nominal signal containment = 95%",
        fontsize=17,
        y=0.985,
    )
    handles, labels = flat_axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.91),
            ncol=min(5, len(labels)),
            frameon=False,
        )
    # endif
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.86))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def make_cut_stats(
    shape: np.ndarray,
    variable: VariableConfig,
    low: float,
    high: float,
    containment: float,
    mode: str,
    iteration: int,
    score: float,
) -> Dict[str, object]:
    mean, std = distribution_moments(shape, variable)
    return {
        "mean": mean,
        "std": std,
        "cut_low": low,
        "cut_high": high,
        "quantile": containment,
        "mode": mode,
        "iteration": iteration + 1,
        "discrimination_score": score,
    }



def nuisance_boundary_flags(
    result: FitResult,
    variable: VariableConfig,
    max_shift_bins: float,
    max_smear_bins: float,
) -> Tuple[str, ...]:
    """Identify nuisance parameters that are effectively on a fit boundary."""

    flags: List[str] = []
    if not result.success:
        flags.append("fit_failed")
        return tuple(flags)
    # endif

    bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
    if is_log_morph_variable(variable):
        if abs(result.shift) >= 0.695:
            flags.append("shift_bound")
        # endif
        if result.sigma_add >= 0.995:
            flags.append("smear_bound")
        # endif
    else:
        if abs(result.shift) >= 0.98 * max_shift_bins * bin_width:
            flags.append("shift_bound")
        # endif
        if result.sigma_add >= 0.98 * max_smear_bins * bin_width:
            flags.append("smear_bound")
        # endif
        if (
            is_asymmetric_additive_variable(variable)
            and math.isfinite(result.sigma_right)
            and result.sigma_right <= 1.0e-8
        ):
            flags.append("right_width_zero")
        # endif
    # endif
    return tuple(flags)


def iterative_signal_shape(
    raw_signal_shape: np.ndarray,
    variable: VariableConfig,
    fitted_result: FitResult,
) -> Tuple[np.ndarray, str]:
    """Use the template-fit nuisance calibration for iterative cut construction.

    The raw reconstructed signal template is used only when the fitted
    calibration is genuinely unusable: missing/nonfinite parameters or a
    failed template transformation. Reduced deviance is intentionally not used
    as a rejection criterion because the very large event samples make even
    visually adequate fits yield large D/ndf values.
    """

    raw = np.asarray(raw_signal_shape, dtype=np.float64)
    usable = (
        fitted_result.success
        and math.isfinite(fitted_result.shift)
        and math.isfinite(fitted_result.sigma_add)
        and (
            not is_asymmetric_additive_variable(variable)
            or math.isfinite(fitted_result.sigma_right)
        )
    )
    if not usable:
        return raw, "raw-MC fallback"
    # endif

    transformed = transform_dvcs_shape(
        raw,
        variable,
        fitted_result.shift,
        fitted_result.sigma_add,
        (
            fitted_result.sigma_right
            if is_asymmetric_additive_variable(variable)
            else None
        ),
    )
    if transformed is None:
        return raw, "raw-MC fallback"
    # endif

    return transformed, "fitted nuisance"


def build_frozen_two_template_result(
    data_histogram: np.ndarray,
    signal_shape: np.ndarray,
    background_shape: np.ndarray,
    fraction: float,
    frozen_result: FitResult,
) -> FitResult:
    """Build a display model using frozen nuisance calibration and a fixed fraction."""

    data = np.asarray(data_histogram, dtype=np.float64)
    signal = np.asarray(signal_shape, dtype=np.float64)
    background = np.asarray(background_shape, dtype=np.float64)
    total_shape = np.clip(
        (1.0 - fraction) * signal + fraction * background,
        0.0,
        None,
    )
    normalization = float(np.sum(total_shape))
    data_total = float(np.sum(data))
    if normalization <= 0.0 or data_total <= 0.0:
        return FitResult(False, "invalid frozen-template model", data_total=data_total)
    # endif

    scale = data_total / normalization
    model = scale * total_shape
    dvcs_component = scale * (1.0 - fraction) * signal
    pi0_component = scale * fraction * background
    return FitResult(
        success=True,
        message="frozen nuisance calibration",
        f_pi0=fraction,
        shift=frozen_result.shift,
        sigma_add=frozen_result.sigma_add,
        sigma_right=frozen_result.sigma_right,
        deviance=poisson_deviance(data, model),
        ndf=max(1, data.size - 1),
        data_total=data_total,
        model_counts=model,
        dvcs_component_counts=dvcs_component,
        pi0_component_counts=pi0_component,
        transformed_dvcs_shape=signal,
        fit_mask=np.ones(data.size, dtype=bool),
        morph_label=frozen_result.morph_label,
    )


def diagnostic_fraction_with_frozen_shapes(
    data_histograms: Mapping[str, np.ndarray],
    signal_shapes: Mapping[str, np.ndarray],
    background_shapes: Mapping[str, np.ndarray],
    variables: Sequence[str],
    propagated_fraction: float,
) -> float:
    """Profile only f_pi0 with frozen shapes and a stabilizing propagated-f prior."""

    usable = [
        branch
        for branch in variables
        if branch in data_histograms
        and branch in signal_shapes
        and branch in background_shapes
        and float(np.sum(data_histograms[branch])) > 0.0
    ]
    if not usable:
        return math.nan
    # endif

    prior_width = max(0.05, 0.35 * max(propagated_fraction, 0.05))

    def objective(candidate: float) -> float:
        value = 0.0
        for branch in usable:
            data = np.asarray(data_histograms[branch], dtype=np.float64)
            total_shape = np.clip(
                (1.0 - candidate) * signal_shapes[branch]
                + candidate * background_shapes[branch],
                0.0,
                None,
            )
            normalization = float(np.sum(total_shape))
            if normalization <= 0.0:
                return 1.0e100
            # endif
            model = float(np.sum(data)) * total_shape / normalization
            value += 0.5 * poisson_deviance(data, model)
        # endfor
        value += 0.5 * ((candidate - propagated_fraction) / prior_width) ** 2
        return value

    result = minimize_scalar(
        objective,
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1.0e-5, "maxiter": 200},
    )
    if not result.success or not math.isfinite(float(result.x)):
        return math.nan
    # endif
    return float(result.x)


def run_dvcs_iterative_cuts(
    period: PeriodConfig,
    topology: TopologyConfig,
    arrays: Mapping[str, Mapping[str, np.ndarray]],
    max_shift_bins: float,
    max_smear_bins: float,
    min_counts: int,
    fraction_variables: Sequence[str],
    shift_prior_bins: float,
    smear_prior_bins: float,
    use_nuisance_penalties: bool,
    core_containment: float,
    fraction_containment: float,
    pi0_calibration: Mapping[str, Tuple[float, float]],
    containments: Mapping[str, float],
    outside_overshoot_penalty_weight: float,
    emiss2_mean_order_penalty_weight: float,
) -> Tuple[List[IterativeCutStep], Dict[str, Dict[str, Dict[str, object]]]]:
    """Develop matched-containment data and MC cuts with one frozen initial calibration and propagated contamination."""

    data_arrays = arrays["data"]
    signal_arrays = arrays["signal"]
    background_arrays = arrays["background"]
    masks = {
        "data": np.ones(len(next(iter(data_arrays.values()))), dtype=bool),
        "signal": np.ones(len(next(iter(signal_arrays.values()))), dtype=bool),
        "background": np.ones(len(next(iter(background_arrays.values()))), dtype=bool),
    }

    initial_data_hists, initial_data_common = common_population_histograms(
        data_arrays,
        fraction_variables,
        VARIABLES,
        masks["data"],
    )
    initial_signal_hists, initial_signal_common = common_population_histograms(
        signal_arrays,
        fraction_variables,
        VARIABLES,
        masks["signal"],
    )
    initial_background_hists, initial_background_common = common_population_histograms(
        background_arrays,
        fraction_variables,
        VARIABLES,
        masks["background"],
    )
    initial_summary = fit_shared_two_templates(
        initial_data_hists,
        initial_signal_hists,
        initial_background_hists,
        topology,
        max_shift_bins,
        max_smear_bins,
        min_counts,
        fraction_variables,
        shift_prior_bins,
        smear_prior_bins,
        use_nuisance_penalties,
        core_containment,
        fraction_containment,
        pi0_calibration,
        outside_overshoot_penalty_weight,
        emiss2_mean_order_penalty_weight,
    )
    if not initial_summary.success or initial_summary.variable_results is None:
        raise RuntimeError(
            f"Initial iterative DVCS calibration failed for "
            f"{period.label} {topology.label}: {initial_summary.message}"
        )
    # endif

    frozen_results = initial_summary.variable_results
    propagated_fraction = float(initial_summary.f_pi0)
    remaining = [variable.branch for variable in VARIABLES]
    steps: List[IterativeCutStep] = []
    json_variants: Dict[str, Dict[str, Dict[str, object]]] = {
        name: {"data": {}, "mc": {}}
        for name in containments
    }

    for iteration in range(len(VARIABLES)):
        data_hists = arrays_to_histograms(data_arrays, masks["data"], VARIABLES)
        signal_hists = arrays_to_histograms(signal_arrays, masks["signal"], VARIABLES)
        background_hists = arrays_to_histograms(
            background_arrays,
            masks["background"],
            VARIABLES,
        )

        current_signal_shapes: Dict[str, np.ndarray] = {}
        current_background_shapes: Dict[str, np.ndarray] = {}
        candidates = []

        for variable in VARIABLES:
            branch = variable.branch
            if branch not in remaining:
                continue
            # endif

            frozen = frozen_results[branch]
            raw_signal_shape = normalized_shape(signal_hists[branch])
            background_shape = normalized_shape(background_hists[branch])
            if (
                not frozen.success
                or raw_signal_shape is None
                or background_shape is None
            ):
                continue
            # endif

            transformed_signal, calibration_source = iterative_signal_shape(
                raw_signal_shape,
                variable,
                frozen,
            )

            current_signal_shapes[branch] = transformed_signal
            current_background_shapes[branch] = background_shape

            boundaries: Dict[str, Dict[str, Tuple[float, float]]] = {}
            # Preserve the same signal survival probability in data and MC:
            # data boundaries come from the morphed DVCS signal template, while
            # reconstructed-MC boundaries come from the raw DVCS signal template.
            for name, containment in containments.items():
                data_low, data_high, _ = containment_window_from_shape(
                    transformed_signal,
                    variable,
                    containment,
                )
                mc_low, mc_high, _ = containment_window_from_shape(
                    raw_signal_shape,
                    variable,
                    containment,
                )
                boundaries[name] = {
                    "data": (data_low, data_high),
                    "mc": (mc_low, mc_high),
                }
            # endfor

            nominal_mc_low, nominal_mc_high = boundaries["nominal"]["mc"]
            signal_values = signal_arrays[branch][masks["signal"]]
            background_values = background_arrays[branch][masks["background"]]
            signal_eff = (
                float(np.mean(apply_window(signal_values, nominal_mc_low, nominal_mc_high)))
                if signal_values.size
                else 0.0
            )
            background_eff = (
                float(np.mean(apply_window(background_values, nominal_mc_low, nominal_mc_high)))
                if background_values.size
                else 1.0
            )

            nominal_data_low, nominal_data_high = boundaries["nominal"]["data"]
            data_values = data_arrays[branch][masks["data"]]
            data_eff = (
                float(np.mean(apply_window(data_values, nominal_data_low, nominal_data_high)))
                if data_values.size
                else 0.0
            )

            # Higher score means stronger pi0 rejection at comparable signal retention.
            score = signal_eff - background_eff
            candidates.append(
                (
                    score,
                    variable,
                    boundaries,
                    data_eff,
                    signal_eff,
                    background_eff,
                    transformed_signal,
                    raw_signal_shape,
                    background_shape,
                    calibration_source,
                )
            )
        # endfor

        if not candidates:
            raise RuntimeError(
                f"No valid iterative cut candidate for {period.label} "
                f"{topology.label} at step {iteration + 1}"
            )
        # endif

        candidates.sort(key=lambda item: (-item[0], item[1].branch))
        (
            score,
            selected_variable,
            boundaries,
            data_eff,
            signal_eff,
            background_eff,
            transformed_signal,
            raw_signal_shape,
            background_shape,
            calibration_source,
        ) = candidates[0]
        branch = selected_variable.branch
        frozen = frozen_results[branch]

        before_histograms = {
            "data": data_hists[branch],
            "signal": signal_hists[branch],
            "background": background_hists[branch],
        }
        display_result = build_frozen_two_template_result(
            data_hists[branch],
            transformed_signal,
            background_shape,
            propagated_fraction,
            frozen,
        )

        denominator = (
            (1.0 - propagated_fraction) * signal_eff
            + propagated_fraction * background_eff
        )
        raw_propagated = (
            propagated_fraction * background_eff / denominator
            if denominator > 0.0
            else propagated_fraction
        )
        # The authoritative reported contamination is monotonic. A late weak
        # variable is not allowed to manufacture an increase in pi0 contamination.
        propagated_after = min(propagated_fraction, raw_propagated)

        nominal_data_low, nominal_data_high = boundaries["nominal"]["data"]
        nominal_mc_low, nominal_mc_high = boundaries["nominal"]["mc"]
        masks["data"] &= apply_window(
            data_arrays[branch],
            nominal_data_low,
            nominal_data_high,
        )
        masks["signal"] &= apply_window(
            signal_arrays[branch],
            nominal_mc_low,
            nominal_mc_high,
        )
        masks["background"] &= apply_window(
            background_arrays[branch],
            nominal_mc_low,
            nominal_mc_high,
        )

        after_histograms = {
            "data": arrays_to_histograms(
                data_arrays, masks["data"], [selected_variable]
            )[branch],
            "signal": arrays_to_histograms(
                signal_arrays, masks["signal"], [selected_variable]
            )[branch],
            "background": arrays_to_histograms(
                background_arrays, masks["background"], [selected_variable]
            )[branch],
        }

        remaining_after = [
            candidate for candidate in remaining if candidate != branch
        ]
        diagnostic_fraction = math.nan

        steps.append(
            IterativeCutStep(
                iteration=iteration,
                variable=branch,
                score=score,
                data_efficiency=data_eff,
                dvcs_mc_efficiency=signal_eff,
                pi0_mc_efficiency=background_eff,
                f_pi0_before=propagated_fraction,
                f_pi0_after=propagated_after,
                f_pi0_raw_propagated=raw_propagated,
                f_pi0_diagnostic=diagnostic_fraction,
                boundaries=boundaries,
                before_histograms=before_histograms,
                after_histograms=after_histograms,
                fit_result=display_result,
                calibration_source=calibration_source,
            )
        )

        for name, containment in containments.items():
            data_low, data_high = boundaries[name]["data"]
            mc_low, mc_high = boundaries[name]["mc"]
            mode = (
                "upper_quantile"
                if is_lower_bounded_morph_variable(selected_variable)
                else "quantile_window"
            )
            json_variants[name]["data"][branch] = make_cut_stats(
                transformed_signal,
                selected_variable,
                data_low,
                data_high,
                containment,
                mode,
                iteration,
                score,
            )
            json_variants[name]["mc"][branch] = make_cut_stats(
                raw_signal_shape,
                selected_variable,
                mc_low,
                mc_high,
                containment,
                mode,
                iteration,
                score,
            )
        # endfor

        propagated_fraction = propagated_after
        remaining.remove(branch)
    # endfor

    return steps, json_variants

def run_pi0_iterative_cuts(
    period: PeriodConfig,
    topology: TopologyConfig,
    arrays: Mapping[str, Mapping[str, np.ndarray]],
    dvcs_order: Sequence[str],
    max_shift_bins: float,
    max_smear_bins: float,
    min_counts: int,
    shift_prior_bins: float,
    smear_prior_bins: float,
    use_nuisance_penalties: bool,
    core_containment: float,
    containments: Mapping[str, float],
) -> Tuple[List[IterativeCutStep], Dict[str, Dict[str, Dict[str, object]]]]:
    """Apply the DVCS-selected order with one frozen direct-pi0 calibration."""

    data_arrays = arrays["data"]
    signal_arrays = arrays["signal"]
    masks = {
        "data": np.ones(len(next(iter(data_arrays.values()))), dtype=bool),
        "signal": np.ones(len(next(iter(signal_arrays.values()))), dtype=bool),
    }
    pi0_lookup = {variable.branch: variable for variable in PI0_VARIABLES}
    mapped_order = [
        "theta_pi0_pi0" if branch == "theta_gamma_gamma" else branch
        for branch in dvcs_order
    ]

    initial_data_hists = arrays_to_histograms(
        data_arrays,
        masks["data"],
        PI0_VARIABLES,
    )
    initial_signal_hists = arrays_to_histograms(
        signal_arrays,
        masks["signal"],
        PI0_VARIABLES,
    )
    frozen_results: Dict[str, FitResult] = {}
    for variable in PI0_VARIABLES:
        result = fit_single_template(
            initial_data_hists[variable.branch],
            initial_signal_hists[variable.branch],
            variable,
            topology,
            max_shift_bins,
            max_smear_bins,
            min_counts,
            shift_prior_bins,
            smear_prior_bins,
            use_nuisance_penalties,
            core_containment,
        )
        # Finite returned parameters are usable even when L-BFGS-B reports a
        # boundary/status warning. Fall back to the raw template otherwise.
        if (
            result.transformed_dvcs_shape is None
            or not math.isfinite(result.shift)
            or not math.isfinite(result.sigma_add)
        ):
            raw_shape = normalized_shape(initial_signal_hists[variable.branch])
            result = FitResult(
                success=raw_shape is not None,
                message="raw-template fallback",
                shift=0.0,
                sigma_add=0.0,
                sigma_right=0.0 if is_asymmetric_additive_variable(variable) else math.nan,
                data_total=float(np.sum(initial_data_hists[variable.branch])),
                transformed_dvcs_shape=raw_shape,
                fit_mask=np.ones(variable.bins, dtype=bool),
                morph_label=(
                    "asymmetric-additive"
                    if is_asymmetric_additive_variable(variable)
                    else "lower-log-space"
                    if is_lower_bounded_morph_variable(variable)
                    else "logit-space"
                    if is_logit_morph_variable(variable)
                    else "upper-log-space"
                    if is_upper_bounded_morph_variable(variable)
                    else "additive"
                ),
            )
        # endif
        if (
            result.transformed_dvcs_shape is not None
            and math.isfinite(result.shift)
            and math.isfinite(result.sigma_add)
        ):
            # A finite L-BFGS-B boundary solution is usable even when SciPy's
            # status flag is false.
            result.success = True
        # endif
        frozen_results[variable.branch] = result
    # endfor

    steps: List[IterativeCutStep] = []
    json_variants: Dict[str, Dict[str, Dict[str, object]]] = {
        name: {"data": {}, "mc": {}}
        for name in containments
    }

    for iteration, branch in enumerate(mapped_order):
        variable = pi0_lookup[branch]
        data_hist = arrays_to_histograms(
            data_arrays, masks["data"], [variable]
        )[branch]
        signal_hist = arrays_to_histograms(
            signal_arrays, masks["signal"], [variable]
        )[branch]
        raw_shape = normalized_shape(signal_hist)
        frozen = frozen_results[branch]
        if raw_shape is None:
            raise RuntimeError(
                f"Empty direct-pi0 signal template for {period.label} "
                f"{topology.label} {branch}"
            )
        # endif
        if not frozen.success:
            frozen = FitResult(
                success=True,
                message="raw-template fallback",
                shift=0.0,
                sigma_add=0.0,
                sigma_right=(
                    0.0 if is_asymmetric_additive_variable(variable) else math.nan
                ),
                data_total=float(np.sum(data_hist)),
                transformed_dvcs_shape=raw_shape,
                fit_mask=np.ones(variable.bins, dtype=bool),
                morph_label=(
                    "asymmetric-additive"
                    if is_asymmetric_additive_variable(variable)
                    else "lower-log-space"
                    if is_lower_bounded_morph_variable(variable)
                    else "logit-space"
                    if is_logit_morph_variable(variable)
                    else "upper-log-space"
                    if is_upper_bounded_morph_variable(variable)
                    else "additive"
                ),
            )
        # endif

        transformed = transform_dvcs_shape(
            raw_shape,
            variable,
            frozen.shift,
            frozen.sigma_add,
            (
                frozen.sigma_right
                if is_asymmetric_additive_variable(variable)
                and math.isfinite(frozen.sigma_right)
                else None
            ),
        )
        if transformed is None:
            transformed = raw_shape
        # endif

        boundaries: Dict[str, Dict[str, Tuple[float, float]]] = {}
        # Use matched containment: morphed-template boundaries for data and
        # raw reconstructed-MC boundaries for MC.
        for name, containment in containments.items():
            data_low, data_high, _ = containment_window_from_shape(
                transformed,
                variable,
                containment,
            )
            mc_low, mc_high, _ = containment_window_from_shape(
                raw_shape,
                variable,
                containment,
            )
            boundaries[name] = {
                "data": (data_low, data_high),
                "mc": (mc_low, mc_high),
            }
        # endfor

        data_low, data_high = boundaries["nominal"]["data"]
        mc_low, mc_high = boundaries["nominal"]["mc"]
        current_data_values = data_arrays[branch][masks["data"]]
        current_signal_values = signal_arrays[branch][masks["signal"]]
        data_eff = (
            float(np.mean(apply_window(current_data_values, data_low, data_high)))
            if current_data_values.size
            else 0.0
        )
        signal_eff = (
            float(np.mean(apply_window(current_signal_values, mc_low, mc_high)))
            if current_signal_values.size
            else 0.0
        )

        display_model = FitResult(
            success=True,
            message="frozen direct-pi0 calibration",
            shift=frozen.shift,
            sigma_add=frozen.sigma_add,
            sigma_right=frozen.sigma_right,
            data_total=float(np.sum(data_hist)),
            model_counts=float(np.sum(data_hist)) * transformed,
            dvcs_component_counts=float(np.sum(data_hist)) * transformed,
            transformed_dvcs_shape=transformed,
            fit_mask=np.ones(variable.bins, dtype=bool),
            morph_label=frozen.morph_label,
        )

        masks["data"] &= apply_window(data_arrays[branch], data_low, data_high)
        masks["signal"] &= apply_window(signal_arrays[branch], mc_low, mc_high)
        after_histograms = {
            "data": arrays_to_histograms(
                data_arrays, masks["data"], [variable]
            )[branch],
            "signal": arrays_to_histograms(
                signal_arrays, masks["signal"], [variable]
            )[branch],
        }

        steps.append(
            IterativeCutStep(
                iteration=iteration,
                variable=branch,
                score=0.0,
                data_efficiency=data_eff,
                dvcs_mc_efficiency=signal_eff,
                pi0_mc_efficiency=0.0,
                f_pi0_before=0.0,
                f_pi0_after=0.0,
                f_pi0_raw_propagated=0.0,
                f_pi0_diagnostic=math.nan,
                boundaries=boundaries,
                before_histograms={"data": data_hist, "signal": signal_hist},
                after_histograms=after_histograms,
                fit_result=display_model,
                calibration_source=(
                    "raw-MC fallback"
                    if frozen.message == "raw-template fallback"
                    else "fitted nuisance"
                ),
            )
        )

        for name, containment in containments.items():
            mode = (
                "upper_quantile"
                if is_lower_bounded_morph_variable(variable)
                else "quantile_window"
            )
            json_variants[name]["data"][branch] = make_cut_stats(
                transformed,
                variable,
                boundaries[name]["data"][0],
                boundaries[name]["data"][1],
                containment,
                mode,
                iteration,
                0.0,
            )
            json_variants[name]["mc"][branch] = make_cut_stats(
                raw_shape,
                variable,
                boundaries[name]["mc"][0],
                boundaries[name]["mc"][1],
                containment,
                mode,
                iteration,
                0.0,
            )
        # endfor
    # endfor

    return steps, json_variants

def develop_iterative_cuts_for_period(
    period: PeriodConfig,
    topologies: Sequence[TopologyConfig],
    output_dir: Path,
    step_size: str,
    max_shift_bins: float,
    max_smear_bins: float,
    min_counts: int,
    fraction_variables: Sequence[str],
    shift_prior_bins: float,
    smear_prior_bins: float,
    use_nuisance_penalties: bool,
    dvcs_core_containment: float,
    dvcs_fraction_containment: float,
    pi0_core_containment: float,
    pi0_calibrations: Mapping[str, Mapping[str, Tuple[float, float]]],
    nominal_containment: float,
    loose_containment: float,
    tight_containment: float,
    outside_overshoot_penalty_weight: float,
    emiss2_mean_order_penalty_weight: float,
    dpi: int,
    preloaded_dvcs_arrays: Optional[
        Mapping[str, Mapping[str, Mapping[str, np.ndarray]]]
    ] = None,
) -> Tuple[Dict[str, Dict[str, object]], List[Dict[str, object]]]:
    containments = {
        "nominal": nominal_containment,
        "loose": loose_containment,
        "tight": tight_containment,
    }
    log(f"Loading event arrays for iterative cuts: {period.label}")
    dvcs_arrays = (
        dict(preloaded_dvcs_arrays)
        if preloaded_dvcs_arrays is not None
        else {
            "data": load_selected_arrays_for_file(
                period.data_file, topologies, step_size, VARIABLES
            ),
            "signal": load_selected_arrays_for_file(
                period.dvcs_mc_file, topologies, step_size, VARIABLES
            ),
            "background": load_selected_arrays_for_file(
                period.pi0_as_dvcs_mc_file, topologies, step_size, VARIABLES
            ),
        }
    )
    pi0_arrays = {
        "data": load_selected_arrays_for_file(period.eppi0_data_file, topologies, step_size, PI0_VARIABLES),
        "signal": load_selected_arrays_for_file(period.eppi0_mc_file, topologies, step_size, PI0_VARIABLES),
    }
    blocks: Dict[str, Dict[str, object]] = {name: {} for name in containments}
    rows: List[Dict[str, object]] = []
    for topology in topologies:
        dvcs_topology_arrays = {
            sample: values[topology.key]
            for sample, values in dvcs_arrays.items()
        }
        dvcs_steps, dvcs_json = run_dvcs_iterative_cuts(
            period,
            topology,
            dvcs_topology_arrays,
            max_shift_bins,
            max_smear_bins,
            min_counts,
            fraction_variables,
            shift_prior_bins,
            smear_prior_bins,
            use_nuisance_penalties,
            dvcs_core_containment,
            dvcs_fraction_containment,
            pi0_calibrations.get(topology.key, {}),
            containments,
            outside_overshoot_penalty_weight,
            emiss2_mean_order_penalty_weight,
        )
        dvcs_order = [step.variable for step in dvcs_steps]
        pi0_topology_arrays = {
            sample: values[topology.key]
            for sample, values in pi0_arrays.items()
        }
        pi0_steps, pi0_json = run_pi0_iterative_cuts(
            period,
            topology,
            pi0_topology_arrays,
            dvcs_order,
            max_shift_bins,
            max_smear_bins,
            min_counts,
            shift_prior_bins,
            smear_prior_bins,
            use_nuisance_penalties,
            pi0_core_containment,
            containments,
        )
        for channel, steps in (("dvcs", dvcs_steps), ("pi0", pi0_steps)):
            development_path = (
                output_dir / ("dvcs_channel" if channel == "dvcs" else "pi0_channel")
                / "iterative_cut_development"
                / f"iterative_cut_development_{period.key}_{topology.key.lower()}.png"
            )
            summary_path = (
                output_dir / ("dvcs_channel" if channel == "dvcs" else "pi0_channel")
                / "iterative_cut_summary"
                / f"iterative_cut_summary_{period.key}_{topology.key.lower()}.png"
            )
            draw_iterative_cut_canvas(
                development_path, period, topology, channel, steps, False, dpi
            )
            draw_iterative_cut_canvas(
                summary_path, period, topology, channel, steps, True, dpi
            )
            log(f"Wrote {development_path}")
            log(f"Wrote {summary_path}")
            for step in steps:
                rows.append({
                    "channel": channel,
                    "period": period.key,
                    "period_label": period.label,
                    "topology": topology.key,
                    "iteration": step.iteration + 1,
                    "variable": step.variable,
                    "score": step.score,
                    "f_pi0_before": step.f_pi0_before,
                    "f_pi0_after": step.f_pi0_after,
                    "f_pi0_raw_propagated": step.f_pi0_raw_propagated,
                    "f_pi0_diagnostic": step.f_pi0_diagnostic,
                    "calibration_source": step.calibration_source,
                    "data_efficiency": step.data_efficiency,
                    "signal_mc_efficiency": step.dvcs_mc_efficiency,
                    "pi0_mc_efficiency": step.pi0_mc_efficiency,
                    "nominal_data_low": step.boundaries["nominal"]["data"][0],
                    "nominal_data_high": step.boundaries["nominal"]["data"][1],
                    "nominal_mc_low": step.boundaries["nominal"]["mc"][0],
                    "nominal_mc_high": step.boundaries["nominal"]["mc"][1],
                    "loose_data_low": step.boundaries["loose"]["data"][0],
                    "loose_data_high": step.boundaries["loose"]["data"][1],
                    "loose_mc_low": step.boundaries["loose"]["mc"][0],
                    "loose_mc_high": step.boundaries["loose"]["mc"][1],
                    "tight_data_low": step.boundaries["tight"]["data"][0],
                    "tight_data_high": step.boundaries["tight"]["data"][1],
                    "tight_mc_low": step.boundaries["tight"]["mc"][0],
                    "tight_mc_high": step.boundaries["tight"]["mc"][1],
                })
            # endfor
        # endfor
        for name in containments:
            blocks[name][cut_json_key("dvcs", period, topology)] = dvcs_json[name]
            blocks[name][cut_json_key("pi0", period, topology)] = pi0_json[name]
        # endfor
    # endfor
    return blocks, rows


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    # endwith


def process_period_worker(
    period: PeriodConfig,
    topologies: Sequence[TopologyConfig],
    output_dir_string: str,
    step_size: str,
    max_shift_bins: float,
    max_smear_bins: float,
    fit_min_counts: int,
    fraction_variables: Sequence[str],
    shift_prior_bins: float,
    smear_prior_bins: float,
    use_nuisance_penalties: bool,
    pi0_core_containment: float,
    dvcs_core_containment: float,
    dvcs_fraction_containment: float,
    log_y: bool,
    dpi: int,
    cut_nominal_containment: float,
    cut_loose_containment: float,
    cut_tight_containment: float,
    outside_overshoot_penalty_weight: float,
    emiss2_mean_order_penalty_weight: float,
    run_iterative_cuts: bool,
) -> Tuple[
    str,
    List[Dict[str, object]],
    List[Dict[str, object]],
    Dict[str, Dict[str, object]],
    List[Dict[str, object]],
]:
    """Process one complete run period inside a worker process."""

    output_dir = Path(output_dir_string)
    log(f"[worker {os.getpid()}] Starting {period.label}")

    pi0_rows, pi0_calibrations = process_pi0_period(
        period,
        topologies,
        output_dir,
        step_size,
        max_shift_bins,
        max_smear_bins,
        fit_min_counts,
        shift_prior_bins,
        smear_prior_bins,
        use_nuisance_penalties,
        pi0_core_containment,
        log_y,
        dpi,
    )

    fraction_population_arrays = {
        "data": load_selected_arrays_for_file(
            period.data_file, topologies, step_size, VARIABLES
        ),
        "signal": load_selected_arrays_for_file(
            period.dvcs_mc_file, topologies, step_size, VARIABLES
        ),
        "background": load_selected_arrays_for_file(
            period.pi0_as_dvcs_mc_file, topologies, step_size, VARIABLES
        ),
    }

    dvcs_rows = process_period(
        period,
        topologies,
        output_dir,
        step_size,
        max_shift_bins,
        max_smear_bins,
        fit_min_counts,
        fraction_variables,
        shift_prior_bins,
        smear_prior_bins,
        use_nuisance_penalties,
        dvcs_core_containment,
        dvcs_fraction_containment,
        outside_overshoot_penalty_weight,
        emiss2_mean_order_penalty_weight,
        log_y,
        dpi,
        pi0_calibrations,
        fraction_population_arrays,
    )

    iterative_blocks: Dict[str, Dict[str, object]] = {
        "nominal": {}, "loose": {}, "tight": {}
    }
    iterative_rows: List[Dict[str, object]] = []
    if run_iterative_cuts:
        iterative_blocks, iterative_rows = develop_iterative_cuts_for_period(
            period,
            topologies,
            output_dir,
            step_size,
            max_shift_bins,
            max_smear_bins,
            fit_min_counts,
            fraction_variables,
            shift_prior_bins,
            smear_prior_bins,
            use_nuisance_penalties,
            dvcs_core_containment,
            dvcs_fraction_containment,
            pi0_core_containment,
            pi0_calibrations,
            cut_nominal_containment,
            cut_loose_containment,
            cut_tight_containment,
            outside_overshoot_penalty_weight,
            emiss2_mean_order_penalty_weight,
            dpi,
            fraction_population_arrays,
        )
    # endif

    log(f"[worker {os.getpid()}] Finished {period.label}")
    return period.key, dvcs_rows, pi0_rows, iterative_blocks, iterative_rows


def main() -> int:
    args = parse_args()
    if gaussian_filter1d is None or minimize is None or minimize_scalar is None or ndtr is None:
        raise RuntimeError("scipy is required for the shared two-template nuisance fit.")
    # endif
    if args.max_shift_bins <= 0.0 or args.max_smear_bins <= 0.0:
        raise ValueError("--max-shift-bins and --max-smear-bins must be positive.")
    # endif
    if args.workers <= 0:
        raise ValueError("--workers must be positive.")
    # endif
    if args.outside_overshoot_penalty < 0.0:
        raise ValueError("--outside-overshoot-penalty must be nonnegative.")
    # endif
    if args.emiss2_mean_order_penalty < 0.0:
        raise ValueError("--emiss2-mean-order-penalty must be nonnegative.")
    # endif
    if not (
        0.50 <= args.cut_tight_containment
        < args.cut_nominal_containment
        < args.cut_loose_containment
        <= 1.0
    ):
        raise ValueError(
            "Iterative containments must satisfy 0.50 <= tight < nominal < loose <= 1.0."
        )
    # endif
    if not (0.50 <= args.dvcs_core_containment <= 1.0):
        raise ValueError("--dvcs-core-containment must satisfy 0.50 <= value <= 1.0.")
    # endif
    if not (args.dvcs_core_containment <= args.dvcs_fraction_containment <= 1.0):
        raise ValueError(
            "--dvcs-fraction-containment must be at least the DVCS core "
            "containment and no greater than 1.0."
        )
    # endif
    if not (0.50 <= args.pi0_core_containment <= 1.0):
        raise ValueError("--pi0-core-containment must satisfy 0.50 <= value <= 1.0.")
    # endif

    periods = selected_periods(args.period)
    topologies = selected_topologies(args.topology)
    output_dir = Path(args.output_dir)
    all_rows: List[Dict[str, object]] = []
    all_pi0_rows: List[Dict[str, object]] = []
    fraction_variables = args.fraction_variable or ["Delta_phi", "theta_gamma_gamma", "pTmiss"]
    if args.shift_prior_bins <= 0.0 or args.smear_prior_bins <= 0.0:
        raise ValueError("--shift-prior-bins and --smear-prior-bins must be positive.")
    # endif

    if not args.no_clean_output:
        import shutil

        stale_paths = (
            output_dir / "dvcs_channel" / "iterative_cut_development",
            output_dir / "dvcs_channel" / "iterative_cut_summary",
            output_dir / "pi0_channel" / "iterative_cut_development",
            output_dir / "pi0_channel" / "iterative_cut_summary",
            output_dir / "iterative_cuts",
        )
        for stale_path in stale_paths:
            if stale_path.exists():
                shutil.rmtree(stale_path)
            # endif
        # endfor
        log("Removed stale iterative-cut outputs")
    # endif

    log(f"ROOT I/O backend: {io_backend()}")
    log(f"Producing DVCS shape and core-fit canvases plus direct-pi0 shape and core-fit canvases")
    log(
        "Selection: topology, (-t1) < 1.0, open_angle_ep2 > 5 deg and "
        "distinct FD sectors; no hard exclusivity cuts are applied."
    )
    log(
        f"Fit: shared f_pi0 uses the common in-range population of "
        f"{fraction_variables}. "
        f"Fraction-region containment={100.0 * args.dvcs_fraction_containment:.0f}%; "
        f"nuisance-region containment={100.0 * args.dvcs_core_containment:.0f}%. "
        "An asymmetric additive core morph is used for theta; ordinary additive core morphs are used for symmetric variables, lower-edge log morphs for pTmiss/theta_gamma_gamma/theta_pi0_pi0 and a bounded logit morph for z2 and an upper-edge log morph for xF2. "
        "DVCS shift priors are centered on the corresponding direct-pi0 core calibrations. "
        "Core-region fitting is the default; Gaussian nuisance penalties are "
        + ("enabled" if not args.disable_nuisance_penalties else "disabled")
        + f"; outside overshoot penalty={args.outside_overshoot_penalty:g}"
        + f"; Emiss2 mean-order penalty={args.emiss2_mean_order_penalty:g}"
    )

    worker_count = min(args.workers, 5, len(periods))
    log(
        f"Processing {len(periods)} run period(s) with "
        f"{worker_count} worker process(es); hard cap = 5"
    )

    period_results: Dict[
        str,
        Tuple[
            List[Dict[str, object]],
            List[Dict[str, object]],
            Dict[str, Dict[str, object]],
            List[Dict[str, object]],
        ],
    ] = {}

    if worker_count == 1:
        for period in periods:
            period_key, dvcs_rows, pi0_rows, cut_blocks, cut_rows = process_period_worker(
                period,
                topologies,
                str(output_dir),
                args.step_size,
                args.max_shift_bins,
                args.max_smear_bins,
                args.fit_min_counts,
                fraction_variables,
                args.shift_prior_bins,
                args.smear_prior_bins,
                not args.disable_nuisance_penalties,
                args.pi0_core_containment,
                args.dvcs_core_containment,
                args.dvcs_fraction_containment,
                args.log_y,
                args.dpi,
                args.cut_nominal_containment,
                args.cut_loose_containment,
                args.cut_tight_containment,
                args.outside_overshoot_penalty,
                args.emiss2_mean_order_penalty,
                not args.skip_iterative_cuts,
            )
            period_results[period_key] = (dvcs_rows, pi0_rows, cut_blocks, cut_rows)
        # endfor
    else:
        spawn_context = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=spawn_context,
        ) as executor:
            future_to_period = {
                executor.submit(
                    process_period_worker,
                    period,
                    topologies,
                    str(output_dir),
                    args.step_size,
                    args.max_shift_bins,
                    args.max_smear_bins,
                    args.fit_min_counts,
                    fraction_variables,
                    args.shift_prior_bins,
                    args.smear_prior_bins,
                    not args.disable_nuisance_penalties,
                    args.pi0_core_containment,
                    args.dvcs_core_containment,
                    args.dvcs_fraction_containment,
                    args.log_y,
                    args.dpi,
                    args.cut_nominal_containment,
                    args.cut_loose_containment,
                    args.cut_tight_containment,
                    args.outside_overshoot_penalty,
                    args.emiss2_mean_order_penalty,
                    not args.skip_iterative_cuts,
                ): period
                for period in periods
            }

            for future in concurrent.futures.as_completed(future_to_period):
                period = future_to_period[future]
                try:
                    period_key, dvcs_rows, pi0_rows, cut_blocks, cut_rows = future.result()
                except Exception as exc:
                    raise RuntimeError(
                        f"Worker failed for {period.label}: {exc}"
                    ) from exc
                # endtry

                period_results[period_key] = (dvcs_rows, pi0_rows, cut_blocks, cut_rows)
                log(f"Collected results for {period.label}")
            # endfor
        # endwith
    # endif

    combined_cut_blocks: Dict[str, Dict[str, object]] = {
        "nominal": {}, "loose": {}, "tight": {}
    }
    iterative_rows_all: List[Dict[str, object]] = []

    # Preserve the configured run-period order in the combined outputs.
    for period in periods:
        dvcs_rows, pi0_rows, cut_blocks, cut_rows = period_results[period.key]
        all_rows.extend(dvcs_rows)
        all_pi0_rows.extend(pi0_rows)
        iterative_rows_all.extend(cut_rows)
        for name in combined_cut_blocks:
            combined_cut_blocks[name].update(cut_blocks.get(name, {}))
        # endfor
    # endfor

    csv_path = output_dir / "dvcs_channel" / "template_fits" / "fit_results.csv"
    write_results_csv(csv_path, all_rows)
    log(f"Wrote {csv_path}")
    pi0_csv_path = output_dir / "pi0_channel" / "template_fits" / "fit_results.csv"
    write_results_csv(pi0_csv_path, all_pi0_rows)
    log(f"Wrote {pi0_csv_path}")
    if not args.skip_iterative_cuts:
        iterative_csv = output_dir / "iterative_cuts" / "iterative_cut_results.csv"
        write_results_csv(iterative_csv, iterative_rows_all)
        log(f"Wrote {iterative_csv}")

        json_dir = output_dir / "iterative_cuts" / "jsons"
        nominal_path = json_dir / "combined_cuts_95.json"
        loose_path = json_dir / "combined_cuts_99.json"
        tight_path = json_dir / "combined_cuts_90.json"
        compatibility_path = json_dir / "combined_cuts.json"
        write_json(nominal_path, combined_cut_blocks["nominal"])
        write_json(loose_path, combined_cut_blocks["loose"])
        write_json(tight_path, combined_cut_blocks["tight"])
        write_json(compatibility_path, combined_cut_blocks["nominal"])
        log(f"Wrote {nominal_path}")
        log(f"Wrote {loose_path}")
        log(f"Wrote {tight_path}")
        log(f"Wrote {compatibility_path}")
    # endif
    log("All requested shape plots and fits completed")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    # endtry
#endif
