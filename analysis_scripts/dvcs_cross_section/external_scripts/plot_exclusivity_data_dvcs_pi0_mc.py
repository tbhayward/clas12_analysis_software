#!/usr/bin/env python3
"""
plot_exclusivity_data_dvcs_pi0_mc.py

For each of the five run periods and three detector topologies, compare the
eight exclusivity/kinematic-variable shapes and perform selected-discriminator two-template fits while keeping all eight
projections as diagnostics:

    data = (1-f_pi0) * shifted-and-smeared DVCS MC
         + f_pi0     * eppi0 MC reconstructed as DVCS.

The detector topology plus the minimal global preselection are required:

    (-t1) < 1.0
    open_angle_ep2 > 5 degrees
    all reconstructed FD particles occupy distinct FD sectors

No hard exclusivity cuts are applied. The DVCS nuisance parameters are fitted
inside an MC-defined signal core, while the shared pi0 fraction is determined in
a broader MC-defined discriminator region. Full distributions remain visible as
validation. The fit uses raw binned data counts and a Poisson likelihood.
The total expected count is fixed to the observed number of entries inside
the plotted range, so the floated parameters describe shape only:

    f_pi0        : one fraction shared by the selected discriminator histograms
    shift        : variable-specific DVCS-template shift or log-scale shift
    sigma_add    : variable-specific additive or log-space smearing

The eppi0 template is fixed. Signed variables use additive shift plus Gaussian
smearing. Positive-definite pTmiss and theta_gamma_gamma use log-space shift
and smearing. All non-driver variables are validation projections: their DVCS nuisance
parameters are profiled at the fitted fraction, but they do not determine
f_pi0. Optional Gaussian nuisance penalties discourage extreme template
shifts and broadenings.

The script also compares reconstructed eppi0 data directly with reconstructed eppi0 MC using theta_pi0_pi0 in place of theta_gamma_gamma. The nuisance fits use xF, xF2, z2 and theta in place of the missing-mass projections. The shape-comparison canvases retain those new variables and additionally show Mx2, Mx2_1, Mx2_2 and pT in a third row. The direct eppi0 nuisance fit is restricted to an MC-defined signal core so that out-of-core backgrounds and tails remain diagnostics rather than forcing excessive template morphing.

Outputs:

  * one unit-area 3x4 shape-comparison canvas for each period/topology combination
  * one 4x2 DVCS-core two-template-fit canvas for each period/topology combination
  * fit_results.csv containing all fitted parameters and diagnostics

Dependencies: Python 3, numpy, matplotlib, scipy and either uproot or PyROOT.
"""

import argparse
import csv
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
    VariableConfig("theta_gamma_gamma", r"$\theta_{\gamma\gamma}$ (rad)", 100, 0.0, 2.0),
    VariableConfig("pTmiss", r"$p_{T}^{\mathrm{miss}}$ (GeV)", 100, 0.0, 0.3),
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
        r"$\theta_2(\gamma)$ (rad)",
        100,
        0.0,
        0.8,
        aliases=("theta2", "theta_2"),
    ),
)


PI0_VARIABLES: Tuple[VariableConfig, ...] = (
    VariableConfig("Delta_phi", r"$\Delta\phi$ (rad)", 100, 2.84159, 3.44159),
    VariableConfig("theta_pi0_pi0", r"$\theta_{\pi^0\pi^0}$ (rad)", 100, 0.0, 2.0),
    VariableConfig("pTmiss", r"$p_{T}^{\mathrm{miss}}$ (GeV)", 100, 0.0, 0.3),
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
        r"$\theta_2(\pi^0)$ (rad)",
        100,
        0.0,
        0.8,
        aliases=("theta2", "theta_2"),
    ),
)


# The nuisance/template fits continue to use the eight variables above.
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
            "Broader DVCS-MC containment used by the selected discriminator "
            "histograms to determine the shared pi0 fraction (default: 0.95)."
        ),
    )
    parser.add_argument(
        "--pi0-core-containment",
        type=float,
        default=0.90,
        help=(
            "MC signal containment used by the direct eppi0 data--MC core fits "
            "(default: 0.90). Signed variables use equal tails; positive-definite "
            "variables use the interval from the lower boundary to this quantile."
        ),
    )
    parser.add_argument(
        "--fraction-variable",
        choices=[v.branch for v in VARIABLES],
        action="append",
        help=(
            "Variable used to determine the shared pi0 fraction. May be supplied "
            "repeatedly. Default: pTmiss and theta_gamma_gamma."
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


def is_positive_morph_variable(variable: VariableConfig) -> bool:
    return variable.branch in {
        "pTmiss",
        "theta_gamma_gamma",
        "theta_pi0_pi0",
        "theta",
        "z2",
    }


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


def transform_dvcs_shape(
    base_shape: np.ndarray,
    variable: VariableConfig,
    shift: float,
    sigma_add: float,
) -> Optional[np.ndarray]:
    if is_positive_morph_variable(variable):
        return transform_positive_shape(base_shape, variable, shift, sigma_add)
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
        if is_positive_morph_variable(variable):
            center = calibrated_shift_center(variable)
            return [(max(-0.70, center - 0.40), min(0.70, center + 0.40)), (0.0, 1.00)]
        # endif
        bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
        center = calibrated_shift_center(variable)
        half_range = max_shift_bins * bin_width
        return [(center - half_range, center + half_range), (0.0, max_smear_bins * bin_width)]

    def nuisance_start(variable: VariableConfig) -> np.ndarray:
        center = calibrated_shift_center(variable)
        if is_positive_morph_variable(variable):
            return np.asarray([center, 0.10], dtype=np.float64)
        # endif
        bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
        return np.asarray([center, 2.0 * bin_width], dtype=np.float64)

    def nuisance_penalty(variable: VariableConfig, shift: float, sigma_add: float) -> float:
        if not use_nuisance_penalties:
            return 0.0
        # endif
        shift_center = calibrated_shift_center(variable)
        if is_positive_morph_variable(variable):
            shift_width = 0.20
            smear_width = 0.40
        else:
            bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
            shift_width = max(shift_prior_bins * bin_width, 1.0e-12)
            smear_width = max(smear_prior_bins * bin_width, 1.0e-12)
        # endif
        return 0.5 * ((shift - shift_center) / shift_width) ** 2 + 0.5 * (sigma_add / smear_width) ** 2

    def build_variable_model(
        variable: VariableConfig,
        fraction: float,
        shift: float,
        sigma_add: float,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        info = prepared[variable.branch]
        transformed = transform_dvcs_shape(info["dvcs_shape"], variable, shift, sigma_add)
        if transformed is None:
            return None
        # endif
        pi0_shape = info["pi0_shape"]
        total_shape = np.clip((1.0 - fraction) * transformed + fraction * pi0_shape, 0.0, None)
        normalization = float(np.sum(total_shape))
        if normalization <= 0.0:
            return None
        # endif
        data_total = float(info["data_total"])
        model = data_total * total_shape / normalization
        dvcs_component = data_total * (1.0 - fraction) * transformed / normalization
        pi0_component = data_total * fraction * pi0_shape / normalization
        return model, dvcs_component, pi0_component, transformed

    def objective_for_mask(
        variable: VariableConfig,
        fraction: float,
        nuisance: np.ndarray,
        mask_name: str,
        include_penalty: bool,
    ) -> float:
        shift = float(nuisance[0])
        sigma_add = float(nuisance[1])
        built = build_variable_model(variable, fraction, shift, sigma_add)
        if built is None:
            return 1.0e100
        # endif
        info = prepared[variable.branch]
        mask = np.asarray(info[mask_name], dtype=bool)
        value = 0.5 * poisson_deviance(info["data"][mask], built[0][mask])
        if include_penalty:
            value += nuisance_penalty(variable, shift, sigma_add)
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
        built = build_variable_model(variable, best_fraction, shift, sigma_add)
        if built is None:
            results[variable.branch] = FitResult(False, "invalid fitted template")
            continue
        # endif
        model, dvcs_component, pi0_component, transformed = built
        display_mask = np.asarray(
            info["fraction_mask"] if variable.branch in requested_set else info["core_mask"],
            dtype=bool,
        )
        variable_deviance = poisson_deviance(info["data"][display_mask], model[display_mask])
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
            deviance=variable_deviance,
            ndf=max(0, used_bins - 2),
            data_total=float(info["data_total"]),
            model_counts=model,
            dvcs_component_counts=dvcs_component,
            pi0_component_counts=pi0_component,
            transformed_dvcs_shape=transformed,
            fit_mask=display_mask,
            morph_label="log-space" if is_positive_morph_variable(variable) else "additive",
            excluded_data_counts=excluded_data,
            excluded_model_counts=excluded_model,
            excluded_excess_counts=excluded_excess,
            excluded_excess_fraction=excluded_fraction,
        )
    # endfor

    n_driver_parameters = 1 + 2 * len(fraction_variables)
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
                rf"$\sigma={result.sigma_add:.4g}$ ({result.morph_label})" + "\n"
                f"$D/ndf={result.deviance:.1f}/{result.ndf}$"
            )
            if result.excluded_data_counts > 0.0:
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
        f"fraction drivers: {', '.join(shared_summary.fraction_variables)}; "
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
    Positive-definite variables retain the interval from the lower plot edge
    through the requested cumulative containment. Existing variable-specific
    edge masks are also respected.
    """
    counts = np.asarray(mc_counts, dtype=np.float64)
    total = float(np.sum(counts))
    base_mask = fit_mask_for_variable(variable, topology)
    if total <= 0.0 or not math.isfinite(total):
        return base_mask
    # endif

    cumulative = np.cumsum(counts) / total
    if is_positive_morph_variable(variable):
        lower_index = 0
        upper_index = int(np.searchsorted(cumulative, containment, side="left"))
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
    """Fit the MC-defined exclusive-pi0 signal core with shift and smearing."""

    data = np.asarray(data_counts, dtype=np.float64)
    mc_shape = normalized_shape(mc_counts)
    data_total = float(np.sum(data))
    if data_total < min_counts or mc_shape is None:
        return FitResult(False, "insufficient counts", data_total=data_total)
    # endif

    mask = mc_signal_containment_mask(mc_counts, variable, topology, core_containment)
    core_data_total = float(np.sum(data[mask]))
    if core_data_total < min_counts:
        return FitResult(False, "insufficient counts in MC-defined signal core", data_total=data_total, fit_mask=mask)
    # endif
    bin_width = (variable.xmax - variable.xmin) / variable.bins
    positive = is_positive_morph_variable(variable)
    if positive:
        bounds = [(-0.8, 0.8), (0.0, 1.2)]
        x0 = np.asarray([0.0, 0.10], dtype=np.float64)
    else:
        bounds = [(-max_shift_bins * bin_width, max_shift_bins * bin_width),
                  (0.0, max_smear_bins * bin_width)]
        x0 = np.asarray([0.0, 0.5 * bin_width], dtype=np.float64)
    # endif

    def build(params: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        transformed = transform_dvcs_shape(mc_shape, variable, float(params[0]), float(params[1]))
        norm = float(np.sum(transformed[mask]))
        if norm <= 0.0:
            return np.zeros_like(data), transformed
        # endif
        # Normalize to the observed yield inside the fitted signal core. The
        # tails remain a genuine prediction rather than forcing the nuisance
        # parameters to account for every out-of-core data event.
        model = core_data_total * transformed / norm
        return model, transformed

    def objective(params: np.ndarray) -> float:
        model, _ = build(params)
        value = poisson_deviance(data[mask], np.clip(model[mask], 1.0e-12, None))
        if use_nuisance_penalties:
            if positive:
                value += (float(params[0]) / 0.20) ** 2 + (float(params[1]) / 0.40) ** 2
            else:
                value += (float(params[0]) / (shift_prior_bins * bin_width)) ** 2
                value += (float(params[1]) / (smear_prior_bins * bin_width)) ** 2
            # endif
        # endif
        return float(value)

    result = minimize(objective, x0, method="L-BFGS-B", bounds=bounds)
    model, transformed = build(np.asarray(result.x, dtype=np.float64))
    nfit = int(np.count_nonzero(mask))
    return FitResult(
        success=bool(result.success),
        message=f"MC-defined {100.0 * core_containment:.1f}% signal-core fit",
        shift=float(result.x[0]),
        sigma_add=float(result.x[1]),
        deviance=poisson_deviance(data[mask], np.clip(model[mask], 1.0e-12, None)),
        ndf=max(1, nfit - 2),
        data_total=data_total,
        model_counts=model,
        dvcs_component_counts=model,
        transformed_dvcs_shape=transformed,
        fit_mask=mask,
        morph_label="log-space" if positive else "additive",
        excluded_data_counts=float(np.sum(data[~mask])),
        excluded_model_counts=float(np.sum(model[~mask])),
        excluded_excess_counts=float(np.sum(data[~mask]) - np.sum(model[~mask])),
        excluded_excess_fraction=(
            float(np.sum(data[~mask]) - np.sum(model[~mask])) / data_total
            if data_total > 0.0 else 0.0
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
            rf"$\sigma_{{add}}={result.sigma_add:.4g}$ ({result.morph_label})" + "\n" +
            rf"$D/ndf={result.deviance:.1f}/{result.ndf}$" + "\n" +
            f"outside-core excess={result.excluded_excess_counts:.0f}",
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
            additive = not is_positive_morph_variable(variable)
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
                "core_containment": pi0_core_containment,
                "fit_bins_used": int(np.count_nonzero(result.fit_mask)) if result.fit_mask is not None else 0,
                "outside_core_data_counts": result.excluded_data_counts,
                "outside_core_model_counts": result.excluded_model_counts,
                "outside_core_excess_counts": result.excluded_excess_counts,
                "shift_or_log_shift": result.shift,
                "shift_bins": result.shift / bin_width if result.success and additive else math.nan,
                "sigma_or_log_sigma": result.sigma_add,
                "sigma_bins": result.sigma_add / bin_width if result.success and additive else math.nan,
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
    log_y: bool,
    dpi: int,
    pi0_core_calibrations: Optional[Mapping[str, Mapping[str, Tuple[float, float]]]] = None,
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

        shared_summary = fit_shared_two_templates(
            data_hists[topology.key],
            dvcs_hists[topology.key],
            pi0_hists[topology.key],
            topology,
            max_shift_bins, max_smear_bins, fit_min_counts,
            fraction_variables, shift_prior_bins, smear_prior_bins,
            use_nuisance_penalties,
            dvcs_core_containment,
            dvcs_fraction_containment,
            (pi0_core_calibrations or {}).get(topology.key, {}),
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
                data_hists[topology.key], dvcs_hists[topology.key],
                pi0_hists[topology.key], topology,
                max_shift_bins, max_smear_bins, fit_min_counts,
                [branch], shift_prior_bins, smear_prior_bins,
                use_nuisance_penalties,
                dvcs_core_containment,
                dvcs_fraction_containment,
                (pi0_core_calibrations or {}).get(topology.key, {}),
            )
        # endfor

        for variable in VARIABLES:
            result = fit_results[variable.branch]
            bin_width = (variable.xmax - variable.xmin) / variable.bins
            additive = not is_positive_morph_variable(variable)
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
                "global_poisson_deviance": shared_summary.deviance,
                "global_ndf": shared_summary.ndf,
                "data_entries_in_range": result.data_total,
                "dvcs_mc_entries_in_range": float(np.sum(dvcs_hists[topology.key][variable.branch])),
                "pi0_mc_entries_in_range": float(np.sum(pi0_hists[topology.key][variable.branch])),
                "morph_type": result.morph_label,
                "shift_or_log_shift": result.shift,
                "shift_err": result.shift_err,
                "shift_bins": result.shift / bin_width if result.success and additive else math.nan,
                "sigma_or_log_sigma": result.sigma_add,
                "sigma_err": result.sigma_add_err,
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


def main() -> int:
    args = parse_args()
    if gaussian_filter1d is None or minimize is None or minimize_scalar is None or ndtr is None:
        raise RuntimeError("scipy is required for the shared two-template nuisance fit.")
    # endif
    if args.max_shift_bins <= 0.0 or args.max_smear_bins <= 0.0:
        raise ValueError("--max-shift-bins and --max-smear-bins must be positive.")
    # endif
    if not (0.50 <= args.dvcs_core_containment < 1.0):
        raise ValueError("--dvcs-core-containment must satisfy 0.50 <= value < 1.0.")
    # endif
    if not (args.dvcs_core_containment <= args.dvcs_fraction_containment < 1.0):
        raise ValueError(
            "--dvcs-fraction-containment must be at least the DVCS core containment and below 1.0."
        )
    # endif
    if not (0.50 <= args.pi0_core_containment < 1.0):
        raise ValueError("--pi0-core-containment must satisfy 0.50 <= value < 1.0.")
    # endif

    periods = selected_periods(args.period)
    topologies = selected_topologies(args.topology)
    output_dir = Path(args.output_dir)
    all_rows: List[Dict[str, object]] = []
    all_pi0_rows: List[Dict[str, object]] = []
    fraction_variables = args.fraction_variable or ["pTmiss", "theta_gamma_gamma"]
    if args.shift_prior_bins <= 0.0 or args.smear_prior_bins <= 0.0:
        raise ValueError("--shift-prior-bins and --smear-prior-bins must be positive.")
    # endif

    log(f"ROOT I/O backend: {io_backend()}")
    log(f"Producing DVCS shape and core-fit canvases plus direct-pi0 shape and core-fit canvases")
    log(
        "Selection: topology, (-t1) < 1.0, open_angle_ep2 > 5 deg and "
        "distinct FD sectors; no hard exclusivity cuts are applied."
    )
    log(
        f"Fit: f_pi0 determined only by {fraction_variables} in the "
        f"{100.0 * args.dvcs_fraction_containment:.0f}% DVCS-MC regions; nuisance morphs use "
        f"the {100.0 * args.dvcs_core_containment:.0f}% DVCS cores. All other variables are validation projections. "
        "Additive morphs are used for signed variables and log-space morphs for pTmiss, theta_gamma_gamma and theta_pi0_pi0. "
        "DVCS shift priors are centered on the corresponding direct-pi0 core calibrations. "
        "Gaussian nuisance penalties are " + ("enabled" if not args.disable_nuisance_penalties else "disabled")
    )

    for period in periods:
        pi0_rows, pi0_calibrations = process_pi0_period(
            period, topologies, output_dir, args.step_size,
            args.max_shift_bins, args.max_smear_bins, args.fit_min_counts,
            args.shift_prior_bins, args.smear_prior_bins,
            not args.disable_nuisance_penalties, args.pi0_core_containment,
            args.log_y, args.dpi,
        )
        all_pi0_rows.extend(pi0_rows)

        all_rows.extend(process_period(
            period, topologies, output_dir, args.step_size,
            args.max_shift_bins, args.max_smear_bins, args.fit_min_counts,
            fraction_variables, args.shift_prior_bins, args.smear_prior_bins,
            not args.disable_nuisance_penalties,
            args.dvcs_core_containment, args.dvcs_fraction_containment,
            args.log_y, args.dpi, pi0_calibrations,
        ))
    # endfor

    csv_path = output_dir / "dvcs_channel" / "template_fits" / "fit_results.csv"
    write_results_csv(csv_path, all_rows)
    log(f"Wrote {csv_path}")
    pi0_csv_path = output_dir / "pi0_channel" / "template_fits" / "fit_results.csv"
    write_results_csv(pi0_csv_path, all_pi0_rows)
    log(f"Wrote {pi0_csv_path}")
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
