#!/usr/bin/env python3
"""
plot_exclusivity_data_dvcs_pi0_mc.py

For each of the five run periods and three detector topologies, compare the
eight exclusivity-variable shapes and perform one simultaneous two-template
fit across all eight projections:

    data = (1-f_pi0) * shifted-and-smeared DVCS MC
         + f_pi0     * eppi0 MC reconstructed as DVCS.

The detector topology plus the minimal global preselection are required:

    (-t1) < 1.0
    open_angle_ep2 > 5 degrees
    all reconstructed FD particles occupy distinct FD sectors

No exclusivity cuts are applied. The fit uses the raw binned data counts and a Poisson likelihood.
The total expected count is fixed to the observed number of entries inside
the plotted range, so the floated parameters describe shape only:

    f_pi0        : one fraction shared by all eight histograms
    shift        : variable-specific DVCS-template shift or log-scale shift
    sigma_add    : variable-specific additive or log-space smearing

The eppi0 template is fixed. Signed variables use additive shift plus Gaussian
smearing. Positive-definite pTmiss and theta_gamma_gamma use log-space shift
and smearing. The clearly unmodelled right side of Mx2_1 is displayed but
excluded from the fit objective, and xF edge bins are excluded to reduce
range-clipping bias.

Outputs:

  * one unit-area 4x2 shape-comparison canvas for each period/topology combination
  * one 4x2 two-template-fit canvas for each period/topology combination
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


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        key="fa18_inb",
        label="Fa18 Inb",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_50nA_10604MeV_epgamma.root",
    ),
    PeriodConfig(
        key="fa18_out",
        label="Fa18 Out",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_50nA_10604MeV_epgamma.root",
    ),
    PeriodConfig(
        key="sp19_inb",
        label="Sp19 Inb",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_50nA_10200MeV_epgamma.root",
    ),
    PeriodConfig(
        key="sp18_inb",
        label="Sp18 Inb",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_50nA_10594MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_inb_50nA_10594MeV_epgamma.root",
    ),
    PeriodConfig(
        key="sp18_out",
        label="Sp18 Out",
        data_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root",
        dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_45nA_10594MeV.root",
        pi0_as_dvcs_mc_file="/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_out_45nA_10594MeV_epgamma.root",
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
    VariableConfig("xF", r"$x_F$", 100, -0.5, 0.2),
    VariableConfig("Emiss2", r"$E_{\mathrm{miss}}^{2}$ (GeV$^2$)", 100, -1.0, 2.0),
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
)


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


def resolve_variable_branches(path: str) -> Dict[str, str]:
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

    for variable in VARIABLES:
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


def empty_histograms() -> Dict[str, np.ndarray]:
    return {
        variable.branch: np.zeros(variable.bins, dtype=np.float64)
        for variable in VARIABLES
    }


def fill_histograms_uproot(
    path: str,
    topologies: Sequence[TopologyConfig],
    step_size: str,
    apply_mx2_1_cut: bool = False,
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, int]]:
    """Read one ROOT file once and fill all requested topology/variable histograms."""

    require_input_file(path)
    resolved = resolve_variable_branches(path)

    expressions = [
        "detector1", "detector2", "t1", "open_angle_ep2",
        "e_phi", "p1_phi", "p2_phi",
    ] + sorted(set(resolved.values()))
    histograms = {topology.key: empty_histograms() for topology in topologies}
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

            for variable in VARIABLES:
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
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, int]]:
    del step_size  # PyROOT performs a direct TTree event loop.

    require_input_file(path)
    resolved = resolve_variable_branches(path)
    histograms = {topology.key: empty_histograms() for topology in topologies}
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

            for variable in VARIABLES:
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
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, int]]:
    backend = io_backend()

    if backend == "uproot":
        return fill_histograms_uproot(path, topologies, step_size, apply_mx2_1_cut)
    # endif

    return fill_histograms_pyroot(path, topologies, step_size, apply_mx2_1_cut)


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
    return variable.branch in {"pTmiss", "theta_gamma_gamma"}


def mx2_1_upper_cut(topology: TopologyConfig) -> float:
    """Topology-dependent hard upper cut used for the second-stage plots and fits."""
    return 0.30 if topology.detector2 == 1 else 0.40


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
    integrated Gaussian probability. This avoids the hard interpolation cutoff
    that previously made some morphed Mx2_2 curves stop near -0.5 GeV^2.
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
) -> SharedFitSummary:
    """Fit one shared f_pi0 and variable-specific DVCS morph nuisances."""
    prepared: Dict[str, Dict[str, object]] = {}
    active_variables: List[VariableConfig] = []

    for variable in VARIABLES:
        data = np.asarray(data_histograms[variable.branch], dtype=np.float64)
        dvcs_shape = normalized_shape(dvcs_histograms[variable.branch])
        pi0_shape = normalized_shape(pi0_histograms[variable.branch])
        data_total = float(np.sum(data))
        if data_total < min_counts or dvcs_shape is None or pi0_shape is None:
            continue
        # endif
        prepared[variable.branch] = {
            "data": data,
            "data_total": data_total,
            "dvcs_shape": dvcs_shape,
            "pi0_shape": pi0_shape,
            "mask": fit_mask_for_variable(variable, topology),
        }
        active_variables.append(variable)
    # endfor

    if not active_variables:
        return SharedFitSummary(False, "no variables have sufficient data and MC")
    # endif

    bounds: List[Tuple[float, float]] = [(0.0, 1.0)]
    starts_by_variable: List[Tuple[float, float]] = []
    for variable in active_variables:
        if is_positive_morph_variable(variable):
            # Multiplicative scale exp(shift), with Gaussian width in log space.
            bounds.extend([(-0.70, 0.70), (0.0, 1.00)])
            starts_by_variable.append((0.0, 0.10))
        else:
            bin_width = (variable.xmax - variable.xmin) / float(variable.bins)
            bounds.extend([
                (-max_shift_bins * bin_width, max_shift_bins * bin_width),
                (0.0, max_smear_bins * bin_width),
            ])
            starts_by_variable.append((0.0, 2.0 * bin_width))
        # endif
    # endfor

    def build_variable_model(
        variable: VariableConfig,
        fraction: float,
        shift: float,
        sigma_add: float,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        info = prepared[variable.branch]
        transformed = transform_dvcs_shape(
            info["dvcs_shape"], variable, shift, sigma_add
        )
        if transformed is None:
            return None
        # endif
        pi0_shape = info["pi0_shape"]
        total_shape = (1.0 - fraction) * transformed + fraction * pi0_shape
        total_shape = np.clip(total_shape, 0.0, None)
        normalization = float(np.sum(total_shape))
        if normalization <= 0.0:
            return None
        # endif
        data_total = float(info["data_total"])
        model = data_total * total_shape / normalization
        dvcs_component = data_total * (1.0 - fraction) * transformed / normalization
        pi0_component = data_total * fraction * pi0_shape / normalization
        return model, dvcs_component, pi0_component, transformed

    def objective(parameters: np.ndarray) -> float:
        fraction = float(parameters[0])
        total_nll = 0.0
        offset = 1
        for variable in active_variables:
            shift = float(parameters[offset])
            sigma_add = float(parameters[offset + 1])
            offset += 2
            built = build_variable_model(variable, fraction, shift, sigma_add)
            if built is None:
                return 1.0e100
            # endif
            model = built[0]
            info = prepared[variable.branch]
            mask = info["mask"]
            total_nll += 0.5 * poisson_deviance(info["data"][mask], model[mask])
        # endfor
        return total_nll

    def variable_objective(
        variable: VariableConfig,
        fraction: float,
        nuisance: np.ndarray,
    ) -> float:
        built = build_variable_model(
            variable, fraction, float(nuisance[0]), float(nuisance[1])
        )
        if built is None:
            return 1.0e100
        # endif
        info = prepared[variable.branch]
        mask = info["mask"]
        return 0.5 * poisson_deviance(info["data"][mask], built[0][mask])

    best_value = math.inf
    best_parameters: Optional[np.ndarray] = None
    best_message = "coordinate-profile fit did not run"

    # The likelihood is separable by variable once f_pi0 is fixed. Alternating
    # variable-wise nuisance fits with a one-dimensional shared-fraction update
    # is much faster and more stable than a 17-parameter numerical-gradient fit.
    for initial_fraction in (0.15, 0.40, 0.70):
        fraction = initial_fraction
        nuisances = np.asarray(starts_by_variable, dtype=np.float64)
        previous_value = math.inf

        for iteration in range(10):
            for index, variable in enumerate(active_variables):
                nuisance_result = minimize(
                    lambda values, v=variable, f=fraction: variable_objective(v, f, values),
                    nuisances[index],
                    method="L-BFGS-B",
                    bounds=bounds[1 + 2 * index : 1 + 2 * index + 2],
                    options={"maxiter": 400, "ftol": 1.0e-9},
                )
                if nuisance_result.success and np.all(np.isfinite(nuisance_result.x)):
                    nuisances[index] = nuisance_result.x
                # endif
            # endfor

            def fraction_objective(candidate_fraction: float) -> float:
                total = 0.0
                for index, variable in enumerate(active_variables):
                    total += variable_objective(variable, candidate_fraction, nuisances[index])
                # endfor
                return total

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

        parameters = [fraction]
        for nuisance in nuisances:
            parameters.extend([float(nuisance[0]), float(nuisance[1])])
        # endfor
        parameter_array = np.asarray(parameters, dtype=np.float64)
        value = objective(parameter_array)
        if value < best_value:
            best_value = value
            best_parameters = parameter_array
            best_message = "coordinate-profile fit converged"
        # endif
    # endfor

    if best_parameters is None or not np.all(np.isfinite(best_parameters)):
        return SharedFitSummary(False, best_message)
    # endif

    class BestFit:
        pass

    best = BestFit()
    best.x = best_parameters
    best.fun = best_value
    best.message = best_message

    errors = np.full(len(best.x), np.nan)
    try:
        fraction_step = 1.0e-3
        f0 = float(best.x[0])
        if fraction_step < f0 < 1.0 - fraction_step:
            left = best.x.copy()
            right = best.x.copy()
            left[0] -= fraction_step
            right[0] += fraction_step
            curvature = (objective(left) - 2.0 * objective(best.x) + objective(right)) / (fraction_step ** 2)
            if curvature > 0.0:
                errors[0] = 1.0 / math.sqrt(curvature)
            # endif
        # endif
    except Exception:
        pass
    # endtry

    fraction = float(best.x[0])
    results: Dict[str, FitResult] = {}
    total_deviance = 0.0
    total_used_bins = 0
    offset = 1
    for variable in VARIABLES:
        if variable.branch not in prepared:
            results[variable.branch] = FitResult(False, "insufficient data or empty template")
            continue
        # endif

        shift = float(best.x[offset])
        sigma_add = float(best.x[offset + 1])
        shift_err = float(errors[offset])
        sigma_err = float(errors[offset + 1])
        offset += 2
        built = build_variable_model(variable, fraction, shift, sigma_add)
        info = prepared[variable.branch]
        if built is None:
            results[variable.branch] = FitResult(False, "invalid fitted template")
            continue
        # endif
        model, dvcs_component, pi0_component, transformed = built
        mask = np.asarray(info["mask"], dtype=bool)
        variable_deviance = poisson_deviance(info["data"][mask], model[mask])
        used_bins = int(np.count_nonzero(mask))
        total_deviance += variable_deviance
        total_used_bins += used_bins
        excluded = ~mask
        excluded_data = float(np.sum(info["data"][excluded]))
        excluded_model = float(np.sum(model[excluded]))
        excluded_excess = excluded_data - excluded_model
        excluded_fraction = (
            excluded_excess / float(info["data_total"])
            if float(info["data_total"]) > 0.0 else 0.0
        )
        results[variable.branch] = FitResult(
            success=True,
            message=str(best.message),
            f_pi0=fraction,
            f_pi0_err=float(errors[0]),
            shift=shift,
            shift_err=shift_err,
            sigma_add=sigma_add,
            sigma_add_err=sigma_err,
            deviance=variable_deviance,
            ndf=max(0, used_bins - 2),
            data_total=float(info["data_total"]),
            model_counts=model,
            dvcs_component_counts=dvcs_component,
            pi0_component_counts=pi0_component,
            transformed_dvcs_shape=transformed,
            fit_mask=mask,
            morph_label="log-space" if is_positive_morph_variable(variable) else "additive",
            excluded_data_counts=excluded_data,
            excluded_model_counts=excluded_model,
            excluded_excess_counts=excluded_excess,
            excluded_excess_fraction=excluded_fraction,
        )
    # endfor

    n_parameters = 1 + 2 * len(active_variables)
    return SharedFitSummary(
        success=True,
        message=str(best.message),
        f_pi0=fraction,
        f_pi0_err=float(errors[0]),
        deviance=total_deviance,
        ndf=max(0, total_used_bins - n_parameters),
        variable_results=results,
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
    fig, axes = plt.subplots(2, 4, figsize=(18.0, 8.8))
    flat_axes = axes.ravel()

    for axis, variable in zip(flat_axes, VARIABLES):
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
    fig.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955),
        ncol=len(labels), frameon=False,
    )
    fig.suptitle(
        f"Exclusivity shapes after {selection_label}: {period.label}, {topology.label}\n"
        f"selected entries: data={selected_counts['data']:,}, "
        f"DVCS MC={selected_counts['dvcs_mc']:,}, "
        rf"$e\pi^0$ MC as DVCS={selected_counts['pi0_mc']:,}",
        fontsize=15, y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))
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
    fig, axes = plt.subplots(2, 4, figsize=(18.0, 8.8))
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
                rf"shared $f_{{\pi^0}}={result.f_pi0:.3f}$" + "\n"
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
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=len(labels), frameon=False)
    fig.suptitle(
        f"Two-template fits after $M_{{x1}}^2<{mx2_1_upper_cut(topology):.1f}$ GeV$^2$: {period.label}, {topology.label}\n"
        f"topology-selected entries: data={selected_counts['data']:,}, "
        f"DVCS MC={selected_counts['dvcs_mc']:,}, "
        rf"$e\pi^0$ MC as DVCS={selected_counts['pi0_mc']:,}",
        fontsize=15, y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)



def process_period(
    period: PeriodConfig,
    topologies: Sequence[TopologyConfig],
    output_dir: Path,
    step_size: str,
    max_shift_bins: float,
    max_smear_bins: float,
    fit_min_counts: int,
    log_y: bool,
    dpi: int,
) -> List[Dict[str, object]]:
    log(f"Starting period {period.label}")

    # Stage 1: minimal global preselection only.
    data_uncut, data_uncut_counts = fill_histograms_for_file(
        period.data_file, topologies, step_size, apply_mx2_1_cut=False
    )
    dvcs_uncut, dvcs_uncut_counts = fill_histograms_for_file(
        period.dvcs_mc_file, topologies, step_size, apply_mx2_1_cut=False
    )
    pi0_uncut, pi0_uncut_counts = fill_histograms_for_file(
        period.pi0_as_dvcs_mc_file, topologies, step_size, apply_mx2_1_cut=False
    )

    # Stage 2 and fit input: apply the topology-dependent Mx2_1 upper cut.
    data_cut, data_cut_counts = fill_histograms_for_file(
        period.data_file, topologies, step_size, apply_mx2_1_cut=True
    )
    dvcs_cut, dvcs_cut_counts = fill_histograms_for_file(
        period.dvcs_mc_file, topologies, step_size, apply_mx2_1_cut=True
    )
    pi0_cut, pi0_cut_counts = fill_histograms_for_file(
        period.pi0_as_dvcs_mc_file, topologies, step_size, apply_mx2_1_cut=True
    )

    rows: List[Dict[str, object]] = []

    for topology in topologies:
        uncut_selected = {
            "data": data_uncut_counts[topology.key],
            "dvcs_mc": dvcs_uncut_counts[topology.key],
            "pi0_mc": pi0_uncut_counts[topology.key],
        }
        uncut_path = output_dir / "shape_comparisons" / f"exclusivity_shapes_{period.key}_{topology.key.lower()}.png"
        draw_shape_canvas(
            uncut_path, period, topology, data_uncut[topology.key],
            dvcs_uncut[topology.key], pi0_uncut[topology.key],
            uncut_selected, log_y, dpi, "minimal preselection",
        )
        log(f"Wrote {uncut_path}")

        cut_value = mx2_1_upper_cut(topology)
        cut_selected = {
            "data": data_cut_counts[topology.key],
            "dvcs_mc": dvcs_cut_counts[topology.key],
            "pi0_mc": pi0_cut_counts[topology.key],
        }
        cut_shape_path = output_dir / "shape_comparisons_mx2_1_cut" / f"exclusivity_shapes_mx2_1_cut_{period.key}_{topology.key.lower()}.png"
        draw_shape_canvas(
            cut_shape_path, period, topology, data_cut[topology.key],
            dvcs_cut[topology.key], pi0_cut[topology.key],
            cut_selected, log_y, dpi,
            rf"minimal preselection and $M_{{x1}}^2<{cut_value:.1f}$ GeV$^2$",
        )
        log(f"Wrote {cut_shape_path}")

        shared_summary = fit_shared_two_templates(
            data_cut[topology.key],
            dvcs_cut[topology.key],
            pi0_cut[topology.key],
            topology,
            max_shift_bins, max_smear_bins, fit_min_counts,
        )
        if not shared_summary.success or shared_summary.variable_results is None:
            raise RuntimeError(
                f"Shared fit failed for {period.label} {topology.label}: {shared_summary.message}"
            )
        # endif
        fit_results = shared_summary.variable_results

        for variable in VARIABLES:
            result = fit_results[variable.branch]
            bin_width = (variable.xmax - variable.xmin) / variable.bins
            additive = not is_positive_morph_variable(variable)
            rows.append({
                "period": period.key, "period_label": period.label,
                "topology": topology.key, "variable": variable.branch,
                "mx2_1_upper_cut": cut_value,
                "success": int(result.success), "message": result.message,
                "shared_f_pi0": shared_summary.f_pi0,
                "shared_f_pi0_err": shared_summary.f_pi0_err,
                "global_poisson_deviance": shared_summary.deviance,
                "global_ndf": shared_summary.ndf,
                "data_entries_in_range": result.data_total,
                "dvcs_mc_entries_in_range": float(np.sum(dvcs_cut[topology.key][variable.branch])),
                "pi0_mc_entries_in_range": float(np.sum(pi0_cut[topology.key][variable.branch])),
                "morph_type": result.morph_label,
                "shift_or_log_shift": result.shift, "shift_err": result.shift_err,
                "shift_bins": result.shift / bin_width if result.success and additive else math.nan,
                "sigma_or_log_sigma": result.sigma_add, "sigma_err": result.sigma_add_err,
                "sigma_bins": result.sigma_add / bin_width if result.success and additive else math.nan,
                "fit_bins_used": int(np.count_nonzero(result.fit_mask)) if result.fit_mask is not None else 0,
                "excluded_data_counts": result.excluded_data_counts,
                "excluded_model_counts": result.excluded_model_counts,
                "excluded_excess_counts": result.excluded_excess_counts,
                "excluded_excess_fraction_of_data": result.excluded_excess_fraction,
                "variable_poisson_deviance": result.deviance, "variable_ndf": result.ndf,
                "variable_deviance_per_ndf": result.deviance / result.ndf if result.success and result.ndf > 0 else math.nan,
            })
        # endfor

        fit_path = output_dir / "template_fits" / f"exclusivity_template_fit_{period.key}_{topology.key.lower()}.png"
        draw_fit_canvas(
            fit_path, period, topology, data_cut[topology.key],
            dvcs_cut[topology.key], pi0_cut[topology.key],
            cut_selected, fit_results, shared_summary, log_y, dpi,
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

    periods = selected_periods(args.period)
    topologies = selected_topologies(args.topology)
    output_dir = Path(args.output_dir)
    all_rows: List[Dict[str, object]] = []

    log(f"ROOT I/O backend: {io_backend()}")
    log(f"Producing {3 * len(periods) * len(topologies)} canvases: uncut shapes, Mx2_1-cut shapes and template fits")
    log(
        "Selection: topology, (-t1) < 1.0, open_angle_ep2 > 5 deg, "
        "distinct FD sectors; no exclusivity cuts"
    )
    log(
        "Fit: one shared f_pi0 per period/topology; additive morphs for signed "
        "variables, log-space morphs for pTmiss and theta_gamma_gamma; "
        "Mx2_1 hard cut applied before second-stage shapes and fits; xF edge bins excluded; "
        "masked excess reported separately"
    )

    for period in periods:
        all_rows.extend(process_period(
            period, topologies, output_dir, args.step_size,
            args.max_shift_bins, args.max_smear_bins, args.fit_min_counts,
            args.log_y, args.dpi,
        ))
    # endfor

    csv_path = output_dir / "template_fits" / "fit_results.csv"
    write_results_csv(csv_path, all_rows)
    log(f"Wrote {csv_path}")
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
