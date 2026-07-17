#!/usr/bin/env python3
"""
plot_exclusivity_data_dvcs_pi0_mc_v1.py

Compare uncut DVCS-candidate exclusivity-variable shapes for:

  1. DVCS data (epgamma skim)
  2. reconstructed DVCS signal MC
  3. reconstructed eppi0 MC misidentified and processed as DVCS

The script produces one 4x2 canvas for every run-period/topology combination:

  5 run periods x 3 topologies = 15 PNG files

Only the requested detector topology is enforced. No global cuts and no
exclusivity cuts are applied by this script.

Default normalization:

  Each sample is independently normalized to unit integral inside the plotted
  range. This is the appropriate first comparison of the SHAPES. It does not
  assume a physical relative normalization between DVCS MC and eppi0 background
  MC.

Optional mixed template:

  Pass --pi0-fraction F, with 0 <= F <= 1, to additionally draw

      (1-F) * normalized DVCS MC + F * normalized eppi0-as-DVCS MC.

  This F is a shape-mixture fraction inside the plotted range. It is not
  automatically the final contamination fraction used by the cross-section
  analysis.

Dependencies:

  Python 3, numpy, matplotlib and either uproot or PyROOT

Example:

  python3 external_scripts/plot_exclusivity_data_dvcs_pi0_mc_v1.py

  python3 external_scripts/plot_exclusivity_data_dvcs_pi0_mc_v1.py \
      --pi0-fraction 0.15 --log-y
"""

import argparse
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
    VariableConfig("xF", r"$x_F$", 100, -0.4, 0.2),
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
        100,
        0.0,
        3.0,
        aliases=("Mx2_egamma", "Mx2_gamma", "Mx2_pi0", "Mx2_x2"),
    ),
)


SAMPLE_LABELS = {
    "data": "DVCS data",
    "dvcs_mc": "DVCS MC",
    "pi0_mc": r"$e\pi^0$ MC as DVCS",
}


SAMPLE_COLORS = {
    "data": "black",
    "dvcs_mc": "tab:blue",
    "pi0_mc": "tab:red",
    "mixture": "tab:green",
}


def log(message: str) -> None:
    now = time.strftime("%H:%M:%S")
    print(f"[{now}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot uncut exclusivity-variable shapes for DVCS data, reconstructed "
            "DVCS MC and eppi0 MC reconstructed as DVCS."
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
        "--pi0-fraction",
        type=float,
        default=None,
        help=(
            "Optionally draw a mixed normalized-MC template with this eppi0 fraction. "
            "Must be between 0 and 1."
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

    required_topology = {"detector1", "detector2"}
    missing_topology = sorted(required_topology - branches)
    if missing_topology:
        raise KeyError(
            f"Missing topology branches in {path}: {', '.join(missing_topology)}"
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


def empty_histograms() -> Dict[str, np.ndarray]:
    return {
        variable.branch: np.zeros(variable.bins, dtype=np.float64)
        for variable in VARIABLES
    }


def fill_histograms_uproot(
    path: str,
    topologies: Sequence[TopologyConfig],
    step_size: str,
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, int]]:
    """Read one ROOT file once and fill all requested topology/variable histograms."""

    require_input_file(path)
    resolved = resolve_variable_branches(path)

    expressions = ["detector1", "detector2"] + sorted(set(resolved.values()))
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

        for topology in topologies:
            mask = (
                (detector1 == topology.detector1)
                & (detector2 == topology.detector2)
            )
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
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, int]]:
    backend = io_backend()

    if backend == "uproot":
        return fill_histograms_uproot(path, topologies, step_size)
    # endif

    return fill_histograms_pyroot(path, topologies, step_size)


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


def draw_canvas(
    output_path: Path,
    period: PeriodConfig,
    topology: TopologyConfig,
    data_histograms: Mapping[str, np.ndarray],
    dvcs_histograms: Mapping[str, np.ndarray],
    pi0_histograms: Mapping[str, np.ndarray],
    selected_counts: Mapping[str, int],
    pi0_fraction: Optional[float],
    log_y: bool,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(18.0, 8.8))
    flat_axes = axes.ravel()

    for axis, variable in zip(flat_axes, VARIABLES):
        edges, centers = bin_geometry(variable)

        data_density, data_error = normalize_density(
            data_histograms[variable.branch], variable
        )
        dvcs_density, _ = normalize_density(
            dvcs_histograms[variable.branch], variable
        )
        pi0_density, _ = normalize_density(
            pi0_histograms[variable.branch], variable
        )

        axis.errorbar(
            centers,
            data_density,
            yerr=data_error,
            fmt="o",
            markersize=2.4,
            linewidth=0.8,
            capsize=0.0,
            color=SAMPLE_COLORS["data"],
            label=SAMPLE_LABELS["data"],
            zorder=4,
        )
        axis.stairs(
            dvcs_density,
            edges,
            color=SAMPLE_COLORS["dvcs_mc"],
            linewidth=1.8,
            label=SAMPLE_LABELS["dvcs_mc"],
            zorder=3,
        )
        axis.stairs(
            pi0_density,
            edges,
            color=SAMPLE_COLORS["pi0_mc"],
            linewidth=1.8,
            label=SAMPLE_LABELS["pi0_mc"],
            zorder=2,
        )

        mixture = None
        if pi0_fraction is not None:
            mixture = (1.0 - pi0_fraction) * dvcs_density + pi0_fraction * pi0_density
            axis.stairs(
                mixture,
                edges,
                color=SAMPLE_COLORS["mixture"],
                linewidth=1.8,
                linestyle="--",
                label=rf"MC mixture ($f_{{\pi^0}}={pi0_fraction:.3f}$)",
                zorder=5,
            )
        # endif

        axis.set_xlim(variable.xmin, variable.xmax)
        axis.set_xlabel(variable.label)
        axis.set_ylabel("unit-normalized density")
        axis.grid(axis="y", alpha=0.25)

        if log_y:
            arrays_for_floor = [data_density, dvcs_density, pi0_density]
            if mixture is not None:
                arrays_for_floor.append(mixture)
            # endif

            floor = positive_y_floor(*arrays_for_floor)
            if floor is not None:
                axis.set_yscale("log")
                axis.set_ylim(bottom=floor)
            # endif
        else:
            axis.set_ylim(bottom=0.0)
        # endif
    # endfor

    handles, labels = flat_axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=len(labels),
        frameon=False,
    )

    fig.suptitle(
        f"Exclusivity shapes before global/exclusivity cuts: "
        f"{period.label}, {topology.label}\n"
        f"topology-selected entries: data={selected_counts['data']:,}, "
        f"DVCS MC={selected_counts['dvcs_mc']:,}, "
        f"$e\\pi^0$ MC as DVCS={selected_counts['pi0_mc']:,}",
        fontsize=15,
        y=0.995,
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
    pi0_fraction: Optional[float],
    log_y: bool,
    dpi: int,
) -> None:
    log(f"Starting period {period.label}")

    data_hists, data_counts = fill_histograms_for_file(
        period.data_file, topologies, step_size
    )
    dvcs_hists, dvcs_counts = fill_histograms_for_file(
        period.dvcs_mc_file, topologies, step_size
    )
    pi0_hists, pi0_counts = fill_histograms_for_file(
        period.pi0_as_dvcs_mc_file, topologies, step_size
    )

    for topology in topologies:
        output_path = output_dir / (
            f"exclusivity_shapes_{period.key}_{topology.key.lower()}.png"
        )

        draw_canvas(
            output_path=output_path,
            period=period,
            topology=topology,
            data_histograms=data_hists[topology.key],
            dvcs_histograms=dvcs_hists[topology.key],
            pi0_histograms=pi0_hists[topology.key],
            selected_counts={
                "data": data_counts[topology.key],
                "dvcs_mc": dvcs_counts[topology.key],
                "pi0_mc": pi0_counts[topology.key],
            },
            pi0_fraction=pi0_fraction,
            log_y=log_y,
            dpi=dpi,
        )
        log(f"Wrote {output_path}")
    # endfor


def main() -> int:
    args = parse_args()

    if args.pi0_fraction is not None and not (0.0 <= args.pi0_fraction <= 1.0):
        raise ValueError("--pi0-fraction must be between 0 and 1 inclusive.")
    # endif

    periods = selected_periods(args.period)
    topologies = selected_topologies(args.topology)
    output_dir = Path(args.output_dir)

    log(f"ROOT I/O backend: {io_backend()}")
    log(
        f"Producing {len(periods) * len(topologies)} canvases "
        f"for {len(periods)} periods and {len(topologies)} topologies"
    )
    log("Selection: detector topology only; no global or exclusivity cuts")
    log("Normalization: each sample has unit integral inside each plotted range")

    if args.pi0_fraction is not None:
        log(f"Additional mixed template enabled with f_pi0={args.pi0_fraction:.6g}")
    # endif

    for period in periods:
        process_period(
            period=period,
            topologies=topologies,
            output_dir=output_dir,
            step_size=args.step_size,
            pi0_fraction=args.pi0_fraction,
            log_y=args.log_y,
            dpi=args.dpi,
        )
    # endfor

    log("All requested canvases completed")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    # endtry
#endif
