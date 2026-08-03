#!/usr/bin/env python3
"""
plot_first_100k_dvcsgen_particle_kinematics_v3_import_fix.py

Compare the first 100,000 entries from the generated DVCSGEN and AAOGEN ROOT
files.

DVCSGEN:
    p1_p, p1_theta  -> generated proton
    p2_p, p2_theta  -> generated real photon

AAOGEN:
    p1_p, p1_theta  -> generated proton
    photon 1 and photon 2 are reconstructed with the exact
    reconstruct_native_eppi0_photons() implementation imported from
    derive_photon_efficiency_scale_factors.py.

No independent isotropic pi0 decay is generated in this script.

The imported reconstruction uses:
    e_theta, e_phi
    pi0_p, pi0_theta, pi0_phi
    open_angle_egamma1, open_angle_egamma2
    Mh_gammagamma
    detector_gamma1, detector_gamma2

and returns the best-ranked physical daughter solution after the same analytic
energy solution, transverse-cone construction, detector-compatibility test,
four-vector closure calculation, and diphoton-mass ranking used by the main
photon-efficiency analysis.

Dependencies:
    uproot, numpy, matplotlib
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import uproot


DEFAULT_DVCSGEN = Path(
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
    "dvcsgen/gen_dvcsgen_rga_fa18_out_50nA_10604MeV.root"
)
DEFAULT_AAOGEN = Path(
    "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/"
    "hipo_files/gen_aaogen_norad_fa18_out_50nA_10604MeV.root"
)
DEFAULT_MAIN_SCRIPT = Path(__file__).with_name(
    "derive_photon_efficiency_scale_factors.py"
)
DEFAULT_OUTPUT = Path(
    "first_100k_dvcsgen_vs_aaogen_reconstructed_photons.png"
)
DEFAULT_TREE = "PhysicsEvents"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Overlay DVCSGEN proton/photon kinematics with AAOGEN proton and "
            "pi0 daughter photons reconstructed by the main efficiency script."
        )
    )
    parser.add_argument(
        "--dvcsgen",
        type=Path,
        default=DEFAULT_DVCSGEN,
    )
    parser.add_argument(
        "--aaogen",
        type=Path,
        default=DEFAULT_AAOGEN,
    )
    parser.add_argument(
        "--main-script",
        type=Path,
        default=DEFAULT_MAIN_SCRIPT,
        help=(
            "Path to derive_photon_efficiency_scale_factors.py. The exact "
            "reconstruct_native_eppi0_photons() function is imported from it."
        ),
    )
    parser.add_argument(
        "--tree",
        default=DEFAULT_TREE,
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=100_000,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
    )

    # These defaults are copied from the main script's current CLI defaults.
    # They may be overridden here for an exact configuration comparison.
    parser.add_argument("--ft-theta-min", type=float, default=2.5)
    parser.add_argument("--ft-theta-max", type=float, default=5.0)
    parser.add_argument("--fd-theta-min", type=float, default=5.0)
    parser.add_argument("--fd-theta-max", type=float, default=35.0)
    parser.add_argument(
        "--pass-photon-closure-max",
        type=float,
        default=0.03,
    )
    parser.add_argument(
        "--cone-energy-denominator-min",
        type=float,
        default=1.0e-5,
    )
    parser.add_argument(
        "--cone-min-photon-energy",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--cone-mass-residual-weight",
        type=float,
        default=0.05,
    )
    return parser.parse_args()


def import_main_module(path: Path) -> ModuleType:
    if not path.exists():
        raise FileNotFoundError(
            f"Main photon-efficiency script does not exist: {path}"
        )
    # endif

    spec = importlib.util.spec_from_file_location(
        "_photon_efficiency_main_for_quick_plot",
        path,
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create import specification for {path}")
    # endif

    module = importlib.util.module_from_spec(spec)

    # Register the dynamically created module before executing it. Python 3.13
    # dataclasses inspect sys.modules while class decorators are being applied.
    # Without this registration, importing the main analysis script fails with
    # an AttributeError inside dataclasses.py.
    sys.modules[spec.name] = module

    try:
        spec.loader.exec_module(module)
    except Exception:
        # Do not leave a partially initialized module cached after a failed
        # import.
        sys.modules.pop(spec.name, None)
        raise
    # endtry

    required = "reconstruct_native_eppi0_photons"
    if not hasattr(module, required):
        raise AttributeError(
            f"{path} does not define {required}()."
        )
    # endif
    return module


def require_tree(
    root_file: uproot.ReadOnlyDirectory,
    tree_name: str,
):
    if tree_name not in root_file:
        raise KeyError(
            f"Tree '{tree_name}' was not found. Available keys: "
            + ", ".join(root_file.keys())
        )
    # endif
    return root_file[tree_name]


def resolve_branch(
    available: set[str],
    aliases: Sequence[str],
    logical_name: str,
) -> str:
    for alias in aliases:
        if alias in available:
            return alias
        # endif
    # endfor

    raise KeyError(
        f"Could not resolve {logical_name}. Tried {tuple(aliases)}.\n"
        "Available branches:\n  "
        + ", ".join(sorted(available))
    )


def finite(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    return values[np.isfinite(values)]


def read_dvcsgen(
    path: Path,
    tree_name: str,
    max_entries: int,
) -> Mapping[str, np.ndarray]:
    with uproot.open(path) as root_file:
        tree = require_tree(root_file, tree_name)
        available = set(tree.keys())

        branches = {
            "p1_p": resolve_branch(
                available,
                ("p1_p", "proton_p", "p_p"),
                "DVCSGEN proton momentum",
            ),
            "p1_theta": resolve_branch(
                available,
                ("p1_theta", "proton_theta", "p_theta"),
                "DVCSGEN proton theta",
            ),
            "p2_p": resolve_branch(
                available,
                ("p2_p", "g_p", "gamma_p", "photon_p"),
                "DVCSGEN photon momentum",
            ),
            "p2_theta": resolve_branch(
                available,
                ("p2_theta", "g_theta", "gamma_theta", "photon_theta"),
                "DVCSGEN photon theta",
            ),
        }

        stop = min(max_entries, int(tree.num_entries))
        arrays = tree.arrays(
            list(branches.values()),
            entry_start=0,
            entry_stop=stop,
            library="np",
        )
    # endwith

    return {
        "proton_p": finite(arrays[branches["p1_p"]]),
        "proton_theta_deg": finite(
            np.degrees(arrays[branches["p1_theta"]])
        ),
        "gamma_p": finite(arrays[branches["p2_p"]]),
        "gamma_theta_deg": finite(
            np.degrees(arrays[branches["p2_theta"]])
        ),
        "entries_read": np.asarray([stop], dtype=np.int64),
    }


AAOGEN_ALIASES: Mapping[str, Sequence[str]] = {
    "e_theta": ("e_theta",),
    "e_phi": ("e_phi",),
    "p1_p": ("p1_p", "p_p", "proton_p"),
    "p1_theta": ("p1_theta", "p_theta", "proton_theta"),
    "pi0_p": ("p2_p", "pi0_p"),
    "pi0_theta": ("p2_theta", "pi0_theta"),
    "pi0_phi": ("p2_phi", "pi0_phi"),
    "opening1": ("open_angle_egamma1",),
    "opening2": ("open_angle_egamma2",),
    "mass": ("Mh_gammagamma", "Mh"),
    "detector1": (
        "detector_gamma1",
        "gamma1_detector",
        "detector1_gamma",
    ),
    "detector2": (
        "detector_gamma2",
        "gamma2_detector",
        "detector2_gamma",
    ),
}


def read_aaogen_with_main_reconstruction(
    path: Path,
    tree_name: str,
    max_entries: int,
    main_module: ModuleType,
    reconstruction_args: argparse.Namespace,
) -> Mapping[str, np.ndarray]:
    with uproot.open(path) as root_file:
        tree = require_tree(root_file, tree_name)
        available = set(tree.keys())

        resolved = {
            logical: resolve_branch(available, aliases, logical)
            for logical, aliases in AAOGEN_ALIASES.items()
        }

        stop = min(max_entries, int(tree.num_entries))
        arrays = tree.arrays(
            sorted(set(resolved.values())),
            entry_start=0,
            entry_stop=stop,
            library="np",
        )
    # endwith

    def array(logical: str, dtype=float) -> np.ndarray:
        return np.asarray(arrays[resolved[logical]], dtype=dtype)

    active_mask = np.ones(stop, dtype=bool)

    result = main_module.reconstruct_native_eppi0_photons(
        array("e_theta"),
        array("e_phi"),
        array("pi0_p"),
        array("pi0_theta"),
        array("pi0_phi"),
        np.deg2rad(array("opening1")),
        np.deg2rad(array("opening2")),
        array("mass"),
        array("detector1", int),
        array("detector2", int),
        reconstruction_args,
        active_mask=active_mask,
    )

    (
        gamma1_energy,
        gamma2_energy,
        solution_gamma1_theta,
        _solution_gamma1_phi,
        solution_gamma2_theta,
        _solution_gamma2_phi,
        _solution_closure,
        _solution_score,
        closure,
        ambiguity,
        _transverse_mismatch,
        _energy_fraction_gamma1,
        number_of_solutions,
        input_finite,
        longitudinal_solution,
        transverse_solution,
        detector_solution,
        closure_pass,
    ) = result

    best_gamma1_theta = solution_gamma1_theta[:, 0]
    best_gamma2_theta = solution_gamma2_theta[:, 0]

    valid = (
        input_finite
        & longitudinal_solution
        & transverse_solution
        & detector_solution
        & closure_pass
        & (number_of_solutions >= 1)
        & np.isfinite(gamma1_energy)
        & np.isfinite(gamma2_energy)
        & np.isfinite(best_gamma1_theta)
        & np.isfinite(best_gamma2_theta)
    )

    return {
        "proton_p": finite(array("p1_p")),
        "proton_theta_deg": finite(np.degrees(array("p1_theta"))),
        "gamma1_p": gamma1_energy[valid],
        "gamma1_theta_deg": np.degrees(best_gamma1_theta[valid]),
        "gamma2_p": gamma2_energy[valid],
        "gamma2_theta_deg": np.degrees(best_gamma2_theta[valid]),
        "entries_read": np.asarray([stop], dtype=np.int64),
        "valid_reconstructions": np.asarray(
            [np.count_nonzero(valid)],
            dtype=np.int64,
        ),
        "input_finite_count": np.asarray(
            [np.count_nonzero(input_finite)],
            dtype=np.int64,
        ),
        "longitudinal_count": np.asarray(
            [np.count_nonzero(longitudinal_solution)],
            dtype=np.int64,
        ),
        "transverse_count": np.asarray(
            [np.count_nonzero(transverse_solution)],
            dtype=np.int64,
        ),
        "detector_count": np.asarray(
            [np.count_nonzero(detector_solution)],
            dtype=np.int64,
        ),
        "closure_count": np.asarray(
            [np.count_nonzero(closure_pass)],
            dtype=np.int64,
        ),
        "median_closure": np.asarray(
            [float(np.nanmedian(closure[detector_solution]))],
            dtype=float,
        ),
        "median_ambiguity": np.asarray(
            [
                float(
                    np.nanmedian(
                        ambiguity[
                            detector_solution
                            & np.isfinite(ambiguity)
                        ]
                    )
                )
            ],
            dtype=float,
        ),
    }


def unit_area_histogram(
    values: np.ndarray,
    bins: int,
    plot_range: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    counts, edges = np.histogram(
        values,
        bins=bins,
        range=plot_range,
    )
    total = float(np.sum(counts))
    if total <= 0.0:
        return np.zeros_like(counts, dtype=float), edges
    # endif
    return counts.astype(float) / total, edges


def main() -> int:
    args = parse_args()

    if args.max_entries <= 0:
        raise ValueError("--max-entries must be positive.")
    # endif
    for path in (args.dvcsgen, args.aaogen, args.main_script):
        if not path.exists():
            raise FileNotFoundError(path)
        # endif
    # endfor

    main_module = import_main_module(args.main_script)

    # The imported function accesses only these reconstruction settings.
    reconstruction_args = SimpleNamespace(
        ft_theta_min=args.ft_theta_min,
        ft_theta_max=args.ft_theta_max,
        fd_theta_min=args.fd_theta_min,
        fd_theta_max=args.fd_theta_max,
        pass_photon_closure_max=args.pass_photon_closure_max,
        cone_energy_denominator_min=args.cone_energy_denominator_min,
        cone_min_photon_energy=args.cone_min_photon_energy,
        cone_mass_residual_weight=args.cone_mass_residual_weight,
    )

    dvcs = read_dvcsgen(
        args.dvcsgen,
        args.tree,
        args.max_entries,
    )
    aaogen = read_aaogen_with_main_reconstruction(
        args.aaogen,
        args.tree,
        args.max_entries,
        main_module,
        reconstruction_args,
    )

    panels = (
        (
            r"Proton momentum (GeV)",
            (0.0, 5.0),
            (
                ("DVCSGEN proton", dvcs["proton_p"], "tab:blue", "-"),
                ("AAOGEN proton", aaogen["proton_p"], "tab:red", "-"),
            ),
        ),
        (
            r"Proton $\theta$ (deg)",
            (0.0, 70.0),
            (
                (
                    "DVCSGEN proton",
                    dvcs["proton_theta_deg"],
                    "tab:blue",
                    "-",
                ),
                (
                    "AAOGEN proton",
                    aaogen["proton_theta_deg"],
                    "tab:red",
                    "-",
                ),
            ),
        ),
        (
            r"Photon momentum (GeV)",
            (0.0, 10.0),
            (
                ("DVCSGEN photon", dvcs["gamma_p"], "tab:blue", "-"),
                (
                    r"AAOGEN reconstructed $\gamma_1$",
                    aaogen["gamma1_p"],
                    "tab:red",
                    "-",
                ),
                (
                    r"AAOGEN reconstructed $\gamma_2$",
                    aaogen["gamma2_p"],
                    "tab:red",
                    "--",
                ),
            ),
        ),
        (
            r"Photon $\theta$ (deg)",
            (0.0, 40.0),
            (
                (
                    "DVCSGEN photon",
                    dvcs["gamma_theta_deg"],
                    "tab:blue",
                    "-",
                ),
                (
                    r"AAOGEN reconstructed $\gamma_1$",
                    aaogen["gamma1_theta_deg"],
                    "tab:red",
                    "-",
                ),
                (
                    r"AAOGEN reconstructed $\gamma_2$",
                    aaogen["gamma2_theta_deg"],
                    "tab:red",
                    "--",
                ),
            ),
        ),
    )

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for axis, (x_label, plot_range, samples) in zip(
        axes.ravel(),
        panels,
    ):
        for sample_label, values, color, linestyle in samples:
            normalized, edges = unit_area_histogram(
                values,
                120,
                plot_range,
            )
            axis.stairs(
                normalized,
                edges,
                color=color,
                linestyle=linestyle,
                linewidth=1.5,
                label=f"{sample_label} ({values.size:,})",
            )
        # endfor

        axis.set_xlabel(x_label)
        axis.set_ylabel("fraction / bin")
        axis.set_xlim(*plot_range)
        axis.set_ylim(bottom=0.0)
        axis.grid(axis="y", alpha=0.25)
        axis.legend(frameon=False)
    # endfor

    entries = int(aaogen["entries_read"][0])
    valid = int(aaogen["valid_reconstructions"][0])
    fig.suptitle(
        f"First {args.max_entries:,} generated entries\n"
        "AAOGEN photons reconstructed with the exact main-analysis helper: "
        f"{valid:,}/{entries:,} valid closure-passing events",
        fontsize=14,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        args.output,
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"Imported reconstruction from: {args.main_script}")
    print(f"DVCSGEN entries read: {int(dvcs['entries_read'][0]):,}")
    print(f"AAOGEN entries read: {entries:,}")
    print(
        "AAOGEN reconstruction cutflow: "
        f"finite={int(aaogen['input_finite_count'][0]):,}, "
        f"longitudinal={int(aaogen['longitudinal_count'][0]):,}, "
        f"transverse={int(aaogen['transverse_count'][0]):,}, "
        f"detector-compatible={int(aaogen['detector_count'][0]):,}, "
        f"closure-pass={int(aaogen['closure_count'][0]):,}, "
        f"final-valid={valid:,}"
    )
    print(
        f"Median normalized closure: "
        f"{float(aaogen['median_closure'][0]):.6g}"
    )
    print(
        f"Median mirror ambiguity: "
        f"{float(aaogen['median_ambiguity'][0]):.6g}"
    )
    print(f"Wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
