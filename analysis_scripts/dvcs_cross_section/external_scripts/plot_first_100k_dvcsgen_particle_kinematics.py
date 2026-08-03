#!/usr/bin/env python3
"""
plot_first_100k_dvcsgen_particle_kinematics_v4_embedded_reconstruction.py

Compare the first 100,000 entries from the generated DVCSGEN and AAOGEN ROOT
files.

DVCSGEN:
    p1_p, p1_theta  -> generated proton
    p2_p, p2_theta  -> generated real photon

AAOGEN:
    p1_p, p1_theta  -> generated proton
    photon 1 and photon 2 are reconstructed with an exact embedded copy of the current
    reconstruct_native_eppi0_photons() implementation from
    derive_photon_efficiency_scale_factors.py.

No independent isotropic pi0 decay is generated in this script.

The embedded reconstruction uses:
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
import math
from pathlib import Path
from types import SimpleNamespace
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


def _unit_vector(theta: float, phi: float) -> np.ndarray:
    st = math.sin(theta)
    return np.asarray([st * math.cos(phi), st * math.sin(phi), math.cos(theta)], dtype=float)

def _orthonormal_basis_about(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return two unit vectors transverse to a unit axis."""
    axis = np.asarray(axis, dtype=float)
    if abs(axis[2]) < 0.9:
        reference = np.asarray([0.0, 0.0, 1.0], dtype=float)
    else:
        reference = np.asarray([1.0, 0.0, 0.0], dtype=float)
    u = np.cross(reference, axis)
    norm_u = float(np.linalg.norm(u))
    if norm_u <= 1.0e-14:
        reference = np.asarray([0.0, 1.0, 0.0], dtype=float)
        u = np.cross(reference, axis)
        norm_u = float(np.linalg.norm(u))
    u /= norm_u
    v = np.cross(axis, u)
    v /= float(np.linalg.norm(v))
    return u, v

def _detector_compatible(theta_rad: float, detector: int, args: argparse.Namespace) -> bool:
    theta_deg = math.degrees(theta_rad)
    if detector == 0:
        return args.ft_theta_min <= theta_deg < args.ft_theta_max
    if detector == 1:
        return args.fd_theta_min <= theta_deg < args.fd_theta_max
    return False

def reconstruct_native_eppi0_photons(
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    pi0_p: np.ndarray,
    pi0_theta: np.ndarray,
    pi0_phi: np.ndarray,
    opening1: np.ndarray,
    opening2: np.ndarray,
    pi0_mass: np.ndarray,
    detector1: np.ndarray,
    detector2: np.ndarray,
    args: argparse.Namespace,
    active_mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, ...]:
    """Recover both detector-compatible cone-mirror solutions when available.

    The first solution index is the nominal best-score solution.  The second is
    the alternate mirror.  Storing both solutions allows the downstream
    ``best``, ``half_weight`` and ``unique_sector`` prescriptions to be applied
    without rerunning the cone reconstruction.
    """
    n_events = len(e_theta)
    g1_E = np.full(n_events, np.nan, dtype=float)
    g2_E = np.full(n_events, np.nan, dtype=float)
    solution_g1_theta = np.full((n_events, 2), np.nan, dtype=float)
    solution_g1_phi = np.full((n_events, 2), np.nan, dtype=float)
    solution_g2_theta = np.full((n_events, 2), np.nan, dtype=float)
    solution_g2_phi = np.full((n_events, 2), np.nan, dtype=float)
    solution_closure = np.full((n_events, 2), np.nan, dtype=float)
    solution_score = np.full((n_events, 2), np.nan, dtype=float)
    closure = np.full(n_events, np.nan, dtype=float)
    ambiguity = np.full(n_events, np.nan, dtype=float)
    transverse_mismatch = np.full(n_events, np.nan, dtype=float)
    energy_fraction_g1 = np.full(n_events, np.nan, dtype=float)
    n_solutions = np.zeros(n_events, dtype=np.int16)
    input_finite = np.zeros(n_events, dtype=bool)
    longitudinal_solution = np.zeros(n_events, dtype=bool)
    transverse_solution = np.zeros(n_events, dtype=bool)
    detector_solution = np.zeros(n_events, dtype=bool)
    closure_pass = np.zeros(n_events, dtype=bool)

    if active_mask is None:
        active_indices = range(n_events)
    else:
        active_indices = np.flatnonzero(np.asarray(active_mask, dtype=bool))

    for i in active_indices:
        values = (
            e_theta[i], e_phi[i], pi0_p[i], pi0_theta[i], pi0_phi[i],
            opening1[i], opening2[i], pi0_mass[i], detector1[i], detector2[i],
        )
        if not all(math.isfinite(float(v)) for v in values):
            continue
        if pi0_p[i] <= 0.0 or pi0_mass[i] <= 0.0:
            continue
        if not (0.0 < opening1[i] < math.pi and 0.0 < opening2[i] < math.pi):
            continue
        input_finite[i] = True

        ehat = _unit_vector(float(e_theta[i]), float(e_phi[i]))
        npi = _unit_vector(float(pi0_theta[i]), float(pi0_phi[i]))
        pvec = float(pi0_p[i]) * npi
        epi = math.sqrt(float(pi0_p[i]) ** 2 + float(pi0_mass[i]) ** 2)
        p_parallel = float(np.dot(pvec, ehat))

        c1 = math.cos(float(opening1[i]))
        c2 = math.cos(float(opening2[i]))
        denom = c1 - c2
        if abs(denom) <= args.cone_energy_denominator_min:
            continue
        E1 = (p_parallel - epi * c2) / denom
        E2 = epi - E1
        if E1 <= args.cone_min_photon_energy or E2 <= args.cone_min_photon_energy:
            continue
        longitudinal_solution[i] = True
        energy_fraction_g1[i] = E1 / epi
        g1_E[i], g2_E[i] = E1, E2

        u, v = _orthonormal_basis_about(ehat)
        p_perp_vec = pvec - p_parallel * ehat
        px = float(np.dot(p_perp_vec, u))
        py = float(np.dot(p_perp_vec, v))
        R = math.hypot(px, py)
        psi = math.atan2(py, px) if R > 1.0e-14 else 0.0
        A = E1 * math.sin(float(opening1[i]))
        B = E2 * math.sin(float(opening2[i]))
        if A <= 1.0e-14 or B <= 1.0e-14:
            continue

        lower = abs(A - B)
        upper = A + B
        target_R = min(max(R, lower), upper)
        transverse_mismatch[i] = abs(R - target_R) / max(epi, 1.0e-12)
        if target_R <= 1.0e-14:
            delta1 = 0.0
            delta2 = math.pi if A >= B else 0.0
        else:
            cos_delta1 = (target_R * target_R + A * A - B * B) / (2.0 * target_R * A)
            cos_delta2 = (target_R * target_R + B * B - A * A) / (2.0 * target_R * B)
            delta1 = math.acos(float(np.clip(cos_delta1, -1.0, 1.0)))
            delta2 = math.acos(float(np.clip(cos_delta2, -1.0, 1.0)))
        transverse_solution[i] = True

        candidates: List[Tuple[float, float, float, float, float, float, float]] = []
        for sign in (+1.0, -1.0):
            beta1 = psi + sign * delta1
            beta2 = psi - sign * delta2
            n1 = c1 * ehat + math.sin(float(opening1[i])) * (math.cos(beta1) * u + math.sin(beta1) * v)
            n2 = c2 * ehat + math.sin(float(opening2[i])) * (math.cos(beta2) * u + math.sin(beta2) * v)
            n1 /= float(np.linalg.norm(n1))
            n2 /= float(np.linalg.norm(n2))
            th1 = math.acos(float(np.clip(n1[2], -1.0, 1.0)))
            ph1 = math.atan2(float(n1[1]), float(n1[0]))
            th2 = math.acos(float(np.clip(n2[2], -1.0, 1.0)))
            ph2 = math.atan2(float(n2[1]), float(n2[0]))
            if not (_detector_compatible(th1, int(detector1[i]), args) and
                    _detector_compatible(th2, int(detector2[i]), args)):
                continue
            four_residual = np.asarray([E1 + E2 - epi, *(E1 * n1 + E2 * n2 - pvec)], dtype=float)
            residual = float(np.linalg.norm(four_residual) / max(epi, 1.0e-12))
            mass2 = 2.0 * E1 * E2 * (1.0 - float(np.dot(n1, n2)))
            mass_residual = abs(mass2 - float(pi0_mass[i]) ** 2) / max(float(pi0_mass[i]) ** 2, 1.0e-12)
            score = math.hypot(residual, args.cone_mass_residual_weight * mass_residual)
            candidates.append((score, residual, mass_residual, th1, ph1, th2, ph2))

        candidates.sort(key=lambda item: item[0])
        n_solutions[i] = min(len(candidates), 2)
        if not candidates:
            continue
        detector_solution[i] = True
        for j, candidate in enumerate(candidates[:2]):
            solution_score[i, j] = candidate[0]
            solution_closure[i, j] = candidate[1]
            solution_g1_theta[i, j], solution_g1_phi[i, j] = candidate[3], candidate[4]
            solution_g2_theta[i, j], solution_g2_phi[i, j] = candidate[5], candidate[6]
        closure[i] = solution_closure[i, 0]
        ambiguity[i] = (solution_score[i, 1] - solution_score[i, 0]
                        if n_solutions[i] > 1 else math.inf)
        closure_pass[i] = solution_closure[i, 0] <= args.pass_photon_closure_max

    return (
        g1_E, g2_E, solution_g1_theta, solution_g1_phi,
        solution_g2_theta, solution_g2_phi, solution_closure, solution_score,
        closure, ambiguity, transverse_mismatch, energy_fraction_g1,
        n_solutions, input_finite, longitudinal_solution, transverse_solution,
        detector_solution, closure_pass,
    )

def read_aaogen_with_main_reconstruction(
    path: Path,
    tree_name: str,
    max_entries: int,
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

    result = reconstruct_native_eppi0_photons(
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
    for path in (args.dvcsgen, args.aaogen):
        if not path.exists():
            raise FileNotFoundError(path)
        # endif
    # endfor


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
        "AAOGEN photons reconstructed with the embedded main-analysis helper: "
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

    print("Using embedded reconstruction copied from the current main analysis script.")
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
