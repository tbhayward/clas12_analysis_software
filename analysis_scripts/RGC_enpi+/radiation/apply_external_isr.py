#!/usr/bin/env python3
"""
apply_external_isr_v4.py

Apply one additional, deterministic external-ISR draw to an existing
internal-ISR RGC e pi+ ROOT tree.

The input tree is expected to contain measured/corrected electron and pion
four-vector variables plus the already sampled internal-ISR photon energy.
Every existing branch is preserved. Beam-energy-dependent derived branches
are replaced consistently, and external-radiation provenance branches are
added.

Required packages:
    python3 -m pip install numpy awkward uproot matplotlib

Example:
    python3 apply_external_isr_v2.py \
        rgc_fa22_inb_NH3_epi+_ISR_mom_corrections.root \
        --geometry-json external_radiation_geometry_v1.json

Important:
  * The script requires an internal ISR photon-energy branch named either
    Egamma_internal or Egamma. It will not guess this energy from the other
    branches.
  * Deterministic sampling uses runnum, the absolute tree-entry index, the
    input-file identity, and --seed-ensemble. This remains reproducible for an
    unchanged input file even when evnum is unavailable.
  * Geometry records marked verified=false are rejected unless
    --allow-unverified-geometry is supplied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot


TREE_DEFAULT = "PhysicsEvents"
STEP_SIZE_DEFAULT = "250 MB"
MODEL_VERSION = "external-isr-v5"

ELECTRON_MASS_GEV = 0.00051099895
PROTON_MASS_GEV = 0.9382720813
PION_MASS_GEV = 0.13957039
NEUTRON_MASS_GEV = 0.9395654133
TWO_PI = 2.0 * math.pi

RECALCULATED_BRANCHES = (
    "Q2", "W", "Mx2", "x", "y", "t", "tmin", "tprime",
    "z", "xF", "pT", "phi",
    "DepA", "DepB", "DepC", "DepV", "DepW",
)

REQUIRED_MEASURED_BRANCHES = (
    "runnum",
    "vz_e",
    "e_p", "e_theta", "e_phi",
    "p_p", "p_theta", "p_phi",
)

ADDED_BRANCH_DTYPES: dict[str, np.dtype] = {
    "Egamma_external": np.dtype(np.float64),
    "Egamma_total": np.dtype(np.float64),
    "external_radiator_thickness": np.dtype(np.float64),
    "external_radiation_uniform": np.dtype(np.float64),
    "effective_beam_energy_externalISR": np.dtype(np.float64),
    "external_radiation_seed_ensemble": np.dtype(np.int32),
    "external_radiation_entry": np.dtype(np.int64),
}


@dataclass(frozen=True)
class FileIdentity:
    period: str
    target: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply deterministic external ISR to an RGC ISR ROOT tree."
    )
    parser.add_argument("input_root", type=Path)
    parser.add_argument("--geometry-json", required=True, type=Path)
    parser.add_argument("--tree-name", default=TREE_DEFAULT)
    parser.add_argument("--output-file", type=Path)
    parser.add_argument("--validation-dir", type=Path)
    parser.add_argument("--period", choices=("su22", "fa22", "sp23"))
    parser.add_argument("--target", choices=("NH3", "C", "CH2", "He", "ET"))
    parser.add_argument("--seed-ensemble", type=int, default=0)
    parser.add_argument("--step-size", default=STEP_SIZE_DEFAULT)
    parser.add_argument("--allow-unverified-geometry", action="store_true")
    parser.add_argument(
        "--external-reference-energy",
        choices=("remaining_after_internal", "nominal"),
        default="remaining_after_internal",
        help=(
            "Energy multiplying exp((r-1)/t). The sequential default uses "
            "the energy remaining after the existing internal-ISR draw."
        ),
    )
    parser.add_argument(
        "--keep-nonfinite",
        action="store_true",
        help="Write events with non-finite recalculated branches instead of failing.",
    )
    return parser.parse_args()


def beam_energy_from_run(runnum: np.ndarray) -> np.ndarray:
    run = np.asarray(runnum, dtype=np.int64)
    energy = np.full(run.shape, np.nan, dtype=np.float64)
    energy[(run >= 16042) & (run <= 17065)] = 10.5473
    energy[(run >= 17067) & (run <= 17716)] = 10.5563
    energy[(run >= 17717) & (run <= 17811)] = 10.5593
    return energy


def infer_period_from_filename(path: Path) -> str | None:
    name = path.name.lower()
    matches = [key for key in ("su22", "fa22", "sp23") if key in name]
    return matches[0] if len(matches) == 1 else None


def infer_target_from_filename(path: Path) -> str | None:
    name = path.name
    ordered = ("NH3", "CH2", "He", "ET", "C")
    matches = []
    for target in ordered:
        token = f"_{target}_"
        if token.lower() in name.lower():
            matches.append(target)
    return matches[0] if len(matches) == 1 else None


def infer_period_from_runs(runs: np.ndarray) -> str:
    run = np.asarray(runs, dtype=np.int64)
    labels = np.full(run.shape, "", dtype="U8")
    labels[(run >= 16042) & (run <= 16788)] = "su22"
    labels[(run >= 16843) & (run <= 17408)] = "fa22"
    labels[(run >= 17477) & (run <= 17811)] = "sp23"
    unique = sorted(set(labels.tolist()) - {""})
    if len(unique) != 1:
        raise RuntimeError(
            f"Could not identify one period from runnum; found {unique}."
        )
    if np.any(labels == ""):
        bad = np.unique(run[labels == ""])[:20]
        raise RuntimeError(f"Unrecognized RGC run numbers, first values: {bad}")
    return unique[0]


def resolve_identity(
    path: Path,
    run_sample: np.ndarray,
    period_override: str | None,
    target_override: str | None,
) -> FileIdentity:
    period_runs = infer_period_from_runs(run_sample)
    period_name = infer_period_from_filename(path)
    period = period_override or period_runs
    if period != period_runs:
        raise RuntimeError(
            f"Period override '{period}' disagrees with runnum period '{period_runs}'."
        )
    if period_name is not None and period_name != period:
        raise RuntimeError(
            f"Filename period '{period_name}' disagrees with runnum period '{period}'."
        )

    target_name = infer_target_from_filename(path)
    target = target_override or target_name
    if target is None:
        raise RuntimeError(
            "Could not infer target from filename. Supply --target explicitly."
        )
    if target_override and target_name and target_override != target_name:
        raise RuntimeError(
            f"Target override '{target_override}' disagrees with filename "
            f"target '{target_name}'."
        )
    return FileIdentity(period=period, target=target)


def default_output_path(input_path: Path) -> Path:
    stem = input_path.stem
    if "_ISR_" in stem:
        stem = stem.replace("_ISR_", "_ISR_externalISR_", 1)
    elif stem.endswith("_ISR"):
        stem += "_externalISR"
    else:
        stem += "_externalISR"
    return input_path.with_name(stem + input_path.suffix)


def file_fingerprint(path: Path) -> str:
    stat_info = path.stat()
    payload = (
        f"{path.resolve()}|{stat_info.st_size}|{stat_info.st_mtime_ns}"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_geometry(
    path: Path,
    identity: FileIdentity,
    allow_unverified: bool,
) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if int(payload.get("schema_version", -1)) != 1:
        raise RuntimeError("Geometry JSON must have schema_version=1.")

    try:
        record = payload["periods"][identity.period]["targets"][identity.target]
    except KeyError as exc:
        raise KeyError(
            f"No geometry for period={identity.period}, target={identity.target}"
        ) from exc

    if not bool(record.get("verified", False)) and not allow_unverified:
        raise RuntimeError(
            f"Geometry for {identity.period}/{identity.target} is marked "
            "verified=false. Correct the JSON or pass "
            "--allow-unverified-geometry for a diagnostic-only run."
        )
    return record


def component_thickness(
    vz: np.ndarray,
    component: dict[str, Any],
) -> np.ndarray:
    kind = str(component["kind"])
    x0_cm = float(component["radiation_length_cm"])
    if not np.isfinite(x0_cm) or x0_cm <= 0.0:
        raise RuntimeError(
            f"Invalid radiation_length_cm for {component.get('name')}: {x0_cm}"
        )

    if kind == "fixed":
        length = np.full(vz.shape, float(component["length_cm"]))
    elif kind == "vertex_dependent":
        z_up = float(component["upstream_face_z_cm"])
        z_down = float(component["downstream_face_z_cm"])
        if z_down <= z_up:
            raise RuntimeError(
                f"Invalid z interval for {component.get('name')}: "
                f"{z_up} to {z_down}"
            )
        direction = str(component.get("beam_direction", "+z"))
        if direction == "+z":
            length = np.clip(vz - z_up, 0.0, z_down - z_up)
        elif direction == "-z":
            length = np.clip(z_down - vz, 0.0, z_down - z_up)
        else:
            raise RuntimeError(
                f"Unsupported beam_direction '{direction}' for "
                f"{component.get('name')}"
            )
        length *= float(component.get("fill_fraction", 1.0))
    else:
        raise RuntimeError(
            f"Unsupported component kind '{kind}' for {component.get('name')}"
        )

    scale = float(component.get("thickness_scale", 1.0))
    return scale * length / x0_cm


def calculate_t_before(vz: np.ndarray, geometry: dict[str, Any]) -> np.ndarray:
    total = np.zeros(vz.shape, dtype=np.float64)
    components = geometry.get("components")
    if not isinstance(components, list) or not components:
        raise RuntimeError("Geometry record must contain a non-empty components list.")
    for component in components:
        total += component_thickness(vz, component)
    if np.any(~np.isfinite(total)) or np.any(total < 0.0):
        bad = total[(~np.isfinite(total)) | (total < 0.0)][:20]
        raise RuntimeError(f"Negative/nonfinite t_before values: {bad}")
    # A zero value is physically allowed for a vertex at or upstream of the
    # first modeled material boundary.  Such an event receives zero sampled
    # external-radiation loss below, corresponding to the t -> 0 limit.
    return total


def splitmix64(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=np.uint64)
    x = x + np.uint64(0x9E3779B97F4A7C15)
    z = x.copy()
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return z ^ (z >> np.uint64(31))


def deterministic_uniform(
    runnum: np.ndarray,
    absolute_entry: np.ndarray,
    fingerprint: str,
    seed_ensemble: int,
) -> np.ndarray:
    file_word = int(fingerprint[:16], 16)
    ensemble_word = (int(seed_ensemble) & 0xFFFFFFFF) << 32
    key = (
        np.asarray(absolute_entry, dtype=np.uint64)
        ^ (np.asarray(runnum, dtype=np.uint64) << np.uint64(32))
        ^ np.uint64(file_word)
        ^ np.uint64(ensemble_word)
    )
    hashed = splitmix64(key)
    # Open interval (0,1), using the upper 53 random bits.
    return ((hashed >> np.uint64(11)).astype(np.float64) + 0.5) / 2.0**53


def momentum_components(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    transverse = momentum * np.sin(theta)
    return (
        transverse * np.cos(phi),
        transverse * np.sin(phi),
        momentum * np.cos(theta),
    )


def calculate_mesonic_tmin(q2: np.ndarray, w: np.ndarray) -> np.ndarray:
    w2 = w * w
    threshold = PION_MASS_GEV + NEUTRON_MASS_GEV
    valid = (
        np.isfinite(q2) & np.isfinite(w) & (q2 > 0.0)
        & (w >= threshold) & (w2 > 0.0)
    )
    result = np.full(w.shape, np.nan, dtype=np.float64)
    if not np.any(valid):
        return result

    wv = w[valid]
    w2v = w2[valid]
    q2v = q2[valid]
    q0_star = (w2v - PROTON_MASS_GEV**2 - q2v) / (2.0 * wv)
    q_star = np.sqrt(np.maximum(q0_star**2 + q2v, 0.0))
    epi_star = (
        w2v + PION_MASS_GEV**2 - NEUTRON_MASS_GEV**2
    ) / (2.0 * wv)
    ppi_star = np.sqrt(np.maximum(epi_star**2 - PION_MASS_GEV**2, 0.0))
    result[valid] = (
        -q2v + PION_MASS_GEV**2
        - 2.0 * (q0_star * epi_star - q_star * ppi_star)
    )
    return result


def calculate_trento_phi(
    q_vec: np.ndarray,
    electron_vec: np.ndarray,
    pion_vec: np.ndarray,
) -> np.ndarray:
    q_norm = np.linalg.norm(q_vec, axis=1)
    qhat = np.divide(
        q_vec,
        q_norm[:, None],
        out=np.full_like(q_vec, np.nan),
        where=q_norm[:, None] > 0.0,
    )
    lepton_normal = np.cross(qhat, electron_vec)
    hadron_normal = np.cross(qhat, pion_vec)
    lepton_norm = np.linalg.norm(lepton_normal, axis=1)
    hadron_norm = np.linalg.norm(hadron_normal, axis=1)
    valid = (
        np.isfinite(qhat).all(axis=1)
        & (lepton_norm > 0.0) & (hadron_norm > 0.0)
    )
    result = np.full(q_norm.shape, np.nan, dtype=np.float64)
    if not np.any(valid):
        return result
    n_l = lepton_normal[valid] / lepton_norm[valid, None]
    n_h = hadron_normal[valid] / hadron_norm[valid, None]
    cosine = np.clip(np.einsum("ij,ij->i", n_l, n_h), -1.0, 1.0)
    unsigned = np.arccos(cosine)
    sine_sign = np.einsum(
        "ij,ij->i",
        np.cross(electron_vec[valid], pion_vec[valid]),
        qhat[valid],
    )
    result[valid] = np.where(sine_sign < 0.0, TWO_PI - unsigned, unsigned)
    return result


def lorentz_boost(
    energy: np.ndarray,
    momentum: np.ndarray,
    beta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    beta2 = np.einsum("ij,ij->i", beta, beta)
    valid = beta2 < 1.0
    gamma = np.full(beta2.shape, np.nan)
    gamma[valid] = 1.0 / np.sqrt(1.0 - beta2[valid])
    bp = np.einsum("ij,ij->i", beta, momentum)
    factor = np.zeros(beta2.shape)
    nonzero = valid & (beta2 > 0.0)
    factor[nonzero] = (
        ((gamma[nonzero] - 1.0) * bp[nonzero] / beta2[nonzero])
        - gamma[nonzero] * energy[nonzero]
    )
    boosted_p = momentum + factor[:, None] * beta
    boosted_e = gamma * (energy - bp)
    return boosted_e, boosted_p


def calculate_kinematics(
    beam_energy: np.ndarray,
    e_p: np.ndarray,
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    p_p: np.ndarray,
    p_theta: np.ndarray,
    p_phi: np.ndarray,
) -> dict[str, np.ndarray]:
    e_px, e_py, e_pz = momentum_components(e_p, e_theta, e_phi)
    pi_px, pi_py, pi_pz = momentum_components(p_p, p_theta, p_phi)

    electron_energy = np.sqrt(e_p**2 + ELECTRON_MASS_GEV**2)
    pion_energy = np.sqrt(p_p**2 + PION_MASS_GEV**2)
    beam_momentum = np.sqrt(
        np.maximum(beam_energy**2 - ELECTRON_MASS_GEV**2, 0.0)
    )

    electron_vec = np.column_stack((e_px, e_py, e_pz))
    pion_vec = np.column_stack((pi_px, pi_py, pi_pz))
    q_energy = beam_energy - electron_energy
    q_vec = np.column_stack((-e_px, -e_py, beam_momentum - e_pz))
    q_three_sq = np.einsum("ij,ij->i", q_vec, q_vec)
    q2 = q_three_sq - q_energy**2

    w2 = (PROTON_MASS_GEV + q_energy)**2 - q_three_sq
    w = np.sqrt(np.maximum(w2, 0.0))

    missing_energy = PROTON_MASS_GEV + q_energy - pion_energy
    missing_vec = q_vec - pion_vec
    mx2 = missing_energy**2 - np.einsum("ij,ij->i", missing_vec, missing_vec)

    with np.errstate(divide="ignore", invalid="ignore"):
        x = q2 / (2.0 * PROTON_MASS_GEV * q_energy)
        y = q_energy / beam_energy
        z = pion_energy / q_energy

    delta_energy = q_energy - pion_energy
    delta_vec = q_vec - pion_vec
    t = delta_energy**2 - np.einsum("ij,ij->i", delta_vec, delta_vec)
    tmin = calculate_mesonic_tmin(q2, w)
    tprime = t - tmin

    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = 2.0 * PROTON_MASS_GEV * x / np.sqrt(q2)
        gamma2 = gamma**2
        common = 1.0 - y - y**2 * gamma2 / 4.0
        dep_a = (
            1.0 - y + y**2 / 2.0 + y**2 * gamma2 / 4.0
        ) / (1.0 + gamma2)
        dep_b = common / (1.0 + gamma2)
        dep_c = y * (1.0 - y / 2.0) / np.sqrt(1.0 + gamma2)
        common_sqrt = np.sqrt(np.maximum(common, 0.0))
        dep_v = (2.0 - y) * common_sqrt / (1.0 + gamma2)
        dep_w = y * common_sqrt / np.sqrt(1.0 + gamma2)

    # gamma*-proton center-of-mass boost.
    gn_energy = q_energy + PROTON_MASS_GEV
    beta = np.divide(
        q_vec,
        gn_energy[:, None],
        out=np.full_like(q_vec, np.nan),
        where=gn_energy[:, None] != 0.0,
    )
    pion_cm_e, pion_cm_p = lorentz_boost(pion_energy, pion_vec, beta)
    q_cm_e, q_cm_p = lorentz_boost(q_energy, q_vec, beta)
    q_cm_mag = np.linalg.norm(q_cm_p, axis=1)
    qhat_cm = np.divide(
        q_cm_p,
        q_cm_mag[:, None],
        out=np.full_like(q_cm_p, np.nan),
        where=q_cm_mag[:, None] > 0.0,
    )
    longitudinal = np.einsum("ij,ij->i", pion_cm_p, qhat_cm)
    transverse_vec = pion_cm_p - longitudinal[:, None] * qhat_cm
    pt = np.linalg.norm(transverse_vec, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        xf = 2.0 * longitudinal / w

    phi = calculate_trento_phi(q_vec, electron_vec, pion_vec)

    return {
        "Q2": q2, "W": w, "Mx2": mx2, "x": x, "y": y,
        "t": t, "tmin": tmin, "tprime": tprime,
        "z": z, "xF": xf, "pT": pt, "phi": phi,
        "DepA": dep_a, "DepB": dep_b, "DepC": dep_c,
        "DepV": dep_v, "DepW": dep_w,
    }


def find_tree(path: Path, requested: str) -> str:
    with uproot.open(path) as root_file:
        direct = {str(key).split(";")[0]: str(key) for key in root_file.keys()}
        if requested in direct:
            return direct[requested]
        candidates = []
        for key in root_file.keys():
            obj = root_file[key]
            if hasattr(obj, "arrays") and hasattr(obj, "num_entries"):
                candidates.append(str(key))
    if len(candidates) == 1:
        return candidates[0]
    raise KeyError(
        f"Tree '{requested}' not found; tree-like objects are {candidates}."
    )


def resolve_internal_egamma_branch(branches: set[str]) -> str:
    if "Egamma_internal" in branches:
        return "Egamma_internal"
    if "Egamma" in branches:
        return "Egamma"
    raise KeyError(
        "No internal ISR photon-energy branch found. Expected "
        "'Egamma_internal' or 'Egamma'. The supplied momentum-correction "
        "workflow may have dropped Egamma; regenerate/preserve it before "
        "running this transformer."
    )


def histogram_accumulator() -> dict[str, list[np.ndarray]]:
    names = (
        "vz_e", "t_before", "egamma_internal", "egamma_external",
        "egamma_total", "effective_beam", "mx2_before", "mx2_after",
        "delta_mx2", "delta_x", "delta_tprime",
    )
    return {name: [] for name in names}


def append_sample(
    accumulator: dict[str, list[np.ndarray]],
    values: dict[str, np.ndarray],
    maximum_total: int = 2_000_000,
) -> None:
    existing = sum(len(part) for part in accumulator["vz_e"])
    remaining = maximum_total - existing
    if remaining <= 0:
        return
    count = min(remaining, len(next(iter(values.values()))))
    if count <= 0:
        return
    # Deterministic thinning by regularly spaced indices.
    indices = np.linspace(0, len(next(iter(values.values()))) - 1, count).astype(int)
    for name, array in values.items():
        accumulator[name].append(np.asarray(array)[indices])


def concatenate_samples(acc: dict[str, list[np.ndarray]]) -> dict[str, np.ndarray]:
    return {
        name: np.concatenate(parts) if parts else np.array([], dtype=float)
        for name, parts in acc.items()
    }


def save_histogram(
    values: np.ndarray,
    path: Path,
    xlabel: str,
    bins: int = 120,
    log_y: bool = False,
) -> None:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(finite, bins=bins, histtype="step", linewidth=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Events")
    if log_y:
        ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_2d(
    x: np.ndarray,
    y: np.ndarray,
    path: Path,
    xlabel: str,
    ylabel: str,
    bins: int = 120,
) -> None:
    mask = np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    hist = ax.hist2d(x[mask], y[mask], bins=bins)
    fig.colorbar(hist[3], ax=ax, label="Events")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_validation_plots(samples: dict[str, np.ndarray], output_dir: Path) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    files: list[str] = []

    one_d = (
        ("t_before", "external_radiator_thickness.png", r"$t_{\mathrm{before}}$ (radiation lengths)", False),
        ("egamma_external", "external_photon_energy.png", r"$E_{\gamma}^{\mathrm{external}}$ (GeV)", True),
        ("egamma_internal", "internal_photon_energy.png", r"$E_{\gamma}^{\mathrm{internal}}$ (GeV)", True),
        ("egamma_total", "total_photon_energy.png", r"$E_{\gamma}^{\mathrm{total}}$ (GeV)", True),
        ("effective_beam", "effective_beam_energy.png", r"$E_{\mathrm{beam}}^{\mathrm{effective}}$ (GeV)", False),
        ("delta_mx2", "delta_mx2.png", r"$M_X^2(\mathrm{int+ext})-M_X^2(\mathrm{int})$ (GeV$^2$)", False),
        ("delta_x", "delta_x.png", r"$x_B(\mathrm{int+ext})-x_B(\mathrm{int})$", False),
        ("delta_tprime", "delta_tprime.png", r"$t'(\mathrm{int+ext})-t'(\mathrm{int})$ (GeV$^2$)", False),
    )
    for key, filename, xlabel, log_y in one_d:
        save_histogram(samples[key], output_dir / filename, xlabel, log_y=log_y)
        files.append(filename)

    two_d = (
        ("vz_e", "t_before", "thickness_vs_vz.png", r"$v_z^e$ (cm)", r"$t_{\mathrm{before}}$ (radiation lengths)"),
        ("vz_e", "egamma_external", "external_energy_vs_vz.png", r"$v_z^e$ (cm)", r"$E_{\gamma}^{\mathrm{external}}$ (GeV)"),
        ("mx2_before", "mx2_after", "mx2_before_vs_after.png", r"$M_X^2$ internal ISR (GeV$^2$)", r"$M_X^2$ internal + external ISR (GeV$^2$)"),
    )
    for xkey, ykey, filename, xlabel, ylabel in two_d:
        save_2d(samples[xkey], samples[ykey], output_dir / filename, xlabel, ylabel)
        files.append(filename)

    # Overlay missing-mass distributions.
    before = samples["mx2_before"]
    after = samples["mx2_after"]
    finite = np.concatenate((before[np.isfinite(before)], after[np.isfinite(after)]))
    if finite.size:
        low, high = np.quantile(finite, [0.001, 0.999])
        bins = np.linspace(low, high, 161)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(before[np.isfinite(before)], bins=bins, histtype="step", linewidth=1.5, label="Internal ISR")
        ax.hist(after[np.isfinite(after)], bins=bins, histtype="step", linewidth=1.5, label="Internal + external ISR")
        ax.set_xlabel(r"$M_X^2$ (GeV$^2$)")
        ax.set_ylabel("Events")
        ax.legend()
        fig.tight_layout()
        filename = "mx2_overlay.png"
        fig.savefig(output_dir / filename, dpi=180)
        plt.close(fig)
        files.append(filename)

    return files


def quantile_summary(values: np.ndarray) -> dict[str, float]:
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return {}
    qs = np.quantile(finite, [0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0])
    return {
        "count": int(finite.size),
        "mean": float(np.mean(finite)),
        "minimum": float(qs[0]),
        "q01": float(qs[1]),
        "q05": float(qs[2]),
        "median": float(qs[3]),
        "q95": float(qs[4]),
        "q99": float(qs[5]),
        "maximum": float(qs[6]),
    }


def main() -> int:
    args = parse_args()
    input_path = args.input_root.resolve()
    geometry_path = args.geometry_json.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if not geometry_path.is_file():
        raise FileNotFoundError(geometry_path)

    tree_name = find_tree(input_path, args.tree_name)
    with uproot.open(input_path) as root_file:
        tree = root_file[tree_name]
        branches = set(tree.keys())
        missing = sorted(set(REQUIRED_MEASURED_BRANCHES) - branches)
        if missing:
            raise KeyError("Input tree is missing: " + ", ".join(missing))
        egamma_branch = resolve_internal_egamma_branch(branches)
        run_sample = tree["runnum"].array(
            entry_stop=min(tree.num_entries, 1_000_000), library="np"
        )
        input_entries = int(tree.num_entries)
        input_branch_dtypes: dict[str, np.dtype] = {}
        unsupported_branches: list[str] = []
        for branch_name in tree.keys():
            branch = tree[branch_name]
            try:
                input_branch_dtypes[branch_name] = np.dtype(
                    branch.interpretation.numpy_dtype
                )
            except Exception:
                unsupported_branches.append(branch_name)
        if unsupported_branches:
            raise TypeError(
                "The input tree contains branches whose NumPy dtype could not "
                "be determined: " + ", ".join(sorted(unsupported_branches))
            )

    identity = resolve_identity(
        input_path, run_sample, args.period, args.target
    )
    geometry = load_geometry(
        geometry_path, identity, args.allow_unverified_geometry
    )

    output_path = (
        args.output_file.resolve()
        if args.output_file
        else default_output_path(input_path).resolve()
    )
    validation_dir = (
        args.validation_dir.resolve()
        if args.validation_dir
        else output_path.with_suffix("").with_name(
            output_path.stem + "_validation"
        )
    )
    provenance_path = output_path.with_suffix(".provenance.json")

    for path in (output_path, provenance_path):
        if path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing output: {path}"
            )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fingerprint = file_fingerprint(input_path)
    read_branches = sorted(
        set(REQUIRED_MEASURED_BRANCHES)
        | {egamma_branch}
        | {name for name in RECALCULATED_BRANCHES if name in branches}
    )

    output_schema: dict[str, Any] = {}
    for name, dtype in input_branch_dtypes.items():
        if name in RECALCULATED_BRANCHES:
            output_schema[name] = np.dtype(np.float64)
        else:
            output_schema[name] = dtype
    for name, dtype in ADDED_BRANCH_DTYPES.items():
        if name in output_schema:
            raise RuntimeError(
                f"Input already contains provenance branch '{name}'."
            )
        output_schema[name] = dtype

    accumulator = histogram_accumulator()
    absolute_start = 0
    written_entries = 0
    invalid_events = 0
    endpoint_capped_events = 0

    with uproot.recreate(output_path, compression=uproot.ZLIB(9)) as output_file:
        output_tree = output_file.mktree(tree_name, output_schema)

        for arrays in uproot.iterate(
            f"{input_path}:{tree_name}",
            expressions=None,
            step_size=args.step_size,
            library="ak",
        ):
            count = len(arrays)
            raw = {name: ak.to_numpy(arrays[name]) for name in arrays.fields}

            runnum = np.asarray(raw["runnum"], dtype=np.int64)
            vz = np.asarray(raw["vz_e"], dtype=np.float64)
            nominal_beam = beam_energy_from_run(runnum)
            internal_egamma = np.asarray(raw[egamma_branch], dtype=np.float64)

            if np.any(~np.isfinite(nominal_beam)):
                bad = np.unique(runnum[~np.isfinite(nominal_beam)])[:20]
                raise RuntimeError(f"No beam energy for run numbers: {bad}")
            if np.any(~np.isfinite(internal_egamma)) or np.any(internal_egamma < 0):
                raise RuntimeError("Internal Egamma contains negative/nonfinite values.")

            t_before = calculate_t_before(vz, geometry)
            entries = np.arange(
                absolute_start, absolute_start + count, dtype=np.int64
            )
            uniform = deterministic_uniform(
                runnum, entries, fingerprint, args.seed_ensemble
            )
            remaining_after_internal = nominal_beam - internal_egamma
            if np.any(remaining_after_internal <= ELECTRON_MASS_GEV):
                raise RuntimeError(
                    "Existing internal Egamma leaves nonpositive beam energy."
                )

            reference = (
                remaining_after_internal
                if args.external_reference_energy == "remaining_after_internal"
                else nominal_beam
            )
            external_egamma = np.zeros_like(reference, dtype=np.float64)
            positive_thickness = t_before > 0.0
            external_egamma[positive_thickness] = reference[positive_thickness] * np.exp(
                (uniform[positive_thickness] - 1.0)
                / t_before[positive_thickness]
            )

            # Mathematically, for the sequential default, exp((u-1)/t) is
            # strictly smaller than one because u is sampled on the open
            # interval (0,1).  For extremely small (1-u)/t, however, exp can
            # round to exactly 1.0 in float64.  That would consume all of the
            # remaining beam energy and produce E_eff <= m_e for a tiny number
            # of events.  Apply only the physical/numerical endpoint cap here:
            # leave the effective incident electron energy strictly above m_e.
            max_external_egamma = np.nextafter(
                remaining_after_internal - ELECTRON_MASS_GEV,
                -np.inf,
            )
            endpoint_capped = external_egamma > max_external_egamma
            external_egamma = np.minimum(external_egamma, max_external_egamma)

            endpoint_capped_events += int(np.count_nonzero(endpoint_capped))

            total_egamma = internal_egamma + external_egamma
            effective_beam = nominal_beam - total_egamma

            valid_beam = effective_beam > ELECTRON_MASS_GEV
            if not np.all(valid_beam):
                first = entries[~valid_beam][:20]
                raise RuntimeError(
                    f"{np.count_nonzero(~valid_beam)} events still have "
                    "nonpositive effective beam energy after the endpoint cap; "
                    f"first absolute entries: {first}."
                )

            before = {
                name: np.asarray(raw[name], dtype=np.float64)
                for name in RECALCULATED_BRANCHES
                if name in raw
            }
            after = calculate_kinematics(
                effective_beam,
                np.asarray(raw["e_p"], dtype=np.float64),
                np.asarray(raw["e_theta"], dtype=np.float64),
                np.asarray(raw["e_phi"], dtype=np.float64),
                np.asarray(raw["p_p"], dtype=np.float64),
                np.asarray(raw["p_theta"], dtype=np.float64),
                np.asarray(raw["p_phi"], dtype=np.float64),
            )

            finite_mask = np.ones(count, dtype=bool)
            for values in after.values():
                finite_mask &= np.isfinite(values)
            invalid = int(np.count_nonzero(~finite_mask))
            invalid_events += invalid

            # A sufficiently hard added external-radiation draw can make the
            # already measured final state kinematically impossible at the new
            # effective beam energy (for example E_e' >= E_beam, Q2 <= 0, or
            # an undefined gamma*-p center-of-mass boost).  Such events cannot
            # belong to the reradiated sample.  Remove them rather than writing
            # NaNs or aborting the entire production.  This is the event-loss
            # component of the additional-radiation response.
            write_mask = finite_mask
            if args.keep_nonfinite:
                write_mask = np.ones(count, dtype=bool)

            output: dict[str, Any] = {
                name: np.asarray(values)[write_mask]
                for name, values in raw.items()
            }
            for name, values in after.items():
                if name in branches:
                    output[name] = np.asarray(values, dtype=np.float64)[write_mask]
            output.update({
                "Egamma_external": external_egamma[write_mask],
                "Egamma_total": total_egamma[write_mask],
                "external_radiator_thickness": t_before[write_mask],
                "external_radiation_uniform": uniform[write_mask],
                "effective_beam_energy_externalISR": effective_beam[write_mask],
                "external_radiation_seed_ensemble": np.full(
                    int(np.count_nonzero(write_mask)), args.seed_ensemble, dtype=np.int32
                ),
                "external_radiation_entry": entries[write_mask],
            })

            output_tree.extend(output)
            written_entries += int(np.count_nonzero(write_mask))

            append_sample(accumulator, {
                "vz_e": vz[finite_mask],
                "t_before": t_before[finite_mask],
                "egamma_internal": internal_egamma[finite_mask],
                "egamma_external": external_egamma[finite_mask],
                "egamma_total": total_egamma[finite_mask],
                "effective_beam": effective_beam[finite_mask],
                "mx2_before": before.get("Mx2", np.full(count, np.nan))[finite_mask],
                "mx2_after": after["Mx2"][finite_mask],
                "delta_mx2": (after["Mx2"] - before.get("Mx2", np.full(count, np.nan)))[finite_mask],
                "delta_x": (after["x"] - before.get("x", np.full(count, np.nan)))[finite_mask],
                "delta_tprime": (after["tprime"] - before.get("tprime", np.full(count, np.nan)))[finite_mask],
            })
            absolute_start += count

    samples = concatenate_samples(accumulator)
    validation_files = make_validation_plots(samples, validation_dir)

    provenance = {
        "schema_version": 1,
        "model_version": MODEL_VERSION,
        "input_file": str(input_path),
        "input_file_fingerprint": fingerprint,
        "input_tree": tree_name,
        "output_file": str(output_path),
        "output_tree": tree_name,
        "geometry_json": str(geometry_path),
        "period": identity.period,
        "target": identity.target,
        "geometry_verified": bool(geometry.get("verified", False)),
        "internal_egamma_branch": egamma_branch,
        "external_reference_energy": args.external_reference_energy,
        "sampling_formula": (
            "Egamma_external = E_reference * exp((uniform - 1) / t_before)"
        ),
        "seed_ensemble": args.seed_ensemble,
        "deterministic_key": (
            "runnum, absolute tree-entry index, input-file fingerprint, "
            "seed ensemble"
        ),
        "input_entries": input_entries,
        "output_entries": written_entries,
        "removed_nonphysical_recalculated_events": invalid_events,
        "external_energy_endpoint_capped_events": endpoint_capped_events,
        "recalculated_branches": [
            name for name in RECALCULATED_BRANCHES if name in branches
        ],
        "added_branches": list(ADDED_BRANCH_DTYPES),
        "validation_directory": str(validation_dir),
        "validation_files": validation_files,
        "summaries": {
            "vz_e_cm": quantile_summary(samples["vz_e"]),
            "t_before_radiation_lengths": quantile_summary(samples["t_before"]),
            "Egamma_internal_GeV": quantile_summary(samples["egamma_internal"]),
            "Egamma_external_GeV": quantile_summary(samples["egamma_external"]),
            "Egamma_total_GeV": quantile_summary(samples["egamma_total"]),
            "effective_beam_energy_GeV": quantile_summary(samples["effective_beam"]),
            "delta_Mx2_GeV2": quantile_summary(samples["delta_mx2"]),
            "delta_x": quantile_summary(samples["delta_x"]),
            "delta_tprime_GeV2": quantile_summary(samples["delta_tprime"]),
        },
        "geometry_record": geometry,
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"Input        : {input_path}")
    print(f"Output       : {output_path}")
    print(f"Provenance   : {provenance_path}")
    print(f"Validation   : {validation_dir}")
    print(f"Period/target: {identity.period}/{identity.target}")
    print(f"Input entries : {input_entries}")
    print(f"Output entries: {written_entries}")
    print(f"Removed       : {invalid_events}")
    print(f"Endpoint caps : {endpoint_capped_events}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"\nFATAL ERROR:\n{exc}", file=sys.stderr)
        raise SystemExit(1)
