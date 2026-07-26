#!/usr/bin/env python3
"""
extract_structure_function_ratios_v5.py

Initial standalone event-level asymmetry extraction for the RGC exclusive
e p -> e' n pi+ analysis.

Run this program from

    RGC_enpi+/asymmetry_extraction/

with sibling analysis directories

    ../channel_selection/
    ../dilution_factor/
    ../momentum_corrections/

Nominal statistical model
-------------------------
The three RGC run periods are fitted simultaneously in every (xB, -tprime)
kinematic bin.  The physics parameters are common to Su22, Fa22, and Sp23,
while the likelihood uses the period-dependent beam polarization, the signed
run-by-run target polarization, the run/helicity accumulated charges, and the
period/bin-dependent dilution factors.

The nominal estimator is an unbinned conditional likelihood.  For an event
with accepted kinematics x in period p, the probability of its observed
(run, beam-helicity) state is

    P(r,h | x,p) =
        Q(r,h) g(x; r,h,p) /
        sum_{r' in p, h'=+/-1} Q(r',h') g(x; r',h',p).

The sum in the denominator is restricted to the same run period.  This permits
the three periods to constrain one common set of structure-function ratios
without assuming that the detector acceptance is identical between periods.
Within one period, helicity- and target-state-independent acceptance is assumed.

The fitted structure-function ratios are

    u1   = F_UU^{cos(phi)}     / F_UU
    u2   = F_UU^{cos(2phi)}    / F_UU
    lu1  = F_LU^{sin(phi)}     / F_UU
    ul1  = F_UL^{sin(phi)}     / F_UU
    ul2  = F_UL^{sin(2phi)}    / F_UU
    ll0  = F_LL                / F_UU
    ll1  = F_LL^{cos(phi)}     / F_UU

The unpolarized cosine modulations u1 and u2 float in the nominal fit.

Target-axis treatments
----------------------
Four fits are performed in every kinematic bin.

  nominal
    P_L = P_t cos(theta_gamma), P_T = 0.

  no_projection
    P_L = P_t, P_T = 0.

  longitudinal_scaled_transverse_plus
    P_L = P_t cos(theta_gamma), P_T = P_t sin(theta_gamma).  The transverse
    stress amplitudes are fixed from the signed nominal longitudinal amplitudes
    in the same bin using the documented harmonic mapping.

  longitudinal_scaled_transverse_minus
    Same geometry and magnitudes as longitudinal_scaled_transverse_plus, but
    every mapped transverse amplitude has the opposite sign.

These longitudinal-scale variants are loose, data-driven stress tests.  They
are not physical models asserting equality between longitudinal and transverse
amplitudes.  For each fitted observable, the target-axis systematic is the
largest absolute displacement from nominal among no_projection and the two
opposite-sign longitudinal-scaled transverse fits.

The transverse stress terms use the standard one-hadron harmonic/depolarization
mapping evaluated for a beam-axis target, phi_S = 0 with the signed target
polarization carrying the target orientation:

  UT:
    A_UT^{sin(phi-phi_S)}                 coefficient 1
    A_UT^{sin(phi+phi_S)}                 coefficient DepB/DepA
    A_UT^{sin(3phi-phi_S)}                coefficient DepB/DepA
    A_UT^{sin(phi_S)}                     coefficient DepV/DepA
    A_UT^{sin(2phi-phi_S)}                coefficient DepV/DepA

  LT:
    A_LT^{cos(phi-phi_S)}                 coefficient DepC/DepA
    A_LT^{cos(phi_S)}                     coefficient DepW/DepA
    A_LT^{cos(2phi-phi_S)}                coefficient DepW/DepA

At phi_S = 0, sin(phi_S) vanishes.  Several remaining transverse terms are
degenerate with longitudinal harmonics; they are therefore fixed rather than
fitted and are used only to estimate the leakage envelope.

Dilution-factor uncertainty
---------------------------
The recommended nominal dilution factor, the average of Methods 1 and 2, is
read from the production dilution-factor JSON.  Its bootstrap statistical
uncertainty is included through one Gaussian-constrained dilution nuisance
parameter per run period and kinematic bin.  The correlated dilution-model
scale uncertainty is intentionally not included here and will be imposed later.

Beam- and target-polarization uncertainties are intentionally ignored by this
program.  Their statistical and systematic components will be combined later
into correlated polarization scale systematics.

Parallelism
-----------
ROOT input is read once in a serial preprocessing pass and written to a compact
cache.  Independent kinematic-bin fits are then distributed across at most
eight worker processes.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

try:
    from iminuit import Minuit
except ImportError as exc:
    raise RuntimeError(
        "This program requires iminuit. Install it with: python -m pip install iminuit"
    ) from exc
# endtry


# =============================================================================
# Fixed analysis definitions
# =============================================================================

PERIODS: tuple[str, ...] = ("su22", "fa22", "sp23")
PERIOD_LABELS: dict[str, str] = {
    "su22": "Su22",
    "fa22": "Fa22",
    "sp23": "Sp23",
}
PERIOD_INDEX: dict[str, int] = {
    period: index for index, period in enumerate(PERIODS)
}

BEAM_POLARIZATION: dict[str, float] = {
    "su22": 0.8384,
    "fa22": 0.8372,
    "sp23": 0.8040,
}

BEAM_ENERGY_GEV: dict[str, float] = {
    "su22": 10.5473,
    "fa22": 10.5563,
    "sp23": 10.5593,
}

XB_BINS: tuple[tuple[float, float], ...] = (
    (0.10, 0.25),
    (0.25, 0.35),
    (0.35, 0.45),
    (0.45, 0.60),
)

MINUS_TPRIME_BINS_GEV2: tuple[tuple[float, float], ...] = (
    (0.05, 0.25),
    (0.25, 0.45),
    (0.45, 0.65),
    (0.65, 0.85),
    (0.85, 1.05),
    (1.05, 1.25),
)

NUMBER_OF_BINS = len(XB_BINS) * len(MINUS_TPRIME_BINS_GEV2)
MAXIMUM_WORKERS = 8
DEFAULT_TREE_NAME = "PhysicsEvents"
DEFAULT_CHUNK_SIZE = "250 MB"
DEFAULT_OUTPUT_DIR = Path("output/asymmetry_extraction")
DEFAULT_CACHE_PATH = DEFAULT_OUTPUT_DIR / "cache/selected_events_v5.npz"

DEFAULT_RUN_INFO_CSV = Path("clas12_run_info.csv")
DEFAULT_CUT_JSON = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "final_carbon_assisted_cuts/tables/"
    "final_carbon_assisted_mx2_cuts.json"
)

DEFAULT_DILUTION_DIR = Path(
    "../dilution_factor/output/dilution_factor_determination"
)

PAPER_VERSIONS_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
    "paper_versions"
)

DEFAULT_INPUTS: dict[str, Path] = {
    period: PAPER_VERSIONS_DIR
    / f"rgc_{period}_inb_NH3_epi+_mom_corrections.root"
    for period in PERIODS
}

FIT_VARIANTS: tuple[str, ...] = (
    "nominal",
    "no_projection",
    "longitudinal_scaled_transverse_plus",
    "longitudinal_scaled_transverse_minus",
)


PHYSICS_PARAMETERS: tuple[str, ...] = (
    "u1",
    "u2",
    "lu1",
    "ul1",
    "ul2",
    "ll0",
    "ll1",
)

PARAMETER_INITIAL_VALUES: dict[str, float] = {
    "u1": 0.0,
    "u2": 0.0,
    "lu1": 0.05,
    "ul1": 0.0,
    "ul2": 0.0,
    "ll0": 0.0,
    "ll1": 0.0,
}

PARAMETER_LIMITS: dict[str, tuple[float, float]] = {
    name: (-1.5, 1.5) for name in PHYSICS_PARAMETERS
}


PARAMETER_LABELS: dict[str, str] = {
    "u1": r"$F_{UU}^{\cos\phi}/F_{UU}$",
    "u2": r"$F_{UU}^{\cos2\phi}/F_{UU}$",
    "lu1": r"$F_{LU}^{\sin\phi}/F_{UU}$",
    "ul1": r"$F_{UL}^{\sin\phi}/F_{UU}$",
    "ul2": r"$F_{UL}^{\sin2\phi}/F_{UU}$",
    "ll0": r"$F_{LL}/F_{UU}$",
    "ll1": r"$F_{LL}^{\cos\phi}/F_{UU}$",
}

PARAMETER_Y_LIMITS: dict[str, tuple[float, float] | None] = {
    "u1": None,
    "u2": None,
    "lu1": (-0.4, 0.4),
    "ul1": (-0.4, 0.4),
    "ul2": (-0.4, 0.4),
    "ll0": (-1.0, 1.0),
    "ll1": (-1.0, 1.0),
}

PROBABILITY_FLOOR = 1.0e-300
CROSS_SECTION_FLOOR = 1.0e-10
INVALID_NLL = 1.0e100

BRANCH_ALIASES: dict[str, tuple[str, ...]] = {
    "runnum": ("runnum", "run", "run_number", "RunNumber"),
    "helicity": ("helicity", "hel", "beam_helicity"),
    "xB": ("xB", "x", "xb", "x_b"),
    "tprime": ("tprime", "t_prime", "tp", "tPrime"),
    "Mx2": (
        "Mx2",
        "mx2",
        "Mx2_epi",
        "Mx2_epip",
        "missing_mass_squared",
        "missing_mass2",
    ),
    "phi": ("phi", "phi1", "phi_h", "trento_phi"),
    "Q2": ("Q2", "q2"),
    "e_p": ("e_p", "p_e", "electron_p"),
    "e_theta": ("e_theta", "theta_e", "electron_theta"),
    "DepA": ("DepA", "depA"),
    "DepB": ("DepB", "depB"),
    "DepC": ("DepC", "depC"),
    "DepV": ("DepV", "depV"),
    "DepW": ("DepW", "depW"),
}


# =============================================================================
# Data containers
# =============================================================================

@dataclass(frozen=True)
class RunRecord:
    period: str
    run: int
    charge_plus: float
    charge_minus: float
    target_polarization: float
    target_polarization_uncertainty: float


@dataclass(frozen=True)
class DilutionRecord:
    period: str
    bin_number: int
    x_index: int
    t_index: int
    value: float
    stat_uncertainty: float


@dataclass(frozen=True)
class CutRecord:
    period: str
    bin_number: int
    x_index: int
    t_index: int
    low_gev2: float
    high_gev2: float


# =============================================================================
# General utilities
# =============================================================================

def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    # endif
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    # endif
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    # endif
    if isinstance(value, np.integer):
        return int(value)
    # endif
    if isinstance(value, np.floating):
        value = float(value)
    # endif
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    # endif
    if isinstance(value, Path):
        return str(value)
    # endif
    return value


def write_json(path: Path, payload: Any) -> None:
    ensure_directory(path.parent)
    path.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
        # endfor
    # endwith
    return digest.hexdigest()


def combined_bin_number(x_index: int, t_index: int) -> int:
    return x_index * len(MINUS_TPRIME_BINS_GEV2) + t_index + 1


def bin_indices(
    x_values: np.ndarray,
    minus_tprime_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_index = np.full(x_values.shape, -1, dtype=np.int16)
    t_index = np.full(minus_tprime_values.shape, -1, dtype=np.int16)

    for index, (low, high) in enumerate(XB_BINS):
        mask = (x_values >= low) & (x_values < high)
        x_index[mask] = index
    # endfor

    for index, (low, high) in enumerate(MINUS_TPRIME_BINS_GEV2):
        mask = (minus_tprime_values >= low) & (minus_tprime_values < high)
        t_index[mask] = index
    # endfor

    valid = (x_index >= 0) & (t_index >= 0)
    number = np.full(x_values.shape, -1, dtype=np.int16)
    number[valid] = (
        x_index[valid] * len(MINUS_TPRIME_BINS_GEV2)
        + t_index[valid]
        + 1
    )
    return x_index, t_index, number


def resolve_branch(
    tree: uproot.behaviors.TTree.TTree,
    aliases: Iterable[str],
) -> str:
    available = {str(key).split(";")[0] for key in tree.keys()}
    for alias in aliases:
        if alias in available:
            return alias
        # endif
    # endfor
    raise RuntimeError(
        f"None of the aliases {tuple(aliases)} was found. "
        f"Available branches include {sorted(available)[:120]}"
    )


def detect_angle_unit(values: np.ndarray, branch_name: str) -> str:
    finite = np.abs(values[np.isfinite(values)])
    if finite.size == 0:
        raise RuntimeError(f"Cannot determine units for empty branch {branch_name}.")
    # endif
    percentile = float(np.percentile(finite, 99.5))
    return "degrees" if percentile > 2.0 * math.pi + 0.25 else "radians"


def angle_to_radians(values: np.ndarray, branch_name: str) -> np.ndarray:
    unit = detect_angle_unit(values, branch_name)
    if unit == "degrees":
        return np.deg2rad(values)
    # endif
    return values


# =============================================================================
# Run information
# =============================================================================

def parse_run_info_csv(path: Path) -> dict[int, RunRecord]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing run-information CSV: {path}")
    # endif

    period_lookup = {
        "Su22": "su22",
        "Fa22": "fa22",
        "Sp23": "sp23",
    }

    records: dict[int, RunRecord] = {}
    current_period: str | None = None
    current_target: str | None = None

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip()
        if not line:
            continue
        # endif

        if line.startswith("#"):
            header = line[1:].strip().split()
            current_period = None
            current_target = None
            if len(header) >= 3 and header[0] == "RGC":
                current_period = period_lookup.get(header[1])
                current_target = header[2]
            # endif
            continue
        # endif

        if current_period not in PERIODS or current_target != "NH3":
            continue
        # endif

        fields = [field.strip() for field in line.split(",")]
        if len(fields) < 6:
            raise RuntimeError(
                f"Malformed NH3 run-info row at line {line_number}: {raw_line}"
            )
        # endif

        try:
            run = int(fields[0])
            charge_plus = float(fields[2])
            charge_minus = float(fields[3])
            target_polarization = float(fields[4])
            target_uncertainty = float(fields[-1])
        except ValueError as exc:
            raise RuntimeError(
                f"Invalid numeric field at line {line_number}: {raw_line}"
            ) from exc
        # endtry

        if run in records:
            raise RuntimeError(f"Duplicate NH3 run {run} in {path}.")
        # endif
        if not (
            math.isfinite(charge_plus)
            and math.isfinite(charge_minus)
            and charge_plus >= 0.0
            and charge_minus >= 0.0
        ):
            raise RuntimeError(
                f"Run {run} has invalid helicity charges: "
                f"Q+={charge_plus}, Q-={charge_minus}. Charges must be "
                "finite and nonnegative."
            )
        # endif

        disabled_by_qa = charge_plus == 0.0 and charge_minus == 0.0
        only_one_helicity_disabled = (charge_plus == 0.0) != (charge_minus == 0.0)
        if only_one_helicity_disabled:
            raise RuntimeError(
                f"Run {run} has only one zero helicity charge: "
                f"Q+={charge_plus}, Q-={charge_minus}. The current likelihood "
                "expects either two positive helicity charges or a fully "
                "QA-disabled run with both charges set to zero."
            )
        # endif

        # Fully QA-disabled runs are intentionally retained in the CSV with
        # Q+ = Q- = 0.  They are valid bookkeeping rows but must never enter
        # the likelihood.  A nonzero finite target polarization is required
        # only for active runs.
        if not disabled_by_qa and (
            not math.isfinite(target_polarization)
            or target_polarization == 0.0
        ):
            raise RuntimeError(
                f"Active run {run} has invalid target polarization "
                f"{target_polarization}."
            )
        # endif

        records[run] = RunRecord(
            period=current_period,
            run=run,
            charge_plus=charge_plus,
            charge_minus=charge_minus,
            target_polarization=target_polarization,
            target_polarization_uncertainty=target_uncertainty,
        )
    # endfor

    for period in PERIODS:
        if not any(
            record.period == period
            and record.charge_plus > 0.0
            and record.charge_minus > 0.0
            for record in records.values()
        ):
            raise RuntimeError(
                f"No active NH3 run records with positive Q+ and Q- were "
                f"found for {period}."
            )
        # endif
    # endfor

    return records


def run_state_arrays(
    run_records: Mapping[int, RunRecord],
) -> dict[str, dict[str, np.ndarray]]:
    result: dict[str, dict[str, np.ndarray]] = {}
    for period in PERIODS:
        period_records = sorted(
            (
                record
                for record in run_records.values()
                if (
                    record.period == period
                    and record.charge_plus > 0.0
                    and record.charge_minus > 0.0
                )
            ),
            key=lambda record: record.run,
        )
        result[period] = {
            "run": np.asarray([record.run for record in period_records], dtype=np.int32),
            "pt": np.asarray(
                [record.target_polarization for record in period_records],
                dtype=np.float64,
            ),
            "q_plus": np.asarray(
                [record.charge_plus for record in period_records],
                dtype=np.float64,
            ),
            "q_minus": np.asarray(
                [record.charge_minus for record in period_records],
                dtype=np.float64,
            ),
        }
    # endfor
    return result


# =============================================================================
# Channel-selection cuts
# =============================================================================

def load_nominal_cuts(path: Path) -> dict[tuple[str, int], CutRecord]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing channel-selection cut JSON: {path}")
    # endif

    payload = json.loads(path.read_text(encoding="utf-8"))
    period_payload = payload.get("periods")
    if not isinstance(period_payload, dict):
        raise RuntimeError(
            f"Cut JSON {path} does not contain a 'periods' mapping."
        )
    # endif

    cuts: dict[tuple[str, int], CutRecord] = {}
    for period in PERIODS:
        rows = period_payload.get(period)
        if not isinstance(rows, list):
            raise RuntimeError(f"Cut JSON has no row list for {period}.")
        # endif

        for row in rows:
            x_index = int(row["x_index"])
            t_index = int(row["t_index"])
            bin_number = int(
                row.get("bin_number", combined_bin_number(x_index, t_index))
            )
            nominal = row.get("nominal")
            if not isinstance(nominal, list) or len(nominal) != 2:
                raise RuntimeError(
                    f"Missing nominal cut interval for {period}, bin {bin_number}."
                )
            # endif
            low, high = float(nominal[0]), float(nominal[1])
            if not (
                math.isfinite(low)
                and math.isfinite(high)
                and high > low
            ):
                raise RuntimeError(
                    f"Invalid nominal cut for {period}, bin {bin_number}: "
                    f"[{low}, {high}]."
                )
            # endif
            cuts[(period, bin_number)] = CutRecord(
                period=period,
                bin_number=bin_number,
                x_index=x_index,
                t_index=t_index,
                low_gev2=low,
                high_gev2=high,
            )
        # endfor
    # endfor

    expected = len(PERIODS) * NUMBER_OF_BINS
    if len(cuts) != expected:
        raise RuntimeError(
            f"Expected {expected} period/bin cuts, found {len(cuts)}."
        )
    # endif
    return cuts


# =============================================================================
# Dilution factors
# =============================================================================

def find_default_dilution_json(directory: Path) -> Path:
    if not directory.is_dir():
        raise FileNotFoundError(
            f"Missing dilution-factor output directory: {directory}"
        )
    # endif

    preferred = [
        directory / "production/tables/dilution_factors_production.json",
        directory / "production/dilution_factors_production.json",
        directory / "analysis/dilution_factors_production.json",
    ]
    for path in preferred:
        if path.is_file():
            return path
        # endif
    # endfor

    candidates = sorted(
        directory.rglob("dilution_factors_production*.json"),
        key=lambda path: (path.stat().st_mtime, path.name),
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            "Could not locate dilution_factors_production*.json under "
            f"{directory}."
        )
    # endif
    return candidates[0]


def _extract_recommended_record(cut_payload: Mapping[str, Any]) -> tuple[float, float]:
    record = cut_payload.get("recommended")
    if not isinstance(record, Mapping):
        raise RuntimeError("Dilution cut payload has no 'recommended' record.")
    # endif

    value_keys = ("value", "central", "recommended_value")
    uncertainty_keys = (
        "stat_uncertainty",
        "statistical_uncertainty",
        "recommended_stat_uncertainty",
    )

    value = None
    for key in value_keys:
        if key in record:
            value = float(record[key])
            break
        # endif
    # endfor

    uncertainty = None
    for key in uncertainty_keys:
        if key in record:
            uncertainty = float(record[key])
            break
        # endif
    # endfor

    if value is None or uncertainty is None:
        raise RuntimeError(
            "Recommended dilution record lacks value/stat_uncertainty fields: "
            f"{dict(record)}"
        )
    # endif
    return value, uncertainty


def load_dilution_factors(path: Path) -> dict[tuple[str, int], DilutionRecord]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing dilution-factor JSON: {path}")
    # endif

    payload = json.loads(path.read_text(encoding="utf-8"))
    periods = payload.get("periods")
    if not isinstance(periods, dict):
        raise RuntimeError(
            f"Dilution JSON {path} does not contain a 'periods' mapping."
        )
    # endif

    records: dict[tuple[str, int], DilutionRecord] = {}
    for period in PERIODS:
        period_payload = periods.get(period)
        if not isinstance(period_payload, Mapping):
            raise RuntimeError(f"Dilution JSON has no payload for {period}.")
        # endif
        bins = period_payload.get("bins")
        if not isinstance(bins, list):
            raise RuntimeError(f"Dilution JSON has no bin list for {period}.")
        # endif

        for row in bins:
            x_index = int(row["x_index"])
            t_index = int(row["t_index"])
            bin_number = int(
                row.get("bin_number", combined_bin_number(x_index, t_index))
            )
            cuts = row.get("cuts")
            if not isinstance(cuts, Mapping) or "nominal" not in cuts:
                raise RuntimeError(
                    f"Dilution JSON lacks nominal cut for {period}, "
                    f"bin {bin_number}."
                )
            # endif
            value, uncertainty = _extract_recommended_record(cuts["nominal"])
            if not (
                math.isfinite(value)
                and math.isfinite(uncertainty)
                and value > 0.0
                and uncertainty >= 0.0
            ):
                raise RuntimeError(
                    f"Invalid dilution factor for {period}, bin {bin_number}: "
                    f"{value} +/- {uncertainty}."
                )
            # endif
            records[(period, bin_number)] = DilutionRecord(
                period=period,
                bin_number=bin_number,
                x_index=x_index,
                t_index=t_index,
                value=value,
                stat_uncertainty=uncertainty,
            )
        # endfor
    # endfor

    expected = len(PERIODS) * NUMBER_OF_BINS
    if len(records) != expected:
        raise RuntimeError(
            f"Expected {expected} dilution records, found {len(records)}."
        )
    # endif
    return records


# =============================================================================
# Input ROOT handling and event cache
# =============================================================================

def parse_input_override(text: str) -> tuple[str, Path]:
    try:
        period_text, path_text = text.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Input overrides must have the form PERIOD=/path/file.root"
        ) from exc
    # endtry
    period = period_text.strip().lower()
    if period not in PERIODS:
        raise argparse.ArgumentTypeError(
            f"Unknown period {period!r}; expected one of {PERIODS}."
        )
    # endif
    return period, Path(path_text).expanduser()


def resolve_tree_branches(
    path: Path,
    tree_name: str,
) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing ROOT input: {path}")
    # endif

    with uproot.open(path) as root_file:
        if tree_name not in root_file:
            available = [str(key).split(";")[0] for key in root_file.keys()]
            raise RuntimeError(
                f"Tree {tree_name!r} not found in {path}. "
                f"Available keys: {available}"
            )
        # endif
        tree = root_file[tree_name]
        return {
            logical_name: resolve_branch(tree, aliases)
            for logical_name, aliases in BRANCH_ALIASES.items()
        }
    # endwith


def compute_theta_gamma(
    electron_momentum_gev: np.ndarray,
    electron_theta_rad: np.ndarray,
    beam_energy_gev: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the virtual-photon angle relative to the incident beam.

    With the incident electron along +z,

        q_T = p_e' sin(theta_e),
        q_z = E_beam - p_e' cos(theta_e),

        sin(theta_gamma) = q_T / |q|,
        cos(theta_gamma) = q_z / |q|.

    This is algebraically equivalent to the analysis-note expression for
    theta_gamma and guarantees sin^2 + cos^2 = 1 numerically.
    """
    q_transverse = electron_momentum_gev * np.sin(electron_theta_rad)
    q_longitudinal = (
        beam_energy_gev
        - electron_momentum_gev * np.cos(electron_theta_rad)
    )
    q_magnitude = np.hypot(q_transverse, q_longitudinal)
    valid = np.isfinite(q_magnitude) & (q_magnitude > 0.0)

    sin_theta = np.full(q_magnitude.shape, np.nan, dtype=np.float64)
    cos_theta = np.full(q_magnitude.shape, np.nan, dtype=np.float64)
    sin_theta[valid] = q_transverse[valid] / q_magnitude[valid]
    cos_theta[valid] = q_longitudinal[valid] / q_magnitude[valid]

    sin_theta = np.clip(sin_theta, -1.0, 1.0)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return sin_theta, cos_theta


def build_event_cache(
    input_paths: Mapping[str, Path],
    tree_name: str,
    chunk_size: str,
    run_records: Mapping[int, RunRecord],
    cuts: Mapping[tuple[str, int], CutRecord],
    cache_path: Path,
) -> dict[str, Any]:
    ensure_directory(cache_path.parent)

    collected: dict[str, list[np.ndarray]] = {
        "period_index": [],
        "runnum": [],
        "helicity": [],
        "bin_number": [],
        "xB": [],
        "minus_tprime": [],
        "Mx2": [],
        "phi": [],
        "Q2": [],
        "sin_theta_gamma": [],
        "cos_theta_gamma": [],
        "rB": [],
        "rC": [],
        "rV": [],
        "rW": [],
    }

    period_statistics: dict[str, dict[str, int]] = {}

    for period in PERIODS:
        path = input_paths[period].expanduser().resolve()
        branches = resolve_tree_branches(path, tree_name)
        expressions = list(dict.fromkeys(branches.values()))
        total_seen = 0
        total_selected = 0
        total_zero_helicity = 0
        total_missing_run = 0
        total_zero_charge_run_events = 0
        missing_run_event_counts: dict[int, int] = {}
        zero_charge_run_event_counts: dict[int, int] = {}

        source = f"{path}:{tree_name}"
        for arrays in uproot.iterate(
            source,
            expressions=expressions,
            step_size=chunk_size,
            library="np",
        ):
            n_chunk = len(arrays[branches["runnum"]])
            total_seen += n_chunk

            runnum = np.asarray(
                arrays[branches["runnum"]],
                dtype=np.int64,
            )
            helicity = np.asarray(
                arrays[branches["helicity"]],
                dtype=np.int8,
            )
            xB = np.asarray(arrays[branches["xB"]], dtype=np.float64)
            minus_tprime = -np.asarray(
                arrays[branches["tprime"]],
                dtype=np.float64,
            )
            mx2 = np.asarray(arrays[branches["Mx2"]], dtype=np.float64)
            phi_raw = np.asarray(arrays[branches["phi"]], dtype=np.float64)
            phi = angle_to_radians(phi_raw, branches["phi"])
            phi = np.mod(phi, 2.0 * math.pi)
            q2 = np.asarray(arrays[branches["Q2"]], dtype=np.float64)
            e_p = np.asarray(arrays[branches["e_p"]], dtype=np.float64)
            e_theta_raw = np.asarray(
                arrays[branches["e_theta"]],
                dtype=np.float64,
            )
            e_theta = angle_to_radians(e_theta_raw, branches["e_theta"])

            dep_a = np.asarray(arrays[branches["DepA"]], dtype=np.float64)
            dep_b = np.asarray(arrays[branches["DepB"]], dtype=np.float64)
            dep_c = np.asarray(arrays[branches["DepC"]], dtype=np.float64)
            dep_v = np.asarray(arrays[branches["DepV"]], dtype=np.float64)
            dep_w = np.asarray(arrays[branches["DepW"]], dtype=np.float64)

            sin_theta_gamma, cos_theta_gamma = compute_theta_gamma(
                e_p,
                e_theta,
                BEAM_ENERGY_GEV[period],
            )

            _, _, bin_number = bin_indices(xB, minus_tprime)

            known_run = np.fromiter(
                (
                    int(run) in run_records
                    and run_records[int(run)].period == period
                    for run in runnum
                ),
                count=runnum.size,
                dtype=bool,
            )
            active_run = np.fromiter(
                (
                    int(run) in run_records
                    and run_records[int(run)].period == period
                    and run_records[int(run)].charge_plus > 0.0
                    and run_records[int(run)].charge_minus > 0.0
                    for run in runnum
                ),
                count=runnum.size,
                dtype=bool,
            )
            zero_charge_run = known_run & ~active_run

            missing_mask = ~known_run
            total_missing_run += int(np.count_nonzero(missing_mask))
            total_zero_charge_run_events += int(
                np.count_nonzero(zero_charge_run)
            )
            total_zero_helicity += int(np.count_nonzero(helicity == 0))

            if np.any(missing_mask):
                missing_runs, missing_counts = np.unique(
                    runnum[missing_mask],
                    return_counts=True,
                )
                for missing_run, missing_count in zip(
                    missing_runs,
                    missing_counts,
                ):
                    key = int(missing_run)
                    missing_run_event_counts[key] = (
                        missing_run_event_counts.get(key, 0)
                        + int(missing_count)
                    )
                # endfor
            # endif

            if np.any(zero_charge_run):
                disabled_runs, disabled_counts = np.unique(
                    runnum[zero_charge_run],
                    return_counts=True,
                )
                for disabled_run, disabled_count in zip(
                    disabled_runs,
                    disabled_counts,
                ):
                    key = int(disabled_run)
                    zero_charge_run_event_counts[key] = (
                        zero_charge_run_event_counts.get(key, 0)
                        + int(disabled_count)
                    )
                # endfor
            # endif

            finite = (
                np.isfinite(xB)
                & np.isfinite(minus_tprime)
                & np.isfinite(mx2)
                & np.isfinite(phi)
                & np.isfinite(q2)
                & np.isfinite(sin_theta_gamma)
                & np.isfinite(cos_theta_gamma)
                & np.isfinite(dep_a)
                & np.isfinite(dep_b)
                & np.isfinite(dep_c)
                & np.isfinite(dep_v)
                & np.isfinite(dep_w)
            )
            base = (
                finite
                & active_run
                & np.isin(helicity, (-1, 1))
                & (bin_number >= 1)
                & (dep_a > 0.0)
            )

            selected = np.zeros(base.shape, dtype=bool)
            for current_bin in range(1, NUMBER_OF_BINS + 1):
                cut = cuts[(period, current_bin)]
                selected |= (
                    base
                    & (bin_number == current_bin)
                    & (mx2 >= cut.low_gev2)
                    & (mx2 < cut.high_gev2)
                )
            # endfor

            if not np.any(selected):
                continue
            # endif

            total_selected += int(np.count_nonzero(selected))
            collected["period_index"].append(
                np.full(
                    int(np.count_nonzero(selected)),
                    PERIOD_INDEX[period],
                    dtype=np.int8,
                )
            )
            collected["runnum"].append(runnum[selected].astype(np.int32))
            collected["helicity"].append(helicity[selected].astype(np.int8))
            collected["bin_number"].append(
                bin_number[selected].astype(np.int16)
            )
            collected["xB"].append(xB[selected])
            collected["minus_tprime"].append(minus_tprime[selected])
            collected["Mx2"].append(mx2[selected])
            collected["phi"].append(phi[selected])
            collected["Q2"].append(q2[selected])
            collected["sin_theta_gamma"].append(sin_theta_gamma[selected])
            collected["cos_theta_gamma"].append(cos_theta_gamma[selected])
            collected["rB"].append(dep_b[selected] / dep_a[selected])
            collected["rC"].append(dep_c[selected] / dep_a[selected])
            collected["rV"].append(dep_v[selected] / dep_a[selected])
            collected["rW"].append(dep_w[selected] / dep_a[selected])
        # endfor

        if missing_run_event_counts:
            details = ", ".join(
                f"{run}: {count:,} events"
                for run, count in sorted(missing_run_event_counts.items())
            )
            raise RuntimeError(
                f"{period} ROOT input contains events from runs absent from "
                f"the NH3 section of the run-information CSV: {details}."
            )
        # endif

        if zero_charge_run_event_counts:
            details = ", ".join(
                f"{run}: {count:,} events"
                for run, count in sorted(
                    zero_charge_run_event_counts.items()
                )
            )
            raise RuntimeError(
                f"{period} ROOT input contains events from QA-disabled runs "
                f"whose CSV charges are Q+=Q-=0: {details}. Zero-charge rows "
                "are allowed in clas12_run_info.csv, but no event from those "
                "runs may enter the ROOT input used for asymmetry extraction."
            )
        # endif

        period_statistics[period] = {
            "events_seen": total_seen,
            "events_selected": total_selected,
            "helicity_zero_events_ignored": total_zero_helicity,
            "events_with_missing_run_info": total_missing_run,
            "events_from_zero_charge_runs": total_zero_charge_run_events,
            "active_runs_in_likelihood": sum(
                1
                for record in run_records.values()
                if (
                    record.period == period
                    and record.charge_plus > 0.0
                    and record.charge_minus > 0.0
                )
            ),
        }
    # endfor

    if not collected["runnum"]:
        raise RuntimeError("No events survived the nominal selection.")
    # endif

    cache = {
        name: np.concatenate(chunks)
        for name, chunks in collected.items()
    }
    np.savez_compressed(cache_path, **cache)

    return {
        "cache_path": str(cache_path.resolve()),
        "number_of_selected_events": int(cache["runnum"].size),
        "period_statistics": period_statistics,
    }


def load_event_cache(path: Path) -> dict[str, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing event cache: {path}")
    # endif
    with np.load(path, allow_pickle=False) as payload:
        result = {
            key: np.asarray(payload[key])
            for key in payload.files
        }
    # endwith

    sizes = {array.shape[0] for array in result.values()}
    if len(sizes) != 1:
        raise RuntimeError(
            f"Cache arrays have inconsistent lengths: {sorted(sizes)}"
        )
    # endif
    return result


# =============================================================================
# Likelihood
# =============================================================================

def longitudinal_scaled_transverse_terms(
    phi: np.ndarray,
    r_b: np.ndarray,
    r_c: np.ndarray,
    r_v: np.ndarray,
    r_w: np.ndarray,
    scales: Mapping[str, float] | None,
    sign: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the transverse-target stress-model angular factors.

    The scales are taken from the nominal fitted longitudinal amplitudes in the
    same kinematic bin.  This is a deliberately loose, data-driven estimate of
    the possible transverse leakage; it is not an assertion that the transverse
    amplitudes equal the longitudinal amplitudes.

    Mapping used:
      UT sin(phi-phi_S)       <- UL sin(phi)
      UT sin(phi+phi_S)       <- UL sin(phi)
      UT sin(3phi-phi_S)      <- UL sin(2phi)
      UT sin(phi_S)           <- UL sin(phi) [vanishes at phi_S = 0]
      UT sin(2phi-phi_S)      <- UL sin(2phi)

      LT cos(phi-phi_S)       <- LL cos(phi)
      LT cos(phi_S)           <- LL constant
      LT cos(2phi-phi_S)      <- LL cos(phi)

    The plus stress test preserves the observed nominal signs.  The minus
    stress test reverses every mapped transverse amplitude coherently.
    """
    if scales is None:
        zeros = np.zeros(phi.shape, dtype=np.float64)
        return zeros, zeros
    # endif

    ul1_scale = sign * float(scales["ul1"])
    ul2_scale = sign * float(scales["ul2"])
    ll0_scale = sign * float(scales["ll0"])
    ll1_scale = sign * float(scales["ll1"])

    # Beam-axis target: phi_S = 0.  The sin(phi_S) term vanishes.
    ut = (
        ul1_scale * np.sin(phi)
        + r_b * ul1_scale * np.sin(phi)
        + r_b * ul2_scale * np.sin(3.0 * phi)
        + r_v * ul2_scale * np.sin(2.0 * phi)
    )
    lt = (
        r_c * ll1_scale * np.cos(phi)
        + r_w * ll0_scale
        + r_w * ll1_scale * np.cos(2.0 * phi)
    )
    return ut, lt


def evaluate_cross_section_factor(
    *,
    variant: str,
    phi: np.ndarray,
    r_b: np.ndarray,
    r_c: np.ndarray,
    r_v: np.ndarray,
    r_w: np.ndarray,
    sin_theta_gamma: np.ndarray,
    cos_theta_gamma: np.ndarray,
    helicity: np.ndarray | float,
    beam_polarization: float,
    target_polarization: np.ndarray | float,
    dilution: float,
    transverse_scales: Mapping[str, float] | None,
    u1: float,
    u2: float,
    lu1: float,
    ul1: float,
    ul2: float,
    ll0: float,
    ll1: float,
) -> np.ndarray:
    h = np.asarray(helicity, dtype=np.float64)
    pt = np.asarray(target_polarization, dtype=np.float64)

    if variant == "no_projection":
        p_longitudinal = pt
        p_transverse = np.zeros_like(pt, dtype=np.float64)
    else:
        p_longitudinal = pt * cos_theta_gamma
        p_transverse = pt * sin_theta_gamma
    # endif

    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cos_2phi = np.cos(2.0 * phi)
    sin_2phi = np.sin(2.0 * phi)

    unpolarized = 1.0 + r_v * u1 * cos_phi + r_b * u2 * cos_2phi
    beam_spin = h * beam_polarization * r_w * lu1 * sin_phi
    target_longitudinal = (
        dilution
        * p_longitudinal
        * (r_v * ul1 * sin_phi + r_b * ul2 * sin_2phi)
    )
    double_longitudinal = (
        h
        * beam_polarization
        * dilution
        * p_longitudinal
        * (r_c * ll0 + r_w * ll1 * cos_phi)
    )

    if variant in (
        "longitudinal_scaled_transverse_plus",
        "longitudinal_scaled_transverse_minus",
    ):
        transverse_sign = (
            1.0
            if variant == "longitudinal_scaled_transverse_plus"
            else -1.0
        )
        ut_fixed, lt_fixed = longitudinal_scaled_transverse_terms(
            phi,
            r_b,
            r_c,
            r_v,
            r_w,
            transverse_scales,
            transverse_sign,
        )
    else:
        ut_fixed = np.zeros(phi.shape, dtype=np.float64)
        lt_fixed = np.zeros(phi.shape, dtype=np.float64)
    # endif

    target_transverse = dilution * p_transverse * ut_fixed
    double_transverse = (
        h * beam_polarization * dilution * p_transverse * lt_fixed
    )

    return (
        unpolarized
        + beam_spin
        + target_longitudinal
        + double_longitudinal
        + target_transverse
        + double_transverse
    )


def make_bin_nll(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    bin_number: int,
    variant: str,
    active_periods: tuple[str, ...] = PERIODS,
    transverse_scales: Mapping[str, float] | None = None,
):
    mask = events["bin_number"] == bin_number
    if len(active_periods) != len(PERIODS):
        allowed_indices = np.asarray(
            [PERIOD_INDEX[period] for period in active_periods],
            dtype=np.int8,
        )
        mask &= np.isin(events["period_index"], allowed_indices)
    # endif
    if not np.any(mask):
        raise RuntimeError(
            f"Bin {bin_number} has no selected events for {active_periods}."
        )
    # endif

    period_index = events["period_index"][mask].astype(np.int8, copy=False)
    runnum = events["runnum"][mask].astype(np.int32, copy=False)
    helicity = events["helicity"][mask].astype(np.int8, copy=False)
    phi = events["phi"][mask].astype(np.float64, copy=False)
    sin_theta = events["sin_theta_gamma"][mask].astype(np.float64, copy=False)
    cos_theta = events["cos_theta_gamma"][mask].astype(np.float64, copy=False)
    r_b = events["rB"][mask].astype(np.float64, copy=False)
    r_c = events["rC"][mask].astype(np.float64, copy=False)
    r_v = events["rV"][mask].astype(np.float64, copy=False)
    r_w = events["rW"][mask].astype(np.float64, copy=False)

    period_event_indices = {
        period: np.flatnonzero(period_index == PERIOD_INDEX[period])
        for period in active_periods
    }

    run_lookup = {
        period: {
            int(run): index
            for index, run in enumerate(run_states[period]["run"])
        }
        for period in active_periods
    }

    observed_run_indices: dict[str, np.ndarray] = {}
    for period in active_periods:
        indices = period_event_indices[period]
        observed_run_indices[period] = np.fromiter(
            (run_lookup[period][int(run)] for run in runnum[indices]),
            count=indices.size,
            dtype=np.int32,
        )
    # endfor

    dilution_central = {
        period: dilution_records[(period, bin_number)].value
        for period in active_periods
    }
    dilution_sigma = {
        period: dilution_records[(period, bin_number)].stat_uncertainty
        for period in active_periods
    }

    # Precompute all parameter-independent event arrays and run-state charge
    # moments once.  The previous implementation recomputed several of these
    # quantities at every Minuit function call.
    period_data: dict[str, dict[str, Any]] = {}
    for period in active_periods:
        indices = period_event_indices[period]
        event_phi = phi[indices]
        state = run_states[period]
        state_pt = state["pt"]
        state_q_plus = state["q_plus"]
        state_q_minus = state["q_minus"]
        observed_state_index = observed_run_indices[period]
        event_h = helicity[indices].astype(np.float64, copy=False)

        period_data[period] = {
            "phi": event_phi,
            "sin_phi": np.sin(event_phi),
            "cos_phi": np.cos(event_phi),
            "sin_2phi": np.sin(2.0 * event_phi),
            "cos_2phi": np.cos(2.0 * event_phi),
            "sin_3phi": np.sin(3.0 * event_phi),
            "r_b": r_b[indices],
            "r_c": r_c[indices],
            "r_v": r_v[indices],
            "r_w": r_w[indices],
            "sin_theta": sin_theta[indices],
            "cos_theta": cos_theta[indices],
            "h": event_h,
            "observed_pt": state_pt[observed_state_index],
            "observed_charge": np.where(
                event_h > 0.0,
                state_q_plus[observed_state_index],
                state_q_minus[observed_state_index],
            ),
            "charge_sum": float(np.sum(state_q_plus + state_q_minus)),
            "helicity_charge_sum": float(
                np.sum(state_q_plus - state_q_minus)
            ),
            "target_charge_sum": float(
                np.sum(state_pt * (state_q_plus + state_q_minus))
            ),
            "helicity_target_charge_sum": float(
                np.sum(state_pt * (state_q_plus - state_q_minus))
            ),
        }
    # endfor

    def nll(
        u1: float,
        u2: float,
        lu1: float,
        ul1: float,
        ul2: float,
        ll0: float,
        ll1: float,
        f_su22: float,
        f_fa22: float,
        f_sp23: float,
    ) -> float:
        parameter_values = (u1, u2, lu1, ul1, ul2, ll0, ll1)
        if not all(math.isfinite(value) for value in parameter_values):
            return INVALID_NLL
        # endif

        dilution_by_period = {
            "su22": f_su22,
            "fa22": f_fa22,
            "sp23": f_sp23,
        }
        total = 0.0

        for period in active_periods:
            dilution = dilution_by_period[period]
            sigma = dilution_sigma[period]
            central = dilution_central[period]
            if not math.isfinite(dilution) or dilution <= 0.0:
                return INVALID_NLL
            # endif
            if sigma > 0.0:
                total += 0.5 * ((dilution - central) / sigma) ** 2
            elif abs(dilution - central) > 1.0e-14:
                return INVALID_NLL
            # endif

            data = period_data[period]
            event_phi = data["phi"]
            event_r_b = data["r_b"]
            event_r_c = data["r_c"]
            event_r_v = data["r_v"]
            event_r_w = data["r_w"]
            event_sin = data["sin_theta"]
            event_cos = data["cos_theta"]
            event_h = data["h"]

            numerator_factor = evaluate_cross_section_factor(
                variant=variant,
                phi=event_phi,
                r_b=event_r_b,
                r_c=event_r_c,
                r_v=event_r_v,
                r_w=event_r_w,
                sin_theta_gamma=event_sin,
                cos_theta_gamma=event_cos,
                helicity=event_h,
                beam_polarization=BEAM_POLARIZATION[period],
                target_polarization=data["observed_pt"],
                dilution=dilution,
                transverse_scales=transverse_scales,
                u1=u1,
                u2=u2,
                lu1=lu1,
                ul1=ul1,
                ul2=ul2,
                ll0=ll0,
                ll1=ll1,
            )
            if (
                np.any(~np.isfinite(numerator_factor))
                or np.any(numerator_factor <= CROSS_SECTION_FLOOR)
            ):
                return INVALID_NLL
            # endif

            unpolarized = (
                1.0
                + event_r_v * u1 * data["cos_phi"]
                + event_r_b * u2 * data["cos_2phi"]
            )
            beam_coefficient = (
                BEAM_POLARIZATION[period]
                * event_r_w
                * lu1
                * data["sin_phi"]
            )

            if variant == "no_projection":
                longitudinal_geometry = np.ones_like(event_phi)
                transverse_geometry = np.zeros_like(event_phi)
            else:
                longitudinal_geometry = event_cos
                transverse_geometry = event_sin
            # endif

            target_longitudinal_coefficient = (
                dilution
                * longitudinal_geometry
                * (
                    event_r_v * ul1 * data["sin_phi"]
                    + event_r_b * ul2 * data["sin_2phi"]
                )
            )
            double_longitudinal_coefficient = (
                BEAM_POLARIZATION[period]
                * dilution
                * longitudinal_geometry
                * (
                    event_r_c * ll0
                    + event_r_w * ll1 * data["cos_phi"]
                )
            )

            if variant in (
                "longitudinal_scaled_transverse_plus",
                "longitudinal_scaled_transverse_minus",
            ):
                scales = transverse_scales or {}
                transverse_sign = (
                    1.0
                    if variant == "longitudinal_scaled_transverse_plus"
                    else -1.0
                )
                ul1_scale = transverse_sign * float(scales["ul1"])
                ul2_scale = transverse_sign * float(scales["ul2"])
                ll0_scale = transverse_sign * float(scales["ll0"])
                ll1_scale = transverse_sign * float(scales["ll1"])
                ut_fixed = (
                    ul1_scale * data["sin_phi"]
                    + event_r_b * ul1_scale * data["sin_phi"]
                    + event_r_b * ul2_scale * data["sin_3phi"]
                    + event_r_v * ul2_scale * data["sin_2phi"]
                )
                lt_fixed = (
                    event_r_c * ll1_scale * data["cos_phi"]
                    + event_r_w * ll0_scale
                    + event_r_w * ll1_scale * data["cos_2phi"]
                )
                target_transverse_coefficient = (
                    dilution * transverse_geometry * ut_fixed
                )
                double_transverse_coefficient = (
                    BEAM_POLARIZATION[period]
                    * dilution
                    * transverse_geometry
                    * lt_fixed
                )
            else:
                target_transverse_coefficient = 0.0
                double_transverse_coefficient = 0.0
            # endif

            target_coefficient = (
                target_longitudinal_coefficient
                + target_transverse_coefficient
            )
            double_coefficient = (
                double_longitudinal_coefficient
                + double_transverse_coefficient
            )

            denominator = (
                data["charge_sum"] * unpolarized
                + data["helicity_charge_sum"] * beam_coefficient
                + data["target_charge_sum"] * target_coefficient
                + data["helicity_target_charge_sum"] * double_coefficient
            )
            if (
                np.any(~np.isfinite(denominator))
                or np.any(denominator <= CROSS_SECTION_FLOOR)
            ):
                return INVALID_NLL
            # endif

            numerator = data["observed_charge"] * numerator_factor
            probability = numerator / denominator
            if (
                np.any(~np.isfinite(probability))
                or np.any(probability <= 0.0)
                or np.any(probability > 1.0 + 1.0e-10)
            ):
                return INVALID_NLL
            # endif
            total -= float(
                np.sum(np.log(np.maximum(probability, PROBABILITY_FLOOR)))
            )
        # endfor

        return total

    metadata = {
        "active_periods": list(active_periods),
        "number_of_events": int(np.count_nonzero(mask)),
        "events_by_period": {
            period: int(period_event_indices.get(period, np.empty(0)).size)
            for period in PERIODS
        },
        "dilution_central": {
            period: (
                dilution_records[(period, bin_number)].value
                if period in active_periods else None
            )
            for period in PERIODS
        },
        "dilution_stat_uncertainty": {
            period: (
                dilution_records[(period, bin_number)].stat_uncertainty
                if period in active_periods else None
            )
            for period in PERIODS
        },
        "mean_sin_theta_gamma": float(np.mean(sin_theta)),
        "rms_sin_theta_gamma": float(np.std(sin_theta, ddof=1))
        if sin_theta.size > 1 else 0.0,
        "mean_cos_theta_gamma": float(np.mean(cos_theta)),
        "rms_cos_theta_gamma": float(np.std(cos_theta, ddof=1))
        if cos_theta.size > 1 else 0.0,
        "mean_xB": float(np.mean(events["xB"][mask])),
        "mean_minus_tprime_gev2": float(
            np.mean(events["minus_tprime"][mask])
        ),
        "mean_Q2_gev2": float(np.mean(events["Q2"][mask])),
    }
    return nll, metadata


def minuit_result_payload(minuit: Minuit) -> dict[str, Any]:
    parameter_names = list(minuit.parameters)
    values = {
        name: float(minuit.values[name])
        for name in parameter_names
    }
    errors = {
        name: float(minuit.errors[name])
        for name in parameter_names
    }

    covariance = None
    correlation = None
    if minuit.covariance is not None:
        covariance = [
            [
                float(minuit.covariance[name_i, name_j])
                for name_j in parameter_names
            ]
            for name_i in parameter_names
        ]
        diagonal = np.sqrt(
            np.maximum(np.diag(np.asarray(covariance)), 0.0)
        )
        matrix = np.asarray(covariance, dtype=np.float64)
        denominator = np.outer(diagonal, diagonal)
        correlation_array = np.divide(
            matrix,
            denominator,
            out=np.zeros_like(matrix),
            where=denominator > 0.0,
        )
        correlation = correlation_array.tolist()
    # endif

    fmin = minuit.fmin
    return {
        "parameter_order": parameter_names,
        "values": values,
        "errors": errors,
        "covariance": covariance,
        "correlation": correlation,
        "minimum_nll": float(minuit.fval),
        "valid": bool(fmin.is_valid),
        "accurate_covariance": bool(fmin.has_accurate_covar),
        "positive_definite_covariance": bool(fmin.has_posdef_covar),
        "parameters_at_limit": bool(fmin.has_parameters_at_limit),
        "edm": float(fmin.edm),
        "edm_goal": float(fmin.edm_goal),
        "nfcn": int(fmin.nfcn),
    }


def fit_one_variant(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    bin_number: int,
    variant: str,
    active_periods: tuple[str, ...] = PERIODS,
    transverse_scales: Mapping[str, float] | None = None,
    initial_values: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    nll, metadata = make_bin_nll(
        events,
        run_states,
        dilution_records,
        bin_number,
        variant,
        active_periods=active_periods,
        transverse_scales=transverse_scales,
    )

    initial = dict(PARAMETER_INITIAL_VALUES)
    if initial_values is not None:
        for name in PHYSICS_PARAMETERS:
            if name in initial_values:
                initial[name] = float(initial_values[name])
            # endif
        # endfor
    # endif
    initial.update(
        {
            f"f_{period}": dilution_records[(period, bin_number)].value
            for period in PERIODS
        }
    )

    def configured_minuit(start_values: Mapping[str, float]) -> Minuit:
        candidate = Minuit(nll, **dict(start_values))
        candidate.errordef = Minuit.LIKELIHOOD
        candidate.print_level = 0
        candidate.strategy = 1

        for parameter_name in PHYSICS_PARAMETERS:
            candidate.limits[parameter_name] = PARAMETER_LIMITS[
                parameter_name
            ]
        # endfor
        for period in PERIODS:
            nuisance_name = f"f_{period}"
            central = dilution_records[(period, bin_number)].value
            sigma = dilution_records[(period, bin_number)].stat_uncertainty
            width = max(8.0 * sigma, 0.20 * central, 0.02)
            candidate.limits[nuisance_name] = (
                max(1.0e-6, central - width),
                central + width,
            )
            if period not in active_periods or sigma == 0.0:
                candidate.fixed[nuisance_name] = True
            # endif
        # endfor
        return candidate

    start_candidates: list[dict[str, float]] = [dict(initial)]
    zero_start = dict(initial)
    zero_start.update({name: 0.0 for name in PHYSICS_PARAMETERS})
    start_candidates.append(zero_start)

    # A deterministic midpoint start is useful for difficult low-statistics
    # bins without introducing run-to-run randomness.
    midpoint_start = dict(initial)
    midpoint_start.update(
        {
            name: 0.5 * float(initial[name])
            for name in PHYSICS_PARAMETERS
        }
    )
    start_candidates.append(midpoint_start)

    attempted: list[Minuit] = []
    for attempt_index, start_values in enumerate(start_candidates):
        candidate = configured_minuit(start_values)
        candidate.migrad(ncall=50000)
        if not candidate.fmin.is_valid:
            candidate.simplex(ncall=50000)
            candidate.strategy = 2
            candidate.migrad(ncall=100000)
        # endif
        attempted.append(candidate)
        if candidate.fmin.is_valid:
            break
        # endif
    # endfor

    valid_attempts = [
        candidate for candidate in attempted if candidate.fmin.is_valid
    ]
    if valid_attempts:
        minuit = min(valid_attempts, key=lambda candidate: candidate.fval)
    else:
        minuit = min(attempted, key=lambda candidate: candidate.fval)
    # endif
    minuit.hesse()

    result = minuit_result_payload(minuit)
    result.update(
        {
            "bin_number": bin_number,
            "variant": variant,
            "active_periods": list(active_periods),
            "transverse_scales": (
                dict(transverse_scales)
                if transverse_scales is not None else None
            ),
            "metadata": metadata,
        }
    )
    return result


_WORKER_EVENTS: dict[str, np.ndarray] | None = None
_WORKER_RUN_STATES: dict[str, dict[str, np.ndarray]] | None = None
_WORKER_DILUTION_RECORDS: dict[tuple[str, int], DilutionRecord] | None = None


def initialize_fit_worker(
    cache_path_text: str,
    run_state_payload: dict[str, dict[str, list[float] | list[int]]],
    dilution_payload: dict[str, dict[str, dict[str, float | int]]],
) -> None:
    """Load immutable fit inputs once per worker process."""
    global _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS

    _WORKER_EVENTS = load_event_cache(Path(cache_path_text))
    _WORKER_RUN_STATES = {
        period: {
            key: np.asarray(values)
            for key, values in state.items()
        }
        for period, state in run_state_payload.items()
    }
    _WORKER_DILUTION_RECORDS = {
        (period, int(bin_text)): DilutionRecord(
            period=period,
            bin_number=int(bin_text),
            x_index=int(record["x_index"]),
            t_index=int(record["t_index"]),
            value=float(record["value"]),
            stat_uncertainty=float(record["stat_uncertainty"]),
        )
        for period, period_payload in dilution_payload.items()
        for bin_text, record in period_payload.items()
    }


def fit_bin_worker(bin_number: int) -> dict[str, Any]:
    if (
        _WORKER_EVENTS is None
        or _WORKER_RUN_STATES is None
        or _WORKER_DILUTION_RECORDS is None
    ):
        raise RuntimeError("Fit worker was not initialized.")
    # endif

    events = _WORKER_EVENTS
    run_states = _WORKER_RUN_STATES
    dilution_records = _WORKER_DILUTION_RECORDS

    nominal = fit_one_variant(
        events,
        run_states,
        dilution_records,
        bin_number,
        "nominal",
    )
    no_projection = fit_one_variant(
        events,
        run_states,
        dilution_records,
        bin_number,
        "no_projection",
        initial_values=nominal["values"],
    )

    transverse_scales = {
        "ul1": nominal["values"]["ul1"],
        "ul2": nominal["values"]["ul2"],
        "ll0": nominal["values"]["ll0"],
        "ll1": nominal["values"]["ll1"],
    }
    longitudinal_scaled_transverse_plus = fit_one_variant(
        events,
        run_states,
        dilution_records,
        bin_number,
        "longitudinal_scaled_transverse_plus",
        transverse_scales=transverse_scales,
        initial_values=nominal["values"],
    )
    longitudinal_scaled_transverse_minus = fit_one_variant(
        events,
        run_states,
        dilution_records,
        bin_number,
        "longitudinal_scaled_transverse_minus",
        transverse_scales=transverse_scales,
        initial_values=nominal["values"],
    )

    variants = {
        "nominal": nominal,
        "no_projection": no_projection,
        "longitudinal_scaled_transverse_plus": (
            longitudinal_scaled_transverse_plus
        ),
        "longitudinal_scaled_transverse_minus": (
            longitudinal_scaled_transverse_minus
        ),
    }

    systematics: dict[str, float] = {}
    for parameter in PHYSICS_PARAMETERS:
        nominal_value = nominal["values"][parameter]
        systematics[parameter] = max(
            abs(no_projection["values"][parameter] - nominal_value),
            abs(
                longitudinal_scaled_transverse_plus["values"][parameter]
                - nominal_value
            ),
            abs(
                longitudinal_scaled_transverse_minus["values"][parameter]
                - nominal_value
            ),
        )
    # endfor

    # Independent nominal fits for period-consistency plots.  These do not
    # enter the quoted combined result or its target-axis systematic.
    period_fits = {
        period: fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
            initial_values=nominal["values"],
        )
        for period in PERIODS
    }

    # Quantify each period's tension with the simultaneous solution.  For each
    # period, compare its nominal NLL at the combined best-fit point with the
    # independently minimized period-only NLL.  The difference is nonnegative
    # up to minimizer precision and is a compact period-consistency diagnostic.
    period_consistency: dict[str, dict[str, float]] = {}
    for period in PERIODS:
        period_nll, _ = make_bin_nll(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
        )
        combined_nll_for_period = float(
            period_nll(
                **{
                    name: nominal["values"][name]
                    for name in (
                        *PHYSICS_PARAMETERS,
                        "f_su22",
                        "f_fa22",
                        "f_sp23",
                    )
                }
            )
        )
        period_minimum_nll = float(period_fits[period]["minimum_nll"])
        period_consistency[period] = {
            "nll_at_combined_solution": combined_nll_for_period,
            "period_only_minimum_nll": period_minimum_nll,
            "delta_nll": max(
                0.0,
                combined_nll_for_period - period_minimum_nll,
            ),
        }
    # endfor

    return {
        "bin_number": bin_number,
        "variants": variants,
        "period_fits": period_fits,
        "period_consistency": period_consistency,
        "target_axis_systematic": systematics,
        "transverse_scales": transverse_scales,
    }


# =============================================================================
# Outputs
# =============================================================================

def flatten_fit_results(
    results: list[dict[str, Any]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for result in sorted(results, key=lambda item: item["bin_number"]):
        nominal = result["variants"]["nominal"]
        metadata = nominal["metadata"]
        row: dict[str, Any] = {
            "bin_number": result["bin_number"],
            "x_index": (result["bin_number"] - 1)
            // len(MINUS_TPRIME_BINS_GEV2),
            "t_index": (result["bin_number"] - 1)
            % len(MINUS_TPRIME_BINS_GEV2),
            "mean_xB": metadata["mean_xB"],
            "mean_Q2_gev2": metadata["mean_Q2_gev2"],
            "mean_minus_tprime_gev2": metadata[
                "mean_minus_tprime_gev2"
            ],
            "number_of_events": metadata["number_of_events"],
            "mean_sin_theta_gamma": metadata["mean_sin_theta_gamma"],
            "rms_sin_theta_gamma": metadata["rms_sin_theta_gamma"],
            "mean_cos_theta_gamma": metadata["mean_cos_theta_gamma"],
            "rms_cos_theta_gamma": metadata["rms_cos_theta_gamma"],
            "nominal_fit_valid": nominal["valid"],
            "nominal_minimum_nll": nominal["minimum_nll"],
            "nominal_edm": nominal["edm"],
        }

        for period in PERIODS:
            row[f"events_{period}"] = metadata["events_by_period"][period]
            row[f"dilution_{period}"] = metadata[
                "dilution_central"
            ][period]
            row[f"dilution_stat_{period}"] = metadata[
                "dilution_stat_uncertainty"
            ][period]
            row[f"fitted_dilution_{period}"] = nominal[
                "values"
            ][f"f_{period}"]
            row[f"fitted_dilution_error_{period}"] = nominal[
                "errors"
            ][f"f_{period}"]

            period_fit = result["period_fits"][period]
            row[f"period_fit_valid_{period}"] = period_fit["valid"]
            row[f"period_fit_nll_{period}"] = period_fit["minimum_nll"]
            row[f"period_fit_edm_{period}"] = period_fit["edm"]
            row[f"period_delta_nll_{period}"] = result[
                "period_consistency"
            ][period]["delta_nll"]
        # endfor

        for parameter in PHYSICS_PARAMETERS:
            row[parameter] = nominal["values"][parameter]
            row[f"{parameter}_stat"] = nominal["errors"][parameter]
            row[f"{parameter}_target_axis_sys"] = result[
                "target_axis_systematic"
            ][parameter]
            for variant in FIT_VARIANTS:
                row[f"{parameter}_{variant}"] = result[
                    "variants"
                ][variant]["values"][parameter]
            # endfor
            for period in PERIODS:
                period_fit = result["period_fits"][period]
                row[f"{parameter}_{period}"] = period_fit[
                    "values"
                ][parameter]
                row[f"{parameter}_stat_{period}"] = period_fit[
                    "errors"
                ][parameter]
            # endfor
        # endfor
        rows.append(row)
    # endfor
    return pd.DataFrame(rows)


def apply_parameter_y_limits(ax: plt.Axes, parameter: str) -> None:
    limits = PARAMETER_Y_LIMITS[parameter]
    if limits is not None:
        ax.set_ylim(*limits)
    # endif


def plot_parameter_summaries(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """Write the original one-parameter, all-24-bin summary plots."""
    ensure_directory(output_dir)
    paths: list[str] = []
    bins = frame["bin_number"].to_numpy()

    for parameter in PHYSICS_PARAMETERS:
        fig, ax = plt.subplots(figsize=(14, 5.5))
        ax.errorbar(
            bins,
            frame[parameter],
            yerr=frame[f"{parameter}_stat"],
            marker="o",
            linestyle="none",
            capsize=2,
            label="Statistical uncertainty",
        )
        ax.errorbar(
            bins,
            frame[parameter],
            yerr=frame[f"{parameter}_target_axis_sys"],
            marker="none",
            linestyle="none",
            capsize=4,
            linewidth=1.0,
            label="Target-axis envelope",
        )
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("Combined kinematic-bin number")
        ax.set_ylabel(PARAMETER_LABELS[parameter])
        ax.set_xticks(bins)
        apply_parameter_y_limits(ax, parameter)
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        path = output_dir / f"{parameter}_summary_v5.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths


def plot_aggregated_by_x(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """
    Write one seven-panel canvas per xB bin, showing the combined result versus
    mean -tprime.  Statistical and target-axis uncertainties are drawn
    separately.
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        subset = frame.loc[frame["x_index"] == x_index].sort_values("t_index")
        fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True)
        axes_flat = axes.ravel()

        for panel, parameter in enumerate(PHYSICS_PARAMETERS):
            ax = axes_flat[panel]
            x_values = subset["mean_minus_tprime_gev2"]
            ax.errorbar(
                x_values,
                subset[parameter],
                yerr=subset[f"{parameter}_stat"],
                marker="o",
                linestyle="none",
                capsize=2,
                label="Statistical uncertainty",
            )
            ax.errorbar(
                x_values,
                subset[parameter],
                yerr=subset[f"{parameter}_target_axis_sys"],
                marker="none",
                linestyle="none",
                capsize=4,
                linewidth=1.0,
                label="Target-axis envelope",
            )
            ax.axhline(0.0, linewidth=0.8)
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            apply_parameter_y_limits(ax, parameter)
            ax.grid(alpha=0.25)
            if panel >= 6:
                ax.set_xlabel(r"Mean $-t^\prime$ (GeV$^2$)")
            # endif
        # endfor

        axes_flat[7].axis("off")
        axes_flat[8].axis("off")
        handles, labels = axes_flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.96, 0.06))
        fig.suptitle(
            rf"${x_low:.2f} \leq x_B < {x_high:.2f}$",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
        path = output_dir / f"xB_bin_{x_index + 1}_combined_v5.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths


def plot_aggregated_by_period(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """
    Write one seven-panel canvas per xB bin with independent Su22, Fa22, Sp23,
    and simultaneous-combination nominal fits on the same axes.
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        subset = frame.loc[frame["x_index"] == x_index].sort_values("t_index")
        x_values = subset["mean_minus_tprime_gev2"]
        fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True)
        axes_flat = axes.ravel()

        for panel, parameter in enumerate(PHYSICS_PARAMETERS):
            ax = axes_flat[panel]
            ax.errorbar(
                x_values,
                subset[parameter],
                yerr=subset[f"{parameter}_stat"],
                marker="o",
                linestyle="-",
                capsize=2,
                label="Combined",
            )
            for period in PERIODS:
                ax.errorbar(
                    x_values,
                    subset[f"{parameter}_{period}"],
                    yerr=subset[f"{parameter}_stat_{period}"],
                    marker="o",
                    linestyle="none",
                    capsize=2,
                    label=PERIOD_LABELS[period],
                )
            # endfor
            ax.axhline(0.0, linewidth=0.8)
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            apply_parameter_y_limits(ax, parameter)
            ax.grid(alpha=0.25)
            if panel >= 6:
                ax.set_xlabel(r"Mean $-t^\prime$ (GeV$^2$)")
            # endif
        # endfor

        axes_flat[7].axis("off")
        axes_flat[8].axis("off")
        handles, labels = axes_flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.96, 0.06))
        fig.suptitle(
            rf"${x_low:.2f} \leq x_B < {x_high:.2f}$: period consistency",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
        path = output_dir / f"xB_bin_{x_index + 1}_by_period_v5.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths



def plot_period_consistency_heatmap(
    frame: pd.DataFrame,
    output_dir: Path,
) -> str:
    """Plot Delta NLL for each period and combined kinematic bin."""
    ensure_directory(output_dir)
    matrix = np.asarray(
        [
            frame[f"period_delta_nll_{period}"].to_numpy(dtype=float)
            for period in PERIODS
        ],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(15, 3.8))
    image = ax.imshow(matrix, aspect="auto", origin="upper")
    ax.set_xticks(np.arange(NUMBER_OF_BINS))
    ax.set_xticklabels(frame["bin_number"].astype(int).tolist())
    ax.set_yticks(np.arange(len(PERIODS)))
    ax.set_yticklabels([PERIOD_LABELS[period] for period in PERIODS])
    ax.set_xlabel("Combined kinematic-bin number")
    ax.set_ylabel("Run period")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(
        r"$\Delta\mathrm{NLL}_p="
        r"\mathrm{NLL}_p(\widehat{\theta}_{\mathrm{combined}})"
        r"-\mathrm{NLL}_{p,\min}$"
    )
    fig.tight_layout()
    path = output_dir / "period_consistency_delta_nll_v5.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)

def write_latex_table(
    frame: pd.DataFrame,
    path: Path,
) -> None:
    ensure_directory(path.parent)
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{rrrrrrrr}",
        r"\hline",
        (
            r"Bin & $F_{UU}^{\cos\phi}/F_{UU}$ & "
            r"$F_{UU}^{\cos2\phi}/F_{UU}$ & "
            r"$F_{LU}^{\sin\phi}/F_{UU}$ & "
            r"$F_{UL}^{\sin\phi}/F_{UU}$ & "
            r"$F_{UL}^{\sin2\phi}/F_{UU}$ & "
            r"$F_{LL}/F_{UU}$ & "
            r"$F_{LL}^{\cos\phi}/F_{UU}$ \\"
        ),
        r"\hline",
    ]

    for row in frame.itertuples(index=False):
        entries = [str(int(row.bin_number))]
        for parameter in PHYSICS_PARAMETERS:
            value = float(getattr(row, parameter))
            stat = float(getattr(row, f"{parameter}_stat"))
            axis = float(getattr(row, f"{parameter}_target_axis_sys"))
            entries.append(
                rf"${value:.5f}\pm{stat:.5f}\pm{axis:.5f}$"
            )
        # endfor
        lines.append(" & ".join(entries) + r" \\")
    # endfor

    lines.extend(
        [
            r"\hline",
            r"\end{tabular}",
            (
                r"\caption{Nominal simultaneous unbinned-likelihood results. "
                r"The first uncertainty is statistical and includes the "
                r"Gaussian-constrained dilution-factor statistical uncertainty. "
                r"The second is the target-axis treatment envelope. "
                r"Polarization and dilution-model scale uncertainties are not "
                r"included.}"
            ),
            r"\end{table}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# =============================================================================
# Command line
# =============================================================================

def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fit RGC exclusive-pi+ structure-function ratios with a simultaneous "
            "three-period unbinned conditional likelihood."
        )
    )
    parser.add_argument(
        "--tree",
        default=DEFAULT_TREE_NAME,
        help=f"ROOT tree name (default: {DEFAULT_TREE_NAME}).",
    )
    parser.add_argument(
        "--input",
        action="append",
        default=[],
        type=parse_input_override,
        metavar="PERIOD=FILE",
        help=(
            "Override one NH3 ROOT input. Repeat for multiple periods, e.g. "
            "--input su22=/path/file.root."
        ),
    )
    parser.add_argument(
        "--run-info-csv",
        type=Path,
        default=DEFAULT_RUN_INFO_CSV,
        help=f"Run-information CSV (default: {DEFAULT_RUN_INFO_CSV}).",
    )
    parser.add_argument(
        "--cut-json",
        type=Path,
        default=DEFAULT_CUT_JSON,
        help=f"Channel-selection cut JSON (default: {DEFAULT_CUT_JSON}).",
    )
    parser.add_argument(
        "--dilution-json",
        type=Path,
        default=None,
        help=(
            "Production dilution-factor JSON. By default the newest "
            "dilution_factors_production*.json is found under "
            f"{DEFAULT_DILUTION_DIR}."
        ),
    )
    parser.add_argument(
        "--dilution-dir",
        type=Path,
        default=DEFAULT_DILUTION_DIR,
        help=f"Dilution-factor output directory (default: {DEFAULT_DILUTION_DIR}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=DEFAULT_CACHE_PATH,
        help=f"Selected-event cache path (default: {DEFAULT_CACHE_PATH}).",
    )
    parser.add_argument(
        "--reuse-cache",
        action="store_true",
        help="Reuse an existing selected-event cache without rereading ROOT.",
    )
    parser.add_argument(
        "--chunk-size",
        default=DEFAULT_CHUNK_SIZE,
        help=f"uproot iterate chunk size (default: {DEFAULT_CHUNK_SIZE}).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=MAXIMUM_WORKERS,
        help=f"Worker processes, hard-limited to {MAXIMUM_WORKERS}.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip summary plot production.",
    )
    return parser


def main() -> int:
    parser = build_argument_parser()
    args = parser.parse_args()

    workers = max(
        1,
        min(
            int(args.workers),
            MAXIMUM_WORKERS,
            os.cpu_count() or 1,
            NUMBER_OF_BINS,
        ),
    )

    input_paths = {period: Path(path) for period, path in DEFAULT_INPUTS.items()}
    for period, path in args.input:
        input_paths[period] = path
    # endfor

    run_info_path = args.run_info_csv.expanduser().resolve()
    cut_json_path = args.cut_json.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    cache_path = args.cache.expanduser().resolve()

    dilution_json_path = (
        args.dilution_json.expanduser().resolve()
        if args.dilution_json is not None
        else find_default_dilution_json(
            args.dilution_dir.expanduser().resolve()
        ).resolve()
    )

    ensure_directory(output_dir)
    tables_dir = output_dir / "tables"
    json_dir = output_dir / "json"
    covariance_dir = output_dir / "covariance"
    plots_dir = output_dir / "plots"
    all_bins_plots_dir = plots_dir / "all_bins"
    aggregated_plots_dir = plots_dir / "aggregated"
    period_plots_dir = aggregated_plots_dir / "by_period"
    latex_dir = output_dir / "latex"
    for directory in (
        tables_dir,
        json_dir,
        covariance_dir,
        plots_dir,
        all_bins_plots_dir,
        aggregated_plots_dir,
        period_plots_dir,
        latex_dir,
    ):
        ensure_directory(directory)
    # endfor

    print("=" * 78)
    print("RGC exclusive enpi+ structure-function-ratio extraction")
    print("=" * 78)
    print(f"Working directory:    {Path.cwd()}")
    print(f"Tree:                 {args.tree}")
    print(f"Run information:      {run_info_path}")
    print(f"Channel cuts:         {cut_json_path}")
    print(f"Dilution factors:     {dilution_json_path}")
    print(f"Output directory:     {output_dir}")
    print(f"Selected-event cache: {cache_path}")
    print(f"Workers:              {workers} (hard maximum {MAXIMUM_WORKERS})")
    print("Fit policy:           one common PDF/parameter set across all three periods")
    print()

    run_records = parse_run_info_csv(run_info_path)
    run_states = run_state_arrays(run_records)
    cuts = load_nominal_cuts(cut_json_path)
    dilution_records = load_dilution_factors(dilution_json_path)

    cache_summary: dict[str, Any]
    if args.reuse_cache:
        events = load_event_cache(cache_path)
        cache_summary = {
            "cache_path": str(cache_path),
            "number_of_selected_events": int(events["runnum"].size),
            "reused": True,
        }
    else:
        cache_summary = build_event_cache(
            input_paths=input_paths,
            tree_name=args.tree,
            chunk_size=args.chunk_size,
            run_records=run_records,
            cuts=cuts,
            cache_path=cache_path,
        )
        cache_summary["reused"] = False
        events = load_event_cache(cache_path)
    # endif

    run_state_payload = {
        period: {
            key: values.tolist()
            for key, values in state.items()
        }
        for period, state in run_states.items()
    }
    dilution_payload = {
        period: {
            str(bin_number): {
                "x_index": record.x_index,
                "t_index": record.t_index,
                "value": record.value,
                "stat_uncertainty": record.stat_uncertainty,
            }
            for (record_period, bin_number), record in dilution_records.items()
            if record_period == period
        }
        for period in PERIODS
    }

    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=initialize_fit_worker,
        initargs=(
            str(cache_path),
            run_state_payload,
            dilution_payload,
        ),
    ) as executor:
        futures = {
            executor.submit(fit_bin_worker, bin_number): bin_number
            for bin_number in range(1, NUMBER_OF_BINS + 1)
        }
        for future in as_completed(futures):
            bin_number = futures[future]
            result = future.result()
            results.append(result)
            nominal = result["variants"]["nominal"]
            period_validity = ",".join(
                f"{PERIOD_LABELS[period]}="
                f"{result['period_fits'][period]['valid']}"
                for period in PERIODS
            )
            print(
                f"[bin {bin_number:02d}] "
                f"N={nominal['metadata']['number_of_events']:,}; "
                f"combined_valid={nominal['valid']}; "
                f"NLL={nominal['minimum_nll']:.6f}; "
                f"EDM={nominal['edm']:.3e}; "
                f"period_fits({period_validity})"
            )
        # endfor
    # endwith

    results.sort(key=lambda item: item["bin_number"])
    frame = flatten_fit_results(results)
    csv_path = tables_dir / "structure_function_ratios_v5.csv"
    frame.to_csv(csv_path, index=False)

    detailed_json_path = json_dir / "structure_function_ratios_v5.json"
    write_json(
        detailed_json_path,
        {
            "schema_version": 5,
            "analysis": "RGC exclusive enpi+ structure-function-ratio extraction",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "fit_policy": (
                "One common set of structure-function ratios is fitted "
                "simultaneously to Su22, Fa22, and Sp23 in each kinematic bin."
            ),
            "target_axis_systematic_definition": (
                "For each fitted longitudinal observable, the systematic is the "
                "largest absolute displacement from nominal among no_projection, "
                "longitudinal_scaled_transverse_plus, and "
                "longitudinal_scaled_transverse_minus."
            ),
            "transverse_stress_model_warning": (
                "The two transverse stress variants use the signed nominal "
                "longitudinal amplitudes as scales for corresponding transverse "
                "harmonics. One preserves all mapped signs and one reverses them. "
                "They are loose data-driven systematic estimates, not physical "
                "equality assumptions."
            ),
            "polarization_uncertainty_policy": (
                "Beam- and target-polarization uncertainties are not included. "
                "They will be imposed later as correlated scale systematics."
            ),
            "dilution_uncertainty_policy": (
                "Only the recommended dilution-factor bootstrap statistical "
                "uncertainties are included through Gaussian-constrained "
                "period-specific nuisance parameters. The correlated dilution "
                "model scale is not included."
            ),
            "beam_polarization": BEAM_POLARIZATION,
            "beam_energy_gev": BEAM_ENERGY_GEV,
            "cache": cache_summary,
            "inputs": {
                "root_files": {
                    period: str(input_paths[period].expanduser().resolve())
                    for period in PERIODS
                },
                "run_info_csv": str(run_info_path),
                "run_info_csv_sha256": sha256_file(run_info_path),
                "channel_cut_json": str(cut_json_path),
                "channel_cut_json_sha256": sha256_file(cut_json_path),
                "dilution_json": str(dilution_json_path),
                "dilution_json_sha256": sha256_file(dilution_json_path),
            },
            "results": results,
        },
    )

    for result in results:
        bin_number = result["bin_number"]
        nominal = result["variants"]["nominal"]
        if nominal["covariance"] is None:
            continue
        # endif
        np.save(
            covariance_dir / f"bin_{bin_number:02d}_nominal_covariance_v5.npy",
            np.asarray(nominal["covariance"], dtype=np.float64),
        )
        np.save(
            covariance_dir / f"bin_{bin_number:02d}_nominal_correlation_v5.npy",
            np.asarray(nominal["correlation"], dtype=np.float64),
        )
    # endfor

    latex_path = latex_dir / "structure_function_ratios_v5.tex"
    write_latex_table(frame, latex_path)

    plot_paths: dict[str, list[str]] = {
        "all_bins": [],
        "aggregated": [],
        "aggregated_by_period": [],
    }
    if not args.skip_plots:
        plot_paths["all_bins"] = plot_parameter_summaries(
            frame, all_bins_plots_dir
        )
        plot_paths["aggregated"] = plot_aggregated_by_x(
            frame, aggregated_plots_dir
        )
        plot_paths["aggregated_by_period"] = plot_aggregated_by_period(
            frame, period_plots_dir
        )
        plot_paths["aggregated_by_period"].append(
            plot_period_consistency_heatmap(frame, period_plots_dir)
        )
    # endif

    manifest_path = output_dir / "asymmetry_extraction_manifest_v5.json"
    write_json(
        manifest_path,
        {
            "schema_version": 5,
            "script": str(Path(__file__).resolve()),
            "workers": workers,
            "products": {
                "csv": str(csv_path),
                "detailed_json": str(detailed_json_path),
                "latex": str(latex_path),
                "covariance_directory": str(covariance_dir),
                "plots": plot_paths,
                "cache": str(cache_path),
            },
        },
    )

    invalid_bins = [
        int(row.bin_number)
        for row in frame.itertuples(index=False)
        if not bool(row.nominal_fit_valid)
    ]

    print()
    print("Asymmetry extraction complete.")
    print(f"  Selected events: {events['runnum'].size:,}")
    print(f"  CSV:             {csv_path}")
    print(f"  Detailed JSON:   {detailed_json_path}")
    print(f"  LaTeX:           {latex_path}")
    print(f"  Manifest:        {manifest_path}")
    if invalid_bins:
        print(f"  WARNING: invalid nominal fits in bins {invalid_bins}")
        return 2
    # endif
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
