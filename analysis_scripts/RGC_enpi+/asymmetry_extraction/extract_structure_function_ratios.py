#!/usr/bin/env python3
"""
extract_structure_function_ratios.py

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
Three fits are performed in every kinematic bin.

  nominal
    P_L = P_t cos(theta_gamma), P_T = 0.

  no_projection
    P_L = P_t, P_T = 0.

  external_data_informed
    P_L = P_t cos(theta_gamma), P_T = P_t sin(theta_gamma).  Four fixed
    transverse amplitudes are supplied from external HERMES measurements:

      F_UT^{sin(phi)} / F_UU        = -0.117
      F_UT^{sin(2phi)} / F_UU       = +0.109
      F_LT^{cos(0 phi)} / F_UU      = -0.150
      F_LT^{cos(phi)} / F_UU        = -0.150

    The two UT inputs are inverse-total-variance weighted averages of measured
    HERMES exclusive-pi+ amplitudes over 0.07 < x_B < 0.35, omitting the
    lowest-x_B bin. Statistical and systematic uncertainties are combined in
    quadrature for the weights. No factor of two is applied.

    The two LT inputs are rounded estimates from the open highest-z HERMES
    SIDIS pi+ points in 0.84 < z < 1.20. Those points are closest to the
    exclusive z -> 1 limit and are not included in the ordinary x and P_hT
    projections. They are documented as digitized/rounded estimates rather
    than exact tabulated values. No factor of two is applied.

Only the directly corresponding harmonics are populated. There is no
transverse sin(3phi) or cos(2phi) input. The external-data-informed fit is a
controlled leakage study, not a claim that the SIDIS amplitudes equal the
exclusive amplitudes at CLAS12 kinematics.

For each observable a, the target-axis systematic is

    max(
        abs(a_no_projection - a_nominal),
        abs(a_external_data_informed - a_nominal)
    ).


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

PERIOD_COLORS: dict[str, str] = {
    "su22": "tab:orange",
    "fa22": "tab:green",
    "sp23": "tab:red",
}
COMBINED_COLOR = "tab:blue"

VARIANT_LABELS: dict[str, str] = {
    "nominal": r"Nominal: $P_L=P_t\cos\theta_\gamma$",
    "no_projection": r"No projection: $P_L=P_t$",
    "external_data_informed": r"External-data-informed transverse leakage",
}
VARIANT_COLORS: dict[str, str] = {
    "nominal": "tab:blue",
    "no_projection": "tab:purple",
    "external_data_informed": "tab:orange",
}

# External transverse amplitudes used in the leakage fit.
#
# Exclusive UT reference:
#   A. Airapetian et al. (HERMES),
#   "Single-spin azimuthal asymmetry in exclusive electroproduction of pi+
#   mesons on transversely polarized protons,"
#   Phys. Lett. B 682 (2010) 345-350,
#   DOI: 10.1016/j.physletb.2009.11.039,
#   arXiv:0907.2596.
#   Official HERMES ASCII data file: ivana.AUTexclpi-1.dat.
#
# SIDIS LT reference:
#   A. Airapetian et al. (HERMES),
#   "Azimuthal single- and double-spin asymmetries in semi-inclusive
#   deep-inelastic lepton scattering by transversely polarized protons,"
#   JHEP 12 (2020) 010,
#   DOI: 10.1007/JHEP12(2020)010,
#   arXiv:2007.07755.
#
# The UT values below are direct inverse-total-variance weighted averages of
# exact tabulated exclusive-pi+ measurements. The LT values are rounded
# estimates from the open highest-z pi+ points (0.84 < z < 1.20), which are not
# included in the ordinary x and P_hT projections. No factor of two is applied
# to any of the four inputs.
EXTERNAL_TRANSVERSE_INPUTS: dict[str, dict[str, Any]] = {
    "ut_sin_phi": {
        "value": -0.117,
        "measured_central_value": -0.11729214289950972,
        "scale_factor": 1.0,
        "observable": "A_UT,l^{sin(phi-phi_S)}",
        "channel": "exclusive pi+",
        "kinematic_selection": {
            "averaging_policy": (
                "Inverse-total-variance weighted average over the HERMES "
                "x_B bins 0.07-0.10, 0.10-0.15, and 0.15-0.35; the "
                "0.03-0.07 bin is omitted. Statistical and systematic "
                "uncertainties are combined in quadrature."
            ),
            "xB_range_used": [0.07, 0.35],
            "weighted_mean_xB": 0.13848440446760085,
            "weighted_mean_uncertainty": 0.09362561740330531,
            "input_rows": [
                {
                    "xB_bin": [0.07, 0.10],
                    "mean_xB": 0.09,
                    "asymmetry": -0.071,
                    "stat_uncertainty": 0.180,
                    "syst_uncertainty": 0.022,
                },
                {
                    "xB_bin": [0.10, 0.15],
                    "mean_xB": 0.13,
                    "asymmetry": -0.093,
                    "stat_uncertainty": 0.130,
                    "syst_uncertainty": 0.029,
                },
                {
                    "xB_bin": [0.15, 0.35],
                    "mean_xB": 0.21,
                    "asymmetry": -0.219,
                    "stat_uncertainty": 0.180,
                    "syst_uncertainty": 0.065,
                },
            ],
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": (
                "A. Airapetian et al., Phys. Lett. B 682 (2010) 345-350"
            ),
            "title": (
                "Single-spin azimuthal asymmetry in exclusive "
                "electroproduction of pi+ mesons on transversely "
                "polarized protons"
            ),
            "doi": "10.1016/j.physletb.2009.11.039",
            "arxiv": "0907.2596",
            "data_source": "Official HERMES ASCII file ivana.AUTexclpi-1.dat",
            "data_note": (
                "Direct weighted average of exact tabulated exclusive-pi+ "
                "measurements; no doubling factor."
            ),
        },
    },
    "ut_sin_2phi": {
        "value": 0.109,
        "measured_central_value": 0.10930777687298002,
        "scale_factor": 1.0,
        "observable": "A_UT,l^{sin(2phi-phi_S)}",
        "channel": "exclusive pi+",
        "kinematic_selection": {
            "averaging_policy": (
                "Inverse-total-variance weighted average over the HERMES "
                "x_B bins 0.07-0.10, 0.10-0.15, and 0.15-0.35; the "
                "0.03-0.07 bin is omitted. Statistical and systematic "
                "uncertainties are combined in quadrature."
            ),
            "xB_range_used": [0.07, 0.35],
            "weighted_mean_xB": 0.13896698591627263,
            "weighted_mean_uncertainty": 0.09952058355881978,
            "input_rows": [
                {
                    "xB_bin": [0.07, 0.10],
                    "mean_xB": 0.09,
                    "asymmetry": -0.196,
                    "stat_uncertainty": 0.181,
                    "syst_uncertainty": 0.106,
                },
                {
                    "xB_bin": [0.10, 0.15],
                    "mean_xB": 0.13,
                    "asymmetry": 0.178,
                    "stat_uncertainty": 0.132,
                    "syst_uncertainty": 0.024,
                },
                {
                    "xB_bin": [0.15, 0.35],
                    "mean_xB": 0.21,
                    "asymmetry": 0.247,
                    "stat_uncertainty": 0.192,
                    "syst_uncertainty": 0.085,
                },
            ],
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": (
                "A. Airapetian et al., Phys. Lett. B 682 (2010) 345-350"
            ),
            "title": (
                "Single-spin azimuthal asymmetry in exclusive "
                "electroproduction of pi+ mesons on transversely "
                "polarized protons"
            ),
            "doi": "10.1016/j.physletb.2009.11.039",
            "arxiv": "0907.2596",
            "data_source": "Official HERMES ASCII file ivana.AUTexclpi-1.dat",
            "data_note": (
                "Direct weighted average of exact tabulated exclusive-pi+ "
                "measurements; no doubling factor."
            ),
        },
    },
    "lt_cos_0phi": {
        "value": -0.150,
        "measured_central_value": -0.150,
        "scale_factor": 1.0,
        "observable": "2<cos(phi_S)>/sqrt(2 epsilon (1-epsilon))",
        "channel": "SIDIS pi+",
        "kinematic_selection": {
            "z_bin": [0.84, 1.20],
            "x_acceptance": [0.023, 0.600],
            "selection_note": (
                "Open highest-z point in the one-dimensional pi+ z "
                "projection. This point is not included in the ordinary "
                "x and P_hT projections and is used as the closest "
                "available SIDIS proxy for z approximately 1."
            ),
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": "A. Airapetian et al., JHEP 12 (2020) 010",
            "title": (
                "Azimuthal single- and double-spin asymmetries in "
                "semi-inclusive deep-inelastic lepton scattering by "
                "transversely polarized protons"
            ),
            "doi": "10.1007/JHEP12(2020)010",
            "arxiv": "2007.07755",
            "figure": 31,
            "data_note": (
                "Digitized/rounded estimate from the open highest-z pi+ "
                "point; not an exact tabulated value; no doubling factor."
            ),
        },
    },
    "lt_cos_phi": {
        "value": -0.150,
        "measured_central_value": -0.150,
        "scale_factor": 1.0,
        "observable": "2<cos(phi-phi_S)>/sqrt(1-epsilon^2)",
        "channel": "SIDIS pi+",
        "kinematic_selection": {
            "z_bin": [0.84, 1.20],
            "x_acceptance": [0.023, 0.600],
            "selection_note": (
                "Open highest-z point in the one-dimensional pi+ z "
                "projection. This point is not included in the ordinary "
                "x and P_hT projections and is used as the closest "
                "available SIDIS proxy for z approximately 1."
            ),
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": "A. Airapetian et al., JHEP 12 (2020) 010",
            "title": (
                "Azimuthal single- and double-spin asymmetries in "
                "semi-inclusive deep-inelastic lepton scattering by "
                "transversely polarized protons"
            ),
            "doi": "10.1007/JHEP12(2020)010",
            "arxiv": "2007.07755",
            "figure": 21,
            "data_note": (
                "Digitized/rounded estimate from the open highest-z pi+ "
                "point; not an exact tabulated value; no doubling factor."
            ),
        },
    },
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
DEFAULT_CACHE_PATH = DEFAULT_OUTPUT_DIR / "cache/selected_events.npz"

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

DEFAULT_ISR_INPUTS: dict[str, Path] = {
    period: PAPER_VERSIONS_DIR
    / f"rgc_{period}_inb_NH3_epi+_ISR_mom_corrections.root"
    for period in PERIODS
}

DEFAULT_CHANNEL_SELECTION_MANIFEST = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "channel_selection_manifest.json"
)
DEFAULT_ISR_CUT_JSON = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "isr/final_carbon_assisted_cuts/tables/"
    "final_carbon_assisted_mx2_cuts.json"
)
DEFAULT_ISR_DILUTION_JSON = Path(
    "../dilution_factor/output/dilution_factor_determination/"
    "isr/dilution_factors_production.json"
)

FIT_VARIANTS: tuple[str, ...] = (
    "nominal",
    "no_projection",
    "external_data_informed",
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
    "u1": (-1.0, 1.0),
    "u2": (-1.0, 1.0),
    "lu1": (-0.5, 0.5),
    "ul1": (-0.5, 0.5),
    "ul2": (-0.5, 0.5),
    "ll0": (-1.0, 1.0),
    "ll1": (-1.0, 1.0),
}

# Physics-motivated panel layout:
#   top row:    unpolarized cosine modulations
#   middle row: beam- and target-single-spin ratios
#   bottom row: double-spin ratios
AGGREGATED_PANEL_ORDER: tuple[str | None, ...] = (
    "u1", "u2", None,
    "lu1", "ul1", "ul2",
    "ll0", "ll1", None,
)

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

def external_data_informed_transverse_terms(
    phi: np.ndarray,
    r_b: np.ndarray,
    r_c: np.ndarray,
    r_v: np.ndarray,
    r_w: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the fixed external-data-informed transverse leakage terms.

    The restricted harmonic mapping is

      F_UT^{sin(phi)} / F_UU
      F_UT^{sin(2phi)} / F_UU
      F_LT^{cos(0 phi)} / F_UU
      F_LT^{cos(phi)} / F_UU

    with no sin(3phi) or cos(2phi) transverse input.

    The HERMES A_UT,l^{sin(phi-phi_S)} term carries no extra
    depolarization factor in this convention. The
    A_UT,l^{sin(2phi-phi_S)} term carries r_v.
    """
    ut_sin_phi = float(EXTERNAL_TRANSVERSE_INPUTS["ut_sin_phi"]["value"])
    ut_sin_2phi = float(
        EXTERNAL_TRANSVERSE_INPUTS["ut_sin_2phi"]["value"]
    )
    lt_cos_0phi = float(EXTERNAL_TRANSVERSE_INPUTS["lt_cos_0phi"]["value"])
    lt_cos_phi = float(EXTERNAL_TRANSVERSE_INPUTS["lt_cos_phi"]["value"])

    ut = (
        ut_sin_phi * np.sin(phi)
        + r_v * ut_sin_2phi * np.sin(2.0 * phi)
    )
    lt = (
        r_w * lt_cos_0phi
        + r_c * lt_cos_phi * np.cos(phi)
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

    if variant == "external_data_informed":
        ut_fixed, lt_fixed = external_data_informed_transverse_terms(
            phi,
            r_b,
            r_c,
            r_v,
            r_w,
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

            if variant == "external_data_informed":
                ut_sin_phi = float(
                    EXTERNAL_TRANSVERSE_INPUTS["ut_sin_phi"]["value"]
                )
                ut_sin_2phi = float(
                    EXTERNAL_TRANSVERSE_INPUTS["ut_sin_2phi"]["value"]
                )
                lt_cos_0phi = float(
                    EXTERNAL_TRANSVERSE_INPUTS["lt_cos_0phi"]["value"]
                )
                lt_cos_phi = float(
                    EXTERNAL_TRANSVERSE_INPUTS["lt_cos_phi"]["value"]
                )
                ut_fixed = (
                    ut_sin_phi * data["sin_phi"]
                    + event_r_v * ut_sin_2phi * data["sin_2phi"]
                )
                lt_fixed = (
                    event_r_w * lt_cos_0phi
                    + event_r_c * lt_cos_phi * data["cos_phi"]
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
    fixed_physics_parameters: Mapping[str, float] | None = None,
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
    if fixed_physics_parameters is not None:
        for name, value in fixed_physics_parameters.items():
            if name not in PHYSICS_PARAMETERS:
                raise RuntimeError(
                    f"Unknown fixed physics parameter {name!r}."
                )
            # endif
            initial[name] = float(value)
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
            if (
                fixed_physics_parameters is not None
                and parameter_name in fixed_physics_parameters
            ):
                candidate.values[parameter_name] = float(
                    fixed_physics_parameters[parameter_name]
                )
                candidate.fixed[parameter_name] = True
            # endif
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

    sign_reversed = dict(initial)
    for name in PHYSICS_PARAMETERS:
        if (
            fixed_physics_parameters is None
            or name not in fixed_physics_parameters
        ):
            sign_reversed[name] = -float(initial[name])
        # endif
    # endfor
    start_candidates.append(sign_reversed)

    bounded_offset = dict(initial)
    offset_pattern = {
        "u1": 0.20,
        "u2": -0.20,
        "lu1": 0.05,
        "ul1": -0.05,
        "ul2": 0.05,
        "ll0": 0.20,
        "ll1": -0.20,
    }
    for name, offset in offset_pattern.items():
        if (
            fixed_physics_parameters is None
            or name not in fixed_physics_parameters
        ):
            low, high = PARAMETER_LIMITS[name]
            bounded_offset[name] = float(
                np.clip(
                    float(initial[name]) + offset,
                    low + 1.0e-6,
                    high - 1.0e-6,
                )
            )
        # endif
    # endfor
    start_candidates.append(bounded_offset)

    attempted: list[Minuit] = []
    for start_values in start_candidates:
        candidate = configured_minuit(start_values)
        candidate.migrad(ncall=50000)
        if not candidate.fmin.is_valid:
            candidate.simplex(ncall=50000)
            candidate.strategy = 2
            candidate.migrad(ncall=120000)
        # endif
        candidate.hesse()
        attempted.append(candidate)
    # endfor

    def candidate_quality(candidate: Minuit) -> tuple[int, float, float]:
        fmin = candidate.fmin
        trustworthy = (
            fmin.is_valid
            and fmin.has_accurate_covar
            and fmin.has_posdef_covar
            and not fmin.has_parameters_at_limit
            and math.isfinite(float(candidate.fval))
        )
        valid = fmin.is_valid and math.isfinite(float(candidate.fval))
        rank = 0 if trustworthy else (1 if valid else 2)
        edm = (
            float(fmin.edm)
            if math.isfinite(float(fmin.edm))
            else math.inf
        )
        return rank, float(candidate.fval), edm

    minuit = min(attempted, key=candidate_quality)

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
            "fixed_physics_parameters": (
                dict(fixed_physics_parameters)
                if fixed_physics_parameters is not None else None
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


def fit_bin_worker(
    bin_number: int,
    include_target_axis_study: bool = True,
) -> dict[str, Any]:
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
    variants = {"nominal": nominal}
    projection_systematic: dict[str, float | None] = {
        parameter: None for parameter in PHYSICS_PARAMETERS
    }
    external_data_systematic: dict[str, float | None] = {
        parameter: None for parameter in PHYSICS_PARAMETERS
    }
    full_three_fit_spread: dict[str, float | None] = {
        parameter: None for parameter in PHYSICS_PARAMETERS
    }
    systematics: dict[str, float | None] = {
        parameter: None for parameter in PHYSICS_PARAMETERS
    }

    if include_target_axis_study:
        no_projection = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "no_projection",
            initial_values=nominal["values"],
        )
        external_data_informed = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "external_data_informed",
            initial_values=nominal["values"],
        )
        variants.update({
            "no_projection": no_projection,
            "external_data_informed": external_data_informed,
        })
        for parameter in PHYSICS_PARAMETERS:
            nominal_value = nominal["values"][parameter]
            no_projection_value = no_projection["values"][parameter]
            external_value = external_data_informed["values"][parameter]
            projection_systematic[parameter] = abs(
                no_projection_value - nominal_value
            )
            external_data_systematic[parameter] = abs(
                external_value - nominal_value
            )
            systematics[parameter] = max(
                projection_systematic[parameter],
                external_data_systematic[parameter],
            )
            full_three_fit_spread[parameter] = (
                max(nominal_value, no_projection_value, external_value)
                - min(nominal_value, no_projection_value, external_value)
            )
        # endfor
    # endif

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

    period_constraint_fits: dict[
        str,
        dict[str, dict[str, Any]],
    ] = {
        "fix_u1": {},
        "fix_u2": {},
        "fix_u1_u2": {},
    }
    for period in PERIODS:
        period_constraint_fits["fix_u1"][period] = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
            initial_values=period_fits[period]["values"],
            fixed_physics_parameters={
                "u1": nominal["values"]["u1"],
            },
        )
        period_constraint_fits["fix_u2"][period] = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
            initial_values=period_fits[period]["values"],
            fixed_physics_parameters={
                "u2": nominal["values"]["u2"],
            },
        )
        period_constraint_fits["fix_u1_u2"][period] = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
            initial_values=period_fits[period]["values"],
            fixed_physics_parameters={
                "u1": nominal["values"]["u1"],
                "u2": nominal["values"]["u2"],
            },
        )
    # endfor

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
        "period_constraint_fits": period_constraint_fits,
        "period_consistency": period_consistency,
        "target_axis_systematic": systematics,
        "projection_systematic": projection_systematic,
        "external_data_systematic": external_data_systematic,
        "full_three_fit_spread": full_three_fit_spread,
        "external_transverse_inputs": (
            EXTERNAL_TRANSVERSE_INPUTS if include_target_axis_study else None
        ),
        "target_axis_study_performed": include_target_axis_study,
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
            row[f"period_fit_accurate_covariance_{period}"] = period_fit[
                "accurate_covariance"
            ]
            row[f"period_fit_positive_definite_covariance_{period}"] = (
                period_fit["positive_definite_covariance"]
            )
            row[f"period_fit_parameters_at_limit_{period}"] = period_fit[
                "parameters_at_limit"
            ]
            row[f"period_fit_nll_{period}"] = period_fit["minimum_nll"]
            row[f"period_fit_edm_{period}"] = period_fit["edm"]
            row[f"period_delta_nll_{period}"] = result[
                "period_consistency"
            ][period]["delta_nll"]

            for constraint_name in ("fix_u1", "fix_u2", "fix_u1_u2"):
                constrained_fit = result["period_constraint_fits"][
                    constraint_name
                ][period]
                prefix = f"{constraint_name}_period_fit"
                row[f"{prefix}_valid_{period}"] = constrained_fit["valid"]
                row[f"{prefix}_accurate_covariance_{period}"] = (
                    constrained_fit["accurate_covariance"]
                )
                row[
                    f"{prefix}_positive_definite_covariance_{period}"
                ] = constrained_fit["positive_definite_covariance"]
                row[f"{prefix}_parameters_at_limit_{period}"] = (
                    constrained_fit["parameters_at_limit"]
                )
                row[f"{prefix}_edm_{period}"] = constrained_fit["edm"]
            # endfor
        # endfor

        ll0_value = float(nominal["values"]["ll0"])
        ll0_error = float(nominal["errors"]["ll0"])
        row["ll0_central_outside_unit_interval"] = abs(ll0_value) > 1.0
        row["ll0_stat_interval_crosses_unit_boundary"] = (
            ll0_value - ll0_error < -1.0
            or ll0_value + ll0_error > 1.0
        )
        row["ll0_distance_to_nearest_unit_boundary"] = (
            1.0 - abs(ll0_value)
        )

        for parameter in PHYSICS_PARAMETERS:
            row[parameter] = nominal["values"][parameter]
            row[f"{parameter}_stat"] = nominal["errors"][parameter]
            row[f"{parameter}_projection_sys"] = result[
                "projection_systematic"
            ][parameter]
            row[f"{parameter}_external_data_sys"] = result[
                "external_data_systematic"
            ][parameter]
            row[f"{parameter}_three_fit_full_spread"] = result[
                "full_three_fit_spread"
            ][parameter]
            row[f"{parameter}_target_axis_sys"] = result[
                "target_axis_systematic"
            ][parameter]
            for variant in FIT_VARIANTS:
                variant_fit = result["variants"].get(variant)
                row[f"{parameter}_{variant}"] = (
                    variant_fit["values"][parameter]
                    if variant_fit is not None
                    else np.nan
                )
            # endfor
            for period in PERIODS:
                period_fit = result["period_fits"][period]
                row[f"{parameter}_{period}"] = period_fit[
                    "values"
                ][parameter]
                row[f"{parameter}_stat_{period}"] = period_fit[
                    "errors"
                ][parameter]

                for constraint_name in (
                    "fix_u1",
                    "fix_u2",
                    "fix_u1_u2",
                ):
                    constrained_fit = result["period_constraint_fits"][
                        constraint_name
                    ][period]
                    row[
                        f"{parameter}_{constraint_name}_{period}"
                    ] = constrained_fit["values"][parameter]
                    row[
                        f"{parameter}_{constraint_name}_stat_{period}"
                    ] = constrained_fit["errors"][parameter]
                # endfor
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


def period_fit_quality_mask(
    frame: pd.DataFrame,
    period: str,
) -> np.ndarray:
    """Return True only for trustworthy period-only diagnostic fits."""
    return (
        frame[f"period_fit_valid_{period}"].astype(bool).to_numpy()
        & frame[
            f"period_fit_accurate_covariance_{period}"
        ].astype(bool).to_numpy()
        & frame[
            f"period_fit_positive_definite_covariance_{period}"
        ].astype(bool).to_numpy()
        & ~frame[
            f"period_fit_parameters_at_limit_{period}"
        ].astype(bool).to_numpy()
    )


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
            label="Target-axis systematic",
        )
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("Combined kinematic-bin number")
        ax.set_ylabel(PARAMETER_LABELS[parameter])
        ax.set_xticks(bins)
        apply_parameter_y_limits(ax, parameter)
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        path = output_dir / f"{parameter}_summary.png"
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
    Write one physics-ordered 3x3 canvas per xB bin.

    Top row: unpolarized cosine modulations.
    Middle row: LU, UL sin(phi), and UL sin(2phi).
    Bottom row: LL and LL cos(phi).
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        subset = frame.loc[frame["x_index"] == x_index].sort_values("t_index")
        fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True)
        axes_flat = axes.ravel()

        for panel, parameter in enumerate(AGGREGATED_PANEL_ORDER):
            ax = axes_flat[panel]
            if parameter is None:
                ax.axis("off")
                continue
            # endif

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
                label="Target-axis systematic",
            )
            ax.axhline(0.0, linewidth=0.8)
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            apply_parameter_y_limits(ax, parameter)
            ax.grid(alpha=0.25)
            if panel >= 6:
                ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
            # endif
        # endfor

        handles, labels = axes_flat[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="lower right",
            bbox_to_anchor=(0.96, 0.06),
        )
        fig.suptitle(
            rf"${x_low:.2f} \leq x_B < {x_high:.2f}$",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
        path = output_dir / f"xB_bin_{x_index + 1}_combined.png"
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
    Write physics-ordered 3x3 canvases comparing the simultaneous result with
    independent Su22, Fa22, and Sp23 diagnostic fits.

    Invalid or covariance-defective period fits are not drawn as ordinary
    error bars.  Their central values are marked with an x so they cannot be
    mistaken for precise measurements with zero or misleading Hesse errors.
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        subset = frame.loc[frame["x_index"] == x_index].sort_values("t_index")
        x_values = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
        fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True)
        axes_flat = axes.ravel()

        for panel, parameter in enumerate(AGGREGATED_PANEL_ORDER):
            ax = axes_flat[panel]
            if parameter is None:
                ax.axis("off")
                continue
            # endif

            ax.errorbar(
                x_values,
                subset[parameter],
                yerr=subset[f"{parameter}_stat"],
                marker="o",
                linestyle="-",
                capsize=2,
                color=COMBINED_COLOR,
                label="Simultaneous",
            )

            for period in PERIODS:
                quality = period_fit_quality_mask(subset, period)
                values = subset[f"{parameter}_{period}"].to_numpy(dtype=float)
                errors = subset[
                    f"{parameter}_stat_{period}"
                ].to_numpy(dtype=float)

                if np.any(quality):
                    ax.errorbar(
                        x_values[quality],
                        values[quality],
                        yerr=errors[quality],
                        marker="o",
                        linestyle="none",
                        capsize=2,
                        color=PERIOD_COLORS[period],
                        label=PERIOD_LABELS[period],
                    )
                # endif

                invalid = ~quality
                if np.any(invalid):
                    ax.plot(
                        x_values[invalid],
                        values[invalid],
                        marker="x",
                        linestyle="none",
                        markersize=8,
                        markeredgewidth=1.8,
                        color=PERIOD_COLORS[period],
                        label=f"{PERIOD_LABELS[period]} invalid",
                    )
                # endif
            # endfor

            ax.axhline(0.0, linewidth=0.8)
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            apply_parameter_y_limits(ax, parameter)
            ax.grid(alpha=0.25)
            if panel >= 6:
                ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
            # endif
        # endfor

        # Deduplicate legend entries produced in every panel.
        handles: list[Any] = []
        labels: list[str] = []
        for ax in axes_flat:
            panel_handles, panel_labels = ax.get_legend_handles_labels()
            for handle, label in zip(panel_handles, panel_labels):
                if label not in labels:
                    handles.append(handle)
                    labels.append(label)
                # endif
            # endfor
        # endfor
        fig.legend(
            handles,
            labels,
            loc="lower right",
            bbox_to_anchor=(0.96, 0.06),
        )
        fig.suptitle(
            rf"${x_low:.2f} \leq x_B < {x_high:.2f}$: period diagnostics",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
        path = output_dir / f"xB_bin_{x_index + 1}_by_period.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths


def plot_target_axis_variants(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """Compare nominal, no-projection, and external-data-informed T fits."""
    ensure_directory(output_dir)
    paths: list[str] = []
    variants = tuple(VARIANT_LABELS)

    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        subset = frame.loc[frame["x_index"] == x_index].sort_values("t_index")
        x_values = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
        fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True)
        axes_flat = axes.ravel()

        for panel, parameter in enumerate(AGGREGATED_PANEL_ORDER):
            ax = axes_flat[panel]
            if parameter is None:
                ax.axis("off")
                continue
            # endif

            for variant in variants:
                ax.errorbar(
                    x_values,
                    subset[f"{parameter}_{variant}"],
                    yerr=(
                        subset[f"{parameter}_stat"]
                        if variant == "nominal"
                        else None
                    ),
                    marker="o",
                    linestyle="-",
                    linewidth=1.0,
                    capsize=2,
                    color=VARIANT_COLORS[variant],
                    label=VARIANT_LABELS[variant],
                )
            # endfor

            # Draw the quoted target-axis systematic as narrow, non-overlapping
            # gray rectangles centered on the nominal points.
            x_array = np.asarray(x_values, dtype=np.float64)
            y_array = subset[parameter].to_numpy(dtype=np.float64)
            sys_array = subset[
                f"{parameter}_target_axis_sys"
            ].to_numpy(dtype=np.float64)
            if x_array.size > 1:
                separations = np.diff(np.sort(x_array))
                half_width = 0.10 * float(np.min(separations))
            else:
                half_width = 0.015
            # endif
            for point_index, (x_point, y_point, y_sys) in enumerate(
                zip(x_array, y_array, sys_array)
            ):
                ax.fill_between(
                    [x_point - half_width, x_point + half_width],
                    [y_point - y_sys, y_point - y_sys],
                    [y_point + y_sys, y_point + y_sys],
                    color="0.65",
                    alpha=0.45,
                    linewidth=0.0,
                    label=(
                        "Quoted target-axis systematic"
                        if point_index == 0 else None
                    ),
                    zorder=1,
                )
            # endfor

            ax.axhline(0.0, linewidth=0.8)
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            apply_parameter_y_limits(ax, parameter)
            ax.grid(alpha=0.25)
            if panel >= 6:
                ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
            # endif
        # endfor

        handles, labels = axes_flat[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="lower right",
            bbox_to_anchor=(0.97, 0.055),
        )
        fig.suptitle(
            rf"${x_low:.2f} \leq x_B < {x_high:.2f}$: target-axis variants",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
        path = output_dir / f"xB_bin_{x_index + 1}_target_axis_variants.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths


def plot_period_stability(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """
    Compare fully free period fits with diagnostic fits fixing u1 only,
    u2 only, or both u1 and u2 to the simultaneous values.
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    constraint_specs = (
        ("free", "Fully free period fit", "o", None),
        ("fix_u1", r"$u_1$ fixed", "s", "tab:purple"),
        ("fix_u2", r"$u_2$ fixed", "^", "tab:brown"),
        ("fix_u1_u2", r"$u_1,u_2$ fixed", "D", COMBINED_COLOR),
    )

    for period in PERIODS:
        for x_index, (x_low, x_high) in enumerate(XB_BINS):
            subset = frame.loc[
                frame["x_index"] == x_index
            ].sort_values("t_index")
            x_values = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
            fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True)
            axes_flat = axes.ravel()

            for panel, parameter in enumerate(AGGREGATED_PANEL_ORDER):
                ax = axes_flat[panel]
                if parameter is None:
                    ax.axis("off")
                    continue
                # endif

                for constraint_name, label, marker, color in constraint_specs:
                    if constraint_name == "free":
                        values = subset[f"{parameter}_{period}"]
                        errors = subset[f"{parameter}_stat_{period}"]
                        draw_color = PERIOD_COLORS[period]
                    else:
                        values = subset[
                            f"{parameter}_{constraint_name}_{period}"
                        ]
                        errors = subset[
                            f"{parameter}_{constraint_name}_stat_{period}"
                        ]
                        draw_color = color
                    # endif

                    ax.errorbar(
                        x_values,
                        values,
                        yerr=errors,
                        marker=marker,
                        linestyle="none",
                        capsize=2,
                        color=draw_color,
                        label=label,
                    )
                # endfor

                ax.axhline(0.0, linewidth=0.8)
                ax.set_ylabel(PARAMETER_LABELS[parameter])
                apply_parameter_y_limits(ax, parameter)
                ax.grid(alpha=0.25)
                if panel >= 6:
                    ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
                # endif
            # endfor

            handles, labels = axes_flat[0].get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                loc="lower right",
                bbox_to_anchor=(0.97, 0.055),
            )
            fig.suptitle(
                rf"{PERIOD_LABELS[period]}, "
                rf"${x_low:.2f} \leq x_B < {x_high:.2f}$: fit stability",
                y=0.995,
            )
            fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
            path = (
                output_dir
                / f"{period}_xB_bin_{x_index + 1}_stability.png"
            )
            fig.savefig(path, dpi=180)
            plt.close(fig)
            paths.append(str(path))
        # endfor
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
    path = output_dir / "period_consistency_delta_nll.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)

def write_latex_table(
    frame: pd.DataFrame,
    path: Path,
    *,
    include_target_axis_uncertainty: bool,
    sample_variant: str,
) -> None:
    """Write a LaTeX table appropriate for the selected sample variant.

    The nominal extraction includes the target-axis envelope as a second
    uncertainty.  The ISR extraction intentionally omits the target-axis
    study, so its table contains only the statistical uncertainty rather than
    attempting to convert the corresponding ``None`` values to floats.
    """
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
            if include_target_axis_uncertainty:
                axis_raw = getattr(row, f"{parameter}_target_axis_sys")
                axis = float(axis_raw)
                entries.append(
                    rf"${value:.5f}\pm{stat:.5f}\pm{axis:.5f}$"
                )
            else:
                entries.append(rf"${value:.5f}\pm{stat:.5f}$")
            # endif
        # endfor
        lines.append(" & ".join(entries) + r" \\")
    # endfor

    if include_target_axis_uncertainty:
        caption = (
            r"Nominal simultaneous unbinned-likelihood results. "
            r"The first uncertainty is statistical and includes the "
            r"Gaussian-constrained dilution-factor statistical uncertainty. "
            r"The second is the target-axis treatment envelope. "
            r"Polarization and dilution-model scale uncertainties are not "
            r"included."
        )
    else:
        caption = (
            rf"{sample_variant.upper()} simultaneous unbinned-likelihood "
            r"diagnostic results. The quoted uncertainty is statistical and "
            r"includes the Gaussian-constrained dilution-factor statistical "
            r"uncertainty. No target-axis envelope is evaluated for the ISR "
            r"sample. Polarization and dilution-model scale uncertainties "
            r"are not included."
        )
    # endif

    lines.extend(
        [
            r"\hline",
            r"\end{tabular}",
            rf"\caption{{{caption}}}",
            r"\end{table}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")



def write_external_transverse_input_tables(
    csv_path: Path,
    json_path: Path,
) -> None:
    """Write the fixed external leakage inputs and full provenance."""
    rows: list[dict[str, Any]] = []
    for term, payload in EXTERNAL_TRANSVERSE_INPUTS.items():
        reference = payload["reference"]
        kinematics = payload["kinematic_selection"]
        rows.append(
            {
                "term": term,
                "input_value": payload["value"],
                "measured_central_value": payload[
                    "measured_central_value"
                ],
                "scale_factor": payload["scale_factor"],
                "observable": payload["observable"],
                "channel": payload["channel"],
                "kinematic_selection_json": json.dumps(
                    kinematics,
                    sort_keys=True,
                ),
                "citation": reference["citation"],
                "title": reference["title"],
                "doi": reference["doi"],
                "arxiv": reference["arxiv"],
                "figure": reference.get("figure"),
                "data_note": reference["data_note"],
            }
        )
    # endfor
    ensure_directory(csv_path.parent)
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    write_json(
        json_path,
        {
            "policy": (
                "The UT inputs are direct inverse-total-variance weighted "
                "averages of exclusive-pi+ measurements. The LT inputs are "
                "rounded estimates from the open highest-z SIDIS pi+ points. "
                "No factor of two is applied."
            ),
            "inputs": EXTERNAL_TRANSVERSE_INPUTS,
        },
    )



# =============================================================================
# Command line
# =============================================================================

def resolve_isr_cut_json(explicit_path: Path | None) -> Path:
    if explicit_path is not None:
        path = explicit_path.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Explicit ISR cut JSON does not exist: {path}")
        return path
    # endif
    manifest = DEFAULT_CHANNEL_SELECTION_MANIFEST.expanduser().resolve()
    if manifest.is_file():
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        retained = payload.get("retained_analysis_products", {})
        isr_root = retained.get("isr")
        if isr_root:
            candidate = Path(isr_root) / "final_carbon_assisted_cuts/tables/final_carbon_assisted_mx2_cuts.json"
            if candidate.is_file():
                return candidate.resolve()
            # endif
        # endif
    # endif
    path = DEFAULT_ISR_CUT_JSON.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Could not resolve ISR channel-selection JSON: {path}")
    return path


def run_analysis_variant(
    *,
    sample_variant: str,
    input_paths: dict[str, Path],
    run_info_path: Path,
    cut_json_path: Path,
    dilution_json_path: Path,
    output_dir: Path,
    cache_path: Path,
    tree_name: str,
    chunk_size: str,
    workers: int,
    reuse_cache: bool,
    skip_plots: bool,
    include_target_axis_study: bool,
) -> dict[str, Any]:
    tables_dir = output_dir / "tables"
    json_dir = output_dir / "json"
    covariance_dir = output_dir / "covariance"
    plots_dir = output_dir / "plots"
    all_bins_plots_dir = plots_dir / "all_bins"
    aggregated_plots_dir = plots_dir / "aggregated"
    period_plots_dir = aggregated_plots_dir / "by_period"
    target_axis_variants_dir = plots_dir / "target_axis_variants"
    period_stability_dir = plots_dir / "period_stability"
    latex_dir = output_dir / "latex"
    for directory in (
        output_dir, tables_dir, json_dir, covariance_dir, plots_dir,
        all_bins_plots_dir, aggregated_plots_dir, period_plots_dir,
        period_stability_dir, latex_dir,
    ):
        ensure_directory(directory)
    # endfor
    if include_target_axis_study:
        ensure_directory(target_axis_variants_dir)
    # endif

    print("-" * 78)
    print(f"Starting {sample_variant} structure-function-ratio analysis")
    print("-" * 78)
    print(f"Channel cuts:         {cut_json_path}")
    print(f"Dilution factors:     {dilution_json_path}")
    print(f"Output directory:     {output_dir}")
    print(f"Selected-event cache: {cache_path}")
    print(f"Target-axis study:    {include_target_axis_study}")

    run_records = parse_run_info_csv(run_info_path)
    run_states = run_state_arrays(run_records)
    cuts = load_nominal_cuts(cut_json_path)
    dilution_records = load_dilution_factors(dilution_json_path)
    if reuse_cache:
        events = load_event_cache(cache_path)
        cache_summary = {"cache_path": str(cache_path), "number_of_selected_events": int(events["runnum"].size), "reused": True}
    else:
        cache_summary = build_event_cache(input_paths=input_paths, tree_name=tree_name, chunk_size=chunk_size, run_records=run_records, cuts=cuts, cache_path=cache_path)
        cache_summary["reused"] = False
        events = load_event_cache(cache_path)
    # endif

    run_state_payload = {period: {key: values.tolist() for key, values in state.items()} for period, state in run_states.items()}
    dilution_payload = {period: {str(bin_number): {"x_index": record.x_index, "t_index": record.t_index, "value": record.value, "stat_uncertainty": record.stat_uncertainty} for (record_period, bin_number), record in dilution_records.items() if record_period == period} for period in PERIODS}
    results = []
    with ProcessPoolExecutor(max_workers=workers, initializer=initialize_fit_worker, initargs=(str(cache_path), run_state_payload, dilution_payload)) as executor:
        futures = {executor.submit(fit_bin_worker, bin_number, include_target_axis_study): bin_number for bin_number in range(1, NUMBER_OF_BINS + 1)}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            nominal = result["variants"]["nominal"]
            print(f"[{sample_variant} bin {result['bin_number']:02d}] N={nominal['metadata']['number_of_events']:,}; valid={nominal['valid']}; NLL={nominal['minimum_nll']:.6f}; EDM={nominal['edm']:.3e}")
        # endfor
    # endwith
    results.sort(key=lambda item: item["bin_number"])
    frame = flatten_fit_results(results)
    csv_path = tables_dir / "structure_function_ratios.csv"
    frame.to_csv(csv_path, index=False)
    detailed_json_path = json_dir / "structure_function_ratios.json"
    write_json(detailed_json_path, {
        "schema_version": 2,
        "analysis": "RGC exclusive enpi+ structure-function-ratio extraction",
        "sample_variant": sample_variant,
        "diagnostic_only": sample_variant == "isr",
        "target_axis_study_performed": include_target_axis_study,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "beam_polarization": BEAM_POLARIZATION,
        "beam_energy_gev": BEAM_ENERGY_GEV,
        "cache": cache_summary,
        "inputs": {
            "root_files": {period: str(input_paths[period].expanduser().resolve()) for period in PERIODS},
            "run_info_csv": str(run_info_path), "run_info_csv_sha256": sha256_file(run_info_path),
            "channel_cut_json": str(cut_json_path), "channel_cut_json_sha256": sha256_file(cut_json_path),
            "dilution_json": str(dilution_json_path), "dilution_json_sha256": sha256_file(dilution_json_path),
        },
        "results": results,
    })
    for result in results:
        nominal = result["variants"]["nominal"]
        if nominal["covariance"] is not None:
            b = result["bin_number"]
            np.save(covariance_dir / f"bin_{b:02d}_nominal_covariance.npy", np.asarray(nominal["covariance"], dtype=np.float64))
            np.save(covariance_dir / f"bin_{b:02d}_nominal_correlation.npy", np.asarray(nominal["correlation"], dtype=np.float64))
        # endif
    # endfor
    latex_path = latex_dir / "structure_function_ratios.tex"
    write_latex_table(
        frame,
        latex_path,
        include_target_axis_uncertainty=include_target_axis_study,
        sample_variant=sample_variant,
    )
    plot_paths = {"all_bins": [], "aggregated": [], "aggregated_by_period": [], "target_axis_variants": [], "period_stability": []}
    if not skip_plots:
        plot_paths["all_bins"] = plot_parameter_summaries(frame, all_bins_plots_dir)
        plot_paths["aggregated"] = plot_aggregated_by_x(frame, aggregated_plots_dir)
        plot_paths["aggregated_by_period"] = plot_aggregated_by_period(frame, period_plots_dir)
        plot_paths["aggregated_by_period"].append(plot_period_consistency_heatmap(frame, period_plots_dir))
        plot_paths["period_stability"] = plot_period_stability(frame, period_stability_dir)
        if include_target_axis_study:
            plot_paths["target_axis_variants"] = plot_target_axis_variants(frame, target_axis_variants_dir)
        # endif
    # endif
    manifest_path = output_dir / "analysis_variant_manifest.json"
    write_json(manifest_path, {"schema_version": 2, "sample_variant": sample_variant, "diagnostic_only": sample_variant == "isr", "target_axis_study_performed": include_target_axis_study, "products": {"csv": str(csv_path), "detailed_json": str(detailed_json_path), "latex": str(latex_path), "covariance_directory": str(covariance_dir), "plots": plot_paths, "cache": str(cache_path)}})
    invalid_bins = frame.loc[~frame["nominal_fit_valid"].astype(bool), "bin_number"].astype(int).tolist()
    return {"sample_variant": sample_variant, "frame": frame, "results": results, "events": int(events["runnum"].size), "csv": str(csv_path), "json": str(detailed_json_path), "latex": str(latex_path), "manifest": str(manifest_path), "invalid_bins": invalid_bins}


def write_nominal_isr_comparison_products(nominal: pd.DataFrame, isr: pd.DataFrame, output_dir: Path) -> dict[str, str]:
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    ensure_directory(tables_dir); ensure_directory(plots_dir)
    keys = ["bin_number", "x_index", "t_index"]
    keep = keys + ["mean_xB", "mean_Q2_gev2", "mean_minus_tprime_gev2", "number_of_events"] + [item for p in PHYSICS_PARAMETERS for item in (p, f"{p}_stat")]
    merged = nominal[keep].merge(isr[keep], on=keys, suffixes=("_nominal", "_isr"), validate="one_to_one")
    for parameter in PHYSICS_PARAMETERS:
        merged[f"{parameter}_isr_minus_nominal"] = merged[f"{parameter}_isr"] - merged[f"{parameter}_nominal"]
        merged[f"{parameter}_absolute_isr_systematic"] = merged[f"{parameter}_isr_minus_nominal"].abs()
        merged[f"{parameter}_combined_stat_uncertainty"] = np.hypot(merged[f"{parameter}_stat_nominal"], merged[f"{parameter}_stat_isr"])
        merged[f"{parameter}_difference_pull"] = merged[f"{parameter}_isr_minus_nominal"] / merged[f"{parameter}_combined_stat_uncertainty"].replace(0.0, np.nan)
    # endfor
    csv_path = tables_dir / "nominal_vs_isr_structure_function_ratios.csv"
    json_path = tables_dir / "nominal_vs_isr_structure_function_ratios.json"
    merged.to_csv(csv_path, index=False)
    write_json(json_path, {"schema_version": 1, "systematic_convention": "absolute point-by-point ISR systematic is abs(ISR - nominal)", "rows": merged.to_dict(orient="records")})
    plot_paths=[]
    x=merged["bin_number"].to_numpy()
    for parameter in PHYSICS_PARAMETERS:
        fig, axes = plt.subplots(2,1,figsize=(11,8),sharex=True)
        axes[0].errorbar(x, merged[f"{parameter}_nominal"], yerr=merged[f"{parameter}_stat_nominal"], fmt="o", capsize=2, label="Nominal")
        axes[0].errorbar(x, merged[f"{parameter}_isr"], yerr=merged[f"{parameter}_stat_isr"], fmt="s", capsize=2, label="ISR")
        axes[0].set_ylabel(PARAMETER_LABELS[parameter]); axes[0].legend(); axes[0].grid(alpha=0.25)
        axes[1].axhline(0.0, linewidth=1)
        axes[1].errorbar(x, merged[f"{parameter}_isr_minus_nominal"], yerr=merged[f"{parameter}_combined_stat_uncertainty"], fmt="o", capsize=2)
        axes[1].set_xlabel("Combined kinematic-bin number"); axes[1].set_ylabel("ISR - nominal"); axes[1].grid(alpha=0.25)
        fig.tight_layout()
        path=plots_dir / f"nominal_vs_isr_{parameter}.png"; fig.savefig(path,dpi=200); plt.close(fig); plot_paths.append(str(path))
    # endfor
    return {"csv": str(csv_path), "json": str(json_path), "plots_directory": str(plots_dir), "plots": plot_paths}


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fit nominal and ISR RGC exclusive-pi+ structure-function ratios.")
    parser.add_argument("--tree", default=DEFAULT_TREE_NAME)
    parser.add_argument("--input", action="append", default=[], type=parse_input_override, metavar="PERIOD=FILE")
    parser.add_argument("--isr-input", action="append", default=[], type=parse_input_override, metavar="PERIOD=FILE")
    parser.add_argument("--run-info-csv", type=Path, default=DEFAULT_RUN_INFO_CSV)
    parser.add_argument("--cut-json", type=Path, default=DEFAULT_CUT_JSON)
    parser.add_argument("--isr-cut-json", type=Path, default=None)
    parser.add_argument("--dilution-json", type=Path, default=None)
    parser.add_argument("--dilution-dir", type=Path, default=DEFAULT_DILUTION_DIR)
    parser.add_argument("--isr-dilution-json", type=Path, default=DEFAULT_ISR_DILUTION_JSON)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache", type=Path, default=None, help="Legacy nominal-cache override.")
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--disable-isr", action="store_true")
    parser.add_argument("--chunk-size", default=DEFAULT_CHUNK_SIZE)
    parser.add_argument("--workers", type=int, default=MAXIMUM_WORKERS)
    parser.add_argument("--skip-plots", action="store_true")
    return parser


def main() -> int:
    args = build_argument_parser().parse_args()
    workers=max(1,min(int(args.workers),MAXIMUM_WORKERS,os.cpu_count() or 1,NUMBER_OF_BINS))
    root=args.output_dir.expanduser().resolve(); nominal_dir=root/"nominal"; isr_dir=root/"isr"; diagnostics_dir=root/"diagnostics"
    for d in (root,nominal_dir,isr_dir,diagnostics_dir): ensure_directory(d)
    nominal_inputs={p:Path(v) for p,v in DEFAULT_INPUTS.items()}
    for period,path in args.input: nominal_inputs[period]=path
    nominal_dilution=(args.dilution_json.expanduser().resolve() if args.dilution_json else find_default_dilution_json(args.dilution_dir.expanduser().resolve()).resolve())
    nominal_cache=(args.cache.expanduser().resolve() if args.cache else nominal_dir/"cache/selected_events.npz")
    nominal_result=run_analysis_variant(sample_variant="nominal", input_paths=nominal_inputs, run_info_path=args.run_info_csv.expanduser().resolve(), cut_json_path=args.cut_json.expanduser().resolve(), dilution_json_path=nominal_dilution, output_dir=nominal_dir, cache_path=nominal_cache, tree_name=args.tree, chunk_size=args.chunk_size, workers=workers, reuse_cache=args.reuse_cache, skip_plots=args.skip_plots, include_target_axis_study=True)
    isr_result=None; comparison=None
    if not args.disable_isr:
        isr_inputs={p:Path(v) for p,v in DEFAULT_ISR_INPUTS.items()}
        for period,path in args.isr_input: isr_inputs[period]=path
        isr_cut=resolve_isr_cut_json(args.isr_cut_json)
        isr_dilution=args.isr_dilution_json.expanduser().resolve()
        if not isr_dilution.is_file(): raise FileNotFoundError(f"ISR dilution JSON does not exist: {isr_dilution}")
        isr_result=run_analysis_variant(sample_variant="isr", input_paths=isr_inputs, run_info_path=args.run_info_csv.expanduser().resolve(), cut_json_path=isr_cut, dilution_json_path=isr_dilution, output_dir=isr_dir, cache_path=isr_dir/"cache/selected_events.npz", tree_name=args.tree, chunk_size=args.chunk_size, workers=workers, reuse_cache=args.reuse_cache, skip_plots=args.skip_plots, include_target_axis_study=False)
        comparison=write_nominal_isr_comparison_products(nominal_result["frame"],isr_result["frame"],diagnostics_dir)
    # endif
    manifest_path=root/"asymmetry_extraction_manifest.json"
    write_json(manifest_path,{"schema_version":3,"production_policy":"Nominal results are production. ISR results are a separate point-by-point diagnostic. Target-axis leakage is evaluated only for nominal data.","nominal":{k:v for k,v in nominal_result.items() if k not in ("frame","results")},"isr":({k:v for k,v in isr_result.items() if k not in ("frame","results")} if isr_result else None),"nominal_isr_comparison":comparison})
    print("\nStructure-function-ratio study complete.")
    print(f"  Nominal:    {nominal_dir}")
    if isr_result: print(f"  ISR:        {isr_dir}\n  Comparison: {diagnostics_dir}")
    print(f"  Manifest:   {manifest_path}")
    invalid=nominal_result["invalid_bins"] + ([] if isr_result is None else isr_result["invalid_bins"])
    return 2 if invalid else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
