#!/usr/bin/env python3
"""
determine_dilution_factor_v13.py

Determine the RGC exclusive e pi+ dilution factor in the fixed 4 xB by 6
(-tprime) bins using three complementary methods:

  Method 1
    Carbon subtraction with exactly the period-wide raw-count normalization
    used by channel_selection_mx2_fits_v24.py:

        alpha_C = N_NH3(0.00 <= Mx2 < 0.40) /
                  N_C  (0.00 <= Mx2 < 0.40)

        f = (N_NH3 - alpha_C N_C) / N_NH3.

  Method 2
    Direct auxiliary-target dilution factor from Eq. (10) of the standalone
    RGC dilution-factor note, using charge-normalized rates from NH3, C, CH2,
    He-bath (MT), and empty-target/foils-only (F) data.

  Method 3
    Packing fraction from Eq. (11), followed by the nonlinear dilution-factor
    expression in Eq. (14).  A per-bin packing fraction is retained as a
    diagnostic.  The production-style Method-3 calculation uses a period-wide
    packing fraction extracted from all 24 kinematic bins for the same Mx2 cut
    variation, because packing fraction is a target-cell property rather than
    a kinematic observable.

The program is intended to be run from

    RGC_enpi+/dilution_factor/

while the channel-selection outputs are in the sibling directory

    RGC_enpi+/channel_selection/

The default exclusivity-cut JSON is therefore

    ../channel_selection/output/channel_selection_mx2_fit_stability/
        final_carbon_assisted_cuts/tables/
        final_carbon_assisted_mx2_cuts.json

All ROOT inputs default to the finalized momentum-corrected paper_versions
files.  The program uses at most seven worker processes.  It writes complete
JSON and CSV products, compact downstream JSON, statistical covariance and
correlation matrices, and diagnostic plots.

Statistical treatment
---------------------
Counts are modeled as a Poisson point process.  For each period, target, and
kinematic bin, events are classified by their complete membership pattern in
four possibly overlapping selections:

    control, tight, nominal, loose.

The 16 disjoint membership-pattern counts are independently Poisson-resampled.
Summing those replicas reconstructs the selected counts while preserving the
correct statistical covariance for nested or overlapping windows.  The shared
Method-1 carbon scale, the shared period packing fraction, and all three cut
variations are recalculated in every replica.  This yields statistically
consistent uncertainties and full 24-bin covariance matrices.

No dilution factor or packing fraction is clipped to a physical interval.
Out-of-range values are retained and flagged as diagnostics.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET
import math
import os
import shutil
from pathlib import Path
import sys
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot


# =============================================================================
# Fixed analysis definitions
# =============================================================================

PERIODS: tuple[str, ...] = ("su22", "fa22", "sp23")
TARGETS: tuple[str, ...] = ("NH3", "C", "CH2", "He", "ET")
CUT_VARIATIONS: tuple[str, ...] = ("tight", "nominal", "loose")
SELECTION_NAMES: tuple[str, ...] = ("control", "tight", "nominal", "loose")

PERIOD_LABELS: dict[str, str] = {
    "su22": "Su22",
    "fa22": "Fa22",
    "sp23": "Sp23",
}

TARGET_LABELS: dict[str, str] = {
    "NH3": "NH$_3$",
    "C": "C",
    "CH2": "CH$_2$",
    "He": "He bath",
    "ET": "Empty target",
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
MAXIMUM_WORKERS = 7

DEFAULT_TREE_NAME = "PhysicsEvents"
DEFAULT_CONTROL_MIN_GEV2 = 0.0
DEFAULT_CONTROL_MAX_GEV2 = 0.40
DEFAULT_REPLICAS = 10000
DEFAULT_SEED = 7302026

DEFAULT_OUTPUT_DIR = Path("output/dilution_factor_determination")
DEFAULT_EPOCH_DEFINITIONS = Path(__file__).resolve().with_name("epoch_definitions.xlsx")
DEFAULT_RUN_INFO_CSV = Path(__file__).resolve().with_name("clas12_run_info.csv")
DEFAULT_CUT_JSON = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "final_carbon_assisted_cuts/tables/"
    "final_carbon_assisted_mx2_cuts.json"
)

PAPER_VERSIONS_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
    "paper_versions"
)

DEFAULT_INPUTS: dict[str, dict[str, Path]] = {
    period: {
        target: PAPER_VERSIONS_DIR
        / f"rgc_{period}_inb_{target}_epi+_mom_corrections.root"
        for target in TARGETS
    }
    for period in PERIODS
}

DEFAULT_ISR_INPUTS: dict[str, dict[str, Path]] = {
    period: {
        target: PAPER_VERSIONS_DIR
        / f"rgc_{period}_inb_{target}_epi+_ISR_externalISR_mom_corrections.root"
        for target in TARGETS
    }
    for period in PERIODS
}

DEFAULT_CHANNEL_SELECTION_MANIFEST = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "channel_selection_manifest.json"
)

# Relative accumulated-charge fractions supplied for this analysis.  They are
# sufficient because the auxiliary-target equations are homogeneous in the
# charge-normalized rates; the common total charge cancels.
DEFAULT_CHARGE_FRACTIONS: dict[str, dict[str, float]] = {
    "su22": {
        "NH3": 0.7624,
        "C": 0.0747,
        "CH2": 0.0390,
        "He": 0.0249,
        "ET": 0.0989,
    },
    "fa22": {
        "NH3": 0.5839,
        "C": 0.1992,
        "CH2": 0.1802,
        "He": 0.0301,
        "ET": 0.0066,
    },
    "sp23": {
        "NH3": 0.4543,
        "C": 0.1509,
        "CH2": 0.1966,
        "He": 0.1205,
        "ET": 0.0777,
    },
}

BRANCH_ALIASES: dict[str, tuple[str, ...]] = {
    "xB": ("xB", "x", "xb", "x_b"),
    "tprime": ("tprime", "t_prime", "tp", "tPrime"),
    "runnum": ("runnum", "run", "run_number", "RunNumber"),
    "Mx2": (
        "Mx2",
        "mx2",
        "Mx2_epi",
        "Mx2_epip",
        "missing_mass_squared",
        "missing_mass2",
    ),
}

# Numerical coefficients from Eqs. (10), (11), and (14) for 15NH3.
# In the note's notation:
#   A  -> NH3
#   C  -> carbon
#   CH -> CH2
#   MT -> cell/foils immersed in liquid helium
#   F  -> foils-only / empty target
EQ10_EQ11_EQ14_COEFFICIENTS: dict[str, float] = {
    "ch2_carbon_hydrogen_combination": 0.880502,
    "he_hydrogen_combination": 0.16604,
    "empty_hydrogen_combination": -0.285536,
    "ch2_carbon_nuclear_combination": 0.195667,
    "he_nuclear_combination": -0.82744,
    "empty_nuclear_combination": 0.023106,
    "packing_prefactor": 0.50734,
}


# =============================================================================
# Data containers
# =============================================================================

@dataclass(frozen=True)
class DatasetSpec:
    period: str
    target: str
    file_path: str
    tree_name: str


@dataclass(frozen=True)
class LoadedDataset:
    period: str
    target: str
    file_path: str
    tree_name: str
    x_branch: str
    tprime_branch: str
    mx2_branch: str
    run_branch: str
    runnum: np.ndarray
    xB: np.ndarray
    minus_tprime_gev2: np.ndarray
    mx2_gev2: np.ndarray


@dataclass(frozen=True)
class CutEntry:
    period: str
    x_index: int
    t_index: int
    xB_min: float
    xB_max: float
    minus_tprime_min_gev2: float
    minus_tprime_max_gev2: float
    mu_gev2: float
    mu_error_gev2: float
    sigma_gev2: float
    sigma_error_gev2: float
    classification: str
    tight_min_gev2: float
    tight_max_gev2: float
    nominal_min_gev2: float
    nominal_max_gev2: float
    loose_min_gev2: float
    loose_max_gev2: float

    def interval(self, variation: str) -> tuple[float, float]:
        if variation == "tight":
            return self.tight_min_gev2, self.tight_max_gev2
        # endif
        if variation == "nominal":
            return self.nominal_min_gev2, self.nominal_max_gev2
        # endif
        if variation == "loose":
            return self.loose_min_gev2, self.loose_max_gev2
        # endif
        raise KeyError(f"Unknown cut variation: {variation}")


@dataclass(frozen=True)
class CountPayload:
    """Observed counts and disjoint membership-pattern counts."""

    period: str
    target: str
    x_index: int
    t_index: int
    bin_number: int
    total_kinematic_count: int
    selected_counts: tuple[int, int, int, int]
    pattern_counts: tuple[int, ...]


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
    if isinstance(value, (np.integer,)):
        return int(value)
    # endif
    if isinstance(value, (np.floating,)):
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


def parse_input_override(text: str) -> tuple[str, str, Path]:
    try:
        left, path_text = text.split("=", 1)
        period, target = left.split(":", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Input overrides must have the form PERIOD:TARGET=/path/file.root"
        ) from exc
    # endtry
    period = period.strip().lower()
    target_lookup = {item.lower(): item for item in TARGETS}
    if period not in PERIODS:
        raise argparse.ArgumentTypeError(
            f"Unknown period {period!r}; expected one of {PERIODS}."
        )
    # endif
    if target.strip().lower() not in target_lookup:
        raise argparse.ArgumentTypeError(
            f"Unknown target {target!r}; expected one of {TARGETS}."
        )
    # endif
    return period, target_lookup[target.strip().lower()], Path(path_text).expanduser()



def resolve_isr_cut_json(
    explicit_path: Path | None,
    channel_selection_manifest: Path,
) -> Path:
    """
    Resolve the ISR-specific channel-selection cut JSON.

    Preference order:
      1. --isr-exclusivity-json
      2. channel-selection manifest diagnostics/retained ISR manifest
      3. ../channel_selection/.../isr/analysis_variant_manifest.json
    """
    if explicit_path is not None:
        path = explicit_path.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Explicit ISR cut JSON does not exist: {path}")
        return path

    manifest_path = channel_selection_manifest.expanduser().resolve()
    candidate_variant_manifest: Path | None = None
    if manifest_path.is_file():
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        retained = payload.get("retained_analysis_products", {})
        isr_root = retained.get("isr")
        if isr_root:
            candidate_variant_manifest = Path(isr_root) / "analysis_variant_manifest.json"

    if candidate_variant_manifest is None:
        candidate_variant_manifest = (
            manifest_path.parent / "isr" / "analysis_variant_manifest.json"
        )

    if not candidate_variant_manifest.is_file():
        raise FileNotFoundError(
            "Could not locate the ISR channel-selection variant manifest. "
            f"Checked: {candidate_variant_manifest}. Run channel selection with "
            "its ISR workflow first, or pass --isr-exclusivity-json explicitly."
        )

    variant_payload = json.loads(
        candidate_variant_manifest.read_text(encoding="utf-8")
    )
    cut_json_text = variant_payload.get("final_cut_json")
    if not cut_json_text:
        raise RuntimeError(
            f"ISR variant manifest has no final_cut_json: {candidate_variant_manifest}"
        )
    cut_path = Path(cut_json_text).expanduser().resolve()
    if not cut_path.is_file():
        raise FileNotFoundError(
            f"ISR cut JSON listed in the variant manifest does not exist: {cut_path}"
        )
    return cut_path


def validate_fraction_table(fractions: dict[str, dict[str, float]]) -> None:
    for period in PERIODS:
        if period not in fractions:
            raise RuntimeError(f"Missing charge fractions for {period}.")
        # endif
        for target in TARGETS:
            if target not in fractions[period]:
                raise RuntimeError(
                    f"Missing charge fraction for {period}/{target}."
                )
            # endif
            value = float(fractions[period][target])
            if not math.isfinite(value) or value <= 0.0:
                raise RuntimeError(
                    f"Invalid charge fraction for {period}/{target}: {value}"
                )
            # endif
        # endfor
        total = sum(float(fractions[period][target]) for target in TARGETS)
        if abs(total - 1.0) > 5.0e-4:
            raise RuntimeError(
                f"Charge fractions for {period} sum to {total:.8f}, not 1."
            )
        # endif
    # endfor


def load_charge_fractions(path: Path | None) -> dict[str, dict[str, float]]:
    if path is None:
        fractions = {
            period: dict(values)
            for period, values in DEFAULT_CHARGE_FRACTIONS.items()
        }
    else:
        payload = json.loads(path.read_text(encoding="utf-8"))
        source = payload.get("charge_fractions", payload)
        fractions = {
            str(period).lower(): {
                str(target): float(value)
                for target, value in values.items()
            }
            for period, values in source.items()
        }
    # endif
    validate_fraction_table(fractions)
    return fractions


def resolve_branch(tree: uproot.behaviors.TTree.TTree, aliases: Iterable[str]) -> str:
    available = set(str(key).split(";")[0] for key in tree.keys())
    for alias in aliases:
        if alias in available:
            return alias
        # endif
    # endfor
    raise RuntimeError(
        f"None of the branch aliases {tuple(aliases)} exists. "
        f"Available branches include: {sorted(available)[:80]}"
    )


def resolve_tree(
    root_file: uproot.reading.ReadOnlyDirectory,
    requested_name: str,
    path: Path,
):
    """
    Resolve a ROOT tree robustly.

    Some external-ISR files are written with key metadata for which
    ``requested_name in root_file`` is false even though
    ``root_file.keys()`` visibly contains the requested tree. Direct lookup
    still succeeds in those files, so membership testing is not used as the
    gatekeeper.
    """
    lookup_names: list[str] = [requested_name]

    for key in root_file.keys():
        key_text = str(key)
        base_name = key_text.split(";")[0]
        if base_name == requested_name:
            lookup_names.extend((key_text, base_name))
        # endif
    # endfor

    seen: set[str] = set()
    lookup_errors: list[str] = []
    for name in lookup_names:
        if name in seen:
            continue
        # endif
        seen.add(name)
        try:
            obj = root_file[name]
        except Exception as exc:
            lookup_errors.append(f"{name!r}: {exc}")
            continue
        # endtry

        if hasattr(obj, "arrays") and hasattr(obj, "keys"):
            return obj
        # endif
    # endfor

    tree_like: list[tuple[str, Any]] = []
    for key in root_file.keys():
        key_text = str(key)
        try:
            obj = root_file[key_text]
        except Exception:
            continue
        # endtry
        if hasattr(obj, "arrays") and hasattr(obj, "keys"):
            tree_like.append((key_text, obj))
        # endif
    # endfor

    if len(tree_like) == 1:
        return tree_like[0][1]
    # endif

    available = [str(key) for key in root_file.keys()]
    details = "; ".join(lookup_errors[:5])
    raise RuntimeError(
        f"Tree {requested_name!r} not found in {path}. "
        f"Available keys: {available}. "
        f"Direct-lookup diagnostics: {details}"
    )




@dataclass(frozen=True)
class RunChargeRecord:
    """One run entry from clas12_run_info.csv."""

    period: str
    target: str
    run: int
    charge: float


def parse_run_info_csv(path: Path) -> dict[tuple[str, str], dict[int, float]]:
    """
    Parse the sectioned CLAS12 run-information CSV.

    The first numeric field is run number. Accumulated charge is defined as
    the sum of the third and fourth CSV columns, matching the requested
    helicity-positive plus helicity-negative charge convention. Only the RGC
    Su22, Fa22, and Sp23 NH3/C/CH2/He/ET sections are retained.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Missing run-information CSV: {path}")
    # endif

    period_lookup = {
        "Su22": "su22",
        "Fa22": "fa22",
        "Sp23": "sp23",
    }
    target_lookup = {
        "NH3": "NH3",
        "C": "C",
        "CH2": "CH2",
        "He": "He",
        "ET": "ET",
    }

    result: dict[tuple[str, str], dict[int, float]] = {
        (period, target): {}
        for period in PERIODS
        for target in TARGETS
    }
    current_key: tuple[str, str] | None = None

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip()
        if not line:
            continue
        # endif
        if line.startswith("#"):
            current_key = None
            header = line[1:].strip()
            match = re.fullmatch(
                r"RGC\s+(Su22|Fa22|Sp23)\s+(NH3|C|CH2|He|ET)",
                header,
            )
            if match is not None:
                current_key = (
                    period_lookup[match.group(1)],
                    target_lookup[match.group(2)],
                )
            # endif
            continue
        # endif
        if current_key is None:
            continue
        # endif

        fields = [field.strip() for field in line.split(",")]
        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed run-info row at line {line_number}: {raw_line}"
            )
        # endif
        try:
            run = int(fields[0])
            charge_column_3 = float(fields[2])
            charge_column_4 = float(fields[3])
            charge = charge_column_3 + charge_column_4
        except ValueError as exc:
            raise RuntimeError(
                f"Invalid run/charge at line {line_number}: {raw_line}"
            ) from exc
        # endtry
        if not math.isfinite(charge) or charge < 0.0:
            raise RuntimeError(
                f"Invalid accumulated charge at line {line_number}: {charge}"
            )
        # endif
        if run in result[current_key]:
            raise RuntimeError(
                f"Duplicate run {run} in section {current_key}."
            )
        # endif
        result[current_key][run] = charge
    # endfor

    for key, run_map in result.items():
        if not run_map:
            raise RuntimeError(
                f"No run-charge entries found for {key[0]}/{key[1]} in {path}."
            )
        # endif
    # endfor
    return result


def absolute_charges_for_loaded_runs(
    loaded: dict[tuple[str, str], LoadedDataset],
    run_charges: dict[tuple[str, str], dict[int, float]],
) -> tuple[
    dict[str, dict[str, float]],
    dict[tuple[str, str], tuple[int, ...]],
]:
    """
    Sum charge only for unique run numbers that actually appear in ROOT trees.

    Rows present in clas12_run_info.csv but absent from the processed ROOT data
    are ignored. This is intentional because QADB processing may remove an
    entire run, leaving zero good charge and no events in the ROOT tree.
    """
    totals: dict[str, dict[str, float]] = {
        period: {} for period in PERIODS
    }
    used_runs: dict[tuple[str, str], tuple[int, ...]] = {}

    for period in PERIODS:
        for target in TARGETS:
            key = (period, target)
            unique_runs = tuple(
                int(run)
                for run in np.unique(loaded[key].runnum)
            )
            used_runs[key] = unique_runs

            charge_map = run_charges[key]
            missing_runs = [
                run for run in unique_runs
                if run not in charge_map
            ]
            if missing_runs:
                raise RuntimeError(
                    f"Missing accumulated charge for ROOT runs in "
                    f"{period}/{target}: {missing_runs}"
                )
            # endif

            nonpositive_runs = [
                run for run in unique_runs
                if charge_map[run] <= 0.0
            ]
            if nonpositive_runs:
                raise RuntimeError(
                    f"ROOT runs with nonpositive accumulated charge in "
                    f"{period}/{target}: {nonpositive_runs}. These runs "
                    "should not appear in the processed ROOT tree if QADB "
                    "removed all good charge."
                )
            # endif

            total = float(
                sum(charge_map[run] for run in unique_runs)
            )
            if not math.isfinite(total) or total <= 0.0:
                raise RuntimeError(
                    f"Nonpositive total charge for loaded ROOT runs in "
                    f"{period}/{target}: {total}"
                )
            # endif
            totals[period][target] = total
        # endfor
    # endfor

    return totals, used_runs


def normalized_charge_fractions_from_absolute(
    absolute_charges: dict[str, float],
) -> dict[str, float]:
    """Normalize positive absolute target charges to fractions summing to one."""
    total = float(sum(absolute_charges[target] for target in TARGETS))
    if not math.isfinite(total) or total <= 0.0:
        raise RuntimeError("Absolute target-charge total is not positive.")
    # endif
    fractions = {
        target: float(absolute_charges[target]) / total
        for target in TARGETS
    }
    return fractions




@dataclass(frozen=True)
class EpochDefinition:
    period: str
    epoch: str
    runs: tuple[int, ...]
    run_min: int
    run_max: int
    spreadsheet_event_count: int
    exposure_fraction: float
    weighted_polarization: float


def _xlsx_column_index(reference: str) -> int:
    letters = "".join(character for character in reference if character.isalpha())
    value = 0
    for character in letters:
        value = value * 26 + (ord(character.upper()) - ord("A") + 1)
    # endfor
    return value - 1


def read_epoch_spreadsheet(path: Path) -> list[EpochDefinition]:
    """
    Read NH3 epoch definitions from the supplied XLSX using only stdlib XML.

    The spreadsheet supplies epoch boundaries and polarization metadata. Accumulated charge is taken from clas12_run_info.csv.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Missing epoch-definition spreadsheet: {path}")
    # endif

    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("x:si", namespace):
                shared_strings.append(
                    "".join(node.text or "" for node in item.findall(".//x:t", namespace))
                )
            # endfor
        # endif
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    # endwith

    table: list[list[Any]] = []
    for row_node in sheet.findall(".//x:sheetData/x:row", namespace):
        values: dict[int, Any] = {}
        for cell in row_node.findall("x:c", namespace):
            reference = cell.attrib.get("r", "A1")
            column = _xlsx_column_index(reference)
            value_node = cell.find("x:v", namespace)
            if value_node is None:
                value: Any = ""
            else:
                raw = value_node.text or ""
                if cell.attrib.get("t") == "s":
                    value = shared_strings[int(raw)]
                elif cell.attrib.get("t") == "b":
                    value = bool(int(raw))
                else:
                    try:
                        numeric = float(raw)
                        value = int(numeric) if numeric.is_integer() else numeric
                    except ValueError:
                        value = raw
                    # endtry
                # endif
            # endif
            values[column] = value
        # endfor
        if values:
            width = max(values) + 1
            table.append([values.get(index, "") for index in range(width)])
        # endif
    # endfor

    if not table:
        raise RuntimeError(f"No rows found in epoch spreadsheet: {path}")
    # endif
    headers = [str(value).strip() for value in table[0]]
    required = {
        "Run Number",
        "Target",
        "Approx. Beam Energy",
        "Number of Events",
        "Target Polarization",
        "Epoch Number",
    }
    missing = required - set(headers)
    if missing:
        raise RuntimeError(
            f"Epoch spreadsheet is missing columns: {sorted(missing)}"
        )
    # endif
    column = {name: headers.index(name) for name in required}

    energy_to_period = {
        10.5473: "su22",
        10.5563: "fa22",
        10.5593: "sp23",
    }
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in table[1:]:
        if len(row) <= max(column.values()):
            continue
        # endif
        if str(row[column["Target"]]).strip().upper() != "NH3":
            continue
        # endif
        epoch = str(row[column["Epoch Number"]]).strip().strip("()")
        if not epoch or not epoch.upper().startswith("P"):
            continue
        # endif
        energy = float(row[column["Approx. Beam Energy"]])
        period = min(
            energy_to_period,
            key=lambda candidate: abs(candidate - energy),
        )
        if abs(period - energy) > 0.01:
            continue
        # endif
        period_name = energy_to_period[period]
        grouped.setdefault((period_name, epoch), []).append(
            {
                "run": int(row[column["Run Number"]]),
                "events": int(row[column["Number of Events"]]),
                "polarization": float(row[column["Target Polarization"]]),
            }
        )
    # endfor

    period_totals = {
        period: sum(
            item["events"]
            for (candidate_period, _), rows in grouped.items()
            if candidate_period == period
            for item in rows
        )
        for period in PERIODS
    }
    definitions: list[EpochDefinition] = []
    for (period, epoch), rows in sorted(
        grouped.items(),
        key=lambda item: (PERIODS.index(item[0][0]), min(row["run"] for row in item[1])),
    ):
        event_count = sum(row["events"] for row in rows)
        weighted_polarization = safe_divide(
            sum(row["events"] * row["polarization"] for row in rows),
            event_count,
        )
        runs = tuple(sorted(row["run"] for row in rows))
        definitions.append(
            EpochDefinition(
                period=period,
                epoch=epoch,
                runs=runs,
                run_min=min(runs),
                run_max=max(runs),
                spreadsheet_event_count=event_count,
                exposure_fraction=(
                    event_count / period_totals[period]
                    if period_totals[period] > 0
                    else math.nan
                ),
                weighted_polarization=float(weighted_polarization),
            )
        )
    # endfor
    return definitions


def count_integrated_epoch_nh3(
    dataset: LoadedDataset,
    epoch: EpochDefinition,
    cuts: dict[tuple[str, int, int], CutEntry],
    control_min_gev2: float,
    control_max_gev2: float,
) -> tuple[int, int, tuple[int, ...]]:
    """
    Count integrated control and nominal NH3 selections for one epoch.

    Epoch membership is defined by the inclusive run-number range. The actual
    runs are taken from the ROOT runnum branch, rather than inferred from the
    spreadsheet event totals.
    """
    run_mask = (
        (dataset.runnum >= epoch.run_min)
        & (dataset.runnum <= epoch.run_max)
    )
    runs_found = tuple(
        int(run)
        for run in np.unique(dataset.runnum[run_mask])
    )
    control_count = 0
    nominal_count = 0
    for x_index, (x_min, x_max) in enumerate(XB_BINS):
        for t_index, (t_min, t_max) in enumerate(MINUS_TPRIME_BINS_GEV2):
            mask = (
                run_mask
                & (dataset.xB >= x_min)
                & (dataset.xB < x_max)
                & (dataset.minus_tprime_gev2 >= t_min)
                & (dataset.minus_tprime_gev2 < t_max)
            )
            values = dataset.mx2_gev2[mask]
            control_count += int(
                np.count_nonzero(
                    (values >= control_min_gev2)
                    & (values < control_max_gev2)
                )
            )
            low, high = cuts[
                (epoch.period, x_index, t_index)
            ].interval("nominal")
            nominal_count += int(
                np.count_nonzero((values >= low) & (values < high))
            )
        # endfor
    # endfor
    return control_count, nominal_count, runs_found


def calculate_epoch_diagnostics(
    loaded: dict[tuple[str, str], LoadedDataset],
    observed_counts: np.ndarray,
    cuts: dict[tuple[str, int, int], CutEntry],
    definitions: list[EpochDefinition],
    run_charges: dict[tuple[str, str], dict[int, float]],
    control_min_gev2: float,
    control_max_gev2: float,
    replicas: int,
    seed: int,
) -> pd.DataFrame:
    """
    Calculate nominal dilution and packing-fraction diagnostics by NH3 epoch.

    For each epoch:
      * NH3 counts and charge use only NH3 runs in that epoch's run range.
      * The actual NH3 runs are read from the ROOT runnum branch.
      * C, CH2, He, and ET counts and charges remain full-period quantities, using only unique runs actually present in their ROOT trees.
      * Charge fractions are recomputed from those five absolute charges.
      * Method 1 is renormalized directly from epoch NH3 control counts to the
        same full-period carbon control sample; no exposure proxy is used.
    """
    rows: list[dict[str, Any]] = []
    nominal_selection_index = SELECTION_NAMES.index("nominal")
    control_selection_index = SELECTION_NAMES.index("control")
    rng = np.random.default_rng(seed + 9107)
    full_charges, full_period_root_runs = (
        absolute_charges_for_loaded_runs(
            loaded=loaded,
            run_charges=run_charges,
        )
    )

    period_epoch_counter = {period: 0 for period in PERIODS}
    for epoch in definitions:
        epoch_index = period_epoch_counter[epoch.period]
        period_epoch_counter[epoch.period] += 1
        p_index = PERIODS.index(epoch.period)

        (
            nh3_control,
            nh3_nominal,
            root_runs_found,
        ) = count_integrated_epoch_nh3(
            loaded[(epoch.period, "NH3")],
            epoch,
            cuts,
            control_min_gev2,
            control_max_gev2,
        )
        if not root_runs_found:
            continue
        # endif

        nh3_charge_map = run_charges[(epoch.period, "NH3")]
        missing_charge_runs = [
            run for run in root_runs_found
            if run not in nh3_charge_map
        ]
        if missing_charge_runs:
            raise RuntimeError(
                f"Missing NH3 accumulated charge for {epoch.period} "
                f"{epoch.epoch} ROOT runs: {missing_charge_runs}"
            )
        # endif
        nonpositive_charge_runs = [
            run for run in root_runs_found
            if nh3_charge_map[run] <= 0.0
        ]
        if nonpositive_charge_runs:
            raise RuntimeError(
                f"NH3 ROOT runs with nonpositive accumulated charge for "
                f"{epoch.period} {epoch.epoch}: "
                f"{nonpositive_charge_runs}. Such runs should have no "
                "events after QADB processing."
            )
        # endif
        epoch_nh3_charge = float(
            sum(nh3_charge_map[run] for run in root_runs_found)
        )

        period_counts = {
            target: int(
                np.sum(
                    observed_counts[
                        p_index,
                        TARGETS.index(target),
                        :,
                        nominal_selection_index,
                    ]
                )
            )
            for target in TARGETS
        }
        period_control_c = int(
            np.sum(
                observed_counts[
                    p_index,
                    TARGETS.index("C"),
                    :,
                    control_selection_index,
                ]
            )
        )

        # Method 1 uses the epoch NH3 control sample and the unchanged
        # full-period carbon control and nominal samples.
        alpha_epoch = (
            nh3_control / period_control_c
            if period_control_c > 0
            else math.nan
        )
        method1 = (
            1.0
            - alpha_epoch * period_counts["C"] / nh3_nominal
            if nh3_nominal > 0
            else math.nan
        )

        epoch_counts = dict(period_counts)
        epoch_counts["NH3"] = nh3_nominal

        epoch_absolute_charges = dict(full_charges[epoch.period])
        epoch_absolute_charges["NH3"] = epoch_nh3_charge
        epoch_charge_fractions = (
            normalized_charge_fractions_from_absolute(
                epoch_absolute_charges
            )
        )

        scalar_counts = {
            target: np.asarray(value, dtype=np.float64)
            for target, value in epoch_counts.items()
        }
        method2 = float(
            method2_equation10_from_counts(
                scalar_counts,
                epoch_charge_fractions,
            )
        )
        packing_fraction = float(
            packing_fraction_equation11_from_counts(
                scalar_counts,
                epoch_charge_fractions,
            )
        )
        recommended = 0.5 * (method1 + method2)

        # Independent-Poisson epoch diagnostic. The full-period auxiliary
        # samples are redrawn for every epoch because they enter every epoch
        # estimate, while only NH3 changes with epoch.
        draws = min(max(replicas, 500), 10000)
        replica_nh3_control = rng.poisson(
            nh3_control,
            size=draws,
        )
        replica_nh3_nominal = rng.poisson(
            nh3_nominal,
            size=draws,
        )
        replica_period_control_c = rng.poisson(
            period_control_c,
            size=draws,
        )
        replica_counts = {
            target: rng.poisson(
                period_counts[target],
                size=draws,
            ).astype(float)
            for target in TARGETS
        }
        replica_counts["NH3"] = replica_nh3_nominal.astype(float)
        replica_alpha = safe_divide(
            replica_nh3_control,
            replica_period_control_c,
        )
        replica_method1 = 1.0 - safe_divide(
            replica_alpha * replica_counts["C"],
            replica_nh3_nominal,
        )
        replica_method2 = method2_equation10_from_counts(
            replica_counts,
            epoch_charge_fractions,
        )
        replica_pf = packing_fraction_equation11_from_counts(
            replica_counts,
            epoch_charge_fractions,
        )
        replica_recommended = 0.5 * (
            replica_method1 + replica_method2
        )

        rows.append(
            {
                "period": epoch.period,
                "period_label": PERIOD_LABELS[epoch.period],
                "epoch": epoch.epoch,
                "epoch_index_within_period": epoch_index,
                "run_min": epoch.run_min,
                "run_max": epoch.run_max,
                "number_of_root_runs_found": len(root_runs_found),
                "root_runs_found": ",".join(
                    str(run) for run in root_runs_found
                ),
                "weighted_target_polarization": (
                    epoch.weighted_polarization
                ),
                "nh3_accumulated_charge": epoch_nh3_charge,
                "full_period_c_charge": (
                    full_charges[epoch.period]["C"]
                ),
                "full_period_ch2_charge": (
                    full_charges[epoch.period]["CH2"]
                ),
                "full_period_he_charge": (
                    full_charges[epoch.period]["He"]
                ),
                "full_period_et_charge": (
                    full_charges[epoch.period]["ET"]
                ),
                "full_period_c_root_runs": ",".join(
                    str(run)
                    for run in full_period_root_runs[
                        (epoch.period, "C")
                    ]
                ),
                "full_period_ch2_root_runs": ",".join(
                    str(run)
                    for run in full_period_root_runs[
                        (epoch.period, "CH2")
                    ]
                ),
                "full_period_he_root_runs": ",".join(
                    str(run)
                    for run in full_period_root_runs[
                        (epoch.period, "He")
                    ]
                ),
                "full_period_et_root_runs": ",".join(
                    str(run)
                    for run in full_period_root_runs[
                        (epoch.period, "ET")
                    ]
                ),
                "charge_fraction_nh3": (
                    epoch_charge_fractions["NH3"]
                ),
                "charge_fraction_c": (
                    epoch_charge_fractions["C"]
                ),
                "charge_fraction_ch2": (
                    epoch_charge_fractions["CH2"]
                ),
                "charge_fraction_he": (
                    epoch_charge_fractions["He"]
                ),
                "charge_fraction_et": (
                    epoch_charge_fractions["ET"]
                ),
                "nh3_control_count": nh3_control,
                "nh3_nominal_count": nh3_nominal,
                "full_period_c_control_count": period_control_c,
                "full_period_c_nominal_count": (
                    period_counts["C"]
                ),
                "method1_epoch_carbon_scale": alpha_epoch,
                "method1_value": method1,
                "method2_value": method2,
                "recommended_dilution_factor": recommended,
                "recommended_stat_uncertainty": float(
                    np.nanstd(
                        replica_recommended,
                        ddof=1,
                    )
                ),
                "packing_fraction": packing_fraction,
                "packing_fraction_stat_uncertainty": float(
                    np.nanstd(replica_pf, ddof=1)
                ),
                "exposure_note": (
                    "NH3 counts and charge are epoch-specific. C, CH2, He, "
                    "and ET counts and charges are unchanged full-period "
                    "quantities. Charge fractions are recomputed for every "
                    "epoch from those absolute charges."
                ),
            }
        )
    # endfor
    return pd.DataFrame(rows)


def plot_epoch_diagnostics(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """
    Plot epoch diagnostics on one common epoch-number axis from 1 through 24.

    Periods are distinguished by color. Colored dash-dot vertical boundaries
    mark the run-period ranges, and inverse-variance weighted period averages
    are shown as horizontal lines with one-standard-deviation bands.
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    period_colors = {
        "su22": "tab:blue",
        "fa22": "tab:orange",
        "sp23": "tab:green",
    }
    period_epoch_ranges = {
        "su22": (1, 12),
        "fa22": (13, 21),
        "sp23": (22, 24),
    }

    plot_frame = frame.copy()
    plot_frame["epoch_number"] = plot_frame["epoch"].map(
        lambda value: int(
            re.search(r"(\d+)", str(value)).group(1)
        )
    )

    for quantity, error, ylabel, stem, title in (
        (
            "recommended_dilution_factor",
            "recommended_stat_uncertainty",
            "Recommended dilution factor",
            "dilution_factor_by_epoch.png",
            "Recommended dilution factor versus NH$_3$ target/polarization epoch",
        ),
        (
            "packing_fraction",
            "packing_fraction_stat_uncertainty",
            "Packing fraction",
            "packing_fraction_by_epoch.png",
            "Packing fraction versus NH$_3$ target/polarization epoch",
        ),
    ):
        fig, ax = plt.subplots(figsize=(16, 8))

        for period in PERIODS:
            subset = plot_frame[
                plot_frame["period"] == period
            ].sort_values("epoch_number")
            if subset.empty:
                continue
            # endif

            color = period_colors[period]
            x = subset["epoch_number"].to_numpy(dtype=float)
            y = subset[quantity].to_numpy(dtype=float)
            yerr = subset[error].to_numpy(dtype=float)

            ax.errorbar(
                x,
                y,
                yerr=yerr,
                marker="s",
                linestyle="none",
                markersize=5,
                capsize=3,
                color="black",
                ecolor="black",
                label=None,
                zorder=4,
            )

            valid = (
                np.isfinite(y)
                & np.isfinite(yerr)
                & (yerr > 0.0)
            )
            if np.count_nonzero(valid) >= 1:
                weights = 1.0 / yerr[valid] ** 2
                weighted_mean = float(
                    np.sum(weights * y[valid]) / np.sum(weights)
                )
                weighted_error = math.sqrt(
                    1.0 / float(np.sum(weights))
                )
                epoch_min, epoch_max = period_epoch_ranges[period]
                ax.hlines(
                    weighted_mean,
                    epoch_min,
                    epoch_max,
                    color=color,
                    linewidth=1.6,
                    label=(
                        f"{PERIOD_LABELS[period]} weighted mean"
                    ),
                    zorder=2,
                )
                ax.fill_between(
                    [epoch_min, epoch_max],
                    weighted_mean - weighted_error,
                    weighted_mean + weighted_error,
                    color=color,
                    alpha=0.16,
                    zorder=1,
                )
            # endif

            epoch_min, epoch_max = period_epoch_ranges[period]
            ax.axvline(
                epoch_min - 0.5,
                color=color,
                linestyle="-.",
                linewidth=1.4,
                zorder=0,
            )
            ax.axvline(
                epoch_max + 0.5,
                color=color,
                linestyle="-.",
                linewidth=1.4,
                zorder=0,
            )
        # endfor

        ax.set_xlim(0.5, 24.5)
        ax.set_xticks(np.arange(1, 25))
        ax.set_xlabel("NH$_3$ epoch")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.25, linestyle="--")
        ax.legend(ncol=3, loc="best")
        fig.tight_layout()

        path = output_dir / stem
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor

    return paths


# =============================================================================
# Input loading and cut parsing
# =============================================================================

def load_one_dataset(spec: DatasetSpec) -> LoadedDataset:
    path = Path(spec.file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing ROOT input: {path}")
    # endif

    with uproot.open(path) as root_file:
        tree = resolve_tree(root_file, spec.tree_name, path)
        x_branch = resolve_branch(tree, BRANCH_ALIASES["xB"])
        tprime_branch = resolve_branch(tree, BRANCH_ALIASES["tprime"])
        mx2_branch = resolve_branch(tree, BRANCH_ALIASES["Mx2"])
        run_branch = resolve_branch(tree, BRANCH_ALIASES["runnum"])
        arrays = tree.arrays(
            [x_branch, tprime_branch, mx2_branch, run_branch],
            library="np",
        )
    # endwith

    xB = np.asarray(arrays[x_branch], dtype=np.float64)
    minus_tprime = -np.asarray(arrays[tprime_branch], dtype=np.float64)
    mx2 = np.asarray(arrays[mx2_branch], dtype=np.float64)
    runnum = np.asarray(arrays[run_branch], dtype=np.int64)
    finite = (
        np.isfinite(xB)
        & np.isfinite(minus_tprime)
        & np.isfinite(mx2)
    )

    return LoadedDataset(
        period=spec.period,
        target=spec.target,
        file_path=str(path.resolve()),
        tree_name=spec.tree_name,
        x_branch=x_branch,
        tprime_branch=tprime_branch,
        mx2_branch=mx2_branch,
        run_branch=run_branch,
        runnum=runnum[finite],
        xB=xB[finite],
        minus_tprime_gev2=minus_tprime[finite],
        mx2_gev2=mx2[finite],
    )


def load_all_datasets(
    input_paths: dict[str, dict[str, Path]],
    tree_name: str,
    workers: int,
) -> dict[tuple[str, str], LoadedDataset]:
    specs = [
        DatasetSpec(period, target, str(input_paths[period][target]), tree_name)
        for period in PERIODS
        for target in TARGETS
    ]
    loaded: dict[tuple[str, str], LoadedDataset] = {}

    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_spec = {
            executor.submit(load_one_dataset, spec): spec
            for spec in specs
        }
        for completed, future in enumerate(as_completed(future_to_spec), start=1):
            spec = future_to_spec[future]
            dataset = future.result()
            loaded[(dataset.period, dataset.target)] = dataset
            print(
                f"Loaded {completed:2d}/{len(specs)}: "
                f"{PERIOD_LABELS[dataset.period]} {dataset.target} "
                f"({dataset.mx2_gev2.size:,} finite events)"
            )
        # endfor
    # endwith

    return loaded


def load_exclusivity_cuts(path: Path) -> dict[tuple[str, int, int], CutEntry]:
    if not path.is_file():
        raise FileNotFoundError(
            "Exclusivity-cut JSON was not found.  The default assumes this "
            "script is run from RGC_enpi+/dilution_factor and channel selection "
            f"was run in RGC_enpi+/channel_selection.  Missing path: {path}"
        )
    # endif

    payload = json.loads(path.read_text(encoding="utf-8"))
    periods_payload = payload.get("periods")
    if not isinstance(periods_payload, dict):
        raise RuntimeError(
            f"Cut JSON {path} does not contain the expected 'periods' object."
        )
    # endif

    cuts: dict[tuple[str, int, int], CutEntry] = {}
    for period in PERIODS:
        entries = periods_payload.get(period)
        if not isinstance(entries, list):
            raise RuntimeError(f"Cut JSON has no list for period {period}.")
        # endif
        for item in entries:
            tight = item["tight"]
            nominal = item["nominal"]
            loose = item["loose"]
            entry = CutEntry(
                period=period,
                x_index=int(item["x_index"]),
                t_index=int(item["t_index"]),
                xB_min=float(item["xB_min"]),
                xB_max=float(item["xB_max"]),
                minus_tprime_min_gev2=float(item["minus_tprime_min_gev2"]),
                minus_tprime_max_gev2=float(item["minus_tprime_max_gev2"]),
                mu_gev2=float(item["mu_gev2"]),
                mu_error_gev2=float(item["mu_error_gev2"]),
                sigma_gev2=float(item["sigma_gev2"]),
                sigma_error_gev2=float(item["sigma_error_gev2"]),
                classification=str(item.get("classification", "unknown")),
                tight_min_gev2=float(tight[0]),
                tight_max_gev2=float(tight[1]),
                nominal_min_gev2=float(nominal[0]),
                nominal_max_gev2=float(nominal[1]),
                loose_min_gev2=float(loose[0]),
                loose_max_gev2=float(loose[1]),
            )
            cuts[(period, entry.x_index, entry.t_index)] = entry
        # endfor
    # endfor

    expected = {
        (period, x_index, t_index)
        for period in PERIODS
        for x_index in range(len(XB_BINS))
        for t_index in range(len(MINUS_TPRIME_BINS_GEV2))
    }
    missing = sorted(expected - set(cuts))
    extra = sorted(set(cuts) - expected)
    if missing or extra:
        raise RuntimeError(
            f"Cut-table coverage mismatch. Missing={missing}; extra={extra}."
        )
    # endif

    for key, entry in cuts.items():
        expected_x = XB_BINS[entry.x_index]
        expected_t = MINUS_TPRIME_BINS_GEV2[entry.t_index]
        values_match = (
            np.allclose([entry.xB_min, entry.xB_max], expected_x)
            and np.allclose(
                [entry.minus_tprime_min_gev2, entry.minus_tprime_max_gev2],
                expected_t,
            )
        )
        if not values_match:
            raise RuntimeError(
                f"Cut-table bin edges for {key} do not match the fixed analysis binning."
            )
        # endif
        for variation in CUT_VARIATIONS:
            low, high = entry.interval(variation)
            if not (math.isfinite(low) and math.isfinite(high) and low < high):
                raise RuntimeError(f"Invalid {variation} interval for {key}: {(low, high)}")
            # endif
        # endfor
    # endfor

    return cuts


# =============================================================================
# Counting and Poisson-membership representation
# =============================================================================

def membership_pattern_counts(masks: list[np.ndarray]) -> tuple[int, ...]:
    if len(masks) != len(SELECTION_NAMES):
        raise ValueError("Expected exactly four selection masks.")
    # endif
    if not masks:
        return tuple()
    # endif
    pattern = np.zeros(masks[0].size, dtype=np.uint8)
    for bit, mask in enumerate(masks):
        pattern |= np.asarray(mask, dtype=np.uint8) << bit
    # endfor
    counts = np.bincount(pattern, minlength=1 << len(masks))
    return tuple(int(value) for value in counts)


def selection_counts_from_patterns(pattern_counts: np.ndarray) -> np.ndarray:
    """Convert [...,16] membership-pattern counts to [...,4] selections."""
    bits = np.arange(16, dtype=np.uint8)
    membership = np.stack(
        [((bits >> bit) & 1).astype(np.int64) for bit in range(4)],
        axis=1,
    )
    return np.asarray(pattern_counts) @ membership


def count_one_dataset(
    dataset: LoadedDataset,
    cuts: dict[tuple[str, int, int], CutEntry],
    control_min_gev2: float,
    control_max_gev2: float,
) -> list[CountPayload]:
    rows: list[CountPayload] = []
    for x_index, (x_min, x_max) in enumerate(XB_BINS):
        x_mask = (dataset.xB >= x_min) & (dataset.xB < x_max)
        for t_index, (t_min, t_max) in enumerate(MINUS_TPRIME_BINS_GEV2):
            kinematic_mask = (
                x_mask
                & (dataset.minus_tprime_gev2 >= t_min)
                & (dataset.minus_tprime_gev2 < t_max)
            )
            values = dataset.mx2_gev2[kinematic_mask]
            entry = cuts[(dataset.period, x_index, t_index)]
            masks = [
                (values >= control_min_gev2) & (values < control_max_gev2),
            ]
            for variation in CUT_VARIATIONS:
                low, high = entry.interval(variation)
                masks.append((values >= low) & (values < high))
            # endfor
            patterns = membership_pattern_counts(masks)
            selected = selection_counts_from_patterns(
                np.asarray(patterns, dtype=np.int64)
            )
            rows.append(
                CountPayload(
                    period=dataset.period,
                    target=dataset.target,
                    x_index=x_index,
                    t_index=t_index,
                    bin_number=combined_bin_number(x_index, t_index),
                    total_kinematic_count=int(values.size),
                    selected_counts=tuple(int(value) for value in selected),
                    pattern_counts=patterns,
                )
            )
        # endfor
    # endfor
    return rows


def build_count_arrays(
    loaded: dict[tuple[str, str], LoadedDataset],
    cuts: dict[tuple[str, int, int], CutEntry],
    control_min_gev2: float,
    control_max_gev2: float,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """
    Return arrays with shape:

      pattern_counts[period,target,bin,16]
      observed_counts[period,target,bin,4]
    """
    patterns = np.zeros(
        (len(PERIODS), len(TARGETS), NUMBER_OF_BINS, 16),
        dtype=np.int64,
    )
    observed = np.zeros(
        (len(PERIODS), len(TARGETS), NUMBER_OF_BINS, 4),
        dtype=np.int64,
    )
    rows: list[dict[str, Any]] = []

    for p_index, period in enumerate(PERIODS):
        for target_index, target in enumerate(TARGETS):
            payloads = count_one_dataset(
                loaded[(period, target)],
                cuts,
                control_min_gev2,
                control_max_gev2,
            )
            for payload in payloads:
                b = payload.bin_number - 1
                patterns[p_index, target_index, b, :] = payload.pattern_counts
                observed[p_index, target_index, b, :] = payload.selected_counts
                row = {
                    "period": period,
                    "target": target,
                    "x_index": payload.x_index,
                    "t_index": payload.t_index,
                    "bin_number": payload.bin_number,
                    "xB_min": XB_BINS[payload.x_index][0],
                    "xB_max": XB_BINS[payload.x_index][1],
                    "minus_tprime_min_gev2": MINUS_TPRIME_BINS_GEV2[payload.t_index][0],
                    "minus_tprime_max_gev2": MINUS_TPRIME_BINS_GEV2[payload.t_index][1],
                    "total_kinematic_count": payload.total_kinematic_count,
                }
                for selection_index, selection in enumerate(SELECTION_NAMES):
                    row[f"count_{selection}"] = payload.selected_counts[selection_index]
                # endfor
                rows.append(row)
            # endfor
        # endfor
    # endfor

    return patterns, observed, rows


# =============================================================================
# Physics estimators
# =============================================================================

def safe_divide(numerator: np.ndarray | float, denominator: np.ndarray | float) -> np.ndarray:
    numerator_array = np.asarray(numerator, dtype=np.float64)
    denominator_array = np.asarray(denominator, dtype=np.float64)
    result_shape = np.broadcast_shapes(numerator_array.shape, denominator_array.shape)
    result = np.full(result_shape, np.nan, dtype=np.float64)
    np.divide(
        numerator_array,
        denominator_array,
        out=result,
        where=np.isfinite(denominator_array) & (denominator_array != 0.0),
    )
    return result


def charge_normalized_rates(
    selected_counts: np.ndarray,
    period: str,
    charge_fractions: dict[str, dict[str, float]],
) -> dict[str, np.ndarray]:
    return {
        target: np.asarray(selected_counts[target], dtype=np.float64)
        / float(charge_fractions[period][target])
        for target in TARGETS
    }


def method1_from_counts(
    counts_by_target: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Exact channel-selection carbon normalization."""
    nh3 = np.asarray(counts_by_target["NH3"], dtype=np.float64)
    carbon = np.asarray(counts_by_target["C"], dtype=np.float64)
    alpha = safe_divide(np.sum(nh3[..., 0], axis=-1), np.sum(carbon[..., 0], axis=-1))
    expanded_alpha = alpha[..., np.newaxis, np.newaxis]
    dilution = 1.0 - safe_divide(
        expanded_alpha * carbon[..., 1:4],
        nh3[..., 1:4],
    )
    return dilution, alpha


def method2_equation10_from_counts(
    counts: dict[str, np.ndarray],
    charge_fractions_period: dict[str, float],
) -> np.ndarray:
    """Thermal-contraction-corrected standalone dilution expression.

    This is a direct vectorized transcription of calculate_dilution_and_error()
    in calculate_dilution_factors.cpp.  The event counts remain raw counts;
    the accumulated-charge fractions enter explicitly as xA, xC, xCH, xHe,
    and xf.  Dividing every count by its charge fraction first and then applying
    a different reduced formula is not algebraically equivalent.
    """
    nA = np.asarray(counts["NH3"], dtype=np.float64)
    nC = np.asarray(counts["C"], dtype=np.float64)
    nCH = np.asarray(counts["CH2"], dtype=np.float64)
    nMT = np.asarray(counts["He"], dtype=np.float64)
    nf = np.asarray(counts["ET"], dtype=np.float64)

    xA = float(charge_fractions_period["NH3"])
    xC = float(charge_fractions_period["C"])
    xCH = float(charge_fractions_period["CH2"])
    xHe = float(charge_fractions_period["He"])
    xf = float(charge_fractions_period["ET"])

    first = -nMT * xA + nA * xHe
    second = (
        -0.579353 * nMT * xC * xCH * xf
        + (
            nf * xC * xCH
            - 3.50431 * nCH * xC * xf
            + 3.08366 * nC * xCH * xf
        ) * xHe
    )
    denominator = nA * xHe * (
        35.88 * nMT * xC * xCH * xf
        - nf * xC * xCH * xHe
        - 43.3586 * nCH * xC * xf * xHe
        + 8.47866 * nC * xCH * xf * xHe
    )
    return safe_divide(12.3729 * first * second, denominator)


def packing_fraction_equation11_from_counts(
    counts: dict[str, np.ndarray],
    charge_fractions_period: dict[str, float],
) -> np.ndarray:
    """Packing-fraction expression preserved in the supplied C++ source."""
    nA = np.asarray(counts["NH3"], dtype=np.float64)
    nC = np.asarray(counts["C"], dtype=np.float64)
    nCH = np.asarray(counts["CH2"], dtype=np.float64)
    nMT = np.asarray(counts["He"], dtype=np.float64)
    nf = np.asarray(counts["ET"], dtype=np.float64)

    xA = float(charge_fractions_period["NH3"])
    xC = float(charge_fractions_period["C"])
    xCH = float(charge_fractions_period["CH2"])
    xHe = float(charge_fractions_period["He"])
    xf = float(charge_fractions_period["ET"])

    numerator = 0.699832 * (nA / xA - nMT / xHe)
    denominator = (
        1.25055 * nCH / xCH
        - 0.23688 * nC / xC
        - 0.013668 * nf / xf
        - nMT / xHe
    )
    return safe_divide(numerator, denominator)


def integrated_packing_fraction_equation11_from_counts(
    counts: dict[str, np.ndarray],
    charge_fractions_period: dict[str, float],
) -> np.ndarray:
    integrated = {
        target: np.sum(value, axis=-2)
        for target, value in counts.items()
    }
    return packing_fraction_equation11_from_counts(
        integrated, charge_fractions_period
    )


def method3_equation14(
    counts: dict[str, np.ndarray],
    charge_fractions_period: dict[str, float],
    period_packing_fraction: np.ndarray,
) -> np.ndarray:
    """Packing-fraction dilution factor using the thermally corrected model.

    Equation (14) expresses the dilution factor in terms of a packing
    fraction rather than the measured ammonia yield in each kinematic bin.
    The supplied C++ source includes the thermally corrected packing-fraction
    relation and the corresponding direct dilution expression, but not a
    separately printed thermally corrected Eq. (14).  We therefore obtain its
    exact algebraic counterpart by solving the C++ packing-fraction equation
    for the charge-normalized ammonia rate in each bin, using the period-wide
    packing fraction, and inserting that reconstructed ammonia rate into the
    same thermally corrected dilution expression used by Method 2.

    If a bin-by-bin packing fraction were inserted, this construction would
    reproduce Method 2 algebraically.  Using the period-wide packing fraction
    makes Method 3 the intended independent target-property determination.
    """
    rates = {
        target: np.asarray(counts[target], dtype=np.float64)
        / float(charge_fractions_period[target])
        for target in TARGETS
    }

    denominator_pf = (
        1.25055 * rates["CH2"]
        - 0.23688 * rates["C"]
        - 0.013668 * rates["ET"]
        - rates["He"]
    )

    # period_packing_fraction has shape [...,3 cuts]. Insert the bin axis so
    # it broadcasts over the 24 kinematic bins.
    pf = np.asarray(period_packing_fraction, dtype=np.float64)[..., np.newaxis, :]
    reconstructed_ammonia_rate = (
        rates["He"] + safe_divide(pf * denominator_pf, 0.699832)
    )

    reconstructed_counts = dict(counts)
    reconstructed_counts["NH3"] = (
        reconstructed_ammonia_rate * float(charge_fractions_period["NH3"])
    )
    return method2_equation10_from_counts(
        reconstructed_counts, charge_fractions_period
    )

def observed_estimators_for_period(
    observed_period: np.ndarray,
    period: str,
    charge_fractions: dict[str, dict[str, float]],
) -> dict[str, Any]:
    counts = {
        target: observed_period[target_index]
        for target_index, target in enumerate(TARGETS)
    }
    method1, alpha = method1_from_counts(counts)
    selected_counts = {
        target: values[:, 1:4]
        for target, values in counts.items()
    }
    rates = charge_normalized_rates(selected_counts, period, charge_fractions)
    method2 = method2_equation10_from_counts(
        selected_counts, charge_fractions[period]
    )
    pf_bin = packing_fraction_equation11_from_counts(
        selected_counts, charge_fractions[period]
    )
    pf_period = integrated_packing_fraction_equation11_from_counts(
        selected_counts, charge_fractions[period]
    )
    method3 = method3_equation14(
        selected_counts, charge_fractions[period], pf_period
    )
    return {
        "method1": method1,
        "method1_alpha": alpha,
        "method2": method2,
        "packing_fraction_bin": pf_bin,
        "packing_fraction_period": pf_period,
        "method3": method3,
        "rates": rates,
    }


# =============================================================================
# Statistical replicas
# =============================================================================

def bootstrap_period_worker(
    period_index: int,
    pattern_counts_period: np.ndarray,
    charge_fractions_period: dict[str, float],
    replicas: int,
    seed: int,
) -> dict[str, Any]:
    """Generate correlated Poisson replicas for one run period in bounded-memory batches."""
    rng = np.random.default_rng(seed)
    batch_size = min(250, replicas)

    method1_output = np.full((replicas, NUMBER_OF_BINS, 3), np.nan, dtype=np.float64)
    method2_output = np.full_like(method1_output, np.nan)
    method3_output = np.full_like(method1_output, np.nan)
    pf_bin_output = np.full_like(method1_output, np.nan)
    alpha_output = np.full(replicas, np.nan, dtype=np.float64)
    pf_period_output = np.full((replicas, 3), np.nan, dtype=np.float64)

    for start in range(0, replicas, batch_size):
        stop = min(start + batch_size, replicas)
        number_in_batch = stop - start
        # Shape: batch,target,bin,16.  Batching avoids a multi-gigabyte
        # allocation when the default 10,000 replicas are requested.
        pattern_replicas = rng.poisson(
            lam=pattern_counts_period[np.newaxis, ...],
            size=(number_in_batch,) + pattern_counts_period.shape,
        )
        selected = selection_counts_from_patterns(pattern_replicas)
        counts = {
            target: selected[:, target_index, :, :]
            for target_index, target in enumerate(TARGETS)
        }

        method1, alpha = method1_from_counts(counts)
        selected_for_auxiliary = {
            target: value[..., 1:4]
            for target, value in counts.items()
        }
        rates = {
            target: selected_for_auxiliary[target] / charge_fractions_period[target]
            for target in TARGETS
        }
        method2 = method2_equation10_from_counts(
            selected_for_auxiliary, charge_fractions_period
        )
        pf_bin = packing_fraction_equation11_from_counts(
            selected_for_auxiliary, charge_fractions_period
        )
        pf_period = integrated_packing_fraction_equation11_from_counts(
            selected_for_auxiliary, charge_fractions_period
        )
        method3 = method3_equation14(
            selected_for_auxiliary, charge_fractions_period, pf_period
        )

        method1_output[start:stop] = method1
        method2_output[start:stop] = method2
        method3_output[start:stop] = method3
        pf_bin_output[start:stop] = pf_bin
        alpha_output[start:stop] = alpha
        pf_period_output[start:stop] = pf_period
    # endfor

    return {
        "period_index": period_index,
        "method1": method1_output,
        "method1_alpha": alpha_output,
        "method2": method2_output,
        "packing_fraction_bin": pf_bin_output,
        "packing_fraction_period": pf_period_output,
        "method3": method3_output,
    }


def run_bootstrap(
    pattern_counts: np.ndarray,
    charge_fractions: dict[str, dict[str, float]],
    replicas: int,
    seed: int,
    workers: int,
) -> dict[str, dict[str, np.ndarray]]:
    output: dict[str, dict[str, np.ndarray]] = {}
    with ProcessPoolExecutor(max_workers=min(workers, len(PERIODS))) as executor:
        futures = {
            executor.submit(
                bootstrap_period_worker,
                period_index,
                pattern_counts[period_index],
                charge_fractions[period],
                replicas,
                seed + 100003 * period_index,
            ): period
            for period_index, period in enumerate(PERIODS)
        }
        for future in as_completed(futures):
            period = futures[future]
            result = future.result()
            output[period] = {
                key: value
                for key, value in result.items()
                if key != "period_index"
            }
            print(f"Completed {replicas:,} statistical replicas for {PERIOD_LABELS[period]}.")
        # endfor
    # endwith
    return output


def summarize_replicas(
    central: np.ndarray,
    replicas: np.ndarray,
) -> dict[str, np.ndarray]:
    finite = np.isfinite(replicas)
    valid_fraction = np.mean(finite, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        standard_deviation = np.nanstd(replicas, axis=0, ddof=1)
        p16 = np.nanpercentile(replicas, 16.0, axis=0)
        p50 = np.nanpercentile(replicas, 50.0, axis=0)
        p84 = np.nanpercentile(replicas, 84.0, axis=0)
    # endwith
    return {
        "central": central,
        "stat_uncertainty": standard_deviation,
        "p16": p16,
        "p50": p50,
        "p84": p84,
        "valid_replica_fraction": valid_fraction,
    }


def pairwise_covariance(samples: np.ndarray) -> np.ndarray:
    """Pairwise-complete covariance for samples shaped [replica,bin]."""
    n_bins = samples.shape[1]
    covariance = np.full((n_bins, n_bins), np.nan, dtype=np.float64)
    for i in range(n_bins):
        for j in range(i, n_bins):
            valid = np.isfinite(samples[:, i]) & np.isfinite(samples[:, j])
            if np.count_nonzero(valid) >= 2:
                value = np.cov(samples[valid, i], samples[valid, j], ddof=1)[0, 1]
                covariance[i, j] = value
                covariance[j, i] = value
            # endif
        # endfor
    # endfor
    return covariance


def covariance_to_correlation(covariance: np.ndarray) -> np.ndarray:
    sigma = np.sqrt(np.diag(covariance))
    denominator = np.outer(sigma, sigma)
    return safe_divide(covariance, denominator)


# =============================================================================
# Output assembly
# =============================================================================

def quality_flags(value: float, valid_fraction: float) -> list[str]:
    flags: list[str] = []
    if not math.isfinite(value):
        flags.append("nonfinite_central_value")
    else:
        if value < 0.0:
            flags.append("below_zero")
        # endif
        if value > 1.0:
            flags.append("above_one")
        # endif
    # endif
    if not math.isfinite(valid_fraction) or valid_fraction < 0.99:
        flags.append("less_than_99_percent_valid_replicas")
    # endif
    return flags


def scalar_record(summary: dict[str, np.ndarray], index: tuple[int, ...]) -> dict[str, Any]:
    value = float(summary["central"][index])
    valid_fraction = float(summary["valid_replica_fraction"][index])
    return {
        "value": value,
        "stat_uncertainty": float(summary["stat_uncertainty"][index]),
        "bootstrap_p16": float(summary["p16"][index]),
        "bootstrap_p50": float(summary["p50"][index]),
        "bootstrap_p84": float(summary["p84"][index]),
        "valid_replica_fraction": valid_fraction,
        "flags": quality_flags(value, valid_fraction),
    }



def relative_method_half_difference(
    method1: np.ndarray,
    method2: np.ndarray,
) -> np.ndarray:
    """
    Return |f1-f2|/(f1+f2), the relative half-difference about their average.

    Since f_rec=(f1+f2)/2 and delta=(|f1-f2|)/2, delta/f_rec is exactly
    |f1-f2|/(f1+f2).
    """
    denominator = method1 + method2
    result = np.full_like(denominator, np.nan, dtype=float)
    valid = (
        np.isfinite(method1)
        & np.isfinite(method2)
        & np.isfinite(denominator)
        & (denominator > 0.0)
    )
    result[valid] = np.abs(method1[valid] - method2[valid]) / denominator[valid]
    return result


def fit_constant_scale(
    values: np.ndarray,
    uncertainties: np.ndarray,
) -> dict[str, float]:
    """Fit a constant to finite positive-uncertainty values."""
    values = np.asarray(values, dtype=float)
    uncertainties = np.asarray(uncertainties, dtype=float)
    valid = (
        np.isfinite(values)
        & np.isfinite(uncertainties)
        & (uncertainties > 0.0)
    )
    if np.count_nonzero(valid) < 2:
        finite_values = values[np.isfinite(values)]
        value = (
            float(np.nanmean(finite_values))
            if finite_values.size
            else math.nan
        )
        return {
            "value": value,
            "uncertainty": math.nan,
            "chi2": math.nan,
            "ndf": 0,
            "chi2_ndf": math.nan,
            "number_of_bins": int(np.count_nonzero(valid)),
        }
    # endif

    weights = 1.0 / uncertainties[valid] ** 2
    weight_sum = float(np.sum(weights))
    fitted = float(np.sum(weights * values[valid]) / weight_sum)
    fitted_uncertainty = math.sqrt(1.0 / weight_sum)
    chi2 = float(
        np.sum(
            ((values[valid] - fitted) / uncertainties[valid]) ** 2
        )
    )
    ndf = int(np.count_nonzero(valid) - 1)
    return {
        "value": fitted,
        "uncertainty": fitted_uncertainty,
        "chi2": chi2,
        "ndf": ndf,
        "chi2_ndf": chi2 / ndf if ndf > 0 else math.nan,
        "number_of_bins": int(np.count_nonzero(valid)),
    }


def build_method_scale_summary(
    central_method1: np.ndarray,
    central_method2: np.ndarray,
    replica_method1: np.ndarray,
    replica_method2: np.ndarray,
) -> dict[str, Any]:
    """
    Build per-bin relative half-differences and period-wide constant fits.

    Arrays have shape (bin, cut) for central values and (replica, bin, cut)
    for replicas.
    """
    central_relative = relative_method_half_difference(
        central_method1,
        central_method2,
    )
    replica_relative = relative_method_half_difference(
        replica_method1,
        replica_method2,
    )
    replica_uncertainty = np.nanstd(
        replica_relative,
        axis=0,
        ddof=1,
    )

    fits: list[dict[str, float]] = []
    for cut_index, variation in enumerate(CUT_VARIATIONS):
        fit = fit_constant_scale(
            central_relative[:, cut_index],
            replica_uncertainty[:, cut_index],
        )
        fit["cut_variation"] = variation
        fits.append(fit)
    # endfor

    return {
        "relative_half_difference": central_relative,
        "relative_half_difference_stat_uncertainty": replica_uncertainty,
        "constant_fits": fits,
    }


def build_output_tables(
    observed: np.ndarray,
    central_results: dict[str, dict[str, Any]],
    bootstrap_results: dict[str, dict[str, np.ndarray]],
    charge_fractions: dict[str, dict[str, float]],
    cuts: dict[tuple[str, int, int], CutEntry],
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any], dict[str, Any]]:
    flat_rows: list[dict[str, Any]] = []
    master_periods: dict[str, Any] = {}
    compact_periods: dict[str, Any] = {}
    summaries_by_period: dict[str, Any] = {}

    method_scale_summaries = {
        period: build_method_scale_summary(
            central_method1=central_results[period]["method1"],
            central_method2=central_results[period]["method2"],
            replica_method1=bootstrap_results[period]["method1"],
            replica_method2=bootstrap_results[period]["method2"],
        )
        for period in PERIODS
    }
    global_method_scale_by_cut: list[dict[str, Any]] = []
    for cut_index, variation in enumerate(CUT_VARIATIONS):
        period_values = np.asarray(
            [
                method_scale_summaries[period]["constant_fits"][cut_index]["value"]
                for period in PERIODS
            ],
            dtype=float,
        )
        global_method_scale_by_cut.append(
            {
                "cut_variation": variation,
                "fraction": float(np.nanmean(period_values)),
                "percent": 100.0 * float(np.nanmean(period_values)),
                "period_scale_fractions": {
                    period: float(
                        method_scale_summaries[period]["constant_fits"][cut_index]["value"]
                    )
                    for period in PERIODS
                },
                "definition": (
                    "Arithmetic mean of the three period-specific constant-fit "
                    "Method-1/Method-2 relative half-difference scales."
                ),
                "correlation_scope": (
                    "One global multiplicative scale fully correlated across "
                    "all periods and kinematic bins."
                ),
            }
        )
    # endfor

    for p_index, period in enumerate(PERIODS):
        central = central_results[period]
        replicas = bootstrap_results[period]
        recommended_central = 0.5 * (central["method1"] + central["method2"])
        recommended_replicas = 0.5 * (replicas["method1"] + replicas["method2"])
        method_scale_summary = method_scale_summaries[period]
        summaries = {
            "method1": summarize_replicas(central["method1"], replicas["method1"]),
            "method2": summarize_replicas(central["method2"], replicas["method2"]),
            "method3": summarize_replicas(central["method3"], replicas["method3"]),
            "recommended": summarize_replicas(
                recommended_central, recommended_replicas
            ),
            "packing_fraction_bin": summarize_replicas(
                central["packing_fraction_bin"], replicas["packing_fraction_bin"]
            ),
            "packing_fraction_period": summarize_replicas(
                central["packing_fraction_period"], replicas["packing_fraction_period"]
            ),
        }
        alpha_summary = summarize_replicas(
            np.asarray(central["method1_alpha"]),
            np.asarray(replicas["method1_alpha"]),
        )
        summaries["method1_alpha"] = alpha_summary
        summaries["method_scale"] = method_scale_summary
        summaries_by_period[period] = summaries

        period_bins: list[dict[str, Any]] = []
        compact_bins: list[dict[str, Any]] = []
        for b in range(NUMBER_OF_BINS):
            x_index = b // len(MINUS_TPRIME_BINS_GEV2)
            t_index = b % len(MINUS_TPRIME_BINS_GEV2)
            cut_entry = cuts[(period, x_index, t_index)]
            cut_payload: dict[str, Any] = {}
            compact_cut_payload: dict[str, Any] = {}
            for cut_index, variation in enumerate(CUT_VARIATIONS):
                selection_index = cut_index + 1
                counts_payload = {
                    target: int(observed[p_index, target_index, b, selection_index])
                    for target_index, target in enumerate(TARGETS)
                }
                rates_payload = {
                    target: float(central["rates"][target][b, cut_index])
                    for target in TARGETS
                }
                method1_record = scalar_record(summaries["method1"], (b, cut_index))
                method2_record = scalar_record(summaries["method2"], (b, cut_index))
                method3_record = scalar_record(summaries["method3"], (b, cut_index))
                recommended_record = scalar_record(
                    summaries["recommended"], (b, cut_index)
                )
                scale_fit = method_scale_summary["constant_fits"][cut_index]
                global_scale = global_method_scale_by_cut[cut_index]
                scale_fraction = float(global_scale["fraction"])
                local_relative_half_difference = float(
                    method_scale_summary["relative_half_difference"][
                        b,
                        cut_index,
                    ]
                )
                local_relative_half_difference_uncertainty = float(
                    method_scale_summary[
                        "relative_half_difference_stat_uncertainty"
                    ][b, cut_index]
                )
                recommended_record["dilution_model_scale_fraction"] = (
                    scale_fraction
                )
                recommended_record["dilution_model_scale_percent"] = (
                    100.0 * scale_fraction
                )
                recommended_record["dilution_model_scale_uncertainty"] = math.nan
                recommended_record["dilution_model_scale_chi2"] = math.nan
                recommended_record[
                    "period_specific_scale_fit_chi2_diagnostic"
                ] = float(scale_fit["chi2"])
                recommended_record["dilution_model_scale_ndf"] = 0
                recommended_record["dilution_model_scale_chi2_ndf"] = math.nan
                recommended_record[
                    "period_specific_scale_fit_chi2_ndf_diagnostic"
                ] = float(scale_fit["chi2_ndf"])
                recommended_record["dilution_model_systematic"] = (
                    scale_fraction * recommended_record["value"]
                )
                recommended_record[
                    "local_method_relative_half_difference"
                ] = local_relative_half_difference
                recommended_record[
                    "local_method_relative_half_difference_stat_uncertainty"
                ] = local_relative_half_difference_uncertainty
                recommended_record["total_uncertainty_quadrature"] = float(
                    math.hypot(
                        recommended_record["stat_uncertainty"],
                        recommended_record["dilution_model_systematic"],
                    )
                )
                recommended_record["definition"] = (
                    "Average of Method 1 and Method 2. The dilution-model "
                    "uncertainty is one global multiplicative scale, equal to "
                    "the arithmetic mean of the three period-specific constant "
                    "fits to |f1-f2|/(f1+f2). It is fully correlated across all "
                    "periods and kinematic bins."
                )
                pf_bin_record = scalar_record(
                    summaries["packing_fraction_bin"], (b, cut_index)
                )
                pf_period_record = scalar_record(
                    summaries["packing_fraction_period"], (cut_index,)
                )
                low, high = cut_entry.interval(variation)
                cut_payload[variation] = {
                    "mx2_min_gev2": low,
                    "mx2_max_gev2": high,
                    "counts": counts_payload,
                    "charge_normalized_rates_relative_units": rates_payload,
                    "method1_carbon_subtraction": method1_record,
                    "method2_equation10": method2_record,
                    "packing_fraction_equation11_bin_diagnostic": pf_bin_record,
                    "packing_fraction_equation11_period_integrated": pf_period_record,
                    "method3_packing_fraction_constrained": method3_record,
                    "recommended_method1_method2_average": recommended_record,
                }
                compact_cut_payload[variation] = {
                    "method1": method1_record,
                    "method2": method2_record,
                    "method3": method3_record,
                    "recommended": recommended_record,
                    "packing_fraction_used_by_method3": pf_period_record,
                }

                row: dict[str, Any] = {
                    "period": period,
                    "period_label": PERIOD_LABELS[period],
                    "x_index": x_index,
                    "t_index": t_index,
                    "bin_number": b + 1,
                    "xB_min": XB_BINS[x_index][0],
                    "xB_max": XB_BINS[x_index][1],
                    "minus_tprime_min_gev2": MINUS_TPRIME_BINS_GEV2[t_index][0],
                    "minus_tprime_max_gev2": MINUS_TPRIME_BINS_GEV2[t_index][1],
                    "cut_variation": variation,
                    "mx2_min_gev2": low,
                    "mx2_max_gev2": high,
                    "method1_carbon_scale_period": float(alpha_summary["central"]),
                    "method1_carbon_scale_stat_uncertainty": float(
                        alpha_summary["stat_uncertainty"]
                    ),
                }
                for target in TARGETS:
                    row[f"count_{target}"] = counts_payload[target]
                    row[f"relative_charge_fraction_{target}"] = charge_fractions[period][target]
                    row[f"rate_{target}_relative_units"] = rates_payload[target]
                # endfor
                for label, record in (
                    ("method1", method1_record),
                    ("method2", method2_record),
                    ("method3", method3_record),
                    ("recommended", recommended_record),
                    ("packing_fraction_bin", pf_bin_record),
                    ("packing_fraction_period", pf_period_record),
                ):
                    row[f"{label}_value"] = record["value"]
                    row[f"{label}_stat_uncertainty"] = record["stat_uncertainty"]
                    row[f"{label}_bootstrap_p16"] = record["bootstrap_p16"]
                    row[f"{label}_bootstrap_p50"] = record["bootstrap_p50"]
                    row[f"{label}_bootstrap_p84"] = record["bootstrap_p84"]
                    row[f"{label}_valid_replica_fraction"] = record[
                        "valid_replica_fraction"
                    ]
                    row[f"{label}_flags"] = ";".join(record["flags"])
                # endfor
                row[
                    "recommended_dilution_model_scale_fraction"
                ] = recommended_record[
                    "dilution_model_scale_fraction"
                ]
                row[
                    "recommended_dilution_model_scale_percent"
                ] = recommended_record[
                    "dilution_model_scale_percent"
                ]
                row[
                    "recommended_dilution_model_systematic_absolute"
                ] = recommended_record[
                    "dilution_model_systematic"
                ]
                row[
                    "recommended_dilution_model_scale_chi2_ndf"
                ] = recommended_record[
                    "dilution_model_scale_chi2_ndf"
                ]
                row[
                    "recommended_local_method_relative_half_difference"
                ] = recommended_record[
                    "local_method_relative_half_difference"
                ]
                row[
                    "recommended_local_method_relative_half_difference_stat_uncertainty"
                ] = recommended_record[
                    "local_method_relative_half_difference_stat_uncertainty"
                ]
                row["recommended_dilution_model_systematic"] = recommended_record[
                    "dilution_model_systematic"
                ]
                row["recommended_total_uncertainty_quadrature"] = recommended_record[
                    "total_uncertainty_quadrature"
                ]
                flat_rows.append(row)
            # endfor

            bin_record = {
                "x_index": x_index,
                "t_index": t_index,
                "bin_number": b + 1,
                "xB_min": XB_BINS[x_index][0],
                "xB_max": XB_BINS[x_index][1],
                "minus_tprime_min_gev2": MINUS_TPRIME_BINS_GEV2[t_index][0],
                "minus_tprime_max_gev2": MINUS_TPRIME_BINS_GEV2[t_index][1],
                "exclusivity_fit": {
                    "mu_gev2": cut_entry.mu_gev2,
                    "mu_error_gev2": cut_entry.mu_error_gev2,
                    "sigma_gev2": cut_entry.sigma_gev2,
                    "sigma_error_gev2": cut_entry.sigma_error_gev2,
                    "classification": cut_entry.classification,
                },
                "cuts": cut_payload,
            }
            period_bins.append(bin_record)
            compact_bins.append(
                {
                    "x_index": x_index,
                    "t_index": t_index,
                    "bin_number": b + 1,
                    "xB_min": XB_BINS[x_index][0],
                    "xB_max": XB_BINS[x_index][1],
                    "minus_tprime_min_gev2": MINUS_TPRIME_BINS_GEV2[t_index][0],
                    "minus_tprime_max_gev2": MINUS_TPRIME_BINS_GEV2[t_index][1],
                    "cuts": compact_cut_payload,
                }
            )
        # endfor

        method_scale_by_cut = {
            global_fit["cut_variation"]: {
                "fraction": float(global_fit["fraction"]),
                "percent": float(global_fit["percent"]),
                "period_specific_scale_fraction_diagnostic": float(
                    method_scale_summary["constant_fits"][cut_index]["value"]
                ),
                "period_specific_scale_percent_diagnostic": 100.0 * float(
                    method_scale_summary["constant_fits"][cut_index]["value"]
                ),
                "correlation_scope": global_fit["correlation_scope"],
                "definition": global_fit["definition"],
            }
            for cut_index, global_fit in enumerate(global_method_scale_by_cut)
        }
        master_periods[period] = {
            "charge_fractions": charge_fractions[period],
            "method1_period_carbon_scale": scalar_record(alpha_summary, ()),
            "dilution_model_scale_by_cut": method_scale_by_cut,
            "period_packing_fraction_by_cut": {
                variation: scalar_record(
                    summaries["packing_fraction_period"], (cut_index,)
                )
                for cut_index, variation in enumerate(CUT_VARIATIONS)
            },
            "bins": period_bins,
        }
        compact_periods[period] = {
            "charge_fractions": charge_fractions[period],
            "recommended_nominal_method": "average_of_method1_and_method2",
            "dilution_model_scale_by_cut": method_scale_by_cut,
            "period_packing_fraction_by_cut": master_periods[period][
                "period_packing_fraction_by_cut"
            ],
            "bins": compact_bins,
        }
    # endfor

    summaries_by_period["_global_method_scale_by_cut"] = (
        global_method_scale_by_cut
    )
    frame = pd.DataFrame(flat_rows)
    return frame, master_periods, compact_periods, summaries_by_period


def write_covariance_products(
    output_dir: Path,
    bootstrap_results: dict[str, dict[str, np.ndarray]],
) -> dict[str, Any]:
    ensure_directory(output_dir)
    manifest: dict[str, Any] = {}
    for period in PERIODS:
        manifest[period] = {}
        for method_key in (
            "method1", "method2", "method3", "recommended", "packing_fraction_bin"
        ):
            manifest[period][method_key] = {}
            if method_key == "recommended":
                samples = 0.5 * (
                    bootstrap_results[period]["method1"]
                    + bootstrap_results[period]["method2"]
                )
            else:
                samples = bootstrap_results[period][method_key]
            for cut_index, variation in enumerate(CUT_VARIATIONS):
                cut_samples = samples[:, :, cut_index]
                covariance = pairwise_covariance(cut_samples)
                correlation = covariance_to_correlation(covariance)
                stem = f"{method_key}_{period}_{variation}_v12"
                cov_path = output_dir / f"{stem}_covariance.json"
                corr_path = output_dir / f"{stem}_correlation.json"
                write_json(
                    cov_path,
                    {
                        "period": period,
                        "quantity": method_key,
                        "cut_variation": variation,
                        "bin_order": list(range(1, NUMBER_OF_BINS + 1)),
                        "matrix": covariance,
                    },
                )
                write_json(
                    corr_path,
                    {
                        "period": period,
                        "quantity": method_key,
                        "cut_variation": variation,
                        "bin_order": list(range(1, NUMBER_OF_BINS + 1)),
                        "matrix": correlation,
                    },
                )
                manifest[period][method_key][variation] = {
                    "covariance": str(cov_path),
                    "correlation": str(corr_path),
                }
            # endfor
        # endfor
    # endfor
    return manifest


# =============================================================================
# Plotting
# =============================================================================

def plot_three_method_comparison(
    output_dir: Path,
    period: str,
    summaries: dict[str, Any],
) -> str:
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    offsets = {"method1": -0.18, "method2": 0.0, "method3": 0.18}
    labels = {
        "method1": "Method 1: carbon subtraction",
        "method2": "Method 2: direct auxiliary-target subtraction",
        "method3": "Method 3: packing-fraction constrained subtraction",
    }

    fig, axes = plt.subplots(3, 1, figsize=(17, 14), sharex=True)
    for cut_index, variation in enumerate(CUT_VARIATIONS):
        ax = axes[cut_index]
        for method_key in ("method1", "method2", "method3"):
            ax.errorbar(
                bins + offsets[method_key],
                summaries[method_key]["central"][:, cut_index],
                yerr=summaries[method_key]["stat_uncertainty"][:, cut_index],
                marker="o",
                linestyle="none",
                markersize=4,
                capsize=2,
                label=labels[method_key],
            )
        ax.set_ylim(0.1, 0.6)
        ax.set_ylabel("Dilution factor")
        ax.set_title(variation.capitalize())
        ax.grid(alpha=0.25)
        ax.legend(ncol=3)
    axes[-1].set_xlabel("Combined kinematic-bin number")
    axes[-1].set_xticks(bins)
    fig.suptitle(f"{PERIOD_LABELS[period]} dilution-factor method comparison")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    path = output_dir / f"three_method_comparison_{period}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_three_period_comparison(
    output_dir: Path,
    method_key: str,
    summaries_by_period: dict[str, Any],
) -> str:
    """Compare Su22, Fa22, and Sp23 for one dilution-factor method."""
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    offsets = {"su22": -0.18, "fa22": 0.0, "sp23": 0.18}
    method_titles = {
        "method1": "Method 1: carbon subtraction",
        "method2": "Method 2: direct auxiliary-target subtraction",
        "method3": "Method 3: packing-fraction constrained subtraction",
    }

    fig, axes = plt.subplots(3, 1, figsize=(17, 14), sharex=True)
    for cut_index, variation in enumerate(CUT_VARIATIONS):
        ax = axes[cut_index]
        for period in PERIODS:
            summary = summaries_by_period[period][method_key]
            ax.errorbar(
                bins + offsets[period],
                summary["central"][:, cut_index],
                yerr=summary["stat_uncertainty"][:, cut_index],
                marker="o",
                linestyle="none",
                markersize=4,
                capsize=2,
                label=PERIOD_LABELS[period],
            )
        ax.set_ylim(0.1, 0.6)
        ax.set_ylabel("Dilution factor")
        ax.set_title(variation.capitalize())
        ax.grid(alpha=0.25)
        ax.legend(ncol=3)
    axes[-1].set_xlabel("Combined kinematic-bin number")
    axes[-1].set_xticks(bins)
    fig.suptitle(f"Run-period comparison — {method_titles[method_key]}")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    path = output_dir / f"three_period_comparison_{method_key}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_packing_fraction_summary(
    output_dir: Path,
    period: str,
    summaries: dict[str, Any],
) -> str:
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    fig, axes = plt.subplots(3, 1, figsize=(17, 14), sharex=True)
    for cut_index, variation in enumerate(CUT_VARIATIONS):
        ax = axes[cut_index]
        central = summaries["packing_fraction_bin"]["central"][:, cut_index]
        uncertainty = summaries["packing_fraction_bin"]["stat_uncertainty"][:, cut_index]
        period_value = summaries["packing_fraction_period"]["central"][cut_index]
        period_uncertainty = summaries["packing_fraction_period"]["stat_uncertainty"][cut_index]
        ax.errorbar(
            bins,
            central,
            yerr=uncertainty,
            marker="o",
            linestyle="none",
            markersize=4,
            capsize=2,
            label="Per-bin auxiliary-target packing fraction",
        )
        ax.axhline(period_value, linewidth=1.2, label="Period-integrated packing fraction")
        ax.axhspan(
            period_value - period_uncertainty,
            period_value + period_uncertainty,
            alpha=0.15,
        )
        ax.set_ylabel("Packing fraction")
        ax.set_title(variation.capitalize())
        ax.grid(alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel("Combined kinematic-bin number")
    axes[-1].set_xticks(bins)
    fig.suptitle(f"{PERIOD_LABELS[period]} packing-fraction diagnostics")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    path = output_dir / f"packing_fraction_summary_{period}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_method1_control_summary(
    output_dir: Path,
    central_results: dict[str, dict[str, Any]],
    summaries_by_period: dict[str, Any],
) -> str:
    ensure_directory(output_dir)
    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(PERIODS))
    values = np.array([central_results[p]["method1_alpha"] for p in PERIODS], dtype=float)
    errors = np.array(
        [summaries_by_period[p]["method1_alpha"]["stat_uncertainty"] for p in PERIODS],
        dtype=float,
    )
    ax.errorbar(x, values, yerr=errors, marker="o", linestyle="none", capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels([PERIOD_LABELS[p] for p in PERIODS])
    ax.set_ylabel("Raw-count NH$_3$/C normalization")
    ax.set_title(r"Method-1 period-wide carbon normalization" "\n" r"$0.00 \leq M_x^2 < 0.40$ GeV$^2$")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = output_dir / "method1_period_carbon_scales.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)



def plot_nominal_method_comparison(
    output_dir: Path,
    period: str,
    summaries: dict[str, Any],
) -> str:
    """Nominal-cut comparison of the three dilution-factor methods."""
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    offsets = {"method1": -0.18, "method2": 0.0, "method3": 0.18}
    labels = {
        "method1": "Method 1: carbon-template subtraction",
        "method2": "Method 2: direct auxiliary-target subtraction",
        "method3": "Method 3: packing-fraction constrained subtraction",
    }
    nominal_index = CUT_VARIATIONS.index("nominal")
    fig, ax = plt.subplots(figsize=(17, 6.5))
    for method_key in ("method1", "method2", "method3"):
        ax.errorbar(
            bins + offsets[method_key],
            summaries[method_key]["central"][:, nominal_index],
            yerr=summaries[method_key]["stat_uncertainty"][:, nominal_index],
            marker="o",
            linestyle="none",
            markersize=4,
            capsize=2,
            label=labels[method_key],
        )
    ax.set_ylim(0.1, 0.6)
    ax.set_xlabel("Combined kinematic-bin number")
    ax.set_ylabel("Dilution factor")
    ax.set_xticks(bins)
    ax.grid(alpha=0.25)
    ax.legend(ncol=3)
    ax.set_title(f"{PERIOD_LABELS[period]} nominal-cut dilution-factor comparison")
    fig.tight_layout()
    path = output_dir / f"nominal_method_comparison_{period}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_nominal_period_comparison(
    output_dir: Path,
    method_key: str,
    summaries_by_period: dict[str, Any],
) -> str:
    """Nominal-cut comparison of all periods for one method."""
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    offsets = {"su22": -0.18, "fa22": 0.0, "sp23": 0.18}
    method_titles = {
        "method1": "carbon-template subtraction",
        "method2": "direct auxiliary-target subtraction",
        "method3": "packing-fraction constrained subtraction",
    }
    nominal_index = CUT_VARIATIONS.index("nominal")
    fig, ax = plt.subplots(figsize=(17, 6.5))
    for period in PERIODS:
        summary = summaries_by_period[period][method_key]
        ax.errorbar(
            bins + offsets[period],
            summary["central"][:, nominal_index],
            yerr=summary["stat_uncertainty"][:, nominal_index],
            marker="o",
            linestyle="none",
            markersize=4,
            capsize=2,
            label=PERIOD_LABELS[period],
        )
    ax.set_ylim(0.1, 0.6)
    ax.set_xlabel("Combined kinematic-bin number")
    ax.set_ylabel("Dilution factor")
    ax.set_xticks(bins)
    ax.grid(alpha=0.25)
    ax.legend(ncol=3)
    ax.set_title(
        "Nominal-cut run-period comparison — " + method_titles[method_key]
    )
    fig.tight_layout()
    path = output_dir / f"nominal_period_comparison_{method_key}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_nominal_packing_fraction_period_comparison(
    output_dir: Path,
    summaries_by_period: dict[str, Any],
) -> str:
    """Compare nominal per-bin and period-integrated packing fractions."""
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    offsets = {"su22": -0.18, "fa22": 0.0, "sp23": 0.18}
    nominal_index = CUT_VARIATIONS.index("nominal")
    fig, ax = plt.subplots(figsize=(17, 6.5))
    for period in PERIODS:
        per_bin = summaries_by_period[period]["packing_fraction_bin"]
        period_pf = summaries_by_period[period]["packing_fraction_period"]
        ax.errorbar(
            bins + offsets[period],
            per_bin["central"][:, nominal_index],
            yerr=per_bin["stat_uncertainty"][:, nominal_index],
            marker="o",
            linestyle="none",
            markersize=4,
            capsize=2,
            label=f"{PERIOD_LABELS[period]} per-bin",
        )
        value = period_pf["central"][nominal_index]
        error = period_pf["stat_uncertainty"][nominal_index]
        ax.axhline(value, linewidth=1.0, linestyle="--")
        ax.axhspan(value - error, value + error, alpha=0.08)
    ax.set_xlabel("Combined kinematic-bin number")
    ax.set_ylabel("Packing fraction")
    ax.set_xticks(bins)
    ax.grid(alpha=0.25)
    ax.legend(ncol=3)
    ax.set_title("Nominal-cut packing-fraction comparison by run period")
    fig.tight_layout()
    path = output_dir / "nominal_packing_fraction_period_comparison.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_nominal_method_differences(
    output_dir: Path,
    summaries_by_period: dict[str, Any],
    bootstrap_results: dict[str, dict[str, np.ndarray]],
) -> str:
    """
    Plot the Method-1/Method-2 relative half-difference and scale fit.

    The plotted quantity is 100*|f1-f2|/(f1+f2), equal to the fractional
    half-difference about the recommended average.
    """
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    nominal_index = CUT_VARIATIONS.index("nominal")
    fig, axes = plt.subplots(3, 1, figsize=(17, 13), sharex=True)
    global_scale = summaries_by_period[
        "_global_method_scale_by_cut"
    ][nominal_index]
    global_percent = float(global_scale["percent"])

    for ax, period in zip(axes, PERIODS):
        scale_summary = summaries_by_period[period]["method_scale"]
        values = (
            100.0
            * scale_summary["relative_half_difference"][:, nominal_index]
        )
        errors = (
            100.0
            * scale_summary[
                "relative_half_difference_stat_uncertainty"
            ][:, nominal_index]
        )
        fit = scale_summary["constant_fits"][nominal_index]
        fit_percent = 100.0 * float(fit["value"])
        fit_error_percent = 100.0 * float(fit["uncertainty"])

        ax.errorbar(
            bins,
            values,
            yerr=errors,
            marker="o",
            linestyle="none",
            markersize=4,
            capsize=2,
            label=(
                r"$100\,|f_1-f_2|/(f_1+f_2)$"
            ),
        )
        ax.axhline(
            fit_percent,
            linewidth=1.2,
            label=f"Period constant fit: {fit_percent:.2f}%",
        )
        ax.axhline(
            global_percent,
            linewidth=1.0,
            linestyle="--",
            label=f"Adopted global scale: {global_percent:.2f}%",
        )
        if np.isfinite(fit_error_percent):
            ax.axhspan(
                fit_percent - fit_error_percent,
                fit_percent + fit_error_percent,
                alpha=0.16,
                label="Fit uncertainty",
            )
        # endif
        ax.set_ylabel("Relative half-difference (%)")
        ax.set_title(
            f"{PERIOD_LABELS[period]}: "
            f"{fit_percent:.2f}% scale, "
            f"$\\chi^2/\\mathrm{{ndf}}$="
            f"{fit['chi2']:.1f}/{fit['ndf']}"
        )
        ax.grid(alpha=0.25)
        ax.legend(ncol=3)
    # endfor

    axes[-1].set_xlabel("Combined kinematic-bin number")
    axes[-1].set_xticks(bins)
    fig.suptitle(
        "Nominal Method-1/Method-2 relative half-difference "
        f"(adopted global scale: {global_percent:.2f}%)"
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    path = output_dir / "nominal_method_scale_systematic.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_nominal_recommended_dilution(
    output_dir: Path,
    summaries_by_period: dict[str, Any],
) -> str:
    """Plot recommended averages with statistical bars only."""
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    offsets = {"su22": -0.18, "fa22": 0.0, "sp23": 0.18}
    nominal_index = CUT_VARIATIONS.index("nominal")
    global_scale = summaries_by_period[
        "_global_method_scale_by_cut"
    ][nominal_index]
    fig, ax = plt.subplots(figsize=(17, 6.5))

    for period in PERIODS:
        rec = summaries_by_period[period]["recommended"]
        values = rec["central"][:, nominal_index]
        stat = rec["stat_uncertainty"][:, nominal_index]
        x = bins + offsets[period]
        ax.errorbar(
            x,
            values,
            yerr=stat,
            marker="o",
            linestyle="none",
            markersize=4,
            capsize=2,
            label=PERIOD_LABELS[period],
        )
    # endfor

    ax.set_ylim(0.1, 0.6)
    ax.set_xlabel("Combined kinematic-bin number")
    ax.set_ylabel("Recommended dilution factor")
    ax.set_xticks(bins)
    ax.grid(alpha=0.25)
    ax.legend(ncol=3)
    ax.set_title(
        "Nominal recommended dilution factor: Method-1/Method-2 average\n"
        f"error bars: statistical; global correlated dilution-model scale: "
        f"{global_scale['percent']:.2f}%"
    )
    fig.tight_layout()
    path = output_dir / "nominal_recommended_dilution_factor.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def write_plots(
    output_dir: Path,
    central_results: dict[str, dict[str, Any]],
    summaries_by_period: dict[str, Any],
    bootstrap_results: dict[str, dict[str, np.ndarray]],
) -> list[str]:
    """Write point/error-bar plots into one flat plot directory."""
    paths: list[str] = []
    for period in PERIODS:
        paths.append(
            plot_three_method_comparison(output_dir, period, summaries_by_period[period])
        )
        paths.append(
            plot_nominal_method_comparison(output_dir, period, summaries_by_period[period])
        )
        paths.append(
            plot_packing_fraction_summary(output_dir, period, summaries_by_period[period])
        )
    for method_key in ("method1", "method2", "method3"):
        paths.append(
            plot_three_period_comparison(output_dir, method_key, summaries_by_period)
        )
        paths.append(
            plot_nominal_period_comparison(output_dir, method_key, summaries_by_period)
        )
    paths.append(
        plot_nominal_packing_fraction_period_comparison(
            output_dir, summaries_by_period
        )
    )
    paths.append(
        plot_nominal_method_differences(
            output_dir, summaries_by_period, bootstrap_results
        )
    )
    paths.append(plot_nominal_recommended_dilution(output_dir, summaries_by_period))
    return paths




def write_nominal_isr_comparison_products(
    output_dir: Path,
    nominal_central: dict[str, dict[str, np.ndarray]],
    nominal_summaries: dict[str, Any],
    isr_central: dict[str, dict[str, np.ndarray]],
    isr_summaries: dict[str, Any],
) -> dict[str, str]:
    """Write compact nominal-versus-ISR numerical and graphical diagnostics."""
    ensure_directory(output_dir)
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    ensure_directory(tables_dir)
    ensure_directory(plots_dir)

    rows: list[dict[str, Any]] = []
    methods = ("method1", "method2", "method3", "recommended")
    for period in PERIODS:
        for method in methods:
            n_summary = nominal_summaries[period][method]
            i_summary = isr_summaries[period][method]
            for cut_index, cut in enumerate(CUT_VARIATIONS):
                for b in range(NUMBER_OF_BINS):
                    nval = float(n_summary["central"][b, cut_index])
                    ival = float(i_summary["central"][b, cut_index])
                    rows.append(
                        {
                            "period": period,
                            "period_label": PERIOD_LABELS[period],
                            "method": method,
                            "cut": cut,
                            "bin_number": b + 1,
                            "x_index": b // len(MINUS_TPRIME_BINS_GEV2),
                            "t_index": b % len(MINUS_TPRIME_BINS_GEV2),
                            "nominal_value": nval,
                            "isr_value": ival,
                            "isr_minus_nominal": ival - nval,
                            "isr_over_nominal": (
                                ival / nval
                                if math.isfinite(nval) and nval != 0.0
                                else math.nan
                            ),
                            "nominal_stat_uncertainty": float(
                                n_summary["stat_uncertainty"][b, cut_index]
                            ),
                            "isr_stat_uncertainty": float(
                                i_summary["stat_uncertainty"][b, cut_index]
                            ),
                        }
                    )
    frame = pd.DataFrame(rows)
    csv_path = tables_dir / "nominal_vs_isr_dilution_factor_comparison.csv"
    json_path = tables_dir / "nominal_vs_isr_dilution_factor_comparison.json"
    frame.to_csv(csv_path, index=False)
    write_json(
        json_path,
        {
            "definition": (
                "ISR-minus-nominal dilution-factor diagnostics. The two samples "
                "are processed independently with their own channel-selection "
                "cut tables and identical dilution-factor machinery. These "
                "differences are diagnostic only and are not assigned as a "
                "systematic uncertainty here."
            ),
            "rows": frame.to_dict(orient="records"),
        },
    )

    nominal_index = CUT_VARIATIONS.index("nominal")
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    fig, axes = plt.subplots(3, 2, figsize=(18, 12), sharex=True)
    for row_index, period in enumerate(PERIODS):
        n = nominal_summaries[period]["recommended"]
        i = isr_summaries[period]["recommended"]
        nval = n["central"][:, nominal_index]
        ival = i["central"][:, nominal_index]
        nerr = n["stat_uncertainty"][:, nominal_index]
        ierr = i["stat_uncertainty"][:, nominal_index]

        axes[row_index, 0].errorbar(
            bins - 0.10, nval, yerr=nerr, marker="o", linestyle="none",
            capsize=2, label="Nominal sample",
        )
        axes[row_index, 0].errorbar(
            bins + 0.10, ival, yerr=ierr, marker="s", linestyle="none",
            capsize=2, label="ISR sample",
        )
        axes[row_index, 0].set_ylim(0.1, 0.6)
        axes[row_index, 0].set_ylabel("Dilution factor")
        axes[row_index, 0].set_title(
            f"{PERIOD_LABELS[period]} recommended dilution factor"
        )
        axes[row_index, 0].grid(alpha=0.25)
        axes[row_index, 0].legend(loc="best")

        difference = ival - nval
        axes[row_index, 1].axhline(0.0, linewidth=1.0, linestyle="--")
        axes[row_index, 1].plot(
            bins, difference, marker="o", linestyle="none"
        )
        finite = difference[np.isfinite(difference)]
        if finite.size:
            span = max(0.01, 1.25 * float(np.max(np.abs(finite))))
            axes[row_index, 1].set_ylim(-span, span)
        axes[row_index, 1].set_ylabel("$f_{ISR}-f_{nominal}$")
        axes[row_index, 1].set_title(
            f"{PERIOD_LABELS[period]} ISR-minus-nominal"
        )
        axes[row_index, 1].grid(alpha=0.25)

    axes[-1, 0].set_xlabel("Combined kinematic-bin number")
    axes[-1, 1].set_xlabel("Combined kinematic-bin number")
    for ax in axes.flat:
        ax.set_xticks(bins)
    fig.suptitle(
        "Nominal-versus-ISR dilution-factor summary\n"
        "recommended result = average of carbon-template and direct auxiliary-target methods"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    plot_path = plots_dir / "nominal_vs_isr_recommended_summary.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    return {
        "comparison_csv": str(csv_path),
        "comparison_json": str(json_path),
        "summary_plot": str(plot_path),
    }


# =============================================================================
# Command-line interface and main program
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Determine RGC enpi+ dilution factors with three methods."
    )
    parser.add_argument(
        "--tree",
        default=DEFAULT_TREE_NAME,
        help=f"ROOT tree name (default: {DEFAULT_TREE_NAME}).",
    )
    parser.add_argument(
        "--exclusivity-json",
        type=Path,
        default=DEFAULT_CUT_JSON,
        help=(
            "Final carbon-assisted Mx2-cut JSON.  Default is the sibling "
            "channel_selection output path."
        ),
    )
    parser.add_argument(
        "--input",
        action="append",
        default=[],
        metavar="PERIOD:TARGET=PATH",
        help=(
            "Override one ROOT input. May be repeated. Example: "
            "--input su22:NH3=/path/file.root"
        ),
    )
    parser.add_argument(
        "--charge-fractions-json",
        type=Path,
        default=None,
        help="Optional JSON overriding the built-in relative charge fractions.",
    )
    parser.add_argument(
        "--control-min",
        type=float,
        default=DEFAULT_CONTROL_MIN_GEV2,
        help="Method-1 carbon control-window minimum (GeV^2).",
    )
    parser.add_argument(
        "--control-max",
        type=float,
        default=DEFAULT_CONTROL_MAX_GEV2,
        help="Method-1 carbon control-window maximum (GeV^2).",
    )
    parser.add_argument(
        "--replicas",
        type=int,
        default=DEFAULT_REPLICAS,
        help=f"Poisson replicas per period (default: {DEFAULT_REPLICAS}).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Base random seed (default: {DEFAULT_SEED}).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=MAXIMUM_WORKERS,
        help="Worker processes; hard-capped at 7 (default: 7).",
    )
    parser.add_argument(
        "--run-info-csv",
        type=Path,
        default=DEFAULT_RUN_INFO_CSV,
        help=(
            "Sectioned CLAS12 run-information CSV containing accumulated "
            "charge by run (default: clas12_run_info.csv beside the script)."
        ),
    )
    parser.add_argument(
        "--epoch-definitions",
        type=Path,
        default=DEFAULT_EPOCH_DEFINITIONS,
        help=(
            "XLSX run/target/polarization epoch table located beside the "
            "script by default."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Stable output directory (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--isr-input",
        action="append",
        default=[],
        metavar="PERIOD:TARGET=PATH",
        help=(
            "Override one ISR momentum-corrected ROOT input. May be repeated. "
            "Defaults use *_ISR_externalISR_mom_corrections.root for all five targets."
        ),
    )
    parser.add_argument(
        "--isr-exclusivity-json",
        type=Path,
        default=None,
        help=(
            "ISR-specific final channel-selection cut JSON. By default this is "
            "resolved from the sibling channel-selection ISR manifest."
        ),
    )
    parser.add_argument(
        "--channel-selection-manifest",
        type=Path,
        default=DEFAULT_CHANNEL_SELECTION_MANIFEST,
        help=(
            "Channel-selection manifest used to locate the retained ISR cut table."
        ),
    )
    parser.add_argument(
        "--disable-isr",
        action="store_true",
        help="Disable the completely separate ISR dilution-factor diagnostic.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Write numerical products without rendering plots.",
    )
    return parser



def publish_canonical_production_products(
    *,
    master_json_path: Path,
    compact_json_path: Path,
    flat_csv_path: Path,
    count_path: Path,
    configuration_path: Path,
    production_tables_dir: Path,
    production_metadata_dir: Path,
) -> dict[str, str]:
    """Publish stable, version-independent downstream products."""
    ensure_directory(production_tables_dir)
    ensure_directory(production_metadata_dir)

    destinations = {
        "master_json": production_tables_dir / "dilution_factors.json",
        "downstream_json": (
            production_tables_dir / "dilution_factors_production.json"
        ),
        "all_methods_csv": (
            production_tables_dir / "dilution_factors_all_methods.csv"
        ),
        "target_counts_csv": (
            production_tables_dir / "target_counts_all_selections.csv"
        ),
        "configuration_json": (
            production_metadata_dir / "configuration.json"
        ),
    }
    sources = {
        "master_json": master_json_path,
        "downstream_json": compact_json_path,
        "all_methods_csv": flat_csv_path,
        "target_counts_csv": count_path,
        "configuration_json": configuration_path,
    }

    for key, source in sources.items():
        source_resolved = source.resolve()
        destination = destinations[key]
        destination_resolved = destination.resolve()
        if source_resolved != destination_resolved:
            temporary_path = destination.with_suffix(
                destination.suffix + ".tmp"
            )
            shutil.copy2(source_resolved, temporary_path)
            temporary_path.replace(destination)
        # endif
    # endfor

    return {key: str(path) for key, path in destinations.items()}


def main() -> int:
    args = build_parser().parse_args()

    if not (math.isfinite(args.control_min) and math.isfinite(args.control_max)):
        raise RuntimeError("Control-window endpoints must be finite.")
    # endif
    if args.control_min >= args.control_max:
        raise RuntimeError("Control-window minimum must be below the maximum.")
    # endif
    if args.replicas < 2:
        raise RuntimeError("At least two replicas are required.")
    # endif

    workers = max(1, min(int(args.workers), MAXIMUM_WORKERS, os.cpu_count() or 1))
    output_dir = args.output_dir.resolve()
    analysis_dir = output_dir / "analysis"
    tables_dir = analysis_dir / "tables"
    covariance_dir = analysis_dir / "covariance"
    plots_dir = analysis_dir / "plots"
    epoch_dir = analysis_dir / "epochs"
    epoch_plots_dir = epoch_dir / "plots"

    isr_dir = output_dir / "isr"
    isr_tables_dir = isr_dir / "tables"
    isr_covariance_dir = isr_dir / "covariance"
    isr_plots_dir = isr_dir / "plots"
    diagnostics_dir = output_dir / "diagnostics"

    production_dir = output_dir / "production"
    production_tables_dir = production_dir / "tables"
    production_metadata_dir = production_dir / "metadata"

    for directory in (
        output_dir,
        analysis_dir,
        tables_dir,
        covariance_dir,
        plots_dir,
        epoch_dir,
        epoch_plots_dir,
        isr_dir,
        isr_tables_dir,
        isr_covariance_dir,
        isr_plots_dir,
        diagnostics_dir,
        production_dir,
        production_tables_dir,
        production_metadata_dir,
    ):
        ensure_directory(directory)
    # endfor

    input_paths = {
        period: dict(targets)
        for period, targets in DEFAULT_INPUTS.items()
    }
    for override_text in args.input:
        period, target, path = parse_input_override(override_text)
        input_paths[period][target] = path
    # endfor

    charge_fractions = load_charge_fractions(args.charge_fractions_json)
    cut_json = args.exclusivity_json.resolve()
    cuts = load_exclusivity_cuts(cut_json)

    print("=" * 78)
    print("RGC exclusive enpi+ dilution-factor determination")
    print("=" * 78)
    print(f"Working directory:       {Path.cwd()}")
    print(f"Exclusivity JSON:        {cut_json}")
    epoch_path = args.epoch_definitions.resolve()
    run_info_path = args.run_info_csv.resolve()
    print(f"Output root:             {output_dir}")
    print(f"Persistent analysis:     {analysis_dir}")
    print(f"Canonical production:    {production_dir}")
    print(f"Epoch definitions:       {epoch_path}")
    print(f"Run information:         {run_info_path}")
    print(f"Tree:                    {args.tree}")
    print(f"Workers:                 {workers} (hard maximum {MAXIMUM_WORKERS})")
    print(f"Poisson replicas/period: {args.replicas:,}")
    print(
        f"Method-1 control region: [{args.control_min:.6g}, "
        f"{args.control_max:.6g}) GeV^2"
    )
    print()

    loaded = load_all_datasets(input_paths, args.tree, workers)
    pattern_counts, observed_counts, count_rows = build_count_arrays(
        loaded,
        cuts,
        args.control_min,
        args.control_max,
    )
    count_frame = pd.DataFrame(count_rows)
    count_path = tables_dir / "target_counts_all_selections.csv"
    count_frame.to_csv(count_path, index=False)

    central_results = {
        period: observed_estimators_for_period(
            observed_counts[p_index],
            period,
            charge_fractions,
        )
        for p_index, period in enumerate(PERIODS)
    }

    bootstrap_results = run_bootstrap(
        pattern_counts,
        charge_fractions,
        args.replicas,
        args.seed,
        workers,
    )

    flat_frame, master_periods, compact_periods, summaries_by_period = build_output_tables(
        observed_counts,
        central_results,
        bootstrap_results,
        charge_fractions,
        cuts,
    )
    flat_csv_path = tables_dir / "dilution_factors_all_methods.csv"
    flat_frame.to_csv(flat_csv_path, index=False)

    epoch_definitions = read_epoch_spreadsheet(epoch_path)
    run_charges = parse_run_info_csv(run_info_path)
    epoch_frame = calculate_epoch_diagnostics(
        loaded=loaded,
        observed_counts=observed_counts,
        cuts=cuts,
        definitions=epoch_definitions,
        run_charges=run_charges,
        control_min_gev2=args.control_min,
        control_max_gev2=args.control_max,
        replicas=args.replicas,
        seed=args.seed,
    )
    epoch_csv_path = epoch_dir / "nh3_epoch_diagnostics.csv"
    epoch_json_path = epoch_dir / "nh3_epoch_diagnostics.json"
    epoch_frame.to_csv(epoch_csv_path, index=False)
    write_json(
        epoch_json_path,
        {
            "definition": (
                "NH3 epochs are defined by inclusive run ranges from the "
                "spreadsheet. Actual NH3 runs are taken from the ROOT runnum "
                "branch, and their accumulated charge is summed from "
                "clas12_run_info.csv. All C, CH2, He, and ET counts and charges "
                "remain full-period values for every epoch. Charge fractions "
                "are recomputed separately for each epoch."
            ),
            "rows": epoch_frame.to_dict(orient="records"),
        },
    )

    covariance_manifest = write_covariance_products(
        covariance_dir,
        bootstrap_results,
    )

    plot_paths: list[str] = []
    if not args.skip_plots:
        plot_paths = write_plots(
            plots_dir, central_results, summaries_by_period, bootstrap_results
        )
        plot_paths.extend(
            plot_epoch_diagnostics(epoch_frame, epoch_plots_dir)
        )
    # endif

    # -------------------------------------------------------------------------
    # Completely separate ISR diagnostic workflow
    # -------------------------------------------------------------------------
    isr_payload: dict[str, Any] | None = None
    comparison_products: dict[str, str] | None = None
    if not args.disable_isr:
        isr_input_paths = {
            period: dict(targets)
            for period, targets in DEFAULT_ISR_INPUTS.items()
        }
        for override_text in args.isr_input:
            period, target, path = parse_input_override(override_text)
            isr_input_paths[period][target] = path

        isr_cut_json = resolve_isr_cut_json(
            args.isr_exclusivity_json,
            args.channel_selection_manifest,
        )
        isr_cuts = load_exclusivity_cuts(isr_cut_json)

        print()
        print("-" * 78)
        print("Starting completely separate ISR dilution-factor diagnostic")
        print("-" * 78)
        print(f"ISR exclusivity JSON:    {isr_cut_json}")
        print(f"ISR output directory:    {isr_dir}")

        isr_loaded = load_all_datasets(isr_input_paths, args.tree, workers)
        (
            isr_pattern_counts,
            isr_observed_counts,
            isr_count_rows,
        ) = build_count_arrays(
            isr_loaded,
            isr_cuts,
            args.control_min,
            args.control_max,
        )
        isr_count_path = isr_tables_dir / "target_counts_all_selections.csv"
        pd.DataFrame(isr_count_rows).to_csv(isr_count_path, index=False)

        isr_central_results = {
            period: observed_estimators_for_period(
                isr_observed_counts[p_index],
                period,
                charge_fractions,
            )
            for p_index, period in enumerate(PERIODS)
        }
        isr_bootstrap_results = run_bootstrap(
            isr_pattern_counts,
            charge_fractions,
            args.replicas,
            args.seed + 1000003,
            workers,
        )
        (
            isr_flat_frame,
            isr_master_periods,
            isr_compact_periods,
            isr_summaries_by_period,
        ) = build_output_tables(
            isr_observed_counts,
            isr_central_results,
            isr_bootstrap_results,
            charge_fractions,
            isr_cuts,
        )
        isr_flat_csv_path = (
            isr_tables_dir / "dilution_factors_all_methods.csv"
        )
        isr_flat_frame.to_csv(isr_flat_csv_path, index=False)

        isr_covariance_manifest = write_covariance_products(
            isr_covariance_dir,
            isr_bootstrap_results,
        )
        isr_plot_paths: list[str] = []
        if not args.skip_plots:
            isr_plot_paths = write_plots(
                isr_plots_dir,
                isr_central_results,
                isr_summaries_by_period,
                isr_bootstrap_results,
            )

        isr_master_json_path = isr_dir / "dilution_factors.json"
        isr_compact_json_path = isr_dir / "dilution_factors_production.json"
        isr_configuration_path = isr_dir / "configuration.json"
        isr_provenance = {
            "sample_variant": "ISR",
            "diagnostic_only": True,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "tree_name": args.tree,
            "workers_used": workers,
            "bootstrap_replicas_per_period": args.replicas,
            "bootstrap_seed": args.seed + 1000003,
            "exclusivity_json": str(isr_cut_json),
            "exclusivity_json_sha256": sha256_file(isr_cut_json),
            "control_region_gev2": [args.control_min, args.control_max],
            "charge_fractions": charge_fractions,
            "root_inputs": {
                period: {
                    target: {
                        "path": isr_loaded[(period, target)].file_path,
                        "branches": {
                            "xB": isr_loaded[(period, target)].x_branch,
                            "tprime": isr_loaded[(period, target)].tprime_branch,
                            "Mx2": isr_loaded[(period, target)].mx2_branch,
                            "runnum": isr_loaded[(period, target)].run_branch,
                        },
                    }
                    for target in TARGETS
                }
                for period in PERIODS
            },
        }
        write_json(isr_configuration_path, isr_provenance)
        write_json(
            isr_master_json_path,
            {
                "schema_version": 1,
                "analysis": "RGC exclusive enpi+ ISR dilution-factor diagnostic",
                "diagnostic_only": True,
                "binning": {
                    "xB": [list(interval) for interval in XB_BINS],
                    "minus_tprime_gev2": [
                        list(interval) for interval in MINUS_TPRIME_BINS_GEV2
                    ],
                    "combined_bin_order": (
                        "x_index major, t_index minor, one-based bin_number"
                    ),
                },
                "cut_variations": list(CUT_VARIATIONS),
                "provenance": isr_provenance,
                "global_dilution_model_scale_by_cut": isr_summaries_by_period[
                    "_global_method_scale_by_cut"
                ],
                "periods": isr_master_periods,
                "covariance_manifest": isr_covariance_manifest,
                "plot_paths": isr_plot_paths,
            },
        )
        write_json(
            isr_compact_json_path,
            {
                "schema_version": 1,
                "analysis": (
                    "RGC exclusive enpi+ ISR dilution factors; diagnostic only"
                ),
                "diagnostic_only": True,
                "recommended_nominal_method": "average_of_method1_and_method2",
                "periods": isr_compact_periods,
                "source_master_json": str(isr_master_json_path),
                "source_exclusivity_json": str(isr_cut_json),
            },
        )

        comparison_products = write_nominal_isr_comparison_products(
            diagnostics_dir,
            central_results,
            summaries_by_period,
            isr_central_results,
            isr_summaries_by_period,
        )
        isr_payload = {
            "root": str(isr_dir),
            "master_json": str(isr_master_json_path),
            "downstream_json": str(isr_compact_json_path),
            "all_methods_csv": str(isr_flat_csv_path),
            "target_counts_csv": str(isr_count_path),
            "covariance": str(isr_covariance_dir),
            "plots": str(isr_plots_dir),
            "exclusivity_json": str(isr_cut_json),
            "comparison_products": comparison_products,
        }

    provenance = {
        "script": Path(__file__).name,
        "schema_version": 12,
        "output_filename_policy": "stable_version_independent",
        "analysis_products_directory": str(analysis_dir),
        "production_products_directory": str(production_dir),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "working_directory": str(Path.cwd()),
        "tree_name": args.tree,
        "maximum_worker_cap": MAXIMUM_WORKERS,
        "workers_used": workers,
        "bootstrap_replicas_per_period": args.replicas,
        "bootstrap_seed": args.seed,
        "exclusivity_json": str(cut_json),
        "exclusivity_json_sha256": sha256_file(cut_json),
        "control_region_gev2": [args.control_min, args.control_max],
        "epoch_definitions_xlsx": str(epoch_path),
        "epoch_definitions_xlsx_sha256": sha256_file(epoch_path),
        "run_info_csv": str(run_info_path),
        "run_info_csv_sha256": sha256_file(run_info_path),
        "run_info_charge_definition": (
            "Accumulated charge per run is CSV column 3 plus CSV column 4. Only unique run numbers found in the corresponding ROOT tree are summed; unused zero-charge CSV rows are ignored."
        ),
        "epoch_charge_treatment": (
            "NH3 charge is the sum of clas12_run_info.csv charges for the "
            "actual ROOT runnum values inside each epoch range. C, CH2, He, "
            "and ET retain their full-period charges. Fractions are "
            "renormalized for each epoch."
        ),
        "isr_diagnostic_enabled": not args.disable_isr,
        "isr_diagnostic": isr_payload,
        "root_inputs": {
            period: {
                target: {
                    "path": loaded[(period, target)].file_path,
                    "branches": {
                        "xB": loaded[(period, target)].x_branch,
                        "tprime": loaded[(period, target)].tprime_branch,
                        "Mx2": loaded[(period, target)].mx2_branch,
                        "runnum": loaded[(period, target)].run_branch,
                    },
                }
                for target in TARGETS
            }
            for period in PERIODS
        },
        "charge_fractions": charge_fractions,
        "global_dilution_model_scale_by_cut": summaries_by_period[
            "_global_method_scale_by_cut"
        ],
        "equation_coefficients": EQ10_EQ11_EQ14_COEFFICIENTS,
        "method1_definition": (
            "Exact channel-selection raw-count normalization pooled over all "
            "24 bins in the control window; no additional charge normalization."
        ),
        "method2_definition": "Exact thermal-contraction-corrected formula transcribed from calculate_dilution_factors.cpp.",
        "method3_definition": (
            "The period-wide packing fraction is calculated from the thermally "
            "corrected C++ relation. Method 3 reconstructs the binwise ammonia "
            "rate from that packing fraction and inserts it into the same "
            "thermally corrected dilution model; this is the algebraic Eq. (14) "
            "counterpart for the corrected coefficients."
        ),
        "statistical_model": (
            "Independent Poisson replicas of 16 disjoint membership-pattern "
            "counts for control/tight/nominal/loose selections."
        ),
    }

    master_payload = {
        "schema_version": 12,
        "analysis": "RGC exclusive enpi+ dilution-factor determination",
        "status": (
            "All three methods are active. Method 2 uses the exact thermally corrected "
            "direct C++ expression; Method 3 uses the period-wide packing "
            "fraction in its algebraically equivalent nonlinear formulation."
        ),
        "binning": {
            "xB": [list(interval) for interval in XB_BINS],
            "minus_tprime_gev2": [list(interval) for interval in MINUS_TPRIME_BINS_GEV2],
            "combined_bin_order": (
                "x_index major, t_index minor, one-based bin_number"
            ),
        },
        "cut_variations": list(CUT_VARIATIONS),
        "provenance": provenance,
        "global_dilution_model_scale_by_cut": summaries_by_period[
            "_global_method_scale_by_cut"
        ],
        "epoch_diagnostics_csv": str(epoch_csv_path),
        "epoch_diagnostics_json": str(epoch_json_path),
        "periods": master_periods,
        "covariance_manifest": covariance_manifest,
        "plot_paths": plot_paths,
        "isr_diagnostic": isr_payload,
        "nominal_isr_comparison_products": comparison_products,
    }
    master_json_path = analysis_dir / "dilution_factors.json"
    write_json(master_json_path, master_payload)

    compact_payload = {
        "schema_version": 12,
        "analysis": "RGC exclusive enpi+ dilution factors for downstream asymmetries",
        "recommended_nominal_method": "average_of_method1_and_method2",
        "note": (
            "Loose, nominal, and tight values are all carried downstream. "
            "The recommended dilution factor is the average of Methods 1 and 2. "
            "Its statistical uncertainty is propagated with common bootstrap "
            "replicas. The Method-1/Method-2 disagreement is represented by a "
            "single global multiplicative scale systematic obtained by averaging "
            "the three period-specific constant fits to |f1-f2|/(f1+f2); it is "
            "fully correlated across all periods and bins. No exclusivity-cut "
            "systematic is assigned "
            "at this stage."
        ),
        "binning": master_payload["binning"],
        "global_dilution_model_scale_by_cut": summaries_by_period[
            "_global_method_scale_by_cut"
        ],
        "epoch_diagnostics_json": str(epoch_json_path),
        "periods": compact_periods,
        "source_master_json": str(master_json_path),
        "source_exclusivity_json": str(cut_json),
        "source_exclusivity_json_sha256": provenance["exclusivity_json_sha256"],
    }
    compact_json_path = analysis_dir / "dilution_factors_production.json"
    write_json(compact_json_path, compact_payload)

    configuration_path = analysis_dir / "configuration.json"
    write_json(configuration_path, provenance)

    production_products = publish_canonical_production_products(
        master_json_path=master_json_path,
        compact_json_path=compact_json_path,
        flat_csv_path=flat_csv_path,
        count_path=count_path,
        configuration_path=configuration_path,
        production_tables_dir=production_tables_dir,
        production_metadata_dir=production_metadata_dir,
    )

    manifest_path = output_dir / "determine_dilution_factor_manifest.json"
    manifest_payload = {
        "schema_version": 1,
        "analysis": "RGC exclusive enpi+ dilution-factor determination",
        "production_policy": (
            "The complete analysis is retained under analysis/. Downstream "
            "consumers use the stable files under production/."
        ),
        "source_channel_selection_json": str(cut_json),
        "source_channel_selection_json_sha256": provenance[
            "exclusivity_json_sha256"
        ],
        "retained_analysis_products": {
            "root": str(analysis_dir),
            "tables": str(tables_dir),
            "covariance": str(covariance_dir),
            "plots": str(plots_dir),
            "epochs": str(epoch_dir),
            "isr": str(isr_dir) if not args.disable_isr else None,
            "nominal_isr_comparisons": (
                str(diagnostics_dir) if not args.disable_isr else None
            ),
        },
        "production_products": production_products,
        "recommended_nominal_method": "average_of_method1_and_method2",
        "filenames_are_version_independent": True,
    }
    write_json(manifest_path, manifest_payload)

    print()
    print(f"Epoch diagnostics:       {epoch_csv_path}")
    print("Dilution-factor calculation complete.")
    print(f"  Complete analysis:     {analysis_dir}")
    print(f"  Master JSON:           {master_json_path}")
    print(f"  Downstream JSON:       {compact_json_path}")
    print(f"  Flat CSV:              {flat_csv_path}")
    print(f"  Count CSV:             {count_path}")
    print(f"  Covariance directory:  {covariance_dir}")
    if not args.skip_plots:
        print(f"  Plot directory:        {plots_dir}")
    # endif
    if not args.disable_isr:
        print(f"  ISR diagnostic:        {isr_dir}")
        if comparison_products is not None:
            print(
                "  Nominal/ISR summary:   "
                f"{comparison_products['summary_plot']}"
            )
    print(f"  Production directory:  {production_dir}")
    print(
        "  Canonical downstream:  "
        f"{production_products['downstream_json']}"
    )
    print(f"  Manifest:              {manifest_path}")
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
