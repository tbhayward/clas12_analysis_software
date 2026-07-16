#!/usr/bin/env python3
"""
check_rcdb_currents.py

Validate the nominal beam-current assignments used in the RGA DVCS analysis
against the CLAS12 Run Conditions Database (RCDB).

For every run in the analysis current tables, the script reads:

    beam_current_request
    beam_current
    torus_scale

It then compares the analysis assignment against the RCDB requested current.

Intentional manual assignments:

    Fa18 Inb runs 5340--5344 -> 40 nA
    Fa18 Out run 5615       -> 50 nA

Run 3218 is intentionally omitted because it is removed from the final
current-dependence study.

Outputs
-------
rcdb_current_check/
    rcdb_run_conditions.csv
        Full run-by-run RCDB output.

    rcdb_current_discrepancies.csv
        Only missing values, mismatches, manual overrides and torus problems.

    rcdb_current_summary.txt
        Human-readable summary.

    rcdb_current_validation_table.tex
        Compact LaTeX table containing only entries requiring attention.

    rcdb_run_conditions_tables.tex
        Full run-by-run LaTeX tables, divided by run period and split into
        manageable ordinary table environments. No longtable package is used.

Status definitions
------------------
MATCH
    The analysis current equals beam_current_request.

MANUAL_OVERRIDE
    The analysis current differs from beam_current_request, but the run is
    one of the explicitly documented manual reassignments.

MISMATCH
    The analysis current differs from beam_current_request and no documented
    manual override applies.

MISSING_REQUEST
    beam_current_request is absent or cannot be interpreted numerically.

MISSING_BEAM_CURRENT
    beam_current is absent or cannot be interpreted numerically. This does
    not by itself invalidate the nominal-current assignment.

MISSING_TORUS
    torus_scale is absent or cannot be interpreted numerically.

TORUS_MISMATCH
    The RCDB torus polarity is inconsistent with the expected run-period
    polarity.

RCDB_ERROR
    RCDB could not return one or more conditions for the run.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

try:
    import rcdb
except ImportError as exc:
    print(
        "ERROR: The Python rcdb package could not be imported.\n"
        "On the JLab farm, first source the CLAS12 environment:\n\n"
        "  tcsh/csh:\n"
        "    source /group/clas12/environment.csh\n\n"
        "  bash:\n"
        "    source /group/clas12/environment.sh\n\n"
        "Then rerun this script.",
        file=sys.stderr,
    )
    raise SystemExit(1) from exc
#endif


DEFAULT_CONNECTION = "mysql://rcdb@clasdb/rcdb"
DEFAULT_OUTPUT_DIRECTORY = Path("rcdb_current_check")

# Current comparisons are effectively comparisons of integer current settings.
# This tolerance allows harmless floating-point or formatting differences such
# as 49.999999 versus 50.
CURRENT_MATCH_TOLERANCE_NA = 0.25

# RCDB torus_scale may contain a value very close to +/-1 rather than exactly
# +/-1. The sign is what matters for this check.
TORUS_ZERO_TOLERANCE = 1.0e-6

# Delay between runs. The local MySQL query is normally fast, but a small delay
# avoids hammering the server unnecessarily.
DEFAULT_QUERY_DELAY_SECONDS = 0.01

# Number of rows placed in each ordinary LaTeX table. This avoids longtable,
# which conflicts with the user's current REVTeX/caption setup.
DEFAULT_LATEX_ROWS_PER_TABLE = 35


@dataclass(frozen=True)
class RunAssignment:
    """Analysis beam-current assignment for one run."""

    period: str
    run: int
    assigned_current_na: float
    expected_torus_sign: int


@dataclass
class RCDBResult:
    """RCDB values and validation status for one run."""

    period: str
    run: int
    assigned_current_na: float
    expected_torus_sign: int
    beam_current_request_na: Optional[float]
    beam_current_na: Optional[float]
    torus_scale: Optional[float]
    status: str
    notes: str


def parse_run_list(text: str) -> list[int]:
    """
    Convert a comma- or whitespace-separated run list into a list of integers.

    The parser also accepts inclusive ranges written with a hyphen, for example:

        "5340-5344"

    Duplicate run numbers are removed while preserving first appearance.
    """

    normalized = text.replace("\n", " ").replace(",", " ")
    tokens = normalized.split()

    runs: list[int] = []
    seen: set[int] = set()

    for token in tokens:
        cleaned = token.strip().rstrip(".")
        if not cleaned:
            continue
        #endif

        range_match = re.fullmatch(r"(\d+)-(\d+)", cleaned)
        if range_match:
            start = int(range_match.group(1))
            stop = int(range_match.group(2))

            if stop < start:
                raise ValueError(
                    f"Invalid descending run range '{cleaned}'."
                )
            #endif

            for run in range(start, stop + 1):
                if run not in seen:
                    runs.append(run)
                    seen.add(run)
                #endif
            #endfor

            continue
        #endif

        if not re.fullmatch(r"\d+", cleaned):
            raise ValueError(f"Cannot interpret run-list token '{token}'.")
        #endif

        run = int(cleaned)
        if run not in seen:
            runs.append(run)
            seen.add(run)
        #endif
    #endfor

    return runs


# ---------------------------------------------------------------------------
# Final nominal current groups
# ---------------------------------------------------------------------------
#
# Run 3218 is intentionally absent from Sp18 Out because it is removed from
# the current-dependence study.
#
# Runs 5247, 5345, 5418, 5419, 5443, 6616 and 6618 are also absent because
# they are removed from the final analysis/current study.
#
# Runs 5158, 5163, 5181, 5519, 5528 and 5627 remain present because they are
# retained in the updated analysis.
# ---------------------------------------------------------------------------

RUN_GROUPS: dict[str, dict[float, list[int]]] = {
    "Sp18 Inb": {
        35.0: parse_run_list(
            """
            3306, 3307, 3315, 3333, 3353, 3359, 3361, 3363, 3378,
            3379, 3384, 3389, 3390, 3403, 3405, 3406, 3407, 3409,
            3411, 3421, 3422
            """
        ),
        50.0: parse_run_list(
            """
            3429, 3431, 3432, 3433, 3434, 3435, 3436, 3441, 3442,
            3459, 3460, 3461, 3462, 3463, 3464, 3465, 3466, 3467,
            3469, 3480, 3482, 3484, 3485, 3488, 3492, 3493, 3501,
            3506, 3507, 3512, 3513, 3517, 3518, 3519, 3520, 3521,
            3522, 3542, 3699, 3700, 3702, 3705, 3708, 3711, 3719,
            3720, 3722, 3738, 3739, 3741, 3748, 3750, 3752, 3771,
            3789, 3790, 3791, 3795, 3796, 3797, 3798, 3799, 3803,
            3804, 3806, 3816, 4003, 4013, 4014, 4015, 4016, 4017,
            4021, 4022, 4025, 4026, 4028, 4030, 4032, 4033, 4037,
            4038, 4039, 4041, 4044, 4045, 4050, 4053, 4054, 4055,
            4058, 4060, 4061, 4067, 4068, 4069, 4070, 4071, 4073,
            4075, 4078, 4080, 4081, 4082, 4083, 4085, 4089, 4090,
            4091, 4092, 4093, 4094, 4095, 4096, 4097, 4099, 4100,
            4103, 4104, 4110, 4112, 4113, 4114, 4115, 4139, 4143,
            4144, 4147, 4148, 4151, 4152, 4154, 4155, 4156, 4157,
            4158, 4161, 4164, 4165, 4166, 4167, 4169, 4174, 4180,
            4181, 4182, 4184, 4189, 4190, 4192, 4193, 4194, 4201,
            4202, 4203, 4204, 4206, 4208, 4211, 4213, 4218, 4219,
            4220, 4221, 4223, 4224, 4226, 4229, 4243, 4244, 4245,
            4247, 4248, 4250, 4253, 4254, 4256, 4262, 4263, 4264,
            4309, 4311, 4312, 4313, 4314, 4315, 4316, 4320, 4321,
            4322, 4323, 4324
            """
        ),
    },
    "Sp18 Out": {
        30.0: parse_run_list(
            """
            3261, 3262, 3266, 3269, 3270, 3282, 3288
            """
        ),
        45.0: parse_run_list(
            """
            3874, 3875, 3878, 3880, 3881, 3883, 3884, 3885, 3888,
            3889, 3891, 3893, 3898, 3903, 3905, 3907, 3908, 3910,
            3911, 3912, 3913, 3915, 3916, 3917, 3919, 3920, 3921,
            3924, 3926, 3928, 3930, 3932, 3933, 3934, 3936, 3938,
            3939, 3940, 3941, 3943, 3944, 3945, 3946, 3948, 3949,
            3950, 3954, 3959, 3963, 3964, 3969, 3970, 3973, 3975,
            3982, 3985, 3986, 3987
            """
        ),
    },
    "Fa18 Inb": {
        40.0: parse_run_list(
            """
            5335, 5339, 5340, 5341, 5342, 5343, 5344
            """
        ),
        45.0: parse_run_list(
            """
            5032, 5036, 5038, 5039, 5040, 5041, 5043, 5045, 5046,
            5047, 5051, 5052, 5053, 5116, 5117, 5119, 5120, 5124,
            5125, 5126, 5127, 5128, 5129, 5130, 5139, 5153, 5158,
            5159, 5160, 5162, 5163, 5164, 5165, 5166, 5167, 5168,
            5169, 5180, 5181, 5182, 5183, 5190, 5191, 5193, 5195,
            5196, 5197, 5198, 5199, 5200, 5201, 5202, 5203, 5204,
            5205, 5206, 5208, 5211, 5212, 5215, 5216, 5219, 5220,
            5221, 5222, 5223, 5230, 5231, 5232, 5233, 5234, 5235,
            5237, 5238, 5239, 5248, 5249, 5252, 5253, 5257, 5258,
            5259, 5261, 5262, 5303, 5304, 5305, 5306, 5307, 5310,
            5311, 5315, 5317, 5318, 5319, 5320, 5323, 5324, 5333,
            5334, 5336, 5346, 5347, 5349, 5351, 5354, 5355, 5367
            """
        ),
        50.0: parse_run_list(
            """
            5356, 5357, 5358, 5359, 5360, 5361, 5362, 5366
            """
        ),
        55.0: parse_run_list(
            """
            5368, 5369, 5372, 5373, 5374, 5375, 5376, 5377, 5378,
            5379, 5380, 5381, 5382, 5383, 5386, 5390, 5391, 5392,
            5393, 5398, 5400, 5401, 5403, 5404, 5406, 5407
            """
        ),
    },
    "Fa18 Out": {
        20.0: parse_run_list(
            """
            5444
            """
        ),
        40.0: parse_run_list(
            """
            5423, 5424, 5425, 5426, 5428, 5429, 5430, 5432, 5434,
            5435, 5436, 5437, 5438, 5440, 5441, 5442, 5445, 5447,
            5448, 5449, 5450, 5451, 5452, 5453, 5454, 5455, 5460,
            5464, 5465, 5466, 5467, 5468, 5469, 5470, 5471, 5472,
            5473, 5474, 5475, 5476, 5478, 5479, 5480, 5481, 5482,
            5483, 5485, 5486, 5487, 5495, 5496, 5497, 5498, 5499,
            5500, 5504
            """
        ),
        50.0: parse_run_list(
            """
            5505, 5507, 5516, 5517, 5518, 5519, 5520, 5521, 5522,
            5523, 5524, 5525, 5526, 5527, 5528, 5530, 5532, 5533,
            5534, 5535, 5536, 5537, 5538, 5540, 5541, 5543, 5544,
            5545, 5546, 5547, 5548, 5549, 5550, 5551, 5552, 5555,
            5556, 5557, 5558, 5559, 5562, 5567, 5569, 5570, 5571,
            5572, 5573, 5574, 5577, 5578, 5591, 5592, 5594, 5597,
            5598, 5600, 5601, 5602, 5603, 5604, 5606, 5607, 5611,
            5612, 5613, 5614, 5615, 5616, 5617, 5618, 5619, 5621,
            5623, 5624, 5625, 5626, 5627, 5628, 5629, 5630, 5631,
            5632, 5633, 5635, 5637, 5638, 5639, 5641, 5643, 5644,
            5645, 5646, 5647, 5648, 5649, 5650, 5651, 5652, 5654,
            5655, 5656, 5662, 5663, 5664, 5665, 5666
            """
        ),
    },
    "Sp19 Inb": {
        50.0: parse_run_list(
            """
            6619, 6620, 6636, 6637, 6638, 6639, 6640, 6642, 6645,
            6647, 6648, 6650, 6651, 6652, 6654, 6655, 6656, 6657,
            6658, 6660, 6661, 6662, 6663, 6664, 6665, 6666, 6667,
            6668, 6669, 6670, 6672, 6673, 6675, 6676, 6677, 6678,
            6680, 6682, 6683, 6684, 6685, 6687, 6688, 6689, 6691,
            6692, 6693, 6694, 6695, 6696, 6697, 6698, 6699, 6704,
            6705, 6706, 6707, 6708, 6709, 6710, 6711, 6712, 6713,
            6714, 6715, 6716, 6717, 6718, 6719, 6729, 6730, 6731,
            6732, 6733, 6736, 6737, 6738, 6739, 6740, 6741, 6742,
            6743, 6744, 6746, 6747, 6748, 6749, 6750, 6753, 6754,
            6755, 6756, 6757, 6759, 6760, 6762, 6763, 6764, 6765,
            6767, 6768, 6769, 6779, 6780, 6781, 6783
            """
        ),
    },
}


EXPECTED_TORUS_SIGN: dict[str, int] = {
    "Sp18 Inb": -1,
    "Sp18 Out": +1,
    "Fa18 Inb": -1,
    "Fa18 Out": +1,
    "Sp19 Inb": -1,
}


# These are analysis assignments intentionally retained even if the RCDB
# beam_current_request condition reports a different value.
MANUAL_CURRENT_OVERRIDES: dict[int, float] = {
    5340: 40.0,
    5341: 40.0,
    5342: 40.0,
    5343: 40.0,
    5344: 40.0,
    5615: 50.0,
}


def build_assignments() -> list[RunAssignment]:
    """Flatten RUN_GROUPS and validate that each run appears exactly once."""

    assignments: list[RunAssignment] = []
    seen: dict[int, tuple[str, float]] = {}

    for period, current_groups in RUN_GROUPS.items():
        if period not in EXPECTED_TORUS_SIGN:
            raise ValueError(
                f"No expected torus polarity was defined for '{period}'."
            )
        #endif

        expected_torus_sign = EXPECTED_TORUS_SIGN[period]

        for assigned_current_na, runs in current_groups.items():
            for run in runs:
                if run in seen:
                    old_period, old_current = seen[run]
                    raise ValueError(
                        f"Run {run} appears more than once: "
                        f"{old_period} at {old_current:g} nA and "
                        f"{period} at {assigned_current_na:g} nA."
                    )
                #endif

                seen[run] = (period, assigned_current_na)
                assignments.append(
                    RunAssignment(
                        period=period,
                        run=run,
                        assigned_current_na=assigned_current_na,
                        expected_torus_sign=expected_torus_sign,
                    )
                )
            #endfor
        #endfor
    #endfor

    assignments.sort(key=lambda item: item.run)
    return assignments


def condition_value(
    provider: Any,
    run: int,
    condition_name: str,
) -> Any:
    """
    Return one RCDB condition value.

    A missing condition is represented by None. RCDB exceptions are propagated
    so the calling function can record an RCDB_ERROR.
    """

    condition = provider.get_condition(run, condition_name)

    if condition is None:
        return None
    #endif

    return condition.value


def to_float(value: Any) -> Optional[float]:
    """
    Convert an RCDB value to a float.

    This accepts numeric values and strings such as:

        "50"
        "50.0"
        "50 nA"
        "-1.0"

    None, empty strings, NaN and infinite values return None.
    """

    if value is None:
        return None
    #endif

    if isinstance(value, bool):
        return float(value)
    #endif

    if isinstance(value, (int, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    #endif

    text = str(value).strip()
    if not text:
        return None
    #endif

    match = re.search(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        text,
    )

    if match is None:
        return None
    #endif

    try:
        number = float(match.group(0))
    except ValueError:
        return None
    #endtry

    return number if math.isfinite(number) else None


def currents_match(left: float, right: float) -> bool:
    """Return True when two current settings agree within the tolerance."""

    return abs(left - right) <= CURRENT_MATCH_TOLERANCE_NA


def torus_sign(value: Optional[float]) -> Optional[int]:
    """Convert torus_scale to -1, 0 or +1."""

    if value is None:
        return None
    #endif

    if abs(value) <= TORUS_ZERO_TOLERANCE:
        return 0
    #endif

    return +1 if value > 0.0 else -1


def classify_result(
    assignment: RunAssignment,
    beam_current_request_na: Optional[float],
    beam_current_na: Optional[float],
    torus_scale: Optional[float],
) -> tuple[str, str]:
    """Assign the primary validation status and explanatory notes."""

    issues: list[str] = []
    current_status: str

    if beam_current_request_na is None:
        current_status = "MISSING_REQUEST"
        issues.append("beam_current_request is missing or nonnumeric")
    elif currents_match(
        assignment.assigned_current_na,
        beam_current_request_na,
    ):
        current_status = "MATCH"
    elif assignment.run in MANUAL_CURRENT_OVERRIDES:
        manual_value = MANUAL_CURRENT_OVERRIDES[assignment.run]

        if currents_match(assignment.assigned_current_na, manual_value):
            current_status = "MANUAL_OVERRIDE"
            issues.append(
                "intentional analysis reassignment: "
                f"RCDB request {beam_current_request_na:g} nA, "
                f"analysis assignment {assignment.assigned_current_na:g} nA"
            )
        else:
            current_status = "MISMATCH"
            issues.append(
                "run has a manual-override entry, but the table assignment "
                f"{assignment.assigned_current_na:g} nA does not equal the "
                f"documented override {manual_value:g} nA"
            )
        #endif
    else:
        current_status = "MISMATCH"
        issues.append(
            f"RCDB request {beam_current_request_na:g} nA differs from "
            f"analysis assignment {assignment.assigned_current_na:g} nA"
        )
    #endif

    actual_torus_sign = torus_sign(torus_scale)

    if torus_scale is None:
        issues.append("torus_scale is missing or nonnumeric")
    elif actual_torus_sign != assignment.expected_torus_sign:
        issues.append(
            f"torus sign {actual_torus_sign:+d} differs from expected "
            f"{assignment.expected_torus_sign:+d}"
        )
    #endif

    if beam_current_na is None:
        issues.append("beam_current is missing or nonnumeric")
    #endif

    # Keep current mismatches and missing requests as the primary status.
    # Otherwise promote torus or measured-current problems to the status field.
    if current_status in {
        "MATCH",
        "MANUAL_OVERRIDE",
    }:
        if torus_scale is None:
            primary_status = "MISSING_TORUS"
        elif actual_torus_sign != assignment.expected_torus_sign:
            primary_status = "TORUS_MISMATCH"
        elif beam_current_na is None:
            primary_status = "MISSING_BEAM_CURRENT"
        else:
            primary_status = current_status
        #endif
    else:
        primary_status = current_status
    #endif

    notes = "; ".join(issues)
    return primary_status, notes


def query_one_run(
    provider: Any,
    assignment: RunAssignment,
) -> RCDBResult:
    """Read and validate all requested RCDB conditions for one run."""

    try:
        raw_request = condition_value(
            provider,
            assignment.run,
            "beam_current_request",
        )
        raw_current = condition_value(
            provider,
            assignment.run,
            "beam_current",
        )
        raw_torus = condition_value(
            provider,
            assignment.run,
            "torus_scale",
        )

        request_na = to_float(raw_request)
        beam_current_na = to_float(raw_current)
        torus_scale = to_float(raw_torus)

        status, notes = classify_result(
            assignment=assignment,
            beam_current_request_na=request_na,
            beam_current_na=beam_current_na,
            torus_scale=torus_scale,
        )

        return RCDBResult(
            period=assignment.period,
            run=assignment.run,
            assigned_current_na=assignment.assigned_current_na,
            expected_torus_sign=assignment.expected_torus_sign,
            beam_current_request_na=request_na,
            beam_current_na=beam_current_na,
            torus_scale=torus_scale,
            status=status,
            notes=notes,
        )

    except Exception as exc:
        return RCDBResult(
            period=assignment.period,
            run=assignment.run,
            assigned_current_na=assignment.assigned_current_na,
            expected_torus_sign=assignment.expected_torus_sign,
            beam_current_request_na=None,
            beam_current_na=None,
            torus_scale=None,
            status="RCDB_ERROR",
            notes=f"{type(exc).__name__}: {exc}",
        )
    #endtry


def format_optional_number(
    value: Optional[float],
    decimals: int = 4,
) -> str:
    """Format a numeric value for CSV/text output."""

    if value is None:
        return ""
    #endif

    return f"{value:.{decimals}f}"


def format_compact_number(value: Optional[float]) -> str:
    """Format a number compactly for LaTeX tables."""

    if value is None:
        return "--"
    #endif

    if abs(value - round(value)) < 1.0e-9:
        return str(int(round(value)))
    #endif

    return f"{value:.4f}".rstrip("0").rstrip(".")


def latex_escape(text: str) -> str:
    """Escape ordinary text for use in a LaTeX table cell."""

    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }

    return "".join(replacements.get(character, character) for character in text)


def write_csv(
    path: Path,
    results: Iterable[RCDBResult],
) -> None:
    """Write run-validation results to a CSV file."""

    fieldnames = [
        "period",
        "run",
        "assigned_current_na",
        "beam_current_request_na",
        "beam_current_na",
        "torus_scale",
        "expected_torus_sign",
        "status",
        "notes",
    ]

    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()

        for result in results:
            writer.writerow(
                {
                    "period": result.period,
                    "run": result.run,
                    "assigned_current_na": format_optional_number(
                        result.assigned_current_na,
                        decimals=3,
                    ),
                    "beam_current_request_na": format_optional_number(
                        result.beam_current_request_na,
                        decimals=3,
                    ),
                    "beam_current_na": format_optional_number(
                        result.beam_current_na,
                        decimals=6,
                    ),
                    "torus_scale": format_optional_number(
                        result.torus_scale,
                        decimals=6,
                    ),
                    "expected_torus_sign": result.expected_torus_sign,
                    "status": result.status,
                    "notes": result.notes,
                }
            )
        #endfor
    #endwith


def result_requires_attention(result: RCDBResult) -> bool:
    """
    Return True for anything other than an ordinary exact/near-exact match.

    Manual overrides are included because they should be documented even
    though they are intentional.
    """

    return result.status != "MATCH"


def chunked(
    values: list[RCDBResult],
    chunk_size: int,
) -> Iterable[list[RCDBResult]]:
    """Yield consecutive fixed-size chunks."""

    for start in range(0, len(values), chunk_size):
        yield values[start : start + chunk_size]
    #endfor


def write_validation_latex(
    path: Path,
    attention_results: list[RCDBResult],
) -> None:
    """Write a compact LaTeX table of all entries requiring attention."""

    lines: list[str] = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\renewcommand{\arraystretch}{1.25}",
        r"\begin{tabular}{|c|c|c|c|c|c|c|}",
        r"\hline",
        (
            r"Period & Run & Assigned (nA) & RCDB request (nA) "
            r"& RCDB current (nA) & Torus & Status \\"
        ),
        r"\hline",
    ]

    if not attention_results:
        lines.extend(
            [
                (
                    r"\multicolumn{7}{|c|}{All run assignments agree with "
                    r"RCDB and no missing values were found.} \\"
                ),
                r"\hline",
            ]
        )
    else:
        for result in attention_results:
            torus_text = format_compact_number(result.torus_scale)

            lines.extend(
                [
                    (
                        f"{latex_escape(result.period)} & "
                        f"{result.run} & "
                        f"{format_compact_number(result.assigned_current_na)} & "
                        f"{format_compact_number(result.beam_current_request_na)} & "
                        f"{format_compact_number(result.beam_current_na)} & "
                        f"{torus_text} & "
                        rf"\texttt{{{latex_escape(result.status)}}} \\"
                    ),
                    r"\hline",
                ]
            )
        #endfor
    #endif

    lines.extend(
        [
            r"\end{tabular}",
            (
                r"\caption{RCDB validation results requiring attention. "
                r"Entries marked \texttt{MANUAL\_OVERRIDE} correspond to "
                r"intentional analysis assignments that differ from the "
                r"RCDB beam-current request.}"
            ),
            r"\label{tab:rcdb_current_validation}",
            r"\end{table}",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")


def write_full_latex_tables(
    path: Path,
    results: list[RCDBResult],
    rows_per_table: int,
) -> None:
    """
    Write complete run-by-run tables using ordinary table/tabular environments.

    Each period is split into tables with at most rows_per_table entries.
    This deliberately avoids longtable.
    """

    period_order = [
        "Sp18 Inb",
        "Sp18 Out",
        "Fa18 Inb",
        "Fa18 Out",
        "Sp19 Inb",
    ]

    lines: list[str] = [
        "% Automatically generated by check_rcdb_currents.py.",
        "% This file deliberately uses ordinary table environments rather",
        "% than longtable.",
        "",
    ]

    for period in period_order:
        period_results = [
            result for result in results if result.period == period
        ]
        period_results.sort(key=lambda result: result.run)

        number_of_parts = math.ceil(len(period_results) / rows_per_table)

        for part_index, part in enumerate(
            chunked(period_results, rows_per_table),
            start=1,
        ):
            safe_label_period = period.lower().replace(" ", "_")
            continuation = (
                ""
                if number_of_parts == 1
                else f", part {part_index} of {number_of_parts}"
            )

            lines.extend(
                [
                    r"\begin{table}[htbp]",
                    r"\centering",
                    r"\scriptsize",
                    r"\setlength{\tabcolsep}{4pt}",
                    r"\renewcommand{\arraystretch}{1.15}",
                    r"\begin{tabular}{|c|c|c|c|c|c|}",
                    r"\hline",
                    (
                        r"Run & Assigned (nA) & RCDB request (nA) "
                        r"& RCDB current (nA) & Torus & Status \\"
                    ),
                    r"\hline",
                ]
            )

            for result in part:
                lines.extend(
                    [
                        (
                            f"{result.run} & "
                            f"{format_compact_number(result.assigned_current_na)} & "
                            f"{format_compact_number(result.beam_current_request_na)} & "
                            f"{format_compact_number(result.beam_current_na)} & "
                            f"{format_compact_number(result.torus_scale)} & "
                            rf"\texttt{{{latex_escape(result.status)}}} \\"
                        ),
                        r"\hline",
                    ]
                )
            #endfor

            lines.extend(
                [
                    r"\end{tabular}",
                    (
                        rf"\caption{{RCDB beam-current and torus conditions "
                        rf"for the RGA {latex_escape(period)} analysis runs"
                        rf"{continuation}.}}"
                    ),
                    (
                        rf"\label{{tab:rcdb_{safe_label_period}_"
                        rf"conditions_part_{part_index}}}"
                    ),
                    r"\end{table}",
                    "",
                ]
            )
        #endfor
    #endfor

    path.write_text("\n".join(lines), encoding="utf-8")


def write_summary(
    path: Path,
    results: list[RCDBResult],
) -> None:
    """Write a human-readable summary of the validation."""

    status_counts: dict[str, int] = {}

    for result in results:
        status_counts[result.status] = status_counts.get(result.status, 0) + 1
    #endfor

    lines: list[str] = [
        "RCDB beam-current validation summary",
        "=" * 38,
        "",
        f"Total runs checked: {len(results)}",
        "",
        "Status counts:",
    ]

    for status in sorted(status_counts):
        lines.append(f"  {status:24s} {status_counts[status]:5d}")
    #endfor

    lines.extend(
        [
            "",
            "Entries requiring attention:",
            "-" * 28,
        ]
    )

    attention_results = [
        result for result in results if result_requires_attention(result)
    ]

    if not attention_results:
        lines.append("None.")
    else:
        for result in attention_results:
            lines.append(
                f"run {result.run:5d}  "
                f"{result.period:9s}  "
                f"assigned={result.assigned_current_na:5.1f} nA  "
                f"request={format_compact_number(result.beam_current_request_na):>7s}  "
                f"current={format_compact_number(result.beam_current_na):>9s}  "
                f"torus={format_compact_number(result.torus_scale):>6s}  "
                f"status={result.status}"
            )

            if result.notes:
                lines.append(f"    {result.notes}")
            #endif
        #endfor
    #endif

    lines.extend(
        [
            "",
            "Documented manual assignments:",
            "  runs 5340--5344 -> 40 nA",
            "  run 5615        -> 50 nA",
            "",
            "Run 3218 is intentionally absent from this validation because it",
            "is removed from the final current-dependence study.",
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line options."""

    parser = argparse.ArgumentParser(
        description=(
            "Validate DVCS run-current assignments against the CLAS12 RCDB."
        )
    )

    parser.add_argument(
        "--connection",
        default=DEFAULT_CONNECTION,
        help=(
            "RCDB connection string "
            f"(default: {DEFAULT_CONNECTION})"
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIRECTORY,
        help=(
            "Directory for CSV, text and LaTeX outputs "
            f"(default: {DEFAULT_OUTPUT_DIRECTORY})"
        ),
    )

    parser.add_argument(
        "--query-delay",
        type=float,
        default=DEFAULT_QUERY_DELAY_SECONDS,
        help=(
            "Delay in seconds after each run query "
            f"(default: {DEFAULT_QUERY_DELAY_SECONDS})"
        ),
    )

    parser.add_argument(
        "--latex-rows-per-table",
        type=int,
        default=DEFAULT_LATEX_ROWS_PER_TABLE,
        help=(
            "Maximum rows in each ordinary full LaTeX table "
            f"(default: {DEFAULT_LATEX_ROWS_PER_TABLE})"
        ),
    )

    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop immediately if an RCDB query fails.",
    )

    return parser.parse_args()


def main() -> int:
    """Program entry point."""

    args = parse_arguments()

    if args.query_delay < 0.0:
        print("ERROR: --query-delay must be nonnegative.", file=sys.stderr)
        return 2
    #endif

    if args.latex_rows_per_table < 1:
        print(
            "ERROR: --latex-rows-per-table must be at least 1.",
            file=sys.stderr,
        )
        return 2
    #endif

    try:
        assignments = build_assignments()
    except ValueError as exc:
        print(f"ERROR in embedded run assignments: {exc}", file=sys.stderr)
        return 2
    #endtry

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Connecting to RCDB:")
    print(f"  {args.connection}")
    print()
    print(f"Runs to check: {len(assignments)}")
    print()

    try:
        provider = rcdb.RCDBProvider(args.connection)
    except Exception as exc:
        print(
            "ERROR: Could not connect to RCDB.\n"
            f"{type(exc).__name__}: {exc}\n\n"
            "Confirm that you are running on a JLab machine and have sourced "
            "the CLAS12 environment.",
            file=sys.stderr,
        )
        return 1
    #endtry

    results: list[RCDBResult] = []

    for index, assignment in enumerate(assignments, start=1):
        print(
            f"[{index:4d}/{len(assignments):4d}] "
            f"run {assignment.run} "
            f"({assignment.period}, "
            f"assigned {assignment.assigned_current_na:g} nA)",
            end="",
            flush=True,
        )

        result = query_one_run(provider, assignment)
        results.append(result)

        print(
            f" -> request "
            f"{format_compact_number(result.beam_current_request_na)} nA, "
            f"beam current "
            f"{format_compact_number(result.beam_current_na)} nA, "
            f"torus {format_compact_number(result.torus_scale)}, "
            f"{result.status}"
        )

        if args.stop_on_error and result.status == "RCDB_ERROR":
            print(
                "\nStopping because --stop-on-error was requested.",
                file=sys.stderr,
            )
            break
        #endif

        if args.query_delay > 0.0:
            time.sleep(args.query_delay)
        #endif
    #endfor

    results.sort(key=lambda result: result.run)

    attention_results = [
        result for result in results if result_requires_attention(result)
    ]

    full_csv_path = args.output_dir / "rcdb_run_conditions.csv"
    discrepancy_csv_path = (
        args.output_dir / "rcdb_current_discrepancies.csv"
    )
    summary_path = args.output_dir / "rcdb_current_summary.txt"
    validation_tex_path = (
        args.output_dir / "rcdb_current_validation_table.tex"
    )
    full_tex_path = (
        args.output_dir / "rcdb_run_conditions_tables.tex"
    )

    write_csv(full_csv_path, results)
    write_csv(discrepancy_csv_path, attention_results)
    write_summary(summary_path, results)
    write_validation_latex(validation_tex_path, attention_results)
    write_full_latex_tables(
        full_tex_path,
        results,
        rows_per_table=args.latex_rows_per_table,
    )

    print()
    print("Validation complete.")
    print()
    print("Generated files:")
    print(f"  {full_csv_path}")
    print(f"  {discrepancy_csv_path}")
    print(f"  {summary_path}")
    print(f"  {validation_tex_path}")
    print(f"  {full_tex_path}")
    print()

    ordinary_matches = sum(
        result.status == "MATCH" for result in results
    )
    manual_overrides = sum(
        result.status == "MANUAL_OVERRIDE" for result in results
    )
    other_attention = len(results) - ordinary_matches - manual_overrides

    print(f"Ordinary matches:       {ordinary_matches}")
    print(f"Manual overrides:       {manual_overrides}")
    print(f"Other entries to check: {other_attention}")

    if other_attention > 0:
        print()
        print(
            "WARNING: At least one unexpected mismatch, missing value, "
            "torus problem or RCDB error was found."
        )
        return 3
    #endif

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
#endif