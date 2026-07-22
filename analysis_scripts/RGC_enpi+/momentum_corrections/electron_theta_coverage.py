#!/usr/bin/env python3
"""
Report the percentage of e pi+ events with 8 <= e_theta <= 14 degrees
for the four RGC calibration periods.

Fa22 is stored in one ROOT file. Its two solenoid configurations are separated
using runnum:
    Fa22, solenoid -1: 16843-17183
    Fa22, solenoid +1: 17185-17408

The four logical periods are processed concurrently with four workers.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
import math

import numpy as np
import uproot


ROOT_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+"
)

TREE_NAME = "PhysicsEvents"
THETA_BRANCH = "e_theta"
RUN_BRANCH = "runnum"

THETA_MIN_DEG = 8.0
THETA_MAX_DEG = 14.0

N_WORKERS = 4
STEP_SIZE = "250 MB"

# e_theta is expected to be in radians in the analysis trees.
# Keep "auto" to verify this from the data and convert accordingly.
INPUT_UNITS = "auto"


@dataclass(frozen=True)
class Period:
    label: str
    source_key: str
    run_min: int | None = None
    run_max: int | None = None


PERIODS = (
    Period("Su22", "su22"),
    Period("Fa22, solenoid -1", "fa22", 16843, 17183),
    Period("Fa22, solenoid +1", "fa22", 17185, 17408),
    Period("Sp23", "sp23"),
)


# Several plausible filename forms are accepted. Exactly one file must be found
# for each physical source period.
SOURCE_PATTERNS = {
    "su22": (
        "*su22*NH3*epi+*.root",
        "*Su22*NH3*epi+*.root",
        "*su22*epi+*.root",
    ),
    "fa22": (
        "*fa22*NH3*epi+*.root",
        "*Fa22*NH3*epi+*.root",
        "*fa22*epi+*.root",
    ),
    "sp23": (
        "*sp23*NH3*epi+*.root",
        "*Sp23*NH3*epi+*.root",
        "*sp23*epi+*.root",
    ),
}


def find_source_file(source_key: str) -> Path:
    matches: set[Path] = set()

    for pattern in SOURCE_PATTERNS[source_key]:
        matches.update(ROOT_DIR.glob(pattern))

    matches = {path.resolve() for path in matches if path.is_file()}

    if len(matches) != 1:
        formatted = "\n".join(f"  {path}" for path in sorted(matches))
        patterns = "\n".join(
            f"  {pattern}" for pattern in SOURCE_PATTERNS[source_key]
        )
        raise RuntimeError(
            f"Expected exactly one {source_key} e pi+ ROOT file in\n"
            f"  {ROOT_DIR}\n"
            f"using patterns:\n{patterns}\n"
            f"but found {len(matches)} unique files:\n{formatted}"
        )

    return next(iter(matches))


def infer_units(values: np.ndarray) -> str:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise RuntimeError("Cannot infer e_theta units from an empty sample.")

    q99 = float(np.percentile(np.abs(finite), 99.0))
    return "rad" if q99 <= 2.0 * math.pi else "deg"


def to_degrees(values: np.ndarray, units: str) -> np.ndarray:
    if units == "rad":
        return np.degrees(values)
    if units == "deg":
        return values
    raise ValueError(f"Unsupported angular units: {units}")


def analyze_period(period: Period, source_path: str) -> dict:
    path = Path(source_path)
    expressions = [THETA_BRANCH]

    if period.run_min is not None:
        expressions.append(RUN_BRANCH)

    total_selected = 0
    inside = 0
    rejected_by_run = 0
    nonfinite = 0

    raw_min = math.inf
    raw_max = -math.inf
    deg_min = math.inf
    deg_max = -math.inf

    diagnostic_samples: list[np.ndarray] = []
    diagnostic_count = 0
    units: str | None = None

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=STEP_SIZE,
        library="np",
    ):
        theta_raw = np.asarray(arrays[THETA_BRANCH], dtype=np.float64)

        if period.run_min is not None:
            runnum = np.asarray(arrays[RUN_BRANCH], dtype=np.int64)
            run_mask = (
                (runnum >= period.run_min)
                & (runnum <= period.run_max)
            )
            rejected_by_run += int(np.count_nonzero(~run_mask))
            theta_raw = theta_raw[run_mask]

        finite_mask = np.isfinite(theta_raw)
        nonfinite += int(np.count_nonzero(~finite_mask))
        theta_raw = theta_raw[finite_mask]

        if theta_raw.size == 0:
            continue

        if units is None:
            units = (
                infer_units(theta_raw[: min(theta_raw.size, 1_000_000)])
                if INPUT_UNITS == "auto"
                else INPUT_UNITS
            )

        theta_deg = to_degrees(theta_raw, units)

        total_selected += int(theta_deg.size)
        inside += int(
            np.count_nonzero(
                (theta_deg >= THETA_MIN_DEG)
                & (theta_deg <= THETA_MAX_DEG)
            )
        )

        raw_min = min(raw_min, float(np.min(theta_raw)))
        raw_max = max(raw_max, float(np.max(theta_raw)))
        deg_min = min(deg_min, float(np.min(theta_deg)))
        deg_max = max(deg_max, float(np.max(theta_deg)))

        if diagnostic_count < 1_000_000:
            take = min(theta_deg.size, 1_000_000 - diagnostic_count)
            diagnostic_samples.append(theta_deg[:take])
            diagnostic_count += take

    if total_selected == 0:
        raise RuntimeError(
            f"{period.label}: no finite events remained after period selection."
        )

    sample = np.concatenate(diagnostic_samples)
    q01, q50, q99 = np.percentile(sample, [1.0, 50.0, 99.0])

    return {
        "label": period.label,
        "file": str(path),
        "run_min": period.run_min,
        "run_max": period.run_max,
        "units": units,
        "inside": inside,
        "total": total_selected,
        "percent": 100.0 * inside / total_selected,
        "rejected_by_run": rejected_by_run,
        "nonfinite": nonfinite,
        "raw_min": raw_min,
        "raw_max": raw_max,
        "deg_min": deg_min,
        "deg_max": deg_max,
        "q01": float(q01),
        "q50": float(q50),
        "q99": float(q99),
    }


def main() -> None:
    source_files = {
        key: find_source_file(key)
        for key in SOURCE_PATTERNS
    }

    results: dict[str, dict] = {}

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        futures = {
            executor.submit(
                analyze_period,
                period,
                str(source_files[period.source_key]),
            ): period.label
            for period in PERIODS
        }

        for future in as_completed(futures):
            label = futures[future]
            results[label] = future.result()

    print(
        f"\nFraction of PhysicsEvents entries with "
        f"{THETA_MIN_DEG:g} <= e_theta <= {THETA_MAX_DEG:g} degrees\n"
    )
    print(
        f"{'Period':<24} {'Inside':>14} {'Total':>14} "
        f"{'Percent':>12} {'Input':>8}"
    )
    print("-" * 78)

    combined_inside = 0
    combined_total = 0

    for period in PERIODS:
        result = results[period.label]
        combined_inside += result["inside"]
        combined_total += result["total"]

        print(
            f"{period.label:<24} "
            f"{result['inside']:>14,d} "
            f"{result['total']:>14,d} "
            f"{result['percent']:>11.3f}% "
            f"{result['units']:>8}"
        )

    print("-" * 78)
    print(
        f"{'Combined':<24} "
        f"{combined_inside:>14,d} "
        f"{combined_total:>14,d} "
        f"{100.0 * combined_inside / combined_total:>11.3f}%"
    )

    print("\nDiagnostics")
    print("-----------")

    for period in PERIODS:
        result = results[period.label]

        print(f"\n{period.label}")
        print(f"  File: {result['file']}")

        if result["run_min"] is not None:
            print(
                f"  Run selection: "
                f"{result['run_min']}-{result['run_max']}"
            )
            print(
                f"  Entries outside run range: "
                f"{result['rejected_by_run']:,}"
            )

        print(f"  Inferred input units: {result['units']}")
        print(
            f"  Raw e_theta range: "
            f"{result['raw_min']:.6g} to {result['raw_max']:.6g}"
        )
        print(
            f"  Converted degree range: "
            f"{result['deg_min']:.3f} to {result['deg_max']:.3f}"
        )
        print(
            f"  Degree percentiles (1%, 50%, 99%): "
            f"{result['q01']:.3f}, "
            f"{result['q50']:.3f}, "
            f"{result['q99']:.3f}"
        )
        print(f"  Nonfinite selected entries: {result['nonfinite']:,}")


if __name__ == "__main__":
    main()
