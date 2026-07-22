#!/usr/bin/env python3
"""
Report the percentage of e pi+ events with 8 <= e_theta <= 14 degrees
for the four RGC calibration periods.

The input files and run-period definitions match the momentum-correction
scripts. Fa22 uses one ROOT file and is split by runnum into the two solenoid
configurations.
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

INPUT_FILES = {
    "su22": ROOT_DIR / "rgc_su22_inb_NH3_epi+.root",
    "fa22": ROOT_DIR / "rgc_fa22_inb_NH3_epi+.root",
    "sp23": ROOT_DIR / "rgc_sp23_inb_NH3_epi+.root",
}

TREE_NAME = "PhysicsEvents"
THETA_BRANCH = "e_theta"
RUN_BRANCH = "runnum"

THETA_MIN_DEG = 8.0
THETA_MAX_DEG = 14.0

N_WORKERS = 4
STEP_SIZE = "250 MB"
INPUT_UNITS = "auto"


@dataclass(frozen=True)
class Period:
    label: str
    source_key: str
    run_min: int
    run_max: int


PERIODS = (
    Period("Su22", "su22", 16042, 16788),
    Period("Fa22, solenoid -1", "fa22", 16843, 17183),
    Period("Fa22, solenoid +1", "fa22", 17185, 17408),
    Period("Sp23", "sp23", 17477, 17811),
)


def validate_inputs() -> None:
    missing = [path for path in INPUT_FILES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing required input file(s):\n"
            + "\n".join(f"  {path}" for path in missing)
        )


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


def analyze_period(period: Period) -> dict:
    path = INPUT_FILES[period.source_key]
    source = f"{path}:{TREE_NAME}"

    total_in_file = 0
    total_in_period = 0
    inside = 0
    nonfinite = 0

    run_seen_min = math.inf
    run_seen_max = -math.inf
    raw_min = math.inf
    raw_max = -math.inf
    deg_min = math.inf
    deg_max = -math.inf

    diagnostic_samples = []
    diagnostic_count = 0
    units = None

    for arrays in uproot.iterate(
        source,
        expressions=[RUN_BRANCH, THETA_BRANCH],
        step_size=STEP_SIZE,
        library="np",
    ):
        runnum = np.asarray(arrays[RUN_BRANCH], dtype=np.int64)
        theta_raw = np.asarray(arrays[THETA_BRANCH], dtype=np.float64)

        if runnum.size != theta_raw.size:
            raise RuntimeError(
                f"{period.label}: runnum and e_theta lengths differ."
            )

        total_in_file += int(runnum.size)

        if runnum.size:
            run_seen_min = min(run_seen_min, int(np.min(runnum)))
            run_seen_max = max(run_seen_max, int(np.max(runnum)))

        period_mask = (
            (runnum >= period.run_min)
            & (runnum <= period.run_max)
        )
        theta_raw = theta_raw[period_mask]
        total_in_period += int(theta_raw.size)

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

    finite_total = total_in_period - nonfinite

    if total_in_period == 0:
        raise RuntimeError(
            f"{period.label}: no events found in run range "
            f"{period.run_min}-{period.run_max}. File run range was "
            f"{run_seen_min}-{run_seen_max}."
        )

    if finite_total == 0:
        raise RuntimeError(
            f"{period.label}: all selected e_theta values are nonfinite."
        )

    sample = np.concatenate(diagnostic_samples)
    q01, q50, q99 = np.percentile(sample, [1.0, 50.0, 99.0])

    return {
        "label": period.label,
        "file": str(path),
        "run_min": period.run_min,
        "run_max": period.run_max,
        "file_entries": total_in_file,
        "selected_entries": total_in_period,
        "finite_entries": finite_total,
        "inside": inside,
        "percent": 100.0 * inside / finite_total,
        "nonfinite": nonfinite,
        "units": units,
        "file_run_min": int(run_seen_min),
        "file_run_max": int(run_seen_max),
        "raw_min": raw_min,
        "raw_max": raw_max,
        "deg_min": deg_min,
        "deg_max": deg_max,
        "q01": float(q01),
        "q50": float(q50),
        "q99": float(q99),
    }


def main() -> None:
    validate_inputs()

    results = {}

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        future_map = {
            executor.submit(analyze_period, period): period.label
            for period in PERIODS
        }

        for future in as_completed(future_map):
            label = future_map[future]
            results[label] = future.result()

    print(
        f"\nFraction of PhysicsEvents entries with "
        f"{THETA_MIN_DEG:g} <= e_theta <= {THETA_MAX_DEG:g} degrees\n"
    )
    print(
        f"{'Period':<24} {'Inside':>14} {'Finite total':>14} "
        f"{'Percent':>12} {'Input':>8}"
    )
    print("-" * 82)

    combined_inside = 0
    combined_total = 0

    for period in PERIODS:
        result = results[period.label]
        combined_inside += result["inside"]
        combined_total += result["finite_entries"]

        print(
            f"{period.label:<24} "
            f"{result['inside']:>14,d} "
            f"{result['finite_entries']:>14,d} "
            f"{result['percent']:>11.3f}% "
            f"{result['units']:>8}"
        )

    print("-" * 82)
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
        print(
            f"  Requested run range: "
            f"{result['run_min']}-{result['run_max']}"
        )
        print(
            f"  Full file run range: "
            f"{result['file_run_min']}-{result['file_run_max']}"
        )
        print(f"  Full file entries: {result['file_entries']:,}")
        print(
            f"  Entries in requested run range: "
            f"{result['selected_entries']:,}"
        )
        print(f"  Nonfinite e_theta entries: {result['nonfinite']:,}")
        print(f"  Inferred e_theta units: {result['units']}")
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


if __name__ == "__main__":
    main()
