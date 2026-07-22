#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import math

import numpy as np
import uproot


ROOT_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+"
)

PERIOD_PATTERNS = {
    "Su22 Inb": "*su22*inb*NH3_epi+.root",
    "Fa22 Inb": "*fa22*inb*NH3_epi+.root",
    "Fa22 Out": "*fa22*out*NH3_epi+.root",
    "Sp23 Inb": "*sp23*inb*NH3_epi+.root",
}

TREE_NAME = "PhysicsEvents"
BRANCH_NAME = "e_theta"

THETA_MIN_DEG = 8.0
THETA_MAX_DEG = 14.0

N_WORKERS = 4
STEP_SIZE = "250 MB"

# "auto" detects radians versus degrees from the branch values.
# Set explicitly to "rad" or "deg" if desired.
INPUT_UNITS = "auto"


def find_one_file(pattern: str) -> Path:
    matches = sorted(ROOT_DIR.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one file matching {pattern!r}, "
            f"but found {len(matches)}:\n"
            + "\n".join(f"  {path}" for path in matches)
        )
    return matches[0]


def infer_units(sample: np.ndarray) -> str:
    """
    Infer whether e_theta is stored in radians or degrees.

    CLAS12 electron polar angles are normally O(0.1--0.5) in radians
    or O(5--30) in degrees. A 99th percentile below 2*pi is therefore
    treated as radians.
    """
    finite = sample[np.isfinite(sample)]
    if finite.size == 0:
        raise RuntimeError("Cannot infer units from an empty/nonfinite sample.")

    q99 = float(np.percentile(np.abs(finite), 99.0))
    return "rad" if q99 <= 2.0 * math.pi else "deg"


def convert_to_degrees(theta: np.ndarray, units: str) -> np.ndarray:
    if units == "rad":
        return np.degrees(theta)
    if units == "deg":
        return theta
    raise ValueError(f"Unsupported units: {units!r}")


def analyze_period(period: str, path_str: str) -> dict:
    path = Path(path_str)
    source = f"{path}:{TREE_NAME}"

    total_finite = 0
    inside = 0

    raw_min = math.inf
    raw_max = -math.inf
    deg_min = math.inf
    deg_max = -math.inf

    sample_chunks = []
    units = None

    for arrays in uproot.iterate(
        source,
        expressions=[BRANCH_NAME],
        step_size=STEP_SIZE,
        library="np",
    ):
        raw = np.asarray(arrays[BRANCH_NAME], dtype=np.float64)
        raw = raw[np.isfinite(raw)]

        if raw.size == 0:
            continue

        if units is None:
            if INPUT_UNITS == "auto":
                units = infer_units(raw[: min(raw.size, 1_000_000)])
            else:
                units = INPUT_UNITS

        theta_deg = convert_to_degrees(raw, units)

        total_finite += raw.size
        inside += int(
            np.count_nonzero(
                (theta_deg >= THETA_MIN_DEG)
                & (theta_deg <= THETA_MAX_DEG)
            )
        )

        raw_min = min(raw_min, float(np.min(raw)))
        raw_max = max(raw_max, float(np.max(raw)))
        deg_min = min(deg_min, float(np.min(theta_deg)))
        deg_max = max(deg_max, float(np.max(theta_deg)))

        # Keep a bounded diagnostic sample for percentiles.
        if sum(chunk.size for chunk in sample_chunks) < 1_000_000:
            remaining = 1_000_000 - sum(chunk.size for chunk in sample_chunks)
            sample_chunks.append(theta_deg[:remaining])

    if total_finite == 0:
        raise RuntimeError(f"{period}: no finite {BRANCH_NAME} values found.")

    diagnostic_sample = np.concatenate(sample_chunks)
    q01, q50, q99 = np.percentile(diagnostic_sample, [1.0, 50.0, 99.0])

    return {
        "period": period,
        "path": str(path),
        "units": units,
        "inside": inside,
        "total": total_finite,
        "percent": 100.0 * inside / total_finite,
        "raw_min": raw_min,
        "raw_max": raw_max,
        "deg_min": deg_min,
        "deg_max": deg_max,
        "deg_q01": float(q01),
        "deg_q50": float(q50),
        "deg_q99": float(q99),
    }


def main() -> None:
    jobs = []
    for period, pattern in PERIOD_PATTERNS.items():
        path = find_one_file(pattern)
        jobs.append((period, path))

    results = {}

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        future_map = {
            executor.submit(analyze_period, period, str(path)): period
            for period, path in jobs
        }

        for future in as_completed(future_map):
            period = future_map[future]
            results[period] = future.result()

    print(
        f"\nFraction of PhysicsEvents entries with "
        f"{THETA_MIN_DEG:g} <= e_theta <= {THETA_MAX_DEG:g} degrees\n"
    )
    print(f"{'Period':<12} {'Inside':>14} {'Total':>14} {'Percent':>12} {'Input':>8}")
    print("-" * 68)

    combined_inside = 0
    combined_total = 0

    for period in PERIOD_PATTERNS:
        result = results[period]
        combined_inside += result["inside"]
        combined_total += result["total"]

        print(
            f"{period:<12} "
            f"{result['inside']:>14,d} "
            f"{result['total']:>14,d} "
            f"{result['percent']:>11.3f}% "
            f"{result['units']:>8}"
        )

    combined_percent = 100.0 * combined_inside / combined_total

    print("-" * 68)
    print(
        f"{'Combined':<12} "
        f"{combined_inside:>14,d} "
        f"{combined_total:>14,d} "
        f"{combined_percent:>11.3f}%"
    )

    print("\nDiagnostics")
    print("-----------")
    for period in PERIOD_PATTERNS:
        r = results[period]
        print(f"\n{period}")
        print(f"  File: {r['path']}")
        print(f"  Inferred input units: {r['units']}")
        print(
            f"  Raw range: {r['raw_min']:.6g} to {r['raw_max']:.6g}"
        )
        print(
            f"  Degree range: {r['deg_min']:.3f} to {r['deg_max']:.3f}"
        )
        print(
            "  Degree percentiles "
            f"(1%, 50%, 99%): "
            f"{r['deg_q01']:.3f}, "
            f"{r['deg_q50']:.3f}, "
            f"{r['deg_q99']:.3f}"
        )


if __name__ == "__main__":
    main()
