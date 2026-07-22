#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import uproot

ROOT_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+"
)

# Adjust only the Fa22-out pattern if its filename uses a different spelling.
PERIOD_PATTERNS = {
    "Su22 Inb": "*su22*inb*NH3_epi+.root",
    "Fa22 Inb": "*fa22*inb*NH3_epi+.root",
    "Fa22 Out": "*fa22*out*NH3_epi+.root",
    "Sp23 Inb": "*sp23*inb*NH3_epi+.root",
}

TREE_NAME = "PhysicsEvents"
BRANCH_NAME = "e_theta"
THETA_MIN = 8.0
THETA_MAX = 14.0


def find_one_file(pattern: str) -> Path:
    matches = sorted(ROOT_DIR.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one file matching {pattern!r}, "
            f"but found {len(matches)}:\n"
            + "\n".join(f"  {path}" for path in matches)
        )
    return matches[0]


def count_events(path: Path) -> tuple[int, int]:
    total = 0
    inside = 0

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=[BRANCH_NAME],
        step_size="250 MB",
        library="np",
    ):
        theta = np.asarray(arrays[BRANCH_NAME], dtype=np.float64)
        theta = theta[np.isfinite(theta)]

        total += theta.size
        inside += np.count_nonzero(
            (theta >= THETA_MIN) & (theta <= THETA_MAX)
        )

    return total, inside


def main() -> None:
    combined_total = 0
    combined_inside = 0

    print(
        f"\nFraction of PhysicsEvents entries with "
        f"{THETA_MIN:g} <= e_theta <= {THETA_MAX:g} degrees\n"
    )
    print(f"{'Period':<12} {'Inside':>14} {'Total':>14} {'Percent':>12}")
    print("-" * 56)

    for period, pattern in PERIOD_PATTERNS.items():
        path = find_one_file(pattern)
        total, inside = count_events(path)
        percent = 100.0 * inside / total if total else float("nan")

        combined_total += total
        combined_inside += inside

        print(
            f"{period:<12} "
            f"{inside:>14,d} "
            f"{total:>14,d} "
            f"{percent:>11.3f}%"
        )

    combined_percent = (
        100.0 * combined_inside / combined_total
        if combined_total
        else float("nan")
    )

    print("-" * 56)
    print(
        f"{'Combined':<12} "
        f"{combined_inside:>14,d} "
        f"{combined_total:>14,d} "
        f"{combined_percent:>11.3f}%"
    )


if __name__ == "__main__":
    main()
