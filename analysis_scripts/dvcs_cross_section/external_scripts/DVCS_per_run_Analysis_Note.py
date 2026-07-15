#!/usr/bin/env python3

"""Create analysis-note-quality DVCS events-per-charge plots by run.

The script reads one ROOT file for each RGA period, obtains the event count per
run from PhysicsEvents/runnum, and normalizes those counts using the Faraday Cup
charge in import/Analysis_Note_charges.txt. Charges are read in nC and converted
to microcoulombs before calculating Events / microC.

Unlike DVCS_per_run.py, removed runs are defined explicitly below. They are
excluded from the per-current mean and displayed with an x marker. The script
performs a complete preflight validation and exits with a detailed error if any
ROOT-tree run lacks either a positive charge or a current assignment.
"""

import argparse
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


TREE_NAME = "PhysicsEvents"
ROOT_DIRECTORY = Path("/volatile/clas12/thayward/DVCS_per_run_for_analysis_note")
DEFAULT_CHARGE_FILE = Path("import/Analysis_Note_charges.txt")
DEFAULT_OUTPUT_DIR = Path("output/dvcs_per_run_analysis_note")

PERIOD_FILES = [
    ("rga_fa18_inb", ROOT_DIRECTORY / "rga_fa18_inb.root"),
    ("rga_fa18_out", ROOT_DIRECTORY / "rga_fa18_out.root"),
    ("rga_sp19_inb", ROOT_DIRECTORY / "rga_sp19_inb.root"),
    ("rga_sp18_inb", ROOT_DIRECTORY / "rga_sp18_inb.root"),
    ("rga_sp18_out", ROOT_DIRECTORY / "rga_sp18_out.root"),
]

PERIOD_DISPLAY_NAMES = {
    "rga_fa18_inb": "RGA Fa18 Inb",
    "rga_fa18_out": "RGA Fa18 Out",
    "rga_sp19_inb": "RGA Sp19 Inb",
    "rga_sp18_inb": "RGA Sp18 Inb",
    "rga_sp18_out": "RGA Sp18 Out",
}

REMOVED_RUNS = {
    "rga_fa18_inb": {5247, 5345},
    "rga_fa18_out": set(),
    "rga_sp19_inb": {6616, 6618},
    "rga_sp18_inb": {
        3355, 3404, 3408, 3449, 3490, 3499, 3500, 3505, 3508, 3526,
        3527, 3528, 3529, 3530, 3531, 3532, 3533, 3534, 3535, 3536,
        3538, 3540, 3544, 3545, 3547, 3548, 3698, 3709, 3712, 3736,
        3793, 3800, 3801, 3807, 3808, 3809, 3810, 3813, 3814, 3815,
        3817, 3867, 3877, 3879, 3882, 3923, 3927, 3929, 3947, 3951,
        3953, 3965, 3967, 3968, 4018, 4059, 4142, 4145, 4146, 4159,
        4160, 4162, 4163, 4176, 4209, 4227, 4246, 4252, 4325,
    },
    "rga_sp18_out": {
        3267, 3867, 3877, 3879, 3882, 3923, 3927, 3929, 3947, 3951,
        3953, 3965, 3967, 3968,
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Make publication-quality DVCS events/microC plots by run."
    )
    parser.add_argument(
        "--charge-file",
        type=Path,
        default=DEFAULT_CHARGE_FILE,
        help=(
            "Charge table containing run number in column 1 and Faraday Cup "
            "charge in nC in column 2. Default: import/Analysis_Note_charges.txt"
        ),
    )
    parser.add_argument(
        "--root-dir",
        type=Path,
        default=ROOT_DIRECTORY,
        help=(
            "Directory containing rga_fa18_inb.root, rga_fa18_out.root, "
            "rga_sp18_inb.root, rga_sp18_out.root and rga_sp19_inb.root."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory in which PNG plots are written.",
    )
    return parser.parse_args()
#enddef


FA18_INB_CURRENT = {
    # 5 nA
    5418: 5, 5419: 5,

    # 40 nA
    5335: 40, 5339: 40, 5341: 40,
    5340: 40, 5342: 40, 5343: 40, 5344: 40,

    # 45 nA
    5032: 45, 5036: 45, 5038: 45, 5039: 45, 5040: 45, 5041: 45,
    5043: 45, 5045: 45, 5052: 45, 5053: 45, 5116: 45, 5117: 45,
    5119: 45, 5120: 45, 5124: 45, 5125: 45, 5126: 45, 5127: 45,
    5139: 45, 5153: 45, 5158: 45, 5162: 45, 5163: 45, 5164: 45,
    5181: 45, 5191: 45, 5193: 45, 5195: 45, 5196: 45, 5197: 45,
    5198: 45, 5199: 45, 5200: 45, 5201: 45, 5202: 45, 5203: 45,
    5204: 45, 5205: 45, 5206: 45, 5208: 45, 5211: 45, 5212: 45,
    5215: 45, 5216: 45, 5219: 45, 5220: 45, 5221: 45, 5222: 45,
    5223: 45, 5230: 45, 5231: 45, 5232: 45, 5233: 45, 5234: 45,
    5235: 45, 5237: 45, 5238: 45, 5247: 45, 5248: 45, 5249: 45,
    5252: 45, 5253: 45, 5257: 45, 5258: 45, 5259: 45, 5261: 45,
    5262: 45, 5303: 45, 5304: 45, 5305: 45, 5306: 45, 5307: 45,
    5310: 45, 5311: 45, 5315: 45, 5317: 45, 5318: 45, 5319: 45,
    5320: 45, 5323: 45, 5324: 45, 5333: 45, 5334: 45, 5346: 45,
    5347: 45, 5349: 45, 5351: 45, 5354: 45, 5355: 45, 5367: 45,

    # Extra Fa18 Inb 45 nA runs
    5046: 45, 5047: 45, 5051: 45,
    5128: 45, 5129: 45, 5130: 45,
    5159: 45, 5160: 45,
    5165: 45, 5166: 45, 5167: 45, 5168: 45, 5169: 45,
    5180: 45, 5182: 45, 5183: 45,
    5190: 45,
    5239: 45,
    5336: 45,

    # 50 nA
    5356: 50, 5357: 50, 5358: 50, 5359: 50, 5360: 50, 5361: 50,
    5362: 50, 5366: 50,

    # 55 nA
    5368: 55, 5369: 55, 5372: 55, 5373: 55, 5374: 55, 5375: 55,
    5376: 55, 5377: 55, 5378: 55, 5379: 55, 5380: 55, 5381: 55,
    5382: 55, 5383: 55, 5386: 55, 5390: 55, 5391: 55, 5392: 55,
    5393: 55, 5398: 55, 5400: 55, 5401: 55, 5403: 55, 5404: 55,
    5406: 55, 5407: 55,
}

FA18_OUT_CURRENT = {
    # 5 nA
    5443: 5,

    # 20 nA
    5444: 20,

    # 40 nA
    5423: 40, 5424: 40, 5425: 40, 5426: 40, 5428: 40, 5429: 40,
    5430: 40, 5432: 40, 5434: 40, 5435: 40, 5436: 40, 5437: 40,
    5438: 40, 5440: 40, 5441: 40, 5442: 40, 5445: 40, 5447: 40,
    5449: 40, 5450: 40, 5451: 40, 5452: 40, 5453: 40, 5454: 40,
    5455: 40, 5460: 40, 5464: 40, 5465: 40, 5466: 40, 5467: 40,
    5468: 40, 5469: 40, 5470: 40, 5471: 40, 5472: 40, 5473: 40,
    5474: 40, 5475: 40, 5476: 40, 5478: 40, 5479: 40, 5480: 40,
    5481: 40, 5482: 40, 5483: 40, 5485: 40, 5486: 40, 5487: 40,
    5497: 40, 5498: 40, 5499: 40, 5500: 40, 5504: 40,

    # Extra Fa18 Out 40 nA runs
    5448: 40, 5495: 40, 5496: 40,

    # 50 nA
    5507: 50, 5516: 50, 5517: 50, 5518: 50, 5519: 50, 5520: 50,
    5521: 50, 5522: 50, 5523: 50, 5524: 50, 5525: 50, 5526: 50,
    5527: 50, 5528: 50, 5530: 50, 5532: 50, 5533: 50, 5534: 50,
    5535: 50, 5536: 50, 5537: 50, 5538: 50, 5540: 50, 5541: 50,
    5543: 50, 5544: 50, 5545: 50, 5546: 50, 5547: 50, 5548: 50,
    5549: 50, 5550: 50, 5551: 50, 5552: 50, 5555: 50, 5556: 50,
    5557: 50, 5558: 50, 5559: 50, 5562: 50, 5569: 50, 5570: 50,
    5571: 50, 5572: 50, 5573: 50, 5574: 50, 5577: 50, 5578: 50,
    5591: 50, 5592: 50, 5594: 50, 5597: 50, 5598: 50, 5600: 50,
    5601: 50, 5602: 50, 5603: 50, 5604: 50, 5606: 50, 5607: 50,
    5611: 50, 5612: 50, 5613: 50, 5614: 50, 5616: 50, 5618: 50,
    5619: 50, 5624: 50, 5625: 50, 5626: 50, 5627: 50, 5628: 50,
    5629: 50, 5630: 50, 5631: 50, 5632: 50, 5633: 50, 5635: 50,
    5637: 50, 5638: 50, 5639: 50, 5641: 50, 5643: 50, 5644: 50,
    5645: 50, 5646: 50, 5647: 50, 5648: 50, 5649: 50, 5650: 50,
    5651: 50, 5652: 50, 5654: 50, 5655: 50, 5656: 50, 5662: 50,
    5663: 50, 5664: 50, 5665: 50, 5666: 50,
    5615: 50,

    # Extra Fa18 Out 50 nA runs
    5505: 50, 5567: 50, 5617: 50, 5621: 50, 5623: 50,

    # Extra Fa18 Out 5 nA run
    5610: 5,
}

TRIGGER_RANGES_SP18 = [
    ("trigger v2-v7", 3031, 3495),
    ("trigger v8", 3495, 3517),
    ("trigger v9", 3517, 3548),
    ("trigger v10", 3709, 3722),
    ("trigger v11-v2_3", 3722, 4325),
]


def resolve_current(period_label, runnum):
    """
    Resolve beam current (nA) for a given period and run number.

    Returns (True, current_nA) if known, else (False, None).
    """
    label = period_label.lower()

    # Fa18 Inb
    if label == "rga_fa18_inb":
        if runnum in FA18_INB_CURRENT:
            return True, FA18_INB_CURRENT[runnum]
        #endif
        return False, None
    #endif

    # Fa18 Out
    if label == "rga_fa18_out":
        if runnum in FA18_OUT_CURRENT:
            return True, FA18_OUT_CURRENT[runnum]
        #endif
        return False, None
    #endif

    # Sp18 Out: 30 nA (3211-3293), 45 nA (3867-3987)
    if label == "rga_sp18_out":
        if 3211 <= runnum <= 3293:
            return True, 30
        #endif
        if 3867 <= runnum <= 3987:
            return True, 45
        #endif
        return False, None
    #endif

    # Sp18 Inb: overrides, then ranges
    if label == "rga_sp18_inb":
        if runnum == 3418:
            return True, 70
        #endif
        if runnum == 3421 or runnum == 3422:
            return True, 35
        #endif
        if runnum == 3429:
            return True, 50
        #endif

        # 35 nA from 3306-3411
        if 3306 <= runnum <= 3411:
            return True, 35
        #endif

        # 50 nA from 3431-4325
        if 3431 <= runnum <= 4325:
            return True, 50
        #endif

        return False, None
    #endif

    # Sp19 Inb: all 50 nA except one 5 nA run
    if label == "rga_sp19_inb":
        if runnum == 6616:
            return True, 5
        #endif
        if runnum == 6618:
            return True, 10
        #endif
        return True, 50
    #endif

    # Other periods: unknown
    return False, None
#enddef



def load_charge_map(charge_file):
    if not charge_file.is_file():
        raise RuntimeError(
            f"Charge file does not exist: {charge_file}\n"
            "Expected a comma-separated table with run number in column 1 and "
            "Faraday Cup charge in nC in column 2."
        )
    #endif

    charge_map = {}

    with charge_file.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue
            #endif

            fields = [field.strip() for field in line.split(",")]

            if len(fields) < 2:
                raise RuntimeError(
                    f"Malformed charge row at {charge_file}:{line_number}: {raw_line.rstrip()}\n"
                    "Every non-comment row must contain at least two comma-separated columns."
                )
            #endif

            try:
                run_number = int(fields[0])
                charge_nC = float(fields[1])
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid run number or charge at {charge_file}:{line_number}: "
                    f"{raw_line.rstrip()}"
                ) from exc
            #endtry

            if run_number in charge_map:
                raise RuntimeError(
                    f"Duplicate charge entry for run {run_number} in {charge_file} "
                    f"(second occurrence at line {line_number})."
                )
            #endif

            charge_map[run_number] = charge_nC / 1000.0
        #endfor
    #endwith

    if not charge_map:
        raise RuntimeError(f"No charge records were loaded from {charge_file}.")
    #endif

    return charge_map
#enddef


def period_files(root_dir):
    return [
        ("rga_fa18_inb", root_dir / "rga_fa18_inb.root"),
        ("rga_fa18_out", root_dir / "rga_fa18_out.root"),
        ("rga_sp19_inb", root_dir / "rga_sp19_inb.root"),
        ("rga_sp18_inb", root_dir / "rga_sp18_inb.root"),
        ("rga_sp18_out", root_dir / "rga_sp18_out.root"),
    ]
#enddef


def read_run_counts(period_label, root_path):
    if not root_path.is_file():
        raise RuntimeError(
            f"Missing ROOT file for {PERIOD_DISPLAY_NAMES[period_label]}: {root_path}\n"
            "Create the analysis-note ROOT files first. The expected filenames are "
            "rga_fa18_inb.root, rga_fa18_out.root, rga_sp18_inb.root, "
            "rga_sp18_out.root and rga_sp19_inb.root."
        )
    #endif

    try:
        with uproot.open(root_path) as root_file:
            if TREE_NAME not in root_file:
                raise RuntimeError(
                    f"Tree '{TREE_NAME}' is missing from {root_path}."
                )
            #endif

            tree = root_file[TREE_NAME]

            if "runnum" not in tree:
                raise RuntimeError(
                    f"Branch 'runnum' is missing from {root_path}:{TREE_NAME}."
                )
            #endif

            run_numbers = tree["runnum"].array(library="np")
    except RuntimeError:
        raise
    except Exception as exc:
        raise RuntimeError(f"Failed to read {root_path}: {exc}") from exc
    #endtry

    unique_runs, event_counts = np.unique(run_numbers, return_counts=True)
    return {int(run): int(count) for run, count in zip(unique_runs, event_counts)}
#enddef


def validate_inputs(all_run_counts, charge_map):
    problems = []

    for period_label, run_counts in all_run_counts.items():
        for run_number in sorted(run_counts):
            if run_number not in charge_map:
                problems.append(
                    f"{PERIOD_DISPLAY_NAMES[period_label]} run {run_number}: "
                    "no entry in Analysis_Note_charges.txt"
                )
                continue
            #endif

            charge_uC = charge_map[run_number]

            if not np.isfinite(charge_uC) or charge_uC <= 0.0:
                problems.append(
                    f"{PERIOD_DISPLAY_NAMES[period_label]} run {run_number}: "
                    f"charge is {charge_uC * 1000.0:g} nC; a positive finite charge is required"
                )
            #endif

            current_known, _ = resolve_current(period_label, run_number)

            if not current_known:
                problems.append(
                    f"{PERIOD_DISPLAY_NAMES[period_label]} run {run_number}: "
                    "beam current is not assigned in resolve_current()"
                )
            #endif
        #endfor
    #endfor

    if problems:
        preview = "\n".join(f"  - {problem}" for problem in problems)
        raise RuntimeError(
            "Input validation failed. No plots were made because every run in every "
            "ROOT tree must have both a positive charge and an explicit current assignment.\n"
            f"{preview}\n"
            "Add the missing charge record(s) to Analysis_Note_charges.txt or add the "
            "missing current assignment(s) to resolve_current(), then rerun."
        )
    #endif
#enddef


def compute_standardized_ylim(per_current_stats):
    y_min_raw = float("inf")
    y_max_raw = float("-inf")

    for stats in per_current_stats.values():
        values = stats["values"]
        errors = stats["errors"]

        if len(values) == 0:
            continue
        #endif

        y_min_raw = min(y_min_raw, float(np.min(values - errors)))
        y_max_raw = max(y_max_raw, float(np.max(values + errors)))
    #endfor

    if not np.isfinite(y_min_raw) or not np.isfinite(y_max_raw):
        return 0.0, 1.0
    #endif

    span = y_max_raw - y_min_raw
    padding = 0.06 * span if span > 0.0 else max(0.05 * abs(y_max_raw), 1.0)
    return y_min_raw - padding, y_max_raw + padding
#enddef


def add_trigger_annotations(ax, x_min, x_max, y_min, y_max):
    boundaries = set()

    for _, start_run, end_run in TRIGGER_RANGES_SP18:
        if x_min <= start_run <= x_max:
            boundaries.add(start_run)
        #endif
        if x_min <= end_run <= x_max:
            boundaries.add(end_run)
        #endif
    #endfor

    for boundary in sorted(boundaries):
        ax.axvline(boundary, color="0.35", linestyle="--", linewidth=1.0, zorder=1)
    #endfor

    y_text = y_max - 0.055 * (y_max - y_min)

    for label, start_run, end_run in TRIGGER_RANGES_SP18:
        visible_start = max(start_run, x_min)
        visible_end = min(end_run, x_max)

        if visible_end < visible_start:
            continue
        #endif

        ax.text(
            0.5 * (visible_start + visible_end),
            y_text,
            label,
            ha="center",
            va="center",
            fontsize=9,
            color="0.25",
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": "white",
                "alpha": 0.80,
                "edgecolor": "none",
            },
        )
    #endfor
#enddef


def build_period_statistics(period_label, run_counts, charge_map):
    grouped = defaultdict(lambda: {"runs": [], "values": [], "errors": [], "removed": []})
    removed_for_period = REMOVED_RUNS[period_label]

    for run_number, event_count in sorted(run_counts.items()):
        _, current_nA = resolve_current(period_label, run_number)
        charge_uC = charge_map[run_number]
        value = event_count / charge_uC
        error = math.sqrt(event_count) / charge_uC if event_count > 0 else 0.0

        grouped[current_nA]["runs"].append(run_number)
        grouped[current_nA]["values"].append(value)
        grouped[current_nA]["errors"].append(error)
        grouped[current_nA]["removed"].append(run_number in removed_for_period)
    #endfor

    statistics = {}

    for current_nA, raw in grouped.items():
        runs = np.asarray(raw["runs"], dtype=int)
        values = np.asarray(raw["values"], dtype=float)
        errors = np.asarray(raw["errors"], dtype=float)
        removed = np.asarray(raw["removed"], dtype=bool)
        accepted = ~removed

        if not np.any(accepted):
            raise RuntimeError(
                f"All {current_nA} nA runs in {PERIOD_DISPLAY_NAMES[period_label]} "
                "are marked as removed, so an accepted-run mean cannot be calculated."
            )
        #endif

        mean = float(np.mean(values[accepted]))
        spread = (
            float(np.std(values[accepted], ddof=1))
            if np.count_nonzero(accepted) > 1
            else 0.0
        )

        statistics[current_nA] = {
            "runs": runs,
            "values": values,
            "errors": errors,
            "removed": removed,
            "accepted": accepted,
            "mean": mean,
            "spread": spread,
        }
    #endfor

    return statistics
#enddef


def make_period_plot(period_label, per_current_stats, output_dir, trigger_version=False):
    fig, ax = plt.subplots(figsize=(12, 6.5))
    color_cycle = [f"C{index}" for index in range(10)]
    sorted_currents = sorted(per_current_stats)

    all_runs = np.concatenate([per_current_stats[current]["runs"] for current in sorted_currents])
    x_min = float(np.min(all_runs))
    x_max = float(np.max(all_runs))

    for index, current_nA in enumerate(sorted_currents):
        stats = per_current_stats[current_nA]
        color = color_cycle[index % len(color_cycle)]
        accepted = stats["accepted"]
        removed = stats["removed"]

        if np.any(accepted):
            ax.errorbar(
                stats["runs"][accepted],
                stats["values"][accepted],
                yerr=stats["errors"][accepted],
                fmt="o",
                markersize=4.5,
                linestyle="none",
                markerfacecolor=color,
                markeredgecolor=color,
                ecolor=color,
                elinewidth=0.9,
                capsize=0,
                label=f"{current_nA} nA",
                zorder=3,
            )
        #endif

        if np.any(removed):
            ax.errorbar(
                stats["runs"][removed],
                stats["values"][removed],
                yerr=stats["errors"][removed],
                fmt="x",
                markersize=7.0,
                markeredgewidth=1.7,
                linestyle="none",
                color=color,
                ecolor=color,
                elinewidth=0.9,
                capsize=0,
                zorder=4,
            )
        #endif

        ax.axhline(stats["mean"], color=color, linestyle="--", linewidth=1.2, zorder=2)
    #endfor

    y_min, y_max = compute_standardized_ylim(per_current_stats)
    ax.set_ylim(y_min, y_max)

    if trigger_version:
        add_trigger_annotations(ax, x_min, x_max, y_min, y_max)
    #endif

    display_name = PERIOD_DISPLAY_NAMES[period_label]
    ax.set_xlabel("Run number", fontsize=13)
    ax.set_ylabel(r"DVCS events / $\mu$C", fontsize=13)
    ax.set_title(f"{display_name}", fontsize=14)
    ax.grid(True, alpha=0.25)
    ax.tick_params(axis="both", labelsize=11)

    handles, labels = ax.get_legend_handles_labels()

    if any(np.any(stats["removed"]) for stats in per_current_stats.values()):
        from matplotlib.lines import Line2D
        handles.append(
            Line2D(
                [0], [0], marker="x", linestyle="none", color="0.2",
                markersize=7, markeredgewidth=1.7,
            )
        )
        labels.append("Removed run")
    #endif

    ax.legend(handles, labels, title="Beam current", fontsize=10, title_fontsize=10)
    fig.tight_layout()

    suffix = "_trigger" if trigger_version else ""
    output_path = output_dir / f"dvcs_per_uC_{period_label}{suffix}.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}")
#enddef


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        charge_map = load_charge_map(args.charge_file)
        all_run_counts = {}

        for period_label, root_path in period_files(args.root_dir):
            all_run_counts[period_label] = read_run_counts(period_label, root_path)
        #endfor

        validate_inputs(all_run_counts, charge_map)

        print("=" * 88)
        print("DVCS per-run analysis-note plots")
        print(f"Charge file: {args.charge_file}")
        print("Charge conversion: nC / 1000 = microC")
        print(f"ROOT directory: {args.root_dir}")
        print(f"Output directory: {args.output_dir}")
        print("=" * 88)

        for period_label, _ in period_files(args.root_dir):
            statistics = build_period_statistics(
                period_label,
                all_run_counts[period_label],
                charge_map,
            )

            print(f"\n{PERIOD_DISPLAY_NAMES[period_label]}")

            for current_nA in sorted(statistics):
                stats = statistics[current_nA]
                n_accepted = int(np.count_nonzero(stats["accepted"]))
                n_removed = int(np.count_nonzero(stats["removed"]))
                print(
                    f"  {current_nA:>2} nA: mean = {stats['mean']:.6f} events/microC, "
                    f"spread = {stats['spread']:.6f}, accepted = {n_accepted}, "
                    f"removed = {n_removed}"
                )
            #endfor

            present_removed = sorted(
                run
                for run in REMOVED_RUNS[period_label]
                if run in all_run_counts[period_label]
            )
            absent_removed = sorted(REMOVED_RUNS[period_label] - set(all_run_counts[period_label]))

            if present_removed:
                print("  Removed runs present in ROOT tree: " + ", ".join(map(str, present_removed)))
            #endif

            if absent_removed:
                print(
                    "  Listed removed runs absent from this ROOT tree: "
                    + ", ".join(map(str, absent_removed))
                )
            #endif

            make_period_plot(period_label, statistics, args.output_dir, trigger_version=False)

            if period_label in {"rga_sp18_inb", "rga_sp18_out"}:
                make_period_plot(period_label, statistics, args.output_dir, trigger_version=True)
            #endif
        #endfor

    except RuntimeError as exc:
        print(f"\nFATAL ERROR:\n{exc}", file=sys.stderr)
        return 1
    #endtry

    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif
