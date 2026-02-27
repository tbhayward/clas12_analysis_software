#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ft_photon_misid_mc_matching_pid_summary_all_mc.py

Purpose:
  Diagnose potential FT "photon" (reco pid=22) mis-ID as electrons by looking at
  MC truth matching PID for forward-angle candidates.

Selection (per file / per period):
  - reconstructed particle_pid == 22
  - theta < THETA_MAX_DEG   (branch: theta, assumed degrees)

Summaries:
  - Always report:
      mc_matching_pid == 22
      mc_matching_pid == 11
      other (everything else)
  - NEW (requested):
      ALSO break out any OTHER individual mc_matching_pid values whose fraction
      exceeds a threshold (default 1%) within that period, and also for the combined total.

Parallel processing:
  - All five MC files processed in parallel
  - HARD LIMIT: at most 5 workers

Run (tcsh, single line):
  python3 ft_photon_misid_mc_matching_pid_summary_all_mc.py --theta-max 6.0 --threshold 0.01

Dependencies:
  pip install uproot numpy
"""

import os
import sys
import argparse
import traceback
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import uproot
from concurrent.futures import ProcessPoolExecutor, as_completed


# -----------------------------------------------------------------------------
# USER SETTINGS (STATIC)
# -----------------------------------------------------------------------------

TREE_NAME = "PhysicsEvents"
MAX_WORKERS = 5

MC_FILES = {
    "Sp18 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_sp18_inb.root",
    "Sp18 Out": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_sp18_out.root",
    # "Fa18 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_fa18_inb.root",
    "Fa18 Inb": "/scratch/thayward/mc_truth_matching_test1.root",
    "Fa18 Out": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_fa18_out.root",
    "Sp19 Inb": "/volatile/clas12/thayward/dvcs_temp/calibration_dvcsgen_rga_sp19_inb.root",
}


# -----------------------------------------------------------------------------
# INTERNALS
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class PidSummary:
    n_total: int
    n_22: int
    n_11: int
    n_other: int
    # full distribution of mc_matching_pid counts for the selected sample
    pid_counts: Dict[int, int]


def fatal(msg: str) -> None:
    sys.stderr.write(f"\nFATAL: {msg}\n")
    sys.exit(1)
    #endif


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Summarize mc_matching_pid for reconstructed photons with theta < threshold, over all 5 MC files."
    )
    p.add_argument(
        "--theta-max",
        type=float,
        default=6.0,
        help="Max theta in degrees for selection (default: 6.0).",
    )
    p.add_argument(
        "--tree",
        type=str,
        default=TREE_NAME,
        help="TTree name (default: PhysicsEvents).",
    )
    p.add_argument(
        "--max-workers",
        type=int,
        default=5,
        help="Max parallel workers (hard limit 5).",
    )
    p.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help="Fraction threshold for listing additional mc_matching_pid categories (default: 0.01 = 1%%).",
    )
    return p.parse_args()
    #endif


def validate_inputs(theta_max: float, tree_name: str, max_workers: int, threshold: float) -> None:
    if theta_max <= 0.0:
        fatal(f"--theta-max must be > 0. Got {theta_max}")
    #endif
    if max_workers <= 0 or max_workers > 5:
        fatal(f"--max-workers must be in [1,5]. Got {max_workers}")
    #endif
    if not tree_name:
        fatal("--tree must be a non-empty string.")
    #endif
    if threshold <= 0.0 or threshold >= 1.0:
        fatal(f"--threshold must be in (0,1). Got {threshold}")
    #endif

    for label, path in MC_FILES.items():
        if not os.path.isfile(path):
            fatal(f"Missing MC file for '{label}': {path}")
        #endif
    #endfor
    #enddef


def summarize_one_file(args: Tuple[str, str, str, float]) -> Tuple[str, PidSummary]:
    """
    Worker: open one MC ROOT file, apply selection, summarize mc_matching_pid.
    Returns (period_label, PidSummary).
    """
    period_label, root_path, tree_name, theta_max = args

    try:
        with uproot.open(root_path) as f:
            if tree_name not in f:
                raise KeyError(f"TTree '{tree_name}' not found. Keys: {list(f.keys())}")
            #endif

            tree = f[tree_name]
            required = ["particle_pid", "theta", "mc_matching_pid"]
            for br in required:
                if br not in tree.keys():
                    raise KeyError(
                        f"Required branch '{br}' not found in tree '{tree_name}'. "
                        f"Available branches (first 80): {list(tree.keys())[:80]}"
                    )
                #endif
            #endfor

            arr = tree.arrays(required, library="np")
            pid = arr["particle_pid"]
            theta = arr["theta"]
            mc_pid = arr["mc_matching_pid"]

            mask = np.isfinite(theta) & np.isfinite(pid) & np.isfinite(mc_pid)
            pid = pid[mask].astype(np.int64)
            theta = theta[mask].astype(np.float64)
            mc_pid = mc_pid[mask].astype(np.int64)

            sel = (pid == 22) & (theta < theta_max)
            mc_sel = mc_pid[sel]

            n_total = int(mc_sel.size)

            if n_total > 0:
                uniq, cnt = np.unique(mc_sel, return_counts=True)
                pid_counts = {int(u): int(c) for (u, c) in zip(uniq, cnt)}
            else:
                pid_counts = {}
            #endif

            n_22 = int(pid_counts.get(22, 0))
            n_11 = int(pid_counts.get(11, 0))
            n_other = int(n_total - n_22 - n_11)

            summary = PidSummary(
                n_total=n_total,
                n_22=n_22,
                n_11=n_11,
                n_other=n_other,
                pid_counts=pid_counts,
            )
            return period_label, summary

    except Exception as e:
        tb = traceback.format_exc()
        raise RuntimeError(
            f"Failed processing MC file:\n"
            f"  period = {period_label}\n"
            f"  path   = {root_path}\n"
            f"  tree   = {tree_name}\n"
            f"  error  = {repr(e)}\n\n{tb}"
        )
    #enddef


def frac(n: int, d: int) -> float:
    if d <= 0:
        return 0.0
    #endif
    return float(n) / float(d)
    #enddef


def print_extra_categories(pid_counts: Dict[int, int], n_total: int, threshold: float) -> None:
    """
    Print any mc_matching_pid categories (excluding 22 and 11) whose fraction >= threshold.
    """
    if n_total <= 0:
        return
    #endif

    extras = []
    for pid_val, n in pid_counts.items():
        if pid_val == 22 or pid_val == 11:
            continue
        #endif
        f = float(n) / float(n_total)
        if f >= threshold:
            extras.append((pid_val, n, f))
        #endif
    #endfor

    extras.sort(key=lambda t: (-t[2], -t[1], t[0]))

    if len(extras) == 0:
        return
    #endif

    print(f"  Additional mc_matching_pid categories with frac >= {100.0*threshold:.2f}%:")
    for pid_val, n, f in extras:
        print(f"    pid == {pid_val}: {n}   frac = {100.0*f:.4f}%")
    #endfor
    #enddef


def merge_pid_counts(dst: Dict[int, int], src: Dict[int, int]) -> None:
    for k, v in src.items():
        dst[k] = int(dst.get(k, 0) + int(v))
    #endfor
    #enddef


def main() -> None:
    args = parse_args()
    theta_max = float(args.theta_max)
    tree_name = str(args.tree)
    max_workers = int(args.max_workers)
    threshold = float(args.threshold)

    validate_inputs(theta_max, tree_name, max_workers, threshold)

    items = sorted(MC_FILES.items(), key=lambda kv: kv[0])
    tasks = [(label, path, tree_name, theta_max) for (label, path) in items]

    n_workers = min(max_workers, len(tasks), MAX_WORKERS)
    if n_workers < 1:
        fatal("No tasks to run.")
    #endif

    results: Dict[str, PidSummary] = {}

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(summarize_one_file, t) for t in tasks]
        for fut in as_completed(futures):
            period_label, summary = fut.result()
            if period_label in results:
                fatal(f"Duplicate result for period '{period_label}'")
            #endif
            results[period_label] = summary
        #endfor
    #endwith

    for period_label in MC_FILES.keys():
        if period_label not in results:
            fatal(f"Missing result for period '{period_label}'")
        #endif
    #endfor

    periods_ordered = ["Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"]

    print("")
    print("------------------------------------------------------------")
    print("FT photon mis-ID diagnostic using mc_matching_pid (MC only)")
    print("Selection:")
    print("  reconstructed particle_pid == 22")
    print(f"  theta < {theta_max:.3f} (deg)")
    print("------------------------------------------------------------")
    print("Per-period summary:")
    print("")

    tot_N = 0
    tot_22 = 0
    tot_11 = 0
    tot_other = 0
    tot_pid_counts: Dict[int, int] = {}

    for period_label in periods_ordered:
        s = results[period_label]
        merge_pid_counts(tot_pid_counts, s.pid_counts)

        tot_N += s.n_total
        tot_22 += s.n_22
        tot_11 += s.n_11
        tot_other += s.n_other

        print(f"{period_label}:")
        print(f"  Total selected: {s.n_total}")
        print(f"    mc_matching_pid == 22: {s.n_22}   frac = {100.0*frac(s.n_22, s.n_total):.4f}%")
        print(f"    mc_matching_pid == 11: {s.n_11}   frac = {100.0*frac(s.n_11, s.n_total):.4f}%")
        print(f"    other:                 {s.n_other}   frac = {100.0*frac(s.n_other, s.n_total):.4f}%")
        print_extra_categories(s.pid_counts, s.n_total, threshold)
        print("")
    #endfor

    print("------------------------------------------------------------")
    print("Combined (all 5 MC files):")
    print(f"  Total selected: {tot_N}")
    print(f"    mc_matching_pid == 22: {tot_22}   frac = {100.0*frac(tot_22, tot_N):.4f}%")
    print(f"    mc_matching_pid == 11: {tot_11}   frac = {100.0*frac(tot_11, tot_N):.4f}%")
    print(f"    other:                 {tot_other}   frac = {100.0*frac(tot_other, tot_N):.4f}%")
    print_extra_categories(tot_pid_counts, tot_N, threshold)
    print("------------------------------------------------------------")
    print("")
#enddef


if __name__ == "__main__":
    main()
#endif