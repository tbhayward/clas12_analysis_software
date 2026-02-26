#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ft_photon_misid_mc_matching_pid_summary.py

Purpose:
  Scan a ROOT TTree and diagnose potential FT "photon" (pid=22) mis-ID as electrons
  by looking at the MC truth-matching PID.

Selection:
  - reconstructed particle_pid == 22   (photons)
  - theta < 6.0 degrees               (branch: theta)

For the selected candidates, summarize mc_matching_pid into:
  - true photons: mc_matching_pid == 22
  - true electrons: mc_matching_pid == 11
  - other: everything else (including 0 / -9999 / etc)

Print:
  - total selected N
  - counts for 22 / 11 / other
  - fractions relative to total selected

Usage (tcsh):
  python3 ft_photon_misid_mc_matching_pid_summary.py --input /path/to/file.root --tree PhysicsEvents --theta-max 6.0

Dependencies:
  pip install uproot numpy
"""

import argparse
import os
import sys

import numpy as np
import uproot


def fatal(msg: str) -> None:
    sys.stderr.write(f"\nFATAL: {msg}\n")
    sys.exit(1)
    #endif


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Summarize mc_matching_pid for reconstructed photons with theta < threshold."
    )
    p.add_argument("--input", required=True, help="Input ROOT file path.")
    p.add_argument("--tree", default="PhysicsEvents", help="TTree name (default: PhysicsEvents).")
    p.add_argument("--theta-max", type=float, default=6.0, help="Max theta in degrees (default: 6.0).")
    return p.parse_args()
    #endif


def main() -> None:
    args = parse_args()

    inpath = args.input
    tname = args.tree
    theta_max = float(args.theta_max)

    if not os.path.isfile(inpath):
        fatal(f"Input file not found: {inpath}")
    #endif
    if theta_max <= 0.0:
        fatal(f"--theta-max must be > 0. Got {theta_max}")
    #endif

    with uproot.open(inpath) as f:
        if tname not in f:
            fatal(f"TTree '{tname}' not found in file. Available keys: {list(f.keys())}")
        #endif

        tree = f[tname]
        required = ["particle_pid", "theta", "mc_matching_pid"]
        for br in required:
            if br not in tree.keys():
                fatal(
                    f"Required branch '{br}' not found in TTree '{tname}'. "
                    f"Available branches (first 80): {list(tree.keys())[:80]}"
                )
            #endif
        #endfor

        arr = tree.arrays(required, library="np")

        pid = arr["particle_pid"]
        theta = arr["theta"]
        mc_pid = arr["mc_matching_pid"]

        # Basic finiteness checks (keep deterministic)
        mask = np.isfinite(theta) & np.isfinite(pid) & np.isfinite(mc_pid)
        pid = pid[mask].astype(np.int64)
        theta = theta[mask].astype(np.float64)
        mc_pid = mc_pid[mask].astype(np.int64)

        # Reco photons with theta < threshold
        sel = (pid == 22) & (theta < theta_max)
        mc_sel = mc_pid[sel]

        n_tot = int(mc_sel.size)
        if n_tot == 0:
            print("")
            print("No entries matched selection:")
            print(f"  particle_pid == 22 and theta < {theta_max:.3f} (deg)")
            print("")
            return
        #endif

        n_true_photon = int(np.count_nonzero(mc_sel == 22))
        n_true_electron = int(np.count_nonzero(mc_sel == 11))
        n_other = int(n_tot - n_true_photon - n_true_electron)

        f_true_photon = float(n_true_photon) / float(n_tot)
        f_true_electron = float(n_true_electron) / float(n_tot)
        f_other = float(n_other) / float(n_tot)

        print("")
        print("------------------------------------------------------------")
        print("FT photon mis-ID diagnostic using mc_matching_pid")
        print("Selection:")
        print(f"  reconstructed particle_pid == 22")
        print(f"  theta < {theta_max:.3f} (deg)")
        print("------------------------------------------------------------")
        print(f"Total selected: {n_tot}")
        print("")
        print("mc_matching_pid breakdown:")
        print(f"  pid == 22 (true photon):   {n_true_photon}   fraction = {100.0*f_true_photon:.4f}%")
        print(f"  pid == 11 (true electron): {n_true_electron}   fraction = {100.0*f_true_electron:.4f}%")
        print(f"  other:                     {n_other}   fraction = {100.0*f_other:.4f}%")
        print("------------------------------------------------------------")
        print("")
    #endwith
#enddef


if __name__ == "__main__":
    main()
#endif