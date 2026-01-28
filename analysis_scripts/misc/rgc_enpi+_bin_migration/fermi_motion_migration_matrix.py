#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mx2_mc_mx2_fraction_table.py

Read a ROOT file at runtime, load the PhysicsEvents TTree, select events with:
  0.81 <= Mx2 <= 1.00

Then compute the percentage distribution of mc_Mx2 values in the bins:
  [-10, 0.81]
  [0.81, 1.00]
  [1.00, 1.15]
  [1.15, 1.30]
  [1.30, 1.45]
  [1.45, 1.60]
  [1.60, 20]

Prints a table to the terminal.

Usage (tcsh):
  python3 mx2_mc_mx2_fraction_table.py /path/to/file.root
"""

import sys
import math

import ROOT


def die(msg):
    raise SystemExit(f"FATAL: {msg}")


def format_pct(x):
    # x is fraction in [0,1]
    return f"{100.0 * x:8.3f}%"


def main():
    if len(sys.argv) != 2:
        die("Expected exactly 1 argument: /path/to/input.root")

    inpath = sys.argv[1]

    ROOT.gROOT.SetBatch(True)

    f = ROOT.TFile.Open(inpath, "READ")
    if not f or f.IsZombie():
        die(f"Could not open input ROOT file: {inpath}")

    t = f.Get("PhysicsEvents")
    if not t:
        die("Tree 'PhysicsEvents' not found in file.")

    # Check branches exist (fail fast).
    if not (t.GetBranch("Mx2") or t.GetLeaf("Mx2")):
        die("Branch/leaf 'Mx2' not found in PhysicsEvents.")
    if not (t.GetBranch("mc_Mx2") or t.GetLeaf("mc_Mx2")):
        die("Branch/leaf 'mc_Mx2' not found in PhysicsEvents.")

    # Define selection and mc_Mx2 bin edges.
    mx2_sel_min = 0.81
    mx2_sel_max = 1.00

    # Bin edges in ascending order.
    # We will treat them as half-open [lo, hi) except the last bin which is [lo, hi].
    edges = [-10.0, 0.81, 1.00, 1.15, 1.30, 1.45, 1.60, 20.0]

    # Counters.
    total_selected = 0
    counts = [0] * (len(edges) - 1)

    # Loop.
    n_entries = int(t.GetEntries())
    for i in range(n_entries):
        t.GetEntry(i)

        mx2 = float(getattr(t, "Mx2"))
        if mx2 < mx2_sel_min or mx2 > mx2_sel_max:
            continue

        mc_mx2 = float(getattr(t, "mc_Mx2"))
        total_selected += 1

        # Fill mc_Mx2 into bins:
        # bins k correspond to [edges[k], edges[k+1]) except last includes upper edge.
        placed = False
        for k in range(len(counts)):
            lo = edges[k]
            hi = edges[k + 1]

            if k < len(counts) - 1:
                if (mc_mx2 >= lo) and (mc_mx2 < hi):
                    counts[k] += 1
                    placed = True
                    break
                #endif
            else:
                if (mc_mx2 >= lo) and (mc_mx2 <= hi):
                    counts[k] += 1
                    placed = True
                    break
                #endif
            #endif
        #endfor

        # If mc_Mx2 is outside the requested binning, we ignore it but warn at end.
        # (Keeps behavior explicit/deterministic.)
        #endif
    #endfor

    if total_selected == 0:
        die(f"No entries passed the selection {mx2_sel_min} <= Mx2 <= {mx2_sel_max}.")

    # Print table.
    print("")
    print("------------------------------------------------------------")
    print(f"Input file: {inpath}")
    print("Tree: PhysicsEvents")
    print(f"Selection: {mx2_sel_min:.2f} <= Mx2 <= {mx2_sel_max:.2f}")
    print(f"Selected entries: {total_selected}")
    print("------------------------------------------------------------")
    print("")
    print(f"{'mc_Mx2 bin':>20}   {'count':>12}   {'percent':>12}")
    print(f"{'-'*20}   {'-'*12}   {'-'*12}")

    for k in range(len(counts)):
        lo = edges[k]
        hi = edges[k + 1]
        label = f"[{lo:.2f}, {hi:.2f}]"
        frac = counts[k] / float(total_selected)
        print(f"{label:>20}   {counts[k]:12d}   {format_pct(frac):>12}")
    #endfor

    # Basic consistency check.
    sum_counts = sum(counts)
    if sum_counts != total_selected:
        missing = total_selected - sum_counts
        print("")
        print(f"WARNING: {missing} selected events had mc_Mx2 outside [{edges[0]}, {edges[-1]}] and were not counted.")
    #endif

    print("")
    print("Done.")


if __name__ == "__main__":
    main()
#endif