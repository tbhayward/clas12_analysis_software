#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mx2_mc_mx2_fraction_table.py

Read a ROOT file at runtime, load the PhysicsEvents TTree, select events with:
  0.81 <= Mx2   <= 1.00
  0.10 <= xB    <= 0.60
  -1.0 <= tprime <= 0.0

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
import ROOT


def die(msg):
    raise SystemExit(f"FATAL: {msg}")


def format_pct(x):
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

    # Fail-fast branch checks.
    required = ["Mx2", "mc_Mx2", "x", "tprime"]
    missing = []
    for br in required:
        if not (t.GetBranch(br) or t.GetLeaf(br)):
            missing.append(br)
        #endif
    #endfor
    if len(missing) > 0:
        die(f"Missing required branch/leaf(s) in PhysicsEvents: {', '.join(missing)}")

    # Selections
    mx2_sel_min = 0.81
    mx2_sel_max = 1.00

    xb_sel_min = 0.10
    xb_sel_max = 0.60

    tprime_sel_min = -1.0
    tprime_sel_max = 0.0

    # mc_Mx2 bin edges (ascending).
    edges = [-10.0, 0.81, 1.00, 1.15, 1.30, 1.45, 1.60, 20.0]
    counts = [0] * (len(edges) - 1)

    total_selected = 0
    total_outside_mc_range = 0

    n_entries = int(t.GetEntries())
    for i in range(n_entries):
        t.GetEntry(i)

        mx2 = float(getattr(t, "Mx2"))
        if mx2 < mx2_sel_min or mx2 > mx2_sel_max:
            continue
        #endif

        xb = float(getattr(t, "x"))
        if xb < xb_sel_min or xb > xb_sel_max:
            continue
        #endif

        tp = float(getattr(t, "tprime"))
        if tp < tprime_sel_min or tp > tprime_sel_max:
            continue
        #endif

        mc_mx2 = float(getattr(t, "mc_Mx2"))

        total_selected += 1

        placed = False
        for k in range(len(counts)):
            lo = edges[k]
            hi = edges[k + 1]

            # Half-open bins [lo, hi), except last bin includes upper edge.
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

        if not placed:
            total_outside_mc_range += 1
        #endif
    #endfor

    if total_selected == 0:
        die(
            "No entries passed the selection: "
            f"{mx2_sel_min:.2f} <= Mx2 <= {mx2_sel_max:.2f}, "
            f"{xb_sel_min:.2f} <= xB <= {xb_sel_max:.2f}, "
            f"{tprime_sel_min:.2f} <= tprime <= {tprime_sel_max:.2f}."
        )

    # Print table.
    print("")
    print("------------------------------------------------------------")
    print(f"Input file: {inpath}")
    print("Tree: PhysicsEvents")
    print("Selection:")
    print(f"  {mx2_sel_min:.2f} <= Mx2    <= {mx2_sel_max:.2f}")
    print(f"  {xb_sel_min:.2f} <= xB     <= {xb_sel_max:.2f}")
    print(f"  {tprime_sel_min:.2f} <= tprime <= {tprime_sel_max:.2f}")
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

    if total_outside_mc_range > 0:
        print("")
        print(
            f"WARNING: {total_outside_mc_range} selected events had mc_Mx2 outside "
            f"[{edges[0]}, {edges[-1]}] and were not counted."
        )
    #endif

    print("")
    print("Done.")


if __name__ == "__main__":
    main()
#endif