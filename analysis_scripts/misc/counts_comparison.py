#!/usr/bin/env python3

import ROOT
import sys

def count_events(filename):
    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        print("Error: could not open file:", filename)
        sys.exit(1)
    #endif

    tree = f.Get("PhysicsEvents")
    if not tree:
        print("Error: tree 'PhysicsEvents' not found in", filename)
        sys.exit(1)
    #endif

    total = tree.GetEntries()
    passed = 0

    for entry in tree:
        if entry.Mx2 < 1.07:
            passed += 1
        #endif
    #endfor

    fraction = passed / total if total > 0 else 0

    return total, passed, fraction
#enddef


if len(sys.argv) != 3:
    print("Usage: python count_mx2_events.py file1.root file2.root")
    sys.exit(1)
#endif

file1 = sys.argv[1]
file2 = sys.argv[2]

total1, pass1, frac1 = count_events(file1)
total2, pass2, frac2 = count_events(file2)

ratio = pass2 / pass1 if pass1 > 0 else 0

print("--------------------------------------------------")
print("File 1:", file1)
print("  total events:", total1)
print("  Mx2 < 1.07:", pass1)
print("  fraction:", frac1)
print("")
print("File 2:", file2)
print("  total events:", total2)
print("  Mx2 < 1.07:", pass2)
print("  fraction:", frac2)
print("")
print("Ratio (file2 / file1) using counts:", ratio)
print("--------------------------------------------------")