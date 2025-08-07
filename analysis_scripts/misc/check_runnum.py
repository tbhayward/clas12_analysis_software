#!/usr/bin/env python3
import sys
import uproot
import numpy as np
from collections import Counter

CSV_PATH = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

def load_csv_runnums(path):
    runnums = set()
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            try:
                runnums.add(int(parts[0]))
            except ValueError:
                pass
    return runnums

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file.root>")
        sys.exit(1)

    rootfile = sys.argv[1]
    csv_runnums = load_csv_runnums(CSV_PATH)

    with uproot.open(rootfile) as f:
        tree = f["PhysicsEvents"]
        data = tree.arrays(["runnum"], library="np")
        runnums = data["runnum"]

    # count how many entries of each runnum
    counts = Counter(int(r) for r in runnums)
    missing = [r for r in counts if r not in csv_runnums]

    if missing:
        print("Run numbers in ROOT tree not found in CSV:")
        for r in sorted(missing):
            print(f"  Run {r}: {counts[r]} entries")
        sys.exit(1)
    else:
        print("All runnum entries in the ROOT file are present in the CSV.")
        sys.exit(0)

if __name__ == "__main__":
    main()