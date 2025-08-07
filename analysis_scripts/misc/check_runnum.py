#!/usr/bin/env python3
import sys
import uproot
import numpy as np

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
                # skip lines that don't have an integer in column 0
                pass
    return runnums

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file.root>")
        sys.exit(1)

    rootfile = sys.argv[1]
    # load allowed runnums from CSV
    csv_runnums = load_csv_runnums(CSV_PATH)

    # open the ROOT file and grab the runnum branch
    with uproot.open(rootfile) as f:
        tree = f["PhysicsEvents"]
        runnums = tree.array("runnum", library="np")

    unique_runnums = np.unique(runnums)
    missing = [int(r) for r in unique_runnums if int(r) not in csv_runnums]

    if missing:
        print("The following runnum values appear in the ROOT tree but are NOT in the CSV:")
        for r in missing:
            print(r)
        sys.exit(1)
    else:
        print("All runnum entries in the ROOT file are present in the CSV.")
        sys.exit(0)

if __name__ == "__main__":
    main()