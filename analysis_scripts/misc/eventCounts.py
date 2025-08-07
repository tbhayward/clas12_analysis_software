#!/usr/bin/env python3
import sys
import uproot
import numpy as np

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file.root>")
        sys.exit(1)

    infile = sys.argv[1]
    # Open the file and get the tree
    with uproot.open(infile) as f:
        tree = f["PhysicsEvents"]
        # Load only the needed branches
        arr = tree.arrays(["t1", "pTmiss", "open_angle_ep2"], library="np")

    t1         = arr["t1"]
    pTmiss     = arr["pTmiss"]
    open_ang   = arr["open_angle_ep2"]

    # Apply cuts: -t1 < 1  ↔  (-t1) < 1
    cut1 = (-t1 < 1)
    #            pTmiss < 0.2
    cut2 = (pTmiss < 0.2)
    #            open_angle_ep2 > 5
    cut3 = (open_ang > 5)

    # Combined mask
    mask = cut1 & cut2 & cut3

    n_selected = np.count_nonzero(mask)
    print(f"Total events passing cuts: {n_selected}")

if __name__ == "__main__":
    main()