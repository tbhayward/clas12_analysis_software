#!/usr/bin/env python3
import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <base_file.root> <recalc_file.root>")
        sys.exit(1)

    base_file, recalc_file = sys.argv[1], sys.argv[2]

    # Open the ROOT files and extract the 't' branch
    with uproot.open(base_file) as f1:
        t_base = f1["PhysicsEvents"].arrays("t", library="np")["t"]
    with uproot.open(recalc_file) as f2:
        t_recalc = f2["PhysicsEvents"].arrays("t", library="np")["t"]

    # Define histogram bins from -10 to 1
    bins = np.linspace(-10, 1, 111)  # 0.1-wide bins

    # Plot
    plt.figure(figsize=(8, 6))
    plt.hist(t_base,   bins=bins, histtype='step', label='base', color='blue')
    plt.hist(t_recalc, bins=bins, histtype='step', label='recalculated', color='red')

    plt.xlabel(r'$t$')
    plt.ylabel('Counts')
    plt.legend(loc='upper right')
    plt.tight_layout()

    # Save
    outpath = "/u/home/thayward/t_example.pdf"
    plt.savefig(outpath)
    print(f"Saved histogram to {outpath}")

if __name__ == "__main__":
    main()