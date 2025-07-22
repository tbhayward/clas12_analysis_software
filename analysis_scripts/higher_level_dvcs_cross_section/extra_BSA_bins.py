#!/usr/bin/env python3
"""
extract_BSA_bins.py

Usage:
    python extract_BSA_bins.py imports/rga_prl_bsa.txt

Reads the BSA file, splits it into bins by φ (when φ decreases you start a new bin),
and prints how many bins were found plus ⟨Q²⟩, ⟨x_B⟩, ⟨t⟩ for each.
"""

import sys
import numpy as np

def main():
    if len(sys.argv) != 2:
        print("Usage: python extract_BSA_bins.py <path/to/rga_prl_bsa.txt>")
        sys.exit(1)

    datafile = sys.argv[1]

    bins = []  # will hold one dict per φ-bin
    current = {'phi': [], 'Q2': [], 'xB': [], 't': []}
    prev_phi = -1.0

    with open(datafile) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            phi, Q2, xB, t, Eb, A, sigA = map(float, line.split())

            # if φ dropped compared to previous, start a new bin
            if current['phi'] and phi < prev_phi:
                bins.append(current)
                current = {'phi': [], 'Q2': [], 'xB': [], 't': []}

            current['phi'].append(phi)
            current['Q2'].append(Q2)
            current['xB'].append(xB)
            current['t'].append(t)
            prev_phi = phi

    # don't forget the last bin
    if current['phi']:
        bins.append(current)

    print(f"Found {len(bins)} φ-bins in {datafile}\n")
    for i, b in enumerate(bins, start=1):
        mean_Q2 = np.mean(b['Q2'])
        mean_xB = np.mean(b['xB'])
        mean_t  = np.mean(b['t'])
        print(f"Bin {i:2d}: ⟨Q²⟩ = {mean_Q2:.3f} GeV², "
              f"⟨x_B⟩ = {mean_xB:.4f}, ⟨t⟩ = {mean_t:.3f} GeV²")

if __name__ == '__main__':
    main()