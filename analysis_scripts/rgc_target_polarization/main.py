#!/usr/bin/env python3
"""
main.py

Main script to load ROOT trees and compute dilution factors, ALL, and Pt.
"""

import data_loader
import calculate_dilution_factor


def main():
    """
    Entry point: define xB bins, load trees, compute dilution factors.
    """
    # Define xB bin edges for analysis
    xB_bins = [0.0, 0.14, 0.24, 0.34, 0.44, 0.54, 0.64, 1.00]

    # Load PhysicsEvents trees for all periods
    trees = data_loader.load_root_trees()

    # Calculate and save dilution-factor CSV + PDF
    calculate_dilution_factor.calculate_and_save(trees, xB_bins)

    # TODO: compute ALL(xB_bins) and solve for Pt using Df from CSV
    # end TODO
#end main

if __name__ == "__main__":
    main()
#endif
