#!/usr/bin/env python3
"""
main.py

Main script to load ROOT trees and run the temporary dilution-factor routine,
ALL, and Pt calculations.
"""

import data_loader
import calculate_dilution_factor as cdf
import plot_normalized_yields as pny


def main():
    """
    Entry point: define xB bins, load trees, compute dilution factors.
    """
    # Define xB bin edges for analysis
    xB_bins = [0.0, 0.14, 0.24, 0.34, 0.44, 0.54, 0.64, 1.00]

    # Load PhysicsEvents trees for all periods
    trees = data_loader.load_root_trees()

    # plot normalized yields for NH3, C, CH2, He, ET
    pny.plot_normalized_yields(trees, xB_bins)

    # Use the temporary manual Df routine
    cdf.calculate_dilution_factor_temp(trees, xB_bins)

    # TODO: compute ALL(xB_bins) and solve for Pt using Df from CSV
    # end TODO


if __name__ == "__main__":
    main()