#ifndef CROSS_SECTION_RUN_PERIOD_CONSISTENCY_H
#define CROSS_SECTION_RUN_PERIOD_CONSISTENCY_H

#include <string>

// -----------------------------------------------------------------------------
// Cross-check of Hayward pass-2 normed cross sections across compatible
// run-period sets.
//
// Available comparison sets are selected by helicity:
//   unpol: Fa18 Inb, Fa18 Out, Sp19 Inb, Sp18 Inb, Sp18 Out
//   pos:   Fa18 Inb, Fa18 Out, Sp19 Inb
//   neg:   Fa18 Inb, Fa18 Out, Sp19 Inb
//
// Sp18 Inb/Out are intentionally excluded from pos/neg comparisons because
// their helicity-resolved accumulated charge is not available. No separate
// 10.2 GeV label is used; Sp19 Inb is the 10.2 GeV dataset.
//
// For each selected helicity state, this module:
//   1) computes a per-cell chi2 and reduced chi2 with respect to the simple
//      arithmetic mean across the available run periods at each phi bin,
//   2) overlays the selected period-dependent phi distributions as ratios to
//      that mean in each (xB, Q2, -t) cell,
//   3) writes a detailed text summary of those chi2 values,
//   4) produces a histogram of the reduced-chi2 distribution and fits it with
//      a reduced-chi2 probability-density function.
//
// Input CSV:
//   output/csvs/dvcs_pass2_analysis.csv
//
// Output directory structure under output_base_dir:
//   unpol/
//     overlays/period_consistency_xB_<ix>.png
//     reduced_chi2_distribution.png
//     reduced_chi2_summary.txt
//   pos/
//     ...
//   neg/
//     ...
// -----------------------------------------------------------------------------

void plot_cross_section_run_period_consistency(const std::string& hayward_csv_path,
                                               const std::string& output_base_dir);

#endif // CROSS_SECTION_RUN_PERIOD_CONSISTENCY_H
