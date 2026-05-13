#ifndef CROSS_SECTION_HAYWARD_CROSS_CHECK_H
#define CROSS_SECTION_HAYWARD_CROSS_CHECK_H

#include <string>

// -----------------------------------------------------------------------------
// Cross-check of experimental cross sections *within* Hayward pass-2,
// comparing different run periods / helicities.
//
// Reads only dvcs_pass2_analysis.csv and produces overlay + ratio plots
// in subdirectories under output_base_dir:
//
//  Subdirs (examples):
//    Fa18Inb_vs_Fa18Out_unpol
//    Fa18Inb_vs_Sp18Inb_unpol
//    Fa18Out_vs_Sp18Out_unpol
//    Sp18Inb_vs_Sp18Out_unpol
//    Fa18_vs_Sp18_unpol
//    Fa18Inb_vs_Fa18Out_pos
//    Fa18Inb_vs_Fa18Out_neg
//
// Sp18 is intentionally unpolarized-only. No Sp18 pos/neg comparison is
// produced or required.
//
//  Each subdir contains:
//    cross_section_counts_xB_<ix>.png
//    cross_section_ratio_xB_<ix>.png
//
// -----------------------------------------------------------------------------

void plot_cross_section_hayward_cross_checks(const std::string& hayward_csv_path,
                                             const std::string& output_base_dir);

#endif // CROSS_SECTION_HAYWARD_CROSS_CHECK_H