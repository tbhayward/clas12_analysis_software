#ifndef RAW_SIGNAL_CROSS_CHECK_H
#define RAW_SIGNAL_CROSS_CHECK_H

#include <string>
#include <vector>

struct LeeRow; // from load_csv.h

// Make raw-yield comparison plots (Hayward vs. Lee) and the ratio (Hayward/Lee).
// - Reads "output/jsons/total_counts.json" master (by default) and finds groups
//   whose names contain "fa18_inb" and "fa18_out" (case-insensitive).
// - Uses your CSV rows (already filtered valid==1) to define the grid and match bins.
//
// output_base_dir should be "output/cross_check/lee/raw_yield" (we'll create it).
//
// It produces, for each xB index:
//   raw_counts_fa18_inb_xB_<ix>.png
//   raw_ratio_fa18_inb_xB_<ix>.png
//   raw_counts_fa18_out_xB_<ix>.png
//   raw_ratio_fa18_out_xB_<ix>.png
//
void plot_raw_yield_cross_checks(const std::vector<LeeRow>& rows,
                                 const std::string& total_counts_master_json,
                                 const std::string& output_base_dir);

#endif // RAW_SIGNAL_CROSS_CHECK_H