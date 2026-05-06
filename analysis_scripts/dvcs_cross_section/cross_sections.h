#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include <map>
#include <string>

// Simple triple used throughout:
//   value = main quantity
//   stat  = statistical component
//   sys   = systematic component (often 0 here)
struct Triple {
    double value;
    double stat;
    double sys;
};

// Map from period/group label to its luminosity triple.
// Convention for LumiMap entries:
//   value = total (unpolarized) accumulated charge
//   stat  = + helicity accumulated charge
//   sys   = - helicity accumulated charge
using LumiMap = std::map<std::string, Triple>;

// Build luminosity map from RGA text files in imports/integrated_luminosity/.
// Uses the conventions described in cross_sections.cpp comments.
LumiMap build_lumi_map();

// Update dvcs_pass2_analysis.csv:
//   - Fill integrated luminosity columns using lumi_map.
//   - Compute cross sections from already-corrected, already-unfolded yields:
//       acceptance corrected yield, ep->epg, exp, <label>, <helicity>
//   - Import Frad/Fbin/bin_volume directly from imports/all_bin_v3.csv.
//   - Write those imported values into both 10.6 GeV and 10.2 GeV CSV columns.
//   - Do not apply imports/efficiency.json or any additional normalization.
// Returns true on success, false on any fatal CSV or I/O problem.
bool compute_cross_sections(const std::string &csv_main,
                            const LumiMap &lumi_map);

// Same as above, but lets the caller explicitly choose the Lee/pass-1 CSV.
// The expected source columns are:
//   Frad
//   Fbin
//   bin_volume
// matched by the pass-2 CSV "bin index" column.
bool compute_cross_sections(const std::string &csv_main,
                            const LumiMap &lumi_map,
                            const std::string &lee_csv_path);

// Plot cross sections vs phi for a given label:
//   - Reads cross section columns for that label from csv_main.
//   - Groups rows by xB range, then lays out subpads in Q2 x |t| grid.
//   - Uses log scale on Y.
//   - Uses a single legend at the top of the canvas.
//   - Overlays BH/KM/VGG theory curves if a xs_phi_all.json exists for
//     this label under theory_json_root.
// Saves PNGs under out_root_dir/<PeriodDir>/.
//
// Returns true on success. If the required cross section columns for this
// particular label do not exist, it prints a message and returns true
// (nothing to plot for that label).
bool plot_cross_sections_for_label(const std::string &csv_main,
                                   const std::string &label,
                                   const std::string &theory_json_root,
                                   const std::string &out_root_dir);


// Regenerate BH / KM / VGG theory curves vs phi and write:
//   output/jsons/cross_sections/10.6_GeV/xs_phi_all.json
//   output/jsons/cross_sections/10.2_GeV/xs_phi_all.json
//
// This function overwrites any existing files. It does NOT modify the CSV.
bool regenerate_theory_jsons(const std::string &csv_main,
                             const std::string &theory_json_root);

#endif // CROSS_SECTIONS_H