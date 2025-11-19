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
//   - Compute cross sections for all label/helicity combinations that have
//     both yield and cross section columns in the CSV.
//   - Compute BH/KM/VGG theory curves for 10.6 GeV and 10.2 GeV and write
//     a single xs_phi_all.json per label in output/jsons/cross_sections/.
// Returns true on success, false on any fatal CSV or I/O problem.
bool compute_cross_sections(const std::string &csv_main,
                            const LumiMap &lumi_map);

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

#endif // CROSS_SECTIONS_H