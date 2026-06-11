#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include <map>
#include <string>

// Simple triple used throughout:
//   value = main quantity
//   stat  = statistical component
//   sys   = systematic component
struct Triple {
    double value;
    double stat;
    double sys;
};

// Map from period/group label to its luminosity triple.
//
// Convention for LumiMap entries:
//   value = total unpolarized accumulated charge
//   stat  = positive-helicity accumulated charge
//   sys   = negative-helicity accumulated charge
using LumiMap = std::map<std::string, Triple>;

struct LumiBuildOptions {
    // Existing charge-column behavior:
    //
    //   true:
    //     use column 2 of the integrated-luminosity import files for all
    //     unpolarized cross-section luminosities.
    //
    //   false:
    //     use the legacy mixed convention:
    //       Sp18 Inb / Sp18 Out             -> column 2
    //       Fa18 Inb / Fa18 Out / Sp19 Inb -> column 3 + column 4
    //
    // This option is ignored if
    // use_columns_3_to_5_charge_sum_scaled_for_unpolarized is true.
    bool use_second_column_charge_for_all_unpolarized = true;

    // New optional unpolarized luminosity mode:
    //
    //   true:
    //     use
    //
    //       L_unpol = columns_3_to_5_charge_sum_scale
    //                 * (column 3 + column 4 + column 5)
    //
    //     for all unpolarized cross-section luminosities.
    //
    //   false:
    //     fall back to use_second_column_charge_for_all_unpolarized / legacy mode.
    //
    // This option has priority over use_second_column_charge_for_all_unpolarized.
    bool use_columns_3_to_5_charge_sum_scaled_for_unpolarized = false;

    // Scale applied in the optional columns-3-to-5 mode.
    // For a 2.5% upward scale, use 1.025.
    double columns_3_to_5_charge_sum_scale = 1.025;
};

// Build luminosity map from imports/integrated_luminosity/.
// Default behavior remains:
//   unpolarized luminosity uses column 2 for all periods.
LumiMap build_lumi_map();

// Build luminosity map with explicit charge-column convention control.
LumiMap build_lumi_map(const LumiBuildOptions& options);

// Update dvcs_pass2_analysis.csv:
//   - Fill integrated luminosity columns using lumi_map.
//   - Compute cross sections from already-corrected, already-unfolded yields:
//       acceptance corrected yield, ep->epg, exp, <label>, <helicity>
//   - Import Frad/Fbin directly from imports/all_bin_v3.csv.
//   - Write those imported Frad/Fbin values into both 10.6 GeV and 10.2 GeV CSV columns.
//   - Read the phase-space-allowed bin_volume from the pass-2 CSV columns
//     filled by bin_volume.cpp.
//   - Do not apply imports/efficiency.json or any additional normalization.
bool compute_cross_sections(const std::string& csv_main,
                            const LumiMap& lumi_map);

// Same as above, but lets the caller explicitly choose the Lee/pass-1 CSV.
bool compute_cross_sections(const std::string& csv_main,
                            const LumiMap& lumi_map,
                            const std::string& lee_csv_path);

// Plot cross sections vs phi for a given label.
bool plot_cross_sections_for_label(const std::string& csv_main,
                                   const std::string& label,
                                   const std::string& theory_json_root,
                                   const std::string& out_root_dir);

// Regenerate BH / KM / VGG theory curves vs phi and write:
//   output/jsons/cross_sections/10.6_GeV/xs_phi_all.json
//   output/jsons/cross_sections/10.2_GeV/xs_phi_all.json
bool regenerate_theory_jsons(const std::string& csv_main,
                             const std::string& theory_json_root);

#endif // CROSS_SECTIONS_H