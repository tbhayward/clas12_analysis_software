#ifndef CURRENT_DEPENDENCE_H
#define CURRENT_DEPENDENCE_H

#include <map>
#include <string>

class TTree;

struct CurrentDependenceOptions {
    // Path to the run charge CSV.
    // Expected columns:
    //   column 1 = run number
    //   column 2 = total RUN::Scaler-like charge
    //   column 3 = positive-helicity charge
    //   column 4 = negative-helicity charge
    //   column 5 = auxiliary scaler component, used only by the optional
    //              columns-3-to-5 scaled mode
    std::string charge_csv_path = "imports/integrated_luminosity/global.csv";

    // Combined 3-sigma cut JSON produced upstream.
    std::string combined_cuts_json = "output/jsons/combined_cuts.json";

    // Output directory for current-dependence diagnostics.
    std::string output_dir = "output/dvcs_current_dependence";

    // Existing override:
    //   false -> compute current-efficiency factors normally
    //   true  -> write all current-efficiency factors as (1,0)
    bool override_to_unity = false;

    // Existing charge-column mode:
    //
    //   true:
    //     use column 2 of the charge import file for all unpolarized data
    //     charge normalization.
    //
    //   false:
    //     use the legacy mixed behavior:
    //       Sp18 Inb / Sp18 Out             -> column 2
    //       Fa18 Inb / Fa18 Out / Sp19 Inb -> column 3 + column 4
    //
    // Spring 2018 is always forced to column 2 in legacy mode.
    //
    // This option is ignored if
    // use_columns_3_to_5_charge_sum_scaled_for_unpolarized is true.
    bool use_second_column_charge_for_all_unpolarized = true;

    // New optional charge-column mode:
    //
    //   true:
    //     use
    //
    //       charge_unpol = columns_3_to_5_charge_sum_scale
    //                      * (column 3 + column 4 + column 5)
    //
    //     for all unpolarized current-dependence data charge normalization.
    //
    //   false:
    //     fall back to use_second_column_charge_for_all_unpolarized / legacy mode.
    //
    // This option has priority over use_second_column_charge_for_all_unpolarized.
    bool use_columns_3_to_5_charge_sum_scaled_for_unpolarized = false;

    // Default scale for the optional columns-3-to-5 mode.
    double columns_3_to_5_charge_sum_scale = 1.025;

    // Default behavior:
    //   true  -> write Sp19 Inb current-efficiency factors using Fa18 Inb
    //            because the Sp19 5 nA luminosity-scan point has suspect
    //            Faraday Cup charge.
    //   false -> use the directly fitted Sp19 Inb current-dependence factor.
    bool use_fa18_inb_current_efficiency_for_sp19_inb = true;

    // Maximum number of OpenMP workers used internally.
    int max_workers = 5;
};

bool update_current_dependence_factors_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& dvcsGenMcTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const CurrentDependenceOptions& options);

#endif // CURRENT_DEPENDENCE_H