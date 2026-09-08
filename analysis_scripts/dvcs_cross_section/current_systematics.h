#ifndef CURRENT_SYSTEMATICS_H
#define CURRENT_SYSTEMATICS_H

#include <string>

struct CurrentSystematicsOptions {
    std::string nuisance_response_csv =
        "output/current_systematics/current_yield_nuisance_responses.csv";
    std::string transfer_shape_csv =
        "output/dvcs_current_dependence/analysis_note/current_systematics/"
        "fa18_sp19_50nA_shape_distance.csv";
    std::string output_dir =
        "output/current_systematics/analysis_note";
    bool write_csv_columns = true;
};

bool evaluate_current_dependence_systematics(
    const std::string& analysis_csv,
    const CurrentSystematicsOptions& options = CurrentSystematicsOptions());

#endif
