#ifndef BIN_CENTERING_SYSTEMATIC_PLOTS_H
#define BIN_CENTERING_SYSTEMATIC_PLOTS_H
#include <string>
bool make_bin_centering_systematic_analysis_note_plots(
    const std::string& pass1_systematics_csv,
    const std::string& correction_csv,
    const std::string& output_dir = "output/systematics/analysis_note/bin_centering");
#endif
