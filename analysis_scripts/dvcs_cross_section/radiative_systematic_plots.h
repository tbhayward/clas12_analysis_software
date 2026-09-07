#ifndef RADIATIVE_SYSTEMATIC_PLOTS_H
#define RADIATIVE_SYSTEMATIC_PLOTS_H

#include <string>

bool make_radiative_systematic_analysis_note_plots(
    const std::string& csv_path,
    const std::string& output_dir =
        "output/systematics/analysis_note/radiative_corrections");

#endif
