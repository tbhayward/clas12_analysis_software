#ifndef PLOT_COMPARISON_H
#define PLOT_COMPARISON_H

#include <string>
#include <vector>
#include <map>

void plot_comparison(
    const std::string& output_dir,
    const std::string& bin_boundaries_csv,
    const std::string& previous_data_csv,
    const std::string& new_data_csv
);

#endif // PLOT_COMPARISON_H