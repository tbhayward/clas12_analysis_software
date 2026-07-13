#ifndef SYSTEMATIC_PROJECTION_PLOTS_H
#define SYSTEMATIC_PROJECTION_PLOTS_H

#include <string>

bool make_systematic_projection_plots(
    const std::string& csv_path,
    const std::string& output_dir = "output/systematics/point_to_point_projections");

#endif
