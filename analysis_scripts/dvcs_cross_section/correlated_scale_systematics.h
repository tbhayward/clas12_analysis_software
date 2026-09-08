#ifndef CORRELATED_SCALE_SYSTEMATICS_H
#define CORRELATED_SCALE_SYSTEMATICS_H

#include <string>

struct CorrelatedScaleSystematicsOptions {
    double theta_bin_width_deg = 4.0;
    int min_ratio_points_per_period = 25;
    double uncorrelated_normalization_fraction = 0.0476;
    std::string output_dir =
        "output/systematics/correlated_scale";
};

bool correlated_scale_systematics(
    const std::string& csv_path,
    const CorrelatedScaleSystematicsOptions& options =
        CorrelatedScaleSystematicsOptions());

#endif
