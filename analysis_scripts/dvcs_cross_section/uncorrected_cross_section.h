#pragma once
#include <string>
#include <map>
#include <vector>

struct Binning {
    double xBmin, xBmax;
    double Q2min, Q2max;
    double tmin, tmax;
};

void compute_uncorrected_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,
    const std::string& unfolded_counts_dir,
    const std::string& luminosity_dir,
    const std::string& output_dir);