// pi0_corrected_counts.h
#pragma once
#include <string>
#include <vector>
struct Binning;

void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,
    const std::string& contamination_dir,
    const std::string& out_root_dir
);