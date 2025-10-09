#ifndef UNCORRECTED_CROSS_SECTION_H
#define UNCORRECTED_CROSS_SECTION_H

#include <string>
#include <vector>

// Use the project's existing Binning struct
#include "load_binning_scheme.h"

/// Compute uncorrected cross sections (per helicity) by dividing unfolded yields
/// by bin volume and integrated luminosity, and make plots/JSON per energy.
void compute_uncorrected_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,   // e.g. "output/jsons" or "output/bin_volume"
    const std::string& unfolded_counts_dir,   // directory with unfolded_<period>.json (your existing outputs)
    const std::string& luminosity_dir,        // "imports/integrated_luminosity"
    const std::string& output_dir             // "output/uncorrected_cross_section"
);

#endif // UNCORRECTED_CROSS_SECTION_H