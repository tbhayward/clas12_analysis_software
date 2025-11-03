#ifndef UNCORRECTED_CROSS_SECTION_H
#define UNCORRECTED_CROSS_SECTION_H

#include <string>
#include <vector>
#include "load_binning_scheme.h"

// Per-helicity uncorrected cross sections (existing):
// Divides unfolded helicity-resolved yields by bin volume and integrated luminosity.
// Writes JSON and per-energy plots.
void compute_uncorrected_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,   // e.g. "output/jsons"
    const std::string& unfolded_counts_dir,   // e.g. "output/jsons"
    const std::string& luminosity_dir,        // e.g. "imports/integrated_luminosity"
    const std::string& output_dir             // e.g. "output/uncorrected_cross_section"
);

// NEW: Unpolarized comparison for Sp18 out vs Fa18 out.
// Uses second column (total) luminosities and sums helicity yields.
// Makes comparison plots only (no JSON).
void compare_unpolarized_cross_sections_sp18out_vs_fa18out(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,   // e.g. "output/jsons"
    const std::string& unfolded_counts_dir,   // e.g. "output/jsons"
    const std::string& luminosity_dir,        // e.g. "imports/integrated_luminosity"
    const std::string& output_dir             // e.g. "output/uncorrected_cross_section"
);

#endif // UNCORRECTED_CROSS_SECTION_H