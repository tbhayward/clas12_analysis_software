#pragma once
#include <string>
#include <vector>

// Use the project's existing Binning struct
#include "load_binning_scheme.h"

// Multiply uncorrected cross sections by (Born/Rad) factors and plot.
// Inputs:
//   - binning_scheme: same scheme used everywhere
//   - uncorrected_xsec_json_dir: e.g. "output/uncorrected_cross_section/jsons"
//       files named: "uncorrected_xsec_<E>.json" (E in {"10.59","10.60","10.2"})
//   - radcorr_json_dir: e.g. "output/jsons"
//       files named: "radiative_corrections_group_<E>.json"
//   - output_dir: e.g. "output/rad_corrected_cross_section" (will create subdirs)
void compute_rad_corrected_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& uncorrected_xsec_json_dir,
    const std::string& radcorr_json_dir,
    const std::string& output_dir);