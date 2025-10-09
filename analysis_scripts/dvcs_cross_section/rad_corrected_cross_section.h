#pragma once
#include <string>
#include <vector>

// If Binning is declared elsewhere (e.g. load_binning_scheme.h), just ensure
// that header is included before this one in your translation unit.
// We only need the type for plotting grid layout.
struct Binning {
    double xBmin=0.0, xBmax=0.0;
    double Q2min=0.0, Q2max=0.0;
    double tmin=0.0,  tmax=0.0;
};

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