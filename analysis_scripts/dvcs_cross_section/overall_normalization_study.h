#ifndef OVERALL_NORMALIZATION_STUDY_H
#define OVERALL_NORMALIZATION_STUDY_H

#include <string>

// Runtime choice for the variable used to write norm, <label> after the BH-edge fit.
enum class OverallNormXAxis {
    PTheta = 0,
    XB = 1
};

struct OverallNormalizationOptions {
    // If true, skip all BH calculations and plots and write (or overwrite)
    // norm, <label> = 1.00 in every CSV row.
    bool override_to_unity = true;

    // If true, use all edge-window points. If false, use the closest-to-edge
    // point per (xB,Q2,|t|) bin, still subject to the edge-window cut.
    bool use_all_points_within_edge_window = true;

    // If true, reject points exactly at d_edge = 0.
    bool require_positive_dedge = true;

    // Edge-window definition in degrees:
    //   d_edge = min(|phi - 0|, |phi - 360|)
    double max_dedge_for_normalization_deg = 10.0;

    // Variable used for the fitted normalization written to norm, <label>.
    OverallNormXAxis norm_x_axis = OverallNormXAxis::XB;

    // Output directory for diagnostic BH-normalization plots.
    std::string output_dir = "output/normalization_study";
};

// Updates norm, <label> in the CSV.
//
// In override mode, writes norm, <label> = 1.00 everywhere and returns without
// requiring model/cross-section columns.
//
// In study mode, reads the already-computed cross-section columns:
//   cross sections, ep->epg, exp, <label>, <helicity>
// compares the edge-window points against the BH model, produces diagnostic
// plots, fits xs/BH versus the chosen x variable, and writes norm, <label>.
bool update_overall_normalization_study_csv(const std::string &csv_path,
                                            const std::string &label,
                                            const std::string &helicity,
                                            const OverallNormalizationOptions &options);

// Backward-compatible wrapper. This intentionally defaults to unity override,
// matching the production-safe behavior.
bool print_bh_normalization_study(const std::string &csv_path,
                                  const std::string &label,
                                  const std::string &helicity);

#endif // OVERALL_NORMALIZATION_STUDY_H
