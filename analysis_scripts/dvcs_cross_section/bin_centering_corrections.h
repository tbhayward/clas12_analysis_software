#ifndef BIN_CENTERING_CORRECTIONS_H
#define BIN_CENTERING_CORRECTIONS_H

#include <string>
#include <vector>

#include "rad_corrected_cross_section.h"  // for Binning definition (or use your binning header if preferred)
#include "model_predictions.h"            // Helicity, ModelPaths, vgg_xs(...), km15_xs(...)

// Model selection for bin-centering corrections
enum class ModelChoice {
    Both,      // Use both models and average them
    VGGOnly,   // Use only VGG model
    KM15Only   // Use only KM15 model
};

/**
 * Compute bin-centered cross sections (no plotting).
 *
 * Inputs:
 *   binning_scheme        : your 3D bin layout
 *   radcorr_xsec_json_dir : directory containing "rad_corrected_xsec_<E>.json" (before-BC input)
 *   output_dir            : base output directory; JSONs are written to "<output_dir>/jsons"
 *   n_steps               : sub-binning granularity per dimension (>=2 recommended)
 *   paths                 : paths to models (VGG, KM15) if needed
 *   vgg_globalfit         : forward --globalfit to dvcsgen for VGG calls
 *   model_choice          : which model(s) to use for bin-centering corrections
 *
 * Output:
 *   JSON files:  <output_dir>/jsons/bin_centered_xsec_<E>.json
 *
 * Note:
 *   This function does not make plots anymore. Use plot_bin_centered_cross_sections(...)
 *   below to generate figures from the saved JSONs.
 */
void compute_bin_centered_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& radcorr_xsec_json_dir,
    const std::string& output_dir,
    int n_steps,
    const ModelPaths& paths,
    bool vgg_globalfit,
    ModelChoice model_choice = ModelChoice::Both
);

/**
 * Plot bin-centered cross sections from precomputed JSONs (fast; no recompute).
 *
 * Inputs:
 *   binning_scheme        : your 3D bin layout
 *   radcorr_xsec_json_dir : directory containing "rad_corrected_xsec_<E>.json" (for overlay: before BC)
 *   bincenter_json_dir    : directory containing "bin_centered_xsec_<E>.json" (after BC; output of compute step)
 *   plots_output_dir      : directory to save PNGs (no JSON written here)
 *
 * Output:
 *   Plots:  <plots_output_dir>/bin_centered_xsec_<E>_xB_<i>.png
 *           <plots_output_dir>/bin_centered_xsec_<E>_xB_<i>_overlay.png
 */
void plot_bin_centered_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& radcorr_xsec_json_dir,
    const std::string& bincenter_json_dir,
    const std::string& plots_output_dir
);

#endif  // BIN_CENTERING_CORRECTIONS_H