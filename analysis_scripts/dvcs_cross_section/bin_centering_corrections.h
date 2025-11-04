#ifndef BIN_CENTERING_CORRECTIONS_H
#define BIN_CENTERING_CORRECTIONS_H

#include <string>
#include <vector>

#include "rad_corrected_cross_section.h"  // for Binning definition and shared helpers, if needed
#include "model_predictions.h"            // Helicity, ModelPaths, vgg_xs(...), km15_xs(...)

// Model selection for bin-centering corrections
enum class ModelChoice {
    Both,      // Use both models and average them (original behavior)
    VGGOnly,   // Use only VGG model
    KM15Only   // Use only KM15 model
};

// Compute bin-centered cross sections starting from the radiatively corrected JSONs.
//
// Inputs:
//   binning_scheme        : your 3D bin layout
//   radcorr_xsec_json_dir : directory that contains "rad_corrected_xsec_<E>.json"
//   output_dir            : base output directory; we create "<output_dir>/jsons" and "<output_dir>/plots"
//   n_steps               : sub-binning granularity per dimension (>=2 recommended; 5 matches your Python)
//   paths                 : paths to models (VGG dvcsgen, KM15 CLI)
//   vgg_globalfit         : forward --globalfit to dvcsgen for VGG calls
//   model_choice          : which model(s) to use for bin-centering corrections
//
// Output:
//   JSON files:  <output_dir>/jsons/bin_centered_xsec_<E>.json
//   Plots:       <output_dir>/plots/bin_centered_xsec_<E>_xB_<i>.png and "..._overlay.png"
//
void compute_bin_centered_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& radcorr_xsec_json_dir,
    const std::string& output_dir,
    int n_steps,
    const ModelPaths& paths,
    bool vgg_globalfit,
    ModelChoice model_choice = ModelChoice::Both
);

#endif  // BIN_CENTERING_CORRECTIONS_H