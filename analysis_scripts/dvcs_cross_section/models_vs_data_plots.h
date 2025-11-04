#ifndef MODELS_VS_DATA_PLOTS_H
#define MODELS_VS_DATA_PLOTS_H

#include <string>
#include <vector>

#include "rad_corrected_cross_section.h"  // Binning definition
#include "model_predictions.h"            // Helicity, ModelPaths, vgg_xs(...), km15_xs(...), vgg_bh_only(...)

/*
 * 1) Compute and save model predictions to JSON for later plotting.
 *
 * Inputs:
 *   binning_scheme         : your 3D bin layout
 *   bincenter_json_dir     : directory with bin_centered_xsec_<E>.json (used to know which (ix,iQ,it) exist)
 *   predictions_output_dir : base output directory; JSONs go to <dir>/jsons/model_predictions_<E>.json
 *   paths                  : model paths for VGG/KM15/BH
 *   vgg_globalfit          : pass-through flag for VGG
 *   phi_dense              : number of phi samples (>=2). Suggest 361 for 1-degree steps.
 */
void compute_model_predictions(
    const std::vector<Binning>& binning_scheme,
    const std::string& bincenter_json_dir,
    const std::string& predictions_output_dir,
    const ModelPaths& paths,
    bool vgg_globalfit,
    int phi_dense = 361
);

/*
 * 2) Plot models vs. bin-centered data using saved predictions.
 *
 * Inputs:
 *   binning_scheme       : your 3D bin layout
 *   bincenter_json_dir   : directory with bin_centered_xsec_<E>.json (data points)
 *   predictions_json_dir : directory with model_predictions_<E>.json (precomputed curves)
 *   plots_output_dir     : base output directory (we make "<plots_output_dir>/plots")
 *   phi_dense_unused     : kept for API symmetry; predictions supply their own phi grid
 */
void plot_models_vs_bincentered(
    const std::vector<Binning>& binning_scheme,
    const std::string& bincenter_json_dir,
    const std::string& predictions_json_dir,
    const std::string& plots_output_dir,
    int phi_dense_unused = 0
);

#endif  // MODELS_VS_DATA_PLOTS_H