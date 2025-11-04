#ifndef MODELS_VS_DATA_PLOTS_H
#define MODELS_VS_DATA_PLOTS_H

#include <string>
#include <vector>

#include "rad_corrected_cross_section.h"  // for Binning definition and shared helpers
#include "model_predictions.h"            // Helicity, ModelPaths, vgg_xs(...), km15_xs(...), bh_xs(...)

/*
 * Make model–data comparison plots for bin-centered cross sections:
 *   • Data points (bin-centered):  + helicity = blue circles,  − helicity = red squares
 *   • KM15 model curves:           + helicity dashed blue,     − helicity dashed red
 *   • VGG model curves:            + helicity dotted blue,     − helicity dotted red
 *   • BH exact curve (helicity-indep.): dashed black
 *
 * Inputs:
 *   binning_scheme         : your 3D bin layout
 *   bincenter_json_dir     : directory with bin_centered_xsec_<E>.json
 *   output_dir             : base output directory (we make "<output_dir>/plots")
 *   paths                  : model paths for VGG/KM15/BH
 *   vgg_globalfit          : pass-through flag for VGG (e.g., --globalfit)
 *   phi_dense              : number of phi samples (>= 73 recommended). Default 361 (0..360 by 1 deg)
 *
 * Notes:
 *   - Energies attempted: 10.59, 10.60, 10.2 (missing files are skipped with a note).
 *   - ROOT drawing happens single-threaded; model sampling can be multithreaded if desired.
 *   - We do not modify any JSONs; this only plots.
 */
void plot_models_vs_bincentered(
    const std::vector<Binning>& binning_scheme,
    const std::string& bincenter_json_dir,
    const std::string& output_dir,
    const ModelPaths& paths,
    bool vgg_globalfit,
    int phi_dense = 361  // 1-degree steps: 0..360 inclusive
);

#endif  // MODELS_VS_DATA_PLOTS_H