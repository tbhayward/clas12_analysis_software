#ifndef BIN_CENTERING_CORRECTIONS_H
#define BIN_CENTERING_CORRECTIONS_H

#include <string>
#include "model_predictions.h"

// Choice of which models to use for Fbin construction.
enum class ModelChoice {
    Both = 0,
    VGGOnly,
    KM15Only
};

// Compute bin-centering corrections Fbin for each phi-binned CSV row and
// write them into:
//   "Fbin, 10.6 GeV"
//   "Fbin, 10.2 GeV"
//
// Inputs:
//   csv_path      : dvcs_pass2_analysis.csv (absolute or relative path)
//   n_steps       : sub-binning per dimension (>=2). Total model calls
//                   per row per model is n_steps^4.
//   paths         : ModelPaths (dvcsgen dir, km15_cli path); empty -> env/defaults.
//   vgg_globalfit : whether to pass --globalfit flag to vgg_xs.
//   model_choice  : Both, VGGOnly, or KM15Only.
//
// Returns true on success, false on failure.
bool update_bin_centering_corrections_csv(
    const std::string& csv_path,
    int n_steps,
    const ModelPaths& paths,
    bool vgg_globalfit,
    ModelChoice model_choice);

// Debug plots: Fbin vs phi for each (xB,Q2,t) bin, for 10.6 GeV and 10.2 GeV.
// Reads:
//   - xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax
//   - xBavg, 10.6 GeV / xBavg, Sp19 Inb
//   - phiavg, 10.6 GeV / phiavg, Sp19 Inb
//   - Fbin, 10.6 GeV / Fbin, 10.2 GeV   (triples "(value, stat, sys)")
//
// out_root_dir is typically "output/bin_volume". This function will create
// subdirs "10.60" and "10.2" under out_root_dir if needed.
void plot_bin_centering_fbin_vs_phi(
    const std::string& csv_path,
    const std::string& out_root_dir);

#endif // BIN_CENTERING_CORRECTIONS_H