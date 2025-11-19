#ifndef BIN_CENTERING_CORRECTIONS_H
#define BIN_CENTERING_CORRECTIONS_H

#include <string>
#include "model_predictions.h"

// Compute and write bin-centering corrections Fbin into the analysis CSV.
// - Fbin, 10.6 GeV   (combined 10.6 group: Sp18/Sp18/Fa18/Fa18)
// - Fbin, 10.2 GeV   (Sp19 Inb)
//
// Binning is taken directly from dvcs_pass2_analysis.csv (Lee-style).
// Only rows with non-empty "xBavg, 10.6 GeV" or "xBavg, Sp19 Inb" are used,
// and only phi-binned rows (phimax - phimin < 359).
//
// Fbin is defined via model predictions:
//   - For each CSV bin, we compute
//       Fbin^model = sigma(center) / <sigma>_bin
//     where <sigma>_bin is the average of the model over a uniform
//     sub-grid in (xB, Q2, |t|, phi).
//   - We do this for BOTH KM15 and VGG (helicity unpolarized).
//   - Final Fbin:
//       Fbin = 0.5 * (Fbin_KM15 + Fbin_VGG)
//       stat = 0.0
//       sys  = std(Fbin_KM15, Fbin_VGG)  (sample std for N=2)
//
// Results are written as "(value, stat, sys)" triples into the
// existing columns:
//   "Fbin, 10.6 GeV"
//   "Fbin, 10.2 GeV"
//
// Parallelization:
//   - Uses OpenMP across bins, hard-capped to 5 threads maximum.

bool update_bin_centering_corrections_csv(
    const std::string& csv_path,
    const std::string& out_root_dir,
    int n_substeps,
    const ModelPaths& model_paths,
    bool vgg_globalfit
);

#endif // BIN_CENTERING_CORRECTIONS_H