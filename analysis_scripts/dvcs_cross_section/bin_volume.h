#ifndef BIN_VOLUME_H
#define BIN_VOLUME_H

#include <map>
#include <string>
#include <vector>

#include "load_binning_scheme.h"  // Binning struct

class TTree; // (unused now, but kept so the signature matches your current main)

/// Computes kinematic bin volumes (independent of generator) and
/// produces per-energy JSON + plots:
///   JSON:  output/jsons/bin_volume_<energy>.json       (energy in {"10.59","10.60","10.2"})
///   Plots: output/bin_volume/<energy>/plot_bin_volume_<energy>_xB_<ix>.png
///
/// NOTE: genMcTrees is no longer used (kept in signature for compatibility).
void compute_and_plot_bin_volume(
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& genMcTrees, // unused
    const std::string& out_root_dir);

#endif // BIN_VOLUME_H