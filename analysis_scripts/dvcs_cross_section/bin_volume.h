#ifndef BIN_VOLUME_H
#define BIN_VOLUME_H

#include <map>
#include <string>
#include <vector>

#include "load_binning_scheme.h" // struct Binning { double xBmin,xBmax,Q2min,Q2max,tmin,tmax; }
class TTree;

// Computes generator-based bin "volume" vs φ for each (xB,Q²,t) bin.
// Grouping is per beam energy:
//   "10.59": gen trees {sp18_inb_gen, sp18_out_gen}
//   "10.60": gen trees {fa18_inb_gen, fa18_out_gen}
//   "10.2" : gen trees {sp19_inb_gen}
// Output:
//   JSON:  output/jsons/bin_volume_<energy>.json
//   Plots: output/bin_volume/<energy>/plot_bin_volume_<energy>_xB_<ix>.png
void compute_and_plot_bin_volume(
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::string& out_root_dir);

#endif // BIN_VOLUME_H