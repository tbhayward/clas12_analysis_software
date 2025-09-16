#ifndef TOTAL_COUNTS_H
#define TOTAL_COUNTS_H

#include "load_binning_scheme.h"  // for Binning definition
#include <TTree.h>
#include <map>
#include <string>
#include <vector>

// Compute total DVCS counts in (xB, Q2, |t|, phi) bins after cuts,
// split by beam helicity (+1/-1). Cuts are loaded from
// output/jsons/combined_cuts.json (format written by exclusivity_cuts.cpp).
//
// periods: names like "DVCS_Sp18_inb", "DVCS_Fa18_out", etc.
// topologies: strings "(FD,FD)", "(CD,FD)", "(CD,FT)" used to map detector1/2.
// binning_scheme: list of 4D bin definitions (same CSV you use elsewhere).
// dataTrees: map run-tag -> TTree*, where run-tag is lowercase like "sp18_inb".
// cuts_json_path: typically "output/jsons/combined_cuts.json".
// output_json_path: where to write the counts JSON.
void compute_total_counts(const std::vector<std::string>& periods,
  const std::vector<std::string>& topologies, const std::vector<Binning>& binning_scheme,
  const std::map<std::string, TTree*>& dataTrees, const std::string& cuts_json_path,
  const std::string& output_json_path);

#endif // TOTAL_COUNTS_H