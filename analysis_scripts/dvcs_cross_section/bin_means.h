#ifndef BIN_MEANS_H
#define BIN_MEANS_H

#include "load_binning_scheme.h"

#include <TTree.h>

#include <map>
#include <string>
#include <vector>

// Compute global bin-averaged kinematics (xB, Q2, |t|, phi) over multiple DVCS periods & topologies.
// - dvcs_periods examples: "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb", "DVCS_Sp18_out", "DVCS_Sp18_inb", "DVCS_Fa18_inb_supp"
// - topologies examples: "(FD,FD)", "(CD,FD)", "(CD,FT)"
// - analysis_type currently unused but kept for parity with your python
// - binning_scheme from load_binning_scheme()
// - output_json path for the resulting JSON
// - dataTrees map holding DVCS data trees keyed by run tags like "fa18_inb", "sp18_out", "fa18_inb_supp", "sp19_inb", etc.
void calculate_bin_means(
    const std::vector<std::string>& dvcs_periods,
    const std::vector<std::string>& topologies,
    const std::string& analysis_type,
    const std::vector<Binning>& binning_scheme,
    const std::string& output_json,
    const std::map<std::string, TTree*>& dataTrees
);

#endif // BIN_MEANS_H