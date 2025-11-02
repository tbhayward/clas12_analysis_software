#ifndef BSA_H
#define BSA_H

#include <map>
#include <string>
#include <vector>

class TTree;

// Forward declaration of your binning struct (as in your codebase)
struct Binning;

// Legacy entry point (kept for backward compatibility).
// Calls the 7-argument variant with an empty radcorr_xsec_json_dir.
void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>&     binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& pi0_corrected_counts_json_path,
    const std::string& out_root_dir);

// New entry point that also accepts the directory containing
// rad_corrected_xsec_<E>.json files for the overlay plots.
void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>&     binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& pi0_corrected_counts_json_path,
    const std::string& out_root_dir,
    const std::string& radcorr_xsec_json_dir);

#endif // BSA_H