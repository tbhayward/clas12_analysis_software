#ifndef RADIATIVE_CORRECTIONS_H
#define RADIATIVE_CORRECTIONS_H

#include <map>
#include <string>
#include <vector>
class TTree;

struct Binning {
    double xBmin, xBmax;
    double Q2min, Q2max;
    double tmin,  tmax; // |t| bounds (positive)
};

// Compute per-φ radiative-correction factors R_C(φ) = A_rad(φ)/A_born(φ)
// using reconstructed MC (with MC-side 3σ exclusivity cuts) over generated MC.
void compute_radiative_corrections(
    const std::vector<std::string>& periods,                        // e.g. "DVCS_Fa18_inb", ...
    const std::vector<std::string>& topologies,                     // {"(FD,FD)","(CD,FD)","(CD,FT)"}
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& genMcTrees_norad,          // keys: "sp18_inb_gen", ...
    const std::map<std::string, TTree*>& recMcTrees_norad,          // keys: "sp18_inb_rec", ...
    const std::map<std::string, TTree*>& genMcTrees_rad,            // keys: "sp18_inb_gen_rad", ...
    const std::map<std::string, TTree*>& recMcTrees_rad,            // keys: "sp18_inb_rec_rad", ...
    const std::string& combined_cuts_json_path,                     // "output/jsons/combined_cuts.json"
    const std::string& out_root_dir                                 // "output"
);

#endif // RADIATIVE_CORRECTIONS_H