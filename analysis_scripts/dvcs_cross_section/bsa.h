#ifndef BSA_H
#define BSA_H

#include <map>
#include <string>

class TTree;

struct BSAOptions {
    // Existing pass-2 analysis CSV.
    std::string csv_path = "output/csvs/dvcs_pass2_analysis.csv";

    // Combined 3-sigma exclusivity cuts produced by exclusivity_cuts.cpp.
    std::string combined_cuts_json = "output/jsons/combined_cuts.json";

    // Output root. The driver writes plots under <output_root>/bsa_plots and
    // JSON summaries under <output_root>/jsons/BSA_counts.
    std::string output_root = "output";

    // Beam polarizations used for the direct count asymmetry correction.
    double beam_pol_sp18_inb = 0.8882;
    double beam_pol_sp18_out = 0.8882;
    double beam_pol_fa18_inb = 0.8592;
    double beam_pol_fa18_out = 0.8922;
    double beam_pol_sp19_inb = 0.8453;

    // If true, subtract the helicity-separated measured ep->eppi0 contribution
    // using the existing bin-by-bin pi0 leakage scale factors derived from the
    // contamination-ratio and normalized-yield CSV columns:
    //   S+ = G+ - f_pi0 * P+
    //   S- = G- - f_pi0 * P-
    // If false, compute the raw ep->epgamma count asymmetry only.
    bool enable_pi0_subtraction = true;

    // If true, create xB-matrix canvases. Each canvas corresponds to one xB bin,
    // with rows in Q2 and columns in |t|. Each subplot shows A_LU(phi).
    bool make_plots = true;

    // Maximum number of independent data trees processed concurrently.
    // Internally capped at seven workers.
    int max_workers = 7;
};

// Recompute beam-spin asymmetries directly from helicity-split measured counts
// after the standard global cuts and data-derived 3-sigma exclusivity cuts.
//
// The function counts both:
//   G+/- = measured ep->epgamma counts
//   P+/- = measured ep->eppi0 counts
//
// and, when enable_pi0_subtraction is true, computes:
//   S+ = G+ - f_pi0 * P+
//   S- = G- - f_pi0 * P-
//
// The BSA columns filled are:
//   BSA, counts, Fa18 Inb
//   BSA, counts, Fa18 Out
//   BSA, counts, Sp19 Inb
//   BSA, counts, Sp18 Inb
//   BSA, counts, Sp18 Out
//   BSA, counts, Fa18
//   BSA, counts, Sp18
//   BSA, counts, 10.6 GeV
//
// Each filled CSV cell has the usual tuple form:
//   (value,stat_unc,0)
bool update_bsa_counts_csv(const std::map<std::string, TTree*>& dvcsDataTrees,
                           const std::map<std::string, TTree*>& eppi0DataTrees,
                           const BSAOptions& options = BSAOptions());

#endif // BSA_H
