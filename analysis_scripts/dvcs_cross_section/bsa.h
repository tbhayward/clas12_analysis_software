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

    // If true, create one phi-dependence canvas per populated xB,Q2,|t| cell.
    bool make_plots = true;

    // Placeholder for future parallelization. The current implementation is
    // intentionally serial because ROOT branch binding is global-state-heavy and
    // the BSA loop is cheap compared with the other production stages.
    int max_workers = 1;
};

// Recompute beam-spin asymmetries directly from helicity-split DVCS data counts
// after the standard global cuts and DVCS 3-sigma exclusivity cuts. This fills:
//   BSA, counts, Fa18 Inb
//   BSA, counts, Fa18 Out
//   BSA, counts, Sp19 Inb
//   BSA, counts, Sp18 Inb
//   BSA, counts, Sp18 Out
//   BSA, counts, Fa18
//   BSA, counts, Sp18
//   BSA, counts, 10.6 GeV
// Each filled CSV cell has the usual tuple form: (value,stat_unc,0).
bool update_bsa_counts_csv(const std::map<std::string, TTree*>& dvcsDataTrees,
                           const BSAOptions& options = BSAOptions());

#endif // BSA_H
