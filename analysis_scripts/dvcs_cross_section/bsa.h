#ifndef BSA_H
#define BSA_H

#include <map>
#include <string>

class TTree;

struct BSAOptions {
    // Existing pass-2 analysis CSV.
    std::string csv_path = "output/csvs/dvcs_pass2_analysis.csv";

    // Canonical nominal (95% containment) production exclusivity cuts
    // produced by plot_exclusivity_data_dvcs_pi0_mc.py.
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

    // Pass-2 pi0 leakage uncertainty established by the contamination study.
    // The BSA systematic is evaluated by repeating the helicity-yield subtraction
    // with f_pi0 scaled by (1 +/- pi0_leakage_relative_uncertainty).
    double pi0_leakage_relative_uncertainty = 0.10;

    // Produce independent FD-photon (FD_FD + CD_FD) and FT-photon (CD_FT)
    // BSA extractions. These are diagnostics only and do not replace the nominal
    // all-topology BSA written to the main CSV.
    bool make_photon_topology_study = true;

    // False-asymmetry study. For each replica, every selected event keeps its
    // run and kinematics but its helicity sign is independently randomized with
    // a deterministic hash. This destroys a physical helicity correlation while
    // retaining the actual run/acceptance population.
    bool make_helicity_scrambling_study = true;
    int helicity_scramble_replicas = 100;
    unsigned long long helicity_scramble_seed = 20260915ULL;

    // Maximum number of independent data trees processed concurrently.
    // Internally capped at seven workers.
    int max_workers = 7;
};

// Recompute beam-spin asymmetries directly from helicity-split measured counts
// after the standard global cuts and Python-optimized nominal production
// exclusivity cuts. Supported production variables are Delta_phi, theta,
// theta_gamma_gamma, pTmiss, Emiss2, Mx2 and Mx2_2. For direct ep->eppi0
// trees, the logical theta_gamma_gamma cut is evaluated with theta_pi0_pi0.
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

// Document the + / - helicity accumulated-charge balance using the same
// authoritative global.csv and final run selections as cross_sections.cpp.
// Periods without usable helicity-resolved Faraday-cup charge (Sp18) are
// reported as unavailable rather than inferred.
bool write_bsa_helicity_charge_balance(
    const std::string& output_csv = "output/bsa_studies/helicity_charge_balance.csv",
    const std::string& charge_csv = "imports/integrated_luminosity/global.csv");

#endif // BSA_H
