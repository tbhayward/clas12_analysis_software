#ifndef EPPI0_NORMALIZATION_H
#define EPPI0_NORMALIZATION_H

#include <map>
#include <string>

class TTree;

struct Eppi0NormalizationOptions {
    std::string charge_csv_path = "imports/integrated_luminosity/global.csv";
    // Canonical nominal (95% containment) Python-optimized production cuts.
    std::string combined_cuts_json = "output/jsons/combined_cuts.json";
    std::string output_dir = "output/data_mc_normalization";
    std::string normalization_json_path = "imports/eppi0_aao_normalization_inputs.json";
    // Use the same event-by-event regional current-response calibration as total_counts.cpp.
    std::string current_response_model_json =
        "output/dvcs_current_dependence/calibration/current_response_model.json";
    bool override_to_unity = false;
    // Remove stale diagnostic PNG/text outputs before a new study run.
    bool clean_output_dir = true;

    // Study mode: derive the pi0 DATA/AAOGEN normalization, write fits,
    // diagnostic plots, and machine-readable summaries, but do NOT overwrite
    // the production normalized-yield columns. This is the default while the
    // absolute pi0-model and photon-efficiency systematics are being validated.
    bool write_normalized_yields = false;

    // Write one compact CSV containing the period/region integrated ratios and
    // fitted cubic coefficients. This is intended for fast external studies.
    bool write_summary_csv = true;
    std::string summary_csv_path =
        "output/data_mc_normalization/eppi0_normalization_summary.csv";

    // Export the accepted reconstructed-AAOgen population used by this study
    // as a sparse fine-grained (xB,Q2,-t,photon-topology) weighted grid.  This
    // is intended for external pi0-model validation/folding without creating a
    // potentially enormous event-level CSV.  The grid is filled in the same MC
    // loop, after the same eppi0 selection, with the same event_norm and
    // event-by-event current-response weight used by the normalization study.
    bool write_accepted_mc_population = true;
    std::string accepted_mc_population_csv_path =
        "output/data_mc_normalization/accepted_aao_population.csv";

    int max_workers = 5;
};

/**
 * update_eppi0_normalization_csv
 *
 * Uses eppi0 DATA, reconstructed AAOGEN MC, and period normalization metadata
 * from options.normalization_json_path to determine period-dependent regional
 * cubic normalization functions
 *
 *   R_pi0(theta_p) = N_data(theta_p) / N_MC(theta_p)
 *
 * after applying the same event-by-event regional current-response calibration
 * used by the production total-counts analysis. Krishna/Neupane proton-efficiency
 * weights are intentionally NOT applied in this normalization diagnostic.
 *
 * The seven cubic fits, one for each p1 FD sector and one CD region, are written to:
 *
 *   eppi0 cross-section normalization cubic, ep->eppi0, data_over_mc, <region>, <period>
 *
 * as:
 *
 *   (a0,a1,a2,a3)
 *
 * with theta_p in degrees and
 *
 *   R_pi0(theta_p) = a0 + a1 theta + a2 theta^2 + a3 theta^3.
 *
 * If options.write_normalized_yields is true, the function then loops over DVCS
 * and eppi0 DATA events and fills:
 *
 *   normalized raw yield, ep->epg,   <topo>, exp, <period>, <helicity>
 *   normalized raw yield, ep->eppi0, <topo>, exp, <period>, <helicity>
 *
 * using event weight:
 *
 *   w = 1 / (current_efficiency_exp_factor * R_pi0(theta_p)).
 *
 * If override_to_unity is true, all cubic cells are written as
 * (1,0,0,0), all scalar normalization factor cells as (1,0), and the
 * normalized raw yields are filled using only the current-efficiency factor.
 */
bool update_eppi0_normalization_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const Eppi0NormalizationOptions& options);

#endif // EPPI0_NORMALIZATION_H