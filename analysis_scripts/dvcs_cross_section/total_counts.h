#ifndef TOTAL_COUNTS_H
#define TOTAL_COUNTS_H

#include <map>
#include <string>

class TTree;

struct TotalCountsOptions {
    // The combined-cuts JSON is expected to be produced by
    // plot_exclusivity_data_dvcs_pi0_mc.py. Reconstructed samples apply only
    // the production variables accepted for their period/topology. Mx2 is
    // mandatory; any accepted cut whose ROOT branch is missing is fatal.
    // If true, DVCS generated/reconstructed MC count columns are filled from
    // no-background dvcsgen files instead of the production-current/background
    // overlaid files. This affects only ep->epg MC counts; data, ep->eppi0 MC,
    // and ep->eppi0->epg background MC are unchanged.
    bool use_nobkg_dvcs_mc_counts = false;

    // Controls the large per-period/topology total-count canvases. The compact
    // analysis-note summaries contain the production information we retain, so
    // these verbose canvases are disabled by default. They can still be enabled
    // for a dedicated debugging run.
    bool make_plots = false;

    // Compact nominal-run summaries intended for the analysis note: raw-yield
    // totals by period/topology, xB projections, and machine-readable tables.
    // Automatic cut-variation jobs should disable these together with make_plots.
    bool make_note_outputs = true;

    // Apply current-dependent reconstruction corrections event-by-event while
    // preserving the original unit-weight raw totals. The response model is
    // produced by current_dependence.cpp before this stage.  Regional DATA
    // responses are used by default; when the calibration JSON contains an
    // optional centered polar-angle term (nominally Sp18 Out ep->epgamma
    // versus theta_e), it is evaluated for each event.  The ordinary event
    // counting uncertainty contains only sum(w^2); fitted response-parameter
    // uncertainties are written separately as correlated calibration nuisances.
    bool apply_event_level_current_correction = true;
    std::string current_response_model_json = "output/dvcs_current_dependence/calibration/current_response_model.json";

    // Final-production guard: the reviewed current-dependence study requires
    // the centered linear e_theta term for Sp18 Out ep->epg DATA. When true,
    // total_counts fails rather than silently reverting to regional-only if
    // that angular model is absent or malformed in the calibration JSON.
    bool require_sp18_out_epg_e_theta_current_model = true;

    // Misidentified ep->eppi0->epg reconstructed MC follows the ep->epg
    // regional response because its reconstructed topology is epgamma.
    bool use_epg_mc_current_factor_for_eppi0_bkg = true;

    // Write the signed one-standard-deviation response of each current-corrected
    // yield to every fitted current-response calibration parameter.  These
    // responses are consumed by current_systematics.cpp and preserve the
    // correlated nuisance structure across physics bins.
    bool write_current_nuisance_responses = true;
    std::string current_nuisance_response_csv =
        "output/current_systematics/current_yield_nuisance_responses.csv";
};

/**
 * update_total_counts_csv
 *
 * Updates the pass-2 analysis CSV with:
 *
 * DATA:
 *   raw yield, ep->epg,   <topo>, exp, <period>, <helicity>
 *   raw yield, ep->eppi0, <topo>, exp, <period>, <helicity>
 *
 * MC:
 *   generated yield, ep->epg, mc, <period>
 *   reconstructed yield, ep->epg, mc, <period>
 *   reconstructed yield, ep->epg, <topo>, mc, <period>
 *
 *   generated yield, ep->eppi0, mc, <period>
 *   reconstructed yield, ep->eppi0, mc, <period>
 *   reconstructed yield, ep->eppi0, <topo>, mc, <period>
 *
 *   reconstructed yield, ep->eppi0->epg, mc, <period>
 *   reconstructed yield, ep->eppi0->epg, <topo>, mc, <period>
 *
 * The nominal workflow also fills the normalized DATA and reconstructed-current-
 * corrected MC columns event-by-event using the regional current-response model
 * produced by current_dependence.cpp. The original raw/generated/reconstructed
 * unit-weight columns remain unchanged and continue to feed the raw-yield note
 * outputs.
 *
 * Return:
 *   true on success, false on fatal error.
 */
bool update_total_counts_csv(const std::string& csv_path,
                             const std::map<std::string, TTree*>& dvcsDataTrees,
                             const std::map<std::string, TTree*>& eppi0DataTrees,
                             const std::map<std::string, TTree*>& dvcsGenMcTrees,
                             const std::map<std::string, TTree*>& dvcsRecMcTrees,
                             const std::map<std::string, TTree*>& eppi0GenMcTrees,
                             const std::map<std::string, TTree*>& eppi0RecMcTrees,
                             const std::map<std::string, TTree*>& eppi0BkgTrees,
                             const std::string& combined_cuts_json,
                             const std::string& out_root_dir,
                             int max_workers,
                             const TotalCountsOptions& options = TotalCountsOptions(),
                             const std::map<std::string, TTree*>& dvcsNoBkgGenMcTrees = std::map<std::string, TTree*>(),
                             const std::map<std::string, TTree*>& dvcsNoBkgRecMcTrees = std::map<std::string, TTree*>());


struct AcceptanceReweightingOptions {
    // Input/output locations.
    std::string combined_cuts_json = "output/jsons/combined_cuts.json";
    std::string current_response_model_json =
        "output/dvcs_current_dependence/calibration/current_response_model.json";
    std::string output_dir = "output/systematics/acceptance_reweighting";

    // Iterative data-driven reweighting controls.
    int max_iterations = 5;
    double damping_power = 0.50;
    double per_iteration_weight_min = 0.50;
    double per_iteration_weight_max = 2.00;
    double cumulative_weight_min = 0.05;
    double cumulative_weight_max = 20.0;
    double convergence_shape_distance = 0.010;
    double convergence_improvement = 0.001;
    int minimum_entries_per_fine_bin = 50;

    // The fine 1D axes are constructed by splitting every nominal analysis
    // interval in half.  This is intentionally finer than the production bins:
    // a weight that is constant across one production bin would cancel exactly
    // in N_rec/N_gen and could not test acceptance-model dependence.
    int subdivisions_per_nominal_interval = 2;

    // Approximate background subtraction when deriving the DATA target shape:
    // each selected DATA event is multiplied by (1-c_pi0) for its production
    // four-dimensional bin.  The pass-2 contamination column is already
    // available in the nominal CSV.
    bool apply_pi0_signal_fraction_to_data = true;

    // Apply the already-finalized event-level current correction to DATA and
    // reconstructed MC before fitting the shape weights.
    bool apply_current_correction = true;

    // Optional pure-BH stress test.  It is disabled by default and is NOT used
    // in the candidate acceptance systematic.  The production candidate is
    // determined solely from the nominal -> DATA-reweighted acceptance change.
    bool enable_bh_reweighting = false;
    bool build_bh_grid_if_missing = false;
    std::string bh_grid_script =
        "external_scripts/build_acceptance_bh_reweight_grid.py";
    std::string bh_grid_csv =
        "output/systematics/acceptance_reweighting/bh_subcell_grid.csv";
    int bh_grid_workers = 7;

    // First-pass safety: write the candidate pass-2 acceptance systematic and
    // all diagnostics, but do not overwrite Syst. err (Acceptance) until the
    // diagnostic result has been reviewed.
    bool install_candidate_as_production_systematic = false;
};

// Re-evaluate the acceptance model dependence using the existing pass-2
// dvcsgen generated/reconstructed MC and nominal DATA trees.
//
// The study constructs:
//   (1) nominal acceptance A0,
//   (2) acceptance Adata after iterative DATA-driven shape reweighting.
//
// The candidate model uncertainty is |Adata-A0|/A0.  A pure-BH stress test may
// still be enabled explicitly, but it is not used in the candidate systematic.
// The study also performs synthetic-MC transfer-closure tests with known smooth
// weights over several distortion amplitudes.  These quantify the residual from
// learning a weight at reconstructed level and applying it to the generated
// denominator.  Because the present trees do not provide a generated->reconstructed
// event map, this is explicitly NOT claimed as a full detector-response closure.
// The conservative candidate combines the DATA-driven excursion with the
// row-level transfer-closure residual.  It is NOT installed into
// Syst. err (Acceptance) unless explicitly requested.
bool run_acceptance_reweighting_study(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsGenMcTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const AcceptanceReweightingOptions& options =
        AcceptanceReweightingOptions());

#endif // TOTAL_COUNTS_H