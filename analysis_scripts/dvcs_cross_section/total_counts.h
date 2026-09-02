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
    // versus theta_e), it is evaluated for each event and its common-gradient
    // calibration uncertainty is propagated as a correlated nuisance.
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

#endif // TOTAL_COUNTS_H