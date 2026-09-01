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

    // Controls the large per-period/topology total-count canvases. Keep true
    // for the nominal production workflow; automatic cut variations set this
    // false because only their numerical counts are needed.
    bool make_plots = true;

    // Compact nominal-run summaries intended for the analysis note: raw-yield
    // totals by period/topology, xB projections, and machine-readable tables.
    // Automatic cut-variation jobs should disable these together with make_plots.
    bool make_note_outputs = true;

    // Apply current-dependent reconstruction corrections event-by-event while
    // preserving the original unit-weight raw totals. The response model is
    // produced by current_dependence.cpp before this stage.
    bool apply_event_level_current_correction = true;
    std::string current_response_model_json = "output/dvcs_current_dependence/calibration/current_response_model.json";

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