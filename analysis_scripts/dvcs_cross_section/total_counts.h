#ifndef TOTAL_COUNTS_H
#define TOTAL_COUNTS_H

#include <map>
#include <string>

class TTree;

struct TotalCountsOptions {
    // If true, DVCS generated/reconstructed MC count columns are filled from
    // no-background dvcsgen files instead of the production-current/background
    // overlaid files. This affects only ep->epg MC counts; data, ep->eppi0 MC,
    // and ep->eppi0->epg background MC are unchanged.
    bool use_nobkg_dvcs_mc_counts = false;

    // Controls the large per-period/topology total-count canvases. Keep true
    // for the nominal production workflow; automatic cut variations set this
    // false because only their numerical counts are needed.
    bool make_plots = true;
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
 * The corresponding "reconstructed current corrected yield" columns are
 * created by initialize_pass2_csv but intentionally not filled here. They should
 * be filled later after current-efficiency correction factors have been
 * determined.
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