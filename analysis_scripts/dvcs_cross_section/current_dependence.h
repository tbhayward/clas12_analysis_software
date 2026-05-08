#ifndef CURRENT_DEPENDENCE_H
#define CURRENT_DEPENDENCE_H

#include <map>
#include <string>

class TTree;

struct CurrentDependenceOptions {
    std::string charge_csv_path = "imports/integrated_luminosity/global.csv";
    std::string combined_cuts_json = "output/jsons/combined_cuts.json";
    std::string output_dir = "output/dvcs_current_dependence";
    bool override_to_unity = false;
    int max_workers = 5;
};

/**
 * update_current_dependence_factors_csv
 *
 * Performs the current-dependence study adapted from dvcs_current_dependence.py.
 *
 * For each channel and period, it determines:
 *
 *   DATA:
 *     weighted_data_rel = event-weighted fitted data response at the actual
 *                         current mixture divided by fitted zero-current response
 *
 *   MC:
 *     mc_ref_rel = fitted MC efficiency at the reference current divided by
 *                  fitted zero-current MC efficiency
 *
 * and writes:
 *
 *   current efficiency factor, ep->epg,   exp, <period> = (weighted_data_rel,stat)
 *   current efficiency factor, ep->epg,   mc,  <period> = (mc_ref_rel,stat)
 *   current efficiency factor, ep->eppi0, exp, <period> = (weighted_data_rel,stat)
 *   current efficiency factor, ep->eppi0, mc,  <period> = (mc_ref_rel,stat)
 *
 * It also applies the saved MC current-efficiency factors to the reconstructed
 * MC yield columns and writes the current-corrected reconstructed MC columns
 * needed by downstream eppi0 normalization, pi0-contamination, acceptance, and
 * unfolding modules:
 *
 *   reconstructed current corrected yield, ep->epg, mc, <period>
 *   reconstructed current corrected yield, ep->epg, <topology>, mc, <period>
 *   reconstructed current corrected yield, ep->eppi0, mc, <period>
 *   reconstructed current corrected yield, ep->eppi0, <topology>, mc, <period>
 *   reconstructed current corrected yield, ep->eppi0->epg, mc, <period>
 *   reconstructed current corrected yield, ep->eppi0->epg, <topology>, mc, <period>
 *
 * The correction applied is N_rec,current-corrected = N_rec / f_current^MC.
 *
 * If options.override_to_unity is true, no ROOT loops are performed, all
 * current-efficiency factors are written as (1,0), and the current-corrected MC
 * reconstructed-yield columns are copied from the uncorrected reconstructed-yield
 * columns.
 */
bool update_current_dependence_factors_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& dvcsGenMcTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const CurrentDependenceOptions& options);

#endif // CURRENT_DEPENDENCE_H