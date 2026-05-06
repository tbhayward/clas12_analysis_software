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
 * It determines the DVCS data and MC current factors directly.
 *
 * It also determines the eppi0 DATA current factor directly from the eppi0
 * data luminosity/current scan. Because there is no eppi0 luminosity-scan MC,
 * the eppi0 MC factor is derived from:
 *
 *   eppi0_mc_factor = eppi0_data_factor * (dvcs_mc_factor / dvcs_data_factor)
 *
 * and the statistical uncertainty is propagated from the eppi0 data factor,
 * the DVCS MC factor, and the DVCS data factor.
 *
 * It writes:
 *
 *   current efficiency factor, ep->epg,   exp, <period> = (dvcs_data_factor,stat)
 *   current efficiency factor, ep->epg,   mc,  <period> = (dvcs_mc_factor,stat)
 *   current efficiency factor, ep->eppi0, exp, <period> = (eppi0_data_factor,stat)
 *   current efficiency factor, ep->eppi0, mc,  <period> = (derived_eppi0_mc_factor,stat)
 *
 * If options.override_to_unity is true, no ROOT loops are performed and all
 * current-efficiency factors are written as (1,0).
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