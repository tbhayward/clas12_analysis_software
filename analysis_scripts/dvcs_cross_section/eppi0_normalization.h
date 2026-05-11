#ifndef EPPI0_NORMALIZATION_H
#define EPPI0_NORMALIZATION_H

#include <map>
#include <string>

class TTree;

struct Eppi0NormalizationOptions {
    std::string charge_csv_path = "imports/integrated_luminosity/global.csv";
    std::string combined_cuts_json = "output/jsons/combined_cuts.json";
    std::string output_dir = "output/data_mc_normalization";
    bool override_to_unity = false;
    int max_workers = 5;

    // Optional diagnostic/testing override. Leave empty for normal production.
    // Keys must be the display period labels, e.g. "Fa18 Inb".
    // If a period is present with value > 0, that value is used instead of
    // gen_tree->GetEntries() for N_gen in the AAOGEN normalization.
    std::map<std::string, double> n_gen_override_by_period;
};

/**
 * update_eppi0_normalization_csv
 *
 * Uses eppi0 DATA and AAOGEN MC to determine a period-dependent quartic
 * normalization function
 *
 *   R_pi0(theta_p) = N_data(theta_p) / N_MC(theta_p)
 *
 * after applying the current-efficiency factors already saved in the CSV.
 *
 * The quartic coefficients are written to:
 *
 *   eppi0 cross-section normalization polynomial, ep->eppi0, data_over_mc, <period>
 *
 * as:
 *
 *   (a0,a1,a2,a3,a4)
 *
 * with theta_p in degrees and
 *
 *   R_pi0(theta_p) = a0 + a1 theta + a2 theta^2 + a3 theta^3 + a4 theta^4.
 *
 * Then the function loops over DVCS and eppi0 DATA events and fills:
 *
 *   normalized raw yield, ep->epg,   <topo>, exp, <period>, <helicity>
 *   normalized raw yield, ep->eppi0, <topo>, exp, <period>, <helicity>
 *
 * using event weight:
 *
 *   w = 1 / (current_efficiency_exp_factor * R_pi0(theta_p)).
 *
 * If override_to_unity is true, all polynomial cells are written as
 * (1,0,0,0,0), all scalar normalization factor cells as (1,0), and the
 * normalized raw yields are filled using only the current-efficiency factor.
 */
bool update_eppi0_normalization_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const Eppi0NormalizationOptions& options);

#endif // EPPI0_NORMALIZATION_H
