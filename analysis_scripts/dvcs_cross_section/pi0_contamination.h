#pragma once
#include <map>
#include <string>
class TTree;

// Overall (helicity-averaged) pi0 contamination calculator.
//
// c = (N_bkg_mc / N_dvcs_csv) * (N_pi0_data / N_pi0_reco_mc)
// with Poisson errors, per CSV row (xB, Q2, |t|, phi).
//
// Inputs:
//   - eppi0DataTrees:  keyed by "<DVCS_tree_key>_eppi0"      (e.g. "DVCS_Fa18_inb_eppi0")
//   - eppi0RecMcTrees: keyed by "<DVCS_tree_key>_rec_mc"     (e.g. "DVCS_Fa18_inb_rec_mc")
//   - eppi0BkgTrees:   keyed by "<DVCS_tree_key>_bkg"        (e.g. "DVCS_Fa18_inb_bkg")
//   - combined_cuts_json: output/jsons/combined_cuts.json
//   - dvcs_csv_path:       output/csvs/dvcs_pass2_analysis.csv
//   - out_root_dir:        "output"
//   - max_workers:         <=0 means auto (omp_get_max_threads)
//
// Behavior:
//   - Periods to run are auto-inferred from CANONICAL_PERIODS() when the three trees exist.
//   - Topologies are always aggregated over (FD, FD), (CD, FD), (CD, FT).
//   - Uses only the CSV for binning and DVCS counts (no external binning scheme).
//   - Writes per-row contamination ratios into: "contamination ratio, <Period Display>".
//
// Output:
//   - Plots: <out_root_dir>/contamination_plots/<PeriodDir>/plot_contamination_<PeriodDir>_xB_<idx>.png
//
// Returns true on success, false on a non-fatal I/O error (fatal issues exit).
bool compute_pi0_contamination_overall(
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::string& combined_cuts_json,
    const std::string& dvcs_csv_path,
    const std::string& out_root_dir,
    int max_workers);