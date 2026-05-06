// yield_totals.h
#ifndef YIELD_TOTALS_H
#define YIELD_TOTALS_H

#include <map>
#include <string>

class TTree;

// Compute run-period/current summaries for the final normalized data yields.
//
// The two reported quantities are:
//
//   1) DVCS normalized pi0-subtracted counts:
//        sum over accepted ep->epgamma DATA events of
//          [1 / (current_efficiency_factor(epg,exp,period) * R_pi0(p1_theta))]
//          * [1 - contamination_ratio(row,period)]
//
//   2) eppi0 normalized counts:
//        sum over accepted ep->eppi0 DATA events of
//          [1 / (current_efficiency_factor(eppi0,exp,period) * R_pi0(p1_theta))]
//
// Events are grouped by period and beam current using runnum -> current maps.
// Global cuts, topology selection, and channel-specific 3-sigma cuts are applied
// before an event contributes to a bin/current total.
//
// The old reconstructed-MC argument was removed because this report is now a
// data-current diagnostic for the normalized analysis counts.
bool compute_yield_totals(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dvcsDataTrees,
                          const std::map<std::string, TTree*>& eppi0DataTrees,
                          const std::string& combined_cuts_json,
                          const std::string& output_txt);

#endif // YIELD_TOTALS_H
