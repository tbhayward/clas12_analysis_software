#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <map>
#include <string>

class TTree;

// Update acceptance columns in the CSV and produce acceptance plots.
//
// This refactored version does not loop over ROOT trees. It uses the MC yield
// columns already written by total_counts.cpp:
//
//   generated yield, ep->epg, mc, <period>
//   reconstructed current corrected yield, ep->epg, <topology>, mc, <period>
//
// Acceptance is written as:
//
//   acceptance, <period> = (value,stat,sys)
//
// with sys = 0 for now.
//
// The TTree maps and JSON paths are retained in the signature for interface
// compatibility with the older analysis driver, but are intentionally unused.
bool update_acceptance_csv(const std::string& csv_path,
                           const std::map<std::string, TTree*>& genMcTrees,
                           const std::map<std::string, TTree*>& recMcTrees,
                           const std::string& combined_cuts_json,
                           const std::string& global_cuts_json,
                           const std::string& out_root_dir);

#endif // ACCEPTANCE_H