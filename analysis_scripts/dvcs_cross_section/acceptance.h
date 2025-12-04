#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <map>
#include <string>

class TTree;

// Update acceptance columns in the CSV and produce acceptance plots.
// Uses non-radiative DVCS MC (generated and reconstructed) and applies:
//   - hard global DVCS cuts on t1, open_angle_ep2, and pTmiss,
//   - run-number–dependent global cuts loaded from global_cuts_json,
//   - topology-dependent 3-sigma cuts loaded from combined_cuts_json.
//
// Returns true on success, false on failure (after printing a diagnostic).
bool update_acceptance_csv(const std::string& csv_path,
                           const std::map<std::string, TTree*>& genMcTrees,
                           const std::map<std::string, TTree*>& recMcTrees,
                           const std::string& combined_cuts_json,
                           const std::string& global_cuts_json,
                           const std::string& out_root_dir);

#endif // ACCEPTANCE_H