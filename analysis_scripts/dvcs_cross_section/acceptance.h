#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <map>
#include <string>

class TTree;

// Update acceptance columns in the CSV and produce acceptance plots.
// Uses non-radiative DVCS MC (generated and reconstructed) and applies
// global DVCS cuts plus 3 sigma cuts loaded from cuts_json.
// Returns true on success, false on failure (after printing a diagnostic).
bool update_acceptance_csv(const std::string& csv_path,
                           const std::map<std::string, TTree*>& genMcTrees,
                           const std::map<std::string, TTree*>& recMcTrees,
                           const std::string& cuts_json,
                           const std::string& out_root_dir);

#endif // ACCEPTANCE_H