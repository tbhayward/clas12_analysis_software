#ifndef TOTAL_COUNTS_H
#define TOTAL_COUNTS_H

#include <map>
#include <string>

class TTree;

bool update_total_counts_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,
    const std::string& global_cuts_json,
    const std::string& out_root_dir,
    int max_workers);

#endif // TOTAL_COUNTS_H