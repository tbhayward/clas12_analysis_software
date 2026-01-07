#ifndef TOTAL_COUNTS_H
#define TOTAL_COUNTS_H

#include <map>
#include <string>

// Forward decls
class TTree;

// Update the CSV with raw yields (pos/neg/unpol) after exclusivity cuts,
// and write per-topology canvases into output/total_counts_plots/<Label>/.
// Periods and topologies are taken internally from periods.h and hard-coded
// to the three standard topologies (FD,FD), (CD,FD), (CD,FT).
bool update_total_counts_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,
    const std::string& out_root_dir,
    int max_workers);

#endif // TOTAL_COUNTS_H