#ifndef TOTAL_COUNTS_H
#define TOTAL_COUNTS_H

#include <map>
#include <string>
#include <vector>
#include <TTree.h>

// Update the master CSV with raw yields (pos/neg/unpol) per period AND per topology,
// and the combined groups (Fa18, Sp18, 10.6 GeV) per topology.
// Also makes per-topology canvases under output/total_counts_plots/<label>/<topo>/.
// Returns true on success.
bool update_total_counts_csv(
    const std::string& csv_path,
    const std::vector<std::string>& periods,                 // canonical tree keys (e.g. "DVCS_Fa18_inb")
    const std::vector<std::string>& topologies,              // {"(FD, FD)","(CD, FD)","(CD, FT)"}
    const std::map<std::string, TTree*>& dvcsDataTrees,      // keys are canonical tree keys
    const std::string& combined_cuts_json,                   // same combined cuts JSON used before
    const std::string& out_root_dir,                         // e.g. "output"
    int max_workers                                          // e.g. 5
);

#endif