#ifndef BIN_MEANS_H
#define BIN_MEANS_H

#include <TTree.h>
#include <map>
#include <string>

// Update the bin-mean columns in-place inside a CSV.
// - csv_path: usually "output/csvs/dvcs_pass2_analysis.csv"
// - dataTrees: map of canonical keys to TTree*, e.g. 
//     "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb",
//     "DVCS_Sp18_inb", "DVCS_Sp18_out"
// - max_workers: number of OpenMP threads to use (up to 5 suggested)
bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers);

#endif // BIN_MEANS_H