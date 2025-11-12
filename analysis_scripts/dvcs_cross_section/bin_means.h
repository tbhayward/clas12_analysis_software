#ifndef BIN_MEANS_H
#define BIN_MEANS_H

#include <map>
#include <string>
class TTree;

// Update the grouped-averages columns for each period and combined groups in the CSV.
// Returns true on success, false on non-fatal I/O failures.
// (Schema problems will print a message and std::exit to avoid a half-written CSV.)
bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers);

#endif // BIN_MEANS_H