#ifndef BIN_MEANS_H
#define BIN_MEANS_H

#include <map>
#include <string>
class TTree;

// Update grouped averages using canonical nominal Python-optimized production cuts.
// Returns true on success, false on non-fatal I/O failures.
// (Schema problems will print a message and std::exit to avoid a half-written CSV.)
struct BinMeansOptions {
    // Write compact, analysis-note-oriented diagnostics after the numerical
    // means have been filled. These are separate from the production CSV and
    // are intended only for documentation/QA.
    bool make_note_outputs = true;
    std::string note_output_dir = "output/bin_means/analysis_note";
};

bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers,
                          const BinMeansOptions& options = BinMeansOptions());

#endif // BIN_MEANS_H