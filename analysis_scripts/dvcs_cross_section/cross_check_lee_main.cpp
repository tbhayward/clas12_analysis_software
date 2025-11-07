#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "load_csv.h"
#include "raw_signal_cross_check.h"
#include "pi0_contamination_cross_check.h"


namespace fs = std::filesystem;

static inline void mkoutdirs() {
    const std::string base = "output/cross_check/lee";
    const char* subdirs[] = {
        "raw_yield", "pi0_contamination", "signal_yield",
        "rad_corrections", "acceptance", "unfolding",
        "bin_volume", "cross_sections"
    };
    fs::create_directories(base);
    for (auto s : subdirs) fs::create_directories(fs::path(base) / s);
}

int main() {
    std::cout << "[lee] cross_check_lee starting...\n";

    // Hard-coded locations as requested
    const std::string all_bin_v3 = "imports/all_bin_v3.csv";
    const std::string full_acc   = "imports/full_acc.csv";

    mkoutdirs();

    int acc_ok = 0, acc_bad = 0;
    auto rows = load_lee_csvs(all_bin_v3, full_acc, acc_ok, acc_bad);

    std::cout << "[lee] Loaded rows: " << rows.size() << "\n";
    std::cout << "[lee] First pass complete.\n";

    // // NEW: make the raw-yield comparison/ratio canvases
    // plot_raw_yield_cross_checks(rows,
    //     "output/jsons/total_counts.json",
    //     "output/cross_check/lee/raw_yield");

    // NEW: pi0 contamination comparisons
    plot_pi0_contam_cross_checks(rows,
        "output/jsons/pi0_contamination_combined.json",
        "output/cross_check/lee/pi0_contamination");

    return 0;
}