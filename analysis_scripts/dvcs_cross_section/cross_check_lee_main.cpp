// cross_check_lee_main.cpp
#include "load_csv.h"

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static void make_output_dirs(const std::string& base) {
    std::vector<std::string> sub = {
        "raw_yield","pi0_contamination","signal_yield",
        "rad_corrections","acceptance","unfolding",
        "bin_volume","cross_sections"
    };
    fs::create_directories(base);
    for (const auto& s : sub) fs::create_directories(fs::path(base) / s);
}

int main() {
    // Hard-coded CSV locations
    const std::string all_csv = "imports/all_bin_v3.csv";
    const std::string acc_csv = "imports/full_acc.csv";
    const std::string out_dir = "output/cross_check/lee";

    make_output_dirs(out_dir);

    std::cout << "[lee] Loading CSVs from imports/ ...\n";
    LeeData data = load_lee_csvs(all_csv, acc_csv, /*verbose=*/true);
    std::cout << "[lee] Loaded rows: " << data.rows.size() << "\n";
    std::cout << "[lee] First pass complete. Add plotting next.\n";
    return 0;
}