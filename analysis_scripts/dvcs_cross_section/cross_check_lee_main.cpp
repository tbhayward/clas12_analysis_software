#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "load_csv.h"
#include "raw_signal_cross_check.h"
#include "pi0_contamination_cross_check.h"
#include "rad_correction_cross_check.h"
#include "unfolded_yields_cross_check.h"
#include "bin_centering_cross_check.h"
#include "bin_volume_cross_check.h"
#include "cross_section_cross_check.h"
#include "acceptance_cross_check.h"
#include "cross_section_hayward_cross_check.h"

namespace fs = std::filesystem;

static inline void mkoutdirs() {
    const std::string base = "output/cross_check/lee";
    const char* subdirs[] = {
        "raw_yield", "pi0_contamination", "signal_yield",
        "rad_corrections", "acceptance", "unfolding",
        "bin_volume", "bin_centering", "cross_sections",
        "cross_sections_hayward"
    };
    fs::create_directories(base);
    for (auto s : subdirs) fs::create_directories(fs::path(base) / s);
}

int main() {
    std::cout << "[lee] cross_check_lee starting...\n";

    // Hard-coded CSV locations
    const std::string lee_csv     = "imports/all_bin_v3.csv";
    const std::string hayward_csv = "output/csvs/dvcs_pass2_analysis.csv";

    mkoutdirs();

    // // Raw yield cross-check (reads directly from CSVs)
    // plot_raw_yield_cross_checks(
    //     lee_csv,                                  // Lee pass-1 CSV
    //     hayward_csv,                              // Hayward pass-2 CSV
    //     "output/cross_check/lee/raw_yield"        // output directory
    // );

    // // Pi0 contamination cross-check (reads directly from CSVs)
    // plot_pi0_contam_cross_checks(
    //     lee_csv,                                  // Lee pass-1 CSV
    //     hayward_csv,                              // Hayward pass-2 CSV
    //     "output/cross_check/lee/pi0_contamination" // output directory
    // );

    // // Radiative correction (Frad) cross-check
    // plot_rad_correction_cross_checks(
    //     lee_csv,
    //     hayward_csv,
    //     "output/cross_check/lee/rad_corrections"
    // );

    // // Acceptance cross-check (Lee vs Fa18 Inb/Out)
    // plot_acceptance_cross_checks(
    //     lee_csv,
    //     hayward_csv,
    //     "output/cross_check/lee/acceptance"
    // );

    // // Unfolded acceptance-corrected yield cross-check
    // plot_unfolded_yields_cross_checks(
    //     lee_csv,
    //     hayward_csv,
    //     "output/cross_check/lee/unfolding"
    // );

    // // Bin-centering correction cross-check (Fbin vs bin_volume)
    // plot_bin_centering_cross_checks(
    //     lee_csv,
    //     hayward_csv,
    //     "output/cross_check/lee/bin_centering"
    // );

    // // Bin volume cross-check
    // plot_bin_volume_cross_checks(
    //     lee_csv,
    //     hayward_csv,
    //     "output/cross_check/lee/bin_volume"
    // );

    // Cross section cross-check
    plot_cross_section_cross_checks(
        lee_csv,
        hayward_csv,
        "output/cross_check/lee/cross_sections"
    );

    // Cross section cross-check *between* Hayward run periods
    plot_cross_section_hayward_cross_checks(
        hayward_csv,
        "output/cross_check/lee/cross_sections_hayward"
    );

    std::cout << "[lee] cross_check_lee complete.\n";
    return 0;
}