#include "write_csv.h"
#include <fstream>
#include <iostream>

void write_csv(const std::string& filename, const std::vector<UnfoldingData>& unfolding_data) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Write the header
    file << "Bin,xB_min,xB_max,xB_avg,Q2_min,Q2_max,Q2_avg,t_min,t_max,t_avg,"
         << "phi_min,phi_max,"
         << "raw_yield_FD_FD_Fa18Inb,raw_yield_CD_FD_Fa18Inb,raw_yield_CD_FT_Fa18Inb,raw_yield_combined_Fa18Inb,acceptance_Fa18Inb,unfolded_yield_Fa18Inb,"
         << "raw_yield_FD_FD_Fa18Out,raw_yield_CD_FD_Fa18Out,raw_yield_CD_FT_Fa18Out,raw_yield_combined_Fa18Out,acceptance_Fa18Out,unfolded_yield_Fa18Out,"
         << "raw_yield_FD_FD_Sp19Inb,raw_yield_CD_FD_Sp19Inb,raw_yield_CD_FT_Sp19Inb,raw_yield_combined_Sp19Inb,acceptance_Sp19Inb,unfolded_yield_Sp19Inb"
         << std::endl;

    // Loop through the unfolding data and write each entry to the file
    for (const auto& data : unfolding_data) {
        for (size_t i = 0; i < data.phi_min.size(); ++i) {
            // Write the binning information
            file << data.bin_number << ","
                 << data.xB_min << "," << data.xB_max << "," << data.xB_avg << ","
                 << data.Q2_min << "," << data.Q2_max << "," << data.Q2_avg << ","
                 << data.t_min << "," << data.t_max << "," << data.t_avg << ","
                 << data.phi_min[i] << "," << data.phi_max[i] << ",";

            // For each period, write raw yields, acceptance, and unfolded yields
            for (int period = 0; period < 3; ++period) {
                for (size_t topo_idx = 0; topo_idx < 4; ++topo_idx) {
                    // Directly access the raw yields for each topology and period
                    file << data.raw_yields[period][topo_idx][i] << ",";
                }

                // Write acceptance and unfolded yield for the current period
                file << data.acceptance[period][i] << ",";
                file << data.unfolded_yields[period][i];

                // Add a comma for the next period (except for the last one)
                if (period < 2) {
                    file << ",";
                }
            }

            file << std::endl;
        }
    }

    file.close();
    std::cout << "CSV file written to " << filename << std::endl;
}