#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>  // For size_t
#include "write_csv.h"
#include "plot_unfolding.h"  // Assuming UnfoldingData is defined here

// Helper function to write a vector to CSV format
void write_vector_to_csv(std::ofstream& file, const std::vector<double>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        file << vec[i];
        if (i < vec.size() - 1) {
            file << ",";
        }
    }
}

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
            // Write the binning information for this phi bin
            file << data.bin_number << ","
                 << data.xB_min << "," << data.xB_max << "," << data.xB_avg << ","
                 << data.Q2_min << "," << data.Q2_max << "," << data.Q2_avg << ","
                 << data.t_min << "," << data.t_max << "," << data.t_avg << ","
                 << data.phi_min[i] << "," << data.phi_max[i] << ",";

            // Write raw yields, acceptance, and unfolded yields for each period
            for (size_t topo_idx = 0; topo_idx < 4; ++topo_idx) {
                size_t index = topo_idx * data.phi_min.size() + i;
                file << data.raw_yields_Fa18Inb[index] << "," 
                     << data.raw_yields_Fa18Out[index] << "," 
                     << data.raw_yields_Sp19Inb[index] << ",";
            }

            file << data.acceptance_Fa18Inb[i] << "," << data.unfolded_yields_Fa18Inb[i] << ","
                 << data.acceptance_Fa18Out[i] << "," << data.unfolded_yields_Fa18Out[i] << ","
                 << data.acceptance_Sp19Inb[i] << "," << data.unfolded_yields_Sp19Inb[i];

            // End the line after all periods' data is written
            file << std::endl;
        }
    }

    file.close();
    std::cout << "CSV file written to " << filename << std::endl;
}