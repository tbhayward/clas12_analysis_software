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
            // Write the binning information
            file << data.bin_number << ","
                 << data.xB_min << "," << data.xB_max << "," << data.xB_avg << ","
                 << data.Q2_min << "," << data.Q2_max << "," << data.Q2_avg << ","
                 << data.t_min << "," << data.t_max << "," << data.t_avg << ","
                 << data.phi_min[i] << "," << data.phi_max[i] << ",";

            // Write data for each period (Fa18Inb, Fa18Out, Sp19Inb)
            for (int period = 0; period < 3; ++period) {
                // Write raw yields for each topology (FD,FD), (CD,FD), (CD,FT), and combined
                for (size_t topo_idx = 0; topo_idx < 4; ++topo_idx) {
                    size_t index = topo_idx * data.phi_min.size() + i;

                    // Debug info for size checks
                    std::cout << "Debug Info: period: " << period 
                              << ", raw_yields size: " << data.raw_yields[period].size() 
                              << ", expected size: " << 4 * data.phi_min.size() 
                              << ", access index: " << index 
                              << std::endl;

                    if (index < data.raw_yields[period].size()) {
                        file << data.raw_yields[period][index] << ",";
                    } else {
                        std::cerr << "Out of bounds access detected at index " << index << std::endl;
                        file << "0,";  // Writing a placeholder value if out-of-bounds
                    }
                }

                // Write acceptance and unfolded yields for the current phi bin
                file << data.acceptance[period][i] << ","; // Only write the current phi bin
                file << data.unfolded_yields[period][i];   // Only write the current phi bin

                if (period < 2) {
                    file << ",";  // Add a comma for the next period (except the last one)
                }
            }

            file << std::endl;
        }
    }

    file.close();
    std::cout << "CSV file written to " << filename << std::endl;
}