#include <fstream>
#include <iostream>
#include "write_csv.h"

void write_csv(const std::string& filename, const std::vector<UnfoldingData>& unfolding_data) {
    std::ofstream csv_file(filename);

    // Check if the file opened successfully
    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write the CSV header
    csv_file << "bin_number,xB_min,xB_max,xB_avg,Q2_min,Q2_max,Q2_avg,t_min,t_max,t_avg,"
             << "phi_min,phi_max,"
             << "raw_yield_FD_FD_1,raw_yield_CD_FD_1,raw_yield_CD_FT_1,raw_yield_combined_1,"
             << "acceptance_1,unfolded_yield_1,"
             << "raw_yield_FD_FD_2,raw_yield_CD_FD_2,raw_yield_CD_FT_2,raw_yield_combined_2,"
             << "acceptance_2,unfolded_yield_2,"
             << "raw_yield_FD_FD_3,raw_yield_CD_FD_3,raw_yield_CD_FT_3,raw_yield_combined_3,"
             << "acceptance_3,unfolded_yield_3\n";

    // Loop over the unfolding data and write each entry to the CSV file
    for (const auto& data : unfolding_data) {
        for (size_t phi_idx = 0; phi_idx < data.phi_min.size(); ++phi_idx) {
            csv_file << data.bin_number << ","
                     << data.xB_min << "," << data.xB_max << "," << data.xB_avg << ","
                     << data.Q2_min << "," << data.Q2_max << "," << data.Q2_avg << ","
                     << data.t_min << "," << data.t_max << "," << data.t_avg << ","
                     << data.phi_min[phi_idx] << "," << data.phi_max[phi_idx] << ",";

            // For each run period, add the raw yields, acceptance, and unfolded yield
            for (int period_idx = 0; period_idx < 3; ++period_idx) {
                csv_file << data.raw_yields[period_idx][0] << ","  // raw_yield_FD_FD
                         << data.raw_yields[period_idx][1] << ","  // raw_yield_CD_FD
                         << data.raw_yields[period_idx][2] << ","  // raw_yield_CD_FT
                         << data.raw_yields[period_idx][3] << ","  // raw_yield_combined
                         << data.acceptance[period_idx][0] << ","  // acceptance
                         << data.unfolded_yields[period_idx][0] << ",";  // unfolded_yield
            }
            csv_file << "\n";  // End of line for this phi bin
        }
    }

    // Close the file after writing
    csv_file.close();
    std::cout << "CSV file saved successfully: " << filename << std::endl;
}