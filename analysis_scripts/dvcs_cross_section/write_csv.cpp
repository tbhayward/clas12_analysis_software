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

// Assuming you have a structure that collects yields for different run periods
struct YieldData {
    double xB_min, xB_max, xB_avg, Q2_min, Q2_max, Q2_avg, t_min, t_max, t_avg, phi_min, phi_max;
    double raw_yield_FD_FD_Fa18Inb, raw_yield_CD_FD_Fa18Inb, raw_yield_CD_FT_Fa18Inb, raw_yield_combined_Fa18Inb, acceptance_Fa18Inb, unfolded_yield_Fa18Inb;
    double raw_yield_FD_FD_Fa18Out, raw_yield_CD_FD_Fa18Out, raw_yield_CD_FT_Fa18Out, raw_yield_combined_Fa18Out, acceptance_Fa18Out, unfolded_yield_Fa18Out;
    double raw_yield_FD_FD_Sp19Inb, raw_yield_CD_FD_Sp19Inb, raw_yield_CD_FT_Sp19Inb, raw_yield_combined_Sp19Inb, acceptance_Sp19Inb, unfolded_yield_Sp19Inb;
};

// This function accumulates data for each unique bin and avoids writing the same bin multiple times
void write_csv(std::vector<YieldData> &yields, const std::string &filename) {
    // Open the output CSV file
    std::ofstream outfile(filename);

    // Write the header
    outfile << "Bin,xB_min,xB_max,xB_avg,Q2_min,Q2_max,Q2_avg,t_min,t_max,t_avg,phi_min,phi_max,"
            << "raw_yield_FD_FD_Fa18Inb,raw_yield_CD_FD_Fa18Inb,raw_yield_CD_FT_Fa18Inb,raw_yield_combined_Fa18Inb,acceptance_Fa18Inb,unfolded_yield_Fa18Inb,"
            << "raw_yield_FD_FD_Fa18Out,raw_yield_CD_FD_Fa18Out,raw_yield_CD_FT_Fa18Out,raw_yield_combined_Fa18Out,acceptance_Fa18Out,unfolded_yield_Fa18Out,"
            << "raw_yield_FD_FD_Sp19Inb,raw_yield_CD_FD_Sp19Inb,raw_yield_CD_FT_Sp19Inb,raw_yield_combined_Sp19Inb,acceptance_Sp19Inb,unfolded_yield_Sp19Inb\n";

    // Write the data rows
    for (const auto &data : yields) {
        outfile << data.xB_min << "," << data.xB_max << "," << data.xB_avg << ","
                << data.Q2_min << "," << data.Q2_max << "," << data.Q2_avg << ","
                << data.t_min << "," << data.t_max << "," << data.t_avg << ","
                << data.phi_min << "," << data.phi_max << ","
                // Fa18Inb
                << data.raw_yield_FD_FD_Fa18Inb << "," << data.raw_yield_CD_FD_Fa18Inb << "," << data.raw_yield_CD_FT_Fa18Inb << "," << data.raw_yield_combined_Fa18Inb << ","
                << data.acceptance_Fa18Inb << "," << data.unfolded_yield_Fa18Inb << ","
                // Fa18Out
                << data.raw_yield_FD_FD_Fa18Out << "," << data.raw_yield_CD_FD_Fa18Out << "," << data.raw_yield_CD_FT_Fa18Out << "," << data.raw_yield_combined_Fa18Out << ","
                << data.acceptance_Fa18Out << "," << data.unfolded_yield_Fa18Out << ","
                // Sp19Inb
                << data.raw_yield_FD_FD_Sp19Inb << "," << data.raw_yield_CD_FD_Sp19Inb << "," << data.raw_yield_CD_FT_Sp19Inb << "," << data.raw_yield_combined_Sp19Inb << ","
                << data.acceptance_Sp19Inb << "," << data.unfolded_yield_Sp19Inb << "\n";
    }

    // Close the file
    outfile.close();
}