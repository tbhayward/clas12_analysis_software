// write_csv.cpp

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>  // For size_t
#include <map>
#include "write_csv.h"
#include "unfolding_data.h"

void write_csv(const std::string& filename, const std::map<std::string, std::vector<UnfoldingData>>& topology_unfolding_data) {
    // Open the file
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Define period names
    std::vector<std::string> period_names = {"Fa18Inb", "Fa18Out", "Sp19Inb"};
    const int n_periods = 3;

    // Define topologies
    std::vector<std::string> topologies = {"FD_FD", "CD_FD", "CD_FT"};
    std::string combined_topology = "combined";

    // Write the header
    file << "Bin,xB_min,xB_max,xB_avg,Q2_min,Q2_max,Q2_avg,t_min,t_max,t_avg,phi_min,phi_max,phi_avg,";

    // Column headers for DVCS data (periods 0-2)
    for (int period = 0; period < n_periods; ++period) {
        std::string period_name = period_names[period];

        for (const auto& topology : topologies) {
            file << "ep->e'pgamma raw_yield_" << topology << "_" << period_name << ",";
        }

        file << "ep->e'pgamma raw_yield_combined_" << period_name << ",";
        file << "ep->e'pgamma acceptance_" << period_name << ",";
        file << "ep->e'pgamma unfolded_yield_" << period_name << ",";
        file << "ep->e'pgamma unfolded_yield_uncertainty_" << period_name << ",";  // Added column for uncertainty

        // Insert contamination_fraction, contamination_uncertainty, and signal_yield
        file << "ep->e'pgamma contamination_fraction_" << period_name << ",";
        file << "ep->e'pgamma contamination_uncertainty_" << period_name << ",";
        file << "ep->e'pgamma signal_yield_" << period_name << ",";
    }

    // Column headers for eppi0 data (periods 3-5)
    for (int period = 0; period < n_periods; ++period) {
        std::string period_name = period_names[period];

        for (const auto& topology : topologies) {
            file << "ep->e'ppi0 raw_yield_" << topology << "_" << period_name << ",";
        }

        file << "ep->e'ppi0 raw_yield_combined_" << period_name << ",";
        file << "ep->e'ppi0 acceptance_" << period_name << ",";
        file << "ep->e'ppi0 unfolded_yield_" << period_name << ",";
        file << "ep->e'ppi0 unfolded_yield_uncertainty_" << period_name;  // Added column for uncertainty

        if (period < n_periods - 1) {
            file << ",";
        }
    }

    file << std::endl;

    // Get the combined data
    const auto& combined_data = topology_unfolding_data.at(combined_topology);

    // Iterate over the data
    for (size_t group_idx = 0; group_idx < combined_data.size(); ++group_idx) {
        const auto& data = combined_data[group_idx];

        for (size_t phi_idx = 0; phi_idx < data.phi_min.size(); ++phi_idx) {
            // Write binning information
            file << data.bin_number << ","
                 << data.xB_min << "," << data.xB_max << "," << data.xB_avg << ","
                 << data.Q2_min << "," << data.Q2_max << "," << data.Q2_avg << ","
                 << data.t_min << "," << data.t_max << "," << data.t_avg << ","
                 << data.phi_min[phi_idx] << "," << data.phi_max[phi_idx] << ","
                 << data.phi_avg[phi_idx] << ",";

            // DVCS data (periods 0-2)
            for (int period = 0; period < n_periods; ++period) {
                // For each topology
                for (const auto& topology : topologies) {
                    const auto& topology_data = topology_unfolding_data.at(topology);
                    int raw_yield = topology_data[group_idx].raw_yields.at(topology)[period][phi_idx];
                    file << raw_yield << ",";
                }

                // Combined raw yield
                int combined_raw_yield = data.combined_raw_yields[period][phi_idx];
                file << combined_raw_yield << ",";

                // Acceptance, unfolded yield, and unfolded yield uncertainty
                double acceptance = data.acceptance[period][phi_idx];
                double unfolded_yield = data.unfolded_yields[period][phi_idx];
                double unfolded_yield_uncertainty = data.unfolded_yield_uncertainty[period][phi_idx];  // Retrieve uncertainty
                file << acceptance << "," << unfolded_yield << "," << unfolded_yield_uncertainty << ",";  // Write uncertainty

                // Contamination fraction, contamination uncertainty, and signal yield
                double contamination_fraction = data.contamination_fraction[period][phi_idx];
                double contamination_uncertainty = data.contamination_uncertainty[period][phi_idx];
                double signal_yield = data.signal_yield[period][phi_idx];
                file << contamination_fraction << "," << contamination_uncertainty << "," << signal_yield << ",";
            }

            // eppi0 data (periods 3-5)
            for (int period = 0; period < n_periods; ++period) {
                int eppi0_period = period + n_periods;  // Adjust index for eppi0 periods

                // For each topology
                for (const auto& topology : topologies) {
                    const auto& topology_data = topology_unfolding_data.at(topology);
                    int raw_yield = topology_data[group_idx].raw_yields.at(topology)[eppi0_period][phi_idx];
                    file << raw_yield << ",";
                }

                // Combined raw yield
                int combined_raw_yield = data.combined_raw_yields[eppi0_period][phi_idx];
                file << combined_raw_yield << ",";

                // Acceptance, unfolded yield, and unfolded yield uncertainty
                double acceptance = data.acceptance[eppi0_period][phi_idx];
                double unfolded_yield = data.unfolded_yields[eppi0_period][phi_idx];
                double unfolded_yield_uncertainty = data.unfolded_yield_uncertainty[eppi0_period][phi_idx];  // Retrieve uncertainty
                file << acceptance << "," << unfolded_yield << "," << unfolded_yield_uncertainty;  // Write uncertainty

                if (period < n_periods - 1) {
                    file << ",";
                }
            }

            file << std::endl;
        }
    }

    file.close();
    std::cout << "CSV file written to " << filename << std::endl;
}