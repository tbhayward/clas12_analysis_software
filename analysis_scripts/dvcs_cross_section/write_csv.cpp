// write_csv.cpp

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>  // For size_t
#include <map>
#include <cmath>    // For M_PI, std::sqrt, std::isnan, std::isfinite
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
        file << "ep->e'pgamma unfolded_yield_uncertainty_" << period_name << ",";  // Unfolded yield uncertainty

        // Insert contamination_fraction, contamination_uncertainty, signal_yield, and signal_yield_uncertainty
        file << "ep->e'pgamma contamination_fraction_" << period_name << ",";
        file << "ep->e'pgamma contamination_uncertainty_" << period_name << ",";
        file << "ep->e'pgamma signal_yield_" << period_name << ",";
        file << "ep->e'pgamma signal_yield_uncertainty_" << period_name << ",";  // Added column for signal yield uncertainty
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
        file << "ep->e'ppi0 unfolded_yield_uncertainty_" << period_name;  // Unfolded yield uncertainty

        if (period < n_periods - 1) {
            file << ",";
        }
    }

    // Add new columns to the header
    file << ",bin_volume,fall_cross_section,fall_cross_section_stat_uncertainty,fall_cross_section_sys_uncertainty" << std::endl;

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
                double unfolded_yield_uncertainty = data.unfolded_yield_uncertainty[period][phi_idx];

                // Check for NaN or Inf and set to zero if found
                if (!std::isfinite(acceptance)) acceptance = 0.0;
                if (!std::isfinite(unfolded_yield)) unfolded_yield = 0.0;
                if (!std::isfinite(unfolded_yield_uncertainty)) unfolded_yield_uncertainty = 0.0;

                file << acceptance << "," << unfolded_yield << "," << unfolded_yield_uncertainty << ",";

                // Contamination fraction, contamination uncertainty, signal yield, and signal yield uncertainty
                double contamination_fraction = data.contamination_fraction[period][phi_idx];
                double contamination_uncertainty = data.contamination_uncertainty[period][phi_idx];
                double signal_yield = data.signal_yield[period][phi_idx];
                double signal_yield_uncertainty = data.signal_yield_uncertainty[period][phi_idx];

                // Check for NaN or Inf and set to zero if found
                if (!std::isfinite(contamination_fraction)) contamination_fraction = 0.0;
                if (!std::isfinite(contamination_uncertainty)) contamination_uncertainty = 0.0;
                if (!std::isfinite(signal_yield)) signal_yield = 0.0;
                if (!std::isfinite(signal_yield_uncertainty)) signal_yield_uncertainty = 0.0;

                file << contamination_fraction << "," << contamination_uncertainty << "," << signal_yield << "," << signal_yield_uncertainty << ",";
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
                double unfolded_yield_uncertainty = data.unfolded_yield_uncertainty[eppi0_period][phi_idx];

                // Check for NaN or Inf and set to zero if found
                if (!std::isfinite(acceptance)) acceptance = 0.0;
                if (!std::isfinite(unfolded_yield)) unfolded_yield = 0.0;
                if (!std::isfinite(unfolded_yield_uncertainty)) unfolded_yield_uncertainty = 0.0;

                file << acceptance << "," << unfolded_yield << "," << unfolded_yield_uncertainty;

                if (period < n_periods - 1) {
                    file << ",";
                }
            }

            // Calculate bin_volume
            double xB_min = data.xB_min;
            double xB_max = data.xB_max;
            double Q2_min = data.Q2_min;
            double Q2_max = data.Q2_max;
            double t_min = data.t_min;
            double t_max = data.t_max;
            double phi_min_deg = data.phi_min[phi_idx];
            double phi_max_deg = data.phi_max[phi_idx];
            double phi_min_rad = phi_min_deg * M_PI / 180.0;
            double phi_max_rad = phi_max_deg * M_PI / 180.0;

            // Adjust for possible wrap-around in phi
            double delta_phi_rad = phi_max_rad - phi_min_rad;
            if (delta_phi_rad < 0) {
                delta_phi_rad += 2 * M_PI;
            }

            // Check for NaN or Inf in bin boundaries
            if (!std::isfinite(delta_phi_rad)) delta_phi_rad = 0.0;
            if (!std::isfinite(xB_max) || !std::isfinite(xB_min)) xB_max = xB_min = 0.0;
            if (!std::isfinite(Q2_max) || !std::isfinite(Q2_min)) Q2_max = Q2_min = 0.0;
            if (!std::isfinite(t_max) || !std::isfinite(t_min)) t_max = t_min = 0.0;

            // Bin volume calculation
            double bin_volume = (xB_max - xB_min) * (Q2_max - Q2_min) * (t_max - t_min) * delta_phi_rad;

            // Check for NaN or negative bin volume
            if (!std::isfinite(bin_volume) || bin_volume <= 0.0) bin_volume = 0.0;

            // Get signal_yield for Fa18Inb (period 0) and Fa18Out (period 1)
            double signal_yield_Fa18Inb = data.signal_yield[0][phi_idx];
            double signal_yield_Fa18Out = data.signal_yield[1][phi_idx];

            // Check for NaN or Inf and set to zero if found
            if (!std::isfinite(signal_yield_Fa18Inb)) signal_yield_Fa18Inb = 0.0;
            if (!std::isfinite(signal_yield_Fa18Out)) signal_yield_Fa18Out = 0.0;

            double combined_signal_yield = signal_yield_Fa18Inb + signal_yield_Fa18Out;

            // Get uncertainties
            double uncertainty_Fa18Inb = data.signal_yield_uncertainty[0][phi_idx];
            double uncertainty_Fa18Out = data.signal_yield_uncertainty[1][phi_idx];

            // Check for NaN or Inf and set to zero if found
            if (!std::isfinite(uncertainty_Fa18Inb)) uncertainty_Fa18Inb = 0.0;
            if (!std::isfinite(uncertainty_Fa18Out)) uncertainty_Fa18Out = 0.0;

            double combined_uncertainty = std::sqrt(uncertainty_Fa18Inb * uncertainty_Fa18Inb + uncertainty_Fa18Out * uncertainty_Fa18Out);

            // Compute fall_cross_section
            double denominator = 8.213e7 * bin_volume; // fa18 integrated luminosity
            double fall_cross_section = 0.0;
            double fall_cross_section_stat_uncertainty = 0.0;

            if (denominator > 0.0) {
                fall_cross_section = combined_signal_yield / denominator;
                fall_cross_section_stat_uncertainty = combined_uncertainty / denominator;
            } else {
                // Handle division by zero
                fall_cross_section = 0.0;
                fall_cross_section_stat_uncertainty = 0.0;
            }

            // Compute fall_cross_section_sys_uncertainty (50% of cross section)
            double fall_cross_section_sys_uncertainty = 1.0 * std::abs(fall_cross_section);

            // Write these values to the file
            file << "," << bin_volume << "," << fall_cross_section << "," << fall_cross_section_stat_uncertainty << "," << fall_cross_section_sys_uncertainty << std::endl;
        }
    }

    file.close();
    std::cout << "CSV file written to " << filename << std::endl;
}