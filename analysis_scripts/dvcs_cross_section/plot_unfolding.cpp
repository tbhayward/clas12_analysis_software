// plot_unfolding.cpp

#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "bin_boundaries.h"
#include "kinematic_cuts.h"
#include "bin_helpers.h"
#include "unfolding_data.h"

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Helper function to check if phi is within bin, accounting for wrap-around
static bool phi_in_bin(double phi_deg, double phi_low, double phi_high) {
    if (phi_low <= phi_high) {
        return phi_deg >= phi_low && phi_deg <= phi_high;
    } else {
        // Wrap-around bin (e.g., from 330 to 30 degrees)
        return phi_deg >= phi_low || phi_deg <= phi_high;
    }
}

// Function implementation
std::vector<UnfoldingData> plot_unfolding(const std::string& output_dir,
                                          int xB_bin,
                                          const std::vector<BinBoundary>& bin_boundaries,
                                          std::vector<TTreeReader>& data_readers,
                                          std::vector<TTreeReader>& mc_gen_readers,
                                          std::vector<TTreeReader>& mc_rec_readers,
                                          std::vector<TTreeReader>& eppi0_readers,
                                          std::vector<TTreeReader>& mc_gen_aaogen_readers,
                                          std::vector<TTreeReader>& mc_rec_aaogen_readers) {
    // Vector to hold all unfolding data
    std::vector<UnfoldingData> all_unfolding_data;

    // Precompute relevant bins for the given xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);
    if (relevant_bins.empty()) {
        std::cerr << "Error: No relevant bins found for xB_bin: " << xB_bin << std::endl;
        return all_unfolding_data;
    }

    // Group the bins by (Q2_low, t_low)
    std::map<std::pair<double, double>, std::vector<int>> bin_groups;
    for (int idx : relevant_bins) {
        const auto& bin = bin_boundaries[idx];
        auto key = std::make_pair(bin.Q2_low, bin.t_low);
        bin_groups[key].push_back(idx);
    }

    if (bin_groups.empty()) {
        std::cerr << "Error: No bin groups created. Exiting function." << std::endl;
        return all_unfolding_data;
    }

    // Names of the periods
    std::vector<std::string> period_names = {"Fa18Inb", "Fa18Out", "Sp19Inb"};
    const int n_periods = 3; // Number of periods for data and MC

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    // Determine number of columns and rows for the canvas
    int n_bins = bin_groups.size();
    int n_columns = std::ceil(std::sqrt(n_bins));
    int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

    gStyle->SetOptStat(0);

    // Initialize UnfoldingData structures and histograms
    for (const auto& group : bin_groups) {
        const auto& idx_list = group.second; // List of indices in bin_boundaries

        // Create an UnfoldingData instance
        UnfoldingData unfolding_data;

        // Initialize bin ranges and averages
        const auto& bin_example = bin_boundaries[idx_list[0]];
        unfolding_data.bin_number = xB_bin;
        unfolding_data.xB_min = bin_example.xB_low;
        unfolding_data.xB_max = bin_example.xB_high;
        unfolding_data.xB_avg = bin_example.xB_avg;
        unfolding_data.Q2_min = bin_example.Q2_low;
        unfolding_data.Q2_max = bin_example.Q2_high;
        unfolding_data.Q2_avg = bin_example.Q2_avg;
        unfolding_data.t_min = bin_example.t_low;
        unfolding_data.t_max = bin_example.t_high;
        unfolding_data.t_avg = bin_example.t_avg;

        // Collect phi bins
        for (int idx : idx_list) {
            const auto& bin = bin_boundaries[idx];
            unfolding_data.phi_min.push_back(bin.phi_low);
            unfolding_data.phi_max.push_back(bin.phi_high);
        }

        size_t n_phi_bins = unfolding_data.phi_min.size();
        size_t total_periods = 6; // 0-2: DVCS, 3-5: eppi0

        // Initialize counts and acceptance vectors
        unfolding_data.raw_yields.resize(total_periods, std::vector<int>(n_phi_bins, 0));
        unfolding_data.acceptance.resize(total_periods, std::vector<double>(n_phi_bins, 0.0));
        unfolding_data.unfolded_yields.resize(total_periods, std::vector<double>(n_phi_bins, 0.0));
        // No changes to UnfoldingData struct; mc_gen_counts and mc_rec_counts will be local variables

        // Initialize contamination_ratio and contamination_error
        unfolding_data.contamination_ratio.resize(total_periods, std::vector<double>(n_phi_bins, 0.0));
        unfolding_data.contamination_error.resize(total_periods, std::vector<double>(n_phi_bins, 0.0));

        // Store the unfolding_data instance
        all_unfolding_data.push_back(unfolding_data);
    }

    // Initialize local variables to store mc_gen_counts and mc_rec_counts
    size_t total_periods = 6; // 0-2: DVCS, 3-5: eppi0
    size_t n_groups = all_unfolding_data.size();
    // mc_gen_counts[period][group_idx][phi_idx]
    std::vector<std::vector<std::vector<double>>> mc_gen_counts(total_periods, std::vector<std::vector<double>>(n_groups));
    std::vector<std::vector<std::vector<double>>> mc_rec_counts(total_periods, std::vector<std::vector<double>>(n_groups));

    // Initialize mc_gen_counts and mc_rec_counts vectors
    for (size_t period = 0; period < total_periods; ++period) {
        for (size_t group_idx = 0; group_idx < n_groups; ++group_idx) {
            size_t n_phi_bins = all_unfolding_data[group_idx].phi_min.size();
            mc_gen_counts[period][group_idx].resize(n_phi_bins, 0.0);
            mc_rec_counts[period][group_idx].resize(n_phi_bins, 0.0);
        }
    }

    // Process data readers (periods 0-2)
    for (int period = 0; period < n_periods; ++period) {
        TTreeReader& data_reader = data_readers[period];
        data_reader.Restart();

        // Set theta_variable_name based on period
        std::string theta_variable_name = "theta_gamma_gamma"; // For DVCS data

        // Initialize TTreeReaderValues
        TTreeReaderValue<double> phi_data(data_reader, "phi2");
        TTreeReaderValue<double> xB_data(data_reader, "x");
        TTreeReaderValue<double> Q2_data(data_reader, "Q2");
        TTreeReaderValue<double> t1_data(data_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_data(data_reader, "open_angle_ep2");
        TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");
        TTreeReaderValue<double> theta_neutral_neutral_data(data_reader, theta_variable_name.c_str());

        while (data_reader.Next()) {
            double phi_deg = *phi_data * RAD_TO_DEG;
            phi_deg = std::fmod(phi_deg + 360.0, 360.0); // Ensure phi in [0, 360)
            double xB_value = *xB_data;
            double Q2_value = *Q2_data;
            double t_abs = std::abs(*t1_data);

            // Apply kinematic cuts
            if (!apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data))
                continue;

            // Check if xB is within the xB_bin range
            const auto& bin_example = bin_boundaries[relevant_bins[0]];
            if (xB_value < bin_example.xB_low || xB_value > bin_example.xB_high)
                continue;

            // Find the bin group
            for (size_t group_idx = 0; group_idx < all_unfolding_data.size(); ++group_idx) {
                UnfoldingData& unfolding_data = all_unfolding_data[group_idx];

                if (Q2_value >= unfolding_data.Q2_min && Q2_value <= unfolding_data.Q2_max &&
                    t_abs >= unfolding_data.t_min && t_abs <= unfolding_data.t_max) {

                    // Find the phi bin within this group
                    for (size_t phi_idx = 0; phi_idx < unfolding_data.phi_min.size(); ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            // Increment raw_yields for this period
                            unfolding_data.raw_yields[period][phi_idx] += 1;
                            break; // Exit phi bin loop
                        }
                    }
                    break; // Exit bin group loop
                }
            }
        }
    }

    // Process eppi0 readers (periods 3-5)
    for (int period = 0; period < n_periods; ++period) {
        TTreeReader& eppi0_reader = eppi0_readers[period];
        eppi0_reader.Restart();

        // Set theta_variable_name based on period
        std::string theta_variable_name = "theta_pi0_pi0"; // For eppi0 data

        // Initialize TTreeReaderValues
        TTreeReaderValue<double> phi_data(eppi0_reader, "phi2");
        TTreeReaderValue<double> xB_data(eppi0_reader, "x");
        TTreeReaderValue<double> Q2_data(eppi0_reader, "Q2");
        TTreeReaderValue<double> t1_data(eppi0_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_data(eppi0_reader, "open_angle_ep2");
        TTreeReaderValue<double> Emiss2_data(eppi0_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_data(eppi0_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_data(eppi0_reader, "pTmiss");
        TTreeReaderValue<double> theta_neutral_neutral_data(eppi0_reader, theta_variable_name.c_str());

        while (eppi0_reader.Next()) {
            double phi_deg = *phi_data * RAD_TO_DEG;
            phi_deg = std::fmod(phi_deg + 360.0, 360.0); // Ensure phi in [0, 360)
            double xB_value = *xB_data;
            double Q2_value = *Q2_data;
            double t_abs = std::abs(*t1_data);

            // Apply kinematic cuts
            if (!apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data))
                continue;

            // Check if xB is within the xB_bin range
            const auto& bin_example = bin_boundaries[relevant_bins[0]];
            if (xB_value < bin_example.xB_low || xB_value > bin_example.xB_high)
                continue;

            // Find the bin group
            for (size_t group_idx = 0; group_idx < all_unfolding_data.size(); ++group_idx) {
                UnfoldingData& unfolding_data = all_unfolding_data[group_idx];

                if (Q2_value >= unfolding_data.Q2_min && Q2_value <= unfolding_data.Q2_max &&
                    t_abs >= unfolding_data.t_min && t_abs <= unfolding_data.t_max) {

                    // Find the phi bin within this group
                    for (size_t phi_idx = 0; phi_idx < unfolding_data.phi_min.size(); ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            // Increment raw_yields for this period (period + 3)
                            unfolding_data.raw_yields[period + 3][phi_idx] += 1;
                            break; // Exit phi bin loop
                        }
                    }
                    break; // Exit bin group loop
                }
            }
        }
    }

    // Process mc_gen_readers (periods 0-2)
    for (int period = 0; period < n_periods; ++period) {
        TTreeReader& mc_gen_reader = mc_gen_readers[period];
        mc_gen_reader.Restart();

        // Set theta_variable_name based on period
        std::string theta_variable_name = "theta_gamma_gamma"; // For DVCS MC data

        // Initialize TTreeReaderValues
        TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi2");
        TTreeReaderValue<double> xB_mc_gen(mc_gen_reader, "x");
        TTreeReaderValue<double> Q2_mc_gen(mc_gen_reader, "Q2");
        TTreeReaderValue<double> t1_mc_gen(mc_gen_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_mc_gen(mc_gen_reader, "open_angle_ep2");
        TTreeReaderValue<double> Emiss2_mc_gen(mc_gen_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_mc_gen(mc_gen_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_mc_gen(mc_gen_reader, "pTmiss");
        TTreeReaderValue<double> theta_neutral_neutral_mc_gen(mc_gen_reader, theta_variable_name.c_str());

        while (mc_gen_reader.Next()) {
            double phi_deg = *phi_mc_gen * RAD_TO_DEG;
            phi_deg = std::fmod(phi_deg + 360.0, 360.0); // Ensure phi in [0, 360)
            double xB_value = *xB_mc_gen;
            double Q2_value = *Q2_mc_gen;
            double t_abs = std::abs(*t1_mc_gen);

            // Apply kinematic cuts
            if (!apply_kinematic_cuts(*t1_mc_gen, *open_angle_ep2_mc_gen, *theta_neutral_neutral_mc_gen, *Emiss2_mc_gen, *Mx2_1_mc_gen, *pTmiss_mc_gen))
                continue;

            // Check if xB is within the xB_bin range
            const auto& bin_example = bin_boundaries[relevant_bins[0]];
            if (xB_value < bin_example.xB_low || xB_value > bin_example.xB_high)
                continue;

            // Find the bin group
            for (size_t group_idx = 0; group_idx < n_groups; ++group_idx) {
                UnfoldingData& unfolding_data = all_unfolding_data[group_idx];

                if (Q2_value >= unfolding_data.Q2_min && Q2_value <= unfolding_data.Q2_max &&
                    t_abs >= unfolding_data.t_min && t_abs <= unfolding_data.t_max) {

                    // Find the phi bin within this group
                    for (size_t phi_idx = 0; phi_idx < unfolding_data.phi_min.size(); ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            // Increment mc_gen_counts for this period and group
                            mc_gen_counts[period][group_idx][phi_idx] += 1.0;
                            break; // Exit phi bin loop
                        }
                    }
                    break; // Exit bin group loop
                }
            }
        }
    }

    // Process mc_rec_readers (periods 0-2)
    for (int period = 0; period < n_periods; ++period) {
        TTreeReader& mc_rec_reader = mc_rec_readers[period];
        mc_rec_reader.Restart();

        // Set theta_variable_name based on period
        std::string theta_variable_name = "theta_gamma_gamma"; // For DVCS MC data

        // Initialize TTreeReaderValues
        TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi2");
        TTreeReaderValue<double> xB_mc_rec(mc_rec_reader, "x");
        TTreeReaderValue<double> Q2_mc_rec(mc_rec_reader, "Q2");
        TTreeReaderValue<double> t1_mc_rec(mc_rec_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_mc_rec(mc_rec_reader, "open_angle_ep2");
        TTreeReaderValue<double> Emiss2_mc_rec(mc_rec_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_mc_rec(mc_rec_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_mc_rec(mc_rec_reader, "pTmiss");
        TTreeReaderValue<double> theta_neutral_neutral_mc_rec(mc_rec_reader, theta_variable_name.c_str());

        while (mc_rec_reader.Next()) {
            double phi_deg = *phi_mc_rec * RAD_TO_DEG;
            phi_deg = std::fmod(phi_deg + 360.0, 360.0); // Ensure phi in [0, 360)
            double xB_value = *xB_mc_rec;
            double Q2_value = *Q2_mc_rec;
            double t_abs = std::abs(*t1_mc_rec);

            // Apply kinematic cuts
            if (!apply_kinematic_cuts(*t1_mc_rec, *open_angle_ep2_mc_rec, *theta_neutral_neutral_mc_rec, *Emiss2_mc_rec, *Mx2_1_mc_rec, *pTmiss_mc_rec))
                continue;

            // Check if xB is within the xB_bin range
            const auto& bin_example = bin_boundaries[relevant_bins[0]];
            if (xB_value < bin_example.xB_low || xB_value > bin_example.xB_high)
                continue;

            // Find the bin group
            for (size_t group_idx = 0; group_idx < n_groups; ++group_idx) {
                UnfoldingData& unfolding_data = all_unfolding_data[group_idx];

                if (Q2_value >= unfolding_data.Q2_min && Q2_value <= unfolding_data.Q2_max &&
                    t_abs >= unfolding_data.t_min && t_abs <= unfolding_data.t_max) {

                    // Find the phi bin within this group
                    for (size_t phi_idx = 0; phi_idx < unfolding_data.phi_min.size(); ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            // Increment mc_rec_counts for this period and group
                            mc_rec_counts[period][group_idx][phi_idx] += 1.0;
                            break; // Exit phi bin loop
                        }
                    }
                    break; // Exit bin group loop
                }
            }
        }
    }

    // Process mc_gen_aaogen_readers and mc_rec_aaogen_readers for eppi0 (periods 3-5)
    for (int period = 0; period < n_periods; ++period) {
        int eppi0_period = period + 3; // Adjust period index for eppi0 data

        // mc_gen_aaogen_readers
        TTreeReader& mc_gen_reader = mc_gen_aaogen_readers[period];
        mc_gen_reader.Restart();

        // Set theta_variable_name based on period
        std::string theta_variable_name = "theta_pi0_pi0"; // For eppi0 MC data

        // Initialize TTreeReaderValues
        TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi2");
        TTreeReaderValue<double> xB_mc_gen(mc_gen_reader, "x");
        TTreeReaderValue<double> Q2_mc_gen(mc_gen_reader, "Q2");
        TTreeReaderValue<double> t1_mc_gen(mc_gen_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_mc_gen(mc_gen_reader, "open_angle_ep2");
        TTreeReaderValue<double> Emiss2_mc_gen(mc_gen_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_mc_gen(mc_gen_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_mc_gen(mc_gen_reader, "pTmiss");
        TTreeReaderValue<double> theta_neutral_neutral_mc_gen(mc_gen_reader, theta_variable_name.c_str());

        while (mc_gen_reader.Next()) {
            double phi_deg = *phi_mc_gen * RAD_TO_DEG;
            phi_deg = std::fmod(phi_deg + 360.0, 360.0); // Ensure phi in [0, 360)
            double xB_value = *xB_mc_gen;
            double Q2_value = *Q2_mc_gen;
            double t_abs = std::abs(*t1_mc_gen);

            // Apply kinematic cuts
            if (!apply_kinematic_cuts(*t1_mc_gen, *open_angle_ep2_mc_gen, *theta_neutral_neutral_mc_gen, *Emiss2_mc_gen, *Mx2_1_mc_gen, *pTmiss_mc_gen))
                continue;

            // Check if xB is within the xB_bin range
            const auto& bin_example = bin_boundaries[relevant_bins[0]];
            if (xB_value < bin_example.xB_low || xB_value > bin_example.xB_high)
                continue;

            // Find the bin group
            for (size_t group_idx = 0; group_idx < n_groups; ++group_idx) {
                UnfoldingData& unfolding_data = all_unfolding_data[group_idx];

                if (Q2_value >= unfolding_data.Q2_min && Q2_value <= unfolding_data.Q2_max &&
                    t_abs >= unfolding_data.t_min && t_abs <= unfolding_data.t_max) {

                    // Find the phi bin within this group
                    for (size_t phi_idx = 0; phi_idx < unfolding_data.phi_min.size(); ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            // Increment mc_gen_counts for this eppi0 period and group
                            mc_gen_counts[eppi0_period][group_idx][phi_idx] += 1.0;
                            break; // Exit phi bin loop
                        }
                    }
                    break; // Exit bin group loop
                }
            }
        }

        // mc_rec_aaogen_readers
        TTreeReader& mc_rec_reader = mc_rec_aaogen_readers[period];
        mc_rec_reader.Restart();

        // Set theta_variable_name based on period
        // Already set to "theta_pi0_pi0" above

        // Initialize TTreeReaderValues
        TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi2");
        TTreeReaderValue<double> xB_mc_rec(mc_rec_reader, "x");
        TTreeReaderValue<double> Q2_mc_rec(mc_rec_reader, "Q2");
        TTreeReaderValue<double> t1_mc_rec(mc_rec_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_mc_rec(mc_rec_reader, "open_angle_ep2");
        TTreeReaderValue<double> Emiss2_mc_rec(mc_rec_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_mc_rec(mc_rec_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_mc_rec(mc_rec_reader, "pTmiss");
        TTreeReaderValue<double> theta_neutral_neutral_mc_rec(mc_rec_reader, theta_variable_name.c_str());

        while (mc_rec_reader.Next()) {
            double phi_deg = *phi_mc_rec * RAD_TO_DEG;
            phi_deg = std::fmod(phi_deg + 360.0, 360.0); // Ensure phi in [0, 360)
            double xB_value = *xB_mc_rec;
            double Q2_value = *Q2_mc_rec;
            double t_abs = std::abs(*t1_mc_rec);

            // Apply kinematic cuts
            if (!apply_kinematic_cuts(*t1_mc_rec, *open_angle_ep2_mc_rec, *theta_neutral_neutral_mc_rec, *Emiss2_mc_rec, *Mx2_1_mc_rec, *pTmiss_mc_rec))
                continue;

            // Check if xB is within the xB_bin range
            const auto& bin_example = bin_boundaries[relevant_bins[0]];
            if (xB_value < bin_example.xB_low || xB_value > bin_example.xB_high)
                continue;

            // Find the bin group
            for (size_t group_idx = 0; group_idx < n_groups; ++group_idx) {
                UnfoldingData& unfolding_data = all_unfolding_data[group_idx];

                if (Q2_value >= unfolding_data.Q2_min && Q2_value <= unfolding_data.Q2_max &&
                    t_abs >= unfolding_data.t_min && t_abs <= unfolding_data.t_max) {

                    // Find the phi bin within this group
                    for (size_t phi_idx = 0; phi_idx < unfolding_data.phi_min.size(); ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            // Increment mc_rec_counts for this eppi0 period and group
                            mc_rec_counts[eppi0_period][group_idx][phi_idx] += 1.0;
                            break; // Exit phi bin loop
                        }
                    }
                    break; // Exit bin group loop
                }
            }
        }
    }

    // Compute acceptances and unfolded yields
    for (size_t group_idx = 0; group_idx < n_groups; ++group_idx) {
        UnfoldingData& unfolding_data = all_unfolding_data[group_idx];
        size_t n_phi_bins = unfolding_data.phi_min.size();

        for (size_t period = 0; period < total_periods; ++period) {
            for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                double mc_gen_count = mc_gen_counts[period][group_idx][phi_idx];
                double mc_rec_count = mc_rec_counts[period][group_idx][phi_idx];

                if (mc_gen_count > 0) {
                    double acceptance = mc_rec_count / mc_gen_count;
                    unfolding_data.acceptance[period][phi_idx] = acceptance;
                } else {
                    unfolding_data.acceptance[period][phi_idx] = 0.0;
                }

                // Unfolded yield
                double raw_yield = unfolding_data.raw_yields[period][phi_idx];
                if (unfolding_data.acceptance[period][phi_idx] > 0) {
                    unfolding_data.unfolded_yields[period][phi_idx] = raw_yield / unfolding_data.acceptance[period][phi_idx];
                } else {
                    unfolding_data.unfolded_yields[period][phi_idx] = 0.0;
                }
            }
        }
    }

    // Plotting code for DVCS data (periods 0-2)
    // Create canvases for each period and plot data vs unfolded yields
    for (int period = 0; period < n_periods; ++period) {
        TCanvas* canvas = new TCanvas(Form("c_unfolding_period_%d", period), "Unfolding Results", canvas_width, canvas_height);
        canvas->Divide(n_columns, n_rows);

        int pad_idx = 1;

        for (size_t group_idx = 0; group_idx < all_unfolding_data.size(); ++group_idx) {
            UnfoldingData& unfolding_data = all_unfolding_data[group_idx];
            size_t n_phi_bins = unfolding_data.phi_min.size();

            std::vector<double> phi_centers(n_phi_bins);
            std::vector<double> phi_widths(n_phi_bins);
            std::vector<double> raw_counts(n_phi_bins);
            std::vector<double> raw_errors(n_phi_bins);
            std::vector<double> unfolded_counts(n_phi_bins);
            std::vector<double> unfolded_errors(n_phi_bins);

            for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                phi_centers[phi_idx] = (unfolding_data.phi_min[phi_idx] + unfolding_data.phi_max[phi_idx]) / 2.0;
                phi_widths[phi_idx] = (unfolding_data.phi_max[phi_idx] - unfolding_data.phi_min[phi_idx]) / 2.0;

                raw_counts[phi_idx] = unfolding_data.raw_yields[period][phi_idx];
                raw_errors[phi_idx] = std::sqrt(unfolding_data.raw_yields[period][phi_idx]);

                unfolded_counts[phi_idx] = unfolding_data.unfolded_yields[period][phi_idx];
                unfolded_errors[phi_idx] = raw_errors[phi_idx] / (unfolding_data.acceptance[period][phi_idx] > 0 ? unfolding_data.acceptance[period][phi_idx] : 1);
            }

            // Create TGraphErrors
            TGraphErrors* graph_raw = new TGraphErrors(n_phi_bins, &phi_centers[0], &raw_counts[0], &phi_widths[0], &raw_errors[0]);
            TGraphErrors* graph_unfolded = new TGraphErrors(n_phi_bins, &phi_centers[0], &unfolded_counts[0], &phi_widths[0], &unfolded_errors[0]);

            // Set styles
            graph_raw->SetMarkerColor(kBlue);
            graph_raw->SetMarkerStyle(20);
            graph_raw->SetLineColor(kBlue);

            graph_unfolded->SetMarkerColor(kRed);
            graph_unfolded->SetMarkerStyle(21);
            graph_unfolded->SetLineColor(kRed);

            // Draw graphs
            canvas->cd(pad_idx);
            TPad* pad = (TPad*)gPad;
            pad->SetLeftMargin(0.15);
            pad->SetBottomMargin(0.15);

            // Create a frame histogram for axes
            double phi_min = 0.0;
            double phi_max = 360.0;
            double count_max = std::max(*std::max_element(raw_counts.begin(), raw_counts.end()),
                                        *std::max_element(unfolded_counts.begin(), unfolded_counts.end())) * 1.2;
            double count_min = 0.0;

            TH1F* frame = pad->DrawFrame(phi_min, count_min, phi_max, count_max);
            frame->GetXaxis()->SetTitle("#phi [deg]");
            frame->GetYaxis()->SetTitle("Counts");

            graph_raw->Draw("P SAME");
            graph_unfolded->Draw("P SAME");

            // Add legend
            TLegend* legend = new TLegend(0.55, 0.65, 0.9, 0.85);
            legend->AddEntry(graph_raw, "Raw Data", "lep");
            legend->AddEntry(graph_unfolded, "Unfolded Data", "lep");
            legend->SetTextSize(0.04);
            legend->Draw();

            // Add title with averages
            std::string title = Form("DVCS, %s, x_{B}=%.3f, Q^{2}=%.3f, -t=%.3f",
                                     period_names[period].c_str(),
                                     unfolding_data.xB_avg,
                                     unfolding_data.Q2_avg,
                                     unfolding_data.t_avg);
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5

            ++pad_idx;
        }

        // Save canvas
        std::string filename = output_dir + "/unfolding_results_DVCS_" + period_names[period] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        canvas->SaveAs(filename.c_str());

        delete canvas;
    }

    // Create plots for eppi0 data (periods 3-5)
    for (int period = 0; period < n_periods; ++period) {
        int eppi0_period = period + 3; // Adjust the period index for eppi0 data
        TCanvas* canvas = new TCanvas(Form("c_unfolding_eppi0_period_%d", eppi0_period), "Unfolding Results eppi0", canvas_width, canvas_height);
        canvas->Divide(n_columns, n_rows);

        int pad_idx = 1;

        for (size_t group_idx = 0; group_idx < all_unfolding_data.size(); ++group_idx) {
            UnfoldingData& unfolding_data = all_unfolding_data[group_idx];
            size_t n_phi_bins = unfolding_data.phi_min.size();

            std::vector<double> phi_centers(n_phi_bins);
            std::vector<double> phi_widths(n_phi_bins);
            std::vector<double> raw_counts(n_phi_bins);
            std::vector<double> raw_errors(n_phi_bins);
            std::vector<double> unfolded_counts(n_phi_bins);
            std::vector<double> unfolded_errors(n_phi_bins);

            for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                phi_centers[phi_idx] = (unfolding_data.phi_min[phi_idx] + unfolding_data.phi_max[phi_idx]) / 2.0;
                phi_widths[phi_idx] = (unfolding_data.phi_max[phi_idx] - unfolding_data.phi_min[phi_idx]) / 2.0;

                raw_counts[phi_idx] = unfolding_data.raw_yields[eppi0_period][phi_idx];
                raw_errors[phi_idx] = std::sqrt(unfolding_data.raw_yields[eppi0_period][phi_idx]);

                unfolded_counts[phi_idx] = unfolding_data.unfolded_yields[eppi0_period][phi_idx];
                unfolded_errors[phi_idx] = raw_errors[phi_idx] / (unfolding_data.acceptance[eppi0_period][phi_idx] > 0 ? unfolding_data.acceptance[eppi0_period][phi_idx] : 1);
            }

            // Create TGraphErrors
            TGraphErrors* graph_raw = new TGraphErrors(n_phi_bins, &phi_centers[0], &raw_counts[0], &phi_widths[0], &raw_errors[0]);
            TGraphErrors* graph_unfolded = new TGraphErrors(n_phi_bins, &phi_centers[0], &unfolded_counts[0], &phi_widths[0], &unfolded_errors[0]);

            // Set styles
            graph_raw->SetMarkerColor(kBlue);
            graph_raw->SetMarkerStyle(20);
            graph_raw->SetLineColor(kBlue);

            graph_unfolded->SetMarkerColor(kRed);
            graph_unfolded->SetMarkerStyle(21);
            graph_unfolded->SetLineColor(kRed);

            // Draw graphs
            canvas->cd(pad_idx);
            TPad* pad = (TPad*)gPad;
            pad->SetLeftMargin(0.15);
            pad->SetBottomMargin(0.15);

            // Create a frame histogram for axes
            double phi_min = 0.0;
            double phi_max = 360.0;
            double count_max = std::max(*std::max_element(raw_counts.begin(), raw_counts.end()),
                                        *std::max_element(unfolded_counts.begin(), unfolded_counts.end())) * 1.2;
            double count_min = 0.0;

            TH1F* frame = pad->DrawFrame(phi_min, count_min, phi_max, count_max);
            frame->GetXaxis()->SetTitle("#phi [deg]");
            frame->GetYaxis()->SetTitle("Counts");

            graph_raw->Draw("P SAME");
            graph_unfolded->Draw("P SAME");

            // Add legend
            TLegend* legend = new TLegend(0.55, 0.65, 0.9, 0.85);
            legend->AddEntry(graph_raw, "Raw Data eppi0", "lep");
            legend->AddEntry(graph_unfolded, "Unfolded Data eppi0", "lep");
            legend->SetTextSize(0.04);
            legend->Draw();

            // Add title with averages
            std::string title = Form("eppi0, %s, x_{B}=%.3f, Q^{2}=%.3f, -t=%.3f",
                                     period_names[period].c_str(),
                                     unfolding_data.xB_avg,
                                     unfolding_data.Q2_avg,
                                     unfolding_data.t_avg);
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5

            ++pad_idx;
        }

        // Save canvas
        std::string filename = output_dir + "/unfolding_results_eppi0_" + period_names[period] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        canvas->SaveAs(filename.c_str());

        delete canvas;
    }

    // Return the unfolding data
    return all_unfolding_data;
}