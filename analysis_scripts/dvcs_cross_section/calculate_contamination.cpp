#include "calculate_contamination.h"
#include "kinematic_cuts.h"
#include "bin_helpers.h"
#include "unfolding_data.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <map>

// Constants
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

void calculate_contamination(
    const std::string& base_output_dir,
    int xB_bin,
    const std::vector<BinBoundary>& bin_boundaries,
    std::vector<TTreeReader>& data_readers,
    std::vector<TTreeReader>& eppi0_readers,
    std::vector<TTreeReader>& mc_rec_aaogen_readers,
    std::vector<TTreeReader>& mc_rec_eppi0_bkg_readers,
    std::map<std::string, std::vector<UnfoldingData>>& bin_data
) {
    // Only process the "combined" topology
    if (bin_data.find("combined") == bin_data.end()) {
        std::cerr << "Error: 'combined' topology not found in bin_data." << std::endl;
        return;
    }

    auto& unfolding_data_vec = bin_data["combined"];

    const int n_periods = 3; // Periods 0-2

    // Loop over each UnfoldingData instance
    for (auto& unfolding_data : unfolding_data_vec) {
        size_t n_phi_bins = unfolding_data.phi_min.size();

        // Initialize contamination_fraction and signal_yield
        unfolding_data.contamination_fraction.resize(n_periods, std::vector<double>(n_phi_bins, 0.0));
        unfolding_data.signal_yield.resize(n_periods, std::vector<double>(n_phi_bins, 0.0));

        // For each period
        for (int period = 0; period < n_periods; ++period) {
            // Initialize counters
            std::vector<double> N_DVCS_data(n_phi_bins, 0.0);
            std::vector<double> N_eppi0_data(n_phi_bins, 0.0);
            std::vector<double> N_eppi0_sim(n_phi_bins, 0.0);
            std::vector<double> N_eppi0_misID_sim(n_phi_bins, 0.0);

            // **Process DVCS data**
            {
                TTreeReader& reader = data_readers[period];
                reader.Restart();

                // Readers for data
                TTreeReaderValue<double> phi_data(reader, "phi2");
                TTreeReaderValue<double> xB_data(reader, "x");
                TTreeReaderValue<double> Q2_data(reader, "Q2");
                TTreeReaderValue<double> t1_data(reader, "t1");
                TTreeReaderValue<double> open_angle_ep2_data(reader, "open_angle_ep2");
                TTreeReaderValue<double> Emiss2_data(reader, "Emiss2");
                TTreeReaderValue<double> Mx2_1_data(reader, "Mx2_1");
                TTreeReaderValue<double> pTmiss_data(reader, "pTmiss");
                TTreeReaderValue<double> theta_neutral_neutral_data(reader, "theta_gamma_gamma");

                while (reader.Next()) {
                    double phi_deg = *phi_data * RAD_TO_DEG;
                    phi_deg = std::fmod(phi_deg + 360.0, 360.0);

                    double xB_value = *xB_data;
                    double Q2_value = *Q2_data;
                    double t_abs = std::abs(*t1_data);

                    // Apply kinematic cuts
                    if (!apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data,
                                              *Emiss2_data, *Mx2_1_data, *pTmiss_data))
                        continue;

                    // Check if xB is within the xB_bin range
                    if (xB_value < unfolding_data.xB_min || xB_value > unfolding_data.xB_max)
                        continue;

                    // Check if Q2 and t are within the bin ranges
                    if (Q2_value < unfolding_data.Q2_min || Q2_value > unfolding_data.Q2_max ||
                        t_abs < unfolding_data.t_min || t_abs > unfolding_data.t_max)
                        continue;

                    // Find the phi bin
                    for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            N_DVCS_data[phi_idx] += 1.0;
                            break;
                        }
                    }
                }
            }

            // **Process eppi0 data**
            {
                TTreeReader& reader = eppi0_readers[period];
                reader.Restart();

                // Readers for eppi0 data
                TTreeReaderValue<double> phi_data(reader, "phi2");
                TTreeReaderValue<double> xB_data(reader, "x");
                TTreeReaderValue<double> Q2_data(reader, "Q2");
                TTreeReaderValue<double> t1_data(reader, "t1");
                TTreeReaderValue<double> open_angle_ep2_data(reader, "open_angle_ep2");
                TTreeReaderValue<double> Emiss2_data(reader, "Emiss2");
                TTreeReaderValue<double> Mx2_1_data(reader, "Mx2_1");
                TTreeReaderValue<double> pTmiss_data(reader, "pTmiss");
                TTreeReaderValue<double> theta_neutral_neutral_data(reader, "theta_pi0_pi0");

                while (reader.Next()) {
                    double phi_deg = *phi_data * RAD_TO_DEG;
                    phi_deg = std::fmod(phi_deg + 360.0, 360.0);

                    double xB_value = *xB_data;
                    double Q2_value = *Q2_data;
                    double t_abs = std::abs(*t1_data);

                    // Apply kinematic cuts
                    if (!apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data,
                                              *Emiss2_data, *Mx2_1_data, *pTmiss_data))
                        continue;

                    // Check if xB is within the xB_bin range
                    if (xB_value < unfolding_data.xB_min || xB_value > unfolding_data.xB_max)
                        continue;

                    // Check if Q2 and t are within the bin ranges
                    if (Q2_value < unfolding_data.Q2_min || Q2_value > unfolding_data.Q2_max ||
                        t_abs < unfolding_data.t_min || t_abs > unfolding_data.t_max)
                        continue;

                    // Find the phi bin
                    for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            N_eppi0_data[phi_idx] += 1.0;
                            break;
                        }
                    }
                }
            }

            // **Process eppi0 simulation**
            {
                TTreeReader& reader = mc_rec_aaogen_readers[period];
                reader.Restart();

                // Readers for eppi0 simulation
                TTreeReaderValue<double> phi_data(reader, "phi2");
                TTreeReaderValue<double> xB_data(reader, "x");
                TTreeReaderValue<double> Q2_data(reader, "Q2");
                TTreeReaderValue<double> t1_data(reader, "t1");
                TTreeReaderValue<double> open_angle_ep2_data(reader, "open_angle_ep2");
                TTreeReaderValue<double> Emiss2_data(reader, "Emiss2");
                TTreeReaderValue<double> Mx2_1_data(reader, "Mx2_1");
                TTreeReaderValue<double> pTmiss_data(reader, "pTmiss");
                TTreeReaderValue<double> theta_neutral_neutral_data(reader, "theta_pi0_pi0");

                while (reader.Next()) {
                    double phi_deg = *phi_data * RAD_TO_DEG;
                    phi_deg = std::fmod(phi_deg + 360.0, 360.0);

                    double xB_value = *xB_data;
                    double Q2_value = *Q2_data;
                    double t_abs = std::abs(*t1_data);

                    // Apply kinematic cuts
                    if (!apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data,
                                              *Emiss2_data, *Mx2_1_data, *pTmiss_data))
                        continue;

                    // Check if xB is within the xB_bin range
                    if (xB_value < unfolding_data.xB_min || xB_value > unfolding_data.xB_max)
                        continue;

                    // Check if Q2 and t are within the bin ranges
                    if (Q2_value < unfolding_data.Q2_min || Q2_value > unfolding_data.Q2_max ||
                        t_abs < unfolding_data.t_min || t_abs > unfolding_data.t_max)
                        continue;

                    // Find the phi bin
                    for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            N_eppi0_sim[phi_idx] += 1.0;
                            break;
                        }
                    }
                }
            }

            // **Process misidentified eppi0 in simulation**
            {
                TTreeReader& reader = mc_rec_eppi0_bkg_readers[period];
                reader.Restart();

                // Readers for misidentified eppi0 simulation
                TTreeReaderValue<double> phi_data(reader, "phi2");
                TTreeReaderValue<double> xB_data(reader, "x");
                TTreeReaderValue<double> Q2_data(reader, "Q2");
                TTreeReaderValue<double> t1_data(reader, "t1");
                TTreeReaderValue<double> open_angle_ep2_data(reader, "open_angle_ep2");
                TTreeReaderValue<double> Emiss2_data(reader, "Emiss2");
                TTreeReaderValue<double> Mx2_1_data(reader, "Mx2_1");
                TTreeReaderValue<double> pTmiss_data(reader, "pTmiss");
                TTreeReaderValue<double> theta_neutral_neutral_data(reader, "theta_gamma_gamma");

                while (reader.Next()) {
                    double phi_deg = *phi_data * RAD_TO_DEG;
                    phi_deg = std::fmod(phi_deg + 360.0, 360.0);

                    double xB_value = *xB_data;
                    double Q2_value = *Q2_data;
                    double t_abs = std::abs(*t1_data);

                    // Apply kinematic cuts
                    if (!apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data,
                                              *Emiss2_data, *Mx2_1_data, *pTmiss_data))
                        continue;

                    // Check if xB is within the xB_bin range
                    if (xB_value < unfolding_data.xB_min || xB_value > unfolding_data.xB_max)
                        continue;

                    // Check if Q2 and t are within the bin ranges
                    if (Q2_value < unfolding_data.Q2_min || Q2_value > unfolding_data.Q2_max ||
                        t_abs < unfolding_data.t_min || t_abs > unfolding_data.t_max)
                        continue;

                    // Find the phi bin
                    for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                        double phi_low = unfolding_data.phi_min[phi_idx];
                        double phi_high = unfolding_data.phi_max[phi_idx];

                        if (phi_in_bin(phi_deg, phi_low, phi_high)) {
                            N_eppi0_misID_sim[phi_idx] += 1.0;
                            break;
                        }
                    }
                }
            }

            // **Calculate contamination fractions and correct unfolded yields**
            for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                double N_DVCS = N_DVCS_data[phi_idx];
                double N_eppi0_data_count = N_eppi0_data[phi_idx];
                double N_eppi0_sim_count = N_eppi0_sim[phi_idx];
                double N_eppi0_misID_sim_count = N_eppi0_misID_sim[phi_idx];

                double contamination = 0.0;

                if (N_DVCS > 0.0 && N_eppi0_sim_count > 0.0) {
                    double ratio = N_eppi0_data_count / N_eppi0_sim_count;
                    contamination = (N_eppi0_misID_sim_count * ratio) / N_DVCS;
                    unfolding_data.contamination_fraction[period][phi_idx] = contamination;
                } else {
                    unfolding_data.contamination_fraction[period][phi_idx] = 0.0;
                }

                // **Correct unfolded yield to get signal yield**
                double unfolded_yield = unfolding_data.unfolded_yields[period][phi_idx];
                double signal_yield = unfolded_yield * (1.0 - contamination);
                unfolding_data.signal_yield[period][phi_idx] = signal_yield;
            }
        }
    }
}