// calculate_contamination.cpp

#include "calculate_contamination.h"
#include "kinematic_cuts.h"
#include "bin_helpers.h"
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Constants
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Function implementation
void calculate_contamination(const std::string& base_output_dir,
                             int xB_bin,
                             const std::vector<BinBoundary>& bin_boundaries,
                             std::vector<TTreeReader>& data_readers,
                             std::vector<TTreeReader>& eppi0_readers,
                             std::vector<TTreeReader>& mc_rec_aaogen_readers,
                             std::vector<TTreeReader>& mc_rec_eppi0_bkg_readers,
                             std::vector<UnfoldingData>& unfolding_data) {
    // Set style to remove stat boxes
    gStyle->SetOptStat(0);

    // Number of periods
    const int n_periods = 3;

    // Names of the periods
    std::vector<std::string> period_names = {"Fa18Inb", "Fa18Out", "Sp19Inb"};

    // Number of phi bins
    const int n_phi_bins = 24;

    // Determine the relevant Q²-t bins for the given xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);
    int n_Q2t_bins = relevant_bins.size();

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    // Determine number of rows and columns for the canvases
    int n_columns = (xB_bin == 3 || xB_bin == 4) ? 5 : std::sqrt(next_perfect_square(n_Q2t_bins));
    int n_rows = std::ceil(static_cast<double>(n_Q2t_bins) / n_columns);

    // Initialize histograms for each period and bin
    std::vector<std::vector<TH1D*>> h_DVCS_data(n_periods, std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_eppi0_data(n_periods, std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_eppi0_sim(n_periods, std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_eppi0_misID_sim(n_periods, std::vector<TH1D*>(n_Q2t_bins));

    // Create histograms
    for (int period = 0; period < n_periods; ++period) {
        for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
            int bin_idx = relevant_bins[idx];
            const auto& bin = bin_boundaries[bin_idx];

            // Create title string with the period name and bin information
            std::string title = Form("%s, <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f",
                                     period_names[period].c_str(), bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));

            // Histograms for DVCS data
            h_DVCS_data[period][idx] = new TH1D(Form("h_DVCS_data_%d_%zu", period, idx),
                                                title.c_str(), n_phi_bins, 0, 360);

            // Histograms for eppi0 data
            h_eppi0_data[period][idx] = new TH1D(Form("h_eppi0_data_%d_%zu", period, idx),
                                                 title.c_str(), n_phi_bins, 0, 360);

            // Histograms for eppi0 simulation
            h_eppi0_sim[period][idx] = new TH1D(Form("h_eppi0_sim_%d_%zu", period, idx),
                                                title.c_str(), n_phi_bins, 0, 360);

            // Histograms for misidentified eppi0 in simulation
            h_eppi0_misID_sim[period][idx] = new TH1D(Form("h_eppi0_misID_sim_%d_%zu", period, idx),
                                                      title.c_str(), n_phi_bins, 0, 360);

            // Set axis labels and format for histograms
            std::vector<TH1D*> histograms = {h_DVCS_data[period][idx], h_eppi0_data[period][idx],
                                             h_eppi0_sim[period][idx], h_eppi0_misID_sim[period][idx]};

            for (auto& hist : histograms) {
                hist->GetXaxis()->SetTitle("#phi [deg]");
                hist->GetYaxis()->SetTitle("Counts");
                hist->GetXaxis()->SetLabelSize(0.05);
                hist->GetYaxis()->SetLabelSize(0.05);
                hist->GetXaxis()->SetTitleSize(0.06);
                hist->GetYaxis()->SetTitleSize(0.06);
                hist->SetMarkerStyle(20);
                hist->SetMarkerSize(1.2);
                // Removed SetDrawOption; specify draw option when calling Draw()
            }
        }
    }

    // Loop over each period to fill histograms
    for (int period = 0; period < n_periods; ++period) {
        // Get TTreeReaders for the current period
        TTreeReader& data_reader = data_readers[period];
        TTreeReader& eppi0_data_reader = eppi0_readers[period];
        TTreeReader& eppi0_sim_reader = mc_rec_aaogen_readers[period];
        TTreeReader& eppi0_misID_sim_reader = mc_rec_eppi0_bkg_readers[period];

        // Define TTreeReaderValues for variables
        // For data_reader
        TTreeReaderValue<double> phi_data(data_reader, "phi");
        TTreeReaderValue<double> xB_data(data_reader, "x");
        TTreeReaderValue<double> Q2_data(data_reader, "Q2");
        TTreeReaderValue<double> t1_data(data_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_data(data_reader, "open_angle_ep2");
        TTreeReaderValue<double> theta_neutral_neutral_data(data_reader, "theta_gamma_gamma");
        TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");

        // For eppi0_data_reader
        TTreeReaderValue<double> phi_eppi0(eppi0_data_reader, "phi");
        TTreeReaderValue<double> xB_eppi0(eppi0_data_reader, "x");
        TTreeReaderValue<double> Q2_eppi0(eppi0_data_reader, "Q2");
        TTreeReaderValue<double> t1_eppi0(eppi0_data_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_eppi0(eppi0_data_reader, "open_angle_ep2");
        TTreeReaderValue<double> theta_neutral_neutral_eppi0(eppi0_data_reader, "theta_pi0_pi0");
        TTreeReaderValue<double> Emiss2_eppi0(eppi0_data_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_eppi0(eppi0_data_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_eppi0(eppi0_data_reader, "pTmiss");

        // For eppi0_sim_reader
        TTreeReaderValue<double> phi_eppi0_sim(eppi0_sim_reader, "phi");
        TTreeReaderValue<double> xB_eppi0_sim(eppi0_sim_reader, "x");
        TTreeReaderValue<double> Q2_eppi0_sim(eppi0_sim_reader, "Q2");
        TTreeReaderValue<double> t1_eppi0_sim(eppi0_sim_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_eppi0_sim(eppi0_sim_reader, "open_angle_ep2");
        TTreeReaderValue<double> theta_neutral_neutral_eppi0_sim(eppi0_sim_reader, "theta_pi0_pi0");
        TTreeReaderValue<double> Emiss2_eppi0_sim(eppi0_sim_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_eppi0_sim(eppi0_sim_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_eppi0_sim(eppi0_sim_reader, "pTmiss");

        // For eppi0_misID_sim_reader
        TTreeReaderValue<double> phi_eppi0_misID_sim(eppi0_misID_sim_reader, "phi");
        TTreeReaderValue<double> xB_eppi0_misID_sim(eppi0_misID_sim_reader, "x");
        TTreeReaderValue<double> Q2_eppi0_misID_sim(eppi0_misID_sim_reader, "Q2");
        TTreeReaderValue<double> t1_eppi0_misID_sim(eppi0_misID_sim_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_eppi0_misID_sim(eppi0_misID_sim_reader, "open_angle_ep2");
        TTreeReaderValue<double> theta_neutral_neutral_eppi0_misID_sim(eppi0_misID_sim_reader, "theta_gamma_gamma");
        TTreeReaderValue<double> Emiss2_eppi0_misID_sim(eppi0_misID_sim_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_eppi0_misID_sim(eppi0_misID_sim_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_eppi0_misID_sim(eppi0_misID_sim_reader, "pTmiss");

        // Loop over data_reader to fill h_DVCS_data
        while (data_reader.Next()) {
            double phi_deg = *phi_data * RAD_TO_DEG;

            for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
                const auto& bin = bin_boundaries[relevant_bins[idx]];

                if ((*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                     *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                     std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high) &&
                    apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data,
                                         *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {

                    h_DVCS_data[period][idx]->Fill(phi_deg);
                }
            }
        }
        data_reader.Restart();

        // Loop over eppi0_data_reader to fill h_eppi0_data
        while (eppi0_data_reader.Next()) {
            double phi_deg = *phi_eppi0 * RAD_TO_DEG;

            for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
                const auto& bin = bin_boundaries[relevant_bins[idx]];

                if ((*xB_eppi0 >= bin.xB_low && *xB_eppi0 <= bin.xB_high &&
                     *Q2_eppi0 >= bin.Q2_low && *Q2_eppi0 <= bin.Q2_high &&
                     std::abs(*t1_eppi0) >= bin.t_low && std::abs(*t1_eppi0) <= bin.t_high) &&
                    apply_kinematic_cuts(*t1_eppi0, *open_angle_ep2_eppi0, *theta_neutral_neutral_eppi0,
                                         *Emiss2_eppi0, *Mx2_1_eppi0, *pTmiss_eppi0)) {

                    h_eppi0_data[period][idx]->Fill(phi_deg);
                }
            }
        }
        eppi0_data_reader.Restart();

        // Loop over eppi0_sim_reader to fill h_eppi0_sim
        while (eppi0_sim_reader.Next()) {
            double phi_deg = *phi_eppi0_sim * RAD_TO_DEG;

            for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
                const auto& bin = bin_boundaries[relevant_bins[idx]];

                if ((*xB_eppi0_sim >= bin.xB_low && *xB_eppi0_sim <= bin.xB_high &&
                     *Q2_eppi0_sim >= bin.Q2_low && *Q2_eppi0_sim <= bin.Q2_high &&
                     std::abs(*t1_eppi0_sim) >= bin.t_low && std::abs(*t1_eppi0_sim) <= bin.t_high) &&
                    apply_kinematic_cuts(*t1_eppi0_sim, *open_angle_ep2_eppi0_sim, *theta_neutral_neutral_eppi0_sim,
                                         *Emiss2_eppi0_sim, *Mx2_1_eppi0_sim, *pTmiss_eppi0_sim)) {

                    h_eppi0_sim[period][idx]->Fill(phi_deg);
                }
            }
        }
        eppi0_sim_reader.Restart();

        // Loop over eppi0_misID_sim_reader to fill h_eppi0_misID_sim
        while (eppi0_misID_sim_reader.Next()) {
            double phi_deg = *phi_eppi0_misID_sim * RAD_TO_DEG;

            for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
                const auto& bin = bin_boundaries[relevant_bins[idx]];

                if ((*xB_eppi0_misID_sim >= bin.xB_low && *xB_eppi0_misID_sim <= bin.xB_high &&
                     *Q2_eppi0_misID_sim >= bin.Q2_low && *Q2_eppi0_misID_sim <= bin.Q2_high &&
                     std::abs(*t1_eppi0_misID_sim) >= bin.t_low && std::abs(*t1_eppi0_misID_sim) <= bin.t_high) &&
                    apply_kinematic_cuts(*t1_eppi0_misID_sim, *open_angle_ep2_eppi0_misID_sim, *theta_neutral_neutral_eppi0_misID_sim,
                                         *Emiss2_eppi0_misID_sim, *Mx2_1_eppi0_misID_sim, *pTmiss_eppi0_misID_sim)) {

                    h_eppi0_misID_sim[period][idx]->Fill(phi_deg);
                }
            }
        }
        eppi0_misID_sim_reader.Restart();
    }

    // Calculate contamination ratio and update UnfoldingData
    for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
        for (int phi_bin = 1; phi_bin <= n_phi_bins; ++phi_bin) {
            for (int period = 0; period < n_periods; ++period) {
                // Retrieve counts
                double N_DVCS_in_data = h_DVCS_data[period][idx]->GetBinContent(phi_bin);
                double N_eppi0_in_data = h_eppi0_data[period][idx]->GetBinContent(phi_bin);
                double N_eppi0_in_sim = h_eppi0_sim[period][idx]->GetBinContent(phi_bin);
                double N_eppi0_misID_in_sim = h_eppi0_misID_sim[period][idx]->GetBinContent(phi_bin);

                double contamination = 0.0;
                double sigma_contamination = 0.0;

                if (N_DVCS_in_data > 0 && N_eppi0_in_sim > 0) {
                    double ratio = N_eppi0_in_data / N_eppi0_in_sim;
                    contamination = (N_eppi0_misID_in_sim * ratio) / N_DVCS_in_data;

                    // Compute partial derivatives
                    double dC_dN_DVCS_in_data = -contamination / N_DVCS_in_data;
                    double dC_dN_eppi0_in_data = N_eppi0_misID_in_sim / (N_eppi0_in_sim * N_DVCS_in_data);
                    double dC_dN_eppi0_in_sim = - (N_eppi0_misID_in_sim * N_eppi0_in_data) / (N_eppi0_in_sim * N_eppi0_in_sim * N_DVCS_in_data);
                    double dC_dN_eppi0_misID_in_sim = ratio / N_DVCS_in_data;

                    // Uncertainties on counts (assuming Poisson)
                    double sigma_N_DVCS_in_data = std::sqrt(N_DVCS_in_data);
                    double sigma_N_eppi0_in_data = std::sqrt(N_eppi0_in_data);
                    double sigma_N_eppi0_in_sim = std::sqrt(N_eppi0_in_sim);
                    double sigma_N_eppi0_misID_in_sim = std::sqrt(N_eppi0_misID_in_sim);

                    // Compute sigma_C^2
                    double sigma_C_squared = 
                        std::pow(dC_dN_DVCS_in_data * sigma_N_DVCS_in_data, 2) +
                        std::pow(dC_dN_eppi0_in_data * sigma_N_eppi0_in_data, 2) +
                        std::pow(dC_dN_eppi0_in_sim * sigma_N_eppi0_in_sim, 2) +
                        std::pow(dC_dN_eppi0_misID_in_sim * sigma_N_eppi0_misID_in_sim, 2);

                    sigma_contamination = std::sqrt(sigma_C_squared);
                }

                // Initialize contamination_ratio and contamination_error if not already done
                if (unfolding_data[idx].contamination_ratio.size() == 0) {
                    unfolding_data[idx].contamination_ratio.resize(n_periods, std::vector<double>(n_phi_bins, 0.0));
                    unfolding_data[idx].contamination_error.resize(n_periods, std::vector<double>(n_phi_bins, 0.0));
                }

                unfolding_data[idx].contamination_ratio[period][phi_bin - 1] = contamination;
                unfolding_data[idx].contamination_error[period][phi_bin - 1] = sigma_contamination;
            }
        }
    }

    // Plot contamination ratios
    for (int period = 0; period < n_periods; ++period) {
        TCanvas* c_contamination = new TCanvas(Form("c_contamination_%d", period),
                                               Form("Contamination Ratio %s", period_names[period].c_str()),
                                               canvas_width, canvas_height);

        c_contamination->Divide(n_columns, n_rows);

        // Vector to hold histograms to prevent them from being deleted
        std::vector<TH1D*> h_contamination_vec;

        for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
            int bin_idx = relevant_bins[idx];
            const auto& bin = bin_boundaries[bin_idx];

            c_contamination->cd(idx + 1);

            // Access the current pad
            TPad* pad = (TPad*)c_contamination->cd(idx + 1);

            // Adjust the margins of the pad
            pad->SetLeftMargin(0.18);    // Increase left margin
            pad->SetBottomMargin(0.18);  // Increase bottom margin
            pad->SetRightMargin(0.05);   // Optional: adjust right margin
            pad->SetTopMargin(0.08);     // Optional: adjust top margin

            // Create the histogram
            TH1D* h_contamination = new TH1D(Form("h_contamination_%d_%zu", period, idx),
                                             Form("%s, <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f",
                                                  period_names[period].c_str(), bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg)),
                                             n_phi_bins, 0, 360);

            // Set axis labels and styles
            h_contamination->GetXaxis()->SetTitle("#phi [deg]");
            h_contamination->GetYaxis()->SetTitle("Contamination Ratio");
            h_contamination->GetXaxis()->SetLabelSize(0.05);
            h_contamination->GetYaxis()->SetLabelSize(0.05);
            h_contamination->GetXaxis()->SetTitleSize(0.06);
            h_contamination->GetYaxis()->SetTitleSize(0.06);
            h_contamination->SetMarkerStyle(20);
            h_contamination->SetMarkerSize(1.2);

            // Fill histogram with data and errors
            for (int phi_bin = 1; phi_bin <= n_phi_bins; ++phi_bin) {
                double contamination = unfolding_data[idx].contamination_ratio[period][phi_bin - 1];
                double sigma_contamination = unfolding_data[idx].contamination_error[period][phi_bin - 1];

                h_contamination->SetBinContent(phi_bin, contamination);
                h_contamination->SetBinError(phi_bin, sigma_contamination);
            }

            h_contamination->Draw("E1");

            // Store the histogram to prevent it from being deleted
            h_contamination_vec.push_back(h_contamination);
        }

        // Save plots into the "contamination_plots" directory
        std::string contamination_plots_dir = base_output_dir + "/contamination_plots";
        std::string output_filename = contamination_plots_dir + "/contamination_ratio_" +
                                      period_names[period] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        c_contamination->SaveAs(output_filename.c_str());

        // Clean up histograms after saving the canvas
        for (auto hist : h_contamination_vec) {
            delete hist;
        }
        h_contamination_vec.clear();

        delete c_contamination;
    }

    // Clean up histograms
    for (int period = 0; period < n_periods; ++period) {
        for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
            delete h_DVCS_data[period][idx];
            delete h_eppi0_data[period][idx];
            delete h_eppi0_sim[period][idx];
            delete h_eppi0_misID_sim[period][idx];
        }
    }
}