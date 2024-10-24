// calculate_contamination.cpp

#include "calculate_contamination.h"
#include "kinematic_cuts.h"
#include "bin_helpers.h"
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <cmath>

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
    // Number of periods
    const int n_periods = 3;

    // Names of the periods
    std::vector<std::string> period_names = {"Fa18Inb", "Fa18Out", "Sp19Inb"};

    // Number of phi bins
    const int n_phi_bins = 24;

    // Determine the relevant Q²-t bins for the given xB_bin
    std::vector<int> relevant_bins = get_bin_indices_for_xB(bin_boundaries, xB_bin);
    int n_Q2t_bins = relevant_bins.size();

    // Initialize histograms for each period and bin
    std::vector<std::vector<TH1D*>> h_DVCS_data(n_periods, std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_eppi0_data(n_periods, std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_eppi0_sim(n_periods, std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_eppi0_misID_sim(n_periods, std::vector<TH1D*>(n_Q2t_bins));

    // Create histograms
    for (int period = 0; period < n_periods; ++period) {
        for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
            int bin_idx = relevant_bins[idx];

            // Histograms for DVCS data
            h_DVCS_data[period][idx] = new TH1D(Form("h_DVCS_data_%d_%zu", period, idx),
                                                Form("DVCS Data %s Bin %d", period_names[period].c_str(), bin_idx),
                                                n_phi_bins, 0, 360);

            // Histograms for eppi0 data
            h_eppi0_data[period][idx] = new TH1D(Form("h_eppi0_data_%d_%zu", period, idx),
                                                 Form("eppi0 Data %s Bin %d", period_names[period].c_str(), bin_idx),
                                                 n_phi_bins, 0, 360);

            // Histograms for eppi0 simulation
            h_eppi0_sim[period][idx] = new TH1D(Form("h_eppi0_sim_%d_%zu", period, idx),
                                                Form("eppi0 Sim %s Bin %d", period_names[period].c_str(), bin_idx),
                                                n_phi_bins, 0, 360);

            // Histograms for misidentified eppi0 in simulation
            h_eppi0_misID_sim[period][idx] = new TH1D(Form("h_eppi0_misID_sim_%d_%zu", period, idx),
                                                      Form("eppi0 MisID Sim %s Bin %d", period_names[period].c_str(), bin_idx),
                                                      n_phi_bins, 0, 360);
        }
    }

    // Loop over each period to fill histograms
    for (int period = 0; period < n_periods; ++period) {
        // Get TTreeReaders for the current period
        TTreeReader& data_reader = data_readers[period];
        TTreeReader& eppi0_data_reader = eppi0_readers[period];
        TTreeReader& eppi0_sim_reader = mc_rec_aaogen_readers[period];
        TTreeReader& eppi0_misID_sim_reader = mc_rec_eppi0_bkg_readers[period];

        // Define TTreeReaderValues for variables (similar to your existing code)
        // For data_reader
        TTreeReaderValue<int> detector1_data(data_reader, "detector1");
        TTreeReaderValue<int> detector2_data(data_reader, "detector2");
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
        TTreeReaderValue<int> detector1_eppi0(eppi0_data_reader, "detector1");
        TTreeReaderValue<int> detector2_eppi0(eppi0_data_reader, "detector2");
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
                int bin_idx = relevant_bins[idx];
                const auto& bin = bin_boundaries[bin_idx];

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
                int bin_idx = relevant_bins[idx];
                const auto& bin = bin_boundaries[bin_idx];

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
                int bin_idx = relevant_bins[idx];
                const auto& bin = bin_boundaries[bin_idx];

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
                int bin_idx = relevant_bins[idx];
                const auto& bin = bin_boundaries[bin_idx];

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
        int bin_idx = relevant_bins[idx];

        for (int phi_bin = 1; phi_bin <= n_phi_bins; ++phi_bin) {
            for (int period = 0; period < n_periods; ++period) {
                double N_DVCS_in_data = h_DVCS_data[period][idx]->GetBinContent(phi_bin);
                double N_eppi0_in_data = h_eppi0_data[period][idx]->GetBinContent(phi_bin);
                double N_eppi0_in_sim = h_eppi0_sim[period][idx]->GetBinContent(phi_bin);
                double N_eppi0_misID_in_sim = h_eppi0_misID_sim[period][idx]->GetBinContent(phi_bin);

                double contamination = 0.0;
                if (N_DVCS_in_data > 0 && N_eppi0_in_sim > 0) {
                    double ratio = N_eppi0_in_data / N_eppi0_in_sim;
                    contamination = (N_eppi0_misID_in_sim * ratio) / N_DVCS_in_data;
                }

                // Update UnfoldingData
                // Ensure contamination_ratio is initialized
                if (unfolding_data[idx].contamination_ratio.size() == 0) {
                    unfolding_data[idx].contamination_ratio.resize(n_periods, std::vector<double>(n_phi_bins, 0.0));
                }
                unfolding_data[idx].contamination_ratio[period][phi_bin - 1] = contamination;
            }
        }
    }

    // Plot contamination ratios
    for (int period = 0; period < n_periods; ++period) {
        TCanvas* c_contamination = new TCanvas(Form("c_contamination_%d", period),
                                               Form("Contamination Ratio %s", period_names[period].c_str()),
                                               1200, 800);

        int n_columns = 2;
        int n_rows = (n_Q2t_bins + 1) / 2;

        c_contamination->Divide(n_columns, n_rows);

        for (size_t idx = 0; idx < relevant_bins.size(); ++idx) {
            c_contamination->cd(idx + 1);

            TH1D* h_contamination = new TH1D(Form("h_contamination_%d_%zu", period, idx),
                                             Form("Contamination %s Bin %d", period_names[period].c_str(), idx),
                                             n_phi_bins, 0, 360);

            for (int phi_bin = 1; phi_bin <= n_phi_bins; ++phi_bin) {
                double contamination = unfolding_data[idx].contamination_ratio[period][phi_bin - 1];
                h_contamination->SetBinContent(phi_bin, contamination);
            }

            h_contamination->GetXaxis()->SetTitle("#phi [deg]");
            h_contamination->GetYaxis()->SetTitle("Contamination Ratio");
            h_contamination->Draw("HIST");

            delete h_contamination;
        }

        std::string output_filename = base_output_dir + "/contamination_ratio_" +
                                      period_names[period] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        c_contamination->SaveAs(output_filename.c_str());
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