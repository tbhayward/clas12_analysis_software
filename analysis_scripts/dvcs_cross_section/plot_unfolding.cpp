#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <string>
#include <vector>
#include "bin_boundaries.h"
#include "kinematic_cuts.h"
#include "bin_helpers.h"
#include "plot_unfolding.h"

std::vector<UnfoldingData> plot_unfolding(const std::string& base_output_dir, 
                                          const std::string& analysisType, 
                                          int xB_bin,
                                          const std::vector<BinBoundary>& bin_boundaries, 
                                          std::vector<TTreeReader>& data_readers,  // Pass by reference
                                          std::vector<TTreeReader>& mc_gen_readers,  // Pass by reference
                                          std::vector<TTreeReader>& mc_rec_readers, // Pass by reference
                                          std::vector<TTreeReader>& eppi0_readers,  // Pass by reference
                                          std::vector<TTreeReader>& mc_gen_aaogen_readers,  // Pass by reference
                                          std::vector<TTreeReader>& mc_rec_aaogen_readers) {  // Pass by reference
    // Set style to remove stat boxes
    gStyle->SetOptStat(0);

    // List of topologies and combined option
    std::vector<std::string> topologies = {"(FD,FD)", "(CD,FD)", "(CD,FT)", "combined"};
    std::string channel_dir_dvcs = "/dvcs";
    std::string channel_dir_eppi0 = "/eppi0";

    // Vector of run period names
    std::vector<std::string> period_names = {"Fa18Inb", "Fa18Out", "Sp19Inb"};

    // Vector to store all the unfolding data
    std::vector<UnfoldingData> all_unfolding_data;
    
    const std::string& period_name_0 = period_names[0];
    const std::string& period_name_1 = period_names[1];
    const std::string& period_name_2 = period_names[2];

    // Construct the output directory for this run period
    std::string output_dir_0 = base_output_dir + "/unfolded" + channel_dir_dvcs + "/" + period_name_0;
    std::string output_dir_1 = base_output_dir + "/unfolded" + channel_dir_dvcs + "/" + period_name_1;
    std::string output_dir_2 = base_output_dir + "/unfolded" + channel_dir_dvcs + "/" + period_name_2;
    std::string output_dir_3 = base_output_dir + "/unfolded" + channel_dir_eppi0 + "/" + period_name_0;
    std::string output_dir_4 = base_output_dir + "/unfolded" + channel_dir_eppi0 + "/" + period_name_1;
    std::string output_dir_5 = base_output_dir + "/unfolded" + channel_dir_eppi0 + "/" + period_name_2;

    // Use references to access the TTreeReader instances directly from the vectors
    TTreeReader& data_reader_0 = data_readers[0];
    TTreeReader& data_reader_1 = data_readers[1];
    TTreeReader& data_reader_2 = data_readers[2];
    TTreeReader& mc_gen_reader_0 = mc_gen_readers[0];
    TTreeReader& mc_gen_reader_1 = mc_gen_readers[1];
    TTreeReader& mc_gen_reader_2 = mc_gen_readers[2];
    TTreeReader& mc_rec_reader_0 = mc_rec_readers[0];
    TTreeReader& mc_rec_reader_1 = mc_rec_readers[1];
    TTreeReader& mc_rec_reader_2 = mc_rec_readers[2];

    TTreeReader& eppi0_reader_0 = eppi0_readers[0];
    TTreeReader& eppi0_reader_1 = eppi0_readers[1];
    TTreeReader& eppi0_reader_2 = eppi0_readers[2];
    TTreeReader& mc_gen_aaogen_reader_0 = mc_gen_aaogen_readers[0];
    TTreeReader& mc_gen_aaogen_reader_1 = mc_gen_aaogen_readers[1];
    TTreeReader& mc_gen_aaogen_reader_2 = mc_gen_aaogen_readers[2];
    TTreeReader& mc_rec_aaogen_reader_0 = mc_rec_aaogen_readers[0];
    TTreeReader& mc_rec_aaogen_reader_1 = mc_rec_aaogen_readers[1];
    TTreeReader& mc_rec_aaogen_reader_2 = mc_rec_aaogen_readers[2];

    // Precompute the relevant bins for the xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);
    int n_Q2t_bins = relevant_bins.size();
    std::cout << "Current xB_bin = " << xB_bin << ", Number of Q2t bins: " << n_Q2t_bins << std::endl;

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    // Determine number of rows and columns for the canvases
    int n_columns = (xB_bin == 3 || xB_bin == 4) ? 5 : std::sqrt(next_perfect_square(n_Q2t_bins));
    int n_rows = std::ceil(static_cast<double>(n_Q2t_bins) / n_columns);

    // Create histograms for data, rec MC, gen MC, and acceptance (for combined only)
    std::vector<std::vector<TH1D*>> h_data_histograms_0(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_data_histograms_1(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_data_histograms_2(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<TH1D*> h_mc_gen_histograms_0(n_Q2t_bins);
    std::vector<TH1D*> h_mc_gen_histograms_1(n_Q2t_bins);
    std::vector<TH1D*> h_mc_gen_histograms_2(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_histograms_0(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_histograms_1(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_histograms_2(n_Q2t_bins);
    std::vector<TH1D*> h_acceptance_histograms_0(n_Q2t_bins);
    std::vector<TH1D*> h_acceptance_histograms_1(n_Q2t_bins);
    std::vector<TH1D*> h_acceptance_histograms_2(n_Q2t_bins);

    std::vector<std::vector<TH1D*>> h_eppi0_histograms_0(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_eppi0_histograms_1(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_eppi0_histograms_2(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<TH1D*> h_mc_gen_aaogen_histograms_0(n_Q2t_bins);
    std::vector<TH1D*> h_mc_gen_aaogen_histograms_1(n_Q2t_bins);
    std::vector<TH1D*> h_mc_gen_aaogen_histograms_2(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_aaogen_histograms_0(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_aaogen_histograms_1(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_aaogen_histograms_2(n_Q2t_bins);
    std::vector<TH1D*> h_acceptance_eppi0_histograms_0(n_Q2t_bins);
    std::vector<TH1D*> h_acceptance_eppi0_histograms_1(n_Q2t_bins);
    std::vector<TH1D*> h_acceptance_eppi0_histograms_2(n_Q2t_bins);

    // Initialize histograms and phi bins
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            // Create title string with the channel, period name, and bin information
            std::string title_0 = Form("%s %s, <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     "dvcs", period_name_0.c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));
            std::string title_1 = Form("%s %s, <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     "dvcs", period_name_1.c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));
            std::string title_2 = Form("%s %s, <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     "dvcs", period_name_2.c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));
            std::string title_3 = Form("%s %s, <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     "eppi0", period_name_0.c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));
            std::string title_4 = Form("%s %s, <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     "eppi0", period_name_1.c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));
            std::string title_5 = Form("%s %s, <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     "eppi0", period_name_2.c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));

            // Create histogram for the data yields with the appropriate title
            h_data_histograms_0[topo_idx][idx] = new TH1D(Form("h_data_0_%zu_%d", topo_idx, idx), title_0.c_str(), 24, 0, 360);
            h_data_histograms_1[topo_idx][idx] = new TH1D(Form("h_data_1_%zu_%d", topo_idx, idx), title_1.c_str(), 24, 0, 360);
            h_data_histograms_2[topo_idx][idx] = new TH1D(Form("h_data_2_%zu_%d", topo_idx, idx), title_2.c_str(), 24, 0, 360);

            h_eppi0_histograms_0[topo_idx][idx] = new TH1D(Form("h_eppi0_0_%zu_%d", topo_idx, idx), title_0.c_str(), 24, 0, 360);
            h_eppi0_histograms_1[topo_idx][idx] = new TH1D(Form("h_eppi0_1_%zu_%d", topo_idx, idx), title_1.c_str(), 24, 0, 360);
            h_eppi0_histograms_2[topo_idx][idx] = new TH1D(Form("h_eppi0_2_%zu_%d", topo_idx, idx), title_2.c_str(), 24, 0, 360);

            // Set axis labels and format for data histograms
            h_data_histograms_0[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_data_histograms_0[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_data_histograms_0[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_data_histograms_0[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_data_histograms_0[topo_idx][idx]->GetXaxis()->SetTitle("#phi");

            h_data_histograms_1[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_data_histograms_1[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_data_histograms_1[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_data_histograms_1[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_data_histograms_1[topo_idx][idx]->GetXaxis()->SetTitle("#phi");

            h_data_histograms_2[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_data_histograms_2[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_data_histograms_2[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_data_histograms_2[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_data_histograms_2[topo_idx][idx]->GetXaxis()->SetTitle("#phi");

            h_eppi0_histograms_0[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_eppi0_histograms_0[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_eppi0_histograms_0[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_eppi0_histograms_0[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_eppi0_histograms_0[topo_idx][idx]->GetXaxis()->SetTitle("#phi");

            h_eppi0_histograms_1[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_eppi0_histograms_1[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_eppi0_histograms_1[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_eppi0_histograms_1[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_eppi0_histograms_1[topo_idx][idx]->GetXaxis()->SetTitle("#phi");

            h_eppi0_histograms_2[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_eppi0_histograms_2[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_eppi0_histograms_2[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_eppi0_histograms_2[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_eppi0_histograms_2[topo_idx][idx]->GetXaxis()->SetTitle("#phi");

            // Set Y-axis title based on whether it's a "Raw Yield" or "Unfolded Yield"
            if (topo_idx == 3) {
                h_data_histograms_0[topo_idx][idx]->GetYaxis()->SetTitle("Unfolded Yield");
                h_data_histograms_1[topo_idx][idx]->GetYaxis()->SetTitle("Unfolded Yield");
                h_data_histograms_2[topo_idx][idx]->GetYaxis()->SetTitle("Unfolded Yield");
            } else {
                h_data_histograms_0[topo_idx][idx]->GetYaxis()->SetTitle("Raw Yield");
                h_data_histograms_1[topo_idx][idx]->GetYaxis()->SetTitle("Raw Yield");
                h_data_histograms_2[topo_idx][idx]->GetYaxis()->SetTitle("Raw Yield");
            }

            if (topo_idx == 3) {
                h_eppi0_histograms_0[topo_idx][idx]->GetYaxis()->SetTitle("Unfolded Yield");
                h_eppi0_histograms_1[topo_idx][idx]->GetYaxis()->SetTitle("Unfolded Yield");
                h_eppi0_histograms_2[topo_idx][idx]->GetYaxis()->SetTitle("Unfolded Yield");
            } else {
                h_eppi0_histograms_0[topo_idx][idx]->GetYaxis()->SetTitle("Raw Yield");
                h_eppi0_histograms_1[topo_idx][idx]->GetYaxis()->SetTitle("Raw Yield");
                h_eppi0_histograms_2[topo_idx][idx]->GetYaxis()->SetTitle("Raw Yield");
            }

            // Set markers and draw vertical error bars without horizontal error bars
            h_data_histograms_0[topo_idx][idx]->SetMarkerStyle(20);
            h_data_histograms_0[topo_idx][idx]->SetMarkerSize(1.2);
            h_data_histograms_0[topo_idx][idx]->SetDrawOption("P E1");
            h_data_histograms_1[topo_idx][idx]->SetMarkerStyle(20);
            h_data_histograms_1[topo_idx][idx]->SetMarkerSize(1.2);
            h_data_histograms_1[topo_idx][idx]->SetDrawOption("P E1");
            h_data_histograms_2[topo_idx][idx]->SetMarkerStyle(20);
            h_data_histograms_2[topo_idx][idx]->SetMarkerSize(1.2);
            h_data_histograms_2[topo_idx][idx]->SetDrawOption("P E1");

            h_eppi0_histograms_0[topo_idx][idx]->SetMarkerStyle(20);
            h_eppi0_histograms_0[topo_idx][idx]->SetMarkerSize(1.2);
            h_eppi0_histograms_0[topo_idx][idx]->SetDrawOption("P E1");
            h_eppi0_histograms_1[topo_idx][idx]->SetMarkerStyle(20);
            h_eppi0_histograms_1[topo_idx][idx]->SetMarkerSize(1.2);
            h_eppi0_histograms_1[topo_idx][idx]->SetDrawOption("P E1");
            h_eppi0_histograms_2[topo_idx][idx]->SetMarkerStyle(20);
            h_eppi0_histograms_2[topo_idx][idx]->SetMarkerSize(1.2);
            h_eppi0_histograms_2[topo_idx][idx]->SetDrawOption("P E1");

            // Create histograms for MC and acceptance if it's the combined histogram
            if (topo_idx == 3) {
                h_mc_gen_histograms_0[idx] = new TH1D(Form("h_mc_gen_combined_0_%d", idx), title_0.c_str(), 24, 0, 360);
                h_mc_rec_histograms_0[idx] = new TH1D(Form("h_mc_rec_combined_0_%d", idx), title_0.c_str(), 24, 0, 360);
                h_acceptance_histograms_0[idx] = new TH1D(Form("h_acceptance_combined_0_%d", idx), title_0.c_str(), 24, 0, 360);

                h_mc_gen_histograms_1[idx] = new TH1D(Form("h_mc_gen_combined_1_%d", idx), title_1.c_str(), 24, 0, 360);
                h_mc_rec_histograms_1[idx] = new TH1D(Form("h_mc_rec_combined_1_%d", idx), title_1.c_str(), 24, 0, 360);
                h_acceptance_histograms_1[idx] = new TH1D(Form("h_acceptance_combined_1_%d", idx), title_1.c_str(), 24, 0, 360);

                h_mc_gen_histograms_2[idx] = new TH1D(Form("h_mc_gen_combined_2_%d", idx), title_2.c_str(), 24, 0, 360);
                h_mc_rec_histograms_2[idx] = new TH1D(Form("h_mc_rec_combined_2_%d", idx), title_2.c_str(), 24, 0, 360);
                h_acceptance_histograms_2[idx] = new TH1D(Form("h_acceptance_combined_2_%d", idx), title_2.c_str(), 24, 0, 360);
            }

            if (topo_idx == 3) {
                h_mc_gen_aaogen_histograms_0[idx] = new TH1D(Form("h_mc_gen_aaogen_combined_0_%d", idx), title_0.c_str(), 24, 0, 360);
                h_mc_rec_aaogen_histograms_0[idx] = new TH1D(Form("h_mc_rec_aaogen_combined_0_%d", idx), title_0.c_str(), 24, 0, 360);
                h_acceptance_eppi0_histograms_0[idx] = new TH1D(Form("h_acceptance_eppi0_combined_0_%d", idx), title_0.c_str(), 24, 0, 360);

                h_mc_gen_aaogen_histograms_1[idx] = new TH1D(Form("h_mc_gen_aaogen_combined_1_%d", idx), title_1.c_str(), 24, 0, 360);
                h_mc_rec_aaogen_histograms_1[idx] = new TH1D(Form("h_mc_rec_aaogen_combined_1_%d", idx), title_1.c_str(), 24, 0, 360);
                h_acceptance_eppi0_histograms_1[idx] = new TH1D(Form("h_acceptance_eppi0_combined_1_%d", idx), title_1.c_str(), 24, 0, 360);

                h_mc_gen_aaogen_histograms_2[idx] = new TH1D(Form("h_mc_gen_aaogen_combined_2_%d", idx), title_2.c_str(), 24, 0, 360);
                h_mc_rec_aaogen_histograms_2[idx] = new TH1D(Form("h_mc_rec_aaogen_combined_2_%d", idx), title_2.c_str(), 24, 0, 360);
                h_acceptance_eppi0_histograms_2[idx] = new TH1D(Form("h_acceptance_eppi0_combined_2_%d", idx), title_2.c_str(), 24, 0, 360);
            }
        }
    }

    // Readers for necessary branches in all datasets (data, mc_gen, mc_rec)
    TTreeReaderValue<int> detector1_data_0(data_reader_0, "detector1");
    TTreeReaderValue<int> detector2_data_0(data_reader_0, "detector2");
    TTreeReaderValue<double> phi_data_0(data_reader_0, "phi");
    TTreeReaderValue<double> xB_data_0(data_reader_0, "x");
    TTreeReaderValue<double> Q2_data_0(data_reader_0, "Q2");
    TTreeReaderValue<double> t1_data_0(data_reader_0, "t1");
    TTreeReaderValue<double> open_angle_ep2_data_0(data_reader_0, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_data_0(data_reader_0, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_data_0(data_reader_0, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data_0(data_reader_0, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data_0(data_reader_0, "pTmiss");

    // Readers for necessary branches in all datasets (data, mc_gen, mc_rec)
    TTreeReaderValue<int> detector1_eppi0_0(eppi0_reader_0, "detector1");
    TTreeReaderValue<int> detector2_eppi0_0(eppi0_reader_0, "detector2");
    TTreeReaderValue<double> phi_eppi0_0(eppi0_reader_0, "phi");
    TTreeReaderValue<double> xB_eppi0_0(eppi0_reader_0, "x");
    TTreeReaderValue<double> Q2_eppi0_0(eppi0_reader_0, "Q2");
    TTreeReaderValue<double> t1_eppi0_0(eppi0_reader_0, "t1");
    TTreeReaderValue<double> open_angle_ep2_eppi0_0(eppi0_reader_0, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_eppi0_0(eppi0_reader_0, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_eppi0_0(eppi0_reader_0, "Emiss2");
    TTreeReaderValue<double> Mx2_1_eppi0_0(eppi0_reader_0, "Mx2_1");
    TTreeReaderValue<double> pTmiss_eppi0_0(eppi0_reader_0, "pTmiss");

    TTreeReaderValue<double> phi_mc_gen_0(mc_gen_reader_0, "phi");
    TTreeReaderValue<double> xB_mc_gen_0(mc_gen_reader_0, "x");
    TTreeReaderValue<double> Q2_mc_gen_0(mc_gen_reader_0, "Q2");
    TTreeReaderValue<double> t1_mc_gen_0(mc_gen_reader_0, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen_0(mc_gen_reader_0, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen_0(mc_gen_reader_0, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_mc_gen_0(mc_gen_reader_0, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen_0(mc_gen_reader_0, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen_0(mc_gen_reader_0, "pTmiss");

    TTreeReaderValue<double> phi_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "phi");
    TTreeReaderValue<double> xB_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "x");
    TTreeReaderValue<double> Q2_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "Q2");
    TTreeReaderValue<double> t1_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen_aaogen_0(mc_gen_aaogen_reader_0, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec_0(mc_rec_reader_0, "phi");
    TTreeReaderValue<double> xB_mc_rec_0(mc_rec_reader_0, "x");
    TTreeReaderValue<double> Q2_mc_rec_0(mc_rec_reader_0, "Q2");
    TTreeReaderValue<double> t1_mc_rec_0(mc_rec_reader_0, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec_0(mc_rec_reader_0, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec_0(mc_rec_reader_0, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_mc_rec_0(mc_rec_reader_0, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec_0(mc_rec_reader_0, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec_0(mc_rec_reader_0, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "phi");
    TTreeReaderValue<double> xB_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "x");
    TTreeReaderValue<double> Q2_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "Q2");
    TTreeReaderValue<double> t1_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec_aaogen_0(mc_rec_aaogen_reader_0, "pTmiss");




    // Readers for necessary branches in all datasets (data, mc_gen, mc_rec)
    TTreeReaderValue<int> detector1_data_1(data_reader_1, "detector1");
    TTreeReaderValue<int> detector2_data_1(data_reader_1, "detector2");
    TTreeReaderValue<double> phi_data_1(data_reader_1, "phi");
    TTreeReaderValue<double> xB_data_1(data_reader_1, "x");
    TTreeReaderValue<double> Q2_data_1(data_reader_1, "Q2");
    TTreeReaderValue<double> t1_data_1(data_reader_1, "t1");
    TTreeReaderValue<double> open_angle_ep2_data_1(data_reader_1, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_data_1(data_reader_1, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_data_1(data_reader_1, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data_1(data_reader_1, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data_1(data_reader_1, "pTmiss");

    // Readers for necessary branches in all datasets (data, mc_gen, mc_rec)
    TTreeReaderValue<int> detector1_eppi0_1(eppi0_reader_1, "detector1");
    TTreeReaderValue<int> detector2_eppi0_1(eppi0_reader_1, "detector2");
    TTreeReaderValue<double> phi_eppi0_1(eppi0_reader_1, "phi");
    TTreeReaderValue<double> xB_eppi0_1(eppi0_reader_1, "x");
    TTreeReaderValue<double> Q2_eppi0_1(eppi0_reader_1, "Q2");
    TTreeReaderValue<double> t1_eppi0_1(eppi0_reader_1, "t1");
    TTreeReaderValue<double> open_angle_ep2_eppi0_1(eppi0_reader_1, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_eppi0_1(eppi0_reader_1, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_eppi0_1(eppi0_reader_1, "Emiss2");
    TTreeReaderValue<double> Mx2_1_eppi0_1(eppi0_reader_1, "Mx2_1");
    TTreeReaderValue<double> pTmiss_eppi0_1(eppi0_reader_1, "pTmiss");

    TTreeReaderValue<double> phi_mc_gen_1(mc_gen_reader_1, "phi");
    TTreeReaderValue<double> xB_mc_gen_1(mc_gen_reader_1, "x");
    TTreeReaderValue<double> Q2_mc_gen_1(mc_gen_reader_1, "Q2");
    TTreeReaderValue<double> t1_mc_gen_1(mc_gen_reader_1, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen_1(mc_gen_reader_1, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen_1(mc_gen_reader_1, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_mc_gen_1(mc_gen_reader_1, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen_1(mc_gen_reader_1, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen_1(mc_gen_reader_1, "pTmiss");

    TTreeReaderValue<double> phi_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "phi");
    TTreeReaderValue<double> xB_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "x");
    TTreeReaderValue<double> Q2_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "Q2");
    TTreeReaderValue<double> t1_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen_aaogen_1(mc_gen_aaogen_reader_1, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec_1(mc_rec_reader_1, "phi");
    TTreeReaderValue<double> xB_mc_rec_1(mc_rec_reader_1, "x");
    TTreeReaderValue<double> Q2_mc_rec_1(mc_rec_reader_1, "Q2");
    TTreeReaderValue<double> t1_mc_rec_1(mc_rec_reader_1, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec_1(mc_rec_reader_1, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec_1(mc_rec_reader_1, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_mc_rec_1(mc_rec_reader_1, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec_1(mc_rec_reader_1, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec_1(mc_rec_reader_1, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "phi");
    TTreeReaderValue<double> xB_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "x");
    TTreeReaderValue<double> Q2_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "Q2");
    TTreeReaderValue<double> t1_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec_aaogen_1(mc_rec_aaogen_reader_1, "pTmiss");


    // Readers for necessary branches in all datasets (data, mc_gen, mc_rec)
    TTreeReaderValue<int> detector1_data_2(data_reader_2, "detector1");
    TTreeReaderValue<int> detector2_data_2(data_reader_2, "detector2");
    TTreeReaderValue<double> phi_data_2(data_reader_2, "phi");
    TTreeReaderValue<double> xB_data_2(data_reader_2, "x");
    TTreeReaderValue<double> Q2_data_2(data_reader_2, "Q2");
    TTreeReaderValue<double> t1_data_2(data_reader_2, "t1");
    TTreeReaderValue<double> open_angle_ep2_data_2(data_reader_2, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_data_2(data_reader_2, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_data_2(data_reader_2, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data_2(data_reader_2, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data_2(data_reader_2, "pTmiss");

    // Readers for necessary branches in all datasets (data, mc_gen, mc_rec)
    TTreeReaderValue<int> detector1_eppi0_2(eppi0_reader_2, "detector1");
    TTreeReaderValue<int> detector2_eppi0_2(eppi0_reader_2, "detector2");
    TTreeReaderValue<double> phi_eppi0_2(eppi0_reader_2, "phi");
    TTreeReaderValue<double> xB_eppi0_2(eppi0_reader_2, "x");
    TTreeReaderValue<double> Q2_eppi0_2(eppi0_reader_2, "Q2");
    TTreeReaderValue<double> t1_eppi0_2(eppi0_reader_2, "t1");
    TTreeReaderValue<double> open_angle_ep2_eppi0_2(eppi0_reader_2, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_eppi0_2(eppi0_reader_2, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_eppi0_2(eppi0_reader_2, "Emiss2");
    TTreeReaderValue<double> Mx2_1_eppi0_2(eppi0_reader_2, "Mx2_1");
    TTreeReaderValue<double> pTmiss_eppi0_2(eppi0_reader_2, "pTmiss");

    TTreeReaderValue<double> phi_mc_gen_2(mc_gen_reader_2, "phi");
    TTreeReaderValue<double> xB_mc_gen_2(mc_gen_reader_2, "x");
    TTreeReaderValue<double> Q2_mc_gen_2(mc_gen_reader_2, "Q2");
    TTreeReaderValue<double> t1_mc_gen_2(mc_gen_reader_2, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen_2(mc_gen_reader_2, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen_2(mc_gen_reader_2, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_mc_gen_2(mc_gen_reader_2, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen_2(mc_gen_reader_2, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen_2(mc_gen_reader_2, "pTmiss");

    TTreeReaderValue<double> phi_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "phi");
    TTreeReaderValue<double> xB_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "x");
    TTreeReaderValue<double> Q2_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "Q2");
    TTreeReaderValue<double> t1_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen_aaogen_2(mc_gen_aaogen_reader_2, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec_2(mc_rec_reader_2, "phi");
    TTreeReaderValue<double> xB_mc_rec_2(mc_rec_reader_2, "x");
    TTreeReaderValue<double> Q2_mc_rec_2(mc_rec_reader_2, "Q2");
    TTreeReaderValue<double> t1_mc_rec_2(mc_rec_reader_2, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec_2(mc_rec_reader_2, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec_2(mc_rec_reader_2, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_mc_rec_2(mc_rec_reader_2, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec_2(mc_rec_reader_2, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec_2(mc_rec_reader_2, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "phi");
    TTreeReaderValue<double> xB_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "x");
    TTreeReaderValue<double> Q2_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "Q2");
    TTreeReaderValue<double> t1_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "theta_pi0_pi0");
    TTreeReaderValue<double> Emiss2_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec_aaogen_2(mc_rec_aaogen_reader_2, "pTmiss");

    // Print before starting loops
    std::cout << "starting dvcs data for fa18 inb" << std::endl;
    // Loop over the data reader and fill histograms
    while (data_reader_0.Next()) {
        double phi_deg = *phi_data_0 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_data_0 >= bin.xB_low && *xB_data_0 <= bin.xB_high &&
                *Q2_data_0 >= bin.Q2_low && *Q2_data_0 <= bin.Q2_high &&
                std::abs(*t1_data_0) >= bin.t_low && std::abs(*t1_data_0) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_data_0, *open_angle_ep2_data_0, *theta_neutral_neutral_data_0, *Emiss2_data_0, *Mx2_1_data_0, *pTmiss_data_0)) {

                // Fill data histograms based on detector topologies
                if (*detector1_data_0 == 1 && *detector2_data_0 == 1) {  // (FD,FD)
                    h_data_histograms_0[0][idx]->Fill(phi_deg);
                } else if (*detector1_data_0 == 2 && *detector2_data_0 == 1) {  // (CD,FD)
                    h_data_histograms_0[1][idx]->Fill(phi_deg);
                } else if (*detector1_data_0 == 2 && *detector2_data_0 == 0) {  // (CD,FT)
                    h_data_histograms_0[2][idx]->Fill(phi_deg);
                }

                // Also fill the combined histogram
                h_data_histograms_0[3][idx]->Fill(phi_deg);
            }
        }
    }

    // Print before starting loops
    std::cout << "starting dvcs data for fa18 out" << std::endl;
    // Loop over the data reader and fill histograms
    while (data_reader_1.Next()) {
        double phi_deg = *phi_data_1 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_data_1 >= bin.xB_low && *xB_data_1 <= bin.xB_high &&
                *Q2_data_1 >= bin.Q2_low && *Q2_data_1 <= bin.Q2_high &&
                std::abs(*t1_data_1) >= bin.t_low && std::abs(*t1_data_1) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_data_1, *open_angle_ep2_data_1, *theta_neutral_neutral_data_1, *Emiss2_data_1, *Mx2_1_data_1, *pTmiss_data_1)) {

                // Fill data histograms based on detector topologies
                if (*detector1_data_1 == 1 && *detector2_data_1 == 1) {  // (FD,FD)
                    h_data_histograms_1[0][idx]->Fill(phi_deg);
                } else if (*detector1_data_1 == 2 && *detector2_data_1 == 1) {  // (CD,FD)
                    h_data_histograms_1[1][idx]->Fill(phi_deg);
                } else if (*detector1_data_1 == 2 && *detector2_data_1 == 0) {  // (CD,FT)
                    h_data_histograms_1[2][idx]->Fill(phi_deg);
                }

                // Also fill the combined histogram
                h_data_histograms_1[3][idx]->Fill(phi_deg);
            }
        }
    }

    // Print before starting loops
    std::cout << "starting dvcs data for sp19 inb" << std::endl;
    // Loop over the data reader and fill histograms
    while (data_reader_2.Next()) {
        double phi_deg = *phi_data_2 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_data_2 >= bin.xB_low && *xB_data_2 <= bin.xB_high &&
                *Q2_data_2 >= bin.Q2_low && *Q2_data_2 <= bin.Q2_high &&
                std::abs(*t1_data_2) >= bin.t_low && std::abs(*t1_data_2) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_data_2, *open_angle_ep2_data_2, *theta_neutral_neutral_data_2, *Emiss2_data_2, *Mx2_1_data_2, *pTmiss_data_2)) {

                // Fill data histograms based on detector topologies
                if (*detector1_data_2 == 1 && *detector2_data_2 == 1) {  // (FD,FD)
                    h_data_histograms_2[0][idx]->Fill(phi_deg);
                } else if (*detector1_data_2 == 2 && *detector2_data_2 == 1) {  // (CD,FD)
                    h_data_histograms_2[1][idx]->Fill(phi_deg);
                } else if (*detector1_data_2 == 2 && *detector2_data_2 == 0) {  // (CD,FT)
                    h_data_histograms_2[2][idx]->Fill(phi_deg);
                }

                // Also fill the combined histogram
                h_data_histograms_2[3][idx]->Fill(phi_deg);
            }
        }
    }



    // Print before starting loops
    std::cout << "starting eppi0 data for fa18 inb" << std::endl;
    // Loop over the data reader and fill histograms
    while (eppi0_reader_0.Next()) {
        double phi_deg = *phi_eppi0_0 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_eppi0_0 >= bin.xB_low && *xB_eppi0_0 <= bin.xB_high &&
                *Q2_eppi0_0 >= bin.Q2_low && *Q2_eppi0_0 <= bin.Q2_high &&
                std::abs(*t1_eppi0_0) >= bin.t_low && std::abs(*t1_eppi0_0) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_eppi0_0, *open_angle_ep2_eppi0_0, *theta_neutral_neutral_eppi0_0, *Emiss2_eppi0_0, *Mx2_1_eppi0_0, *pTmiss_eppi0_0)) {

                // Fill data histograms based on detector topologies
                if (*detector1_eppi0_0 == 1 && *detector2_eppi0_0 == 1) {  // (FD,FD)
                    h_eppi0_histograms_0[0][idx]->Fill(phi_deg);
                } else if (*detector1_eppi0_0 == 2 && *detector2_eppi0_0 == 1) {  // (CD,FD)
                    h_eppi0_histograms_0[1][idx]->Fill(phi_deg);
                } else if (*detector1_eppi0_0 == 2 && *detector2_eppi0_0 == 0) {  // (CD,FT)
                    h_eppi0_histograms_0[2][idx]->Fill(phi_deg);
                }

                // Also fill the combined histogram
                h_eppi0_histograms_0[3][idx]->Fill(phi_deg);
            }
        }
    }

    // Print before starting loops
    std::cout << "starting eppi0 data for fa18 out" << std::endl;
    // Loop over the data reader and fill histograms
    while (eppi0_reader_1.Next()) {
        double phi_deg = *phi_eppi0_1 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_eppi0_1 >= bin.xB_low && *xB_eppi0_1 <= bin.xB_high &&
                *Q2_eppi0_1 >= bin.Q2_low && *Q2_eppi0_1 <= bin.Q2_high &&
                std::abs(*t1_eppi0_1) >= bin.t_low && std::abs(*t1_eppi0_1) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_eppi0_1, *open_angle_ep2_eppi0_1, *theta_neutral_neutral_eppi0_1, *Emiss2_eppi0_1, *Mx2_1_eppi0_1, *pTmiss_eppi0_1)) {

                // Fill data histograms based on detector topologies
                if (*detector1_eppi0_1 == 1 && *detector2_eppi0_1 == 1) {  // (FD,FD)
                    h_eppi0_histograms_1[0][idx]->Fill(phi_deg);
                } else if (*detector1_eppi0_1 == 2 && *detector2_eppi0_1 == 1) {  // (CD,FD)
                    h_eppi0_histograms_1[1][idx]->Fill(phi_deg);
                } else if (*detector1_eppi0_1 == 2 && *detector2_eppi0_1 == 0) {  // (CD,FT)
                    h_eppi0_histograms_1[2][idx]->Fill(phi_deg);
                }

                // Also fill the combined histogram
                h_eppi0_histograms_1[3][idx]->Fill(phi_deg);
            }
        }
    }

    // Print before starting loops
    std::cout << "starting eppi0 data for sp19 inb" << std::endl;
    // Loop over the data reader and fill histograms
    while (eppi0_reader_2.Next()) {
        double phi_deg = *phi_eppi0_2 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_eppi0_2 >= bin.xB_low && *xB_eppi0_2 <= bin.xB_high &&
                *Q2_eppi0_2 >= bin.Q2_low && *Q2_eppi0_2 <= bin.Q2_high &&
                std::abs(*t1_eppi0_2) >= bin.t_low && std::abs(*t1_eppi0_2) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_eppi0_2, *open_angle_ep2_eppi0_2, *theta_neutral_neutral_eppi0_2, *Emiss2_eppi0_2, *Mx2_1_eppi0_2, *pTmiss_eppi0_2)) {

                // Fill data histograms based on detector topologies
                if (*detector1_eppi0_2 == 1 && *detector2_eppi0_2 == 1) {  // (FD,FD)
                    h_eppi0_histograms_2[0][idx]->Fill(phi_deg);
                } else if (*detector1_eppi0_2 == 2 && *detector2_eppi0_2 == 1) {  // (CD,FD)
                    h_eppi0_histograms_2[1][idx]->Fill(phi_deg);
                } else if (*detector1_eppi0_2 == 2 && *detector2_eppi0_2 == 0) {  // (CD,FT)
                    h_eppi0_histograms_2[2][idx]->Fill(phi_deg);
                }

                // Also fill the combined histogram
                h_eppi0_histograms_2[3][idx]->Fill(phi_deg);
            }
        }
    }





    std::cout << "starting dvcsgen gen mc for fa18 inb" << std::endl;
    // Loop over the MC generated reader and fill histograms
    while (mc_gen_reader_0.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen_0 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_gen_0 >= bin.xB_low && *xB_mc_gen_0 <= bin.xB_high &&
                *Q2_mc_gen_0 >= bin.Q2_low && *Q2_mc_gen_0 <= bin.Q2_high &&
                std::abs(*t1_mc_gen_0) >= bin.t_low && std::abs(*t1_mc_gen_0) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_gen_0, *open_angle_ep2_mc_gen_0, *theta_neutral_neutral_mc_gen_0, *Emiss2_mc_gen_0, *Mx2_1_mc_gen_0, *pTmiss_mc_gen_0)) {

                // Fill generated MC histogram (combined only)
                h_mc_gen_histograms_0[idx]->Fill(phi_mc_gen_deg);
            }
        }
    }

    std::cout << "starting dvcsgen gen mc for fa18 out" << std::endl;
    // Loop over the MC generated reader and fill histograms
    while (mc_gen_reader_1.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen_1 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_gen_1 >= bin.xB_low && *xB_mc_gen_1 <= bin.xB_high &&
                *Q2_mc_gen_1 >= bin.Q2_low && *Q2_mc_gen_1 <= bin.Q2_high &&
                std::abs(*t1_mc_gen_1) >= bin.t_low && std::abs(*t1_mc_gen_1) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_gen_1, *open_angle_ep2_mc_gen_1, *theta_neutral_neutral_mc_gen_1, *Emiss2_mc_gen_1, *Mx2_1_mc_gen_1, *pTmiss_mc_gen_1)) {

                // Fill generated MC histogram (combined only)
                h_mc_gen_histograms_1[idx]->Fill(phi_mc_gen_deg);
            }
        }
    }

    std::cout << "starting dvcsgen gen mc for sp19 inb" << std::endl;
    // Loop over the MC generated reader and fill histograms
    while (mc_gen_reader_2.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen_2 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_gen_2 >= bin.xB_low && *xB_mc_gen_2 <= bin.xB_high &&
                *Q2_mc_gen_2 >= bin.Q2_low && *Q2_mc_gen_2 <= bin.Q2_high &&
                std::abs(*t1_mc_gen_2) >= bin.t_low && std::abs(*t1_mc_gen_2) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_gen_2, *open_angle_ep2_mc_gen_2, *theta_neutral_neutral_mc_gen_2, *Emiss2_mc_gen_2, *Mx2_1_mc_gen_2, *pTmiss_mc_gen_2)) {

                // Fill generated MC histogram (combined only)
                h_mc_gen_histograms_2[idx]->Fill(phi_mc_gen_deg);
            }
        }
    }




    std::cout << "starting aaogen gen mc for fa18 inb" << std::endl;
    // Loop over the MC generated reader and fill histograms
    while (mc_gen_aaogen_reader_0.Next()) {
        double phi_mc_gen_aaogen_deg = *phi_mc_gen_aaogen_0 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_gen_aaogen_0 >= bin.xB_low && *xB_mc_gen_aaogen_0 <= bin.xB_high &&
                *Q2_mc_gen_aaogen_0 >= bin.Q2_low && *Q2_mc_gen_aaogen_0 <= bin.Q2_high &&
                std::abs(*t1_mc_gen_aaogen_0) >= bin.t_low && std::abs(*t1_mc_gen_aaogen_0) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_gen_aaogen_0, *open_angle_ep2_mc_gen_aaogen_0, *theta_neutral_neutral_mc_gen_aaogen_0, *Emiss2_mc_gen_aaogen_0, *Mx2_1_mc_gen_aaogen_0, *pTmiss_mc_gen_aaogen_0)) {

                // Fill generated MC histogram (combined only)
                h_mc_gen_aaogen_histograms_0[idx]->Fill(phi_mc_gen_aaogen_deg);
            }
        }
    }

    std::cout << "starting aaogen gen mc for fa18 out" << std::endl;
    // Loop over the MC generated reader and fill histograms
    while (mc_gen_aaogen_reader_1.Next()) {
        double phi_mc_gen_aaogen_deg = *phi_mc_gen_aaogen_1 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_gen_aaogen_1 >= bin.xB_low && *xB_mc_gen_aaogen_1 <= bin.xB_high &&
                *Q2_mc_gen_aaogen_1 >= bin.Q2_low && *Q2_mc_gen_aaogen_1 <= bin.Q2_high &&
                std::abs(*t1_mc_gen_aaogen_1) >= bin.t_low && std::abs(*t1_mc_gen_aaogen_1) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_gen_aaogen_1, *open_angle_ep2_mc_gen_aaogen_1, *theta_neutral_neutral_mc_gen_aaogen_1, *Emiss2_mc_gen_aaogen_1, *Mx2_1_mc_gen_aaogen_1, *pTmiss_mc_gen_aaogen_1)) {

                // Fill generated MC histogram (combined only)
                h_mc_gen_aaogen_histograms_1[idx]->Fill(phi_mc_gen_aaogen_deg);
            }
        }
    }

    std::cout << "starting aaogen gen mc for sp19 inb" << std::endl;
    // Loop over the MC generated reader and fill histograms
    while (mc_gen_aaogen_reader_2.Next()) {
        double phi_mc_gen_aaogen_deg = *phi_mc_gen_aaogen_2 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_gen_aaogen_2 >= bin.xB_low && *xB_mc_gen_aaogen_2 <= bin.xB_high &&
                *Q2_mc_gen_aaogen_2 >= bin.Q2_low && *Q2_mc_gen_aaogen_2 <= bin.Q2_high &&
                std::abs(*t1_mc_gen_aaogen_2) >= bin.t_low && std::abs(*t1_mc_gen_aaogen_2) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_gen_aaogen_2, *open_angle_ep2_mc_gen_aaogen_2, *theta_neutral_neutral_mc_gen_aaogen_2, *Emiss2_mc_gen_aaogen_2, *Mx2_1_mc_gen_aaogen_2, *pTmiss_mc_gen_aaogen_2)) {

                // Fill generated MC histogram (combined only)
                h_mc_gen_aaogen_histograms_2[idx]->Fill(phi_mc_gen_aaogen_deg);
            }
        }
    }






    std::cout << "starting dvcsgen rec mc for fa18 inb" << std::endl;
    // Loop over the MC reconstructed reader and fill histograms
    while (mc_rec_reader_0.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec_0 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_rec_0 >= bin.xB_low && *xB_mc_rec_0 <= bin.xB_high &&
                *Q2_mc_rec_0 >= bin.Q2_low && *Q2_mc_rec_0 <= bin.Q2_high &&
                std::abs(*t1_mc_rec_0) >= bin.t_low && std::abs(*t1_mc_rec_0) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_rec_0, *open_angle_ep2_mc_rec_0, *theta_neutral_neutral_mc_rec_0, *Emiss2_mc_rec_0, *Mx2_1_mc_rec_0, *pTmiss_mc_rec_0)) {

                // Fill reconstructed MC histogram (combined only)
                h_mc_rec_histograms_0[idx]->Fill(phi_mc_rec_deg);
            }
        }
    }

    std::cout << "starting dvcsgen rec mc for fa18 out" << std::endl;
    // Loop over the MC reconstructed reader and fill histograms
    while (mc_rec_reader_1.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec_1 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_rec_1 >= bin.xB_low && *xB_mc_rec_1 <= bin.xB_high &&
                *Q2_mc_rec_1 >= bin.Q2_low && *Q2_mc_rec_1 <= bin.Q2_high &&
                std::abs(*t1_mc_rec_1) >= bin.t_low && std::abs(*t1_mc_rec_1) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_rec_1, *open_angle_ep2_mc_rec_1, *theta_neutral_neutral_mc_rec_1, *Emiss2_mc_rec_1, *Mx2_1_mc_rec_1, *pTmiss_mc_rec_1)) {

                // Fill reconstructed MC histogram (combined only)
                h_mc_rec_histograms_1[idx]->Fill(phi_mc_rec_deg);
            }
        }
    }

    std::cout << "starting dvcsgen rec mc for sp19 inb" << std::endl;
    // Loop over the MC reconstructed reader and fill histograms
    while (mc_rec_reader_2.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec_2 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_rec_2 >= bin.xB_low && *xB_mc_rec_2 <= bin.xB_high &&
                *Q2_mc_rec_2 >= bin.Q2_low && *Q2_mc_rec_2 <= bin.Q2_high &&
                std::abs(*t1_mc_rec_2) >= bin.t_low && std::abs(*t1_mc_rec_2) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_rec_2, *open_angle_ep2_mc_rec_2, *theta_neutral_neutral_mc_rec_2, *Emiss2_mc_rec_2, *Mx2_1_mc_rec_2, *pTmiss_mc_rec_2)) {

                // Fill reconstructed MC histogram (combined only)
                h_mc_rec_histograms_2[idx]->Fill(phi_mc_rec_deg);
            }
        }
    }




    std::cout << "starting aaogen rec mc for fa18 inb" << std::endl;
    // Loop over the MC reconstructed reader and fill histograms
    while (mc_rec_aaogen_reader_0.Next()) {
        double phi_mc_rec_aaogen_deg = *phi_mc_rec_aaogen_0 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_rec_aaogen_0 >= bin.xB_low && *xB_mc_rec_aaogen_0 <= bin.xB_high &&
                *Q2_mc_rec_aaogen_0 >= bin.Q2_low && *Q2_mc_rec_aaogen_0 <= bin.Q2_high &&
                std::abs(*t1_mc_rec_aaogen_0) >= bin.t_low && std::abs(*t1_mc_rec_aaogen_0) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_rec_aaogen_0, *open_angle_ep2_mc_rec_aaogen_0, *theta_neutral_neutral_mc_rec_aaogen_0, *Emiss2_mc_rec_aaogen_0, *Mx2_1_mc_rec_aaogen_0, *pTmiss_mc_rec_aaogen_0)) {

                // Fill reconstructed MC histogram (combined only)
                h_mc_rec_aaogen_histograms_0[idx]->Fill(phi_mc_rec_aaogen_deg);
            }
        }
    }

    std::cout << "starting aaogen rec mc for fa18 out" << std::endl;
    // Loop over the MC reconstructed reader and fill histograms
    while (mc_rec_aaogen_reader_1.Next()) {
        double phi_mc_rec_aaogen_deg = *phi_mc_rec_aaogen_1 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_rec_aaogen_1 >= bin.xB_low && *xB_mc_rec_aaogen_1 <= bin.xB_high &&
                *Q2_mc_rec_aaogen_1 >= bin.Q2_low && *Q2_mc_rec_aaogen_1 <= bin.Q2_high &&
                std::abs(*t1_mc_rec_aaogen_1) >= bin.t_low && std::abs(*t1_mc_rec_aaogen_1) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_rec_aaogen_1, *open_angle_ep2_mc_rec_aaogen_1, *theta_neutral_neutral_mc_rec_aaogen_1, *Emiss2_mc_rec_aaogen_1, *Mx2_1_mc_rec_aaogen_1, *pTmiss_mc_rec_aaogen_1)) {

                // Fill reconstructed MC histogram (combined only)
                h_mc_rec_aaogen_histograms_1[idx]->Fill(phi_mc_rec_aaogen_deg);
            }
        }
    }

    std::cout << "starting aaogen rec mc for sp19 inb" << std::endl;
    // Loop over the MC reconstructed reader and fill histograms
    while (mc_rec_aaogen_reader_2.Next()) {
        double phi_mc_rec_aaogen_deg = *phi_mc_rec_aaogen_2 * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_rec_aaogen_2 >= bin.xB_low && *xB_mc_rec_aaogen_2 <= bin.xB_high &&
                *Q2_mc_rec_aaogen_2 >= bin.Q2_low && *Q2_mc_rec_aaogen_2 <= bin.Q2_high &&
                std::abs(*t1_mc_rec_aaogen_2) >= bin.t_low && std::abs(*t1_mc_rec_aaogen_2) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_rec_aaogen_2, *open_angle_ep2_mc_rec_aaogen_2, *theta_neutral_neutral_mc_rec_aaogen_2, *Emiss2_mc_rec_aaogen_2, *Mx2_1_mc_rec_aaogen_2, *pTmiss_mc_rec_aaogen_2)) {

                // Fill reconstructed MC histogram (combined only)
                h_mc_rec_aaogen_histograms_2[idx]->Fill(phi_mc_rec_aaogen_deg);
            }
        }
    }


    // Compute acceptance for combined
    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        if (h_mc_gen_histograms_0[idx]->Integral() > 0) {
            h_acceptance_histograms_0[idx]->Divide(h_mc_rec_histograms_0[idx], h_mc_gen_histograms_0[idx], 1, 1, "B");
        }

        if (h_mc_gen_histograms_1[idx]->Integral() > 0) {
            h_acceptance_histograms_1[idx]->Divide(h_mc_rec_histograms_1[idx], h_mc_gen_histograms_1[idx], 1, 1, "B");
        }

        if (h_mc_gen_histograms_2[idx]->Integral() > 0) {
            h_acceptance_histograms_2[idx]->Divide(h_mc_rec_histograms_2[idx], h_mc_gen_histograms_2[idx], 1, 1, "B");
        }

        if (h_mc_gen_aaogen_histograms_0[idx]->Integral() > 0) {
            h_acceptance_eppi0_histograms_0[idx]->Divide(h_mc_rec_aaogen_histograms_0[idx], h_mc_gen_aaogen_histograms_0[idx], 1, 1, "B");
        }

        if (h_mc_gen_aaogen_histograms_1[idx]->Integral() > 0) {
            h_acceptance_eppi0_histograms_1[idx]->Divide(h_mc_rec_aaogen_histograms_1[idx], h_mc_gen_aaogen_histograms_1[idx], 1, 1, "B");
        }

        if (h_mc_gen_aaogen_histograms_2[idx]->Integral() > 0) {
            h_acceptance_eppi0_histograms_2[idx]->Divide(h_mc_rec_aaogen_histograms_2[idx], h_mc_gen_aaogen_histograms_2[idx], 1, 1, "B");
        }
    }

    // Gather and store unfolding data
    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        const auto& bin = bin_boundaries[relevant_bins[idx]];
        UnfoldingData unfolding_data;
        
        unfolding_data.bin_number = idx;
        unfolding_data.xB_min = bin.xB_low;
        unfolding_data.xB_max = bin.xB_high;
        unfolding_data.xB_avg = bin.xB_avg;
        unfolding_data.Q2_min = bin.Q2_low;
        unfolding_data.Q2_max = bin.Q2_high;
        unfolding_data.Q2_avg = bin.Q2_avg;
        unfolding_data.t_min = bin.t_low;
        unfolding_data.t_max = bin.t_high;
        unfolding_data.t_avg = bin.t_avg;

        // Resize vectors for each period
        unfolding_data.raw_yields.resize(period_names.size()*2, std::vector<int>(topologies.size() * 24));  // Topology x Phi bins
        unfolding_data.acceptance.resize(period_names.size()*2, std::vector<double>(24));  // Acceptance for each phi_bin
        unfolding_data.unfolded_yields.resize(period_names.size()*2, std::vector<double>(24));  // Unfolded yield for each phi_bin

        // For each phi bin (24 bins from 0 to 360)
        for (int phi_bin = 1; phi_bin <= 24; ++phi_bin) {
            double phi_min = (phi_bin - 1) * 15.0;
            double phi_max = phi_bin * 15.0;
            unfolding_data.phi_min.push_back(phi_min);
            unfolding_data.phi_max.push_back(phi_max);

            // Get the raw yield for each topology
            for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
                int raw_yield_0 = h_data_histograms_0[topo_idx][idx]->GetBinContent(phi_bin);
                int raw_yield_1 = h_data_histograms_1[topo_idx][idx]->GetBinContent(phi_bin);
                int raw_yield_2 = h_data_histograms_2[topo_idx][idx]->GetBinContent(phi_bin);
                unfolding_data.raw_yields[0][topo_idx * 24 + (phi_bin - 1)] = raw_yield_0;
                unfolding_data.raw_yields[1][topo_idx * 24 + (phi_bin - 1)] = raw_yield_1;
                unfolding_data.raw_yields[2][topo_idx * 24 + (phi_bin - 1)] = raw_yield_2;

                int raw_yield_3 = h_eppi0_histograms_0[topo_idx][idx]->GetBinContent(phi_bin);
                int raw_yield_4 = h_eppi0_histograms_1[topo_idx][idx]->GetBinContent(phi_bin);
                int raw_yield_5 = h_eppi0_histograms_2[topo_idx][idx]->GetBinContent(phi_bin);
                unfolding_data.raw_yields[3][topo_idx * 24 + (phi_bin - 1)] = raw_yield_3;
                unfolding_data.raw_yields[4][topo_idx * 24 + (phi_bin - 1)] = raw_yield_4;
                unfolding_data.raw_yields[5][topo_idx * 24 + (phi_bin - 1)] = raw_yield_5;
            }

            // Store acceptance for combined topology (topo_idx == 3)
            double acceptance_value_0 = h_acceptance_histograms_0[idx]->GetBinContent(phi_bin);
            double acceptance_value_1 = h_acceptance_histograms_1[idx]->GetBinContent(phi_bin);
            double acceptance_value_2 = h_acceptance_histograms_2[idx]->GetBinContent(phi_bin);
            unfolding_data.acceptance[0][phi_bin - 1] = acceptance_value_0;
            unfolding_data.acceptance[1][phi_bin - 1] = acceptance_value_1;
            unfolding_data.acceptance[2][phi_bin - 1] = acceptance_value_2;

            double acceptance_eppi0_value_0 = h_acceptance_eppi0_histograms_0[idx]->GetBinContent(phi_bin);
            double acceptance_eppi0_value_1 = h_acceptance_eppi0_histograms_1[idx]->GetBinContent(phi_bin);
            double acceptance_eppi0_value_2 = h_acceptance_eppi0_histograms_2[idx]->GetBinContent(phi_bin);
            unfolding_data.acceptance[3][phi_bin - 1] = acceptance_eppi0_value_0;
            unfolding_data.acceptance[4][phi_bin - 1] = acceptance_eppi0_value_1;
            unfolding_data.acceptance[5][phi_bin - 1] = acceptance_eppi0_value_2;

            // Store unfolded yield for combined topology (topo_idx == 3)
            // For period 0
            if (acceptance_value_0 > 0) {
                double unfolded_yield_0 = h_data_histograms_0[3][idx]->GetBinContent(phi_bin) / acceptance_value_0;
                unfolding_data.unfolded_yields[0][phi_bin - 1] = unfolded_yield_0;
            } else {
                unfolding_data.unfolded_yields[0][phi_bin - 1] = 0.0;
            }

            // For period 1
            if (acceptance_value_1 > 0) {
                double unfolded_yield_1 = h_data_histograms_1[3][idx]->GetBinContent(phi_bin) / acceptance_value_1;
                unfolding_data.unfolded_yields[1][phi_bin - 1] = unfolded_yield_1;
            } else {
                unfolding_data.unfolded_yields[1][phi_bin - 1] = 0.0;
            }

            // For period 2
            if (acceptance_value_2 > 0) {
                double unfolded_yield_2 = h_data_histograms_2[3][idx]->GetBinContent(phi_bin) / acceptance_value_2;
                unfolding_data.unfolded_yields[2][phi_bin - 1] = unfolded_yield_2;
            } else {
                unfolding_data.unfolded_yields[2][phi_bin - 1] = 0.0;
            }


            if (acceptance_eppi0_value_0 > 0) {
                double unfolded_eppi0_yield_0 = h_eppi0_histograms_0[3][idx]->GetBinContent(phi_bin) / acceptance_eppi0_value_0;
                unfolding_data.unfolded_yields[3][phi_bin - 1] = unfolded_eppi0_yield_0;
            } else {
                unfolding_data.unfolded_yields[3][phi_bin - 1] = 0.0;
            }

            // For period 1
            if (acceptance_eppi0_value_1 > 0) {
                double unfolded_eppi0_yield_1 = h_eppi0_histograms_1[3][idx]->GetBinContent(phi_bin) / acceptance_eppi0_value_1;
                unfolding_data.unfolded_yields[4][phi_bin - 1] = unfolded_eppi0_yield_1;
            } else {
                unfolding_data.unfolded_yields[4][phi_bin - 1] = 0.0;
            }

            // For period 2
            if (acceptance_eppi0_value_2 > 0) {
                double unfolded_eppi0_yield_2 = h_eppi0_histograms_2[3][idx]->GetBinContent(phi_bin) / acceptance_eppi0_value_2;
                unfolding_data.unfolded_yields[5][phi_bin - 1] = unfolded_eppi0_yield_2;
            } else {
                unfolding_data.unfolded_yields[5][phi_bin - 1] = 0.0;
            }
        }

        // Add the unfolding data for this bin to the overall results
        all_unfolding_data.push_back(unfolding_data);
    }

    // Plot and save histograms (data divided by acceptance for combined)
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        TCanvas* canvas_yield_0 = new TCanvas(Form("c_yield_0_%zu", topo_idx), "Yields", canvas_width, canvas_height);
        TCanvas* canvas_yield_1 = new TCanvas(Form("c_yield_1_%zu", topo_idx), "Yields", canvas_width, canvas_height);
        TCanvas* canvas_yield_2 = new TCanvas(Form("c_yield_2_%zu", topo_idx), "Yields", canvas_width, canvas_height);
        canvas_yield_0->Divide(n_columns, n_rows);
        canvas_yield_1->Divide(n_columns, n_rows);
        canvas_yield_2->Divide(n_columns, n_rows);

        TCanvas* canvas_yield_3 = new TCanvas(Form("c_yield_3_%zu", topo_idx), "Yields", canvas_width, canvas_height);
        TCanvas* canvas_yield_4 = new TCanvas(Form("c_yield_4_%zu", topo_idx), "Yields", canvas_width, canvas_height);
        TCanvas* canvas_yield_5 = new TCanvas(Form("c_yield_5_%zu", topo_idx), "Yields", canvas_width, canvas_height);
        canvas_yield_3->Divide(n_columns, n_rows);
        canvas_yield_4->Divide(n_columns, n_rows);
        canvas_yield_5->Divide(n_columns, n_rows);

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            // For canvas_yield_0
            TPad* pad_yield_0 = (TPad*)canvas_yield_0->cd(idx + 1);
            pad_yield_0->SetLeftMargin(0.2);
            pad_yield_0->SetBottomMargin(0.15);

            // For canvas_yield_1
            TPad* pad_yield_1 = (TPad*)canvas_yield_1->cd(idx + 1);
            pad_yield_1->SetLeftMargin(0.2);
            pad_yield_1->SetBottomMargin(0.15);

            // For canvas_yield_2
            TPad* pad_yield_2 = (TPad*)canvas_yield_2->cd(idx + 1);
            pad_yield_2->SetLeftMargin(0.2);
            pad_yield_2->SetBottomMargin(0.15);

            // For canvas_yield_3
            TPad* pad_yield_3 = (TPad*)canvas_yield_3->cd(idx + 1);
            pad_yield_3->SetLeftMargin(0.2);
            pad_yield_3->SetBottomMargin(0.15);

            // For canvas_yield_4
            TPad* pad_yield_4 = (TPad*)canvas_yield_4->cd(idx + 1);
            pad_yield_4->SetLeftMargin(0.2);
            pad_yield_4->SetBottomMargin(0.15);

            // For canvas_yield_5
            TPad* pad_yield_5 = (TPad*)canvas_yield_5->cd(idx + 1);
            pad_yield_5->SetLeftMargin(0.2);
            pad_yield_5->SetBottomMargin(0.15);

            // If this is the combined topology, we need to unfold by dividing data by acceptance
            if (topo_idx == 3) {  // Combined histograms
                int nBins = h_data_histograms_0[topo_idx][idx]->GetNbinsX();
                for (int bin = 1; bin <= nBins; ++bin) {
                    // For period 0
                    double acceptance_value_0 = h_acceptance_histograms_0[idx]->GetBinContent(bin);
                    if (acceptance_value_0 > 0) {
                        double data_value_0 = h_data_histograms_0[topo_idx][idx]->GetBinContent(bin);
                        double data_error_0 = h_data_histograms_0[topo_idx][idx]->GetBinError(bin);
                        double unfolded_value_0 = data_value_0 / acceptance_value_0;
                        double unfolded_error_0 = data_error_0 / acceptance_value_0;
                        h_data_histograms_0[topo_idx][idx]->SetBinContent(bin, unfolded_value_0);
                        h_data_histograms_0[topo_idx][idx]->SetBinError(bin, unfolded_error_0);
                    } else {
                        h_data_histograms_0[topo_idx][idx]->SetBinContent(bin, 0);
                        h_data_histograms_0[topo_idx][idx]->SetBinError(bin, 0);
                    }
                    // For period 1
                    double acceptance_value_1 = h_acceptance_histograms_1[idx]->GetBinContent(bin);
                    if (acceptance_value_1 > 0) {
                        double data_value_1 = h_data_histograms_1[topo_idx][idx]->GetBinContent(bin);
                        double data_error_1 = h_data_histograms_1[topo_idx][idx]->GetBinError(bin);
                        double unfolded_value_1 = data_value_1 / acceptance_value_1;
                        double unfolded_error_1 = data_error_1 / acceptance_value_1;
                        h_data_histograms_1[topo_idx][idx]->SetBinContent(bin, unfolded_value_1);
                        h_data_histograms_1[topo_idx][idx]->SetBinError(bin, unfolded_error_1);
                    } else {
                        h_data_histograms_1[topo_idx][idx]->SetBinContent(bin, 0);
                        h_data_histograms_1[topo_idx][idx]->SetBinError(bin, 0);
                    }
                    // For period 2
                    double acceptance_value_2 = h_acceptance_histograms_2[idx]->GetBinContent(bin);
                    if (acceptance_value_2 > 0) {
                        double data_value_2 = h_data_histograms_2[topo_idx][idx]->GetBinContent(bin);
                        double data_error_2 = h_data_histograms_2[topo_idx][idx]->GetBinError(bin);
                        double unfolded_value_2 = data_value_2 / acceptance_value_2;
                        double unfolded_error_2 = data_error_2 / acceptance_value_2;
                        h_data_histograms_2[topo_idx][idx]->SetBinContent(bin, unfolded_value_2);
                        h_data_histograms_2[topo_idx][idx]->SetBinError(bin, unfolded_error_2);
                    } else {
                        h_data_histograms_2[topo_idx][idx]->SetBinContent(bin, 0);
                        h_data_histograms_2[topo_idx][idx]->SetBinError(bin, 0);
                    }
                    // For period 3
                    double acceptance_value_3 = h_acceptance_eppi0_histograms_0[idx]->GetBinContent(bin);
                    if (acceptance_value_3 > 0) {
                        double data_value_3 = h_eppi0_histograms_0[topo_idx][idx]->GetBinContent(bin);
                        double data_error_3 = h_eppi0_histograms_0[topo_idx][idx]->GetBinError(bin);
                        double unfolded_value_3 = data_value_3 / acceptance_value_3;
                        double unfolded_error_3 = data_error_3 / acceptance_value_3;
                        h_eppi0_histograms_0[topo_idx][idx]->SetBinContent(bin, unfolded_value_3);
                        h_eppi0_histograms_0[topo_idx][idx]->SetBinError(bin, unfolded_error_3);
                    } else {
                        h_eppi0_histograms_0[topo_idx][idx]->SetBinContent(bin, 0);
                        h_eppi0_histograms_0[topo_idx][idx]->SetBinError(bin, 0);
                    }
                    // For period 4
                    double acceptance_value_4 = h_acceptance_eppi0_histograms_1[idx]->GetBinContent(bin);
                    if (acceptance_value_4 > 0) {
                        double data_value_4 = h_eppi0_histograms_1[topo_idx][idx]->GetBinContent(bin);
                        double data_error_4 = h_eppi0_histograms_1[topo_idx][idx]->GetBinError(bin);
                        double unfolded_value_4 = data_value_4 / acceptance_value_4;
                        double unfolded_error_4 = data_error_4 / acceptance_value_4;
                        h_eppi0_histograms_1[topo_idx][idx]->SetBinContent(bin, unfolded_value_4);
                        h_eppi0_histograms_1[topo_idx][idx]->SetBinError(bin, unfolded_error_4);
                    } else {
                        h_eppi0_histograms_1[topo_idx][idx]->SetBinContent(bin, 0);
                        h_eppi0_histograms_1[topo_idx][idx]->SetBinError(bin, 0);
                    }
                    // For period 3
                    double acceptance_value_5 = h_acceptance_eppi0_histograms_2[idx]->GetBinContent(bin);
                    if (acceptance_value_5 > 0) {
                        double data_value_5 = h_eppi0_histograms_2[topo_idx][idx]->GetBinContent(bin);
                        double data_error_5 = h_eppi0_histograms_2[topo_idx][idx]->GetBinError(bin);
                        double unfolded_value_5 = data_value_5 / acceptance_value_5;
                        double unfolded_error_5 = data_error_5 / acceptance_value_5;
                        h_eppi0_histograms_2[topo_idx][idx]->SetBinContent(bin, unfolded_value_5);
                        h_eppi0_histograms_2[topo_idx][idx]->SetBinError(bin, unfolded_error_5);
                    } else {
                        h_eppi0_histograms_2[topo_idx][idx]->SetBinContent(bin, 0);
                        h_eppi0_histograms_2[topo_idx][idx]->SetBinError(bin, 0);
                    }
                }
                // Adjust y-axis range for each histogram
                double max_value_0 = h_data_histograms_0[topo_idx][idx]->GetMaximum();
                double max_value_1 = h_data_histograms_1[topo_idx][idx]->GetMaximum();
                double max_value_2 = h_data_histograms_2[topo_idx][idx]->GetMaximum();
                h_data_histograms_0[topo_idx][idx]->GetYaxis()->SetRangeUser(0, 1.1 * max_value_0);
                h_data_histograms_1[topo_idx][idx]->GetYaxis()->SetRangeUser(0, 1.1 * max_value_1);
                h_data_histograms_2[topo_idx][idx]->GetYaxis()->SetRangeUser(0, 1.1 * max_value_2);
                // Adjust y-axis range for each histogram
                double max_value_3 = h_eppi0_histograms_0[topo_idx][idx]->GetMaximum();
                double max_value_4 = h_eppi0_histograms_1[topo_idx][idx]->GetMaximum();
                double max_value_5 = h_eppi0_histograms_2[topo_idx][idx]->GetMaximum();
                h_eppi0_histograms_0[topo_idx][idx]->GetYaxis()->SetRangeUser(0, 1.1 * max_value_3);
                h_eppi0_histograms_1[topo_idx][idx]->GetYaxis()->SetRangeUser(0, 1.1 * max_value_4);
                h_eppi0_histograms_2[topo_idx][idx]->GetYaxis()->SetRangeUser(0, 1.1 * max_value_5);
            }

            // Draw the histograms on their respective canvases and pads
            pad_yield_0->cd();
            h_data_histograms_0[topo_idx][idx]->Draw("E1");

            pad_yield_1->cd();
            h_data_histograms_1[topo_idx][idx]->Draw("E1");

            pad_yield_2->cd();
            h_data_histograms_2[topo_idx][idx]->Draw("E1");

            // Draw the histograms on their respective canvases and pads
            pad_yield_3->cd();
            h_eppi0_histograms_0[topo_idx][idx]->Draw("E1");

            pad_yield_4->cd();
            h_eppi0_histograms_1[topo_idx][idx]->Draw("E1");

            pad_yield_5->cd();
            h_eppi0_histograms_2[topo_idx][idx]->Draw("E1");
        }

        // Save the canvases with unique filenames
        std::string filename_yield_0 = output_dir_0 + "/yields/yields_dvcs_" + period_names[0] + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        std::string filename_yield_1 = output_dir_1 + "/yields/yields_dvcs_" + period_names[1] + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        std::string filename_yield_2 = output_dir_2 + "/yields/yields_dvcs_" + period_names[2] + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        // Save the canvases with unique filenames
        std::string filename_yield_3 = output_dir_3 + "/yields/yields_eppi0_" + period_names[0] + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        std::string filename_yield_4 = output_dir_4 + "/yields/yields_eppi0_" + period_names[1] + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        std::string filename_yield_5 = output_dir_5 + "/yields/yields_eppi0_" + period_names[2] + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";

        canvas_yield_0->SaveAs(filename_yield_0.c_str());
        canvas_yield_1->SaveAs(filename_yield_1.c_str());
        canvas_yield_2->SaveAs(filename_yield_2.c_str());

        canvas_yield_3->SaveAs(filename_yield_3.c_str());
        canvas_yield_4->SaveAs(filename_yield_4.c_str());
        canvas_yield_5->SaveAs(filename_yield_5.c_str());

        delete canvas_yield_3;
        delete canvas_yield_4;
        delete canvas_yield_5;
    }

    // Plot and save acceptance for combined only for each period separately

    // Create canvases for each period
    TCanvas* canvas_acceptance_0 = new TCanvas("c_acceptance_combined_0", "Acceptance dvcs Combined Period 0", canvas_width, canvas_height);
    TCanvas* canvas_acceptance_1 = new TCanvas("c_acceptance_combined_1", "Acceptance dvcs Combined Period 1", canvas_width, canvas_height);
    TCanvas* canvas_acceptance_2 = new TCanvas("c_acceptance_combined_2", "Acceptance dvcs Combined Period 2", canvas_width, canvas_height);
    // Create canvases for each period
    TCanvas* canvas_acceptance_3 = new TCanvas("c_acceptance_combined_3", "Acceptance eppi0 Combined Period 0", canvas_width, canvas_height);
    TCanvas* canvas_acceptance_4 = new TCanvas("c_acceptance_combined_4", "Acceptance eppi0 Combined Period 1", canvas_width, canvas_height);
    TCanvas* canvas_acceptance_5 = new TCanvas("c_acceptance_combined_5", "Acceptance eppi0 Combined Period 2", canvas_width, canvas_height);

    canvas_acceptance_0->Divide(n_columns, n_rows);
    canvas_acceptance_1->Divide(n_columns, n_rows);
    canvas_acceptance_2->Divide(n_columns, n_rows);
    canvas_acceptance_3->Divide(n_columns, n_rows);
    canvas_acceptance_4->Divide(n_columns, n_rows);
    canvas_acceptance_5->Divide(n_columns, n_rows);

    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        // For canvas_acceptance_0
        TPad* pad_acceptance_0 = (TPad*)canvas_acceptance_0->cd(idx + 1);
        pad_acceptance_0->SetLeftMargin(0.20);
        pad_acceptance_0->SetBottomMargin(0.15);

        // Set the axis labels and styles for acceptance histograms for period 0
        h_acceptance_histograms_0[idx]->GetXaxis()->SetTitle("#phi");
        h_acceptance_histograms_0[idx]->GetYaxis()->SetTitle("Acceptance");
        h_acceptance_histograms_0[idx]->GetXaxis()->SetLabelSize(0.05);
        h_acceptance_histograms_0[idx]->GetYaxis()->SetLabelSize(0.05);
        h_acceptance_histograms_0[idx]->GetXaxis()->SetTitleSize(0.06);
        h_acceptance_histograms_0[idx]->GetYaxis()->SetTitleSize(0.06);

        h_acceptance_histograms_0[idx]->Draw("E1");

        // For canvas_acceptance_1
        TPad* pad_acceptance_1 = (TPad*)canvas_acceptance_1->cd(idx + 1);
        pad_acceptance_1->SetLeftMargin(0.20);
        pad_acceptance_1->SetBottomMargin(0.15);

        // Set the axis labels and styles for acceptance histograms for period 1
        h_acceptance_histograms_1[idx]->GetXaxis()->SetTitle("#phi");
        h_acceptance_histograms_1[idx]->GetYaxis()->SetTitle("Acceptance");
        h_acceptance_histograms_1[idx]->GetXaxis()->SetLabelSize(0.05);
        h_acceptance_histograms_1[idx]->GetYaxis()->SetLabelSize(0.05);
        h_acceptance_histograms_1[idx]->GetXaxis()->SetTitleSize(0.06);
        h_acceptance_histograms_1[idx]->GetYaxis()->SetTitleSize(0.06);

        h_acceptance_histograms_1[idx]->Draw("E1");

        // For canvas_acceptance_2
        TPad* pad_acceptance_2 = (TPad*)canvas_acceptance_2->cd(idx + 1);
        pad_acceptance_2->SetLeftMargin(0.20);
        pad_acceptance_2->SetBottomMargin(0.15);

        // Set the axis labels and styles for acceptance histograms for period 2
        h_acceptance_histograms_2[idx]->GetXaxis()->SetTitle("#phi");
        h_acceptance_histograms_2[idx]->GetYaxis()->SetTitle("Acceptance");
        h_acceptance_histograms_2[idx]->GetXaxis()->SetLabelSize(0.05);
        h_acceptance_histograms_2[idx]->GetYaxis()->SetLabelSize(0.05);
        h_acceptance_histograms_2[idx]->GetXaxis()->SetTitleSize(0.06);
        h_acceptance_histograms_2[idx]->GetYaxis()->SetTitleSize(0.06);

        h_acceptance_histograms_2[idx]->Draw("E1");

        // For canvas_acceptance_3
        TPad* pad_acceptance_3 = (TPad*)canvas_acceptance_3->cd(idx + 1);
        pad_acceptance_3->SetLeftMargin(0.20);
        pad_acceptance_3->SetBottomMargin(0.15);

        // Set the axis labels and styles for acceptance histograms for period 2
        h_acceptance_eppi0_histograms_0[idx]->GetXaxis()->SetTitle("#phi");
        h_acceptance_eppi0_histograms_0[idx]->GetYaxis()->SetTitle("Acceptance");
        h_acceptance_eppi0_histograms_0[idx]->GetXaxis()->SetLabelSize(0.05);
        h_acceptance_eppi0_histograms_0[idx]->GetYaxis()->SetLabelSize(0.05);
        h_acceptance_eppi0_histograms_0[idx]->GetXaxis()->SetTitleSize(0.06);
        h_acceptance_eppi0_histograms_0[idx]->GetYaxis()->SetTitleSize(0.06);

        h_acceptance_eppi0_histograms_0[idx]->Draw("E1");

        // For canvas_acceptance_4
        TPad* pad_acceptance_4 = (TPad*)canvas_acceptance_4->cd(idx + 1);
        pad_acceptance_4->SetLeftMargin(0.20);
        pad_acceptance_4->SetBottomMargin(0.15);

        // Set the axis labels and styles for acceptance histograms for period 2
        h_acceptance_eppi0_histograms_1[idx]->GetXaxis()->SetTitle("#phi");
        h_acceptance_eppi0_histograms_1[idx]->GetYaxis()->SetTitle("Acceptance");
        h_acceptance_eppi0_histograms_1[idx]->GetXaxis()->SetLabelSize(0.05);
        h_acceptance_eppi0_histograms_1[idx]->GetYaxis()->SetLabelSize(0.05);
        h_acceptance_eppi0_histograms_1[idx]->GetXaxis()->SetTitleSize(0.06);
        h_acceptance_eppi0_histograms_1[idx]->GetYaxis()->SetTitleSize(0.06);

        h_acceptance_eppi0_histograms_1[idx]->Draw("E1");

        // For canvas_acceptance_5
        TPad* pad_acceptance_5 = (TPad*)canvas_acceptance_5->cd(idx + 1);
        pad_acceptance_5->SetLeftMargin(0.20);
        pad_acceptance_5->SetBottomMargin(0.15);

        // Set the axis labels and styles for acceptance histograms for period 2
        h_acceptance_eppi0_histograms_2[idx]->GetXaxis()->SetTitle("#phi");
        h_acceptance_eppi0_histograms_2[idx]->GetYaxis()->SetTitle("Acceptance");
        h_acceptance_eppi0_histograms_2[idx]->GetXaxis()->SetLabelSize(0.05);
        h_acceptance_eppi0_histograms_2[idx]->GetYaxis()->SetLabelSize(0.05);
        h_acceptance_eppi0_histograms_2[idx]->GetXaxis()->SetTitleSize(0.06);
        h_acceptance_eppi0_histograms_2[idx]->GetYaxis()->SetTitleSize(0.06);

        h_acceptance_eppi0_histograms_2[idx]->Draw("E1");
    }

    // Save each canvas with a unique filename
    std::string filename_acceptance_0 = output_dir_0 + "/acceptances/acceptances_combined_dvcs_" + period_names[0] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    std::string filename_acceptance_1 = output_dir_1 + "/acceptances/acceptances_combined_dvcs_" + period_names[1] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    std::string filename_acceptance_2 = output_dir_2 + "/acceptances/acceptances_combined_dvcs_" + period_names[2] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    // Save each canvas with a unique filename
    std::string filename_acceptance_3 = output_dir_3 + "/acceptances/acceptances_combined_eppi0_" + period_names[0] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    std::string filename_acceptance_4 = output_dir_4 + "/acceptances/acceptances_combined_eppi0_" + period_names[1] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    std::string filename_acceptance_5 = output_dir_5 + "/acceptances/acceptances_combined_eppi0_" + period_names[2] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";

    canvas_acceptance_0->SaveAs(filename_acceptance_0.c_str());
    canvas_acceptance_1->SaveAs(filename_acceptance_1.c_str());
    canvas_acceptance_2->SaveAs(filename_acceptance_2.c_str());
    canvas_acceptance_3->SaveAs(filename_acceptance_3.c_str());
    canvas_acceptance_4->SaveAs(filename_acceptance_4.c_str());
    canvas_acceptance_5->SaveAs(filename_acceptance_5.c_str());

    delete canvas_acceptance_0;
    delete canvas_acceptance_1;
    delete canvas_acceptance_2;
    delete canvas_acceptance_3;
    delete canvas_acceptance_4;
    delete canvas_acceptance_5;

    // Clean up histograms
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_data_histograms_0[topo_idx]) delete h;
    }
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_data_histograms_1[topo_idx]) delete h;
    }
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_data_histograms_2[topo_idx]) delete h;
    }
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_eppi0_histograms_0[topo_idx]) delete h;
    }
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_eppi0_histograms_1[topo_idx]) delete h;
    }
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_eppi0_histograms_2[topo_idx]) delete h;
    }

    for (auto& h : h_mc_gen_histograms_0) delete h;
    for (auto& h : h_mc_gen_histograms_1) delete h;
    for (auto& h : h_mc_gen_histograms_2) delete h;
    for (auto& h : h_mc_gen_aaogen_histograms_0) delete h;
    for (auto& h : h_mc_gen_aaogen_histograms_1) delete h;
    for (auto& h : h_mc_gen_aaogen_histograms_2) delete h;
    for (auto& h : h_mc_rec_histograms_0) delete h;
    for (auto& h : h_mc_rec_histograms_1) delete h;
    for (auto& h : h_mc_rec_histograms_2) delete h;
    for (auto& h : h_mc_rec_aaogen_histograms_0) delete h;
    for (auto& h : h_mc_rec_aaogen_histograms_0) delete h;
    for (auto& h : h_mc_rec_aaogen_histograms_0) delete h;
    for (auto& h : h_acceptance_histograms_0) delete h;
    for (auto& h : h_acceptance_histograms_1) delete h;
    for (auto& h : h_acceptance_histograms_2) delete h;
    for (auto& h : h_acceptance_eppi0_histograms_0) delete h;
    for (auto& h : h_acceptance_eppi0_histograms_1) delete h;
    for (auto& h : h_acceptance_eppi0_histograms_2) delete h;

    // Reset the readers after each iteration
    data_readers[0].Restart(); data_readers[1].Restart(); data_readers[2].Restart();
    mc_gen_readers[0].Restart(); mc_gen_readers[1].Restart(); mc_gen_readers[2].Restart();
    mc_rec_readers[0].Restart(); mc_rec_readers[1].Restart(); mc_rec_readers[2].Restart();
    eppi0_readers[0].Restart(); eppi0_readers[1].Restart(); eppi0_readers[2].Restart();
    mc_gen_aaogen_readers[0].Restart(); mc_gen_aaogen_readers[1].Restart(); mc_gen_aaogen_readers[2].Restart();
    mc_rec_aaogen_readers[0].Restart(); mc_rec_aaogen_readers[1].Restart(); mc_rec_aaogen_readers[2].Restart();

    // Return all the gathered unfolding data
    return all_unfolding_data;
}