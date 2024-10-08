#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLatex.h>  // For adding text in each subplot
#include <cmath>     // For conversion from radians to degrees
#include <string>
#include <vector>
#include "bin_boundaries.h"  // Include the header where BinBoundary is defined
#include <algorithm> // For remove_if
#include <cctype> // For isspace
#include "kinematic_cuts.h"

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Helper function to remove spaces and parentheses
std::string clean_bin_label(const std::string& label) {
    std::string clean_label = label;
    clean_label.erase(std::remove_if(clean_label.begin(), clean_label.end(), [](unsigned char c) {
        return std::isspace(c) || c == '(' || c == ')';
    }), clean_label.end());
    return clean_label;
}

// Find the next perfect square greater than or equal to the given number
int next_perfect_square(int n) {
    int square_root = std::ceil(std::sqrt(n));
    return square_root * square_root;
}

// Precompute the relevant bin indices based on xB_bin
std::vector<int> precompute_relevant_bins(int xB_bin, const std::vector<BinBoundary>& bin_boundaries) {
    std::vector<int> relevant_bins;
    for (size_t bin_idx = 0; bin_idx < bin_boundaries.size(); ++bin_idx) {
        std::string bin_label = clean_bin_label(bin_boundaries[bin_idx].bin_label);
        size_t first_comma = bin_label.find(',');

        if (first_comma != std::string::npos) {
            int xB_label = std::stoi(bin_label.substr(0, first_comma));
            if (xB_label == xB_bin) {
                relevant_bins.push_back(bin_idx);
            }
        }
    }
    return relevant_bins;
}

// Plot function for DVCS data/MC comparison
void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, const std::vector<BinBoundary>& bin_boundaries, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader) {

    // Precompute the relevant bins for the xB_bin to avoid redundant string parsing
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);

    // Count the number of Q²-t bins for the current xB bin
    int n_Q2t_bins = relevant_bins.size();
    std::cout << "Current xB_bin = " << xB_bin << ", Number of Q2t bins: " << n_Q2t_bins << std::endl;

    int n_subplots = next_perfect_square(n_Q2t_bins);
    int n_columns = std::sqrt(n_subplots);
    int n_rows = n_columns;

    TCanvas* canvas = new TCanvas("c1", "Data vs MC", 1200, 800);
    canvas->Divide(n_columns, n_rows);

    gStyle->SetOptStat(0);

    std::vector<TH1D*> h_data_histograms(n_Q2t_bins);
    std::vector<TH1D*> h_mc_gen_histograms(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_histograms(n_Q2t_bins);

    // Create histograms only for the relevant bins
    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        const auto& bin = bin_boundaries[relevant_bins[idx]];
        std::string title = Form("x_{B}: %.2f-%.2f, Q^{2}: %.2f-%.2f, |t|: %.2f-%.2f",
                                 bin.xB_low, bin.xB_high, bin.Q2_low, bin.Q2_high, std::abs(bin.t_low), std::abs(bin.t_high));

        h_data_histograms[idx] = new TH1D(Form("h_data_%d", idx), title.c_str(), 24, 0, 360);
        h_mc_gen_histograms[idx] = new TH1D(Form("h_mc_gen_%d", idx), "gen", 24, 0, 360);
        h_mc_rec_histograms[idx] = new TH1D(Form("h_mc_rec_%d", idx), "rec", 24, 0, 360);

        h_data_histograms[idx]->GetXaxis()->SetTitle("#phi");
        h_data_histograms[idx]->GetYaxis()->SetTitle("Normalized Counts");

        h_mc_gen_histograms[idx]->GetXaxis()->SetTitle("#phi");
        h_mc_gen_histograms[idx]->GetYaxis()->SetTitle("Normalized Counts");

        h_mc_rec_histograms[idx]->GetXaxis()->SetTitle("#phi");
        h_mc_rec_histograms[idx]->GetYaxis()->SetTitle("Normalized Counts");
    }

    std::cout << "There are " << n_Q2t_bins << " bins created." << std::endl;

    // Restart the readers before looping over data
    data_reader.Restart();
    mc_gen_reader.Restart();
    mc_rec_reader.Restart();

    // Readers for data and MC
    TTreeReaderValue<double> phi_data(data_reader, "phi");
    TTreeReaderValue<double> xB_data(data_reader, "x");
    TTreeReaderValue<double> Q2_data(data_reader, "Q2");
    TTreeReaderValue<double> t1_data(data_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_data(data_reader, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_data(data_reader, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");

    TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi");
    TTreeReaderValue<double> xB_mc_gen(mc_gen_reader, "x");
    TTreeReaderValue<double> Q2_mc_gen(mc_gen_reader, "Q2");
    TTreeReaderValue<double> t1_mc_gen(mc_gen_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen(mc_gen_reader, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen(mc_gen_reader, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_mc_gen(mc_gen_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen(mc_gen_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen(mc_gen_reader, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi");
    TTreeReaderValue<double> xB_mc_rec(mc_rec_reader, "x");
    TTreeReaderValue<double> Q2_mc_rec(mc_rec_reader, "Q2");
    TTreeReaderValue<double> t1_mc_rec(mc_rec_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec(mc_rec_reader, "open_angle_ep2");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec(mc_rec_reader, "theta_gamma_gamma");
    TTreeReaderValue<double> Emiss2_mc_rec(mc_rec_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec(mc_rec_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec(mc_rec_reader, "pTmiss");

    // Data histograms filling
    std::cout << "Started data " << std::endl;
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if (*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high &&
                apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {

                h_data_histograms[idx]->Fill(phi_deg);
                break;
            }
        }
    }

    // MC-generated histograms filling
    std::cout << "Started mc gen " << std::endl;
    while (mc_gen_reader.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen * RAD_TO_DEG;
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if (*xB_mc_gen >= bin.xB_low && *xB_mc_gen <= bin.xB_high &&
                *Q2_mc_gen >= bin.Q2_low && *Q2_mc_gen <= bin.Q2_high &&
                std::abs(*t1_mc_gen) >= bin.t_low && std::abs(*t1_mc_gen) <= bin.t_high &&
                apply_kinematic_cuts(*t1_mc_gen, *open_angle_ep2_mc_gen, *theta_neutral_neutral_mc_gen, *Emiss2_mc_gen, *Mx2_1_mc_gen, *pTmiss_mc_gen)) {

                h_mc_gen_histograms[idx]->Fill(phi_mc_gen_deg);
                break;
            }
        }
    }

    // MC-reconstructed histograms filling
    std::cout << "Started mc rec " << std::endl;
    while (mc_rec_reader.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec * RAD_TO_DEG;
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if (*xB_mc_rec >= bin.xB_low && *xB_mc_rec <= bin.xB_high &&
                *Q2_mc_rec >= bin.Q2_low && *Q2_mc_rec <= bin.Q2_high &&
                std::abs(*t1_mc_rec) >= bin.t_low && std::abs(*t1_mc_rec) <= bin.t_high &&
                apply_kinematic_cuts(*t1_mc_rec, *open_angle_ep2_mc_rec, *theta_neutral_neutral_mc_rec, *Emiss2_mc_rec, *Mx2_1_mc_rec, *pTmiss_mc_rec)) {

                h_mc_rec_histograms[idx]->Fill(phi_mc_rec_deg);
                break;
            }
        }
    }

    // Normalize histograms and plot in each subplot
    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        canvas->cd(idx + 1);

        TH1D* h_data = h_data_histograms[idx];
        TH1D* h_mc_gen = h_mc_gen_histograms[idx];
        TH1D* h_mc_rec = h_mc_rec_histograms[idx];

        if (h_data->Integral() > 0) h_data->Scale(1.0 / h_data->Integral());
        if (h_mc_gen->Integral() > 0) h_mc_gen->Scale(1.0 / h_mc_gen->Integral());
        if (h_mc_rec->Integral() > 0) h_mc_rec->Scale(1.0 / h_mc_rec->Integral());

        double max_value = std::max({h_data->GetMaximum(), h_mc_gen->GetMaximum(), h_mc_rec->GetMaximum()});

        h_data->SetMaximum(1.35 * max_value);
        h_data->SetMinimum(0);

        h_data->SetMarkerColor(kBlue);
        h_data->SetMarkerStyle(20);
        h_data->SetLineColor(kBlue);

        h_mc_gen->SetLineColor(kBlack);
        h_mc_gen->SetLineStyle(2);

        h_mc_rec->SetMarkerColor(kRed);
        h_mc_rec->SetMarkerStyle(22);
        h_mc_rec->SetLineColor(kRed);

        h_data->Draw("E1");
        h_mc_rec->Draw("E1 SAME");
        h_mc_gen->Draw("HIST SAME");

        TLegend* legend = new TLegend(0.575, 0.6, 0.9, 0.9);
        legend->AddEntry(h_data, "Data", "lep");
        legend->AddEntry(h_mc_rec, "Reconstructed MC", "lep");
        legend->AddEntry(h_mc_gen, "Generated MC", "l");
        legend->SetTextSize(0.04);
        legend->Draw();
    }

    // Save canvas to the output directory
    std::string filename = output_dir + "/phi_data_mc_comparison_xB_bin_" + std::to_string(xB_bin) + ".png";
    canvas->SaveAs(filename.c_str());

    // Clean up histograms and canvas
    for (auto& h : h_data_histograms) delete h;
    for (auto& h : h_mc_gen_histograms) delete h;
    for (auto& h : h_mc_rec_histograms) delete h;
    delete canvas;
}