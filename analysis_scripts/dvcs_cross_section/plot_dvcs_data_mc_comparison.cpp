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

// Helper function to find the closest lower power of 10 for a given value
double closest_lower_power_of_ten(double value) {
    if (value <= 0) return 0.1;  // Safety check for non-positive values
    return std::pow(10, std::floor(std::log10(value)));
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

void plot_dvcs_data_mc_comparison(const std::string& output_dir, const std::string& analysisType, 
                                  const std::string& topology, int xB_bin,
                                  const std::vector<BinBoundary>& bin_boundaries, 
                                  TTreeReader& data_reader, 
                                  TTreeReader& mc_gen_reader, 
                                  TTreeReader& mc_rec_reader) {

    // Precompute the relevant bins for the xB_bin to avoid redundant string parsing
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);

    // Count the number of Q²-t bins for the current xB bin
    int n_Q2t_bins = relevant_bins.size();
    std::cout << "Current xB_bin = " << xB_bin << ", Number of Q2t bins: " << n_Q2t_bins << std::endl;

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    // For xB_bin 3 and 4, use a 6x5 layout (30 subplots), otherwise use next perfect square logic
    int n_subplots;
    int n_columns;
    int n_rows;

    if (xB_bin == 3 || xB_bin == 4) {
        // Special case for xB_bin 3 and 4, use 6x5 layout
        n_subplots = 30;  // Fixed to 30 bins
        n_columns = 5;    // 5 columns
        n_rows = 6;       // 6 rows
    } else {
        // Default case: use perfect square logic
        n_subplots = next_perfect_square(n_Q2t_bins);
        n_columns = std::sqrt(n_subplots);
        n_rows = n_columns;
    }

    TCanvas* canvas = new TCanvas("c1", "Data vs MC", canvas_width, canvas_height);
    canvas->Divide(n_columns, n_rows);

    gStyle->SetOptStat(0);

    std::vector<TH1D*> h_data_histograms(n_Q2t_bins);
    std::vector<TH1D*> h_mc_gen_histograms(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_histograms(n_Q2t_bins);

    // Create histograms only for the relevant bins
    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        const auto& bin = bin_boundaries[relevant_bins[idx]];

        // Update title to use the analysisType and topology
        std::string title = Form("%s, %s: x_{B} avg: %.2f, Q^{2} avg: %.2f, -t avg: %.2f", 
                                 analysisType.c_str(), 
                                 topology.c_str(), 
                                 bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));

        h_data_histograms[idx] = new TH1D(Form("h_data_%d", idx), title.c_str(), 24, 0, 360);
        h_mc_gen_histograms[idx] = new TH1D(Form("h_mc_gen_%d", idx), "gen", 24, 0, 360);
        h_mc_rec_histograms[idx] = new TH1D(Form("h_mc_rec_%d", idx), "rec", 24, 0, 360);

        // Increase axis label font size slightly
        h_data_histograms[idx]->GetXaxis()->SetLabelSize(0.05);
        h_data_histograms[idx]->GetYaxis()->SetLabelSize(0.05);

        h_mc_gen_histograms[idx]->GetXaxis()->SetLabelSize(0.05);
        h_mc_gen_histograms[idx]->GetYaxis()->SetLabelSize(0.05);

        h_mc_rec_histograms[idx]->GetXaxis()->SetLabelSize(0.05);
        h_mc_rec_histograms[idx]->GetYaxis()->SetLabelSize(0.05);

        // Increase axis title font size
        h_data_histograms[idx]->GetXaxis()->SetTitleSize(0.06);
        h_data_histograms[idx]->GetYaxis()->SetTitleSize(0.06);

        h_mc_gen_histograms[idx]->GetXaxis()->SetTitleSize(0.06);
        h_mc_gen_histograms[idx]->GetYaxis()->SetTitleSize(0.06);

        h_mc_rec_histograms[idx]->GetXaxis()->SetTitleSize(0.06);
        h_mc_rec_histograms[idx]->GetYaxis()->SetTitleSize(0.06);

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
    TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");

    // Handle theta_neutral_neutral based on analysis type (dvcs or eppi0)
    TTreeReaderValue<double>* theta_neutral_neutral_data;
    TTreeReaderValue<double>* theta_neutral_neutral_mc_gen;
    TTreeReaderValue<double>* theta_neutral_neutral_mc_rec;
    if (analysisType == "dvcs") {
        theta_neutral_neutral_data = new TTreeReaderValue<double>(data_reader, "theta_gamma_gamma");
        theta_neutral_neutral_mc_gen = new TTreeReaderValue<double>(mc_gen_reader, "theta_gamma_gamma");
        theta_neutral_neutral_mc_rec = new TTreeReaderValue<double>(mc_rec_reader, "theta_gamma_gamma");
    } else {
        theta_neutral_neutral_data = new TTreeReaderValue<double>(data_reader, "theta_pi0_pi0");
        theta_neutral_neutral_mc_gen = new TTreeReaderValue<double>(mc_gen_reader, "theta_pi0_pi0");
        theta_neutral_neutral_mc_rec = new TTreeReaderValue<double>(mc_rec_reader, "theta_pi0_pi0");
    }

    // Readers for detector status variables
    TTreeReaderValue<int> detector1_data(data_reader, "detector1");
    TTreeReaderValue<int> detector2_data(data_reader, "detector2");
    TTreeReaderValue<int> detector1_mc(mc_rec_reader, "detector1");
    TTreeReaderValue<int> detector2_mc(mc_rec_reader, "detector2");

    TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi");
    TTreeReaderValue<double> xB_mc_gen(mc_gen_reader, "x");
    TTreeReaderValue<double> Q2_mc_gen(mc_gen_reader, "Q2");
    TTreeReaderValue<double> t1_mc_gen(mc_gen_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen(mc_gen_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_mc_gen(mc_gen_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen(mc_gen_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen(mc_gen_reader, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi");
    TTreeReaderValue<double> xB_mc_rec(mc_rec_reader, "x");
    TTreeReaderValue<double> Q2_mc_rec(mc_rec_reader, "Q2");
    TTreeReaderValue<double> t1_mc_rec(mc_rec_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec(mc_rec_reader, "open_angle_ep2");
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
                ((topology == "(FD,FD)" && *detector1_data == 1 && *detector2_data == 1) ||
                 (topology == "(CD,FD)" && *detector1_data == 2 && *detector2_data == 1) ||
                 (topology == "(CD,FT)" && *detector1_data == 2 && *detector2_data == 0)) &&
                apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, **theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {

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
                ((topology == "(FD,FD)" && *detector1_mc == 1 && *detector2_mc == 1) ||
                 (topology == "(CD,FD)" && *detector1_mc == 2 && *detector2_mc == 1) ||
                 (topology == "(CD,FT)" && *detector1_mc == 2 && *detector2_mc == 0)) &&
                apply_kinematic_cuts(*t1_mc_gen, *open_angle_ep2_mc_gen, **theta_neutral_neutral_mc_gen, *Emiss2_mc_gen, *Mx2_1_mc_gen, *pTmiss_mc_gen)) {

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
                ((topology == "(FD,FD)" && *detector1_mc == 1 && *detector2_mc == 1) ||
                 (topology == "(CD,FD)" && *detector1_mc == 2 && *detector2_mc == 1) ||
                 (topology == "(CD,FT)" && *detector1_mc == 2 && *detector2_mc == 0)) &&
                apply_kinematic_cuts(*t1_mc_rec, *open_angle_ep2_mc_rec, **theta_neutral_neutral_mc_rec, *Emiss2_mc_rec, *Mx2_1_mc_rec, *pTmiss_mc_rec)) {

                h_mc_rec_histograms[idx]->Fill(phi_mc_rec_deg);
                break;
            }
        }
    }

    // Normalize histograms and plot in each subplot
    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        TPad* pad = (TPad*)canvas->cd(idx + 1);  // Get the current pad (subplot)
    
        // Set the margins for the subplot (adjust these values as needed)
        pad->SetLeftMargin(0.15);   // Add space to the left
        pad->SetBottomMargin(0.15); // Add space below

        TH1D* h_data = h_data_histograms[idx];
        TH1D* h_mc_gen = h_mc_gen_histograms[idx];
        TH1D* h_mc_rec = h_mc_rec_histograms[idx];

        // Get integrals before scaling
        double integral_data = h_data->Integral();
        double integral_mc_rec = h_mc_rec->Integral();
        double integral_mc_gen = h_mc_gen->Integral();

        // Scale reconstructed MC to have the same integral as data
        if (integral_mc_rec > 0 && integral_data > 0) {
            h_mc_rec->Scale(integral_data / integral_mc_rec);
        }

        // Now scale generated MC using the ratio of its original integral to the original reconstructed MC integral
        if (integral_mc_gen > 0 && integral_mc_rec > 0) {
            h_mc_gen->Scale((integral_mc_gen / integral_mc_rec) * (integral_data / integral_mc_rec));
        }

        // Find the minimum non-zero value in each histogram (for this subplot)
        double min_data = h_data->GetMinimum(0);  // Get minimum non-zero value
        double min_mc_gen = h_mc_gen->GetMinimum(0);  // Get minimum non-zero value
        double min_mc_rec = h_mc_rec->GetMinimum(0);  // Get minimum non-zero value

        // Find the smallest non-zero value across the three histograms
        double local_min = std::numeric_limits<double>::max();
        if (min_data > 0) local_min = std::min(local_min, min_data);
        if (min_mc_gen > 0) local_min = std::min(local_min, min_mc_gen);
        if (min_mc_rec > 0) local_min = std::min(local_min, min_mc_rec);

        // Find the closest lower power of 10 for the log scale minimum
        double min_y_log = closest_lower_power_of_ten(local_min);

        // Find the maximum bin content across all histograms for plotting purposes
        double max_value = std::max({h_data->GetMaximum(), h_mc_gen->GetMaximum(), h_mc_rec->GetMaximum()});

        // Set the maximum for the histograms
        h_data->SetMaximum(1.35 * max_value);

        // Set the log scale minimum value for this subplot
        h_data->SetMinimum(min_y_log);
        h_mc_gen->SetMinimum(min_y_log);
        h_mc_rec->SetMinimum(min_y_log);

        // Set log scale for the canvas
        canvas->cd(idx + 1)->SetLogy();

        // Set colors and styles for plotting
        h_data->SetMarkerColor(kBlue);
        h_data->SetMarkerStyle(20);
        h_data->SetLineColor(kBlue);

        h_mc_gen->SetLineColor(kBlack);
        h_mc_gen->SetLineStyle(2);

        h_mc_rec->SetMarkerColor(kRed);
        h_mc_rec->SetMarkerStyle(22);
        h_mc_rec->SetLineColor(kRed);

        // Draw histograms
        h_data->Draw("E1");
        h_mc_rec->Draw("E1 SAME");
        h_mc_gen->Draw("HIST SAME");

        // Add legend
        TLegend* legend = new TLegend(0.575, 0.45, 0.9, 0.75);
        legend->AddEntry(h_data, "Data", "lep");
        legend->AddEntry(h_mc_rec, "Reconstructed MC", "lep");
        legend->AddEntry(h_mc_gen, "Generated MC", "l");
        legend->SetTextSize(0.04);
        legend->Draw();
    }

    // Save canvas to the output directory with the analysisType and topology in the filename
    std::string filename = output_dir + "/phi_data_mc_comparison_" + analysisType + "_" + topology + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    canvas->SaveAs(filename.c_str());

    // Clean up histograms and canvas
    for (auto& h : h_data_histograms) delete h;
    for (auto& h : h_mc_gen_histograms) delete h;
    for (auto& h : h_mc_rec_histograms) delete h;
    delete canvas;

    // Clean up dynamically allocated memory
    delete theta_neutral_neutral_data;
    delete theta_neutral_neutral_mc;
}