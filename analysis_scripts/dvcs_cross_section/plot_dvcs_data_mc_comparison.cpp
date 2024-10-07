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

int count_Q2t_bins_for_xB(int xB_bin, const std::vector<BinBoundary>& bin_boundaries) {
    int n_Q2t_bins = 0;

    for (const auto& bin : bin_boundaries) {
        // Clean the bin label to remove spaces and parentheses
        std::string bin_label = clean_bin_label(bin.bin_label);

        // Parse the xB_bin part of the label from the cleaned bin_label
        try {
            // Now split the label by commas and extract the first part
            size_t first_comma = bin_label.find(',');
            if (first_comma != std::string::npos) {
                int xB_label = std::stoi(bin_label.substr(0, first_comma)); // Extract xB_bin

                // Check if the extracted xB_bin matches the current xB_bin
                if (xB_label == xB_bin) {
                    n_Q2t_bins++;
                }
            } else {
                std::cerr << "Error: Unexpected bin_label format: " << bin_label << std::endl;
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error parsing xB_bin from label: '" << bin_label << "' -> " << e.what() << std::endl;
        }
    }

    return n_Q2t_bins;
}

// Plot function for DVCS data/MC comparison
void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, const std::vector<BinBoundary>& bin_boundaries, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader) {

    // Count the number of Q²-t bins for the current xB bin
    int n_Q2t_bins = count_Q2t_bins_for_xB(xB_bin, bin_boundaries);
    std::cout << "Current xB_bin = " << xB_bin << ", Number of Q2t bins: " << n_Q2t_bins << std::endl;

    // Determine the next largest perfect square to create square canvas
    int n_subplots = next_perfect_square(n_Q2t_bins);
    int n_columns = std::sqrt(n_subplots);  // Number of columns for a square layout
    int n_rows = n_columns;

    // Create canvas with dynamic subdivision
    TCanvas* canvas = new TCanvas("c1", "Data vs MC", 1200, 800);
    canvas->Divide(n_columns, n_rows);  // Dynamically divide the canvas

    // Disable stat boxes globally
    gStyle->SetOptStat(0);

    // Create histograms for each Q²-t bin (only for the required number of bins)
    std::vector<TH1D*> h_data_histograms(n_Q2t_bins);
    std::vector<TH1D*> h_mc_gen_histograms(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_histograms(n_Q2t_bins);

    // Initialize histograms
    int histogram_idx = 0;
    for (const auto& bin : bin_boundaries) {
        std::string bin_label = clean_bin_label(bin.bin_label);
        size_t first_comma = bin_label.find(',');

        if (first_comma != std::string::npos) {
            int xB_label = std::stoi(bin_label.substr(0, first_comma)); // Extract xB_bin from bin label
            if (xB_label == xB_bin) {
                if (histogram_idx >= n_Q2t_bins) break;

                // Create a title string based on the kinematic constraints
                std::string title = Form("x_{B}: %.2f-%.2f, Q^{2}: %.2f-%.2f, |t|: %.2f-%.2f",
                                         bin.xB_low, bin.xB_high,
                                         bin.Q2_low, bin.Q2_high,
                                         std::abs(bin.t_low), std::abs(bin.t_high));

                h_data_histograms[histogram_idx] = new TH1D(Form("h_data_%d", histogram_idx), title.c_str(), 24, 0, 360);
                h_mc_gen_histograms[histogram_idx] = new TH1D(Form("h_mc_gen_%d", histogram_idx), "gen", 24, 0, 360);
                h_mc_rec_histograms[histogram_idx] = new TH1D(Form("h_mc_rec_%d", histogram_idx), "rec", 24, 0, 360);

                // Set axis labels
                h_data_histograms[histogram_idx]->GetXaxis()->SetTitle("#phi");
                h_data_histograms[histogram_idx]->GetYaxis()->SetTitle("Normalized Counts");

                h_mc_gen_histograms[histogram_idx]->GetXaxis()->SetTitle("#phi");
                h_mc_gen_histograms[histogram_idx]->GetYaxis()->SetTitle("Normalized Counts");

                h_mc_rec_histograms[histogram_idx]->GetXaxis()->SetTitle("#phi");
                h_mc_rec_histograms[histogram_idx]->GetYaxis()->SetTitle("Normalized Counts");

                histogram_idx++;
            }
        }
    }

    std::cout << "There are " << histogram_idx << " bins created." << std::endl;

    // Now, fill the histograms by looping over the entries in the tree

    // Fill the data histograms
    data_reader.Restart();
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;  // Convert phi to degrees

        histogram_idx = 0;  // Reset the histogram index for each entry
        for (const auto& bin : bin_boundaries) {
            std::string bin_label = clean_bin_label(bin.bin_label);
            size_t first_comma = bin_label.find(',');

            if (first_comma != std::string::npos) {
                int xB_label = std::stoi(bin_label.substr(0, first_comma)); // Extract xB_bin from bin label

                if (xB_label == xB_bin) {
                    if (*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                        *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                        std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high) {
                        h_data_histograms[histogram_idx]->Fill(phi_deg);
                    }
                    histogram_idx++;
                }
            }
        }
    }

    // Fill the MC-generated histograms
    mc_gen_reader.Restart();
    while (mc_gen_reader.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen * RAD_TO_DEG;

        histogram_idx = 0;
        for (const auto& bin : bin_boundaries) {
            std::string bin_label = clean_bin_label(bin.bin_label);
            size_t first_comma = bin_label.find(',');

            if (first_comma != std::string::npos) {
                int xB_label = std::stoi(bin_label.substr(0, first_comma));

                if (xB_label == xB_bin) {
                    h_mc_gen_histograms[histogram_idx]->Fill(phi_mc_gen_deg);
                    histogram_idx++;
                }
            }
        }
    }

    // Fill the MC-reconstructed histograms
    mc_rec_reader.Restart();
    while (mc_rec_reader.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec * RAD_TO_DEG;

        histogram_idx = 0;
        for (const auto& bin : bin_boundaries) {
            std::string bin_label = clean_bin_label(bin.bin_label);
            size_t first_comma = bin_label.find(',');

            if (first_comma != std::string::npos) {
                int xB_label = std::stoi(bin_label.substr(0, first_comma));

                if (xB_label == xB_bin) {
                    h_mc_rec_histograms[histogram_idx]->Fill(phi_mc_rec_deg);
                    histogram_idx++;
                }
            }
        }
    }

    // Normalize histograms and plot in each subplot
    histogram_idx = 0;
    for (int subplot_idx = 1; subplot_idx <= n_Q2t_bins; ++subplot_idx) {
        canvas->cd(subplot_idx);

        TH1D* h_data = h_data_histograms[histogram_idx];
        TH1D* h_mc_gen = h_mc_gen_histograms[histogram_idx];
        TH1D* h_mc_rec = h_mc_rec_histograms[histogram_idx];

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

        histogram_idx++;
    }

    std::string filename = output_dir + "/phi_data_mc_comparison_xB_bin_" + std::to_string(xB_bin) + ".png";
    canvas->SaveAs(filename.c_str());

    for (auto& h : h_data_histograms) delete h;
    for (auto& h : h_mc_gen_histograms) delete h;
    for (auto& h : h_mc_rec_histograms) delete h;
    delete canvas;
}