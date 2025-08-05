#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLatex.h>
#include <cmath>
#include <string>
#include <vector>
#include "all_bin_data.h" // To access AllBinData structure

// Helper function to precompute relevant bins
std::vector<int> precompute_relevant_bins(int xB_bin, const std::vector<AllBinData>& all_bin_data) {
    std::vector<int> relevant_bins;
    for (size_t bin_idx = 0; bin_idx < all_bin_data.size(); ++bin_idx) {
        if (all_bin_data[bin_idx].bin_number == xB_bin) {
            relevant_bins.push_back(bin_idx);
        }
    }
    return relevant_bins;
}

// Function to plot yield comparison
void plot_yield_comparison(const std::string& output_dir, int xB_bin, const std::vector<AllBinData>& all_bin_data) {

    // Precompute relevant bins for the xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, all_bin_data);
    int n_Q2t_bins = relevant_bins.size();
    std::cout << "Current xB_bin = " << xB_bin << ", Number of Q2t bins: " << n_Q2t_bins << std::endl;

    // Set canvas size and layout (using the same layout logic as before)
    TCanvas* canvas = new TCanvas("c1", "Yield Comparison", 1200, 800);
    int n_columns = 5; // Adjust as necessary
    int n_rows = std::ceil(static_cast<double>(n_Q2t_bins) / n_columns);
    canvas->Divide(n_columns, n_rows);

    // Style settings
    gStyle->SetOptStat(0);

    // Loop over all relevant bins and create a histogram for each bin
    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        const AllBinData& bin = all_bin_data[relevant_bins[idx]];

        // Set up the bin title
        std::string title = Form("xB: %.4f, Q^2: %.4f, t: %.4f, Phi [%.4f, %.4f]",
                                 bin.xB_avg, bin.Q2_avg, bin.t_avg, bin.phi_min, bin.phi_max);

        // Create histograms
        TH1D* h_my_data = new TH1D(Form("h_my_data_%d", idx), title.c_str(), 24, bin.phi_min, bin.phi_max);
        TH1D* h_colleague_data = new TH1D(Form("h_colleague_data_%d", idx), title.c_str(), 24, bin.phi_min, bin.phi_max);

        // Fill histograms with my data and my colleague's data (replace with actual data)
        h_my_data->SetBinContent(1, bin.yield_epg_FD_FD_inb);    // Example for your data (you can loop or handle as needed)
        h_colleague_data->SetBinContent(1, bin.yield_eppi0_FD_FD_inb); // Example for colleague's data

        // Customize histogram aesthetics
        h_my_data->SetMarkerColor(kBlue);
        h_my_data->SetMarkerStyle(20);
        h_my_data->SetLineColor(kBlue);

        h_colleague_data->SetMarkerColor(kRed);
        h_colleague_data->SetMarkerStyle(22);
        h_colleague_data->SetLineColor(kRed);

        // Plot histograms
        canvas->cd(idx + 1);  // Go to the next subplot
        h_my_data->Draw("E1");
        h_colleague_data->Draw("E1 SAME");

        // Add a legend
        TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend->AddEntry(h_my_data, "My Data", "lep");
        legend->AddEntry(h_colleague_data, "Colleague Data", "lep");
        legend->Draw();
    }

    // Save canvas to the output directory
    std::string filename = output_dir + "/yield_comparison_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    canvas->SaveAs(filename.c_str());

    // Clean up memory
    delete canvas;
}