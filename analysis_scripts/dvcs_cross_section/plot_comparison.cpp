#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <filesystem>

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TSystem.h>

#include "bin_boundaries.h"
#include "unfolding_data.h"
#include "plot_comparison.h"

struct PreviousUnfoldingData {
    int bin_number;
    double xB_min, xB_max, xB_avg;
    double Q2_min, Q2_max, Q2_avg;
    double t_min, t_max, t_avg;
    double phi_min, phi_max, phi_avg;
    double unfolded_yield_dvcs_Fa18Inb;
    double unfolded_yield_dvcs_Fa18Out;
    double unfolded_yield_eppi0_Fa18Inb;
    double unfolded_yield_eppi0_Fa18Out;
};

// Function to read previous unfolding data from CSV
// Function to read unfolding data from CSV
std::vector<UnfoldingData> read_unfolding_data(const std::string& filename, bool is_previous_data) {
    std::vector<UnfoldingData> data_vector;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open data file " << filename << std::endl;
        return data_vector;
    }

    std::string line;
    // Read the header line
    if (!std::getline(file, line)) {
        std::cerr << "Error: Could not read header line from file " << filename << std::endl;
        return data_vector;
    }

    // Parse the header to get column indices
    std::istringstream header_stream(line);
    std::string header_cell;
    std::vector<std::string> headers;
    while (std::getline(header_stream, header_cell, ',')) {
        headers.push_back(header_cell);
    }

    // Map column names to indices
    std::map<std::string, int> column_indices;
    for (size_t i = 0; i < headers.size(); ++i) {
        column_indices[headers[i]] = i;
    }

    // Now, set up the column names based on whether we're reading previous data or new data
    int idx_bin_number;
    int idx_xBmin, idx_xBmax, idx_xBavg;
    int idx_Q2min, idx_Q2max, idx_Q2avg;
    int idx_t_abs_min, idx_t_abs_max, idx_t_abs_avg;
    int idx_phimin, idx_phimax, idx_phiavg;

    int idx_unfolded_yield_dvcs_Fa18Inb;
    int idx_unfolded_yield_dvcs_Fa18Out;
    int idx_unfolded_yield_eppi0_Fa18Inb;
    int idx_unfolded_yield_eppi0_Fa18Out;

    if (is_previous_data) {
        // For the previous data CSV, use the column names from the old CSV
        idx_bin_number = column_indices["Bin Name"];
        idx_xBmin = column_indices["xBmin"];
        idx_xBmax = column_indices["xBmax"];
        idx_xBavg = column_indices["xBavg"];
        idx_Q2min = column_indices["Q2min"];
        idx_Q2max = column_indices["Q2max"];
        idx_Q2avg = column_indices["Q2avg"];
        idx_t_abs_min = column_indices["t_abs_min"];
        idx_t_abs_max = column_indices["t_abs_max"];
        idx_t_abs_avg = column_indices["t_abs_avg"];
        idx_phimin = column_indices["phimin"];
        idx_phimax = column_indices["phimax"];
        idx_phiavg = column_indices["phiavg"];

        // Unfolded yields for DVCS
        idx_unfolded_yield_dvcs_Fa18Inb = column_indices["acceptance corrected yield, ep->epg, exp"];
        idx_unfolded_yield_dvcs_Fa18Out = column_indices["acceptance corrected yield, ep->epg, exp, outbending"]; // Adjust if necessary

        // For eppi0 unfolded yields, we might need to adjust the column names based on what's available
        idx_unfolded_yield_eppi0_Fa18Inb = column_indices["acceptance corrected yield, ep->eppi0, exp"]; // Adjust if necessary
        idx_unfolded_yield_eppi0_Fa18Out = column_indices["acceptance corrected yield, ep->eppi0, exp, outbending"]; // Adjust if necessary
    } else {
        // For the new data CSV, use the column names from the new CSV
        idx_bin_number = column_indices["Bin"];
        idx_xBmin = column_indices["xB_min"];
        idx_xBmax = column_indices["xB_max"];
        idx_xBavg = column_indices["xB_avg"];
        idx_Q2min = column_indices["Q2_min"];
        idx_Q2max = column_indices["Q2_max"];
        idx_Q2avg = column_indices["Q2_avg"];
        idx_t_abs_min = column_indices["t_min"];
        idx_t_abs_max = column_indices["t_max"];
        idx_t_abs_avg = column_indices["t_avg"];
        idx_phimin = column_indices["phi_min"];
        idx_phimax = column_indices["phi_max"];
        idx_phiavg = column_indices["phi_avg"];

        // Unfolded yields for DVCS
        idx_unfolded_yield_dvcs_Fa18Inb = column_indices["ep->e'pgamma unfolded_yield_Fa18Inb"];
        idx_unfolded_yield_dvcs_Fa18Out = column_indices["ep->e'pgamma unfolded_yield_Fa18Out"];

        // Unfolded yields for eppi0
        idx_unfolded_yield_eppi0_Fa18Inb = column_indices["ep->e'ppi0 unfolded_yield_Fa18Inb"];
        idx_unfolded_yield_eppi0_Fa18Out = column_indices["ep->e'ppi0 unfolded_yield_Fa18Out"];
    }

    // Read the data lines
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        std::istringstream line_stream(line);
        std::string cell;
        std::vector<std::string> cells;
        while (std::getline(line_stream, cell, ',')) {
            cells.push_back(cell);
        }

        if (cells.size() < headers.size()) {
            // Handle the case where some cells might be empty
            cells.resize(headers.size());
        }

        UnfoldingData data;

        data.bin_number = std::stoi(cells[idx_bin_number]);
        data.xB_min = std::stod(cells[idx_xBmin]);
        data.xB_max = std::stod(cells[idx_xBmax]);
        data.xB_avg = std::stod(cells[idx_xBavg]);
        data.Q2_min = std::stod(cells[idx_Q2min]);
        data.Q2_max = std::stod(cells[idx_Q2max]);
        data.Q2_avg = std::stod(cells[idx_Q2avg]);
        data.t_min = std::stod(cells[idx_t_abs_min]);
        data.t_max = std::stod(cells[idx_t_abs_max]);
        data.t_avg = std::stod(cells[idx_t_abs_avg]);
        data.phi_min = std::stod(cells[idx_phimin]);
        data.phi_max = std::stod(cells[idx_phimax]);
        data.phi_avg = std::stod(cells[idx_phiavg]);

        // Unfolded yields
        // Handle possible missing data gracefully
        if (cells[idx_unfolded_yield_dvcs_Fa18Inb].empty()) {
            data.unfolded_yield_dvcs_Fa18Inb = 0.0;
        } else {
            data.unfolded_yield_dvcs_Fa18Inb = std::stod(cells[idx_unfolded_yield_dvcs_Fa18Inb]);
        }

        if (cells[idx_unfolded_yield_dvcs_Fa18Out].empty()) {
            data.unfolded_yield_dvcs_Fa18Out = 0.0;
        } else {
            data.unfolded_yield_dvcs_Fa18Out = std::stod(cells[idx_unfolded_yield_dvcs_Fa18Out]);
        }

        if (cells[idx_unfolded_yield_eppi0_Fa18Inb].empty()) {
            data.unfolded_yield_eppi0_Fa18Inb = 0.0;
        } else {
            data.unfolded_yield_eppi0_Fa18Inb = std::stod(cells[idx_unfolded_yield_eppi0_Fa18Inb]);
        }

        if (cells[idx_unfolded_yield_eppi0_Fa18Out].empty()) {
            data.unfolded_yield_eppi0_Fa18Out = 0.0;
        } else {
            data.unfolded_yield_eppi0_Fa18Out = std::stod(cells[idx_unfolded_yield_eppi0_Fa18Out]);
        }

        data_vector.push_back(data);
    }

    return data_vector;
}

void plot_comparison(
    const std::string& output_dir,
    const std::vector<BinBoundary>& bin_boundaries,
    const std::string& previous_data_csv,
    const std::map<std::string, std::vector<UnfoldingData>>& all_unfolding_data
) {
    // Read previous data
    std::vector<PreviousUnfoldingData> previous_data = read_previous_unfolding_data(previous_data_csv);
    if (previous_data.empty()) {
        std::cerr << "Error: No previous data read from file." << std::endl;
        return;
    }

    // Names of the periods
    std::vector<std::string> period_names = {"Fa18Inb", "Fa18Out"};
    const int n_periods = 2; // Only Fa18Inb and Fa18Out

    // Analysis types
    std::vector<std::string> analysis_types = {"dvcs", "eppi0"};

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    gStyle->SetOptStat(0);
    std::cout << "before getting all unfolding data" << std::endl;

    // Debug: print keys in all_unfolding_data
    std::cout << "Keys in all_unfolding_data: ";
    for (const auto& kv : all_unfolding_data) {
        std::cout << kv.first << " ";
    }
    std::cout << std::endl;

    // Attempt to access "combined" key
    const auto& combined_unfolding_data = all_unfolding_data.at("combined");
    std::cout << "after getting all unfolding data" << std::endl;
    // Create a map to quickly access previous data by bin identifiers
    // We'll use a combination of bin_number and phi_avg as the key
    std::map<std::pair<int, double>, PreviousUnfoldingData> previous_data_map;
    for (const auto& prev_data : previous_data) {
        previous_data_map[{prev_data.bin_number, prev_data.phi_avg}] = prev_data;
    }

    // Group the data by xB bins
    std::map<int, std::vector<size_t>> xB_bin_indices; // Map from xB bin number to vector of indices

    for (size_t idx = 0; idx < combined_unfolding_data.size(); ++idx) {
        int xB_bin = combined_unfolding_data[idx].bin_number;
        xB_bin_indices[xB_bin].push_back(idx);
    }

    // For each xB bin
    for (const auto& [xB_bin, indices] : xB_bin_indices) {
        // For each period (Fa18Inb, Fa18Out)
        for (int period = 0; period < n_periods; ++period) {
            std::string period_name = period_names[period];

            // For each analysis type (dvcs, eppi0)
            for (const auto& analysis_type : analysis_types) {
                // Create output directories
                std::string output_subdir = output_dir + "/cross_check/" + analysis_type + "/" + period_name;
                gSystem->mkdir(output_subdir.c_str(), true);

                size_t n_bins = indices.size();
                int n_columns = std::ceil(std::sqrt(n_bins));
                int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

                // Create the canvas
                TCanvas* canvas = new TCanvas(Form("c_comparison_%s_%s_xBbin_%d", analysis_type.c_str(), period_name.c_str(), xB_bin),
                                              "Comparison Plots", canvas_width, canvas_height);
                canvas->Divide(n_columns, n_rows);

                for (size_t idx = 0; idx < indices.size(); ++idx) {
                    size_t data_idx = indices[idx];
                    const UnfoldingData& unfolding_data = combined_unfolding_data[data_idx];
                    size_t n_phi_bins = unfolding_data.phi_min.size();

                    // Prepare data for plotting
                    std::vector<double> phi_centers(n_phi_bins);
                    std::vector<double> phi_widths(n_phi_bins);
                    std::vector<double> unfolded_counts_new(n_phi_bins);
                    std::vector<double> unfolded_errors_new(n_phi_bins);
                    std::vector<double> unfolded_counts_prev(n_phi_bins);

                    for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                        phi_centers[phi_idx] = unfolding_data.phi_avg[phi_idx];
                        phi_widths[phi_idx] = (unfolding_data.phi_max[phi_idx] - unfolding_data.phi_min[phi_idx]) / 2.0;

                        // New data
                        if (analysis_type == "dvcs") {
                            unfolded_counts_new[phi_idx] = unfolding_data.unfolded_yields[period][phi_idx];
                            double raw_yield = unfolding_data.combined_raw_yields[period][phi_idx];
                            unfolded_errors_new[phi_idx] = (unfolding_data.acceptance[period][phi_idx] > 0)
                                                               ? std::sqrt(raw_yield) / unfolding_data.acceptance[period][phi_idx]
                                                               : 0.0;
                        } else if (analysis_type == "eppi0") {
                            int eppi0_period = period + 3; // Adjust period index for eppi0 data
                            unfolded_counts_new[phi_idx] = unfolding_data.unfolded_yields[eppi0_period][phi_idx];
                            double raw_yield = unfolding_data.combined_raw_yields[eppi0_period][phi_idx];
                            unfolded_errors_new[phi_idx] = (unfolding_data.acceptance[eppi0_period][phi_idx] > 0)
                                                               ? std::sqrt(raw_yield) / unfolding_data.acceptance[eppi0_period][phi_idx]
                                                               : 0.0;
                        }

                        // Previous data
                        auto key = std::make_pair(unfolding_data.bin_number, phi_centers[phi_idx]);
                        if (previous_data_map.find(key) != previous_data_map.end()) {
                            const PreviousUnfoldingData& prev_data = previous_data_map[key];
                            if (analysis_type == "dvcs") {
                                if (period == 0) { // Fa18Inb
                                    unfolded_counts_prev[phi_idx] = prev_data.unfolded_yield_dvcs_Fa18Inb;
                                } else if (period == 1) { // Fa18Out
                                    unfolded_counts_prev[phi_idx] = prev_data.unfolded_yield_dvcs_Fa18Out;
                                }
                            } else if (analysis_type == "eppi0") {
                                if (period == 0) { // Fa18Inb
                                    unfolded_counts_prev[phi_idx] = prev_data.unfolded_yield_eppi0_Fa18Inb;
                                } else if (period == 1) { // Fa18Out
                                    unfolded_counts_prev[phi_idx] = prev_data.unfolded_yield_eppi0_Fa18Out;
                                }
                            }
                        } else {
                            unfolded_counts_prev[phi_idx] = 0.0; // No matching bin in previous data
                        }
                    }

                    // Move to the appropriate pad
                    canvas->cd(idx + 1);

                    // Create TGraphErrors for new data
                    TGraphErrors* graph_new = new TGraphErrors(n_phi_bins, &phi_centers[0], &unfolded_counts_new[0], &phi_widths[0], &unfolded_errors_new[0]);

                    // Create TGraph for previous data (assuming no errors)
                    TGraph* graph_prev = new TGraph(n_phi_bins, &phi_centers[0], &unfolded_counts_prev[0]);

                    // Set styles
                    graph_new->SetMarkerColor(kBlue);
                    graph_new->SetMarkerStyle(20);
                    graph_new->SetLineColor(kBlue);

                    graph_prev->SetMarkerColor(kRed);
                    graph_prev->SetMarkerStyle(24);
                    graph_prev->SetLineColor(kRed);

                    // Draw graphs
                    TPad* pad = (TPad*)gPad;
                    pad->SetLeftMargin(0.15);
                    pad->SetBottomMargin(0.15);

                    // Create a frame histogram for axes
                    double phi_min = 0.0;
                    double phi_max = 360.0;
                    double max_new = *std::max_element(unfolded_counts_new.begin(), unfolded_counts_new.end());
                    double max_prev = *std::max_element(unfolded_counts_prev.begin(), unfolded_counts_prev.end());
                    double count_max = std::max(max_new, max_prev) * 1.2;
                    double count_min = 0.0;

                    TH1F* frame = pad->DrawFrame(phi_min, count_min, phi_max, count_max);
                    frame->GetXaxis()->SetTitle("#phi [deg]");
                    frame->GetYaxis()->SetTitle("Unfolded Counts");

                    graph_new->Draw("P SAME");
                    graph_prev->Draw("P SAME");

                    // Add legend
                    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
                    legend->AddEntry(graph_new, "Hayward, pass-2", "p");
                    legend->AddEntry(graph_prev, "Lee, pass-1", "p");
                    legend->Draw();

                    // Add title with averages
                    std::string title = Form("%s, %s, x_{B}=%.3f, Q^{2}=%.3f, -t=%.3f",
                                             analysis_type.c_str(),
                                             period_name.c_str(),
                                             unfolding_data.xB_avg,
                                             unfolding_data.Q2_avg,
                                             unfolding_data.t_avg);
                    TLatex latex;
                    latex.SetNDC();
                    latex.SetTextSize(0.04);
                    latex.SetTextAlign(22); // Center alignment
                    latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5

                    // Clean up graphs
                    // delete graph_new;
                    // delete graph_prev;
                }

                // Save the canvas
                std::string filename = output_subdir + "/comparison_" + analysis_type + "_" + period_name +
                                       "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
                canvas->SaveAs(filename.c_str());

                // Clean up
                delete canvas;
            }
        }
    }
}