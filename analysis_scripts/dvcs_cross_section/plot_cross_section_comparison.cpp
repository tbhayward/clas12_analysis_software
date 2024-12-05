// plot_cross_section_comparison.cpp

#include "plot_cross_section_comparison.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <tuple>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>

// ROOT includes for plotting
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>

// Helper function to ensure a directory exists
void ensure_directory_exists(const std::string &path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        mkdir(path.c_str(), 0777);
    }
}

// Structure to hold bin data
struct CrossSectionBinData {
    int bin_number;
    double xB_min, xB_max, xB_avg;
    double Q2_min, Q2_max, Q2_avg;
    double t_min, t_max, t_avg;
    double phi_min, phi_max, phi_avg;
    double cross_section;
    double stat_uncertainty;
    double sys_uncertainty;
};

// Type alias for a vector of CrossSectionBinData
typedef std::vector<CrossSectionBinData> CrossSectionData;

// Helper function to read data from the first CSV file format
CrossSectionData read_first_csv(const std::string &filename) {
    CrossSectionData data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return data;
    }

    std::string line;
    // Read header
    if (!std::getline(file, line)) {
        std::cerr << "Error: Cannot read header from file " << filename << std::endl;
        return data;
    }

    // Split header
    std::vector<std::string> headers;
    std::istringstream ss_header(line);
    std::string header;
    while (std::getline(ss_header, header, ',')) {
        headers.push_back(header);
    }

    // Map header names to indices
    std::map<std::string, int> header_indices;
    for (size_t i = 0; i < headers.size(); ++i) {
        header_indices[headers[i]] = i;
    }

    // Required columns
    std::string bin_col = "Bin";
    std::string xB_avg_col = "xB_avg";
    std::string xB_min_col = "xB_min";
    std::string xB_max_col = "xB_max";
    std::string Q2_avg_col = "Q2_avg";
    std::string Q2_min_col = "Q2_min";
    std::string Q2_max_col = "Q2_max";
    std::string t_avg_col = "t_avg";
    std::string t_min_col = "t_min";
    std::string t_max_col = "t_max";
    std::string phi_avg_col = "phi_avg";
    std::string phi_min_col = "phi_min";
    std::string phi_max_col = "phi_max";
    std::string cross_section_col = "cross sections, ep->epg, exp";
    std::string stat_unc_col = "cross sections, ep->epg, exp, stat. unc.";
    std::string sys_unc_col = "cross sections, ep->epg, exp, syst. unc. (up)";

    // Find indices of required columns
    std::vector<std::string> required_cols = {bin_col, xB_avg_col, xB_min_col, xB_max_col,
                                              Q2_avg_col, Q2_min_col, Q2_max_col,
                                              t_avg_col, t_min_col, t_max_col,
                                              phi_avg_col, phi_min_col, phi_max_col,
                                              cross_section_col};
    for (const auto &col : required_cols) {
        if (header_indices.find(col) == header_indices.end()) {
            std::cerr << "Error: Column " << col << " not found in " << filename << std::endl;
            return data;
        }
    }

    // Since 'stat_unc_col' is 4 columns after 'cross_section_col', find its index
    int cross_section_idx = header_indices[cross_section_col];
    int stat_unc_idx = cross_section_idx + 4;
    if (stat_unc_idx >= headers.size()) {
        std::cerr << "Error: Stat. unc. column index out of range in " << filename << std::endl;
        return data;
    }
    int sys_unc_idx = stat_unc_idx + 1;
    if (sys_unc_idx >= headers.size()) {
        std::cerr << "Error: Syst. unc. column index out of range in " << filename << std::endl;
        return data;
    }

    // Read data lines
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        if (tokens.size() <= sys_unc_idx) continue; // Skip incomplete lines

        CrossSectionBinData entry;
        entry.bin_number = std::stoi(tokens[header_indices[bin_col]]);
        entry.xB_avg = std::stod(tokens[header_indices[xB_avg_col]]);
        entry.xB_min = std::stod(tokens[header_indices[xB_min_col]]);
        entry.xB_max = std::stod(tokens[header_indices[xB_max_col]]);
        entry.Q2_avg = std::stod(tokens[header_indices[Q2_avg_col]]);
        entry.Q2_min = std::stod(tokens[header_indices[Q2_min_col]]);
        entry.Q2_max = std::stod(tokens[header_indices[Q2_max_col]]);
        entry.t_avg = std::stod(tokens[header_indices[t_avg_col]]);
        entry.t_min = std::stod(tokens[header_indices[t_min_col]]);
        entry.t_max = std::stod(tokens[header_indices[t_max_col]]);
        entry.phi_avg = std::stod(tokens[header_indices[phi_avg_col]]);
        entry.phi_min = std::stod(tokens[header_indices[phi_min_col]]);
        entry.phi_max = std::stod(tokens[header_indices[phi_max_col]]);
        entry.cross_section = std::stod(tokens[cross_section_idx]);
        entry.stat_uncertainty = std::stod(tokens[stat_unc_idx]);
        entry.sys_uncertainty = std::stod(tokens[sys_unc_idx]);

        data.push_back(entry);
    }

    file.close();
    return data;
}

// Function to read the second CSV file
CrossSectionData read_second_csv(const std::string &filename) {
    CrossSectionData data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return data;
    }

    std::string line;
    // Read header
    if (!std::getline(file, line)) {
        std::cerr << "Error: Cannot read header from file " << filename << std::endl;
        return data;
    }

    // Split header
    std::vector<std::string> headers;
    std::istringstream ss_header(line);
    std::string header;
    while (std::getline(ss_header, header, ',')) {
        headers.push_back(header);
    }

    // Map header names to indices
    std::map<std::string, int> header_indices;
    for (size_t i = 0; i < headers.size(); ++i) {
        header_indices[headers[i]] = i;
    }

    // Required columns
    std::string bin_col = "Bin";
    std::string xB_avg_col = "xB_avg";
    std::string xB_min_col = "xB_min";
    std::string xB_max_col = "xB_max";
    std::string Q2_avg_col = "Q2_avg";
    std::string Q2_min_col = "Q2_min";
    std::string Q2_max_col = "Q2_max";
    std::string t_avg_col = "t_avg";
    std::string t_min_col = "t_min";
    std::string t_max_col = "t_max";
    std::string phi_avg_col = "phi_avg";
    std::string phi_min_col = "phi_min";
    std::string phi_max_col = "phi_max";
    std::string cross_section_col = "fall_cross_section";
    std::string stat_unc_col = "fall_cross_section_stat_uncertainty";
    std::string sys_unc_col = "fall_cross_section_sys_uncertainty";

    // Check if required columns are present
    std::vector<std::string> required_cols = {bin_col, xB_avg_col, xB_min_col, xB_max_col,
                                              Q2_avg_col, Q2_min_col, Q2_max_col,
                                              t_avg_col, t_min_col, t_max_col,
                                              phi_avg_col, phi_min_col, phi_max_col,
                                              cross_section_col, stat_unc_col, sys_unc_col};
    for (const auto &col : required_cols) {
        if (header_indices.find(col) == header_indices.end()) {
            std::cerr << "Error: Column " << col << " not found in " << filename << std::endl;
            return data;
        }
    }

    // Read data lines
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        CrossSectionBinData entry;
        entry.bin_number = std::stoi(tokens[header_indices[bin_col]]);
        entry.xB_avg = std::stod(tokens[header_indices[xB_avg_col]]);
        entry.xB_min = std::stod(tokens[header_indices[xB_min_col]]);
        entry.xB_max = std::stod(tokens[header_indices[xB_max_col]]);
        entry.Q2_avg = std::stod(tokens[header_indices[Q2_avg_col]]);
        entry.Q2_min = std::stod(tokens[header_indices[Q2_min_col]]);
        entry.Q2_max = std::stod(tokens[header_indices[Q2_max_col]]);
        entry.t_avg = std::stod(tokens[header_indices[t_avg_col]]);
        entry.t_min = std::stod(tokens[header_indices[t_min_col]]);
        entry.t_max = std::stod(tokens[header_indices[t_max_col]]);
        entry.phi_avg = std::stod(tokens[header_indices[phi_avg_col]]);
        entry.phi_min = std::stod(tokens[header_indices[phi_min_col]]);
        entry.phi_max = std::stod(tokens[header_indices[phi_max_col]]);
        entry.cross_section = std::stod(tokens[header_indices[cross_section_col]]);
        entry.stat_uncertainty = std::stod(tokens[header_indices[stat_unc_col]]);
        entry.sys_uncertainty = std::stod(tokens[header_indices[sys_unc_col]]);

        data.push_back(entry);
    }

    file.close();
    return data;
}

// Function to find unique xB bins in the data
std::vector<std::pair<double, double>> find_unique_xB_bins(const CrossSectionData &data) {
    std::set<std::pair<double, double>> unique_bins;
    for (const auto &bin : data) {
        unique_bins.emplace(bin.xB_min, bin.xB_max);
    }
    return std::vector<std::pair<double, double>>(unique_bins.begin(), unique_bins.end());
}

// Function to filter data for a specific xB bin
CrossSectionData filter_data_by_xB(const CrossSectionData &data, const std::pair<double, double> &xB_range) {
    CrossSectionData filtered_data;
    for (const auto &bin : data) {
        if (bin.xB_min == xB_range.first && bin.xB_max == xB_range.second) {
            filtered_data.push_back(bin);
        }
    }
    return filtered_data;
}

// Plotting function for each xB bin with comparison between datasets
void plot_for_xB_bin(const CrossSectionData &data_first, const CrossSectionData &data_second, int xB_index, const std::string &output_dir) {
    // Step 1: Identify unique (Q2min, Q2max, tmin, tmax) bins in the first dataset
    std::map<std::tuple<double, double, double, double>, CrossSectionData> qt_bins_first, qt_bins_second;

    // Populate qt_bins_first and qt_bins_second
    for (const auto &bin : data_first) {
        auto key = std::make_tuple(bin.Q2_min, bin.Q2_max, bin.t_min, bin.t_max);
        qt_bins_first[key].push_back(bin);
    }
    for (const auto &bin : data_second) {
        auto key = std::make_tuple(bin.Q2_min, bin.Q2_max, bin.t_min, bin.t_max);
        qt_bins_second[key].push_back(bin);
    }

    // Step 2: Determine the grid size for subplots based on the first dataset's bins
    int num_plots = qt_bins_first.size();
    int grid_size = std::ceil(std::sqrt(num_plots));

    // Adjust grid size for canvas layout (if needed)
    if ((xB_index == 3 || xB_index == 4) && grid_size * (grid_size - 1) >= num_plots) {
        grid_size -= 1;
    }

    // Step 3: Set up the canvas
    TCanvas canvas("canvas", "Cross Section Cross Check", 1200, 1200);
    canvas.Divide(grid_size, grid_size, 0.02, 0.02); // Small padding between pads

    int plot_index = 1;
    for (const auto &[qt_key, bins_first] : qt_bins_first) {
        canvas.cd(plot_index);
        gPad->SetLeftMargin(0.15);  // Adds padding to the left
        gPad->SetBottomMargin(0.15); // Adds padding to the bottom
        gPad->SetLogy();  // Set log scale for the y-axis

        // Prepare vectors for the first dataset
        std::vector<double> phi_values_first, cs_first, phi_errors_first, stat_err_first, sys_err_first;
        for (const auto &bin : bins_first) {
            phi_values_first.push_back(bin.phi_avg);
            cs_first.push_back(bin.cross_section);
            stat_err_first.push_back(bin.stat_uncertainty);
            sys_err_first.push_back(bin.sys_uncertainty);
            double phi_bin_width = (bin.phi_max - bin.phi_min) / 2.0; // Calculate half-width of the phi bin
            phi_errors_first.push_back(phi_bin_width);
        }

        // Create TGraphErrors for statistical uncertainties (first dataset)
        TGraphErrors *graph_first_stat = new TGraphErrors(phi_values_first.size(),
                                                          &phi_values_first[0], &cs_first[0],
                                                          &phi_errors_first[0], &stat_err_first[0]);
        graph_first_stat->SetMarkerStyle(20);
        graph_first_stat->SetMarkerColor(kBlue);
        graph_first_stat->SetLineColor(kBlue);

        // Create TGraphAsymmErrors for systematic uncertainties (first dataset)
        std::vector<double> sys_err_first_down = sys_err_first;
        std::vector<double> sys_err_first_up = sys_err_first;

        TGraphAsymmErrors *graph_first_sys = new TGraphAsymmErrors(phi_values_first.size(),
                                                                   &phi_values_first[0], &cs_first[0],
                                                                   &phi_errors_first[0], &phi_errors_first[0],
                                                                   &sys_err_first_down[0], &sys_err_first_up[0]);
        graph_first_sys->SetFillColorAlpha(kBlue, 0.35);
        graph_first_sys->SetLineColor(kBlue);

        // Prepare and plot the second dataset if there’s a matching key
        auto it_second = qt_bins_second.find(qt_key);
        if (it_second != qt_bins_second.end()) {
            const auto &bins_second = it_second->second;

            std::vector<double> phi_values_second, cs_second, phi_errors_second, stat_err_second, sys_err_second;

            for (const auto &bin : bins_second) {
                phi_values_second.push_back(bin.phi_avg);
                cs_second.push_back(bin.cross_section);
                stat_err_second.push_back(bin.stat_uncertainty);
                sys_err_second.push_back(bin.sys_uncertainty);
                double phi_bin_width = (bin.phi_max - bin.phi_min) / 2.0;
                phi_errors_second.push_back(phi_bin_width);
            }

            // Create TGraphErrors for statistical uncertainties (second dataset)
            TGraphErrors *graph_second_stat = new TGraphErrors(phi_values_second.size(),
                                                               &phi_values_second[0], &cs_second[0],
                                                               &phi_errors_second[0], &stat_err_second[0]);
            graph_second_stat->SetMarkerStyle(21);
            graph_second_stat->SetMarkerColor(kRed);
            graph_second_stat->SetLineColor(kRed);

            // Create TGraphAsymmErrors for systematic uncertainties (second dataset)
            std::vector<double> sys_err_second_down = sys_err_second;
            std::vector<double> sys_err_second_up = sys_err_second;

            TGraphAsymmErrors *graph_second_sys = new TGraphAsymmErrors(phi_values_second.size(),
                                                                        &phi_values_second[0], &cs_second[0],
                                                                        &phi_errors_second[0], &phi_errors_second[0],
                                                                        &sys_err_second_down[0], &sys_err_second_up[0]);
            graph_second_sys->SetFillColorAlpha(kRed, 0.35);
            graph_second_sys->SetLineColor(kRed);

            // Determine y-axis range
            double y_min = 0.1;
            double y_max = 10;

            // Create a frame for the plot
            TH1F *frame = gPad->DrawFrame(0, y_min, 360, y_max);
            frame->GetXaxis()->SetTitle("#phi [deg]");
            frame->GetYaxis()->SetTitle("d#sigma/dx_{B}dQ^{2}d|t|d#phi (nb/GeV^{4})");

            // Draw systematic uncertainties as bands
            graph_first_sys->Draw("E3 SAME");
            graph_second_sys->Draw("E3 SAME");

            // Draw statistical uncertainties as error bars
            graph_first_stat->Draw("P SAME");
            graph_second_stat->Draw("P SAME");

            // Create and draw legend
            TLegend* legend = new TLegend(0.55, 0.65, 0.9, 0.85);
            legend->SetTextSize(0.03);
            legend->AddEntry(graph_first_stat, "First CSV", "lep");
            legend->AddEntry(graph_second_stat, "Second CSV", "lep");
            legend->Draw();

            // Add title with bin information
            double xB_avg = bins_first[0].xB_avg;
            double Q2_avg = bins_first[0].Q2_avg;
            double t_avg = bins_first[0].t_avg;
            std::string title = Form("x_{B} = %.2f, Q^{2} = %.2f, -t = %.2f", xB_avg, Q2_avg, t_avg);
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5

        } else {
            // If no matching (Q2, t) bin in the second dataset, plot only the first dataset
            double y_min = 0.1;
            double y_max = 10;

            // Create a frame for the plot
            TH1F *frame = gPad->DrawFrame(0, y_min, 360, y_max);
            frame->GetXaxis()->SetTitle("#phi [deg]");
            frame->GetYaxis()->SetTitle("d#sigma/dx_{B}dQ^{2}d|t|d#phi (nb/GeV^{4})");

            // Draw systematic uncertainties as bands
            graph_first_sys->Draw("E3 SAME");

            // Draw statistical uncertainties as error bars
            graph_first_stat->Draw("P SAME");

            // Create and draw legend
            TLegend* legend = new TLegend(0.55, 0.65, 0.9, 0.85);
            legend->SetTextSize(0.03);
            legend->AddEntry(graph_first_stat, "First CSV", "lep");
            legend->Draw();

            // Add title with bin information
            double xB_avg = bins_first[0].xB_avg;
            double Q2_avg = bins_first[0].Q2_avg;
            double t_avg = bins_first[0].t_avg;
            std::string title = Form("x_{B} = %.2f, Q^{2} = %.2f, -t = %.2f", xB_avg, Q2_avg, t_avg);
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5
        }

        plot_index++;
    }

    // Save the canvas
    std::string save_path = Form("output/cross_section_cross_check/cross_section_cross_check_xB_%d.pdf", xB_index);
    canvas.SaveAs(save_path.c_str());
}

// Main function to control the import, processing, and plotting
void plot_cross_section_comparison(const std::string &first_csv_file,
                                   const std::string &second_csv_file) {
    ensure_directory_exists("output/cross_section_cross_check/");

    // Read data from CSV files
    CrossSectionData data_first_csv = read_first_csv(first_csv_file);
    CrossSectionData data_second_csv = read_second_csv(second_csv_file);

    // Find unique xB bins
    auto unique_xB_bins = find_unique_xB_bins(data_first_csv);

    int xB_index = 0;
    for (const auto &xB_range : unique_xB_bins) {
        // Filter data for the current xB bin range in both datasets
        auto filtered_data_first = filter_data_by_xB(data_first_csv, xB_range);
        auto filtered_data_second = filter_data_by_xB(data_second_csv, xB_range);

        // Pass both datasets to the plotting function
        plot_for_xB_bin(filtered_data_first, filtered_data_second, xB_index, output_dir);
        xB_index++;
    }
}