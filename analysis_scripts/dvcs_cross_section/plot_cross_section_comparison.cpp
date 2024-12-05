// plot_cross_section_comparison.cpp

#include "plot_cross_section_comparison.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>
#include <set>
#include <tuple>
#include <algorithm>

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TColor.h>

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

typedef std::vector<CrossSectionBinData> CrossSectionData;

CrossSectionData read_first_csv(const std::string& filename);
CrossSectionData read_second_csv(const std::string& filename);

void plot_cross_section_comparison(const std::string& first_csv_file,
                                   const std::string& second_csv_file,
                                   const std::string& output_dir) {
    // Read data from CSV files
    CrossSectionData data_first_csv = read_first_csv(first_csv_file);
    CrossSectionData data_second_csv = read_second_csv(second_csv_file);

    // Organize data by bins
    // Map from (bin_number) to vector of CrossSectionBinData for each phi bin
    typedef int BinKey; // bin_number
    std::map<BinKey, std::vector<CrossSectionBinData>> bins_first_csv;
    std::map<BinKey, std::vector<CrossSectionBinData>> bins_second_csv;

    // Organize first CSV data
    for (const auto& entry : data_first_csv) {
        BinKey key = entry.bin_number;
        bins_first_csv[key].push_back(entry);
    }

    // Organize second CSV data
    for (const auto& entry : data_second_csv) {
        BinKey key = entry.bin_number;
        bins_second_csv[key].push_back(entry);
    }

    // Now, for each bin, we need to plot the data
    // We need to determine the number of canvases and subplots

    // Collect all unique bins
    std::set<BinKey> all_bins;
    for (const auto& item : bins_first_csv) {
        all_bins.insert(item.first);
    }
    for (const auto& item : bins_second_csv) {
        all_bins.insert(item.first);
    }

    // Determine number of bins
    int n_bins = all_bins.size();

    // Determine number of columns and rows for the canvas
    int n_columns = std::ceil(std::sqrt(n_bins));
    int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    TCanvas* canvas = new TCanvas("c1", "Cross Section Comparison", canvas_width, canvas_height);
    canvas->Divide(n_columns, n_rows);

    gStyle->SetOptStat(0);

    int pad_idx = 1;

    for (const auto& bin_key : all_bins) {
        canvas->cd(pad_idx);
        TPad* pad = (TPad*)gPad;
        pad->SetLeftMargin(0.15);
        pad->SetBottomMargin(0.15);
        pad->SetLogy();

        // Get data for this bin
        std::vector<CrossSectionBinData> data_first = bins_first_csv[bin_key];
        std::vector<CrossSectionBinData> data_second = bins_second_csv[bin_key];

        // Prepare data for plotting
        // For each dataset, extract phi, cross_section, stat_uncertainty, sys_uncertainty

        // Sort data by phi_avg
        auto sort_by_phi = [](const CrossSectionBinData& a, const CrossSectionBinData& b) {
            return a.phi_avg < b.phi_avg;
        };
        std::sort(data_first.begin(), data_first.end(), sort_by_phi);
        std::sort(data_second.begin(), data_second.end(), sort_by_phi);

        // Prepare arrays for plotting
        int n_points_first = data_first.size();
        int n_points_second = data_second.size();

        std::vector<double> phi_first(n_points_first);
        std::vector<double> cs_first(n_points_first);
        std::vector<double> stat_err_first(n_points_first);
        std::vector<double> sys_err_first(n_points_first);

        std::vector<double> phi_second(n_points_second);
        std::vector<double> cs_second(n_points_second);
        std::vector<double> stat_err_second(n_points_second);
        std::vector<double> sys_err_second(n_points_second);

        for (int i = 0; i < n_points_first; ++i) {
            phi_first[i] = data_first[i].phi_avg;
            cs_first[i] = data_first[i].cross_section;
            stat_err_first[i] = data_first[i].stat_uncertainty;
            sys_err_first[i] = data_first[i].sys_uncertainty;
        }

        for (int i = 0; i < n_points_second; ++i) {
            phi_second[i] = data_second[i].phi_avg;
            cs_second[i] = data_second[i].cross_section;
            stat_err_second[i] = data_second[i].stat_uncertainty;
            sys_err_second[i] = data_second[i].sys_uncertainty;
        }

        // Create TGraphErrors for statistical uncertainties
        TGraphErrors* graph_first_stat = new TGraphErrors(n_points_first, &phi_first[0], &cs_first[0], nullptr, &stat_err_first[0]);
        TGraphErrors* graph_second_stat = new TGraphErrors(n_points_second, &phi_second[0], &cs_second[0], nullptr, &stat_err_second[0]);

        // Create TGraphAsymmErrors for systematic uncertainties
        // Since the systematic uncertainties are symmetric, we can use the same value for both up and down
        std::vector<double> sys_err_first_up(n_points_first);
        std::vector<double> sys_err_first_down(n_points_first);
        std::vector<double> sys_err_second_up(n_points_second);
        std::vector<double> sys_err_second_down(n_points_second);

        for (int i = 0; i < n_points_first; ++i) {
            sys_err_first_up[i] = sys_err_first[i];
            sys_err_first_down[i] = sys_err_first[i];
        }
        for (int i = 0; i < n_points_second; ++i) {
            sys_err_second_up[i] = sys_err_second[i];
            sys_err_second_down[i] = sys_err_second[i];
        }

        TGraphAsymmErrors* graph_first_sys = new TGraphAsymmErrors(n_points_first, &phi_first[0], &cs_first[0], nullptr, nullptr, &sys_err_first_down[0], &sys_err_first_up[0]);
        TGraphAsymmErrors* graph_second_sys = new TGraphAsymmErrors(n_points_second, &phi_second[0], &cs_second[0], nullptr, nullptr, &sys_err_second_down[0], &sys_err_second_up[0]);

        // Set styles
        graph_first_stat->SetMarkerColor(kBlue);
        graph_first_stat->SetMarkerStyle(20);
        graph_first_stat->SetLineColor(kBlue);

        graph_second_stat->SetMarkerColor(kRed);
        graph_second_stat->SetMarkerStyle(21);
        graph_second_stat->SetLineColor(kRed);

        graph_first_sys->SetFillColorAlpha(kBlue, 0.35);
        graph_first_sys->SetLineColor(kBlue);

        graph_second_sys->SetFillColorAlpha(kRed, 0.35);
        graph_second_sys->SetLineColor(kRed);

        // Determine y-axis range
        double y_min = 0.1;
        double y_max = 10;

        TH1F* frame = pad->DrawFrame(0, y_min, 360, y_max);
        frame->GetXaxis()->SetTitle("#phi [deg]");
        frame->GetYaxis()->SetTitle("d#sigma/dx_{B}dQ^{2}d|t|d#phi (nb/GeV^{4})");

        // Draw systematic uncertainties as bands
        graph_first_sys->Draw("E3 SAME");
        graph_second_sys->Draw("E3 SAME");

        // Draw statistical uncertainties as error bars
        graph_first_stat->Draw("P SAME");
        graph_second_stat->Draw("P SAME");

        // Add legend
        TLegend* legend = new TLegend(0.15, 0.65, 0.5, 0.85);
        legend->AddEntry(graph_first_stat, "ep->epg (First CSV)", "lep");
        legend->AddEntry(graph_second_stat, "ep->epg (Second CSV)", "lep");
        legend->SetTextSize(0.04);
        legend->Draw();

        // Add title with bin information
        int bin_number = bin_key;

        // Get average values from first dataset if available, else from second
        double xB_avg = (data_first.size() > 0) ? data_first[0].xB_avg : data_second[0].xB_avg;
        double Q2_avg = (data_first.size() > 0) ? data_first[0].Q2_avg : data_second[0].Q2_avg;
        double t_avg = (data_first.size() > 0) ? data_first[0].t_avg : data_second[0].t_avg;

        std::string title = Form("Bin %d: x_{B}=%.3f, Q^{2}=%.3f, |t|=%.3f", bin_number, xB_avg, Q2_avg, t_avg);
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextAlign(22); // Center alignment
        latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5

        ++pad_idx;
    }

    // Save canvas
    std::string filename = output_dir + "/cross_section_cross_check.pdf";
    canvas->SaveAs(filename.c_str());

    // Clean up
    delete canvas;
}

// Function to read the first CSV file
CrossSectionData read_first_csv(const std::string& filename) {
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
    std::string Q2_avg_col = "Q2_avg";
    std::string t_avg_col = "t_avg";
    std::string phi_avg_col = "phi_avg";
    std::string cross_section_col = "cross sections, ep->epg, exp";
    std::string stat_unc_col = "cross sections, ep->epg, exp, stat. unc.";
    std::string sys_unc_col = "cross sections, ep->epg, exp, syst. unc. (up)";

    // Find indices of required columns
    std::vector<std::string> required_cols = {bin_col, xB_avg_col, Q2_avg_col, t_avg_col, phi_avg_col, cross_section_col};
    for (const auto& col : required_cols) {
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
        entry.Q2_avg = std::stod(tokens[header_indices[Q2_avg_col]]);
        entry.t_avg = std::stod(tokens[header_indices[t_avg_col]]);
        entry.phi_avg = std::stod(tokens[header_indices[phi_avg_col]]);
        entry.cross_section = std::stod(tokens[cross_section_idx]);
        entry.stat_uncertainty = std::stod(tokens[stat_unc_idx]);
        entry.sys_uncertainty = std::stod(tokens[sys_unc_idx]);

        data.push_back(entry);
    }

    file.close();
    return data;
}

// Function to read the second CSV file
CrossSectionData read_second_csv(const std::string& filename) {
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
    std::string Q2_avg_col = "Q2_avg";
    std::string t_avg_col = "t_avg";
    std::string phi_avg_col = "phi_avg";
    std::string cross_section_col = "fall_cross_section";
    std::string stat_unc_col = "fall_cross_section_stat_uncertainty";
    std::string sys_unc_col = "fall_cross_section_sys_uncertainty";

    // Check if required columns are present
    std::vector<std::string> required_cols = {bin_col, xB_avg_col, Q2_avg_col, t_avg_col, phi_avg_col, cross_section_col, stat_unc_col, sys_unc_col};
    for (const auto& col : required_cols) {
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
        entry.Q2_avg = std::stod(tokens[header_indices[Q2_avg_col]]);
        entry.t_avg = std::stod(tokens[header_indices[t_avg_col]]);
        entry.phi_avg = std::stod(tokens[header_indices[phi_avg_col]]);
        entry.cross_section = std::stod(tokens[header_indices[cross_section_col]]);
        entry.stat_uncertainty = std::stod(tokens[header_indices[stat_unc_col]]);
        entry.sys_uncertainty = std::stod(tokens[header_indices[sys_unc_col]]);

        data.push_back(entry);
    }

    file.close();
    return data;
}