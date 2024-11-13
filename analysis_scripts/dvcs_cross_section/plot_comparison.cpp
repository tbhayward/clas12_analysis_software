#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <filesystem>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TSystem.h>
#include "bin_boundaries.h"

// Function to read unfolded yield data from CSV
std::vector<std::vector<double>> read_unfolded_yields(const std::string& filename, const std::vector<std::string>& columns) {
    std::vector<std::vector<double>> yield_data(columns.size());
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return yield_data;
    }

    // Read header and map columns
    std::string line;
    std::getline(file, line); // Header line
    std::istringstream header_stream(line);
    std::map<std::string, int> column_indices;
    std::string header_cell;
    int col_idx = 0;
    
    while (std::getline(header_stream, header_cell, ',')) {
        column_indices[header_cell] = col_idx++;
    }
    
    // Validate that required columns are in the CSV
    for (const auto& col : columns) {
        if (column_indices.find(col) == column_indices.end()) {
            std::cerr << "Error: Column " << col << " not found in CSV file." << std::endl;
            return yield_data;
        }
    }
    
    // Read each row and extract relevant columns
    while (std::getline(file, line)) {
        std::istringstream line_stream(line);
        std::string cell;
        int current_idx = 0;
        std::vector<std::string> cells(col_idx); // Total columns from header

        while (std::getline(line_stream, cell, ',')) {
            cells[current_idx++] = cell;
        }

        // Populate yield data for required columns
        for (size_t i = 0; i < columns.size(); ++i) {
            yield_data[i].push_back(std::stod(cells[column_indices[columns[i]]]));
        }
    }

    return yield_data;
}

void plot_comparison(const std::string& output_dir, const std::vector<BinBoundary>& bin_boundaries,
                     const std::string& csv_fa18_inb, const std::string& csv_fa18_out) {
    // Define the columns of interest for each CSV file
    std::vector<std::string> columns_fa18_inb = {
        "ep->e'pgamma unfolded_yield_Fa18Inb",
        "ep->e'ppi0 unfolded_yield_Fa18Inb"
    };
    std::vector<std::string> columns_fa18_out = {
        "ep->e'pgamma unfolded_yield_Fa18Out",
        "ep->e'ppi0 unfolded_yield_Fa18Out"
    };

    // Read unfolded yield data from both CSV files
    std::vector<std::vector<double>> yields_fa18_inb = read_unfolded_yields(csv_fa18_inb, columns_fa18_inb);
    std::vector<std::vector<double>> yields_fa18_out = read_unfolded_yields(csv_fa18_out, columns_fa18_out);

    if (yields_fa18_inb[0].empty() || yields_fa18_out[0].empty()) {
        std::cerr << "Error: Unfolded yield data could not be read." << std::endl;
        return;
    }

    // Start plotting setup
    gStyle->SetOptStat(0);
    int canvas_width = 1200;
    int canvas_height = 800;
    TCanvas* canvas = new TCanvas("c_comparison", "Comparison Plots", canvas_width, canvas_height);
    canvas->Divide(2, 1);

    // DVCS Comparison Plot (Inb vs Out)
    canvas->cd(1);
    TGraphErrors* graph_dvcs_inb = new TGraphErrors(bin_boundaries.size());
    TGraphErrors* graph_dvcs_out = new TGraphErrors(bin_boundaries.size());

    for (size_t i = 0; i < bin_boundaries.size(); ++i) {
        double xB_avg = bin_boundaries[i].xB_avg;
        graph_dvcs_inb->SetPoint(i, xB_avg, yields_fa18_inb[0][i]);
        graph_dvcs_out->SetPoint(i, xB_avg, yields_fa18_out[0][i]);
    }

    graph_dvcs_inb->SetMarkerColor(kBlue);
    graph_dvcs_inb->SetMarkerStyle(20);
    graph_dvcs_inb->SetLineColor(kBlue);
    graph_dvcs_inb->SetTitle("DVCS Yield Comparison");
    graph_dvcs_inb->GetXaxis()->SetTitle("x_{B}");
    graph_dvcs_inb->GetYaxis()->SetTitle("Unfolded Yield");
    graph_dvcs_inb->Draw("AP");

    graph_dvcs_out->SetMarkerColor(kRed);
    graph_dvcs_out->SetMarkerStyle(24);
    graph_dvcs_out->SetLineColor(kRed);
    graph_dvcs_out->Draw("P SAME");

    TLegend* legend_dvcs = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_dvcs->AddEntry(graph_dvcs_inb, "Fa18 Inbending", "p");
    legend_dvcs->AddEntry(graph_dvcs_out, "Fa18 Outbending", "p");
    legend_dvcs->Draw();

    // Eppi0 Comparison Plot (Inb vs Out)
    canvas->cd(2);
    TGraphErrors* graph_eppi0_inb = new TGraphErrors(bin_boundaries.size());
    TGraphErrors* graph_eppi0_out = new TGraphErrors(bin_boundaries.size());

    for (size_t i = 0; i < bin_boundaries.size(); ++i) {
        double xB_avg = bin_boundaries[i].xB_avg;
        graph_eppi0_inb->SetPoint(i, xB_avg, yields_fa18_inb[1][i]);
        graph_eppi0_out->SetPoint(i, xB_avg, yields_fa18_out[1][i]);
    }

    graph_eppi0_inb->SetMarkerColor(kBlue);
    graph_eppi0_inb->SetMarkerStyle(20);
    graph_eppi0_inb->SetLineColor(kBlue);
    graph_eppi0_inb->SetTitle("Eppi0 Yield Comparison");
    graph_eppi0_inb->GetXaxis()->SetTitle("x_{B}");
    graph_eppi0_inb->GetYaxis()->SetTitle("Unfolded Yield");
    graph_eppi0_inb->Draw("AP");

    graph_eppi0_out->SetMarkerColor(kRed);
    graph_eppi0_out->SetMarkerStyle(24);
    graph_eppi0_out->SetLineColor(kRed);
    graph_eppi0_out->Draw("P SAME");

    TLegend* legend_eppi0 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_eppi0->AddEntry(graph_eppi0_inb, "Fa18 Inbending", "p");
    legend_eppi0->AddEntry(graph_eppi0_out, "Fa18 Outbending", "p");
    legend_eppi0->Draw();

    // Save canvas to file
    std::string filename = output_dir + "/comparison_yields.pdf";
    canvas->SaveAs(filename.c_str());

    // Clean up
    delete graph_dvcs_inb;
    delete graph_dvcs_out;
    delete graph_eppi0_inb;
    delete graph_eppi0_out;
    delete canvas;
}