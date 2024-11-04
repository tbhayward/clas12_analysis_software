// rad_study.cpp
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <sys/stat.h> // For mkdir

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cout << "Usage: rad_study data_file.root input_file.root ISR_or_FSR_string" << std::endl;
        return 1;
    }

    std::string data_filename = argv[1];
    std::string input_filename = argv[2];
    std::string rad_type = argv[3];

    // Open the ROOT files
    TFile* data_file = TFile::Open(data_filename.c_str());
    TFile* input_file = TFile::Open(input_filename.c_str());

    if (!data_file || !input_file) {
        std::cout << "Error opening input files." << std::endl;
        return 1;
    }

    // Get the TTrees
    TTree* data_tree = (TTree*)data_file->Get("PhysicsEvents");
    TTree* input_tree = (TTree*)input_file->Get("PhysicsEvents");

    if (!data_tree || !input_tree) {
        std::cout << "Error getting TTrees from files." << std::endl;
        return 1;
    }

    // Create output directory
    mkdir("output", 0755);
    mkdir("output/rad_study", 0755);

    // Set style options
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.05);

    // Variables to plot
    const int nVars = 11;
    std::string varNames[nVars] = {"e_p", "e_theta", "Mx2", "phi", "pT", "Q2", "W", "x", "xF", "xi", "z"};
    double varMin[nVars] = {1, 0, -1, 0, 0, 0, 2, 0, -1, 0, 0};
    double varMax[nVars] = {9, 35, 4, 2*TMath::Pi(), 1.2, 10, 4, 0.7, 1, 1, 1};
    std::string xLabels[nVars] = {"e_{p} (GeV)", "e_{#theta}", "M_{x}^{2} (GeV^{2})", "#phi", "P_{T} (GeV)", "Q^{2} (GeV^{2})", "W (GeV)", "x_{B}", "x_{F}", "#xi", "z"};
    int nBins[nVars] = {80, 70, 100, 64, 60, 100, 100, 70, 100, 100, 100};

    // Selection criteria
    std::string selection = "Mx2 > 1.8225";

    for (int i = 0; i < nVars; ++i) {
        // Create histograms
        std::string hist_name_data = "h_data_" + varNames[i];
        std::string hist_name_input = "h_input_" + varNames[i];

        TH1D* h_data = new TH1D(hist_name_data.c_str(), "", nBins[i], varMin[i], varMax[i]);
        TH1D* h_input = new TH1D(hist_name_input.c_str(), "", nBins[i], varMin[i], varMax[i]);

        // Apply selection except for Mx2
        std::string data_selection = (varNames[i] == "Mx2") ? "" : selection;
        std::string input_selection = data_selection;

        // Project data onto histograms
        data_tree->Project(hist_name_data.c_str(), varNames[i].c_str(), data_selection.c_str());
        input_tree->Project(hist_name_input.c_str(), varNames[i].c_str(), input_selection.c_str());

        // Create ratio histogram
        TH1D* h_ratio = (TH1D*)h_input->Clone(("h_ratio_" + varNames[i]).c_str());
        h_ratio->Divide(h_data);

        // Create canvas and draw
        TCanvas* c = new TCanvas("c", "c", 800, 600);
        c->SetGrid();

        // Adjust margins to prevent cutting off labels
        c->SetLeftMargin(0.15);
        c->SetBottomMargin(0.15);

        h_ratio->SetTitle("");
        h_ratio->GetXaxis()->SetTitle(xLabels[i].c_str());
        h_ratio->GetYaxis()->SetTitle(("ratio (" + rad_type + "/data)").c_str());
        h_ratio->GetXaxis()->SetLabelSize(0.04);
        h_ratio->GetYaxis()->SetLabelSize(0.04);
        h_ratio->GetXaxis()->SetTitleSize(0.05);
        h_ratio->GetYaxis()->SetTitleSize(0.05);
        h_ratio->GetXaxis()->SetTitleOffset(1.2);
        h_ratio->GetYaxis()->SetTitleOffset(1.4);
        h_ratio->SetLineColor(kBlue);
        h_ratio->SetLineWidth(2);

        // Set y-axis range from 0.7 to 1.3
        h_ratio->GetYaxis()->SetRangeUser(0.4, 1.3);

        h_ratio->Draw("HIST");

        // Save the plot
        std::string output_filename = "output/rad_study/" + varNames[i] + "_" + rad_type + ".png";
        c->SaveAs(output_filename.c_str());

        // Clean up
        delete h_data;
        delete h_input;
        delete h_ratio;
        delete c;
    }

    // Close files
    data_file->Close();
    input_file->Close();

    return 0;
}