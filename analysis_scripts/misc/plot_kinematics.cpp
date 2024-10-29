#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h> // For mkdir

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TAxis.h"

int main(int argc, char** argv) {

    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " input_root_file" << std::endl;
        return 1;
    }

    std::string inputFileName = argv[1];

    TFile *file = TFile::Open(inputFileName.c_str(), "READ");
    if (!file || !file->IsOpen()) {
        std::cerr << "Error opening file " << inputFileName << std::endl;
        return 1;
    }

    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: TTree PhysicsEvents not found in file " << inputFileName << std::endl;
        return 1;
    }

    // Define cuts
    std::string cuts = "open_angle_ep2 > 5 && theta_gamma_gamma < 0.7 && Emiss2 < 1 && pTmiss < 0.15";

    // Define bin edges
    std::vector<double> xB_bin_edges = {0.062, 0.090, 0.118, 0.155, 0.204, 0.268, 0.357, 0.446, 0.581};
    std::vector<double> Q2_bin_edges = {1.000, 1.200, 1.456, 1.912, 2.510, 3.295, 4.326, 5.761};
    std::vector<double> t_bin_edges = {0.110, 0.150, 0.250, 0.400, 0.600, 0.800, 1.000};

    // Create histograms
    TH2D *hQ2_vs_x = new TH2D("hQ2_vs_x", "Q^{2} vs x_{B};x_{B};Q^{2} (GeV^{2})", 100, 0.05, 0.6, 100, 0.6, 8.0);
    TH1D *h_t = new TH1D("h_t", "-t distribution;-t (GeV^{2});Counts", 100, 0, 1.0);
    TH1D *h_phi = new TH1D("h_phi", "#phi distribution;#phi (degrees);Counts", 360, 0, 360);

    // Fill histograms
    tree->Draw("Q2:x>>hQ2_vs_x", cuts.c_str(), "COLZ");
    tree->Draw("(-(t1)>>h_t", cuts.c_str());
    tree->Draw("(180/3.14159)*phi>>h_phi", cuts.c_str());

    // Set style
    gStyle->SetOptStat(0);

    // Create output directory if it doesn't exist
    struct stat info;
    if (stat("output", &info) != 0) {
        mkdir("output", 0777);
    }

    // Draw hQ2_vs_x
    TCanvas *c1 = new TCanvas("c1", "Q^{2} vs x_{B}", 800, 600);
    hQ2_vs_x->Draw("COLZ");

    // Draw bin lines for xB
    for (size_t i = 0; i < xB_bin_edges.size(); ++i) {
        double x = xB_bin_edges[i];
        TLine *line = new TLine(x, 0.6, x, 8.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("same");
    }

    // Draw bin lines for Q2
    for (size_t i = 0; i < Q2_bin_edges.size(); ++i) {
        double y = Q2_bin_edges[i];
        TLine *line = new TLine(0.05, y, 0.6, y);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("same");
    }

    // Save the plot
    c1->SaveAs("output/Q2_vs_x.pdf");

    // Draw h_t
    TCanvas *c2 = new TCanvas("c2", "-t distribution", 800, 600);
    h_t->Draw();

    // Draw vertical lines at t bin edges
    for (size_t i = 0; i < t_bin_edges.size(); ++i) {
        double x = t_bin_edges[i];
        TLine *line = new TLine(x, 0, x, h_t->GetMaximum());
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("same");
    }

    // Save the plot
    c2->SaveAs("output/t_distribution.pdf");

    // Draw h_phi
    TCanvas *c3 = new TCanvas("c3", "#phi distribution", 800, 600);
    h_phi->Draw();

    // Draw vertical lines every 15 degrees
    for (int phi = 0; phi <= 360; phi += 15) {
        TLine *line = new TLine(phi, 1400, phi, 9550);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("same");
    }

    // Save the plot
    c3->SaveAs("output/phi_distribution.pdf");

    // Close the ROOT file
    file->Close();

    return 0;
}