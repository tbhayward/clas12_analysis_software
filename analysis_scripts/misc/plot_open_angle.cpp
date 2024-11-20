#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>

void plot_open_angle() {
    // Open the input ROOT file
    TFile *file = TFile::Open("/volatile/clas12/thayward/compass_rad_workshop/rga_sp19_inb_egammaX.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file!" << std::endl;
        return;
    }

    // Access the TTree
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: Cannot find PhysicsEvents tree!" << std::endl;
        file->Close();
        return;
    }

    // Create a histogram for the branch 'open_angle'
    TH1D *hist = new TH1D("hist", ";#theta_{e'#gamma};Normalized Counts", 100, 0, 180); // Adjust range/bins as needed
    tree->Draw("open_angle>>hist", "", "goff");

    // Normalize the histogram
    if (hist->Integral() > 0)
        hist->Scale(1.0 / hist->Integral());

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Open Angle Plot", 800, 600);
    canvas->SetLeftMargin(0.125); // Left margin padding
    gStyle->SetOptStat(0); // Disable the stat box

    // Draw the histogram
    hist->SetLineColor(kBlack); // Set histogram color to black
    hist->SetLineWidth(2); // Set line width
    hist->Draw("HIST");

    // Save the plot
    canvas->SaveAs("output/FSR.pdf");

    // Clean up
    delete canvas;
    delete hist;
    file->Close();
    delete file;

    std::cout << "Plot saved to output/FSR.pdf" << std::endl;
}

// Main function for standalone compilation
int main() {
    plot_open_angle();
    return 0;
}