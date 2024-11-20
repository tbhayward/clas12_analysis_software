#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>

void plot_open_angle() {
    // Open the input ROOT file
    TFile *file = TFile::Open("/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root");
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
    TH1D *hist = new TH1D("hist", ";#theta_{e'#gamma};Normalized Counts", 150, 0, 40); // Adjust range/bins as needed
    tree->Draw("open_angle>>hist", "", "goff");

    // Normalize the histogram
    if (hist->Integral() > 0)
        hist->Scale(1.0 / hist->Integral());

    // Set the y-axis range explicitly
    hist->SetMinimum(1e-3); // Minimum y-axis value
    hist->SetMaximum(1e-1); // Maximum y-axis value

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Open Angle Plot", 800, 600);
    canvas->SetLeftMargin(0.125); // Left margin padding
    gStyle->SetOptStat(0); // Disable the stat box

    // Enable logarithmic scale for y-axis
    canvas->SetLogy();

    // Draw the histogram
    hist->SetLineColor(kBlack); // Set histogram color to black
    hist->SetLineWidth(2); // Set line width
    hist->Draw("HIST");

    // Save the plot
    canvas->SaveAs("output/FSR.pdf");

    // Access the branch 'open_angle' for looping
    Double_t open_angle_ep2 = 0; // Change to Double_t to ensure compatibility
    tree->SetBranchAddress("open_angle_ep2", &open_angle_ep2);

    std::cout << "Open angle values:" << std::endl;

    // Loop through the tree entries and print the values of 'open_angle'
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (open_angle_ep2 != 0) { // Debug: Check if the value is non-zero
            std::cout << "Entry " << i << ": " << open_angle_ep2 << std::endl;
        } else {
            std::cerr << "Entry " << i << ": open_angle is zero or uninitialized!" << std::endl;
        }
    }

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