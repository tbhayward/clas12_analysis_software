#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>

void misIDPlot() {
    // Open the ROOT file and get the tree
    TFile *file = new TFile("/scratch/thayward/ratios/epi-X_inb.root");
    TTree *tree = (TTree*)file->Get("PhysicsEvents");

    // Create a histogram for the fraction calculation
    TH1F *hFraction = new TH1F("hFraction", ";p (GeV);k^{-} #rightarrow #pi^{-}", 15, 0, 7);

    // Loop over the tree and fill the histogram
    double p_p;
    int matching_p1_pid;
    tree->SetBranchAddress("p_p", &p_p);
    tree->SetBranchAddress("matching_p1_pid", &matching_p1_pid);

    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (matching_p1_pid == 321) {
            hFraction->Fill(p_p);
        }
    }

    // Normalize the histogram to get the fraction
    hFraction->Scale(1.0 / tree->GetEntries());

    // Style the plot
    gStyle->SetOptStat(0);
    hFraction->GetXaxis()->CenterTitle();
    hFraction->GetXaxis()->SetLabelSize(0.04);
    hFraction->GetYaxis()->SetLabelSize(0.04);
    hFraction->SetMarkerStyle(20);
    hFraction->SetMarkerSize(1.2);

    // Draw the plot
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetLeftMargin(0.15);
    hFraction->Draw("E");

    // Save the canvas as a PNG file
    c1->SaveAs("output/epi-X_misid.png");
}

// This allows the script to be run standalone from the ROOT interpreter
void misID() {
    misIDPlot();
}
