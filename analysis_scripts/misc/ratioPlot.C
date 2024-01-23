#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>

void ratioPlot() {
    // Open the ROOT files and get the trees
    TFile *f1 = new TFile("/scratch/thayward/ratios/epi-X_out.root");
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");

    TFile *f2 = new TFile("/scratch/thayward/ratios/ek-X_out.root");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    // Create histograms for pi- and k- counts
    TH1F *h1 = new TH1F("h1", "pi- counts", 50, 0, 6);
    TH1F *h2 = new TH1F("h2", "k- counts", 50, 0, 6);

    // Fill histograms with the branch variable p_p
    tree1->Draw("p_p>>h1");
    tree2->Draw("p_p>>h2");

    // Create a ratio histogram
    TH1F *hRatio = (TH1F*)h1->Clone("hRatio");
    hRatio->Divide(h2);

    // Style the plot
    gStyle->SetOptStat(0);
    hRatio->SetTitle("Negatives Outbending;p (GeV);#pi^{-} /k^{-}");
    hRatio->GetXaxis()->CenterTitle();
    hRatio->GetXaxis()->SetLabelSize(0.04);
    hRatio->GetYaxis()->SetLabelSize(0.04);

    // Draw the ratio plot
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    hRatio->Draw("E");

    // Save the canvas as a PNG file
    c1->SaveAs("output/negative_out_ratio.png");

    // Save the canvas as a .C file
    c1->SaveAs("ratio.C");
}

// This allows the script to be run standalone from the ROOT interpreter
void ratio() {
    ratioPlot();
}
