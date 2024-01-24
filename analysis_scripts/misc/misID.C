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
    TH1F *hFraction = new TH1F("hFraction", ";p (GeV);% k^{-} #rightarrow #pi^{-}", 12, 0, 7);

    // Loop over the tree and fill the histogram
    double p_p;
    int matching_p1_pid;
    tree->SetBranchAddress("p_p", &p_p);
    tree->SetBranchAddress("matching_p1_pid", &matching_p1_pid);

    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (matching_p1_pid == -321) {
            hFraction->Fill(p_p);
        }
    }

    // Normalize the histogram to get the fraction
    hFraction->Scale(100.0 / tree->GetEntries());

    // Create a TGraphErrors from the histogram with only vertical error bars
    TGraphErrors *graph = new TGraphErrors();
    for (int i = 1; i <= hFraction->GetNbinsX(); ++i) {
        graph->SetPoint(i-1, hFraction->GetBinCenter(i), hFraction->GetBinContent(i));
        graph->SetPointError(i-1, 0, hFraction->GetBinError(i)); // Set horizontal error to 0
    }

    // Style the graph with markers
    graph->SetTitle(";p (GeV);% k^{-} #rightarrow #pi^{-}"); // Setting title and axis labels
    graph->SetMarkerStyle(20);  // Style 20 is a filled circle
    graph->SetMarkerSize(1.2);  // Adjust the size as needed

    // Draw the plot using TGraphErrors
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetLeftMargin(0.15);
    graph->Draw("AP"); // "AP" to draw the graph with markers and lines

    // Set axis styles
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetLabelSize(0.06);
    graph->GetYaxis()->SetLabelSize(0.06);
    graph->GetYaxis()->SetTitleOffset(1.5); // Adjust if necessary

    // Save the canvas as a PNG file
    c1->SaveAs("output/epi-X_misid.png");
}

// This allows the script to be run standalone from the ROOT interpreter
void misID() {
    misIDPlot();
}
