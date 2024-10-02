// Importing necessary ROOT libraries
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void plot_zeta_vs_xi() {
    // Open the ROOT files
    TFile *file1 = TFile::Open("/home/thayward/zeta1.root");
    TFile *file2 = TFile::Open("/home/thayward/zeta2.root");
    TFile *file3 = TFile::Open("/home/thayward/zeta3.root");

    // Get the trees from the files
    TTree *tree1 = (TTree*)file1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)file2->Get("PhysicsEvents");
    TTree *tree3 = (TTree*)file3->Get("PhysicsEvents");

    // Create histograms for zeta
    TH1F *hist1 = new TH1F("hist1", "Zeta Distribution", 100, 0, 1);
    TH1F *hist2 = new TH1F("hist2", "Xi Distribution", 100, 0, 1);
    TH1F *hist3 = new TH1F("hist2", "Xi_h Distribution", 100, 0, 1);

    // Set up the branch variables
    double zeta, xF;
    tree1->SetBranchAddress("zeta", &zeta);
    tree1->SetBranchAddress("xF", &xF);
    tree2->SetBranchAddress("zeta", &zeta);
    tree2->SetBranchAddress("xF", &xF);
    tree3->SetBranchAddress("zeta", &zeta);
    tree3->SetBranchAddress("xF", &xF);

    // Fill the histograms for tree1 (zeta1.root)
    Long64_t nentries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nentries1; i++) {
        tree1->GetEntry(i);
        if (xF < 0) hist1->Fill(zeta);
    }

    // Fill the histograms for tree2 (zeta2.root)
    Long64_t nentries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nentries2; i++) {
        tree2->GetEntry(i);
        if (xF < 0) hist2->Fill(zeta);
    }

    // Fill the histograms for tree3 (zeta3.root)
    Long64_t nentries3 = tree3->GetEntries();
    for (Long64_t i = 0; i < nentries3; i++) {
        tree3->GetEntry(i);
        if (xF < 0) hist3->Fill(xi);
    }

    // Create a canvas to draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Zeta vs Xi", 800, 600);

    // Set histogram styles
    hist1->SetLineColor(kRed);
    hist2->SetLineColor(kBlue);
    hist2->SetLineColor(kGreen);
    hist1->SetLineWidth(2);
    hist2->SetLineWidth(2);
    hist3->SetLineWidth(2);

    // Axis labels and remove stats box
    hist1->GetXaxis()->SetTitle("");
    hist1->GetYaxis()->SetTitle("Counts");
    hist1->SetStats(0);  // Remove statistics box
    hist2->SetStats(0);  // Remove statistics box
    hist3->SetStats(0);  // Remove statistics box

    // Draw the histograms
    hist1->Draw();
    hist2->Draw("SAME");
    hist3->Draw("SAME");

    // Create a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, "#zeta", "l");
    legend->AddEntry(hist2, "#xi", "l");
    legend->AddEntry(hist3, "#xi_h", "l");
    legend->Draw();

    // Save the plot as a PNG file
    c1->SaveAs("/home/thayward/zeta_vs_xi.png");

    // Clean up
    delete c1;
    delete hist1;
    delete hist2;
    delete hist3;
    file1->Close();
    file2->Close();
    file3->Close();
}