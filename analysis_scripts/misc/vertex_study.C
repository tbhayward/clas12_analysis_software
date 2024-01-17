#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

void vertex_study() {
    // Create a canvas to draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.15);

    // Open the ROOT files and get the TTrees
    TFile *file1 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root");
    TTree *tree1 = (TTree*)file1->Get("PhysicsEvents");
    TFile *file2 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root");
    TTree *tree2 = (TTree*)file2->Get("PhysicsEvents");
    TFile *file3 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root");
    TTree *tree3 = (TTree*)file3->Get("PhysicsEvents");

    // Create histograms
    TH1D *h1 = new TH1D("h1", "RGA Fa18 Inb;v_{z};Normalized counts", 100, -15, 5);
    TH1D *h2 = new TH1D("h2", "", 100, -15, 5);
    TH1D *h3 = new TH1D("h3", "", 100, -15, 5);

    // Set histogram colors
    h1->SetLineColor(kBlack);
    h2->SetLineColor(kRed);
    h3->SetLineColor(kBlue);

    // Fill the histograms
    tree1->Draw("vz_e>>h1");
    tree2->Draw("vz_p>>h2");
    tree3->Draw("vz_p>>h3");

    // Normalize histograms
    h1->Scale(1.0 / h1->Integral());
    h2->Scale(1.0 / h2->Integral());
    h3->Scale(1.0 / h3->Integral());

    // Draw the histograms on the same canvas
    h1->Draw("HIST");
    h2->Draw("HIST same");
    h3->Draw("HIST same");

    // Add a legend
    TLegend *leg = new TLegend(0.8, 0.7, 0.9, 0.9);
    leg->AddEntry(h1, "e^{-}", "l");
    leg->AddEntry(h2, "#pi^{-}", "l");
    leg->AddEntry(h3, "k^{-}", "l");
    leg->Draw();

    // Remove the stat box from all histograms
    h1->SetStats(0);
    h2->SetStats(0);
    h3->SetStats(0);

    // Save the histogram as an image
    c1->SaveAs("output/neg_vertices.png");

    // Close the files
    file1->Close();
    file2->Close();
    file3->Close();
}
