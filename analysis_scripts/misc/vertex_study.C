#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>

void vertex_study() {
    // Create a canvas to draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.15);

    // Open the first ROOT file
    TFile *file1 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root");
    TTree *tree1 = (TTree*)file1->Get("PhysicsEvents");

    // Open the second ROOT file
    TFile *file2 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root");
    TTree *tree2 = (TTree*)file2->Get("PhysicsEvents");

    // Open the third ROOT file
    TFile *file3 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root");
    TTree *tree3 = (TTree*)file3->Get("PhysicsEvents");

    // Create histograms
    TH1D *h1 = new TH1D("h1", "v_{z} distribution;v_{z};counts", 100, -15, 5);
    TH1D *h2 = new TH1D("h2", "v_{z} distribution", 100, -15, 5);
    TH1D *h3 = new TH1D("h3", "v_{z} distribution", 100, -15, 5);

    // Set histogram colors
    h1->SetLineColor(kBlack);
    h2->SetLineColor(kRed);
    h3->SetLineColor(kBlue);

    // Fill the histograms
    tree1->Draw("vz_e>>h1");
    tree2->Draw("vz_p>>h2");
    tree3->Draw("vz_p>>h3");

    // Draw the histograms on the same canvas
    h1->Draw();
    h2->Draw("same");
    h3->Draw("same");

    // Customize histogram appearance
    h1->GetXaxis()->CenterTitle(true);
    h1->GetYaxis()->CenterTitle(true);
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleSize(0.05);

    // Add a legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h1, "e^{-}", "l");
    leg->AddEntry(h2, "#pi^{-}", "l");
    leg->AddEntry(h3, "k^{-}", "l");
    leg->Draw();

    // Add plot label
    TText *label = new TText(0.2, 0.85, "RGA Fa18 Inb");
    label->SetTextSize(0.04);
    label->Draw();

    // Remove the stat box from all histograms
    h1->SetStats(0);
    h2->SetStats(0);
    h3->SetStats(0);

    // Save the histogram as an image in the specified output directory
    c1->SaveAs("output/neg_vertices.png");

    // Close the files
    file1->Close();
    file2->Close();
    file3->Close();
}
