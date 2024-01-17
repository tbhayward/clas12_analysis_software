#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>

void vertex_study() {
    // Create a canvas and divide it into two subpads
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 1200);
    c1->Divide(1, 2); // 1 column, 2 rows

    // First subpad (Inb data)
    c1->cd(1);
    TPad *pad1 = new TPad("pad1", "Pad1", 0.0, 0.5, 1.0, 1.0);
    pad1->SetBottomMargin(0.15);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();

    // Open the ROOT files for Inb data and get the TTrees
    TFile *file1 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root");
    TTree *tree1 = (TTree*)file1->Get("PhysicsEvents");
    TFile *file2 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root");
    TTree *tree2 = (TTree*)file2->Get("PhysicsEvents");
    TFile *file3 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root");
    TTree *tree3 = (TTree*)file3->Get("PhysicsEvents");

    // Create histograms for Inb data
    TH1D *h1 = new TH1D("h1", "RGA Fa18 Inb;v_{z} (cm);Normalized counts", 100, -10, 5);
    TH1D *h2 = new TH1D("h2", "", 100, -10, 5);
    TH1D *h3 = new TH1D("h3", "", 100, -10, 5);

    // Set histogram colors for Inb data
    h1->SetLineColor(kBlack);
    h2->SetLineColor(kRed);
    h3->SetLineColor(kBlue);

    // Fill the histograms for Inb data
    tree1->Draw("vz_e>>h1");
    tree2->Draw("vz_p>>h2");
    tree3->Draw("vz_p>>h3");

    // Normalize histograms for Inb data
    h1->Scale(1.0 / h1->Integral());
    h2->Scale(1.0 / h2->Integral());
    h3->Scale(1.0 / h3->Integral());

    // Draw the histograms on the first subpad
    h1->Draw("HIST");
    h2->Draw("HIST same");
    h3->Draw("HIST same");

    // Add a legend for Inb data
    TLegend *leg1 = new TLegend(0.85, 0.7, 0.9, 0.9);
    leg1->AddEntry(h1, "e^{-}", "l");
    leg1->AddEntry(h2, "#pi^{-}", "l");
    leg1->AddEntry(h3, "k^{-}", "l");
    leg1->Draw();

    // Remove the stat box from all Inb histograms
    h1->SetStats(0);
    h2->SetStats(0);
    h3->SetStats(0);

    // Second subpad (Out data)
    c1->cd(2);
    TPad *pad2 = new TPad("pad2", "Pad2", 0.0, 0.0, 1.0, 0.5);
    pad2->SetBottomMargin(0.15);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();

    // Open the ROOT files for Out data and get the TTrees
    TFile *file4 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_eX.root");
    TTree *tree4 = (TTree*)file4->Get("PhysicsEvents");
    TFile *file5 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_epi-X.root");
    TTree *tree5 = (TTree*)file5->Get("PhysicsEvents");
    TFile *file6 = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_ek-X.root");
    TTree *tree6 = (TTree*)file6->Get("PhysicsEvents");

    // Create histograms for Out data
    TH1D *h4 = new TH1D("h4", "RGA Fa18 Out;v_{z} (cm);Normalized counts", 100, -10, 5);
    TH1D *h5 = new TH1D("h5", "", 100, -10, 5);
    TH1D *h6 = new TH1D("h6", "", 100, -10, 5);

    // Set histogram colors for Out data
    h4->SetLineColor(kBlack);
    h5->SetLineColor(kRed);
    h6->SetLineColor(kBlue);

    // Fill the histograms for Out data
    tree4->Draw("vz_e>>h4");
    tree5->Draw("vz_p>>h5");
    tree6->Draw("vz_p>>h6");

    // Normalize histograms for Out data
    h4->Scale(1.0 / h4->Integral());
    h5->Scale(1.0 / h5->Integral());
    h6->Scale(1.0 / h6->Integral());

    // Draw the histograms on the second subpad
    h4->Draw("HIST");
    h5->Draw("HIST same");
    h6->Draw("HIST same");

    // Add a legend for Out data
    TLegend *leg2 = new TLegend(0.85, 0.7, 0.9, 0.9);
    leg2->AddEntry(h4, "e^{-}", "l");
    leg2->AddEntry(h5, "#pi^{-}", "l");
    leg2->AddEntry(h6, "k^{-}", "l");
    leg2->Draw();

    // Remove the stat box from all Out histograms
    h4->SetStats(0);
    h5->SetStats(0);
    h6->SetStats(0);

    // Save the entire canvas with both subpads
    c1->cd();
    c1->SaveAs("output/all_vertices.png");

    // Close all files
    file1->Close();
    file2->Close();
    file3->Close();
    file4->Close();
    file5->Close();
    file6->Close();
}
