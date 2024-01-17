#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <string>

void vertex_study() {
    // Create a canvas and divide it into two subpads (two columns, one row)
    TCanvas *c1 = new TCanvas("c1", "Canvas", 1600, 800);
    c1->Divide(2, 1); // 2 columns, 1 row

    // Increase font size globally
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetLegendTextSize(0.04);

    // Process "Inb" data
    c1->cd(1);
    TPad *pad1 = new TPad("pad1", "Inb Pad", 0, 0, 1, 1);
    pad1->SetBottomMargin(0.15);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();

    // Open files and get TTrees for "Inb" data
    TFile *file_eX_inb = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root");
    TTree *tree_eX_inb = (TTree*)file_eX_inb->Get("PhysicsEvents");
    TFile *file_epi_X_inb = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root");
    TTree *tree_epi_X_inb = (TTree*)file_epi_X_inb->Get("PhysicsEvents");
    TFile *file_ek_X_inb = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root");
    TTree *tree_ek_X_inb = (TTree*)file_ek_X_inb->Get("PhysicsEvents");

    // Create histograms for "Inb" data
    TH1D *h_eX_inb = new TH1D("h_eX_inb", "RGA Fa18 Inb;v_{z} (cm);Normalized counts", 100, -10, 5);
    TH1D *h_epi_X_inb = new TH1D("h_epi_X_inb", "", 100, -10, 5);
    TH1D *h_ek_X_inb = new TH1D("h_ek_X_inb", "", 100, -10, 5);

    h_eX_inb->SetLineColor(kBlack);
    h_epi_X_inb->SetLineColor(kRed);
    h_ek_X_inb->SetLineColor(kBlue);

    // Fill histograms for "Inb" data
    tree_eX_inb->Draw("vz_e>>h_eX_inb", "", "goff");
    tree_epi_X_inb->Draw("vz_p>>h_epi_X_inb", "", "goff");
    tree_ek_X_inb->Draw("vz_p>>h_ek_X_inb", "", "goff");

    // Normalize histograms for "Inb" data
    h_eX_inb->Scale(1.0 / h_eX_inb->Integral());
    h_epi_X_inb->Scale(1.0 / h_epi_X_inb->Integral());
    h_ek_X_inb->Scale(1.0 / h_ek_X_inb->Integral());

    // Draw histograms for "Inb" data
    h_eX_inb->Draw("HIST");
    h_epi_X_inb->Draw("HIST same");
    h_ek_X_inb->Draw("HIST same");

    // Add legend for "Inb" data
    TLegend *leg_inb = new TLegend(0.85, 0.7, 0.9, 0.9);
    leg_inb->AddEntry(h_eX_inb, "e^{-}", "l");
    leg_inb->AddEntry(h_epi_X_inb, "#pi^{-}", "l");
    leg_inb->AddEntry(h_ek_X_inb, "k^{-}", "l");
    leg_inb->Draw();

    // Process "Out" data
    c1->cd(2);
    TPad *pad2 = new TPad("pad2", "Out Pad", 0, 0, 1, 1);
    pad2->SetBottomMargin(0.15);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();

    // Open files and get TTrees for "Out" data
    TFile *file_eX_out = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_eX.root");
    TTree *tree_eX_out = (TTree*)file_eX_out->Get("PhysicsEvents");
    TFile *file_epi_X_out = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_epi-X.root");
    TTree *tree_epi_X_out = (TTree*)file_epi_X_out->Get("PhysicsEvents");
    TFile *file_ek_X_out = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_ek-X.root");
    TTree *tree_ek_X_out = (TTree*)file_ek_X_out->Get("PhysicsEvents");

    // Create histograms for "Out" data
    TH1D *h_eX_out = new TH1D("h_eX_out", "RGA Fa18 Out;v_{z} (cm);Normalized counts", 100, -10, 5);
    TH1D *h_epi_X_out = new TH1D("h_epi_X_out", "", 100, -10, 5);
    TH1D *h_ek_X_out = new TH1D("h_ek_X_out", "", 100, -10, 5);

    h_eX_out->SetLineColor(kBlack);
    h_epi_X_out->SetLineColor(kRed);
    h_ek_X_out->SetLineColor(kBlue);

    // Fill histograms for "Out" data
    tree_eX_out->Draw("vz_e>>h_eX_out", "", "goff");
    tree_epi_X_out->Draw("vz_p>>h_epi_X_out", "", "goff");
    tree_ek_X_out->Draw("vz_p>>h_ek_X_out", "", "goff");

    // Normalize histograms for "Out" data
    h_eX_out->Scale(1.0 / h_eX_out->Integral());
    h_epi_X_out->Scale(1.0 / h_epi_X_out->Integral());
    h_ek_X_out->Scale(1.0 / h_ek_X_out->Integral());

    // Draw histograms for "Out" data
    h_eX_out->Draw("HIST");
    h_epi_X_out->Draw("HIST same");
    h_ek_X_out->Draw("HIST same");

    // Add legend for "Out" data
    TLegend *leg_out = new TLegend(0.85, 0.7, 0.9, 0.9);
    leg_out->AddEntry(h_eX_out, "e^{-}", "l");
    leg_out->AddEntry(h_epi_X_out, "#pi^{-}", "l");
    leg_out->AddEntry(h_ek_X_out, "k^{-}", "l");
    leg_out->Draw();

    // Close files
    file_eX_inb->Close();
    file_epi_X_inb->Close();
    file_ek_X_inb->Close();
    file_eX_out->Close();
    file_epi_X_out->Close();
    file_ek_X_out->Close();

    // Save the canvas
    c1->SaveAs("output/all_vertices.png");
}
