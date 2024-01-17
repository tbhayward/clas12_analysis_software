#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <vector>
#include <string>

void vertex_study() {
    // List of file paths and titles
    std::vector<std::string> eX_files = {
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_eX.root"
        // Add more file paths as needed
    };
    std::vector<std::string> epi_X_files = {
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_epi-X.root"
        // Add more file paths as needed
    };
    std::vector<std::string> ek_X_files = {
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_ek-X.root"
        // Add more file paths as needed
    };
    std::vector<std::string> titles = {"RGA Fa18 Inb", "RGA Fa18 Out"}; // Add more titles as needed

    // Create a canvas and divide it into two subpads (two columns, one row)
    TCanvas *c1 = new TCanvas("c1", "Canvas", 1600, 800);
    c1->Divide(2, 1); // 2 columns, 1 row

    // Increase font size globally
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetLegendTextSize(0.04);

    // Processing for each dataset
    for (size_t i = 0; i < eX_files.size(); ++i) {
        c1->cd(i+1);
        TPad *pad = new TPad(Form("pad%d", i), Form("Pad%d", i), 0, 0, 1, 1);
        pad->SetBottomMargin(0.15);
        pad->SetLeftMargin(0.15);
        pad->Draw();
        pad->cd();

        // Open files and get TTrees
        TFile *file_eX = TFile::Open(eX_files[i].c_str());
        TTree *tree_eX = (TTree*)file_eX->Get("PhysicsEvents");
        TFile *file_epi_X = TFile::Open(epi_X_files[i].c_str());
        TTree *tree_epi_X = (TTree*)file_epi_X->Get("PhysicsEvents");
        TFile *file_ek_X = TFile::Open(ek_X_files[i].c_str());
        TTree *tree_ek_X = (TTree*)file_ek_X->Get("PhysicsEvents");

        // Create histograms
        TH1D *h_eX = new TH1D(Form("h_eX_%d", i), Form("%s;v_{z} (cm);Normalized counts", titles[i].c_str()), 100, -10, 5);
        TH1D *h_epi_X = new TH1D(Form("h_epi_X_%d", i), "", 100, -10, 5);
        TH1D *h_ek_X = new TH1D(Form("h_ek_X_%d", i), "", 100, -10, 5);

        // Set colors
        h_eX->SetLineColor(kBlack);
        h_epi_X->SetLineColor(kRed);
        h_ek_X->SetLineColor(kBlue);

        // Fill histograms
        tree_eX->Draw(Form("vz_e>>h_eX_%d", i));
        tree_epi_X->Draw(Form("vz_p>>h_epi_X_%d", i));
        tree_ek_X->Draw(Form("vz_p>>h_ek_X_%d", i));

        // Normalize histograms
        h_eX->Scale(1.0 / h_eX->Integral());
        h_epi_X->Scale(1.0 / h_epi_X->Integral());
        h_ek_X->Scale(1.0 / h_ek_X->Integral());

        // Draw histograms
        h_eX->Draw("HIST");
        h_epi_X->Draw("HIST same");
        h_ek_X->Draw("HIST same");

        // Add legend
        TLegend *leg = new TLegend(0.85, 0.7, 0.9, 0.9);
        leg->AddEntry(h_eX, "e^{-}", "l");
        leg->AddEntry(h_epi_X, "#pi^{-}", "l");
        leg->AddEntry(h_ek_X, "k^{-}", "l");
        leg->Draw();

        // Close files
        file_eX->Close();
        file_epi_X->Close();
        file_ek_X->Close();
    }

    // Save the canvas
    c1->SaveAs("output/all_vertices.png");
}

