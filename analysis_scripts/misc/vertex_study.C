// Include necessary ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>

// Function to draw histograms for a given panel
void DrawHistogramsForPanel(const char* file_eX, const char* file_epiX, const char* file_ekX, TPad* pad, const char* title) {
    // Load files and trees
    TFile* fileEX = new TFile(file_eX);
    TTree* treeEX = (TTree*)fileEX->Get("PhysicsEvents");
    TFile* fileEpiX = new TFile(file_epiX);
    TTree* treeEpiX = (TTree*)fileEpiX->Get("PhysicsEvents");
    TFile* fileEkX = new TFile(file_ekX);
    TTree* treeEkX = (TTree*)fileEkX->Get("PhysicsEvents");

    // Create histograms
    TH1F* h_eX = new TH1F("h_eX", title, 100, -15, 10);
    h_eX->SetLineColor(kBlack);
    TH1F* h_epiX = new TH1F("h_epiX", title, 100, -15, 10);
    h_epiX->SetLineColor(kRed);
    TH1F* h_ekX = new TH1F("h_ekX", title, 100, -15, 10);
    h_ekX->SetLineColor(kBlue);

    // Fill histograms
    treeEX->Draw("vz_e>>h_eX");
    treeEpiX->Draw("vz_p>>h_epiX");
    treeEkX->Draw("vz_p>>h_ekX");

    // Normalize histograms
    h_eX->Scale(1.0 / h_eX->Integral());
    h_epiX->Scale(1.0 / h_epiX->Integral());
    h_ekX->Scale(1.0 / h_ekX->Integral());

    // Draw histograms on the pad
    pad->cd();
    h_eX->Draw("HIST");
    h_epiX->Draw("HIST SAME");
    h_ekX->Draw("HIST SAME");

    // Add a legend
    TLegend* legend = new TLegend(0.7, 0.75, 0.35, 0.9);
    legend->SetTextSize(0.04); // Increase font size
    legend->AddEntry(h_eX, "e^{-}", "l");
    legend->AddEntry(h_epiX, "#pi^{-}", "l");
    legend->AddEntry(h_ekX, "k^{-}", "l");
    legend->Draw();

    // Style the histograms
    h_eX->GetXaxis()->SetTitle("v_{z} (cm)");
    h_eX->GetYaxis()->SetTitle("normalized counts");
    h_eX->GetXaxis()->CenterTitle();
    h_eX->GetYaxis()->CenterTitle();
    h_eX->GetXaxis()->SetTitleSize(0.05);
    h_eX->GetYaxis()->SetTitleSize(0.05);

    // Remove the stat box
    gStyle->SetOptStat(0);
}

void vertex_study() {
    // Create a canvas with multiple pads
    TCanvas *c1 = new TCanvas("c1", "Vertex Study", 1200, 800);
    c1->Divide(3, 2); // 3 columns, 2 rows

    // RGA Fa18 Inb
    DrawHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root",
        (TPad*)c1->cd(1),
        "RGA Fa18 Inb"
    );

    // RGA Fa18 Out
    DrawHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_ek-X.root",
        (TPad*)c1->cd(2),
        "RGA Fa18 Out"
    );

    // RGA Sp19 Inb
    DrawHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_ek-X.root",
        (TPad*)c1->cd(3),
        "RGA Sp19 Inb"
    );

    // RGB Sp19 Inb
    DrawHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_ek-X.root",
        (TPad*)c1->cd(4),
        "RGB Sp19 Inb"
    );

    // RGB Fa19 Out
    DrawHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_ek-X.root",
        (TPad*)c1->cd(5),
        "RGB Fa19 Out"
    );

    // RGB Sp20 Inb
    DrawHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_ek-X.root",
        (TPad*)c1->cd(6),
        "RGB Sp20 Inb"
    );

    // Save the canvas as a PNG file
    c1->SaveAs("output/neg_vz.png");
}

int main() {
    vertex_study();
    return 0;
}
