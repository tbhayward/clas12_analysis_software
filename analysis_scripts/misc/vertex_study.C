// Include necessary ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>

// Function to draw histograms for a given panel
void DrawNegHistogramsForPanel(const char* file_eX, const char* file_epiX, const char* file_ekX, TPad* pad, const char* title) {
    pad->SetLeftMargin(0.15); // Increase left margin
    pad->SetBottomMargin(0.15); // Increase bottom margin

    // Load files and trees
    TFile* fileEX = new TFile(file_eX);
    TTree* treeEX = (TTree*)fileEX->Get("PhysicsEvents");
    TFile* fileEpiX = new TFile(file_epiX);
    TTree* treeEpiX = (TTree*)fileEpiX->Get("PhysicsEvents");
    TFile* fileEkX = new TFile(file_ekX);
    TTree* treeEkX = (TTree*)fileEkX->Get("PhysicsEvents");

    // Create histograms
    TH1F* h_eX = new TH1F("h_eX", title, 100, -15, 1);
    h_eX->SetLineColor(kBlack);
    TH1F* h_epiX = new TH1F("h_epiX", title, 100, -15, 1);
    h_epiX->SetLineColor(kRed);
    TH1F* h_ekX = new TH1F("h_ekX", title, 100, -15, 1);
    h_ekX->SetLineColor(kBlue);

    // Fill histograms
    treeEX->Draw("vz_e>>h_eX");
    treeEpiX->Draw("vz_p>>h_epiX");
    treeEkX->Draw("vz_p>>h_ekX");

    // Normalize histograms
    h_eX->Scale(1.0 / h_eX->Integral());
    h_epiX->Scale(1.0 / h_epiX->Integral());
    h_ekX->Scale(1.0 / h_ekX->Integral());

    // Find the maximum value among the histograms
    double maxVal = h_eX->GetMaximum();

    // Set y-axis to 20% higher than the largest maximum
    double maxYAxis = maxVal * 1.3;
    h_eX->SetMaximum(maxYAxis);
    h_epiX->SetMaximum(maxYAxis); // This might be redundant but ensures consistency
    h_ekX->SetMaximum(maxYAxis); // This might be redundant but ensures consistency

    // Draw histograms on the pad
    pad->cd();
    h_eX->Draw("HIST");
    // h_epiX->Draw("HIST SAME");
    // h_ekX->Draw("HIST SAME");

    // Calculate mean and standard deviation for each histogram
    double meanEpiX = h_eX->GetMean();
    double stdEpiX = h_eX->GetStdDev();

    // Add a legend with mean and std
    TLegend* legend = new TLegend(0.4, 0.9, 0.9, 0.8);
    legend->SetTextSize(0.04); // Increase font size
    TString legendEntryEpiX = Form("e^{-}, #mu = %.2f, #sigma = %.2f", meanEpiX, stdEpiX);
    legend->AddEntry(h_eX, legendEntryEpiX, "l");
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

void DrawDiffNegHistogramsForPanel(const char* file_epiX, const char* file_epX, 
        TPad* pad, const char* title, double vz_min, double vz_max) {
    pad->SetLeftMargin(0.15); // Increase left margin
    pad->SetBottomMargin(0.15); // Increase bottom margin

    // Load files and trees
    TFile* fileEpiX = new TFile(file_epiX);
    TTree* treeEpiX = (TTree*)fileEpiX->Get("PhysicsEvents");
    TFile* fileEpX = new TFile(file_epX);
    TTree* treeEpX = (TTree*)fileEpX->Get("PhysicsEvents");

    // Create histograms for vz_e - vz_p
    TH1F* h_diffEpiX = new TH1F("h_diffEpiX", title, 100, -8, 8);
    h_diffEpiX->SetLineColor(kRed);
    TH1F* h_diffEpX = new TH1F("h_diffEpX", title, 100, -8, 8);
    h_diffEpX->SetLineColor(kBlue);

    // Construct the draw command with condition
    TString drawCommandEpiX = Form("vz_e-vz_p>>h_diffEpiX", vz_min, vz_max);
    TString drawCommandEpX = Form("vz_e-vz_p>>h_diffEpX", vz_min, vz_max);
    TString condition = Form("%f < vz_e && vz_e < %f", vz_min, vz_max);
    // Fill histograms with the difference vz_e - vz_p, with the condition
    treeEpiX->Draw(drawCommandEpiX + " " + condition);
    treeEpX->Draw(drawCommandEpX + " " + condition);

    // Normalize histograms
    h_diffEpiX->Scale(1.0 / h_diffEpiX->Integral());
    h_diffEpX->Scale(1.0 / h_diffEpX->Integral());

    // Find the maximum value among the histograms
    double maxValEpiX = h_diffEpiX->GetMaximum();
    double maxValEpX = h_diffEpX->GetMaximum();
    double maxVal = TMath::Max(maxValEpiX, maxValEpX);

    // Set y-axis to 20% higher than the largest maximum
    double maxYAxis = maxVal * 1.3;
    h_diffEpiX->SetMaximum(maxYAxis);
    h_diffEpX->SetMaximum(maxYAxis); // This might be redundant but ensures consistency

    // Draw histograms on the pad
    pad->cd();
    h_diffEpiX->Draw("HIST");
    h_diffEpX->Draw("HIST SAME");

    // Calculate mean and standard deviation for each histogram
    double meanEpiX = h_diffEpiX->GetMean();
    double stdEpiX = h_diffEpiX->GetStdDev();
    double meanEpX = h_diffEpX->GetMean();
    double stdEpX = h_diffEpX->GetStdDev();

    // Add a legend with mean and std
    TLegend* legend = new TLegend(0.4, 0.9, 0.9, 0.8);
    legend->SetTextSize(0.04); // Increase font size
    TString legendEntryEpiX = Form("#pi^{-}, #mu = %.2f, #sigma = %.2f", meanEpiX, stdEpiX);
    TString legendEntryEpX = Form("k^{-}, #mu = %.2f, #sigma = %.2f", meanEpX, stdEpX);
    legend->AddEntry(h_diffEpiX, legendEntryEpiX, "l");
    legend->AddEntry(h_diffEpX, legendEntryEpX, "l");
    legend->Draw();

    // Style the histograms
    h_diffEpiX->GetXaxis()->SetTitle("v_{z_{e}} - v_{z_{h}} (cm)");
    h_diffEpiX->GetYaxis()->SetTitle("normalized counts");
    h_diffEpiX->GetXaxis()->CenterTitle();
    h_diffEpiX->GetYaxis()->CenterTitle();
    h_diffEpiX->GetXaxis()->SetTitleSize(0.05);
    h_diffEpiX->GetYaxis()->SetTitleSize(0.05);

    // Remove the stat box
    gStyle->SetOptStat(0);
}

void DrawDiffPosHistogramsForPanel(const char* file_epiX, const char* file_epX, 
        TPad* pad, const char* title, double vz_min, double vz_max) {
    pad->SetLeftMargin(0.15); // Increase left margin
    pad->SetBottomMargin(0.15); // Increase bottom margin

    // Load files and trees
    TFile* fileEpiX = new TFile(file_epiX);
    TTree* treeEpiX = (TTree*)fileEpiX->Get("PhysicsEvents");
    TFile* fileEpX = new TFile(file_epX);
    TTree* treeEpX = (TTree*)fileEpX->Get("PhysicsEvents");

    // Create histograms for vz_e - vz_p
    TH1F* h_diffEpiX = new TH1F("h_diffEpiX", title, 100, -8, 8);
    h_diffEpiX->SetLineColor(kRed);
    TH1F* h_diffEpX = new TH1F("h_diffEpX", title, 100, -8, 8);
    h_diffEpX->SetLineColor(kBlue);

    // Construct the draw command with condition
    TString drawCommandEpiX = Form("vz_e-vz_p>>h_diffEpiX", vz_min, vz_max);
    TString drawCommandEpX = Form("vz_e-vz_p>>h_diffEpX", vz_min, vz_max);
    TString condition = Form("%f < vz_e && vz_e < %f", vz_min, vz_max);
    // Fill histograms with the difference vz_e - vz_p, with the condition
    treeEpiX->Draw(drawCommandEpiX + " " + condition);
    treeEpX->Draw(drawCommandEpX + " " + condition);

    // Normalize histograms
    h_diffEpiX->Scale(1.0 / h_diffEpiX->Integral());
    h_diffEpX->Scale(1.0 / h_diffEpX->Integral());

    // Find the maximum value among the histograms
    double maxValEpiX = h_diffEpiX->GetMaximum();
    double maxValEpX = h_diffEpX->GetMaximum();
    double maxVal = TMath::Max(maxValEpiX, maxValEpX);

    // Set y-axis to 20% higher than the largest maximum
    double maxYAxis = maxVal * 1.3;
    h_diffEpiX->SetMaximum(maxYAxis);
    h_diffEpX->SetMaximum(maxYAxis); // This might be redundant but ensures consistency

    // Draw histograms on the pad
    pad->cd();
    h_diffEpiX->Draw("HIST");
    h_diffEpX->Draw("HIST SAME");

    // Calculate mean and standard deviation for each histogram
    double meanEpiX = h_diffEpiX->GetMean();
    double stdEpiX = h_diffEpiX->GetStdDev();
    double meanEpX = h_diffEpX->GetMean();
    double stdEpX = h_diffEpX->GetStdDev();

    // Add a legend with mean and std
    TLegend* legend = new TLegend(0.4, 0.9, 0.9, 0.8);
    legend->SetTextSize(0.04); // Increase font size
    TString legendEntryEpiX = Form("#pi^{+}, #mu = %.2f, #sigma = %.2f", meanEpiX, stdEpiX);
    TString legendEntryEpX = Form("p, #mu = %.2f, #sigma = %.2f", meanEpX, stdEpX);
    legend->AddEntry(h_diffEpiX, legendEntryEpiX, "l");
    legend->AddEntry(h_diffEpX, legendEntryEpX, "l");
    legend->Draw();

    // Style the histograms
    h_diffEpiX->GetXaxis()->SetTitle("v_{z_{e}} - v_{z_{h}} (cm)");
    h_diffEpiX->GetYaxis()->SetTitle("normalized counts");
    h_diffEpiX->GetXaxis()->CenterTitle();
    h_diffEpiX->GetYaxis()->CenterTitle();
    h_diffEpiX->GetXaxis()->SetTitleSize(0.05);
    h_diffEpiX->GetYaxis()->SetTitleSize(0.05);

    // Remove the stat box
    gStyle->SetOptStat(0);
}


void vertex_study() {
    // Create a canvas with multiple pads
    TCanvas *c_neg = new TCanvas("c_neg", "Vertex Study", 1200, 800);
    c_neg->Divide(3, 2); // 3 columns, 2 rows
    // RGA Fa18 Inb
    DrawNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root",
        (TPad*)c_neg->cd(1), "RGA Fa18 Inb"
    );
    // RGA Fa18 Out
    DrawNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_ek-X.root",
        (TPad*)c_neg->cd(2), "RGA Fa18 Out"
    );
    // RGA Sp19 Inb
    DrawNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_ek-X.root",
        (TPad*)c_neg->cd(3), "RGA Sp19 Inb"
    );
    // RGB Sp19 Inb
    DrawNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_ek-X.root",
        (TPad*)c_neg->cd(4), "RGB Sp19 Inb"
    );
    // RGB Fa19 Out
    DrawNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_ek-X.root",
        (TPad*)c_neg->cd(5), "RGB Fa19 Out"
    );
    // RGB Sp20 Inb
    DrawNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_eX.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_ek-X.root",
        (TPad*)c_neg->cd(6), "RGB Sp20 Inb"
    );
    // Save the canvas as a PNG file
    c_neg->SaveAs("output/neg_vz.png");

    // Create a canvas with multiple pads
    TCanvas *c_neg_diff = new TCanvas("c_neg", "Vertex Study", 1200, 800);
    c_neg_diff->Divide(3, 2); // 3 columns, 2 rows
    // RGA Fa18 Inb
    DrawDiffNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root",
        (TPad*)c_neg_diff->cd(1), "RGA Fa18 Inb", -6, 1
    );
    // RGA Fa18 Out
    DrawDiffNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_ek-X.root",
        (TPad*)c_neg_diff->cd(2), "RGA Fa18 Out", -7, 0
    );
    // RGA Sp19 Inb
    DrawDiffNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_ek-X.root",
        (TPad*)c_neg_diff->cd(3), "RGA Sp19 Inb", -6, 1
    );
    // RGB Sp19 Inb
    DrawDiffNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_ek-X.root",
        (TPad*)c_neg_diff->cd(4), "RGB Sp19 Inb", -6, 1
    );
    // RGB Fa19 Out
    DrawDiffNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_ek-X.root",
        (TPad*)c_neg_diff->cd(5), "RGB Fa19 Out", -7, 0
    );
    // RGB Sp20 Inb
    DrawDiffNegHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_epi-X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_ek-X.root",
        (TPad*)c_neg_diff->cd(6), "RGB Sp20 Inb", -6, 1
    );
    // Save the canvas as a PNG file
    c_neg_diff->SaveAs("output/neg_diff_vz.png");

    // Create a canvas with multiple pads
    TCanvas *c_pos_diff = new TCanvas("c_neg", "Vertex Study", 1200, 800);
    c_pos_diff->Divide(3, 2); // 3 columns, 2 rows
    // RGA Fa18 Inb
    DrawDiffPosHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi+X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epX.root",
        (TPad*)c_pos_diff->cd(1), "RGA Fa18 Inb", -6, 1
    );
    // RGA Fa18 Out
    DrawDiffPosHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_epi+X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/fa18_out/rga_fa18_out_epX.root",
        (TPad*)c_pos_diff->cd(2), "RGA Fa18 Out", -7, 0
    );
    // RGA Sp19 Inb
    DrawDiffPosHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_epi+X.root",
        "/volatile/clas12/thayward/vertex_studies/rga/sp19_inb/rga_sp19_inb_epX.root",
        (TPad*)c_pos_diff->cd(3), "RGA Sp19 Inb", -6, 1
    );
    // RGB Sp19 Inb
    DrawDiffPosHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_epi+X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp19_inb/rgb_sp19_inb_epX.root",
        (TPad*)c_pos_diff->cd(4), "RGB Sp19 Inb", -6, 1
    );
    // RGB Fa19 Out
    DrawDiffPosHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_epi+X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/fa19_out/rgb_fa19_out_epX.root",
        (TPad*)c_pos_diff->cd(5), "RGB Fa19 Out", -7, 0
    );
    // RGB Sp20 Inb
    DrawDiffPosHistogramsForPanel(
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_epi+X.root",
        "/volatile/clas12/thayward/vertex_studies/rgb/sp20_inb/rgb_sp20_inb_epX.root",
        (TPad*)c_pos_diff->cd(6), "RGB Sp20 Inb", -6, 1
    );
    // Save the canvas as a PNG file
    c_pos_diff->SaveAs("output/pos_diff_vz.png");
}

int main() {
    vertex_study();
    return 0;
}
