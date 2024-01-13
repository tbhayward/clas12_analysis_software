#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <iostream>

const std::string output_dir = "output/rgc_ready_for_cooking_plots/";

// Function to create a histogram for a given dataset
TH1D* createHistogram(TTree* tree, const char* name, const char* title, 
        const char* variable, const char* cut, double norm) {
    TH1D* hist = new TH1D(name, title, 80, -4, 4);
    tree->Draw((std::string(variable) + ">>" + name).c_str(), cut, "goff");
    hist->Scale(1.0 / norm);
    hist->SetStats(kFALSE);
    return hist;
}

void rgc_preparation() {
    
    const char* files[] = {
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_eX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epi+X.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epi+pi-X.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_eX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epi+X.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epi+pi-X.root"
    };

    // Open files and get trees
    TFile* filesOpened[8];
    TTree* trees[8];
    for (int i = 0; i < 8; i++) {
        filesOpened[i] = new TFile(files[i]);
        trees[i] = (TTree*)filesOpened[i]->Get("PhysicsEvents");
    }

    double rga_H2_norm = 159661.55+145813.73;
    double rgc_pos_NH3_norm = 19355.9+19392.53+21683.25+21621.178;
    double rgc_neg_NH3_norm = 21282.264+21217.414+21303.576+21297.766;
    double rgc_NH3_norm = rgc_pos_NH3_norm+rgc_neg_NH3_norm;
    double rgc_C_norm = 8883.014+8834.256;

    // Create canvas with 4x2 pads
    TCanvas *c1 = new TCanvas("c1", "Data Analysis", 1600, 1200);
    c1->Divide(2, 4);

    // Histograms for each dataset
    TH1D* hists[8];

    // Set style
    gStyle->SetOptStat(0);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(.5);
    gStyle->SetLabelSize(0.04, "XY");
    gStyle->SetTitleSize(0.04, "XY");

    const char* titles[] = {"eX", "e#pi^{+}X", "epX", "e#pi^{+}#pi^{-}X"};
    const char* variables[] = {"Mx", "Mx", "Mx", "Mx"};
    const char* cuts[] = {"runnum != 16297", "runnum == 16297"};
    double norms[] = {rga_norm, rgc_pos_norm, rgc_neg_norm, rgc_carbon_norm};

    for (int i = 0; i < 4; i++) {
        // Left column plots (NH3, C, H2 distributions)
        c1->cd(i * 2 + 1);
        TPad *pad1 = (TPad*)c1->GetPad(i * 2 + 1);
        pad1->SetBottomMargin(0.15);
        pad1->SetLeftMargin(0.15);

        hists[i] = createHistogram(trees[i], (std::string("h_rg") + titles[i]).c_str(), titles[i], variables[i], cuts[0], norms[i]);
        hists[i]->SetLineColor(kBlue);
        hists[i]->GetXaxis()->SetTitle("M_{X} (GeV)");
        hists[i]->GetYaxis()->SetTitle("Counts / nC");
        hists[i]->Draw();

        hists[i + 4] = createHistogram(trees[i + 4], (std::string("h_rgc") + titles[i]).c_str(), "", variables[i], cuts[1], norms[i + 4]);
        hists[i + 4]->SetLineColor(kGreen);
        hists[i + 4]->Draw("same");

        // Add legend
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(hists[i], "NH3", "l");
        leg->AddEntry(hists[i + 4], "C", "l");
        leg->Draw();

        // Right column plots (NH3/C ratio)
        c1->cd(i * 2 + 2);
        TPad *pad2 = (TPad*)c1->GetPad(i * 2 + 2);
        pad2->SetBottomMargin(0.15);
        pad2->SetLeftMargin(0.15);

        TH1D* ratioHist = (TH1D*)hists[i]->Clone((std::string("ratio_") + titles[i]).c_str());
        ratioHist->Divide(hists[i + 4]);
        ratioHist->SetLineColor(kRed);
        ratioHist->GetXaxis()->SetTitle("M_{X} (GeV)");
        ratioHist->GetYaxis()->SetTitle("NH_{3} / C");
        ratioHist->Draw();
        // Add label for ratio plots
        TLatex *label = new TLatex();
        label->SetTextSize(0.04);
        label->DrawLatexNDC(0.7, 0.8, "NH3/C Ratio");
    }

    // Save the canvas as "normalizations.png"
    std::string final_output = output_dir + "normalizations.png";
    c1->SaveAs(final_output.c_str());

    // Clean up
    for (int i = 0; i < 8; i++) {
        delete filesOpened[i];
    }
    delete c1;

}

int main() {
    rgc_preparation();
    return 0;
}
