#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <iostream>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h> // For TMath::Max

const std::string output_dir = "output/rgc_ready_for_cooking_plots/";

// Function to create a histogram for a given dataset
TH1D* createHistogram(TTree* tree, const char* name, const char* title, 
                      const char* variable, const char* cut, double norm) {
    TH1D* hist = new TH1D(name, title, 100, -2, 4);
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

    TCanvas *c1 = new TCanvas("c1", "Data Analysis", 1600, 1200);
    c1->Divide(2, 4);

    TH1D* hists[12]; // 4 plots * 3 histograms per plot

    gStyle->SetOptStat(0);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(.5);
    gStyle->SetLabelSize(0.05, "XY");
    gStyle->SetTitleSize(0.05, "XY");

    const char* titles[] = {"eX", "e#pi^{+}X", "epX", "e#pi^{+}#pi^{-}X"};
    const char* variables[] = {"Mx", "Mx", "Mx", "Mx"};
    const char* cuts_NH3[] = {"runnum != 16297", "runnum != 16297", "runnum != 16297", "runnum != 16297"};
    const char* cuts_C[] = {"runnum == 16297", "runnum == 16297", "runnum == 16297", "runnum == 16297"};

    for (int i = 0; i < 4; i++) {
        c1->cd(i * 2 + 1);
        TPad* pad1 = (TPad*)c1->GetPad(i * 2 + 1);
        pad1->SetBottomMargin(0.15);
        pad1->SetLeftMargin(0.15);

        // Adjust variable and range for the eX plot
        const char* var = (i == 0) ? "x" : "Mx"; // Use "x" for the first plot
        double xMin = (i == 0) ? 0.0 : -2.0; // Set minimum x to 0 for the first plot
        double xMax = (i == 0) ? 1.0 : 4.0; // Set maximum x to 1 for the first plot

        // Creating H2 histogram (red)
        hists[i] = createHistogram(trees[i], (std::string("h_rga_") + titles[i]).c_str(), titles[i], variables[i], "", rga_H2_norm);
        hists[i]->SetLineColor(kRed);
        hists[i]->GetXaxis()->SetTitle("M_{X} (GeV)");
        hists[i]->GetYaxis()->SetTitle("Counts / nC");
        hists[i]->GetXaxis()->SetTitleSize(0.05);
        hists[i]->GetYaxis()->SetTitleSize(0.05);
        hists[i]->Draw();

        // Creating NH3 histogram (blue)
        hists[i + 4] = createHistogram(trees[i + 4], (std::string("h_rgc_nh3_") + titles[i]).c_str(), "", variables[i], cuts_NH3[i], rgc_NH3_norm);
        hists[i + 4]->SetLineColor(kBlue);
        hists[i + 4]->Draw("same");

        // Creating C histogram (green)
        hists[i + 8] = createHistogram(trees[i + 4], (std::string("h_rgc_c_") + titles[i]).c_str(), "", variables[i], cuts_C[i], rgc_C_norm);
        hists[i + 8]->SetLineColor(kGreen);
        hists[i + 8]->Draw("same");

        // Determine the maximum value of the histograms
        double maxVal = TMath::Max(hists[i]->GetMaximum(), hists[i + 4]->GetMaximum());
        maxVal = TMath::Max(maxVal, hists[i + 8]->GetMaximum());
        double newMax = maxVal * 1.2;  // 20% higher than the maximum value
        hists[i]->SetMaximum(newMax);
        hists[i + 4]->SetMaximum(newMax);
        hists[i + 8]->SetMaximum(newMax);

        // Add legend in the top left
        TLegend* leg = new TLegend(0.15, 0.7, 0.35, 0.9); // Adjusted coordinates for top left
        leg->AddEntry(hists[i], "H2", "l");
        leg->AddEntry(hists[i + 4], "NH3", "l");
        leg->AddEntry(hists[i + 8], "C", "l");
        leg->Draw();

        // Set x-axis label for eX plots
        const char* xAxisTitle = (i == 0) ? "x_{B}" : "M_{X} (GeV)";
        hists[i]->GetXaxis()->SetTitle(xAxisTitle);
        hists[i + 4]->GetXaxis()->SetTitle(xAxisTitle);
        hists[i + 8]->GetXaxis()->SetTitle(xAxisTitle);

        // Right column plots (NH3/C ratio)
        c1->cd(i * 2 + 2);
        TPad* pad2 = (TPad*)c1->GetPad(i * 2 + 2);
        pad2->SetBottomMargin(0.15);
        pad2->SetLeftMargin(0.15);

        TH1D* ratioHist = (TH1D*)hists[i + 4]->Clone((std::string("ratio_") + titles[i]).c_str());
        ratioHist->Divide(hists[i + 8]);
        ratioHist->SetLineColor(kBlack);  // Changed to black
        ratioHist->GetXaxis()->SetTitle("M_{X} (GeV)");
        ratioHist->GetYaxis()->SetTitle("NH_{3} / C");
        ratioHist->GetXaxis()->SetTitleSize(0.05);
        ratioHist->GetYaxis()->SetTitleSize(0.05);
        ratioHist->Draw();

        // Add label for ratio plots
        TLatex *ratioLabel = new TLatex();
        ratioLabel->SetTextSize(0.04);
        ratioLabel->DrawLatexNDC(0.7, 0.8, "NH3/C Ratio");
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


