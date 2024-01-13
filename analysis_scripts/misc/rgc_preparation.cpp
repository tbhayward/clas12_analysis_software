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

TH1D* createHistogram(TTree* tree, const char* name, const char* title, const char* variable, 
    const char* cut, double norm, double xMin, double xMax) {
    TH1D* hist = new TH1D(name, title, 100, xMin, xMax);
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

        // Define histogram ranges for each channel
        double xMin, xMax;
        if (i == 0) {        // eX
            xMin = 1.0; xMax = 4.0;
        } else if (i == 1) { // epi+X
            xMin = 0.0; xMax = 3.5;
        } else if (i == 2) { // epX
            xMin = -2.0; xMax = 3.0;
        } else {             // epi+pi-X
            xMin = 0.0; xMax = 3.0;
        }

        // Creating H2 histogram (red)
        hists[i] = createHistogram(trees[i], (std::string("h_rga_") + titles[i]).c_str(), titles[i], variables[i], "", rga_H2_norm, xMin, xMax);
        hists[i]->SetLineColor(kRed);
        hists[i]->GetXaxis()->SetTitle("M_{X} (GeV)");
        hists[i]->GetYaxis()->SetTitle("Counts / nC");
        hists[i]->GetXaxis()->SetTitleSize(0.05);
        hists[i]->GetYaxis()->SetTitleSize(0.05);
        hists[i]->Draw();

        // Creating NH3 histogram (blue)
        hists[i + 4] = createHistogram(trees[i + 4], (std::string("h_rgc_nh3_") + titles[i]).c_str(), "", variables[i], cuts_NH3[i], rgc_NH3_norm, xMin, xMax);
        hists[i + 4]->SetLineColor(kBlue);
        hists[i + 4]->Draw("same");

        // Creating C histogram (green)
        hists[i + 8] = createHistogram(trees[i + 4], (std::string("h_rgc_c_") + titles[i]).c_str(), "", variables[i], cuts_C[i], rgc_C_norm, xMin, xMax);
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

        // Right column plots (NH3/C ratio)
        c1->cd(i * 2 + 2);
        TPad* pad2 = (TPad*)c1->GetPad(i * 2 + 2);
        pad2->SetBottomMargin(0.15);
        pad2->SetLeftMargin(0.15);

        TH1D* ratioHist = (TH1D*)hists[i + 4]->Clone((std::string("ratio_") + titles[i]).c_str());
        ratioHist->Divide(hists[i + 8]);
        ratioHist->SetLineColor(kBlack);  // Changed to black
        // Set y-axis range from 0.5 to 1.2
        ratioHist->SetMinimum(0.5);
        ratioHist->SetMaximum(1.2);
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


