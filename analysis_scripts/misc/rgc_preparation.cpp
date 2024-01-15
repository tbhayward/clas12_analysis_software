#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <iostream>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TF1.h>         // For TF1 class
#include <TGraphErrors.h> // For TGraphErrors class
#include <TMath.h> // For TMath::Max

const std::string output_dir = "output/rgc_ready_for_cooking_plots/";

TH1D* createHistogram(TTree* tree, const char* name, const char* title, const char* variable, 
    const char* cut, double norm, double xMin, double xMax) {
    TH1D* hist = new TH1D(name, "", 30, xMin, xMax);
    tree->Draw((std::string(variable) + ">>" + name).c_str(), cut, "goff");
    hist->Scale(1.0 / norm);
    hist->SetStats(kFALSE);
    return hist;
}

void rgc_preparation() {
    const char* files[] = {
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_elastic_eX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epi+X.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epi+pi-X.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_elastic_eX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epi+X.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epi+pi-X.root"
    };

    const char* newFilesEX[] = {
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_eX.root",
        "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_eX.root"
    };

    TFile* filesOpened[8];
    TTree* trees[8];
    for (int i = 0; i < 8; i++) {
        filesOpened[i] = new TFile(files[i]);
        trees[i] = (TTree*)filesOpened[i]->Get("PhysicsEvents");
    }

    double rga_H2_norm = 159661.55+145813.73;
    // double rga_H2_norm = 53381.99+41401.77;
    // double rgc_pos_NH3_norm = 19355.9+19392.53+21683.25+21621.178;
    double rgc_pos_NH3_norm = 40926.816+45619.707;
    // double rgc_neg_NH3_norm = 21282.264+21217.414+21303.576+21297.766;
    double rgc_neg_NH3_norm = 44842.426+44852.715;
    double rgc_NH3_norm = rgc_pos_NH3_norm+rgc_neg_NH3_norm;
    double rgc_C_norm = 8883.014+8834.256;

    // // Compute normalization factors based on the number of entries under specific conditions
    // double rga_H2_norm = trees[0]->GetEntries();
    // double rgc_pos_NH3_norm = trees[4]->GetEntries("runnum == 16320 || runnum == 16327");
    // double rgc_neg_NH3_norm = trees[4]->GetEntries("runnum == 16346 || runnum == 16353");
    // double rgc_NH3_norm = rgc_pos_NH3_norm+rgc_neg_NH3_norm;
    // double rgc_C_norm = trees[4]->GetEntries("runnum == 16297");

    TCanvas *c1 = new TCanvas("c1", "Data Analysis", 2200, 1200);
    c1->Divide(3, 4);

    TH1D* hists[12]; // 4 plots * 3 histograms per plot
    TH1D* xHists[12]; // 4 plots * 3 histograms per plot for "x"
    double xBMin = 0;
    double xBMax = 0.7;

    gStyle->SetOptStat(0);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(.8);
    gStyle->SetLabelSize(0.08, "XY");
    gStyle->SetTitleSize(0.08, "XY");

    const char* titles[] = {"eX", "e#pi^{+}X", "epX", "e#pi^{+}#pi^{-}X"};
    const char* variables[] = {"Mx", "Mx", "Mx", "Mx"};
    const char* cuts_NH3[] = {"runnum != 16297", "runnum != 16297", "runnum != 16297", "runnum != 16297"};
    const char* cuts_C[] = {"runnum == 16297", "runnum == 16297", "runnum == 16297", "runnum == 16297"};

    for (int i = 0; i < 4; i++) {
        c1->cd(i * 3 + 1);
        TPad* pad1 = (TPad*)c1->GetPad(i * 3 + 1);
        pad1->SetBottomMargin(0.20);
        pad1->SetLeftMargin(0.20);


        // Define histogram ranges for each channel
        double xMin, xMax;
        if (i == 0) {        // eX
            xMin = 0.4; xMax = 1.3;
        } else if (i == 1) { // epi+X
            xMin = 0.4; xMax = 1.3;
        } else if (i == 2) { // epX
            xMin = -1.0; xMax = 1.0;
        } else {             // epi+pi-X
            xMin = 0.4; xMax = 1.3;
        }

        // Creating H2 histogram (red)
        hists[i] = createHistogram(trees[i], (std::string("h_rga_") + titles[i]).c_str(), titles[i], variables[i], "", rga_H2_norm, xMin, xMax);
        hists[i]->SetLineColor(kRed);
        hists[i]->GetXaxis()->SetTitle("M_{X} (GeV)");
        hists[i]->GetYaxis()->SetTitle("Counts / nC");
        hists[i]->GetXaxis()->SetTitleSize(0.08);
        hists[i]->GetYaxis()->SetTitleSize(0.08);
        hists[i]->Draw();

        // Creating NH3 histogram (blue)
        hists[i + 4] = createHistogram(trees[i + 4], (std::string("h_rgc_nh3_") + titles[i]).c_str(), "", variables[i], cuts_NH3[i], rgc_NH3_norm, xMin, xMax);
        hists[i + 4]->SetLineColor(kBlue);
        hists[i + 4]->Draw("same");

        // Creating C histogram (green)
        hists[i + 8] = createHistogram(trees[i + 4], (std::string("h_rgc_c_") + titles[i]).c_str(), "", variables[i], cuts_C[i], rgc_C_norm, xMin, xMax);
        hists[i + 8]->SetLineColor(kGreen);
        hists[i + 8]->Draw("same");

        // Create a title using TLatex at the top of each pad
        TLatex *title = new TLatex();
        title->SetTextSize(0.08); // Adjust text size as needed
        title->SetTextAlign(22); // Center alignment
        title->DrawLatexNDC(0.5, 0.98, titles[i]); // Draw title at the top center of the pad

        // Determine the maximum value of the histograms
        double maxVal = TMath::Max(hists[i]->GetMaximum(), hists[i + 4]->GetMaximum());
        maxVal = TMath::Max(maxVal, hists[i + 8]->GetMaximum());
        double newMax = maxVal * 1.2;  // 20% higher than the maximum value
        hists[i]->SetMaximum(newMax);
        hists[i + 4]->SetMaximum(newMax);
        hists[i + 8]->SetMaximum(newMax);

        // Add legend in the top left
        TLegend* leg = new TLegend(0.20, 0.7, 0.35, 0.9); // Adjusted coordinates for top left
        leg->AddEntry(hists[i], "H2 (RGA)", "l");
        leg->AddEntry(hists[i + 4], "NH3 (RGC)", "l");
        leg->AddEntry(hists[i + 8], "C (RGC)", "l");
        leg->Draw();

        /* ~~~~~~~~~~~~~~~~ SECOND COLUMN ~~~~~~~~~~~~~~~~ */

        // Right column plots (NH3/C ratio)
        c1->cd(i * 3 + 2);
        TPad* pad2 = (TPad*)c1->GetPad(i * 3 + 2);
        pad2->SetBottomMargin(0.20);
        pad2->SetLeftMargin(0.15);

        TH1D* ratioHist = (TH1D*)hists[i + 4]->Clone((std::string("ratio_") + titles[i]).c_str());
        ratioHist->Divide(hists[i + 8]);
        ratioHist->SetLineColor(kBlack);  // Changed to black
        // Set y-axis range from 0.5 to 1.2
        ratioHist->SetMinimum(0.5);
        ratioHist->SetMaximum(1.5);
        ratioHist->GetXaxis()->SetTitle("M_{X} (GeV)");
        ratioHist->GetYaxis()->SetTitle("NH_{3} / C");
        ratioHist->GetXaxis()->SetTitleSize(0.08);
        ratioHist->GetYaxis()->SetTitleSize(0.08);
        ratioHist->Draw();

        // Define fit range based on the channel
        double fitMin, fitMax;
        if (i == 0 || i == 1 || i == 3) { // eX, epi+X, epi+pi-X
            fitMin = 0.4; fitMax = 0.8;
        } else if (i == 2) { // epX
            fitMin = -1.0; fitMax = -0.2;
        }

        // Perform a linear fit on the histogram
        TF1 *fitFunc = new TF1("fitFunc", "pol0", fitMin, fitMax);
        ratioHist->Fit(fitFunc, "RQ"); // "R" for range, "Q" for quiet mode

        // Draw the fit function over the full range
        fitFunc->SetLineColor(kRed);
        fitFunc->SetRange(fitMax,5);
        fitFunc->SetLineStyle(3); // Dashed line
        fitFunc->Draw("same");

        // Draw the fit function again over the fit range with a solid line
        TF1 *fitFuncSolid = (TF1*)fitFunc->Clone();
        fitFuncSolid->SetRange(fitMin, fitMax);
        fitFuncSolid->SetLineStyle(1); // Solid line
        fitFuncSolid->Draw("same");

        // Add updated label for ratio plots
        double normalization = fitFunc->GetParameter(0); // Get the constant value
        double normalizationError = fitFunc->GetParError(0); // Get the error on the constant
        char label[100];
        sprintf(label, "s (normalization) = %.3f #pm %.3f", 
            normalization, normalizationError); // Format the label with uncertainty

        TLatex *ratioLabel = new TLatex();
        ratioLabel->SetTextSize(0.06);
        ratioLabel->DrawLatexNDC(0.6, 0.8, label);

        /* ~~~~~~~~~~~~~~~~ THIRD COLUMN ~~~~~~~~~~~~~~~~ */

        // Determine which file to use for the third column
        const char* fileH2 = (i == 0) ? newFilesEX[0] : files[i];
        const char* fileNH3 = (i == 0) ? newFilesEX[1] : files[i + 4];
        const char* fileC = (i == 0) ? newFilesEX[1] : files[i + 4];

        // Open the new files for H2 and NH3 histograms for the third column (eX case)
        TFile* fileH2Opened = new TFile(fileH2);
        TFile* fileNH3Opened = new TFile(fileNH3);
        TTree* treeH2 = (TTree*)fileH2Opened->Get("PhysicsEvents");
        TTree* treeNH3 = (TTree*)fileNH3Opened->Get("PhysicsEvents");

        // Right column plots (H2 and scaled difference)
        c1->cd(i * 3 + 3); // Adjust to access the third column
        TPad* pad3 = (TPad*)c1->GetPad(i * 3 + 3);
        pad3->SetBottomMargin(0.25);
        pad3->SetLeftMargin(0.2);

        // Create histograms for "x" for H2, NH3, and C
        xHists[i] = createHistogram(treeH2, 
            (std::string("x_hist_h2_") + std::to_string(i)).c_str(), "", "x", "", 
            rga_H2_norm, xBMin, xBMax);
        xHists[i + 4] = createHistogram(treeNH3, 
            (std::string("x_hist_nh3_") + std::to_string(i)).c_str(), "", "x", 
            cuts_NH3[i], rgc_NH3_norm, xBMin, xBMax);
        xHists[i + 8] = createHistogram(trees[i + 4], 
            (std::string("x_hist_c_") + std::to_string(i)).c_str(), "", "x", 
            cuts_C[i], rgc_C_norm, xBMin, xBMax);


        // Plot H2 histogram for "x" (red)
        TH1D* xH2Hist = (TH1D*)xHists[i]->Clone();
        xH2Hist->SetLineColor(kRed);
        xH2Hist->Draw();

        // Create and plot the difference histogram for "x" (black)
        TH1D* xDiffHist = (TH1D*)xHists[i + 4]->Clone(); 
        xDiffHist->Scale(1/normalization);
        TH1D* scaledCHist = (TH1D*)xHists[i + 8]->Clone(); 
        scaledCHist->Scale(-1);
        xDiffHist->Add(scaledCHist);
        xDiffHist->SetLineColor(kBlack);
        xDiffHist->Draw("same");

        // // Create and plot the difference histogram for "x" (black)
        // TH1D* xDiffHist = (TH1D*)xHists[i + 4]->Clone(); 
        // TH1D* scaledCHist = (TH1D*)xHists[i + 8]->Clone(); 
        // scaledCHist->Scale(-normalization);
        // xDiffHist->Add(scaledCHist);
        // xDiffHist->SetLineColor(kBlack);
        // xDiffHist->Draw("same");

        // Set y-axis title for xDiffHist
        xDiffHist->GetYaxis()->SetTitle("NH_{3} - s*C (Counts/nC)");
        xDiffHist->GetYaxis()->SetTitleSize(0.08);

        // Find the maximum value between xH2Hist and xDiffHist
        double maxValThirdCol = TMath::Max(xH2Hist->GetMaximum(), xDiffHist->GetMaximum());
        double newMaxThirdCol = maxValThirdCol * 1.20; // 20% higher than the max value
        xH2Hist->SetMaximum(newMaxThirdCol); // Set max for H2 histogram
        xDiffHist->SetMaximum(newMaxThirdCol); // Set max for difference histogram

        // Set axis titles
        xH2Hist->GetXaxis()->SetTitle("x_{B}");
        xH2Hist->GetXaxis()->SetTitleSize(0.08);
        xDiffHist->GetXaxis()->SetTitle("x_{B}");
        xDiffHist->GetXaxis()->SetTitleSize(0.08);

        // Draw histograms again to update axis settings
        xH2Hist->Draw();
        xDiffHist->Draw("same");

        // Add legend for the third column at top right
        TLegend* thirdColLeg = new TLegend(0.70, 0.75, 0.90, 0.90); // Top right coordinates
        thirdColLeg->AddEntry(xH2Hist, "H2 (RGA)", "l");
        thirdColLeg->AddEntry(xDiffHist, "NH_{3} - s*C (RGC)", "l");
        thirdColLeg->Draw();
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


