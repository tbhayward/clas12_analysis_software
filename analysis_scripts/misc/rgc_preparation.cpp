#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <iostream>
#include <string>
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

    double rga_H2_norm = 51.79738;
    // double rga_H2_norm = 53381.99+41401.77;
    // double rgc_pos_NH3_norm = 19355.9+19392.53+21683.25+21621.178;
    double rgc_pos_NH3_norm = 45.7595+41.050555;
    // double rgc_neg_NH3_norm = 21282.264+21217.414+21303.576+21297.766;
    double rgc_neg_NH3_norm = 44.970203+44.98406;
    double rgc_NH3_norm = rgc_pos_NH3_norm+rgc_neg_NH3_norm;
    double rgc_C_norm = 18.91757;

    // double rga_H2_norm = 159661.55;
    // double rgc_NH3_norm = 41392.934 + 43299.863;
    // double rgc_C_norm =  43098.254;

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

    gStyle->SetOptStat(0);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(.8);
    gStyle->SetLabelSize(0.08, "XY");
    gStyle->SetTitleSize(0.08, "XY");

    const char* titles[] = {"eX", "e#pi^{+}X", "epX", "e#pi^{+}#pi^{-}X"};
    const char* variables[] = {"Mx", "Mx", "Mx", "Mx"};
    const char* cuts_NH3[] = {"runnum != 16297", "runnum != 16297", "runnum != 16297", "runnum != 16297"};
    const char* cuts_C[] = {"runnum == 16297", "runnum == 16297", "runnum == 16297", "runnum == 16297"};

    // const char* cuts_NH3[] = {"runnum != 16293", "runnum != 16293", "runnum != 16293", "runnum != 16293"};
    // const char* cuts_C[] = {"runnum == 16293", "runnum == 16293", "runnum == 16293", "runnum == 16293"};


    for (int i = 0; i < 4; i++) {
        c1->cd(i * 3 + 1);
        TPad* pad1 = (TPad*)c1->GetPad(i * 3 + 1);
        pad1->SetBottomMargin(0.20);
        pad1->SetLeftMargin(0.20);


        // Define histogram ranges for each channel
        double xMin, xMax;
        if (i == 0) {        // eX
            xMin = 0.3; xMax = 1.1;
        } else if (i == 1) { // epi+X
            xMin = 0.3; xMax = 1.1;
        } else if (i == 2) { // epX
            xMin = -1.0; xMax = 1.0;
        } else {             // epi+pi-X
            xMin = 0.3; xMax = 1.1;
        }

        // Creating H2 histogram (red)
        hists[i] = createHistogram(trees[i], (std::string("h_rga_") + titles[i]).c_str(), titles[i], variables[i], "", rga_H2_norm, xMin, xMax);
        hists[i]->SetLineColor(kRed);
        hists[i]->GetXaxis()->SetTitle("M_{X} (GeV)");
        hists[i]->GetYaxis()->SetTitle("Counts / mC");
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
            fitMin = 0.3; fitMax = 0.75;
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
        c1->cd(i * 3 + 3); // Access the third column for each row
        TPad* pad3 = (TPad*)c1->GetPad(i * 3 + 3);
        pad3->SetBottomMargin(0.20);
        pad3->SetLeftMargin(0.20);

        // Modify cuts with additional constraints
        const char* additionalCuts = "W > 0 && y < 0.75";

        // For NH3 histograms
        std::string nh3CutsWithConstraints = std::string(cuts_NH3[i]) + " && " + additionalCuts;

        // For Carbon histograms
        std::string cCutsWithConstraints = std::string(cuts_C[i]) + " && " + additionalCuts;

        // For H2 histograms, as there are no initial cuts, just use the additional constraints
        std::string h2CutsWithConstraints = additionalCuts;

        std::string var = "x";
        if (i == 0) {        // eX
            xMin = 0.0; xMax = 10; var = "Q2";
        } else if (i == 1) { // epi+X
            xMin = 0.0; xMax = 1.2; var = "pT";
        } else if (i == 2) { // epX
            xMin = -1.0; xMax = 1.0; var = "xF";
        } else {             // epi+pi-X
            xMin = 0.0; xMax = 1.5; var = "Mh";
        }

        // Creating a new NH3 histogram for Mx in the third column
        TH1D* nh3HistThirdCol = createHistogram(trees[i + 4],
            (std::string("nh3_third_col_") + titles[i]).c_str(), titles[i],
            var.c_str(), nh3CutsWithConstraints.c_str(), rgc_NH3_norm, xMin, xMax);

        // Creating a new Carbon histogram for Mx in the third column
        TH1D* cHistThirdCol = createHistogram(trees[i + 4], 
            (std::string("c_third_col_") + titles[i]).c_str(), titles[i],
        var.c_str(), cCutsWithConstraints.c_str(), rgc_C_norm, xMin, xMax);

        // Scale the Carbon histogram by the normalization factor
        cHistThirdCol->Scale(normalization);

        // Subtract the scaled Carbon histogram from the NH3 histogram
        nh3HistThirdCol->Add(cHistThirdCol, -1); // The second argument "-1" is for subtraction

        // Set line color for the resulting histogram
        nh3HistThirdCol->SetLineColor(kBlack);

        // Set titles and sizes for axes
        nh3HistThirdCol->GetXaxis()->SetTitle("x_{B}");
        nh3HistThirdCol->GetYaxis()->SetTitle("NH_{3} - s*C (Counts/mC)");
        nh3HistThirdCol->GetXaxis()->SetTitleSize(0.08);
        nh3HistThirdCol->GetYaxis()->SetTitleSize(0.08);

        // integrate histogram
        nh3HistThirdCol->Scale(1/nh3HistThirdCol->Integral());

        // Draw the resulting histogram
        nh3HistThirdCol->Draw();

        // Creating a new H2 histogram for Mx in the third column
        TH1D* h2HistThirdCol = createHistogram(trees[i], 
            (std::string("h2_third_col_") + titles[i]).c_str(), titles[i], 
            var.c_str(), h2CutsWithConstraints.c_str(), rga_H2_norm, xMin, xMax);
        h2HistThirdCol->SetLineColor(kRed); // Set line color to red
        h2HistThirdCol->GetXaxis()->SetTitle("x_{B}");
        h2HistThirdCol->GetYaxis()->SetTitle("Counts / mC");
        h2HistThirdCol->GetXaxis()->SetTitleSize(0.08);
        h2HistThirdCol->GetYaxis()->SetTitleSize(0.08);
        
        // integrate histogram
        h2HistThirdCol->Scale(1/h2HistThirdCol->Integral());

        // Overlay the H2 histogram on top of the NH3 minus Carbon histogram
        h2HistThirdCol->Draw("same");



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


