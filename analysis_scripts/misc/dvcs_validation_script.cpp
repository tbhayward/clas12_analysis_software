#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
#include <TGraphErrors.h>
#include <TLine.h>

void plot_dvcs_energy_loss_validation(const char* file1, const char* file2, const char* titleSuffix) {
    // Open the ROOT files
    TFile *f1 = new TFile(file1);
    TFile *f2 = new TFile(file2);

    // Access the "PhysicsEvents" trees
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    // Prepare variables for reading the branches
    Double_t p1_theta1, p1_theta2;
    Double_t Mx2_1_1, Mx2_1_2;
    Double_t eta2_1, eta2_2;
    Double_t t1_1, t1_2;
    Double_t theta_gamma_gamma_1, theta_gamma_gamma_2;
    Double_t Emiss2_1, Emiss2_2;
    Double_t pTmiss_1, pTmiss_2;

    // Set branch addresses
    tree1->SetBranchAddress("p1_theta", &p1_theta1);
    tree1->SetBranchAddress("Mx2_1", &Mx2_1_1);
    tree1->SetBranchAddress("eta2", &eta2_1);
    tree1->SetBranchAddress("t1", &t1_1);
    tree1->SetBranchAddress("theta_gamma_gamma", &theta_gamma_gamma_1);
    tree1->SetBranchAddress("Emiss2", &Emiss2_1);
    tree1->SetBranchAddress("pTmiss", &pTmiss_1);

    tree2->SetBranchAddress("p1_theta", &p1_theta2);
    tree2->SetBranchAddress("Mx2_1", &Mx2_1_2);
    tree2->SetBranchAddress("eta2", &eta2_2);
    tree2->SetBranchAddress("t1", &t1_2);
    tree2->SetBranchAddress("theta_gamma_gamma", &theta_gamma_gamma_2);
    tree2->SetBranchAddress("Emiss2", &Emiss2_2);
    tree2->SetBranchAddress("pTmiss", &pTmiss_2);

    // Create canvas and divide it into 3x4 subplots
    TCanvas *c1 = new TCanvas("c1", "DVCS Energy Loss Validation", 1200, 900);
    c1->Divide(4, 3);

    // Set the theta bins and corresponding histogram ranges
    const int nBins = 10; // 10 theta bins + 1 fully integrated case
    Double_t thetaBins[nBins + 1] = {5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65}; // 10 equally spaced bins

    TH1D *h1[nBins + 1]; // +1 for the fully integrated case
    TH1D *h2[nBins + 1]; // +1 for the fully integrated case

    Double_t mu1_values[nBins], sigma1_values[nBins];
    Double_t mu2_values[nBins], sigma2_values[nBins];
    Double_t theta_sum[nBins] = {0.0};
    Int_t theta_count[nBins] = {0};
    Double_t theta_mean[nBins] = {0.0};

    // Create histograms for each theta bin with 50 bins
    h1[0] = new TH1D("h1_integrated", Form("Integrated #theta [5, 65] %s", titleSuffix), 35, -0.3, 0.3);
    h2[0] = new TH1D("h2_integrated", Form("Integrated #theta [5, 65] %s", titleSuffix), 35, -0.3, 0.3);

    for (int i = 0; i < nBins; ++i) {
        h1[i + 1] = new TH1D(Form("h1_%d", i), Form("#theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 35, -0.3, 0.3);
        h2[i + 1] = new TH1D(Form("h2_%d", i), Form("#theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 35, -0.3, 0.3);
    }

    // Fill the histograms with the cuts applied and calculate theta means
    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t j = 0; j < nEntries1; ++j) {
        tree1->GetEntry(j);
        Double_t thetaDeg1 = p1_theta1 * (180.0 / TMath::Pi()); // Convert to degrees
        if (thetaDeg1 >= 5 && thetaDeg1 < 65 
            && eta2_1 < 0 
            && t1_1 > -2 
            && theta_gamma_gamma_1 < 0.6
            && Emiss2_1 < 0.5 
            && pTmiss_1 < 0.125
                ) {
            h1[0]->Fill(Mx2_1_1); 
            // h1[0]->Fill(theta_gamma_gamma_1); 
        }
        for (int i = 0; i < nBins; ++i) {
            if (thetaDeg1 >= thetaBins[i] && thetaDeg1 < thetaBins[i + 1] 
                && eta2_1 < 0 
                && t1_1 > -2 
                && theta_gamma_gamma_1 < 0.6
                && Emiss2_1 < 0.5 
                && pTmiss_1 < 0.125
                ) {
                h1[i + 1]->Fill(Mx2_1_1);
                // h1[i + 1]->Fill(theta_gamma_gamma_1);
                theta_sum[i] += thetaDeg1;
                theta_count[i]++;
            }
        }
    }

    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t j = 0; j < nEntries2; ++j) {
        tree2->GetEntry(j);
        Double_t thetaDeg2 = p1_theta2 * (180.0 / TMath::Pi()); // Convert to degrees
        if (thetaDeg2 >= 5 && thetaDeg2 < 65 
            && eta2_2 < 0 
            && t1_2 > -2 
            && theta_gamma_gamma_2 < 0.6
            && Emiss2_2 < 0.5 
            && pTmiss_2 < 0.125
                ) {
            h2[0]->Fill(Mx2_1_2); // Fully integrated case
            // h2[0]->Fill(theta_gamma_gamma_2); // Fully integrated case
        }
        for (int i = 0; i < nBins; ++i) {
            if (thetaDeg2 >= thetaBins[i] && thetaDeg2 < thetaBins[i + 1] 
                && eta2_2 < 0 
                && t1_2 > -2 
                && theta_gamma_gamma_2 < 0.6 
                && Emiss2_2 < 0.5 
                && pTmiss_2 < 0.125
                ) {
                h2[i + 1]->Fill(Mx2_1_2);
                // h2[i + 1]->Fill(theta_gamma_gamma_2);
                theta_sum[i] += thetaDeg2;
                theta_count[i]++;
            }
        }
    }

    // Calculate mean theta values
    for (int i = 0; i < nBins; ++i) {
        if (theta_count[i] > 0) {
            theta_mean[i] = theta_sum[i] / theta_count[i];
        } else {
            theta_mean[i] = 0.5 * (thetaBins[i] + thetaBins[i + 1]); // Default to midpoint if no entries
        }
    }

    // Draw the fully integrated case in the first subplot
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);

    Double_t maxVal1 = h1[0]->GetMaximum();
    Double_t maxVal2 = h2[0]->GetMaximum();
    Double_t maxVal = std::max(maxVal1, maxVal2);
    h1[0]->SetMaximum(1.7 * maxVal);
    h1[0]->SetMinimum(0);
    h2[0]->SetMaximum(1.7 * maxVal);
    h2[0]->SetMinimum(0);

    h1[0]->SetMarkerStyle(20);
    h1[0]->SetMarkerSize(0.8);
    h1[0]->SetMarkerColor(kBlack);
    h1[0]->SetStats(0);
    h1[0]->Draw("E");

    h2[0]->SetMarkerStyle(21);
    h2[0]->SetMarkerSize(0.8);
    h2[0]->SetMarkerColor(kRed);
    h2[0]->SetStats(0);
    h2[0]->Draw("E SAME");

    TF1 *fit1_int = new TF1("fit1_integrated", "gaus(0) + pol1(3)", -0.3, 0.3);
    TF1 *fit2_int = new TF1("fit2_integrated", "gaus(0) + pol1(3)", -0.3, 0.3);

    fit1_int->SetParameters(0.8 * maxVal1, 0, 0.2);
    fit1_int->SetParLimits(1, -0.15, 0.15);
    fit1_int->SetParLimits(2, 0, 0.3);
    fit2_int->SetParameters(0.8 * maxVal2, 0, 0.2);
    fit2_int->SetParLimits(1, -0.15, 0.15);
    fit2_int->SetParLimits(2, 0, 0.3);

    fit1_int->SetLineWidth(1);
    fit2_int->SetLineWidth(1);
    h1[0]->Fit(fit1_int, "Q");
    h2[0]->Fit(fit2_int, "Q");

    fit1_int->SetLineColor(kBlack);
    fit1_int->Draw("SAME");
    fit2_int->SetLineColor(kRed);
    fit2_int->Draw("SAME");

    TLegend *legend_int = new TLegend(0.25, 0.75, 0.9, 0.9);
    legend_int->SetTextSize(0.03);
    legend_int->AddEntry(h1[0], Form("Uncorrected: #mu=%.3f, #sigma=%.3f", fit1_int->GetParameter(1), fit1_int->GetParameter(2)), "lep");
    legend_int->AddEntry(h2[0], Form("Corrected: #mu=%.3f, #sigma=%.3f", fit2_int->GetParameter(1), fit2_int->GetParameter(2)), "lep");
    legend_int->Draw();

    h1[0]->GetXaxis()->SetTitle("M_{xp}^{2} (GeV^{2})");
    h1[0]->GetYaxis()->SetTitle("Counts");

    // Draw the remaining subplots for the binned cases
    for (int i = 1; i <= nBins; ++i) {
        c1->cd(i + 1)->SetLeftMargin(0.15);
        c1->cd(i + 1)->SetBottomMargin(0.15);

        maxVal1 = h1[i]->GetMaximum();
        maxVal2 = h2[i]->GetMaximum();
        maxVal = std::max(maxVal1, maxVal2);
        h1[i]->SetMaximum(1.7 * maxVal);
        h1[i]->SetMinimum(0);
        h2[i]->SetMaximum(1.7 * maxVal);
        h2[i]->SetMinimum(0);

        h1[i]->SetMarkerStyle(20);
        h1[i]->SetMarkerSize(0.8);
        h1[i]->SetMarkerColor(kBlack);
        h1[i]->SetStats(0);
        h1[i]->Draw("E");

        h2[i]->SetMarkerStyle(21);
        h2[i]->SetMarkerSize(0.8);
        h2[i]->SetMarkerColor(kRed);
        h2[i]->SetStats(0);
        h2[i]->Draw("E SAME");

        TF1 *fit1 = new TF1(Form("fit1_%d", i), "gaus(0) + pol1(3)", -0.3, 0.3);
        TF1 *fit2 = new TF1(Form("fit2_%d", i), "gaus(0) + pol1(3)", -0.3, 0.3);

        fit1->SetParameters(0.8 * maxVal1, 0, 0.2);
        fit1->SetParLimits(1, -0.15, 0.15);
        fit1->SetParLimits(2, 0, 0.3);
        fit2->SetParameters(0.8 * maxVal2, 0, 0.2);
        fit2->SetParLimits(1, -0.15, 0.15);
        fit2->SetParLimits(2, 0, 0.3);

        fit1->SetLineWidth(1);
        fit2->SetLineWidth(1);
        h1[i]->Fit(fit1, "Q");
        h2[i]->Fit(fit2, "Q");

        fit1->SetLineColor(kBlack);
        fit1->Draw("SAME");
        fit2->SetLineColor(kRed);
        fit2->Draw("SAME");

        mu1_values[i - 1] = fit1->GetParameter(1);
        sigma1_values[i - 1] = fit1->GetParameter(2);
        mu2_values[i - 1] = fit2->GetParameter(1);
        sigma2_values[i - 1] = fit2->GetParameter(2);

        TLegend *legend = new TLegend(0.25, 0.75, 0.9, 0.9); // Adjusted the legend position to the top right corner
        legend->SetTextSize(0.03); // Decrease the font size in the legend
        legend->AddEntry(h1[i], Form("Uncorrected: #mu=%.3f, #sigma=%.3f", mu1_values[i - 1], sigma1_values[i - 1]), "lep");
        legend->AddEntry(h2[i], Form("Corrected: #mu=%.3f, #sigma=%.3f", mu2_values[i - 1], sigma2_values[i - 1]), "lep");
        legend->Draw();

        h1[i]->GetXaxis()->SetTitle("M_{xp}^{2} (GeV^{2})");
        h1[i]->GetYaxis()->SetTitle("Counts");
    }

    // Plot TGraphErrors in the last subplot (12th pad)
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);

    TGraphErrors *gr1 = new TGraphErrors(nBins, theta_mean, mu1_values, 0, sigma1_values);
    TGraphErrors *gr2 = new TGraphErrors(nBins, theta_mean, mu2_values, 0, sigma2_values);

    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.8);
    gr1->SetMarkerColor(kBlack);
    gr1->Draw("AP");

    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(0.8);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("P SAME");

    // Add dashed gray line at y = 0
    TLine *line = new TLine(0, 0, 70, 0);
    line->SetLineColor(kGray);
    line->SetLineStyle(2); // Dashed line
    line->Draw("SAME");

    // Customize the last subplot
    gr1->SetTitle(""); // Remove the "Graph" title
    gr1->GetXaxis()->SetTitle("#theta"); // Set x-axis label
    gr1->GetYaxis()->SetTitle("#mu (GeV^{2})"); // Set y-axis label
    gr1->GetXaxis()->SetTitleSize(0.05); // Match font size with other plots
    gr1->GetYaxis()->SetTitleSize(0.05); // Match font size with other plots
    gr1->GetXaxis()->SetLimits(0, 70); // Set x-axis range
    gr1->GetYaxis()->SetRangeUser(-0.20, 0.20); // Set y-axis range

    // Add legend to the last plot, positioned in the top right
    TLegend *legend12 = new TLegend(0.6, 0.75, 0.9, 0.9); // Adjusted for horizontal and vertical size
    legend12->SetTextSize(0.03);
    legend12->AddEntry(gr1, "Uncorrected", "lep");
    legend12->AddEntry(gr2, "Corrected", "lep");
    legend12->Draw();

    // Save the canvas as a PDF
    TString outputFileName = "output/dvcs_" + TString(titleSuffix) + "_energy_loss_validation.pdf";
    c1->SaveAs(outputFileName);

    // Clean up
    delete c1;
    delete f1;
    delete f2;
} 

void plot_rho0_energy_loss_validation(const char* file1, const char* file2, const char* titleSuffix) {
    // Open the ROOT files
    TFile *f1 = new TFile(file1);
    TFile *f2 = new TFile(file2);

    // Access the "PhysicsEvents" trees
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    // Prepare variables for reading the branches
    Double_t p1_theta1, p1_theta2;
    Double_t Mx1_1, Mx1_2;
    Double_t Mx2_1, Mx2_2;

    // Set branch addresses
    tree1->SetBranchAddress("p_theta", &p1_theta1);
    tree1->SetBranchAddress("Mx", &Mx1_1);

    tree2->SetBranchAddress("p_theta", &p1_theta2);
    tree2->SetBranchAddress("Mx", &Mx1_2);

    // Create canvas and divide it into 3x4 subplots
    TCanvas *c1 = new TCanvas("c1", "Rho0 Energy Loss Validation", 1200, 900);
    c1->Divide(4, 3);

    // Set the theta bins and corresponding histogram ranges
    const int nBins = 11;
    Double_t thetaBins[nBins + 1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 65}; // 5 to 65 degrees

    TH1D *h1[nBins];
    TH1D *h2[nBins];

    Double_t mu1_values[nBins], sigma1_values[nBins];
    Double_t mu2_values[nBins], sigma2_values[nBins];
    Double_t theta_mean[nBins];

    // Rho0 mass in GeV/c^2
    const Double_t rho0_mass = 0.77526;

    for (int i = 0; i < nBins; ++i) {
        TPad *pad = (TPad*)c1->cd(i + 1);
        pad->SetLeftMargin(0.15); // Add padding to the left of each subplot
        pad->SetBottomMargin(0.15); // Add padding to the bottom of each subplot

        // Create histograms for each theta bin with 50 bins
        h1[i] = new TH1D(Form("h1_%d", i), Form("Mx1 for #theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 30, 0.55, 0.95);
        h2[i] = new TH1D(Form("h2_%d", i), Form("Mx1 for #theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 30, 0.55, 0.95);

        // Set text sizes
        h1[i]->GetXaxis()->SetTitleSize(0.05);
        h1[i]->GetYaxis()->SetTitleSize(0.05);
        h1[i]->GetXaxis()->SetLabelSize(0.04);
        h1[i]->GetYaxis()->SetLabelSize(0.04);

        h2[i]->GetXaxis()->SetTitleSize(0.05);
        h2[i]->GetYaxis()->SetTitleSize(0.05);
        h2[i]->GetXaxis()->SetLabelSize(0.04);
        h2[i]->GetYaxis()->SetLabelSize(0.04);

        // Fill the histograms without applying additional cuts
        Long64_t nEntries1 = tree1->GetEntries();
        // for (Long64_t j = 0; j < nEntries1; ++j) {
        for (Long64_t j = 0; j < 5000000; ++j) {
            tree1->GetEntry(j);
            Double_t thetaDeg1 = p1_theta1 * (180.0 / TMath::Pi()); // Convert to degrees
            if (thetaDeg1 >= thetaBins[i] && thetaDeg1 < thetaBins[i + 1]) {
                h1[i]->Fill(Mx1_1);
            }
        }

        Long64_t nEntries2 = tree2->GetEntries();
        // for (Long64_t j = 0; j < nEntries2; ++j) {
        for (Long64_t j = 0; j < 5000000; ++j) {    
            tree2->GetEntry(j);
            Double_t thetaDeg2 = p1_theta2 * (180.0 / TMath::Pi()); // Convert to degrees
            if (thetaDeg2 >= thetaBins[i] && thetaDeg2 < thetaBins[i + 1]) {
                h2[i]->Fill(Mx1_2);
            }
        }

        // Set the y-axis range to 0 to 1.7 times the maximum bin value
        Double_t maxVal1 = h1[i]->GetMaximum();
        Double_t maxVal2 = h2[i]->GetMaximum();
        Double_t maxVal = std::max(maxVal1, maxVal2);
        h1[i]->SetMaximum(1.7 * maxVal);
        h1[i]->SetMinimum(0);
        h2[i]->SetMaximum(1.7 * maxVal);
        h2[i]->SetMinimum(0);

        // Draw histograms as points with error bars
        h1[i]->SetMarkerStyle(20);
        h1[i]->SetMarkerSize(0.8); // Make the points smaller
        h1[i]->SetMarkerColor(kBlack);
        h1[i]->SetStats(0); // Remove stat box
        h1[i]->Draw("E");
        h2[i]->SetMarkerStyle(21);
        h2[i]->SetMarkerSize(0.8); // Make the points smaller
        h2[i]->SetMarkerColor(kRed);
        h2[i]->SetStats(0); // Remove stat box
        h2[i]->Draw("E SAME");

        // Fit histograms to Gaussian plus quadratic background
        TF1 *fit1 = new TF1(Form("fit1_%d", i), "gaus(0) + pol3(3)", 0.5, 1.3);
        TF1 *fit2 = new TF1(Form("fit2_%d", i), "gaus(0) + pol3 (3)", 0.5, 1.3);

        // Set initial parameter guesses and limits
        Double_t amplitudeGuess1 = 0.8 * maxVal1;
        Double_t amplitudeGuess2 = 0.8 * maxVal2;
        fit1->SetParameters(amplitudeGuess1, rho0_mass, 0.1);
        // fit1->SetParLimits(1, rho0_mass - 0.15, rho0_mass + 0.15); // mu limits around rho0 mass
        // fit1->SetParLimits(2, 0, 0.3);      // sigma limit
        fit2->SetParameters(amplitudeGuess2, rho0_mass, 0.1);
        // fit2->SetParLimits(1, rho0_mass - 0.05, rho0_mass + 0.05); // mu limits around rho0 mass
        // fit2->SetParLimits(2, 0, 0.3);      // sigma l imit

        fit1->SetLineWidth(1); // Make the line thinner
        fit2->SetLineWidth(1); // Make the line thinner
        h1[i]->Fit(fit1, "Q");
        h2[i]->Fit(fit2, "Q");

        // Draw fitted functions
        fit1->SetLineColor(kBlack);
        fit1->Draw("SAME");
        fit2->SetLineColor(kRed);
        fit2->Draw("SAME");

        // Get fit parameters (only from the Gaussian)
        mu1_values[i] = fit1->GetParameter(1);
        sigma1_values[i] = fit1->GetParameter(2);
        mu2_values[i] = fit2->GetParameter(1);
        sigma2_values[i] = fit2->GetParameter(2);

        // Calculate the mean theta value for the bin
        theta_mean[i] = 0.5 * (thetaBins[i] + thetaBins[i + 1]);

        // Add legend with mu and sigma values in the top right corner
        TLegend *legend = new TLegend(0.25, 0.75, 0.9, 0.9); // Adjusted the legend position to the top right corner
        legend->SetTextSize(0.03); // Decrease the font size in the legend
        legend->AddEntry(h1[i], Form("Uncorrected: #mu=%.3f, #sigma=%.3f", mu1_values[i], sigma1_values[i]), "lep");
        legend->AddEntry(h2[i], Form("Corrected: #mu=%.3f, #sigma=%.3f", mu2_values[i], sigma2_values[i]), "lep");
        legend->Draw();

        // Label the axes
        h1[i]->GetXaxis()->SetTitle("Mx1 GeV/c^2");
        h1[i]->GetYaxis()->SetTitle("Counts");
    }

    // Plot TGraphErrors in the last subplot (12th pad)
    TPad *pad12 = (TPad*) c1->cd(12);
    pad12->SetLeftMargin(0.15);
    pad12->SetBottomMargin(0.15);

    TGraphErrors *gr1 = new TGraphErrors(nBins, theta_mean, mu1_values, 0, sigma1_values);
    TGraphErrors *gr2 = new TGraphErrors(nBins, theta_mean, mu2_values, 0, sigma2_values);

    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.8);
    gr1->SetMarkerColor(kBlack);
    gr1->Draw("AP");

    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(0.8);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("P SAME");

    // Add dashed gray line at y = 0
    TLine *line = new TLine(thetaBins[0], rho0_mass, thetaBins[nBins], rho0_mass);
    line->SetLineColor(kGray);
    line->SetLineStyle(2); // Dashed line
    line->Draw("SAME");

    // Customize the last subplot
    gr1->SetTitle(""); // Remove the "Graph" title
    gr1->GetXaxis()->SetTitle("#theta"); // Set x-axis label
    gr1->GetYaxis()->SetTitle("#mu");    // Set y-axis label

    // Add legend to the last plot, positioned in the top right
    TLegend *legend12 = new TLegend(0.6, 0.75, 0.9, 0.9); // Adjusted for horizontal and vertical size
    legend12->SetTextSize(0.03);
    legend12->AddEntry(gr1, "Uncorrected", "lep");
    legend12->AddEntry(gr2, "Corrected", "lep");
    legend12->Draw();

    // Save the canvas as a PDF
    c1->SaveAs("output/rho0_energy_loss_validation.pdf");

    // Clean up
    delete c1;
    delete f1;
    delete f2;
}

void plot_elastic_energy_loss_validation(const char* file1, const char* file2, const char* titleSuffix) {
    // Open the ROOT files
    TFile *f1 = new TFile(file1);
    TFile *f2 = new TFile(file2);

    // Access the "PhysicsEvents" trees
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    // Prepare variables for reading the branches
    Double_t p1_theta1, p1_theta2;
    Double_t Mx2_1, Mx2_2;
    Double_t W_1, W_2;

    // Set branch addresses
    tree1->SetBranchAddress("p_theta", &p1_theta1);
    tree1->SetBranchAddress("Mx2", &Mx2_1);
    tree1->SetBranchAddress("W", &W_1);

    tree2->SetBranchAddress("p_theta", &p1_theta2);
    tree2->SetBranchAddress("Mx2", &Mx2_2);
    tree2->SetBranchAddress("W", &W_2);

    // Create canvas and divide it into 3x4 subplots
    TCanvas *c1 = new TCanvas("c1", "Elastic Energy Loss Validation", 1200, 900);
    c1->Divide(4, 3);

    // Set the theta bins and corresponding histogram ranges
    const int nBins = 11; // 10 theta bins + 1 fully integrated case
    // Double_t thetaBins[nBins + 1] = {5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65}; // 10 equally spaced bins
    Double_t thetaBins[nBins + 1] = {5, 31, 33, 35, 37, 39, 41, 43, 46, 49, 65}; // 10 equally spaced bins

    TH1D *h1[nBins + 1]; // +1 for the fully integrated case
    TH1D *h2[nBins + 1]; // +1 for the fully integrated case

    Double_t mu1_values[nBins], sigma1_values[nBins];
    Double_t mu2_values[nBins], sigma2_values[nBins];
    Double_t theta_sum[nBins] = {0.0};
    Int_t theta_count[nBins] = {0};
    Double_t theta_mean[nBins] = {0.0};

    // Create histograms for each theta bin with 50 bins
    h1[0] = new TH1D("h1_integrated", Form("Integrated #theta [5, 65] %s", titleSuffix), 25, -0.005, 0.005);
    h2[0] = new TH1D("h2_integrated", Form("Integrated #theta [5, 65] %s", titleSuffix), 25, -0.005, 0.005);

    for (int i = 0; i < nBins; ++i) {
        h1[i + 1] = new TH1D(Form("h1_%d", i), Form("#theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 25, -0.005, 0.005);
        h2[i + 1] = new TH1D(Form("h2_%d", i), Form("#theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 25, -0.005, 0.005);
    }

    // Fill the histograms and calculate theta means
    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t j = 0; j < nEntries1; ++j) {
        tree1->GetEntry(j);
        Double_t thetaDeg1 = p1_theta1 * (180.0 / TMath::Pi()); // Convert to degrees
        if (thetaDeg1 >= 5 && thetaDeg1 < 65 && W_1 < 1.0) {
            h1[0]->Fill(Mx2_1); // Fully integrated case
        }
        for (int i = 0; i < nBins; ++i) {
            if (thetaDeg1 >= thetaBins[i] && thetaDeg1 < thetaBins[i + 1] && W_1 < 1.0) {
                h1[i + 1]->Fill(Mx2_1);
                theta_sum[i] += thetaDeg1;
                theta_count[i]++;
            }
        }
    }

    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t j = 0; j < nEntries2; ++j) {
        tree2->GetEntry(j);
        Double_t thetaDeg2 = p1_theta2 * (180.0 / TMath::Pi()); // Convert to degrees
        if (thetaDeg2 >= 5 && thetaDeg2 < 65 && W_2 < 1.0) {
            h2[0]->Fill(Mx2_2); // Fully integrated case
        }
        for (int i = 0; i < nBins; ++i) {
            if (thetaDeg2 >= thetaBins[i] && thetaDeg2 < thetaBins[i + 1] && W_2 < 1.0) {
                h2[i + 1]->Fill(Mx2_2);
                theta_sum[i] += thetaDeg2;
                theta_count[i]++;
            }
        }
    }

    // Calculate mean theta values
    for (int i = 0; i < nBins; ++i) {
        if (theta_count[i] > 0) {
            theta_mean[i] = theta_sum[i] / theta_count[i];
        } else {
            theta_mean[i] = 0.5 * (thetaBins[i] + thetaBins[i + 1]); // Default to midpoint if no entries
        }
    }

    // Draw the fully integrated case in the first subplot
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);

    Double_t maxVal1 = h1[0]->GetMaximum();
    Double_t maxVal2 = h2[0]->GetMaximum();
    Double_t maxVal = std::max(maxVal1, maxVal2);
    h1[0]->SetMaximum(1.7 * maxVal);
    h1[0]->SetMinimum(0);
    h2[0]->SetMaximum(1.7 * maxVal);
    h2[0]->SetMinimum(0);

    h1[0]->SetMarkerStyle(20);
    h1[0]->SetMarkerSize(0.8);
    h1[0]->SetMarkerColor(kBlack);
    h1[0]->SetStats(0);
    h1[0]->Draw("E");

    h2[0]->SetMarkerStyle(21);
    h2[0]->SetMarkerSize(0.8);
    h2[0]->SetMarkerColor(kRed);
    h2[0]->SetStats(0);
    h2[0]->Draw("E SAME");

    // Define the Gaussian + quadratic background fit functions
    TF1 *fit1_int = new TF1("fit1_integrated", "gaus(0) + pol4(3)", -0.005, 0.005);
    TF1 *fit2_int = new TF1("fit2_integrated", "gaus(0) + pol4(3)", -0.005, 0.005);

    // Set the initial guesses and parameter limits for the fit
    fit1_int->SetParameters(0.5 * maxVal1, 0, 0.0025);  // Initial guesses: amplitude, mu, sigma
    fit1_int->SetParLimits(1, -0.0025, 0.0025);  // mu limits
    fit1_int->SetParLimits(2, 0, 0.005);         // sigma limits

    fit2_int->SetParameters(0.5 * maxVal2, 0, 0.0025);  // Initial guesses: amplitude, mu, sigma
    fit2_int->SetParLimits(1, -0.0025, 0.0025);  // mu limits
    fit2_int->SetParLimits(2, 0, 0.005);         // sigma limits

    fit1_int->SetLineWidth(1);
    fit2_int->SetLineWidth(1);
    h1[0]->Fit(fit1_int, "Q");
    h2[0]->Fit(fit2_int, "Q");

    fit1_int->SetLineColor(kBlack);
    fit1_int->Draw("SAME");
    fit2_int->SetLineColor(kRed);
    fit2_int->Draw("SAME");

    TLegend *legend_int = new TLegend(0.25, 0.75, 0.9, 0.9);
    legend_int->SetTextSize(0.03);
    legend_int->AddEntry(h1[0], Form("Uncorrected: #mu=%.4f, #sigma=%.4f", fit1_int->GetParameter(1), fit1_int->GetParameter(2)), "lep");
    legend_int->AddEntry(h2[0], Form("Corrected: #mu=%.4f, #sigma=%.4f", fit2_int->GetParameter(1), fit2_int->GetParameter(2)), "lep");
    legend_int->Draw();

    h1[0]->GetXaxis()->SetTitle("M_{x2} (GeV^{2})");
    h1[0]->GetYaxis()->SetTitle("Counts");

    // Draw the remaining subplots for the binned cases
    for (int i = 1; i <= nBins; ++i) {
        c1->cd(i + 1)->SetLeftMargin(0.15);
        c1->cd(i + 1)->SetBottomMargin(0.15);

        maxVal1 = h1[i]->GetMaximum();
        maxVal2 = h2[i]->GetMaximum();
        maxVal = std::max(maxVal1, maxVal2);
        h1[i]->SetMaximum(1.7 * maxVal);
        h1[i]->SetMinimum(0);
        h2[i]->SetMaximum(1.7 * maxVal);
        h2[i]->SetMinimum(0);

        h1[i]->SetMarkerStyle(20);
        h1[i]->SetMarkerSize(0.8);
        h1[i]->SetMarkerColor(kBlack);
        h1[i]->SetStats(0);
        h1[i]->Draw("E");

        h2[i]->SetMarkerStyle(21);
        h2[i]->SetMarkerSize(0.8);
        h2[i]->SetMarkerColor(kRed);
        h2[i]->SetStats(0);
        h2[i]->Draw("E SAME");

        // Define the Gaussian + quadratic background fit functions for each bin
        TF1 *fit1 = new TF1(Form("fit1_%d", i), "gaus(0) + pol4(3)", -0.005, 0.005);
        TF1 *fit2 = new TF1(Form("fit2_%d", i), "gaus(0) + pol4(3)", -0.005, 0.005);

        // Set the initial guesses and parameter limits for the fit
        fit1->SetParameters(0.5 * maxVal1, 0, 0.0025);  // Initial guesses: amplitude, mu, sigma
        fit1->SetParLimits(1, -0.0025, 0.0025);  // mu limits
        fit1->SetParLimits(2, 0, 0.005);         // sigma limits

        fit2->SetParameters(0.5 * maxVal2, 0, 0.0025);  // Initial guesses: amplitude, mu, sigma
        fit2->SetParLimits(1, -0.0025, 0.0025);  // mu limits
        fit2->SetParLimits(2, 0, 0.005);         // sigma limits

        fit1->SetLineWidth(1);
        fit2->SetLineWidth(1);
        h1[i]->Fit(fit1, "Q");
        h2[i]->Fit(fit2, "Q");

        fit1->SetLineColor(kBlack);
        fit1->Draw("SAME");
        fit2->SetLineColor(kRed);
        fit2->Draw("SAME");

        mu1_values[i - 1] = fit1->GetParameter(1);
        sigma1_values[i - 1] = fit1->GetParameter(2);
        mu2_values[i - 1] = fit2->GetParameter(1);
        sigma2_values[i - 1] = fit2->GetParameter(2);

        TLegend *legend = new TLegend(0.25, 0.75, 0.9, 0.9); // Adjusted the legend position to the top right corner
        legend->SetTextSize(0.03); // Decrease the font size in the legend
        legend->AddEntry(h1[i], Form("Uncorrected: #mu=%.4f, #sigma=%.4f", mu1_values[i - 1], sigma1_values[i - 1]), "lep");
        legend->AddEntry(h2[i], Form("Corrected: #mu=%.4f, #sigma=%.4f", mu2_values[i - 1], sigma2_values[i - 1]), "lep");
        legend->Draw();

        h1[i]->GetXaxis()->SetTitle("M_{x2} (GeV^{2})");
        h1[i]->GetYaxis()->SetTitle("Counts");
    }

    // Plot TGraphErrors in the last subplot (12th pad)
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);

    TGraphErrors *gr1 = new TGraphErrors(nBins, theta_mean, mu1_values, 0, sigma1_values);
    TGraphErrors *gr2 = new TGraphErrors(nBins, theta_mean, mu2_values, 0, sigma2_values);

    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.8);
    gr1->SetMarkerColor(kBlack);
    gr1->Draw("AP");

    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(0.8);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("P SAME");

    // Add dashed gray line at y = 0
    TLine *line = new TLine(0, 0, 70, 0);
    line->SetLineColor(kGray);
    line->SetLineStyle(2); // Dashed line
    line->Draw("SAME");

    // Customize the last subplot
    gr1->SetTitle(""); // Remove the "Graph" title
    gr1->GetXaxis()->SetTitle("#theta"); // Set x-axis label
    gr1->GetYaxis()->SetTitle("#mu (GeV^{2})"); // Set y-axis label
    gr1->GetXaxis()->SetTitleSize(0.05); // Match font size with other plots
    gr1->GetYaxis()->SetTitleSize(0.05); // Match font size with other plots
    gr1->GetXaxis()->SetLimits(25, 60); // Set x-axis range
    gr1->GetYaxis()->SetRangeUser(-0.005, 0.005); // Set y-axis range

    // Add legend to the last plot, positioned in the top right
    TLegend *legend12 = new TLegend(0.6, 0.75, 0.9, 0.9); // Adjusted for horizontal and vertical size
    legend12->SetTextSize(0.03);
    legend12->AddEntry(gr1, "Uncorrected", "lep");
    legend12->AddEntry(gr2, "Corrected", "lep");
    legend12->Draw();

    // Save the canvas as a PDF
    TString outputFileName = "output/elastic_" + TString(titleSuffix) + "_energy_loss_validation.pdf";
    c1->SaveAs(outputFileName);

    // Clean up
    delete c1;
    delete f1;
    delete f2;
}

void plot_pi0_energy_loss_validation(const char* file1, const char* file2, const char* titleSuffix) {
    // Open the ROOT files
    TFile *f1 = new TFile(file1);
    TFile *f2 = new TFile(file2);

    // Access the "PhysicsEvents" trees
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    // Prepare variables for reading the branches
    Double_t p1_theta1, p1_theta2;
    Double_t Mx1_1, Mx1_2;

    // Set branch addresses
    tree1->SetBranchAddress("p1_theta", &p1_theta1);
    tree1->SetBranchAddress("Mx", &Mx1_1);

    tree2->SetBranchAddress("p1_theta", &p1_theta2);
    tree2->SetBranchAddress("Mx", &Mx1_2);

    // Create canvas and divide it into 3x4 subplots
    TCanvas *c1 = new TCanvas("c1", "pi0 Energy Loss Validation", 1200, 900);
    c1->Divide(4, 3);

    // Set the theta bins and corresponding histogram ranges
    const int nBins = 10; // 10 theta bins + 1 fully integrated case
    Double_t thetaBins[nBins + 1] = {5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65};  // 5 to 65 degrees

    TH1D *h1[nBins + 1]; // +1 for the fully integrated case
    TH1D *h2[nBins + 1]; // +1 for the fully integrated case

    Double_t mu1_values[nBins], sigma1_values[nBins];
    Double_t mu2_values[nBins], sigma2_values[nBins];
    Double_t theta_sum[nBins] = {0.0};
    Int_t theta_count[nBins] = {0};
    Double_t theta_mean[nBins] = {0.0};

    // Create histograms for each theta bin with 50 bins
    h1[0] = new TH1D("h1_integrated", Form("Integrated #theta [5, 65] %s", titleSuffix), 30, -1.0, 2.22);
    h2[0] = new TH1D("h2_integrated", Form("Integrated #theta [5, 65] %s", titleSuffix), 30, -1.0, 2.22);

    for (int i = 0; i < nBins; ++i) {
        h1[i + 1] = new TH1D(Form("h1_%d", i), Form("Mx1 for #theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 30, -1.0, 2.22);
        h2[i + 1] = new TH1D(Form("h2_%d", i), Form("Mx1 for #theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 30, -1.0, 2.22);
    }

    // Fill the histograms and calculate theta means
    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t j = 0; j < nEntries1; ++j) {
        tree1->GetEntry(j);
        Double_t thetaDeg1 = p1_theta1 * (180.0 / TMath::Pi()); // Convert to degrees
        if (thetaDeg1 >= 5 && thetaDeg1 < 65) {
            h1[0]->Fill(Mx1_1); // Fully integrated case
        }
        for (int i = 0; i < nBins; ++i) {
            if (thetaDeg1 >= thetaBins[i] && thetaDeg1 < thetaBins[i + 1]) {
                h1[i + 1]->Fill(Mx1_1);
                theta_sum[i] += thetaDeg1;
                theta_count[i]++;
            }
        }
    }

    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t j = 0; j < nEntries2; ++j) {
        tree2->GetEntry(j);
        Double_t thetaDeg2 = p1_theta2 * (180.0 / TMath::Pi()); // Convert to degrees
        if (thetaDeg2 >= 5 && thetaDeg2 < 65) {
            h2[0]->Fill(Mx1_2); // Fully integrated case
        }
        for (int i = 0; i < nBins; ++i) {
            if (thetaDeg2 >= thetaBins[i] && thetaDeg2 < thetaBins[i + 1]) {
                h2[i + 1]->Fill(Mx1_2);
                theta_sum[i] += thetaDeg2;
                theta_count[i]++;
            }
        }
    }

    // Calculate mean theta values
    for (int i = 0; i < nBins; ++i) {
        if (theta_count[i] > 0) {
            theta_mean[i] = theta_sum[i] / theta_count[i];
        } else {
            theta_mean[i] = 0.5 * (thetaBins[i] + thetaBins[i + 1]); // Default to midpoint if no entries
        }
    }

    // Draw the fully integrated case in the first subplot
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);

    Double_t maxVal1 = h1[0]->GetMaximum();
    Double_t maxVal2 = h2[0]->GetMaximum();
    Double_t maxVal = std::max(maxVal1, maxVal2);
    h1[0]->SetMaximum(1.7 * maxVal);
    h1[0]->SetMinimum(0);
    h2[0]->SetMaximum(1.7 * maxVal);
    h2[0]->SetMinimum(0);

    h1[0]->SetMarkerStyle(20);
    h1[0]->SetMarkerSize(0.8);
    h1[0]->SetMarkerColor(kBlack);
    h1[0]->SetStats(0);
    h1[0]->Draw("E");

    h2[0]->SetMarkerStyle(21);
    h2[0]->SetMarkerSize(0.8);
    h2[0]->SetMarkerColor(kRed);
    h2[0]->SetStats(0);
    h2[0]->Draw("E SAME");

    TF1 *fit1_int = new TF1("fit1_integrated", "gaus(0) + pol3(3)", 0.0, 0.22);
    TF1 *fit2_int = new TF1("fit2_integrated", "gaus(0) + pol3(3)", 0.0, 0.22);

    fit1_int->SetParameters(0.8 * maxVal1, 0.135, 0.01); // Initial guesses for amplitude, mean, sigma
    fit2_int->SetParameters(0.8 * maxVal2, 0.135, 0.01);

    fit1_int->SetLineWidth(1);
    fit2_int->SetLineWidth(1);
    h1[0]->Fit(fit1_int, "Q");
    h2[0]->Fit(fit2_int, "Q");

    fit1_int->SetLineColor(kBlack);
    fit1_int->Draw("SAME");
    fit2_int->SetLineColor(kRed);
    fit2_int->Draw("SAME");

    TLegend *legend_int = new TLegend(0.25, 0.75, 0.9, 0.9);
    legend_int->SetTextSize(0.03);
    legend_int->AddEntry(h1[0], Form("Uncorrected: #mu=%.3f, #sigma=%.3f", fit1_int->GetParameter(1), fit1_int->GetParameter(2)), "lep");
    legend_int->AddEntry(h2[0], Form("Corrected: #mu=%.3f, #sigma=%.3f", fit2_int->GetParameter(1), fit2_int->GetParameter(2)), "lep");
    legend_int->Draw();

    h1[0]->GetXaxis()->SetTitle("M_{xp} (GeV)");
    h1[0]->GetYaxis()->SetTitle("Counts");

    // Draw the remaining subplots for the binned cases
    for (int i = 1; i <= nBins; ++i) {
        c1->cd(i + 1)->SetLeftMargin(0.15);
        c1->cd(i + 1)->SetBottomMargin(0.15);

        maxVal1 = h1[i]->GetMaximum();
        maxVal2 = h2[i]->GetMaximum();
        maxVal = std::max(maxVal1, maxVal2);
        h1[i]->SetMaximum(1.7 * maxVal);
        h1[i]->SetMinimum(0);
        h2[i]->SetMaximum(1.7 * maxVal);
        h2[i]->SetMinimum(0);

        h1[i]->SetMarkerStyle(20);
        h1[i]->SetMarkerSize(0.8);
        h1[i]->SetMarkerColor(kBlack);
        h1[i]->SetStats(0);
        h1[i]->Draw("E");

        h2[i]->SetMarkerStyle(21);
        h2[i]->SetMarkerSize(0.8);
        h2[i]->SetMarkerColor(kRed);
        h2[i]->SetStats(0);
        h2[i]->Draw("E SAME");

        TF1 *fit1 = new TF1(Form("fit1_%d", i), "gaus(0) + pol3(3)", 0.0, 0.22);
        TF1 *fit2 = new TF1(Form("fit2_%d", i), "gaus(0) + pol3(3)", 0.0, 0.22);

        fit1->SetParameters(0.8 * maxVal1, 0.135, 0.01);
        fit2->SetParameters(0.8 * maxVal2, 0.135, 0.01);

        fit1->SetLineWidth(1);
        fit2->SetLineWidth(1);
        h1[i]->Fit(fit1, "Q");
        h2[i]->Fit(fit2, "Q");

        fit1->SetLineColor(kBlack);
        fit1->Draw("SAME");
        fit2->SetLineColor(kRed);
        fit2->Draw("SAME");

        mu1_values[i - 1] = fit1->GetParameter(1);
        sigma1_values[i - 1] = fit1->GetParameter(2);
        mu2_values[i - 1] = fit2->GetParameter(1);
        sigma2_values[i - 1] = fit2->GetParameter(2);

        TLegend *legend = new TLegend(0.25, 0.75, 0.9, 0.9); // Adjusted the legend position to the top right corner
        legend->SetTextSize(0.03); // Decrease the font size in the legend
        legend->AddEntry(h1[i], Form("Uncorrected: #mu=%.3f, #sigma=%.3f", mu1_values[i - 1], sigma1_values[i - 1]), "lep");
        legend->AddEntry(h2[i], Form("Corrected: #mu=%.3f, #sigma=%.3f", mu2_values[i - 1], sigma2_values[i - 1]), "lep");
        legend->Draw();

        h1[i]->GetXaxis()->SetTitle("M_{xp} (GeV)");
        h1[i]->GetYaxis()->SetTitle("Counts");
    }

    // Plot TGraphErrors in the last subplot (12th pad)
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);

    TGraphErrors *gr1 = new TGraphErrors(nBins, theta_mean, mu1_values, 0, sigma1_values);
    TGraphErrors *gr2 = new TGraphErrors(nBins, theta_mean, mu2_values, 0, sigma2_values);

    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.8);
    gr1->SetMarkerColor(kBlack);
    gr1->Draw("AP");

    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(0.8);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("P SAME");

    // Add dashed gray line at y = pi0_mass
    TLine *line = new TLine(0, 0.135, 70, 0.135);
    line->SetLineColor(kGray);
    line->SetLineStyle(2); // Dashed line
    line->Draw("SAME");

    // Customize the last subplot
    gr1->SetTitle(""); // Remove the "Graph" title
    gr1->GetXaxis()->SetTitle("#theta"); // Set x-axis label
    gr1->GetYaxis()->SetTitle("#mu (GeV)");    // Set y-axis label
    gr1->GetXaxis()->SetTitleSize(0.05); // Match font size with other plots
    gr1->GetYaxis()->SetTitleSize(0.05); // Match font size with other plots
    gr1->GetXaxis()->SetLimits(5, 65); // Set x-axis range
    gr1->GetYaxis()->SetRangeUser(0.13, 0.14); // Set y-axis range

    // Add legend to the last plot, positioned in the top right
    TLegend *legend12 = new TLegend(0.6, 0.75, 0.9, 0.9); // Adjusted for horizontal and vertical size
    legend12->SetTextSize(0.03);
    legend12->AddEntry(gr1, "Uncorrected", "lep");
    legend12->AddEntry(gr2, "Corrected", "lep");
    legend12->Draw();

    // Save the canvas as a PDF
    TString outputFileName = "output/pi0_" + TString(titleSuffix) + "_energy_loss_validation.pdf";
    c1->SaveAs(outputFileName);

    // Clean up
    delete c1;
    delete f1;
    delete f2;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <file1.root> <file2.root> <titleSuffix>" << std::endl;
        return 1;
    }

    plot_dvcs_energy_loss_validation(argv[1], argv[2], argv[3]);
    // plot_rho0_energy_loss_validation(argv[1], argv[2], argv[3]);
    // plot_pi0_energy_loss_validation(argv[1], argv[2], argv[3]);
    // plot_elastic_energy_loss_validation(argv[1], argv[2], argv[3]);
    return 0;
}