#include <TMath.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

void generateData(double B, double C, long unsigned int N, std::vector<double>& phiVec) {
    TRandom3 rand(0); // Seed for reproducibility
    double A = 1.0; // Assuming A is constant
    long unsigned int accepted = 0;
    while (accepted < N) {
        double phi = rand.Uniform(0, 2*TMath::Pi());
        double value = A * (1 + B * TMath::Cos(phi) + C * TMath::Cos(2 * phi));
        double randVal = rand.Uniform(0, A * (1 + fabs(B) + fabs(C))); // Max possible value

        if (randVal < value) {
            phiVec.push_back(phi);
            accepted++;
        }
    }
}

void plotForExclusion(const std::vector<double>& phiVec, double B, double C, int canvasIndex, double exclusionFraction, TCanvas* masterCanvas) {
    double fitRangeMin = 2 * TMath::Pi() * exclusionFraction / 2.0;
    double fitRangeMax = 2 * TMath::Pi() * (1 - exclusionFraction / 2.0);
    double binWidth = (2*TMath::Pi()) / 24;

    TGraphErrors *graphIncluded = new TGraphErrors(); // Graph for included data points
    TGraphErrors *graphExcluded = new TGraphErrors(); // Graph for excluded data points
    graphExcluded->SetMarkerStyle(4); // Open circles for excluded points

    int binIndexIncluded = 0, binIndexExcluded = 0;
    for (double phi : phiVec) {
        double yValue = 1.0; // Assuming a constant value for illustration, replace with your actual y value or bin count logic
        double error = sqrt(yValue); // Example error calculation
        if (phi >= fitRangeMin && phi <= fitRangeMax) {
            graphIncluded->SetPoint(binIndexIncluded, phi, yValue);
            graphIncluded->SetPointError(binIndexIncluded, 0, error);
            binIndexIncluded++;
        } else {
            graphExcluded->SetPoint(binIndexExcluded, phi, yValue);
            graphExcluded->SetPointError(binIndexExcluded, 0, error);
            binIndexExcluded++;
        }
    }

    // Fitting only the included data points
    TF1 *fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", fitRangeMin, fitRangeMax);
    fitFunc->SetParameters(1, B, C); // Example parameters

    masterCanvas->cd(canvasIndex);
    graphIncluded->Draw("AP");
    graphIncluded->Fit(fitFunc, "Q");
    graphExcluded->Draw("P SAME");

    // Setting up the plot
    graphIncluded->GetXaxis()->SetTitle("#phi");
    graphIncluded->GetYaxis()->SetTitle("Counts");
    graphIncluded->GetXaxis()->SetLimits(0, 2*TMath::Pi());
    // Adjust y-axis range as needed
    // graphIncluded->GetYaxis()->SetRangeUser(0, maxYValue);

    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("SAME");

    // TPaveText for parameters
    TPaveText *pt = new TPaveText(0.1, 0.65, 0.5, 0.9, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(Form("Exclusion: %.0f%%", exclusionFraction * 100));
    pt->AddText(Form("A_{UU}^{cos#phi} = %.3f #pm %.3f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
    pt->AddText(Form("A_{UU}^{cos2#phi} = %.3f #pm %.3f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    pt->AddText(Form("#chi^{2}/ndf = %.3f", chi2 / ndf));
    pt->Draw();
}


int acceptanceStudy(double B, double C) {
    const long unsigned int N = 1e5;
    std::vector<double> phiVec;
    generateData(B, C, N, phiVec);

    TCanvas *masterCanvas = new TCanvas("masterCanvas", "Phi Distribution Fits", 1200, 800);
    masterCanvas->Divide(3, 2);

    // Loop over different exclusions
    for (int i = 0; i <= 5; ++i) {
        double exclusionFraction = i * 0.05; // 0%, 5%, 10%, 15%, 20%, 25%
        plotForExclusion(phiVec, B, C, i + 1, exclusionFraction, masterCanvas);
    }

    std::ostringstream filename;
    filename << "output/acceptance_study_B=" << B << "_C=" << C << ".png";
    masterCanvas->SaveAs(filename.str().c_str());

    delete masterCanvas;
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <B> <C>" << std::endl;
        return 1;
    }
    double B = atof(argv[1]);
    double C = atof(argv[2]);

    return acceptanceStudy(B, C);
}
