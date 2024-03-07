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

void plotForExclusion(const std::vector<double>& phiVec, double B, double C, int canvasIndex, int binsToExclude, TCanvas* masterCanvas) {
    double binWidth = (2 * TMath::Pi()) / 24;
    TGraphErrors *graphIncluded = new TGraphErrors();
    TGraphErrors *graphExcluded = new TGraphErrors();
    graphExcluded->SetMarkerStyle(4); // Open circles for excluded points

    // Initialize bin counts
    std::vector<int> binCountsIncluded(24, 0), binCountsExcluded(24, 0);

    for (double phi : phiVec) {
        int binIndex = static_cast<int>(phi / binWidth);
        if (binIndex >= binsToExclude && binIndex < 24 - binsToExclude) {
            binCountsIncluded[binIndex]++;
        } else {
            binCountsExcluded[binIndex]++;
        }
    }

    // Set points for included and excluded graphs
    for (int i = 0; i < 24; ++i) {
        double binCenter = binWidth * i + binWidth / 2;
        if (binCountsIncluded[i] > 0) {
            graphIncluded->SetPoint(graphIncluded->GetN(), binCenter, binCountsIncluded[i]);
            graphIncluded->SetPointError(graphIncluded->GetN() - 1, 0, TMath::Sqrt(binCountsIncluded[i]));
        }
        if (binCountsExcluded[i] > 0) {
            graphExcluded->SetPoint(graphExcluded->GetN(), binCenter, binCountsExcluded[i]);
            graphExcluded->SetPointError(graphExcluded->GetN() - 1, 0, TMath::Sqrt(binCountsExcluded[i]));
        }
    }

    double exclusionPercentage = binsToExclude * 100.0 / 12.0; // 12.0 because 24 bins total, and we exclude from both sides

    masterCanvas->cd(canvasIndex);
    graphIncluded->SetMarkerStyle(20); // Filled circles for included points
    graphIncluded->Draw("AP");
    graphExcluded->Draw("P SAME");
    graphExcluded->GetXaxis()->SetLimits(0, TMath::TwoPi());
    graphExcluded->GetYaxis()->SetLimits(3000, 7000);

    // Perform the fit within the limited range
	TF1 *fitFuncLimited = new TF1("fitFuncLimited", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", binWidth * binsToExclude, 2*TMath::Pi() - binWidth * binsToExclude);
	fitFuncLimited->SetParameters(1, B, C);
	graphIncluded->Fit(fitFuncLimited, "Q R");

	// Define a new TF1 that uses the parameters from the fit but extends across the full range
	TF1 *fitFuncFullRange = new TF1("fitFuncFullRange", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
	fitFuncFullRange->SetParameters(fitFuncLimited->GetParameters());
	fitFuncFullRange->SetLineColor(kRed);

	masterCanvas->cd(canvasIndex);
	graphIncluded->Draw("AP");
	graphExcluded->Draw("P SAME");
	fitFuncFullRange->Draw("SAME");

    // Customize axis limits and labels
    graphIncluded->GetXaxis()->SetTitle("#phi");
    graphIncluded->GetYaxis()->SetTitle("Counts");
    graphIncluded->GetXaxis()->SetLimits(0, TMath::TwoPi());
    graphIncluded->GetYaxis()->SetLimits(3000, 7000);

    // Adding TPaveText for fit parameters and chi2/ndf
    TPaveText *pt = new TPaveText(0.1, 0.65, 0.5, 0.9, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(Form("Exclusion: %.1f%%", exclusionPercentage));
    pt->AddText(Form("A_{UU}^{cos#phi} = %.3f #pm %.3f", fitFuncLimited->GetParameter(1), fitFuncLimited->GetParError(1)));
    pt->AddText(Form("A_{UU}^{cos2#phi} = %.3f #pm %.3f", fitFuncLimited->GetParameter(2), fitFuncLimited->GetParError(2)));
    double chi2 = fitFuncLimited->GetChisquare();
    double ndf = fitFuncLimited->GetNDF();
    pt->AddText(Form("#chi^{2}/ndf = %.3f", chi2 / ndf));
    pt->Draw();
}

int acceptanceStudy(double B, double C) {
    const long unsigned int N = 1e5; // Increased number of points
    std::vector<double> phiVec;
    generateData(B, C, N, phiVec);

    TCanvas *masterCanvas = new TCanvas("masterCanvas", "Phi Distribution Fits", 1200, 800);
    masterCanvas->Divide(3, 2);

    // Loop over different numbers of bins to exclude
    for (int i = 0; i <= 5; ++i) {
        plotForExclusion(phiVec, B, C, i + 1, i, masterCanvas);
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
