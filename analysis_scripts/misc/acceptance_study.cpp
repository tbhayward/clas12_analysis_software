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

// New function to create the plot for a given exclusion percentage
void createPlotWithExclusion(const std::vector<double>& phiVec, double B, double C, int nBins, double exclusionPercent, int canvasIndex, TCanvas* masterCanvas) {
    // Calculate exclusion based on percentage
    long unsigned int N = phiVec.size();
    int excludeFromEachSide = N * (exclusionPercent / 100.0) / 2;
    std::vector<double> includedPhi, excludedPhi;

    // Sort phi values to apply exclusion symmetrically
    std::vector<double> sortedPhi = phiVec;
    std::sort(sortedPhi.begin(), sortedPhi.end());
    for (long unsigned int i = excludeFromEachSide; i < N - excludeFromEachSide; ++i) {
        includedPhi.push_back(sortedPhi[i]);
    }
    for (long unsigned int i = 0; i < excludeFromEachSide; ++i) {
        excludedPhi.push_back(sortedPhi[i]);
        excludedPhi.push_back(sortedPhi[N - 1 - i]);
    }

    double binWidth = (2*TMath::Pi()) / nBins;
    TGraphErrors *graphIncluded = new TGraphErrors();
    TGraphErrors *graphExcluded = new TGraphErrors();
    graphExcluded->SetMarkerStyle(4); // Empty circles for excluded points

    // Fill the graphs for included and excluded phi values
    std::vector<int> binCountsIncluded(nBins, 0), binCountsExcluded(nBins, 0);
    for (double phi : includedPhi) {
        int binIndex = static_cast<int>(phi / binWidth);
        binCountsIncluded[binIndex]++;
    }
    for (double phi : excludedPhi) {
        int binIndex = static_cast<int>(phi / binWidth);
        binCountsExcluded[binIndex]++;
    }

    for (int i = 0; i < nBins; ++i) {
        double binCenter = binWidth * i + binWidth / 2;
        if (binCountsIncluded[i] > 0) {
            double binError = sqrt(binCountsIncluded[i]); // Statistical error
            graphIncluded->SetPoint(graphIncluded->GetN(), binCenter, binCountsIncluded[i]);
            graphIncluded->SetPointError(graphIncluded->GetN() - 1, 0, binError);
        }
        if (binCountsExcluded[i] > 0) {
            graphExcluded->SetPoint(graphExcluded->GetN(), binCenter, binCountsExcluded[i]);
        }
    }

    // Create and set up the fit function
    TF1 *fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
    fitFunc->SetParameters(1, B, C);

    // Drawing
    masterCanvas->cd(canvasIndex);
    graphIncluded->SetMarkerStyle(20); // Filled circles for included points
    graphIncluded->Draw("AP");
    graphIncluded->Fit(fitFunc, "Q");
    graphExcluded->Draw("P SAME");
    fitFunc->Draw("SAME");

    // Customize axes and fit function appearance
    graphIncluded->GetXaxis()->SetLimits(0, TMath::TwoPi());
    graphIncluded->GetYaxis()->SetRangeUser(0, *std::max_element(binCountsIncluded.begin(), binCountsIncluded.end()) * 1.2);
    graphIncluded->GetXaxis()->SetTitle("#phi");
    graphIncluded->GetYaxis()->SetTitle("Counts");
    fitFunc->SetLineColor(kRed);

    // Optionally, add TPaveText or other annotations as needed
    // For example, showing the exclusion percentage on the plot
    TPaveText *pt = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
    pt->SetFillColor(0);
    pt->AddText(Form("Exclusion: %.0f%%", exclusionPercent));
    pt->Draw();
}

int acceptanceStudy(double B, double C) {
    const long unsigned int N = 1e4; // Number of points to generate
    std::vector<double> phiVec;
    generateData(B, C, N, phiVec);

    const int nBins = 24;
    TCanvas *masterCanvas = new TCanvas("masterCanvas", "Phi Distribution Fits", 1200, 800);
    masterCanvas->Divide(3, 2); // 3x2 grid

    // Loop over exclusions: 0%, 5%, 10%, 15%, 20%, 25%
    for (int i = 0; i <= 5; ++i) {
        double exclusionPercent = i * 5.0;
        createPlotWithExclusion(phiVec, B, C, nBins, exclusionPercent, i+1, masterCanvas);
    }

    // Save the canvas with all plots
    std::ostringstream filename;
    filename << "output/acceptance_study_all_B=" << B << "_C=" << C << ".png";
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
