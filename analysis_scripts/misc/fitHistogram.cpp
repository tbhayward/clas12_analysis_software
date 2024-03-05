#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>

// Define the fitting function
double fitFunction(double *x, double *par) {
    double linearPart = par[0] + par[1] * x[0];
    double gaussianPart = par[2] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[3]) / par[4], 2));
    return linearPart + gaussianPart;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <ROOT file path>" << std::endl;
        return 1;
    }

    // Load the ROOT file
    TFile *file = TFile::Open(argv[1]);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }

    // Retrieve the histogram
    TH1F *hist = (TH1F*)file->Get("canvas"); // Ensure the histogram name matches
    if (!hist) {
        std::cerr << "Histogram not found" << std::endl;
        return 1;
    }

    // Create the fit function
    TF1 *fitFunc = new TF1("fitFunc", fitFunction, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 5);
    fitFunc->SetParameters(1, 1, 1, hist->GetMean(), hist->GetStdDev()); // Initial guesses
    fitFunc->SetParNames("LinearConst", "LinearSlope", "GaussAmplitude", "GaussMean", "GaussSigma");

    // Constrain the amplitude of the Gaussian to be > 0
    fitFunc->SetParLimits(2, 0.0001, 1e6); // Assuming amplitude is the third parameter

    // Fit the histogram
    hist->Fit(fitFunc, "R");

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Fit Histogram", 800, 600);
    hist->Draw();
    fitFunc->Draw("same");

    // Save the canvas as a PNG image
    c1->SaveAs("/Users/tbhayward/Desktop/fit_distro.png");

    // Cleanup
    delete c1;
    file->Close();
    delete file;

    return 0;
}
