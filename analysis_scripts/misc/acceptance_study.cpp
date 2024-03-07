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

void generateData(double B, double C, int N, double* phi, double* values) {
    TRandom3 rand(0); // Seed for reproducibility
    for (int i = 0; i < N; ++i) {
        phi[i] = rand.Uniform(0, 2*TMath::Pi());
        double A = 1.0; // Assuming A is constant
        values[i] = rand.Gaus(A * (1 + B * TMath::Cos(phi[i]) + C * TMath::Cos(2 * phi[i])), 0.0); // Adding some noise
    }
}

int acceptanceStudy(double B, double C) {
    const long unsigned int N = 5e5; // Number of points
    double phi[N], values[N];
    generateData(B, C, N, phi, values);
    std::cout << "Data generated." << std::endl;

    // Binning data
    const int nBins = 24;
    double binWidth = (2*TMath::Pi()) / nBins;
    TGraphErrors *graph = new TGraphErrors();
    int binCounts[nBins] = {0};
    double binValues[nBins] = {0.0};
    double binErrors[nBins] = {0.0};
    for (int i = 0; i < N; ++i) {
        int binIndex = static_cast<int>(phi[i] / binWidth);
        binCounts[binIndex]++;
        binValues[binIndex]+=values[i];
    }
    for (int i = 0; i < nBins; ++i) {
        if (binCounts[i] > 0) {
            binValues[i] = binCounts[i]; // Average value for the bin
            binErrors[i] = sqrt(binCounts[i]); // Statistical error
            graph->SetPoint(i, binWidth*i + binWidth/2, binValues[i]);
            graph->SetPointError(i, 0, binErrors[i]);
        }
    }

    // Fitting
    TF1 *fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
    fitFunc->SetParameters(1, B, C); // Initial parameters guess
    graph->Fit(fitFunc, "Q"); // Quiet mode

    // Creating canvas and drawing
    TCanvas *c = new TCanvas("c", "Fitting with Chi2 Minimization", 800, 600);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlack);
    graph->Draw("APE");
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");

    // Adding TPaveText for fit parameters
    TPaveText *pt = new TPaveText(0.1, 0.7, 0.5, 0.9, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(Form("A = %.3f #pm %.3f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
    pt->AddText(Form("B = %.3f #pm %.3f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
    pt->AddText(Form("C = %.3f #pm %.3f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));
    pt->Draw();

    // Saving the plot
    std::ostringstream filename;
    filename << "output/acceptance_study_B=" << B << "_C=" << C << ".png";
    c->SaveAs(filename.str().c_str());

    delete c;
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
