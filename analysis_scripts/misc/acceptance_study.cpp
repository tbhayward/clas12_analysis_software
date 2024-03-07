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

int acceptanceStudy(double B, double C) {
    const long unsigned int N = 1e5; // Target number of points
    std::vector<double> phiVec;
    generateData(B, C, N, phiVec);
    std::cout << "Data generated. Number of accepted phi values: " << phiVec.size() << std::endl;

    // Binning data
    const int nBins = 24;
    double binWidth = (2*TMath::Pi()) / nBins;
    TGraphErrors *graph = new TGraphErrors();
    int binCounts[nBins] = {0};

    for (double phi : phiVec) {
        int binIndex = static_cast<int>(phi / binWidth);
        binCounts[binIndex]++;
    }

    for (int i = 0; i < nBins; ++i) {
        double binCenter = binWidth * i + binWidth / 2;
        double binError = sqrt(binCounts[i]); // Statistical error
        if (binCounts[i] > 0) {
            graph->SetPoint(graph->GetN(), binCenter, binCounts[i]);
            graph->SetPointError(graph->GetN() - 1, 0, binError);
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
    graph->GetXaxis()->SetLimits(0, TMath::TwoPi()); // Set x-axis range to 0 to 2pi
    graph->GetYaxis()->SetRangeUser(*std::min_element(binCounts, binCounts + nBins) * 0.8, *std::max_element(binCounts, binCounts + nBins) * 1.2); 
    graph->GetXaxis()->SetTitle("#phi");
    graph->GetYaxis()->SetTitle("Counts");

    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");

    // Adding TPaveText for fit parameters
	TPaveText *pt = new TPaveText(0.1, 0.7, 0.5, 0.9, "NDC");
	pt->SetFillColor(0);
	pt->SetTextAlign(12);
	// Remove or comment out the line for A
	// pt->AddText(Form("A = %.3f #pm %.3f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
	pt->AddText(Form("B = %.3f #pm %.3f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
	pt->AddText(Form("C = %.3f #pm %.3f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));
	// Add a line for chi2/ndf
	double chi2 = fitFunc->GetChisquare();
	double ndf = fitFunc->GetNDF(); // number of degrees of freedom
	pt->AddText(Form("#chi^{2}/ndf = %.3f", chi2 / ndf));
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
