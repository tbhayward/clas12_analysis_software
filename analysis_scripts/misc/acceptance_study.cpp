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
#include <utility> // For std::pair

// Global structure to hold deviations. 
// Each pair corresponds to an exclusion case, holding two vectors for deviations in B and C respectively.
std::vector<std::pair<std::vector<double>, std::vector<double>>> deviationsForCases(6); 

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

void setupCanvas(TCanvas* masterCanvas) {
    masterCanvas->Divide(3, 2, 0, 0); // Set gaps to 0

    for (int i = 1; i <= 6; ++i) {
        masterCanvas->cd(i);
        TPad *pad = (TPad*)gPad;
        
        pad->SetTopMargin((i <= 3) ? 0.02 : 0.00); // Smaller margin for top row
        pad->SetBottomMargin((i > 3) ? 0.15 : 0.00); // Larger bottom margin only for bottom row
        pad->SetLeftMargin(((i-1) % 3 == 0) ? 0.15 : 0.00); // Larger left margin only for first column
        pad->SetRightMargin(((i-1) % 3 == 2) ? 0.01 : 0.00); // Uniform right margin
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

    double exclusionPercentage = (binsToExclude / 24.0) * 100.0; // Update formula to accurately represent the exclusion percentage

    masterCanvas->cd(canvasIndex);
    graphIncluded->SetMarkerStyle(20); // Filled circles for included points
    graphIncluded->Draw("AP");
    graphExcluded->Draw("P SAME");
    graphExcluded->GetXaxis()->SetLimits(0, TMath::TwoPi());
    graphExcluded->GetYaxis()->SetLimits(13000, 17000);

    // Find the maximum value of the distribution for the initial guess of [0]
	double maxY = 0;
	for (int i = 0; i < graphIncluded->GetN(); ++i) {
	    double x, y;
	    graphIncluded->GetPoint(i, x, y);
	    if (y > maxY) maxY = y;
	}

	// Perform the fit within the limited range
	TF1 *fitFuncLimited = new TF1("fitFuncLimited", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", binWidth * binsToExclude, 2*TMath::Pi() - binWidth * binsToExclude);

	// Set initial parameters
	fitFuncLimited->SetParameter(0, maxY); // Initial guess for [0] is the max Y value
	fitFuncLimited->SetParameter(1, B);    // Initial guess for [1]
	fitFuncLimited->SetParameter(2, C);    // Initial guess for [2]

	// Set parameter limits for [1] and [2]
	fitFuncLimited->SetParLimits(1, -1, 1); // Constrain [1] to be between -1 and 1
	fitFuncLimited->SetParLimits(2, -1, 1); // Constrain [2] to be between -1 and 1

	// Fit the graph with these settings
	graphIncluded->Fit(fitFuncLimited, "Q R");

	// Extract fitted parameters and their errors
    double fittedB = fitFuncLimited->GetParameter(1);
    double fittedC = fitFuncLimited->GetParameter(2);
    double errB = fitFuncLimited->GetParError(1);
    double errC = fitFuncLimited->GetParError(2);

    // Calculate deviations in sigma
    double deviationSigmaB = (fittedB - B) / errB;
    double deviationSigmaC = (fittedC - C) / errC;

    // Store deviations
    deviationsForCases[canvasIndex-1].first.push_back(deviationSigmaB); // Store deviation for B
    deviationsForCases[canvasIndex-1].second.push_back(deviationSigmaC); // Store deviation for C

	// Define a new TF1 that uses the parameters from the fit but extends across the full range
	TF1 *fitFuncFullRange = new TF1("fitFuncFullRange", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
	fitFuncFullRange->SetParameters(fitFuncLimited->GetParameters());
	fitFuncFullRange->SetLineColor(kRed);

	masterCanvas->cd(canvasIndex);
	graphIncluded->Draw("AP");
	graphIncluded->GetYaxis()->SetRangeUser(0.8*fitFuncLimited->GetParameter(0), 1.3 *fitFuncLimited->GetParameter(0));
	graphExcluded->Draw("P SAME");
	fitFuncFullRange->Draw("SAME");

    // Customize axis limits and labels
    graphIncluded->GetXaxis()->SetTitle("#phi");
    graphIncluded->GetYaxis()->SetTitle("Counts");
    graphIncluded->GetXaxis()->SetTitleSize(0.05);
    graphIncluded->GetYaxis()->SetTitleSize(0.05);
    if (canvasIndex<5) {
    	graphIncluded->GetXaxis()->SetLimits(0, TMath::TwoPi());
    } else {
    	graphIncluded->GetXaxis()->SetLimits(0.001, TMath::TwoPi());
    }

    // Calculate the start position for TPaveText based on the column, directly here
    double textStartX = (canvasIndex % 3 == 1) ? 0.15 : 0.0; // Example: 0.15 for left column, adjust 0.12 for middle/right columns
    TPaveText *pt = new TPaveText(textStartX, 0.75, textStartX + 0.45, 1.0, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.04); // Adjust the text size if needed
    pt->AddText(Form("Exclusion: %.1f%%", exclusionPercentage));
	pt->AddText(Form("A_{UU}^{cos#phi} = %.3f #pm %.3f", fitFuncLimited->GetParameter(1), fitFuncLimited->GetParError(1)));
	pt->AddText(Form("A_{UU}^{cos2#phi} = %.3f #pm %.3f", fitFuncLimited->GetParameter(2), fitFuncLimited->GetParError(2)));
	double chi2 = fitFuncLimited->GetChisquare();
	double ndf = fitFuncLimited->GetNDF();
	pt->AddText(Form("#chi^{2}/ndf = %.3f", chi2 / ndf));
	pt->Draw();


    // After drawing the graphs and fit function, adjust axis titles visibility if necessary
    bool isLeftColumn = (canvasIndex % 3 == 1); // Adjust based on your actual layout
    bool isBottomRow = (canvasIndex > 3); // Adjust based on your actual layout

    if (!isLeftColumn) {
        graphIncluded->GetYaxis()->SetTitle("");
    }
    if (!isBottomRow) {
        graphIncluded->GetXaxis()->SetTitle("");
    }
}

void acceptanceStudy(double B, double C, int iterations) {
    for (int loop = 0; loop < iterations; ++loop) {
        std::cout << "Starting loop " << loop + 1 << " of " << iterations << std::endl;

        std::vector<double> phiVec;
        generateData(B, C, 1e5, phiVec); // Generate fresh data for each iteration

        // Create the example plot only during the first iteration
        TCanvas *masterCanvas = nullptr;
        masterCanvas = new TCanvas("masterCanvas", "Phi Distribution Fits", 1200, 800);
        setupCanvas(masterCanvas);

        int exclusionSteps[] = {0, 1, 2, 3, 4, 5};
        for (int i = 0; i < sizeof(exclusionSteps)/sizeof(exclusionSteps[0]); ++i) {
            int binsToExclude = exclusionSteps[i];
            plotForExclusion(phiVec, B, C, i + 1, binsToExclude, masterCanvas);
        }

        if (loop == 0) {
            std::ostringstream filename;
            filename << "output/acceptance_study_B=" << B << "_C=" << C << ".png";
            masterCanvas->SaveAs(filename.str().c_str());
        }
        delete masterCanvas;
    }
}


int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <B> <C>" << std::endl;
        return 1;
    }
    double B = atof(argv[1]);
    double C = atof(argv[2]);

    // Run the acceptance study n times
    acceptanceStudy(B, C, 10);
}
