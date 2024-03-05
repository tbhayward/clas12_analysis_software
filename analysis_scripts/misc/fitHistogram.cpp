#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>

int main(int argc, char **argv) {

    TFile *file = TFile::Open(argv[1]);

    TH1F *hist = (TH1F*)file->Get("canvas");

    TF1 *fitFunc = new TF1("fitFunc", "pol1 + gaus(2)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Set parameter names
    fitFunc->SetParNames("b", "m", "A", "mu", "sigma");
    // Initial parameter guesses
    fitFunc->SetParameters(0, 0.1, hist->GetMaximum(), 0.94, 0.1); 
    // fitFunc->SetParLimits(2, 0.0001, 1e10); 
    hist->Fit(fitFunc, "R");

    TCanvas *c1 = new TCanvas("c1", "Histogram Fit", 800, 600);
    hist->Draw();
    fitFunc->Draw("same");

    // Save the canvas as a PNG image
    c1->SaveAs("/u/home/thayward/fit_distro.png");

    // Cleanup
    delete c1;
    file->Close();
    delete file;

    return 0;
}
