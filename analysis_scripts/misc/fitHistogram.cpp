#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h" 
#include "TPaveText.h"
#include <iostream>

int main(int argc, char **argv) {

    TFile *file = TFile::Open(argv[1]);

    TH1F *hist = (TH1F*)file->Get("canvas");

    TF1 *fitFunc = new TF1("fitFunc", "pol1 + gaus(2)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Set parameter names
    fitFunc->SetParNames("b", "m", "A", "mu", "sigma");
    // Initial parameter guesses
    fitFunc->SetParameters(0, 0.1, hist->GetMaximum(), 0.94, 0.1); 
    fitFunc->SetParLimits(4, 0.0000, 1e10); // positive sigma
    hist->Fit(fitFunc, "R");
    hist->SetStats(kFALSE);

    TCanvas *c1 = new TCanvas("c1", "Histogram Fit", 800, 600);
    hist->Draw();
    fitFunc->Draw("same");

    // Create a TPaveText to display fitted parameters
    TPaveText *pt = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC"); // NDC coordinates
    pt->SetBorderSize(1);
    pt->SetFillStyle(1001);
    pt->SetTextAlign(12);
    pt->AddText(Form("b = %.3f #pm %.3f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
    pt->AddText(Form("m = %.3f #pm %.3f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
    pt->AddText(Form("A = %.3f #pm %.3f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));
    pt->AddText(Form("mu = %.3f #pm %.3f", fitFunc->GetParameter(3), fitFunc->GetParError(3)));
    pt->AddText(Form("sigma = %.3f #pm %.3f", fitFunc->GetParameter(4), fitFunc->GetParError(4)));
    pt->Draw();

    // Save the canvas as a PNG image
    c1->SaveAs("/u/home/thayward/fit_distro.png");

    // Cleanup
    delete c1;
    file->Close();
    delete file;

    return 0;
}
