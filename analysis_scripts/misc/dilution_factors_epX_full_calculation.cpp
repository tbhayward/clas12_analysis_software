#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TF1.h>
#include <iostream>

void plot_dilution_factor(const char* variable_name, const char* x_title, double x_min, double x_max, int n_bins, TCanvas* canvas, int pad) {
    // Create a dummy histogram for the X axis range
    TH1D *h_dummy = new TH1D(Form("h_%s_dummy", variable_name), "", n_bins, x_min, x_max);
    
    // Placeholder for the dilution factor
    TGraphErrors *gr_dilution = new TGraphErrors();
    for (int i = 1; i <= n_bins; ++i) {
        double bin_center = h_dummy->GetBinCenter(i);
        gr_dilution->SetPoint(i - 1, bin_center, 0.175); // Placeholder value
        gr_dilution->SetPointError(i - 1, 0, 0.01); // Placeholder error
    }
    delete h_dummy;

    // Draw on the canvas
    canvas->cd(pad);
    gPad->SetLeftMargin(0.15);
    gr_dilution->SetTitle(Form(";%s;D_{f}", x_title));
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->Draw("AP");
    gr_dilution->GetXaxis()->SetRangeUser(x_min, x_max);
    gr_dilution->GetYaxis()->SetRangeUser(0.1, 0.3);

    // Fit to a third-degree polynomial (optional)
    TF1 *fit_poly = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2", x_min, x_max);
    gr_dilution->Fit(fit_poly, "RQ");
    fit_poly->SetLineColor(kRed);
    fit_poly->Draw("SAME");

    // Retrieve and display chi2/ndf in the top left
    double chi2 = fit_poly->GetChisquare();
    int ndf = fit_poly->GetNDF();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2 / ndf));

    // Optional: Add fit parameters box (this is kept simple here)
    TPaveText *pt = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
    pt->SetBorderSize(1);
    pt->SetFillStyle(1001); // Solid fill style
    pt->SetFillColor(kWhite); // White background
    pt->AddText(Form("Placeholder fit"));
    pt->Draw();
}

void one_dimensional(TFile* nh3, TFile* c, TFile* ch, TFile* he, TFile* empty) {
    // Create a canvas and divide it into 1 row and 3 columns
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1600, 600);
    c1->Divide(3, 1);

    // Plot for x-Bjorken
    plot_dilution_factor("x", "x_{B} (GeV)", 0.06, 0.6, 50, c1, 1);

    // Plot for transverse momentum
    plot_dilution_factor("pT", "P_{T} (GeV)", 0, 1.0, 50, c1, 2);

    // Plot for x-Feynman
    plot_dilution_factor("xF", "x_{F} (GeV)", -0.8, 0.5, 50, c1, 3);

    // Save the canvas as a PNG file
    c1->SaveAs("output/one_dimensional_placeholder.png");

    // Clean up
    delete c1;
}


int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <NH3 ROOT file> <Carbon ROOT file> <CH2 ROOT file> <Helium ROOT file> <Empty ROOT file>" << std::endl;
        return 1;
    }

    // Open the ROOT files
    TFile *nh3 = TFile::Open(argv[1]);
    TFile *c = TFile::Open(argv[2]);
    TFile *ch = TFile::Open(argv[3]);
    TFile *he = TFile::Open(argv[4]);
    TFile *empty = TFile::Open(argv[5]);

    // Check if files opened successfully
    if (!nh3 || nh3->IsZombie() || !c || c->IsZombie() || !ch || ch->IsZombie() || !he || he->IsZombie() || !empty || empty->IsZombie()) {
        std::cerr << "Error opening one or more files!" << std::endl;
        if (nh3) nh3->Close();
        if (c) c->Close();
        if (ch) ch->Close();
        if (he) he->Close();
        if (empty) empty->Close();
        return 2;
    }

    // Call the one-dimensional function
    one_dimensional(nh3, c, ch, he, empty);

    // Safely close the ROOT files
    nh3->Close();
    c->Close();
    ch->Close();
    he->Close();
    empty->Close();

    return 0;
}