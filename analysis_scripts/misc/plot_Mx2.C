#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLatex.h>

void plot_Mx2() {
    // Open the ROOT file
    TFile *file = TFile::Open("/Users/tbhayward/Desktop/updated_script_test/rga_fa18_inb_epX.root");
    if (!file || file->IsZombie()) {
        printf("Error: Could not open file.\n");
        return;
    }

    // Load the tree
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        printf("Error: Could not load tree.\n");
        file->Close();
        return;
    }

    // Create a histogram for Mx2 with 150 bins between -0.4 and 0.4
    TH1F *h_Mx2 = new TH1F("h_Mx2", "RGA Fa18 Inb", 150, -0.4, 0.4);

    // Draw the Mx2 branch into the histogram
    tree->Draw("Mx2 >> h_Mx2");

    // Set axis labels
    h_Mx2->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h_Mx2->GetYaxis()->SetTitle("Counts");

    // Set histogram style to data points with vertical error bars
    h_Mx2->SetMarkerStyle(20); // Circle marker
    h_Mx2->SetMarkerSize(0.7); // Slightly smaller marker size
    h_Mx2->SetMarkerColor(kBlack); // Marker color black
    h_Mx2->SetLineColor(kBlack);   // Line color black

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Mx^{2} Plot", 800, 600);

    // Remove the stat box
    gStyle->SetOptStat(0);

    // Set additional spacing on the left to avoid cutting off labels
    c1->SetLeftMargin(0.15);

    // Draw the histogram with error bars (E1 option)
    h_Mx2->Draw("E1");

    // Define a fitting function: quadratic + single gaussian
    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[4])**2)", -0.4, 0.4);

    // Set initial parameters for the fit
    fitFunc->SetParameters(100, 0, 50, 0, 0.1); // Initial guesses for params

    // Set parameter names for clarity
    fitFunc->SetParNames("const", "linear", "Gaussian_amp", "Gaussian_mean", "Gaussian_sigma");

    // Perform the fit
    h_Mx2->Fit("fitFunc");

    // Extract fit parameters for Gaussian
    double mu = fitFunc->GetParameter(3);       // Mean of Gaussian
    double sigma = fitFunc->GetParameter(4);    // Sigma of Gaussian

    // Create a LaTeX object to display the mu and sigma values
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC(); // Use normalized coordinates
    latex.DrawLatex(0.4, 0.3, Form("#mu = %.3f (GeV^{2}), #sigma = %.3f (GeV^{2})", mu, sigma));

    // Save the plot
    c1->SaveAs("/Users/tbhayward/Desktop/rga_Mx2_plot.png");

    // Clean up
    file->Close();
    delete c1;
    delete h_Mx2;
}