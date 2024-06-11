#include <iostream>
#include <string>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TLorentzVector.h>

void dilution_factors_rho(const char* nh3_file, const char* c_file) {
    // Open the ROOT files
    TFile *nh3 = TFile::Open(nh3_file);
    TFile *carbon = TFile::Open(c_file);
    if (!nh3 || nh3->IsZombie() || !carbon || carbon->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    // Get the PhysicsEvents trees
    TTree *tree_nh3;
    TTree *tree_carbon;
    nh3->GetObject("PhysicsEvents", tree_nh3);
    carbon->GetObject("PhysicsEvents", tree_carbon);
    if (!tree_nh3 || !tree_carbon) {
        std::cerr << "Error: PhysicsEvents tree not found!" << std::endl;
        nh3->Close();
        carbon->Close();
        return;
    }

    // Define the variables for trimming
    Double_t Mx23, z23, p1_p, p1_phi, p1_theta, p2_p, p2_theta, p2_phi, p3_theta;
    tree_nh3->SetBranchAddress("Mx23", &Mx23);
    tree_nh3->SetBranchAddress("z23", &z23);
    tree_nh3->SetBranchAddress("p1_p", &p1_p);
    tree_nh3->SetBranchAddress("p1_phi", &p1_phi);
    tree_nh3->SetBranchAddress("p1_theta", &p1_theta);
    tree_nh3->SetBranchAddress("p2_p", &p2_p);
    tree_nh3->SetBranchAddress("p2_theta", &p2_theta);
    tree_nh3->SetBranchAddress("p2_phi", &p2_phi);
    tree_nh3->SetBranchAddress("p3_theta", &p3_theta);

    // Create new trees with cuts applied
    TTree *new_tree_nh3 = tree_nh3->CloneTree(0);
    TTree *new_tree_carbon = tree_carbon->CloneTree(0);

    Long64_t nentries = tree_nh3->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree_nh3->GetEntry(i);

        // Apply the cuts
        if (fabs(Mx23 - 0.95) < 0.06 && z23 > 0.9) {
            // Calculate the calculated p3_theta
            TLorentzVector p1, p2, p3;
            p1.SetMagThetaPhi(p1_p, p1_theta, p1_phi); // Use provided theta and phi directly
            p2.SetMagThetaPhi(p2_p, p2_theta, p2_phi); // Use provided theta and phi directly

            // Calculate the missing four-momentum
            TLorentzVector system = p1 + p2;
            p3 = -system; // Missing momentum

            Double_t calculated_p3_theta = p3.Theta();
            std::cout << p3_theta << " " << calculated_p3_theta << " " << (p3_theta - calculated_p3_theta) << std::endl;

            // Apply the cut on p3_theta
            if (fabs(p3_theta - calculated_p3_theta) < 0.05) { // Adjusted tolerance to 0.05
                new_tree_nh3->Fill();
            }
        }
    }

    // Repeat the same process for carbon tree
    tree_carbon->SetBranchAddress("Mx23", &Mx23);
    tree_carbon->SetBranchAddress("z23", &z23);
    tree_carbon->SetBranchAddress("p1_p", &p1_p);
    tree_carbon->SetBranchAddress("p1_phi", &p1_phi);
    tree_carbon->SetBranchAddress("p1_theta", &p1_theta);
    tree_carbon->SetBranchAddress("p2_p", &p2_p);
    tree_carbon->SetBranchAddress("p2_theta", &p2_theta);
    tree_carbon->SetBranchAddress("p2_phi", &p2_phi);
    tree_carbon->SetBranchAddress("p3_theta", &p3_theta);

    nentries = tree_carbon->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree_carbon->GetEntry(i);

        // Apply the cuts
        if (fabs(Mx23 - 0.95) < 0.06 && z23 > 0.9) {
            // Calculate the calculated p3_theta
            TLorentzVector p1, p2, p3;
            p1.SetMagThetaPhi(p1_p, p1_theta, p1_phi); // Use provided theta and phi directly
            p2.SetMagThetaPhi(p2_p, p2_theta, p2_phi); // Use provided theta and phi directly

            // Calculate the missing four-momentum
            TLorentzVector system = p1 + p2;
            p3 = -system; // Missing momentum

            Double_t calculated_p3_theta = p3.Theta();
            std::cout << p3_theta << " " << calculated_p3_theta << " " << (p3_theta - calculated_p3_theta) << std::endl;

            // Apply the cut on p3_theta
            if (fabs(p3_theta - calculated_p3_theta) < 0.05) { // Adjusted tolerance to 0.05
                new_tree_carbon->Fill();
            }
        }
    }

    // Use the new trees for further analysis
    tree_nh3 = new_tree_nh3;
    tree_carbon = new_tree_carbon;

    // Create histograms for xF2
    // TH1D *h_xF2_nh3 = new TH1D("h_xF2_nh3", "x_{F2} Distribution; x_{F2}; Counts", 100, -2.5, 1);
    // TH1D *h_xF2_carbon = new TH1D("h_xF2_carbon", "x_{F2} Distribution; x_{F2}; Counts", 100, -2.5, 1);
    TH1D *h_xF2_nh3 = new TH1D("h_xF2_nh3", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -1, 2);
    TH1D *h_xF2_carbon = new TH1D("h_xF2_carbon", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -1, 2);

    // Fill the histograms
    // tree_nh3->Draw("xF2>>h_xF2_nh3");
    // tree_carbon->Draw("xF2>>h_xF2_carbon");
    tree_nh3->Draw("Mx>>h_xF2_nh3");
    tree_carbon->Draw("Mx>>h_xF2_carbon");
    // Fill the histograms with cuts
    // tree_nh3->Draw("Mx>>h_xF2_nh3", "helicity > 0 && target_pol < 0");
    // tree_carbon->Draw("Mx>>h_xF2_carbon", "helicity > 0");
    // Scale the carbon histogram by 1/2
    // h_xF2_carbon->Scale(0.5);

    // Create canvas and divide it into four panels
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1200, 1200);
    c1->Divide(2, 3);

    // First panel: plot xF2 histograms
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy(); // Log scale to better see differences
    h_xF2_nh3->SetLineColor(kBlue);
    h_xF2_carbon->SetLineColor(kRed);
    h_xF2_nh3->Draw();
    h_xF2_carbon->Draw("SAME");

    // Add legend
    TLegend *leg = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg->AddEntry(h_xF2_nh3, "NH_{3}", "l");
    leg->AddEntry(h_xF2_carbon, "C", "l");
    leg->Draw();

    // Remove statboxes
    h_xF2_nh3->SetStats(0);
    h_xF2_carbon->SetStats(0);

    // Second panel: ratio of NH3 to Carbon counts
    c1->cd(2);
    gPad->SetLeftMargin(0.15);
    TGraphErrors *gr_ratio = new TGraphErrors();
    for (int i = 1; i <= h_xF2_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_xF2_nh3->GetBinContent(i);
        double c_counts = h_xF2_carbon->GetBinContent(i);
        if (c_counts > 0) {
            double ratio = nh3_counts / c_counts;
            double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
            gr_ratio->SetPoint(i - 1, h_xF2_nh3->GetBinCenter(i), ratio);
            gr_ratio->SetPointError(i - 1, 0, error);
        }
    }
    // gr_ratio->SetTitle("NH_{3} to Carbon Ratio; x_{F2}; Ratio");
    gr_ratio->SetTitle("NH_{3} to Carbon Ratio; M_{x} (GeV); Ratio");
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->Draw("AP");

    // Fit the data from -2.5 to -1 to a constant
    TF1 *fit_const = new TF1("fit_const", "[0]", -1, -0.5);
    gr_ratio->Fit(fit_const, "R");
    fit_const->SetLineColor(kRed);
    fit_const->Draw("SAME");

    // Add fit constant value and uncertainty
    double fit_value = fit_const->GetParameter(0);
    double fit_error = fit_const->GetParError(0);
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.85, Form("Fit Const, s = %.3f #pm %.3f", fit_value, fit_error));

    // Third panel: Mh23 histograms scaled by the fit constant
    c1->cd(3);
    gPad->SetLeftMargin(0.15);
    TH1D *h_Mh23_nh3 = new TH1D("h_Mh23_nh3", "P_{1T}P_{2T} Distribution; P_{1T}P_{2T} (GeV); Counts", 5, 0, 0.5);
    TH1D *h_Mh23_carbon = new TH1D("h_Mh23_carbon", "P_{1T}P_{2T} Distribution; P_{1T}P_{2T} (GeV); Counts", 5, 0, 0.5);
    tree_nh3->Draw("Mh23>>h_Mh23_nh3"); 
    tree_carbon->Draw("Mh23>>h_Mh23_carbon");

    double scale_factor = fit_const->GetParameter(0);
    h_Mh23_carbon->Scale(scale_factor);

    h_Mh23_nh3->SetLineColor(kBlue);
    h_Mh23_carbon->SetLineColor(kRed);
    h_Mh23_nh3->Draw();
    h_Mh23_carbon->Draw("SAME");

    // Add legend
    TLegend *leg_Mh23 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg_Mh23->SetTextSize(0.04);
    leg_Mh23->AddEntry(h_Mh23_nh3, "NH_{3}", "l");
    leg_Mh23->AddEntry(h_Mh23_carbon, "s*C", "l");
    leg_Mh23->Draw();

    // Remove statboxes
    h_Mh23_nh3->SetStats(0);
    h_Mh23_carbon->SetStats(0);

    // Fourth panel: (NH3 - Carbon) / NH3 with fit to a third-degree polynomial
    c1->cd(4);
    gPad->SetLeftMargin(0.15);
    TGraphErrors *gr_dilution = new TGraphErrors();
    for (int i = 1; i <= h_Mh23_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_Mh23_nh3->GetBinContent(i);
        double c_counts = h_Mh23_carbon->GetBinContent(i);
        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;
            double error = std::sqrt((c_counts / nh3_counts) * (c_counts / nh3_counts) / nh3_counts + c_counts / (nh3_counts * nh3_counts));
            gr_dilution->SetPoint(i - 1, h_Mh23_nh3->GetBinCenter(i), dilution);
            gr_dilution->SetPointError(i - 1, 0, error);
        }
    }
    gr_dilution->SetTitle("Dilution Factor; M_{#pi#pi} (GeV); (NH3 - Carbon) / NH3");
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->Draw("AP");

    // Fit to a third-degree polynomial
    TF1 *fit_poly = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2", 0, 0.5);
    gr_dilution->Fit(fit_poly, "R");
    fit_poly->SetLineColor(kRed);
    fit_poly->Draw("SAME");

    // Add fit parameters box
    TPaveText *pt = new TPaveText(0.7, 0.7, 0.9, 0.9, "brNDC");
    pt->SetBorderSize(1);
    pt->SetFillStyle(1001); // Solid fill style
    pt->SetFillColor(kWhite); // White background
    pt->AddText(Form("p0 = %.3f", fit_poly->GetParameter(0)));
    pt->AddText(Form("p1 = %.3f", fit_poly->GetParameter(1)));
    pt->AddText(Form("p2 = %.3f", fit_poly->GetParameter(2)));
    pt->Draw();

    // Fifth panel: xB histograms scaled by the fit constant
    c1->cd(5);
    gPad->SetLeftMargin(0.15);
    TH1D *h_xB_nh3 = new TH1D("h_xB_nh3", "x_{B} Distribution; x_{B}; Counts", 50, 0.08, 0.58);
    TH1D *h_xB_carbon = new TH1D("h_xB_carbon", "x_{B} Distribution; x_{B}; Counts", 50, 0.08, 0.58);
    tree_nh3->Draw("x>>h_xB_nh3");
    tree_carbon->Draw("x>>h_xB_carbon");

    h_xB_carbon->Scale(scale_factor);

    h_xB_nh3->SetLineColor(kBlue);
    h_xB_carbon->SetLineColor(kRed);
    h_xB_nh3->Draw();
    h_xB_carbon->Draw("SAME");

    // Add legend
    TLegend *leg_xB = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg_xB->SetTextSize(0.04);
    leg_xB->AddEntry(h_xB_nh3, "NH_{3}", "l");
    leg_xB->AddEntry(h_xB_carbon, "s*C", "l");
    leg_xB->Draw();

    // Remove statboxes
    h_xB_nh3->SetStats(0);
    h_xB_carbon->SetStats(0);

    // Sixth panel: (NH3 - Carbon) / NH3 with fit to a third-degree polynomial for xB
    c1->cd(6);
    gPad->SetLeftMargin(0.15);
    TGraphErrors *gr_dilution_xB = new TGraphErrors();
    for (int i = 1; i <= h_xB_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_xB_nh3->GetBinContent(i);
        double c_counts = h_xB_carbon->GetBinContent(i);
        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;
            double error = std::sqrt((c_counts / nh3_counts) * (c_counts / nh3_counts) / nh3_counts + c_counts / (nh3_counts * nh3_counts));
            gr_dilution_xB->SetPoint(i - 1, h_xB_nh3->GetBinCenter(i), dilution);
            gr_dilution_xB->SetPointError(i - 1, 0, error);
        }
    }
    gr_dilution_xB->SetTitle("Dilution Factor; x_{B}; (NH3 - Carbon) / NH3");
    gr_dilution_xB->SetMarkerStyle(20);
    gr_dilution_xB->Draw("AP");

    std::vector<std::string> output;
    output.push_back("{");

    for (int i = 1; i <= h_xB_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_xB_nh3->GetBinContent(i);
        double c_counts = h_xB_carbon->GetBinContent(i);
        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;
            double error = std::sqrt((c_counts / nh3_counts) * (c_counts / nh3_counts) / nh3_counts + c_counts / (nh3_counts * nh3_counts));
            std::ostringstream ss;
            ss << "{" << h_xB_nh3->GetBinCenter(i) << ", " << dilution << ", " << error << "}";
            output.push_back(ss.str());
        }
    }

    output.push_back("}");

    for (size_t i = 0; i < output.size(); ++i) {
        std::cout << output[i];
        if (i != output.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;

    // Fit to a third-degree polynomial
    TF1 *fit_poly_xB = new TF1("fit_poly_xB", "[0] + [1]*x + [2]*x^2", 0.08, 0.58);
    gr_dilution_xB->Fit(fit_poly_xB, "R");
    fit_poly_xB->SetLineColor(kRed);
    fit_poly_xB->Draw("SAME");

    // Add fit parameters box
    TPaveText *pt_xB = new TPaveText(0.7, 0.7, 0.9, 0.9, "brNDC");
    pt_xB->SetBorderSize(1);
    pt_xB->SetFillStyle(1001); // Solid fill style
    pt_xB->SetFillColor(kWhite); // White background
    pt_xB->AddText(Form("p0 = %.3f", fit_poly_xB->GetParameter(0)));
    pt_xB->AddText(Form("p1 = %.3f", fit_poly_xB->GetParameter(1)));
    pt_xB->AddText(Form("p2 = %.3f", fit_poly_xB->GetParameter(2)));
    pt_xB->Draw();

    // Save the canvas
    c1->SaveAs("dilution_factors.pdf");

    // Clean up
    nh3->Close();
    carbon->Close();
    delete h_xF2_nh3;
    delete h_xF2_carbon;
    delete h_Mh23_nh3;
    delete h_Mh23_carbon;
    delete gr_ratio;
    delete fit_const;
    delete gr_dilution;
    delete fit_poly;
    delete pt;
    delete c1;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <NH3 ROOT file> <Carbon ROOT file>" << std::endl;
        return 1;
    }

    dilution_factors_rho(argv[1], argv[2]);
    return 0;
}