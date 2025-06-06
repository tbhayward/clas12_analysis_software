#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMath.h>
#include <iostream>
#include <algorithm>

void plot_dvcs_energy_loss_validation(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* file4,
    const char* titleSuffix) {
    const int nFiles = 4;
    const char* files[nFiles] = { file1, file2, file3, file4 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "Timothy's",
        "Krishna's",
        "Mariana's"
    };

    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = (TTree*)f[i]->Get("PhysicsEvents");
    }

    Double_t p1_theta[nFiles], Mx2_1[nFiles];
    Double_t eta2[nFiles], t1[nFiles], theta_gamma_gamma[nFiles];
    Double_t Emiss2[nFiles], pTmiss[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("p1_theta",          &p1_theta[i]);
        tree[i]->SetBranchAddress("Mx2_1",             &Mx2_1[i]);
        tree[i]->SetBranchAddress("eta2",              &eta2[i]);
        tree[i]->SetBranchAddress("t1",                &t1[i]);
        tree[i]->SetBranchAddress("theta_gamma_gamma", &theta_gamma_gamma[i]);
        tree[i]->SetBranchAddress("Emiss2",            &Emiss2[i]);
        tree[i]->SetBranchAddress("pTmiss",            &pTmiss[i]);
    }

    TCanvas* c1 = new TCanvas("c1", "DVCS Energy Loss Validation", 1200, 900);
    c1->Divide(4, 3);

    const int    nBins        = 10;
    Double_t     thetaBins[nBins+1] = {5,15,20,25,30,35,40,45,50,60,100};
    const int    nbMx2_hi     = 35;
    const int    nbMx2_lo     = nbMx2_hi/2;
    const Double_t mx2_min    = -0.3, mx2_max = +0.3;

    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];
    Double_t theta_sum[nBins]     = {0};
    Int_t    theta_count[nBins]   = {0};
    Double_t theta_mean[nBins]    = {0};
    Double_t mx2_sum[nFiles][nBins]   = {{0}};
    Int_t    mx2_count[nFiles][nBins] = {{0}};
    Double_t mx2_mean[nFiles][nBins]  = {{0}};

    // create histograms
    for (int i = 0; i < nFiles; ++i) {
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [5,65] %s", titleSuffix),
            nbMx2_hi, mx2_min, mx2_max
        );
        for (int b = 0; b < nBins; ++b) {
            int nb = (b < nBins/2 ? nbMx2_lo : nbMx2_hi);
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.0f,%.0f] %s",
                     thetaBins[b], thetaBins[b+1],
                     titleSuffix),
                nb, mx2_min, mx2_max
            );
        }
    }

    // fill & accumulate
    for (int i = 0; i < nFiles; ++i) {
        Long64_t nEntries = tree[i]->GetEntries();
        for (Long64_t j = 0; j < nEntries; ++j) {
            tree[i]->GetEntry(j);
            Double_t thetaDeg = p1_theta[i] * 180.0 / TMath::Pi();
            if (thetaDeg >= 5 && thetaDeg < 65 &&
                eta2[i] < 0 &&
                t1[i]   > -2 &&
                theta_gamma_gamma[i] < 0.6 &&
                Emiss2[i] < 0.5 &&
                pTmiss[i] < 0.125) {

                h[i][0]->Fill(Mx2_1[i]);
                for (int b = 0; b < nBins; ++b) {
                    if (thetaDeg >= thetaBins[b] && thetaDeg < thetaBins[b+1]) {
                        h[i][b+1]->Fill(Mx2_1[i]);
                        theta_sum[b]   += thetaDeg;
                        theta_count[b] += 1;
                        mx2_sum[i][b]   += Mx2_1[i];
                        mx2_count[i][b] += 1;
                    }
                }
            }
        }
    }

    // compute means
    for (int b = 0; b < nBins; ++b) {
        theta_mean[b] = (theta_count[b] > 0
                         ? theta_sum[b] / theta_count[b]
                         : NAN);
        for (int i = 0; i < nFiles; ++i) {
            mx2_mean[i][b] = (mx2_count[i][b] > 0
                              ? mx2_sum[i][b] / mx2_count[i][b]
                              : NAN);
        }
    }

    // integrated pad
    c1->cd(1);
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    Double_t globalMax = 0;
    for (int i = 0; i < nFiles; ++i) {
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    }
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        h[i][0]->SetStats(0);
        if (i == 0) h[i][0]->Draw("E");
        else        h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(Form("fitInt%d", i),
                            "gaus(0)+pol1(3)",
                            mx2_min, mx2_max);
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(),  // A init
            0.0,                          // mu init
            0.1                           // sigma init
        );
        fitInt[i]->SetParLimits(1, -0.15, 0.15);
        fitInt[i]->SetParLimits(2,  0.0, 0.3);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }
    TLegend* legInt = new TLegend(0.25, 0.75, 0.9, 0.9);
    legInt->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        legInt->AddEntry(
            h[i][0],
            Form("%s: #mu=%.3f, #sigma=%.3f",
                 corrLabels[i],
                 fitInt[i]->GetParameter(1),
                 fitInt[i]->GetParameter(2)),
            "lep"
        );
    }
    legInt->Draw();
    h[0][0]->GetXaxis()->SetTitle("M_{xp}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // theta‐binned pads
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1);
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        Double_t binMax = 0;
        for (int i = 0; i < nFiles; ++i) {
            binMax = std::max(binMax, h[i][b]->GetMaximum());
        }
        for (int i = 0; i < nFiles; ++i) {
            h[i][b]->SetMaximum(1.7 * binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20 + i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack + i);
            h[i][b]->SetStats(0);
            if (i == 0) h[i][b]->Draw("E");
            else        h[i][b]->Draw("E SAME");

            TF1* fbin = new TF1(Form("fitBin%d_%d", i, b),
                                "gaus(0)+pol1(3)",
                                mx2_min, mx2_max);
            fbin->SetParameters(
                0.8 * h[i][b]->GetMaximum(),
                0.0,
                0.1
            );
            fbin->SetParLimits(1, -0.15, 0.15);
            fbin->SetParLimits(2,  0.0, 0.3);
            fbin->SetLineColor(kBlack + i);
            fbin->SetLineWidth(1);
            h[i][b]->Fit(fbin, "Q");
            fbin->Draw("SAME");
            mu[i][b-1]    = fbin->GetParameter(1);
            sigma[i][b-1] = fbin->GetParameter(2);
        }
        TLegend* legB = new TLegend(0.25, 0.75, 0.9, 0.9);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("M_{xp}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // final pad: mu vs mean theta
    c1->cd(12);
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);
    TGraph* gr[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        std::vector<double> xs, ys;
        for (int b = 0; b < nBins; ++b) {
            if (theta_count[b] > 0) {
                xs.push_back(theta_mean[b]);
                ys.push_back(mu[i][b]);
            }
        }
        gr[i] = new TGraph(xs.size(), xs.data(), ys.data());
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack + i);
        if (i == 0) gr[i]->Draw("AP");
        else        gr[i]->Draw("P SAME");
    }
    TLine* line = new TLine(0, 0, 90, 0);
    line->SetLineColor(kGray);
    line->SetLineStyle(2);
    line->Draw("SAME");
    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(0, 90);
    gr[0]->GetYaxis()->SetRangeUser(-0.1, 0.1);

    TLegend* leg12 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        leg12->AddEntry(gr[i], corrLabels[i], "lep");
    }
    leg12->Draw();

    TString outname = TString::Format("output/dvcs_%s_energy_loss_validation.pdf", titleSuffix);
    c1->SaveAs(outname);

    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        f[i]->Close();
        delete f[i];
    }
}


void plot_eppi0_energy_loss_validation(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* file4,
    const char* titleSuffix) {
    const int nFiles = 4;
    const char* files[nFiles] = { file1, file2, file3, file4 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "Timothy's",
        "Krishna's",
        "Mariana's"
    };

    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = (TTree*)f[i]->Get("PhysicsEvents");
    }

    Double_t p1_theta[nFiles], Mx2_1[nFiles];
    Double_t eta2[nFiles], t1[nFiles], theta_gamma_gamma[nFiles];
    Double_t Emiss2[nFiles], pTmiss[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("p1_theta",          &p1_theta[i]);
        tree[i]->SetBranchAddress("Mx2_1",             &Mx2_1[i]);
        tree[i]->SetBranchAddress("eta2",              &eta2[i]);
        tree[i]->SetBranchAddress("t1",                &t1[i]);
        tree[i]->SetBranchAddress("theta_pi0_pi0", &theta_gamma_gamma[i]);
        tree[i]->SetBranchAddress("Emiss2",            &Emiss2[i]);
        tree[i]->SetBranchAddress("pTmiss",            &pTmiss[i]);
    }

    TCanvas* c1 = new TCanvas("c1", "DVCS Energy Loss Validation", 1200, 900);
    c1->Divide(4, 3);

    const int    nBins        = 10;
    Double_t     thetaBins[nBins+1] = {5,15,20,25,30,35,40,45,50,60,100};
    const int    nbMx2_hi     = 35;
    const int    nbMx2_lo     = nbMx2_hi/2;
    const Double_t mx2_min    = -0.3, mx2_max = +0.3;

    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];
    Double_t theta_sum[nBins]     = {0};
    Int_t    theta_count[nBins]   = {0};
    Double_t theta_mean[nBins]    = {0};
    Double_t mx2_sum[nFiles][nBins]   = {{0}};
    Int_t    mx2_count[nFiles][nBins] = {{0}};
    Double_t mx2_mean[nFiles][nBins]  = {{0}};

    // create histograms
    for (int i = 0; i < nFiles; ++i) {
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [5,65] %s", titleSuffix),
            nbMx2_hi, mx2_min, mx2_max
        );
        for (int b = 0; b < nBins; ++b) {
            int nb = (b < nBins/2 ? nbMx2_lo : nbMx2_hi);
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.0f,%.0f] %s",
                     thetaBins[b], thetaBins[b+1],
                     titleSuffix),
                nb, mx2_min, mx2_max
            );
        }
    }

    // fill & accumulate
    for (int i = 0; i < nFiles; ++i) {
        Long64_t nEntries = tree[i]->GetEntries();
        for (Long64_t j = 0; j < nEntries; ++j) {
            tree[i]->GetEntry(j);
            Double_t thetaDeg = p1_theta[i] * 180.0 / TMath::Pi();
            if (thetaDeg >= 5 && thetaDeg < 65 &&
                eta2[i] < 0 &&
                t1[i]   > -2 &&
                theta_gamma_gamma[i] < 0.6 &&
                Emiss2[i] < 0.5 &&
                pTmiss[i] < 0.125) {

                h[i][0]->Fill(Mx2_1[i]);
                for (int b = 0; b < nBins; ++b) {
                    if (thetaDeg >= thetaBins[b] && thetaDeg < thetaBins[b+1]) {
                        h[i][b+1]->Fill(Mx2_1[i]);
                        theta_sum[b]   += thetaDeg;
                        theta_count[b] += 1;
                        mx2_sum[i][b]   += Mx2_1[i];
                        mx2_count[i][b] += 1;
                    }
                }
            }
        }
    }

    // compute means
    for (int b = 0; b < nBins; ++b) {
        theta_mean[b] = (theta_count[b] > 0
                         ? theta_sum[b] / theta_count[b]
                         : NAN);
        for (int i = 0; i < nFiles; ++i) {
            mx2_mean[i][b] = (mx2_count[i][b] > 0
                              ? mx2_sum[i][b] / mx2_count[i][b]
                              : NAN);
        }
    }

    // integrated pad
    c1->cd(1);
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    Double_t globalMax = 0;
    for (int i = 0; i < nFiles; ++i) {
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    }
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        h[i][0]->SetStats(0);
        if (i == 0) h[i][0]->Draw("E");
        else        h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(Form("fitInt%d", i),
                            "gaus(0)+pol1(3)",
                            mx2_min, mx2_max);
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(),  // A init
            0.13497*0.13497,                          // mu init
            0.1                           // sigma init
        );
        fitInt[i]->SetParLimits(1, -0.15, 0.2);
        fitInt[i]->SetParLimits(2,  0.0, 0.3);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }
    TLegend* legInt = new TLegend(0.25, 0.75, 0.9, 0.9);
    legInt->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        legInt->AddEntry(
            h[i][0],
            Form("%s: #mu=%.3f, #sigma=%.3f",
                 corrLabels[i],
                 fitInt[i]->GetParameter(1),
                 fitInt[i]->GetParameter(2)),
            "lep"
        );
    }
    legInt->Draw();
    h[0][0]->GetXaxis()->SetTitle("M_{xp}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // theta‐binned pads
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1);
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        Double_t binMax = 0;
        for (int i = 0; i < nFiles; ++i) {
            binMax = std::max(binMax, h[i][b]->GetMaximum());
        }
        for (int i = 0; i < nFiles; ++i) {
            h[i][b]->SetMaximum(1.7 * binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20 + i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack + i);
            h[i][b]->SetStats(0);
            if (i == 0) h[i][b]->Draw("E");
            else        h[i][b]->Draw("E SAME");

            TF1* fbin = new TF1(Form("fitBin%d_%d", i, b),
                                "gaus(0)+pol1(3)",
                                mx2_min, mx2_max);
            fbin->SetParameters(
                0.8 * h[i][b]->GetMaximum(),
                0.13497*0.13497,
                0.1
            );
            fbin->SetParLimits(1, -0.15, 0.2);
            fbin->SetParLimits(2,  0.0, 0.3);
            fbin->SetLineColor(kBlack + i);
            fbin->SetLineWidth(1);
            h[i][b]->Fit(fbin, "Q");
            fbin->Draw("SAME");
            mu[i][b-1]    = fbin->GetParameter(1);
            sigma[i][b-1] = fbin->GetParameter(2);
        }
        TLegend* legB = new TLegend(0.25, 0.75, 0.9, 0.9);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("M_{xp}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // final pad: mu vs mean theta
    c1->cd(12);
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);
    TGraph* gr[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        std::vector<double> xs, ys;
        for (int b = 0; b < nBins; ++b) {
            if (theta_count[b] > 0) {
                xs.push_back(theta_mean[b]);
                ys.push_back(mu[i][b]);
            }
        }
        gr[i] = new TGraph(xs.size(), xs.data(), ys.data());
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack + i);
        if (i == 0) gr[i]->Draw("AP");
        else        gr[i]->Draw("P SAME");
    }
    TLine* line = new TLine(0, 0.13497*0.13497, 90, 0.13497*0.13497);
    line->SetLineColor(kGray);
    line->SetLineStyle(2);
    line->Draw("SAME");
    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(0, 90);
    gr[0]->GetYaxis()->SetRangeUser(-0.1, 0.1);

    TLegend* leg12 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        leg12->AddEntry(gr[i], corrLabels[i], "lep");
    }
    leg12->Draw();

    TString outname = TString::Format("output/eppi0_%s_energy_loss_validation.pdf", titleSuffix);
    c1->SaveAs(outname);

    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        f[i]->Close();
        delete f[i];
    }
}

void plot_mx2_comparison_elastic(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* file4,
    const char* titleSuffix) {
    const int nFiles = 4;
    const char* files[nFiles] = { file1, file2, file3, file4 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "Timothy's",
        "Krishna's",
        "Mariana's"
    };

    // 1) Open files & trees
    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = f[i] ? (TTree*)f[i]->Get("PhysicsEvents") : nullptr;
    }

    // 2) Branch variables: p_theta and W
    Double_t p_theta[nFiles], W[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        if (!tree[i]) continue;
        tree[i]->SetBranchAddress("p_theta", &p_theta[i]);
        tree[i]->SetBranchAddress("Mx2",       &W[i]);
    }

    // 3) θ–bins: 10 evenly spaced from 15 to 70 deg
    const int nBins = 10;
    Double_t thetaBins[nBins+1];
    for (int b = 0; b <= nBins; ++b) {
        thetaBins[b] = 15.0 + b * (70.0-15.0)/nBins;
    }

    // 4) Canvas & grid
    TCanvas* c1 = new TCanvas("c1","Elastic W Comparison",1200,900);
    c1->Divide(4,3);

    // 5) Histogram parameters for W in [0.7,1.1]
    const int    nbW     = 35;
    const double W_min   = -0.04, W_max = 0.04;

    // 6) Allocate histograms and fit objects
    TH1D*   h[nFiles][nBins+1];
    TF1*    fitInt[nFiles];
    TF1*    fitBin[nFiles][nBins];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];

    // 7) Create histograms
    for (int i = 0; i < nFiles; ++i) {
        // integrated over all θ
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [15,70] %s", titleSuffix),
            nbW, W_min, W_max
        );
        // per-θ bins
        for (int b = 0; b < nBins; ++b) {
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.1f,%.1f] %s",
                     thetaBins[b], thetaBins[b+1],
                     titleSuffix),
                nbW, W_min, W_max
            );
        }
    }

    // 8) Fill histograms (no W cut now)
    for (int i = 0; i < nFiles; ++i) {
        if (!tree[i]) continue;
        Long64_t N = tree[i]->GetEntries();
        for (Long64_t ev = 0; ev < N; ++ev) {
            tree[i]->GetEntry(ev);
            double θ = p_theta[i] * 180.0 / TMath::Pi();
            // require θ in [15,70]
            if (θ < 15.0 || θ >= 70.0) continue;
            // integrated
            h[i][0]->Fill(W[i]);
            // per-θ
            for (int b = 0; b < nBins; ++b) {
                if (θ >= thetaBins[b] && θ < thetaBins[b+1]) {
                    h[i][b+1]->Fill(W[i]);
                }
            }
        }
    }

    // 9) Draw integrated pad with Gaussian+quadratic fit
    c1->cd(1);
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    double globalMax = 0;
    for (int i = 0; i < nFiles; ++i)
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        h[i][0]->SetStats(0);
        if (i == 0) h[i][0]->Draw("E");
        else        h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(
            Form("fitInt%d", i),
            "gaus(0)+pol2(3)",
            W_min, W_max
        );
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(),  // amplitude
            0.0,                        // μ init ≃ proton mass
            0.005,                         // σ init
            0, 0, 0                       // polynomial terms
        );
        fitInt[i]->SetParLimits(1, -0.02, 0.02);
        fitInt[i]->SetParLimits(2, 0.0,  0.1);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }
    {
        TLegend* leg = new TLegend(0.25,0.75,0.9,0.9);
        leg->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            leg->AddEntry(
                h[i][0],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     fitInt[i]->GetParameter(1),
                     fitInt[i]->GetParameter(2)),
                "lep"
            );
        }
        leg->Draw();
    }
    h[0][0]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // 10) Draw θ‐binned pads with fits
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1);
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        double binMax = 0;
        for (int i = 0; i < nFiles; ++i)
            binMax = std::max(binMax, h[i][b]->GetMaximum());
        for (int i = 0; i < nFiles; ++i) {
            h[i][b]->SetMaximum(1.7 * binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20 + i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack + i);
            h[i][b]->SetStats(0);
            if (i == 0) h[i][b]->Draw("E");
            else        h[i][b]->Draw("E SAME");

            fitBin[i][b-1] = new TF1(
                Form("fitBin%d_%d", i, b),
                "gaus(0)+pol2(3)",
                W_min, W_max
            );
            fitBin[i][b-1]->SetParameters(
                0.8 * h[i][b]->GetMaximum(),
                0.0,
                0.05,
                0, 0, 0
            );
            fitBin[i][b-1]->SetParLimits(1, -0.02, 0.02);
            fitBin[i][b-1]->SetParLimits(2, 0.0,  0.1);
            fitBin[i][b-1]->SetLineColor(kBlack + i);
            fitBin[i][b-1]->SetLineWidth(1);
            h[i][b]->Fit(fitBin[i][b-1], "Q");
            fitBin[i][b-1]->Draw("SAME");

            mu[i][b-1]    = fitBin[i][b-1]->GetParameter(1);
            sigma[i][b-1] = fitBin[i][b-1]->GetParameter(2);
        }
        TLegend* legB = new TLegend(0.25,0.75,0.9,0.9);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
        h[0][b]->GetXaxis()->SetRangeUser(W_min, W_max);
    }

    // 11) Save & clean up
    c1->SaveAs(Form("output/W_elastic_comparison_%s.pdf", titleSuffix));
    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        if (f[i]) { f[i]->Close(); delete f[i]; }
    }
}

void plot_two_pions(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* file4,
    const char* titleSuffix) {
    const int nFiles = 4;
    const char* files[nFiles] = { file1, file2, file3, file4 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "Timothy's",
        "Krishna's",
        "Mariana's"
    };

    // 1) Open files & trees
    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = f[i] ? static_cast<TTree*>(f[i]->Get("PhysicsEvents")) : nullptr;
        if (!tree[i]) {
            std::cerr << "[two_pions] ERROR opening PhysicsEvents in "
                      << files[i] << "\n";
            return;
        }
    }

    // 2) Branches: missing mass² and proton θ
    Double_t Mx2[nFiles], p1_theta[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("Mx2_1",   &Mx2[i]);
        tree[i]->SetBranchAddress("p1_theta",&p1_theta[i]);
    }

    // 3) θ‐bin edges: [5–10],[10–15],…,[50–100]
    const int nBins = 10;
    Double_t thetaBins[nBins+1] = {5,10,15,20,25,30,35,40,45,50,100};

    // 4) Canvas & grid
    TCanvas* c1 = new TCanvas("c1","Two-Pion M_{x}^{2} Validation",1200,900);
    c1->Divide(4,3);

    // 5) Histogram range
    const int    nbMx2   = 25;
    const double mx2_min = 0.35, mx2_max = 0.95;

    // 6) Storage
    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    TF1*     fitBin[nFiles][nBins];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];
    Double_t theta_sum[nBins]   = {0};
    Int_t    theta_count[nBins] = {0};
    Double_t theta_mean[nBins]  = {0};

    // 7) Create histograms
    for (int i = 0; i < nFiles; ++i) {
        // integrated over all θ
        h[i][0] = new TH1D(
            Form("h%d_int",i),
            Form("Integrated #theta [5,100] %s", titleSuffix),
            nbMx2, mx2_min, mx2_max
        );
        h[i][0]->SetStats(0);

        // per‐θ bins
        for (int b = 0; b < nBins; ++b) {
            h[i][b+1] = new TH1D(
                Form("h%d_%d",i,b),
                Form("#theta [%.0f,%.0f] %s", thetaBins[b], thetaBins[b+1], titleSuffix),
                nbMx2, mx2_min, mx2_max
            );
            h[i][b+1]->SetStats(0);
        }
    }

    // 8) Fill & accumulate θ
    for (int i = 0; i < nFiles; ++i) {
        Long64_t N = tree[i]->GetEntries();
        for (Long64_t ev=0; ev<N; ++ev) {
            tree[i]->GetEntry(ev);
            double θ = p1_theta[i]*180.0/TMath::Pi();
            h[i][0]->Fill(Mx2[i]);
            for (int b = 0; b < nBins; ++b) {
                if (θ >= thetaBins[b] && θ < thetaBins[b+1]) {
                    h[i][b+1]->Fill(Mx2[i]);
                    theta_sum[b]   += θ;
                    theta_count[b] += 1;
                }
            }
        }
    }

    // 9) Compute mean θ per bin
    for (int b = 0; b < nBins; ++b) {
        theta_mean[b] = (theta_count[b] > 0
                         ? theta_sum[b] / theta_count[b]
                         : NAN);
    }

    // 10) Draw & fit integrated
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    double globalMax = 0;
    for (int i = 0; i < nFiles; ++i)
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7*globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20+i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack+i);
        if (i==0) h[i][0]->Draw("E");
        else      h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(
            Form("fitInt%d",i),
            "gaus(0)+pol2(3)",
            mx2_min, mx2_max
        );
        // μ₀ ≃ m_{π±}² ≃ (0.1396)² ≃ 0.0195, σ₀ = 0.20
        fitInt[i]->SetParameters(
            0.8*h[i][0]->GetMaximum(),
            0.775*0.775, 0.05,
            0,0,0
        );
        fitInt[i]->SetParLimits(1, 0.4, 0.8);
        fitInt[i]->SetParLimits(2, 0.00, 0.30);
        fitInt[i]->SetLineColor(kBlack+i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i],"Q");
        fitInt[i]->Draw("SAME");
    }
    {
        TLegend* leg = new TLegend(0.25,0.75,0.9,0.9);
        leg->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            leg->AddEntry(
                h[i][0],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     fitInt[i]->GetParameter(1),
                     fitInt[i]->GetParameter(2)),
                "lep"
            );
        }
        leg->Draw();
    }
    h[0][0]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // 11) θ‐binned pads
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        double binMax = 0;
        for (int i = 0; i < nFiles; ++i)
            binMax = std::max(binMax, h[i][b]->GetMaximum());
        for (int i = 0; i < nFiles; ++i) {
            h[i][b]->SetMaximum(1.7*binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20+i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack+i);
            if (i==0) h[i][b]->Draw("E");
            else      h[i][b]->Draw("E SAME");

            fitBin[i][b-1] = new TF1(
                Form("fitBin%d_%d",i,b),
                "gaus(0)+pol2(3)",
                mx2_min, mx2_max
            );
            fitBin[i][b-1]->SetParameters(
                0.8*h[i][b]->GetMaximum(),
                0.775*0.775, 0.05,
                0,0,0
            );
            fitBin[i][b-1]->SetParLimits(1, 0.4, 0.8);
            fitBin[i][b-1]->SetParLimits(2, 0.00, 0.30);
            fitBin[i][b-1]->SetLineColor(kBlack+i);
            fitBin[i][b-1]->SetLineWidth(1);
            h[i][b]->Fit(fitBin[i][b-1],"Q");
            fitBin[i][b-1]->Draw("SAME");

            mu[i][b-1]    = fitBin[i][b-1]->GetParameter(1);
            sigma[i][b-1] = fitBin[i][b-1]->GetParameter(2);
        }
        TLegend* legB = new TLegend(0.25,0.75,0.9,0.9);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // 12) Final pad: μ vs. <θ>
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);
    TGraph* gr[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        std::vector<double> xs, ys;
        for (int b = 0; b < nBins; ++b) {
            if (theta_count[b] > 0) {
                xs.push_back(theta_mean[b]);
                ys.push_back(mu[i][b]);
            }
        }
        gr[i] = new TGraph(xs.size(), xs.data(), ys.data());
        gr[i]->SetMarkerStyle(20+i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack+i);
        if (i==0) gr[i]->Draw("AP");
        else      gr[i]->Draw("P SAME");
    }
    // dashed zero line
    TLine* line = new TLine(thetaBins[0],0.775*0.775,thetaBins[nBins],0.775*0.775);
    line->SetLineColor(kGray);
    line->SetLineStyle(2);
    line->Draw("SAME");

    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(thetaBins[0], thetaBins[nBins]);
    gr[0]->GetYaxis()->SetRangeUser(0.5, 0.7);

    TLegend* leg12 = new TLegend(0.6,0.75,0.9,0.9);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        leg12->AddEntry(gr[i], corrLabels[i], "lep");
    }
    leg12->Draw();

    // 13) Save & cleanup
    c1->SaveAs(Form("output/two_pions_%s.pdf", titleSuffix));
    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        f[i]->Close();
        delete f[i];
    }
}

void plot_dvcs_sebastian_energy_loss_validation(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* titleSuffix) {
    // 0) Setup
    const int nFiles = 3;
    const char* files[nFiles] = { file1, file2, file3 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "SA e^{-} + #gamma",
        "SA e^{-} + #gamma, TBH p"
    };

    // 1) Open files & trees
    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = f[i]
                  ? static_cast<TTree*>(f[i]->Get("PhysicsEvents"))
                  : nullptr;
        if (!tree[i]) {
            std::cerr << "[sebastian] ERROR opening PhysicsEvents in "
                      << files[i] << "\n";
            return;
        }
    }

    // 2) Branch variables
    Double_t p2_theta[nFiles], Mx2_1[nFiles];
    Double_t eta2[nFiles], t1[nFiles], theta_gamma_gamma[nFiles];
    Double_t Emiss2[nFiles], pTmiss[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("p2_theta",           &p2_theta[i]);
        tree[i]->SetBranchAddress("Mx2_1",              &Mx2_1[i]);
        tree[i]->SetBranchAddress("eta2",               &eta2[i]);
        tree[i]->SetBranchAddress("t1",                 &t1[i]);
        tree[i]->SetBranchAddress("theta_gamma_gamma",  &theta_gamma_gamma[i]);
        tree[i]->SetBranchAddress("Emiss2",             &Emiss2[i]);
        tree[i]->SetBranchAddress("pTmiss",             &pTmiss[i]);
    }

    // 3) Canvas & layout
    TCanvas* c1 = new TCanvas(
        "c1",
        "DVCS Energy Loss Validation (Sebastian)",
        1200, 900
    );
    c1->Divide(4, 3);

    // 4) θ bins: 5 → 32 deg, 10 equal bins
    const int    nBins = 10;
    Double_t     thetaBins[nBins+1];
    for (int b = 0; b <= nBins; ++b) {
        thetaBins[b] = 0.0 + b * (32.0 - 0.0) / nBins;
    }

    // 5) Histogram parameters
    const int    nbHi     = 35;
    const int    nbLo     = nbHi/2;  // 17
    const double mx2_min  = -0.3;
    const double mx2_max  = +0.3;

    // 6) Storage
    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];
    Double_t theta_sum[nBins]   = {0};
    Int_t    theta_count[nBins] = {0};
    Double_t theta_mean[nBins]  = {0};

    // 7) Create histograms (no stats box)
    for (int i = 0; i < nFiles; ++i) {
        // integrated
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [5,32] %s", titleSuffix),
            nbHi, mx2_min, mx2_max
        );
        h[i][0]->SetStats(false);

        // per-θ slices
        for (int b = 0; b < nBins; ++b) {
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.1f,%.1f] %s",
                     thetaBins[b], thetaBins[b+1], titleSuffix),
                (b < nBins/2 ? nbLo : nbHi),
                mx2_min, mx2_max
            );
            h[i][b+1]->SetStats(false);
        }
    }

    // 8) Fill & accumulate (with kinematic cuts)
    for (int i = 0; i < nFiles; ++i) {
        Long64_t N = tree[i]->GetEntries();
        for (Long64_t ev = 0; ev < N; ++ev) {
            tree[i]->GetEntry(ev);
            double θ = p2_theta[i] * 180.0 / TMath::Pi();
            if (
                // θ >= 5.0 && 
                θ < 32.0 &&
                eta2[i] <  0    &&
                t1[i]   > -2    &&
                theta_gamma_gamma[i] < 0.6 &&
                Emiss2[i] < 0.5 &&
                pTmiss[i] < 0.125)
            {
                h[i][0]->Fill(Mx2_1[i]);
                for (int b = 0; b < nBins; ++b) {
                    if (θ >= thetaBins[b] && θ < thetaBins[b+1]) {
                        h[i][b+1]->Fill(Mx2_1[i]);
                        theta_sum[b]   += θ;
                        theta_count[b] += 1;
                    }
                }
            }
        }
    }

    // 9) Compute <θ>
    for (int b = 0; b < nBins; ++b) {
        theta_mean[b] = (theta_count[b] > 0)
                       ? theta_sum[b] / theta_count[b]
                       : NAN;
    }

    // 10) Draw & fit integrated pad (1)
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    double globalMax = 0;
    for (int i = 0; i < nFiles; ++i)
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        if (i == 0) h[i][0]->Draw("E");
        else        h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(
            Form("fitInt%d", i),
            "gaus(0)+pol1(3)",
            mx2_min, mx2_max
        );
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(),  // A
            0.0,                          // μ init
            0.1                           // σ init
        );
        fitInt[i]->SetParLimits(1, -0.15, 0.15);
        fitInt[i]->SetParLimits(2,  0.0, 0.3);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }
    // integrated legend
    TLegend* legInt = new TLegend(0.20, 0.75, 0.95, 0.90);
    legInt->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        legInt->AddEntry(
            h[i][0],
            Form("%s: #mu=%.3f, #sigma=%.3f",
                 corrLabels[i],
                 fitInt[i]->GetParameter(1),
                 fitInt[i]->GetParameter(2)),
            "lep"
        );
    }
    legInt->Draw();
    h[0][0]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // 11) θ‐binned pads & fits (2–11)
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        double binMax = 0;
        for (int i = 0; i < nFiles; ++i)
            binMax = std::max(binMax, h[i][b]->GetMaximum());

        for (int i = 0; i < nFiles; ++i) {
            h[i][b]->SetMaximum(1.7 * binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20 + i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack + i);
            if (i == 0) h[i][b]->Draw("E");
            else        h[i][b]->Draw("E SAME");

            TF1* fbin = new TF1(
                Form("fitBin%d_%d", i, b),
                "gaus(0)+pol1(3)",
                mx2_min, mx2_max
            );
            fbin->SetParameters(
                0.8 * h[i][b]->GetMaximum(),
                0.0,
                0.1
            );
            fbin->SetParLimits(1, -0.15, 0.15);
            fbin->SetParLimits(2,  0.0, 0.3);
            fbin->SetLineColor(kBlack + i);
            fbin->SetLineWidth(1);
            h[i][b]->Fit(fbin, "Q");
            fbin->Draw("SAME");

            mu[i][b-1]    = fbin->GetParameter(1);
            sigma[i][b-1] = fbin->GetParameter(2);
        }
        // per-θ legend
        TLegend* legB = new TLegend(0.20, 0.75, 0.95, 0.90);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // 12) final pad: μ vs. <θ>
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);
    TGraph* gr[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        std::vector<double> xs, ys;
        for (int b = 0; b < nBins; ++b) {
            if (theta_count[b] > 0) {
                xs.push_back(theta_mean[b]);
                ys.push_back(mu[i][b]);
            }
        }
        gr[i] = new TGraph(xs.size(), xs.data(), ys.data());
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack + i);
        if (i == 0) gr[i]->Draw("AP");
        else        gr[i]->Draw("P SAME");
    }
    TLine* zero = new TLine(5, 0, 32, 0);
    zero->SetLineColor(kGray);
    zero->SetLineStyle(2);
    zero->Draw("SAME");

    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(5, 32);
    gr[0]->GetYaxis()->SetRangeUser(-0.1, 0.1);

    // legend on final pad
    TLegend* leg12 = new TLegend(0.20, 0.75, 0.95, 0.90);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        leg12->AddEntry(gr[i], corrLabels[i], "lep");
    }
    leg12->Draw();

    // 13) Save & cleanup
    TString outname = TString::Format(
        "output/dvcs_sebastian_%s_energy_loss_validation.pdf",
        titleSuffix
    );
    c1->SaveAs(outname);

    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        f[i]->Close();
        delete f[i];
    }
}

void plot_dvcs_sebastian_energy_loss_Emiss2_validation(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* titleSuffix) {
    // 0) Setup
    const int nFiles = 3;
    const char* files[nFiles] = { file1, file2, file3 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "SA e^{-} + #gamma",
        "SA e^{-} + #gamma, TBH p"
    };

    // 1) Open files & trees
    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = f[i]
                  ? static_cast<TTree*>(f[i]->Get("PhysicsEvents"))
                  : nullptr;
        if (!tree[i]) {
            std::cerr << "[sebastian] ERROR opening PhysicsEvents in "
                      << files[i] << "\n";
            return;
        }
    }

    // 2) Branch variables
    Double_t p2_theta[nFiles], Mx2_1[nFiles];
    Double_t eta2[nFiles], t1[nFiles], theta_gamma_gamma[nFiles];
    Double_t Emiss2[nFiles], pTmiss[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("p2_theta",           &p2_theta[i]);
        tree[i]->SetBranchAddress("Mx2_1",              &Mx2_1[i]);
        tree[i]->SetBranchAddress("eta2",               &eta2[i]);
        tree[i]->SetBranchAddress("t1",                 &t1[i]);
        tree[i]->SetBranchAddress("theta_gamma_gamma",  &theta_gamma_gamma[i]);
        tree[i]->SetBranchAddress("Emiss2",             &Emiss2[i]);
        tree[i]->SetBranchAddress("pTmiss",             &pTmiss[i]);
    }

    // 3) Canvas & layout
    TCanvas* c1 = new TCanvas(
        "c1",
        "DVCS Energy Loss Validation (Sebastian)",
        1200, 900
    );
    c1->Divide(4, 3);

    // 4) θ bins: 5 → 32 deg, 10 equal bins
    const int    nBins = 10;
    Double_t     thetaBins[nBins+1];
    for (int b = 0; b <= nBins; ++b) {
        thetaBins[b] = 5.0 + b * (32.0 - 5.0) / nBins;
    }

    // 5) Histogram parameters
    const int    nbHi     = 35;
    const int    nbLo     = nbHi/2;  // 17
    const double mx2_min  = -1;
    const double mx2_max  = +2;

    // 6) Storage
    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];
    Double_t theta_sum[nBins]   = {0};
    Int_t    theta_count[nBins] = {0};
    Double_t theta_mean[nBins]  = {0};

    // 7) Create histograms (no stats box)
    for (int i = 0; i < nFiles; ++i) {
        // integrated
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [5,32] %s", titleSuffix),
            nbHi, mx2_min, mx2_max
        );
        h[i][0]->SetStats(false);

        // per-θ slices
        for (int b = 0; b < nBins; ++b) {
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.1f,%.1f] %s",
                     thetaBins[b], thetaBins[b+1], titleSuffix),
                (b < nBins/2 ? nbLo : nbHi),
                mx2_min, mx2_max
            );
            h[i][b+1]->SetStats(false);
        }
    }

    // 8) Fill & accumulate (with kinematic cuts)
    for (int i = 0; i < nFiles; ++i) {
        Long64_t N = tree[i]->GetEntries();
        for (Long64_t ev = 0; ev < N; ++ev) {
            tree[i]->GetEntry(ev);
            double θ = p2_theta[i] * 180.0 / TMath::Pi();
            if (θ >= 5.0 && θ < 32.0 &&
                eta2[i] <  0    &&
                t1[i]   > -2 &&
                theta_gamma_gamma[i] < 0.6 
                // &&
                // Emiss2[i] < 0.5 &&
                // pTmiss[i] < 0.125
                )
            {
                h[i][0]->Fill(Emiss2[i]);
                for (int b = 0; b < nBins; ++b) {
                    if (θ >= thetaBins[b] && θ < thetaBins[b+1]) {
                        h[i][b+1]->Fill(Emiss2[i]);
                        theta_sum[b]   += θ;
                        theta_count[b] += 1;
                    }
                }
            }
        }
    }

    // 9) Compute <θ>
    for (int b = 0; b < nBins; ++b) {
        theta_mean[b] = (theta_count[b] > 0)
                       ? theta_sum[b] / theta_count[b]
                       : NAN;
    }

    // 10) Draw & fit integrated pad (1)
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    double globalMax = 0;
    for (int i = 0; i < nFiles; ++i)
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        if (i == 0) h[i][0]->Draw("E");
        else        h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(
            Form("fitInt%d", i),
            "gaus(0)+pol1(3)",
            mx2_min, mx2_max
        );
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(),  // A
            0.0,                          // μ init
            0.1                           // σ init
        );
        fitInt[i]->SetParLimits(1, -1, 1);
        fitInt[i]->SetParLimits(2,  0.0, 0.8);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }
    // integrated legend
    TLegend* legInt = new TLegend(0.20, 0.75, 0.95, 0.90);
    legInt->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        legInt->AddEntry(
            h[i][0],
            Form("%s: #mu=%.3f, #sigma=%.3f",
                 corrLabels[i],
                 fitInt[i]->GetParameter(1),
                 fitInt[i]->GetParameter(2)),
            "lep"
        );
    }
    legInt->Draw();
    h[0][0]->GetXaxis()->SetTitle("E_{miss}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // 11) θ‐binned pads & fits (2–11)
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        double binMax = 0;
        for (int i = 0; i < nFiles; ++i)
            binMax = std::max(binMax, h[i][b]->GetMaximum());

        for (int i = 0; i < nFiles; ++i) {
            h[i][b]->SetMaximum(1.7 * binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20 + i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack + i);
            if (i == 0) h[i][b]->Draw("E");
            else        h[i][b]->Draw("E SAME");

            TF1* fbin = new TF1(
                Form("fitBin%d_%d", i, b),
                "gaus(0)+pol1(3)",
                mx2_min, mx2_max
            );
            fbin->SetParameters(
                0.8 * h[i][b]->GetMaximum(),
                0.0,
                0.1
            );
            fbin->SetParLimits(1, -1, 1);
            fbin->SetParLimits(2,  0.0, 0.8);
            fbin->SetLineColor(kBlack + i);
            fbin->SetLineWidth(1);
            h[i][b]->Fit(fbin, "Q");
            fbin->Draw("SAME");

            mu[i][b-1]    = fbin->GetParameter(1);
            sigma[i][b-1] = fbin->GetParameter(2);
        }
        // per-θ legend
        TLegend* legB = new TLegend(0.20, 0.75, 0.95, 0.90);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("E_{miss}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // 12) final pad: μ vs. <θ>
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);
    TGraph* gr[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        std::vector<double> xs, ys;
        for (int b = 0; b < nBins; ++b) {
            if (theta_count[b] > 0) {
                xs.push_back(theta_mean[b]);
                ys.push_back(mu[i][b]);
            }
        }
        gr[i] = new TGraph(xs.size(), xs.data(), ys.data());
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack + i);
        if (i == 0) gr[i]->Draw("AP");
        else        gr[i]->Draw("P SAME");
    }
    TLine* zero = new TLine(5, 0, 32, 0);
    zero->SetLineColor(kGray);
    zero->SetLineStyle(2);
    zero->Draw("SAME");

    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(5, 32);
    gr[0]->GetYaxis()->SetRangeUser(-0.1, 0.1);

    // legend on final pad
    TLegend* leg12 = new TLegend(0.20, 0.75, 0.95, 0.90);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        leg12->AddEntry(gr[i], corrLabels[i], "lep");
    }
    leg12->Draw();

    // 13) Save & cleanup
    TString outname = TString::Format(
        "output/dvcs_sebastian_%s_energy_loss_Emiss2_validation.pdf",
        titleSuffix
    );
    c1->SaveAs(outname);

    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        f[i]->Close();
        delete f[i];
    }
}

void plot_eppi0_sebastian_energy_loss_validation(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* titleSuffix,
    bool plotPi0Mass /* = false */) {
    const int nFiles = 3;
    const char* files[nFiles] = { file1, file2, file3 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "SA #gamma",
        "SA e^{-} + #gamma"
    };

    // choose branch name, range, axis‐label and output filename
    const char* branchName  = plotPi0Mass ? "Mh_gammagamma" : "Mx2_2";
    const double rng_min    = plotPi0Mass ? 0.11 : 0.40;
    const double rng_max    = plotPi0Mass ? 0.16 : 1.40;
    const char* xAxisTitle  = plotPi0Mass
                              ? "M_{#gamma#gamma} (GeV)"
                              : "M_{x (e#pi^{0})}^{2} (GeV^{2})";
    const char* outPrefix   = plotPi0Mass
                              ? "output/eppi0_sebastian_pi0mass_%s_energy_loss_validation.pdf"
                              : "output/eppi0_sebastian_%s_energy_loss_validation.pdf";

    // how many datasets to draw
    const int nDraw = plotPi0Mass ? 2 : nFiles;

    std::cout << "[sebastian] plotPi0Mass=" << plotPi0Mass
              << "  branch=" << branchName << "\n";

    // 1) Open files & trees
    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = f[i]
                  ? static_cast<TTree*>(f[i]->Get("PhysicsEvents"))
                  : nullptr;
        if (!tree[i]) {
            std::cerr << "[sebastian] ERROR: cannot open PhysicsEvents in "
                      << files[i] << "\n";
            return;
        }
    }

    // 2) Branch variables
    Double_t p1_theta[nFiles], p2_theta[nFiles], mass[nFiles];
    Double_t eta2[nFiles], t1[nFiles], theta_pi0_pi0[nFiles];
    Double_t Emiss2[nFiles], pTmiss[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("p1_theta",      &p1_theta[i]);
        tree[i]->SetBranchAddress("p2_theta",      &p2_theta[i]);
        tree[i]->SetBranchAddress(branchName,      &mass[i]);
        tree[i]->SetBranchAddress("eta2",          &eta2[i]);
        tree[i]->SetBranchAddress("t1",            &t1[i]);
        tree[i]->SetBranchAddress("theta_pi0_pi0", &theta_pi0_pi0[i]);
        tree[i]->SetBranchAddress("Emiss2",        &Emiss2[i]);
        tree[i]->SetBranchAddress("pTmiss",        &pTmiss[i]);
    }

    // 3) Canvas & layout
    TCanvas* c1 = new TCanvas(
        "c1",
        "eppi0 Energy Loss Validation (Sebastian)",
        1200, 900
    );
    c1->Divide(4, 3);

    // 4) θ‐bins
    const int nBins = 10;
    Double_t thetaBins[nBins+1];
    if (!plotPi0Mass) {
        // original bins: [5,10],[10,15],…,[50,100]
        // Double_t def[11] = {5,10,15,20,25,30,35,40,45,50,100};
        // for (int b = 0; b <= nBins; ++b) thetaBins[b] = def[b];
        // equal 10 bins from 5 to 32 deg
        for (int b = 0; b <= nBins; ++b)
            thetaBins[b] = 0.0 + b * (32.0 - 0.0) / nBins;
    } else {
        // equal 10 bins from 5 to 32 deg
        for (int b = 0; b <= nBins; ++b)
            thetaBins[b] = 0.0 + b * (32.0 - 0.0) / nBins;
    }
    const int nbHi = 35/2, nbLo = 35/2;

    // 5) Storage
    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];
    Double_t theta_sum[nBins]   = {0};
    Int_t    theta_count[nBins] = {0};
    Double_t theta_mean[nBins]  = {0};

    // 6) Create histograms (no stats)
    for (int i = 0; i < nFiles; ++i) {
        // integrated over all θ
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [%.0f,%.0f] %s",
                 thetaBins[0], thetaBins[nBins], titleSuffix),
            nbHi, rng_min, rng_max
        );
        h[i][0]->SetStats(false);

        // per‐θ slices
        for (int b = 0; b < nBins; ++b) {
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.1f,%.1f] %s",
                     thetaBins[b], thetaBins[b+1], titleSuffix),
                (b < nBins/2 ? nbLo : nbHi),
                rng_min, rng_max
            );
            h[i][b+1]->SetStats(false);
        }
    }

    // 7) Fill & accumulate
    for (int i = 0; i < nFiles; ++i) {
        Long64_t N = tree[i]->GetEntries();
        for (Long64_t ev = 0; ev < N; ++ev) {
            tree[i]->GetEntry(ev);
            // choose which theta to use
            Double_t θ = (plotPi0Mass ? p2_theta[i] : p1_theta[i])
                         * 180.0 / TMath::Pi();
            if (θ >= thetaBins[0] && θ < thetaBins[nBins] &&
                eta2[i] < 0 &&
                t1[i]   > -1.79 &&
                theta_pi0_pi0[i] < 0.6 &&
                Emiss2[i] < 0.5 &&
                pTmiss[i] < 0.125)
            {
                h[i][0]->Fill(mass[i]);
                for (int b = 0; b < nBins; ++b) {
                    if (θ >= thetaBins[b] && θ < thetaBins[b+1]) {
                        h[i][b+1]->Fill(mass[i]);
                        theta_sum[b]   += θ;
                        theta_count[b] += 1;
                    }
                }
            }
        }
    }

    // 8) Compute <θ>
    for (int b = 0; b < nBins; ++b) {
        theta_mean[b] = theta_count[b] > 0
                       ? theta_sum[b] / theta_count[b]
                       : NAN;
    }

    // 9) Draw & fit integrated pad (1)
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    Double_t globalMax = 0;
    for (int i = 0; i < nFiles; ++i)
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    for (int i = 0; i < nDraw; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        if (i == 0) h[i][0]->Draw("E");
        else        h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(
            Form("fitInt%d", i),
            "gaus(0)+pol0(3)",
            rng_min, rng_max
        );
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(),
            (plotPi0Mass ? 0.135 : 0.88),
            0.1
        );
        fitInt[i]->SetParLimits(1,
            plotPi0Mass ? 0.11 : 0.75,
            plotPi0Mass ? 0.16 : 1.00
        );
        fitInt[i]->SetParLimits(2, 0.00, 0.02);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }

    // 10) Integrated legend (wider)
    TLegend* legInt = new TLegend(0.20, 0.75, 0.95, 0.90);
    legInt->SetTextSize(0.03);
    for (int i = 0; i < nDraw; ++i) {
        legInt->AddEntry(
            h[i][0],
            Form("%s: #mu=%.3f, #sigma=%.3f",
                 corrLabels[i],
                 fitInt[i]->GetParameter(1),
                 fitInt[i]->GetParameter(2)),
            "lep"
        );
    }
    legInt->Draw();

    h[0][0]->GetXaxis()->SetTitle(xAxisTitle);
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // 11) θ‐binned pads (2–11)
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        Double_t binMax = 0;
        for (int i = 0; i < nFiles; ++i)
            binMax = std::max(binMax, h[i][b]->GetMaximum());

        for (int i = 0; i < nDraw; ++i) {
            h[i][b]->SetMaximum(1.7 * binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20 + i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack + i);
            if (i == 0) h[i][b]->Draw("E");
            else        h[i][b]->Draw("E SAME");

            TF1* fbin = new TF1(
                Form("fitBin%d_%d", i, b),
                "gaus(0)+pol0(3)",
                rng_min, rng_max
            );
            fbin->SetParameters(
                0.8 * h[i][b]->GetMaximum(),
                (plotPi0Mass ? 0.135 : 0.88),
                0.1
            );
            fbin->SetParLimits(1,
                plotPi0Mass ? 0.11 : 0.75,
                plotPi0Mass ? 0.16 : 1.00
            );
            fbin->SetParLimits(2, 0.00, 0.02);
            fbin->SetLineColor(kBlack + i);
            fbin->SetLineWidth(1);
            h[i][b]->Fit(fbin, "Q");
            fbin->Draw("SAME");

            mu[i][b-1]    = fbin->GetParameter(1);
            sigma[i][b-1] = fbin->GetParameter(2);
        }

        TLegend* legB = new TLegend(0.20, 0.75, 0.95, 0.90);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nDraw; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();

        h[0][b]->GetXaxis()->SetTitle(xAxisTitle);
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // 12) Final pad: μ vs. <θ>
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);
    TGraph* gr[nFiles];
    for (int i = 0; i < nDraw; ++i) {
        std::vector<double> xs, ys;
        for (int b = 0; b < nBins; ++b) {
            if (theta_count[b] > 0) {
                xs.push_back(theta_mean[b]);
                ys.push_back(mu[i][b]);
            }
        }
        gr[i] = new TGraph(xs.size(), xs.data(), ys.data());
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack + i);
        if (i == 0) gr[i]->Draw("AP");
        else        gr[i]->Draw("P SAME");
    }

    const double pi0Mass = 0.135;
    // zero line
    TLine* zero = new TLine(0, 0.937*0.937, 35, 0.937*0.937);
    zero->SetLineColor(kGray);
    zero->SetLineStyle(2);
    zero->Draw("SAME");

    // π0‐mass line
    TLine* pi0line = new TLine(0, pi0Mass, 35, pi0Mass);
    pi0line->SetLineColor(kGray);
    pi0line->SetLineStyle(2);
    pi0line->Draw("SAME");

    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV)");
    gr[0]->GetXaxis()->SetLimits(
        plotPi0Mass ? 0 : 0,
        plotPi0Mass ? 35 : 35
    );
    gr[0]->GetYaxis()->SetRangeUser(
        plotPi0Mass ? 0.131 : 0.5,
        plotPi0Mass ? 0.139 : 1.3 
    );

    TLegend* leg12 = new TLegend(0.20, 0.75, 0.95, 0.90);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nDraw; ++i) {
        leg12->AddEntry(gr[i], corrLabels[i], "lep");
    }
    leg12->Draw();

    // 13) Save & cleanup
    TString outname = TString::Format(outPrefix, titleSuffix);
    c1->SaveAs(outname);

    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        f[i]->Close();
        delete f[i];
    }
}

void plot_exclusive_pip_energy_loss_validation(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* titleSuffix) {
    // labels for our three datasets
    const int nFiles = 3;
    const char* files[nFiles] = { file1, file2, file3 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "e^{-} Corrections",
        "e^{-} + #pi Corrections"
    };

    // 1) open files & trees
    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = f[i] ? (TTree*)f[i]->Get("PhysicsEvents") : nullptr;
        if (!tree[i]) {
            std::cerr << "[exclusive π⁺] ERROR opening PhysicsEvents in "
                      << files[i] << "\n";
            return;
        }
    }

    // 2) branch variables: use p_theta now
    Double_t p_theta[nFiles], Mx2[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("p_theta", &p_theta[i]);
        tree[i]->SetBranchAddress("Mx2",     &Mx2[i]);
    }

    // 3) make canvas & divide
    TCanvas* c1 = new TCanvas(
        "c1",
        "exclusive #pi^{+} Energy Loss Validation",
        1200, 900
    );
    c1->Divide(4, 3);

    // 4) θ-bins from 5 to 40° in 10 equal bins
    const int nBins = 10;
    Double_t thetaBins[nBins+1];
    for (int b = 0; b <= nBins; ++b) {
        thetaBins[b] = 5.0 + b * (40.0 - 5.0) / nBins;
    }

    // 5) Mx2 histogram params
    const int    nbHi  = 35;
    const int    nbLo  = nbHi/2; // 17
    const double xMin  = 0.5;
    const double xMax  = 1.2;

    // 6) storage
    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];
    Double_t theta_sum[nBins]   = {0};
    Int_t    theta_count[nBins] = {0};
    Double_t theta_mean[nBins]  = {0};

    // 7) create histograms (no stats box)
    for (int i = 0; i < nFiles; ++i) {
        // integrated over θ
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [5,40] %s (exclusive #pi^{+})", titleSuffix),
            nbHi, xMin, xMax
        );
        h[i][0]->SetStats(false);

        // per-θ slices
        for (int b = 0; b < nBins; ++b) {
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.1f,%.1f] %s (exclusive #pi^{+})",
                     thetaBins[b], thetaBins[b+1], titleSuffix),
                (b < nBins/2 ? nbLo : nbHi),
                xMin, xMax
            );
            h[i][b+1]->SetStats(false);
        }
    }

    // 8) fill & accumulate θ sums
    for (int i = 0; i < nFiles; ++i) {
        Long64_t N = tree[i]->GetEntries();
        for (Long64_t ev = 0; ev < N; ++ev) {
            tree[i]->GetEntry(ev);
            double θ = p_theta[i] * 180.0 / TMath::Pi();
            if (θ < thetaBins[0] || θ >= thetaBins[nBins]) continue;

            h[i][0]->Fill(Mx2[i]);
            for (int b = 0; b < nBins; ++b) {
                if (θ >= thetaBins[b] && θ < thetaBins[b+1]) {
                    h[i][b+1]->Fill(Mx2[i]);
                    theta_sum[b]   += θ;
                    theta_count[b] += 1;
                }
            }
        }
    }

    // 9) compute <θ>
    for (int b = 0; b < nBins; ++b) {
        theta_mean[b] = (theta_count[b] > 0)
                       ? theta_sum[b] / theta_count[b]
                       : NAN;
    }

    // 10) draw & fit integrated pad (1)
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    double globalMax = 0;
    for (int i = 0; i < nFiles; ++i)
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        if (i == 0) h[i][0]->Draw("E");
        else        h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(
            Form("fitInt%d", i),
            "gaus(0)+pol1(3)",
            xMin, xMax
        );
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(),
            (xMin + xMax)/2.0,  // μ init
            0.1                 // σ init
        );
        fitInt[i]->SetParLimits(1, xMin, xMax);
        fitInt[i]->SetParLimits(2, 0.0, (xMax - xMin)/2.0);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }

    // 11) integrated legend
    TLegend* legInt = new TLegend(0.15, 0.75, 0.9, 0.9);
    legInt->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        legInt->AddEntry(
            h[i][0],
            Form("%s: #mu=%.3f, #sigma=%.3f",
                 corrLabels[i],
                 fitInt[i]->GetParameter(1),
                 fitInt[i]->GetParameter(2)),
            "lep"
        );
    }
    legInt->Draw();
    h[0][0]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // 12) θ-binned pads & fits (2–11)
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        double binMax = 0;
        for (int i = 0; i < nFiles; ++i)
            binMax = std::max(binMax, h[i][b]->GetMaximum());
        for (int i = 0; i < nFiles; ++i) {
            h[i][b]->SetMaximum(1.7 * binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20 + i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack + i);
            if (i == 0) h[i][b]->Draw("E");
            else        h[i][b]->Draw("E SAME");

            TF1* fbin = new TF1(
                Form("fitBin%d_%d", i, b),
                "gaus(0)+pol1(3)",
                xMin, xMax
            );
            fbin->SetParameters(
                0.8 * h[i][b]->GetMaximum(),
                (xMin + xMax)/2.0,
                0.1
            );
            fbin->SetParLimits(1, xMin, xMax);
            fbin->SetParLimits(2, 0.0, (xMax - xMin)/2.0);
            fbin->SetLineColor(kBlack + i);
            fbin->SetLineWidth(1);
            h[i][b]->Fit(fbin, "Q");
            fbin->Draw("SAME");
            mu[i][b-1]    = fbin->GetParameter(1);
            sigma[i][b-1] = fbin->GetParameter(2);
        }

        TLegend* legB = new TLegend(0.15, 0.75, 0.9, 0.9);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // 13) final pad: μ vs. <θ> (pad 12)
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);
    TGraph* gr[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        std::vector<double> xs, ys;
        for (int b = 0; b < nBins; ++b) {
            if (theta_count[b] > 0) {
                xs.push_back(theta_mean[b]);
                ys.push_back(mu[i][b]);
            }
        }
        gr[i] = new TGraph(xs.size(), xs.data(), ys.data());
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack + i);
        if (i == 0) gr[i]->Draw("AP");
        else        gr[i]->Draw("P SAME");
    }

    // dashed zero line
    TLine* zero = new TLine(5, 0, 40, 0);
    zero->SetLineColor(kGray);
    zero->SetLineStyle(2);
    zero->Draw("SAME");

    // horizontal line at neutron mass² ≃ (0.9396 GeV)²
    const double neutron2 = 0.9396*0.9396;
    TLine* neutronLine = new TLine(5, neutron2, 40, neutron2);
    neutronLine->SetLineColor(kGray);
    neutronLine->SetLineStyle(2);
    neutronLine->Draw("SAME");

    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(5, 40);
    gr[0]->GetYaxis()->SetRangeUser(0.75, 1.05);

    TLegend* leg12 = new TLegend(0.4, 0.75, 0.9, 0.9);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        leg12->AddEntry(gr[i], corrLabels[i], "lep");
    }
    leg12->Draw();

    // 14) save & cleanup
    TString outname = TString::Format(
        "output/exclusive_pip_%s_energy_loss_validation.pdf",
        titleSuffix
    );
    c1->SaveAs(outname);

    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        f[i]->Close();
        delete f[i];
    }
}

void plot_exclusive_twopion_energy_loss_validation(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* titleSuffix) {
    // labels for our three datasets
    const int nFiles = 3;
    const char* files[nFiles] = { file1, file2, file3 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "e^{-} Corrections",
        "e^{-} + #pi Corrections"
    };

    // 1) open files & trees
    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = f[i] ? (TTree*)f[i]->Get("PhysicsEvents") : nullptr;
        if (!tree[i]) {
            std::cerr << "[exclusive π⁺π⁻] ERROR opening PhysicsEvents in "
                      << files[i] << "\n";
            return;
        }
    }

    // 2) branch variables: now p2_theta and Mx2_23
    Double_t p2_theta[nFiles], Mx2_23[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("p2_theta", &p2_theta[i]);
        tree[i]->SetBranchAddress("Mx2_23",   &Mx2_23[i]);
    }

    // 3) make canvas & divide
    TCanvas* c1 = new TCanvas(
        "c1",
        "exclusive #pi^{+}#pi^{-} Energy Loss Validation",
        1200, 900
    );
    c1->Divide(4, 3);

    // 4) θ-bins from 5 to 40° in 10 equal bins
    const int nBins = 10;
    Double_t thetaBins[nBins+1];
    for (int b = 0; b <= nBins; ++b) {
        thetaBins[b] = 5.0 + b * (40.0 - 5.0) / nBins;
    }

    // 5) Mx2_23 histogram params
    const int    nbHi  = 35;
    const int    nbLo  = nbHi/2; // 17
    const double xMin  = 0.5;
    const double xMax  = 1.2;

    // 6) storage
    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    Double_t mu[nFiles][nBins], sigma[nFiles][nBins];
    Double_t theta_sum[nBins]   = {0};
    Int_t    theta_count[nBins] = {0};
    Double_t theta_mean[nBins]  = {0};

    // 7) create histograms (no stats box)
    for (int i = 0; i < nFiles; ++i) {
        // integrated over θ
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [5,40] %s (exclusive #pi^{+}#pi^{-})", titleSuffix),
            nbHi, xMin, xMax
        );
        h[i][0]->SetStats(false);

        // per-θ slices
        for (int b = 0; b < nBins; ++b) {
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.1f,%.1f] %s (exclusive #pi^{+}#pi^{-})",
                     thetaBins[b], thetaBins[b+1], titleSuffix),
                (b < nBins/2 ? nbLo : nbHi),
                xMin, xMax
            );
            h[i][b+1]->SetStats(false);
        }
    }

    // 8) fill & accumulate θ sums
    for (int i = 0; i < nFiles; ++i) {
        Long64_t N = tree[i]->GetEntries();
        for (Long64_t ev = 0; ev < N; ++ev) {
            tree[i]->GetEntry(ev);
            double θ = p2_theta[i] * 180.0 / TMath::Pi();
            if (θ < thetaBins[0] || θ >= thetaBins[nBins]) continue;

            h[i][0]->Fill(Mx2_23[i]);
            for (int b = 0; b < nBins; ++b) {
                if (θ >= thetaBins[b] && θ < thetaBins[b+1]) {
                    h[i][b+1]->Fill(Mx2_23[i]);
                    theta_sum[b]   += θ;
                    theta_count[b] += 1;
                }
            }
        }
    }

    // 9) compute <θ>
    for (int b = 0; b < nBins; ++b) {
        theta_mean[b] = (theta_count[b] > 0)
                       ? theta_sum[b] / theta_count[b]
                       : NAN;
    }

    // 10) draw & fit integrated pad (1)
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    double globalMax = 0;
    for (int i = 0; i < nFiles; ++i)
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        if (i == 0) h[i][0]->Draw("E");
        else        h[i][0]->Draw("E SAME");

        fitInt[i] = new TF1(
            Form("fitInt%d", i),
            "gaus(0)+pol1(3)",
            xMin, xMax
        );
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(),
            (xMin + xMax)/2.0,  // μ init
            0.1                 // σ init
        );
        fitInt[i]->SetParLimits(1, xMin, xMax);
        fitInt[i]->SetParLimits(2, 0.0, (xMax - xMin)/2.0);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }

    // 11) integrated legend
    TLegend* legInt = new TLegend(0.20, 0.75, 0.9, 0.9);
    legInt->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        legInt->AddEntry(
            h[i][0],
            Form("%s: #mu=%.3f, #sigma=%.3f",
                 corrLabels[i],
                 fitInt[i]->GetParameter(1),
                 fitInt[i]->GetParameter(2)),
            "lep"
        );
    }
    legInt->Draw();
    h[0][0]->GetXaxis()->SetTitle("M_{23}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // 12) θ-binned pads & fits (2–11)
    for (int b = 1; b <= nBins; ++b) {
        c1->cd(b+1)->SetLeftMargin(0.15);
        c1->cd(b+1)->SetBottomMargin(0.15);
        double binMax = 0;
        for (int i = 0; i < nFiles; ++i)
            binMax = std::max(binMax, h[i][b]->GetMaximum());
        for (int i = 0; i < nFiles; ++i) {
            h[i][b]->SetMaximum(1.7 * binMax);
            h[i][b]->SetMinimum(0);
            h[i][b]->SetMarkerStyle(20 + i);
            h[i][b]->SetMarkerSize(0.8);
            h[i][b]->SetMarkerColor(kBlack + i);
            if (i == 0) h[i][b]->Draw("E");
            else        h[i][b]->Draw("E SAME");

            TF1* fbin = new TF1(
                Form("fitBin%d_%d", i, b),
                "gaus(0)+pol1(3)",
                xMin, xMax
            );
            fbin->SetParameters(
                0.8 * h[i][b]->GetMaximum(),
                (xMin + xMax)/2.0,
                0.1
            );
            fbin->SetParLimits(1, xMin, xMax);
            fbin->SetParLimits(2, 0.0, (xMax - xMin)/2.0);
            fbin->SetLineColor(kBlack + i);
            fbin->SetLineWidth(1);
            h[i][b]->Fit(fbin, "Q");
            fbin->Draw("SAME");
            mu[i][b-1]    = fbin->GetParameter(1);
            sigma[i][b-1] = fbin->GetParameter(2);
        }

        TLegend* legB = new TLegend(0.20, 0.75, 0.9, 0.9);
        legB->SetTextSize(0.03);
        for (int i = 0; i < nFiles; ++i) {
            legB->AddEntry(
                h[i][b],
                Form("%s: #mu=%.3f, #sigma=%.3f",
                     corrLabels[i],
                     mu[i][b-1],
                     sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("M_{23}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // 13) final pad: μ vs. <θ> (pad 12)
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);
    TGraph* gr[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        std::vector<double> xs, ys;
        for (int b = 0; b < nBins; ++b) {
            if (theta_count[b] > 0) {
                xs.push_back(theta_mean[b]);
                ys.push_back(mu[i][b]);
            }
        }
        gr[i] = new TGraph(xs.size(), xs.data(), ys.data());
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack + i);
        if (i == 0) gr[i]->Draw("AP");
        else        gr[i]->Draw("P SAME");
    }

    // dashed zero line
    TLine* zero = new TLine(5, 0, 40, 0);
    zero->SetLineColor(kGray);
    zero->SetLineStyle(2);
    zero->Draw("SAME");

    // horizontal line at proton mass² ≃ (0.9383 GeV)²
    const double proton2 = 0.9383*0.9383;
    TLine* protonLine = new TLine(5, proton2, 40, proton2);
    protonLine->SetLineColor(kGray);
    protonLine->SetLineStyle(2);
    protonLine->Draw("SAME");

    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(5, 40);
    gr[0]->GetYaxis()->SetRangeUser(0.75, 1.05);

    TLegend* leg12 = new TLegend(0.35, 0.75, 0.9, 0.9);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        leg12->AddEntry(gr[i], corrLabels[i], "lep");
    }
    leg12->Draw();

    // 14) save & cleanup
    TString outname = TString::Format(
        "output/exclusive_twopion_%s_energy_loss_validation.pdf",
        titleSuffix
    );
    c1->SaveAs(outname);

    delete c1;
    for (int i = 0; i < nFiles; ++i) {
        f[i]->Close();
        delete f[i];
    }
}

int main(int argc, char** argv) {
    std::cout << std::flush; //test
    if (!(argc == 4 || argc == 5 || argc == 6)) {
        std::cerr << "Usage: " << argv[0]
                  << " <file1.root> <file2.root> <file3.root> (<file4.root>) <titleSuffix>\n";
        return 1;
    }

    // plot_dvcs_energy_loss_validation(
    //     argv[1], argv[2], argv[3], argv[4], argv[5]
    // );

    // plot_eppi0_energy_loss_validation(
    //     argv[1], argv[2], argv[3], argv[4], argv[5]
    // );

    // plot_mx2_comparison_elastic(
    //     argv[1], argv[2], argv[3], argv[4], argv[5]
    // );

    // plot_two_pions(
    //     argv[1], argv[2], argv[3], argv[4], argv[5]
    // );

    /////////

    plot_dvcs_sebastian_energy_loss_validation(
        argv[1], argv[2], argv[3], argv[4]
    );

    // plot_dvcs_sebastian_energy_loss_Emiss2_validation(
    //     argv[1], argv[2], argv[3], argv[4]
    // );

    // plot_eppi0_sebastian_energy_loss_validation(
    //     argv[1], argv[2], argv[3], argv[4], false
    // );

    // plot_eppi0_sebastian_energy_loss_validation(
    //     argv[1], argv[2], argv[3], argv[4], true
    // );

    /////////

    // plot_exclusive_pip_energy_loss_validation(
    //     argv[1], argv[2], argv[3], argv[4]
    // );

    // plot_exclusive_twopion_energy_loss_validation(
    //     argv[1], argv[2], argv[3], argv[4]
    // );


    return 0;
}