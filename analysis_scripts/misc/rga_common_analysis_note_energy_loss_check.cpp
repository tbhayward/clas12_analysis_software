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
    const double W_min   = -0.1, W_max = 0.1;

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
        fitInt[i]->SetParLimits(1, 0.8, 1.0);
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
            fitBin[i][b-1]->SetParLimits(1, 0.8, 1.0);
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
        tree[i]->SetBranchAddress("Mx2_13",   &Mx2[i]);
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
    const double mx2_min = -0.12, mx2_max = 0.20;

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
            0.0195, 0.05,
            0,0,0
        );
        fitInt[i]->SetParLimits(1, 0.00, 0.05);
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
                0.0195, 0.05,
                0,0,0
            );
            fitBin[i][b-1]->SetParLimits(1, 0.00, 0.05);
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
    TLine* line = new TLine(thetaBins[0],0,thetaBins[nBins],0);
    line->SetLineColor(kGray);
    line->SetLineStyle(2);
    line->Draw("SAME");

    gr[0]->GetXaxis()->SetTitle("#theta (deg)");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(thetaBins[0], thetaBins[nBins]);
    gr[0]->GetYaxis()->SetRangeUser(-0.01, 0.04);

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

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <file1.root> <file2.root> <file3.root> <file4.root> <titleSuffix>\n";
        return 1;
    }

    // plot_dvcs_energy_loss_validation(
    //     argv[1], argv[2], argv[3], argv[4], argv[5]
    // );

    // plot_eppi0_energy_loss_validation(
    //     argv[1], argv[2], argv[3], argv[4], argv[5]
    // );

    plot_mx2_comparison_elastic(
        argv[1], argv[2], argv[3], argv[4], argv[5]
    );

    // plot_two_pions(
    //     argv[1], argv[2], argv[3], argv[4], argv[5]
    // );


    return 0;
}