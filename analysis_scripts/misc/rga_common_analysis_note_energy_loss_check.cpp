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
    const char* titleSuffix
) {
    const int nFiles = 4;
    const char* files[nFiles] = { file1, file2, file3, file4 };
    TFile*    f[nFiles];
    TTree*    tree[nFiles];

    // open files & trees
    for (int i = 0; i < nFiles; ++i) {
        f[i]    = TFile::Open(files[i]);
        tree[i] = (TTree*)f[i]->Get("PhysicsEvents");
    }

    // branch variables
    Double_t p1_theta[nFiles], Mx2_1[nFiles];
    Double_t eta2[nFiles], t1[nFiles], theta_gamma_gamma[nFiles];
    Double_t Emiss2[nFiles], pTmiss[nFiles];

    // set branch addresses
    for (int i = 0; i < nFiles; ++i) {
        tree[i]->SetBranchAddress("p1_theta",          &p1_theta[i]);
        tree[i]->SetBranchAddress("Mx2_1",             &Mx2_1[i]);
        tree[i]->SetBranchAddress("eta2",              &eta2[i]);
        tree[i]->SetBranchAddress("t1",                &t1[i]);
        tree[i]->SetBranchAddress("theta_gamma_gamma", &theta_gamma_gamma[i]);
        tree[i]->SetBranchAddress("Emiss2",            &Emiss2[i]);
        tree[i]->SetBranchAddress("pTmiss",            &pTmiss[i]);
    }

    // canvas and pads
    TCanvas* c1 = new TCanvas("c1", "DVCS Energy Loss Validation", 1200, 900);
    c1->Divide(4, 3);

    const int    nBins        = 10;
    Double_t     thetaBins[nBins+1] = {5,15,20,25,30,35,40,45,50,60,100};
    const int    nbMx2        = 35;
    const Double_t mx2_min    = -0.3, mx2_max = +0.3;

    // histograms, fit functions, and theta accumulators
    TH1D*    h[nFiles][nBins+1];
    TF1*     fitInt[nFiles];
    Double_t mu[nFiles][nBins];
    Double_t sigma[nFiles][nBins];
    Double_t theta_sum[nBins]   = {0};
    Int_t    theta_count[nBins] = {0};
    Double_t theta_mean[nBins]  = {0};

    // create histograms for each file and bin
    for (int i = 0; i < nFiles; ++i) {
        // integrated
        h[i][0] = new TH1D(
            Form("h%d_int", i),
            Form("Integrated #theta [5,65] %s (file %d)", titleSuffix, i+1),
            nbMx2, mx2_min, mx2_max
        );
        // theta bins
        for (int b = 0; b < nBins; ++b) {
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.0f,%.0f] %s (file %d)",
                     thetaBins[b], thetaBins[b+1], titleSuffix, i+1),
                nbMx2, mx2_min, mx2_max
            );
        }
    }

    // fill histograms with kinematic cuts
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
                pTmiss[i] < 0.125
            ) {
                h[i][0]->Fill(Mx2_1[i]);
                for (int b = 0; b < nBins; ++b) {
                    if (thetaDeg >= thetaBins[b] && thetaDeg < thetaBins[b+1]) {
                        h[i][b+1]->Fill(Mx2_1[i]);
                        theta_sum[b]   += thetaDeg;
                        theta_count[b] += 1;
                    }
                }
            }
        }
    }

    // compute mean theta per bin
    for (int b = 0; b < nBins; ++b) {
        if (theta_count[b] > 0) {
            theta_mean[b] = theta_sum[b] / theta_count[b];
        } else {
            theta_mean[b] = 0.5 * (thetaBins[b] + thetaBins[b+1]);
        }
    }

    // integrated pad (pad 1)
    c1->cd(1);
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);
    Double_t globalMax = 0;
    for (int i = 0; i < nFiles; ++i) {
        globalMax = std::max(globalMax, h[i][0]->GetMaximum());
    }
    // draw & fit integrated
    for (int i = 0; i < nFiles; ++i) {
        h[i][0]->SetMaximum(1.7 * globalMax);
        h[i][0]->SetMinimum(0);
        h[i][0]->SetMarkerStyle(20 + i);
        h[i][0]->SetMarkerSize(0.8);
        h[i][0]->SetMarkerColor(kBlack + i);
        h[i][0]->SetStats(0);
        if (i == 0) {
            h[i][0]->Draw("E");
        } else {
            h[i][0]->Draw("E SAME");
        }

        fitInt[i] = new TF1(
            Form("fitInt%d", i),
            "gaus(0)+pol1(3)",
            mx2_min, mx2_max
        );
        fitInt[i]->SetParameters(
            0.8 * h[i][0]->GetMaximum(), 0, 0.2
        );
        fitInt[i]->SetParLimits(1, -0.15, 0.15);
        fitInt[i]->SetParLimits(2,  0.0, 0.3);
        fitInt[i]->SetLineColor(kBlack + i);
        fitInt[i]->SetLineWidth(1);
        h[i][0]->Fit(fitInt[i], "Q");
        fitInt[i]->Draw("SAME");
    }

    // legend for integrated
    TLegend* legInt = new TLegend(0.25, 0.75, 0.9, 0.9);
    legInt->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        legInt->AddEntry(
            h[i][0],
            Form("File %d: #mu=%.3f, #sigma=%.3f",
                 i+1,
                 fitInt[i]->GetParameter(1),
                 fitInt[i]->GetParameter(2)),
            "lep"
        );
    }
    legInt->Draw();
    h[0][0]->GetXaxis()->SetTitle("M_{xp}^{2} (GeV^{2})");
    h[0][0]->GetYaxis()->SetTitle("Counts");

    // theta‐binned pads (2–11)
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
            if (i == 0) {
                h[i][b]->Draw("E");
            } else {
                h[i][b]->Draw("E SAME");
            }

            TF1* fbin = new TF1(
                Form("fitBin%d_%d", i, b),
                "gaus(0)+pol1(3)",
                mx2_min, mx2_max
            );
            fbin->SetParameters(
                0.8 * h[i][b]->GetMaximum(), 0, 0.2
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
                Form("File %d: #mu=%.3f, #sigma=%.3f",
                     i+1, mu[i][b-1], sigma[i][b-1]),
                "lep"
            );
        }
        legB->Draw();
        h[0][b]->GetXaxis()->SetTitle("M_{xp}^{2} (GeV^{2})");
        h[0][b]->GetYaxis()->SetTitle("Counts");
    }

    // final pad: mu vs theta with no error bars
    c1->cd(12);
    c1->cd(12)->SetLeftMargin(0.20);
    c1->cd(12)->SetBottomMargin(0.15);

    TGraph* gr[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        gr[i] = new TGraph(nBins, theta_mean, mu[i]);
        gr[i]->SetMarkerStyle(20 + i);
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlack + i);
        if (i == 0) {
            gr[i]->Draw("AP");
        } else {
            gr[i]->Draw("P SAME");
        }
    }

    TLine* line = new TLine(0, 0, 100, 0);
    line->SetLineColor(kGray);
    line->SetLineStyle(2);
    line->Draw("SAME");

    gr[0]->GetXaxis()->SetTitle("#theta");
    gr[0]->GetYaxis()->SetTitle("#mu (GeV^{2})");
    gr[0]->GetXaxis()->SetLimits(0, 100);
    gr[0]->GetYaxis()->SetRangeUser(-0.20, 0.20);

    TLegend* leg12 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg12->SetTextSize(0.03);
    for (int i = 0; i < nFiles; ++i) {
        leg12->AddEntry(gr[i], Form("File %d", i+1), "lep");
    }
    leg12->Draw();

    // save and clean up
    TString outname = TString::Format("output/dvcs_%s_energy_loss_validation.pdf", titleSuffix);
    c1->SaveAs(outname);

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
    plot_dvcs_energy_loss_validation(
        argv[1], argv[2], argv[3], argv[4], argv[5]
    );
    return 0;
}