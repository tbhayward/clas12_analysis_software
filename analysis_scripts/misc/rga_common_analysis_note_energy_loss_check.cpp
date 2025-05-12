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
            Form("Integrated #theta [5,65] %s (%s)", titleSuffix, corrLabels[i]),
            nbMx2_hi, mx2_min, mx2_max
        );
        for (int b = 0; b < nBins; ++b) {
            int nb = (b < nBins/2 ? nbMx2_lo : nbMx2_hi);
            h[i][b+1] = new TH1D(
                Form("h%d_%d", i, b),
                Form("#theta [%.0f,%.0f] %s (%s)",
                     thetaBins[b], thetaBins[b+1],
                     titleSuffix, corrLabels[i]),
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

void plot_mx2_comparison_elastic(
    const char* file1,
    const char* file2,
    const char* file3,
    const char* file4,
    const char* titleSuffix
) {
    const int nFiles = 4;
    const char* files[nFiles] = { file1, file2, file3, file4 };
    const char* corrLabels[nFiles] = {
        "No Corrections",
        "Timothy's",
        "Krishna's",
        "Mariana's"
    };

    // --- open files and get trees, with null checks
    TFile* f[nFiles];
    TTree* tree[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        f[i] = TFile::Open(files[i]);
        if (!f[i] || f[i]->IsZombie()) {
            std::cerr << "Error: could not open file " << files[i] << "\n";
            tree[i] = nullptr;
            continue;
        }
        tree[i] = dynamic_cast<TTree*>(f[i]->Get("PhysicsEvents"));
        if (!tree[i]) {
            std::cerr << "Error: no PhysicsEvents tree in " << files[i] << "\n";
        }
    }

    // --- set up branches only for valid trees
    Double_t p_theta[nFiles], Mx2[nFiles];
    Int_t    det[nFiles];
    for (int i = 0; i < nFiles; ++i) {
        if (!tree[i]) continue;
        tree[i]->SetBranchAddress("p_theta",  &p_theta[i]);
        tree[i]->SetBranchAddress("Mx2",      &Mx2[i]);
        tree[i]->SetBranchAddress("detector", &det[i]);
    }

    // --- detector definitions
    struct DetConfig { const char* name; std::vector<double> bins; std::vector<std::string> labels; };
    DetConfig dets[2] = {
        { "Forward", {0,8,11,14,17,20,23,26,29,32,35,38,41,80}, {} },
        { "Central",{0,36,39,42,45,48,51,54,57,180},           {} }
    };
    for (auto &dc : dets) {
        dc.labels.push_back("All θ");
        for (size_t i = 1; i+1 < dc.bins.size(); ++i)
            dc.labels.push_back(
                std::to_string((int)dc.bins[i]) + "-" +
                std::to_string((int)dc.bins[i+1])
            );
    }

    // --- histogram parameters
    const Double_t mx2_min = 0.5, mx2_max = 1.2;
    const int    nbMx2_hi  = 50, nbMx2_lo = nbMx2_hi/2;

    // --- loop over detectors
    for (int d = 0; d < 2; ++d) {
        auto &DC = dets[d];
        int detNum = d+1;
        int nPlots = DC.labels.size();
        int nCols  = 4;
        int nRows  = ((nPlots-1) + nCols-1)/nCols + 1;

        TCanvas* c = new TCanvas(
            Form("c_%s", DC.name),
            Form("Elastic Mx² Comparison — %s", DC.name),
            1200, 300*nRows
        );
        c->Divide(nCols, nRows);

        // --- create histograms
        static const int MAXP = 16;
        TH1D* h[nFiles][MAXP];
        for (int i = 0; i < nFiles; ++i) {
            h[i][0] = new TH1D(
                Form("h_%s_int_%d", DC.name, i),
                Form("Integrated %s (%s)", DC.name, corrLabels[i]),
                nbMx2_hi, mx2_min, mx2_max
            );
            for (int b=0; b<nPlots-1; ++b) {
                int nb = (b < (nPlots-1)/2 ? nbMx2_lo : nbMx2_hi);
                h[i][b+1] = new TH1D(
                    Form("h_%s_%d_%d", DC.name, b, i),
                    Form("θ[%s] %s (%s)",
                         DC.labels[b+1].c_str(), DC.name, corrLabels[i]),
                    nb, mx2_min, mx2_max
                );
            }
        }

        // --- fill them
        for (int i = 0; i < nFiles; ++i) {
            if (!tree[i]) continue;
            Long64_t N = tree[i]->GetEntries();
            for (Long64_t evt=0; evt<N; ++evt) {
                tree[i]->GetEntry(evt);
                if (det[i] != detNum) continue;
                double θdeg = p_theta[i]*180.0/TMath::Pi();
                h[i][0]->Fill(Mx2[i]);
                for (size_t b=0; b+1<DC.bins.size(); ++b) {
                    if (θdeg >= DC.bins[b] && θdeg < DC.bins[b+1])
                        h[i][b+1]->Fill(Mx2[i]);
                }
            }
        }

        // --- draw integrated
        c->cd(1)->SetLeftMargin(0.15)->SetBottomMargin(0.15);
        double gmax=0;
        for (int i=0; i<nFiles; ++i) if (h[i][0]) gmax = std::max(gmax, h[i][0]->GetMaximum());
        for (int i=0; i<nFiles; ++i) {
            if (!h[i][0]) continue;
            h[i][0]->SetMaximum(1.7*gmax);
            h[i][0]->SetLineColor(kBlack+i);
            h[i][0]->SetLineStyle(i);
            h[i][0]->Draw(i==0?"HIST":"HIST SAME");
        }
        auto legInt = new TLegend(0.6,0.7,0.9,0.9);
        for (int i=0; i<nFiles; ++i) legInt->AddEntry(h[i][0], corrLabels[i], "l");
        legInt->Draw();
        h[0][0]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h[0][0]->GetYaxis()->SetTitle("Counts");

        // --- draw theta slices
        for (int p=1; p<nPlots; ++p) {
            int pad = p+1;
            c->cd(pad)->SetLeftMargin(0.15)->SetBottomMargin(0.15);
            double bmax=0;
            for (int i=0; i<nFiles; ++i) if (h[i][p]) bmax = std::max(bmax, h[i][p]->GetMaximum());
            for (int i=0; i<nFiles; ++i) {
                if (!h[i][p]) continue;
                h[i][p]->SetMaximum(1.7*bmax);
                h[i][p]->SetLineColor(kBlack+i);
                h[i][p]->SetLineStyle(i);
                h[i][p]->Draw(i==0?"HIST":"HIST SAME");
            }
            if (p < nCols) {
                auto leg = new TLegend(0.6,0.7,0.9,0.9);
                for (int i=0; i<nFiles; ++i) leg->AddEntry(h[i][p], corrLabels[i], "l");
                leg->Draw();
            }
            h[0][p]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
            h[0][p]->GetYaxis()->SetTitle("Counts");
            h[0][p]->GetXaxis()->SetRangeUser(mx2_min, mx2_max);
            c->cd(pad)->SetTitle(DC.labels[p].c_str());
        }

        // --- save
        TString out = TString::Format(
            "output/Mx2_elastic_comparison_%s_%s.pdf",
            DC.name, titleSuffix
        );
        c->SaveAs(out);
        delete c;
    }

    // --- cleanup
    for (int i = 0; i < nFiles; ++i) {
        if (f[i]) { f[i]->Close(); delete f[i]; }
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

    plot_mx2_comparison_elastic(
        argv[1], argv[2], argv[3], argv[4], argv[5]
    );


    return 0;
}