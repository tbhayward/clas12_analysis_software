/*
 * File: calorimeter_comparison.cpp
 *
 * Compile with:
 *   g++ calorimeter_comparison.cpp -o calorimeter_comparison \
 *       `root-config --cflags --libs`
 *
 * Run with:
 *   ./calorimeter_comparison [Nevents] [dataFile] [mcFile]
 *
 * What it does:
 *   1) Creates output/ and output/cal/ if they do not exist.
 *   2) Loops over data and MC, filling photon (PID=22) and electron (PID=11)
 *      calorimeter hit-position histograms for PCal, ECin, and ECout.
 *   3) Makes unnormalized 2D plots (data vs. MC) with all three layers on one canvas.
 *   4) Normalizes each histogram and computes data/MC ratio.
 *   5) Creates a 3-panel ratio plot with mean (μ) and stddev (σ) annotated per layer.
 *   6) Creates a second 3-panel plot drawing points: blue for 0.5≤ratio≤2, red otherwise,
 *      with no histogram color scale (just axes).
 */

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TString.h"

// ----------------------------------------------------------------------------
// Unify color scale for three TH2D histograms
// ----------------------------------------------------------------------------
void SetSame2DScale(TH2D* a, TH2D* b, TH2D* c) {
    double mn = 1e9, mx = -1e9;
    for (auto* h : {a,b,c}) {
        mn = std::min(mn, h->GetMinimum());
        mx = std::max(mx, h->GetMaximum());
    }
    for (auto* h : {a,b,c}) {
        h->SetMinimum(mn);
        h->SetMaximum(mx);
    }
}

int main(int argc, char** argv) {
    // 1) Create output directories
    gSystem->mkdir("output",    kTRUE);
    gSystem->mkdir("output/cal",kTRUE);

    // 2) Parse arguments
    Long64_t maxEvents = -1;
    if (argc > 1) {
        maxEvents = std::stoll(argv[1]);
        if (maxEvents == 0) maxEvents = -1;
    }
    bool useData = (argc > 2);
    bool useMC   = (argc > 3);
    std::string dataFile = useData ? argv[2] : "";
    std::string mcFile   = useMC   ? argv[3] : "";

    // 3) Set up TChains
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    if (useData) dataCh.Add(dataFile.c_str()); else dataCh.Add("/work/.../cal_data.root");
    if (useMC)   mcCh.Add(mcFile.c_str());   else mcCh.Add("/work/.../cal_mc.root");

    // 4) Branch addresses
    Int_t    pid;
    Double_t x1,y1,x4,y4,x7,y7;
    dataCh.SetBranchAddress("particle_pid", &pid);
    dataCh.SetBranchAddress("cal_x_1",      &x1);
    dataCh.SetBranchAddress("cal_y_1",      &y1);
    dataCh.SetBranchAddress("cal_x_4",      &x4);
    dataCh.SetBranchAddress("cal_y_4",      &y4);
    dataCh.SetBranchAddress("cal_x_7",      &x7);
    dataCh.SetBranchAddress("cal_y_7",      &y7);
    mcCh .SetBranchAddress("particle_pid", &pid);
    mcCh .SetBranchAddress("cal_x_1",      &x1);
    mcCh .SetBranchAddress("cal_y_1",      &y1);
    mcCh .SetBranchAddress("cal_x_4",      &x4);
    mcCh .SetBranchAddress("cal_y_4",      &y4);
    mcCh .SetBranchAddress("cal_x_7",      &x7);
    mcCh .SetBranchAddress("cal_y_7",      &y7);

    // 5) Configuration
    const int NB = 200;
    const double xmin = -450, xmax = 450, ymin = -450, ymax = 450;
    std::vector<std::pair<int,std::string>> species = {{22,"photon"},{11,"electron"}};
    std::vector<std::string> layers = {"PCal","ECin","ECout"};

    for (auto& sp : species) {
        int pidVal = sp.first;
        std::string label = sp.second;

        // Histograms for data, MC
        std::vector<TH2D*> hD(3), hM(3), hR(3);
        for (int i = 0; i < 3; ++i) {
            hD[i] = new TH2D(Form("hD_%s_%s",label.c_str(),layers[i].c_str()),
                "; x (cm); y (cm)", NB, xmin, xmax, NB, ymin, ymax);
            hM[i] = (TH2D*)hD[i]->Clone(Form("hM_%s_%s",label.c_str(),layers[i].c_str()));
        }

        // 6) Fill data
        Long64_t nD = dataCh.GetEntries(); if (maxEvents>0 && maxEvents<nD) nD = maxEvents;
        for (Long64_t e = 0; e < nD; ++e) {
            dataCh.GetEntry(e);
            if (pid != pidVal) continue;
            double xs[3] = {x1, x4, x7};
            double ys[3] = {y1, y4, y7};
            for (int i = 0; i < 3; ++i) {
                if (xs[i] != -9999 && ys[i] != -9999)
                    hD[i]->Fill(xs[i], ys[i]);
            }
        }
        // 7) Fill MC
        Long64_t nM = mcCh.GetEntries(); if (maxEvents>0 && maxEvents<nM) nM = maxEvents;
        for (Long64_t e = 0; e < nM; ++e) {
            mcCh.GetEntry(e);
            if (pid != pidVal) continue;
            double xs[3] = {x1, x4, x7};
            double ys[3] = {y1, y4, y7};
            for (int i = 0; i < 3; ++i) {
                if (xs[i] != -9999 && ys[i] != -9999)
                    hM[i]->Fill(xs[i], ys[i]);
            }
        }

        gStyle->SetOptStat(0);

        // 8) Unnormalized plots
        { TCanvas c1("c_data_uncut","Data Uncut",1800,600); c1.Divide(3,1);
          SetSame2DScale(hD[0],hD[1],hD[2]);
          for(int i=0;i<3;i++){ c1.cd(i+1); gPad->SetLogz(); hD[i]->Draw("COLZ"); }
          c1.SaveAs(Form("output/cal/data_uncut_%s.png",label.c_str())); }
        { TCanvas c2("c_mc_uncut","MC Uncut",1800,600); c2.Divide(3,1);
          SetSame2DScale(hM[0],hM[1],hM[2]);
          for(int i=0;i<3;i++){ c2.cd(i+1); gPad->SetLogz(); hM[i]->Draw("COLZ"); }
          c2.SaveAs(Form("output/cal/mc_uncut_%s.png",label.c_str())); }

        // 9) Normalize & ratio
        for(int i=0;i<3;i++){
            double id=hD[i]->Integral(); if(id>0)hD[i]->Scale(1/id);
            double im=hM[i]->Integral(); if(im>0)hM[i]->Scale(1/im);
            hR[i] = (TH2D*)hD[i]->Clone(Form("hR_%s_%s",label.c_str(),layers[i].c_str()));
            hR[i]->Divide(hM[i]);
        }
        // 10) Stats
        std::vector<double> mu(3), sigma(3);
        for(int i=0;i<3;i++){
            double sum=0,sum2=0; int cnt=0;
            for(int ix=1; ix<=NB; ix++) for(int iy=1; iy<=NB; iy++){
                double v=hR[i]->GetBinContent(ix,iy);
                if(v==0) continue;
                sum+=v; sum2+=v*v; cnt++;
            }
            mu[i] = cnt?sum/cnt:0;
            sigma[i] = cnt?sqrt(sum2/cnt - mu[i]*mu[i]):0;
        }

        // 11) Ratio canvas
        { TCanvas c3("c_ratio","Data/MC Ratio",1800,600); c3.Divide(3,1);
          SetSame2DScale(hR[0],hR[1],hR[2]);
          for(int i=0;i<3;i++){
            c3.cd(i+1); gPad->SetLogz(); hR[i]->Draw("COLZ");
            TLatex t; t.SetNDC(); t.SetTextAlign(31); t.SetTextSize(0.03);
            t.DrawLatex(0.95,0.95,Form("μ=%.3f",mu[i]));
            t.DrawLatex(0.95,0.90,Form("σ=%.3f",sigma[i]));
          }
          c3.SaveAs(Form("output/cal/ratio_%s.png",label.c_str())); }

        // 12) Scatter outliers (red/blue dots) with just axes
        { TCanvas c4("c_outliers","Outliers",1800,600); c4.Divide(3,1);
          for(int i=0;i<3;i++){
            c4.cd(i+1);
            gPad->SetLogz(0);
            // draw only frame axes
            hR[i]->Draw("axis");
            // prepare graphs
            TGraph gIn, gOut;
            gIn.SetMarkerStyle(20); gIn.SetMarkerColor(kBlue);
            gOut.SetMarkerStyle(20); gOut.SetMarkerColor(kRed);
            // fill points
            for(int ix=1; ix<=NB; ix++){
              double px = hR[i]->GetXaxis()->GetBinCenter(ix);
              for(int iy=1; iy<=NB; iy++){
                double py = hR[i]->GetYaxis()->GetBinCenter(iy);
                double v = hR[i]->GetBinContent(ix,iy);
                if(v==0) continue;
                if(v>2.0 || v<0.5) gOut.SetPoint(gOut.GetN(), px, py);
                else               gIn.SetPoint(gIn.GetN(),  px, py);
              }
            }
            // draw dots
            gIn.Draw("P");
            gOut.Draw("P same");
          }
          c4.SaveAs(Form("output/cal/ratio_%s_outliers.png",label.c_str())); }
    }
    return 0;
}
