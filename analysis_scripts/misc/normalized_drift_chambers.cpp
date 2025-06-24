/*****************************************************************************
 * File: normalized_drift_chambers_comparison.cpp
 *
 * Compile with:
 *   g++ normalized_drift_chambers_comparison.cpp -o normalized_drift_chambers_comparison \
 *       `root-config --cflags --libs`
 *
 * Run with:
 *   ./normalized_drift_chambers_comparison [Nevents] [dataFile] [mcFile] [torus]
 *
 *   torus: +1 for outbending, -1 for inbending
 *
 * What it does:
 *   1) Creates output/ and output/dc/ if they do not exist.
 *   2) Parses torus to apply DC fiducial cuts.
 *   3) Loops over data and MC, filling per-region 2D histograms for
 *      electrons (PID=11) and protons (PID=2212), applying fiducial cuts.
 *   4) Normalizes each histogram and computes data/MC ratio.
 *   5) Draws “COLZ” canvases for Data, MC, Ratio (with mean/std in legend).
 *   6) Draws outlier maps with 4-color palette:
 *        • 0.5–2     → blue
 *        • 1/3–0.5 & 2–3   → light red (orange)
 *        • 1/5–1/3 & 3–5   → dark red
 *        • <1/5 or >5 → pink
 *****************************************************************************/

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <utility>
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TLegend.h"

// ----------------------------------------------------------------------------
// DC fiducial cut based on torus polarity and track theta (deg)
// ----------------------------------------------------------------------------
bool dc_fiducial_cut(int torus, double theta_deg,
                     double edge1, double edge2, double edge3) {
    bool inbending  = (torus < 0);
    bool outbending = (torus > 0);
    if (inbending) {
        if (theta_deg < 10.0) {
            return edge1 > 10.0 && edge2 > 10.0 && edge3 > 10.0;
        } else {
            return edge1 > 3.0 && edge2 > 3.0 && edge3 > 10.0;
        }
    }
    else if (outbending) {
        return edge1 > 3.0 && edge2 > 3.0 && edge3 > 10.0;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Ensure same color scale across three TH2D histograms
// ----------------------------------------------------------------------------
void SetSame2DScale(TH2D* h1, TH2D* h2, TH2D* h3) {
    double mn =  1e9, mx = -1e9;
    for (auto* h : {h1,h2,h3}) {
        mn = std::min(mn, h->GetMinimum());
        mx = std::max(mx, h->GetMaximum());
    }
    for (auto* h : {h1,h2,h3}) {
        h->SetMinimum(mn);
        h->SetMaximum(mx);
    }
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " [Nevents] [dataFile] [mcFile] [torus]\n"
                  << "  torus: +1 for outbending, -1 for inbending" << std::endl;
        return 1;
    }

    // create output directories
    gSystem->mkdir("output",    kTRUE);
    gSystem->mkdir("output/dc", kTRUE);

    // parse args
    Long64_t maxEvents = std::stoll(argv[1]);
    if (maxEvents == 0) maxEvents = -1;
    std::string dataFile = argv[2];
    std::string mcFile   = argv[3];
    int torus = std::stoi(argv[4]);

    // setup chains
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    dataCh.Add(dataFile.c_str());
    mcCh .Add(mcFile  .c_str());

    // branches
    Int_t    pid;
    Double_t x6,y6,x18,y18,x36,y36;
    Double_t edge6,edge18,edge36;
    Double_t theta;
    dataCh.SetBranchAddress("particle_pid", &pid);
    dataCh.SetBranchAddress("traj_x_6",     &x6);
    dataCh.SetBranchAddress("traj_y_6",     &y6);
    dataCh.SetBranchAddress("traj_x_18",    &x18);
    dataCh.SetBranchAddress("traj_y_18",    &y18);
    dataCh.SetBranchAddress("traj_x_36",    &x36);
    dataCh.SetBranchAddress("traj_y_36",    &y36);
    dataCh.SetBranchAddress("traj_edge_6",  &edge6);
    dataCh.SetBranchAddress("traj_edge_18", &edge18);
    dataCh.SetBranchAddress("traj_edge_36", &edge36);
    dataCh.SetBranchAddress("theta",        &theta);
    mcCh .SetBranchAddress("particle_pid", &pid);
    mcCh .SetBranchAddress("traj_x_6",     &x6);
    mcCh .SetBranchAddress("traj_y_6",     &y6);
    mcCh .SetBranchAddress("traj_x_18",    &x18);
    mcCh .SetBranchAddress("traj_y_18",    &y18);
    mcCh .SetBranchAddress("traj_x_36",    &x36);
    mcCh .SetBranchAddress("traj_y_36",    &y36);
    mcCh .SetBranchAddress("traj_edge_6",  &edge6);
    mcCh .SetBranchAddress("traj_edge_18", &edge18);
    mcCh .SetBranchAddress("traj_edge_36", &edge36);
    mcCh .SetBranchAddress("theta",        &theta);

    // constants
    const int NB2=200;
    const double xmins[3]={-180,-280,-450}, xmaxs[3]={180,280,450};
    std::vector<std::pair<int,std::string>> species = {{11,"electron"},{2212,"proton"}};

    // book histograms
    std::vector<std::vector<TH2D*>> hData(species.size()), hMC(species.size()), hRatio(species.size());
    for(size_t s=0;s<species.size();++s){
        auto &lab=species[s].second;
        for(int r=0;r<3;++r){
            hData[s].push_back(new TH2D(Form("hD_%s_r%d", lab.c_str(), r+1),
                Form("%s Data R%d; x; y", lab.c_str(), r+1),
                NB2, xmins[r], xmaxs[r], NB2, xmins[r], xmaxs[r]));
            hMC  [s].push_back((TH2D*)hData[s][r]->Clone(Form("hM_%s_r%d", lab.c_str(), r+1)));
        }
    }

    // fill data
    Long64_t nD=dataCh.GetEntries(); if(maxEvents>0 && maxEvents<nD) nD=maxEvents;
    for(Long64_t i=0;i<nD;++i){ dataCh.GetEntry(i);
        if (!dc_fiducial_cut(torus, theta*180.0/M_PI, edge6, edge18, edge36)) continue;
        for(size_t s=0;s<species.size();++s){
            if(pid!=species[s].first) continue;
            double xs[3]={x6,x18,x36}, ys[3]={y6,y18,y36};
            for(int r=0;r<3;++r) if(xs[r]!=-9999) hData[s][r]->Fill(xs[r], ys[r]);
        }
    }

    // fill MC
    Long64_t nM=mcCh.GetEntries(); if(maxEvents>0 && maxEvents<nM) nM=maxEvents;
    for(Long64_t i=0;i<nM;++i){ mcCh.GetEntry(i);
        if (!dc_fiducial_cut(torus, theta*180.0/M_PI, edge6, edge18, edge36)) continue;
        for(size_t s=0;s<species.size();++s){
            if(pid!=species[s].first) continue;
            double xs[3]={x6,x18,x36}, ys[3]={y6,y18,y36};
            for(int r=0;r<3;++r) if(xs[r]!=-9999) hMC[s][r]->Fill(xs[r], ys[r]);
        }
    }

    // normalize, ratio, stats
    std::vector<std::vector<double>> meanV(species.size(), std::vector<double>(3)), sigmaV=meanV;
    for(size_t s=0;s<species.size();++s){
        for(int r=0;r<3;++r){
            double id=hData[s][r]->Integral(); if(id>0) hData[s][r]->Scale(1./id);
            double im=hMC  [s][r]->Integral(); if(im>0) hMC  [s][r]->Scale(1./im);
            hRatio[s].push_back((TH2D*)hData[s][r]->Clone(Form("hR_%s_r%d", species[s].second.c_str(), r+1)));
            hRatio[s][r]->Divide(hMC[s][r]);
            int cnt=0; double sum=0,sum2=0;
            for(int ix=1;ix<=NB2;++ix) for(int iy=1;iy<=NB2;++iy){
                double v=hRatio[s][r]->GetBinContent(ix,iy);
                if(v<=0) continue;
                sum+=v; sum2+=v*v; cnt++;
            }
            if(cnt>0){
                meanV[s][r]=sum/cnt;
                sigmaV[s][r]=std::sqrt(sum2/cnt - meanV[s][r]*meanV[s][r]);
            }
        }
    }

    gStyle->SetOptStat(0);
    // draw Data, MC, Ratio
    for(size_t s=0;s<species.size();++s){
        auto &lab=species[s].second;

        // Data
        TCanvas cD(Form("cD_%s", lab.c_str()), Form("%s Data", lab.c_str()), 1800,600);
        cD.Divide(3,1);
        SetSame2DScale(hData[s][0],hData[s][1],hData[s][2]);
        for(int r=0;r<3;++r){
            cD.cd(r+1);
            gPad->SetLeftMargin(.15);
            gPad->SetRightMargin(.15);
            gPad->SetLogz();
            hData[s][r]->Draw("COLZ");
        }
        cD.SaveAs(Form("output/dc/data_%s.png", lab.c_str()));

        // MC
        TCanvas cM(Form("cM_%s", lab.c_str()), Form("%s MC", lab.c_str()), 1800,600);
        cM.Divide(3,1);
        SetSame2DScale(hMC[s][0],hMC[s][1],hMC[s][2]);
        for(int r=0;r<3;++r){
            cM.cd(r+1);
            gPad->SetLeftMargin(.15);
            gPad->SetRightMargin(.15);
            gPad->SetLogz();
            hMC[s][r]->Draw("COLZ");
        }
        cM.SaveAs(Form("output/dc/mc_%s.png", lab.c_str()));

        // Ratio
        TCanvas cR(Form("cR_%s", lab.c_str()), Form("%s Data/MC Ratio", lab.c_str()), 1800,600);
        cR.Divide(3,1);
        SetSame2DScale(hRatio[s][0],hRatio[s][1],hRatio[s][2]);
        for(int r=0;r<3;++r){
            cR.cd(r+1);
            gPad->SetLeftMargin(.15);
            gPad->SetRightMargin(.15);
            gPad->SetLogz();
            hRatio[s][r]->Draw("COLZ");
            TLegend leg(0.6,0.7,0.9,0.9);
            leg.SetFillColor(kWhite);
            leg.SetBorderSize(1);
            leg.SetTextSize(0.03);
            leg.AddEntry((TObject*)0, Form("Mean=%.3f", meanV[s][r]), "");
            leg.AddEntry((TObject*)0, Form("StdDev=%.3f", sigmaV[s][r]), "");
            leg.Draw();
        }
        cR.SaveAs(Form("output/dc/ratio_%s.png", lab.c_str()));
    }

    // four-color palette for outliers
    Int_t pal[4] = { kBlue, kOrange, kRed, kPink };
    gStyle->SetPalette(4, pal);

    // draw outlier maps with new scaling
    for(size_t s=0;s<species.size();++s){
        auto &lab=species[s].second;
        TCanvas cO(Form("cO_%s", lab.c_str()), Form("%s Ratio Outliers", lab.c_str()), 1800,600);
        cO.Divide(3,1);
        for(int r=0;r<3;++r){
            cO.cd(r+1);
            gPad->SetLeftMargin(.15);
            gPad->SetRightMargin(.15);
            gPad->SetLogz(0);

            // clone & reset for map
            TH2D* m = (TH2D*)hRatio[s][r]->Clone(Form("m_%s_r%d", lab.c_str(), r+1));
            m->Reset();

            // fill with category integers
            for(int ix=1; ix<=NB2; ++ix){
                for(int iy=1; iy<=NB2; ++iy){
                    double v = hRatio[s][r]->GetBinContent(ix, iy);
                    if(v <= 0) continue;
                    int binCat = 4; // default: <1/5 or >5 (pink)
                    if(v >= 0.5 && v <= 2.0) {
                        binCat = 1; // blue
                    }
                    else if((v >= 1.0/3.0 && v < 0.5) || (v > 2.0   && v <= 3.0)) {
                        binCat = 2; // light red (orange)
                    }
                    else if((v >= 1.0/5.0 && v < 1.0/3.0) || (v > 3.0 && v <= 5.0)) {
                        binCat = 3; // dark red
                    }
                    m->SetBinContent(ix, iy, binCat);
                }
            }

            // set discrete contour levels
            m->SetContour(4);
            m->SetContourLevel(0, 1);
            m->SetContourLevel(1, 2);
            m->SetContourLevel(2, 3);
            m->SetContourLevel(3, 4);

            // draw
            m->Draw("COLZ");
        }
        cO.SaveAs(Form("output/dc/outliers_%s.png", lab.c_str()));
    }

    return 0;
}