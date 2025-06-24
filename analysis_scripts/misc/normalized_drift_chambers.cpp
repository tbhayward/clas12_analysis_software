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
 *   2) Parses torus to apply original and strict DC fiducial cuts.
 *   3) Loops over data and MC, filling per-region 2D histograms for
 *      electrons (PID=11) and protons (PID=2212):
 *        • no cuts
 *        • original cuts
 *        • strict cuts
 *   4) Normalizes each histogram and computes data/MC ratio for
 *      both original- and strict-cut sets.
 *   5) Draws “COLZ” canvases for Data, MC, and two ratio versions
 *      (original and strict) with mean/stddev and survival-percentage legends.
 *   6) Draws outlier maps (4-color) for both ratio versions.
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
// Original DC fiducial cut based on torus polarity and track theta (deg)
// ----------------------------------------------------------------------------
bool dc_fiducial_cut(int torus, double theta_deg,
                     double edge1, double edge2, double edge3) {
    bool inbending  = (torus < 0);
    bool outbending = (torus > 0);
    if (inbending) {
        if (theta_deg < 10.0)
            return edge1 > 10.0 && edge2 > 10.0 && edge3 > 10.0;
        else
            return edge1 > 3.0 && edge2 > 3.0 && edge3 > 10.0;
    }
    else if (outbending) {
        return edge1 > 3.0 && edge2 > 3.0 && edge3 > 10.0;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Strict DC fiducial cut: adds radial restrictions per region
// ----------------------------------------------------------------------------
bool dc_fiducial_cut_strict(int torus, int pid, double theta_deg,
                            double edge1, double edge2, double edge3,
                            double x6, double x18, double x36,
                            double y6, double y18, double y36) {
    double r6  = std::hypot(x6,  y6);
    double r18 = std::hypot(x18, y18);
    double r36 = std::hypot(x36, y36);

    // inbending, electron
    if (torus < 0 && pid == 11) {
        bool pass = (theta_deg < 10.0
                     ? (edge1>10 && edge2>10 && edge3>14)
                     : (edge1>3  && edge2>3  && edge3>14));
        if (!pass) return false;
        if (r6<50 || r18<70 || r36<60) return false;
        return true;
    }
    // inbending, proton
    if (torus < 0 && pid == 2212) {
        bool pass = (theta_deg < 10.0
                     ? (edge1>10 && edge2>10 && edge3>10)
                     : (edge1>3  && edge2>3  && edge3>10));
        if (!pass) return false;
        if (r6<60  || r6>150  ||
            r18<120|| r18>210 ||
            r36<260|| r36>380 ) return false;
        return true;
    }
    // outbending, electron
    if (torus > 0 && pid == 11) {
        bool pass = edge1>3 && edge2>3 && edge3>14;
        if (!pass) return false;
        if (r6<35  || r6>60  ||
            r18<78 || r18>98 ||
            r36<100|| r36>180) return false;
        return true;
    }
    // outbending, proton
    if (torus > 0 && pid == 2212) {
        bool pass = edge1>3 && edge2>3 && edge3>14;
        if (!pass) return false;
        if (r6<85  || r6>145 ||
            r18<120|| r18>190 ||
            r36<90 || r36>250) return false;
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Ensure same color scale across three TH2D histograms
// ----------------------------------------------------------------------------
void SetSame2DScale(TH2D* h1, TH2D* h2, TH2D* h3) {
    double mn = 1e9, mx = -1e9;
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
                  << "  torus: +1 outbending, -1 inbending\n";
        return 1;
    }

    gSystem->mkdir("output",    kTRUE);
    gSystem->mkdir("output/dc", kTRUE);

    Long64_t maxEvents = std::stoll(argv[1]);
    if (maxEvents == 0) maxEvents = -1;
    std::string dataFile = argv[2], mcFile = argv[3];
    int torus = std::stoi(argv[4]);

    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    dataCh.Add(dataFile.c_str());
    mcCh .Add(mcFile .c_str());

    Int_t    pid;
    Double_t x6,y6,x18,y18,x36,y36;
    Double_t edge6,edge18,edge36, theta;
    dataCh.SetBranchAddress("particle_pid",&pid);
    dataCh.SetBranchAddress("traj_x_6",    &x6);
    dataCh.SetBranchAddress("traj_y_6",    &y6);
    dataCh.SetBranchAddress("traj_x_18",   &x18);
    dataCh.SetBranchAddress("traj_y_18",   &y18);
    dataCh.SetBranchAddress("traj_x_36",   &x36);
    dataCh.SetBranchAddress("traj_y_36",   &y36);
    dataCh.SetBranchAddress("traj_edge_6", &edge6);
    dataCh.SetBranchAddress("traj_edge_18",&edge18);
    dataCh.SetBranchAddress("traj_edge_36",&edge36);
    dataCh.SetBranchAddress("theta",       &theta);
    mcCh .SetBranchAddress("particle_pid",&pid);
    mcCh .SetBranchAddress("traj_x_6",    &x6);
    mcCh .SetBranchAddress("traj_y_6",    &y6);
    mcCh .SetBranchAddress("traj_x_18",   &x18);
    mcCh .SetBranchAddress("traj_y_18",   &y18);
    mcCh .SetBranchAddress("traj_x_36",   &x36);
    mcCh .SetBranchAddress("traj_y_36",   &y36);
    mcCh .SetBranchAddress("traj_edge_6", &edge6);
    mcCh .SetBranchAddress("traj_edge_18",&edge18);
    mcCh .SetBranchAddress("traj_edge_36",&edge36);
    mcCh .SetBranchAddress("theta",       &theta);

    const int NB2 = 200;
    const double xmins[3] = {-180,-280,-450}, xmaxs[3] = {180,280,450};
    std::vector<std::pair<int,std::string>> species = {{11,"electron"},{2212,"proton"}};

    std::vector<std::vector<TH2D*>> 
      hDataOrig(species.size()), hMCOrig(species.size()), hRatioOrig(species.size()),
      hDataStr (species.size()), hMCStr  (species.size()), hRatioStr (species.size());

    for (size_t s=0; s<species.size(); ++s) {
      auto &lab=species[s].second;
      for (int r=0; r<3; ++r) {
        hDataOrig[s].push_back(new TH2D(
          Form("hD_orig_%s_r%d", lab.c_str(),r+1),
          Form("%s Data Orig R%d; x; y",lab.c_str(),r+1),
          NB2,xmins[r],xmaxs[r],NB2,xmins[r],xmaxs[r]));
        hMCOrig  [s].push_back((TH2D*)hDataOrig[s][r]->Clone(
          Form("hM_orig_%s_r%d", lab.c_str(),r+1)));
        hDataStr [s].push_back(new TH2D(
          Form("hD_strict_%s_r%d", lab.c_str(),r+1),
          Form("%s Data Strict R%d; x; y",lab.c_str(),r+1),
          NB2,xmins[r],xmaxs[r],NB2,xmins[r],xmaxs[r]));
        hMCStr   [s].push_back((TH2D*)hDataStr[s][r]->Clone(
          Form("hM_strict_%s_r%d", lab.c_str(),r+1)));
      }
    }

    std::vector<std::vector<Long64_t>> 
      nNoCut(species.size(), std::vector<Long64_t>(3,0)),
      nOrig (species.size(), std::vector<Long64_t>(3,0)),
      nStr  (species.size(), std::vector<Long64_t>(3,0));

    Long64_t nD = dataCh.GetEntries();
    if (maxEvents>0 && maxEvents<nD) nD=maxEvents;
    for (Long64_t i=0;i<nD;++i) {
      dataCh.GetEntry(i);
      for (size_t s=0;s<species.size();++s) {
        if (pid!=species[s].first) continue;
        double xs[3]={x6,x18,x36}, ys[3]={y6,y18,y36};
        for (int r=0;r<3;++r) {
          if (xs[r]==-9999||ys[r]==-9999) continue;
          nNoCut[s][r]++;
          if (dc_fiducial_cut(torus, theta*180./M_PI, edge6,edge18,edge36)) {
            hDataOrig[s][r]->Fill(xs[r], ys[r]);
            nOrig[s][r]++;
          }
          if (dc_fiducial_cut_strict(torus, pid, theta*180./M_PI,
                                     edge6,edge18,edge36,
                                     x6,x18,x36, y6,y18,y36)) {
            hDataStr[s][r]->Fill(xs[r], ys[r]);
            nStr[s][r]++;
          }
        }
      }
    }

    Long64_t nM = mcCh.GetEntries();
    if (maxEvents>0 && maxEvents<nM) nM=maxEvents;
    for (Long64_t i=0;i<nM;++i) {
      mcCh.GetEntry(i);
      for (size_t s=0;s<species.size();++s) {
        if (pid!=species[s].first) continue;
        double xs[3]={x6,x18,x36}, ys[3]={y6,y18,y36};
        for (int r=0;r<3;++r) {
          if (xs[r]==-9999||ys[r]==-9999) continue;
          if (dc_fiducial_cut(torus, theta*180./M_PI, edge6,edge18,edge36))
            hMCOrig[s][r]->Fill(xs[r], ys[r]);
          if (dc_fiducial_cut_strict(torus, pid, theta*180./M_PI,
                                     edge6,edge18,edge36,
                                     x6,x18,x36, y6,y18,y36))
            hMCStr  [s][r]->Fill(xs[r], ys[r]);
        }
      }
    }

    std::vector<std::vector<double>> meanOrig(species.size(), std::vector<double>(3)),
                                sigmaOrig(meanOrig),
                                meanStr (meanOrig),
                                sigmaStr(meanOrig);

    for (size_t s=0;s<species.size();++s) {
      for (int r=0;r<3;++r) {
        double dI = hDataOrig[s][r]->Integral();
        if (dI>0) hDataOrig[s][r]->Scale(1./dI);
        double mI = hMCOrig  [s][r]->Integral();
        if (mI>0) hMCOrig[s][r]->Scale(1./mI);
        hRatioOrig[s].push_back((TH2D*)hDataOrig[s][r]->Clone(
          Form("hR_orig_%s_r%d", species[s].second.c_str(), r+1)
        ));
        hRatioOrig[s][r]->Divide(hMCOrig[s][r]);
        int cnt=0; double sum=0,sum2=0;
        for (int ix=1; ix<=NB2; ++ix) for (int iy=1; iy<=NB2; ++iy){
          double v=hRatioOrig[s][r]->GetBinContent(ix,iy);
          if (v>0) { sum+=v; sum2+=v*v; cnt++; }
        }
        if (cnt>0) {
          meanOrig[s][r]  = sum/cnt;
          sigmaOrig[s][r] = std::sqrt(sum2/cnt - meanOrig[s][r]*meanOrig[s][r]);
        }
        double dI2 = hDataStr[s][r]->Integral();
        if (dI2>0) hDataStr[s][r]->Scale(1./dI2);
        double mI2 = hMCStr  [s][r]->Integral();
        if (mI2>0) hMCStr[s][r]->Scale(1./mI2);
        hRatioStr[s].push_back((TH2D*)hDataStr[s][r]->Clone(
          Form("hR_strict_%s_r%d", species[s].second.c_str(), r+1)
        ));
        hRatioStr[s][r]->Divide(hMCStr[s][r]);
        cnt=0; sum=0; sum2=0;
        for (int ix=1; ix<=NB2; ++ix) for (int iy=1; iy<=NB2; ++iy){
          double v=hRatioStr[s][r]->GetBinContent(ix,iy);
          if (v>0) { sum+=v; sum2+=v*v; cnt++; }
        }
        if (cnt>0) {
          meanStr [s][r]  = sum/cnt;
          sigmaStr[s][r] = std::sqrt(sum2/cnt - meanStr[s][r]*meanStr[s][r]);
        }
      }
    }

    gStyle->SetOptStat(0);

    // draw original ratio
    for (size_t s=0;s<species.size();++s) {
      auto &lab=species[s].second;
      TCanvas cR(Form("cR_orig_%s", lab.c_str()),
                 Form("%s Ratio Orig", lab.c_str()), 1800,600);
      cR.Divide(3,1);
      SetSame2DScale(hRatioOrig[s][0],hRatioOrig[s][1],hRatioOrig[s][2]);
      for (int r=0;r<3;++r) {
        cR.cd(r+1);
        gPad->SetLeftMargin(.15); gPad->SetRightMargin(.15); gPad->SetLogz();
        hRatioOrig[s][r]->Draw("COLZ");
        gPad->Update();
        double pO = nNoCut[s][r] ? 100.0*nOrig[s][r]/nNoCut[s][r] : 0;
        double pS = nNoCut[s][r] ? 100.0*nStr [s][r]/nNoCut[s][r] : 0;
        TLegend legF(0.15,0.75,0.45,0.90);
        legF.SetFillColor(kWhite); legF.SetBorderSize(1); legF.SetTextSize(0.03);
        legF.AddEntry((TObject*)0,"No cuts: 100%","");  
        legF.AddEntry((TObject*)0,Form("Orig cuts: %.1f%%",pO),"");  
        legF.AddEntry((TObject*)0,Form("Strict cuts: %.1f%%",pS),"");  
        legF.Draw();
        TLegend legS(0.60,0.75,0.90,0.90);
        legS.SetFillColor(kWhite); legS.SetBorderSize(1); legS.SetTextSize(0.03);
        legS.AddEntry((TObject*)0,Form("Mean=%.3f", meanOrig[s][r]),"");
        legS.AddEntry((TObject*)0,Form("Std=%.3f",  sigmaOrig[s][r]),"");
        legS.Draw();
      }
      cR.SaveAs(Form("output/dc/ratio_orig_%s.png", lab.c_str()));
    }

    // draw strict ratio
    for (size_t s=0;s<species.size();++s) {
      auto &lab=species[s].second;
      TCanvas cR(Form("cR_strict_%s", lab.c_str()),
                 Form("%s Ratio Strict", lab.c_str()), 1800,600);
      cR.Divide(3,1);
      SetSame2DScale(hRatioStr[s][0],hRatioStr[s][1],hRatioStr[s][2]);
      for (int r=0;r<3;++r) {
        cR.cd(r+1);
        gPad->SetLeftMargin(.15); gPad->SetRightMargin(.15); gPad->SetLogz();
        hRatioStr[s][r]->Draw("COLZ");
        gPad->Update();
        double pO = nNoCut[s][r] ? 100.0*nOrig[s][r]/nNoCut[s][r] : 0;
        double pS = nNoCut[s][r] ? 100.0*nStr [s][r]/nNoCut[s][r] : 0;
        TLegend legF(0.15,0.75,0.45,0.90);
        legF.SetFillColor(kWhite); legF.SetBorderSize(1); legF.SetTextSize(0.03);
        legF.AddEntry((TObject*)0,"No cuts: 100%","");
        legF.AddEntry((TObject*)0,Form("Orig cuts: %.1f%%",pO),"");
        legF.AddEntry((TObject*)0,Form("Strict cuts: %.1f%%",pS),"");
        legF.Draw();
        TLegend legS(0.60,0.75,0.90,0.90);
        legS.SetFillColor(kWhite); legS.SetBorderSize(1); legS.SetTextSize(0.03);
        legS.AddEntry((TObject*)0,Form("Mean=%.3f", meanStr [s][r]),"");
        legS.AddEntry((TObject*)0,Form("Std=%.3f",  sigmaStr[s][r]),"");
        legS.Draw();
      }
      cR.SaveAs(Form("output/dc/ratio_strict_%s.png", lab.c_str()));
    }

    // four-color palette for outliers
    Int_t pal[4] = { kBlue, kOrange, kRed, kPink };
    gStyle->SetPalette(4,pal);

    // draw outlier maps for orig and strict
    for (size_t s=0;s<species.size();++s) {
      auto &lab=species[s].second;
      // orig outliers
      {
        TCanvas cO(Form("cO_orig_%s", lab.c_str()),
                   Form("%s Outliers Orig", lab.c_str()), 1800,600);
        cO.Divide(3,1);
        for (int r=0;r<3;++r) {
          cO.cd(r+1);
          gPad->SetLeftMargin(.15); gPad->SetRightMargin(.15); gPad->SetLogz(0);
          TH2D* m=(TH2D*)hRatioOrig[s][r]->Clone();
          m->Reset();
          for (int ix=1; ix<=NB2; ++ix) for (int iy=1; iy<=NB2; ++iy) {
            double v=hRatioOrig[s][r]->GetBinContent(ix,iy);
            if (v<=0) continue;
            int cat=4;
            if      (v>=0.5&&v<=2.0)                          cat=1;
            else if ((v>=1.0/3.0&&v<0.5)||(v>2.0&&v<=3.0))     cat=2;
            else if ((v>=1.0/5.0&&v<1.0/3.0)||(v>3.0&&v<=5.0)) cat=3;
            m->SetBinContent(ix,iy,cat);
          }
          m->SetContour(4);
          for(int i=0;i<4;++i) m->SetContourLevel(i,i+1);
          m->Draw("COLZ");
        }
        cO.SaveAs(Form("output/dc/outliers_orig_%s.png", lab.c_str()));
      }
      // strict outliers
      {
        TCanvas cO(Form("cO_strict_%s", lab.c_str()),
                   Form("%s Outliers Strict", lab.c_str()), 1800,600);
        cO.Divide(3,1);
        for (int r=0;r<3;++r) {
          cO.cd(r+1);
          gPad->SetLeftMargin(.15); gPad->SetRightMargin(.15); gPad->SetLogz(0);
          TH2D* m=(TH2D*)hRatioStr[s][r]->Clone();
          m->Reset();
          for (int ix=1; ix<=NB2; ++ix) for (int iy=1; iy<=NB2; ++iy) {
            double v=hRatioStr[s][r]->GetBinContent(ix,iy);
            if (v<=0) continue;
            int cat=4;
            if      (v>=0.5&&v<=2.0)                          cat=1;
            else if ((v>=1.0/3.0&&v<0.5)||(v>2.0&&v<=3.0))     cat=2;
            else if ((v>=1.0/5.0&&v<1.0/3.0)||(v>3.0&&v<=5.0)) cat=3;
            m->SetBinContent(ix,iy,cat);
          }
          m->SetContour(4);
          for(int i=0;i<4;++i) m->SetContourLevel(i,i+1);
          m->Draw("COLZ");
        }
        cO.SaveAs(Form("output/dc/outliers_strict_%s.png", lab.c_str()));
      }
    }

    return 0;
}