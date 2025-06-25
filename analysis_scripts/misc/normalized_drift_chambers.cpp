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
 *   4) Prints survival rates to console for electrons and protons.
 *   5) Normalizes each histogram and computes data/MC ratio for
 *      both original- and strict-cut sets.
 *   6) Draws “COLZ” canvases for:
 *        • Data (original cuts)
 *        • MC   (original cuts)
 *        • Data/MC ratio (original vs strict) with survival% in title
 *        • Data/MC ratio strict  with survival% in title
 *   7) Draws outlier maps (4-color) for both ratio versions.
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
    bool inbending  = (torus < 0), outbending = (torus > 0);
    if (inbending) {
        if (theta_deg < 10.0)
            return edge1>10 && edge2>10 && edge3>10;
        else
            return edge1>3 && edge2>3 && edge3>10;
    } else if (outbending) {
        return edge1>3 && edge2>3 && edge3>10;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Strict DC fiducial cut with additional radial constraints
// ----------------------------------------------------------------------------
bool dc_fiducial_cut_strict(int torus, int pid, double theta_deg,
                            double edge1, double edge2, double edge3,
                            double x6, double x18, double x36,
                            double y6, double y18, double y36) {
    double r6  = std::hypot(x6,  y6);
    double r18 = std::hypot(x18, y18);
    double r36 = std::hypot(x36, y36);

    if (torus < 0 && pid == 11) {
        bool pass = (theta_deg < 10.0
                     ? (edge1>14 && edge2>14 && edge3>14)
                     : (edge1>6  && edge2>6  && edge3>14));
        if (!pass || r6<50 || r18<70 || r36<70 || r36>260) return false;
        return true;
    }
    if (torus < 0 && pid == 2212) {
        bool pass = (theta_deg < 10.0
                     ? (edge1>10 && edge2>10 && edge3>10)
                     : (edge1>3  && edge2>3  && edge3>10));
        if (!pass || r6<60 || r6>150 ||
            r18<120|| r18>210 ||
            r36<280|| r36>380 ) return false;
        return true;
    }
    if (torus > 0 && pid == 11) {
        bool pass = (edge1>3 && edge2>3 && edge3>14);
        if (!pass || r6<35 || r6>60 ||
            r18<78|| r18>98 ||
            r36<100||r36>180) return false;
        return true;
    }
    if (torus > 0 && pid == 2212) {
        bool pass = (edge1>3 && edge2>3 && edge3>14);
        if (!pass || r6<85 || r6>145 ||
            r18<120|| r18>190 ||
            r36<90|| r36>250) return false;
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Match color scales across three TH2Ds
// ----------------------------------------------------------------------------
void SetSame2DScale(TH2D* a, TH2D* b, TH2D* c) {
    double mn = 1e9, mx = -1e9;
    for (auto*h : {a,b,c}) {
        mn = std::min(mn, h->GetMinimum());
        mx = std::max(mx, h->GetMaximum());
    }
    for (auto*h : {a,b,c}) {
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

    // create output dirs
    gSystem->mkdir("output",    kTRUE);
    gSystem->mkdir("output/dc", kTRUE);

    // parse args
    Long64_t maxE = std::stoll(argv[1]);
    if (maxE == 0) maxE = -1;
    std::string dataFile = argv[2];
    std::string mcFile   = argv[3];
    int torus = std::stoi(argv[4]);

    // chains
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    dataCh.Add(dataFile.c_str());
    mcCh.Add(mcFile.c_str());

    // branches
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

    // constants
    const int NB2 = 200;
    const double xmins[3] = {-180,-280,-450}, xmaxs[3] = {180,280,450};
    std::vector<std::pair<int,std::string>> species = {{11,"electron"},{2212,"proton"}};

    // histograms
    std::vector<std::vector<TH2D*>> 
      Dorig(species.size()), Morig(species.size()), Rorig(species.size()),
      Dstr (species.size()), Mstr (species.size()), Rstr (species.size());

    for (size_t s=0; s<species.size(); ++s) {
        auto &lab = species[s].second;
        for (int r=0; r<3; ++r) {
            Dorig[s].push_back(new TH2D(
              Form("Dorig_%s_r%d", lab.c_str(), r+1),
              Form("%s Data Orig R%d; x; y", lab.c_str(), r+1),
              NB2, xmins[r], xmaxs[r], NB2, xmins[r], xmaxs[r]));
            Morig[s].push_back((TH2D*)Dorig[s][r]->Clone(
              Form("Morig_%s_r%d", lab.c_str(), r+1)));
            Dstr[s].push_back(new TH2D(
              Form("Dstr_%s_r%d", lab.c_str(), r+1),
              Form("%s Data Strict R%d; x; y", lab.c_str(), r+1),
              NB2, xmins[r], xmaxs[r], NB2, xmins[r], xmaxs[r]));
            Mstr[s].push_back((TH2D*)Dstr[s][r]->Clone(
              Form("Mstr_%s_r%d", lab.c_str(), r+1)));
        }
    }

    // survival counters
    std::vector<std::vector<Long64_t>> 
      nNo(species.size(), std::vector<Long64_t>(3,0)),
      nOr(species.size(), std::vector<Long64_t>(3,0)),
      nSt(species.size(), std::vector<Long64_t>(3,0));

    // fill data
    Long64_t nD = dataCh.GetEntries();
    if (maxE>0 && maxE<nD) nD = maxE;
    for (Long64_t i=0; i<nD; ++i) {
        dataCh.GetEntry(i);
        for (size_t s=0; s<species.size(); ++s) {
            if (pid != species[s].first) continue;
            double xs[3] = {x6, x18, x36}, ys[3] = {y6, y18, y36};
            for (int r=0; r<3; ++r) {
                if (xs[r]==-9999 || ys[r]==-9999) continue;
                nNo[s][r]++;
                if (dc_fiducial_cut(torus, theta*180./M_PI, edge6, edge18, edge36)) {
                    Dorig[s][r]->Fill(xs[r], ys[r]);
                    nOr[s][r]++;
                }
                if (dc_fiducial_cut_strict(torus, pid, theta*180./M_PI,
                                           edge6, edge18, edge36,
                                           x6, x18, x36, y6, y18, y36)) {
                    Dstr[s][r]->Fill(xs[r], ys[r]);
                    nSt[s][r]++;
                }
            }
        }
    }

    // print survival rates
    std::cout << "=== Survival rates (data) ===\n";
    for (size_t s=0; s<species.size(); ++s) {
        const auto &lab = species[s].second;
        Long64_t totNo=0, totOr=0, totSt=0;
        for (int r=0; r<3; ++r) {
            totNo += nNo[s][r];
            totOr += nOr[s][r];
            totSt += nSt[s][r];
        }
        double pOr = totNo ? 100.0*totOr/totNo : 0.0;
        double pSt = totNo ? 100.0*totSt/totNo : 0.0;
        std::cout << lab << ": original = " << pOr << "%, strict = " << pSt << "%\n";
    }
    std::cout << "=============================\n";

    // fill MC
    Long64_t nM = mcCh.GetEntries();
    if (maxE>0 && maxE<nM) nM = maxE;
    for (Long64_t i=0; i<nM; ++i) {
        mcCh.GetEntry(i);
        for (size_t s=0; s<species.size(); ++s) {
            if (pid != species[s].first) continue;
            double xs[3] = {x6, x18, x36}, ys[3] = {y6, y18, y36};
            for (int r=0; r<3; ++r) {
                if (xs[r]==-9999 || ys[r]==-9999) continue;
                if (dc_fiducial_cut(torus, theta*180./M_PI, edge6, edge18, edge36))
                    Morig[s][r]->Fill(xs[r], ys[r]);
                if (dc_fiducial_cut_strict(torus, pid, theta*180./M_PI,
                                           edge6, edge18, edge36,
                                           x6, x18, x36, y6, y18, y36))
                    Mstr[s][r]->Fill(xs[r], ys[r]);
            }
        }
    }

    // normalize and build ratios & stats
    std::vector<std::vector<double>> mO(species.size(), std::vector<double>(3)),
                                sO(mO), mS(mO), sS(mO);
    for (size_t s=0; s<species.size(); ++s) {
        for (int r=0; r<3; ++r) {
            // original
            double dI = Dorig[s][r]->Integral();
            if (dI>0) Dorig[s][r]->Scale(1./dI);
            double mI = Morig[s][r]->Integral();
            if (mI>0) Morig[s][r]->Scale(1./mI);
            Rorig[s].push_back((TH2D*)Dorig[s][r]->Clone(
                Form("Rorig_%s_r%d", species[s].second.c_str(), r+1)));
            Rorig[s][r]->Divide(Morig[s][r]);
            int cnt=0; double sum=0, sum2=0;
            for (int ix=1; ix<=NB2; ++ix) for (int iy=1; iy<=NB2; ++iy) {
                double v = Rorig[s][r]->GetBinContent(ix,iy);
                if (v>0) { sum+=v; sum2+=v*v; cnt++; }
            }
            if (cnt>0) {
                mO[s][r] = sum/cnt;
                sO[s][r] = std::sqrt(sum2/cnt - mO[s][r]*mO[s][r]);
            }
            // strict
            double dI2 = Dstr[s][r]->Integral();
            if (dI2>0) Dstr[s][r]->Scale(1./dI2);
            double mI2 = Mstr[s][r]->Integral();
            if (mI2>0) Mstr[s][r]->Scale(1./mI2);
            Rstr[s].push_back((TH2D*)Dstr[s][r]->Clone(
                Form("Rstr_%s_r%d", species[s].second.c_str(), r+1)));
            Rstr[s][r]->Divide(Mstr[s][r]);
            cnt=0; sum=0; sum2=0;
            for (int ix=1; ix<=NB2; ++ix) for (int iy=1; iy<=NB2; ++iy) {
                double v = Rstr[s][r]->GetBinContent(ix,iy);
                if (v>0) { sum+=v; sum2+=v*v; cnt++; }
            }
            if (cnt>0) {
                mS[s][r] = sum/cnt;
                sS[s][r] = std::sqrt(sum2/cnt - mS[s][r]*mS[s][r]);
            }
        }
    }

    gStyle->SetOptStat(0);

    // draw Data & MC histograms (original cuts)
    for (size_t s=0; s<species.size(); ++s) {
        auto &lab = species[s].second;
        // Data original
        {
            TCanvas c(Form("cD_%s", lab.c_str()),
                     Form("%s Data (orig)", lab.c_str()), 1800, 600);
            c.Divide(3,1);
            SetSame2DScale(Dorig[s][0], Dorig[s][1], Dorig[s][2]);
            for (int r=0; r<3; ++r) {
                c.cd(r+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                Dorig[s][r]->Draw("COLZ");
            }
            c.SaveAs(Form("output/dc/data_orig_%s.png", lab.c_str()));
        }
        // MC original
        {
            TCanvas c(Form("cM_%s", lab.c_str()),
                     Form("%s MC (orig)", lab.c_str()), 1800, 600);
            c.Divide(3,1);
            SetSame2DScale(Morig[s][0], Morig[s][1], Morig[s][2]);
            for (int r=0; r<3; ++r) {
                c.cd(r+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                Morig[s][r]->Draw("COLZ");
            }
            c.SaveAs(Form("output/dc/mc_orig_%s.png", lab.c_str()));
        }
    }

    // four-color palette for outliers
    Int_t pal[4] = { kBlue, kOrange, kRed, kPink };
    gStyle->SetPalette(4, pal);

    // draw original ratio with survival % in title
    for (size_t s=0; s<species.size(); ++s) {
        auto &lab = species[s].second;
        TCanvas c(Form("RatioOrig_%s", lab.c_str()),
                 Form("%s Ratio Orig", lab.c_str()), 1800, 600);
        c.Divide(3,1);
        SetSame2DScale(Rorig[s][0], Rorig[s][1], Rorig[s][2]);
        for (int r=0; r<3; ++r) {
            c.cd(r+1);
            double pOr = nNo[s][r]?100.0*nOr[s][r]/nNo[s][r]:0.0;
            double pSt = nNo[s][r]?100.0*nSt[s][r]/nNo[s][r]:0.0;
            Rorig[s][r]->SetTitle(
              Form("R%d: Orig %.1f%%, Strict %.1f%%; x; y", r+1, pOr, pSt)
            );
            gPad->SetLeftMargin(.15);
            gPad->SetRightMargin(.15);
            gPad->SetLogz();
            Rorig[s][r]->Draw("COLZ");
        }
        c.SaveAs(Form("output/dc/ratio_orig_%s.png", lab.c_str()));
    }

    // draw strict ratio with survival % in title
    for (size_t s=0; s<species.size(); ++s) {
        auto &lab = species[s].second;
        TCanvas c(Form("RatioStr_%s", lab.c_str()),
                 Form("%s Ratio Strict", lab.c_str()), 1800, 600);
        c.Divide(3,1);
        SetSame2DScale(Rstr[s][0], Rstr[s][1], Rstr[s][2]);
        for (int r=0; r<3; ++r) {
            c.cd(r+1);
            double pOr = nNo[s][r]?100.0*nOr[s][r]/nNo[s][r]:0.0;
            double pSt = nNo[s][r]?100.0*nSt[s][r]/nNo[s][r]:0.0;
            Rstr[s][r]->SetTitle(
              Form("R%d: Orig %.1f%%, Strict %.1f%%; x; y", r+1, pOr, pSt)
            );
            gPad->SetLeftMargin(.15);
            gPad->SetRightMargin(.15);
            gPad->SetLogz();
            Rstr[s][r]->Draw("COLZ");
        }
        c.SaveAs(Form("output/dc/ratio_strict_%s.png", lab.c_str()));
    }

    // draw outliers (orig then strict)
    for (size_t s=0; s<species.size(); ++s) {
        auto &lab = species[s].second;
        // orig
        {
            TCanvas c(Form("OutlOrig_%s", lab.c_str()),
                     Form("%s Outliers Orig", lab.c_str()),
                     1800,600);
            c.Divide(3,1);
            for (int r=0; r<3; ++r) {
                c.cd(r+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                // ensure linear z
                gPad->SetLogz(0);

                TH2D* m = (TH2D*)Rorig[s][r]->Clone();
                m->Reset();

                // fill categorical bins 1–4
                for (int ix=1; ix<=NB2; ++ix) {
                    for (int iy=1; iy<=NB2; ++iy) {
                        double v = Rorig[s][r]->GetBinContent(ix,iy);
                        if (v <= 0) continue;
                        int cat = 4;
                        if      (v >= 0.5  && v <= 2.0)                        cat = 1;
                        else if ((v >= 1./3 && v <  0.5) || (v > 2.0   && v <= 3.0)) cat = 2;
                        else if ((v >= 1./5 && v < 1./3) || (v > 3.0   && v <= 5.0)) cat = 3;
                        m->SetBinContent(ix, iy, cat);
                    }
                }

                // force z-range 1–4 with one bin per category
                m->SetMinimum(1.);
                m->SetMaximum(4.);
                m->GetZaxis()->SetRangeUser(1.,4.);
                m->GetZaxis()->SetNdivisions(4, kFALSE);

                m->SetContour(4);
                for (int i=0; i<4; ++i) m->SetContourLevel(i, i+1);

                m->Draw("COLZ");
            }
            c.SaveAs(Form("output/dc/outliers_orig_%s.png", lab.c_str()));
        }
        // strict
        {
            TCanvas c(Form("OutlStr_%s", lab.c_str()),
                     Form("%s Outliers Strict", lab.c_str()),
                     1800,600);
            c.Divide(3,1);
            for (int r=0; r<3; ++r) {
                c.cd(r+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz(0);

                TH2D* m = (TH2D*)Rstr[s][r]->Clone();
                m->Reset();

                for (int ix=1; ix<=NB2; ++ix) {
                    for (int iy=1; iy<=NB2; ++iy) {
                        double v = Rstr[s][r]->GetBinContent(ix,iy);
                        if (v <= 0) continue;
                        int cat = 4;
                        if      (v >= 0.5  && v <= 2.0)                        cat = 1;
                        else if ((v >= 1./3 && v <  0.5) || (v > 2.0   && v <= 3.0)) cat = 2;
                        else if ((v >= 1./5 && v < 1./3) || (v > 3.0   && v <= 5.0)) cat = 3;
                        m->SetBinContent(ix, iy, cat);
                    }
                }

                m->SetMinimum(1.);
                m->SetMaximum(4.);
                m->GetZaxis()->SetRangeUser(1.,4.);
                m->GetZaxis()->SetNdivisions(4, kFALSE);

                m->SetContour(4);
                for (int i=0; i<4; ++i) m->SetContourLevel(i, i+1);

                m->Draw("COLZ");
            }
            c.SaveAs(Form("output/dc/outliers_strict_%s.png", lab.c_str()));
        }
    }

    return 0;
}