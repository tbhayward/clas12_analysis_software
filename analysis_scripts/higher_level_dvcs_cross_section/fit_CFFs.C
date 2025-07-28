/**
 * fit_CFFs.C
 * ──────────
 * program to fit DVCS Compton Form Factors (CFFs) imaginary then real parts
 *
 * Strategies:
 *   1) fit only Im–parts → BSA data
 *   2) two-step: (a) fit Im–parts → BSA, (b) fit renormReal → xsec
 *
 * Usage:
 *   ./fit_CFFs --strategy <1|2> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>
 *             [--constraint <0|1>] [--input <BSA_file>] [--plot-fits]
 *
 * Compile:
 *   g++ -O2 fit_CFFs.C `root-config --cflags --libs` -lMinuit -o fit_CFFs
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <map>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TStyle.h"

// pull in full BMK_DVCS + CFF code, with globals
#include "DVCS_xsec.C"

// extern flags & renormalizations
extern bool   hasH, hasHt, hasE, hasEt;
extern double renormImag, renormReal;

// ----------------------------------------------------------------------------
// GPD–H defaults (in DVCS_xsec.C)
extern double r_H;
extern double n_H;
extern double alpha0_H;
extern double alpha1_H;
extern double b_H;
extern double beta_H;
extern double M2_H;
extern double P_H;

// GPD–Htilde
extern double r_Ht;
extern double n_Ht;
extern double alpha0_Ht;
extern double alpha1_Ht;
extern double b_Ht;
extern double M2_Ht;
extern double P_Ht;

// GPD–E
extern double r_E;
extern double n_E;
extern double alpha0_E;
extern double alpha1_E;
extern double b_E;
extern double M2_E;
extern double P_E;

// GPD–Etilde
extern double r_Et;
extern double n_Et;
extern double alpha0_Et;
extern double alpha1_Et;
extern double b_Et;
extern double M2_Et;
extern double P_Et;

// control flags
static int   gStrategy    = 0;
static int   gStage       = 1;
static int   gConstraint  = 0;
static bool  gPlotBinFits = false;
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char* gXsFile  = "imports/rga_pass1_xsec_2018.txt";

// raw data + binned observables
struct DataPoint { double phi,Q2,xB,t,Eb,A,sigA; };
static std::vector<DataPoint> bsaData, xsData;
static std::vector<double> bin_xB, bin_Q2, bin_t, bin_Eb;
static std::vector<double> bin_A, bin_dA, bin_chi2;
static std::vector<int>    bin_M;
static int     Nbins          = 0;
static double  reducedAmpChi2 = 0.0;

// forward-declare helper
void PlotBinFit(int ibin, const std::string &ts);

// parse_args: --strategy, -H, -Ht, -E, -Et, [--constraint], [--input], [--plot-fits]
void parse_args(int argc, char** argv){
    static struct option opts[] = {
        {"strategy",   required_argument, nullptr, 's'},
        {"H",          required_argument, nullptr, 'h'},
        {"Ht",         required_argument, nullptr, 't'},
        {"E",          required_argument, nullptr, 'e'},
        {"Et",         required_argument, nullptr, 'x'},
        {"constraint", required_argument, nullptr, 'C'},
        {"input",      required_argument, nullptr, 'i'},
        {"plot-fits",  no_argument,       nullptr, 'p'},
        {nullptr,0,nullptr,0}
    };
    int c;
    while((c=getopt_long(argc,argv,"s:h:t:e:x:C:i:p",opts,nullptr))!=-1){
        switch(c){
          case 's': gStrategy   = std::atoi(optarg); break;
          case 'h': hasH        = std::atoi(optarg); break;
          case 't': hasHt       = std::atoi(optarg); break;
          case 'e': hasE        = std::atoi(optarg); break;
          case 'x': hasEt       = std::atoi(optarg); break;
          case 'C': gConstraint = std::atoi(optarg); break;
          case 'i': gBsaFile    = std::string(optarg); break;
          case 'p': gPlotBinFits= true;                break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy<1|2> -H<0|1> -Ht<0|1>"
                     <<" -E<0|1> -Et<0|1> [--constraint<0|1>]"
                     <<" [--input <BSA_file>] [--plot-fits]\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2){
        std::cerr<<"Invalid strategy\n"; std::exit(1);
    }
}

// Load BSA + XSC
void LoadData(){
    auto read=[&](const char* fn, auto &v){
        std::ifstream in(fn);
        if(!in){std::cerr<<"ERROR: cannot open "<<fn<<"\n"; std::exit(1);}
        std::string line;
        while(std::getline(in,line)){
            if(line.empty()||line[0]=='#') continue;
            std::istringstream iss(line);
            DataPoint d;
            iss>>d.phi>>d.Q2>>d.xB>>d.t>>d.Eb>>d.A>>d.sigA;
            if(gConstraint==1 && (-d.t/d.Q2)>=0.2) continue;
            v.push_back(d);
        }
    };
    read(gBsaFile.c_str(), bsaData);
    read(gXsFile,          xsData);
}

// Bin & sinφ fit → bin_A,bin_dA,bin_chi2
void BinBsaData(){
    bin_xB.clear(); bin_Q2.clear(); bin_t.clear(); bin_Eb.clear();
    bin_A.clear(); bin_dA.clear(); bin_chi2.clear(); bin_M.clear();
    if(bsaData.empty()) return;
    size_t start=0;
    for(size_t i=1;i<=bsaData.size();++i){
        bool newbin=(i==bsaData.size()||bsaData[i].phi<bsaData[i-1].phi);
        if(newbin){
            auto pts=std::vector<DataPoint>(bsaData.begin()+start,
                                            bsaData.begin()+i);
            bin_M.push_back(pts.size());
            double SwA=0,Sw2=0;
            for(auto &d:pts){
                double s=std::sin(d.phi*TMath::Pi()/180.);
                double w=1./(d.sigA*d.sigA);
                SwA+=w*d.A*s; Sw2+=w*s*s;
            }
            double A_bin=SwA/Sw2, dA_bin=1./std::sqrt(Sw2);
            bin_A.push_back(A_bin); bin_dA.push_back(dA_bin);
            double chi2=0;
            for(auto &d:pts){
                double s=std::sin(d.phi*TMath::Pi()/180.);
                double diff=d.A - A_bin*s;
                chi2+=diff*diff/(d.sigA*d.sigA);
            }
            bin_chi2.push_back(chi2);
            double sumw=0,Sx=0,Sq=0,St=0,Se=0;
            for(auto &d:pts){
                double w=1./(d.sigA*d.sigA);
                sumw+=w; Sx+=w*d.xB; Sq+=w*d.Q2; St+=w*d.t; Se+=w*d.Eb;
            }
            bin_xB.push_back(Sx/sumw);
            bin_Q2.push_back(Sq/sumw);
            bin_t .push_back(St/sumw);
            bin_Eb.push_back(Se/sumw);
            start=i;
        }
    }
    Nbins=bin_A.size();
    double totChi2=std::accumulate(bin_chi2.begin(),bin_chi2.end(),0.0);
    int totDof=std::accumulate(bin_M.begin(),bin_M.end(),0)-Nbins;
    reducedAmpChi2 = totDof>0? totChi2/totDof : 0.0;
}

// prepare list of Im-fitting parameters
static std::vector<std::string> parNamesIm;
void build_par_list(){
    parNamesIm.clear();
    if(hasH )  parNamesIm.insert(parNamesIm.end(),
               // {"r_H","alpha0_H","alpha1_H","b_H","M2_H","P_H"});
                {"r_H","alpha0_H","alpha1_H", "b_H", "beta_H"});
    if(hasHt)  parNamesIm.insert(parNamesIm.end(),
               {"r_Ht","alpha0_Ht","alpha1_Ht","b_Ht","M2_Ht","P_Ht"});
    if(hasE )  parNamesIm.insert(parNamesIm.end(),
               {"r_E","alpha0_E","alpha1_E","b_E","M2_E","P_E"});
    if(hasEt)  parNamesIm.insert(parNamesIm.end(),
               {"r_Et","alpha0_Et","alpha1_Et","b_Et","M2_Et","P_Et"});
}

// χ² function for Im-fit and RenormReal-fit
void fcn(int&, double*, double &f, double *par, int){
    int ip=0;
    if(gStage==1){
        if(hasH ){
          r_H       = par[ip++];
          alpha0_H  = par[ip++];
          alpha1_H  = par[ip++];
          b_H       = par[ip++];
          M2_H      = par[ip++];
          P_H       = par[ip++];
        }
        if(hasHt){
          r_Ht      = par[ip++];
          alpha0_Ht = par[ip++];
          alpha1_Ht = par[ip++];
          b_Ht      = par[ip++];
          M2_Ht     = par[ip++];
          P_Ht      = par[ip++];
        }
        if(hasE ){
          r_E       = par[ip++];
          alpha0_E  = par[ip++];
          alpha1_E  = par[ip++];
          b_E       = par[ip++];
          M2_E      = par[ip++];
          P_E       = par[ip++];
        }
        if(hasEt){
          r_Et      = par[ip++];
          alpha0_Et = par[ip++];
          alpha1_Et = par[ip++];
          b_Et      = par[ip++];
          M2_Et     = par[ip++];
          P_Et      = par[ip++];
        }
        double chi2=0;
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0,
                          bin_Eb[k],bin_xB[k],bin_Q2[k],bin_t[k],0.0);
            double modelA = dvcs.s1_I()/dvcs.c0_BH();
            double r      = (bin_A[k] - modelA)/bin_dA[k];
            chi2 += r*r;
        }
        f = chi2;
    } else {
        renormReal = par[ip++];
        double chi2=0;
        for(auto &d: xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double m=dvcs.CrossSection();
            double r=(d.A - renormReal*m)/d.sigA;
            chi2+=r*r;
        }
        f = chi2;
    }
}

int main(int argc, char** argv) {
    parse_args(argc, argv);

    // timestamp
    time_t now = time(nullptr);
    char tb[32];
    strftime(tb, sizeof(tb), "%Y%m%d_%H%M%S", localtime(&now));
    std::cout << "\n=== Strategy=" << gStrategy
              << "  H=" << hasH << " Ht=" << hasHt
              << "  E=" << hasE  << " Et=" << hasEt
              << "  constraint=" << gConstraint
              << "  input=" << gBsaFile
              << "  plot-fits=" << (gPlotBinFits ? "ON" : "OFF")
              << "  ts=" << tb << " ===\n";

    LoadData();
    BinBsaData();
    std::cout << " BSA bins=" << Nbins << " (raw=" << bsaData.size() << ")\n";
    std::cout << " Reduced χ² per amp-fit = " << reducedAmpChi2 << "\n\n";

    if (gPlotBinFits) {
        gStyle->SetOptStat(0);
        for (int ib = 0; ib < Nbins; ++ib) PlotBinFit(ib, tb);
        std::cout << " Wrote plots to output/plots/binned_fits/\n\n";
    }

    // ─── Stage 1: Im-fit ─────────────────────────────────────────────────────────
    gStage = 1;
    build_par_list();
    int nim = parNamesIm.size();
    std::vector<double> imVal(nim), imErr(nim);
    double chi2_im, edm, errdef;
    int nv, nx, ic, ndf_im;
    {
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);

        // ── define all Im–fit parameters ─────────
        for (int i = 0; i < nim; ++i) {
            const auto &nm = parNamesIm[i];
            double step = 0.01;
            double init = 0.0, lo = -1e3, hi = 1e3;

            if      (nm == "r_H")       init = r_H;
            else if (nm == "alpha0_H")  { init = alpha0_H; }
            else if (nm == "alpha1_H")  init = alpha1_H;
            else if (nm == "b_H")       { init = b_H; lo = 0.0; }
            else if (nm == "beta_H")    { init = beta_H;  }

            else if (nm == "r_Ht")      init = r_Ht;
            else if (nm == "alpha0_Ht") { init = alpha0_Ht; lo = 0.2; hi = 1.0; }
            else if (nm == "alpha1_Ht") init = alpha1_Ht;
            else if (nm == "b_Ht")      init = b_Ht;
            else if (nm == "M2_Ht")     { init = M2_Ht;    lo = 0.0; hi = 2.0; }
            else if (nm == "P_Ht")      { init = P_Ht;     lo = 0.0; hi = 5.0; }

            else if (nm == "r_E")       init = r_E;
            else if (nm == "alpha0_E")  { init = alpha0_E; lo = 0.2; hi = 1.0; }
            else if (nm == "alpha1_E")  init = alpha1_E;
            else if (nm == "b_E")       init = b_E;
            else if (nm == "M2_E")      { init = M2_E;     lo = 0.0; hi = 2.0; }
            else if (nm == "P_E")       { init = P_E;      lo = 0.0; hi = 5.0; }

            else if (nm == "r_Et")      init = r_Et;
            else if (nm == "alpha0_Et") { init = alpha0_Et;lo = 0.2; hi = 1.0; }
            else if (nm == "alpha1_Et") init = alpha1_Et;
            else if (nm == "b_Et")      init = b_Et;
            else if (nm == "M2_Et")     { init = M2_Et;    lo = 0.0; hi = 2.0; }
            else if (nm == "P_Et")      { init = P_Et;     lo = 0.0; hi = 5.0; }

            minu.DefineParameter(i, nm.c_str(), init, step, lo, hi);
        }

        std::cout << " Stage1: fitting Im-CFF parameters…\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_im, edm, errdef, nv, nx, ic);
        for (int i = 0; i < nim; ++i) {
            minu.GetParameter(i, imVal[i], imErr[i]);
        }
        ndf_im = Nbins - nim;
    }

    std::map<std::string,double> valMap, errMap;
    for (int i = 0; i < nim; ++i) {
        valMap[parNamesIm[i]] = imVal[i];
        errMap[parNamesIm[i]] = imErr[i];
    }

    // ─── Stage 2: RenormReal-fit ────────────────────────────────────────────────
    if (gStrategy == 2) {
        gStage = 2;
        double chi2_re, edm2, errdef2;
        int nv2, nx2, ic2, ndf_re;
        double reVal, reErr;
        {
            TMinuit m2(1);
            m2.SetPrintLevel(1);
            m2.SetFCN(fcn);
            m2.DefineParameter(0, "renormReal", renormReal, 0.01, -1e3, 1e3);
            std::cout << " Stage2: fitting renormReal...\n";
            m2.Migrad();
            m2.Command("HESSE");
            m2.mnstat(chi2_re, edm2, errdef2, nv2, nx2, ic2);
            m2.GetParameter(0, reVal, reErr);
            ndf_re = xsData.size() - 1;
        }
        valMap["renormReal"] = reVal;
        errMap["renormReal"] = reErr;
        chi2_im = chi2_re;
        ndf_im  = ndf_re;
    }

    // ─── Write out results ───────────────────────────────────────────────────────
    std::vector<std::string> outNames = parNamesIm;
    if (gStrategy == 2) outNames.push_back("renormReal");
    system("mkdir -p output/fit_results");
    std::ofstream fout(std::string("output/fit_results/fit_results_") + tb + ".txt");
    fout << "# fit_CFFs results\n"
         << "timestamp   " << tb << "\n"
         << "strategy    " << gStrategy << "\n"
         << "constraint  " << gConstraint << "\n"
         << "input       " << gBsaFile << "\n"
         << "H " << hasH << "  Ht " << hasHt
         << "  E " << hasE << "  Et " << hasEt << "\n"
         << "# parameters:";
    for (auto &n : outNames) fout << " " << n;
    fout << "\n# values:\n";
    for (auto &n : outNames) fout << valMap[n] << " ";
    fout << "\n# errors:\n";
    for (auto &n : outNames) fout << errMap[n] << " ";
    fout << "\n# chi2 ndf chi2/ndf\n"
         << chi2_im << " " << ndf_im << " " << (chi2_im / ndf_im) << "\n";
    fout.close();

    std::cout << "\n--- Fit Results ---\n";
    for (auto &n : outNames) {
        std::cout << " " << n << " = " << valMap[n]
                  << " ± " << errMap[n] << "\n";
    }
    std::cout << " χ²/ndf = " << chi2_im << "/" << ndf_im
              << " = " << (chi2_im / ndf_im) << "\n";
    std::cout << " Reduced χ² per amp-fit = " << reducedAmpChi2 << "\n";

    return 0;
}

// ──────────────────────────────────────────────────────────────────────────────
void PlotBinFit(int ibin, const std::string &ts) {
    if (!gPlotBinFits) return;

    // regroup raw points into φ-bins
    std::vector<std::vector<DataPoint>> bins;
    size_t start = 0;
    for (size_t i = 1; i <= bsaData.size(); ++i) {
        bool newbin = (i == bsaData.size() || bsaData[i].phi < bsaData[i-1].phi);
        if (newbin) {
            bins.emplace_back(bsaData.begin() + start, bsaData.begin() + i);
            start = i;
        }
    }
    if (ibin < 0 || ibin >= int(bins.size())) return;
    auto &dp = bins[ibin];
    int n = dp.size();
    if (n < 2) return;

    // ensure directory
    system("mkdir -p output/plots/binned_fits");

    // disable default stat box & grid
    gStyle->SetOptStat(0);

    // graph
    TGraphErrors *gr = new TGraphErrors(n);
    for (int i = 0; i < n; ++i) {
        gr->SetPoint(i, dp[i].phi, dp[i].A);
        gr->SetPointError(i, 0.0, dp[i].sigA);
    }
    gr->SetMarkerStyle(20);

    // fit with offset C + A sinφ/(1+B cosφ)
    TF1 *f1 = new TF1(Form("f_bin%d", ibin),
        "[0] + [1]*sin(x*TMath::Pi()/180)/(1+[2]*cos(x*TMath::Pi()/180))",
        0, 360);
    f1->SetParameter(0, 0.0);
    f1->SetParameter(1, bin_A[ibin]);
    f1->SetParameter(2, 0.0);
    f1->SetParLimits(1, -1.0, 1.0);
    f1->SetParLimits(2, -1.0, 1.0);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    gr->Fit(f1, "RQN");

    // get chi2 and ndf
    double chi2 = f1->GetChisquare();
    double ndf  = f1->GetNDF();

    // canvas + frame
    TCanvas *c = new TCanvas(Form("c_bin%d", ibin), "", 600, 500);
    TH1F *frame = new TH1F(Form("frame%d", ibin), "", 360, 0, 360);
    frame->SetMinimum(-0.6);
    frame->SetMaximum(0.6);
    frame->GetXaxis()->SetTitle("#phi (deg)");
    frame->GetYaxis()->SetTitle("A_{LU}");
    frame->Draw("AXIS");
    frame->GetXaxis()->SetLimits(0, 360);
    frame->GetXaxis()->SetRangeUser(0, 360);
    gPad->Modified(); gPad->Update();

    gr->Draw("P same");
    f1->Draw("L same");

    TLegend *leg = new TLegend(0.60, 0.75, 0.90, 0.90);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(0);
    leg->AddEntry(gr, "data", "p");
    leg->AddEntry(f1,
        Form("C = %.3f +/- %.3f", f1->GetParameter(0), f1->GetParError(0)),
        "l");
    leg->AddEntry(f1,
        Form("A = %.3f +/- %.3f", f1->GetParameter(1), f1->GetParError(1)),
        "l");
    leg->AddEntry(f1,
        Form("B = %.3f +/- %.3f", f1->GetParameter(2), f1->GetParError(2)),
        "l");
    leg->AddEntry((TObject*)0,
        Form("#chi^{2}/ndf = %.2f", chi2/ndf),
        "");
    leg->Draw();

    c->SaveAs(Form("output/plots/binned_fits/BinFit_%s_bin%d.pdf", ts.c_str(), ibin));

    delete leg;
    delete frame;
    delete c;
    delete gr;
    delete f1;
}