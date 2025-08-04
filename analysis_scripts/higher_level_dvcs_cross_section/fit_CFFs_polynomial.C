/**
 * fit_CFFs_polynomial.cpp
 * ───────────────────────
 * program to fit DVCS Compton Form Factors (CFFs) imaginary then real parts
 * using a polynomial+exponential ansatz for the Im parts and dispersion for Re.
 *
 * Strategies:
 *   1) fit only Im–parts → BSA data (extracted sine amplitudes)
 *   2) two-step: (a) fit Im–parts → BSA, (b) fit renormReal → xsec
 *
 * Usage:
 *   ./fit_CFFs --strategy <1|2> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>
 *             [--constraint <0|1>] [--input <BSA_file>] [--plot-fits]
 *
 *   constraint=0 : no cuts
 *   constraint=1 : per-point (-t/Q2>=0.2 or -t>=1) AND per-bin
 *                  (χ²_red>=2 or |B|>=0.9) cuts applied to the sine fits
 *
 * Compile:
 *   g++ -O2 fit_CFFs_polynomial.cpp `root-config --cflags --libs` -lMinuit -o fit_CFFs_polynomial
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
#include <functional>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TObject.h"

// pull in the polynomial DVCS machinery (defines BMK_DVCS, GetIm*/GetRe*, and all globals)
#include "DVCS_xsec_polynomial.C"

// -----------------------------------------------------------------------------
// externs from DVCS_xsec_polynomial.cpp so we can drive the fit
extern bool hasH, hasHt, hasE, hasEt;
extern double renormImag, renormReal;

// Imaginary part polynomial parameters
extern double A_ImH, B_ImH, c1_ImH, c2_ImH, c3_ImH, d1_ImH, d2_ImH, d3_ImH;
extern double A_ImHt, B_ImHt, c1_ImHt, c2_ImHt, c3_ImHt, d1_ImHt, d2_ImHt, d3_ImHt;
extern double A_ImE, B_ImE, c1_ImE, c2_ImE, c3_ImE, d1_ImE, d2_ImE, d3_ImE;
extern double A_ImEt, B_ImEt, c1_ImEt, c2_ImEt, c3_ImEt, d1_ImEt, d2_ImEt, d3_ImEt;

// Real-part subtraction constants (D-term-like)
extern double C0_H, MD2_H, lambda_H;
extern double C0_Ht, MD2_Ht, lambda_Ht;
extern double C0_E, MD2_E, lambda_E;
extern double C0_Et, MD2_Et, lambda_Et;

// Forward declarations for CFF access (should already be visible via include)
double GetImH(double xi, double t);
double GetImHt(double xi, double t);
double GetImE(double xi, double t);
double GetImEt(double xi, double t);
double GetReH(double xi, double t);
double GetReHt(double xi, double t);
double GetReE(double xi, double t);
double GetReEt(double xi, double t);

// -----------------------------------------------------------------------------
// control flags & defaults
static int   gStrategy    = 0;
static int   gStage       = 1;
static int   gConstraint  = 0;  // 0 or 1
static bool  gPlotBinFits = false;
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char* gXsFile  = "imports/rga_pass1_xsec_2018.txt";

// raw data + binned observables
struct DataPoint { double phi,Q2,xB,t,Eb,A,sigA; };
static std::vector<DataPoint> bsaData, xsData;
static std::vector<double>    bin_xB, bin_Q2, bin_t, bin_Eb;
static std::vector<double>    bin_A, bin_dA;
static std::vector<int>       bin_M;
static std::vector<double>    bin_redChi2; // per-bin sine reduced chi2
static int     Nbins          = 0;
static double  reducedAmpChi2 = 0.0;

// helper forward
void PlotBinFit(int ibin, const std::string &ts);

// parse command-line arguments
void parse_args(int argc, char** argv){
    static struct option opts[] = {
      {"strategy",   required_argument,nullptr,'s'},
      {"H",          required_argument,nullptr,'h'},
      {"Ht",         required_argument,nullptr,'t'},
      {"E",          required_argument,nullptr,'e'},
      {"Et",         required_argument,nullptr,'x'},
      {"constraint", required_argument,nullptr,'C'},
      {"input",      required_argument,nullptr,'i'},
      {"plot-fits",  no_argument,      nullptr,'p'},
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
          case 'i': gBsaFile    = optarg;            break;
          case 'p': gPlotBinFits= true;              break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy<1|2> -H<0|1> -Ht<0|1>"
                     <<" -E<0|1> -Et<0|1> [--constraint<0|1>]"
                     <<" [--input <BSA_file>] [--plot-fits]\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2 || gConstraint<0||gConstraint>1){
        std::cerr<<"Invalid strategy or constraint\n";
        std::exit(1);
    }
}

// Load BSA + XSC, apply per-point cuts when constraint==1
void LoadData(){
    auto read=[&](const char* fn, auto &v){
        std::ifstream in(fn);
        if(!in){ std::cerr<<"ERROR: cannot open "<<fn<<"\n"; std::exit(1); }
        std::string line;
        while(std::getline(in,line)){
            if(line.empty()||line[0]=='#') continue;
            std::istringstream iss(line);
            DataPoint d; iss>>d.phi>>d.Q2>>d.xB>>d.t>>d.Eb>>d.A>>d.sigA;
            if(gConstraint==1){
              if((-d.t/d.Q2) >= 0.2) continue;
              if((-d.t)     >= 1.0) continue;
            }
            v.push_back(d);
        }
    };
    read(gBsaFile.c_str(), bsaData);
    read(gXsFile,          xsData);
}

// Bin & sinφ fit → bin_A,bin_dA; then drop bins failing the per-bin cuts
void BinBsaData(){
    bin_xB.clear(); bin_Q2.clear(); bin_t.clear(); bin_Eb.clear();
    bin_A.clear(); bin_dA.clear(); bin_M.clear(); bin_redChi2.clear();
    if(bsaData.empty()) return;

    size_t start=0;
    for(size_t i=1;i<=bsaData.size();++i){
        bool newbin = (i==bsaData.size() || bsaData[i].phi < bsaData[i-1].phi);
        if(!newbin) continue;
        std::vector<DataPoint> pts(bsaData.begin()+start, bsaData.begin()+i);
        start = i;
        int M = pts.size();
        if(M<3) continue; // need at least 3 points

        // weighted estimate for A and its uncertainty
        double SwA=0, Sw2=0;
        for(auto &d: pts){
            double s = std::sin(d.phi * TMath::Pi()/180.);
            double w = 1.0/(d.sigA * d.sigA);
            SwA += w * d.A * s;
            Sw2 += w * s * s;
        }
        double A_bin = SwA / Sw2;
        double dA_bin = 1.0 / std::sqrt(Sw2);

        // perform three-parameter fit: C + A sinφ/(1 + B cosφ)
        TF1 ftmp("ftmp",
                 "[0] + [1]*sin(x*TMath::Pi()/180.)/(1+[2]*cos(x*TMath::Pi()/180.))",
                 0, 360);
        ftmp.SetParameter(0, 0.0);     // C
        ftmp.SetParameter(1, A_bin);   // A
        ftmp.SetParameter(2, 0.0);     // B
        TGraphErrors gr(M);
        for(int j=0;j<M;++j){
            gr.SetPoint(j, pts[j].phi, pts[j].A);
            gr.SetPointError(j, 0, pts[j].sigA);
        }
        gr.SetMarkerStyle(20);
        gr.Fit(&ftmp, "Q0R"); // quiet, no range restriction beyond [0,360]

        double Bfit    = ftmp.GetParameter(2);
        double chi2    = ftmp.GetChisquare();
        double redchi2 = (M>1) ? chi2 / (M - 1) : 1e9;

        // apply per-bin cuts after the fit
        if(gConstraint==1){
          if(redchi2 >= 2.0 || std::fabs(Bfit) >= 0.9) continue;
        }

        // store surviving bin (use weighted A_bin, not fitted A to preserve original definition)
        bin_M.push_back(M);
        bin_A.push_back(A_bin);
        bin_dA.push_back(dA_bin);
        bin_redChi2.push_back(redchi2);

        double sumw=0, Sx=0, Sq=0, St=0, Se=0;
        for(auto &d: pts){
            double w = 1.0/(d.sigA * d.sigA);
            sumw += w; Sx += w * d.xB; Sq += w * d.Q2; St += w * d.t; Se += w * d.Eb;
        }
        bin_xB.push_back(Sx/sumw);
        bin_Q2.push_back(Sq/sumw);
        bin_t .push_back(St/sumw);
        bin_Eb.push_back(Se/sumw);
    }

    Nbins = bin_A.size();

    // compute overall reduced χ² per amp-fit: sum((A/dA)^2) over total degrees of freedom
    double totChi2=0; int totDof=0;
    for(int k=0;k<Nbins;++k){
        totChi2 += std::pow(bin_A[k]/bin_dA[k],2);
        totDof += (bin_M[k]-1);
    }
    reducedAmpChi2 = totDof>0 ? totChi2 / totDof : 0.0;
}

// build list of parameters to float in Im-fit (polynomial coefficients) depending on active CFFs
static std::vector<std::string> parNamesIm;
void build_par_list(){
    parNamesIm.clear();
    if(hasH){
        parNamesIm.insert(parNamesIm.end(),
            {"A_ImH","B_ImH",
             "c1_ImH","c2_ImH","c3_ImH",
             "d1_ImH","d2_ImH","d3_ImH"});
    }
    if(hasHt){
        parNamesIm.insert(parNamesIm.end(),
            {"A_ImHt","B_ImHt",
             "c1_ImHt","c2_ImHt","c3_ImHt",
             "d1_ImHt","d2_ImHt","d3_ImHt"});
    }
    if(hasE){
        parNamesIm.insert(parNamesIm.end(),
            {"A_ImE","B_ImE",
             "c1_ImE","c2_ImE","c3_ImE",
             "d1_ImE","d2_ImE","d3_ImE"});
    }
    if(hasEt){
        parNamesIm.insert(parNamesIm.end(),
            {"A_ImEt","B_ImEt",
             "c1_ImEt","c2_ImEt","c3_ImEt",
             "d1_ImEt","d2_ImEt","d3_ImEt"});
    }
}

// χ² function used by TMinuit
void fcn(int&, double*, double &f, double *par, int){
    int ip=0;
    if(gStage==1){
        // assign polynomial parameters in the same order as build_par_list()
        if(hasH){
            A_ImH   = par[ip++];
            B_ImH   = par[ip++];
            c1_ImH  = par[ip++];
            c2_ImH  = par[ip++];
            c3_ImH  = par[ip++];
            d1_ImH  = par[ip++];
            d2_ImH  = par[ip++];
            d3_ImH  = par[ip++];
        }
        if(hasHt){
            A_ImHt  = par[ip++];
            B_ImHt  = par[ip++];
            c1_ImHt = par[ip++];
            c2_ImHt = par[ip++];
            c3_ImHt = par[ip++];
            d1_ImHt = par[ip++];
            d2_ImHt = par[ip++];
            d3_ImHt = par[ip++];
        }
        if(hasE){
            A_ImE   = par[ip++];
            B_ImE   = par[ip++];
            c1_ImE  = par[ip++];
            c2_ImE  = par[ip++];
            c3_ImE  = par[ip++];
            d1_ImE  = par[ip++];
            d2_ImE  = par[ip++];
            d3_ImE  = par[ip++];
        }
        if(hasEt){
            A_ImEt  = par[ip++];
            B_ImEt  = par[ip++];
            c1_ImEt = par[ip++];
            c2_ImEt = par[ip++];
            c3_ImEt = par[ip++];
            d1_ImEt = par[ip++];
            d2_ImEt = par[ip++];
            d3_ImEt = par[ip++];
        }

        double chi2=0;
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0,
                          bin_Eb[k], bin_xB[k], bin_Q2[k], bin_t[k], 0.0);
            double modelA = dvcs.s1_I()/dvcs.c0_BH();
            double resid  = (bin_A[k] - modelA)/bin_dA[k];
            chi2 += resid*resid;
        }
        f = chi2;
    } else {
        // Stage 2: fit renormReal to cross-section data
        renormReal = par[ip++];
        double chi2=0;
        for(auto &d: xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double m     = dvcs.CrossSection();
            double resid = (d.A - renormReal*m)/d.sigA;
            chi2 += resid*resid;
        }
        f = chi2;
    }
}

int main(int argc, char** argv) {
    parse_args(argc, argv);

    // timestamp
    time_t now = time(nullptr);
    char tb[32];
    strftime(tb,sizeof(tb),"%Y%m%d_%H%M%S",localtime(&now));
    std::cout<<"\n=== Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<"  E="<<hasE<<" Et="<<hasEt
             <<"  constraint="<<gConstraint
             <<"  input="<<gBsaFile
             <<"  plot-fits="<<(gPlotBinFits?"ON":"OFF")
             <<"  ts="<<tb<<" ===\n\n";

    LoadData();
    BinBsaData();

    // print binned points
    std::cout<<"Data points entering Im-fit:\n";
    for(int k=0;k<Nbins;++k){
      std::cout<<" Point "<<(k+1)
               <<", Eb="<<bin_Eb[k]
               <<", xB="<<bin_xB[k]
               <<", Q2="<<bin_Q2[k]
               <<", -t="<<-bin_t[k]
               <<", Amp="<<bin_A[k]
               <<"+/-"<<bin_dA[k]
               <<", χ²_red="<<bin_redChi2[k]
               <<"\n";
    }

    double avgBinChi2 = bin_redChi2.empty()
                      ? 0.0
                      : std::accumulate(bin_redChi2.begin(),bin_redChi2.end(),0.0)
                        / bin_redChi2.size();
    std::cout<<"\n BSA bins="<<Nbins
             <<" (raw="<<bsaData.size()<<")\n"
             <<" Average χ² per bin-fit = "<<avgBinChi2<<"\n"
             <<" Reduced χ² per amp-fit = "<<reducedAmpChi2<<"\n\n";

    if(gPlotBinFits){
      gStyle->SetOptStat(0);
      for(int ib=0; ib<Nbins; ++ib) PlotBinFit(ib,tb);
      std::cout<<" Wrote plots to output/plots/binned_fits/\n\n";
    }

    // ─── Stage 1: Im-fit ─────────────────────────────────────────────────────────
    gStage = 1;
    build_par_list();
    int nim = parNamesIm.size();
    std::vector<double> imVal(nim), imErr(nim);
    double chi2_im=0, edm=0, errdef=0;
    int nv=0,nx=0,ic=0, ndf_im=0;
    {
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);

        for (int i = 0; i < nim; ++i) {
            const auto &nm = parNamesIm[i];
            double init = 1.0;
            double step = 0.01;
            double lo = -10.0, hi = 10.0;

            // sensible starting values from current globals
            if (nm == "A_ImH")       init = A_ImH, lo = 0.0, hi = 10.0;
            else if (nm == "B_ImH")  init = B_ImH, lo = 0.0, hi = 10.0;
            else if (nm == "c1_ImH") init = c1_ImH;
            else if (nm == "c2_ImH") init = c2_ImH;
            else if (nm == "c3_ImH") init = c3_ImH;
            else if (nm == "d1_ImH") init = d1_ImH;
            else if (nm == "d2_ImH") init = d2_ImH;
            else if (nm == "d3_ImH") init = d3_ImH;

            else if (nm == "A_ImHt")       init = A_ImHt, lo = 0.0, hi = 10.0;
            else if (nm == "B_ImHt")       init = B_ImHt, lo = 0.0, hi = 10.0;
            else if (nm == "c1_ImHt")      init = c1_ImHt;
            else if (nm == "c2_ImHt")      init = c2_ImHt;
            else if (nm == "c3_ImHt")      init = c3_ImHt;
            else if (nm == "d1_ImHt")      init = d1_ImHt;
            else if (nm == "d2_ImHt")      init = d2_ImHt;
            else if (nm == "d3_ImHt")      init = d3_ImHt;

            else if (nm == "A_ImE")       init = A_ImE, lo = 0.0, hi = 10.0;
            else if (nm == "B_ImE")       init = B_ImE, lo = 0.0, hi = 10.0;
            else if (nm == "c1_ImE")      init = c1_ImE;
            else if (nm == "c2_ImE")      init = c2_ImE;
            else if (nm == "c3_ImE")      init = c3_ImE;
            else if (nm == "d1_ImE")      init = d1_ImE;
            else if (nm == "d2_ImE")      init = d2_ImE;
            else if (nm == "d3_ImE")      init = d3_ImE;

            else if (nm == "A_ImEt")       init = A_ImEt, lo = 0.0, hi = 10.0;
            else if (nm == "B_ImEt")       init = B_ImEt, lo = 0.0, hi = 10.0;
            else if (nm == "c1_ImEt")      init = c1_ImEt;
            else if (nm == "c2_ImEt")      init = c2_ImEt;
            else if (nm == "c3_ImEt")      init = c3_ImEt;
            else if (nm == "d1_ImEt")      init = d1_ImEt;
            else if (nm == "d2_ImEt")      init = d2_ImEt;
            else if (nm == "d3_ImEt")      init = d3_ImEt;

            minu.DefineParameter(i, nm.c_str(), init, step, lo, hi);
        }

        std::cout<<" Stage1: fitting Im-CFF parameters…\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_im,edm,errdef,nv,nx,ic);
        for (int i = 0; i < nim; ++i) {
            minu.GetParameter(i, imVal[i], imErr[i]);
        }
        ndf_im = Nbins - nim;
    }

    // collect Im-fit results
    std::map<std::string,double> valMap, errMap;
    for (int i = 0; i < nim; ++i) {
        valMap[parNamesIm[i]] = imVal[i];
        errMap[parNamesIm[i]] = imErr[i];
    }

    // ─── Stage 2: renormReal-fit ────────────────────────────────────────────────
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
            std::cout<<" Stage2: fitting renormReal…\n";
            m2.Migrad();
            m2.Command("HESSE");
            m2.mnstat(chi2_re,edm2,errdef2,nv2,nx2,ic2);
            m2.GetParameter(0, reVal, reErr);
            ndf_re = xsData.size() - 1;
        }
        valMap["renormReal"] = reVal;
        errMap["renormReal"] = reErr;
        chi2_im = chi2_re;
        ndf_im = ndf_re;
    }

    // ─── Output results ─────────────────────────────────────────────────────────
    std::vector<std::string> outNames = parNamesIm;
    if (gStrategy == 2) outNames.push_back("renormReal");
    system("mkdir -p output/fit_results");
    std::ofstream fout(std::string("output/fit_results/fit_results_") + tb + ".txt");
    fout<<"# fit_CFFs_polynomial results\n"
        <<"timestamp   "<<tb<<"\n"
        <<"strategy    "<<gStrategy<<"\n"
        <<"constraint  "<<gConstraint<<"\n"
        <<"input       "<<gBsaFile<<"\n"
        <<"H "<<hasH<<" Ht "<<hasHt<<" E "<<hasE<<" Et "<<hasEt<<"\n"
        <<"# parameters:";
    for(auto &n: outNames) fout<<" "<<n;
    fout<<"\n# values:\n";
    for(auto &n: outNames) fout<<valMap[n]<<" ";
    fout<<"\n# errors:\n";
    for(auto &n: outNames) fout<<errMap[n]<<" ";
    fout<<"\n# chi2 ndf chi2/ndf\n"
        <<chi2_im<<" "<<ndf_im<<" "<<(ndf_im>0?chi2_im/ndf_im:0)<<"\n";
    fout.close();

    // print to stdout
    std::cout<<"\n--- Fit Results ---\n";
    for(auto &n: outNames){
      std::cout<<" "<<n<<" = "<<valMap[n]
               <<" ± "<<errMap[n]<<"\n";
    }
    std::cout<<" χ²/ndf = "<<chi2_im<<"/"<<ndf_im
             <<" = "<<(ndf_im>0?chi2_im/ndf_im:0)<<"\n";
    std::cout<<" Average χ² per bin-fit = "<<avgBinChi2<<"\n";
    std::cout<<" Reduced χ² per amp-fit = "<<reducedAmpChi2<<"\n";

    return 0;
}

// ──────────────────────────────────────────────────────────────────────────────
void PlotBinFit(int ibin, const std::string &ts) {
    if (!gPlotBinFits) return;

    // regroup raw points into φ-bins
    std::vector<std::vector<DataPoint>> bins;
    size_t start=0;
    for(size_t i=1;i<=bsaData.size();++i){
      bool newbin=(i==bsaData.size()||bsaData[i].phi<bsaData[i-1].phi);
      if(newbin){
        bins.emplace_back(bsaData.begin()+start,bsaData.begin()+i);
        start=i;
      }
    }
    if(ibin<0||ibin>=int(bins.size())) return;
    auto &dp=bins[ibin];
    int n=dp.size(); if(n<3) return;   // skip if insufficient points

    system("mkdir -p output/plots/binned_fits");
    gStyle->SetOptStat(0);

    // graph with errors and visible marker
    TGraphErrors *gr=new TGraphErrors(n);
    for(int i=0;i<n;++i){
      gr->SetPoint(i, dp[i].phi, dp[i].A);
      gr->SetPointError(i, 0, dp[i].sigA);
    }
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);

    // perform full three-parameter fit for display: C + A sinφ/(1+B cosφ)
    TF1 *f1=new TF1(Form("f_bin%d",ibin),
      "[0] + [1]*sin(x*TMath::Pi()/180.)/(1+[2]*cos(x*TMath::Pi()/180.))",
      0,360);
    f1->SetParameter(0, 0.0);
    f1->SetParameter(1, bin_A[ibin]);
    f1->SetParameter(2, 0.0);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    gr->Fit(f1,"RQN");

    double chi2=f1->GetChisquare(), ndf=f1->GetNDF();

    TCanvas *c=new TCanvas(Form("c_bin%d",ibin),"",600,500);
    TH1F *frame=new TH1F(Form("frame%d",ibin),"",360,0,360);
    frame->SetMinimum(-0.6); frame->SetMaximum(0.6);
    frame->GetXaxis()->SetTitle("#phi (deg)");
    frame->GetYaxis()->SetTitle("A_{LU}");
    frame->Draw("AXIS");
    frame->GetXaxis()->SetLimits(0,360);
    frame->GetXaxis()->SetRangeUser(0,360);
    gPad->Modified(); gPad->Update();

    gr->Draw("P same");
    f1->Draw("L same");

    TLegend *leg=new TLegend(0.60,0.75,0.90,0.90);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(0);
    leg->AddEntry(gr,"data","p");
    leg->AddEntry(f1,
        Form("C = %.3f ± %.3f", f1->GetParameter(0), f1->GetParError(0)),
        "l");
    leg->AddEntry(f1,
        Form("A = %.3f ± %.3f", f1->GetParameter(1), f1->GetParError(1)),
        "l");
    leg->AddEntry(f1,
        Form("B = %.3f ± %.3f", f1->GetParameter(2), f1->GetParError(2)),
        "l");
    leg->AddEntry((TObject*)0,
        Form("#chi^{2}/ndf = %.2f", (ndf>0?chi2/ndf:0.)),
        "");
    leg->Draw();

    c->SaveAs(Form("output/plots/binned_fits/BinFit_%s_bin%d.pdf",ts.c_str(),ibin));

    delete leg; delete frame; delete c; delete gr; delete f1;
}