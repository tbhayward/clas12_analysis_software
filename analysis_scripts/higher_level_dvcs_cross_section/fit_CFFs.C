/**
 * fit_CFFs.C
 * ──────────
 * program to fit DVCS Compton Form Factors (CFFs) imaginary then real parts
 *
 * Strategies:
 *   1) fit only Im–parts -> BSA data
 *   2) two-step: (a) fit Im–parts -> BSA, (b) fit only renormReal -> xsec
 *
 * Usage:
 *   ./fit_CFFs --strategy <1|2> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>
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
#include <cstdlib>
#include <ctime>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"

// pull in full BMK_DVCS + CFF code, with globals
#include "DVCS_xsec.C"

// extern flags & renormalizations
extern bool   hasH, hasHt, hasE, hasEt;
extern double renormImag, renormReal;

// -------------------------------------------------------------------------------------------------
//   Compton Form Factor (CFF) models—imaginary parts
// -------------------------------------------------------------------------------------------------
extern double alpha0_H,  alpha1_H,  n_H,   b_H,   M2_H,  P_H;
extern double alpha0_Ht, alpha1_Ht, n_Ht,  b_Ht,  M2_Ht, P_Ht;
extern double alpha0_E,  alpha1_E,  n_E,   b_E,   M2_E,  P_E;
extern double alpha0_Et, alpha1_Et, n_Et,  b_Et,  M2_Et, P_Et;

// ──────────────────────────────────────────────────────────────────────────────
// global controls (must be before LoadData)
static int  gStrategy    = 0;       // 1 or 2
static int  gStage       = 1;       // 1=Im fit, 2=Re fit
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char* gXsFile    = "imports/rga_pass1_xsec_2018.txt";

// ──────────────────────────────────────────────────────────────────────────────
// Data structures & loading
struct DataPoint { double phi,Q2,xB,t,Eb,A,sigA; };
static std::vector<DataPoint> bsaData, xsData;

// Binned observables:
static std::vector<double> bin_phi, bin_xB, bin_Q2, bin_t, bin_Eb;
static std::vector<double> bin_A, bin_dA;
static int Nbins = 0;

void LoadData(){
    auto read = [&](const char* fn, std::vector<DataPoint>& v){
        std::ifstream in(fn);
        if(!in){ std::cerr<<"ERROR: cannot open "<<fn<<"\n"; std::exit(1); }
        std::string line;
        while(std::getline(in,line)){
            if(line.empty()||line[0]=='#') continue;
            std::istringstream iss(line);
            DataPoint d;
            iss>>d.phi>>d.Q2>>d.xB>>d.t>>d.Eb>>d.A>>d.sigA;
            v.push_back(d);
        }
    };
    read(gBsaFile.c_str(), bsaData);
    read(gXsFile,          xsData);
}

// ──────────────────────────────────────────────────────────────────────────────
// Bin the BSA data by contiguous (xB, Q2, t) blocks and extract
// weighted sin-phi amplitudes A_bin ± dA_bin
void BinBsaData(){
    bin_A.clear(); bin_dA.clear();
    bin_xB.clear(); bin_Q2.clear(); bin_t.clear(); bin_Eb.clear();
    if(bsaData.empty()) return;

    size_t start=0;
    auto flush_bin = [&](size_t end){
        double sum_w_sinA=0, sum_w_sin2=0;
        double sum_xB=0, sum_Q2=0, sum_t=0, sum_Eb=0;
        int M = end - start;
        for(size_t i=start;i<end;++i){
            auto &d = bsaData[i];
            double phi_rad = d.phi * TMath::Pi()/180.0;
            double s = sin(phi_rad);
            double w = 1.0/(d.sigA*d.sigA);
            sum_w_sinA += w * d.A * s;
            sum_w_sin2  += w * s*s;
            sum_xB += d.xB; sum_Q2 += d.Q2; sum_t += d.t; sum_Eb += d.Eb;
        }
        double Abin = sum_w_sinA / sum_w_sin2;
        double dA   = 1.0/std::sqrt(sum_w_sin2);
        bin_A.push_back(Abin);
        bin_dA.push_back(dA);
        bin_xB.push_back(sum_xB/M);
        bin_Q2.push_back(sum_Q2/M);
        bin_t.push_back(sum_t/M);
        bin_Eb.push_back(sum_Eb/M);
    };

    for(size_t i=1;i<=bsaData.size();++i){
        bool newbin = (i==bsaData.size()
                   || bsaData[i].xB != bsaData[i-1].xB
                   || bsaData[i].Q2 != bsaData[i-1].Q2
                   || bsaData[i].t  != bsaData[i-1].t);
        if(newbin){
            flush_bin(i);
            start = i;
        }
    }
    Nbins = bin_A.size();
}

// ──────────────────────────────────────────────────────────────────────────────
// parse_args(): read --strategy, -H, -Ht, -E, -Et
void parse_args(int argc, char** argv){
    static struct option long_opts[] = {
        {"strategy", required_argument, nullptr, 's'},
        {"H",        required_argument, nullptr, 'h'},
        {"Ht",       required_argument, nullptr, 't'},
        {"E",        required_argument, nullptr, 'e'},
        {"Et",       required_argument, nullptr, 'x'},
        {nullptr,0,nullptr,0}
    };
    int c;
    while((c = getopt_long(argc, argv, "s:h:t:e:x:", long_opts, nullptr)) != -1){
        switch(c){
          case 's': gStrategy = std::atoi(optarg);   break;
          case 'h': hasH       = std::atoi(optarg);   break;
          case 't': hasHt      = std::atoi(optarg);   break;
          case 'e': hasE       = std::atoi(optarg);   break;
          case 'x': hasEt      = std::atoi(optarg);   break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy <1|2> -H <0|1> -Ht <0|1>"
                     <<" -E <0|1> -Et <0|1>\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2){
        std::cerr<<"Invalid strategy "<<gStrategy<<"\n";
        std::exit(1);
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// build_par_list(): list of Im–part parameters for Minuit
static std::vector<std::string> parNamesIm;
void build_par_list(){
    parNamesIm.clear();
    parNamesIm.push_back("renormImag");
    if(hasH){
        parNamesIm.insert(parNamesIm.end(),
            {"alpha0_H","alpha1_H","n_H","b_H","M2_H","P_H"});
    }
    if(hasHt){
        parNamesIm.insert(parNamesIm.end(),
            {"alpha0_Ht","alpha1_Ht","n_Ht","b_Ht","M2_Ht","P_Ht"});
    }
    if(hasE){
        parNamesIm.insert(parNamesIm.end(),
            {"alpha0_E","alpha1_E","n_E","b_E","M2_E","P_E"});
    }
    if(hasEt){
        parNamesIm.insert(parNamesIm.end(),
            {"alpha0_Et","alpha1_Et","n_Et","b_Et","M2_Et","P_Et"});
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// fcn(): Minuit’s χ², unpack Im–parameters only, then compare
//      A_bin  vs  s1_I/c0_BH  for each bin
void fcn(int& /*npar*/, double* /*grad*/, double &f, double *par, int /*iflag*/){
    int ip = 0;
    // unpack Im–parameters
    renormImag = par[ip++];
    if(hasH){
      alpha0_H = par[ip++]; alpha1_H = par[ip++];
      n_H      = par[ip++]; b_H      = par[ip++];
      M2_H     = par[ip++]; P_H      = par[ip++];
    }
    if(hasHt){
      alpha0_Ht= par[ip++]; alpha1_Ht= par[ip++];
      n_Ht     = par[ip++]; b_Ht     = par[ip++];
      M2_Ht    = par[ip++]; P_Ht     = par[ip++];
    }
    if(hasE){
      alpha0_E = par[ip++]; alpha1_E = par[ip++];
      n_E      = par[ip++]; b_E      = par[ip++];
      M2_E     = par[ip++]; P_E      = par[ip++];
    }
    if(hasEt){
      alpha0_Et= par[ip++]; alpha1_Et= par[ip++];
      n_Et     = par[ip++]; b_Et     = par[ip++];
      M2_Et    = par[ip++]; P_Et     = par[ip++];
    }

    double chi2 = 0;
    // only Im‐fit (strategy 1 or first stage of 2)
    for(int k=0;k<Nbins;++k){
        BMK_DVCS dvcs(-1,1,0,
                      bin_Eb[k], bin_xB[k], bin_Q2[k], bin_t[k], 0.0);
        double num    = dvcs.s1_I();   // sinφ interference ∝ Im CFF
        double den    = dvcs.c0_BH();  // BH c₀ only
        double modelA = num/den;
        double resid  = (bin_A[k] - modelA)/bin_dA[k];
        chi2 += resid*resid;
    }
    f = chi2;
}

// ──────────────────────────────────────────────────────────────────────────────
// main(): orchestrate fits
int main(int argc, char** argv){
    parse_args(argc, argv);

    std::cout<<"\n=== Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<" E="<<hasE<<" Et="<<hasEt<<" ===\n";
    LoadData();
    BinBsaData();
    std::cout<<" BSA bins = "<<Nbins<<"  (from "
             <<bsaData.size()<<" raw points)\n\n";

    // ─── Stage 1: Im fit ────────────────────────────────────────────────────────
    gStage = 1;
    build_par_list();
    int nim = parNamesIm.size();
    std::vector<double> imVal(nim), imErr(nim);
    double chi2_im, edm, errdef; int nv,nx,ic, ndf_im;

    {
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);

        // define each parameter with its true default as the start value
        for(int i=0;i<nim;++i){
            const auto &nm = parNamesIm[i];
            double init, step, mn, mx;

            if(nm=="renormImag") {
                init = renormImag; step=0.01; mn=0;   mx=10;
            }
            else if(nm=="alpha0_H") {
                init = alpha0_H;  step=0.001; mn=-5;  mx=5;
            }
            else if(nm=="alpha1_H") {
                init = alpha1_H;  step=0.001; mn=-10; mx=10;
            }
            else if(nm=="n_H") {
                init = n_H;       step=0.001; mn=0;   mx=10;
            }
            else if(nm=="b_H") {
                init = b_H;       step=0.001; mn=0;   mx=10;
            }
            else if(nm=="M2_H") {
                init = M2_H;      step=0.001; mn=0;   mx=10;
            }
            else if(nm=="P_H") {
                init = P_H;       step=0.001; mn=0;   mx=10;
            }
            else if(nm=="alpha0_Ht") {
                init = alpha0_Ht; step=0.001; mn=-5;  mx=5;
            }
            else if(nm=="alpha1_Ht") {
                init = alpha1_Ht; step=0.001; mn=-10; mx=10;
            }
            else if(nm=="n_Ht") {
                init = n_Ht;      step=0.001; mn=0;   mx=10;
            }
            else if(nm=="b_Ht") {
                init = b_Ht;      step=0.001; mn=0;   mx=10;
            }
            else if(nm=="M2_Ht") {
                init = M2_Ht;     step=0.001; mn=0;   mx=10;
            }
            else if(nm=="P_Ht") {
                init = P_Ht;      step=0.001; mn=0;   mx=10;
            }
            else if(nm=="alpha0_E") {
                init = alpha0_E;  step=0.001; mn=-5;  mx=5;
            }
            else if(nm=="alpha1_E") {
                init = alpha1_E;  step=0.001; mn=-10; mx=10;
            }
            else if(nm=="n_E") {
                init = n_E;       step=0.001; mn=0;   mx=10;
            }
            else if(nm=="b_E") {
                init = b_E;       step=0.001; mn=0;   mx=10;
            }
            else if(nm=="M2_E") {
                init = M2_E;      step=0.001; mn=0;   mx=10;
            }
            else if(nm=="P_E") {
                init = P_E;       step=0.001; mn=0;   mx=10;
            }
            else if(nm=="alpha0_Et") {
                init = alpha0_Et; step=0.001; mn=-5;  mx=5;
            }
            else if(nm=="alpha1_Et") {
                init = alpha1_Et; step=0.001; mn=-10; mx=10;
            }
            else if(nm=="n_Et") {
                init = n_Et;      step=0.001; mn=0;   mx=10;
            }
            else if(nm=="b_Et") {
                init = b_Et;      step=0.001; mn=0;   mx=10;
            }
            else if(nm=="M2_Et") {
                init = M2_Et;     step=0.001; mn=0;   mx=10;
            }
            else if(nm=="P_Et") {
                init = P_Et;      step=0.001; mn=0;   mx=10;
            }
            else {
                // fallback
                init = 0; step=0.1; mn=-1e3; mx=1e3;
            }

            minu.DefineParameter(i, nm.c_str(), init, step, mn, mx);
        }

        std::cout<<"Stage1: fitting Im→BSA bins...\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_im, edm, errdef, nv, nx, ic);
        for(int i=0;i<nim;++i) minu.GetParameter(i, imVal[i], imErr[i]);
        ndf_im = Nbins - nim;
    }

    std::cout<<"\n--- Im Results ---\n";
    for(int i=0;i<nim;++i)
        std::cout<<" "<<parNamesIm[i]
                 <<" = "<<imVal[i]
                 <<" ± "<<imErr[i]<<"\n";
    std::cout<<" χ²/ndf = "<<chi2_im<<"/"<<ndf_im
             <<" = "<<(chi2_im/ndf_im)<<"\n\n";

    // ─── Stage 2: renormReal → xsec, if requested ───────────────────────────────
    if(gStrategy==2){
        gStage = 2;
        double chi2_re, ndf_re;
        double reVal, reErr;
        {
            TMinuit m2(1);
            m2.SetPrintLevel(1);
            m2.SetFCN(fcn);
            m2.DefineParameter(0,"renormReal",renormReal,0.01,0,10);
            std::cout<<"Stage2: fitting renormReal→xsec...\n";
            m2.Migrad();
            m2.Command("HESSE");
            m2.mnstat(chi2_re, edm, errdef, nv, nx, ic);
            m2.GetParameter(0, reVal, reErr);
            ndf_re = xsData.size() - 1;
        }
        std::cout<<"\n--- Re Results ---\n";
        std::cout<<" renormReal = "<<reVal
                 <<" ± "<<reErr<<"\n";
        std::cout<<" χ²/ndf = "<<chi2_re<<"/"<<ndf_re
                 <<" = "<<(chi2_re/ndf_re)<<"\n\n";
    }

    return 0;
}