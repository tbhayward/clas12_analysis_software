/**
 * fit_CFFs.C
 * ──────────
 * program to fit DVCS Compton Form Factors (CFFs)
 *
 * Strategies:
 *   1) fit only Im–parts -> BSA data
 *   2) fit only renormReal -> xsec data
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

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"

// bring in DVCS_xsec.C (defines BMK_DVCS, all GetIm*/GetRe*, and globals)
#include "DVCS_xsec.C"

// flags for which CFFs to include
extern bool   hasH, hasHt, hasE, hasEt;
// overall normalization factors
extern double renormImag, renormReal;

// shape parameters (all defined in DVCS_xsec.C)
extern double alpha0_H, alpha1_H, n_H,   b_H,   Mm2_H,  P_H;
extern double alpha0_Ht,alpha1_Ht,n_Ht, b_Ht,  Mm2_Ht, P_Ht;
extern double alpha0_E, alpha1_E, n_E,   b_E,   Mm2_E,  P_E;
extern double alpha0_Et,alpha1_Et,n_Et, b_Et,  Mm2_Et, P_Et;

// ──────────────────────────────────────────────────────────────────────────────
// which strategy to run
static int gStrategy = 0;  // 1 or 2

// one data‐point: φ, Q², xB, t, Eb, observable A, uncertainty σA
struct DataPoint { double phi, Q2, xB, t, Eb, A, sigA; };
static std::vector<DataPoint> bsaData, xsData;

// load BSA + xsec files
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
    read("imports/rga_prl_bsa.txt",       bsaData);
    read("imports/rga_pass1_xsec_2018.txt", xsData);
}

// parse command‐line
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
    while((c=getopt_long(argc,argv,"s:h:t:e:x:",long_opts,nullptr))!=-1){
        switch(c){
          case 's': gStrategy = std::atoi(optarg); break;
          case 'h': hasH       = std::atoi(optarg); break;
          case 't': hasHt      = std::atoi(optarg); break;
          case 'e': hasE       = std::atoi(optarg); break;
          case 'x': hasEt      = std::atoi(optarg); break;
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

// for strategy 1: we’ll fit {renormImag, α0,α1,n,b,Mm2,P} for each selected Im–CFF
static std::vector<std::string> parNames;

// build parameter list for the “Im → BSA” fit
void build_par_list(){
    parNames.clear();
    parNames.push_back("renormImag");
    if(hasH)   parNames.insert(parNames.end(),{"alpha0_H","alpha1_H","n_H","b_H","Mm2_H","P_H"});
    if(hasHt)  parNames.insert(parNames.end(),{"alpha0_Ht","alpha1_Ht","n_Ht","b_Ht","Mm2_Ht","P_Ht"});
    if(hasE)   parNames.insert(parNames.end(),{"alpha0_E","alpha1_E","n_E","b_E","Mm2_E","P_E"});
    if(hasEt)  parNames.insert(parNames.end(),{"alpha0_Et","alpha1_Et","n_Et","b_Et","Mm2_Et","P_Et"});
}

// χ² function for Minuit
void fcn(int& /*npar*/, double* /*grad*/, double &f, double *par, int /*iflag*/){
    int ip=0;
    if(gStrategy==1){
        renormImag = par[ip++];
        if(hasH){   alpha0_H=par[ip++]; alpha1_H=par[ip++]; n_H=par[ip++]; b_H=par[ip++]; Mm2_H=par[ip++]; P_H=par[ip++]; }
        if(hasHt){  alpha0_Ht=par[ip++];alpha1_Ht=par[ip++];n_Ht=par[ip++];b_Ht=par[ip++];Mm2_Ht=par[ip++];P_Ht=par[ip++]; }
        if(hasE){   alpha0_E=par[ip++]; alpha1_E=par[ip++]; n_E=par[ip++]; b_E=par[ip++]; Mm2_E=par[ip++]; P_E=par[ip++]; }
        if(hasEt){  alpha0_Et=par[ip++];alpha1_Et=par[ip++];n_Et=par[ip++];b_Et=par[ip++];Mm2_Et=par[ip++];P_Et=par[ip++]; }
        double chi2=0;
        for(auto &d: bsaData){
            BMK_DVCS dvcs(-1,1,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double r=(d.A - dvcs.BSA())/d.sigA;
            chi2 += r*r;
        }
        f = chi2;
    }
    else { // gStrategy==2
        renormReal = par[ip++];
        double chi2=0;
        for(auto &d: xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double r=(d.A - dvcs.CrossSection())/d.sigA;
            chi2 += r*r;
        }
        f = chi2;
    }
}

int main(int argc, char** argv){
    parse_args(argc,argv);
    std::cout<<"\n=== Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<" E="<<hasE<<" Et="<<hasEt<<" ===\n";
    LoadData();
    std::cout<<" Loaded "<<bsaData.size()<<" BSA and "
             <<xsData.size()<<" xsec points\n\n";

    if(gStrategy==1){
        // --- Strategy 1: fit Im → BSA ---
        build_par_list();
        int npar = parNames.size();
        std::vector<double> val(npar), err(npar);
        double chi2, edm, errdef; int nv,nx,ic;
        TMinuit minu(npar);
        minu.SetPrintLevel(1);
        minu.SetErrorDef(1.0);
        minu.SetFCN(fcn);
        for(int i=0;i<npar;++i){
            double init=1.0, step=0.01;
            const auto &nm = parNames[i];
            if(nm=="alpha0_H"||nm=="alpha0_Ht"||nm=="alpha0_E"||nm=="alpha0_Et") init=0.43;
            if(nm=="alpha1_H"||nm=="alpha1_Ht"||nm=="alpha1_E"||nm=="alpha1_Et") init=0.85;
            if(nm=="n_H")   init=n_H;
            if(nm=="n_Ht")  init=n_Ht;
            if(nm=="n_E")   init=n_E;
            if(nm=="n_Et")  init=n_Et;
            minu.DefineParameter(i, nm.c_str(), init, step, -1e3, 1e3);
        }
        std::cout<<"Stage1: fitting Im→BSA...\n";
        minu.Migrad();
        minu.mnstat(chi2,edm,errdef,nv,nx,ic);
        for(int i=0;i<npar;++i) minu.GetParameter(i, val[i], err[i]);
        int ndf = int(bsaData.size()) - npar;
        std::cout<<"\n--- Im Results ---\n";
        for(int i=0;i<npar;++i)
            std::cout<<" "<<parNames[i]<<" = "<<val[i]<<" ± "<<err[i]<<"\n";
        std::cout<<" χ²/ndf = "<<chi2<<"/"<<ndf
                 <<" = "<<(chi2/ndf)<<"\n\n";
    }
    else {
        // --- Strategy 2: fit renormReal → xsec ---
        std::vector<double> val(1), err(1);
        double chi2, edm, errdef; int nv,nx,ic;
        TMinuit minu(1);
        minu.SetPrintLevel(1);
        minu.SetErrorDef(1.0);
        minu.SetFCN(fcn);
        minu.DefineParameter(0, "renormReal", renormReal, 0.1, 0, 10);
        std::cout<<"Strategy 2: fitting renormReal→xsec...\n";
        minu.Migrad();
        minu.mnstat(chi2,edm,errdef,nv,nx,ic);
        minu.GetParameter(0, val[0], err[0]);
        int ndf = int(xsData.size()) - 1;
        std::cout<<"\n--- Re Results ---\n";
        std::cout<<" renormReal = "<<val[0]<<" ± "<<err[0]<<"\n";
        std::cout<<" χ²/ndf = "<<chi2<<"/"<<ndf
                 <<" = "<<(chi2/ndf)<<"\n\n";
    }

    return 0;
}