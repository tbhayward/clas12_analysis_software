// fit_CFFs.C
// ──────────
// program to fit DVCS Compton Form Factors (CFFs) imaginary then real parts
//
// Strategies:
//   1) fit only Im‐parts -> BSA data
//   2) two‐step: (a) fit Im‐parts -> BSA, (b) fit Re‐parts -> xsec
//
// Usage:
//   ./fit_CFFs -strategy <1|2> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1> [--input <BSA_file>]
//
// Compile:
//   g++ -O2 fit_CFFs.C `root-config --cflags --libs` -lMinuit -o fit_CFFs

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
//   Compton Form Factor (CFF) models—imaginary parts (r_* removed, factors absorbed)
// -------------------------------------------------------------------------------------------------
extern double alpha0_H, alpha1_H, n_H, b_H, M2_H, P_H;
extern double alpha0_Ht, alpha1_Ht, n_Ht, b_Ht, M2_Ht, P_Ht;
extern double alpha0_E, alpha1_E, n_E, b_E, M2_E, P_E;
extern double alpha0_Et, alpha1_Et, n_Et, b_Et, M2_Et, P_Et;

// -------------------------------------------------------------------------------------------------
//   Compton Form Factor (CFF) models—real parts (dispersion‐relation subtraction)
// -------------------------------------------------------------------------------------------------
extern double C0_H,    MD2_H,    lambda_H;
extern double C0_Ht,   MD2_Ht,   lambda_Ht;
extern double C0_E,    MD2_E,    lambda_E;
extern double C0_Et,   MD2_Et,   lambda_Et;

// real‐ansatz globals (subtraction constants)
extern double A_H, B_H, C_H, D_H, E_H;     // note: unused in DR version
extern double A_Ht, B_Ht;                  // unused
extern double A_E, B_E, C_E, D_E;          // unused
extern double A_Et, B_Et;                  // unused

// ──────────────────────────────────────────────────────────────────────────────
// Data structures & loading
struct DataPoint { double phi,Q2,xB,t,Eb,A,sigA; };
static std::vector<DataPoint> bsaData, xsData;

void LoadData(){
    auto read = [&](const char* fn, std::vector<DataPoint>& v){
        std::ifstream in(fn);
        if(!in){ std::cerr<<"ERROR: open "<<fn<<"\n"; std::exit(1); }
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
// global controls
static int  gStrategy    = 0;       // 1 or 2
static int  gStage       = 1;       // 1=Im fit, 2=Re fit
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char* gXsFile = "imports/rga_pass1_xsec_2018.txt";

// ──────────────────────────────────────────────────────────────────────────────
// parse_args(): read -strategy, -H, -Ht, -E, -Et, [--input <BSA_file>]
// ──────────────────────────────────────────────────────────────────────────────
void parse_args(int argc, char** argv){
    static struct option long_opts[] = {
        {"strategy", required_argument, nullptr, 's'},
        {"H",        required_argument, nullptr, 'h'},
        {"Ht",       required_argument, nullptr, 't'},
        {"E",        required_argument, nullptr, 'e'},
        {"Et",       required_argument, nullptr, 'x'},
        {"input",    required_argument, nullptr, 'i'},
        {nullptr,0,nullptr,0}
    };
    int c;
    while((c = getopt_long(argc, argv, "s:h:t:e:x:i:", long_opts, nullptr)) != -1){
        switch(c){
          case 's': gStrategy = std::atoi(optarg);      break;
          case 'h': hasH       = std::atoi(optarg);      break;
          case 't': hasHt      = std::atoi(optarg);      break;
          case 'e': hasE       = std::atoi(optarg);      break;
          case 'x': hasEt      = std::atoi(optarg);      break;
          case 'i': gBsaFile   = std::string(optarg);    break;  // override BSA file
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy <1|2> -H <0|1> -Ht <0|1>"
                     <<" -E <0|1> -Et <0|1> [--input <BSA_file>]\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2){
        std::cerr<<"Invalid strategy "<<gStrategy<<"\n";
        std::exit(1);
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// build_par_list(): list of parameters for Minuit
static std::vector<std::string> parNamesIm, parNamesRe;
void build_par_list(){
    // for gStage: 1→Im, 2→Re
    if(gStage==1){
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
    else if(gStage==2){
        parNamesRe.clear();
        if(hasH){
            parNamesRe.insert(parNamesRe.end(),
                {"C0_H","MD2_H","lambda_H"});
        }
        if(hasHt){
            parNamesRe.insert(parNamesRe.end(),
                {"C0_Ht","MD2_Ht","lambda_Ht"});
        }
        if(hasE){
            parNamesRe.insert(parNamesRe.end(),
                {"C0_E","MD2_E","lambda_E"});
        }
        if(hasEt){
            parNamesRe.insert(parNamesRe.end(),
                {"C0_Et","MD2_Et","lambda_Et"});
        }
        parNamesRe.push_back("renormReal");
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// fcn(): Minuit’s χ², unpack parameters
void fcn(int& /*npar*/, double* /*grad*/, double &f, double *par, int /*iflag*/){
    int ip=0;
    if(gStage==1){
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
    }
    else if(gStage==2){
        if(hasH){
            C0_H      = par[ip++]; MD2_H     = par[ip++];
            lambda_H = par[ip++];
        }
        if(hasHt){
            C0_Ht     = par[ip++]; MD2_Ht    = par[ip++];
            lambda_Ht= par[ip++];
        }
        if(hasE){
            C0_E      = par[ip++]; MD2_E     = par[ip++];
            lambda_E = par[ip++];
        }
        if(hasEt){
            C0_Et     = par[ip++]; MD2_Et    = par[ip++];
            lambda_Et= par[ip++];
        }
        renormReal = par[ip++];
    }

    // compute χ²
    double chi2=0;
    if(gStage==1){
        for(auto &d: bsaData){
            BMK_DVCS dvcs(-1,1,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double r=(d.A - dvcs.BSA())/d.sigA;
            chi2 += r*r;
        }
    } else {
        for(auto &d: xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double r=(d.A - dvcs.CrossSection())/d.sigA;
            chi2 += r*r;
        }
    }
    f = chi2;
}

// ──────────────────────────────────────────────────────────────────────────────
// main(): orchestrate fits
int main(int argc, char** argv){
    parse_args(argc, argv);
    std::cout<<"\n=== Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<" E="<<hasE<<" Et="<<hasEt
             <<"  BSA="<<gBsaFile<<" XSEC="<<gXsFile<<" ===\n";
    LoadData();
    std::cout<<" Loaded "<<bsaData.size()<<" BSA and "
             <<xsData.size()<<" xsec points\n\n";

    std::vector<std::string> imNames, reNames;
    std::vector<double>     imVal,  imErr,  reVal, reErr;
    double chi2_im=0, chi2_re=0;
    int    ndf_im=0, ndf_re=0;

    // --- Stage 1: Im fit always for both strategies ---
    gStage = 1;
    build_par_list();
    int nim = parNamesIm.size();
    imNames = parNamesIm;
    imVal.assign(nim,0); imErr.assign(nim,0);
    {
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);
        for(int i=0;i<nim;++i){
            double init=1.0, step=0.01;
            auto &nm = imNames[i];
            if(nm=="alpha0_H"||nm=="alpha0_Ht"||nm=="alpha0_E")
                init = 0.43;
            if(nm=="alpha1_H"||nm=="alpha1_Ht"||nm=="alpha1_E")
                init = 0.85;
            if(nm=="n_H")   init = 2.43;
            if(nm=="n_Ht")  init = 1.68;
            if(nm=="n_E")   init = 1.215;
            minu.DefineParameter(i, nm.c_str(), init, step, -1e3, 1e3);
        }
        std::cout<<"Stage1: fitting Im→BSA...\n";
        minu.Migrad();
        double edm, errdef; int nv,nx,ic;
        for(int i=0;i<nim;++i) minu.GetParameter(i, imVal[i], imErr[i]);
        minu.mnstat(chi2_im, edm, errdef, nv, nx, ic);
        ndf_im = int(bsaData.size()) - nim;
        std::cout<<"\n--- Im Results ---\n";
        for(int i=0;i<nim;++i)
            std::cout<<" "<<imNames[i]<<" = "<<imVal[i]<<" ± "<<imErr[i]<<"\n";
        std::cout<<" χ²/ndf = "<<chi2_im<<"/"<<ndf_im
                 <<" = "<<(chi2_im/ndf_im)<<"\n\n";
    }

    if(gStrategy==2){
        // --- Stage 2: Re fit, hold Im fixed ---
        gStage = 2;
        build_par_list();
        int nre = parNamesRe.size();
        reNames = parNamesRe;
        reVal.assign(nre,0); reErr.assign(nre,0);
        TMinuit minu(nre);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);
        for(int i=0;i<nre;++i){
            double init=1.0, step=0.01;
            minu.DefineParameter(i, reNames[i].c_str(), init, step, -1e3, 1e3);
        }
        std::cout<<"Stage2: fitting Re→xsec (Im held fixed)...\n";
        minu.Migrad();
        double edm,errdef; int nv,nx,ic;
        for(int i=0;i<nre;++i) minu.GetParameter(i, reVal[i], reErr[i]);
        minu.mnstat(chi2_re, edm, errdef, nv, nx, ic);
        ndf_re = int(xsData.size()) - nre;
        std::cout<<"\n--- Re Results ---\n";
        for(int i=0;i<nre;++i)
            std::cout<<" "<<reNames[i]<<" = "<<reVal[i]<<" ± "<<reErr[i]<<"\n";
        std::cout<<" χ²/ndf = "<<chi2_re<<"/"<<ndf_re
                 <<" = "<<(chi2_re/ndf_re)<<"\n\n";
    }

    // --- Write combined results ---
    system("mkdir -p output/fit_results/");
    time_t now = time(nullptr);
    char buf[32];
    strftime(buf,sizeof(buf),"%Y%m%d_%H%M%S",localtime(&now));
    std::string outfn = std::string("output/fit_results/fit_results_") + buf + ".txt";
    std::ofstream out(outfn);
    out<<"# fit_CFFs results\n";
    out<<"timestamp   "<<buf<<"\n";
    out<<"strategy    "<<gStrategy<<"\n";
    out<<"H "<<hasH<<"  Ht "<<hasHt<<"  E "<<hasE<<"  Et "<<hasEt<<"\n";
    out<<"# parameters:";
    for(auto &n: imNames) out<<" "<<n;
    for(auto &n: reNames) out<<" "<<n;
    out<<"\n# values:\n";
    for(double v: imVal) out<<v<<" ";
    for(double v: reVal) out<<v<<" ";
    out<<"\n# errors:\n";
    for(double e: imErr) out<<e<<" ";
    for(double e: reErr) out<<e<<" ";
    out<<"\n# chi2 ndf chi2/ndf\n";
    if(gStrategy==1){
        out<<chi2_im<<" "<<ndf_im<<" "<<(chi2_im/ndf_im)<<"\n";
    } else {
        out<<chi2_re<<" "<<ndf_re<<" "<<(chi2_re/ndf_re)<<"\n";
    }
    out.close();
    std::cout<<"Wrote fit results to "<<outfn<<"\n";
    return 0;
}