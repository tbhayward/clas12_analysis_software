/**
 * fit_CFFs.C
 * ──────────
 * program to fit DVCS Compton Form Factors (CFFs) imaginary then real parts
 *
 * Strategies:
 *   1) fit only Im–parts -> BSA data
 *   2) two-step: (a) fit Im–parts -> BSA, (b) fit Re–parts -> xsec
 *
 * Usage:
 *   ./fit_CFFs --strategy <1|2> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1> [--input <BSA_file>]
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
extern double alpha0_H,  alpha1_H,  n_H,   b_H,   M2_H,  P_H;
extern double alpha0_Ht, alpha1_Ht, n_Ht,  b_Ht,  M2_Ht, P_Ht;
extern double alpha0_E,  alpha1_E,  n_E,   b_E,   M2_E,  P_E;
extern double alpha0_Et, alpha1_Et, n_Et,  b_Et,  M2_Et, P_Et;

// -------------------------------------------------------------------------------------------------
//   Compton Form Factor (CFF) models—real parts (dispersion‐relation subtraction)
// -------------------------------------------------------------------------------------------------
extern double C0_H,    MD2_H,    lambda_H;
extern double C0_Ht,   MD2_Ht,   lambda_Ht;
extern double C0_E,    MD2_E,    lambda_E;
extern double C0_Et,   MD2_Et,   lambda_Et;

// ──────────────────────────────────────────────────────────────────────────────
// global controls (must be before LoadData)
static int  gStrategy    = 0;       // 1 or 2
static int  gStage       = 1;       // 1=Im fit, 2=Re fit
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char*  gXsFile   = "imports/rga_pass1_xsec_2018.txt";

// ──────────────────────────────────────────────────────────────────────────────
// Data structures & loading
struct DataPoint { double phi,Q2,xB,t,Eb,A,sigA; };
static std::vector<DataPoint> bsaData, xsData;

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
// parse_args(): read --strategy, -H, -Ht, -E, -Et, [--input <BSA_file>]
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
          case 's': gStrategy = std::atoi(optarg);   break;
          case 'h': hasH       = std::atoi(optarg);   break;
          case 't': hasHt      = std::atoi(optarg);   break;
          case 'e': hasE       = std::atoi(optarg);   break;
          case 'x': hasEt      = std::atoi(optarg);   break;
          case 'i': gBsaFile   = std::string(optarg); break;
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
    else { // gStage==2
        if(hasH){
            C0_H      = par[ip++]; MD2_H     = par[ip++];
            lambda_H  = par[ip++];
        }
        if(hasHt){
            C0_Ht     = par[ip++]; MD2_Ht    = par[ip++];
            lambda_Ht = par[ip++];
        }
        if(hasE){
            C0_E      = par[ip++]; MD2_E     = par[ip++];
            lambda_E  = par[ip++];
        }
        if(hasEt){
            C0_Et     = par[ip++]; MD2_Et    = par[ip++];
            lambda_Et = par[ip++];
        }
        renormReal = par[ip++];
    }

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

    // ─── Stage 1: Im fit ────────────────────────────────────────────────────────
    gStage = 1;
    build_par_list();
    int nim = parNamesIm.size();
    std::vector<double> imVal(nim,0), imErr(nim,0);
    double chi2_im, edm, errdef; int nv,nx,ic, ndf_im;
    {
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);
        for(int i=0;i<nim;++i){
            const auto &nm = parNamesIm[i];
            double init = 1.0, step = 0.01;
            if(nm=="renormImag")   init = 1.0;
            else if(nm=="alpha0_H") init = alpha0_H;
            else if(nm=="alpha1_H") init = alpha1_H;
            else if(nm=="n_H")      init = n_H;
            else if(nm=="b_H")      init = b_H;
            else if(nm=="M2_H")     init = M2_H;
            else if(nm=="P_H")      init = P_H;
            else if(nm=="alpha0_Ht") init = alpha0_Ht;
            else if(nm=="alpha1_Ht") init = alpha1_Ht;
            else if(nm=="n_Ht")      init = n_Ht;
            else if(nm=="b_Ht")      init = b_Ht;
            else if(nm=="M2_Ht")     init = M2_Ht;
            else if(nm=="P_Ht")      init = P_Ht;
            else if(nm=="alpha0_E")  init = alpha0_E;
            else if(nm=="alpha1_E")  init = alpha1_E;
            else if(nm=="n_E")       init = n_E;
            else if(nm=="b_E")       init = b_E;
            else if(nm=="M2_E")      init = M2_E;
            else if(nm=="P_E")       init = P_E;
            else if(nm=="alpha0_Et") init = alpha0_Et;
            else if(nm=="alpha1_Et") init = alpha1_Et;
            else if(nm=="n_Et")      init = n_Et;
            else if(nm=="b_Et")      init = b_Et;
            else if(nm=="M2_Et")     init = M2_Et;
            else if(nm=="P_Et")      init = P_Et;
            minu.DefineParameter(i, nm.c_str(), init, step, -1e3, 1e3);
        }
        std::cout<<"Stage1: fitting Im→BSA...\n";
        minu.Migrad();
        minu.mnstat(chi2_im, edm, errdef, nv, nx, ic);
        for(int i=0;i<nim;++i) minu.GetParameter(i, imVal[i], imErr[i]);
        ndf_im = int(bsaData.size()) - nim;
    }
    std::cout<<"\n--- Im Results ---\n";
    for(int i=0;i<nim;++i)
        std::cout<<" "<<parNamesIm[i]<<" = "<<imVal[i]<<" ± "<<imErr[i]<<"\n";
    std::cout<<" χ²/ndf = "<<chi2_im<<"/"<<ndf_im
             <<" = "<<(chi2_im/ndf_im)<<"\n\n";

    // ─── Stage 2: Re fit (strategy 2) ────────────────────────────────────────────
    if(gStrategy==2){
        gStage = 2;
        build_par_list();
        int nre = parNamesRe.size();
        std::vector<double> reVal(nre,0), reErr(nre,0);
        double chi2_re; int ndf_re;
        {
            TMinuit minu2(nre);
            minu2.SetPrintLevel(1);
            minu2.SetFCN(fcn);
            for(int i=0;i<nre;++i){
                const auto &nm = parNamesRe[i];
                double init = 1.0, step = 0.01;
                if(nm=="C0_H")       init = C0_H;
                else if(nm=="MD2_H")    init = MD2_H;
                else if(nm=="lambda_H") init = lambda_H;
                else if(nm=="C0_Ht")    init = C0_Ht;
                else if(nm=="MD2_Ht")   init = MD2_Ht;
                else if(nm=="lambda_Ht")init = lambda_Ht;
                else if(nm=="C0_E")     init = C0_E;
                else if(nm=="MD2_E")    init = MD2_E;
                else if(nm=="lambda_E") init = lambda_E;
                else if(nm=="C0_Et")    init = C0_Et;
                else if(nm=="MD2_Et")   init = MD2_Et;
                else if(nm=="lambda_Et")init = lambda_Et;
                else if(nm=="renormReal") init = renormReal;
                minu2.DefineParameter(i, nm.c_str(), init, step, -1e3, 1e3);
            }
            std::cout<<"Stage2: fitting Re→xsec (Im fixed)...\n";
            minu2.Migrad();
            minu2.mnstat(chi2_re, edm, errdef, nv, nx, ic);
            for(int i=0;i<nre;++i) minu2.GetParameter(i, reVal[i], reErr[i]);
            ndf_re = int(xsData.size()) - nre;
        }
        std::cout<<"\n--- Re Results ---\n";
        for(int i=0;i<nre;++i)
            std::cout<<" "<<parNamesRe[i]<<" = "<<reVal[i]<<" ± "<<reErr[i]<<"\n";
        std::cout<<" χ²/ndf = "<<chi2_re<<"/"<<ndf_re
                 <<" = "<<(chi2_re/ndf_re)<<"\n\n";

        // write both Im & Re results
        system("mkdir -p output/fit_results/");
        time_t now = time(nullptr);
        char buf[32];
        strftime(buf,sizeof(buf),"%Y%m%d_%H%M%S",localtime(&now));
        std::ofstream out(std::string("output/fit_results/fit_results_")+buf+".txt");
        out<<"# fit_CFFs results\n";
        out<<"timestamp "<<buf<<"\n";
        out<<"strategy "<<gStrategy<<"\n";
        out<<"H "<<hasH<<" Ht "<<hasHt<<" E "<<hasE<<" Et "<<hasEt<<"\n";
        out<<"# parameters:";
        for(auto &n: parNamesIm) out<<" "<<n;
        for(auto &n: parNamesRe) out<<" "<<n;
        out<<"\n# values:\n";
        for(auto &v: imVal) out<<v<<" ";
        for(auto &v: reVal) out<<v<<" ";
        out<<"\n# errors:\n";
        for(auto &e: imErr) out<<e<<" ";
        for(auto &e: reErr) out<<e<<" ";
        out<<"\n";
        out.close();
        std::cout<<"Wrote fit results to output/fit_results/fit_results_"<<buf<<".txt\n";
    }

    return 0;
}