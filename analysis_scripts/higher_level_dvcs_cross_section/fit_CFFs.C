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

// imaginary‐part model parameters
extern double alpha0_H,  alpha1_H,  n_H,   b_H,   M2_H,  P_H;
extern double alpha0_Ht, alpha1_Ht, n_Ht,  b_Ht,  M2_Ht, P_Ht;
extern double alpha0_E,  alpha1_E,  n_E,   b_E,   M2_E,  P_E;
extern double alpha0_Et, alpha1_Et, n_Et,  b_Et,  M2_Et, P_Et;

// ──────────────────────────────────────────────────────────────────────────────
// control flags (must be before LoadData)
static int  gStrategy = 0;    // 1 or 2
static int  gStage    = 1;    // 1 = Im‐fit, 2 = Re‐fit
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char* gXsFile   = "imports/rga_pass1_xsec_2018.txt";

// ──────────────────────────────────────────────────────────────────────────────
// raw data + binned observables
struct DataPoint { double phi, Q2, xB, t, Eb, A, sigA; };
static std::vector<DataPoint> bsaData, xsData;
static std::vector<double> bin_xB, bin_Q2, bin_t, bin_Eb;
static std::vector<double> bin_A, bin_dA;
static int Nbins = 0;

// ──────────────────────────────────────────────────────────────────────────────
// Load raw BSA and XSC data
void LoadData(){
    auto read = [&](const char* fn, auto &v){
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
// Bin the BSA data by phi‐drop bins of (xB,Q2,t) and extract
// weighted sinφ amplitudes A_bin ± dA_bin
void BinBsaData(){
    bin_xB.clear(); bin_Q2.clear(); bin_t.clear(); bin_Eb.clear();
    bin_A.clear();  bin_dA.clear();
    if(bsaData.empty()) return;

    size_t start = 0;
    auto flush = [&](size_t end){
        double SwA=0, Sw2=0, Sx=0, Sq=0, St=0, Se=0;
        int M = end - start;
        for(size_t i=start;i<end;++i){
            auto &d = bsaData[i];
            double rad = d.phi * TMath::Pi()/180.0;
            double s = std::sin(rad), w = 1.0/(d.sigA*d.sigA);
            SwA += w * d.A * s;
            Sw2 += w * s*s;
            Sx += d.xB; Sq += d.Q2; St += d.t; Se += d.Eb;
        }
        bin_A .push_back(SwA/Sw2);
        bin_dA.push_back(1.0/std::sqrt(Sw2));
        bin_xB.push_back(Sx/M);
        bin_Q2.push_back(Sq/M);
        bin_t .push_back(St/M);
        bin_Eb.push_back(Se/M);
    };

    for(size_t i=1;i<=bsaData.size();++i){
        bool newbin = (i==bsaData.size() || bsaData[i].phi < bsaData[i-1].phi);
        if(newbin){
            flush(i);
            start = i;
        }
    }
    Nbins = bin_A.size();
}

// ──────────────────────────────────────────────────────────────────────────────
// parse_args(): read --strategy, -H, -Ht, -E, -Et
void parse_args(int argc, char** argv){
    static struct option opts[] = {
        {"strategy", required_argument, nullptr, 's'},
        {"H",        required_argument, nullptr, 'h'},
        {"Ht",       required_argument, nullptr, 't'},
        {"E",        required_argument, nullptr, 'e'},
        {"Et",       required_argument, nullptr, 'x'},
        {nullptr,0,nullptr,0}
    };
    int c;
    while((c=getopt_long(argc,argv,"s:h:t:e:x:",opts,nullptr))!=-1){
        switch(c){
          case 's': gStrategy=std::atoi(optarg); break;
          case 'h': hasH     =std::atoi(optarg); break;
          case 't': hasHt    =std::atoi(optarg); break;
          case 'e': hasE     =std::atoi(optarg); break;
          case 'x': hasEt    =std::atoi(optarg); break;
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
// build_par_list(): which Im‐parameters to fit (only n_*, b_*, M2_* ≥0)
static std::vector<std::string> parNamesIm;
void build_par_list(){
    parNamesIm.clear();
    if(hasH )  parNamesIm.insert(parNamesIm.end(), {"n_H","b_H","M2_H"});
    if(hasHt)  parNamesIm.insert(parNamesIm.end(), {"n_Ht","b_Ht","M2_Ht"});
    if(hasE )  parNamesIm.insert(parNamesIm.end(), {"n_E","b_E","M2_E"});
    if(hasEt)  parNamesIm.insert(parNamesIm.end(), {"n_Et","b_Et","M2_Et"});
    // renormImag, alpha0_*, alpha1_* and P_* are now fixed constants
}

// ──────────────────────────────────────────────────────────────────────────────
// fcn(): Minuit’s χ², either Im‐fit (stage1) or renormReal‐fit (stage2)
void fcn(int&, double*, double &f, double *par, int){
    int ip = 0;
    if(gStage==1){
        // ** force fixed values **
        renormImag = 1.0;
        alpha0_H   = alpha0_Ht   = alpha0_E   = alpha0_Et   = 0.43;
        alpha1_H   = alpha1_Ht   = alpha1_E   = alpha1_Et   = 0.85;
        P_H = P_Ht = P_E = P_Et = 1.0;

        // now read only the floated n, b, M2 in exactly this order
        if(hasH ){
          n_H  = par[ip++];
          b_H  = par[ip++];
          M2_H = par[ip++];
        }
        if(hasHt){
          n_Ht  = par[ip++];
          b_Ht  = par[ip++];
          M2_Ht = par[ip++];
        }
        if(hasE ){
          n_E  = par[ip++];
          b_E  = par[ip++];
          M2_E = par[ip++];
        }
        if(hasEt){
          n_Et  = par[ip++];
          b_Et  = par[ip++];
          M2_Et = par[ip++];
        }

        double chi2 = 0;
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0,
                          bin_Eb[k], bin_xB[k], bin_Q2[k], bin_t[k], 0.0);
            double modelA = dvcs.s1_I()/dvcs.c0_BH();
            double r = (bin_A[k] - modelA)/bin_dA[k];
            chi2 += r*r;
        }
        f = chi2;
    }
    else {
        // renormReal fit in stage 2
        renormReal = par[ip++];
        double chi2 = 0;
        for(auto &d: xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double m  = dvcs.CrossSection();
            double r  = (d.A - renormReal*m)/d.sigA;
            chi2 += r*r;
        }
        f = chi2;
    }
}

// ──────────────────────────────────────────────────────────────────────────────
int main(int argc, char** argv){
    parse_args(argc,argv);

    std::cout<<"\n=== Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<" E="<<hasE<<" Et="<<hasEt<<" ===\n";

    LoadData();
    BinBsaData();
    std::cout<<" BSA bins = "<<Nbins<<"  (from "
             <<bsaData.size()<<" raw points)\n\n";

    // enforce all P_* = 1.0 up front
    P_H = P_Ht = P_E = P_Et = 1.0;

    // prepare result vectors
    std::vector<std::string> outNames;
    std::vector<double>      outVal, outErr;
    double finalChi2=0;
    int    finalNdf =0;

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

        // define each Im parameter with bounds [0,1e3]
        for(int i=0;i<nim;++i){
            minu.DefineParameter(i,
                                parNamesIm[i].c_str(),
                                1.0,    // init
                                0.01,   // step
                                0.0,    // lower
                                1e3     // upper
            );
        }

        std::cout<<"Stage1: fitting Im–CFFs (n,b,M2 ≥0; fixed renormImag, α₀, α₁, P)...\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_im,edm,errdef,nv,nx,ic);
        for(int i=0;i<nim;++i) minu.GetParameter(i,imVal[i],imErr[i]);
        ndf_im = Nbins - nim;
    }

    // store Im‐fit results
    outNames  = parNamesIm;
    outVal    = imVal;
    outErr    = imErr;
    finalChi2 = chi2_im;
    finalNdf  = ndf_im;

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
        double chi2_re, edm2, errdef2; int nv2,nx2,ic2, ndf_re;
        double reVal, reErr;

        {
            TMinuit m2(1);
            m2.SetPrintLevel(1);
            m2.SetFCN(fcn);
            m2.DefineParameter(0,"renormReal",1.0,0.01,0.0,1e3);
            std::cout<<"Stage2: fitting renormReal (≥0)...\n";
            m2.Migrad();
            m2.Command("HESSE");
            m2.mnstat(chi2_re,edm2,errdef2,nv2,nx2,ic2);
            m2.GetParameter(0,reVal,reErr);
            ndf_re = xsData.size() - 1;
        }

        outNames .push_back("renormReal");
        outVal   .push_back(reVal);
        outErr   .push_back(reErr);
        finalChi2 = chi2_re;
        finalNdf  = ndf_re;

        std::cout<<"\n--- Re Results ---\n";
        std::cout<<" renormReal = "<<reVal
                 <<" ± "<<reErr<<"\n";
        std::cout<<" χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                 <<" = "<<(finalChi2/finalNdf)<<"\n\n";
    }

    // ─── Write results to timestamped file ───────────────────────────────────────
    system("mkdir -p output/fit_results/");
    time_t now = time(nullptr);
    char tb[32];
    strftime(tb,sizeof(tb),"%Y%m%d_%H%M%S",localtime(&now));
    std::string fname = "output/fit_results/fit_results_" + std::string(tb) + ".txt";
    std::ofstream out(fname);
    out<<"# fit_CFFs results\n";
    out<<"timestamp   "<<tb<<"\n";
    out<<"strategy    "<<gStrategy<<"\n";
    out<<"H "<<hasH<<"  Ht "<<hasHt<<"  E "<<hasE<<"  Et "<<hasEt<<"\n";
    out<<"# parameters:";
    for(auto &n: outNames) out<<" "<<n;
    out<<"\n# values:\n";
    for(auto &v: outVal) out<<v<<" ";
    out<<"\n# errors:\n";
    for(auto &e: outErr) out<<e<<" ";
    out<<"\n# chi2 ndf chi2/ndf\n";
    out<<finalChi2<<" "<<finalNdf<<" "<<(finalChi2/finalNdf)<<"\n";
    out.close();

    std::cout<<"Wrote fit results to "<<fname<<"\n";
    return 0;
}