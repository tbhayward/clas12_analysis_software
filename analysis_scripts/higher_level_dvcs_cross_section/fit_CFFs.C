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
 *   ./fit_CFFs --strategy <1|2> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1> [--constraint <0|1>]
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
#include <map>               // for valMap / errMap

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"

// pull in full BMK_DVCS + CFF code, with globals
#include "DVCS_xsec.C"

// extern flags & renormalizations
extern bool   hasH, hasHt, hasE, hasEt;
extern double renormImag, renormReal;

// imaginary-part model parameters
extern double alpha0_H,  alpha1_H,  n_H,   b_H,   M2_H,  P_H;
extern double alpha0_Ht, alpha1_Ht, n_Ht,  b_Ht,  M2_Ht, P_Ht;
extern double alpha0_E,  alpha1_E,  n_E,   b_E,   M2_E,  P_E;
extern double alpha0_Et, alpha1_Et, n_Et,  b_Et,  M2_Et, P_Et;

// ──────────────────────────────────────────────────────────────────────────────
// control flags
static int  gStrategy   = 0;    // 1 or 2
static int  gStage      = 1;    // 1 = Im-fit, 2 = Re-fit
static int  gConstraint = 0;    // 0 = no cut, 1 = apply -t/Q2 < 0.2
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char* gXsFile  = "imports/rga_pass1_xsec_2018.txt";

// ──────────────────────────────────────────────────────────────────────────────
// raw data + binned observables
struct DataPoint { double phi, Q2, xB, t, Eb, A, sigA; };
static std::vector<DataPoint> bsaData, xsData;
static std::vector<double> bin_xB, bin_Q2, bin_t, bin_Eb;
static std::vector<double> bin_A, bin_dA;
static int Nbins = 0;

// ──────────────────────────────────────────────────────────────────────────────
// Load raw BSA and XSC data, optionally applying the -t/Q2 constraint
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
            if(gConstraint==1){
                // skip any single point with -t/Q2 ≥ 0.2
                if( (-d.t / d.Q2) >= 0.2 ) continue;
            }
            v.push_back(d);
        }
    };
    read(gBsaFile.c_str(), bsaData);
    read(gXsFile,          xsData);
}

// ──────────────────────────────────────────────────────────────────────────────
// Bin the BSA data by phi-drop bins and extract sinφ amplitudes
void BinBsaData(){
    bin_xB.clear(); bin_Q2.clear(); bin_t.clear(); bin_Eb.clear();
    bin_A .clear(); bin_dA.clear();
    if(bsaData.empty()) return;

    size_t start = 0;
    auto flush = [&](size_t end){
        double SwA=0, Sw2=0, Sx=0, Sq=0, St=0, Se=0;
        int M = end - start;
        for(size_t i=start;i<end;++i){
            auto &d = bsaData[i];
            double rad = d.phi * TMath::Pi()/180.0;
            double s   = std::sin(rad);
            double w   = 1.0/(d.sigA*d.sigA);
            SwA += w*d.A*s;
            Sw2 += w*s*s;
            Sx  += d.xB;  Sq += d.Q2;
            St  += d.t;   Se += d.Eb;
        }
        bin_A .push_back(SwA/Sw2);
        bin_dA.push_back(1.0/std::sqrt(Sw2));
        bin_xB.push_back(Sx/M);
        bin_Q2.push_back(Sq/M);
        bin_t .push_back(St/M);
        bin_Eb.push_back(Se/M);
    };

    for(size_t i=1;i<=bsaData.size();++i){
        bool newbin = (i==bsaData.size()
                       || bsaData[i].phi < bsaData[i-1].phi);
        if(newbin){
            flush(i);
            start = i;
        }
    }
    Nbins = bin_A.size();
}

// ──────────────────────────────────────────────────────────────────────────────
// parse_args(): read --strategy, -H, -Ht, -E, -Et, [--constraint]
void parse_args(int argc, char** argv){
    static struct option opts[] = {
        {"strategy",   required_argument, nullptr, 's'},
        {"H",          required_argument, nullptr, 'h'},
        {"Ht",         required_argument, nullptr, 't'},
        {"E",          required_argument, nullptr, 'e'},
        {"Et",         required_argument, nullptr, 'x'},
        {"constraint", required_argument, nullptr, 'C'},
        {nullptr,0,nullptr,0}
    };
    int c;
    while((c=getopt_long(argc,argv,"s:h:t:e:x:C:",opts,nullptr))!=-1){
        switch(c){
          case 's': gStrategy   = std::atoi(optarg); break;
          case 'h': hasH        = std::atoi(optarg); break;
          case 't': hasHt       = std::atoi(optarg); break;
          case 'e': hasE        = std::atoi(optarg); break;
          case 'x': hasEt       = std::atoi(optarg); break;
          case 'C': gConstraint = std::atoi(optarg); break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy <1|2> -H <0|1> -Ht <0|1>"
                     <<" -E <0|1> -Et <0|1> [--constraint <0|1>]\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2){
        std::cerr<<"Invalid strategy "<<gStrategy<<"\n"; std::exit(1);
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// build_par_list(): which Im-parameters to fit
static std::vector<std::string> parNamesIm;
void build_par_list(){
    parNamesIm.clear();
    parNamesIm.push_back("renormImag");
    if(hasH )  parNamesIm.insert(parNamesIm.end(),
                   {"alpha0_H","alpha1_H","n_H","b_H","M2_H","P_H"});
    if(hasHt)  parNamesIm.insert(parNamesIm.end(),
                   {"alpha0_Ht","alpha1_Ht","n_Ht","b_Ht","M2_Ht","P_Ht"});
    if(hasE )  parNamesIm.insert(parNamesIm.end(),
                   {"alpha0_E","alpha1_E","n_E","b_E","M2_E","P_E"});
    if(hasEt)  parNamesIm.insert(parNamesIm.end(),
                   {"alpha0_Et","alpha1_Et","n_Et","b_Et","M2_Et","P_Et"});
}

// ──────────────────────────────────────────────────────────────────────────────
// fcn(): Minuit’s χ² for Im-fit (gStage=1) or renormReal-fit (gStage=2)
void fcn(int&, double*, double &f, double *par, int){
    int ip=0;
    if(gStage==1){
        renormImag = par[ip++];
        if(hasH ){
          alpha0_H = par[ip++]; alpha1_H = par[ip++];
          n_H      = par[ip++]; b_H      = par[ip++];
          M2_H     = par[ip++]; P_H      = par[ip++];
        }
        if(hasHt){
          alpha0_Ht = par[ip++]; alpha1_Ht = par[ip++];
          n_Ht      = par[ip++]; b_Ht      = par[ip++];
          M2_Ht     = par[ip++]; P_Ht      = par[ip++];
        }
        if(hasE ){
          alpha0_E = par[ip++]; alpha1_E = par[ip++];
          n_E      = par[ip++]; b_E      = par[ip++];
          M2_E     = par[ip++]; P_E      = par[ip++];
        }
        if(hasEt){
          alpha0_Et = par[ip++]; alpha1_Et = par[ip++];
          n_Et      = par[ip++]; b_Et      = par[ip++];
          M2_Et     = par[ip++]; P_Et      = par[ip++];
        }
        double chi2=0;
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0,
                          bin_Eb[k],bin_xB[k],bin_Q2[k],bin_t[k],0.0);
            double modelA = dvcs.s1_I()/dvcs.c0_BH();
            double r       = (bin_A[k]-modelA)/bin_dA[k];
            chi2 += r*r;
        }
        f = chi2;
    }
    else {
        renormReal = par[ip++];
        double chi2=0;
        for(auto &d: xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double m = dvcs.CrossSection();
            double r = (d.A - renormReal*m)/d.sigA;
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
             <<" E="<<hasE<<" Et="<<hasEt
             <<"  constraint="<<gConstraint<<" ===\n";

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

        // define each Im parameter, enforce M2_* ≥ 0
        for(int i=0;i<nim;++i){
            const auto &nm = parNamesIm[i];
            double init=0, step=0.01;
            // pull starting values from extern globals
            #define GETINIT(NAME) if(nm==#NAME) init = NAME;
            GETINIT(renormImag)
            GETINIT(alpha0_H)   GETINIT(alpha1_H)
            GETINIT(n_H)        GETINIT(b_H)
            GETINIT(M2_H)       GETINIT(P_H)
            GETINIT(alpha0_Ht)  GETINIT(alpha1_Ht)
            GETINIT(n_Ht)       GETINIT(b_Ht)
            GETINIT(M2_Ht)      GETINIT(P_Ht)
            GETINIT(alpha0_E)   GETINIT(alpha1_E)
            GETINIT(n_E)        GETINIT(b_E)
            GETINIT(M2_E)       GETINIT(P_E)
            GETINIT(alpha0_Et)  GETINIT(alpha1_Et)
            GETINIT(n_Et)       GETINIT(b_Et)
            GETINIT(M2_Et)      GETINIT(P_Et)
            #undef GETINIT

            double lower=-1e3, upper=1e3;
            if(nm.rfind("M2_",0)==0) lower = 0.0;   // M2 must stay ≥0

            minu.DefineParameter(i, nm.c_str(), init, step, lower, upper);
            // renormImag remains fixed? If yes, uncomment next line:
            // if(nm=="renormImag") minu.FixParameter(i);
        }

        std::cout<<"Stage1: fitting Im–CFFs (all shape + P floated; renormImag floated)...\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_im,edm,errdef,nv,nx,ic);
        for(int i=0;i<nim;++i)
            minu.GetParameter(i, imVal[i], imErr[i]);
        ndf_im = Nbins - nim;
    }

    // ─── Collect all Im‐parameters (float or fixed) into maps ─────────────
    std::map<std::string,double> valMap, errMap;
    for(int i=0;i<nim;++i){
        const auto &name = parNamesIm[i];
        valMap[name] = imVal[i];
        errMap[name] = imErr[i];
    }
    auto FIXED = [&](auto& var, const char* nm){
        if(!valMap.count(nm)){
            valMap[nm] = var;
            errMap[nm] = 0.0;
        }
    };
    FIXED(renormImag,  "renormImag");
    FIXED(alpha0_H,    "alpha0_H");    FIXED(alpha1_H,    "alpha1_H");
    FIXED(n_H,         "n_H");         FIXED(b_H,         "b_H");
    FIXED(M2_H,        "M2_H");        FIXED(P_H,         "P_H");
    FIXED(alpha0_Ht,   "alpha0_Ht");   FIXED(alpha1_Ht,   "alpha1_Ht");
    FIXED(n_Ht,        "n_Ht");        FIXED(b_Ht,        "b_Ht");
    FIXED(M2_Ht,       "M2_Ht");       FIXED(P_Ht,        "P_Ht");
    FIXED(alpha0_E,    "alpha0_E");    FIXED(alpha1_E,    "alpha1_E");
    FIXED(n_E,         "n_E");         FIXED(b_E,         "b_E");
    FIXED(M2_E,        "M2_E");        FIXED(P_E,         "P_E");
    FIXED(alpha0_Et,   "alpha0_Et");   FIXED(alpha1_Et,   "alpha1_Et");
    FIXED(n_Et,        "n_Et");        FIXED(b_Et,        "b_Et");
    FIXED(M2_Et,       "M2_Et");       FIXED(P_Et,        "P_Et");

    double finalChi2 = chi2_im;
    int    finalNdf  = ndf_im;

    std::cout<<"\n--- Im Results ---\n";
    for(auto &kv : valMap){
        std::cout<<" "<<kv.first
                 <<" = "<<kv.second
                 <<" ± "<<errMap[kv.first]<<"\n";
    }
    std::cout<<" χ²/ndf = "<<finalChi2<<"/"<<finalNdf
             <<" = "<<(finalChi2/finalNdf)<<"\n\n";

    // ─── Stage 2: renormReal → xsec, if requested ───────────────────────────────
    double reVal=0, reErr=0;
    if(gStrategy==2){
        gStage = 2;
        double chi2_re, edm2, errdef2; int nv2,nx2,ic2, ndf_re;

        {
            TMinuit m2(1);
            m2.SetPrintLevel(1);
            m2.SetFCN(fcn);
            m2.DefineParameter(0,"renormReal",renormReal,0.01,-1e3,1e3);
            std::cout<<"Stage2: fitting renormReal...\n";
            m2.Migrad();
            m2.Command("HESSE");
            m2.mnstat(chi2_re,edm2,errdef2,nv2,nx2,ic2);
            m2.GetParameter(0,reVal,reErr);
            ndf_re = xsData.size() - 1;
        }
        valMap["renormReal"] = reVal;   errMap["renormReal"] = reErr;
        finalChi2 = chi2_re;          finalNdf = ndf_re;

        std::cout<<"\n--- Re Results ---\n";
        std::cout<<" renormReal = "<<reVal
                 <<" ± "<<reErr<<"\n";
        std::cout<<" χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                 <<" = "<<(finalChi2/finalNdf)<<"\n\n";
    }

    // ─── Write full-parameter results to timestamped file ───────────────────────
    system("mkdir -p output/fit_results/");
    time_t now = time(nullptr);
    char tb[32];
    strftime(tb,sizeof(tb),"%Y%m%d_%H%M%S",localtime(&now));
    std::string fname = "output/fit_results/fit_results_"+std::string(tb)+".txt";
    std::ofstream fout(fname);
    fout<<"# fit_CFFs results\n";
    fout<<"timestamp   "<<tb<<"\n";
    fout<<"strategy    "<<gStrategy<<"\n";
    fout<<"constraint  "<<gConstraint<<"\n";
    fout<<"H "<<hasH<<"  Ht "<<hasHt<<"  E "<<hasE<<"  Et "<<hasEt<<"\n";
    fout<<"# parameters:";
    for(auto &kv: valMap) fout<<" "<<kv.first;
    fout<<"\n# values:\n";
    for(auto &kv: valMap)   fout<<kv.second<<" ";
    fout<<"\n# errors:\n";
    for(auto &kv: errMap)   fout<<kv.second<<" ";
    fout<<"\n# chi2 ndf chi2/ndf\n";
    fout<<finalChi2<<" "<<finalNdf<<" "<<(finalChi2/finalNdf)<<"\n";
    fout.close();

    std::cout<<"Wrote fit results to "<<fname<<"\n";
    return 0;
}