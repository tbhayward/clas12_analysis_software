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
 *             [--constraint <0|1>] [--input <BSA_file>]
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
#include <map>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"

// pull in full BMK_DVCS + CFF code, with globals
#include "DVCS_xsec.C"

// extern flags & renormalizations
extern bool   hasH, hasHt, hasE, hasEt;
extern double renormImag, renormReal;

// imaginary-part model parameters
extern double r_H,      alpha0_H,  alpha1_H,  n_H,   b_H,   M2_H,  P_H;
extern double r_Ht,     alpha0_Ht, alpha1_Ht, n_Ht,  b_Ht,  M2_Ht, P_Ht;
extern double r_E,      alpha0_E,  alpha1_E,  n_E,   b_E,   M2_E,  P_E;
extern double r_Et,     alpha0_Et, alpha1_Et, n_Et,  b_Et,  M2_Et, P_Et;

// ──────────────────────────────────────────────────────────────────────────────
// control flags
static int  gStrategy   = 0;    // 1 or 2
static int  gStage      = 1;    // 1 = Im-fit, 2 = Re-fit
static int  gConstraint = 0;    // 0 = no cut, 1 = apply -t/Q2<0.2
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
// temporaries for bin-wise φ–fit
static std::vector<double> _phis, _As, _sigAs;

// Minuit FCN for A0,B fit: model = A0 * sinφ / (1 + B cosφ)
void fcn_bin(int&, double*, double &f, double *par, int){
    double A0 = par[0], B = par[1];
    double chi2 = 0;
    for(size_t i=0; i<_phis.size(); ++i){
        double model = A0 * std::sin(_phis[i]) / (1 + B * std::cos(_phis[i]));
        double diff  = _As[i] - model;
        chi2 += diff*diff / (_sigAs[i]*_sigAs[i]);
    }
    f = chi2;
}

// ──────────────────────────────────────────────────────────────────────────────
// Load raw BSA + XSC, apply constraint if requested
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
            if(gConstraint==1 && (-d.t/d.Q2) >= 0.2) continue;
            v.push_back(d);
        }
    };
    read(gBsaFile.c_str(), bsaData);
    read(gXsFile,          xsData);
}

// ──────────────────────────────────────────────────────────────────────────────
// Bin BSA by φ‐drop and **fit** A₀ sinφ/(1+B cosφ) to extract A₀ ± σA₀
void BinBsaData(){
    bin_xB.clear(); bin_Q2.clear(); bin_t.clear(); bin_Eb.clear();
    bin_A .clear(); bin_dA.clear();
    if(bsaData.empty()) return;

    size_t start = 0;
    for(size_t i=1; i<=bsaData.size(); ++i){
        bool newbin = (i==bsaData.size() || bsaData[i].phi < bsaData[i-1].phi);
        if(!newbin) continue;

        // collect this bin's raw points
        int M = i - start;
        double Sx=0,Sq=0,St=0,Se=0;
        _phis.clear(); _As.clear(); _sigAs.clear();
        for(size_t j=start; j< i; ++j){
            auto &d = bsaData[j];
            Sx += d.xB; Sq += d.Q2; St += d.t; Se += d.Eb;
            double rad = d.phi * TMath::Pi()/180.0;
            _phis.push_back(rad);
            _As.push_back(d.A);
            _sigAs.push_back(d.sigA);
        }

        //AVERAGED kinematics
        bin_xB.push_back(Sx/M);
        bin_Q2.push_back(Sq/M);
        bin_t .push_back(St/M);
        bin_Eb.push_back(Se/M);

        // now fit A₀,B via Minuit
        TMinuit minu(2);
        minu.SetPrintLevel(-1);
        minu.SetFCN(fcn_bin);

        // initial guesses: A₀ ~ max(A), B~0
        double A0_init = 0;
        for(double a:_As) if(std::abs(a)>A0_init) A0_init=std::abs(a);
        minu.DefineParameter(0, "A0", A0_init, 0.01, 0.0, 1e3);
        minu.DefineParameter(1, "B",  0.0,      0.1, -10.0,10.0);

        minu.Migrad();
        minu.Command("HESSE");

        double A0_fit, errA0, B_fit, errB;
        minu.GetParameter(0, A0_fit, errA0);
        // (we ignore B_fit here)

        bin_A .push_back(A0_fit);
        bin_dA.push_back(errA0);

        start = i;
    }
    Nbins = bin_A.size();
}

// ──────────────────────────────────────────────────────────────────────────────
// parse_args: --strategy, -H, -Ht, -E, -Et, [--constraint], [--input]
void parse_args(int argc,char**argv){
    static struct option opts[]={
        {"strategy",   required_argument,nullptr,'s'},
        {"H",          required_argument,nullptr,'h'},
        {"Ht",         required_argument,nullptr,'t'},
        {"E",          required_argument,nullptr,'e'},
        {"Et",         required_argument,nullptr,'x'},
        {"constraint", required_argument,nullptr,'C'},
        {"input",      required_argument,nullptr,'i'},
        {nullptr,0,nullptr,0}
    };
    int c;
    while((c=getopt_long(argc,argv,"s:h:t:e:x:C:i:",opts,nullptr))!=-1){
        switch(c){
          case 's': gStrategy   = atoi(optarg);           break;
          case 'h': hasH        = atoi(optarg);           break;
          case 't': hasHt       = atoi(optarg);           break;
          case 'e': hasE        = atoi(optarg);           break;
          case 'x': hasEt       = atoi(optarg);           break;
          case 'C': gConstraint = atoi(optarg);           break;
          case 'i': gBsaFile    = std::string(optarg);    break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy<1|2> -H<0|1> -Ht<0|1>"
                     <<" -E<0|1> -Et<0|1> [--constraint<0|1>]"
                     <<" [--input <BSA_file>]\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2){
        std::cerr<<"Invalid strategy\n"; std::exit(1);
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// build which Im‐parameters to fit
static std::vector<std::string> parNamesIm;
void build_par_list(){
    parNamesIm.clear();
    parNamesIm.push_back("renormImag");
    if(hasH )  parNamesIm.insert(parNamesIm.end(),
                   {"r_H","alpha0_H","alpha1_H","n_H","b_H","M2_H","P_H"});
    if(hasHt)  parNamesIm.insert(parNamesIm.end(),
                   {"r_Ht","alpha0_Ht","alpha1_Ht","n_Ht","b_Ht","M2_Ht","P_Ht"});
    if(hasE )  parNamesIm.insert(parNamesIm.end(),
                   {"r_E","alpha0_E","alpha1_E","n_E","b_E","M2_E","P_E"});
    if(hasEt)  parNamesIm.insert(parNamesIm.end(),
                   {"r_Et","alpha0_Et","alpha1_Et","n_Et","b_Et","M2_Et","P_Et"});
}

// ──────────────────────────────────────────────────────────────────────────────
// χ² function: Im‐fit (gStage=1) or renormReal‐fit (gStage=2)
void fcn(int&, double*, double &f, double *par, int){
    int ip=0;
    if(gStage==1){
        renormImag = par[ip++];
        if(hasH ){
          r_H      = par[ip++];
          alpha0_H = par[ip++]; alpha1_H = par[ip++];
          n_H      = par[ip++]; b_H      = par[ip++];
          M2_H     = par[ip++]; P_H      = par[ip++];
        }
        if(hasHt){
          r_Ht     = par[ip++];
          alpha0_Ht= par[ip++]; alpha1_Ht= par[ip++];
          n_Ht     = par[ip++]; b_Ht     = par[ip++];
          M2_Ht    = par[ip++]; P_Ht     = par[ip++];
        }
        if(hasE ){
          r_E      = par[ip++];
          alpha0_E = par[ip++]; alpha1_E = par[ip++];
          n_E      = par[ip++]; b_E      = par[ip++];
          M2_E     = par[ip++]; P_E      = par[ip++];
        }
        if(hasEt){
          r_Et     = par[ip++];
          alpha0_Et= par[ip++]; alpha1_Et= par[ip++];
          n_Et     = par[ip++]; b_Et     = par[ip++];
          M2_Et    = par[ip++]; P_Et     = par[ip++];
        }
        double chi2=0;
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0,
                          bin_Eb[k],bin_xB[k],bin_Q2[k],bin_t[k],0.0);
            double mA = dvcs.s1_I()/dvcs.c0_BH();
            double r  = (bin_A[k]-mA)/bin_dA[k];
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
int main(int argc,char**argv){
    parse_args(argc,argv);
    std::cout<<"\n=== Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<"  E="<<hasE<<" Et="<<hasEt
             <<"  constraint="<<gConstraint
             <<"  input="<<gBsaFile<<" ===\n";

    LoadData();
    BinBsaData();
    std::cout<<" BSA bins="<<Nbins<<" (raw "<<bsaData.size()<<")\n\n";

    // ─── Stage 1: Im‐fit ─────────────────────────────────────────────────────────
    gStage=1;
    build_par_list();
    int nim=parNamesIm.size();
    std::vector<double> imVal(nim), imErr(nim);
    double chi2_im,edm,errdef; int nv,nx,ic,ndf_im;

    {
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);

        for(int i=0;i<nim;++i){
            const auto &nm=parNamesIm[i];
            double init=0, step=0.01;
            #define GETINIT(NAME) if(nm==#NAME) init=NAME;
            GETINIT(renormImag)
            GETINIT(r_H)        GETINIT(alpha0_H)
            GETINIT(alpha1_H)   GETINIT(n_H)
            GETINIT(b_H)        GETINIT(M2_H)
            GETINIT(P_H)
            GETINIT(r_Ht)       GETINIT(alpha0_Ht)
            GETINIT(alpha1_Ht)  GETINIT(n_Ht)
            GETINIT(b_Ht)       GETINIT(M2_Ht)
            GETINIT(P_Ht)
            GETINIT(r_E)        GETINIT(alpha0_E)
            GETINIT(alpha1_E)   GETINIT(n_E)
            GETINIT(b_E)        GETINIT(M2_E)
            GETINIT(P_E)
            GETINIT(r_Et)       GETINIT(alpha0_Et)
            GETINIT(alpha1_Et)  GETINIT(n_Et)
            GETINIT(b_Et)       GETINIT(M2_Et)
            GETINIT(P_Et)
            #undef GETINIT

            double lo=-1e3, hi=1e3;
            if(nm.rfind("M2_",0)==0 || nm.rfind("r_",0)==0) lo=0.0;
            minu.DefineParameter(i,nm.c_str(),init,step,lo,hi);

            // fix the ones we never float
            if(nm=="renormImag"
            || nm.rfind("alpha0_",0)==0
            || nm.rfind("alpha1_",0)==0
            || nm.rfind("n_",0)==0
            || nm.rfind("P_",0)==0)
              minu.FixParameter(i);
        }

        std::cout<<"Stage1: fitting r_*, b_*, M2_* (all others fixed)...\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_im,edm,errdef,nv,nx,ic);
        for(int i=0;i<nim;++i)
            minu.GetParameter(i,imVal[i],imErr[i]);
        ndf_im = Nbins - nim;
    }

    // collect into maps
    std::map<std::string,double> valMap, errMap;
    for(int i=0;i<nim;++i){
        valMap[parNamesIm[i]] = imVal[i];
        errMap[parNamesIm[i]] = imErr[i];
    }

    // ─── Stage 2: renormReal-fit ────────────────────────────────────────────────
    if(gStrategy==2){
        gStage=2;
        double chi2_re,edm2,errdef2; int nv2,nx2,ic2,ndf_re;
        double reVal,reErr;
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
        valMap["renormReal"]=reVal;
        errMap["renormReal"]=reErr;
        chi2_im = chi2_re;
        ndf_im  = ndf_re;
    }

    // prepare output names
    std::vector<std::string> outNames = parNamesIm;
    if(gStrategy==2) outNames.push_back("renormReal");

    // write results file
    time_t now=time(nullptr);
    char tb[32];
    strftime(tb,sizeof(tb),"%Y%m%d_%H%M%S",localtime(&now));
    std::string fname="output/fit_results/fit_results_"+std::string(tb)+".txt";
    std::ofstream fout(fname);

    fout<<"# fit_CFFs results\n";
    fout<<"timestamp   "<<tb<<"\n";
    fout<<"strategy    "<<gStrategy<<"\n";
    fout<<"constraint  "<<gConstraint<<"\n";
    fout<<"input       "<<gBsaFile<<"\n";
    fout<<"H "<<hasH<<"  Ht "<<hasHt<<"  E "<<hasE<<"  Et "<<hasEt<<"\n";
    fout<<"# parameters:";
    for(auto &n: outNames) fout<<" "<<n;
    fout<<"\n# values:\n";
    for(auto &n: outNames) fout<<valMap[n]<<" ";
    fout<<"\n# errors:\n";
    for(auto &n: outNames) fout<<errMap[n]<<" ";
    fout<<"\n# chi2 ndf chi2/ndf\n";
    fout<<chi2_im<<" "<<ndf_im<<" "<<(chi2_im/ndf_im)<<"\n";
    fout.close();

    std::cout<<"Wrote fit results to "<<fname<<"\n\n";

    // echo to stdout
    std::cout<<"--- Fit Results ---\n";
    for(auto &n: outNames){
      std::cout<<" "<<n<<" = "<<valMap[n]<<" ± "<<errMap[n]<<"\n";
    }
    std::cout<<" χ²/ndf = "<<chi2_im<<"/"<<ndf_im
             <<" = "<<(chi2_im/ndf_im)<<"\n";

    return 0;
}