// fit_CFFs.C
// ──────────
// program to fit arbitrary subsets of the DVCS Compton Form Factors (CFFs)
// under three possible strategies:
//   1) fit only Im‐parts → BSA data
//   2) fit Im‐parts + renormReal → BSA + xsec simultaneously
//   3) two‐step: (a) fit Im‐parts→BSA, (b) fit renormReal→xsec
//
// Usage:
//   ./fit_CFFs -strategy <1|2|3> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>
//
// e.g.
//   ./fit_CFFs -strategy 2 -H 1 -Ht 1 -E 0 -Et 0
//
// Compile:
//   g++ -O2 fit_CFFs.C `root-config --cflags --libs` -lMinuit -o fit_CFFs
// Run as above.

#include <cmath>
#include <iostream>
#include <vector>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"

// pull in your full BMK_DVCS + CFF code, with globals
#include "DVCS_xsec.C"

// ──────────────────────────────────────────────────────────────────────────────
// globals controlling which CFFs to include & strategy/stage
static int  gStrategy = 0;   // 1, 2 or 3
static int  gStage    = 1;   // for Strategy 3: stage=1 (Im→BSA) or 2 (Re→xsec)
bool hasH  = false;          // fit H?
bool hasHt = false;          // fit Ĥ?
bool hasE  = false;          // fit E?
bool hasEt = false;          // fit Ė?

// we will dynamically build our parameter list here:
static std::vector<std::string> parNames;
static int idx_H   = -1,
           idx_Ht  = -1,
           idx_E   = -1,
           idx_Et  = -1,
           idx_R   = -1;

// -----------------------------------------------------------------------------
// parse_args(): read -strategy, -H, -Ht, -E, -Et
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
    while((c = getopt_long(argc,argv,"s:h:t:e:x:",long_opts,nullptr))!=-1){
        switch(c){
          case 's': gStrategy = std::atoi(optarg); break;
          case 'h': hasH       = std::atoi(optarg); break;
          case 't': hasHt      = std::atoi(optarg); break;
          case 'e': hasE       = std::atoi(optarg); break;
          case 'x': hasEt      = std::atoi(optarg); break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" -strategy <1|2|3> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>3){
        std::cerr<<"Invalid strategy "<<gStrategy<<"\n"; std::exit(1);
    }
}

// -----------------------------------------------------------------------------
// Build the list of parameters (names & initial values) in the order we will
// hand them to Minuit.  We record the starting index of each block.
// -----------------------------------------------------------------------------
void build_par_list(){
    parNames.clear();
    // 0: always fit a global renormImag
    parNames.push_back("renormImag");
    // 1–7: H ansatz
    if(hasH){
        idx_H = parNames.size();
        parNames.push_back("r_H");
        parNames.push_back("alpha0_H");
        parNames.push_back("alpha1_H");
        parNames.push_back("n_H");
        parNames.push_back("b_H");
        parNames.push_back("Mm2_H");
        parNames.push_back("P_H");
    }
    // next 7: Ht ansatz
    if(hasHt){
        idx_Ht = parNames.size();
        parNames.push_back("r_Ht");
        parNames.push_back("alpha0_Ht");
        parNames.push_back("alpha1_Ht");
        parNames.push_back("n_Ht");
        parNames.push_back("b_Ht");
        parNames.push_back("Mm2_Ht");
        parNames.push_back("P_Ht");
    }
    // next 7: E ansatz
    if(hasE){
        idx_E = parNames.size();
        parNames.push_back("r_E");
        parNames.push_back("alpha0_E");
        parNames.push_back("alpha1_E");
        parNames.push_back("n_E");
        parNames.push_back("b_E");
        parNames.push_back("Mm2_E");
        parNames.push_back("P_E");
    }
    // Et currently no shape parameters, so we just toggle it on/off
    // finally: renormReal if doing a combined BSA+xsec fit or stage2
    if((gStrategy==2) || (gStrategy==3 && gStage==2)){
        idx_R = parNames.size();
        parNames.push_back("renormReal");
    }
}

// -----------------------------------------------------------------------------
// fcn(): Minuit’s χ² function.  Unpack *only* those parameters that are turned on.
// -----------------------------------------------------------------------------
void fcn(int &npar, double *grad, double &f, double *par, int iflag){
    int ip=0;
    // 1) renormImag always
    renormImag = par[ip++];
    // 2) H block
    if(hasH){
        r_H      = par[ip++];
        alpha0_H = par[ip++];
        alpha1_H = par[ip++];
        n_H      = par[ip++];
        b_H      = par[ip++];
        Mm2_H    = par[ip++];
        P_H      = par[ip++];
    }
    // 3) Ht block
    if(hasHt){
        r_Ht      = par[ip++];
        alpha0_Ht = par[ip++];
        alpha1_Ht = par[ip++];
        n_Ht      = par[ip++];
        b_Ht      = par[ip++];
        Mm2_Ht    = par[ip++];
        P_Ht      = par[ip++];
    }
    // 4) E block
    if(hasE){
        r_E      = par[ip++];
        alpha0_E = par[ip++];
        alpha1_E = par[ip++];
        n_E      = par[ip++];
        b_E      = par[ip++];
        Mm2_E    = par[ip++];
        P_E      = par[ip++];
    }
    // 5) renormReal for combined fits
    if(idx_R>=0){
        renormReal = par[ip++];
    }

    // switch on only the CFFs we want
    hasH  = hasH;
    hasHt = hasHt;
    hasE  = hasE;
    hasEt = hasEt;

    double chi2 = 0;

    // 6) BSA?
    bool doBSA = (gStrategy==1||gStrategy==2|| (gStrategy==3&&gStage==1));
    if(doBSA){
        for(auto &d : bsaData){
            BMK_DVCS dvcs(-1,1,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double mA = dvcs.BSA();
            double r  = (d.A - mA)/d.sigA;
            chi2 += r*r;
        }
    }

    // 7) xsec?
    bool doXS = (gStrategy==2|| (gStrategy==3&&gStage==2));
    if(doXS){
        for(auto &d : xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double mA = dvcs.CrossSection();
            double r  = (d.A - mA)/d.sigA;
            chi2 += r*r;
        }
    }

    f = chi2;
}

// -----------------------------------------------------------------------------
// main(): parse args, build parameter list, run Minuit fits, write results
// -----------------------------------------------------------------------------
int main(int argc, char** argv){
    parse_args(argc, argv);

    std::cout<<"\n=== fit_CFFs Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<" E="<<hasE<<" Et="<<hasEt<<" ===\n";

    LoadData();
    std::cout<<"  loaded "<<bsaData.size()<<" BSA and "
             <<xsData.size()<<" xsec points\n\n";

    // first stage, build parNames based on flags & gStage
    build_par_list();
    int npar = parNames.size();

    // storage for results
    std::vector<double> val(npar), err(npar);
    double finalChi2=0; int finalNdf=0;

    // ─── Strategies 1 & 2 ──────────────────────────────────────────────
    if(gStrategy==1 || gStrategy==2){
        TMinuit minuit(npar);
        minuit.SetPrintLevel(1);
        minuit.SetFCN(fcn);

        // define all parameters with reasonable defaults & steps
        for(int i=0;i<npar;++i){
            double init = 1.0, step=0.1, mn=0, mx=10;
            // customize some defaults:
            auto &name = parNames[i];
            if(name=="alpha0_H"||name=="alpha0_Ht"||name=="alpha0_E") init=0.43, step=0.05, mn=-5, mx=5;
            if(name=="alpha1_H"||name=="alpha1_Ht"||name=="alpha1_E") init=0.85, step=0.05, mn=-10, mx=10;
            if(name.find("n_")==0)  init=1.35, step=0.05;
            if(name.find("b_")==0)  init=0.4,  step=0.05;
            if(name.find("Mm2_")==0){init=0.64; step=0.05;}
            if(name.find("P_")==0)  init=1.0,  step=0.05;
            if(name=="r_Ht")        init=7.0,  step=1.0;
            if(name=="r_H" || name=="r_E") init=0.9, step=0.1;
            if(name=="renormReal")  init=1.0,  step=0.1;

            minuit.DefineParameter(i, name.c_str(), init, step, mn, mx);
        }

        std::cout<<"Running Migrad()…\n";
        minuit.Migrad();

        // retrieve
        for(int i=0;i<npar;++i)
            minuit.GetParameter(i,val[i],err[i]);

        // χ² stats
        double edm, errdef; int nv, nx, ic;
        minuit.mnstat(finalChi2,edm,errdef,nv,nx,ic);
        finalNdf = (gStrategy==1
                   ? int(bsaData.size()) - npar
                   : int(bsaData.size()+xsData.size()) - npar);

        // print to console
        std::cout<<"\n=== Results ===\n";
        for(int i=0;i<npar;++i){
            std::cout<<" "<<parNames[i]<<" = "<<val[i]<<" ± "<<err[i]<<"\n";
        }
        std::cout<<" χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                 <<" = "<<(finalChi2/finalNdf)<<"\n\n";
    }

    // ─── Strategy 3: two‐step ────────────────────────────────────────────
    else {
        // -- Stage 1: fit all Im‐blocks exactly as above but WITHOUT renormReal
        gStage = 1;
        build_par_list();
        // strip out renormReal if it snuck in
        if(idx_R>=0){
            parNames.resize(idx_R);
            idx_R = -1;
        }
        npar = parNames.size();
        {
            TMinuit minuit(npar);
            minuit.SetPrintLevel(1);
            minuit.SetFCN(fcn);
            for(int i=0;i<npar;++i){
                // same DefineParameter loop as before...
                double init=1,step=0.1,mn=0,mx=10;
                auto &name = parNames[i];
                if(name=="alpha0_H"||name=="alpha0_Ht"||name=="alpha0_E") init=0.43,step=0.05,mn=-5,mx=5;
                if(name=="alpha1_H"||name=="alpha1_Ht"||name=="alpha1_E") init=0.85,step=0.05,mn=-10,mx=10;
                if(name.find("n_")==0)  init=1.35,step=0.05;
                if(name.find("b_")==0)  init=0.4, step=0.05;
                if(name.find("Mm2_")==0){init=0.64;step=0.05;}
                if(name.find("P_")==0)  init=1.0, step=0.05;
                if(name=="r_Ht")        init=7.0, step=1.0;
                if(name=="r_H"||name=="r_E") init=0.9, step=0.1;
                minuit.DefineParameter(i,name.c_str(),init,step,mn,mx);
            }
            std::cout<<"Stage 1: fitting Im→BSA…\n";
            minuit.Migrad();

            for(int i=0;i<npar;++i)
                minuit.GetParameter(i,val[i],err[i]);
            double chi2,edm,errdef; int nv,nx,ic;
            minuit.mnstat(chi2,edm,errdef,nv,nx,ic);
            finalChi2 = chi2;
            finalNdf  = int(bsaData.size()) - npar;
            std::cout<<"\nStage 1 χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                     <<" = "<<(finalChi2/finalNdf)<<"\n\n";
        }

        // -- Stage 2: only renormReal → xsec
        gStage = 2;
        parNames.clear();
        parNames.push_back("renormReal");
        npar = 1;
        {
            TMinuit minuit(npar);
            minuit.SetPrintLevel(1);
            minuit.SetFCN(fcn);
            minuit.DefineParameter(0,"renormReal",1.0,0.1,0,10);
            std::cout<<"Stage 2: fitting renormReal→xsec…\n";
            minuit.Migrad();
            minuit.GetParameter(0,val[0],err[0]);
            double chi2,edm,errdef; int nv,nx,ic;
            minuit.mnstat(chi2,edm,errdef,nv,nx,ic);
            finalChi2 = chi2;
            finalNdf  = int(xsData.size()) - 1;
            std::cout<<"\nStage 2 χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                     <<" = "<<(finalChi2/finalNdf)<<"\n\n";
        }
    }

    // ─── Write results file ────────────────────────────────────────────────
    system("mkdir -p output/fit_results/");
    time_t now = time(nullptr);
    char tb[32];
    strftime(tb,sizeof(tb),"%Y%m%d_%H%M%S",localtime(&now));
    std::string fname = std::string("output/fit_results/fit_results_")+tb+".txt";
    std::ofstream out(fname);
    out<<"# fit_CFFs results\n";
    out<<"timestamp   "<<tb<<"\n";
    out<<"strategy    "<<gStrategy<<"\n";
    out<<"H  "<<hasH<<"  Ht "<<hasHt<<"  E "<<hasE<<"  Et "<<hasEt<<"\n";
    out<<"# parameters: ";
    for(auto &n:parNames) out<<n<<" ";
    out<<"\n# values:\n";
    for(double v: val) out<<v<<" "; out<<"\n";
    out<<"# errors:\n";
    for(double e: err) out<<e<<" "; out<<"\n";
    out<<"# chi2 ndf chi2/ndf\n";
    out<<finalChi2<<" "<<finalNdf<<" "<<(finalChi2/finalNdf)<<"\n";
    out.close();

    std::cout<<"Wrote fit results to "<<fname<<"\n";
    return 0;
}