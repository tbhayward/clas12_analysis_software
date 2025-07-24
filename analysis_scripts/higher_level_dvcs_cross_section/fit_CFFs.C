// fit_CFFs.C
// ──────────
// program to fit arbitrary subsets of the DVCS Compton Form Factors (CFFs)
// under four possible strategies:
//   1) fit only Im‐parts -> BSA data
//   2) fit only Re‐parts (new real‐ansatz parameters + renormReal) -> xsec data
//   3) fit Im‐parts + Re‐parts -> BSA + xsec simultaneously
//   4) two‐step: (a) fit Im‐parts -> BSA, (b) fit Re‐parts -> xsec
//
// Usage:
//   ./fit_CFFs -strategy <1|2|3|4> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>
//
// e.g.
//   ./fit_CFFs -strategy 3 -H 1 -Ht 1 -E 0 -Et 0
//
// Compile:
//   g++ -O2 fit_CFFs.C `root-config --cflags --libs` -lMinuit -o fit_CFFs
// Run as above.

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

// all of the CFF shape‐parameter globals (Im parts)
extern double r_H,      alpha0_H,  alpha1_H,  n_H,   b_H,   M2_H,  P_H;
extern double r_Ht,     alpha0_Ht, alpha1_Ht, n_Ht,  b_Ht,  M2_Ht, P_Ht;
extern double r_E,      alpha0_E,  alpha1_E,  n_E,   b_E,   M2_E,  P_E;
extern double r_Et,     alpha0_Et, alpha1_Et, n_Et,  b_Et,  M2_Et, P_Et;

// real‐ansatz globals (Re parts)
extern double A_H, B_H, C_H, D_H, E_H;
extern double A_Ht, B_Ht;
extern double A_E, B_E, C_E, D_E;
extern double A_Et, B_Et;

// ──────────────────────────────────────────────────────────────────────────────
// globals controlling strategy & stage:
static int  gStrategy = 0;   // 1, 2, 3 or 4
static int  gStage    = 1;   // for Strategy 4: stage=1 (Im->BSA) or 2 (Re->xsec)

// Data structures & loading
struct DataPoint {
    double phi, Q2, xB, t, Eb, A, sigA;
};
static std::vector<DataPoint> bsaData, xsData;

void LoadData(){
    auto readFile = [&](const char* fname,
                        std::vector<DataPoint>& vec)
    {
        std::ifstream in(fname);
        if(!in){ std::cerr<<"ERROR: cannot open "<<fname<<"\n"; std::exit(1);}        
        std::string line;
        while(std::getline(in,line)){
            if(line.empty()||line[0]=='#') continue;
            std::istringstream iss(line);
            DataPoint d;
            iss>>d.phi>>d.Q2>>d.xB>>d.t
               >>d.Eb>>d.A>>d.sigA;
            vec.push_back(d);
        }
    };
    readFile("imports/rga_prl_bsa.txt",       bsaData);
    readFile("imports/rga_pass1_xsec_2018.txt", xsData);
}

// ──────────────────────────────────────────────────────────────────────────────
// dynamic parameter list
static std::vector<std::string> parNames;
static int idx_H=-1, idx_Ht=-1, idx_E=-1, idx_Et=-1, idx_R=-1;

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
                     <<" -strategy <1|2|3|4> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>4){
        std::cerr<<"Invalid strategy "<<gStrategy<<"\n";
        std::exit(1);
    }
}

// build_par_list(): construct the list of Minuit parameters based on flags & strategy
void build_par_list(){
    parNames.clear();
    idx_H = idx_Ht = idx_E = idx_Et = idx_R = -1;

    // imaginary‐fit stage: BSA
    bool doImag = (gStrategy==1) || (gStrategy==3) || (gStrategy==4 && gStage==1);
    if(doImag){
        parNames.push_back("renormImag");
        if(hasH){ idx_H = parNames.size();
            parNames.insert(parNames.end(),
                {"r_H","alpha0_H","alpha1_H","n_H","b_H","M2_H","P_H"});
        }
        if(hasHt){ idx_Ht = parNames.size();
            parNames.insert(parNames.end(),
                {"r_Ht","alpha0_Ht","alpha1_Ht","n_Ht","b_Ht","M2_Ht","P_Ht"});
        }
        if(hasE){ idx_E = parNames.size();
            parNames.insert(parNames.end(),
                {"r_E","alpha0_E","alpha1_E","n_E","b_E","M2_E","P_E"});
        }
        if(hasEt){ idx_Et = parNames.size();
            parNames.insert(parNames.end(),
                {"r_Et","alpha0_Et","alpha1_Et","n_Et","b_Et","M2_Et","P_Et"});
        }
    }

    // real‐fit stage: xsec
    bool doReal = (gStrategy==2) || (gStrategy==3) || (gStrategy==4 && gStage==2);
    if(doReal){
        if(hasH){
            parNames.insert(parNames.end(),
                {"A_H","B_H","C_H","D_H","E_H"});
        }
        if(hasHt){
            parNames.insert(parNames.end(), {"A_Ht","B_Ht"});
        }
        if(hasE){
            parNames.insert(parNames.end(),
                {"A_E","B_E","C_E","D_E"});
        }
        if(hasEt){
            parNames.insert(parNames.end(), {"A_Et","B_Et"});
        }
        idx_R = parNames.size();
        parNames.push_back("renormReal");
    }
}

// fcn(): Minuit’s chi2.  Unpack par[] into the global CFF parameters.
void fcn(int & /*npar*/, double* /*grad*/, double &f, double *par, int /*iflag*/){
    int ip = 0;
    // imaginary stage
    bool doImag = (gStrategy==1) || (gStrategy==3) || (gStrategy==4 && gStage==1);
    if(doImag){
        renormImag = par[ip++];
        if(hasH){ r_H=par[ip++]; alpha0_H=par[ip++]; alpha1_H=par[ip++];
                  n_H=par[ip++]; b_H=par[ip++]; M2_H=par[ip++]; P_H=par[ip++]; }
        if(hasHt){ r_Ht=par[ip++]; alpha0_Ht=par[ip++]; alpha1_Ht=par[ip++];
                   n_Ht=par[ip++]; b_Ht=par[ip++]; M2_Ht=par[ip++]; P_Ht=par[ip++]; }
        if(hasE){ r_E=par[ip++]; alpha0_E=par[ip++]; alpha1_E=par[ip++];
                  n_E=par[ip++]; b_E=par[ip++]; M2_E=par[ip++]; P_E=par[ip++]; }
        if(hasEt){ r_Et=par[ip++]; alpha0_Et=par[ip++]; alpha1_Et=par[ip++];
                   n_Et=par[ip++]; b_Et=par[ip++]; M2_Et=par[ip++]; P_Et=par[ip++]; }
    }

    // real stage
    bool doReal = (gStrategy==2) || (gStrategy==3) || (gStrategy==4 && gStage==2);
    if(doReal){
        if(hasH){ A_H=par[ip++]; B_H=par[ip++]; C_H=par[ip++]; D_H=par[ip++]; E_H=par[ip++]; }
        if(hasHt){ A_Ht=par[ip++]; B_Ht=par[ip++]; }
        if(hasE){ A_E=par[ip++]; B_E=par[ip++]; C_E=par[ip++]; D_E=par[ip++]; }
        if(hasEt){ A_Et=par[ip++]; B_Et=par[ip++]; }
        renormReal = par[ip++];
    }

    // compute chi2
    double chi2 = 0;
    bool doBSA = (gStrategy==1) || (gStrategy==3) || (gStrategy==4 && gStage==1);
    if(doBSA){
        for(auto &d : bsaData){
            BMK_DVCS dvcs(-1,1,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double mA = dvcs.BSA(); double r = (d.A - mA)/d.sigA;
            chi2 += r*r;
        }
    }
    bool doXS = (gStrategy==2) || (gStrategy==3) || (gStrategy==4 && gStage==2);
    if(doXS){
        for(auto &d : xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double mA = dvcs.CrossSection(); double r = (d.A - mA)/d.sigA;
            chi2 += r*r;
        }
    }
    f = chi2;
}

// main(): orchestrate parsing, fitting, and I/O
int main(int argc, char** argv){
    parse_args(argc, argv);
    std::cout<<"\n=== fit_CFFs Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<" E="<<hasE<<" Et="<<hasEt<<" ===\n";
    LoadData();
    std::cout<<"  loaded "<<bsaData.size()
             <<" BSA and "<<xsData.size()
             <<" xsec points\n\n";

    std::vector<double> val, err;
    double finalChi2=0; int finalNdf=0;

    // Strategy 2: fit real parts + renormReal -> xsec only
    if(gStrategy==2){
        gStage = 2;
        build_par_list();
        int npar = parNames.size();
        val.assign(npar,0); err.assign(npar,0);
        TMinuit minuit(npar); minuit.SetPrintLevel(1); minuit.SetFCN(fcn);
        for(int i=0;i<npar;++i){
            auto &nm = parNames[i];
            double init=1.0, step=0.01;
            if(nm=="A_H") init = A_H;
            else if(nm=="B_H") init = B_H;
            else if(nm=="C_H") init = C_H;
            else if(nm=="D_H") init = D_H;
            else if(nm=="E_H") init = E_H;
            else if(nm=="A_Ht") init = A_Ht;
            else if(nm=="B_Ht") init = B_Ht;
            else if(nm=="A_E") init = A_E;
            else if(nm=="B_E") init = B_E;
            else if(nm=="C_E") init = C_E;
            else if(nm=="D_E") init = D_E;
            else if(nm=="A_Et") init = A_Et;
            else if(nm=="B_Et") init = B_Et;
            else if(nm=="renormReal") init = 1.0;
            // six-arg DefineParameter: [-1e6, +1e6] bounds
            minuit.DefineParameter(i, nm.c_str(), init, step, -1e6, 1e6);
        }
        std::cout<<"Strategy 2: fitting Re parts + renormReal->xsec only...\n";
        minuit.Migrad();
        for(int i=0;i<npar;++i) minuit.GetParameter(i,val[i],err[i]);
        double edm,errdef; int nv,nx,ic; double chi2;
        minuit.mnstat(chi2,edm,errdef,nv,nx,ic);
        finalChi2 = chi2; finalNdf = int(xsData.size()) - npar;
        std::cout<<"\n=== Results ===\n";
        for(int i=0;i<npar;++i)
            std::cout<<" "<<parNames[i]<<" = "<<val[i]<<" ± "<<err[i]<<"\n";
        std::cout<<" χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                 <<" = "<<(finalChi2/finalNdf)<<"\n\n";
    }
    // Strategies 1 & 3
    else if(gStrategy==1 || gStrategy==3){
        bool sim = (gStrategy==3);
        std::cout<<"Strategy "<<(sim?3:1)
                 <<": "<<(sim?"fitting Im+Re -> BSA+xsec simultaneously...":"fitting Im->BSA only...")<<"\n";
        gStage = 1;
        build_par_list(); int npar = parNames.size();
        val.assign(npar,0); err.assign(npar,0);
        TMinuit minuit(npar); minuit.SetPrintLevel(1); minuit.SetFCN(fcn);
        for(int i=0;i<npar;++i){
            double init=1.0, step=0.01;
            auto &name = parNames[i];
            if(name.find("alpha0_")==0)     init=0.43;
            if(name.find("alpha1_")==0)     init=0.85;
            if(name.find("n_")==0)          init=1.35;
            if(name.find("b_")==0)          init=0.4;
            if(name.find("M2_")==0)         init=0.64;
            if(name.find("P_")==0)          init=1.0;
            if(name=="r_Ht")                init=7.0;
            if(name=="renormReal")          init=1.0;
            minuit.DefineParameter(i, name.c_str(), init, step, -1e6, 1e6);
        }
        minuit.Migrad();
        for(int i=0;i<npar;++i) minuit.GetParameter(i,val[i],err[i]);
        double edm,errdef; int nv,nx,ic;
        minuit.mnstat(finalChi2,edm,errdef,nv,nx,ic);
        finalNdf = (sim
                   ? int(bsaData.size()+xsData.size()) - npar
                   : int(bsaData.size()) - npar);
        std::cout<<"\n=== Results ===\n";
        for(int i=0;i<npar;++i)
            std::cout<<" "<<parNames[i]<<" = "<<val[i]<<" ± "<<err[i]<<"\n";
        std::cout<<" χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                 <<" = "<<(finalChi2/finalNdf)<<"\n\n";
    }
    // Strategy 4: two-step
    else {
        std::cout<<"Strategy 4: two-step fit: (1) Im->BSA, (2) Re->xsec...\n";
        // stage 1
        gStage = 1; build_par_list(); if(idx_R>=0) parNames.resize(idx_R);
        int n1 = parNames.size(); val.assign(n1,0); err.assign(n1,0);
        {
            TMinuit m1(n1); m1.SetPrintLevel(1); m1.SetFCN(fcn);
            for(int i=0;i<n1;++i){
                auto &nm = parNames[i];
                double init=1.0, step=0.01;
                if(nm.find("alpha0_")==0) init=0.43;
                if(nm.find("alpha1_")==0) init=0.85;
                if(nm.find("n_")==0)     init=1.35;
                if(nm.find("b_")==0)     init=0.4;
                if(nm.find("M2_")==0)    init=0.64;
                if(nm.find("P_")==0)     init=1.0;
                if(nm=="r_Ht")           init=7.0;
                m1.DefineParameter(i, nm.c_str(), init, step, -1e6, 1e6);
            }
            std::cout<<"Stage 1: fitting Im->BSA...\n";
            m1.Migrad(); double edm,errdef; int nv,nx,ic; double chi2;
            m1.mnstat(chi2,edm,errdef,nv,nx,ic);
            finalChi2=chi2; finalNdf=int(bsaData.size())-n1;
            std::cout<<"\nStage 1 χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                     <<" = "<<(finalChi2/finalNdf)<<"\n\n";
        }
        // stage 2
        gStage = 2;
        build_par_list();
        int n2 = parNames.size(); val.assign(n2,0); err.assign(n2,0);
        {
            TMinuit m2(n2); m2.SetPrintLevel(1); m2.SetFCN(fcn);
            for(int i=0;i<n2;++i){
                auto &nm = parNames[i];
                double init=1.0, step=0.01;
                if(nm=="A_H") init=A_H;
                else if(nm=="B_H") init=B_H;
                else if(nm=="C_H") init=C_H;
                else if(nm=="D_H") init=D_H;
                else if(nm=="E_H") init=E_H;
                else if(nm=="A_Ht") init=A_Ht;
                else if(nm=="B_Ht") init=B_Ht;
                else if(nm=="A_E") init=A_E;
                else if(nm=="B_E") init=B_E;
                else if(nm=="C_E") init=C_E;
                else if(nm=="D_E") init=D_E;
                else if(nm=="A_Et") init=A_Et;
                else if(nm=="B_Et") init=B_Et;
                else if(nm=="renormReal") init=1.0;
                m2.DefineParameter(i, nm.c_str(), init, step, -1e6, 1e6);
            }
            std::cout<<"Stage 2: fitting Re parts + renormReal->xsec...\n";
            m2.Migrad();
            for(int i=0;i<n2;++i) m2.GetParameter(i,val[i],err[i]);
            double edm,errdef; int nv,nx,ic; double chi2;
            m2.mnstat(chi2,edm,errdef,nv,nx,ic);
            finalChi2=chi2; finalNdf=int(xsData.size())-n2;
            std::cout<<"\nStage 2 Results:\n";
            for(int i=0;i<n2;++i)
                std::cout<<" "<<parNames[i]<<" = "<<val[i]<<" ± "<<err[i]<<"\n";
            std::cout<<" χ²/ndf = "<<finalChi2<<"/"<<finalNdf
                     <<" = "<<(finalChi2/finalNdf)<<"\n\n";
        }
    }

    // write results file
    system("mkdir -p output/fit_results/");
    time_t now = time(nullptr);
    char tb[32];
    strftime(tb,sizeof(tb),"%Y%m%d_%H%M%S",localtime(&now));
    std::string fname = std::string("output/fit_results/fit_results_") + tb + ".txt";
    std::ofstream out(fname);
    out<<"# fit_CFFs results\n";
    out<<"timestamp   "<<tb<<"\n";
    out<<"strategy    "<<gStrategy<<"\n";
    out<<"H "<<hasH<<"  Ht "<<hasHt<<"  E "<<hasE<<"  Et "<<hasEt<<"\n";
    out<<"# parameters:";
    for(auto &n:parNames) out<<" "<<n;
    out<<"\n# values:\n";
    for(double v:val) out<<v<<" "; out<<"\n";
    out<<"# errors:\n";
    for(double e:err) out<<e<<" "; out<<"\n";
    out<<"# chi2 ndf chi2/ndf\n";
    out<<finalChi2<<" "<<finalNdf<<" "<<(finalChi2/finalNdf)<<"\n";
    out.close();

    std::cout<<"Wrote fit results to "<<fname<<"\n";
    return 0;
}