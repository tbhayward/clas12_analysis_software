// FitH.C
// ──────
// Stand‐alone program to fit the GPD‐H ansatz under three possible strategies:
//   1) fit only ImH parameters to BSA data
//   2) fit ImH parameters + renormReal to BSA + xsec simultaneously
//   3) two‐step: (a) fit ImH→BSA, (b) fit renormReal→xsec
//
// Usage:  ./FitH <strategy>    where <strategy> ∈ {1,2,3}
//
// Compile:
//   g++ -O2 FitH.C `root-config --cflags --libs` -lMinuit -o FitH
// Run, e.g.:
//   ./FitH 3

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <ctime>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"

// pull in BMK_DVCS implementation (with globals
//   renormImag, alpha0, alpha1, n_val, b_val, Mm2_val, P_val, renormReal
//   and the functions GetImH, GetReH, CrossSection(), BSA(), etc.)
#include "DVCS_xsec.C"

// ──────────────────────────────────────────────────────────────────────────────
// globals controlling strategy & stage:
static int gStrategy = 0;   // 1, 2 or 3
static int gStage    = 1;   // for Strategy 3: stage=1 (ImH) or 2 (ReH)

// struct for each measurement
struct DataPoint {
    double phi, Q2, xB, t, Eb, A, sigA;
};

// two data sets
static std::vector<DataPoint> bsaData, xsData;

// -----------------------------------------------------------------------------
// LoadData(): read imports/*.txt into bsaData & xsData
void LoadData(){
    auto readFile = [&](const char* fname,
                        std::vector<DataPoint>& vec)
    {
        std::ifstream in(fname);
        if(!in){
            std::cerr<<"ERROR: cannot open "<<fname<<"\n"; std::exit(1);
        }
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

// -----------------------------------------------------------------------------
// fcn(): Minuit’s χ² function.  Depending on gStrategy/gStage it
//         • unpacks par[] → globals
//         • loops over BSA and/or xsec
void fcn(int& /*npar*/, double* /*grad*/, double& f, double* par, int /*iflag*/){
    // 1) unpack parameters
    if(gStrategy==1 || gStrategy==2
    || (gStrategy==3 && gStage==1)){
        renormImag = par[0];
        alpha0     = par[1];
        alpha1     = par[2];
        n_val      = par[3];
        b_val      = par[4];
        Mm2_val    = par[5];
        P_val      = par[6];
    }
    if(gStrategy==2){
        renormReal = par[7];
    } else if(gStrategy==3 && gStage==2){
        renormReal = par[0];
    }

    // only H‐term
    hasH = true;
    hasHt = hasE = hasEt = false;

    double chi2 = 0;

    // 2) BSA?
    bool doBSA = ( gStrategy==1 || gStrategy==2
                || (gStrategy==3 && gStage==1) );
    if(doBSA){
        for(auto &d : bsaData){
            BMK_DVCS dvcs(-1,1,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double mA = dvcs.BSA();
            double r  = (d.A - mA)/d.sigA;
            chi2 += r*r;
        }
    }

    // 3) xsec?
    bool doXS = ( gStrategy==2
               || (gStrategy==3 && gStage==2) );
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
// main(): parse strategy, run fits, save results + χ²/ndf to output/fit_results_*.txt
int main(int argc, char** argv){
    if(argc<2){
        std::cerr<<"Usage: "<<argv[0]<<" <strategy 1|2|3>\n";
        return 1;
    }
    gStrategy = std::atoi(argv[1]);
    if(gStrategy<1||gStrategy>3){
        std::cerr<<"Invalid strategy "<<gStrategy<<"\n"; return 1;
    }

    std::cout<<"\n=== FitH Strategy="<<gStrategy<<" ===\n";
    LoadData();
    std::cout<<"  loaded "<<bsaData.size()
             <<" BSA and "<<xsData.size()<<" xsec points\n\n";

    // storage for final results
    double valIm[7], errIm[7];
    double valReal=renormReal, errReal=0;
    double finalChi2=0;
    int    finalNdf=0;

    // ────────────────────────────────
    // Strategy 1 or 2: single Minuit run
    if(gStrategy==1 || gStrategy==2){
        int npar = (gStrategy==1?7:8);
        TMinuit minuit(npar);
        minuit.SetPrintLevel(1);
        minuit.SetFCN(fcn);

        // ImH params 0–6
        minuit.DefineParameter(0,"renormImag",1.0,0.1,0,10);
        minuit.DefineParameter(1,"alpha0",    0.43,0.05,-5,5);
        minuit.DefineParameter(2,"alpha1",    0.85,0.05,-10,10);
        minuit.DefineParameter(3,"n_val",     1.35,0.05,0,10);
        minuit.DefineParameter(4,"b_val",     0.40,0.05,0,10);
        minuit.DefineParameter(5,"Mm2_val",   0.64,0.05,0,10);
        minuit.DefineParameter(6,"P_val",     1.00,0.05,0,10);
        if(gStrategy==2){
            minuit.DefineParameter(7,"renormReal",1.0,0.1,0,10);
        }

        std::cout<<"Running Migrad()…\n";
        minuit.Migrad();

        // retrieve
        double v[8], e[8];
        for(int i=0;i<npar;++i) minuit.GetParameter(i,v[i],e[i]);
        valReal = (gStrategy==2? v[7] : renormReal);
        errReal = (gStrategy==2? e[7] : 0);
        for(int i=0;i<7;++i){
            valIm[i]=v[i];
            errIm[i]=e[i];
        }

        minuit.mnstat(finalChi2, /*edm*/finalChi2, /*errdef*/finalChi2,
                      /*nvpar*/finalNdf, /*nparx*/finalNdf,
                      /*icstat*/finalNdf);
        // above mnstat args overwritten, so re-call properly:
        double edm, errdef; int nv, nx, ic;
        minuit.mnstat(finalChi2, edm, errdef, nv, nx, ic);
        finalNdf = (gStrategy==1
                  ? int(bsaData.size())-7
                  : int(bsaData.size()+xsData.size())-8);

        // console
        std::cout<<"\n=== Results ===\n"
                 <<" renormImag="<<valIm[0]<<"±"<<errIm[0]<<"\n"
                 <<" alpha0    ="<<valIm[1]<<"±"<<errIm[1]<<"\n"
                 <<" alpha1    ="<<valIm[2]<<"±"<<errIm[2]<<"\n"
                 <<" n_val     ="<<valIm[3]<<"±"<<errIm[3]<<"\n"
                 <<" b_val     ="<<valIm[4]<<"±"<<errIm[4]<<"\n"
                 <<" Mm2_val   ="<<valIm[5]<<"±"<<errIm[5]<<"\n"
                 <<" P_val     ="<<valIm[6]<<"±"<<errIm[6]<<"\n";
        if(gStrategy==2)
            std::cout<<" renormReal="<<valReal<<"±"<<errReal<<"\n";
        std::cout<<" χ²/ndf   ="<<finalChi2<<"/"<<finalNdf
                 <<" ="<<(finalChi2/finalNdf)<<"\n\n";
    }

    // Strategy 3: two‐step
    else {
        // — Stage 1: ImH→BSA
        gStage = 1;
        {
            TMinuit minuit(7);
            minuit.SetPrintLevel(1);
            minuit.SetFCN(fcn);
            minuit.DefineParameter(0,"renormImag",1.0,0.1,0,10);
            minuit.DefineParameter(1,"alpha0",    0.43,0.05,-5,5);
            minuit.DefineParameter(2,"alpha1",    0.85,0.05,-10,10);
            minuit.DefineParameter(3,"n_val",     1.35,0.05,0,10);
            minuit.DefineParameter(4,"b_val",     0.40,0.05,0,10);
            minuit.DefineParameter(5,"Mm2_val",   0.64,0.05,0,10);
            minuit.DefineParameter(6,"P_val",     1.00,0.05,0,10);

            std::cout<<"Stage 1: fitting ImH→BSA…\n";
            minuit.Migrad();

            double v[7], e[7];
            for(int i=0;i<7;++i) minuit.GetParameter(i,v[i],e[i]);
            for(int i=0;i<7;++i){
                valIm[i]=v[i];
                errIm[i]=e[i];
            }

            double chi2, edm, errdef; int nv, nx, ic;
            minuit.mnstat(chi2, edm, errdef, nv, nx, ic);
            finalChi2 = chi2;
            finalNdf  = int(bsaData.size()) - 7;

            std::cout<<"\nStage 1 results:\n"
                     <<" renormImag="<<v[0]<<"±"<<e[0]<<", "
                     <<"alpha0="   <<v[1]<<"±"<<e[1]<<", "
                     <<"alpha1="   <<v[2]<<"±"<<e[2]<<"\n"
                     <<" n_val="    <<v[3]<<"±"<<e[3]<<", "
                     <<"b_val="    <<v[4]<<"±"<<e[4]<<", "
                     <<"Mm2_val="  <<v[5]<<"±"<<e[5]<<", "
                     <<"P_val="    <<v[6]<<"±"<<e[6]<<"\n"
                     <<" χ²/ndf=" <<finalChi2<<"/"<<finalNdf
                     <<"="<<(finalChi2/finalNdf)<<"\n\n";
        }

        // — Stage 2: renormReal→xsec
        gStage = 2;
        {
            TMinuit minuit(1);
            minuit.SetPrintLevel(1);
            minuit.SetFCN(fcn);
            minuit.DefineParameter(0,"renormReal",1.0,0.1,0,10);

            std::cout<<"Stage 2: fitting renormReal→xsec…\n";
            minuit.Migrad();

            double rR,eR;
            minuit.GetParameter(0,rR,eR);
            valReal = rR; errReal = eR;

            double chi2, edm, errdef; int nv, nx, ic;
            minuit.mnstat(chi2, edm, errdef, nv, nx, ic);
            finalChi2 = chi2;
            finalNdf  = int(xsData.size()) - 1;

            std::cout<<"\nStage 2 results:\n"
                     <<" renormReal="<<rR<<"±"<<eR<<"\n"
                     <<" χ²/ndf   ="<<finalChi2<<"/"<<finalNdf
                     <<"="<<(finalChi2/finalNdf)<<"\n\n";
        }
    }

    // make output dir
    system("mkdir -p output");

    // time‐stamp
    time_t now = time(nullptr);
    tm*    lt  = localtime(&now);
    char   tb[32];
    strftime(tb,sizeof(tb),"%Y%m%d_%H%M%S",lt);
    std::string fname = std::string("output/fit_results_")+tb+".txt";
    std::ofstream out(fname);

    out<<"# FitH results\n";
    out<<"timestamp    "<<tb<<"\n";
    out<<"strategy     "<<gStrategy<<"\n";
    out<<"# values: renormImag alpha0 alpha1 n_val b_val Mm2_val P_val renormReal\n";
    out<<valIm[0]<<" "<<valIm[1]<<" "<<valIm[2]<<" "
       <<valIm[3]<<" "<<valIm[4]<<" "<<valIm[5]<<" "
       <<valIm[6]<<" "<<valReal<<"\n";
    out<<"# errors:\n";
    out<<errIm[0]<<" "<<errIm[1]<<" "<<errIm[2]<<" "
       <<errIm[3]<<" "<<errIm[4]<<" "<<errIm[5]<<" "
       <<errIm[6]<<" "<<errReal<<"\n";
    out<<"# chi2 ndof chi2/ndof\n";
    out<<finalChi2<<" "<<finalNdf<<" "<<(finalChi2/finalNdf)<<"\n";
    out.close();

    std::cout<<"Wrote fit results to "<<fname<<"\n";
    return 0;
}