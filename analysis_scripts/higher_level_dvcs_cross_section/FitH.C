// FitH.C
// ──────
// Stand‐alone program to fit the GPD‐H ansatz under three possible strategies:
//   1) fit only ImH parameters to BSA data
//   2) fit ImH parameters + renormReal to BSA + xsec simultaneously
//   3) two‐step: (a) fit ImH→BSA, (b) fit renormReal→xsec
//
// Usage:
//   ./FitH <strategy>
// where <strategy> in {1,2,3} as above.
//
// Compile with:
//   g++ -O2 FitH.C `root-config --cflags --libs` -lMinuit -o FitH
// Run with, e.g.:
//   ./FitH 3

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

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
    readFile("imports/rga_prl_bsa.txt",      bsaData);
    readFile("imports/rga_pass1_xsec_2018.txt", xsData);
}

// -----------------------------------------------------------------------------
// fcn(): Minuit’s chi2 function.  Depending on gStrategy/gStage it
//         either copies the correct subset of par[] → global parameters
//         or loops over BSA and/or xsec
//
void fcn(int& /*npar*/, double* /*grad*/, double& f, double* par, int /*iflag*/){
    // 1) unpack parameters
    if(gStrategy==1 || gStrategy==2 || (gStrategy==3 && gStage==1)){
        // always these first seven: ImH ansatz
        renormImag = par[0];
        alpha0     = par[1];
        alpha1     = par[2];
        n_val      = par[3];
        b_val      = par[4];
        Mm2_val    = par[5];
        P_val      = par[6];
    }
    if(gStrategy==2){
        // simultaneous fit: 8th parameter is renormReal
        renormReal = par[7];
    } else if(gStrategy==3 && gStage==2){
        // two‐step: now only renormReal is floating
        renormReal = par[0];
    }

    // always use only H‐term
    hasH = true;
    hasHt = hasE = hasEt = false;

    double chi2 = 0;

    // 2) BSA loop
    bool doBSA = ( gStrategy==1 || gStrategy==2 || (gStrategy==3 && gStage==1) );
    if(doBSA){
        for(auto &d : bsaData){
            BMK_DVCS dvcs(-1, 1, 0, d.Eb, d.xB, d.Q2, d.t, d.phi);
            double mA = dvcs.BSA();
            double r  = (d.A - mA)/d.sigA; // A is the ALU and sigA is stat uncertainty
            chi2 += r*r;
        }
    }

    // 3) xsec loop
    bool doXS = ( gStrategy==2 || (gStrategy==3 && gStage==2) );
    if(doXS){
        for(auto &d : xsData){
            BMK_DVCS dvcs(-1, 0, 0, d.Eb, d.xB, d.Q2, d.t, d.phi);
            double mA = dvcs.CrossSection();
            double r  = (d.A - mA)/d.sigA; // A is the cross section and sigA is stat uncertainty
            chi2 += r*r;
        }
    }

    f = chi2;
}

// -----------------------------------------------------------------------------
// main(): parse strategy, run the appropriate fits, report results
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

    // load all data
    LoadData();
    std::cout<<"  loaded "<<bsaData.size() <<" BSA and "<<xsData.size() <<" xsec points\n\n";

    // Strategy 1 or 2: single Minuit run
    if(gStrategy==1 || gStrategy==2){
        // choose #params
        int npar = (gStrategy==1? 7 : 8);
        TMinuit minuit(npar);
        minuit.SetPrintLevel(1);
        minuit.SetFCN(fcn);

        // define params 0–6 (ImH)
        minuit.DefineParameter(0,"renormImag",1.0,0.1,0,10);
        minuit.DefineParameter(1,"alpha0",    0.43,0.05,-5,5);
        minuit.DefineParameter(2,"alpha1",    0.85,0.05,-10,10);
        minuit.DefineParameter(3,"n_val",     1.35,0.05,0,10);
        minuit.DefineParameter(4,"b_val",     0.40,0.05,0,10);
        minuit.DefineParameter(5,"Mm2_val",   0.64,0.05,0,10);
        minuit.DefineParameter(6,"P_val",     1.00,0.05,0,10);
        if(gStrategy==2){
            // simultaneous ReH
            minuit.DefineParameter(7,"renormReal",1.0,0.1,0,10);
        }

        std::cout<<"Running Migrad()…\n";
        minuit.Migrad();

        // retrieve
        double val[8], err[8];
        for(int i=0;i<npar;++i){
            minuit.GetParameter(i,val[i],err[i]);
        }
        double chi2,edm,errdef; int nv,nx,ic;
        minuit.mnstat(chi2,edm,errdef,nv,nx,ic);
        int ndof = (gStrategy==1 ? int(bsaData.size())-7 : int(bsaData.size()+xsData.size())-8);

        // print out
        std::cout<<"\n=== Results ===\n";
        std::cout<<" renormImag = "<<val[0]<<" ± "<<err[0]<<"\n";
        std::cout<<" alpha0     = "<<val[1]<<" ± "<<err[1]<<"\n";
        std::cout<<" alpha1     = "<<val[2]<<" ± "<<err[2]<<"\n";
        std::cout<<" n_val      = "<<val[3]<<" ± "<<err[3]<<"\n";
        std::cout<<" b_val      = "<<val[4]<<" ± "<<err[4]<<"\n";
        std::cout<<" Mm2_val    = "<<val[5]<<" ± "<<err[5]<<"\n";
        std::cout<<" P_val      = "<<val[6]<<" ± "<<err[6]<<"\n";
        if(gStrategy==2){
            std::cout<<" renormReal = "<<val[7]<<" ± "<<err[7]<<"\n";
        }
        std::cout<<" χ²/ndof    = "<<chi2<<"/"<<ndof <<" = "<<(chi2/ndof)<<"\n\n";
    }
    // Strategy 3: two‐step
    else {
        // --- stage 1: ImH → BSA
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

            // get & print
            double val[7], err[7];
            for(int i=0;i<7;++i) minuit.GetParameter(i,val[i],err[i]);
            double chi2,edm,errdef; int nv,nx,ic;
            minuit.mnstat(chi2,edm,errdef,nv,nx,ic);
            int ndof = int(bsaData.size())-7;
            std::cout<<"\nStage 1 results:\n";
            std::cout<<" renormImag="<<val[0]<<"±"<<err[0]<<", "
                     <<"alpha0="<<val[1]<<"±"<<err[1]<<", "
                     <<"alpha1="<<val[2]<<"±"<<err[2]<<"\n"
                     <<" n_val="<<val[3]<<"±"<<err[3]<<", "
                     <<"b_val="<<val[4]<<"±"<<err[4]<<", "
                     <<"Mm2_val="<<val[5]<<"±"<<err[5]<<", "
                     <<"P_val="<<val[6]<<"±"<<err[6]<<"\n"
                     <<" χ²/ndof="<<chi2<<"/"<<ndof
                     <<"="<<(chi2/ndof)<<"\n\n";
        }

        // --- stage 2: renormReal → xsec
        gStage = 2;
        {
            TMinuit minuit(1);
            minuit.SetPrintLevel(1);
            minuit.SetFCN(fcn);
            minuit.DefineParameter(0,"renormReal",1.0,0.1,0,10);

            std::cout<<"Stage 2: fitting renormReal→xsec…\n";
            minuit.Migrad();

            double rR, eR;
            minuit.GetParameter(0,rR,eR);
            double chi2,edm,errdef; int nv,nx,ic;
            minuit.mnstat(chi2,edm,errdef,nv,nx,ic);
            int ndof = int(xsData.size())-1;
            std::cout<<"\nStage 2 results:\n" <<" renormReal="<<rR<<"±"<<eR<<"\n"
                     <<" χ²/ndof="<<chi2<<"/"<<ndof <<"="<<(chi2/ndof)<<"\n\n";
        }
    }

    return 0;
}