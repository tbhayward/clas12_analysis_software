// FitImH.cpp
// ──────────
// Stand-alone program to fit eight parameters of the GPD-H ansatz:
//   0) renormImag (global scale of Im H)
//   1) alpha0     (intercept of α(t) = α0 + α1 t)
//   2) alpha1     (slope   of α(t))
//   3) n_val      (ξ-exponent)
//   4) b_val      (power of (1−ξ)/(1+ξ) factor)
//   5) Mm2_val    (mass scale in t-dependence)
//   6) P_val      (power of t-dependence factor)
//   7) renormReal (global scale of Re H)
// Uses BMK_DVCS in DVCS_xsec.C and two data files (BSA & xsec).
// Minimization via ROOT’s TMinuit.

// Compile with:
//   g++ -O2 FitImH.cpp `root-config --cflags --libs` -lMinuit -o FitImH
// Run with:
//   ./FitImH

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"

// pull in the full BMK_DVCS implementation (with globals alpha0, alpha1, etc.)
#include "DVCS_xsec.C"

// -----------------------------------------------------------------------------
// small struct for one measurement
struct DataPoint {
    double phi, Q2, xB, t, Eb, A, sigA;
};

// globals holding our two data sets
static std::vector<DataPoint> bsaData, xsData;

// -----------------------------------------------------------------------------
// LoadData(): populate bsaData & xsData from two imports/*.txt files
void LoadData() {
    auto readFile = [&](const char* fname,
                        std::vector<DataPoint>& vec)
    {
        std::ifstream in(fname);
        if(!in) {
            std::cerr << "ERROR: cannot open " << fname << "\n";
            std::exit(1);
        }
        std::string line;
        while(std::getline(in, line)) {
            if(line.empty() || line[0]=='#') continue;
            std::istringstream iss(line);
            DataPoint d;
            iss >> d.phi >> d.Q2 >> d.xB >> d.t
                >> d.Eb  >> d.A  >> d.sigA;
            vec.push_back(d);
        }
    };

    readFile("imports/rga_prl_bsa.txt",       bsaData);
    readFile("imports/rga_pass1_xsec_2018.txt", xsData);
}

// -----------------------------------------------------------------------------
// fcn(): TMinuit cost function → total χ² over both data sets
void fcn(int & /*npar*/, double * /*grad*/, double &f,
         double *par, int /*iflag*/) {
    // copy Minuit parameters into globals
    renormImag = par[0];
    alpha0     = par[1];
    alpha1     = par[2];
    n_val      = par[3];
    b_val      = par[4];
    Mm2_val    = par[5];
    P_val      = par[6];
    renormReal = par[7];  // new

    hasH  = true;
    hasHt = false;
    hasE  = false;
    hasEt = false;

    double chi2 = 0;

    // BSA data
    for(auto &d : bsaData) {
        BMK_DVCS dvcs(-1, 1, 0, d.Eb, d.xB, d.Q2, d.t, d.phi);
        double modelA = dvcs.BSA();
        double res    = (d.A - modelA)/d.sigA;
        chi2 += res*res;
    }

    // unpolarized xsec data
    for(auto &d : xsData) {
        BMK_DVCS dvcs(-1, 0, 0, d.Eb, d.xB, d.Q2, d.t, d.phi);
        double modelA = dvcs.CrossSection();
        double res    = (d.A - modelA)/d.sigA;
        chi2 += res*res;
    }

    f = chi2;
}

// -----------------------------------------------------------------------------
// main(): configure Minuit, run Migrad(), and report the fit
int main(int /*argc*/, char** /*argv*/) {
    std::cout << "\n=== FitImH Stand-alone (8 parameters) ===\n";

    // 1) load the data
    LoadData();
    std::cout << "Read " << bsaData.size()
              << " BSA points and "
              << xsData.size() << " xsec points.\n\n";

    // 2) set up TMinuit with 8 parameters
    TMinuit minuit(8);
    minuit.SetPrintLevel(1);
    minuit.SetFCN(fcn);

    // define parameters: index, name, start, step, lower, upper
    minuit.DefineParameter(0, "renormImag", 1.0, 0.1,   0.0,  10.0);
    minuit.DefineParameter(1, "alpha0",     0.43,0.05,  -5.0,   5.0);
    minuit.DefineParameter(2, "alpha1",     0.85,0.05, -10.0,  10.0);
    minuit.DefineParameter(3, "n_val",      1.35,0.05,   0.0,  10.0);
    minuit.DefineParameter(4, "b_val",      0.40,0.05,   0.0,  10.0);
    minuit.DefineParameter(5, "Mm2_val",    0.64,0.05,   0.0,  10.0);
    minuit.DefineParameter(6, "P_val",      1.00,0.05,   0.0,  10.0);
    minuit.DefineParameter(7, "renormReal", 1.0, 0.1,   0.0,  10.0);

    // 3) run the minimizer
    std::cout << "Running Migrad() with 8 parameters...\n";
    minuit.Migrad();

    // 4) extract & print results
    double renormI, errI;
    double a0, err0, a1, err1;
    double nv, errn, bv, errb;
    double m2, errm2, Pv, errP;
    double renormR, errR;

    minuit.GetParameter(0, renormI, errI);
    minuit.GetParameter(1, a0,      err0);
    minuit.GetParameter(2, a1,      err1);
    minuit.GetParameter(3, nv,      errn);
    minuit.GetParameter(4, bv,      errb);
    minuit.GetParameter(5, m2,      errm2);
    minuit.GetParameter(6, Pv,      errP);
    minuit.GetParameter(7, renormR, errR);

    double chi2, edm, errdef;
    int nvpar, nparx, icstat;
    minuit.mnstat(chi2, edm, errdef, nvpar, nparx, icstat);
    int ndof = static_cast<int>(bsaData.size() + xsData.size()) - 8;

    std::cout << "\n=== Fit Results ===\n"
              << " renormImag = " << renormI
              << " ± " << errI << "\n"
              << " alpha0     = " << a0
              << " ± " << err0 << "\n"
              << " alpha1     = " << a1
              << " ± " << err1 << "\n"
              << " n_val      = " << nv
              << " ± " << errn << "\n"
              << " b_val      = " << bv
              << " ± " << errb << "\n"
              << " Mm2_val    = " << m2
              << " ± " << errm2<< "\n"
              << " P_val      = " << Pv
              << " ± " << errP << "\n"
              << " renormReal = " << renormR
              << " ± " << errR << "\n"
              << " χ²/ndof    = " << chi2
              << "/" << ndof << " = " << (chi2/ndof)
              << "\n\n";

    return 0;
}