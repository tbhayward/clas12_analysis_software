// FitImH.cpp
// ──────────
// Stand-alone program to fit seven parameters of the Im H valence ansatz:
//   0) renormImag (global scale)
//   1) alpha0     (intercept of α(t) = α0 + α1 t)
//   2) alpha1     (slope   of α(t))
//   3) n_val      (exponent power of ξ‐dependence)
//   4) b_val      (power of (1−ξ)/(1+ξ) factor)
//   5) Mm2_val    (mass scale in the t‐dependence)
//   6) P_val      (power of the t‐dependence factor)
// Uses BMK_DVCS in DVCS_xsec.C and two data files (BSA & xsec).
// Minimization via ROOT’s TMinuit.
//
// Compile with:
//   g++ -O2 FitImH.cpp `root-config --cflags --libs` -lMinuit -o FitImH
//
// Run with:
//   ./FitImH

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"

// pull in your full BMK_DVCS implementation (with globals alpha0, alpha1, n_val, etc. declared)
#include "DVCS_xsec.C"

// -----------------------------------------------------------------------------
// small struct for one measurement
struct DataPoint {
    double phi, Q2, xB, t, Eb, A, sigA;
};

// globals holding our two data sets
static std::vector<DataPoint> bsaData, xsData;

// -----------------------------------------------------------------------------
// LoadData(): populate bsaData & xsData from the two imports/*.txt files.
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
            // skip blank lines or headers beginning with '#'
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
void fcn(int & /*npar*/,
         double * /*grad*/,
         double &f,
         double *par,
         int /*iflag*/)
{
    // copy Minuit parameters into our global model switches:
    renormImag = par[0];
    alpha0     = par[1];
    alpha1     = par[2];
    n_val      = par[3];
    b_val      = par[4];
    Mm2_val    = par[5];
    P_val      = par[6];

    // only include the H‐term
    hasH  = true;
    hasHt = false;
    hasE  = false;
    hasEt = false;

    double chi2 = 0;

    // beam‐spin asymmetry data
    for(auto &d : bsaData) {
        BMK_DVCS dvcs(-1, 1, 0,
                     d.Eb, d.xB, d.Q2, d.t, d.phi);
        double modelA = dvcs.BSA();
        double res    = (d.A - modelA)/d.sigA;
        chi2 += res*res;
    }

    // // unpolarized cross‐section data
    // for(auto &d : xsData) {
    //     BMK_DVCS dvcs(-1, 0, 0,
    //                  d.Eb, d.xB, d.Q2, d.t, d.phi);
    //     double modelA = dvcs.CrossSection();
    //     double res    = (d.A - modelA)/d.sigA;
    //     chi2 += res*res;
    // }

    f = chi2;
}

// -----------------------------------------------------------------------------
// main(): configure Minuit, run Migrad(), and report the fit
int main(int /*argc*/, char** /*argv*/) {
    std::cout << "\n=== FitImH Stand‐alone (7 parameters) ===\n";

    // 1) load the data
    LoadData();
    std::cout << "Read " << bsaData.size()
              << " BSA points and "
              << xsData.size()
              << " xsec points.\n\n";

    // 2) set up TMinuit with 7 parameters
    TMinuit minuit(7);
    minuit.SetPrintLevel(1);     // 0 = silent, 1+ shows progress
    minuit.SetFCN(fcn);

    // define each parameter: name, start, step, lower, upper
    minuit.DefineParameter(0, "renormImag", 1.0, 0.1,   0.0,  10.0);
    minuit.DefineParameter(1, "alpha0",     0.43,0.05,  -5.0,   5.0);
    minuit.DefineParameter(2, "alpha1",     0.85,0.05, -10.0,  10.0);
    minuit.DefineParameter(3, "n_val",      1.35,0.05,   0.0,  10.0);
    minuit.DefineParameter(4, "b_val",      0.40,0.05,   0.0,  10.0);
    minuit.DefineParameter(5, "Mm2_val",    0.64,0.05,   0.0,  10.0);
    minuit.DefineParameter(6, "P_val",      1.00,0.05,   0.0,  10.0);

    // 3) run the minimizer
    std::cout << "Running Migrad() with 7 parameters...\n";
    minuit.Migrad();

    // 4) extract & print results
    double renormFit, errRenorm;
    double a0fit, erra0, a1fit, erra1;
    double nfit, errn, bfit, errb;
    double m2fit, errm2, Pfit, errP;

    minuit.GetParameter(0, renormFit, errRenorm);
    minuit.GetParameter(1, a0fit,     erra0);
    minuit.GetParameter(2, a1fit,     erra1);
    minuit.GetParameter(3, nfit,      errn);
    minuit.GetParameter(4, bfit,      errb);
    minuit.GetParameter(5, m2fit,     errm2);
    minuit.GetParameter(6, Pfit,      errP);

    double chi2, edm, errdef;
    int nvpar, nparx, icstat;
    minuit.mnstat(chi2, edm, errdef, nvpar, nparx, icstat);

    int ndof = static_cast<int>(bsaData.size() + xsData.size()) - 7;

    std::cout << "\n=== Fit Results ===\n"
              << " renormImag = " << renormFit
              << " ± " << errRenorm << "\n"
              << " alpha0     = " << a0fit
              << " ± " << erra0    << "\n"
              << " alpha1     = " << a1fit
              << " ± " << erra1    << "\n"
              << " n_val      = " << nfit
              << " ± " << errn     << "\n"
              << " b_val      = " << bfit
              << " ± " << errb     << "\n"
              << " Mm2_val    = " << m2fit
              << " ± " << errm2    << "\n"
              << " P_val      = " << Pfit
              << " ± " << errP     << "\n"
              << " χ²/ndof    = " << chi2
              << "/" << ndof     << " = "
              << (chi2/ndof)  << "\n\n";

    return 0;
}