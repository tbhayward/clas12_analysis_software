// FitImH.cpp
// ──────────
// Stand‐alone program to fit:
//   1) renormImag   (global scale of Im H)
//   2) alpha0       (intercept of α(t))
//   3) alpha1       (slope of  α(t))
// using BMK_DVCS from DVCS_xsec.C and two data files.
// Uses ROOT’s TMinuit for χ² minimization.
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

// include the full DVCS_xsec.C implementation (with your alpha0/alpha1 edits)
#include "DVCS_xsec.C"

// -----------------------------------------------------------------------------
// struct to hold one data point
struct DataPoint {
    double phi, Q2, xB, t, Eb, A, sigA;
};

// globals for the two data sets
static std::vector<DataPoint> bsaData, xsData;

// -----------------------------------------------------------------------------
// LoadData(): read the two text files under imports/, skip lines starting with '#'.
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

    readFile("imports/rga_prl_bsa.txt",      bsaData);
    readFile("imports/rga_pass1_xsec_2018.txt", xsData);
}

// -----------------------------------------------------------------------------
// fcn(): TMinuit cost function → total χ² for (renormImag, alpha0, alpha1)
void fcn(int & /*npar*/,
         double * /*grad*/,
         double &f,
         double *par,
         int /*iflag*/)
{
    // copy Minuit params into our globals
    renormImag = par[0];
    alpha0     = par[1];
    alpha1     = par[2];

    // only H-term active
    hasH  = true;
    hasHt = false;
    hasE  = false;
    hasEt = false;

    double chi2 = 0;

    // BSA points
    for(auto &d : bsaData) {
        BMK_DVCS dvcs(-1, 1, 0, d.Eb, d.xB, d.Q2, d.t, d.phi);
        double modelA = dvcs.BSA();
        double res    = (d.A - modelA)/d.sigA;
        chi2 += res*res;
    }

    // unpolarized cross‐section points
    for(auto &d : xsData) {
        BMK_DVCS dvcs(-1, 0, 0, d.Eb, d.xB, d.Q2, d.t, d.phi);
        double modelA = dvcs.CrossSection();
        double res    = (d.A - modelA)/d.sigA;
        chi2 += res*res;
    }

    f = chi2;
}

// -----------------------------------------------------------------------------
// main(): set up Minuit, run Migrad(), report results.
int main(int /*argc*/, char** /*argv*/) {
    std::cout << "\n=== FitImH Stand‐alone (3 parameters) ===\n";

    // 1) load data
    LoadData();
    std::cout << "Read " << bsaData.size()
              << " BSA points and "
              << xsData.size()
              << " xsec points.\n\n";

    // 2) configure TMinuit(3)
    TMinuit minuit(3);
    minuit.SetPrintLevel(1);
    minuit.SetFCN(fcn);

    // param 0 = renormImag
    minuit.DefineParameter(0, "renormImag", 1.0, 0.1, 0.0,  10.0);
    // param 1 = alpha0 (intercept)
    minuit.DefineParameter(1, "alpha0",     0.43, 0.05, -5.0, 5.0);
    // param 2 = alpha1 (slope)
    minuit.DefineParameter(2, "alpha1",     0.85, 0.05, -10.0, 10.0);

    // 3) minimize
    std::cout << "Running Migrad() with 3 parameters...\n";
    minuit.Migrad();

    // 4) retrieve & print
    double renormFit, errRenorm;
    double a0fit,   erra0;
    double a1fit,   erra1;
    minuit.GetParameter(0, renormFit, errRenorm);
    minuit.GetParameter(1, a0fit,     erra0);
    minuit.GetParameter(2, a1fit,     erra1);

    double chi2, edm, errdef;
    int nvpar, nparx, icstat;
    minuit.mnstat(chi2, edm, errdef, nvpar, nparx, icstat);
    int ndof = static_cast<int>(bsaData.size() + xsData.size()) - 3;

    std::cout << "\n=== Fit Results ===\n"
              << " renormImag = " << renormFit
              << " ± "         << errRenorm << "\n"
              << " alpha0     = " << a0fit
              << " ± "         << erra0    << "\n"
              << " alpha1     = " << a1fit
              << " ± "         << erra1    << "\n"
              << " χ²/ndof    = " << chi2
              << "/" << ndof   << " = "
              << (chi2/ndof)  << "\n\n";

    return 0;
}