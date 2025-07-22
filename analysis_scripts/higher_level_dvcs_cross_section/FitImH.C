// FitImH.cpp
// ──────────
// Stand-alone program to fit the single overall scale of Im H (renormImag)
// using your BMK_DVCS class (DVCS_xsec.C) and two data files under imports/.
// Uses ROOT’s TMinuit for χ² minimization.
//
// Compile with:
//   g++ FitImH.cpp DVCS_xsec.C $(root-config --cflags --libs) -o FitImH
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

// include your full class & functions
#include "DVCS_xsec.C"

// -----------------------------------------------------------------------------
// Simple struct to hold one data point
struct DataPoint {
    double phi, Q2, xB, t, Eb, A, sigA;
};

// global containers for BSA and cross-section data
static std::vector<DataPoint> bsaData, xsData;

// -----------------------------------------------------------------------------
// LoadData(): read both text files, skipping lines beginning with '#'.
void LoadData() {
    auto readFile = [&](const char* fname, std::vector<DataPoint>& vec){
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

    readFile("imports/rga_prl_bsa.txt",    bsaData);
    readFile("imports/rga_pass1_xsec_2018.txt", xsData);
}

// -----------------------------------------------------------------------------
// fcn(): the TMinuit cost function, computes total χ² for a given renormImag.
void fcn(int & /*npar*/, double * /*grad*/, double &f, double *par, int /*iflag*/) {
    // par[0] = renormImag
    renormImag = par[0];

    // switch on only the H‐term
    hasH  = true;
    hasHt = false;
    hasE  = false;
    hasEt = false;

    double chi2 = 0;

    // loop over BSA points
    for(auto &d : bsaData) {
        BMK_DVCS dvcs(-1, 1, 0, d.Eb, d.xB, d.Q2, d.t, d.phi);
        double modelA = dvcs.BSA();
        double res     = (d.A - modelA)/d.sigA;
        chi2 += res*res;
    }

    // loop over unpolarized cross-section points
    for(auto &d : xsData) {
        BMK_DVCS dvcs(-1, 0, 0, d.Eb, d.xB, d.Q2, d.t, d.phi);
        double modelA = dvcs.CrossSection();
        double res     = (d.A - modelA)/d.sigA;
        chi2 += res*res;
    }

    f = chi2;
}

// -----------------------------------------------------------------------------
// main(): entry point
int main(int argc, char** argv) {
    std::cout << "\n=== FitImH Stand-alone ===\n";

    // 1) load data
    LoadData();
    std::cout << "Read " << bsaData.size()
              << " BSA points and "
              << xsData.size()
              << " cross-section points.\n\n";

    // 2) set up TMinuit for 1 parameter
    TMinuit minuit(1);
    minuit.SetPrintLevel(1);     // 0 = quiet, 1+ = status
    minuit.SetFCN(fcn);

    // define parameter 0 = renormImag, start=1.0, step=0.1, bounds [0,10]
    minuit.DefineParameter(0, "renormImag", 1.0, 0.1, 0.0, 10.0);

    // 3) run Migrad()
    std::cout << "Running Migrad()...\n";
    minuit.Migrad();

    // 4) retrieve & print results
    double renormFit, errRenorm;
    minuit.GetParameter(0, renormFit, errRenorm);

    double chi2, edm, errdef;
    int nvpar, nparx, icstat;
    minuit.mnstat(chi2, edm, errdef, nvpar, nparx, icstat);
    int ndof = static_cast<int>(bsaData.size() + xsData.size()) - 1;

    std::cout << "\n=== Fit Results ===\n"
              << " renormImag = " << renormFit
              << " ± " << errRenorm << "\n"
              << " χ²/ndof    = "
              << chi2 << "/" << ndof
              << " = " << (chi2/ndof) << "\n\n";

    return 0;
}