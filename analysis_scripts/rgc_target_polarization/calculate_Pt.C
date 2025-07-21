// plot_Pt_by_run.C

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <Rtypes.h>      // for Int_t, Double_t
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>

// ---------- CONFIGURATION ----------
// 0 = all three periods, 1 = RGC_Su22 only, 2 = RGC_Fa22 only, 3 = RGC_Sp23 only
const int runMode = 1;
// testRun: 0 means process all runs, >0 will restrict to that single run number
const int testRun = 16137;  // set to your run of interest, or 0 to do all

// xB bin edges
static const std::vector<double> xB_bins = {
    0.00, 0.14, 0.24, 0.34, 0.44, 0.54, 0.64, 1.00
};

// beam polarization (fraction) & stat uncertainty
static const std::map<std::string,double> Pb = {
    {"RGC_Su22", 0.8384},
    {"RGC_Fa22", 0.8372},
    {"RGC_Sp23", 0.8040}
};
static const std::map<std::string,double> sigma_Pb = {
    {"RGC_Su22", 0.0086},
    {"RGC_Fa22", 0.0045},
    {"RGC_Sp23", 0.0061}
};

// A_LL models
double ALL_GRV(double x) {
    return 0.00823729 + 1.62853*x
         - 1.38493*x*x + 1.07047*x*x*x
         - 0.747653*x*x*x*x;
}
double ALL_ABD(double x) {
    return 0.0558035 + 1.23137*x
         - 1.05596*x*x + 1.95783*x*x*x
         - 1.22263*x*x*x*x;
}

// path to per-run charge & target-polarity CSV
const char* RUNINFO =
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "asymmetry_extraction/imports/clas12_run_info.csv";

// hard-coded NH3 ROOT file paths per period
static const std::map<std::string,std::string> filePaths = {
    {"RGC_Su22", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_NH3_eX.root"},
    {"RGC_Fa22", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_NH3_eX.root"},
    {"RGC_Sp23", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_NH3_eX.root"}
};

// choose periods to process
std::vector<std::string> selectPeriods() {
    if      (runMode==1) return {"RGC_Su22"};
    else if (runMode==2) return {"RGC_Fa22"};
    else if (runMode==3) return {"RGC_Sp23"};
    else                 return {"RGC_Su22","RGC_Fa22","RGC_Sp23"};
}

int main() {
    auto periods = selectPeriods();
    std::cout << "[Start] processing:";
    for (auto& p : periods) std::cout << " " << p;
    std::cout << "\n\n";

    // --- 1) Read runinfo CSV, build chargeMap, chargePlusMap, chargeMinusMap, signMap ---
    std::ifstream runinfo(RUNINFO);
    if (!runinfo) {
        std::cerr << "[Error] cannot open " << RUNINFO << "\n";
        return 1;
    }
    std::map<int,double> chargeMap, chargePlusMap, chargeMinusMap;
    std::map<int,int>    signMap;
    {
        std::string line;
        while (std::getline(runinfo, line)) {
            if (line.empty() || line[0]=='#') continue;
            std::stringstream ss(line);
            int run; double chTotal, chPlus, chMinus, pol_s, pol_e; char comma;
            ss >> run >> comma
               >> chTotal  >> comma
               >> chPlus   >> comma
               >> chMinus  >> comma
               >> pol_s    >> comma
               >> pol_e;
            chargeMap[run]      = chTotal;
            chargePlusMap[run]  = chPlus;
            chargeMinusMap[run] = chMinus;
            signMap[run]        = (pol_s > 0 ? +1 : -1);
        }
    }
    std::cout << "[Loaded] " << chargeMap.size()
              << " runs with total/± charge & sign info\n\n";

    // --- 2) Load dilution factors and x_mean from CSV ---
    std::ifstream dfcsv("output/dilution_factor.csv");
    if (!dfcsv) {
        std::cerr << "[Error] cannot open output/dilution_factor.csv\n";
        return 1;
    }
    std::vector<double> xMean, Df, sDf;
    {
        std::string hdr; std::getline(dfcsv, hdr);
        std::string ln;
        while (std::getline(dfcsv, ln) && xMean.size() < xB_bins.size()-1) {
            std::stringstream ss(ln);
            double xm, df, dfs; char comma;
            ss >> xm >> comma >> df >> comma >> dfs;
            xMean.push_back(xm);
            Df   .push_back(df);
            sDf  .push_back(dfs);
        }
    }
    std::cout << "[Loaded] " << xMean.size()
              << " x_mean & dilution entries\n\n";

    // --- 3) Prepare output ---
    gSystem->mkdir("output", true);
    std::ofstream out("output/Pt_by_run.txt");
    out << "Run\tPt_GRV\tsigma_GRV\tPt_ABD\tsigma_ABD\n";

    // --- 4) Single-pass per-period loop with early exit for testRun ---
    const size_t nBins = xB_bins.size() - 1;
    for (auto& period : periods) {
        std::cout << "[Period] " << period << "\n";

        // open NH3 file & tree
        TFile f(filePaths.at(period).c_str());
        TTree* tree = (TTree*)f.Get("PhysicsEvents");
        if (!tree) {
            std::cerr << "[Error] PhysicsEvents not found in "
                      << filePaths.at(period) << "\n";
            f.Close();
            continue;
        }

        // build run-index map and allocate counters
        std::map<int,size_t> runToIdx;
        std::vector<int> runs;
        for (auto& kv : chargeMap) {
            runs.push_back(kv.first);
            runToIdx[kv.first] = runs.size() - 1;
        }
        const size_t nRuns = runs.size();
        std::vector<std::vector<long>> Np(nRuns, std::vector<long>(nBins, 0));
        std::vector<std::vector<long>> Nm(nRuns, std::vector<long>(nBins, 0));

        // set branches (including y)
        Int_t    runnum;
        Double_t x, y;
        Int_t    helicity;
        tree->SetBranchAddress("runnum",   &runnum);
        tree->SetBranchAddress("x",        &x);
        tree->SetBranchAddress("y",        &y);
        tree->SetBranchAddress("helicity", &helicity);

        Long64_t N = tree->GetEntries();
        bool inTest = false;
        std::cout << "  [Scan] " << N << " events\n";
        for (Long64_t i = 0; i < N; ++i) {
            tree->GetEntry(i);
            if (i % 10000000 == 0 && i > 0)
                std::cout << "    event " << i << "/" << N << "\n";

            // apply y-cut
            if (y >= 0.75) continue;

            // early exit logic for single-run test
            if (testRun > 0) {
                if (!inTest) {
                    if (runnum == testRun) {
                        inTest = true;
                        std::cout << "    [Entering run " << testRun << "]\n";
                    } else {
                        continue;
                    }
                } else if (runnum != testRun) {
                    std::cout << "    [Leaving run " << testRun << "]\n";
                    break;
                }
            }

            auto it = runToIdx.find(runnum);
            if (it == runToIdx.end()) continue;
            size_t ridx = it->second;
            int sgn = signMap[runnum];

            // ** flip helicity sign here **
            int hel = -helicity;

            int bin = std::upper_bound(xB_bins.begin(), xB_bins.end(), x)
                    - xB_bins.begin() - 1;
            if (bin < 0 || bin >= (int)nBins) continue;

            if (hel * sgn > 0) Np[ridx][bin]++;
            else               Nm[ridx][bin]++;
        }
        f.Close();
        std::cout << "  [Done] filled histograms for " << period << "\n\n";

        // compute Pt per run
        for (size_t i = 0; i < nRuns; ++i) {
            int run = runs[i];
            if (testRun > 0 && run != testRun) continue;

            std::vector<double> xv, yg, ye_g, ya, ye_a;
            std::cout << "  [Compute] Run " << run << "\n";
            double cp = chargePlusMap[run];
            double cm = chargeMinusMap[run];

            for (size_t b = 0; b < nBins; ++b) {
                long raw_p = Np[i][b], raw_m = Nm[i][b];
                if (cp <= 0 || cm <= 0) {
                    std::cout << "    bin " << b << ": missing charge info, skip\n";
                    continue;
                }
                double p = raw_p / cp;
                double m = raw_m / cm;
                double S = p + m;
                if (S < 1e-12) {
                    std::cout << "    bin " << b << ": S≈0, skip\n";
                    continue;
                }
                double Δ  = p - m;
                double xm = xMean[b];
                double df = Df[b], s_df = sDf[b];
                double pb = Pb.at(period), s_pb = sigma_Pb.at(period);

                double a_grv = ALL_GRV(xm);
                double a_abd = ALL_ABD(xm);

                double Ag   = a_grv * df * pb;
                double Aa   = a_abd * df * pb;
                double Pt_g = Δ / (Ag * S);
                double Pt_a = Δ / (Aa * S);

                std::cout
                  << "    bin " << b
                  << ": Np=" << raw_p
                  << ", Nm=" << raw_m
                  << ", normNp=" << p
                  << ", normNm=" << m
                  << ", Df=" << df
                  << ", Pb=" << pb
                  << ", A_GRV=" << a_grv
                  << ", A_ABD=" << a_abd
                  << ", Pt_GRV_bin=" << Pt_g
                  << ", Pt_ABD_bin=" << Pt_a << "\n";

                // propagate stats using sqrt(N)/Q
                double var_p = raw_p / (cp*cp);
                double var_m = raw_m / (cm*cm);

                double dPg_p = (S - Δ)/(Ag*S*S);
                double dPg_n = -(S + Δ)/(Ag*S*S);
                double var_g = dPg_p*dPg_p * var_p
                             + dPg_n*dPg_n * var_m;
                double err_g = std::sqrt(var_g);
                err_g = std::sqrt(err_g*err_g
                    + std::pow(Pt_g*s_df/df,2)
                    + std::pow(Pt_g*s_pb/pb,2));

                double dPa_p = (S - Δ)/(Aa*S*S);
                double dPa_n = -(S + Δ)/(Aa*S*S);
                double var_a = dPa_p*dPa_p * var_p
                             + dPa_n*dPa_n * var_m;
                double err_a = std::sqrt(var_a);
                err_a = std::sqrt(err_a*err_a
                    + std::pow(Pt_a*s_df/df,2)
                    + std::pow(Pt_a*s_pb/pb,2));

                xv .push_back(xm);
                yg .push_back(Pt_g); ye_g.push_back(err_g);
                ya .push_back(Pt_a); ye_a.push_back(err_a);
            }
            if (xv.empty()) {
                std::cout << "    [No valid bins for run " << run << "]\n\n";
                continue;
            }

            // fit GRV constant
            TGraphErrors g_grv(xv.size(), &xv[0], &yg[0], nullptr, &ye_g[0]);
            TF1 fit0("fit0","[0]",0,1);
            g_grv.Fit(&fit0,"Q");
            double Pt_grv = fit0.GetParameter(0), s_grv = fit0.GetParError(0);

            // fit ABD constant
            TGraphErrors g_abd(xv.size(), &xv[0], &ya[0], nullptr, &ye_a[0]);
            TF1 fit1("fit1","[0]",0,1);
            g_abd.Fit(&fit1,"Q");
            double Pt_abd = fit1.GetParameter(0), s_abd = fit1.GetParError(0);

            out << run << "\t"
                << Pt_grv << "\t" << s_grv << "\t"
                << Pt_abd << "\t" << s_abd << "\n";
            std::cout << "    -> Fit Pt_GRV=" << Pt_grv << "±" << s_grv
                      << ", Pt_ABD=" << Pt_abd << "±" << s_abd << "\n\n";
        }
    }

    out.close();
    std::cout << "[Done] wrote output/Pt_by_run.txt\n";
    return 0;
}