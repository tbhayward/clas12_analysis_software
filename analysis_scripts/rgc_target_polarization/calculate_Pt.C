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
#include <iomanip>       // for std::fixed, std::setprecision
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>

// ---------- CONFIGURATION ----------
// 0 = all three periods, 1 = RGC_Su22 only, 2 = RGC_Fa22 only, 3 = RGC_Sp23 only
const int runMode = 1;
// testRun: 0 means process all runs, >0 will restrict to that single run number
const int testRun = 0;  // set to your run of interest, or 0 to do all

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
    return 0.00823729 + 1.62853*x - 1.38493*x*x + 1.07047*x*x*x - 0.747653*x*x*x*x;
}
double ALL_ABD(double x) {
    return 0.0558035 + 1.23137*x - 1.05596*x*x + 1.95783*x*x*x - 1.22263*x*x*x*x;
}

// path to per-run charge & target‐polarity CSV
const char* RUNINFO =
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "asymmetry_extraction/imports/clas12_run_info.csv";

// hard‐coded NH3 ROOT file paths per period
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

    // --- 1) Read runinfo CSV, build chargeMap, chargePlusMap, chargeMinusMap, signMap, targetPolMap ---
    std::ifstream runinfo(RUNINFO);
    if (!runinfo) {
        std::cerr << "[Error] cannot open " << RUNINFO << "\n";
        return 1;
    }
    std::map<int,double> chargeMap, chargePlusMap, chargeMinusMap;
    std::map<int,int>    signMap;
    std::map<int,double> targetPolMap;
    {
        std::string line;
        while (std::getline(runinfo, line)) {
            if (line.empty() || line[0]=='#') continue;
            std::stringstream ss(line);
            int run; double chTotal, chPlus, chMinus, pol_s, pol_e; char comma;
            ss >> run >> comma >> chTotal  >> comma >> chPlus   >> comma >> chMinus  >> comma >>
                pol_s >> comma >> pol_e;
            chargeMap[run]      = chTotal;
            chargePlusMap[run]  = chPlus;
            chargeMinusMap[run] = chMinus;
            signMap[run]        = (pol_s > 0 ? +1 : -1);
            targetPolMap[run]   = pol_s;
            // // Debug: print loaded polarity sign
            // std::cout << "[Debug] run " << run << "  pol_s=" << pol_s << "  signMap=" <<
            //     signMap[run] << "\n";
        }
    }
    std::cout << "[Loaded] " << chargeMap.size()
              << " runs with total/± charge & sign/pol info\n\n";

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
    out << "Run\tPt_GRV\tsigma_GRV\tsys_GRV\tPt_ABD\tsigma_ABD\tsys_ABD\tPt_avg\tavg_sig\tavg_sys\n";

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
        std::vector<std::vector<long>>   Np(nRuns, std::vector<long>(nBins, 0));
        std::vector<std::vector<long>>   Nm(nRuns, std::vector<long>(nBins, 0));
        std::vector<std::vector<double>> depA_sum(nRuns, std::vector<double>(nBins, 0.0));
        std::vector<std::vector<double>> depC_sum(nRuns, std::vector<double>(nBins, 0.0));
        std::vector<std::vector<long>>   dep_cnt  (nRuns, std::vector<long>(nBins, 0));

        // set branches 
        Int_t    runnum;
        Double_t x, y;
        Int_t    helicity;
        Double_t depA, depC;
        tree->SetBranchAddress("runnum",   &runnum);
        tree->SetBranchAddress("x",        &x);
        tree->SetBranchAddress("y",        &y);
        tree->SetBranchAddress("helicity", &helicity);
        tree->SetBranchAddress("DepA", &depA);
        tree->SetBranchAddress("DepC", &depC);

        Long64_t N = tree->GetEntries();
        bool inTest = false;
        std::cout << "  [Scan] " << N << " events\n";
        for (Long64_t i = 0; i < N; ++i) {
            tree->GetEntry(i);
            if (i % 10000000 == 0 && i > 0)
                std::cout << "    event " << i << "/" << N << "\n";

            // skip zero‐helicity entirely
            if (helicity == 0) continue;

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

            // flip helicity sign and combine with target‐sign (current convention: no flip)
            int hel = helicity;
            int sgn = signMap[runnum];

            int bin = std::upper_bound(xB_bins.begin(), xB_bins.end(), x)
                    - xB_bins.begin() - 1;
            if (bin < 0 || bin >= (int)nBins) continue;

            // fill yields
            if (hel * sgn > 0) Np[ridx][bin]++;
            else               Nm[ridx][bin]++;

            // accumulate depolarization averages (over all accepted events)
            depA_sum[ridx][bin] += depA;
            depC_sum[ridx][bin] += depC;
            dep_cnt  [ridx][bin] += 1;
        }
        f.Close();
        std::cout << "  [Done] filled histograms for " << period << "\n\n";

        for (size_t i = 0; i < runs.size(); ++i) {
            int run = runs[i];
            double cp = chargePlusMap[run];
            double cm = chargeMinusMap[run];
            double pb = Pb.at(period);
            double targPol = targetPolMap[run];

            std::vector<double> xv, yg, ye_g, ya, ye_a;
            std::vector<double> grv_bins, abd_bins;  // store per-bin Pt for model sys calc

            for (size_t b = 0; b < nBins - 1; ++b) {  // skip the last xB bin
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
                double delta = p - m;
                double asym  = delta / S;

                // depolarization averages for this run/bin
                double dep_mean_A = 0.0, dep_mean_C = 0.0, depRatio = 1.0;
                if (dep_cnt[i][b] > 0) {
                    dep_mean_A = depA_sum[i][b] / (double)dep_cnt[i][b];
                    dep_mean_C = depC_sum[i][b] / (double)dep_cnt[i][b];
                    if (std::abs(dep_mean_C) > 0.0) depRatio = dep_mean_A / dep_mean_C;
                }

                // model A_LL and effective scale
                double xm    = xMean[b];
                double df    = Df[b], s_df = sDf[b];
                double a_grv = ALL_GRV(xm);
                double a_abd = ALL_ABD(xm);

                double Ag    = a_grv * df * pb;
                double Aa    = a_abd * df * pb;

                // depolarization-corrected Pt from asymmetry
                double Pt_g  = (asym * depRatio) / Ag;
                double Pt_a  = (asym * depRatio) / Aa;

                // print bin details
                std::cout << std::fixed << std::setprecision(3)
                          << "    bin " << b
                          << ": Np="    << raw_p
                          << ", Nm="    << raw_m
                          << ", normNp="<< p
                          << ", normNm="<< m
                          << ", asym="  << asym
                          << ", Df="    << df
                          << ", Pb="    << pb
                          << ", DepA="  << dep_mean_A
                          << ", DepC="  << dep_mean_C
                          << ", DepA/DepC=" << depRatio
                          << ", asym*(DepA/DepC)/(Df*Pb)=" << (asym*depRatio)/(df*pb)
                          << ", pol_tgt="<< targPol
                          << ", A_GRV=" << a_grv
                          << ", A_ABD=" << a_abd
                          << ", Pt_GRV_bin="<< Pt_g
                          << ", Pt_ABD_bin="<< Pt_a
                          << "\n";

                // propagate stats
                double var_p = raw_p / (cp*cp);
                double var_m = raw_m / (cm*cm);

                double dPg_p = (S - delta)/(Ag*S*S);
                double dPg_n = -(S + delta)/(Ag*S*S);
                double var_g = (depRatio*depRatio) * (dPg_p*dPg_p * var_p
                                                    + dPg_n*dPg_n * var_m);
                double err_g = std::sqrt(var_g);
                err_g = std::sqrt(err_g*err_g
                    + std::pow(Pt_g*s_df/df,2)
                    + std::pow(Pt_g*sigma_Pb.at(period)/pb,2));

                double dPa_p = (S - delta)/(Aa*S*S);
                double dPa_n = -(S + delta)/(Aa*S*S);
                double var_a = (depRatio*depRatio) * (dPa_p*dPa_p * var_p
                                                    + dPa_n*dPa_n * var_m);
                double err_a = std::sqrt(var_a);
                err_a = std::sqrt(err_a*err_a
                    + std::pow(Pt_a*s_df/df,2)
                    + std::pow(Pt_a*sigma_Pb.at(period)/pb,2));

                xv .push_back(xm);
                yg .push_back(Pt_g); ye_g.push_back(err_g);
                ya .push_back(Pt_a); ye_a.push_back(err_a);

                grv_bins.push_back(Pt_g);
                abd_bins.push_back(Pt_a);
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

            // enforce sign from signMap
            double signCorr = signMap.count(run) ? signMap.at(run) : 1.0;
            Pt_grv *= signCorr;
            Pt_abd *= signCorr;

            // per-model bin-to-bin systematics (sample stddev over bins)
            auto sample_std = [](const std::vector<double>& v)->double {
                if (v.size() < 2) return 0.0;
                double mean = 0.0;
                for (double val : v) mean += val;
                mean /= (double)v.size();
                double ss = 0.0;
                for (double val : v) {
                    double d = val - mean;
                    ss += d*d;
                }
                return std::sqrt(ss / (double)(v.size()-1));
            };
            double sys_grv = sample_std(grv_bins);
            double sys_abd = sample_std(abd_bins);

            // extra output columns
            double Pt_avg  = 0.5 * (Pt_grv + Pt_abd);
            double avg_sig = std::max(s_grv, s_abd);

            // avg_sys = sqrt( max(sys_GRV, sys_ABD)^2 + std(sys_GRV, sys_ABD)^2 )
            double sys_mean = 0.5 * (sys_grv + sys_abd);
            double sys_std  = std::sqrt( ((sys_grv - sys_mean) * (sys_grv - sys_mean) +
                                          (sys_abd - sys_mean) * (sys_abd - sys_mean)) / 1.0 ); // N-1=1 for 2 samples
            double avg_sys  = std::sqrt(std::pow(std::max(sys_grv, sys_abd), 2) + std::pow(sys_std, 2));

            // write to output in requested order
            out << run << "\t"
                << std::fixed << std::setprecision(3)
                << Pt_grv << "\t" << s_grv << "\t" << sys_grv << "\t"
                << Pt_abd << "\t" << s_abd << "\t" << sys_abd << "\t"
                << Pt_avg << "\t" << avg_sig << "\t" << avg_sys << "\n";

            std::cout << std::fixed << std::setprecision(3)
                      << "    -> Fit Pt_GRV=" << Pt_grv << "±" << s_grv
                      << ", Pt_ABD=" << Pt_abd << "±" << s_abd
                      << ", Pt_avg=" << Pt_avg
                      << ", avg_sig=" << avg_sig
                      << ", avg_sys=" << avg_sys
                      << ", sys_GRV=" << sys_grv
                      << ", sys_ABD=" << sys_abd << "\n\n";
        }
    }

    out.close();
    std::cout << "[Done] wrote output/Pt_by_run.txt\n";
    return 0;
}