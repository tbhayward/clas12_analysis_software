// plot_Pt_by_run.C
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <Rtypes.h>      // For Int_t, Double_t
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <cmath>

// ---------- CONFIGURATION ----------
// 0 = process all three periods
// 1 = RGC_Su22 only
// 2 = RGC_Fa22 only
// 3 = RGC_Sp23 only
const int runMode = 1;

// xB bin edges
static const std::vector<double> xB_bins = {0.0,0.14,0.24,0.34,0.44,0.54,0.64,1.00};

// beam polarization (fraction) & stat uncertainty
static const std::map<std::string,double> Pb = {
  {"RGC_Su22", 0.8384}, {"RGC_Fa22", 0.8372}, {"RGC_Sp23", 0.8040}
};
static const std::map<std::string,double> sigma_Pb = {
  {"RGC_Su22", 0.0086}, {"RGC_Fa22", 0.0045}, {"RGC_Sp23", 0.0061}
};

// models for A_LL(x)
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

// path to per-run charge & target‐polarity CSV
const char* RUNINFO =
  "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
  "asymmetry_extraction/imports/clas12_run_info.csv";

// hard‐coded ROOT file paths for each period and target
std::map<std::string,std::map<std::string,std::string>> filePaths = {
  {"RGC_Su22", {
      {"NH3", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_NH3_eX.root"},
      {"C",   "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_C_eX.root"},
      {"CH2", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_CH2_eX.root"},
      {"He",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_He_eX.root"},
      {"ET",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_ET_eX.root"}
  }},
  {"RGC_Fa22", {
      {"NH3", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_NH3_eX.root"},
      {"C",   "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_C_eX.root"},
      {"CH2", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_CH2_eX.root"},
      {"He",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_He_eX.root"},
      {"ET",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_ET_eX.root"}
  }},
  {"RGC_Sp23", {
      {"NH3", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_NH3_eX.root"},
      {"C",   "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_C_eX.root"},
      {"CH2", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_CH2_eX.root"},
      {"He",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_He_eX.root"},
      {"ET",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_ET_eX.root"}
  }}
};

// helper to pick which periods to run
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
  std::cout << std::endl;

  // --- read runinfo CSV ---
  std::ifstream runinfo(RUNINFO);
  if (!runinfo) {
    std::cerr << "[Error] cannot open " << RUNINFO << std::endl;
    return 1;
  }
  std::map<int,double> chargeMap;
  std::map<int,int>   signMap;
  {
    std::string line;
    while (std::getline(runinfo, line)) {
      if (line.empty() || line[0]=='#') continue;
      std::stringstream ss(line);
      int run; double ch, d1, d2, pol_s, pol_e; char comma;
      ss >> run >> comma >> ch >> comma >> d1 >> comma >> d2 >> comma >> pol_s >> comma >> pol_e;
      chargeMap[run] = ch;
      signMap[run]   = (pol_s > 0 ? +1 : -1);
    }
  }
  std::cout << "[Loaded] " << chargeMap.size() << " runs with charge info\n";

  // --- load dilution CSV ---
  std::ifstream dfcsv("output/dilution_factor.csv");
  if (!dfcsv) {
    std::cerr << "[Error] cannot open output/dilution_factor.csv\n";
    return 1;
  }
  std::vector<double> Df, sDf;
  { std::string hdr; std::getline(dfcsv, hdr);
    std::string ln;
    while (std::getline(dfcsv, ln) && Df.size() < xB_bins.size()-1) {
      std::stringstream ss(ln);
      double xm, df, dfs; char comma;
      ss >> xm >> comma >> df >> comma >> dfs;
      Df.push_back(df);
      sDf.push_back(dfs);
    }
  }

  // --- prepare output log ---
  gSystem->mkdir("output", true);
  std::ofstream out("output/Pt_by_run.txt");
  out << "Run\tPt_GRV\tσ_GRV\tPt_ABD\tσ_ABD\n";

  // --- loop over selected periods ---
  for (auto& period : periods) {
    std::cout << "[Period] " << period << std::endl;

    // open NH3 file (structure shared)
    const auto& path = filePaths[period]["NH3"];
    TFile f(path.c_str());
    TTree* tree = (TTree*)f.Get("PhysicsEvents");
    if (!tree) {
      std::cerr << "[Error] no PhysicsEvents in " << path << std::endl;
      f.Close();
      continue;
    }

    // prepare per-run counts per bin
    struct Counts { std::vector<long> Np, Nm; };
    std::map<int,Counts> cnt;
    for (auto& kv : chargeMap) {
      cnt[kv.first] = { std::vector<long>(xB_bins.size()-1,0),
                        std::vector<long>(xB_bins.size()-1,0) };
    }

    // set branches with correct types
    Int_t   runnum;
    Double_t x;
    Int_t   helicity;
    tree->SetBranchAddress("runnum",   &runnum);
    tree->SetBranchAddress("x",        &x);
    tree->SetBranchAddress("helicity", &helicity);

    Long64_t N = tree->GetEntries();
    std::cout << "  [Scan] " << N << " events\n";
    for (Long64_t i = 0; i < N; ++i) {
      tree->GetEntry(i);
      auto it = cnt.find(runnum);
      if (it == cnt.end()) continue;
      int sgn = signMap[runnum];
      int b = std::upper_bound(xB_bins.begin(), xB_bins.end(), x)
            - xB_bins.begin() - 1;
      if (b < 0 || b >= (int)xB_bins.size()-1) continue;
      if (helicity * sgn > 0) it->second.Np[b]++;
      else                    it->second.Nm[b]++;
    }

    // compute Pt per run
    for (auto& kv : cnt) {
      int run = kv.first;
      auto& C = kv.second;
      std::vector<double> xv, yg, ye_ga, ya, ye_ab;
      for (size_t b = 0; b < C.Np.size(); ++b) {
        long Np = C.Np[b], Nm = C.Nm[b];
        long S  = Np + Nm;
        if (S < 1) continue;
        double Δ  = double(Np) - double(Nm);
        double xm = 0.5 * (xB_bins[b] + xB_bins[b+1]);
        double df = Df[b], s_df = sDf[b];
        double pb = Pb.at(period), s_pb = sigma_Pb.at(period);

        double Ag = ALL_GRV(xm) * df * pb;
        double Aa = ALL_ABD(xm) * df * pb;

        double Pt_g = Δ / (Ag * S);
        double Pt_a = Δ / (Aa * S);

        // stat error propagation (Poisson + Pb, Df)
        double dPg_p = (S - Δ) / (Ag * S * S);
        double dPg_n = -(S + Δ) / (Ag * S * S);
        double varg = dPg_p*dPg_p * Np + dPg_n*dPg_n * Nm;
        double eg   = std::sqrt(varg);
        eg = std::sqrt( eg*eg
                      + std::pow(Pt_g * s_df/df, 2)
                      + std::pow(Pt_g * s_pb/pb, 2) );

        double dPa_p = (S - Δ) / (Aa * S * S);
        double dPa_n = -(S + Δ) / (Aa * S * S);
        double vara = dPa_p*dPa_p * Np + dPa_n*dPa_n * Nm;
        double ea   = std::sqrt(vara);
        ea = std::sqrt( ea*ea
                      + std::pow(Pt_a * s_df/df, 2)
                      + std::pow(Pt_a * s_pb/pb, 2) );

        xv.push_back(xm);
        yg.push_back(Pt_g); ye_ga.push_back(eg);
        ya.push_back(Pt_a); ye_ab.push_back(ea);
      }
      if (xv.empty()) continue;

      // fit GRV
      TGraphErrors g_grv(xv.size(), &xv[0], &yg[0], nullptr, &ye_ga[0]);
      TF1 f0("f0","[0]",0,1);
      g_grv.Fit(&f0,"Q");
      double Pt_grv = f0.GetParameter(0), s_grv = f0.GetParError(0);

      // fit ABD
      TGraphErrors g_abd(xv.size(), &xv[0], &ya[0], nullptr, &ye_ab[0]);
      TF1 f1("f1","[0]",0,1);
      g_abd.Fit(&f1,"Q");
      double Pt_abd = f1.GetParameter(0), s_abd = f1.GetParError(0);

      out << run << "\t"
          << Pt_grv << "\t" << s_grv << "\t"
          << Pt_abd << "\t" << s_abd << "\n";
      std::cout << "    Run=" << run
                << " Pt_GRV=" << Pt_grv << "±" << s_grv
                << " Pt_ABD=" << Pt_abd << "±" << s_abd << "\n";
    }

    f.Close();
  }

  out.close();
  std::cout << "[Done] wrote output/Pt_by_run.txt" << std::endl;
  return 0;
}