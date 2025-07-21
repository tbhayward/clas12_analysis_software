// plot_Pt_by_run.C
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <cmath>

// hard‐coded xB bin edges
static const std::vector<double> xB_bins = {0.0,0.14,0.24,0.34,0.44,0.54,0.64,1.00};

// beam polarizations (fraction) & their STAT uncertainties
static const std::map<std::string,double> Pb = {
  {"RGC_Su22", 0.8384}, {"RGC_Fa22", 0.8372}, {"RGC_Sp23", 0.8040}
};
static const std::map<std::string,double> sigma_Pb = {
  {"RGC_Su22", 0.0086}, {"RGC_Fa22", 0.0045}, {"RGC_Sp23", 0.0061}
};

// model A_LL(x)
double ALL_GRV(double x) {
  return 0.00823729 + 1.62853*x - 1.38493*x*x
         + 1.07047*x*x*x - 0.747653*x*x*x*x;
}
double ALL_ABD(double x) {
  return 0.0558035 + 1.23137*x - 1.05596*x*x
         + 1.95783*x*x*x - 1.22263*x*x*x*x;
}

int main(int argc, char** argv) {
  if (argc!=2) {
    std::cerr << "Usage: ./plot_Pt_by_run <RGC_Su22|RGC_Fa22|RGC_Sp23>\n";
    return 1;
  }
  std::string period = argv[1];
  // map of ROOT files (only NH3 used to grab the tree structure)
  std::map<std::string,std::map<std::string,std::string>> filePaths = {
    {"RGC_Su22", {{"NH3","/work/.../rgc_su22_inb_NH3_eX.root"}}},
    {"RGC_Fa22", {{"NH3","/work/.../rgc_fa22_inb_NH3_eX.root"}}},
    {"RGC_Sp23", {{"NH3","/work/.../rgc_sp23_inb_NH3_eX.root"}}}
  };
  // open tree
  TFile f(filePaths[period]["NH3"].c_str());
  if (!f.IsOpen()) { std::cerr<<"Cannot open file\n"; return 1; }
  TTree* tree = (TTree*)f.Get("PhysicsEvents");
  if (!tree) { std::cerr<<"No PhysicsEvents\n"; return 1; }

  // --- read runinfo CSV ---
  std::ifstream runinfo("/u/home/thayward/.../clas12_run_info.csv");
  std::map<int,double> chargeMap;
  std::map<int,int>    signMap;
  {
    std::string line;
    while (std::getline(runinfo,line)) {
      if (line.empty()||line[0]=='#') continue;
      std::stringstream ss(line);
      int run; double ch, dummy1,dummy2; double pol_sign, pol_err;
      char comma;
      ss>>run>>comma>>ch>>comma>>dummy1>>comma>>dummy2>>comma>>pol_sign>>comma>>pol_err;
      chargeMap[run] = ch;
      signMap[run]   = (pol_sign>0?+1:-1);
    }
  }

  // --- read dilution factors from CSV (single vector of bins) ---
  std::ifstream dfcsv("output/dilution_factor.csv");
  std::vector<double> Df, sDf;
  { std::string header; std::getline(dfcsv,header);
    std::string ln; while(std::getline(dfcsv,ln) && (int)Df.size()<((int)xB_bins.size()-1)) {
      std::stringstream ss(ln);
      double xmean; double df, dfs;
      ss>>xmean; ss.ignore(1); ss>>df; ss.ignore(1); ss>>dfs;
      Df .push_back(df);
      sDf.push_back(dfs);
    }
  }

  // prepare output
  gSystem->mkdir("output",true);
  std::ofstream out("output/Pt_by_run.txt");
  out<<"Run\tPt_GRV\tσ_GRV\tPt_ABD\tσ_ABD\n";

  // set up tree branches
  Int_t   runnum;
  Float_t x, helicity;
  tree->SetBranchAddress("runnum",&runnum);
  tree->SetBranchAddress("x",     &x);
  tree->SetBranchAddress("helicity",&helicity);

  // allocate per‐run per‐bin counters
  struct Counts { std::vector<long> Np, Nm; };
  std::map<int,Counts> cnt;
  for (auto& [run,_]: chargeMap) {
    cnt[run] = { std::vector<long>(xB_bins.size()-1,0),
                 std::vector<long>(xB_bins.size()-1,0) };
  }

  // single pass over all events
  Long64_t N = tree->GetEntries();
  for (Long64_t i=0; i<N; ++i) {
    tree->GetEntry(i);
    auto it = cnt.find(runnum);
    if (it==cnt.end()) continue;
    int sgn = signMap[runnum];
    int bin = std::upper_bound(xB_bins.begin(),xB_bins.end(),x) 
            - xB_bins.begin() - 1;
    if (bin<0 || bin>= (int)xB_bins.size()-1) continue;
    if (helicity*sgn > 0) it->second.Np[bin]++;
    else                  it->second.Nm[bin]++;
  }

  // for each run, build Pt(x) points, fit to constant
  for (auto& [run,counts]: cnt) {
    // prepare graph for GRV
    std::vector<double> xv, yv_grv, yerr_grv, yv_abd, yerr_abd;
    for (int b=0; b<(int)xB_bins.size()-1; ++b) {
      long Np = counts.Np[b], Nm = counts.Nm[b];
      long S  = Np + Nm;
      if (S<1) continue;
      double Δ = (double)Np - (double)Nm;
      double xm = 0.5*(xB_bins[b]+xB_bins[b+1]);
      double df = Df[b], s_df = sDf[b];
      double pb = Pb.at(period), s_pb = sigma_Pb.at(period);
      double Ag = ALL_GRV(xm)*df*pb;
      double Aa = ALL_ABD(xm)*df*pb;
      double Pt_g = Δ/(Ag * S);
      double Pt_a = Δ/(Aa * S);
      // stat error from Poisson + df & pb stat
      double A = Δ, T = S;
      double dPg_p = (T - A)/(Ag*T*T), dPg_n = -(T + A)/(Ag*T*T);
      double varg = dPg_p*dPg_p * Np + dPg_n*dPg_n * Nm;
      double errg = std::sqrt(varg);
      errg = std::sqrt(errg*errg + (Pt_g*s_df/df)*(Pt_g*s_df/df)
                              + (Pt_g*s_pb/pb)*(Pt_g*s_pb/pb));
      double dPa_p = (T - A)/(Aa*T*T), dPa_n = -(T + A)/(Aa*T*T);
      double vara = dPa_p*dPa_p * Np + dPa_n*dPa_n * Nm;
      double erra = std::sqrt(vara);
      erra = std::sqrt(erra*erra + (Pt_a*s_df/df)*(Pt_a*s_df/df)
                              + (Pt_a*s_pb/pb)*(Pt_a*s_pb/pb));
      xv.push_back(xm);
      yv_grv.push_back(Pt_g); yerr_grv.push_back(errg);
      yv_abd.push_back(Pt_a); yerr_abd.push_back(erra);
    }
    if (xv.empty()) continue;

    // fit GRV
    TGraphErrors g_grv(xv.size(), &xv[0], &yv_grv[0],
                       nullptr, &yerr_grv[0]);
    TF1 f0("f0","[0]",0,1);
    g_grv.Fit(&f0,"Q");
    double Pt_grv = f0.GetParameter(0),
           s_grv  = f0.GetParError(0);

    // fit ABD
    TGraphErrors g_abd(xv.size(), &xv[0], &yv_abd[0],
                       nullptr, &yerr_abd[0]);
    TF1 f1("f1","[0]",0,1);
    g_abd.Fit(&f1,"Q");
    double Pt_abd = f1.GetParameter(0),
           s_abd  = f1.GetParError(0);

    out<<run<<"\t"
       <<Pt_grv<<"\t"<<s_grv<<"\t"
       <<Pt_abd<<"\t"<<s_abd<<"\n";
  }

  out.close();
  std::cout<<"[Done] wrote output/Pt_by_run.txt\n";
  return 0;
}