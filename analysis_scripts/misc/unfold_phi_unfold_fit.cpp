// unfold_phi_unfold_fit.cpp
// Exclusive pi+ unfolding and cos(n*phi) fits (phi plotted in degrees).
// Produces AUU^{cos phi} and AUU^{cos2 phi} per x_B bin (raw from fit),
// and scaled FUU ratios vs x_B using depolarization factors.
// Text outputs mirror your "chi2Fits_*" style: propertyNameAUUcosphi = {...}; etc.
//
// Build (csh):
//   g++ -O2 -std=c++17 unfold_phi_unfold_fit.cpp `root-config --cflags --libs` -o unfold_phi_unfold_fit
//
// Run:
//   ./unfold_phi_unfold_fit data.root gen.root rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]
//        [--x x] [--phi phi] [--t t] [--mx2 Mx2] [--fid fiducial_status]
//        [--DepA DepA] [--DepB DepB] [--DepV DepV] [--DepW DepW]
//        [--phibins 24] [--no-fid-cut] [--debug N] [--list-branches]
//
// Outputs per property (to output/enpi+/):
//   - PDF  : <property>_unfolded.pdf
//   - TEXT : <property>_unfolded_fit_arrays.txt
//            propertyNameAUUcosphi = { {xB, val, err}, ... };
//            propertyNameAUUcos2phi = { {xB, val, err}, ... };

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSystem.h"
#include "TTree.h"

struct BranchNames {
  std::string x   = "x";               // x_B
  std::string phi = "phi";             // radians in input
  std::string t   = "t";               // signed t
  std::string mx2 = "Mx2";
  std::string fid = "fiducial_status";
  // Depolarization factors (present on DATA)
  std::string DepA = "DepA";
  std::string DepB = "DepB";
  std::string DepV = "DepV";
  std::string DepW = "DepW";
};

struct Config {
  std::string data_path;
  std::string gen_path;
  std::string rec_path;
  std::string which_property = "all";
  BranchNames bn;
  int  phi_nbins      = 24;
  bool apply_fid_cut  = true;
  int  debug_print    = 0;
  bool list_branches  = false;
};

static void print_usage(const char* argv0) {
  std::cout << "Usage:\n  " << argv0
            << " data.root gen.root rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]\n"
            << "    [--x x] [--phi phi] [--t t] [--mx2 Mx2] [--fid fiducial_status]\n"
            << "    [--DepA DepA] [--DepB DepB] [--DepV DepV] [--DepW DepW]\n"
            << "    [--phibins 24] [--no-fid-cut] [--debug N] [--list-branches]\n";
}

static bool parse_args(int argc, char** argv, Config& cfg) {
  if (argc < 5) { print_usage(argv[0]); return false; }
  cfg.data_path = argv[1];
  cfg.gen_path  = argv[2];
  cfg.rec_path  = argv[3];
  cfg.which_property = argv[4];

  for (int i=5;i<argc;i++) {
    std::string a = argv[i];
    auto eat = [&](const char* key, std::string& out) {
      size_t eq = a.find('=');
      if (a == key) {
        if (i+1 < argc) { out = argv[++i]; return true; }
      } else if (eq != std::string::npos && a.substr(0,eq) == key) {
        out = a.substr(eq+1); return true;
      }
      return false;
    };
    auto eat_int = [&](const char* key, int& out) {
      size_t eq = a.find('=');
      if (a == key) {
        if (i+1 < argc) { out = std::atoi(argv[++i]); return true; }
      } else if (eq != std::string::npos && a.substr(0,eq) == key) {
        out = std::atoi(a.substr(eq+1).c_str()); return true;
      }
      return false;
    };

    if (eat("--x",   cfg.bn.x)) continue;
    if (eat("--phi", cfg.bn.phi)) continue;
    if (eat("--t",   cfg.bn.t)) continue;
    if (eat("--mx2", cfg.bn.mx2)) continue;
    if (eat("--fid", cfg.bn.fid)) continue;
    if (eat("--DepA", cfg.bn.DepA)) continue;
    if (eat("--DepB", cfg.bn.DepB)) continue;
    if (eat("--DepV", cfg.bn.DepV)) continue;
    if (eat("--DepW", cfg.bn.DepW)) continue;
    if (eat_int("--phibins", cfg.phi_nbins)) continue;
    if (a == "--no-fid-cut") { cfg.apply_fid_cut = false; continue; }
    if (eat_int("--debug", cfg.debug_print)) continue;
    if (a == "--list-branches") { cfg.list_branches = true; continue; }

    std::cerr << "Unknown argument: " << a << "\n";
  }
  return true;
}

static std::vector<std::string> all_props() {
  return {"enpi", "enpiLowt", "enpiMidt", "enpiHight"};
}
static std::vector<std::string> properties_to_run(const std::string& which) {
  if (which == "all") return all_props();
  return {which};
}

static std::vector<double> x_edges() {
  return {0.06, 0.14, 0.22, 0.30, 0.40, 0.50, 0.60};
}
static int xbin_index(double x, const std::vector<double>& e) {
  for (size_t i=0;i+1<e.size();++i) if (x >= e[i] && x < e[i+1]) return (int)i;
  return -1;
}

static inline bool pass_common(double t, double mx2) {
  return (std::fabs(t) <= 1.0) && (std::fabs(t) >= 0.0) && (mx2 > 0.80) && (mx2 < 1.00);
}
static inline bool pass_prop_window(const std::string& prop, double t) {
  const double at = std::fabs(t);
  if (prop == "enpi")      return true;
  if (prop == "enpiLowt")  return (at >= 0.00 && at <= 0.30);
  if (prop == "enpiMidt")  return (at >= 0.30 && at <= 0.70);
  if (prop == "enpiHight") return (at >= 0.70);
  return false;
}

static std::string prop_title(const std::string& p) {
  if (p == "enpi")      return "ep -> e' n #pi^{+}";
  if (p == "enpiLowt")  return "ep -> e' n #pi^{+}, low |t|";
  if (p == "enpiMidt")  return "ep -> e' n #pi^{+}, mid |t|";
  if (p == "enpiHight") return "ep -> e' n #pi^{+}, high |t|";
  return p;
}
static std::string prop_tlabel(const std::string& prop) {
  if (prop == "enpi")      return "-t \\in [0.00, 1.00]";
  if (prop == "enpiLowt")  return "-t \\in [0.00, 0.30]";
  if (prop == "enpiMidt")  return "-t \\in [0.30, 0.70]";
  if (prop == "enpiHight") return "-t \\in [0.70, 1.00]";
  return "";
}

static void ensure_dir(const std::string& d) {
  if (gSystem->AccessPathName(d.c_str())) gSystem->mkdir(d.c_str(), kTRUE);
}

static void list_tree_branches(TTree* tr, const char* label) {
  if (!tr) { std::cout << "[" << label << "] Tree is null.\n"; return; }
  std::cout << "\n[" << label << "] Branches on " << tr->GetName() << ":\n";
  TObjArray* bl = tr->GetListOfBranches();
  if (!bl) { std::cout << "  (no branches?)\n"; return; }
  for (int i=0;i<bl->GetEntries();++i) {
    TBranch* br = dynamic_cast<TBranch*>(bl->At(i));
    if (!br) continue;
    std::cout << "  - " << br->GetName() << "\n";
  }
  std::cout.flush();
}

// ------------ Binding ------------
struct BranchHandles {
  double x=0, phi=0, t=0, mx2=0;
  int fid=0;
  bool has_fid=false;
  // Dep factors on DATA
  double DepA=0, DepB=0, DepV=0, DepW=0;
  bool has_DepA=false, has_DepB=false, has_DepV=false, has_DepW=false;
};

static bool bind_tree(TTree* tr, const BranchNames& bn, BranchHandles& bh, const char* who) {
  if (!tr) return false;
  bool ok = true;
  ok &= (tr->GetBranch(bn.x.c_str())  != nullptr);
  ok &= (tr->GetBranch(bn.phi.c_str())!= nullptr);
  ok &= (tr->GetBranch(bn.t.c_str())  != nullptr);
  ok &= (tr->GetBranch(bn.mx2.c_str())!= nullptr);
  if (!ok) {
    std::cerr << "ERROR: missing one of required branches on " << who << " (need: "
              << bn.x << ", " << bn.phi << ", " << bn.t << ", " << bn.mx2 << ")\n";
    return false;
  }
  tr->SetBranchAddress(bn.x.c_str(),   &bh.x);
  tr->SetBranchAddress(bn.phi.c_str(), &bh.phi);
  tr->SetBranchAddress(bn.t.c_str(),   &bh.t);
  tr->SetBranchAddress(bn.mx2.c_str(), &bh.mx2);

  bh.has_fid = (tr->GetBranch(bn.fid.c_str()) != nullptr);
  if (bh.has_fid) tr->SetBranchAddress(bn.fid.c_str(), &bh.fid);

  // Optional dep branches (on DATA)
  bh.has_DepA = (tr->GetBranch(bn.DepA.c_str()) != nullptr);
  bh.has_DepB = (tr->GetBranch(bn.DepB.c_str()) != nullptr);
  bh.has_DepV = (tr->GetBranch(bn.DepV.c_str()) != nullptr);
  bh.has_DepW = (tr->GetBranch(bn.DepW.c_str()) != nullptr);
  if (bh.has_DepA) tr->SetBranchAddress(bn.DepA.c_str(), &bh.DepA);
  if (bh.has_DepB) tr->SetBranchAddress(bn.DepB.c_str(), &bh.DepB);
  if (bh.has_DepV) tr->SetBranchAddress(bn.DepV.c_str(), &bh.DepV);
  if (bh.has_DepW) tr->SetBranchAddress(bn.DepW.c_str(), &bh.DepW);

  return true;
}

// ------------ Hists & dep means ------------
struct HSet {
  std::vector<std::unique_ptr<TH1D>> D; // data
  std::vector<std::unique_ptr<TH1D>> G; // gen
  std::vector<std::unique_ptr<TH1D>> R; // rec
};
static std::map<std::string,HSet> make_hist_map(int phibins) {
  std::map<std::string,HSet> m;
  const auto props = all_props();
  const auto xe = x_edges();
  const int NX = (int)xe.size()-1;
  const double pmin = 0.0, pmax = 360.0;

  for (const auto& p : props) {
    HSet hs;
    hs.D.reserve(NX); hs.G.reserve(NX); hs.R.reserve(NX);
    for (int i=0;i<NX;i++) {
      auto mk = [&](const std::string& tag)->std::unique_ptr<TH1D> {
        std::string name = tag + "_" + p + "_xbin" + std::to_string(i);
        auto h = std::make_unique<TH1D>(name.c_str(), "", phibins, pmin, pmax);
        h->Sumw2();
        return h;
      };
      hs.D.emplace_back(mk("hD"));
      hs.G.emplace_back(mk("hG"));
      hs.R.emplace_back(mk("hR"));
    }
    m.emplace(p, std::move(hs));
  }
  return m;
}

struct DepMeans {
  std::vector<double> sumA, sumB, sumV, sumW;
  std::vector<long long> count;
};
static std::map<std::string, DepMeans> make_dep_map() {
  std::map<std::string, DepMeans> dm;
  const auto props = all_props();
  const int NX = (int)x_edges().size()-1;
  for (const auto& p : props) {
    DepMeans d;
    d.sumA.assign(NX, 0.0);
    d.sumB.assign(NX, 0.0);
    d.sumV.assign(NX, 0.0);
    d.sumW.assign(NX, 0.0);
    d.count.assign(NX, 0);
    dm.emplace(p, std::move(d));
  }
  return dm;
}

// ------------ Debug counters ------------
struct Counters {
  Long64_t total=0, pass_common_cnt=0, pass_fid_cnt=0;
  Long64_t fail_abs_t=0, fail_mx2=0, fail_no_fid_branch=0, fail_fid_lt_100=0, fail_xbin=0;
  std::map<std::string, Long64_t> pass_prop_cnt;
};

static void debug_event_print(const char* label, int idx, const BranchHandles& b,
                              const std::vector<std::string>& props, const std::vector<double>& xe) {
  std::cout << "[" << label << " evt " << idx << "] "
            << "x_B=" << b.x << "  t=" << b.t << "  Mx2=" << b.mx2
            << "  fid=" << (b.has_fid? b.fid : -9999)
            << "  phi(rad)=" << b.phi
            << "  phi(deg)=" << (b.phi*180.0/TMath::Pi()) << "\n";
  for (const auto& p : props) {
    bool in_prop = pass_prop_window(p, b.t);
    int ib = xbin_index(b.x, xe);
    std::cout << "    prop=" << p << "  pass=" << (in_prop?"Y":"N") << "  xbin=" << ib << "\n";
  }
  std::cout.flush();
}

// ------------ Loop & fill ------------
static void loop_tree_fill(
  TTree* tr, const BranchNames& bn, bool apply_fid, bool count_fid,
  std::map<std::string,HSet>& H, std::map<std::string,DepMeans>* depPtr,
  Counters& C, int debugN, const char* dbg_label)
{
  if (!tr) return;
  BranchHandles b;
  if (!bind_tree(tr, bn, b, dbg_label)) {
    std::cerr << "ERROR: failed to bind branches on " << dbg_label << "\n";
    return;
  }
  if (apply_fid && !b.has_fid) {
    std::cerr << "WARNING: " << dbg_label << " missing fiducial branch '" << bn.fid << "'. Fiducial cut cannot be applied.\n";
  }

  const auto xe = x_edges();
  const auto props = all_props();
  const Long64_t nent = tr->GetEntries();
  const int debug_limit = std::max(0, debugN);

  for (Long64_t i=0;i<nent;i++) {
    tr->GetEntry(i);
    C.total++;

    if (debug_limit > 0 && i < debug_limit) {
      debug_event_print(dbg_label, (int)i, b, props, xe);
    }

    const bool ok_abs_t = (std::fabs(b.t) <= 1.0) && (std::fabs(b.t) >= 0.0);
    const bool ok_mx2   = (b.mx2 > 0.80) && (b.mx2 < 1.00);
    if (!ok_abs_t) C.fail_abs_t++;
    if (!ok_mx2)   C.fail_mx2++;
    if (!(ok_abs_t && ok_mx2)) continue;
    C.pass_common_cnt++;

    if (apply_fid) {
      if (!b.has_fid) { C.fail_no_fid_branch++; continue; }
      if (b.fid < 100) { C.fail_fid_lt_100++; continue; }
    }
    if (count_fid) C.pass_fid_cnt++;

    int ib = xbin_index(b.x, xe);
    if (ib < 0) { C.fail_xbin++; continue; }

    double phideg = b.phi * 180.0 / TMath::Pi();
    phideg = std::fmod(phideg, 360.0);
    if (phideg < 0) phideg += 360.0;

    for (const auto& p : props) {
      if (!pass_prop_window(p, b.t)) continue;
      C.pass_prop_cnt[p]++;

      if (std::string(dbg_label) == "DATA") {
        H[p].D[ib]->Fill(phideg);
        if (depPtr) {
          auto& dep = depPtr->at(p);
          if (b.has_DepA) dep.sumA[ib] += b.DepA;
          if (b.has_DepB) dep.sumB[ib] += b.DepB;
          if (b.has_DepV) dep.sumV[ib] += b.DepV;
          if (b.has_DepW) dep.sumW[ib] += b.DepW;
          dep.count[ib] += 1;
        }
      } else if (std::string(dbg_label) == "GEN") {
        H[p].G[ib]->Fill(phideg);
      } else if (std::string(dbg_label) == "REC") {
        H[p].R[ib]->Fill(phideg);
      }
    }

    if ((i>0) && (i % 10000000 == 0)) {
      std::cout << "[" << dbg_label << "] processed " << i << " / " << nent << " entries..." << std::endl;
    }
  }
}

// ------------ Fit ------------
struct FitResult {
  double C=0, dC=0;
  double A=0, dA=0;
  double B=0, dB=0;
  double chi2=0; int ndf=0;
  int npoints=0;
};

static FitResult make_unfold_graph_and_fit(
  TH1D* hD, TH1D* hG, TH1D* hR, TGraphErrors*& g, TF1*& ffit)
{
  FitResult fr;
  if (!hD || !hG || !hR) { g=nullptr; ffit=nullptr; return fr; }

  const int nb = hD->GetNbinsX();
  std::vector<double> X, Y, EY;
  X.reserve(nb); Y.reserve(nb); EY.reserve(nb);

  for (int b=1; b<=nb; ++b) {
    const double D = hD->GetBinContent(b);
    const double G = hG->GetBinContent(b);
    const double R = hR->GetBinContent(b);
    const double xc = hD->GetBinCenter(b); // degrees

    if (R <= 0.0 || G <= 0.0) continue;

    const double U = (D * G) / R;
    const double term1 = (G*G)/(R*R) * D;
    const double term2 = (D*D)/(R*R) * G;
    const double term3 = (D*D*G*G)/(R*R*R);
    const double sig = std::sqrt(std::max(0.0, term1 + term2 + term3));

    X.push_back(xc); Y.push_back(U); EY.push_back(sig);
  }

  fr.npoints = (int)X.size();
  if (fr.npoints == 0) { g=nullptr; ffit=nullptr; return fr; }

  g = new TGraphErrors(fr.npoints);
  for (int i=0;i<fr.npoints;i++) {
    g->SetPoint(i, X[i], Y[i]);
    g->SetPointError(i, 0.0, EY[i]);
  }
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  g->SetLineWidth(1);

  const double xmin = 0.0, xmax = 360.0;
  ffit = new TF1("fitCAB_deg",
                 "[0]*(1 + [1]*cos(TMath::Pi()/180.0*x) + [2]*cos(2.0*TMath::Pi()/180.0*x))",
                 xmin, xmax);
  ffit->SetParNames("C","A","B");
  double ymean = 0.0; for (double v : Y) ymean += v; ymean /= (double)Y.size();
  ffit->SetParameters( (ymean>0?ymean:1.0), 0.0, 0.0 );
  ffit->SetParLimits(0, 0.0, 1e12);
  ffit->SetParLimits(1, -2.0, 2.0);
  ffit->SetParLimits(2, -2.0, 2.0);

  g->Fit(ffit, "Q"); // quiet

  fr.C  = ffit->GetParameter(0);
  fr.dC = ffit->GetParError(0);
  fr.A  = ffit->GetParameter(1);
  fr.dA = ffit->GetParError(1);
  fr.B  = ffit->GetParameter(2);
  fr.dB = ffit->GetParError(2);
  fr.chi2 = ffit->GetChisquare();
  fr.ndf  = ffit->GetNDF();
  return fr;
}

// ------------ Save text in your required style ------------
static void save_arrays_style(
  const std::string& prop,
  const std::vector<double>& xcenters,
  const std::vector<double>& Fcos, const std::vector<double>& dFcos,
  const std::vector<double>& Fcos2, const std::vector<double>& dFcos2)
{
  ensure_dir("output/enpi+");
  std::string out = "output/enpi+/" + prop + "_unfolded_fit_arrays.txt";
  std::ofstream ofs(out);
  ofs << std::fixed << std::setprecision(9);

  // propertyNameAUUcosphi = { {xB, val, err}, ... };
  ofs << prop << "AUUcosphi = {";
  for (size_t i=0;i<xcenters.size();++i) {
    ofs << "{" << xcenters[i] << ", " << Fcos[i] << ", " << dFcos[i] << "}";
    if (i + 1 < xcenters.size()) ofs << ", ";
  }
  ofs << "};\n";

  // propertyNameAUUcos2phi = { {xB, val, err}, ... };
  ofs << prop << "AUUcos2phi = {";
  for (size_t i=0;i<xcenters.size();++i) {
    ofs << "{" << xcenters[i] << ", " << Fcos2[i] << ", " << dFcos2[i] << "}";
    if (i + 1 < xcenters.size()) ofs << ", ";
  }
  ofs << "};\n";

  ofs.close();
  std::cout << "Wrote arrays: " << out << "\n";
}

// ------------ Draw/save per property ------------
static void draw_property_and_save(
  const std::string& prop,
  const std::map<std::string,HSet>& H,
  const std::map<std::string,DepMeans>& depMap)
{
  const auto xe = x_edges();
  const int NX = (int)xe.size()-1;

  // Global style (more left & bottom margin so nothing clips)
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadBottomMargin(0.13); // tiny extra bottom padding
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  TCanvas* c = new TCanvas(("c_"+prop).c_str(), ("Unfolded phi fits: "+prop).c_str(), 1200, 900);
  c->Divide(3,3);

  std::vector<double> xcenters, xlos, xhis;
  std::vector<double> Avec, dAvec, Bvec, dBvec, Cvec, dCvec;
  std::vector<double> meanDepA(NX,0.0), meanDepB(NX,0.0), meanDepV(NX,0.0);

  const auto& dep = depMap.at(prop);

  for (int ib=0; ib<NX; ++ib) {
    c->cd(ib+1);

    TH1D* hD = H.at(prop).D[ib].get();
    TH1D* hG = H.at(prop).G[ib].get();
    TH1D* hR = H.at(prop).R[ib].get();

    TGraphErrors* g = nullptr;
    TF1* fit = nullptr;
    FitResult fr = make_unfold_graph_and_fit(hD,hG,hR,g,fit);

    const double xl = xe[ib], xh = xe[ib+1];
    std::string titleCore = prop_title(prop);
    std::string tlabel = prop_tlabel(prop);
    std::string title = titleCore + Form(", x_{B} in [%.2f, %.2f), %s", xl, xh, tlabel.c_str());

    // compute mean dep factors for this bin from DATA
    if (dep.count[ib] > 0) {
      meanDepA[ib] = dep.sumA[ib] / (double)dep.count[ib];
      meanDepB[ib] = dep.sumB[ib] / (double)dep.count[ib];
      meanDepV[ib] = dep.sumV[ib] / (double)dep.count[ib];
    } else {
      meanDepA[ib] = meanDepB[ib] = meanDepV[ib] = 0.0;
    }

    if (g) {
      g->SetTitle((title + ";#phi (deg);Unfolded yield").c_str());
      g->Draw("AP");
      if (fit) fit->Draw("LSAME");

      // Legend: AUU^{cos phi}, AUU^{cos2 phi} (raw A,B) – small, bordered, inside pad
      auto leg = new TLegend(0.56, 0.72, 0.94, 0.92);
      leg->SetBorderSize(1);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.035);
      leg->AddEntry((TObject*)0, Form("A_{UU}^{#cos#phi} = %.3g #pm %.3g", fr.A, fr.dA), "");
      leg->AddEntry((TObject*)0, Form("A_{UU}^{#cos2#phi} = %.3g #pm %.3g", fr.B, fr.dB), "");
      leg->Draw();

      xcenters.push_back(0.5*(xl+xh));
      xlos.push_back(xl);
      xhis.push_back(xh);
      Avec.push_back(fr.A); dAvec.push_back(fr.dA);
      Bvec.push_back(fr.B); dBvec.push_back(fr.dB);
      Cvec.push_back(fr.C); dCvec.push_back(fr.dC);
    } else {
      TH1D* frame = new TH1D("frame","",10,0,360);
      frame->SetTitle((title + ";#phi (deg);Unfolded yield").c_str());
      frame->SetMinimum(0); frame->SetMaximum(1);
      frame->Draw("AXIS");
      auto leg = new TLegend(0.56, 0.82, 0.94, 0.92);
      leg->SetBorderSize(1); leg->SetFillStyle(0); leg->SetTextSize(0.035);
      leg->AddEntry((TObject*)0, "No valid unfolded points", "");
      leg->Draw();
    }
  }

  // Build scaled FUU ratios from means: Fcos = (DepA/DepV) * A   and   Fcos2 = (DepA/DepB) * B
  std::vector<double> Fcos, dFcos, Fcos2, dFcos2;
  Fcos.reserve(xcenters.size()); dFcos.reserve(xcenters.size());
  Fcos2.reserve(xcenters.size()); dFcos2.reserve(xcenters.size());
  for (size_t i=0;i<xcenters.size();++i) {
    const int ib = (int)i; // aligns with push order
    const double a = Avec[i], da = dAvec[i];
    const double b = Bvec[i], db = dBvec[i];
    const double mA = meanDepA[ib], mB = meanDepB[ib], mV = meanDepV[ib];

    double f1=0, df1=0, f2=0, df2=0;
    if (mV != 0.0) { f1 = (mA/mV) * a; df1 = (mA/mV) * da; }
    if (mB != 0.0) { f2 = (mA/mB) * b; df2 = (mA/mB) * db; }

    Fcos.push_back(f1);  dFcos.push_back(df1);
    Fcos2.push_back(f2); dFcos2.push_back(df2);
  }

  // Bottom row: pad 8 — x_B dependence of the scaled FUU ratios (no title, top-right legend)
  c->cd(7); gPad->Clear();
  c->cd(8);
  if (!xcenters.empty()) {
    auto gF1 = new TGraphErrors((int)xcenters.size());
    auto gF2 = new TGraphErrors((int)xcenters.size());
    for (int i=0;i<(int)xcenters.size();++i) {
      gF1->SetPoint(i, xcenters[i], Fcos[i]);  gF1->SetPointError(i, 0.0, dFcos[i]);
      gF2->SetPoint(i, xcenters[i], Fcos2[i]); gF2->SetPointError(i, 0.0, dFcos2[i]);
    }
    gF1->SetMarkerStyle(20); gF1->SetMarkerSize(1.0);
    gF2->SetMarkerStyle(21); gF2->SetMarkerSize(1.0);

    // No title, just axes
    gF1->SetTitle(";x_{B};Ratio");
    gF1->GetYaxis()->SetRangeUser(-1.0, 1.0);
    gF1->Draw("AP");
    gF2->Draw("P SAME");

    // Legend in top-right, inside pad
    auto legAB = new TLegend(0.62, 0.74, 0.94, 0.92);
    legAB->SetBorderSize(1); legAB->SetFillStyle(0); legAB->SetTextSize(0.035);
    legAB->AddEntry(gF1, "F_{UU}^{\\cos\\phi}/F_{UU}", "p");
    legAB->AddEntry(gF2, "F_{UU}^{\\cos2\\phi}/F_{UU}", "p");
    legAB->Draw();
  } else {
    TH1D* frame = new TH1D("frameAB",";x_{B};Ratio",10,0.06,0.60);
    frame->SetMinimum(-1.0); frame->SetMaximum(1.0);
    frame->Draw("AXIS");
    auto leg = new TLegend(0.62, 0.80, 0.94, 0.92);
    leg->SetBorderSize(1); leg->SetFillStyle(0); leg->SetTextSize(0.035);
    leg->AddEntry((TObject*)0, "No fit results collected", "");
    leg->Draw();
  }
  c->cd(9); gPad->Clear();

  ensure_dir("output/enpi+");
  const std::string outpdf = "output/enpi+/" + prop + "_unfolded.pdf";
  c->SaveAs(outpdf.c_str());
  std::cout << "Saved: " << outpdf << "\n";

  // Save arrays in the exact style you require (scaled values)
  save_arrays_style(prop, xcenters, Fcos, dFcos, Fcos2, dFcos2);
}

// ------------ Counters ------------
static void print_counters(const char* label, const Counters& C, const std::vector<std::string>& props) {
  std::cout << "\n[" << label << " counters]\n";
  std::cout << "  total read              : " << C.total << "\n";
  std::cout << "  pass common cuts        : " << C.pass_common_cnt << "\n";
  std::cout << "  pass fiducial (>=100)   : " << C.pass_fid_cnt << "\n";
  std::cout << "  fail |t|>1 or <0        : " << C.fail_abs_t << "\n";
  std::cout << "  fail Mx2 not (0.80,1.00): " << C.fail_mx2 << "\n";
  std::cout << "  fail no fid branch      : " << C.fail_no_fid_branch << "\n";
  std::cout << "  fail fid<100            : " << C.fail_fid_lt_100 << "\n";
  std::cout << "  fail x outside bins     : " << C.fail_xbin << "\n";
  for (const auto& p : props) {
    auto it = C.pass_prop_cnt.find(p);
    long long v = (it==C.pass_prop_cnt.end()? 0 : it->second);
    std::cout << "  pass property " << p << "     : " << v << "\n";
  }
}

// ------------ main ------------
int main(int argc, char** argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) return 1;

  std::unique_ptr<TFile> fD(TFile::Open(cfg.data_path.c_str(), "READ"));
  std::unique_ptr<TFile> fG(TFile::Open(cfg.gen_path.c_str(),  "READ"));
  std::unique_ptr<TFile> fR(TFile::Open(cfg.rec_path.c_str(),  "READ"));
  if (!fD || fD->IsZombie() || !fG || fG->IsZombie() || !fR || fR->IsZombie()) {
    std::cerr << "ERROR: could not open one or more input files.\n";
    return 2;
  }

  TTree* tD = (TTree*)fD->Get("PhysicsEvents");
  TTree* tG = (TTree*)fG->Get("PhysicsEvents");
  TTree* tR = (TTree*)fR->Get("PhysicsEvents");
  if (!tD || !tG || !tR) {
    std::cerr << "ERROR: missing PhysicsEvents tree in one or more files.\n";
    return 3;
  }

  if (cfg.list_branches) {
    list_tree_branches(tD, "DATA");
    list_tree_branches(tG, "GEN");
    list_tree_branches(tR, "REC");
  }

  auto H = make_hist_map(cfg.phi_nbins);
  auto Dep = make_dep_map();

  Counters cD, cG, cR;
  loop_tree_fill(tD, cfg.bn, cfg.apply_fid_cut, /*count_fid=*/true,  H, &Dep, cD, cfg.debug_print, "DATA");
  loop_tree_fill(tG, cfg.bn, /*apply_fid=*/false, /*count_fid=*/false, H, nullptr, cG, cfg.debug_print, "GEN");
  loop_tree_fill(tR, cfg.bn, cfg.apply_fid_cut, /*count_fid=*/true,  H, nullptr, cR, cfg.debug_print, "REC");

  const auto props = properties_to_run(cfg.which_property);

  // Totals sanity & detailed counters
  for (const auto& p : props) {
    const auto& hs = H.at(p);
    long long dsum=0, gsum=0, rsum=0;
    for (size_t i=0;i<hs.D.size();++i) {
      dsum += (long long)hs.D[i]->GetEntries();
      gsum += (long long)hs.G[i]->GetEntries();
      rsum += (long long)hs.R[i]->GetEntries();
    }
    std::cout << "[check] property " << p << " totals: "
              << "D=" << dsum << "  G=" << gsum << "  R=" << rsum << "\n";
  }
  print_counters("DATA", cD, props);
  print_counters("GEN",  cG, props);
  print_counters("REC",  cR, props);

  // Draw and save per property
  for (const auto& p : props) {
    draw_property_and_save(p, H, Dep);
  }

  return 0;
}