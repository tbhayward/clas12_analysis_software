// unfold_phi_unfold_fit.cpp
// Exclusive pi+ unfolding and cos(n*phi) fits (phi in degrees).
// Optional sine terms are toggled by --fit-sin (adds sin(phi) and sin(2phi)).
// Two canvases per property (enpi, enpiLowt, enpiMidt, enpiHight):
//   1) Main 2x3:
//        pads 1–4 : unfolded-yield points + fit for the four x_B bins
//                    (Y range now forced to [0, 2*max(points)]).
//                    Legend moved to top-right and made taller.
//        pad 5    : FUU ratios vs x_B (points only)
//        pad 6    : fitted amplitudes vs x_B (points only)
//   2) “Sines” 2x3:
//        pads 1–4 : DATA histograms of sin(phi) and sin(2phi) overlaid,
//                    legend bigger and pinned to top-right, shows <sinφ> and <sin2φ>.
//        pads 5–6 : left empty.
//
// Binning:
//   x_B edges: 0.10, 0.30, 0.40, 0.50, 0.60 (4 bins)
//   |t| windows:
//     enpi      : 0.10 ≤ |t| ≤ 1.20
//     enpiLowt  : 0.10 ≤ |t| ≤ 0.4667
//     enpiMidt  : 0.4667 ≤ |t| ≤ 0.8333
//     enpiHight : 0.8333 ≤ |t| ≤ 1.20
//
// Build:
//   g++ -O2 -std=c++17 unfold_phi_unfold_fit.cpp `root-config --cflags --libs` -o unfold_phi_unfold_fit
//
// Run:
//   ./unfold_phi_unfold_fit data.root gen.root rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]
//        [--x x] [--phi phi] [--t t] [--mx2 Mx2] [--fid fiducial_status]
//        [--DepA DepA] [--DepB DepB] [--DepV DepV] [--DepW DepW]
//        [--phibins 24] [--no-fid-cut] [--fit-sin] [--debug N] [--list-branches]
//
// Outputs to output/enpi+/:
//   <property>_unfolded.pdf
//   <property>_sines.pdf
//   <property>_unfolded_fit_arrays.txt
//     propAUUcosphi   = { {xB, val, err}, ... };
//     propAUUcos2phi  = { {xB, val, err}, ... };
//     (if --fit-sin) propAUUsinphi   = { {xB, val, err}, ... };
//     (if --fit-sin) propAUUsin2phi  = { {xB, val, err}, ... };

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
#include "TPave.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TLine.h"

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
  bool fit_sin        = false;   // --fit-sin turns on BOTH sin(phi) and sin(2phi)
  int  debug_print    = 0;
  bool list_branches  = false;
};

static void print_usage(const char* argv0) {
  std::cout << "Usage:\n  " << argv0
            << " data.root gen.root rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]\n"
            << "    [--x x] [--phi phi] [--t t] [--mx2 Mx2] [--fid fiducial_status]\n"
            << "    [--DepA DepA] [--DepB DepB] [--DepV DepV] [--DepW DepW]\n"
            << "    [--phibins 24] [--no-fid-cut] [--fit-sin] [--debug N] [--list-branches]\n";
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
    if (a == "--fit-sin")    { cfg.fit_sin = true; continue; }
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

// ---- x bins: 0.10, 0.30, 0.40, 0.50, 0.60 (4 bins)
static std::vector<double> x_edges() { return {0.10, 0.30, 0.40, 0.50, 0.60}; }
static int xbin_index(double x, const std::vector<double>& e) {
  for (size_t i=0;i+1<e.size();++i) if (x >= e[i] && x < e[i+1]) return (int)i;
  return -1;
}

// ---- |t| windows
static inline bool pass_prop_window(const std::string& prop, double t) {
  const double at = std::fabs(t);
  if (prop == "enpi")      return (at >= 0.10 && at <= 1.20);
  if (prop == "enpiLowt")  return (at >= 0.10 && at <= 0.4667);
  if (prop == "enpiMidt")  return (at >= 0.4667 && at <= 0.8333);
  if (prop == "enpiHight") return (at >= 0.8333 && at <= 1.20);
  return false;
}
static std::string prop_title(const std::string& p) {
  if (p == "enpi")      return "ep -> e' n #pi^{+}";
  if (p == "enpiLowt")  return "ep -> e' n #pi^{+}, low |t|";
  if (p == "enpiMidt")  return "ep -> e' n #pi^{+}, mid |t|";
  if (p == "enpiHight") return "ep -> e' n #pi^{+}, high |t|";
  return p;
}
static std::string prop_tlabel(const std::string& p) {
  if (p == "enpi")      return "-t #in [0.10, 1.20]";
  if (p == "enpiLowt")  return "-t #in [0.10, 0.4667]";
  if (p == "enpiMidt")  return "-t #in [0.4667, 0.8333]";
  if (p == "enpiHight") return "-t #in [0.8333, 1.20]";
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
    std::cerr << "ERROR: missing required branches on " << who << "\n";
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
  Long64_t fail_xbin=0;
  std::map<std::string, Long64_t> pass_prop_cnt;
};

static void debug_event_print(const char* label, int idx, const BranchHandles& b,
                              const std::vector<std::string>& props, const std::vector<double>& xe) {
  std::cout << "[" << label << " evt " << idx << "] "
            << "x_{B}=" << b.x << "  t=" << b.t << "  Mx2=" << b.mx2
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
  TTree* tr, const BranchNames& bn, bool apply_fid, bool /*count_fid*/,
  std::map<std::string,HSet>& H, std::map<std::string,DepMeans>* depPtr,
  Counters& C, int debugN, const char* dbg_label)
{
  if (!tr) return;
  BranchHandles b;
  if (!bind_tree(tr, bn, b, dbg_label)) {
    std::cerr << "ERROR: failed to bind branches on " << dbg_label << "\n";
    return;
  }
  const auto xe = x_edges();
  const auto props = all_props();
  const Long64_t nent = tr->GetEntries();
  const int debug_limit = std::max(0, debugN);

  for (Long64_t i=0;i<nent;i++) {
    tr->GetEntry(i);
    if (debug_limit > 0 && i < debug_limit) debug_event_print(dbg_label, (int)i, b, props, xe);

    // Common cuts: Mx2 and optional fiducial; (t selection handled per property later)
    const bool ok_mx2 = (b.mx2 > 0.80) && (b.mx2 < 1.00);
    if (!ok_mx2) continue;
    if (apply_fid) { if (!b.has_fid || b.fid < 100) continue; }
    C.pass_common_cnt++;

    int ib = xbin_index(b.x, xe);
    if (ib < 0) { C.fail_xbin++; continue; }

    double phideg = b.phi * 180.0 / TMath::Pi();
    phideg = std::fmod(phideg, 360.0); if (phideg < 0) phideg += 360.0;

    for (const auto& p : props) {
      if (!pass_prop_window(p, b.t)) continue;
      C.pass_prop_cnt[p]++;

      if (std::string(dbg_label) == "DATA") {
        H[p].D[ib]->Fill(phideg);
        if (depPtr) {
          auto& dep = depPtr->at(p);
          dep.sumA[ib] += b.DepA;
          dep.sumB[ib] += b.DepB;
          dep.sumV[ib] += b.DepV;
          dep.sumW[ib] += b.DepW;
          dep.count[ib] += 1;
        }
      } else if (std::string(dbg_label) == "GEN") {
        H[p].G[ib]->Fill(phideg);
      } else if (std::string(dbg_label) == "REC") {
        H[p].R[ib]->Fill(phideg);
      }
    }
  }
  C.total = tr->GetEntries();
}

// ------------ Fit ------------
struct FitResult {
  double C=0, dC=0;
  double A=0, dA=0;      // cosφ
  double B=0, dB=0;      // cos2φ
  double D=0, dD=0;      // sinφ   (only if --fit-sin)
  double E=0, dE=0;      // sin2φ  (only if --fit-sin)
  double chi2=0; int ndf=0;
  int npoints=0;
};

static FitResult make_unfold_graph_and_fit(
  TH1D* hD, TH1D* hG, TH1D* hR, TGraphErrors*& g, TF1*& ffit, bool fit_sin)
{
  FitResult fr;
  if (!hD || !hG || !hR) { g=nullptr; ffit=nullptr; return fr; }

  const int nb = hD->GetNbinsX();
  std::vector<double> X, Y, EY;
  X.reserve(nb); Y.reserve(nb); EY.reserve(nb);

  for (int b=1; b<=nb; ++b) {
    const double Dd = hD->GetBinContent(b);
    const double G  = hG->GetBinContent(b);
    const double R  = hR->GetBinContent(b);
    const double xc = hD->GetBinCenter(b); // degrees
    if (R <= 0.0 || G <= 0.0) continue;

    const double U = (Dd * G) / R;
    const double term1 = (G*G)/(R*R) * Dd;
    const double term2 = (Dd*Dd)/(R*R) * G;
    const double term3 = (Dd*Dd*G*G)/(R*R*R);
    const double sig = std::sqrt(std::max(0.0, term1 + term2 + term3));

    X.push_back(xc); Y.push_back(U); EY.push_back(sig);
  }

  fr.npoints = (int)X.size();
  if (fr.npoints == 0) { g=nullptr; ffit=nullptr; return fr; }

  g = new TGraphErrors(fr.npoints);
  for (int i=0;i<fr.npoints;i++) { g->SetPoint(i, X[i], Y[i]); g->SetPointError(i, 0.0, EY[i]); }
  g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->SetLineWidth(1);

  const double xmin = 0.0, xmax = 360.0;
  if (!fit_sin) {
    ffit = new TF1("fitCAB_deg",
                   "[0]*(1 + [1]*cos(TMath::Pi()/180.0*x) + [2]*cos(2.0*TMath::Pi()/180.0*x))",
                   xmin, xmax);
    ffit->SetParNames("C","A","B");
    ffit->SetParameters(1.0, 0.0, 0.0);
  } else {
    ffit = new TF1("fitCABDE_deg",
                   "[0]*(1 + [1]*cos(TMath::Pi()/180.0*x) + [2]*cos(2.0*TMath::Pi()/180.0*x) + [3]*sin(TMath::Pi()/180.0*x) + [4]*sin(2.0*TMath::Pi()/180.0*x))",
                   xmin, xmax);
    ffit->SetParNames("C","A","B","D","E");
    ffit->SetParameters(1.0, 0.0, 0.0, 0.0, 0.0);
  }
  // Reasonable bounds
  ffit->SetParLimits(0, 0.0, 1e12);
  for (int ip=1; ip<=4; ++ip) if (ip < ffit->GetNpar()) ffit->SetParLimits(ip, -2.0, 2.0);

  g->Fit(ffit, "Q"); // quiet

  fr.C  = ffit->GetParameter(0); fr.dC = ffit->GetParError(0);
  fr.A  = ffit->GetParameter(1); fr.dA = ffit->GetParError(1);
  fr.B  = ffit->GetParameter(2); fr.dB = ffit->GetParError(2);
  if (fit_sin) {
    fr.D = ffit->GetParameter(3); fr.dD = ffit->GetParError(3);
    fr.E = ffit->GetParameter(4); fr.dE = ffit->GetParError(4);
  }
  fr.chi2 = ffit->GetChisquare(); fr.ndf  = ffit->GetNDF();

  return fr;
}

// ---------- Legend helpers (white background, top-right, adaptive height) ----------
static void DrawFitLegendBoxTR(const FitResult& fr, bool fit_sin) {
  // Lines: A, B, plus D,E if fit_sin
  int nLines = 2 + (fit_sin ? 2 : 0);
  double x2=0.93, x1=0.55;
  double y2=0.92;
  double lineH = 0.08;               // taller lines than before
  double padH  = 0.03;
  double boxH  = padH + nLines*lineH + 0.04;
  double y1    = y2 - boxH;
  if (y1 < 0.55) { y1 = 0.55; }      // keep it inside
  TPave* box = new TPave(x1,y1,x2,y2,1,"NDC");
  box->SetFillStyle(1001); box->SetFillColor(kWhite);
  box->SetLineColor(kBlack); box->SetLineWidth(2);
  box->Draw("same");

  TLatex lat;
  lat.SetNDC(); lat.SetTextColor(kBlack); lat.SetTextSize(0.040); lat.SetTextAlign(13);

  double xtext = x1 + 0.03;
  double yline = y2 - 0.12;
  lat.DrawLatex(xtext, yline, Form("A_{UU}^{cos#phi} = %.3g  #pm %.3g",  fr.A, fr.dA));
  yline -= lineH;
  lat.DrawLatex(xtext, yline, Form("A_{UU}^{cos2#phi} = %.3g  #pm %.3g", fr.B, fr.dB));
  if (fit_sin) {
    yline -= lineH;
    lat.DrawLatex(xtext, yline, Form("A_{UU}^{sin#phi} = %.3g  #pm %.3g",  fr.D, fr.dD));
    yline -= lineH;
    lat.DrawLatex(xtext, yline, Form("A_{UU}^{sin2#phi} = %.3g  #pm %.3g", fr.E, fr.dE));
  }
}

static void DrawBottomLegendBox(double x1, double y1, double x2, double y2,
                                TGraphErrors* gF1, TGraphErrors* gF2,
                                bool fit_sin, TGraphErrors* gFs, TGraphErrors* gFs2)
{
  TPave* box = new TPave(x1,y1,x2,y2,1,"NDC");
  box->SetFillStyle(1001); box->SetFillColor(kWhite);
  box->SetLineColor(kBlack); box->SetLineWidth(2);
  box->Draw("same");

  TLatex lat; lat.SetNDC(); lat.SetTextColor(kBlack); lat.SetTextSize(0.034); lat.SetTextAlign(13);
  const double dx = 0.05*(x2-x1);
  const double mdx = 0.02*(x2-x1);
  double xtext = x1 + dx + mdx;
  double xmark = x1 + dx;
  double yline = y2 - 0.20*(y2-y1);

  if (gF1) {
    TMarker* m1 = new TMarker(xmark, yline, gF1->GetMarkerStyle()); m1->SetNDC(true);
    m1->SetMarkerColor(gF1->GetMarkerColor()); m1->SetMarkerSize(gF1->GetMarkerSize()); m1->Draw();
    lat.DrawLatex(xtext, yline, "F_{UU}^{cos#phi}/F_{UU}");
  }
  yline -= 0.25*(y2-y1);
  if (gF2) {
    TMarker* m2 = new TMarker(xmark, yline, gF2->GetMarkerStyle()); m2->SetNDC(true);
    m2->SetMarkerColor(gF2->GetMarkerColor()); m2->SetMarkerSize(gF2->GetMarkerSize()); m2->Draw();
    lat.DrawLatex(xtext, yline, "F_{UU}^{cos2#phi}/F_{UU}");
  }
  if (fit_sin && gFs) {
    yline -= 0.25*(y2-y1);
    TMarker* m3 = new TMarker(xmark, yline, gFs->GetMarkerStyle()); m3->SetNDC(true);
    m3->SetMarkerColor(gFs->GetMarkerColor()); m3->SetMarkerSize(gFs->GetMarkerSize()); m3->Draw();
    lat.DrawLatex(xtext, yline, "F_{UU}^{sin#phi}/F_{UU}");
  }
  if (fit_sin && gFs2) {
    yline -= 0.25*(y2-y1);
    TMarker* m4 = new TMarker(xmark, yline, gFs2->GetMarkerStyle()); m4->SetNDC(true);
    m4->SetMarkerColor(gFs2->GetMarkerColor()); m4->SetMarkerSize(gFs2->GetMarkerSize()); m4->Draw();
    lat.DrawLatex(xtext, yline, "F_{UU}^{sin2#phi}/F_{UU}");
  }
}

// ------------ Save text arrays ------------
static void save_arrays_style(
  const std::string& prop,
  const std::vector<double>& xcenters,
  const std::vector<double>& Fcos, const std::vector<double>& dFcos,
  const std::vector<double>& Fcos2, const std::vector<double>& dFcos2,
  bool fit_sin,
  const std::vector<double>& Fsin, const std::vector<double>& dFsin,
  const std::vector<double>& Fsin2, const std::vector<double>& dFsin2) {
  ensure_dir("output/enpi+");
  std::string out = "output/enpi+/" + prop + "_unfolded_fit_arrays.txt";
  std::ofstream ofs(out);
  ofs << std::fixed << std::setprecision(9);

  ofs << prop << "AUUcosphi = {";
  for (size_t i=0;i<xcenters.size();++i) {
    ofs << "{" << xcenters[i] << ", " << Fcos[i] << ", " << dFcos[i] << "}";
    if (i + 1 < xcenters.size()) ofs << ", ";
  }
  ofs << "};\n";

  ofs << prop << "AUUcos2phi = {";
  for (size_t i=0;i<xcenters.size();++i) {
    ofs << "{" << xcenters[i] << ", " << Fcos2[i] << ", " << dFcos2[i] << "}";
    if (i + 1 < xcenters.size()) ofs << ", ";
  }
  ofs << "};\n";

  if (fit_sin) {
    ofs << prop << "AUUsinphi = {";
    for (size_t i=0;i<xcenters.size();++i) {
      ofs << "{" << xcenters[i] << ", " << Fsin[i] << ", " << dFsin[i] << "}";
      if (i + 1 < xcenters.size()) ofs << ", ";
    }
    ofs << "};\n";

    ofs << prop << "AUUsin2phi = {";
    for (size_t i=0;i<xcenters.size();++i) {
      ofs << "{" << xcenters[i] << ", " << Fsin2[i] << ", " << dFsin2[i] << "}";
      if (i + 1 < xcenters.size()) ofs << ", ";
    }
    ofs << "};\n";
  }

  ofs.close();
  std::cout << "Wrote arrays: " << out << "\n";
}

// ---- Sines canvas helpers ----
struct SineHists {
  std::unique_ptr<TH1D> hSin;
  std::unique_ptr<TH1D> hSin2;
  double meanSin=0, errSin=0;
  double meanSin2=0, errSin2=0;
};
static std::map<std::string, std::vector<SineHists>>
build_sine_hists(TTree* tD, const BranchNames& bn, bool apply_fid) {
  std::map<std::string, std::vector<SineHists>> out;
  if (!tD) return out;
  BranchHandles b;
  if (!bind_tree(tD, bn, b, "DATA")) return out;
  const auto props = all_props();
  const auto xe = x_edges();
  const int NX = (int)xe.size()-1;

  for (const auto& p : props) {
    std::vector<SineHists> v(NX);
    for (int i=0;i<NX;i++) {
      v[i].hSin  = std::make_unique<TH1D>(Form("hSin_%s_x%d",  p.c_str(), i),  ";sin#phi;Events",  100, -1.0, 1.0);
      v[i].hSin2 = std::make_unique<TH1D>(Form("hSin2_%s_x%d", p.c_str(), i),  ";sin2#phi;Events", 100, -1.0, 1.0);
    }
    out.emplace(p, std::move(v));
  }

  const Long64_t nent = tD->GetEntries();
  std::vector<double> sum1(props.size()*NX,0.0), sum1sq(props.size()*NX,0.0);
  std::vector<double> sum2(props.size()*NX,0.0), sum2sq(props.size()*NX,0.0);
  std::vector<long long> cnt(props.size()*NX,0);

  for (Long64_t i=0;i<nent;i++) {
    tD->GetEntry(i);
    const bool ok_mx2 = (b.mx2 > 0.80) && (b.mx2 < 1.00);
    if (!ok_mx2) continue;
    if (apply_fid) { if (!b.has_fid || b.fid < 100) continue; }

    int ib = xbin_index(b.x, xe);
    if (ib < 0) continue;

    const double s1 = std::sin(b.phi);
    const double s2 = std::sin(2.0*b.phi);

    for (size_t ip=0; ip<props.size(); ++ip) {
      const auto& p = props[ip];
      if (!pass_prop_window(p, b.t)) continue;
      auto& rec = out[p][ib];
      rec.hSin->Fill(s1);
      rec.hSin2->Fill(s2);
      const size_t idx = ip*NX + ib;
      sum1[idx]  += s1; sum1sq[idx] += s1*s1;
      sum2[idx]  += s2; sum2sq[idx] += s2*s2;
      cnt[idx]   += 1;
    }
  }

  for (size_t ip=0; ip<props.size(); ++ip) {
    const auto& p = props[ip];
    for (int ib=0; ib<NX; ++ib) {
      const size_t idx = ip*NX + ib;
      auto& rec = out[p][ib];
      const double N = (double)cnt[idx];
      if (N > 0) {
        rec.meanSin  = sum1[idx]/N;
        rec.meanSin2 = sum2[idx]/N;
        const double var1 = std::max(0.0, sum1sq[idx]/N - rec.meanSin*rec.meanSin);
        const double var2 = std::max(0.0, sum2sq[idx]/N - rec.meanSin2*rec.meanSin2);
        rec.errSin  = std::sqrt(var1 / N);
        rec.errSin2 = std::sqrt(var2 / N);
      }
    }
  }
  return out;
}

static void draw_property_and_save(
  const std::string& prop,
  const std::map<std::string,HSet>& H,
  const std::map<std::string,DepMeans>& depMap,
  bool fit_sin)
{
  const auto xe = x_edges();
  const int NX = (int)xe.size()-1;

  // Global style
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  TCanvas* c = new TCanvas(("c_"+prop).c_str(), ("Unfolded phi fits: "+prop).c_str(), 1500, 1000);
  c->Divide(3,2); // 2 rows x 3 cols

  std::vector<double> xcenters;
  std::vector<double> Avec, dAvec, Bvec, dBvec, Dvec, dDvec, Evec, dEvec;
  std::vector<double> meanDepA(NX,0.0), meanDepB(NX,0.0), meanDepV(NX,0.0), meanDepW(NX,0.0);
  const auto& dep = depMap.at(prop);

  // Pads 1–4: unfolded φ and fits per x_B bin
  for (int ib=0; ib<NX; ++ib) {
    c->cd(ib+1);
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.06);
    gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.12);

    TH1D* hD = H.at(prop).D[ib].get();
    TH1D* hG = H.at(prop).G[ib].get();
    TH1D* hR = H.at(prop).R[ib].get();

    TGraphErrors* g = nullptr; TF1* fit = nullptr;
    FitResult fr = make_unfold_graph_and_fit(hD,hG,hR,g,fit,fit_sin);

    const double xl = xe[ib], xh = xe[ib+1];
    std::string title = prop_title(prop) + Form(", x_{B} #in [%.2f, %.2f), %s", xl, xh, prop_tlabel(prop).c_str());

    if (dep.count[ib] > 0) {
      meanDepA[ib] = dep.sumA[ib] / (double)dep.count[ib];
      meanDepB[ib] = dep.sumB[ib] / (double)dep.count[ib];
      meanDepV[ib] = dep.sumV[ib] / (double)dep.count[ib];
      meanDepW[ib] = dep.sumW[ib] / (double)dep.count[ib];
    }

    if (g) {
      g->SetTitle((title + ";#phi (deg);Unfolded yield").c_str());

      // --- NEW: force y-range to [0, 2*max(points)] to leave room for tall legend
      double ymax = -1e300, xx, yy;
      for (int ip=0; ip<g->GetN(); ++ip) { g->GetPoint(ip, xx, yy); ymax = std::max(ymax, yy); }
      if (!(ymax > 0)) ymax = 1.0;
      g->GetYaxis()->SetRangeUser(0.0, 2.0*ymax);

      g->Draw("AP"); if (fit) fit->Draw("LSAME");

      // Tall, top-right legend
      DrawFitLegendBoxTR(fr, fit_sin);
      gPad->Modified(); gPad->Update();

      xcenters.push_back(0.5*(xl+xh));
      Avec.push_back(fr.A); dAvec.push_back(fr.dA);
      Bvec.push_back(fr.B); dBvec.push_back(fr.dB);
      if (fit_sin) { Dvec.push_back(fr.D); dDvec.push_back(fr.dD); Evec.push_back(fr.E); dEvec.push_back(fr.dE); }
    } else {
      TH1D* frame = new TH1D("frame","",10,0,360);
      frame->SetTitle((title + ";#phi (deg);Unfolded yield").c_str());
      frame->SetMinimum(0); frame->SetMaximum(1);
      frame->Draw("AXIS");
      DrawFitLegendBoxTR(FitResult{}, fit_sin);
      gPad->Modified(); gPad->Update();
    }
  }

  // Scale to FUU ratios using depolarization means.
  // Mapping: cosφ -> A/V, cos2φ -> A/B, sinφ -> A/W, sin2φ -> A/B (proxy; update if a dedicated dep is added)
  std::vector<double> Fcos, dFcos, Fcos2, dFcos2, Fsin, dFsin, Fsin2, dFsin2;
  for (size_t i=0;i<xcenters.size();++i) {
    const int ib = (int)i;
    const double mA = meanDepA[ib], mB = meanDepB[ib], mV = meanDepV[ib], mW = meanDepW[ib];

    double f1=0, df1=0; if (mV!=0.0) { f1=(mA/mV)*Avec[i]; df1=(mA/mV)*dAvec[i]; }
    Fcos.push_back(f1); dFcos.push_back(df1);

    double f2=0, df2=0; if (mB!=0.0) { f2=(mA/mB)*Bvec[i]; df2=(mA/mB)*dBvec[i]; }
    Fcos2.push_back(f2); dFcos2.push_back(df2);

    if (fit_sin) {
      double fs=0, dfs=0; if (mW!=0.0) { fs=(mA/mW)*Dvec[i]; dfs=(mA/mW)*dDvec[i]; }
      Fsin.push_back(fs); dFsin.push_back(dfs);

      // No dedicated dep branch for sin2ϕ in current inputs; use B as a proxy.
      double fs2=0, dfs2=0; if (mB!=0.0) { fs2=(mA/mB)*Evec[i]; dfs2=(mA/mB)*dEvec[i]; }
      Fsin2.push_back(fs2); dFsin2.push_back(dfs2);
    }
  }

  // Pad 5: FUU ratios vs x_B (points only)
  c->cd(5);
  gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.06);
  gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.12);

  TGraphErrors *gF1=nullptr, *gF2=nullptr, *gFs=nullptr, *gFs2=nullptr;
  if (!xcenters.empty()) {
    const int N = (int)xcenters.size();
    gF1 = new TGraphErrors(N);
    gF2 = new TGraphErrors(N);
    if (fit_sin) { gFs = new TGraphErrors(N); gFs2 = new TGraphErrors(N); }

    for (int i=0;i<N;++i) {
      gF1->SetPoint(i, xcenters[i] - 0.004, Fcos[i]);   gF1->SetPointError(i, 0.0, dFcos[i]);
      gF2->SetPoint(i, xcenters[i] + 0.004, Fcos2[i]);  gF2->SetPointError(i, 0.0, dFcos2[i]);
      if (fit_sin) {
        gFs ->SetPoint(i, xcenters[i], Fsin[i]);   gFs ->SetPointError(i, 0.0, dFsin[i]);
        gFs2->SetPoint(i, xcenters[i], Fsin2[i]);  gFs2->SetPointError(i, 0.0, dFsin2[i]);
      }
    }

    gF1->SetMarkerStyle(20); gF1->SetMarkerSize(1.0); gF1->SetMarkerColor(kRed);
    gF2->SetMarkerStyle(21); gF2->SetMarkerSize(1.0); gF2->SetMarkerColor(kBlue);
    if (fit_sin) {
      gFs ->SetMarkerStyle(22); gFs ->SetMarkerSize(1.0); gFs ->SetMarkerColor(kGreen+2);
      gFs2->SetMarkerStyle(23); gFs2->SetMarkerSize(1.0); gFs2->SetMarkerColor(kMagenta+2);
    }

    gF1->SetTitle(";x_{B};Ratio");
    gF1->GetYaxis()->SetRangeUser(-1.0, 1.0);
    gF1->Draw("AP");
    gF1->GetXaxis()->SetLimits(0.10, 0.60);
    gF2->Draw("P SAME");
    if (fit_sin) { gFs->Draw("P SAME"); gFs2->Draw("P SAME"); }

    TLine* line0_mid = new TLine(0.10, 0.0, 0.60, 0.0);
    line0_mid->SetLineStyle(2); line0_mid->SetLineWidth(1); line0_mid->SetLineColor(kBlack);
    line0_mid->Draw("SAME");

    if (fit_sin) {
      DrawBottomLegendBox(0.58, 0.62, 0.94, 0.90, gF1, gF2, true, gFs, gFs2);
    } else {
      DrawBottomLegendBox(0.64, 0.70, 0.94, 0.88, gF1, gF2, false, nullptr, nullptr);
    }
    gPad->Modified(); gPad->Update();
  } else {
    TH1D* frame = new TH1D("frameAB",";x_{B};Ratio",10,0.10,0.60);
    frame->SetMinimum(-1.0); frame->SetMaximum(1.0);
    frame->Draw("AXIS");
  }

  // Pad 6: Amplitudes vs x_B (points only)
  c->cd(6);
  gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.06);
  gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.12);

  if (!xcenters.empty()) {
    const int N = (int)xcenters.size();
    auto gA = new TGraphErrors(N);
    auto gB = new TGraphErrors(N);
    TGraphErrors* gD=nullptr; TGraphErrors* gE=nullptr;

    for (int i=0;i<N;++i) {
      gA->SetPoint(i, xcenters[i] - 0.004, Avec[i]); gA->SetPointError(i, 0.0, dAvec[i]);
      gB->SetPoint(i, xcenters[i] + 0.004, Bvec[i]); gB->SetPointError(i, 0.0, dBvec[i]);
      if (fit_sin) {
        if (!gD) gD = new TGraphErrors(N);
        if (!gE) gE = new TGraphErrors(N);
        gD->SetPoint(i, xcenters[i], Dvec[i]); gD->SetPointError(i, 0.0, dDvec[i]);
        gE->SetPoint(i, xcenters[i], Evec[i]); gE->SetPointError(i, 0.0, dEvec[i]);
      }
    }
    gA->SetMarkerColor(kRed);   gA->SetMarkerStyle(20); gA->SetMarkerSize(1.0);
    gB->SetMarkerColor(kBlue);  gB->SetMarkerStyle(21); gB->SetMarkerSize(1.0);
    if (fit_sin) { gD->SetMarkerColor(kGreen+2); gD->SetMarkerStyle(22); gD->SetMarkerSize(1.0);
                   gE->SetMarkerColor(kMagenta+2); gE->SetMarkerStyle(23); gE->SetMarkerSize(1.0); }

    gA->SetTitle(";x_{B};Amplitude");
    gA->GetYaxis()->SetRangeUser(-1.0, 1.0);
    gA->Draw("AP");
    gA->GetXaxis()->SetLimits(0.10, 0.60);
    gB->Draw("P SAME");
    if (fit_sin) { gD->Draw("P SAME"); gE->Draw("P SAME"); }

    TLine* line0_amp = new TLine(0.10, 0.0, 0.60, 0.0);
    line0_amp->SetLineStyle(2); line0_amp->SetLineWidth(1); line0_amp->SetLineColor(kBlack);
    line0_amp->Draw("SAME");

    auto legAmp = new TLegend(0.55,0.66,0.94,0.90);
    legAmp->SetFillStyle(1001); legAmp->SetFillColor(kWhite);
    legAmp->SetBorderSize(1); legAmp->SetTextSize(0.034);
    legAmp->AddEntry(gA,"A_{UU}^{cos#phi}","p");
    legAmp->AddEntry(gB,"A_{UU}^{cos2#phi}","p");
    if (fit_sin) { legAmp->AddEntry(gD,"A_{UU}^{sin#phi}","p"); legAmp->AddEntry(gE,"A_{UU}^{sin2#phi}","p"); }
    legAmp->Draw();
  } else {
    TH1D* frame = new TH1D("frameAmp",";x_{B};Amplitude",10,0.10,0.60);
    frame->SetMinimum(-1.0); frame->SetMaximum(1.0);
    frame->Draw("AXIS");
  }

  ensure_dir("output/enpi+");
  const std::string outpdf = "output/enpi+/" + prop + "_unfolded.pdf";
  c->SaveAs(outpdf.c_str());
  std::cout << "Saved: " << outpdf << "\n";

  // Save arrays
  save_arrays_style(prop, xcenters, Fcos, dFcos, Fcos2, dFcos2, fit_sin, Fsin, dFsin, Fsin2, dFsin2);
}

// ---- Draw the “sines” canvas (DATA-only histograms of sinφ and sin2φ) ----
static void draw_sines_canvas(
  const std::string& prop,
  const std::map<std::string, std::vector<SineHists>>& S)
{
  const auto xe = x_edges();
  const int NX = (int)xe.size()-1;

  TCanvas* c = new TCanvas(("cSines_"+prop).c_str(), ("Sine distributions: "+prop).c_str(), 1500, 1000);
  c->Divide(3,2);

  for (int ib=0; ib<NX; ++ib) {
    c->cd(ib+1);
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.06);
    gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.12);

    auto& sh = S.at(prop)[ib];
    // Use same bins & titles already set. Style:
    sh.hSin ->SetLineColor(kBlue+1);  sh.hSin ->SetLineWidth(3);
    sh.hSin2->SetLineColor(kRed+1);   sh.hSin2->SetLineWidth(3);

    // Auto y max across both and give headroom
    double ymax = std::max(sh.hSin->GetMaximum(), sh.hSin2->GetMaximum());
    if (!(ymax > 0)) ymax = 1.0;
    sh.hSin->SetMaximum(1.15*ymax);
    sh.hSin2->SetMaximum(1.15*ymax);

    const double xl = xe[ib], xh = xe[ib+1];
    sh.hSin->SetTitle((prop_title(prop) + Form(", x_{B} #in [%.2f, %.2f), %s", xl, xh, prop_tlabel(prop).c_str())).c_str());
    sh.hSin->GetYaxis()->SetTitle("Events");

    sh.hSin->Draw("HIST");
    sh.hSin2->Draw("HIST SAME");

    // Bigger legend, top-right, with means
    double x2=0.90, x1=0.58, y2=0.89, y1=0.68; // bigger box
    TPave* box = new TPave(x1,y1,x2,y2,1,"NDC");
    box->SetFillStyle(1001); box->SetFillColor(kWhite);
    box->SetLineColor(kBlack); box->SetLineWidth(2);
    box->Draw("same");

    TLatex lat;
    lat.SetNDC(); lat.SetTextColor(kBlack); lat.SetTextSize(0.042); lat.SetTextAlign(13);
    double xtext = x1 + 0.03;
    double yline = y2 - 0.10;
    lat.SetTextColor(kBlue+1);
    lat.DrawLatex(xtext, yline, Form("#LTsin#phi#GT = %.4f  #pm %.4f",  sh.meanSin,  sh.errSin));
    yline -= 0.18;
    lat.SetTextColor(kRed+1);
    lat.DrawLatex(xtext, yline, Form("#LTsin2#phi#GT = %.4f  #pm %.4f", sh.meanSin2, sh.errSin2));
  }

  ensure_dir("output/enpi+");
  const std::string outpdf = "output/enpi+/" + prop + "_sines.pdf";
  c->SaveAs(outpdf.c_str());
  std::cout << "Saved: " << outpdf << "\n";
}

// ------------ Counters & main ------------
static void print_counters(const char* label, const Counters& C, const std::vector<std::string>& props) {
  std::cout << "\n[" << label << " counters]\n";
  std::cout << "  total read              : " << C.total << "\n";
  std::cout << "  pass common+fid cuts    : " << C.pass_common_cnt << "\n";
  std::cout << "  fail x outside bins     : " << C.fail_xbin << "\n";
  for (const auto& p : props) {
    auto it = C.pass_prop_cnt.find(p);
    long long v = (it==C.pass_prop_cnt.end()? 0 : it->second);
    std::cout << "  pass property " << p << "     : " << v << "\n";
  }
}

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

  auto H   = make_hist_map(cfg.phi_nbins);
  auto Dep = make_dep_map();

  Counters cD, cG, cR;
  loop_tree_fill(tD, cfg.bn, cfg.apply_fid_cut, /*count_fid=*/true,  H, &Dep, cD, cfg.debug_print, "DATA");
  loop_tree_fill(tG, cfg.bn, /*apply_fid=*/false, /*count_fid=*/false, H, nullptr, cG, cfg.debug_print, "GEN");
  loop_tree_fill(tR, cfg.bn, cfg.apply_fid_cut, /*count_fid=*/true,  H, nullptr, cR, cfg.debug_print, "REC");

  const auto props = properties_to_run(cfg.which_property);

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

  // Build sines once from DATA for the “sines” canvases
  auto SH = build_sine_hists(tD, cfg.bn, cfg.apply_fid_cut);

  for (const auto& p : props) {
    draw_property_and_save(p, H, Dep, cfg.fit_sin);
    draw_sines_canvas(p, SH);
  }

  return 0;
}