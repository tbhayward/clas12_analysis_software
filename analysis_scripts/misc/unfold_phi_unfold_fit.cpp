// unfold_phi_unfold_fit.cpp
// Exclusive pi+ unfolding and cos(n*phi) fits (phi in degrees).
// --fit-sin toggles *both* sin(phi) and sin(2phi) in the fit and in the
//   bottom-row summary plots.
//
// Layout per property (2x3 canvas):
//   [1] xB bin 0 unfolded (points) + CAB[+D+E] fit (red)
//   [2] xB bin 1 unfolded + fit
//   [3] xB bin 2 unfolded + fit
//   [4] xB bin 3 unfolded + fit
//   [5] (BLANK by request — structure-function ratio plot removed)
//   [6] Amplitudes vs xB (points only), Y in [-0.5, 0.5], horizontal legend
//
// A second 2x3 canvas per property shows, for each xB bin (first 4 pads),
//   two overlaid histograms: sin(phi) (blue) and sin(2phi) (red).
//   The legend (top-right, inside) lists <sin(phi)> and <sin(2phi)> with
//   standard errors. Legend text is normal-weight and small; legend widened
//   back out and aligned to the top-right corner of the subplot. Legend border
//   thickness matches the subplot frame thickness.
//
// Per–xB-bin unfolded-yield pads auto-scale to [0, 2*max(points)] to leave room
// for legends. Those legends are drawn as white TPaves fully inside each pad,
// horizontally anchored to the upper-right and with the text block vertically
// centered (nudged slightly upward so it never hugs the lower edge).
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
// Text output (to output/enpi+/):
//   <property>_unfolded_fit_arrays.txt
//     propAUUcosphi  = { {xB, val, err}, ... };
//     propAUUcos2phi = { {xB, val, err}, ... };
//     propAUUsinphi  = { ... }             (only if --fit-sin)
//     propAUUsin2phi = { ... }             (only if --fit-sin)

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

// ------------------ CLI/Config ------------------

struct BranchNames {
  std::string x   = "x";
  std::string phi = "phi";   // radians
  std::string t   = "t";     // signed
  std::string mx2 = "Mx2";
  std::string fid = "fiducial_status";
  // Depolarization (on DATA)
  std::string DepA = "DepA";
  std::string DepB = "DepB";
  std::string DepV = "DepV";
  std::string DepW = "DepW";
};

struct Config {
  std::string data_path, gen_path, rec_path;
  std::string which_property = "all";
  BranchNames bn;
  int  phi_nbins      = 24;
  bool apply_fid_cut  = true;
  bool fit_sin        = false; // toggles sinφ & sin2φ
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

// ------------------ binning & helper text ------------------

static std::vector<std::string> all_props() {
  return {"enpi", "enpiLowt", "enpiMidt", "enpiHight"};
}
static std::vector<std::string> properties_to_run(const std::string& which) {
  if (which == "all") return all_props();
  return {which};
}

// four xB bins
static std::vector<double> x_edges() { return {0.10, 0.30, 0.40, 0.50, 0.60}; }
static int xbin_index(double x, const std::vector<double>& e) {
  for (size_t i=0;i+1<e.size();++i) if (x >= e[i] && x < e[i+1]) return (int)i;
  return -1;
}

// absolute-t windows per property
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
static std::string prop_tlabel(const std::string& prop) {
  if (prop == "enpi")      return "-t #in [0.10, 1.20]";
  if (prop == "enpiLowt")  return "-t #in [0.10, 0.4667]";
  if (prop == "enpiMidt")  return "-t #in [0.4667, 0.8333]";
  if (prop == "enpiHight") return "-t #in [0.8333, 1.20]";
  return "";
}
static void ensure_dir(const std::string& d) {
  if (gSystem->AccessPathName(d.c_str())) gSystem->mkdir(d.c_str(), kTRUE);
}

// ------------------ binding ------------------

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

// ------------------ hist sets & dep means ------------------

struct HSet {
  std::vector<std::unique_ptr<TH1D>> D; // data (phi deg)
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
  const int NX = (int)x_edges().size()-1;
  for (const auto& p : all_props()) {
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

// ------------------ counters/debug (minimal) ------------------

struct Counters {
  Long64_t total=0, pass_common_cnt=0, fail_xbin=0;
  std::map<std::string, Long64_t> pass_prop_cnt;
};
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

// ------------------ event loop to fill histograms ------------------

static void loop_tree_fill(
  TTree* tr, const BranchNames& bn, bool apply_fid,
  std::map<std::string,HSet>& H,
  std::map<std::string,DepMeans>* depPtr,
  Counters& C, int /*debugN*/, const char* dbg_label)
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
  C.total = nent;

  for (Long64_t i=0;i<nent;i++) {
    tr->GetEntry(i);

    // common cuts
    const bool ok_abs_t = (std::fabs(b.t) >= 0.10 && std::fabs(b.t) <= 1.20);
    const bool ok_mx2   = (b.mx2 > 0.80) && (b.mx2 < 1.00);
    if (!(ok_abs_t && ok_mx2)) continue;
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
}

// ------------------ fitting ------------------

struct FitResult {
  double C=0, dC=0;   // overall scale
  double A=0, dA=0;   // cosφ
  double B=0, dB=0;   // cos2φ
  double D=0, dD=0;   // sinφ (if --fit-sin)
  double E=0, dE=0;   // sin2φ (if --fit-sin)
  double chi2=0; int ndf=0; int npoints=0;
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
  } else {
    ffit = new TF1("fitCABDE_deg",
                   "[0]*(1 + [1]*cos(TMath::Pi()/180.0*x) + [2]*cos(2.0*TMath::Pi()/180.0*x) + [3]*sin(TMath::Pi()/180.0*x) + [4]*sin(2.0*TMath::Pi()/180.0*x))",
                   xmin, xmax);
    ffit->SetParNames("C","A","B","D","E");
  }
  double ymean = 0.0; for (double v : Y) ymean += v; ymean = (fr.npoints>0 ? ymean/fr.npoints : 1.0);
  if (ymean<=0) ymean = 1.0;
  if (!fit_sin) {
    ffit->SetParameters(ymean, 0.0, 0.0);
  } else {
    ffit->SetParameters(ymean, 0.0, 0.0, 0.0, 0.0);
  }
  ffit->SetParLimits(0, 0.0, 1e12);
  for (int ip=1; ip<= (fit_sin?4:2); ++ip) ffit->SetParLimits(ip, -2.0, 2.0);

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

// ------------- legend helpers -------------

// Unfolded-yield pad legend (white box in upper-right, widened again, border = frame thickness)
static void DrawUnfoldLegendBox_UR(double x1, double y1, double x2, double y2,
                                   const FitResult& fr, bool fit_sin) {
  TPave* box = new TPave(x1,y1,x2,y2,1,"NDC");
  box->SetFillStyle(1001); box->SetFillColor(kWhite);
  box->SetLineColor(kBlack); box->SetLineWidth(1); // match frame thickness
  box->Draw("same");

  // Build text lines (omit pure AUU term intentionally)
  std::vector<std::string> lines;
  {
    std::ostringstream s1; s1<<std::fixed<<std::setprecision(4)
      << "A_{UU}^{cos#phi}  =  " << std::setw(7) << fr.A << "   #pm  " << std::setw(7) << fr.dA;
    std::ostringstream s2; s2<<std::fixed<<std::setprecision(4)
      << "A_{UU}^{cos2#phi} =  " << std::setw(7) << fr.B << "   #pm  " << std::setw(7) << fr.dB;
    lines.push_back(s1.str()); lines.push_back(s2.str());
    if (fit_sin) {
      std::ostringstream s3; s3<<std::fixed<<std::setprecision(4)
        << "A_{UU}^{sin#phi}  =  " << std::setw(7) << fr.D << "   #pm  " << std::setw(7) << fr.dD;
      std::ostringstream s4; s4<<std::fixed<<std::setprecision(4)
        << "A_{UU}^{sin2#phi} =  " << std::setw(7) << fr.E << "   #pm  " << std::setw(7) << fr.dE;
      lines.push_back(s3.str()); lines.push_back(s4.str());
    }
  }

  TLatex lat; lat.SetNDC(); lat.SetTextColor(kBlack); lat.SetTextSize(0.032); lat.SetTextAlign(13);

  const double boxH = (y2-y1);
  const double lineStep = 0.23*boxH;              // spacing per line
  const double textBlockH = lineStep * (double)lines.size();
  // Vertical centering, nudged a bit upward to avoid lower-edge crowding
  double ystart = y1 + 0.5*(boxH - textBlockH) + 0.2*boxH;
  double xtext  = x1 + 0.05*(x2-x1);

  for (size_t i=0;i<lines.size();++i) {
    lat.DrawLatex(xtext, ystart + (double)(lines.size()-1-i)*lineStep, lines[i].c_str());
  }
}

// Horizontal legend for bottom amplitude pad (kept NARROWER and placed to the right)
static void DrawHorizontalLegendAmplitude(bool haveSin, TGraphErrors* gCos, TGraphErrors* gCos2,
                                          TGraphErrors* gSin, TGraphErrors* gSin2) {
  auto L = new TLegend(0.50, 0.78, 0.94, 0.92); // tuned narrow legend (these were "perfect")
  L->SetNColumns(haveSin ? 4 : 2);
  L->SetBorderSize(1);
  L->SetFillStyle(1001);
  L->SetFillColor(kWhite);
  L->SetTextSize(0.030);

  L->AddEntry(gCos,  "A_{UU}^{cos#phi}", "p");
  L->AddEntry(gCos2, "A_{UU}^{cos2#phi}","p");
  if (haveSin) {
    L->AddEntry(gSin,  "A_{UU}^{sin#phi}", "p");
    L->AddEntry(gSin2, "A_{UU}^{sin2#phi}","p");
  }
  L->Draw();
}

// ------------------ save arrays ------------------

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

// ------------------ draw one property ------------------

static void draw_property_and_save(
  const std::string& prop,
  const std::map<std::string,HSet>& H,
  const std::map<std::string,DepMeans>& depMap,
  bool fit_sin)
{
  const auto xe = x_edges();
  const int NX = (int)xe.size()-1;

  // Global style & margins
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.22);   // extra padding so Y labels never clip
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetTitleSize(0.050, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  TCanvas* c = new TCanvas(("c_"+prop).c_str(), ("Unfolded phi fits: "+prop).c_str(), 1500, 1000);
  c->Divide(3,2); // 2 rows x 3 cols (pads 1..6)

  std::vector<double> xcenters;
  std::vector<double> Avec, dAvec, Bvec, dBvec, Dvec, dDvec, Evec, dEvec;
  std::vector<double> meanDepA(NX,0.0), meanDepB(NX,0.0), meanDepV(NX,0.0), meanDepW(NX,0.0);

  const auto& dep = depMap.at(prop);

  // First four pads: per-xB-bin unfolded φ and fit
  for (int ib=0; ib<NX; ++ib) {
    c->cd(ib+1);
    gPad->SetLeftMargin(0.22);
    gPad->SetRightMargin(0.06);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.12);

    TH1D* hD = H.at(prop).D[ib].get();
    TH1D* hG = H.at(prop).G[ib].get();
    TH1D* hR = H.at(prop).R[ib].get();

    TGraphErrors* g = nullptr; TF1* fit = nullptr;
    FitResult fr = make_unfold_graph_and_fit(hD,hG,hR,g,fit,fit_sin);

    const double xl = xe[ib], xh = xe[ib+1];
    std::string ttl = prop_title(prop) + Form(", x_{B} #in [%.2f, %.2f),  %s", xl, xh, prop_tlabel(prop).c_str());

    if (dep.count[ib] > 0) {
      meanDepA[ib] = dep.sumA[ib] / (double)dep.count[ib];
      meanDepB[ib] = dep.sumB[ib] / (double)dep.count[ib];
      meanDepV[ib] = dep.sumV[ib] / (double)dep.count[ib];
      meanDepW[ib] = dep.sumW[ib] / (double)dep.count[ib];
    }

    if (g) {
      g->SetTitle((ttl + ";#phi (deg);Unfolded yield").c_str());

      // y-range from points: [0, 2.0*max]
      double ymin=1e300,ymax=-1e300,xx,yy;
      for (int ip=0; ip<g->GetN(); ++ip) { g->GetPoint(ip, xx, yy); ymin = std::min(ymin,yy); ymax = std::max(ymax,yy); }
      if (!(ymax>0)) { ymin = 0.0; ymax = 1.0; }
      const double lo = 0.0;
      const double hi = 2.0*ymax;
      g->GetYaxis()->SetRangeUser(lo, hi);

      g->Draw("AP"); if (fit) fit->Draw("LSAME");

      // Legend LAST (white background), inside UR — widened again, aligned to UR
      DrawUnfoldLegendBox_UR(0.46, 0.62, 0.94, 0.88, fr, fit_sin);
      gPad->Modified(); gPad->Update();

      xcenters.push_back(0.5*(xl+xh));
      Avec.push_back(fr.A); dAvec.push_back(fr.dA);
      Bvec.push_back(fr.B); dBvec.push_back(fr.dB);
      if (fit_sin) { Dvec.push_back(fr.D); dDvec.push_back(fr.dD); Evec.push_back(fr.E); dEvec.push_back(fr.dE); }
    } else {
      TH1D* frame = new TH1D("frame","",10,0,360);
      frame->SetTitle((ttl + ";#phi (deg);Unfolded yield").c_str());
      frame->SetMinimum(0); frame->SetMaximum(1);
      frame->Draw("AXIS");
      DrawUnfoldLegendBox_UR(0.46, 0.66, 0.94, 0.92, FitResult{}, fit_sin);
      gPad->Modified(); gPad->Update();
    }
  }

  // Scale to FUU ratios using depolarization means
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
      double fs2=0, dfs2=0; if (mW!=0.0) { fs2=(mA/mW)*Evec[i]; dfs2=(mA/mW)*dEvec[i]; }
      Fsin.push_back(fs); dFsin.push_back(dfs);
      Fsin2.push_back(fs2); dFsin2.push_back(dfs2);
    }
  }

  // Pad 5 (middle bottom): BLANK (structure-function ratio plot removed)
  c->cd(5); gPad->Clear();

  // Pad 6 (bottom right): amplitudes vs xB (points only, horizontal legend)
  c->cd(6);
  gPad->SetLeftMargin(0.22); gPad->SetRightMargin(0.06);
  gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.12);

  if (!xcenters.empty()) {
    const int N = (int)xcenters.size();
    auto gA  = new TGraphErrors(N);
    auto gB  = new TGraphErrors(N);
    TGraphErrors *gD=nullptr, *gE=nullptr;

    for (int i=0;i<N;i++) {
      gA->SetPoint(i, xcenters[i] - 0.004, Avec[i]); gA->SetPointError(i, 0.0, dAvec[i]);
      gB->SetPoint(i, xcenters[i] + 0.004, Bvec[i]); gB->SetPointError(i, 0.0, dBvec[i]);
      if (fit_sin) {
        if (!gD) gD = new TGraphErrors(N);
        if (!gE) gE = new TGraphErrors(N);
        gD->SetPoint(i, xcenters[i] - 0.012, Dvec[i]); gD->SetPointError(i, 0.0, dDvec[i]);
        gE->SetPoint(i, xcenters[i] + 0.012, Evec[i]); gE->SetPointError(i, 0.0, dEvec[i]);
      }
    }
    gA->SetMarkerColor(kRed);   gA->SetMarkerStyle(20); gA->SetMarkerSize(1.0);
    gB->SetMarkerColor(kBlue);  gB->SetMarkerStyle(21); gB->SetMarkerSize(1.0);
    if (fit_sin) {
      gD->SetMarkerColor(kGreen+2); gD->SetMarkerStyle(22); gD->SetMarkerSize(1.0);
      gE->SetMarkerColor(kMagenta+2); gE->SetMarkerStyle(23); gE->SetMarkerSize(1.0);
    }

    gA->SetTitle(";x_{B};Amplitude");
    gA->GetYaxis()->SetRangeUser(-0.5, 0.5);
    gA->Draw("AP");
    gA->GetXaxis()->SetLimits(0.10, 0.60);
    gB->Draw("P SAME");
    if (fit_sin) { gD->Draw("P SAME"); gE->Draw("P SAME"); }

    auto line0 = new TLine(0.10, 0.0, 0.60, 0.0);
    line0->SetLineStyle(2); line0->SetLineWidth(1); line0->SetLineColor(kBlack);
    line0->Draw("SAME");

    DrawHorizontalLegendAmplitude(fit_sin, gA, gB, (fit_sin?gD:nullptr), (fit_sin?gE:nullptr));
    gPad->Modified(); gPad->Update();
  } else {
    TH1D* frame = new TH1D("frameAmp",";x_{B};Amplitude",10,0.10,0.60);
    frame->SetMinimum(-0.5); frame->SetMaximum(0.5);
    frame->Draw("AXIS");
  }

  ensure_dir("output/enpi+");
  const std::string outpdf = "output/enpi+/" + prop + "_unfolded.pdf";
  c->SaveAs(outpdf.c_str());
  std::cout << "Saved: " << outpdf << "\n";

  // Save arrays (scaled values)
  save_arrays_style(prop, xcenters, Fcos, dFcos, Fcos2, dFcos2, fit_sin, Fsin, dFsin, Fsin2, dFsin2);

  delete c;
}

// ------------------ sine-moments canvas (sinφ & sin2φ histos) ------------------

struct Moments {
  double m1=0, e1=0; // <sinφ>, std error
  double m2=0, e2=0; // <sin2φ>, std error
};

static Moments compute_sine_moments(const std::vector<double>& s1, const std::vector<double>& s2) {
  Moments M;
  const size_t N = s1.size();
  if (N==0) return M;
  double sum1=0, sum2=0, sum1sq=0, sum2sq=0;
  for (size_t i=0;i<N;i++){ sum1+=s1[i]; sum2+=s2[i]; sum1sq+=s1[i]*s1[i]; sum2sq+=s2[i]*s2[i]; }
  M.m1 = sum1 / (double)N;
  M.m2 = sum2 / (double)N;
  double var1 = (N>1? (sum1sq/(double)N - M.m1*M.m1) * (double)N/((double)N-1) : 0.0);
  double var2 = (N>1? (sum2sq/(double)N - M.m2*M.m2) * (double)N/((double)N-1) : 0.0);
  M.e1 = (N>0? std::sqrt(std::max(0.0, var1)) / std::sqrt((double)N) : 0.0);
  M.e2 = (N>0? std::sqrt(std::max(0.0, var2)) / std::sqrt((double)N) : 0.0);
  return M;
}

static void draw_sin_moments_canvas(
  TTree* tD, const BranchNames& bn, const std::string& prop, int nbins_display = 100)
{
  if (!tD) return;

  // global style (match top canvas margins)
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.22);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetTitleSize(0.050, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  auto xe = x_edges();
  const int NX = (int)xe.size()-1;

  // Set up reading
  BranchHandles b;
  if (!bind_tree(tD, bn, b, "DATA_moments")) return;

  // Per xB bin collections
  std::vector<std::unique_ptr<TH1D>> hSin(NX), hSin2(NX);
  std::vector<std::vector<double>> vecSin(NX), vecSin2(NX);
  for (int i=0;i<NX;i++) {
    hSin[i]  = std::make_unique<TH1D>(Form("hSin_%s_%d",  prop.c_str(), i), "", nbins_display, -1.0, 1.0);
    hSin2[i] = std::make_unique<TH1D>(Form("hSin2_%s_%d", prop.c_str(), i), "", nbins_display, -1.0, 1.0);
    hSin[i]->Sumw2(); hSin2[i]->Sumw2();
  }

  // Fill (DATA only) with property window & common cuts
  const Long64_t nent = tD->GetEntries();
  for (Long64_t i=0;i<nent;i++) {
    tD->GetEntry(i);
    const bool ok_abs_t = (std::fabs(b.t) >= 0.10 && std::fabs(b.t) <= 1.20);
    const bool ok_mx2   = (b.mx2 > 0.80) && (b.mx2 < 1.00);
    if (!(ok_abs_t && ok_mx2)) continue;
    if (!b.has_fid || b.fid < 100) continue;
    if (!pass_prop_window(prop, b.t)) continue;

    int ib = xbin_index(b.x, xe);
    if (ib < 0) continue;

    double s1 = std::sin(b.phi);
    double s2 = std::sin(2.0*b.phi);
    hSin[ib]->Fill(s1); hSin2[ib]->Fill(s2);
    vecSin[ib].push_back(s1); vecSin2[ib].push_back(s2);
  }

  // Draw 2x3 canvas: first 4 pads are the histos
  TCanvas* c = new TCanvas(("c_sinmom_"+prop).c_str(), ("Sine moments: "+prop).c_str(), 1500, 1000);
  c->Divide(3,2);
  for (int ib=0; ib<NX; ++ib) {
    c->cd(ib+1);
    gPad->SetLeftMargin(0.22); gPad->SetRightMargin(0.06);
    gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.12);

    // styles
    hSin[ib]->SetLineColor(kBlue);    hSin[ib]->SetLineWidth(3);
    hSin2[ib]->SetLineColor(kRed+1);  hSin2[ib]->SetLineWidth(3);

    // axis titles
    const double xl = xe[ib], xh = xe[ib+1];
    hSin[ib]->SetTitle(Form("%s,  x_{B} #in [%.2f, %.2f),  %s; sin#phi; Events",
                            prop_title(prop).c_str(), xl, xh, prop_tlabel(prop).c_str()));

    // autoscale y to max of both
    double ymax = std::max(hSin[ib]->GetMaximum(), hSin2[ib]->GetMaximum());
    hSin[ib]->SetMaximum(1.10*ymax);
    hSin[ib]->SetMinimum(0.0);
    hSin[ib]->Draw("HIST");
    hSin2[ib]->Draw("HIST SAME");

    // compute means
    Moments M = compute_sine_moments(vecSin[ib], vecSin2[ib]);

    // widened legend box, aligned to UR, border thickness = frame (1)
    const double x1=0.54, y1=0.66, x2=0.94, y2=0.92;
    TPave* box = new TPave(x1,y1,x2,y2,1,"NDC");
    box->SetFillStyle(1001); box->SetFillColor(kWhite);
    box->SetLineColor(kBlack); box->SetLineWidth(1);
    box->Draw("same");

    TLatex lat; lat.SetNDC(); lat.SetTextColor(kBlack);
    lat.SetTextSize(0.030); lat.SetTextAlign(13);
    const double boxH = (y2-y1);
    const double step = 0.40*boxH;
    double ys = y1 + 0.5*(boxH - 2.0*step) + 0.06*boxH; // centered, a hair up
    double xt = x1 + 0.06*(x2-x1);

    lat.SetTextColor(kBlue);
    lat.DrawLatex(xt, ys+step, Form("<sin#phi>  =  %.4f  #pm  %.4f", M.m1, M.e1));
    lat.SetTextColor(kRed+1);
    lat.DrawLatex(xt, ys,      Form("<sin2#phi> =  %.4f  #pm  %.4f", M.m2, M.e2));
  }

  // leave pads 5 and 6 empty on this canvas
  c->cd(5); gPad->Clear();
  c->cd(6); gPad->Clear();

  ensure_dir("output/enpi+");
  const std::string outpdf = "output/enpi+/" + prop + "_sinmoments.pdf";
  c->SaveAs(outpdf.c_str());
  std::cout << "Saved: " << outpdf << "\n";
  delete c;
}

// ------------------ main ------------------

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
    auto list_tree_branches = [](TTree* tr, const char* label){
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
    };
    list_tree_branches(tD, "DATA");
    list_tree_branches(tG, "GEN");
    list_tree_branches(tR, "REC");
  }

  auto H   = make_hist_map(cfg.phi_nbins);
  auto Dep = make_dep_map();

  Counters cD, cG, cR;
  loop_tree_fill(tD, cfg.bn, cfg.apply_fid_cut, H, &Dep, cD, cfg.debug_print, "DATA");
  loop_tree_fill(tG, cfg.bn, /*apply_fid=*/false, H, nullptr, cG, cfg.debug_print, "GEN");
  loop_tree_fill(tR, cfg.bn, cfg.apply_fid_cut, H, nullptr, cR, cfg.debug_print, "REC");

  const auto props = properties_to_run(cfg.which_property);
  print_counters("DATA", cD, props);

  for (const auto& p : props) {
    // main unfolded + amplitudes canvas
    draw_property_and_save(p, H, Dep, cfg.fit_sin);
    // sine-moment overlaid hist canvas
    draw_sin_moments_canvas(tD, cfg.bn, p);
  }

  return 0;
}