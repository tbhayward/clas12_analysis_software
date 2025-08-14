// unfold_phi_unfold_fit.cpp
// Exclusive pi+ unfolding and cos(n*phi) fits (phi plotted in degrees).
// Build (csh):
//   g++ -O2 -std=c++17 unfold_phi_unfold_fit.cpp `root-config --cflags --libs` -o unfold_phi_unfold_fit
// Run:
//   ./unfold_phi_unfold_fit data.root gen.root rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]
//                                [--x x] [--phi phi] [--t t] [--mx2 Mx2] [--fid fiducial_status]
//                                [--phibins 24] [--no-fid-cut] [--debug N] [--list-branches]
//
// Outputs per property:
//   PDF:  output/enpi+/<property>_unfolded.pdf
//   TXT:  output/enpi+/<property>_unfolded_fit_results.txt
//
// Notes:
// - Reads each tree ONCE per dataset; fills all four properties simultaneously.
// - Input phi assumed in radians. Plots and fits in degrees (0..360).
// - Unfolded yield per phi bin = D * G / R, with Poisson error propagation.
// - Fits C * (1 + A*cos(pi/180*x) + B*cos(2*pi/180*x)) on unfolded points.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSystem.h"
#include "TTree.h"

struct BranchNames {
  std::string x   = "x";
  std::string phi = "phi";
  std::string t   = "t";
  std::string mx2 = "Mx2";
  std::string fid = "fiducial_status";
};

struct Config {
  std::string data_path;
  std::string gen_path;
  std::string rec_path;
  std::string which_property = "all";
  BranchNames bn;
  int  phi_nbins      = 24;
  bool apply_fid_cut  = true;
  int  debug_print    = 0;  // print this many raw events per dataset
  bool list_branches  = false;
};

// ---------- CLI ----------
static void print_usage(const char* argv0) {
  std::cout << "Usage:\n  " << argv0
            << " data.root gen.root rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]\n"
            << "    [--x x] [--phi phi] [--t t] [--mx2 Mx2] [--fid fiducial_status]\n"
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
    if (eat_int("--phibins", cfg.phi_nbins)) continue;
    if (a == "--no-fid-cut") { cfg.apply_fid_cut = false; continue; }
    if (eat_int("--debug", cfg.debug_print)) continue;
    if (a == "--list-branches") { cfg.list_branches = true; continue; }

    std::cerr << "Unknown argument: " << a << "\n";
  }
  return true;
}

// ---------- Helpers ----------
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
  if (prop == "enpi")      return true;
  if (prop == "enpiLowt")  return (std::fabs(t) >= 0.00 && std::fabs(t) <= 0.30);
  if (prop == "enpiMidt")  return (std::fabs(t) >= 0.30 && std::fabs(t) <= 0.70);
  if (prop == "enpiHight") return (std::fabs(t) >= 0.70); // bounded by common
  return false;
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

// ---------- Binding ----------
struct BranchHandles {
  double x=0, phi=0, t=0, mx2=0;
  int fid=0;
  bool has_fid=false;
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
  return true;
}

// ---------- Hist containers ----------
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

  const double pmin = 0.0;
  const double pmax = 360.0;

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

// ---------- Debug counters ----------
struct Counters {
  Long64_t total=0, pass_common_cnt=0, pass_fid_cnt=0;
  Long64_t fail_abs_t=0, fail_mx2=0, fail_no_fid_branch=0, fail_fid_lt_100=0, fail_xbin=-0;
  std::map<std::string, Long64_t> pass_prop_cnt;
};

static void debug_event_print(const char* label, int idx, const BranchHandles& b,
                              const std::vector<std::string>& props, const std::vector<double>& xe) {
  std::cout << "[" << label << " evt " << idx << "] "
            << "x=" << b.x << "  t=" << b.t << "  Mx2=" << b.mx2
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

// ---------- Fill loops ----------
static void loop_tree_fill(
  TTree* tr, const BranchNames& bn, bool apply_fid, bool count_fid,
  std::map<std::string,HSet>& H, Counters& C, int debugN, const char* dbg_label)
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

    // Unconditional early prints for first N entries:
    if (debug_limit > 0 && i < debug_limit) {
      debug_event_print(dbg_label, (int)i, b, props, xe);
    }

    // Track reasons for failures:
    bool ok_abs_t = (std::fabs(b.t) <= 1.0) && (std::fabs(b.t) >= 0.0);
    bool ok_mx2   = (b.mx2 > 0.80) && (b.mx2 < 1.00);
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
      if (std::string(dbg_label) == "DATA") H[p].D[ib]->Fill(phideg);
      else if (std::string(dbg_label) == "GEN") H[p].G[ib]->Fill(phideg);
      else if (std::string(dbg_label) == "REC") H[p].R[ib]->Fill(phideg);
    }

    // Also print every ~10M entries to show progress on huge files
    if ((i>0) && (i % 10000000 == 0)) {
      std::cout << "[" << dbg_label << "] processed " << i << " / " << nent << " entries..." << std::endl;
    }
  }
}

// ---------- Fit ----------
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

    // U = D*G/R
    const double U = (D * G) / R;

    // Poisson propagation: Var(U) = (G/R)^2*D + (D/R)^2*G + (D*G/R^2)^2*R
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

// ---------- I/O ----------
static std::string prop_title(const std::string& p) {
  if (p == "enpi")      return "ep -> e' n #pi^{+}";
  if (p == "enpiLowt")  return "ep -> e' n #pi^{+}, low |t|";
  if (p == "enpiMidt")  return "ep -> e' n #pi^{+}, mid |t|";
  if (p == "enpiHight") return "ep -> e' n #pi^{+}, high |t|";
  return p;
}

static void save_results_txt(
  const std::string& prop,
  const std::vector<double>& xcenters,
  const std::vector<double>& xlo,
  const std::vector<double>& xhi,
  const std::vector<double>& Avec, const std::vector<double>& dAvec,
  const std::vector<double>& Bvec, const std::vector<double>& dBvec,
  const std::vector<double>& Cvec, const std::vector<double>& dCvec,
  const std::vector<double>& chi2, const std::vector<int>& ndf)
{
  ensure_dir("output/enpi+");
  std::string out = "output/enpi+/" + prop + "_unfolded_fit_results.txt";
  std::ofstream ofs(out);
  ofs << "# property: " << prop << "\n";
  ofs << "# columns: i  x_center  x_lo  x_hi  C  dC  A  dA  B  dB  chi2  ndf\n";
  for (size_t i=0;i<xcenters.size();++i) {
    ofs << i << " "
        << xcenters[i] << " "
        << xlo[i] << " "
        << xhi[i] << " "
        << Cvec[i] << " " << dCvec[i] << " "
        << Avec[i] << " " << dAvec[i] << " "
        << Bvec[i] << " " << dBvec[i] << " "
        << chi2[i] << " " << ndf[i] << "\n";
  }
  ofs.close();
  std::cout << "Wrote results: " << out << "\n";
}

static void draw_property_and_save(
  const std::string& prop,
  const std::map<std::string,HSet>& H)
{
  const auto xe = x_edges();
  const int NX = (int)xe.size()-1;

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.125);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  TCanvas* c = new TCanvas(("c_"+prop).c_str(), ("Unfolded phi fits: "+prop).c_str(), 1200, 900);
  c->Divide(3,3);

  std::vector<double> xcenters, xlos, xhis;
  std::vector<double> Avec, dAvec, Bvec, dBvec, Cvec, dCvec, Chi2;
  std::vector<int>    NDF;

  for (int ib=0; ib<NX; ++ib) {
    c->cd(ib+1);

    TH1D* hD = H.at(prop).D[ib].get();
    TH1D* hG = H.at(prop).G[ib].get();
    TH1D* hR = H.at(prop).R[ib].get();

    TGraphErrors* g = nullptr;
    TF1* fit = nullptr;
    FitResult fr = make_unfold_graph_and_fit(hD,hG,hR,g,fit);

    const double xl = xe[ib], xh = xe[ib+1];
    std::string title = prop_title(prop) + Form(", x in [%.2f, %.2f)", xl, xh);

    if (g) {
      g->SetTitle((title + ";phi (deg);Unfolded yield").c_str());
      g->Draw("AP");
      if (fit) fit->Draw("LSAME");

      auto leg = new TLegend(0.55, 0.58, 0.94, 0.94);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry((TObject*)0, Form("C = %.3g +/- %.3g", fr.C, fr.dC), "");
      leg->AddEntry((TObject*)0, Form("A = %.3g +/- %.3g", fr.A, fr.dA), "");
      leg->AddEntry((TObject*)0, Form("B = %.3g +/- %.3g", fr.B, fr.dB), "");
      leg->AddEntry((TObject*)0, Form("chi2/ndf = %.1f/%d", fr.chi2, fr.ndf), "");
      leg->Draw();

      xcenters.push_back(0.5*(xl+xh));
      xlos.push_back(xl);
      xhis.push_back(xh);
      Avec.push_back(fr.A); dAvec.push_back(fr.dA);
      Bvec.push_back(fr.B); dBvec.push_back(fr.dB);
      Cvec.push_back(fr.C); dCvec.push_back(fr.dC);
      Chi2.push_back(fr.chi2); NDF.push_back(fr.ndf);
    } else {
      TH1D* frame = new TH1D("frame","",10,0,360);
      frame->SetTitle((title + ";phi (deg);Unfolded yield").c_str());
      frame->SetMinimum(0); frame->SetMaximum(1);
      frame->Draw("AXIS");
      auto leg = new TLegend(0.55, 0.70, 0.94, 0.90);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      leg->AddEntry((TObject*)0, "No valid unfolded points", "");
      leg->Draw();
    }
  }

  // Bottom row: pad 8 shows A(x), B(x)
  c->cd(7); gPad->Clear();
  c->cd(8);
  if (!xcenters.empty()) {
    auto gA = new TGraphErrors((int)xcenters.size());
    auto gB = new TGraphErrors((int)xcenters.size());
    for (int i=0;i<(int)xcenters.size();++i) {
      gA->SetPoint(i, xcenters[i], Avec[i]); gA->SetPointError(i, 0.0, dAvec[i]);
      gB->SetPoint(i, xcenters[i], Bvec[i]); gB->SetPointError(i, 0.0, dBvec[i]);
    }
    gA->SetMarkerStyle(20); gA->SetMarkerSize(1.0);
    gB->SetMarkerStyle(21); gB->SetMarkerSize(1.0);

    gA->SetTitle("A(x) and B(x);x;Coefficient");
    gA->GetYaxis()->SetRangeUser(-1.0, 1.0);
    gA->Draw("AP");
    gB->Draw("P SAME");

    auto legAB = new TLegend(0.15, 0.75, 0.45, 0.92);
    legAB->SetBorderSize(0); legAB->SetFillStyle(0);
    legAB->AddEntry(gA, "A(x) from fits", "p");
    legAB->AddEntry(gB, "B(x) from fits", "p");
    legAB->Draw();
  } else {
    TH1D* frame = new TH1D("frameAB","A(x) and B(x);x;Coefficient",10,0.06,0.60);
    frame->SetMinimum(-1.0); frame->SetMaximum(1.0);
    frame->Draw("AXIS");
    auto leg = new TLegend(0.2, 0.8, 0.6, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry((TObject*)0, "No fit results collected", "");
    leg->Draw();
  }
  c->cd(9); gPad->Clear();

  ensure_dir("output/enpi+");
  const std::string outpdf = "output/enpi+/" + prop + "_unfolded.pdf";
  c->SaveAs(outpdf.c_str());
  std::cout << "Saved: " << outpdf << "\n";

  save_results_txt(prop, xcenters, xlos, xhis, Avec, dAvec, Bvec, dBvec, Cvec, dCvec, Chi2, NDF);
}

// ---------- Final reporting ----------
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

// ---------- main ----------
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

  Counters cD, cG, cR;
  loop_tree_fill(tD, cfg.bn, cfg.apply_fid_cut, /*count_fid=*/true,  H, cD, cfg.debug_print, "DATA");
  loop_tree_fill(tG, cfg.bn, /*apply_fid=*/false, /*count_fid=*/false, H, cG, cfg.debug_print, "GEN");
  loop_tree_fill(tR, cfg.bn, cfg.apply_fid_cut, /*count_fid=*/true,  H, cR, cfg.debug_print, "REC");

  const auto props = properties_to_run(cfg.which_property);

  // Totals sanity
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

  // Detailed counters
  print_counters("DATA", cD, props);
  print_counters("GEN",  cG, props);
  print_counters("REC",  cR, props);

  // Draw and save
  for (const auto& p : props) {
    draw_property_and_save(p, H);
  }

  return 0;
}