// compare_data_mc.cpp
// Compare reconstructed DATA vs MC in your enpi(-t) selections,
// split into the same xB bins as the unfolding code, with overlay+ratio.
//
// • One canvas per (property, xB-bin): a compact 4×2 grid of “cells”.
//   Each cell is internally split into top (overlay, normalized) and bottom (ratio).
// • Variables shown (7 total):
//     1) e_p        [2,9] GeV
//     2) e_theta    [0,40] deg   (stored in tree as radians -> converted to deg)
//     3) p_p        [1,4] GeV
//     4) p_theta    [0,80] deg   (stored in tree as radians -> converted to deg)
//     5) Q2         [1,8] (GeV/c)^2
//     6) xB         [0,1]
//     7) -t         [0,1.2]  (uses fabs(t))
// • Colors: DATA black, MC red.
// • All hists normalized to their own integrals (no magic scale factors).
// • Fast single read per file & property; no copying vectors<unique_ptr>.
//
// Build:
//   g++ -O2 -std=c++17 compare_data_mc.cpp `root-config --cflags --libs` -o compare_data_mc
//
// Run (example with your branch names):
//   ./compare_data_mc data.root mc.root all \
//     --x x --t t --mx2 Mx2 --fid fiducial_status \
//     --ep e_p --eth e_theta --pp p_p --pth p_theta --Q2 Q2
//
// List branches if unsure:
//   ./compare_data_mc data.root mc.root all --list-branches
//
// Outputs:
//   output/compare/compare_<property>_xbin<idx>.pdf   (idx=0..3)
//
// Properties:
//   enpi       : |t| in [0.10,1.20]
//   enpiLowt   : |t| in [0.10,0.4667]
//   enpiMidt   : |t| in [0.4667,0.8333]
//   enpiHight  : |t| in [0.8333,1.20]
//
// Notes:
// - Fiducial cut (fid>=100) applied if the branch exists; otherwise skipped.
// - Mx2 in (0.80,1.00) always applied.
// - Trees are assumed to be named "PhysicsEvents" in both files.

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
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"

// ------------------ CLI / Config ------------------

struct BranchNames {
  // physics
  std::string x   = "x";
  std::string t   = "t";
  std::string mx2 = "Mx2";
  std::string fid = "fiducial_status";   // optional
  // kinematics for plots
  std::string ep  = "e_p";
  std::string eth = "e_theta";           // radians in tree
  std::string pp  = "p_p";
  std::string pth = "p_theta";           // radians in tree
  std::string Q2  = "Q2";
};

struct Config {
  std::string data_path;
  std::string mc_path;
  std::string which_property = "all";
  BranchNames bn;
  bool list_branches = false;
};

static void print_usage(const char* argv0) {
  std::cout << "Usage:\n  " << argv0
            << " data.root mc.root [enpi|enpiLowt|enpiMidt|enpiHight|all]\n"
            << "   [--x x] [--t t] [--mx2 Mx2] [--fid fiducial_status]\n"
            << "   [--ep e_p] [--eth e_theta] [--pp p_p] [--pth p_theta] [--Q2 Q2]\n"
            << "   [--list-branches]\n";
}

static bool parse_args(int argc, char** argv, Config& cfg) {
  if (argc < 4) { print_usage(argv[0]); return false; }
  cfg.data_path = argv[1];
  cfg.mc_path   = argv[2];
  cfg.which_property = argv[3];

  for (int i=4;i<argc;i++) {
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
    if (eat("--x", cfg.bn.x)) continue;
    if (eat("--t", cfg.bn.t)) continue;
    if (eat("--mx2", cfg.bn.mx2)) continue;
    if (eat("--fid", cfg.bn.fid)) continue;
    if (eat("--ep", cfg.bn.ep)) continue;
    if (eat("--eth", cfg.bn.eth)) continue;
    if (eat("--pp", cfg.bn.pp)) continue;
    if (eat("--pth", cfg.bn.pth)) continue;
    if (eat("--Q2", cfg.bn.Q2)) continue;
    if (a == "--list-branches") { cfg.list_branches = true; continue; }

    std::cerr << "Unknown argument: " << a << "\n";
  }
  return true;
}

// ------------------ Helpers: binning & props ------------------

static std::vector<std::string> all_props() {
  return {"enpi", "enpiLowt", "enpiMidt", "enpiHight"};
}
static std::vector<std::string> properties_to_run(const std::string& which) {
  if (which == "all") return all_props();
  return {which};
}

// xB bin edges (4 bins)
static std::vector<double> x_edges() { return {0.10, 0.30, 0.40, 0.50, 0.60}; }
static int xbin_index(double x, const std::vector<double>& e) {
  for (size_t i=0;i+1<e.size();++i) if (x >= e[i] && x < e[i+1]) return (int)i;
  return -1;
}

// |t| windows per property
static inline bool pass_prop_window(const std::string& prop, double t) {
  const double at = std::fabs(t);
  if (prop == "enpi")      return (at >= 0.10 && at <= 1.20);
  if (prop == "enpiLowt")  return (at >= 0.10 && at <= 0.4667);
  if (prop == "enpiMidt")  return (at >= 0.4667 && at <= 0.8333);
  if (prop == "enpiHight") return (at >= 0.8333 && at <= 1.20);
  return false;
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

// ------------------ Branch binding ------------------

struct BranchHandles {
  double x=0, t=0, mx2=0;
  int fid=0; bool has_fid=false;

  double ep=0, eth=0, pp=0, pth=0, Q2=0;

  bool has_x=true, has_t=true, has_mx2=true;
  bool has_ep=true, has_eth=true, has_pp=true, has_pth=true, has_Q2=true;
};

static bool bind_tree(TTree* tr, const BranchNames& bn, BranchHandles& b, const std::string& who) {
  if (!tr) return false;
  auto missing = [&](const char* nm){ std::cerr << "  missing: " << nm << "\n"; };

  // Required for selection
  if (!tr->GetBranch(bn.x.c_str()))  { missing(bn.x.c_str());  b.has_x=false; }
  if (!tr->GetBranch(bn.t.c_str()))  { missing(bn.t.c_str());  b.has_t=false; }
  if (!tr->GetBranch(bn.mx2.c_str())){ missing(bn.mx2.c_str()); b.has_mx2=false; }

  // Variables to plot
  if (!tr->GetBranch(bn.ep.c_str()))  { missing(bn.ep.c_str());  b.has_ep=false; }
  if (!tr->GetBranch(bn.eth.c_str())) { missing(bn.eth.c_str()); b.has_eth=false; }
  if (!tr->GetBranch(bn.pp.c_str()))  { missing(bn.pp.c_str());  b.has_pp=false; }
  if (!tr->GetBranch(bn.pth.c_str())) { missing(bn.pth.c_str()); b.has_pth=false; }
  if (!tr->GetBranch(bn.Q2.c_str()))  { missing(bn.Q2.c_str());  b.has_Q2=false; }

  bool ok = (b.has_x && b.has_t && b.has_mx2 && b.has_ep && b.has_eth && b.has_pp && b.has_pth && b.has_Q2);
  if (!ok) {
    std::cerr << "ERROR: missing required branches on " << who << "\n";
    // still bind what exists to avoid crashes
  }

  if (b.has_x)   tr->SetBranchAddress(bn.x.c_str(),  &b.x);
  if (b.has_t)   tr->SetBranchAddress(bn.t.c_str(),  &b.t);
  if (b.has_mx2) tr->SetBranchAddress(bn.mx2.c_str(),&b.mx2);

  b.has_fid = (tr->GetBranch(bn.fid.c_str()) != nullptr);
  if (b.has_fid) tr->SetBranchAddress(bn.fid.c_str(), &b.fid);

  if (b.has_ep)  tr->SetBranchAddress(bn.ep.c_str(),  &b.ep);
  if (b.has_eth) tr->SetBranchAddress(bn.eth.c_str(), &b.eth);
  if (b.has_pp)  tr->SetBranchAddress(bn.pp.c_str(),  &b.pp);
  if (b.has_pth) tr->SetBranchAddress(bn.pth.c_str(), &b.pth);
  if (b.has_Q2)  tr->SetBranchAddress(bn.Q2.c_str(),  &b.Q2);

  return ok;
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

// ------------------ Variables & histograms ------------------

enum VarID { EP=0, ETH_DEG=1, PP=2, PTH_DEG=3, Q2VAR=4, XB=5, MT=6, NVAR=7 };

struct VarInfo {
  const char* name;
  const char* xtitle;
  int    nbins;
  double xmin, xmax;
};

static VarInfo varInfo[NVAR] = {
  {"e_p",      "e_{p}  (GeV)",                35,  2.0,  9.0},   // EP
  {"e_theta",  "e_{#theta}  (deg)",           40,  0.0, 40.0},   // ETH_DEG (from rad)
  {"p_p",      "p_{p}  (GeV)",                30,  1.0,  4.0},   // PP
  {"p_theta",  "p_{#theta}  (deg)",           40,  0.0, 80.0},   // PTH_DEG (from rad)
  {"Q2",       "Q^{2}  ((GeV/c)^{2})",        35,  1.0,  8.0},   // Q2VAR
  {"xB",       "x_{B}",                        40,  0.0,  1.0},   // XB
  {"mt",       "-t  (GeV^{2})",               40,  0.0,  1.2}    // MT (|t|)
};

struct HistGrid {
  // h[var][xb]
  std::vector<std::vector<TH1D*>> h;
  HistGrid() { h.assign(NVAR, std::vector<TH1D*>(4,nullptr)); }
  ~HistGrid() {
    for (int v=0; v<NVAR; ++v)
      for (int xb=0; xb<4; ++xb)
        if (h[v][xb]) delete h[v][xb];
  }
};

static TH1D* make_hist(const std::string& tag, VarID v, int xb) {
  auto& I = varInfo[v];
  std::string hname = "h_" + std::string(I.name) + "_" + tag + Form("_xb%d", xb);
  auto* h = new TH1D(hname.c_str(), "", I.nbins, I.xmin, I.xmax);
  h->Sumw2();
  h->SetDirectory(nullptr); // avoid ROOT directory warnings
  return h;
}

// ------------------ Filling ------------------

static void build_hists_one_pass(
  TTree* tr, const BranchNames& bn, const std::string& prop,
  HistGrid& HG, bool isData)
{
  if (!tr) return;
  BranchHandles b;
  bool ok = bind_tree(tr, bn, b, std::string(isData?"DATA":"MC")+"("+prop+")");

  // Make hists if not yet
  const auto xe = x_edges();
  for (int xb=0; xb<4; ++xb) {
    if (!HG.h[EP][xb])      HG.h[EP][xb]      = make_hist(std::string(isData?"D":"M")+"_"+prop, EP, xb);
    if (!HG.h[ETH_DEG][xb]) HG.h[ETH_DEG][xb] = make_hist(std::string(isData?"D":"M")+"_"+prop, ETH_DEG, xb);
    if (!HG.h[PP][xb])      HG.h[PP][xb]      = make_hist(std::string(isData?"D":"M")+"_"+prop, PP, xb);
    if (!HG.h[PTH_DEG][xb]) HG.h[PTH_DEG][xb] = make_hist(std::string(isData?"D":"M")+"_"+prop, PTH_DEG, xb);
    if (!HG.h[Q2VAR][xb])   HG.h[Q2VAR][xb]   = make_hist(std::string(isData?"D":"M")+"_"+prop, Q2VAR, xb);
    if (!HG.h[XB][xb])      HG.h[XB][xb]      = make_hist(std::string(isData?"D":"M")+"_"+prop, XB, xb);
    if (!HG.h[MT][xb])      HG.h[MT][xb]      = make_hist(std::string(isData?"D":"M")+"_"+prop, MT, xb);
  }

  if (!ok) {
    std::cerr << "Proceeding with whatever branches are available on "
              << (isData?"DATA":"MC") << "(" << prop << ").\n";
  }

  const Long64_t nent = tr->GetEntries();
  for (Long64_t i=0;i<nent;i++) {
    tr->GetEntry(i);

    if (b.has_x && b.has_t && b.has_mx2) {
      if (!(std::fabs(b.t) >= 0.10 && std::fabs(b.t) <= 1.20)) continue;
      if (!pass_prop_window(prop, b.t)) continue;
      if (!(b.mx2 > 0.80 && b.mx2 < 1.00)) continue;
    } else {
      continue;
    }

    // fiducial if present
    if (b.has_fid && b.fid < 100) continue;

    int xb = xbin_index(b.x, xe);
    if (xb < 0) continue;

    // Fill variables (respect missing-branch flags)
    if (b.has_ep)  HG.h[EP][xb]->Fill(b.ep);
    if (b.has_eth) HG.h[ETH_DEG][xb]->Fill(b.eth * 180.0/TMath::Pi());
    if (b.has_pp)  HG.h[PP][xb]->Fill(b.pp);
    if (b.has_pth) HG.h[PTH_DEG][xb]->Fill(b.pth * 180.0/TMath::Pi());
    if (b.has_Q2)  HG.h[Q2VAR][xb]->Fill(b.Q2);
    if (b.has_x)   HG.h[XB][xb]->Fill(b.x);
    if (b.has_t)   HG.h[MT][xb]->Fill(std::fabs(b.t));
  }
}

// ------------------ Drawing (overlay + ratio in one cell) ------------------

struct Colors {
  int data = kBlack;
  int mc   = kRed+1;
};

static void style_hist_line(TH1D* h, int color, int width=2) {
  h->SetLineColor(color);
  h->SetLineWidth(width);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}

static double safe_integral(const TH1D* h) {
  if (!h) return 0.0;
  double sum = 0.0;
  for (int i=1;i<=h->GetNbinsX();++i) sum += h->GetBinContent(i);
  return sum;
}

static TH1D* clone_detached(const TH1D* src, const std::string& newname) {
  auto* c = (TH1D*)src->Clone(newname.c_str());
  c->SetDirectory(nullptr);
  return c;
}

static void normalize_in_place(TH1D* h) {
  double I = safe_integral(h);
  if (I > 0) h->Scale(1.0/I);
}

static void compute_ratio(const TH1D* hD, const TH1D* hM, TH1D* hR) {
  const int nb = hR->GetNbinsX();
  for (int i=1;i<=nb;++i) {
    const double d = hD->GetBinContent(i);
    const double m = hM->GetBinContent(i);
    double r=0, er=0;
    if (m > 0) {
      r = d/m;
      double ed = hD->GetBinError(i);
      double em = hM->GetBinError(i);
      er = std::sqrt( (ed*ed)/(m*m) + (d*d*em*em)/(m*m*m*m) );
    }
    hR->SetBinContent(i, r);
    hR->SetBinError(i, er);
  }
}

// Draw overlay+ratio inside the current pad (we split it into two subpads)
static void draw_dist_ratio_cell(TH1D* hD, TH1D* hM,
                                 const std::string& xtitle, const std::string& ytitle,
                                 double xmin, double xmax,
                                 const Colors& col)
{
  // Create subpads
  TPad* pTop = new TPad("pTop","",0.00,0.32,1.00,1.00);  // ~68% height
  TPad* pBot = new TPad("pBot","",0.00,0.00,1.00,0.32);  // ~32% height
  pTop->SetBottomMargin(0.10);
  pTop->SetLeftMargin(0.16); pTop->SetRightMargin(0.05); pTop->SetTopMargin(0.08);
  pBot->SetTopMargin(0.02);  pBot->SetBottomMargin(0.35);
  pBot->SetLeftMargin(0.16); pBot->SetRightMargin(0.05);
  pTop->Draw(); pBot->Draw();

  // --- TOP: overlay (normalized) ---
  pTop->cd();
  TH1D* dN = nullptr; TH1D* mN = nullptr;
  if (hD) dN = clone_detached(hD, std::string(hD->GetName())+"_N");
  if (hM) mN = clone_detached(hM, std::string(hM->GetName())+"_N");

  if (dN) normalize_in_place(dN);
  if (mN) normalize_in_place(mN);

  // Axis frame from data if exists, else MC
  TH1D* frame = nullptr;
  if (dN && safe_integral(dN) > 0) frame = dN;
  else if (mN && safe_integral(mN) > 0) frame = mN;

  if (frame) {
    frame->SetTitle( (";"+xtitle+";"+ytitle).c_str() );
    frame->GetXaxis()->SetRangeUser(xmin, xmax);

    double ymax = 0.0;
    if (dN) ymax = std::max(ymax, dN->GetMaximum());
    if (mN) ymax = std::max(ymax, mN->GetMaximum());
    frame->SetMaximum(1.15*ymax);
    frame->SetMinimum(0.0);

    style_hist_line(frame, (frame==dN? col.data : col.mc), 2);
    frame->Draw("HIST");
    if (dN && frame!=dN) { style_hist_line(dN, col.data, 2); dN->Draw("HIST SAME"); }
    if (mN && frame!=mN) { style_hist_line(mN, col.mc,   2); mN->Draw("HIST SAME"); }
  } else {
    TH1D* ax = new TH1D("ax","",10,xmin,xmax);
    ax->SetDirectory(nullptr);
    ax->SetTitle( (";"+xtitle+";"+ytitle).c_str() );
    ax->SetMinimum(0.0); ax->SetMaximum(1.0);
    ax->Draw("AXIS");
  }

  // --- BOTTOM: ratio (Data/MC of normalized) ---
  pBot->cd();
  TH1D* r = nullptr;
  if (dN && mN) {
    r = clone_detached(dN, std::string(dN->GetName()) + "_ratio"); // <-- FIX: use helper (c_str handled)
    compute_ratio(dN, mN, r);
    r->SetTitle( (";"+xtitle+";Data/MC").c_str() );
    r->GetXaxis()->SetRangeUser(xmin, xmax);
    r->GetYaxis()->SetTitleSize(0.11);
    r->GetYaxis()->SetLabelSize(0.10);
    r->GetXaxis()->SetTitleSize(0.12);
    r->GetXaxis()->SetLabelSize(0.10);

    double rmin=+1e9, rmax=-1e9;
    for (int i=1;i<=r->GetNbinsX();++i) {
      double v = r->GetBinContent(i);
      if (v<=0) continue;
      rmin = std::min(rmin,v); rmax = std::max(rmax,v);
    }
    if (rmin>rmax) { rmin=0.5; rmax=1.5; }
    double lo = std::min(0.5, rmin*0.9);
    double hi = std::max(1.5, rmax*1.1);
    r->SetMinimum(lo); r->SetMaximum(hi);

    style_hist_line(r, kBlack, 2);
    r->Draw("HIST");

    TLine* l1 = new TLine(xmin,1.0,xmax,1.0);
    l1->SetLineStyle(2); l1->SetLineWidth(1); l1->SetLineColor(kBlack);
    l1->Draw("SAME");
  } else {
    TH1D* ax = new TH1D("rax","",10,xmin,xmax);
    ax->SetDirectory(nullptr);
    ax->SetTitle( (";"+xtitle+";Data/MC").c_str() );
    ax->SetMinimum(0.5); ax->SetMaximum(1.5);
    ax->GetYaxis()->SetTitleSize(0.11);
    ax->GetYaxis()->SetLabelSize(0.10);
    ax->GetXaxis()->SetTitleSize(0.12);
    ax->GetXaxis()->SetLabelSize(0.10);
    ax->Draw("AXIS");
    TLine* l1 = new TLine(xmin,1.0,xmax,1.0);
    l1->SetLineStyle(2); l1->SetLineWidth(1); l1->SetLineColor(kBlack);
    l1->Draw("SAME");
  }

  if (dN) delete dN;
  if (mN) delete mN;
  if (r)  delete r;
}

// ------------------ Canvas per x-bin ------------------

static std::string prop_tlabel(const std::string&); // fwd decl

static void draw_global_legend_and_labels(
  const std::string& prop, int xbin, const std::vector<double>& xe, const Colors& col)
{
  gPad->Clear();
  gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.12);

  TH1D* h = new TH1D("hHost","",10,0,1);
  h->SetDirectory(nullptr);
  h->SetMinimum(0); h->SetMaximum(1);
  h->SetTitle("; ; ");
  h->Draw("AXIS");

  auto L = new TLegend(0.15, 0.65, 0.85, 0.88);
  L->SetFillStyle(1001); L->SetFillColor(kWhite);
  L->SetBorderSize(1);   L->SetTextSize(0.035);
  TH1D* hD = new TH1D("hDummyData","",1,0,1);  hD->SetDirectory(nullptr); style_hist_line(hD, col.data, 2);
  TH1D* hM = new TH1D("hDummyMC","",1,0,1);    hM->SetDirectory(nullptr); style_hist_line(hM, col.mc,   2);
  L->AddEntry(hD, "DATA (norm. to integral)", "l");
  L->AddEntry(hM, "MC (norm. to integral)",   "l");
  L->Draw();

  TLatex lat; lat.SetNDC();
  lat.SetTextSize(0.040); lat.SetTextAlign(13);
  std::string tlabel = prop_tlabel(prop);
  std::string xblabel = Form("x_{B} #in [%.2f, %.2f)", xe[xbin], xe[xbin+1]);

  lat.DrawLatex(0.15, 0.58, Form("Property: %s", prop.c_str()));
  lat.DrawLatex(0.15, 0.52, tlabel.c_str());
  lat.DrawLatex(0.15, 0.46, xblabel.c_str());

  delete hD; delete hM;
}

static std::vector<VarID> default_var_order() {
  return { EP, ETH_DEG, PP, PTH_DEG, Q2VAR, XB, MT };
}

static void draw_canvas_for_xbin(
  const std::string& prop, int xbin,
  const HistGrid& D, const HistGrid& M)
{
  ensure_dir("output/compare");
  auto xe = x_edges();

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas(Form("c_%s_xb%d", prop.c_str(), xbin),
                           Form("Compare DATA vs MC: %s, xB bin %d", prop.c_str(), xbin),
                           1600, 900);
  c->Divide(4,2); // 8 cells; each cell has overlay+ratio inside

  Colors col; // black data, red MC
  auto vars = default_var_order();

  for (size_t iv=0; iv<vars.size(); ++iv) {
    int cell = (int)iv + 1;
    c->cd(cell);

    VarID v = vars[iv];
    TH1D* hD = D.h[v][xbin];
    TH1D* hM = M.h[v][xbin];

    auto& I = varInfo[v];
    std::string xtitle = I.xtitle;
    std::string ytitle = "Normalized counts";

    draw_dist_ratio_cell(hD, hM, xtitle, ytitle, I.xmin, I.xmax, col);
  }

  c->cd(8);
  draw_global_legend_and_labels(prop, xbin, xe, col);

  std::string outpdf = Form("output/compare/compare_%s_xbin%d.pdf", prop.c_str(), xbin);
  c->Print(outpdf.c_str(), "pdf");
  std::cout << "Saved: " << outpdf << "\n";
  delete c;
}

// ------------------ main ------------------

int main(int argc, char** argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) return 1;

  std::unique_ptr<TFile> fD(TFile::Open(cfg.data_path.c_str(), "READ"));
  std::unique_ptr<TFile> fM(TFile::Open(cfg.mc_path.c_str(),   "READ"));
  if (!fD || fD->IsZombie() || !fM || fM->IsZombie()) {
    std::cerr << "ERROR: could not open one or both input files.\n";
    return 2;
  }

  TTree* tD = (TTree*)fD->Get("PhysicsEvents");
  TTree* tM = (TTree*)fM->Get("PhysicsEvents");
  if (!tD || !tM) {
    std::cerr << "ERROR: missing PhysicsEvents tree in one or both files.\n";
    if (tD) std::cerr << "  DATA has it; MC missing.\n";
    if (tM) std::cerr << "  MC has it; DATA missing.\n";
    return 3;
  }

  if (cfg.list_branches) {
    list_tree_branches(tD, "DATA");
    list_tree_branches(tM, "MC");
  }

  auto props = properties_to_run(cfg.which_property);

  for (const auto& p : props) {
    HistGrid HD, HM;
    build_hists_one_pass(tD, cfg.bn, p, HD, /*isData=*/true);
    build_hists_one_pass(tM, cfg.bn, p, HM, /*isData=*/false);

    for (int xb=0; xb<4; ++xb) {
      draw_canvas_for_xbin(p, xb, HD, HM);
    }
  }

  return 0;
}