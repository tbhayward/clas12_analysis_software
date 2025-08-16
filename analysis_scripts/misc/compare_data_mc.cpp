// compare_data_mc.cpp
// Overlay DATA (black) vs MC (red) for key reconstructed observables,
// normalized HIST-BY-HIST to unit area (shape comparison), plus a ratio panel.
// Uses same kinematic binning you’ve been using:
//
//   x bins: [0.10,0.30), [0.30,0.40), [0.40,0.50), [0.50,0.60)
//   t-slices ("property"):
//     enpi      : 0.10 <= |t| <= 1.20
//     enpiLowt  : 0.10 <= |t| <= 0.4667
//     enpiMidt  : 0.4667 <= |t| <= 0.8333
//     enpiHight : 0.8333 <= |t| <= 1.20
//
// Cuts (DATA & MC): fiducial_status >= 100 (if branch exists); 0.80 < Mx2 < 1.00
//
// Variables & fixed display ranges (all normalized independently to unit area):
//   e_p      : 2 -> 9  GeV
//   e_theta  : 0 -> 40 deg (branch is radians; converted here)
//   p_p      : 1 -> 4  GeV
//   p_theta  : 0 -> 80 deg (branch is radians; converted here)
//   Q2       : 1 -> 8
//   xB       : 0 -> 1
//   -t       : 0 -> 1.2  (we plot |t|)
//
// Output: one canvas per (property, x-bin), 4x2 cells. Each cell contains
// an upper overlay (data vs mc) and a lower ratio (Data/MC).
//
// Build:
//   g++ -O2 -std=c++17 compare_data_mc.cpp `root-config --cflags --libs` -o compare_data_mc
//
// Run:
//   ./compare_data_mc data_rec.root mc_rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]
//        [--x x] [--t t] [--mx2 Mx2] [--fid fiducial_status]
//        [--ep e_p] [--eth e_theta] [--pp p_p] [--pth p_theta] [--Q2 Q2]
//        [--list-branches]

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
#include <string>
#include <vector>

#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

// ------------- CLI -------------

struct BranchNames {
  // shared kinematics/cuts
  std::string x   = "x";                // xB
  std::string t   = "t";                // signed t
  std::string mx2 = "Mx2";
  std::string fid = "fiducial_status";
  // reconstructed observables
  std::string ep  = "e_p";              // GeV
  std::string eth = "e_theta";          // radians in tree
  std::string pp  = "p_p";              // GeV
  std::string pth = "p_theta";          // radians in tree
  std::string Q2  = "Q2";
};

struct Config {
  std::string data_path, mc_path;
  std::string which_property = "all";
  BranchNames bn;
  bool list_branches = false;
};

static void print_usage(const char* argv0) {
  std::cout << "Usage:\n  " << argv0
            << " data_rec.root mc_rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]\n"
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
    if (eat("--x",   cfg.bn.x)) continue;
    if (eat("--t",   cfg.bn.t)) continue;
    if (eat("--mx2", cfg.bn.mx2)) continue;
    if (eat("--fid", cfg.bn.fid)) continue;
    if (eat("--ep",  cfg.bn.ep))  continue;
    if (eat("--eth", cfg.bn.eth)) continue;
    if (eat("--pp",  cfg.bn.pp))  continue;
    if (eat("--pth", cfg.bn.pth)) continue;
    if (eat("--Q2",  cfg.bn.Q2))  continue;
    if (a == "--list-branches") { cfg.list_branches = true; continue; }
    std::cerr << "Unknown argument: " << a << "\n";
  }
  return true;
}

// ------------- binning & selections -------------

static std::vector<std::string> all_props() {
  return {"enpi", "enpiLowt", "enpiMidt", "enpiHight"};
}
static std::vector<std::string> properties_to_run(const std::string& which) {
  if (which == "all") return all_props();
  return {which};
}
static std::vector<double> x_edges() { return {0.10, 0.30, 0.40, 0.50, 0.60}; }
static int xbin_index(double x, const std::vector<double>& e) {
  for (size_t i=0;i+1<e.size();++i) if (x >= e[i] && x < e[i+1]) return (int)i;
  return -1;
}
static inline bool pass_prop_window(const std::string& prop, double t) {
  const double at = std::fabs(t);
  if (prop == "enpi")      return (at >= 0.10 && at <= 1.20);
  if (prop == "enpiLowt")  return (at >= 0.10 && at <= 0.4667);
  if (prop == "enpiMidt")  return (at >= 0.4667 && at <= 0.8333);
  if (prop == "enpiHight") return (at >= 0.8333 && at <= 1.20);
  return false;
}
static void ensure_dir(const std::string& d) {
  if (gSystem->AccessPathName(d.c_str())) gSystem->mkdir(d.c_str(), kTRUE);
}

// ------------- binding -------------

struct Handles {
  double x=0, t=0, mx2=0;
  int fid=0; bool has_fid=false;
  double ep=0, eth=0, pp=0, pth=0, Q2=0;
};

static bool bind_tree(TTree* tr, const BranchNames& bn, Handles& h, const char* who) {
  if (!tr) return false;
  auto req = [&](const std::string& b){ return tr->GetBranch(b.c_str()) != nullptr; };
  if (!(req(bn.x) && req(bn.t) && req(bn.mx2) && req(bn.ep) && req(bn.eth) &&
        req(bn.pp) && req(bn.pth) && req(bn.Q2))) {
    std::cerr << "ERROR: missing required branches on " << who << "\n";
    return false;
  }
  tr->SetBranchAddress(bn.x.c_str(),   &h.x);
  tr->SetBranchAddress(bn.t.c_str(),   &h.t);
  tr->SetBranchAddress(bn.mx2.c_str(), &h.mx2);
  tr->SetBranchAddress(bn.ep.c_str(),  &h.ep);
  tr->SetBranchAddress(bn.eth.c_str(), &h.eth);
  tr->SetBranchAddress(bn.pp.c_str(),  &h.pp);
  tr->SetBranchAddress(bn.pth.c_str(), &h.pth);
  tr->SetBranchAddress(bn.Q2.c_str(),  &h.Q2);
  h.has_fid = (tr->GetBranch(bn.fid.c_str()) != nullptr);
  if (h.has_fid) tr->SetBranchAddress(bn.fid.c_str(), &h.fid);
  return true;
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

// ------------- histogram helpers -------------

static void normalize_to_area(TH1* h) {
  if (!h) return;
  double integral = h->Integral();
  if (integral > 0) h->Scale(1.0 / integral);
}

struct VarDef {
  std::string name;         // internal label
  std::string title;        // top pad title (left of / above legend)
  std::string xtitle;       // x-axis label
  int    nbins;
  double xmin, xmax;
  bool   use_degrees;       // true for *_theta vars (radian->degree)
};

static std::vector<VarDef> var_list() {
  // Bin counts chosen to be smooth but not noisy
  return {
    {"e_p",     "e_{p}",          "e_{p} (GeV)",          70, 2.0,  9.0,  false},
    {"e_theta", "e_{#theta}",     "e_{#theta} (deg)",     80, 0.0, 40.0,  true },
    {"p_p",     "p_{p}",          "p_{p} (GeV)",          60, 1.0,  4.0,  false},
    {"p_theta", "p_{#theta}",     "p_{#theta} (deg)",     80, 0.0, 80.0,  true },
    {"Q2",      "Q^{2}",          "Q^{2} (GeV^{2})",      70, 1.0,  8.0,  false},
    {"xB",      "x_{B}",          "x_{B}",                60, 0.0,  1.0,  false},
    {"t",       "-t",             "-t (GeV^{2})",         60, 0.0,  1.2,  false}
  };
}

// Extract a value from Handles based on VarDef
static inline double get_value(const Handles& h, const VarDef& v) {
  if (v.name == "e_p")     return h.ep;
  if (v.name == "e_theta") return h.eth * 180.0 / TMath::Pi();
  if (v.name == "p_p")     return h.pp;
  if (v.name == "p_theta") return h.pth * 180.0 / TMath::Pi();
  if (v.name == "Q2")      return h.Q2;
  if (v.name == "xB")      return h.x;
  if (v.name == "t")       return std::fabs(h.t);
  return 0.0;
}

// Create and fill a normalized hist for a given (tree, property, xbin, var)
static std::unique_ptr<TH1D> make_hist_for(
  TTree* tr, const BranchNames& bn, const std::string& prop, int xbin,
  const VarDef& var, const std::vector<double>& xedges)
{
  auto h = std::make_unique<TH1D>(
    Form("h_%s_prop_%s_xbin_%d", var.name.c_str(), prop.c_str(), xbin),
    "", var.nbins, var.xmin, var.xmax);
  h->Sumw2();

  Handles hd; if (!bind_tree(tr, bn, hd, "TREE")) return h;

  const Long64_t N = tr->GetEntries();
  for (Long64_t i=0;i<N;i++) {
    tr->GetEntry(i);
    // cuts
    if (hd.has_fid && hd.fid < 100) continue;
    if (!(hd.mx2 > 0.80 && hd.mx2 < 1.00)) continue;
    if (!pass_prop_window(prop, hd.t)) continue;

    int ib = xbin_index(hd.x, xedges);
    if (ib != xbin) continue;

    double val = get_value(hd, var);
    if (val >= var.xmin && val < var.xmax) h->Fill(val);
  }

  normalize_to_area(h.get());
  return h;
}

// ------------- draw split cell (overlay + ratio) -------------

struct Colors { int data=1; int mc=2; }; // kBlack=1, kRed=2

static void draw_one_cell(
  TPad* parent, TH1D* hD, TH1D* hM,
  const std::string& title, const std::string& xtitle,
  double xmin, double xmax, const Colors& col)
{
  if (!parent) return;
  parent->cd();
  parent->Clear();

  // Subpads
  // top for distributions, bottom for ratio
  TPad* pTop = new TPad("pTop","", 0.00, 0.35, 1.00, 0.98);
  TPad* pBot = new TPad("pBot","", 0.00, 0.08, 1.00, 0.32);
  pTop->SetBottomMargin(0.04);
  pTop->SetTopMargin(0.12);
  pTop->SetLeftMargin(0.16);
  pTop->SetRightMargin(0.04);

  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.40);
  pBot->SetLeftMargin(0.16);
  pBot->SetRightMargin(0.04);

  pTop->SetTicks(1,1); pBot->SetTicks(1,1);
  pTop->Draw(); pBot->Draw();

  // --------- TOP: overlay ---------
  pTop->cd();

  // Style: DATA black markers+line; MC red line
  if (hD) { hD->SetLineColor(col.data); hD->SetMarkerColor(col.data); hD->SetMarkerStyle(20); hD->SetMarkerSize(0.7); hD->SetLineWidth(2); }
  if (hM) { hM->SetLineColor(col.mc);   hM->SetLineWidth(2); }

  // Axis frame
  TH1D* frame = (TH1D*)hD->Clone("frame_overlay");
  frame->Reset();
  frame->SetTitle(Form("%s;%s;Normalized entries", title.c_str(), xtitle.c_str()));
  frame->GetXaxis()->SetRangeUser(xmin, xmax);

  // y-range auto from max of both
  double ymax = 1.0;
  if (hD) ymax = std::max(ymax, hD->GetMaximum());
  if (hM) ymax = std::max(ymax, hM->GetMaximum());
  frame->SetMinimum(0.0);
  frame->SetMaximum(ymax * 1.25);
  frame->Draw("AXIS");

  if (hM) hM->Draw("HIST SAME");
  if (hD) hD->Draw("E SAME");

  // Legend (UR inside)
  auto L = new TLegend(0.62, 0.72, 0.94, 0.90);
  L->SetBorderSize(1);
  L->SetFillStyle(1001);
  L->SetFillColor(kWhite);
  L->SetTextSize(0.030);
  L->AddEntry(hD, "Data", "lep");
  L->AddEntry(hM, "MC",   "l");
  L->Draw();

  // --------- BOTTOM: ratio (Data / MC) ---------
  pBot->cd();
  TH1D* r = nullptr;
  if (hD && hM) {
    r = (TH1D*)hD->Clone("ratio");
    // Prevent div-by-zero: set bins where MC=0 to 0 with zero error
    for (int b=1;b<=r->GetNbinsX();++b) {
      double mc = hM->GetBinContent(b);
      if (mc<=0) { r->SetBinContent(b, 0.0); r->SetBinError(b, 0.0); }
    }
    r->Divide(hM); // uses stored errors (propagation)
    r->SetLineColor(kBlack);
    r->SetMarkerColor(kBlack);
    r->SetMarkerStyle(20);
    r->SetMarkerSize(0.6);
  }

  TH1D* fr = (TH1D*)hD->Clone("frame_ratio");
  fr->Reset();
  fr->SetTitle(Form(";%s;Data/MC", xtitle.c_str()));
  fr->GetXaxis()->SetRangeUser(xmin, xmax);
  fr->SetMinimum(0.5);
  fr->SetMaximum(1.5);
  fr->Draw("AXIS");
  if (r) r->Draw("E SAME");

  // y=1 line
  TLine* l1 = new TLine(xmin, 1.0, xmax, 1.0);
  l1->SetLineStyle(2);
  l1->SetLineWidth(1);
  l1->SetLineColor(kBlack);
  l1->Draw("SAME");
}

// ------------- one canvas for a given (property, xbin) -------------

static void draw_canvas_for_bin(
  TTree* tD, TTree* tM, const BranchNames& bn,
  const std::string& prop, int xbin)
{
  auto vars = var_list();
  const auto xe = x_edges();

  // Canvas: 4 columns x 2 rows = 8 cells (7 variables + 1 blank)
  TCanvas* c = new TCanvas(
    Form("c_cmp_%s_xbin%d", prop.c_str(), xbin),
    Form("DATA vs MC | %s | x_{B} in [%.2f,%.2f)", prop.c_str(), xe[xbin], xe[xbin+1]),
    1600, 900);
  c->Divide(4,2);

  gStyle->SetOptStat(0);

  Colors colors; colors.data = kBlack; colors.mc = kRed;

  for (size_t iv=0; iv<vars.size(); ++iv) {
    c->cd((int)iv+1);
    TPad* cell = (TPad*)gPad;

    auto hD = make_hist_for(tD, bn, prop, xbin, vars[iv], xe);
    auto hM = make_hist_for(tM, bn, prop, xbin, vars[iv], xe);

    draw_one_cell(cell, hD.get(), hM.get(),
                  vars[iv].title, vars[iv].xtitle,
                  vars[iv].xmin, vars[iv].xmax, colors);
  }

  // leave last pad (8th) blank
  c->cd(8); gPad->Clear();

  ensure_dir("output/compare");
  std::string out = Form("output/compare/compare_%s_xbin%d.pdf", prop.c_str(), xbin);
  c->SaveAs(out.c_str());
  std::cout << "Saved: " << out << "\n";

  delete c;
}

// ------------- main -------------

int main(int argc, char** argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) return 1;

  std::unique_ptr<TFile> fD(TFile::Open(cfg.data_path.c_str(), "READ"));
  std::unique_ptr<TFile> fM(TFile::Open(cfg.mc_path.c_str(),   "READ"));
  if (!fD || fD->IsZombie() || !fM || fM->IsZombie()) {
    std::cerr << "ERROR: could not open input files.\n";
    return 2;
  }

  TTree* tD = (TTree*)fD->Get("PhysicsEvents");
  TTree* tM = (TTree*)fM->Get("PhysicsEvents");
  if (!tD || !tM) {
    std::cerr << "ERROR: missing PhysicsEvents tree in one or both files.\n";
    return 3;
  }

  if (cfg.list_branches) {
    list_tree_branches(tD, "DATA");
    list_tree_branches(tM, "MC");
  }

  auto props = properties_to_run(cfg.which_property);
  const auto xe = x_edges();
  for (const auto& p : props) {
    for (int xb=0; xb<(int)xe.size()-1; ++xb) {
      draw_canvas_for_bin(tD, tM, cfg.bn, p, xb);
    }
  }
  return 0;
}