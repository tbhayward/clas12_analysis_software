// compare_data_mc.cpp
// DATA vs MC comparisons with distributions and ratio (DATA/MC)
// - Uses your xB and |t| binning from the unfolding script
// - One-pass fill per (tree, property): fast
// - DATA shown in black, MC in red
// - Angles in the tree are in radians; plotted in degrees
// - Each histogram is normalized to its own integral (shape comparison)
// - Per property and per xB bin, produce one 4x4 canvas:
//     * Top 2 rows (pads 1..8): 7 variables (1 blank pad reserved for legend)
//     * Bottom 2 rows (pads 9..16): corresponding 7 ratio pads (DATA/MC), aligned
//
// Variables & ranges:
//   e_p      : [2, 9] GeV
//   e_theta  : [0, 40] deg  (input radians -> deg)
//   p_p      : [1, 4] GeV
//   p_theta  : [0, 80] deg  (input radians -> deg)
//   Q2       : [1, 8] GeV^2
//   xB       : [0, 1]
//   -t       : [0, 1.2]
//
// Build:
//   g++ -O2 -std=c++17 compare_data_mc.cpp `root-config --cflags --libs` -o compare_data_mc
//
// Run:
//   ./compare_data_mc data.root mc.root [enpi|enpiLowt|enpiMidt|enpiHight|all]
//     [--x x] [--t t] [--mx2 Mx2] [--fid fiducial_status]
//     [--ep ep] [--eth e_theta] [--pp pp] [--pth p_theta] [--Q2 Q2]
//     [--list-branches] [--debug N]
//
// Output:
//   output/compare/compare_<property>_xbin<i>.pdf  (i=0..3)
//
// Notes:
// - Cuts applied: fid>=100 (if available), Mx2 in (0.80,1.00), |t| per property window,
//                 xB in the specified x bin on each canvas.
// - Ratios use normalized shapes: DATA_norm / MC_norm, with a dashed y=1 line.
// - If MC is zero in a bin, ratio content & error are set to 0 in that bin.

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
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

// ------------------ CLI/Config ------------------

struct BranchNames {
  std::string x   = "x";
  std::string t   = "t";     // signed
  std::string mx2 = "Mx2";
  std::string fid = "fiducial_status";

  // Kinematics for comparison (defaults; override via CLI if needed)
  std::string ep  = "ep";    // electron momentum [GeV]
  std::string eth = "eth";   // electron polar angle [rad]
  std::string pp  = "pp";    // proton momentum [GeV]
  std::string pth = "pth";   // proton polar angle [rad]
  std::string Q2  = "Q2";    // Q^2 [GeV^2]
};

struct Config {
  std::string data_path, mc_path;
  std::string which_property = "all";
  BranchNames bn;
  bool list_branches = false;
  int  debug_print   = 0;
};

static void print_usage(const char* argv0) {
  std::cout << "Usage:\n  " << argv0
            << " data.root mc.root [enpi|enpiLowt|enpiMidt|enpiHight|all]\n"
            << "    [--x x] [--t t] [--mx2 Mx2] [--fid fiducial_status]\n"
            << "    [--ep ep] [--eth e_theta] [--pp pp] [--pth p_theta] [--Q2 Q2]\n"
            << "    [--list-branches] [--debug N]\n";
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
    if (eat("--t",   cfg.bn.t)) continue;
    if (eat("--mx2", cfg.bn.mx2)) continue;
    if (eat("--fid", cfg.bn.fid)) continue;
    if (eat("--ep",  cfg.bn.ep)) continue;
    if (eat("--eth", cfg.bn.eth)) continue;
    if (eat("--pp",  cfg.bn.pp)) continue;
    if (eat("--pth", cfg.bn.pth)) continue;
    if (eat("--Q2",  cfg.bn.Q2)) continue;

    if (a == "--list-branches") { cfg.list_branches = true; continue; }
    if (eat_int("--debug", cfg.debug_print)) continue;

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

struct Handles {
  double x=0, t=0, mx2=0, ep=0, eth=0, pp=0, pth=0, Q2=0;
  int fid=0;
  bool has_fid=false;
};

static bool bind_tree(TTree* tr, const BranchNames& bn, Handles& h, const char* who) {
  if (!tr) return false;
  bool ok = true;
  ok &= (tr->GetBranch(bn.x.c_str())  != nullptr);
  ok &= (tr->GetBranch(bn.t.c_str())  != nullptr);
  ok &= (tr->GetBranch(bn.mx2.c_str())!= nullptr);
  ok &= (tr->GetBranch(bn.ep.c_str()) != nullptr);
  ok &= (tr->GetBranch(bn.eth.c_str())!= nullptr);
  ok &= (tr->GetBranch(bn.pp.c_str()) != nullptr);
  ok &= (tr->GetBranch(bn.pth.c_str())!= nullptr);
  ok &= (tr->GetBranch(bn.Q2.c_str()) != nullptr);
  if (!ok) {
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

// ------------------ variable defs ------------------

struct VarDef {
  std::string name;
  std::string title;   // for pad title
  std::string xtitle;  // axis
  int    nbins;
  double xmin, xmax;
};

static std::vector<VarDef> var_list() {
  return {
    {"e_p",     "e_{p}",             "e_{p}  [GeV]",         70,  2.0,  9.0},
    {"e_theta", "#theta_{e} (deg)",  "#theta_{e}  [deg]",    80,  0.0, 40.0},
    {"p_p",     "p_{p}",             "p_{p}  [GeV]",         60,  1.0,  4.0},
    {"p_theta", "#theta_{p} (deg)",  "#theta_{p}  [deg]",    80,  0.0, 80.0},
    {"Q2",      "Q^{2}",             "Q^{2}  [GeV^{2}]",     70,  1.0,  8.0},
    {"xB",      "x_{B}",             "x_{B}",                 50,  0.0,  1.0},
    {"mt",      "-t",                "-t  [GeV^{2}]",        60,  0.0,  1.2}
  };
}

// ------------------ histogram grid & fast one-pass fill ------------------

using HistGrid = std::vector<std::vector<std::unique_ptr<TH1D>>>; // [ivar][xbin]

static HistGrid build_hists_one_pass(
  TTree* tr, const BranchNames& bn, const std::string& prop,
  const std::vector<VarDef>& vars, const std::vector<double>& xedges)
{
  const int NV = (int)vars.size();
  const int NX = (int)xedges.size()-1;

  HistGrid H(NV, std::vector<std::unique_ptr<TH1D>>(NX));
  for (int v=0; v<NV; ++v) {
    for (int xb=0; xb<NX; ++xb) {
      H[v][xb] = std::make_unique<TH1D>(
        Form("h_%s_%s_xb%d", vars[v].name.c_str(), prop.c_str(), xb),
        "", vars[v].nbins, vars[v].xmin, vars[v].xmax);
      H[v][xb]->Sumw2();
    }
  }

  Handles h; if (!bind_tree(tr, bn, h, (prop+"(FAST)").c_str())) return H;

  const Long64_t N = tr->GetEntries();
  for (Long64_t i=0;i<N;i++) {
    tr->GetEntry(i);

    // Common selection
    if (h.has_fid && h.fid < 100) continue;
    if (!(h.mx2 > 0.80 && h.mx2 < 1.00)) continue;
    if (!pass_prop_window(prop, h.t)) continue;

    int xb = xbin_index(h.x, xedges);
    if (xb < 0 || xb >= NX) continue;

    // Derived quantities
    const double e_theta_deg = h.eth * 180.0 / TMath::Pi();
    const double p_theta_deg = h.pth * 180.0 / TMath::Pi();
    const double tabs        = std::fabs(h.t);

    // Fill all variables (respect plotting ranges)
    for (int v=0; v<NV; ++v) {
      double val = 0.0;
      if      (vars[v].name == "e_p")     val = h.ep;
      else if (vars[v].name == "e_theta") val = e_theta_deg;
      else if (vars[v].name == "p_p")     val = h.pp;
      else if (vars[v].name == "p_theta") val = p_theta_deg;
      else if (vars[v].name == "Q2")      val = h.Q2;
      else if (vars[v].name == "xB")      val = h.x;
      else if (vars[v].name == "mt")      val = tabs;

      if (val >= vars[v].xmin && val < vars[v].xmax) H[v][xb]->Fill(val);
    }
  }

  // Normalize each histogram to its own integral (shape-only)
  for (int v=0; v<NV; ++v)
    for (int xb=0; xb<NX; ++xb) {
      double I = H[v][xb]->Integral();
      if (I > 0) H[v][xb]->Scale(1.0/I);
    }

  return H;
}

// Safe ratio = data/mc (per-bin; handle mc==0)
static std::unique_ptr<TH1D> make_ratio(const TH1D* hD, const TH1D* hM, const char* name) {
  auto r = std::make_unique<TH1D>(*hD);
  r->SetName(name);
  const int nb = r->GetNbinsX();
  for (int b=1;b<=nb;++b) {
    double d  = hD->GetBinContent(b);
    double dd = hD->GetBinError(b);
    double m  = hM->GetBinContent(b);
    double md = hM->GetBinError(b);
    if (m <= 0) { r->SetBinContent(b, 0.0); r->SetBinError(b, 0.0); continue; }
    double val = d/m;
    // simple error propagation: (d/m)*sqrt( (dd/d)^2 + (md/m)^2 )
    double err = 0.0;
    if (d>0 && m>0) err = val*std::sqrt( (dd>0? (dd*dd)/(d*d):0.0) + (md>0? (md*md)/(m*m):0.0) );
    r->SetBinContent(b, val);
    r->SetBinError(b, err);
  }
  return r;
}

// ------------------ drawing helpers ------------------

struct Colors { int data=kBlack; int mc=kRed; };

static void style_axes(TH1* h, const char* xtitle, const char* ytitle, double ytitle_offset=1.45) {
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitleOffset(ytitle_offset);
}

static void draw_distributions_cell(TH1D* hD, TH1D* hM,
                                    const std::string& title,
                                    const std::string& xtitle,
                                    double xmin, double xmax,
                                    const Colors& colors)
{
  if (!hD || !hM) return;

  // Style
  hD->SetLineColor(colors.data); hD->SetLineWidth(3);
  hM->SetLineColor(colors.mc);   hM->SetLineWidth(3);
  hD->SetStats(0); hM->SetStats(0);

  // Axes/titles
  style_axes(hD, xtitle.c_str(), "Normalized counts");
  hD->SetTitle(title.c_str());
  hD->GetXaxis()->SetLimits(xmin, xmax);

  // Y range: a bit of headroom over both
  double ymax = std::max(hD->GetMaximum(), hM->GetMaximum());
  if (ymax <= 0) ymax = 1.0;
  hD->SetMinimum(0.0);
  hD->SetMaximum(1.15*ymax);

  hD->Draw("HIST");
  hM->Draw("HIST SAME");
}

static void draw_ratio_cell(TH1D* r, const std::string& xtitle, double xmin, double xmax) {
  if (!r) return;
  style_axes(r, xtitle.c_str(), "DATA / MC", 1.55);
  r->SetLineColor(kBlack); r->SetMarkerColor(kBlack);
  r->SetMarkerStyle(20); r->SetMarkerSize(0.8);
  r->SetLineWidth(2);

  r->GetXaxis()->SetLimits(xmin, xmax);
  r->SetMinimum(0.5);
  r->SetMaximum(1.5);
  r->Draw("E1");

  // unity line
  TLine* l1 = new TLine(xmin, 1.0, xmax, 1.0);
  l1->SetLineStyle(2);
  l1->SetLineColor(kGray+2);
  l1->Draw("SAME");
}

static void draw_global_legend() {
  auto L = new TLegend(0.68, 0.74, 0.93, 0.92);
  L->SetBorderSize(1);
  L->SetFillStyle(1001);
  L->SetFillColor(kWhite);
  L->SetTextSize(0.035);
  // dummy lines
  auto hD = new TH1D("hDummyD","",1,0,1);
  auto hM = new TH1D("hDummyM","",1,0,1);
  hD->SetLineColor(kBlack); hD->SetLineWidth(3);
  hM->SetLineColor(kRed);   hM->SetLineWidth(3);
  L->AddEntry(hD, "DATA (norm)", "l");
  L->AddEntry(hM, "MC (norm)",   "l");
  L->Draw();
}

// Map variable index -> pad index (4x4 grid)
// Dist pads in {1..8}, ratio pads = distPad + 8
static std::vector<int> dist_pad_positions() {
  // Fill across rows: 1..4 (row1), 5..8 (row2)
  // Use 7 slots, leave pad 8 for legend.
  return {1,2,3,4,5,6,7}; // pad 8 reserved
}

// ------------------ draw per property & x-bin ------------------

static void draw_canvas_for_xbin(
  const std::string& prop, int xbin,
  const HistGrid& HD, const HistGrid& HM,
  const std::vector<VarDef>& vars,
  const std::vector<double>& xedges)
{
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetTitleSize(0.050, "XYZ");
  gStyle->SetLabelSize(0.040, "XYZ");

  const double xl = xedges[xbin], xh = xedges[xbin+1];
  TCanvas* c = new TCanvas(
    Form("c_cmp_%s_xbin%d", prop.c_str(), xbin),
    Form("DATA vs MC | %s | x_{B} in [%.2f, %.2f),  %s",
         prop_title(prop).c_str(), xl, xh, prop_tlabel(prop).c_str()),
    1600, 1200);
  c->Divide(4,4);

  Colors colors; colors.data = kBlack; colors.mc = kRed;

  auto pads = dist_pad_positions(); // size 7, pad 8 kept for legend

  for (size_t iv=0; iv<vars.size(); ++iv) {
    int pDist = pads[iv];
    int pRatio = pDist + 8;

    // Draw distributions
    c->cd(pDist);
    {
      gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.10);
      TH1D* hD = HD[iv][xbin].get();
      TH1D* hM = HM[iv][xbin].get();
      draw_distributions_cell(hD, hM, vars[iv].title, vars[iv].xtitle,
                              vars[iv].xmin, vars[iv].xmax, colors);
    }

    // Draw ratio
    c->cd(pRatio);
    {
      gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.10);
      auto r = make_ratio(HD[iv][xbin].get(), HM[iv][xbin].get(),
                          Form("r_%s_%s_xb%d", prop.c_str(), vars[iv].name.c_str(), xbin).c_str());
      draw_ratio_cell(r.get(), vars[iv].xtitle, vars[iv].xmin, vars[iv].xmax);
    }
  }

  // Pad 8 (top-right of the distributions grid) reserved for a clean legend
  c->cd(8);
  gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.10); gPad->SetTopMargin(0.10);
  gPad->Clear();
  draw_global_legend();

  // Save
  ensure_dir("output/compare");
  std::string out = Form("output/compare/compare_%s_xbin%d.pdf", prop.c_str(), xbin);
  c->SaveAs(out.c_str());
  std::cout << "Saved: " << out << "\n";
  delete c;
}

// ------------------ Orchestration ------------------

int main(int argc, char** argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) return 1;

  std::unique_ptr<TFile> fD(TFile::Open(cfg.data_path.c_str(), "READ"));
  std::unique_ptr<TFile> fM(TFile::Open(cfg.mc_path.c_str(),   "READ"));
  if (!fD || fD->IsZombie() || !fM || fM->IsZombie()) {
    std::cerr << "ERROR: could not open one or both input files.\n";
    return 2;
  }

  // Heuristic tree name (same as unfolding script)
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

  const auto props = properties_to_run(cfg.which_property);
  const auto vars  = var_list();
  const auto xe    = x_edges();
  const int  NX    = (int)xe.size()-1;

  // Build once per (file, property) — fast
  for (const auto& p : props) {
    HistGrid HD = build_hists_one_pass(tD, cfg.bn, p, vars, xe);
    HistGrid HM = build_hists_one_pass(tM, cfg.bn, p, vars, xe);

    for (int xb=0; xb<NX; ++xb) {
      draw_canvas_for_xbin(p, xb, HD, HM, vars, xe);
    }
  }

  return 0;
}