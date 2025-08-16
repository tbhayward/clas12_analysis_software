// compare_data_mc.cpp
// Compare reconstructed DATA vs reconstructed MC (REC) in xB bins and |t| windows.
// For each property (enpi, enpiLowt, enpiMidt, enpiHight), produce a multi-page PDF
// with 4 pages (one per xB bin). Each page is a 3x3 grid of variables; each cell
// contains an upper distribution overlay (DATA points, MC histogram) and a lower
// Data/MC ratio panel.
//
// Variables (7):
//  1) e_p       (electron momentum)        default branch "e_p"
//  2) e_theta   (electron polar angle)     "e_theta"
//  3) p_p       (pi+ momentum)             "p_p"
//  4) p_theta   (pi+ polar angle)          "p_theta"
//  5) Q2        (GeV^2)                    "Q2"
//  6) xB        (uses branch name in Config.bn.x, default "x")
//  7) -t        (absolute of branch Config.bn.t)
//
// Distributions use common binning for DATA and MC; bottom pad shows Data/MC ratio.
// MC is shape-normalized to Data per (variable, xBin, property) by default
// (disable with --no-norm).
//
// Cuts (as in your unfolding code):
//  - Mx2 in (0.80, 1.00)
//  - |t| in [0.10, 1.20]  (plus the property-specific sub-window)
//  - fiducial_status >= 100 (unless --no-fid-cut)
//  - xB in the 4 fixed bins: [0.10,0.30), [0.30,0.40), [0.40,0.50), [0.50,0.60)
//
// Output: output/enpi+/compare_<property>.pdf (multi-page).
//
// Build:
//   g++ -O2 -std=c++17 compare_data_mc.cpp `root-config --cflags --libs` -o compare_data_mc
//
// Run:
//   ./compare_data_mc data.root rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]
//       [--x x] [--t t] [--mx2 Mx2] [--fid fiducial_status]
//       [--ep e_p] [--etheta e_theta] [--pp p_p] [--ptheta p_theta] [--Q2 Q2]
//       [--nbins 50] [--no-fid-cut] [--no-norm] [--ratio 0.5 1.5]
//
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
#include "TLatex.h"

// ------------------ Config / CLI ------------------
struct BranchNames {
  std::string x   = "x";
  std::string t   = "t";
  std::string mx2 = "Mx2";
  std::string fid = "fiducial_status";
  // extra kinematics
  std::string e_p      = "e_p";
  std::string e_theta  = "e_theta";
  std::string p_p      = "p_p";
  std::string p_theta  = "p_theta";
  std::string Q2       = "Q2";
};

struct Config {
  std::string data_path;
  std::string rec_path;
  std::string which_property = "all";
  BranchNames bn;
  bool apply_fid_cut = true;
  bool do_shape_norm = true;
  int  nbins = 50;
  double ratio_min = 0.5, ratio_max = 1.5;
};

static void ensure_dir(const std::string& d) {
  if (gSystem->AccessPathName(d.c_str())) gSystem->mkdir(d.c_str(), kTRUE);
}

static void print_usage(const char* argv0) {
  std::cout << "Usage:\n  " << argv0
            << " data.root rec.root [enpi|enpiLowt|enpiMidt|enpiHight|all]\n"
            << "    [--x x] [--t t] [--mx2 Mx2] [--fid fiducial_status]\n"
            << "    [--ep e_p] [--etheta e_theta] [--pp p_p] [--ptheta p_theta] [--Q2 Q2]\n"
            << "    [--nbins 50] [--no-fid-cut] [--no-norm] [--ratio 0.5 1.5]\n";
}

static bool parse_args(int argc, char** argv, Config& cfg) {
  if (argc < 4) { print_usage(argv[0]); return false; }
  cfg.data_path = argv[1];
  cfg.rec_path  = argv[2];
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
    auto eat_dbl2 = [&](const char* key, double& a1, double& a2) {
      if (a == key && i+2 < argc) { a1 = std::atof(argv[++i]); a2 = std::atof(argv[++i]); return true; }
      return false;
    };
    if (eat("--x", cfg.bn.x)) continue;
    if (eat("--t", cfg.bn.t)) continue;
    if (eat("--mx2", cfg.bn.mx2)) continue;
    if (eat("--fid", cfg.bn.fid)) continue;
    if (eat("--ep", cfg.bn.e_p)) continue;
    if (eat("--etheta", cfg.bn.e_theta)) continue;
    if (eat("--pp", cfg.bn.p_p)) continue;
    if (eat("--ptheta", cfg.bn.p_theta)) continue;
    if (eat("--Q2", cfg.bn.Q2)) continue;
    if (eat_int("--nbins", cfg.nbins)) continue;
    if (a == "--no-fid-cut") { cfg.apply_fid_cut=false; continue; }
    if (a == "--no-norm") { cfg.do_shape_norm=false; continue; }
    if (eat_dbl2("--ratio", cfg.ratio_min, cfg.ratio_max)) continue;
    std::cerr << "Unknown arg: " << a << "\n";
  }
  return true;
}

// ------------------ Binning & property windows ------------------
static std::vector<double> x_edges() { return {0.10,0.30,0.40,0.50,0.60}; }
static int xbin_index(double x, const std::vector<double>& e) {
  for (size_t i=0;i+1<e.size();++i) if (x >= e[i] && x < e[i+1]) return (int)i;
  return -1;
}
static inline bool in_common_cuts(double t, double mx2, int fid, bool has_fid, bool apply_fid) {
  const double at = std::fabs(t);
  if (!(at >= 0.10 && at <= 1.20)) return false;
  if (!(mx2 > 0.80 && mx2 < 1.00)) return false;
  if (apply_fid) { if (!has_fid || fid < 100) return false; }
  return true;
}
static inline bool pass_prop_window(const std::string& prop, double t) {
  const double at = std::fabs(t);
  if (prop == "enpi")      return (at >= 0.10 && at <= 1.20);
  if (prop == "enpiLowt")  return (at >= 0.10 && at <= 0.4667);
  if (prop == "enpiMidt")  return (at >= 0.4667 && at <= 0.8333);
  if (prop == "enpiHight") return (at >= 0.8333 && at <= 1.20);
  return false;
}
static std::vector<std::string> all_props() {
  return {"enpi","enpiLowt","enpiMidt","enpiHight"};
}
static std::vector<std::string> props_to_run(const std::string& which) {
  if (which=="all") return all_props();
  return {which};
}
static std::string prop_title(const std::string& p) {
  if (p == "enpi")      return "ep #rightarrow e' n #pi^{+}";
  if (p == "enpiLowt")  return "low |t|";
  if (p == "enpiMidt")  return "mid |t|";
  if (p == "enpiHight") return "high |t|";
  return p;
}
static std::string prop_tlabel(const std::string& prop) {
  if (prop == "enpi")      return "-t #in [0.10, 1.20]";
  if (prop == "enpiLowt")  return "-t #in [0.10, 0.4667]";
  if (prop == "enpiMidt")  return "-t #in [0.4667, 0.8333]";
  if (prop == "enpiHight") return "-t #in [0.8333, 1.20]";
  return "";
}

// ------------------ Branch binding ------------------
struct BranchHandles {
  // basic
  double x=0, t=0, mx2=0; int fid=0;
  bool has_fid=false;
  // extras
  double e_p=0, e_theta=0, p_p=0, p_theta=0, Q2=0;
  bool has_e_p=false, has_e_theta=false, has_p_p=false, has_p_theta=false, has_Q2=false;
};

static bool bind_tree(TTree* tr, const BranchNames& bn, BranchHandles& b) {
  if (!tr) return false;
  auto has = [&](const std::string& name){ return tr->GetBranch(name.c_str()) != nullptr; };
  if (!has(bn.x) || !has(bn.t) || !has(bn.mx2)) {
    std::cerr << "ERROR: missing required branches x/t/Mx2 on tree " << tr->GetName() << "\n";
    return false;
  }
  tr->SetBranchAddress(bn.x.c_str(), &b.x);
  tr->SetBranchAddress(bn.t.c_str(), &b.t);
  tr->SetBranchAddress(bn.mx2.c_str(), &b.mx2);
  b.has_fid = has(bn.fid);
  if (b.has_fid) tr->SetBranchAddress(bn.fid.c_str(), &b.fid);

  b.has_e_p     = has(bn.e_p);     if (b.has_e_p)     tr->SetBranchAddress(bn.e_p.c_str(), &b.e_p);
  b.has_e_theta = has(bn.e_theta); if (b.has_e_theta) tr->SetBranchAddress(bn.e_theta.c_str(), &b.e_theta);
  b.has_p_p     = has(bn.p_p);     if (b.has_p_p)     tr->SetBranchAddress(bn.p_p.c_str(), &b.p_p);
  b.has_p_theta = has(bn.p_theta); if (b.has_p_theta) tr->SetBranchAddress(bn.p_theta.c_str(), &b.p_theta);
  b.has_Q2      = has(bn.Q2);      if (b.has_Q2)      tr->SetBranchAddress(bn.Q2.c_str(), &b.Q2);
  return true;
}

// ------------------ Variable descriptors ------------------
struct VarDesc {
  std::string key;        // internal id
  std::string xtitle;     // axis label
  int nbins;
  double vmin, vmax;
  bool available = true;  // set false if branch missing
};

static std::vector<VarDesc> make_var_descs(const BranchHandles& b, int nbins) {
  std::vector<VarDesc> V;

  // Tuned, conservative default ranges; adjust here if your kinematics differ.
  V.push_back({"e_p",     "e momentum p_{e} (GeV)",            nbins, 0.0, 12.0,  b.has_e_p});
  V.push_back({"e_theta", "e polar angle #theta_{e} (deg)",    nbins, 0.0, 40.0,  b.has_e_theta});
  V.push_back({"p_p",     "#pi^{+} momentum p_{#pi} (GeV)",    nbins, 0.0, 6.0,   b.has_p_p});
  V.push_back({"p_theta", "#pi^{+} polar #theta_{#pi} (deg)",  nbins, 0.0, 90.0,  b.has_p_theta});
  V.push_back({"Q2",      "Q^{2} (GeV^{2})",                   nbins, 0.0, 10.0,  b.has_Q2});
  V.push_back({"xB",      "x_{B}",                             nbins, 0.10, 0.60,  true});
  V.push_back({"t",       "-t (GeV^{2})",                      nbins, 0.10, 1.20,  true});
  return V;
}

static inline double get_var_value(const VarDesc& vd, const BranchHandles& b) {
  if (vd.key=="e_p")     return b.e_p;
  if (vd.key=="e_theta") return b.e_theta;
  if (vd.key=="p_p")     return b.p_p;
  if (vd.key=="p_theta") return b.p_theta;
  if (vd.key=="Q2")      return b.Q2;
  if (vd.key=="xB")      return b.x;
  if (vd.key=="t")       return std::fabs(b.t);
  return 0.0;
}

// ------------------ Hist containers ------------------
struct HistPair { std::unique_ptr<TH1D> hD, hM; };
using VarHists = std::vector<HistPair>;              // one per variable
using XbinVarHists = std::vector<VarHists>;          // size = NX x Nvar

static XbinVarHists make_hists_for_property(const std::vector<VarDesc>& vars,
                                            int nx, const std::string& prop) {
  XbinVarHists hv(nx);
  for (int ib=0; ib<nx; ++ib) {
    hv[ib].reserve(vars.size());
    for (size_t iv=0; iv<vars.size(); ++iv) {
      const auto& v = vars[iv];
      // If unavailable, still create tiny placeholder to keep drawing flow simple
      auto hD = std::make_unique<TH1D>(Form("hD_%s_%s_ib%d", prop.c_str(), v.key.c_str(), ib),
                                       "", v.nbins, v.vmin, v.vmax);
      auto hM = std::make_unique<TH1D>(Form("hM_%s_%s_ib%d", prop.c_str(), v.key.c_str(), ib),
                                       "", v.nbins, v.vmin, v.vmax);
      hD->Sumw2(); hM->Sumw2();
      hv[ib].push_back({std::move(hD), std::move(hM)});
    }
  }
  return hv;
}

// ------------------ Fill loops ------------------
static void fill_tree_into_hists(TTree* tr, const BranchNames& bn, const std::string& prop,
                                 bool apply_fid, const std::vector<VarDesc>& vars,
                                 XbinVarHists& hv, bool isData)
{
  BranchHandles b;
  if (!bind_tree(tr, bn, b)) return;

  const auto xe = x_edges();
  const Long64_t nent = tr->GetEntries();
  for (Long64_t i=0; i<nent; ++i) {
    tr->GetEntry(i);
    if (!in_common_cuts(b.t, b.mx2, b.fid, b.has_fid, apply_fid)) continue;
    if (!pass_prop_window(prop, b.t)) continue;
    int ib = xbin_index(b.x, xe);
    if (ib < 0) continue;

    for (size_t iv=0; iv<vars.size(); ++iv) {
      const auto& v = vars[iv];
      if (!v.available) continue;
      const double val = get_var_value(v, b);
      if (std::isnan(val)) continue;
      if (isData) hv[ib][iv].hD->Fill(val);
      else        hv[ib][iv].hM->Fill(val);
    }
  }
}

// ------------------ Drawing helpers ------------------
struct Colors {
  int dataMarker=20, dataColor=kBlack;
  int mcLineColor=kAzure+2, mcFillColor=0, mcLineWidth=2;
};

static void draw_dist_ratio_cell(TH1D* hD, TH1D* hM, const std::string& xtitle,
                                 const std::string& titleTop,
                                 double ratioMin, double ratioMax,
                                 bool do_norm, const Colors& C)
{
  // Split pad
  TPad* pTop = new TPad("pTop","",0,0.32,1,1);  pTop->SetNumber(1);
  TPad* pBot = new TPad("pBot","",0,0,1,0.32);  pBot->SetNumber(2);
  pTop->SetBottomMargin(0.02); pTop->SetTopMargin(0.08);
  pTop->SetLeftMargin(0.18);   pTop->SetRightMargin(0.06);
  pBot->SetTopMargin(0.02);    pBot->SetBottomMargin(0.38);
  pBot->SetLeftMargin(0.18);   pBot->SetRightMargin(0.06);
  pTop->Draw(); pBot->Draw();

  // normalize MC to Data if requested
  double sf = 1.0;
  if (do_norm) {
    double iD = hD->Integral(0, hD->GetNbinsX()+1);
    double iM = hM->Integral(0, hM->GetNbinsX()+1);
    if (iM>0) { sf = iD / iM; hM->Scale(sf); }
  }

  // --- TOP: overlay ---
  pTop->cd();
  hM->SetLineColor(C.mcLineColor);
  hM->SetLineWidth(C.mcLineWidth);
  hM->SetFillStyle(0);
  hM->SetTitle(titleTop.c_str());
  hM->GetXaxis()->SetLabelSize(0); // hide x labels on top
  hM->GetYaxis()->SetTitle("Events");
  hM->GetYaxis()->SetTitleSize(0.06);
  hM->GetYaxis()->SetLabelSize(0.05);
  hM->GetYaxis()->SetTitleOffset(1.25);

  // autoscale y by union
  double ymax = std::max(hM->GetMaximum(), hD->GetMaximum());
  hM->SetMaximum(1.25*ymax);
  hM->Draw("HIST");

  auto gD = (TH1D*)hD->Clone(Form("%s_pts", hD->GetName()));
  gD->SetMarkerStyle(C.dataMarker);
  gD->SetMarkerColor(C.dataColor);
  gD->SetLineColor(C.dataColor);
  gD->Draw("E1 SAME");

  // legend (UR)
  auto L = new TLegend(0.58,0.72,0.94,0.92);
  L->SetBorderSize(1); L->SetFillStyle(1001); L->SetFillColor(kWhite); L->SetTextSize(0.045);
  L->AddEntry(gD,"DATA","lep");
  L->AddEntry(hM, (do_norm?Form("MC (norm x %.3g)", sf):"MC").c_str(), "l");
  L->Draw();

  // --- BOTTOM: ratio ---
  pBot->cd();
  TH1D* r = (TH1D*)hD->Clone(Form("%s_ratio", hD->GetName()));
  r->Divide(hM);
  r->SetTitle("");
  r->GetYaxis()->SetTitle("Data/MC");
  r->GetXaxis()->SetTitle(xtitle.c_str());
  r->GetXaxis()->SetTitleSize(0.13);
  r->GetYaxis()->SetTitleSize(0.12);
  r->GetXaxis()->SetLabelSize(0.12);
  r->GetYaxis()->SetLabelSize(0.11);
  r->GetYaxis()->SetTitleOffset(0.50);
  r->GetYaxis()->SetNdivisions(505);
  r->SetMarkerStyle(C.dataMarker);
  r->SetMarkerColor(C.dataColor);
  r->SetLineColor(C.dataColor);
  r->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  r->Draw("E1");

  TLine* l1 = new TLine(r->GetXaxis()->GetXmin(),1.0,r->GetXaxis()->GetXmax(),1.0);
  l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw("SAME");
}

// ------------------ One property -> multipage PDF ------------------
static void draw_property_pdf(const std::string& prop,
                              const std::vector<VarDesc>& vars,
                              const XbinVarHists& hv,
                              const Config& cfg)
{
  ensure_dir("output/enpi+");
  std::string pdf = "output/enpi+/compare_" + prop + ".pdf";

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas(("c_"+prop).c_str(), ("Data vs MC: "+prop).c_str(), 1800, 1200);
  c->Divide(3,3);

  // open multipage
  c->Print((pdf + "[").c_str());

  auto xe = x_edges();
  const int NX = (int)xe.size()-1;
  Colors colors;

  for (int ib=0; ib<NX; ++ib) {
    c->cd(0); // ensure clean
    // draw 7 variables into first 7 pads
    int pad = 1;
    for (size_t iv=0; iv<vars.size(); ++iv) {
      c->cd(pad++);
      if (!vars[iv].available) {
        // nice placeholder
        gPad->Clear();
        TLatex t; t.SetNDC(); t.SetTextSize(0.05);
        t.DrawLatex(0.20,0.55,Form("%s unavailable on tree", vars[iv].key.c_str()));
        continue;
      }
      // compose small title with xB, t-window
      std::string ttl = Form("x_{B} #in [%.2f, %.2f),   %s",
                             xe[ib], xe[ib+1], prop_tlabel(prop).c_str());
      draw_dist_ratio_cell(hv[ib][iv].hD.get(), hv[ib][iv].hM.get(),
                           vars[iv].xtitle, ttl,
                           cfg.ratio_min, cfg.ratio_max,
                           cfg.do_shape_norm, colors);
    }
    // leave pad 8 and 9 blank or label
    for (; pad<=9; ++pad) {
      c->cd(pad);
      gPad->Clear();
      if (pad==9) {
        TLatex t; t.SetNDC(); t.SetTextSize(0.060); t.SetTextAlign(22);
        t.DrawLatex(0.50,0.55,Form("%s — %s", prop_title(prop).c_str(), prop_tlabel(prop).c_str()));
        TLatex xb; xb.SetNDC(); xb.SetTextSize(0.050); xb.SetTextAlign(22);
        xb.DrawLatex(0.50,0.42,Form("x_{B} bin: [%.2f, %.2f)", xe[ib], xe[ib+1]));
      }
    }
    c->Print(pdf.c_str());
  }

  // close multipage
  c->Print((pdf + "]").c_str());
  std::cout << "Saved: " << pdf << "\n";
  delete c;
}

// ------------------ main ------------------
int main(int argc, char** argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) return 1;

  std::unique_ptr<TFile> fD(TFile::Open(cfg.data_path.c_str(), "READ"));
  std::unique_ptr<TFile> fM(TFile::Open(cfg.rec_path.c_str(),  "READ")); // REC MC
  if (!fD || fD->IsZombie() || !fM || fM->IsZombie()) {
    std::cerr << "ERROR: could not open input files.\n";
    return 2;
  }
  TTree* tD = (TTree*)fD->Get("PhysicsEvents");
  TTree* tM = (TTree*)fM->Get("PhysicsEvents");
  if (!tD || !tM) {
    std::cerr << "ERROR: missing 'PhysicsEvents' tree in one or both files.\n";
    return 3;
  }

  // Probe branches to decide availability / set ranges
  BranchHandles probe;
  if (!bind_tree(tD, cfg.bn, probe)) return 4; // use DATA to probe
  auto varDescs = make_var_descs(probe, cfg.nbins);

  // For each property, build hist sets and fill
  const auto props = props_to_run(cfg.which_property);
  const int NX = (int)x_edges().size()-1;

  for (const auto& prop : props) {
    auto hv = make_hists_for_property(varDescs, NX, prop);
    // fill DATA
    fill_tree_into_hists(tD, cfg.bn, prop, cfg.apply_fid_cut, varDescs, hv, /*isData=*/true);
    // fill MC (REC)
    fill_tree_into_hists(tM, cfg.bn, prop, cfg.apply_fid_cut, varDescs, hv, /*isData=*/false);

    // style (global)
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");

    // draw to one multi-page PDF per property
    draw_property_pdf(prop, varDescs, hv, cfg);
  }

  return 0;
}