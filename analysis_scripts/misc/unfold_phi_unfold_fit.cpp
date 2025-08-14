// unfold_phi_unfold_fit.cpp
// Build & run (csh):
//   g++ -O2 -std=c++17 unfold_phi_unfold_fit.cpp `root-config --cflags --libs` -o unfold_phi_unfold_fit
//   ./unfold_phi_unfold_fit data.root gen.root recmc.root [enpi|enpiLowt|enpiMidt|enpiHight|all] [--x x] [--phi phi] [--t t] [--mx2 Mx2] [--phibins 24] [--phi-deg]
// Output: output/enpi+/PROPERTY_unfolded.pdf

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

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
  std::string x  = "x";
  std::string phi = "phi";
  std::string t   = "t";
  std::string mx2 = "Mx2";
  std::string fid = "fiducial_status";
};

struct FitResult {
  double C=0, dC=0;
  double A=0, dA=0;
  double B=0, dB=0;
  double chi2=0; int ndf=0;
};

static inline bool in_range(double v, double lo, double hi) {
  return (v >= lo) && (v < hi);
}

static inline double wrap_phi_0_2pi(double phi) {
  const double twopi = 2.0 * TMath::Pi();
  phi = std::fmod(phi, twopi);
  if (phi < 0) phi += twopi;
  return phi;
}

struct Cuts {
  // physics cuts common to all four properties:
  std::cout << t << " " << mx2 << std::endl;
  static bool pass_common(double t, double mx2) {
    return (std::fabs(t) <= 1.0) && (std::fabs(t) >= 0.0) && (mx2 > 0.80) && (mx2 < 1.00);
  }
  static bool pass_property_window(const std::string& prop, double t) {
    if (prop == "enpi")       return true; // just the common cuts (including |t|<=1)
    if (prop == "enpiLowt")   return (std::fabs(t) >= 0.00 && std::fabs(t) <= 0.30);
    if (prop == "enpiMidt")   return (std::fabs(t) >= 0.30 && std::fabs(t) <= 0.70);
    if (prop == "enpiHight")  return (std::fabs(t) >= 0.70); // up to 1.0 bounded by common cut
    return false;
  }
};

struct BranchHandles {
  // Pointers bound directly to TTree branches
  double x=0, phi=0, t=0, mx2=0;
  int fid=0;
  bool has_fid=false;
};

static bool bind_branches(TTree* tr, const BranchNames& bn, BranchHandles& bh, bool want_fid) {
  if (!tr) return false;
  bool ok = true;
  ok &= (tr->GetBranch(bn.x.c_str())  != nullptr);
  ok &= (tr->GetBranch(bn.phi.c_str())!= nullptr);
  ok &= (tr->GetBranch(bn.t.c_str())  != nullptr);
  ok &= (tr->GetBranch(bn.mx2.c_str())!= nullptr);
  if (!ok) return false;

  tr->SetBranchAddress(bn.x.c_str(),   &bh.x);
  tr->SetBranchAddress(bn.phi.c_str(), &bh.phi);
  tr->SetBranchAddress(bn.t.c_str(),   &bh.t);
  tr->SetBranchAddress(bn.mx2.c_str(), &bh.mx2);

  bh.has_fid = (tr->GetBranch(bn.fid.c_str()) != nullptr);
  if (want_fid && bh.has_fid) tr->SetBranchAddress(bn.fid.c_str(), &bh.fid);
  return true;
}

static std::string property_title(const std::string& prop) {
  if (prop == "enpi")      return "ep -> e' n #pi^{+}";
  if (prop == "enpiLowt")  return "ep -> e' n #pi^{+}, low |t|";
  if (prop == "enpiMidt")  return "ep -> e' n #pi^{+}, mid |t|";
  if (prop == "enpiHight") return "ep -> e' n #pi^{+}, high |t|";
  return prop;
}

static void ensure_dir(const std::string& d) {
  if (gSystem->AccessPathName(d.c_str())) {
    gSystem->mkdir(d.c_str(), kTRUE);
  }
}

struct Config {
  std::string data_path;
  std::string gen_path;
  std::string rec_path;
  std::string which_property = "all"; // or one of: enpi, enpiLowt, enpiMidt, enpiHight
  BranchNames bn;
  int phi_nbins = 24;
  bool phi_in_degrees = false;
};

static void print_usage(const char* argv0) {
  std::cout << "Usage:\n  " << argv0 << " data.root gen.root recmc.root [enpi|enpiLowt|enpiMidt|enpiHight|all]"
            << " [--x x] [--phi phi] [--t t] [--mx2 Mx2] [--phibins 24] [--phi-deg]\n";
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

    if (eat("--x", cfg.bn.x)) continue;
    if (eat("--phi", cfg.bn.phi)) continue;
    if (eat("--t", cfg.bn.t)) continue;
    if (eat("--mx2", cfg.bn.mx2)) continue;
    if (eat_int("--phibins", cfg.phi_nbins)) continue;
    if (a == "--phi-deg") { cfg.phi_in_degrees = true; continue; }

    std::cerr << "Unknown argument: " << a << "\n";
  }
  return true;
}

static std::vector<std::string> properties_to_run(const std::string& which) {
  if (which == "all") return {"enpi", "enpiLowt", "enpiMidt", "enpiHight"};
  return {which};
}

static std::vector<double> x_edges() {
  // 6 x bins: [0.06-0.14), [0.14-0.22), [0.22-0.30), [0.30-0.40), [0.40-0.50), [0.50-0.60)
  return {0.06, 0.14, 0.22, 0.30, 0.40, 0.50, 0.60};
}

static int find_bin(double x, const std::vector<double>& edges) {
  for (size_t i=0;i+1<edges.size();++i) {
    if (in_range(x, edges[i], edges[i+1])) return static_cast<int>(i);
  }
  return -1;
}

static std::vector<std::unique_ptr<TH1D>>
make_phi_harray(const std::string& prefix, int nbins, double pmin, double pmax, int nbx) {
  std::vector<std::unique_ptr<TH1D>> v; v.reserve(nbx);
  for (int i=0;i<nbx;i++) {
    std::string name = prefix + Form("_bin%d", i);
    auto h = std::make_unique<TH1D>(name.c_str(), "", nbins, pmin, pmax);
    h->Sumw2();
    v.emplace_back(std::move(h));
  }
  return v;
}

static void fill_histos_per_file(
  TTree* tr, const BranchNames& bn, const std::string& prop,
  const std::vector<double>& xedges, std::vector<std::unique_ptr<TH1D>>& Hphi,
  int phi_nbins, double phi_min, double phi_max,
  bool apply_fiducial_cut, bool phi_is_deg)
{
  if (!tr) return;
  BranchHandles bh;
  if (!bind_branches(tr, bn, bh, apply_fiducial_cut)) {
    std::cerr << "Error: could not bind branches (check names). Tree: " << tr->GetName() << "\n";
    return;
  }
  const Long64_t nent = tr->GetEntries();
  for (Long64_t i=0;i<nent;i++) {
    tr->GetEntry(i);

    if (!Cuts::pass_common(bh.t, bh.mx2)) continue;
    if (!Cuts::pass_property_window(prop, bh.t)) continue;
    if (apply_fiducial_cut) {
      if (!bh.has_fid) continue;
      if (bh.fid < 100) continue;
    }

    int ib = find_bin(bh.x, xedges);
    if (ib < 0) continue;

    double phi = bh.phi;
    if (phi_is_deg || std::fabs(phi) > 3.5) {
      phi = phi * TMath::Pi() / 180.0;
    }
    phi = wrap_phi_0_2pi(phi);

    Hphi[ib]->Fill(phi);
  }
}

// Build unfolded TGraphErrors for one x bin and fit C*(1 + A*cos + B*cos2)
static FitResult build_graph_and_fit_from_histos(
  TH1D* hD, TH1D* hG, TH1D* hR, TGraphErrors*& out_graph, TF1*& out_fit)
{
  FitResult fr;
  if (!hD || !hG || !hR) return fr;

  const int nb = hD->GetNbinsX();
  std::vector<double> vx, vy, vex, vey; vx.reserve(nb); vy.reserve(nb); vex.assign(nb, 0.0); vey.reserve(nb);

  for (int b=1; b<=nb; ++b) {
    const double D = hD->GetBinContent(b);
    const double G = hG->GetBinContent(b);
    const double R = hR->GetBinContent(b);

    // Skip points where we cannot evaluate the ratio
    if (R <= 0.0 || G <= 0.0) continue;

    const double f = (D * G) / R;

    // Error propagation (Poisson on D,G,R):
    // f = D * G / R
    // Var(f) = (G/R)^2 Var(D) + (D/R)^2 Var(G) + (D*G/R^2)^2 Var(R)
    // with Var(D)=D, Var(G)=G, Var(R)=R
    const double term1 = (G*G)/(R*R) * D;
    const double term2 = (D*D)/(R*R) * G;
    const double term3 = (D*D*G*G)/(R*R*R); // (D*G/R^2)^2 * R = D^2 G^2 / R^3
    const double sigma = std::sqrt(std::max(0.0, term1 + term2 + term3));

    const double xc = hD->GetBinCenter(b);
    vx.push_back(xc);
    vy.push_back(f);
    vey.push_back(sigma);
  }

  if (vx.empty()) {
    out_graph = nullptr; out_fit = nullptr; return fr;
  }

  out_graph = new TGraphErrors((int)vx.size(), vx.data(), vy.data(), nullptr, vey.data());
  out_graph->SetMarkerStyle(20);
  out_graph->SetMarkerSize(0.8);
  out_graph->SetLineWidth(1);

  // Fit function: C * (1 + A*cos(x) + B*cos(2*x))
  const double xmin = hD->GetXaxis()->GetXmin();
  const double xmax = hD->GetXaxis()->GetXmax();
  out_fit = new TF1("fitCAB", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", xmin, xmax);
  out_fit->SetParNames("C","A","B");

  // Initial guesses
  double ymean = 0.0; for (double v : vy) ymean += v; ymean /= (double)vy.size();
  out_fit->SetParameters( (ymean>0?ymean:1.0), 0.0, 0.0 );

  // Reasonable bounds to stabilize (optional)
  out_fit->SetParLimits(0, 0.0, 1e12); // C >= 0
  out_fit->SetParLimits(1, -2.0, 2.0);
  out_fit->SetParLimits(2, -2.0, 2.0);

  // Weighted chi2 fit using errors from TGraphErrors
  out_graph->Fit(out_fit, "Q"); // quiet; remove Q if you want fit printouts

  fr.C  = out_fit->GetParameter(0);
  fr.dC = out_fit->GetParError(0);
  fr.A  = out_fit->GetParameter(1);
  fr.dA = out_fit->GetParError(1);
  fr.B  = out_fit->GetParameter(2);
  fr.dB = out_fit->GetParError(2);
  fr.chi2 = out_fit->GetChisquare();
  fr.ndf  = out_fit->GetNDF();

  return fr;
}

static void draw_one_property(
  const Config& cfg,
  const std::string& prop,
  TFile* fData, TFile* fGen, TFile* fRec)
{
  std::cout << "Processing property: " << prop << "\n";

  TTree* tData = (TTree*)fData->Get("PhysicsEvents");
  TTree* tGen  = (TTree*)fGen->Get("PhysicsEvents");
  TTree* tRec  = (TTree*)fRec->Get("PhysicsEvents");
  if (!tData || !tGen || !tRec) {
    std::cerr << "Error: could not load PhysicsEvents from one or more files.\n";
    return;
  }

  const auto xedges = x_edges();
  const int NX = (int)xedges.size()-1;

  // Phi axis
  const int    PHIB = cfg.phi_nbins;
  const double PHIM = 0.0;
  const double PHIX = 2.0 * TMath::Pi();

  // Build phi hist arrays for each x bin
  auto hD = make_phi_harray("hD", PHIB, PHIM, PHIX, NX);
  auto hG = make_phi_harray("hG", PHIB, PHIM, PHIX, NX);
  auto hR = make_phi_harray("hR", PHIB, PHIM, PHIX, NX);

  // Fill histograms
  fill_histos_per_file(tData, cfg.bn, prop, xedges, hD, PHIB, PHIM, PHIX, /*apply_fiducial_cut=*/true,  cfg.phi_in_degrees);
  fill_histos_per_file(tGen,  cfg.bn, prop, xedges, hG, PHIB, PHIM, PHIX, /*apply_fiducial_cut=*/false, cfg.phi_in_degrees);
  fill_histos_per_file(tRec,  cfg.bn, prop, xedges, hR, PHIB, PHIM, PHIX, /*apply_fiducial_cut=*/true,  cfg.phi_in_degrees);

  // Style tweaks
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.125);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  // Canvas 3x3
  auto c = new TCanvas(("c_"+prop).c_str(), ("Unfolded phi fits: "+prop).c_str(), 1200, 900);
  c->Divide(3,3);

  std::vector<double> xcenters; xcenters.reserve(NX);
  std::vector<double> Avec, dAvec, Bvec, dBvec; Avec.reserve(NX); dAvec.reserve(NX); Bvec.reserve(NX); dBvec.reserve(NX);

  for (int ib=0; ib<NX; ++ib) {
    c->cd(ib+1);

    TGraphErrors* g = nullptr;
    TF1* fit = nullptr;
    FitResult fr = build_graph_and_fit_from_histos(hD[ib].get(), hG[ib].get(), hR[ib].get(), g, fit);

    // Title for this pad
    const double xl = xedges[ib], xh = xedges[ib+1];
    std::string title = property_title(prop) + Form(", x in [%.2f, %.2f)", xl, xh);
    if (g) {
      g->SetTitle((title + ";#phi (rad);Unfolded counts").c_str());
      g->Draw("AP");
      if (fit) fit->Draw("LSAME");

      // legend with C, A, B and uncertainties, plus chi2/ndf
      auto leg = new TLegend(0.55, 0.60, 0.94, 0.94);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry((TObject*)0, Form("C = %.3g #pm %.3g", fr.C, fr.dC), "");
      leg->AddEntry((TObject*)0, Form("A = %.3g #pm %.3g", fr.A, fr.dA), "");
      leg->AddEntry((TObject*)0, Form("B = %.3g #pm %.3g", fr.B, fr.dB), "");
      leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.1f/%d", fr.chi2, fr.ndf), "");
      leg->Draw();

      // store A(x), B(x)
      xcenters.push_back(0.5*(xl+xh));
      Avec.push_back(fr.A); dAvec.push_back(fr.dA);
      Bvec.push_back(fr.B); dBvec.push_back(fr.dB);
    } else {
      // Empty pad labeling
      TH1D* frame = new TH1D("frame","",10,0,2*TMath::Pi());
      frame->SetTitle((title + ";#phi (rad);Unfolded counts").c_str());
      frame->SetMinimum(0); frame->SetMaximum(1);
      frame->Draw("AXIS");
      auto leg = new TLegend(0.55, 0.60, 0.94, 0.90);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      leg->AddEntry((TObject*)0, "No valid points", "");
      leg->Draw();
    }
  }

  // Bottom row, middle pad (pad 8): A(x) and B(x)
  c->cd(7); gPad->Clear(); // left bottom: empty on purpose
  c->cd(8);

  TGraphErrors* gA = nullptr;
  TGraphErrors* gB = nullptr;

  if (!xcenters.empty()) {
    gA = new TGraphErrors((int)xcenters.size());
    gB = new TGraphErrors((int)xcenters.size());
    for (int i=0;i<(int)xcenters.size();++i) {
      gA->SetPoint(i, xcenters[i], Avec[i]);
      gA->SetPointError(i, 0.0, dAvec[i]);
      gB->SetPoint(i, xcenters[i], Bvec[i]);
      gB->SetPointError(i, 0.0, dBvec[i]);
    }
    gA->SetMarkerStyle(20);
    gB->SetMarkerStyle(21);
    gA->SetMarkerSize(1.0);
    gB->SetMarkerSize(1.0);

    // Draw A first to set axes
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

  c->cd(9); gPad->Clear(); // right bottom: empty on purpose

  // Save
  ensure_dir("output/enpi+");
  const std::string outpdf = "output/enpi+/" + prop + "_unfolded.pdf";
  c->SaveAs(outpdf.c_str());
  std::cout << "Saved: " << outpdf << "\n";
}

int main(int argc, char** argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) return 1;

  std::unique_ptr<TFile> fData(TFile::Open(cfg.data_path.c_str(), "READ"));
  std::unique_ptr<TFile> fGen (TFile::Open(cfg.gen_path.c_str(),  "READ"));
  std::unique_ptr<TFile> fRec (TFile::Open(cfg.rec_path.c_str(),  "READ"));
  if (!fData || fData->IsZombie() || !fGen || fGen->IsZombie() || !fRec || fRec->IsZombie()) {
    std::cerr << "Error opening input files.\n";
    return 2;
  }

  auto props = properties_to_run(cfg.which_property);
  for (const auto& p : props) {
    draw_one_property(cfg, p, fData.get(), fGen.get(), fRec.get());
  }

  return 0;
}