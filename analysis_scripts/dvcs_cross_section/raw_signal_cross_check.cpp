#include "raw_signal_cross_check.h"
#include "load_csv.h"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include <TGaxis.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

// ---------------- small utilities ----------------
static inline void info(const std::string& s) { std::cout << "[cross] " << s << std::endl; }
static inline void warn(const std::string& s) { std::cout << "[cross][warn] " << s << std::endl; }

static inline std::string slower(std::string s) {
    for (auto& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

static constexpr int N_PHI_BINS = 12;
static inline int phiBinFromDeg(double phi_deg) {
    // robust binning for 12 bins on [0,360)
    if (!std::isfinite(phi_deg)) return -1;
    double p = fmod(phi_deg, 360.0);
    if (p < 0.0) p += 360.0;
    const double w = 360.0 / double(N_PHI_BINS);
    int ip = (int)std::floor(p / w);
    if (ip < 0) ip = 0;
    if (ip >= N_PHI_BINS) ip = N_PHI_BINS - 1;
    return ip;
}

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> v(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) v[i] = (i + 0.5)*step;
    return v;
}

static inline void degreeTicks(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

// ---------------- load total_counts master and pick Fa18 groups ----------------
static bool load_master_groups(const std::string& master_path,
                               json& group_inb, json& group_out,
                               std::string& name_inb, std::string& name_out)
{
    std::ifstream f(master_path);
    if (!f.is_open()) {
        warn("Cannot open total-counts master: " + master_path);
        return false;
    }
    json J; f >> J;
    if (!J.contains("groups")) {
        warn("Master JSON has no 'groups' key.");
        return false;
    }
    const json& G = J["groups"];
    std::string found_inb, found_out;
    json jinb, jout;

    for (auto it = G.begin(); it != G.end(); ++it) {
        std::string k = it.key();
        std::string kl = slower(k);
        if (kl.find("fa18_inb") != std::string::npos || kl.find("fa18inb") != std::string::npos) {
            jinb = it.value();
            found_inb = k;
        }
        if (kl.find("fa18_out") != std::string::npos || kl.find("fa18out") != std::string::npos) {
            jout = it.value();
            found_out = k;
        }
    }
    if (!jinb.contains("bins") || !jout.contains("bins")) {
        warn("Could not find both Fa18 inb/out groups inside master JSON.");
        return false;
    }
    group_inb = jinb["bins"];
    group_out = jout["bins"];
    name_inb  = found_inb;
    name_out  = found_out;
    return true;
}

// ---------------- build indices from CSV (Lee) ----------------
struct AxisSets {
    std::vector<std::pair<double,double>> xB;
    std::map<int, std::vector<std::pair<double,double>>> Q2_by_ix;
    std::map<int, std::vector<std::pair<double,double>>> t_by_ix;
};

static AxisSets build_axes_from_rows(const std::vector<LeeRow>& rows) {
    std::set<std::pair<double,double>> xbset;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> q2set_by_xb;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> tset_by_xb;

    for (const auto& r : rows) {
        std::pair<double,double> xb(r.xBmin, r.xBmax);
        xbset.insert(xb);
        q2set_by_xb[xb].insert(std::make_pair(r.Q2min, r.Q2max));
        tset_by_xb[xb].insert(std::make_pair(r.tmin,  r.tmax));
    }
    AxisSets ax;
    ax.xB.assign(xbset.begin(), xbset.end());
    for (int ix = 0; ix < (int)ax.xB.size(); ++ix) {
        const auto& xb = ax.xB[ix];
        auto& q2set = q2set_by_xb[xb];
        auto& tset  = tset_by_xb[xb];
        ax.Q2_by_ix[ix] = std::vector<std::pair<double,double>>(q2set.begin(), q2set.end());
        ax.t_by_ix[ix]  = std::vector<std::pair<double,double>>(tset.begin(),  tset.end());
    }
    return ax;
}

static inline int find_index(const std::pair<double,double>& r,
                             const std::vector<std::pair<double,double>>& v)
{
    for (int i=0;i<(int)v.size();++i) if (v[i]==r) return i;
    return -1;
}

// ---------------- collect per-bin Lee counts into index space ----------------
struct LeePerBin {
    // [ix][iQ][it][ip] -> counts (separate for inb/out)
    std::map<std::tuple<int,int,int,int>, double> inb;
    std::map<std::tuple<int,int,int,int>, double> out;
};

static LeePerBin map_lee_to_indices(const std::vector<LeeRow>& rows,
                                    const AxisSets& ax)
{
    LeePerBin m;
    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        int ix = find_index(xb, ax.xB); if (ix<0) continue;
        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);

        int iQ = find_index(std::make_pair(r.Q2min,r.Q2max), Q2s);
        int it = find_index(std::make_pair(r.tmin, r.tmax),  Ts);
        if (iQ<0 || it<0) continue;

        int ip = phiBinFromDeg(r.phiavg);
        if (ip<0) continue;

        m.inb[{ix,iQ,it,ip}] += r.raw_inb_sum; // sum in case of duplicates
        m.out[{ix,iQ,it,ip}] += r.raw_out_sum;
    }
    return m;
}

// ---------------- extract our counts from JSON group ----------------
static std::map<std::tuple<int,int,int,int>, double>
extract_my_counts(const json& bins_object)
{
    std::map<std::tuple<int,int,int,int>, double> out;
    for (auto it = bins_object.begin(); it != bins_object.end(); ++it) {
        const std::string& key = it.key(); // "(ix,iQ,it,ip)"
        int ix=0,iQ=0,itb=0,ip=0;
        if (std::sscanf(key.c_str(), "(%d,%d,%d,%d)", &ix,&iQ,&itb,&ip) != 4) continue;
        const json& cell = it.value();
        double total = 0.0;
        if (cell.contains("total")) {
            try { total = cell["total"].get<double>(); } catch(...) { total = 0.0; }
        } else if (cell.contains("helicity")) {
            // fallback
            const json& h = cell["helicity"];
            double p = 0.0, m = 0.0;
            try { p = h["+1"].get<double>(); } catch(...) {}
            try { m = h["-1"].get<double>(); } catch(...) {}
            total = p + m;
        }
        out[{ix,iQ,itb,ip}] = total;
    }
    return out;
}

// ---------------- plotting core ----------------
static void draw_degree_ticks_here(double labelSize) {
    degreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), labelSize);
}

static void draw_one_canvas(const std::string& title,
                            const std::vector<std::pair<double,double>>& Q2s,
                            const std::vector<std::pair<double,double>>& Ts,
                            const std::vector<double>& phiC,
                            // series A and B: vectors over ip=0..11
                            const std::function<bool(int,int,int,std::vector<double>&)>& fillA,
                            const std::function<bool(int,int,int,std::vector<double>&)>& fillB,
                            const std::string& out_png,
                            bool draw_ratio_only)
{
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows==0 || ncols==0) return;

    const int W = 320*ncols + 220;  // match style from your recent plots
    const int H = 260*nrows + 260;

    TCanvas* c = new TCanvas(("c_"+out_png).c_str(), "", W, H);
    c->cd();

    // top band
    TPad* pTop = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
    pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
    pTop->cd();

    TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.40);
    head.DrawLatex(0.50, 0.68, title.c_str());

    // legend (single, compact)
    // series A: ours, series B: Lee
    TGraph a; a.SetMarkerStyle(20); a.SetMarkerColor(kBlack); a.SetLineColor(kBlack); a.SetLineWidth(3);
    TGraph b; b.SetMarkerStyle(24); b.SetMarkerColor(kOrange+7); b.SetLineColor(kOrange+7); b.SetLineWidth(3);
    TLegend* legTop = new TLegend(0.08, 0.10, 0.92, 0.56);
    legTop->SetNColumns(2);
    legTop->SetBorderSize(0);
    legTop->SetFillStyle(0);
    legTop->SetTextFont(42);
    legTop->SetTextSize(0.22);
    legTop->AddEntry(&a, draw_ratio_only ? "Hayward/Lee" : "Hayward (this work)", "lp");
    legTop->AddEntry(&b, draw_ratio_only ? "y = 1"       : "Lee (legacy)",       "lp");
    legTop->Draw();

    // grid
    c->cd();
    TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.90);
    pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.00, 0.00);

    // for each panel
    std::vector<double> A(N_PHI_BINS, 0.0), Bv(N_PHI_BINS, 0.0), ex(N_PHI_BINS, 0.0);
    for (int r=0; r<nrows; ++r) {
        for (int ccol=0; ccol<ncols; ++ccol) {
            pGrid->cd(r*ncols + ccol + 1);
            gPad->SetGrid(1,1);
            gPad->SetTicks(1,1);
            gPad->SetTopMargin(0.06);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.18);   // extra padding
            gPad->SetRightMargin(0.08);

            bool haveA = fillA(ccol, r, /*ignored*/0, A);
            bool haveB = fillB(ccol, r, /*ignored*/0, Bv);

            double ymin=0.0, ymax=1.0;
            if (!draw_ratio_only) {
                // counts mode: set a nice headroom based on both
                double local_max = 0.0;
                for (int ip=0; ip<N_PHI_BINS; ++ip) local_max = std::max(local_max, std::max(A[ip], Bv[ip]));
                ymax = std::max(1.0, local_max*1.25);
            } else {
                // ratio mode: show around 1
                ymin = 0.5; ymax = 1.5;
            }

            TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle(draw_ratio_only ? "Hayward / Lee" : "Raw counts");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.050);
            frame->GetXaxis()->SetLabelSize(0.050);
            frame->GetXaxis()->SetTitleOffset(1.10);
            frame->GetYaxis()->SetTitleOffset(1.55);

            draw_degree_ticks_here(0.050);

            // per-panel label
            TLatex lab; lab.SetNDC(); lab.SetTextSize(0.070); lab.SetTextAlign(11); lab.SetTextFont(42);
            lab.DrawLatex(0.15, 0.94,
                Form("Q^{2} #in [%.2g, %.2g],  -t #in [%.2g, %.2g]",
                     Q2s[ccol].first, Q2s[ccol].second, Ts[r].first, Ts[r].second));

            if (draw_ratio_only) {
                // ratio points
                std::vector<double> R(N_PHI_BINS, 0.0);
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const double denom = Bv[ip];
                    R[ip] = (denom > 0.0) ? (A[ip]/denom) : NAN;
                }
                TGraph* gr = new TGraph(N_PHI_BINS,
                                        const_cast<double*>(phiC.data()),
                                        const_cast<double*>(R.data()));
                gr->SetMarkerStyle(20);
                gr->SetMarkerColor(kBlack);
                gr->SetLineColor(kBlack);
                gr->SetLineWidth(2);
                gr->Draw("LP SAME");

                TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                one->SetLineStyle(2);
                one->SetLineWidth(2);
                one->SetLineColor(kOrange+7);
                one->Draw("SAME");
            } else {
                // A (ours)
                TGraph* gA = new TGraph(N_PHI_BINS,
                                        const_cast<double*>(phiC.data()),
                                        const_cast<double*>(A.data()));
                gA->SetMarkerStyle(20);
                gA->SetMarkerColor(kBlack);
                gA->SetLineColor(kBlack);
                gA->SetLineWidth(3);
                gA->Draw("LP SAME");

                // B (Lee)
                TGraph* gB = new TGraph(N_PHI_BINS,
                                        const_cast<double*>(phiC.data()),
                                        const_cast<double*>(Bv.data()));
                gB->SetMarkerStyle(24);
                gB->SetMarkerColor(kOrange+7);
                gB->SetLineColor(kOrange+7);
                gB->SetLineWidth(3);
                gB->Draw("LP SAME");
            }
        }
    }

    c->SaveAs(out_png.c_str());
    delete c;
}

// ---------------- driver ----------------
void plot_raw_yield_cross_checks(const std::vector<LeeRow>& rows,
                                 const std::string& total_counts_master_json,
                                 const std::string& output_base_dir)
{
    fs::create_directories(output_base_dir);

    // style (mirror your recent choices)
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);

    // 1) Load master and find Fa18 groups
    json bins_inb, bins_out; std::string name_inb, name_out;
    if (!load_master_groups(total_counts_master_json, bins_inb, bins_out, name_inb, name_out)) {
        warn("Cannot proceed without Fa18 inb/out groups.");
        return;
    }
    info("Using groups: " + name_inb + " and " + name_out);

    // 2) Build axes from CSV rows (defines which panels to draw)
    AxisSets ax = build_axes_from_rows(rows);
    const auto phiC = phiCentersDeg();

    // 3) Map Lee rows into index space
    LeePerBin lee = map_lee_to_indices(rows, ax);

    // 4) Extract our totals per (ix,iQ,it,ip) from JSON groups
    auto ours_inb = extract_my_counts(bins_inb);
    auto ours_out = extract_my_counts(bins_out);

    // Helper to fill series for a given ix, iQ, it: return false if empty
    auto fill_series = [&](const std::map<std::tuple<int,int,int,int>, double>& ours,
                           const std::map<std::tuple<int,int,int,int>, double>& lee_side,
                           int ix, int iQ, int it,
                           std::vector<double>& A, std::vector<double>& B)->bool
    {
        bool any=false;
        for (int ip=0; ip<N_PHI_BINS; ++ip) {
            auto k = std::make_tuple(ix,iQ,it,ip);
            double a = 0.0, b = 0.0;
            auto itA = ours.find(k);
            if (itA != ours.end()) a = itA->second;
            auto itB = lee_side.find(k);
            if (itB != lee_side.end()) b = itB->second;
            A[ip] = a;
            B[ip] = b;
            any |= (a>0.0 || b>0.0);
        }
        return any;
    };

    // For each xB slice, draw 2 canvases for inb and 2 for out
    for (int ix=0; ix<(int)ax.xB.size(); ++ix) {
        const auto& Q2s = ax.Q2_by_ix[ix];
        const auto& Ts  = ax.t_by_ix[ix];
        if (Q2s.empty() || Ts.empty()) continue;

        // INBENDING
        {
            const std::string title_counts =
                Form("Raw yields: %s   x_{B} #in [%.3g, %.3g]", name_inb.c_str(), ax.xB[ix].first, ax.xB[ix].second);
            const std::string title_ratio  =
                Form("Raw yields ratio (Hayward/Lee): %s   x_{B} #in [%.3g, %.3g]", name_inb.c_str(), ax.xB[ix].first, ax.xB[ix].second);

            auto fillA_inb = [&](int iQcol, int irow, int, std::vector<double>& A)->bool {
                // A = ours
                std::vector<double> tmp(N_PHI_BINS,0.0);
                return fill_series(ours_inb, lee.inb, ix, iQcol, irow, A, tmp);
            };
            auto fillB_inb = [&](int iQcol, int irow, int, std::vector<double>& B)->bool {
                // B = lee
                std::vector<double> tmp(N_PHI_BINS,0.0);
                return fill_series(ours_inb, lee.inb, ix, iQcol, irow, tmp, B);
            };

            const std::string f_counts = (fs::path(output_base_dir)/Form("raw_counts_fa18_inb_xB_%d.png", ix)).string();
            const std::string f_ratio  = (fs::path(output_base_dir)/Form("raw_ratio_fa18_inb_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, phiC, fillA_inb, fillB_inb, f_counts, /*ratio=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, phiC, fillA_inb, fillB_inb, f_ratio,  /*ratio=*/true);
            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }

        // OUTBENDING
        {
            const std::string title_counts =
                Form("Raw yields: %s   x_{B} #in [%.3g, %.3g]", name_out.c_str(), ax.xB[ix].first, ax.xB[ix].second);
            const std::string title_ratio  =
                Form("Raw yields ratio (Hayward/Lee): %s   x_{B} #in [%.3g, %.3g]", name_out.c_str(), ax.xB[ix].first, ax.xB[ix].second);

            auto fillA_out = [&](int iQcol, int irow, int, std::vector<double>& A)->bool {
                std::vector<double> tmp(N_PHI_BINS,0.0);
                return fill_series(ours_out, lee.out, ix, iQcol, irow, A, tmp);
            };
            auto fillB_out = [&](int iQcol, int irow, int, std::vector<double>& B)->bool {
                std::vector<double> tmp(N_PHI_BINS,0.0);
                return fill_series(ours_out, lee.out, ix, iQcol, irow, tmp, B);
            };

            const std::string f_counts = (fs::path(output_base_dir)/Form("raw_counts_fa18_out_xB_%d.png", ix)).string();
            const std::string f_ratio  = (fs::path(output_base_dir)/Form("raw_ratio_fa18_out_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, phiC, fillA_out, fillB_out, f_counts, /*ratio=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, phiC, fillA_out, fillB_out, f_ratio,  /*ratio=*/true);
            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }
    }
}