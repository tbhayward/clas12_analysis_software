// raw_signal_cross_check.cpp
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
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cctype>
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

// ---------- small utilities ----------
static inline void info(const std::string& s) { std::cout << "[cross] " << s << std::endl; }
static inline void warn(const std::string& s) { std::cout << "[cross][warn] " << s << std::endl; }
static inline std::string slower(std::string s) { for (auto& c : s) c = (char)std::tolower((unsigned char)c); return s; }

static constexpr int N_PHI_BINS = 12;

static inline int phiBinFromDeg(double phi_deg) {
    if (!std::isfinite(phi_deg)) return -1;
    double p = std::fmod(phi_deg, 360.0);
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

// Prefer exact “fa18_inb/out”; otherwise accept keys with “fa18” and (“inb”/“out”) but NOT “supp”.
static std::string choose_fa18_key(const std::vector<std::string>& keys, bool want_inb) {
    const std::string exact = want_inb ? "fa18_inb" : "fa18_out";
    for (const auto& k : keys) if (slower(k) == exact) return k;
    std::vector<std::string> cands;
    for (const auto& k : keys) {
        const std::string kl = slower(k);
        const bool ok = kl.find("fa18")!=std::string::npos &&
                        kl.find("supp")==std::string::npos &&
                        ((want_inb && kl.find("inb")!=std::string::npos) ||
                         (!want_inb && kl.find("out")!=std::string::npos));
        if (ok) cands.push_back(k);
    }
    if (!cands.empty()) return cands.front();
    return std::string();
}

// ---------- load total_counts master and pick Fa18 groups ----------
static bool load_master_groups(const std::string& master_path,
                               json& group_inb, json& group_out,
                               std::string& name_inb, std::string& name_out) {
    std::ifstream f(master_path);
    if (!f.is_open()) { warn("Cannot open total-counts master: " + master_path); return false; }
    json J; f >> J;
    if (!J.contains("groups")) { warn("Master JSON has no 'groups' key."); return false; }
    const json& G = J["groups"];

    std::vector<std::string> names; names.reserve(G.size());
    for (auto it = G.begin(); it != G.end(); ++it) names.push_back(it.key());

    name_inb = choose_fa18_key(names, /*want_inb=*/true);
    name_out = choose_fa18_key(names, /*want_inb=*/false);
    if (name_inb.empty() || name_out.empty()) {
        warn("Could not find both Fa18_inb and Fa18_out (excluding any *_supp).");
        return false;
    }
    const json& jinb = G.at(name_inb);
    const json& jout = G.at(name_out);
    if (!jinb.contains("bins") || !jout.contains("bins")) { warn("Chosen groups lack 'bins'."); return false; }
    group_inb = jinb["bins"];
    group_out = jout["bins"];
    return true;
}

// ---------- axes from CSV (Lee) ----------
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
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        xbset.insert(xb);
        q2set_by_xb[xb].insert({r.Q2min, r.Q2max});
        tset_by_xb[xb].insert({r.tmin,  r.tmax});
    }

    AxisSets ax;
    ax.xB.assign(xbset.begin(), xbset.end());
    for (int ix=0; ix<(int)ax.xB.size(); ++ix) {
        const auto& xb = ax.xB[ix];
        ax.Q2_by_ix[ix] = { q2set_by_xb[xb].begin(), q2set_by_xb[xb].end() };
        ax.t_by_ix[ix]  = {  tset_by_xb[xb].begin(),  tset_by_xb[xb].end()  };
    }
    return ax;
}

static inline int find_index(const std::pair<double,double>& r,
                             const std::vector<std::pair<double,double>>& v) {
    for (int i=0;i<(int)v.size();++i) if (v[i]==r) return i; return -1;
}

// ---------- Lee counts in index space ----------
struct LeePerBin {
    std::map<std::tuple<int,int,int,int>, double> inb;
    std::map<std::tuple<int,int,int,int>, double> out;
};

static LeePerBin map_lee_to_indices(const std::vector<LeeRow>& rows,
                                    const AxisSets& ax) {
    LeePerBin m;
    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        const int ix  = find_index(xb, ax.xB);
        if (ix<0) continue;
        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);

        const int iQ = find_index({r.Q2min,r.Q2max}, Q2s);
        const int it = find_index({r.tmin, r.tmax},  Ts);
        if (iQ<0 || it<0) continue;

        const int ip = phiBinFromDeg(r.phiavg);
        if (ip<0) continue;

        m.inb[{ix,iQ,it,ip}] += r.raw_inb_sum;
        m.out[{ix,iQ,it,ip}] += r.raw_out_sum;
    }
    return m;
}

// ---------- extract our counts ----------
static std::map<std::tuple<int,int,int,int>, double>
extract_my_counts(const json& bins_object) {
    std::map<std::tuple<int,int,int,int>, double> out;
    for (auto it = bins_object.begin(); it != bins_object.end(); ++it) {
        int ix=0,iQ=0,itb=0,ip=0;
        if (std::sscanf(it.key().c_str(), "(%d,%d,%d,%d)", &ix,&iQ,&itb,&ip) != 4) continue;
        const json& cell = it.value();
        double total = 0.0;
        if (cell.contains("total")) {
            try { total = cell["total"].get<double>(); } catch (...) { total = 0.0; }
        } else if (cell.contains("helicity")) {
            const json& h = cell["helicity"];
            double p=0.0, m=0.0;
            try { p = h["+1"].get<double>(); } catch (...) {}
            try { m = h["-1"].get<double>(); } catch (...) {}
            total = p + m;
        }
        out[{ix,iQ,itb,ip}] = total;
    }
    return out;
}

// ---------- plotting helpers ----------
static void draw_degree_ticks_here(double labelSize) {
    degreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), labelSize);
}

// Markers + error bars only (no connecting line)
static TGraphErrors* graph_pe1(const std::vector<double>& X,
                               const std::vector<double>& Y,
                               const std::vector<double>& EY,
                               int markerStyle, int color) {
    std::vector<double> ex(X.size(), 0.0);
    auto* g = new TGraphErrors((int)X.size(),
                               const_cast<double*>(X.data()),
                               const_cast<double*>(Y.data()),
                               ex.data(),
                               const_cast<double*>(EY.data()));
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);
    g->Draw("PE1 SAME"); // points + error bars, no connecting line
    return g;
}

static void draw_one_canvas(const std::string& title,
                            const std::vector<std::pair<double,double>>& Q2s,
                            const std::vector<std::pair<double,double>>& Ts,
                            const std::vector<double>& phiC,
                            const std::function<bool(int,int,std::vector<double>&,std::vector<double>&)>& fillBoth,
                            const std::string& out_png,
                            bool draw_ratio_only) {
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows==0 || ncols==0) return;

    const int W = 320*ncols + 220;
    const int H = 260*nrows + 260;

    TCanvas* c = new TCanvas(("c_"+out_png).c_str(), "", W, H);
    c->cd();

    // Top band (smaller title)
    TPad* pTop = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
    pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
    pTop->cd();

    TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42);
    head.SetTextSize(0.18); // ~half the previous size
    head.DrawLatex(0.50, 0.65, title.c_str());

    // Legend
    c->cd();
    TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.90);
    pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
    pGrid->cd();

    // Place legend in the top band for consistency
    pTop->cd();
    TLegend* legTop = new TLegend(0.08, 0.10, 0.92, 0.56);
    legTop->SetNColumns(2);
    legTop->SetBorderSize(0);
    legTop->SetFillStyle(0);
    legTop->SetTextFont(42);
    legTop->SetTextSize(0.22);
    if (draw_ratio_only) {
        TGraph dummyRatio; dummyRatio.SetMarkerStyle(20); dummyRatio.SetMarkerColor(kBlack);
        TGraph y1;         y1.SetLineStyle(2); y1.SetLineColor(kOrange+7);
        legTop->AddEntry(&dummyRatio, "Hayward/Lee", "p");
        legTop->AddEntry(&y1,         "y = 1",       "l");
    } else {
        TGraph h; h.SetMarkerStyle(20); h.SetMarkerColor(kBlack);
        TGraph l; l.SetMarkerStyle(24); l.SetMarkerColor(kOrange+7);
        legTop->AddEntry(&h, "Hayward (pass-2)", "p");
        legTop->AddEntry(&l, "Lee (pass-1)",     "p");
    }
    legTop->Draw();

    // Grid
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.00, 0.00);

    // Panels
    for (int r=0; r<nrows; ++r) {
        for (int ccol=0; ccol<ncols; ++ccol) {
            pGrid->cd(r*ncols + ccol + 1);
            gPad->SetGrid(1,1);
            gPad->SetTicks(1,1);
            gPad->SetTopMargin(0.10);   // a touch more space at top
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.08);

            std::vector<double> A(N_PHI_BINS,0.0), B(N_PHI_BINS,0.0);
            (void)fillBoth(ccol, r, A, B);

            double ymin=0.0, ymax=1.0;
            TH1* frame = nullptr;

            if (!draw_ratio_only) {
                double local_max = 0.0;
                for (int ip=0; ip<N_PHI_BINS; ++ip) local_max = std::max(local_max, std::max(A[ip], B[ip]));
                ymax = std::max(1.0, local_max*1.25);
                frame = gPad->DrawFrame(0.0, 0.0, 360.0, ymax);
            } else {
                ymin = 0.5; ymax = 1.5;
                frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
            }

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

            // Panel label (nudged UP a bit; was ~0.92)
            TLatex lab; lab.SetNDC(); lab.SetTextSize(0.070); lab.SetTextAlign(11); lab.SetTextFont(42);
            lab.DrawLatex(0.15, 0.95,
                Form("Q^{2} #in [%.2g, %.2g],  -t #in [%.2g, %.2g]",
                     Q2s[ccol].first, Q2s[ccol].second, Ts[r].first, Ts[r].second));

            const int black  = kBlack;
            const int orange = kOrange+7;

            if (draw_ratio_only) {
                std::vector<double> x, y, ey;
                x.reserve(N_PHI_BINS); y.reserve(N_PHI_BINS); ey.reserve(N_PHI_BINS);
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const double NL = B[ip];
                    const double NH = A[ip];
                    if (NL <= 0.0) continue; // undefined ratio
                    const double R  = (NH <= 0.0) ? 0.0 : NH/NL;
                    double eR = 0.0;
                    if (NH > 0.0) eR = R * std::sqrt(1.0/NH + 1.0/NL); // Poisson prop
                    x.push_back(phiC[ip]);
                    y.push_back(R);
                    ey.push_back(eR);
                }
                graph_pe1(x, y, ey, 20, black);

                TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                one->SetLineStyle(2);
                one->SetLineWidth(2);
                one->SetLineColor(orange);
                one->Draw("SAME");
            } else {
                std::vector<double> eyA(N_PHI_BINS, 0.0), eyB(N_PHI_BINS, 0.0);
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    eyA[ip] = (A[ip] > 0.0) ? std::sqrt(A[ip]) : 0.0;
                    eyB[ip] = (B[ip] > 0.0) ? std::sqrt(B[ip]) : 0.0;
                }
                graph_pe1(phiC, A, eyA, 20, black);   // Hayward (pass-2)
                graph_pe1(phiC, B, eyB, 24, orange);  // Lee (pass-1)
            }
        }
    }

    c->SaveAs(out_png.c_str());
    delete c;
}

// ---------- driver ----------
void plot_raw_yield_cross_checks(const std::vector<LeeRow>& rows,
                                 const std::string& total_counts_master_json,
                                 const std::string& output_base_dir) {
    fs::create_directories(output_base_dir);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);

    json bins_inb, bins_out; std::string name_inb, name_out;
    if (!load_master_groups(total_counts_master_json, bins_inb, bins_out, name_inb, name_out)) {
        warn("Cannot proceed without Fa18 inb/out groups.");
        return;
    }
    info("Using groups: " + name_inb + " and " + name_out);

    AxisSets ax = build_axes_from_rows(rows);
    const auto phiC = phiCentersDeg();

    LeePerBin lee = map_lee_to_indices(rows, ax);
    auto ours_inb = extract_my_counts(bins_inb);
    auto ours_out = extract_my_counts(bins_out);

    auto make_fillBoth = [&](const std::map<std::tuple<int,int,int,int>, double>& ours,
                             const std::map<std::tuple<int,int,int,int>, double>& lee_side,
                             int ix) {
        return [&,ix](int iQcol, int irow, std::vector<double>& A, std::vector<double>& B)->bool {
            bool any=false;
            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                auto k = std::make_tuple(ix,iQcol,irow,ip);
                double a = 0.0, b = 0.0;
                auto itA = ours.find(k); if (itA != ours.end()) a = itA->second;
                auto itB = lee_side.find(k); if (itB != lee_side.end()) b = itB->second;
                A[ip] = a; B[ip] = b; any |= (a>0.0 || b>0.0);
            }
            return any;
        };
    };

    for (int ix=0; ix<(int)ax.xB.size(); ++ix) {
        const auto& Q2s = ax.Q2_by_ix[ix];
        const auto& Ts  = ax.t_by_ix[ix];
        if (Q2s.empty() || Ts.empty()) continue;

        // INB
        {
            const std::string title_counts =
                Form("Raw yields: %s   x_{B} #in [%.3g, %.3g]", name_inb.c_str(), ax.xB[ix].first, ax.xB[ix].second);
            const std::string title_ratio  =
                Form("Raw yields ratio (Hayward/Lee): %s   x_{B} #in [%.3g, %.3g]", name_inb.c_str(), ax.xB[ix].first, ax.xB[ix].second);

            auto fillBoth_inb = make_fillBoth(ours_inb, lee.inb, ix);

            const std::string f_counts = (fs::path(output_base_dir)/Form("raw_counts_fa18_inb_xB_%d.png", ix)).string();
            const std::string f_ratio  = (fs::path(output_base_dir)/Form("raw_ratio_fa18_inb_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, phiC, fillBoth_inb, f_counts, /*ratio=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, phiC, fillBoth_inb, f_ratio,  /*ratio=*/true);
            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }

        // OUT
        {
            const std::string title_counts =
                Form("Raw yields: %s   x_{B} #in [%.3g, %.3g]", name_out.c_str(), ax.xB[ix].first, ax.xB[ix].second);
            const std::string title_ratio  =
                Form("Raw yields ratio (Hayward/Lee): %s   x_{B} #in [%.3g, %.3g]", name_out.c_str(), ax.xB[ix].first, ax.xB[ix].second);

            auto fillBoth_out = make_fillBoth(ours_out, lee.out, ix);

            const std::string f_counts = (fs::path(output_base_dir)/Form("raw_counts_fa18_out_xB_%d.png", ix)).string();
            const std::string f_ratio  = (fs::path(output_base_dir)/Form("raw_ratio_fa18_out_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, phiC, fillBoth_out, f_counts, /*ratio=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, phiC, fillBoth_out, f_ratio,  /*ratio=*/true);
            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }
    }
}