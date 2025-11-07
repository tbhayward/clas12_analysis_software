// pi0_contamination_cross_check.cpp
#include "pi0_contamination_cross_check.h"
#include "load_csv.h"  // LeeRow + loader/field names

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <TCanvas.h>
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
#include <tuple>
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

static void draw_degree_ticks_here(double labelSize) {
    degreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), labelSize);
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

// ---------- Lee contamination mapped to index space ----------
struct LeeContamPerBin {
    // store arithmetic mean if duplicates occur (rare)
    std::map<std::tuple<int,int,int,int>, double> inb;
    std::map<std::tuple<int,int,int,int>, double> out;
};

static LeeContamPerBin map_lee_contam_to_indices(const std::vector<LeeRow>& rows,
                                                 const AxisSets& ax) {
    LeeContamPerBin m;

    // Track counts to average if duplicates
    std::map<std::tuple<int,int,int,int>, std::pair<double,int>> tmp_inb;
    std::map<std::tuple<int,int,int,int>, std::pair<double,int>> tmp_out;

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

        auto key = std::make_tuple(ix,iQ,it,ip);

        // colleague's CSV provides single helicity-independent contamination per side
        if (std::isfinite(r.contam_inb) && r.contam_inb > 0.0) {
            auto& a = tmp_inb[key];
            a.first += r.contam_inb;
            a.second += 1;
        }
        if (std::isfinite(r.contam_out) && r.contam_out > 0.0) {
            auto& b = tmp_out[key];
            b.first += r.contam_out;
            b.second += 1;
        }
    }

    for (const auto& kv : tmp_inb) m.inb[kv.first] = kv.second.first / std::max(1, kv.second.second);
    for (const auto& kv : tmp_out) m.out[kv.first] = kv.second.first / std::max(1, kv.second.second);
    return m;
}

// ---------- load our pi0 contamination (helicity-averaged) from combined JSON ----------
struct OurContam {
    // value and error (helicity-averaged with N_data weights)
    std::map<std::tuple<int,int,int,int>, std::pair<double,double>> inb;
    std::map<std::tuple<int,int,int,int>, std::pair<double,double>> out;
};

static bool load_our_contam_groups(const std::string& combined_path,
                                   json& bins_inb, json& bins_out,
                                   std::string& name_inb, std::string& name_out) {
    std::ifstream f(combined_path);
    if (!f.is_open()) { warn("Cannot open pi0 combined JSON: " + combined_path); return false; }
    json J; f >> J;

    if (!J.contains("periods")) { warn("Combined pi0 JSON has no 'periods' key."); return false; }
    const json& P = J["periods"];

    std::vector<std::string> names; names.reserve(P.size());
    for (auto it = P.begin(); it != P.end(); ++it) names.push_back(it.key());

    name_inb = choose_fa18_key(names, /*want_inb=*/true);
    name_out = choose_fa18_key(names, /*want_inb=*/false);
    if (name_inb.empty() || name_out.empty()) {
        warn("Could not find both fa18_inb and fa18_out in combined pi0 JSON (excluding any *_supp).");
        return false;
    }
    const json& jinb = P.at(name_inb);
    const json& jout = P.at(name_out);
    if (!jinb.contains("bins") || !jout.contains("bins")) { warn("Chosen periods lack 'bins'."); return false; }
    bins_inb = jinb["bins"];
    bins_out = jout["bins"];
    return true;
}

static std::map<std::tuple<int,int,int,int>, std::pair<double,double>>
extract_our_helicity_averaged(const json& bins_object) {
    // returns map of (ix,iQ,it,ip) -> (c_avg, err_avg)
    std::map<std::tuple<int,int,int,int>, std::pair<double,double>> out;

    for (auto it = bins_object.begin(); it != bins_object.end(); ++it) {
        int ix=0,iQ=0,itb=0,ip=0;
        if (std::sscanf(it.key().c_str(), "(%d,%d,%d,%d)", &ix,&iQ,&itb,&ip) != 4) continue;

        const json& cell = it.value();
        // Required fields in our contamination JSON (written by your analysis):
        //  N_data.helicity.+1, N_data.helicity.-1
        //  contamination.+1.value, contamination.+1.err, contamination.-1.value, contamination.-1.err

        long long Np = 0, Nm = 0;
        double cp = 0.0, cm = 0.0, ep = 0.0, em = 0.0;

        try { Np = cell.at("N_data").at("helicity").at("+1").get<long long>(); } catch (...) {}
        try { Nm = cell.at("N_data").at("helicity").at("-1").get<long long>(); } catch (...) {}
        try { cp = cell.at("contamination").at("+1").at("value").get<double>(); } catch (...) {}
        try { cm = cell.at("contamination").at("-1").at("value").get<double>(); } catch (...) {}
        try { ep = cell.at("contamination").at("+1").at("err").get<double>();   } catch (...) {}
        try { em = cell.at("contamination").at("-1").at("err").get<double>();   } catch (...) {}

        const double Ntot = double(Np + Nm);
        if (Ntot <= 0.0) { out[{ix,iQ,itb,ip}] = {0.0, 0.0}; continue; }

        const double wp = double(Np) / Ntot;
        const double wm = double(Nm) / Ntot;

        const double c_avg = wp*cp + wm*cm;
        // simple weighted error combination (treating c_plus and c_minus uncorrelated)
        const double e_avg = std::sqrt( (wp*ep)*(wp*ep) + (wm*em)*(wm*em) );

        out[{ix,iQ,itb,ip}] = {c_avg, e_avg};
    }
    return out;
}

// ---------- plotting helpers ----------
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
    g->Draw("PE1 SAME");
    return g;
}

static std::string safe_canvas_name(const std::string& out_png) {
    return fs::path(out_png).filename().string();
}

static double compute_canvas_ymax(bool ratio_mode,
                                  const std::vector<std::pair<double,double>>& Q2s,
                                  const std::vector<std::pair<double,double>>& Ts,
                                  const std::function<bool(int,int,std::vector<double>&,std::vector<double>&,std::vector<double>&)>& fillBoth)
{
    // fillBoth returns (ours, ours_err, lee)
    double ymax = 0.0;
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    std::vector<double> phiC = phiCentersDeg();

    for (int r=0; r<nrows; ++r) {
        for (int ccol=0; ccol<ncols; ++ccol) {
            std::vector<double> A(N_PHI_BINS,0.0), EA(N_PHI_BINS,0.0), B(N_PHI_BINS,0.0);
            (void)fillBoth(ccol, r, A, EA, B);

            if (!ratio_mode) {
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    ymax = std::max(ymax, A[ip] + EA[ip]);
                    ymax = std::max(ymax, B[ip]); // Lee has no errors
                }
            } else {
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const double NL = B[ip];
                    const double NH = A[ip];
                    const double eA = EA[ip];
                    if (NL <= 0.0) continue; // undefined ratio
                    const double R  = (NH <= 0.0) ? 0.0 : NH/NL;
                    double eR = 0.0;
                    if (NL > 0.0) {
                        // eR = eA / NL  == R * (eA / NH) if NH>0
                        eR = (NH > 0.0) ? (R * (eA / NH)) : (eA / NL);
                    }
                    ymax = std::max(ymax, R + eR);
                }
            }
        }
    }

    if (ymax <= 0.0) ymax = ratio_mode ? 1.0 : 0.10;
    if (ratio_mode)  ymax = std::max(ymax, 1.0);
    return ymax;
}

static void draw_one_canvas(const std::string& title,
                            const std::vector<std::pair<double,double>>& Q2s,
                            const std::vector<std::pair<double,double>>& Ts,
                            const std::function<bool(int,int,std::vector<double>&,std::vector<double>&,std::vector<double>&)>& fillBoth,
                            const std::string& out_png,
                            bool draw_ratio_only) {
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows==0 || ncols==0) return;

    const double canvas_ymax = compute_canvas_ymax(draw_ratio_only, Q2s, Ts, fillBoth);

    const int W = 320*ncols + 220;
    const int H = 260*nrows + 260;

    const std::string cname = safe_canvas_name(out_png);
    TCanvas* c = new TCanvas(cname.c_str(), cname.c_str(), W, H);
    c->cd();

    // Top band
    TPad* pTop = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
    pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
    pTop->cd();

    TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42);
    head.SetTextSize(0.18);
    head.DrawLatex(0.50, 0.65, title.c_str());

    // Legend
    std::vector<TObject*> legend_keepalive;
    TLegend* legTop = new TLegend(0.08, 0.10, 0.92, 0.56);
    legTop->SetNColumns(2);
    legTop->SetBorderSize(0);
    legTop->SetFillStyle(0);
    legTop->SetTextFont(42);
    legTop->SetTextSize(0.22);
    if (draw_ratio_only) {
        auto* mRatio = new TMarker(0.0, 0.0, 20);
        mRatio->SetMarkerColor(kBlack);
        auto* lnY1   = new TLine(0.0, 0.0, 1.0, 0.0);
        lnY1->SetLineStyle(2);
        lnY1->SetLineWidth(2);
        lnY1->SetLineColor(kOrange+7);
        legend_keepalive.push_back(mRatio);
        legend_keepalive.push_back(lnY1);
        legTop->AddEntry(mRatio, "Hayward/Lee", "p");
        legTop->AddEntry(lnY1,   "y = 1",       "l");
    } else {
        auto* mH = new TMarker(0.0, 0.0, 20);
        mH->SetMarkerColor(kBlack);
        auto* mL = new TMarker(0.0, 0.0, 24);
        mL->SetMarkerColor(kOrange+7);
        legend_keepalive.push_back(mH);
        legend_keepalive.push_back(mL);
        legTop->AddEntry(mH, "Hayward (pass-2)", "p");
        legTop->AddEntry(mL, "Lee (pass-1)",     "p");
    }
    legTop->Draw();

    // Grid
    c->cd();
    TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.90);
    pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.00, 0.00);

    const auto phiC = phiCentersDeg();

    for (int r=0; r<nrows; ++r) {
        for (int ccol=0; ccol<ncols; ++ccol) {
            pGrid->cd(r*ncols + ccol + 1);
            gPad->SetGrid(1,1);
            gPad->SetTicks(1,1);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.08);

            std::vector<double> A(N_PHI_BINS,0.0), EA(N_PHI_BINS,0.0), B(N_PHI_BINS,0.0);
            (void)fillBoth(ccol, r, A, EA, B);

            TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, canvas_ymax);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle(draw_ratio_only ? "Hayward / Lee" : "#pi^{0} contamination");
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

            TLatex lab; lab.SetNDC(); lab.SetTextSize(0.070); lab.SetTextAlign(11); lab.SetTextFont(42);
            lab.DrawLatex(0.15, 0.92,
                Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                     Q2s[ccol].first, Q2s[ccol].second, Ts[r].first, Ts[r].second));

            const int black  = kBlack;
            const int orange = kOrange+7;

            if (draw_ratio_only) {
                std::vector<double> x, y, ey;
                x.reserve(N_PHI_BINS); y.reserve(N_PHI_BINS); ey.reserve(N_PHI_BINS);
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const double L = B[ip];     // Lee
                    const double H = A[ip];     // ours
                    const double eH = EA[ip];   // our error
                    if (L <= 0.0) continue;     // undefined
                    const double R  = (H <= 0.0) ? 0.0 : H/L;
                    double eR = 0.0;
                    if (H > 0.0) eR = R * (eH / H); // treat Lee as exact here
                    else if (eH > 0.0) eR = eH / L;
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
                // counts mode: plot ours (with error bars) and Lee (no errors)
                std::vector<double> eLee(N_PHI_BINS, 0.0);
                graph_pe1(phiC, A,  EA,   20, black);
                graph_pe1(phiC, B, eLee,  24, orange);
            }
        }
    }

    c->SaveAs(out_png.c_str());
    delete legTop;
    delete c;
}

// ---------- driver ----------
void plot_pi0_contam_cross_checks(const std::vector<LeeRow>& rows,
                                  const std::string& pi0_combined_json,
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

    // Load our contamination bins from combined JSON (pick Fa18 inb/out)
    json bins_inb, bins_out; std::string name_inb, name_out;
    if (!load_our_contam_groups(pi0_combined_json, bins_inb, bins_out, name_inb, name_out)) {
        warn("Cannot proceed without Fa18 inb/out in pi0 combined JSON.");
        return;
    }
    info("Using pi0 periods: " + name_inb + " and " + name_out);

    // Build axes from Lee rows and map Lee contamination
    AxisSets ax = build_axes_from_rows(rows);
    LeeContamPerBin lee = map_lee_contam_to_indices(rows, ax);

    // Extract our helicity-averaged contamination and errors
    auto ours_inb = extract_our_helicity_averaged(bins_inb);
    auto ours_out = extract_our_helicity_averaged(bins_out);

    // Helper lambda to bind fetches into plotting
    auto make_fillBoth = [&](const std::map<std::tuple<int,int,int,int>, std::pair<double,double>>& ours,
                             const std::map<std::tuple<int,int,int,int>, double>& lee_side,
                             int ix) {
        return [&,ix](int iQcol, int irow,
                      std::vector<double>& A,    // our value
                      std::vector<double>& EA,   // our error
                      std::vector<double>& B     // Lee value
                      )->bool {
            bool any=false;
            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                auto k = std::make_tuple(ix,iQcol,irow,ip);
                double a=0.0, ea=0.0, b=0.0;
                auto itA = ours.find(k); if (itA != ours.end()) { a = itA->second.first; ea = itA->second.second; }
                auto itB = lee_side.find(k); if (itB != lee_side.end()) { b = itB->second; }
                A[ip]  = a;
                EA[ip] = ea;
                B[ip]  = b;
                any |= (a>0.0 || b>0.0);
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
                Form("#pi^{0} contamination: %s   x_{B} #in [%.3g, %.3g]",
                     name_inb.c_str(), ax.xB[ix].first, ax.xB[ix].second);
            const std::string title_ratio  =
                Form("#pi^{0} contamination ratio (Hayward/Lee): %s   x_{B} #in [%.3g, %.3g]",
                     name_inb.c_str(), ax.xB[ix].first, ax.xB[ix].second);

            auto fillBoth_inb = make_fillBoth(ours_inb, lee.inb, ix);

            const std::string f_counts = (fs::path(output_base_dir)/Form("pi0_counts_fa18_inb_xB_%d.png", ix)).string();
            const std::string f_ratio  = (fs::path(output_base_dir)/Form("pi0_ratio_fa18_inb_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, fillBoth_inb, f_counts, /*ratio=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, fillBoth_inb, f_ratio,  /*ratio=*/true);
            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }

        // OUT
        {
            const std::string title_counts =
                Form("#pi^{0} contamination: %s   x_{B} #in [%.3g, %.3g]",
                     name_out.c_str(), ax.xB[ix].first, ax.xB[ix].second);
            const std::string title_ratio  =
                Form("#pi^{0} contamination ratio (Hayward/Lee): %s   x_{B} #in [%.3g, %.3g]",
                     name_out.c_str(), ax.xB[ix].first, ax.xB[ix].second);

            auto fillBoth_out = make_fillBoth(ours_out, lee.out, ix);

            const std::string f_counts = (fs::path(output_base_dir)/Form("pi0_counts_fa18_out_xB_%d.png", ix)).string();
            const std::string f_ratio  = (fs::path(output_base_dir)/Form("pi0_ratio_fa18_out_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, fillBoth_out, f_counts, /*ratio=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, fillBoth_out, f_ratio,  /*ratio=*/true);
            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }
    }
}