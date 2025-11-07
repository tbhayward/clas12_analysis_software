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
#include <TMarker.h>
#include <TString.h>
#include <cstdio>

namespace fs = std::filesystem;

// ---------- small utilities ----------
static inline void info(const std::string& s) { std::cout << "[cross] " << s << std::endl; }
static inline void warn(const std::string& s) { std::cout << "[cross][warn] " << s << std::endl; }
static inline std::string slower(std::string s) { for (auto& c : s) c = (char)std::tolower((unsigned char)c); return s; }

static inline bool approx_equal(double a, double b, double abs = 1e-3, double rel = 5e-4) {
    double diff = std::fabs(a - b);
    if (diff <= abs) return true;
    double m = std::max(std::fabs(a), std::fabs(b));
    return diff <= rel * (m > 0.0 ? m : 1.0);
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

static inline int find_index_fuzzy(const std::pair<double,double>& r,
                                   const std::vector<std::pair<double,double>>& v) {
    for (int i=0;i<(int)v.size();++i) {
        if (approx_equal(r.first, v[i].first) && approx_equal(r.second, v[i].second)) return i;
    }
    return -1;
}

// ---------- Lee contamination gathered by numeric bin ----------
using Key3 = std::tuple<int,int,int>; // (ix,iQ,it)

struct LeeSeries {
    // phi (deg) and contamination value (no error)
    std::vector<double> phi;
    std::vector<double> val;
};

struct LeeContamByCell {
    std::map<Key3, LeeSeries> inb;
    std::map<Key3, LeeSeries> out;
};

static LeeContamByCell collect_lee_series(const std::vector<LeeRow>& rows, const AxisSets& ax) {
    LeeContamByCell out;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        const int ix  = find_index_fuzzy(xb, ax.xB);
        if (ix<0) continue;

        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);

        const int iQ = find_index_fuzzy({r.Q2min,r.Q2max}, Q2s);
        const int it = find_index_fuzzy({r.tmin, r.tmax},  Ts);
        if (iQ<0 || it<0) continue;

        const double phi = r.phiavg;
        const Key3 k{ix,iQ,it};

        if (std::isfinite(r.contam_inb) && r.contam_inb>0.0) {
            auto& s = out.inb[k];
            s.phi.push_back(phi);
            s.val.push_back(r.contam_inb);
        }
        if (std::isfinite(r.contam_out) && r.contam_out>0.0) {
            auto& s = out.out[k];
            s.phi.push_back(phi);
            s.val.push_back(r.contam_out);
        }
    }

    // sort by phi for nicer lines
    auto sorter = [](LeeSeries& s) {
        std::vector<int> idx(s.phi.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b){ return s.phi[a] < s.phi[b]; });
        std::vector<double> p, v;
        p.reserve(idx.size()); v.reserve(idx.size());
        for (int i : idx) { p.push_back(s.phi[i]); v.push_back(s.val[i]); }
        s.phi.swap(p); s.val.swap(v);
    };
    for (auto& kv : out.inb) sorter(kv.second);
    for (auto& kv : out.out) sorter(kv.second);
    return out;
}

// ---------- read our JSON, map to Lee bins by numeric ranges ----------
static bool load_our_periods(const std::string& combined_path,
                             json& period_inb, json& period_out,
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
    period_inb = P.at(name_inb);
    period_out = P.at(name_out);
    if (!period_inb.contains("bins") || !period_out.contains("bins")) {
        warn("Chosen periods lack 'bins'."); return false;
    }
    return true;
}

static std::vector<std::pair<double,double>> read_pair_vec(const json& arr) {
    std::vector<std::pair<double,double>> out;
    for (const auto& x : arr) {
        if (x.is_array() && x.size()==2) out.emplace_back(x[0].get<double>(), x[1].get<double>());
    }
    return out;
}
static std::map<int, std::vector<std::pair<double,double>>> read_pair_map_by_ix(const json& obj) {
    std::map<int, std::vector<std::pair<double,double>>> out;
    for (auto it = obj.begin(); it != obj.end(); ++it) {
        int ix = -1;
        try { ix = std::stoi(it.key()); } catch (...) { continue; }
        out[ix] = read_pair_vec(it.value());
    }
    return out;
}

struct OurAxesFromJSON {
    bool ok = false;
    std::vector<std::pair<double,double>> xb;
    std::map<int, std::vector<std::pair<double,double>>> q2_by_ix;
    std::map<int, std::vector<std::pair<double,double>>> t_by_ix;
    int nphi_default = 12;
};

static OurAxesFromJSON try_read_our_axes(const json& period) {
    OurAxesFromJSON A;
    if (!period.contains("axes")) return A;
    const json& ax = period["axes"];

    if (ax.contains("xB") && ax["xB"].is_array()) A.xb = read_pair_vec(ax["xB"]);
    if (ax.contains("Q2_by_ix") && ax["Q2_by_ix"].is_object()) A.q2_by_ix = read_pair_map_by_ix(ax["Q2_by_ix"]);
    if (ax.contains("t_by_ix")  && ax["t_by_ix"].is_object())  A.t_by_ix  = read_pair_map_by_ix(ax["t_by_ix"]);
    if (ax.contains("N_phi")) {
        try { A.nphi_default = ax["N_phi"].get<int>(); } catch (...) {}
    } else if (ax.contains("phi_bins")) {
        try { A.nphi_default = ax["phi_bins"].get<int>(); } catch (...) {}
    }
    A.ok = (!A.xb.empty());
    return A;
}

static double read_phi_center_deg(const json& cell, int ip, int nphi_fallback) {
    // Try common keys for a stored phi center
    const char* keys[] = { "phi_center_deg", "phi_center", "phi_avg_deg", "phi_deg", "phi" };
    for (const char* k : keys) {
        auto it = cell.find(k);
        if (it != cell.end() && it->is_number()) return it->get<double>();
    }
    // Reconstruct from ip
    const double w = 360.0 / double(std::max(1, nphi_fallback));
    return (ip + 0.5) * w;
}

struct OurSeries {
    // our φ, value, and error
    std::vector<double> phi;
    std::vector<double> val;
    std::vector<double> err;
};
struct OurContamByCell {
    std::map<Key3, OurSeries> inb;
    std::map<Key3, OurSeries> out;
};

static bool get_bin_key(const std::string& key, int& ix, int& iQ, int& it, int& ip) {
    // keys like "(ix,iQ,it,ip)"
    return std::sscanf(key.c_str(), "(%d,%d,%d,%d)", &ix,&iQ,&it,&ip) == 4;
}

static OurContamByCell collect_our_series(const json& period, const AxisSets& lee_ax) {
    OurContamByCell out;
    const json& bins = period["bins"];
    const OurAxesFromJSON A = try_read_our_axes(period);

    if (!A.ok) warn("No 'axes' found in period JSON; falling back to index-based mapping (may misalign).");

    for (auto it = bins.begin(); it != bins.end(); ++it) {
        int ix=-1,iQ=-1,itb=-1,ip=-1;
        if (!get_bin_key(it.key(), ix,iQ,itb,ip)) continue;
        const json& cell = it.value();

        // helicity-weighted average and error
        long long Np = 0, Nm = 0;
        double cp = 0.0, cm = 0.0, ep = 0.0, em = 0.0;
        try { Np = cell.at("N_data").at("helicity").at("+1").get<long long>(); } catch (...) {}
        try { Nm = cell.at("N_data").at("helicity").at("-1").get<long long>(); } catch (...) {}
        try { cp = cell.at("contamination").at("+1").at("value").get<double>(); } catch (...) {}
        try { cm = cell.at("contamination").at("-1").at("value").get<double>(); } catch (...) {}
        try { ep = cell.at("contamination").at("+1").at("err").get<double>();   } catch (...) {}
        try { em = cell.at("contamination").at("-1").at("err").get<double>();   } catch (...) {}

        const double Ntot = double(Np + Nm);
        if (Ntot <= 0.0) continue;
        const double wp = double(Np) / Ntot;
        const double wm = double(Nm) / Ntot;
        const double c_avg = wp*cp + wm*cm;
        const double e_avg = std::sqrt((wp*ep)*(wp*ep) + (wm*em)*(wm*em));

        // determine numeric ranges for this JSON bin
        std::pair<double,double> xbR{0,0}, q2R{0,0}, tR{0,0};
        bool ranges_ok = false;

        if (A.ok) {
            if (ix >= 0 && ix < (int)A.xb.size()) {
                xbR = A.xb[ix];
                auto itQ = A.q2_by_ix.find(ix);
                auto itT = A.t_by_ix.find(ix);
                if (itQ != A.q2_by_ix.end() && itT != A.t_by_ix.end()) {
                    const auto& Q2s = itQ->second;
                    const auto& Ts  = itT->second;
                    if (iQ >= 0 && iQ < (int)Q2s.size() &&
                        itb >= 0 && itb < (int)Ts.size()) {
                        q2R = Q2s[iQ];
                        tR  = Ts[itb];
                        ranges_ok = true;
                    }
                }
            }
        }
        if (!ranges_ok) {
            // try per-cell ranges, if present
            try {
                if (cell.contains("ranges")) {
                    const json& R = cell["ranges"];
                    auto rx = R.at("xB");
                    auto rq = R.at("Q2");
                    auto rt = R.at("t");
                    xbR = { rx[0].get<double>(), rx[1].get<double>() };
                    q2R = { rq[0].get<double>(), rq[1].get<double>() };
                    tR  = { rt[0].get<double>(), rt[1].get<double>() };
                    ranges_ok = true;
                }
            } catch (...) {}
        }
        if (!ranges_ok) {
            // last resort: skip (we do not trust index-only mapping)
            continue;
        }

        // map to Lee grid
        const int lix  = find_index_fuzzy(xbR, lee_ax.xB);
        if (lix < 0) continue;
        const auto& Q2sL = lee_ax.Q2_by_ix.at(lix);
        const auto& TsL  = lee_ax.t_by_ix.at(lix);
        const int liQ = find_index_fuzzy(q2R, Q2sL);
        const int lit = find_index_fuzzy(tR,  TsL);
        if (liQ < 0 || lit < 0) continue;

        // φ center
        const int nphi_fallback = A.ok ? A.nphi_default : 12;
        const double phi_deg = read_phi_center_deg(cell, ip, nphi_fallback);

        // stash
        const Key3 k{lix,liQ,lit};
        // choose inb/out container at call site, so just store both here
        // (we will push into the specific side later)
        // To keep code simple, return both series for this function's caller.
        // We'll let the caller decide which side (inb/out) is wanted.
        // To avoid duplication, we store both side values here:
        // Use tags in the cell if present, else assume same c_avg applies to both helicities but same side.
        // In this code path we don't know the side; the caller will pass which side map to fill.
        // So we return a neutral structure and populate side-specific series in the wrapper.
        // For simplicity, we store into both and let the caller pick.
        // (No-op if not used.)
        // We'll use a tiny wrapper below to pick the side.
        // For now, push into a temporary neutral store:
        // We'll just reuse the same function and fill at caller time.
        // To keep this function side-agnostic, we return nothing per-side here.
        // (See wrapper below.)
        // However, to keep code linear, we will directly fill into both 'inb' and 'out' here; the
        // caller will select which one to read when plotting.
        out.inb[k].phi.push_back(phi_deg);
        out.inb[k].val.push_back(c_avg);
        out.inb[k].err.push_back(e_avg);

        out.out[k].phi.push_back(phi_deg);
        out.out[k].val.push_back(c_avg);
        out.out[k].err.push_back(e_avg);
    }

    // sort by φ
    auto sorter = [](OurSeries& s) {
        std::vector<int> idx(s.phi.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b){ return s.phi[a] < s.phi[b]; });
        std::vector<double> p, v, e;
        p.reserve(idx.size()); v.reserve(idx.size()); e.reserve(idx.size());
        for (int i : idx) { p.push_back(s.phi[i]); v.push_back(s.val[i]); e.push_back(s.err[i]); }
        s.phi.swap(p); s.val.swap(v); s.err.swap(e);
    };
    for (auto& kv : out.inb) sorter(kv.second);
    for (auto& kv : out.out) sorter(kv.second);
    return out;
}

// ---------- plotting helpers ----------
static TGraphErrors* graph_xyerr(const std::vector<double>& X,
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

static TGraphErrors* graph_xy(const std::vector<double>& X,
                              const std::vector<double>& Y,
                              int markerStyle, int color) {
    std::vector<double> ex(X.size(), 0.0), ey(X.size(), 0.0);
    auto* g = new TGraphErrors((int)X.size(),
                               const_cast<double*>(X.data()),
                               const_cast<double*>(Y.data()),
                               ex.data(), ey.data());
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);
    g->Draw("PE1 SAME");
    return g;
}

static double ymax_from_series(const std::vector<double>& y, const std::vector<double>* ey=nullptr) {
    double m = 0.0;
    for (size_t i=0;i<y.size();++i) {
        double v = y[i];
        if (ey) v += std::max(0.0, (*ey)[i]);
        if (v > m) m = v;
    }
    return m;
}

static double compute_canvas_ymax_dynamic(
        bool ratio_mode,
        const std::vector<std::pair<double,double>>& Q2s,
        const std::vector<std::pair<double,double>>& Ts,
        const std::function<bool(int,int,
                                 std::vector<double>&,std::vector<double>&,std::vector<double>&, // ours X,Y,E
                                 std::vector<double>&,std::vector<double>&                         // lee X,Y
                                 )>& fetch) {
    double ymax = 0.0;
    for (int r=0; r<(int)Ts.size(); ++r) {
        for (int c=0; c<(int)Q2s.size(); ++c) {
            std::vector<double> xH,yH,eH,xL,yL;
            (void)fetch(c,r,xH,yH,eH,xL,yL);
            if (!ratio_mode) {
                ymax = std::max(ymax, ymax_from_series(yH, &eH));
                ymax = std::max(ymax, ymax_from_series(yL, nullptr));
            } else {
                // ratio at our φ points by nearest Lee φ within tol
                const double tol = 20.0; // deg
                std::vector<double> R, eR;
                for (size_t i=0;i<xH.size();++i) {
                    double phi = xH[i];
                    // find nearest Lee φ
                    double best_d = 1e9; int jbest = -1;
                    for (size_t j=0;j<xL.size();++j) {
                        double d = std::fabs(xL[j]-phi);
                        if (d < best_d) { best_d = d; jbest = (int)j; }
                    }
                    if (jbest>=0 && best_d <= tol && yL[jbest] > 0.0) {
                        double r  = (yH[i] <= 0.0) ? 0.0 : yH[i] / yL[jbest];
                        double er = 0.0;
                        if (yH[i] > 0.0) er = r * (eH[i] / yH[i]);
                        else if (eH[i] > 0.0) er = eH[i] / yL[jbest];
                        R.push_back(r); eR.push_back(er);
                    }
                }
                ymax = std::max(ymax, ymax_from_series(R, &eR));
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
                            const std::function<bool(int,int,
                                                     std::vector<double>&,std::vector<double>&,std::vector<double>&,
                                                     std::vector<double>&,std::vector<double>&)>& fetch,
                            const std::string& out_png,
                            bool draw_ratio_only) {
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows==0 || ncols==0) return;

    const double canvas_ymax = compute_canvas_ymax_dynamic(draw_ratio_only, Q2s, Ts, fetch);

    const int W = 320*ncols + 220;
    const int H = 260*nrows + 260;

    const std::string cname = fs::path(out_png).filename().string();
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
    std::vector<TObject*> keep;
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
        keep.push_back(mRatio); keep.push_back(lnY1);
        legTop->AddEntry(mRatio, "Hayward/Lee (at Hayward #phi)", "p");
        legTop->AddEntry(lnY1,   "y = 1",                          "l");
    } else {
        auto* mH = new TMarker(0.0, 0.0, 20);
        mH->SetMarkerColor(kBlack);
        auto* mL = new TMarker(0.0, 0.0, 24);
        mL->SetMarkerColor(kOrange+7);
        keep.push_back(mH); keep.push_back(mL);
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

    for (int r=0; r<nrows; ++r) {
        for (int col=0; col<ncols; ++col) {
            pGrid->cd(r*ncols + col + 1);
            gPad->SetGrid(1,1);
            gPad->SetTicks(1,1);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.08);

            std::vector<double> xH,yH,eH,xL,yL;
            (void)fetch(col, r, xH, yH, eH, xL, yL);

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
                     Q2s[col].first, Q2s[col].second, Ts[r].first, Ts[r].second));

            const int black  = kBlack;
            const int orange = kOrange+7;

            if (draw_ratio_only) {
                // build ratio at our phi
                const double tol = 20.0; // deg
                std::vector<double> xr, yr, er;
                xr.reserve(yH.size()); yr.reserve(yH.size()); er.reserve(yH.size());
                for (size_t i=0;i<xH.size();++i) {
                    double phi = xH[i];
                    // nearest Lee phi
                    double best_d = 1e9; int jbest = -1;
                    for (size_t j=0;j<xL.size();++j) {
                        double d = std::fabs(xL[j]-phi);
                        if (d < best_d) { best_d = d; jbest = (int)j; }
                    }
                    if (jbest>=0 && best_d <= tol && yL[jbest] > 0.0) {
                        double r  = (yH[i] <= 0.0) ? 0.0 : yH[i] / yL[jbest];
                        double er_i = 0.0;
                        if (yH[i] > 0.0) er_i = r * (eH[i] / yH[i]);
                        else if (eH[i] > 0.0) er_i = eH[i] / yL[jbest];
                        xr.push_back(phi);
                        yr.push_back(r);
                        er.push_back(er_i);
                    }
                }
                graph_xyerr(xr, yr, er, 20, black);
                TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                one->SetLineStyle(2);
                one->SetLineWidth(2);
                one->SetLineColor(orange);
                one->Draw("SAME");
            } else {
                graph_xyerr(xH, yH, eH, 20, black);
                graph_xy(xL, yL, 24, orange);
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

    // Load our contamination periods (Fa18 inb/out)
    json pinb, pout; std::string name_inb, name_out;
    if (!load_our_periods(pi0_combined_json, pinb, pout, name_inb, name_out)) {
        warn("Cannot proceed without Fa18 inb/out in pi0 combined JSON.");
        return;
    }
    info("Using pi0 periods: " + name_inb + " and " + name_out);

    // Build axes from Lee rows and gather series
    AxisSets ax = build_axes_from_rows(rows);
    LeeContamByCell lee = collect_lee_series(rows, ax);

    // Collect our series and map them onto Lee grid by numeric ranges
    OurContamByCell ours_inb = collect_our_series(pinb, ax);
    OurContamByCell ours_out = collect_our_series(pout, ax);

    auto make_fetch = [&](bool inb_side,
                          int ix)->std::function<bool(int,int,
                                                      std::vector<double>&,std::vector<double>&,std::vector<double>&,
                                                      std::vector<double>&,std::vector<double>&)> {
        return [&,ix,inb_side](int iQcol, int irow,
                               std::vector<double>& xH,
                               std::vector<double>& yH,
                               std::vector<double>& eH,
                               std::vector<double>& xL,
                               std::vector<double>& yL)->bool {
            const Key3 k{ix,iQcol,irow};

            // ours
            {
                const auto& M = inb_side ? ours_inb.inb : ours_out.out;
                auto it = M.find(k);
                if (it != M.end()) {
                    xH = it->second.phi;
                    yH = it->second.val;
                    eH = it->second.err;
                } else { xH.clear(); yH.clear(); eH.clear(); }
            }

            // lee
            {
                const auto& ML = inb_side ? lee.inb : lee.out;
                auto it = ML.find(k);
                if (it != ML.end()) {
                    xL = it->second.phi;
                    yL = it->second.val;
                } else { xL.clear(); yL.clear(); }
            }

            return (!xH.empty() || !xL.empty());
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

            auto fetch_inb = make_fetch(true, ix);

            const std::string f_counts = (fs::path(output_base_dir)/Form("pi0_counts_fa18_inb_xB_%d.png", ix)).string();
            const std::string f_ratio  = (fs::path(output_base_dir)/Form("pi0_ratio_fa18_inb_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, fetch_inb, f_counts, /*ratio=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, fetch_inb, f_ratio,  /*ratio=*/true);
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

            auto fetch_out = make_fetch(false, ix);

            const std::string f_counts = (fs::path(output_base_dir)/Form("pi0_counts_fa18_out_xB_%d.png", ix)).string();
            const std::string f_ratio  = (fs::path(output_base_dir)/Form("pi0_ratio_fa18_out_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, fetch_out, f_counts, /*ratio=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, fetch_out, f_ratio,  /*ratio=*/true);
            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }
    }
}