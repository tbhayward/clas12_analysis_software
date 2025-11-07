// pi0_contamination_cross_check.cpp  (ordered-index remap + robust matching + deep debug)
// -----------------------------------------------------------------------------
// - Reads Lee's CSV rows (load_csv.h).
// - Loads your combined pi0 JSON and extracts (ix,iQ,it,ip) indices and hel-avg values.
// - Rebuilds YOUR numeric bin edges from imports/integrated_bin_v2.csv (load_binning_scheme.h)
//   but PRESERVES original index order (no sorting), so JSON indices map 1:1.
// - Remaps your (ix,iQ,it) to Lee's panels by numeric comparison with tolerance.
// - Plots your points at the 12 fixed phi centers implied by ip.
// - Splits Lee’s CSV into inb/out series so sides are compared consistently.
// - Adds extensive debugging controlled by env vars:
//     PI0X_DEBUG=1                -> verbose mapping diagnostics
//     PI0X_ACCEPT_NEAREST=1       -> if strict match fails, accept nearest index within soft tol
//     PI0X_SOFT_REL=0.01          -> relative soft tol for nearest accept (default 1e-2)
// -----------------------------------------------------------------------------

#include "pi0_contamination_cross_check.h"
#include "load_csv.h"              // LeeRow + CSV loader
#include "load_binning_scheme.h"   // Binning + loader for imports/integrated_bin_v2.csv

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
#include <TMarker.h>
#include <TString.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

// ---------- global debug toggles ----------
static bool   g_debug = false;
static bool   g_accept_nearest = false;
static double g_soft_rel = 1e-2;

// ---------- small utilities ----------
static inline void info(const std::string& s) { std::cout << "[cross] " << s << std::endl; }
static inline void warn(const std::string& s) { std::cout << "[cross][warn] " << s << std::endl; }
static inline void dbg (const std::string& s) { if (g_debug) std::cout << "[cross][dbg] " << s << std::endl; }

static inline std::string slower(std::string s){
    for (auto& c: s) c = (char)std::tolower((unsigned char)c);
    return s;
}

static inline std::string fmt_pair(const std::pair<double,double>& p){
    std::ostringstream os; os.setf(std::ios::fixed); os.precision(6);
    os << "[" << p.first << "," << p.second << "]";
    return os.str();
}

static inline std::string fmt_pairs(const std::vector<std::pair<double,double>>& v){
    std::ostringstream os; os << "{";
    for (size_t i=0;i<v.size();++i){ if (i) os << ", "; os << fmt_pair(v[i]); }
    os << "}";
    return os.str();
}

// A bit looser than before to accommodate CSV vs JSON float formatting
static inline bool approx_equal(double a, double b, double abs=1e-4, double rel=2e-3){
    double d = std::fabs(a-b);
    if (d <= abs) return true;
    double m = std::max(std::fabs(a), std::fabs(b));
    return d <= rel * (m>0.0 ? m : 1.0);
}

static inline void degreeTicks(double xmin, double ymin, double xmax, double labelSize){
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}
static void draw_degree_ticks_here(double labelSize){
    degreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), labelSize);
}

// pick Fa18 inb/out from combined JSON period names (avoid *_supp)
static std::string choose_fa18_key(const std::vector<std::string>& keys, bool want_inb){
    const std::string exact = want_inb ? "fa18_inb" : "fa18_out";
    for (const auto& k : keys) if (slower(k)==exact) return k;
    for (const auto& k : keys){
        const std::string kl = slower(k);
        const bool ok = kl.find("fa18")!=std::string::npos &&
                        kl.find("supp")==std::string::npos &&
                        ((want_inb && kl.find("inb")!=std::string::npos) ||
                         (!want_inb && kl.find("out")!=std::string::npos));
        if (ok) return k;
    }
    return std::string();
}

// ---------- build Lee axes from his CSV (order doesn't matter here) ----------
struct AxisSets {
    std::vector<std::pair<double,double>> xB;
    std::map<int, std::vector<std::pair<double,double>>> Q2_by_ix;
    std::map<int, std::vector<std::pair<double,double>>> t_by_ix;
};

static AxisSets build_axes_from_rows(const std::vector<LeeRow>& rows){
    std::set<std::pair<double,double>> xbset;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> q2_by_xb;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> t_by_xb;

    for (const auto& r : rows){
        auto xb = std::make_pair(r.xBmin, r.xBmax);
        xbset.insert(xb);
        q2_by_xb[xb].insert({r.Q2min, r.Q2max});
        t_by_xb[xb].insert({r.tmin,   r.tmax});
    }

    AxisSets ax;
    ax.xB.assign(xbset.begin(), xbset.end());
    for (int ix=0; ix<(int)ax.xB.size(); ++ix){
        const auto& xb = ax.xB[ix];
        ax.Q2_by_ix[ix] = { q2_by_xb[xb].begin(), q2_by_xb[xb].end() };
        ax.t_by_ix[ix]  = {  t_by_xb[xb].begin(),  t_by_xb[xb].end()  };
    }
    return ax;
}

// nearest index helper for debugging and optional soft acceptance
static inline int nearest_index(const std::pair<double,double>& r,
                                const std::vector<std::pair<double,double>>& v,
                                double& best_absdiff_frac_out){
    int best = -1;
    double best_rel = std::numeric_limits<double>::max();
    const double span = std::max(1e-12, std::fabs(r.second - r.first));
    for (int i=0;i<(int)v.size();++i){
        double d1 = std::fabs(r.first  - v[i].first);
        double d2 = std::fabs(r.second - v[i].second);
        double rel = std::max(d1, d2) / span; // relative to r's span
        if (rel < best_rel){ best_rel = rel; best = i; }
    }
    best_absdiff_frac_out = best_rel;
    return best;
}

static inline int find_index_close(const std::pair<double,double>& r,
                                   const std::vector<std::pair<double,double>>& v){
    for (int i=0;i<(int)v.size();++i){
        if (approx_equal(r.first, v[i].first) && approx_equal(r.second, v[i].second)) return i;
    }
    return -1;
}

// ---------- rebuild YOUR numeric axes from integrated_bin_v2.csv (ORDERED) ----------
struct MyAxes {
    std::vector<std::pair<double,double>> xB; // ordered first-appearance
    std::map<int, std::vector<std::pair<double,double>>> Q2_by_ix; // per-ix ordered first-appearance
    std::map<int, std::vector<std::pair<double,double>>> t_by_ix;  // per-ix ordered first-appearance
};

static MyAxes build_my_axes_ordered(const std::vector<Binning>& scheme){
    // preserve order of first appearance for xB
    std::vector<std::pair<double,double>> xB_order;
    auto seen_xb = std::set<std::pair<double,double>>();
    xB_order.reserve(16);

    for (const auto& b : scheme){
        std::pair<double,double> xb = {b.xBmin, b.xBmax};
        if (!seen_xb.count(xb)){
            seen_xb.insert(xb);
            xB_order.push_back(xb);
        }
    }

    MyAxes ax;
    ax.xB = xB_order;

    // For each ix (by order), collect Q2 and |t| ranges in order first seen at that ix
    for (int ix=0; ix<(int)ax.xB.size(); ++ix){
        const auto xb = ax.xB[ix];
        std::vector<std::pair<double,double>> q_order, t_order;
        std::set<std::pair<double,double>> q_seen, t_seen;

        for (const auto& b : scheme){
            if (std::make_pair(b.xBmin,b.xBmax) != xb) continue;

            std::pair<double,double> q = {b.Q2min, b.Q2max};
            std::pair<double,double> t = {std::fabs(b.tmin), std::fabs(b.tmax)};

            if (!q_seen.count(q)){ q_seen.insert(q); q_order.push_back(q); }
            if (!t_seen.count(t)){ t_seen.insert(t); t_order.push_back(t); }
        }
        ax.Q2_by_ix[ix] = q_order;
        ax.t_by_ix[ix]  = t_order;
    }
    return ax;
}

// ---------- read our JSON, pick Fa18 inb/out ----------
static bool load_our_periods(const std::string& combined_path,
                             json& bins_inb, json& bins_out,
                             std::string& name_inb, std::string& name_out){
    std::ifstream f(combined_path);
    if (!f.is_open()){ warn("Cannot open: " + combined_path); return false; }
    json J; f >> J;
    if (!J.contains("periods")){ warn("No 'periods' in combined JSON"); return false; }

    const json& P = J["periods"];
    std::vector<std::string> names; names.reserve(P.size());
    for (auto it = P.begin(); it != P.end(); ++it) names.push_back(it.key());

    name_inb = choose_fa18_key(names, true);
    name_out = choose_fa18_key(names, false);
    if (name_inb.empty() || name_out.empty()){
        warn("Could not locate both Fa18 inb and Fa18 out in combined JSON"); return false;
    }
    if (!P.at(name_inb).contains("bins") || !P.at(name_out).contains("bins")){
        warn("Selected periods lack 'bins'"); return false;
    }
    bins_inb = P.at(name_inb)["bins"];
    bins_out = P.at(name_out)["bins"];
    return true;
}

// ---------- parse our cells (helicity-avg), then REMAP onto Lee grid ----------
static constexpr int N_PHI_BINS = 12;

static inline double phi_center_from_ip(int ip){
    const double w = 360.0 / double(N_PHI_BINS);
    if (ip < 0) ip = 0;
    if (ip >= N_PHI_BINS) ip = N_PHI_BINS - 1;
    return (ip + 0.5) * w;
}

using Key3 = std::tuple<int,int,int>; // (ix_L, iQ_L, it_L)

struct Series {
    std::vector<double> phi;   // deg
    std::vector<double> val;   // contamination
    std::vector<double> err;   // our uncertainty
};

static void accumulate_point(std::map<Key3, Series>& dst, const Key3& k,
                             double phi_deg, double y, double ey){
    auto& s = dst[k];
    s.phi.push_back(phi_deg);
    s.val.push_back(y);
    s.err.push_back(ey);
}

static inline std::string fmt_series(const Series& s){
    std::ostringstream os; os << "N=" << s.phi.size() << " phi=[";
    for (size_t i=0;i<s.phi.size();++i){ if(i) os<<","; os<<s.phi[i]; }
    os << "] val=[";
    for (size_t i=0;i<s.val.size();++i){ if(i) os<<","; os<<s.val[i]; }
    os << "] err=[";
    for (size_t i=0;i<s.err.size();++i){ if(i) os<<","; os<<s.err[i]; }
    os << "]";
    return os.str();
}

// average helper when N_data missing
static inline bool build_helavg(double cp, double cm, double ep, double em,
                                long long Np, long long Nm,
                                double& c_avg, double& e_avg){
    const double Ntot = double(Np + Nm);
    if (Ntot > 0.0){
        const double wp = (Np>0 ? double(Np)/Ntot : 0.0);
        const double wm = (Nm>0 ? double(Nm)/Ntot : 0.0);
        c_avg = wp*cp + wm*cm;
        e_avg = std::sqrt((wp*ep)*(wp*ep) + (wm*em)*(wm*em));
        return true;
    }
    bool have_p = std::isfinite(cp);
    bool have_m = std::isfinite(cm);
    if (have_p && have_m){
        c_avg = 0.5*(cp+cm);
        e_avg = 0.5*std::sqrt(ep*ep + em*em);
        return true;
    }
    if (have_p){ c_avg = cp; e_avg = ep; return true; }
    if (have_m){ c_avg = cm; e_avg = em; return true; }
    return false;
}

static void collect_and_remap_ours(const json& bins_object,
                                   const MyAxes& myAx,
                                   const AxisSets& leeAx,
                                   std::map<Key3, Series>& out_side,
                                   const char* side_tag){
    size_t kept = 0, dropped_bad_index = 0, dropped_nomatch = 0, dropped_empty = 0;
    std::map<int,size_t> kept_by_ixL;

    const int lee_last_ix = (int)leeAx.xB.size() - 1;

    size_t iter_count = 0;
    for (auto it = bins_object.begin(); it != bins_object.end(); ++it, ++iter_count){
        // key is "(ix,iQ,it,ip)"
        int ix=-1,iQ=-1,itb=-1,ip=-1;
        if (std::sscanf(it.key().c_str(), "(%d,%d,%d,%d)", &ix,&iQ,&itb,&ip) != 4) continue;
        const json& cell = it.value();

        long long Np=0, Nm=0;
        double cp=NAN, cm=NAN, ep=0.0, em=0.0;
        try { Np = cell.at("N_data").at("helicity").at("+1").get<long long>(); } catch(...) {}
        try { Nm = cell.at("N_data").at("helicity").at("-1").get<long long>(); } catch(...) {}
        try { cp = cell.at("contamination").at("+1").at("value").get<double>(); } catch(...) {}
        try { cm = cell.at("contamination").at("-1").at("value").get<double>(); } catch(...) {}
        try { ep = cell.at("contamination").at("+1").at("err").get<double>(); } catch(...) {}
        try { em = cell.at("contamination").at("-1").at("err").get<double>(); } catch(...) {}

        double c_avg=0.0, e_avg=0.0;
        if (!build_helavg(cp, cm, ep, em, Np, Nm, c_avg, e_avg)){
            dropped_empty++;
            if (g_debug){
                std::ostringstream os; os << side_tag << " drop(empty) key=" << it.key();
                dbg(os.str());
            }
            continue;
        }

        // numeric ranges from YOUR CSV axes, preserving index order
        if (ix < 0 || ix >= (int)myAx.xB.size()){
            dropped_bad_index++;
            if (g_debug){
                std::ostringstream os; os << side_tag << " drop(bad ix) key=" << it.key()
                                          << " myAx.xB.size=" << myAx.xB.size();
                dbg(os.str());
            }
            continue;
        }
        const auto xbR = myAx.xB[ix];

        const auto itQv = myAx.Q2_by_ix.find(ix);
        const auto itTv = myAx.t_by_ix.find(ix);
        if (itQv==myAx.Q2_by_ix.end() || itTv==myAx.t_by_ix.end()){
            dropped_bad_index++;
            if (g_debug){
                std::ostringstream os; os << side_tag << " drop(missing Q/t vectors) ix=" << ix
                                          << " key=" << it.key();
                dbg(os.str());
            }
            continue;
        }
        const auto& Q2s_my = itQv->second;
        const auto& Ts_my  = itTv->second;
        if (iQ < 0 || iQ >= (int)Q2s_my.size()){
            dropped_bad_index++;
            if (g_debug){
                std::ostringstream os; os << side_tag << " drop(bad iQ) ix=" << ix
                                          << " iQ=" << iQ << " Q2s_my.size=" << Q2s_my.size()
                                          << " key=" << it.key();
                dbg(os.str());
            }
            continue;
        }
        if (itb < 0 || itb >= (int)Ts_my.size()){
            dropped_bad_index++;
            if (g_debug){
                std::ostringstream os; os << side_tag << " drop(bad it) ix=" << ix
                                          << " it=" << itb << " Ts_my.size=" << Ts_my.size()
                                          << " key=" << it.key();
                dbg(os.str());
            }
            continue;
        }
        const auto q2R = Q2s_my[iQ];
        const auto tR  = Ts_my[itb];

        // map to Lee index space by numeric compare
        int ixL = find_index_close(xbR, leeAx.xB);
        if (ixL < 0){
            double frac=-1.0; int near_ix = nearest_index(xbR, leeAx.xB, frac);
            std::ostringstream os;
            os << side_tag << " xb no-strict-match my=" << fmt_pair(xbR)
               << " nearest Lee ix=" << near_ix << " Lee=" << (near_ix>=0?fmt_pair(leeAx.xB[near_ix]):std::string("NA"))
               << " frac_diff=" << frac;
            dbg(os.str());
            if (g_accept_nearest && near_ix>=0 && frac <= g_soft_rel){
                ixL = near_ix;
                dbg(std::string(side_tag) + " xb accepting nearest due to PI0X_ACCEPT_NEAREST");
            } else {
                dropped_nomatch++;
                continue;
            }
        }
        const auto& Q2s_L = leeAx.Q2_by_ix.at(ixL);
        const auto& Ts_L  = leeAx.t_by_ix.at(ixL);
        int iQL = find_index_close(q2R, Q2s_L);
        if (iQL < 0){
            double frac=-1.0; int near_iQ = nearest_index(q2R, Q2s_L, frac);
            std::ostringstream os;
            os << side_tag << " Q2 no-strict-match ixL=" << ixL
               << " my=" << fmt_pair(q2R)
               << " nearest iQ=" << near_iQ
               << " Lee=" << (near_iQ>=0?fmt_pair(Q2s_L[near_iQ]):std::string("NA"))
               << " frac_diff=" << frac;
            dbg(os.str());
            if (g_accept_nearest && near_iQ>=0 && frac <= g_soft_rel){
                iQL = near_iQ;
                dbg(std::string(side_tag) + " Q2 accepting nearest due to PI0X_ACCEPT_NEAREST");
            } else {
                dropped_nomatch++;
                continue;
            }
        }
        int itL = find_index_close(tR, Ts_L);
        if (itL < 0){
            double frac=-1.0; int near_it = nearest_index(tR, Ts_L, frac);
            std::ostringstream os;
            os << side_tag << " t no-strict-match ixL=" << ixL
               << " my=" << fmt_pair(tR)
               << " nearest it=" << near_it
               << " Lee=" << (near_it>=0?fmt_pair(Ts_L[near_it]):std::string("NA"))
               << " frac_diff=" << frac;
            dbg(os.str());
            if (g_accept_nearest && near_it>=0 && frac <= g_soft_rel){
                itL = near_it;
                dbg(std::string(side_tag) + " t accepting nearest due to PI0X_ACCEPT_NEAREST");
            } else {
                dropped_nomatch++;
                continue;
            }
        }

        const double phi_deg = phi_center_from_ip(ip);
        accumulate_point(out_side, Key3{ixL,iQL,itL}, phi_deg, c_avg, e_avg);
        kept++;
        kept_by_ixL[ixL]++;

        // Focused debug for last xB Lee panel
        if (g_debug && ixL == lee_last_ix){
            std::ostringstream os;
            os << side_tag << " mapped-to-last-xB: key=" << it.key()
               << " myXb=" << fmt_pair(xbR)
               << " LeeXb=" << fmt_pair(leeAx.xB[ixL])
               << " myQ2=" << fmt_pair(q2R)
               << " LeeQ2=" << fmt_pair(leeAx.Q2_by_ix.at(ixL)[iQL])
               << " myt="  << fmt_pair(tR)
               << " Leet=" << fmt_pair(leeAx.t_by_ix.at(ixL)[itL])
               << " phi=" << phi_deg
               << " c=" << c_avg << " e=" << e_avg;
            dbg(os.str());
        }
    }

    // sort by phi for nicer lines
    auto sort_by_phi = [](Series& s){
        std::vector<int> idx(s.phi.size());
        for (int i=0;i<(int)idx.size();++i) idx[i]=i;
        std::sort(idx.begin(), idx.end(), [&](int a,int b){ return s.phi[a] < s.phi[b]; });
        std::vector<double> p,v,e; p.reserve(idx.size()); v.reserve(idx.size()); e.reserve(idx.size());
        for (int i: idx){ p.push_back(s.phi[i]); v.push_back(s.val[i]); e.push_back(s.err[i]); }
        s.phi.swap(p); s.val.swap(v); s.err.swap(e);
    };
    for (auto& kv : out_side) sort_by_phi(kv.second);

    // summary
    {
        std::ostringstream os;
        os << side_tag << " collect_and_remap_ours: kept=" << kept
           << " drop_bad_index=" << dropped_bad_index
           << " drop_nomatch=" << dropped_nomatch
           << " drop_empty=" << dropped_empty
           << " unique_cells=" << out_side.size();
        info(os.str());
    }

    if (g_debug){
        // dump how many points per Lee xB panel we produced
        std::ostringstream os;
        os << side_tag << " kept_by_Lee_ix:";
        for (int i=0;i<(int)leeAx.xB.size();++i){
            size_t n = kept_by_ixL.count(i) ? kept_by_ixL[i] : 0;
            os << " ix=" << i << "(" << fmt_pair(leeAx.xB[i]) << ")=" << n;
        }
        dbg(os.str());
    }
}

// ---------- gather Lee series (split inb/out) ----------
struct LeeSeries { std::vector<double> phi; std::vector<double> val; };

static void collect_lee_split(const std::vector<LeeRow>& rows, const AxisSets& ax,
                              std::map<Key3, LeeSeries>& lee_inb,
                              std::map<Key3, LeeSeries>& lee_out){
    for (const auto& r : rows){
        int ix = find_index_close({r.xBmin,r.xBmax}, ax.xB);
        if (ix<0) continue;
        int iQ = find_index_close({r.Q2min,r.Q2max}, ax.Q2_by_ix.at(ix));
        int it = find_index_close({r.tmin, r.tmax }, ax.t_by_ix.at(ix));
        if (iQ<0 || it<0) continue;
        const double phi = r.phiavg;

        if (std::isfinite(r.contam_inb) && r.contam_inb>0.0){
            auto& s = lee_inb[{ix,iQ,it}]; s.phi.push_back(phi); s.val.push_back(r.contam_inb);
        }
        if (std::isfinite(r.contam_out) && r.contam_out>0.0){
            auto& s = lee_out[{ix,iQ,it}]; s.phi.push_back(phi); s.val.push_back(r.contam_out);
        }
    }
    auto sort_series = [](std::map<Key3, LeeSeries>& M){
        for (auto& kv : M){
            auto& s = kv.second;
            std::vector<int> idx(s.phi.size());
            for (int i=0;i<(int)idx.size();++i) idx[i]=i;
            std::sort(idx.begin(), idx.end(), [&](int a,int b){ return s.phi[a] < s.phi[b]; });
            std::vector<double> p,v; p.reserve(idx.size()); v.reserve(idx.size());
            for (int i: idx){ p.push_back(s.phi[i]); v.push_back(s.val[i]); }
            s.phi.swap(p); s.val.swap(v);
        }
    };
    sort_series(lee_inb);
    sort_series(lee_out);

    if (g_debug){
        std::ostringstream os1, os2;
        size_t cnt1=0, cnt2=0;
        for (auto& kv : lee_inb) cnt1 += kv.second.phi.size();
        for (auto& kv : lee_out) cnt2 += kv.second.phi.size();
        os1 << "Lee inb series: cells=" << lee_inb.size() << " total points=" << cnt1;
        os2 << "Lee out series: cells=" << lee_out.size() << " total points=" << cnt2;
        dbg(os1.str());
        dbg(os2.str());
    }
}

// ---------- plotting helpers ----------
static double ymax_series(const std::vector<double>& y, const std::vector<double>* ey=nullptr){
    double m=0.0; for (size_t i=0;i<y.size();++i){ double v=y[i]; if (ey) v += std::max(0.0, (*ey)[i]); m = std::max(m, v); } return m;
}

static double compute_canvas_ymax(bool ratio_mode,
                                  const std::vector<std::pair<double,double>>& Q2s,
                                  const std::vector<std::pair<double,double>>& Ts,
                                  const std::function<void(int,int,
                                                           std::vector<double>&,std::vector<double>&,std::vector<double>&,
                                                           std::vector<double>&,std::vector<double>&)>& fetch){
    double ymax = 0.0;
    for (int r=0; r<(int)Ts.size(); ++r){
        for (int c=0; c<(int)Q2s.size(); ++c){
            std::vector<double> xH,yH,eH,xL,yL;
            fetch(c,r,xH,yH,eH,xL,yL);
            if (!ratio_mode){
                ymax = std::max(ymax, ymax_series(yH, &eH));
                ymax = std::max(ymax, ymax_series(yL, nullptr));
            } else {
                const double tol = 20.0;
                std::vector<double> R, eR;
                for (size_t i=0;i<xH.size();++i){
                    double phi = xH[i];
                    double best=1e9; int jbest=-1;
                    for (size_t j=0;j<xL.size();++j){
                        double d = std::fabs(xL[j]-phi);
                        if (d<best){ best=d; jbest=(int)j; }
                    }
                    if (jbest>=0 && best<=tol && yL[jbest]>0.0){
                        double r = (yH[i]<=0.0) ? 0.0 : yH[i]/yL[jbest];
                        double er=0.0;
                        if (yH[i]>0.0) er = r*(eH[i]/yH[i]);
                        else if (eH[i]>0.0) er = eH[i]/yL[jbest];
                        R.push_back(r); eR.push_back(er);
                    }
                }
                ymax = std::max(ymax, ymax_series(R, &eR));
            }
        }
    }
    if (ymax<=0.0) ymax = ratio_mode ? 1.0 : 0.10;
    if (ratio_mode) ymax = std::max(ymax, 1.0);
    return ymax;
}

static TGraphErrors* draw_points(const std::vector<double>& X,
                                 const std::vector<double>& Y,
                                 const std::vector<double>& EY,
                                 int mstyle, int color){
    std::vector<double> EX(X.size(), 0.0);
    auto* g = new TGraphErrors((int)X.size(),
                               const_cast<double*>(X.data()),
                               const_cast<double*>(Y.data()),
                               EX.data(),
                               const_cast<double*>(EY.data()));
    g->SetMarkerStyle(mstyle);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);
    g->Draw("PE1 SAME");
    return g;
}

static void draw_canvas(const std::string& title,
                        const std::vector<std::pair<double,double>>& Q2s,
                        const std::vector<std::pair<double,double>>& Ts,
                        const std::function<void(int,int,
                                                 std::vector<double>&,std::vector<double>&,std::vector<double>&,
                                                 std::vector<double>&,std::vector<double>&)>& fetch,
                        const std::string& out_png,
                        bool ratio_mode){
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows==0 || ncols==0) return;

    const double ymax = compute_canvas_ymax(ratio_mode, Q2s, Ts, fetch);

    const int W = 320*ncols + 220;
    const int H = 260*nrows + 260;

    TCanvas* c = new TCanvas(fs::path(out_png).filename().string().c_str(), "", W, H);
    c->cd();

    TPad* pTop = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
    pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
    pTop->cd();
    TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.18);
    head.DrawLatex(0.50, 0.65, title.c_str());

    c->cd();
    TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.90);
    pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
    pGrid->cd(); pGrid->Divide(ncols, nrows, 0.00, 0.00);

    for (int r=0; r<nrows; ++r){
        for (int col=0; col<ncols; ++col){
            pGrid->cd(r*ncols + col + 1);
            gPad->SetGrid(1,1);
            gPad->SetTicks(1,1);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.125);
            gPad->SetRightMargin(0.08);

            std::vector<double> xH,yH,eH,xL,yL;
            fetch(col,r,xH,yH,eH,xL,yL);

            TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, ymax);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle(ratio_mode ? "Hayward / Lee" : "#pi^{0} contamination");
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
                Form("Q^{2} in [%.2g, %.2g], -t in [%.2g, %.2g]",
                     Q2s[col].first, Q2s[col].second, Ts[r].first, Ts[r].second));

            if (ratio_mode){
                const double tol = 20.0;
                std::vector<double> xr, yr, er;
                for (size_t i=0;i<xH.size();++i){
                    double best=1e9; int jbest=-1;
                    for (size_t j=0;j<xL.size();++j){
                        double d = std::fabs(xL[j]-xH[i]);
                        if (d<best){ best=d; jbest=(int)j; }
                    }
                    if (jbest>=0 && best<=tol && yL[jbest]>0.0){
                        double r = (yH[i]<=0.0) ? 0.0 : yH[i]/yL[jbest];
                        double er_i = 0.0;
                        if (yH[i]>0.0) er_i = r*(eH[i]/yH[i]);
                        else if (eH[i]>0.0) er_i = eH[i]/yL[jbest];
                        xr.push_back(xH[i]); yr.push_back(r); er.push_back(er_i);
                    }
                }
                draw_points(xr, yr, er, 20, kBlack);
                TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                one->SetLineStyle(2); one->SetLineWidth(2); one->SetLineColor(kOrange+7); one->Draw("SAME");
            } else {
                draw_points(xH, yH, eH, 20, kBlack);        // ours
                std::vector<double> e0(xL.size(),0.0);
                draw_points(xL, yL, e0, 24, kOrange+7);     // Lee
            }
        }
    }

    c->SaveAs(out_png.c_str());
    delete c;
}

// ---------- driver ----------
void plot_pi0_contam_cross_checks(const std::vector<LeeRow>& rows,
                                  const std::string& pi0_combined_json,
                                  const std::string& output_base_dir){
    // read debug envs
    if (const char* s = std::getenv("PI0X_DEBUG"))             g_debug = (s[0] != '0');
    if (const char* s = std::getenv("PI0X_ACCEPT_NEAREST"))    g_accept_nearest = (s[0] != '0');
    if (const char* s = std::getenv("PI0X_SOFT_REL"))          g_soft_rel = std::atof(s);

    std::ostringstream boot;
    boot << "debug=" << g_debug
         << " accept_nearest=" << g_accept_nearest
         << " soft_rel=" << g_soft_rel;
    info(boot.str());

    fs::create_directories(output_base_dir);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetTextFont(42);

    // 1) Build Lee axes and gather Lee series (split by side)
    AxisSets leeAx = build_axes_from_rows(rows);

    info("Lee xB panels = " + std::to_string(leeAx.xB.size()));
    for (int ix=0; ix<(int)leeAx.xB.size(); ++ix){
        dbg("Lee ix="+std::to_string(ix)+" xB="+fmt_pair(leeAx.xB[ix])+
            " Q2="+fmt_pairs(leeAx.Q2_by_ix[ix])+
            " t="+fmt_pairs(leeAx.t_by_ix[ix]));
    }

    std::map<Key3, LeeSeries> lee_inb, lee_out;
    collect_lee_split(rows, leeAx, lee_inb, lee_out);

    // 2) Load our combined JSON, pick Fa18 inb/out
    json bins_inb, bins_out; std::string name_inb, name_out;
    if (!load_our_periods(pi0_combined_json, bins_inb, bins_out, name_inb, name_out)) return;
    info("Using periods: " + name_inb + " and " + name_out);

    // 3) Rebuild YOUR axes from CSV (ORDER PRESERVED)
    const std::string my_csv = "imports/integrated_bin_v2.csv";
    auto my_scheme = load_binning_scheme(my_csv);
    if (my_scheme.empty()){ warn("Your binning CSV parsed empty: "+my_csv); return; }
    MyAxes myAx = build_my_axes_ordered(my_scheme);
    info("My xB count = " + std::to_string(myAx.xB.size()));
    for (int ix=0; ix<(int)myAx.xB.size(); ++ix){
        dbg("My ix="+std::to_string(ix)+" xB="+fmt_pair(myAx.xB[ix])+
            " Q2="+fmt_pairs(myAx.Q2_by_ix[ix])+
            " t="+fmt_pairs(myAx.t_by_ix[ix]));
    }

    // quick check: print last panel ranges side-by-side
    if (!leeAx.xB.empty()){
        int L = (int)leeAx.xB.size()-1;
        dbg(std::string("Lee last xB=")+fmt_pair(leeAx.xB[L])+
            " Q2="+fmt_pairs(leeAx.Q2_by_ix[L])+
            " t="+fmt_pairs(leeAx.t_by_ix[L]));
    }
    if (!myAx.xB.empty()){
        int M = (int)myAx.xB.size()-1;
        dbg(std::string("My last xB=")+fmt_pair(myAx.xB[M])+
            " Q2="+fmt_pairs(myAx.Q2_by_ix[M])+
            " t="+fmt_pairs(myAx.t_by_ix[M]));
    }

    // 4) Collect our points and REMAP (index -> numeric -> Lee indices)
    std::map<Key3, Series> ours_inb, ours_out;
    collect_and_remap_ours(bins_inb, myAx, leeAx, ours_inb, "INB");
    collect_and_remap_ours(bins_out, myAx, leeAx, ours_out, "OUT");

    // Sanity dump: if last Lee panel is empty for ours, try to log any nearby candidates
    const int lee_last_ix = (int)leeAx.xB.size() - 1;
    auto debug_empty_panel = [&](int ix, const std::map<Key3, Series>& ours_map, const char* tag){
        bool any=false;
        for (const auto& kv : ours_map){
            int ixl, iQ, it; std::tie(ixl,iQ,it)=kv.first;
            if (ixl==ix && !kv.second.phi.empty()){ any=true; break; }
        }
        if (!any){
            std::ostringstream os;
            os << "No points for " << tag << " at Lee ix=" << ix
               << " xB=" << fmt_pair(leeAx.xB[ix])
               << " (dumping first 5 cells overall for " << tag << ")";
            warn(os.str());
            int c=0;
            for (const auto& kv : ours_map){
                if (c>=5) break;
                std::ostringstream d;
                int ixl, iQ, it; std::tie(ixl,iQ,it)=kv.first;
                d << "  cell ix=" << ixl << " iQ=" << iQ << " it=" << it
                  << " series: " << fmt_series(kv.second);
                dbg(d.str());
                ++c;
            }
        }
    };

    debug_empty_panel(lee_last_ix, ours_inb, "ours_inb");
    debug_empty_panel(lee_last_ix, ours_out, "ours_out");

    // 5) Plot per Lee xB panel
    for (int ix=0; ix<(int)leeAx.xB.size(); ++ix){
        const auto& Q2s = leeAx.Q2_by_ix.at(ix);
        const auto& Ts  = leeAx.t_by_ix.at(ix);
        if (Q2s.empty() || Ts.empty()) continue;

        auto fetch_inb = [&](int iQ, int it,
                             std::vector<double>& xH, std::vector<double>& yH, std::vector<double>& eH,
                             std::vector<double>& xL, std::vector<double>& yL){
            const Key3 k{ix,iQ,it};
            auto itH = ours_inb.find(k);
            if (itH != ours_inb.end()){ xH = itH->second.phi; yH = itH->second.val; eH = itH->second.err; }
            else { xH.clear(); yH.clear(); eH.clear(); }
            auto itL = lee_inb.find(k);
            if (itL != lee_inb.end()){ xL = itL->second.phi; yL = itL->second.val; }
            else { xL.clear(); yL.clear(); }
        };

        auto fetch_out = [&](int iQ, int it,
                             std::vector<double>& xH, std::vector<double>& yH, std::vector<double>& eH,
                             std::vector<double>& xL, std::vector<double>& yL){
            const Key3 k{ix,iQ,it};
            auto itH = ours_out.find(k);
            if (itH != ours_out.end()){ xH = itH->second.phi; yH = itH->second.val; eH = itH->second.err; }
            else { xH.clear(); yH.clear(); eH.clear(); }
            auto itL = lee_out.find(k);
            if (itL != lee_out.end()){ xL = itL->second.phi; yL = itL->second.val; }
            else { xL.clear(); yL.clear(); }
        };

        const std::string title_inb_counts =
            Form("#pi^{0} contamination: %s   x_{B} in [%.3g, %.3g]",
                 name_inb.c_str(), leeAx.xB[ix].first, leeAx.xB[ix].second);
        const std::string title_inb_ratio =
            Form("#pi^{0} contamination ratio (Hayward/Lee): %s   x_{B} in [%.3g, %.3g]",
                 name_inb.c_str(), leeAx.xB[ix].first, leeAx.xB[ix].second);

        const std::string out_counts_inb = (fs::path(output_base_dir)/Form("pi0_counts_fa18_inb_xB_%d.png", ix)).string();
        const std::string out_ratio_inb  = (fs::path(output_base_dir)/Form("pi0_ratio_fa18_inb_xB_%d.png",  ix)).string();

        draw_canvas(title_inb_counts, Q2s, Ts, fetch_inb, out_counts_inb, /*ratio=*/false);
        draw_canvas(title_inb_ratio,  Q2s, Ts, fetch_inb, out_ratio_inb,  /*ratio=*/true);
        info("Saved: " + out_counts_inb);
        info("Saved: " + out_ratio_inb);

        const std::string title_out_counts =
            Form("#pi^{0} contamination: %s   x_{B} in [%.3g, %.3g]",
                 name_out.c_str(), leeAx.xB[ix].first, leeAx.xB[ix].second);
        const std::string title_out_ratio =
            Form("#pi^{0} contamination ratio (Hayward/Lee): %s   x_{B} in [%.3g, %.3g]",
                 name_out.c_str(), leeAx.xB[ix].first, leeAx.xB[ix].second);

        const std::string out_counts_out = (fs::path(output_base_dir)/Form("pi0_counts_fa18_out_xB_%d.png", ix)).string();
        const std::string out_ratio_out  = (fs::path(output_base_dir)/Form("pi0_ratio_fa18_out_xB_%d.png",  ix)).string();

        draw_canvas(title_out_counts, Q2s, Ts, fetch_out, out_counts_out, /*ratio=*/false);
        draw_canvas(title_out_ratio,  Q2s, Ts, fetch_out, out_ratio_out,  /*ratio=*/true);
        info("Saved: " + out_counts_out);
        info("Saved: " + out_ratio_out);
    }
}