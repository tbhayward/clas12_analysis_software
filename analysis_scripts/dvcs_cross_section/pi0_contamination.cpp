// pi0_contamination.cpp
// ------------------------------------------------------------
// Overall (helicity-averaged) pi0 contamination estimator (STRICT, fail-fast)
//   c = (N_bkg_mc / N_dvcs_csv) * (N_pi0_data / N_pi0_reco_mc)
// with Poisson errors. DVCS counts come from CSV; eppi0 counts are re-counted.
//
// Cuts policy:
//   - Global kinematic gates: open_angle_ep2 > 5 deg, |t1| < 1.0, pTmiss <= 0.20
//   - Topology must be one of: (FD, FD), (CD, FD), (CD, FT)  [hard-coded]
//   - 3sigma exclusivity windows:
//       * eppi0 DATA counts:   use eppi0_*_*  "data" block
//       * eppi0 RECO MC:       use eppi0_*_*  "mc"   block
//       * pi0->DVCS BKG MC:    use DVCS_*_*   "mc"   block
//
// Binning & phi handling:
//   - All bins come from the Lee CSV rows (no fixed N_PHI_BINS).
//   - Each CSV row defines xBmin/xBmax, Q2min/Q2max, t_abs_min/max, phimin/phimax (or phiavg fallback).
//   - DVCS counts are the sum over (FD, FD), (CD, FD), (CD, FT) for the period, helicity-aggregated (unpol).
//
// Output:
//   - Plots only: <out_root_dir>/contamination_plots/<PeriodDir>/plot_contamination_<PeriodDir>_xB_<idx>.png
// ------------------------------------------------------------

#include "pi0_contamination.h"
#include "periods.h"            // PeriodDef, CANONICAL_PERIODS()

#include <nlohmann/json.hpp>

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include <TROOT.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

[[noreturn]] static void fatal(const std::string& msg) {
    std::cerr << "[pi0_contamination][FATAL] " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

static inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return (char)std::tolower(c); });
    return s;
}

static inline std::string trim(std::string s) {
    size_t a=0, b=s.size();
    while (a<b && std::isspace((unsigned char)s[a])) ++a;
    while (b>a && std::isspace((unsigned char)s[b-1])) --b;
    return s.substr(a, b-a);
}

static inline std::string period_dir_label(std::string L) {
    for (char& c : L) if (c==' ') c = '_';
    return L;
}

static inline bool passesTopology(int det1, int det2, const std::string& topo) {
    if (topo == "(FD, FD)") return (det1==1 && det2==1);
    if (topo == "(CD, FD)") return (det1==2 && det2==1);
    if (topo == "(CD, FT)") return (det1==2 && det2==0);
    return false;
}

static inline bool passesGlobal(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

static inline double wrap_deg(double d) {
    double w = std::fmod(d, 360.0);
    if (w < 0.0) w += 360.0;
    if (w >= 360.0) w = std::nextafter(360.0, 0.0);
    return w;
}

// -------- CSV helpers --------
struct CsvRow {
    double xb_min=0, xb_max=0, q2_min=0, q2_max=0, tab_min=0, tab_max=0;
    double phimin=0, phimax=360.0, phiavg=std::numeric_limits<double>::quiet_NaN();
    long long n_dvcs_csv=0;
};

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out; out.reserve(64);
    std::string cur; cur.reserve(line.size());
    bool inq=false;
    for (char c : line) {
        if (c=='"') { inq=!inq; continue; }
        if (c==',' && !inq) { out.emplace_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    out.emplace_back(cur);
    return out;
}

static int find_col_exact(const std::vector<std::string>& hdr, const std::string& name) {
    for (int i=0;i<(int)hdr.size();++i) if (hdr[i]==name) return i;
    return -1;
}

static int find_col_anycase(const std::vector<std::string>& hdr, const std::vector<std::string>& aliases) {
    for (int i=0;i<(int)hdr.size();++i) {
        std::string h = to_lower(trim(hdr[i]));
        for (const auto& a : aliases) {
            if (h == to_lower(trim(a))) return i;
        }
    }
    return -1;
}

static int find_col_contains_lower(const std::vector<std::string>& hdr, const std::string& needle_lower) {
    for (int i=0;i<(int)hdr.size();++i) {
        std::string h = to_lower(hdr[i]);
        if (h.find(needle_lower) != std::string::npos) return i;
    }
    return -1;
}

static long long to_ll(const std::string& s) {
    if (s.empty()) return 0;
    try { return (long long)std::llround(std::stod(s)); } catch (...) { return 0; }
}
static double to_d(const std::string& s, double def=std::numeric_limits<double>::quiet_NaN()) {
    if (s.empty()) return def;
    try { return std::stod(s); } catch (...) { return def; }
}

static std::string dvcs_unpol_col(const std::string& topo, const std::string& period_display) {
    std::ostringstream os;
    os << "raw yield, ep->epg, " << topo << ", exp, " << period_display << ", unpol";
    return os.str();
}

static std::vector<CsvRow> load_csv_rows_and_dvcs_counts(
    const std::string& csv_path,
    const std::string& period_display)
{
    std::ifstream ifs(csv_path);
    if (!ifs) fatal("Cannot open DVCS CSV: " + csv_path);

    std::string header_line;
    if (!std::getline(ifs, header_line)) fatal("Empty DVCS CSV: " + csv_path);
    auto hdr = split_csv_line(header_line);

    const int c_xb_min = find_col_anycase(hdr, {"xBmin","xbmin"});
    const int c_xb_max = find_col_anycase(hdr, {"xBmax","xbmax"});
    const int c_q2_min = find_col_anycase(hdr, {"Q2min","q2min"});
    const int c_q2_max = find_col_anycase(hdr, {"Q2max","q2max"});
    const int c_tab_min= find_col_anycase(hdr, {"t_abs_min","tabmin","|t|_min","tmin_abs"});
    const int c_tab_max= find_col_anycase(hdr, {"t_abs_max","tabmax","|t|_max","tmax_abs"});

    int c_phi_min = find_col_anycase(hdr, {"phimin","phi_min"});
    int c_phi_max = find_col_anycase(hdr, {"phimax","phi_max"});
    int c_phi_avg = -1;
    {
        // prefer "phiavg, <Period Display>"
        const std::string needle = "phiavg, " + to_lower(period_display);
        for (int i=0;i<(int)hdr.size();++i) if (to_lower(hdr[i]) == needle) { c_phi_avg=i; break; }
        if (c_phi_avg < 0) c_phi_avg = find_col_contains_lower(hdr, "phiavg");
    }

    if (c_xb_min<0||c_xb_max<0||c_q2_min<0||c_q2_max<0||c_tab_min<0||c_tab_max<0) {
        fatal("CSV missing one or more bin-edge columns (xB/Q2/|t|).");
    }
    // We can live without phimin/phimax if phiavg exists (we create a narrow bin)
    if ((c_phi_min<0 || c_phi_max<0) && c_phi_avg<0) {
        fatal("CSV has neither phimin/phimax nor phiavg.");
    }

    // fixed topologies for aggregation
    const std::vector<std::string> topos = {"(FD, FD)","(CD, FD)","(CD, FT)"};
    std::vector<int> dvcs_cols; dvcs_cols.reserve(topos.size());
    for (const auto& topo : topos) {
        const std::string cname = dvcs_unpol_col(topo, period_display);
        int idx = find_col_exact(hdr, cname);
        if (idx < 0) fatal("DVCS CSV missing column: \"" + cname + "\"");
        dvcs_cols.push_back(idx);
    }

    std::vector<CsvRow> rows;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        auto f = split_csv_line(line);
        if ((int)f.size() < (int)hdr.size()) f.resize(hdr.size());

        CsvRow r;
        r.xb_min  = to_d(f[c_xb_min]);
        r.xb_max  = to_d(f[c_xb_max]);
        r.q2_min  = to_d(f[c_q2_min]);
        r.q2_max  = to_d(f[c_q2_max]);
        r.tab_min = to_d(f[c_tab_min]);
        r.tab_max = to_d(f[c_tab_max]);

        if (c_phi_min >= 0 && c_phi_max >= 0) {
            r.phimin = to_d(f[c_phi_min], 0.0);
            r.phimax = to_d(f[c_phi_max], 360.0);
            if (c_phi_avg >= 0) r.phiavg = to_d(f[c_phi_avg], std::numeric_limits<double>::quiet_NaN());
        } else {
            double pav = to_d(f[c_phi_avg], std::numeric_limits<double>::quiet_NaN());
            if (!std::isfinite(pav)) fatal("CSV phiavg value is missing or invalid.");
            r.phiavg = pav;
            r.phimin = pav - 15.0;
            r.phimax = pav + 15.0;
        }

        long long sum_unpol = 0;
        for (int c : dvcs_cols) sum_unpol += to_ll(f[c]);
        r.n_dvcs_csv = sum_unpol;

        rows.push_back(r);
    }
    if (rows.empty()) fatal("DVCS CSV has no data rows.");
    return rows;
}

static bool row_accepts(const CsvRow& r, double xB, double Q2, double t_abs, double phi_deg) {
    if (!(xB  >= r.xb_min  && xB  < r.xb_max )) return false;
    if (!(Q2  >= r.q2_min  && Q2  < r.q2_max )) return false;
    if (!(t_abs >= r.tab_min && t_abs < r.tab_max)) return false;

    double p = wrap_deg(phi_deg);
    double a = wrap_deg(r.phimin), b = wrap_deg(r.phimax);
    if (a <= b) return (p >= a && p < b);
    return (p >= a || p < b); // wrap-around bin
}

// -------- Cuts loader --------
struct Stats { double mean=0.0; double std=0.0; };
using VarCutMap = std::map<std::string, Stats>;   // var -> {mean,std}
struct CutPair { VarCutMap data, mc; };           // per block: "data":{}, "mc":{}
using CutTable = std::map<std::string, CutPair>;  // key: "DVCS_Fa18_Inb_FD_FD", "eppi0_Fa18_Inb_FD_FD"

static inline std::string to_cased_period_key(const std::string& tree_key) {
    // "DVCS_Fa18_inb" -> "Fa18_Inb"
    if (tree_key.rfind("DVCS_", 0) == 0) {
        std::string tail = tree_key.substr(5);
        if (tail.size() >= 4) {
            if (tail.rfind("_inb") == tail.size()-4) return tail.substr(0, tail.size()-4) + "_Inb";
            if (tail.rfind("_out") == tail.size()-4) return tail.substr(0, tail.size()-4) + "_Out";
        }
        return tail;
    }
    return tree_key;
}

static CutTable load_cuts(const std::string& path) {
    std::ifstream in(path);
    if (!in) fatal("Cannot open cuts JSON: " + path);
    nlohmann::json j; in >> j;

    CutTable out;
    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key(); // e.g. "eppi0_Fa18_Inb_FD_FD"
        const auto& blk = it.value();
        CutPair cp;

        auto fill_one = [&](const char* which, VarCutMap& dst){
            if (!blk.contains(which) || !blk[which].is_object()) return;
            for (auto vit = blk[which].begin(); vit != blk[which].end(); ++vit) {
                const std::string vname = vit.key();
                const auto& vs = vit.value();
                if (!vs.contains("mean") || !vs.contains("std")) continue;
                Stats s;
                s.mean = vs["mean"].get<double>();
                s.std  = vs["std"].get<double>();
                if (std::isfinite(s.std) && s.std > 0.0) dst[vname] = s;
            }
        };

        fill_one("data", cp.data);
        fill_one("mc",   cp.mc);
        if (!cp.data.empty() || !cp.mc.empty()) out[key] = std::move(cp);
    }

    if (out.empty()) fatal("Cuts JSON contains no DVCS_* or eppi0_* blocks with numeric mean/std.");
    return out;
}

static inline bool within3(double v, const Stats& s) {
    return (v >= s.mean - 3.0*s.std) && (v <= s.mean + 3.0*s.std);
}

static bool passes_cuts(const VarCutMap& cuts, const std::map<std::string,double>& vals) {
    for (const auto& kv : cuts) {
        auto it = vals.find(kv.first);
        if (it == vals.end()) {
            fatal("Cuts reference variable '" + kv.first + "' not available in tree.");
        }
        if (!within3(it->second, kv.second)) return false;
    }
    return true;
}

// -------- Branch binders --------
struct BinderCommon {
    int detector1=0, detector2=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0;

    // optional variables that may appear in cuts
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_gamma_gamma=0, xF=0;
    double theta_pi0_pi0=0; // appears in eppi0 blocks

    std::set<std::string> have;

    static void bindI(TTree* t, const char* n, int* a) {
        if (!t->GetBranch(n)) fatal(std::string("Missing int branch '") + n + "'");
        t->SetBranchAddress(n, a);
    }
    static void bindD(TTree* t, const char* n, double* a) {
        if (!t->GetBranch(n)) fatal(std::string("Missing double branch '") + n + "'");
        t->SetBranchAddress(n, a);
    }

    void bind_core(TTree* t) {
        bindI(t, "detector1", &detector1);
        bindI(t, "detector2", &detector2);
        bindD(t, "t1",       &t1);
        bindD(t, "open_angle_ep2", &open_angle_ep2);
        bindD(t, "pTmiss",   &pTmiss);
        bindD(t, "x",        &x);
        bindD(t, "Q2",       &Q2);
        bindD(t, "phi2",     &phi2);

        auto tryD = [&](const char* n, double* a){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); have.insert(n); }
        };
        tryD("Emiss2", &Emiss2);
        tryD("Mx2",    &Mx2);
        tryD("Mx2_1",  &Mx2_1);
        tryD("Mx2_2",  &Mx2_2);
        tryD("theta_gamma_gamma", &theta_gamma_gamma);
        tryD("xF",     &xF);
        tryD("theta_pi0_pi0", &theta_pi0_pi0);
    }

    std::map<std::string,double> cut_vals() const {
        std::map<std::string,double> v;
        v["pTmiss"] = pTmiss;
        if (have.count("Emiss2"))            v["Emiss2"] = Emiss2;
        if (have.count("Mx2"))               v["Mx2"]    = Mx2;
        if (have.count("Mx2_1"))             v["Mx2_1"]  = Mx2_1;
        if (have.count("Mx2_2"))             v["Mx2_2"]  = Mx2_2;
        if (have.count("theta_gamma_gamma")) v["theta_gamma_gamma"] = theta_gamma_gamma;
        if (have.count("xF"))                v["xF"] = xF;
        if (have.count("theta_pi0_pi0"))     v["theta_pi0_pi0"] = theta_pi0_pi0;
        return v;
    }
};

struct BinderEppi0Data : public BinderCommon {
    int helicity=0; // not used for aggregation, but we bind to be consistent
    void bind(TTree* t) {
        bind_core(t);
        if (!t->GetBranch("helicity")) fatal("Missing 'helicity' in eppi0 DATA tree.");
        t->SetBranchAddress("helicity", &helicity);
    }
};

struct BinderMC : public BinderCommon {
    void bind(TTree* t) { bind_core(t); }
};

// -------- Period helpers --------
template <typename MapT>
static TTree* getTreeOrDie(const MapT& m, const std::string& key, const char* human) {
    auto it = m.find(key);
    if (it != m.end() && it->second) return it->second;
    std::ostringstream oss;
    oss << "Missing " << human << " tree under key \"" << key << "\". Available: {";
    bool first=true; for (const auto& kv : m) { if (!first) oss << ", "; first=false; oss << "\"" << kv.first << "\""; }
    oss << "}";
    fatal(oss.str());
    return nullptr;
}

// -------- Plotting --------
struct RowCounts {
    long long n_dvcs_csv=0;
    long long n_pi0_data=0;
    long long n_pi0_reco=0;
    long long n_pi0_bkg=0;
    double    phi_center=std::numeric_limits<double>::quiet_NaN();
};

static void ensure_dir(const std::string& d) {
    namespace fs=std::filesystem;
    std::error_code ec;
    fs::create_directories(d, ec);
    if (ec) fatal("Cannot create directory: " + d + " (" + ec.message() + ")");
}

static void plot_period(
    const std::string& period_dir,
    const std::vector<CsvRow>& rows,
    const std::vector<RowCounts>& cnts,
    const std::string& out_root_dir)
{
    // build unique lists for layout from CSV rows
    auto uniq_ranges = [](const std::vector<CsvRow>& V, char which){
        std::set<std::pair<double,double>> st;
        for (const auto& r : V) {
            if (which=='x') st.emplace(r.xb_min,  r.xb_max);
            if (which=='Q') st.emplace(r.q2_min,  r.q2_max);
            if (which=='t') st.emplace(r.tab_min, r.tab_max);
        }
        return std::vector<std::pair<double,double>>(st.begin(), st.end());
    };
    const auto Xs = uniq_ranges(rows,'x');
    const auto Qs = uniq_ranges(rows,'Q');
    const auto Ts = uniq_ranges(rows,'t');

    auto find_idx = [](const std::pair<double,double>& r, const std::vector<std::pair<double,double>>& V){
        for (int i=0;i<(int)V.size();++i) if (V[i]==r) return i; return -1;
    };

    struct Key { int ix, iQ, it; };
    auto lessK = [](const Key& a, const Key& b){
        if (a.ix != b.ix) return a.ix < b.ix;
        if (a.iQ != b.iQ) return a.iQ < b.iQ;
        return a.it < b.it;
    };
    std::map<Key, std::vector<int>, decltype(lessK)> by_k(lessK);

    for (int r=0; r<(int)rows.size(); ++r) {
        Key K{find_idx({rows[r].xb_min,rows[r].xb_max}, Xs),
              find_idx({rows[r].q2_min,rows[r].q2_max}, Qs),
              find_idx({rows[r].tab_min,rows[r].tab_max}, Ts)};
        if (K.ix>=0 && K.iQ>=0 && K.it>=0) by_k[K].push_back(r);
    }

    const std::string plot_dir = (std::filesystem::path(out_root_dir) / "contamination_plots" / period_dir).string();
    ensure_dir(plot_dir);

    for (int ix=0; ix<(int)Xs.size(); ++ix) {
        // find which Q2,t exist for this xB
        std::set<int> Qset, Tset;
        for (const auto& kv : by_k) if (kv.first.ix == ix) { Qset.insert(kv.first.iQ); Tset.insert(kv.first.it); }
        if (Qset.empty() || Tset.empty()) continue;
        std::vector<int> vQ(Qset.begin(), Qset.end());
        std::vector<int> vT(Tset.begin(), Tset.end());

        const int nrows = (int)vQ.size();
        const int ncols = (int)vT.size();
        const int W = 260*ncols + 140;
        const int H = 230*nrows + 180;

        TCanvas* c = new TCanvas(Form("c_contam_%s_xB%d", period_dir.c_str(), ix), "", W, H);
        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->Draw();
        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.0, 1.0, 0.90);
        pGrid->SetFillStyle(0); pGrid->Draw(); pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextSize(0.085);
        head.DrawLatex(0.5, 0.52,
            Form("#pi^{0} contamination vs #phi   |   %s   |   xB in (%.3g, %.3g)",
                 period_dir.c_str(), Xs[ix].first, Xs[ix].second));

        for (int rr=0; rr<nrows; ++rr) {
            for (int cc=0; cc<ncols; ++cc) {
                pGrid->cd(rr*ncols+cc+1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.18);
                gPad->SetBottomMargin(0.14);
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.06);

                const Key K{ix, vQ[rr], vT[cc]};
                auto it = by_k.find(K);
                if (it == by_k.end()) {
                    TH1* fr = gPad->DrawFrame(0.0, 0.0, 360.0, 0.1);
                    fr->GetXaxis()->SetTitle("#phi (deg)");
                    fr->GetYaxis()->SetTitle("#pi^{0} contamination");
                    continue;
                }

                struct Pt { double x, y, ey; };
                std::vector<Pt> pts; pts.reserve(it->second.size());
                double ymax = 0.0;

                for (int ridx : it->second) {
                    const auto& R = rows[ridx];
                    const auto& C = cnts[ridx];

                    const double Nb = (double)C.n_pi0_bkg;
                    const double Nd = (double)C.n_dvcs_csv;
                    const double Ne = (double)C.n_pi0_data;
                    const double Nr = (double)C.n_pi0_reco;

                    double val=0.0, err=0.0;
                    if (Nb>0.0 && Nd>0.0 && Ne>0.0 && Nr>0.0) {
                        val = (Nb/Nd) * (Ne/Nr);
                        auto rel = [](double n){ return (n>0.0) ? 1.0/std::sqrt(n) : 0.0; };
                        const double rel2 = std::pow(rel(Nb),2) + std::pow(rel(Nd),2)
                                          + std::pow(rel(Ne),2) + std::pow(rel(Nr),2);
                        err = val * std::sqrt(rel2);
                    }

                    double xc = std::isfinite(C.phi_center) ? C.phi_center : 0.5*(R.phimin + R.phimax);
                    xc = wrap_deg(xc);

                    pts.push_back({xc, val, err});
                    ymax = std::max(ymax, val+err);
                }

                std::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b){ return a.x < b.x; });
                std::vector<double> X, Y, eX, eY;
                X.reserve(pts.size()); Y.reserve(pts.size()); eX.assign(pts.size(),0.0); eY.reserve(pts.size());
                for (const auto& p : pts) { X.push_back(p.x); Y.push_back(p.y); eY.push_back(p.ey); }

                const double yhi = std::min(1.0, (ymax>0.0 ? ymax*1.25 : 0.10));
                TH1* fr = gPad->DrawFrame(0.0, 0.0, 360.0, yhi);
                fr->GetXaxis()->SetTitle("#phi (deg)");
                fr->GetYaxis()->SetTitle("#pi^{0} contamination");
                fr->GetXaxis()->CenterTitle();
                fr->GetYaxis()->CenterTitle();
                fr->GetXaxis()->SetNdivisions(505);

                TGraphErrors* gr = new TGraphErrors((int)X.size(), X.data(), Y.data(), eX.data(), eY.data());
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.9);
                gr->SetLineWidth(2);
                gr->Draw("P SAME");

                TLatex lab; lab.SetNDC(); lab.SetTextAlign(13); lab.SetTextSize(0.045);
                lab.DrawLatex(0.14, 0.96,
                    Form("Q2 in (%.3g, %.3g),   |t| in (%.3g, %.3g)",
                         Qs[vQ[rr]].first, Qs[vQ[rr]].second, Ts[vT[cc]].first, Ts[vT[cc]].second));
            }
        }

        const std::string fpath = (std::filesystem::path(plot_dir) /
            ("plot_contamination_" + period_dir + "_xB_" + std::to_string(ix) + ".png")).string();
        c->SaveAs(fpath.c_str());
        delete c;
    }
}

// -------- Core --------
} // end anon ns

bool compute_pi0_contamination_overall(
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::string& combined_cuts_json,
    const std::string& dvcs_csv_path,
    const std::string& out_root_dir)
{
    // fixed topologies we always aggregate over
    const std::vector<std::string> TOPO_LIST = {"(FD, FD)","(CD, FD)","(CD, FT)"};

    // load cuts once
    const CutTable cuts = load_cuts(combined_cuts_json);

    // loop periods that we can actually run (have all three trees)
    for (const auto& P : CANONICAL_PERIODS()) {
        const std::string tree_key = P.tree_key;          // e.g. "DVCS_Fa18_inb"
        const std::string display  = P.label;             // e.g. "Fa18 Inb"
        const std::string period_cased = to_cased_period_key(tree_key); // e.g. "Fa18_Inb"
        const std::string period_dir = period_dir_label(display);       // "Fa18_Inb"

        const std::string key_data = tree_key + "_eppi0";
        const std::string key_rec  = tree_key + "_rec_mc";
        const std::string key_bkg  = tree_key + "_bkg";

        // only run if all needed trees exist
        auto itD = eppi0DataTrees.find(key_data);
        auto itR = eppi0RecMcTrees.find(key_rec);
        auto itB = eppi0BkgTrees.find(key_bkg);
        if (itD == eppi0DataTrees.end() || !itD->second ||
            itR == eppi0RecMcTrees.end() || !itR->second ||
            itB == eppi0BkgTrees.end()   || !itB->second) {
            // silently skip periods that do not have the trio of trees
            continue;
        }

        // load CSV rows and aggregate DVCS counts for this period
        const std::vector<CsvRow> csv_rows = load_csv_rows_and_dvcs_counts(dvcs_csv_path, display);
        std::vector<RowCounts> counts(csv_rows.size());
        for (size_t i=0;i<csv_rows.size();++i) {
            counts[i].n_dvcs_csv = csv_rows[i].n_dvcs_csv;
            counts[i].phi_center  = std::isfinite(csv_rows[i].phiavg)
                                  ? csv_rows[i].phiavg
                                  : 0.5*(csv_rows[i].phimin + csv_rows[i].phimax);
        }

        // prepare topology-key helpers for cuts
        auto topo_to_key = [](const std::string& s)->std::string {
            if (s=="(FD, FD)") return "FD_FD";
            if (s=="(CD, FD)") return "CD_FD";
            if (s=="(CD, FT)") return "CD_FT";
            fatal("Unknown topology string: " + s);
            return "";
        };

        auto cuts_for = [&](const char* which, const std::string& topoKey)->const CutPair& {
            // which in {"DVCS","eppi0"}, key is like "DVCS_Fa18_Inb_FD_FD"
            const std::string k = std::string(which) + "_" + period_cased + "_" + topoKey;
            auto it = cuts.find(k);
            if (it == cuts.end()) fatal("Missing cuts block: " + k);
            return it->second;
        };

        // eppi0 DATA counting
        {
            TTree* t = itD->second;
            BinderEppi0Data b; b.bind(t);
            const Long64_t N = t->GetEntries();
            for (Long64_t i=0;i<N;++i) {
                t->GetEntry(i);
                if (!passesGlobal(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

                std::string matched_topo;
                for (const auto& topoStr : TOPO_LIST) {
                    if (passesTopology(b.detector1, b.detector2, topoStr)) { matched_topo = topoStr; break; }
                }
                if (matched_topo.empty()) continue;

                const std::string topoKey = topo_to_key(matched_topo);
                const auto& CP = cuts_for("eppi0", topoKey);
                if (!passes_cuts(CP.data, b.cut_vals())) continue;

                const double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1);
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / M_PI);
                for (size_t r=0;r<csv_rows.size();++r) {
                    if (row_accepts(csv_rows[r], xB, Q2, tt, phi_deg)) { counts[r].n_pi0_data++; break; }
                }
            }
        }

        // eppi0 RECO MC counting
        {
            TTree* t = itR->second;
            BinderMC b; b.bind(t);
            const Long64_t N = t->GetEntries();
            for (Long64_t i=0;i<N;++i) {
                t->GetEntry(i);
                if (!passesGlobal(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

                std::string matched_topo;
                for (const auto& topoStr : TOPO_LIST) {
                    if (passesTopology(b.detector1, b.detector2, topoStr)) { matched_topo = topoStr; break; }
                }
                if (matched_topo.empty()) continue;

                const std::string topoKey = topo_to_key(matched_topo);
                const auto& CP = cuts_for("eppi0", topoKey);
                if (!passes_cuts(CP.mc, b.cut_vals())) continue;

                const double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1);
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / M_PI);
                for (size_t r=0;r<csv_rows.size();++r) {
                    if (row_accepts(csv_rows[r], xB, Q2, tt, phi_deg)) { counts[r].n_pi0_reco++; break; }
                }
            }
        }

        // pi0->DVCS background MC counting (DVCS cuts, mc)
        {
            TTree* t = itB->second;
            BinderMC b; b.bind(t);
            const Long64_t N = t->GetEntries();
            for (Long64_t i=0;i<N;++i) {
                t->GetEntry(i);
                if (!passesGlobal(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

                std::string matched_topo;
                for (const auto& topoStr : TOPO_LIST) {
                    if (passesTopology(b.detector1, b.detector2, topoStr)) { matched_topo = topoStr; break; }
                }
                if (matched_topo.empty()) continue;

                const std::string topoKey = topo_to_key(matched_topo);
                const auto& CP = cuts_for("DVCS", topoKey); // important: DVCS "mc" cuts here
                if (!passes_cuts(CP.mc, b.cut_vals())) continue;

                const double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1);
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / M_PI);
                for (size_t r=0;r<csv_rows.size();++r) {
                    if (row_accepts(csv_rows[r], xB, Q2, tt, phi_deg)) { counts[r].n_pi0_bkg++; break; }
                }
            }
        }

        // plots
        plot_period(period_dir, csv_rows, counts, out_root_dir);
        std::cout << "[pi0_contamination] Plotted " << display << " (" << period_dir << ")\n";
    }

    return true;
}