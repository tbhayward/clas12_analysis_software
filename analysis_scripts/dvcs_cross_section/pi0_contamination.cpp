// pi0_contamination.cpp
// ------------------------------------------------------------
// pi0 contamination estimator with tuple CSV writes: "(value, stat, sys)"
// Strict CSV policy: do not add columns; resolve existing columns via aliases.
// Uses label_aliases.h (same helpers as total_counts.cpp).
// ------------------------------------------------------------

#include "pi0_contamination.h"
#include "periods.h"
#include "label_aliases.h"  // provides topology_aliases(), period_aliases(), helicity_aliases(), find_col_alias()

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
#include <TString.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

static constexpr double PI_CONST = 3.14159265358979323846;

static inline bool is_debug() { return std::getenv("PI0_CONTAM_DEBUG") != nullptr; }
static inline bool trace_rows_env() { return std::getenv("PI0_CONTAM_TRACE") != nullptr; }

static int resolve_threads(int max_workers) {
    int threads = 1;
#ifdef _OPENMP
    threads = (max_workers <= 0) ? omp_get_max_threads() : max_workers;
    if (const char* env = std::getenv("OMP_NUM_THREADS")) {
        int env_n = std::max(1, std::atoi(env));
        threads = std::min(threads, env_n);
    }
    if (threads > 5) threads = 5;
    if (threads < 1) threads = 1;
#endif
    return threads;
}

[[noreturn]] static void fatal(const std::string& msg) {
    std::cerr << "[pi0_contamination][FATAL] " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

// ---------- CSV doc ----------
struct CsvDoc {
    std::vector<std::string> header;
    std::map<std::string,int> index;
    std::vector<std::vector<std::string>> rows;

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
    static void write_field(std::ostream& os, const std::string& s) {
        bool needq = s.find(',') != std::string::npos || s.find('"') != std::string::npos;
        if (!needq) { os << s; return; }
        os << '"';
        for (char ch : s) { if (ch=='"') os << "\"\""; else os << ch; }
        os << '"';
    }
    static double toD(const std::string& s) {
        if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
        char* e=nullptr; double v = std::strtod(s.c_str(), &e);
        return (e==s.c_str()) ? std::numeric_limits<double>::quiet_NaN() : v;
    }
    bool load(const std::string& path) {
        std::ifstream fin(path);
        if (!fin.is_open()) { std::cerr << "[pi0_contamination] cannot open CSV: " << path << "\n"; return false; }
        std::string line;
        if (!std::getline(fin, line)) { std::cerr << "[pi0_contamination] empty CSV: " << path << "\n"; return false; }
        header = split_csv_line(line);
        index.clear(); for (int i=0;i<(int)header.size();++i) index[header[i]] = i;
        rows.clear();
        while (std::getline(fin, line)) if (!line.empty()) rows.push_back(split_csv_line(line));
        for (auto& r : rows) r.resize(header.size());
        return true;
    }
    bool save_atomic(const std::string& path) const {
        const std::string tmp = path + ".tmp";
        {
            std::ofstream fout(tmp);
            if (!fout.is_open()) { std::cerr << "[pi0_contamination] cannot write CSV tmp: " << tmp << "\n"; return false; }
            for (size_t i=0;i<header.size();++i) { write_field(fout, header[i]); if (i+1<header.size()) fout<<','; }
            fout << "\n";
            for (const auto& row : rows) {
                for (size_t i=0;i<row.size();++i) { write_field(fout, row[i]); if (i+1<row.size()) fout<<','; }
                fout << "\n";
            }
        }
        std::error_code ec;
        std::filesystem::rename(tmp, path, ec);
        if (ec) {
            std::remove(path.c_str());
            std::filesystem::rename(tmp, path, ec);
            if (ec) { std::cerr << "[pi0_contamination] ERROR: atomic rename failed (" << ec.message() << ")\n"; return false; }
        }
        return true;
    }
    int col(const std::string& name) const {
        auto it=index.find(name); return (it==index.end())? -1 : it->second;
    }
    double as_double(int r, int c) const {
        if (r<0 || r>=(int)rows.size() || c<0 || c>=(int)header.size()) return std::numeric_limits<double>::quiet_NaN();
        return toD(rows[r][c]);
    }
    void set_string(int r, int c, const std::string& s) {
        if (r<0 || r>=(int)rows.size() || c<0 || c>=(int)header.size()) return;
        rows[r][c] = s;
    }
};

// tuple formatter
static inline std::string tuple_string(double v, double sv, double sy, int prec=8) {
    std::ostringstream oss; oss.setf(std::ios::fixed);
    oss << "(" << std::setprecision(prec) << v
        << ", " << std::setprecision(prec) << sv
        << ", " << std::setprecision(prec) << sy << ")";
    return oss.str();
}

static void require_columns_or_die(const CsvDoc& csv, const std::vector<std::string>& names, const std::string& ctx) {
    std::vector<std::string> missing;
    for (const auto& n : names) if (csv.col(n) < 0) missing.push_back(n);
    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "Missing required CSV columns for " << ctx << ":\n";
        for (const auto& n : missing) oss << "  - " << n << "\n";
        fatal(oss.str());
    }
}

// ---- title-casing helper (deterministic)
static std::string title_case_token(const std::string& tok) {
    if (tok.empty()) return tok;
    std::string out = tok;
    out[0] = (char)std::toupper((unsigned char)out[0]);
    for (size_t i=1;i<out.size();++i) out[i] = (char)std::tolower((unsigned char)out[i]);
    return out;
}
static std::string to_title_space(const std::string& p) {
    std::string out; out.reserve(p.size()+2);
    std::string cur;
    for (char c : p) {
        if (c=='_') { if (!cur.empty()) { out += title_case_token(cur); out += ' '; cur.clear(); } }
        else cur.push_back(c);
    }
    if (!cur.empty()) { out += title_case_token(cur); }
    return out;
}

// -------- alias-aware column resolvers --------
static int col_alias(const CsvDoc& csv, const std::vector<std::string>& candidates) {
    for (const auto& n : candidates) {
        int c = csv.col(n);
        if (c >= 0) return c;
    }
    return -1;
}

static int dvcs_unpol_col_alias(const CsvDoc& csv, const std::string& topo_label, const std::string& period_display) {
    std::vector<std::string> names;
    const auto tops = topology_aliases(topo_label);
    const auto pers = period_aliases(period_display);
    const auto hels = helicity_aliases("unpol");
    for (const auto& t : tops) {
        for (const auto& p : pers) {
            for (const auto& h : hels) {
                std::ostringstream os;
                os << "raw yield, ep->epg, " << t << ", exp, " << p << ", " << h;
                names.push_back(os.str());
            }
        }
    }
    return col_alias(csv, names);
}

static int avg_col_alias(const CsvDoc& csv, const std::string& base, const std::string& period_display) {
    std::vector<std::string> names;
    const auto pers = period_aliases(period_display);
    for (const auto& p : pers) {
        names.push_back(base + ", " + p);
    }
    return col_alias(csv, names);
}

// -------- CSV row materialization (alias-aware) --------
struct CsvRow {
    int row_index=-1;
    double xb_min=0, xb_max=0, q2_min=0, q2_max=0, tab_min=0, tab_max=0;
    double phimin=0, phimax=360.0, phiavg=std::numeric_limits<double>::quiet_NaN();
    double xb_avg=std::numeric_limits<double>::quiet_NaN();
    double q2_avg=std::numeric_limits<double>::quiet_NaN();
    double tab_avg=std::numeric_limits<double>::quiet_NaN();
    long long n_dvcs_csv=0;
};

static std::vector<CsvRow> materialize_rows_for_period(
    const CsvDoc& csv,
    const std::string& period_display)
{
    // base scalars (no aliases expected for these names)
    const int c_xb_min   = csv.col("xBmin");
    const int c_xb_max   = csv.col("xBmax");
    const int c_q2_min   = csv.col("Q2min");
    const int c_q2_max   = csv.col("Q2max");
    const int c_tab_min  = csv.col("t_abs_min");
    const int c_tab_max  = csv.col("t_abs_max");
    const int c_phi_min  = csv.col("phimin");
    const int c_phi_max  = csv.col("phimax");

    const std::vector<std::string> need_base = {
        "xBmin","xBmax","Q2min","Q2max","t_abs_min","t_abs_max","phimin","phimax"
    };
    require_columns_or_die(csv, need_base, "base ranges");

    // alias-aware averages that are period-tagged
    const int c_phi_avg  = avg_col_alias(csv, "phiavg", period_display);
    const int c_xb_avg   = avg_col_alias(csv, "xBavg", period_display);
    const int c_q2_avg   = avg_col_alias(csv, "Q2avg", period_display);
    const int c_tab_avg  = avg_col_alias(csv, "t_abs_avg", period_display);
    if (c_phi_avg < 0 || c_xb_avg < 0 || c_q2_avg < 0 || c_tab_avg < 0) {
        fatal("Missing one of period-tagged average columns (phiavg,xBavg,Q2avg,t_abs_avg) for period aliases of \"" + period_display + "\".");
    }

    // alias-aware DVCS counts per topology
    const int c_fd_fd    = dvcs_unpol_col_alias(csv, "(FD, FD)", period_display);
    const int c_cd_fd    = dvcs_unpol_col_alias(csv, "(CD, FD)", period_display);
    const int c_cd_ft    = dvcs_unpol_col_alias(csv, "(CD, FT)", period_display);
    if (c_fd_fd < 0 && c_cd_fd < 0 && c_cd_ft < 0) {
        fatal("Could not resolve any DVCS unpolarized columns via aliases for period \"" + period_display + "\".");
    }

    std::vector<CsvRow> rows;
    rows.reserve(csv.rows.size());
    for (int r=0; r<(int)csv.rows.size(); ++r) {
        CsvRow cr;
        cr.row_index = r;
        cr.xb_min  = csv.as_double(r, c_xb_min);
        cr.xb_max  = csv.as_double(r, c_xb_max);
        cr.q2_min  = csv.as_double(r, c_q2_min);
        cr.q2_max  = csv.as_double(r, c_q2_max);
        cr.tab_min = csv.as_double(r, c_tab_min);
        cr.tab_max = csv.as_double(r, c_tab_max);

        cr.phimin  = csv.as_double(r, c_phi_min);
        cr.phimax  = csv.as_double(r, c_phi_max);
        cr.phiavg  = csv.as_double(r, c_phi_avg);

        cr.xb_avg  = csv.as_double(r, c_xb_avg);
        cr.q2_avg  = csv.as_double(r, c_q2_avg);
        cr.tab_avg = csv.as_double(r, c_tab_avg);

        long long sum_unpol = 0;
        if (c_fd_fd >= 0) {
            double v = csv.as_double(r, c_fd_fd);
            if (std::isfinite(v) && v > 0) sum_unpol += (long long)std::llround(v);
        }
        if (c_cd_fd >= 0) {
            double v = csv.as_double(r, c_cd_fd);
            if (std::isfinite(v) && v > 0) sum_unpol += (long long)std::llround(v);
        }
        if (c_cd_ft >= 0) {
            double v = csv.as_double(r, c_cd_ft);
            if (std::isfinite(v) && v > 0) sum_unpol += (long long)std::llround(v);
        }
        cr.n_dvcs_csv = sum_unpol;

        rows.push_back(cr);
    }
    if (rows.empty()) fatal("DVCS CSV has no data rows.");
    return rows;
}

static inline double wrap_deg(double d) {
    if (!std::isfinite(d)) return std::numeric_limits<double>::quiet_NaN();
    double w = std::fmod(d, 360.0);
    if (w < 0.0) w += 360.0;
    if (w >= 360.0) w = std::nextafter(360.0, 0.0);
    return w;
}

static bool row_accepts(const CsvRow& r, double xB, double Q2, double t_abs, double phi_deg) {
    if (!(xB  >= r.xb_min  && xB  < r.xb_max )) return false;
    if (!(Q2  >= r.q2_min  && Q2  < r.q2_max )) return false;
    if (!(t_abs >= r.tab_min && t_abs < r.tab_max)) return false;

    double p = wrap_deg(phi_deg);
    double a = wrap_deg(r.phimin), b = wrap_deg(r.phimax);
    if (a <= b) return (p >= a && p < b);
    return (p >= a || p < b);
}

// -------- Cuts loader --------
struct Stats { double mean=0.0; double std=0.0; };
using VarCutMap = std::map<std::string, Stats>;
struct CutPair { VarCutMap data, mc; };
using CutTable = std::map<std::string, CutPair>;

static inline std::string to_cased_period_key(const std::string& tree_key) {
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
        const std::string key = it.key();
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
    if (out.empty()) fatal("Cuts JSON contains no usable blocks with numeric mean/std.");
    return out;
}
static inline bool within3(double v, const Stats& s) {
    return (v >= s.mean - 3.0*s.std) && (v <= s.mean + 3.0*s.std);
}
static bool passes_cuts(const VarCutMap& cuts, const std::map<std::string,double>& vals) {
    for (const auto& kv : cuts) {
        auto it = vals.find(kv.first);
        if (it == vals.end()) fatal("Cuts reference variable '" + kv.first + "' not available in tree.");
        if (!within3(it->second, kv.second)) return false;
    }
    return true;
}

// -------- Branch binders --------
struct BinderCommon {
    int detector1=0, detector2=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0;

    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_gamma_gamma=0, xF=0;
    double theta_pi0_pi0=0;

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
    int helicity=0;
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

static bool dir_exists(const std::string& d) {
    std::error_code ec;
    return std::filesystem::exists(d, ec);
}

static void ensure_dir(const std::string& d) {
    std::error_code ec;
    std::filesystem::create_directories(d, ec);
}

static inline bool is_fin(double v) { return std::isfinite(v); }

static double avg_if_finite_or_midbin(double avg, double lo, double hi) {
    if (is_fin(avg)) return avg;
    if (std::isfinite(lo) && std::isfinite(hi)) return 0.5*(lo+hi);
    return std::numeric_limits<double>::quiet_NaN();
}

struct CellMeans {
    double xb=std::numeric_limits<double>::quiet_NaN();
    double q2=std::numeric_limits<double>::quiet_NaN();
    double tab=std::numeric_limits<double>::quiet_NaN();
};

// Compute cell means across the set of row indices
static CellMeans compute_cell_means(const std::vector<int>& ridxs, const std::vector<CsvRow>& rows, double xb_lo, double xb_hi) {
    double xb_m=0, q2_m=0, t_m=0; int nxb=0,nq2=0,nt=0;
    for (int r : ridxs) {
        if (is_fin(rows[r].xb_avg))  { xb_m += rows[r].xb_avg;  ++nxb; }
        if (is_fin(rows[r].q2_avg))  { q2_m += rows[r].q2_avg;  ++nq2; }
        if (is_fin(rows[r].tab_avg)) { t_m  += rows[r].tab_avg; ++nt;  }
    }
    CellMeans cm;
    cm.xb  = (nxb>0) ? (xb_m/nxb) : avg_if_finite_or_midbin(std::numeric_limits<double>::quiet_NaN(), xb_lo, xb_hi);
    cm.q2  = (nq2>0) ? (q2_m/nq2) : std::numeric_limits<double>::quiet_NaN();
    cm.tab = (nt >0) ? (t_m /nt ) : std::numeric_limits<double>::quiet_NaN();
    return cm;
}

struct Key { int iQ, it; };

static void plot_period(
    const std::string& period_dir,
    const std::vector<CsvRow>& rows,
    const std::vector<RowCounts>& cnts,
    const std::string& out_root_dir,
    const std::vector<int>& rows_in_slice,
    double xb_lo, double xb_hi,
    int slice_index)
{
    auto uniq_ranges = [](const std::vector<int>& idxs, const std::vector<CsvRow>& V, char which){
        std::set<std::pair<double,double>> st;
        for (int r : idxs) {
            if (which=='Q') st.emplace(V[r].q2_min,  V[r].q2_max);
            if (which=='t') st.emplace(V[r].tab_min, V[r].tab_max);
        }
        return std::vector<std::pair<double,double>>(st.begin(), st.end());
    };
    const auto Qs = uniq_ranges(rows_in_slice, rows, 'Q');
    const auto Ts = uniq_ranges(rows_in_slice, rows, 't');

    auto find_idx = [](const std::pair<double,double>& r, const std::vector<std::pair<double,double>>& V){
        for (int i=0;i<(int)V.size();++i) if (V[i]==r) return i; return -1;
    };

    auto lessK = [](const Key& a, const Key& b){
        if (a.iQ != b.iQ) return a.iQ < b.iQ;
        return a.it < b.it;
    };
    std::map<Key, std::vector<int>, decltype(lessK)> by_k(lessK);

    for (int ridx : rows_in_slice) {
        Key K{find_idx({rows[ridx].q2_min,rows[ridx].q2_max}, Qs),
              find_idx({rows[ridx].tab_min,rows[ridx].tab_max}, Ts)};
        if (K.iQ>=0 && K.it>=0) by_k[K].push_back(ridx);
    }

    const std::string plot_dir = (std::filesystem::path(out_root_dir) / "contamination_plots" / period_dir).string();
    ensure_dir(plot_dir);

    const int nrows = (int)Qs.size();
    const int ncols = (int)Ts.size();
    const int W = 300*ncols + 160;
    const int H = 260*nrows + 260;

    TCanvas* c = new TCanvas(Form("c_contam_%s_xB_%d", period_dir.c_str(), slice_index), "", W, H);
    TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
    pTop->SetFillStyle(0); pTop->Draw();
    TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.0, 1.0, 0.90);
    pGrid->SetFillStyle(0); pGrid->Draw(); pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

    // Canvas title: PERIOD   xB in [lo, hi]
    pTop->cd();
    TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextSize(0.16);
    head.DrawLatex(0.5, 0.50,
        Form("%s   xB in [%.4f, %.4f]", period_dir.c_str(), xb_lo, xb_hi));

    for (int rr=0; rr<nrows; ++rr) {
        for (int cc=0; cc<ncols; ++cc) {
            pGrid->cd(rr*ncols+cc+1);
            gPad->SetGrid(1,1);
            gPad->SetTopMargin(0.24);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.160);
            gPad->SetRightMargin(0.07);

            TH1* fr = gPad->DrawFrame(0.0, 0.0, 360.0, 1.0);
            fr->GetXaxis()->SetTitle("#phi (deg)");
            fr->GetYaxis()->SetTitle("pi0 contamination");
            fr->GetXaxis()->CenterTitle();
            fr->GetYaxis()->CenterTitle();
            fr->GetXaxis()->SetNdivisions(505);
            fr->GetXaxis()->SetTitleSize(0.070);
            fr->GetYaxis()->SetTitleSize(0.070);
            fr->GetXaxis()->SetLabelSize(0.060);
            fr->GetYaxis()->SetLabelSize(0.060);

            auto it = by_k.find(Key{rr, cc});
            if (it == by_k.end()) continue;

            struct Pt { double x, y, ey; };
            std::vector<Pt> pts; pts.reserve(it->second.size());

            for (int ridx : it->second) {
                const auto& R = rows[ridx];
                const auto& C = cnts[ridx];

                const double Nb = (double)C.n_pi0_bkg;
                const double Nd = (double)C.n_dvcs_csv;
                const double Ne = (double)C.n_pi0_data;
                const double Nr = (double)C.n_pi0_reco;

                double val=std::numeric_limits<double>::quiet_NaN();
                double err=0.0;
                if (Nb>0.0 && Nd>0.0 && Ne>0.0 && Nr>0.0) {
                    val = (Nb/Nd) * (Ne/Nr);
                    auto rel = [](double n){ return (n>0.0) ? 1.0/std::sqrt(n) : 0.0; };
                    const double rel2 = std::pow(rel(Nb),2) + std::pow(rel(Nd),2)
                                      + std::pow(rel(Ne),2) + std::pow(rel(Nr),2);
                    err = (std::isfinite(val) ? val : 0.0) * std::sqrt(rel2);
                } else {
                    // optional trace to diagnose why a bin is empty
                    if (trace_rows_env()) {
                        std::cout << "[pi0_contamination][TRACE_PLOT] skip row="<<R.row_index
                                  << " {Nd="<< (long long)Nd << ", Nr="<<(long long)Nr
                                  << ", Nb="<<(long long)Nb << ", Ne="<<(long long)Ne << "}\n";
                    }
                }

                if (std::isfinite(val)) {
                    double xc = std::isfinite(C.phi_center) ? C.phi_center : R.phiavg;
                    xc = wrap_deg(xc);
                    pts.push_back({xc, val, err});
                }
            }

            std::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b){ return a.x < b.x; });
            if (!pts.empty()) {
                std::vector<double> X, Y, eX, eY;
                X.reserve(pts.size()); Y.reserve(pts.size()); eY.reserve(pts.size());
                for (const auto& p : pts) { X.push_back(p.x); Y.push_back(p.y); eY.push_back(p.ey); }
                eX.assign(pts.size(), 0.0);
                TGraphErrors* gr = new TGraphErrors((int)X.size(), X.data(), Y.data(), eX.data(), eY.data());
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.7);
                gr->SetLineWidth(1);
                gr->Draw("P SAME");
            }

            const auto cm = compute_cell_means(it->second, rows, xb_lo, xb_hi);
            TLatex lab; lab.SetNDC(); lab.SetTextAlign(13); lab.SetTextSize(0.055);
            lab.DrawLatex(0.12, 0.83,
                Form("<xB>=%.4f   <Q^{2}>=%.4f (GeV^{2})   <-t>=%.4f (GeV^{2})",
                     cm.xb, cm.q2, cm.tab));
        }
    }

    const std::string plot_dir = (std::filesystem::path(out_root_dir) / "contamination_plots" / period_dir).string();
    const std::string fpath = (std::filesystem::path(plot_dir) /
        ("plot_contamination_" + period_dir + "_xB_" + std::to_string(slice_index) + ".png")).string();
    c->SaveAs(fpath.c_str());
    delete c;
}

// -------- Diagnostics --------
struct BlockDbg {
    long long entries=0, pass_global=0, pass_topo=0, pass_cuts=0, matched=0;
    std::array<long long,3> by_topo{{0,0,0}};
    void print(const std::string& tag, const std::string& period) const {
        std::cout << "[pi0_contamination]["<<period<<"] " << tag
                  << " entries="<<entries
                  << " pass_global="<<pass_global
                  << " pass_topo="<<pass_topo
                  << " pass_cuts="<<pass_cuts
                  << " matched_rows="<<matched
                  << " topo_counts=("<<by_topo[0]<<","<<by_topo[1]<<","<<by_topo[2]<<")\n";
    }
};

} // end anon ns

bool compute_pi0_contamination_overall(
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::string& combined_cuts_json,
    const std::string& dvcs_csv_path,
    const std::string& out_root_dir,
    int max_workers)
{
    gROOT->SetBatch(kTRUE);
    ROOT::EnableThreadSafety();

    CsvDoc csv;
    if (!csv.load(dvcs_csv_path)) return false;

    const CutTable cuts = load_cuts(combined_cuts_json);

    struct Job { std::string tree_key, display, period_dir; };
    std::vector<Job> jobs;
    for (const auto& P : CANONICAL_PERIODS()) {
        const std::string tree_key = P.tree_key;
        const std::string key_data = tree_key + "_eppi0";
        const std::string key_rec  = tree_key + "_rec_mc";
        const std::string key_bkg  = tree_key + "_bkg";
        if (eppi0DataTrees.count(key_data) && eppi0DataTrees.at(key_data) &&
            eppi0RecMcTrees.count(key_rec) && eppi0RecMcTrees.at(key_rec) &&
            eppi0BkgTrees.count(key_bkg)   && eppi0BkgTrees.at(key_bkg)) {
            jobs.push_back({tree_key, P.label, P.label});
        }
    }
    if (jobs.empty()) {
        std::cerr << "[pi0_contamination] Nothing to do: no periods with all three trees present.\n";
        return true;
    }

    // Ensure the output contamination column exists (do not add)
    for (const auto& J : jobs) {
        const auto per_aliases = period_aliases(J.display);
        bool ok=false;
        for (const auto& p : per_aliases) {
            const std::string colname = "contamination ratio, " + p;
            if (csv.col(colname) >= 0) { ok=true; break; }
        }
        if (!ok) fatal("CSV missing output contamination column for period aliases of \"" + J.display + "\".");
    }

    std::mutex csv_mtx;

    int threads = resolve_threads(max_workers);
#ifdef _OPENMP
    std::cout << "[pi0_contamination] Using " << threads << " threads (hard cap = 5)\n";
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
    for (int ip=0; ip<(int)jobs.size(); ++ip) {
        const auto& J = jobs[ip];
        const std::string period_cased = to_cased_period_key(J.tree_key);

        std::string period_dir = J.display;
        for (char& c : period_dir) if (c==' ') c='_';

        std::vector<CsvRow> rows = materialize_rows_for_period(csv, J.display);
        std::vector<RowCounts> counts(rows.size());
        for (size_t i=0;i<rows.size();++i) {
            counts[i].n_dvcs_csv = rows[i].n_dvcs_csv;
            counts[i].phi_center = rows[i].phiavg;
        }

        auto topo_to_key = [](const std::string& s)->std::string {
            if (s=="(FD, FD)") return "FD_FD";
            if (s=="(CD, FD)") return "CD_FD";
            if (s=="(CD, FT)") return "CD_FT";
            fatal("Unknown topology string: " + s);
            return "";
        };
        auto cuts_for = [&](const char* which, const std::string& topoKey)->const CutPair& {
            const std::string k = std::string(which) + "_" + period_cased + "_" + topoKey;
            auto it = cuts.find(k);
            if (it == cuts.end()) fatal("Missing cuts block: " + k);
            return it->second;
        };

        TTree* tD = getTreeOrDie(eppi0DataTrees, J.tree_key + "_eppi0", "eppi0 DATA");
        TTree* tR = getTreeOrDie(eppi0RecMcTrees, J.tree_key + "_rec_mc", "eppi0 RECO MC");
        TTree* tB = getTreeOrDie(eppi0BkgTrees,   J.tree_key + "_bkg",    "BKG->DVCS MC");

        BlockDbg dbgD, dbgR, dbgB;

        // eppi0 DATA
        {
            BinderEppi0Data b; b.bind(tD);
            const Long64_t N = tD->GetEntries(); dbgD.entries = (long long)N;
            for (Long64_t i=0;i<N;++i) {
                tD->GetEntry(i);
                if (!(b.open_angle_ep2 > 5.0)) continue; dbgD.pass_global++;
                if (!((-b.t1) < 1.0)) continue;
                if (!(b.pTmiss <= 0.20)) continue;

                std::string matched_topo;
                int topo_idx=-1;
                if (b.detector1==1 && b.detector2==1) { matched_topo="(FD, FD)"; topo_idx=0; }
                else if (b.detector1==2 && b.detector2==1) { matched_topo="(CD, FD)"; topo_idx=1; }
                else if (b.detector1==2 && b.detector2==0) { matched_topo="(CD, FT)"; topo_idx=2; }
                else { continue; }
                dbgD.pass_topo++; dbgD.by_topo[(topo_idx<0?0:topo_idx)]++;

                const std::string topoKey = topo_to_key(matched_topo);
                const auto& CP = cuts_for("eppi0", topoKey);
                if (!passes_cuts(CP.data, b.cut_vals())) continue; dbgD.pass_cuts++;

                const double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1);
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / PI_CONST);
                bool matched=false;
                for (size_t r=0;r<rows.size();++r) {
                    if (row_accepts(rows[r], xB, Q2, tt, phi_deg)) { counts[r].n_pi0_data++; matched=true; break; }
                }
                if (matched) dbgD.matched++;
            }
        }

        // eppi0 RECO MC
        {
            BinderMC b; b.bind(tR);
            const Long64_t N = tR->GetEntries(); dbgR.entries = (long long)N;
            for (Long64_t i=0;i<N;++i) {
                tR->GetEntry(i);
                if (!(b.open_angle_ep2 > 5.0)) continue; dbgR.pass_global++;
                if (!((-b.t1) < 1.0)) continue;
                if (!(b.pTmiss <= 0.20)) continue;

                std::string matched_topo;
                int topo_idx=-1;
                if (b.detector1==1 && b.detector2==1) { matched_topo="(FD, FD)"; topo_idx=0; }
                else if (b.detector1==2 && b.detector2==1) { matched_topo="(CD, FD)"; topo_idx=1; }
                else if (b.detector1==2 && b.detector2==0) { matched_topo="(CD, FT)"; topo_idx=2; }
                else { continue; }
                dbgR.pass_topo++; dbgR.by_topo[(topo_idx<0?0:topo_idx)]++;

                const std::string topoKey = topo_to_key(matched_topo);
                const auto& CP = cuts_for("eppi0", topoKey);
                if (!passes_cuts(CP.mc, b.cut_vals())) continue; dbgR.pass_cuts++;

                const double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1);
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / PI_CONST);
                bool matched=false;
                for (size_t r=0;r<rows.size();++r) {
                    if (row_accepts(rows[r], xB, Q2, tt, phi_deg)) { counts[r].n_pi0_reco++; matched=true; break; }
                }
                if (matched) dbgR.matched++;
            }
        }

        // pi0->DVCS BKG MC (apply DVCS mc cuts)
        {
            BinderMC b; b.bind(tB);
            const Long64_t N = tB->GetEntries(); dbgB.entries = (long long)N;
            for (Long64_t i=0;i<N;++i) {
                tB->GetEntry(i);
                if (!(b.open_angle_ep2 > 5.0)) continue; dbgB.pass_global++;
                if (!((-b.t1) < 1.0)) continue;
                if (!(b.pTmiss <= 0.20)) continue;

                std::string matched_topo;
                int topo_idx=-1;
                if (b.detector1==1 && b.detector2==1) { matched_topo="(FD, FD)"; topo_idx=0; }
                else if (b.detector1==2 && b.detector2==1) { matched_topo="(CD, FD)"; topo_idx=1; }
                else if (b.detector1==2 && b.detector2==0) { matched_topo="(CD, FT)"; topo_idx=2; }
                else { continue; }
                dbgB.pass_topo++; dbgB.by_topo[(topo_idx<0?0:topo_idx)]++;

                const std::string topoKey = topo_to_key(matched_topo);
                const auto& CP = cuts_for("DVCS", topoKey);
                if (!passes_cuts(CP.mc, b.cut_vals())) continue; dbgB.pass_cuts++;

                const double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1);
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / PI_CONST);
                bool matched=false;
                for (size_t r=0;r<rows.size();++r) {
                    if (row_accepts(rows[r], xB, Q2, tt, phi_deg)) { counts[r].n_pi0_bkg++; matched=true; break; }
                }
                if (matched) dbgB.matched++;
            }
        }

        long long Nd_sum=0, Ne_sum=0, Nr_sum=0, Nb_sum=0;
        int rows_with_dvcs=0, rows_with_vals=0;
        for (size_t r=0;r<rows.size();++r) {
            Nd_sum += counts[r].n_dvcs_csv;
            Ne_sum += counts[r].n_pi0_data;
            Nr_sum += counts[r].n_pi0_reco;
            Nb_sum += counts[r].n_pi0_bkg;
            if (counts[r].n_dvcs_csv > 0) rows_with_dvcs++;
            if (counts[r].n_pi0_data>0 || counts[r].n_pi0_reco>0 || counts[r].n_pi0_bkg>0) rows_with_vals++;
        }
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            std::cout << "[pi0_contamination] Period " << J.display
                      << " rows="<<rows.size()
                      << " DVCS_sum="<<Nd_sum<<" rows_with_dvcs="<<rows_with_dvcs
                      << " data_sum="<<Ne_sum<<" reco_sum="<<Nr_sum<<" bkg_sum="<<Nb_sum
                      << " rows_with_any_count="<<rows_with_vals << "\n";
            dbgD.print("DATA", J.display);
            dbgR.print("RECO_MC", J.display);
            dbgB.print("BKG_MC", J.display);
        }

        // Where Ne>0 but any of Nd/Nr/Nb==0, print a brief row-level reason
        if (trace_rows_env()) {
            for (size_t i=0;i<rows.size();++i) {
                const auto& R = rows[i];
                const auto& C = counts[i];
                if (C.n_pi0_data>0 && (C.n_dvcs_csv==0 || C.n_pi0_reco==0 || C.n_pi0_bkg==0)) {
                    std::cout << "[pi0_contamination][TRACE_MISS] row="<<R.row_index
                              << " Ne>0 but Nd="<<C.n_dvcs_csv
                              << " Nr="<<C.n_pi0_reco
                              << " Nb="<<C.n_pi0_bkg << "\n";
                }
            }
        }

        // Write contamination tuples
        const auto per_aliases = period_aliases(J.display);
        int c_contam = -1;
        for (const auto& p : per_aliases) {
            const std::string colname = "contamination ratio, " + p;
            c_contam = csv.col(colname);
            if (c_contam >= 0) break;
        }
        if (c_contam < 0) fatal("CSV missing output contamination column for period \"" + J.display + "\".");

        size_t wrote=0;
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            std::lock_guard<std::mutex> lock(csv_mtx);
            for (size_t i=0;i<rows.size();++i) {
                const auto& R = rows[i];
                const auto& C = counts[i];

                const double Nb = (double)C.n_pi0_bkg;
                const double Nd = (double)C.n_dvcs_csv;
                const double Ne = (double)C.n_pi0_data;
                const double Nr = (double)C.n_pi0_reco;

                if (Nb>0.0 && Nd>0.0 && Ne>0.0 && Nr>0.0) {
                    const double val = (Nb/Nd) * (Ne/Nr);
                    auto rel = [](double n){ return (n>0.0) ? 1.0/std::sqrt(n) : 0.0; };
                    const double rel2 = std::pow(rel(Nb),2) + std::pow(rel(Nd),2)
                                      + std::pow(rel(Ne),2) + std::pow(rel(Nr),2);
                    const double estat = val * std::sqrt(rel2);
                    const double esys  = 0.0;
                    csv.set_string(R.row_index, c_contam, tuple_string(val, estat, esys, 8));
                    ++wrote;
                } else {
                    // leave blank; optional noisy trace is handled above
                }
            }
            if (!csv.save_atomic(dvcs_csv_path)) {
                fatal("Failed to save updated CSV: " + dvcs_csv_path);
            }
        }

#ifdef _OPENMP
        #pragma omp critical
#endif
        std::cout << "[pi0_contamination] Wrote " << wrote
                  << " tuple values to contamination column for period " << J.display << " (saved)\n";

        // Group rows by xB slice and plot
        std::map<std::pair<double,double>, std::vector<int>> slice_rows;
        for (int r=0; r<(int)rows.size(); ++r) {
            slice_rows[{rows[r].xb_min, rows[r].xb_max}].push_back(r);
        }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            int slice_counter_local = 0;
            for (const auto& kv : slice_rows) {
                const auto& idxs = kv.second;
                const double xb_lo = kv.first.first;
                const double xb_hi = kv.first.second;

                ++slice_counter_local;
                plot_period(period_dir, rows, counts, out_root_dir, idxs,
                            xb_lo, xb_hi, slice_counter_local);
            }
        }
    } // end ip

    return true;
}