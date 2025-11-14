// pi0_contamination.cpp
// ------------------------------------------------------------
// Overall pi0 contamination estimator (STRICT: no aliases, no fallbacks)
// Writes ONLY to an existing header column named exactly:
//   "contamination ratio, <Period Label>"
// For Fa18 Inb this is: "contamination ratio, Fa18 Inb"
//
// If that exact column is missing for a period, we fail fast.
//
// Cuts policy (strict, from project conventions):
//   - open_angle_ep2 > 5 deg, (-t1) < 1.0, pTmiss <= 0.20
//   - 3sigma exclusivity cuts loaded from combined_cuts_json
//   - Topologies: (FD, FD), (CD, FD), (CD, FT)
//
// CSV policy (strict):
//   - Never add columns. Never rename. Fail fast if target column absent.
//
// Plot aesthetics: margins consistent with project (left=0.125).
// ------------------------------------------------------------

#include "pi0_contamination.h"
#include "periods.h" // CANONICAL_PERIODS(): vector of { tree_key, label }

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

#include <algorithm>
#include <array>
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

static inline bool is_debug() { return std::getenv("PI0_CONTAM_DEBUG") != nullptr; }
static inline bool trace_rows() { return std::getenv("PI0_CONTAM_TRACE") != nullptr; }

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
    void set_double(int r, int c, double v) {
        if (r<0 || r>=(int)rows.size() || c<0 || c>=(int)header.size()) return;
        std::ostringstream oss; oss.setf(std::ios::fixed); oss<<std::setprecision(8)<<v;
        rows[r][c]=oss.str();
    }
};

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

// ---- title-casing helper (deterministic) ----
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

// -------- CSV row materialization (strict) --------
struct CsvRow {
    int row_index=-1;
    double xb_min=0, xb_max=0, q2_min=0, q2_max=0, tab_min=0, tab_max=0;
    double phimin=0, phimax=360.0, phiavg=std::numeric_limits<double>::quiet_NaN();
    double xb_avg=std::numeric_limits<double>::quiet_NaN();
    double q2_avg=std::numeric_limits<double>::quiet_NaN();
    double tab_avg=std::numeric_limits<double>::quiet_NaN();
    long long n_dvcs_csv=0;
};

static std::string dvcs_unpol_col(const std::string& topo, const std::string& period_display) {
    std::ostringstream os;
    os << "raw yield, ep->epg, " << topo << ", exp, " << period_display << ", unpol";
    return os.str();
}

static std::vector<CsvRow> materialize_rows_for_period_strict(
    const CsvDoc& csv,
    const std::string& period_display)
{
    const std::string disp = to_title_space(period_display);

    std::vector<std::string> required = {
        "xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max",
        "phimin", "phimax",
        "phiavg, " + disp,
        "xBavg, " + disp,
        "Q2avg, " + disp,
        "t_abs_avg, " + disp,
        dvcs_unpol_col("(FD, FD)", disp),
        dvcs_unpol_col("(CD, FD)", disp),
        dvcs_unpol_col("(CD, FT)", disp)
    };
    require_columns_or_die(csv, required, "period \"" + disp + "\"");

    const int c_xb_min   = csv.col("xBmin");
    const int c_xb_max   = csv.col("xBmax");
    const int c_q2_min   = csv.col("Q2min");
    const int c_q2_max   = csv.col("Q2max");
    const int c_tab_min  = csv.col("t_abs_min");
    const int c_tab_max  = csv.col("t_abs_max");
    const int c_phi_min  = csv.col("phimin");
    const int c_phi_max  = csv.col("phimax");
    const int c_phi_avg  = csv.col("phiavg, " + disp);
    const int c_xb_avg   = csv.col("xBavg, " + disp);
    const int c_q2_avg   = csv.col("Q2avg, " + disp);
    const int c_tab_avg  = csv.col("t_abs_avg, " + disp);

    const int c_fd_fd    = csv.col(dvcs_unpol_col("(FD, FD)", disp));
    const int c_cd_fd    = csv.col(dvcs_unpol_col("(CD, FD)", disp));
    const int c_cd_ft    = csv.col(dvcs_unpol_col("(CD, FT)", disp));

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
        double vfd = csv.as_double(r, c_fd_fd);
        double vcf = csv.as_double(r, c_cd_fd);
        double vct = csv.as_double(r, c_cd_ft);
        if (std::isfinite(vfd) && vfd > 0) sum_unpol += (long long)std::llround(vfd);
        if (std::isfinite(vcf) && vcf > 0) sum_unpol += (long long)std::llround(vcf);
        if (std::isfinite(vct) && vct > 0) sum_unpol += (long long)std::llround(vct);
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
    if (out.empty()) fatal("Cuts JSON contains no DVCS_* or eppi0_* blocks with numeric mean/std.");
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

static void plot_period(
    const std::string& period_dir,
    const std::vector<CsvRow>& rows,
    const std::vector<RowCounts>& cnts,
    const std::string& out_root_dir,
    const std::vector<int>& rows_in_slice,
    double xb_mean_slice,
    double q2_mean_slice,
    double tab_mean_slice,
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

    struct Key { int iQ, it; }; // declare ONCE here
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
    if (!dir_exists(plot_dir)) {
        std::cerr << "[pi0_contamination] WARNING: plot directory missing: " << plot_dir
                  << " (did you run makeOutputDirs()?)\n";
    }

    const int nrows = (int)Qs.size();
    const int ncols = (int)Ts.size();
    const int W = 260*ncols + 140;
    const int H = 230*nrows + 200;

    TCanvas* c = new TCanvas(Form("c_contam_%s_xB%d", period_dir.c_str(), slice_index), "", W, H);
    TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
    pTop->SetFillStyle(0); pTop->Draw();
    TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.0, 1.0, 0.90);
    pGrid->SetFillStyle(0); pGrid->Draw(); pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

    pTop->cd();
    TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextSize(0.135);
    head.DrawLatex(0.5, 0.50,
        Form("pi0 contamination vs phi   |   %s   |   <xB>=%.4f   <Q^{2}>=%.4f (GeV^{2})   <|t|>=%.4f (GeV^{2})",
             period_dir.c_str(), xb_mean_slice, q2_mean_slice, tab_mean_slice));

    for (int rr=0; rr<nrows; ++rr) {
        for (int cc=0; cc<ncols; ++cc) {
            pGrid->cd(rr*ncols+cc+1);
            gPad->SetGrid(1,1);
            gPad->SetTopMargin(0.18);
            gPad->SetBottomMargin(0.14);
            gPad->SetLeftMargin(0.125);
            gPad->SetRightMargin(0.06);

            Key K{rr, cc};
            auto it = by_k.find(K);

            TH1* fr = gPad->DrawFrame(0.0, 0.0, 360.0, 1.0);
            fr->GetXaxis()->SetTitle("phi (deg)");
            fr->GetYaxis()->SetTitle("pi0 contamination");
            fr->GetXaxis()->CenterTitle();
            fr->GetYaxis()->CenterTitle();
            fr->GetXaxis()->SetNdivisions(505);

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
                X.reserve(pts.size()); Y.reserve(pts.size()); eX.assign(pts.size(),0.0); eY.reserve(pts.size());
                for (const auto& p : pts) { X.push_back(p.x); Y.push_back(p.y); eY.push_back(p.ey); }
                TGraphErrors* gr = new TGraphErrors((int)X.size(), X.data(), Y.data(), eX.data(), eY.data());
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.9);
                gr->SetLineWidth(2);
                gr->Draw("P SAME");
            }

            TLatex lab; lab.SetNDC(); lab.SetTextAlign(13); lab.SetTextSize(0.045);
            lab.DrawLatex(0.14, 0.96, "Q^{2} and |t| cell");
        }
    }

    const std::string fpath = (std::filesystem::path(out_root_dir) /
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

// ------------ Core ------------
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
    namespace fs=std::filesystem;

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

    // Preflight: require exact output column to exist (no creation, no aliasing).
    for (const auto& J : jobs) {
        const std::string colname = "contamination ratio, " + J.display;
        const int idx = csv.col(colname);
        if (idx < 0) {
            fatal("CSV missing output column: \"" + colname + "\"");
        }
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
        const std::string period_dir   = [] (const std::string& s){
            std::string out = s; for (char& c : out) if (c==' ') c='_'; return out;
        }(J.display);

        std::vector<CsvRow> rows = materialize_rows_for_period_strict(csv, J.display);
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

        // eppi0 DATA counting
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
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / M_PI);
                bool matched=false;
                for (size_t r=0;r<rows.size();++r) {
                    if (row_accepts(rows[r], xB, Q2, tt, phi_deg)) { counts[r].n_pi0_data++; matched=true; break; }
                }
                if (matched) dbgD.matched++;
            }
        }

        // eppi0 RECO MC counting
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
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / M_PI);
                bool matched=false;
                for (size_t r=0;r<rows.size();++r) {
                    if (row_accepts(rows[r], xB, Q2, tt, phi_deg)) { counts[r].n_pi0_reco++; matched=true; break; }
                }
                if (matched) dbgR.matched++;
            }
        }

        // pi0->DVCS background MC counting (apply DVCS cuts, mc)
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
                const double phi_deg = wrap_deg(b.phi2 * 180.0 / M_PI);
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

        // Write contamination ratios into the EXACT CSV output column (no creation, no aliasing)
        const std::string contam_col = "contamination ratio, " + J.display;
        const int c_contam = csv.col(contam_col);
        if (c_contam < 0) {
            fatal("CSV missing output column: \"" + contam_col + "\"");
        } else {
            size_t wrote=0;
            std::lock_guard<std::mutex> lock(csv_mtx);
            for (size_t i=0;i<rows.size();++i) {
                const double Nb = (double)counts[i].n_pi0_bkg;
                const double Nd = (double)counts[i].n_dvcs_csv;
                const double Ne = (double)counts[i].n_pi0_data;
                const double Nr = (double)counts[i].n_pi0_reco;
                if (Nb>0.0 && Nd>0.0 && Ne>0.0 && Nr>0.0) {
                    const double val = (Nb/Nd)*(Ne/Nr);
                    csv.set_double(rows[i].row_index, c_contam, val);
                    ++wrote;
                }
            }
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                std::cout << "[pi0_contamination] Wrote " << wrote
                          << " values to column \"" << contam_col
                          << "\" for period " << J.display << "\n";
            }
        }

        // Per-xB slice plots (compute indices per slice for this period)
        std::map<std::pair<double,double>, std::vector<int>> slice_rows;
        for (int r=0; r<(int)rows.size(); ++r) {
            slice_rows[{rows[r].xb_min, rows[r].xb_max}].push_back(r);
        }
        int slice_idx=0;
        for (const auto& kv : slice_rows) {
            const auto& idxs = kv.second;
            double xb_m=0, q2_m=0, t_m=0; int n=0;
            for (int ridx : idxs) {
                xb_m += rows[ridx].xb_avg;
                q2_m += rows[ridx].q2_avg;
                t_m  += rows[ridx].tab_avg;
                ++n;
            }
            if (n>0) { xb_m/=n; q2_m/=n; t_m/=n; }
            plot_period(period_dir, rows, counts, out_root_dir, idxs, xb_m, q2_m, t_m, ++slice_idx);
        }
    } // end ip periods

    // Save CSV (atomic)
    if (!csv.save_atomic(dvcs_csv_path)) {
        fatal("Failed to save updated CSV: " + dvcs_csv_path);
    }

    return true;
}