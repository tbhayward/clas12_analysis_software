// bin_means.cpp
// Progress-enabled, thread-safe binding for ROOT. Single pass per TTree.
// Uses ONLY phi2 (radians) -> converts to degrees for CSV phimin/phimax checks.
// Skips any bin with fewer than MIN_BIN_COUNT entries.

#include "bin_means.h"
#include "periods.h"
#include "load_binning_scheme.h"
#include "global_cuts.h"

#include <TTree.h>
#include <TROOT.h>          // ROOT::EnableThreadSafety
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

static constexpr double PI      = 3.14159265358979323846;
static constexpr double RAD2DEG = 180.0 / PI;
static constexpr double DEG2RAD = PI / 180.0;

// Minimum entries required to write a bin's averages
static constexpr long long MIN_BIN_COUNT = 1;

static const std::string kCutsJSON = "output/jsons/combined_cuts.json";

// ---------------- printing helpers ----------------
static inline void print_banner(const std::string& msg) {
    std::cout << "[bin_means] " << msg << std::endl;
}
static inline void print_status_singleline(const std::string& tag,
                                           double pct,
                                           long long i,
                                           long long N,
                                           long long pass_global,
                                           long long pass_3sig,
                                           long long used_rows,
                                           double elapsed_s) {
    double rate = (elapsed_s > 0.0) ? (double)i / elapsed_s : 0.0;
    std::cout << "[bin_means][" << tag << "] "
              << std::fixed << std::setprecision(1)
              << pct << "%  "
              << i << "/" << N
              << "  global=" << pass_global
              << "  excl=" << pass_3sig
              << "  used=" << used_rows
              << "  rate=" << std::setprecision(2) << rate << " ev/s"
              << std::endl;
}

// ---------------- labels and mapping ----------------
enum class Topology { FD_FD, CD_FD, CD_FT };

static inline const char* topo_tag(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }
    return "FD_FD";
}

struct PeriodKeyMap {
    std::string json_tag;      // e.g. "DVCS_Fa18_Out" (matches combined_cuts.json keys)
    std::string csv_label;     // e.g. "Fa18 Out"      (matches CSV column labels)
    std::string period_label;  // e.g. "fa18_out"      (matches global_cuts period labels)
};

static inline PeriodKeyMap period_key_to_tags(const std::string& k) {
    std::string s; s.reserve(k.size());
    for (char c : k) s.push_back((char)std::tolower((unsigned char)c));
    auto has = [&](const char* sub){ return s.find(sub) != std::string::npos; };

    // Deterministic mapping only. No fallbacks.
    if (has("fa18") && has("inb")) return {"DVCS_Fa18_Inb", "Fa18 Inb", "fa18_inb"};
    if (has("fa18") && has("out")) return {"DVCS_Fa18_Out", "Fa18 Out", "fa18_out"};
    if (has("sp19") && has("inb")) return {"DVCS_Sp19_Inb", "Sp19 Inb", "sp19_inb"};
    if (has("sp18") && has("inb")) return {"DVCS_Sp18_Inb", "Sp18 Inb", "sp18_inb"};
    if (has("sp18") && has("out")) return {"DVCS_Sp18_Out", "Sp18 Out", "sp18_out"};

    // Fail-fast: if we cannot map it, do not proceed.
    std::ostringstream ss;
    ss << "[bin_means] FATAL: cannot map period key '" << k
       << "' into {json_tag,csv_label,period_label}.";
    std::cerr << ss.str() << std::endl;
    std::exit(EXIT_FAILURE);

    return {"", "", ""};
}

static inline const std::vector<std::string>& csv_period_labels() {
    static const std::vector<std::string> v = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    return v;
}
static inline const std::vector<std::string>& csv_group_labels() {
    static const std::vector<std::string> v = { "Fa18", "Sp18", "10.6 GeV" };
    return v;
}

// Theta column label sets (must match initialize_pass2_csv.cpp schema exactly)
static inline const std::vector<std::string>& e_theta_labels() {
    static const std::vector<std::string> v = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out", "Fa18", "Sp18", "10.6 GeV"
    };
    return v;
}
static inline const std::vector<std::string>& p_theta_labels() {
    static const std::vector<std::string> v = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out", "Fa18", "Sp18", "10.6 GeV"
    };
    return v;
}
static inline const std::vector<std::string>& g_theta_labels() {
    static const std::vector<std::string> v = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out", "Fa18", "Sp18", "10.6 GeV"
    };
    return v;
}

// ---------------- CSV I/O ----------------
struct CSV {
    std::vector<std::string> header;
    std::unordered_map<std::string,int> index;
    std::vector<std::vector<std::string>> rows;
};

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool inq = false;
    for (char c : line) {
        if (c == '"') { inq = !inq; continue; }
        if (c == ',' && !inq) { out.push_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

static bool load_csv(const std::string& path, CSV& csv) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::cerr << "[bin_means] FATAL: cannot open CSV: " << path << std::endl;
        return false;
    }
    std::string line;
    if (!std::getline(fin, line)) {
        std::cerr << "[bin_means] FATAL: empty CSV: " << path << std::endl;
        return false;
    }
    csv.header = split_csv_line(line);
    csv.index.clear();
    for (int i = 0; i < (int)csv.header.size(); ++i) csv.index[csv.header[i]] = i;

    csv.rows.clear();
    while (std::getline(fin, line)) {
        if (!line.empty()) csv.rows.push_back(split_csv_line(line));
    }
    return true;
}

static bool write_csv(const std::string& path, const CSV& csv) {
    std::ofstream fout(path);
    if (!fout.is_open()) {
        std::cerr << "[bin_means] FATAL: cannot write CSV: " << path << std::endl;
        return false;
    }
    // header
    for (size_t i = 0; i < csv.header.size(); ++i) {
        const std::string& s = csv.header[i];
        bool needq = s.find(',') != std::string::npos || s.find('"') != std::string::npos;
        if (needq) {
            fout << '"';
            for (char ch : s) fout << (ch == '"' ? "\"\"" : std::string(1, ch));
            fout << '"';
        } else {
            fout << s;
        }
        if (i + 1 < csv.header.size()) fout << ',';
    }
    fout << "\n";
    // rows
    for (const auto& row : csv.rows) {
        for (size_t i = 0; i < row.size(); ++i) {
            const std::string& s = row[i];
            bool needq = s.find(',') != std::string::npos || s.find('"') != std::string::npos;
            if (needq) {
                fout << '"';
                for (char ch : s) fout << (ch == '"' ? "\"\"" : std::string(1, ch));
                fout << '"';
            } else {
                fout << s;
            }
            if (i + 1 < row.size()) fout << ',';
        }
        fout << "\n";
    }
    return true;
}

static int col(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);
    if (it == csv.index.end()) {
        std::cerr << "[bin_means] FATAL: column missing: " << name << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return it->second;
}

static inline std::string col_xBavg(const std::string& lab)  { return "xBavg, "      + lab; }
static inline std::string col_Q2avg(const std::string& lab)  { return "Q2avg, "      + lab; }
static inline std::string col_tavg (const std::string& lab)  { return "t_abs_avg, "  + lab; }
static inline std::string col_phiavg(const std::string& lab) { return "phiavg, "     + lab; }

static inline std::string col_e_theta(const std::string& lab) { return "e_theta, " + lab; }
static inline std::string col_p_theta(const std::string& lab) { return "p_theta, " + lab; }
static inline std::string col_g_theta(const std::string& lab) { return "g_theta, " + lab; }

// ---------------- accumulators ----------------
struct Accum {
    double sx  = 0.0, sQ  = 0.0, st  = 0.0, sp  = 0.0;
    double seT = 0.0, spT = 0.0, sgT = 0.0; // theta sums in degrees
    long long n = 0;

    void add(double x, double Q2, double tabs, double phi_deg,
             double e_theta_deg, double p_theta_deg, double g_theta_deg) {
        sx  += x;
        sQ  += Q2;
        st  += tabs;
        sp  += phi_deg;
        seT += e_theta_deg;
        spT += p_theta_deg;
        sgT += g_theta_deg;
        ++n;
    }

    double mx()  const { return n ? sx  / n : std::numeric_limits<double>::quiet_NaN(); }
    double mQ()  const { return n ? sQ  / n : std::numeric_limits<double>::quiet_NaN(); }
    double mt()  const { return n ? st  / n : std::numeric_limits<double>::quiet_NaN(); }
    double mp()  const { return n ? sp  / n : std::numeric_limits<double>::quiet_NaN(); }
    double meT() const { return n ? seT / n : std::numeric_limits<double>::quiet_NaN(); }
    double mpT() const { return n ? spT / n : std::numeric_limits<double>::quiet_NaN(); }
    double mgT() const { return n ? sgT / n : std::numeric_limits<double>::quiet_NaN(); }
};

struct RowWin {
    double xBmin = 0.0, xBmax = 0.0;
    double Q2min = 0.0, Q2max = 0.0;
    double tmin = 0.0, tmax = 0.0;
    double pmin = 0.0, pmax = 0.0;
    bool valid = false;
};

struct Interval {
    double min = 0.0;
    double max = 0.0;
};

struct BinLookup {
    std::vector<RowWin> rows;
    std::vector<int> valid_rows;
    std::vector<Interval> x_bins;
    std::vector<Interval> q2_bins;
    std::vector<Interval> t_bins;
    std::map<std::tuple<int, int, int>, std::vector<int>> rows_by_kin_bin;
    bool indexed = true;
};

static int interval_index(const std::vector<Interval>& bins, double value) {
    auto it = std::upper_bound(
        bins.begin(), bins.end(), value,
        [](double v, const Interval& bin) { return v < bin.min; });
    if (it == bins.begin()) return -1;
    --it;
    return (value >= it->min && value < it->max)
        ? static_cast<int>(it - bins.begin()) : -1;
}

// ---------- ROOT binding serialized to avoid Cling races ----------
static std::mutex g_root_bind_mutex;

// ---------------- branch binder ----------------
struct BranchBinder {
    // run info
    int runnum = 0;              bool has_runnum = false;

    // topology info
    int detector1 = 0;           bool has_det1 = false;
    int detector2 = 0;           bool has_det2 = false;

    // global/exclusivity inputs
    double t1 = 0.0;             bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open = false; // degrees
    double pTmiss = 0.0;         bool has_pT = false;
    double Delta_phi = 0.0;       bool has_Delta_phi = false;
    double theta = 0.0;           bool has_theta = false;
    double Emiss2 = 0.0;         bool has_Em2 = false;
    double Mx2 = 0.0;            bool has_Mx2 = false;
    double Mx2_2 = 0.0;          bool has_Mx2_2 = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gg = false;

    // Inputs used by global_cuts when dvcsgen ycol cut is enabled
    double e_p = 0.0;            bool has_e_p = false;
    double e_theta = 0.0;        bool has_e_theta = false; // radians (also used for bin means)
    double e_phi = 0.0;          bool has_e_phi = false;

    double p2_p = 0.0;           bool has_p2_p = false;
    double p2_theta = 0.0;       bool has_p2_theta = false; // radians (also used for bin means: photon theta)
    double p2_phi = 0.0;         bool has_p2_phi = false;

    // proton theta/phi (branch names p1_theta/p1_phi, radians)
    double p1_theta = 0.0;       bool has_p1_theta = false;
    double p1_phi = 0.0;         bool has_p1_phi = false;

    // binning vars
    double x = 0.0;              bool has_x = false;
    double Q2 = 0.0;             bool has_Q2 = false;
    double phi2_rad = 0.0;       bool has_phi2 = false;  // radians

    void bind(TTree* t) {
        if (!t) return;

        // Serialize ROOT calls to avoid cling races.
        std::lock_guard<std::mutex> lock(g_root_bind_mutex);

        // Read ONLY what we need; avoid broken/unused branches like Delta_eta, DepV, etc.
        t->SetBranchStatus("*", 0);

        auto enable = [&](const char* n){
            if (t->GetBranch(n)) t->SetBranchStatus(n, 1);
        };

        // Exclusivity/global cuts
        enable("runnum");
        enable("detector1");
        enable("detector2");
        enable("t1");
        enable("open_angle_ep2");
        enable("pTmiss");
        enable("Emiss2");
        enable("Mx2");
        enable("Mx2_2");
        enable("Delta_phi");
        enable("theta");
        enable("theta_gamma_gamma");

        // Used by dvcsgen ycol cut (and now also by bin means theta averages)
        enable("e_p");
        enable("e_theta");
        enable("e_phi");
        enable("p1_phi");
        enable("p2_p");
        enable("p2_theta");
        enable("p2_phi");

        // proton theta for bin means
        enable("p1_theta");

        // Binning variables
        enable("x");
        enable("Q2");
        enable("phi2"); // radians

        // Keep ROOT's tree cache off here to avoid oddities across multiple threads/files.
        t->SetCacheSize(0);

        auto bindD = [&](const char* n, double* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };
        auto bindI = [&](const char* n, int* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };

        bindI("runnum",    &runnum,    has_runnum);
        bindI("detector1", &detector1, has_det1);
        bindI("detector2", &detector2, has_det2);

        bindD("t1",                &t1,                has_t1);
        bindD("open_angle_ep2",    &open_angle_ep2,    has_open);
        bindD("pTmiss",            &pTmiss,            has_pT);
        bindD("Emiss2",            &Emiss2,            has_Em2);
        bindD("Mx2",               &Mx2,               has_Mx2);
        bindD("Mx2_2",             &Mx2_2,             has_Mx2_2);
        bindD("Delta_phi",         &Delta_phi,         has_Delta_phi);
        bindD("theta",             &theta,             has_theta);
        bindD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gg);

        bindD("e_p",     &e_p,     has_e_p);
        bindD("e_theta", &e_theta, has_e_theta);
        bindD("e_phi",   &e_phi,   has_e_phi);
        bindD("p1_phi",  &p1_phi,  has_p1_phi);

        bindD("p2_p",     &p2_p,     has_p2_p);
        bindD("p2_theta", &p2_theta, has_p2_theta);
        bindD("p2_phi",   &p2_phi,   has_p2_phi);

        bindD("p1_theta", &p1_theta, has_p1_theta);

        bindD("x",    &x,    has_x);
        bindD("Q2",   &Q2,   has_Q2);
        bindD("phi2", &phi2_rad, has_phi2);
    }

    bool readyForCuts() const {
        return has_t1 && has_open && has_pT;
    }
    bool readyForAverages() const {
        return has_x && has_Q2 && has_phi2;
    }
    bool readyForThetaMeans() const {
        // Required for filling e_theta / p_theta / g_theta (in degrees) into CSV.
        // - e_theta uses branch e_theta
        // - p_theta uses branch p1_theta
        // - g_theta uses branch p2_theta
        return has_e_theta && has_p1_theta && has_p2_theta;
    }

    double phi_deg() const {
        if (!has_phi2) return std::numeric_limits<double>::quiet_NaN();
        double deg = std::fmod(phi2_rad * RAD2DEG, 360.0);
        if (deg < 0.0) deg += 360.0;
        if (deg >= 360.0) deg = std::nextafter(360.0, 0.0);
        return deg;
    }

    double e_theta_deg() const {
        if (!has_e_theta) return std::numeric_limits<double>::quiet_NaN();
        return e_theta * RAD2DEG;
    }
    double p_theta_deg() const {
        if (!has_p1_theta) return std::numeric_limits<double>::quiet_NaN();
        return p1_theta * RAD2DEG;
    }
    double g_theta_deg() const {
        if (!has_p2_theta) return std::numeric_limits<double>::quiet_NaN();
        return p2_theta * RAD2DEG;
    }

    int topology_index() const {
        if (!has_det1 || !has_det2) return -1;
        if (detector1 == 1 && detector2 == 1) return (int)Topology::FD_FD;
        if (detector1 == 2 && detector2 == 1) return (int)Topology::CD_FD;
        if (detector1 == 2 && detector2 == 0) return (int)Topology::CD_FT;
        return -1;
    }
};

// ---------------- cuts ----------------
// ---------------- cuts ----------------
static inline bool passes_global(const BranchBinder& b, const std::string& period_label) {
    if (!b.readyForCuts()) return false;

    if (b.has_runnum && is_excluded_run(b.runnum)) return false;

    const GlobalCutConfig& cfg = default_global_cuts();

    if (cfg.enable_topology_filter || global_cuts_require_sector_phi(cfg) || cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_det1 && b.has_det2)) {
            std::cerr << "[bin_means] FATAL: topology/sector/global-ycol selection requires detector1/detector2.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    if (global_cuts_require_sector_phi(cfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            std::cerr << "[bin_means] FATAL: sector selection requires e_phi, p1_phi, and p2_phi branches.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    if (global_cuts_require_auxiliary_kinematics(cfg)) {
        if (!(b.has_e_theta && b.has_e_phi &&
              b.has_p1_theta && b.has_p1_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            std::cerr << "[bin_means] FATAL: auxiliary fiducial cuts require e_theta, e_phi, p1_theta, p1_phi, p2_p, p2_theta, p2_phi branches.\n";
            std::exit(EXIT_FAILURE);
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p1_theta, b.p1_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            std::cerr << "[bin_means] FATAL: dvcsgen ycol cut requires e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi.\n";
            std::exit(EXIT_FAILURE);
        }

        if (global_cuts_require_sector_phi(cfg)) {
            return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                      b.detector1, b.detector2,
                                      period_label,
                                      b.e_p, b.e_theta, b.e_phi, b.p1_theta, b.p1_phi,
                                      b.p2_p, b.p2_theta, b.p2_phi,
                                      cfg);
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    if (global_cuts_require_sector_phi(cfg)) {
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_phi, b.p1_phi, b.p2_phi,
                                  cfg);
    }

    if (cfg.enable_topology_filter) {
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  cfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss, cfg);
} //endfor

struct ProductionCut {
    double low = -std::numeric_limits<double>::infinity();
    double high = std::numeric_limits<double>::infinity();
};
using CutMap = std::unordered_map<std::string, ProductionCut>;
static std::unordered_map<std::string, CutMap> g_cut_cache;
static std::once_flag g_cut_once;

static const std::set<std::string>& supported_production_variables() {
    static const std::set<std::string> vars = {
        "Delta_phi", "theta", "theta_gamma_gamma", "pTmiss",
        "Emiss2", "Mx2", "Mx2_2"
    };
    return vars;
}

static ProductionCut parse_production_cut(const nlohmann::json& obj,
                                          const std::string& key,
                                          const std::string& variable) {
    if (!obj.is_object()) {
        std::cerr << "[bin_means] FATAL: malformed cut object "
                  << key << "/" << variable << std::endl;
        std::exit(EXIT_FAILURE);
    }
    ProductionCut c;
    if (obj.contains("cut_low") && !obj["cut_low"].is_null())
        c.low = obj["cut_low"].get<double>();
    if (obj.contains("cut_high") && !obj["cut_high"].is_null())
        c.high = obj["cut_high"].get<double>();
    if (!std::isfinite(c.low) && !std::isfinite(c.high) &&
        obj.contains("mean") && obj.contains("std")) {
        const double mean=obj["mean"].get<double>();
        const double sigma=obj["std"].get<double>();
        if (std::isfinite(mean) && std::isfinite(sigma) && sigma>0.0) {
            c.low=mean-3.0*sigma;
            c.high=mean+3.0*sigma;
        }
    }
    if (!(std::isfinite(c.low)||std::isfinite(c.high)) || c.low>c.high) {
        std::cerr << "[bin_means] FATAL: invalid cut "
                  << key << "/" << variable << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return c;
}

static void load_production_cuts_once() {
    std::ifstream fin(kCutsJSON);
    if (!fin.is_open()) {
        std::cerr << "[bin_means] FATAL: cannot open " << kCutsJSON << std::endl;
        std::exit(EXIT_FAILURE);
    }
    nlohmann::json j;
    try { fin >> j; }
    catch (const std::exception& e) {
        std::cerr << "[bin_means] FATAL: failed parsing " << kCutsJSON
                  << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (!j.is_object()) {
        std::cerr << "[bin_means] FATAL: cuts JSON is not an object." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    const std::vector<std::string> periods={
        "DVCS_Fa18_Inb","DVCS_Fa18_Out","DVCS_Sp19_Inb",
        "DVCS_Sp18_Inb","DVCS_Sp18_Out"
    };
    const std::vector<std::string> topologies={"FD_FD","CD_FD","CD_FT"};
    for (const auto& period:periods) {
        for (const auto& topology:topologies) {
            const std::string key=period+"_"+topology;
            if (!j.contains(key) || !j[key].is_object() ||
                !j[key].contains("data") || !j[key]["data"].is_object()) {
                std::cerr << "[bin_means] FATAL: missing data cut block "
                          << key << std::endl;
                std::exit(EXIT_FAILURE);
            }
            CutMap cuts;
            for (auto it=j[key]["data"].begin();it!=j[key]["data"].end();++it) {
                if (supported_production_variables().count(it.key())==0) {
                    std::cerr << "[bin_means] FATAL: unsupported production variable "
                              << it.key() << " in " << key << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                cuts.emplace(it.key(),parse_production_cut(it.value(),key,it.key()));
            }
            if (cuts.find("Mx2")==cuts.end()) {
                std::cerr << "[bin_means] FATAL: missing mandatory Mx2 in "
                          << key << std::endl;
                std::exit(EXIT_FAILURE);
            }
            g_cut_cache.emplace(key,std::move(cuts));
        }
    }
}

static const CutMap& resolve_production_cuts(const std::string& period_json_tag,
                                             Topology topo) {
    std::call_once(g_cut_once,load_production_cuts_once);
    const std::string key=period_json_tag+"_"+topo_tag(topo);
    auto it=g_cut_cache.find(key);
    if (it==g_cut_cache.end()) {
        std::cerr << "[bin_means] FATAL: missing cached cuts " << key << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return it->second;
}

static bool production_value(const BranchBinder& b,const std::string& v,
                             bool& has,double& value) {
    if (v=="Delta_phi") {has=b.has_Delta_phi; value=b.Delta_phi;}
    else if (v=="theta") {has=b.has_theta; value=b.theta;}
    else if (v=="theta_gamma_gamma") {has=b.has_theta_gg; value=b.theta_gamma_gamma;}
    else if (v=="pTmiss") {has=b.has_pT; value=b.pTmiss;}
    else if (v=="Emiss2") {has=b.has_Em2; value=b.Emiss2;}
    else if (v=="Mx2") {has=b.has_Mx2; value=b.Mx2;}
    else if (v=="Mx2_2") {has=b.has_Mx2_2; value=b.Mx2_2;}
    else return false;
    return true;
}

static bool passes_production_cuts(const CutMap& cuts,const BranchBinder& b) {
    for (const auto& kv:cuts) {
        bool has=false; double value=0.0;
        if (!production_value(b,kv.first,has,value)) {
            std::cerr << "[bin_means] FATAL: unsupported cached variable "
                      << kv.first << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (!has) {
            std::cerr << "[bin_means] FATAL: missing ROOT branch required by cut "
                      << kv.first << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (!std::isfinite(value) || value<kv.second.low || value>kv.second.high)
            return false;
    }
    return true;
}

static inline bool in_range(double v, double a, double b) { return (v >= a) && (v < b); }
static bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (pmax_deg > pmin_deg) return in_range(phi_deg, pmin_deg, pmax_deg);
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg); // wrap-around
}
static bool row_accepts_kin(const BranchBinder& b,
                            double xBmin, double xBmax,
                            double Q2min, double Q2max,
                            double tmin, double tmax) {
    const double tabs = std::fabs(b.t1);
    if (!in_range(b.x,  xBmin, xBmax)) return false;
    if (!in_range(b.Q2, Q2min, Q2max)) return false;
    if (!in_range(tabs, tmin, tmax))   return false;
    return true;
}

static BinLookup build_bin_lookup(const CSV& csv) {
    const int c_xBmin = col(csv, "xBmin"), c_xBmax = col(csv, "xBmax");
    const int c_Q2min = col(csv, "Q2min"), c_Q2max = col(csv, "Q2max");
    const int c_tmin = col(csv, "t_abs_min"), c_tmax = col(csv, "t_abs_max");
    const int c_pmin = col(csv, "phimin"), c_pmax = col(csv, "phimax");
    const int c_valid = col(csv, "valid bin");
    auto toD = [](const std::string& text) {
        if (text.empty()) return std::numeric_limits<double>::quiet_NaN();
        char* end = nullptr;
        const double value = std::strtod(text.c_str(), &end);
        return end == text.c_str() ? std::numeric_limits<double>::quiet_NaN() : value;
    };
    BinLookup lookup;
    lookup.rows.reserve(csv.rows.size());
    for (const auto& row : csv.rows) {
        RowWin w;
        w.xBmin=toD(row[c_xBmin]); w.xBmax=toD(row[c_xBmax]);
        w.Q2min=toD(row[c_Q2min]); w.Q2max=toD(row[c_Q2max]);
        w.tmin=toD(row[c_tmin]); w.tmax=toD(row[c_tmax]);
        w.pmin=toD(row[c_pmin]); w.pmax=toD(row[c_pmax]);
        const std::string& v=row[c_valid];
        w.valid=(v=="1" || v=="1.0" || v=="true" || v=="TRUE");
        lookup.rows.push_back(w);
        if (w.valid) {
            lookup.valid_rows.push_back(static_cast<int>(lookup.rows.size()) - 1);
            lookup.x_bins.push_back({w.xBmin,w.xBmax});
            lookup.q2_bins.push_back({w.Q2min,w.Q2max});
            lookup.t_bins.push_back({w.tmin,w.tmax});
        }
    }
    auto normalize=[](std::vector<Interval>& bins) {
        std::sort(bins.begin(),bins.end(),[](const Interval&a,const Interval&b){
            return a.min<b.min || (a.min==b.min && a.max<b.max);
        });
        bins.erase(std::unique(bins.begin(),bins.end(),[](const Interval&a,const Interval&b){
            return a.min==b.min && a.max==b.max;
        }),bins.end());
    };
    normalize(lookup.x_bins); normalize(lookup.q2_bins); normalize(lookup.t_bins);

    auto non_overlapping = [](const std::vector<Interval>& bins) {
        for (std::size_t i = 1; i < bins.size(); ++i) {
            if (bins[i].min < bins[i - 1].max) return false;
        }
        return true;
    };
    lookup.indexed = non_overlapping(lookup.x_bins) &&
                     non_overlapping(lookup.q2_bins) &&
                     non_overlapping(lookup.t_bins);
    if (!lookup.indexed) {
        print_banner("CSV kinematic intervals overlap; using exact linear row scan fallback");
        return lookup;
    }

    for (int r=0;r<static_cast<int>(lookup.rows.size());++r) {
        const RowWin& w=lookup.rows[static_cast<std::size_t>(r)];
        if (!w.valid) continue;
        const int ix=interval_index(lookup.x_bins,w.xBmin);
        const int iq=interval_index(lookup.q2_bins,w.Q2min);
        const int it=interval_index(lookup.t_bins,w.tmin);
        if (ix<0 || iq<0 || it<0) {
            std::cerr << "[bin_means] FATAL: failed to index a valid CSV bin row.\n";
            std::exit(EXIT_FAILURE);
        }
        lookup.rows_by_kin_bin[std::make_tuple(ix,iq,it)].push_back(r);
    }
    return lookup;
}

// ---------------- per period work ----------------
struct PeriodResult {
    std::unordered_map<int, Accum> per_row; // row index -> Accum
};

static PeriodResult process_period(const std::string& period_key, TTree* tree, const BinLookup& lookup) {
    PeriodResult R;
    if (!tree) return R;

    PeriodKeyMap tags = period_key_to_tags(period_key);
    BranchBinder b; b.bind(tree);

    if (!b.readyForCuts() || !b.readyForAverages()) {
        std::cerr << "[bin_means] FATAL: Tree for '" << period_key
                  << "' missing branches (t1/open_angle_ep2/pTmiss/x/Q2/phi2)." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (!b.has_det1 || !b.has_det2) {
        std::cerr << "[bin_means] FATAL: Tree for '" << period_key
                  << "' missing detector1/detector2 for topology determination." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Required for theta means (in degrees) written to CSV
    if (!b.readyForThetaMeans()) {
        std::cerr << "[bin_means] FATAL: Tree for '" << period_key
                  << "' missing theta branches required for theta means. Required: "
                  << "e_theta (electron), p1_theta (proton), p2_theta (photon)."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    const Long64_t N = tree->GetEntries();
    print_banner("Processing period " + period_key + " with " + std::to_string((long long)N) + " entries");

    const bool dbg = (std::getenv("BINMEANS_DEBUG") != nullptr);

    const std::array<const CutMap*, 3> production_cuts = {{
        &resolve_production_cuts(tags.json_tag, Topology::FD_FD),
        &resolve_production_cuts(tags.json_tag, Topology::CD_FD),
        &resolve_production_cuts(tags.json_tag, Topology::CD_FT)
    }};

    long long n_pass_global = 0;
    long long n_pass_excl   = 0;
    long long n_used_rows   = 0;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (dbg && i < 3) {
            std::cout << "[bin_means][" << period_key << "] sample "
                      << i << " t1=" << b.t1
                      << " open_angle_ep2=" << b.open_angle_ep2
                      << " pTmiss=" << b.pTmiss
                      << " x=" << b.x << " Q2=" << b.Q2
                      << " phi2(rad)=" << b.phi2_rad
                      << " theta_gamma_gamma=" << b.theta_gamma_gamma
                      << " e_theta(rad)=" << b.e_theta
                      << " p1_theta(rad)=" << b.p1_theta
                      << " p2_theta(rad)=" << b.p2_theta
                      << " det1=" << (b.has_det1 ? b.detector1 : -1)
                      << " det2=" << (b.has_det2 ? b.detector2 : -1)
                      << " runnum=" << (b.has_runnum ? b.runnum : -1)
                      << std::endl;
        }

        if (!passes_global(b, tags.period_label)) continue;
        ++n_pass_global;

        int topo_idx = b.topology_index();
        if (topo_idx < 0 || topo_idx > 2) continue;
        if (!passes_production_cuts(*production_cuts[static_cast<std::size_t>(topo_idx)], b)) continue;
        ++n_pass_excl;

        const double phi_deg = b.phi_deg();

        const double e_theta_deg = b.e_theta_deg();
        const double p_theta_deg = b.p_theta_deg();
        const double g_theta_deg = b.g_theta_deg();

        const double tabs = std::fabs(b.t1);
        bool used_any_row = false;

        auto test_rows = [&](const std::vector<int>& candidate_rows) {
            for (int r : candidate_rows) {
                const RowWin& w = lookup.rows[static_cast<std::size_t>(r)];
                if (!row_accepts_kin(b, w.xBmin, w.xBmax,
                                     w.Q2min, w.Q2max,
                                     w.tmin, w.tmax)) continue;
                if (!row_accepts_phi(phi_deg, w.pmin, w.pmax)) continue;
                R.per_row[r].add(b.x, b.Q2, tabs, phi_deg,
                                 e_theta_deg, p_theta_deg, g_theta_deg);
                used_any_row = true;
            }
        };

        if (lookup.indexed) {
            const int ix = interval_index(lookup.x_bins, b.x);
            const int iq = interval_index(lookup.q2_bins, b.Q2);
            const int it = interval_index(lookup.t_bins, tabs);
            if (ix < 0 || iq < 0 || it < 0) continue;

            const auto candidates = lookup.rows_by_kin_bin.find(std::make_tuple(ix, iq, it));
            if (candidates == lookup.rows_by_kin_bin.end()) continue;
            test_rows(candidates->second);
        } else {
            test_rows(lookup.valid_rows);
        }

        if (used_any_row) ++n_used_rows;
    }

    print_banner("Finished " + period_key +
                 "  global_pass=" + std::to_string(n_pass_global) +
                 "  sig_pass="    + std::to_string(n_pass_excl) +
                 "  used="        + std::to_string(n_used_rows));

    return R;
}

// ---------------- combine groups and write ----------------
static inline std::string fmt8(double v) {
    if (!std::isfinite(v)) return std::string();
    std::ostringstream oss; oss << std::fixed << std::setprecision(8) << v;
    return oss.str();
}

static void fill_combined_groups(CSV& csv,
                                 const std::unordered_map<std::string, std::unordered_map<int, Accum>>& per_period) {
    auto get = [&](const std::string& lab)->const std::unordered_map<int, Accum>* {
        auto it = per_period.find(lab);
        return (it == per_period.end()) ? nullptr : &it->second;
    };

    const auto* Fa18Inb = get("Fa18 Inb");
    const auto* Fa18Out = get("Fa18 Out");
    const auto* Sp18Inb = get("Sp18 Inb");
    const auto* Sp18Out = get("Sp18 Out");
    const auto* Sp19Inb = get("Sp19 Inb"); (void)Sp19Inb;

    const int c_xB_Fa18   = col(csv, col_xBavg("Fa18"));
    const int c_Q2_Fa18   = col(csv, col_Q2avg("Fa18"));
    const int c_t_Fa18    = col(csv, col_tavg ("Fa18"));
    const int c_phi_Fa18  = col(csv, col_phiavg("Fa18"));

    const int c_xB_Sp18   = col(csv, col_xBavg("Sp18"));
    const int c_Q2_Sp18   = col(csv, col_Q2avg("Sp18"));
    const int c_t_Sp18    = col(csv, col_tavg ("Sp18"));
    const int c_phi_Sp18  = col(csv, col_phiavg("Sp18"));

    const int c_xB_106    = col(csv, col_xBavg("10.6 GeV"));
    const int c_Q2_106    = col(csv, col_Q2avg("10.6 GeV"));
    const int c_t_106     = col(csv, col_tavg ("10.6 GeV"));
    const int c_phi_106   = col(csv, col_phiavg("10.6 GeV"));

    // theta group columns (all must exist per schema; fail fast if missing)
    const int c_e_Fa18   = col(csv, col_e_theta("Fa18"));
    const int c_p_Fa18   = col(csv, col_p_theta("Fa18"));
    const int c_g_Fa18   = col(csv, col_g_theta("Fa18"));

    const int c_e_Sp18   = col(csv, col_e_theta("Sp18"));
    const int c_p_Sp18   = col(csv, col_p_theta("Sp18"));
    const int c_g_Sp18   = col(csv, col_g_theta("Sp18"));

    const int c_e_106    = col(csv, col_e_theta("10.6 GeV"));
    const int c_p_106    = col(csv, col_p_theta("10.6 GeV"));
    const int c_g_106    = col(csv, col_g_theta("10.6 GeV"));

    auto combine = [&](const std::vector<const std::unordered_map<int, Accum>*>& parts, int row)->Accum {
        Accum a;
        for (auto* p : parts) {
            if (!p) continue;
            auto it = p->find(row);
            if (it == p->end()) continue;
            const Accum& r = it->second;
            a.sx  += r.sx;
            a.sQ  += r.sQ;
            a.st  += r.st;
            a.sp  += r.sp;
            a.seT += r.seT;
            a.spT += r.spT;
            a.sgT += r.sgT;
            a.n   += r.n;
        }
        return a;
    };

    print_banner("Combining period means into Fa18, Sp18, 10.6 GeV groups (skip if n < 10)");
    long long wrote_Fa18 = 0, wrote_Sp18 = 0, wrote_106 = 0;
    long long skip_Fa18 = 0, skip_Sp18 = 0, skip_106 = 0;

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        {
            Accum a = combine({Fa18Inb, Fa18Out}, r);
            if (a.n >= MIN_BIN_COUNT) {
                csv.rows[r][c_xB_Fa18]  = fmt8(a.mx());
                csv.rows[r][c_Q2_Fa18]  = fmt8(a.mQ());
                csv.rows[r][c_t_Fa18]   = fmt8(a.mt());
                csv.rows[r][c_phi_Fa18] = fmt8(a.mp());

                csv.rows[r][c_e_Fa18]   = fmt8(a.meT());
                csv.rows[r][c_p_Fa18]   = fmt8(a.mpT());
                csv.rows[r][c_g_Fa18]   = fmt8(a.mgT());

                ++wrote_Fa18;
            } else {
                csv.rows[r][c_xB_Fa18].clear();
                csv.rows[r][c_Q2_Fa18].clear();
                csv.rows[r][c_t_Fa18].clear();
                csv.rows[r][c_phi_Fa18].clear();

                csv.rows[r][c_e_Fa18].clear();
                csv.rows[r][c_p_Fa18].clear();
                csv.rows[r][c_g_Fa18].clear();

                ++skip_Fa18;
            }
        }
        {
            Accum a = combine({Sp18Inb, Sp18Out}, r);
            if (a.n >= MIN_BIN_COUNT) {
                csv.rows[r][c_xB_Sp18]  = fmt8(a.mx());
                csv.rows[r][c_Q2_Sp18]  = fmt8(a.mQ());
                csv.rows[r][c_t_Sp18]   = fmt8(a.mt());
                csv.rows[r][c_phi_Sp18] = fmt8(a.mp());

                csv.rows[r][c_e_Sp18]   = fmt8(a.meT());
                csv.rows[r][c_p_Sp18]   = fmt8(a.mpT());
                csv.rows[r][c_g_Sp18]   = fmt8(a.mgT());

                ++wrote_Sp18;
            } else {
                csv.rows[r][c_xB_Sp18].clear();
                csv.rows[r][c_Q2_Sp18].clear();
                csv.rows[r][c_t_Sp18].clear();
                csv.rows[r][c_phi_Sp18].clear();

                csv.rows[r][c_e_Sp18].clear();
                csv.rows[r][c_p_Sp18].clear();
                csv.rows[r][c_g_Sp18].clear();

                ++skip_Sp18;
            }
        }
        {
            Accum a = combine({Fa18Inb, Fa18Out, Sp18Inb, Sp18Out}, r);
            if (a.n >= MIN_BIN_COUNT) {
                csv.rows[r][c_xB_106]  = fmt8(a.mx());
                csv.rows[r][c_Q2_106]  = fmt8(a.mQ());
                csv.rows[r][c_t_106]   = fmt8(a.mt());
                csv.rows[r][c_phi_106] = fmt8(a.mp());

                csv.rows[r][c_e_106]   = fmt8(a.meT());
                csv.rows[r][c_p_106]   = fmt8(a.mpT());
                csv.rows[r][c_g_106]   = fmt8(a.mgT());

                ++wrote_106;
            } else {
                csv.rows[r][c_xB_106].clear();
                csv.rows[r][c_Q2_106].clear();
                csv.rows[r][c_t_106].clear();
                csv.rows[r][c_phi_106].clear();

                csv.rows[r][c_e_106].clear();
                csv.rows[r][c_p_106].clear();
                csv.rows[r][c_g_106].clear();

                ++skip_106;
            }
        }
    }

    print_banner("Group write summary: Fa18 wrote=" + std::to_string(wrote_Fa18) +
                 " skipped=" + std::to_string(skip_Fa18) +
                 " | Sp18 wrote=" + std::to_string(wrote_Sp18) +
                 " skipped=" + std::to_string(skip_Sp18) +
                 " | 10.6 wrote=" + std::to_string(wrote_106) +
                 " skipped=" + std::to_string(skip_106));
}


// ---------------- analysis-note diagnostics ----------------
static inline bool finite_cell(const std::string& s, double& value) {
    if (s.empty()) return false;
    char* end = nullptr;
    value = std::strtod(s.c_str(), &end);
    return end != s.c_str() && std::isfinite(value);
}

static inline double interval_center_for_phi(double lo, double hi) {
    if (hi > lo) return 0.5 * (lo + hi);
    // Wrapped interval, e.g. 330--30 degrees.
    double width = (360.0 - lo) + hi;
    double c = lo + 0.5 * width;
    if (c >= 360.0) c -= 360.0;
    return c;
}

static void write_analysis_note_bin_means_csv(
    const CSV& csv,
    const BinLookup& lookup,
    const std::unordered_map<std::string, std::unordered_map<int, Accum>>& per_period,
    const std::string& outdir) {

    gSystem->mkdir(outdir.c_str(), true);
    const std::string path = outdir + "/kinematic_means_by_bin.csv";
    std::ofstream out(path);
    if (!out.is_open()) {
        std::cerr << "[bin_means] WARNING: cannot write analysis-note CSV: " << path << std::endl;
        return;
    }

    const int cx106 = col(csv, col_xBavg("10.6 GeV"));
    const int cq106 = col(csv, col_Q2avg("10.6 GeV"));
    const int ct106 = col(csv, col_tavg("10.6 GeV"));
    const int cp106 = col(csv, col_phiavg("10.6 GeV"));
    const int ce106 = col(csv, col_e_theta("10.6 GeV"));
    const int cpr106 = col(csv, col_p_theta("10.6 GeV"));
    const int cg106 = col(csv, col_g_theta("10.6 GeV"));

    const int cx19 = col(csv, col_xBavg("Sp19 Inb"));
    const int cq19 = col(csv, col_Q2avg("Sp19 Inb"));
    const int ct19 = col(csv, col_tavg("Sp19 Inb"));
    const int cp19 = col(csv, col_phiavg("Sp19 Inb"));

    out << "row,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax,"
        << "xBavg_10p6,Q2avg_10p6,tavg_10p6,phiavg_10p6,e_theta_10p6,p_theta_10p6,g_theta_10p6,"
        << "selected_events_10p6,xBavg_sp19,Q2avg_sp19,tavg_sp19,phiavg_sp19,selected_events_sp19\n";

    auto n_for = [&](const std::string& period, int row)->long long {
        auto ip = per_period.find(period);
        if (ip == per_period.end()) return 0;
        auto ir = ip->second.find(row);
        return (ir == ip->second.end()) ? 0 : ir->second.n;
    };

    for (int r = 0; r < static_cast<int>(lookup.rows.size()); ++r) {
        const RowWin& w = lookup.rows[static_cast<std::size_t>(r)];
        if (!w.valid) continue;
        const long long n106 = n_for("Fa18 Inb", r) + n_for("Fa18 Out", r) +
                               n_for("Sp18 Inb", r) + n_for("Sp18 Out", r);
        const long long n19 = n_for("Sp19 Inb", r);
        out << r << ',' << w.xBmin << ',' << w.xBmax << ','
            << w.Q2min << ',' << w.Q2max << ',' << w.tmin << ',' << w.tmax << ','
            << w.pmin << ',' << w.pmax << ','
            << csv.rows[r][cx106] << ',' << csv.rows[r][cq106] << ',' << csv.rows[r][ct106] << ',' << csv.rows[r][cp106] << ','
            << csv.rows[r][ce106] << ',' << csv.rows[r][cpr106] << ',' << csv.rows[r][cg106] << ',' << n106 << ','
            << csv.rows[r][cx19] << ',' << csv.rows[r][cq19] << ',' << csv.rows[r][ct19] << ',' << csv.rows[r][cp19] << ',' << n19 << '\n';
    }
    std::cout << "[bin_means] Analysis-note bin-mean table written to: " << path << std::endl;
}

static void make_analysis_note_mean_vs_center_plot(const CSV& csv,
                                                    const BinLookup& lookup,
                                                    const std::string& outdir) {
    gSystem->mkdir(outdir.c_str(), true);
    gStyle->SetOptStat(0);

    const int c_x = col(csv, col_xBavg("10.6 GeV"));
    const int c_q = col(csv, col_Q2avg("10.6 GeV"));
    const int c_t = col(csv, col_tavg("10.6 GeV"));
    const int c_p = col(csv, col_phiavg("10.6 GeV"));

    std::array<std::vector<double>,4> centers, means;
    for (int r = 0; r < static_cast<int>(lookup.rows.size()); ++r) {
        const RowWin& w = lookup.rows[static_cast<std::size_t>(r)];
        if (!w.valid) continue;
        const std::array<double,4> xcenter = {
            0.5*(w.xBmin+w.xBmax), 0.5*(w.Q2min+w.Q2max),
            0.5*(w.tmin+w.tmax), interval_center_for_phi(w.pmin,w.pmax)
        };
        const std::array<int,4> cc = {c_x,c_q,c_t,c_p};
        for (int k=0;k<4;++k) {
            double m=0.0;
            if (!finite_cell(csv.rows[r][cc[k]], m)) continue;
            centers[k].push_back(xcenter[k]);
            means[k].push_back(m);
        }
    }

    TCanvas c("c_note_bin_means", "", 1500, 1200);
    c.Divide(2,2,0.012,0.012);
    const std::array<const char*,4> xt = {"bin center x_{B}", "bin center Q^{2} (GeV^{2})", "bin center -t (GeV^{2})", "bin center #phi (deg)"};
    const std::array<const char*,4> yt = {"#LT x_{B} #GT", "#LT Q^{2} #GT (GeV^{2})", "#LT -t #GT (GeV^{2})", "#LT #phi #GT (deg)"};

    for (int k=0;k<4;++k) {
        c.cd(k+1);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.08);
        if (centers[k].empty()) continue;
        double xmin=*std::min_element(centers[k].begin(),centers[k].end());
        double xmax=*std::max_element(centers[k].begin(),centers[k].end());
        double ymin=*std::min_element(means[k].begin(),means[k].end());
        double ymax=*std::max_element(means[k].begin(),means[k].end());
        double lo=std::min(xmin,ymin), hi=std::max(xmax,ymax), pad=0.06*(hi-lo>0?hi-lo:1.0);
        TH1D frame((std::string("frame_means_")+std::to_string(k)).c_str(), "", 10, lo-pad, hi+pad);
        frame.SetMinimum(lo-pad); frame.SetMaximum(hi+pad);
        frame.GetXaxis()->SetTitle(xt[k]); frame.GetYaxis()->SetTitle(yt[k]);
        frame.GetXaxis()->SetTitleSize(0.050); frame.GetYaxis()->SetTitleSize(0.050);
        frame.GetXaxis()->SetLabelSize(0.043); frame.GetYaxis()->SetLabelSize(0.043);
        frame.DrawCopy("AXIS");
        TLine diag(lo-pad,lo-pad,hi+pad,hi+pad); diag.SetLineStyle(2); diag.SetLineWidth(2); diag.DrawClone("SAME");
        TGraph gr(static_cast<int>(centers[k].size()), centers[k].data(), means[k].data());
        gr.SetMarkerStyle(20); gr.SetMarkerSize(0.55); gr.DrawClone("P SAME");
    }
    c.cd(0);
    const std::string path=outdir+"/kinematic_means_vs_bin_centers.png";
    c.SaveAs(path.c_str());
    std::cout << "[bin_means] Analysis-note mean-vs-center plot written to: " << path << std::endl;
}

static void make_analysis_note_kinematic_coverage_plot(const CSV& csv,
                                                       const BinLookup& lookup,
                                                       const std::string& outdir) {
    gSystem->mkdir(outdir.c_str(), true);
    const int cx=col(csv,col_xBavg("10.6 GeV"));
    const int cq=col(csv,col_Q2avg("10.6 GeV"));
    const int ct=col(csv,col_tavg("10.6 GeV"));
    std::vector<double> x,q,t;
    std::set<std::tuple<long long,long long,long long>> seen;
    for (int r=0;r<static_cast<int>(lookup.rows.size());++r) {
        const RowWin& w=lookup.rows[static_cast<std::size_t>(r)]; if(!w.valid) continue;
        double xv,qv,tv; if(!finite_cell(csv.rows[r][cx],xv)||!finite_cell(csv.rows[r][cq],qv)||!finite_cell(csv.rows[r][ct],tv)) continue;
        auto key=std::make_tuple((long long)std::llround(xv*1e5),(long long)std::llround(qv*1e5),(long long)std::llround(tv*1e5));
        if(!seen.insert(key).second) continue;
        x.push_back(xv); q.push_back(qv); t.push_back(tv);
    }
    TCanvas c("c_note_coverage","",1500,650); c.Divide(2,1,0.012,0.012);
    c.cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.08);
    if(!x.empty()) { TGraph g((int)x.size(),x.data(),q.data()); g.SetTitle(""); g.SetMarkerStyle(20); g.SetMarkerSize(0.65); g.GetXaxis()->SetTitle("#LT x_{B} #GT"); g.GetYaxis()->SetTitle("#LT Q^{2} #GT (GeV^{2})"); g.DrawClone("AP"); }
    c.cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.08);
    if(!x.empty()) { TGraph g((int)x.size(),x.data(),t.data()); g.SetTitle(""); g.SetMarkerStyle(20); g.SetMarkerSize(0.65); g.GetXaxis()->SetTitle("#LT x_{B} #GT"); g.GetYaxis()->SetTitle("#LT -t #GT (GeV^{2})"); g.DrawClone("AP"); }
    const std::string path=outdir+"/kinematic_mean_coverage.png"; c.SaveAs(path.c_str());
    std::cout << "[bin_means] Analysis-note kinematic-coverage plot written to: " << path << std::endl;
}

bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers,
                          const BinMeansOptions& options) {
    // Make ROOT set up its global mutexes.
    ROOT::EnableThreadSafety();

    CSV csv;
    if (!load_csv(csv_path, csv)) return false;

    const BinLookup bin_lookup = build_bin_lookup(csv);

    // Existing mean columns
    std::unordered_map<std::string,int> cxB, cQ2, ct, cphi;
    for (const auto& lab : csv_period_labels()) {
        cxB[lab]  = col(csv, col_xBavg(lab));
        cQ2[lab]  = col(csv, col_Q2avg(lab));
        ct[lab]   = col(csv, col_tavg (lab));
        cphi[lab] = col(csv, col_phiavg(lab));
    }

    // Theta mean columns
    std::unordered_map<std::string,int> ceT, cpT, cgT;

    for (const auto& lab : e_theta_labels()) {
        ceT[lab] = col(csv, col_e_theta(lab));
    }
    for (const auto& lab : p_theta_labels()) {
        cpT[lab] = col(csv, col_p_theta(lab));
    }
    for (const auto& lab : g_theta_labels()) {
        cgT[lab] = col(csv, col_g_theta(lab));
    }

    std::vector<std::string> period_keys;
    for (const auto& P : CANONICAL_PERIODS()) {
        auto it = dataTrees.find(P.tree_key);
        if (it != dataTrees.end() && it->second) period_keys.push_back(P.tree_key);
    }
    if (period_keys.empty()) {
        std::cerr << "[bin_means] FATAL: no DVCS trees available in dataTrees." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    print_banner("Will process periods:");
    for (const auto& k : period_keys) std::cout << "  - " << k << std::endl;

    std::vector<PeriodResult> results(period_keys.size());

    const int nth = std::max(1, std::min(7, max_workers));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
    for (int i = 0; i < (int)period_keys.size(); ++i) {
        const std::string& pk = period_keys[i];
        auto it = dataTrees.find(pk);
        TTree* t = (it == dataTrees.end()) ? nullptr : it->second;
        results[i] = process_period(pk, t, bin_lookup);
    }

    std::unordered_map<std::string, std::unordered_map<int, Accum>> per_period_rows;

    for (int i = 0; i < (int)period_keys.size(); ++i) {
        const std::string& pk = period_keys[i];
        PeriodKeyMap tags = period_key_to_tags(pk);
        if (tags.csv_label.empty()) continue;

        const auto& perrow = results[i].per_row;
        per_period_rows[tags.csv_label] = perrow;

        long long wrote = 0, skipped = 0;
        for (const auto& kv : perrow) {
            int r = kv.first;
            const Accum& a = kv.second;
            if (a.n >= MIN_BIN_COUNT) {
                // Existing means
                csv.rows[r][cxB[tags.csv_label]]  = fmt8(a.mx());
                csv.rows[r][cQ2[tags.csv_label]]  = fmt8(a.mQ());
                csv.rows[r][ct [tags.csv_label]]  = fmt8(a.mt());
                csv.rows[r][cphi[tags.csv_label]] = fmt8(a.mp());

                // Theta means (degrees): all must exist now
                auto itE = ceT.find(tags.csv_label);
                if (itE == ceT.end()) {
                    std::cerr << "[bin_means] FATAL: expected column missing for e_theta label: "
                              << tags.csv_label << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                auto itP = cpT.find(tags.csv_label);
                if (itP == cpT.end()) {
                    std::cerr << "[bin_means] FATAL: expected column missing for p_theta label: "
                              << tags.csv_label << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                auto itG = cgT.find(tags.csv_label);
                if (itG == cgT.end()) {
                    std::cerr << "[bin_means] FATAL: expected column missing for g_theta label: "
                              << tags.csv_label << std::endl;
                    std::exit(EXIT_FAILURE);
                }

                csv.rows[r][itE->second] = fmt8(a.meT());
                csv.rows[r][itP->second] = fmt8(a.mpT());
                csv.rows[r][itG->second] = fmt8(a.mgT());

                ++wrote;
            } else {
                // Existing means
                csv.rows[r][cxB[tags.csv_label]].clear();
                csv.rows[r][cQ2[tags.csv_label]].clear();
                csv.rows[r][ct [tags.csv_label]].clear();
                csv.rows[r][cphi[tags.csv_label]].clear();

                // Theta means
                auto itE = ceT.find(tags.csv_label);
                if (itE == ceT.end()) {
                    std::cerr << "[bin_means] FATAL: expected column missing for e_theta label: "
                              << tags.csv_label << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                auto itP = cpT.find(tags.csv_label);
                if (itP == cpT.end()) {
                    std::cerr << "[bin_means] FATAL: expected column missing for p_theta label: "
                              << tags.csv_label << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                auto itG = cgT.find(tags.csv_label);
                if (itG == cgT.end()) {
                    std::cerr << "[bin_means] FATAL: expected column missing for g_theta label: "
                              << tags.csv_label << std::endl;
                    std::exit(EXIT_FAILURE);
                }

                csv.rows[r][itE->second].clear();
                csv.rows[r][itP->second].clear();
                csv.rows[r][itG->second].clear();

                ++skipped;
            }
        }
        print_banner("Wrote " + std::to_string(wrote) +
                     " row means (n >= " + std::to_string(MIN_BIN_COUNT) + ") and skipped " +
                     std::to_string(skipped) + " for " + tags.csv_label);
    }

    fill_combined_groups(csv, per_period_rows);

    if (options.make_note_outputs) {
        write_analysis_note_bin_means_csv(csv, bin_lookup, per_period_rows, options.note_output_dir);
        make_analysis_note_mean_vs_center_plot(csv, bin_lookup, options.note_output_dir);
        make_analysis_note_kinematic_coverage_plot(csv, bin_lookup, options.note_output_dir);
    }

    if (!write_csv(csv_path, csv)) return false;
    print_banner(std::string("Updated bin means in: ") + csv_path);
    return true;
}