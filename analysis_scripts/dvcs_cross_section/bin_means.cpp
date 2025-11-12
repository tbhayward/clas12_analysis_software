// bin_means.cpp
// Progress-enabled version. Single pass through each TTree per period.
// Uses ONLY phi2 (radians) -> converts to degrees for CSV phimin/phimax checks.

#include "bin_means.h"
#include "periods.h"
#include "load_binning_scheme.h"

#include <TTree.h>
#include <TStopwatch.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
                                           double elapsed_s)
{
    double rate = (elapsed_s > 0.0) ? (double)i / elapsed_s : 0.0;
    std::cout << "[bin_means][" << tag << "] "
              << std::fixed << std::setprecision(1)
              << pct << "%  "
              << i << "/" << N
              << "  global=" << pass_global
              << "  sig=" << pass_3sig
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
    std::string json_tag;   // e.g. "DVCS_Fa18_Inb"
    std::string csv_label;  // e.g. "Fa18 Inb"
};

static inline PeriodKeyMap period_key_to_tags(const std::string& k) {
    std::string s; s.reserve(k.size());
    for (char c : k) s.push_back((char)std::tolower((unsigned char)c));
    auto has = [&](const char* sub){ return s.find(sub) != std::string::npos; };

    if (has("fa18") && has("inb")) return {"DVCS_Fa18_Inb", "Fa18 Inb"};
    if (has("fa18") && has("out")) return {"DVCS_Fa18_Out", "Fa18 Out"};
    if (has("sp19") && has("inb")) return {"DVCS_Sp19_Inb", "Sp19 Inb"};
    if (has("sp18") && has("inb")) return {"DVCS_Sp18_Inb", "Sp18 Inb"};
    if (has("sp18") && has("out")) return {"DVCS_Sp18_Out", "Sp18 Out"};
    return {"", ""};
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

static inline std::string col_xBavg(const std::string& lab)  { return "xBavg, "     + lab; }
static inline std::string col_Q2avg(const std::string& lab)  { return "Q2avg, "     + lab; }
static inline std::string col_tavg (const std::string& lab)  { return "t_abs_avg, " + lab; }
static inline std::string col_phiavg(const std::string& lab) { return "phiavg, "    + lab; }

// ---------------- accumulators ----------------
struct Accum {
    double sx = 0.0, sQ = 0.0, st = 0.0, sp = 0.0;
    long long n = 0;
    void add(double x, double Q2, double tabs, double phi_deg) { sx += x; sQ += Q2; st += tabs; sp += phi_deg; ++n; }
    double mx() const { return n ? sx / n : std::numeric_limits<double>::quiet_NaN(); }
    double mQ() const { return n ? sQ / n : std::numeric_limits<double>::quiet_NaN(); }
    double mt() const { return n ? st / n : std::numeric_limits<double>::quiet_NaN(); }
    double mp() const { return n ? sp / n : std::numeric_limits<double>::quiet_NaN(); }
};

// ---------------- branch binder ----------------
struct BranchBinder {
    // global/exclusivity inputs
    double t1 = 0.0;               bool has_t1 = false;
    double open_angle_ep2 = 0.0;   bool has_open = false; // degrees
    double pTmiss = 0.0;           bool has_pT = false;
    double Emiss2 = 0.0;           bool has_Em2 = false;
    double Mx2 = 0.0;              bool has_Mx2 = false;
    double Mx2_1 = 0.0;            bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;            bool has_Mx2_2 = false;
    double xF = 0.0;               bool has_xF = false;

    // binning vars
    double x = 0.0;                bool has_x = false;
    double Q2 = 0.0;               bool has_Q2 = false;
    double phi2_rad = 0.0;         bool has_phi2 = false;  // radians

    void bind(TTree* t) {
        if (!t) return;
        auto bindD = [&](const char* n, double* a, bool& f){ if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };

        bindD("t1", &t1, has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_open);
        bindD("pTmiss", &pTmiss, has_pT);
        bindD("Emiss2", &Emiss2, has_Em2);
        bindD("Mx2", &Mx2, has_Mx2);
        bindD("Mx2_1", &Mx2_1, has_Mx2_1);
        bindD("Mx2_2", &Mx2_2, has_Mx2_2);
        bindD("xF", &xF, has_xF);

        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);
        bindD("phi2", &phi2_rad, has_phi2); // ONLY phi2 (radians)
    }

    bool readyForCuts() const {
        return has_t1 && has_open && has_pT;
    }
    bool readyForAverages() const {
        return has_x && has_Q2 && has_phi2;
    }
    double phi_deg() const {
        if (!has_phi2) return std::numeric_limits<double>::quiet_NaN();
        double deg = std::fmod(phi2_rad * RAD2DEG, 360.0);
        if (deg < 0.0) deg += 360.0;
        if (deg >= 360.0) deg = std::nextafter(360.0, 0.0);
        return deg;
    }
};

// ---------------- cuts ----------------
static inline bool passes_global(const BranchBinder& b) {
    if (!b.readyForCuts()) return false;
    if (b.open_angle_ep2 <= 5.0) return false;    // degrees
    if ((-b.t1) > 1.0)          return false;     // |t| < 1.0
    if (b.pTmiss > 0.20)        return false;     // GeV
    return true;
}

struct Sigmas { double mean = std::numeric_limits<double>::quiet_NaN(); double std = std::numeric_limits<double>::quiet_NaN(); };
using VarMap = std::unordered_map<std::string, Sigmas>;
static std::unordered_map<std::string, VarMap> g_sigma_cache;
static std::once_flag g_sigma_once;

static double parse_number(const std::string& s, size_t& i) {
    while (i < s.size() && !(std::isdigit((unsigned char)s[i]) || s[i]=='-' || s[i]=='+')) ++i;
    size_t st = i;
    while (i < s.size() && (std::isdigit((unsigned char)s[i]) || s[i]=='.' || s[i]=='e' || s[i]=='E' || s[i]=='+' || s[i]=='-')) ++i;
    if (st == i) return std::numeric_limits<double>::quiet_NaN();
    return std::strtod(s.c_str() + st, nullptr);
}
static Sigmas extract_mean_std(const std::string& block, const std::string& var) {
    Sigmas out;
    size_t data_pos = block.find("\"data\"");
    if (data_pos == std::string::npos) return out;
    size_t var_pos = block.find("\"" + var + "\"", data_pos);
    if (var_pos == std::string::npos) return out;
    size_t mpos = block.find("\"mean\"", var_pos);
    size_t spos = block.find("\"std\"",  var_pos);
    if (mpos != std::string::npos) { size_t i = mpos; out.mean = parse_number(block, i); }
    if (spos != std::string::npos) { size_t i = spos; out.std  = parse_number(block, i); }
    return out;
}
static std::string extract_object(const std::string& full, const std::string& key) {
    const std::string qk = "\"" + key + "\"";
    size_t p = full.find(qk);
    if (p == std::string::npos) return {};
    size_t brace = full.find('{', p);
    if (brace == std::string::npos) return {};
    int depth = 0;
    size_t i = brace;
    for (; i < full.size(); ++i) {
        if (full[i] == '{') ++depth;
        else if (full[i] == '}') { --depth; if (depth == 0) { ++i; break; } }
    }
    if (i <= brace) return {};
    return full.substr(brace, i - brace);
}
static void load_sigmas_once() {
    std::ifstream fin(kCutsJSON);
    if (!fin.is_open()) {
        std::cerr << "[bin_means] WARNING: cannot open " << kCutsJSON << " — 3-sigma cuts skipped." << std::endl;
        return;
    }
    std::string text((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

    const std::vector<std::string> P = { "DVCS_Fa18_Inb", "DVCS_Fa18_Out", "DVCS_Sp19_Inb",
                                         "DVCS_Sp18_Inb", "DVCS_Sp18_Out" };
    const std::vector<std::string> T = { "FD_FD", "CD_FD", "CD_FT" };
    const std::vector<std::string> VARS = { "Emiss2", "Mx2", "Mx2_1", "Mx2_2", "pTmiss", "theta_gamma_gamma", "xF" };

    for (const auto& p : P) {
        for (const auto& t : T) {
            const std::string key = p + "_" + t;
            const std::string obj = extract_object(text, key);
            if (obj.empty()) continue;
            VarMap vm;
            for (const auto& v : VARS) {
                Sigmas s = extract_mean_std(obj, v);
                if (std::isfinite(s.std)) vm[v] = s;
            }
            if (!vm.empty()) g_sigma_cache.emplace(key, std::move(vm));
        }
    }
}
enum class CutMode { TwoSided, UpperOnly };
static CutMode var_mode(const std::string& var) {
    if (var == "Emiss2" || var == "pTmiss" || var == "theta_gamma_gamma") return CutMode::UpperOnly;
    return CutMode::TwoSided;
}
static bool passes_3sigma_for_topo(const std::string& period_json_tag,
                                   Topology topo,
                                   const BranchBinder& b) {
    std::call_once(g_sigma_once, load_sigmas_once);
    if (period_json_tag.empty()) return true;

    const std::string key = period_json_tag + "_" + std::string(topo_tag(topo));
    auto it = g_sigma_cache.find(key);
    if (it == g_sigma_cache.end()) return true;

    const VarMap& vm = it->second;
    auto check = [&](const char* var, bool has_val, double val)->bool {
        auto itv = vm.find(var);
        if (itv == vm.end()) return true;
        if (!has_val) return true;
        const Sigmas s = itv->second;
        if (!std::isfinite(s.mean) || !std::isfinite(s.std)) return true;

        const CutMode mode = var_mode(var);
        if (mode == CutMode::UpperOnly) {
            return (val <= s.mean + 3.0 * s.std);
        } else {
            const double d = std::abs(val - s.mean);
            return (d <= 3.0 * s.std);
        }
    };

    const bool ok_theta = check("theta_gamma_gamma", b.has_open, b.open_angle_ep2 * DEG2RAD);
    const bool ok_pT    = check("pTmiss",           b.has_pT,   b.pTmiss);
    const bool ok_Em2   = check("Emiss2",           b.has_Em2,  b.Emiss2);
    const bool ok_Mx2   = check("Mx2",              b.has_Mx2,  b.Mx2);
    const bool ok_Mx21  = check("Mx2_1",            b.has_Mx2_1,b.Mx2_1);
    const bool ok_Mx22  = check("Mx2_2",            b.has_Mx2_2,b.Mx2_2);
    const bool ok_xF    = check("xF",               b.has_xF,   b.xF);

    return ok_theta && ok_pT && ok_Em2 && ok_Mx2 && ok_Mx21 && ok_Mx22 && ok_xF;
}
// Accept if ANY topology gate passes.
static bool passes_all_exclusivity(const std::string& period_json_tag, const BranchBinder& b) {
    return passes_3sigma_for_topo(period_json_tag, Topology::FD_FD, b)
        || passes_3sigma_for_topo(period_json_tag, Topology::CD_FD, b)
        || passes_3sigma_for_topo(period_json_tag, Topology::CD_FT, b);
}

static inline bool in_range(double v, double a, double b) { return (v >= a) && (v < b); }
static bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (pmax_deg > pmin_deg) return in_range(phi_deg, pmin_deg, pmax_deg);
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg); // wrap-around case
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

// ---------------- per period work ----------------
struct PeriodResult {
    std::unordered_map<int, Accum> per_row; // row index -> Accum
};

static PeriodResult process_period(const std::string& period_key, TTree* tree, const CSV& csv) {
    PeriodResult R;
    if (!tree) return R;

    PeriodKeyMap tags = period_key_to_tags(period_key);
    BranchBinder b; b.bind(tree);

    if (!b.readyForCuts() || !b.readyForAverages()) {
        std::cerr << "[bin_means] FATAL: Tree for '" << period_key
                  << "' missing branches (t1/open_angle_ep2/pTmiss/x/Q2/phi2)." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Fetch needed CSV columns once.
    const int c_xBmin = col(csv, "xBmin");
    const int c_xBmax = col(csv, "xBmax");
    const int c_Q2min = col(csv, "Q2min");
    const int c_Q2max = col(csv, "Q2max");
    const int c_tmin  = col(csv, "t_abs_min");
    const int c_tmax  = col(csv, "t_abs_max");
    const int c_pmin  = col(csv, "phimin");
    const int c_pmax  = col(csv, "phimax");
    const int c_valid = col(csv, "valid bin");

    auto toD = [](const std::string& s)->double {
        if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
        char* e = nullptr; double v = std::strtod(s.c_str(), &e);
        return (e == s.c_str()) ? std::numeric_limits<double>::quiet_NaN() : v;
    };

    // Pre-cache parsed CSV row windows to avoid string->double repeatedly.
    struct RowWin {
        double xBmin, xBmax, Q2min, Q2max, tmin, tmax, pmin, pmax;
        bool valid;
    };
    std::vector<RowWin> rows; rows.reserve(csv.rows.size());
    for (const auto& row : csv.rows) {
        RowWin w;
        w.xBmin = toD(row[c_xBmin]); w.xBmax = toD(row[c_xBmax]);
        w.Q2min = toD(row[c_Q2min]); w.Q2max = toD(row[c_Q2max]);
        w.tmin  = toD(row[c_tmin]);  w.tmax  = toD(row[c_tmax]);
        w.pmin  = toD(row[c_pmin]);  w.pmax  = toD(row[c_pmax]);
        const std::string& vs = row[c_valid];
        w.valid = (vs == "1" || vs == "1.0" || vs == "true" || vs == "TRUE");
        rows.push_back(w);
    }

    const Long64_t N = tree->GetEntries();
    print_banner("Processing period " + period_key + " with " + std::to_string((long long)N) + " entries");

    long long n_pass_global = 0;
    long long n_pass_3sig   = 0;
    long long n_used_rows   = 0;

    TStopwatch sw; sw.Start();

    // Progress cadence: print about every 2% or every 1M events, whichever comes first.
    const Long64_t cadence_by_pct = std::max<Long64_t>( (Long64_t) (0.02 * (double)N), 1 );
    const Long64_t cadence_by_abs = 1000000;
    const Long64_t cadence = std::min(cadence_by_pct, cadence_by_abs);

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i); // single pass through the tree

        if (!passes_global(b)) continue;
        ++n_pass_global;

        if (!passes_all_exclusivity(tags.json_tag, b)) continue;
        ++n_pass_3sig;

        const double phi_deg = b.phi_deg();
        bool used_any_row = false;

        // Test this event against CSV rows (no additional tree reads).
        for (int r = 0; r < (int)rows.size(); ++r) {
            const RowWin& w = rows[r];
            if (!w.valid) continue;
            if (!row_accepts_kin(b, w.xBmin, w.xBmax, w.Q2min, w.Q2max, w.tmin, w.tmax)) continue;
            if (!row_accepts_phi(phi_deg, w.pmin, w.pmax)) continue;

            R.per_row[r].add(b.x, b.Q2, std::fabs(b.t1), phi_deg);
            used_any_row = true;
        }
        if (used_any_row) ++n_used_rows;

        if (i == 0 || (i % cadence) == 0 || i + 1 == N) {
            double pct = (N > 0) ? (100.0 * (double)i / (double)N) : 100.0;
            print_status_singleline(period_key, pct, (long long)i, (long long)N,
                                    n_pass_global, n_pass_3sig, n_used_rows, sw.RealTime());
            sw.Continue(); // reset interval timing
        }
    }

    TStopwatch swtot; // simple call to show total not needed; we already printed periodic rates.
    print_banner("Finished " + period_key +
                 "  global_pass=" + std::to_string(n_pass_global) +
                 "  sig_pass="    + std::to_string(n_pass_3sig) +
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
    const auto* Sp19Inb = get("Sp19 Inb"); (void)Sp19Inb; // excluded from 10.6 GeV combo

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

    auto combine = [&](const std::vector<const std::unordered_map<int, Accum>*>& parts, int row)->Accum {
        Accum a;
        for (auto* p : parts) {
            if (!p) continue;
            auto it = p->find(row);
            if (it == p->end()) continue;
            const Accum& r = it->second;
            a.sx += r.sx; a.sQ += r.sQ; a.st += r.st; a.sp += r.sp; a.n += r.n;
        }
        return a;
    };

    print_banner("Combining period means into Fa18, Sp18, 10.6 GeV groups");
    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        {
            Accum a = combine({Fa18Inb, Fa18Out}, r);
            if (a.n > 0) {
                csv.rows[r][c_xB_Fa18]  = fmt8(a.mx());
                csv.rows[r][c_Q2_Fa18]  = fmt8(a.mQ());
                csv.rows[r][c_t_Fa18]   = fmt8(a.mt());
                csv.rows[r][c_phi_Fa18] = fmt8(a.mp());
            }
        }
        {
            Accum a = combine({Sp18Inb, Sp18Out}, r);
            if (a.n > 0) {
                csv.rows[r][c_xB_Sp18]  = fmt8(a.mx());
                csv.rows[r][c_Q2_Sp18]  = fmt8(a.mQ());
                csv.rows[r][c_t_Sp18]   = fmt8(a.mt());
                csv.rows[r][c_phi_Sp18] = fmt8(a.mp());
            }
        }
        {
            Accum a = combine({Fa18Inb, Fa18Out, Sp18Inb, Sp18Out}, r); // 10.6 GeV only
            if (a.n > 0) {
                csv.rows[r][c_xB_106]  = fmt8(a.mx());
                csv.rows[r][c_Q2_106]  = fmt8(a.mQ());
                csv.rows[r][c_t_106]   = fmt8(a.mt());
                csv.rows[r][c_phi_106] = fmt8(a.mp());
            }
        }
    }
}

bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers) {
    CSV csv;
    if (!load_csv(csv_path, csv)) return false;

    // Resolve period columns for output once.
    std::unordered_map<std::string,int> cxB, cQ2, ct, cphi;
    for (const auto& lab : csv_period_labels()) {
        cxB[lab]  = col(csv, col_xBavg(lab));
        cQ2[lab]  = col(csv, col_Q2avg(lab));
        ct[lab]   = col(csv, col_tavg (lab));
        cphi[lab] = col(csv, col_phiavg(lab));
    }

    // Determine available period keys.
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

    // Process each period (parallel across periods if requested).
    std::vector<PeriodResult> results(period_keys.size());

    const int nth = std::max(1, std::min(5, max_workers));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
    for (int i = 0; i < (int)period_keys.size(); ++i) {
        const std::string& pk = period_keys[i];
        auto it = dataTrees.find(pk);
        TTree* t = (it == dataTrees.end()) ? nullptr : it->second;
        results[i] = process_period(pk, t, csv);
    }

    // Write period means to CSV, and collect for group combos.
    std::unordered_map<std::string, std::unordered_map<int, Accum>> per_period_rows; // label -> row -> Accum

    for (int i = 0; i < (int)period_keys.size(); ++i) {
        const std::string& pk = period_keys[i];
        PeriodKeyMap tags = period_key_to_tags(pk);
        if (tags.csv_label.empty()) continue;

        const auto& perrow = results[i].per_row;
        per_period_rows[tags.csv_label] = perrow;

        long long filled = 0;
        for (const auto& kv : perrow) {
            int r = kv.first;
            const Accum& a = kv.second;
            if (a.n <= 0) continue;
            ++filled;

            csv.rows[r][cxB[tags.csv_label]]  = fmt8(a.mx());
            csv.rows[r][cQ2[tags.csv_label]]  = fmt8(a.mQ());
            csv.rows[r][ct [tags.csv_label]]  = fmt8(a.mt());
            csv.rows[r][cphi[tags.csv_label]] = fmt8(a.mp());
        }
        print_banner("Wrote " + std::to_string(filled) + " row means for " + tags.csv_label);
    }

    // Combined groups
    fill_combined_groups(csv, per_period_rows);

    if (!write_csv(csv_path, csv)) return false;
    print_banner(std::string("Updated bin means in: ") + csv_path);
    return true;
}