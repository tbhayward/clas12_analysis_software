// bin_means.cpp
// Compute per-bin means (xB, Q2, |t|, phi) per period, then fill combined groups,
// using global cuts AND all derived 3-sigma exclusivity cuts from combined_cuts.json.

#include "bin_means.h"
#include "periods.h"
#include "load_binning_scheme.h"

#include <TTree.h>
#include <TChain.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cctype>
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

// ---------------- constants ----------------
static constexpr double PI = 3.14159265358979323846;
static constexpr double DEG2RAD = PI / 180.0;
static const std::string kCutsJSON = "output/jsons/combined_cuts.json";

// We accept all topologies for the averages; we only use them to look up the 3-sigma gates.
enum class Topology { FD_FD, CD_FD, CD_FT };

static inline const char* topo_tag(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }
    return "FD_FD";
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

// Map your tree key to the JSON block and CSV label.
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

// ---------------- CSV utils ----------------
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
            for (char ch : s) {
                if (ch == '"') fout << "\"\"";
                else           fout << ch;
            }
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
                for (char ch : s) {
                    if (ch == '"') fout << "\"\"";
                    else           fout << ch;
                }
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

static inline std::string col_xBavg(const std::string& group)  { return "xBavg, "    + group; }
static inline std::string col_Q2avg(const std::string& group)  { return "Q2avg, "    + group; }
static inline std::string col_tavg (const std::string& group)  { return "t_abs_avg, "+ group; }
static inline std::string col_phiavg(const std::string& group) { return "phiavg, "   + group; }

// ---------------- tuple hasher ----------------
struct TupleHash {
    template <class T> static inline void hc(std::size_t& s, const T& v) {
        std::hash<T> h; s ^= h(v) + 0x9e3779b97f4a7c15ULL + (s<<6) + (s>>2);
    }
    template <class... Ts> std::size_t operator()(const std::tuple<Ts...>& t) const {
        std::size_t seed = 0;
        std::apply([&](const Ts&... elems){ (hc(seed, elems), ...); }, t);
        return seed;
    }
};

// ---------------- branches ----------------
struct BranchBinder {
    // det flags (bound in case you want later logic)
    int detector1 = 0; bool has_d1 = false;
    int detector2 = 0; bool has_d2 = false;

    // exclusivity & globals
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
    double phi2 = 0.0;             bool has_phi = false;  // degrees

    void bind(TTree* t) {
        if (!t) return;

        // If this is a TChain, force-load the first tree so branches are discoverable.
        if (auto ch = dynamic_cast<TChain*>(t)) {
            ch->LoadTree(0);   // does not require addresses set
        }

        auto bindI = [&](const char* n, int* a, bool& f) {
            // Try to bind regardless; mark present if a branch or leaf with this name exists
            t->SetBranchAddress(n, a);
            f = (t->GetBranch(n) != nullptr) || (t->GetLeaf(n) != nullptr);
        };
        auto bindD = [&](const char* n, double* a, bool& f) {
            t->SetBranchAddress(n, a);
            f = (t->GetBranch(n) != nullptr) || (t->GetLeaf(n) != nullptr);
        };

        bindI("detector1", &detector1, has_d1);
        bindI("detector2", &detector2, has_d2);

        bindD("t1",               &t1,               has_t1);
        bindD("open_angle_ep2",   &open_angle_ep2,   has_open);   // degrees
        bindD("pTmiss",           &pTmiss,           has_pT);
        bindD("Emiss2",           &Emiss2,           has_Em2);
        bindD("Mx2",              &Mx2,              has_Mx2);
        bindD("Mx2_1",            &Mx2_1,            has_Mx2_1);
        bindD("Mx2_2",            &Mx2_2,            has_Mx2_2);
        bindD("xF",               &xF,               has_xF);

        bindD("x",                &x,                has_x);
        bindD("Q2",               &Q2,               has_Q2);
        bindD("phi2",             &phi2,             has_phi);    // degrees
    }

    bool readyForCuts() const {
        // require core globals; the 3-sigma variables are optional and applied if present
        return has_t1 && has_open && has_pT;
    }
    bool readyForAverages() const {
        return has_x && has_Q2 && has_phi;
    }
};

// ---------------- global cuts ----------------
static inline bool passes_global(const BranchBinder& b) {
    if (!b.readyForCuts()) return false;
    if (b.open_angle_ep2 <= 5.0) return false;    // deg
    if ((-b.t1) > 1.0)          return false;     // |t| < 1.0
    if (b.pTmiss > 0.20)        return false;     // GeV
    return true;
}

// ---------------- 3-sigma cuts loader ----------------
struct Sigmas { double mean = std::numeric_limits<double>::quiet_NaN(); double std = std::numeric_limits<double>::quiet_NaN(); };
using VarMap = std::unordered_map<std::string, Sigmas>;           // var -> {mean,std}
static std::unordered_map<std::string, VarMap> g_sigma_cache;     // "DVCS_Fa18_Inb_FD_FD" -> map
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

    const std::vector<std::string> VARS = {
        "Emiss2", "Mx2", "Mx2_1", "Mx2_2", "pTmiss", "theta_gamma_gamma", "xF"
    };

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
    return CutMode::TwoSided; // Mx2, Mx2_1, Mx2_2, xF
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
        if (itv == vm.end()) return true;          // no derived cut for this var -> pass
        if (!has_val) return true;                 // branch not available -> pass
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

    // theta: JSON stored in radians; branch is degrees
    const bool ok_theta = check("theta_gamma_gamma", b.has_open, b.open_angle_ep2 * DEG2RAD);
    const bool ok_pT    = check("pTmiss",           b.has_pT,   b.pTmiss);
    const bool ok_Em2   = check("Emiss2",           b.has_Em2,  b.Emiss2);
    const bool ok_Mx2   = check("Mx2",              b.has_Mx2,  b.Mx2);
    const bool ok_Mx21  = check("Mx2_1",            b.has_Mx2_1,b.Mx2_1);
    const bool ok_Mx22  = check("Mx2_2",            b.has_Mx2_2,b.Mx2_2);
    const bool ok_xF    = check("xF",               b.has_xF,   b.xF);

    return ok_theta && ok_pT && ok_Em2 && ok_Mx2 && ok_Mx21 && ok_Mx22 && ok_xF;
}

// Accept if ANY topology’s 3-sigma gate passes (accept-all-topologies for averages).
static bool passes_all_exclusivity(const std::string& period_json_tag, const BranchBinder& b) {
    return passes_3sigma_for_topo(period_json_tag, Topology::FD_FD, b)
        || passes_3sigma_for_topo(period_json_tag, Topology::CD_FD, b)
        || passes_3sigma_for_topo(period_json_tag, Topology::CD_FT, b);
}

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

static bool in_range(double v, double a, double b) { return (v >= a) && (v < b); }

static bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (!std::isfinite(phi_deg) || !std::isfinite(pmin_deg) || !std::isfinite(pmax_deg)) return false;
    if (pmax_deg > pmin_deg) return in_range(phi_deg, pmin_deg, pmax_deg);
    // wrap-around safety
    return phi_deg >= pmin_deg || phi_deg < pmax_deg;
}

static bool row_accepts_kin(const BranchBinder& b,
                            double xBmin, double xBmax,
                            double Q2min, double Q2max,
                            double tmin, double tmax) {
    const double tabs = std::fabs(b.t1);
    if (!in_range(b.x,   xBmin, xBmax))   return false;
    if (!in_range(b.Q2,  Q2min, Q2max))   return false;
    if (!in_range(tabs,  tmin,  tmax))    return false;
    return true;
}

// ---------------- per-period processing ----------------
struct PeriodResult {
    std::unordered_map<int, Accum> per_row; // row index -> Accum
};

static PeriodResult process_period(const std::string& period_key, TTree* tree, const CSV& csv) {
    PeriodResult R;
    if (!tree) return R;

    PeriodKeyMap tags = period_key_to_tags(period_key);
    if (tags.json_tag.empty() || tags.csv_label.empty()) {
        std::cerr << "[bin_means] WARNING: cannot map period key '" << period_key
                  << "' to JSON/CSV labels; proceeding without 3-sigma cuts for this period."
                  << std::endl;
    }

    BranchBinder b; b.bind(tree);
    if (!b.readyForCuts() || !b.readyForAverages()) {
        std::cerr << "[bin_means] FATAL: Tree for '" << period_key
                  << "' missing branches (t1/open_angle_ep2/pTmiss/x/Q2/phi2)." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // CSV columns we need
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

    const Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        if (!passes_global(b)) continue;

        // 3-sigma (any topology OK)
        if (!passes_all_exclusivity(tags.json_tag, b)) continue;

        for (int r = 0; r < (int)csv.rows.size(); ++r) {
            const auto& row = csv.rows[r];

            // valid row?
            bool valid = true;
            if (c_valid < (int)row.size()) {
                const std::string& vs = row[c_valid];
                valid = (vs == "1" || vs == "1.0" || vs == "true" || vs == "TRUE");
            }
            if (!valid) continue;

            const double xBmin  = toD(row[c_xBmin]);
            const double xBmax  = toD(row[c_xBmax]);
            const double Q2min  = toD(row[c_Q2min]);
            const double Q2max  = toD(row[c_Q2max]);
            const double tmin   = toD(row[c_tmin]);
            const double tmax   = toD(row[c_tmax]);
            const double pmin   = toD(row[c_pmin]);
            const double pmax   = toD(row[c_pmax]);

            if (!row_accepts_kin(b, xBmin, xBmax, Q2min, Q2max, tmin, tmax)) continue;
            if (!row_accepts_phi(b.phi2, pmin, pmax)) continue;

            R.per_row[r].add(b.x, b.Q2, std::fabs(b.t1), b.phi2);
        }
    }

    return R;
}

// ---------------- main entry ----------------
static inline std::string fmt8(double v) {
    if (!std::isfinite(v)) return std::string();
    std::ostringstream oss; oss << std::fixed << std::setprecision(8) << v;
    return oss.str();
}

static void fill_combined_groups(CSV& csv,
                                 const std::unordered_map<std::string, std::unordered_map<int, Accum>>& per_period) {
    // Build helpers: map CSV label -> per-row Accum (already computed)
    auto get = [&](const std::string& lab)->const std::unordered_map<int, Accum>* {
        auto it = per_period.find(lab);
        return (it == per_period.end()) ? nullptr : &it->second;
    };

    const auto* Fa18Inb = get("Fa18 Inb");
    const auto* Fa18Out = get("Fa18 Out");
    const auto* Sp18Inb = get("Sp18 Inb");
    const auto* Sp18Out = get("Sp18 Out");
    const auto* Sp19Inb = get("Sp19 Inb");

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

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        // Fa18 = Inb + Out
        {
            Accum a = combine({Fa18Inb, Fa18Out}, r);
            if (a.n > 0) {
                csv.rows[r][c_xB_Fa18]  = fmt8(a.mx());
                csv.rows[r][c_Q2_Fa18]  = fmt8(a.mQ());
                csv.rows[r][c_t_Fa18]   = fmt8(a.mt());
                csv.rows[r][c_phi_Fa18] = fmt8(a.mp());
            }
        }
        // Sp18 = Inb + Out
        {
            Accum a = combine({Sp18Inb, Sp18Out}, r);
            if (a.n > 0) {
                csv.rows[r][c_xB_Sp18]  = fmt8(a.mx());
                csv.rows[r][c_Q2_Sp18]  = fmt8(a.mQ());
                csv.rows[r][c_t_Sp18]   = fmt8(a.mt());
                csv.rows[r][c_phi_Sp18] = fmt8(a.mp());
            }
        }
        // 10.6 GeV = Fa18 + Sp18 (exclude Sp19 10.2 GeV)
        {
            Accum a = combine({Fa18Inb, Fa18Out, Sp18Inb, Sp18Out}, r);
            if (a.n > 0) {
                csv.rows[r][c_xB_106]  = fmt8(a.mx());
                csv.rows[r][c_Q2_106]  = fmt8(a.mQ());
                csv.rows[r][c_t_106]   = fmt8(a.mt());
                csv.rows[r][c_phi_106] = fmt8(a.mp());
            }
        }
        (void)Sp19Inb; // intentionally unused for 10.6 GeV group
    }
}

bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers) {
    CSV csv;
    if (!load_csv(csv_path, csv)) return false;

    // Resolve per-period CSV columns for writing period means
    std::unordered_map<std::string,int> cxB, cQ2, ct, cphi;
    for (const auto& lab : csv_period_labels()) {
        cxB[lab]  = col(csv, col_xBavg(lab));
        cQ2[lab]  = col(csv, col_Q2avg(lab));
        ct[lab]   = col(csv, col_tavg (lab));
        cphi[lab] = col(csv, col_phiavg(lab));
    }

    // Determine available period keys from dataTrees
    std::vector<std::string> period_keys;
    for (const auto& P : CANONICAL_PERIODS()) {
        auto it = dataTrees.find(P.tree_key);
        if (it != dataTrees.end() && it->second) period_keys.push_back(P.tree_key);
    }
    if (period_keys.empty()) {
        std::cerr << "[bin_means] FATAL: no DVCS trees available in dataTrees." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Process periods (optionally in parallel across periods)
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

    // Write period results into CSV cells
    std::unordered_map<std::string, std::unordered_map<int, Accum>> per_period_rows; // label -> row -> Accum

    for (int i = 0; i < (int)period_keys.size(); ++i) {
        const std::string& pk = period_keys[i];
        PeriodKeyMap tags = period_key_to_tags(pk);
        if (tags.csv_label.empty()) continue; // skip unknown

        const auto& perrow = results[i].per_row;
        per_period_rows[tags.csv_label] = perrow;

        for (const auto& kv : perrow) {
            int r = kv.first;
            const Accum& a = kv.second;
            if (a.n <= 0) continue;

            csv.rows[r][cxB[tags.csv_label]]  = fmt8(a.mx());
            csv.rows[r][cQ2[tags.csv_label]]  = fmt8(a.mQ());
            csv.rows[r][ct [tags.csv_label]]  = fmt8(a.mt());
            csv.rows[r][cphi[tags.csv_label]] = fmt8(a.mp());
        }
    }

    // Fill combined groups from the period maps (no TTree reread)
    fill_combined_groups(csv, per_period_rows);

    if (!write_csv(csv_path, csv)) return false;
    std::cout << "[bin_means] Updated bin means in: " << csv_path << std::endl;
    return true;
}