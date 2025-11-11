#include "bin_means.h"

#include <TTree.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// Strategy
// --------
// 1) Load the CSV and identify valid rows (valid bin == 1). Build bin ranges
//    (xB, Q2, |t|, phi in degrees) from Sangbaek's CSV and a map from
//    (ix, iQ, it, ip) -> row_index.
// 2) Build the five periods and their tree keys from dataTrees:
//      "Fa18 Inb" -> "DVCS_Fa18_inb",
//      "Fa18 Out" -> "DVCS_Fa18_out",
//      "Sp19 Inb" -> "DVCS_Sp19_inb",
//      "Sp18 Inb" -> "DVCS_Sp18_inb",
//      "Sp18 Out" -> "DVCS_Sp18_out".
// 3) Parallelize ACROSS TREES (max_workers up to 5). Each worker scans its
//    TTree exactly once, applies BOTH global and 3-sigma exclusivity cuts,
//    accepts all topologies, and accumulates per-bin sums and counts.
// 4) After all workers finish, compute the combined groups (Fa18, Sp18,
//    10.6 GeV) by merging the per-period sums and counts (weighted means).
// 5) Write the four averages back into the CSV columns. Column names are
//    "base, Group" (comma + space), e.g. "xBavg, Fa18 Inb".
// ============================================================================

static const char* kCutsJSON = "output/jsons/combined_cuts.json";

// ---------------- CSV helpers ----------------
static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());
    bool in_quotes = false;
    for (char c : line) {
        if (c == '"') { in_quotes = !in_quotes; continue; }
        if (c == ',' && !in_quotes) { out.push_back(cur); cur.clear(); }
        else { cur.push_back(c); }
    }
    out.push_back(cur);
    return out;
}

static std::string join_csv_row(const std::vector<std::string>& fields) {
    std::ostringstream oss;
    for (size_t i = 0; i < fields.size(); ++i) {
        const std::string& s = fields[i];
        bool need_quotes = (s.find(',') != std::string::npos) || (s.find('"') != std::string::npos);
        if (need_quotes) {
            oss << '"';
            for (char c : s) { if (c == '"') oss << "\"\""; else oss << c; }
            oss << '"';
        } else {
            oss << s;
        }
        if (i + 1 < fields.size()) oss << ',';
    }
    return oss.str();
}

static std::unordered_map<std::string,int>
build_header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string,int> m;
    for (int i = 0; i < (int)header.size(); ++i) m[header[i]] = i;
    return m;
}

static std::string get_col(const std::vector<std::string>& row,
                           const std::unordered_map<std::string,int>& idx,
                           const std::string& name) {
    auto it = idx.find(name);
    if (it == idx.end()) return std::string();
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return std::string();
    return row[j];
}

static inline int ToInt(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

static inline double ToDouble(const std::string& s) {
    if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return std::numeric_limits<double>::quiet_NaN();
    return v;
}

// ---------------- Groups and period wiring ----------------
static const std::vector<std::string> kPeriodGroups = {
    "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
};

static const std::vector<std::string> kAllGroups = {
    "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out",
    "Fa18", "Sp18", "10.6 GeV"
};

// Map CSV display name -> canonical tree key used in dataTrees.
static const std::map<std::string,std::string> kDisplayToTree = {
    {"Fa18 Inb", "DVCS_Fa18_inb"},
    {"Fa18 Out", "DVCS_Fa18_out"},
    {"Sp19 Inb", "DVCS_Sp19_inb"},
    {"Sp18 Inb", "DVCS_Sp18_inb"},
    {"Sp18 Out", "DVCS_Sp18_out"}
};

// Combined membership (only display period names).
static std::map<std::string, std::vector<std::string>> build_combined_groups() {
    std::map<std::string, std::vector<std::string>> m;
    m["Fa18"]     = {"Fa18 Inb", "Fa18 Out"};
    m["Sp18"]     = {"Sp18 Inb", "Sp18 Out"};
    m["10.6 GeV"] = {"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"};
    return m;
}

// For 3-sigma JSON lookups: treeKey -> JSON base block name.
static const std::map<std::string,std::string> kTreeToJsonBase = {
    {"DVCS_Fa18_inb", "DVCS_Fa18_Inb"},
    {"DVCS_Fa18_out", "DVCS_Fa18_Out"},
    {"DVCS_Sp19_inb", "DVCS_Sp19_Inb"},
    {"DVCS_Sp18_inb", "DVCS_Sp18_Inb"},
    {"DVCS_Sp18_out", "DVCS_Sp18_Out"}
};

static inline std::string topo_string_from_det(int d1, int d2) {
    if (d1 == 1 && d2 == 1) return "FD_FD";
    if (d1 == 2 && d2 == 1) return "CD_FD";
    if (d1 == 2 && d2 == 0) return "CD_FT";
    return std::string();
}

// ---------------- Branch binder ----------------
struct BranchBinder {
    int    detector1 = 0; bool has_detector1 = false;
    int    detector2 = 0; bool has_detector2 = false;

    double t1 = 0.0;             bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open_angle_ep2 = false;
    double pTmiss = 0.0;         bool has_pTmiss = false;

    double x = 0.0;              bool has_x = false;
    double Q2 = 0.0;             bool has_Q2 = false;
    double phi2 = 0.0;           bool has_phi2 = false;
    double Delta_phi = 0.0;      bool has_Delta_phi = false;

    double Emiss2 = 0.0;         bool has_Emiss2 = false;
    double Mx2 = 0.0;            bool has_Mx2 = false;
    double theta_gamma_gamma = 0.0; bool has_tgg = false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI = [&](const char* n, int* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };
        auto bindD = [&](const char* n, double* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };

        bindI("detector1", &detector1, has_detector1);
        bindI("detector2", &detector2, has_detector2);

        bindD("t1", &t1, has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bindD("pTmiss", &pTmiss, has_pTmiss);

        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);
        bindD("phi2", &phi2, has_phi2);
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);

        bindD("Emiss2", &Emiss2, has_Emiss2);
        bindD("Mx2", &Mx2, has_Mx2);
        bindD("theta_gamma_gamma", &theta_gamma_gamma, has_tgg);
    }

    // Phi is ALWAYS in degrees here.
    double phi_deg() const {
        if (has_phi2) return phi2;
        if (has_Delta_phi) return Delta_phi;
        return std::numeric_limits<double>::quiet_NaN();
    }
};

// ---------------- Cuts ----------------
static inline bool passes_simple_global(const BranchBinder& b) {
    if (!(b.has_t1 && b.has_open_angle_ep2 && b.has_pTmiss)) return false;
    if (!(b.open_angle_ep2 > 5.0)) return false;   // deg
    if ((-b.t1) > 1.0) return false;               // |t| <= 1.0
    if (b.pTmiss > 0.20) return false;             // (GeV)
    return true;
}

struct MeanStd { double mean = 0.0; double std = 0.0; bool ok = false; };
struct CutTable {
    std::map<std::pair<std::string,std::string>, MeanStd> data; // {fullkey,var} -> mean,std
    bool loaded = false;
};

static bool parse_next_number(const std::string& s, size_t start, double& out, size_t& endpos) {
    size_t i = start;
    while (i < s.size() && std::isspace((unsigned char)s[i])) ++i;
    if (i >= s.size()) return false;
    if (s[i] == '+' || s[i] == '-') ++i;
    size_t j = i;
    bool dot = false, exp = false;
    while (j < s.size()) {
        char c = s[j];
        if (std::isdigit((unsigned char)c)) { ++j; continue; }
        if (!dot && c == '.') { dot = true; ++j; continue; }
        if (!exp && (c == 'e' || c == 'E')) {
            exp = true; ++j;
            if (j < s.size() && (s[j] == '+' || s[j] == '-')) ++j;
            continue;
        }
        break;
    }
    if (j == i) return false;
    out = std::strtod(s.substr(i, j - i).c_str(), nullptr);
    endpos = j;
    return true;
}

static MeanStd extract_mean_std(const std::string& json,
                                const std::string& key_prefix,
                                const std::string& var) {
    MeanStd ms;
    std::string anchor = "\"" + key_prefix + "\"";
    size_t p = json.find(anchor);
    if (p == std::string::npos) return ms;

    size_t q = json.find("\"data\"", p);
    if (q == std::string::npos) return ms;

    std::string varkey = "\"" + var + "\"";
    size_t r = json.find(varkey, q);
    if (r == std::string::npos) return ms;

    size_t mpos = json.find("\"mean\"", r);
    if (mpos == std::string::npos) return ms;
    mpos = json.find(':', mpos);
    if (mpos == std::string::npos) return ms;

    double mean_val = 0.0, std_val = 0.0;
    size_t endpos = mpos + 1;
    if (!parse_next_number(json, endpos, mean_val, endpos)) return ms;

    size_t spos = json.find("\"std\"", r);
    if (spos == std::string::npos) return ms;
    spos = json.find(':', spos);
    if (spos == std::string::npos) return ms;

    size_t endpos2 = spos + 1;
    if (!parse_next_number(json, endpos2, std_val, endpos2)) return ms;

    ms.mean = mean_val;
    ms.std  = std_val;
    ms.ok   = true;
    return ms;
}

static CutTable load_cuts_json_once() {
    CutTable table;
    std::ifstream fin(kCutsJSON);
    if (!fin.is_open()) {
        std::cerr << "[bin_means] Warning: cannot open " << kCutsJSON
                  << ". 3-sigma cuts will be skipped (global cuts still applied)." << std::endl;
        return table;
    }
    std::ostringstream oss; oss << fin.rdbuf();
    std::string js = oss.str();

    const char* vars[4] = {"Emiss2", "Mx2", "pTmiss", "theta_gamma_gamma"};
    std::vector<std::string> jsonBases;
    for (auto& kv : kTreeToJsonBase) jsonBases.push_back(kv.second);
    const char* topos[3] = {"FD_FD","CD_FD","CD_FT"};

    for (const auto& base : jsonBases) {
        for (const char* topo : topos) {
            std::string key = base + std::string("_") + topo;
            for (const char* v : vars) {
                table.data[{key, v}] = extract_mean_std(js, key, v);
            }
        }
    }
    table.loaded = true;
    return table;
}

static const CutTable& get_cuts_table() {
    static CutTable table = load_cuts_json_once();
    return table;
}

static inline bool pass_3sigma_one(const MeanStd& ms, double value, bool one_sided) {
    if (!ms.ok || !std::isfinite(value)) return true;
    if (ms.std <= 0.0) return true;
    if (one_sided) return (value <= (ms.mean + 3.0 * ms.std));
    return (std::fabs(value - ms.mean) <= 3.0 * ms.std);
}

static MeanStd get_ms(const std::string& fullkey, const std::string& var) {
    const CutTable& ct = get_cuts_table();
    auto it = ct.data.find({fullkey, var});
    if (it == ct.data.end()) return MeanStd{};
    return it->second;
}

static bool passes_3sigma_all(const std::string& tree_key, int det1, int det2, const BranchBinder& b) {
    const CutTable& ct = get_cuts_table();
    if (!ct.loaded) return true;

    auto itBase = kTreeToJsonBase.find(tree_key);
    if (itBase == kTreeToJsonBase.end()) return true;

    std::string topo = topo_string_from_det(det1, det2);
    if (topo.empty()) return true;

    std::string full = itBase->second + "_" + topo;

    const MeanStd ms_Emiss2 = get_ms(full, "Emiss2");
    const MeanStd ms_Mx2    = get_ms(full, "Mx2");
    const MeanStd ms_pT     = get_ms(full, "pTmiss");
    const MeanStd ms_tgg    = get_ms(full, "theta_gamma_gamma");

    bool ok = true;
    if (b.has_Emiss2)   ok = ok && pass_3sigma_one(ms_Emiss2, b.Emiss2, false);
    if (b.has_Mx2)      ok = ok && pass_3sigma_one(ms_Mx2,    b.Mx2,    false);
    if (b.has_pTmiss)   ok = ok && pass_3sigma_one(ms_pT,     b.pTmiss, true);
    if (b.has_tgg)      ok = ok && pass_3sigma_one(ms_tgg,    b.theta_gamma_gamma, true);
    return ok;
}

// ---------------- Bin helpers (phi in degrees) ----------------
struct Ranges {
    std::vector<std::pair<double,double>> xb;
    std::vector<std::pair<double,double>> q2;
    std::vector<std::pair<double,double>> tabs;
    std::vector<std::pair<double,double>> phi; // degrees
};

static int find_bin(double v, const std::vector<std::pair<double,double>>& rs) {
    for (int i = 0; i < (int)rs.size(); ++i) {
        const auto& r = rs[i];
        if (r.first <= v && v < r.second) return i;
        if (v == r.second && i+1 == (int)rs.size()) return i; // include last-edge
    }
    return -1;
}

static double wrap_360(double d) {
    double w = std::fmod(d, 360.0);
    if (w < 0.0) w += 360.0;
    return w;
}

static void build_ranges_and_index(
    const std::vector<std::vector<std::string>>& rows,
    const std::unordered_map<std::string,int>& H,
    Ranges& rng,
    std::map<std::tuple<int,int,int,int>, int>& tuple_to_row)
{
    std::set<std::pair<double,double>> sxb, sq2, st, sphi;

    for (int row_index = 0; row_index < (int)rows.size(); ++row_index) {
        const auto& r = rows[row_index];
        double xbmin = ToDouble(get_col(r, H, "xBmin"));
        double xbmax = ToDouble(get_col(r, H, "xBmax"));
        double q2min = ToDouble(get_col(r, H, "Q2min"));
        double q2max = ToDouble(get_col(r, H, "Q2max"));
        double tmin  = ToDouble(get_col(r, H, "t_abs_min"));
        double tmax  = ToDouble(get_col(r, H, "t_abs_max"));
        double phimin = ToDouble(get_col(r, H, "phimin"));
        double phimax = ToDouble(get_col(r, H, "phimax"));

        sxb.emplace(xbmin, xbmax);
        sq2.emplace(q2min, q2max);
        st.emplace(tmin, tmax);
        sphi.emplace(phimin, phimax);
    }

    rng.xb.assign(sxb.begin(), sxb.end());
    rng.q2.assign(sq2.begin(), sq2.end());
    rng.tabs.assign(st.begin(), st.end());
    rng.phi.assign(sphi.begin(), sphi.end());

    for (int row_index = 0; row_index < (int)rows.size(); ++row_index) {
        const auto& r = rows[row_index];
        double xbmin = ToDouble(get_col(r, H, "xBmin"));
        double xbmax = ToDouble(get_col(r, H, "xBmax"));
        double q2min = ToDouble(get_col(r, H, "Q2min"));
        double q2max = ToDouble(get_col(r, H, "Q2max"));
        double tmin  = ToDouble(get_col(r, H, "t_abs_min"));
        double tmax  = ToDouble(get_col(r, H, "t_abs_max"));
        double phimin = ToDouble(get_col(r, H, "phimin"));
        double phimax = ToDouble(get_col(r, H, "phimax"));

        int ix = find_bin(xbmin + 0.5 * (xbmax - xbmin), rng.xb);
        int iQ = find_bin(q2min + 0.5 * (q2max - q2min), rng.q2);
        int it = find_bin(tmin  + 0.5 * (tmax  - tmin ), rng.tabs);

        double midphi = wrap_360(phimin + 0.5 * (phimax - phimin));
        int ip = find_bin(midphi, rng.phi);

        if (ix >= 0 && iQ >= 0 && it >= 0 && ip >= 0) {
            tuple_to_row[{ix,iQ,it,ip}] = row_index;
        }
    }
}

// ---------------- Accumulators ----------------
struct Accum {
    double sum_x = 0.0;
    double sum_Q2 = 0.0;
    double sum_ta = 0.0;
    double sum_phi = 0.0; // degrees
    long long n = 0;
};

struct PeriodAccum {
    std::vector<Accum> bins; // one per valid CSV row
};

// ---------------- Column name resolution (comma style) ----------------
// Preferred header form is: "base, Group" (comma + space).
// For robustness, we also accept a couple of legacy aliases.
static int find_avg_col(const std::unordered_map<std::string,int>& H,
                        const std::string& base, const std::string& group) {
    const std::string primary = base + ", " + group;   // matches initialize_pass2_csv
    auto it = H.find(primary);
    if (it != H.end()) return it->second;

    // Fallbacks (if someone imports an older CSV)
    const std::string alt1 = base + " (" + group + ")";
    const std::string alt2 = base + "," + group;       // comma no space
    const std::string alt3 = base + " " + group;       // space no comma

    auto try1 = H.find(alt1); if (try1 != H.end()) return try1->second;
    auto try2 = H.find(alt2); if (try2 != H.end()) return try2->second;
    auto try3 = H.find(alt3); if (try3 != H.end()) return try3->second;

    std::cerr << "[bin_means] FATAL: column missing: " << primary
              << " (also tried \"" << alt1 << "\", \"" << alt2 << "\", \"" << alt3 << "\")"
              << std::endl;
    std::exit(EXIT_FAILURE);
}

// ---------------- Write means into CSV rows ----------------
static void write_means_into_rows(const std::string& group_name,
                                  const std::vector<Accum>& bins,
                                  std::vector<std::vector<std::string>>& rows,
                                  const std::unordered_map<std::string,int>& H)
{
    const int XB  = find_avg_col(H, "xBavg",     group_name);
    const int Q2  = find_avg_col(H, "Q2avg",     group_name);
    const int TA  = find_avg_col(H, "t_abs_avg", group_name);
    const int PHI = find_avg_col(H, "phiavg",    group_name);

    std::ostringstream oss;
    oss.setf(std::ios::fixed);

    for (int r = 0; r < (int)rows.size(); ++r) {
        const Accum& A = (r < (int)bins.size()) ? bins[r] : Accum{};
        if (A.n <= 0) {
            rows[r][XB].clear();
            rows[r][Q2].clear();
            rows[r][TA].clear();
            rows[r][PHI].clear();
            continue;
        }
        double xb_mean  = A.sum_x   / (double)A.n;
        double q2_mean  = A.sum_Q2  / (double)A.n;
        double ta_mean  = A.sum_ta  / (double)A.n;
        double phi_mean = wrap_360(A.sum_phi / (double)A.n); // degrees

        oss.str(""); oss << std::setprecision(8) << xb_mean;  rows[r][XB]  = oss.str();
        oss.str(""); oss << std::setprecision(8) << q2_mean;  rows[r][Q2]  = oss.str();
        oss.str(""); oss << std::setprecision(8) << ta_mean;  rows[r][TA]  = oss.str();
        oss.str(""); oss << std::setprecision(8) << phi_mean; rows[r][PHI] = oss.str();
    }
}

// ---------------- Public API ----------------
bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers)
{
    // Load CSV
    std::ifstream fin(csv_path);
    if (!fin.is_open()) {
        std::cerr << "[bin_means] ERROR: cannot open CSV: " << csv_path << std::endl;
        return false;
    }
    std::string header_line;
    if (!std::getline(fin, header_line)) {
        std::cerr << "[bin_means] ERROR: CSV is empty: " << csv_path << std::endl;
        return false;
    }
    std::vector<std::string> header = split_csv_line(header_line);
    auto H = build_header_index(header);

    std::vector<std::vector<std::string>> rows_all;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) { rows_all.push_back(std::vector<std::string>{}); continue; }
        rows_all.push_back(split_csv_line(line));
    }
    fin.close();

    int col_valid = (H.count("valid bin") ? H.at("valid bin") : -1);
    std::vector<int> valid_row_indices;
    for (int i = 0; i < (int)rows_all.size(); ++i) {
        const auto& r = rows_all[i];
        if (col_valid < 0 || col_valid >= (int)r.size()) continue;
        if (ToInt(r[col_valid]) == 1) valid_row_indices.push_back(i);
    }

    if (valid_row_indices.empty()) {
        std::cerr << "[bin_means] WARNING: no valid rows (valid bin == 1). Nothing to do." << std::endl;
        return true;
    }

    std::vector<std::vector<std::string>> valid_rows;
    valid_rows.reserve(valid_row_indices.size());
    for (int idx : valid_row_indices) valid_rows.push_back(rows_all[idx]);

    // Bin ranges and row index mapping
    struct Ranges rng;
    std::map<std::tuple<int,int,int,int>, int> tuple_to_row;
    build_ranges_and_index(valid_rows, H, rng, tuple_to_row);

    // Prepare per-period accumulators (only the five periods are scanned).
    const int n_valid = (int)valid_rows.size();
    std::map<std::string, PeriodAccum> period_acc; // key: display period name

    for (const auto& disp : kPeriodGroups) {
        period_acc[disp].bins.assign(n_valid, Accum{});
    }

    // Build list of trees we actually have.
    struct PeriodTree {
        std::string display;   // e.g., "Fa18 Inb"
        std::string tree_key;  // e.g., "DVCS_Fa18_inb"
    };
    std::vector<PeriodTree> trees_to_scan;
    trees_to_scan.reserve(kPeriodGroups.size());
    for (const auto& disp : kPeriodGroups) {
        auto itKey = kDisplayToTree.find(disp);
        if (itKey == kDisplayToTree.end()) continue;
        auto itTree = dataTrees.find(itKey->second);
        if (itTree == dataTrees.end() || !itTree->second) continue;
        trees_to_scan.push_back({disp, itKey->second});
    }

    if (max_workers <= 0) max_workers = 1;

    // Parallelize across trees; each TTree is read by at most one thread.
#ifdef _OPENMP
    #pragma omp parallel for num_threads(std::min<int>(max_workers, (int)trees_to_scan.size())) schedule(dynamic)
#endif
    for (int itree = 0; itree < (int)trees_to_scan.size(); ++itree) {
        const PeriodTree PT = trees_to_scan[itree];
        TTree* t = dataTrees.at(PT.tree_key);

        // Local accumulator for this period (thread-local).
        std::vector<Accum> local_bins(n_valid);

        BranchBinder b; b.bind(t);
        const Long64_t nent = t->GetEntries();

        for (Long64_t i = 0; i < nent; ++i) {
            t->GetEntry(i);

            // Apply BOTH cut sets.
            if (!passes_simple_global(b)) continue;

            int d1 = b.has_detector1 ? b.detector1 : 0;
            int d2 = b.has_detector2 ? b.detector2 : 0;
            if (!passes_3sigma_all(PT.tree_key, d1, d2, b)) continue;

            if (!(b.has_x && b.has_Q2 && b.has_t1)) continue;

            double xb  = b.x;
            double Q2  = b.Q2;
            double tab = std::fabs(b.t1);
            double phi_deg = b.phi_deg();

            if (!std::isfinite(xb) || !std::isfinite(Q2) || !std::isfinite(tab) || !std::isfinite(phi_deg)) continue;

            phi_deg = wrap_360(phi_deg);

            int ix = find_bin(xb,  rng.xb);
            int iQ = find_bin(Q2,  rng.q2);
            int it = find_bin(tab, rng.tabs);
            int ip = find_bin(phi_deg, rng.phi);
            if (ix < 0 || iQ < 0 || it < 0 || ip < 0) continue;

            auto rit = tuple_to_row.find({ix,iQ,it,ip});
            if (rit == tuple_to_row.end()) continue;
            int row_index = rit->second;

            Accum& A = local_bins[row_index];
            A.sum_x   += xb;
            A.sum_Q2  += Q2;
            A.sum_ta  += tab;
            A.sum_phi += phi_deg;
            A.n       += 1;
        }

        // Merge into the period accumulator.
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            auto& dst = period_acc[PT.display].bins;
            for (int r = 0; r < n_valid; ++r) {
                dst[r].sum_x   += local_bins[r].sum_x;
                dst[r].sum_Q2  += local_bins[r].sum_Q2;
                dst[r].sum_ta  += local_bins[r].sum_ta;
                dst[r].sum_phi += local_bins[r].sum_phi;
                dst[r].n       += local_bins[r].n;
            }
        }
    }

    // Build combined groups from period sums (weighted means via summed counts).
    std::map<std::string, std::vector<Accum>> all_group_bins;
    // Start by copying the five periods straight through.
    for (const auto& disp : kPeriodGroups) {
        all_group_bins[disp] = period_acc[disp].bins;
    }

    // Now the combined ones.
    auto combined = build_combined_groups();
    for (const auto& kv : combined) {
        const std::string& cname = kv.first;
        const auto& members = kv.second;
        std::vector<Accum> bins(n_valid);
        for (const auto& mdisp : members) {
            const auto& src = period_acc[mdisp].bins;
            for (int r = 0; r < n_valid; ++r) {
                bins[r].sum_x   += src[r].sum_x;
                bins[r].sum_Q2  += src[r].sum_Q2;
                bins[r].sum_ta  += src[r].sum_ta;
                bins[r].sum_phi += src[r].sum_phi;
                bins[r].n       += src[r].n;
            }
        }
        all_group_bins[cname] = std::move(bins);
    }

    // Map valid-rows back into full CSV after writing means for each group.
    std::vector<std::vector<std::string>> working_valid = valid_rows;
    for (const auto& gname : kAllGroups) {
        auto it = all_group_bins.find(gname);
        if (it == all_group_bins.end()) continue;
        write_means_into_rows(gname, it->second, working_valid, H);
    }
    for (int k = 0; k < (int)working_valid.size(); ++k) {
        rows_all[valid_row_indices[k]] = working_valid[k];
    }

    // Write atomically
    std::string tmp_path = csv_path + ".tmp.binmeans";
    {
        std::ofstream fout(tmp_path);
        if (!fout.is_open()) {
            std::cerr << "[bin_means] ERROR: cannot open temp CSV: " << tmp_path << std::endl;
            return false;
        }
        fout << join_csv_row(header) << "\n";
        for (const auto& r : rows_all) {
            if (r.empty()) { fout << "\n"; continue; }
            fout << join_csv_row(r) << "\n";
        }
    }
    if (std::rename(tmp_path.c_str(), csv_path.c_str()) != 0) {
        std::perror("[bin_means] rename failed");
        return false;
    }

    std::cout << "[bin_means] Updated bin means in: " << csv_path << std::endl;
    return true;
}