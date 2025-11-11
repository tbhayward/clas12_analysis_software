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

// -----------------------------------------------------------------------------
// Overview:
//  - Read CSV bins (xB, Q2, |t|, phi) and build bin lookups.
//  - For each group (5 periods + Fa18 + Sp18 + 10.6 GeV):
//      * Scan relevant TTree(s).
//      * Apply BOTH simple global cuts AND 3-sigma exclusivity cuts.
//      * Accept ALL topologies; topology is only used to select 3-sigma thresholds.
//      * Accumulate x, Q2, |t|, phi and counts per CSV row.
//      * Compute means and write them back into the group's four average columns.
//  - Write to a temp file, then rename to csv_path.
// Notes:
//  - Phi is ALWAYS in degrees here (both CSV bin edges and tree branches).
//  - Exclusivity 3-sigma thresholds are loaded from output/jsons/combined_cuts.json.
// -----------------------------------------------------------------------------

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
            for (char c : s) oss << (c == '"' ? "\"\"" : std::string(1, c));
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

// ---------------- Group wiring ----------------
static const std::vector<std::string> kGroups = {
    "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out",
    "Fa18", "Sp18", "10.6 GeV"
};

static const std::map<std::string,std::string> kDisplayToTree = {
    {"Fa18 Inb", "DVCS_Fa18_inb"},
    {"Fa18 Out", "DVCS_Fa18_out"},
    {"Sp19 Inb", "DVCS_Sp19_inb"},
    {"Sp18 Inb", "DVCS_Sp18_inb"},
    {"Sp18 Out", "DVCS_Sp18_out"}
};

static std::map<std::string, std::vector<std::string>> build_combined_groups() {
    std::map<std::string, std::vector<std::string>> m;
    m["Fa18"]     = {"Fa18 Inb", "Fa18 Out"};
    m["Sp18"]     = {"Sp18 Inb", "Sp18 Out"};
    m["10.6 GeV"] = {"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"};
    return m;
}

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

    // Phi always in DEGREES now.
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
    // key: "<JSONBASE>_<TOPO>", var: "Emiss2" / "Mx2" / "pTmiss" / "theta_gamma_gamma"
    std::map<std::pair<std::string,std::string>, MeanStd> data;
    bool loaded = false;
};

// Tiny numeric parser
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

static bool passes_3sigma_all(const std::string& tree_key, int det1, int det2, const BranchBinder& b) {
    const CutTable& ct = get_cuts_table();
    if (!ct.loaded) return true;

    auto itBase = kTreeToJsonBase.find(tree_key);
    if (itBase == kTreeToJsonBase.end()) return true;

    std::string topo = topo_string_from_det(det1, det2);
    if (topo.empty()) return true;

    std::string full = itBase->second + "_" + topo;

    const MeanStd ms_Emiss2 = ct.data.at({full, "Emiss2"});
    const MeanStd ms_Mx2    = ct.data.at({full, "Mx2"});
    const MeanStd ms_pT     = ct.data.at({full, "pTmiss"});
    const MeanStd ms_tgg    = ct.data.at({full, "theta_gamma_gamma"});

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
        if (v == r.second && i+1 == (int)rs.size()) return i;
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

struct GroupAccum {
    std::vector<Accum> bins; // one per CSV row
};

// ---------------- Scan one group ----------------
static void scan_group_fill(const std::string& group_name,
                            const std::vector<std::string>& group_display_periods,
                            const std::map<std::string, TTree*>& dataTrees,
                            const Ranges& rng,
                            const std::map<std::tuple<int,int,int,int>, int>& tuple_to_row,
                            GroupAccum& out)
{
    int max_row_index = -1;
    for (const auto& kv : tuple_to_row) if (kv.second > max_row_index) max_row_index = kv.second;
    int rows_count = max_row_index + 1;
    out.bins.assign(rows_count, Accum{});

    std::vector<std::pair<std::string,TTree*>> trees;
    for (const auto& disp : group_display_periods) {
        auto itT = kDisplayToTree.find(disp);
        if (itT == kDisplayToTree.end()) continue;
        auto itMap = dataTrees.find(itT->second);
        if (itMap != dataTrees.end() && itMap->second) trees.emplace_back(itT->second, itMap->second);
    }

    for (const auto& kv : trees) {
        const std::string& tree_key = kv.first;
        TTree* t = kv.second;
        BranchBinder b; b.bind(t);
        const Long64_t nent = t->GetEntries();

        for (Long64_t i = 0; i < nent; ++i) {
            t->GetEntry(i);

            if (!passes_simple_global(b)) continue;

            int d1 = b.has_detector1 ? b.detector1 : 0;
            int d2 = b.has_detector2 ? b.detector2 : 0;
            if (!passes_3sigma_all(tree_key, d1, d2, b)) continue;

            if (!(b.has_x && b.has_Q2 && b.has_t1)) continue;
            double xb  = b.x;
            double Q2  = b.Q2;
            double tab = std::fabs(b.t1);
            if (!std::isfinite(xb) || !std::isfinite(Q2) || !std::isfinite(tab)) continue;

            double phi_deg = b.phi_deg();
            if (!std::isfinite(phi_deg)) continue;
            phi_deg = wrap_360(phi_deg);

            int ix = find_bin(xb,  rng.xb);
            int iQ = find_bin(Q2,  rng.q2);
            int it = find_bin(tab, rng.tabs);
            int ip = find_bin(phi_deg, rng.phi);
            if (ix < 0 || iQ < 0 || it < 0 || ip < 0) continue;

            auto itRow = tuple_to_row.find({ix,iQ,it,ip});
            if (itRow == tuple_to_row.end()) continue;
            int row_index = itRow->second;
            if (row_index < 0 || row_index >= (int)out.bins.size()) continue;

            Accum& A = out.bins[row_index];
            A.sum_x   += xb;
            A.sum_Q2  += Q2;
            A.sum_ta  += tab;
            A.sum_phi += phi_deg;
            A.n       += 1;
        }
    }
}

// ---------------- Write means into CSV rows ----------------
static void write_means_into_rows(const std::string& group_name,
                                  const GroupAccum& gacc,
                                  std::vector<std::vector<std::string>>& rows,
                                  const std::unordered_map<std::string,int>& H)
{
    auto need = [&](const std::string& col)->int {
        auto it = H.find(col);
        if (it == H.end()) {
            std::cerr << "[bin_means] FATAL: column missing: " << col << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return it->second;
    };

    const std::string xb_col   = "xBavg ("     + group_name + ")";
    const std::string q2_col   = "Q2avg ("     + group_name + ")";
    const std::string t_col    = "t_abs_avg (" + group_name + ")";
    const std::string phi_col  = "phiavg ("    + group_name + ")";

    const int XB  = need(xb_col);
    const int Q2  = need(q2_col);
    const int TA  = need(t_col);
    const int PHI = need(phi_col);

    std::ostringstream oss;
    oss.setf(std::ios::fixed);

    for (int r = 0; r < (int)rows.size(); ++r) {
        const Accum& A = (r < (int)gacc.bins.size()) ? gacc.bins[r] : Accum{};
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
        double phi_mean = wrap_360(A.sum_phi / (double)A.n);

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

    Ranges rng;
    std::map<std::tuple<int,int,int,int>, int> tuple_to_row;
    build_ranges_and_index(valid_rows, H, rng, tuple_to_row);

    auto combined = build_combined_groups();

    struct GroupPlan {
        std::string name;
        std::vector<std::string> periods; // display names
        GroupAccum accum;
    };
    std::vector<GroupPlan> plans;

    // 5 individual periods we know about
    for (const auto& g : kGroups) {
        if (kDisplayToTree.count(g)) {
            plans.push_back(GroupPlan{g, std::vector<std::string>{g}, GroupAccum{}});
        }
    }
    // Combined groups
    for (const auto& kv : combined) {
        plans.push_back(GroupPlan{kv.first, kv.second, GroupAccum{}});
    }

    if (max_workers <= 0) max_workers = 1;
#ifdef _OPENMP
    #pragma omp parallel for num_threads(max_workers) schedule(dynamic)
#endif
    for (int ig = 0; ig < (int)plans.size(); ++ig) {
        scan_group_fill(plans[ig].name,
                        plans[ig].periods,
                        dataTrees,
                        rng,
                        tuple_to_row,
                        plans[ig].accum);
    }

    // Map valid-rows back into full CSV
    std::vector<std::vector<std::string>> working_valid = valid_rows;
    for (const auto& gp : plans) {
        write_means_into_rows(gp.name, gp.accum, working_valid, H);
    }
    for (int k = 0; k < (int)working_valid.size(); ++k) {
        rows_all[valid_row_indices[k]] = working_valid[k];
    }

    // Temp write then replace
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