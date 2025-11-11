// bin_means.cpp
//
// Computes per-row (bin) means of xB, Q2, |t|, phi and writes them
// back into the analysis CSV, using BOTH global simple cuts and
// 3-sigma exclusivity windows loaded from output/jsons/combined_cuts.json.
// We accept all topologies for the averages; topology is only used to
// select the proper 3-sigma windows from the JSON.
//
// Dependencies: nlohmann::json single header.

#include "bin_means.h"
#include "periods.h"              // CANONICAL_PERIODS()
#include "load_binning_scheme.h"  // Binning (only for signature parity)

#include <TTree.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cctype>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef NLOHMANN_JSON_HPP
  #include <nlohmann/json.hpp>
#endif

using json = nlohmann::json;

// ---------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------

static const char* kCombinedCutsPath = "output/jsons/combined_cuts.json";

// CSV column names for averages we will fill
static const std::vector<std::string> kPeriodsPretty = {
    "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
};
static const std::vector<std::string> kCombinedPretty = {
    "Fa18", "Sp18", "10.6 GeV"
};

// Build the set of all average-column names we will write per row
static std::vector<std::string> avgColumns_for(const char* base) {
    std::vector<std::string> cols;
    cols.reserve(8);
    for (auto& p : kPeriodsPretty)  cols.emplace_back(std::string(base) + ", " + p);
    for (auto& g : kCombinedPretty) cols.emplace_back(std::string(base) + ", " + g);
    return cols;
}

static const std::vector<std::string> kCols_xB = avgColumns_for("xBavg");
static const std::vector<std::string> kCols_Q2 = avgColumns_for("Q2avg");
static const std::vector<std::string> kCols_t  = avgColumns_for("t_abs_avg");
static const std::vector<std::string> kCols_ph = avgColumns_for("phiavg");

// ---------------------------------------------------------------------
// Minimal CSV utilities (quote-safe split/join)
// ---------------------------------------------------------------------

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    out.reserve(64);
    std::string cur;
    cur.reserve(line.size());
    bool in_quotes = false;
    for (char c : line) {
        if (c == '"') { in_quotes = !in_quotes; continue; }
        if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
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

static inline std::string get_col(const std::vector<std::string>& row,
                                  const std::unordered_map<std::string,int>& idx,
                                  const std::string& name) {
    auto it = idx.find(name);
    if (it == idx.end()) return std::string();
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return std::string();
    return row[j];
}

static inline void put_col(std::vector<std::string>& row,
                           const std::unordered_map<std::string,int>& idx,
                           const std::string& name,
                           const std::string& value) {
    auto it = idx.find(name);
    if (it == idx.end()) return;
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return;
    row[j] = value;
}

static inline double to_double_default(const std::string& s, double dv) {
    if (s.empty()) return dv;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return dv;
    return v;
}

static inline int to_int_default(const std::string& s, int dv) {
    if (s.empty()) return dv;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return dv;
    return (int)v;
}

// ---------------------------------------------------------------------
// Topology and branch binder
// ---------------------------------------------------------------------

enum class Topology { FD_FD, CD_FD, CD_FT, Unknown };

static inline Topology classify_topology(int detector1, int detector2) {
    if (detector1 == 1 && detector2 == 1) return Topology::FD_FD;
    if (detector1 == 2 && detector2 == 1) return Topology::CD_FD;
    if (detector1 == 2 && detector2 == 0) return Topology::CD_FT;
    return Topology::Unknown;
}

static inline const char* topo_to_json_suffix(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
        default: return "FD_FD"; // harmless fallback
    }
}

struct BranchBinder {
    // topology classification inputs
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    // simple-cut branches
    double t1 = 0.0;             bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open = false;
    double pTmiss = 0.0;         bool has_pTmiss = false;

    // kinematic to average
    double x = 0.0;              bool has_x = false;
    double Q2 = 0.0;             bool has_Q2 = false;
    double phi2 = 0.0;           bool has_phi2 = false;
    double Delta_phi = 0.0;      bool has_Delta_phi = false;

    // exclusivity JSON vars (bind if present; skip checks if missing)
    double Emiss2 = 0.0;         bool has_Emiss2 = false;
    double Mx2 = 0.0;            bool has_Mx2 = false;
    double Mx2_1 = 0.0;          bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;          bool has_Mx2_2 = false;
    double theta_gamma_gamma = 0.0; bool has_tgg = false;
    double xF = 0.0;               bool has_xF = false;
    // pi0 only (harmless if absent for DVCS trees)
    double theta_pi0_pi0 = 0.0;    bool has_tpi0 = false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI = [&](const char* n, int* a, bool& f) {
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };
        auto bindD = [&](const char* n, double* a, bool& f) {
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };

        bindI("detector1", &detector1, has_detector1);
        bindI("detector2", &detector2, has_detector2);

        bindD("t1", &t1, has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_open);
        bindD("pTmiss", &pTmiss, has_pTmiss);

        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);
        bindD("phi2", &phi2, has_phi2);
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);

        bindD("Emiss2", &Emiss2, has_Emiss2);
        bindD("Mx2", &Mx2, has_Mx2);
        bindD("Mx2_1", &Mx2_1, has_Mx2_1);
        bindD("Mx2_2", &Mx2_2, has_Mx2_2);
        bindD("theta_gamma_gamma", &theta_gamma_gamma, has_tgg);
        bindD("xF", &xF, has_xF);

        bindD("theta_pi0_pi0", &theta_pi0_pi0, has_tpi0);
    }

    inline double phi() const {
        if (has_phi2) return phi2;
        if (has_Delta_phi) return Delta_phi;
        return std::numeric_limits<double>::quiet_NaN();
    }

    inline bool ready_simple() const {
        return has_detector1 && has_detector2 && has_t1 && has_open && has_pTmiss;
    }
    inline bool ready_binning() const {
        return has_x && has_Q2 && (has_phi2 || has_Delta_phi);
    }
};

// ---------------------------------------------------------------------
// Simple global cuts (always applied)
// ---------------------------------------------------------------------

static inline bool passes_simple_global(const BranchBinder& b) {
    // Same spirit as earlier simple cuts
    if (!b.ready_simple()) return false;

    if (b.open_angle_ep2 <= 5.0) return false;
    if ((-b.t1) > 1.0)          return false;
    if (b.pTmiss > 0.20)        return false;

    return true;
}

// ---------------------------------------------------------------------
// Load and query JSON 3-sigma cuts
// ---------------------------------------------------------------------

// Cache JSON content once
static const json& get_cuts_json() {
    static std::once_flag flag;
    static json J;
    std::call_once(flag, [](){
        std::ifstream fin(kCombinedCutsPath);
        if (!fin.is_open()) {
            std::cerr << "[bin_means] WARNING: Could not open " << kCombinedCutsPath
                      << " — 3-sigma exclusivity cuts will be skipped.\n";
            return;
        }
        try {
            fin >> J;
        } catch (const std::exception& e) {
            std::cerr << "[bin_means] WARNING: Failed to parse combined cuts JSON: "
                      << e.what() << " — skipping 3-sigma.\n";
            J = json(); // empty
        }
    });
    return J;
}

// Normalize period_key into the key style used in the JSON.
// Example: "DVCS_Fa18_inb" -> "DVCS_Fa18_Inb"
//          "DVCS_Sp18_out" -> "DVCS_Sp18_Out"
static std::string normalize_period_key_for_json(const std::string& period_key) {
    std::string s = period_key;
    // Make sure prefix "DVCS_" is capitalized as-is.
    // Replace "_inb" and "_out" with capitalized.
    auto pos = s.rfind("_inb");
    if (pos != std::string::npos) s.replace(pos, 4, "_Inb");
    pos = s.rfind("_out");
    if (pos != std::string::npos) s.replace(pos, 4, "_Out");
    // Common variants: "DVCS_fa18_inb" -> "DVCS_Fa18_Inb"
    // Capitalize Fa/Sp prefix if present
    // Scan for "DVCS_" then capitalize next two letters
    if (s.rfind("DVCS_", 0) == 0 && s.size() >= 7) {
        // After DVCS_, expect something like Fa18 or Sp19
        if (std::islower((unsigned char)s[5])) s[5] = std::toupper((unsigned char)s[5]);
        if (std::islower((unsigned char)s[6])) s[6] = std::toupper((unsigned char)s[6]);
    }
    return s;
}

// Build full JSON key: e.g. "DVCS_Fa18_Inb_CD_FD"
static std::string json_key_for(const std::string& period_key, Topology topo) {
    std::string base = normalize_period_key_for_json(period_key);
    base += "_";
    base += topo_to_json_suffix(topo);
    return base;
}

// Check a single variable against mean +- 3*std if present
static inline bool check_3sigma_one(const json& node, const char* var_name, double value) {
    if (!node.contains(var_name)) return true; // no constraint -> pass
    const auto& v = node[var_name];
    if (!v.contains("mean") || !v.contains("std")) return true;
    double mu = v["mean"].is_number() ? (double)v["mean"] : 0.0;
    double sd = v["std"].is_number()  ? (double)v["std"]  : 0.0;
    if (sd <= 0.0) return true; // ill-defined -> pass
    double lo = mu - 3.0 * sd;
    double hi = mu + 3.0 * sd;
    return (value >= lo && value <= hi);
}

// Apply JSON 3-sigma cuts for DATA side.
// If JSON missing or the specific key missing, this returns true (non-fatal).
static bool passes_3sigma_data(const std::string& period_key, Topology topo, const BranchBinder& b) {
    const json& J = get_cuts_json();
    if (J.is_null() || J.empty()) return true;

    std::string key = json_key_for(period_key, topo);
    auto it = J.find(key);
    if (it == J.end()) {
        // Key not found — try Unknown fallback: accept.
        return true;
    }
    const json& node = *it;
    if (!node.contains("data")) return true;
    const json& D = node["data"];

    // Only check variables that are bound in the tree
    if (b.has_Emiss2 && !check_3sigma_one(D, "Emiss2", b.Emiss2)) return false;
    if (b.has_Mx2   && !check_3sigma_one(D, "Mx2",   b.Mx2))     return false;
    if (b.has_Mx2_1 && !check_3sigma_one(D, "Mx2_1", b.Mx2_1))   return false;
    if (b.has_Mx2_2 && !check_3sigma_one(D, "Mx2_2", b.Mx2_2))   return false;
    if (b.has_pTmiss && !check_3sigma_one(D, "pTmiss", b.pTmiss)) return false;
    if (b.has_tgg   && !check_3sigma_one(D, "theta_gamma_gamma", b.theta_gamma_gamma)) return false;
    if (b.has_xF    && !check_3sigma_one(D, "xF", b.xF))         return false;
    if (b.has_tpi0  && !check_3sigma_one(D, "theta_pi0_pi0", b.theta_pi0_pi0)) return false;

    return true;
}

// Both sets of cuts together (always enforced)
static inline bool passes_all_cuts(const std::string& period_key, const BranchBinder& b) {
    if (!passes_simple_global(b)) return false;

    Topology topo = classify_topology(b.detector1, b.detector2);
    // If topology is unknown, try all three and accept if any passes JSON cut.
    if (topo == Topology::Unknown) {
        if (passes_3sigma_data(period_key, Topology::FD_FD, b)) return true;
        if (passes_3sigma_data(period_key, Topology::CD_FD, b)) return true;
        if (passes_3sigma_data(period_key, Topology::CD_FT, b)) return true;
        return false;
    }
    return passes_3sigma_data(period_key, topo, b);
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

static inline double wrap_to_2pi(double a) {
    if (!std::isfinite(a)) return 0.0;
    double w = std::fmod(a, 2.0 * M_PI);
    if (w < 0.0) w += 2.0 * M_PI;
    return w;
}

struct Accum {
    double sx = 0.0, sQ = 0.0, st = 0.0, sp = 0.0;
    long long n = 0;
    void add(double x, double Q2, double tabs, double phi) {
        sx += x; sQ += Q2; st += tabs; sp += wrap_to_2pi(phi); ++n;
    }
    std::tuple<double,double,double,double> means() const {
        if (n <= 0) return {0.0,0.0,0.0,0.0};
        double inv = 1.0 / (double)n;
        return {sx*inv, sQ*inv, st*inv, sp*inv};
    }
};

// Check if event is inside this CSV row's bin box
static inline bool in_row_bin(double xB, double Q2, double tabs, double phi,
                              double xBmin, double xBmax,
                              double Q2min, double Q2max,
                              double tmin_abs, double tmax_abs,
                              double phimin, double phimax) {
    if (!(xB >= xBmin && xB < xBmax)) return false;
    if (!(Q2 >= Q2min && Q2 < Q2max)) return false;
    if (!(tabs >= tmin_abs && tabs < tmax_abs)) return false;
    // phi in [0, 2pi), phimin/phimax from CSV are assumed in radians
    double wphi = wrap_to_2pi(phi);
    double pmin = wrap_to_2pi(phimin);
    double pmax = wrap_to_2pi(phimax);
    if (pmin <= pmax) {
        return (wphi >= pmin && wphi < pmax);
    } else {
        // Wrapped interval, e.g. [5.5, 0.2)
        return (wphi >= pmin || wphi < pmax);
    }
}

// Period groups for combined averages
static inline bool is_Fa18(const std::string& pk) {
    return pk.find("Fa18") != std::string::npos;
}
static inline bool is_Sp18(const std::string& pk) {
    return pk.find("Sp18") != std::string::npos;
}
static inline bool is_106(const std::string& pk) {
    // 2018 data: Fa18 and Sp18 (10.6 GeV running)
    return is_Fa18(pk) || is_Sp18(pk);
}

// ---------------------------------------------------------------------
// Per-period worker: returns vector<Accum> of size rows.size()
// ---------------------------------------------------------------------

struct WorkerResult {
    std::vector<Accum> per_row; // same size as csv rows
};

static WorkerResult run_period_worker(const std::string& period_key,
                                      TTree* tree,
                                      const std::vector<std::vector<std::string>>& rows,
                                      const std::unordered_map<std::string,int>& Hidx) {
    WorkerResult R;
    R.per_row.resize(rows.size());

    if (!tree) return R;

    BranchBinder b; b.bind(tree);
    if (!b.ready_binning() || !b.ready_simple()) {
        std::cerr << "[bin_means] FATAL: Tree for " << period_key
                  << " missing needed branches for cuts/binning.\n";
        return R;
    }

    const Long64_t nent = tree->GetEntries();
    for (Long64_t i = 0; i < nent; ++i) {
        tree->GetEntry(i);

        if (!passes_all_cuts(period_key, b)) continue;

        double xB = b.x;
        double Q2 = b.Q2;
        double tabs = std::fabs(b.t1);
        double phi = b.phi();
        if (!std::isfinite(xB) || !std::isfinite(Q2) ||
            !std::isfinite(tabs) || !std::isfinite(phi)) continue;

        // Find which CSV row this event belongs to by scanning bins
        // (This is simple and robust; if performance becomes an issue we
        // can add indexing later.)
        for (size_t r = 0; r < rows.size(); ++r) {
            const auto& row = rows[r];
            double xBmin = to_double_default(get_col(row, Hidx, "xBmin"), 0.0);
            double xBmax = to_double_default(get_col(row, Hidx, "xBmax"), 0.0);
            double Q2min = to_double_default(get_col(row, Hidx, "Q2min"), 0.0);
            double Q2max = to_double_default(get_col(row, Hidx, "Q2max"), 0.0);
            double tminA = to_double_default(get_col(row, Hidx, "t_abs_min"), 0.0);
            double tmaxA = to_double_default(get_col(row, Hidx, "t_abs_max"), 0.0);
            double pmin  = to_double_default(get_col(row, Hidx, "phimin"), 0.0);
            double pmax  = to_double_default(get_col(row, Hidx, "phimax"), 0.0);

            if (!in_row_bin(xB, Q2, tabs, phi, xBmin, xBmax, Q2min, Q2max, tminA, tmaxA, pmin, pmax)) continue;

            R.per_row[r].add(xB, Q2, tabs, phi);
        }
    }

    return R;
}

// ---------------------------------------------------------------------
// Main entry
// ---------------------------------------------------------------------

void calculate_bin_means(
    const std::vector<std::string>& dvcs_periods,
    const std::vector<std::string>& /*topologies_unused*/,
    const std::string& /*analysis_type*/,
    const std::vector<Binning>& /*binning_scheme_unused*/,
    const std::string& /*output_json_unused*/,
    const std::map<std::string, TTree*>& dataTrees
) {
    // Load CSV to update
    const std::string csv_inout = "output/csvs/dvcs_pass2_analysis.csv";

    // Read entire CSV
    std::ifstream fin(csv_inout);
    if (!fin.is_open()) {
        std::cerr << "[bin_means] FATAL: Could not open CSV: " << csv_inout << "\n";
        std::exit(EXIT_FAILURE);
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        std::cerr << "[bin_means] FATAL: Empty CSV: " << csv_inout << "\n";
        std::exit(EXIT_FAILURE);
    }
    std::vector<std::string> header = split_csv_line(header_line);
    auto Hidx = build_header_index(header);

    // Read rows
    std::vector<std::vector<std::string>> rows;
    rows.reserve(4096);
    for (std::string line; std::getline(fin, line); ) {
        if (line.empty()) continue;
        rows.emplace_back(split_csv_line(line));
    }
    fin.close();

    // Validate that all needed average columns exist
    auto check_cols = [&](const std::vector<std::string>& cols) {
        for (auto& c : cols) {
            if (!Hidx.count(c)) {
                std::cerr << "[bin_means] FATAL: Missing CSV column: " << c << "\n";
                std::exit(EXIT_FAILURE);
            }
        }
    };
    check_cols(kCols_xB);
    check_cols(kCols_Q2);
    check_cols(kCols_t);
    check_cols(kCols_ph);

    // Launch up to 5 workers (one per period) to accumulate per-row means
    std::vector<std::future<WorkerResult>> futs;
    futs.reserve(dvcs_periods.size());

    for (const auto& pk : dvcs_periods) {
        auto it = dataTrees.find(pk);
        if (it == dataTrees.end() || !it->second) {
            std::cerr << "[bin_means] WARNING: Missing tree for " << pk << " — skipping period.\n";
            continue;
        }
        futs.emplace_back(std::async(std::launch::async, run_period_worker,
                                     pk, it->second, std::cref(rows), std::cref(Hidx)));
    }

    // Collect per-period results; build a map keyed by "pretty" label
    // We assume dvcs_periods are canonical tree keys like "DVCS_Fa18_inb".
    std::unordered_map<std::string, std::vector<Accum>> perPeriod; // pretty -> per-row
    size_t fut_i = 0;
    for (const auto& pk : dvcs_periods) {
        // Map canonical period key to pretty string used in header
        std::string pretty;
        if (pk.find("Fa18_inb") != std::string::npos)      pretty = "Fa18 Inb";
        else if (pk.find("Fa18_out") != std::string::npos) pretty = "Fa18 Out";
        else if (pk.find("Sp19_inb") != std::string::npos) pretty = "Sp19 Inb";
        else if (pk.find("Sp18_inb") != std::string::npos) pretty = "Sp18 Inb";
        else if (pk.find("Sp18_out") != std::string::npos) pretty = "Sp18 Out";
        else {
            // Unrecognized; skip fetch for this pk
            if (fut_i < futs.size()) { (void)futs[fut_i].get(); ++fut_i; }
            continue;
        }

        if (fut_i >= futs.size()) break;
        WorkerResult R = futs[fut_i++].get();
        if (R.per_row.size() != rows.size()) {
            std::cerr << "[bin_means] WARNING: Worker size mismatch for " << pk << "\n";
            R.per_row.resize(rows.size());
        }
        perPeriod[pretty] = std::move(R.per_row);
    }

    // Build combined groups per row by summing corresponding period accumulators
    auto empty_acc_vec = std::vector<Accum>(rows.size());
    auto get_or_empty = [&](const std::string& key) -> const std::vector<Accum>& {
        auto it = perPeriod.find(key);
        if (it == perPeriod.end()) return empty_acc_vec;
        return it->second;
    };

    std::vector<Accum> Fa18(rows.size()), Sp18(rows.size()), G106(rows.size());
    {
        const auto& Fa18Inb = get_or_empty("Fa18 Inb");
        const auto& Fa18Out = get_or_empty("Fa18 Out");
        for (size_t r = 0; r < rows.size(); ++r) {
            Fa18[r].sx = Fa18Inb[r].sx + Fa18Out[r].sx;
            Fa18[r].sQ = Fa18Inb[r].sQ + Fa18Out[r].sQ;
            Fa18[r].st = Fa18Inb[r].st + Fa18Out[r].st;
            Fa18[r].sp = Fa18Inb[r].sp + Fa18Out[r].sp;
            Fa18[r].n  = Fa18Inb[r].n  + Fa18Out[r].n;
        }
    }
    {
        const auto& Sp18Inb = get_or_empty("Sp18 Inb");
        const auto& Sp18Out = get_or_empty("Sp18 Out");
        for (size_t r = 0; r < rows.size(); ++r) {
            Sp18[r].sx = Sp18Inb[r].sx + Sp18Out[r].sx;
            Sp18[r].sQ = Sp18Inb[r].sQ + Sp18Out[r].sQ;
            Sp18[r].st = Sp18Inb[r].st + Sp18Out[r].st;
            Sp18[r].sp = Sp18Inb[r].sp + Sp18Out[r].sp;
            Sp18[r].n  = Sp18Inb[r].n  + Sp18Out[r].n;
        }
    }
    {
        for (size_t r = 0; r < rows.size(); ++r) {
            // 10.6 GeV group = Fa18 + Sp18 (the 2018 running)
            G106[r].sx = Fa18[r].sx + Sp18[r].sx;
            G106[r].sQ = Fa18[r].sQ + Sp18[r].sQ;
            G106[r].st = Fa18[r].st + Sp18[r].st;
            G106[r].sp = Fa18[r].sp + Sp18[r].sp;
            G106[r].n  = Fa18[r].n  + Sp18[r].n;
        }
    }

    // Write averages into the CSV rows
    auto fmt = [](double v) {
        std::ostringstream oss;
        oss.setf(std::ios::fixed); oss << std::setprecision(8) << v;
        return oss.str();
    };

    // Helper to extract mean tuple
    auto means_tuple = [](const Accum& A) {
        return A.means();
    };

    for (size_t r = 0; r < rows.size(); ++r) {
        auto put_all = [&](const std::string& pretty, size_t period_slot) {
            auto itP = perPeriod.find(pretty);
            if (itP == perPeriod.end()) return;
            const Accum& A = itP->second[r];
            auto [mx, mQ, mt, mp] = means_tuple(A);
            put_col(rows[r], Hidx, kCols_xB[period_slot], fmt(mx));
            put_col(rows[r], Hidx, kCols_Q2[period_slot], fmt(mQ));
            put_col(rows[r], Hidx, kCols_t [period_slot], fmt(mt));
            put_col(rows[r], Hidx, kCols_ph[period_slot], fmt(mp));
        };

        // Period slots 0..4
        put_all("Fa18 Inb", 0);
        put_all("Fa18 Out", 1);
        put_all("Sp19 Inb", 2);
        put_all("Sp18 Inb", 3);
        put_all("Sp18 Out", 4);

        // Combined slots 5..7
        {
            auto [mx, mQ, mt, mp] = means_tuple(Fa18[r]);
            put_col(rows[r], Hidx, kCols_xB[5], fmt(mx));
            put_col(rows[r], Hidx, kCols_Q2[5], fmt(mQ));
            put_col(rows[r], Hidx, kCols_t [5], fmt(mt));
            put_col(rows[r], Hidx, kCols_ph[5], fmt(mp));
        }
        {
            auto [mx, mQ, mt, mp] = means_tuple(Sp18[r]);
            put_col(rows[r], Hidx, kCols_xB[6], fmt(mx));
            put_col(rows[r], Hidx, kCols_Q2[6], fmt(mQ));
            put_col(rows[r], Hidx, kCols_t [6], fmt(mt));
            put_col(rows[r], Hidx, kCols_ph[6], fmt(mp));
        }
        {
            auto [mx, mQ, mt, mp] = means_tuple(G106[r]);
            put_col(rows[r], Hidx, kCols_xB[7], fmt(mx));
            put_col(rows[r], Hidx, kCols_Q2[7], fmt(mQ));
            put_col(rows[r], Hidx, kCols_t [7], fmt(mt));
            put_col(rows[r], Hidx, kCols_ph[7], fmt(mp));
        }
    }

    // Overwrite CSV in place (your main.cpp should have already copied a backup)
    std::ofstream fout(csv_inout);
    if (!fout.is_open()) {
        std::cerr << "[bin_means] FATAL: Could not open CSV for writing: " << csv_inout << "\n";
        std::exit(EXIT_FAILURE);
    }
    fout << join_csv_row(header) << "\n";
    for (const auto& row : rows) fout << join_csv_row(row) << "\n";
    fout.close();

    std::cout << "[bin_means] Updated per-row means in " << csv_inout << "\n";
}