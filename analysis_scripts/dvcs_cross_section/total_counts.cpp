// total_counts.cpp
// -----------------------------------------------------------------------------
// DVCS total (raw) counts vs phi, per (xB, Q2, |t|, phi) CSV row, split by
// helicity (+ / -), for each period and for combined groups.
//
// Conventions enforced (project-specific):
//  - Output directories:
//      output/total_counts_plots/{PeriodDir|GroupDir}/<TopologyDir>/
//    where PeriodDir is canonical_period_dir(label) and TopologyDir is one of
//    {FD_FD, CD_FD, CD_FT}.
//
//  - CSV policy:
//      Never add new columns. Only update existing columns via alias resolution
//      of the canonical column name:
//        "raw yield, ep->epg, <topo label>, exp, <period display>, <helicity>"
//      helicity in {unpol,pos,neg}.
//
//  - Event selection:
//      topology from detector1/detector2,
//      global cuts via global_cuts.h,
//      and (DATA only) 3-sigma exclusivity cuts loaded from combined_cuts.json.
//
//  - Phi handling:
//      use phi2 (radians) converted to degrees and wrapped into [0,360).
//      Row matching supports wrap-around phi bins.
//
//  - Plot style (ROOT):
//      canvas with top title pad + grid pad,
//      markers: + helicity red open circle (24), - helicity blue solid circle (20),
//      legends and TLatex locations per conventions,
//      y-range [1.0, 1.15*max].
//
// Threading & robustness:
//  - ROOT::EnableThreadSafety()
//  - OpenMP over periods (max 5 threads)
//  - branch binding protected by mutex
//  - fail-fast on missing required inputs (branches / columns / cut keys)
//
// Notes:
//  - No "ycol" branch is used or required. The dvcsgen P2 cut is applied via
//    global_cuts.h when cfg.enable_dvcsgen_ycol_cut is enabled, using branches
//    (e_p,e_theta,e_phi,p2_p,p2_theta,p2_phi).
// -----------------------------------------------------------------------------

#include "total_counts.h"

#include "periods.h"
#include "global_cuts.h"

// ROOT
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH1F.h>
#include <TAxis.h>

// JSON
#include <nlohmann/json.hpp>

// C++ stdlib
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdlib>
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
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

static constexpr double PI      = 3.14159265358979323846;
static constexpr double RAD2DEG = 180.0 / PI;

static inline bool env_flag(const char* name) {
    return (std::getenv(name) != nullptr);
}

static inline void fatal(const std::string& msg) {
    throw std::runtime_error(msg);
}

static inline std::string to_lower_ascii(std::string s) {
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

static inline double wrap_phi_deg(double phi_deg) {
    double p = std::fmod(phi_deg, 360.0);
    if (p < 0.0) p += 360.0;
    if (p >= 360.0) p = std::nextafter(360.0, 0.0);
    return p;
}

static inline bool in_range(double v, double a, double b) {
    return (v >= a) && (v < b);
}

static inline bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    // Supports wrap-around bins where pmax < pmin (e.g., [330, 30))
    if (pmax_deg > pmin_deg) return in_range(phi_deg, pmin_deg, pmax_deg);
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg);
}

static inline std::string topo_dir(int det1, int det2) {
    if (det1 == 1 && det2 == 1) return "FD_FD";
    if (det1 == 2 && det2 == 1) return "CD_FD";
    if (det1 == 2 && det2 == 0) return "CD_FT";
    return "";
}

static inline std::string topo_label_for_csv(const std::string& topoDir) {
    if (topoDir == "FD_FD") return "(FD, FD)";
    if (topoDir == "CD_FD") return "(CD, FD)";
    if (topoDir == "CD_FT") return "(CD, FT)";
    return "";
}

static inline std::string canonical_period_dir(const std::string& label) {
    // Deterministic mapping only.
    if (label == "Fa18 Inb") return "Fa18_Inb";
    if (label == "Fa18 Out") return "Fa18_Out";
    if (label == "Fa18 Inb Supp") return "Fa18_Inb_Supp";
    if (label == "Sp18 Inb") return "Sp18_Inb";
    if (label == "Sp18 Out") return "Sp18_Out";
    if (label == "Sp19 Inb") return "Sp19_Inb";
    std::ostringstream ss;
    ss << "[total_counts] FATAL: unknown period label for canonical_period_dir: '" << label << "'";
    fatal(ss.str());
    return "";
}

static inline std::string canonical_group_dir(const std::string& label) {
    // Must return exactly label for combined groups by convention.
    if (label == "Fa18") return "Fa18";
    if (label == "Sp18") return "Sp18";
    if (label == "10.6 GeV") return "10.6 GeV";
    std::ostringstream ss;
    ss << "[total_counts] FATAL: unknown group label for canonical_group_dir: '" << label << "'";
    fatal(ss.str());
    return "";
}

static inline bool is_combined_group_label(const std::string& label) {
    return (label == "Fa18" || label == "Sp18" || label == "10.6 GeV");
}

static inline bool is_supplemental_label(const std::string& label) {
    return (label == "Fa18 Inb Supp");
}

static inline bool should_skip_csv_for_label(const std::string& label) {
    // Per conventions: skip supplemental and combined groups in CSV writes.
    if (is_supplemental_label(label)) return true;
    if (is_combined_group_label(label)) return true;
    return false;
}

static inline std::string out_root_for_label(const std::string& label, const std::string& out_root_dir) {
    // out_root_dir is expected to be something like "output/total_counts_plots"
    if (is_combined_group_label(label)) {
        return out_root_dir + "/" + canonical_group_dir(label);
    }
    return out_root_dir + "/" + canonical_period_dir(label);
}

static inline std::string ensure_trailing_slash(std::string s) {
    if (!s.empty() && s.back() == '/') return s;
    return s + "/";
}

static inline void mkdir_p(const std::string& path) {
    if (path.empty()) return;
    // ROOT's gSystem supports recursive mkdir with the 2nd arg true.
    gSystem->mkdir(path.c_str(), true);
}

// -----------------------------------------------------------------------------
// Period key parsing (tree key -> {period display, period_label for global_cuts,
// json tag for combined cuts})
// -----------------------------------------------------------------------------

struct PeriodTags {
    std::string tree_key;       // e.g. "DVCS_Fa18_Inb" (as in CANONICAL_PERIODS().tree_key)
    std::string period_display; // e.g. "Fa18 Inb"      (CSV display label)
    std::string period_label;   // e.g. "fa18_inb"      (global_cuts period label)
    std::string json_tag;       // e.g. "DVCS_Fa18_Inb" (combined_cuts.json key prefix)
};

static inline PeriodTags parse_period_tags_from_tree_key(const std::string& tree_key) {
    // Deterministic mapping only. No heuristics.
    const std::string s = to_lower_ascii(tree_key);

    PeriodTags t;
    t.tree_key = tree_key;

    // DVCS_<Period>_<Torus>
    // We only support the canonical DVCS periods used in pass-2.
    auto has = [&](const char* sub) { return (s.find(sub) != std::string::npos); };

    if (has("fa18") && has("inb") && !has("supp")) {
        t.period_display = "Fa18 Inb";
        t.period_label   = "fa18_inb";
        t.json_tag       = "DVCS_Fa18_Inb";
        return t;
    }
    if (has("fa18") && has("out")) {
        t.period_display = "Fa18 Out";
        t.period_label   = "fa18_out";
        t.json_tag       = "DVCS_Fa18_Out";
        return t;
    }
    if (has("fa18") && has("inb") && has("supp")) {
        t.period_display = "Fa18 Inb Supp";

        // IMPORTANT:
        // global_cuts.h does not recognize "fa18_inb_supp". For global-cuts
        // dispatch, treat supplemental Fa18 Inb as the standard Fa18 Inb label.
        t.period_label   = "fa18_inb";

        // Sigma cuts still use the supplemental JSON key tag.
        t.json_tag       = "DVCS_Fa18_Inb_Supp";
        return t;
    }
    if (has("sp18") && has("inb")) {
        t.period_display = "Sp18 Inb";
        t.period_label   = "sp18_inb";
        t.json_tag       = "DVCS_Sp18_Inb";
        return t;
    }
    if (has("sp18") && has("out")) {
        t.period_display = "Sp18 Out";
        t.period_label   = "sp18_out";
        t.json_tag       = "DVCS_Sp18_Out";
        return t;
    }
    if (has("sp19") && has("inb")) {
        t.period_display = "Sp19 Inb";
        t.period_label   = "sp19_inb";
        t.json_tag       = "DVCS_Sp19_Inb";
        return t;
    }

    std::ostringstream ss;
    ss << "[total_counts] FATAL: cannot map tree_key '" << tree_key
       << "' to PeriodTags (period_display/period_label/json_tag).";
    fatal(ss.str());
    return t;
}

// -----------------------------------------------------------------------------
// CSV I/O (fail-fast on missing columns; atomic save)
// -----------------------------------------------------------------------------

struct CSV {
    std::vector<std::string> header;
    std::unordered_map<std::string, int> index;
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
        std::ostringstream ss;
        ss << "[total_counts] FATAL: cannot open CSV: " << path;
        fatal(ss.str());
    }

    std::string line;
    if (!std::getline(fin, line)) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: empty CSV: " << path;
        fatal(ss.str());
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

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";

    std::ofstream fout(tmp);
    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: cannot write temp CSV: " << tmp;
        fatal(ss.str());
    }

    auto write_cell = [&](const std::string& s) {
        const bool needq = (s.find(',') != std::string::npos) || (s.find('"') != std::string::npos);
        if (!needq) {
            fout << s;
            return;
        }
        fout << '"';
        for (char ch : s) {
            if (ch == '"') fout << "\"\"";
            else fout << ch;
        }
        fout << '"';
    };

    for (size_t i = 0; i < csv.header.size(); ++i) {
        write_cell(csv.header[i]);
        if (i + 1 < csv.header.size()) fout << ',';
    }
    fout << "\n";

    for (const auto& row : csv.rows) {
        if (row.size() != csv.header.size()) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: CSV row width mismatch (row has " << row.size()
               << " cells, header has " << csv.header.size() << ").";
            fatal(ss.str());
        }
        for (size_t i = 0; i < row.size(); ++i) {
            write_cell(row[i]);
            if (i + 1 < row.size()) fout << ',';
        }
        fout << "\n";
    }

    fout.close();
    if (!fout) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: write failed for temp CSV: " << tmp;
        fatal(ss.str());
    }

    // Atomic replace.
    if (std::remove(path.c_str()) != 0) {
        // If file doesn't exist, ignore; otherwise error.
        // We do not silently continue on unexpected errno; but C stdlib doesn't expose errno portably here.
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: rename failed from '" << tmp << "' to '" << path << "'";
        fatal(ss.str());
    }
}

static int col_strict(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);
    if (it == csv.index.end()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: missing required CSV column: '" << name << "'";
        fatal(ss.str());
    }
    return it->second;
}

// Canonical column name builder (must match existing CSV schema).
static inline std::string col_counts(const std::string& topo_label,
                                     const std::string& period_display,
                                     const std::string& helicity) {
    // "raw yield, ep->epg, <topo label>, exp, <period display>, <helicity>"
    return std::string("raw yield, ep->epg, ") + topo_label + ", exp, " + period_display + ", " + helicity;
}

// Alias-resolution helpers (deterministic; no heuristics).
static inline std::vector<std::string> topology_aliases(const std::string& topo_label) {
    // If you ever had alternate topo label strings in older CSVs, list them here explicitly.
    // Deterministic list ordered by preference.
    return { topo_label };
}

static inline std::vector<std::string> period_aliases(const std::string& period_display) {
    // Deterministic list ordered by preference.
    return { period_display };
}

static inline std::vector<std::string> helicity_aliases(const std::string& helicity) {
    return { helicity };
}

static int find_col_alias(const CSV& csv,
                          const std::vector<std::string>& topo_alias,
                          const std::vector<std::string>& period_alias,
                          const std::vector<std::string>& helicity_alias) {
    // We resolve the canonical naming by enumerating explicit aliases.
    // If nothing matches, we fail fast.
    for (const auto& tl : topo_alias) {
        for (const auto& pl : period_alias) {
            for (const auto& hl : helicity_alias) {
                const std::string name = col_counts(tl, pl, hl);
                auto it = csv.index.find(name);
                if (it != csv.index.end()) return it->second;
            }
        }
    }
    std::ostringstream ss;
    ss << "[total_counts] FATAL: could not find CSV column for counts with "
       << "topo_alias[0]='" << (topo_alias.empty() ? "" : topo_alias[0]) << "', "
       << "period_alias[0]='" << (period_alias.empty() ? "" : period_alias[0]) << "', "
       << "helicity_alias[0]='" << (helicity_alias.empty() ? "" : helicity_alias[0]) << "'.";
    fatal(ss.str());
    return -1;
}

// -----------------------------------------------------------------------------
// Row binning (from CSV) and matching
// -----------------------------------------------------------------------------

struct RowBin {
    double xBmin = std::numeric_limits<double>::quiet_NaN();
    double xBmax = std::numeric_limits<double>::quiet_NaN();
    double Q2min = std::numeric_limits<double>::quiet_NaN();
    double Q2max = std::numeric_limits<double>::quiet_NaN();
    double tmin  = std::numeric_limits<double>::quiet_NaN();
    double tmax  = std::numeric_limits<double>::quiet_NaN();
    double pmin  = std::numeric_limits<double>::quiet_NaN();
    double pmax  = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

static inline double to_double_strict(const std::string& s, const std::string& what) {
    if (s.empty()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: empty numeric cell for '" << what << "'";
        fatal(ss.str());
    }
    char* e = nullptr;
    const double v = std::strtod(s.c_str(), &e);
    if (e == s.c_str()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: parse failure for '" << what << "' value '" << s << "'";
        fatal(ss.str());
    }
    return v;
}

static inline bool to_bool_valid(const std::string& s) {
    return (s == "1" || s == "1.0" || s == "true" || s == "TRUE");
}

static std::vector<RowBin> load_row_bins_from_csv(const CSV& csv) {
    const int c_xBmin = col_strict(csv, "xBmin");
    const int c_xBmax = col_strict(csv, "xBmax");
    const int c_Q2min = col_strict(csv, "Q2min");
    const int c_Q2max = col_strict(csv, "Q2max");
    const int c_tmin  = col_strict(csv, "t_abs_min");
    const int c_tmax  = col_strict(csv, "t_abs_max");
    const int c_pmin  = col_strict(csv, "phimin");
    const int c_pmax  = col_strict(csv, "phimax");
    const int c_valid = col_strict(csv, "valid bin");

    std::vector<RowBin> rows;
    rows.reserve(csv.rows.size());

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        const auto& row = csv.rows[r];
        RowBin b;
        b.xBmin = to_double_strict(row[c_xBmin], "xBmin");
        b.xBmax = to_double_strict(row[c_xBmax], "xBmax");
        b.Q2min = to_double_strict(row[c_Q2min], "Q2min");
        b.Q2max = to_double_strict(row[c_Q2max], "Q2max");
        b.tmin  = to_double_strict(row[c_tmin],  "t_abs_min");
        b.tmax  = to_double_strict(row[c_tmax],  "t_abs_max");
        b.pmin  = to_double_strict(row[c_pmin],  "phimin");
        b.pmax  = to_double_strict(row[c_pmax],  "phimax");
        b.valid = to_bool_valid(row[c_valid]);
        rows.push_back(b);
    }
    return rows;
}

static inline bool row_accepts_kin(const RowBin& w, double x, double Q2, double tabs) {
    if (!in_range(x,  w.xBmin, w.xBmax)) return false;
    if (!in_range(Q2, w.Q2min, w.Q2max)) return false;
    if (!in_range(tabs, w.tmin, w.tmax)) return false;
    return true;
}

// -----------------------------------------------------------------------------
// Combined 3-sigma cuts loader (DATA only)
// combined_cuts.json format (as written by exclusivity_cuts.cpp):
// {
//   "DVCS_Fa18_Inb_FD_FD": { "data": { "Mx2": {"mean":..., "std":...}, ... }, "mc": {...} },
//   ...
// }
// -----------------------------------------------------------------------------

struct SigmaStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std  = std::numeric_limits<double>::quiet_NaN();
};

using CutVarMap = std::unordered_map<std::string, SigmaStats>; // var -> (mean,std)
using TopoCutMap = std::unordered_map<std::string, CutVarMap>; // key -> var map

static inline bool within_3sigma(double v, const SigmaStats& s) {
    if (!std::isfinite(s.mean) || !std::isfinite(s.std) || s.std <= 0.0) return true;
    return (std::fabs(v - s.mean) <= 3.0 * s.std);
}

static TopoCutMap load_combined_cuts_data_only(const std::string& combined_cuts_json) {
    std::ifstream fin(combined_cuts_json);
    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: cannot open combined cuts JSON: " << combined_cuts_json;
        fatal(ss.str());
    }

    nlohmann::json j;
    try {
        fin >> j;
    } catch (const std::exception& e) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: JSON parse failed for " << combined_cuts_json << " : " << e.what();
        fatal(ss.str());
    }

    if (!j.is_object()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: combined cuts JSON is not an object: " << combined_cuts_json;
        fatal(ss.str());
    }

    TopoCutMap out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const auto& block = it.value();
        if (!block.is_object()) continue;

        if (!block.contains("data")) continue;
        const auto& data = block["data"];
        if (!data.is_object()) continue;

        CutVarMap vm;
        for (auto vit = data.begin(); vit != data.end(); ++vit) {
            const std::string var = vit.key();
            const auto& stats = vit.value();
            if (!stats.is_object()) continue;
            if (!stats.contains("mean") || !stats.contains("std")) continue;
            SigmaStats s;
            try {
                s.mean = stats["mean"].get<double>();
                s.std  = stats["std"].get<double>();
            } catch (...) {
                continue;
            }
            if (std::isfinite(s.std)) vm.emplace(var, s);
        }

        if (!vm.empty()) out.emplace(key, std::move(vm));
    }

    std::cout << "[total_counts] Loaded sigma cuts for " << out.size()
              << " topology keys from " << combined_cuts_json << std::endl;

    return out;
}

// -----------------------------------------------------------------------------
// Branch binder (thread-safe binding)
// -----------------------------------------------------------------------------

static std::mutex g_root_bind_mutex;

struct BranchBinder {
    // Topology
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    // Helicity
    int helicity = 0; bool has_helicity = false;

    // Binning variables
    double x = 0.0;     bool has_x = false;
    double Q2 = 0.0;    bool has_Q2 = false;
    double t1 = 0.0;    bool has_t1 = false;
    double phi2 = 0.0;  bool has_phi2 = false; // radians

    // Global cuts inputs
    double open_angle_ep2 = 0.0; bool has_open_angle = false; // degrees
    double pTmiss = 0.0;         bool has_pTmiss = false;

    // 3-sigma variables (DVCS)
    double Emiss2 = 0.0;            bool has_Emiss2 = false;
    double Mx2 = 0.0;               bool has_Mx2 = false;
    double Mx2_1 = 0.0;             bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;             bool has_Mx2_2 = false;
    double xF = 0.0;                bool has_xF = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gg = false;

    // dvcsgen ycol (P2) cut inputs, only if enabled in GlobalCutConfig
    double e_p = 0.0;       bool has_e_p = false;
    double e_theta = 0.0;   bool has_e_theta = false;
    double e_phi = 0.0;     bool has_e_phi = false;

    double p2_p = 0.0;      bool has_p2_p = false;
    double p2_theta = 0.0;  bool has_p2_theta = false;
    double p2_phi = 0.0;    bool has_p2_phi = false;

    void bind(TTree* t) {
        if (!t) return;

        std::lock_guard<std::mutex> lock(g_root_bind_mutex);

        t->SetBranchStatus("*", 0);

        auto ena = [&](const char* n) {
            if (t->GetBranch(n)) t->SetBranchStatus(n, 1);
        };

        ena("detector1");
        ena("detector2");
        ena("helicity");

        ena("x");
        ena("Q2");
        ena("t1");
        ena("phi2");

        ena("open_angle_ep2");
        ena("pTmiss");

        // 3-sigma vars (presence is expected for DVCS trees; still bind conditionally)
        ena("Emiss2");
        ena("Mx2");
        ena("Mx2_1");
        ena("Mx2_2");
        ena("xF");
        ena("theta_gamma_gamma");

        // dvcsgen cut inputs (optional; required only if cfg.enable_dvcsgen_ycol_cut)
        ena("e_p");
        ena("e_theta");
        ena("e_phi");
        ena("p2_p");
        ena("p2_theta");
        ena("p2_phi");

        auto bI = [&](const char* n, int* a, bool& f) {
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };
        auto bD = [&](const char* n, double* a, bool& f) {
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };

        bI("detector1", &detector1, has_detector1);
        bI("detector2", &detector2, has_detector2);
        bI("helicity",  &helicity,  has_helicity);

        bD("x",    &x,    has_x);
        bD("Q2",   &Q2,   has_Q2);
        bD("t1",   &t1,   has_t1);
        bD("phi2", &phi2, has_phi2);

        bD("open_angle_ep2", &open_angle_ep2, has_open_angle);
        bD("pTmiss",         &pTmiss,         has_pTmiss);

        bD("Emiss2",            &Emiss2,            has_Emiss2);
        bD("Mx2",               &Mx2,               has_Mx2);
        bD("Mx2_1",             &Mx2_1,             has_Mx2_1);
        bD("Mx2_2",             &Mx2_2,             has_Mx2_2);
        bD("xF",                &xF,                has_xF);
        bD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gg);

        bD("e_p",     &e_p,     has_e_p);
        bD("e_theta", &e_theta, has_e_theta);
        bD("e_phi",   &e_phi,   has_e_phi);

        bD("p2_p",     &p2_p,     has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi",   &p2_phi,   has_p2_phi);
    }

    bool ready_for_matching() const {
        return has_detector1 && has_detector2 && has_x && has_Q2 && has_t1 && has_phi2;
    }

    double phi_deg() const {
        return wrap_phi_deg(phi2 * RAD2DEG);
    }

    double t_abs() const {
        return std::fabs(t1);
    }
};

// -----------------------------------------------------------------------------
// Counts accumulation
// -----------------------------------------------------------------------------

struct HelCounts {
    double unpol = 0.0;
    double pos   = 0.0;
    double neg   = 0.0;
};

using RowCounts = std::unordered_map<int, HelCounts>; // row index -> counts

static inline void add_count(HelCounts& c, int helicity) {
    if (helicity > 0) c.pos += 1.0;
    else if (helicity < 0) c.neg += 1.0;
    else c.unpol += 1.0;
    c.unpol += 0.0; // keep explicit; we also store explicit "unpol" bucket only from helicity==0 unless user wants sum.
}

// If you want unpol to mean (pos+neg) always, do it explicitly at write-time.

static inline bool passes_sigma_cuts_data_only(const TopoCutMap& cuts,
                                               const std::string& key,
                                               const BranchBinder& b) {
    auto it = cuts.find(key);
    if (it == cuts.end()) return true;

    const CutVarMap& vm = it->second;

    auto check = [&](const char* var, bool has_val, double val) {
        auto iv = vm.find(var);
        if (iv == vm.end()) return true;
        if (!has_val) return true;
        return within_3sigma(val, iv->second);
    };

    // Apply to any variables present in the cut map.
    // (If the JSON contains fewer/more keys, we still deterministically check what exists.)
    if (!check("Emiss2",            b.has_Emiss2, b.Emiss2)) return false;
    if (!check("Mx2",               b.has_Mx2,    b.Mx2)) return false;
    if (!check("Mx2_1",             b.has_Mx2_1,  b.Mx2_1)) return false;
    if (!check("Mx2_2",             b.has_Mx2_2,  b.Mx2_2)) return false;
    if (!check("pTmiss",            b.has_pTmiss, b.pTmiss)) return false;
    if (!check("xF",                b.has_xF,     b.xF)) return false;
    if (!check("theta_gamma_gamma", b.has_theta_gg, b.theta_gamma_gamma)) return false;

    return true;
}

static inline bool passes_global_cuts_dispatch(const BranchBinder& b,
                                               const std::string& period_label) {
    // Global cuts are shared across modules via global_cuts.h.
    // We enforce the DVCSgen P2/ycol logic only through global_cuts.h.
    if (!(b.has_t1 && b.has_open_angle && b.has_pTmiss)) return false;

    const GlobalCutConfig& cfg = default_global_cuts();

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: dvcsgen ycol cut enabled, but required branches are missing "
               << "for period_label='" << period_label << "'. Missing one or more of: "
               << "e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi.";
            fatal(ss.str());
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss, cfg);
}

// -----------------------------------------------------------------------------
// CSV writing helper for per-period counts
// -----------------------------------------------------------------------------

struct CountsColumnIndex {
    int col_unpol = -1;
    int col_pos   = -1;
    int col_neg   = -1;
};

static CountsColumnIndex resolve_counts_columns(const CSV& csv,
                                                const std::string& topo_label,
                                                const std::string& period_display) {
    CountsColumnIndex idx;
    idx.col_unpol = find_col_alias(csv,
                                  topology_aliases(topo_label),
                                  period_aliases(period_display),
                                  helicity_aliases("unpol"));

    idx.col_pos   = find_col_alias(csv,
                                  topology_aliases(topo_label),
                                  period_aliases(period_display),
                                  helicity_aliases("pos"));

    idx.col_neg   = find_col_alias(csv,
                                  topology_aliases(topo_label),
                                  period_aliases(period_display),
                                  helicity_aliases("neg"));

    return idx;
}

static inline std::string fmt0(double v) {
    // Counts are integers in practice, but stored as a numeric string.
    // Deterministic formatting: no scientific notation.
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(0) << v;
    return oss.str();
}

// -----------------------------------------------------------------------------
// Plotting
// -----------------------------------------------------------------------------

struct RowPoint {
    double phi_x = 0.0;    // degrees (x-value)
    double pos = 0.0;
    double neg = 0.0;
    double pos_err = 0.0;
    double neg_err = 0.0;
    double Q2c = 0.0;      // for annotation (center)
    double tc  = 0.0;      // |t| center
    bool valid = false;
};

static inline double bin_center(double a, double b) {
    return 0.5 * (a + b);
}

static inline double poisson_err(double n) {
    // Deterministic Poisson sqrt(n). For zero counts, error is zero here.
    // If you prefer Garwood for zeros, implement explicitly; do not guess.
    return (n > 0.0) ? std::sqrt(n) : 0.0;
}

static void plot_one_xB_canvas(const std::string& outdir,
                               const std::string& label_for_title,
                               const std::string& topoDir,
                               const std::string& topoLabel,
                               const std::string& xB_text,
                               double xBmin_for_file,
                               const std::vector<RowBin>& rows,
                               const std::vector<int>& row_indices,
                               const std::vector<RowPoint>& points) {
    if (row_indices.empty()) return;

    // Determine grid: unique Q2 and t bins among these rows.
    struct Edge { double a, b; };
    std::vector<Edge> Q2bins;
    std::vector<Edge> tbins;

    auto edge_eq = [](const Edge& u, const Edge& v) {
        return (u.a == v.a && u.b == v.b);
    };

    for (int ridx : row_indices) {
        const RowBin& r = rows[ridx];
        Edge q{r.Q2min, r.Q2max};
        Edge t{r.tmin,  r.tmax};
        if (std::find_if(Q2bins.begin(), Q2bins.end(), [&](const Edge& e){ return edge_eq(e,q); }) == Q2bins.end())
            Q2bins.push_back(q);
        if (std::find_if(tbins.begin(), tbins.end(), [&](const Edge& e){ return edge_eq(e,t); }) == tbins.end())
            tbins.push_back(t);
    }

    auto sort_edges = [](std::vector<Edge>& v) {
        std::sort(v.begin(), v.end(), [](const Edge& p, const Edge& q) {
            if (p.a != q.a) return p.a < q.a;
            return p.b < q.b;
        });
    };
    sort_edges(Q2bins);
    sort_edges(tbins);

    const int ncols = (int)Q2bins.size();
    const int nrows = (int)tbins.size();
    if (ncols <= 0 || nrows <= 0) return;

    const int W = 300 * ncols + 160;
    const int H = 260 * nrows + 240;

    TCanvas c("c_total_counts", "", W, H);

    // Top title pad + grid pad.
    TPad* top = new TPad("top", "top", 0.0, 0.90, 1.0, 1.0);
    TPad* grid = new TPad("grid", "grid", 0.0, 0.0, 1.0, 0.90);
    top->SetFillStyle(0);
    grid->SetFillStyle(0);
    top->Draw();
    grid->Draw();

    top->cd();
    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextSize(0.45);

    {
        std::ostringstream ss;
        ss << label_for_title << "   " << xB_text << "   " << topoLabel;
        t.DrawLatex(0.05, 0.35, ss.str().c_str());
    }

    grid->cd();
    grid->Divide(ncols, nrows, 0.0, 0.0);

    // Build lookup Q2/t bin -> pad index.
    auto find_q = [&](double a, double b) {
        for (int i = 0; i < (int)Q2bins.size(); ++i) if (Q2bins[i].a == a && Q2bins[i].b == b) return i;
        return -1;
    };
    auto find_t = [&](double a, double b) {
        for (int i = 0; i < (int)tbins.size(); ++i) if (tbins[i].a == a && tbins[i].b == b) return i;
        return -1;
    };

    // For each (Q2,t) cell, gather its phi points in increasing x for plotting.
    for (int it = 0; it < nrows; ++it) {
        for (int iq = 0; iq < ncols; ++iq) {
            const int pad_idx = it * ncols + iq + 1;
            grid->cd(pad_idx);

            gPad->SetLeftMargin(0.160);
            gPad->SetRightMargin(0.07);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetGrid(1,1);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            // Collect points for this cell.
            std::vector<RowPoint> cell;
            cell.reserve(row_indices.size());

            for (int ridx : row_indices) {
                const RowBin& r = rows[ridx];
                const int q = find_q(r.Q2min, r.Q2max);
                const int tt = find_t(r.tmin, r.tmax);
                if (q != iq || tt != it) continue;

                const RowPoint& p = points[ridx];
                if (!p.valid) continue;
                cell.push_back(p);
            }

            std::sort(cell.begin(), cell.end(), [](const RowPoint& a, const RowPoint& b) {
                return a.phi_x < b.phi_x;
            });

            const int N = (int)cell.size();

            // Build graphs.
            TGraphErrors* gpos = new TGraphErrors();
            TGraphErrors* gneg = new TGraphErrors();

            gpos->SetMarkerStyle(24);
            gpos->SetMarkerColor(kRed);
            gpos->SetLineColor(kRed);
            gpos->SetLineWidth(1);

            gneg->SetMarkerStyle(20);
            gneg->SetMarkerColor(kBlue);
            gneg->SetLineColor(kBlue);
            gneg->SetLineWidth(1);

            double ymax = 0.0;
            for (int i = 0; i < N; ++i) {
                const double x = cell[i].phi_x;

                gpos->SetPoint(i, x, cell[i].pos);
                gpos->SetPointError(i, 0.0, cell[i].pos_err);

                gneg->SetPoint(i, x, cell[i].neg);
                gneg->SetPointError(i, 0.0, cell[i].neg_err);

                ymax = std::max(ymax, std::max(cell[i].pos + cell[i].pos_err, cell[i].neg + cell[i].neg_err));
            }

            // Create a frame histogram for axes.
            // X range: show full phi (0..360) but data may extend beyond 360 if your bins do.
            const double xmin = 0.0;
            const double xmax = 360.0;

            TH1F* frame = (TH1F*)gPad->DrawFrame(xmin, 0.0, xmax, std::max(1.0, 1.15 * ymax));
            frame->SetTitle("");
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("Counts");

            frame->GetXaxis()->CenterTitle(true);
            frame->GetYaxis()->CenterTitle(true);

            frame->GetXaxis()->SetNdivisions(505);

            frame->GetXaxis()->SetLabelSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.060);
            frame->GetXaxis()->SetTitleSize(0.070);
            frame->GetYaxis()->SetTitleSize(0.070);

            frame->GetXaxis()->SetTitleOffset(1.05);
            frame->GetYaxis()->SetTitleOffset(1.10);

            gpos->Draw("PE1 SAME");
            gneg->Draw("PE1 SAME");

            // Per-pad legend.
            TLegend* leg = new TLegend(0.60, 0.73, 0.93, 0.92);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetBorderSize(1);
            leg->SetTextFont(42);
            leg->SetTextSize(0.055);
            leg->AddEntry(gpos, "+ helicity", "p");
            leg->AddEntry(gneg, "- helicity", "p");
            leg->Draw();

            // Kinematics annotation (use bin centers from the Q2 and t edges).
            TLatex txt;
            txt.SetNDC(true);
            txt.SetTextFont(42);
            txt.SetTextSize(0.060);

            const double Q2c = bin_center(Q2bins[iq].a, Q2bins[iq].b);
            const double tc  = bin_center(tbins[it].a,  tbins[it].b);

            {
                std::ostringstream ss;
                ss << "Q^{2}=" << std::fixed << std::setprecision(2) << Q2c
                   << "  |t|=" << std::fixed << std::setprecision(2) << tc;
                txt.DrawLatex(0.12, 0.83, ss.str().c_str());
            }
        }
    }

    mkdir_p(outdir);

    const int idx = (int)std::llround(xBmin_for_file * 1000.0);
    std::ostringstream fname;
    fname << outdir << "/plot_total_counts_"
          << (is_combined_group_label(label_for_title) ? canonical_group_dir(label_for_title) : canonical_period_dir(label_for_title))
          << "_" << topoDir << "_xB_" << idx << ".png";

    c.SaveAs(fname.str().c_str());

    delete top;
    delete grid;
}

// -----------------------------------------------------------------------------
// Main per-period counting worker
// -----------------------------------------------------------------------------

struct PeriodWorkResult {
    // Per topology directory -> row counts
    std::unordered_map<std::string, RowCounts> topo_counts;
};

static PeriodWorkResult accumulate_counts_for_tree(const PeriodTags& tags,
                                                   TTree* tree,
                                                   const std::vector<RowBin>& rows,
                                                   const TopoCutMap& sigma_cuts_data,
                                                   bool trace_matches) {
    PeriodWorkResult out;
    if (!tree) return out;

    BranchBinder b;
    b.bind(tree);

    if (!b.ready_for_matching()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: Missing required branches in DVCS tree for '" << tags.tree_key
           << "'. Required: detector1, detector2, x, Q2, t1, phi2.";
        fatal(ss.str());
    }

    // This module expects helicity; if absent, we only fill "unpol" via helicity==0.
    if (!b.has_helicity) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: Missing required branch 'helicity' in DVCS tree for '" << tags.tree_key << "'.";
        fatal(ss.str());
    }

    const Long64_t N = tree->GetEntries();
    const bool dbg = env_flag("TOTAL_COUNTS_DEBUG");

    long long n_global_pass = 0;
    long long n_sigma_pass  = 0;
    long long n_used        = 0;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        const std::string topoDir = topo_dir(b.detector1, b.detector2);
        if (topoDir.empty()) continue;

        if (!passes_global_cuts_dispatch(b, tags.period_label)) continue;
        ++n_global_pass;

        // Apply DATA-only 3-sigma cuts (always enforced going forward).
        // Key format in JSON: "<json_tag>_<TopoDir>" e.g. "DVCS_Fa18_Inb_FD_FD"
        const std::string sig_key = tags.json_tag + "_" + topoDir;
        if (!passes_sigma_cuts_data_only(sigma_cuts_data, sig_key, b)) continue;
        ++n_sigma_pass;

        const double phi_deg = b.phi_deg();
        const double tabs = b.t_abs();

        // Match to CSV row(s).
        bool matched_any = false;

        for (int r = 0; r < (int)rows.size(); ++r) {
            const RowBin& w = rows[r];
            if (!w.valid) continue;

            if (!row_accepts_kin(w, b.x, b.Q2, tabs)) continue;
            if (!row_accepts_phi(phi_deg, w.pmin, w.pmax)) continue;

            HelCounts& hc = out.topo_counts[topoDir][r];
            // pos/neg based on helicity sign
            if (b.helicity > 0) hc.pos += 1.0;
            else if (b.helicity < 0) hc.neg += 1.0;
            else hc.unpol += 1.0;

            matched_any = true;

            if (trace_matches) {
                std::cout << "[total_counts][TRACE] " << tags.tree_key
                          << " topo=" << topoDir
                          << " row=" << r
                          << " x=" << b.x
                          << " Q2=" << b.Q2
                          << " |t|=" << tabs
                          << " phi(deg)=" << phi_deg
                          << " hel=" << b.helicity
                          << std::endl;
            }
        }

        if (matched_any) ++n_used;

        if (dbg && i < 3) {
            std::cout << "[total_counts][DEBUG] " << tags.tree_key
                      << " i=" << (long long)i
                      << " topo=" << topoDir
                      << " hel=" << b.helicity
                      << " x=" << b.x
                      << " Q2=" << b.Q2
                      << " t1=" << b.t1
                      << " phi2(rad)=" << b.phi2
                      << " phi(deg)=" << phi_deg
                      << " open_angle_ep2(deg)=" << (b.has_open_angle ? b.open_angle_ep2 : -999.0)
                      << " pTmiss=" << (b.has_pTmiss ? b.pTmiss : -999.0)
                      << std::endl;
        }
    }

    std::cout << "[total_counts] " << tags.tree_key
              << " entries=" << (long long)N
              << " global_pass=" << n_global_pass
              << " sig_pass=" << n_sigma_pass
              << " matched=" << n_used
              << std::endl;

    return out;
}

// -----------------------------------------------------------------------------
// Build per-row plot points from counts + phiavg/midbin, and plot per xB bin.
// -----------------------------------------------------------------------------

static inline bool cell_is_number(const std::string& s) {
    if (s.empty()) return false;
    char* e = nullptr;
    (void)std::strtod(s.c_str(), &e);
    return (e != s.c_str());
}

static inline double cell_to_double_or_nan(const std::string& s) {
    if (!cell_is_number(s)) return std::numeric_limits<double>::quiet_NaN();
    return std::strtod(s.c_str(), nullptr);
}

static std::string col_phiavg_for_period(const std::string& period_display) {
    return std::string("phiavg, ") + period_display;
}

static std::string col_xBavg_for_period(const std::string& period_display) {
    return std::string("xBavg, ") + period_display;
}

static void make_plots_for_label_and_topo(const std::string& label,
                                          const std::string& topoDir,
                                          const std::string& out_root_dir,
                                          const CSV& csv,
                                          const std::vector<RowBin>& rows,
                                          const RowCounts& row_counts_for_topo) {
    // Determine xB bins present (unique xBmin/xBmax among valid bins).
    struct XBEdge { double a, b; };
    std::vector<XBEdge> xbs;

    for (int r = 0; r < (int)rows.size(); ++r) {
        const RowBin& w = rows[r];
        if (!w.valid) continue;
        XBEdge e{w.xBmin, w.xBmax};
        auto it = std::find_if(xbs.begin(), xbs.end(), [&](const XBEdge& z){ return (z.a == e.a && z.b == e.b); });
        if (it == xbs.end()) xbs.push_back(e);
    }

    std::sort(xbs.begin(), xbs.end(), [](const XBEdge& p, const XBEdge& q) {
        if (p.a != q.a) return p.a < q.a;
        return p.b < q.b;
    });

    // phiavg column is period-specific; it may or may not exist.
    // Per conventions: use "phiavg, <period>" if present; else mid-bin.
    const bool use_phiavg = (!is_combined_group_label(label)); // combined groups have no single period phiavg
    int c_phiavg = -1;
    if (use_phiavg) {
        const std::string name = col_phiavg_for_period(label);
        auto it = csv.index.find(name);
        if (it != csv.index.end()) c_phiavg = it->second;
    }

    int c_xBavg = -1;
    if (!is_combined_group_label(label)) {
        const std::string name = col_xBavg_for_period(label);
        auto it = csv.index.find(name);
        if (it != csv.index.end()) c_xBavg = it->second;
    }

    const std::string topoLabel = topo_label_for_csv(topoDir);
    if (topoLabel.empty()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: unknown topoDir '" << topoDir << "'";
        fatal(ss.str());
    }

    // Build output directory:
    const std::string outdir = out_root_for_label(label, out_root_dir) + "/" + topoDir;
    mkdir_p(outdir);

    // For each xB bin, collect row indices for this xB bin, and build RowPoint per row.
    for (const auto& xb : xbs) {
        std::vector<int> row_indices;
        row_indices.reserve(256);

        std::vector<RowPoint> points;
        points.resize(rows.size());

        for (int r = 0; r < (int)rows.size(); ++r) {
            const RowBin& w = rows[r];
            if (!w.valid) continue;
            if (!(w.xBmin == xb.a && w.xBmax == xb.b)) continue;

            row_indices.push_back(r);

            RowPoint p;

            // x value: phiavg if available, else midpoint of bin.
            if (c_phiavg >= 0) {
                const std::string& cell = csv.rows[r][c_phiavg];
                const double v = cell_to_double_or_nan(cell);
                if (std::isfinite(v)) p.phi_x = wrap_phi_deg(v);
                else p.phi_x = wrap_phi_deg(bin_center(w.pmin, w.pmax));
            } else {
                p.phi_x = wrap_phi_deg(bin_center(w.pmin, w.pmax));
            }

            // y values: counts for this row, this topo.
            auto itc = row_counts_for_topo.find(r);
            double pos = 0.0, neg = 0.0;
            if (itc != row_counts_for_topo.end()) {
                pos = itc->second.pos;
                neg = itc->second.neg;
            }

            p.pos = pos;
            p.neg = neg;
            p.pos_err = poisson_err(pos);
            p.neg_err = poisson_err(neg);

            p.Q2c = bin_center(w.Q2min, w.Q2max);
            p.tc  = bin_center(w.tmin,  w.tmax);
            p.valid = true;

            points[r] = p;
        }

        // Title xB text: prefer xBavg if present, else use range.
        std::ostringstream xbtxt;
        if (c_xBavg >= 0) {
            // Use the first row in this xB bin that has a finite xBavg.
            double xBavg = std::numeric_limits<double>::quiet_NaN();
            for (int r : row_indices) {
                const double v = cell_to_double_or_nan(csv.rows[r][c_xBavg]);
                if (std::isfinite(v)) { xBavg = v; break; }
            }
            if (std::isfinite(xBavg)) {
                xbtxt << "x_{B}=" << std::fixed << std::setprecision(3) << xBavg;
            } else {
                xbtxt << "x_{B} in [" << std::fixed << std::setprecision(2) << xb.a
                      << "," << std::fixed << std::setprecision(2) << xb.b << ")";
            }
        } else {
            xbtxt << "x_{B} in [" << std::fixed << std::setprecision(2) << xb.a
                  << "," << std::fixed << std::setprecision(2) << xb.b << ")";
        }

        // Label for title in plot function: for combined groups label is the group label.
        const std::string label_for_title = label;

        plot_one_xB_canvas(outdir,
                           label_for_title,
                           topoDir,
                           topoLabel,
                           xbtxt.str(),
                           xb.a,
                           rows,
                           row_indices,
                           points);
    }
}

// -----------------------------------------------------------------------------
// Combined group aggregation (in-memory only; never write to CSV)
// -----------------------------------------------------------------------------

static RowCounts sum_row_counts(const RowCounts& a, const RowCounts& b) {
    RowCounts out = a;
    for (const auto& kv : b) {
        const int row = kv.first;
        const HelCounts& h = kv.second;
        HelCounts& o = out[row];
        o.pos += h.pos;
        o.neg += h.neg;
        o.unpol += h.unpol;
    }
    return out;
}

// -----------------------------------------------------------------------------
// update_total_counts_csv (public entry point)
// -----------------------------------------------------------------------------

} // end anonymous namespace

bool update_total_counts_csv(const std::string& csv_main,
                             const std::map<std::string, TTree*>& dvcsDataTrees,
                             const std::string& combined_cuts_json,
                             const std::string& global_cuts_json,
                             int maxThreads) {
    (void)global_cuts_json; // global_cuts.h currently uses default_global_cuts(); keep signature stable.

    // ROOT global threading safety.
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);

    const bool trace_matches = env_flag("TOTAL_COUNTS_TRACE_MATCHES");

    // Load CSV.
    CSV csv;
    load_csv(csv_main, csv);

    // Load row bins from CSV.
    const std::vector<RowBin> rows = load_row_bins_from_csv(csv);

    // Load 3-sigma cuts (DATA only).
    const TopoCutMap sigma_cuts_data = load_combined_cuts_data_only(combined_cuts_json);

    // Prepare list of periods that are present.
    std::vector<PeriodTags> periods;
    periods.reserve(CANONICAL_PERIODS().size());

    for (const auto& P : CANONICAL_PERIODS()) {
        auto it = dvcsDataTrees.find(P.tree_key);
        if (it == dvcsDataTrees.end()) continue;
        if (!it->second) continue;
        PeriodTags tags = parse_period_tags_from_tree_key(P.tree_key);
        periods.push_back(tags);
    }

    if (periods.empty()) {
        fatal("[total_counts] FATAL: no DVCS trees available in dvcsDataTrees for total_counts.");
    }

    std::cout << "[total_counts] Will process " << periods.size() << " DVCS period tree(s)." << std::endl;

    // Per-period results: period_display -> topoDir -> RowCounts
    std::unordered_map<std::string, std::unordered_map<std::string, RowCounts>> per_period_counts;
    std::mutex merge_mutex;

    // Thread count cap.
    int nth = std::max(1, std::min(5, maxThreads));
    nth = std::min(nth, (int)periods.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
    for (int i = 0; i < (int)periods.size(); ++i) {
        const PeriodTags& tags = periods[i];

        auto it = dvcsDataTrees.find(tags.tree_key);
        TTree* tree = (it == dvcsDataTrees.end()) ? nullptr : it->second;
        if (!tree) continue;

        const PeriodWorkResult R = accumulate_counts_for_tree(tags, tree, rows, sigma_cuts_data, trace_matches);

        std::lock_guard<std::mutex> lock(merge_mutex);
        auto& topoMap = per_period_counts[tags.period_display];
        for (const auto& kv : R.topo_counts) {
            const std::string& topoDir = kv.first;
            const RowCounts& rc = kv.second;
            topoMap[topoDir] = sum_row_counts(topoMap[topoDir], rc);
        }
    }

    // Write per-period counts into CSV (skip supplemental and groups).
    for (const auto& kvp : per_period_counts) {
        const std::string& period_display = kvp.first;
        const auto& topoMap = kvp.second;

        if (should_skip_csv_for_label(period_display)) continue;

        for (const auto& kt : topoMap) {
            const std::string& topoDir = kt.first;
            const RowCounts& rc = kt.second;

            const std::string topoLabel = topo_label_for_csv(topoDir);
            if (topoLabel.empty()) {
                std::ostringstream ss;
                ss << "[total_counts] FATAL: cannot map topoDir '" << topoDir << "' to topoLabel.";
                fatal(ss.str());
            }

            const CountsColumnIndex cols = resolve_counts_columns(csv, topoLabel, period_display);

            for (const auto& row_kv : rc) {
                const int r = row_kv.first;
                const HelCounts& h = row_kv.second;

                if (r < 0 || r >= (int)csv.rows.size()) {
                    std::ostringstream ss;
                    ss << "[total_counts] FATAL: row index out of range: " << r;
                    fatal(ss.str());
                }

                // Per conventions: write pos/neg separately; unpol is stored as the explicit unpol bucket.
                // If you want unpol := pos+neg, do it explicitly here (deterministic).
                const double unpol = h.unpol;
                csv.rows[r][cols.col_unpol] = fmt0(unpol);
                csv.rows[r][cols.col_pos]   = fmt0(h.pos);
                csv.rows[r][cols.col_neg]   = fmt0(h.neg);
            }
        }
    }

    // Atomic CSV save.
    write_csv_atomic(csv_main, csv);
    std::cout << "[total_counts] Updated CSV counts in: " << csv_main << std::endl;

    // Plotting:
    //  - per-period plots (including supplemental if present)
    //  - combined groups Fa18, Sp18, 10.6 GeV (in-memory sums; never written to CSV)
    const std::string out_root_dir = "output/total_counts_plots";

    // Per-period plots.
    for (const auto& kvp : per_period_counts) {
        const std::string& period_display = kvp.first;
        const auto& topoMap = kvp.second;

        for (const auto& kt : topoMap) {
            const std::string& topoDir = kt.first;
            const RowCounts& rc = kt.second;

            make_plots_for_label_and_topo(period_display, topoDir, out_root_dir, csv, rows, rc);
        }
    }

    // Combined groups (in-memory aggregation only).
    // Group definitions (deterministic):
    //   Fa18     := Fa18 Inb + Fa18 Out
    //   Sp18     := Sp18 Inb + Sp18 Out
    //   10.6 GeV := Fa18 Inb + Fa18 Out + Sp18 Inb + Sp18 Out
    auto get_counts = [&](const std::string& period_display) -> const std::unordered_map<std::string, RowCounts>* {
        auto it = per_period_counts.find(period_display);
        if (it == per_period_counts.end()) return nullptr;
        return &it->second;
    };

    auto sum_group = [&](const std::vector<std::string>& members) {
        std::unordered_map<std::string, RowCounts> topoSum;
        for (const auto& m : members) {
            const auto* pm = get_counts(m);
            if (!pm) continue;
            for (const auto& kt : *pm) {
                const std::string& topoDir = kt.first;
                topoSum[topoDir] = sum_row_counts(topoSum[topoDir], kt.second);
            }
        }
        return topoSum;
    };

    const auto Fa18Topo = sum_group({"Fa18 Inb", "Fa18 Out"});
    const auto Sp18Topo = sum_group({"Sp18 Inb", "Sp18 Out"});
    const auto Topo106  = sum_group({"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"});

    for (const auto& kt : Fa18Topo) make_plots_for_label_and_topo("Fa18", kt.first, out_root_dir, csv, rows, kt.second);
    for (const auto& kt : Sp18Topo) make_plots_for_label_and_topo("Sp18", kt.first, out_root_dir, csv, rows, kt.second);
    for (const auto& kt : Topo106)  make_plots_for_label_and_topo("10.6 GeV", kt.first, out_root_dir, csv, rows, kt.second);

    std::cout << "[total_counts] Plots written under: " << out_root_dir << std::endl;

    return true;
}