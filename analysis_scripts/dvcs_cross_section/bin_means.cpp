#include "bin_means.h"

// Optional period helper if you have it; not required here.
// #include "periods.h"

// Optional 3-sigma exclusivity if you have it; we guard with __has_include.
#if __has_include("exclusivity_cuts.h")
  #include "exclusivity_cuts.h"
  #define HAVE_EXCL_3SIGMA 1
#else
  #define HAVE_EXCL_3SIGMA 0
#endif

#include <TTree.h>
#include <TMath.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// ------------------------------
// Minimal CSV utilities
// ------------------------------
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
            for (char ch : s) { if (ch == '"') oss << "\"\""; else oss << ch; }
            oss << '"';
        } else {
            oss << s;
        }
        if (i + 1 < fields.size()) oss << ',';
    }
    return oss.str();
}

static std::unordered_map<std::string, int>
build_header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> m;
    for (int i = 0; i < (int)header.size(); ++i) m[header[i]] = i;
    return m;
}

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    return s;
}

static bool as_double(const std::string& s, double& out) {
    if (s.empty()) return false;
    char* endp = nullptr;
    out = std::strtod(s.c_str(), &endp);
    return (endp != s.c_str());
}

static int as_int_def(const std::string& s, int def) {
    if (s.empty()) return def;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return def;
    return (int)v;
}

// ------------------------------
// Period labeling helpers
// ------------------------------
static std::string csv_label_for_period_key(std::string key) {
    std::string k = to_lower(key);
    if (k.find("fa18") != std::string::npos && k.find("inb") != std::string::npos) return "Fa18 Inb";
    if (k.find("fa18") != std::string::npos && k.find("out") != std::string::npos) return "Fa18 Out";
    if (k.find("sp18") != std::string::npos && k.find("inb") != std::string::npos) return "Sp18 Inb";
    if (k.find("sp18") != std::string::npos && k.find("out") != std::string::npos) return "Sp18 Out";
    if (k.find("sp19") != std::string::npos && k.find("inb") != std::string::npos) return "Sp19 Inb";
    // Fallback: return original key if we do not recognize it.
    return key;
}

static bool is_10p6_label(const std::string& csv_label) {
    return (csv_label == "Fa18 Inb" || csv_label == "Fa18 Out" ||
            csv_label == "Sp18 Inb" || csv_label == "Sp18 Out");
}

// ------------------------------
// Topology helpers
// ------------------------------
enum class Topology {
    FD_FD,
    CD_FD,
    CD_FT
};

static const std::vector<Topology>& all_topos() {
    static const std::vector<Topology> T = {Topology::FD_FD, Topology::CD_FD, Topology::CD_FT};
    return T;
}

static std::string topo_to_csv_title(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "(FD, FD)";
        case Topology::CD_FD: return "(CD, FD)";
        case Topology::CD_FT: return "(CD, FT)";
    }
    return "(FD, FD)";
}

static bool event_matches_topology(int detector1, int detector2, Topology t) {
    switch (t) {
        case Topology::FD_FD: return (detector1 == 1 && detector2 == 1);
        case Topology::CD_FD: return (detector1 == 2 && detector2 == 1);
        case Topology::CD_FT: return (detector1 == 2 && detector2 == 0);
    }
    return false;
}

// Determine per-row allowed topologies by inspecting raw-yield columns.
// If exactly one topology has any non-empty "unpol" yield across any of the 5 periods,
// we lock to that topology. Otherwise we accept all.
static std::vector<Topology> detect_row_topos(const std::vector<std::string>& row,
                                              const std::vector<std::string>& header) {
    auto H = build_header_index(header);
    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    std::set<Topology> seen;

    for (Topology t : all_topos()) {
        std::string topo_name = topo_to_csv_title(t);
        for (const auto& per : periods) {
            std::ostringstream col;
            col << "raw yield, ep->epg, " << topo_name << ", exp, " << per << ", unpol";
            auto it = H.find(col.str());
            if (it == H.end()) continue;
            int idx = it->second;
            if (idx < 0 || idx >= (int)row.size()) continue;
            const std::string& cell = row[idx];
            if (!cell.empty()) {
                // Non-empty is enough; optionally check numeric > 0.0
                double v = 0.0;
                if (as_double(cell, v)) {
                    if (v > 0.0) seen.insert(t);
                } else {
                    // Any text means "present"; accept it.
                    seen.insert(t);
                }
            }
        }
    }

    if (seen.size() == 1) return { *seen.begin() };
    // Ambiguous or none: accept all
    return all_topos();
}

// ------------------------------
// Branch binder
// ------------------------------
struct BranchBinder {
    // topology
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    // global+exclusivity cuts
    double t1 = 0.0;             bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open = false;
    double pTmiss = 0.0;         bool has_pTmiss = false;

    // binning vars
    double x = 0.0;              bool has_x = false;
    double Q2 = 0.0;             bool has_Q2 = false;
    double phi2 = 0.0;           bool has_phi2 = false;
    double Delta_phi = 0.0;      bool has_Delta = false;

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
        bindD("open_angle_ep2", &open_angle_ep2, has_open);
        bindD("pTmiss", &pTmiss, has_pTmiss);

        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);
        bindD("phi2", &phi2, has_phi2);
        bindD("Delta_phi", &Delta_phi, has_Delta);
    }

    double phi_rad() const {
        if (has_phi2) return phi2;
        if (has_Delta) return Delta_phi;
        return std::numeric_limits<double>::quiet_NaN();
    }

    bool ready_for_cuts() const {
        return has_detector1 && has_detector2 && has_t1 && has_open && has_pTmiss;
    }
    bool ready_for_bins() const {
        return has_x && has_Q2 && (has_phi2 || has_Delta);
    }
};

// ------------------------------
// Cuts
// ------------------------------
static bool apply_global_cuts(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

static bool apply_excl_3sigma(const std::string& period_key,
                              Topology topo,
                              const BranchBinder& b) {
#if HAVE_EXCL_3SIGMA
    // If you have a real implementation, wire it here.
    // Example signature:
    // return exclusivity_cuts::passes3sigma(period_key, topo_to_csv_title(topo),
    //                                       b.t1, b.open_angle_ep2, b.pTmiss /*, ...*/);
    return exclusivity_cuts::passes3sigma(period_key, topo_to_csv_title(topo),
                                          b.t1, b.open_angle_ep2, b.pTmiss);
#else
    (void)period_key; (void)topo; (void)b;
    // Fallback to "true": only global cuts are applied.
    return true;
#endif
}

// ------------------------------
// Bin structures built from CSV
// ------------------------------
struct RowDef {
    int row_index = -1;          // CSV row index
    double xBmin = 0.0, xBmax = 0.0;
    double Q2min = 0.0, Q2max = 0.0;
    double tmin  = 0.0, tmax  = 0.0;     // absolute t
    double phimin = 0.0, phimax = 0.0;   // in CSV units (deg or rad)
    std::vector<Topology> allowed_topos; // detected for this row
    bool valid = true;
};

struct CellKey {
    int ix = -1, iQ = -1, it = -1;
    bool operator<(const CellKey& o) const {
        if (ix != o.ix) return ix < o.ix;
        if (iQ != o.iQ) return iQ < o.iQ;
        return it < o.it;
    }
};

static bool in_range_inclusive_low_exclusive_high(double v, double a, double b) {
    return (v >= a && v < b);
}

static bool phi_in_interval(double phi, double a, double b, double period) {
    // Handles wrap-around intervals; phi, a, b in same units in [0,period)
    if (a <= b) {
        return (phi >= a && phi < b);
    } else {
        // wrap across boundary
        return (phi >= a || phi < b);
    }
}

// ------------------------------
// Accumulators
// ------------------------------
struct Accum {
    double sum_xB = 0.0;
    double sum_Q2 = 0.0;
    double sum_t  = 0.0;
    // Circular mean for phi:
    double sum_c = 0.0; // sum cos(phi)
    double sum_s = 0.0; // sum sin(phi)
    long long n = 0;
    void add(double xB, double Q2, double tabs, double phi, bool phi_is_deg) {
        sum_xB += xB;
        sum_Q2 += Q2;
        sum_t  += tabs;
        double ang = phi_is_deg ? (phi * M_PI / 180.0) : phi;
        sum_c += std::cos(ang);
        sum_s += std::sin(ang);
        ++n;
    }
    void mean(double& m_xB, double& m_Q2, double& m_t, double& m_phi, bool phi_out_deg) const {
        if (n <= 0) { m_xB = m_Q2 = m_t = m_phi = std::numeric_limits<double>::quiet_NaN(); return; }
        m_xB = sum_xB / (double)n;
        m_Q2 = sum_Q2 / (double)n;
        m_t  = sum_t  / (double)n;
        double ang = std::atan2(sum_s, sum_c); // [-pi,pi]
        if (ang < 0.0) ang += 2.0 * M_PI;
        m_phi = phi_out_deg ? (ang * 180.0 / M_PI) : ang;
    }
};

struct PeriodResults {
    // One accumulator per CSV row.
    std::vector<Accum> per_row;
};

// ------------------------------
// Worker: process one period
// ------------------------------
static PeriodResults process_period(const std::string& period_key,
                                    TTree* tree,
                                    const std::vector<RowDef>& rows,
                                    const std::vector<double>& x_edges_min,
                                    const std::vector<double>& x_edges_max,
                                    const std::vector<double>& Q_edges_min,
                                    const std::vector<double>& Q_edges_max,
                                    const std::vector<double>& t_edges_min,
                                    const std::vector<double>& t_edges_max,
                                    const std::map<CellKey, std::vector<int>>& cell_to_rows,
                                    bool csv_phi_is_deg) {
    PeriodResults R;
    R.per_row.assign(rows.size(), Accum{});

    if (!tree) return R;

    BranchBinder b; b.bind(tree);
    if (!b.ready_for_cuts() || !b.ready_for_bins()) {
        std::cerr << "[bin_means] FATAL: Tree for '" << period_key
                  << "' lacks required branches.\n";
        return R;
    }

    const Long64_t nent = tree->GetEntries();
    const double phi_period = csv_phi_is_deg ? 360.0 : (2.0 * M_PI);

    auto find_index = [](double v, const std::vector<double>& mins,
                         const std::vector<double>& maxs) -> int {
        for (int i = 0; i < (int)mins.size(); ++i) {
            if (v >= mins[i] && v < maxs[i]) return i;
        }
        return -1;
    };

    for (Long64_t i = 0; i < nent; ++i) {
        tree->GetEntry(i);

        // quick global cuts
        if (!apply_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

        const double xB = b.x;
        const double Q2 = b.Q2;
        const double tabs = std::fabs(b.t1);
        double phi_rad = b.phi_rad();
        if (!std::isfinite(xB) || !std::isfinite(Q2) ||
            !std::isfinite(tabs) || !std::isfinite(phi_rad)) continue;

        // Convert to CSV phi units and wrap
        double phi_csv = csv_phi_is_deg ? (phi_rad * 180.0 / M_PI) : phi_rad;
        // wrap into [0,period)
        phi_csv = std::fmod(phi_csv, phi_period);
        if (phi_csv < 0.0) phi_csv += phi_period;

        // locate (ix,iQ,it)
        int ix = find_index(xB, x_edges_min, x_edges_max);
        if (ix < 0) continue;
        int iQ = find_index(Q2, Q_edges_min, Q_edges_max);
        if (iQ < 0) continue;
        int it = find_index(tabs, t_edges_min, t_edges_max);
        if (it < 0) continue;

        CellKey ck{ix, iQ, it};
        auto it_rows = cell_to_rows.find(ck);
        if (it_rows == cell_to_rows.end()) continue;

        // Try each candidate phi row in this cell
        for (int row_idx : it_rows->second) {
            const RowDef& rd = rows[row_idx];

            // Topology filter
            bool topo_ok = false;
            for (Topology t : rd.allowed_topos) {
                if (event_matches_topology(b.detector1, b.detector2, t)) {
                    // 3-sigma exclusivity for this topology (if available)
                    if (apply_excl_3sigma(period_key, t, b)) {
                        topo_ok = true;
                        break;
                    }
                }
            }
            if (!topo_ok) continue;

            // Phi interval test (with wrap support)
            if (!phi_in_interval(phi_csv, rd.phimin, rd.phimax, phi_period)) continue;

            // Accumulate
            R.per_row[row_idx].add(xB, Q2, tabs, phi_csv, csv_phi_is_deg);
            // Note: do not break; if multiple overlapping phi rows existed (should not),
            // we would add more than once. Binning is assumed disjoint.
            break;
        }
    }

    return R;
}

// ------------------------------
// Build bin grid and mapping from CSV
// ------------------------------
struct GridBuild {
    std::vector<RowDef> rows;
    std::vector<double> x_edges_min, x_edges_max;
    std::vector<double> Q_edges_min, Q_edges_max;
    std::vector<double> t_edges_min, t_edges_max;
    std::map<CellKey, std::vector<int>> cell_to_rows;
    bool csv_phi_is_deg = true;
};

static bool parse_double_or_die(const std::string& name,
                                const std::vector<std::string>& row,
                                const std::unordered_map<std::string,int>& H,
                                double& out) {
    auto it = H.find(name);
    if (it == H.end()) return false;
    if (!as_double(row[it->second], out)) return false;
    return true;
}

static GridBuild build_grid_from_csv(const std::vector<std::vector<std::string>>& table,
                                     const std::vector<std::string>& header) {
    GridBuild G;

    auto H = build_header_index(header);

    // Detect phi unit: if any phimax > ~7.0, assume degrees
    bool any_phimax_gt7 = false;
    for (size_t r = 1; r < table.size(); ++r) {
        const auto& row = table[r];
        auto it = H.find("phimax");
        if (it == H.end() || it->second >= (int)row.size()) continue;
        double v = 0.0;
        if (as_double(row[it->second], v) && v > 7.0) { any_phimax_gt7 = true; break; }
    }
    G.csv_phi_is_deg = any_phimax_gt7;

    // Collect unique edges in xB, Q2, t_abs
    std::set<std::pair<double,double>> x_ranges, Q_ranges, t_ranges;

    for (size_t r = 1; r < table.size(); ++r) {
        const auto& row = table[r];
        int valid = 0;
        auto it_valid = H.find("valid bin");
        if (it_valid != H.end() && it_valid->second < (int)row.size()) {
            valid = as_int_def(row[it_valid->second], 0);
        }
        if (valid != 1) continue;

        RowDef rd;
        rd.row_index = (int)r;

        if (!parse_double_or_die("xBmin", row, H, rd.xBmin)) continue;
        if (!parse_double_or_die("xBmax", row, H, rd.xBmax)) continue;
        if (!parse_double_or_die("Q2min", row, H, rd.Q2min)) continue;
        if (!parse_double_or_die("Q2max", row, H, rd.Q2max)) continue;
        if (!parse_double_or_die("t_abs_min", row, H, rd.tmin)) continue;
        if (!parse_double_or_die("t_abs_max", row, H, rd.tmax)) continue;
        if (!parse_double_or_die("phimin", row, H, rd.phimin)) continue;
        if (!parse_double_or_die("phimax", row, H, rd.phimax)) continue;

        // Normalize phi into [0,period)
        const double phi_period = G.csv_phi_is_deg ? 360.0 : (2.0 * M_PI);
        auto wrap = [&](double a)->double {
            double w = std::fmod(a, phi_period);
            if (w < 0.0) w += phi_period;
            return w;
        };
        rd.phimin = wrap(rd.phimin);
        rd.phimax = wrap(rd.phimax);

        // Detect allowed topologies from row contents
        rd.allowed_topos = detect_row_topos(row, header);
        G.rows.push_back(rd);

        x_ranges.emplace(rd.xBmin, rd.xBmax);
        Q_ranges.emplace(rd.Q2min, rd.Q2max);
        t_ranges.emplace(rd.tmin,  rd.tmax);
    }

    // Expand edges into vectors (preserving CSV order as seen)
    for (const auto& p : x_ranges) { G.x_edges_min.push_back(p.first);  G.x_edges_max.push_back(p.second); }
    for (const auto& p : Q_ranges) { G.Q_edges_min.push_back(p.first);  G.Q_edges_max.push_back(p.second); }
    for (const auto& p : t_ranges) { G.t_edges_min.push_back(p.first);  G.t_edges_max.push_back(p.second); }

    // Build mapping from (ix,iQ,it) to all phi rows
    auto find_index = [](double a, double b, const std::vector<double>& mins,
                         const std::vector<double>& maxs) -> int {
        for (int i = 0; i < (int)mins.size(); ++i) {
            if (a == mins[i] && b == maxs[i]) return i;
        }
        return -1;
    };

    for (const RowDef& rd : G.rows) {
        int ix = find_index(rd.xBmin, rd.xBmax, G.x_edges_min, G.x_edges_max);
        int iQ = find_index(rd.Q2min, rd.Q2max, G.Q_edges_min, G.Q_edges_max);
        int it = find_index(rd.tmin,  rd.tmax,  G.t_edges_min, G.t_edges_max);
        if (ix < 0 || iQ < 0 || it < 0) continue;
        G.cell_to_rows[CellKey{ix,iQ,it}].push_back(rd.row_index);
    }

    return G;
}

// ------------------------------
// Column index helpers for writing results
// ------------------------------
struct AvgColIdx {
    // For each variable, eight columns in order:
    //   Fa18 Inb, Fa18 Out, Sp19 Inb, Sp18 Inb, Sp18 Out, Fa18, Sp18, 10.6 GeV
    int xBavg[8]  = {-1,-1,-1,-1,-1,-1,-1,-1};
    int Q2avg[8]  = {-1,-1,-1,-1,-1,-1,-1,-1};
    int tabsavg[8]= {-1,-1,-1,-1,-1,-1,-1,-1};
    int phiavg[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
};

static int find_header_index(const std::vector<std::string>& header,
                             const std::string& name) {
    for (int i = 0; i < (int)header.size(); ++i) {
        if (header[i] == name) return i;
    }
    return -1;
}

static bool locate_avg_columns(const std::vector<std::string>& header, AvgColIdx& idx) {
    static const char* groups[8] = {
        "Fa18 Inb","Fa18 Out","Sp19 Inb","Sp18 Inb","Sp18 Out","Fa18","Sp18","10.6 GeV"
    };
    bool ok = true;
    for (int g = 0; g < 8; ++g) {
        {
            std::ostringstream n; n << "xBavg, " << groups[g];
            idx.xBavg[g] = find_header_index(header, n.str());
            ok = ok && (idx.xBavg[g] >= 0);
        }
        {
            std::ostringstream n; n << "Q2avg, " << groups[g];
            idx.Q2avg[g] = find_header_index(header, n.str());
            ok = ok && (idx.Q2avg[g] >= 0);
        }
        {
            std::ostringstream n; n << "t_abs_avg, " << groups[g];
            idx.tabsavg[g] = find_header_index(header, n.str());
            ok = ok && (idx.tabsavg[g] >= 0);
        }
        {
            std::ostringstream n; n << "phiavg, " << groups[g];
            idx.phiavg[g] = find_header_index(header, n.str());
            ok = ok && (idx.phiavg[g] >= 0);
        }
    }
    return ok;
}

// ------------------------------
// Main driver
// ------------------------------
bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers) {
    // 1) Load CSV into memory
    std::ifstream fin(csv_path);
    if (!fin.is_open()) {
        std::cerr << "[bin_means] ERROR: Could not open CSV: " << csv_path << "\n";
        return false;
    }
    std::vector<std::vector<std::string>> table;
    std::string line;
    while (std::getline(fin, line)) table.push_back(split_csv_line(line));
    fin.close();
    if (table.empty()) {
        std::cerr << "[bin_means] ERROR: Empty CSV: " << csv_path << "\n";
        return false;
    }
    const std::vector<std::string> header = table[0];

    // 2) Build grid from CSV rows
    GridBuild G = build_grid_from_csv(table, header);
    if (G.rows.empty()) {
        std::cerr << "[bin_means] ERROR: No valid rows found in CSV.\n";
        return false;
    }

    // 3) Prepare per-period workers
    //    We expect keys like DVCS_Fa18_inb, DVCS_Fa18_out, DVCS_Sp19_inb, DVCS_Sp18_inb, DVCS_Sp18_out.
    struct WorkItem { std::string key; TTree* tree; };
    std::vector<WorkItem> work;
    work.reserve(dataTrees.size());
    for (const auto& kv : dataTrees) {
        if (kv.second) work.push_back(WorkItem{kv.first, kv.second});
    }

    // Limit workers to 5
    if (max_workers <= 0) max_workers = 1;
    if (max_workers > 5)  max_workers = 5;

    // Launch one task per period; if more than max_workers, run some in-band after futures complete.
    std::vector<std::future<PeriodResults>> futures;
    futures.reserve(work.size());

    auto submit = [&](const WorkItem& w) {
        return std::async(std::launch::async,
            [&G, &w]() -> PeriodResults {
                return process_period(w.key, w.tree,
                                      G.rows,
                                      G.x_edges_min, G.x_edges_max,
                                      G.Q_edges_min, G.Q_edges_max,
                                      G.t_edges_min, G.t_edges_max,
                                      G.cell_to_rows,
                                      G.csv_phi_is_deg);
            });
    };

    size_t next = 0;
    while (next < work.size()) {
        // Fill up to max_workers outstanding
        while (next < work.size() && futures.size() < (size_t)max_workers) {
            futures.push_back(submit(work[next]));
            ++next;
        }
        // Wait for the first future to complete to free a slot
        if (!futures.empty()) {
            futures.front().wait();
            // Compact: move completed future to a results array
            // but we can just rotate vector
            std::rotate(futures.begin(), futures.begin() + 1, futures.end());
        }
    }

    // At this point, all tasks have been submitted; wait for all to finish and collect.
    std::vector<PeriodResults> per_period_results;
    per_period_results.reserve(work.size());
    for (auto& f : futures) per_period_results.push_back(f.get());

    // 4) Build maps from period label to results, keeping order:
    //    idx order: 0=Fa18 Inb, 1=Fa18 Out, 2=Sp19 Inb, 3=Sp18 Inb, 4=Sp18 Out
    struct LabeledResults { std::string label; PeriodResults R; };
    std::vector<LabeledResults> ordered(5); // we will fill by label
    auto place = [&](const std::string& label, PeriodResults&& R) {
        if (label == "Fa18 Inb") ordered[0] = {label, std::move(R)};
        else if (label == "Fa18 Out") ordered[1] = {label, std::move(R)};
        else if (label == "Sp19 Inb") ordered[2] = {label, std::move(R)};
        else if (label == "Sp18 Inb") ordered[3] = {label, std::move(R)};
        else if (label == "Sp18 Out") ordered[4] = {label, std::move(R)};
        else {
            // Unknown label; ignore silently
        }
    };

    for (size_t i = 0; i < work.size(); ++i) {
        std::string lbl = csv_label_for_period_key(work[i].key);
        place(lbl, std::move(per_period_results[i]));
    }

    // 5) Locate output columns for averages
    AvgColIdx Aidx;
    if (!locate_avg_columns(header, Aidx)) {
        std::cerr << "[bin_means] ERROR: Could not find grouped-average columns in CSV header.\n";
        return false;
    }

    auto write_cell = [&](std::vector<std::string>& row, int col, double v) {
        if (col < 0 || col >= (int)row.size()) return;
        if (!std::isfinite(v)) { row[col].clear(); return; }
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(8) << v;
        row[col] = ss.str();
    };

    // 6) For each row: compute means for 5 periods and 3 combined groups
    //    groups order for columns:
    //    0 Fa18 Inb, 1 Fa18 Out, 2 Sp19 Inb, 3 Sp18 Inb, 4 Sp18 Out, 5 Fa18, 6 Sp18, 7 10.6 GeV
    for (const RowDef& rd : G.rows) {
        // Gather per-period accumulators for this row
        Accum per[5];
        bool have[5] = {false,false,false,false,false};

        auto pull = [&](int slot, const std::string& label) {
            if (ordered[slot].label != label) return;
            if ((int)ordered[slot].R.per_row.size() <= rd.row_index) return;
            per[slot] = ordered[slot].R.per_row[rd.row_index];
            have[slot] = true;
        };

        pull(0, "Fa18 Inb");
        pull(1, "Fa18 Out");
        pull(2, "Sp19 Inb");
        pull(3, "Sp18 Inb");
        pull(4, "Sp18 Out");

        auto combine = [&](const std::vector<int>& idxs) -> Accum {
            Accum out{};
            for (int j : idxs) {
                if (!have[j]) continue;
                out.sum_xB += per[j].sum_xB;
                out.sum_Q2 += per[j].sum_Q2;
                out.sum_t  += per[j].sum_t;
                out.sum_c  += per[j].sum_c;
                out.sum_s  += per[j].sum_s;
                out.n      += per[j].n;
            }
            return out;
        };

        Accum fa18 = combine({0,1});
        Accum sp18 = combine({3,4});
        Accum g106 = combine({0,1,3,4}); // only 10.6 GeV periods

        // Compute means and write into CSV cells
        auto write_group = [&](int gslot, const Accum& acc) {
            double mx = 0, mQ = 0, mt = 0, mphi = 0;
            acc.mean(mx, mQ, mt, mphi, G.csv_phi_is_deg);
            write_cell(table[rd.row_index], Aidx.xBavg[gslot], mx);
            write_cell(table[rd.row_index], Aidx.Q2avg[gslot], mQ);
            write_cell(table[rd.row_index], Aidx.tabsavg[gslot], mt);
            write_cell(table[rd.row_index], Aidx.phiavg[gslot], mphi);
        };

        // Period groups (0..4)
        for (int g = 0; g < 5; ++g) write_group(g, per[g]);
        // Combined groups
        write_group(5, fa18);
        write_group(6, sp18);
        write_group(7, g106);
    }

    // 7) Write CSV back to disk (in place)
    std::ofstream fout(csv_path);
    if (!fout.is_open()) {
        std::cerr << "[bin_means] ERROR: Could not open CSV for writing: " << csv_path << "\n";
        return false;
    }
    fout << join_csv_row(header) << "\n";
    for (size_t r = 1; r < table.size(); ++r) {
        fout << join_csv_row(table[r]) << "\n";
    }
    fout.close();

    std::cout << "[bin_means] Updated grouped bin-averaged kinematics in " << csv_path << "\n";
    return true;
}