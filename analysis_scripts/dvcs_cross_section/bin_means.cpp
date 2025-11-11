#include "bin_means.h"
#include "exclusivity_cuts.h"   // passes3sigma(...) + Topology enum + topo_to_csv_title(...)
#include "global_cuts.h"        // global_cuts::passes(...)
#include "periods.h"            // CANONICAL_PERIODS() if needed

#include <TTree.h>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// ------------------------------
// Helpers: CSV minimal parsing
// ------------------------------
static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur; cur.reserve(line.size());
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
        bool need_quotes = s.find(',') != std::string::npos || s.find('"') != std::string::npos;
        if (need_quotes) {
            oss << '"';
            for (char ch : s) oss << (ch == '"' ? "\"\"" : std::string(1, ch));
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

static bool get_double(const std::vector<std::string>& row,
                       const std::unordered_map<std::string,int>& H,
                       const std::string& name, double& out) {
    auto it = H.find(name);
    if (it == H.end()) return false;
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return false;
    if (row[j].empty()) return false;
    char* endp = nullptr;
    out = std::strtod(row[j].c_str(), &endp);
    return endp != row[j].c_str();
}

static int get_int(const std::vector<std::string>& row,
                   const std::unordered_map<std::string,int>& H,
                   const std::string& name, int& out) {
    auto it = H.find(name);
    if (it == H.end()) return 0;
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return 0;
    if (row[j].empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(row[j].c_str(), &endp, 10);
    if (endp == row[j].c_str()) return 0;
    out = (int)v;
    return 1;
}

static inline double wrap_deg_360(double phi_deg) {
    if (!std::isfinite(phi_deg)) return std::numeric_limits<double>::quiet_NaN();
    double w = std::fmod(phi_deg, 360.0);
    if (w < 0.0) w += 360.0;
    if (w >= 360.0) w = std::nextafter(360.0, 0.0);
    return w;
}

// ------------------------------
// Branch binder (minimal)
// ------------------------------
struct BranchBinder {
    int    detector1 = 0; bool has_detector1 = false;
    int    detector2 = 0; bool has_detector2 = false;

    double t1 = 0.0;             bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open_angle_ep2 = false;
    double pTmiss = 0.0;         bool has_pTmiss = false;

    double x = 0.0;              bool has_x = false;
    double Q2 = 0.0;             bool has_Q2 = false;
    double phi2 = 0.0;           bool has_phi2 = false; // degrees preferred
    double Delta_phi = 0.0;      bool has_Delta_phi = false; // fallback in degrees

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
    }

    inline double phi_deg() const {
        if (has_phi2) return phi2;
        if (has_Delta_phi) return Delta_phi;
        return std::numeric_limits<double>::quiet_NaN();
    }

    inline bool readyForCuts() const {
        return has_t1 && has_open_angle_ep2 && has_pTmiss;
    }
    inline bool readyForVars() const {
        return has_x && has_Q2 && (has_phi2 || has_Delta_phi);
    }
};

// ------------------------------
// Bin grid assembled from CSV
// ------------------------------
struct RowDef {
    int    row_index = -1; // "bin index" from CSV
    double xBmin=0, xBmax=0;
    double Q2min=0, Q2max=0;
    double tminAbs=0, tmaxAbs=0;
    double phimin=0, phimax=0; // degrees

    // CSV column indices for write-back:
    int col_xBavg[8]  = {-1,-1,-1,-1,-1,-1,-1,-1};
    int col_Q2avg[8]  = {-1,-1,-1,-1,-1,-1,-1,-1};
    int col_tavg[8]   = {-1,-1,-1,-1,-1,-1,-1,-1};
    int col_phiavg[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
};

static const std::vector<std::string>& group_labels() {
    static const std::vector<std::string> labs = {
        "Fa18 Inb","Fa18 Out","Sp19 Inb","Sp18 Inb","Sp18 Out","Fa18","Sp18","10.6 GeV"
    };
    return labs;
}
enum GroupIdx {
    G_Fa18_Inb=0, G_Fa18_Out=1, G_Sp19_Inb=2, G_Sp18_Inb=3, G_Sp18_Out=4, G_Fa18=5, G_Sp18=6, G_10p6=7
};

// map canonical period key -> group index
static int period_to_group(const std::string& key) {
    if (key == "DVCS_Fa18_inb") return G_Fa18_Inb;
    if (key == "DVCS_Fa18_out") return G_Fa18_Out;
    if (key == "DVCS_Sp19_inb") return G_Sp19_Inb;
    if (key == "DVCS_Sp18_inb") return G_Sp18_Inb;
    if (key == "DVCS_Sp18_out") return G_Sp18_Out;
    return -1;
}

// ------------------------------
// Fast 4D bin indexer based on unique ranges
// ------------------------------
static std::vector<std::pair<double,double>> unique_ranges_from_rows(
    const std::vector<RowDef>& rows, char which)
{
    std::set<std::pair<double,double>> s;
    for (const auto& r : rows) {
        if (which=='x') s.emplace(r.xBmin, r.xBmax);
        else if (which=='Q') s.emplace(r.Q2min, r.Q2max);
        else if (which=='t') s.emplace(r.tminAbs, r.tmaxAbs);
        else if (which=='p') s.emplace(r.phimin, r.phimax);
    }
    return {s.begin(), s.end()};
}

static int find_idx(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < (int)ranges.size(); ++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

// map 4-tuple -> row slot
using Key4 = std::tuple<int,int,int,int>;

// ------------------------------
// Accumulator
// ------------------------------
struct Accum {
    double sum_x=0, sum_Q2=0, sum_t=0, sum_phi=0;
    long long n=0;
    inline void add(double xB, double Q2, double tabs, double phi_deg) {
        sum_x   += xB;
        sum_Q2  += Q2;
        sum_t   += tabs;
        sum_phi += wrap_deg_360(phi_deg);
        ++n;
    }
};

// ------------------------------
// Cut wrapper: BOTH global and 3-sigma, accept ANY topology
// ------------------------------
static bool passes_all_cuts(const std::string& period_key, const BranchBinder& b) {
    if (!global_cuts::passes(b.t1, b.open_angle_ep2, b.pTmiss)) return false;
    // accept if ANY topology's 3-sigma passes:
    static const Topology topos[] = {Topology::FD_FD, Topology::CD_FD, Topology::CD_FT};
    bool ok3 = false;
    for (Topology tp : topos) {
        if (passes3sigma(period_key, topo_to_csv_title(tp), b.t1, b.open_angle_ep2, b.pTmiss)) {
            ok3 = true; break;
        }
    }
    return ok3;
}

// ------------------------------
// Main
// ------------------------------
bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int /*worker_count*/)
{
    // 1) Load CSV
    std::ifstream fin(csv_path);
    if (!fin.is_open()) {
        std::cerr << "[bin_means] FATAL: cannot open CSV: " << csv_path << std::endl;
        return false;
    }
    std::string header_line;
    if (!std::getline(fin, header_line)) {
        std::cerr << "[bin_means] FATAL: empty CSV: " << csv_path << std::endl;
        return false;
    }
    std::vector<std::string> header = split_csv_line(header_line);
    auto H = build_header_index(header);

    // verify required bin-def columns
    const char* req[] = {"bin index","xBmin","xBmax","Q2min","Q2max","t_abs_min","t_abs_max","phimin","phimax"};
    for (auto* n : req) {
        if (!H.count(n)) { std::cerr << "[bin_means] FATAL: CSV missing column: " << n << std::endl; return false; }
    }

    // 2) Slurp all rows and build RowDef list
    std::vector<std::vector<std::string>> rows_csv;
    rows_csv.reserve(4000);
    std::string line;
    while (std::getline(fin, line)) {
        if (!line.empty()) rows_csv.push_back(split_csv_line(line));
    }
    fin.close();

    std::vector<RowDef> rows;
    rows.reserve(rows_csv.size());
    for (size_t r = 0; r < rows_csv.size(); ++r) {
        const auto& cols = rows_csv[r];
        RowDef R;
        // bin index
        int idx = -1;
        get_int(cols, H, "bin index", idx);
        R.row_index = idx;

        // edges
        get_double(cols, H, "xBmin", R.xBmin); get_double(cols, H, "xBmax", R.xBmax);
        get_double(cols, H, "Q2min", R.Q2min); get_double(cols, H, "Q2max", R.Q2max);
        get_double(cols, H, "t_abs_min", R.tminAbs); get_double(cols, H, "t_abs_max", R.tmaxAbs);
        get_double(cols, H, "phimin", R.phimin); get_double(cols, H, "phimax", R.phimax);

        // write-back columns
        const auto& groups = group_labels();
        for (int g = 0; g < (int)groups.size(); ++g) {
            // build exact header names
            std::string s1 = "xBavg, "     + groups[g];
            std::string s2 = "Q2avg, "     + groups[g];
            std::string s3 = "t_abs_avg, " + groups[g];
            std::string s4 = "phiavg, "    + groups[g];
            auto it1 = H.find(s1); auto it2 = H.find(s2);
            auto it3 = H.find(s3); auto it4 = H.find(s4);
            if (it1 == H.end() || it2 == H.end() || it3 == H.end() || it4 == H.end()) {
                std::cerr << "[bin_means] FATAL: column missing: " 
                          << (it1==H.end()?s1: it2==H.end()?s2: it3==H.end()?s3:s4) << std::endl;
                return false;
            }
            R.col_xBavg[g]  = it1->second;
            R.col_Q2avg[g]  = it2->second;
            R.col_tavg[g]   = it3->second;
            R.col_phiavg[g] = it4->second;
        }
        rows.push_back(R);
    }

    // 3) Build fast 4D index
    const auto x_ranges  = unique_ranges_from_rows(rows, 'x');
    const auto Q_ranges  = unique_ranges_from_rows(rows, 'Q');
    const auto t_ranges  = unique_ranges_from_rows(rows, 't');
    const auto p_ranges  = unique_ranges_from_rows(rows, 'p');

    std::unordered_map<Key4, int> lut; // -> row slot
    lut.reserve(rows.size()*2);
    for (int i = 0; i < (int)rows.size(); ++i) {
        int ix = find_idx(rows[i].xBmin, x_ranges);
        int iQ = find_idx(rows[i].Q2min, Q_ranges);
        int it = find_idx(rows[i].tminAbs, t_ranges);
        int ip = find_idx(rows[i].phimin,  p_ranges);
        if (ix<0 || iQ<0 || it<0 || ip<0) {
            std::cerr << "[bin_means] WARNING: could not index row i="<<i<<" bin="<<rows[i].row_index<<"\n";
            continue;
        }
        lut.emplace(Key4{ix,iQ,it,ip}, i);
    }

    // 4) Per-period accumulators (store to compute combined later)
    const int R = (int)rows.size();
    std::vector<std::vector<Accum>> per_group_acc(5, std::vector<Accum>(R)); // 5 periods only
    auto add_event = [&](int gidx, int row_slot, double xB, double Q2, double tabs, double phi){
        if (gidx<0 || gidx>=5) return;
        if (row_slot<0 || row_slot>=R) return;
        per_group_acc[gidx][row_slot].add(xB,Q2,tabs,phi);
    };

    // 5) Loop each period ONCE (single-thread for ROOT safety)
    for (const auto& kv : dataTrees) {
        const std::string period_key = kv.first;
        int gidx = period_to_group(period_key); // only 5 period groups here
        if (gidx < 0) continue;                 // ignore other keys

        TTree* t = kv.second;
        if (!t) { std::cerr << "[bin_means] WARNING: null tree for " << period_key << "\n"; continue; }

        BranchBinder b; b.bind(t);
        if (!b.readyForCuts() || !b.readyForVars()) {
            std::cerr << "[bin_means] FATAL: tree '"<<period_key<<"' missing required branches.\n";
            return false;
        }

        const Long64_t nent = t->GetEntries();
        for (Long64_t i = 0; i < nent; ++i) {
            t->GetEntry(i);
            if (!passes_all_cuts(period_key, b)) continue;

            double xB  = b.x;
            double Q2  = b.Q2;
            double tab = std::fabs(b.t1);
            double phi = wrap_deg_360(b.phi_deg());
            if (!std::isfinite(xB) || !std::isfinite(Q2) || !std::isfinite(tab) || !std::isfinite(phi)) continue;

            int ix = find_idx(xB,  x_ranges);
            int iQ = find_idx(Q2,  Q_ranges);
            int it = find_idx(tab, t_ranges);
            int ip = find_idx(phi, p_ranges);
            if (ix<0 || iQ<0 || it<0 || ip<0) continue;

            auto itL = lut.find(Key4{ix,iQ,it,ip});
            if (itL == lut.end()) continue;
            int row_slot = itL->second;

            add_event(gidx, row_slot, xB, Q2, tab, phi);
        }
    }

    // 6) Write per-period means into CSV rows
    auto write_means_for_group = [&](int gidx_period, int g_write_slot) {
        for (int i = 0; i < R; ++i) {
            const auto& a = per_group_acc[gidx_period][i];
            if (a.n <= 0) continue;
            rows_csv[i][ rows[i].col_xBavg[g_write_slot]  ] = (static_cast<std::ostringstream&>(std::ostringstream() << std::setprecision(8) << std::fixed << (a.sum_x / a.n))).str();
            rows_csv[i][ rows[i].col_Q2avg[g_write_slot]  ] = (static_cast<std::ostringstream&>(std::ostringstream() << std::setprecision(8) << std::fixed << (a.sum_Q2 / a.n))).str();
            rows_csv[i][ rows[i].col_tavg[g_write_slot]   ] = (static_cast<std::ostringstream&>(std::ostringstream() << std::setprecision(8) << std::fixed << (a.sum_t / a.n))).str();
            rows_csv[i][ rows[i].col_phiavg[g_write_slot] ] = (static_cast<std::ostringstream&>(std::ostringstream() << std::setprecision(8) << std::fixed << (a.sum_phi / a.n))).str();
        }
    };

    // per-period (exact slot == group index)
    write_means_for_group(G_Fa18_Inb, G_Fa18_Inb);
    write_means_for_group(G_Fa18_Out, G_Fa18_Out);
    write_means_for_group(G_Sp19_Inb, G_Sp19_Inb);
    write_means_for_group(G_Sp18_Inb, G_Sp18_Inb);
    write_means_for_group(G_Sp18_Out, G_Sp18_Out);

    // 7) Combined groups built from per-period sums
    auto combine_and_write = [&](const std::vector<int>& src_gidx, int dst_slot) {
        std::vector<Accum> combined(R);
        for (int i = 0; i < R; ++i) {
            for (int g : src_gidx) {
                combined[i].sum_x   += per_group_acc[g][i].sum_x;
                combined[i].sum_Q2  += per_group_acc[g][i].sum_Q2;
                combined[i].sum_t   += per_group_acc[g][i].sum_t;
                combined[i].sum_phi += per_group_acc[g][i].sum_phi;
                combined[i].n       += per_group_acc[g][i].n;
            }
            if (combined[i].n > 0) {
                rows_csv[i][ rows[i].col_xBavg[dst_slot]  ] = (static_cast<std::ostringstream&>(std::ostringstream() << std::setprecision(8) << std::fixed << (combined[i].sum_x   / combined[i].n))).str();
                rows_csv[i][ rows[i].col_Q2avg[dst_slot]  ] = (static_cast<std::ostringstream&>(std::ostringstream() << std::setprecision(8) << std::fixed << (combined[i].sum_Q2  / combined[i].n))).str();
                rows_csv[i][ rows[i].col_tavg[dst_slot]   ] = (static_cast<std::ostringstream&>(std::ostringstream() << std::setprecision(8) << std::fixed << (combined[i].sum_t   / combined[i].n))).str();
                rows_csv[i][ rows[i].col_phiavg[dst_slot] ] = (static_cast<std::ostringstream&>(std::ostringstream() << std::setprecision(8) << std::fixed << (combined[i].sum_phi / combined[i].n))).str();
            }
        }
    };

    // Fa18 = Fa18 Inb + Fa18 Out
    combine_and_write({G_Fa18_Inb, G_Fa18_Out}, G_Fa18);

    // Sp18 = Sp18 Inb + Sp18 Out
    combine_and_write({G_Sp18_Inb, G_Sp18_Out}, G_Sp18);

    // 10.6 GeV = Fa18 Inb + Fa18 Out + Sp18 Inb + Sp18 Out  (exclude Sp19 Inb)
    combine_and_write({G_Fa18_Inb, G_Fa18_Out, G_Sp18_Inb, G_Sp18_Out}, G_10p6);

    // 8) Write CSV back
    std::ofstream fout(csv_path);
    if (!fout.is_open()) {
        std::cerr << "[bin_means] FATAL: cannot re-open CSV for write: " << csv_path << std::endl;
        return false;
    }
    fout << join_csv_row(header) << "\n";
    for (const auto& row : rows_csv) fout << join_csv_row(row) << "\n";
    fout.close();

    std::cout << "[bin_means] Updated bin means in: " << csv_path << std::endl;
    return true;
}