// acceptance.cpp
// DVCS acceptance from non-radiative MC:
//   - acceptance, <period> = N_rec / N_gen per (xB, Q2, |t|, phi) bin
//   - N_rec from reconstructed MC with:
//       * hard global DVCS exclusivity cuts (t1, open_angle_ep2, pTmiss),
//       * global run-blacklist and global cuts from global_cuts.h,
//       * AND 3-sigma cuts loaded from combined_cuts.json, applied
//         topology-by-topology (FD_FD, CD_FD, CD_FT) via detector1/2.
//   - N_gen from generated MC without exclusivity cuts
//   - Periods: Fa18 Inb, Fa18 Out, Sp19 Inb, Sp18 Inb, Sp18 Out
//   - Acceptance is stored as "(value, stat, sys)" triple
//     with stat the binomial uncertainty on N_rec / N_gen and
//     sys = 0 for now.
//   - Produces per-period acceptance vs phi plots in the usual
//     xB by (Q2, |t|) canvas layout under output/acceptance/<PeriodDir>/.

#include "acceptance.h"
#include "global_cuts.h"   // <- reuse the same GLOBAL cuts as total_counts.cpp

#include <nlohmann/json.hpp>

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TError.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cstdio>

namespace {

using nlohmann::json;

// ---------------- CSV helper ----------------

struct CsvDoc {
    std::vector<std::string> header;
    std::unordered_map<std::string,int> index;
    std::vector<std::vector<std::string>> rows;

    static std::vector<std::string> split_csv_line(const std::string& line) {
        std::vector<std::string> out;
        std::string cur;
        bool inq = false;
        for (char c : line) {
            if (c == '"') {
                inq = !inq;
                continue;
            }
            if (c == ',' && !inq) {
                out.push_back(cur);
                cur.clear();
            } else {
                cur.push_back(c);
            }
        }
        out.push_back(cur);
        return out;
    }

    static void write_field(std::ostream& os, const std::string& s) {
        bool needq = s.find(',') != std::string::npos || s.find('"') != std::string::npos;
        if (!needq) {
            os << s;
            return;
        }
        os << '"';
        for (char ch : s) {
            if (ch == '"') {
                os << "\"\"";
            } else {
                os << ch;
            }
        }
        os << '"';
    }

    static double to_double(const std::string& s) {
        if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
        char* e = nullptr;
        double v = std::strtod(s.c_str(), &e);
        if (e == s.c_str()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return v;
    }

    bool load(const std::string& path) {
        std::ifstream fin(path);
        if (!fin.is_open()) {
            std::cerr << "[acceptance] ERROR: cannot open CSV: " << path << "\n";
            return false;
        }
        std::string line;
        if (!std::getline(fin, line)) {
            std::cerr << "[acceptance] ERROR: empty CSV: " << path << "\n";
            return false;
        }
        header = split_csv_line(line);
        index.clear();
        for (int i = 0; i < (int)header.size(); ++i) {
            index[header[i]] = i;
        }
        rows.clear();
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            rows.push_back(split_csv_line(line));
        }
        for (auto& r : rows) {
            r.resize(header.size());
        }
        return true;
    }

    bool save_atomic(const std::string& path) const {
        const std::string tmp = path + ".tmp";
        {
            std::ofstream fout(tmp);
            if (!fout.is_open()) {
                std::cerr << "[acceptance] ERROR: cannot write CSV tmp: " << tmp << "\n";
                return false;
            }
            for (size_t i = 0; i < header.size(); ++i) {
                write_field(fout, header[i]);
                if (i + 1 < header.size()) fout << ',';
            }
            fout << "\n";
            for (const auto& row : rows) {
                for (size_t i = 0; i < row.size(); ++i) {
                    write_field(fout, row[i]);
                    if (i + 1 < row.size()) fout << ',';
                }
                fout << "\n";
            }
        }
        std::error_code ec;
        std::filesystem::rename(tmp, path, ec);
        if (ec) {
            std::remove(path.c_str());
            std::filesystem::rename(tmp, path, ec);
            if (ec) {
                std::cerr << "[acceptance] ERROR: atomic rename failed ("
                          << ec.message() << ")\n";
                return false;
            }
        }
        return true;
    }

    int nrows() const { return (int)rows.size(); }

    int col_index(const std::string& name) const {
        auto it = index.find(name);
        return (it == index.end()) ? -1 : it->second;
    }

    double as_double(int r, int c) const {
        if (r < 0 || r >= nrows()) return std::numeric_limits<double>::quiet_NaN();
        if (c < 0 || c >= (int)header.size()) return std::numeric_limits<double>::quiet_NaN();
        return to_double(rows[r][c]);
    }
};

// ------------- basic helpers -------------

static inline double PI() { return 3.14159265358979323846; }
static inline double RAD2DEG(double r) { return r * 180.0 / PI(); }

static inline double wrap_deg(double phi) {
    double x = std::fmod(phi, 360.0);
    if (x < 0.0) x += 360.0;
    return x;
}

static inline std::string to_lower_nospace(std::string s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == ' ' || c == '\t' || c == '_') continue;
        out.push_back((char)std::tolower((unsigned char)c));
    }
    return out;
}

static std::string canonical_period_dir(const std::string& L) {
    const std::string k = to_lower_nospace(L);
    if (k == "fa18inb") return "Fa18_Inb";
    if (k == "fa18out") return "Fa18_Out";
    if (k == "fa18inbsupp") return "Fa18_Inb_Supp";
    if (k == "sp18inb") return "Sp18_Inb";
    if (k == "sp18out") return "Sp18_Out";
    if (k == "sp19inb") return "Sp19_Inb";
    std::string s = L;
    for (char& c : s) {
        if (c == ' ') c = '_';
    }
    return s;
}

static inline std::string period_dir_for_label(const std::string& L) {
    return canonical_period_dir(L);
}

static void ensure_dir(const std::string& p) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (!fs::exists(p)) {
        fs::create_directories(p, ec);
        if (ec) {
            std::cerr << "[acceptance] FATAL: could not create directory: "
                      << p << " (" << ec.message() << ")\n";
            std::exit(EXIT_FAILURE);
        }
    }
}

static double safe_mean(const std::vector<double>& v) {
    double s = 0.0;
    int n = 0;
    for (double x : v) {
        if (std::isfinite(x)) {
            s += x;
            ++n;
        }
    }
    return n ? s / n : std::numeric_limits<double>::quiet_NaN();
}

// ------------- triple formatting / parsing "(value, stat, sys)" -------------

static std::string format_triple(double v, double s_stat, double s_sys) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << "("
        << std::setprecision(8) << v << ", "
        << std::setprecision(8) << s_stat << ", "
        << std::setprecision(8) << s_sys
        << ")";
    return oss.str();
}

static bool parse_triple(const std::string& s,
                         double& value,
                         double& stat,
                         double& sys)
{
    value = std::numeric_limits<double>::quiet_NaN();
    stat  = std::numeric_limits<double>::quiet_NaN();
    sys   = std::numeric_limits<double>::quiet_NaN();

    std::string trimmed;
    trimmed.reserve(s.size());
    for (char c : s) {
        if (!std::isspace((unsigned char)c)) trimmed.push_back(c);
    }
    if (trimmed.empty()) return false;

    if (trimmed.front() == '(' && trimmed.back() == ')') {
        trimmed = trimmed.substr(1, trimmed.size() - 2);
    }

    std::vector<std::string> parts;
    std::string cur;
    for (char c : trimmed) {
        if (c == ',') {
            parts.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    parts.push_back(cur);

    if (parts.size() != 3) {
        return false;
    }

    char* e1 = nullptr;
    char* e2 = nullptr;
    char* e3 = nullptr;

    value = std::strtod(parts[0].c_str(), &e1);
    stat  = std::strtod(parts[1].c_str(), &e2);
    sys   = std::strtod(parts[2].c_str(), &e3);

    if (e1 == parts[0].c_str()) return false;
    if (e2 == parts[1].c_str()) return false;
    if (e3 == parts[2].c_str()) return false;

    return true;
}

// ------------- periods and tag mapping -------------

static const std::vector<std::string> kPeriods = {
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out"
};

struct McTags {
    std::string genTag;
    std::string recTag;
};

static std::map<std::string, McTags> build_mc_tag_map() {
    std::map<std::string, McTags> m;
    m["Fa18 Inb"] = {"DVCS_Fa18_inb_gen", "DVCS_Fa18_inb_rec"};
    m["Fa18 Out"] = {"DVCS_Fa18_out_gen", "DVCS_Fa18_out_rec"};
    m["Sp19 Inb"] = {"DVCS_Sp19_inb_gen", "DVCS_Sp19_inb_rec"};
    m["Sp18 Inb"] = {"DVCS_Sp18_inb_gen", "DVCS_Sp18_inb_rec"};
    m["Sp18 Out"] = {"DVCS_Sp18_out_gen", "DVCS_Sp18_out_rec"};
    return m;
}

// ----------------- 3-sigma cuts from JSON (topology-dependent) -----------------

struct SigmaCut {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std  = std::numeric_limits<double>::quiet_NaN();
};

using VarCutMap  = std::unordered_map<std::string, SigmaCut>;      // var -> {mean,std}
using TopoCutMap = std::unordered_map<std::string, VarCutMap>;     // "DVCS_<PeriodDir>_<TopoDir>" -> VarCutMap

static inline bool within_3sigma(double val, const SigmaCut& sc) {
    if (!std::isfinite(val) || !std::isfinite(sc.mean) || !std::isfinite(sc.std) || sc.std <= 0.0) {
        return true;
    }
    return std::fabs(val - sc.mean) <= 3.0 * sc.std;
}

static TopoCutMap load_sigma_cuts_data(const std::string& path) {
    TopoCutMap out;

    if (path.empty()) {
        std::cerr << "[acceptance] FATAL: combined cuts JSON path is empty.\n";
        std::exit(EXIT_FAILURE);
    }

    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::cerr << "[acceptance] FATAL: could not open combined cuts JSON: " << path << "\n";
        std::exit(EXIT_FAILURE);
    }

    json j;
    try {
        fin >> j;
    } catch (const std::exception& e) {
        std::cerr << "[acceptance] FATAL: failed to parse combined cuts JSON ("
                  << e.what() << ")\n";
        std::exit(EXIT_FAILURE);
    }

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key(); // e.g. "DVCS_Fa18_Inb_FD_FD" or maybe pi0 keys
        const json& block = it.value();

        if (!block.contains("data") || !block["data"].is_object()) continue;

        VarCutMap m;
        for (auto vit = block["data"].begin(); vit != block["data"].end(); ++vit) {
            const std::string vname = vit.key(); // Emiss2, Mx2, ...
            const json& vs = vit.value();
            SigmaCut sc;
            if (vs.contains("mean") && vs["mean"].is_number()) sc.mean = vs["mean"].get<double>();
            if (vs.contains("std")  && vs["std"].is_number())  sc.std  = vs["std"].get<double>();
            if (std::isfinite(sc.std) && sc.std > 0.0) {
                m[vname] = sc;
            }
        }
        if (!m.empty()) out[key] = std::move(m);
    }

    std::cout << "[acceptance] Loaded sigma cuts for " << out.size()
              << " topology keys from " << path << "\n";
    return out;
}

// ---------------- topology (FD_FD, CD_FD, CD_FT) ----------------

enum class Topology { FD_FD = 0, CD_FD = 1, CD_FT = 2 };

static inline const char* topo_dir(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }
    return "FD_FD";
}

struct TopologyResolver {
    int detector1 = 0;
    int detector2 = 0;
    bool have1 = false;
    bool have2 = false;

    void enable_and_bind(TTree* t, const std::string& period_label) {
        t->SetBranchStatus("*", 1);

        have1 = (t->GetBranch("detector1") != nullptr);
        have2 = (t->GetBranch("detector2") != nullptr);
        if (!(have1 && have2)) {
            std::cerr << "[acceptance] FATAL: Missing detector1/detector2 in rec MC tree for period "
                      << period_label << ".\n";
            std::exit(EXIT_FAILURE);
        }

        t->SetBranchAddress("detector1", &detector1);
        t->SetBranchAddress("detector2", &detector2);
    }

    int index() const {
        if (detector1 == 1 && detector2 == 1) return (int)Topology::FD_FD;
        if (detector1 == 2 && detector2 == 1) return (int)Topology::CD_FD;
        if (detector1 == 2 && detector2 == 0) return (int)Topology::CD_FT;
        return -1;
    }
};

// ------------- bin specification per CSV row -------------

struct RowBin {
    double xbmin, xbmax;
    double q2min, q2max;
    double tmin,  tmax;
    double phimin, phimax;
    int    row_index;   // index into CsvDoc::rows
};

// We only build RowBin entries for "standard" phi bins,
// NOT for phi-integrated bins (phi width >= 360 degrees).
static std::vector<RowBin> build_row_bins(const CsvDoc& csv) {
    const int c_xb_min  = csv.col_index("xBmin");
    const int c_xb_max  = csv.col_index("xBmax");
    const int c_q2_min  = csv.col_index("Q2min");
    const int c_q2_max  = csv.col_index("Q2max");
    const int c_tab_min = csv.col_index("t_abs_min");
    const int c_tab_max = csv.col_index("t_abs_max");
    const int c_phi_min = csv.col_index("phimin");
    const int c_phi_max = csv.col_index("phimax");

    if (c_xb_min < 0 || c_xb_max < 0 ||
        c_q2_min < 0 || c_q2_max < 0 ||
        c_tab_min < 0 || c_tab_max < 0 ||
        c_phi_min < 0 || c_phi_max < 0) {
        std::cerr << "[acceptance] FATAL: missing one or more bin-edge columns "
                  << "(xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax)\n";
        std::exit(EXIT_FAILURE);
    }

    std::vector<RowBin> bins;
    const int NR = csv.nrows();
    for (int r = 0; r < NR; ++r) {
        RowBin b;
        b.xbmin  = csv.as_double(r, c_xb_min);
        b.xbmax  = csv.as_double(r, c_xb_max);
        b.q2min  = csv.as_double(r, c_q2_min);
        b.q2max  = csv.as_double(r, c_q2_max);
        b.tmin   = csv.as_double(r, c_tab_min);
        b.tmax   = csv.as_double(r, c_tab_max);
        b.phimin = csv.as_double(r, c_phi_min);
        b.phimax = csv.as_double(r, c_phi_max);
        b.row_index = r;

        const double phi_width = b.phimax - b.phimin;
        if (!std::isfinite(phi_width)) {
            continue;
        }

        if (std::fabs(phi_width) >= 359.0) {
            continue;
        }

        if (!(b.xbmax > b.xbmin &&
              b.q2max > b.q2min &&
              b.tmax  > b.tmin  &&
              b.phimax > b.phimin)) {
            continue;
        }

        bins.push_back(b);
    }

    std::cout << "[acceptance] Built " << bins.size()
              << " RowBin entries out of " << NR << " CSV rows "
              << "(skipped phi-integrated or invalid bins).\n";

    if (bins.empty()) {
        std::cerr << "[acceptance] FATAL: no usable bins found in CSV.\n";
        std::exit(EXIT_FAILURE);
    }

    return bins;
}

// For each period, mark which CSV rows actually have data for that period.
static std::map<std::string, std::vector<bool>>
build_row_has_data(const CsvDoc& csv)
{
    const int NR = csv.nrows();
    std::map<std::string, std::vector<bool>> has_data;

    for (const auto& per : kPeriods) {
        std::string cname = "xBavg, " + per;
        int c = csv.col_index(cname);
        if (c < 0) {
            std::cerr << "[acceptance] FATAL: missing column '" << cname
                      << "' needed to decide where to compute acceptance.\n";
            std::exit(EXIT_FAILURE);
        }

        std::vector<bool> flags(NR, false);
        for (int r = 0; r < NR; ++r) {
            const std::string& cell = csv.rows[r][c];
            if (cell.empty()) {
                flags[r] = false;
            } else {
                double v = CsvDoc::to_double(cell);
                flags[r] = std::isfinite(v);
            }
        }
        has_data[per] = flags;
    }

    return has_data;
}

// Find which CSV row an event belongs to, given xB, Q2, |t|, phiDeg.
static int find_row_for_event(double xB,
                              double Q2,
                              double tAbs,
                              double phiDeg,
                              const std::vector<RowBin>& bins)
{
    for (const auto& b : bins) {
        if (!(xB >= b.xbmin && xB < b.xbmax)) continue;
        if (!(Q2 >= b.q2min && Q2 < b.q2max)) continue;
        if (!(tAbs >= b.tmin && tAbs < b.tmax)) continue;
        if (!(phiDeg >= b.phimin && phiDeg < b.phimax)) continue;
        return b.row_index;
    }
    return -1;
}

// ------------- exclusivity for reconstructed MC (global + 3-sigma) -------------

static bool rec_passes_exclusivity(double t1,
                                   double open_angle_ep2_deg,
                                   int runnum,
                                   const std::map<std::string,double>& recVars,
                                   const VarCutMap* cutsForTopo)
{
    auto itP = recVars.find("pTmiss");
    if (itP == recVars.end()) {
        std::cerr << "[acceptance] FATAL: internal error: pTmiss not found in recVars map.\n";
        std::exit(EXIT_FAILURE);
    }
    const double pTmiss = itP->second;

    // Global run blacklist + global cuts (same as total_counts.cpp)
    if (is_excluded_run(runnum)) return false;
    if (!passes_global_cuts(t1, open_angle_ep2_deg, pTmiss)) return false;

    // Topology-dependent 3-sigma cuts from combined_cuts.json
    if (cutsForTopo) {
        for (const auto& kv : *cutsForTopo) {
            const std::string& vname = kv.first;
            const SigmaCut&    sc    = kv.second;

            auto iv = recVars.find(vname);
            if (iv == recVars.end()) {
                std::cerr << "[acceptance] FATAL: missing branch value for variable '"
                          << vname << "' while applying 3-sigma cuts.\n";
                std::exit(EXIT_FAILURE);
            }
            const double val = iv->second;
            if (!within_3sigma(val, sc)) return false;
        }
    }

    return true;
}

// Build the set of branch names we actually need for the 3-sigma cuts
// plus pTmiss for the global exclusivity.
static std::set<std::string>
build_cut_branches_for_period(const std::string& period_label,
                              const TopoCutMap& sigmaCuts)
{
    std::set<std::string> vars;
    const std::string period_dir = period_dir_for_label(period_label);

    // Topology-dependent 3-sigma variables
    for (int topo_idx = 0; topo_idx <= 2; ++topo_idx) {
        const Topology T = (Topology)topo_idx;
        const std::string key =
            std::string("DVCS_") + period_dir + "_" + topo_dir(T);

        auto it = sigmaCuts.find(key);
        if (it == sigmaCuts.end()) continue;

        for (const auto& kv : it->second) {
            vars.insert(kv.first);
        }
    }

    // Always include pTmiss for the hard global cut
    vars.insert("pTmiss");
    return vars;
}

// ------------- MC counting per period -------------

static void accumulate_counts_for_period(const std::string& period_label,
                                         const CsvDoc& csv,
                                         const std::vector<RowBin>& bins,
                                         const std::vector<bool>& row_has_data,
                                         TTree* genTree,
                                         TTree* recTree,
                                         const TopoCutMap& sigmaCuts,
                                         std::vector<double>& gen_counts,
                                         std::vector<double>& rec_counts)
{
    if (!genTree || !recTree) {
        std::cerr << "[acceptance] FATAL: null TTree pointer for period "
                  << period_label << ".\n";
        std::exit(EXIT_FAILURE);
    }

    const int NR = csv.nrows();
    gen_counts.assign(NR, 0.0);
    rec_counts.assign(NR, 0.0);

    const char* br_x      = "x";
    const char* br_Q2     = "Q2";
    const char* br_t1     = "t1";
    const char* br_phi2   = "phi2";
    const char* br_oa_ep2 = "open_angle_ep2";
    const char* br_runnum = "runnum";

    // Required branches in generated MC
    if (!genTree->GetBranch(br_x) ||
        !genTree->GetBranch(br_Q2) ||
        !genTree->GetBranch(br_t1) ||
        !genTree->GetBranch(br_phi2)) {
        std::cerr << "[acceptance] FATAL: missing one or more branches in genTree for period "
                  << period_label << " (expected: x, Q2, t1, phi2).\n";
        std::exit(EXIT_FAILURE);
    }

    // Required branches in reconstructed MC (kinematics + topology + exclusivity vars + runnum)
    if (!recTree->GetBranch(br_x) ||
        !recTree->GetBranch(br_Q2) ||
        !recTree->GetBranch(br_t1) ||
        !recTree->GetBranch(br_phi2) ||
        !recTree->GetBranch(br_oa_ep2) ||
        !recTree->GetBranch(br_runnum)) {
        std::cerr << "[acceptance] FATAL: missing one or more kinematic/exclusivity/run branches in recTree for period "
                  << period_label
                  << " (expected: x, Q2, t1, phi2, open_angle_ep2, runnum).\n";
        std::exit(EXIT_FAILURE);
    }

    // Bind topology resolver (detector1/detector2) for rec MC
    TopologyResolver topo;
    topo.enable_and_bind(recTree, period_label);

    // Bind generated MC branches
    double g_x    = 0.0;
    double g_Q2   = 0.0;
    double g_t1   = 0.0;
    double g_phi2 = 0.0;

    genTree->SetBranchAddress(br_x,    &g_x);
    genTree->SetBranchAddress(br_Q2,   &g_Q2);
    genTree->SetBranchAddress(br_t1,   &g_t1);
    genTree->SetBranchAddress(br_phi2, &g_phi2);

    const Long64_t Ngen = genTree->GetEntries();
    Long64_t used_gen = 0;

    int next_pct_gen = 10;

    for (Long64_t i = 0; i < Ngen; ++i) {
        genTree->GetEntry(i);

        const double xB   = g_x;
        const double Q2   = g_Q2;
        const double tAbs = std::fabs(g_t1);
        const double phiD = wrap_deg(RAD2DEG(g_phi2));

        int row = find_row_for_event(xB, Q2, tAbs, phiD, bins);
        if (row < 0 || row >= NR) continue;
        if (!row_has_data[row]) continue;

        gen_counts[row] += 1.0;
        ++used_gen;

        if (Ngen > 0 && next_pct_gen <= 100) {
            double pct = 100.0 * (double)(i + 1) / (double)Ngen;
            while (pct >= (double)next_pct_gen && next_pct_gen <= 100) {
                std::cout << "[acceptance] Period " << period_label
                          << " gen progress: " << (double)next_pct_gen << "% ("
                          << (long long)(i + 1) << "/" << (long long)Ngen << ")\n";
                next_pct_gen += 10;
            }
        }
    }

    std::cout << "[acceptance] Period " << period_label
              << " gen MC: total entries = " << Ngen
              << " ; binned (with period-data flag) = " << used_gen << "\n";

    // Bind reconstructed MC branches
    double r_x      = 0.0;
    double r_Q2     = 0.0;
    double r_t1     = 0.0;
    double r_phi2   = 0.0;
    double r_oa     = 0.0;
    int    r_run    = 0;

    recTree->SetBranchAddress(br_x,      &r_x);
    recTree->SetBranchAddress(br_Q2,     &r_Q2);
    recTree->SetBranchAddress(br_t1,     &r_t1);
    recTree->SetBranchAddress(br_phi2,   &r_phi2);
    recTree->SetBranchAddress(br_oa_ep2, &r_oa);
    recTree->SetBranchAddress(br_runnum, &r_run);

    // Build the set of branches needed for 3-sigma cuts + pTmiss
    const std::set<std::string> cutBranches =
        build_cut_branches_for_period(period_label, sigmaCuts);

    // Map for recVars
    std::map<std::string,double> recVars;
    recVars.clear();
    for (const auto& vname : cutBranches) {
        if (!recTree->GetBranch(vname.c_str())) {
            std::cerr << "[acceptance] FATAL: recTree for period " << period_label
                      << " is missing branch '" << vname
                      << "' required for exclusivity / 3-sigma cuts.\n";
            std::exit(EXIT_FAILURE);
        }
        recVars[vname] = 0.0;
        recTree->SetBranchAddress(vname.c_str(), &recVars[vname]);
    }

    const Long64_t Nrec = recTree->GetEntries();
    Long64_t used_rec = 0;
    Long64_t passed_excl = 0;

    int next_pct_rec = 10;

    const std::string period_dir = period_dir_for_label(period_label);

    for (Long64_t i = 0; i < Nrec; ++i) {
        recTree->GetEntry(i);

        const int topo_idx = topo.index();
        if (topo_idx < 0 || topo_idx > 2) {
            if (Nrec > 0 && next_pct_rec <= 100) {
                double pct = 100.0 * (double)(i + 1) / (double)Nrec;
                while (pct >= (double)next_pct_rec && next_pct_rec <= 100) {
                    std::cout << "[acceptance] Period " << period_label
                              << " rec progress: " << (double)next_pct_rec << "% ("
                              << (long long)(i + 1) << "/" << (long long)Nrec << ")\n";
                    next_pct_rec += 10;
                }
            }
            continue;
        }
        const Topology T = (Topology)topo_idx;
        const std::string topo_key =
            std::string("DVCS_") + period_dir + "_" + topo_dir(T);

        const VarCutMap* cutsForTopo = nullptr;
        auto itTopoCuts = sigmaCuts.find(topo_key);
        if (itTopoCuts != sigmaCuts.end()) {
            cutsForTopo = &itTopoCuts->second;
        }

        if (!rec_passes_exclusivity(r_t1,
                                    r_oa,
                                    r_run,
                                    recVars,
                                    cutsForTopo)) {
            if (Nrec > 0 && next_pct_rec <= 100) {
                double pct = 100.0 * (double)(i + 1) / (double)Nrec;
                while (pct >= (double)next_pct_rec && next_pct_rec <= 100) {
                    std::cout << "[acceptance] Period " << period_label
                              << " rec progress: " << (double)next_pct_rec << "% ("
                              << (long long)(i + 1) << "/" << (long long)Nrec << ")\n";
                    next_pct_rec += 10;
                }
            }
            continue;
        }
        ++passed_excl;

        const double xB   = r_x;
        const double Q2   = r_Q2;
        const double tAbs = std::fabs(r_t1);
        const double phiD = wrap_deg(RAD2DEG(r_phi2));

        int row = find_row_for_event(xB, Q2, tAbs, phiD, bins);
        if (row < 0 || row >= NR) {
            if (Nrec > 0 && next_pct_rec <= 100) {
                double pct = 100.0 * (double)(i + 1) / (double)Nrec;
                while (pct >= (double)next_pct_rec && next_pct_rec <= 100) {
                    std::cout << "[acceptance] Period " << period_label
                              << " rec progress: " << (double)next_pct_rec << "% ("
                              << (long long)(i + 1) << "/" << (long long)Nrec << ")\n";
                    next_pct_rec += 10;
                }
            }
            continue;
        }

        if (!row_has_data[row]) {
            if (Nrec > 0 && next_pct_rec <= 100) {
                double pct = 100.0 * (double)(i + 1) / (double)Nrec;
                while (pct >= (double)next_pct_rec && next_pct_rec <= 100) {
                    std::cout << "[acceptance] Period " << period_label
                              << " rec progress: " << (double)next_pct_rec << "% ("
                              << (long long)(i + 1) << "/" << (long long)Nrec << ")\n";
                    next_pct_rec += 10;
                }
            }
            continue;
        }

        rec_counts[row] += 1.0;
        ++used_rec;

        if (Nrec > 0 && next_pct_rec <= 100) {
            double pct = 100.0 * (double)(i + 1) / (double)Nrec;
            while (pct >= (double)next_pct_rec && next_pct_rec <= 100) {
                std::cout << "[acceptance] Period " << period_label
                          << " rec progress: " << (double)next_pct_rec << "% ("
                          << (long long)(i + 1) << "/" << (long long)Nrec << ")\n";
                next_pct_rec += 10;
            }
        }
    }

    std::cout << "[acceptance] Period " << period_label
              << " rec MC: total entries = " << Nrec
              << " ; passed global+3sigma exclusivity = " << passed_excl
              << " ; binned (with period-data flag) = " << used_rec << "\n";
}

// ------------- CSV filling (triples) -------------

static bool fill_acceptance_columns(CsvDoc& csv,
                                    const std::map<std::string, std::vector<double>>& gen_all,
                                    const std::map<std::string, std::vector<double>>& rec_all,
                                    const std::map<std::string, std::vector<bool>>& row_has_data)
{
    const int NR = csv.nrows();

    std::map<std::string,int> acc_idx;
    for (const auto& per : kPeriods) {
        const std::string cname = "acceptance, " + per;
        int idx = csv.col_index(cname);
        if (idx < 0) {
            std::cerr << "[acceptance] FATAL: missing acceptance column '"
                      << cname << "' in CSV header.\n";
            return false;
        }
        acc_idx[per] = idx;
    }

    std::size_t cells_written = 0;

    for (const auto& per : kPeriods) {
        auto itG = gen_all.find(per);
        auto itR = rec_all.find(per);
        auto itD = row_has_data.find(per);

        if (itG == gen_all.end() || itR == rec_all.end() || itD == row_has_data.end()) {
            std::cerr << "[acceptance] FATAL: internal error, missing counts or flags for period "
                      << per << ".\n";
            return false;
        }

        const std::vector<double>& gen = itG->second;
        const std::vector<double>& rec = itR->second;
        const std::vector<bool>&   has = itD->second;

        if ((int)gen.size() != NR || (int)rec.size() != NR || (int)has.size() != NR) {
            std::cerr << "[acceptance] FATAL: size mismatch for counts/flags vectors in period "
                      << per << ".\n";
            return false;
        }

        const int c_acc = acc_idx[per];

        double period_gen_sum = 0.0;
        double period_rec_sum = 0.0;

        for (int r = 0; r < NR; ++r) {
            if (!has[r]) {
                continue;
            }

            const double Ng = gen[r];
            const double Nr = rec[r];

            period_gen_sum += Ng;
            period_rec_sum += Nr;

            if (Ng <= 0.0) {
                continue;
            }

            const double acc = (Nr > 0.0) ? (Nr / Ng) : 0.0;

            const double var = acc * (1.0 - acc) / Ng;
            const double acc_stat = (var > 0.0) ? std::sqrt(var) : 0.0;

            csv.rows[r][c_acc] = format_triple(acc, acc_stat, 0.0);
            ++cells_written;
        }

        std::cout << "[acceptance] Period " << per
                  << " summary: total Ng = " << period_gen_sum
                  << " ; total Nr = " << period_rec_sum << "\n";
    }

    std::cout << "[acceptance] Filled acceptance triple columns; cells written: "
              << cells_written << "\n";
    return true;
}

// ------------- plotting -------------

struct CellData {
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> EX;
    std::vector<double> EY;
    std::vector<double> q2means;
    std::vector<double> tmeans;
};

static void draw_acceptance_canvases(const std::string& period_label,
                                     const CsvDoc& csv,
                                     const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const int c_xb_min  = csv.col_index("xBmin");
    const int c_xb_max  = csv.col_index("xBmax");
    const int c_q2_min  = csv.col_index("Q2min");
    const int c_q2_max  = csv.col_index("Q2max");
    const int c_tab_min = csv.col_index("t_abs_min");
    const int c_tab_max = csv.col_index("t_abs_max");
    const int c_phi_min = csv.col_index("phimin");
    const int c_phi_max = csv.col_index("phimax");

    if (c_xb_min < 0 || c_xb_max < 0 || c_q2_min < 0 || c_q2_max < 0 ||
        c_tab_min < 0 || c_tab_max < 0 || c_phi_min < 0 || c_phi_max < 0) {
        std::cerr << "[acceptance] FATAL: missing bin-edge columns needed for plotting.\n";
        std::exit(EXIT_FAILURE);
    }

    int c_phiavg = csv.col_index("phiavg, " + period_label);
    int c_q2avg  = csv.col_index("Q2avg, " + period_label);
    int c_tabavg = csv.col_index("t_abs_avg, " + period_label);
    int c_xbavg  = csv.col_index("xBavg, " + period_label);

    const std::string acc_col_name = "acceptance, " + period_label;
    const int c_acc = csv.col_index(acc_col_name);
    if (c_acc < 0) {
        std::cerr << "[acceptance] FATAL: missing column '" << acc_col_name
                  << "' for plotting.\n";
        std::exit(EXIT_FAILURE);
    }

    std::set<std::pair<double,double>> xb_set;
    for (int r = 0; r < csv.nrows(); ++r) {
        xb_set.emplace(csv.as_double(r, c_xb_min), csv.as_double(r, c_xb_max));
    }

    const double head_size = 0.14;
    const double label_sz  = 0.050;
    const double title_sz  = 0.045;

    const std::string period_dir = canonical_period_dir(period_label);
    const std::string base_dir =
        (fs::path(out_root_dir) / "acceptance" / period_dir).string();
    ensure_dir(base_dir);

    for (auto xb : xb_set) {
        std::set<std::pair<double,double>> q2set;
        std::set<std::pair<double,double>> tset_all;

        for (int r = 0; r < csv.nrows(); ++r) {
            const double xbmin = csv.as_double(r, c_xb_min);
            const double xbmax = csv.as_double(r, c_xb_max);
            if (std::fabs(xbmin - xb.first) < 1e-9 &&
                std::fabs(xbmax - xb.second) < 1e-9) {
                q2set.emplace(csv.as_double(r, c_q2_min),
                              csv.as_double(r, c_q2_max));
            }
        }
        for (auto q2r : q2set) {
            for (int r = 0; r < csv.nrows(); ++r) {
                const double xbmin = csv.as_double(r, c_xb_min);
                const double xbmax = csv.as_double(r, c_xb_max);
                const double q2min = csv.as_double(r, c_q2_min);
                const double q2max = csv.as_double(r, c_q2_max);
                if (std::fabs(xbmin - xb.first) < 1e-9 &&
                    std::fabs(xbmax - xb.second) < 1e-9 &&
                    std::fabs(q2min - q2r.first) < 1e-9 &&
                    std::fabs(q2max - q2r.second) < 1e-9) {
                    tset_all.emplace(csv.as_double(r, c_tab_min),
                                     csv.as_double(r, c_tab_max));
                }
            }
        }

        std::vector<std::pair<double,double>> Q2s(q2set.begin(), q2set.end());
        std::vector<std::pair<double,double>> Ts (tset_all.begin(), tset_all.end());
        if (Q2s.empty() || Ts.empty()) continue;

        const int nrows = (int)Q2s.size();
        const int ncols = (int)Ts.size();
        const int W = 300 * ncols + 160;
        const int H = 260 * nrows + 240;

        std::vector<CellData> cells(nrows * ncols);
        double canvas_max = 0.0;

        std::vector<double> xbmeans;
        for (int r = 0; r < csv.nrows(); ++r) {
            const double xbmin = csv.as_double(r, c_xb_min);
            const double xbmax = csv.as_double(r, c_xb_max);
            if (std::fabs(xbmin - xb.first) < 1e-9 &&
                std::fabs(xbmax - xb.second) < 1e-9) {
                if (c_xbavg >= 0) {
                    xbmeans.push_back(csv.as_double(r, c_xbavg));
                } else {
                    xbmeans.push_back(0.5 * (xb.first + xb.second));
                }
            }
        }
        const double xb_mean_for_title = safe_mean(xbmeans);

        for (int rr = 0; rr < nrows; ++rr) {
            for (int cc = 0; cc < ncols; ++cc) {
                const auto& qpair = Q2s[rr];
                const auto& tpair = Ts[cc];

                std::vector<int> rows_for_cell;
                for (int r = 0; r < csv.nrows(); ++r) {
                    const double xbmin = csv.as_double(r, c_xb_min);
                    const double xbmax = csv.as_double(r, c_xb_max);
                    const double q2min = csv.as_double(r, c_q2_min);
                    const double q2max = csv.as_double(r, c_q2_max);
                    const double tmin  = csv.as_double(r, c_tab_min);
                    const double tmax  = csv.as_double(r, c_tab_max);
                    if (std::fabs(xbmin - xb.first) < 1e-9 &&
                        std::fabs(xbmax - xb.second) < 1e-9 &&
                        std::fabs(q2min - qpair.first) < 1e-9 &&
                        std::fabs(q2max - qpair.second) < 1e-9 &&
                        std::fabs(tmin  - tpair.first) < 1e-9 &&
                        std::fabs(tmax  - tpair.second) < 1e-9) {
                        rows_for_cell.push_back(r);
                    }
                }

                std::sort(rows_for_cell.begin(), rows_for_cell.end(),
                          [&](int a, int b) {
                              return csv.as_double(a, c_phi_min) <
                                     csv.as_double(b, c_phi_min);
                          });

                CellData& C = cells[rr * ncols + cc];
                C.X.reserve(rows_for_cell.size());
                C.EX.assign(rows_for_cell.size(), 0.0);

                for (int r : rows_for_cell) {
                    const double pmin = csv.as_double(r, c_phi_min);
                    const double pmax = csv.as_double(r, c_phi_max);
                    double xphi = 0.5 * (pmin + pmax);
                    if (c_phiavg >= 0) {
                        const double pav = csv.as_double(r, c_phiavg);
                        if (std::isfinite(pav) && pav > 0.0 && pav < 360.0) {
                            xphi = pav;
                        }
                    }

                    const std::string& acc_cell = csv.rows[r][c_acc];
                    if (acc_cell.empty()) {
                        continue;
                    }
                    double acc_val = 0.0;
                    double acc_stat = 0.0;
                    double acc_sys = 0.0;
                    if (!parse_triple(acc_cell, acc_val, acc_stat, acc_sys)) {
                        double tmp = CsvDoc::to_double(acc_cell);
                        if (!std::isfinite(tmp)) continue;
                        acc_val  = tmp;
                        acc_stat = 0.0;
                    }

                    C.X.push_back(xphi);
                    C.Y.push_back(acc_val);
                    C.EY.push_back(acc_stat);
                    C.EX.push_back(0.0);

                    const double q2m = (c_q2avg >= 0)
                        ? csv.as_double(r, c_q2avg)
                        : 0.5 * (qpair.first + qpair.second);
                    const double tm  = (c_tabavg >= 0)
                        ? csv.as_double(r, c_tabavg)
                        : 0.5 * (tpair.first + tpair.second);
                    C.q2means.push_back(q2m);
                    C.tmeans.push_back(tm);

                    if (std::isfinite(acc_val)) {
                        canvas_max = std::max(canvas_max, acc_val);
                    }
                }
            }
        }

        std::string cname = "c_acc_" + period_dir + "_xB_" +
                            std::to_string((int)std::round(xb.first * 1000.0));
        TCanvas* c = new TCanvas(cname.c_str(), "", W, H);

        TPad* pTop  = new TPad("pTop", "pTop", 0.0, 0.86, 1.0, 1.0);
        TPad* pGrid = new TPad("pGrid","pGrid",0.0, 0.00, 1.0, 0.86);
        pTop->SetFillStyle(0);  pTop->Draw();
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        TLatex head;
        head.SetNDC();
        head.SetTextSize(head_size);
        head.SetTextAlign(22);
        head.DrawLatex(0.5, 0.58,
                       Form("%s   <xB>=%.2f",
                            period_label.c_str(),
                            xb_mean_for_title));

        const double y_lo = 0.0;
        const double y_hi = std::max(1.0, canvas_max * 1.15);

        for (int rr = 0; rr < nrows; ++rr) {
            for (int cc = 0; cc < ncols; ++cc) {
                pGrid->cd(rr * ncols + cc + 1);
                gPad->SetGrid(1, 1);
                gPad->SetTopMargin(0.24);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.07);
                gPad->SetTickx(1);
                gPad->SetTicky(1);

                TH1* frame = gPad->DrawFrame(0.0, y_lo, 360.0, y_hi);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Acceptance");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(title_sz);
                frame->GetYaxis()->SetTitleSize(title_sz);
                frame->GetXaxis()->SetLabelSize(label_sz);
                frame->GetYaxis()->SetLabelSize(label_sz);
                frame->GetYaxis()->SetTitleOffset(1.25);
                frame->GetXaxis()->SetTitleOffset(1.05);

                const CellData& C = cells[rr * ncols + cc];

                if (C.X.empty()) {
                    continue;
                }

                TGraphErrors* gacc = new TGraphErrors(
                    (int)C.X.size(),
                    (double*)C.X.data(),
                    (double*)C.Y.data(),
                    (double*)C.EX.data(),
                    (double*)C.EY.data());

                gacc->SetMarkerStyle(20);
                gacc->SetMarkerColor(kBlack);
                gacc->SetLineColor(kBlack);
                gacc->SetLineWidth(1);

                gacc->Draw("PE1 SAME");

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.058);
                lab.SetTextAlign(13);
                const double q2m = safe_mean(C.q2means);
                const double tm  = safe_mean(C.tmeans);
                lab.DrawLatex(0.16, 0.88,
                              Form("Q^{2}=%.2f  |t|=%.2f", q2m, tm));
            }
        }

        const std::string fpath =
            (fs::path(base_dir) /
             ("plot_acceptance_" +
              period_dir + "_xB_" +
              std::to_string((int)std::round(xb.first * 1000.0)) +
              ".png")).string();

        c->SaveAs(fpath.c_str());
        delete c;
    }
}

} // end anonymous namespace

// ------------- public entry point -------------

bool update_acceptance_csv(const std::string& csv_path,
                           const std::map<std::string, TTree*>& genMcTrees,
                           const std::map<std::string, TTree*>& recMcTrees,
                           const std::string& combined_cuts_json,
                           const std::string& global_cuts_json, // unused, kept for API symmetry
                           const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    (void)global_cuts_json; // global cuts handled via global_cuts.h

    const std::string csv_abs = fs::absolute(csv_path).string();
    std::error_code ec;
    const uintmax_t size_before =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[acceptance] CSV: " << csv_abs
              << " (size=" << size_before << ")\n";

    CsvDoc csv;
    if (!csv.load(csv_path)) {
        std::cerr << "[acceptance] ERROR: failed to load CSV.\n";
        return false;
    }

    std::vector<RowBin> bins = build_row_bins(csv);
    std::map<std::string, std::vector<bool>> has_data = build_row_has_data(csv);

    const TopoCutMap sigmaCuts   = load_sigma_cuts_data(combined_cuts_json);

    const auto tagMap = build_mc_tag_map();
    std::map<std::string, std::vector<double>> gen_all;
    std::map<std::string, std::vector<double>> rec_all;

    for (const auto& per : kPeriods) {
        auto itTag = tagMap.find(per);
        if (itTag == tagMap.end()) {
            std::cerr << "[acceptance] FATAL: no MC tag mapping for period "
                      << per << ".\n";
            return false;
        }
        const McTags& tags = itTag->second;

        auto itG = genMcTrees.find(tags.genTag);
        auto itR = recMcTrees.find(tags.recTag);
        if (itG == genMcTrees.end() || itR == recMcTrees.end()) {
            std::cerr << "[acceptance] FATAL: missing MC tree(s) for period "
                      << per << " (expected tags: gen="
                      << tags.genTag << " rec=" << tags.recTag << ").\n";
            return false;
        }

        const std::vector<bool>& flags = has_data[per];

        std::vector<double> gen_counts;
        std::vector<double> rec_counts;
        accumulate_counts_for_period(per,
                                     csv,
                                     bins,
                                     flags,
                                     itG->second,
                                     itR->second,
                                     sigmaCuts,
                                     gen_counts,
                                     rec_counts);
        gen_all[per] = gen_counts;
        rec_all[per] = rec_counts;
    }

    if (!fill_acceptance_columns(csv, gen_all, rec_all, has_data)) {
        std::cerr << "[acceptance] ERROR: fill_acceptance_columns failed.\n";
        return false;
    }

    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[acceptance] ERROR: failed to save updated CSV.\n";
        return false;
    }

    const uintmax_t size_after =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[acceptance] Updated CSV: " << csv_abs
              << " (size " << size_before << " -> " << size_after << ")\n";

    for (const auto& per : kPeriods) {
        draw_acceptance_canvases(per, csv, out_root_dir);
    }

    std::cout << "[acceptance] Acceptance plotting finished.\n";
    return true;
}