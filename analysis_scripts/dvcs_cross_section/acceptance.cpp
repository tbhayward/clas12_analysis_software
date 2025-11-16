// acceptance.cpp
// DVCS acceptance from non-radiative MC:
//   - acceptance, <period> = N_rec / N_gen per (xB, Q2, |t|, phi) bin
//   - N_rec from reconstructed MC with global DVCS exclusivity cuts
//     AND 3 sigma cuts loaded from combined_cuts.json
//   - N_gen from generated MC without exclusivity cuts
//   - Periods: Fa18 Inb, Fa18 Out, Sp19 Inb, Sp18 Inb, Sp18 Out
//   - Produces per-period acceptance vs phi plots in the usual
//     xB by (Q2, |t|) canvas layout under output/acceptance/<PeriodDir>/.

#include "acceptance.h"

#include <nlohmann/json.hpp>

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TLegend.h>
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

// periods
static const std::vector<std::string> kPeriods = {
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out"
};

// Map period label -> MC tree tags
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

// 3 sigma config for DVCS MC (pTmiss)
struct ThreeSigmaConfig {
    bool have_pTmiss;
    double mean_pTmiss;
    double sigma_pTmiss;
};

// Recursively search JSON for a single pTmiss node with "mean" and "sigma".
static void find_pTmiss_node(const json& j,
                             bool& found,
                             double& mean,
                             double& sigma,
                             int& count)
{
    if (found && count > 1) return;

    if (j.is_object()) {
        auto it = j.find("pTmiss");
        if (it != j.end() && it->is_object()) {
            const json& node = *it;
            if (node.contains("mean") && node.contains("sigma") &&
                node["mean"].is_number() && node["sigma"].is_number()) {
                ++count;
                mean  = node["mean"].get<double>();
                sigma = node["sigma"].get<double>();
                found = true;
            }
        }
        for (auto it2 = j.begin(); it2 != j.end(); ++it2) {
            find_pTmiss_node(it2.value(), found, mean, sigma, count);
        }
    } else if (j.is_array()) {
        for (const auto& el : j) {
            find_pTmiss_node(el, found, mean, sigma, count);
        }
    }
}

static ThreeSigmaConfig load_three_sigma_config(const std::string& cuts_json_path) {
    ThreeSigmaConfig cfg;
    cfg.have_pTmiss   = false;
    cfg.mean_pTmiss   = 0.0;
    cfg.sigma_pTmiss  = 0.0;

    std::ifstream fin(cuts_json_path);
    if (!fin.is_open()) {
        std::cerr << "[acceptance] FATAL: cannot open cuts JSON: "
                  << cuts_json_path << "\n";
        std::exit(EXIT_FAILURE);
    }

    json j;
    try {
        fin >> j;
    } catch (const std::exception& e) {
        std::cerr << "[acceptance] FATAL: failed to parse cuts JSON ("
                  << e.what() << ")\n";
        std::exit(EXIT_FAILURE);
    }

    bool found = false;
    double mean = 0.0;
    double sigma = 0.0;
    int count = 0;
    find_pTmiss_node(j, found, mean, sigma, count);

    if (!found) {
        std::cerr << "[acceptance] FATAL: no pTmiss entry with mean and sigma "
                  << "found in cuts JSON " << cuts_json_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (count > 1) {
        std::cerr << "[acceptance] FATAL: multiple pTmiss entries with mean and sigma "
                  << "found in cuts JSON " << cuts_json_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (!(sigma > 0.0)) {
        std::cerr << "[acceptance] FATAL: non-positive sigma for pTmiss in cuts JSON\n";
        std::exit(EXIT_FAILURE);
    }

    cfg.have_pTmiss  = true;
    cfg.mean_pTmiss  = mean;
    cfg.sigma_pTmiss = sigma;

    std::cout << "[acceptance] Loaded 3 sigma pTmiss cuts: mean="
              << cfg.mean_pTmiss << " sigma=" << cfg.sigma_pTmiss
              << " (3 sigma window will be applied)\n";

    return cfg;
}

// ---------- binning helpers ----------

struct BinKey {
    int ix;
    int iq;
    int it;
    int ip;

    bool operator<(const BinKey& o) const {
        if (ix != o.ix) return ix < o.ix;
        if (iq != o.iq) return iq < o.iq;
        if (it != o.it) return it < o.it;
        return ip < o.ip;
    }
};

struct Binning {
    std::vector<double> xb_edges;
    std::vector<double> q2_edges;
    std::vector<double> t_edges;
    std::vector<double> phi_edges;
    std::map<BinKey,int> key_to_row;
};

static void build_axis_edges(const CsvDoc& csv,
                             int c_min,
                             int c_max,
                             std::vector<double>& edges)
{
    std::set<double> s;
    const int NR = csv.nrows();
    for (int r = 0; r < NR; ++r) {
        const double vmin = csv.as_double(r, c_min);
        const double vmax = csv.as_double(r, c_max);
        if (std::isfinite(vmin)) s.insert(vmin);
        if (std::isfinite(vmax)) s.insert(vmax);
    }
    edges.assign(s.begin(), s.end());
    if (edges.size() < 2) {
        std::cerr << "[acceptance] FATAL: less than two edges in axis.\n";
        std::exit(EXIT_FAILURE);
    }
}

static int find_interval_index(double vmin,
                               double vmax,
                               const std::vector<double>& edges)
{
    int idx_min = -1;
    int idx_max = -1;
    const int N = (int)edges.size();
    for (int i = 0; i < N; ++i) {
        if (std::fabs(edges[i] - vmin) < 1e-9) idx_min = i;
        if (std::fabs(edges[i] - vmax) < 1e-9) idx_max = i;
    }
    if (idx_min < 0 || idx_max < 0) return -1;
    if (idx_max != idx_min + 1) return -1;
    return idx_min;
}

static int find_value_index(double v,
                            const std::vector<double>& edges)
{
    if (!std::isfinite(v)) return -1;
    const int N = (int)edges.size();
    if (v < edges.front() || v > edges.back()) return -1;
    if (std::fabs(v - edges.back()) < 1e-9) return N - 2;
    auto it = std::upper_bound(edges.begin(), edges.end(), v);
    if (it == edges.begin() || it == edges.end()) return -1;
    int idx = (int)(it - edges.begin()) - 1;
    if (idx < 0 || idx >= N - 1) return -1;
    return idx;
}

static Binning build_binning_from_csv(const CsvDoc& csv) {
    Binning B;

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

    build_axis_edges(csv, c_xb_min,  c_xb_max,  B.xb_edges);
    build_axis_edges(csv, c_q2_min,  c_q2_max,  B.q2_edges);
    build_axis_edges(csv, c_tab_min, c_tab_max, B.t_edges);
    build_axis_edges(csv, c_phi_min, c_phi_max, B.phi_edges);

    const int NR = csv.nrows();
    for (int r = 0; r < NR; ++r) {
        const double xbmin  = csv.as_double(r, c_xb_min);
        const double xbmax  = csv.as_double(r, c_xb_max);
        const double q2min  = csv.as_double(r, c_q2_min);
        const double q2max  = csv.as_double(r, c_q2_max);
        const double tmin   = csv.as_double(r, c_tab_min);
        const double tmax   = csv.as_double(r, c_tab_max);
        const double phimin = csv.as_double(r, c_phi_min);
        const double phimax = csv.as_double(r, c_phi_max);

        const int ix = find_interval_index(xbmin, xbmax, B.xb_edges);
        const int iq = find_interval_index(q2min, q2max, B.q2_edges);
        const int it = find_interval_index(tmin,  tmax,  B.t_edges);
        const int ip = find_interval_index(phimin,phimax,B.phi_edges);

        if (ix < 0 || iq < 0 || it < 0 || ip < 0) {
            std::cerr << "[acceptance] FATAL: could not map CSV row " << r
                      << " to a unique 4D bin (xB,Q2,|t|,phi).\n";
            std::exit(EXIT_FAILURE);
        }

        BinKey key;
        key.ix = ix;
        key.iq = iq;
        key.it = it;
        key.ip = ip;

        auto itExisting = B.key_to_row.find(key);
        if (itExisting != B.key_to_row.end()) {
            std::cerr << "[acceptance] FATAL: duplicate bin key for rows "
                      << itExisting->second << " and " << r << ".\n";
            std::exit(EXIT_FAILURE);
        }
        B.key_to_row[key] = r;
    }

    return B;
}

// ---------- exclusivity for reconstructed MC (global + 3 sigma) ----------

static bool rec_passes_exclusivity(double t1,
                                   double open_angle_ep2_deg,
                                   double pTmiss,
                                   const ThreeSigmaConfig& cfg)
{
    const double t_abs = std::fabs(t1);
    if (!std::isfinite(t_abs) || t_abs <= 0.0 || t_abs >= 1.0) return false;
    if (!std::isfinite(open_angle_ep2_deg) || open_angle_ep2_deg <= 5.0) return false;
    if (!std::isfinite(pTmiss)) return false;
    if (pTmiss > 0.20) return false;

    if (cfg.have_pTmiss) {
        const double dev = pTmiss - cfg.mean_pTmiss;
        if (std::fabs(dev) > 3.0 * cfg.sigma_pTmiss) return false;
    }

    return true;
}

// ---------- MC counting ----------

static void fill_mc_counts_for_period(const std::string& period_label,
                                      const Binning& binning,
                                      const CsvDoc& csv,
                                      TTree* genTree,
                                      TTree* recTree,
                                      const ThreeSigmaConfig& cfg,
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
    const char* br_pTmiss = "pTmiss";

    if (!genTree->GetBranch(br_x) ||
        !genTree->GetBranch(br_Q2) ||
        !genTree->GetBranch(br_t1) ||
        !genTree->GetBranch(br_phi2)) {
        std::cerr << "[acceptance] FATAL: missing one or more branches in genTree for period "
                  << period_label << " (expected: x, Q2, t1, phi2).\n";
        std::exit(EXIT_FAILURE);
    }

    if (!recTree->GetBranch(br_x) ||
        !recTree->GetBranch(br_Q2) ||
        !recTree->GetBranch(br_t1) ||
        !recTree->GetBranch(br_phi2) ||
        !recTree->GetBranch(br_oa_ep2) ||
        !recTree->GetBranch(br_pTmiss)) {
        std::cerr << "[acceptance] FATAL: missing one or more branches in recTree for period "
                  << period_label
                  << " (expected: x, Q2, t1, phi2, open_angle_ep2, pTmiss).\n";
        std::exit(EXIT_FAILURE);
    }

    // generated MC
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

    for (Long64_t i = 0; i < Ngen; ++i) {
        genTree->GetEntry(i);

        const double xB   = g_x;
        const double Q2   = g_Q2;
        const double tAbs = std::fabs(g_t1);
        const double phiD = wrap_deg(RAD2DEG(g_phi2));

        const int ix = find_value_index(xB,   binning.xb_edges);
        const int iq = find_value_index(Q2,   binning.q2_edges);
        const int it = find_value_index(tAbs, binning.t_edges);
        const int ip = find_value_index(phiD, binning.phi_edges);

        if (ix < 0 || iq < 0 || it < 0 || ip < 0) {
            continue;
        }

        BinKey key;
        key.ix = ix;
        key.iq = iq;
        key.it = it;
        key.ip = ip;

        auto itRow = binning.key_to_row.find(key);
        if (itRow == binning.key_to_row.end()) {
            continue;
        }

        const int row = itRow->second;
        if (row < 0 || row >= NR) continue;

        gen_counts[row] += 1.0;
        ++used_gen;
    }

    std::cout << "[acceptance] Period " << period_label
              << " gen MC: total entries = " << Ngen
              << " ; binned = " << used_gen << "\n";

    // reconstructed MC
    double r_x      = 0.0;
    double r_Q2     = 0.0;
    double r_t1     = 0.0;
    double r_phi2   = 0.0;
    double r_oa     = 0.0;
    double r_pTmiss = 0.0;

    recTree->SetBranchAddress(br_x,      &r_x);
    recTree->SetBranchAddress(br_Q2,     &r_Q2);
    recTree->SetBranchAddress(br_t1,     &r_t1);
    recTree->SetBranchAddress(br_phi2,   &r_phi2);
    recTree->SetBranchAddress(br_oa_ep2, &r_oa);
    recTree->SetBranchAddress(br_pTmiss, &r_pTmiss);

    const Long64_t Nrec = recTree->GetEntries();
    Long64_t used_rec = 0;
    Long64_t passed_excl = 0;

    for (Long64_t i = 0; i < Nrec; ++i) {
        recTree->GetEntry(i);

        if (!rec_passes_exclusivity(r_t1, r_oa, r_pTmiss, cfg)) {
            continue;
        }
        ++passed_excl;

        const double xB   = r_x;
        const double Q2   = r_Q2;
        const double tAbs = std::fabs(r_t1);
        const double phiD = wrap_deg(RAD2DEG(r_phi2));

        const int ix = find_value_index(xB,   binning.xb_edges);
        const int iq = find_value_index(Q2,   binning.q2_edges);
        const int it = find_value_index(tAbs, binning.t_edges);
        const int ip = find_value_index(phiD, binning.phi_edges);

        if (ix < 0 || iq < 0 || it < 0 || ip < 0) {
            continue;
        }

        BinKey key;
        key.ix = ix;
        key.iq = iq;
        key.it = it;
        key.ip = ip;

        auto itRow = binning.key_to_row.find(key);
        if (itRow == binning.key_to_row.end()) {
            continue;
        }

        const int row = itRow->second;
        if (row < 0 || row >= NR) continue;

        rec_counts[row] += 1.0;
        ++used_rec;
    }

    std::cout << "[acceptance] Period " << period_label
              << " rec MC: total entries = " << Nrec
              << " ; passed exclusivity = " << passed_excl
              << " ; binned = " << used_rec << "\n";
}

// ---------- CSV update and plotting ----------

static bool fill_acceptance_columns(CsvDoc& csv,
                                    const std::map<std::string, std::vector<double>>& gen_all,
                                    const std::map<std::string, std::vector<double>>& rec_all)
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
        if (itG == gen_all.end() || itR == rec_all.end()) {
            std::cerr << "[acceptance] FATAL: internal error, missing counts for period "
                      << per << ".\n";
            return false;
        }
        const std::vector<double>& gen = itG->second;
        const std::vector<double>& rec = itR->second;
        if ((int)gen.size() != NR || (int)rec.size() != NR) {
            std::cerr << "[acceptance] FATAL: size mismatch for counts vectors in period "
                      << per << ".\n";
            return false;
        }

        const int c_acc = acc_idx[per];
        double period_gen_sum = 0.0;
        double period_rec_sum = 0.0;

        for (int r = 0; r < NR; ++r) {
            const double Ng = gen[r];
            const double Nr = rec[r];

            period_gen_sum += Ng;
            period_rec_sum += Nr;

            double acc = 0.0;
            if (Ng > 0.0) {
                acc = Nr / Ng;
            } else {
                acc = 0.0;
            }

            std::ostringstream oss;
            oss.setf(std::ios::fixed);
            oss << std::setprecision(8) << acc;
            csv.rows[r][c_acc] = oss.str();
            ++cells_written;
        }

        std::cout << "[acceptance] Period " << per
                  << " summary: total Ng = " << period_gen_sum
                  << " ; total Nr = " << period_rec_sum << "\n";
    }

    std::cout << "[acceptance] Filled acceptance columns; cells written: "
              << cells_written << "\n";
    return true;
}

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
                C.EY.assign(rows_for_cell.size(), 0.0);

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
                    C.X.push_back(xphi);

                    const double acc = csv.as_double(r, c_acc);
                    C.Y.push_back(acc);

                    const double q2m = (c_q2avg >= 0)
                        ? csv.as_double(r, c_q2avg)
                        : 0.5 * (qpair.first + qpair.second);
                    const double tm  = (c_tabavg >= 0)
                        ? csv.as_double(r, c_tabavg)
                        : 0.5 * (tpair.first + tpair.second);
                    C.q2means.push_back(q2m);
                    C.tmeans.push_back(tm);

                    if (std::isfinite(acc)) {
                        canvas_max = std::max(canvas_max, acc);
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

bool update_acceptance_csv(const std::string& csv_path,
                           const std::map<std::string, TTree*>& genMcTrees,
                           const std::map<std::string, TTree*>& recMcTrees,
                           const std::string& cuts_json,
                           const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

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

    // Build binning from CSV once
    Binning binning = build_binning_from_csv(csv);

    // Load 3 sigma config (pTmiss) once and use for all periods
    ThreeSigmaConfig cfg = load_three_sigma_config(cuts_json);

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

        std::vector<double> gen_counts;
        std::vector<double> rec_counts;
        fill_mc_counts_for_period(per, binning, csv,
                                  itG->second, itR->second,
                                  cfg,
                                  gen_counts, rec_counts);
        gen_all[per] = gen_counts;
        rec_all[per] = rec_counts;
    }

    if (!fill_acceptance_columns(csv, gen_all, rec_all)) {
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