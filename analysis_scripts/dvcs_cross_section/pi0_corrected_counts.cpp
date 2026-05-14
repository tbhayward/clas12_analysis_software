// pi0_corrected_counts.cpp
// Pi0-corrected DVCS signal yields:
//   - Fill "signal yield, ep->epg, exp, <period>, <hel>" columns
//     with "(value, stat, sys)" triples.
//   - Uses normalized raw DVCS yields, not raw unnormalized counts.
//   - Sp18 Inb and Sp18 Out are unpolarized-only because those periods do not
//     have valid helicity-resolved accumulated charge.
//   - Produce per-period yield plots vs phi for helicity-resolved periods.

#include "pi0_corrected_counts.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TLegend.h>

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
        const bool needq =
            s.find(',') != std::string::npos ||
            s.find('"') != std::string::npos;

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
        if (s.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        char* e = nullptr;
        const double v = std::strtod(s.c_str(), &e);

        if (e == s.c_str()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        return v;
    }

    bool load(const std::string& path) {
        std::ifstream fin(path);

        if (!fin.is_open()) {
            std::cerr << "[pi0_corrected] ERROR: cannot open CSV: "
                      << path << "\n";
            return false;
        }

        std::string line;

        if (!std::getline(fin, line)) {
            std::cerr << "[pi0_corrected] ERROR: empty CSV: "
                      << path << "\n";
            return false;
        }

        header = split_csv_line(line);

        index.clear();

        for (int i = 0; i < (int)header.size(); ++i) {
            index[header[i]] = i;
        }

        rows.clear();

        while (std::getline(fin, line)) {
            if (line.empty()) {
                continue;
            }

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
                std::cerr << "[pi0_corrected] ERROR: cannot write CSV tmp: "
                          << tmp << "\n";
                return false;
            }

            for (size_t i = 0; i < header.size(); ++i) {
                write_field(fout, header[i]);

                if (i + 1 < header.size()) {
                    fout << ',';
                }
            }

            fout << "\n";

            for (const auto& row : rows) {
                for (size_t i = 0; i < row.size(); ++i) {
                    write_field(fout, row[i]);

                    if (i + 1 < row.size()) {
                        fout << ',';
                    }
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
                std::cerr << "[pi0_corrected] ERROR: atomic rename failed ("
                          << ec.message() << ")\n";
                return false;
            }
        }

        return true;
    }

    int nrows() const {
        return (int)rows.size();
    }

    bool has_col(const std::string& name) const {
        return index.find(name) != index.end();
    }

    int col_index(const std::string& name) const {
        auto it = index.find(name);
        return (it == index.end()) ? -1 : it->second;
    }

    double as_double(int r, int c) const {
        if (r < 0 || r >= nrows()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (c < 0 || c >= (int)header.size()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        return to_double(rows[r][c]);
    }
};

static inline double PI() {
    return 3.14159265358979323846;
}

static inline double RAD2DEG(double r) {
    return r * 180.0 / PI();
}

static inline std::string to_lower_nospace(std::string s) {
    std::string out;
    out.reserve(s.size());

    for (char c : s) {
        if (c == ' ' || c == '\t' || c == '_') {
            continue;
        }

        out.push_back((char)std::tolower((unsigned char)c));
    }

    return out;
}

static std::string canonical_period_dir(const std::string& L) {
    const std::string k = to_lower_nospace(L);

    if (k == "fa18inb") {
        return "Fa18_Inb";
    }

    if (k == "fa18out") {
        return "Fa18_Out";
    }

    if (k == "fa18inbsupp") {
        return "Fa18_Inb_Supp";
    }

    if (k == "sp18inb") {
        return "Sp18_Inb";
    }

    if (k == "sp18out") {
        return "Sp18_Out";
    }

    if (k == "sp19inb") {
        return "Sp19_Inb";
    }

    std::string s = L;

    for (char& c : s) {
        if (c == ' ') {
            c = '_';
        }
    }

    return s;
}

static void ensure_dir(const std::string& p) {
    namespace fs = std::filesystem;

    std::error_code ec;

    if (!fs::exists(p)) {
        fs::create_directories(p, ec);

        if (ec) {
            std::cerr << "[pi0_corrected] FATAL: could not create directory: "
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

static bool parse_tuple3(const std::string& s,
                         double& value,
                         double& stat,
                         double& sys) {
    value = std::numeric_limits<double>::quiet_NaN();
    stat  = std::numeric_limits<double>::quiet_NaN();
    sys   = std::numeric_limits<double>::quiet_NaN();

    std::string trimmed;
    trimmed.reserve(s.size());

    for (char c : s) {
        if (!std::isspace((unsigned char)c)) {
            trimmed.push_back(c);
        }
    }

    if (trimmed.empty()) {
        return false;
    }

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
        std::cerr << "[pi0_corrected] ERROR: tuple has "
                  << parts.size() << " fields (expected 3): '"
                  << s << "'\n";
        return false;
    }

    char* e1 = nullptr;
    char* e2 = nullptr;
    char* e3 = nullptr;

    value = std::strtod(parts[0].c_str(), &e1);
    stat  = std::strtod(parts[1].c_str(), &e2);
    sys   = std::strtod(parts[2].c_str(), &e3);

    if (e1 == parts[0].c_str()) {
        return false;
    }

    if (e2 == parts[1].c_str()) {
        return false;
    }

    if (e3 == parts[2].c_str()) {
        return false;
    }

    return true;
}

static bool parse_tuple3_or_empty(const std::string& s,
                                  double& value,
                                  double& stat,
                                  double& sys) {
    value = 0.0;
    stat = 0.0;
    sys = 0.0;

    if (s.empty()) {
        return true;
    }

    return parse_tuple3(s, value, stat, sys);
}

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

static void add_tuple_in_quadrature(double v,
                                    double stat,
                                    double sys,
                                    double& sum_v,
                                    double& sum_stat_var,
                                    double& sum_sys_var) {
    sum_v += v;
    sum_stat_var += stat * stat;
    sum_sys_var += sys * sys;
}

// S = (1 - c) * N_norm
// Var_stat(S) = (1 - c)^2 Var_stat(N_norm) + N_norm^2 Var_stat(c)
static void compute_signal_and_stat(double norm_yield,
                                    double norm_yield_stat,
                                    double c_val,
                                    double c_stat,
                                    double& S,
                                    double& S_stat) {
    if (!std::isfinite(norm_yield) || norm_yield <= 0.0) {
        S = 0.0;
        S_stat = 0.0;
        return;
    }

    if (!std::isfinite(norm_yield_stat) || norm_yield_stat < 0.0) {
        norm_yield_stat = 0.0;
    }

    if (!std::isfinite(c_val)) {
        c_val = 0.0;
    }

    if (!std::isfinite(c_stat) || c_stat < 0.0) {
        c_stat = 0.0;
    }

    const double one_minus_c = 1.0 - c_val;

    const double var =
        one_minus_c * one_minus_c * norm_yield_stat * norm_yield_stat +
        norm_yield * norm_yield * c_stat * c_stat;

    S = one_minus_c * norm_yield;
    S_stat = std::sqrt(std::max(0.0, var));
}

static const std::vector<std::string> kPeriods = {
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out"
};

static const std::vector<std::string> kHelicities = {
    "unpol",
    "pos",
    "neg"
};

static const std::vector<std::string> kUnpolarizedOnlyHelicities = {
    "unpol"
};

static bool has_helicity_resolved_data_columns(const std::string& period) {
    if (period == "Sp18 Inb" || period == "Sp18 Out") {
        return false;
    }

    return true;
}

static const std::vector<std::string>& helicities_for_period(const std::string& period) {
    if (!has_helicity_resolved_data_columns(period)) {
        return kUnpolarizedOnlyHelicities;
    }

    return kHelicities;
}

static const std::vector<std::string> kTopos = {
    "(FD, FD)",
    "(CD, FD)",
    "(CD, FT)"
};

struct CellData {
    std::vector<double> X;
    std::vector<double> Yp;
    std::vector<double> Ym;
    std::vector<double> EX;
    std::vector<double> EYp;
    std::vector<double> EYm;
    std::vector<double> q2means;
    std::vector<double> tmeans;
};

static bool fill_signal_yields(CsvDoc& csv) {
    std::map<std::string,int> cont_idx;

    for (const auto& per : kPeriods) {
        const std::string cname = "contamination ratio, " + per;
        const int idx = csv.col_index(cname);

        if (idx < 0) {
            std::cerr << "[pi0_corrected] FATAL: missing column '"
                      << cname << "'. Run pi0_contamination first.\n";
            return false;
        }

        cont_idx[per] = idx;
    }

    std::map<std::string, std::map<std::string,int>> sig_idx;

    for (const auto& per : kPeriods) {
        for (const auto& hel : helicities_for_period(per)) {
            std::ostringstream name;
            name << "signal yield, ep->epg, exp, " << per << ", " << hel;

            const int idx = csv.col_index(name.str());

            if (idx < 0) {
                std::cerr << "[pi0_corrected] FATAL: missing signal-yield column '"
                          << name.str() << "' in CSV header.\n";
                return false;
            }

            sig_idx[per][hel] = idx;
        }
    }

    std::map<std::string, std::map<std::string, std::map<std::string,int>>> raw_idx;

    for (const auto& per : kPeriods) {
        for (const auto& topo : kTopos) {
            for (const auto& hel : helicities_for_period(per)) {
                std::ostringstream nm;

                nm << "normalized raw yield, ep->epg, " << topo
                   << ", exp, " << per << ", " << hel;

                const int idx = csv.col_index(nm.str());

                if (idx < 0) {
                    std::cerr << "[pi0_corrected] FATAL: missing normalized-raw-yield column '"
                              << nm.str() << "' in CSV header.\n";
                    return false;
                }

                raw_idx[per][topo][hel] = idx;
            }
        }
    }

    const int NR = csv.nrows();
    std::size_t cells_written = 0;

    for (int r = 0; r < NR; ++r) {
        for (const auto& per : kPeriods) {
            const int c_cont = cont_idx[per];
            const std::string& cs = csv.rows[r][c_cont];

            double c_val  = 0.0;
            double c_stat = 0.0;
            double c_sys  = 0.0;

            if (!cs.empty()) {
                double cv = 0.0;
                double cs_stat = 0.0;
                double cs_sys  = 0.0;

                const bool ok = parse_tuple3(cs, cv, cs_stat, cs_sys);

                if (!ok) {
                    std::cerr << "[pi0_corrected] FATAL: failed to parse contamination '"
                              << cs << "' for period " << per
                              << " row " << r << "\n";
                    return false;
                }

                c_val  = cv;
                c_stat = cs_stat;
                c_sys  = cs_sys;
            } else {
                c_val  = 0.0;
                c_stat = 0.0;
                c_sys  = 0.0;
            }

            std::map<std::string,double> raw_sum;
            std::map<std::string,double> raw_stat_var;
            std::map<std::string,double> raw_sys_var;

            for (const auto& hel : helicities_for_period(per)) {
                raw_sum[hel] = 0.0;
                raw_stat_var[hel] = 0.0;
                raw_sys_var[hel] = 0.0;
            }

            for (const auto& topo : kTopos) {
                for (const auto& hel : helicities_for_period(per)) {
                    const int c_raw = raw_idx[per][topo][hel];
                    const std::string& sraw = csv.rows[r][c_raw];

                    double v = 0.0;
                    double v_stat = 0.0;
                    double v_sys = 0.0;

                    if (!parse_tuple3_or_empty(sraw, v, v_stat, v_sys)) {
                        std::cerr << "[pi0_corrected] FATAL: failed to parse normalized-raw-yield triple in '"
                                  << "normalized raw yield, ep->epg, " << topo
                                  << ", exp, " << per << ", " << hel
                                  << "' at row " << r << " value '"
                                  << sraw << "'\n";
                        return false;
                    }

                    raw_sum[hel] += v;
                    raw_stat_var[hel] += v_stat * v_stat;
                    raw_sys_var[hel] += v_sys * v_sys;
                }
            }

            for (const auto& hel : helicities_for_period(per)) {
                const double norm = raw_sum[hel];
                const double norm_stat =
                    std::sqrt(std::max(0.0, raw_stat_var[hel]));

                double S = 0.0;
                double S_stat = 0.0;

                compute_signal_and_stat(norm,
                                        norm_stat,
                                        c_val,
                                        c_stat,
                                        S,
                                        S_stat);

                const int c_sig = sig_idx[per][hel];
                csv.rows[r][c_sig] = format_triple(S, S_stat, 0.0);

                ++cells_written;
            }
        }
    }

    std::cout << "[pi0_corrected] Filled signal-yield triple columns; cells written: "
              << cells_written << "\n";

    return true;
}

static void draw_signal_yield_canvases(const std::string& period_label,
                                       const CsvDoc& csv,
                                       const std::string& out_root_dir) {
    namespace fs = std::filesystem;

    if (!has_helicity_resolved_data_columns(period_label)) {
        std::cout << "[pi0_corrected] Skipping helicity-separated signal-yield plots for "
                  << period_label
                  << " because this period is unpolarized-only.\n";
        return;
    }

    const int c_xb_min  = csv.col_index("xBmin");
    const int c_xb_max  = csv.col_index("xBmax");
    const int c_q2_min  = csv.col_index("Q2min");
    const int c_q2_max  = csv.col_index("Q2max");
    const int c_tab_min = csv.col_index("t_abs_min");
    const int c_tab_max = csv.col_index("t_abs_max");
    const int c_phi_min = csv.col_index("phimin");
    const int c_phi_max = csv.col_index("phimax");

    if (c_xb_min < 0 ||
        c_xb_max < 0 ||
        c_q2_min < 0 ||
        c_q2_max < 0 ||
        c_tab_min < 0 ||
        c_tab_max < 0 ||
        c_phi_min < 0 ||
        c_phi_max < 0) {
        std::cerr << "[pi0_corrected] FATAL: missing bin-edge columns needed for plotting.\n";
        std::exit(EXIT_FAILURE);
    }

    const int c_phiavg = csv.col_index("phiavg, " + period_label);
    const int c_q2avg  = csv.col_index("Q2avg, " + period_label);
    const int c_tabavg = csv.col_index("t_abs_avg, " + period_label);
    const int c_xbavg  = csv.col_index("xBavg, " + period_label);

    const std::string cont_col_name = "contamination ratio, " + period_label;
    const int c_contam = csv.col_index(cont_col_name);

    if (c_contam < 0) {
        std::cerr << "[pi0_corrected] FATAL: missing column '"
                  << cont_col_name
                  << "'. Did you run pi0_contamination?\n";
        std::exit(EXIT_FAILURE);
    }

    std::map<std::string, std::map<std::string,int>> raw_idx;

    for (const auto& topo : kTopos) {
        for (const auto& hel : kHelicities) {
            std::ostringstream name;

            name << "normalized raw yield, ep->epg, " << topo
                 << ", exp, " << period_label
                 << ", " << hel;

            const int idx = csv.col_index(name.str());

            if (idx < 0) {
                std::cerr << "[pi0_corrected] FATAL: missing normalized-raw-yield column: '"
                          << name.str() << "'\n";
                std::exit(EXIT_FAILURE);
            }

            raw_idx[topo][hel] = idx;
        }
    }

    std::set<std::pair<double,double>> xb_set;

    for (int r = 0; r < csv.nrows(); ++r) {
        xb_set.emplace(csv.as_double(r, c_xb_min),
                       csv.as_double(r, c_xb_max));
    }

    const double head_size = 0.14;
    const double label_sz  = 0.050;
    const double title_sz  = 0.045;
    const double leg_txt   = 0.050;

    const std::string period_dir = canonical_period_dir(period_label);
    const std::string base_dir =
        (fs::path(out_root_dir) / "signal_yield_plots" / period_dir).string();

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
        std::vector<std::pair<double,double>> Ts(tset_all.begin(), tset_all.end());

        if (Q2s.empty() || Ts.empty()) {
            continue;
        }

        const int nrows = (int)Q2s.size();
        const int ncols = (int)Ts.size();
        const int W = 300 * ncols + 160;
        const int H = 260 * nrows + 240;

        std::vector<CellData> cells(nrows * ncols);
        double canvas_max = 1.0;

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
                        std::fabs(tmin - tpair.first) < 1e-9 &&
                        std::fabs(tmax - tpair.second) < 1e-9) {
                        rows_for_cell.push_back(r);
                    }
                }

                std::sort(rows_for_cell.begin(),
                          rows_for_cell.end(),
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

                    C.X.push_back(xphi);

                    const std::string& cs = csv.rows[r][c_contam];

                    double c_val  = 0.0;
                    double c_stat = 0.0;
                    double c_sys  = 0.0;

                    if (!cs.empty()) {
                        double cv = 0.0;
                        double cs_stat = 0.0;
                        double cs_sys  = 0.0;

                        const bool have_cont = parse_tuple3(cs,
                                                            cv,
                                                            cs_stat,
                                                            cs_sys);

                        if (!have_cont) {
                            std::cerr << "[pi0_corrected] FATAL: failed to parse contamination '"
                                      << cs << "' for period " << period_label
                                      << " row " << r << "\n";
                            std::exit(EXIT_FAILURE);
                        }

                        c_val  = cv;
                        c_stat = cs_stat;
                        c_sys  = cs_sys;
                    } else {
                        c_val  = 0.0;
                        c_stat = 0.0;
                        c_sys  = 0.0;
                    }

                    double norm_pos = 0.0;
                    double norm_neg = 0.0;
                    double norm_pos_stat_var = 0.0;
                    double norm_neg_stat_var = 0.0;
                    double norm_pos_sys_var = 0.0;
                    double norm_neg_sys_var = 0.0;

                    for (const auto& topo : kTopos) {
                        const int c_pos = raw_idx[topo].at("pos");
                        const int c_neg = raw_idx[topo].at("neg");

                        const std::string& s_pos = csv.rows[r][c_pos];
                        const std::string& s_neg = csv.rows[r][c_neg];

                        double vpos = 0.0;
                        double spos = 0.0;
                        double ypos = 0.0;

                        if (!parse_tuple3_or_empty(s_pos, vpos, spos, ypos)) {
                            std::cerr << "[pi0_corrected] FATAL: failed to parse normalized pos yield triple in '"
                                      << "normalized raw yield, ep->epg, " << topo
                                      << ", exp, " << period_label << ", pos"
                                      << "' at row " << r << " value '"
                                      << s_pos << "'\n";
                            std::exit(EXIT_FAILURE);
                        }

                        add_tuple_in_quadrature(vpos,
                                                spos,
                                                ypos,
                                                norm_pos,
                                                norm_pos_stat_var,
                                                norm_pos_sys_var);

                        double vneg = 0.0;
                        double sneg = 0.0;
                        double yneg = 0.0;

                        if (!parse_tuple3_or_empty(s_neg, vneg, sneg, yneg)) {
                            std::cerr << "[pi0_corrected] FATAL: failed to parse normalized neg yield triple in '"
                                      << "normalized raw yield, ep->epg, " << topo
                                      << ", exp, " << period_label << ", neg"
                                      << "' at row " << r << " value '"
                                      << s_neg << "'\n";
                            std::exit(EXIT_FAILURE);
                        }

                        add_tuple_in_quadrature(vneg,
                                                sneg,
                                                yneg,
                                                norm_neg,
                                                norm_neg_stat_var,
                                                norm_neg_sys_var);
                    }

                    double S_pos = 0.0;
                    double S_pos_stat = 0.0;
                    double S_neg = 0.0;
                    double S_neg_stat = 0.0;

                    compute_signal_and_stat(norm_pos,
                                            std::sqrt(std::max(0.0, norm_pos_stat_var)),
                                            c_val,
                                            c_stat,
                                            S_pos,
                                            S_pos_stat);

                    compute_signal_and_stat(norm_neg,
                                            std::sqrt(std::max(0.0, norm_neg_stat_var)),
                                            c_val,
                                            c_stat,
                                            S_neg,
                                            S_neg_stat);

                    C.Yp.push_back(S_pos);
                    C.Ym.push_back(S_neg);
                    C.EYp.push_back(S_pos_stat);
                    C.EYm.push_back(S_neg_stat);

                    const double q2m =
                        (c_q2avg >= 0)
                        ? csv.as_double(r, c_q2avg)
                        : 0.5 * (qpair.first + qpair.second);

                    const double tm =
                        (c_tabavg >= 0)
                        ? csv.as_double(r, c_tabavg)
                        : 0.5 * (tpair.first + tpair.second);

                    C.q2means.push_back(q2m);
                    C.tmeans.push_back(tm);

                    if (std::isfinite(S_pos)) {
                        canvas_max = std::max(canvas_max, S_pos);
                    }

                    if (std::isfinite(S_neg)) {
                        canvas_max = std::max(canvas_max, S_neg);
                    }
                }
            }
        }

        const std::string cname =
            "c_sig_" + period_dir + "_xB_" +
            std::to_string((int)std::round(xb.first * 1000.0));

        TCanvas* c = new TCanvas(cname.c_str(), "", W, H);

        TPad* pTop  = new TPad("pTop", "pTop", 0.0, 0.86, 1.0, 1.0);
        TPad* pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.86);

        pTop->SetFillStyle(0);
        pTop->Draw();

        pGrid->SetFillStyle(0);
        pGrid->Draw();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();

        TLatex head;
        head.SetNDC();
        head.SetTextSize(head_size);
        head.SetTextAlign(22);

        head.DrawLatex(0.5,
                       0.58,
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
                frame->GetYaxis()->SetTitle("Signal yield");
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

                TGraphErrors* gp = new TGraphErrors((int)C.X.size(),
                                                    (double*)C.X.data(),
                                                    (double*)C.Yp.data(),
                                                    (double*)C.EX.data(),
                                                    (double*)C.EYp.data());

                TGraphErrors* gm = new TGraphErrors((int)C.X.size(),
                                                    (double*)C.X.data(),
                                                    (double*)C.Ym.data(),
                                                    (double*)C.EX.data(),
                                                    (double*)C.EYm.data());

                gp->SetMarkerStyle(24);
                gp->SetMarkerColor(kRed);
                gp->SetLineColor(kRed);
                gp->SetLineWidth(1);

                gm->SetMarkerStyle(20);
                gm->SetMarkerColor(kBlue);
                gm->SetLineColor(kBlue);
                gm->SetLineWidth(1);

                gp->Draw("PE1 SAME");
                gm->Draw("PE1 SAME");

                TLegend* leg = new TLegend(0.60, 0.72, 0.93, 0.92);

                leg->SetBorderSize(1);
                leg->SetLineColor(kBlack);
                leg->SetFillStyle(1001);
                leg->SetFillColor(kWhite);
                leg->SetTextSize(leg_txt);
                leg->AddEntry(gp, "+ helicity", "PE");
                leg->AddEntry(gm, "- helicity", "PE");
                leg->Draw();

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.058);
                lab.SetTextAlign(13);

                const double q2m = safe_mean(C.q2means);
                const double tm  = safe_mean(C.tmeans);

                lab.DrawLatex(0.16,
                              0.88,
                              Form("Q^{2}=%.2f  |t|=%.2f", q2m, tm));
            }
        }

        const std::string fpath =
            (fs::path(base_dir) /
             ("plot_signal_yield_" +
              period_dir + "_xB_" +
              std::to_string((int)std::round(xb.first * 1000.0)) +
              ".png")).string();

        c->SaveAs(fpath.c_str());

        delete c;
    }
}

} // end anonymous namespace

bool update_pi0_corrected_counts_csv(const std::string& csv_path,
                                     const std::string& out_root_dir) {
    namespace fs = std::filesystem;

    const std::string csv_abs = fs::absolute(csv_path).string();

    std::error_code ec;

    const uintmax_t size_before =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[pi0_corrected] CSV: " << csv_abs
              << " (size=" << size_before << ")\n";

    CsvDoc csv;

    if (!csv.load(csv_path)) {
        return false;
    }

    if (!fill_signal_yields(csv)) {
        std::cerr << "[pi0_corrected] ERROR: fill_signal_yields failed.\n";
        return false;
    }

    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[pi0_corrected] ERROR: failed to save updated CSV.\n";
        return false;
    }

    const uintmax_t size_after =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[pi0_corrected] Updated CSV: " << csv_abs
              << " (size " << size_before << " -> " << size_after << ")\n";

    for (const auto& per : kPeriods) {
        draw_signal_yield_canvases(per, csv, out_root_dir);
    }

    std::cout << "[pi0_corrected] Signal-yield plotting finished.\n";

    return true;
}