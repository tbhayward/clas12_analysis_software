// unfolding.cpp
// Acceptance-corrected (unfolded) DVCS signal yields:
//   - Uses existing CSV columns:
//       signal yield, ep->epg, exp, <period>, <hel>  (triples)
//       acceptance, <period>                         (triples)
//   - Fills (for periods and combined groups):
//       acceptance corrected yield, ep->epg, exp, <label>, <hel>  (triples)
//     where <label> is one of:
//       Fa18 Inb, Fa18 Out, Sp19 Inb, Sp18 Inb, Sp18 Out,
//       Fa18, Sp18, 2018 (10.6 GeV)
//   - Periods are corrected bin by bin using:
//       U = S / A
//       Var(U) = (dU/dS)^2 Var(S) + (dU/dA)^2 Var(A)
//               = (1/A^2) sigma_S^2 + (S^2 / A^4) sigma_A^2
//   - Combined groups are simple sums of unfolded yields across member periods,
//     with variances added in quadrature.
//   - Produces xB-sliced canvases vs phi for pos/neg helicities (unpol skipped)
//     into output/unfolded_yields/<PeriodOrGroupDir>/.

#include "unfolding.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TH1D.h>
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
            std::cerr << "[unfolding] ERROR: cannot open CSV: " << path << "\n";
            return false;
        }
        std::string line;
        if (!std::getline(fin, line)) {
            std::cerr << "[unfolding] ERROR: empty CSV: " << path << "\n";
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
                std::cerr << "[unfolding] ERROR: cannot write CSV tmp: " << tmp << "\n";
                return false;
            }
            // header
            for (size_t i = 0; i < header.size(); ++i) {
                write_field(fout, header[i]);
                if (i + 1 < header.size()) fout << ',';
            }
            fout << "\n";
            // rows
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
                std::cerr << "[unfolding] ERROR: atomic rename failed ("
                          << ec.message() << ")\n";
                return false;
            }
        }
        return true;
    }

    int nrows() const { return (int)rows.size(); }

    bool has_col(const std::string& name) const {
        return index.find(name) != index.end();
    }

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

// basic helpers
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

// For periods and groups (Fa18, Sp18, 2018 (10.6 GeV)) -> directory name
static std::string canonical_dir_for_label(const std::string& L) {
    const std::string k = to_lower_nospace(L);
    if (k == "fa18inb" || k == "fa18out" ||
        k == "fa18inbsupp" || k == "sp18inb" ||
        k == "sp18out" || k == "sp19inb") {
        return canonical_period_dir(L);
    }
    if (k == "fa18") return "Fa18";
    if (k == "sp18") return "Sp18";
    // CSV label is "2018 (10.6 GeV)", mean columns use "10.6 GeV",
    // output dir should be "10.6_GeV"
    if (L == "2018 (10.6 GeV)" || k == "2018(10.6gev)") return "10.6_GeV";

    // fallback: spaces to underscore
    std::string s = L;
    for (char& c : s) {
        if (c == ' ') c = '_';
    }
    return s;
}

// For bin-mean columns we need the name that actually appears in the CSV header.
// Periods: use the label itself (Fa18 Inb, Fa18 Out, Sp18 Inb, Sp18 Out, Sp19 Inb)
// Groups:
//   Fa18  -> "Fa18"
//   Sp18  -> "Sp18"
//   2018 (10.6 GeV) -> "10.6 GeV"
static std::string mean_label_for_label(const std::string& L) {
    if (L == "2018 (10.6 GeV)") return "10.6 GeV";
    return L;
}

static void ensure_dir(const std::string& p) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(p, ec);
    if (ec) {
        std::cerr << "[unfolding] FATAL: could not create directory: "
                  << p << " (" << ec.message() << ")\n";
        std::exit(EXIT_FAILURE);
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

// Fixed triple parser for "(value, stat, sys)".
// Parses directly from the CSV cell without constructing a whitespace-stripped
// copy, token strings, or a temporary vector.
static bool parse_triple(const std::string& s,
                         double& value,
                         double& stat,
                         double& sys)
{
    value = std::numeric_limits<double>::quiet_NaN();
    stat  = std::numeric_limits<double>::quiet_NaN();
    sys   = std::numeric_limits<double>::quiet_NaN();

    const char* p = s.c_str();
    const char* const end = p + s.size();

    auto skip_space = [&]() {
        while (p < end && std::isspace(static_cast<unsigned char>(*p))) {
            ++p;
        }
    };

    auto parse_number = [&](double& out) -> bool {
        skip_space();
        char* parsed_end = nullptr;
        out = std::strtod(p, &parsed_end);
        if (parsed_end == p) {
            return false;
        }
        p = parsed_end;
        skip_space();
        return true;
    };

    skip_space();
    if (p == end) {
        return false;
    }

    const bool parenthesized = (*p == '(');
    if (parenthesized) {
        ++p;
    }

    if (!parse_number(value)) return false;
    if (p >= end || *p != ',') return false;
    ++p;

    if (!parse_number(stat)) return false;
    if (p >= end || *p != ',') return false;
    ++p;

    if (!parse_number(sys)) return false;

    if (parenthesized) {
        if (p >= end || *p != ')') return false;
        ++p;
        skip_space();
    }

    return p == end;
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

// U = S/A.
// Statistical propagation:
//   Var_stat(U) = (1/A^2) Var_stat(S) + (S^2/A^4) Var_stat(A)
// Systematic propagation:
//   Var_sys(U)  = (1/A^2) Var_sys(S)  + (S^2/A^4) Var_sys(A)
static void compute_unfolded(double S, double S_stat, double S_sys,
                             double A, double A_stat, double A_sys,
                             double& U, double& U_stat, double& U_sys)
{
    if (!std::isfinite(S) || !std::isfinite(S_stat) || !std::isfinite(S_sys) ||
        !std::isfinite(A) || !std::isfinite(A_stat) || !std::isfinite(A_sys) ||
        A <= 0.0) {
        U = std::numeric_limits<double>::quiet_NaN();
        U_stat = std::numeric_limits<double>::quiet_NaN();
        U_sys = std::numeric_limits<double>::quiet_NaN();
        return;
    }

    const double invA  = 1.0 / A;
    const double invA2 = invA * invA;
    const double invA4 = invA2 * invA2;

    const double varS_stat = S_stat * S_stat;
    const double varA_stat = A_stat * A_stat;
    const double varS_sys  = S_sys * S_sys;
    const double varA_sys  = A_sys * A_sys;

    double varU_stat = varS_stat * invA2 + (S * S) * varA_stat * invA4;
    double varU_sys  = varS_sys  * invA2 + (S * S) * varA_sys  * invA4;

    if (varU_stat < 0.0 && std::fabs(varU_stat) < 1.0e-15) varU_stat = 0.0;
    if (varU_sys  < 0.0 && std::fabs(varU_sys)  < 1.0e-15) varU_sys  = 0.0;

    U = S * invA;
    U_stat = (varU_stat > 0.0) ? std::sqrt(varU_stat) : 0.0;
    U_sys  = (varU_sys  > 0.0) ? std::sqrt(varU_sys)  : 0.0;
}

// periods, helicities, groups
static const std::vector<std::string> kPeriods = {
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out"
};

static const std::vector<std::string> kGroups = {
    "Fa18",
    "Sp18",
    "2018 (10.6 GeV)"
};

static const std::map<std::string,std::vector<std::string>> kGroupMembers = {
    {"Fa18", {"Fa18 Inb", "Fa18 Out"}},
    {"Sp18", {"Sp18 Inb", "Sp18 Out"}},
    {"2018 (10.6 GeV)", {"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"}}
};

static const std::vector<std::string> kHelicities = {
    "unpol",
    "pos",
    "neg"
};

static const std::vector<std::string> kUnpolarizedOnlyHelicities = {
    "unpol"
};

static bool has_helicity_resolved_columns(const std::string& label) {
    if (label == "Sp18 Inb" || label == "Sp18 Out" ||
        label == "Sp18" || label == "2018 (10.6 GeV)") {
        return false;
    }

    return true;
}

static const std::vector<std::string>& helicities_for_label(const std::string& label) {
    if (has_helicity_resolved_columns(label)) {
        return kHelicities;
    }

    return kUnpolarizedOnlyHelicities;
}

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

// Core CSV update:
//   - fills all "acceptance corrected yield, ep->epg, exp, <label>, <hel>" columns
//   - builds in-memory unfolded values/uncertainties for plotting.
static bool fill_unfolded_yields(
    CsvDoc& csv,
    std::map<std::string, std::map<std::string,std::vector<double>>>& unfolded_val,
    std::map<std::string, std::map<std::string,std::vector<double>>>& unfolded_stat,
    std::map<std::string, std::map<std::string,std::vector<double>>>& unfolded_sys)
{
    const int NR = csv.nrows();

    // indices for signal yield per period/helicity
    std::map<std::string, std::map<std::string,int>> sig_idx;
    for (const auto& per : kPeriods) {
        for (const auto& hel : helicities_for_label(per)) {
            std::ostringstream name;
            name << "signal yield, ep->epg, exp, " << per << ", " << hel;
            int idx = csv.col_index(name.str());
            if (idx < 0) {
                std::cerr << "[unfolding] FATAL: missing signal-yield column '"
                          << name.str() << "'. Did you run pi0_corrected_counts?\n";
                return false;
            }
            sig_idx[per][hel] = idx;
        }
    }

    // indices for acceptance triples per period
    std::map<std::string,int> acc_idx;
    for (const auto& per : kPeriods) {
        std::string cname = "acceptance, " + per;
        int idx = csv.col_index(cname);
        if (idx < 0) {
            std::cerr << "[unfolding] FATAL: missing acceptance column '"
                      << cname << "'. Did you run acceptance?\n";
            return false;
        }
        acc_idx[per] = idx;
    }

    // indices for unfolded-yield columns per label (periods + groups) and helicity
    std::map<std::string, std::map<std::string,int>> unf_idx;
    std::vector<std::string> allLabels;
    allLabels.insert(allLabels.end(), kPeriods.begin(), kPeriods.end());
    allLabels.insert(allLabels.end(), kGroups.begin(),  kGroups.end());

    for (const auto& lab : allLabels) {
        for (const auto& hel : helicities_for_label(lab)) {
            std::ostringstream nm;
            nm << "acceptance corrected yield, ep->epg, exp, "
               << lab << ", " << hel;
            int idx = csv.col_index(nm.str());
            if (idx < 0) {
                std::cerr << "[unfolding] FATAL: missing unfolded-yield column '"
                          << nm.str() << "' in CSV header.\n";
                return false;
            }
            unf_idx[lab][hel] = idx;
        }
    }

    // initialize in-memory storage: unfolded_val[label][hel][row]
    for (const auto& lab : allLabels) {
        for (const auto& hel : helicities_for_label(lab)) {
            unfolded_val[lab][hel].assign(NR,
                std::numeric_limits<double>::quiet_NaN());
            unfolded_stat[lab][hel].assign(NR,
                std::numeric_limits<double>::quiet_NaN());
            unfolded_sys[lab][hel].assign(NR,
                std::numeric_limits<double>::quiet_NaN());
        }
    }

    std::size_t cells_written = 0;

    // First: per-period unfolding S/A
    for (const auto& per : kPeriods) {
        const int c_acc = acc_idx[per];

        for (int r = 0; r < NR; ++r) {
            const std::string& sA = csv.rows[r][c_acc];

            if (sA.empty()) {
                // No acceptance for this bin/period: clear unfolded cells
                for (const auto& hel : helicities_for_label(per)) {
                    csv.rows[r][unf_idx[per][hel]].clear();
                    unfolded_val[per][hel][r]  =
                        std::numeric_limits<double>::quiet_NaN();
                    unfolded_stat[per][hel][r] =
                        std::numeric_limits<double>::quiet_NaN();
                    unfolded_sys[per][hel][r] =
                        std::numeric_limits<double>::quiet_NaN();
                }
                continue;
            }

            double A_val  = 0.0;
            double A_stat = 0.0;
            double A_sys  = 0.0;
            if (!parse_triple(sA, A_val, A_stat, A_sys)) {
                std::cerr << "[unfolding] FATAL: failed to parse acceptance '"
                          << sA << "' for period " << per
                          << " row " << r << "\n";
                return false;
            }

            if (!std::isfinite(A_val) || A_val <= 0.0) {
                // Non-physical or zero acceptance: clear unfolded cells
                for (const auto& hel : helicities_for_label(per)) {
                    csv.rows[r][unf_idx[per][hel]].clear();
                    unfolded_val[per][hel][r]  =
                        std::numeric_limits<double>::quiet_NaN();
                    unfolded_stat[per][hel][r] =
                        std::numeric_limits<double>::quiet_NaN();
                    unfolded_sys[per][hel][r] =
                        std::numeric_limits<double>::quiet_NaN();
                }
                continue;
            }

            for (const auto& hel : helicities_for_label(per)) {
                int c_sig = sig_idx[per][hel];
                const std::string& sS = csv.rows[r][c_sig];

                if (sS.empty()) {
                    // No signal yield: clear unfolded cell for this hel
                    csv.rows[r][unf_idx[per][hel]].clear();
                    unfolded_val[per][hel][r]  =
                        std::numeric_limits<double>::quiet_NaN();
                    unfolded_stat[per][hel][r] =
                        std::numeric_limits<double>::quiet_NaN();
                    unfolded_sys[per][hel][r] =
                        std::numeric_limits<double>::quiet_NaN();
                    continue;
                }

                double S_val  = 0.0;
                double S_stat = 0.0;
                double S_sys  = 0.0;
                if (!parse_triple(sS, S_val, S_stat, S_sys)) {
                    std::cerr << "[unfolding] FATAL: failed to parse signal-yield '"
                              << sS << "' for period " << per
                              << " helicity " << hel
                              << " row " << r << "\n";
                    return false;
                }

                double U = 0.0;
                double U_stat = 0.0;
                double U_sys = 0.0;
                compute_unfolded(S_val, S_stat, S_sys,
                                 A_val, A_stat, A_sys,
                                 U, U_stat, U_sys);

                unfolded_val[per][hel][r]  = U;
                unfolded_stat[per][hel][r] = U_stat;
                unfolded_sys[per][hel][r] = U_sys;

                csv.rows[r][unf_idx[per][hel]] = format_triple(U, U_stat, U_sys);
                ++cells_written;
            }
        }
    }

    // Next: combined groups = sums over member-period unfolded yields
    for (const auto& grp : kGroups) {
        auto itMembers = kGroupMembers.find(grp);
        if (itMembers == kGroupMembers.end()) {
            std::cerr << "[unfolding] FATAL: no member-period list for group '"
                      << grp << "'.\n";
            return false;
        }
        const std::vector<std::string>& members = itMembers->second;

        for (int r = 0; r < NR; ++r) {
            for (const auto& hel : helicities_for_label(grp)) {
                double sum_val = 0.0;
                double sum_stat_var = 0.0;
                double sum_sys_var = 0.0;
                int    count   = 0;

                for (const auto& per : members) {
                    double v = unfolded_val[per][hel][r];
                    double s = unfolded_stat[per][hel][r];
                    double y = unfolded_sys[per][hel][r];
                    if (!std::isfinite(v) || !std::isfinite(s) || !std::isfinite(y)) {
                        continue;
                    }
                    sum_val += v;
                    sum_stat_var += s * s;
                    sum_sys_var += y * y;
                    ++count;
                }

                if (count == 0) {
                    // No contributing periods in this bin: clear cell
                    csv.rows[r][unf_idx[grp][hel]].clear();
                    unfolded_val[grp][hel][r]  =
                        std::numeric_limits<double>::quiet_NaN();
                    unfolded_stat[grp][hel][r] =
                        std::numeric_limits<double>::quiet_NaN();
                    unfolded_sys[grp][hel][r] =
                        std::numeric_limits<double>::quiet_NaN();
                    continue;
                }

                const double U   = sum_val;
                const double U_s = (sum_stat_var > 0.0) ? std::sqrt(sum_stat_var) : 0.0;
                const double U_y = (sum_sys_var > 0.0) ? std::sqrt(sum_sys_var) : 0.0;

                unfolded_val[grp][hel][r]  = U;
                unfolded_stat[grp][hel][r] = U_s;
                unfolded_sys[grp][hel][r] = U_y;

                csv.rows[r][unf_idx[grp][hel]] = format_triple(U, U_s, U_y);
                ++cells_written;
            }
        }
    }

    std::cout << "[unfolding] Filled acceptance-corrected yield triple columns; cells written: "
              << cells_written << "\n";

    return true;
}

// Draw per-label unfolded-yield canvases (pos/neg helicities only)
static void draw_unfolded_canvases(
    const std::string& label,
    const CsvDoc& csv,
    const std::map<std::string,std::vector<double>>& val_by_hel,
    const std::map<std::string,std::vector<double>>& stat_by_hel,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const int NR = csv.nrows();

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
        std::cerr << "[unfolding] FATAL: missing bin-edge columns needed for plotting.\n";
        std::exit(EXIT_FAILURE);
    }

    // For mean values, use the actual label that appears in the header.
    const std::string mean_label = mean_label_for_label(label);

    int c_phiavg = csv.col_index("phiavg, " + mean_label);
    int c_q2avg  = csv.col_index("Q2avg, " + mean_label);
    int c_tabavg = csv.col_index("t_abs_avg, " + mean_label);
    int c_xbavg  = csv.col_index("xBavg, " + mean_label);

    if (c_xbavg < 0) {
        std::cerr << "[unfolding] FATAL: missing xBavg column for label '"
                  << label << "' (expected 'xBavg, " << mean_label << "').\n";
        std::exit(EXIT_FAILURE);
    }

    // We only plot pos/neg helicities
    auto itPosVal = val_by_hel.find("pos");
    auto itNegVal = val_by_hel.find("neg");
    auto itPosErr = stat_by_hel.find("pos");
    auto itNegErr = stat_by_hel.find("neg");

    if (itPosVal == val_by_hel.end() || itNegVal == val_by_hel.end() ||
        itPosErr == stat_by_hel.end() || itNegErr == stat_by_hel.end()) {
        std::cerr << "[unfolding] FATAL: internal error: missing pos/neg arrays for label '"
                  << label << "'.\n";
        std::exit(EXIT_FAILURE);
    }

    const std::vector<double>& posVal = itPosVal->second;
    const std::vector<double>& negVal = itNegVal->second;
    const std::vector<double>& posErr = itPosErr->second;
    const std::vector<double>& negErr = itNegErr->second;

    if ((int)posVal.size() != NR || (int)negVal.size() != NR ||
        (int)posErr.size() != NR || (int)negErr.size() != NR) {
        std::cerr << "[unfolding] FATAL: size mismatch in val/stat arrays for label '"
                  << label << "'.\n";
        std::exit(EXIT_FAILURE);
    }

    // distinct xB bins across the CSV
    std::set<std::pair<double,double>> xb_set;
    for (int r = 0; r < NR; ++r) {
        xb_set.emplace(csv.as_double(r, c_xb_min), csv.as_double(r, c_xb_max));
    }

    const double head_size = 0.14;
    const double label_sz  = 0.050;
    const double title_sz  = 0.045;
    const double leg_txt   = 0.050;

    const std::string dir_name  = canonical_dir_for_label(label);
    const std::string base_dir  =
        (fs::path(out_root_dir) / "unfolded_yields" / dir_name).string();
    ensure_dir(base_dir);

    for (auto xb : xb_set) {
        // group Q2 and t_abs bins within this xB bin
        std::set<std::pair<double,double>> q2set;
        std::set<std::pair<double,double>> tset_all;

        for (int r = 0; r < NR; ++r) {
            const double xbmin = csv.as_double(r, c_xb_min);
            const double xbmax = csv.as_double(r, c_xb_max);
            if (std::fabs(xbmin - xb.first) < 1e-9 &&
                std::fabs(xbmax - xb.second) < 1e-9) {
                q2set.emplace(csv.as_double(r, c_q2_min),
                              csv.as_double(r, c_q2_max));
            }
        }
        for (auto q2r : q2set) {
            for (int r = 0; r < NR; ++r) {
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
        double canvas_max = 1.0;
        double canvas_min = std::numeric_limits<double>::infinity();

        std::vector<double> xbmeans;
        for (int r = 0; r < NR; ++r) {
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

        // build per-cell datasets
        for (int rr = 0; rr < nrows; ++rr) {
            for (int cc = 0; cc < ncols; ++cc) {
                const auto& qpair = Q2s[rr];
                const auto& tpair = Ts[cc];

                std::vector<int> rows_for_cell;
                for (int r = 0; r < NR; ++r) {
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

                    const double Yp  = posVal[r];
                    const double Ym  = negVal[r];
                    const double Ep  = posErr[r];
                    const double Em  = negErr[r];

                    // Only use bins where both pos and neg are defined
                    if (!std::isfinite(Yp) || !std::isfinite(Ym) ||
                        !std::isfinite(Ep) || !std::isfinite(Em)) {
                        continue;
                    }

                    C.X.push_back(xphi);
                    C.Yp.push_back(Yp);
                    C.Ym.push_back(Ym);
                    C.EYp.push_back(Ep);
                    C.EYm.push_back(Em);
                    C.EX.push_back(0.0);

                    const double q2m = (c_q2avg >= 0)
                        ? csv.as_double(r, c_q2avg)
                        : 0.5 * (qpair.first + qpair.second);
                    const double tm  = (c_tabavg >= 0)
                        ? csv.as_double(r, c_tabavg)
                        : 0.5 * (tpair.first + tpair.second);
                    C.q2means.push_back(q2m);
                    C.tmeans.push_back(tm);

                    if (std::isfinite(Yp) && Yp > 0.0) {
                        canvas_max = std::max(canvas_max, Yp);
                        canvas_min = std::min(canvas_min, Yp);
                    }
                    if (std::isfinite(Ym) && Ym > 0.0) {
                        canvas_max = std::max(canvas_max, Ym);
                        canvas_min = std::min(canvas_min, Ym);
                    }
                }
            }
        }

        std::string cname = "c_unf_" + dir_name + "_xB_" +
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
                            label.c_str(),
                            xb_mean_for_title));

        // Common log-scale y-range per canvas
        double y_lo = 1.0;
        double y_hi = 1.0;
        if (std::isfinite(canvas_min) && canvas_min > 0.0) {
            y_lo = std::max(1.0, 0.5 * canvas_min);
        }
        if (canvas_max > 0.0) {
            y_hi = canvas_max * 1.15;
        }
        if (y_hi <= y_lo) {
            y_hi = y_lo * 10.0;
        }

        for (int rr = 0; rr < nrows; ++rr) {
            for (int cc = 0; cc < ncols; ++cc) {
                pGrid->cd(rr * ncols + cc + 1);
                gPad->SetGrid(1, 1);
                gPad->SetTopMargin(0.24);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.20);  // larger left margin to avoid label clipping
                gPad->SetRightMargin(0.07);
                gPad->SetTickx(1);
                gPad->SetTicky(1);
                gPad->SetLogy(1);           // log-scale y axis

                TH1* frame = gPad->DrawFrame(0.0, y_lo, 360.0, y_hi);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Unfolded yield");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(title_sz);
                frame->GetYaxis()->SetTitleSize(title_sz);
                frame->GetXaxis()->SetLabelSize(label_sz);
                frame->GetYaxis()->SetLabelSize(label_sz);
                frame->GetYaxis()->SetTitleOffset(1.8); // push label further left
                frame->GetXaxis()->SetTitleOffset(1.05);

                const CellData& C = cells[rr * ncols + cc];

                if (C.X.empty()) {
                    continue;
                }

                TGraphErrors* gp = new TGraphErrors(
                    (int)C.X.size(),
                    (double*)C.X.data(),
                    (double*)C.Yp.data(),
                    (double*)C.EX.data(),
                    (double*)C.EYp.data());

                TGraphErrors* gm = new TGraphErrors(
                    (int)C.X.size(),
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
                lab.DrawLatex(0.16, 0.88,
                              Form("Q^{2}=%.2f  |t|=%.2f", q2m, tm));
            }
        }

        const std::string fpath =
            (fs::path(base_dir) /
             ("plot_unfolded_yield_" +
              dir_name + "_xB_" +
              std::to_string((int)std::round(xb.first * 1000.0)) +
              ".png")).string();

        c->SaveAs(fpath.c_str());
        delete c;
    }
}


// -----------------------------------------------------------------------------
// Analysis-note summary outputs for the acceptance correction.
//
// The production calculation in this file is a bin-by-bin correction U=S/A.
// These plots intentionally use the phrase "acceptance-corrected yield" rather
// than implying a response-matrix unfolding that is not performed here.
// -----------------------------------------------------------------------------

static int note_period_color_unf(int ip) {
    const int colors[5] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kOrange+7};
    return colors[ip];
}

static double note_quantile_sorted_unf(const std::vector<double>& v, double q) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (q <= 0.0) return v.front();
    if (q >= 1.0) return v.back();
    const double x=q*(double)(v.size()-1);
    const size_t i0=(size_t)std::floor(x), i1=std::min(i0+1,v.size()-1);
    const double f=x-(double)i0;
    return v[i0]*(1.0-f)+v[i1]*f;
}

static bool note_parse_triple_cell(const CsvDoc& csv, int row, int col,
                                   double& v, double& s, double& y) {
    if (row<0 || row>=csv.nrows() || col<0 || col>=(int)csv.header.size()) return false;
    const std::string& cell=csv.rows[row][col];
    if (cell.empty()) return false;
    return parse_triple(cell,v,s,y);
}

static void write_unfolding_note_summary(
    const CsvDoc& csv,
    const std::map<std::string,std::map<std::string,std::vector<double>>>& unfolded_val,
    const std::map<std::string,std::map<std::string,std::vector<double>>>& unfolded_stat,
    const std::string& note_dir) {

    std::ofstream f((std::filesystem::path(note_dir)/"acceptance_correction_summary.csv").string());
    f << "period,n_bins,median_inverse_acceptance,p16_inverse_acceptance,p84_inverse_acceptance,median_signal_relative_stat,median_acceptance_relative_stat,median_corrected_relative_stat\n";
    f << std::setprecision(10);

    for (const auto& per:kPeriods) {
        const int ca=csv.col_index("acceptance, "+per);
        const int cs=csv.col_index("signal yield, ep->epg, exp, "+per+", unpol");
        std::vector<double> inva, rsig, racc, rout;
        for (int r=0;r<csv.nrows();++r) {
            double A=0,As=0,Ay=0,S=0,Ss=0,Sy=0;
            if (!note_parse_triple_cell(csv,r,ca,A,As,Ay) || !(A>0)) continue;
            inva.push_back(1.0/A);
            if (std::isfinite(As) && As>=0) racc.push_back(As/A);
            if (note_parse_triple_cell(csv,r,cs,S,Ss,Sy) && S!=0.0) rsig.push_back(std::fabs(Ss/S));
            const double U=unfolded_val.at(per).at("unpol")[r];
            const double Us=unfolded_stat.at(per).at("unpol")[r];
            if (std::isfinite(U) && U!=0.0 && std::isfinite(Us)) rout.push_back(std::fabs(Us/U));
        }
        std::sort(inva.begin(),inva.end()); std::sort(rsig.begin(),rsig.end());
        std::sort(racc.begin(),racc.end()); std::sort(rout.begin(),rout.end());
        f << per << ',' << inva.size() << ','
          << note_quantile_sorted_unf(inva,.50) << ',' << note_quantile_sorted_unf(inva,.16) << ',' << note_quantile_sorted_unf(inva,.84) << ','
          << note_quantile_sorted_unf(rsig,.50) << ',' << note_quantile_sorted_unf(racc,.50) << ',' << note_quantile_sorted_unf(rout,.50) << '\n';
    }
}

static void draw_unfolding_note_factor_distribution(const CsvDoc& csv,
                                                    const std::string& note_dir) {
    TCanvas c("c_note_accfac_dist","",1100,700);
    c.SetLeftMargin(0.11); c.SetRightMargin(0.035); c.SetBottomMargin(0.12); c.SetTopMargin(0.10); c.SetTicks(1,1);
    std::vector<TH1D*> hs; double ymax=0.0;
    for(int ip=0;ip<5;++ip){
        auto* h=new TH1D(Form("h_note_accfac_%d",ip),"",50,0,100);
        int ca=csv.col_index("acceptance, "+kPeriods[ip]);
        for(int r=0;r<csv.nrows();++r){double A=0,s=0,y=0;if(note_parse_triple_cell(csv,r,ca,A,s,y)&&A>0){double f=1.0/A;if(f<100.0)h->Fill(f);}}
        if(h->Integral()>0)h->Scale(1.0/h->Integral());
        h->SetLineColor(note_period_color_unf(ip));h->SetMarkerColor(note_period_color_unf(ip));h->SetMarkerStyle(20+ip);h->SetMarkerSize(.75);h->SetLineWidth(2);
        ymax=std::max(ymax,h->GetMaximum());hs.push_back(h);
    }
    TH1F frame("h_note_accfac_frame","",100,0,100);frame.SetMinimum(0);frame.SetMaximum(1.2*ymax);frame.GetXaxis()->SetTitle("Acceptance correction factor 1/A");frame.GetYaxis()->SetTitle("Fraction of populated analysis bins");frame.GetYaxis()->SetTitleOffset(1.15);frame.Draw();
    for(auto* h:hs)h->Draw("PE SAME");
    TLegend leg(0.72,0.60,0.94,0.86);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.036);for(int ip=0;ip<5;++ip)leg.AddEntry(hs[ip],kPeriods[ip].c_str(),"pe");leg.Draw();
    TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.048);t.DrawLatex(.11,.93,"Distribution of the acceptance correction factor");
    TLatex n;n.SetNDC();n.SetTextFont(42);n.SetTextSize(.030);n.DrawLatex(.12,.84,"Shown: 1/A < 100; the low-acceptance tail is retained in the numerical summary");
    c.SaveAs((std::filesystem::path(note_dir)/"acceptance_correction_factor_distribution.png").string().c_str());
    for(auto* h:hs)delete h;
}

static void draw_unfolding_note_uncertainty_summary(const CsvDoc& csv,
                                                    const std::string& note_dir) {
    TCanvas c("c_note_acc_unc","",1100,700);c.SetLeftMargin(.11);c.SetRightMargin(.035);c.SetBottomMargin(.16);c.SetTopMargin(.10);c.SetGridy();c.SetTicks(1,1);
    TH1F frame("h_note_acc_unc","",5,.5,5.5);frame.SetMinimum(0);frame.SetMaximum(.20);frame.GetYaxis()->SetTitle("Median fractional statistical uncertainty");frame.GetYaxis()->SetTitleOffset(1.10);frame.GetXaxis()->SetLabelSize(.043);for(int ip=0;ip<5;++ip)frame.GetXaxis()->SetBinLabel(ip+1,kPeriods[ip].c_str());frame.Draw();
    TGraphErrors gs,ga; gs.SetMarkerStyle(24);gs.SetMarkerSize(1.25);gs.SetMarkerColor(kGray+2);gs.SetLineColor(kGray+2); ga.SetMarkerStyle(20);ga.SetMarkerSize(1.25);ga.SetMarkerColor(kBlue+1);ga.SetLineColor(kBlue+1);
    for(int ip=0;ip<5;++ip){std::vector<double> rs,ra;int cs=csv.col_index("signal yield, ep->epg, exp, "+kPeriods[ip]+", unpol"),ca=csv.col_index("acceptance, "+kPeriods[ip]);for(int r=0;r<csv.nrows();++r){double S=0,Ss=0,Sy=0,A=0,As=0,Ay=0;if(note_parse_triple_cell(csv,r,cs,S,Ss,Sy)&&S!=0)rs.push_back(std::fabs(Ss/S));if(note_parse_triple_cell(csv,r,ca,A,As,Ay)&&A>0)ra.push_back(As/A);}std::sort(rs.begin(),rs.end());std::sort(ra.begin(),ra.end());gs.SetPoint(ip,ip+1-.08,note_quantile_sorted_unf(rs,.50));ga.SetPoint(ip,ip+1+.08,note_quantile_sorted_unf(ra,.50));}
    gs.Draw("P SAME");ga.Draw("P SAME");TLegend leg(.66,.68,.94,.84);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.036);leg.AddEntry(&gs,"Signal yield","p");leg.AddEntry(&ga,"Acceptance","p");leg.Draw();
    TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.048);t.DrawLatex(.11,.93,"Statistical precision entering the acceptance correction");
    c.SaveAs((std::filesystem::path(note_dir)/"acceptance_correction_uncertainty_summary.png").string().c_str());
}

static void draw_unfolding_note_phi_example(
    const CsvDoc& csv,
    const std::map<std::string,std::map<std::string,std::vector<double>>>& unfolded_val,
    const std::map<std::string,std::map<std::string,std::vector<double>>>& unfolded_stat,
    const std::string& note_dir) {

    const int cx0=csv.col_index("xBmin"),cx1=csv.col_index("xBmax"),cq0=csv.col_index("Q2min"),cq1=csv.col_index("Q2max"),ct0=csv.col_index("t_abs_min"),ct1=csv.col_index("t_abs_max"),cp0=csv.col_index("phimin"),cp1=csv.col_index("phimax");
    const double tx0=.204,tx1=.268;
    struct Key{double q0,q1,t0,t1;};std::vector<Key> keys;
    for(int r=0;r<csv.nrows();++r){if(csv.as_double(r,cx0)!=tx0||csv.as_double(r,cx1)!=tx1)continue;Key k{csv.as_double(r,cq0),csv.as_double(r,cq1),csv.as_double(r,ct0),csv.as_double(r,ct1)};bool seen=false;for(const auto& z:keys)if(z.q0==k.q0&&z.q1==k.q1&&z.t0==k.t0&&z.t1==k.t1){seen=true;break;}if(!seen)keys.push_back(k);}
    if(keys.empty())return;
    int best=-1,totalbest=-1;Key sel=keys.front();
    for(const auto& k:keys){int mn=100000,total=0;for(int ip=0;ip<5;++ip){int cs=csv.col_index("signal yield, ep->epg, exp, "+kPeriods[ip]+", unpol");int n=0;for(int r=0;r<csv.nrows();++r){if(csv.as_double(r,cx0)!=tx0||csv.as_double(r,cx1)!=tx1||csv.as_double(r,cq0)!=k.q0||csv.as_double(r,cq1)!=k.q1||csv.as_double(r,ct0)!=k.t0||csv.as_double(r,ct1)!=k.t1)continue;double S=0,Ss=0,Sy=0;if(note_parse_triple_cell(csv,r,cs,S,Ss,Sy)&&S>0&&std::isfinite(unfolded_val.at(kPeriods[ip]).at("unpol")[r])&&unfolded_val.at(kPeriods[ip]).at("unpol")[r]>0)++n;}mn=std::min(mn,n);total+=n;}if(mn>best||(mn==best&&total>totalbest)){best=mn;totalbest=total;sel=k;}}

    TCanvas c("c_note_acc_corr_phi","",1200,760);c.Divide(3,2,.001,.001);
    for(int ip=0;ip<5;++ip){c.cd(ip+1);gPad->SetLeftMargin(.15);gPad->SetRightMargin(.04);gPad->SetBottomMargin(.14);gPad->SetTopMargin(.16);gPad->SetGridy();gPad->SetTicks(1,1);gPad->SetLogy();int cs=csv.col_index("signal yield, ep->epg, exp, "+kPeriods[ip]+", unpol");int cpa=csv.col_index("phiavg, "+kPeriods[ip]);std::vector<double>x,ys,yu,es,eu;double ymin=1e300,ymax=0;
        for(int r=0;r<csv.nrows();++r){if(csv.as_double(r,cx0)!=tx0||csv.as_double(r,cx1)!=tx1||csv.as_double(r,cq0)!=sel.q0||csv.as_double(r,cq1)!=sel.q1||csv.as_double(r,ct0)!=sel.t0||csv.as_double(r,ct1)!=sel.t1)continue;double S=0,Ss=0,Sy=0;if(!note_parse_triple_cell(csv,r,cs,S,Ss,Sy)||!(S>0))continue;double U=unfolded_val.at(kPeriods[ip]).at("unpol")[r],Us=unfolded_stat.at(kPeriods[ip]).at("unpol")[r];if(!(U>0)||!std::isfinite(Us))continue;double phi=.5*(csv.as_double(r,cp0)+csv.as_double(r,cp1));if(cpa>=0){double z=csv.as_double(r,cpa);if(std::isfinite(z)&&z>0&&z<360)phi=z;}x.push_back(phi);ys.push_back(S);es.push_back(Ss);yu.push_back(U);eu.push_back(Us);ymin=std::min(ymin,std::min(S,U));ymax=std::max(ymax,std::max(S+Ss,U+Us));}
        if(!(ymin>0)||!(ymax>0))continue;TH1* frame=gPad->DrawFrame(0,std::max(1e-3,.55*ymin),360,1.35*ymax);frame->GetXaxis()->SetTitle("#phi (deg)");frame->GetYaxis()->SetTitle("Yield");frame->GetYaxis()->SetTitleOffset(1.25);frame->GetXaxis()->SetTitleSize(.055);frame->GetYaxis()->SetTitleSize(.055);frame->GetXaxis()->SetLabelSize(.045);frame->GetYaxis()->SetLabelSize(.045);std::vector<double>ex(x.size(),0.0);auto* gs=new TGraphErrors((int)x.size(),x.data(),ys.data(),ex.data(),es.data());auto* gu=new TGraphErrors((int)x.size(),x.data(),yu.data(),ex.data(),eu.data());gs->SetMarkerStyle(24);gs->SetMarkerColor(kGray+2);gs->SetLineColor(kGray+2);gu->SetMarkerStyle(20);gu->SetMarkerColor(note_period_color_unf(ip));gu->SetLineColor(note_period_color_unf(ip));gs->Draw("PE SAME");gu->Draw("PE SAME");auto* leg=new TLegend(.54,.67,.94,.84);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.038);leg->AddEntry(gs,"Signal yield","pe");leg->AddEntry(gu,"Acceptance corrected","pe");leg->Draw();TLatex lab;lab.SetNDC();lab.SetTextFont(42);lab.SetTextSize(.050);lab.DrawLatex(.15,.91,kPeriods[ip].c_str());}
    c.cd(6);gPad->SetFillStyle(0);TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.065);t.DrawLatex(.08,.84,"Representative acceptance correction");std::ostringstream kin;kin<<std::fixed<<std::setprecision(3)<<tx0<<" < x_{B} < "<<tx1<<"\\n"<<sel.q0<<" < Q^{2} < "<<sel.q1<<" GeV^{2}\\n"<<sel.t0<<" < |t| < "<<sel.t1<<" GeV^{2}";TLatex l;l.SetNDC();l.SetTextFont(42);l.SetTextSize(.050);l.DrawLatex(.08,.66,Form("%.3f < x_{B} < %.3f",tx0,tx1));l.DrawLatex(.08,.53,Form("%.3f < Q^{2} < %.3f GeV^{2}",sel.q0,sel.q1));l.DrawLatex(.08,.40,Form("%.3f < |t| < %.3f GeV^{2}",sel.t0,sel.t1));l.SetTextSize(.042);l.DrawLatex(.08,.22,"U=S/A; error bars are statistical");
    c.SaveAs((std::filesystem::path(note_dir)/"acceptance_corrected_yield_phi_example.png").string().c_str());
}

static void write_unfolding_analysis_note_outputs(
    const CsvDoc& csv,
    const std::map<std::string,std::map<std::string,std::vector<double>>>& unfolded_val,
    const std::map<std::string,std::map<std::string,std::vector<double>>>& unfolded_stat,
    const std::string& out_root_dir) {
    const std::string note_dir=(std::filesystem::path(out_root_dir)/"unfolded_yields"/"analysis_note").string();
    ensure_dir(note_dir);
    write_unfolding_note_summary(csv,unfolded_val,unfolded_stat,note_dir);
    draw_unfolding_note_factor_distribution(csv,note_dir);
    draw_unfolding_note_uncertainty_summary(csv,note_dir);
    draw_unfolding_note_phi_example(csv,unfolded_val,unfolded_stat,note_dir);
}

} // end anonymous namespace

bool update_unfolded_yields_csv(const std::string& csv_path,
                                const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const std::string csv_abs = fs::absolute(csv_path).string();
    std::error_code ec;
    const uintmax_t size_before =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[unfolding] CSV: " << csv_abs
              << " (size=" << size_before << ")\n";

    CsvDoc csv;
    if (!csv.load(csv_path)) {
        return false;
    }

    // in-memory unfolded yields for plotting
    std::map<std::string, std::map<std::string,std::vector<double>>> unfolded_val;
    std::map<std::string, std::map<std::string,std::vector<double>>> unfolded_stat;
    std::map<std::string, std::map<std::string,std::vector<double>>> unfolded_sys;

    if (!fill_unfolded_yields(csv, unfolded_val, unfolded_stat, unfolded_sys)) {
        std::cerr << "[unfolding] ERROR: fill_unfolded_yields failed.\n";
        return false;
    }

    // save CSV
    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[unfolding] ERROR: failed to save updated CSV.\n";
        return false;
    }

    const uintmax_t size_after =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[unfolding] Updated CSV: " << csv_abs
              << " (size " << size_before << " -> " << size_after << ")\n";

    // plotting: periods and combined groups
    std::vector<std::string> labels;
    labels.insert(labels.end(), kPeriods.begin(), kPeriods.end());
    labels.insert(labels.end(), kGroups.begin(),  kGroups.end());

    for (const auto& lab : labels) {
        if (!has_helicity_resolved_columns(lab)) {
            std::cout << "[unfolding] Skipping pos/neg unfolded-yield plots for unpolarized-only label: "
                      << lab << "\n";
            continue;
        }

        draw_unfolded_canvases(lab, csv,
                               unfolded_val.at(lab),
                               unfolded_stat.at(lab),
                               out_root_dir);
    }

    write_unfolding_analysis_note_outputs(csv, unfolded_val, unfolded_stat, out_root_dir);
    std::cout << "[unfolding] Analysis-note summary outputs written under: "
              << out_root_dir << "/unfolded_yields/analysis_note\n";

    std::cout << "[unfolding] Unfolded-yield plotting finished.\n";
    return true;
}