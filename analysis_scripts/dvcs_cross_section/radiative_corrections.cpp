// radiative_corrections.cpp
// Frad from generated Born and radiative MC, written into CSV as
// "(value, stat, sys)" triples in the columns:
//   "Frad, 10.6 GeV"   (Sp18 Inb + Sp18 Out + Fa18 Inb + Fa18 Out,
//                       with Fa18 Inb optionally skippable via a boolean)
//   "Frad, 10.2 GeV"   (Sp19 Inb)
//
// Binning comes from dvcs_pass2_analysis.csv (Lee-style).
//
// NEW BEHAVIOR:
//   - We now apply the *DVCS global kinematic cuts* from global_cuts.h
//     to the generated Born and Radiative MC before binning and
//     normalization:
//         * (-t1) < t1_abs_max
//         * open_angle_ep2 > open_angle_min_deg
//         * pTmiss <= pTmiss_max
//     via passes_global_cuts(t1, open_angle_ep2, pTmiss).
//
//   - The global normalizations N_born and N_rad used in
//       Frad_i = (a_i / N_born) / (b_i / N_rad)
//     are now defined as the total number of generated events that:
//         * pass the DVCS global cuts, and
//         * fall into any DVCS analysis bin (xB, Q2, |t|, phi).
//
// Frad is defined as Born/Rad with global MC normalizations per energy group:
//   Frad_i = (a_i / N_born) / (b_i / N_rad) = a_i * N_rad / (b_i * N_born)
//
// Also produces Frad vs phi canvases per beam energy and xB bin under
//   output/radiative_correction_plots/10.60
//   output/radiative_correction_plots/10.2
//
// Parallelization:
//   - 4 worker threads:
//       (10.6 GeV, Born), (10.6 GeV, Rad),
//       (10.2 GeV, Born), (10.2 GeV, Rad)
//   - Each worker loops over its MC trees, applies DVCS global cuts,
//     accumulates row-wise counts and a global N_total (post-cut, in-bin),
//     then returns. CSV writes and plotting occur only after all workers
//     complete.

#include "radiative_corrections.h"
#include "global_cuts.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
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
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cstdio>

namespace {

struct CsvDoc {
    std::vector<std::string> header;
    std::unordered_map<std::string,int> index;
    std::vector<std::vector<std::string> > rows;

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
        char* e = 0;
        double v = std::strtod(s.c_str(), &e);
        if (e == s.c_str()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return v;
    }

    bool load(const std::string& path) {
        std::ifstream fin(path.c_str());
        if (!fin.is_open()) {
            std::cerr << "[radcorr] ERROR: cannot open CSV: " << path << "\n";
            return false;
        }
        std::string line;
        if (!std::getline(fin, line)) {
            std::cerr << "[radcorr] ERROR: empty CSV: " << path << "\n";
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
        for (std::size_t r = 0; r < rows.size(); ++r) {
            rows[r].resize(header.size());
        }
        return true;
    }

    bool save_atomic(const std::string& path) const {
        const std::string tmp = path + ".tmp";
        {
            std::ofstream fout(tmp.c_str());
            if (!fout.is_open()) {
                std::cerr << "[radcorr] ERROR: cannot write CSV tmp: " << tmp << "\n";
                return false;
            }
            for (std::size_t i = 0; i < header.size(); ++i) {
                write_field(fout, header[i]);
                if (i + 1 < header.size()) fout << ',';
            }
            fout << "\n";
            for (std::size_t r = 0; r < rows.size(); ++r) {
                const std::vector<std::string>& row = rows[r];
                for (std::size_t i = 0; i < row.size(); ++i) {
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
                std::cerr << "[radcorr] ERROR: atomic rename failed ("
                          << ec.message() << ")\n";
                return false;
            }
        }
        return true;
    }

    int nrows() const { return (int)rows.size(); }

    int col_index(const std::string& name) const {
        std::unordered_map<std::string,int>::const_iterator it = index.find(name);
        return (it == index.end()) ? -1 : it->second;
    }

    double as_double(int r, int c) const {
        if (r < 0 || r >= nrows()) return std::numeric_limits<double>::quiet_NaN();
        if (c < 0 || c >= (int)header.size()) return std::numeric_limits<double>::quiet_NaN();
        return to_double(rows[r][c]);
    }
};

// basic helpers

static inline double PI() { return 3.14159265358979323846; }
static inline double RAD2DEG(double r) { return r * 180.0 / PI(); }

static inline double wrap_deg(double phi) {
    double x = std::fmod(phi, 360.0);
    if (x < 0.0) x += 360.0;
    return x;
}

static double safe_mean(const std::vector<double>& v) {
    double s = 0.0;
    int n = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        double x = v[i];
        if (std::isfinite(x)) {
            s += x;
            ++n;
        }
    }
    return n ? s / n : std::numeric_limits<double>::quiet_NaN();
}

// "(value, stat, sys)" formatting/parsing

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
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (!std::isspace((unsigned char)c)) trimmed.push_back(c);
    }
    if (trimmed.empty()) return false;

    if (trimmed.front() == '(' && trimmed.back() == ')') {
        trimmed = trimmed.substr(1, trimmed.size() - 2);
    }

    std::vector<std::string> parts;
    std::string cur;
    for (std::size_t i = 0; i < trimmed.size(); ++i) {
        char c = trimmed[i];
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

    char* e1 = 0;
    char* e2 = 0;
    char* e3 = 0;

    value = std::strtod(parts[0].c_str(), &e1);
    stat  = std::strtod(parts[1].c_str(), &e2);
    sys   = std::strtod(parts[2].c_str(), &e3);

    if (e1 == parts[0].c_str()) return false;
    if (e2 == parts[1].c_str()) return false;
    if (e3 == parts[2].c_str()) return false;

    return true;
}

// binning from CSV: one RowBin per phi-binned row

struct RowBin {
    double xbmin, xbmax;
    double q2min, q2max;
    double tmin,  tmax;
    double phimin, phimax;
    int    row_index;
};

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
        std::cerr << "[radcorr] FATAL: missing one or more bin-edge columns "
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
        if (!std::isfinite(phi_width)) continue;

        // Skip phi-integrated bins (we only define Frad for phi-binned rows)
        if (std::fabs(phi_width) >= 359.0) continue;

        if (!(b.xbmax > b.xbmin &&
              b.q2max > b.q2min &&
              b.tmax  > b.tmin  &&
              b.phimax > b.phimin)) {
            continue;
        }

        bins.push_back(b);
    }

    std::cout << "[radcorr] Built " << bins.size()
              << " RowBin entries out of " << NR << " CSV rows "
              << "(skipped phi-integrated or invalid bins).\n";

    if (bins.empty()) {
        std::cerr << "[radcorr] FATAL: no usable bins found in CSV.\n";
        std::exit(EXIT_FAILURE);
    }

    return bins;
}

// For each energy group, mark which CSV rows actually have data.
static std::vector<bool>
build_row_has_data(const CsvDoc& csv, const std::string& xbavg_col_name)
{
    const int NR = csv.nrows();
    int c = csv.col_index(xbavg_col_name);
    if (c < 0) {
        std::cerr << "[radcorr] FATAL: missing column '" << xbavg_col_name
                  << "' needed to decide where to compute Frad.\n";
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
    return flags;
}

// Find which CSV row an event belongs to, given kinematics in degrees.
static int find_row_for_event(double xB,
                              double Q2,
                              double tAbs,
                              double phiDeg,
                              const std::vector<RowBin>& bins)
{
    for (std::size_t i = 0; i < bins.size(); ++i) {
        const RowBin& b = bins[i];
        if (!(xB   >= b.xbmin  && xB   < b.xbmax))  continue;
        if (!(Q2   >= b.q2min  && Q2   < b.q2max))  continue;
        if (!(tAbs >= b.tmin   && tAbs < b.tmax))   continue;
        if (!(phiDeg >= b.phimin && phiDeg < b.phimax)) continue;
        return b.row_index;
    }
    return -1;
}

// MC accumulation with DVCS global cuts and progress prints.
// Returns N_accepted = total number of entries that:
//   - pass passes_global_cuts, and
//   - fall into any DVCS analysis bin.
static double accumulate_counts_for_tree(const std::string& group_label,
                                         const std::string& tree_label,
                                         TTree* tree,
                                         const std::vector<RowBin>& bins,
                                         const std::vector<bool>& row_has_data,
                                         std::vector<double>& counts)
{
    if (!tree) {
        std::cerr << "[radcorr] FATAL: null TTree pointer for "
                  << group_label << " (" << tree_label << ").\n";
        std::exit(EXIT_FAILURE);
    }

    const int NR = (int)row_has_data.size();
    if ((int)counts.size() != NR) {
        std::cerr << "[radcorr] FATAL: counts vector size mismatch in accumulate_counts_for_tree "
                  << "for " << group_label << " (" << tree_label << ").\n";
        std::exit(EXIT_FAILURE);
    }

    const char* br_x        = "x";
    const char* br_Q2       = "Q2";
    const char* br_t1       = "t1";
    const char* br_phi2     = "phi2";
    const char* br_oa_ep2   = "open_angle_ep2";
    const char* br_pTmiss   = "pTmiss";

    // Required branches in generated MC:
    //   - x, Q2, t1, phi2 for binning
    //   - open_angle_ep2, pTmiss for DVCS global cuts
    if (!tree->GetBranch(br_x) ||
        !tree->GetBranch(br_Q2) ||
        !tree->GetBranch(br_t1) ||
        !tree->GetBranch(br_phi2) ||
        !tree->GetBranch(br_oa_ep2) ||
        !tree->GetBranch(br_pTmiss)) {
        std::cerr << "[radcorr] FATAL: missing one or more branches in tree for "
                  << group_label << " (" << tree_label
                  << ") (expected: x, Q2, t1, phi2, open_angle_ep2, pTmiss).\n";
        std::exit(EXIT_FAILURE);
    }

    double g_x            = 0.0;
    double g_Q2           = 0.0;
    double g_t1           = 0.0;
    double g_phi2         = 0.0;
    double g_open_angle   = 0.0;
    double g_pTmiss       = 0.0;

    tree->SetBranchAddress(br_x,          &g_x);
    tree->SetBranchAddress(br_Q2,         &g_Q2);
    tree->SetBranchAddress(br_t1,         &g_t1);
    tree->SetBranchAddress(br_phi2,       &g_phi2);
    tree->SetBranchAddress(br_oa_ep2,     &g_open_angle);
    tree->SetBranchAddress(br_pTmiss,     &g_pTmiss);

    const Long64_t N = tree->GetEntries();
    Long64_t used = 0;

    int next_pct = 10;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        // Apply DVCS global kinematic cuts at generator level:
        //   (-t1) < t1_abs_max
        //   open_angle_ep2 > open_angle_min_deg
        //   pTmiss <= pTmiss_max
        if (!passes_global_cuts(g_t1, g_open_angle, g_pTmiss)) {
            // Still count toward progress percentage, but not toward N_accepted.
            if (N > 0 && next_pct <= 100) {
                double pct = 100.0 * (double)(i + 1) / (double)N;
                while (pct >= (double)next_pct && next_pct <= 100) {
                    std::cout << "[radcorr] Group " << group_label
                              << ", tree " << tree_label
                              << " progress: " << (double)next_pct << "% ("
                              << (long long)(i + 1) << "/"
                              << (long long)N << ")\n";
                    next_pct += 10;
                }
            }
            continue;
        }

        const double xB   = g_x;
        const double Q2   = g_Q2;
        const double tAbs = std::fabs(g_t1);
        const double phiD = wrap_deg(RAD2DEG(g_phi2));

        int row = find_row_for_event(xB, Q2, tAbs, phiD, bins);
        if (row < 0 || row >= NR) {
            if (N > 0 && next_pct <= 100) {
                double pct = 100.0 * (double)(i + 1) / (double)N;
                while (pct >= (double)next_pct && next_pct <= 100) {
                    std::cout << "[radcorr] Group " << group_label
                              << ", tree " << tree_label
                              << " progress: " << (double)next_pct << "% ("
                              << (long long)(i + 1) << "/"
                              << (long long)N << ")\n";
                    next_pct += 10;
                }
            }
            continue;
        }

        // NOTE: we *still* do not gate on row_has_data[row] here.
        // Generator MC should be counted for any analysis bin, regardless of
        // whether that energy group has reconstructed data in that bin.

        counts[row] += 1.0;
        ++used;

        if (N > 0 && next_pct <= 100) {
            double pct = 100.0 * (double)(i + 1) / (double)N;
            while (pct >= (double)next_pct && next_pct <= 100) {
                std::cout << "[radcorr] Group " << group_label
                          << ", tree " << tree_label
                          << " progress: " << (double)next_pct << "% ("
                          << (long long)(i + 1) << "/"
                          << (long long)N << ")\n";
                next_pct += 10;
            }
        }
    }

    std::cout << "[radcorr] Group " << group_label
              << ", tree " << tree_label
              << ": total entries = " << (long long)N
              << " ; passed DVCS global cuts and binned (any DVCS bin) = "
              << (long long)used << "\n";

    // Return N_accepted (post-cuts, in-bin) for use as N_born / N_rad
    // global normalizations at the group level.
    return (double)used;
}

// Fill one Frad column for a given energy group using global N_born, N_rad.
static bool fill_frad_for_group(const std::string& group_label,
                                const std::string& frad_col_name,
                                CsvDoc& csv,
                                const std::vector<bool>& row_has_data,
                                const std::vector<double>& born_counts,
                                const std::vector<double>& rad_counts,
                                double N_born,
                                double N_rad)
{
    const int NR = csv.nrows();
    if ((int)row_has_data.size() != NR ||
        (int)born_counts.size()   != NR ||
        (int)rad_counts.size()    != NR) {
        std::cerr << "[radcorr] FATAL: size mismatch in fill_frad_for_group for "
                  << group_label << ".\n";
        return false;
    }

    if (N_born <= 0.0 || N_rad <= 0.0) {
        std::cerr << "[radcorr] FATAL: non-positive global N_born or N_rad in fill_frad_for_group for "
                  << group_label << " (N_born=" << N_born
                  << ", N_rad=" << N_rad << ").\n";
        return false;
    }

    const int c_frad = csv.col_index(frad_col_name);
    if (c_frad < 0) {
        std::cerr << "[radcorr] FATAL: missing Frad column '"
                  << frad_col_name << "' in CSV header.\n";
        return false;
    }

    // Column indices for kinematic bin edges (for debug printing).
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
        std::cerr << "[radcorr] FATAL: missing bin-edge columns in fill_frad_for_group for "
                  << group_label << ".\n";
        return false;
    }

    std::size_t cells_written = 0;
    int debug_printed = 0;
    const int debug_limit = 20;

    std::cout << "[radcorr] Debug: per-bin Frad construction for group "
              << group_label << " (up to " << debug_limit << " bins):\n";
    std::cout << std::fixed << std::setprecision(4);

    for (int r = 0; r < NR; ++r) {
        if (!row_has_data[r]) continue;

        const double a = born_counts[r];
        const double b = rad_counts[r];

        if (a <= 0.0 || b <= 0.0) {
            // No usable information in this row
            continue;
        }

        const double RC = (a * N_rad) / (b * N_born); // Born/Rad with global norms

        if (!std::isfinite(RC)) {
            // Do not store NaN or inf
            continue;
        }

        double rel_var = 0.0;
        const double a_eff = std::max(a, 1.0);
        const double b_eff = std::max(b, 1.0);
        rel_var = 1.0 / a_eff + 1.0 / b_eff;
        double sRC = 0.0;
        if (rel_var > 0.0) {
            sRC = std::fabs(RC) * std::sqrt(rel_var);
        }

        // Debug print for the first debug_limit bins where we actually write Frad.
        if (debug_printed < debug_limit) {
            const double xbmin  = csv.as_double(r, c_xb_min);
            const double xbmax  = csv.as_double(r, c_xb_max);
            const double q2min  = csv.as_double(r, c_q2_min);
            const double q2max  = csv.as_double(r, c_q2_max);
            const double tmin   = csv.as_double(r, c_tab_min);
            const double tmax   = csv.as_double(r, c_tab_max);
            const double phimin = csv.as_double(r, c_phi_min);
            const double phimax = csv.as_double(r, c_phi_max);

            std::cout << "[radcorr]   row " << r
                      << " : xB=[" << xbmin << "," << xbmax << "]"
                      << " Q2=[" << q2min << "," << q2max << "]"
                      << " |t|=[" << tmin << "," << tmax << "]"
                      << " phi=[" << phimin << "," << phimax << "]"
                      << " ; a=" << a
                      << " b=" << b
                      << " N_born=" << N_born
                      << " N_rad=" << N_rad
                      << " => Frad=" << RC
                      << " +/- " << sRC
                      << "\n";
            ++debug_printed;
        }

        csv.rows[r][c_frad] = format_triple(RC, sRC, 0.0);
        ++cells_written;
    }

    std::cout << "[radcorr] Group " << group_label
              << ": filled Frad column '" << frad_col_name
              << "' with global N_born=" << N_born
              << " N_rad=" << N_rad
              << " ; cells written: " << cells_written << "\n";

    return true;
}

// plotting

struct CellData {
    std::vector<double> X;    // phi position (deg)
    std::vector<double> Y;    // Frad
    std::vector<double> EXL;  // left x error (deg)
    std::vector<double> EXH;  // right x error (deg)
    std::vector<double> EY;   // y error (stat, symmetric)
    std::vector<double> q2means;
    std::vector<double> tmeans;
};

static void ensure_dir(const std::string& p) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (!fs::exists(p)) {
        fs::create_directories(p, ec);
        if (ec) {
            std::cerr << "[radcorr] FATAL: could not create directory: "
                      << p << " (" << ec.message() << ")\n";
            std::exit(EXIT_FAILURE);
        }
    }
}

static void draw_radiative_canvases(const std::string& group_label,
                                    const std::string& xbavg_col,
                                    const std::string& phiavg_col,
                                    const std::string& q2avg_col,
                                    const std::string& tabavg_col,
                                    const std::string& frad_col_name,
                                    const CsvDoc& csv,
                                    const std::vector<bool>& row_has_data,
                                    const std::string& out_root_dir,
                                    const std::string& energy_dir)
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
        std::cerr << "[radcorr] FATAL: missing bin-edge columns needed for plotting.\n";
        std::exit(EXIT_FAILURE);
    }

    const int c_phiavg = csv.col_index(phiavg_col);
    const int c_q2avg  = csv.col_index(q2avg_col);
    const int c_tabavg = csv.col_index(tabavg_col);
    const int c_xbavg  = csv.col_index(xbavg_col);

    const int c_frad = csv.col_index(frad_col_name);
    if (c_frad < 0) {
        std::cerr << "[radcorr] FATAL: missing column '" << frad_col_name
                  << "' for plotting.\n";
        std::exit(EXIT_FAILURE);
    }

    std::set<std::pair<double,double> > xb_set;
    for (int r = 0; r < csv.nrows(); ++r) {
        if (!row_has_data[r]) continue;
        const double xbmin = csv.as_double(r, c_xb_min);
        const double xbmax = csv.as_double(r, c_xb_max);
        xb_set.insert(std::make_pair(xbmin, xbmax));
    }

    const double head_size = 0.14;
    const double label_sz  = 0.050;
    const double title_sz  = 0.045;

    const std::string base_dir =
        (fs::path(out_root_dir) / "radiative_correction_plots" / energy_dir).string();
    ensure_dir(base_dir);

    for (std::set<std::pair<double,double> >::const_iterator itX = xb_set.begin();
         itX != xb_set.end(); ++itX) {

        const std::pair<double,double> xb = *itX;

        std::set<std::pair<double,double> > q2set;
        std::set<std::pair<double,double> > tset_all;

        for (int r = 0; r < csv.nrows(); ++r) {
            if (!row_has_data[r]) continue;
            const double xbmin = csv.as_double(r, c_xb_min);
            const double xbmax = csv.as_double(r, c_xb_max);
            if (std::fabs(xbmin - xb.first) < 1e-9 &&
                std::fabs(xbmax - xb.second) < 1e-9) {
                q2set.insert(std::make_pair(csv.as_double(r, c_q2_min),
                                            csv.as_double(r, c_q2_max)));
            }
        }
        for (std::set<std::pair<double,double> >::const_iterator itQ = q2set.begin();
             itQ != q2set.end(); ++itQ) {
            const std::pair<double,double> q2r = *itQ;
            for (int r = 0; r < csv.nrows(); ++r) {
                if (!row_has_data[r]) continue;
                const double xbmin = csv.as_double(r, c_xb_min);
                const double xbmax = csv.as_double(r, c_xb_max);
                const double q2min = csv.as_double(r, c_q2_min);
                const double q2max = csv.as_double(r, c_q2_max);
                if (std::fabs(xbmin - xb.first) < 1e-9 &&
                    std::fabs(xbmax - xb.second) < 1e-9 &&
                    std::fabs(q2min - q2r.first) < 1e-9 &&
                    std::fabs(q2max - q2r.second) < 1e-9) {
                    tset_all.insert(std::make_pair(csv.as_double(r, c_tab_min),
                                                   csv.as_double(r, c_tab_max)));
                }
            }
        }

        std::vector<std::pair<double,double> > Q2s(q2set.begin(), q2set.end());
        std::vector<std::pair<double,double> > Ts (tset_all.begin(), tset_all.end());
        if (Q2s.empty() || Ts.empty()) continue;

        const int nrows = (int)Q2s.size();
        const int ncols = (int)Ts.size();
        const int W = 300 * ncols + 160;
        const int H = 260 * nrows + 240;

        std::vector<CellData> cells(nrows * ncols);
        double canvas_max = -std::numeric_limits<double>::infinity();
        double canvas_min =  std::numeric_limits<double>::infinity();

        std::vector<double> xbmeans;
        for (int r = 0; r < csv.nrows(); ++r) {
            if (!row_has_data[r]) continue;
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
                const std::pair<double,double>& qpair = Q2s[rr];
                const std::pair<double,double>& tpair = Ts[cc];

                std::vector<int> rows_for_cell;
                for (int r = 0; r < csv.nrows(); ++r) {
                    if (!row_has_data[r]) continue;
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
                C.EXL.reserve(rows_for_cell.size());
                C.EXH.reserve(rows_for_cell.size());
                C.EY.reserve(rows_for_cell.size());

                for (std::size_t k = 0; k < rows_for_cell.size(); ++k) {
                    int r = rows_for_cell[k];

                    const double pmin = csv.as_double(r, c_phi_min);
                    const double pmax = csv.as_double(r, c_phi_max);
                    double xphi = 0.5 * (pmin + pmax);
                    if (c_phiavg >= 0) {
                        const double pav = csv.as_double(r, c_phiavg);
                        if (std::isfinite(pav) && pav > 0.0 && pav < 360.0) {
                            xphi = pav;
                        }
                    }

                    const std::string& frad_cell = csv.rows[r][c_frad];
                    if (frad_cell.empty()) continue;

                    double frad_val = 0.0;
                    double frad_stat = 0.0;
                    double frad_sys = 0.0;
                    if (!parse_triple(frad_cell, frad_val, frad_stat, frad_sys)) {
                        double tmp = CsvDoc::to_double(frad_cell);
                        if (!std::isfinite(tmp)) continue;
                        frad_val  = tmp;
                        frad_stat = 0.0;
                    }

                    if (!std::isfinite(frad_val)) continue;

                    // Horizontal errors: asymmetric, to show full bin width around phiavg.
                    double exl = xphi - pmin;
                    double exh = pmax - xphi;
                    if (exl < 0.0) exl = 0.0;
                    if (exh < 0.0) exh = 0.0;

                    C.X.push_back(xphi);
                    C.Y.push_back(frad_val);
                    C.EY.push_back(frad_stat);
                    C.EXL.push_back(exl);
                    C.EXH.push_back(exh);

                    const double q2m = (c_q2avg >= 0)
                        ? csv.as_double(r, c_q2avg)
                        : 0.5 * (qpair.first + qpair.second);
                    const double tm  = (c_tabavg >= 0)
                        ? csv.as_double(r, c_tabavg)
                        : 0.5 * (tpair.first + tpair.second);
                    C.q2means.push_back(q2m);
                    C.tmeans.push_back(tm);

                    if (frad_val > canvas_max) canvas_max = frad_val;
                    if (frad_val < canvas_min) canvas_min = frad_val;
                }
            }
        }

        // If no points for this xB, skip.
        if (!std::isfinite(canvas_max) || !std::isfinite(canvas_min)) {
            continue;
        }

        // Build common y-range for all pads on this xB canvas.
        double y_min = canvas_min;
        double y_max = canvas_max;

        if (y_min == y_max) {
            double base = (y_min == 0.0) ? 1.0 : std::fabs(y_min);
            y_min -= 0.1 * base;
            y_max += 0.1 * base;
        }

        double span = y_max - y_min;
        if (span <= 0.0) span = 1.0;
        double pad = 0.15 * span;
        double y_lo = y_min - pad;
        double y_hi = y_max + pad;

        std::ostringstream cname;
        cname << "c_frad_" << energy_dir << "_xB_"
              << (int)std::round(xb.first * 1000.0);
        TCanvas* c = new TCanvas(cname.str().c_str(), "", W, H);

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
                       Form("Frad   %s   <xB>=%.2f",
                            group_label.c_str(),
                            xb_mean_for_title));

        for (int rr = 0; rr < nrows; ++rr) {
            for (int cc = 0; cc < ncols; ++cc) {
                pGrid->cd(rr * ncols + cc + 1);
                gPad->SetGrid(1, 1);
                gPad->SetTopMargin(0.24);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.160);
                gPad->SetRightMargin(0.07);
                gPad->SetTickx(1);
                gPad->SetTicky(1);

                TH1* frame = gPad->DrawFrame(0.0, y_lo, 360.0, y_hi);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("F_{rad}");
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

                TGraphAsymmErrors* gfrad = new TGraphAsymmErrors(
                    (int)C.X.size(),
                    (double*)C.X.data(),
                    (double*)C.Y.data(),
                    (double*)C.EXL.data(),
                    (double*)C.EXH.data(),
                    (double*)C.EY.data(),   // EY low
                    (double*)C.EY.data()    // EY high (symmetric)
                );

                gfrad->SetMarkerStyle(20);
                gfrad->SetMarkerColor(kBlack);
                gfrad->SetLineColor(kBlack);
                gfrad->SetLineWidth(1);

                gfrad->Draw("PE1 SAME");

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
             ("plot_frad_" + energy_dir + "_xB_" +
              std::to_string((int)std::round(xb.first * 1000.0)) +
              ".png")).string();

        c->SaveAs(fpath.c_str());
        delete c;
    }
}

// group definition

struct McPair {
    std::string bornTag;
    std::string radTag;
    std::string period_label; // just for print messages
};

struct RCGroup {
    std::string label;        // "10.6 GeV" or "10.2 GeV"
    std::string xbavg_col;    // "xBavg, 10.6 GeV" or "xBavg, Sp19 Inb"
    std::string phiavg_col;   // "phiavg, 10.6 GeV" or "phiavg, Sp19 Inb"
    std::string q2avg_col;    // "Q2avg, 10.6 GeV" or "Q2avg, Sp19 Inb"
    std::string tabavg_col;   // "t_abs_avg, 10.6 GeV" or "t_abs_avg, Sp19 Inb"
    std::string frad_col;     // "Frad, 10.6 GeV" or "Frad, 10.2 GeV"
    std::string energy_dir;   // "10.60" or "10.2"
    std::vector<McPair> mc_pairs;
};

struct GroupAccum {
    std::vector<double> born_counts;
    std::vector<double> rad_counts;
    double N_born;
    double N_rad;
    GroupAccum() : N_born(0.0), N_rad(0.0) {}
};

// Worker function for a single (group, type) pair.
// isBorn == true  -> use genMcTrees (Born), fill born_counts and N_born.
// isBorn == false -> use radGenMcTrees (Rad), fill rad_counts and N_rad.
static void accumulate_group_type(const RCGroup& G,
                                  const std::vector<RowBin>& bins,
                                  const std::vector<bool>& row_has_data,
                                  const std::map<std::string, TTree*>& treeMap,
                                  bool isBorn,
                                  std::vector<double>& counts_out,
                                  double& Ntotal_out)
{
    const std::string type_label = isBorn ? "Born" : "Rad";
    const int NR = (int)row_has_data.size();

    counts_out.assign(NR, 0.0);
    Ntotal_out = 0.0;

    for (std::size_t ip = 0; ip < G.mc_pairs.size(); ++ip) {
        const McPair& P = G.mc_pairs[ip];
        const std::string tag = isBorn ? P.bornTag : P.radTag;

        std::map<std::string,TTree*>::const_iterator itT = treeMap.find(tag);
        if (itT == treeMap.end()) {
            std::cerr << "[radcorr] FATAL: missing MC tree with tag '"
                      << tag << "' for group " << G.label
                      << " (" << type_label << ", period " << P.period_label << ").\n";
            std::exit(EXIT_FAILURE);
        }

        TTree* tree = itT->second;
        const std::string tree_label = P.period_label + " " + type_label;

        double N_here = accumulate_counts_for_tree(G.label,
                                                   tree_label,
                                                   tree,
                                                   bins,
                                                   row_has_data,
                                                   counts_out);
        Ntotal_out += N_here;
    }

    std::cout << "[radcorr] Group " << G.label << " " << type_label
              << ": global N_total (all periods, post-cuts, in-bin) = "
              << Ntotal_out << "\n";
}

} // end anonymous namespace

// public entry point

bool update_radiative_corrections_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& radGenMcTrees,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    // -------------------------------------------------------------------------
    // Temporary toggle:
    //   - Set SKIP_FA18_INB_10P6 = true to *omit* the Fa18 Inb period from the
    //     10.6 GeV Frad calculation (both Born and Rad).
    //   - When your Fa18 Inb MC is fixed, set this back to false to restore
    //     the full four-period combination.
    // -------------------------------------------------------------------------
    const bool SKIP_FA18_INB_10P6 = true;

    const std::string csv_abs = fs::absolute(csv_path).string();
    std::error_code ec;
    const uintmax_t size_before =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[radcorr] CSV: " << csv_abs
              << " (size=" << size_before << ")\n";

    CsvDoc csv;
    if (!csv.load(csv_path)) {
        std::cerr << "[radcorr] ERROR: failed to load CSV.\n";
        return false;
    }

    std::vector<RowBin> bins = build_row_bins(csv);

    // Define groups: 10.6 GeV (4 periods) and 10.2 GeV (Sp19 Inb only)
    RCGroup g10p6;
    g10p6.label      = "10.6 GeV";
    g10p6.xbavg_col  = "xBavg, 10.6 GeV";
    g10p6.phiavg_col = "phiavg, 10.6 GeV";
    g10p6.q2avg_col  = "Q2avg, 10.6 GeV";
    g10p6.tabavg_col = "t_abs_avg, 10.6 GeV";
    g10p6.frad_col   = "Frad, 10.6 GeV";
    g10p6.energy_dir = "10.60";
    g10p6.mc_pairs.push_back(McPair{"DVCS_Sp18_inb_gen",  "DVCS_Sp18_inb_gen_rad",  "Sp18 Inb"});
    g10p6.mc_pairs.push_back(McPair{"DVCS_Sp18_out_gen",  "DVCS_Sp18_out_gen_rad",  "Sp18 Out"});

    if (SKIP_FA18_INB_10P6) {
        std::cout << "[radcorr] NOTE: SKIP_FA18_INB_10P6=true, omitting Fa18 Inb "
                  << "from 10.6 GeV Frad (using Sp18 Inb, Sp18 Out, Fa18 Out only).\n";
    } else {
        g10p6.mc_pairs.push_back(McPair{"DVCS_Fa18_inb_gen",  "DVCS_Fa18_inb_gen_rad",  "Fa18 Inb"});
    }

    g10p6.mc_pairs.push_back(McPair{"DVCS_Fa18_out_gen",  "DVCS_Fa18_out_gen_rad",  "Fa18 Out"});

    RCGroup g10p2;
    g10p2.label      = "10.2 GeV";
    g10p2.xbavg_col  = "xBavg, Sp19 Inb";
    g10p2.phiavg_col = "phiavg, Sp19 Inb";
    g10p2.q2avg_col  = "Q2avg, Sp19 Inb";
    g10p2.tabavg_col = "t_abs_avg, Sp19 Inb";
    g10p2.frad_col   = "Frad, 10.2 GeV";
    g10p2.energy_dir = "10.2";
    g10p2.mc_pairs.push_back(McPair{"DVCS_Sp19_inb_gen",  "DVCS_Sp19_inb_gen_rad",  "Sp19 Inb"});

    std::vector<RCGroup> groups;
    groups.push_back(g10p6);
    groups.push_back(g10p2);

    // Pre-check: ensure all required trees exist before launching threads.
    for (std::size_t ig = 0; ig < groups.size(); ++ig) {
        const RCGroup& G = groups[ig];
        for (std::size_t ip = 0; ip < G.mc_pairs.size(); ++ip) {
            const McPair& P = G.mc_pairs[ip];
            if (genMcTrees.find(P.bornTag) == genMcTrees.end()) {
                std::cerr << "[radcorr] FATAL: missing Born MC tree with tag '"
                          << P.bornTag << "' for group " << G.label
                          << " (period " << P.period_label << ").\n";
                return false;
            }
            if (radGenMcTrees.find(P.radTag) == radGenMcTrees.end()) {
                std::cerr << "[radcorr] FATAL: missing Rad MC tree with tag '"
                          << P.radTag << "' for group " << G.label
                          << " (period " << P.period_label << ").\n";
                return false;
            }
        }
    }

    // Enable ROOT thread safety for concurrent TTree access.
    ROOT::EnableThreadSafety();

    // Build row_has_data for each group (decides where we compute Frad).
    std::vector<bool> row_has_data_10p6 = build_row_has_data(csv, g10p6.xbavg_col);
    std::vector<bool> row_has_data_10p2 = build_row_has_data(csv, g10p2.xbavg_col);

    const int NR = csv.nrows();
    (void)NR;

    GroupAccum acc10p6;
    GroupAccum acc10p2;

    // Launch 4 worker threads:
    //   - 10.6 GeV Born
    //   - 10.6 GeV Rad
    //   - 10.2 GeV Born
    //   - 10.2 GeV Rad

    std::thread t_10p6_born(
        accumulate_group_type,
        std::cref(g10p6),
        std::cref(bins),
        std::cref(row_has_data_10p6),
        std::cref(genMcTrees),
        true,
        std::ref(acc10p6.born_counts),
        std::ref(acc10p6.N_born)
    );

    std::thread t_10p6_rad(
        accumulate_group_type,
        std::cref(g10p6),
        std::cref(bins),
        std::cref(row_has_data_10p6),
        std::cref(radGenMcTrees),
        false,
        std::ref(acc10p6.rad_counts),
        std::ref(acc10p6.N_rad)
    );

    std::thread t_10p2_born(
        accumulate_group_type,
        std::cref(g10p2),
        std::cref(bins),
        std::cref(row_has_data_10p2),
        std::cref(genMcTrees),
        true,
        std::ref(acc10p2.born_counts),
        std::ref(acc10p2.N_born)
    );

    std::thread t_10p2_rad(
        accumulate_group_type,
        std::cref(g10p2),
        std::cref(bins),
        std::cref(row_has_data_10p2),
        std::cref(radGenMcTrees),
        false,
        std::ref(acc10p2.rad_counts),
        std::ref(acc10p2.N_rad)
    );

    // Wait for all workers to finish before touching CSV or producing plots.
    t_10p6_born.join();
    t_10p6_rad.join();
    t_10p2_born.join();
    t_10p2_rad.join();

    // Now fill Frad columns and draw plots in the main thread.

    if (!fill_frad_for_group(g10p6.label,
                             g10p6.frad_col,
                             csv,
                             row_has_data_10p6,
                             acc10p6.born_counts,
                             acc10p6.rad_counts,
                             acc10p6.N_born,
                             acc10p6.N_rad)) {
        std::cerr << "[radcorr] ERROR: fill_frad_for_group failed for "
                  << g10p6.label << ".\n";
        return false;
    }

    if (!fill_frad_for_group(g10p2.label,
                             g10p2.frad_col,
                             csv,
                             row_has_data_10p2,
                             acc10p2.born_counts,
                             acc10p2.rad_counts,
                             acc10p2.N_born,
                             acc10p2.N_rad)) {
        std::cerr << "[radcorr] ERROR: fill_frad_for_group failed for "
                  << g10p2.label << ".\n";
        return false;
    }

    // Draw canvases (sequential, no threads).
    draw_radiative_canvases(g10p6.label,
                            g10p6.xbavg_col,
                            g10p6.phiavg_col,
                            g10p6.q2avg_col,
                            g10p6.tabavg_col,
                            g10p6.frad_col,
                            csv,
                            row_has_data_10p6,
                            out_root_dir,
                            g10p6.energy_dir);

    draw_radiative_canvases(g10p2.label,
                            g10p2.xbavg_col,
                            g10p2.phiavg_col,
                            g10p2.q2avg_col,
                            g10p2.tabavg_col,
                            g10p2.frad_col,
                            csv,
                            row_has_data_10p2,
                            out_root_dir,
                            g10p2.energy_dir);

    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[radcorr] ERROR: failed to save updated CSV.\n";
        return false;
    }

    const uintmax_t size_after =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[radcorr] Updated CSV: " << csv_abs
              << " (size " << size_before << " -> " << size_after << ")\n";
    std::cout << "[radcorr] Radiative corrections (Frad) complete.\n";

    return true;
}