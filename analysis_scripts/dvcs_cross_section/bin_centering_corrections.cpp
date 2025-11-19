#include "bin_centering_corrections.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// ROOT plotting includes for debug plots
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TGaxis.h>
#include <TLine.h>

namespace {

namespace fs = std::filesystem;

// ------------------------- CSV helper -------------------------

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
        bool needq = (s.find(',') != std::string::npos) ||
                     (s.find('"') != std::string::npos);
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
        char* e = 0;
        double v = std::strtod(s.c_str(), &e);
        if (e == s.c_str()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return v;
    }

    bool load(const std::string& path) {
        std::cout << "[bincenter] Loading CSV: " << path << "\n";
        std::ifstream fin(path.c_str());
        if (!fin.is_open()) {
            std::cerr << "[bincenter] ERROR: cannot open CSV: " << path << "\n";
            return false;
        }
        std::string line;
        if (!std::getline(fin, line)) {
            std::cerr << "[bincenter] ERROR: empty CSV: " << path << "\n";
            return false;
        }
        header = split_csv_line(line);
        index.clear();
        for (int i = 0; i < (int)header.size(); ++i) {
            index[header[i]] = i;
        }
        std::cout << "[bincenter] Header has " << header.size() << " columns.\n";
        rows.clear();
        int n_read = 0;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            rows.push_back(split_csv_line(line));
            ++n_read;
        }
        for (std::size_t r = 0; r < rows.size(); ++r) {
            rows[r].resize(header.size());
        }
        std::cout << "[bincenter] Loaded " << n_read << " data rows.\n";
        return true;
    }

    bool save_atomic(const std::string& path) const {
        std::cout << "[bincenter] Saving updated CSV (atomic write) to: " << path << "\n";
        const std::string tmp = path + ".tmp";
        {
            std::ofstream fout(tmp.c_str());
            if (!fout.is_open()) {
                std::cerr << "[bincenter] ERROR: cannot write CSV tmp: " << tmp << "\n";
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
        fs::rename(tmp, path, ec);
        if (ec) {
            std::cerr << "[bincenter] First rename attempt failed (" << ec.message()
                      << "), trying to remove target and retry.\n";
            std::remove(path.c_str());
            fs::rename(tmp, path, ec);
            if (ec) {
                std::cerr << "[bincenter] ERROR: atomic rename failed ("
                          << ec.message() << ")\n";
                return false;
            }
        }
        std::cout << "[bincenter] CSV save complete.\n";
        return true;
    }

    int nrows() const {
        return (int)rows.size();
    }

    int col_index(const std::string& name) const {
        std::unordered_map<std::string,int>::const_iterator it = index.find(name);
        if (it == index.end()) return -1;
        return it->second;
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

// "(value, stat, sys)" triple formatter.
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

// Parse "(value, stat, sys)" into three doubles.
// Returns true on success.
static bool parse_triple(const std::string& s, double& v, double& s_stat, double& s_sys) {
    std::string t = s;
    // Trim
    while (!t.empty() && std::isspace(static_cast<unsigned char>(t.front()))) t.erase(t.begin());
    while (!t.empty() && std::isspace(static_cast<unsigned char>(t.back())))  t.pop_back();
    if (t.size() < 5 || t.front() != '(' || t.back() != ')') {
        return false;
    }
    t = t.substr(1, t.size() - 2); // remove '(' and ')'

    std::vector<std::string> parts;
    std::string cur;
    for (std::size_t i = 0; i < t.size(); ++i) {
        char c = t[i];
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

    auto to_d = [](const std::string& s_in, double& out) -> bool {
        std::string u = s_in;
        while (!u.empty() && std::isspace(static_cast<unsigned char>(u.front()))) u.erase(u.begin());
        while (!u.empty() && std::isspace(static_cast<unsigned char>(u.back())))  u.pop_back();
        if (u.empty()) return false;
        char* e = 0;
        double v = std::strtod(u.c_str(), &e);
        if (e == u.c_str()) return false;
        out = v;
        return true;
    };

    if (!to_d(parts[0], v))      return false;
    if (!to_d(parts[1], s_stat)) return false;
    if (!to_d(parts[2], s_sys))  return false;
    return true;
}

static inline bool finite_pos(double x) {
    return std::isfinite(x) && x > 0.0;
}

// Flags for which rows are actually populated for a given group,
// based on a specific xBavg column.
static std::vector<bool>
build_row_has_data(const CsvDoc& csv, const std::string& xbavg_col_name) {
    const int NR = csv.nrows();
    int c = csv.col_index(xbavg_col_name);
    if (c < 0) {
        std::cerr << "[bincenter] FATAL: missing column '" << xbavg_col_name
                  << "' needed to decide where to compute Fbin.\n";
        std::exit(EXIT_FAILURE);
    }

    std::vector<bool> flags(NR, false);
    int count_true = 0;
    for (int r = 0; r < NR; ++r) {
        const std::string& cell = csv.rows[r][c];
        if (cell.empty()) {
            flags[r] = false;
        } else {
            double v = CsvDoc::to_double(cell);
            flags[r] = std::isfinite(v);
        }
        if (flags[r]) ++count_true;
    }
    std::cout << "[bincenter] Rows with non-empty '" << xbavg_col_name
              << "': " << count_true << " / " << NR << "\n";
    return flags;
}

// Helper struct for result storage per row.
struct FbinResult {
    double value;
    double stat;
    double sys;
    bool   valid;
    FbinResult() : value(1.0), stat(0.0), sys(0.0), valid(false) {}
};

// Compute Fbin for a single row and energy group, given:
//
//   reference point (bin averages):
//      xB_center, Q2_center, tpos_center, phi_center_deg
//
//   bin edges:
//      xBmin,xBmax, Q2min,Q2max, t_abs_min,t_abs_max, phimin,phimax
//
//   energy, model setup, and sub-grid size.
static FbinResult compute_fbin_for_row(
    double xB_center,
    double Q2_center,
    double tpos_center,
    double phi_center_deg,
    double xBmin,
    double xBmax,
    double Q2min,
    double Q2max,
    double t_abs_min,
    double t_abs_max,
    double phimin,
    double phimax,
    double Ebeam,
    int n_steps,
    const ModelPaths& paths,
    bool vgg_globalfit,
    ModelChoice model_choice)
{
    FbinResult res;

    // Basic sanity checks on bin widths.
    const double dx = xBmax   - xBmin;
    const double dQ = Q2max   - Q2min;
    const double dt = t_abs_max - t_abs_min;
    const double dP = phimax  - phimin;

    if (!(dx > 0.0 && dQ > 0.0 && dt > 0.0 && dP > 0.0)) {
        return res;
    }

    if (!(finite_pos(xB_center) &&
          finite_pos(Q2_center) &&
          finite_pos(tpos_center))) {
        return res;
    }
    if (!std::isfinite(phi_center_deg)) {
        return res;
    }

    // Sub-binning size (per dimension). Keep at least 2.
    const int nx = std::max(2, n_steps);
    const int nQ = std::max(2, n_steps);
    const int nt = std::max(2, n_steps);
    const int np = std::max(2, n_steps);

    const double dx_step = dx / double(nx);
    const double dQ_step = dQ / double(nQ);
    const double dt_step = dt / double(nt);
    const double dP_step = dP / double(np);

    // Center-point model values.
    double center_km15 = 0.0;
    double center_vgg  = 0.0;

    const bool want_km15 = (model_choice == ModelChoice::Both ||
                            model_choice == ModelChoice::KM15Only);
    const bool want_vgg  = (model_choice == ModelChoice::Both ||
                            model_choice == ModelChoice::VGGOnly);

    if (want_km15) {
        center_km15 = km15_xs(
            xB_center, Q2_center, tpos_center, phi_center_deg,
            Ebeam, Helicity::Unpol, paths);
    }
    if (want_vgg) {
        center_vgg = vgg_xs(
            xB_center, Q2_center, tpos_center, phi_center_deg,
            Ebeam, Helicity::Unpol, paths, vgg_globalfit);
    }

    const bool need_km15 = want_km15 && finite_pos(center_km15);
    const bool need_vgg  = want_vgg  && finite_pos(center_vgg);

    if (!need_km15 && !need_vgg) {
        return res;
    }

    // Integrate model over the bin using midpoints in each dimension.
    double sum_km15 = 0.0;
    double sum_vgg  = 0.0;
    int    cnt_km15 = 0;
    int    cnt_vgg  = 0;

    for (int ix = 0; ix < nx; ++ix) {
        const double xB = xBmin + (ix + 0.5) * dx_step;
        for (int iQ = 0; iQ < nQ; ++iQ) {
            const double Q2 = Q2min + (iQ + 0.5) * dQ_step;
            for (int it = 0; it < nt; ++it) {
                const double tpos = t_abs_min + (it + 0.5) * dt_step;
                for (int ip = 0; ip < np; ++ip) {
                    const double phi_deg = phimin + (ip + 0.5) * dP_step;

                    if (need_km15) {
                        const double vk = km15_xs(
                            xB, Q2, tpos, phi_deg, Ebeam,
                            Helicity::Unpol, paths);
                        if (finite_pos(vk)) {
                            sum_km15 += vk;
                            ++cnt_km15;
                        }
                    }

                    if (need_vgg) {
                        const double vv = vgg_xs(
                            xB, Q2, tpos, phi_deg, Ebeam,
                            Helicity::Unpol, paths, vgg_globalfit);
                        if (finite_pos(vv)) {
                            sum_vgg += vv;
                            ++cnt_vgg;
                        }
                    }
                }
            }
        }
    }

    double f_km15 = 0.0;
    double f_vgg  = 0.0;

    if (need_km15 && cnt_km15 > 0) {
        const double avg_km15 = sum_km15 / double(cnt_km15);
        if (finite_pos(avg_km15)) {
            f_km15 = center_km15 / avg_km15;
        }
    }
    if (need_vgg && cnt_vgg > 0) {
        const double avg_vgg = sum_vgg / double(cnt_vgg);
        if (finite_pos(avg_vgg)) {
            f_vgg = center_vgg / avg_vgg;
        }
    }

    double f_final = 1.0;
    double s_sys   = 0.0;

    switch (model_choice) {
        case ModelChoice::Both:
            if (finite_pos(f_km15) && finite_pos(f_vgg)) {
                f_final = 0.5 * (f_km15 + f_vgg);
                // Systematic = absolute difference between the two models.
                s_sys   = std::fabs(f_km15 - f_vgg);
            } else if (finite_pos(f_km15)) {
                f_final = f_km15;
                s_sys   = 1.0; // flag: only one model available
            } else if (finite_pos(f_vgg)) {
                f_final = f_vgg;
                s_sys   = 1.0; // flag: only one model available
            } else {
                return res;
            }
            break;

        case ModelChoice::VGGOnly:
            if (!finite_pos(f_vgg)) {
                return res;
            }
            f_final = f_vgg;
            s_sys   = 0.0;
            break;

        case ModelChoice::KM15Only:
            if (!finite_pos(f_km15)) {
                return res;
            }
            f_final = f_km15;
            s_sys   = 0.0;
            break;
    }

    if (!std::isfinite(f_final) || f_final <= 0.0) {
        return res;
    }

    res.value = f_final;
    res.stat  = 0.0;    // no statistical uncertainty here; pure model factor
    res.sys   = s_sys;
    res.valid = true;
    return res;
}

// Utility to gather unique ranges (min,max) from CSV columns.
static std::vector<std::pair<double,double> >
gather_unique_ranges(const CsvDoc& csv, int c_min, int c_max) {
    std::set<std::pair<double,double> > s;
    const int NR = csv.nrows();
    for (int r = 0; r < NR; ++r) {
        double lo = csv.as_double(r, c_min);
        double hi = csv.as_double(r, c_max);
        if (!std::isfinite(lo) || !std::isfinite(hi)) continue;
        if (hi <= lo) continue;
        s.insert(std::make_pair(lo, hi));
    }
    return std::vector<std::pair<double,double> >(s.begin(), s.end());
}

// Simple helper for degree tick marks under the main frame axis.
static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

} // anonymous namespace

// -------------------------------------------------------------
// Public: compute Fbin and fill CSV
// -------------------------------------------------------------
bool update_bin_centering_corrections_csv(
    const std::string& csv_path,
    int n_steps,
    const ModelPaths& paths,
    bool vgg_globalfit,
    ModelChoice model_choice)
{
    std::cout << "============================================================\n";
    std::cout << "[bincenter] Starting bin-centering corrections.\n";
    std::cout << "[bincenter] Input CSV: " << csv_path << "\n";
    std::cout << "[bincenter] n_steps (per dimension) = " << n_steps
              << " (total model calls per row per model ~ n_steps^4)\n";
    std::cout << "[bincenter] Model choice = "
              << (model_choice == ModelChoice::Both ? "Both"
                  : (model_choice == ModelChoice::VGGOnly ? "VGGOnly" : "KM15Only"))
              << ", vgg_globalfit = " << (vgg_globalfit ? "true" : "false") << "\n";

    const std::string csv_abs = fs::absolute(csv_path).string();
    std::error_code ec;
    const uintmax_t size_before =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[bincenter] CSV absolute path: " << csv_abs
              << " (size = " << size_before << " bytes)\n";

    if (n_steps < 2) {
        std::cerr << "[bincenter] FATAL: n_steps must be >= 2 (got "
                  << n_steps << ")\n";
        return false;
    }

    CsvDoc csv;
    if (!csv.load(csv_path)) {
        std::cerr << "[bincenter] ERROR: failed to load CSV.\n";
        return false;
    }

    const int NR = csv.nrows();
    if (NR <= 0) {
        std::cerr << "[bincenter] FATAL: CSV has no data rows.\n";
        return false;
    }
    std::cout << "[bincenter] Total rows to inspect: " << NR << "\n";

    // Bin-edge columns (must exist).
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
        std::cerr << "[bincenter] FATAL: missing one or more bin-edge columns "
                  << "(xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax)\n";
        return false;
    }

    // Group-specific average columns and Fbin columns.

    // 10.6 GeV group (combined periods).
    const std::string col_xbavg_10p6  = "xBavg, 10.6 GeV";
    const std::string col_q2avg_10p6  = "Q2avg, 10.6 GeV";
    const std::string col_tabavg_10p6 = "t_abs_avg, 10.6 GeV";
    const std::string col_phiavg_10p6 = "phiavg, 10.6 GeV";
    const std::string col_fbin_10p6   = "Fbin, 10.6 GeV";

    const int c_xbavg_10p6  = csv.col_index(col_xbavg_10p6);
    const int c_q2avg_10p6  = csv.col_index(col_q2avg_10p6);
    const int c_tabavg_10p6 = csv.col_index(col_tabavg_10p6);
    const int c_phiavg_10p6 = csv.col_index(col_phiavg_10p6);
    const int c_fbin_10p6   = csv.col_index(col_fbin_10p6);

    if (c_xbavg_10p6 < 0 || c_q2avg_10p6 < 0 ||
        c_tabavg_10p6 < 0 || c_phiavg_10p6 < 0 ||
        c_fbin_10p6   < 0) {
        std::cerr << "[bincenter] FATAL: missing one or more columns for 10.6 GeV group:\n"
                  << "  required: \"" << col_xbavg_10p6  << "\", "
                  << "\"" << col_q2avg_10p6  << "\", "
                  << "\"" << col_tabavg_10p6 << "\", "
                  << "\"" << col_phiavg_10p6 << "\", "
                  << "\"" << col_fbin_10p6   << "\"\n";
        return false;
    }

    // 10.2 GeV group (Sp19 Inb only).
    const std::string col_xbavg_10p2  = "xBavg, Sp19 Inb";
    const std::string col_q2avg_10p2  = "Q2avg, Sp19 Inb";
    const std::string col_tabavg_10p2 = "t_abs_avg, Sp19 Inb";
    const std::string col_phiavg_10p2 = "phiavg, Sp19 Inb";
    const std::string col_fbin_10p2   = "Fbin, 10.2 GeV";

    const int c_xbavg_10p2  = csv.col_index(col_xbavg_10p2);
    const int c_q2avg_10p2  = csv.col_index(col_q2avg_10p2);
    const int c_tabavg_10p2 = csv.col_index(col_tabavg_10p2);
    const int c_phiavg_10p2 = csv.col_index(col_phiavg_10p2);
    const int c_fbin_10p2   = csv.col_index(col_fbin_10p2);

    if (c_xbavg_10p2 < 0 || c_q2avg_10p2 < 0 ||
        c_tabavg_10p2 < 0 || c_phiavg_10p2 < 0 ||
        c_fbin_10p2   < 0) {
        std::cerr << "[bincenter] FATAL: missing one or more columns for 10.2 GeV group:\n"
                  << "  required: \"" << col_xbavg_10p2  << "\", "
                  << "\"" << col_q2avg_10p2  << "\", "
                  << "\"" << col_tabavg_10p2 << "\", "
                  << "\"" << col_phiavg_10p2 << "\", "
                  << "\"" << col_fbin_10p2   << "\"\n";
        return false;
    }

    std::cout << "[bincenter] All required columns found. Building row masks...\n";

    // Determine which rows actually carry data for each group,
    // based on whether xBavg column is populated.
    std::vector<bool> row_has_data_10p6 = build_row_has_data(csv, col_xbavg_10p6);
    std::vector<bool> row_has_data_10p2 = build_row_has_data(csv, col_xbavg_10p2);

    // Storage for results.
    std::vector<FbinResult> res_10p6(NR);
    std::vector<FbinResult> res_10p2(NR);

#ifdef _OPENMP
    int hard_cap = 5;
    int want     = omp_get_max_threads();
    int use      = std::min(hard_cap, std::max(1, want));
    omp_set_num_threads(use);
    std::cout << "[bincenter] OpenMP enabled with " << use
              << " worker(s) (hard-capped at 5)\n";
#else
    std::cout << "[bincenter] OpenMP not available; running single-threaded\n";
#endif

    const double Ebeam_10p6 = 10.6;
    const double Ebeam_10p2 = 10.2;

    std::cout << "[bincenter] Beam energies: 10.6 GeV and 10.2 GeV.\n";
    std::cout << "[bincenter] Beginning Fbin computation over all rows...\n";

    std::atomic<int> processed_rows(0);

    // Parallel loop over rows. For each row, we may compute 10.6, 10.2, both, or neither.
#pragma omp parallel for schedule(dynamic)
    for (int r = 0; r < NR; ++r) {
        // Bin edges (shared across groups).
        const double xbmin  = csv.as_double(r, c_xb_min);
        const double xbmax  = csv.as_double(r, c_xb_max);
        const double q2min  = csv.as_double(r, c_q2_min);
        const double q2max  = csv.as_double(r, c_q2_max);
        const double tmin   = csv.as_double(r, c_tab_min);
        const double tmax   = csv.as_double(r, c_tab_max);
        const double phimin = csv.as_double(r, c_phi_min);
        const double phimax = csv.as_double(r, c_phi_max);

        // 10.6 GeV group.
        if (row_has_data_10p6[r]) {
            const double xb_c   = csv.as_double(r, c_xbavg_10p6);
            const double q2_c   = csv.as_double(r, c_q2avg_10p6);
            const double tpos_c = csv.as_double(r, c_tabavg_10p6);
            const double phi_c  = csv.as_double(r, c_phiavg_10p6);

            FbinResult fr = compute_fbin_for_row(
                xb_c, q2_c, tpos_c, phi_c,
                xbmin, xbmax,
                q2min, q2max,
                tmin, tmax,
                phimin, phimax,
                Ebeam_10p6,
                n_steps,
                paths,
                vgg_globalfit,
                model_choice);

            res_10p6[r] = fr;
        }

        // 10.2 GeV group (Sp19 Inb).
        if (row_has_data_10p2[r]) {
            const double xb_c   = csv.as_double(r, c_xbavg_10p2);
            const double q2_c   = csv.as_double(r, c_q2avg_10p2);
            const double tpos_c = csv.as_double(r, c_tabavg_10p2);
            const double phi_c  = csv.as_double(r, c_phiavg_10p2);

            FbinResult fr = compute_fbin_for_row(
                xb_c, q2_c, tpos_c, phi_c,
                xbmin, xbmax,
                q2min, q2max,
                tmin, tmax,
                phimin, phimax,
                Ebeam_10p2,
                n_steps,
                paths,
                vgg_globalfit,
                model_choice);

            res_10p2[r] = fr;
        }

        int done = processed_rows.fetch_add(1, std::memory_order_relaxed) + 1;
        if (done % 50 == 0 || done == NR) {
#pragma omp critical
            {
                std::cout << "[bincenter] Progress: " << done << " / " << NR
                          << " rows processed.\n";
            }
        }
    } // end parallel for

    std::cout << "[bincenter] Finished Fbin model evaluations for all rows.\n";
    std::cout << "[bincenter] Writing results back into CSV columns...\n";

    // Fill CSV with Fbin triples.
    std::size_t write_10p6 = 0;
    std::size_t write_10p2 = 0;

    for (int r = 0; r < NR; ++r) {
        if (row_has_data_10p6[r] && res_10p6[r].valid) {
            csv.rows[r][c_fbin_10p6] =
                format_triple(res_10p6[r].value, res_10p6[r].stat, res_10p6[r].sys);
            ++write_10p6;
        }
        if (row_has_data_10p2[r] && res_10p2[r].valid) {
            csv.rows[r][c_fbin_10p2] =
                format_triple(res_10p2[r].value, res_10p2[r].stat, res_10p2[r].sys);
            ++write_10p2;
        }
    }

    std::cout << "[bincenter] Rows with valid Fbin, 10.6 GeV: " << write_10p6 << "\n";
    std::cout << "[bincenter] Rows with valid Fbin, 10.2 GeV: " << write_10p2 << "\n";

    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[bincenter] ERROR: failed to save updated CSV.\n";
        return false;
    }

    const uintmax_t size_after =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[bincenter] Updated CSV: " << csv_abs
              << " (size " << size_before << " -> " << size_after << " bytes)\n";
    std::cout << "[bincenter] Bin-centering corrections complete.\n";
    std::cout << "============================================================\n";

    return true;
}

// -------------------------------------------------------------
// Debug plots: Fbin vs phi
// -------------------------------------------------------------
void plot_bin_centering_fbin_vs_phi(
    const std::string& csv_path,
    const std::string& out_root_dir)
{
    std::cout << "============================================================\n";
    std::cout << "[bincenter-plot] Starting Fbin vs phi debug plot generation.\n";
    std::cout << "[bincenter-plot] CSV: " << csv_path << "\n";
    std::cout << "[bincenter-plot] Output root directory: " << out_root_dir << "\n";

    CsvDoc csv;
    if (!csv.load(csv_path)) {
        std::cerr << "[bincenter-plot] ERROR: failed to load CSV: "
                  << csv_path << "\n";
        return;
    }

    const int NR = csv.nrows();
    if (NR <= 0) {
        std::cerr << "[bincenter-plot] CSV has no data rows.\n";
        return;
    }

    // Columns we will use.
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
        std::cerr << "[bincenter-plot] FATAL: missing bin-edge columns.\n";
        return;
    }

    // Columns for 10.6 GeV.
    const std::string col_xbavg_10p6  = "xBavg, 10.6 GeV";
    const std::string col_phiavg_10p6 = "phiavg, 10.6 GeV";
    const std::string col_fbin_10p6   = "Fbin, 10.6 GeV";

    const int c_xbavg_10p6  = csv.col_index(col_xbavg_10p6);
    const int c_phiavg_10p6 = csv.col_index(col_phiavg_10p6);
    const int c_fbin_10p6   = csv.col_index(col_fbin_10p6);

    if (c_xbavg_10p6 < 0 || c_phiavg_10p6 < 0 || c_fbin_10p6 < 0) {
        std::cerr << "[bincenter-plot] WARNING: missing some 10.6 GeV columns; "
                  << "10.6 Fbin vs phi plots will be skipped.\n";
    }

    // Columns for 10.2 GeV (Sp19 Inb).
    const std::string col_xbavg_10p2  = "xBavg, Sp19 Inb";
    const std::string col_phiavg_10p2 = "phiavg, Sp19 Inb";
    const std::string col_fbin_10p2   = "Fbin, 10.2 GeV";

    const int c_xbavg_10p2  = csv.col_index(col_xbavg_10p2);
    const int c_phiavg_10p2 = csv.col_index(col_phiavg_10p2);
    const int c_fbin_10p2   = csv.col_index(col_fbin_10p2);

    if (c_xbavg_10p2 < 0 || c_phiavg_10p2 < 0 || c_fbin_10p2 < 0) {
        std::cerr << "[bincenter-plot] WARNING: missing some 10.2 GeV columns; "
                  << "10.2 Fbin vs phi plots will be skipped.\n";
    }

    // Build dir structure: out_root_dir/10.60 and out_root_dir/10.2.
    const fs::path base(out_root_dir);
    const fs::path dir_10p6 = base / "10.60";
    const fs::path dir_10p2 = base / "10.2";

    std::error_code ec;
    fs::create_directories(dir_10p6, ec);
    fs::create_directories(dir_10p2, ec);

    if (ec) {
        std::cerr << "[bincenter-plot] WARNING: error while creating directories: "
                  << ec.message() << "\n";
    }

    // Unique ranges for xB, Q2, t.
    std::vector<std::pair<double,double> > xb_ranges  = gather_unique_ranges(csv, c_xb_min,  c_xb_max);
    std::vector<std::pair<double,double> > q2_ranges  = gather_unique_ranges(csv, c_q2_min,  c_q2_max);
    std::vector<std::pair<double,double> > t_ranges   = gather_unique_ranges(csv, c_tab_min, c_tab_max);

    std::cout << "[bincenter-plot] Found "
              << xb_ranges.size() << " xB ranges, "
              << q2_ranges.size() << " Q2 ranges, "
              << t_ranges.size()  << " t ranges.\n";

    if (xb_ranges.empty() || q2_ranges.empty() || t_ranges.empty()) {
        std::cerr << "[bincenter-plot] No non-empty (xB,Q2,t) ranges found.\n";
        return;
    }

    // Global ROOT style tweaks (no stats, tick marks on both axes).
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);

    auto build_row_has_data_local = [&](const std::string& col_name)->std::vector<bool> {
        if (csv.col_index(col_name) < 0) {
            std::cout << "[bincenter-plot] Column '" << col_name
                      << "' not found; treating all rows as empty for this group.\n";
            return std::vector<bool>(NR, false);
        }
        return ::build_row_has_data(csv, col_name);
    };

    std::vector<bool> row_has_data_10p6 = build_row_has_data_local(col_xbavg_10p6);
    std::vector<bool> row_has_data_10p2 = build_row_has_data_local(col_xbavg_10p2);

    // Make a helper lambda that, for a given energy group, makes canvases
    // for all xB ranges and saves them under a given directory.
    auto make_plots_for_energy =
        [&](const std::string& energy_label,
            const fs::path& out_dir,
            int c_xbavg,
            int c_phiavg,
            int c_fbin,
            const std::vector<bool>& row_has_data)
    {
        if (c_xbavg < 0 || c_phiavg < 0 || c_fbin < 0) {
            std::cout << "[bincenter-plot] Skipping energy " << energy_label
                      << " GeV: required columns not present.\n";
            return;
        }

        std::cout << "[bincenter-plot] ------------------------------------------------------------\n";
        std::cout << "[bincenter-plot] Building Fbin vs phi plots for energy " << energy_label
                  << " GeV.\n";

        int canvases_made = 0;
        int pads_with_data = 0;

        // For each xB range, we build a canvas of Q2 (columns) vs t (rows).
        for (std::size_t ix = 0; ix < xb_ranges.size(); ++ix) {
            const std::pair<double,double>& xb = xb_ranges[ix];

            // Build slices of Q2 and t that actually have data for this xB
            // and this energy group.
            std::set<std::pair<double,double> > q2_slice_set;
            std::set<std::pair<double,double> > t_slice_set;

            for (int r = 0; r < NR; ++r) {
                if (!row_has_data[r]) continue;
                double xbmin = csv.as_double(r, c_xb_min);
                double xbmax = csv.as_double(r, c_xb_max);
                if (!std::isfinite(xbmin) || !std::isfinite(xbmax)) continue;
                if (xbmin != xb.first || xbmax != xb.second) continue;

                double q2min = csv.as_double(r, c_q2_min);
                double q2max = csv.as_double(r, c_q2_max);
                double tmin  = csv.as_double(r, c_tab_min);
                double tmax  = csv.as_double(r, c_tab_max);

                if (std::isfinite(q2min) && std::isfinite(q2max) && q2max > q2min) {
                    q2_slice_set.insert(std::make_pair(q2min, q2max));
                }
                if (std::isfinite(tmin) && std::isfinite(tmax) && tmax > tmin) {
                    t_slice_set.insert(std::make_pair(tmin, tmax));
                }
            }

            if (q2_slice_set.empty() || t_slice_set.empty()) {
                // No data at this xB for this energy group.
                continue;
            }

            std::vector<std::pair<double,double> > q2_slice(
                q2_slice_set.begin(), q2_slice_set.end());
            std::vector<std::pair<double,double> > t_slice(
                t_slice_set.begin(), t_slice_set.end());

            const int ncols = (int)q2_slice.size();
            const int nrows = (int)t_slice.size();

            std::cout << "[bincenter-plot] Energy " << energy_label
                      << " GeV, xB bin index " << ix
                      << " with xB in [" << xb.first << ", " << xb.second
                      << "]: Q2 bins = " << ncols
                      << ", t bins = " << nrows << ".\n";

            // Determine y-range over all Fbin values in this xB bin.
            double fmin = std::numeric_limits<double>::infinity();
            double fmax = 0.0;

            for (int r = 0; r < NR; ++r) {
                if (!row_has_data[r]) continue;

                double xbmin = csv.as_double(r, c_xb_min);
                double xbmax = csv.as_double(r, c_xb_max);
                if (!std::isfinite(xbmin) || !std::isfinite(xbmax)) continue;
                if (xbmin != xb.first || xbmax != xb.second) continue;

                const std::string& fstr = csv.rows[r][c_fbin];
                if (fstr.empty()) continue;

                double val, s_stat, s_sys;
                if (!parse_triple(fstr, val, s_stat, s_sys)) continue;
                if (!std::isfinite(val) || val <= 0.0) continue;

                fmin = std::min(fmin, val);
                fmax = std::max(fmax, val);
            }

            if (!std::isfinite(fmin) || fmax <= 0.0) {
                std::cout << "[bincenter-plot]   No valid Fbin values in this xB bin; skipping canvas.\n";
                continue;
            }

            // Add a bit of padding; clamp to a sensible range around 1.
            double ymin = std::max(0.5, 0.95 * fmin);
            double ymax = std::min(2.0, 1.05 * fmax);
            if (ymax <= ymin) {
                ymin = 0.8;
                ymax = 1.2;
            }

            std::cout << "[bincenter-plot]   Global Fbin range in this xB bin: ["
                      << fmin << ", " << fmax << "], frame y-range = ["
                      << ymin << ", " << ymax << "].\n";

            // Canvas dimensions: consistent with your other modules.
            const int W = 300 * ncols + 160;
            const int H = 260 * nrows + 240;

            std::ostringstream cname;
            cname << "c_Fbin_phi_" << energy_label << "_xB" << ix;
            TCanvas* c = new TCanvas(cname.str().c_str(),
                                     cname.str().c_str(), W, H);

            TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
            pTop->SetFillStyle(0);
            pTop->SetBorderSize(0);
            pTop->Draw();

            TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.90);
            pGrid->SetFillStyle(0);
            pGrid->SetBorderSize(0);
            pGrid->Draw();
            pGrid->cd();
            pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

            // Top title.
            pTop->cd();
            TLatex head;
            head.SetNDC();
            head.SetTextAlign(22);
            head.SetTextFont(42);
            head.SetTextSize(0.16);
            std::ostringstream tit;
            tit << "F_{bin} vs #phi   " << energy_label
                << " GeV   x_{B} #in ["
                << std::setprecision(3) << xb.first << ", "
                << std::setprecision(3) << xb.second << "]";
            head.DrawLatex(0.5, 0.55, tit.str().c_str());

            // Now fill each pad.
            for (int ir = 0; ir < nrows; ++ir) {
                const std::pair<double,double>& tR = t_slice[ir];
                for (int ic = 0; ic < ncols; ++ic) {
                    const std::pair<double,double>& qR = q2_slice[ic];

                    pGrid->cd(ir * ncols + ic + 1);
                    gPad->SetGrid(1,1);
                    gPad->SetTopMargin(0.16);
                    gPad->SetBottomMargin(0.18);
                    gPad->SetLeftMargin(0.16);
                    gPad->SetRightMargin(0.07);

                    TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
                    frame->GetXaxis()->SetTitle("#phi (deg)");
                    frame->GetYaxis()->SetTitle("F_{bin}");
                    frame->GetXaxis()->CenterTitle();
                    frame->GetYaxis()->CenterTitle();
                    frame->GetXaxis()->SetNdivisions(505);
                    frame->GetXaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetTitleSize(0.060);
                    frame->GetXaxis()->SetLabelSize(0.050);
                    frame->GetYaxis()->SetLabelSize(0.050);
                    frame->GetXaxis()->SetTitleOffset(1.15);
                    frame->GetYaxis()->SetTitleOffset(1.30);

                    // Nice degree ticks underneath.
                    drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(),
                                    gPad->GetUxmax(), 0.050);

                    // Collect Fbin points for this (xB,Q2,t) bin.
                    std::vector<double> vx, vy, exl, exh, eyl, eyh;

                    for (int r = 0; r < NR; ++r) {
                        if (!row_has_data[r]) continue;

                        double xbmin = csv.as_double(r, c_xb_min);
                        double xbmax = csv.as_double(r, c_xb_max);
                        double q2min = csv.as_double(r, c_q2_min);
                        double q2max = csv.as_double(r, c_q2_max);
                        double tmin  = csv.as_double(r, c_tab_min);
                        double tmax  = csv.as_double(r, c_tab_max);

                        if (!std::isfinite(xbmin) || !std::isfinite(xbmax) ||
                            !std::isfinite(q2min) || !std::isfinite(q2max) ||
                            !std::isfinite(tmin)  || !std::isfinite(tmax)) {
                            continue;
                        }

                        if (xbmin != xb.first || xbmax != xb.second) continue;
                        if (q2min != qR.first || q2max != qR.second) continue;
                        if (tmin  != tR.first || tmax  != tR.second) continue;

                        const std::string& fstr = csv.rows[r][c_fbin];
                        if (fstr.empty()) continue;

                        double val, s_stat, s_sys;
                        if (!parse_triple(fstr, val, s_stat, s_sys)) continue;
                        if (!std::isfinite(val) || val <= 0.0) continue;

                        double phiavg = csv.as_double(r, c_phiavg);
                        double phimin = csv.as_double(r, c_phi_min);
                        double phimax = csv.as_double(r, c_phi_max);
                        if (!std::isfinite(phiavg) ||
                            !std::isfinite(phimin) ||
                            !std::isfinite(phimax)) {
                            continue;
                        }

                        double ex_low  = std::max(0.0, phiavg - phimin);
                        double ex_high = std::max(0.0, phimax - phiavg);
                        double ey      = val;
                        double ey_err  = std::max(0.0, s_sys); // use sys as vertical error

                        vx.push_back(phiavg);
                        vy.push_back(ey);
                        exl.push_back(ex_low);
                        exh.push_back(ex_high);
                        eyl.push_back(ey_err);
                        eyh.push_back(ey_err);
                    }

                    if (vx.empty()) {
                        // No data for this pad; leave blank.
                        continue;
                    }

                    ++pads_with_data;

                    TGraphAsymmErrors* gr =
                        new TGraphAsymmErrors((int)vx.size(),
                                              vx.data(), vy.data(),
                                              exl.data(), exh.data(),
                                              eyl.data(), eyh.data());
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerSize(1.0);
                    gr->SetLineWidth(2);
                    gr->SetLineColor(kBlue+1);
                    gr->SetMarkerColor(kBlue+1);
                    gr->Draw("PE1 SAME");

                    // Reference line at Fbin = 1
                    TLine* line1 = new TLine(0.0, 1.0, 360.0, 1.0);
                    line1->SetLineStyle(2);
                    line1->SetLineColor(kGray+2);
                    line1->Draw("SAME");

                    // Panel label with Q2 and -t ranges.
                    TLatex lab;
                    lab.SetNDC();
                    lab.SetTextFont(42);
                    lab.SetTextSize(0.060);
                    lab.SetTextAlign(11);
                    std::ostringstream txt;
                    txt << "Q^{2} #in ["
                        << std::setprecision(3) << qR.first << ", "
                        << std::setprecision(3) << qR.second << "],  "
                        << "-t #in ["
                        << std::setprecision(3) << tR.first << ", "
                        << std::setprecision(3) << tR.second << "]";
                    lab.DrawLatex(0.14, 0.88, txt.str().c_str());
                } // ic
            } // ir

            // Save canvas.
            std::ostringstream fname;
            fname << "Fbin_vs_phi_" << energy_label << "_xB_" << ix << ".png";
            const fs::path outP = out_dir / fname.str();
            c->SaveAs(outP.string().c_str());
            ++canvases_made;
            std::cout << "[bincenter-plot]   Saved canvas: " << outP.string() << "\n";
            delete c;
        } // ix

        std::cout << "[bincenter-plot] Finished energy " << energy_label
                  << " GeV: " << canvases_made << " canvases with "
                  << pads_with_data << " non-empty pads.\n";
    };

    // Make plots for 10.6 GeV and 10.2 GeV.
    make_plots_for_energy("10.6", dir_10p6,
                          c_xbavg_10p6, c_phiavg_10p6, c_fbin_10p6,
                          row_has_data_10p6);

    make_plots_for_energy("10.2", dir_10p2,
                          c_xbavg_10p2, c_phiavg_10p2, c_fbin_10p2,
                          row_has_data_10p2);

    std::cout << "[bincenter-plot] Finished Fbin vs phi debug plots into "
              << out_root_dir << "/10.60 and /10.2\n";
    std::cout << "============================================================\n";
}