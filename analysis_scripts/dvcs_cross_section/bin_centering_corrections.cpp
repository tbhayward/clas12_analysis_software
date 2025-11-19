// bin_centering_corrections.cpp
// --------------------------------------------
// Bin-centering corrections Fbin written into dvcs_pass2_analysis.csv
//
// Columns filled:
//   "Fbin, 10.6 GeV"   (combined 10.6 group: Sp18/Fa18 periods)
//   "Fbin, 10.2 GeV"   (Sp19 Inb)
//
// Logic:
//   * Bins are taken from dvcs_pass2_analysis.csv (Lee-style) via bin-edge
//     columns: xBmin/xBmax, Q2min/Q2max, t_abs_min/t_abs_max, phimin/phimax.
//   * Only phi-binned rows are used (phimax - phimin < 359).
//   * We use "xBavg, 10.6 GeV" and "xBavg, Sp19 Inb" as masks to decide
//     which rows belong to each energy group.
//   * For each energy group and row, we compute Fbin using model predictions:
//       - Two models: KM15 and VGG (Helicity::Unpol).
//       - Define bin center (xB_c, Q2_c, t_pos_c, phi_c) at midpoints
//         of the bin edges.
//       - Compute sigma_center^model at the center.
//       - Compute <sigma>_bin^model as the average over a uniform
//         sub-grid of size n_substeps^4 in (xB, Q2, |t|, phi).
//       - Fbin^model = sigma_center^model / <sigma>_bin^model.
//       - Final Fbin:
//             Fbin = 0.5 * (Fbin_KM15 + Fbin_VGG)
//             stat = 0.0
//             sys  = std(Fbin_KM15, Fbin_VGG)
//         where the std is the sample standard deviation for N=2.
//   * If only one model is usable (the other fails or returns non-positive/
//     non-finite values), we fall back to that single model with sys = 0.
//     If both fail, the row is left unchanged.
//
// Parallelization:
//   * OpenMP parallel-for over (group,row) tasks, with omp_set_num_threads
//     hard-capped at 5 workers.
//
// CSV policy:
//   * We DO NOT create or rename columns.
//   * We require the columns:
//       xBmin, xBmax, Q2min, Q2max, t_abs_min, t_abs_max, phimin, phimax
//       xBavg, 10.6 GeV
//       xBavg, Sp19 Inb
//       Fbin, 10.6 GeV
//       Fbin, 10.2 GeV
//     and we fail fast if any are missing.
//   * Updated CSV is written atomically using a tmp file and rename.

#include "bin_centering_corrections.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

namespace fs = std::filesystem;

// -----------------------------
// Simple CSV document helper
// -----------------------------
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
    } //enddef

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
    } //enddef

    static double to_double(const std::string& s) {
        if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
        char* e = 0;
        double v = std::strtod(s.c_str(), &e);
        if (e == s.c_str()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return v;
    } //enddef

    bool load(const std::string& path) {
        std::ifstream fin(path.c_str());
        if (!fin.is_open()) {
            std::cerr << "[fbin] ERROR: cannot open CSV: " << path << "\n";
            return false;
        }
        std::string line;
        if (!std::getline(fin, line)) {
            std::cerr << "[fbin] ERROR: empty CSV: " << path << "\n";
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
    } //enddef

    bool save_atomic(const std::string& path) const {
        const std::string tmp = path + ".tmp";
        {
            std::ofstream fout(tmp.c_str());
            if (!fout.is_open()) {
                std::cerr << "[fbin] ERROR: cannot write CSV tmp: " << tmp << "\n";
                return false;
            }
            // header
            for (std::size_t i = 0; i < header.size(); ++i) {
                write_field(fout, header[i]);
                if (i + 1 < header.size()) fout << ',';
            }
            fout << "\n";
            // rows
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
            std::remove(path.c_str());
            fs::rename(tmp, path, ec);
            if (ec) {
                std::cerr << "[fbin] ERROR: atomic rename failed ("
                          << ec.message() << ")\n";
                return false;
            }
        }
        return true;
    } //enddef

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
}; //endstruct CsvDoc

// "(value, stat, sys)" formatting
static std::string format_triple(double v, double s_stat, double s_sys) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << "("
        << std::setprecision(8) << v << ", "
        << std::setprecision(8) << s_stat << ", "
        << std::setprecision(8) << s_sys
        << ")";
    return oss.str();
} //enddef

static inline bool finite_pos(double x) {
    return std::isfinite(x) && x > 0.0;
} //enddef

// --------------------------------------
// Binning from CSV: one RowBin per row
// --------------------------------------
struct RowBin {
    double xbmin, xbmax;
    double q2min, q2max;
    double tmin,  tmax;   // positive |t|
    double phimin, phimax;
    int    row_index;
}; //endstruct

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
        std::cerr << "[fbin] FATAL: missing one or more bin-edge columns "
                  << "(xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax)\n";
        std::exit(EXIT_FAILURE);
    }

    std::vector<RowBin> bins;
    const int NR = csv.nrows();
    bins.reserve(NR);

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

        // Skip phi-integrated bins; Fbin is defined only for phi-binned rows.
        if (std::fabs(phi_width) >= 359.0) continue;

        if (!(b.xbmax > b.xbmin &&
              b.q2max > b.q2min &&
              b.tmax  > b.tmin  &&
              b.phimax > b.phimin)) {
            continue;
        }

        bins.push_back(b);
    }

    std::cout << "[fbin] Built " << bins.size()
              << " RowBin entries out of " << csv.nrows()
              << " CSV rows (skipped phi-integrated or invalid bins).\n";

    if (bins.empty()) {
        std::cerr << "[fbin] FATAL: no usable bins found in CSV.\n";
        std::exit(EXIT_FAILURE);
    }

    return bins;
} //enddef

// For each energy group, mark which CSV rows should get Fbin.
// We only care about whether the gating column has a finite numeric value.
static std::vector<bool>
build_row_has_data(const CsvDoc& csv, const std::string& xbavg_col_name)
{
    const int NR = csv.nrows();
    int c = csv.col_index(xbavg_col_name);
    if (c < 0) {
        std::cerr << "[fbin] FATAL: missing column '" << xbavg_col_name
                  << "' needed to decide where to compute Fbin.\n";
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
} //enddef

// ---------------------------
// Model evaluation per row
// ---------------------------
//
// For a single row (one DVCS bin) and one beam energy, compute
//   Fbin_KM15 and Fbin_VGG via sub-binning in xB,Q2,|t|,phi.
//
// Inputs:
//   b         : RowBin with edges in xB, Q2, |t|, phi (deg)
//   Ebeam     : beam energy (GeV)
//   n_steps   : sub-binning steps per dimension (>= 1)
//   paths     : ModelPaths for dvcsgen/km15
//   vgg_globalfit : flag passed through to vgg_xs
//
// Outputs:
//   fbin_mean_out : combined Fbin
//   fbin_sys_out  : systematic (std between models)
// Returns:
//   true if at least one model produced a valid Fbin; false otherwise.
static bool compute_fbin_for_row(
    const RowBin& b,
    double Ebeam,
    int n_steps,
    const ModelPaths& paths,
    bool vgg_globalfit,
    double& fbin_mean_out,
    double& fbin_sys_out
) {
    fbin_mean_out = std::numeric_limits<double>::quiet_NaN();
    fbin_sys_out  = std::numeric_limits<double>::quiet_NaN();

    // Bin centers
    const double xb_c   = 0.5 * (b.xbmin  + b.xbmax);
    const double Q2_c   = 0.5 * (b.q2min  + b.q2max);
    const double tpos_c = 0.5 * (b.tmin   + b.tmax);    // positive |t|
    const double phi_c  = 0.5 * (b.phimin + b.phimax);  // mid-bin center (deg)

    if (!(xb_c > 0.0 && Q2_c > 0.0 && tpos_c > 0.0)) {
        return false;
    }

    // Center cross sections (unpolarized)
    double km15_center = km15_xs(xb_c, Q2_c, tpos_c, phi_c, Ebeam,
                                 Helicity::Unpol, paths);
    double vgg_center  = vgg_xs (xb_c, Q2_c, tpos_c, phi_c, Ebeam,
                                 Helicity::Unpol, paths, vgg_globalfit);

    const bool have_center_km15 = finite_pos(km15_center);
    const bool have_center_vgg  = finite_pos(vgg_center);

    if (!have_center_km15 && !have_center_vgg) {
        return false;
    }

    const int nx   = std::max(1, n_steps);
    const int nQ   = std::max(1, n_steps);
    const int nt   = std::max(1, n_steps);
    const int nphi = std::max(1, n_steps);

    const double dx   = (b.xbmax  - b.xbmin ) / double(nx);
    const double dQ   = (b.q2max  - b.q2min ) / double(nQ);
    const double dt   = (b.tmax   - b.tmin  ) / double(nt);
    const double dphi = (b.phimax - b.phimin) / double(nphi);

    if (!(dx > 0.0 && dQ > 0.0 && dt > 0.0 && dphi > 0.0)) {
        return false;
    }

    double sum_km15 = 0.0;
    double sum_vgg  = 0.0;
    int    cnt_k    = 0;
    int    cnt_v    = 0;

    for (int ix = 0; ix < nx; ++ix) {
        const double xB = b.xbmin + (ix + 0.5) * dx;
        for (int iQ = 0; iQ < nQ; ++iQ) {
            const double Q2 = b.q2min + (iQ + 0.5) * dQ;
            for (int it = 0; it < nt; ++it) {
                const double tpos = b.tmin + (it + 0.5) * dt;
                for (int ip = 0; ip < nphi; ++ip) {
                    const double phi = b.phimin + (ip + 0.5) * dphi;

                    if (have_center_km15) {
                        double v = km15_xs(xB, Q2, tpos, phi, Ebeam,
                                           Helicity::Unpol, paths);
                        if (finite_pos(v)) {
                            sum_km15 += v;
                            ++cnt_k;
                        }
                    }
                    if (have_center_vgg) {
                        double v = vgg_xs(xB, Q2, tpos, phi, Ebeam,
                                          Helicity::Unpol, paths, vgg_globalfit);
                        if (finite_pos(v)) {
                            sum_vgg += v;
                            ++cnt_v;
                        }
                    }
                } //endfor ip
            } //endfor it
        } //endfor iQ
    } //endfor ix

    double f_km15 = std::numeric_limits<double>::quiet_NaN();
    double f_vgg  = std::numeric_limits<double>::quiet_NaN();
    bool   ok_k   = false;
    bool   ok_v   = false;

    if (have_center_km15 && cnt_k > 0) {
        const double avg_k = sum_km15 / double(cnt_k);
        if (finite_pos(avg_k)) {
            f_km15 = km15_center / avg_k;
            if (finite_pos(f_km15)) ok_k = true;
        }
    }

    if (have_center_vgg && cnt_v > 0) {
        const double avg_v = sum_vgg / double(cnt_v);
        if (finite_pos(avg_v)) {
            f_vgg = vgg_center / avg_v;
            if (finite_pos(f_vgg)) ok_v = true;
        }
    }

    if (!ok_k && !ok_v) {
        return false;
    }

    // Combine models:
    // If both are valid:
    //   Fbin = mean, sys = std(Fbin_KM15, Fbin_VGG) with N=2.
    // If only one is valid:
    //   Fbin = that model, sys = 0.
    if (ok_k && ok_v) {
        const double mean = 0.5 * (f_km15 + f_vgg);
        const double diff = f_km15 - f_vgg;
        const double var  = (diff * diff) / 2.0; // sample variance for N=2
        const double stdv = std::sqrt(var);

        fbin_mean_out = mean;
        fbin_sys_out  = stdv;
    } else if (ok_k) {
        fbin_mean_out = f_km15;
        fbin_sys_out  = 0.0;
    } else { // ok_v only
        fbin_mean_out = f_vgg;
        fbin_sys_out  = 0.0;
    }

    return finite_pos(fbin_mean_out);
} //enddef

// -------------------------------------
// Group and fill helpers
// -------------------------------------
struct FbinGroup {
    std::string label;      // "10.6 GeV" or "10.2 GeV"
    std::string xbavg_col;  // "xBavg, 10.6 GeV" or "xBavg, Sp19 Inb"
    std::string fbin_col;   // "Fbin, 10.6 GeV" or "Fbin, 10.2 GeV"
    double      Ebeam;      // 10.6 or 10.2
}; //endstruct

struct Task {
    int group_index;   // index into groups vector
    int bin_index;     // index into RowBin vector
}; //endstruct

static bool fill_fbin_for_group(const FbinGroup& G,
                                const CsvDoc& csv,
                                const std::vector<bool>& row_has_data,
                                const std::vector<double>& fbin_vals,
                                const std::vector<double>& fbin_sys,
                                CsvDoc& csv_out)
{
    const int NR = csv.nrows();
    if ((int)row_has_data.size() != NR ||
        (int)fbin_vals.size()    != NR ||
        (int)fbin_sys.size()     != NR) {
        std::cerr << "[fbin] FATAL: size mismatch in fill_fbin_for_group for "
                  << G.label << ".\n";
        return false;
    }

    const int c_fbin = csv.col_index(G.fbin_col);
    if (c_fbin < 0) {
        std::cerr << "[fbin] FATAL: missing Fbin column '"
                  << G.fbin_col << "' in CSV header.\n";
        return false;
    }

    std::size_t written = 0;

    for (int r = 0; r < NR; ++r) {
        if (!row_has_data[r]) continue;

        const double v  = fbin_vals[r];
        const double es = fbin_sys[r];

        if (!std::isfinite(v) || v <= 0.0) {
            continue;
        }

        // Stat uncertainty is 0.0 (model-based); sys is between-model std.
        csv_out.rows[r][c_fbin] = format_triple(v, 0.0, std::max(0.0, es));
        ++written;
    }

    std::cout << "[fbin] Group " << G.label
              << ": filled column '" << G.fbin_col
              << "' ; cells written: " << written << "\n";

    return true;
} //enddef

} // anonymous namespace

// --------------------------------------------------------
// Public entry point
// --------------------------------------------------------
bool update_bin_centering_corrections_csv(
    const std::string& csv_path,
    const std::string& out_root_dir,
    int n_substeps,
    const ModelPaths& model_paths,
    bool vgg_globalfit
) {
    (void)out_root_dir; // currently unused; kept for interface symmetry

    const std::string csv_abs = fs::absolute(csv_path).string();
    std::error_code ec;
    const uintmax_t size_before =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[fbin] CSV: " << csv_abs
              << " (size=" << size_before << ")\n";

    CsvDoc csv;
    if (!csv.load(csv_path)) {
        std::cerr << "[fbin] ERROR: failed to load CSV.\n";
        return false;
    }

    // Build row bins (phi-binned only)
    std::vector<RowBin> bins = build_row_bins(csv);
    const int NR = csv.nrows();
    (void)NR;

    // Define groups: 10.6 GeV (combined) and 10.2 GeV (Sp19 Inb)
    FbinGroup g10p6;
    g10p6.label     = "10.6 GeV";
    g10p6.xbavg_col = "xBavg, 10.6 GeV";
    g10p6.fbin_col  = "Fbin, 10.6 GeV";
    g10p6.Ebeam     = 10.6;

    FbinGroup g10p2;
    g10p2.label     = "10.2 GeV";
    g10p2.xbavg_col = "xBavg, Sp19 Inb";
    g10p2.fbin_col  = "Fbin, 10.2 GeV";
    g10p2.Ebeam     = 10.2;

    std::vector<FbinGroup> groups;
    groups.push_back(g10p6);
    groups.push_back(g10p2);

    // Build row_has_data masks for each group
    std::vector<bool> row_has_data_10p6 = build_row_has_data(csv, g10p6.xbavg_col);
    std::vector<bool> row_has_data_10p2 = build_row_has_data(csv, g10p2.xbavg_col);

    // Sanity check that the Fbin columns exist up front (fail fast if not).
    if (csv.col_index(g10p6.fbin_col) < 0) {
        std::cerr << "[fbin] FATAL: missing column '" << g10p6.fbin_col
                  << "' in CSV header.\n";
        return false;
    }
    if (csv.col_index(g10p2.fbin_col) < 0) {
        std::cerr << "[fbin] FATAL: missing column '" << g10p2.fbin_col
                  << "' in CSV header.\n";
        return false;
    }

    // Prepare storage for Fbin values per group and row.
    std::vector<double> fbin_vals_10p6(NR, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> fbin_sys_10p6 (NR, 0.0);
    std::vector<double> fbin_vals_10p2(NR, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> fbin_sys_10p2 (NR, 0.0);

    // Build task list: (group, bin)
    std::vector<Task> tasks;
    tasks.reserve(bins.size() * groups.size());

    for (std::size_t ib = 0; ib < bins.size(); ++ib) {
        const RowBin& b = bins[ib];
        if (b.row_index < 0 || b.row_index >= NR) continue;

        if (row_has_data_10p6[b.row_index]) {
            Task t;
            t.group_index = 0;
            t.bin_index   = (int)ib;
            tasks.push_back(t);
        }
        if (row_has_data_10p2[b.row_index]) {
            Task t;
            t.group_index = 1;
            t.bin_index   = (int)ib;
            tasks.push_back(t);
        }
    }

    std::cout << "[fbin] Total tasks (group,row) to compute: "
              << tasks.size() << "\n";

    if (tasks.empty()) {
        std::cout << "[fbin] No rows require Fbin; nothing to do.\n";
        return true;
    }

    // Configure OpenMP (hard cap at 5 workers)
#ifdef _OPENMP
    {
        int want = omp_get_max_threads();
        int use  = std::min(5, std::max(1, want));
        omp_set_num_threads(use);
        std::cout << "[fbin] OpenMP enabled with " << use
                  << " worker(s) (hard-capped at 5)\n";
    }
#else
    std::cout << "[fbin] OpenMP not available; running single-threaded\n";
#endif

    const int n_steps = std::max(1, n_substeps);

    // Parallel loop over tasks
#pragma omp parallel for schedule(dynamic)
    for (std::size_t it = 0; it < tasks.size(); ++it) {
        const Task& T = tasks[it];
        const int gidx = T.group_index;
        const int bidx = T.bin_index;

        const RowBin& b = bins[bidx];
        double fval = std::numeric_limits<double>::quiet_NaN();
        double fsys = std::numeric_limits<double>::quiet_NaN();

        const FbinGroup& G = groups[gidx];

        bool ok = compute_fbin_for_row(
            b,
            G.Ebeam,
            n_steps,
            model_paths,
            vgg_globalfit,
            fval,
            fsys
        );

        if (!ok) continue;

        if (gidx == 0) {
            fbin_vals_10p6[b.row_index] = fval;
            fbin_sys_10p6 [b.row_index] = fsys;
        } else if (gidx == 1) {
            fbin_vals_10p2[b.row_index] = fval;
            fbin_sys_10p2 [b.row_index] = fsys;
        }
    } //end parallel for

    // Fill CSV columns for each group (sequential)
    if (!fill_fbin_for_group(g10p6,
                             csv,
                             row_has_data_10p6,
                             fbin_vals_10p6,
                             fbin_sys_10p6,
                             csv)) {
        std::cerr << "[fbin] ERROR: fill_fbin_for_group failed for "
                  << g10p6.label << ".\n";
        return false;
    }

    if (!fill_fbin_for_group(g10p2,
                             csv,
                             row_has_data_10p2,
                             fbin_vals_10p2,
                             fbin_sys_10p2,
                             csv)) {
        std::cerr << "[fbin] ERROR: fill_fbin_for_group failed for "
                  << g10p2.label << ".\n";
        return false;
    }

    // Save updated CSV atomically
    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[fbin] ERROR: failed to save updated CSV.\n";
        return false;
    }

    const uintmax_t size_after =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[fbin] Updated CSV: " << csv_abs
              << " (size " << size_before << " -> " << size_after << ")\n";
    std::cout << "[fbin] Bin-centering corrections (Fbin) complete.\n";

    return true;
} //enddef