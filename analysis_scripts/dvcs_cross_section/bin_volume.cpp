// bin_volume.cpp
// Kinematic bin volume calculator, written into CSV as "(value, stat, sys)"
// triples in the columns:
//   "bin_volume, 10.6 GeV"          phase-space-allowed volume, 10.6 GeV
//   "bin_volume, 10.2 GeV"          phase-space-allowed volume, 10.2 GeV
//   "cubic bin_volume, 10.6 GeV"    rectangular/cubic 4D volume, 10.6 GeV
//   "cubic bin_volume, 10.2 GeV"    rectangular/cubic 4D volume, 10.2 GeV
//
// Binning comes from dvcs_pass2_analysis.csv (Lee-style).
// No MC trees are used. Volumes are purely kinematic:
//
//   V_bin = ∫_{xBmin}^{xBmax} ∫_{Q2min}^{Q2max}
//           ∫_{|t|min}^{|t|max} ∫_{phi_min}^{phi_max} mask(xB,Q2,t;E) dxB dQ2 d|t| dphi
//
// where mask enforces:
//   - t > t_min(xB, Q2) (DVCS kinematics, t is negative)
//   - 0.19 < y < 0.8113919276, with y = Q2 / (2 Mp xB E_beam)
//   - W > 2.0, with W^2 = Mp^2 + Q2 (1/xB - 1)
//
// Volumes are computed via a deterministic 3D grid in (xB, Q2, t) and
// scaled by the geometric extent in phi.
//
// The analysis-note output also contains a read-only diagnostic comparing the
// current t boundary against the standard DVCS t_min expression.  This
// diagnostic does NOT alter the production bin-volume calculation.
//
// Also produces bin-volume vs phi canvases per beam energy and xB bin under
//   output/bin_volume/10.60
//   output/bin_volume/10.2
//
// All valid phi-binned rows in the pass-2 CSV are filled for both beam
// energies. Phi-integrated rows (phi width ~360 deg) are skipped.

#include "bin_volume.h"

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TROOT.h>
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
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

static inline double PI() { return 3.14159265358979323846; }

// ---------------- style bootstrap ----------------
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);
        const int rf = 42;
        gStyle->SetTitleFont(rf, "XYZ");
        gStyle->SetLabelFont(rf, "XYZ");
        gStyle->SetTextFont(rf);
    }
} _style_guard;

// ---------------- helpers: CSV ----------------

struct CsvDoc {
    std::vector<std::string> header;
    std::map<std::string,int> index;
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
            std::cerr << "[binvol] ERROR: cannot open CSV: " << path << "\n";
            return false;
        }
        std::string line;
        if (!std::getline(fin, line)) {
            std::cerr << "[binvol] ERROR: empty CSV: " << path << "\n";
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
                std::cerr << "[binvol] ERROR: cannot write CSV tmp: " << tmp << "\n";
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
                std::cerr << "[binvol] ERROR: atomic rename failed ("
                          << ec.message() << ")\n";
                return false;
            }
        }
        return true;
    }

    int nrows() const { return (int)rows.size(); }

    int col_index(const std::string& name) const {
        std::map<std::string,int>::const_iterator it = index.find(name);
        return (it == index.end()) ? -1 : it->second;
    }

    double as_double(int r, int c) const {
        if (r < 0 || r >= nrows()) return std::numeric_limits<double>::quiet_NaN();
        if (c < 0 || c >= (int)header.size()) return std::numeric_limits<double>::quiet_NaN();
        return to_double(rows[r][c]);
    }
};

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

// ---------------- deterministic 3D (xB,Q2,t) volume under kinematic masks ----------------

static double calculate_cubic_bin_volume(double xB_min, double xB_max,
                                      double Q2_min, double Q2_max,
                                      double t_abs_min, double t_abs_max,
                                      double phi_min, double phi_max)
{
    return
        (xB_max - xB_min) *
        (Q2_max - Q2_min) *
        (t_abs_max - t_abs_min) *
        (phi_max - phi_min);
}

static double calculate_phase_space_allowed_fraction(double xB_min, double xB_max,
                                                     double Q2_min, double Q2_max,
                                                     double t_abs_min, double t_abs_max,
                                                     double E_beam)
{
    constexpr int    n_steps = 10;
    constexpr double Mp      = 0.938272; // Proton mass (GeV)
    int valid_count = 0;

    // Convert |t| edges to physical t<0 for the loop.
    const double t_phys_min = -t_abs_max;
    const double t_phys_max = -t_abs_min;

    const double dxB = (xB_max - xB_min) / n_steps;
    const double dQ2 = (Q2_max - Q2_min) / n_steps;
    const double dt  = (t_phys_max - t_phys_min) / n_steps;

    for (int i = 0; i < n_steps; ++i) {
        const double xB = xB_min + (i + 0.5) * dxB;

        for (int j = 0; j < n_steps; ++j) {
            const double Q2 = Q2_min + (j + 0.5) * dQ2;

            // These quantities do not depend on t, so evaluate them once per
            // (xB,Q2) grid point rather than once for all ten t samples.
            const double y  = Q2 / (2.0 * Mp * xB * E_beam);
            const double W2 = Mp * Mp + Q2 * (1.0 / xB - 1.0);
            const double W  = (W2 > 0.0) ? std::sqrt(W2) : 0.0;
            if (!(y > 0.19 && y < 0.8113919276 && W > 2.0)) {
                continue;
            }

            // DVCS t_min(xB,Q2) (negative).
            const double sqrt_term = std::sqrt(1.0 + (4.0 * Mp * Mp * xB * xB) / Q2);
            const double t_min_val = -Q2 * (1.0 - xB) * (1.0 - xB) /
                                     (xB * (1.0 + sqrt_term));

            for (int k = 0; k < n_steps; ++k) {
                const double t = t_phys_min + (k + 0.5) * dt;
                if (t > t_min_val) {
                    ++valid_count;
                }
            }
        }
    }

    return static_cast<double>(valid_count) /
           static_cast<double>(n_steps * n_steps * n_steps);
}


// -----------------------------------------------------------------------------
// Diagnostic comparison of the currently implemented t boundary with the
// standard DVCS t_min expression.  These helpers are deliberately separate from
// the production volume calculation so this study cannot alter the CSV values.
// -----------------------------------------------------------------------------

static double binvol_current_t_boundary(double xB, double Q2)
{
    constexpr double Mp = 0.938272;
    const double sqrt_term =
        std::sqrt(1.0 + (4.0 * Mp * Mp * xB * xB) / Q2);

    return -Q2 * (1.0 - xB) * (1.0 - xB) /
           (xB * (1.0 + sqrt_term));
}

static double binvol_standard_dvcs_tmin(double xB, double Q2)
{
    constexpr double Mp = 0.938272;

    const double eps2 = 4.0 * Mp * Mp * xB * xB / Q2;
    const double root = std::sqrt(1.0 + eps2);

    const double minus_tmin =
        Q2 *
        (2.0 * (1.0 - xB) * (1.0 - root) + eps2) /
        (4.0 * xB * (1.0 - xB) + eps2);

    return -minus_tmin;
}

enum class BinvolTBoundaryMode {
    CurrentImplementation,
    StandardDVCS
};

static double calculate_phase_space_allowed_fraction_t_diagnostic(
    double xB_min, double xB_max,
    double Q2_min, double Q2_max,
    double t_abs_min, double t_abs_max,
    double E_beam,
    BinvolTBoundaryMode mode)
{
    constexpr int n_steps = 10;
    constexpr double Mp = 0.938272;

    int valid_count = 0;

    const double t_phys_min = -t_abs_max;
    const double t_phys_max = -t_abs_min;

    const double dxB = (xB_max - xB_min) / n_steps;
    const double dQ2 = (Q2_max - Q2_min) / n_steps;
    const double dt  = (t_phys_max - t_phys_min) / n_steps;

    for (int i = 0; i < n_steps; ++i) {
        const double xB = xB_min + (i + 0.5) * dxB;

        for (int j = 0; j < n_steps; ++j) {
            const double Q2 = Q2_min + (j + 0.5) * dQ2;

            // Keep every non-t cut exactly the same as the production
            // bin-volume calculation.  The only changed ingredient in this
            // diagnostic is the t-boundary definition/inequality.
            const double y  = Q2 / (2.0 * Mp * xB * E_beam);
            const double W2 = Mp * Mp + Q2 * (1.0 / xB - 1.0);
            const double W  = (W2 > 0.0) ? std::sqrt(W2) : 0.0;

            if (!(y > 0.19 && y < 0.8113919276 && W > 2.0)) {
                continue;
            }

            const double boundary =
                (mode == BinvolTBoundaryMode::CurrentImplementation)
                    ? binvol_current_t_boundary(xB, Q2)
                    : binvol_standard_dvcs_tmin(xB, Q2);

            for (int k = 0; k < n_steps; ++k) {
                const double t = t_phys_min + (k + 0.5) * dt;

                bool accepted = false;

                if (mode == BinvolTBoundaryMode::CurrentImplementation) {
                    // Exactly the condition used by the current production code.
                    accepted = (t > boundary);
                } else {
                    // Standard DVCS physical condition:
                    // |t| >= |t_min|  <=>  t <= t_min for negative t.
                    accepted = (t <= boundary);
                }

                if (accepted) {
                    ++valid_count;
                }
            }
        }
    }

    return static_cast<double>(valid_count) /
           static_cast<double>(n_steps * n_steps * n_steps);
}

// ---------------- helpers for deciding which rows to use ----------------

// For each energy group, mark which CSV rows actually have data based on xBavg column.
static std::vector<bool>
build_row_has_data(const CsvDoc& csv, const std::string& xbavg_col_name)
{
    const int NR = csv.nrows();
    int c = csv.col_index(xbavg_col_name);
    if (c < 0) {
        std::cerr << "[binvol] FATAL: missing column '" << xbavg_col_name
                  << "' needed to decide where to compute bin volume.\n";
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

// ---------------- plotting helpers ----------------

static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize){
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

struct CellData {
    std::vector<double> X;    // phi position (deg)
    std::vector<double> Y;    // bin volume
    std::vector<double> EXL;  // left x error (deg)
    std::vector<double> EXH;  // right x error (deg)
    std::vector<double> EY;   // y error (stat, symmetric; here zero)
};

struct VolumeGroup {
    std::string label;        // "10.6 GeV" or "10.2 GeV" for titles
    std::string xbavg_col;    // "xBavg, 10.6 GeV" or "xBavg, Sp19 Inb"
    std::string phiavg_col;   // "phiavg, 10.6 GeV" or "phiavg, Sp19 Inb"
    std::string q2avg_col;    // "Q2avg, 10.6 GeV" or "Q2avg, Sp19 Inb"
    std::string tabavg_col;   // "t_abs_avg, 10.6 GeV" or "t_abs_avg, Sp19 Inb"
    std::string binvol_col;         // "bin_volume, 10.6 GeV" or "bin_volume, 10.2 GeV"
    std::string cubic_binvol_col;   // "cubic bin_volume, 10.6 GeV" or "cubic bin_volume, 10.2 GeV"
    std::string energy_dir;         // "10.60" or "10.2" for directory names
    double Ebeam;             // beam energy (GeV)
};

struct PhaseSpaceCellKey {
    double xbmin;
    double xbmax;
    double q2min;
    double q2max;
    double tmin;
    double tmax;

    bool operator<(const PhaseSpaceCellKey& other) const {
        return std::tie(xbmin, xbmax, q2min, q2max, tmin, tmax) <
               std::tie(other.xbmin, other.xbmax,
                        other.q2min, other.q2max,
                        other.tmin, other.tmax);
    }
};

// Compute bin volumes for a given group and fill its CSV column.
static void compute_bin_volumes_for_group(const VolumeGroup& G,
                                          CsvDoc& csv,
                                          const std::vector<bool>& row_has_data)
{
    const int NR = csv.nrows();
    if ((int)row_has_data.size() != NR) {
        std::cerr << "[binvol] FATAL: row_has_data size mismatch for group "
                  << G.label << ".\n";
        std::exit(EXIT_FAILURE);
    }

    const int c_xb_min  = csv.col_index("xBmin");
    const int c_xb_max  = csv.col_index("xBmax");
    const int c_q2_min  = csv.col_index("Q2min");
    const int c_q2_max  = csv.col_index("Q2max");
    const int c_tab_min = csv.col_index("t_abs_min");
    const int c_tab_max = csv.col_index("t_abs_max");
    const int c_phi_min = csv.col_index("phimin");
    const int c_phi_max = csv.col_index("phimax");
    const int c_binvol       = csv.col_index(G.binvol_col);
    const int c_cubic_binvol = csv.col_index(G.cubic_binvol_col);

    if (c_xb_min < 0 || c_xb_max < 0 ||
        c_q2_min < 0 || c_q2_max < 0 ||
        c_tab_min < 0 || c_tab_max < 0 ||
        c_phi_min < 0 || c_phi_max < 0) {
        std::cerr << "[binvol] FATAL: missing one or more bin-edge columns "
                  << "(xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax)\n";
        std::exit(EXIT_FAILURE);
    }

    if (c_binvol < 0) {
        std::cerr << "[binvol] FATAL: missing phase-space bin volume column '"
                  << G.binvol_col << "' in CSV header.\n";
        std::exit(EXIT_FAILURE);
    }

    if (c_cubic_binvol < 0) {
        std::cerr << "[binvol] FATAL: missing cubic bin volume column '"
                  << G.cubic_binvol_col << "' in CSV header.\n";
        std::exit(EXIT_FAILURE);
    }

    int bins_written = 0;

    // All phi rows belonging to the same (xB,Q2,|t|) cell have the same
    // phase-space-allowed fraction. Cache the deterministic 10^3-grid result
    // once per 3D cell and apply the row-specific phi width afterward.
    std::map<PhaseSpaceCellKey, double> allowed_fraction_cache;

    for (int r = 0; r < NR; ++r) {
        if (!row_has_data[r]) continue;

        const double xbmin  = csv.as_double(r, c_xb_min);
        const double xbmax  = csv.as_double(r, c_xb_max);
        const double q2min  = csv.as_double(r, c_q2_min);
        const double q2max  = csv.as_double(r, c_q2_max);
        const double tmin   = csv.as_double(r, c_tab_min);
        const double tmax   = csv.as_double(r, c_tab_max);
        const double phimin = csv.as_double(r, c_phi_min);
        const double phimax = csv.as_double(r, c_phi_max);

        if (!std::isfinite(xbmin)  || !std::isfinite(xbmax) ||
            !std::isfinite(q2min)  || !std::isfinite(q2max) ||
            !std::isfinite(tmin)   || !std::isfinite(tmax)  ||
            !std::isfinite(phimin) || !std::isfinite(phimax)) {
            continue;
        }

        const double phi_width = phimax - phimin;
        // Skip phi-integrated bins (we only define bin volume for phi-binned rows).
        if (std::fabs(phi_width) >= 359.0) {
            continue;
        }

        if (!(xbmax > xbmin &&
              q2max > q2min &&
              tmax  > tmin)) {
            continue;
        }

        const double phi_min_rad = phimin * PI() / 180.0;
        const double phi_max_rad = phimax * PI() / 180.0;

        const double V_cubic = calculate_cubic_bin_volume(
            xbmin, xbmax,
            q2min, q2max,
            tmin,  tmax,
            phi_min_rad, phi_max_rad
        );

        const PhaseSpaceCellKey key = {
            xbmin, xbmax, q2min, q2max, tmin, tmax
        };

        double allowed_fraction = 0.0;
        std::map<PhaseSpaceCellKey, double>::const_iterator fit =
            allowed_fraction_cache.find(key);
        if (fit == allowed_fraction_cache.end()) {
            allowed_fraction = calculate_phase_space_allowed_fraction(
                xbmin, xbmax,
                q2min, q2max,
                tmin,  tmax,
                G.Ebeam
            );
            allowed_fraction_cache.emplace(key, allowed_fraction);
        } else {
            allowed_fraction = fit->second;
        }

        const double V_allowed = V_cubic * allowed_fraction;

        if (!std::isfinite(V_cubic) || !std::isfinite(V_allowed)) {
            continue;
        }

        double Vcubic_clamped = V_cubic;
        if (Vcubic_clamped < 0.0) Vcubic_clamped = 0.0;

        double Vallowed_clamped = V_allowed;
        if (Vallowed_clamped < 0.0) Vallowed_clamped = 0.0;

        csv.rows[r][c_binvol] = format_triple(Vallowed_clamped, 0.0, 0.0);
        csv.rows[r][c_cubic_binvol] = format_triple(Vcubic_clamped, 0.0, 0.0);
        ++bins_written;
    }

    std::cout << "[binvol] Group " << G.label
              << ": filled columns '" << G.binvol_col
              << "' and '" << G.cubic_binvol_col
              << "' for " << bins_written << " phi-binned rows.\n";
}

// Draw bin-volume vs phi canvases for a given group.
static void draw_bin_volume_canvases(const VolumeGroup& G,
                                     const CsvDoc& csv,
                                     const std::vector<bool>& row_has_data,
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

    const int c_phiavg  = csv.col_index(G.phiavg_col); // may be -1
    const int c_binvol  = csv.col_index(G.binvol_col);

    if (c_xb_min < 0 || c_xb_max < 0 ||
        c_q2_min < 0 || c_q2_max < 0 ||
        c_tab_min < 0 || c_tab_max < 0 ||
        c_phi_min < 0 || c_phi_max < 0) {
        std::cerr << "[binvol] FATAL: missing bin-edge columns needed for plotting.\n";
        std::exit(EXIT_FAILURE);
    }

    if (c_binvol < 0) {
        std::cerr << "[binvol] FATAL: missing bin volume column '"
                  << G.binvol_col << "' for plotting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Unique xB bins for this group
    std::set<std::pair<double,double> > xb_set;
    for (int r = 0; r < NR; ++r) {
        if (!row_has_data[r]) continue;

        const double xbmin  = csv.as_double(r, c_xb_min);
        const double xbmax  = csv.as_double(r, c_xb_max);
        const double phimin = csv.as_double(r, c_phi_min);
        const double phimax = csv.as_double(r, c_phi_max);
        if (!std::isfinite(xbmin) || !std::isfinite(xbmax) ||
            !std::isfinite(phimin) || !std::isfinite(phimax)) {
            continue;
        }
        const double phi_width = phimax - phimin;
        if (std::fabs(phi_width) >= 359.0) continue; // skip phi-integrated

        xb_set.insert(std::make_pair(xbmin, xbmax));
    }

    if (xb_set.empty()) {
        std::cout << "[binvol] Group " << G.label
                  << ": no phi-binned rows with data; no plots produced.\n";
        return;
    }

    const double head_size = 0.20;
    const double label_sz  = 0.048;
    const double title_sz  = 0.060;

    const std::string base_dir =
        (fs::path(out_root_dir) / "bin_volume" / G.energy_dir).string();
    // Directories are assumed to already exist from makeOutputDirs().

    for (std::set<std::pair<double,double> >::const_iterator itX = xb_set.begin();
         itX != xb_set.end(); ++itX) {

        const std::pair<double,double> xb = *itX;

        // Unique Q2 and t ranges for this xB and group
        std::set<std::pair<double,double> > q2set;
        std::set<std::pair<double,double> > tset;

        for (int r = 0; r < NR; ++r) {
            if (!row_has_data[r]) continue;

            const double xbmin  = csv.as_double(r, c_xb_min);
            const double xbmax  = csv.as_double(r, c_xb_max);
            const double q2min  = csv.as_double(r, c_q2_min);
            const double q2max  = csv.as_double(r, c_q2_max);
            const double tmin   = csv.as_double(r, c_tab_min);
            const double tmax   = csv.as_double(r, c_tab_max);
            const double phimin = csv.as_double(r, c_phi_min);
            const double phimax = csv.as_double(r, c_phi_max);

            if (!std::isfinite(xbmin) || !std::isfinite(xbmax) ||
                !std::isfinite(q2min) || !std::isfinite(q2max) ||
                !std::isfinite(tmin)  || !std::isfinite(tmax)  ||
                !std::isfinite(phimin) || !std::isfinite(phimax)) {
                continue;
            }

            const double phi_width = phimax - phimin;
            if (std::fabs(phi_width) >= 359.0) continue;

            if (std::fabs(xbmin - xb.first) < 1e-9 &&
                std::fabs(xbmax - xb.second) < 1e-9) {
                q2set.insert(std::make_pair(q2min, q2max));
                tset.insert(std::make_pair(tmin, tmax));
            }
        }

        std::vector<std::pair<double,double> > Q2s(q2set.begin(), q2set.end());
        std::vector<std::pair<double,double> > Ts (tset.begin(),  tset.end());
        if (Q2s.empty() || Ts.empty()) continue;

        // Old bin_volume layout: rows = t bins, cols = Q2 bins
        const int nrows = (int)Ts.size();
        const int ncols = (int)Q2s.size();

        const int W = 280 * ncols + 160;
        const int H = 240 * nrows + 170;

        std::ostringstream cname;
        cname << "c_binvol_" << G.energy_dir << "_xB_"
              << (int)std::round(xb.first * 1000.0);
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0);
        pTop->SetBorderSize(0);
        pTop->Draw();

        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0);
        pGrid->SetBorderSize(0);
        pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // Title
        pTop->cd();
        TLatex head;
        head.SetNDC();
        head.SetTextAlign(22);
        head.SetTextFont(42);
        head.SetTextSize(head_size);
        std::ostringstream tit;
        tit << "Bin Volume (kinematic)  " << G.label
            << "   x_{B} #in [" << std::setprecision(2) << xb.first
            << ", " << std::setprecision(2) << xb.second << "]";
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        // Panels
        for (int rrow = 0; rrow < nrows; ++rrow) {
            const std::pair<double,double>& tpair = Ts[rrow];

            for (int ccol = 0; ccol < ncols; ++ccol) {
                const std::pair<double,double>& qpair = Q2s[ccol];

                pGrid->cd(rrow * ncols + ccol + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 0.0); // y autoscale
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetLabelSize(0.0001); // hide default numeric labels
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("Kinematic bin volume");
                ax->CenterTitle();
                ay->CenterTitle();
                ax->SetNdivisions(505);
                ax->SetTitleSize(title_sz);
                ay->SetTitleSize(title_sz);
                ay->SetLabelSize(label_sz);
                ax->SetTitleOffset(1.25);
                ay->SetTitleOffset(1.35);

                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), label_sz);

                std::vector<int> rows_for_cell;
                rows_for_cell.reserve(16);

                // Collect all rows for this (xB,Q2,t) cell
                for (int r = 0; r < NR; ++r) {
                    if (!row_has_data[r]) continue;

                    const double xbmin  = csv.as_double(r, c_xb_min);
                    const double xbmax  = csv.as_double(r, c_xb_max);
                    const double q2min  = csv.as_double(r, c_q2_min);
                    const double q2max  = csv.as_double(r, c_q2_max);
                    const double tmin   = csv.as_double(r, c_tab_min);
                    const double tmax   = csv.as_double(r, c_tab_max);
                    const double phimin = csv.as_double(r, c_phi_min);
                    const double phimax = csv.as_double(r, c_phi_max);

                    if (!std::isfinite(xbmin) || !std::isfinite(xbmax) ||
                        !std::isfinite(q2min) || !std::isfinite(q2max) ||
                        !std::isfinite(tmin)  || !std::isfinite(tmax)  ||
                        !std::isfinite(phimin) || !std::isfinite(phimax)) {
                        continue;
                    }

                    const double phi_width = phimax - phimin;
                    if (std::fabs(phi_width) >= 359.0) continue;

                    if (std::fabs(xbmin - xb.first) < 1e-9 &&
                        std::fabs(xbmax - xb.second) < 1e-9 &&
                        std::fabs(q2min - qpair.first) < 1e-9 &&
                        std::fabs(q2max - qpair.second) < 1e-9 &&
                        std::fabs(tmin  - tpair.first) < 1e-9 &&
                        std::fabs(tmax  - tpair.second) < 1e-9) {
                        rows_for_cell.push_back(r);
                    }
                }

                if (rows_for_cell.empty()) {
                    // Leave pad blank
                    continue;
                }

                std::sort(rows_for_cell.begin(), rows_for_cell.end(),
                          [&](int a, int b) {
                              return csv.as_double(a, c_phi_min) <
                                     csv.as_double(b, c_phi_min);
                          });

                CellData C;
                C.X.reserve(rows_for_cell.size());
                C.EXL.reserve(rows_for_cell.size());
                C.EXH.reserve(rows_for_cell.size());
                C.EY.reserve(rows_for_cell.size());
                C.Y.reserve(rows_for_cell.size());

                double ymax = 0.0;

                for (std::size_t k = 0; k < rows_for_cell.size(); ++k) {
                    int r = rows_for_cell[k];

                    const double pmin = csv.as_double(r, c_phi_min);
                    const double pmax = csv.as_double(r, c_phi_max);
                    if (!std::isfinite(pmin) || !std::isfinite(pmax)) continue;

                    double xphi = 0.5 * (pmin + pmax);
                    if (c_phiavg >= 0) {
                        const double pav = csv.as_double(r, c_phiavg);
                        if (std::isfinite(pav) && pav > 0.0 && pav < 360.0) {
                            xphi = pav;
                        }
                    }

                    const std::string& cell = csv.rows[r][c_binvol];
                    if (cell.empty()) continue;

                    double V = 0.0;
                    double s_stat = 0.0;
                    double s_sys  = 0.0;
                    if (!parse_triple(cell, V, s_stat, s_sys)) {
                        double tmp = CsvDoc::to_double(cell);
                        if (!std::isfinite(tmp)) continue;
                        V = tmp;
                    }

                    if (!std::isfinite(V)) continue;

                    double exl = xphi - pmin;
                    double exh = pmax - xphi;
                    if (exl < 0.0) exl = 0.0;
                    if (exh < 0.0) exh = 0.0;

                    C.X.push_back(xphi);
                    C.Y.push_back(V);
                    C.EXL.push_back(exl);
                    C.EXH.push_back(exh);
                    C.EY.push_back(0.0); // deterministic volume, no stat error

                    if (V > ymax) ymax = V;
                }

                if (C.X.empty()) {
                    continue;
                }

                if (ymax <= 0.0) ymax = 1.0;
                frame->GetYaxis()->SetRangeUser(0.0, ymax * 1.20);

                TGraphAsymmErrors* gvol = new TGraphAsymmErrors(
                    (int)C.X.size(),
                    (double*)C.X.data(),
                    (double*)C.Y.data(),
                    (double*)C.EXL.data(),
                    (double*)C.EXH.data(),
                    (double*)C.EY.data(),
                    (double*)C.EY.data()
                );

                gvol->SetMarkerStyle(20);
                gvol->SetMarkerSize(1.0);
                gvol->SetLineWidth(2);

                gvol->Draw("PE1 SAME");

                // annotate Q2 and -t (edges, like original bin_volume)
                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.040);
                lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         qpair.first, qpair.second,
                         tpair.first, tpair.second));
            }
        }

        const std::string fpath =
            (fs::path(base_dir) /
             ("plot_bin_volume_" + G.energy_dir + "_xB_" +
              std::to_string((int)std::round(xb.first * 1000.0)) +
              ".png")).string();

        c->SaveAs(fpath.c_str());
        delete c;
    }
}

// -----------------------------------------------------------------------------
// Analysis-note diagnostics
// -----------------------------------------------------------------------------

static double note_quantile(std::vector<double> v, double q)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    if (v.size() == 1) return v[0];
    const double x = q * double(v.size() - 1);
    const std::size_t i = static_cast<std::size_t>(std::floor(x));
    const std::size_t j = std::min(i + 1, v.size() - 1);
    const double f = x - double(i);
    return v[i] * (1.0 - f) + v[j] * f;
}

static bool note_same_edge(double a, double b)
{
    return std::fabs(a - b) < 1e-10;
}

struct NoteVolumePoint {
    double xbmin=0, xbmax=0;
    double q2min=0, q2max=0;
    double tmin=0, tmax=0;
    double phimin=0, phimax=0;
    double allowed=0;
    double cubic=0;
    double fraction=0;
};

static std::vector<NoteVolumePoint>
collect_note_volume_points(const CsvDoc& csv, const VolumeGroup& G)
{
    const int c_x0=csv.col_index("xBmin");
    const int c_x1=csv.col_index("xBmax");
    const int c_q0=csv.col_index("Q2min");
    const int c_q1=csv.col_index("Q2max");
    const int c_t0=csv.col_index("t_abs_min");
    const int c_t1=csv.col_index("t_abs_max");
    const int c_p0=csv.col_index("phimin");
    const int c_p1=csv.col_index("phimax");
    const int c_v =csv.col_index(G.binvol_col);
    const int c_c =csv.col_index(G.cubic_binvol_col);

    std::vector<NoteVolumePoint> out;
    if(c_x0<0||c_x1<0||c_q0<0||c_q1<0||c_t0<0||c_t1<0||
       c_p0<0||c_p1<0||c_v<0||c_c<0) return out;

    for(int r=0;r<csv.nrows();++r){
        const double p0=csv.as_double(r,c_p0);
        const double p1=csv.as_double(r,c_p1);
        if(!std::isfinite(p0)||!std::isfinite(p1)) continue;
        if(std::fabs(p1-p0) > 100.0) continue; // skip phi-integrated rows

        double v=0,vs=0,vy=0,c=0,cs=0,cy=0;
        if(!parse_triple(csv.rows[r][c_v],v,vs,vy)){
            v=CsvDoc::to_double(csv.rows[r][c_v]);
        }
        if(!parse_triple(csv.rows[r][c_c],c,cs,cy)){
            c=CsvDoc::to_double(csv.rows[r][c_c]);
        }
        if(!std::isfinite(v)||!std::isfinite(c)||c<=0||v<0) continue;

        NoteVolumePoint p;
        p.xbmin=csv.as_double(r,c_x0); p.xbmax=csv.as_double(r,c_x1);
        p.q2min=csv.as_double(r,c_q0); p.q2max=csv.as_double(r,c_q1);
        p.tmin=csv.as_double(r,c_t0); p.tmax=csv.as_double(r,c_t1);
        p.phimin=p0; p.phimax=p1;
        p.allowed=v; p.cubic=c;
        p.fraction=std::max(0.0,std::min(1.0,v/c));
        if(std::isfinite(p.xbmin)&&std::isfinite(p.xbmax)&&
           std::isfinite(p.q2min)&&std::isfinite(p.q2max)&&
           std::isfinite(p.tmin)&&std::isfinite(p.tmax)) out.push_back(p);
    }
    return out;
}

static void write_bin_volume_analysis_note_outputs(
    const CsvDoc& csv,
    const VolumeGroup& g10p6,
    const VolumeGroup& g10p2,
    const std::string& out_root_dir)
{
    namespace fs=std::filesystem;
    const fs::path note_dir=fs::path(out_root_dir)/"bin_volume"/"analysis_note";
    std::error_code ec;
    fs::create_directories(note_dir,ec);
    if(ec){
        std::cerr<<"[binvol-note] WARNING: cannot create "<<note_dir
                 <<": "<<ec.message()<<"\n";
        return;
    }

    const std::vector<NoteVolumePoint> p106=collect_note_volume_points(csv,g10p6);
    const std::vector<NoteVolumePoint> p102=collect_note_volume_points(csv,g10p2);
    if(p106.empty()||p102.empty()){
        std::cerr<<"[binvol-note] WARNING: no populated bin-volume points found.\n";
        return;
    }

    auto unique_bulk=[](const std::vector<NoteVolumePoint>& pts){
        std::vector<NoteVolumePoint> out;
        std::set<std::tuple<double,double,double,double,double,double>> seen;
        for(const auto&p:pts){
            auto key=std::make_tuple(p.xbmin,p.xbmax,p.q2min,p.q2max,p.tmin,p.tmax);
            if(seen.insert(key).second) out.push_back(p);
        }
        return out;
    };

    // Since the phase-space fraction does not depend on phi, summarize one
    // entry per (xB,Q2,|t|) cell rather than counting each phi bin 24 times.
    const std::vector<NoteVolumePoint> bulk106=unique_bulk(p106);
    const std::vector<NoteVolumePoint> bulk102=unique_bulk(p102);


    // ---------------------------------------------------------------------
    // t-boundary diagnostic.
    //
    // This does not feed back into the production volume calculation.  It
    // compares the current t_min_val implementation with the standard DVCS
    // t_min expression while holding the y and W selections fixed.
    // ---------------------------------------------------------------------
    {
        struct TDiagPoint {
            double xb=0.0;
            double q2=0.0;
            double current_t=0.0;
            double standard_t=0.0;
            double frac_current_106=0.0;
            double frac_standard_106=0.0;
            double frac_current_102=0.0;
            double frac_standard_102=0.0;
        };

        std::vector<TDiagPoint> diag;
        diag.reserve(bulk106.size());

        for (const auto& p : bulk106) {
            const double xb = 0.5 * (p.xbmin + p.xbmax);
            const double q2 = 0.5 * (p.q2min + p.q2max);

            TDiagPoint d;
            d.xb = xb;
            d.q2 = q2;
            d.current_t = binvol_current_t_boundary(xb, q2);
            d.standard_t = binvol_standard_dvcs_tmin(xb, q2);

            d.frac_current_106 =
                calculate_phase_space_allowed_fraction_t_diagnostic(
                    p.xbmin,p.xbmax,p.q2min,p.q2max,p.tmin,p.tmax,
                    g10p6.Ebeam,
                    BinvolTBoundaryMode::CurrentImplementation);

            d.frac_standard_106 =
                calculate_phase_space_allowed_fraction_t_diagnostic(
                    p.xbmin,p.xbmax,p.q2min,p.q2max,p.tmin,p.tmax,
                    g10p6.Ebeam,
                    BinvolTBoundaryMode::StandardDVCS);

            d.frac_current_102 =
                calculate_phase_space_allowed_fraction_t_diagnostic(
                    p.xbmin,p.xbmax,p.q2min,p.q2max,p.tmin,p.tmax,
                    g10p2.Ebeam,
                    BinvolTBoundaryMode::CurrentImplementation);

            d.frac_standard_102 =
                calculate_phase_space_allowed_fraction_t_diagnostic(
                    p.xbmin,p.xbmax,p.q2min,p.q2max,p.tmin,p.tmax,
                    g10p2.Ebeam,
                    BinvolTBoundaryMode::StandardDVCS);

            diag.push_back(d);
        }

        if (!diag.empty()) {
            // Full machine-readable diagnostic.
            {
                std::ofstream o((note_dir/"tmin_diagnostic_bins.csv").string());
                o<<"xB_center,Q2_center,current_t_boundary,standard_dvcs_tmin,"
                    "allowed_fraction_current_10p6,allowed_fraction_standard_10p6,"
                    "allowed_fraction_current_10p2,allowed_fraction_standard_10p2\n";
                o<<std::setprecision(10);
                for(const auto&d:diag){
                    o<<d.xb<<","<<d.q2<<","
                     <<d.current_t<<","<<d.standard_t<<","
                     <<d.frac_current_106<<","<<d.frac_standard_106<<","
                     <<d.frac_current_102<<","<<d.frac_standard_102<<"\n";
                }
            }

            // Compact numerical summary of the impact on the allowed fraction.
            std::vector<double> delta106,delta102,ratio_t;
            int same106=0,same102=0;
            int diff_gt_001_106=0,diff_gt_001_102=0;
            int diff_gt_005_106=0,diff_gt_005_102=0;

            for(const auto&d:diag){
                const double a106=std::fabs(d.frac_current_106-d.frac_standard_106);
                const double a102=std::fabs(d.frac_current_102-d.frac_standard_102);

                delta106.push_back(a106);
                delta102.push_back(a102);

                if(a106<1e-12) ++same106;
                if(a102<1e-12) ++same102;
                if(a106>0.01) ++diff_gt_001_106;
                if(a102>0.01) ++diff_gt_001_102;
                if(a106>0.05) ++diff_gt_005_106;
                if(a102>0.05) ++diff_gt_005_102;

                if(std::fabs(d.standard_t)>0.0){
                    ratio_t.push_back(std::fabs(d.current_t/d.standard_t));
                }
            }

            {
                std::ofstream o((note_dir/"tmin_diagnostic_summary.csv").string());
                o<<"quantity,value\n";
                o<<std::setprecision(10);
                o<<"cells,"<<diag.size()<<"\n";
                o<<"median_abs_current_over_standard_t_boundary,"
                 <<note_quantile(ratio_t,.50)<<"\n";
                o<<"p16_abs_current_over_standard_t_boundary,"
                 <<note_quantile(ratio_t,.16)<<"\n";
                o<<"p84_abs_current_over_standard_t_boundary,"
                 <<note_quantile(ratio_t,.84)<<"\n";

                o<<"median_abs_allowed_fraction_difference_10p6,"
                 <<note_quantile(delta106,.50)<<"\n";
                o<<"p84_abs_allowed_fraction_difference_10p6,"
                 <<note_quantile(delta106,.84)<<"\n";
                o<<"max_abs_allowed_fraction_difference_10p6,"
                 <<*std::max_element(delta106.begin(),delta106.end())<<"\n";
                o<<"identical_fraction_cells_10p6,"<<same106<<"\n";
                o<<"cells_difference_gt_0p01_10p6,"<<diff_gt_001_106<<"\n";
                o<<"cells_difference_gt_0p05_10p6,"<<diff_gt_005_106<<"\n";

                o<<"median_abs_allowed_fraction_difference_10p2,"
                 <<note_quantile(delta102,.50)<<"\n";
                o<<"p84_abs_allowed_fraction_difference_10p2,"
                 <<note_quantile(delta102,.84)<<"\n";
                o<<"max_abs_allowed_fraction_difference_10p2,"
                 <<*std::max_element(delta102.begin(),delta102.end())<<"\n";
                o<<"identical_fraction_cells_10p2,"<<same102<<"\n";
                o<<"cells_difference_gt_0p01_10p2,"<<diff_gt_001_102<<"\n";
                o<<"cells_difference_gt_0p05_10p2,"<<diff_gt_005_102<<"\n";
            }

            // Choose three actual Q2-bin centers spanning the available range.
            std::vector<double> q2vals;
            for(const auto&d:diag){
                bool seen=false;
                for(double q:q2vals){
                    if(std::fabs(q-d.q2)<1e-10){ seen=true; break; }
                }
                if(!seen) q2vals.push_back(d.q2);
            }
            std::sort(q2vals.begin(),q2vals.end());

            std::vector<double> chosen_q2;
            if(!q2vals.empty()){
                chosen_q2.push_back(q2vals.front());
                if(q2vals.size()>2) chosen_q2.push_back(q2vals[q2vals.size()/2]);
                if(q2vals.size()>1) chosen_q2.push_back(q2vals.back());
            }

            TCanvas c("c_tmin_diag","",1500,650);
            c.Divide(2,1,0.004,0.004);

            // Left panel: directly compare the two boundary expressions.
            c.cd(1);
            gPad->SetLeftMargin(.15);
            gPad->SetRightMargin(.035);
            gPad->SetBottomMargin(.14);
            gPad->SetTopMargin(.12);
            gPad->SetGridy();
            gPad->SetTicks(1,1);

            double max_abs_boundary=0.0;
            for(const auto&d:diag){
                max_abs_boundary=std::max(max_abs_boundary,std::fabs(d.current_t));
                max_abs_boundary=std::max(max_abs_boundary,std::fabs(d.standard_t));
            }
            if(max_abs_boundary<=0) max_abs_boundary=1.0;

            TH1F f1("h_tmin_diag_formula","",100,.05,.60);
            f1.SetMinimum(0.0);
            f1.SetMaximum(1.08*max_abs_boundary);
            f1.GetXaxis()->SetTitle("x_{B}");
            f1.GetYaxis()->SetTitle("|t boundary| (GeV^{2})");
            f1.GetXaxis()->SetTitleSize(.050);
            f1.GetYaxis()->SetTitleSize(.050);
            f1.GetXaxis()->SetLabelSize(.041);
            f1.GetYaxis()->SetLabelSize(.041);
            f1.GetYaxis()->SetTitleOffset(1.35);
            f1.Draw();

            std::vector<TGraph*> current_graphs;
            std::vector<TGraph*> standard_graphs;
            const int cols[3]={kBlue+1,kGreen+2,kMagenta+1};

            for(std::size_t iq=0;iq<chosen_q2.size() && iq<3;++iq){
                auto* gc=new TGraph();
                auto* gs=new TGraph();
                int nc=0,ns=0;

                std::vector<TDiagPoint> linepts;
                for(const auto&d:diag){
                    if(std::fabs(d.q2-chosen_q2[iq])<1e-10) linepts.push_back(d);
                }
                std::sort(linepts.begin(),linepts.end(),
                    [](const TDiagPoint&a,const TDiagPoint&b){ return a.xb<b.xb; });

                for(const auto&d:linepts){
                    gc->SetPoint(nc++,d.xb,std::fabs(d.current_t));
                    gs->SetPoint(ns++,d.xb,std::fabs(d.standard_t));
                }

                gc->SetMarkerStyle(24);
                gc->SetMarkerSize(1.0);
                gc->SetLineStyle(2);
                gc->SetLineWidth(2);
                gc->SetMarkerColor(cols[iq]);
                gc->SetLineColor(cols[iq]);

                gs->SetMarkerStyle(20);
                gs->SetMarkerSize(1.0);
                gs->SetLineStyle(1);
                gs->SetLineWidth(2);
                gs->SetMarkerColor(cols[iq]);
                gs->SetLineColor(cols[iq]);

                gc->Draw("PL SAME");
                gs->Draw("PL SAME");

                current_graphs.push_back(gc);
                standard_graphs.push_back(gs);
            }

            TLegend leg1(.52,.57,.93,.86);
            leg1.SetBorderSize(0);
            leg1.SetFillStyle(0);
            leg1.SetTextSize(.030);
            for(std::size_t iq=0;iq<chosen_q2.size() && iq<3;++iq){
                std::ostringstream s1,s2;
                s1<<std::fixed<<std::setprecision(2)
                  <<"Current, Q^{2}="<<chosen_q2[iq];
                s2<<std::fixed<<std::setprecision(2)
                  <<"Standard, Q^{2}="<<chosen_q2[iq];
                leg1.AddEntry(current_graphs[iq],s1.str().c_str(),"pl");
                leg1.AddEntry(standard_graphs[iq],s2.str().c_str(),"pl");
            }
            leg1.Draw();

            TLatex l1;
            l1.SetNDC(); l1.SetTextFont(42); l1.SetTextSize(.043);
            l1.DrawLatex(.15,.93,"Boundary expressions at representative Q^{2}");

            // Right panel: the actual quantity that matters for this analysis,
            // the resulting phase-space fraction in every populated cell.
            c.cd(2);
            gPad->SetLeftMargin(.14);
            gPad->SetRightMargin(.04);
            gPad->SetBottomMargin(.14);
            gPad->SetTopMargin(.12);
            gPad->SetGrid();
            gPad->SetTicks(1,1);

            TH1F f2("h_tmin_diag_fraction","",100,0.0,1.02);
            f2.SetMinimum(0.0);
            f2.SetMaximum(1.02);
            f2.GetXaxis()->SetTitle("Allowed fraction: current implementation");
            f2.GetYaxis()->SetTitle("Allowed fraction: standard DVCS t_{min}");
            f2.GetXaxis()->SetTitleSize(.047);
            f2.GetYaxis()->SetTitleSize(.047);
            f2.GetXaxis()->SetLabelSize(.041);
            f2.GetYaxis()->SetLabelSize(.041);
            f2.GetYaxis()->SetTitleOffset(1.25);
            f2.Draw();

            TLine diagline(0.0,0.0,1.0,1.0);
            diagline.SetLineStyle(2);
            diagline.SetLineWidth(2);
            diagline.SetLineColor(kGray+2);
            diagline.Draw();

            TGraph g106,g102;
            int n106=0,n102=0;
            for(const auto&d:diag){
                g106.SetPoint(n106++,d.frac_current_106,d.frac_standard_106);
                g102.SetPoint(n102++,d.frac_current_102,d.frac_standard_102);
            }

            g106.SetMarkerStyle(20);
            g106.SetMarkerSize(.80);
            g106.SetMarkerColor(kBlue+1);
            g106.SetLineColor(kBlue+1);

            g102.SetMarkerStyle(24);
            g102.SetMarkerSize(.80);
            g102.SetMarkerColor(kRed+1);
            g102.SetLineColor(kRed+1);

            g106.Draw("P SAME");
            g102.Draw("P SAME");

            TLegend leg2(.63,.18,.92,.29);
            leg2.SetBorderSize(0);
            leg2.SetFillStyle(0);
            leg2.SetTextSize(.033);
            leg2.AddEntry(&g106,"10.6 GeV","p");
            leg2.AddEntry(&g102,"10.2 GeV","p");
            leg2.Draw();

            TLatex l2;
            l2.SetNDC(); l2.SetTextFont(42); l2.SetTextSize(.043);
            l2.DrawLatex(.14,.93,"Impact on the calculated bin-volume fraction");

            c.cd(0);
            TLatex title;
            title.SetNDC();
            title.SetTextFont(42);
            title.SetTextAlign(22);
            title.SetTextSize(.032);
            title.DrawLatex(.50,.985,
                "Diagnostic comparison of the current t boundary and standard DVCS t_{min}");

            c.SaveAs((note_dir/"tmin_diagnostic_comparison.png").string().c_str());

            for(auto*g:current_graphs) delete g;
            for(auto*g:standard_graphs) delete g;
        }
    }


    auto fractions=[](const std::vector<NoteVolumePoint>& pts){
        std::vector<double> v;
        for(const auto&p:pts) if(std::isfinite(p.fraction)) v.push_back(p.fraction);
        return v;
    };
    const std::vector<double> f106=fractions(bulk106);
    const std::vector<double> f102=fractions(bulk102);

    // Machine-readable global summary (for tables and QA).
    {
        std::ofstream o((note_dir/"bin_volume_summary.csv").string());
        o<<"beam_energy_GeV,bulk_cells,p16_phase_space_fraction,median_phase_space_fraction,p84_phase_space_fraction,min_phase_space_fraction,max_phase_space_fraction\n";
        auto row=[&](double E,const std::vector<double>& v){
            o<<std::fixed<<std::setprecision(6)
             <<E<<","<<v.size()<<","
             <<note_quantile(v,.16)<<","<<note_quantile(v,.50)<<","<<note_quantile(v,.84)<<","
             <<*std::min_element(v.begin(),v.end())<<","<<*std::max_element(v.begin(),v.end())<<"\n";
        };
        row(g10p6.Ebeam,f106);
        row(g10p2.Ebeam,f102);
    }

    // xB-resolved summary CSV (for tables and QA).
    std::vector<std::pair<double,double>> xedges;
    for(const auto&p:bulk106){
        std::pair<double,double> e(p.xbmin,p.xbmax);
        if(std::find(xedges.begin(),xedges.end(),e)==xedges.end()) xedges.push_back(e);
    }
    std::sort(xedges.begin(),xedges.end());

    {
        std::ofstream o((note_dir/"bin_volume_xB_summary.csv").string());
        o<<"beam_energy_GeV,xBmin,xBmax,cells,p16_phase_space_fraction,median_phase_space_fraction,p84_phase_space_fraction\n";
        auto write_energy=[&](const VolumeGroup&G,const std::vector<NoteVolumePoint>& pts){
            for(const auto&e:xedges){
                std::vector<double> v;
                for(const auto&p:pts){
                    if(note_same_edge(p.xbmin,e.first)&&note_same_edge(p.xbmax,e.second)) v.push_back(p.fraction);
                }
                if(v.empty()) continue;
                o<<std::fixed<<std::setprecision(6)
                 <<G.Ebeam<<","<<e.first<<","<<e.second<<","<<v.size()<<","
                 <<note_quantile(v,.16)<<","<<note_quantile(v,.50)<<","<<note_quantile(v,.84)<<"\n";
            }
        };
        write_energy(g10p6,bulk106);
        write_energy(g10p2,bulk102);
    }

    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // 1) Representative xB slices at fixed (Q2,|t|).
    //
    // These slices show how the y, W, and t_min boundaries clip a row of
    // neighboring xB bins.  The allowed fraction is independent of phi for a
    // fixed (xB,Q2,|t|) cell because the same Delta-phi multiplies the nominal
    // and physically allowed four-dimensional volumes.
    {
        struct SliceKey {
            double q0=0,q1=0,t0=0,t1=0;
        };

        auto same_slice=[](const SliceKey&a,const SliceKey&b){
            return note_same_edge(a.q0,b.q0)&&note_same_edge(a.q1,b.q1)&&
                   note_same_edge(a.t0,b.t0)&&note_same_edge(a.t1,b.t1);
        };

        std::vector<SliceKey> keys;
        for(const auto&p:bulk106){
            SliceKey k{p.q2min,p.q2max,p.tmin,p.tmax};
            bool seen=false;
            for(const auto&x:keys){ if(same_slice(k,x)){seen=true;break;} }
            if(!seen) keys.push_back(k);
        }

        auto xbins_for=[&](const std::vector<NoteVolumePoint>&pts,const SliceKey&k){
            std::vector<std::pair<double,double>> xs;
            for(const auto&p:pts){
                if(!note_same_edge(p.q2min,k.q0)||!note_same_edge(p.q2max,k.q1)||
                   !note_same_edge(p.tmin,k.t0)||!note_same_edge(p.tmax,k.t1)) continue;
                std::pair<double,double> e(p.xbmin,p.xbmax);
                if(std::find(xs.begin(),xs.end(),e)==xs.end()) xs.push_back(e);
            }
            std::sort(xs.begin(),xs.end());
            return xs;
        };

        // For each Q2 interval, retain the |t| row with the largest common xB
        // coverage between the two beam energies.  Then select low-, mid-, and
        // high-Q2 examples.  This keeps the figure deterministic and avoids
        // hand-picking an unusually dramatic row.
        std::map<std::pair<double,double>,std::pair<SliceKey,int>> best_by_q;
        for(const auto&k:keys){
            const auto x106=xbins_for(bulk106,k);
            const auto x102=xbins_for(bulk102,k);
            int common=0;
            for(const auto&e:x106) if(std::find(x102.begin(),x102.end(),e)!=x102.end()) ++common;
            auto q=std::make_pair(k.q0,k.q1);
            auto it=best_by_q.find(q);
            if(it==best_by_q.end()||common>it->second.second) best_by_q[q]={k,common};
        }

        std::vector<SliceKey> candidates;
        for(const auto&kv:best_by_q) if(kv.second.second>=3) candidates.push_back(kv.second.first);
        std::sort(candidates.begin(),candidates.end(),[](const SliceKey&a,const SliceKey&b){
            return .5*(a.q0+a.q1)<.5*(b.q0+b.q1);
        });

        std::vector<SliceKey> chosen;
        if(candidates.size()>=3){
            chosen.push_back(candidates.front());
            chosen.push_back(candidates[candidates.size()/2]);
            chosen.push_back(candidates.back());
        } else {
            chosen=candidates;
        }

        if(!chosen.empty()){
            std::ofstream o((note_dir/"bin_volume_representative_xB_slices.csv").string());
            o<<"panel,beam_energy_GeV,Q2min,Q2max,tmin,tmax,xBmin,xBmax,xBcenter,phase_space_fraction\n";
            o<<std::setprecision(10);

            TCanvas c("c_note_binvol_xb_slices","",1500,590);
            c.Divide((int)chosen.size(),1,.002,.002);

            for(std::size_t ip=0;ip<chosen.size();++ip){
                c.cd((int)ip+1);
                gPad->SetLeftMargin(ip==0?.16:.12);
                gPad->SetRightMargin(.035);
                gPad->SetBottomMargin(.16);
                gPad->SetTopMargin(.20);
                gPad->SetGridy();
                gPad->SetTicks(1,1);

                TH1F* frame=new TH1F(Form("h_note_binvol_slice_%zu",ip),"",100,.05,.60);
                frame->SetMinimum(0.0);
                frame->SetMaximum(1.08);
                frame->GetXaxis()->SetTitle("x_{B}");
                frame->GetYaxis()->SetTitle(ip==0?"Physically allowed fraction, V_{bin}/V_{cubic}":"");
                frame->GetXaxis()->SetTitleSize(.055);
                frame->GetYaxis()->SetTitleSize(.050);
                frame->GetXaxis()->SetLabelSize(.045);
                frame->GetYaxis()->SetLabelSize(.043);
                frame->GetYaxis()->SetTitleOffset(ip==0?1.45:1.10);
                frame->Draw();

                TGraphAsymmErrors* g106=new TGraphAsymmErrors();
                TGraphAsymmErrors* g102=new TGraphAsymmErrors();
                g106->SetMarkerStyle(20); g106->SetMarkerSize(1.25);
                g106->SetMarkerColor(kBlue+1); g106->SetLineColor(kBlue+1); g106->SetLineWidth(2);
                g102->SetMarkerStyle(21); g102->SetMarkerSize(1.15);
                g102->SetMarkerColor(kRed+1); g102->SetLineColor(kRed+1); g102->SetLineWidth(2);

                auto fill=[&](TGraphAsymmErrors* g,const std::vector<NoteVolumePoint>&pts,
                              double E,double xoff){
                    std::vector<NoteVolumePoint> row;
                    for(const auto&p:pts){
                        if(!note_same_edge(p.q2min,chosen[ip].q0)||!note_same_edge(p.q2max,chosen[ip].q1)||
                           !note_same_edge(p.tmin,chosen[ip].t0)||!note_same_edge(p.tmax,chosen[ip].t1)) continue;
                        row.push_back(p);
                    }
                    std::sort(row.begin(),row.end(),[](const NoteVolumePoint&a,const NoteVolumePoint&b){
                        return a.xbmin<b.xbmin;
                    });
                    int n=0;
                    for(const auto&p:row){
                        const double xc=.5*(p.xbmin+p.xbmax);
                        const double exl=xc-p.xbmin;
                        const double exh=p.xbmax-xc;
                        g->SetPoint(n,xc+xoff,p.fraction);
                        g->SetPointError(n,exl,exh,0,0);
                        o<<(ip+1)<<","<<E<<","<<chosen[ip].q0<<","<<chosen[ip].q1<<","<<chosen[ip].t0<<","<<chosen[ip].t1<<","<<p.xbmin<<","<<p.xbmax<<","<<xc<<","<<p.fraction<<"\n";
                        ++n;
                    }
                };

                fill(g106,bulk106,g10p6.Ebeam,-.0015);
                fill(g102,bulk102,g10p2.Ebeam,+.0015);
                g106->Draw("PL SAME");
                g102->Draw("PL SAME");

                TLine* one=new TLine(.05,1.0,.60,1.0);
                one->SetLineStyle(2); one->SetLineWidth(2); one->SetLineColor(kGray+2); one->Draw();

                TLatex lab;
                lab.SetNDC(); lab.SetTextFont(42); lab.SetTextSize(.038);
                lab.DrawLatex(ip==0?.17:.13,.88,
                    Form("%.3f < Q^{2} < %.3f GeV^{2}",chosen[ip].q0,chosen[ip].q1));
                lab.DrawLatex(ip==0?.17:.13,.82,
                    Form("%.3f < |t| < %.3f GeV^{2}",chosen[ip].t0,chosen[ip].t1));

                if(ip==chosen.size()-1){
                    TLegend* leg=new TLegend(.55,.61,.92,.76);
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(.037);
                    leg->AddEntry(g106,"10.6 GeV","pl");
                    leg->AddEntry(g102,"10.2 GeV","pl");
                    leg->Draw();
                }
            }

            c.cd(0);
            TLatex title;
            title.SetNDC(); title.SetTextFont(42); title.SetTextAlign(22); title.SetTextSize(.034);
            title.DrawLatex(.50,.975,"Representative phase-space clipping across rows of x_{B} bins");

            c.SaveAs((note_dir/"bin_volume_representative_xB_slices.png").string().c_str());
        }
    }


    // 2) Representative phi dependence of the raw four-dimensional volume.
    //
    // The physical mask itself depends only on (xB,Q2,|t|), but the analysis
    // uses unequal phi-bin widths.  Therefore both the nominal rectangular
    // volume and the physically allowed volume vary with phi through Delta-phi.
    // Their ratio remains constant for a fixed (xB,Q2,|t|) cell.
    {
        struct PhiCell {
            double x0=0,x1=0,q0=0,q1=0,t0=0,t1=0;
        };

        auto same_cell=[](const NoteVolumePoint&p,const PhiCell&c){
            return note_same_edge(p.xbmin,c.x0)&&note_same_edge(p.xbmax,c.x1)&&
                   note_same_edge(p.q2min,c.q0)&&note_same_edge(p.q2max,c.q1)&&
                   note_same_edge(p.tmin,c.t0)&&note_same_edge(p.tmax,c.t1);
        };

        std::vector<PhiCell> cells;
        for(const auto&p:p106){
            PhiCell c{p.xbmin,p.xbmax,p.q2min,p.q2max,p.tmin,p.tmax};
            bool seen=false;
            for(const auto&x:cells){
                if(note_same_edge(c.x0,x.x0)&&note_same_edge(c.x1,x.x1)&&
                   note_same_edge(c.q0,x.q0)&&note_same_edge(c.q1,x.q1)&&
                   note_same_edge(c.t0,x.t0)&&note_same_edge(c.t1,x.t1)){
                    seen=true; break;
                }
            }
            if(!seen) cells.push_back(c);
        }

        // Prefer a common cell with full phi coverage and moderate clipping.
        // A target allowed fraction near 0.7 makes the nominal and allowed
        // curves visibly distinct without choosing an extreme corner bin.
        const double target_fraction=0.70;
        int best_common_phi=-1;
        double best_score=1e9;
        PhiCell best;
        bool have_best=false;

        for(const auto&c:cells){
            std::vector<double> f106_cell,f102_cell;
            std::set<std::pair<double,double>> phi106,phi102;

            for(const auto&p:p106){
                if(!same_cell(p,c) || !std::isfinite(p.fraction)) continue;
                f106_cell.push_back(p.fraction);
                phi106.insert({p.phimin,p.phimax});
            }
            for(const auto&p:p102){
                if(!same_cell(p,c) || !std::isfinite(p.fraction)) continue;
                f102_cell.push_back(p.fraction);
                phi102.insert({p.phimin,p.phimax});
            }
            if(f106_cell.empty()||f102_cell.empty()) continue;

            int common_phi=0;
            for(const auto&e:phi106) if(phi102.count(e)) ++common_phi;
            if(common_phi<4) continue;

            const double med106=note_quantile(f106_cell,.50);
            const double med102=note_quantile(f102_cell,.50);
            const double score=std::fabs(.5*(med106+med102)-target_fraction);

            if(common_phi>best_common_phi ||
               (common_phi==best_common_phi && score<best_score)){
                best_common_phi=common_phi;
                best_score=score;
                best=c;
                have_best=true;
            }
        }

        if(have_best){
            struct PhiPoint {
                double phi=0, cubic=0, allowed106=0, allowed102=0;
                double frac106=0, frac102=0;
            };
            std::vector<PhiPoint> pts;

            for(const auto&p:p106){
                if(!same_cell(p,best) || p.cubic<=0 || p.allowed<=0) continue;

                const NoteVolumePoint* match102=nullptr;
                for(const auto&q:p102){
                    if(!same_cell(q,best)) continue;
                    if(note_same_edge(q.phimin,p.phimin)&&note_same_edge(q.phimax,p.phimax)){
                        match102=&q; break;
                    }
                }
                if(!match102 || match102->allowed<=0) continue;

                PhiPoint x;
                x.phi=.5*(p.phimin+p.phimax);
                x.cubic=p.cubic;
                x.allowed106=p.allowed;
                x.allowed102=match102->allowed;
                x.frac106=p.allowed/p.cubic;
                x.frac102=match102->allowed/match102->cubic;
                pts.push_back(x);
            }

            std::sort(pts.begin(),pts.end(),[](const PhiPoint&a,const PhiPoint&b){
                return a.phi<b.phi;
            });

            if(!pts.empty()){
                double vmax=0.0;
                for(const auto&p:pts) vmax=std::max(vmax,p.cubic);

                // Choose a display scale automatically so the ordinate is O(1).
                double scale=1.0;
                int exponent=0;
                if(vmax>0){
                    exponent=(int)std::floor(std::log10(vmax));
                    scale=std::pow(10.0,-exponent);
                }

                TCanvas c("c_note_binvol_phi","",1150,720);
                c.SetLeftMargin(.135);
                c.SetRightMargin(.035);
                c.SetBottomMargin(.14);
                c.SetTopMargin(.095);
                c.SetGridy();
                c.SetTicks(1,1);

                TH1F frame("h_note_binvol_phi","",100,0,360);
                frame.SetMinimum(0.0);
                frame.SetMaximum(1.18*vmax*scale);
                frame.GetXaxis()->SetTitle("#phi (deg)");

                std::ostringstream ytitle;
                ytitle<<"Four-dimensional bin volume";
                if(exponent!=0) ytitle<<" (#times10^{"<<exponent<<"})";
                frame.GetYaxis()->SetTitle(ytitle.str().c_str());
                frame.GetXaxis()->SetTitleSize(.050);
                frame.GetYaxis()->SetTitleSize(.050);
                frame.GetXaxis()->SetLabelSize(.041);
                frame.GetYaxis()->SetLabelSize(.041);
                frame.GetYaxis()->SetTitleOffset(1.30);
                frame.Draw();

                TGraph gcubic,g106,g102;
                gcubic.SetMarkerStyle(24);
                gcubic.SetMarkerSize(1.20);
                gcubic.SetMarkerColor(kGray+2);
                gcubic.SetLineColor(kGray+2);
                gcubic.SetLineWidth(2);

                g106.SetMarkerStyle(20);
                g106.SetMarkerSize(1.20);
                g106.SetMarkerColor(kBlue+1);
                g106.SetLineColor(kBlue+1);
                g106.SetLineWidth(2);

                g102.SetMarkerStyle(21);
                g102.SetMarkerSize(1.10);
                g102.SetMarkerColor(kRed+1);
                g102.SetLineColor(kRed+1);
                g102.SetLineWidth(2);

                for(std::size_t i=0;i<pts.size();++i){
                    gcubic.SetPoint((int)i,pts[i].phi,pts[i].cubic*scale);
                    g106.SetPoint((int)i,pts[i].phi,pts[i].allowed106*scale);
                    g102.SetPoint((int)i,pts[i].phi,pts[i].allowed102*scale);
                }

                gcubic.Draw("PL SAME");
                g106.Draw("PL SAME");
                g102.Draw("PL SAME");

                TLegend leg(.62,.69,.93,.86);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextSize(.033);
                leg.AddEntry(&gcubic,"Nominal rectangular volume","pl");
                leg.AddEntry(&g106,"Allowed volume, 10.6 GeV","pl");
                leg.AddEntry(&g102,"Allowed volume, 10.2 GeV","pl");
                leg.Draw();

                TLatex title;
                title.SetNDC();
                title.SetTextFont(42);
                title.SetTextSize(.047);
                title.DrawLatex(.135,.925,
                    "Representative #phi dependence of the four-dimensional bin volume");

                TLatex note;
                note.SetNDC();
                note.SetTextFont(42);
                note.SetTextSize(.027);
                note.DrawLatex(.145,.845,
                    Form("%.3f < x_{B} < %.3f,  %.3f < Q^{2} < %.3f GeV^{2},  %.3f < |t| < %.3f GeV^{2}",
                         best.x0,best.x1,best.q0,best.q1,best.t0,best.t1));

                std::ofstream o((note_dir/"bin_volume_phi_example.csv").string());
                o<<"xBmin,xBmax,Q2min,Q2max,tmin,tmax,phi,cubic_volume,allowed_volume_10p6,allowed_fraction_10p6,allowed_volume_10p2,allowed_fraction_10p2\n";
                o<<std::setprecision(10);
                for(const auto&p:pts){
                    o<<best.x0<<","<<best.x1<<","
                     <<best.q0<<","<<best.q1<<","
                     <<best.t0<<","<<best.t1<<","
                     <<p.phi<<","<<p.cubic<<","
                     <<p.allowed106<<","<<p.frac106<<","
                     <<p.allowed102<<","<<p.frac102<<"\n";
                }

                c.SaveAs((note_dir/"bin_volume_phi_example.png").string().c_str());
            }
        }
    }

    // 3) Representative xB-Q2 map for a fixed |t| interval.  Choose the
    // interval with the largest number of common cells in the 10.6 GeV table.
    std::map<std::pair<double,double>,int> tcounts;
    for(const auto&p:bulk106) ++tcounts[{p.tmin,p.tmax}];
    std::pair<double,double> best_t=tcounts.begin()->first;
    int best_n=-1;
    for(const auto&kv:tcounts){
        if(kv.second>best_n){best_n=kv.second; best_t=kv.first;}
    }

    auto make_map=[&](const std::vector<NoteVolumePoint>&pts,
                      const std::string& hname,
                      const std::string& title_text,
                      const std::string& outname)
    {
        std::vector<double> xbounds,qbounds;
        for(const auto&p:pts){
            if(!note_same_edge(p.tmin,best_t.first)||!note_same_edge(p.tmax,best_t.second)) continue;
            xbounds.push_back(p.xbmin); xbounds.push_back(p.xbmax);
            qbounds.push_back(p.q2min); qbounds.push_back(p.q2max);
        }
        std::sort(xbounds.begin(),xbounds.end());
        xbounds.erase(std::unique(xbounds.begin(),xbounds.end(),
                     [](double a,double b){return note_same_edge(a,b);}),xbounds.end());
        std::sort(qbounds.begin(),qbounds.end());
        qbounds.erase(std::unique(qbounds.begin(),qbounds.end(),
                     [](double a,double b){return note_same_edge(a,b);}),qbounds.end());
        if(xbounds.size()<2||qbounds.size()<2) return;

        TH2D h(hname.c_str(),"",
               int(xbounds.size()-1),xbounds.data(),
               int(qbounds.size()-1),qbounds.data());

        for(const auto&p:pts){
            if(!note_same_edge(p.tmin,best_t.first)||!note_same_edge(p.tmax,best_t.second)) continue;
            const double xc=.5*(p.xbmin+p.xbmax);
            const double qc=.5*(p.q2min+p.q2max);
            h.SetBinContent(h.FindBin(xc,qc),p.fraction);
        }

        TCanvas c((hname+"_canvas").c_str(),"",1050,760);
        c.SetLeftMargin(.13); c.SetRightMargin(.14);
        c.SetBottomMargin(.13); c.SetTopMargin(.10);
        c.SetTicks(1,1);

        h.SetMinimum(0.0); h.SetMaximum(1.0);
        h.GetXaxis()->SetTitle("x_{B}");
        h.GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
        h.GetZaxis()->SetTitle("V_{bin}/V_{cubic}");
        h.GetXaxis()->SetTitleSize(.048); h.GetYaxis()->SetTitleSize(.048);
        h.GetZaxis()->SetTitleSize(.043);
        h.GetXaxis()->SetLabelSize(.039); h.GetYaxis()->SetLabelSize(.039);
        h.GetZaxis()->SetLabelSize(.036);
        h.GetYaxis()->SetTitleOffset(1.15);
        h.GetZaxis()->SetTitleOffset(1.20);
        h.Draw("COLZ TEXT");

        TLatex title;
        title.SetNDC(); title.SetTextFont(42); title.SetTextSize(.045);
        title.DrawLatex(.13,.935,title_text.c_str());

        TLatex note;
        note.SetNDC(); note.SetTextFont(42); note.SetTextSize(.030);
        note.DrawLatex(.14,.865,
            Form("%.3f < |t| < %.3f GeV^{2}",best_t.first,best_t.second));

        c.SaveAs((note_dir/outname).string().c_str());
    };

    make_map(bulk106,"h_note_binvol_map_106",
             "Physical phase-space coverage of the nominal bins at 10.6 GeV",
             "bin_volume_phase_space_map_10p6.png");
    make_map(bulk102,"h_note_binvol_map_102",
             "Physical phase-space coverage of the nominal bins at 10.2 GeV",
             "bin_volume_phase_space_map_10p2.png");

    std::cout<<"[binvol-note] Wrote analysis-note outputs to "<<note_dir<<"\n";
}


} // end anonymous namespace

// =====================================================================
// Public driver
// =====================================================================

bool update_bin_volume_csv(const std::string& csv_path,
                           const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const std::string csv_abs = fs::absolute(csv_path).string();
    std::error_code ec;
    const uintmax_t size_before =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[binvol] CSV: " << csv_abs
              << " (size=" << size_before << ")\n";

    CsvDoc csv;
    if (!csv.load(csv_path)) {
        std::cerr << "[binvol] ERROR: failed to load CSV.\n";
        return false;
    }

    const int NR = csv.nrows();
    std::cout << "[binvol] Loaded CSV with " << NR << " data rows.\n";

    // Define groups: 10.6 GeV (4 periods) and 10.2 GeV (Sp19 Inb only)
    VolumeGroup g10p6;
    g10p6.label      = "10.6 GeV";
    g10p6.xbavg_col  = "xBavg, 10.6 GeV";
    g10p6.phiavg_col = "phiavg, 10.6 GeV";
    g10p6.q2avg_col  = "Q2avg, 10.6 GeV";
    g10p6.tabavg_col = "t_abs_avg, 10.6 GeV";
    g10p6.binvol_col       = "bin_volume, 10.6 GeV";
    g10p6.cubic_binvol_col = "cubic bin_volume, 10.6 GeV";
    g10p6.energy_dir       = "10.60";
    g10p6.Ebeam      = 10.60;

    VolumeGroup g10p2;
    g10p2.label      = "10.2 GeV";
    g10p2.xbavg_col  = "xBavg, Sp19 Inb";
    g10p2.phiavg_col = "phiavg, Sp19 Inb";
    g10p2.q2avg_col  = "Q2avg, Sp19 Inb";
    g10p2.tabavg_col = "t_abs_avg, Sp19 Inb";
    g10p2.binvol_col       = "bin_volume, 10.2 GeV";
    g10p2.cubic_binvol_col = "cubic bin_volume, 10.2 GeV";
    g10p2.energy_dir       = "10.2";
    g10p2.Ebeam      = 10.20;

    // Compute bin volumes for all phi-binned rows for both energies.
    // Do not gate this on xBavg/data availability: cross_sections.cpp needs
    // deterministic phase-space volumes wherever a yield later appears.
    std::vector<bool> row_has_data_10p6(NR, true);
    std::vector<bool> row_has_data_10p2(NR, true);

    // Compute volumes (CSV-only; no trees needed).
    compute_bin_volumes_for_group(g10p6, csv, row_has_data_10p6);
    compute_bin_volumes_for_group(g10p2, csv, row_has_data_10p2);

    // Draw canvases (sequential, no threads).
    draw_bin_volume_canvases(g10p6, csv, row_has_data_10p6, out_root_dir);
    draw_bin_volume_canvases(g10p2, csv, row_has_data_10p2, out_root_dir);

    // Compact, note-quality diagnostics.
    write_bin_volume_analysis_note_outputs(csv, g10p6, g10p2, out_root_dir);

    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[binvol] ERROR: failed to save updated CSV.\n";
        return false;
    }

    const uintmax_t size_after =
        fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[binvol] Updated CSV: " << csv_abs
              << " (size " << size_before << " -> " << size_after << ")\n";
    std::cout << "[binvol] Bin-volume computation complete.\n";

    return true;
}