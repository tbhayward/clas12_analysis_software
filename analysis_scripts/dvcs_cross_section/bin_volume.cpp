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
// Also produces bin-volume vs phi canvases per beam energy and xB bin under
//   output/bin_volume/10.60
//   output/bin_volume/10.2
//
// All valid phi-binned rows in the pass-2 CSV are filled for both beam
// energies. Phi-integrated rows (phi width ~360 deg) are skipped.

#include "bin_volume.h"

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
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

static double calculate_phase_space_allowed_bin_volume(double xB_min, double xB_max,
                                                       double Q2_min, double Q2_max,
                                                       double t_abs_min, double t_abs_max,   // positive |t| edges
                                                       double phi_min, double phi_max,       // radians
                                                       double E_beam)
{
    constexpr int    n_steps = 10;
    constexpr double Mp      = 0.938272; // Proton mass (GeV)
    int valid_count = 0;

    // Convert |t| edges to physical t<0 for the loop
    const double t_phys_min = -t_abs_max;
    const double t_phys_max = -t_abs_min;

    const double dxB = (xB_max - xB_min) / n_steps;
    const double dQ2 = (Q2_max - Q2_min) / n_steps;
    const double dt  = (t_phys_max - t_phys_min) / n_steps;

    for (int i = 0; i < n_steps; ++i) {
        const double xB = xB_min + (i + 0.5) * dxB;

        for (int j = 0; j < n_steps; ++j) {
            const double Q2 = Q2_min + (j + 0.5) * dQ2;

            // DVCS t_min(xB,Q2) (negative)
            const double sqrt_term = std::sqrt(1.0 + (4.0 * Mp * Mp * xB * xB) / Q2);
            const double t_min_val = -Q2 * (1.0 - xB) * (1.0 - xB) / (xB * (1.0 + sqrt_term));

            for (int k = 0; k < n_steps; ++k) {
                const double t = t_phys_min + (k + 0.5) * dt;

                // y and W cuts
                const double y  = Q2 / (2.0 * Mp * xB * E_beam);
                const double W2 = Mp * Mp + Q2 * (1.0 / xB - 1.0);
                const double W  = (W2 > 0.0) ? std::sqrt(W2) : 0.0;

                if (t > t_min_val && y > 0.19 && y < 0.8113919276 && W > 2.0) {
                    ++valid_count;
                }
            }
        }
    }

    const double fraction = (double)valid_count / (double)(n_steps * n_steps * n_steps);
    const double cubic_volume = calculate_cubic_bin_volume(
        xB_min, xB_max,
        Q2_min, Q2_max,
        t_abs_min, t_abs_max,
        phi_min, phi_max
    );

    return cubic_volume * fraction;
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

        const double V_allowed = calculate_phase_space_allowed_bin_volume(
            xbmin, xbmax,
            q2min, q2max,
            tmin,  tmax,
            phi_min_rad, phi_max_rad,
            G.Ebeam
        );

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