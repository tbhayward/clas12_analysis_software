// radiative_corrections.cpp
//
// Compute radiative correction factors Frad from generator MC and write
// them into dvcs_pass2_analysis.csv as three-tuples "(value, stat, sys)"
// in the columns:
//
//   "Frad, 10.6 GeV"
//   "Frad, 10.2 GeV"
//
// using the Lee CSV binning (xBmin/xBmax, Q2min/Q2max, t_abs_min/t_abs_max, phimin/phimax).
//
// Beam energies:
//   10.6 GeV: DVCS_Sp18_inb, DVCS_Sp18_out, DVCS_Fa18_inb, DVCS_Fa18_out
//   10.2 GeV: DVCS_Sp19_inb
//
// The 10.6 GeV factors are filled only when "xBavg, 10.6 GeV" is non-empty.
// The 10.2 GeV factors are filled only when "xBavg, Sp19 Inb" is non-empty.

#include "radiative_corrections.h"

#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TH1.h>
#include <TLatex.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

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

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

// -------------------------- small helpers --------------------------

static inline std::string trim(const std::string& s) {
    size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    size_t j = s.size();
    while (j > i && std::isspace(static_cast<unsigned char>(s[j - 1]))) --j;
    return s.substr(i, j - i);
}

static inline bool is_empty(const std::string& s) {
    return trim(s).empty();
}

static inline int find_column(const std::vector<std::string>& header,
                              const std::string& name) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) return static_cast<int>(i);
    }
    return -1;
}

static inline double parse_double_or_nan(const std::string& s) {
    std::string t = trim(s);
    if (t.empty()) return std::numeric_limits<double>::quiet_NaN();
    try {
        return std::stod(t);
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

static inline int range_index(const std::pair<double,double>& r,
                              std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i) {
        if (ranges[i].first == r.first && ranges[i].second == r.second) return i;
    }
    ranges.push_back(r);
    return static_cast<int>(ranges.size()) - 1;
}

static inline int findBin1D(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int k = 0; k < static_cast<int>(ranges.size()); ++k) {
        const double lo = ranges[k].first;
        const double hi = ranges[k].second;
        if (v >= lo && v < hi) return k;
    }
    // allow right edge to land in the last bin
    if (!ranges.empty()) {
        const int last = static_cast<int>(ranges.size()) - 1;
        if (v == ranges[last].second) return last;
    }
    return -1;
}

// ------------------------ CSV read / write ------------------------

struct CsvTable {
    std::vector<std::string> header;
    std::vector< std::vector<std::string> > rows;
};

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string field;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (in_quotes) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    field.push_back('"');
                    ++i;
                } else {
                    in_quotes = false;
                }
            } else {
                field.push_back(c);
            }
        } else {
            if (c == ',') {
                out.push_back(field);
                field.clear();
            } else if (c == '"') {
                in_quotes = true;
            } else {
                field.push_back(c);
            }
        }
    }
    out.push_back(field);
    return out;
}

static bool read_csv(const std::string& path, CsvTable& table) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[radcorr] ERROR: cannot open CSV \"" << path << "\" for reading.\n";
        return false;
    }
    std::string line;
    if (!std::getline(ifs, line)) {
        std::cerr << "[radcorr] ERROR: CSV \"" << path << "\" appears to be empty.\n";
        return false;
    }
    table.header = split_csv_line(line);

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        table.rows.push_back(split_csv_line(line));
    }

    // Normalize row width
    const size_t ncols = table.header.size();
    for (auto& r : table.rows) {
        if (r.size() < ncols) r.resize(ncols);
        else if (r.size() > ncols) r.resize(ncols);
    }
    return true;
}

static std::string csv_escape(const std::string& field) {
    bool need_quotes = false;
    for (char c : field) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            need_quotes = true;
            break;
        }
    }
    if (!need_quotes) return field;
    std::string out;
    out.push_back('"');
    for (char c : field) {
        if (c == '"') out.push_back('"');
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

static bool write_csv_atomic(const std::string& path, const CsvTable& table) {
    namespace fs = std::filesystem;
    fs::path p(path);
    fs::path tmp = p;
    tmp += ".tmp";

    {
        std::ofstream ofs(tmp.string());
        if (!ofs) {
            std::cerr << "[radcorr] ERROR: cannot open temp CSV \"" << tmp.string()
                      << "\" for writing.\n";
            return false;
        }
        // header
        for (size_t i = 0; i < table.header.size(); ++i) {
            if (i) ofs << ",";
            ofs << csv_escape(table.header[i]);
        }
        ofs << "\n";
        // rows
        for (const auto& row : table.rows) {
            for (size_t i = 0; i < row.size(); ++i) {
                if (i) ofs << ",";
                ofs << csv_escape(row[i]);
            }
            ofs << "\n";
        }
    }

    std::error_code ec;
    fs::rename(tmp, p, ec);
    if (ec) {
        std::cerr << "[radcorr] ERROR: failed to replace \"" << path
                  << "\" with temp file \"" << tmp.string() << "\": " << ec.message() << "\n";
        return false;
    }
    return true;
}

// ------------------- MC accumulation and RC -------------------

struct BranchesGen {
    double x = 0.0;
    double Q2 = 0.0;
    double t1 = 0.0;
    double phi2 = 0.0;
    bool has_x   = false;
    bool has_Q2  = false;
    bool has_t1  = false;
    bool has_phi = false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindD = [&](const char* name, double* addr, bool& flag) {
            if (t->GetBranch(name)) {
                t->SetBranchAddress(name, addr);
                flag = true;
            }
        };
        bindD("x",   &x,   has_x);
        bindD("Q2",  &Q2,  has_Q2);
        bindD("t1",  &t1,  has_t1);
        bindD("phi2",&phi2,has_phi);
    }

    bool ready() const {
        return has_x && has_Q2 && has_t1 && has_phi;
    }
};

struct CountsStore {
    std::map< std::tuple<int,int,int,int>, double > born_phi;
    std::map< std::tuple<int,int,int,int>, double > rad_phi;
    std::map< std::tuple<int,int,int>,     double > born_tot;
    std::map< std::tuple<int,int,int>,     double > rad_tot;
};

static void accumulate_generated(
    TTree* t,
    bool isBorn,
    const std::vector<std::pair<double,double>>& xB_ranges,
    const std::vector<std::pair<double,double>>& Q2_ranges,
    const std::vector<std::pair<double,double>>& t_ranges,
    CountsStore& acc)
{
    if (!t) return;
    BranchesGen b;
    b.bind(t);
    if (!b.ready()) {
        std::cerr << "[radcorr] WARNING: generated tree is missing one or more branches "
                  << "(x, Q2, t1, phi2); this tree will be skipped.\n";
        return;
    }

    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);

        const double xB = b.x;
        const double Q2 = b.Q2;
        const double tt = std::fabs(b.t1);

        // phi2 is in radians
        double phi = b.phi2;
        double w   = std::fmod(phi, TWO_PI);
        if (w < 0.0) w += TWO_PI;
        const double width = TWO_PI / static_cast<double>(N_PHI_BINS);
        int iphi = static_cast<int>(std::floor(w / width));
        if (iphi < 0) iphi = 0;
        if (iphi >= N_PHI_BINS) iphi = N_PHI_BINS - 1;

        int ix = findBin1D(xB, xB_ranges);
        int iQ = findBin1D(Q2, Q2_ranges);
        int it = findBin1D(tt,  t_ranges);
        if (ix < 0 || iQ < 0 || it < 0) continue;

        auto key3 = std::make_tuple(ix, iQ, it);
        auto key4 = std::make_tuple(ix, iQ, it, iphi);

        if (isBorn) {
            acc.born_phi[key4] += 1.0;
            acc.born_tot[key3] += 1.0;
        } else {
            acc.rad_phi[key4]  += 1.0;
            acc.rad_tot[key3]  += 1.0;
        }
    }
}

static std::vector<double> phi_centers_deg() {
    std::vector<double> v(N_PHI_BINS);
    const double step = 360.0 / static_cast<double>(N_PHI_BINS);
    for (int i = 0; i < N_PHI_BINS; ++i) {
        v[i] = (i + 0.5) * step;
    }
    return v;
}

struct PhiArrays {
    std::vector<double> phi_deg;
    std::vector<double> rc;
    std::vector<double> rc_err;
    double A_born = 0.0;
    double B_rad  = 0.0;
};

using RCPerCell = std::map< std::tuple<int,int,int>, PhiArrays >;

static RCPerCell compute_rc_per_cell(
    const CountsStore& acc,
    const std::vector<std::pair<double,double>>& xB_ranges,
    const std::vector<std::pair<double,double>>& Q2_ranges,
    const std::vector<std::pair<double,double>>& t_ranges)
{
    RCPerCell out;
    const auto PHI_DEG = phi_centers_deg();

    for (int ix = 0; ix < static_cast<int>(xB_ranges.size()); ++ix) {
        for (int iQ = 0; iQ < static_cast<int>(Q2_ranges.size()); ++iQ) {
            for (int it = 0; it < static_cast<int>(t_ranges.size()); ++it) {
                auto key3 = std::make_tuple(ix, iQ, it);
                double A = 0.0;
                double B = 0.0;

                auto itA = acc.born_tot.find(key3);
                if (itA != acc.born_tot.end()) A = itA->second;
                auto itB = acc.rad_tot.find(key3);
                if (itB != acc.rad_tot.end()) B = itB->second;

                PhiArrays pa;
                pa.phi_deg = PHI_DEG;
                pa.rc.assign(N_PHI_BINS, 1.0);
                pa.rc_err.assign(N_PHI_BINS, 0.0);
                pa.A_born = A;
                pa.B_rad  = B;

                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    auto key4 = std::make_tuple(ix, iQ, it, ip);
                    double a = 0.0;
                    double b = 0.0;
                    auto ia = acc.born_phi.find(key4);
                    if (ia != acc.born_phi.end()) a = ia->second;
                    auto ib = acc.rad_phi.find(key4);
                    if (ib != acc.rad_phi.end()) b = ib->second;

                    double RC  = 1.0;
                    double sRC = 0.0;

                    if (A > 0.0 && B > 0.0 && a > 0.0 && b > 0.0) {
                        RC = (a * B) / (b * A);
                        double aSafe = std::max(a, 1.0);
                        double ASafe = std::max(A, 1.0);
                        double bSafe = std::max(b, 1.0);
                        double BSafe = std::max(B, 1.0);
                        sRC = RC * std::sqrt(1.0 / aSafe + 1.0 / ASafe
                                           + 1.0 / bSafe + 1.0 / BSafe);
                    } else if (A > 0.0 && B > 0.0) {
                        RC  = 1.0;
                        sRC = 0.0;
                    }

                    if (!std::isfinite(RC)) RC = 1.0;
                    RC = std::max(0.0, std::min(RC, 2.0));

                    pa.rc[ip]     = RC;
                    pa.rc_err[ip] = sRC;
                }

                out[key3] = pa;
            }
        }
    }

    return out;
}

// ---------------------- plotting helpers ----------------------

struct RowBinKey {
    int ix   = -1;
    int iQ   = -1;
    int it   = -1;
    int iphi = -1;
};

static void plot_group_rc(
    const std::string& energy_label_for_title,   // "10.60" or "10.2"
    const std::vector<std::pair<double,double>>& xB_ranges,
    const std::vector<std::pair<double,double>>& Q2_ranges,
    const std::vector<std::pair<double,double>>& t_ranges,
    const std::vector<RowBinKey>& row_keys,
    const RCPerCell& rcPerCell,
    const std::string& out_dir_plots)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(out_dir_plots, ec);

    // For each xB range, determine which Q2 and t indices are actually used
    const int nx = static_cast<int>(xB_ranges.size());

    for (int ix = 0; ix < nx; ++ix) {
        std::set<int> usedQ;
        std::set<int> usedT;
        for (const auto& rk : row_keys) {
            if (rk.ix == ix) {
                if (rk.iQ >= 0) usedQ.insert(rk.iQ);
                if (rk.it >= 0) usedT.insert(rk.it);
            }
        }
        if (usedQ.empty() || usedT.empty()) continue;

        std::vector<int> qIdx(usedQ.begin(), usedQ.end());
        std::vector<int> tIdx(usedT.begin(), usedT.end());

        // sort Q2 and t by their low edge
        std::sort(qIdx.begin(), qIdx.end(),
                  [&](int a, int b){ return Q2_ranges[a].first < Q2_ranges[b].first; });
        std::sort(tIdx.begin(), tIdx.end(),
                  [&](int a, int b){ return t_ranges[a].first < t_ranges[b].first; });

        const int ncols = static_cast<int>(qIdx.size());
        const int nrows = static_cast<int>(tIdx.size());

        const int W = 280 * ncols + 160;
        const int H = 240 * nrows + 170;

        std::ostringstream cname;
        cname << "c_rc_E" << energy_label_for_title << "_xB" << ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        TPad* pTop  = new TPad("pTop", "pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0);
        pTop->SetBorderSize(0);
        pTop->Draw();

        TPad* pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0);
        pGrid->SetBorderSize(0);
        pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        TLatex head;
        head.SetNDC();
        head.SetTextAlign(22);
        head.SetTextFont(42);
        head.SetTextSize(0.36);
        std::ostringstream tit;
        tit << "Radiative correction  E_{beam}=" << energy_label_for_title
            << " GeV   x_{B} #in [" << std::setprecision(2) << std::fixed
            << xB_ranges[ix].first << ", " << xB_ranges[ix].second << "]";
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        for (int r = 0; r < nrows; ++r) {
            const int it = tIdx[r];
            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int iQ = qIdx[ccol];
                pGrid->cd(r * ncols + ccol + 1);

                gPad->SetGrid(1, 1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 2.0);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetTitle("#phi (deg)");
                ay->SetTitle("Born/Rad");
                ax->CenterTitle();
                ay->CenterTitle();

                ax->SetNdivisions(505);
                ax->SetTitleSize(0.060);
                ay->SetTitleSize(0.060);
                ay->SetLabelSize(0.048);
                ax->SetTitleOffset(1.25);
                ay->SetTitleOffset(1.35);

                auto itCell = rcPerCell.find(std::make_tuple(ix, iQ, it));
                if (itCell == rcPerCell.end()) continue;
                const auto& pa = itCell->second;

                int npt = static_cast<int>(pa.phi_deg.size());
                std::vector<double> ex(npt, 0.0);
                TGraphErrors* gr = new TGraphErrors(
                    npt,
                    const_cast<double*>(pa.phi_deg.data()),
                    const_cast<double*>(pa.rc.data()),
                    ex.data(),
                    const_cast<double*>(pa.rc_err.data())
                );
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(1.0);
                gr->SetLineWidth(1);
                gr->Draw("P SAME");

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.040);
                lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   |t| #in [%.2g, %.2g]",
                        Q2_ranges[iQ].first, Q2_ranges[iQ].second,
                        t_ranges[it].first,  t_ranges[it].second));
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/rc_E" << energy_label_for_title << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ---------------------- main driver ----------------------

struct EnergyConfig {
    std::string energy_label_column;   // used only for messages
    std::string xBavg_col_name;        // gate column (non-empty => fill Frad)
    std::string frad_col_name;         // output column "(val,stat,sys)"
    std::string plot_dir_suffix;       // "10.60" or "10.2"
    std::vector<std::string> period_tags; // DVCS_Sp18_inb, etc.
};

} // anon

bool update_radiative_corrections_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& radGenMcTrees,
    const std::string& out_root_dir)
{
    CsvTable table;
    if (!read_csv(csv_path, table)) {
        return false;
    }

    // Required columns
    std::vector<std::string> missing;

    auto require_col = [&](const std::string& name) -> int {
        int idx = find_column(table.header, name);
        if (idx < 0) missing.push_back(name);
        return idx;
    };

    int c_xb_min   = require_col("xBmin");
    int c_xb_max   = require_col("xBmax");
    int c_q2_min   = require_col("Q2min");
    int c_q2_max   = require_col("Q2max");
    int c_t_min    = require_col("t_abs_min");
    int c_t_max    = require_col("t_abs_max");
    int c_phi_min  = require_col("phimin");
    int c_phi_max  = require_col("phimax");

    int c_xbavg_106 = require_col("xBavg, 10.6 GeV");
    int c_xbavg_sp19 = require_col("xBavg, Sp19 Inb");

    int c_frad_106 = require_col("Frad, 10.6 GeV");
    int c_frad_102 = require_col("Frad, 10.2 GeV");

    if (!missing.empty()) {
        std::cerr << "[radcorr][FATAL] CSV \"" << csv_path
                  << "\" is missing required columns:\n";
        for (const auto& name : missing) {
            std::cerr << "  - " << name << "\n";
        }
        return false;
    }

    const size_t nrows = table.rows.size();
    if (nrows == 0) {
        std::cerr << "[radcorr][FATAL] CSV \"" << csv_path
                  << "\" has no data rows.\n";
        return false;
    }

    // Build bin ranges and row->(ix,iQ,it,iphi) mapping
    std::vector<std::pair<double,double>> xB_ranges;
    std::vector<std::pair<double,double>> Q2_ranges;
    std::vector<std::pair<double,double>> t_ranges;
    std::vector<RowBinKey> row_keys(nrows);

    for (size_t i = 0; i < nrows; ++i) {
        const auto& row = table.rows[i];
        double xbmin  = parse_double_or_nan(row[c_xb_min]);
        double xbmax  = parse_double_or_nan(row[c_xb_max]);
        double q2min  = parse_double_or_nan(row[c_q2_min]);
        double q2max  = parse_double_or_nan(row[c_q2_max]);
        double tmin   = parse_double_or_nan(row[c_t_min]);
        double tmax   = parse_double_or_nan(row[c_t_max]);
        double phimin = parse_double_or_nan(row[c_phi_min]);
        double phimax = parse_double_or_nan(row[c_phi_max]);

        if (!std::isfinite(xbmin) || !std::isfinite(xbmax) ||
            !std::isfinite(q2min) || !std::isfinite(q2max) ||
            !std::isfinite(tmin)  || !std::isfinite(tmax)  ||
            !std::isfinite(phimin)|| !std::isfinite(phimax)) {
            row_keys[i] = RowBinKey();
            continue;
        }

        int ix = range_index(std::make_pair(xbmin, xbmax), xB_ranges);
        int iQ = range_index(std::make_pair(q2min, q2max), Q2_ranges);
        int it = range_index(std::make_pair(tmin,  tmax),  t_ranges);

        // map phimin/phimax to phi bin index
        const double center = 0.5 * (phimin + phimax);
        const double step   = 360.0 / static_cast<double>(N_PHI_BINS);
        int iphi = static_cast<int>(std::floor(center / step));
        if (iphi < 0) iphi = 0;
        if (iphi >= N_PHI_BINS) iphi = N_PHI_BINS - 1;

        row_keys[i].ix   = ix;
        row_keys[i].iQ   = iQ;
        row_keys[i].it   = it;
        row_keys[i].iphi = iphi;
    }

    // Define energy groups
    std::vector<EnergyConfig> energies;
    energies.push_back(EnergyConfig{
        "10.6 GeV",
        "xBavg, 10.6 GeV",
        "Frad, 10.6 GeV",
        "10.60",
        {"DVCS_Sp18_inb", "DVCS_Sp18_out", "DVCS_Fa18_inb", "DVCS_Fa18_out"}
    });
    energies.push_back(EnergyConfig{
        "10.2 GeV",
        "xBavg, Sp19 Inb",
        "Frad, 10.2 GeV",
        "10.2",
        {"DVCS_Sp19_inb"}
    });

    // Precompute which rows are active for each energy
    std::map<std::string, std::vector<bool> > row_active;
    for (const auto& E : energies) {
        int idx_xbavg = find_column(table.header, E.xBavg_col_name);
        if (idx_xbavg < 0) {
            std::cerr << "[radcorr][FATAL] Missing xBavg gating column \""
                      << E.xBavg_col_name << "\" for " << E.energy_label_column << ".\n";
            return false;
        }
        std::vector<bool> active(nrows, false);
        for (size_t i = 0; i < nrows; ++i) {
            const std::string val = table.rows[i][idx_xbavg];
            if (!is_empty(val)) active[i] = true;
        }
        row_active[E.energy_label_column] = std::move(active);
    }

    // For each energy, accumulate MC, compute RC per cell, and then fill CSV
    namespace fs = std::filesystem;
    fs::path plot_root = fs::path(out_root_dir) / "radiative_correction_plots";

    // We keep the RCPerCell per energy so we can also make plots.
    std::map<std::string, RCPerCell> rc_by_energy;

    for (const auto& E : energies) {
        CountsStore acc;

        // Accumulate over all periods in this energy group
        for (const auto& tag : E.period_tags) {
            const std::string bornKey = tag + "_gen";
            const std::string radKey  = tag + "_gen_rad";

            auto itB = genMcTrees.find(bornKey);
            auto itR = radGenMcTrees.find(radKey);

            if (itB == genMcTrees.end() || itR == radGenMcTrees.end()) {
                std::cerr << "[radcorr][FATAL] Missing generator MC trees for period \""
                          << tag << "\" in energy group " << E.energy_label_column
                          << ". Expected keys \"" << bornKey << "\" and \""
                          << radKey << "\".\n";
                return false;
            }

            TTree* tBorn = itB->second;
            TTree* tRad  = itR->second;

            accumulate_generated(tBorn, true,
                                 xB_ranges, Q2_ranges, t_ranges,
                                 acc);
            accumulate_generated(tRad, false,
                                 xB_ranges, Q2_ranges, t_ranges,
                                 acc);
        }

        RCPerCell rc = compute_rc_per_cell(acc, xB_ranges, Q2_ranges, t_ranges);
        rc_by_energy[E.energy_label_column] = rc;

        // Diagnostic plots
        {
            std::error_code ec;
            fs::create_directories(plot_root / E.plot_dir_suffix, ec);
            plot_group_rc(E.plot_dir_suffix,
                          xB_ranges, Q2_ranges, t_ranges,
                          row_keys,
                          rc,
                          (plot_root / E.plot_dir_suffix).string());
        }

        std::cout << "[radcorr] Computed radiative corrections for "
                  << E.energy_label_column << ".\n";
    }

    // Now fill the CSV columns with "(value, stat, 0.0)"
    for (const auto& E : energies) {
        RCPerCell& rc = rc_by_energy[E.energy_label_column];

        int idx_frad = find_column(table.header, E.frad_col_name);
        if (idx_frad < 0) {
            std::cerr << "[radcorr][FATAL] Missing output column \""
                      << E.frad_col_name << "\" while filling "
                      << E.energy_label_column << " radiative corrections.\n";
            return false;
        }

        const std::vector<bool>& active = row_active[E.energy_label_column];

        for (size_t i = 0; i < nrows; ++i) {
            if (!active[i]) continue;

            const RowBinKey& rk = row_keys[i];
            if (rk.ix < 0 || rk.iQ < 0 || rk.it < 0 || rk.iphi < 0) {
                std::cerr << "[radcorr][FATAL] Row " << i
                          << " has invalid bin indices while trying to fill "
                          << E.frad_col_name << ".\n";
                return false;
            }

            auto itCell = rc.find(std::make_tuple(rk.ix, rk.iQ, rk.it));
            if (itCell == rc.end()) {
                std::cerr << "[radcorr][FATAL] No RC cell found for row " << i
                          << " in " << E.energy_label_column
                          << " (ix=" << rk.ix << ", iQ=" << rk.iQ
                          << ", it=" << rk.it << ").\n";
                return false;
            }
            const PhiArrays& pa = itCell->second;
            if (rk.iphi < 0 || rk.iphi >= static_cast<int>(pa.rc.size()) ||
                rk.iphi >= static_cast<int>(pa.rc_err.size())) {
                std::cerr << "[radcorr][FATAL] Invalid phi index " << rk.iphi
                          << " for row " << i << " in " << E.energy_label_column << ".\n";
                return false;
            }

            double val  = pa.rc[rk.iphi];
            double stat = pa.rc_err[rk.iphi];
            if (!std::isfinite(val))  val  = 1.0;
            if (!std::isfinite(stat)) stat = 0.0;

            std::ostringstream buf;
            buf << "("
                << std::setprecision(6) << std::fixed << val
                << ", "
                << std::setprecision(6) << std::fixed << stat
                << ", 0.0)";
            table.rows[i][idx_frad] = buf.str();
        }
    }

    if (!write_csv_atomic(csv_path, table)) {
        return false;
    }

    std::cout << "[radcorr] Updated radiative correction columns in \"" << csv_path << "\".\n";
    return true;
}