// norm_scenario_comparison.cpp
// -----------------------------------------------------------------------------
// Standalone comparison plotter for three normalized DVCS cross-section CSVs.
//
// Use case:
//   1) CSV A: eppi0 normalization = 1, BH normalization = 1
//   2) CSV B: eppi0 normalization floats, BH normalization = 1
//   3) CSV C: eppi0 normalization = 1, BH normalization floats
//
// This program reads all three CSVs at runtime and makes:
//   - normed cross-section comparison plots with all three data scenarios and
//     BH/KM15/VGG predictions overlaid;
//   - BH-only ratio plots, data / BH, for the same three scenarios.
//
// It does not modify any CSV. It only reads the existing columns:
//   normed cross sections, ep->epg, exp, <label>, <helicity>
// and the kinematic-bin columns xBmin/xBmax, Q2min/Q2max, t_abs_min/t_abs_max,
// phimin/phimax, and phiavg, <label>.
//
// Place in:
//   /u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/
//
// Compile through the Makefile target added below.
// -----------------------------------------------------------------------------

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TH1.h>
#include <TSystem.h>
#include <TLine.h>

#include <nlohmann/json.hpp>

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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace {

using Range = std::pair<double, double>;
using QTKey = std::pair<Range, Range>;

struct Triple {
    double value = 0.0;
    double stat = 0.0;
    double sys = 0.0;
};

struct Point {
    double phi = 0.0;
    double xs = 0.0;
    double err = 0.0;
};

struct BinData {
    std::vector<Point> points[3];
    size_t theory_row = 0;
    bool have_theory_row = false;
};

struct XBGroup {
    Range xb = {0.0, 0.0};
    std::map<QTKey, BinData> bins;
};

struct CsvData {
    std::string path;
    std::string scenario_name;
    std::map<Range, XBGroup> xb_groups;
};

struct TheoryCurves {
    std::vector<double> phi_deg;
    std::vector<double> bh[3];
    std::vector<double> km[3];
    std::vector<double> vgg[3];
};

static const std::vector<std::string> kHelNames = {"unpol", "pos", "neg"};
static const std::vector<std::string> kHelTitles = {"unpolarized", "+ helicity", "- helicity"};
static const std::vector<std::string> kDefaultLabels = {
    "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out",
    "Fa18", "Sp18", "10.6 GeV", "10.2 GeV"
};

static std::string trim(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;

    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;

    return s.substr(b, e - b);
}

static std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> out;
    std::string field;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

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
            if (c == '"') {
                in_quotes = true;
            } else if (c == ',') {
                out.push_back(field);
                field.clear();
            } else {
                field.push_back(c);
            }
        }
    }

    out.push_back(field);
    return out;
}

static int find_col_optional(const std::vector<std::string> &header,
                             const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(header[i]) == target) return (int)i;
    }

    return -1;
}

static int find_col_required(const std::vector<std::string> &header,
                             const std::string &target) {
    const int idx = find_col_optional(header, target);
    if (idx < 0) throw std::runtime_error("Missing required column: " + target);
    return idx;
}

static double parse_double_or_nan(const std::string &cell) {
    const std::string s = trim(cell);
    if (s.empty()) return std::numeric_limits<double>::quiet_NaN();

    char *endptr = nullptr;
    const double v = std::strtod(s.c_str(), &endptr);
    if (endptr == s.c_str()) return std::numeric_limits<double>::quiet_NaN();
    return v;
}

static Triple parse_tuple3(const std::string &cell_raw) {
    Triple out;
    std::string s = trim(cell_raw);
    if (s.empty()) return out;

    if (!s.empty() && s.front() == '(' && s.back() == ')') {
        s = trim(s.substr(1, s.size() - 2));
    }

    if (s.empty()) return out;

    std::vector<std::string> parts;
    std::string token;

    for (char c : s) {
        if (c == ',') {
            parts.push_back(trim(token));
            token.clear();
        } else {
            token.push_back(c);
        }
    }
    parts.push_back(trim(token));

    auto to_double = [](const std::string &x) -> double {
        if (x.empty()) return 0.0;
        return std::atof(x.c_str());
    };

    if (parts.size() > 0U) out.value = to_double(parts[0]);
    if (parts.size() > 1U) out.stat = to_double(parts[1]);
    if (parts.size() > 2U) out.sys = to_double(parts[2]);
    return out;
}

static bool good_positive(double x) {
    return std::isfinite(x) && x > 0.0;
}

static std::string canonical_dir(const std::string &label) {
    if (label == "Fa18 Inb") return "Fa18_Inb";
    if (label == "Fa18 Out") return "Fa18_Out";
    if (label == "Fa18 Inb Supp") return "Fa18_Inb_Supp";
    if (label == "Sp18 Inb") return "Sp18_Inb";
    if (label == "Sp18 Out") return "Sp18_Out";
    if (label == "Sp19 Inb") return "Sp19_Inb";
    if (label == "Fa18") return "Fa18";
    if (label == "Sp18") return "Sp18";
    if (label == "10.6 GeV") return "10.6_GeV";
    if (label == "10.2 GeV") return "10.2_GeV";

    std::string out = label;
    std::replace(out.begin(), out.end(), ' ', '_');
    return out;
}

static std::string theory_energy_label_for(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") return "10.2 GeV";
    return "10.6 GeV";
}

static void ensure_dir(const fs::path &p) {
    std::error_code ec;
    fs::create_directories(p, ec);
    if (!fs::exists(p, ec)) {
        throw std::runtime_error("Could not create output directory: " + p.string());
    }
}

static bool column_has_nonzero_data(const std::vector<std::vector<std::string>> &rows,
                                    int col) {
    if (col < 0) return false;

    for (const auto &row : rows) {
        if ((size_t)col >= row.size()) continue;
        const Triple t = parse_tuple3(row[(size_t)col]);
        if (good_positive(t.value)) return true;
    }

    return false;
}

static CsvData read_csv_for_label(const std::string &path,
                                  const std::string &scenario_name,
                                  const std::string &label,
                                  std::vector<bool> *helicity_available) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Could not open CSV: " + path);

    std::string line;
    if (!std::getline(fin, line)) throw std::runtime_error("CSV is empty: " + path);

    const std::vector<std::string> header = split_csv_line(line);

    const int c_xb_min = find_col_required(header, "xBmin");
    const int c_xb_max = find_col_required(header, "xBmax");
    const int c_q2_min = find_col_required(header, "Q2min");
    const int c_q2_max = find_col_required(header, "Q2max");
    const int c_t_min = find_col_required(header, "t_abs_min");
    const int c_t_max = find_col_required(header, "t_abs_max");
    const int c_phi_min = find_col_required(header, "phimin");
    const int c_phi_max = find_col_required(header, "phimax");

    int c_phi_avg = find_col_optional(header, "phiavg, " + label);
    if (c_phi_avg < 0) c_phi_avg = -1;

    int c_xs[3] = {-1, -1, -1};
    for (int ih = 0; ih < 3; ++ih) {
        c_xs[ih] = find_col_optional(
            header,
            "normed cross sections, ep->epg, exp, " + label + ", " + kHelNames[(size_t)ih]
        );
    }

    std::vector<std::vector<std::string>> rows;
    rows.reserve(8192);
    while (std::getline(fin, line)) {
        if (trim(line).empty()) continue;
        rows.push_back(split_csv_line(line));
    }

    if (helicity_available) {
        helicity_available->assign(3, false);
        for (int ih = 0; ih < 3; ++ih) {
            (*helicity_available)[(size_t)ih] = column_has_nonzero_data(rows, c_xs[ih]);
        }
    }

    CsvData out;
    out.path = path;
    out.scenario_name = scenario_name;

    for (size_t irow = 0; irow < rows.size(); ++irow) {
        const auto &row = rows[irow];

        auto get = [&](int col) -> std::string {
            if (col < 0 || (size_t)col >= row.size()) return "";
            return row[(size_t)col];
        };

        const double xb_min = parse_double_or_nan(get(c_xb_min));
        const double xb_max = parse_double_or_nan(get(c_xb_max));
        const double q2_min = parse_double_or_nan(get(c_q2_min));
        const double q2_max = parse_double_or_nan(get(c_q2_max));
        const double t_min = parse_double_or_nan(get(c_t_min));
        const double t_max = parse_double_or_nan(get(c_t_max));

        if (!std::isfinite(xb_min) || !std::isfinite(xb_max) ||
            !std::isfinite(q2_min) || !std::isfinite(q2_max) ||
            !std::isfinite(t_min) || !std::isfinite(t_max)) {
            continue;
        }

        double phi = std::numeric_limits<double>::quiet_NaN();
        if (c_phi_avg >= 0) phi = parse_double_or_nan(get(c_phi_avg));

        if (!std::isfinite(phi)) {
            const double phi_min = parse_double_or_nan(get(c_phi_min));
            const double phi_max = parse_double_or_nan(get(c_phi_max));
            if (!std::isfinite(phi_min) || !std::isfinite(phi_max)) continue;
            phi = 0.5 * (phi_min + phi_max);
        }

        Range xb = {xb_min, xb_max};
        Range q2 = {q2_min, q2_max};
        Range tt = {t_min, t_max};
        QTKey key = {q2, tt};

        XBGroup &g = out.xb_groups[xb];
        g.xb = xb;

        BinData &bd = g.bins[key];
        bd.theory_row = irow;
        bd.have_theory_row = true;

        for (int ih = 0; ih < 3; ++ih) {
            if (c_xs[ih] < 0) continue;
            const Triple xs = parse_tuple3(get(c_xs[ih]));
            if (!good_positive(xs.value)) continue;

            Point p;
            p.phi = phi;
            p.xs = xs.value;
            p.err = std::max(0.0, xs.stat);
            bd.points[ih].push_back(p);
        }
    }

    for (auto &xbit : out.xb_groups) {
        for (auto &bit : xbit.second.bins) {
            for (int ih = 0; ih < 3; ++ih) {
                std::sort(bit.second.points[ih].begin(), bit.second.points[ih].end(),
                          [](const Point &a, const Point &b) {
                              return a.phi < b.phi;
                          });
            }
        }
    }

    return out;
}

static std::vector<Range> collect_q2_ranges(const std::vector<CsvData> &csvs,
                                            const Range &xb) {
    std::set<Range> s;
    for (const auto &csv : csvs) {
        auto itx = csv.xb_groups.find(xb);
        if (itx == csv.xb_groups.end()) continue;
        for (const auto &kv : itx->second.bins) s.insert(kv.first.first);
    }
    return std::vector<Range>(s.begin(), s.end());
}

static std::vector<Range> collect_t_ranges(const std::vector<CsvData> &csvs,
                                           const Range &xb) {
    std::set<Range> s;
    for (const auto &csv : csvs) {
        auto itx = csv.xb_groups.find(xb);
        if (itx == csv.xb_groups.end()) continue;
        for (const auto &kv : itx->second.bins) s.insert(kv.first.second);
    }
    return std::vector<Range>(s.begin(), s.end());
}

static std::vector<Range> collect_xb_ranges(const std::vector<CsvData> &csvs) {
    std::set<Range> s;
    for (const auto &csv : csvs) {
        for (const auto &kv : csv.xb_groups) s.insert(kv.first);
    }
    return std::vector<Range>(s.begin(), s.end());
}

static std::map<size_t, TheoryCurves> load_theory_for_label(const std::string &label,
                                                            const std::string &theory_root) {
    std::map<size_t, TheoryCurves> out;

    const std::string energy_label = theory_energy_label_for(label);
    const fs::path file = fs::path(theory_root) / canonical_dir(energy_label) / "xs_phi_all.json";

    if (!fs::exists(file)) {
        std::cerr << "[norm_scenario_comparison] WARNING: no theory JSON at "
                  << file.string() << "\n";
        return out;
    }

    std::ifstream fin(file);
    if (!fin) {
        std::cerr << "[norm_scenario_comparison] WARNING: could not open theory JSON at "
                  << file.string() << "\n";
        return out;
    }

    json j;
    fin >> j;

    const std::vector<double> phi_deg = j.value("phi_deg", std::vector<double>{});
    if (phi_deg.empty() || !j.contains("rows") || !j["rows"].is_object()) {
        std::cerr << "[norm_scenario_comparison] WARNING: malformed theory JSON at "
                  << file.string() << "\n";
        return out;
    }

    for (auto it = j["rows"].begin(); it != j["rows"].end(); ++it) {
        size_t row_index = 0;
        try {
            row_index = (size_t)std::stoul(it.key());
        } catch (...) {
            continue;
        }

        const json &cell = it.value();
        TheoryCurves tc;
        tc.phi_deg = phi_deg;

        tc.bh[0] = cell["BH"].value("unpol", std::vector<double>{});
        tc.bh[1] = cell["BH"].value("pos", std::vector<double>{});
        tc.bh[2] = cell["BH"].value("neg", std::vector<double>{});

        tc.km[0] = cell["KM"].value("unpol", std::vector<double>{});
        tc.km[1] = cell["KM"].value("pos", std::vector<double>{});
        tc.km[2] = cell["KM"].value("neg", std::vector<double>{});

        tc.vgg[0] = cell["VGG"].value("unpol", std::vector<double>{});
        tc.vgg[1] = cell["VGG"].value("pos", std::vector<double>{});
        tc.vgg[2] = cell["VGG"].value("neg", std::vector<double>{});

        out[row_index] = std::move(tc);
    }

    std::cout << "[norm_scenario_comparison] Loaded theory label=" << label
              << " energy=" << energy_label << " rows=" << out.size()
              << " file=" << file.string() << "\n";

    return out;
}

static double interpolate_curve(const std::vector<double> &xs,
                                const std::vector<double> &ys,
                                double x) {
    if (xs.empty() || xs.size() != ys.size()) return std::numeric_limits<double>::quiet_NaN();
    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();

    auto hi = std::lower_bound(xs.begin(), xs.end(), x);
    if (hi == xs.end()) return ys.back();
    const size_t ihi = (size_t)std::distance(xs.begin(), hi);
    if (ihi == 0U) return ys.front();
    const size_t ilo = ihi - 1U;

    const double x0 = xs[ilo];
    const double x1 = xs[ihi];
    const double y0 = ys[ilo];
    const double y1 = ys[ihi];
    if (x1 == x0) return y0;

    const double f = (x - x0) / (x1 - x0);
    return y0 + f * (y1 - y0);
}

static TGraphErrors *make_point_graph(const std::vector<Point> &points,
                                       int color,
                                       int marker_style) {
    TGraphErrors *g = new TGraphErrors();
    g->SetMarkerStyle(marker_style);
    g->SetMarkerSize(0.8);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);

    int ip = 0;
    for (const auto &p : points) {
        if (!good_positive(p.xs)) continue;
        g->SetPoint(ip, p.phi, p.xs);
        g->SetPointError(ip, 0.0, std::max(0.0, p.err));
        ++ip;
    }

    return g;
}

static TGraph *make_curve_graph(const std::vector<double> &phi,
                                const std::vector<double> &ys,
                                int color,
                                int style,
                                int width) {
    TGraph *g = new TGraph();
    g->SetLineColor(color);
    g->SetLineStyle(style);
    g->SetLineWidth(width);

    if (phi.size() != ys.size()) return g;

    int ip = 0;
    for (size_t i = 0; i < phi.size(); ++i) {
        if (!good_positive(ys[i])) continue;
        g->SetPoint(ip, phi[i], ys[i]);
        ++ip;
    }

    return g;
}

static TGraphErrors *make_bh_ratio_graph(const std::vector<Point> &points,
                                         const TheoryCurves &tc,
                                         int ih,
                                         int color,
                                         int marker_style) {
    TGraphErrors *g = new TGraphErrors();
    g->SetMarkerStyle(marker_style);
    g->SetMarkerSize(0.8);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);

    if (tc.phi_deg.size() != tc.bh[ih].size()) return g;

    int ip = 0;
    for (const auto &p : points) {
        if (!good_positive(p.xs)) continue;
        const double model = interpolate_curve(tc.phi_deg, tc.bh[ih], p.phi);
        if (!good_positive(model)) continue;

        g->SetPoint(ip, p.phi, p.xs / model);
        g->SetPointError(ip, 0.0, std::max(0.0, p.err / model));
        ++ip;
    }

    return g;
}

static void draw_empty_axes(double ymin,
                            double ymax,
                            const std::string &ytitle,
                            bool logy) {
    TH1D *frame = new TH1D("frame", "", 1, 0.0, 360.0);
    frame->SetDirectory(nullptr);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->GetXaxis()->SetTitle("#phi (deg)");
    frame->GetYaxis()->SetTitle(ytitle.c_str());
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetLabelSize(0.060);
    frame->GetYaxis()->SetLabelSize(0.060);
    frame->GetXaxis()->SetTitleSize(0.070);
    frame->GetYaxis()->SetTitleSize(0.070);
    frame->GetXaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetTitleOffset(1.00);
    frame->Draw("AXIS");

    if (logy) gPad->SetLogy(true);
}

static void draw_bin_title(const Range &q2, const Range &tt) {
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.060);
    lat.SetTextAlign(12);

    std::ostringstream oss;
    oss << "Q^{2}=(" << std::fixed << std::setprecision(2)
        << q2.first << ", " << q2.second << ")  |t|=("
        << tt.first << ", " << tt.second << ")";

    lat.DrawLatex(0.12, 0.90, oss.str().c_str());
}

static std::pair<double, double> cross_section_yrange(const std::vector<const BinData*> &bins,
                                                      const std::vector<const TheoryCurves*> &theories,
                                                      int ih) {
    double ymin = std::numeric_limits<double>::max();
    double ymax = 0.0;

    for (const BinData *bd : bins) {
        if (!bd) continue;
        for (const auto &p : bd->points[ih]) {
            if (!good_positive(p.xs)) continue;
            ymin = std::min(ymin, std::max(1e-30, p.xs - p.err));
            ymax = std::max(ymax, p.xs + p.err);
        }
    }

    for (const TheoryCurves *tc : theories) {
        if (!tc) continue;
        const std::vector<const std::vector<double>*> curves = {&tc->bh[ih], &tc->km[ih], &tc->vgg[ih]};
        for (const auto *curve : curves) {
            for (double y : *curve) {
                if (!good_positive(y)) continue;
                ymin = std::min(ymin, y);
                ymax = std::max(ymax, y);
            }
        }
    }

    if (!good_positive(ymax) || ymin == std::numeric_limits<double>::max()) return {1e-4, 1.0};

    ymin *= 0.50;
    ymax *= 2.00;
    if (!good_positive(ymin)) ymin = 1e-4;
    if (ymax <= ymin) ymax = ymin * 10.0;
    return {ymin, ymax};
}

static void make_cross_section_canvas(const std::string &label,
                                      const Range &xb,
                                      int xb_index,
                                      int ih,
                                      const std::vector<CsvData> &csvs,
                                      const std::map<size_t, TheoryCurves> &theory,
                                      const fs::path &outdir) {
    const std::vector<Range> q2_ranges = collect_q2_ranges(csvs, xb);
    const std::vector<Range> t_ranges = collect_t_ranges(csvs, xb);
    if (q2_ranges.empty() || t_ranges.empty()) return;

    const int ncols = (int)q2_ranges.size();
    const int nrows = (int)t_ranges.size();
    const int W = 300 * ncols + 160;
    const int H = 260 * nrows + 240;

    TCanvas *c = new TCanvas("c_norm_scenario_xs", "c_norm_scenario_xs", W, H);

    TPad *pTop = new TPad("pTop", "pTop", 0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad *pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.78);
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
    head.SetTextSize(0.060);

    std::ostringstream title;
    title << "Normed cross sections comparison, ep #rightarrow ep#gamma   "
          << label << "   x_{B} in (" << std::fixed << std::setprecision(3)
          << xb.first << ", " << xb.second << ")   (" << kHelTitles[(size_t)ih] << ")";
    head.DrawLatex(0.5, 0.86, title.str().c_str());

    const int scenario_colors[3] = {kBlack, kRed + 1, kBlue + 1};
    const int scenario_markers[3] = {20, 24, 25};

    TGraphErrors dummy_s0;
    TGraphErrors dummy_s1;
    TGraphErrors dummy_s2;
    TGraph dummy_bh;
    TGraph dummy_km;
    TGraph dummy_vgg;

    TGraphErrors *dummy_s[3] = {&dummy_s0, &dummy_s1, &dummy_s2};
    for (int is = 0; is < 3; ++is) {
        dummy_s[is]->SetMarkerStyle(scenario_markers[is]);
        dummy_s[is]->SetMarkerColor(scenario_colors[is]);
        dummy_s[is]->SetLineColor(scenario_colors[is]);
        dummy_s[is]->SetLineWidth(1);
    }

    dummy_bh.SetLineColor(kGreen + 2);
    dummy_bh.SetLineStyle(2);
    dummy_bh.SetLineWidth(2);
    dummy_km.SetLineColor(kMagenta + 1);
    dummy_km.SetLineStyle(1);
    dummy_km.SetLineWidth(2);
    dummy_vgg.SetLineColor(kOrange + 7);
    dummy_vgg.SetLineStyle(1);
    dummy_vgg.SetLineWidth(2);

    TLegend *leg = new TLegend(0.12, 0.12, 0.88, 0.72);
    leg->SetNColumns(3);
    leg->SetBorderSize(1);
    leg->SetLineColor(kBlack);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045);
    leg->AddEntry(&dummy_s0, csvs[0].scenario_name.c_str(), "lep");
    leg->AddEntry(&dummy_s1, csvs[1].scenario_name.c_str(), "lep");
    leg->AddEntry(&dummy_s2, csvs[2].scenario_name.c_str(), "lep");
    leg->AddEntry(&dummy_bh, "BH", "l");
    leg->AddEntry(&dummy_km, "KM15", "l");
    leg->AddEntry(&dummy_vgg, "VGG", "l");
    leg->Draw();

    for (int ir = 0; ir < nrows; ++ir) {
        for (int ic = 0; ic < ncols; ++ic) {
            const Range &tt = t_ranges[(size_t)ir];
            const Range &q2 = q2_ranges[(size_t)ic];
            const QTKey key = {q2, tt};

            pGrid->cd(ir * ncols + ic + 1);
            gPad->SetGrid(1, 1);
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.07);
            gPad->SetLogy(true);

            std::vector<const BinData*> bin_ptrs;
            std::vector<const TheoryCurves*> tc_ptrs;
            const TheoryCurves *first_tc = nullptr;

            for (const auto &csv : csvs) {
                const BinData *bd = nullptr;
                auto itx = csv.xb_groups.find(xb);
                if (itx != csv.xb_groups.end()) {
                    auto itb = itx->second.bins.find(key);
                    if (itb != itx->second.bins.end()) bd = &(itb->second);
                }
                bin_ptrs.push_back(bd);

                if (bd && bd->have_theory_row) {
                    auto itt = theory.find(bd->theory_row);
                    if (itt != theory.end()) {
                        tc_ptrs.push_back(&(itt->second));
                        if (!first_tc) first_tc = &(itt->second);
                    }
                }
            }

            const auto yr = cross_section_yrange(bin_ptrs, tc_ptrs, ih);
            draw_empty_axes(yr.first, yr.second,
                            "d^{4}#sigma_{norm} / (dx_{B} dQ^{2} d|t| d#phi)",
                            true);
            draw_bin_title(q2, tt);

            if (first_tc) {
                TGraph *gbh = make_curve_graph(first_tc->phi_deg, first_tc->bh[ih], kGreen + 2, 2, 2);
                TGraph *gkm = make_curve_graph(first_tc->phi_deg, first_tc->km[ih], kMagenta + 1, 1, 2);
                TGraph *gvgg = make_curve_graph(first_tc->phi_deg, first_tc->vgg[ih], kOrange + 7, 1, 2);
                gbh->Draw("L SAME");
                gkm->Draw("L SAME");
                gvgg->Draw("L SAME");
            }

            for (int is = 0; is < 3; ++is) {
                if (!bin_ptrs[(size_t)is]) continue;
                TGraphErrors *gp = make_point_graph(bin_ptrs[(size_t)is]->points[ih],
                                                    scenario_colors[is],
                                                    scenario_markers[is]);
                gp->Draw("PE1 SAME");
            }
        }
    }

    ensure_dir(outdir);
    const fs::path out = outdir / ("norm_scenario_compare_" + kHelNames[(size_t)ih] + "_" +
                                  canonical_dir(label) + "_xB_" + std::to_string(xb_index) + ".png");
    c->SaveAs(out.string().c_str());
    delete c;
}

static void make_bh_ratio_canvas(const std::string &label,
                                 const Range &xb,
                                 int xb_index,
                                 int ih,
                                 const std::vector<CsvData> &csvs,
                                 const std::map<size_t, TheoryCurves> &theory,
                                 const fs::path &outdir) {
    const std::vector<Range> q2_ranges = collect_q2_ranges(csvs, xb);
    const std::vector<Range> t_ranges = collect_t_ranges(csvs, xb);
    if (q2_ranges.empty() || t_ranges.empty()) return;

    const int ncols = (int)q2_ranges.size();
    const int nrows = (int)t_ranges.size();
    const int W = 300 * ncols + 160;
    const int H = 260 * nrows + 240;

    TCanvas *c = new TCanvas("c_norm_scenario_bh_ratio", "c_norm_scenario_bh_ratio", W, H);

    TPad *pTop = new TPad("pTop", "pTop", 0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad *pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.78);
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
    head.SetTextSize(0.060);

    std::ostringstream title;
    title << "Normed cross-section / BH ratios, ep #rightarrow ep#gamma   "
          << label << "   x_{B} in (" << std::fixed << std::setprecision(3)
          << xb.first << ", " << xb.second << ")   (" << kHelTitles[(size_t)ih] << ")";
    head.DrawLatex(0.5, 0.86, title.str().c_str());

    const int scenario_colors[3] = {kBlack, kRed + 1, kBlue + 1};
    const int scenario_markers[3] = {20, 24, 25};

    TGraphErrors dummy_s0;
    TGraphErrors dummy_s1;
    TGraphErrors dummy_s2;
    TGraphErrors *dummy_s[3] = {&dummy_s0, &dummy_s1, &dummy_s2};
    for (int is = 0; is < 3; ++is) {
        dummy_s[is]->SetMarkerStyle(scenario_markers[is]);
        dummy_s[is]->SetMarkerColor(scenario_colors[is]);
        dummy_s[is]->SetLineColor(scenario_colors[is]);
        dummy_s[is]->SetLineWidth(1);
    }

    TLegend *leg = new TLegend(0.18, 0.18, 0.82, 0.72);
    leg->SetNColumns(1);
    leg->SetBorderSize(1);
    leg->SetLineColor(kBlack);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045);
    leg->AddEntry(&dummy_s0, csvs[0].scenario_name.c_str(), "lep");
    leg->AddEntry(&dummy_s1, csvs[1].scenario_name.c_str(), "lep");
    leg->AddEntry(&dummy_s2, csvs[2].scenario_name.c_str(), "lep");
    leg->Draw();

    for (int ir = 0; ir < nrows; ++ir) {
        for (int ic = 0; ic < ncols; ++ic) {
            const Range &tt = t_ranges[(size_t)ir];
            const Range &q2 = q2_ranges[(size_t)ic];
            const QTKey key = {q2, tt};

            pGrid->cd(ir * ncols + ic + 1);
            gPad->SetGrid(1, 1);
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.07);
            gPad->SetLogy(false);

            draw_empty_axes(0.0, 2.0, "data / BH", false);
            draw_bin_title(q2, tt);

            TLine line_one(0.0, 1.0, 360.0, 1.0);
            line_one.SetLineColor(kBlack);
            line_one.SetLineStyle(1);
            line_one.SetLineWidth(1);
            line_one.Draw("SAME");

            for (int is = 0; is < 3; ++is) {
                const BinData *bd = nullptr;
                auto itx = csvs[(size_t)is].xb_groups.find(xb);
                if (itx != csvs[(size_t)is].xb_groups.end()) {
                    auto itb = itx->second.bins.find(key);
                    if (itb != itx->second.bins.end()) bd = &(itb->second);
                }

                if (!bd || !bd->have_theory_row) continue;

                auto itt = theory.find(bd->theory_row);
                if (itt == theory.end()) continue;

                TGraphErrors *gr = make_bh_ratio_graph(bd->points[ih], itt->second, ih,
                                                       scenario_colors[is],
                                                       scenario_markers[is]);
                gr->Draw("PE1 SAME");
            }
        }
    }

    ensure_dir(outdir);
    const fs::path out = outdir / ("norm_scenario_bh_ratio_" + kHelNames[(size_t)ih] + "_" +
                                  canonical_dir(label) + "_xB_" + std::to_string(xb_index) + ".png");
    c->SaveAs(out.string().c_str());
    delete c;
}

static std::vector<std::string> parse_label_list(int argc,
                                                 char **argv,
                                                 int first_label_arg) {
    if (argc <= first_label_arg) return kDefaultLabels;

    std::vector<std::string> labels;
    for (int i = first_label_arg; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--label" || arg == "--labels") continue;

        size_t start = 0;
        while (start <= arg.size()) {
            const size_t comma = arg.find(',', start);
            const std::string token = trim(arg.substr(start, comma == std::string::npos ? std::string::npos : comma - start));
            if (!token.empty()) labels.push_back(token);
            if (comma == std::string::npos) break;
            start = comma + 1U;
        }
    }

    if (labels.empty()) return kDefaultLabels;
    return labels;
}

static void print_usage(const char *prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " <csv_unity> <csv_eppi0_float> <csv_bh_float> [theory_json_root] [out_root] [labels...]\n\n"
        << "Defaults:\n"
        << "  theory_json_root = output/jsons/cross_sections\n"
        << "  out_root         = output/norm_scenario_comparison\n\n"
        << "Examples:\n"
        << "  " << prog << " output/csvs/dvcs_pass2_unity.csv output/csvs/dvcs_pass2_eppi0_float.csv output/csvs/dvcs_pass2_bh_float.csv\n\n"
        << "  " << prog << " csv1.csv csv2.csv csv3.csv output/jsons/cross_sections output/norm_scenario_comparison \"10.6 GeV\"\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    try {
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetLegendFont(42);
        TH1::AddDirectory(kFALSE);

        const std::string csv_unity = argv[1];
        const std::string csv_eppi0 = argv[2];
        const std::string csv_bh = argv[3];

        const std::string theory_root = (argc >= 5) ? argv[4] : "output/jsons/cross_sections";
        const std::string out_root = (argc >= 6) ? argv[5] : "output/norm_scenario_comparison";

        const std::vector<std::string> labels = parse_label_list(argc, argv, 6);

        ensure_dir(out_root);

        for (const std::string &label : labels) {
            std::cout << "[norm_scenario_comparison] Processing label: " << label << "\n";

            std::vector<bool> avail0;
            std::vector<bool> avail1;
            std::vector<bool> avail2;

            std::vector<CsvData> csvs;
            csvs.push_back(read_csv_for_label(csv_unity, "e#pi^{0}=1, BH=1", label, &avail0));
            csvs.push_back(read_csv_for_label(csv_eppi0, "e#pi^{0} float, BH=1", label, &avail1));
            csvs.push_back(read_csv_for_label(csv_bh, "e#pi^{0}=1, BH float", label, &avail2));

            std::vector<bool> available(3, false);
            for (int ih = 0; ih < 3; ++ih) {
                available[(size_t)ih] = avail0[(size_t)ih] || avail1[(size_t)ih] || avail2[(size_t)ih];
            }

            if (!available[0] && !available[1] && !available[2]) {
                std::cout << "[norm_scenario_comparison]   no normed cross-section data found for "
                          << label << "; skipping.\n";
                continue;
            }

            const std::map<size_t, TheoryCurves> theory = load_theory_for_label(label, theory_root);
            if (theory.empty()) {
                std::cout << "[norm_scenario_comparison]   no theory loaded for " << label
                          << "; plots will still be attempted but overlays/ratios may be empty.\n";
            }

            const fs::path outdir = fs::path(out_root) / canonical_dir(label);
            ensure_dir(outdir);

            const std::vector<Range> xb_ranges = collect_xb_ranges(csvs);
            for (size_t ixb = 0; ixb < xb_ranges.size(); ++ixb) {
                for (int ih = 0; ih < 3; ++ih) {
                    if (!available[(size_t)ih]) continue;
                    make_cross_section_canvas(label, xb_ranges[ixb], (int)ixb, ih, csvs, theory, outdir);
                    make_bh_ratio_canvas(label, xb_ranges[ixb], (int)ixb, ih, csvs, theory, outdir);
                }
            }

            std::cout << "[norm_scenario_comparison]   saved plots under "
                      << outdir.string() << "\n";
        }

        std::cout << "[norm_scenario_comparison] Done. Output root: " << out_root << "\n";
    } catch (const std::exception &e) {
        std::cerr << "[norm_scenario_comparison] FATAL: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
