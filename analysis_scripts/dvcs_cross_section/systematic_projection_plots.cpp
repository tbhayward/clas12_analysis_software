#include "systematic_projection_plots.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TPad.h>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string> > rows;
    std::unordered_map<std::string, int> index;
};

struct Component {
    std::string column;
    std::string label;
    int marker = 20;
    int color = 1;
    int line_style = 1;
};

struct BinAccumulator {
    double x_sum = 0.0;
    int x_count = 0;
    std::vector<double> abs_sum;
    std::vector<double> rel_sum;
    std::vector<int> abs_count;
    std::vector<int> rel_count;
};

static std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;
    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
    return s.substr(b, e - b);
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

static CsvTable read_csv(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Could not open CSV: " + path);
    CsvTable t;
    std::string line;
    if (!std::getline(fin, line)) throw std::runtime_error("Empty CSV: " + path);
    t.header = split_csv_line(line);
    for (int i = 0; i < (int)t.header.size(); ++i) t.index[t.header[(size_t)i]] = i;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        auto row = split_csv_line(line);
        row.resize(t.header.size());
        t.rows.push_back(std::move(row));
    }
    return t;
}

static double number(const std::string& s) {
    const std::string x = trim(s);
    if (x.empty()) return std::numeric_limits<double>::quiet_NaN();
    char* end = nullptr;
    const double v = std::strtod(x.c_str(), &end);
    return end == x.c_str() ? std::numeric_limits<double>::quiet_NaN() : v;
}

static std::string cell(const CsvTable& t, const std::vector<std::string>& row,
                        const std::string& col) {
    auto it = t.index.find(col);
    if (it == t.index.end()) return "";
    return row[(size_t)it->second];
}

static bool parse_tuple_first(const std::string& raw, double& value) {
    std::string s = trim(raw);
    while (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
        s = trim(s.substr(1, s.size() - 2));
    }
    if (s.size() < 3 || s.front() != '(' || s.back() != ')') return false;
    s = s.substr(1, s.size() - 2);
    const size_t comma = s.find(',');
    const std::string first = comma == std::string::npos ? s : s.substr(0, comma);
    value = number(first);
    return std::isfinite(value);
}

static std::vector<Component> components() {
    return {
        {"Syst. err (pi0 subtraction)", "#pi^{0} subtraction", 20, kBlue + 1, 1},
        {"Syst. err (Acceptance)", "Acceptance", 21, kRed + 1, 1},
        {"Syst.err (Frad)", "F_{rad}", 22, kGreen + 2, 1},
        {"Syst.err (Fbin)", "F_{bin}", 23, kMagenta + 1, 1},
        {"Syst. err (exclusivity cuts)", "Exclusivity cuts", 24, kOrange + 7, 1},
        {"Syst. err (fiducial cuts)", "Fiducial cuts", 25, kCyan + 2, 1},
        {"Syst. err (point-to-point total)", "Total", 29, kBlack, 1}
    };
}

struct VariableSpec {
    std::string name;
    std::string x_title;
    std::string min_col;
    std::string max_col;
    std::string mean_col;
    double grouping_width = 0.0;
};

static std::vector<VariableSpec> variables() {
    return {
        {"xB", "x_{B}", "xBmin", "xBmax", "", 0.0},
        {"Q2", "Q^{2} (GeV^{2})", "Q2min", "Q2max", "", 0.0},
        {"t", "|t| (GeV^{2})", "t_abs_min", "t_abs_max", "", 0.0},
        {"phi", "#phi (deg)", "phimin", "phimax", "", 0.0},
        {"e_theta", "#theta_{e} (deg)", "", "", "e_theta, 10.6 GeV", 1.0},
        {"p_theta", "#theta_{p} (deg)", "", "", "p_theta, 10.6 GeV", 1.0},
        {"g_theta", "#theta_{#gamma} (deg)", "", "", "g_theta, 10.6 GeV", 1.0}
    };
}

static bool row_x(const CsvTable& t, const std::vector<std::string>& row,
                  const VariableSpec& v, double& x, std::string& key) {
    if (!v.mean_col.empty()) {
        x = number(cell(t, row, v.mean_col));
        if (!std::isfinite(x)) return false;
        const double center = v.grouping_width > 0.0
            ? std::floor(x / v.grouping_width) * v.grouping_width + 0.5 * v.grouping_width
            : x;
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(6) << center;
        key = ss.str();
        x = center;
        return true;
    }
    const double lo = number(cell(t, row, v.min_col));
    const double hi = number(cell(t, row, v.max_col));
    if (!std::isfinite(lo) || !std::isfinite(hi)) return false;
    x = 0.5 * (lo + hi);
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << lo << '|' << hi;
    key = ss.str();
    return true;
}

static void style_pad(TPad* p) {
    p->SetLeftMargin(0.105);
    p->SetRightMargin(0.025);
    p->SetBottomMargin(0.135);
    p->SetTopMargin(0.155);
    p->SetGridx(false);
    p->SetGridy(true);
    p->SetTicks(1, 1);
    p->SetLogy(true);
}

static double log_upper_bound(double maximum, double fallback) {
    if (!(maximum > 0.0) || !std::isfinite(maximum)) return fallback;
    return std::pow(10.0, std::ceil(std::log10(maximum * 1.35)));
}

static bool make_one(const CsvTable& table, const VariableSpec& var,
                     const std::string& output_dir) {
    const auto comps = components();
    std::map<std::string, BinAccumulator> bins;

    const std::string xs_col =
        "normed cross sections, ep->epg, exp, 10.6 GeV, unpol";

    for (const auto& row : table.rows) {
        double xs = std::numeric_limits<double>::quiet_NaN();
        if (!parse_tuple_first(cell(table, row, xs_col), xs) || xs == 0.0) continue;

        double x = 0.0;
        std::string key;
        if (!row_x(table, row, var, x, key)) continue;

        auto& b = bins[key];
        if (b.abs_sum.empty()) {
            b.abs_sum.assign(comps.size(), 0.0);
            b.rel_sum.assign(comps.size(), 0.0);
            b.abs_count.assign(comps.size(), 0);
            b.rel_count.assign(comps.size(), 0);
        }
        b.x_sum += x;
        ++b.x_count;

        for (size_t j = 0; j < comps.size(); ++j) {
            const double s = number(cell(table, row, comps[j].column));
            if (!std::isfinite(s) || s < 0.0) continue;
            b.abs_sum[j] += s;
            ++b.abs_count[j];
            b.rel_sum[j] += 100.0 * s / std::fabs(xs);
            ++b.rel_count[j];
        }
    }

    struct Point { double x; std::vector<double> a, r; };
    std::vector<Point> points;
    for (const auto& kv : bins) {
        const auto& b = kv.second;
        if (b.x_count <= 0) continue;
        Point p;
        p.x = b.x_sum / b.x_count;
        p.a.resize(comps.size(), std::numeric_limits<double>::quiet_NaN());
        p.r.resize(comps.size(), std::numeric_limits<double>::quiet_NaN());
        for (size_t j = 0; j < comps.size(); ++j) {
            if (b.abs_count[j] > 0) p.a[j] = b.abs_sum[j] / b.abs_count[j];
            if (b.rel_count[j] > 0) p.r[j] = b.rel_sum[j] / b.rel_count[j];
        }
        points.push_back(std::move(p));
    }
    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });
    if (points.empty()) return false;

    TCanvas c(("c_syst_" + var.name).c_str(), "", 1500, 1100);
    c.Divide(1, 2, 0.0, 0.015);
    std::vector<std::unique_ptr<TGraphErrors> > graphs_abs;
    std::vector<std::unique_ptr<TGraphErrors> > graphs_rel;

    auto draw_panel = [&](int pad_index, bool relative) {
        TPad* pad = static_cast<TPad*>(c.cd(pad_index));
        style_pad(pad);

        const double ymin = relative ? 1.0e-4 : 1.0e-6;
        double observed_max = 0.0;
        for (const auto& p : points) {
            const auto& values = relative ? p.r : p.a;
            for (double y : values) {
                if (std::isfinite(y) && y > 0.0) observed_max = std::max(observed_max, y);
            }
        }

        // The relative panel is intentionally capped at 100% so a handful of
        // pathological low-cross-section bins do not destroy the useful scale.
        const double ymax = relative ? 100.0 : log_upper_bound(observed_max, 1.0);

        TLegend* leg = new TLegend(0.12, 0.835, 0.88, 0.955);
        leg->SetBorderSize(1);
        leg->SetLineColor(kGray + 1);
        leg->SetFillColor(kWhite);
        leg->SetFillStyle(1001);
        leg->SetTextFont(42);
        leg->SetTextSize(0.031);
        leg->SetNColumns(4);
        leg->SetMargin(0.18);
        leg->SetColumnSeparation(0.08);

        bool first = true;
        for (size_t j = 0; j < comps.size(); ++j) {
            std::vector<double> x, y, ex, ey;
            for (const auto& p : points) {
                const double yy = relative ? p.r[j] : p.a[j];
                if (!std::isfinite(yy) || yy <= 0.0) continue;
                x.push_back(p.x);
                y.push_back(yy);
                ex.push_back(0.0);
                ey.push_back(0.0);
            }
            if (x.empty()) continue;

            auto g = std::make_unique<TGraphErrors>(
                static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
            g->SetMarkerStyle(comps[j].marker);
            g->SetMarkerColor(comps[j].color);
            g->SetLineColor(comps[j].color);
            g->SetLineStyle(comps[j].line_style);
            g->SetMarkerSize(comps[j].label == "Total" ? 1.15 : 0.90);
            g->SetLineWidth(comps[j].label == "Total" ? 4 : 2);
            g->SetTitle("");

            if (first) {
                g->GetXaxis()->SetTitle(var.x_title.c_str());
                g->GetYaxis()->SetTitle(relative
                    ? "Mean point-to-point systematic / |#sigma| (%)"
                    : "Mean point-to-point systematic (nb/GeV^{4})");
                g->GetXaxis()->SetTitleFont(42);
                g->GetYaxis()->SetTitleFont(42);
                g->GetXaxis()->SetLabelFont(42);
                g->GetYaxis()->SetLabelFont(42);
                g->GetXaxis()->SetTitleSize(0.047);
                g->GetYaxis()->SetTitleSize(0.047);
                g->GetXaxis()->SetLabelSize(0.039);
                g->GetYaxis()->SetLabelSize(0.039);
                g->GetYaxis()->SetTitleOffset(1.02);
                g->GetXaxis()->SetTitleOffset(1.08);
                g->GetYaxis()->SetRangeUser(ymin, ymax);
                g->GetYaxis()->SetMoreLogLabels(true);
                g->GetYaxis()->SetNoExponent(false);
                g->Draw("ALP");
                first = false;
            } else {
                g->Draw("LP SAME");
            }

            leg->AddEntry(g.get(), comps[j].label.c_str(), "lp");
            if (relative) graphs_rel.push_back(std::move(g));
            else graphs_abs.push_back(std::move(g));
        }

        leg->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.038);
        latex.SetTextAlign(31);
        latex.DrawLatex(0.975, 0.975,
            relative ? "Relative average systematic size" : "Absolute average systematic size");

        if (relative) {
            latex.SetTextAlign(11);
            latex.SetTextSize(0.027);
            latex.SetTextColor(kGray + 2);
            latex.DrawLatex(0.115, 0.785, "Values above 100% are clipped by the display range");
            latex.SetTextColor(kBlack);
        }

        pad->RedrawAxis();
    };

    draw_panel(1, false);
    draw_panel(2, true);
    c.SaveAs((output_dir + "/point_to_point_systematics_vs_" + var.name + ".png").c_str());

    std::ofstream out(output_dir + "/point_to_point_systematics_vs_" + var.name + ".csv");
    out << "x";
    for (const auto& comp : comps) out << ',' << comp.label << "_absolute," << comp.label << "_percent";
    out << '\n';
    for (const auto& p : points) {
        out << std::setprecision(12) << p.x;
        for (size_t j = 0; j < comps.size(); ++j) out << ',' << p.a[j] << ',' << p.r[j];
        out << '\n';
    }
    return true;
}

} // namespace

bool make_systematic_projection_plots(const std::string& csv_path,
                                      const std::string& output_dir) {
    try {
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);
        fs::create_directories(output_dir);
        const CsvTable table = read_csv(csv_path);

        for (const auto& comp : components()) {
            if (table.index.find(comp.column) == table.index.end()) {
                throw std::runtime_error("Missing systematic column: " + comp.column);
            }
        }

        int made = 0;
        for (const auto& var : variables()) {
            if (make_one(table, var, output_dir)) ++made;
        }
        std::cout << "[systematic-projections] Made " << made
                  << " projection canvases in " << output_dir << "\n";
        return made == (int)variables().size();
    } catch (const std::exception& e) {
        std::cerr << "[systematic-projections] ERROR: " << e.what() << "\n";
        return false;
    }
}
