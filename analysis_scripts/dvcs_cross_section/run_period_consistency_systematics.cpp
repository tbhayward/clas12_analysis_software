#include "run_period_consistency_systematics.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TH1D.h>
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
#include <set>
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

struct TupleValue {
    bool ok = false;
    double value = 0.0;
    double stat = 0.0;
};

struct ColumnInput {
    std::string source_period;
    std::string column;
};

struct PeriodInput {
    std::string period;
    std::vector<ColumnInput> columns;
};

struct ConsistencyCase {
    std::string label;
    std::string file_tag;
    std::string avg_label;
    std::string helicity;
    std::string reference_column;
    std::vector<PeriodInput> inputs;
};

struct PointPack {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ey;
};

struct VariableConfig {
    std::string key;
    std::string column_prefix;
    std::string title;
    int summary_bins = 8;
};

struct PeriodRatioAccumulator {
    double sum_w = 0.0;
    double sum_wr = 0.0;
    int n = 0;
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

static CsvTable read_csv_or_throw(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Could not open CSV: " + path);
    CsvTable table;
    std::string line;
    if (!std::getline(fin, line)) throw std::runtime_error("CSV is empty: " + path);
    table.header = split_csv_line(line);
    std::set<std::string> seen;
    std::set<std::string> dup;
    for (int i = 0; i < (int)table.header.size(); ++i) {
        const std::string name = trim(table.header[(size_t)i]);
        if (!seen.insert(name).second) dup.insert(name);
        table.index[name] = i;
    }
    if (!dup.empty()) {
        std::ostringstream msg;
        msg << "CSV header contains duplicate columns:";
        for (const auto& d : dup) msg << "\n  - " << d;
        throw std::runtime_error(msg.str());
    }
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = split_csv_line(line);
        if (row.size() < table.header.size()) row.resize(table.header.size());
        if (row.size() != table.header.size()) throw std::runtime_error("CSV row/header size mismatch: " + path);
        table.rows.push_back(std::move(row));
    }
    return table;
}

static void require_columns(const CsvTable& table,
                            const std::vector<std::string>& required,
                            const std::string& context) {
    std::vector<std::string> missing;
    for (const auto& name : required) {
        if (table.index.find(name) == table.index.end()) missing.push_back(name);
    }
    if (!missing.empty()) {
        std::ostringstream msg;
        msg << "Missing required columns for " << context << ":";
        for (const auto& name : missing) msg << "\n  - " << name;
        throw std::runtime_error(msg.str());
    }
}

static bool parse_double(const std::string& raw, double& out) {
    const std::string s = trim(raw);
    if (s.empty()) return false;
    char* end = nullptr;
    out = std::strtod(s.c_str(), &end);
    if (end == s.c_str()) return false;
    while (end && *end != '\0') {
        if (!std::isspace((unsigned char)(*end))) return false;
        ++end;
    }
    return std::isfinite(out);
}

static bool parse_numeric_or_tuple_first(const std::string& raw, double& out) {
    if (parse_double(raw, out)) return true;
    std::string s = trim(raw);
    if (s.empty()) return false;
    if (s.front() == '(' && s.back() == ')') s = trim(s.substr(1, s.size() - 2));
    const std::vector<std::string> fields = split_csv_line(s);
    if (fields.empty()) return false;
    return parse_double(fields[0], out);
}

static TupleValue parse_tuple_value(const std::string& raw) {
    TupleValue out;
    std::string s = trim(raw);
    if (s.empty()) return out;
    if (s.front() == '(' && s.back() == ')') s = trim(s.substr(1, s.size() - 2));
    const std::vector<std::string> fields = split_csv_line(s);
    if (fields.size() < 2) return out;
    double value = 0.0;
    double stat = 0.0;
    if (!parse_double(fields[0], value)) return out;
    if (!parse_double(fields[1], stat)) return out;
    if (!std::isfinite(value) || !std::isfinite(stat) || stat <= 0.0) return out;
    out.ok = true;
    out.value = value;
    out.stat = stat;
    return out;
}

static double get_numeric_or_nan(const CsvTable& table,
                                 const std::vector<std::string>& row,
                                 const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) throw std::runtime_error("Missing required column: " + column);
    double value = 0.0;
    if (!parse_numeric_or_tuple_first(row[(size_t)it->second], value)) return std::numeric_limits<double>::quiet_NaN();
    return value;
}

static TupleValue get_tuple(const CsvTable& table,
                            const std::vector<std::string>& row,
                            const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) throw std::runtime_error("Missing required tuple column: " + column);
    return parse_tuple_value(row[(size_t)it->second]);
}

static std::string cross_section_column(const std::string& label, const std::string& helicity) {
    return "normed cross sections, ep->epg, exp, " + label + ", " + helicity;
}

static std::string avg_column(const std::string& variable_prefix, const std::string& avg_label) {
    return variable_prefix + ", " + avg_label;
}

static ColumnInput column_input(const std::string& source_period, const std::string& column) {
    return ColumnInput{source_period, column};
}

static PeriodInput single_input(const std::string& label, const std::string& source_period, const std::string& column) {
    PeriodInput input;
    input.period = label;
    input.columns.push_back(column_input(source_period, column));
    return input;
}

static PeriodInput grouped_input(const std::string& label, const std::vector<ColumnInput>& columns) {
    PeriodInput input;
    input.period = label;
    input.columns = columns;
    return input;
}

static std::vector<ConsistencyCase> consistency_cases() {
    return {
        {
            "10.6 GeV unpol", "10p6_GeV_unpol", "10.6 GeV", "unpol",
            cross_section_column("10.6 GeV", "unpol"),
            {
                single_input("Fa18 Inb", "Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                single_input("Fa18 Out", "Fa18 Out", cross_section_column("Fa18 Out", "unpol")),
                single_input("Sp18 Inb", "Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
                single_input("Sp18 Out", "Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
            }
        },
        {
            "Fa18 unpol", "Fa18_unpol", "Fa18", "unpol",
            cross_section_column("Fa18", "unpol"),
            {
                single_input("Fa18 Inb", "Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                single_input("Fa18 Out", "Fa18 Out", cross_section_column("Fa18 Out", "unpol"))
            }
        },
        {
            "Fa18 pos", "Fa18_pos", "Fa18", "pos",
            cross_section_column("Fa18", "pos"),
            {
                single_input("Fa18 Inb", "Fa18 Inb", cross_section_column("Fa18 Inb", "pos")),
                single_input("Fa18 Out", "Fa18 Out", cross_section_column("Fa18 Out", "pos"))
            }
        },
        {
            "Fa18 neg", "Fa18_neg", "Fa18", "neg",
            cross_section_column("Fa18", "neg"),
            {
                single_input("Fa18 Inb", "Fa18 Inb", cross_section_column("Fa18 Inb", "neg")),
                single_input("Fa18 Out", "Fa18 Out", cross_section_column("Fa18 Out", "neg"))
            }
        },
        {
            "Sp18 unpol", "Sp18_unpol", "Sp18", "unpol",
            cross_section_column("Sp18", "unpol"),
            {
                single_input("Sp18 Inb", "Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
                single_input("Sp18 Out", "Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
            }
        },
        {
            "Fa18 vs Sp18 unpol", "Fa18_vs_Sp18_unpol", "10.6 GeV", "unpol",
            cross_section_column("10.6 GeV", "unpol"),
            {
                single_input("Fa18", "", cross_section_column("Fa18", "unpol")),
                single_input("Sp18", "", cross_section_column("Sp18", "unpol"))
            }
        },
        {
            "Inb vs Out unpol", "Inb_vs_Out_unpol", "10.6 GeV", "unpol",
            cross_section_column("10.6 GeV", "unpol"),
            {
                grouped_input("Inb", {
                    column_input("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                    column_input("Sp18 Inb", cross_section_column("Sp18 Inb", "unpol"))
                }),
                grouped_input("Out", {
                    column_input("Fa18 Out", cross_section_column("Fa18 Out", "unpol")),
                    column_input("Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
                })
            }
        }
    };
}

static std::vector<VariableConfig> variable_configs() {
    return {
        {"xB", "xBavg", "x_{B}", 8},
        {"Q2", "Q2avg", "Q^{2} (GeV^{2})", 8},
        {"t", "t_abs_avg", "-t (GeV^{2})", 8},
        {"phi", "phiavg", "#phi (deg)", 12}
    };
}

static void validate_schema(const CsvTable& table,
                            const std::vector<ConsistencyCase>& cases,
                            const std::vector<VariableConfig>& vars) {
    std::vector<std::string> required;
    for (const auto& c : cases) {
        required.push_back(c.reference_column);
        for (const auto& input : c.inputs) {
            for (const auto& col : input.columns) required.push_back(col.column);
        }
        for (const auto& var : vars) required.push_back(avg_column(var.column_prefix, c.avg_label));
    }
    require_columns(table, required, "run-period consistency systematics");
}

static bool weighted_mean(const std::vector<TupleValue>& values, TupleValue& out) {
    double sum_w = 0.0;
    double sum_wx = 0.0;
    for (const auto& v : values) {
        if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) continue;
        const double w = 1.0 / (v.stat * v.stat);
        sum_w += w;
        sum_wx += w * v.value;
    }
    if (sum_w <= 0.0) return false;
    out.ok = true;
    out.value = sum_wx / sum_w;
    out.stat = 1.0 / std::sqrt(sum_w);
    return std::isfinite(out.value) && std::isfinite(out.stat) && out.stat > 0.0;
}

static TupleValue scaled_tuple(const TupleValue& input,
                               const std::string& source_period,
                               bool use_scaling,
                               const std::map<std::string, double>& scale_by_period) {
    if (!input.ok || !use_scaling || source_period.empty()) return input;
    const auto it = scale_by_period.find(source_period);
    if (it == scale_by_period.end()) return input;
    const double scale = it->second;
    if (!std::isfinite(scale) || std::abs(scale) <= 0.0) return input;
    TupleValue out = input;
    out.value = input.value / scale;
    out.stat = input.stat / std::abs(scale);
    return out;
}

static TupleValue evaluate_input_for_row(const CsvTable& table,
                                         const std::vector<std::string>& row,
                                         const PeriodInput& input,
                                         bool use_scaling,
                                         const std::map<std::string, double>& scale_by_period) {
    std::vector<TupleValue> values;
    for (const auto& col : input.columns) {
        TupleValue v = get_tuple(table, row, col.column);
        v = scaled_tuple(v, col.source_period, use_scaling, scale_by_period);
        if (v.ok) values.push_back(v);
    }
    TupleValue out;
    weighted_mean(values, out);
    return out;
}

static double compute_chi2_reduced_against_ref(const std::vector<TupleValue>& values, double ref_value) {
    int n = 0;
    double chi2 = 0.0;
    for (const auto& v : values) {
        if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) continue;
        const double residual = v.value - ref_value;
        chi2 += residual * residual / (v.stat * v.stat);
        ++n;
    }
    const double ndf = (double)n - 1.0;
    if (ndf <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return chi2 / ndf;
}

static std::vector<std::string> base_periods() {
    return {"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"};
}

static std::map<std::string, double> compute_10p6_unpol_scale_factors(const CsvTable& table) {
    const std::vector<std::string> periods = base_periods();
    std::map<std::string, PeriodRatioAccumulator> acc_by_period;
    for (const auto& period : periods) acc_by_period[period] = PeriodRatioAccumulator();

    for (const auto& row : table.rows) {
        const TupleValue ref = get_tuple(table, row, cross_section_column("10.6 GeV", "unpol"));
        if (!ref.ok || !std::isfinite(ref.value) || std::abs(ref.value) <= 0.0) continue;
        for (const auto& period : periods) {
            const TupleValue v = get_tuple(table, row, cross_section_column(period, "unpol"));
            if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) continue;
            const double ratio = v.value / ref.value;
            const double ratio_stat = std::abs(v.stat / ref.value);
            if (!std::isfinite(ratio) || !std::isfinite(ratio_stat) || ratio_stat <= 0.0) continue;
            const double w = 1.0 / (ratio_stat * ratio_stat);
            PeriodRatioAccumulator& acc = acc_by_period[period];
            acc.sum_w += w;
            acc.sum_wr += w * ratio;
            acc.n += 1;
        }
    }

    std::map<std::string, double> scale_by_period;
    std::cout << "[run-period-consistency] 10.6 GeV unpol scale factors used for scaled plots:\n";
    for (const auto& period : periods) {
        const PeriodRatioAccumulator& acc = acc_by_period[period];
        if (acc.sum_w <= 0.0 || acc.n <= 0) throw std::runtime_error("Could not compute scale factor for period " + period);
        const double mean_ratio = acc.sum_wr / acc.sum_w;
        const double mean_ratio_stat = 1.0 / std::sqrt(acc.sum_w);
        if (!std::isfinite(mean_ratio) || std::abs(mean_ratio) <= 0.0) throw std::runtime_error("Invalid scale factor for period " + period);
        scale_by_period[period] = mean_ratio;
        std::cout << "  " << period << " scale = " << std::setprecision(10) << mean_ratio
                  << " +/- " << mean_ratio_stat << "   n=" << acc.n << "\n";
    }
    return scale_by_period;
}

static std::vector<double> choose_x_range(const std::vector<double>& xs) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    for (const double x : xs) {
        if (!std::isfinite(x)) continue;
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
    }
    if (!std::isfinite(xmin) || !std::isfinite(xmax)) return {0.0, 1.0};
    if (xmin == xmax) {
        const double pad = (std::abs(xmin) > 0.0) ? 0.05 * std::abs(xmin) : 1.0;
        return {xmin - pad, xmax + pad};
    }
    const double pad = 0.07 * (xmax - xmin);
    return {xmin - pad, xmax + pad};
}

static double percentile(std::vector<double> values, double q) {
    std::vector<double> clean;
    for (const double v : values) if (std::isfinite(v)) clean.push_back(v);
    if (clean.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(clean.begin(), clean.end());
    if (q <= 0.0) return clean.front();
    if (q >= 1.0) return clean.back();
    const double pos = q * (double)(clean.size() - 1);
    const size_t lo = (size_t)std::floor(pos);
    const size_t hi = (size_t)std::ceil(pos);
    const double frac = pos - (double)lo;
    if (lo == hi) return clean[lo];
    return clean[lo] * (1.0 - frac) + clean[hi] * frac;
}

static std::vector<double> choose_chi2_y_range(const PointPack& pack) {
    std::vector<double> ys;
    for (const double y : pack.y) if (std::isfinite(y) && y >= 0.0) ys.push_back(y);
    if (ys.empty()) return {0.0, 5.0};
    double ymax = std::max(5.0, 1.25 * percentile(ys, 0.95));
    if (!std::isfinite(ymax) || ymax <= 0.0) ymax = 5.0;
    return {0.0, ymax};
}

static std::vector<double> choose_pull_y_range(const std::map<std::string, PointPack>& by_period) {
    double ymin = -5.0;
    double ymax = 5.0;
    std::vector<double> vals;
    for (const auto& kv : by_period) for (const double y : kv.second.y) if (std::isfinite(y)) vals.push_back(y);
    if (!vals.empty()) {
        ymin = std::min(-5.0, 1.25 * percentile(vals, 0.02));
        ymax = std::max(5.0, 1.25 * percentile(vals, 0.98));
    }
    return {ymin, ymax};
}

static void set_pad_style() {
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.16);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
}

static TGraph* make_marker_graph(const PointPack& pack, int color, int marker_style, double marker_size) {
    TGraph* graph = new TGraph((int)pack.x.size());
    for (int i = 0; i < (int)pack.x.size(); ++i) graph->SetPoint(i, pack.x[(size_t)i], pack.y[(size_t)i]);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetMarkerStyle(marker_style);
    graph->SetMarkerSize(marker_size);
    graph->SetLineWidth(1);
    return graph;
}

static TH1D* make_frame(const std::string& name,
                        const std::string& xtitle,
                        const std::string& ytitle,
                        double xmin, double xmax, double ymin, double ymax) {
    TH1D* frame = new TH1D(name.c_str(), "", 100, xmin, xmax);
    frame->SetDirectory(nullptr);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->GetXaxis()->SetTitle(xtitle.c_str());
    frame->GetYaxis()->SetTitle(ytitle.c_str());
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetTitleSize(0.060);
    frame->GetYaxis()->SetTitleSize(0.060);
    frame->GetXaxis()->SetTitleOffset(1.10);
    frame->GetYaxis()->SetTitleOffset(1.20);
    frame->GetXaxis()->SetLabelSize(0.047);
    frame->GetYaxis()->SetLabelSize(0.047);
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(505);
    return frame;
}

static void draw_subpad_label(const std::string& text) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.058);
    latex.SetTextAlign(22);
    latex.DrawLatex(0.50, 0.940, text.c_str());
}

static void draw_top_title(TPad* title_pad, const std::string& title_text) {
    title_pad->cd();
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.34);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.50, title_text.c_str());
}

static void make_title_and_grid_pads(TCanvas& canvas, TPad*& title_pad, TPad*& grid_pad) {
    canvas.cd();
    title_pad = new TPad("title_pad", "title_pad", 0.0, 0.925, 1.0, 1.0);
    title_pad->SetFillColor(kWhite);
    title_pad->SetFrameFillColor(kWhite);
    title_pad->SetBorderMode(0);
    title_pad->SetMargin(0.0, 0.0, 0.0, 0.0);
    title_pad->Draw();
    canvas.cd();
    grid_pad = new TPad("grid_pad", "grid_pad", 0.0, 0.0, 1.0, 0.925);
    grid_pad->SetFillColor(kWhite);
    grid_pad->SetFrameFillColor(kWhite);
    grid_pad->SetBorderMode(0);
    grid_pad->SetMargin(0.02, 0.02, 0.02, 0.02);
    grid_pad->Draw();
    grid_pad->Divide(2, 2, 0.001, 0.001);
}

static std::string plot_prefix(bool use_scaling) {
    return use_scaling ? "scaled" : "unscaled";
}

static std::string title_prefix(bool use_scaling) {
    return use_scaling ? "Scaled " : "Unscaled ";
}

static void draw_chi2_canvas(const fs::path& out_dir,
                             const ConsistencyCase& c,
                             const std::map<std::string, PointPack>& chi2_by_var,
                             bool use_scaling) {
    const std::string prefix = plot_prefix(use_scaling);
    TCanvas canvas(("c_chi2_" + c.file_tag + "_" + prefix).c_str(), "c_chi2", 1600, 1200);
    canvas.SetFillColor(kWhite);
    TPad* title_pad = nullptr;
    TPad* grid_pad = nullptr;
    make_title_and_grid_pads(canvas, title_pad, grid_pad);
    draw_top_title(title_pad, title_prefix(use_scaling) + "run-period consistency: " + c.label);
    const std::vector<VariableConfig> vars = variable_configs();
    std::vector<std::unique_ptr<TH1D> > frames;
    std::vector<std::unique_ptr<TGraph> > graphs;
    std::vector<std::unique_ptr<TLine> > lines;
    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        grid_pad->cd(ivar + 1);
        set_pad_style();
        PointPack pack;
        const auto it = chi2_by_var.find(vars[(size_t)ivar].key);
        if (it != chi2_by_var.end()) pack = it->second;
        const std::vector<double> xr = choose_x_range(pack.x);
        const std::vector<double> yr = choose_chi2_y_range(pack);
        std::unique_ptr<TH1D> frame(make_frame("frame_chi2_" + c.file_tag + "_" + prefix + "_" + vars[(size_t)ivar].key,
                                               vars[(size_t)ivar].title, "#chi^{2}/ndf", xr[0], xr[1], yr[0], yr[1]));
        frame->Draw("AXIS");
        std::unique_ptr<TLine> line(new TLine(xr[0], 1.0, xr[1], 1.0));
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kRed + 1);
        line->Draw("SAME");
        if (!pack.x.empty()) {
            std::unique_ptr<TGraph> graph(make_marker_graph(pack, kBlack, 20, 0.55));
            graph->Draw("P SAME");
            graphs.push_back(std::move(graph));
        }
        draw_subpad_label(vars[(size_t)ivar].title);
        frames.push_back(std::move(frame));
        lines.push_back(std::move(line));
    }
    canvas.cd();
    canvas.Modified();
    canvas.Update();
    canvas.SaveAs((out_dir / (prefix + "_reduced_chi2.pdf")).string().c_str());
}

static void draw_multi_period_canvas(const fs::path& out_dir,
                                     const ConsistencyCase& c,
                                     const std::map<std::string, std::map<std::string, PointPack> >& by_var_period,
                                     bool use_scaling,
                                     const std::string& kind,
                                     const std::string& ytitle,
                                     double ymin_fixed,
                                     double ymax_fixed,
                                     bool fixed_y) {
    const std::string prefix = plot_prefix(use_scaling);
    TCanvas canvas(("c_" + kind + "_" + c.file_tag + "_" + prefix).c_str(), "c_multi", 1600, 1200);
    canvas.SetFillColor(kWhite);
    TPad* title_pad = nullptr;
    TPad* grid_pad = nullptr;
    make_title_and_grid_pads(canvas, title_pad, grid_pad);
    draw_top_title(title_pad, title_prefix(use_scaling) + kind + ": " + c.label);
    const std::vector<VariableConfig> vars = variable_configs();
    const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7};
    const std::vector<int> markers = {20, 24, 21, 25, 22, 26};
    std::vector<std::unique_ptr<TH1D> > frames;
    std::vector<std::unique_ptr<TGraph> > graphs;
    std::vector<std::unique_ptr<TLine> > lines;
    std::vector<std::unique_ptr<TLegend> > legends;

    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        grid_pad->cd(ivar + 1);
        set_pad_style();
        std::vector<double> all_x;
        const auto it_var = by_var_period.find(vars[(size_t)ivar].key);
        if (it_var != by_var_period.end()) {
            for (const auto& input : c.inputs) {
                const auto it_pack = it_var->second.find(input.period);
                if (it_pack == it_var->second.end()) continue;
                for (const double x : it_pack->second.x) all_x.push_back(x);
            }
        }
        const std::vector<double> xr = choose_x_range(all_x);
        double ymin = ymin_fixed;
        double ymax = ymax_fixed;
        if (!fixed_y && it_var != by_var_period.end()) {
            const std::vector<double> pr = choose_pull_y_range(it_var->second);
            ymin = pr[0];
            ymax = pr[1];
        }
        std::unique_ptr<TH1D> frame(make_frame("frame_" + kind + "_" + c.file_tag + "_" + prefix + "_" + vars[(size_t)ivar].key,
                                               vars[(size_t)ivar].title, ytitle, xr[0], xr[1], ymin, ymax));
        frame->Draw("AXIS");
        std::unique_ptr<TLine> line(new TLine(xr[0], (kind == "pulls" ? 0.0 : 1.0), xr[1], (kind == "pulls" ? 0.0 : 1.0)));
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kRed + 1);
        line->Draw("SAME");
        std::unique_ptr<TLegend> legend(new TLegend(0.54, 0.68, 0.93, 0.89));
        legend->SetBorderSize(1);
        legend->SetFillStyle(1001);
        legend->SetFillColor(kWhite);
        legend->SetTextFont(42);
        legend->SetTextSize(0.040);
        int iper = 0;
        if (it_var != by_var_period.end()) {
            for (const auto& input : c.inputs) {
                const auto it_pack = it_var->second.find(input.period);
                if (it_pack == it_var->second.end() || it_pack->second.x.empty()) {
                    ++iper;
                    continue;
                }
                const int color = colors[(size_t)(iper % (int)colors.size())];
                const int marker = markers[(size_t)(iper % (int)markers.size())];
                std::unique_ptr<TGraph> graph(make_marker_graph(it_pack->second, color, marker, 0.48));
                graph->Draw("P SAME");
                legend->AddEntry(graph.get(), input.period.c_str(), "p");
                graphs.push_back(std::move(graph));
                ++iper;
            }
        }
        legend->Draw();
        draw_subpad_label(vars[(size_t)ivar].title);
        frames.push_back(std::move(frame));
        lines.push_back(std::move(line));
        legends.push_back(std::move(legend));
    }
    canvas.cd();
    canvas.Modified();
    canvas.Update();
    canvas.SaveAs((out_dir / (prefix + "_" + kind + ".pdf")).string().c_str());
}

static void fill_case_outputs(const CsvTable& table,
                              const ConsistencyCase& c,
                              const std::map<std::string, double>& scale_by_period,
                              bool use_scaling,
                              std::map<std::string, PointPack>& chi2_by_var,
                              std::map<std::string, std::map<std::string, PointPack> >& ratios_by_var_period,
                              std::map<std::string, std::map<std::string, PointPack> >& pulls_by_var_period) {
    const std::vector<VariableConfig> vars = variable_configs();
    for (const auto& row : table.rows) {
        const TupleValue ref = get_tuple(table, row, c.reference_column);
        if (!ref.ok || !std::isfinite(ref.value) || std::abs(ref.value) <= 0.0) continue;

        std::vector<TupleValue> values;
        values.reserve(c.inputs.size());
        for (const auto& input : c.inputs) {
            values.push_back(evaluate_input_for_row(table, row, input, use_scaling, scale_by_period));
        }
        const double chi2_ndf = compute_chi2_reduced_against_ref(values, ref.value);
        if (!std::isfinite(chi2_ndf)) continue;

        for (const auto& var : vars) {
            const double x = get_numeric_or_nan(table, row, avg_column(var.column_prefix, c.avg_label));
            if (!std::isfinite(x)) continue;
            chi2_by_var[var.key].x.push_back(x);
            chi2_by_var[var.key].y.push_back(chi2_ndf);
            chi2_by_var[var.key].ey.push_back(0.0);

            for (size_t iper = 0; iper < c.inputs.size(); ++iper) {
                const TupleValue& v = values[iper];
                if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) continue;
                const double ratio = v.value / ref.value;
                const double pull = (v.value - ref.value) / v.stat;
                if (std::isfinite(ratio)) {
                    PointPack& pack = ratios_by_var_period[var.key][c.inputs[iper].period];
                    pack.x.push_back(x);
                    pack.y.push_back(ratio);
                    pack.ey.push_back(0.0);
                }
                if (std::isfinite(pull)) {
                    PointPack& pack = pulls_by_var_period[var.key][c.inputs[iper].period];
                    pack.x.push_back(x);
                    pack.y.push_back(pull);
                    pack.ey.push_back(0.0);
                }
            }
        }
    }
}

static double mean_vec(const std::vector<double>& v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    double sum = 0.0;
    for (const double x : v) sum += x;
    return sum / (double)v.size();
}

static double rms_vec(const std::vector<double>& v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double m = mean_vec(v);
    double sum = 0.0;
    for (const double x : v) sum += (x - m) * (x - m);
    return std::sqrt(sum / (double)v.size());
}

static void write_binned_chi2_summary(const fs::path& path,
                                      const ConsistencyCase& c,
                                      bool use_scaling,
                                      const std::map<std::string, PointPack>& chi2_by_var) {
    std::ofstream fout(path);
    if (!fout.is_open()) throw std::runtime_error("Could not write chi2 summary: " + path.string());
    fout << "case,scaled,variable,bin index,x min,x max,n,mean chi2ndf,rms chi2ndf\n";
    std::cout << "[run-period-consistency] Binned chi2/ndf summary for " << (use_scaling ? "scaled " : "") << c.label << "\n";
    for (const auto& var : variable_configs()) {
        const auto it = chi2_by_var.find(var.key);
        if (it == chi2_by_var.end() || it->second.x.empty()) continue;
        const PointPack& pack = it->second;
        const std::vector<double> xr = choose_x_range(pack.x);
        const double xmin = xr[0];
        const double xmax = xr[1];
        const int nb = var.summary_bins;
        std::vector<std::vector<double> > bins((size_t)nb);
        for (size_t i = 0; i < pack.x.size(); ++i) {
            const double x = pack.x[i];
            const double y = pack.y[i];
            if (!std::isfinite(x) || !std::isfinite(y)) continue;
            int ib = (int)std::floor((x - xmin) / (xmax - xmin) * (double)nb);
            if (ib < 0) ib = 0;
            if (ib >= nb) ib = nb - 1;
            bins[(size_t)ib].push_back(y);
        }
        for (int ib = 0; ib < nb; ++ib) {
            const double lo = xmin + (xmax - xmin) * (double)ib / (double)nb;
            const double hi = xmin + (xmax - xmin) * (double)(ib + 1) / (double)nb;
            const double m = mean_vec(bins[(size_t)ib]);
            const double r = rms_vec(bins[(size_t)ib]);
            fout << c.label << "," << (use_scaling ? 1 : 0) << "," << var.key << "," << ib << ","
                 << std::setprecision(12) << lo << "," << hi << "," << bins[(size_t)ib].size() << ","
                 << m << "," << r << "\n";
            std::cout << "    " << var.key << " bin " << ib << " (" << std::setprecision(5) << lo << ", " << hi
                      << "): n=" << bins[(size_t)ib].size() << " mean=" << std::setprecision(8) << m
                      << " rms=" << r << "\n";
        }
    }
}

static void write_pull_summary(const fs::path& path,
                               const ConsistencyCase& c,
                               bool use_scaling,
                               const std::map<std::string, std::map<std::string, PointPack> >& pulls_by_var_period) {
    std::ofstream fout(path);
    if (!fout.is_open()) throw std::runtime_error("Could not write pull summary: " + path.string());
    fout << "case,scaled,variable,period,n,mean pull,rms pull\n";
    std::cout << "[run-period-consistency] Pull summary for " << (use_scaling ? "scaled " : "") << c.label << "\n";
    for (const auto& var : variable_configs()) {
        const auto it_var = pulls_by_var_period.find(var.key);
        if (it_var == pulls_by_var_period.end()) continue;
        for (const auto& kv : it_var->second) {
            const double m = mean_vec(kv.second.y);
            const double r = rms_vec(kv.second.y);
            fout << c.label << "," << (use_scaling ? 1 : 0) << "," << var.key << "," << kv.first << ","
                 << kv.second.y.size() << "," << std::setprecision(12) << m << "," << r << "\n";
            if (var.key == "phi") {
                std::cout << "    pull " << std::left << std::setw(10) << kv.first << std::right
                          << ": n=" << kv.second.y.size() << " mean=" << std::setprecision(8) << m
                          << " rms=" << r << "\n";
            }
        }
    }
}

static void make_plots_for_case(const fs::path& root_out_dir,
                                const CsvTable& table,
                                const ConsistencyCase& c,
                                const std::map<std::string, double>& scale_by_period,
                                bool use_scaling) {
    const fs::path out_dir = root_out_dir / c.file_tag;
    fs::create_directories(out_dir);

    std::map<std::string, PointPack> chi2_by_var;
    std::map<std::string, std::map<std::string, PointPack> > ratios_by_var_period;
    std::map<std::string, std::map<std::string, PointPack> > pulls_by_var_period;
    fill_case_outputs(table, c, scale_by_period, use_scaling, chi2_by_var, ratios_by_var_period, pulls_by_var_period);

    draw_chi2_canvas(out_dir, c, chi2_by_var, use_scaling);
    draw_multi_period_canvas(out_dir, c, ratios_by_var_period, use_scaling, "period_ratios", "#sigma_{i}/#sigma_{CSV ref}", 0.0, 2.0, true);
    draw_multi_period_canvas(out_dir, c, pulls_by_var_period, use_scaling, "pulls", "(#sigma_{i}-#sigma_{CSV ref})/#delta#sigma_{i}", -5.0, 5.0, false);
    write_binned_chi2_summary(out_dir / (plot_prefix(use_scaling) + "_binned_chi2_summary.csv"), c, use_scaling, chi2_by_var);
    write_pull_summary(out_dir / (plot_prefix(use_scaling) + "_pull_summary.csv"), c, use_scaling, pulls_by_var_period);

    size_t npoints = 0;
    const auto it = chi2_by_var.find("phi");
    if (it != chi2_by_var.end()) npoints = it->second.x.size();
    std::cout << "[run-period-consistency] Completed " << (use_scaling ? "scaled " : "")
              << c.label << " with " << npoints << " valid bins.\n";
}

} // namespace

bool run_period_consistency_systematics(const std::string& csv_path,
                                        const std::string& output_root_dir) {
    try {
        const fs::path out_dir = fs::path(output_root_dir) / "run_period_consistency";
        fs::create_directories(out_dir);
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTitleFont(42, "");
        gStyle->SetTextFont(42);
        gStyle->SetEndErrorSize(0);

        const CsvTable table = read_csv_or_throw(csv_path);
        const std::vector<ConsistencyCase> cases = consistency_cases();
        const std::vector<VariableConfig> vars = variable_configs();
        validate_schema(table, cases, vars);
        const std::map<std::string, double> scale_by_period = compute_10p6_unpol_scale_factors(table);

        std::cout << "[run-period-consistency] CSV rows loaded: " << table.rows.size() << "\n";
        std::cout << "[run-period-consistency] Output directory: " << out_dir.string() << "\n";
        for (const auto& c : cases) {
            make_plots_for_case(out_dir, table, c, scale_by_period, false);
            make_plots_for_case(out_dir, table, c, scale_by_period, true);
        }
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[run-period-consistency] FATAL: " << e.what() << "\n";
        return false;
    }
}
