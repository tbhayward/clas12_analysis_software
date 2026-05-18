#include "run_period_consistency_systematics.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TSystem.h>
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

struct PeriodInput {
    std::string period;
    std::string column;
};

struct ConsistencyCase {
    std::string label;
    std::string file_tag;
    std::string avg_label;
    std::string helicity;
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
};

static std::string trim(const std::string& s) {
    const size_t first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }

    const size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());

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

static std::unordered_map<std::string, int>
build_header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> idx;

    for (int i = 0; i < (int)header.size(); ++i) {
        idx[header[(size_t)i]] = i;
    }

    return idx;
}

static void assert_no_duplicate_columns(const std::vector<std::string>& header) {
    std::set<std::string> seen;
    std::set<std::string> duplicates;

    for (const auto& name : header) {
        if (!seen.insert(name).second) {
            duplicates.insert(name);
        }
    }

    if (!duplicates.empty()) {
        std::ostringstream msg;
        msg << "CSV header contains duplicate columns:";
        for (const auto& name : duplicates) {
            msg << "\n  - " << name;
        }
        throw std::runtime_error(msg.str());
    }
}

static CsvTable read_csv_or_throw(const std::string& path) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        throw std::runtime_error("Could not open CSV: " + path);
    }

    CsvTable table;

    std::string line;
    if (!std::getline(fin, line)) {
        throw std::runtime_error("CSV is empty: " + path);
    }

    table.header = split_csv_line(line);
    assert_no_duplicate_columns(table.header);
    table.index = build_header_index(table.header);

    while (std::getline(fin, line)) {
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> row = split_csv_line(line);

        if (row.size() < table.header.size()) {
            row.resize(table.header.size());
        }

        if (row.size() != table.header.size()) {
            std::ostringstream msg;
            msg << "CSV row has " << row.size()
                << " columns, but header has " << table.header.size()
                << " columns. File: " << path;
            throw std::runtime_error(msg.str());
        }

        table.rows.push_back(std::move(row));
    }

    return table;
}

static void require_columns(const CsvTable& table,
                            const std::vector<std::string>& required,
                            const std::string& context) {
    std::vector<std::string> missing;

    for (const auto& name : required) {
        if (table.index.find(name) == table.index.end()) {
            missing.push_back(name);
        }
    }

    if (!missing.empty()) {
        std::ostringstream msg;
        msg << "Missing required columns for " << context << ":";
        for (const auto& name : missing) {
            msg << "\n  - " << name;
        }
        throw std::runtime_error(msg.str());
    }
}

static bool parse_double(const std::string& raw, double& out) {
    const std::string s = trim(raw);
    if (s.empty()) {
        return false;
    }

    char* end = nullptr;
    out = std::strtod(s.c_str(), &end);

    if (end == s.c_str()) {
        return false;
    }

    while (end && *end != '\0') {
        if (!std::isspace((unsigned char)(*end))) {
            return false;
        }
        ++end;
    }

    return std::isfinite(out);
}

static TupleValue parse_tuple_value(const std::string& raw) {
    TupleValue out;
    std::string s = trim(raw);

    if (s.empty()) {
        return out;
    }

    if (s.front() == '(' && s.back() == ')') {
        s = s.substr(1, s.size() - 2);
    }

    const std::vector<std::string> fields = split_csv_line(s);
    if (fields.empty()) {
        return out;
    }

    double value = 0.0;
    if (!parse_double(fields[0], value)) {
        return out;
    }

    double stat = 0.0;
    if (fields.size() >= 2) {
        if (!parse_double(fields[1], stat)) {
            return out;
        }
    } else {
        return out;
    }

    if (!std::isfinite(value) || !std::isfinite(stat) || stat <= 0.0) {
        return out;
    }

    out.ok = true;
    out.value = value;
    out.stat = stat;
    return out;
}

static double get_required_double(const CsvTable& table,
                                  const std::vector<std::string>& row,
                                  const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) {
        throw std::runtime_error("Missing required column: " + column);
    }

    double value = 0.0;
    if (!parse_double(row[(size_t)it->second], value)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return value;
}

static TupleValue get_tuple(const CsvTable& table,
                            const std::vector<std::string>& row,
                            const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) {
        throw std::runtime_error("Missing required tuple column: " + column);
    }

    return parse_tuple_value(row[(size_t)it->second]);
}

static std::string cross_section_column(const std::string& period,
                                        const std::string& helicity) {
    return "normed cross sections, ep->epg, exp, " + period + ", " + helicity;
}

static std::string avg_column(const std::string& variable_prefix,
                              const std::string& avg_label) {
    return variable_prefix + ", " + avg_label;
}

static std::vector<ConsistencyCase> consistency_cases() {
    return {
        {
            "10.6 GeV unpol",
            "10p6_GeV_unpol",
            "10.6 GeV",
            "unpol",
            {
                {"Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")},
                {"Fa18 Out", cross_section_column("Fa18 Out", "unpol")},
                {"Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")},
                {"Sp18 Out", cross_section_column("Sp18 Out", "unpol")}
            }
        },
        {
            "Fa18 unpol",
            "Fa18_unpol",
            "Fa18",
            "unpol",
            {
                {"Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")},
                {"Fa18 Out", cross_section_column("Fa18 Out", "unpol")}
            }
        },
        {
            "Fa18 pos",
            "Fa18_pos",
            "Fa18",
            "pos",
            {
                {"Fa18 Inb", cross_section_column("Fa18 Inb", "pos")},
                {"Fa18 Out", cross_section_column("Fa18 Out", "pos")}
            }
        },
        {
            "Fa18 neg",
            "Fa18_neg",
            "Fa18",
            "neg",
            {
                {"Fa18 Inb", cross_section_column("Fa18 Inb", "neg")},
                {"Fa18 Out", cross_section_column("Fa18 Out", "neg")}
            }
        },
        {
            "Sp18 unpol",
            "Sp18_unpol",
            "Sp18",
            "unpol",
            {
                {"Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")},
                {"Sp18 Out", cross_section_column("Sp18 Out", "unpol")}
            }
        }
    };
}

static std::vector<VariableConfig> variable_configs() {
    return {
        {"xB", "xBavg", "x_{B}"},
        {"Q2", "Q2avg", "Q^{2} (GeV^{2})"},
        {"t", "t_abs_avg", "-t (GeV^{2})"},
        {"phi", "phiavg", "#phi (deg)"}
    };
}

static void validate_schema(const CsvTable& table,
                            const std::vector<ConsistencyCase>& cases,
                            const std::vector<VariableConfig>& vars) {
    std::vector<std::string> required;

    for (const auto& c : cases) {
        for (const auto& input : c.inputs) {
            required.push_back(input.column);
        }

        for (const auto& var : vars) {
            required.push_back(avg_column(var.column_prefix, c.avg_label));
        }
    }

    require_columns(table, required, "run-period consistency systematics");
}

static bool compute_weighted_mean(const std::vector<TupleValue>& values,
                                  double& mean,
                                  double& mean_stat) {
    double sum_w = 0.0;
    double sum_wx = 0.0;

    for (const auto& v : values) {
        if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) {
            return false;
        }

        const double w = 1.0 / (v.stat * v.stat);
        sum_w += w;
        sum_wx += w * v.value;
    }

    if (sum_w <= 0.0) {
        return false;
    }

    mean = sum_wx / sum_w;
    mean_stat = 1.0 / std::sqrt(sum_w);
    return std::isfinite(mean) && std::isfinite(mean_stat) && mean_stat > 0.0;
}

static double compute_chi2_reduced(const std::vector<TupleValue>& values,
                                   double mean) {
    if (values.size() < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double chi2 = 0.0;
    for (const auto& v : values) {
        const double residual = v.value - mean;
        chi2 += residual * residual / (v.stat * v.stat);
    }

    const double ndf = (double)values.size() - 1.0;
    if (ndf <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return chi2 / ndf;
}

static std::vector<double> choose_x_range(const std::vector<double>& xs) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();

    for (const double x : xs) {
        if (!std::isfinite(x)) {
            continue;
        }
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
    }

    if (!std::isfinite(xmin) || !std::isfinite(xmax)) {
        return {0.0, 1.0};
    }

    if (xmin == xmax) {
        const double pad = (std::abs(xmin) > 0.0) ? 0.05 * std::abs(xmin) : 1.0;
        return {xmin - pad, xmax + pad};
    }

    const double pad = 0.07 * (xmax - xmin);
    return {xmin - pad, xmax + pad};
}

static std::vector<double> choose_y_range(const std::vector<PointPack>& packs,
                                          double floor_value,
                                          bool include_one) {
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    if (include_one) {
        ymin = std::min(ymin, 1.0);
        ymax = std::max(ymax, 1.0);
    }

    for (const auto& pack : packs) {
        for (size_t i = 0; i < pack.y.size(); ++i) {
            const double y = pack.y[i];
            if (!std::isfinite(y)) {
                continue;
            }

            const double ey = (i < pack.ey.size() && std::isfinite(pack.ey[i])) ? pack.ey[i] : 0.0;
            ymin = std::min(ymin, y - ey);
            ymax = std::max(ymax, y + ey);
        }
    }

    if (!std::isfinite(ymin) || !std::isfinite(ymax)) {
        return {floor_value, floor_value + 1.0};
    }

    ymin = std::min(ymin, floor_value);

    if (ymin == ymax) {
        const double pad = (std::abs(ymin) > 0.0) ? 0.20 * std::abs(ymin) : 1.0;
        return {ymin - pad, ymax + pad};
    }

    const double pad = 0.12 * (ymax - ymin);
    return {std::max(floor_value, ymin - pad), ymax + pad};
}

static void set_pad_style() {
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.12);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicks(1, 1);
}

static TGraphErrors* make_graph(const PointPack& pack,
                                int color,
                                int marker_style) {
    TGraphErrors* g = new TGraphErrors((int)pack.x.size());

    for (int i = 0; i < (int)pack.x.size(); ++i) {
        g->SetPoint(i, pack.x[(size_t)i], pack.y[(size_t)i]);
        const double ey = (i < (int)pack.ey.size()) ? pack.ey[(size_t)i] : 0.0;
        g->SetPointError(i, 0.0, ey);
    }

    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerStyle(marker_style);
    g->SetMarkerSize(0.65);
    g->SetLineWidth(1);
    return g;
}

static void draw_subpad_title(const std::string& text) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.070);
    latex.SetTextAlign(22);
    latex.DrawLatex(0.50, 0.965, text.c_str());
}

static void draw_chi2_canvas(const std::string& out_dir,
                             const ConsistencyCase& c,
                             const std::map<std::string, PointPack>& chi2_by_var) {
    TCanvas canvas(("c_chi2_" + c.file_tag).c_str(),
                   ("Run-period reduced chi2: " + c.label).c_str(),
                   1500,
                   1100);
    canvas.Divide(2, 2);

    const std::vector<VariableConfig> vars = variable_configs();

    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        canvas.cd(ivar + 1);
        set_pad_style();

        const auto it = chi2_by_var.find(vars[(size_t)ivar].key);
        PointPack pack;
        if (it != chi2_by_var.end()) {
            pack = it->second;
        }

        const std::vector<double> xr = choose_x_range(pack.x);
        std::vector<PointPack> ypacks = {pack};
        const std::vector<double> yr = choose_y_range(ypacks, 0.0, true);

        TH1D frame(("frame_chi2_" + c.file_tag + "_" + vars[(size_t)ivar].key).c_str(),
                   "",
                   100,
                   xr[0],
                   xr[1]);
        frame.SetMinimum(yr[0]);
        frame.SetMaximum(yr[1]);
        frame.GetXaxis()->SetTitle(vars[(size_t)ivar].title.c_str());
        frame.GetYaxis()->SetTitle("#chi^{2}/ndf");
        frame.GetXaxis()->CenterTitle();
        frame.GetYaxis()->CenterTitle();
        frame.GetXaxis()->SetTitleSize(0.060);
        frame.GetYaxis()->SetTitleSize(0.060);
        frame.GetXaxis()->SetLabelSize(0.050);
        frame.GetYaxis()->SetLabelSize(0.050);
        frame.Draw("AXIS");

        TLine line(xr[0], 1.0, xr[1], 1.0);
        line.SetLineStyle(2);
        line.SetLineWidth(1);
        line.SetLineColor(kRed + 1);
        line.Draw("SAME");

        if (!pack.x.empty()) {
            TGraphErrors* graph = make_graph(pack, kBlack, 20);
            graph->Draw("P SAME");
        }

        draw_subpad_title(vars[(size_t)ivar].title);
    }

    canvas.cd();
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.028);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.985, ("Run-period consistency: " + c.label).c_str());

    const std::string png = out_dir + "/" + c.file_tag + "_reduced_chi2.png";
    const std::string pdf = out_dir + "/" + c.file_tag + "_reduced_chi2.pdf";
    canvas.SaveAs(png.c_str());
    canvas.SaveAs(pdf.c_str());
}

static void draw_ratio_canvas(const std::string& out_dir,
                              const ConsistencyCase& c,
                              const std::map<std::string, std::map<std::string, PointPack> >& ratios_by_var_period) {
    TCanvas canvas(("c_ratio_" + c.file_tag).c_str(),
                   ("Run-period ratios: " + c.label).c_str(),
                   1500,
                   1100);
    canvas.Divide(2, 2);

    const std::vector<VariableConfig> vars = variable_configs();
    const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7};
    const std::vector<int> markers = {20, 24, 21, 25, 22, 26};

    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        canvas.cd(ivar + 1);
        set_pad_style();

        std::vector<double> all_x;
        std::vector<PointPack> ypacks;

        const auto it_var = ratios_by_var_period.find(vars[(size_t)ivar].key);
        if (it_var != ratios_by_var_period.end()) {
            for (const auto& input : c.inputs) {
                const auto it_pack = it_var->second.find(input.period);
                if (it_pack == it_var->second.end()) {
                    continue;
                }
                ypacks.push_back(it_pack->second);
                for (const double x : it_pack->second.x) {
                    all_x.push_back(x);
                }
            }
        }

        const std::vector<double> xr = choose_x_range(all_x);
        const std::vector<double> yr = choose_y_range(ypacks, 0.0, true);

        TH1D frame(("frame_ratio_" + c.file_tag + "_" + vars[(size_t)ivar].key).c_str(),
                   "",
                   100,
                   xr[0],
                   xr[1]);
        frame.SetMinimum(yr[0]);
        frame.SetMaximum(yr[1]);
        frame.GetXaxis()->SetTitle(vars[(size_t)ivar].title.c_str());
        frame.GetYaxis()->SetTitle("#sigma_{i}/#bar{#sigma}");
        frame.GetXaxis()->CenterTitle();
        frame.GetYaxis()->CenterTitle();
        frame.GetXaxis()->SetTitleSize(0.060);
        frame.GetYaxis()->SetTitleSize(0.060);
        frame.GetXaxis()->SetLabelSize(0.050);
        frame.GetYaxis()->SetLabelSize(0.050);
        frame.Draw("AXIS");

        TLine line(xr[0], 1.0, xr[1], 1.0);
        line.SetLineStyle(2);
        line.SetLineWidth(1);
        line.SetLineColor(kRed + 1);
        line.Draw("SAME");

        TLegend legend(0.56, 0.68, 0.93, 0.88);
        legend.SetBorderSize(1);
        legend.SetFillStyle(1001);
        legend.SetFillColor(kWhite);
        legend.SetTextSize(0.045);

        int iper = 0;
        if (it_var != ratios_by_var_period.end()) {
            for (const auto& input : c.inputs) {
                const auto it_pack = it_var->second.find(input.period);
                if (it_pack == it_var->second.end()) {
                    continue;
                }

                const PointPack& pack = it_pack->second;
                if (pack.x.empty()) {
                    ++iper;
                    continue;
                }

                const int color = colors[(size_t)(iper % (int)colors.size())];
                const int marker = markers[(size_t)(iper % (int)markers.size())];
                TGraphErrors* graph = make_graph(pack, color, marker);
                graph->Draw("P SAME");
                legend.AddEntry(graph, input.period.c_str(), "pe");
                ++iper;
            }
        }

        legend.Draw();
        draw_subpad_title(vars[(size_t)ivar].title);
    }

    canvas.cd();
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.028);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.985, ("Run-period ratios: " + c.label).c_str());

    const std::string png = out_dir + "/" + c.file_tag + "_period_ratios.png";
    const std::string pdf = out_dir + "/" + c.file_tag + "_period_ratios.pdf";
    canvas.SaveAs(png.c_str());
    canvas.SaveAs(pdf.c_str());
}

static void fill_case_outputs(const CsvTable& table,
                              const ConsistencyCase& c,
                              std::map<std::string, PointPack>& chi2_by_var,
                              std::map<std::string, std::map<std::string, PointPack> >& ratios_by_var_period) {
    const std::vector<VariableConfig> vars = variable_configs();

    for (const auto& row : table.rows) {
        std::vector<TupleValue> values;
        values.reserve(c.inputs.size());

        for (const auto& input : c.inputs) {
            values.push_back(get_tuple(table, row, input.column));
        }

        double mean = 0.0;
        double mean_stat = 0.0;
        if (!compute_weighted_mean(values, mean, mean_stat)) {
            continue;
        }

        if (!std::isfinite(mean) || std::abs(mean) <= 0.0) {
            continue;
        }

        const double chi2_ndf = compute_chi2_reduced(values, mean);
        if (!std::isfinite(chi2_ndf)) {
            continue;
        }

        for (const auto& var : vars) {
            const std::string col = avg_column(var.column_prefix, c.avg_label);
            const double x = get_required_double(table, row, col);
            if (!std::isfinite(x)) {
                continue;
            }

            chi2_by_var[var.key].x.push_back(x);
            chi2_by_var[var.key].y.push_back(chi2_ndf);
            chi2_by_var[var.key].ey.push_back(0.0);

            for (size_t iper = 0; iper < c.inputs.size(); ++iper) {
                const TupleValue& v = values[iper];
                const double ratio = v.value / mean;
                const double ratio_stat = std::abs(v.stat / mean);

                if (!std::isfinite(ratio) || !std::isfinite(ratio_stat)) {
                    continue;
                }

                PointPack& pack = ratios_by_var_period[var.key][c.inputs[iper].period];
                pack.x.push_back(x);
                pack.y.push_back(ratio);
                pack.ey.push_back(ratio_stat);
            }
        }
    }
}

} // namespace

bool run_period_consistency_systematics(const std::string& csv_path,
                                        const std::string& output_root_dir,
                                        bool make_all_cases) {
    try {
        const std::string out_dir = output_root_dir + "/run_period_consistency";
        fs::create_directories(out_dir);

        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");

        const CsvTable table = read_csv_or_throw(csv_path);
        const std::vector<ConsistencyCase> cases = consistency_cases();
        const std::vector<VariableConfig> vars = variable_configs();
        validate_schema(table, cases, vars);

        std::cout << "[run-period-consistency] CSV rows loaded: " << table.rows.size() << "\n";
        std::cout << "[run-period-consistency] Output directory: " << out_dir << "\n";

        for (const auto& c : cases) {
            std::map<std::string, PointPack> chi2_by_var;
            std::map<std::string, std::map<std::string, PointPack> > ratios_by_var_period;

            fill_case_outputs(table, c, chi2_by_var, ratios_by_var_period);
            draw_chi2_canvas(out_dir, c, chi2_by_var);
            draw_ratio_canvas(out_dir, c, ratios_by_var_period);

            size_t npoints = 0;
            const auto it = chi2_by_var.find("phi");
            if (it != chi2_by_var.end()) {
                npoints = it->second.x.size();
            }

            std::cout << "[run-period-consistency] Completed " << c.label
                      << " with " << npoints << " valid bins.\n";
        }

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[run-period-consistency] FATAL: " << e.what() << "\n";
        return false;
    }
}
