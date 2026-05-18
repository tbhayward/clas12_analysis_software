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

struct VariableGroup {
    std::string label;
    std::string file_tag;
    int ncols = 2;
    int nrows = 2;
    std::vector<VariableConfig> variables;
};

struct PeriodRatioAccumulator {
    double sum_w = 0.0;
    double sum_wr = 0.0;
    int n = 0;
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
        idx[trim(header[(size_t)i])] = i;
    }

    return idx;
}

static void assert_no_duplicate_columns(const std::vector<std::string>& header) {
    std::set<std::string> seen;
    std::set<std::string> duplicates;

    for (const auto& raw_name : header) {
        const std::string name = trim(raw_name);

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

static bool parse_double(const std::string& raw,
                         double& out) {
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

static bool parse_numeric_or_tuple_first(const std::string& raw,
                                         double& out) {
    if (parse_double(raw, out)) {
        return true;
    }

    std::string s = trim(raw);
    if (s.empty()) {
        return false;
    }

    if (s.front() == '(' && s.back() == ')') {
        s = trim(s.substr(1, s.size() - 2));
    }

    const std::vector<std::string> fields = split_csv_line(s);
    if (fields.empty()) {
        return false;
    }

    return parse_double(fields[0], out);
}

static TupleValue parse_tuple_value(const std::string& raw) {
    TupleValue out;
    std::string s = trim(raw);

    if (s.empty()) {
        return out;
    }

    if (s.front() == '(' && s.back() == ')') {
        s = trim(s.substr(1, s.size() - 2));
    }

    const std::vector<std::string> fields = split_csv_line(s);
    if (fields.size() < 2) {
        return out;
    }

    double value = 0.0;
    double stat = 0.0;

    if (!parse_double(fields[0], value)) {
        return out;
    }

    if (!parse_double(fields[1], stat)) {
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

static double get_numeric_or_nan(const CsvTable& table,
                                 const std::vector<std::string>& row,
                                 const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) {
        throw std::runtime_error("Missing required numeric column: " + column);
    }

    double value = 0.0;
    if (!parse_numeric_or_tuple_first(row[(size_t)it->second], value)) {
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

static ColumnInput column_input(const std::string& source_period,
                                const std::string& column) {
    ColumnInput out;
    out.source_period = source_period;
    out.column = column;
    return out;
}

static PeriodInput single_input(const std::string& label,
                                const std::string& source_period,
                                const std::string& column) {
    PeriodInput out;
    out.period = label;
    out.columns.push_back(column_input(source_period, column));
    return out;
}

static PeriodInput grouped_input(const std::string& label,
                                 const std::vector<ColumnInput>& columns) {
    PeriodInput out;
    out.period = label;
    out.columns = columns;
    return out;
}

static std::vector<ConsistencyCase> consistency_cases() {
    return {
        {
            "10.6 GeV unpol",
            "10p6_GeV_unpol",
            "10.6 GeV",
            "unpol",
            {
                single_input("Fa18 Inb", "Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                single_input("Fa18 Out", "Fa18 Out", cross_section_column("Fa18 Out", "unpol")),
                single_input("Sp18 Inb", "Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
                single_input("Sp18 Out", "Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
            }
        },
        {
            "Fa18 unpol",
            "Fa18_unpol",
            "Fa18",
            "unpol",
            {
                single_input("Fa18 Inb", "Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                single_input("Fa18 Out", "Fa18 Out", cross_section_column("Fa18 Out", "unpol"))
            }
        },
        {
            "Fa18 pos",
            "Fa18_pos",
            "Fa18",
            "pos",
            {
                single_input("Fa18 Inb", "Fa18 Inb", cross_section_column("Fa18 Inb", "pos")),
                single_input("Fa18 Out", "Fa18 Out", cross_section_column("Fa18 Out", "pos"))
            }
        },
        {
            "Fa18 neg",
            "Fa18_neg",
            "Fa18",
            "neg",
            {
                single_input("Fa18 Inb", "Fa18 Inb", cross_section_column("Fa18 Inb", "neg")),
                single_input("Fa18 Out", "Fa18 Out", cross_section_column("Fa18 Out", "neg"))
            }
        },
        {
            "Sp18 unpol",
            "Sp18_unpol",
            "Sp18",
            "unpol",
            {
                single_input("Sp18 Inb", "Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
                single_input("Sp18 Out", "Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
            }
        },
        {
            "Fa18 vs Sp18 unpol",
            "Fa18_vs_Sp18_unpol",
            "10.6 GeV",
            "unpol",
            {
                grouped_input(
                    "Fa18",
                    {
                        column_input("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                        column_input("Fa18 Out", cross_section_column("Fa18 Out", "unpol"))
                    }
                ),
                grouped_input(
                    "Sp18",
                    {
                        column_input("Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
                        column_input("Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
                    }
                )
            }
        },
        {
            "Inb vs Out unpol",
            "Inb_vs_Out_unpol",
            "10.6 GeV",
            "unpol",
            {
                grouped_input(
                    "Inb",
                    {
                        column_input("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                        column_input("Sp18 Inb", cross_section_column("Sp18 Inb", "unpol"))
                    }
                ),
                grouped_input(
                    "Out",
                    {
                        column_input("Fa18 Out", cross_section_column("Fa18 Out", "unpol")),
                        column_input("Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
                    }
                )
            }
        }
    };
}

static std::vector<VariableGroup> variable_groups() {
    return {
        {
            "physics kinematics",
            "physics_kinematics",
            2,
            2,
            {
                {"xB", "xBavg", "x_{B}"},
                {"Q2", "Q2avg", "Q^{2} (GeV^{2})"},
                {"t", "t_abs_avg", "-t (GeV^{2})"},
                {"phi", "phiavg", "#phi (deg)"}
            }
        },
        {
            "polar angles",
            "polar_angles",
            3,
            1,
            {
                {"e_theta", "e_theta", "#theta_{e} (deg)"},
                {"p_theta", "p_theta", "#theta_{p} (deg)"},
                {"g_theta", "g_theta", "#theta_{#gamma} (deg)"}
            }
        }
    };
}

static std::vector<VariableConfig> all_variable_configs() {
    std::vector<VariableConfig> out;

    for (const auto& group : variable_groups()) {
        for (const auto& var : group.variables) {
            out.push_back(var);
        }
    }

    return out;
}

static void validate_schema(const CsvTable& table,
                            const std::vector<ConsistencyCase>& cases) {
    std::vector<std::string> required;
    const std::vector<VariableConfig> vars = all_variable_configs();

    for (const auto& c : cases) {
        for (const auto& input : c.inputs) {
            for (const auto& col : input.columns) {
                required.push_back(col.column);
            }
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
            continue;
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

static TupleValue scaled_tuple(const TupleValue& input,
                               const std::string& source_period,
                               bool use_scaling,
                               const std::map<std::string, double>& scale_by_period) {
    if (!input.ok) {
        return input;
    }

    if (!use_scaling) {
        return input;
    }

    if (source_period.empty()) {
        return input;
    }

    const auto it = scale_by_period.find(source_period);
    if (it == scale_by_period.end()) {
        return input;
    }

    const double scale = it->second;
    if (!std::isfinite(scale) || std::abs(scale) <= 0.0) {
        return input;
    }

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
    values.reserve(input.columns.size());

    for (const auto& col : input.columns) {
        TupleValue v = get_tuple(table, row, col.column);
        v = scaled_tuple(v, col.source_period, use_scaling, scale_by_period);

        if (v.ok) {
            values.push_back(v);
        }
    }

    TupleValue out;

    double mean = 0.0;
    double mean_stat = 0.0;

    if (!compute_weighted_mean(values, mean, mean_stat)) {
        return out;
    }

    out.ok = true;
    out.value = mean;
    out.stat = mean_stat;

    return out;
}

static double compute_chi2_reduced(const std::vector<TupleValue>& values,
                                   double mean) {
    int n = 0;
    double chi2 = 0.0;

    for (const auto& v : values) {
        if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) {
            continue;
        }

        const double residual = v.value - mean;
        chi2 += residual * residual / (v.stat * v.stat);
        ++n;
    }

    const double ndf = (double)n - 1.0;
    if (ndf <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return chi2 / ndf;
}

static std::vector<std::string> base_periods() {
    return {
        "Fa18 Inb",
        "Fa18 Out",
        "Sp18 Inb",
        "Sp18 Out"
    };
}

static std::map<std::string, double>
compute_10p6_unpol_scale_factors(const CsvTable& table) {
    const std::vector<std::string> periods = base_periods();

    std::map<std::string, PeriodRatioAccumulator> acc_by_period;
    for (const auto& period : periods) {
        acc_by_period[period] = PeriodRatioAccumulator();
    }

    for (const auto& row : table.rows) {
        std::vector<TupleValue> values;
        values.reserve(periods.size());

        for (const auto& period : periods) {
            values.push_back(get_tuple(table, row, cross_section_column(period, "unpol")));
        }

        double mean = 0.0;
        double mean_stat = 0.0;

        if (!compute_weighted_mean(values, mean, mean_stat)) {
            continue;
        }

        if (!std::isfinite(mean) || std::abs(mean) <= 0.0) {
            continue;
        }

        for (size_t i = 0; i < periods.size(); ++i) {
            const TupleValue& v = values[i];

            if (!v.ok || !std::isfinite(v.value) || v.stat <= 0.0) {
                continue;
            }

            const double ratio = v.value / mean;
            const double ratio_stat = std::abs(v.stat / mean);

            if (!std::isfinite(ratio) || !std::isfinite(ratio_stat) || ratio_stat <= 0.0) {
                continue;
            }

            const double w = 1.0 / (ratio_stat * ratio_stat);

            PeriodRatioAccumulator& acc = acc_by_period[periods[i]];
            acc.sum_w += w;
            acc.sum_wr += w * ratio;
            acc.n += 1;
        }
    }

    std::map<std::string, double> scale_by_period;

    std::cout << "[run-period-consistency] 10.6 GeV unpol scale factors used for scaled plots:\n";

    for (const auto& period : periods) {
        const PeriodRatioAccumulator& acc = acc_by_period[period];

        if (acc.sum_w <= 0.0 || acc.n <= 0) {
            std::ostringstream msg;
            msg << "Could not compute 10.6 GeV unpol scale factor for period " << period;
            throw std::runtime_error(msg.str());
        }

        const double mean_ratio = acc.sum_wr / acc.sum_w;
        const double mean_ratio_stat = 1.0 / std::sqrt(acc.sum_w);

        if (!std::isfinite(mean_ratio) || std::abs(mean_ratio) <= 0.0) {
            std::ostringstream msg;
            msg << "Invalid 10.6 GeV unpol scale factor for period " << period;
            throw std::runtime_error(msg.str());
        }

        scale_by_period[period] = mean_ratio;

        std::cout << "  " << period
                  << " scale = " << std::setprecision(10) << mean_ratio
                  << " +/- " << mean_ratio_stat
                  << "   n=" << acc.n << "\n";
    }

    return scale_by_period;
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
                                          double reference_value,
                                          bool include_reference,
                                          bool use_errors) {
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    if (include_reference) {
        ymin = std::min(ymin, reference_value);
        ymax = std::max(ymax, reference_value);
    }

    for (const auto& pack : packs) {
        for (size_t i = 0; i < pack.y.size(); ++i) {
            const double y = pack.y[i];

            if (!std::isfinite(y)) {
                continue;
            }

            double ey = 0.0;
            if (use_errors && i < pack.ey.size() && std::isfinite(pack.ey[i])) {
                ey = std::abs(pack.ey[i]);
            }

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
    gPad->SetTopMargin(0.14);
    gPad->SetBottomMargin(0.17);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
}

static TGraphErrors* make_graph(const PointPack& pack,
                                int color,
                                int marker_style,
                                bool use_errors) {
    TGraphErrors* graph = new TGraphErrors((int)pack.x.size());

    for (int i = 0; i < (int)pack.x.size(); ++i) {
        graph->SetPoint(i, pack.x[(size_t)i], pack.y[(size_t)i]);

        const double ey =
            (use_errors && i < (int)pack.ey.size() && std::isfinite(pack.ey[(size_t)i]))
            ? pack.ey[(size_t)i]
            : 0.0;

        graph->SetPointError(i, 0.0, ey);
    }

    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetMarkerStyle(marker_style);
    graph->SetMarkerSize(0.55);
    graph->SetLineWidth(1);

    return graph;
}

static TH1D* make_frame(const std::string& name,
                        const std::string& xtitle,
                        const std::string& ytitle,
                        double xmin,
                        double xmax,
                        double ymin,
                        double ymax) {
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
    frame->GetXaxis()->SetLabelSize(0.050);
    frame->GetYaxis()->SetLabelSize(0.050);

    frame->GetXaxis()->SetTitleOffset(1.10);
    frame->GetYaxis()->SetTitleOffset(1.20);

    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(505);

    return frame;
}

static void draw_subpad_title(const std::string& text) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.064);
    latex.SetTextAlign(22);
    latex.DrawLatex(0.50, 0.950, text.c_str());
}

static void draw_top_title(TPad* title_pad,
                           const std::string& title_text) {
    title_pad->cd();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.28);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.55, title_text.c_str());
}

static void make_title_and_grid_pads(TCanvas& canvas,
                                     TPad*& title_pad,
                                     TPad*& grid_pad,
                                     int ncols,
                                     int nrows) {
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
    grid_pad->Divide(ncols, nrows, 0.001, 0.001);
}

static std::string scaled_tag(bool use_scaling) {
    return use_scaling ? "scaled" : "unscaled";
}

static std::string scaled_title(bool use_scaling) {
    return use_scaling ? "Scaled " : "Unscaled ";
}

static void draw_chi2_canvas(const fs::path& out_dir,
                             const ConsistencyCase& c,
                             const VariableGroup& group,
                             const std::map<std::string, PointPack>& chi2_by_var,
                             bool use_scaling) {
    const int width = 760 * group.ncols;
    const int height = 560 * group.nrows + 120;

    TCanvas canvas(("c_chi2_" + c.file_tag + "_" + group.file_tag + "_" + scaled_tag(use_scaling)).c_str(),
                   "c_chi2",
                   width,
                   height);
    canvas.SetFillColor(kWhite);

    TPad* title_pad = nullptr;
    TPad* grid_pad = nullptr;
    make_title_and_grid_pads(canvas, title_pad, grid_pad, group.ncols, group.nrows);
    draw_top_title(title_pad, scaled_title(use_scaling) + "run-period consistency: " + c.label + " (" + group.label + ")");

    std::vector<std::unique_ptr<TH1D> > frames;
    std::vector<std::unique_ptr<TGraphErrors> > graphs;
    std::vector<std::unique_ptr<TLine> > lines;

    for (int ivar = 0; ivar < (int)group.variables.size(); ++ivar) {
        grid_pad->cd(ivar + 1);
        set_pad_style();

        const VariableConfig& var = group.variables[(size_t)ivar];

        PointPack pack;
        const auto it = chi2_by_var.find(var.key);
        if (it != chi2_by_var.end()) {
            pack = it->second;
        }

        const std::vector<double> xr = choose_x_range(pack.x);
        const std::vector<PointPack> ypacks = {pack};
        const std::vector<double> yr = choose_y_range(ypacks, 0.0, 1.0, true, false);

        std::unique_ptr<TH1D> frame(
            make_frame("frame_chi2_" + c.file_tag + "_" + group.file_tag + "_" + scaled_tag(use_scaling) + "_" + var.key,
                       var.title,
                       "#chi^{2}/ndf",
                       xr[0],
                       xr[1],
                       yr[0],
                       yr[1])
        );
        frame->Draw("AXIS");

        std::unique_ptr<TLine> line(new TLine(xr[0], 1.0, xr[1], 1.0));
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kRed + 1);
        line->Draw("SAME");

        if (!pack.x.empty()) {
            std::unique_ptr<TGraphErrors> graph(make_graph(pack, kBlack, 20, false));
            graph->Draw("P SAME");
            graphs.push_back(std::move(graph));
        }

        draw_subpad_title(var.title);

        frames.push_back(std::move(frame));
        lines.push_back(std::move(line));
    }

    canvas.cd();
    canvas.Modified();
    canvas.Update();

    const fs::path png =
        out_dir / (scaled_tag(use_scaling) + "_reduced_chi2_" + group.file_tag + ".png");

    canvas.SaveAs(png.string().c_str());
}

static void draw_ratio_canvas(const fs::path& out_dir,
                              const ConsistencyCase& c,
                              const VariableGroup& group,
                              const std::map<std::string, std::map<std::string, PointPack> >& ratios_by_var_period,
                              bool use_scaling,
                              bool use_errors) {
    const int width = 760 * group.ncols;
    const int height = 560 * group.nrows + 120;

    const std::string error_tag = use_errors ? "_with_stat_errors" : "";

    TCanvas canvas(("c_ratio_" + c.file_tag + "_" + group.file_tag + "_" + scaled_tag(use_scaling) + error_tag).c_str(),
                   "c_ratio",
                   width,
                   height);
    canvas.SetFillColor(kWhite);

    TPad* title_pad = nullptr;
    TPad* grid_pad = nullptr;
    make_title_and_grid_pads(canvas, title_pad, grid_pad, group.ncols, group.nrows);

    draw_top_title(title_pad,
                   scaled_title(use_scaling) + "run-period ratios: " + c.label +
                   " (" + group.label + ")" + (use_errors ? " with stat errors" : ""));

    const std::vector<int> colors = {
        kBlack,
        kRed + 1,
        kBlue + 1,
        kGreen + 2,
        kMagenta + 1,
        kOrange + 7
    };

    const std::vector<int> markers = {
        20,
        24,
        21,
        25,
        22,
        26
    };

    std::vector<std::unique_ptr<TH1D> > frames;
    std::vector<std::unique_ptr<TGraphErrors> > graphs;
    std::vector<std::unique_ptr<TLine> > lines;
    std::vector<std::unique_ptr<TLegend> > legends;

    for (int ivar = 0; ivar < (int)group.variables.size(); ++ivar) {
        grid_pad->cd(ivar + 1);
        set_pad_style();

        const VariableConfig& var = group.variables[(size_t)ivar];

        std::vector<double> all_x;
        std::vector<PointPack> ypacks;

        const auto it_var = ratios_by_var_period.find(var.key);
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
        const std::vector<double> yr = choose_y_range(ypacks, 0.0, 1.0, true, use_errors);

        std::unique_ptr<TH1D> frame(
            make_frame("frame_ratio_" + c.file_tag + "_" + group.file_tag + "_" + scaled_tag(use_scaling) + error_tag + "_" + var.key,
                       var.title,
                       "#sigma_{i}/#bar{#sigma}",
                       xr[0],
                       xr[1],
                       yr[0],
                       yr[1])
        );
        frame->Draw("AXIS");

        std::unique_ptr<TLine> line(new TLine(xr[0], 1.0, xr[1], 1.0));
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kRed + 1);
        line->Draw("SAME");

        std::unique_ptr<TLegend> legend(new TLegend(0.55, 0.68, 0.93, 0.88));
        legend->SetBorderSize(1);
        legend->SetFillStyle(1001);
        legend->SetFillColor(kWhite);
        legend->SetTextFont(42);
        legend->SetTextSize(0.040);

        int iper = 0;

        if (it_var != ratios_by_var_period.end()) {
            for (const auto& input : c.inputs) {
                const auto it_pack = it_var->second.find(input.period);
                if (it_pack == it_var->second.end()) {
                    ++iper;
                    continue;
                }

                const PointPack& pack = it_pack->second;
                if (pack.x.empty()) {
                    ++iper;
                    continue;
                }

                const int color = colors[(size_t)(iper % (int)colors.size())];
                const int marker = markers[(size_t)(iper % (int)markers.size())];

                std::unique_ptr<TGraphErrors> graph(make_graph(pack, color, marker, use_errors));
                graph->Draw("P SAME");
                legend->AddEntry(graph.get(), input.period.c_str(), "p");
                graphs.push_back(std::move(graph));

                ++iper;
            }
        }

        legend->Draw();
        draw_subpad_title(var.title);

        frames.push_back(std::move(frame));
        lines.push_back(std::move(line));
        legends.push_back(std::move(legend));
    }

    canvas.cd();
    canvas.Modified();
    canvas.Update();

    const fs::path png =
        out_dir / (scaled_tag(use_scaling) + "_period_ratios_" + group.file_tag + error_tag + ".png");

    canvas.SaveAs(png.string().c_str());
}

static void fill_case_outputs(const CsvTable& table,
                              const ConsistencyCase& c,
                              const std::map<std::string, double>& scale_by_period,
                              bool use_scaling,
                              std::map<std::string, PointPack>& chi2_by_var,
                              std::map<std::string, std::map<std::string, PointPack> >& ratios_by_var_period) {
    const std::vector<VariableConfig> vars = all_variable_configs();

    for (const auto& row : table.rows) {
        std::vector<TupleValue> values;
        values.reserve(c.inputs.size());

        for (const auto& input : c.inputs) {
            const TupleValue v = evaluate_input_for_row(table, row, input, use_scaling, scale_by_period);
            if (v.ok) {
                values.push_back(v);
            } else {
                values.push_back(TupleValue());
            }
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
            const double x = get_numeric_or_nan(table, row, col);

            if (!std::isfinite(x)) {
                continue;
            }

            chi2_by_var[var.key].x.push_back(x);
            chi2_by_var[var.key].y.push_back(chi2_ndf);
            chi2_by_var[var.key].ey.push_back(0.0);

            for (size_t iper = 0; iper < c.inputs.size(); ++iper) {
                const TupleValue& v = values[iper];

                if (!v.ok || !std::isfinite(v.value) || v.stat <= 0.0) {
                    continue;
                }

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

static void make_plots_for_case(const fs::path& root_out_dir,
                                const CsvTable& table,
                                const ConsistencyCase& c,
                                const std::map<std::string, double>& scale_by_period,
                                bool use_scaling) {
    const fs::path out_dir = root_out_dir / c.file_tag;
    fs::create_directories(out_dir);

    std::map<std::string, PointPack> chi2_by_var;
    std::map<std::string, std::map<std::string, PointPack> > ratios_by_var_period;

    fill_case_outputs(table,
                      c,
                      scale_by_period,
                      use_scaling,
                      chi2_by_var,
                      ratios_by_var_period);

    for (const auto& group : variable_groups()) {
        draw_chi2_canvas(out_dir,
                         c,
                         group,
                         chi2_by_var,
                         use_scaling);

        draw_ratio_canvas(out_dir,
                          c,
                          group,
                          ratios_by_var_period,
                          use_scaling,
                          false);

        draw_ratio_canvas(out_dir,
                          c,
                          group,
                          ratios_by_var_period,
                          use_scaling,
                          true);
    }

    size_t npoints = 0;
    const auto it = chi2_by_var.find("phi");
    if (it != chi2_by_var.end()) {
        npoints = it->second.x.size();
    }

    std::cout << "[run-period-consistency] Completed "
              << (use_scaling ? "scaled " : "unscaled ")
              << c.label
              << " with " << npoints << " valid bins.\n";
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

        validate_schema(table, cases);

        const std::map<std::string, double> scale_by_period =
            compute_10p6_unpol_scale_factors(table);

        std::cout << "[run-period-consistency] CSV rows loaded: "
                  << table.rows.size() << "\n";

        std::cout << "[run-period-consistency] Output directory: "
                  << out_dir.string() << "\n";

        for (const auto& c : cases) {
            make_plots_for_case(out_dir,
                                table,
                                c,
                                scale_by_period,
                                false);

            make_plots_for_case(out_dir,
                                table,
                                c,
                                scale_by_period,
                                true);
        }

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[run-period-consistency] FATAL: " << e.what() << "\n";
        return false;
    }
}