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
    std::vector<PeriodInput> inputs;
};

struct PointPack {
    std::vector<double> x;
    std::vector<double> y;
};

struct VariableConfig {
    std::string key;
    std::string column_prefix;
    std::string title;
};

struct PeriodRatioAccumulator {
    double sum_w = 0.0;
    double sum_wr = 0.0;
    int n = 0;
};

static std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) {
        ++b;
    }

    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) {
        --e;
    }

    return s.substr(b, e - b);
}

static std::vector<std::string> split_csv_line(const std::string& line) {
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

static std::unordered_map<std::string, int>
build_header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> out;

    for (int i = 0; i < (int)header.size(); ++i) {
        out[trim(header[(size_t)i])] = i;
    }

    return out;
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
    if (fields.size() < 2U) {
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

static std::string cross_section_column(const std::string& label,
                                        const std::string& helicity) {
    return "normed cross sections, ep->epg, exp, " + label + ", " + helicity;
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

static std::vector<ConsistencyCase> all_consistency_cases() {
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

static std::vector<VariableConfig> physics_variable_configs() {
    return {
        {"xB",  "xBavg",      "x_{B}"},
        {"Q2",  "Q2avg",      "Q^{2} (GeV^{2})"},
        {"t",   "t_abs_avg",  "-t (GeV^{2})"},
        {"phi", "phiavg",     "#phi (deg)"}
    };
}

static std::vector<VariableConfig> angle_variable_configs() {
    return {
        {"e_theta", "e_theta", "#theta_{e} (deg)"},
        {"p_theta", "p_theta", "#theta_{p} (deg)"},
        {"g_theta", "g_theta", "#theta_{#gamma} (deg)"}
    };
}

static std::vector<VariableConfig> all_variable_configs() {
    std::vector<VariableConfig> out = physics_variable_configs();
    const std::vector<VariableConfig> angles = angle_variable_configs();

    for (const auto& var : angles) {
        out.push_back(var);
    }

    return out;
}

static void validate_schema(const CsvTable& table,
                            const std::vector<ConsistencyCase>& cases) {
    std::vector<std::string> required;

    for (const auto& c : cases) {
        for (const auto& input : c.inputs) {
            for (const auto& col : input.columns) {
                required.push_back(col.column);
            }
        }

        for (const auto& var : all_variable_configs()) {
            required.push_back(avg_column(var.column_prefix, c.avg_label));
        }
    }

    require_columns(table, required, "run-period consistency systematics");
}

static bool compute_weighted_mean_skip_invalid(const std::vector<TupleValue>& values,
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

static bool compute_leave_one_out_mean(const std::vector<TupleValue>& values,
                                       size_t leave_out_index,
                                       double& mean,
                                       double& mean_stat) {
    double sum_w = 0.0;
    double sum_wx = 0.0;

    for (size_t i = 0; i < values.size(); ++i) {
        if (i == leave_out_index) {
            continue;
        }

        const TupleValue& v = values[i];

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

    if (!compute_weighted_mean_skip_invalid(values, mean, mean_stat)) {
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

static double compute_leave_one_out_chi2_reduced(const std::vector<TupleValue>& values) {
    int n = 0;
    double chi2 = 0.0;

    for (size_t i = 0; i < values.size(); ++i) {
        const TupleValue& v = values[i];

        if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) {
            continue;
        }

        double loo_mean = 0.0;
        double loo_mean_stat = 0.0;

        if (!compute_leave_one_out_mean(values, i, loo_mean, loo_mean_stat)) {
            continue;
        }

        const double var = v.stat * v.stat + loo_mean_stat * loo_mean_stat;
        if (!std::isfinite(var) || var <= 0.0) {
            continue;
        }

        const double residual = v.value - loo_mean;
        chi2 += residual * residual / var;
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
            const TupleValue v = get_tuple(table, row, cross_section_column(period, "unpol"));
            values.push_back(v);
        }

        double mean = 0.0;
        double mean_stat = 0.0;

        if (!compute_weighted_mean_skip_invalid(values, mean, mean_stat)) {
            continue;
        }

        if (!std::isfinite(mean) || std::abs(mean) <= 0.0) {
            continue;
        }

        for (size_t iper = 0; iper < periods.size(); ++iper) {
            const TupleValue& v = values[iper];

            if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) {
                continue;
            }

            const double ratio = v.value / mean;
            const double ratio_stat = std::abs(v.stat / mean);

            if (!std::isfinite(ratio) || !std::isfinite(ratio_stat) || ratio_stat <= 0.0) {
                continue;
            }

            const double w = 1.0 / (ratio_stat * ratio_stat);

            PeriodRatioAccumulator& acc = acc_by_period[periods[iper]];
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
            throw std::runtime_error("Could not compute 10.6 GeV unpol scale factor for period " + period);
        }

        const double mean_ratio = acc.sum_wr / acc.sum_w;
        const double mean_ratio_stat = 1.0 / std::sqrt(acc.sum_w);

        if (!std::isfinite(mean_ratio) || std::abs(mean_ratio) <= 0.0) {
            throw std::runtime_error("Invalid 10.6 GeV unpol scale factor for period " + period);
        }

        scale_by_period[period] = mean_ratio;

        std::cout << "  " << std::left << std::setw(10) << period
                  << " scale = " << std::right << std::setprecision(10)
                  << mean_ratio << " +/- " << mean_ratio_stat
                  << "   n=" << acc.n << "\n";
    }

    return scale_by_period;
}

static double percentile(std::vector<double> values,
                         double q) {
    std::vector<double> clean;

    for (const double v : values) {
        if (std::isfinite(v)) {
            clean.push_back(v);
        }
    }

    if (clean.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::sort(clean.begin(), clean.end());

    if (q <= 0.0) {
        return clean.front();
    }

    if (q >= 1.0) {
        return clean.back();
    }

    const double pos = q * (double)(clean.size() - 1);
    const size_t lo = (size_t)std::floor(pos);
    const size_t hi = (size_t)std::ceil(pos);
    const double frac = pos - (double)lo;

    if (lo == hi) {
        return clean[lo];
    }

    return clean[lo] * (1.0 - frac) + clean[hi] * frac;
}

static std::vector<double> choose_x_range_from_packs(const std::vector<PointPack>& packs) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();

    for (const auto& pack : packs) {
        for (const double x : pack.x) {
            if (!std::isfinite(x)) {
                continue;
            }

            xmin = std::min(xmin, x);
            xmax = std::max(xmax, x);
        }
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

static std::vector<double> choose_chi2_y_range(const PointPack& pack) {
    std::vector<double> ys;

    for (const double y : pack.y) {
        if (std::isfinite(y) && y >= 0.0) {
            ys.push_back(y);
        }
    }

    if (ys.empty()) {
        return {0.0, 5.0};
    }

    double ymax = std::max(5.0, 1.25 * percentile(ys, 0.95));

    if (!std::isfinite(ymax) || ymax <= 0.0) {
        ymax = 5.0;
    }

    return {0.0, ymax};
}

static std::vector<double> choose_pull_y_range(const std::vector<PointPack>& packs) {
    std::vector<double> ys;

    for (const auto& pack : packs) {
        for (const double y : pack.y) {
            if (std::isfinite(y)) {
                ys.push_back(y);
            }
        }
    }

    if (ys.empty()) {
        return {-5.0, 5.0};
    }

    double ymin = std::min(-5.0, 1.25 * percentile(ys, 0.02));
    double ymax = std::max(5.0, 1.25 * percentile(ys, 0.98));

    if (!std::isfinite(ymin) || !std::isfinite(ymax) || ymin >= ymax) {
        return {-5.0, 5.0};
    }

    return {ymin, ymax};
}

static void set_pad_style() {
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.12);
    gPad->SetBottomMargin(0.16);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
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

    frame->GetXaxis()->SetLabelSize(0.047);
    frame->GetYaxis()->SetLabelSize(0.047);

    frame->GetXaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetTitleOffset(1.20);

    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(505);

    return frame;
}

static TGraph* make_marker_graph(const PointPack& pack,
                                 int color,
                                 int marker_style,
                                 double marker_size) {
    TGraph* graph = new TGraph((int)pack.x.size());

    for (int i = 0; i < (int)pack.x.size(); ++i) {
        graph->SetPoint(i, pack.x[(size_t)i], pack.y[(size_t)i]);
    }

    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetMarkerStyle(marker_style);
    graph->SetMarkerSize(marker_size);
    graph->SetLineWidth(1);

    return graph;
}

static void draw_subpad_label(const std::string& text) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.060);
    latex.SetTextAlign(22);
    latex.DrawLatex(0.50, 0.955, text.c_str());
}

static void draw_top_title(TPad* title_pad,
                           const std::string& title_text) {
    title_pad->cd();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.30);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.50, title_text.c_str());
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

static std::string plot_prefix(bool use_scaling) {
    return use_scaling ? "scaled" : "unscaled";
}

static std::string title_prefix(bool use_scaling) {
    return use_scaling ? "Scaled " : "Unscaled ";
}

static std::string canvas_tag_for_vars(const std::string& var_group_tag) {
    if (var_group_tag == "physics") {
        return "physics_kinematics";
    }

    if (var_group_tag == "angles") {
        return "polar_angles";
    }

    return var_group_tag;
}

static void draw_chi2_canvas(const fs::path& out_dir,
                             const ConsistencyCase& c,
                             const std::map<std::string, PointPack>& chi2_by_var,
                             const std::vector<VariableConfig>& vars,
                             const std::string& var_group_tag,
                             bool use_scaling,
                             const std::string& diagnostic_tag,
                             const std::string& diagnostic_title) {
    const int ncols = (var_group_tag == "angles") ? 3 : 2;
    const int nrows = (var_group_tag == "angles") ? 1 : 2;
    const int width = (var_group_tag == "angles") ? 1800 : 1600;
    const int height = (var_group_tag == "angles") ? 720 : 1200;

    const std::string prefix = plot_prefix(use_scaling);
    const std::string tag = canvas_tag_for_vars(var_group_tag);

    TCanvas canvas(("c_chi2_" + c.file_tag + "_" + prefix + "_" + diagnostic_tag + "_" + tag).c_str(),
                   ("Run-period reduced chi2: " + c.label).c_str(),
                   width,
                   height);
    canvas.SetFillColor(kWhite);

    TPad* title_pad = nullptr;
    TPad* grid_pad = nullptr;
    make_title_and_grid_pads(canvas, title_pad, grid_pad, ncols, nrows);
    draw_top_title(title_pad, title_prefix(use_scaling) + diagnostic_title + ": " + c.label + " (" + tag + ")");

    std::vector<std::unique_ptr<TH1D> > frames;
    std::vector<std::unique_ptr<TGraph> > graphs;
    std::vector<std::unique_ptr<TLine> > lines;

    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        grid_pad->cd(ivar + 1);
        set_pad_style();

        PointPack pack;
        const auto it = chi2_by_var.find(vars[(size_t)ivar].key);
        if (it != chi2_by_var.end()) {
            pack = it->second;
        }

        const std::vector<PointPack> xpacks = {pack};
        const std::vector<double> xr = choose_x_range_from_packs(xpacks);
        const std::vector<double> yr = choose_chi2_y_range(pack);

        std::unique_ptr<TH1D> frame(
            make_frame("frame_chi2_" + c.file_tag + "_" + prefix + "_" + diagnostic_tag + "_" + tag + "_" + vars[(size_t)ivar].key,
                       vars[(size_t)ivar].title,
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

    const fs::path out_path = out_dir / (prefix + "_" + diagnostic_tag + "_reduced_chi2_" + tag + ".png");
    canvas.SaveAs(out_path.string().c_str());
}

static void draw_multi_period_canvas(const fs::path& out_dir,
                                     const ConsistencyCase& c,
                                     const std::map<std::string, std::map<std::string, PointPack> >& by_var_period,
                                     const std::vector<VariableConfig>& vars,
                                     const std::string& var_group_tag,
                                     bool use_scaling,
                                     const std::string& kind,
                                     const std::string& ytitle,
                                     double ymin_fixed,
                                     double ymax_fixed,
                                     bool fixed_y,
                                     const std::string& diagnostic_title) {
    const int ncols = (var_group_tag == "angles") ? 3 : 2;
    const int nrows = (var_group_tag == "angles") ? 1 : 2;
    const int width = (var_group_tag == "angles") ? 1800 : 1600;
    const int height = (var_group_tag == "angles") ? 720 : 1200;

    const std::string prefix = plot_prefix(use_scaling);
    const std::string tag = canvas_tag_for_vars(var_group_tag);

    TCanvas canvas(("c_" + kind + "_" + c.file_tag + "_" + prefix + "_" + tag).c_str(),
                   ("Run-period " + kind + ": " + c.label).c_str(),
                   width,
                   height);
    canvas.SetFillColor(kWhite);

    TPad* title_pad = nullptr;
    TPad* grid_pad = nullptr;
    make_title_and_grid_pads(canvas, title_pad, grid_pad, ncols, nrows);
    draw_top_title(title_pad, title_prefix(use_scaling) + diagnostic_title + ": " + c.label + " (" + tag + ")");

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
    std::vector<std::unique_ptr<TGraph> > graphs;
    std::vector<std::unique_ptr<TLine> > lines;
    std::vector<std::unique_ptr<TLegend> > legends;

    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        grid_pad->cd(ivar + 1);
        set_pad_style();

        std::vector<PointPack> xpacks;
        std::vector<PointPack> ypacks;

        const auto it_var = by_var_period.find(vars[(size_t)ivar].key);
        if (it_var != by_var_period.end()) {
            for (const auto& input : c.inputs) {
                const auto it_pack = it_var->second.find(input.period);
                if (it_pack == it_var->second.end()) {
                    continue;
                }

                xpacks.push_back(it_pack->second);
                ypacks.push_back(it_pack->second);
            }
        }

        const std::vector<double> xr = choose_x_range_from_packs(xpacks);

        double ymin = ymin_fixed;
        double ymax = ymax_fixed;
        if (!fixed_y) {
            const std::vector<double> pr = choose_pull_y_range(ypacks);
            ymin = pr[0];
            ymax = pr[1];
        }

        std::unique_ptr<TH1D> frame(
            make_frame("frame_" + kind + "_" + c.file_tag + "_" + prefix + "_" + tag + "_" + vars[(size_t)ivar].key,
                       vars[(size_t)ivar].title,
                       ytitle,
                       xr[0],
                       xr[1],
                       ymin,
                       ymax)
        );
        frame->Draw("AXIS");

        const double ref_y = (kind.find("pulls") != std::string::npos) ? 0.0 : 1.0;

        std::unique_ptr<TLine> line(new TLine(xr[0], ref_y, xr[1], ref_y));
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

    const fs::path out_path = out_dir / (prefix + "_" + kind + "_" + tag + ".png");
    canvas.SaveAs(out_path.string().c_str());
}

static double mean_vec(const std::vector<double>& values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum = 0.0;
    for (const double v : values) {
        sum += v;
    }

    return sum / (double)values.size();
}

static double rms_vec(const std::vector<double>& values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double mean = mean_vec(values);
    double sum = 0.0;

    for (const double v : values) {
        const double d = v - mean;
        sum += d * d;
    }

    return std::sqrt(sum / (double)values.size());
}

static void print_summary_block(const std::string& label,
                                const std::map<std::string, PointPack>& chi2_by_var,
                                const std::map<std::string, std::map<std::string, PointPack> >& pulls_by_var_period) {
    std::cout << label << "\n";

    const std::vector<VariableConfig> vars = all_variable_configs();

    for (const auto& var : vars) {
        const auto it = chi2_by_var.find(var.key);
        if (it == chi2_by_var.end()) {
            continue;
        }

        std::cout << "    chi2/ndf vs " << std::left << std::setw(8) << var.key
                  << ": n=" << std::right << it->second.y.size()
                  << " mean=" << std::setprecision(8) << mean_vec(it->second.y)
                  << " rms=" << rms_vec(it->second.y) << "\n";
    }

    const auto it_phi = pulls_by_var_period.find("phi");
    if (it_phi != pulls_by_var_period.end()) {
        for (const auto& kv : it_phi->second) {
            std::cout << "    pull " << std::left << std::setw(10) << kv.first
                      << ": n=" << std::right << kv.second.y.size()
                      << " mean=" << std::setprecision(8) << mean_vec(kv.second.y)
                      << " rms=" << rms_vec(kv.second.y) << "\n";
        }
    }
}

static void fill_case_outputs(const CsvTable& table,
                              const ConsistencyCase& c,
                              const std::map<std::string, double>& scale_by_period,
                              bool use_scaling,
                              std::map<std::string, PointPack>& chi2_by_var,
                              std::map<std::string, PointPack>& loo_chi2_by_var,
                              std::map<std::string, std::map<std::string, PointPack> >& ratios_by_var_period,
                              std::map<std::string, std::map<std::string, PointPack> >& pulls_by_var_period,
                              std::map<std::string, std::map<std::string, PointPack> >& loo_ratios_by_var_period,
                              std::map<std::string, std::map<std::string, PointPack> >& loo_pulls_by_var_period) {
    const std::vector<VariableConfig> vars = all_variable_configs();

    for (const auto& row : table.rows) {
        std::vector<TupleValue> values;
        values.reserve(c.inputs.size());

        for (const auto& input : c.inputs) {
            const TupleValue v = evaluate_input_for_row(table, row, input, use_scaling, scale_by_period);
            values.push_back(v);
        }

        double mean = 0.0;
        double mean_stat = 0.0;

        if (!compute_weighted_mean_skip_invalid(values, mean, mean_stat)) {
            continue;
        }

        if (!std::isfinite(mean) || std::abs(mean) <= 0.0) {
            continue;
        }

        const double chi2_ndf = compute_chi2_reduced(values, mean);
        const double loo_chi2_ndf = compute_leave_one_out_chi2_reduced(values);

        if (!std::isfinite(chi2_ndf) && !std::isfinite(loo_chi2_ndf)) {
            continue;
        }

        for (const auto& var : vars) {
            const std::string col = avg_column(var.column_prefix, c.avg_label);
            const double x = get_numeric_or_nan(table, row, col);

            if (!std::isfinite(x)) {
                continue;
            }

            if (std::isfinite(chi2_ndf)) {
                chi2_by_var[var.key].x.push_back(x);
                chi2_by_var[var.key].y.push_back(chi2_ndf);
            }

            if (std::isfinite(loo_chi2_ndf)) {
                loo_chi2_by_var[var.key].x.push_back(x);
                loo_chi2_by_var[var.key].y.push_back(loo_chi2_ndf);
            }

            for (size_t iper = 0; iper < c.inputs.size(); ++iper) {
                const TupleValue& v = values[iper];

                if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) {
                    continue;
                }

                const double ratio = v.value / mean;
                const double pull = (v.value - mean) / v.stat;

                if (std::isfinite(ratio)) {
                    PointPack& pack = ratios_by_var_period[var.key][c.inputs[iper].period];
                    pack.x.push_back(x);
                    pack.y.push_back(ratio);
                }

                if (std::isfinite(pull)) {
                    PointPack& pack = pulls_by_var_period[var.key][c.inputs[iper].period];
                    pack.x.push_back(x);
                    pack.y.push_back(pull);
                }

                double loo_mean = 0.0;
                double loo_mean_stat = 0.0;

                if (!compute_leave_one_out_mean(values, iper, loo_mean, loo_mean_stat)) {
                    continue;
                }

                if (!std::isfinite(loo_mean) || std::abs(loo_mean) <= 0.0) {
                    continue;
                }

                const double loo_ratio = v.value / loo_mean;
                const double loo_var = v.stat * v.stat + loo_mean_stat * loo_mean_stat;
                const double loo_pull = (loo_var > 0.0 && std::isfinite(loo_var))
                                      ? (v.value - loo_mean) / std::sqrt(loo_var)
                                      : std::numeric_limits<double>::quiet_NaN();

                if (std::isfinite(loo_ratio)) {
                    PointPack& pack = loo_ratios_by_var_period[var.key][c.inputs[iper].period];
                    pack.x.push_back(x);
                    pack.y.push_back(loo_ratio);
                }

                if (std::isfinite(loo_pull)) {
                    PointPack& pack = loo_pulls_by_var_period[var.key][c.inputs[iper].period];
                    pack.x.push_back(x);
                    pack.y.push_back(loo_pull);
                }
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
    std::map<std::string, PointPack> loo_chi2_by_var;

    std::map<std::string, std::map<std::string, PointPack> > ratios_by_var_period;
    std::map<std::string, std::map<std::string, PointPack> > pulls_by_var_period;

    std::map<std::string, std::map<std::string, PointPack> > loo_ratios_by_var_period;
    std::map<std::string, std::map<std::string, PointPack> > loo_pulls_by_var_period;

    fill_case_outputs(table,
                      c,
                      scale_by_period,
                      use_scaling,
                      chi2_by_var,
                      loo_chi2_by_var,
                      ratios_by_var_period,
                      pulls_by_var_period,
                      loo_ratios_by_var_period,
                      loo_pulls_by_var_period);

    const std::vector<VariableConfig> physics_vars = physics_variable_configs();
    const std::vector<VariableConfig> angle_vars = angle_variable_configs();

    draw_chi2_canvas(out_dir,
                     c,
                     chi2_by_var,
                     physics_vars,
                     "physics",
                     use_scaling,
                     "all_mean",
                     "all-mean #chi^{2}/ndf");

    draw_chi2_canvas(out_dir,
                     c,
                     chi2_by_var,
                     angle_vars,
                     "angles",
                     use_scaling,
                     "all_mean",
                     "all-mean #chi^{2}/ndf");

    draw_chi2_canvas(out_dir,
                     c,
                     loo_chi2_by_var,
                     physics_vars,
                     "physics",
                     use_scaling,
                     "leave_one_out",
                     "leave-one-out #chi^{2}/ndf");

    draw_chi2_canvas(out_dir,
                     c,
                     loo_chi2_by_var,
                     angle_vars,
                     "angles",
                     use_scaling,
                     "leave_one_out",
                     "leave-one-out #chi^{2}/ndf");

    draw_multi_period_canvas(out_dir,
                             c,
                             ratios_by_var_period,
                             physics_vars,
                             "physics",
                             use_scaling,
                             "all_mean_period_ratios",
                             "#sigma_{i}/#bar{#sigma}_{all}",
                             0.0,
                             3.0,
                             true,
                             "all-mean period ratios");

    draw_multi_period_canvas(out_dir,
                             c,
                             ratios_by_var_period,
                             angle_vars,
                             "angles",
                             use_scaling,
                             "all_mean_period_ratios",
                             "#sigma_{i}/#bar{#sigma}_{all}",
                             0.0,
                             3.0,
                             true,
                             "all-mean period ratios");

    draw_multi_period_canvas(out_dir,
                             c,
                             loo_ratios_by_var_period,
                             physics_vars,
                             "physics",
                             use_scaling,
                             "leave_one_out_period_ratios",
                             "#sigma_{i}/#bar{#sigma}_{-i}",
                             0.0,
                             3.0,
                             true,
                             "leave-one-out period ratios");

    draw_multi_period_canvas(out_dir,
                             c,
                             loo_ratios_by_var_period,
                             angle_vars,
                             "angles",
                             use_scaling,
                             "leave_one_out_period_ratios",
                             "#sigma_{i}/#bar{#sigma}_{-i}",
                             0.0,
                             3.0,
                             true,
                             "leave-one-out period ratios");

    draw_multi_period_canvas(out_dir,
                             c,
                             pulls_by_var_period,
                             physics_vars,
                             "physics",
                             use_scaling,
                             "all_mean_pulls",
                             "(#sigma_{i}-#bar{#sigma}_{all})/#delta#sigma_{i}",
                             -5.0,
                             5.0,
                             false,
                             "all-mean pulls");

    draw_multi_period_canvas(out_dir,
                             c,
                             pulls_by_var_period,
                             angle_vars,
                             "angles",
                             use_scaling,
                             "all_mean_pulls",
                             "(#sigma_{i}-#bar{#sigma}_{all})/#delta#sigma_{i}",
                             -5.0,
                             5.0,
                             false,
                             "all-mean pulls");

    draw_multi_period_canvas(out_dir,
                             c,
                             loo_pulls_by_var_period,
                             physics_vars,
                             "physics",
                             use_scaling,
                             "leave_one_out_pulls",
                             "(#sigma_{i}-#bar{#sigma}_{-i})/#sqrt{#delta#sigma_{i}^{2}+#delta#bar{#sigma}_{-i}^{2}}",
                             -5.0,
                             5.0,
                             false,
                             "leave-one-out pulls");

    draw_multi_period_canvas(out_dir,
                             c,
                             loo_pulls_by_var_period,
                             angle_vars,
                             "angles",
                             use_scaling,
                             "leave_one_out_pulls",
                             "(#sigma_{i}-#bar{#sigma}_{-i})/#sqrt{#delta#sigma_{i}^{2}+#delta#bar{#sigma}_{-i}^{2}}",
                             -5.0,
                             5.0,
                             false,
                             "leave-one-out pulls");

    print_summary_block(
        "[run-period-consistency] All-mean summary for " +
        std::string(use_scaling ? "scaled " : "unscaled ") + c.label,
        chi2_by_var,
        pulls_by_var_period
    );

    print_summary_block(
        "[run-period-consistency] Leave-one-out summary for " +
        std::string(use_scaling ? "scaled " : "unscaled ") + c.label,
        loo_chi2_by_var,
        loo_pulls_by_var_period
    );

    size_t npoints = 0;
    const auto it = chi2_by_var.find("phi");
    if (it != chi2_by_var.end()) {
        npoints = it->second.x.size();
    }

    size_t loo_npoints = 0;
    const auto it_loo = loo_chi2_by_var.find("phi");
    if (it_loo != loo_chi2_by_var.end()) {
        loo_npoints = it_loo->second.x.size();
    }

    std::cout << "[run-period-consistency] Completed "
              << (use_scaling ? "scaled " : "unscaled ")
              << c.label
              << " with " << npoints
              << " all-mean valid bins and "
              << loo_npoints
              << " leave-one-out valid bins.\n";
}

} // namespace

bool run_period_consistency_systematics(const std::string& csv_path,
                                        const std::string& output_root_dir,
                                        bool make_all_cases) {
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
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);

        const CsvTable table = read_csv_or_throw(csv_path);
        const std::vector<ConsistencyCase> all_cases = all_consistency_cases();

        std::vector<ConsistencyCase> cases_to_run;
        for (const auto& c : all_cases) {
            if (make_all_cases || c.file_tag == "10p6_GeV_unpol") {
                cases_to_run.push_back(c);
            }
        }

        validate_schema(table, cases_to_run);

        const std::map<std::string, double> scale_by_period =
            compute_10p6_unpol_scale_factors(table);

        std::cout << "[run-period-consistency] CSV rows loaded: "
                  << table.rows.size() << "\n";
        std::cout << "[run-period-consistency] Output directory: "
                  << out_dir.string() << "\n";
        std::cout << "[run-period-consistency] make_all_cases = "
                  << (make_all_cases ? "true" : "false") << "\n";

        for (const auto& c : cases_to_run) {
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
        std::cerr << "[run-period-consistency] FATAL: "
                  << e.what() << "\n";
        return false;
    }
}