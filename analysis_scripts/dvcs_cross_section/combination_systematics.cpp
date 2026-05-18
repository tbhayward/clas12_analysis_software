#include "combination_systematics.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

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

enum class CentralValueMode {
    StatWeighted,
    EqualPeriod,
    RandomEffects
};

static constexpr CentralValueMode kFillOutputMode = CentralValueMode::StatWeighted;
static constexpr int kMaxPolynomialOrder = 5;

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string> > rows;
    std::unordered_map<std::string, int> index;
};

struct TupleValue {
    bool ok = false;
    double value = 0.0;
    double stat = 0.0;
    double sys = 0.0;
};

struct PeriodInput {
    std::string period;
    std::vector<std::string> columns;
};

struct CombinationCase {
    std::string label;
    std::string output_column;
    bool fill_output = true;
    std::vector<PeriodInput> inputs;
};

struct PeriodRatioAccumulator {
    double sum_w = 0.0;
    double sum_wr = 0.0;
    int n = 0;
};

struct PeriodRatioSummary {
    std::string period;

    int n = 0;
    double mean_ratio = 0.0;
    double mean_ratio_stat = 0.0;

    int loo_n = 0;
    double loo_mean_ratio = 0.0;
    double loo_mean_ratio_stat = 0.0;
};

struct CombinationResult {
    std::string label;
    std::string central_value_mode;
    std::string output_column;
    bool fill_output = false;

    int valid_bins = 0;
    int ratio_points = 0;
    int loo_ratio_points = 0;

    double s_obs = 0.0;
    double s_stat_exp = 0.0;
    double s_comb = 0.0;

    double loo_s_obs = 0.0;
    double loo_s_stat_exp = 0.0;
    double loo_s_comb = 0.0;

    std::vector<PeriodRatioSummary> period_summaries;
};

struct VariableConfig {
    std::string key;
    std::string column_prefix;
    std::string title;
};

struct RatioPoint {
    std::string case_label;
    std::string central_mode;
    std::string reference_type;
    std::string period;
    std::string variable_key;

    double x = 0.0;
    double y = 0.0;
    double ey = 0.0;
};

struct FitResultSummary {
    std::string case_label;
    std::string central_mode;
    std::string reference_type;
    std::string period;
    std::string variable_key;
    int order = 0;
    int n = 0;
    int ndf = 0;
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    double chi2_ndf = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> params;
    std::vector<double> errors;
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

static std::string sanitize_for_path(const std::string& s) {
    std::string out;
    out.reserve(s.size());

    for (const char c : s) {
        if (std::isalnum((unsigned char)c)) {
            out.push_back(c);
        } else if (c == '+' || c == '-') {
            out.push_back(c);
        } else {
            out.push_back('_');
        }
    }

    while (out.find("__") != std::string::npos) {
        const size_t pos = out.find("__");
        out.replace(pos, 2, "_");
    }

    if (!out.empty() && out.front() == '_') {
        out.erase(out.begin());
    }

    if (!out.empty() && out.back() == '_') {
        out.pop_back();
    }

    if (out.empty()) {
        out = "unnamed";
    }

    return out;
}

static std::string central_value_mode_name(CentralValueMode mode) {
    if (mode == CentralValueMode::StatWeighted) {
        return "stat_weighted";
    }

    if (mode == CentralValueMode::EqualPeriod) {
        return "equal_period";
    }

    if (mode == CentralValueMode::RandomEffects) {
        return "random_effects";
    }

    return "unknown";
}

static std::vector<CentralValueMode> central_value_modes() {
    return {
        CentralValueMode::StatWeighted,
        CentralValueMode::EqualPeriod,
        CentralValueMode::RandomEffects
    };
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

static std::string join_csv_row(const std::vector<std::string>& fields) {
    std::ostringstream oss;

    for (size_t i = 0; i < fields.size(); ++i) {
        const std::string& s = fields[i];

        const bool quote =
            s.find(',') != std::string::npos ||
            s.find('"') != std::string::npos ||
            s.find('\n') != std::string::npos ||
            s.find('\r') != std::string::npos;

        if (quote) {
            oss << '"';

            for (const char ch : s) {
                if (ch == '"') {
                    oss << "\"\"";
                } else {
                    oss << ch;
                }
            }

            oss << '"';
        } else {
            oss << s;
        }

        if (i + 1 < fields.size()) {
            oss << ',';
        }
    }

    return oss.str();
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

    std::set<std::string> seen;
    std::set<std::string> dup;

    for (int i = 0; i < (int)table.header.size(); ++i) {
        const std::string name = trim(table.header[(size_t)i]);

        if (!seen.insert(name).second) {
            dup.insert(name);
        }

        table.index[name] = i;
    }

    if (!dup.empty()) {
        std::ostringstream msg;
        msg << "CSV header contains duplicate columns:";

        for (const auto& d : dup) {
            msg << "\n  - " << d;
        }

        throw std::runtime_error(msg.str());
    }

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
                << " columns, but header has "
                << table.header.size()
                << " columns. File: " << path;

            throw std::runtime_error(msg.str());
        }

        table.rows.push_back(std::move(row));
    }

    return table;
}

static void write_csv_or_throw(const std::string& path,
                               const CsvTable& table) {
    const std::string tmp_path = path + ".tmp";
    std::ofstream fout(tmp_path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open temporary CSV for writing: " + tmp_path);
    }

    fout << join_csv_row(table.header) << "\n";

    for (const auto& row : table.rows) {
        if (row.size() != table.header.size()) {
            throw std::runtime_error("Internal CSV row/header size mismatch.");
        }

        fout << join_csv_row(row) << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing temporary CSV: " + tmp_path);
    }

    fs::rename(tmp_path, path);
}

static void require_columns(const CsvTable& table,
                            const std::vector<std::string>& required,
                            const std::string& context) {
    std::vector<std::string> missing;

    for (const auto& col : required) {
        if (table.index.find(col) == table.index.end()) {
            missing.push_back(col);
        }
    }

    if (!missing.empty()) {
        std::ostringstream msg;
        msg << "Missing required columns for " << context << ":";

        for (const auto& col : missing) {
            msg << "\n  - " << col;
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
    double sys = 0.0;

    if (!parse_double(fields[0], value)) {
        return out;
    }

    if (!parse_double(fields[1], stat)) {
        return out;
    }

    if (fields.size() > 2U) {
        parse_double(fields[2], sys);
    }

    if (!std::isfinite(value) || !std::isfinite(stat) || stat <= 0.0) {
        return out;
    }

    out.ok = true;
    out.value = value;
    out.stat = stat;
    out.sys = sys;
    return out;
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

static std::vector<TupleValue> valid_values_only(const std::vector<TupleValue>& values) {
    std::vector<TupleValue> out;

    for (const auto& v : values) {
        if (v.ok &&
            std::isfinite(v.value) &&
            std::isfinite(v.stat) &&
            v.stat > 0.0) {
            out.push_back(v);
        }
    }

    return out;
}

static bool combine_stat_weighted(const std::vector<TupleValue>& values,
                                  TupleValue& out) {
    const std::vector<TupleValue> valid = valid_values_only(values);

    if (valid.empty()) {
        return false;
    }

    double sum_w = 0.0;
    double sum_wx = 0.0;

    for (const auto& v : valid) {
        const double w = 1.0 / (v.stat * v.stat);
        sum_w += w;
        sum_wx += w * v.value;
    }

    if (sum_w <= 0.0) {
        return false;
    }

    out.ok = true;
    out.value = sum_wx / sum_w;
    out.stat = 1.0 / std::sqrt(sum_w);
    out.sys = 0.0;

    return std::isfinite(out.value) &&
           std::isfinite(out.stat) &&
           out.stat > 0.0;
}

static bool combine_equal_period(const std::vector<TupleValue>& values,
                                 TupleValue& out) {
    const std::vector<TupleValue> valid = valid_values_only(values);

    if (valid.empty()) {
        return false;
    }

    double sum_x = 0.0;
    double sum_stat2 = 0.0;

    for (const auto& v : valid) {
        sum_x += v.value;
        sum_stat2 += v.stat * v.stat;
    }

    const double n = (double)valid.size();

    out.ok = true;
    out.value = sum_x / n;
    out.stat = std::sqrt(sum_stat2) / n;
    out.sys = 0.0;

    return std::isfinite(out.value) &&
           std::isfinite(out.stat) &&
           out.stat > 0.0;
}

static double dersimonian_laird_tau2(const std::vector<TupleValue>& valid) {
    if (valid.size() < 2U) {
        return 0.0;
    }

    double sum_w = 0.0;
    double sum_w2 = 0.0;
    double sum_wx = 0.0;

    for (const auto& v : valid) {
        const double w = 1.0 / (v.stat * v.stat);
        sum_w += w;
        sum_w2 += w * w;
        sum_wx += w * v.value;
    }

    if (sum_w <= 0.0) {
        return 0.0;
    }

    const double mean = sum_wx / sum_w;

    double q = 0.0;
    for (const auto& v : valid) {
        const double w = 1.0 / (v.stat * v.stat);
        const double residual = v.value - mean;
        q += w * residual * residual;
    }

    const double ndf = (double)valid.size() - 1.0;
    const double c = sum_w - (sum_w2 / sum_w);

    if (c <= 0.0) {
        return 0.0;
    }

    const double tau2 = (q - ndf) / c;

    if (!std::isfinite(tau2) || tau2 <= 0.0) {
        return 0.0;
    }

    return tau2;
}

static bool combine_random_effects(const std::vector<TupleValue>& values,
                                   TupleValue& out) {
    const std::vector<TupleValue> valid = valid_values_only(values);

    if (valid.empty()) {
        return false;
    }

    const double tau2 = dersimonian_laird_tau2(valid);

    double sum_w = 0.0;
    double sum_wx = 0.0;

    for (const auto& v : valid) {
        const double var = v.stat * v.stat + tau2;

        if (!std::isfinite(var) || var <= 0.0) {
            continue;
        }

        const double w = 1.0 / var;
        sum_w += w;
        sum_wx += w * v.value;
    }

    if (sum_w <= 0.0) {
        return false;
    }

    out.ok = true;
    out.value = sum_wx / sum_w;
    out.stat = 1.0 / std::sqrt(sum_w);
    out.sys = std::sqrt(std::max(0.0, tau2));

    return std::isfinite(out.value) &&
           std::isfinite(out.stat) &&
           out.stat > 0.0;
}

static bool combine_values(const std::vector<TupleValue>& values,
                           CentralValueMode mode,
                           TupleValue& out) {
    if (mode == CentralValueMode::StatWeighted) {
        return combine_stat_weighted(values, out);
    }

    if (mode == CentralValueMode::EqualPeriod) {
        return combine_equal_period(values, out);
    }

    if (mode == CentralValueMode::RandomEffects) {
        return combine_random_effects(values, out);
    }

    return false;
}

static TupleValue get_input_value(const CsvTable& table,
                                  const std::vector<std::string>& row,
                                  const PeriodInput& input,
                                  CentralValueMode mode) {
    std::vector<TupleValue> values;
    values.reserve(input.columns.size());

    for (const auto& col : input.columns) {
        values.push_back(get_tuple(table, row, col));
    }

    TupleValue out;

    if (!combine_values(values, mode, out)) {
        return TupleValue();
    }

    return out;
}

static std::string cross_section_column(const std::string& label,
                                        const std::string& helicity) {
    return "normed cross sections, ep->epg, exp, " + label + ", " + helicity;
}

static std::string combination_column(const std::string& label,
                                      const std::string& helicity) {
    return cross_section_column(label, helicity) + ", combination sys";
}

static std::string avg_column(const std::string& variable_prefix,
                              const std::string& avg_label) {
    return variable_prefix + ", " + avg_label;
}

static PeriodInput input_single(const std::string& label,
                                const std::string& col) {
    PeriodInput out;
    out.period = label;
    out.columns.push_back(col);
    return out;
}

static PeriodInput input_group(const std::string& label,
                               const std::vector<std::string>& cols) {
    PeriodInput out;
    out.period = label;
    out.columns = cols;
    return out;
}

static std::vector<CombinationCase> combination_cases() {
    return {
        {
            "10.6 GeV unpol",
            combination_column("10.6 GeV", "unpol"),
            true,
            {
                input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                input_single("Fa18 Out", cross_section_column("Fa18 Out", "unpol")),
                input_single("Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
                input_single("Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
            }
        },
        {
            "Fa18 unpol",
            combination_column("Fa18", "unpol"),
            true,
            {
                input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
                input_single("Fa18 Out", cross_section_column("Fa18 Out", "unpol"))
            }
        },
        {
            "Fa18 pos",
            combination_column("Fa18", "pos"),
            true,
            {
                input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "pos")),
                input_single("Fa18 Out", cross_section_column("Fa18 Out", "pos"))
            }
        },
        {
            "Fa18 neg",
            combination_column("Fa18", "neg"),
            true,
            {
                input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "neg")),
                input_single("Fa18 Out", cross_section_column("Fa18 Out", "neg"))
            }
        },
        {
            "Sp18 unpol",
            combination_column("Sp18", "unpol"),
            true,
            {
                input_single("Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
                input_single("Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
            }
        },
        {
            "Fa18 vs Sp18 unpol",
            "",
            false,
            {
                input_group("Fa18", {
                    cross_section_column("Fa18 Inb", "unpol"),
                    cross_section_column("Fa18 Out", "unpol")
                }),
                input_group("Sp18", {
                    cross_section_column("Sp18 Inb", "unpol"),
                    cross_section_column("Sp18 Out", "unpol")
                })
            }
        },
        {
            "Inb vs Out unpol",
            "",
            false,
            {
                input_group("Inb", {
                    cross_section_column("Fa18 Inb", "unpol"),
                    cross_section_column("Sp18 Inb", "unpol")
                }),
                input_group("Out", {
                    cross_section_column("Fa18 Out", "unpol"),
                    cross_section_column("Sp18 Out", "unpol")
                })
            }
        }
    };
}

static std::vector<std::string> sp19_proxy_output_columns() {
    return {
        combination_column("Sp19 Inb", "unpol")
    };
}

static std::vector<VariableConfig> fit_variable_configs() {
    return {
        {"xB",      "xBavg",     "x_{B}"},
        {"Q2",      "Q2avg",     "Q^{2} (GeV^{2})"},
        {"t",       "t_abs_avg", "-t (GeV^{2})"},
        {"phi",     "phiavg",    "#phi (deg)"},
        {"e_theta", "e_theta",   "#theta_{e} (deg)"},
        {"p_theta", "p_theta",   "#theta_{p} (deg)"},
        {"g_theta", "g_theta",   "#theta_{#gamma} (deg)"}
    };
}

static void validate_schema(const CsvTable& table,
                            const std::vector<CombinationCase>& cases) {
    std::vector<std::string> required;

    for (const auto& c : cases) {
        if (c.fill_output) {
            required.push_back(c.output_column);
        }

        for (const auto& input : c.inputs) {
            for (const auto& col : input.columns) {
                required.push_back(col);
            }
        }

        for (const auto& var : fit_variable_configs()) {
            required.push_back(avg_column(var.column_prefix, "10.6 GeV"));
        }
    }

    for (const auto& col : sp19_proxy_output_columns()) {
        required.push_back(col);
    }

    require_columns(table, required, "combination systematics");
}

static std::string format_double(double value) {
    if (!std::isfinite(value)) {
        return "";
    }

    std::ostringstream oss;
    oss << std::setprecision(10) << value;
    return oss.str();
}

static std::string csv_escape_field(const std::string& s) {
    const bool need_quotes =
        s.find(',') != std::string::npos ||
        s.find('"') != std::string::npos ||
        s.find('\n') != std::string::npos ||
        s.find('\r') != std::string::npos;

    if (!need_quotes) {
        return s;
    }

    std::ostringstream oss;
    oss << '"';

    for (const char ch : s) {
        if (ch == '"') {
            oss << "\"\"";
        } else {
            oss << ch;
        }
    }

    oss << '"';
    return oss.str();
}

static bool get_leave_one_out_reference(const std::vector<TupleValue>& input_values,
                                        size_t leave_out_index,
                                        CentralValueMode mode,
                                        TupleValue& ref_out) {
    std::vector<TupleValue> loo_values;
    loo_values.reserve(input_values.size());

    for (size_t i = 0; i < input_values.size(); ++i) {
        if (i == leave_out_index) {
            continue;
        }

        loo_values.push_back(input_values[i]);
    }

    return combine_values(loo_values, mode, ref_out);
}

static void append_ratio_points_for_row(const CsvTable& table,
                                        const std::vector<std::string>& row,
                                        const CombinationCase& c,
                                        CentralValueMode mode,
                                        const std::vector<TupleValue>& input_values,
                                        const TupleValue& ref,
                                        std::vector<RatioPoint>& ratio_points,
                                        std::vector<RatioPoint>& loo_ratio_points) {
    for (const auto& var : fit_variable_configs()) {
        const std::string col = avg_column(var.column_prefix, "10.6 GeV");
        const double x = get_numeric_or_nan(table, row, col);

        if (!std::isfinite(x)) {
            continue;
        }

        for (size_t i = 0; i < c.inputs.size(); ++i) {
            const TupleValue& v = input_values[i];

            if (!v.ok ||
                !std::isfinite(v.value) ||
                !std::isfinite(v.stat) ||
                v.stat <= 0.0 ||
                !ref.ok ||
                !std::isfinite(ref.value) ||
                std::abs(ref.value) <= 0.0) {
                continue;
            }

            const double ratio = v.value / ref.value;
            const double ratio_stat = std::abs(v.stat / ref.value);

            if (std::isfinite(ratio) &&
                std::isfinite(ratio_stat) &&
                ratio_stat > 0.0) {
                RatioPoint p;
                p.case_label = c.label;
                p.central_mode = central_value_mode_name(mode);
                p.reference_type = "all_mean";
                p.period = c.inputs[i].period;
                p.variable_key = var.key;
                p.x = x;
                p.y = ratio;
                p.ey = ratio_stat;
                ratio_points.push_back(p);
            }

            TupleValue loo_ref;
            if (get_leave_one_out_reference(input_values, i, mode, loo_ref) &&
                loo_ref.ok &&
                std::isfinite(loo_ref.value) &&
                std::isfinite(loo_ref.stat) &&
                std::abs(loo_ref.value) > 0.0 &&
                std::abs(v.value) > 0.0) {
                const double loo_ratio = v.value / loo_ref.value;
                const double loo_ratio_stat = std::abs(loo_ratio) * std::sqrt(
                    (v.stat / v.value) * (v.stat / v.value) +
                    (loo_ref.stat / loo_ref.value) * (loo_ref.stat / loo_ref.value)
                );

                if (std::isfinite(loo_ratio) &&
                    std::isfinite(loo_ratio_stat) &&
                    loo_ratio_stat > 0.0) {
                    RatioPoint p;
                    p.case_label = c.label;
                    p.central_mode = central_value_mode_name(mode);
                    p.reference_type = "leave_one_out";
                    p.period = c.inputs[i].period;
                    p.variable_key = var.key;
                    p.x = x;
                    p.y = loo_ratio;
                    p.ey = loo_ratio_stat;
                    loo_ratio_points.push_back(p);
                }
            }
        }
    }
}

static CombinationResult evaluate_case(const CsvTable& table,
                                       const CombinationCase& c,
                                       CentralValueMode mode,
                                       std::vector<RatioPoint>* ratio_points,
                                       std::vector<RatioPoint>* loo_ratio_points) {
    CombinationResult result;
    result.label = c.label;
    result.central_value_mode = central_value_mode_name(mode);
    result.output_column = c.output_column;
    result.fill_output = c.fill_output;

    std::map<std::string, PeriodRatioAccumulator> acc_by_period;
    std::map<std::string, PeriodRatioAccumulator> loo_acc_by_period;

    for (const auto& input : c.inputs) {
        acc_by_period[input.period] = PeriodRatioAccumulator();
        loo_acc_by_period[input.period] = PeriodRatioAccumulator();
    }

    for (const auto& row : table.rows) {
        std::vector<TupleValue> input_values;
        input_values.reserve(c.inputs.size());

        for (const auto& input : c.inputs) {
            input_values.push_back(get_input_value(table, row, input, mode));
        }

        TupleValue ref;
        if (!combine_values(input_values, mode, ref)) {
            continue;
        }

        if (!ref.ok || !std::isfinite(ref.value) || std::abs(ref.value) <= 0.0) {
            continue;
        }

        if (ratio_points != nullptr && loo_ratio_points != nullptr) {
            append_ratio_points_for_row(table,
                                        row,
                                        c,
                                        mode,
                                        input_values,
                                        ref,
                                        *ratio_points,
                                        *loo_ratio_points);
        }

        int n_valid_in_bin = 0;

        for (size_t i = 0; i < c.inputs.size(); ++i) {
            const TupleValue& v = input_values[i];

            if (!v.ok || !std::isfinite(v.value) || v.stat <= 0.0) {
                continue;
            }

            const double ratio = v.value / ref.value;
            const double ratio_stat = std::abs(v.stat / ref.value);

            if (std::isfinite(ratio) &&
                std::isfinite(ratio_stat) &&
                ratio_stat > 0.0) {
                const double w = 1.0 / (ratio_stat * ratio_stat);

                PeriodRatioAccumulator& acc = acc_by_period[c.inputs[i].period];
                acc.sum_w += w;
                acc.sum_wr += w * ratio;
                acc.n += 1;

                ++result.ratio_points;
                ++n_valid_in_bin;
            }

            TupleValue loo_ref;
            if (get_leave_one_out_reference(input_values, i, mode, loo_ref) &&
                loo_ref.ok &&
                std::isfinite(loo_ref.value) &&
                std::isfinite(loo_ref.stat) &&
                std::abs(loo_ref.value) > 0.0 &&
                std::abs(v.value) > 0.0) {
                const double loo_ratio = v.value / loo_ref.value;

                const double loo_ratio_stat = std::abs(loo_ratio) * std::sqrt(
                    (v.stat / v.value) * (v.stat / v.value) +
                    (loo_ref.stat / loo_ref.value) * (loo_ref.stat / loo_ref.value)
                );

                if (std::isfinite(loo_ratio) &&
                    std::isfinite(loo_ratio_stat) &&
                    loo_ratio_stat > 0.0) {
                    const double loo_w = 1.0 / (loo_ratio_stat * loo_ratio_stat);

                    PeriodRatioAccumulator& loo_acc = loo_acc_by_period[c.inputs[i].period];
                    loo_acc.sum_w += loo_w;
                    loo_acc.sum_wr += loo_w * loo_ratio;
                    loo_acc.n += 1;

                    ++result.loo_ratio_points;
                }
            }
        }

        if (n_valid_in_bin >= 2) {
            ++result.valid_bins;
        }
    }

    double sum_obs2 = 0.0;
    double sum_stat2 = 0.0;
    int n_period = 0;

    double loo_sum_obs2 = 0.0;
    double loo_sum_stat2 = 0.0;
    int loo_n_period = 0;

    for (const auto& input : c.inputs) {
        PeriodRatioSummary p;
        p.period = input.period;

        const PeriodRatioAccumulator& acc = acc_by_period[input.period];
        if (acc.sum_w > 0.0 && acc.n > 0) {
            p.n = acc.n;
            p.mean_ratio = acc.sum_wr / acc.sum_w;
            p.mean_ratio_stat = 1.0 / std::sqrt(acc.sum_w);

            if (std::isfinite(p.mean_ratio) &&
                std::isfinite(p.mean_ratio_stat) &&
                p.mean_ratio_stat > 0.0) {
                const double residual = p.mean_ratio - 1.0;
                sum_obs2 += residual * residual;
                sum_stat2 += p.mean_ratio_stat * p.mean_ratio_stat;
                ++n_period;
            }
        }

        const PeriodRatioAccumulator& loo_acc = loo_acc_by_period[input.period];
        if (loo_acc.sum_w > 0.0 && loo_acc.n > 0) {
            p.loo_n = loo_acc.n;
            p.loo_mean_ratio = loo_acc.sum_wr / loo_acc.sum_w;
            p.loo_mean_ratio_stat = 1.0 / std::sqrt(loo_acc.sum_w);

            if (std::isfinite(p.loo_mean_ratio) &&
                std::isfinite(p.loo_mean_ratio_stat) &&
                p.loo_mean_ratio_stat > 0.0) {
                const double loo_residual = p.loo_mean_ratio - 1.0;
                loo_sum_obs2 += loo_residual * loo_residual;
                loo_sum_stat2 += p.loo_mean_ratio_stat * p.loo_mean_ratio_stat;
                ++loo_n_period;
            }
        }

        if (p.n > 0 || p.loo_n > 0) {
            result.period_summaries.push_back(p);
        }
    }

    if (n_period > 0) {
        result.s_obs = std::sqrt(std::max(0.0, sum_obs2 / (double)n_period));
        result.s_stat_exp = std::sqrt(std::max(0.0, sum_stat2 / (double)n_period));
        result.s_comb = std::sqrt(std::max(
            0.0,
            result.s_obs * result.s_obs -
            result.s_stat_exp * result.s_stat_exp
        ));
    }

    if (loo_n_period > 0) {
        result.loo_s_obs = std::sqrt(std::max(0.0, loo_sum_obs2 / (double)loo_n_period));
        result.loo_s_stat_exp = std::sqrt(std::max(0.0, loo_sum_stat2 / (double)loo_n_period));
        result.loo_s_comb = std::sqrt(std::max(
            0.0,
            result.loo_s_obs * result.loo_s_obs -
            result.loo_s_stat_exp * result.loo_s_stat_exp
        ));
    }

    return result;
}

static void fill_output_column(CsvTable& table,
                               const std::string& column,
                               double value) {
    const auto it = table.index.find(column);

    if (it == table.index.end()) {
        throw std::runtime_error("Missing output column: " + column);
    }

    const std::string value_string = format_double(value);
    const size_t icol = (size_t)it->second;

    for (auto& row : table.rows) {
        row[icol] = value_string;
    }
}

static void write_summary_csv(const std::string& path,
                              const std::vector<CombinationResult>& results) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open summary CSV: " + path);
    }

    fout << "case,central value mode,output column,fill output,valid bins,ratio points,loo ratio points,"
         << "s_obs_period,s_stat_period,s_comb,s_comb percent,"
         << "loo_s_obs_period,loo_s_stat_period,loo_s_comb,loo_s_comb percent,"
         << "period,period ratio points,period mean ratio,period mean ratio stat,"
         << "period loo ratio points,period loo mean ratio,period loo mean ratio stat\n";

    for (const auto& r : results) {
        if (r.period_summaries.empty()) {
            fout << csv_escape_field(r.label) << ","
                 << csv_escape_field(r.central_value_mode) << ","
                 << csv_escape_field(r.output_column) << ","
                 << (r.fill_output ? "true" : "false") << ","
                 << r.valid_bins << ","
                 << r.ratio_points << ","
                 << r.loo_ratio_points << ","
                 << format_double(r.s_obs) << ","
                 << format_double(r.s_stat_exp) << ","
                 << format_double(r.s_comb) << ","
                 << format_double(100.0 * r.s_comb) << ","
                 << format_double(r.loo_s_obs) << ","
                 << format_double(r.loo_s_stat_exp) << ","
                 << format_double(r.loo_s_comb) << ","
                 << format_double(100.0 * r.loo_s_comb) << ","
                 << ",,,,,,\n";
            continue;
        }

        for (const auto& p : r.period_summaries) {
            fout << csv_escape_field(r.label) << ","
                 << csv_escape_field(r.central_value_mode) << ","
                 << csv_escape_field(r.output_column) << ","
                 << (r.fill_output ? "true" : "false") << ","
                 << r.valid_bins << ","
                 << r.ratio_points << ","
                 << r.loo_ratio_points << ","
                 << format_double(r.s_obs) << ","
                 << format_double(r.s_stat_exp) << ","
                 << format_double(r.s_comb) << ","
                 << format_double(100.0 * r.s_comb) << ","
                 << format_double(r.loo_s_obs) << ","
                 << format_double(r.loo_s_stat_exp) << ","
                 << format_double(r.loo_s_comb) << ","
                 << format_double(100.0 * r.loo_s_comb) << ","
                 << csv_escape_field(p.period) << ","
                 << p.n << ","
                 << format_double(p.mean_ratio) << ","
                 << format_double(p.mean_ratio_stat) << ","
                 << p.loo_n << ","
                 << format_double(p.loo_mean_ratio) << ","
                 << format_double(p.loo_mean_ratio_stat) << "\n";
        }
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing summary CSV: " + path);
    }
}

static void print_summary_table(const std::vector<CombinationResult>& results) {
    std::cout << "\n[combination-systematics] Summary\n";

    std::cout << std::left
              << std::setw(30) << "case"
              << std::setw(18) << "central"
              << std::right
              << std::setw(10) << "bins"
              << std::setw(12) << "ratios"
              << std::setw(12) << "loo"
              << std::setw(14) << "s_obs"
              << std::setw(14) << "s_stat"
              << std::setw(14) << "s_comb"
              << std::setw(14) << "percent"
              << std::setw(14) << "loo_comb"
              << std::setw(14) << "loo_percent"
              << "\n";

    for (const auto& r : results) {
        std::cout << std::left
                  << std::setw(30) << r.label
                  << std::setw(18) << r.central_value_mode
                  << std::right
                  << std::setw(10) << r.valid_bins
                  << std::setw(12) << r.ratio_points
                  << std::setw(12) << r.loo_ratio_points
                  << std::setw(14) << std::setprecision(6) << r.s_obs
                  << std::setw(14) << std::setprecision(6) << r.s_stat_exp
                  << std::setw(14) << std::setprecision(6) << r.s_comb
                  << std::setw(14) << std::setprecision(6) << 100.0 * r.s_comb
                  << std::setw(14) << std::setprecision(6) << r.loo_s_comb
                  << std::setw(14) << std::setprecision(6) << 100.0 * r.loo_s_comb
                  << "\n";

        for (const auto& p : r.period_summaries) {
            std::cout << "    period mean ratio "
                      << std::left << std::setw(10) << p.period
                      << " = " << std::right
                      << std::setprecision(8) << p.mean_ratio
                      << " +/- " << std::setprecision(8) << p.mean_ratio_stat
                      << "   n=" << p.n
                      << "\n";

            std::cout << "    leave-one-out ratio "
                      << std::left << std::setw(10) << p.period
                      << " = " << std::right
                      << std::setprecision(8) << p.loo_mean_ratio
                      << " +/- " << std::setprecision(8) << p.loo_mean_ratio_stat
                      << "   n=" << p.loo_n
                      << "\n";
        }
    }

    std::cout << "\n";
}

static bool is_fill_mode(CentralValueMode mode) {
    return mode == kFillOutputMode;
}

static std::vector<double> choose_x_range(const std::vector<RatioPoint>& points) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();

    for (const auto& p : points) {
        if (!std::isfinite(p.x)) {
            continue;
        }

        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
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

static std::vector<double> choose_y_range(const std::vector<RatioPoint>& points) {
    double ymin = 1.0;
    double ymax = 1.0;

    for (const auto& p : points) {
        if (!std::isfinite(p.y)) {
            continue;
        }

        const double ey = (std::isfinite(p.ey) && p.ey > 0.0) ? p.ey : 0.0;

        ymin = std::min(ymin, p.y - ey);
        ymax = std::max(ymax, p.y + ey);
    }

    if (!std::isfinite(ymin) || !std::isfinite(ymax)) {
        return {0.0, 2.0};
    }

    ymin = std::max(0.0, ymin - 0.12 * (ymax - ymin));
    ymax = ymax + 0.12 * (ymax - ymin);

    if (ymax < 2.0) {
        ymax = 2.0;
    }

    if (ymin >= ymax) {
        ymin = 0.0;
        ymax = 2.0;
    }

    return {ymin, ymax};
}

static const VariableConfig* find_variable_config(const std::string& key) {
    static const std::vector<VariableConfig> vars = fit_variable_configs();

    for (const auto& var : vars) {
        if (var.key == key) {
            return &var;
        }
    }

    return nullptr;
}

static std::vector<RatioPoint> filter_ratio_points(const std::vector<RatioPoint>& points,
                                                   const std::string& case_label,
                                                   const std::string& central_mode,
                                                   const std::string& reference_type,
                                                   const std::string& period,
                                                   const std::string& variable_key) {
    std::vector<RatioPoint> out;

    for (const auto& p : points) {
        if (p.case_label != case_label) {
            continue;
        }

        if (p.central_mode != central_mode) {
            continue;
        }

        if (p.reference_type != reference_type) {
            continue;
        }

        if (p.period != period) {
            continue;
        }

        if (p.variable_key != variable_key) {
            continue;
        }

        if (!std::isfinite(p.x) ||
            !std::isfinite(p.y) ||
            !std::isfinite(p.ey) ||
            p.ey <= 0.0) {
            continue;
        }

        out.push_back(p);
    }

    std::sort(out.begin(), out.end(),
              [](const RatioPoint& a, const RatioPoint& b) {
                  return a.x < b.x;
              });

    return out;
}

static TGraphErrors* make_ratio_graph(const std::vector<RatioPoint>& points) {
    if (points.empty()) {
        return nullptr;
    }

    TGraphErrors* graph = new TGraphErrors((int)points.size());

    for (int i = 0; i < (int)points.size(); ++i) {
        graph->SetPoint(i, points[(size_t)i].x, points[(size_t)i].y);
        graph->SetPointError(i, 0.0, points[(size_t)i].ey);
    }

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.55);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlack);
    graph->SetLineWidth(1);

    return graph;
}

static void set_plot_style() {
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
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleOffset(1.10);
    frame->GetYaxis()->SetTitleOffset(1.25);

    return frame;
}

static std::string polynomial_formula(int order) {
    if (order < 0 || order > kMaxPolynomialOrder) {
        throw std::runtime_error("Unsupported polynomial order.");
    }

    std::ostringstream oss;
    oss << "pol" << order;
    return oss.str();
}

static FitResultSummary draw_and_fit_one_subset(const std::vector<RatioPoint>& points,
                                                const fs::path& out_dir,
                                                const std::string& case_label,
                                                const std::string& central_mode,
                                                const std::string& reference_type,
                                                const std::string& period,
                                                const VariableConfig& var,
                                                int order) {
    fs::create_directories(out_dir);

    FitResultSummary summary;
    summary.case_label = case_label;
    summary.central_mode = central_mode;
    summary.reference_type = reference_type;
    summary.period = period;
    summary.variable_key = var.key;
    summary.order = order;
    summary.n = (int)points.size();

    const std::vector<double> xr = choose_x_range(points);
    const std::vector<double> yr = choose_y_range(points);

    const std::string canvas_name =
        "c_fit_" + sanitize_for_path(case_label) + "_" +
        sanitize_for_path(central_mode) + "_" +
        sanitize_for_path(reference_type) + "_" +
        sanitize_for_path(period) + "_" +
        sanitize_for_path(var.key) + "_p" + std::to_string(order);

    TCanvas canvas(canvas_name.c_str(), canvas_name.c_str(), 1200, 900);
    canvas.SetFillColor(kWhite);
    canvas.cd();
    set_plot_style();

    std::unique_ptr<TH1D> frame(
        make_frame("frame_" + canvas_name,
                   var.title,
                   "#sigma_{i}/#bar{#sigma}",
                   xr[0],
                   xr[1],
                   yr[0],
                   yr[1])
    );

    frame->Draw("AXIS");

    TLine unity(xr[0], 1.0, xr[1], 1.0);
    unity.SetLineColor(kRed + 1);
    unity.SetLineStyle(2);
    unity.SetLineWidth(2);
    unity.Draw("SAME");

    std::unique_ptr<TGraphErrors> graph(make_ratio_graph(points));

    TF1* fit = nullptr;

    if (graph) {
        graph->Draw("P SAME");

        const std::string formula = polynomial_formula(order);
        fit = new TF1(("fit_" + canvas_name).c_str(), formula.c_str(), xr[0], xr[1]);
        fit->SetLineColor(kBlue + 1);
        fit->SetLineWidth(3);

        if ((int)points.size() >= order + 1) {
            TFitResultPtr fit_ptr = graph->Fit(fit, "QRS");

            fit->Draw("SAME");

            summary.chi2 = fit->GetChisquare();
            summary.ndf = fit->GetNDF();

            if (summary.ndf > 0) {
                summary.chi2_ndf = summary.chi2 / (double)summary.ndf;
            }

            for (int ip = 0; ip <= order; ++ip) {
                summary.params.push_back(fit->GetParameter(ip));
                summary.errors.push_back(fit->GetParError(ip));
            }
        } else {
            summary.ndf = (int)points.size() - (order + 1);
        }
    }

    TLegend legend(0.50, 0.68, 0.93, 0.88);
    legend.SetBorderSize(1);
    legend.SetFillStyle(1001);
    legend.SetFillColor(kWhite);
    legend.SetTextFont(42);
    legend.SetTextSize(0.032);

    if (graph) {
        legend.AddEntry(graph.get(), period.c_str(), "pe");
    }

    if (fit && std::isfinite(summary.chi2_ndf)) {
        std::ostringstream fit_label;
        fit_label << "p" << order << ", #chi^{2}/ndf = "
                  << std::fixed << std::setprecision(2)
                  << summary.chi2_ndf;
        legend.AddEntry(fit, fit_label.str().c_str(), "l");
    } else if (fit) {
        std::ostringstream fit_label;
        fit_label << "p" << order << ", fit unavailable";
        legend.AddEntry(fit, fit_label.str().c_str(), "l");
    }

    legend.Draw();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.034);
    title.SetTextAlign(22);

    std::ostringstream title_text;
    title_text << case_label << "   " << central_mode
               << "   " << reference_type
               << "   " << period
               << "   p" << order;
    title.DrawLatex(0.50, 0.965, title_text.str().c_str());

    canvas.Modified();
    canvas.Update();

    const fs::path out_path =
        out_dir /
        ("fit_" + sanitize_for_path(case_label) + "_" +
         sanitize_for_path(central_mode) + "_" +
         sanitize_for_path(reference_type) + "_" +
         sanitize_for_path(period) + "_" +
         sanitize_for_path(var.key) + "_p" +
         std::to_string(order) + ".png");

    canvas.SaveAs(out_path.string().c_str());

    delete fit;

    return summary;
}

static void write_fit_summary_csv(const fs::path& path,
                                  const std::vector<FitResultSummary>& summaries) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open fit summary CSV: " + path.string());
    }

    fout << "case,central mode,reference type,period,variable,polynomial order,n,ndf,chi2,chi2ndf";

    for (int ip = 0; ip <= kMaxPolynomialOrder; ++ip) {
        fout << ",p" << ip << ",p" << ip << " err";
    }

    fout << "\n";

    for (const auto& s : summaries) {
        fout << csv_escape_field(s.case_label) << ","
             << csv_escape_field(s.central_mode) << ","
             << csv_escape_field(s.reference_type) << ","
             << csv_escape_field(s.period) << ","
             << csv_escape_field(s.variable_key) << ","
             << s.order << ","
             << s.n << ","
             << s.ndf << ","
             << format_double(s.chi2) << ","
             << format_double(s.chi2_ndf);

        for (int ip = 0; ip <= kMaxPolynomialOrder; ++ip) {
            if (ip < (int)s.params.size()) {
                fout << "," << format_double(s.params[(size_t)ip])
                     << "," << format_double(s.errors[(size_t)ip]);
            } else {
                fout << ",,";
            }
        }

        fout << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing fit summary CSV: " + path.string());
    }
}

static void make_kinematic_fit_plots(const std::vector<RatioPoint>& all_ratio_points,
                                     const std::vector<RatioPoint>& all_loo_ratio_points,
                                     const std::vector<CombinationCase>& cases,
                                     const fs::path& out_dir) {
    const fs::path fit_root = out_dir / "kinematic_dependence_fits";
    fs::create_directories(fit_root);

    std::vector<FitResultSummary> summaries;

    const std::vector<VariableConfig> variables = fit_variable_configs();

    for (const auto& c : cases) {
        for (const auto& mode : central_value_modes()) {
            const std::string mode_name = central_value_mode_name(mode);

            for (const std::string& reference_type : {"all_mean", "leave_one_out"}) {
                const std::vector<RatioPoint>& source_points =
                    (reference_type == "all_mean") ? all_ratio_points : all_loo_ratio_points;

                for (const auto& input : c.inputs) {
                    for (const auto& var : variables) {
                        const std::vector<RatioPoint> points =
                            filter_ratio_points(source_points,
                                                c.label,
                                                mode_name,
                                                reference_type,
                                                input.period,
                                                var.key);

                        if (points.empty()) {
                            continue;
                        }

                        for (int order = 0; order <= kMaxPolynomialOrder; ++order) {
                            const fs::path subset_dir =
                                fit_root /
                                sanitize_for_path(c.label) /
                                sanitize_for_path(mode_name) /
                                sanitize_for_path(reference_type) /
                                sanitize_for_path(input.period) /
                                sanitize_for_path(var.key) /
                                ("p" + std::to_string(order));

                            FitResultSummary s =
                                draw_and_fit_one_subset(points,
                                                        subset_dir,
                                                        c.label,
                                                        mode_name,
                                                        reference_type,
                                                        input.period,
                                                        var,
                                                        order);

                            summaries.push_back(s);

                            std::cout << "[combination-systematics][fit] "
                                      << c.label << " "
                                      << mode_name << " "
                                      << reference_type << " "
                                      << input.period << " "
                                      << var.key << " p" << order
                                      << " n=" << s.n
                                      << " chi2/ndf=" << std::setprecision(8)
                                      << s.chi2_ndf
                                      << "\n";
                        }
                    }
                }
            }
        }
    }

    const fs::path summary_path = fit_root / "kinematic_dependence_fit_summary.csv";
    write_fit_summary_csv(summary_path, summaries);

    std::cout << "[combination-systematics] Wrote kinematic fit summary CSV: "
              << summary_path.string() << "\n";
}

} // namespace

bool combination_systematics(const std::string& csv_path,
                             const std::string& output_root_dir) {
    try {
        const std::string out_dir_string = output_root_dir + "/combination_systematics";
        const fs::path out_dir(out_dir_string);
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

        CsvTable table = read_csv_or_throw(csv_path);
        const std::vector<CombinationCase> cases = combination_cases();

        validate_schema(table, cases);

        std::cout << "[combination-systematics] CSV rows loaded: "
                  << table.rows.size() << "\n";
        std::cout << "[combination-systematics] Output directory: "
                  << out_dir_string << "\n";
        std::cout << "[combination-systematics] CSV output columns filled using central mode: "
                  << central_value_mode_name(kFillOutputMode) << "\n";

        std::vector<CombinationResult> results;
        results.reserve(cases.size() * central_value_modes().size() + 1);

        std::vector<RatioPoint> all_ratio_points;
        std::vector<RatioPoint> all_loo_ratio_points;

        double ten6_unpol_s_comb_for_fill =
            std::numeric_limits<double>::quiet_NaN();

        for (const auto& mode : central_value_modes()) {
            for (const auto& c : cases) {
                CombinationResult result =
                    evaluate_case(table,
                                  c,
                                  mode,
                                  &all_ratio_points,
                                  &all_loo_ratio_points);

                if (is_fill_mode(mode) && c.fill_output) {
                    fill_output_column(table, result.output_column, result.s_comb);
                }

                if (is_fill_mode(mode) && c.label == "10.6 GeV unpol") {
                    ten6_unpol_s_comb_for_fill = result.s_comb;
                }

                results.push_back(result);
            }
        }

        if (!std::isfinite(ten6_unpol_s_comb_for_fill)) {
            throw std::runtime_error("Could not determine 10.6 GeV unpol s_comb for Sp19 proxy fill.");
        }

        for (const auto& col : sp19_proxy_output_columns()) {
            fill_output_column(table, col, ten6_unpol_s_comb_for_fill);
        }

        CombinationResult sp19;
        sp19.label = "Sp19 Inb unpol";
        sp19.central_value_mode = central_value_mode_name(kFillOutputMode);
        sp19.output_column = combination_column("Sp19 Inb", "unpol");
        sp19.fill_output = true;
        sp19.s_obs = ten6_unpol_s_comb_for_fill;
        sp19.s_comb = ten6_unpol_s_comb_for_fill;
        sp19.loo_s_obs = ten6_unpol_s_comb_for_fill;
        sp19.loo_s_comb = ten6_unpol_s_comb_for_fill;
        results.push_back(sp19);

        const std::string summary_csv =
            out_dir_string + "/combination_systematics_summary.csv";

        write_summary_csv(summary_csv, results);
        print_summary_table(results);

        make_kinematic_fit_plots(all_ratio_points,
                                 all_loo_ratio_points,
                                 cases,
                                 out_dir);

        write_csv_or_throw(csv_path, table);

        std::cout << "[combination-systematics] Wrote summary CSV: "
                  << summary_csv << "\n";
        std::cout << "[combination-systematics] Updated CSV: "
                  << csv_path << "\n";

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[combination-systematics] FATAL: "
                  << e.what() << "\n";
        return false;
    }
}