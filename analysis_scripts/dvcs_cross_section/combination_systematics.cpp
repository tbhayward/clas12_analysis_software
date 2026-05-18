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

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
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
static constexpr double kMinChi2NdfImprovement = 0.01;
static constexpr int kMaxFitWorkers = 5;

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
    std::string plain_title;
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

    bool accepted = false;
    bool attempted = false;

    double chi2 = std::numeric_limits<double>::quiet_NaN();
    double chi2_ndf = std::numeric_limits<double>::quiet_NaN();
    double improvement_from_previous = std::numeric_limits<double>::quiet_NaN();

    std::vector<double> params;
    std::vector<double> errors;
};

struct FitTask {
    std::string case_label;
    std::string central_mode;
    std::string reference_type;
    std::string period;
    VariableConfig variable;
    std::vector<RatioPoint> points;
};

struct FitTaskResult {
    FitTask task;
    std::vector<FitResultSummary> attempted;
    std::vector<FitResultSummary> accepted;
};

struct ScaleReferencePoint {
    double theta = 0.0;
    double s_obs = 0.0;
    double s_stat = 0.0;
    double s_comb = 0.0;
};

static std::mutex gPrintMutex;

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
        {"xB",      "xBavg",     "x_{B}",                  "xB"},
        {"Q2",      "Q2avg",     "Q^{2} (GeV^{2})",        "Q2"},
        {"t",       "t_abs_avg", "-t (GeV^{2})",           "t"},
        {"phi",     "phiavg",    "#phi (deg)",             "phi"},
        {"e_theta", "e_theta",   "#theta_{e} (deg)",       "e_theta"},
        {"p_theta", "p_theta",   "#theta_{p} (deg)",       "p_theta"},
        {"g_theta", "g_theta",   "#theta_{#gamma} (deg)",  "g_theta"}
    };
}

static std::vector<std::string> ten6_periods() {
    return {
        "Fa18 Inb",
        "Fa18 Out",
        "Sp18 Inb",
        "Sp18 Out"
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
    }

    for (const auto& var : fit_variable_configs()) {
        required.push_back(avg_column(var.column_prefix, "10.6 GeV"));
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
    if (c.label != "10.6 GeV unpol") {
        return;
    }

    if (mode != CentralValueMode::StatWeighted) {
        return;
    }

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

static std::vector<double> choose_x_range_from_points(const std::vector<RatioPoint>& points) {
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

static bool solve_linear_system(std::vector<std::vector<double> > a,
                                std::vector<double> b,
                                std::vector<double>& x) {
    const int n = (int)b.size();

    if ((int)a.size() != n) {
        return false;
    }

    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double best = std::abs(a[(size_t)pivot][(size_t)col]);

        for (int row = col + 1; row < n; ++row) {
            const double v = std::abs(a[(size_t)row][(size_t)col]);
            if (v > best) {
                best = v;
                pivot = row;
            }
        }

        if (!std::isfinite(best) || best <= 1.0e-30) {
            return false;
        }

        if (pivot != col) {
            std::swap(a[(size_t)pivot], a[(size_t)col]);
            std::swap(b[(size_t)pivot], b[(size_t)col]);
        }

        const double diag = a[(size_t)col][(size_t)col];

        for (int j = col; j < n; ++j) {
            a[(size_t)col][(size_t)j] /= diag;
        }
        b[(size_t)col] /= diag;

        for (int row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }

            const double factor = a[(size_t)row][(size_t)col];
            if (factor == 0.0) {
                continue;
            }

            for (int j = col; j < n; ++j) {
                a[(size_t)row][(size_t)j] -= factor * a[(size_t)col][(size_t)j];
            }
            b[(size_t)row] -= factor * b[(size_t)col];
        }
    }

    x = b;

    for (const double v : x) {
        if (!std::isfinite(v)) {
            return false;
        }
    }

    return true;
}

static bool invert_matrix(const std::vector<std::vector<double> >& a,
                          std::vector<std::vector<double> >& inv) {
    const int n = (int)a.size();

    if (n <= 0) {
        return false;
    }

    inv.assign((size_t)n, std::vector<double>((size_t)n, 0.0));

    for (int col = 0; col < n; ++col) {
        std::vector<double> rhs((size_t)n, 0.0);
        rhs[(size_t)col] = 1.0;

        std::vector<double> sol;
        if (!solve_linear_system(a, rhs, sol)) {
            return false;
        }

        for (int row = 0; row < n; ++row) {
            inv[(size_t)row][(size_t)col] = sol[(size_t)row];
        }
    }

    return true;
}

static double eval_poly(const std::vector<double>& params,
                        double x) {
    double y = 0.0;
    double pow_x = 1.0;

    for (const double p : params) {
        y += p * pow_x;
        pow_x *= x;
    }

    return y;
}

static double eval_poly_error_diag_only(const std::vector<double>& errors,
                                        double x) {
    double var = 0.0;
    double pow_x = 1.0;

    for (const double e : errors) {
        if (std::isfinite(e)) {
            var += e * e * pow_x * pow_x;
        }

        pow_x *= x;
    }

    if (!std::isfinite(var) || var < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(var);
}

static FitResultSummary weighted_polynomial_fit(const FitTask& task,
                                                int order,
                                                double previous_chi2_ndf) {
    FitResultSummary out;
    out.case_label = task.case_label;
    out.central_mode = task.central_mode;
    out.reference_type = task.reference_type;
    out.period = task.period;
    out.variable_key = task.variable.key;
    out.order = order;
    out.n = (int)task.points.size();
    out.attempted = true;

    const int npar = order + 1;
    out.ndf = out.n - npar;

    if (out.n < npar || out.ndf <= 0) {
        return out;
    }

    std::vector<std::vector<double> > normal((size_t)npar, std::vector<double>((size_t)npar, 0.0));
    std::vector<double> rhs((size_t)npar, 0.0);

    for (const auto& p : task.points) {
        if (!std::isfinite(p.x) ||
            !std::isfinite(p.y) ||
            !std::isfinite(p.ey) ||
            p.ey <= 0.0) {
            continue;
        }

        const double w = 1.0 / (p.ey * p.ey);
        std::vector<double> xp((size_t)npar, 1.0);

        for (int i = 1; i < npar; ++i) {
            xp[(size_t)i] = xp[(size_t)(i - 1)] * p.x;
        }

        for (int i = 0; i < npar; ++i) {
            rhs[(size_t)i] += w * xp[(size_t)i] * p.y;

            for (int j = 0; j < npar; ++j) {
                normal[(size_t)i][(size_t)j] += w * xp[(size_t)i] * xp[(size_t)j];
            }
        }
    }

    std::vector<double> params;
    if (!solve_linear_system(normal, rhs, params)) {
        return out;
    }

    out.params = params;

    std::vector<std::vector<double> > cov;
    if (invert_matrix(normal, cov)) {
        out.errors.assign((size_t)npar, 0.0);

        for (int i = 0; i < npar; ++i) {
            const double v = cov[(size_t)i][(size_t)i];
            if (std::isfinite(v) && v >= 0.0) {
                out.errors[(size_t)i] = std::sqrt(v);
            }
        }
    } else {
        out.errors.assign((size_t)npar, std::numeric_limits<double>::quiet_NaN());
    }

    double chi2 = 0.0;
    int n_used = 0;

    for (const auto& p : task.points) {
        if (!std::isfinite(p.x) ||
            !std::isfinite(p.y) ||
            !std::isfinite(p.ey) ||
            p.ey <= 0.0) {
            continue;
        }

        const double residual = p.y - eval_poly(out.params, p.x);
        chi2 += residual * residual / (p.ey * p.ey);
        ++n_used;
    }

    out.n = n_used;
    out.ndf = n_used - npar;
    out.chi2 = chi2;

    if (out.ndf > 0) {
        out.chi2_ndf = out.chi2 / (double)out.ndf;
    }

    if (order == 0) {
        out.improvement_from_previous = std::numeric_limits<double>::quiet_NaN();
        out.accepted = std::isfinite(out.chi2_ndf);
    } else if (std::isfinite(previous_chi2_ndf) && std::isfinite(out.chi2_ndf)) {
        out.improvement_from_previous = previous_chi2_ndf - out.chi2_ndf;
        out.accepted = out.improvement_from_previous >= kMinChi2NdfImprovement;
    }

    return out;
}

static FitTaskResult run_iterative_fit_task(const FitTask& task) {
    FitTaskResult result;
    result.task = task;

    double previous_chi2_ndf = std::numeric_limits<double>::quiet_NaN();

    for (int order = 0; order <= kMaxPolynomialOrder; ++order) {
        FitResultSummary fit = weighted_polynomial_fit(task, order, previous_chi2_ndf);
        result.attempted.push_back(fit);

        {
            std::lock_guard<std::mutex> lock(gPrintMutex);

            std::cout << "[combination-systematics][fit] "
                      << fit.case_label << " "
                      << fit.central_mode << " "
                      << fit.reference_type << " "
                      << fit.period << " "
                      << fit.variable_key << " p" << fit.order
                      << " n=" << fit.n
                      << " chi2/ndf=" << std::setprecision(8)
                      << fit.chi2_ndf;

            if (order > 0) {
                std::cout << " improvement=" << std::setprecision(8)
                          << fit.improvement_from_previous;
            }

            std::cout << (fit.accepted ? " accepted" : " rejected")
                      << "\n";
        }

        if (!fit.accepted) {
            break;
        }

        previous_chi2_ndf = fit.chi2_ndf;
        result.accepted.push_back(fit);
    }

    return result;
}

static std::vector<FitTaskResult> run_fit_tasks_parallel(const std::vector<FitTask>& tasks) {
    std::vector<FitTaskResult> results;
    results.reserve(tasks.size());

    const unsigned int hardware_threads = std::max(1U, std::thread::hardware_concurrency());
    const int n_workers = std::max(1, std::min(kMaxFitWorkers, (int)hardware_threads));

    std::cout << "[combination-systematics] Running kinematic-dependence fits with "
              << n_workers << " worker(s).\n";

    size_t next_task = 0;

    while (next_task < tasks.size()) {
        std::vector<std::future<FitTaskResult> > futures;

        for (int i = 0; i < n_workers && next_task < tasks.size(); ++i) {
            const FitTask task = tasks[next_task];
            futures.push_back(std::async(std::launch::async,
                                         [task]() {
                                             return run_iterative_fit_task(task);
                                         }));
            ++next_task;
        }

        for (auto& fut : futures) {
            results.push_back(fut.get());
        }
    }

    return results;
}

static void write_fit_summary_csv(const fs::path& path,
                                  const std::vector<FitTaskResult>& task_results) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open fit summary CSV: " + path.string());
    }

    fout << "case,central mode,reference type,period,variable,polynomial order,n,ndf,"
         << "chi2,chi2ndf,improvement from previous,attempted,accepted";

    for (int ip = 0; ip <= kMaxPolynomialOrder; ++ip) {
        fout << ",p" << ip << ",p" << ip << " err";
    }

    fout << "\n";

    for (const auto& task_result : task_results) {
        for (const auto& s : task_result.attempted) {
            fout << csv_escape_field(s.case_label) << ","
                 << csv_escape_field(s.central_mode) << ","
                 << csv_escape_field(s.reference_type) << ","
                 << csv_escape_field(s.period) << ","
                 << csv_escape_field(s.variable_key) << ","
                 << s.order << ","
                 << s.n << ","
                 << s.ndf << ","
                 << format_double(s.chi2) << ","
                 << format_double(s.chi2_ndf) << ","
                 << format_double(s.improvement_from_previous) << ","
                 << (s.attempted ? "true" : "false") << ","
                 << (s.accepted ? "true" : "false");

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
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing fit summary CSV: " + path.string());
    }
}

static TGraphErrors* make_ratio_graph(const std::vector<RatioPoint>& points,
                                      int color,
                                      int marker_style) {
    if (points.empty()) {
        return nullptr;
    }

    TGraphErrors* graph = new TGraphErrors((int)points.size());

    for (int i = 0; i < (int)points.size(); ++i) {
        graph->SetPoint(i, points[(size_t)i].x, points[(size_t)i].y);
        graph->SetPointError(i, 0.0, points[(size_t)i].ey);
    }

    graph->SetMarkerStyle(marker_style);
    graph->SetMarkerSize(0.34);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetLineWidth(1);

    return graph;
}

static void set_plot_style() {
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.13);
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
    frame->GetXaxis()->SetTitleSize(0.060);
    frame->GetYaxis()->SetTitleSize(0.060);
    frame->GetXaxis()->SetLabelSize(0.050);
    frame->GetYaxis()->SetLabelSize(0.050);
    frame->GetXaxis()->SetTitleOffset(1.10);
    frame->GetYaxis()->SetTitleOffset(1.20);

    return frame;
}

static std::vector<double> get_combined_x_range_for_canvas(const std::vector<FitTaskResult>& period_results) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();

    for (const auto& r : period_results) {
        const std::vector<double> xr = choose_x_range_from_points(r.task.points);
        xmin = std::min(xmin, xr[0]);
        xmax = std::max(xmax, xr[1]);
    }

    if (!std::isfinite(xmin) || !std::isfinite(xmax) || xmin >= xmax) {
        return {0.0, 1.0};
    }

    return {xmin, xmax};
}

static std::string fit_function_formula(int order) {
    std::ostringstream oss;

    for (int i = 0; i <= order; ++i) {
        if (i > 0) {
            oss << " + ";
        }

        oss << "[" << i << "]";

        if (i == 1) {
            oss << "*x";
        } else if (i > 1) {
            oss << "*TMath::Power(x," << i << ")";
        }
    }

    return oss.str();
}

static TF1* make_root_fit_function(const FitResultSummary& summary,
                                   double xmin,
                                   double xmax,
                                   int color) {
    const std::string name =
        "f_" + sanitize_for_path(summary.reference_type) + "_" +
        sanitize_for_path(summary.period) + "_" +
        sanitize_for_path(summary.variable_key) + "_p" +
        std::to_string(summary.order);

    TF1* f = new TF1(name.c_str(), fit_function_formula(summary.order).c_str(), xmin, xmax);

    for (int i = 0; i < (int)summary.params.size(); ++i) {
        f->SetParameter(i, summary.params[(size_t)i]);
    }

    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(1 + summary.order);

    return f;
}

static void draw_kinematic_canvas_for_variable(const fs::path& out_dir,
                                               const std::string& reference_type,
                                               const VariableConfig& variable,
                                               const std::vector<FitTaskResult>& all_results) {
    std::vector<FitTaskResult> period_results;

    for (const auto& r : all_results) {
        if (r.task.reference_type != reference_type) {
            continue;
        }

        if (r.task.variable.key != variable.key) {
            continue;
        }

        period_results.push_back(r);
    }

    if (period_results.empty()) {
        return;
    }

    std::map<std::string, FitTaskResult> by_period;
    for (const auto& r : period_results) {
        by_period[r.task.period] = r;
    }

    const std::vector<double> xr = get_combined_x_range_for_canvas(period_results);

    const fs::path ref_dir = out_dir / "kinematic_dependence_fits" / sanitize_for_path(reference_type);
    fs::create_directories(ref_dir);

    const std::string canvas_name =
        "c_kinematic_dependence_" +
        sanitize_for_path(reference_type) + "_" +
        sanitize_for_path(variable.key);

    TCanvas canvas(canvas_name.c_str(), canvas_name.c_str(), 1700, 1300);
    canvas.SetFillColor(kWhite);
    canvas.Divide(2, 2, 0.001, 0.001);

    const std::vector<int> fit_colors = {
        kBlue + 1,
        kRed + 1,
        kGreen + 2,
        kMagenta + 1,
        kOrange + 7,
        kCyan + 2
    };

    const std::vector<std::string> periods = ten6_periods();

    std::vector<std::unique_ptr<TH1D> > frames;
    std::vector<std::unique_ptr<TGraphErrors> > graphs;
    std::vector<std::unique_ptr<TLine> > unity_lines;
    std::vector<std::unique_ptr<TLegend> > legends;
    std::vector<std::unique_ptr<TLatex> > labels;
    std::vector<std::unique_ptr<TF1> > functions;

    for (int iper = 0; iper < (int)periods.size(); ++iper) {
        canvas.cd(iper + 1);
        set_plot_style();

        const std::string& period = periods[(size_t)iper];

        std::unique_ptr<TH1D> frame(
            make_frame("frame_" + canvas_name + "_" + sanitize_for_path(period),
                       variable.title,
                       "#sigma_{i}/#bar{#sigma}",
                       xr[0],
                       xr[1],
                       0.0,
                       2.0)
        );
        frame->Draw("AXIS");

        std::unique_ptr<TLine> unity(new TLine(xr[0], 1.0, xr[1], 1.0));
        unity->SetLineColor(kRed + 1);
        unity->SetLineStyle(2);
        unity->SetLineWidth(1);
        unity->Draw("SAME");

        std::unique_ptr<TLegend> legend(new TLegend(0.56, 0.735, 0.94, 0.895));
        legend->SetBorderSize(1);
        legend->SetFillStyle(1001);
        legend->SetFillColor(kWhite);
        legend->SetTextFont(42);
        legend->SetTextSize(0.019);
        legend->SetTextAlign(32);
        legend->SetMargin(0.10);
        legend->SetEntrySeparation(0.05);

        const auto it_period = by_period.find(period);

        if (it_period != by_period.end()) {
            std::unique_ptr<TGraphErrors> graph(make_ratio_graph(it_period->second.task.points, kBlack, 20));

            if (graph) {
                graph->Draw("P SAME");
                legend->AddEntry(graph.get(), period.c_str(), "pe");
                graphs.push_back(std::move(graph));
            }

            for (int ifit = 0; ifit < (int)it_period->second.accepted.size(); ++ifit) {
                const FitResultSummary& fit = it_period->second.accepted[(size_t)ifit];
                const int color = fit_colors[(size_t)(ifit % (int)fit_colors.size())];

                std::unique_ptr<TF1> func(make_root_fit_function(fit, xr[0], xr[1], color));
                func->Draw("SAME");

                std::ostringstream label;
                label << "p" << fit.order
                      << ", #chi^{2}/ndf = "
                      << std::fixed << std::setprecision(2)
                      << fit.chi2_ndf;

                legend->AddEntry(func.get(), label.str().c_str(), "l");
                functions.push_back(std::move(func));
            }
        }

        legend->Draw();

        std::unique_ptr<TLatex> sublabel(new TLatex());
        sublabel->SetNDC();
        sublabel->SetTextFont(42);
        sublabel->SetTextSize(0.055);
        sublabel->SetTextAlign(22);
        sublabel->DrawLatex(0.50, 0.945, period.c_str());

        frames.push_back(std::move(frame));
        unity_lines.push_back(std::move(unity));
        legends.push_back(std::move(legend));
        labels.push_back(std::move(sublabel));
    }

    canvas.cd();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.018);
    title.SetTextAlign(22);

    std::ostringstream title_text;
    title_text << "10.6 GeV unpol, stat weighted, "
               << reference_type
               << ", kinematic dependence vs "
               << variable.plain_title;

    title.DrawLatex(0.50, 0.992, title_text.str().c_str());

    canvas.Modified();
    canvas.Update();

    const fs::path out_path =
        ref_dir /
        ("kinematic_dependence_" +
         sanitize_for_path(reference_type) + "_" +
         sanitize_for_path(variable.key) + ".png");

    canvas.SaveAs(out_path.string().c_str());
}

static std::vector<FitTask> build_fit_tasks(const std::vector<RatioPoint>& all_ratio_points,
                                            const std::vector<RatioPoint>& all_loo_ratio_points) {
    std::vector<FitTask> tasks;

    const std::string case_label = "10.6 GeV unpol";
    const std::string central_mode = central_value_mode_name(CentralValueMode::StatWeighted);
    const std::vector<std::string> reference_types = {
        "all_mean",
        "leave_one_out"
    };

    for (const auto& reference_type : reference_types) {
        const std::vector<RatioPoint>& source_points =
            (reference_type == "all_mean") ? all_ratio_points : all_loo_ratio_points;

        for (const auto& period : ten6_periods()) {
            for (const auto& variable : fit_variable_configs()) {
                std::vector<RatioPoint> points =
                    filter_ratio_points(source_points,
                                        case_label,
                                        central_mode,
                                        reference_type,
                                        period,
                                        variable.key);

                if (points.empty()) {
                    continue;
                }

                FitTask task;
                task.case_label = case_label;
                task.central_mode = central_mode;
                task.reference_type = reference_type;
                task.period = period;
                task.variable = variable;
                task.points = std::move(points);

                tasks.push_back(std::move(task));
            }
        }
    }

    return tasks;
}

static std::map<std::string, FitResultSummary>
extract_final_e_theta_all_mean_fits(const std::vector<FitTaskResult>& fit_results) {
    std::map<std::string, FitResultSummary> out;

    for (const auto& result : fit_results) {
        if (result.task.case_label != "10.6 GeV unpol") {
            continue;
        }

        if (result.task.central_mode != "stat_weighted") {
            continue;
        }

        if (result.task.reference_type != "all_mean") {
            continue;
        }

        if (result.task.variable.key != "e_theta") {
            continue;
        }

        if (result.accepted.empty()) {
            continue;
        }

        out[result.task.period] = result.accepted.back();
    }

    for (const auto& period : ten6_periods()) {
        if (out.find(period) == out.end()) {
            throw std::runtime_error("Missing final accepted all_mean e_theta fit for period " + period);
        }
    }

    return out;
}

static std::vector<ScaleReferencePoint>
compute_e_theta_scale_systematic_reference(const std::map<std::string, FitResultSummary>& fits_by_period) {
    std::vector<ScaleReferencePoint> out;

    for (int itheta = 1; itheta <= 30; ++itheta) {
        const double theta = (double)itheta;

        std::vector<double> residuals;
        std::vector<double> stat2_values;

        for (const auto& period : ten6_periods()) {
            const auto it = fits_by_period.find(period);
            if (it == fits_by_period.end()) {
                throw std::runtime_error("Missing e_theta fit for period " + period);
            }

            const FitResultSummary& fit = it->second;
            const double scale = eval_poly(fit.params, theta);
            const double scale_stat = eval_poly_error_diag_only(fit.errors, theta);

            if (std::isfinite(scale)) {
                residuals.push_back(scale - 1.0);
            }

            if (std::isfinite(scale_stat)) {
                stat2_values.push_back(scale_stat * scale_stat);
            }
        }

        if (residuals.empty()) {
            continue;
        }

        double obs2 = 0.0;
        for (const double r : residuals) {
            obs2 += r * r;
        }
        obs2 /= (double)residuals.size();

        double stat2 = 0.0;
        if (!stat2_values.empty()) {
            for (const double v : stat2_values) {
                stat2 += v;
            }
            stat2 /= (double)stat2_values.size();
        }

        ScaleReferencePoint p;
        p.theta = theta;
        p.s_obs = std::sqrt(std::max(0.0, obs2));
        p.s_stat = std::sqrt(std::max(0.0, stat2));
        p.s_comb = std::sqrt(std::max(0.0, obs2 - stat2));

        out.push_back(p);
    }

    return out;
}

static void write_e_theta_scale_reference_csv(const fs::path& path,
                                              const std::map<std::string, FitResultSummary>& fits_by_period,
                                              const std::vector<ScaleReferencePoint>& reference_points) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open e_theta scale reference CSV: " + path.string());
    }

    fout << "theta_e_deg";

    for (const auto& period : ten6_periods()) {
        fout << "," << csv_escape_field(period) << " scale"
             << "," << csv_escape_field(period) << " scale stat";
    }

    fout << ",s_obs,s_stat,s_comb,s_comb percent\n";

    for (const auto& ref : reference_points) {
        fout << format_double(ref.theta);

        for (const auto& period : ten6_periods()) {
            const FitResultSummary& fit = fits_by_period.at(period);
            const double scale = eval_poly(fit.params, ref.theta);
            const double scale_stat = eval_poly_error_diag_only(fit.errors, ref.theta);

            fout << "," << format_double(scale)
                 << "," << format_double(scale_stat);
        }

        fout << "," << format_double(ref.s_obs)
             << "," << format_double(ref.s_stat)
             << "," << format_double(ref.s_comb)
             << "," << format_double(100.0 * ref.s_comb)
             << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing e_theta scale reference CSV: " + path.string());
    }
}

static void print_e_theta_scale_reference(const std::map<std::string, FitResultSummary>& fits_by_period,
                                          const std::vector<ScaleReferencePoint>& reference_points) {
    std::cout << "\n[combination-systematics] e_theta-dependent scale systematic reference\n";
    std::cout << std::left
              << std::setw(10) << "theta"
              << std::right;

    for (const auto& period : ten6_periods()) {
        std::cout << std::setw(16) << (sanitize_for_path(period) + "_f");
    }

    std::cout << std::setw(14) << "s_obs"
              << std::setw(14) << "s_stat"
              << std::setw(14) << "s_comb"
              << std::setw(14) << "percent"
              << "\n";

    for (const auto& ref : reference_points) {
        std::cout << std::left
                  << std::setw(10) << std::setprecision(4) << ref.theta
                  << std::right;

        for (const auto& period : ten6_periods()) {
            const FitResultSummary& fit = fits_by_period.at(period);
            const double scale = eval_poly(fit.params, ref.theta);
            std::cout << std::setw(16) << std::setprecision(8) << scale;
        }

        std::cout << std::setw(14) << std::setprecision(8) << ref.s_obs
                  << std::setw(14) << std::setprecision(8) << ref.s_stat
                  << std::setw(14) << std::setprecision(8) << ref.s_comb
                  << std::setw(14) << std::setprecision(8) << 100.0 * ref.s_comb
                  << "\n";
    }

    std::cout << "\n";
}

static void make_kinematic_fit_plots(const std::vector<RatioPoint>& all_ratio_points,
                                     const std::vector<RatioPoint>& all_loo_ratio_points,
                                     const fs::path& out_dir) {
    const fs::path fit_root = out_dir / "kinematic_dependence_fits";
    fs::create_directories(fit_root);

    std::vector<FitTask> tasks = build_fit_tasks(all_ratio_points, all_loo_ratio_points);

    std::cout << "[combination-systematics] Kinematic-dependence fit tasks: "
              << tasks.size() << "\n";

    std::vector<FitTaskResult> fit_results = run_fit_tasks_parallel(tasks);

    const fs::path summary_path = fit_root / "kinematic_dependence_fit_summary.csv";
    write_fit_summary_csv(summary_path, fit_results);

    for (const std::string& reference_type : {"all_mean", "leave_one_out"}) {
        for (const auto& variable : fit_variable_configs()) {
            draw_kinematic_canvas_for_variable(out_dir,
                                               reference_type,
                                               variable,
                                               fit_results);
        }
    }

    const std::map<std::string, FitResultSummary> e_theta_fits =
        extract_final_e_theta_all_mean_fits(fit_results);

    const std::vector<ScaleReferencePoint> reference_points =
        compute_e_theta_scale_systematic_reference(e_theta_fits);

    const fs::path e_theta_scale_path =
        fit_root / "e_theta_scale_systematic_reference.csv";

    write_e_theta_scale_reference_csv(e_theta_scale_path,
                                      e_theta_fits,
                                      reference_points);

    print_e_theta_scale_reference(e_theta_fits,
                                  reference_points);

    std::cout << "[combination-systematics] Wrote kinematic fit summary CSV: "
              << summary_path.string() << "\n";
    std::cout << "[combination-systematics] Wrote e_theta scale systematic reference CSV: "
              << e_theta_scale_path.string() << "\n";
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
        gStyle->SetEndErrorSize(1);
        gStyle->SetErrorX(0.0);
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