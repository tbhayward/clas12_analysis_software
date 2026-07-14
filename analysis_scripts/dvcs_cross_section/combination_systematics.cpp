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

static constexpr int kMaxPolynomialOrder = 5;
static constexpr double kMinChi2NdfImprovement = 0.01;
static constexpr int kMaxFitWorkers = 5;
static constexpr double kExpectedSp19OverFa18InbRatio = 0.95;
static constexpr bool kVerboseFitOrderPrints = false;


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

enum class CombinationMetric {
    RelativeRatio,
    AbsoluteDifference
};

static bool is_absolute_metric(CombinationMetric metric) {
    return metric == CombinationMetric::AbsoluteDifference;
}

static std::string metric_name(CombinationMetric metric) {
    return is_absolute_metric(metric) ? "absolute_difference" : "relative_ratio";
}

static bool label_is_bsa(const std::string& label) {
    return label.rfind("BSA:", 0) == 0;
}

struct CombinationCase {
    std::string label;
    std::string output_column;
    bool fill_output = true;
    CombinationMetric metric = CombinationMetric::RelativeRatio;
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
};

struct CombinationResult {
    std::string label;
    std::string central_value_mode = "stat_weighted";
    std::string output_column;
    bool fill_output = false;
    CombinationMetric metric = CombinationMetric::RelativeRatio;

    int valid_bins = 0;
    int ratio_points = 0;

    double s_obs = 0.0;
    double s_stat_exp = 0.0;
    double s_comb = 0.0;

    double half_width = 0.0;

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

    // Auxiliary coordinate used for multi-dimensional diagnostic studies.
    // For the present workflow this stores the 10.6-GeV average photon polar angle
    // associated with the same row as the fitted ratio point. It is intentionally
    // not used by the ordinary one-dimensional fits.
    double g_theta = std::numeric_limits<double>::quiet_NaN();
};

struct DirectRatioPoint {
    std::string variable_key;
    double x = 0.0;
    double y = 0.0;
    double ey = 0.0;
};

struct DirectRatioSummary {
    std::string label;
    int n = 0;
    double weighted_mean = std::numeric_limits<double>::quiet_NaN();
    double weighted_mean_stat = std::numeric_limits<double>::quiet_NaN();
    double pull_from_unity = std::numeric_limits<double>::quiet_NaN();
    double pull_from_expected = std::numeric_limits<double>::quiet_NaN();
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

    double x_min = std::numeric_limits<double>::quiet_NaN();
    double x_max = std::numeric_limits<double>::quiet_NaN();

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
    int n_valid = 0;
    double mean_scale = std::numeric_limits<double>::quiet_NaN();
    double s_obs = std::numeric_limits<double>::quiet_NaN();
    double s_stat = 0.0;
    double s_comb = std::numeric_limits<double>::quiet_NaN();
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

static std::vector<std::string> sp19_proxy_output_columns();
static std::string cross_section_target_charge_column(const std::string& label,
                                                      const std::string& helicity);
static std::string cross_section_total_scale_column(const std::string& label,
                                                    const std::string& helicity);
static std::string bsa_beam_polarization_column(const std::string& label);
static std::string bsa_total_scale_column(const std::string& label);

static void ensure_column(CsvTable& table,
                          const std::string& column,
                          const std::string& initial_value = "") {
    if (column.empty()) {
        return;
    }

    if (table.index.find(column) != table.index.end()) {
        return;
    }

    table.index[column] = (int)table.header.size();
    table.header.push_back(column);

    for (auto& row : table.rows) {
        row.push_back(initial_value);
    }
}

static void ensure_output_columns(CsvTable& table,
                                  const std::vector<CombinationCase>& cases) {
    for (const auto& c : cases) {
        if (c.fill_output && !c.output_column.empty()) {
            ensure_column(table, c.output_column);
        }
    }

    for (const auto& col : sp19_proxy_output_columns()) {
        ensure_column(table, col);
    }

    const std::vector<std::pair<std::string, std::string> > xs_targets = {
        {"10.6 GeV", "unpol"}, {"Fa18", "unpol"}, {"Fa18", "pos"},
        {"Fa18", "neg"}, {"Sp18", "unpol"}, {"Sp19 Inb", "unpol"}
    };
    for (const auto& target : xs_targets) {
        ensure_column(table, cross_section_target_charge_column(target.first, target.second));
        ensure_column(table, cross_section_total_scale_column(target.first, target.second));
    }

    for (const auto& label : std::vector<std::string>{"10.6 GeV", "Fa18", "Sp18", "Sp19 Inb"}) {
        ensure_column(table, bsa_beam_polarization_column(label));
        ensure_column(table, bsa_total_scale_column(label));
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

    // Some upstream CSV writers produce ordinary tuple cells such as
    //     (0.1, 0.02, 0)
    // while the BSA stage can produce cells that arrive here with literal
    // quote characters still present, for example
    //     "(0.1,0.02,0)"
    // Strip balanced quote layers before stripping the tuple parentheses.
    // This keeps the parser backward compatible with the cross-section tuple
    // columns and makes it robust to the BSA tuple columns.
    bool changed = true;
    while (changed && !s.empty()) {
        changed = false;
        s = trim(s);

        if (s.size() >= 2U &&
            ((s.front() == '"' && s.back() == '"') ||
             (s.front() == '\'' && s.back() == '\''))) {
            s = trim(s.substr(1, s.size() - 2));
            changed = true;
        }
    }

    if (s.size() >= 2U && s.front() == '(' && s.back() == ')') {
        s = trim(s.substr(1, s.size() - 2));
    }

    changed = true;
    while (changed && !s.empty()) {
        changed = false;
        s = trim(s);

        if (s.size() >= 2U &&
            ((s.front() == '"' && s.back() == '"') ||
             (s.front() == '\'' && s.back() == '\''))) {
            s = trim(s.substr(1, s.size() - 2));
            changed = true;
        }
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

static TupleValue get_input_value(const CsvTable& table,
                                  const std::vector<std::string>& row,
                                  const PeriodInput& input) {
    std::vector<TupleValue> values;
    values.reserve(input.columns.size());

    for (const auto& col : input.columns) {
        values.push_back(get_tuple(table, row, col));
    }

    TupleValue out;

    if (!combine_stat_weighted(values, out)) {
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


static std::string cross_section_target_charge_column(const std::string& label,
                                                      const std::string& helicity) {
    return cross_section_column(label, helicity) + ", target thickness and charge sys";
}

static std::string cross_section_total_scale_column(const std::string& label,
                                                    const std::string& helicity) {
    return cross_section_column(label, helicity) + ", total scale sys";
}

static std::string bsa_column(const std::string& label) {
    return "BSA, counts, " + label;
}

static std::string bsa_combination_column(const std::string& label) {
    return bsa_column(label) + ", combination sys";
}

static std::string bsa_beam_polarization_column(const std::string& label) {
    return bsa_column(label) + ", beam polarization sys";
}

static std::string bsa_total_scale_column(const std::string& label) {
    return bsa_column(label) + ", total scale sys";
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

static CombinationCase make_relative_case(const std::string& label,
                                          const std::string& output_column,
                                          bool fill_output,
                                          std::vector<PeriodInput> inputs) {
    CombinationCase c;
    c.label = label;
    c.output_column = output_column;
    c.fill_output = fill_output;
    c.metric = CombinationMetric::RelativeRatio;
    c.inputs = std::move(inputs);
    return c;
}

static CombinationCase make_absolute_case(const std::string& label,
                                          const std::string& output_column,
                                          bool fill_output,
                                          std::vector<PeriodInput> inputs) {
    CombinationCase c;
    c.label = label;
    c.output_column = output_column;
    c.fill_output = fill_output;
    c.metric = CombinationMetric::AbsoluteDifference;
    c.inputs = std::move(inputs);
    return c;
}

static std::vector<CombinationCase> combination_cases() {
    std::vector<CombinationCase> cases;

    cases.push_back(make_relative_case(
        "cross section: 10.6 GeV unpol",
        combination_column("10.6 GeV", "unpol"),
        true,
        {
            input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
            input_single("Fa18 Out", cross_section_column("Fa18 Out", "unpol")),
            input_single("Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
            input_single("Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
        }
    ));

    cases.push_back(make_relative_case(
        "cross section: Fa18 unpol",
        combination_column("Fa18", "unpol"),
        true,
        {
            input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
            input_single("Fa18 Out", cross_section_column("Fa18 Out", "unpol"))
        }
    ));

    cases.push_back(make_relative_case(
        "cross section: Fa18 pos",
        combination_column("Fa18", "pos"),
        true,
        {
            input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "pos")),
            input_single("Fa18 Out", cross_section_column("Fa18 Out", "pos"))
        }
    ));

    cases.push_back(make_relative_case(
        "cross section: Fa18 neg",
        combination_column("Fa18", "neg"),
        true,
        {
            input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "neg")),
            input_single("Fa18 Out", cross_section_column("Fa18 Out", "neg"))
        }
    ));

    cases.push_back(make_relative_case(
        "cross section: Sp18 unpol",
        combination_column("Sp18", "unpol"),
        true,
        {
            input_single("Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")),
            input_single("Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
        }
    ));

    cases.push_back(make_relative_case(
        "cross section: Fa18 vs Sp18 unpol",
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
    ));

    cases.push_back(make_relative_case(
        "cross section: Inb vs Out unpol",
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
    ));

    cases.push_back(make_relative_case(
        "cross section: Fa18 Inb vs Sp18 Inb unpol",
        "",
        false,
        {
            input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
            input_single("Sp18 Inb", cross_section_column("Sp18 Inb", "unpol"))
        }
    ));

    cases.push_back(make_relative_case(
        "cross section: Fa18 Out vs Sp18 Out unpol",
        "",
        false,
        {
            input_single("Fa18 Out", cross_section_column("Fa18 Out", "unpol")),
            input_single("Sp18 Out", cross_section_column("Sp18 Out", "unpol"))
        }
    ));

    cases.push_back(make_relative_case(
        "cross section: Fa18 Inb vs Sp19 Inb unpol",
        "",
        false,
        {
            input_single("Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")),
            input_single("Sp19 Inb", cross_section_column("Sp19 Inb", "unpol"))
        }
    ));

    // BSA is dimensionless and crosses zero, so a ratio-to-mean definition is
    // mathematically ill-conditioned. For BSA cases we use the same structure
    // as the cross-section calculation, but with absolute residuals:
    //     delta_i = A_i - <A>
    //     s_comb = sqrt(RMS(delta_i)^2 - RMS(stat_i)^2)
    // The resulting output-column value is an absolute A_LU uncertainty.
    cases.push_back(make_absolute_case(
        "BSA: 10.6 GeV",
        bsa_combination_column("10.6 GeV"),
        true,
        {
            input_single("Fa18 Inb", bsa_column("Fa18 Inb")),
            input_single("Fa18 Out", bsa_column("Fa18 Out")),
            input_single("Sp18 Inb", bsa_column("Sp18 Inb")),
            input_single("Sp18 Out", bsa_column("Sp18 Out"))
        }
    ));

    cases.push_back(make_absolute_case(
        "BSA: Fa18",
        bsa_combination_column("Fa18"),
        true,
        {
            input_single("Fa18 Inb", bsa_column("Fa18 Inb")),
            input_single("Fa18 Out", bsa_column("Fa18 Out"))
        }
    ));

    cases.push_back(make_absolute_case(
        "BSA: Sp18",
        bsa_combination_column("Sp18"),
        true,
        {
            input_single("Sp18 Inb", bsa_column("Sp18 Inb")),
            input_single("Sp18 Out", bsa_column("Sp18 Out"))
        }
    ));

    cases.push_back(make_absolute_case(
        "BSA: Fa18 vs Sp18",
        "",
        false,
        {
            input_group("Fa18", {
                bsa_column("Fa18 Inb"),
                bsa_column("Fa18 Out")
            }),
            input_group("Sp18", {
                bsa_column("Sp18 Inb"),
                bsa_column("Sp18 Out")
            })
        }
    ));

    cases.push_back(make_absolute_case(
        "BSA: Inb vs Out",
        "",
        false,
        {
            input_group("Inb", {
                bsa_column("Fa18 Inb"),
                bsa_column("Sp18 Inb")
            }),
            input_group("Out", {
                bsa_column("Fa18 Out"),
                bsa_column("Sp18 Out")
            })
        }
    ));

    cases.push_back(make_absolute_case(
        "BSA: Fa18 Inb vs Sp18 Inb",
        "",
        false,
        {
            input_single("Fa18 Inb", bsa_column("Fa18 Inb")),
            input_single("Sp18 Inb", bsa_column("Sp18 Inb"))
        }
    ));

    cases.push_back(make_absolute_case(
        "BSA: Fa18 Out vs Sp18 Out",
        "",
        false,
        {
            input_single("Fa18 Out", bsa_column("Fa18 Out")),
            input_single("Sp18 Out", bsa_column("Sp18 Out"))
        }
    ));

    cases.push_back(make_absolute_case(
        "BSA: Fa18 Inb vs Sp19 Inb",
        "",
        false,
        {
            input_single("Fa18 Inb", bsa_column("Fa18 Inb")),
            input_single("Sp19 Inb", bsa_column("Sp19 Inb"))
        }
    ));

    return cases;
}

static std::vector<std::string> sp19_proxy_output_columns() {
    return {
        combination_column("Sp19 Inb", "unpol"),
        bsa_combination_column("Sp19 Inb")
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

    required.push_back(cross_section_column("Fa18 Inb", "unpol"));
    required.push_back(cross_section_column("Sp19 Inb", "unpol"));

    for (const auto& var : fit_variable_configs()) {
        required.push_back(avg_column(var.column_prefix, "Sp19 Inb"));
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

static bool has_complete_valid_fit_inputs(const CombinationCase& c,
                                          const std::vector<TupleValue>& input_values,
                                          const TupleValue& ref) {
    // The kinematic-dependence fits are intended to characterize the four-period
    // 10.6-GeV spread.  Do not let lower-statistics or missing-period bins enter
    // these fit/reference studies, because they can distort the apparent
    // variable dependence.  The global s_comb calculation itself is left
    // unchanged and can still use bins with two or more valid inputs.
    if (c.label != "cross section: 10.6 GeV unpol") {
        return false;
    }

    if (input_values.size() != c.inputs.size() || c.inputs.size() != 4) {
        return false;
    }

    if (!ref.ok || !std::isfinite(ref.value) || std::abs(ref.value) <= 0.0) {
        return false;
    }

    for (const auto& v : input_values) {
        if (!v.ok ||
            !std::isfinite(v.value) ||
            !std::isfinite(v.stat) ||
            v.stat <= 0.0 ||
            std::abs(v.value) <= 0.0) {
            return false;
        }
    }

    return true;
}

static void append_ratio_points_for_row(const CsvTable& table,
                                        const std::vector<std::string>& row,
                                        const CombinationCase& c,
                                        const std::vector<TupleValue>& input_values,
                                        const TupleValue& ref,
                                        std::vector<RatioPoint>& ratio_points) {
    if (!has_complete_valid_fit_inputs(c, input_values, ref)) {
        return;
    }

    const double g_theta_10p6 = get_numeric_or_nan(table, row, avg_column("g_theta", "10.6 GeV"));

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
                p.central_mode = "stat_weighted";
                p.reference_type = "all_mean";
                p.period = c.inputs[i].period;
                p.variable_key = var.key;
                p.x = x;
                p.y = ratio;
                p.ey = ratio_stat;
                p.g_theta = g_theta_10p6;
                ratio_points.push_back(p);
            }
        }
    }
}

static CombinationResult evaluate_case(const CsvTable& table,
                                       const CombinationCase& c,
                                       std::vector<RatioPoint>* ratio_points) {
    CombinationResult result;
    result.label = c.label;
    result.central_value_mode = "stat_weighted";
    result.output_column = c.output_column;
    result.fill_output = c.fill_output;
    result.metric = c.metric;

    std::map<std::string, PeriodRatioAccumulator> acc_by_period;

    for (const auto& input : c.inputs) {
        acc_by_period[input.period] = PeriodRatioAccumulator();
    }

    for (const auto& row : table.rows) {
        std::vector<TupleValue> input_values;
        input_values.reserve(c.inputs.size());

        for (const auto& input : c.inputs) {
            input_values.push_back(get_input_value(table, row, input));
        }

        TupleValue ref;
        if (!combine_stat_weighted(input_values, ref)) {
            continue;
        }

        if (!ref.ok || !std::isfinite(ref.value)) {
            continue;
        }

        if (!is_absolute_metric(c.metric) && std::abs(ref.value) <= 0.0) {
            continue;
        }

        if (ratio_points != nullptr) {
            append_ratio_points_for_row(table,
                                        row,
                                        c,
                                        input_values,
                                        ref,
                                        *ratio_points);
        }

        int n_valid_in_bin = 0;

        for (size_t i = 0; i < c.inputs.size(); ++i) {
            const TupleValue& v = input_values[i];

            if (!v.ok || !std::isfinite(v.value) || !std::isfinite(v.stat) || v.stat <= 0.0) {
                continue;
            }

            double metric_value = 0.0;
            double metric_stat = 0.0;

            if (is_absolute_metric(c.metric)) {
                metric_value = v.value - ref.value;
                metric_stat = std::abs(v.stat);
            } else {
                metric_value = v.value / ref.value;
                metric_stat = std::abs(v.stat / ref.value);
            }

            if (std::isfinite(metric_value) &&
                std::isfinite(metric_stat) &&
                metric_stat > 0.0) {
                const double w = 1.0 / (metric_stat * metric_stat);

                PeriodRatioAccumulator& acc = acc_by_period[c.inputs[i].period];
                acc.sum_w += w;
                acc.sum_wr += w * metric_value;
                acc.n += 1;

                ++result.ratio_points;
                ++n_valid_in_bin;
            }
        }

        if (n_valid_in_bin >= 2) {
            ++result.valid_bins;
        }
    }

    double sum_obs2 = 0.0;
    double sum_stat2 = 0.0;
    int n_period = 0;

    double min_metric = std::numeric_limits<double>::infinity();
    double max_metric = -std::numeric_limits<double>::infinity();

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
                const double residual = is_absolute_metric(c.metric)
                    ? p.mean_ratio
                    : (p.mean_ratio - 1.0);

                sum_obs2 += residual * residual;
                sum_stat2 += p.mean_ratio_stat * p.mean_ratio_stat;
                ++n_period;

                min_metric = std::min(min_metric, p.mean_ratio);
                max_metric = std::max(max_metric, p.mean_ratio);
            }

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

    if (std::isfinite(min_metric) &&
        std::isfinite(max_metric) &&
        max_metric >= min_metric) {
        result.half_width = 0.5 * (max_metric - min_metric);
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

static void fill_scale_systematic_columns(CsvTable& table) {
    static constexpr double kCrossSectionTargetChargeFraction = 0.0476;
    static constexpr double kBsaBeamPolarizationFraction = 0.0400;

    const std::vector<std::pair<std::string, std::string> > xs_targets = {
        {"10.6 GeV", "unpol"}, {"Fa18", "unpol"}, {"Fa18", "pos"},
        {"Fa18", "neg"}, {"Sp18", "unpol"}, {"Sp19 Inb", "unpol"}
    };

    for (const auto& target : xs_targets) {
        const int i_comb = table.index.at(combination_column(target.first, target.second));
        const int i_component = table.index.at(cross_section_target_charge_column(target.first, target.second));
        const int i_total = table.index.at(cross_section_total_scale_column(target.first, target.second));

        for (auto& row : table.rows) {
            double combination = 0.0;
            if (!parse_double(row[(size_t)i_comb], combination) || !std::isfinite(combination)) {
                row[(size_t)i_component].clear();
                row[(size_t)i_total].clear();
                continue;
            }
            row[(size_t)i_component] = format_double(kCrossSectionTargetChargeFraction);
            row[(size_t)i_total] = format_double(std::hypot(combination, kCrossSectionTargetChargeFraction));
        }
    }

    for (const auto& label : std::vector<std::string>{"10.6 GeV", "Fa18", "Sp18", "Sp19 Inb"}) {
        const int i_bsa = table.index.at(bsa_column(label));
        const int i_comb = table.index.at(bsa_combination_column(label));
        const int i_component = table.index.at(bsa_beam_polarization_column(label));
        const int i_total = table.index.at(bsa_total_scale_column(label));

        for (auto& row : table.rows) {
            const TupleValue bsa = parse_tuple_value(row[(size_t)i_bsa]);
            double combination = 0.0;
            if (!bsa.ok || !parse_double(row[(size_t)i_comb], combination) || !std::isfinite(combination)) {
                row[(size_t)i_component].clear();
                row[(size_t)i_total].clear();
                continue;
            }
            const double beam_abs = kBsaBeamPolarizationFraction * std::abs(bsa.value);
            row[(size_t)i_component] = format_double(beam_abs);
            row[(size_t)i_total] = format_double(std::hypot(combination, beam_abs));
        }
    }
}

static void write_summary_csv(const std::string& path,
                              const std::vector<CombinationResult>& results) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open summary CSV: " + path);
    }

    fout << "case,central value mode,metric,output column,fill output,valid bins,points,"
         << "s_obs,s_stat,s_comb,half_width,"
         << "s_obs_percent_for_relative_cases,s_stat_percent_for_relative_cases,"
         << "s_comb_percent_for_relative_cases,half_width_percent_for_relative_cases,"
         << "period,period points,period mean metric,period mean metric stat\n";

    for (const auto& r : results) {
        if (r.period_summaries.empty()) {
            fout << csv_escape_field(r.label) << ","
                 << csv_escape_field(r.central_value_mode) << ","
                 << csv_escape_field(metric_name(r.metric)) << ","
                 << csv_escape_field(r.output_column) << ","
                 << (r.fill_output ? "true" : "false") << ","
                 << r.valid_bins << ","
                 << r.ratio_points << ","
                 << format_double(r.s_obs) << ","
                 << format_double(r.s_stat_exp) << ","
                 << format_double(r.s_comb) << ","
                 << format_double(r.half_width) << ","
                 << (is_absolute_metric(r.metric) ? "" : format_double(100.0 * r.s_obs)) << ","
                 << (is_absolute_metric(r.metric) ? "" : format_double(100.0 * r.s_stat_exp)) << ","
                 << (is_absolute_metric(r.metric) ? "" : format_double(100.0 * r.s_comb)) << ","
                 << (is_absolute_metric(r.metric) ? "" : format_double(100.0 * r.half_width)) << ",,,,\n";
            continue;
        }

        for (const auto& p : r.period_summaries) {
            fout << csv_escape_field(r.label) << ","
                 << csv_escape_field(r.central_value_mode) << ","
                 << csv_escape_field(metric_name(r.metric)) << ","
                 << csv_escape_field(r.output_column) << ","
                 << (r.fill_output ? "true" : "false") << ","
                 << r.valid_bins << ","
                 << r.ratio_points << ","
                 << format_double(r.s_obs) << ","
                 << format_double(r.s_stat_exp) << ","
                 << format_double(r.s_comb) << ","
                 << format_double(r.half_width) << ","
                 << (is_absolute_metric(r.metric) ? "" : format_double(100.0 * r.s_obs)) << ","
                 << (is_absolute_metric(r.metric) ? "" : format_double(100.0 * r.s_stat_exp)) << ","
                 << (is_absolute_metric(r.metric) ? "" : format_double(100.0 * r.s_comb)) << ","
                 << (is_absolute_metric(r.metric) ? "" : format_double(100.0 * r.half_width)) << ","
                 << csv_escape_field(p.period) << ","
                 << p.n << ","
                 << format_double(p.mean_ratio) << ","
                 << format_double(p.mean_ratio_stat) << "\n";
        }
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing summary CSV: " + path);
    }
}

static void print_summary_table(const std::vector<CombinationResult>& results) {
    auto important = [](const CombinationResult& r) {
        return r.fill_output ||
               r.label == "cross section: Fa18 Inb vs Sp19 Inb unpol" ||
               r.label == "BSA: Fa18 Inb vs Sp19 Inb" ||
               r.label == "cross section: Fa18 vs Sp18 unpol" ||
               r.label == "BSA: Fa18 vs Sp18" ||
               r.label == "cross section: Inb vs Out unpol" ||
               r.label == "BSA: Inb vs Out";
    };

    std::cout << "\n[combination-systematics] Key relative cross-section s_comb results\n";
    std::cout << std::left
              << std::setw(46) << "case"
              << std::right
              << std::setw(9) << "bins"
              << std::setw(10) << "points"
              << std::setw(12) << "s_obs(%)"
              << std::setw(12) << "s_stat(%)"
              << std::setw(12) << "s_comb(%)"
              << std::setw(12) << "half(%)"
              << "\n";

    for (const auto& r : results) {
        if (!important(r) || is_absolute_metric(r.metric)) {
            continue;
        }

        std::cout << std::left
                  << std::setw(46) << r.label
                  << std::right
                  << std::setw(9) << r.valid_bins
                  << std::setw(10) << r.ratio_points
                  << std::setw(12) << std::fixed << std::setprecision(3) << 100.0 * r.s_obs
                  << std::setw(12) << std::fixed << std::setprecision(3) << 100.0 * r.s_stat_exp
                  << std::setw(12) << std::fixed << std::setprecision(3) << 100.0 * r.s_comb
                  << std::setw(12) << std::fixed << std::setprecision(3) << 100.0 * r.half_width
                  << "\n";
    }

    std::cout << "\n[combination-systematics] Key absolute BSA s_comb results\n";
    std::cout << std::left
              << std::setw(46) << "case"
              << std::right
              << std::setw(9) << "bins"
              << std::setw(10) << "points"
              << std::setw(12) << "s_obs"
              << std::setw(12) << "s_stat"
              << std::setw(12) << "s_comb"
              << std::setw(12) << "half"
              << "\n";

    for (const auto& r : results) {
        if (!important(r) || !is_absolute_metric(r.metric)) {
            continue;
        }

        std::cout << std::left
                  << std::setw(46) << r.label
                  << std::right
                  << std::setw(9) << r.valid_bins
                  << std::setw(10) << r.ratio_points
                  << std::setw(12) << std::fixed << std::setprecision(5) << r.s_obs
                  << std::setw(12) << std::fixed << std::setprecision(5) << r.s_stat_exp
                  << std::setw(12) << std::fixed << std::setprecision(5) << r.s_comb
                  << std::setw(12) << std::fixed << std::setprecision(5) << r.half_width
                  << "\n";
    }

    std::cout << "\n";
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
    double x_min = std::numeric_limits<double>::infinity();
    double x_max = -std::numeric_limits<double>::infinity();

    for (const auto& p : task.points) {
        if (!std::isfinite(p.x) ||
            !std::isfinite(p.y) ||
            !std::isfinite(p.ey) ||
            p.ey <= 0.0) {
            continue;
        }

        x_min = std::min(x_min, p.x);
        x_max = std::max(x_max, p.x);

        const double residual = p.y - eval_poly(out.params, p.x);
        chi2 += residual * residual / (p.ey * p.ey);
        ++n_used;
    }

    out.n = n_used;
    out.ndf = n_used - npar;
    out.chi2 = chi2;

    if (std::isfinite(x_min) && std::isfinite(x_max) && x_min <= x_max) {
        out.x_min = x_min;
        out.x_max = x_max;
    }

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

        if (kVerboseFitOrderPrints) {
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
    std::cout << "[combination-systematics] Fit/reference inputs require all four 10.6-GeV "
              << "period values to be valid in each bin.\n";

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
         << "chi2,chi2ndf,x_min,x_max,improvement from previous,attempted,accepted";

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
                 << format_double(s.x_min) << ","
                 << format_double(s.x_max) << ","
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

static DirectRatioSummary summarize_direct_ratio_points(const std::vector<DirectRatioPoint>& points,
                                                        const std::string& label) {
    DirectRatioSummary summary;
    summary.label = label;

    double sum_w = 0.0;
    double sum_wr = 0.0;

    for (const auto& p : points) {
        if (!std::isfinite(p.y) ||
            !std::isfinite(p.ey) ||
            p.ey <= 0.0) {
            continue;
        }

        const double w = 1.0 / (p.ey * p.ey);
        sum_w += w;
        sum_wr += w * p.y;
        summary.n += 1;
    }

    if (sum_w > 0.0 && summary.n > 0) {
        summary.weighted_mean = sum_wr / sum_w;
        summary.weighted_mean_stat = 1.0 / std::sqrt(sum_w);

        if (summary.weighted_mean_stat > 0.0) {
            summary.pull_from_unity =
                (summary.weighted_mean - 1.0) / summary.weighted_mean_stat;

            summary.pull_from_expected =
                (summary.weighted_mean - kExpectedSp19OverFa18InbRatio) /
                summary.weighted_mean_stat;
        }
    }

    return summary;
}

static std::map<std::string, std::vector<DirectRatioPoint> >
build_fa18_inb_sp19_inb_direct_ratio_points(const CsvTable& table) {
    std::map<std::string, std::vector<DirectRatioPoint> > out;

    const std::string fa18_col = cross_section_column("Fa18 Inb", "unpol");
    const std::string sp19_col = cross_section_column("Sp19 Inb", "unpol");

    for (const auto& row : table.rows) {
        const TupleValue fa18 = get_tuple(table, row, fa18_col);
        const TupleValue sp19 = get_tuple(table, row, sp19_col);

        if (!fa18.ok ||
            !sp19.ok ||
            !std::isfinite(fa18.value) ||
            !std::isfinite(sp19.value) ||
            !std::isfinite(fa18.stat) ||
            !std::isfinite(sp19.stat) ||
            fa18.stat <= 0.0 ||
            sp19.stat <= 0.0 ||
            std::abs(fa18.value) <= 0.0 ||
            std::abs(sp19.value) <= 0.0) {
            continue;
        }

        const double ratio = sp19.value / fa18.value;
        const double ratio_stat = std::abs(ratio) * std::sqrt(
            (sp19.stat / sp19.value) * (sp19.stat / sp19.value) +
            (fa18.stat / fa18.value) * (fa18.stat / fa18.value)
        );

        if (!std::isfinite(ratio) ||
            !std::isfinite(ratio_stat) ||
            ratio_stat <= 0.0) {
            continue;
        }

        for (const auto& var : fit_variable_configs()) {
            const std::string x_col = avg_column(var.column_prefix, "Sp19 Inb");
            const double x = get_numeric_or_nan(table, row, x_col);

            if (!std::isfinite(x)) {
                continue;
            }

            DirectRatioPoint p;
            p.variable_key = var.key;
            p.x = x;
            p.y = ratio;
            p.ey = ratio_stat;

            out[var.key].push_back(p);
        }
    }

    for (auto& kv : out) {
        std::sort(kv.second.begin(),
                  kv.second.end(),
                  [](const DirectRatioPoint& a, const DirectRatioPoint& b) {
                      return a.x < b.x;
                  });
    }

    return out;
}

static std::vector<double> choose_direct_ratio_x_range(const std::vector<DirectRatioPoint>& points) {
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

static TGraphErrors* make_direct_ratio_graph(const std::vector<DirectRatioPoint>& points) {
    if (points.empty()) {
        return nullptr;
    }

    TGraphErrors* graph = new TGraphErrors((int)points.size());

    for (int i = 0; i < (int)points.size(); ++i) {
        graph->SetPoint(i, points[(size_t)i].x, points[(size_t)i].y);
        graph->SetPointError(i, 0.0, points[(size_t)i].ey);
    }

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.36);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlack);
    graph->SetLineWidth(1);

    return graph;
}

static void write_fa18_inb_sp19_inb_direct_ratio_csv(
    const fs::path& path,
    const std::map<std::string, std::vector<DirectRatioPoint> >& points_by_var) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open Fa18/Sp19 direct-ratio CSV: " + path.string());
    }

    fout << "variable,x,sp19_over_fa18_inb,sp19_over_fa18_inb_stat\n";

    for (const auto& var : fit_variable_configs()) {
        const auto it = points_by_var.find(var.key);
        if (it == points_by_var.end()) {
            continue;
        }

        for (const auto& p : it->second) {
            fout << csv_escape_field(var.key) << ","
                 << format_double(p.x) << ","
                 << format_double(p.y) << ","
                 << format_double(p.ey) << "\n";
        }
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing Fa18/Sp19 direct-ratio CSV: " + path.string());
    }
}

static void write_fa18_inb_sp19_inb_direct_ratio_summary_csv(
    const fs::path& path,
    const std::map<std::string, std::vector<DirectRatioPoint> >& points_by_var) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open Fa18/Sp19 direct-ratio summary CSV: " + path.string());
    }

    fout << "variable,n,weighted mean,weighted mean stat,pull from 1,pull from expected 0.95,expected ratio\n";

    for (const auto& var : fit_variable_configs()) {
        const auto it = points_by_var.find(var.key);
        if (it == points_by_var.end()) {
            continue;
        }

        const DirectRatioSummary s =
            summarize_direct_ratio_points(it->second, var.key);

        fout << csv_escape_field(var.key) << ","
             << s.n << ","
             << format_double(s.weighted_mean) << ","
             << format_double(s.weighted_mean_stat) << ","
             << format_double(s.pull_from_unity) << ","
             << format_double(s.pull_from_expected) << ","
             << format_double(kExpectedSp19OverFa18InbRatio) << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing Fa18/Sp19 direct-ratio summary CSV: " + path.string());
    }
}

static void print_fa18_inb_sp19_inb_direct_ratio_summary(
    const std::map<std::string, std::vector<DirectRatioPoint> >& points_by_var) {
    std::cout << "\n[combination-systematics] Fa18 Inb vs Sp19 Inb direct diagnostic\n";
    std::cout << "[combination-systematics] Ratio shown is Sp19 Inb / Fa18 Inb. "
              << "BH-only expectation supplied by user: approximately "
              << kExpectedSp19OverFa18InbRatio << ".\n";

    std::cout << std::left
              << std::setw(12) << "variable"
              << std::right
              << std::setw(10) << "n"
              << std::setw(16) << "mean"
              << std::setw(16) << "stat"
              << std::setw(16) << "pull_1"
              << std::setw(16) << "pull_0.95"
              << "\n";

    for (const auto& var : fit_variable_configs()) {
        const auto it = points_by_var.find(var.key);
        if (it == points_by_var.end()) {
            continue;
        }

        const DirectRatioSummary s =
            summarize_direct_ratio_points(it->second, var.key);

        std::cout << std::left
                  << std::setw(12) << var.key
                  << std::right
                  << std::setw(10) << s.n
                  << std::setw(16) << std::setprecision(8) << s.weighted_mean
                  << std::setw(16) << std::setprecision(8) << s.weighted_mean_stat
                  << std::setw(16) << std::setprecision(8) << s.pull_from_unity
                  << std::setw(16) << std::setprecision(8) << s.pull_from_expected
                  << "\n";
    }

    std::cout << "\n";
}

static void draw_fa18_inb_sp19_inb_direct_ratio_canvas(
    const fs::path& out_dir,
    const std::map<std::string, std::vector<DirectRatioPoint> >& points_by_var) {
    const fs::path diag_dir = out_dir / "fa18_inb_vs_sp19_inb";
    fs::create_directories(diag_dir);

    const std::string canvas_name = "c_fa18_inb_vs_sp19_inb_direct_ratio";
    TCanvas canvas(canvas_name.c_str(), canvas_name.c_str(), 1900, 1300);
    canvas.SetFillColor(kWhite);
    canvas.Divide(4, 2, 0.001, 0.001);

    std::vector<std::unique_ptr<TH1D> > frames;
    std::vector<std::unique_ptr<TGraphErrors> > graphs;
    std::vector<std::unique_ptr<TLine> > lines;
    std::vector<std::unique_ptr<TLegend> > legends;
    std::vector<std::unique_ptr<TLatex> > labels;

    const std::vector<VariableConfig> vars = fit_variable_configs();

    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        canvas.cd(ivar + 1);
        set_plot_style();

        std::vector<DirectRatioPoint> points;
        const auto it = points_by_var.find(vars[(size_t)ivar].key);
        if (it != points_by_var.end()) {
            points = it->second;
        }

        const std::vector<double> xr = choose_direct_ratio_x_range(points);

        std::unique_ptr<TH1D> frame(
            make_frame("frame_fa18_sp19_" + vars[(size_t)ivar].key,
                       vars[(size_t)ivar].title,
                       "#sigma_{Sp19 Inb}/#sigma_{Fa18 Inb}",
                       xr[0],
                       xr[1],
                       0.0,
                       2.0)
        );

        frame->Draw("AXIS");

        std::unique_ptr<TLine> unity(new TLine(xr[0], 1.0, xr[1], 1.0));
        unity->SetLineColor(kRed + 1);
        unity->SetLineStyle(2);
        unity->SetLineWidth(2);
        unity->Draw("SAME");

        std::unique_ptr<TLine> expected(new TLine(xr[0],
                                                  kExpectedSp19OverFa18InbRatio,
                                                  xr[1],
                                                  kExpectedSp19OverFa18InbRatio));
        expected->SetLineColor(kBlue + 1);
        expected->SetLineStyle(7);
        expected->SetLineWidth(2);
        expected->Draw("SAME");

        std::unique_ptr<TGraphErrors> graph(make_direct_ratio_graph(points));

        if (graph) {
            graph->Draw("P SAME");
        }

        std::unique_ptr<TLegend> legend(new TLegend(0.48, 0.70, 0.94, 0.89));
        legend->SetBorderSize(1);
        legend->SetFillStyle(1001);
        legend->SetFillColor(kWhite);
        legend->SetTextFont(42);
        legend->SetTextSize(0.030);

        if (graph) {
            legend->AddEntry(graph.get(), "Sp19 Inb / Fa18 Inb", "pe");
        }

        legend->AddEntry(unity.get(), "unity", "l");
        legend->AddEntry(expected.get(), "BH-only expected #approx 0.95", "l");
        legend->Draw();

        std::unique_ptr<TLatex> label(new TLatex());
        label->SetNDC();
        label->SetTextFont(42);
        label->SetTextSize(0.055);
        label->SetTextAlign(22);
        label->DrawLatex(0.50, 0.955, vars[(size_t)ivar].title.c_str());

        frames.push_back(std::move(frame));
        lines.push_back(std::move(unity));
        lines.push_back(std::move(expected));
        graphs.push_back(std::move(graph));
        legends.push_back(std::move(legend));
        labels.push_back(std::move(label));
    }

    canvas.cd(8);
    gPad->SetFillColor(kWhite);

    TLatex info;
    info.SetNDC();
    info.SetTextFont(42);
    info.SetTextSize(0.055);
    info.SetTextAlign(22);
    info.DrawLatex(0.50, 0.78, "Direct energy/timing diagnostic");
    info.DrawLatex(0.50, 0.64, "Sp19 Inb: 10.200 GeV");
    info.DrawLatex(0.50, 0.52, "Fa18 Inb: 10.604 GeV");
    info.DrawLatex(0.50, 0.40, "Same torus orientation");
    info.DrawLatex(0.50, 0.28, "Same reconstruction pass");

    canvas.cd();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.022);
    title.SetTextAlign(22);
    title.DrawLatex(0.50,
                    0.992,
                    "Fa18 Inb vs Sp19 Inb unpolarized direct comparison");

    canvas.Modified();
    canvas.Update();

    const fs::path out_path = diag_dir / "fa18_inb_vs_sp19_inb_direct_ratio.png";
    canvas.SaveAs(out_path.string().c_str());
}

static void make_fa18_inb_sp19_inb_direct_diagnostic(const CsvTable& table,
                                                     const fs::path& out_dir) {
    const fs::path diag_dir = out_dir / "fa18_inb_vs_sp19_inb";
    fs::create_directories(diag_dir);

    const std::map<std::string, std::vector<DirectRatioPoint> > points_by_var =
        build_fa18_inb_sp19_inb_direct_ratio_points(table);

    const fs::path point_csv = diag_dir / "fa18_inb_vs_sp19_inb_direct_ratio_points.csv";
    const fs::path summary_csv = diag_dir / "fa18_inb_vs_sp19_inb_direct_ratio_summary.csv";

    write_fa18_inb_sp19_inb_direct_ratio_csv(point_csv, points_by_var);
    write_fa18_inb_sp19_inb_direct_ratio_summary_csv(summary_csv, points_by_var);
    print_fa18_inb_sp19_inb_direct_ratio_summary(points_by_var);
    draw_fa18_inb_sp19_inb_direct_ratio_canvas(out_dir, points_by_var);

    std::cout << "[combination-systematics] Wrote Fa18 Inb vs Sp19 Inb direct ratio CSV: "
              << point_csv.string() << "\n";
    std::cout << "[combination-systematics] Wrote Fa18 Inb vs Sp19 Inb direct ratio summary CSV: "
              << summary_csv.string() << "\n";
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
    title_text << "cross section: 10.6 GeV unpol, stat weighted, "
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

static std::vector<FitTask> build_fit_tasks(const std::vector<RatioPoint>& all_ratio_points) {
    std::vector<FitTask> tasks;

    const std::string case_label = "cross section: 10.6 GeV unpol";
    const std::string central_mode = "stat_weighted";
    const std::string reference_type = "all_mean";

    for (const auto& period : ten6_periods()) {
        for (const auto& variable : fit_variable_configs()) {
            std::vector<RatioPoint> points =
                filter_ratio_points(all_ratio_points,
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

    return tasks;
}

static std::map<std::string, FitResultSummary>
extract_final_variable_all_mean_fits(const std::vector<FitTaskResult>& fit_results,
                                     const std::string& variable_key) {
    std::map<std::string, FitResultSummary> out;

    for (const auto& result : fit_results) {
        if (result.task.case_label != "cross section: 10.6 GeV unpol") {
            continue;
        }

        if (result.task.central_mode != "stat_weighted") {
            continue;
        }

        if (result.task.reference_type != "all_mean") {
            continue;
        }

        if (result.task.variable.key != variable_key) {
            continue;
        }

        if (result.accepted.empty()) {
            continue;
        }

        out[result.task.period] = result.accepted.back();
    }

    for (const auto& period : ten6_periods()) {
        if (out.find(period) == out.end()) {
            throw std::runtime_error("Missing final accepted all_mean " + variable_key +
                                     " fit for period " + period);
        }
    }

    return out;
}

static bool value_inside_fit_support(const FitResultSummary& fit,
                                     double x) {
    if (!std::isfinite(fit.x_min) || !std::isfinite(fit.x_max)) {
        return false;
    }

    return x >= fit.x_min && x <= fit.x_max;
}

static std::vector<double>
build_reference_grid(const std::map<std::string, FitResultSummary>& fits_by_period,
                     const VariableConfig& variable) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();

    for (const auto& kv : fits_by_period) {
        const FitResultSummary& fit = kv.second;
        if (!std::isfinite(fit.x_min) || !std::isfinite(fit.x_max)) {
            continue;
        }
        xmin = std::min(xmin, fit.x_min);
        xmax = std::max(xmax, fit.x_max);
    }

    std::vector<double> grid;
    if (!std::isfinite(xmin) || !std::isfinite(xmax) || xmax <= xmin) {
        return grid;
    }

    if (variable.column_prefix == "e_theta" ||
        variable.column_prefix == "p_theta" ||
        variable.column_prefix == "g_theta") {
        const int lo = (int)std::ceil(xmin);
        const int hi = (int)std::floor(xmax);
        for (int i = lo; i <= hi; ++i) {
            grid.push_back((double)i);
        }
        return grid;
    }

    if (variable.key == "phi") {
        const double step = 5.0;
        const int n = (int)std::floor((xmax - xmin) / step);
        for (int i = 0; i <= n; ++i) {
            grid.push_back(xmin + step * (double)i);
        }
        if (grid.empty() || grid.back() < xmax) {
            grid.push_back(xmax);
        }
        return grid;
    }

    const int n_grid = 80;
    for (int i = 0; i < n_grid; ++i) {
        const double f = (n_grid == 1) ? 0.0 : (double)i / (double)(n_grid - 1);
        grid.push_back(xmin + f * (xmax - xmin));
    }

    return grid;
}

static std::vector<ScaleReferencePoint>
compute_scale_systematic_reference(const std::map<std::string, FitResultSummary>& fits_by_period,
                                   const VariableConfig& variable) {
    std::vector<ScaleReferencePoint> out;

    const std::vector<double> grid = build_reference_grid(fits_by_period, variable);

    for (const double x : grid) {
        std::vector<double> scales;

        for (const auto& period : ten6_periods()) {
            const auto it = fits_by_period.find(period);
            if (it == fits_by_period.end()) {
                throw std::runtime_error("Missing " + variable.key + " fit for period " + period);
            }

            const FitResultSummary& fit = it->second;
            if (!value_inside_fit_support(fit, x)) {
                continue;
            }

            const double scale = eval_poly(fit.params, x);
            if (std::isfinite(scale) && scale > 0.0) {
                scales.push_back(scale);
            }
        }

        ScaleReferencePoint p;
        p.theta = x;
        p.n_valid = (int)scales.size();

        if (scales.size() >= 2U) {
            double mean_scale = 0.0;
            for (const double scale : scales) {
                mean_scale += scale;
            }
            mean_scale /= (double)scales.size();

            if (std::isfinite(mean_scale) && std::abs(mean_scale) > 0.0) {
                double obs2 = 0.0;
                for (const double scale : scales) {
                    const double rel = scale / mean_scale - 1.0;
                    obs2 += rel * rel;
                }
                obs2 /= (double)scales.size();

                p.mean_scale = mean_scale;
                p.s_obs = std::sqrt(std::max(0.0, obs2));
                p.s_stat = 0.0;
                p.s_comb = p.s_obs;
            }
        }

        out.push_back(p);
    }

    return out;
}

static void write_scale_reference_csv(const fs::path& path,
                                      const VariableConfig& variable,
                                      const std::map<std::string, FitResultSummary>& fits_by_period,
                                      const std::vector<ScaleReferencePoint>& reference_points) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open scale reference CSV: " + path.string());
    }

    fout << csv_escape_field(variable.plain_title);

    for (const auto& period : ten6_periods()) {
        fout << "," << csv_escape_field(period) << " scale"
             << "," << csv_escape_field(period) << " in fit range";
    }

    fout << ",n_valid,mean_scale,s_obs,s_stat,s_comb,s_comb percent\n";

    for (const auto& ref : reference_points) {
        fout << format_double(ref.theta);

        for (const auto& period : ten6_periods()) {
            const FitResultSummary& fit = fits_by_period.at(period);
            const bool in_range = value_inside_fit_support(fit, ref.theta);
            const double scale = in_range ? eval_poly(fit.params, ref.theta)
                                          : std::numeric_limits<double>::quiet_NaN();

            fout << "," << format_double(scale)
                 << "," << (in_range ? "true" : "false");
        }

        fout << "," << ref.n_valid
             << "," << format_double(ref.mean_scale)
             << "," << format_double(ref.s_obs)
             << "," << format_double(ref.s_stat)
             << "," << format_double(ref.s_comb)
             << "," << format_double(100.0 * ref.s_comb)
             << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing scale reference CSV: " + path.string());
    }
}

static void print_scale_reference_summary(const VariableConfig& variable,
                                          const std::vector<ScaleReferencePoint>& reference_points) {
    int n_valid = 0;
    double sum_s_comb = 0.0;
    double max_s_comb = 0.0;

    for (const auto& ref : reference_points) {
        if (ref.n_valid < 2 || !std::isfinite(ref.s_comb)) {
            continue;
        }

        ++n_valid;
        sum_s_comb += ref.s_comb;
        max_s_comb = std::max(max_s_comb, ref.s_comb);
    }

    const double mean_s_comb = n_valid > 0 ? sum_s_comb / (double)n_valid : 0.0;

    std::cout << "[combination-systematics] " << variable.plain_title
              << " scale reference: "
              << n_valid << " usable points, mean s_comb="
              << std::fixed << std::setprecision(3) << 100.0 * mean_s_comb << "%"
              << ", max s_comb="
              << std::fixed << std::setprecision(3) << 100.0 * max_s_comb << "%\n";
}

static void draw_s_obs_canvas_for_variable(const fs::path& fit_root,
                                           const VariableConfig& variable,
                                           const std::vector<ScaleReferencePoint>& reference_points) {
    std::vector<double> x;
    std::vector<double> y;

    for (const auto& ref : reference_points) {
        if (ref.n_valid < 2) {
            continue;
        }

        if (!std::isfinite(ref.theta) ||
            !std::isfinite(ref.s_obs)) {
            continue;
        }

        x.push_back(ref.theta);
        y.push_back(ref.s_obs);
    }

    if (x.empty()) {
        std::cout << "[combination-systematics] No finite " << variable.plain_title
                  << " s_obs points available for plotting.\n";
        return;
    }

    double xmin = *std::min_element(x.begin(), x.end());
    double xmax = *std::max_element(x.begin(), x.end());
    if (!(xmax > xmin)) {
        xmin -= 1.0;
        xmax += 1.0;
    }

    double ymax = 0.0;
    for (const double v : y) {
        if (std::isfinite(v)) {
            ymax = std::max(ymax, v);
        }
    }

    if (!std::isfinite(ymax) || ymax <= 0.0) {
        ymax = 1.0;
    } else {
        ymax *= 1.20;
    }

    const std::string canvas_name = "c_" + sanitize_for_path(variable.key) + "_s_obs_dependence";
    TCanvas canvas(canvas_name.c_str(),
                   (variable.plain_title + " s_obs dependence").c_str(),
                   1100,
                   850);
    canvas.SetFillColor(kWhite);
    set_plot_style();

    std::unique_ptr<TH1D> frame(
        make_frame("frame_" + canvas_name,
                   variable.title,
                   "s_{obs}",
                   xmin,
                   xmax,
                   0.0,
                   ymax)
    );

    frame->Draw("AXIS");

    TGraphErrors graph((int)x.size());

    for (int i = 0; i < (int)x.size(); ++i) {
        graph.SetPoint(i, x[(size_t)i], y[(size_t)i]);
        graph.SetPointError(i, 0.0, 0.0);
    }

    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(0.9);
    graph.SetMarkerColor(kBlack);
    graph.SetLineColor(kBlack);
    graph.SetLineWidth(2);
    graph.Draw("PL SAME");

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.038);
    title.SetTextAlign(22);
    std::ostringstream title_text;
    title_text << variable.plain_title << "-dependent observed spread from run-period fits";
    title.DrawLatex(0.50, 0.955, title_text.str().c_str());

    TLatex note;
    note.SetNDC();
    note.SetTextFont(42);
    note.SetTextSize(0.026);
    note.SetTextAlign(13);
    note.DrawLatex(0.19, 0.86, "s_{obs} = RMS[f_{i}(x)/#LT f(x)#GT - 1]");
    note.DrawLatex(0.19, 0.82, "Only coordinates inside each period fit support are used.");

    canvas.Modified();
    canvas.Update();

    const fs::path out_path = fit_root / (sanitize_for_path(variable.key) + "_s_obs_dependence.png");
    canvas.SaveAs(out_path.string().c_str());
}


struct GammaThetaBinConfig {
    std::string key;
    std::string label;
    std::string title;
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
};

static std::vector<GammaThetaBinConfig> gamma_theta_bins_for_ptheta_study() {
    return {
        {"gtheta_lt_5",   "#theta_{#gamma} < 5#circ",        "gtheta < 5 deg",      -std::numeric_limits<double>::infinity(), 5.0},
        {"gtheta_5_10",   "5#circ #leq #theta_{#gamma} < 10#circ",  "5 <= gtheta < 10 deg",   5.0, 10.0},
        {"gtheta_10_15",  "10#circ #leq #theta_{#gamma} < 15#circ", "10 <= gtheta < 15 deg", 10.0, 15.0},
        {"gtheta_15_20",  "15#circ #leq #theta_{#gamma} < 20#circ", "15 <= gtheta < 20 deg", 15.0, 20.0},
        {"gtheta_ge_20",  "#theta_{#gamma} #geq 20#circ",       "gtheta >= 20 deg",    20.0, std::numeric_limits<double>::infinity()}
    };
}

static bool point_in_gamma_theta_bin(const RatioPoint& p,
                                     const GammaThetaBinConfig& bin) {
    if (!std::isfinite(p.g_theta)) {
        return false;
    }

    return p.g_theta >= bin.min && p.g_theta < bin.max;
}

static VariableConfig ptheta_variable_for_gamma_bin(const GammaThetaBinConfig& bin) {
    VariableConfig out;
    out.key = "p_theta_by_" + bin.key;
    out.column_prefix = "p_theta";
    out.title = "#theta_{p} (deg)";
    out.plain_title = "p_theta, " + bin.title;
    return out;
}

static std::vector<FitTask>
build_p_theta_by_gamma_theta_fit_tasks(const std::vector<RatioPoint>& all_ratio_points,
                                       const std::vector<GammaThetaBinConfig>& bins) {
    std::vector<FitTask> tasks;

    const std::string case_label = "cross section: 10.6 GeV unpol";
    const std::string central_mode = "stat_weighted";
    const std::string reference_type = "all_mean";

    for (const auto& bin : bins) {
        const VariableConfig variable = ptheta_variable_for_gamma_bin(bin);

        for (const auto& period : ten6_periods()) {
            std::vector<RatioPoint> points;

            for (const auto& p : all_ratio_points) {
                if (p.case_label != case_label ||
                    p.central_mode != central_mode ||
                    p.reference_type != reference_type ||
                    p.period != period ||
                    p.variable_key != "p_theta" ||
                    !point_in_gamma_theta_bin(p, bin)) {
                    continue;
                }

                points.push_back(p);
            }

            if (points.empty()) {
                continue;
            }

            std::sort(points.begin(),
                      points.end(),
                      [](const RatioPoint& a, const RatioPoint& b) {
                          return a.x < b.x;
                      });

            FitTask task;
            task.case_label = case_label + ", " + bin.title;
            task.central_mode = central_mode;
            task.reference_type = reference_type;
            task.period = period;
            task.variable = variable;
            task.points = std::move(points);

            tasks.push_back(std::move(task));
        }
    }

    return tasks;
}

static std::map<std::string, FitResultSummary>
extract_final_fits_for_variable_key(const std::vector<FitTaskResult>& fit_results,
                                    const std::string& variable_key) {
    std::map<std::string, FitResultSummary> out;

    for (const auto& result : fit_results) {
        if (result.task.variable.key != variable_key) {
            continue;
        }

        if (result.task.central_mode != "stat_weighted" ||
            result.task.reference_type != "all_mean" ||
            result.accepted.empty()) {
            continue;
        }

        out[result.task.period] = result.accepted.back();
    }

    for (const auto& period : ten6_periods()) {
        if (out.find(period) == out.end()) {
            throw std::runtime_error("Missing final accepted gamma-binned p_theta fit for period " + period +
                                     " in variable " + variable_key);
        }
    }

    return out;
}

static std::vector<double> reference_x_values(const std::vector<ScaleReferencePoint>& reference_points) {
    std::vector<double> x;
    for (const auto& ref : reference_points) {
        if (ref.n_valid >= 2 && std::isfinite(ref.theta) && std::isfinite(ref.s_obs)) {
            x.push_back(ref.theta);
        }
    }
    return x;
}

static std::vector<double> reference_y_values(const std::vector<ScaleReferencePoint>& reference_points) {
    std::vector<double> y;
    for (const auto& ref : reference_points) {
        if (ref.n_valid >= 2 && std::isfinite(ref.theta) && std::isfinite(ref.s_obs)) {
            y.push_back(ref.s_obs);
        }
    }
    return y;
}

static void draw_s_obs_canvas_for_gamma_binned_ptheta(const fs::path& fit_root,
                                                      const GammaThetaBinConfig& bin,
                                                      const VariableConfig& variable,
                                                      const std::vector<ScaleReferencePoint>& reference_points) {
    const std::vector<double> x = reference_x_values(reference_points);
    const std::vector<double> y = reference_y_values(reference_points);

    if (x.empty()) {
        std::cout << "[combination-systematics] No finite p_theta s_obs points for "
                  << bin.title << ".\n";
        return;
    }

    double xmin = *std::min_element(x.begin(), x.end());
    double xmax = *std::max_element(x.begin(), x.end());
    if (!(xmax > xmin)) {
        xmin -= 1.0;
        xmax += 1.0;
    }

    double ymax = 0.0;
    for (const double v : y) {
        if (std::isfinite(v)) {
            ymax = std::max(ymax, v);
        }
    }
    ymax = (std::isfinite(ymax) && ymax > 0.0) ? 1.20 * ymax : 1.0;

    const std::string canvas_name = "c_p_theta_s_obs_dependence_" + sanitize_for_path(bin.key);
    TCanvas canvas(canvas_name.c_str(),
                   ("p_theta s_obs dependence, " + bin.title).c_str(),
                   1100,
                   850);
    canvas.SetFillColor(kWhite);
    set_plot_style();

    std::unique_ptr<TH1D> frame(
        make_frame("frame_" + canvas_name,
                   "#theta_{p} (deg)",
                   "s_{obs}",
                   xmin,
                   xmax,
                   0.0,
                   ymax)
    );

    frame->Draw("AXIS");

    TGraphErrors graph((int)x.size());
    for (int i = 0; i < (int)x.size(); ++i) {
        graph.SetPoint(i, x[(size_t)i], y[(size_t)i]);
        graph.SetPointError(i, 0.0, 0.0);
    }

    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(0.9);
    graph.SetMarkerColor(kBlack);
    graph.SetLineColor(kBlack);
    graph.SetLineWidth(2);
    graph.Draw("PL SAME");

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.036);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.955, "p_{#theta}-dependent observed spread in #theta_{#gamma} bin");
    title.SetTextSize(0.031);
    title.DrawLatex(0.50, 0.915, bin.label.c_str());

    TLatex note;
    note.SetNDC();
    note.SetTextFont(42);
    note.SetTextSize(0.026);
    note.SetTextAlign(13);
    note.DrawLatex(0.19, 0.84, "s_{obs} = RMS[f_{i}(#theta_{p})/#LT f(#theta_{p})#GT - 1]");
    note.DrawLatex(0.19, 0.80, "Only coordinates inside each period fit support are used.");

    canvas.Modified();
    canvas.Update();

    const fs::path out_path = fit_root / ("p_theta_s_obs_dependence_" + sanitize_for_path(bin.key) + ".png");
    canvas.SaveAs(out_path.string().c_str());
}

static int gamma_bin_color(const int ibin) {
    static const int colors[] = {
        kBlue + 1,
        kRed + 1,
        kGreen + 2,
        kMagenta + 1,
        kOrange + 7,
        kCyan + 2,
        kViolet + 1
    };
    const int n = (int)(sizeof(colors) / sizeof(colors[0]));
    return colors[ibin % n];
}

static int gamma_bin_marker(const int ibin) {
    static const int markers[] = {20, 21, 22, 23, 33, 34, 29};
    const int n = (int)(sizeof(markers) / sizeof(markers[0]));
    return markers[ibin % n];
}

static double interpolate_reference_s_obs(const std::vector<ScaleReferencePoint>& reference_points,
                                          const double x) {
    std::vector<std::pair<double, double> > xy;
    xy.reserve(reference_points.size());

    for (const auto& ref : reference_points) {
        if (ref.n_valid >= 2 && std::isfinite(ref.theta) && std::isfinite(ref.s_obs)) {
            xy.emplace_back(ref.theta, ref.s_obs);
        }
    }

    if (xy.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::sort(xy.begin(), xy.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    constexpr double kExactTolerance = 1.0e-9;
    for (const auto& p : xy) {
        if (std::abs(p.first - x) < kExactTolerance) {
            return p.second;
        }
    }

    if (x < xy.front().first || x > xy.back().first) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    for (size_t i = 1; i < xy.size(); ++i) {
        const double x0 = xy[i - 1].first;
        const double y0 = xy[i - 1].second;
        const double x1 = xy[i].first;
        const double y1 = xy[i].second;

        if (x >= x0 && x <= x1 && x1 > x0) {
            const double f = (x - x0) / (x1 - x0);
            return y0 + f * (y1 - y0);
        }
    }

    return std::numeric_limits<double>::quiet_NaN();
}

static int count_ratio_points_near_ptheta_in_gamma_bin(const std::vector<RatioPoint>& all_ratio_points,
                                                       const GammaThetaBinConfig& bin,
                                                       const double theta_p) {
    int count = 0;

    for (const auto& p : all_ratio_points) {
        if (p.case_label != "cross section: 10.6 GeV unpol" ||
            p.central_mode != "stat_weighted" ||
            p.reference_type != "all_mean" ||
            p.variable_key != "p_theta" ||
            !point_in_gamma_theta_bin(p, bin) ||
            !std::isfinite(p.x)) {
            continue;
        }

        if (std::abs(p.x - theta_p) <= 0.5) {
            ++count;
        }
    }

    return count;
}

static std::vector<ScaleReferencePoint>
compute_weighted_recombined_gamma_sliced_ptheta_reference(
    const std::vector<RatioPoint>& all_ratio_points,
    const std::vector<GammaThetaBinConfig>& bins,
    const std::map<std::string, std::vector<ScaleReferencePoint> >& references_by_bin_key,
    const std::vector<ScaleReferencePoint>& inclusive_ptheta_reference_points) {

    std::vector<ScaleReferencePoint> out;

    for (const auto& inclusive_ref : inclusive_ptheta_reference_points) {
        if (inclusive_ref.n_valid < 2 ||
            !std::isfinite(inclusive_ref.theta) ||
            !std::isfinite(inclusive_ref.s_obs)) {
            continue;
        }

        const double theta_p = inclusive_ref.theta;
        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        int contributing_gamma_bins = 0;

        for (const auto& bin : bins) {
            const auto it = references_by_bin_key.find(bin.key);
            if (it == references_by_bin_key.end()) {
                continue;
            }

            const double sliced_s_obs = interpolate_reference_s_obs(it->second, theta_p);
            if (!std::isfinite(sliced_s_obs)) {
                continue;
            }

            int weight = count_ratio_points_near_ptheta_in_gamma_bin(all_ratio_points, bin, theta_p);
            if (weight <= 0) {
                // The fitted curve can still have support at this point even when no
                // actual row lands in the narrow +/-0.5 deg counting window.  Give
                // it unit weight rather than dropping it, but keep the weight small
                // so real populated regions dominate the recombination.
                weight = 1;
            }

            weighted_sum += (double)weight * sliced_s_obs;
            weight_sum += (double)weight;
            ++contributing_gamma_bins;
        }

        if (contributing_gamma_bins <= 0 || weight_sum <= 0.0) {
            continue;
        }

        ScaleReferencePoint p;
        p.theta = theta_p;
        p.n_valid = contributing_gamma_bins;
        p.mean_scale = std::numeric_limits<double>::quiet_NaN();
        p.s_obs = weighted_sum / weight_sum;
        p.s_stat = 0.0;
        p.s_comb = p.s_obs;
        out.push_back(p);
    }

    return out;
}

static void write_recombined_gamma_sliced_reference_csv(
    const fs::path& path,
    const std::vector<GammaThetaBinConfig>& bins,
    const std::map<std::string, std::vector<ScaleReferencePoint> >& references_by_bin_key,
    const std::vector<ScaleReferencePoint>& inclusive_ptheta_reference_points,
    const std::vector<ScaleReferencePoint>& recombined_reference_points) {

    std::ofstream fout(path);
    if (!fout.is_open()) {
        throw std::runtime_error("Could not open recombined gamma-sliced reference CSV: " + path.string());
    }

    fout << "theta_p,inclusive_s_obs,recombined_gamma_sliced_s_obs,recombined_n_gamma_bins";
    for (const auto& bin : bins) {
        fout << "," << csv_escape_field(bin.title + " s_obs");
    }
    fout << "\n";

    std::map<double, ScaleReferencePoint> recombined_by_theta;
    for (const auto& ref : recombined_reference_points) {
        if (std::isfinite(ref.theta)) {
            recombined_by_theta[ref.theta] = ref;
        }
    }

    for (const auto& inclusive_ref : inclusive_ptheta_reference_points) {
        if (inclusive_ref.n_valid < 2 ||
            !std::isfinite(inclusive_ref.theta) ||
            !std::isfinite(inclusive_ref.s_obs)) {
            continue;
        }

        const double theta_p = inclusive_ref.theta;
        const auto it_recombined = recombined_by_theta.find(theta_p);

        fout << format_double(theta_p)
             << "," << format_double(inclusive_ref.s_obs);

        if (it_recombined != recombined_by_theta.end()) {
            fout << "," << format_double(it_recombined->second.s_obs)
                 << "," << it_recombined->second.n_valid;
        } else {
            fout << ",,0";
        }

        for (const auto& bin : bins) {
            const auto it_bin = references_by_bin_key.find(bin.key);
            const double y = (it_bin == references_by_bin_key.end())
                                 ? std::numeric_limits<double>::quiet_NaN()
                                 : interpolate_reference_s_obs(it_bin->second, theta_p);
            fout << "," << format_double(y);
        }

        fout << "\n";
    }
}

static void draw_aggregate_gamma_binned_ptheta_canvas(
    const fs::path& fit_root,
    const std::vector<GammaThetaBinConfig>& bins,
    const std::map<std::string, std::vector<ScaleReferencePoint> >& references_by_bin_key,
    const std::vector<ScaleReferencePoint>& inclusive_ptheta_reference_points,
    const std::vector<ScaleReferencePoint>& recombined_gamma_sliced_reference_points) {

    double global_xmin = std::numeric_limits<double>::infinity();
    double global_xmax = -std::numeric_limits<double>::infinity();
    double global_ymax = 0.0;

    auto scan_xy = [&](const std::vector<ScaleReferencePoint>& refs) {
        const std::vector<double> x = reference_x_values(refs);
        const std::vector<double> y = reference_y_values(refs);
        for (const double v : x) {
            if (std::isfinite(v)) {
                global_xmin = std::min(global_xmin, v);
                global_xmax = std::max(global_xmax, v);
            }
        }
        for (const double v : y) {
            if (std::isfinite(v)) {
                global_ymax = std::max(global_ymax, v);
            }
        }
    };

    scan_xy(inclusive_ptheta_reference_points);
    scan_xy(recombined_gamma_sliced_reference_points);

    for (const auto& bin : bins) {
        const auto it = references_by_bin_key.find(bin.key);
        if (it != references_by_bin_key.end()) {
            scan_xy(it->second);
        }
    }

    if (!std::isfinite(global_xmin) || !std::isfinite(global_xmax) || global_xmax <= global_xmin) {
        std::cout << "[combination-systematics] No finite gamma-binned p_theta s_obs points for overlay canvas.\n";
        return;
    }

    const double xpad = 0.04 * (global_xmax - global_xmin);
    global_xmin -= xpad;
    global_xmax += xpad;
    global_ymax = (std::isfinite(global_ymax) && global_ymax > 0.0) ? 1.18 * global_ymax : 1.0;

    TCanvas canvas("c_p_theta_by_gtheta_s_obs_dependence",
                   "p_theta s_obs dependence by photon-angle bin",
                   1200,
                   900);
    canvas.SetFillColor(kWhite);
    set_plot_style();
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.13);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.04);

    std::unique_ptr<TH1D> frame(
        make_frame("frame_p_theta_by_gtheta_overlay",
                   "#theta_{p} (deg)",
                   "s_{obs}",
                   global_xmin,
                   global_xmax,
                   0.0,
                   global_ymax)
    );
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.038);
    frame->GetYaxis()->SetLabelSize(0.038);
    frame->Draw("AXIS");

    TLegend legend(0.50, 0.56, 0.93, 0.88);
    legend.SetBorderSize(1);
    legend.SetFillColor(kWhite);
    legend.SetTextFont(42);
    legend.SetTextSize(0.024);

    std::vector<std::unique_ptr<TGraphErrors> > graphs;
    graphs.reserve(bins.size() + 2U);

    auto make_graph = [](const std::vector<ScaleReferencePoint>& refs,
                         const int marker,
                         const int color,
                         const int line_width,
                         const int line_style) -> std::unique_ptr<TGraphErrors> {
        const std::vector<double> x = reference_x_values(refs);
        const std::vector<double> y = reference_y_values(refs);
        if (x.empty()) {
            return nullptr;
        }

        std::unique_ptr<TGraphErrors> graph(new TGraphErrors((int)x.size()));
        for (int i = 0; i < (int)x.size(); ++i) {
            graph->SetPoint(i, x[(size_t)i], y[(size_t)i]);
            graph->SetPointError(i, 0.0, 0.0);
        }
        graph->SetMarkerStyle(marker);
        graph->SetMarkerSize(0.85);
        graph->SetMarkerColor(color);
        graph->SetLineColor(color);
        graph->SetLineWidth(line_width);
        graph->SetLineStyle(line_style);
        return graph;
    };

    std::unique_ptr<TGraphErrors> inclusive_graph =
        make_graph(inclusive_ptheta_reference_points, 20, kBlack, 4, 1);
    if (inclusive_graph) {
        inclusive_graph->Draw("PL SAME");
        legend.AddEntry(inclusive_graph.get(), "inclusive #theta_{p} fit", "lp");
        graphs.push_back(std::move(inclusive_graph));
    }

    std::unique_ptr<TGraphErrors> recombined_graph =
        make_graph(recombined_gamma_sliced_reference_points, 24, kGray + 2, 4, 2);
    if (recombined_graph) {
        recombined_graph->Draw("PL SAME");
        legend.AddEntry(recombined_graph.get(), "weighted recombination of #theta_{#gamma} slices", "lp");
        graphs.push_back(std::move(recombined_graph));
    }

    for (int ib = 0; ib < (int)bins.size(); ++ib) {
        const auto& bin = bins[(size_t)ib];
        const auto it = references_by_bin_key.find(bin.key);
        if (it == references_by_bin_key.end()) {
            continue;
        }

        const int color = gamma_bin_color(ib);
        std::unique_ptr<TGraphErrors> graph =
            make_graph(it->second, gamma_bin_marker(ib), color, 2, 1);
        if (!graph) {
            continue;
        }

        graph->Draw("PL SAME");
        legend.AddEntry(graph.get(), bin.label.c_str(), "lp");
        graphs.push_back(std::move(graph));
    }

    legend.Draw();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.034);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.955, "#theta_{p}-dependent observed spread: inclusive vs #theta_{#gamma}-sliced fits");

    TLatex note;
    note.SetNDC();
    note.SetTextFont(42);
    note.SetTextSize(0.021);
    note.SetTextAlign(13);
    note.DrawLatex(0.16, 0.86, "black: original inclusive #theta_{p} fit");
    note.DrawLatex(0.16, 0.82, "gray dashed: population-weighted recombination of sliced curves");
    note.DrawLatex(0.16, 0.78, "colored: separate #theta_{p} fits inside #theta_{#gamma} bins");

    canvas.Modified();
    canvas.Update();

    const fs::path out_path = fit_root / "p_theta_by_gtheta_s_obs_dependence.png";
    canvas.SaveAs(out_path.string().c_str());
}

static void make_gamma_binned_ptheta_scale_reference(
    const std::vector<RatioPoint>& all_ratio_points,
    const fs::path& fit_root,
    const std::vector<ScaleReferencePoint>& inclusive_ptheta_reference_points) {
    const std::vector<GammaThetaBinConfig> bins = gamma_theta_bins_for_ptheta_study();
    const std::vector<FitTask> tasks = build_p_theta_by_gamma_theta_fit_tasks(all_ratio_points, bins);

    if (tasks.empty()) {
        std::cout << "[combination-systematics] No p_theta-by-g_theta fit tasks could be built.\n";
        return;
    }

    std::cout << "[combination-systematics] Gamma-binned p_theta fit tasks: "
              << tasks.size() << "\n";

    const std::vector<FitTaskResult> fit_results = run_fit_tasks_parallel(tasks);

    const fs::path summary_path = fit_root / "p_theta_by_gtheta_fit_summary.csv";
    write_fit_summary_csv(summary_path, fit_results);

    std::map<std::string, std::vector<ScaleReferencePoint> > references_by_bin_key;

    for (const auto& bin : bins) {
        const VariableConfig variable = ptheta_variable_for_gamma_bin(bin);

        try {
            const std::map<std::string, FitResultSummary> fits =
                extract_final_fits_for_variable_key(fit_results, variable.key);

            const std::vector<ScaleReferencePoint> reference_points =
                compute_scale_systematic_reference(fits, variable);

            references_by_bin_key[bin.key] = reference_points;

            const fs::path scale_path =
                fit_root / ("p_theta_scale_systematic_reference_" + sanitize_for_path(bin.key) + ".csv");

            write_scale_reference_csv(scale_path,
                                      variable,
                                      fits,
                                      reference_points);

            std::cout << "[combination-systematics] p_theta scale reference in "
                      << bin.title << ": ";
            print_scale_reference_summary(variable, reference_points);

            draw_s_obs_canvas_for_gamma_binned_ptheta(fit_root,
                                                      bin,
                                                      variable,
                                                      reference_points);
        } catch (const std::exception& e) {
            std::cout << "[combination-systematics] Skipping p_theta scale reference for "
                      << bin.title << ": " << e.what() << "\n";
        }
    }

    const std::vector<ScaleReferencePoint> recombined_gamma_sliced_reference_points =
        compute_weighted_recombined_gamma_sliced_ptheta_reference(
            all_ratio_points,
            bins,
            references_by_bin_key,
            inclusive_ptheta_reference_points);

    const fs::path recombined_path = fit_root / "p_theta_by_gtheta_recombined_reference.csv";
    write_recombined_gamma_sliced_reference_csv(recombined_path,
                                                bins,
                                                references_by_bin_key,
                                                inclusive_ptheta_reference_points,
                                                recombined_gamma_sliced_reference_points);

    draw_aggregate_gamma_binned_ptheta_canvas(fit_root,
                                              bins,
                                              references_by_bin_key,
                                              inclusive_ptheta_reference_points,
                                              recombined_gamma_sliced_reference_points);

    std::cout << "[combination-systematics] Wrote p_theta gamma-sliced recombination CSV: "
              << recombined_path.string() << "\n";

    std::cout << "[combination-systematics] Wrote gamma-binned p_theta fit summary CSV: "
              << summary_path.string() << "\n";
}

static std::vector<ScaleReferencePoint> make_kinematic_fit_plots(const std::vector<RatioPoint>& all_ratio_points,
                                                               const fs::path& out_dir) {
    const fs::path fit_root = out_dir / "kinematic_dependence_fits";
    fs::create_directories(fit_root);

    std::vector<FitTask> tasks = build_fit_tasks(all_ratio_points);

    std::cout << "[combination-systematics] Kinematic-dependence fit tasks: "
              << tasks.size() << "\n";

    std::vector<FitTaskResult> fit_results = run_fit_tasks_parallel(tasks);

    const fs::path summary_path = fit_root / "kinematic_dependence_fit_summary.csv";
    write_fit_summary_csv(summary_path, fit_results);

    for (const auto& variable : fit_variable_configs()) {
        draw_kinematic_canvas_for_variable(out_dir,
                                           "all_mean",
                                           variable,
                                           fit_results);
    }

    std::cout << "[combination-systematics] Scale-reference summaries from fitted period ratios:\n";

    std::vector<ScaleReferencePoint> inclusive_ptheta_reference_points;

    for (const auto& variable : fit_variable_configs()) {
        if (variable.key == "phi") {
            continue;
        }

        const std::map<std::string, FitResultSummary> fits =
            extract_final_variable_all_mean_fits(fit_results, variable.key);

        const std::vector<ScaleReferencePoint> reference_points =
            compute_scale_systematic_reference(fits, variable);

        const fs::path scale_path =
            fit_root / (sanitize_for_path(variable.key) + "_scale_systematic_reference.csv");

        write_scale_reference_csv(scale_path,
                                  variable,
                                  fits,
                                  reference_points);

        print_scale_reference_summary(variable,
                                      reference_points);

        draw_s_obs_canvas_for_variable(fit_root,
                                       variable,
                                       reference_points);

        if (variable.key == "p_theta") {
            inclusive_ptheta_reference_points = reference_points;
        }
    }

    make_gamma_binned_ptheta_scale_reference(all_ratio_points,
                                             fit_root,
                                             inclusive_ptheta_reference_points);

    std::cout << "[combination-systematics] Wrote kinematic fit summary CSV: "
              << summary_path.string() << "\n";
    std::cout << "[combination-systematics] Wrote variable-dependent scale-reference CSVs under: "
              << fit_root.string() << "\n";

    return inclusive_ptheta_reference_points;
}

static double interpolate_reference_s_obs_clamped(const std::vector<ScaleReferencePoint>& reference_points,
                                                  const double x,
                                                  bool& clamped) {
    clamped = false;

    std::vector<std::pair<double, double> > xy;
    xy.reserve(reference_points.size());

    for (const auto& ref : reference_points) {
        if (ref.n_valid >= 2 &&
            std::isfinite(ref.theta) &&
            std::isfinite(ref.s_obs)) {
            xy.emplace_back(ref.theta, ref.s_obs);
        }
    }

    if (xy.empty() || !std::isfinite(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::sort(xy.begin(), xy.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    constexpr double kExactTolerance = 1.0e-9;
    for (const auto& p : xy) {
        if (std::abs(p.first - x) < kExactTolerance) {
            return p.second;
        }
    }

    if (x < xy.front().first) {
        clamped = true;
        return xy.front().second;
    }

    if (x > xy.back().first) {
        clamped = true;
        return xy.back().second;
    }

    for (size_t i = 1; i < xy.size(); ++i) {
        const double x0 = xy[i - 1].first;
        const double y0 = xy[i - 1].second;
        const double x1 = xy[i].first;
        const double y1 = xy[i].second;

        if (x >= x0 && x <= x1 && x1 > x0) {
            const double f = (x - x0) / (x1 - x0);
            return y0 + f * (y1 - y0);
        }
    }

    return std::numeric_limits<double>::quiet_NaN();
}

static int fill_cross_section_combination_column_from_ptheta_reference(
    CsvTable& table,
    const std::string& output_column,
    const std::string& ptheta_average_label,
    const std::vector<ScaleReferencePoint>& ptheta_reference,
    const double fallback_value,
    int& n_interpolated,
    int& n_clamped,
    int& n_fallback) {

    n_interpolated = 0;
    n_clamped = 0;
    n_fallback = 0;

    ensure_column(table, output_column);

    const auto it_out = table.index.find(output_column);
    if (it_out == table.index.end()) {
        throw std::runtime_error("Missing output column after ensure_column: " + output_column);
    }

    const std::string ptheta_column = avg_column("p_theta", ptheta_average_label);
    const auto it_theta = table.index.find(ptheta_column);
    if (it_theta == table.index.end()) {
        throw std::runtime_error("Missing proton-angle column needed for theta_p-dependent systematic fill: " + ptheta_column);
    }

    for (auto& row : table.rows) {
        if (row.size() < table.header.size()) {
            row.resize(table.header.size());
        }

        double theta_p = std::numeric_limits<double>::quiet_NaN();
        parse_numeric_or_tuple_first(row[(size_t)it_theta->second], theta_p);

        bool clamped = false;
        double assigned = interpolate_reference_s_obs_clamped(ptheta_reference,
                                                              theta_p,
                                                              clamped);

        if (std::isfinite(assigned)) {
            ++n_interpolated;
            if (clamped) {
                ++n_clamped;
            }
        } else if (std::isfinite(fallback_value)) {
            assigned = fallback_value;
            ++n_fallback;
        }

        if (std::isfinite(assigned)) {
            row[(size_t)it_out->second] = format_double(assigned);
        }
    }

    return n_interpolated + n_fallback;
}

static void fill_final_cross_section_combination_systematics_from_ptheta_reference(
    CsvTable& table,
    const std::vector<ScaleReferencePoint>& ptheta_reference,
    const double fallback_ten6_s_comb) {

    if (ptheta_reference.empty()) {
        std::cout << "[combination-systematics] theta_p reference is empty; "
                  << "retaining scalar cross-section combination-systematic fills.\n";
        return;
    }

    int n_interp_10p6 = 0;
    int n_clamp_10p6 = 0;
    int n_fallback_10p6 = 0;
    const int n_fill_10p6 = fill_cross_section_combination_column_from_ptheta_reference(
        table,
        combination_column("10.6 GeV", "unpol"),
        "10.6 GeV",
        ptheta_reference,
        fallback_ten6_s_comb,
        n_interp_10p6,
        n_clamp_10p6,
        n_fallback_10p6);

    int n_interp_sp19 = 0;
    int n_clamp_sp19 = 0;
    int n_fallback_sp19 = 0;
    const int n_fill_sp19 = fill_cross_section_combination_column_from_ptheta_reference(
        table,
        combination_column("Sp19 Inb", "unpol"),
        "Sp19 Inb",
        ptheta_reference,
        fallback_ten6_s_comb,
        n_interp_sp19,
        n_clamp_sp19,
        n_fallback_sp19);

    std::cout << "[combination-systematics] Filled final cross-section combination-sys columns "
              << "from the inclusive theta_p scale-reference curve.\n";
    std::cout << "[combination-systematics]   10.6 GeV unpol: "
              << n_fill_10p6 << " rows filled, "
              << n_clamp_10p6 << " endpoint-clamped, "
              << n_fallback_10p6 << " scalar-fallback.\n";
    std::cout << "[combination-systematics]   Sp19 Inb unpol proxy: "
              << n_fill_sp19 << " rows filled, "
              << n_clamp_sp19 << " endpoint-clamped, "
              << n_fallback_sp19 << " scalar-fallback.\n";
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

        ensure_output_columns(table, cases);
        validate_schema(table, cases);

        std::cout << "[combination-systematics] CSV rows loaded: "
                  << table.rows.size() << "\n";
        std::cout << "[combination-systematics] Output directory: "
                  << out_dir_string << "\n";
        std::cout << "[combination-systematics] CSV output columns filled using central mode: stat_weighted\n";
        std::cout << "[combination-systematics] Cross sections use relative ratio spread; BSAs use absolute A_LU spread.\n";

        std::vector<CombinationResult> results;
        results.reserve(cases.size() + 2);

        std::vector<RatioPoint> all_ratio_points;

        double ten6_unpol_s_comb_for_fill =
            std::numeric_limits<double>::quiet_NaN();
        double ten6_bsa_s_comb_for_fill =
            std::numeric_limits<double>::quiet_NaN();

        for (const auto& c : cases) {
            CombinationResult result =
                evaluate_case(table,
                              c,
                              &all_ratio_points);

            if (c.fill_output) {
                fill_output_column(table, result.output_column, result.s_comb);
            }

            if (c.label == "cross section: 10.6 GeV unpol") {
                ten6_unpol_s_comb_for_fill = result.s_comb;
            } else if (c.label == "BSA: 10.6 GeV") {
                ten6_bsa_s_comb_for_fill = result.s_comb;
            }

            results.push_back(result);
        }

        if (!std::isfinite(ten6_unpol_s_comb_for_fill)) {
            throw std::runtime_error("Could not determine cross-section 10.6 GeV unpol s_comb for Sp19 proxy fill.");
        }

        fill_output_column(table,
                           combination_column("Sp19 Inb", "unpol"),
                           ten6_unpol_s_comb_for_fill);

        CombinationResult sp19_xs;
        sp19_xs.label = "cross section: Sp19 Inb unpol proxy";
        sp19_xs.central_value_mode = "stat_weighted";
        sp19_xs.output_column = combination_column("Sp19 Inb", "unpol");
        sp19_xs.fill_output = true;
        sp19_xs.metric = CombinationMetric::RelativeRatio;
        sp19_xs.s_obs = ten6_unpol_s_comb_for_fill;
        sp19_xs.s_comb = ten6_unpol_s_comb_for_fill;
        sp19_xs.half_width = ten6_unpol_s_comb_for_fill;
        results.push_back(sp19_xs);

        if (std::isfinite(ten6_bsa_s_comb_for_fill)) {
            fill_output_column(table,
                               bsa_combination_column("Sp19 Inb"),
                               ten6_bsa_s_comb_for_fill);

            CombinationResult sp19_bsa;
            sp19_bsa.label = "BSA: Sp19 Inb proxy";
            sp19_bsa.central_value_mode = "stat_weighted";
            sp19_bsa.output_column = bsa_combination_column("Sp19 Inb");
            sp19_bsa.fill_output = true;
            sp19_bsa.metric = CombinationMetric::AbsoluteDifference;
            sp19_bsa.s_obs = ten6_bsa_s_comb_for_fill;
            sp19_bsa.s_comb = ten6_bsa_s_comb_for_fill;
            sp19_bsa.half_width = ten6_bsa_s_comb_for_fill;
            results.push_back(sp19_bsa);
        } else {
            std::cout << "[combination-systematics] BSA 10.6 GeV s_comb was not finite; "
                      << "Sp19 BSA proxy combination-sys column was left blank.\n";
        }

        const std::string summary_csv =
            out_dir_string + "/combination_systematics_summary.csv";

        write_summary_csv(summary_csv, results);
        print_summary_table(results);

        const std::vector<ScaleReferencePoint> inclusive_ptheta_reference_points =
            make_kinematic_fit_plots(all_ratio_points,
                                     out_dir);

        fill_final_cross_section_combination_systematics_from_ptheta_reference(
            table,
            inclusive_ptheta_reference_points,
            ten6_unpol_s_comb_for_fill);

        make_fa18_inb_sp19_inb_direct_diagnostic(table,
                                                 out_dir);

        // Populate the independent normalization-scale components and the final
        // scale totals only after the theta_p-dependent combination columns have
        // received their final row-by-row values.
        fill_scale_systematic_columns(table);

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