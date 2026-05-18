#include "combination_systematics.h"

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

enum class CentralValueMode {
    StatWeighted,
    EqualPeriod,
    RandomEffects
};

static constexpr CentralValueMode kFillOutputMode = CentralValueMode::StatWeighted;

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
};

struct CombinationResult {
    std::string label;
    std::string central_value_mode;
    std::string output_column;
    bool fill_output = false;

    int valid_bins = 0;
    int ratio_points = 0;

    double s_obs = 0.0;
    double s_stat_exp = 0.0;
    double s_comb = 0.0;

    std::vector<PeriodRatioSummary> period_summaries;
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

static CombinationResult evaluate_case(const CsvTable& table,
                                       const CombinationCase& c,
                                       CentralValueMode mode) {
    CombinationResult result;
    result.label = c.label;
    result.central_value_mode = central_value_mode_name(mode);
    result.output_column = c.output_column;
    result.fill_output = c.fill_output;

    std::map<std::string, PeriodRatioAccumulator> acc_by_period;

    for (const auto& input : c.inputs) {
        acc_by_period[input.period] = PeriodRatioAccumulator();
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

        int n_valid_in_bin = 0;

        for (size_t i = 0; i < c.inputs.size(); ++i) {
            const TupleValue& v = input_values[i];

            if (!v.ok || !std::isfinite(v.value) || v.stat <= 0.0) {
                continue;
            }

            const double ratio = v.value / ref.value;
            const double ratio_stat = std::abs(v.stat / ref.value);

            if (!std::isfinite(ratio) ||
                !std::isfinite(ratio_stat) ||
                ratio_stat <= 0.0) {
                continue;
            }

            const double w = 1.0 / (ratio_stat * ratio_stat);

            PeriodRatioAccumulator& acc = acc_by_period[c.inputs[i].period];
            acc.sum_w += w;
            acc.sum_wr += w * ratio;
            acc.n += 1;

            ++result.ratio_points;
            ++n_valid_in_bin;
        }

        if (n_valid_in_bin >= 2) {
            ++result.valid_bins;
        }
    }

    double sum_obs2 = 0.0;
    double sum_stat2 = 0.0;
    int n_period = 0;

    for (const auto& input : c.inputs) {
        const PeriodRatioAccumulator& acc = acc_by_period[input.period];

        if (acc.sum_w <= 0.0 || acc.n <= 0) {
            continue;
        }

        PeriodRatioSummary p;
        p.period = input.period;
        p.n = acc.n;
        p.mean_ratio = acc.sum_wr / acc.sum_w;
        p.mean_ratio_stat = 1.0 / std::sqrt(acc.sum_w);

        if (!std::isfinite(p.mean_ratio) ||
            !std::isfinite(p.mean_ratio_stat) ||
            p.mean_ratio_stat <= 0.0) {
            continue;
        }

        const double residual = p.mean_ratio - 1.0;
        sum_obs2 += residual * residual;
        sum_stat2 += p.mean_ratio_stat * p.mean_ratio_stat;
        ++n_period;

        result.period_summaries.push_back(p);
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

    fout << "case,central value mode,output column,fill output,valid bins,ratio points,"
         << "s_obs_period,s_stat_period,s_comb,s_comb percent,"
         << "period,period ratio points,period mean ratio,period mean ratio stat\n";

    for (const auto& r : results) {
        if (r.period_summaries.empty()) {
            fout << csv_escape_field(r.label) << ","
                 << csv_escape_field(r.central_value_mode) << ","
                 << csv_escape_field(r.output_column) << ","
                 << (r.fill_output ? "true" : "false") << ","
                 << r.valid_bins << ","
                 << r.ratio_points << ","
                 << format_double(r.s_obs) << ","
                 << format_double(r.s_stat_exp) << ","
                 << format_double(r.s_comb) << ","
                 << format_double(100.0 * r.s_comb) << ","
                 << ",,,\n";
            continue;
        }

        for (const auto& p : r.period_summaries) {
            fout << csv_escape_field(r.label) << ","
                 << csv_escape_field(r.central_value_mode) << ","
                 << csv_escape_field(r.output_column) << ","
                 << (r.fill_output ? "true" : "false") << ","
                 << r.valid_bins << ","
                 << r.ratio_points << ","
                 << format_double(r.s_obs) << ","
                 << format_double(r.s_stat_exp) << ","
                 << format_double(r.s_comb) << ","
                 << format_double(100.0 * r.s_comb) << ","
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
    std::cout << "\n[combination-systematics] Summary\n";

    std::cout << std::left
              << std::setw(30) << "case"
              << std::setw(18) << "central"
              << std::right
              << std::setw(10) << "bins"
              << std::setw(12) << "ratios"
              << std::setw(14) << "s_obs"
              << std::setw(14) << "s_stat"
              << std::setw(14) << "s_comb"
              << std::setw(14) << "percent"
              << "\n";

    for (const auto& r : results) {
        std::cout << std::left
                  << std::setw(30) << r.label
                  << std::setw(18) << r.central_value_mode
                  << std::right
                  << std::setw(10) << r.valid_bins
                  << std::setw(12) << r.ratio_points
                  << std::setw(14) << std::setprecision(6) << r.s_obs
                  << std::setw(14) << std::setprecision(6) << r.s_stat_exp
                  << std::setw(14) << std::setprecision(6) << r.s_comb
                  << std::setw(14) << std::setprecision(6) << 100.0 * r.s_comb
                  << "\n";

        for (const auto& p : r.period_summaries) {
            std::cout << "    period mean ratio "
                      << std::left << std::setw(10) << p.period
                      << " = " << std::right
                      << std::setprecision(8) << p.mean_ratio
                      << " +/- " << std::setprecision(8) << p.mean_ratio_stat
                      << "   n=" << p.n
                      << "\n";
        }
    }

    std::cout << "\n";
}

static bool is_fill_mode(CentralValueMode mode) {
    return mode == kFillOutputMode;
}

} // namespace

bool combination_systematics(const std::string& csv_path,
                             const std::string& output_root_dir) {
    try {
        const std::string out_dir = output_root_dir + "/combination_systematics";
        fs::create_directories(out_dir);

        CsvTable table = read_csv_or_throw(csv_path);
        const std::vector<CombinationCase> cases = combination_cases();

        validate_schema(table, cases);

        std::cout << "[combination-systematics] CSV rows loaded: "
                  << table.rows.size() << "\n";
        std::cout << "[combination-systematics] Output directory: "
                  << out_dir << "\n";
        std::cout << "[combination-systematics] CSV output columns filled using central mode: "
                  << central_value_mode_name(kFillOutputMode) << "\n";

        std::vector<CombinationResult> results;
        results.reserve(cases.size() * central_value_modes().size() + 1);

        double ten6_unpol_s_comb_for_fill =
            std::numeric_limits<double>::quiet_NaN();

        for (const auto& mode : central_value_modes()) {
            for (const auto& c : cases) {
                CombinationResult result = evaluate_case(table, c, mode);

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
        results.push_back(sp19);

        const std::string summary_csv =
            out_dir + "/combination_systematics_summary.csv";

        write_summary_csv(summary_csv, results);
        print_summary_table(results);
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