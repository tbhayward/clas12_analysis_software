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

struct CombinationCase {
    std::string label;
    std::string output_column;
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
    std::string output_column;
    int valid_bins = 0;
    int ratio_points = 0;
    double s_obs = 0.0;
    double s_stat_exp = 0.0;
    double s_comb = 0.0;
    std::vector<PeriodRatioSummary> period_summaries;
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

static std::string join_csv_row(const std::vector<std::string>& fields) {
    std::ostringstream oss;

    for (size_t i = 0; i < fields.size(); ++i) {
        const std::string& s = fields[i];

        const bool need_quotes =
            (s.find(',') != std::string::npos) ||
            (s.find('"') != std::string::npos) ||
            (s.find('\n') != std::string::npos) ||
            (s.find('\r') != std::string::npos);

        if (need_quotes) {
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
            throw std::runtime_error("Internal CSV row/header size mismatch before writing.");
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

static std::string combination_column(const std::string& label,
                                      const std::string& helicity) {
    return "normed cross sections, ep->epg, exp, " + label + ", " + helicity + ", combination sys";
}

static std::vector<CombinationCase> combination_cases() {
    return {
        {
            "10.6 GeV unpol",
            combination_column("10.6 GeV", "unpol"),
            {
                {"Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")},
                {"Fa18 Out", cross_section_column("Fa18 Out", "unpol")},
                {"Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")},
                {"Sp18 Out", cross_section_column("Sp18 Out", "unpol")}
            }
        },
        {
            "Fa18 unpol",
            combination_column("Fa18", "unpol"),
            {
                {"Fa18 Inb", cross_section_column("Fa18 Inb", "unpol")},
                {"Fa18 Out", cross_section_column("Fa18 Out", "unpol")}
            }
        },
        {
            "Fa18 pos",
            combination_column("Fa18", "pos"),
            {
                {"Fa18 Inb", cross_section_column("Fa18 Inb", "pos")},
                {"Fa18 Out", cross_section_column("Fa18 Out", "pos")}
            }
        },
        {
            "Fa18 neg",
            combination_column("Fa18", "neg"),
            {
                {"Fa18 Inb", cross_section_column("Fa18 Inb", "neg")},
                {"Fa18 Out", cross_section_column("Fa18 Out", "neg")}
            }
        },
        {
            "Sp18 unpol",
            combination_column("Sp18", "unpol"),
            {
                {"Sp18 Inb", cross_section_column("Sp18 Inb", "unpol")},
                {"Sp18 Out", cross_section_column("Sp18 Out", "unpol")}
            }
        },
        {
            "Fa18 vs Sp18 unpol",
            "",
            {
                {"Fa18", cross_section_column("Fa18", "unpol")},
                {"Sp18", cross_section_column("Sp18", "unpol")}
            }
        },
        {
            "Inb vs Out unpol",
            "",
            {
                {"Inb", cross_section_column("Inb", "unpol")},
                {"Out", cross_section_column("Out", "unpol")}
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
        if (!c.output_column.empty()) {
            required.push_back(c.output_column);
        }

        for (const auto& input : c.inputs) {
            required.push_back(input.column);
        }
    }

    for (const auto& col : sp19_proxy_output_columns()) {
        required.push_back(col);
    }

    require_columns(table, required, "combination systematics");
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
        (s.find(',') != std::string::npos) ||
        (s.find('"') != std::string::npos) ||
        (s.find('\n') != std::string::npos) ||
        (s.find('\r') != std::string::npos);

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
                                       const CombinationCase& c) {
    CombinationResult result;
    result.label = c.label;
    result.output_column = c.output_column;

    std::map<std::string, PeriodRatioAccumulator> acc_by_period;

    int n_valid_bins = 0;
    int n_ratio = 0;

    for (const auto& input : c.inputs) {
        acc_by_period[input.period] = PeriodRatioAccumulator();
    }

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

        ++n_valid_bins;

        for (size_t i = 0; i < values.size(); ++i) {
            const TupleValue& v = values[i];

            const double ratio = v.value / mean;
            const double ratio_stat = std::abs(v.stat / mean);

            if (!std::isfinite(ratio) || !std::isfinite(ratio_stat) || ratio_stat <= 0.0) {
                continue;
            }

            const double w = 1.0 / (ratio_stat * ratio_stat);

            PeriodRatioAccumulator& acc = acc_by_period[c.inputs[i].period];
            acc.sum_w += w;
            acc.sum_wr += w * ratio;
            acc.n += 1;

            ++n_ratio;
        }
    }

    result.valid_bins = n_valid_bins;
    result.ratio_points = n_ratio;

    double sum_obs_period2 = 0.0;
    double sum_stat_period2 = 0.0;
    int n_period = 0;

    for (const auto& input : c.inputs) {
        const auto it = acc_by_period.find(input.period);
        if (it == acc_by_period.end()) {
            continue;
        }

        const PeriodRatioAccumulator& acc = it->second;
        if (acc.sum_w <= 0.0 || acc.n <= 0) {
            continue;
        }

        PeriodRatioSummary summary;
        summary.period = input.period;
        summary.n = acc.n;
        summary.mean_ratio = acc.sum_wr / acc.sum_w;
        summary.mean_ratio_stat = 1.0 / std::sqrt(acc.sum_w);

        if (!std::isfinite(summary.mean_ratio) ||
            !std::isfinite(summary.mean_ratio_stat) ||
            summary.mean_ratio_stat <= 0.0) {
            continue;
        }

        const double residual = summary.mean_ratio - 1.0;

        sum_obs_period2 += residual * residual;
        sum_stat_period2 += summary.mean_ratio_stat * summary.mean_ratio_stat;
        ++n_period;

        result.period_summaries.push_back(summary);
    }

    if (n_period > 0) {
        const double s_obs2 = sum_obs_period2 / (double)n_period;
        const double s_stat2 = sum_stat_period2 / (double)n_period;
        const double s_comb2 = std::max(0.0, s_obs2 - s_stat2);

        result.s_obs = std::sqrt(std::max(0.0, s_obs2));
        result.s_stat_exp = std::sqrt(std::max(0.0, s_stat2));
        result.s_comb = std::sqrt(s_comb2);
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
        throw std::runtime_error("Could not open combination-systematics summary CSV: " + path);
    }

    fout << "case,output column,valid bins,ratio points,s_obs_period,s_stat_period,s_comb,s_comb percent,period,period ratio points,period mean ratio,period mean ratio stat\n";

    for (const auto& r : results) {
        if (r.period_summaries.empty()) {
            fout << csv_escape_field(r.label) << ","
                 << csv_escape_field(r.output_column) << ","
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
                 << csv_escape_field(r.output_column) << ","
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
        throw std::runtime_error("Failed while writing combination-systematics summary CSV: " + path);
    }
}

static void print_summary_table(const std::vector<CombinationResult>& results) {
    std::cout << "\n";
    std::cout << "[combination-systematics] Summary\n";
    std::cout << std::left
              << std::setw(24) << "case"
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
                  << std::setw(24) << r.label
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
                      << " = " << std::right << std::setprecision(8) << p.mean_ratio
                      << " +/- " << std::setprecision(8) << p.mean_ratio_stat
                      << "   n=" << p.n
                      << "\n";
        }
    }

    std::cout << "\n";
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

        std::vector<CombinationResult> results;
        results.reserve(cases.size() + 1);

        double ten6_unpol_s_comb = std::numeric_limits<double>::quiet_NaN();

        for (const auto& c : cases) {
            CombinationResult result = evaluate_case(table, c);

            if (!result.output_column.empty()) {
                fill_output_column(table, result.output_column, result.s_comb);
            }

            if (c.label == "10.6 GeV unpol") {
                ten6_unpol_s_comb = result.s_comb;
            }

            results.push_back(result);
        }

        if (!std::isfinite(ten6_unpol_s_comb)) {
            throw std::runtime_error("Could not determine 10.6 GeV unpol s_comb for Sp19 Inb proxy fill.");
        }

        for (const auto& col : sp19_proxy_output_columns()) {
            fill_output_column(table, col, ten6_unpol_s_comb);
        }

        CombinationResult sp19_proxy;
        sp19_proxy.label = "Sp19 Inb unpol";
        sp19_proxy.output_column = combination_column("Sp19 Inb", "unpol");
        sp19_proxy.valid_bins = 0;
        sp19_proxy.ratio_points = 0;
        sp19_proxy.s_obs = ten6_unpol_s_comb;
        sp19_proxy.s_stat_exp = 0.0;
        sp19_proxy.s_comb = ten6_unpol_s_comb;
        results.push_back(sp19_proxy);

        const std::string summary_csv = out_dir + "/combination_systematics_summary.csv";
        write_summary_csv(summary_csv, results);
        print_summary_table(results);

        write_csv_or_throw(csv_path, table);

        std::cout << "[combination-systematics] Wrote summary CSV: "
                  << summary_csv << "\n";
        std::cout << "[combination-systematics] Updated CSV: "
                  << csv_path << "\n";

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[combination-systematics] FATAL: " << e.what() << "\n";
        return false;
    }
}