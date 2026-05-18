#include "combination_point_to_point_systematics.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TH1.h>
#include <TLine.h>
#include <TF1.h>

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

using Range = std::pair<double, double>;
using QTKey = std::pair<Range, Range>;

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
    std::string column;
};

struct PeriodScale {
    std::string period;
    double scale = 1.0;
    double scale_stat = 0.0;
    int n = 0;
};

struct PeriodRatioAccumulator {
    double sum_w = 0.0;
    double sum_wr = 0.0;
    int n = 0;
};

struct BinResult {
    size_t row_index = 0;

    Range xb_range;
    Range q2_range;
    Range t_range;

    int xb_index = -1;

    double xBavg = std::numeric_limits<double>::quiet_NaN();
    double Q2avg = std::numeric_limits<double>::quiet_NaN();
    double t_abs_avg = std::numeric_limits<double>::quiet_NaN();
    double phiavg = std::numeric_limits<double>::quiet_NaN();
    double e_theta = std::numeric_limits<double>::quiet_NaN();
    double p_theta = std::numeric_limits<double>::quiet_NaN();
    double g_theta = std::numeric_limits<double>::quiet_NaN();

    double mean_scaled = std::numeric_limits<double>::quiet_NaN();
    double mean_scaled_stat = std::numeric_limits<double>::quiet_NaN();

    double chi2_ndf_before = std::numeric_limits<double>::quiet_NaN();
    double chi2_ndf_after = std::numeric_limits<double>::quiet_NaN();

    double f_ptp = std::numeric_limits<double>::quiet_NaN();
    double f_ptp_percent = std::numeric_limits<double>::quiet_NaN();
    double f_ptp_stat_percent = std::numeric_limits<double>::quiet_NaN();

    int n_periods = 0;
};

struct PlotPoint {
    double phi = 0.0;
    double y = 0.0;
};

struct BinData {
    std::vector<PlotPoint> points;
};

struct GroupByXB {
    std::map<QTKey, BinData> bins;
    Range xb_range;
    int xb_index = -1;
};

struct VariableScatterConfig {
    std::string key;
    std::string title;
    bool angle_canvas = false;
};

static const std::vector<std::string> kBasePeriods = {
    "Fa18 Inb",
    "Fa18 Out",
    "Sp18 Inb",
    "Sp18 Out"
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

static std::string csv_escape_field(const std::string& s) {
    bool need_quotes = false;

    for (const char c : s) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            need_quotes = true;
            break;
        }
    }

    if (!need_quotes) {
        return s;
    }

    std::ostringstream oss;
    oss << '"';

    for (const char c : s) {
        if (c == '"') {
            oss << "\"\"";
        } else {
            oss << c;
        }
    }

    oss << '"';
    return oss.str();
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
    if (fields.empty()) {
        return out;
    }

    double value = 0.0;
    if (!parse_double(fields[0], value)) {
        return out;
    }

    double stat = 0.0;
    if (fields.size() > 1U) {
        if (!parse_double(fields[1], stat)) {
            return out;
        }
    } else {
        return out;
    }

    double sys = 0.0;
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

static int get_col(const CsvTable& table,
                   const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) {
        throw std::runtime_error("Missing required column: " + column);
    }

    return it->second;
}

static int get_col_optional(const CsvTable& table,
                            const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) {
        return -1;
    }

    return it->second;
}

static double get_double_or_nan(const CsvTable& table,
                                const std::vector<std::string>& row,
                                const std::string& column) {
    const int col = get_col(table, column);

    double value = 0.0;
    if (!parse_double(row[(size_t)col], value)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return value;
}

static TupleValue get_tuple(const CsvTable& table,
                            const std::vector<std::string>& row,
                            const std::string& column) {
    const int col = get_col(table, column);
    return parse_tuple_value(row[(size_t)col]);
}

static std::string cross_section_column(const std::string& period,
                                        const std::string& helicity) {
    return "normed cross sections, ep->epg, exp, " + period + ", " + helicity;
}

static std::vector<PeriodInput> base_inputs() {
    std::vector<PeriodInput> inputs;

    for (const auto& period : kBasePeriods) {
        PeriodInput input;
        input.period = period;
        input.column = cross_section_column(period, "unpol");
        inputs.push_back(input);
    }

    return inputs;
}

static void validate_schema(const CsvTable& table) {
    std::vector<std::string> required = {
        "xBmin",
        "xBmax",
        "Q2min",
        "Q2max",
        "t_abs_min",
        "t_abs_max",
        "phimin",
        "phimax",
        "xBavg, 10.6 GeV",
        "Q2avg, 10.6 GeV",
        "t_abs_avg, 10.6 GeV",
        "phiavg, 10.6 GeV",
        "e_theta, 10.6 GeV",
        "p_theta, 10.6 GeV",
        "g_theta, 10.6 GeV"
    };

    for (const auto& input : base_inputs()) {
        required.push_back(input.column);
    }

    require_columns(table, required, "combination point-to-point systematics");
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

static std::vector<PeriodScale> compute_10p6_unpol_scale_factors(const CsvTable& table) {
    const std::vector<PeriodInput> inputs = base_inputs();

    std::map<std::string, PeriodRatioAccumulator> acc_by_period;
    for (const auto& input : inputs) {
        acc_by_period[input.period] = PeriodRatioAccumulator();
    }

    for (const auto& row : table.rows) {
        std::vector<TupleValue> values;
        values.reserve(inputs.size());

        for (const auto& input : inputs) {
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

        for (size_t i = 0; i < inputs.size(); ++i) {
            const TupleValue& v = values[i];

            if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) {
                continue;
            }

            const double ratio = v.value / mean;
            const double ratio_stat = std::abs(v.stat / mean);

            if (!std::isfinite(ratio) || !std::isfinite(ratio_stat) || ratio_stat <= 0.0) {
                continue;
            }

            const double w = 1.0 / (ratio_stat * ratio_stat);

            PeriodRatioAccumulator& acc = acc_by_period[inputs[i].period];
            acc.sum_w += w;
            acc.sum_wr += w * ratio;
            acc.n += 1;
        }
    }

    std::vector<PeriodScale> scales;
    scales.reserve(inputs.size());

    std::cout << "[combination-ptp] 10.6 GeV unpol scale factors used before point-to-point estimate:\n";

    for (const auto& input : inputs) {
        const PeriodRatioAccumulator& acc = acc_by_period[input.period];

        if (acc.sum_w <= 0.0 || acc.n <= 0) {
            throw std::runtime_error("Could not compute scale factor for period: " + input.period);
        }

        PeriodScale scale;
        scale.period = input.period;
        scale.scale = acc.sum_wr / acc.sum_w;
        scale.scale_stat = 1.0 / std::sqrt(acc.sum_w);
        scale.n = acc.n;

        if (!std::isfinite(scale.scale) || std::abs(scale.scale) <= 0.0) {
            throw std::runtime_error("Invalid scale factor for period: " + input.period);
        }

        std::cout << "  " << scale.period
                  << " scale = " << std::setprecision(10) << scale.scale
                  << " +/- " << scale.scale_stat
                  << "   n=" << scale.n << "\n";

        scales.push_back(scale);
    }

    return scales;
}

static std::map<std::string, double>
scale_map_from_vector(const std::vector<PeriodScale>& scales) {
    std::map<std::string, double> out;

    for (const auto& s : scales) {
        out[s.period] = s.scale;
    }

    return out;
}

static TupleValue apply_scale(const TupleValue& input,
                              double scale) {
    TupleValue out;

    if (!input.ok || !std::isfinite(scale) || std::abs(scale) <= 0.0) {
        return out;
    }

    out.ok = true;
    out.value = input.value / scale;
    out.stat = input.stat / std::abs(scale);
    out.sys = input.sys / std::abs(scale);

    return out;
}

static bool weighted_mean_with_extra_fraction(const std::vector<TupleValue>& values,
                                              double f,
                                              double reference_mean,
                                              double& mean,
                                              double& mean_stat) {
    if (!std::isfinite(f) || f < 0.0 || !std::isfinite(reference_mean)) {
        return false;
    }

    const double extra = f * std::abs(reference_mean);

    double sum_w = 0.0;
    double sum_wx = 0.0;

    for (const auto& v : values) {
        if (!v.ok || !std::isfinite(v.value) || !std::isfinite(v.stat) || v.stat <= 0.0) {
            continue;
        }

        const double var = v.stat * v.stat + extra * extra;
        if (!(var > 0.0) || !std::isfinite(var)) {
            continue;
        }

        const double w = 1.0 / var;
        sum_w += w;
        sum_wx += w * v.value;
    }

    if (sum_w <= 0.0) {
        return false;
    }

    mean = sum_wx / sum_w;
    mean_stat = 1.0 / std::sqrt(sum_w);

    return std::isfinite(mean) && std::isfinite(mean_stat);
}

static double chi2_ndf_with_extra_fraction(const std::vector<TupleValue>& values,
                                           double f,
                                           double reference_mean) {
    int n = 0;

    double mean = 0.0;
    double mean_stat = 0.0;

    if (!weighted_mean_with_extra_fraction(values, f, reference_mean, mean, mean_stat)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double extra = f * std::abs(reference_mean);

    double chi2 = 0.0;

    for (const auto& v : values) {
        if (!v.ok || !std::isfinite(v.value) || !std::isfinite(v.stat) || v.stat <= 0.0) {
            continue;
        }

        const double var = v.stat * v.stat + extra * extra;
        if (!(var > 0.0) || !std::isfinite(var)) {
            continue;
        }

        const double residual = v.value - mean;
        chi2 += residual * residual / var;
        ++n;
    }

    const double ndf = (double)n - 1.0;
    if (ndf <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return chi2 / ndf;
}

static double solve_fraction_for_chi2_unity(const std::vector<TupleValue>& scaled_values,
                                            double reference_mean,
                                            double chi2_initial) {
    if (!std::isfinite(chi2_initial)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (chi2_initial <= 1.0) {
        return 0.0;
    }

    double lo = 0.0;
    double hi = 0.01;

    for (int iter = 0; iter < 80; ++iter) {
        const double chi2_hi = chi2_ndf_with_extra_fraction(scaled_values, hi, reference_mean);

        if (std::isfinite(chi2_hi) && chi2_hi <= 1.0) {
            break;
        }

        hi *= 2.0;

        if (hi > 10.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    for (int iter = 0; iter < 100; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double chi2_mid = chi2_ndf_with_extra_fraction(scaled_values, mid, reference_mean);

        if (!std::isfinite(chi2_mid)) {
            hi = mid;
            continue;
        }

        if (chi2_mid > 1.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return hi;
}

static int optional_xb_index(const CsvTable& table,
                             const std::vector<std::string>& row) {
    const int col = get_col_optional(table, "xB index");
    if (col < 0) {
        return -1;
    }

    double value = 0.0;
    if (!parse_double(row[(size_t)col], value)) {
        return -1;
    }

    return (int)std::llround(value);
}

static BinResult evaluate_row(const CsvTable& table,
                              const std::vector<std::string>& row,
                              size_t row_index,
                              const std::map<std::string, double>& scale_by_period) {
    BinResult result;
    result.row_index = row_index;

    result.xb_range = Range(
        get_double_or_nan(table, row, "xBmin"),
        get_double_or_nan(table, row, "xBmax")
    );

    result.q2_range = Range(
        get_double_or_nan(table, row, "Q2min"),
        get_double_or_nan(table, row, "Q2max")
    );

    result.t_range = Range(
        get_double_or_nan(table, row, "t_abs_min"),
        get_double_or_nan(table, row, "t_abs_max")
    );

    result.xb_index = optional_xb_index(table, row);

    result.xBavg = get_double_or_nan(table, row, "xBavg, 10.6 GeV");
    result.Q2avg = get_double_or_nan(table, row, "Q2avg, 10.6 GeV");
    result.t_abs_avg = get_double_or_nan(table, row, "t_abs_avg, 10.6 GeV");
    result.phiavg = get_double_or_nan(table, row, "phiavg, 10.6 GeV");

    if (!std::isfinite(result.phiavg)) {
        const double phimin = get_double_or_nan(table, row, "phimin");
        const double phimax = get_double_or_nan(table, row, "phimax");

        if (std::isfinite(phimin) && std::isfinite(phimax)) {
            result.phiavg = 0.5 * (phimin + phimax);
        }
    }

    result.e_theta = get_double_or_nan(table, row, "e_theta, 10.6 GeV");
    result.p_theta = get_double_or_nan(table, row, "p_theta, 10.6 GeV");
    result.g_theta = get_double_or_nan(table, row, "g_theta, 10.6 GeV");

    std::vector<TupleValue> scaled_values;

    for (const auto& input : base_inputs()) {
        const TupleValue raw = get_tuple(table, row, input.column);

        const auto it_scale = scale_by_period.find(input.period);
        if (it_scale == scale_by_period.end()) {
            throw std::runtime_error("Missing scale factor for period: " + input.period);
        }

        const TupleValue scaled = apply_scale(raw, it_scale->second);

        if (scaled.ok) {
            scaled_values.push_back(scaled);
        }
    }

    result.n_periods = (int)scaled_values.size();
    if (result.n_periods < 2) {
        return result;
    }

    double mean = 0.0;
    double mean_stat_value = 0.0;

    if (!compute_weighted_mean(scaled_values, mean, mean_stat_value)) {
        return result;
    }

    if (!std::isfinite(mean) || std::abs(mean) <= 0.0) {
        return result;
    }

    result.mean_scaled = mean;
    result.mean_scaled_stat = mean_stat_value;

    result.chi2_ndf_before =
        chi2_ndf_with_extra_fraction(scaled_values, 0.0, result.mean_scaled);

    result.f_ptp =
        solve_fraction_for_chi2_unity(scaled_values,
                                      result.mean_scaled,
                                      result.chi2_ndf_before);

    if (std::isfinite(result.f_ptp)) {
        result.chi2_ndf_after =
            chi2_ndf_with_extra_fraction(scaled_values,
                                         result.f_ptp,
                                         result.mean_scaled);
        result.f_ptp_percent = 100.0 * result.f_ptp;

        if (std::isfinite(result.mean_scaled_stat) && std::abs(result.mean_scaled) > 0.0) {
            result.f_ptp_stat_percent =
                100.0 * result.mean_scaled_stat / std::abs(result.mean_scaled);
        }
    }

    return result;
}

static std::vector<BinResult>
evaluate_all_bins(const CsvTable& table,
                  const std::map<std::string, double>& scale_by_period) {
    std::vector<BinResult> results;
    results.reserve(table.rows.size());

    for (size_t i = 0; i < table.rows.size(); ++i) {
        BinResult result = evaluate_row(table, table.rows[i], i + 1, scale_by_period);

        if (!std::isfinite(result.f_ptp)) {
            continue;
        }

        if (!std::isfinite(result.phiavg)) {
            continue;
        }

        results.push_back(result);
    }

    return results;
}

static std::string format_double(double value) {
    if (!std::isfinite(value)) {
        return "";
    }

    std::ostringstream oss;
    oss << std::setprecision(12) << value;
    return oss.str();
}

static void write_scale_summary_csv(const std::string& path,
                                    const std::vector<PeriodScale>& scales) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open scale summary CSV for writing: " + path);
    }

    fout << "period,scale,scale stat,n\n";

    for (const auto& s : scales) {
        fout << csv_escape_field(s.period) << ","
             << format_double(s.scale) << ","
             << format_double(s.scale_stat) << ","
             << s.n << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing scale summary CSV: " + path);
    }
}

static void write_bin_summary_csv(const std::string& path,
                                  const std::vector<BinResult>& results) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open point-to-point summary CSV for writing: " + path);
    }

    fout << "row index,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,"
         << "xBavg,Q2avg,t_abs_avg,phiavg,e_theta,p_theta,g_theta,n periods,"
         << "mean scaled,mean scaled stat,chi2ndf before,chi2ndf after,"
         << "f_ptp,f_ptp percent,f_ptp stat percent\n";

    for (const auto& r : results) {
        fout << r.row_index << ","
             << format_double(r.xb_range.first) << ","
             << format_double(r.xb_range.second) << ","
             << format_double(r.q2_range.first) << ","
             << format_double(r.q2_range.second) << ","
             << format_double(r.t_range.first) << ","
             << format_double(r.t_range.second) << ","
             << format_double(r.xBavg) << ","
             << format_double(r.Q2avg) << ","
             << format_double(r.t_abs_avg) << ","
             << format_double(r.phiavg) << ","
             << format_double(r.e_theta) << ","
             << format_double(r.p_theta) << ","
             << format_double(r.g_theta) << ","
             << r.n_periods << ","
             << format_double(r.mean_scaled) << ","
             << format_double(r.mean_scaled_stat) << ","
             << format_double(r.chi2_ndf_before) << ","
             << format_double(r.chi2_ndf_after) << ","
             << format_double(r.f_ptp) << ","
             << format_double(r.f_ptp_percent) << ","
             << format_double(r.f_ptp_stat_percent) << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing point-to-point summary CSV: " + path);
    }
}

static std::map<Range, GroupByXB>
group_results_by_xb(const std::vector<BinResult>& results) {
    std::map<Range, GroupByXB> grouped;

    for (const auto& r : results) {
        GroupByXB& group = grouped[r.xb_range];
        group.xb_range = r.xb_range;

        if (group.xb_index < 0 && r.xb_index >= 0) {
            group.xb_index = r.xb_index;
        }

        QTKey key(r.q2_range, r.t_range);

        PlotPoint p;
        p.phi = r.phiavg;
        p.y = r.f_ptp_percent;

        group.bins[key].points.push_back(p);
    }

    return grouped;
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

static double choose_global_ymax_percent(const std::vector<BinResult>& results) {
    std::vector<double> ys;

    for (const auto& r : results) {
        if (std::isfinite(r.f_ptp_percent) && r.f_ptp_percent >= 0.0) {
            ys.push_back(r.f_ptp_percent);
        }
    }

    if (ys.empty()) {
        return 10.0;
    }

    const double p95 = percentile(ys, 0.95);
    const double max_y = *std::max_element(ys.begin(), ys.end());

    double ymax = std::max(5.0, 1.25 * p95);

    if (max_y > ymax && max_y < 2.0 * ymax) {
        ymax = 1.15 * max_y;
    }

    if (!std::isfinite(ymax) || ymax <= 0.0) {
        ymax = 10.0;
    }

    return ymax;
}

static TGraph* make_ptp_graph(const std::vector<PlotPoint>& points) {
    if (points.empty()) {
        return nullptr;
    }

    std::vector<PlotPoint> sorted = points;
    std::sort(sorted.begin(), sorted.end(),
              [](const PlotPoint& a, const PlotPoint& b) {
                  return a.phi < b.phi;
              });

    const int n = (int)sorted.size();
    std::vector<double> x((size_t)n);
    std::vector<double> y((size_t)n);

    for (int i = 0; i < n; ++i) {
        x[(size_t)i] = sorted[(size_t)i].phi;
        y[(size_t)i] = sorted[(size_t)i].y;
    }

    TGraph* graph = new TGraph(n, x.data(), y.data());
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.85);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlack);
    graph->SetLineWidth(1);

    return graph;
}

static void draw_one_xb_canvas(const Range& xb_range,
                               const GroupByXB& group,
                               int xb_counter,
                               double y_max_percent,
                               const fs::path& out_dir) {
    if (group.bins.empty()) {
        return;
    }

    std::set<Range> q2_set;
    std::set<Range> t_set;

    for (const auto& kv : group.bins) {
        q2_set.insert(kv.first.first);
        t_set.insert(kv.first.second);
    }

    std::vector<Range> q2_slice(q2_set.begin(), q2_set.end());
    std::vector<Range> t_slice(t_set.begin(), t_set.end());

    if (q2_slice.empty() || t_slice.empty()) {
        return;
    }

    const int ncols = (int)q2_slice.size();
    const int nrows = (int)t_slice.size();

    const int width = 300 * ncols + 160;
    const int height = 260 * nrows + 240;

    TCanvas* canvas = new TCanvas("c_combination_ptp",
                                  "c_combination_ptp",
                                  width,
                                  height);

    TPad* pTop = new TPad("pTopCombinationPtp", "pTopCombinationPtp",
                          0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad* pGrid = new TPad("pGridCombinationPtp", "pGridCombinationPtp",
                           0.0, 0.00, 1.0, 0.78);
    pGrid->SetFillStyle(0);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(0.060);

    std::ostringstream title;
    title << "Point-to-point run-period combination systematic, 10.6 GeV unpol"
          << "   x_{B} in ("
          << std::fixed << std::setprecision(3)
          << xb_range.first << ", " << xb_range.second << ")";

    head.DrawLatex(0.5, 0.83, title.str().c_str());

    TLatex sub;
    sub.SetNDC();
    sub.SetTextAlign(22);
    sub.SetTextFont(42);
    sub.SetTextSize(0.040);
    sub.DrawLatex(0.5, 0.50,
                  "Global period scale offsets removed first; plotted value is the residual fractional uncertainty needed to make #chi^{2}/ndf = 1");

    for (int r = 0; r < nrows; ++r) {
        const Range& t_range = t_slice[(size_t)r];

        for (int c = 0; c < ncols; ++c) {
            const Range& q2_range = q2_slice[(size_t)c];

            pGrid->cd(r * ncols + c + 1);

            gPad->SetGrid(1, 1);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.07);
            gPad->SetTicks(1, 1);

            TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, y_max_percent);

            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("Combination point-to-point sys (%)");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.070);
            frame->GetYaxis()->SetTitleSize(0.070);
            frame->GetXaxis()->SetLabelSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.060);
            frame->GetXaxis()->SetTitleOffset(1.00);
            frame->GetYaxis()->SetTitleOffset(1.15);

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(0.060);
            lab.SetTextAlign(11);
            lab.SetTextFont(42);

            lab.DrawLatex(
                0.12, 0.83,
                Form("Q^{2}=(%.2f, %.2f)  |t|=(%.2f, %.2f)",
                     q2_range.first, q2_range.second,
                     t_range.first, t_range.second)
            );

            const QTKey key(q2_range, t_range);
            const auto it_bin = group.bins.find(key);

            if (it_bin == group.bins.end()) {
                continue;
            }

            TGraph* graph = make_ptp_graph(it_bin->second.points);
            if (graph) {
                graph->Draw("P SAME");
            }
        }
    }

    canvas->cd();
    canvas->Modified();
    canvas->Update();

    const int xb_name =
        (group.xb_index >= 0 ? group.xb_index : xb_counter);

    std::ostringstream fname;
    fname << "combination_point_to_point_sys_10p6_GeV_unpol_xB_"
          << xb_name << ".pdf";

    const fs::path out_path = out_dir / fname.str();
    canvas->SaveAs(out_path.string().c_str());

    delete canvas;
}

static std::vector<VariableScatterConfig> physics_scatter_variables() {
    return {
        {"xB", "x_{B}", false},
        {"Q2", "Q^{2} (GeV^{2})", false},
        {"t", "-t (GeV^{2})", false},
        {"phi", "#phi (deg)", false}
    };
}

static std::vector<VariableScatterConfig> angle_scatter_variables() {
    return {
        {"e_theta", "#theta_{e} (deg)", true},
        {"p_theta", "#theta_{p} (deg)", true},
        {"g_theta", "#theta_{#gamma} (deg)", true}
    };
}

static double value_for_scatter_variable(const BinResult& r,
                                         const std::string& key) {
    if (key == "xB") return r.xBavg;
    if (key == "Q2") return r.Q2avg;
    if (key == "t") return r.t_abs_avg;
    if (key == "phi") return r.phiavg;
    if (key == "e_theta") return r.e_theta;
    if (key == "p_theta") return r.p_theta;
    if (key == "g_theta") return r.g_theta;

    return std::numeric_limits<double>::quiet_NaN();
}

static std::vector<double> choose_scatter_x_range(const std::vector<double>& xs) {
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

static TGraph* make_scatter_graph(const std::vector<double>& x,
                                  const std::vector<double>& y) {
    if (x.empty() || y.empty() || x.size() != y.size()) {
        return nullptr;
    }

    TGraph* graph = new TGraph((int)x.size());

    for (int i = 0; i < (int)x.size(); ++i) {
        graph->SetPoint(i, x[(size_t)i], y[(size_t)i]);
    }

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.55);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlack);
    graph->SetLineWidth(1);

    return graph;
}

static TGraphErrors* make_scatter_graph_errors(const std::vector<double>& x,
                                               const std::vector<double>& y,
                                               const std::vector<double>& ey) {
    if (x.empty() || y.empty() || ey.empty()) {
        return nullptr;
    }

    if (x.size() != y.size() || x.size() != ey.size()) {
        return nullptr;
    }

    TGraphErrors* graph = new TGraphErrors((int)x.size());

    for (int i = 0; i < (int)x.size(); ++i) {
        graph->SetPoint(i, x[(size_t)i], y[(size_t)i]);
        graph->SetPointError(i, 0.0, ey[(size_t)i]);
    }

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.55);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlack);
    graph->SetLineWidth(1);

    return graph;
}

static std::string poly_formula(int order) {
    if (order < 0 || order > 5) {
        throw std::runtime_error("Invalid polynomial order requested.");
    }

    if (order == 0) return "[0]";
    if (order == 1) return "[0] + [1]*x";
    if (order == 2) return "[0] + [1]*x + [2]*x*x";
    if (order == 3) return "[0] + [1]*x + [2]*x*x + [3]*x*x*x";
    if (order == 4) return "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x";

    return "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x";
}

static int poly_color(int order) {
    if (order == 0) return kRed + 1;
    if (order == 1) return kBlue + 1;
    if (order == 2) return kGreen + 2;
    if (order == 3) return kMagenta + 1;
    if (order == 4) return kOrange + 7;

    return kCyan + 2;
}

static int poly_line_style(int order) {
    if (order == 0) return 2;
    if (order == 1) return 3;
    if (order == 2) return 4;
    if (order == 3) return 5;
    if (order == 4) return 6;

    return 7;
}

static void fit_and_draw_polynomials(TGraphErrors* graph,
                                     const std::string& fit_tag,
                                     const std::string& variable_key,
                                     double xmin,
                                     double xmax,
                                     TLegend* legend,
                                     std::vector<std::unique_ptr<TF1> >& fit_functions) {
    if (!graph) {
        return;
    }

    const int n_points = graph->GetN();
    if (n_points < 2) {
        return;
    }

    std::cout << "[combination-ptp] Polynomial fits for "
              << fit_tag << " vs " << variable_key
              << " with n=" << n_points << "\n";

    for (int order = 0; order <= 5; ++order) {
        if (n_points <= order + 1) {
            std::cout << "    p" << order << ": skipped, not enough points for fit\n";
            continue;
        }

        const std::string name =
            "fit_" + fit_tag + "_" + variable_key + "_p" + std::to_string(order);

        std::unique_ptr<TF1> fit(new TF1(name.c_str(),
                                         poly_formula(order).c_str(),
                                         xmin,
                                         xmax));

        fit->SetLineColor(poly_color(order));
        fit->SetLineStyle(poly_line_style(order));
        fit->SetLineWidth(2);

        graph->Fit(fit.get(), "Q0");

        const double chi2 = fit->GetChisquare();
        const int ndf = fit->GetNDF();
        const double chi2_ndf =
            (ndf > 0 ? chi2 / (double)ndf : std::numeric_limits<double>::quiet_NaN());

        fit->Draw("SAME");

        std::ostringstream leg_label;
        leg_label << "p" << order << " #chi^{2}/ndf=";

        if (std::isfinite(chi2_ndf)) {
            leg_label << std::fixed << std::setprecision(2) << chi2_ndf;
        } else {
            leg_label << "n/a";
        }

        legend->AddEntry(fit.get(), leg_label.str().c_str(), "l");

        std::cout << "    p" << order
                  << ": chi2=" << std::setprecision(10) << chi2
                  << " ndf=" << ndf
                  << " chi2/ndf=" << chi2_ndf;

        for (int ip = 0; ip <= order; ++ip) {
            std::cout << " p" << ip << "=" << fit->GetParameter(ip)
                      << " +/- " << fit->GetParError(ip);
        }

        std::cout << "\n";

        fit_functions.push_back(std::move(fit));
    }
}

static void draw_scatter_canvas(const std::vector<BinResult>& results,
                                const std::vector<VariableScatterConfig>& vars,
                                const fs::path& out_dir,
                                const std::string& file_tag,
                                const std::string& title_tag,
                                int ncols,
                                int nrows,
                                double y_max_percent) {
    if (results.empty()) {
        return;
    }

    const int width = 650 * ncols;
    const int height = 520 * nrows + 140;
    TCanvas* canvas = new TCanvas(("c_ptp_scatter_" + file_tag).c_str(),
                                  ("c_ptp_scatter_" + file_tag).c_str(),
                                  width,
                                  height);

    TPad* pTop = new TPad(("pTop_" + file_tag).c_str(),
                          ("pTop_" + file_tag).c_str(),
                          0.0,
                          0.89,
                          1.0,
                          1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad* pGrid = new TPad(("pGrid_" + file_tag).c_str(),
                           ("pGrid_" + file_tag).c_str(),
                           0.0,
                           0.00,
                           1.0,
                           0.89);
    pGrid->SetFillStyle(0);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(0.075);
    head.DrawLatex(0.5, 0.58, title_tag.c_str());

    std::vector<std::unique_ptr<TGraph> > graphs;

    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        pGrid->cd(ivar + 1);

        gPad->SetGrid(1, 1);
        gPad->SetTopMargin(0.14);
        gPad->SetBottomMargin(0.18);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.07);
        gPad->SetTicks(1, 1);

        std::vector<double> xs;
        std::vector<double> ys;

        for (const auto& r : results) {
            const double x = value_for_scatter_variable(r, vars[(size_t)ivar].key);
            const double y = r.f_ptp_percent;

            if (!std::isfinite(x) || !std::isfinite(y)) {
                continue;
            }

            xs.push_back(x);
            ys.push_back(y);
        }

        const std::vector<double> xr = choose_scatter_x_range(xs);
        TH1* frame = gPad->DrawFrame(xr[0], 0.0, xr[1], y_max_percent);
        frame->GetXaxis()->SetTitle(vars[(size_t)ivar].title.c_str());
        frame->GetYaxis()->SetTitle("Combination point-to-point sys (%)");
        frame->GetXaxis()->CenterTitle();
        frame->GetYaxis()->CenterTitle();
        frame->GetXaxis()->SetNdivisions(505);
        frame->GetXaxis()->SetTitleSize(0.060);
        frame->GetYaxis()->SetTitleSize(0.060);
        frame->GetXaxis()->SetLabelSize(0.050);
        frame->GetYaxis()->SetLabelSize(0.050);
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitleOffset(1.25);

        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(42);
        lab.SetTextSize(0.060);
        lab.SetTextAlign(22);
        lab.DrawLatex(0.50, 0.93, vars[(size_t)ivar].title.c_str());

        TGraph* graph = make_scatter_graph(xs, ys);
        if (graph) {
            graph->Draw("P SAME");
            graphs.emplace_back(graph);
        }
    }

    canvas->cd();
    canvas->Modified();
    canvas->Update();

    canvas->SaveAs((out_dir / ("combination_point_to_point_sys_" + file_tag + ".pdf")).string().c_str());

    delete canvas;
}

static void draw_scatter_canvas_with_errors_and_fits(
    const std::vector<BinResult>& results,
    const std::vector<VariableScatterConfig>& vars,
    const fs::path& out_dir,
    const std::string& file_tag,
    const std::string& title_tag,
    int ncols,
    int nrows,
    double y_max_percent) {
    if (results.empty()) {
        return;
    }

    const int width = 650 * ncols;
    const int height = 520 * nrows + 140;
    TCanvas* canvas = new TCanvas(("c_ptp_scatter_err_fit_" + file_tag).c_str(),
                                  ("c_ptp_scatter_err_fit_" + file_tag).c_str(),
                                  width,
                                  height);

    TPad* pTop = new TPad(("pTop_err_fit_" + file_tag).c_str(),
                          ("pTop_err_fit_" + file_tag).c_str(),
                          0.0,
                          0.89,
                          1.0,
                          1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad* pGrid = new TPad(("pGrid_err_fit_" + file_tag).c_str(),
                           ("pGrid_err_fit_" + file_tag).c_str(),
                           0.0,
                           0.00,
                           1.0,
                           0.89);
    pGrid->SetFillStyle(0);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(0.065);
    head.DrawLatex(0.5, 0.64, title_tag.c_str());

    TLatex sub;
    sub.SetNDC();
    sub.SetTextAlign(22);
    sub.SetTextFont(42);
    sub.SetTextSize(0.038);
    sub.DrawLatex(0.5, 0.30,
                  "Error bars show 100 #times #delta#bar{#sigma}_{scaled}/|#bar{#sigma}_{scaled}|; fits are p0 through p5");

    std::vector<std::unique_ptr<TGraphErrors> > graphs;
    std::vector<std::unique_ptr<TLegend> > legends;
    std::vector<std::unique_ptr<TF1> > fit_functions;

    for (int ivar = 0; ivar < (int)vars.size(); ++ivar) {
        pGrid->cd(ivar + 1);

        gPad->SetGrid(1, 1);
        gPad->SetTopMargin(0.14);
        gPad->SetBottomMargin(0.18);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.07);
        gPad->SetTicks(1, 1);

        std::vector<double> xs;
        std::vector<double> ys;
        std::vector<double> eys;

        for (const auto& r : results) {
            const double x = value_for_scatter_variable(r, vars[(size_t)ivar].key);
            const double y = r.f_ptp_percent;
            double ey = r.f_ptp_stat_percent;

            if (!std::isfinite(x) || !std::isfinite(y)) {
                continue;
            }

            if (!std::isfinite(ey) || ey < 0.0) {
                ey = 0.0;
            }

            xs.push_back(x);
            ys.push_back(y);
            eys.push_back(ey);
        }

        const std::vector<double> xr = choose_scatter_x_range(xs);
        TH1* frame = gPad->DrawFrame(xr[0], 0.0, xr[1], y_max_percent);
        frame->GetXaxis()->SetTitle(vars[(size_t)ivar].title.c_str());
        frame->GetYaxis()->SetTitle("Combination point-to-point sys (%)");
        frame->GetXaxis()->CenterTitle();
        frame->GetYaxis()->CenterTitle();
        frame->GetXaxis()->SetNdivisions(505);
        frame->GetXaxis()->SetTitleSize(0.060);
        frame->GetYaxis()->SetTitleSize(0.060);
        frame->GetXaxis()->SetLabelSize(0.050);
        frame->GetYaxis()->SetLabelSize(0.050);
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitleOffset(1.25);

        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(42);
        lab.SetTextSize(0.060);
        lab.SetTextAlign(22);
        lab.DrawLatex(0.50, 0.93, vars[(size_t)ivar].title.c_str());

        std::unique_ptr<TLegend> legend(new TLegend(0.50, 0.50, 0.93, 0.89));
        legend->SetBorderSize(1);
        legend->SetFillColor(kWhite);
        legend->SetFillStyle(1001);
        legend->SetTextFont(42);
        legend->SetTextSize(0.032);

        TGraphErrors* graph = make_scatter_graph_errors(xs, ys, eys);

        if (graph) {
            graph->Draw("P SAME");
            legend->AddEntry(graph, "points", "lep");

            fit_and_draw_polynomials(graph,
                                     file_tag,
                                     vars[(size_t)ivar].key,
                                     xr[0],
                                     xr[1],
                                     legend.get(),
                                     fit_functions);

            graphs.emplace_back(graph);
        }

        legend->Draw();
        legends.push_back(std::move(legend));
    }

    canvas->cd();
    canvas->Modified();
    canvas->Update();

    canvas->SaveAs((out_dir / ("combination_point_to_point_sys_" +
                               file_tag +
                               "_with_stat_errors_and_fits.pdf")).string().c_str());

    delete canvas;
}

static void make_global_scatter_plots(const std::vector<BinResult>& results,
                                      const fs::path& out_dir) {
    const double y_max_percent = choose_global_ymax_percent(results);

    draw_scatter_canvas(results,
                        physics_scatter_variables(),
                        out_dir,
                        "physics_kinematics",
                        "Point-to-point run-period combination systematic vs physics kinematics",
                        2,
                        2,
                        y_max_percent);

    draw_scatter_canvas(results,
                        angle_scatter_variables(),
                        out_dir,
                        "polar_angles",
                        "Point-to-point run-period combination systematic vs polar angles",
                        3,
                        1,
                        y_max_percent);

    draw_scatter_canvas_with_errors_and_fits(
        results,
        physics_scatter_variables(),
        out_dir,
        "physics_kinematics",
        "Point-to-point run-period combination systematic vs physics kinematics",
        2,
        2,
        y_max_percent);

    draw_scatter_canvas_with_errors_and_fits(
        results,
        angle_scatter_variables(),
        out_dir,
        "polar_angles",
        "Point-to-point run-period combination systematic vs polar angles",
        3,
        1,
        y_max_percent);
}

static void make_ptp_plots(const std::vector<BinResult>& results,
                           const fs::path& out_dir) {
    const std::map<Range, GroupByXB> grouped = group_results_by_xb(results);
    const double y_max_percent = choose_global_ymax_percent(results);

    int xb_counter = 0;

    for (const auto& kv : grouped) {
        draw_one_xb_canvas(kv.first,
                           kv.second,
                           xb_counter,
                           y_max_percent,
                           out_dir);
        ++xb_counter;
    }
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

static void print_global_summary(const std::vector<BinResult>& results) {
    if (results.empty()) {
        std::cout << "[combination-ptp] No valid bins.\n";
        return;
    }

    std::vector<double> fvals;
    std::vector<double> fstat_vals;
    std::vector<double> chi2_before;
    std::vector<double> chi2_after;

    for (const auto& r : results) {
        if (std::isfinite(r.f_ptp_percent)) {
            fvals.push_back(r.f_ptp_percent);
        }

        if (std::isfinite(r.f_ptp_stat_percent)) {
            fstat_vals.push_back(r.f_ptp_stat_percent);
        }

        if (std::isfinite(r.chi2_ndf_before)) {
            chi2_before.push_back(r.chi2_ndf_before);
        }

        if (std::isfinite(r.chi2_ndf_after)) {
            chi2_after.push_back(r.chi2_ndf_after);
        }
    }

    std::cout << "[combination-ptp] Summary for 10.6 GeV unpol residual point-to-point systematic\n";
    std::cout << "  valid bins                  = " << results.size() << "\n";
    std::cout << "  mean chi2/ndf before        = " << std::setprecision(8) << mean_vec(chi2_before) << "\n";
    std::cout << "  mean chi2/ndf after         = " << std::setprecision(8) << mean_vec(chi2_after) << "\n";
    std::cout << "  mean f_ptp percent          = " << std::setprecision(8) << mean_vec(fvals) << "\n";
    std::cout << "  mean f_ptp stat percent     = " << std::setprecision(8) << mean_vec(fstat_vals) << "\n";
    std::cout << "  median f_ptp percent        = " << std::setprecision(8) << percentile(fvals, 0.50) << "\n";
    std::cout << "  p68 f_ptp percent           = " << std::setprecision(8) << percentile(fvals, 0.68) << "\n";
    std::cout << "  p90 f_ptp percent           = " << std::setprecision(8) << percentile(fvals, 0.90) << "\n";
    std::cout << "  p95 f_ptp percent           = " << std::setprecision(8) << percentile(fvals, 0.95) << "\n";
    std::cout << "  max f_ptp percent           = " << std::setprecision(8)
              << (fvals.empty() ? std::numeric_limits<double>::quiet_NaN()
                                : *std::max_element(fvals.begin(), fvals.end()))
              << "\n";
}

} // namespace

bool combination_point_to_point_systematics(const std::string& csv_path,
                                            const std::string& output_root_dir) {
    try {
        const fs::path out_dir =
            fs::path(output_root_dir) / "combination_point_to_point_systematics";

        fs::create_directories(out_dir);

        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTitleFont(42, "");
        gStyle->SetTextFont(42);
        gStyle->SetEndErrorSize(0);

        const CsvTable table = read_csv_or_throw(csv_path);
        validate_schema(table);

        std::cout << "[combination-ptp] CSV rows loaded: "
                  << table.rows.size() << "\n";
        std::cout << "[combination-ptp] Output directory: "
                  << out_dir.string() << "\n";

        const std::vector<PeriodScale> scales =
            compute_10p6_unpol_scale_factors(table);

        const std::map<std::string, double> scale_by_period =
            scale_map_from_vector(scales);

        const std::vector<BinResult> results =
            evaluate_all_bins(table, scale_by_period);

        const fs::path scale_csv =
            out_dir / "period_scale_factors_used.csv";

        const fs::path bin_csv =
            out_dir / "combination_point_to_point_systematics_by_bin.csv";

        write_scale_summary_csv(scale_csv.string(), scales);
        write_bin_summary_csv(bin_csv.string(), results);

        make_ptp_plots(results, out_dir);
        make_global_scatter_plots(results, out_dir);
        print_global_summary(results);

        std::cout << "[combination-ptp] Wrote scale summary: "
                  << scale_csv.string() << "\n";
        std::cout << "[combination-ptp] Wrote bin summary: "
                  << bin_csv.string() << "\n";

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[combination-ptp] FATAL: " << e.what() << "\n";
        return false;
    }
}