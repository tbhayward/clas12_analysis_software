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

static constexpr int kMaxEThetaScaleFitOrder = 5;
static constexpr double kMinEThetaScaleFitImprovement = 0.01;

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

struct RatioPoint {
    std::string period;
    double e_theta = std::numeric_limits<double>::quiet_NaN();
    double ratio = std::numeric_limits<double>::quiet_NaN();
    double ratio_stat = std::numeric_limits<double>::quiet_NaN();
};

struct EThetaScaleFit {
    std::string period;
    int order = -1;
    int n = 0;
    int ndf = 0;
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    double chi2_ndf = std::numeric_limits<double>::quiet_NaN();
    double x_min = std::numeric_limits<double>::quiet_NaN();
    double x_max = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> params;
    std::vector<double> errors;
};

struct ScaleModel {
    std::string tag;
    std::string label;
    bool use_e_theta_fit = false;
    std::vector<PeriodScale> constant_scales;
    std::map<std::string, double> constant_scale_by_period;
    std::map<std::string, EThetaScaleFit> etheta_fit_by_period;
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

    std::cout << "[combination-ptp] 10.6 GeV unpol constant scale factors:\n";

    for (const auto& input : inputs) {
        const PeriodRatioAccumulator& acc = acc_by_period[input.period];

        if (acc.sum_w <= 0.0 || acc.n <= 0) {
            throw std::runtime_error("Could not compute constant scale factor for period: " + input.period);
        }

        PeriodScale scale;
        scale.period = input.period;
        scale.scale = acc.sum_wr / acc.sum_w;
        scale.scale_stat = 1.0 / std::sqrt(acc.sum_w);
        scale.n = acc.n;

        if (!std::isfinite(scale.scale) || std::abs(scale.scale) <= 0.0) {
            throw std::runtime_error("Invalid constant scale factor for period: " + input.period);
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

static std::vector<RatioPoint> build_e_theta_ratio_points(const CsvTable& table) {
    const std::vector<PeriodInput> inputs = base_inputs();
    std::vector<RatioPoint> points;

    for (const auto& row : table.rows) {
        const double e_theta = get_double_or_nan(table, row, "e_theta, 10.6 GeV");

        if (!std::isfinite(e_theta)) {
            continue;
        }

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

            if (!v.ok || !std::isfinite(v.value) || !std::isfinite(v.stat) || v.stat <= 0.0) {
                continue;
            }

            RatioPoint p;
            p.period = inputs[i].period;
            p.e_theta = e_theta;
            p.ratio = v.value / mean;
            p.ratio_stat = std::abs(v.stat / mean);

            if (!std::isfinite(p.ratio) ||
                !std::isfinite(p.ratio_stat) ||
                p.ratio_stat <= 0.0) {
                continue;
            }

            points.push_back(p);
        }
    }

    return points;
}

static std::vector<RatioPoint> filter_ratio_points_for_period(const std::vector<RatioPoint>& points,
                                                              const std::string& period) {
    std::vector<RatioPoint> out;

    for (const auto& p : points) {
        if (p.period != period) {
            continue;
        }

        if (!std::isfinite(p.e_theta) ||
            !std::isfinite(p.ratio) ||
            !std::isfinite(p.ratio_stat) ||
            p.ratio_stat <= 0.0) {
            continue;
        }

        out.push_back(p);
    }

    std::sort(out.begin(), out.end(),
              [](const RatioPoint& a, const RatioPoint& b) {
                  return a.e_theta < b.e_theta;
              });

    return out;
}

static EThetaScaleFit weighted_polynomial_fit_ratio_points(const std::string& period,
                                                           const std::vector<RatioPoint>& points,
                                                           int order,
                                                           double previous_chi2_ndf) {
    EThetaScaleFit out;
    out.period = period;
    out.order = order;
    out.n = (int)points.size();

    const int npar = order + 1;
    out.ndf = out.n - npar;

    if (out.n < npar || out.ndf <= 0) {
        return out;
    }

    std::vector<std::vector<double> > normal((size_t)npar, std::vector<double>((size_t)npar, 0.0));
    std::vector<double> rhs((size_t)npar, 0.0);

    for (const auto& p : points) {
        if (!std::isfinite(p.e_theta) ||
            !std::isfinite(p.ratio) ||
            !std::isfinite(p.ratio_stat) ||
            p.ratio_stat <= 0.0) {
            continue;
        }

        const double w = 1.0 / (p.ratio_stat * p.ratio_stat);
        std::vector<double> xp((size_t)npar, 1.0);

        for (int i = 1; i < npar; ++i) {
            xp[(size_t)i] = xp[(size_t)(i - 1)] * p.e_theta;
        }

        for (int i = 0; i < npar; ++i) {
            rhs[(size_t)i] += w * xp[(size_t)i] * p.ratio;

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
            } else {
                out.errors[(size_t)i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    } else {
        out.errors.assign((size_t)npar, std::numeric_limits<double>::quiet_NaN());
    }

    double chi2 = 0.0;
    int n_used = 0;
    double x_min = std::numeric_limits<double>::infinity();
    double x_max = -std::numeric_limits<double>::infinity();

    for (const auto& p : points) {
        if (!std::isfinite(p.e_theta) ||
            !std::isfinite(p.ratio) ||
            !std::isfinite(p.ratio_stat) ||
            p.ratio_stat <= 0.0) {
            continue;
        }

        x_min = std::min(x_min, p.e_theta);
        x_max = std::max(x_max, p.e_theta);

        const double residual = p.ratio - eval_poly(out.params, p.e_theta);
        chi2 += residual * residual / (p.ratio_stat * p.ratio_stat);
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

    bool accepted = false;

    if (order == 0) {
        accepted = std::isfinite(out.chi2_ndf);
    } else if (std::isfinite(previous_chi2_ndf) && std::isfinite(out.chi2_ndf)) {
        const double improvement = previous_chi2_ndf - out.chi2_ndf;
        accepted = improvement >= kMinEThetaScaleFitImprovement;
    }

    if (!accepted) {
        out.params.clear();
        out.errors.clear();
    }

    return out;
}

static EThetaScaleFit choose_e_theta_scale_fit(const std::string& period,
                                               const std::vector<RatioPoint>& points) {
    double previous_chi2_ndf = std::numeric_limits<double>::quiet_NaN();
    EThetaScaleFit last_accepted;

    for (int order = 0; order <= kMaxEThetaScaleFitOrder; ++order) {
        EThetaScaleFit fit =
            weighted_polynomial_fit_ratio_points(period,
                                                 points,
                                                 order,
                                                 previous_chi2_ndf);

        const bool attempted_ok =
            !fit.params.empty() &&
            std::isfinite(fit.chi2_ndf) &&
            fit.ndf > 0;

        bool accepted = false;

        if (attempted_ok) {
            if (order == 0) {
                accepted = true;
            } else {
                const double improvement = previous_chi2_ndf - fit.chi2_ndf;
                accepted = std::isfinite(improvement) &&
                           improvement >= kMinEThetaScaleFitImprovement;
            }
        }

        std::cout << "[combination-ptp][e-theta-scale-fit] "
                  << period
                  << " p" << order
                  << " n=" << fit.n
                  << " chi2/ndf=" << std::setprecision(8) << fit.chi2_ndf
                  << (accepted ? " accepted" : " rejected")
                  << "\n";

        if (!accepted) {
            break;
        }

        previous_chi2_ndf = fit.chi2_ndf;
        last_accepted = fit;
    }

    if (last_accepted.params.empty()) {
        throw std::runtime_error("Could not determine e_theta scale fit for period: " + period);
    }

    return last_accepted;
}

static std::map<std::string, EThetaScaleFit>
compute_e_theta_scale_fits(const CsvTable& table) {
    const std::vector<RatioPoint> all_points = build_e_theta_ratio_points(table);

    std::map<std::string, EThetaScaleFit> out;

    std::cout << "[combination-ptp] 10.6 GeV unpol e_theta-dependent scale factors:\n";

    for (const auto& period : kBasePeriods) {
        const std::vector<RatioPoint> points =
            filter_ratio_points_for_period(all_points, period);

        if (points.empty()) {
            throw std::runtime_error("No e_theta ratio points available for period: " + period);
        }

        EThetaScaleFit fit = choose_e_theta_scale_fit(period, points);

        std::cout << "  " << period
                  << " selected p" << fit.order
                  << " chi2/ndf=" << std::setprecision(10) << fit.chi2_ndf
                  << " support=(" << fit.x_min << ", " << fit.x_max << ")";

        for (int ip = 0; ip < (int)fit.params.size(); ++ip) {
            std::cout << " p" << ip << "=" << fit.params[(size_t)ip];

            if (ip < (int)fit.errors.size()) {
                std::cout << "+/-" << fit.errors[(size_t)ip];
            }
        }

        std::cout << "\n";

        out[period] = fit;
    }

    return out;
}

static double e_theta_fit_scale_value(const EThetaScaleFit& fit,
                                      double e_theta) {
    if (!std::isfinite(e_theta)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (!std::isfinite(fit.x_min) ||
        !std::isfinite(fit.x_max) ||
        e_theta < fit.x_min ||
        e_theta > fit.x_max) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (fit.params.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double scale = eval_poly(fit.params, e_theta);

    if (!std::isfinite(scale) || std::abs(scale) <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return scale;
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

static void fill_common_bin_metadata(const CsvTable& table,
                                     const std::vector<std::string>& row,
                                     size_t row_index,
                                     BinResult& result) {
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
}

static double scale_value_for_period_and_row(const ScaleModel& model,
                                             const std::string& period,
                                             double e_theta) {
    if (!model.use_e_theta_fit) {
        const auto it = model.constant_scale_by_period.find(period);

        if (it == model.constant_scale_by_period.end()) {
            throw std::runtime_error("Missing constant scale factor for period: " + period);
        }

        return it->second;
    }

    const auto it = model.etheta_fit_by_period.find(period);

    if (it == model.etheta_fit_by_period.end()) {
        throw std::runtime_error("Missing e_theta scale fit for period: " + period);
    }

    return e_theta_fit_scale_value(it->second, e_theta);
}

static BinResult evaluate_row(const CsvTable& table,
                              const std::vector<std::string>& row,
                              size_t row_index,
                              const ScaleModel& model) {
    BinResult result;
    fill_common_bin_metadata(table, row, row_index, result);

    std::vector<TupleValue> scaled_values;

    for (const auto& input : base_inputs()) {
        const TupleValue raw = get_tuple(table, row, input.column);
        const double scale =
            scale_value_for_period_and_row(model, input.period, result.e_theta);

        if (!std::isfinite(scale) || std::abs(scale) <= 0.0) {
            continue;
        }

        const TupleValue scaled = apply_scale(raw, scale);

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
                  const ScaleModel& model) {
    std::vector<BinResult> results;
    results.reserve(table.rows.size());

    for (size_t i = 0; i < table.rows.size(); ++i) {
        BinResult result = evaluate_row(table, table.rows[i], i + 1, model);

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

static void write_constant_scale_summary_csv(const std::string& path,
                                             const std::vector<PeriodScale>& scales) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open constant scale summary CSV for writing: " + path);
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
        throw std::runtime_error("Failed while writing constant scale summary CSV: " + path);
    }
}

static void write_e_theta_fit_summary_csv(const std::string& path,
                                          const std::map<std::string, EThetaScaleFit>& fits) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open e_theta scale fit summary CSV for writing: " + path);
    }

    fout << "period,order,n,ndf,chi2,chi2ndf,e_theta_min,e_theta_max";

    for (int ip = 0; ip <= kMaxEThetaScaleFitOrder; ++ip) {
        fout << ",p" << ip << ",p" << ip << " err";
    }

    fout << "\n";

    for (const auto& kv : fits) {
        const EThetaScaleFit& fit = kv.second;

        fout << csv_escape_field(fit.period) << ","
             << fit.order << ","
             << fit.n << ","
             << fit.ndf << ","
             << format_double(fit.chi2) << ","
             << format_double(fit.chi2_ndf) << ","
             << format_double(fit.x_min) << ","
             << format_double(fit.x_max);

        for (int ip = 0; ip <= kMaxEThetaScaleFitOrder; ++ip) {
            if (ip < (int)fit.params.size()) {
                fout << "," << format_double(fit.params[(size_t)ip]);

                if (ip < (int)fit.errors.size()) {
                    fout << "," << format_double(fit.errors[(size_t)ip]);
                } else {
                    fout << ",";
                }
            } else {
                fout << ",,";
            }
        }

        fout << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing e_theta scale fit summary CSV: " + path);
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
                               const fs::path& out_dir,
                               const ScaleModel& model) {
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

    TCanvas* canvas = new TCanvas(("c_combination_ptp_" + model.tag).c_str(),
                                  ("c_combination_ptp_" + model.tag).c_str(),
                                  width,
                                  height);

    TPad* pTop = new TPad(("pTopCombinationPtp_" + model.tag).c_str(),
                          ("pTopCombinationPtp_" + model.tag).c_str(),
                          0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad* pGrid = new TPad(("pGridCombinationPtp_" + model.tag).c_str(),
                           ("pGridCombinationPtp_" + model.tag).c_str(),
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

    std::string subtitle;
    if (model.use_e_theta_fit) {
        subtitle =
            "e_{#theta}-dependent period scale offsets removed first; plotted value is the residual fractional uncertainty needed to make #chi^{2}/ndf = 1";
    } else {
        subtitle =
            "Constant global period scale offsets removed first; plotted value is the residual fractional uncertainty needed to make #chi^{2}/ndf = 1";
    }

    sub.DrawLatex(0.5, 0.50, subtitle.c_str());

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
    fname << "combination_point_to_point_sys_10p6_GeV_unpol_"
          << model.tag
          << "_xB_"
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
                                      const fs::path& out_dir,
                                      const ScaleModel& model) {
    const double y_max_percent = choose_global_ymax_percent(results);

    draw_scatter_canvas(results,
                        physics_scatter_variables(),
                        out_dir,
                        model.tag + "_physics_kinematics",
                        "Point-to-point run-period combination systematic vs physics kinematics, " + model.label,
                        2,
                        2,
                        y_max_percent);

    draw_scatter_canvas(results,
                        angle_scatter_variables(),
                        out_dir,
                        model.tag + "_polar_angles",
                        "Point-to-point run-period combination systematic vs polar angles, " + model.label,
                        3,
                        1,
                        y_max_percent);

    draw_scatter_canvas_with_errors_and_fits(
        results,
        physics_scatter_variables(),
        out_dir,
        model.tag + "_physics_kinematics",
        "Point-to-point run-period combination systematic vs physics kinematics, " + model.label,
        2,
        2,
        y_max_percent);

    draw_scatter_canvas_with_errors_and_fits(
        results,
        angle_scatter_variables(),
        out_dir,
        model.tag + "_polar_angles",
        "Point-to-point run-period combination systematic vs polar angles, " + model.label,
        3,
        1,
        y_max_percent);
}

static void make_ptp_plots(const std::vector<BinResult>& results,
                           const fs::path& out_dir,
                           const ScaleModel& model) {
    const std::map<Range, GroupByXB> grouped = group_results_by_xb(results);
    const double y_max_percent = choose_global_ymax_percent(results);

    int xb_counter = 0;

    for (const auto& kv : grouped) {
        draw_one_xb_canvas(kv.first,
                           kv.second,
                           xb_counter,
                           y_max_percent,
                           out_dir,
                           model);
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

static void print_global_summary(const std::vector<BinResult>& results,
                                 const ScaleModel& model) {
    if (results.empty()) {
        std::cout << "[combination-ptp] No valid bins for " << model.tag << ".\n";
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

    std::cout << "[combination-ptp] Summary for 10.6 GeV unpol residual point-to-point systematic, "
              << model.label << "\n";
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

static void write_global_summary_csv_line(std::ofstream& fout,
                                          const ScaleModel& model,
                                          const std::vector<BinResult>& results) {
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

    fout << csv_escape_field(model.tag) << ","
         << csv_escape_field(model.label) << ","
         << results.size() << ","
         << format_double(mean_vec(chi2_before)) << ","
         << format_double(mean_vec(chi2_after)) << ","
         << format_double(mean_vec(fvals)) << ","
         << format_double(mean_vec(fstat_vals)) << ","
         << format_double(percentile(fvals, 0.50)) << ","
         << format_double(percentile(fvals, 0.68)) << ","
         << format_double(percentile(fvals, 0.90)) << ","
         << format_double(percentile(fvals, 0.95)) << ","
         << format_double(fvals.empty() ? std::numeric_limits<double>::quiet_NaN()
                                        : *std::max_element(fvals.begin(), fvals.end()))
         << "\n";
}

static void write_global_comparison_summary_csv(const fs::path& path,
                                                const std::vector<std::pair<ScaleModel, std::vector<BinResult> > >& all_results) {
    std::ofstream fout(path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open global comparison summary CSV for writing: " + path.string());
    }

    fout << "model tag,model label,valid bins,mean chi2ndf before,mean chi2ndf after,"
         << "mean f_ptp percent,mean f_ptp stat percent,median f_ptp percent,"
         << "p68 f_ptp percent,p90 f_ptp percent,p95 f_ptp percent,max f_ptp percent\n";

    for (const auto& item : all_results) {
        write_global_summary_csv_line(fout, item.first, item.second);
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing global comparison summary CSV: " + path.string());
    }
}

static void draw_e_theta_scale_fit_canvas(const std::vector<RatioPoint>& all_points,
                                          const std::map<std::string, EThetaScaleFit>& fits,
                                          const fs::path& out_dir) {
    TCanvas* canvas = new TCanvas("c_e_theta_scale_fits",
                                  "c_e_theta_scale_fits",
                                  1700,
                                  1300);

    canvas->Divide(2, 2, 0.001, 0.001);

    std::vector<std::unique_ptr<TGraphErrors> > graphs;
    std::vector<std::unique_ptr<TF1> > funcs;
    std::vector<std::unique_ptr<TLegend> > legends;

    for (int iper = 0; iper < (int)kBasePeriods.size(); ++iper) {
        const std::string& period = kBasePeriods[(size_t)iper];

        canvas->cd(iper + 1);
        gPad->SetGrid(1, 1);
        gPad->SetTopMargin(0.15);
        gPad->SetBottomMargin(0.16);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.07);
        gPad->SetTicks(1, 1);

        const std::vector<RatioPoint> points =
            filter_ratio_points_for_period(all_points, period);

        std::vector<double> xs;
        for (const auto& p : points) {
            xs.push_back(p.e_theta);
        }

        const std::vector<double> xr = choose_scatter_x_range(xs);
        TH1* frame = gPad->DrawFrame(xr[0], 0.0, xr[1], 2.0);
        frame->GetXaxis()->SetTitle("#theta_{e} (deg)");
        frame->GetYaxis()->SetTitle("#sigma_{i}/#bar{#sigma}");
        frame->GetXaxis()->CenterTitle();
        frame->GetYaxis()->CenterTitle();
        frame->GetXaxis()->SetNdivisions(505);
        frame->GetYaxis()->SetNdivisions(505);
        frame->GetXaxis()->SetTitleSize(0.060);
        frame->GetYaxis()->SetTitleSize(0.060);
        frame->GetXaxis()->SetLabelSize(0.050);
        frame->GetYaxis()->SetLabelSize(0.050);
        frame->GetXaxis()->SetTitleOffset(1.10);
        frame->GetYaxis()->SetTitleOffset(1.15);

        TLine unity(xr[0], 1.0, xr[1], 1.0);
        unity.SetLineColor(kRed + 1);
        unity.SetLineStyle(2);
        unity.SetLineWidth(1);
        unity.Draw("SAME");

        std::unique_ptr<TGraphErrors> graph(new TGraphErrors((int)points.size()));

        for (int i = 0; i < (int)points.size(); ++i) {
            graph->SetPoint(i, points[(size_t)i].e_theta, points[(size_t)i].ratio);
            graph->SetPointError(i, 0.0, points[(size_t)i].ratio_stat);
        }

        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.35);
        graph->SetMarkerColor(kBlack);
        graph->SetLineColor(kBlack);
        graph->SetLineWidth(1);
        graph->Draw("P SAME");

        std::unique_ptr<TLegend> legend(new TLegend(0.48, 0.70, 0.93, 0.88));
        legend->SetBorderSize(1);
        legend->SetFillColor(kWhite);
        legend->SetFillStyle(1001);
        legend->SetTextFont(42);
        legend->SetTextSize(0.032);
        legend->AddEntry(graph.get(), period.c_str(), "pe");

        const auto it_fit = fits.find(period);
        if (it_fit == fits.end()) {
            throw std::runtime_error("Missing e_theta fit while drawing fit canvas for period: " + period);
        }

        const EThetaScaleFit& fit = it_fit->second;
        std::unique_ptr<TF1> func(new TF1(("f_e_theta_scale_" + sanitize_for_path(period)).c_str(),
                                          poly_formula(fit.order).c_str(),
                                          fit.x_min,
                                          fit.x_max));

        for (int ip = 0; ip < (int)fit.params.size(); ++ip) {
            func->SetParameter(ip, fit.params[(size_t)ip]);
        }

        func->SetLineColor(kBlue + 1);
        func->SetLineWidth(2);
        func->SetLineStyle(1);
        func->Draw("SAME");

        std::ostringstream fit_label;
        fit_label << "p" << fit.order
                  << ", #chi^{2}/ndf="
                  << std::fixed << std::setprecision(2)
                  << fit.chi2_ndf;
        legend->AddEntry(func.get(), fit_label.str().c_str(), "l");
        legend->Draw();

        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(42);
        lab.SetTextSize(0.055);
        lab.SetTextAlign(22);
        lab.DrawLatex(0.50, 0.94, period.c_str());

        graphs.push_back(std::move(graph));
        funcs.push_back(std::move(func));
        legends.push_back(std::move(legend));
    }

    canvas->cd();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.022);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.992, "10.6 GeV unpol e_{#theta}-dependent period scale fits");

    canvas->Modified();
    canvas->Update();

    const fs::path out_path = out_dir / "e_theta_dependent_period_scale_fits.pdf";
    canvas->SaveAs(out_path.string().c_str());

    delete canvas;
}

static void draw_model_comparison_scatter(const std::vector<BinResult>& constant_results,
                                          const std::vector<BinResult>& etheta_results,
                                          const fs::path& out_dir) {
    std::map<size_t, BinResult> constant_by_row;
    std::map<size_t, BinResult> etheta_by_row;

    for (const auto& r : constant_results) {
        constant_by_row[r.row_index] = r;
    }

    for (const auto& r : etheta_results) {
        etheta_by_row[r.row_index] = r;
    }

    std::vector<double> x_const;
    std::vector<double> y_etheta;

    for (const auto& kv : constant_by_row) {
        const auto it = etheta_by_row.find(kv.first);
        if (it == etheta_by_row.end()) {
            continue;
        }

        if (!std::isfinite(kv.second.f_ptp_percent) ||
            !std::isfinite(it->second.f_ptp_percent)) {
            continue;
        }

        x_const.push_back(kv.second.f_ptp_percent);
        y_etheta.push_back(it->second.f_ptp_percent);
    }

    if (x_const.empty()) {
        return;
    }

    double max_value = 0.0;

    for (size_t i = 0; i < x_const.size(); ++i) {
        max_value = std::max(max_value, x_const[i]);
        max_value = std::max(max_value, y_etheta[i]);
    }

    max_value = std::max(5.0, 1.15 * max_value);

    TCanvas* canvas = new TCanvas("c_ptp_constant_vs_etheta",
                                  "c_ptp_constant_vs_etheta",
                                  900,
                                  800);

    gPad->SetGrid(1, 1);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.06);
    gPad->SetTicks(1, 1);

    TH1* frame = gPad->DrawFrame(0.0, 0.0, max_value, max_value);
    frame->GetXaxis()->SetTitle("P2P sys with constant scale removal (%)");
    frame->GetYaxis()->SetTitle("P2P sys with e_{#theta}-dependent scale removal (%)");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetTitleSize(0.050);
    frame->GetYaxis()->SetTitleSize(0.050);
    frame->GetXaxis()->SetLabelSize(0.043);
    frame->GetYaxis()->SetLabelSize(0.043);
    frame->GetXaxis()->SetTitleOffset(1.20);
    frame->GetYaxis()->SetTitleOffset(1.35);

    TLine unity(0.0, 0.0, max_value, max_value);
    unity.SetLineColor(kRed + 1);
    unity.SetLineStyle(2);
    unity.SetLineWidth(2);
    unity.Draw("SAME");

    TGraph* graph = new TGraph((int)x_const.size());
    for (int i = 0; i < (int)x_const.size(); ++i) {
        graph->SetPoint(i, x_const[(size_t)i], y_etheta[(size_t)i]);
    }

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.55);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlack);
    graph->Draw("P SAME");

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.038);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.955, "Point-to-point systematic comparison by scale-removal method");

    TLegend legend(0.20, 0.76, 0.58, 0.88);
    legend.SetBorderSize(1);
    legend.SetFillColor(kWhite);
    legend.SetFillStyle(1001);
    legend.SetTextFont(42);
    legend.SetTextSize(0.030);
    legend.AddEntry(graph, "matched bins", "p");
    legend.AddEntry(&unity, "same value", "l");
    legend.Draw();

    canvas->Modified();
    canvas->Update();

    const fs::path out_path = out_dir / "comparison_constant_vs_e_theta_ptp_sys.pdf";
    canvas->SaveAs(out_path.string().c_str());

    delete graph;
    delete canvas;
}

static void run_one_model_outputs(const CsvTable& table,
                                  const ScaleModel& model,
                                  const fs::path& model_dir,
                                  std::vector<BinResult>& results_out) {
    fs::create_directories(model_dir);

    const std::vector<BinResult> results =
        evaluate_all_bins(table, model);

    results_out = results;

    const fs::path bin_csv =
        model_dir / "combination_point_to_point_systematics_by_bin.csv";

    write_bin_summary_csv(bin_csv.string(), results);

    make_ptp_plots(results, model_dir, model);
    make_global_scatter_plots(results, model_dir, model);
    print_global_summary(results, model);

    std::cout << "[combination-ptp] Wrote bin summary for " << model.tag
              << ": " << bin_csv.string() << "\n";
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

        const std::vector<PeriodScale> constant_scales =
            compute_10p6_unpol_scale_factors(table);

        const std::map<std::string, double> constant_scale_by_period =
            scale_map_from_vector(constant_scales);

        const std::map<std::string, EThetaScaleFit> etheta_fit_by_period =
            compute_e_theta_scale_fits(table);

        const std::vector<RatioPoint> etheta_ratio_points =
            build_e_theta_ratio_points(table);

        const fs::path constant_scale_csv =
            out_dir / "constant_period_scale_factors_used.csv";

        const fs::path etheta_fit_csv =
            out_dir / "e_theta_dependent_period_scale_fits.csv";

        write_constant_scale_summary_csv(constant_scale_csv.string(), constant_scales);
        write_e_theta_fit_summary_csv(etheta_fit_csv.string(), etheta_fit_by_period);
        draw_e_theta_scale_fit_canvas(etheta_ratio_points, etheta_fit_by_period, out_dir);

        ScaleModel constant_model;
        constant_model.tag = "constant_scale";
        constant_model.label = "constant period-scale removal";
        constant_model.use_e_theta_fit = false;
        constant_model.constant_scales = constant_scales;
        constant_model.constant_scale_by_period = constant_scale_by_period;

        ScaleModel etheta_model;
        etheta_model.tag = "e_theta_scale";
        etheta_model.label = "e_{#theta}-dependent period-scale removal";
        etheta_model.use_e_theta_fit = true;
        etheta_model.constant_scales = constant_scales;
        etheta_model.constant_scale_by_period = constant_scale_by_period;
        etheta_model.etheta_fit_by_period = etheta_fit_by_period;

        std::vector<BinResult> constant_results;
        std::vector<BinResult> etheta_results;

        run_one_model_outputs(table,
                              constant_model,
                              out_dir / "constant_scale",
                              constant_results);

        run_one_model_outputs(table,
                              etheta_model,
                              out_dir / "e_theta_scale",
                              etheta_results);

        draw_model_comparison_scatter(constant_results,
                                      etheta_results,
                                      out_dir);

        std::vector<std::pair<ScaleModel, std::vector<BinResult> > > all_results;
        all_results.push_back(std::make_pair(constant_model, constant_results));
        all_results.push_back(std::make_pair(etheta_model, etheta_results));

        const fs::path global_summary_csv =
            out_dir / "combination_point_to_point_systematics_global_summary.csv";

        write_global_comparison_summary_csv(global_summary_csv,
                                            all_results);

        std::cout << "[combination-ptp] Wrote constant scale summary: "
                  << constant_scale_csv.string() << "\n";
        std::cout << "[combination-ptp] Wrote e_theta scale fit summary: "
                  << etheta_fit_csv.string() << "\n";
        std::cout << "[combination-ptp] Wrote global comparison summary: "
                  << global_summary_csv.string() << "\n";

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[combination-ptp] FATAL: " << e.what() << "\n";
        return false;
    }
}