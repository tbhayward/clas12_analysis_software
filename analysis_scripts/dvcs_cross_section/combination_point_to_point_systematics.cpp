#include "combination_point_to_point_systematics.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TH1.h>

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
    bool ok = false;
    size_t row_index = 0;
    Range xb_range;
    Range q2_range;
    Range t_range;
    int xb_index = -1;
    double xBavg = std::numeric_limits<double>::quiet_NaN();
    double Q2avg = std::numeric_limits<double>::quiet_NaN();
    double t_abs_avg = std::numeric_limits<double>::quiet_NaN();
    double phiavg = std::numeric_limits<double>::quiet_NaN();
    double ref_value = std::numeric_limits<double>::quiet_NaN();
    double ref_stat = std::numeric_limits<double>::quiet_NaN();
    double chi2_ndf_before = std::numeric_limits<double>::quiet_NaN();
    double chi2_ndf_after = std::numeric_limits<double>::quiet_NaN();
    double f_ptp = std::numeric_limits<double>::quiet_NaN();
    double f_ptp_percent = std::numeric_limits<double>::quiet_NaN();
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

static const std::vector<std::string> kBasePeriods = {
    "Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"
};

static std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;
    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
    return s.substr(b, e - b);
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

static CsvTable read_csv_or_throw(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Could not open CSV: " + path);

    CsvTable table;
    std::string line;
    if (!std::getline(fin, line)) throw std::runtime_error("CSV is empty: " + path);

    table.header = split_csv_line(line);
    std::set<std::string> seen;
    std::set<std::string> dup;
    for (int i = 0; i < (int)table.header.size(); ++i) {
        const std::string name = trim(table.header[(size_t)i]);
        if (!seen.insert(name).second) dup.insert(name);
        table.index[name] = i;
    }
    if (!dup.empty()) {
        std::ostringstream msg;
        msg << "CSV header contains duplicate columns:";
        for (const auto& d : dup) msg << "\n  - " << d;
        throw std::runtime_error(msg.str());
    }

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = split_csv_line(line);
        if (row.size() < table.header.size()) row.resize(table.header.size());
        if (row.size() != table.header.size()) {
            std::ostringstream msg;
            msg << "CSV row has " << row.size() << " columns, but header has "
                << table.header.size() << " columns. File: " << path;
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
    for (const auto& col : required) if (table.index.find(col) == table.index.end()) missing.push_back(col);
    if (!missing.empty()) {
        std::ostringstream msg;
        msg << "Missing required columns for " << context << ":";
        for (const auto& col : missing) msg << "\n  - " << col;
        throw std::runtime_error(msg.str());
    }
}

static bool parse_double(const std::string& raw, double& out) {
    const std::string s = trim(raw);
    if (s.empty()) return false;
    char* end = nullptr;
    out = std::strtod(s.c_str(), &end);
    if (end == s.c_str()) return false;
    while (end && *end != '\0') {
        if (!std::isspace((unsigned char)(*end))) return false;
        ++end;
    }
    return std::isfinite(out);
}

static TupleValue parse_tuple_value(const std::string& raw) {
    TupleValue out;
    std::string s = trim(raw);
    if (s.empty()) return out;
    if (s.front() == '(' && s.back() == ')') s = trim(s.substr(1, s.size() - 2));
    const std::vector<std::string> fields = split_csv_line(s);
    if (fields.size() < 2) return out;
    double value = 0.0;
    double stat = 0.0;
    double sys = 0.0;
    if (!parse_double(fields[0], value)) return out;
    if (!parse_double(fields[1], stat)) return out;
    if (fields.size() > 2) parse_double(fields[2], sys);
    if (!std::isfinite(value) || !std::isfinite(stat) || stat <= 0.0) return out;
    out.ok = true;
    out.value = value;
    out.stat = stat;
    out.sys = sys;
    return out;
}

static int get_col(const CsvTable& table, const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) throw std::runtime_error("Missing required column: " + column);
    return it->second;
}

static int get_col_optional(const CsvTable& table, const std::string& column) {
    const auto it = table.index.find(column);
    if (it == table.index.end()) return -1;
    return it->second;
}

static double get_double_or_nan(const CsvTable& table,
                                const std::vector<std::string>& row,
                                const std::string& column) {
    const int col = get_col(table, column);
    double value = 0.0;
    if (!parse_double(row[(size_t)col], value)) return std::numeric_limits<double>::quiet_NaN();
    return value;
}

static TupleValue get_tuple(const CsvTable& table,
                            const std::vector<std::string>& row,
                            const std::string& column) {
    return parse_tuple_value(row[(size_t)get_col(table, column)]);
}

static std::string cross_section_column(const std::string& label,
                                        const std::string& helicity) {
    return "normed cross sections, ep->epg, exp, " + label + ", " + helicity;
}

static std::string reference_column() {
    return cross_section_column("10.6 GeV", "unpol");
}

static std::vector<PeriodInput> base_inputs() {
    std::vector<PeriodInput> inputs;
    for (const auto& period : kBasePeriods) {
        inputs.push_back({period, cross_section_column(period, "unpol")});
    }
    return inputs;
}

static void validate_schema(const CsvTable& table) {
    std::vector<std::string> required = {
        "xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max",
        "phimin", "phimax", "xBavg, 10.6 GeV", "Q2avg, 10.6 GeV",
        "t_abs_avg, 10.6 GeV", "phiavg, 10.6 GeV", reference_column()
    };
    for (const auto& input : base_inputs()) required.push_back(input.column);
    require_columns(table, required, "combination point-to-point systematics");
}

static bool weighted_mean(const std::vector<TupleValue>& values,
                          TupleValue& out) {
    double sum_w = 0.0;
    double sum_wx = 0.0;
    for (const auto& v : values) {
        if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) continue;
        const double w = 1.0 / (v.stat * v.stat);
        sum_w += w;
        sum_wx += w * v.value;
    }
    if (sum_w <= 0.0) return false;
    out.ok = true;
    out.value = sum_wx / sum_w;
    out.stat = 1.0 / std::sqrt(sum_w);
    return std::isfinite(out.value) && std::isfinite(out.stat) && out.stat > 0.0;
}

static std::vector<PeriodScale> compute_10p6_unpol_scale_factors(const CsvTable& table) {
    const std::vector<PeriodInput> inputs = base_inputs();
    std::map<std::string, PeriodRatioAccumulator> acc_by_period;
    for (const auto& input : inputs) acc_by_period[input.period] = PeriodRatioAccumulator();

    for (const auto& row : table.rows) {
        const TupleValue ref = get_tuple(table, row, reference_column());
        if (!ref.ok || std::abs(ref.value) <= 0.0 || !std::isfinite(ref.value)) continue;
        for (const auto& input : inputs) {
            const TupleValue v = get_tuple(table, row, input.column);
            if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) continue;
            const double ratio = v.value / ref.value;
            const double ratio_stat = std::abs(v.stat / ref.value);
            if (!std::isfinite(ratio) || !std::isfinite(ratio_stat) || ratio_stat <= 0.0) continue;
            const double w = 1.0 / (ratio_stat * ratio_stat);
            PeriodRatioAccumulator& acc = acc_by_period[input.period];
            acc.sum_w += w;
            acc.sum_wr += w * ratio;
            acc.n += 1;
        }
    }

    std::vector<PeriodScale> scales;
    std::cout << "[combination-ptp] 10.6 GeV unpol scale factors used before point-to-point estimate:\n";
    for (const auto& input : inputs) {
        const PeriodRatioAccumulator& acc = acc_by_period[input.period];
        if (acc.sum_w <= 0.0 || acc.n <= 0) throw std::runtime_error("Could not compute scale factor for period: " + input.period);
        PeriodScale scale;
        scale.period = input.period;
        scale.scale = acc.sum_wr / acc.sum_w;
        scale.scale_stat = 1.0 / std::sqrt(acc.sum_w);
        scale.n = acc.n;
        if (!std::isfinite(scale.scale) || std::abs(scale.scale) <= 0.0) throw std::runtime_error("Invalid scale factor for period: " + input.period);
        std::cout << "  " << scale.period << " scale = " << std::setprecision(10) << scale.scale
                  << " +/- " << scale.scale_stat << "   n=" << scale.n << "\n";
        scales.push_back(scale);
    }
    return scales;
}

static std::map<std::string, double> scale_map_from_vector(const std::vector<PeriodScale>& scales) {
    std::map<std::string, double> out;
    for (const auto& s : scales) out[s.period] = s.scale;
    return out;
}

static TupleValue apply_scale(const TupleValue& input, double scale) {
    TupleValue out;
    if (!input.ok || !std::isfinite(scale) || std::abs(scale) <= 0.0) return out;
    out.ok = true;
    out.value = input.value / scale;
    out.stat = input.stat / std::abs(scale);
    out.sys = input.sys / std::abs(scale);
    return out;
}

static double chi2_ndf_with_fraction_against_ref(const std::vector<TupleValue>& values,
                                                 double f,
                                                 double ref_value) {
    if (values.size() < 2U || !std::isfinite(ref_value) || std::abs(ref_value) <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double extra = f * std::abs(ref_value);
    double chi2 = 0.0;
    for (const auto& v : values) {
        if (!v.ok || v.stat <= 0.0 || !std::isfinite(v.value)) return std::numeric_limits<double>::quiet_NaN();
        const double var = v.stat * v.stat + extra * extra;
        if (!(var > 0.0) || !std::isfinite(var)) return std::numeric_limits<double>::quiet_NaN();
        const double residual = v.value - ref_value;
        chi2 += residual * residual / var;
    }
    const double ndf = (double)values.size() - 1.0;
    if (ndf <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return chi2 / ndf;
}

static double solve_fraction_for_chi2_unity(const std::vector<TupleValue>& values,
                                            double ref_value,
                                            double chi2_initial) {
    if (!std::isfinite(chi2_initial)) return std::numeric_limits<double>::quiet_NaN();
    if (chi2_initial <= 1.0) return 0.0;

    double lo = 0.0;
    double hi = 0.01;
    for (int iter = 0; iter < 80; ++iter) {
        const double chi2_hi = chi2_ndf_with_fraction_against_ref(values, hi, ref_value);
        if (std::isfinite(chi2_hi) && chi2_hi <= 1.0) break;
        hi *= 2.0;
        if (hi > 10.0) return std::numeric_limits<double>::quiet_NaN();
    }
    for (int iter = 0; iter < 100; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double chi2_mid = chi2_ndf_with_fraction_against_ref(values, mid, ref_value);
        if (!std::isfinite(chi2_mid) || chi2_mid <= 1.0) hi = mid;
        else lo = mid;
    }
    return hi;
}

static int optional_xb_index(const CsvTable& table, const std::vector<std::string>& row) {
    const int col = get_col_optional(table, "xB index");
    if (col < 0) return -1;
    double value = 0.0;
    if (!parse_double(row[(size_t)col], value)) return -1;
    return (int)std::llround(value);
}

static BinResult evaluate_row(const CsvTable& table,
                              const std::vector<std::string>& row,
                              size_t row_index,
                              const std::map<std::string, double>& scale_by_period) {
    BinResult result;
    result.row_index = row_index;
    result.xb_range = Range(get_double_or_nan(table, row, "xBmin"), get_double_or_nan(table, row, "xBmax"));
    result.q2_range = Range(get_double_or_nan(table, row, "Q2min"), get_double_or_nan(table, row, "Q2max"));
    result.t_range = Range(get_double_or_nan(table, row, "t_abs_min"), get_double_or_nan(table, row, "t_abs_max"));
    result.xb_index = optional_xb_index(table, row);
    result.xBavg = get_double_or_nan(table, row, "xBavg, 10.6 GeV");
    result.Q2avg = get_double_or_nan(table, row, "Q2avg, 10.6 GeV");
    result.t_abs_avg = get_double_or_nan(table, row, "t_abs_avg, 10.6 GeV");
    result.phiavg = get_double_or_nan(table, row, "phiavg, 10.6 GeV");
    if (!std::isfinite(result.phiavg)) {
        const double phimin = get_double_or_nan(table, row, "phimin");
        const double phimax = get_double_or_nan(table, row, "phimax");
        if (std::isfinite(phimin) && std::isfinite(phimax)) result.phiavg = 0.5 * (phimin + phimax);
    }

    const TupleValue ref = get_tuple(table, row, reference_column());
    if (!ref.ok || !std::isfinite(ref.value) || std::abs(ref.value) <= 0.0) return result;
    result.ref_value = ref.value;
    result.ref_stat = ref.stat;

    std::vector<TupleValue> scaled_values;
    for (const auto& input : base_inputs()) {
        const TupleValue raw = get_tuple(table, row, input.column);
        const auto it_scale = scale_by_period.find(input.period);
        if (it_scale == scale_by_period.end()) throw std::runtime_error("Missing scale factor for period: " + input.period);
        const TupleValue scaled = apply_scale(raw, it_scale->second);
        if (scaled.ok) scaled_values.push_back(scaled);
    }

    result.n_periods = (int)scaled_values.size();
    if (result.n_periods < 2) return result;

    result.chi2_ndf_before = chi2_ndf_with_fraction_against_ref(scaled_values, 0.0, result.ref_value);
    result.f_ptp = solve_fraction_for_chi2_unity(scaled_values, result.ref_value, result.chi2_ndf_before);
    if (!std::isfinite(result.f_ptp)) return result;
    result.chi2_ndf_after = chi2_ndf_with_fraction_against_ref(scaled_values, result.f_ptp, result.ref_value);
    result.f_ptp_percent = 100.0 * result.f_ptp;
    if (!std::isfinite(result.chi2_ndf_after) || !std::isfinite(result.f_ptp_percent) || !std::isfinite(result.phiavg)) return result;

    result.ok = true;
    return result;
}

static std::vector<BinResult> evaluate_all_bins(const CsvTable& table,
                                                const std::map<std::string, double>& scale_by_period) {
    std::vector<BinResult> results;
    results.reserve(table.rows.size());
    for (size_t i = 0; i < table.rows.size(); ++i) {
        BinResult result = evaluate_row(table, table.rows[i], i + 1, scale_by_period);
        if (result.ok) results.push_back(result);
    }
    return results;
}

static std::string format_double(double value) {
    if (!std::isfinite(value)) return "";
    std::ostringstream oss;
    oss << std::setprecision(12) << value;
    return oss.str();
}

static std::string csv_escape_field(const std::string& s) {
    const bool need_quotes = s.find(',') != std::string::npos || s.find('"') != std::string::npos || s.find('\n') != std::string::npos || s.find('\r') != std::string::npos;
    if (!need_quotes) return s;
    std::ostringstream oss;
    oss << '"';
    for (const char ch : s) {
        if (ch == '"') oss << "\"\"";
        else oss << ch;
    }
    oss << '"';
    return oss.str();
}

static void write_scale_summary_csv(const std::string& path, const std::vector<PeriodScale>& scales) {
    std::ofstream fout(path);
    if (!fout.is_open()) throw std::runtime_error("Could not open scale summary CSV: " + path);
    fout << "period,scale,scale stat,n\n";
    for (const auto& s : scales) {
        fout << csv_escape_field(s.period) << "," << format_double(s.scale) << "," << format_double(s.scale_stat) << "," << s.n << "\n";
    }
    fout.close();
    if (!fout) throw std::runtime_error("Failed while writing scale summary CSV: " + path);
}

static void write_bin_summary_csv(const std::string& path, const std::vector<BinResult>& results) {
    std::ofstream fout(path);
    if (!fout.is_open()) throw std::runtime_error("Could not open point-to-point summary CSV: " + path);
    fout << "row index,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,xBavg,Q2avg,t_abs_avg,phiavg,n periods,reference value,reference stat,chi2ndf before,chi2ndf after,f_ptp,f_ptp percent\n";
    for (const auto& r : results) {
        fout << r.row_index << "," << format_double(r.xb_range.first) << "," << format_double(r.xb_range.second) << ","
             << format_double(r.q2_range.first) << "," << format_double(r.q2_range.second) << ","
             << format_double(r.t_range.first) << "," << format_double(r.t_range.second) << ","
             << format_double(r.xBavg) << "," << format_double(r.Q2avg) << "," << format_double(r.t_abs_avg) << ","
             << format_double(r.phiavg) << "," << r.n_periods << "," << format_double(r.ref_value) << ","
             << format_double(r.ref_stat) << "," << format_double(r.chi2_ndf_before) << "," << format_double(r.chi2_ndf_after) << ","
             << format_double(r.f_ptp) << "," << format_double(r.f_ptp_percent) << "\n";
    }
    fout.close();
    if (!fout) throw std::runtime_error("Failed while writing point-to-point summary CSV: " + path);
}

static std::map<Range, GroupByXB> group_results_by_xb(const std::vector<BinResult>& results) {
    std::map<Range, GroupByXB> grouped;
    for (const auto& r : results) {
        GroupByXB& group = grouped[r.xb_range];
        group.xb_range = r.xb_range;
        if (group.xb_index < 0 && r.xb_index >= 0) group.xb_index = r.xb_index;
        PlotPoint p;
        p.phi = r.phiavg;
        p.y = r.f_ptp_percent;
        group.bins[QTKey(r.q2_range, r.t_range)].points.push_back(p);
    }
    return grouped;
}

static double percentile(std::vector<double> values, double q) {
    std::vector<double> clean;
    for (const double v : values) if (std::isfinite(v)) clean.push_back(v);
    if (clean.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(clean.begin(), clean.end());
    if (q <= 0.0) return clean.front();
    if (q >= 1.0) return clean.back();
    const double pos = q * (double)(clean.size() - 1);
    const size_t lo = (size_t)std::floor(pos);
    const size_t hi = (size_t)std::ceil(pos);
    const double frac = pos - (double)lo;
    if (lo == hi) return clean[lo];
    return clean[lo] * (1.0 - frac) + clean[hi] * frac;
}

static double choose_global_ymax_percent(const std::vector<BinResult>& results) {
    std::vector<double> ys;
    for (const auto& r : results) if (std::isfinite(r.f_ptp_percent) && r.f_ptp_percent >= 0.0) ys.push_back(r.f_ptp_percent);
    if (ys.empty()) return 10.0;
    const double p95 = percentile(ys, 0.95);
    double ymax = std::max(5.0, 1.25 * p95);
    if (!std::isfinite(ymax) || ymax <= 0.0) ymax = 10.0;
    return ymax;
}

static TGraph* make_ptp_graph(const std::vector<PlotPoint>& points) {
    if (points.empty()) return nullptr;
    std::vector<PlotPoint> sorted = points;
    std::sort(sorted.begin(), sorted.end(), [](const PlotPoint& a, const PlotPoint& b) { return a.phi < b.phi; });
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
    if (group.bins.empty()) return;
    std::set<Range> q2_set;
    std::set<Range> t_set;
    for (const auto& kv : group.bins) {
        q2_set.insert(kv.first.first);
        t_set.insert(kv.first.second);
    }
    std::vector<Range> q2_slice(q2_set.begin(), q2_set.end());
    std::vector<Range> t_slice(t_set.begin(), t_set.end());
    if (q2_slice.empty() || t_slice.empty()) return;

    const int ncols = (int)q2_slice.size();
    const int nrows = (int)t_slice.size();
    const int width = 300 * ncols + 160;
    const int height = 260 * nrows + 240;
    TCanvas* canvas = new TCanvas("c_combination_ptp", "c_combination_ptp", width, height);

    TPad* pTop = new TPad("pTopCombinationPtp", "pTopCombinationPtp", 0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();
    TPad* pGrid = new TPad("pGridCombinationPtp", "pGridCombinationPtp", 0.0, 0.00, 1.0, 0.78);
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
    title << "Point-to-point run-period combination systematic, 10.6 GeV unpol   x_{B} in ("
          << std::fixed << std::setprecision(3) << xb_range.first << ", " << xb_range.second << ")";
    head.DrawLatex(0.5, 0.83, title.str().c_str());
    TLatex sub;
    sub.SetNDC();
    sub.SetTextAlign(22);
    sub.SetTextFont(42);
    sub.SetTextSize(0.040);
    sub.DrawLatex(0.5, 0.50, "Reference is the CSV 10.6 GeV combined value; global period scale offsets are removed first");

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
            lab.DrawLatex(0.12, 0.83, Form("Q^{2}=(%.2f, %.2f)  |t|=(%.2f, %.2f)", q2_range.first, q2_range.second, t_range.first, t_range.second));
            const auto it_bin = group.bins.find(QTKey(q2_range, t_range));
            if (it_bin == group.bins.end()) continue;
            TGraph* graph = make_ptp_graph(it_bin->second.points);
            if (graph) graph->Draw("P SAME");
        }
    }
    canvas->cd();
    canvas->Modified();
    canvas->Update();
    const int xb_name = (group.xb_index >= 0 ? group.xb_index : xb_counter);
    std::ostringstream fname;
    fname << "combination_point_to_point_sys_10p6_GeV_unpol_xB_" << xb_name << ".pdf";
    canvas->SaveAs((out_dir / fname.str()).string().c_str());
    delete canvas;
}

static void make_ptp_plots(const std::vector<BinResult>& results, const fs::path& out_dir) {
    const std::map<Range, GroupByXB> grouped = group_results_by_xb(results);
    const double y_max_percent = choose_global_ymax_percent(results);
    int xb_counter = 0;
    for (const auto& kv : grouped) {
        draw_one_xb_canvas(kv.first, kv.second, xb_counter, y_max_percent, out_dir);
        ++xb_counter;
    }
}

static double mean(const std::vector<double>& v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    double sum = 0.0;
    for (const double x : v) sum += x;
    return sum / (double)v.size();
}

static void print_global_summary(const std::vector<BinResult>& results) {
    std::vector<double> fvals;
    std::vector<double> chi2_before;
    std::vector<double> chi2_after;
    for (const auto& r : results) {
        if (std::isfinite(r.f_ptp_percent)) fvals.push_back(r.f_ptp_percent);
        if (std::isfinite(r.chi2_ndf_before)) chi2_before.push_back(r.chi2_ndf_before);
        if (std::isfinite(r.chi2_ndf_after)) chi2_after.push_back(r.chi2_ndf_after);
    }
    std::cout << "[combination-ptp] Summary for 10.6 GeV unpol residual point-to-point systematic\n";
    std::cout << "  valid bins                  = " << results.size() << "\n";
    std::cout << "  mean chi2/ndf before        = " << std::setprecision(8) << mean(chi2_before) << "\n";
    std::cout << "  mean chi2/ndf after         = " << std::setprecision(8) << mean(chi2_after) << "\n";
    std::cout << "  mean f_ptp percent          = " << std::setprecision(8) << mean(fvals) << "\n";
    std::cout << "  median f_ptp percent        = " << std::setprecision(8) << percentile(fvals, 0.50) << "\n";
    std::cout << "  p68 f_ptp percent           = " << std::setprecision(8) << percentile(fvals, 0.68) << "\n";
    std::cout << "  p90 f_ptp percent           = " << std::setprecision(8) << percentile(fvals, 0.90) << "\n";
    std::cout << "  p95 f_ptp percent           = " << std::setprecision(8) << percentile(fvals, 0.95) << "\n";
    std::cout << "  max f_ptp percent           = " << std::setprecision(8)
              << (fvals.empty() ? std::numeric_limits<double>::quiet_NaN() : *std::max_element(fvals.begin(), fvals.end())) << "\n";
}

} // namespace

bool combination_point_to_point_systematics(const std::string& csv_path,
                                            const std::string& output_root_dir) {
    try {
        const fs::path out_dir = fs::path(output_root_dir) / "combination_point_to_point_systematics";
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
        std::cout << "[combination-ptp] CSV rows loaded: " << table.rows.size() << "\n";
        std::cout << "[combination-ptp] Output directory: " << out_dir.string() << "\n";

        const std::vector<PeriodScale> scales = compute_10p6_unpol_scale_factors(table);
        const std::map<std::string, double> scale_by_period = scale_map_from_vector(scales);
        const std::vector<BinResult> results = evaluate_all_bins(table, scale_by_period);

        const fs::path scale_csv = out_dir / "period_scale_factors_used.csv";
        const fs::path bin_csv = out_dir / "combination_point_to_point_systematics_by_bin.csv";
        write_scale_summary_csv(scale_csv.string(), scales);
        write_bin_summary_csv(bin_csv.string(), results);
        make_ptp_plots(results, out_dir);
        print_global_summary(results);

        std::cout << "[combination-ptp] Wrote scale summary: " << scale_csv.string() << "\n";
        std::cout << "[combination-ptp] Wrote bin summary: " << bin_csv.string() << "\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[combination-ptp] FATAL: " << e.what() << "\n";
        return false;
    }
}
