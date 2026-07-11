// sp19_inb_energy_scaling_systematics.cpp
// -----------------------------------------------------------------------------
// Energy-corrected Sp19 Inb normalization diagnostic.
//
// Purpose:
//   Sp19 Inb was taken at about 10.2 GeV, while the main RGA pass-2 DVCS
//   cross-section combination uses the 10.6 GeV periods.  This stage estimates
//   the multiplicative model factor needed to translate the Sp19 Inb unpolarized
//   cross section to the equivalent 10.6 GeV beam-energy point at the same
//   measured bin kinematics.  It then compares the energy-corrected Sp19 Inb
//   result with Fa18 Inb and with the 10.6 GeV combined result as functions of
//   several diagnostic variables.
//
// Model correction:
//   For each CSV row, evaluate BH, VGG and KM15 at E=10.604 GeV and E=10.1998 GeV
//   using the Sp19 Inb average kinematics when available.  The correction factor is
//
//       C_E = average_m [ sigma_m(E=10.604) / sigma_m(E=10.1998) ].
//
//   The RMS spread among the valid model ratios is saved as an energy-correction
//   model systematic diagnostic.
//
// Outputs:
//   - Adds row-level diagnostic columns to the pass-2 CSV.
//   - Writes one consolidated diagnostic CSV under
//       output/systematics/sp19_inb_energy_scaling/
//   - Produces one PNG summary canvas with subplots versus
//     e_theta, p_theta, g_theta, xB, Q2 and |t|.
// -----------------------------------------------------------------------------

#include "sp19_inb_energy_scaling_systematics.h"
#include "model_predictions.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TPad.h>
#include <TBox.h>
#include <TH1D.h>

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
#include <numeric>
#include <atomic>
#include <mutex>
#include <thread>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

static constexpr double kE10p6 = 10.604;
static constexpr double kE10p2 = 10.1998;
static constexpr int kNBinsDiagnostic = 12;
static constexpr int kMaxSp19ScaleWorkers = 6;

static const std::string kColFa18Inb =
    "normed cross sections, ep->epg, exp, Fa18 Inb, unpol";
static const std::string kColSp19Inb =
    "normed cross sections, ep->epg, exp, Sp19 Inb, unpol";
static const std::string kColTenSix =
    "normed cross sections, ep->epg, exp, 10.6 GeV, unpol";

static const std::string kOutEnergyCorrection =
    "Sp19 Inb energy correction factor to 10.6 GeV, model average";
static const std::string kOutEnergyCorrectionRms =
    "Sp19 Inb energy correction factor to 10.6 GeV, model rms";
static const std::string kOutEnergyCorrectionN =
    "Sp19 Inb energy correction factor to 10.6 GeV, model count";
static const std::string kOutSp19Corrected =
    "normed cross sections, ep->epg, exp, Sp19 Inb energy-corrected to 10.6 GeV, unpol";
static const std::string kOutRatioToFa18 =
    "ratio, Sp19 Inb energy-corrected to 10.6 GeV over Fa18 Inb, unpol";
static const std::string kOutRatioToTenSix =
    "ratio, Sp19 Inb energy-corrected to 10.6 GeV over 10.6 GeV, unpol";
static const std::string kOutScaleToTenSix =
    "scale factor, Sp19 Inb energy-corrected to 10.6 GeV to match 10.6 GeV, unpol";
static const std::string kOutScaleToFa18 =
    "scale factor, Sp19 Inb energy-corrected to 10.6 GeV to match Fa18 Inb, unpol";

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

struct RowPoint {
    int row_index = -1;
    double xB = 0.0;
    double Q2 = 0.0;
    double t = 0.0;
    double phi = 0.0;
    double e_theta = 0.0;
    double p_theta = 0.0;
    double g_theta = 0.0;
    double fa18 = 0.0;
    double fa18_stat = 0.0;
    double sp19 = 0.0;
    double sp19_stat = 0.0;
    double ten6 = 0.0;
    double ten6_stat = 0.0;
    double cE = 0.0;
    double cE_rms = 0.0;
    int cE_n = 0;
    double sp19_corr = 0.0;
    double sp19_corr_stat = 0.0;
    double ratio_to_fa18 = 0.0;
    double ratio_to_fa18_stat = 0.0;
    double ratio_to_ten6 = 0.0;
    double ratio_to_ten6_stat = 0.0;
    double scale_to_fa18 = 0.0;
    double scale_to_fa18_stat = 0.0;
    double scale_to_ten6 = 0.0;
    double scale_to_ten6_stat = 0.0;
};

struct BinnedPoint {
    double x = 0.0;
    double ex = 0.0;
    double y = 0.0;

    // Statistical error on the projected weighted mean.  This is usually small
    // and is kept separately so the diagnostic CSV remains quantitatively useful.
    double ey_stat = 0.0;

    // RMS spread of the row-level ratios inside this one-dimensional projection
    // bin.  This is the quantity that makes the projection-composition scatter
    // visible.
    double y_rms = 0.0;

    // Average absolute contribution from the BH/VGG/KM15 energy-correction model
    // spread inside this projection bin.
    double ey_model = 0.0;

    // Error bar used on the summary canvas.  This is intentionally not just the
    // statistical error on the mean; it includes the row-level projected RMS and
    // the model-correction spread so the plot does not overstate the significance
    // of projection wiggles.
    double ey_display = 0.0;

    int n = 0;
};


struct WorkItem {
    int row_index = -1;
    double xB_model = 0.0;
    double Q2_model = 0.0;
    double t_model = 0.0;
    double phi_model = 0.0;
    double xB = 0.0;
    double Q2 = 0.0;
    double t = 0.0;
    double phi = 0.0;
    double e_theta = 0.0;
    double p_theta = 0.0;
    double g_theta = 0.0;
    TupleValue fa18;
    TupleValue sp19;
    TupleValue ten6;
};

struct WorkResult {
    bool ok = false;
    bool missing_model = false;
    WorkItem item;
    RowPoint point;
    double bh_ratio = std::numeric_limits<double>::quiet_NaN();
    double vgg_ratio = std::numeric_limits<double>::quiet_NaN();
    double km15_ratio = std::numeric_limits<double>::quiet_NaN();
};

static std::string trim_copy(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return std::string();
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
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
        const bool quote = s.find(',') != std::string::npos ||
                           s.find('"') != std::string::npos ||
                           s.find('\n') != std::string::npos ||
                           s.find('\r') != std::string::npos;
        if (quote) {
            oss << '"';
            for (const char c : s) {
                if (c == '"') oss << "\"\"";
                else oss << c;
            }
            oss << '"';
        } else {
            oss << s;
        }
        if (i + 1 < fields.size()) oss << ',';
    }
    return oss.str();
}

static std::unordered_map<std::string, int> build_index(const std::vector<std::string>& h) {
    std::unordered_map<std::string, int> idx;
    for (int i = 0; i < (int)h.size(); ++i) idx[h[(size_t)i]] = i;
    return idx;
}

static CsvTable read_csv_or_throw(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Could not open CSV: " + path);
    CsvTable t;
    std::string line;
    if (!std::getline(fin, line)) throw std::runtime_error("CSV is empty: " + path);
    t.header = split_csv_line(line);
    t.index = build_index(t.header);
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        auto row = split_csv_line(line);
        if (row.size() < t.header.size()) row.resize(t.header.size());
        if (row.size() != t.header.size()) {
            std::ostringstream msg;
            msg << "CSV row has " << row.size() << " columns but header has "
                << t.header.size() << ": " << path;
            throw std::runtime_error(msg.str());
        }
        t.rows.push_back(std::move(row));
    }
    return t;
}

static void write_csv_or_throw(const std::string& path, const CsvTable& t) {
    const std::string tmp = path + ".tmp";
    std::ofstream fout(tmp);
    if (!fout.is_open()) throw std::runtime_error("Could not write temporary CSV: " + tmp);
    fout << join_csv_row(t.header) << "\n";
    for (const auto& row : t.rows) fout << join_csv_row(row) << "\n";
    fout.close();
    if (!fout) throw std::runtime_error("Failed writing temporary CSV: " + tmp);
    fs::rename(tmp, path);
}

static bool parse_double(const std::string& s_in, double& v) {
    const std::string s = trim_copy(s_in);
    if (s.empty()) return false;
    try {
        size_t pos = 0;
        v = std::stod(s, &pos);
        if (pos != s.size()) return false;
        return std::isfinite(v);
    } catch (...) {
        return false;
    }
}

static std::string strip_balanced_quotes(std::string s) {
    s = trim_copy(s);
    bool changed = true;
    while (changed && s.size() >= 2) {
        changed = false;
        if ((s.front() == '"' && s.back() == '"') ||
            (s.front() == '\'' && s.back() == '\'')) {
            s = trim_copy(s.substr(1, s.size() - 2));
            changed = true;
        }
    }
    return s;
}

static TupleValue parse_tuple_cell(std::string s) {
    TupleValue tv;
    s = strip_balanced_quotes(s);
    if (s.empty()) return tv;
    if (s.front() == '(' && s.back() == ')') s = s.substr(1, s.size() - 2);
    std::vector<std::string> parts;
    std::string cur;
    std::istringstream iss(s);
    while (std::getline(iss, cur, ',')) parts.push_back(trim_copy(cur));
    if (parts.empty()) return tv;
    double value = 0.0, stat = 0.0, sys = 0.0;
    if (!parse_double(parts[0], value)) return tv;
    if (parts.size() > 1) parse_double(parts[1], stat);
    if (parts.size() > 2) parse_double(parts[2], sys);
    if (!std::isfinite(value)) return tv;
    tv.ok = true;
    tv.value = value;
    tv.stat = std::fabs(stat);
    tv.sys = std::fabs(sys);
    return tv;
}

static std::string format_double(double x, int precision = 10) {
    if (!std::isfinite(x)) return std::string();
    std::ostringstream oss;
    oss << std::setprecision(precision) << x;
    return oss.str();
}

static std::string format_tuple(double value, double stat, double sys) {
    std::ostringstream oss;
    oss << "(" << std::setprecision(12) << value
        << "," << std::setprecision(12) << stat
        << "," << std::setprecision(12) << sys << ")";
    return oss.str();
}

static int require_column(const CsvTable& t, const std::string& c) {
    auto it = t.index.find(c);
    if (it == t.index.end()) throw std::runtime_error("Missing required column: " + c);
    return it->second;
}

static int find_column(const CsvTable& t, const std::string& c) {
    auto it = t.index.find(c);
    return (it == t.index.end()) ? -1 : it->second;
}

static int ensure_column(CsvTable& t, const std::string& c) {
    auto it = t.index.find(c);
    if (it != t.index.end()) return it->second;
    const int idx = (int)t.header.size();
    t.index[c] = idx;
    t.header.push_back(c);
    for (auto& row : t.rows) row.push_back(std::string());
    return idx;
}

static bool get_double_col(const CsvTable& t, const std::vector<std::string>& row,
                           const std::string& c, double& v) {
    const int idx = find_column(t, c);
    if (idx < 0 || idx >= (int)row.size()) return false;
    return parse_double(row[(size_t)idx], v);
}

static bool get_tuple_col(const CsvTable& t, const std::vector<std::string>& row,
                          const std::string& c, TupleValue& tv) {
    const int idx = find_column(t, c);
    if (idx < 0 || idx >= (int)row.size()) return false;
    tv = parse_tuple_cell(row[(size_t)idx]);
    return tv.ok;
}

static double midpoint_from_cols(const CsvTable& t, const std::vector<std::string>& row,
                                 const std::string& lo_col, const std::string& hi_col,
                                 double fallback = std::numeric_limits<double>::quiet_NaN()) {
    double lo = 0.0, hi = 0.0;
    if (get_double_col(t, row, lo_col, lo) && get_double_col(t, row, hi_col, hi)) {
        return 0.5 * (lo + hi);
    }
    return fallback;
}

static double get_period_or_bin_average(const CsvTable& t,
                                        const std::vector<std::string>& row,
                                        const std::string& base,
                                        const std::string& period,
                                        const std::string& min_col,
                                        const std::string& max_col) {
    double v = 0.0;
    if (get_double_col(t, row, base + ", " + period, v)) return v;
    if (get_double_col(t, row, base + ", 10.6 GeV", v)) return v;
    return midpoint_from_cols(t, row, min_col, max_col);
}

static double get_average_pair_or_fallback(const CsvTable& t,
                                           const std::vector<std::string>& row,
                                           const std::string& base,
                                           const std::string& min_col,
                                           const std::string& max_col) {
    double a = 0.0, b = 0.0;
    bool oka = get_double_col(t, row, base + ", Fa18 Inb", a);
    bool okb = get_double_col(t, row, base + ", Sp19 Inb", b);
    if (oka && okb) return 0.5 * (a + b);
    if (oka) return a;
    if (okb) return b;
    double ten = 0.0;
    if (get_double_col(t, row, base + ", 10.6 GeV", ten)) return ten;
    return midpoint_from_cols(t, row, min_col, max_col);
}

static double ratio_stat_error(double ratio, double a, double da, double b, double db) {
    if (!(std::isfinite(ratio) && a > 0.0 && b > 0.0)) return 0.0;
    const double ra = (da > 0.0) ? da / a : 0.0;
    const double rb = (db > 0.0) ? db / b : 0.0;
    return std::fabs(ratio) * std::sqrt(ra * ra + rb * rb);
}

static std::string model_key(const std::string& model,
                             double xB, double Q2, double t, double phi, double ebeam) {
    std::ostringstream oss;
    oss << model << ","
        << std::fixed << std::setprecision(6)
        << xB << "," << Q2 << "," << t << "," << phi << "," << ebeam;
    return oss.str();
}

static std::unordered_map<std::string, double> read_model_cache(const std::string& path) {
    std::unordered_map<std::string, double> cache;
    std::ifstream fin(path);
    if (!fin.is_open()) return cache;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line.rfind("model,", 0) == 0) continue;
        auto parts = split_csv_line(line);
        if (parts.size() < 7) continue;
        const std::string key = parts[0] + "," + parts[1] + "," + parts[2] + "," +
                                parts[3] + "," + parts[4] + "," + parts[5];
        double v = 0.0;
        if (parse_double(parts[6], v)) cache[key] = v;
    }
    return cache;
}

static void write_model_cache(const std::string& path,
                              const std::unordered_map<std::string, double>& cache) {
    std::ofstream fout(path);
    if (!fout.is_open()) throw std::runtime_error("Could not write model cache: " + path);
    fout << "model,xB,Q2,t_abs,phi_deg,Ebeam,xs\n";
    for (const auto& kv : cache) {
        fout << kv.first << "," << std::setprecision(12) << kv.second << "\n";
    }
}

static bool model_xs_uncached(const std::string& model,
                              double xB, double Q2, double t, double phi, double ebeam,
                              double& xs) {
    double val = 0.0;
    try {
        if (model == "BH") {
            val = vgg_bh_only(xB, Q2, t, phi, ebeam, ModelPaths(), false);
        } else if (model == "VGG") {
            val = vgg_xs(xB, Q2, t, phi, ebeam, Helicity::Unpol, ModelPaths(), false);
        } else if (model == "KM15") {
            val = km15_xs(xB, Q2, t, phi, ebeam, Helicity::Unpol, ModelPaths());
        }
    } catch (const std::exception& e) {
        std::cerr << "[sp19-scale] WARNING: " << model << " model failed at "
                  << "xB=" << xB << ", Q2=" << Q2 << ", t=" << t
                  << ", phi=" << phi << ", E=" << ebeam << ": " << e.what() << "\n";
        val = 0.0;
    }

    xs = val;
    return std::isfinite(xs) && xs > 0.0;
}

static bool model_xs_cached_threadsafe(const std::string& model,
                                       double xB, double Q2, double t, double phi, double ebeam,
                                       std::unordered_map<std::string, double>& cache,
                                       std::mutex& cache_mutex,
                                       double& xs) {
    const std::string key = model_key(model, xB, Q2, t, phi, ebeam);
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        auto it = cache.find(key);
        if (it != cache.end()) {
            xs = it->second;
            return std::isfinite(xs) && xs > 0.0;
        }
    }

    double val = 0.0;
    const bool ok = model_xs_uncached(model, xB, Q2, t, phi, ebeam, val);
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        cache[key] = val;
    }
    xs = val;
    return ok;
}

static bool energy_correction_factor_threadsafe(double xB, double Q2, double t, double phi,
                                                std::unordered_map<std::string, double>& cache,
                                                std::mutex& cache_mutex,
                                                double& avg, double& rms, int& nvalid,
                                                double& bh_ratio, double& vgg_ratio, double& km15_ratio) {
    avg = 0.0;
    rms = 0.0;
    nvalid = 0;
    bh_ratio = vgg_ratio = km15_ratio = std::numeric_limits<double>::quiet_NaN();

    std::vector<double> ratios;
    for (const std::string model : {std::string("BH"), std::string("VGG"), std::string("KM15")}) {
        double xs_hi = 0.0, xs_lo = 0.0;
        const bool ok_hi = model_xs_cached_threadsafe(model, xB, Q2, t, phi, kE10p6, cache, cache_mutex, xs_hi);
        const bool ok_lo = model_xs_cached_threadsafe(model, xB, Q2, t, phi, kE10p2, cache, cache_mutex, xs_lo);
        if (ok_hi && ok_lo && xs_lo > 0.0) {
            const double r = xs_hi / xs_lo;
            if (std::isfinite(r) && r > 0.0) {
                ratios.push_back(r);
                if (model == "BH") bh_ratio = r;
                if (model == "VGG") vgg_ratio = r;
                if (model == "KM15") km15_ratio = r;
            }
        }
    }

    if (ratios.empty()) return false;
    avg = std::accumulate(ratios.begin(), ratios.end(), 0.0) / (double)ratios.size();
    double ss = 0.0;
    for (const double r : ratios) ss += (r - avg) * (r - avg);
    rms = std::sqrt(ss / (double)ratios.size());
    nvalid = (int)ratios.size();
    return std::isfinite(avg) && avg > 0.0;
}

static WorkResult process_sp19_scale_item(const WorkItem& item,
                                          std::unordered_map<std::string, double>& model_cache,
                                          std::mutex& cache_mutex) {
    WorkResult res;
    res.item = item;

    double bh_r = 0.0, vgg_r = 0.0, km15_r = 0.0;
    double cE = 0.0, cE_rms = 0.0;
    int cE_n = 0;
    if (!energy_correction_factor_threadsafe(item.xB_model, item.Q2_model, item.t_model, item.phi_model,
                                             model_cache, cache_mutex,
                                             cE, cE_rms, cE_n,
                                             bh_r, vgg_r, km15_r)) {
        res.missing_model = true;
        return res;
    }

    RowPoint p;
    p.row_index = item.row_index;
    p.xB = item.xB;
    p.Q2 = item.Q2;
    p.t = item.t;
    p.phi = item.phi;
    p.e_theta = item.e_theta;
    p.p_theta = item.p_theta;
    p.g_theta = item.g_theta;
    p.fa18 = item.fa18.value;
    p.fa18_stat = item.fa18.stat;
    p.sp19 = item.sp19.value;
    p.sp19_stat = item.sp19.stat;
    p.ten6 = item.ten6.value;
    p.ten6_stat = item.ten6.stat;
    p.cE = cE;
    p.cE_rms = cE_rms;
    p.cE_n = cE_n;
    p.sp19_corr = item.sp19.value * cE;
    p.sp19_corr_stat = item.sp19.stat * cE;
    p.ratio_to_fa18 = p.sp19_corr / item.fa18.value;
    p.ratio_to_fa18_stat = ratio_stat_error(p.ratio_to_fa18, p.sp19_corr, p.sp19_corr_stat, item.fa18.value, item.fa18.stat);
    p.ratio_to_ten6 = p.sp19_corr / item.ten6.value;
    p.ratio_to_ten6_stat = ratio_stat_error(p.ratio_to_ten6, p.sp19_corr, p.sp19_corr_stat, item.ten6.value, item.ten6.stat);
    p.scale_to_fa18 = item.fa18.value / p.sp19_corr;
    p.scale_to_fa18_stat = ratio_stat_error(p.scale_to_fa18, item.fa18.value, item.fa18.stat, p.sp19_corr, p.sp19_corr_stat);
    p.scale_to_ten6 = item.ten6.value / p.sp19_corr;
    p.scale_to_ten6_stat = ratio_stat_error(p.scale_to_ten6, item.ten6.value, item.ten6.stat, p.sp19_corr, p.sp19_corr_stat);

    res.ok = true;
    res.point = p;
    res.bh_ratio = bh_r;
    res.vgg_ratio = vgg_r;
    res.km15_ratio = km15_r;
    return res;
}

static std::vector<BinnedPoint> make_binned_points(const std::vector<RowPoint>& rows,
                                                   const std::string& variable,
                                                   const std::string& quantity) {
    std::vector<double> xs;
    xs.reserve(rows.size());
    auto get_x = [&](const RowPoint& p) {
        if (variable == "xB") return p.xB;
        if (variable == "Q2") return p.Q2;
        if (variable == "t") return p.t;
        if (variable == "e_theta") return p.e_theta;
        if (variable == "p_theta") return p.p_theta;
        if (variable == "g_theta") return p.g_theta;
        return std::numeric_limits<double>::quiet_NaN();
    };
    auto get_y = [&](const RowPoint& p) {
        if (quantity == "ratio_to_fa18") return p.ratio_to_fa18;
        if (quantity == "ratio_to_ten6") return p.ratio_to_ten6;
        if (quantity == "scale_to_fa18") return p.scale_to_fa18;
        if (quantity == "scale_to_ten6") return p.scale_to_ten6;
        return std::numeric_limits<double>::quiet_NaN();
    };
    auto get_ey_stat = [&](const RowPoint& p) {
        if (quantity == "ratio_to_fa18") return p.ratio_to_fa18_stat;
        if (quantity == "ratio_to_ten6") return p.ratio_to_ten6_stat;
        if (quantity == "scale_to_fa18") return p.scale_to_fa18_stat;
        if (quantity == "scale_to_ten6") return p.scale_to_ten6_stat;
        return 0.0;
    };
    auto get_model_abs = [&](const RowPoint& p) {
        const double y = get_y(p);
        const double rel = (p.cE > 0.0 && std::isfinite(p.cE_rms)) ? std::fabs(p.cE_rms / p.cE) : 0.0;
        return (std::isfinite(y) && y > 0.0) ? std::fabs(y) * rel : 0.0;
    };

    for (const auto& p : rows) {
        const double x = get_x(p);
        const double y = get_y(p);
        if (std::isfinite(x) && std::isfinite(y) && y > 0.0) xs.push_back(x);
    }
    if (xs.empty()) return {};
    std::sort(xs.begin(), xs.end());
    double xmin = xs.front();
    double xmax = xs.back();
    if (!(xmax > xmin)) return {};
    const double width = (xmax - xmin) / (double)kNBinsDiagnostic;

    std::vector<BinnedPoint> out;
    for (int ib = 0; ib < kNBinsDiagnostic; ++ib) {
        const double lo = xmin + ib * width;
        const double hi = (ib + 1 == kNBinsDiagnostic) ? xmax + 1e-12 : xmin + (ib + 1) * width;
        double sw = 0.0, swy = 0.0, swx = 0.0;
        int n = 0;
        std::vector<double> y_values;
        std::vector<double> weights;
        double model_abs_sum = 0.0;
        for (const auto& p : rows) {
            const double x = get_x(p);
            const double y = get_y(p);
            double ey_stat = get_ey_stat(p);
            if (!(std::isfinite(x) && std::isfinite(y) && y > 0.0)) continue;
            if (x < lo || x >= hi) continue;
            if (!(std::isfinite(ey_stat) && ey_stat > 0.0)) ey_stat = 1.0;
            const double w = 1.0 / (ey_stat * ey_stat);
            sw += w;
            swy += w * y;
            swx += w * x;
            y_values.push_back(y);
            weights.push_back(w);
            model_abs_sum += get_model_abs(p);
            ++n;
        }
        if (n <= 0 || sw <= 0.0) continue;
        BinnedPoint bp;
        bp.x = swx / sw;
        bp.ex = 0.5 * width;
        bp.y = swy / sw;
        bp.ey_stat = std::sqrt(1.0 / sw);
        bp.n = n;

        // Weighted RMS of the row-level ratios in this projected bin.  This is a
        // spread diagnostic, not an error on the mean.
        double wrms_num = 0.0;
        for (size_t i = 0; i < y_values.size(); ++i) {
            const double d = y_values[i] - bp.y;
            wrms_num += weights[i] * d * d;
        }
        bp.y_rms = (sw > 0.0) ? std::sqrt(wrms_num / sw) : 0.0;
        bp.ey_model = (n > 0) ? model_abs_sum / (double)n : 0.0;
        bp.ey_display = std::sqrt(bp.ey_stat * bp.ey_stat + bp.y_rms * bp.y_rms + bp.ey_model * bp.ey_model);
        out.push_back(bp);
    }
    return out;
}

static double finite_or_nan(double x) {
    return std::isfinite(x) ? x : std::numeric_limits<double>::quiet_NaN();
}

static double graph_y_min_from_central(const std::vector<BinnedPoint>& pts) {
    double out = std::numeric_limits<double>::infinity();
    for (const auto& p : pts) {
        if (!std::isfinite(p.y)) continue;
        const double e = std::isfinite(p.ey_stat) ? p.ey_stat : 0.0;
        out = std::min(out, p.y - e);
    }
    return out;
}

static double graph_y_max_from_central(const std::vector<BinnedPoint>& pts) {
    double out = -std::numeric_limits<double>::infinity();
    for (const auto& p : pts) {
        if (!std::isfinite(p.y)) continue;
        const double e = std::isfinite(p.ey_stat) ? p.ey_stat : 0.0;
        out = std::max(out, p.y + e);
    }
    return out;
}

static double projection_band_half_height(const BinnedPoint& p) {
    const double rms = std::isfinite(p.y_rms) ? p.y_rms : 0.0;
    const double model = std::isfinite(p.ey_model) ? p.ey_model : 0.0;
    return std::sqrt(rms * rms + model * model);
}

static void draw_projection_boxes(const std::vector<BinnedPoint>& pts,
                                  int color,
                                  double alpha,
                                  double xmin,
                                  double xmax) {
    for (const auto& p : pts) {
        if (!(std::isfinite(p.x) && std::isfinite(p.y))) continue;
        const double band = projection_band_half_height(p);
        if (!(std::isfinite(band) && band > 0.0)) continue;
        double half_width = std::isfinite(p.ex) && p.ex > 0.0 ? p.ex : 0.0;
        if (half_width <= 0.0) half_width = 0.005 * std::max(1.0, xmax - xmin);
        const double xlo = std::max(xmin, p.x - half_width);
        const double xhi = std::min(xmax, p.x + half_width);
        if (!(xhi > xlo)) continue;
        TBox* box = new TBox(xlo, p.y - band, xhi, p.y + band);
        box->SetFillColorAlpha(color, alpha);
        box->SetLineColorAlpha(color, 0.0);
        box->Draw("same");
    }
}

static TGraphErrors* make_stat_graph(const std::vector<BinnedPoint>& pts,
                                     int marker,
                                     int color) {
    std::vector<double> x, ex, y, ey;
    x.reserve(pts.size()); ex.reserve(pts.size()); y.reserve(pts.size()); ey.reserve(pts.size());
    for (const auto& p : pts) {
        if (!(std::isfinite(p.x) && std::isfinite(p.y))) continue;
        x.push_back(p.x);
        ex.push_back(0.0);
        y.push_back(p.y);
        ey.push_back(std::isfinite(p.ey_stat) ? p.ey_stat : 0.0);
    }
    if (x.empty()) return nullptr;
    TGraphErrors* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(0.72);
    g->SetLineWidth(2);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    return g;
}

static void draw_sp19_kinematic_summary_canvas(
        const std::string& out_png,
        const std::vector<std::pair<std::string, std::string> >& variables,
        const std::map<std::string, std::vector<BinnedPoint> >& binned_fa18,
        const std::map<std::string, std::vector<BinnedPoint> >& binned_ten6) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    TCanvas c("c_sp19_scale_summary", "c_sp19_scale_summary", 1900, 1180);
    c.Divide(3, 2, 0.006, 0.006);

    for (size_t iv = 0; iv < variables.size(); ++iv) {
        c.cd((int)iv + 1);
        TPad* pad = (TPad*)gPad;
        pad->SetLeftMargin(0.13);
        pad->SetRightMargin(0.04);
        pad->SetBottomMargin(0.14);
        pad->SetTopMargin(0.12);
        pad->SetGridy(true);

        const std::string& key = variables[iv].first;
        const std::string& label = variables[iv].second;
        const auto it_f = binned_fa18.find(key);
        const auto it_t = binned_ten6.find(key);
        const std::vector<BinnedPoint> empty;
        const auto& fa18 = (it_f == binned_fa18.end()) ? empty : it_f->second;
        const auto& ten6 = (it_t == binned_ten6.end()) ? empty : it_t->second;

        if (fa18.empty() && ten6.empty()) {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.050);
            lat.DrawLatex(0.20, 0.52, ("No points for " + key).c_str());
            continue;
        }

        double xmin = std::numeric_limits<double>::infinity();
        double xmax = -std::numeric_limits<double>::infinity();
        auto scan_x = [&](const std::vector<BinnedPoint>& pts) {
            for (const auto& p : pts) {
                if (!std::isfinite(p.x)) continue;
                const double ex = std::isfinite(p.ex) && p.ex > 0.0 ? p.ex : 0.0;
                xmin = std::min(xmin, p.x - ex);
                xmax = std::max(xmax, p.x + ex);
            }
        };
        scan_x(fa18);
        scan_x(ten6);
        if (!std::isfinite(xmin) || !std::isfinite(xmax) || !(xmax > xmin)) {
            xmin = 0.0; xmax = 1.0;
        }
        const double xrange = xmax - xmin;
        xmin -= 0.03 * xrange;
        xmax += 0.03 * xrange;

        // The central values are the main diagnostic.  The shaded boxes show the
        // row-level projection spread and model-correction spread, which can be
        // large because a one-dimensional projection mixes several other
        // kinematic and topology variables.  Use the central-value range for the
        // frame so the trend is readable, while allowing the boxes to extend to
        // the plot boundary if the projection spread is large.
        double ycentral_min = std::min(graph_y_min_from_central(fa18), graph_y_min_from_central(ten6));
        double ycentral_max = std::max(graph_y_max_from_central(fa18), graph_y_max_from_central(ten6));
        if (!std::isfinite(ycentral_min) || !std::isfinite(ycentral_max) || !(ycentral_max > ycentral_min)) {
            ycentral_min = 0.75;
            ycentral_max = 1.05;
        }
        ycentral_min = std::min(ycentral_min, 1.0);
        ycentral_max = std::max(ycentral_max, 1.0);
        double ymin = std::max(0.45, ycentral_min - 0.12);
        double ymax = std::min(1.25, ycentral_max + 0.12);
        if (!(ymax > ymin + 0.20)) {
            const double mid = 0.5 * (ymax + ymin);
            ymin = mid - 0.15;
            ymax = mid + 0.15;
        }
        ymin = std::min(ymin, 0.90);
        ymax = std::max(ymax, 1.10);

        const std::string frame_name = "frame_sp19_scale_" + key + "_" + std::to_string(iv);
        TH1D* frame = new TH1D(frame_name.c_str(), "", 1, xmin, xmax);
        frame->SetDirectory(nullptr);
        frame->SetMinimum(ymin);
        frame->SetMaximum(ymax);
        frame->GetXaxis()->SetTitle(label.c_str());
        frame->GetYaxis()->SetTitle("Ratio");
        frame->GetXaxis()->SetTitleSize(0.050);
        frame->GetYaxis()->SetTitleSize(0.050);
        frame->GetXaxis()->SetLabelSize(0.040);
        frame->GetYaxis()->SetLabelSize(0.040);
        frame->GetYaxis()->SetTitleOffset(1.15);
        frame->Draw("axis");

        TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        draw_projection_boxes(fa18, kBlue + 1, 0.14, xmin, xmax);
        draw_projection_boxes(ten6, kRed + 1, 0.14, xmin, xmax);

        TGraphErrors* g_fa18 = make_stat_graph(fa18, 20, kBlue + 1);
        TGraphErrors* g_ten6 = make_stat_graph(ten6, 21, kRed + 1);
        if (g_fa18) g_fa18->Draw("P same");
        if (g_ten6) g_ten6->Draw("P same");

        frame->Draw("axis same");

        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.048);
        lat.DrawLatex(0.16, 0.925, key.c_str());

        if (iv == 0) {
            TLegend* leg = new TLegend(0.38, 0.66, 0.95, 0.90);
            leg->SetBorderSize(1);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetTextSize(0.030);
            TBox* blue_box = new TBox(0.0, 0.0, 1.0, 1.0);
            blue_box->SetFillColorAlpha(kBlue + 1, 0.18);
            blue_box->SetLineColor(kBlue + 1);
            TBox* red_box = new TBox(0.0, 0.0, 1.0, 1.0);
            red_box->SetFillColorAlpha(kRed + 1, 0.18);
            red_box->SetLineColor(kRed + 1);
            if (g_fa18) leg->AddEntry(g_fa18, "mean: Sp19 corr / Fa18 Inb", "p");
            if (g_ten6) leg->AddEntry(g_ten6, "mean: Sp19 corr / 10.6 GeV", "p");
            leg->AddEntry(blue_box, "blue box: projection RMS #oplus model", "f");
            leg->AddEntry(red_box, "red box: projection RMS #oplus model", "f");
            leg->Draw("same");
        }

        if (iv == 1) {
            TLatex note;
            note.SetNDC(true);
            note.SetTextSize(0.030);
            note.DrawLatex(0.16, 0.83, "markers: weighted projected means");
            note.DrawLatex(0.16, 0.77, "thin bars: stat. error on mean");
            note.DrawLatex(0.16, 0.71, "boxes: row RMS + energy-model spread");
        }
    }

    c.SaveAs(out_png.c_str());
}

static void write_sp19_diagnostics_csv(const std::string& path,
                                       const std::vector<RowPoint>& pts,
                                       const std::vector<std::pair<std::string, std::string> >& variables,
                                       const std::map<std::string, std::vector<BinnedPoint> >& binned_fa18,
                                       const std::map<std::string, std::vector<BinnedPoint> >& binned_ten6,
                                       size_t cache_before,
                                       size_t cache_after,
                                       int n_missing_data,
                                       int n_missing_model,
                                       double mean_cE,
                                       double mean_cE_rms_rel,
                                       const std::pair<double,double>& ratio_fa18,
                                       const std::pair<double,double>& ratio_ten6,
                                       const std::pair<double,double>& scale_fa18,
                                       const std::pair<double,double>& scale_ten6) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("Could not write diagnostics CSV: " + path);

    const std::vector<std::string> header = {
        "record_type", "quantity", "variable", "bin_index",
        "x", "x_half_width", "value", "stat_error_on_mean", "projected_row_rms",
        "energy_model_abs_error", "display_error", "n_points", "row",
        "xB", "Q2", "t_abs", "phi_deg", "e_theta", "p_theta", "g_theta",
        "Fa18_Inb", "Fa18_Inb_stat", "Sp19_Inb", "Sp19_Inb_stat", "Ten6", "Ten6_stat",
        "energy_correction_avg", "energy_correction_rms", "energy_correction_rms_rel", "energy_model_count",
        "Sp19_energy_corrected", "Sp19_energy_corrected_stat",
        "ratio_sp19corr_over_fa18", "ratio_sp19corr_over_fa18_stat",
        "ratio_sp19corr_over_10p6", "ratio_sp19corr_over_10p6_stat",
        "scale_fa18_over_sp19corr", "scale_fa18_over_sp19corr_stat",
        "scale_10p6_over_sp19corr", "scale_10p6_over_sp19corr_stat", "notes"
    };
    f << join_csv_row(header) << "\n";

    auto blank_row = [&]() { return std::vector<std::string>(header.size(), std::string()); };
    auto emit = [&](const std::vector<std::string>& row) { f << join_csv_row(row) << "\n"; };
    auto set = [&](std::vector<std::string>& row, const std::string& col, const std::string& value) {
        auto it = std::find(header.begin(), header.end(), col);
        if (it != header.end()) row[(size_t)std::distance(header.begin(), it)] = value;
    };

    auto write_global = [&](const std::string& quantity, double value, double stat, const std::string& notes) {
        auto row = blank_row();
        set(row, "record_type", "global");
        set(row, "quantity", quantity);
        set(row, "value", format_double(value, 12));
        set(row, "stat_error_on_mean", format_double(stat, 12));
        set(row, "notes", notes);
        emit(row);
    };

    write_global("n_valid_rows", (double)pts.size(), 0.0, "Rows with Fa18 Inb, Sp19 Inb, 10.6 GeV and valid model correction");
    write_global("n_rows_missing_data", (double)n_missing_data, 0.0, "Rows skipped before model correction");
    write_global("n_rows_missing_model", (double)n_missing_model, 0.0, "Rows skipped because all model corrections failed");
    write_global("model_cache_entries_before", (double)cache_before, 0.0, "Model cache size before this run");
    write_global("model_cache_entries_after", (double)cache_after, 0.0, "Model cache size after this run");
    write_global("mean_energy_correction_factor", mean_cE, 0.0, "Average of per-bin BH/VGG/KM15 average factors");
    write_global("mean_relative_model_rms", mean_cE_rms_rel, 0.0, "Average over rows of RMS(BH,VGG,KM15)/mean");
    write_global("weighted_mean_sp19corr_over_fa18", ratio_fa18.first, ratio_fa18.second, "Energy-corrected Sp19 divided by Fa18 Inb");
    write_global("weighted_mean_sp19corr_over_10p6", ratio_ten6.first, ratio_ten6.second, "Energy-corrected Sp19 divided by 10.6 GeV combined");
    write_global("recommended_scale_sp19corr_to_fa18", scale_fa18.first, scale_fa18.second, "Fa18 Inb divided by energy-corrected Sp19");
    write_global("recommended_scale_sp19corr_to_10p6", scale_ten6.first, scale_ten6.second, "10.6 GeV combined divided by energy-corrected Sp19");

    for (const auto& kv : variables) {
        const std::string& v = kv.first;
        auto write_binned = [&](const std::string& quantity, const std::vector<BinnedPoint>& bins) {
            for (size_t ib = 0; ib < bins.size(); ++ib) {
                const auto& p = bins[ib];
                auto row = blank_row();
                set(row, "record_type", "binned");
                set(row, "quantity", quantity);
                set(row, "variable", v);
                set(row, "bin_index", std::to_string(ib));
                set(row, "x", format_double(p.x, 12));
                set(row, "x_half_width", format_double(p.ex, 12));
                set(row, "value", format_double(p.y, 12));
                set(row, "stat_error_on_mean", format_double(p.ey_stat, 12));
                set(row, "projected_row_rms", format_double(p.y_rms, 12));
                set(row, "energy_model_abs_error", format_double(p.ey_model, 12));
                set(row, "display_error", format_double(p.ey_display, 12));
                set(row, "n_points", std::to_string(p.n));
                set(row, "notes", "display_error=sqrt(stat_error_on_mean^2+projected_row_rms^2+energy_model_abs_error^2)");
                emit(row);
            }
        };
        auto itf = binned_fa18.find(v);
        if (itf != binned_fa18.end()) write_binned("ratio_sp19corr_over_fa18", itf->second);
        auto itt = binned_ten6.find(v);
        if (itt != binned_ten6.end()) write_binned("ratio_sp19corr_over_10p6", itt->second);
    }

    for (const auto& p : pts) {
        auto row = blank_row();
        set(row, "record_type", "row");
        set(row, "quantity", "row_values");
        set(row, "row", std::to_string(p.row_index));
        set(row, "xB", format_double(p.xB, 12));
        set(row, "Q2", format_double(p.Q2, 12));
        set(row, "t_abs", format_double(p.t, 12));
        set(row, "phi_deg", format_double(p.phi, 12));
        set(row, "e_theta", format_double(p.e_theta, 12));
        set(row, "p_theta", format_double(p.p_theta, 12));
        set(row, "g_theta", format_double(p.g_theta, 12));
        set(row, "Fa18_Inb", format_double(p.fa18, 12));
        set(row, "Fa18_Inb_stat", format_double(p.fa18_stat, 12));
        set(row, "Sp19_Inb", format_double(p.sp19, 12));
        set(row, "Sp19_Inb_stat", format_double(p.sp19_stat, 12));
        set(row, "Ten6", format_double(p.ten6, 12));
        set(row, "Ten6_stat", format_double(p.ten6_stat, 12));
        set(row, "energy_correction_avg", format_double(p.cE, 12));
        set(row, "energy_correction_rms", format_double(p.cE_rms, 12));
        set(row, "energy_correction_rms_rel", format_double((p.cE > 0.0) ? p.cE_rms / p.cE : 0.0, 12));
        set(row, "energy_model_count", std::to_string(p.cE_n));
        set(row, "Sp19_energy_corrected", format_double(p.sp19_corr, 12));
        set(row, "Sp19_energy_corrected_stat", format_double(p.sp19_corr_stat, 12));
        set(row, "ratio_sp19corr_over_fa18", format_double(p.ratio_to_fa18, 12));
        set(row, "ratio_sp19corr_over_fa18_stat", format_double(p.ratio_to_fa18_stat, 12));
        set(row, "ratio_sp19corr_over_10p6", format_double(p.ratio_to_ten6, 12));
        set(row, "ratio_sp19corr_over_10p6_stat", format_double(p.ratio_to_ten6_stat, 12));
        set(row, "scale_fa18_over_sp19corr", format_double(p.scale_to_fa18, 12));
        set(row, "scale_fa18_over_sp19corr_stat", format_double(p.scale_to_fa18_stat, 12));
        set(row, "scale_10p6_over_sp19corr", format_double(p.scale_to_ten6, 12));
        set(row, "scale_10p6_over_sp19corr_stat", format_double(p.scale_to_ten6_stat, 12));
        emit(row);
    }
}

static std::pair<double,double> weighted_mean_and_error(const std::vector<RowPoint>& pts,
                                                        const std::string& quantity) {
    double sw = 0.0, swy = 0.0;
    for (const auto& p : pts) {
        double y = 0.0, ey = 0.0;
        if (quantity == "scale_to_fa18") { y = p.scale_to_fa18; ey = p.scale_to_fa18_stat; }
        else if (quantity == "scale_to_ten6") { y = p.scale_to_ten6; ey = p.scale_to_ten6_stat; }
        else if (quantity == "ratio_to_fa18") { y = p.ratio_to_fa18; ey = p.ratio_to_fa18_stat; }
        else if (quantity == "ratio_to_ten6") { y = p.ratio_to_ten6; ey = p.ratio_to_ten6_stat; }
        if (!(std::isfinite(y) && y > 0.0)) continue;
        if (!(std::isfinite(ey) && ey > 0.0)) ey = 1.0;
        const double w = 1.0 / (ey * ey);
        sw += w;
        swy += w * y;
    }
    if (sw <= 0.0) return {0.0, 0.0};
    return {swy / sw, std::sqrt(1.0 / sw)};
}

} // namespace

bool sp19_inb_energy_scaling_systematics(const std::string& csv_path,
                                         const std::string& output_root_dir) {
    try {
        const fs::path out_dir = fs::path(output_root_dir) / "sp19_inb_energy_scaling";
        fs::create_directories(out_dir);

        std::cout << "[sp19-scale] Starting Sp19 Inb energy-scaling diagnostic.\n";
        std::cout << "[sp19-scale] CSV: " << csv_path << "\n";
        std::cout << "[sp19-scale] Output directory: " << out_dir.string() << "\n";
        std::cout << "[sp19-scale] Beam-energy correction: " << kE10p2
                  << " GeV -> " << kE10p6 << " GeV using BH, VGG and KM15.\n";

        CsvTable table = read_csv_or_throw(csv_path);
        require_column(table, kColFa18Inb);
        require_column(table, kColSp19Inb);
        require_column(table, kColTenSix);

        const int idx_cE = ensure_column(table, kOutEnergyCorrection);
        const int idx_cE_rms = ensure_column(table, kOutEnergyCorrectionRms);
        const int idx_cE_n = ensure_column(table, kOutEnergyCorrectionN);
        const int idx_sp19_corr = ensure_column(table, kOutSp19Corrected);
        const int idx_ratio_fa18 = ensure_column(table, kOutRatioToFa18);
        const int idx_ratio_ten6 = ensure_column(table, kOutRatioToTenSix);
        const int idx_scale_ten6 = ensure_column(table, kOutScaleToTenSix);
        const int idx_scale_fa18 = ensure_column(table, kOutScaleToFa18);

        const std::string cache_path = (out_dir / "model_energy_cache.csv").string();
        auto model_cache = read_model_cache(cache_path);
        const size_t cache_before = model_cache.size();

        std::vector<RowPoint> points;
        points.reserve(table.rows.size());
        int n_missing_data = 0;
        int n_missing_model = 0;

        std::vector<WorkItem> work_items;
        work_items.reserve(table.rows.size());

        for (int ir = 0; ir < (int)table.rows.size(); ++ir) {
            const auto& row = table.rows[(size_t)ir];

            TupleValue fa18, sp19, ten6;
            const bool ok_fa18 = get_tuple_col(table, row, kColFa18Inb, fa18);
            const bool ok_sp19 = get_tuple_col(table, row, kColSp19Inb, sp19);
            const bool ok_ten6 = get_tuple_col(table, row, kColTenSix, ten6);
            if (!ok_fa18 || !ok_sp19 || !ok_ten6 ||
                fa18.value <= 0.0 || sp19.value <= 0.0 || ten6.value <= 0.0) {
                ++n_missing_data;
                continue;
            }

            const double xB_model = get_period_or_bin_average(table, row, "xBavg", "Sp19 Inb", "xBmin", "xBmax");
            const double Q2_model = get_period_or_bin_average(table, row, "Q2avg", "Sp19 Inb", "Q2min", "Q2max");
            const double t_model = get_period_or_bin_average(table, row, "t_abs_avg", "Sp19 Inb", "t_abs_min", "t_abs_max");
            const double phi_model = get_period_or_bin_average(table, row, "phiavg", "Sp19 Inb", "phimin", "phimax");

            if (!(std::isfinite(xB_model) && std::isfinite(Q2_model) && std::isfinite(t_model) && std::isfinite(phi_model)) ||
                xB_model <= 0.0 || Q2_model <= 0.0 || t_model <= 0.0) {
                ++n_missing_data;
                continue;
            }

            WorkItem item;
            item.row_index = ir;
            item.xB_model = xB_model;
            item.Q2_model = Q2_model;
            item.t_model = t_model;
            item.phi_model = phi_model;
            item.xB = get_average_pair_or_fallback(table, row, "xBavg", "xBmin", "xBmax");
            item.Q2 = get_average_pair_or_fallback(table, row, "Q2avg", "Q2min", "Q2max");
            item.t = get_average_pair_or_fallback(table, row, "t_abs_avg", "t_abs_min", "t_abs_max");
            item.phi = get_average_pair_or_fallback(table, row, "phiavg", "phimin", "phimax");
            item.e_theta = get_average_pair_or_fallback(table, row, "e_theta", "", "");
            item.p_theta = get_average_pair_or_fallback(table, row, "p_theta", "", "");
            item.g_theta = get_average_pair_or_fallback(table, row, "g_theta", "", "");
            item.fa18 = fa18;
            item.sp19 = sp19;
            item.ten6 = ten6;
            work_items.push_back(item);
        }

        int requested_workers = kMaxSp19ScaleWorkers;
        if (const char* env_workers = std::getenv("SP19_SCALE_WORKERS")) {
            try {
                requested_workers = std::stoi(env_workers);
            } catch (...) {
                requested_workers = kMaxSp19ScaleWorkers;
            }
        }
        const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
        const int n_workers = std::max(1, std::min(kMaxSp19ScaleWorkers,
                                  std::min(requested_workers, (int)hw)));
        std::cout << "[sp19-scale] Queued " << work_items.size()
                  << " valid Sp19/Fa18 rows for model energy correction."
                  << " Using " << n_workers << " worker(s)"
                  << " (set SP19_SCALE_WORKERS to override, capped at 6).\n";

        std::vector<WorkResult> results(work_items.size());
        std::atomic<size_t> next_index{0};
        std::atomic<size_t> completed{0};
        std::mutex cache_mutex;
        std::mutex cout_mutex;

        auto worker = [&]() {
            while (true) {
                const size_t i = next_index.fetch_add(1);
                if (i >= work_items.size()) break;
                results[i] = process_sp19_scale_item(work_items[i], model_cache, cache_mutex);
                const size_t done = completed.fetch_add(1) + 1;
                if (done == work_items.size() || done % 50 == 0) {
                    std::lock_guard<std::mutex> lock(cout_mutex);
                    std::cout << "[sp19-scale] Model corrections completed "
                              << done << "/" << work_items.size() << " rows...\n";
                }
            }
        };

        std::vector<std::thread> threads;
        threads.reserve((size_t)n_workers);
        for (int iw = 0; iw < n_workers; ++iw) threads.emplace_back(worker);
        for (auto& th : threads) th.join();

        for (const auto& r : results) {
            if (!r.ok) {
                if (r.missing_model) ++n_missing_model;
                continue;
            }
            const RowPoint& p = r.point;
            auto& row = table.rows[(size_t)p.row_index];
            const double cE = p.cE;
            const double cE_rms = p.cE_rms;
            const int cE_n = p.cE_n;

            row[(size_t)idx_cE] = format_double(cE, 12);
            row[(size_t)idx_cE_rms] = format_double(cE_rms, 12);
            row[(size_t)idx_cE_n] = std::to_string(cE_n);
            row[(size_t)idx_sp19_corr] = format_tuple(p.sp19_corr, p.sp19_corr_stat, 0.0);
            row[(size_t)idx_ratio_fa18] = format_tuple(p.ratio_to_fa18, p.ratio_to_fa18_stat, (cE > 0.0 ? cE_rms / cE * p.ratio_to_fa18 : 0.0));
            row[(size_t)idx_ratio_ten6] = format_tuple(p.ratio_to_ten6, p.ratio_to_ten6_stat, (cE > 0.0 ? cE_rms / cE * p.ratio_to_ten6 : 0.0));
            row[(size_t)idx_scale_fa18] = format_tuple(p.scale_to_fa18, p.scale_to_fa18_stat, (cE > 0.0 ? cE_rms / cE * p.scale_to_fa18 : 0.0));
            row[(size_t)idx_scale_ten6] = format_tuple(p.scale_to_ten6, p.scale_to_ten6_stat, (cE > 0.0 ? cE_rms / cE * p.scale_to_ten6 : 0.0));
            points.push_back(p);
        }

        write_model_cache(cache_path, model_cache);
        write_csv_or_throw(csv_path, table);

        const auto scale_fa18 = weighted_mean_and_error(points, "scale_to_fa18");
        const auto scale_ten6 = weighted_mean_and_error(points, "scale_to_ten6");
        const auto ratio_fa18 = weighted_mean_and_error(points, "ratio_to_fa18");
        const auto ratio_ten6 = weighted_mean_and_error(points, "ratio_to_ten6");

        double mean_cE = 0.0, mean_cE_rms_rel = 0.0;
        if (!points.empty()) {
            for (const auto& p : points) {
                mean_cE += p.cE;
                if (p.cE > 0.0) mean_cE_rms_rel += p.cE_rms / p.cE;
            }
            mean_cE /= (double)points.size();
            mean_cE_rms_rel /= (double)points.size();
        }

        const std::vector<std::pair<std::string, std::string> > vars = {
            {"e_theta", "#theta_{e} (deg)"},
            {"p_theta", "#theta_{p} (deg)"},
            {"g_theta", "#theta_{#gamma} (deg)"},
            {"xB", "x_{B}"},
            {"Q2", "Q^{2} (GeV^{2})"},
            {"t", "-t (GeV^{2})"}
        };

        std::map<std::string, std::vector<BinnedPoint> > binned_fa18;
        std::map<std::string, std::vector<BinnedPoint> > binned_ten6;
        for (const auto& v : vars) {
            binned_fa18[v.first] = make_binned_points(points, v.first, "ratio_to_fa18");
            binned_ten6[v.first] = make_binned_points(points, v.first, "ratio_to_ten6");
        }

        write_sp19_diagnostics_csv((out_dir / "sp19_inb_energy_scaling_diagnostics.csv").string(),
                                   points, vars, binned_fa18, binned_ten6,
                                   cache_before, model_cache.size(),
                                   n_missing_data, n_missing_model,
                                   mean_cE, mean_cE_rms_rel,
                                   ratio_fa18, ratio_ten6,
                                   scale_fa18, scale_ten6);

        draw_sp19_kinematic_summary_canvas(
            (out_dir / "sp19_inb_energy_scaling_kinematic_summary.png").string(),
            vars, binned_fa18, binned_ten6);

        std::cout << "[sp19-scale] Valid comparison rows: " << points.size() << "\n";
        std::cout << "[sp19-scale] Rows skipped for missing data: " << n_missing_data << "\n";
        std::cout << "[sp19-scale] Rows skipped for missing model correction: " << n_missing_model << "\n";
        std::cout << "[sp19-scale] Mean energy correction factor 10.2->10.6: "
                  << std::setprecision(6) << mean_cE
                  << "  (mean model RMS/mean = " << mean_cE_rms_rel << ")\n";
        std::cout << "[sp19-scale] Weighted mean Sp19(corrected)/Fa18 Inb = "
                  << ratio_fa18.first << " +/- " << ratio_fa18.second << "\n";
        std::cout << "[sp19-scale] Weighted mean Sp19(corrected)/10.6 GeV = "
                  << ratio_ten6.first << " +/- " << ratio_ten6.second << "\n";
        std::cout << "[sp19-scale] Recommended scale to match Fa18 Inb = "
                  << scale_fa18.first << " +/- " << scale_fa18.second << "\n";
        std::cout << "[sp19-scale] Recommended scale to match 10.6 GeV combined = "
                  << scale_ten6.first << " +/- " << scale_ten6.second << "\n";
        std::cout << "[sp19-scale] Wrote diagnostics to " << out_dir.string() << "\n";

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[sp19-scale] FATAL: " << e.what() << "\n";
        return false;
    }
}
