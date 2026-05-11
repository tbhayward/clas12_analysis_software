#include "eppi0_normalization.h"

#include "global_cuts.h"

#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TLine.h>
#include <TObject.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

static constexpr double PI = 3.14159265358979323846;
static constexpr double RAD2DEG = 180.0 / PI;
static constexpr double CHARGE_TO_MC_FACTOR = 1.0e-6;
static constexpr double RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC = 1316.875;

static constexpr double RATIO_Y_MIN = 0.0;
static constexpr double RATIO_Y_MAX = 2.0;

static const std::vector<std::string> CSV_PERIOD_ORDER = {
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out"
};

static const std::vector<std::string> WORK_PERIOD_ORDER = {
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb"
};

static const std::vector<std::string> HELICITIES = {
    "unpol",
    "pos",
    "neg"
};

struct ChannelConfig {
    std::string csv_channel;
    std::string cut_prefix;
    std::string title;
};

static ChannelConfig dvcs_config() {
    ChannelConfig cfg;
    cfg.csv_channel = "ep->epg";
    cfg.cut_prefix = "DVCS";
    cfg.title = "ep #rightarrow ep#gamma";
    return cfg;
}

static ChannelConfig eppi0_config() {
    ChannelConfig cfg;
    cfg.csv_channel = "ep->eppi0";
    cfg.cut_prefix = "eppi0";
    cfg.title = "ep #rightarrow ep#pi_{0}";
    return cfg;
}

struct PeriodTags {
    std::string display;
    std::string period_label;
    std::string period_code;
};

struct CubicFit {
    double a[4] = {1.0, 0.0, 0.0, 0.0};
    double ea[4] = {0.0, 0.0, 0.0, 0.0};
    double x_min = 0.0;
    double x_max = 0.0;
    bool valid = false;

    double eval(double x) const {
        return a[0] + x * (a[1] + x * (a[2] + x * a[3]));
    }
};

struct QuarticFit {
    double a[5] = {1.0, 0.0, 0.0, 0.0, 0.0};

    double eval(double x) const {
        return a[0] + x * (a[1] + x * (a[2] + x * (a[3] + x * a[4])));
    }
};

struct RegionNormalization {
    std::string region;
    CubicFit fit;
    double integrated_ratio = 1.0;
    double integrated_ratio_err = 0.0;
};

struct SummaryRatioCurve {
    std::string key;
    std::string title;
    std::string x_title;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ey;
};

struct PeriodNormalization {
    std::string period;
    std::map<std::string, RegionNormalization> regions;
    std::map<std::string, SummaryRatioCurve> summary_ratio_curves;
};

static void fatal(const std::string& msg);

struct AaoPeriodInput {
    double sigma_min_microbarn = 0.0;
    double sigma_max_microbarn = 0.0;
    double sigma_mid_microbarn = 0.0;
    double sigma_half_width_microbarn = 0.0;
    double sigma_relative_uncertainty_from_range = 0.0;
    double generated_events = 0.0;
};

static AaoPeriodInput read_aao_period_input(const std::string& json_path, const std::string& period) {
    std::ifstream in(json_path.c_str());

    if (!in.is_open()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: could not open AAOGEN normalization input JSON: " << json_path;
        fatal(ss.str());
    }

    nlohmann::json j;
    in >> j;

    if (!j.contains("periods") || !j["periods"].is_object()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: AAOGEN normalization input JSON is missing object 'periods': " << json_path;
        fatal(ss.str());
    }

    const nlohmann::json& periods = j["periods"];

    if (!periods.contains(period)) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: AAOGEN normalization input JSON is missing period '" << period << "': " << json_path;
        fatal(ss.str());
    }

    const nlohmann::json& jp = periods[period];

    const std::vector<std::string> required = {
        "sigma_min_microbarn",
        "sigma_max_microbarn",
        "sigma_mid_microbarn",
        "sigma_half_width_microbarn",
        "sigma_relative_uncertainty_from_range",
        "generated_events"
    };

    std::vector<std::string> missing;

    for (const std::string& key : required) {
        if (!jp.contains(key)) {
            missing.push_back(key);
        }
    }

    if (!missing.empty()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: AAOGEN normalization input JSON period '" << period << "' is missing keys:";
        for (const std::string& key : missing) {
            ss << " " << key;
        }
        fatal(ss.str());
    }

    AaoPeriodInput out;
    out.sigma_min_microbarn = jp.at("sigma_min_microbarn").get<double>();
    out.sigma_max_microbarn = jp.at("sigma_max_microbarn").get<double>();
    out.sigma_mid_microbarn = jp.at("sigma_mid_microbarn").get<double>();
    out.sigma_half_width_microbarn = jp.at("sigma_half_width_microbarn").get<double>();
    out.sigma_relative_uncertainty_from_range = jp.at("sigma_relative_uncertainty_from_range").get<double>();
    out.generated_events = jp.at("generated_events").get<double>();

    if (!(out.sigma_mid_microbarn > 0.0)) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: AAOGEN normalization input JSON period '" << period << "' has non-positive sigma_mid_microbarn.";
        fatal(ss.str());
    }

    if (!(out.generated_events > 0.0)) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: AAOGEN normalization input JSON period '" << period << "' has non-positive generated_events.";
        fatal(ss.str());
    }

    return out;
}

static void fatal(const std::string& msg) {
    throw std::runtime_error(msg);
}

static std::string lower_ascii(std::string s) {
    for (char& c : s) {
        c = (char)std::tolower((unsigned char)c);
    }

    return s;
}

static bool has_substr(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

static PeriodTags parse_period_from_key(const std::string& key) {
    const std::string s = lower_ascii(key);
    PeriodTags t;

    if (has_substr(s, "fa18") && has_substr(s, "inb") && has_substr(s, "supp")) {
        t.display = "Fa18 Inb Supp";
        t.period_label = "fa18_inb";
        t.period_code = "Fa18_Inb_Supp";
        return t;
    }

    if (has_substr(s, "fa18") && has_substr(s, "inb")) {
        t.display = "Fa18 Inb";
        t.period_label = "fa18_inb";
        t.period_code = "Fa18_Inb";
        return t;
    }

    if (has_substr(s, "fa18") && has_substr(s, "out")) {
        t.display = "Fa18 Out";
        t.period_label = "fa18_out";
        t.period_code = "Fa18_Out";
        return t;
    }

    if (has_substr(s, "sp19") && has_substr(s, "inb")) {
        t.display = "Sp19 Inb";
        t.period_label = "sp19_inb";
        t.period_code = "Sp19_Inb";
        return t;
    }

    if (has_substr(s, "sp18") && has_substr(s, "inb")) {
        t.display = "Sp18 Inb";
        t.period_label = "sp18_inb";
        t.period_code = "Sp18_Inb";
        return t;
    }

    if (has_substr(s, "sp18") && has_substr(s, "out")) {
        t.display = "Sp18 Out";
        t.period_label = "sp18_out";
        t.period_code = "Sp18_Out";
        return t;
    }

    std::ostringstream ss;
    ss << "[eppi0_norm] FATAL: cannot parse period from tree key: " << key;
    fatal(ss.str());

    return t;
}

static std::string period_dir(const std::string& period) {
    if (period == "Fa18 Inb") return "Fa18_Inb";
    if (period == "Fa18 Out") return "Fa18_Out";
    if (period == "Sp19 Inb") return "Sp19_Inb";
    if (period == "Sp18 Inb") return "Sp18_Inb";
    if (period == "Sp18 Out") return "Sp18_Out";
    if (period == "Fa18 Inb Supp") return "Fa18_Inb_Supp";

    std::ostringstream ss;
    ss << "[eppi0_norm] FATAL: unknown period label: " << period;
    fatal(ss.str());

    return "";
}

static void mkdir_p(const std::string& path) {
    if (!path.empty()) {
        gSystem->mkdir(path.c_str(), true);
    }
}

static double wrap_phi_deg(double v) {
    double p = std::fmod(v, 360.0);

    if (p < 0.0) {
        p += 360.0;
    }

    if (p >= 360.0) {
        p = std::nextafter(360.0, 0.0);
    }

    return p;
}

static bool in_range(double v, double a, double b) {
    return v >= a && v < b;
}

static bool row_accepts_phi(double phi, double pmin, double pmax) {
    if (pmax > pmin) {
        return in_range(phi, pmin, pmax);
    }

    return phi >= pmin || phi < pmax;
}

static std::string topo_dir(int det1, int det2) {
    if (det1 == 1 && det2 == 1) return "FD_FD";
    if (det1 == 2 && det2 == 1) return "CD_FD";
    if (det1 == 2 && det2 == 0) return "CD_FT";

    return "";
}

static std::string topo_label_for_csv(const std::string& topo) {
    if (topo == "FD_FD") return "(FD, FD)";
    if (topo == "CD_FD") return "(CD, FD)";
    if (topo == "CD_FT") return "(CD, FT)";

    return "";
}

// -----------------------------------------------------------------------------
// CSV helpers
// -----------------------------------------------------------------------------

struct CSV {
    std::vector<std::string> header;
    std::unordered_map<std::string, int> index;
    std::vector<std::vector<std::string>> rows;
};

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool inq = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (c == '"') {
            if (inq && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                inq = !inq;
            }
        } else if (c == ',' && !inq) {
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

            for (char ch : s) {
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

static void load_csv_strict(const std::string& path, CSV& csv) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: cannot open CSV: " << path;
        fatal(ss.str());
    }

    std::string line;

    if (!std::getline(fin, line)) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: empty CSV: " << path;
        fatal(ss.str());
    }

    csv.header = split_csv_line(line);
    csv.index.clear();

    for (int i = 0; i < (int)csv.header.size(); ++i) {
        if (csv.index.find(csv.header[i]) != csv.index.end()) {
            std::ostringstream ss;
            ss << "[eppi0_norm] FATAL: duplicate CSV column: " << csv.header[i];
            fatal(ss.str());
        }

        csv.index[csv.header[i]] = i;
    }

    csv.rows.clear();

    while (std::getline(fin, line)) {
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> row = split_csv_line(line);

        if (row.size() < csv.header.size()) {
            row.resize(csv.header.size(), "");
        }

        if (row.size() != csv.header.size()) {
            fatal("[eppi0_norm] FATAL: CSV row width mismatch while loading.");
        }

        csv.rows.push_back(std::move(row));
    }
}

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";

    std::ofstream fout(tmp);

    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: cannot write temp CSV: " << tmp;
        fatal(ss.str());
    }

    fout << join_csv_row(csv.header) << "\n";

    for (const auto& row : csv.rows) {
        fout << join_csv_row(row) << "\n";
    }

    fout.close();

    if (!fout) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: failed writing temp CSV: " << tmp;
        fatal(ss.str());
    }

    (void)std::remove(path.c_str());

    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: failed replacing CSV: " << path;
        fatal(ss.str());
    }
}

static int col_strict(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);

    if (it == csv.index.end()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: missing CSV column: " << name;
        fatal(ss.str());
    }

    return it->second;
}

static double to_double_strict(const std::string& s, const std::string& what) {
    if (s.empty()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: empty numeric cell for " << what;
        fatal(ss.str());
    }

    char* e = nullptr;
    const double v = std::strtod(s.c_str(), &e);

    if (e == s.c_str()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: parse failure for " << what << " value " << s;
        fatal(ss.str());
    }

    return v;
}

static bool valid_cell(const std::string& s) {
    return s == "1" || s == "1.0" || s == "true" || s == "TRUE";
}

static std::string tuple2(double a, double b) {
    std::ostringstream ss;
    ss << std::setprecision(12) << "(" << a << "," << b << ")";
    return ss.str();
}

static std::string tuple4(const CubicFit& p) {
    std::ostringstream ss;
    ss << std::setprecision(12)
       << "(" << p.a[0] << "," << p.a[1] << "," << p.a[2] << "," << p.a[3] << ")";
    return ss.str();
}

static bool parse_tuple_numbers(const std::string& s, std::vector<double>& vals) {
    vals.clear();

    std::string t = s;

    for (char& c : t) {
        if (c == '(' || c == ')' || c == '"') {
            c = ' ';
        }
    }

    std::stringstream ss(t);
    std::string part;

    while (std::getline(ss, part, ',')) {
        char* e = nullptr;
        const double v = std::strtod(part.c_str(), &e);

        if (e == part.c_str()) {
            return false;
        }

        vals.push_back(v);
    }

    return !vals.empty();
}

static double read_current_factor(const CSV& csv,
                                  const ChannelConfig& cfg,
                                  const std::string& sample,
                                  const std::string& period) {
    const std::string name = "current efficiency factor, " + cfg.csv_channel + ", " + sample + ", " + period;
    const int c = col_strict(csv, name);

    if (csv.rows.empty()) {
        fatal("[eppi0_norm] FATAL: CSV has no rows.");
    }

    std::vector<double> vals;

    if (!parse_tuple_numbers(csv.rows.front()[c], vals) || vals.size() < 1) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: cannot parse current-efficiency tuple from column " << name
           << " value '" << csv.rows.front()[c] << "'.";
        fatal(ss.str());
    }

    if (!(std::isfinite(vals[0]) && vals[0] > 0.0)) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: invalid current-efficiency factor in column " << name
           << ": " << vals[0];
        fatal(ss.str());
    }

    return vals[0];
}

static std::vector<std::string> normalization_regions() {
    return {"Sector 1", "Sector 2", "Sector 3", "Sector 4", "Sector 5", "Sector 6", "CD"};
}

static std::string norm_fit_col(const std::string& period, const std::string& region) {
    return "eppi0 cross-section normalization cubic, ep->eppi0, data_over_mc, " + region + ", " + period;
}

static std::string norm_factor_col(const std::string& period, const std::string& region) {
    return "eppi0 cross-section normalization factor, ep->eppi0, data_over_mc, " + region + ", " + period;
}

static std::string normalized_raw_yield_col(const ChannelConfig& cfg,
                                            const std::string& topo_label,
                                            const std::string& period,
                                            const std::string& hel) {
    return "normalized raw yield, " + cfg.csv_channel + ", " + topo_label + ", exp, " + period + ", " + hel;
}

// -----------------------------------------------------------------------------
// Row bins
// -----------------------------------------------------------------------------

struct RowBin {
    double xBmin = 0.0;
    double xBmax = 0.0;
    double Q2min = 0.0;
    double Q2max = 0.0;
    double tmin = 0.0;
    double tmax = 0.0;
    double pmin = 0.0;
    double pmax = 0.0;
    bool valid = false;
};

static std::vector<RowBin> load_rows(const CSV& csv) {
    const int c_xmin = col_strict(csv, "xBmin");
    const int c_xmax = col_strict(csv, "xBmax");
    const int c_qmin = col_strict(csv, "Q2min");
    const int c_qmax = col_strict(csv, "Q2max");
    const int c_tmin = col_strict(csv, "t_abs_min");
    const int c_tmax = col_strict(csv, "t_abs_max");
    const int c_pmin = col_strict(csv, "phimin");
    const int c_pmax = col_strict(csv, "phimax");
    const int c_valid = col_strict(csv, "valid bin");

    std::vector<RowBin> rows;
    rows.reserve(csv.rows.size());

    for (int i = 0; i < (int)csv.rows.size(); ++i) {
        const auto& r = csv.rows[i];

        RowBin b;
        b.xBmin = to_double_strict(r[c_xmin], "xBmin");
        b.xBmax = to_double_strict(r[c_xmax], "xBmax");
        b.Q2min = to_double_strict(r[c_qmin], "Q2min");
        b.Q2max = to_double_strict(r[c_qmax], "Q2max");
        b.tmin = to_double_strict(r[c_tmin], "t_abs_min");
        b.tmax = to_double_strict(r[c_tmax], "t_abs_max");
        b.pmin = to_double_strict(r[c_pmin], "phimin");
        b.pmax = to_double_strict(r[c_pmax], "phimax");
        b.valid = valid_cell(r[c_valid]);

        rows.push_back(b);
    }

    return rows;
}

static bool row_accepts(const RowBin& r, double x, double Q2, double tabs, double phi) {
    if (!r.valid) return false;
    if (!in_range(x, r.xBmin, r.xBmax)) return false;
    if (!in_range(Q2, r.Q2min, r.Q2max)) return false;
    if (!in_range(tabs, r.tmin, r.tmax)) return false;
    if (!row_accepts_phi(phi, r.pmin, r.pmax)) return false;

    return true;
}

// -----------------------------------------------------------------------------
// Charge and sigma cuts
// -----------------------------------------------------------------------------

static std::unordered_map<int, double> read_charge_csv(const std::string& path) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: cannot open charge CSV: " << path;
        fatal(ss.str());
    }

    std::unordered_map<int, double> out;
    std::string line;

    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::vector<std::string> f = split_csv_line(line);

        if (f.size() < 2) {
            continue;
        }

        char* e = nullptr;
        int run = (int)std::strtol(f[0].c_str(), &e, 10);

        if (e == f[0].c_str()) {
            continue;
        }

        e = nullptr;
        double q = std::strtod(f[1].c_str(), &e);

        if (e == f[1].c_str()) {
            continue;
        }

        out[run] = q;
    }

    return out;
}

struct SigmaStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std = std::numeric_limits<double>::quiet_NaN();
};

using CutVarMap = std::unordered_map<std::string, SigmaStats>;
using TopoCutMap = std::unordered_map<std::string, CutVarMap>;

static TopoCutMap load_sigma_cuts(const std::string& path, const std::string& sample_key) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: cannot open combined cuts JSON: " << path;
        fatal(ss.str());
    }

    nlohmann::json j;
    fin >> j;

    TopoCutMap out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const auto& block = it.value();

        if (!block.is_object() || !block.contains(sample_key)) {
            continue;
        }

        const auto& sample = block[sample_key];

        if (!sample.is_object()) {
            continue;
        }

        CutVarMap vm;

        for (auto vit = sample.begin(); vit != sample.end(); ++vit) {
            const auto& obj = vit.value();

            if (!obj.is_object() || !obj.contains("mean") || !obj.contains("std")) {
                continue;
            }

            SigmaStats s;
            s.mean = obj["mean"].get<double>();
            s.std = obj["std"].get<double>();
            vm[vit.key()] = s;
        }

        if (!vm.empty()) {
            out[key] = vm;
        }
    }

    return out;
}

static bool within_3sigma(double v, const SigmaStats& s) {
    if (!std::isfinite(s.mean) || !std::isfinite(s.std) || s.std <= 0.0) {
        return true;
    }

    return std::fabs(v - s.mean) <= 3.0 * s.std;
}

static bool check_sigma(const CutVarMap& vm, const std::string& var, bool has, double val) {
    auto it = vm.find(var);

    if (it == vm.end()) {
        return true;
    }

    if (!has) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: sigma cut requires missing branch " << var;
        fatal(ss.str());
    }

    return within_3sigma(val, it->second);
}

// -----------------------------------------------------------------------------
// Branch binding and cuts
// -----------------------------------------------------------------------------

static std::mutex g_bind_mutex;

struct Branches {
    int runnum = 0; bool has_runnum = false;
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;
    int helicity = 0; bool has_helicity = false;

    double x = 0.0; bool has_x = false;
    double Q2 = 0.0; bool has_Q2 = false;
    double t1 = 0.0; bool has_t1 = false;
    double phi2 = 0.0; bool has_phi2 = false;

    double p1_p = 0.0; bool has_p1_p = false;
    double p1_theta = 0.0; bool has_p1_theta = false;
    double p1_phi = 0.0; bool has_p1_phi = false;
    double p2_p = 0.0; bool has_p2_p = false;
    double p2_theta = 0.0; bool has_p2_theta = false;
    double p2_phi = 0.0; bool has_p2_phi = false;

    double weight = 1.0; bool has_weight = false;

    double open_angle_ep2 = 0.0; bool has_open_angle_ep2 = false;
    double pTmiss = 0.0; bool has_pTmiss = false;
    double Emiss2 = 0.0; bool has_Emiss2 = false;
    double Mx2 = 0.0; bool has_Mx2 = false;
    double Mx2_1 = 0.0; bool has_Mx2_1 = false;
    double Mx2_2 = 0.0; bool has_Mx2_2 = false;
    double xF = 0.0; bool has_xF = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double theta_pi0_pi0 = 0.0; bool has_theta_pi0_pi0 = false;

    double e_p = 0.0; bool has_e_p = false;
    double e_theta = 0.0; bool has_e_theta = false;
    double e_phi = 0.0; bool has_e_phi = false;

    void bind(TTree* t, bool need_weight) {
        std::lock_guard<std::mutex> lock(g_bind_mutex);

        t->SetBranchStatus("*", 0);

        auto ena = [&](const char* name) {
            if (t->GetBranch(name)) {
                t->SetBranchStatus(name, 1);
            }
        };

        ena("runnum");
        ena("detector1");
        ena("detector2");
        ena("helicity");

        ena("x");
        ena("Q2");
        ena("t1");
        ena("phi2");

        ena("p1_p");
        ena("p1_theta");
        ena("p1_phi");
        ena("p2_p");
        ena("p2_theta");
        ena("p2_phi");

        ena("open_angle_ep2");
        ena("pTmiss");
        ena("Emiss2");
        ena("Mx2");
        ena("Mx2_1");
        ena("Mx2_2");
        ena("xF");
        ena("theta_gamma_gamma");
        ena("theta_pi0_pi0");

        ena("e_p");
        ena("e_theta");
        ena("e_phi");

        if (need_weight) {
            ena("weight");
        }

        t->SetCacheSize(0);

        auto bI = [&](const char* name, int* addr, bool& flag) {
            if (t->GetBranch(name)) {
                t->SetBranchAddress(name, addr);
                flag = true;
            }
        };

        auto bD = [&](const char* name, double* addr, bool& flag) {
            if (t->GetBranch(name)) {
                t->SetBranchAddress(name, addr);
                flag = true;
            }
        };

        bI("runnum", &runnum, has_runnum);
        bI("detector1", &detector1, has_detector1);
        bI("detector2", &detector2, has_detector2);
        bI("helicity", &helicity, has_helicity);

        bD("x", &x, has_x);
        bD("Q2", &Q2, has_Q2);
        bD("t1", &t1, has_t1);
        bD("phi2", &phi2, has_phi2);

        bD("p1_p", &p1_p, has_p1_p);
        bD("p1_theta", &p1_theta, has_p1_theta);
        bD("p1_phi", &p1_phi, has_p1_phi);
        bD("p2_p", &p2_p, has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi", &p2_phi, has_p2_phi);

        bD("weight", &weight, has_weight);

        bD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bD("pTmiss", &pTmiss, has_pTmiss);
        bD("Emiss2", &Emiss2, has_Emiss2);
        bD("Mx2", &Mx2, has_Mx2);
        bD("Mx2_1", &Mx2_1, has_Mx2_1);
        bD("Mx2_2", &Mx2_2, has_Mx2_2);
        bD("xF", &xF, has_xF);
        bD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        bD("theta_pi0_pi0", &theta_pi0_pi0, has_theta_pi0_pi0);

        bD("e_p", &e_p, has_e_p);
        bD("e_theta", &e_theta, has_e_theta);
        bD("e_phi", &e_phi, has_e_phi);
    }

    double phi_deg() const {
        return wrap_phi_deg(phi2 * RAD2DEG);
    }

    double t_abs() const {
        return std::fabs(t1);
    }

    double p1_theta_deg() const {
        return p1_theta * RAD2DEG;
    }

    double p1_phi_deg() const {
        return wrap_phi_deg(p1_phi * RAD2DEG);
    }

    double p2_theta_deg() const {
        return p2_theta * RAD2DEG;
    }

    double p2_phi_deg() const {
        return wrap_phi_deg(p2_phi * RAD2DEG);
    }
};

static bool passes_global_dispatch(const Branches& b, const PeriodTags& tags) {
    if (!(b.has_t1 && b.has_open_angle_ep2 && b.has_pTmiss)) {
        return false;
    }

    if (b.has_runnum && is_excluded_run(b.runnum)) {
        return false;
    }

    const GlobalCutConfig& cfg = default_global_cuts();

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_detector1 && b.has_detector2 &&
              b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[eppi0_norm] FATAL: missing branches required by global ycol/topology cuts.");
        }

        return passes_global_cuts(b.t1,
                                  b.open_angle_ep2,
                                  b.pTmiss,
                                  b.detector1,
                                  b.detector2,
                                  tags.period_label,
                                  b.e_p,
                                  b.e_theta,
                                  b.e_phi,
                                  b.p2_p,
                                  b.p2_theta,
                                  b.p2_phi,
                                  cfg);
    }

    if (!(b.has_detector1 && b.has_detector2)) {
        fatal("[eppi0_norm] FATAL: missing detector1/detector2.");
    }

    return passes_global_cuts(b.t1,
                              b.open_angle_ep2,
                              b.pTmiss,
                              b.detector1,
                              b.detector2,
                              cfg);
}

static bool passes_sigma_dispatch(const ChannelConfig& cfg,
                                  const PeriodTags& tags,
                                  const TopoCutMap& cuts,
                                  const Branches& b) {
    if (!(b.has_detector1 && b.has_detector2)) {
        return false;
    }

    const std::string topo = topo_dir(b.detector1, b.detector2);

    if (topo.empty()) {
        return false;
    }

    const std::string key = cfg.cut_prefix + "_" + tags.period_code + "_" + topo;

    auto it = cuts.find(key);

    if (it == cuts.end()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: missing sigma-cut key: " << key;
        fatal(ss.str());
    }

    const CutVarMap& vm = it->second;

    if (!check_sigma(vm, "Emiss2", b.has_Emiss2, b.Emiss2)) return false;
    if (!check_sigma(vm, "Mx2", b.has_Mx2, b.Mx2)) return false;
    if (!check_sigma(vm, "Mx2_1", b.has_Mx2_1, b.Mx2_1)) return false;
    if (!check_sigma(vm, "Mx2_2", b.has_Mx2_2, b.Mx2_2)) return false;
    if (!check_sigma(vm, "pTmiss", b.has_pTmiss, b.pTmiss)) return false;
    if (!check_sigma(vm, "xF", b.has_xF, b.xF)) return false;

    if (cfg.csv_channel == "ep->eppi0") {
        if (!check_sigma(vm, "theta_pi0_pi0", b.has_theta_pi0_pi0, b.theta_pi0_pi0)) return false;
    } else {
        if (!check_sigma(vm, "theta_gamma_gamma", b.has_theta_gamma_gamma, b.theta_gamma_gamma)) return false;
    }

    return true;
}

static bool passes_event_selection(const ChannelConfig& cfg,
                                   const PeriodTags& tags,
                                   const TopoCutMap& cuts,
                                   const Branches& b) {
    if (!passes_global_dispatch(b, tags)) {
        return false;
    }

    if (!passes_sigma_dispatch(cfg, tags, cuts, b)) {
        return false;
    }

    return true;
}

// -----------------------------------------------------------------------------
// Histogram configuration
// -----------------------------------------------------------------------------

struct VarConfig {
    std::string key;
    std::string title;
    std::string x_title;
    int particle = 1;
    int kind = 0; // 0 theta, 1 momentum, 2 phi, 3 kinematics
    int n_panels = 7;
};

static std::vector<VarConfig> variable_configs() {
    return {
        {"p1_theta", "p_{1} #theta", "p_{1} #theta (deg)", 1, 0, 7},
        {"p1_p", "p_{1} momentum", "p_{1} momentum (GeV)", 1, 1, 7},
        {"p1_phi", "p_{1} #phi", "p_{1} #phi (deg)", 1, 2, 2},

        {"p2_theta", "p_{2} #theta", "p_{2} #theta (deg)", 2, 0, 7},
        {"p2_p", "p_{2} momentum", "p_{2} momentum (GeV)", 2, 1, 7},
        {"p2_phi", "p_{2} #phi", "p_{2} #phi (deg)", 2, 2, 2},

        {"x", "x_{B}", "x_{B}", 0, 3, 1},
        {"Q2", "Q^{2}", "Q^{2} (GeV^{2})", 0, 3, 1},
        {"minus_t1", "-t", "-t (GeV^{2})", 0, 3, 1},
        {"Mx2", "M_{X}^{2}", "M_{X}^{2} (GeV^{2})", 0, 3, 1},
        {"Mx2_1", "M_{X}^{2}(ep)", "M_{X}^{2}(ep) (GeV^{2})", 0, 3, 1},
        {"Mx2_2", "M_{X}^{2}(e#gamma)", "M_{X}^{2}(e#gamma) (GeV^{2})", 0, 3, 1},
        {"pTmiss", "p_{T}^{miss}", "p_{T}^{miss} (GeV)", 0, 3, 1}
    };
}

static int sector_from_phi(double phi) {
    const double p = wrap_phi_deg(phi);

    if ((p >= 330.0 && p < 360.0) || (p >= 0.0 && p < 30.0)) return 1;
    if (p >= 30.0 && p < 90.0) return 2;
    if (p >= 90.0 && p < 150.0) return 3;
    if (p >= 150.0 && p < 210.0) return 4;
    if (p >= 210.0 && p < 270.0) return 5;
    if (p >= 270.0 && p < 330.0) return 6;

    return 0;
}

static int panel_index(const VarConfig& v, const Branches& b) {
    if (v.kind == 3) {
        return 0;
    }

    double theta = (v.particle == 1) ? b.p1_theta_deg() : b.p2_theta_deg();
    double phi = (v.particle == 1) ? b.p1_phi_deg() : b.p2_phi_deg();

    if (v.kind == 2) {
        if (v.particle == 1) {
            if (theta >= 0.0 && theta < 40.0) return 0;
            if (theta >= 40.0 && theta < 70.0) return 1;
        } else {
            if (theta >= 5.0 && theta < 40.0) return 0;
            if (theta >= 0.0 && theta < 5.0) return 1;
        }

        return -1;
    }

    if (v.particle == 1) {
        if (theta >= 0.0 && theta < 40.0) {
            int sec = sector_from_phi(phi);
            return (sec >= 1 && sec <= 6) ? sec - 1 : -1;
        }

        if (theta >= 40.0 && theta < 70.0) return 6;

        return -1;
    }

    if (theta >= 5.0 && theta < 40.0) {
        int sec = sector_from_phi(phi);
        return (sec >= 1 && sec <= 6) ? sec - 1 : -1;
    }

    if (theta >= 0.0 && theta < 5.0) return 6;

    return -1;
}

static double value_for_var(const VarConfig& v, const Branches& b) {
    if (v.kind == 3) {
        if (v.key == "x") return b.x;
        if (v.key == "Q2") return b.Q2;
        if (v.key == "minus_t1") return -b.t1;
        if (v.key == "Mx2") return b.Mx2;
        if (v.key == "Mx2_1") return b.Mx2_1;
        if (v.key == "Mx2_2") return b.Mx2_2;
        if (v.key == "pTmiss") return b.pTmiss;

        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: unknown kinematic variable key: " << v.key;
        fatal(ss.str());
    }

    if (v.particle == 1) {
        if (v.kind == 0) return b.p1_theta_deg();
        if (v.kind == 1) return b.p1_p;
        return b.p1_phi_deg();
    }

    if (v.kind == 0) return b.p2_theta_deg();
    if (v.kind == 1) return b.p2_p;

    return b.p2_phi_deg();
}

static int nbins_for_panel(const VarConfig& v, int panel) {
    if (v.kind == 3) {
        if (v.key == "x") return 35;
        if (v.key == "Q2") return 36;
        if (v.key == "minus_t1") return 40;
        if (v.key == "Mx2") return 60;
        if (v.key == "Mx2_1") return 60;
        if (v.key == "Mx2_2") return 60;
        if (v.key == "pTmiss") return 60;

        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: unknown kinematic variable key: " << v.key;
        fatal(ss.str());
    }

    if (v.kind == 0) return (panel < 6) ? 12 : 23;
    if (v.kind == 1) return (panel < 6) ? 12 : 23;

    return (panel == 0) ? 12 : 24;
}

static void range_for_panel(const VarConfig& v, int panel, double& xmin, double& xmax) {
    if (v.kind == 3) {
        if (v.key == "x") {
            xmin = 0.0;
            xmax = 0.7;
            return;
        }

        if (v.key == "Q2") {
            xmin = 1.0;
            xmax = 7.0;
            return;
        }

        if (v.key == "minus_t1") {
            xmin = 0.0;
            xmax = 1.0;
            return;
        }

        if (v.key == "Mx2") {
            xmin = -0.03;
            xmax = 0.03;
            return;
        }

        if (v.key == "Mx2_1") {
            xmin = -1.5;
            xmax = 1.5;
            return;
        }

        if (v.key == "Mx2_2") {
            xmin = 0.0;
            xmax = 3.0;
            return;
        }

        if (v.key == "pTmiss") {
            xmin = 0.0;
            xmax = 0.3;
            return;
        }

        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: unknown kinematic variable key: " << v.key;
        fatal(ss.str());
    }

    if (v.kind == 0) {
        if (v.particle == 1) {
            xmin = (panel < 6) ? 0.0 : 40.0;
            xmax = (panel < 6) ? 40.0 : 70.0;
        } else {
            xmin = (panel < 6) ? 5.0 : 0.0;
            xmax = (panel < 6) ? 40.0 : 5.0;
        }

        return;
    }

    if (v.kind == 1) {
        xmin = (v.particle == 1) ? 0.3 : 2.0;
        xmax = (v.particle == 1) ? 1.3 : 6.0;
        return;
    }

    xmin = 0.0;
    xmax = 360.0;
}

static std::string panel_label(const VarConfig& v, int panel) {
    if (v.kind == 3) {
        return v.title;
    }

    if (v.kind == 2) {
        if (v.particle == 1) return (panel == 0) ? "FD" : "CD";
        return (panel == 0) ? "FD" : "FT";
    }

    if (panel < 6) {
        std::ostringstream ss;
        ss << "Sector " << panel + 1;
        return ss.str();
    }

    if (v.particle == 1) return "CD";

    return "FT";
}

static std::string normalization_region_for_event(const Branches& b) {
    const double theta = b.p1_theta_deg();

    if (theta >= 0.0 && theta < 40.0) {
        const int sec = sector_from_phi(b.p1_phi_deg());

        if (sec >= 1 && sec <= 6) {
            std::ostringstream ss;
            ss << "Sector " << sec;
            return ss.str();
        }

        return "";
    }

    if (theta >= 40.0 && theta < 70.0) {
        return "CD";
    }

    return "";
}

struct HistSet {
    std::map<std::string, std::vector<TH1D*>> hists;
};

static HistSet make_hist_set(const std::string& prefix) {
    HistSet hs;

    for (const VarConfig& v : variable_configs()) {
        std::vector<TH1D*> vec;

        for (int p = 0; p < v.n_panels; ++p) {
            double xmin = 0.0;
            double xmax = 1.0;
            range_for_panel(v, p, xmin, xmax);

            const int nb = nbins_for_panel(v, p);

            std::ostringstream name;
            name << prefix << "_" << v.key << "_panel_" << p;

            TH1D* h = new TH1D(name.str().c_str(), "", nb, xmin, xmax);
            h->Sumw2();
            h->GetXaxis()->SetTitle(v.x_title.c_str());
            h->GetYaxis()->SetTitle("Counts");

            vec.push_back(h);
        }

        hs.hists[v.key] = vec;
    }

    return hs;
}

static void delete_hist_set(HistSet& hs) {
    for (auto& kv : hs.hists) {
        for (TH1D* h : kv.second) {
            delete h;
        }

        kv.second.clear();
    }

    hs.hists.clear();
}

static void fill_hist_set(HistSet& hs, const Branches& b, double w) {
    for (const VarConfig& v : variable_configs()) {
        const int p = panel_index(v, b);

        if (p < 0 || p >= v.n_panels) {
            continue;
        }

        const double value = value_for_var(v, b);
        hs.hists[v.key][p]->Fill(value, w);
    }
}

static void add_hist_set(HistSet& dest, const HistSet& src) {
    for (const auto& kv : src.hists) {
        auto it = dest.hists.find(kv.first);

        if (it == dest.hists.end()) {
            continue;
        }

        for (size_t i = 0; i < kv.second.size(); ++i) {
            it->second[i]->Add(kv.second[i]);
        }
    }
}

static void add_bin_to_compatible_hist(TH1D& target,
                                       double x,
                                       double content,
                                       double error) {
    if (!(content != 0.0 || error != 0.0)) {
        return;
    }

    const int b = target.FindBin(x);

    if (b < 1 || b > target.GetNbinsX()) {
        return;
    }

    const double old_content = target.GetBinContent(b);
    const double old_error = target.GetBinError(b);

    target.SetBinContent(b, old_content + content);
    target.SetBinError(b, std::sqrt(old_error * old_error + error * error));
}

static void add_hist_to_compatible_hist(TH1D& target, const TH1D* source) {
    if (!source) {
        return;
    }

    for (int b = 1; b <= source->GetNbinsX(); ++b) {
        const double x = source->GetBinCenter(b);
        const double content = source->GetBinContent(b);
        const double error = source->GetBinError(b);

        add_bin_to_compatible_hist(target, x, content, error);
    }
}

static void build_p1_theta_fit_histograms_by_region(const HistSet& data_hists,
                                                    const HistSet& mc_hists,
                                                    const std::string& period,
                                                    std::map<std::string, TH1D*>& h_data_fit,
                                                    std::map<std::string, TH1D*>& h_mc_fit) {
    auto itd = data_hists.hists.find("p1_theta");
    auto itm = mc_hists.hists.find("p1_theta");

    if (itd == data_hists.hists.end() || itm == mc_hists.hists.end()) {
        fatal("[eppi0_norm] FATAL: p1_theta histograms are missing while building regional fit histograms.");
    }

    if (itd->second.size() != 7 || itm->second.size() != 7) {
        fatal("[eppi0_norm] FATAL: p1_theta regional histogram count is not 7.");
    }

    const std::vector<std::string> regions = normalization_regions();

    for (int panel = 0; panel < 7; ++panel) {
        const std::string& region = regions[panel];

        h_data_fit[region] = (TH1D*)itd->second[panel]->Clone(
            ("h_data_fit_" + period_dir(period) + "_" + std::to_string(panel)).c_str()
        );

        h_mc_fit[region] = (TH1D*)itm->second[panel]->Clone(
            ("h_mc_fit_" + period_dir(period) + "_" + std::to_string(panel)).c_str()
        );

        h_data_fit[region]->SetDirectory(nullptr);
        h_mc_fit[region]->SetDirectory(nullptr);

        std::cout << "[eppi0_norm] " << period << " " << region
                  << " p1_theta fit histogram integrals: data=" << h_data_fit[region]->Integral()
                  << " mc=" << h_mc_fit[region]->Integral() << std::endl;
    }
}

static CubicFit fit_p1_theta_ratio_region(const std::string& period,
                                          const std::string& region,
                                          const std::string& outdir,
                                          TH1D* h_data,
                                          TH1D* h_mc) {
    CubicFit p;

    if (!h_data || !h_mc) {
        return p;
    }

    TGraphErrors gr;
    int ip = 0;

    double fit_xmin = std::numeric_limits<double>::infinity();
    double fit_xmax = -std::numeric_limits<double>::infinity();

    for (int b = 1; b <= h_data->GetNbinsX(); ++b) {
        const double d = h_data->GetBinContent(b);
        const double m = h_mc->GetBinContent(b);

        if (!(d > 0.0 && m > 0.0)) {
            continue;
        }

        const double x = h_data->GetBinCenter(b);
        const double ed = h_data->GetBinError(b);
        const double em = h_mc->GetBinError(b);
        const double r = d / m;

        double er = 0.0;

        if (ed > 0.0 || em > 0.0) {
            const double rd = (ed > 0.0) ? ed / d : 0.0;
            const double rm = (em > 0.0) ? em / m : 0.0;
            er = std::fabs(r) * std::sqrt(rd * rd + rm * rm);
        }

        if (!(er > 0.0)) {
            er = 1.0;
        }

        gr.SetPoint(ip, x, r);
        gr.SetPointError(ip, 0.0, er);
        ++ip;

        const double low = h_data->GetBinLowEdge(b);
        const double high = h_data->GetBinLowEdge(b + 1);

        fit_xmin = std::min(fit_xmin, low);
        fit_xmax = std::max(fit_xmax, high);
    }

    if (ip < 4 || !(fit_xmax > fit_xmin)) {
        p.a[0] = 1.0;
        p.x_min = (region == "CD") ? 40.0 : 0.0;
        p.x_max = (region == "CD") ? 70.0 : 40.0;
        p.valid = true;

        std::cout << "[eppi0_norm] WARNING: " << period << " " << region
                  << " has fewer than four valid data/MC bins. Using unity cubic fit."
                  << std::endl;

        return p;
    }

    TF1 f("f_ratio_pol3", "pol3", fit_xmin, fit_xmax);
    TFitResultPtr fit_result = gr.Fit(&f, "SQ0", "", fit_xmin, fit_xmax);
    (void)fit_result;

    for (int i = 0; i < 4; ++i) {
        p.a[i] = f.GetParameter(i);
        p.ea[i] = f.GetParError(i);
    }

    p.x_min = fit_xmin;
    p.x_max = fit_xmax;
    p.valid = true;

    const std::string fit_outdir = outdir + "/p1/fits";
    mkdir_p(fit_outdir);

    std::string safe_region = region;

    for (char& c : safe_region) {
        if (c == ' ') {
            c = '_';
        }
    }

    TCanvas c(("c_p1_theta_ratio_fit_" + safe_region).c_str(), "", 1000, 750);
    c.SetGrid(1, 1);
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.08);

    double xframe_min = 0.0;
    double xframe_max = 40.0;

    if (region == "CD") {
        xframe_min = 40.0;
        xframe_max = 70.0;
    }

    TH1D* frame = (TH1D*)gPad->DrawFrame(xframe_min, RATIO_Y_MIN, xframe_max, RATIO_Y_MAX);
    frame->SetTitle((period + "  " + region + "  ep #rightarrow ep#pi_{0}  p_{1} #theta ratio cubic fit").c_str());
    frame->GetXaxis()->SetTitle("p_{1} #theta (deg)");
    frame->GetYaxis()->SetTitle("data / MC");
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);

    TLine* unity_line = new TLine(xframe_min, 1.0, xframe_max, 1.0);
    unity_line->SetLineStyle(2);
    unity_line->SetLineColor(kRed + 1);
    unity_line->SetLineWidth(2);
    unity_line->SetBit(TObject::kCanDelete);
    unity_line->Draw("SAME");

    gr.SetMarkerStyle(20);
    gr.SetMarkerColor(kBlack);
    gr.SetLineColor(kBlack);
    gr.Draw("PE SAME");

    TF1 fdraw("f_ratio_pol3_draw", "pol3", fit_xmin, fit_xmax);

    for (int i = 0; i < 4; ++i) {
        fdraw.SetParameter(i, p.a[i]);
    }

    fdraw.SetLineColor(kRed + 1);
    fdraw.SetLineWidth(2);
    fdraw.Draw("L SAME");

    TLegend leg(0.56, 0.66, 0.90, 0.84);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(1001);
    leg.SetBorderSize(1);
    leg.AddEntry(&gr, "data / MC", "pe");
    leg.AddEntry(&fdraw, "cubic fit", "l");
    leg.Draw();

    c.SaveAs((fit_outdir + "/p1_theta_ratio_cubic_fit_" + safe_region + ".png").c_str());

    return p;
}

static std::vector<TGraphErrors*> make_ratio_graphs(const std::vector<TH1D*>& hd,
                                                    const std::vector<TH1D*>& hm) {
    std::vector<TGraphErrors*> out;

    for (size_t i = 0; i < hd.size(); ++i) {
        TGraphErrors* g = new TGraphErrors();
        int ip = 0;

        for (int b = 1; b <= hd[i]->GetNbinsX(); ++b) {
            const double d = hd[i]->GetBinContent(b);
            const double m = hm[i]->GetBinContent(b);

            if (!(d > 0.0 && m > 0.0)) {
                continue;
            }

            const double ed = hd[i]->GetBinError(b);
            const double em = hm[i]->GetBinError(b);
            const double r = d / m;

            const double er = std::fabs(r) * std::sqrt((ed / d) * (ed / d) + (em / m) * (em / m));

            g->SetPoint(ip, hd[i]->GetBinCenter(b), r);
            g->SetPointError(ip, 0.0, er);
            ++ip;
        }

        g->SetMarkerStyle(20);
        g->SetMarkerColor(kBlack);
        g->SetLineColor(kBlack);

        out.push_back(g);
    }

    return out;
}

static SummaryRatioCurve make_summary_ratio_curve(const std::string& period,
                                                  const VarConfig& v,
                                                  const std::vector<TH1D*>& hd,
                                                  const std::vector<TH1D*>& hm) {
    SummaryRatioCurve curve;
    curve.key = v.key;
    curve.title = v.title;
    curve.x_title = v.x_title;

    if (hd.empty() || hm.empty() || !hd[0] || !hm[0]) {
        return curve;
    }

    TH1D* h_data = hd[0];
    TH1D* h_mc = hm[0];

    const int n_bins = std::min(h_data->GetNbinsX(), h_mc->GetNbinsX());

    for (int b = 1; b <= n_bins; ++b) {
        const double d = h_data->GetBinContent(b);
        const double m = h_mc->GetBinContent(b);

        if (!(d > 0.0 && m > 0.0)) {
            continue;
        }

        const double ed = h_data->GetBinError(b);
        const double em = h_mc->GetBinError(b);
        const double ratio = d / m;

        const double rel_d = (ed > 0.0) ? ed / d : 0.0;
        const double rel_m = (em > 0.0) ? em / m : 0.0;
        const double ratio_err = std::fabs(ratio) * std::sqrt(rel_d * rel_d + rel_m * rel_m);

        curve.x.push_back(h_data->GetBinCenter(b));
        curve.y.push_back(ratio);
        curve.ey.push_back(ratio_err);
    }

    std::cout << "[eppi0_norm] " << period
              << " summary ratio curve " << v.key
              << " points=" << curve.x.size() << std::endl;

    return curve;
}

struct PlotNormInfo {
    double total_charge_mc = 0.0;
    double integrated_luminosity_pb_inv = 0.0;
    double n_gen = 0.0;
    double event_norm = 0.0;
};

static void style_hist_for_plot(TH1D* h, int color, int marker) {
    if (!h) {
        return;
    }

    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.8);
    h->SetStats(false);
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.45);
}

static void style_ratio_frame(TH1D* frame) {
    if (!frame) {
        return;
    }

    frame->SetStats(false);
    frame->SetMinimum(RATIO_Y_MIN);
    frame->SetMaximum(RATIO_Y_MAX);
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetTitleSize(0.050);
    frame->GetYaxis()->SetTitleSize(0.050);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.45);
}

static void draw_panel_title(const std::string& text) {
    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(42);
    latex.SetTextSize(0.050);
    latex.DrawLatex(0.18, 0.925, text.c_str());
}

static void draw_normalization_pad(const std::string& period,
                                   const VarConfig& v,
                                   const PlotNormInfo& norm,
                                   bool ratio_mode) {
    TPad* pad = (TPad*)gPad;

    if (!pad) {
        return;
    }

    pad->Clear();
    pad->SetFillColor(kWhite);
    pad->SetLeftMargin(0.08);
    pad->SetRightMargin(0.04);
    pad->SetTopMargin(0.04);
    pad->SetBottomMargin(0.04);

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(42);
    latex.SetTextSize((v.kind == 2) ? 0.040 : 0.034);

    const double x0 = 0.10;

    latex.DrawLatex(x0, 0.90, period.c_str());

    if (v.kind == 3) {
        latex.DrawLatex(x0, 0.80, "kinematics");
    } else {
        latex.DrawLatex(x0, 0.80, (v.particle == 1) ? "p_{1}" : "p_{2}");
    }

    latex.DrawLatex(x0, 0.70, v.title.c_str());
    latex.DrawLatex(x0, 0.60, ratio_mode ? "Ratio normalization" : "MC normalization");

    std::ostringstream ss;

    ss << std::setprecision(6) << "Q = " << norm.total_charge_mc << " mC";
    latex.DrawLatex(x0, 0.49, ss.str().c_str());

    ss.str("");
    ss.clear();
    ss << std::setprecision(6) << "L_{int} = " << norm.integrated_luminosity_pb_inv << " pb^{-1}";
    latex.DrawLatex(x0, 0.39, ss.str().c_str());

    ss.str("");
    ss.clear();
    ss << std::setprecision(6) << "N_{gen} = " << norm.n_gen;
    latex.DrawLatex(x0, 0.29, ss.str().c_str());

    ss.str("");
    ss.clear();
    ss << std::setprecision(6) << "event weight = " << norm.event_norm;
    latex.DrawLatex(x0, 0.19, ss.str().c_str());

    latex.DrawLatex(x0, 0.09, ratio_mode ? "ratio: data / MC" : "MC fill: const AAOgen norm");
}

static TLegend* make_comparison_legend(TH1D* h_data, TH1D* h_mc) {
    TLegend* leg = new TLegend(0.56, 0.705, 0.90, 0.875);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.032);
    leg->SetMargin(0.22);
    leg->AddEntry(h_data, "data", "l");
    leg->AddEntry(h_mc, "AAOgen MC", "l");

    return leg;
}

static std::pair<double, double> common_y_range(const std::vector<TH1D*>& hd,
                                                const std::vector<TH1D*>& hm,
                                                const std::vector<int>& panels,
                                                bool /*log_y*/) {
    double ymax = 0.0;

    for (int p : panels) {
        if (p < 0 || p >= (int)hd.size() || p >= (int)hm.size()) {
            continue;
        }

        ymax = std::max(ymax, hd[p]->GetMaximum());
        ymax = std::max(ymax, hm[p]->GetMaximum());
    }

    if (!(ymax > 0.0)) {
        ymax = 1.0;
    }

    return {0.0, 1.20 * ymax};
}

static std::pair<double, double> y_range_for_panel(const VarConfig& v,
                                                   const std::vector<TH1D*>& hd,
                                                   const std::vector<TH1D*>& hm,
                                                   int panel,
                                                   bool log_y) {
    if (v.kind == 3) {
        return common_y_range(hd, hm, {0}, log_y);
    }

    if (v.kind == 2) {
        if (panel == 0) return common_y_range(hd, hm, {0}, log_y);
        if (panel == 1) return common_y_range(hd, hm, {1}, log_y);
    } else {
        if (panel >= 0 && panel <= 5) return common_y_range(hd, hm, {0, 1, 2, 3, 4, 5}, log_y);
        if (panel == 6) return common_y_range(hd, hm, {6}, log_y);
    }

    return {0.0, 1.0};
}

static void draw_comparison_panel(const std::string& /*period*/,
                                  const VarConfig& v,
                                  TH1D* h_data,
                                  TH1D* h_mc,
                                  int panel,
                                  const std::pair<double, double>& yr,
                                  bool /*log_y*/) {
    TPad* pad = (TPad*)gPad;

    pad->SetLeftMargin(0.16);
    pad->SetRightMargin(0.06);
    pad->SetTopMargin(0.12);
    pad->SetBottomMargin(0.14);
    pad->SetLogy(false);
    pad->SetGrid(0, 0);

    style_hist_for_plot(h_data, kBlue + 1, 20);
    style_hist_for_plot(h_mc, kRed + 1, 24);

    double xmin = 0.0;
    double xmax = 1.0;
    range_for_panel(v, panel, xmin, xmax);

    TH1D* frame = (TH1D*)pad->DrawFrame(xmin, yr.first, xmax, yr.second);
    frame->SetStats(false);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(v.x_title.c_str());
    frame->GetYaxis()->SetTitle("Counts");
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetTitleSize(0.050);
    frame->GetYaxis()->SetTitleSize(0.050);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.45);

    h_data->Draw("HIST SAME");
    h_mc->Draw("HIST SAME");

    draw_panel_title(panel_label(v, panel));

    TLegend* leg = make_comparison_legend(h_data, h_mc);
    leg->Draw();

    pad->RedrawAxis();
}

static void draw_ratio_panel(const std::string& /*period*/,
                             const VarConfig& v,
                             TGraphErrors* gr,
                             int panel,
                             const std::map<std::string, CubicFit>* fits_for_p1_theta) {
    TPad* pad = (TPad*)gPad;

    pad->SetLeftMargin(0.16);
    pad->SetRightMargin(0.06);
    pad->SetTopMargin(0.12);
    pad->SetBottomMargin(0.14);
    pad->SetGrid(0, 0);

    double xmin = 0.0;
    double xmax = 1.0;
    range_for_panel(v, panel, xmin, xmax);

    TH1D* frame = (TH1D*)pad->DrawFrame(xmin, RATIO_Y_MIN, xmax, RATIO_Y_MAX);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(v.x_title.c_str());
    frame->GetYaxis()->SetTitle("data / MC");
    style_ratio_frame(frame);

    TLine* unity_line = new TLine(xmin, 1.0, xmax, 1.0);
    unity_line->SetLineStyle(2);
    unity_line->SetLineColor(kRed + 1);
    unity_line->SetLineWidth(2);
    unity_line->SetBit(TObject::kCanDelete);
    unity_line->Draw("SAME");

    if (gr) {
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.8);
        gr->SetMarkerColor(kBlack);
        gr->SetLineColor(kBlack);
        gr->Draw("PE SAME");
    }

    if (fits_for_p1_theta && v.key == "p1_theta") {
        const std::string region = panel_label(v, panel);
        auto it_fit = fits_for_p1_theta->find(region);

        if (it_fit != fits_for_p1_theta->end() && it_fit->second.valid) {
            const CubicFit& fit = it_fit->second;
            const double draw_xmin = std::max(xmin, fit.x_min);
            const double draw_xmax = std::min(xmax, fit.x_max);

            if (draw_xmax > draw_xmin) {
                TF1 f("f_cubic_overlay", "pol3", draw_xmin, draw_xmax);

                for (int i = 0; i < 4; ++i) {
                    f.SetParameter(i, fit.a[i]);
                }

                f.SetLineColor(kRed + 1);
                f.SetLineWidth(2);
                f.DrawCopy("L SAME");
            }
        }
    }

    draw_panel_title(panel_label(v, panel));
    pad->RedrawAxis();
}

static void draw_sector_comparison_canvas(const std::string& out_png,
                                          const std::string& period,
                                          const VarConfig& v,
                                          const std::vector<TH1D*>& hd,
                                          const std::vector<TH1D*>& hm,
                                          const PlotNormInfo& norm,
                                          bool log_y) {
    TCanvas canvas(("canvas_" + v.key + (log_y ? "_comparison_log" : "_comparison_linear")).c_str(), "", 1600, 900);
    canvas.Divide(4, 2);

    const int canvas_pad_for_panel[7] = {1, 2, 3, 5, 6, 7, 8};

    for (int panel = 0; panel < 7; ++panel) {
        canvas.cd(canvas_pad_for_panel[panel]);

        const auto yr = y_range_for_panel(v, hd, hm, panel, log_y);

        draw_comparison_panel(period,
                              v,
                              hd[panel],
                              hm[panel],
                              panel,
                              yr,
                              log_y);
    }

    canvas.cd(4);
    draw_normalization_pad(period, v, norm, false);

    canvas.SaveAs(out_png.c_str());
}

static void draw_phi_comparison_canvas(const std::string& out_png,
                                       const std::string& period,
                                       const VarConfig& v,
                                       const std::vector<TH1D*>& hd,
                                       const std::vector<TH1D*>& hm,
                                       const PlotNormInfo& norm,
                                       bool log_y) {
    TCanvas canvas(("canvas_" + v.key + (log_y ? "_comparison_log" : "_comparison_linear")).c_str(), "", 1600, 500);
    canvas.Divide(3, 1);

    for (int panel = 0; panel < 2; ++panel) {
        canvas.cd(panel + 1);

        const auto yr = y_range_for_panel(v, hd, hm, panel, log_y);

        draw_comparison_panel(period,
                              v,
                              hd[panel],
                              hm[panel],
                              panel,
                              yr,
                              log_y);
    }

    canvas.cd(3);
    draw_normalization_pad(period, v, norm, false);

    canvas.SaveAs(out_png.c_str());
}

static void draw_kinematics_comparison_canvas(const std::string& out_png,
                                              const std::string& period,
                                              const VarConfig& v,
                                              const std::vector<TH1D*>& hd,
                                              const std::vector<TH1D*>& hm,
                                              const PlotNormInfo& norm,
                                              bool log_y) {
    TCanvas canvas(("canvas_" + v.key + (log_y ? "_comparison_log" : "_comparison_linear")).c_str(), "", 1100, 750);
    canvas.Divide(2, 1);

    canvas.cd(1);

    const auto yr = y_range_for_panel(v, hd, hm, 0, log_y);

    draw_comparison_panel(period,
                          v,
                          hd[0],
                          hm[0],
                          0,
                          yr,
                          log_y);

    canvas.cd(2);
    draw_normalization_pad(period, v, norm, false);

    canvas.SaveAs(out_png.c_str());
}

static void draw_sector_ratio_canvas(const std::string& out_png,
                                     const std::string& period,
                                     const VarConfig& v,
                                     const std::vector<TGraphErrors*>& gr,
                                     const PlotNormInfo& norm,
                                     const std::map<std::string, CubicFit>* fits_for_p1_theta) {
    TCanvas canvas(("canvas_" + v.key + "_ratio_linear").c_str(), "", 1600, 900);
    canvas.Divide(4, 2);

    const int canvas_pad_for_panel[7] = {1, 2, 3, 5, 6, 7, 8};

    for (int panel = 0; panel < 7; ++panel) {
        canvas.cd(canvas_pad_for_panel[panel]);

        draw_ratio_panel(period,
                         v,
                         gr[panel],
                         panel,
                         fits_for_p1_theta);
    }

    canvas.cd(4);
    draw_normalization_pad(period, v, norm, true);

    canvas.SaveAs(out_png.c_str());
}

static void draw_phi_ratio_canvas(const std::string& out_png,
                                  const std::string& period,
                                  const VarConfig& v,
                                  const std::vector<TGraphErrors*>& gr,
                                  const PlotNormInfo& norm,
                                  const std::map<std::string, CubicFit>* fits_for_p1_theta) {
    TCanvas canvas(("canvas_" + v.key + "_ratio_linear").c_str(), "", 1600, 500);
    canvas.Divide(3, 1);

    for (int panel = 0; panel < 2; ++panel) {
        canvas.cd(panel + 1);

        draw_ratio_panel(period,
                         v,
                         gr[panel],
                         panel,
                         fits_for_p1_theta);
    }

    canvas.cd(3);
    draw_normalization_pad(period, v, norm, true);

    canvas.SaveAs(out_png.c_str());
}

static void draw_kinematics_ratio_canvas(const std::string& out_png,
                                         const std::string& period,
                                         const VarConfig& v,
                                         const std::vector<TGraphErrors*>& gr,
                                         const PlotNormInfo& norm,
                                         const std::map<std::string, CubicFit>* fits_for_p1_theta) {
    TCanvas canvas(("canvas_" + v.key + "_ratio_linear").c_str(), "", 1100, 750);
    canvas.Divide(2, 1);

    canvas.cd(1);

    draw_ratio_panel(period,
                     v,
                     gr[0],
                     0,
                     fits_for_p1_theta);

    canvas.cd(2);
    draw_normalization_pad(period, v, norm, true);

    canvas.SaveAs(out_png.c_str());
}

static void plot_variable(const std::string& outdir,
                          const std::string& period,
                          const VarConfig& v,
                          const std::vector<TH1D*>& hd,
                          const std::vector<TH1D*>& hm,
                          const PlotNormInfo& norm,
                          const std::map<std::string, CubicFit>* fits_for_p1_theta) {
    mkdir_p(outdir);

    std::string subdir;

    if (v.kind == 3) {
        subdir = outdir + "/kinematics";
    } else {
        subdir = outdir + ((v.particle == 1) ? "/p1" : "/p2");
    }

    mkdir_p(subdir);

    const std::string linear_png = subdir + "/" + v.key + "_comparison_linear.png";
    const std::string ratio_png = subdir + "/" + v.key + "_ratio_linear.png";

    if (v.kind == 3) {
        draw_kinematics_comparison_canvas(linear_png,
                                          period,
                                          v,
                                          hd,
                                          hm,
                                          norm,
                                          false);
    } else if (v.kind == 2) {
        draw_phi_comparison_canvas(linear_png,
                                   period,
                                   v,
                                   hd,
                                   hm,
                                   norm,
                                   false);
    } else {
        draw_sector_comparison_canvas(linear_png,
                                      period,
                                      v,
                                      hd,
                                      hm,
                                      norm,
                                      false);
    }

    std::vector<TGraphErrors*> gr = make_ratio_graphs(hd, hm);

    if (v.kind == 3) {
        draw_kinematics_ratio_canvas(ratio_png,
                                     period,
                                     v,
                                     gr,
                                     norm,
                                     fits_for_p1_theta);
    } else if (v.kind == 2) {
        draw_phi_ratio_canvas(ratio_png,
                              period,
                              v,
                              gr,
                              norm,
                              fits_for_p1_theta);
    } else {
        draw_sector_ratio_canvas(ratio_png,
                                 period,
                                 v,
                                 gr,
                                 norm,
                                 fits_for_p1_theta);
    }

    for (TGraphErrors* g : gr) {
        delete g;
    }
}

// -----------------------------------------------------------------------------
// Period data/MC normalization histograms
// -----------------------------------------------------------------------------

struct PeriodHistResult {
    std::string period;
    HistSet data_hists;
    HistSet mc_hists;
    CubicFit poly;
    double integrated_ratio = 1.0;
    double integrated_ratio_err = 0.0;

    PeriodHistResult() :
        data_hists(make_hist_set("data_empty")),
        mc_hists(make_hist_set("mc_empty")) {
    }
};

static double data_charge_for_tree(TTree* tree,
                                   const std::unordered_map<int, double>& charge_map) {
    if (!tree) {
        return 0.0;
    }

    Branches b;
    b.bind(tree, false);

    if (!b.has_runnum) {
        fatal("[eppi0_norm] FATAL: data tree missing runnum.");
    }

    std::set<int> runs;
    const Long64_t N = tree->GetEntries();

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        runs.insert(b.runnum);
    }

    double q = 0.0;

    for (int run : runs) {
        auto it = charge_map.find(run);

        if (it == charge_map.end()) {
            std::ostringstream ss;
            ss << "[eppi0_norm] FATAL: run " << run << " missing from charge CSV.";
            fatal(ss.str());
        }

        q += it->second;
    }

    return q;
}

static void fill_eppi0_data_hists_analysis(const ChannelConfig& cfg,
                                           const PeriodTags& tags,
                                           TTree* tree,
                                           const TopoCutMap& data_cuts,
                                           double current_factor,
                                           HistSet& hists) {
    if (!tree) {
        return;
    }

    Branches b;
    b.bind(tree, false);

    const Long64_t N = tree->GetEntries();
    const double w = 1.0 / current_factor;

    long long n_pass = 0;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (!passes_event_selection(cfg, tags, data_cuts, b)) {
            continue;
        }

        ++n_pass;
        fill_hist_set(hists, b, w);
    }

    std::cout << "[eppi0_norm] analysis DATA " << tags.display
              << " entries=" << (long long)N
              << " pass=" << n_pass
              << " current_factor=" << current_factor
              << std::endl;
}

static void fill_eppi0_mc_hists_analysis(const ChannelConfig& cfg,
                                         const PeriodTags& tags,
                                         TTree* tree,
                                         const TopoCutMap& mc_cuts,
                                         double event_norm,
                                         double current_factor,
                                         HistSet& hists) {
    if (!tree) {
        return;
    }

    Branches b;
    b.bind(tree, true);

    const Long64_t N = tree->GetEntries();

    long long n_pass = 0;
    double sum_weight = 0.0;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (!passes_event_selection(cfg, tags, mc_cuts, b)) {
            continue;
        }

        ++n_pass;

        const double w = event_norm / current_factor;

        fill_hist_set(hists, b, w);
        sum_weight += w;
    }

    std::cout << "[eppi0_norm] analysis MC " << tags.display
              << " entries=" << (long long)N
              << " pass=" << n_pass
              << " sum_weight=" << sum_weight
              << " event_norm=" << event_norm
              << " current_factor=" << current_factor
              << std::endl;
}

static PeriodNormalization run_period_normalization(const std::string& period,
                                                    TTree* data_tree,
                                                    TTree* rec_tree,
                                                    const CSV& csv,
                                                    const std::unordered_map<int, double>& charge_map,
                                                    const TopoCutMap& data_cuts,
                                                    const TopoCutMap& mc_cuts,
                                                    const std::string& output_dir,
                                                    const std::string& normalization_json_path) {
    const ChannelConfig epi = eppi0_config();
    const PeriodTags parsed = parse_period_from_key(period);

    PeriodTags tags;
    tags.display = period;
    tags.period_label = parsed.period_label;
    tags.period_code = parsed.period_code;

    const double data_eff = read_current_factor(csv, epi, "exp", period);
    const double mc_eff = read_current_factor(csv, epi, "mc", period);

    const double q_raw = data_charge_for_tree(data_tree, charge_map);
    const double q_mc = q_raw * CHARGE_TO_MC_FACTOR;
    const double lint = RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC * q_mc;

    const AaoPeriodInput aao_input = read_aao_period_input(normalization_json_path, period);
    const double n_gen = aao_input.generated_events;
    const double sigma_integrated_pb = aao_input.sigma_mid_microbarn * 1.0e6;
    const double expected_generated_yield = lint * sigma_integrated_pb;
    const double event_norm = expected_generated_yield / n_gen;

    PlotNormInfo norm_info;
    norm_info.total_charge_mc = q_mc;
    norm_info.integrated_luminosity_pb_inv = lint;
    norm_info.n_gen = n_gen;
    norm_info.event_norm = event_norm;

    const std::string period_root = output_dir + "/" + period_dir(period);
    mkdir_p(period_root);

    HistSet ana_data = make_hist_set("ana_data_" + period_dir(period));
    HistSet ana_mc = make_hist_set("ana_mc_" + period_dir(period));

    fill_eppi0_data_hists_analysis(epi,
                                   tags,
                                   data_tree,
                                   data_cuts,
                                   data_eff,
                                   ana_data);

    fill_eppi0_mc_hists_analysis(epi,
                                 tags,
                                 rec_tree,
                                 mc_cuts,
                                 event_norm,
                                 mc_eff,
                                 ana_mc);

    std::map<std::string, SummaryRatioCurve> summary_ratio_curves;

    for (const VarConfig& v : variable_configs()) {
        if (v.key == "minus_t1" || v.key == "Q2" || v.key == "x") {
            summary_ratio_curves[v.key] = make_summary_ratio_curve(period,
                                                                    v,
                                                                    ana_data.hists[v.key],
                                                                    ana_mc.hists[v.key]);
        }
    }

    for (const VarConfig& v : variable_configs()) {
        plot_variable(period_root,
                      period,
                      v,
                      ana_data.hists[v.key],
                      ana_mc.hists[v.key],
                      norm_info,
                      nullptr);
    }

    std::map<std::string, TH1D*> h_data_fit;
    std::map<std::string, TH1D*> h_mc_fit;

    build_p1_theta_fit_histograms_by_region(ana_data,
                                            ana_mc,
                                            period,
                                            h_data_fit,
                                            h_mc_fit);

    std::map<std::string, RegionNormalization> region_norms;

    for (const std::string& region : normalization_regions()) {
        CubicFit fit = fit_p1_theta_ratio_region(period,
                                                 region,
                                                 period_root,
                                                 h_data_fit[region],
                                                 h_mc_fit[region]);

        const double data_int_region = h_data_fit[region]->Integral();
        const double mc_int_region = h_mc_fit[region]->Integral();

        double data_err2_region = 0.0;
        double mc_err2_region = 0.0;

        for (int b = 1; b <= h_data_fit[region]->GetNbinsX(); ++b) {
            data_err2_region += h_data_fit[region]->GetBinError(b) * h_data_fit[region]->GetBinError(b);
            mc_err2_region += h_mc_fit[region]->GetBinError(b) * h_mc_fit[region]->GetBinError(b);
        }

        double ratio_region = 1.0;
        double ratio_err_region = 0.0;

        if (data_int_region > 0.0 && mc_int_region > 0.0) {
            ratio_region = data_int_region / mc_int_region;
            ratio_err_region = std::fabs(ratio_region) * std::sqrt(
                data_err2_region / (data_int_region * data_int_region) +
                mc_err2_region / (mc_int_region * mc_int_region)
            );
        }

        RegionNormalization rn;
        rn.region = region;
        rn.fit = fit;
        rn.integrated_ratio = ratio_region;
        rn.integrated_ratio_err = ratio_err_region;

        region_norms[region] = rn;
    }

    for (const VarConfig& v : variable_configs()) {
        if (v.key == "p1_theta") {
            std::map<std::string, CubicFit> fits_for_plot;

            for (const auto& kv : region_norms) {
                fits_for_plot[kv.first] = kv.second.fit;
            }

            plot_variable(period_root,
                          period,
                          v,
                          ana_data.hists[v.key],
                          ana_mc.hists[v.key],
                          norm_info,
                          &fits_for_plot);
        }
    }

    double data_int = 0.0;
    double mc_int = 0.0;
    double data_err2 = 0.0;
    double mc_err2 = 0.0;

    for (const std::string& region : normalization_regions()) {
        data_int += h_data_fit[region]->Integral();
        mc_int += h_mc_fit[region]->Integral();

        for (int b = 1; b <= h_data_fit[region]->GetNbinsX(); ++b) {
            data_err2 += h_data_fit[region]->GetBinError(b) * h_data_fit[region]->GetBinError(b);
            mc_err2 += h_mc_fit[region]->GetBinError(b) * h_mc_fit[region]->GetBinError(b);
        }
    }

    double ratio = 1.0;
    double ratio_err = 0.0;

    if (data_int > 0.0 && mc_int > 0.0) {
        ratio = data_int / mc_int;
        ratio_err = std::fabs(ratio) * std::sqrt(
            data_err2 / (data_int * data_int) +
            mc_err2 / (mc_int * mc_int)
        );
    }

    std::ofstream dbg((period_root + "/normalization_debug_summary.txt").c_str());

    if (dbg.is_open()) {
        dbg << "period " << period << "\n";
        dbg << "data_current_factor " << std::setprecision(12) << data_eff << "\n";
        dbg << "mc_current_factor " << std::setprecision(12) << mc_eff << "\n";
        dbg << "data_charge_nC " << std::setprecision(12) << q_raw << "\n";
        dbg << "charge_to_mc_factor " << std::setprecision(12) << CHARGE_TO_MC_FACTOR << "\n";
        dbg << "integrated_luminosity_pb_inv " << std::setprecision(12) << lint << "\n";
        dbg << "normalization_json_path " << normalization_json_path << "\n";
        dbg << "sigma_min_microbarn " << std::setprecision(12) << aao_input.sigma_min_microbarn << "\n";
        dbg << "sigma_max_microbarn " << std::setprecision(12) << aao_input.sigma_max_microbarn << "\n";
        dbg << "sigma_mid_microbarn " << std::setprecision(12) << aao_input.sigma_mid_microbarn << "\n";
        dbg << "sigma_half_width_microbarn " << std::setprecision(12) << aao_input.sigma_half_width_microbarn << "\n";
        dbg << "sigma_relative_uncertainty_from_range " << std::setprecision(12) << aao_input.sigma_relative_uncertainty_from_range << "\n";
        dbg << "n_gen " << std::setprecision(12) << n_gen << "\n";
        dbg << "expected_generated_yield " << std::setprecision(12) << expected_generated_yield << "\n";
        dbg << "event_norm " << std::setprecision(12) << event_norm << "\n";
        dbg << "analysis_data_integral_p1_theta " << std::setprecision(12) << data_int << "\n";
        dbg << "analysis_mc_integral_p1_theta " << std::setprecision(12) << mc_int << "\n";
        dbg << "analysis_ratio " << std::setprecision(12) << ratio << "\n";
        dbg << "analysis_ratio_stat " << std::setprecision(12) << ratio_err << "\n";
    }

    PeriodNormalization out;
    out.period = period;
    out.regions = region_norms;
    out.summary_ratio_curves = summary_ratio_curves;

    for (auto& kv : h_data_fit) {
        delete kv.second;
    }

    for (auto& kv : h_mc_fit) {
        delete kv.second;
    }

    delete_hist_set(ana_data);
    delete_hist_set(ana_mc);

    std::cout << "[eppi0_norm] " << period
              << " regional_cubic_fits=" << region_norms.size()
              << " integrated_ratio=" << ratio
              << " +/- " << ratio_err
              << " data_eff=" << data_eff
              << " mc_eff=" << mc_eff
              << " sigma_mid_microbarn=" << aao_input.sigma_mid_microbarn
              << " n_gen=" << n_gen
              << " event_norm=" << event_norm
              << std::endl;

    return out;
}

static std::string key_for_period_map(const std::map<std::string, TTree*>& trees,
                                      const std::string& period) {
    for (const auto& kv : trees) {
        PeriodTags tags = parse_period_from_key(kv.first);

        if (tags.display == period) {
            return kv.first;
        }
    }

    return "";
}

static TTree* tree_for_period(const std::map<std::string, TTree*>& trees,
                              const std::string& period,
                              const std::string& label) {
    const std::string key = key_for_period_map(trees, period);

    if (key.empty()) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: missing " << label << " tree for period " << period;
        fatal(ss.str());
    }

    return trees.at(key);
}

// -----------------------------------------------------------------------------
// Normalized raw yield filling
// -----------------------------------------------------------------------------

struct HelCounts {
    double unpol = 0.0;
    double pos = 0.0;
    double neg = 0.0;
};

using RowCounts = std::unordered_map<int, HelCounts>;

static void accumulate_normalized_yields_for_tree(const ChannelConfig& cfg,
                                                  const PeriodTags& tags,
                                                  TTree* tree,
                                                  const std::vector<RowBin>& rows,
                                                  const TopoCutMap& data_cuts,
                                                  double current_factor,
                                                  const PeriodNormalization& norm,
                                                  std::unordered_map<std::string, RowCounts>& out) {
    if (!tree) {
        return;
    }

    Branches b;
    b.bind(tree, false);

    if (!b.has_helicity) {
        fatal("[eppi0_norm] FATAL: data tree missing helicity.");
    }

    const Long64_t N = tree->GetEntries();

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (!passes_event_selection(cfg, tags, data_cuts, b)) {
            continue;
        }

        const std::string topo = topo_dir(b.detector1, b.detector2);

        if (topo.empty()) {
            continue;
        }

        const double theta = b.p1_theta_deg();
        const std::string region = normalization_region_for_event(b);

        if (region.empty()) {
            continue;
        }

        auto it_region = norm.regions.find(region);

        if (it_region == norm.regions.end()) {
            std::ostringstream ss;
            ss << "[eppi0_norm] FATAL: missing normalization region " << region
               << " for period " << norm.period;
            fatal(ss.str());
        }

        const double ratio = it_region->second.fit.eval(theta);

        if (!(std::isfinite(ratio) && ratio > 0.0)) {
            continue;
        }

        const double weight = 1.0 / (current_factor * ratio);
        const double phi = b.phi_deg();
        const double tabs = b.t_abs();

        for (int r = 0; r < (int)rows.size(); ++r) {
            if (!row_accepts(rows[r], b.x, b.Q2, tabs, phi)) {
                continue;
            }

            HelCounts& hc = out[topo][r];

            if (b.helicity > 0) {
                hc.pos += weight;
            } else if (b.helicity < 0) {
                hc.neg += weight;
            } else {
                hc.unpol += weight;
            }
        }
    }
}

static void write_normalized_counts_to_csv(CSV& csv,
                                           const ChannelConfig& cfg,
                                           const std::string& period,
                                           const std::unordered_map<std::string, RowCounts>& topo_counts) {
    for (const auto& kt : topo_counts) {
        const std::string topo_label = topo_label_for_csv(kt.first);

        if (topo_label.empty()) {
            continue;
        }

        const int c_unpol = col_strict(csv, normalized_raw_yield_col(cfg, topo_label, period, "unpol"));
        const int c_pos = col_strict(csv, normalized_raw_yield_col(cfg, topo_label, period, "pos"));
        const int c_neg = col_strict(csv, normalized_raw_yield_col(cfg, topo_label, period, "neg"));

        for (const auto& kr : kt.second) {
            const int row = kr.first;

            if (row < 0 || row >= (int)csv.rows.size()) {
                fatal("[eppi0_norm] FATAL: row index out of range.");
            }

            const HelCounts& h = kr.second;
            const double unpol = h.pos + h.neg + h.unpol;

            std::ostringstream su;
            std::ostringstream sp;
            std::ostringstream sn;

            su << std::fixed << std::setprecision(8) << unpol;
            sp << std::fixed << std::setprecision(8) << h.pos;
            sn << std::fixed << std::setprecision(8) << h.neg;

            csv.rows[row][c_unpol] = su.str();
            csv.rows[row][c_pos] = sp.str();
            csv.rows[row][c_neg] = sn.str();
        }
    }
}

static const PeriodNormalization& find_norm(const std::vector<PeriodNormalization>& norms,
                                            const std::string& period) {
    for (const auto& n : norms) {
        if (n.period == period) {
            return n;
        }
    }

    std::ostringstream ss;
    ss << "[eppi0_norm] FATAL: missing normalization result for " << period;
    fatal(ss.str());

    return norms.front();
}

static void fill_normalized_yields(CSV& csv,
                                   const std::vector<RowBin>& rows,
                                   const std::map<std::string, TTree*>& dvcsDataTrees,
                                   const std::map<std::string, TTree*>& eppi0DataTrees,
                                   const std::vector<PeriodNormalization>& norms,
                                   const TopoCutMap& data_cuts,
                                   int max_workers) {
    struct WorkItem {
        ChannelConfig cfg;
        PeriodTags tags;
        TTree* tree = nullptr;
        double current_factor = 1.0;
        PeriodNormalization norm;
    };

    std::vector<WorkItem> items;

    const ChannelConfig dvcs = dvcs_config();
    const ChannelConfig epi = eppi0_config();

    auto add_items = [&](const ChannelConfig& cfg, const std::map<std::string, TTree*>& trees) {
        for (const auto& kv : trees) {
            PeriodTags tags = parse_period_from_key(kv.first);

            if (tags.display == "Fa18 Inb Supp") {
                continue;
            }

            WorkItem w;
            w.cfg = cfg;
            w.tags = tags;
            w.tree = kv.second;
            w.current_factor = read_current_factor(csv, cfg, "exp", tags.display);
            w.norm = find_norm(norms, tags.display);

            items.push_back(w);
        }
    };

    add_items(dvcs, dvcsDataTrees);
    add_items(epi, eppi0DataTrees);

    std::vector<std::unordered_map<std::string, RowCounts>> results(items.size());

    int nth = std::max(1, std::min(5, max_workers));
    nth = std::min(nth, std::max(1, (int)items.size()));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int i = 0; i < (int)items.size(); ++i) {
        accumulate_normalized_yields_for_tree(items[i].cfg,
                                              items[i].tags,
                                              items[i].tree,
                                              rows,
                                              data_cuts,
                                              items[i].current_factor,
                                              items[i].norm,
                                              results[i]);
    }

    for (int i = 0; i < (int)items.size(); ++i) {
        write_normalized_counts_to_csv(csv,
                                       items[i].cfg,
                                       items[i].tags.display,
                                       results[i]);
    }
}

// -----------------------------------------------------------------------------
// Parent summary plots
// -----------------------------------------------------------------------------

static std::string region_short_label(const std::string& region) {
    if (region == "CD") {
        return "CD";
    }

    return region;
}

static int color_for_period_index(int i) {
    static const int colors[5] = {
        kBlack,
        kRed + 1,
        kBlue + 1,
        kGreen + 2,
        kMagenta + 1
    };

    if (i < 0) {
        return kBlack;
    }

    return colors[i % 5];
}

static const PeriodNormalization* find_norm_ptr(const std::vector<PeriodNormalization>& norms,
                                                const std::string& period) {
    for (const PeriodNormalization& n : norms) {
        if (n.period == period) {
            return &n;
        }
    }

    return nullptr;
}

static const RegionNormalization* find_region_norm_ptr(const PeriodNormalization& norm,
                                                       const std::string& region) {
    auto it = norm.regions.find(region);

    if (it == norm.regions.end()) {
        return nullptr;
    }

    return &it->second;
}

static void draw_fit_canvas_frame(TPad* pad,
                                  const std::string& panel_label,
                                  double x_min,
                                  double x_max,
                                  const std::string& title_suffix) {
    pad->cd();
    pad->SetGrid(1, 1);
    pad->SetLeftMargin(0.16);
    pad->SetRightMargin(0.05);
    pad->SetBottomMargin(0.16);
    pad->SetTopMargin(0.09);

    TH1D* frame = (TH1D*)pad->DrawFrame(x_min, RATIO_Y_MIN, x_max, RATIO_Y_MAX);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("p_{1} #theta (deg)");
    frame->GetYaxis()->SetTitle("data / MC");
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);

    TLine* unity_line = new TLine(x_min, 1.0, x_max, 1.0);
    unity_line->SetLineStyle(2);
    unity_line->SetLineColor(kRed + 1);
    unity_line->SetLineWidth(2);
    unity_line->SetBit(TObject::kCanDelete);
    unity_line->Draw("SAME");

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextSize(0.052);
    latex.SetTextFont(42);
    latex.DrawLatex(0.20, 0.84, (panel_label + title_suffix).c_str());
}

static void save_all_period_p1_theta_fit_summary(const std::string& output_dir,
                                                 const std::vector<PeriodNormalization>& norms) {
    mkdir_p(output_dir);

    TCanvas c("c_eppi0_p1_theta_fit_summary_all_periods", "", 1800, 900);
    c.Divide(4, 2, 0.001, 0.001);

    std::vector<TF1*> owned_funcs;

    const std::vector<std::string> regions = normalization_regions();

    for (int ipanel = 0; ipanel < 7; ++ipanel) {
        const int pad_index = (ipanel < 3) ? (ipanel + 1) : (ipanel + 2);
        TPad* pad = (TPad*)c.cd(pad_index);

        const std::string& region = regions[ipanel];
        const bool is_cd = (region == "CD");

        const double x_min = is_cd ? 40.0 : 0.0;
        const double x_max = is_cd ? 70.0 : 40.0;

        draw_fit_canvas_frame(pad,
                              region_short_label(region),
                              x_min,
                              x_max,
                              "");

        for (int ip = 0; ip < (int)CSV_PERIOD_ORDER.size(); ++ip) {
            const std::string& period = CSV_PERIOD_ORDER[ip];

            const PeriodNormalization* pn = find_norm_ptr(norms, period);

            if (!pn) {
                continue;
            }

            const RegionNormalization* rn = find_region_norm_ptr(*pn, region);

            if (!rn || !rn->fit.valid) {
                continue;
            }

            const CubicFit& fit = rn->fit;

            if (!(fit.x_max > fit.x_min)) {
                continue;
            }

            TF1* f = new TF1(("f_all_period_" + period_dir(period) + "_" + std::to_string(ipanel)).c_str(),
                             "pol3",
                             fit.x_min,
                             fit.x_max);

            for (int k = 0; k < 4; ++k) {
                f->SetParameter(k, fit.a[k]);
            }

            f->SetLineColor(color_for_period_index(ip));
            f->SetLineWidth(3);
            f->Draw("L SAME");

            owned_funcs.push_back(f);
        }
    }

    TPad* info = (TPad*)c.cd(4);
    info->SetLeftMargin(0.05);
    info->SetRightMargin(0.05);
    info->SetTopMargin(0.05);
    info->SetBottomMargin(0.05);
    info->Clear();

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(42);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.08, 0.90, "ep #rightarrow ep#pi_{0}");
    latex.DrawLatex(0.08, 0.82, "p_{1} #theta normalization fits");
    latex.DrawLatex(0.08, 0.74, "Pass-2 cubic fits");
    latex.DrawLatex(0.08, 0.66, "Line ranges: populated data/MC bins only");

    TLegend* legend = new TLegend(0.08, 0.20, 0.88, 0.58);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->SetBorderSize(1);
    legend->SetTextSize(0.045);

    std::vector<TF1*> legend_funcs;

    for (int ip = 0; ip < (int)CSV_PERIOD_ORDER.size(); ++ip) {
        TF1* dummy = new TF1(("f_legend_period_" + std::to_string(ip)).c_str(), "pol0", 0.0, 1.0);
        dummy->SetParameter(0, 1.0);
        dummy->SetLineColor(color_for_period_index(ip));
        dummy->SetLineWidth(3);

        legend_funcs.push_back(dummy);
        legend->AddEntry(dummy, CSV_PERIOD_ORDER[ip].c_str(), "l");
    }

    legend->Draw();

    c.SaveAs((output_dir + "/eppi0_p1_theta_cubic_fits_all_periods.png").c_str());

    for (TF1* f : owned_funcs) {
        delete f;
    }

    for (TF1* f : legend_funcs) {
        delete f;
    }

    delete legend;
}

static bool pass1_quartic_coefficients(const std::string& period,
                                       const std::string& region,
                                       QuarticFit& q) {
    const int region_index = (region == "Sector 1") ? 0 :
                             (region == "Sector 2") ? 1 :
                             (region == "Sector 3") ? 2 :
                             (region == "Sector 4") ? 3 :
                             (region == "Sector 5") ? 4 :
                             (region == "Sector 6") ? 5 :
                             (region == "CD")       ? 6 : -1;

    if (region_index < 0) {
        return false;
    }

    static const double inb[7][5] = {
        { 5.205723e-6, -6.553104e-4,  2.857282e-2, -5.275748e-1,  4.390724 },
        { 1.736551e-5, -1.905288e-3,  7.459240e-2, -1.243871e0,  8.406038 },
        {-3.858605e-6,  4.012865e-4, -1.489836e-2,  2.088062e-1,  8.284200e-2},
        { 8.245591e-6, -8.403340e-4,  2.952225e-2, -4.468835e-1,  3.527616 },
        {-7.027999e-6,  4.804512e-4, -8.624948e-3, -3.522498e-2,  2.414833 },
        {-9.300985e-6,  8.132209e-4, -2.646195e-2,  3.715724e-1, -9.530212e-1},
        { 2.430049e-5, -5.364220e-3,  4.440329e-1, -1.633226e1,  2.258644e2}
    };

    static const double out[7][5] = {
        { 6.193140e-5, -7.809714e-3,  3.609521e-1, -7.247023e0,  5.419993e1},
        {-3.086642e-5,  3.301421e-3, -1.326819e-1,  2.385607e0, -1.532828e1},
        { 1.265237e-5, -2.382045e-3,  1.419987e-1, -3.430274e0,  3.015195e1},
        {-3.351760e-5,  3.303882e-3, -1.176380e-1,  1.776518e0, -8.604291e0},
        { 2.073593e-5, -3.141968e-3,  1.675443e-1, -3.791117e0,  3.185059e1},
        {-1.690791e-4,  1.976225e-2, -8.510287e-1,  1.595090e1, -1.084887e2},
        { 7.497258e-6, -1.734757e-3,  1.505486e-1, -5.804050e0,  8.456443e1}
    };

    const double (*src)[5] = nullptr;

    if (period == "Fa18 Inb") {
        src = inb;
    } else if (period == "Fa18 Out") {
        src = out;
    } else {
        return false;
    }

    q.a[0] = src[region_index][4];
    q.a[1] = src[region_index][3];
    q.a[2] = src[region_index][2];
    q.a[3] = src[region_index][1];
    q.a[4] = src[region_index][0];

    return true;
}

static void save_pass1_pass2_p1_theta_comparison(const std::string& output_dir,
                                                 const std::vector<PeriodNormalization>& norms,
                                                 const std::string& period) {
    mkdir_p(output_dir);

    const PeriodNormalization* pn = find_norm_ptr(norms, period);

    if (!pn) {
        std::cout << "[eppi0_norm] WARNING: skipping pass-1 comparison for missing period "
                  << period << std::endl;
        return;
    }

    TCanvas c(("c_eppi0_p1_theta_pass1_pass2_" + period_dir(period)).c_str(), "", 1800, 900);
    c.Divide(4, 2, 0.001, 0.001);

    std::vector<TF1*> owned_funcs;
    const std::vector<std::string> regions = normalization_regions();

    for (int ipanel = 0; ipanel < 7; ++ipanel) {
        const int pad_index = (ipanel < 3) ? (ipanel + 1) : (ipanel + 2);
        TPad* pad = (TPad*)c.cd(pad_index);

        const std::string& region = regions[ipanel];
        const bool is_cd = (region == "CD");

        const double x_frame_min = is_cd ? 40.0 : 0.0;
        const double x_frame_max = is_cd ? 70.0 : 40.0;

        draw_fit_canvas_frame(pad,
                              region_short_label(region),
                              x_frame_min,
                              x_frame_max,
                              "");

        const RegionNormalization* rn = find_region_norm_ptr(*pn, region);

        if (!rn || !rn->fit.valid) {
            continue;
        }

        const CubicFit& cfit = rn->fit;

        if (!(cfit.x_max > cfit.x_min)) {
            continue;
        }

        TF1* fpass2 = new TF1(("f_pass2_" + period_dir(period) + "_" + std::to_string(ipanel)).c_str(),
                              "pol3",
                              cfit.x_min,
                              cfit.x_max);

        for (int k = 0; k < 4; ++k) {
            fpass2->SetParameter(k, cfit.a[k]);
        }

        fpass2->SetLineColor(kBlue + 1);
        fpass2->SetLineWidth(3);
        fpass2->Draw("L SAME");

        owned_funcs.push_back(fpass2);

        QuarticFit qfit;

        if (pass1_quartic_coefficients(period, region, qfit)) {
            TF1* fpass1 = new TF1(("f_pass1_" + period_dir(period) + "_" + std::to_string(ipanel)).c_str(),
                                  "pol4",
                                  cfit.x_min,
                                  cfit.x_max);

            for (int k = 0; k < 5; ++k) {
                fpass1->SetParameter(k, qfit.a[k]);
            }

            fpass1->SetLineColor(kRed + 1);
            fpass1->SetLineWidth(3);
            fpass1->SetLineStyle(2);
            fpass1->Draw("L SAME");

            owned_funcs.push_back(fpass1);
        }
    }

    TPad* info = (TPad*)c.cd(4);
    info->SetLeftMargin(0.05);
    info->SetRightMargin(0.05);
    info->SetTopMargin(0.05);
    info->SetBottomMargin(0.05);
    info->Clear();

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(42);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.08, 0.90, period.c_str());
    latex.DrawLatex(0.08, 0.82, "ep #rightarrow ep#pi_{0}");
    latex.DrawLatex(0.08, 0.74, "p_{1} #theta normalization");
    latex.DrawLatex(0.08, 0.66, "Pass-2 cubic vs pass-1 quartic");
    latex.DrawLatex(0.08, 0.58, "Both lines over pass-2 populated range");

    TLegend legend(0.08, 0.30, 0.88, 0.50);
    legend.SetFillColor(kWhite);
    legend.SetFillStyle(1001);
    legend.SetBorderSize(1);
    legend.SetTextSize(0.050);

    TF1 leg_pass2("leg_pass2", "pol0", 0.0, 1.0);
    leg_pass2.SetParameter(0, 1.0);
    leg_pass2.SetLineColor(kBlue + 1);
    leg_pass2.SetLineWidth(3);

    TF1 leg_pass1("leg_pass1", "pol0", 0.0, 1.0);
    leg_pass1.SetParameter(0, 1.0);
    leg_pass1.SetLineColor(kRed + 1);
    leg_pass1.SetLineWidth(3);
    leg_pass1.SetLineStyle(2);

    legend.AddEntry(&leg_pass2, "pass-2 cubic", "l");
    legend.AddEntry(&leg_pass1, "pass-1 quartic", "l");
    legend.Draw();

    c.SaveAs((output_dir + "/eppi0_p1_theta_pass1_pass2_" + period_dir(period) + ".png").c_str());

    for (TF1* f : owned_funcs) {
        delete f;
    }
}

static void save_all_period_kinematic_ratio_summary(const std::string& output_dir,
                                                    const std::vector<PeriodNormalization>& norms,
                                                    const std::string& key,
                                                    const std::string& title,
                                                    const std::string& x_title) {
    mkdir_p(output_dir);

    double x_min = 0.0;
    double x_max = 1.0;

    VarConfig frame_var;
    frame_var.key = key;
    frame_var.title = title;
    frame_var.x_title = x_title;
    frame_var.particle = 0;
    frame_var.kind = 3;
    frame_var.n_panels = 1;

    range_for_panel(frame_var, 0, x_min, x_max);

    double y_max = 0.0;

    for (const PeriodNormalization& norm : norms) {
        auto it = norm.summary_ratio_curves.find(key);

        if (it == norm.summary_ratio_curves.end()) {
            continue;
        }

        const SummaryRatioCurve& curve = it->second;

        for (size_t i = 0; i < curve.y.size(); ++i) {
            const double candidate = curve.y[i] + curve.ey[i];

            if (std::isfinite(candidate)) {
                y_max = std::max(y_max, candidate);
            }
        }
    }

    if (!(y_max > 0.0)) {
        y_max = RATIO_Y_MAX;
    }

    y_max = std::max(RATIO_Y_MAX, 1.20 * y_max);

    TCanvas c(("c_eppi0_summary_ratio_" + key).c_str(), "", 1200, 850);
    c.SetGrid(1, 1);
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.08);

    TH1D* frame = (TH1D*)gPad->DrawFrame(x_min, RATIO_Y_MIN, x_max, y_max);
    frame->SetTitle(("ep #rightarrow ep#pi_{0}  all periods  " + title + " ratio").c_str());
    frame->GetXaxis()->SetTitle(x_title.c_str());
    frame->GetYaxis()->SetTitle("data / MC");
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetTitleSize(0.050);
    frame->GetYaxis()->SetTitleSize(0.050);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.25);

    TLine unity_line(x_min, 1.0, x_max, 1.0);
    unity_line.SetLineStyle(2);
    unity_line.SetLineColor(kGray + 2);
    unity_line.SetLineWidth(2);
    unity_line.Draw("SAME");

    std::vector<TGraphErrors*> graphs;

    TLegend legend(0.58, 0.66, 0.90, 0.88);
    legend.SetFillColor(kWhite);
    legend.SetFillStyle(1001);
    legend.SetBorderSize(1);
    legend.SetTextSize(0.035);

    for (int ip = 0; ip < (int)CSV_PERIOD_ORDER.size(); ++ip) {
        const std::string& period = CSV_PERIOD_ORDER[ip];
        const PeriodNormalization* pn = find_norm_ptr(norms, period);

        if (!pn) {
            continue;
        }

        auto it = pn->summary_ratio_curves.find(key);

        if (it == pn->summary_ratio_curves.end()) {
            continue;
        }

        const SummaryRatioCurve& curve = it->second;

        if (curve.x.empty()) {
            continue;
        }

        TGraphErrors* g = new TGraphErrors();

        for (int i = 0; i < (int)curve.x.size(); ++i) {
            g->SetPoint(i, curve.x[i], curve.y[i]);
            g->SetPointError(i, 0.0, curve.ey[i]);
        }

        const int color = color_for_period_index(ip);

        g->SetMarkerStyle(20 + ip);
        g->SetMarkerSize(0.85);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->SetLineWidth(1);
        g->Draw("PE SAME");

        legend.AddEntry(g, period.c_str(), "pe");
        graphs.push_back(g);
    }

    legend.Draw();

    c.SaveAs((output_dir + "/eppi0_summary_ratio_" + key + ".png").c_str());

    for (TGraphErrors* g : graphs) {
        delete g;
    }
}

static void save_all_period_kinematic_ratio_summaries(const std::string& output_dir,
                                                      const std::vector<PeriodNormalization>& norms) {
    if (norms.empty()) {
        return;
    }

    save_all_period_kinematic_ratio_summary(output_dir,
                                            norms,
                                            "minus_t1",
                                            "-t",
                                            "-t (GeV^{2})");

    save_all_period_kinematic_ratio_summary(output_dir,
                                            norms,
                                            "Q2",
                                            "Q^{2}",
                                            "Q^{2} (GeV^{2})");

    save_all_period_kinematic_ratio_summary(output_dir,
                                            norms,
                                            "x",
                                            "x_{B}",
                                            "x_{B}");
}

static void save_parent_fit_summary_plots(const std::string& output_dir,
                                          const std::vector<PeriodNormalization>& norms) {
    if (norms.empty()) {
        return;
    }

    mkdir_p(output_dir);

    save_all_period_p1_theta_fit_summary(output_dir, norms);
    save_all_period_kinematic_ratio_summaries(output_dir, norms);
    save_pass1_pass2_p1_theta_comparison(output_dir, norms, "Fa18 Inb");
    save_pass1_pass2_p1_theta_comparison(output_dir, norms, "Fa18 Out");
}

// -----------------------------------------------------------------------------
// CSV output for normalization fits
// -----------------------------------------------------------------------------

static void write_norms_to_csv(CSV& csv, const std::vector<PeriodNormalization>& norms) {
    for (const std::string& period : CSV_PERIOD_ORDER) {
        const PeriodNormalization& n = find_norm(norms, period);

        for (const std::string& region : normalization_regions()) {
            auto it_region = n.regions.find(region);

            if (it_region == n.regions.end()) {
                std::ostringstream ss;
                ss << "[eppi0_norm] FATAL: missing " << region << " normalization for " << period;
                fatal(ss.str());
            }

            const int c_fit = col_strict(csv, norm_fit_col(period, region));
            const int c_factor = col_strict(csv, norm_factor_col(period, region));

            const std::string fit_cell = tuple4(it_region->second.fit);
            const std::string factor_cell = tuple2(it_region->second.integrated_ratio,
                                                   it_region->second.integrated_ratio_err);

            for (auto& row : csv.rows) {
                row[c_fit] = fit_cell;
                row[c_factor] = factor_cell;
            }
        }
    }
}

static std::vector<PeriodNormalization> unity_norms() {
    std::vector<PeriodNormalization> norms;

    for (const std::string& period : CSV_PERIOD_ORDER) {
        PeriodNormalization n;
        n.period = period;

        for (const std::string& region : normalization_regions()) {
            RegionNormalization rn;
            rn.region = region;
            rn.fit.valid = true;
            rn.fit.a[0] = 1.0;
            rn.fit.x_min = (region == "CD") ? 40.0 : 0.0;
            rn.fit.x_max = (region == "CD") ? 70.0 : 40.0;
            rn.integrated_ratio = 1.0;
            rn.integrated_ratio_err = 0.0;

            n.regions[region] = rn;
        }

        norms.push_back(n);
    }

    return norms;
}

} // namespace

bool update_eppi0_normalization_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const Eppi0NormalizationOptions& options) {
    try {
        ROOT::EnableThreadSafety();
        TH1::AddDirectory(kFALSE);
        gStyle->SetOptStat(0);

        CSV csv;
        load_csv_strict(csv_path, csv);

        const std::vector<RowBin> rows = load_rows(csv);

        const TopoCutMap data_cuts = load_sigma_cuts(options.combined_cuts_json, "data");
        const TopoCutMap mc_cuts = load_sigma_cuts(options.combined_cuts_json, "mc");

        std::vector<PeriodNormalization> norms;

        if (options.override_to_unity) {
            std::cout << "[eppi0_norm] Override enabled: using unity eppi0 normalization polynomial."
                      << std::endl;

            norms = unity_norms();
        } else {
            mkdir_p(options.output_dir);

            const std::unordered_map<int, double> charge_map =
                read_charge_csv(options.charge_csv_path);

            for (const std::string& period : WORK_PERIOD_ORDER) {
                TTree* data_tree = tree_for_period(eppi0DataTrees, period, "eppi0 data");
                TTree* rec_tree = tree_for_period(eppi0RecMcTrees, period, "eppi0 reconstructed MC");

                PeriodNormalization n = run_period_normalization(period,
                                                                 data_tree,
                                                                 rec_tree,
                                                                 csv,
                                                                 charge_map,
                                                                 data_cuts,
                                                                 mc_cuts,
                                                                 options.output_dir,
                                                                 options.normalization_json_path);

                norms.push_back(n);
            }
        }

        if (!options.override_to_unity) {
            save_parent_fit_summary_plots(options.output_dir, norms);
        }

        write_norms_to_csv(csv, norms);

        fill_normalized_yields(csv,
                               rows,
                               dvcsDataTrees,
                               eppi0DataTrees,
                               norms,
                               data_cuts,
                               options.max_workers);

        write_csv_atomic(csv_path, csv);

        std::cout << "[eppi0_norm] Updated eppi0 normalization polynomial columns and normalized raw yields in: "
                  << csv_path << std::endl;

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}