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

static constexpr double RATIO_Y_MIN = 0.2;
static constexpr double RATIO_Y_MAX = 1.5;

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

struct Poly4 {
    double a[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
    double ea[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    bool valid = false;

    double eval(double x) const {
        return a[0] + x * (a[1] + x * (a[2] + x * (a[3] + x * a[4])));
    }
};

struct PeriodNormalization {
    std::string period;
    Poly4 poly;
    double integrated_ratio = 1.0;
    double integrated_ratio_err = 0.0;
};

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
    if (p < 0.0) p += 360.0;
    if (p >= 360.0) p = std::nextafter(360.0, 0.0);
    return p;
}

static bool in_range(double v, double a, double b) {
    return v >= a && v < b;
}

static bool row_accepts_phi(double phi, double pmin, double pmax) {
    if (pmax > pmin) return in_range(phi, pmin, pmax);
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
                if (ch == '"') oss << "\"\"";
                else oss << ch;
            }
            oss << '"';
        } else {
            oss << s;
        }
        if (i + 1 < fields.size()) oss << ',';
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
        if (line.empty()) continue;
        std::vector<std::string> row = split_csv_line(line);
        if (row.size() < csv.header.size()) row.resize(csv.header.size(), "");
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
    for (const auto& row : csv.rows) fout << join_csv_row(row) << "\n";
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

static std::string tuple5(const Poly4& p) {
    std::ostringstream ss;
    ss << std::setprecision(12)
       << "(" << p.a[0] << "," << p.a[1] << "," << p.a[2] << "," << p.a[3] << "," << p.a[4] << ")";
    return ss.str();
}

static bool parse_tuple_numbers(const std::string& s, std::vector<double>& vals) {
    vals.clear();
    std::string t = s;
    for (char& c : t) {
        if (c == '(' || c == ')' || c == '"') c = ' ';
    }

    std::stringstream ss(t);
    std::string part;
    while (std::getline(ss, part, ',')) {
        char* e = nullptr;
        const double v = std::strtod(part.c_str(), &e);
        if (e == part.c_str()) return false;
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
    if (csv.rows.empty()) fatal("[eppi0_norm] FATAL: CSV has no rows.");

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

static std::string norm_poly_col(const std::string& period) {
    return "eppi0 cross-section normalization polynomial, ep->eppi0, data_over_mc, " + period;
}

static std::string norm_factor_col(const std::string& period) {
    return "eppi0 cross-section normalization factor, ep->eppi0, data_over_mc, " + period;
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
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> f = split_csv_line(line);
        if (f.size() < 2) continue;
        char* e = nullptr;
        int run = (int)std::strtol(f[0].c_str(), &e, 10);
        if (e == f[0].c_str()) continue;
        e = nullptr;
        double q = std::strtod(f[1].c_str(), &e);
        if (e == f[1].c_str()) continue;
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
        if (!block.is_object() || !block.contains(sample_key)) continue;
        const auto& sample = block[sample_key];
        if (!sample.is_object()) continue;

        CutVarMap vm;
        for (auto vit = sample.begin(); vit != sample.end(); ++vit) {
            const auto& obj = vit.value();
            if (!obj.is_object() || !obj.contains("mean") || !obj.contains("std")) continue;
            SigmaStats s;
            s.mean = obj["mean"].get<double>();
            s.std = obj["std"].get<double>();
            vm[vit.key()] = s;
        }
        if (!vm.empty()) out[key] = vm;
    }
    return out;
}

static bool within_3sigma(double v, const SigmaStats& s) {
    if (!std::isfinite(s.mean) || !std::isfinite(s.std) || s.std <= 0.0) return true;
    return std::fabs(v - s.mean) <= 3.0 * s.std;
}

static bool check_sigma(const CutVarMap& vm, const std::string& var, bool has, double val) {
    auto it = vm.find(var);
    if (it == vm.end()) return true;
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
            if (t->GetBranch(name)) t->SetBranchStatus(name, 1);
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
        if (need_weight) ena("weight");

        t->SetCacheSize(0);

        auto bI = [&](const char* name, int* addr, bool& flag) {
            if (t->GetBranch(name)) { t->SetBranchAddress(name, addr); flag = true; }
        };
        auto bD = [&](const char* name, double* addr, bool& flag) {
            if (t->GetBranch(name)) { t->SetBranchAddress(name, addr); flag = true; }
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

    double phi_deg() const { return wrap_phi_deg(phi2 * RAD2DEG); }
    double t_abs() const { return std::fabs(t1); }
    double p1_theta_deg() const { return p1_theta * RAD2DEG; }
    double p1_phi_deg() const { return wrap_phi_deg(p1_phi * RAD2DEG); }
    double p2_theta_deg() const { return p2_theta * RAD2DEG; }
    double p2_phi_deg() const { return wrap_phi_deg(p2_phi * RAD2DEG); }
};

static bool passes_global_dispatch(const Branches& b, const PeriodTags& tags) {
    if (!(b.has_t1 && b.has_open_angle_ep2 && b.has_pTmiss)) return false;
    if (b.has_runnum && is_excluded_run(b.runnum)) return false;

    const GlobalCutConfig& cfg = default_global_cuts();

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_detector1 && b.has_detector2 && b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[eppi0_norm] FATAL: missing branches required by global ycol/topology cuts.");
        }
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  tags.period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    if (!(b.has_detector1 && b.has_detector2)) {
        fatal("[eppi0_norm] FATAL: missing detector1/detector2.");
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                              b.detector1, b.detector2, cfg);
}

static bool passes_sigma_dispatch(const ChannelConfig& cfg,
                                  const PeriodTags& tags,
                                  const TopoCutMap& cuts,
                                  const Branches& b) {
    if (!(b.has_detector1 && b.has_detector2)) return false;
    const std::string topo = topo_dir(b.detector1, b.detector2);
    if (topo.empty()) return false;

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
    if (!passes_global_dispatch(b, tags)) return false;
    if (!passes_sigma_dispatch(cfg, tags, cuts, b)) return false;
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
    int kind = 0; // 0 theta, 1 momentum, 2 phi
    int n_panels = 7;
};

static std::vector<VarConfig> variable_configs() {
    return {
        {"p1_theta", "p_{1} #theta", "p_{1} #theta (deg)", 1, 0, 7},
        {"p1_p", "p_{1} momentum", "p_{1} momentum (GeV)", 1, 1, 7},
        {"p1_phi", "p_{1} #phi", "p_{1} #phi (deg)", 1, 2, 2},
        {"p2_theta", "p_{2} #theta", "p_{2} #theta (deg)", 2, 0, 7},
        {"p2_p", "p_{2} momentum", "p_{2} momentum (GeV)", 2, 1, 7},
        {"p2_phi", "p_{2} #phi", "p_{2} #phi (deg)", 2, 2, 2}
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
    if (v.kind == 0) return (panel < 6) ? 12 : 23;
    if (v.kind == 1) return (panel < 6) ? 12 : 23;
    return (panel == 0) ? 12 : 24;
}

static void range_for_panel(const VarConfig& v, int panel, double& xmin, double& xmax) {
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

struct HistSet {
    std::map<std::string, std::vector<TH1D*>> hists;
};

static HistSet make_hist_set(const std::string& prefix) {
    HistSet hs;
    for (const VarConfig& v : variable_configs()) {
        std::vector<TH1D*> vec;
        for (int p = 0; p < v.n_panels; ++p) {
            double xmin = 0.0, xmax = 1.0;
            range_for_panel(v, p, xmin, xmax);
            int nb = nbins_for_panel(v, p);
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
        for (TH1D* h : kv.second) delete h;
        kv.second.clear();
    }
    hs.hists.clear();
}

static void fill_hist_set(HistSet& hs, const Branches& b, double w) {
    for (const VarConfig& v : variable_configs()) {
        int p = panel_index(v, b);
        if (p < 0 || p >= v.n_panels) continue;
        const double value = value_for_var(v, b);
        hs.hists[v.key][p]->Fill(value, w);
    }
}

static void add_hist_set(HistSet& dest, const HistSet& src) {
    for (const auto& kv : src.hists) {
        auto it = dest.hists.find(kv.first);
        if (it == dest.hists.end()) continue;
        for (size_t i = 0; i < kv.second.size(); ++i) {
            it->second[i]->Add(kv.second[i]);
        }
    }
}

// -----------------------------------------------------------------------------
// Fitting and plotting
// -----------------------------------------------------------------------------

static Poly4 fit_p1_theta_ratio(const std::string& period,
                                const std::string& outdir,
                                TH1D* h_data,
                                TH1D* h_mc) {
    Poly4 p;
    if (!h_data || !h_mc) return p;

    TGraphErrors gr;
    int ip = 0;
    for (int b = 1; b <= h_data->GetNbinsX(); ++b) {
        const double d = h_data->GetBinContent(b);
        const double m = h_mc->GetBinContent(b);
        if (!(d > 0.0 && m > 0.0)) continue;
        const double ed = h_data->GetBinError(b);
        const double em = h_mc->GetBinError(b);
        const double r = d / m;
        double er = 0.0;
        if (ed > 0.0 || em > 0.0) {
            const double rd = (ed > 0.0) ? ed / d : 0.0;
            const double rm = (em > 0.0) ? em / m : 0.0;
            er = std::fabs(r) * std::sqrt(rd * rd + rm * rm);
        }
        if (!(er > 0.0)) er = 1.0;
        gr.SetPoint(ip, h_data->GetBinCenter(b), r);
        gr.SetPointError(ip, 0.0, er);
        ++ip;
    }

    if (ip < 5) {
        p.a[0] = 1.0;
        p.valid = true;
        return p;
    }

    TF1 f("f_ratio_pol4", "pol4", 0.0, 70.0);
    TFitResultPtr fit_result = gr.Fit(&f, "SQ0");
    (void)fit_result;

    for (int i = 0; i < 5; ++i) {
        p.a[i] = f.GetParameter(i);
        p.ea[i] = f.GetParError(i);
    }
    p.valid = true;

    mkdir_p(outdir);
    TCanvas c("c_p1_theta_ratio_fit", "", 1000, 750);
    c.SetGrid(1, 1);
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.08);

    TH1D* frame = (TH1D*)gPad->DrawFrame(0.0, RATIO_Y_MIN, 70.0, RATIO_Y_MAX);
    frame->SetTitle((period + "  ep #rightarrow ep#pi_{0}  p_{1} #theta ratio fit").c_str());
    frame->GetXaxis()->SetTitle("p_{1} #theta (deg)");
    frame->GetYaxis()->SetTitle("data / MC");
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);

    gr.SetMarkerStyle(20);
    gr.SetMarkerColor(kBlack);
    gr.SetLineColor(kBlack);
    gr.Draw("PE SAME");
    f.SetLineColor(kRed + 1);
    f.SetLineWidth(2);
    f.Draw("L SAME");

    TLegend leg(0.58, 0.74, 0.90, 0.88);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(1001);
    leg.SetBorderSize(1);
    leg.AddEntry(&gr, "data / MC", "pe");
    leg.AddEntry(&f, "quartic fit", "l");
    leg.Draw();

    c.SaveAs((outdir + "/p1_theta_ratio_quartic_fit.pdf").c_str());
    c.SaveAs((outdir + "/p1_theta_ratio_quartic_fit.png").c_str());

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
            if (!(d > 0.0 && m > 0.0)) continue;
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

static void plot_variable(const std::string& outdir,
                          const std::string& period,
                          const VarConfig& v,
                          const std::vector<TH1D*>& hd,
                          const std::vector<TH1D*>& hm,
                          const Poly4* poly_for_p1_theta) {
    mkdir_p(outdir);

    const int ncols = (v.kind == 2) ? 3 : 4;
    const int nrows = (v.kind == 2) ? 1 : 2;
    const int W = 350 * ncols + 120;
    const int H = 300 * nrows + 120;

    TCanvas c(("c_" + v.key).c_str(), "", W, H);
    c.Divide(ncols, nrows, 0.001, 0.001);

    for (int p = 0; p < v.n_panels; ++p) {
        c.cd(p + 1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.06);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.12);
        gPad->SetGrid(1, 1);

        TH1D* d = hd[p];
        TH1D* m = hm[p];
        d->SetLineColor(kBlack);
        d->SetLineWidth(2);
        m->SetLineColor(kRed + 1);
        m->SetLineWidth(2);
        d->SetMarkerColor(kBlack);
        m->SetMarkerColor(kRed + 1);
        d->SetMarkerStyle(20);
        m->SetMarkerStyle(24);

        double ymax = std::max(d->GetMaximum(), m->GetMaximum());
        if (!(ymax > 0.0)) ymax = 1.0;
        d->SetMaximum(1.3 * ymax);
        d->SetMinimum(0.0);
        d->SetTitle((period + "  " + panel_label(v, p)).c_str());
        d->GetXaxis()->SetTitle(v.x_title.c_str());
        d->GetYaxis()->SetTitle("Counts");
        d->Draw("HIST E");
        m->Draw("HIST E SAME");

        if (p == 0) {
            TLegend* leg = new TLegend(0.55, 0.72, 0.90, 0.88);
            leg->SetFillColor(kWhite);
            leg->SetFillStyle(1001);
            leg->SetBorderSize(1);
            leg->AddEntry(d, "data", "l");
            leg->AddEntry(m, "AAOGEN MC", "l");
            leg->Draw();
        }
    }

    c.SaveAs((outdir + "/" + v.key + "_comparison_linear.pdf").c_str());
    c.SaveAs((outdir + "/" + v.key + "_comparison_linear.png").c_str());

    TCanvas clog(("c_log_" + v.key).c_str(), "", W, H);
    clog.Divide(ncols, nrows, 0.001, 0.001);
    for (int p = 0; p < v.n_panels; ++p) {
        clog.cd(p + 1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.06);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.12);
        gPad->SetGrid(1, 1);
        gPad->SetLogy(1);
        TH1D* d = hd[p];
        TH1D* m = hm[p];
        double ymax = std::max(d->GetMaximum(), m->GetMaximum());
        if (!(ymax > 0.0)) ymax = 1.0;
        d->SetMaximum(30.0 * ymax);
        d->SetMinimum(0.5);
        d->SetTitle((period + "  " + panel_label(v, p)).c_str());
        d->Draw("HIST E");
        m->Draw("HIST E SAME");
        if (p == 0) {
            TLegend* leg = new TLegend(0.55, 0.72, 0.90, 0.88);
            leg->SetFillColor(kWhite);
            leg->SetFillStyle(1001);
            leg->SetBorderSize(1);
            leg->AddEntry(d, "data", "l");
            leg->AddEntry(m, "AAOGEN MC", "l");
            leg->Draw();
        }
    }
    clog.SaveAs((outdir + "/" + v.key + "_comparison_log.pdf").c_str());
    clog.SaveAs((outdir + "/" + v.key + "_comparison_log.png").c_str());

    std::vector<TGraphErrors*> gr = make_ratio_graphs(hd, hm);
    TCanvas cr(("c_ratio_" + v.key).c_str(), "", W, H);
    cr.Divide(ncols, nrows, 0.001, 0.001);

    for (int p = 0; p < v.n_panels; ++p) {
        cr.cd(p + 1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.06);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.12);
        gPad->SetGrid(1, 1);

        double xmin = 0.0, xmax = 1.0;
        range_for_panel(v, p, xmin, xmax);
        TH1D* frame = (TH1D*)gPad->DrawFrame(xmin, RATIO_Y_MIN, xmax, RATIO_Y_MAX);
        frame->SetTitle((period + "  " + panel_label(v, p)).c_str());
        frame->GetXaxis()->SetTitle(v.x_title.c_str());
        frame->GetYaxis()->SetTitle("data / MC");
        gr[p]->Draw("PE SAME");

        if (poly_for_p1_theta && v.key == "p1_theta") {
            TF1 f("f_poly_overlay", "pol4", xmin, xmax);
            for (int i = 0; i < 5; ++i) f.SetParameter(i, poly_for_p1_theta->a[i]);
            f.SetLineColor(kRed + 1);
            f.SetLineWidth(2);
            f.DrawCopy("L SAME");
        }
    }

    cr.SaveAs((outdir + "/" + v.key + "_ratio_linear.pdf").c_str());
    cr.SaveAs((outdir + "/" + v.key + "_ratio_linear.png").c_str());

    for (TGraphErrors* g : gr) delete g;
}

// -----------------------------------------------------------------------------
// Period data/MC normalization histograms
// -----------------------------------------------------------------------------

struct PeriodHistResult {
    std::string period;
    HistSet data_hists;
    HistSet mc_hists;
    Poly4 poly;
    double integrated_ratio = 1.0;
    double integrated_ratio_err = 0.0;

    PeriodHistResult() : data_hists(make_hist_set("data_empty")), mc_hists(make_hist_set("mc_empty")) {}
};

static long long generated_entries(TTree* tree) {
    return tree ? (long long)tree->GetEntries() : 0;
}

static double data_charge_for_tree(TTree* tree,
                                   const std::unordered_map<int, double>& charge_map) {
    if (!tree) return 0.0;
    Branches b;
    b.bind(tree, false);
    if (!b.has_runnum) fatal("[eppi0_norm] FATAL: data tree missing runnum.");

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

static void fill_eppi0_data_hists(const ChannelConfig& cfg,
                                  const PeriodTags& tags,
                                  TTree* tree,
                                  const TopoCutMap& data_cuts,
                                  double current_factor,
                                  HistSet& hists) {
    if (!tree) return;
    Branches b;
    b.bind(tree, false);
    const Long64_t N = tree->GetEntries();
    const double w = 1.0 / current_factor;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        if (!passes_event_selection(cfg, tags, data_cuts, b)) continue;
        fill_hist_set(hists, b, w);
    }
}

static void fill_eppi0_mc_hists(const ChannelConfig& cfg,
                                const PeriodTags& tags,
                                TTree* tree,
                                const TopoCutMap& mc_cuts,
                                double event_norm,
                                double current_factor,
                                HistSet& hists) {
    if (!tree) return;
    Branches b;
    b.bind(tree, true);
    if (!b.has_weight) fatal("[eppi0_norm] FATAL: eppi0 reconstructed MC tree missing weight branch.");

    const Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        if (!passes_event_selection(cfg, tags, mc_cuts, b)) continue;
        const double w = event_norm * b.weight / current_factor;
        fill_hist_set(hists, b, w);
    }
}

static PeriodNormalization run_period_normalization(const std::string& period,
                                                    TTree* data_tree,
                                                    TTree* gen_tree,
                                                    TTree* rec_tree,
                                                    const CSV& csv,
                                                    const std::unordered_map<int, double>& charge_map,
                                                    const TopoCutMap& data_cuts,
                                                    const TopoCutMap& mc_cuts,
                                                    const std::string& output_dir) {
    const ChannelConfig epi = eppi0_config();
    PeriodTags tags;
    tags.display = period;
    tags.period_label = parse_period_from_key(period).period_label;
    tags.period_code = parse_period_from_key(period).period_code;

    const double data_eff = read_current_factor(csv, epi, "exp", period);
    const double mc_eff = read_current_factor(csv, epi, "mc", period);

    const double q_raw = data_charge_for_tree(data_tree, charge_map);
    const double q_mc = q_raw * CHARGE_TO_MC_FACTOR;
    const double lint = RGA_LUMINOSITY_FACTOR_PB_INV_PER_MC * q_mc;
    const double n_gen = (double)generated_entries(gen_tree);
    if (!(n_gen > 0.0)) {
        std::ostringstream ss;
        ss << "[eppi0_norm] FATAL: generated eppi0 MC entries are zero for " << period;
        fatal(ss.str());
    }
    const double norm = lint / n_gen;

    PeriodHistResult hr;
    delete_hist_set(hr.data_hists);
    delete_hist_set(hr.mc_hists);
    hr.period = period;
    hr.data_hists = make_hist_set("data_" + period_dir(period));
    hr.mc_hists = make_hist_set("mc_" + period_dir(period));

    fill_eppi0_data_hists(epi, tags, data_tree, data_cuts, data_eff, hr.data_hists);
    fill_eppi0_mc_hists(epi, tags, rec_tree, mc_cuts, norm, mc_eff, hr.mc_hists);

    const std::string odir = output_dir + "/" + period_dir(period);
    mkdir_p(odir);

    for (const VarConfig& v : variable_configs()) {
        const Poly4* poly_ptr = nullptr;
        if (v.key == "p1_theta") poly_ptr = &hr.poly;
        plot_variable(odir, period, v, hr.data_hists.hists[v.key], hr.mc_hists.hists[v.key], poly_ptr);
    }

    // Fit using all p1_theta panels merged into one 0-70 degree histogram.
    TH1D h_data_fit(("h_data_fit_" + period_dir(period)).c_str(), "", 35, 0.0, 70.0);
    TH1D h_mc_fit(("h_mc_fit_" + period_dir(period)).c_str(), "", 35, 0.0, 70.0);
    h_data_fit.Sumw2();
    h_mc_fit.Sumw2();
    for (TH1D* h : hr.data_hists.hists["p1_theta"]) h_data_fit.Add(h);
    for (TH1D* h : hr.mc_hists.hists["p1_theta"]) h_mc_fit.Add(h);

    Poly4 poly = fit_p1_theta_ratio(period, odir, &h_data_fit, &h_mc_fit);

    // Redraw p1 theta ratio with quartic overlay after fit is known.
    for (const VarConfig& v : variable_configs()) {
        if (v.key == "p1_theta") {
            plot_variable(odir, period, v, hr.data_hists.hists[v.key], hr.mc_hists.hists[v.key], &poly);
        }
    }

    const double data_int = h_data_fit.Integral();
    const double mc_int = h_mc_fit.Integral();
    double data_err2 = 0.0;
    double mc_err2 = 0.0;
    for (int b = 1; b <= h_data_fit.GetNbinsX(); ++b) {
        data_err2 += h_data_fit.GetBinError(b) * h_data_fit.GetBinError(b);
        mc_err2 += h_mc_fit.GetBinError(b) * h_mc_fit.GetBinError(b);
    }

    double ratio = 1.0;
    double ratio_err = 0.0;
    if (data_int > 0.0 && mc_int > 0.0) {
        ratio = data_int / mc_int;
        ratio_err = std::fabs(ratio) * std::sqrt(data_err2 / (data_int * data_int) + mc_err2 / (mc_int * mc_int));
    }

    PeriodNormalization out;
    out.period = period;
    out.poly = poly;
    out.integrated_ratio = ratio;
    out.integrated_ratio_err = ratio_err;

    delete_hist_set(hr.data_hists);
    delete_hist_set(hr.mc_hists);

    std::cout << "[eppi0_norm] " << period
              << " polynomial=" << tuple5(poly)
              << " integrated_ratio=" << ratio
              << " +/- " << ratio_err
              << std::endl;

    return out;
}

static std::string key_for_period_map(const std::map<std::string, TTree*>& trees,
                                      const std::string& period) {
    for (const auto& kv : trees) {
        PeriodTags tags = parse_period_from_key(kv.first);
        if (tags.display == period) return kv.first;
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
                                                  const Poly4& poly,
                                                  std::unordered_map<std::string, RowCounts>& out) {
    if (!tree) return;

    Branches b;
    b.bind(tree, false);
    if (!b.has_helicity) fatal("[eppi0_norm] FATAL: data tree missing helicity.");

    const Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (!passes_event_selection(cfg, tags, data_cuts, b)) continue;

        const std::string topo = topo_dir(b.detector1, b.detector2);
        if (topo.empty()) continue;

        const double theta = b.p1_theta_deg();
        const double ratio = poly.eval(theta);
        if (!(std::isfinite(ratio) && ratio > 0.0)) continue;

        const double weight = 1.0 / (current_factor * ratio);
        const double phi = b.phi_deg();
        const double tabs = b.t_abs();

        for (int r = 0; r < (int)rows.size(); ++r) {
            if (!row_accepts(rows[r], b.x, b.Q2, tabs, phi)) continue;
            HelCounts& hc = out[topo][r];
            if (b.helicity > 0) hc.pos += weight;
            else if (b.helicity < 0) hc.neg += weight;
            else hc.unpol += weight;
        }
    }
}

static void write_normalized_counts_to_csv(CSV& csv,
                                           const ChannelConfig& cfg,
                                           const std::string& period,
                                           const std::unordered_map<std::string, RowCounts>& topo_counts) {
    for (const auto& kt : topo_counts) {
        const std::string topo_label = topo_label_for_csv(kt.first);
        if (topo_label.empty()) continue;

        const int c_unpol = col_strict(csv, normalized_raw_yield_col(cfg, topo_label, period, "unpol"));
        const int c_pos = col_strict(csv, normalized_raw_yield_col(cfg, topo_label, period, "pos"));
        const int c_neg = col_strict(csv, normalized_raw_yield_col(cfg, topo_label, period, "neg"));

        for (const auto& kr : kt.second) {
            const int row = kr.first;
            if (row < 0 || row >= (int)csv.rows.size()) fatal("[eppi0_norm] FATAL: row index out of range.");
            const HelCounts& h = kr.second;
            const double unpol = h.pos + h.neg + h.unpol;
            std::ostringstream su, sp, sn;
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
        if (n.period == period) return n;
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
        Poly4 poly;
    };

    std::vector<WorkItem> items;
    const ChannelConfig dvcs = dvcs_config();
    const ChannelConfig epi = eppi0_config();

    auto add_items = [&](const ChannelConfig& cfg, const std::map<std::string, TTree*>& trees) {
        for (const auto& kv : trees) {
            PeriodTags tags = parse_period_from_key(kv.first);
            if (tags.display == "Fa18 Inb Supp") continue;
            WorkItem w;
            w.cfg = cfg;
            w.tags = tags;
            w.tree = kv.second;
            w.current_factor = read_current_factor(csv, cfg, "exp", tags.display);
            w.poly = find_norm(norms, tags.display).poly;
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
                                             items[i].poly,
                                             results[i]);
    }

    for (int i = 0; i < (int)items.size(); ++i) {
        write_normalized_counts_to_csv(csv,
                                       items[i].cfg,
                                       items[i].tags.display,
                                       results[i]);
    }
}

static void write_norms_to_csv(CSV& csv, const std::vector<PeriodNormalization>& norms) {
    for (const std::string& period : CSV_PERIOD_ORDER) {
        const PeriodNormalization& n = find_norm(norms, period);
        const int c_poly = col_strict(csv, norm_poly_col(period));
        const int c_factor = col_strict(csv, norm_factor_col(period));
        const std::string poly_cell = tuple5(n.poly);
        const std::string factor_cell = tuple2(n.integrated_ratio, n.integrated_ratio_err);
        for (auto& row : csv.rows) {
            row[c_poly] = poly_cell;
            row[c_factor] = factor_cell;
        }
    }
}

static std::vector<PeriodNormalization> unity_norms() {
    std::vector<PeriodNormalization> norms;
    for (const std::string& period : CSV_PERIOD_ORDER) {
        PeriodNormalization n;
        n.period = period;
        n.poly.valid = true;
        n.poly.a[0] = 1.0;
        n.integrated_ratio = 1.0;
        n.integrated_ratio_err = 0.0;
        norms.push_back(n);
    }
    return norms;
}

} // namespace

bool update_eppi0_normalization_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
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
            std::cout << "[eppi0_norm] Override enabled: using unity eppi0 normalization polynomial." << std::endl;
            norms = unity_norms();
        } else {
            mkdir_p(options.output_dir);
            const std::unordered_map<int, double> charge_map = read_charge_csv(options.charge_csv_path);

            for (const std::string& period : WORK_PERIOD_ORDER) {
                TTree* data_tree = tree_for_period(eppi0DataTrees, period, "eppi0 data");
                TTree* gen_tree = tree_for_period(eppi0GenMcTrees, period, "eppi0 generated MC");
                TTree* rec_tree = tree_for_period(eppi0RecMcTrees, period, "eppi0 reconstructed MC");

                PeriodNormalization n = run_period_normalization(period,
                                                                 data_tree,
                                                                 gen_tree,
                                                                 rec_tree,
                                                                 csv,
                                                                 charge_map,
                                                                 data_cuts,
                                                                 mc_cuts,
                                                                 options.output_dir);
                norms.push_back(n);
            }
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
