#include "pi0_subtracted_kinematics.h"

#include "exclusivity_cuts.h"
#include "global_cuts.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <TSystem.h>
#include <TString.h>
#include <TDataType.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
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

static constexpr double kPi = 3.14159265358979323846;
static constexpr double kRad2Deg = 180.0 / kPi;

struct PeriodDef {
    std::string pretty;
    std::string period_label;
    std::string cut_period_code;
    std::string data_key;
    std::string mc_key;
};

struct PlotVar {
    std::string branch;
    std::string title;
    std::string x_label;
    int nbins = 100;
    double xmin = 0.0;
    double xmax = 1.0;
    bool radians_to_degrees = false;
};

enum class DataShapeMode {
    RawPostExclusivity,
    Pi0Subtracted
};

struct PanelSpec {
    std::string label;
    PeriodDef period;
    int detector1_filter = -1;
    int detector2_filter = -1;
    int electron_sector_filter = -1;
};

struct WeightSummary {
    std::string period;
    int n_valid_bins = 0;
    int n_positive_weight_bins = 0;
    int n_zero_weight_bins = 0;
    int n_negative_signal_bins = 0;
    double total_csv_raw_epg = 0.0;
    double total_csv_signal = 0.0;
    double sum_weight = 0.0;
    double sum_weight2 = 0.0;
    double max_weight = 0.0;
};

struct FillSummary {
    std::string mode;
    std::string category;
    std::string panel;
    std::string period;
    std::string sample;
    long long entries = 0;
    long long pass_global = 0;
    long long pass_selector = 0;
    long long pass_exclusivity = 0;
    long long filled_events = 0;
    double filled_weight_sum = 0.0;
};

struct RowBin {
    double xBmin = 0.0;
    double xBmax = 0.0;
    double Q2min = 0.0;
    double Q2max = 0.0;
    double tmin = 0.0;
    double tmax = 0.0;
    double phimin = 0.0;
    double phimax = 0.0;
    bool valid = false;
};

struct CsvInfo {
    std::vector<std::string> header;
    std::unordered_map<std::string, int> col;
    std::vector<std::vector<std::string> > rows;
    std::vector<RowBin> bins;
    std::map<std::string, std::vector<double> > data_weight_by_period_row;
    std::map<std::string, WeightSummary> weight_summary_by_period;
};

static std::mutex g_root_mutex;

static void fatal(const std::string& msg) {
    throw std::runtime_error(msg);
}

static std::string sanitize(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_') out.push_back(c);
        else out.push_back('_');
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

static double parse_double_or_nan(const std::string& s) {
    if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
    char* endp = nullptr;
    const double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return std::numeric_limits<double>::quiet_NaN();
    return v;
}

static bool parse_tuple_numbers(const std::string& cell, std::vector<double>& out) {
    out.clear();
    std::string s;
    s.reserve(cell.size());
    for (char c : cell) {
        if (c == '(' || c == ')' || c == '[' || c == ']' || c == '"') s.push_back(' ');
        else s.push_back(c);
    }

    std::stringstream ss(s);
    std::string part;
    while (std::getline(ss, part, ',')) {
        const double v = parse_double_or_nan(part);
        if (std::isfinite(v)) out.push_back(v);
    }

    return !out.empty();
}

static double first_tuple_value_or_default(const std::string& cell, double def) {
    std::vector<double> vals;
    if (!parse_tuple_numbers(cell, vals)) return def;
    if (vals.empty() || !std::isfinite(vals[0])) return def;
    return vals[0];
}

static bool is_true_valid(const std::string& s) {
    return s == "1" || s == "1.0" || s == "true" || s == "TRUE" || s == "True";
}

static int col_strict(const CsvInfo& csv, const std::string& name) {
    auto it = csv.col.find(name);
    if (it == csv.col.end()) {
        std::ostringstream ss;
        ss << "[pi0_subtracted_kinematics] FATAL: missing required CSV column: " << name;
        fatal(ss.str());
    }
    return it->second;
}

static int col_optional(const CsvInfo& csv, const std::string& name) {
    auto it = csv.col.find(name);
    if (it == csv.col.end()) return -1;
    return it->second;
}

static std::string raw_yield_topo_col(const std::string& topo,
                                      const std::string& period) {
    return "raw yield, ep->epg, " + topo + ", exp, " + period + ", unpol";
}

static std::string signal_yield_col(const std::string& period) {
    return "signal yield, ep->epg, exp, " + period + ", unpol";
}

static const std::vector<std::string>& topo_labels() {
    static const std::vector<std::string> labels = {
        "(FD, FD)",
        "(CD, FD)",
        "(CD, FT)"
    };
    return labels;
}

static std::string topo_dir_name(int detector1, int detector2) {
    if (detector1 == 1 && detector2 == 1) return "FD_FD";
    if (detector1 == 2 && detector2 == 1) return "CD_FD";
    if (detector1 == 2 && detector2 == 0) return "CD_FT";
    return "unknown_topology";
}

static std::string data_mode_name(DataShapeMode mode) {
    if (mode == DataShapeMode::RawPostExclusivity) return "raw_post_exclusivity";
    return "pi0_subtracted";
}

static std::string data_mode_title(DataShapeMode mode) {
    if (mode == DataShapeMode::RawPostExclusivity) return "raw post-exclusivity";
    return "#pi^{0}-subtracted";
}

static double wrap_phi_deg(double phi_deg) {
    double p = std::fmod(phi_deg, 360.0);
    if (p < 0.0) p += 360.0;
    if (p >= 360.0) p = std::nextafter(360.0, 0.0);
    return p;
}

static bool in_range(double v, double a, double b) {
    return v >= a && v < b;
}

static bool row_accepts_phi(double phi_deg, double pmin, double pmax) {
    if (pmax > pmin) return in_range(phi_deg, pmin, pmax);
    return phi_deg >= pmin || phi_deg < pmax;
}

static int find_csv_row_for_event(const CsvInfo& csv,
                                  double xB,
                                  double Q2,
                                  double t_abs,
                                  double phi_deg) {
    for (int i = 0; i < (int)csv.bins.size(); ++i) {
        const RowBin& r = csv.bins[(size_t)i];
        if (!r.valid) continue;
        if (!in_range(xB, r.xBmin, r.xBmax)) continue;
        if (!in_range(Q2, r.Q2min, r.Q2max)) continue;
        if (!in_range(t_abs, r.tmin, r.tmax)) continue;
        if (!row_accepts_phi(phi_deg, r.phimin, r.phimax)) continue;
        return i;
    }
    return -1;
}

static std::vector<PeriodDef> periods() {
    return {
        {"Sp18 Inb", "sp18_inb", "Sp18_Inb", "DVCS_Sp18_inb", "DVCS_Sp18_inb_rec"},
        {"Sp18 Out", "sp18_out", "Sp18_Out", "DVCS_Sp18_out", "DVCS_Sp18_out_rec"},
        {"Fa18 Inb", "fa18_inb", "Fa18_Inb", "DVCS_Fa18_inb", "DVCS_Fa18_inb_rec"},
        {"Fa18 Out", "fa18_out", "Fa18_Out", "DVCS_Fa18_out", "DVCS_Fa18_out_rec"},
        {"Sp19 Inb", "sp19_inb", "Sp19_Inb", "DVCS_Sp19_inb", "DVCS_Sp19_inb_rec"}
    };
}

static std::vector<PlotVar> plot_vars() {
    return {
        {"e_p",       "Electron momentum",        "e p (GeV)",                 100, 1.0, 9.0,  false},
        {"e_theta",   "Electron polar angle",     "e #theta (deg)",            100, 0.0, 45.0, true},
        {"p1_p",      "Proton momentum",          "p p (GeV)",                 100, 0.0, 3.0,  false},
        {"p1_theta",  "Proton polar angle",       "p #theta (deg)",            100, 0.0, 90.0, true},
        {"p2_p",      "Photon momentum",          "#gamma p (GeV)",            100, 1.0, 9.0,  false},
        {"p2_theta",  "Photon polar angle",       "#gamma #theta (deg)",       100, 0.0, 45.0, true}
    };
}

static CsvInfo load_csv_info(const std::string& path) {
    CsvInfo csv;
    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[pi0_subtracted_kinematics] FATAL: cannot open CSV: " << path;
        fatal(ss.str());
    }

    std::string line;
    if (!std::getline(fin, line)) {
        fatal("[pi0_subtracted_kinematics] FATAL: empty CSV.");
    }

    csv.header = split_csv_line(line);
    for (int i = 0; i < (int)csv.header.size(); ++i) {
        if (csv.col.find(csv.header[(size_t)i]) != csv.col.end()) {
            std::ostringstream ss;
            ss << "[pi0_subtracted_kinematics] FATAL: duplicate CSV column: "
               << csv.header[(size_t)i];
            fatal(ss.str());
        }
        csv.col[csv.header[(size_t)i]] = i;
    }

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = split_csv_line(line);
        if (row.size() < csv.header.size()) row.resize(csv.header.size(), "");
        if (row.size() != csv.header.size()) {
            fatal("[pi0_subtracted_kinematics] FATAL: CSV row width mismatch.");
        }
        csv.rows.push_back(std::move(row));
    }

    const int c_xmin = col_strict(csv, "xBmin");
    const int c_xmax = col_strict(csv, "xBmax");
    const int c_qmin = col_strict(csv, "Q2min");
    const int c_qmax = col_strict(csv, "Q2max");
    const int c_tmin = col_strict(csv, "t_abs_min");
    const int c_tmax = col_strict(csv, "t_abs_max");
    const int c_pmin = col_strict(csv, "phimin");
    const int c_pmax = col_strict(csv, "phimax");
    const int c_valid = col_strict(csv, "valid bin");

    csv.bins.reserve(csv.rows.size());
    for (const auto& row : csv.rows) {
        RowBin b;
        b.xBmin = parse_double_or_nan(row[(size_t)c_xmin]);
        b.xBmax = parse_double_or_nan(row[(size_t)c_xmax]);
        b.Q2min = parse_double_or_nan(row[(size_t)c_qmin]);
        b.Q2max = parse_double_or_nan(row[(size_t)c_qmax]);
        b.tmin = parse_double_or_nan(row[(size_t)c_tmin]);
        b.tmax = parse_double_or_nan(row[(size_t)c_tmax]);
        b.phimin = parse_double_or_nan(row[(size_t)c_pmin]);
        b.phimax = parse_double_or_nan(row[(size_t)c_pmax]);
        b.valid = is_true_valid(row[(size_t)c_valid]);
        csv.bins.push_back(b);
    }

    for (const PeriodDef& p : periods()) {
        const int c_sig = col_strict(csv, signal_yield_col(p.pretty));
        std::vector<int> raw_cols;
        for (const std::string& topo : topo_labels()) {
            const int c = col_optional(csv, raw_yield_topo_col(topo, p.pretty));
            if (c >= 0) raw_cols.push_back(c);
        }

        if (raw_cols.empty()) {
            std::ostringstream ss;
            ss << "[pi0_subtracted_kinematics] FATAL: no raw epgamma DATA yield columns found for period "
               << p.pretty;
            fatal(ss.str());
        }

        std::vector<double> weights(csv.rows.size(), 0.0);
        WeightSummary summary;
        summary.period = p.pretty;

        for (size_t ir = 0; ir < csv.rows.size(); ++ir) {
            const auto& row = csv.rows[ir];
            if (!csv.bins[ir].valid) continue;

            summary.n_valid_bins += 1;

            const double signal = first_tuple_value_or_default(row[(size_t)c_sig], 0.0);
            double raw_sum = 0.0;
            for (int c_raw : raw_cols) {
                const double v = first_tuple_value_or_default(row[(size_t)c_raw], 0.0);
                if (std::isfinite(v) && v > 0.0) raw_sum += v;
            }

            if (std::isfinite(raw_sum) && raw_sum > 0.0) {
                summary.total_csv_raw_epg += raw_sum;
            }
            if (std::isfinite(signal)) {
                summary.total_csv_signal += signal;
                if (signal < 0.0) summary.n_negative_signal_bins += 1;
            }

            if (std::isfinite(signal) && signal > 0.0 && raw_sum > 0.0) {
                const double w = signal / raw_sum;
                weights[ir] = w;
                summary.n_positive_weight_bins += 1;
                summary.sum_weight += w;
                summary.sum_weight2 += w * w;
                summary.max_weight = std::max(summary.max_weight, w);
            } else {
                summary.n_zero_weight_bins += 1;
            }
        }

        csv.weight_summary_by_period[p.pretty] = summary;
        csv.data_weight_by_period_row[p.pretty] = std::move(weights);
    }

    std::cout << "[pi0_subtracted_kinematics] Loaded CSV with "
              << csv.rows.size() << " rows and pi0-subtracted signal weights from "
              << path << std::endl;

    return csv;
}

static std::string topo_key(int detector1, int detector2) {
    if (detector1 == 1 && detector2 == 1) return "FD_FD";
    if (detector1 == 2 && detector2 == 1) return "CD_FD";
    if (detector1 == 2 && detector2 == 0) return "CD_FT";
    return "";
}

static bool value_within_cut(double v, const Stats& s) {
    if (!std::isfinite(v)) return false;

    if (std::isfinite(s.cut_low) && std::isfinite(s.cut_high) && s.cut_high > s.cut_low) {
        return v >= s.cut_low && v <= s.cut_high;
    }

    if (s.mode == "upper_quantile") {
        if (!std::isfinite(s.cut_high)) return true;
        return v <= s.cut_high;
    }

    if (!(std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0)) {
        return true;
    }

    return v >= s.mean - 3.0 * s.std && v <= s.mean + 3.0 * s.std;
}

static std::map<std::string, CutDict> load_combined_cuts_json(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[pi0_subtracted_kinematics] FATAL: cannot open combined cuts JSON: " << path;
        fatal(ss.str());
    }

    nlohmann::json j;
    fin >> j;

    std::map<std::string, CutDict> out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        if (!it.value().is_object()) continue;
        CutDict cd;

        auto fill_sample = [&](const char* sample, std::map<std::string, Stats>& target) {
            if (!it.value().contains(sample) || !it.value()[sample].is_object()) return;
            const auto& block = it.value()[sample];
            for (auto vit = block.begin(); vit != block.end(); ++vit) {
                if (!vit.value().is_object()) continue;
                Stats s;
                s.mean = std::numeric_limits<double>::quiet_NaN();
                s.std = std::numeric_limits<double>::quiet_NaN();
                s.cut_low = std::numeric_limits<double>::quiet_NaN();
                s.cut_high = std::numeric_limits<double>::quiet_NaN();
                s.quantile = 0.0;
                s.mode = "symmetric_3sigma";
                if (vit.value().contains("mean")) s.mean = vit.value()["mean"].get<double>();
                if (vit.value().contains("std")) s.std = vit.value()["std"].get<double>();
                if (vit.value().contains("cut_low")) s.cut_low = vit.value()["cut_low"].get<double>();
                if (vit.value().contains("cut_high")) s.cut_high = vit.value()["cut_high"].get<double>();
                if (vit.value().contains("quantile")) s.quantile = vit.value()["quantile"].get<double>();
                if (vit.value().contains("mode")) s.mode = vit.value()["mode"].get<std::string>();
                target[vit.key()] = s;
            }
        };

        fill_sample("data", cd.data);
        fill_sample("mc", cd.mc);
        out[it.key()] = cd;
    }

    std::cout << "[pi0_subtracted_kinematics] Loaded " << out.size()
              << " combined-cut entries from " << path << std::endl;

    return out;
}

struct Branches {
    int runnum = 0; bool has_runnum = false;
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    double x = 0.0; bool has_x = false;
    double Q2 = 0.0; bool has_Q2 = false;
    double phi2 = 0.0; bool has_phi2 = false;

    double t1 = 0.0; bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open_angle_ep2 = false;
    double pTmiss = 0.0; bool has_pTmiss = false;
    double Delta_phi = 0.0; bool has_Delta_phi = false;
    double Emiss2 = 0.0; bool has_Emiss2 = false;
    double Mx2 = 0.0; bool has_Mx2 = false;
    double Mx2_1 = 0.0; bool has_Mx2_1 = false;
    double Mx2_2 = 0.0; bool has_Mx2_2 = false;
    double xF = 0.0; bool has_xF = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;

    double e_p = 0.0; bool has_e_p = false;
    double e_theta = 0.0; bool has_e_theta = false;
    double e_phi = 0.0; bool has_e_phi = false;

    double p1_p = 0.0; bool has_p1_p = false;
    double p1_theta = 0.0; bool has_p1_theta = false;
    double p1_phi = 0.0; bool has_p1_phi = false;

    double p2_p = 0.0; bool has_p2_p = false;
    double p2_theta = 0.0; bool has_p2_theta = false;
    double p2_phi = 0.0; bool has_p2_phi = false;

    void bind(TTree* t) {
        auto ena = [&](const char* n) {
            if (t->GetBranch(n)) t->SetBranchStatus(n, 1);
        };

        ena("runnum"); ena("detector1"); ena("detector2");
        ena("x"); ena("Q2"); ena("phi2");
        ena("t1"); ena("open_angle_ep2"); ena("pTmiss"); ena("Delta_phi");
        ena("Emiss2"); ena("Mx2"); ena("Mx2_1"); ena("Mx2_2"); ena("xF");
        ena("theta_gamma_gamma");
        ena("e_p"); ena("e_theta"); ena("e_phi");
        ena("p1_p"); ena("p1_theta"); ena("p1_phi");
        ena("p2_p"); ena("p2_theta"); ena("p2_phi");

        auto bI = [&](const char* n, int* a, bool& f) {
            if (t->GetBranch(n)) {
                t->SetBranchAddress(n, a, nullptr, nullptr, kInt_t, false);
                f = true;
            }
        };

        auto bD = [&](const char* n, double* a, bool& f) {
            if (t->GetBranch(n)) {
                t->SetBranchAddress(n, a, nullptr, nullptr, kDouble_t, false);
                f = true;
            }
        };

        bI("runnum", &runnum, has_runnum);
        bI("detector1", &detector1, has_detector1);
        bI("detector2", &detector2, has_detector2);

        bD("x", &x, has_x);
        bD("Q2", &Q2, has_Q2);
        bD("phi2", &phi2, has_phi2);
        bD("t1", &t1, has_t1);
        bD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bD("pTmiss", &pTmiss, has_pTmiss);
        bD("Delta_phi", &Delta_phi, has_Delta_phi);
        bD("Emiss2", &Emiss2, has_Emiss2);
        bD("Mx2", &Mx2, has_Mx2);
        bD("Mx2_1", &Mx2_1, has_Mx2_1);
        bD("Mx2_2", &Mx2_2, has_Mx2_2);
        bD("xF", &xF, has_xF);
        bD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);

        bD("e_p", &e_p, has_e_p);
        bD("e_theta", &e_theta, has_e_theta);
        bD("e_phi", &e_phi, has_e_phi);
        bD("p1_p", &p1_p, has_p1_p);
        bD("p1_theta", &p1_theta, has_p1_theta);
        bD("p1_phi", &p1_phi, has_p1_phi);
        bD("p2_p", &p2_p, has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi", &p2_phi, has_p2_phi);
    }

    double plot_value(const PlotVar& v, bool& ok) const {
        ok = true;
        double value = 0.0;

        if (v.branch == "e_p") { if (!has_e_p) ok = false; value = e_p; }
        else if (v.branch == "e_theta") { if (!has_e_theta) ok = false; value = e_theta; }
        else if (v.branch == "p1_p") { if (!has_p1_p) ok = false; value = p1_p; }
        else if (v.branch == "p1_theta") { if (!has_p1_theta) ok = false; value = p1_theta; }
        else if (v.branch == "p2_p") { if (!has_p2_p) ok = false; value = p2_p; }
        else if (v.branch == "p2_theta") { if (!has_p2_theta) ok = false; value = p2_theta; }
        else ok = false;

        if (!ok || !std::isfinite(value)) return 0.0;
        return v.radians_to_degrees ? value * kRad2Deg : value;
    }

    std::map<std::string, double> exclusivity_values() const {
        std::map<std::string, double> m;
        if (has_Delta_phi) m["Delta_phi"] = Delta_phi;
        if (has_pTmiss) m["pTmiss"] = pTmiss;
        if (has_Emiss2) m["Emiss2"] = Emiss2;
        if (has_Mx2) m["Mx2"] = Mx2;
        if (has_Mx2_1) m["Mx2_1"] = Mx2_1;
        if (has_Mx2_2) m["Mx2_2"] = Mx2_2;
        if (has_xF) m["xF"] = xF;
        if (has_theta_gamma_gamma) m["theta_gamma_gamma"] = theta_gamma_gamma;
        return m;
    }
};

static bool passes_global_cuts_for_event(const Branches& b,
                                         const PeriodDef& period) {
    const GlobalCutConfig& cfg = default_global_cuts();

    if (b.has_runnum && is_excluded_run(b.runnum, cfg)) return false;
    if (!(b.has_t1 && b.has_open_angle_ep2)) return false;
    if (cfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (!(b.has_detector1 && b.has_detector2)) return false;

    if (global_cuts_require_sector_phi(cfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            fatal("[pi0_subtracted_kinematics] FATAL: sector cuts require e_phi, p1_phi, and p2_phi.");
        }
    }

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi && b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[pi0_subtracted_kinematics] FATAL: dvcsgen ycol cut requires e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi.");
        }

        if (global_cuts_require_sector_phi(cfg)) {
            return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                      b.detector1, b.detector2,
                                      period.period_label,
                                      b.e_p, b.e_theta, b.e_phi,
                                      b.p1_phi,
                                      b.p2_p, b.p2_theta, b.p2_phi,
                                      cfg);
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period.period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    if (global_cuts_require_sector_phi(cfg)) {
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  b.e_phi, b.p1_phi, b.p2_phi,
                                  cfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                              b.detector1, b.detector2,
                              cfg);
}

static bool passes_exclusivity_cuts(const Branches& b,
                                    const PeriodDef& period,
                                    const std::map<std::string, CutDict>& cuts,
                                    bool use_data_cuts) {
    if (!(b.has_detector1 && b.has_detector2)) return false;
    const std::string topo = topo_key(b.detector1, b.detector2);
    if (topo.empty()) return false;

    const std::string key = "DVCS_" + period.cut_period_code + "_" + topo;
    auto it = cuts.find(key);
    if (it == cuts.end()) {
        std::ostringstream ss;
        ss << "[pi0_subtracted_kinematics] FATAL: missing combined-cut key: " << key;
        fatal(ss.str());
    }

    const std::map<std::string, Stats>& cut_map = use_data_cuts ? it->second.data : it->second.mc;
    const std::map<std::string, double> vals = b.exclusivity_values();

    for (const auto& kv : cut_map) {
        auto vit = vals.find(kv.first);
        if (vit == vals.end()) {
            std::ostringstream ss;
            ss << "[pi0_subtracted_kinematics] FATAL: missing branch needed for exclusivity cut: "
               << kv.first;
            fatal(ss.str());
        }
        if (!value_within_cut(vit->second, kv.second)) return false;
    }

    return true;
}

static double data_weight_for_event(const CsvInfo& csv,
                                    const PeriodDef& period,
                                    const Branches& b) {
    if (!(b.has_x && b.has_Q2 && b.has_t1 && b.has_phi2)) {
        fatal("[pi0_subtracted_kinematics] FATAL: missing x/Q2/t1/phi2 for CSV-bin lookup.");
    }

    const int row = find_csv_row_for_event(csv,
                                           b.x,
                                           b.Q2,
                                           std::fabs(b.t1),
                                           wrap_phi_deg(b.phi2 * kRad2Deg));
    if (row < 0) return 0.0;

    auto it = csv.data_weight_by_period_row.find(period.pretty);
    if (it == csv.data_weight_by_period_row.end()) return 0.0;
    if (row >= (int)it->second.size()) return 0.0;

    const double w = it->second[(size_t)row];
    if (!std::isfinite(w) || w <= 0.0) return 0.0;
    return w;
}

static void require_trees(const std::vector<PeriodDef>& ps,
                          const std::map<std::string, TTree*>& dataTrees,
                          const std::map<std::string, TTree*>& mcTrees) {
    std::vector<std::string> missing;
    for (const PeriodDef& p : ps) {
        auto itd = dataTrees.find(p.data_key);
        if (itd == dataTrees.end() || itd->second == nullptr) missing.push_back("data:" + p.data_key);
        auto itm = mcTrees.find(p.mc_key);
        if (itm == mcTrees.end() || itm->second == nullptr) missing.push_back("mc:" + p.mc_key);
    }

    if (!missing.empty()) {
        std::ostringstream ss;
        ss << "[pi0_subtracted_kinematics] FATAL: missing required trees: ";
        for (size_t i = 0; i < missing.size(); ++i) {
            if (i) ss << ", ";
            ss << missing[i];
        }
        fatal(ss.str());
    }
}

static bool passes_panel_selector(const Branches& b, const PanelSpec& panel) {
    if (panel.detector1_filter >= 0) {
        if (!b.has_detector1 || b.detector1 != panel.detector1_filter) return false;
    }
    if (panel.detector2_filter >= 0) {
        if (!b.has_detector2 || b.detector2 != panel.detector2_filter) return false;
    }
    if (panel.electron_sector_filter >= 0) {
        if (!b.has_e_phi) return false;
        return fd_sector_from_phi_rad(b.e_phi) == panel.electron_sector_filter;
    }
    return true;
}

static double event_weight_for_mode(const CsvInfo& csv,
                                    const PeriodDef& period,
                                    const Branches& b,
                                    DataShapeMode mode) {
    if (mode == DataShapeMode::RawPostExclusivity) {
        return 1.0;
    }
    return data_weight_for_event(csv, period, b);
}

static FillSummary fill_tree_panel(TTree* tree,
                                   const PanelSpec& panel,
                                   const std::vector<PlotVar>& vars,
                                   const std::map<std::string, CutDict>& cuts,
                                   const CsvInfo& csv,
                                   bool is_data,
                                   DataShapeMode mode,
                                   const std::string& category,
                                   std::vector<TH1D*>& hists) {
    if (!tree) fatal("[pi0_subtracted_kinematics] FATAL: null TTree.");

    FillSummary summary;
    summary.mode = data_mode_name(mode);
    summary.category = category;
    summary.panel = panel.label;
    summary.period = panel.period.pretty;
    summary.sample = is_data ? "data" : "mc";

    {
        std::lock_guard<std::mutex> lock(g_root_mutex);
        tree->SetBranchStatus("*", 0);
    }

    Branches b;
    {
        std::lock_guard<std::mutex> lock(g_root_mutex);
        b.bind(tree);
    }

    const Long64_t nentries = tree->GetEntries();
    summary.entries = (long long)nentries;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        if (!passes_global_cuts_for_event(b, panel.period)) continue;
        summary.pass_global += 1;

        if (!passes_panel_selector(b, panel)) continue;
        summary.pass_selector += 1;

        if (!passes_exclusivity_cuts(b, panel.period, cuts, is_data)) continue;
        summary.pass_exclusivity += 1;

        double weight = 1.0;
        if (is_data) {
            weight = event_weight_for_mode(csv, panel.period, b, mode);
            if (!(std::isfinite(weight) && weight > 0.0)) continue;
        }

        summary.filled_events += 1;
        summary.filled_weight_sum += weight;

        for (size_t iv = 0; iv < vars.size(); ++iv) {
            bool ok = false;
            const double value = b.plot_value(vars[iv], ok);
            if (!ok || !std::isfinite(value)) continue;
            hists[iv]->Fill(value, weight);
        }
    }

    {
        std::lock_guard<std::mutex> lock(g_root_mutex);
        tree->SetBranchStatus("*", 1);
    }

    std::cout << "[pi0_subtracted_kinematics] "
              << summary.sample
              << " mode=" << summary.mode
              << " category=" << category
              << " panel=" << panel.label
              << " period=" << panel.period.pretty
              << " entries=" << summary.entries
              << " pass_global=" << summary.pass_global
              << " pass_selector=" << summary.pass_selector
              << " pass_exclusivity=" << summary.pass_exclusivity
              << " filled=" << summary.filled_events
              << " weight_sum=" << summary.filled_weight_sum
              << std::endl;

    return summary;
}

static void normalize_hist(TH1D* h) {
    if (!h) return;
    const double integral = h->Integral();
    if (integral > 0.0) h->Scale(1.0 / integral);
}

static void style_hist(TH1D* h, int color, int marker) {
    if (!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.80);
    h->SetLineWidth(2);
    h->SetStats(0);
}

static void draw_info_pad(const PlotVar& var,
                          DataShapeMode mode,
                          const std::string& category_title) {
    gPad->SetLeftMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->DrawFrame(0.0, 0.0, 1.0, 1.0);

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextAlign(22);
    lat.SetTextSize(0.050);
    lat.DrawLatex(0.50, 0.78, ("DVCS shape: " + var.title).c_str());
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.50, 0.67, category_title.c_str());

    lat.SetTextSize(0.031);
    if (mode == DataShapeMode::Pi0Subtracted) {
        lat.DrawLatex(0.50, 0.56, "DATA: accepted ep#gamma events weighted to CSV #pi^{0}-subtracted signal yield per 4D bin");
        lat.DrawLatex(0.50, 0.47, "This is a bin-weighted signal-shape approximation, not event-level background tagging");
    } else {
        lat.DrawLatex(0.50, 0.56, "DATA: raw accepted ep#gamma candidates after global and topology-dependent exclusivity cuts");
        lat.DrawLatex(0.50, 0.47, "No #pi^{0}-subtraction weights are applied in this diagnostic");
    }
    lat.DrawLatex(0.50, 0.38, "MC: reconstructed DVCS MC with the same global and topology-dependent exclusivity cuts");
    lat.DrawLatex(0.50, 0.29, "Each DATA and MC histogram is normalized to its own integral");
    lat.DrawLatex(0.50, 0.20, "Angles are converted from radians to degrees");

    std::ostringstream range_ss;
    range_ss << "Range: [" << var.xmin << ", " << var.xmax << "], bins=" << var.nbins;
    lat.DrawLatex(0.50, 0.11, range_ss.str().c_str());
}

static void save_canvas_panels(const PlotVar& var,
                               const std::vector<PanelSpec>& panels,
                               const std::vector<TH1D*>& data_hists,
                               const std::vector<TH1D*>& mc_hists,
                               DataShapeMode mode,
                               const std::string& category_title,
                               const std::string& out_dir) {
    double ymax = 0.0;
    for (size_t ip = 0; ip < panels.size(); ++ip) {
        normalize_hist(data_hists[ip]);
        normalize_hist(mc_hists[ip]);
        style_hist(data_hists[ip], kBlue + 1, 20);
        style_hist(mc_hists[ip], kRed + 1, 24);
        ymax = std::max(ymax, data_hists[ip]->GetMaximum());
        ymax = std::max(ymax, mc_hists[ip]->GetMaximum());
    }
    if (!(ymax > 0.0)) ymax = 1.0;

    std::filesystem::create_directories(out_dir);

    TCanvas c(("c_" + data_mode_name(mode) + "_" + sanitize(category_title) + "_" + sanitize(var.branch)).c_str(),
              "",
              2100,
              1200);
    c.Divide(3, 2, 0.002, 0.002);

    const int npanels_to_draw = std::min((int)panels.size(), 6);
    for (int ip = 0; ip < npanels_to_draw; ++ip) {
        c.cd(ip + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.12);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetGrid(1, 1);

        TH1D* hd = data_hists[(size_t)ip];
        TH1D* hm = mc_hists[(size_t)ip];
        hd->SetTitle("");
        hd->SetMaximum(1.25 * ymax);
        hd->GetXaxis()->SetTitle(var.x_label.c_str());
        hd->GetYaxis()->SetTitle("Normalized counts");
        hd->GetXaxis()->CenterTitle(true);
        hd->GetYaxis()->CenterTitle(true);
        hd->GetYaxis()->SetTitleOffset(1.45);

        hd->Draw("HIST");
        hm->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.50, 0.68, 0.92, 0.84);
        leg->SetFillStyle(1001);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->SetTextFont(42);
        leg->SetTextSize(0.034);
        leg->AddEntry(hd, (std::string("Data, ") + data_mode_title(mode)).c_str(), "l");
        leg->AddEntry(hm, "Reconstructed DVCS MC", "l");
        leg->Draw();

        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextFont(42);
        lat.SetTextSize(0.052);
        lat.SetTextAlign(13);
        lat.DrawLatex(0.18, 0.92, panels[(size_t)ip].label.c_str());
    }

    if (panels.size() < 6) {
        c.cd(6);
        draw_info_pad(var, mode, category_title);
    }

    const std::string out_path = out_dir + "/" + sanitize(var.branch) + "_" + data_mode_name(mode) + "_data_vs_rec_mc.png";
    c.SaveAs(out_path.c_str());
}

static TH1D* make_hist(const std::string& name, const PlotVar& v) {
    TH1D* h = new TH1D(name.c_str(), "", v.nbins, v.xmin, v.xmax);
    h->SetDirectory(nullptr);
    h->Sumw2();
    return h;
}

static void write_weight_diagnostics_csv(const CsvInfo& csv, const std::string& diagnostics_dir) {
    std::filesystem::create_directories(diagnostics_dir);
    const std::string path = diagnostics_dir + "/data_pi0_weight_summary.csv";
    std::ofstream fout(path);
    if (!fout.is_open()) {
        fatal("[pi0_subtracted_kinematics] FATAL: cannot write " + path);
    }

    fout << "period,n_valid_bins,n_positive_weight_bins,n_zero_or_unusable_weight_bins,n_negative_signal_bins,"
         << "total_csv_raw_epgamma,total_csv_pi0_subtracted_signal,mean_positive_weight,rms_positive_weight,max_positive_weight\n";

    for (const PeriodDef& p : periods()) {
        auto it = csv.weight_summary_by_period.find(p.pretty);
        if (it == csv.weight_summary_by_period.end()) continue;
        const WeightSummary& s = it->second;
        double mean = std::numeric_limits<double>::quiet_NaN();
        double rms = std::numeric_limits<double>::quiet_NaN();
        if (s.n_positive_weight_bins > 0) {
            mean = s.sum_weight / double(s.n_positive_weight_bins);
            const double mean2 = s.sum_weight2 / double(s.n_positive_weight_bins);
            const double var = std::max(0.0, mean2 - mean * mean);
            rms = std::sqrt(var);
        }

        fout << p.pretty << ","
             << s.n_valid_bins << ","
             << s.n_positive_weight_bins << ","
             << s.n_zero_weight_bins << ","
             << s.n_negative_signal_bins << ","
             << std::setprecision(12) << s.total_csv_raw_epg << ","
             << std::setprecision(12) << s.total_csv_signal << ","
             << std::setprecision(12) << mean << ","
             << std::setprecision(12) << rms << ","
             << std::setprecision(12) << s.max_weight << "\n";
    }
}

static void write_fill_diagnostics_csv(const std::vector<FillSummary>& summaries,
                                       const std::string& diagnostics_dir) {
    std::filesystem::create_directories(diagnostics_dir);
    const std::string path = diagnostics_dir + "/event_fill_summary.csv";
    std::ofstream fout(path);
    if (!fout.is_open()) {
        fatal("[pi0_subtracted_kinematics] FATAL: cannot write " + path);
    }

    fout << "mode,category,panel,period,sample,entries,pass_global,pass_selector,pass_exclusivity,filled_events,filled_weight_sum\n";
    for (const FillSummary& s : summaries) {
        fout << s.mode << ","
             << s.category << ","
             << s.panel << ","
             << s.period << ","
             << s.sample << ","
             << s.entries << ","
             << s.pass_global << ","
             << s.pass_selector << ","
             << s.pass_exclusivity << ","
             << s.filled_events << ","
             << std::setprecision(12) << s.filled_weight_sum << "\n";
    }
}

static std::vector<PanelSpec> inclusive_period_panels(const std::vector<PeriodDef>& ps) {
    std::vector<PanelSpec> panels;
    for (const PeriodDef& p : ps) {
        PanelSpec panel;
        panel.label = p.pretty;
        panel.period = p;
        panels.push_back(panel);
    }
    return panels;
}

static std::vector<PanelSpec> topology_period_panels(const std::vector<PeriodDef>& ps,
                                                     int detector1,
                                                     int detector2) {
    std::vector<PanelSpec> panels;
    for (const PeriodDef& p : ps) {
        PanelSpec panel;
        panel.label = p.pretty;
        panel.period = p;
        panel.detector1_filter = detector1;
        panel.detector2_filter = detector2;
        panels.push_back(panel);
    }
    return panels;
}

static std::vector<PanelSpec> sp18_out_electron_sector_panels(const std::vector<PeriodDef>& ps) {
    auto it = std::find_if(ps.begin(), ps.end(), [](const PeriodDef& p) { return p.pretty == "Sp18 Out"; });
    if (it == ps.end()) {
        fatal("[pi0_subtracted_kinematics] FATAL: cannot find Sp18 Out period definition.");
    }

    std::vector<PanelSpec> panels;
    for (int sector = 1; sector <= 6; ++sector) {
        PanelSpec panel;
        panel.label = "Electron sector " + std::to_string(sector);
        panel.period = *it;
        panel.electron_sector_filter = sector;
        panels.push_back(panel);
    }
    return panels;
}

static void make_panel_plot_set(const std::string& category_key,
                                const std::string& category_title,
                                const std::string& out_dir,
                                DataShapeMode mode,
                                const std::vector<PanelSpec>& panels,
                                const std::vector<PlotVar>& vars,
                                const std::map<std::string, TTree*>& dataTrees,
                                const std::map<std::string, TTree*>& mcTrees,
                                const std::map<std::string, CutDict>& cuts,
                                const CsvInfo& csv,
                                std::vector<FillSummary>& fill_summaries) {
    std::filesystem::create_directories(out_dir);

    std::vector<std::vector<TH1D*> > data_hists(panels.size());
    std::vector<std::vector<TH1D*> > mc_hists(panels.size());

    for (size_t ip = 0; ip < panels.size(); ++ip) {
        for (const PlotVar& v : vars) {
            data_hists[ip].push_back(make_hist("h_data_" + sanitize(category_key) + "_" + data_mode_name(mode) + "_" + sanitize(panels[ip].label) + "_" + sanitize(v.branch), v));
            mc_hists[ip].push_back(make_hist("h_mc_" + sanitize(category_key) + "_" + data_mode_name(mode) + "_" + sanitize(panels[ip].label) + "_" + sanitize(v.branch), v));
        }
    }

    for (size_t ip = 0; ip < panels.size(); ++ip) {
        const PanelSpec& panel = panels[ip];
        auto itd = dataTrees.find(panel.period.data_key);
        auto itm = mcTrees.find(panel.period.mc_key);
        if (itd == dataTrees.end() || itd->second == nullptr || itm == mcTrees.end() || itm->second == nullptr) {
            fatal("[pi0_subtracted_kinematics] FATAL: missing tree while making category " + category_key + " panel " + panel.label);
        }

        fill_summaries.push_back(fill_tree_panel(itd->second,
                                                 panel,
                                                 vars,
                                                 cuts,
                                                 csv,
                                                 true,
                                                 mode,
                                                 category_key,
                                                 data_hists[ip]));
        fill_summaries.push_back(fill_tree_panel(itm->second,
                                                 panel,
                                                 vars,
                                                 cuts,
                                                 csv,
                                                 false,
                                                 mode,
                                                 category_key,
                                                 mc_hists[ip]));
    }

    for (size_t iv = 0; iv < vars.size(); ++iv) {
        std::vector<TH1D*> data_for_var;
        std::vector<TH1D*> mc_for_var;
        data_for_var.reserve(panels.size());
        mc_for_var.reserve(panels.size());
        for (size_t ip = 0; ip < panels.size(); ++ip) {
            data_for_var.push_back(data_hists[ip][iv]);
            mc_for_var.push_back(mc_hists[ip][iv]);
        }
        save_canvas_panels(vars[iv], panels, data_for_var, mc_for_var, mode, category_title, out_dir);
    }

    for (size_t ip = 0; ip < panels.size(); ++ip) {
        for (TH1D* h : data_hists[ip]) delete h;
        for (TH1D* h : mc_hists[ip]) delete h;
    }
}

} // namespace

bool plot_pi0_subtracted_dvcs_kinematics(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& out_dir,
    int max_workers) {

    (void)max_workers;

    try {
        ROOT::EnableThreadSafety();
        TH1::AddDirectory(kFALSE);
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);

        std::filesystem::create_directories(out_dir);
        std::filesystem::create_directories(out_dir + "/diagnostics");

        const CsvInfo csv = load_csv_info(csv_path);
        const std::map<std::string, CutDict> cuts = load_combined_cuts_json(combined_cuts_json);
        const std::vector<PeriodDef> ps = periods();
        const std::vector<PlotVar> vars = plot_vars();

        require_trees(ps, dvcsDataTrees, dvcsRecMcTrees);

        std::vector<FillSummary> fill_summaries;

        write_weight_diagnostics_csv(csv, out_dir + "/diagnostics");

        const std::vector<DataShapeMode> modes = {
            DataShapeMode::Pi0Subtracted,
            DataShapeMode::RawPostExclusivity
        };

        for (DataShapeMode mode : modes) {
            const std::string mode_dir = data_mode_name(mode);

            make_panel_plot_set("inclusive",
                                "Inclusive over topology and electron sector",
                                out_dir + "/inclusive/" + mode_dir,
                                mode,
                                inclusive_period_panels(ps),
                                vars,
                                dvcsDataTrees,
                                dvcsRecMcTrees,
                                cuts,
                                csv,
                                fill_summaries);

            struct TopologyRequest {
                int detector1;
                int detector2;
            };

            const std::vector<TopologyRequest> topologies = {
                {1, 1},
                {2, 1},
                {2, 0}
            };

            for (const TopologyRequest& tr : topologies) {
                const std::string tdir = topo_dir_name(tr.detector1, tr.detector2);
                make_panel_plot_set("topology_" + tdir,
                                    "Topology " + tdir,
                                    out_dir + "/topology/" + mode_dir + "/" + tdir,
                                    mode,
                                    topology_period_panels(ps, tr.detector1, tr.detector2),
                                    vars,
                                    dvcsDataTrees,
                                    dvcsRecMcTrees,
                                    cuts,
                                    csv,
                                    fill_summaries);
            }

            make_panel_plot_set("sp18_out_electron_sector",
                                "Sp18 Out by electron FD sector",
                                out_dir + "/electron_sector_sp18_out/" + mode_dir,
                                mode,
                                sp18_out_electron_sector_panels(ps),
                                vars,
                                dvcsDataTrees,
                                dvcsRecMcTrees,
                                cuts,
                                csv,
                                fill_summaries);
        }

        write_fill_diagnostics_csv(fill_summaries, out_dir + "/diagnostics");

        std::cout << "[pi0_subtracted_kinematics] Wrote inclusive, raw, topology-split, "
                  << "Sp18 Out electron-sector, and diagnostics outputs below "
                  << out_dir << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}
