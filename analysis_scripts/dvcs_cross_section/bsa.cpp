// bsa.cpp
// -----------------------------------------------------------------------------
// Pi0-subtracted, direct count-based beam-spin asymmetry stage for the DVCS
// pass-2 workflow.
//
// This stage intentionally stays at the measured-count level. It does not use
// acceptance unfolding, bin migration unfolding, radiative corrections, or model
// bin-centering factors. It does apply:
//   * the current process-wide global cuts from global_cuts.cpp,
//   * the data-derived channel/topology/period exclusivity windows in
//     output/jsons/combined_cuts.json,
//   * helicity-separated measured ep->eppi0 subtraction using the existing
//     bin-by-bin contamination ratios in the pass-2 CSV.
//
// Per period and CSV row:
//   G+/- = measured ep->epgamma counts after cuts
//   P+/- = measured ep->eppi0 counts after cuts
//   f_pi0 = contamination_ratio * N_norm(epgamma) / N_norm(eppi0)
//   S+/- = G+/- - f_pi0 * P+/-
//   A_LU = (S+ - S-) / [P_beam * (S+ + S-)]
//
// For combined groups with multiple beam polarizations:
//   A_LU = sum_i(S_i+ - S_i-) / sum_i[P_i * (S_i+ + S_i-)]
// -----------------------------------------------------------------------------

#include "bsa.h"
#include "global_cuts.h"

// ROOT
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

// JSON
#include <nlohmann/json.hpp>

// C++ stdlib
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
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

namespace {

static constexpr double PI = 3.14159265358979323846;
static constexpr double RAD2DEG = 180.0 / PI;

struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTextFont(42);
        gStyle->SetFrameLineWidth(2);
    }
} g_style_init;

static std::mutex g_root_bind_mutex;

static inline void fatal(const std::string& msg) {
    throw std::runtime_error(msg);
}

static inline std::string to_lower_ascii(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

static inline std::string sanitize_token(std::string s) {
    for (char& c : s) {
        if (!std::isalnum(static_cast<unsigned char>(c))) {
            c = '_';
        }
    }
    return s;
}

static inline double wrap_phi_deg(double phi_deg) {
    double p = std::fmod(phi_deg, 360.0);
    if (p < 0.0) {
        p += 360.0;
    }
    if (p >= 360.0) {
        p = std::nextafter(360.0, 0.0);
    }
    return p;
}

static inline bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (pmax_deg > pmin_deg) {
        return (phi_deg >= pmin_deg && phi_deg < pmax_deg);
    }
    return (phi_deg >= pmin_deg || phi_deg < pmax_deg);
}

static inline double phi_center_deg(double pmin_deg, double pmax_deg) {
    double lo = wrap_phi_deg(pmin_deg);
    double hi = wrap_phi_deg(pmax_deg);
    if (hi <= lo) {
        hi += 360.0;
    }
    return wrap_phi_deg(0.5 * (lo + hi));
}

static inline double phi_half_width_deg(double pmin_deg, double pmax_deg) {
    double lo = wrap_phi_deg(pmin_deg);
    double hi = wrap_phi_deg(pmax_deg);
    if (hi <= lo) {
        hi += 360.0;
    }
    return 0.5 * (hi - lo);
}

static inline double delta_phi_rad_from_two_phi(double phi_a, double phi_b) {
    double d = std::fmod(phi_a - phi_b, 2.0 * PI);
    if (d <= -PI) {
        d += 2.0 * PI;
    }
    if (d > PI) {
        d -= 2.0 * PI;
    }
    return std::fabs(d);
}

static inline std::string topo_key_from_detectors(int detector1, int detector2) {
    if (detector1 == 1 && detector2 == 1) return "FD_FD";
    if (detector1 == 2 && detector2 == 1) return "CD_FD";
    if (detector1 == 2 && detector2 == 0) return "CD_FT";
    return "";
}

static inline std::string topo_csv_from_key(const std::string& topo_key) {
    if (topo_key == "FD_FD") return "(FD, FD)";
    if (topo_key == "CD_FD") return "(CD, FD)";
    if (topo_key == "CD_FT") return "(CD, FT)";
    return "";
}

static const std::vector<std::string>& topo_keys() {
    static const std::vector<std::string> v = {"FD_FD", "CD_FD", "CD_FT"};
    return v;
}

static inline std::string period_display_from_tree_key(const std::string& tree_key) {
    const std::string s = to_lower_ascii(tree_key);
    const auto has = [&](const char* sub) { return s.find(sub) != std::string::npos; };

    if (has("fa18") && has("inb") && !has("supp")) return "Fa18 Inb";
    if (has("fa18") && has("out")) return "Fa18 Out";
    if (has("sp19") && has("inb")) return "Sp19 Inb";
    if (has("sp18") && has("inb")) return "Sp18 Inb";
    if (has("sp18") && has("out")) return "Sp18 Out";
    return "";
}

static inline std::string period_code_from_display(const std::string& period) {
    if (period == "Fa18 Inb") return "Fa18_Inb";
    if (period == "Fa18 Out") return "Fa18_Out";
    if (period == "Sp19 Inb") return "Sp19_Inb";
    if (period == "Sp18 Inb") return "Sp18_Inb";
    if (period == "Sp18 Out") return "Sp18_Out";
    fatal("[bsa] unknown display period label: " + period);
    return "";
}

static inline std::string period_global_cuts_label_from_display(const std::string& period) {
    if (period == "Fa18 Inb") return "fa18_inb";
    if (period == "Fa18 Out") return "fa18_out";
    if (period == "Sp19 Inb") return "sp19_inb";
    if (period == "Sp18 Inb") return "sp18_inb";
    if (period == "Sp18 Out") return "sp18_out";
    fatal("[bsa] unknown display period label for global cuts: " + period);
    return "";
}

static inline double beam_pol_for_period(const std::string& period, const BSAOptions& opt) {
    if (period == "Sp18 Inb") return opt.beam_pol_sp18_inb;
    if (period == "Sp18 Out") return opt.beam_pol_sp18_out;
    if (period == "Fa18 Inb") return opt.beam_pol_fa18_inb;
    if (period == "Fa18 Out") return opt.beam_pol_fa18_out;
    if (period == "Sp19 Inb") return opt.beam_pol_sp19_inb;
    fatal("[bsa] no beam polarization configured for period: " + period);
    return std::numeric_limits<double>::quiet_NaN();
}

static const std::vector<std::string>& base_period_order() {
    static const std::vector<std::string> v = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    return v;
}

static const std::vector<std::string>& output_group_order() {
    static const std::vector<std::string> v = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out",
        "Fa18", "Sp18", "10.6 GeV"
    };
    return v;
}

static std::vector<std::string> component_periods_for_group(const std::string& group) {
    if (group == "Fa18 Inb") return {"Fa18 Inb"};
    if (group == "Fa18 Out") return {"Fa18 Out"};
    if (group == "Sp19 Inb") return {"Sp19 Inb"};
    if (group == "Sp18 Inb") return {"Sp18 Inb"};
    if (group == "Sp18 Out") return {"Sp18 Out"};
    if (group == "Fa18") return {"Fa18 Inb", "Fa18 Out"};
    if (group == "Sp18") return {"Sp18 Inb", "Sp18 Out"};
    if (group == "10.6 GeV") return {"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"};
    fatal("[bsa] unknown output group: " + group);
    return {};
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

static void load_csv(const std::string& path, CSV& csv) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        fatal("[bsa] cannot open CSV: " + path);
    }

    std::string line;
    if (!std::getline(fin, line)) {
        fatal("[bsa] empty CSV: " + path);
    }

    csv.header = split_csv_line(line);
    csv.index.clear();
    for (int i = 0; i < static_cast<int>(csv.header.size()); ++i) {
        if (csv.index.count(csv.header[i])) {
            fatal("[bsa] duplicate CSV column: " + csv.header[i]);
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
            fatal("[bsa] CSV row width mismatch in " + path);
        }
        csv.rows.push_back(std::move(row));
    }
}

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";
    std::ofstream fout(tmp);
    if (!fout.is_open()) {
        fatal("[bsa] cannot write temporary CSV: " + tmp);
    }

    auto write_cell = [&](const std::string& s) {
        const bool quote =
            s.find(',') != std::string::npos ||
            s.find('"') != std::string::npos ||
            s.find('\n') != std::string::npos ||
            s.find('\r') != std::string::npos;

        if (!quote) {
            fout << s;
            return;
        }

        fout << '"';
        for (char c : s) {
            if (c == '"') {
                fout << "\"\"";
            } else {
                fout << c;
            }
        }
        fout << '"';
    };

    for (size_t i = 0; i < csv.header.size(); ++i) {
        write_cell(csv.header[i]);
        if (i + 1 < csv.header.size()) fout << ',';
    }
    fout << '\n';

    for (const auto& row : csv.rows) {
        for (size_t i = 0; i < row.size(); ++i) {
            write_cell(row[i]);
            if (i + 1 < row.size()) fout << ',';
        }
        fout << '\n';
    }

    fout.close();
    if (!fout) {
        fatal("[bsa] failed while writing temporary CSV: " + tmp);
    }

    (void)std::remove(path.c_str());
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        fatal("[bsa] failed to rename " + tmp + " to " + path);
    }
}

static int col_strict(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);
    if (it == csv.index.end()) {
        fatal("[bsa] missing required CSV column: " + name);
    }
    return it->second;
}

static int col_optional(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);
    if (it == csv.index.end()) {
        return -1;
    }
    return it->second;
}

static double to_double_strict(const std::string& s, const std::string& what) {
    if (s.empty()) {
        fatal("[bsa] empty numeric cell for " + what);
    }
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str()) {
        fatal("[bsa] parse failure for " + what + ": " + s);
    }
    return v;
}

static bool to_bool_valid(const std::string& s) {
    return (s == "1" || s == "1.0" || s == "true" || s == "TRUE");
}

static bool parse_first_number(const std::string& s, double& value) {
    value = std::numeric_limits<double>::quiet_NaN();
    if (s.empty()) {
        return false;
    }

    size_t i = 0;
    while (i < s.size() && !(std::isdigit(static_cast<unsigned char>(s[i])) ||
                             s[i] == '-' || s[i] == '+' || s[i] == '.')) {
        ++i;
    }
    if (i >= s.size()) {
        return false;
    }

    char* end = nullptr;
    const double v = std::strtod(s.c_str() + i, &end);
    if (end == s.c_str() + i || !std::isfinite(v)) {
        return false;
    }

    value = v;
    return true;
}

static double cell_value_or_zero(const CSV& csv, int row_index, const std::string& col_name) {
    const int c = col_optional(csv, col_name);
    if (c < 0) {
        return 0.0;
    }
    double v = 0.0;
    if (!parse_first_number(csv.rows[row_index][c], v)) {
        return 0.0;
    }
    if (!std::isfinite(v)) {
        return 0.0;
    }
    return v;
}

// -----------------------------------------------------------------------------
// Row binning
// -----------------------------------------------------------------------------

struct RowBin {
    double xBmin = std::numeric_limits<double>::quiet_NaN();
    double xBmax = std::numeric_limits<double>::quiet_NaN();
    double Q2min = std::numeric_limits<double>::quiet_NaN();
    double Q2max = std::numeric_limits<double>::quiet_NaN();
    double tmin = std::numeric_limits<double>::quiet_NaN();
    double tmax = std::numeric_limits<double>::quiet_NaN();
    double pmin = std::numeric_limits<double>::quiet_NaN();
    double pmax = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

struct AxisBin {
    double min = 0.0;
    double max = 0.0;
};

struct FastBinning {
    std::vector<AxisBin> xbins;
    std::vector<AxisBin> qbins;
    std::vector<AxisBin> tbins;
    std::vector<std::vector<std::vector<std::vector<int>>>> rows_by_xqt;
};

static std::vector<RowBin> load_row_bins_from_csv(const CSV& csv) {
    const int c_xBmin = col_strict(csv, "xBmin");
    const int c_xBmax = col_strict(csv, "xBmax");
    const int c_Q2min = col_strict(csv, "Q2min");
    const int c_Q2max = col_strict(csv, "Q2max");
    const int c_tmin = col_strict(csv, "t_abs_min");
    const int c_tmax = col_strict(csv, "t_abs_max");
    const int c_pmin = col_strict(csv, "phimin");
    const int c_pmax = col_strict(csv, "phimax");
    const int c_valid = col_strict(csv, "valid bin");

    std::vector<RowBin> rows;
    rows.reserve(csv.rows.size());
    for (int r = 0; r < static_cast<int>(csv.rows.size()); ++r) {
        const auto& row = csv.rows[r];
        RowBin b;
        b.xBmin = to_double_strict(row[c_xBmin], "xBmin");
        b.xBmax = to_double_strict(row[c_xBmax], "xBmax");
        b.Q2min = to_double_strict(row[c_Q2min], "Q2min");
        b.Q2max = to_double_strict(row[c_Q2max], "Q2max");
        b.tmin = to_double_strict(row[c_tmin], "t_abs_min");
        b.tmax = to_double_strict(row[c_tmax], "t_abs_max");
        b.pmin = to_double_strict(row[c_pmin], "phimin");
        b.pmax = to_double_strict(row[c_pmax], "phimax");
        b.valid = to_bool_valid(row[c_valid]);
        rows.push_back(b);
    }
    return rows;
}

static void add_unique_axis_bin(std::vector<AxisBin>& bins, double minv, double maxv) {
    for (const AxisBin& b : bins) {
        if (b.min == minv && b.max == maxv) {
            return;
        }
    }
    bins.push_back({minv, maxv});
}

static void sort_axis_bins(std::vector<AxisBin>& bins) {
    std::sort(bins.begin(), bins.end(), [](const AxisBin& a, const AxisBin& b) {
        if (a.min != b.min) return a.min < b.min;
        return a.max < b.max;
    });
}

static int find_axis_bin_index(const std::vector<AxisBin>& bins, double value) {
    for (int i = 0; i < static_cast<int>(bins.size()); ++i) {
        if (value >= bins[i].min && value < bins[i].max) {
            return i;
        }
    }
    return -1;
}

static int find_axis_bin_exact(const std::vector<AxisBin>& bins, double minv, double maxv) {
    for (int i = 0; i < static_cast<int>(bins.size()); ++i) {
        if (bins[i].min == minv && bins[i].max == maxv) {
            return i;
        }
    }
    return -1;
}

static FastBinning build_fast_binning(const std::vector<RowBin>& rows) {
    FastBinning fb;

    for (const RowBin& r : rows) {
        if (!r.valid) {
            continue;
        }
        add_unique_axis_bin(fb.xbins, r.xBmin, r.xBmax);
        add_unique_axis_bin(fb.qbins, r.Q2min, r.Q2max);
        add_unique_axis_bin(fb.tbins, r.tmin, r.tmax);
    }

    sort_axis_bins(fb.xbins);
    sort_axis_bins(fb.qbins);
    sort_axis_bins(fb.tbins);

    fb.rows_by_xqt.resize(fb.xbins.size());
    for (size_t ix = 0; ix < fb.xbins.size(); ++ix) {
        fb.rows_by_xqt[ix].resize(fb.qbins.size());
        for (size_t iq = 0; iq < fb.qbins.size(); ++iq) {
            fb.rows_by_xqt[ix][iq].resize(fb.tbins.size());
        }
    }

    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        const RowBin& row = rows[r];
        if (!row.valid) {
            continue;
        }

        const int ix = find_axis_bin_exact(fb.xbins, row.xBmin, row.xBmax);
        const int iq = find_axis_bin_exact(fb.qbins, row.Q2min, row.Q2max);
        const int it = find_axis_bin_exact(fb.tbins, row.tmin, row.tmax);
        if (ix < 0 || iq < 0 || it < 0) {
            fatal("[bsa] failed to build fast row-bin lookup");
        }
        fb.rows_by_xqt[ix][iq][it].push_back(r);
    }

    std::cout << "[bsa] Fast bin lookup built with "
              << fb.xbins.size() << " xB bins, "
              << fb.qbins.size() << " Q2 bins, "
              << fb.tbins.size() << " |t| bins.\n";
    return fb;
}

// -----------------------------------------------------------------------------
// Exclusivity cuts
// -----------------------------------------------------------------------------

struct SigmaStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std = std::numeric_limits<double>::quiet_NaN();
    double cut_low = std::numeric_limits<double>::quiet_NaN();
    double cut_high = std::numeric_limits<double>::quiet_NaN();
    double quantile = 0.0;
    std::string mode = "symmetric_3sigma";
};

using CutVarMap = std::unordered_map<std::string, SigmaStats>;
using TopoCutMap = std::unordered_map<std::string, CutVarMap>;

static bool within_cut_window(double v, const SigmaStats& s) {
    if (!std::isfinite(v)) {
        return false;
    }

    if (s.mode == "upper_quantile") {
        if (!std::isfinite(s.cut_high)) {
            return true;
        }
        return v <= s.cut_high;
    }

    double lo = s.cut_low;
    double hi = s.cut_high;
    if (!(std::isfinite(lo) && std::isfinite(hi)) || hi <= lo) {
        if (!std::isfinite(s.mean) || !std::isfinite(s.std) || s.std <= 0.0) {
            return true;
        }
        lo = s.mean - 3.0 * s.std;
        hi = s.mean + 3.0 * s.std;
    }

    return (v >= lo && v <= hi);
}

static TopoCutMap load_combined_cuts(const std::string& combined_cuts_json) {
    std::ifstream fin(combined_cuts_json);
    if (!fin.is_open()) {
        fatal("[bsa] cannot open combined cuts JSON: " + combined_cuts_json);
    }

    nlohmann::json j;
    fin >> j;
    if (!j.is_object()) {
        fatal("[bsa] combined cuts JSON is not an object: " + combined_cuts_json);
    }

    TopoCutMap out;
    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const auto& block = it.value();
        if (!block.is_object() || !block.contains("data") || !block["data"].is_object()) {
            continue;
        }

        CutVarMap vm;
        for (auto vit = block["data"].begin(); vit != block["data"].end(); ++vit) {
            const std::string var = vit.key();
            const auto& stats = vit.value();
            if (!stats.is_object()) {
                continue;
            }

            SigmaStats s;
            try {
                if (stats.contains("mean")) s.mean = stats["mean"].get<double>();
                if (stats.contains("std")) s.std = stats["std"].get<double>();
                if (stats.contains("cut_low")) s.cut_low = stats["cut_low"].get<double>();
                if (stats.contains("cut_high")) s.cut_high = stats["cut_high"].get<double>();
                if (stats.contains("quantile")) s.quantile = stats["quantile"].get<double>();
                if (stats.contains("mode")) s.mode = stats["mode"].get<std::string>();
            } catch (...) {
                continue;
            }

            if (!std::isfinite(s.cut_low) || !std::isfinite(s.cut_high) || s.cut_high <= s.cut_low) {
                if (std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0) {
                    s.cut_low = s.mean - 3.0 * s.std;
                    s.cut_high = s.mean + 3.0 * s.std;
                }
            }

            vm.emplace(var, s);
        }

        if (!vm.empty()) {
            out.emplace(key, std::move(vm));
        }
    }

    std::cout << "[bsa] Loaded " << out.size()
              << " data cut blocks from " << combined_cuts_json << "\n";
    return out;
}

// -----------------------------------------------------------------------------
// Branch binding and event selection
// -----------------------------------------------------------------------------

struct BranchBinder {
    int runnum = 0; bool has_runnum = false;
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;
    int helicity = 0; bool has_helicity = false;

    double x = 0.0; bool has_x = false;
    double Q2 = 0.0; bool has_Q2 = false;
    double t1 = 0.0; bool has_t1 = false;
    double phi2 = 0.0; bool has_phi2 = false;

    double open_angle_ep2 = 0.0; bool has_open_angle_ep2 = false;
    double pTmiss = 0.0; bool has_pTmiss = false;
    double Emiss2 = 0.0; bool has_Emiss2 = false;
    double Mx2 = 0.0; bool has_Mx2 = false;
    double Mx2_1 = 0.0; bool has_Mx2_1 = false;
    double Mx2_2 = 0.0; bool has_Mx2_2 = false;
    double xF = 0.0; bool has_xF = false;
    double Delta_phi = 0.0; bool has_Delta_phi = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double theta_pi0_pi0 = 0.0; bool has_theta_pi0_pi0 = false;

    double e_p = 0.0; bool has_e_p = false;
    double e_theta = 0.0; bool has_e_theta = false;
    double e_phi = 0.0; bool has_e_phi = false;
    double p1_theta = 0.0; bool has_p1_theta = false;
    double p1_phi = 0.0; bool has_p1_phi = false;
    double p2_p = 0.0; bool has_p2_p = false;
    double p2_theta = 0.0; bool has_p2_theta = false;
    double p2_phi = 0.0; bool has_p2_phi = false;

    void bind(TTree* t) {
        if (!t) return;

        std::lock_guard<std::mutex> lock(g_root_bind_mutex);
        t->SetBranchStatus("*", 0);

        auto enable = [&](const char* name) {
            if (t->GetBranch(name)) {
                t->SetBranchStatus(name, 1);
            }
        };

        const char* names[] = {
            "runnum", "detector1", "detector2", "helicity",
            "x", "Q2", "t1", "phi2",
            "open_angle_ep2", "pTmiss", "Emiss2", "Mx2", "Mx2_1", "Mx2_2",
            "xF", "Delta_phi", "theta_gamma_gamma", "theta_pi0_pi0",
            "e_p", "e_theta", "e_phi", "p1_theta", "p1_phi",
            "p2_p", "p2_theta", "p2_phi"
        };
        for (const char* name : names) {
            enable(name);
        }

        t->SetCacheSize(0);

        auto bind_int = [&](const char* name, int* ptr, bool& has) {
            if (t->GetBranch(name)) {
                t->SetBranchAddress(name, ptr);
                has = true;
            }
        };
        auto bind_double = [&](const char* name, double* ptr, bool& has) {
            if (t->GetBranch(name)) {
                t->SetBranchAddress(name, ptr);
                has = true;
            }
        };

        bind_int("runnum", &runnum, has_runnum);
        bind_int("detector1", &detector1, has_detector1);
        bind_int("detector2", &detector2, has_detector2);
        bind_int("helicity", &helicity, has_helicity);

        bind_double("x", &x, has_x);
        bind_double("Q2", &Q2, has_Q2);
        bind_double("t1", &t1, has_t1);
        bind_double("phi2", &phi2, has_phi2);
        bind_double("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bind_double("pTmiss", &pTmiss, has_pTmiss);
        bind_double("Emiss2", &Emiss2, has_Emiss2);
        bind_double("Mx2", &Mx2, has_Mx2);
        bind_double("Mx2_1", &Mx2_1, has_Mx2_1);
        bind_double("Mx2_2", &Mx2_2, has_Mx2_2);
        bind_double("xF", &xF, has_xF);
        bind_double("Delta_phi", &Delta_phi, has_Delta_phi);
        bind_double("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        bind_double("theta_pi0_pi0", &theta_pi0_pi0, has_theta_pi0_pi0);
        bind_double("e_p", &e_p, has_e_p);
        bind_double("e_theta", &e_theta, has_e_theta);
        bind_double("e_phi", &e_phi, has_e_phi);
        bind_double("p1_theta", &p1_theta, has_p1_theta);
        bind_double("p1_phi", &p1_phi, has_p1_phi);
        bind_double("p2_p", &p2_p, has_p2_p);
        bind_double("p2_theta", &p2_theta, has_p2_theta);
        bind_double("p2_phi", &p2_phi, has_p2_phi);
    }

    double phi_deg() const { return wrap_phi_deg(phi2 * RAD2DEG); }
    double t_abs() const { return std::fabs(t1); }

    double delta_phi_value(bool& has_val) const {
        if (has_Delta_phi) {
            has_val = true;
            return Delta_phi;
        }
        if (has_p1_phi && has_p2_phi) {
            has_val = true;
            return delta_phi_rad_from_two_phi(p1_phi, p2_phi);
        }
        has_val = false;
        return 0.0;
    }
};

static bool value_for_sigma_var(const BranchBinder& b,
                                const std::string& var,
                                double& value) {
    if (var == "Emiss2" && b.has_Emiss2) { value = b.Emiss2; return true; }
    if (var == "Mx2" && b.has_Mx2) { value = b.Mx2; return true; }
    if (var == "Mx2_1" && b.has_Mx2_1) { value = b.Mx2_1; return true; }
    if (var == "Mx2_2" && b.has_Mx2_2) { value = b.Mx2_2; return true; }
    if (var == "pTmiss" && b.has_pTmiss) { value = b.pTmiss; return true; }
    if (var == "xF" && b.has_xF) { value = b.xF; return true; }
    if (var == "theta_gamma_gamma" && b.has_theta_gamma_gamma) { value = b.theta_gamma_gamma; return true; }
    if (var == "theta_pi0_pi0" && b.has_theta_pi0_pi0) { value = b.theta_pi0_pi0; return true; }
    if (var == "Delta_phi") {
        bool has_val = false;
        value = b.delta_phi_value(has_val);
        return has_val;
    }
    return false;
}

static bool passes_sigma_cuts(const TopoCutMap& cuts,
                              const std::string& cut_key,
                              const BranchBinder& b) {
    auto it = cuts.find(cut_key);
    if (it == cuts.end()) {
        fatal("[bsa] missing exclusivity cut key in combined_cuts.json: " + cut_key);
    }

    for (const auto& kv : it->second) {
        double v = 0.0;
        if (!value_for_sigma_var(b, kv.first, v)) {
            continue;
        }
        if (!within_cut_window(v, kv.second)) {
            return false;
        }
    }
    return true;
}

static bool passes_global_cuts_dispatch(const BranchBinder& b,
                                        const std::string& period_display) {
    const GlobalCutConfig& cfg = default_global_cuts();
    const std::string period_label = period_global_cuts_label_from_display(period_display);

    if (!(b.has_t1 && b.has_open_angle_ep2)) return false;
    if (cfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (b.has_runnum && is_excluded_run(b.runnum)) return false;

    if (cfg.enable_topology_filter || global_cuts_require_sector_phi(cfg) ||
        cfg.enable_dvcsgen_ycol_cut || global_cuts_require_auxiliary_kinematics(cfg)) {
        if (!(b.has_detector1 && b.has_detector2)) {
            fatal("[bsa] global topology/sector/auxiliary cuts require detector1 and detector2 branches");
        }
    }

    if (global_cuts_require_sector_phi(cfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            fatal("[bsa] sector selection requires e_phi, p1_phi and p2_phi branches");
        }
    }

    if (global_cuts_require_auxiliary_kinematics(cfg)) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p1_theta && b.has_p1_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[bsa] auxiliary fiducial cuts require e/p/gamma kinematic branches");
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p1_theta, b.p1_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[bsa] dvcsgen ycol cut requires e/p2 kinematic branches");
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
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

// -----------------------------------------------------------------------------
// Counts and BSA math
// -----------------------------------------------------------------------------

struct HelCounts {
    double plus = 0.0;
    double minus = 0.0;
};

using RowCounts = std::unordered_map<int, HelCounts>;
using PeriodCounts = std::unordered_map<std::string, RowCounts>;

static void add_event(HelCounts& h, int helicity) {
    if (helicity > 0) {
        h.plus += 1.0;
    } else if (helicity < 0) {
        h.minus += 1.0;
    }
}

static PeriodCounts accumulate_counts(const std::map<std::string, TTree*>& trees,
                                      const std::string& channel_prefix,
                                      const std::vector<RowBin>& rows,
                                      const FastBinning& fast_bins,
                                      const TopoCutMap& sigma_cuts) {
    PeriodCounts out;
    for (const std::string& period : base_period_order()) {
        out[period] = RowCounts();
    }

    for (const auto& kv : trees) {
        const std::string& tree_key = kv.first;
        TTree* tree = kv.second;
        if (!tree) {
            continue;
        }

        const std::string period = period_display_from_tree_key(tree_key);
        if (period.empty()) {
            std::cout << "[bsa] Skipping non-canonical data tree key for "
                      << channel_prefix << ": " << tree_key << "\n";
            continue;
        }

        BranchBinder b;
        b.bind(tree);

        if (!(b.has_detector1 && b.has_detector2 && b.has_helicity &&
              b.has_x && b.has_Q2 && b.has_t1 && b.has_phi2)) {
            fatal("[bsa] tree " + tree_key +
                  " is missing one or more required branches: detector1, detector2, helicity, x, Q2, t1, phi2");
        }

        const std::string period_code = period_code_from_display(period);
        const Long64_t n_entries = tree->GetEntries();

        long long n_topology = 0;
        long long n_global = 0;
        long long n_sigma = 0;
        long long n_matched = 0;

        for (Long64_t i = 0; i < n_entries; ++i) {
            tree->GetEntry(i);

            const std::string topo = topo_key_from_detectors(b.detector1, b.detector2);
            if (topo.empty()) {
                continue;
            }
            ++n_topology;

            if (!passes_global_cuts_dispatch(b, period)) {
                continue;
            }
            ++n_global;

            const std::string cut_key = channel_prefix + "_" + period_code + "_" + topo;
            if (!passes_sigma_cuts(sigma_cuts, cut_key, b)) {
                continue;
            }
            ++n_sigma;

            const int ix = find_axis_bin_index(fast_bins.xbins, b.x);
            if (ix < 0) continue;
            const int iq = find_axis_bin_index(fast_bins.qbins, b.Q2);
            if (iq < 0) continue;
            const int it = find_axis_bin_index(fast_bins.tbins, b.t_abs());
            if (it < 0) continue;

            const double phi_deg = b.phi_deg();
            const std::vector<int>& candidate_rows = fast_bins.rows_by_xqt[ix][iq][it];
            bool matched = false;

            for (int row_index : candidate_rows) {
                const RowBin& rb = rows[row_index];
                if (!row_accepts_phi(phi_deg, rb.pmin, rb.pmax)) {
                    continue;
                }
                add_event(out[period][row_index], b.helicity);
                matched = true;
            }

            if (matched) {
                ++n_matched;
            }
        }

        std::cout << "[bsa] channel=" << channel_prefix
                  << " tree=" << tree_key
                  << " period=" << period
                  << " entries=" << static_cast<long long>(n_entries)
                  << " topology=" << n_topology
                  << " global=" << n_global
                  << " sigma=" << n_sigma
                  << " matched=" << n_matched
                  << "\n";
    }

    return out;
}

struct LeakInfo {
    bool valid = false;
    double contamination_ratio = 0.0;
    double norm_epg = 0.0;
    double norm_eppi0 = 0.0;
    double f_pi0 = 0.0;
};

using PeriodLeakRows = std::unordered_map<std::string, std::vector<LeakInfo>>;

static std::string normalized_raw_col(const std::string& channel,
                                      const std::string& topo_csv,
                                      const std::string& period,
                                      const std::string& helicity) {
    return "normalized raw yield, " + channel + ", " + topo_csv + ", exp, " + period + ", " + helicity;
}

static PeriodLeakRows build_pi0_leakage(const CSV& csv, const BSAOptions& options) {
    PeriodLeakRows leaks;

    for (const std::string& period : base_period_order()) {
        std::vector<LeakInfo> v(csv.rows.size());
        const std::string c_cont = "contamination ratio, " + period;
        const int c_cont_idx = col_optional(csv, c_cont);
        if (c_cont_idx < 0) {
            std::cout << "[bsa] WARNING: missing " << c_cont
                      << "; pi0 subtraction scale will be zero for this period.\n";
        }

        int valid_count = 0;
        for (int r = 0; r < static_cast<int>(csv.rows.size()); ++r) {
            LeakInfo info;

            double c = 0.0;
            if (c_cont_idx >= 0) {
                (void)parse_first_number(csv.rows[r][c_cont_idx], c);
            }
            info.contamination_ratio = c;

            double norm_epg = 0.0;
            double norm_eppi0 = 0.0;
            for (const std::string& topo_key : topo_keys()) {
                const std::string topo_csv = topo_csv_from_key(topo_key);
                norm_epg += cell_value_or_zero(csv, r, normalized_raw_col("ep->epg", topo_csv, period, "unpol"));
                norm_eppi0 += cell_value_or_zero(csv, r, normalized_raw_col("ep->eppi0", topo_csv, period, "unpol"));
            }

            info.norm_epg = norm_epg;
            info.norm_eppi0 = norm_eppi0;

            if (options.enable_pi0_subtraction && c > 0.0 && norm_epg > 0.0 && norm_eppi0 > 0.0) {
                info.f_pi0 = c * norm_epg / norm_eppi0;
                info.valid = std::isfinite(info.f_pi0) && info.f_pi0 >= 0.0;
            } else {
                info.f_pi0 = 0.0;
                info.valid = true;
            }

            if (info.valid) {
                ++valid_count;
            }
            v[r] = info;
        }

        std::cout << "[bsa] pi0 leakage scale factors for " << period
                  << ": valid rows=" << valid_count << "/" << csv.rows.size() << "\n";
        leaks[period] = std::move(v);
    }

    return leaks;
}

struct AsymResult {
    bool valid = false;
    double value = 0.0;
    double stat = 0.0;

    double raw_value = std::numeric_limits<double>::quiet_NaN();
    double raw_stat = std::numeric_limits<double>::quiet_NaN();
    bool raw_valid = false;

    double pi0_value = std::numeric_limits<double>::quiet_NaN();
    double pi0_stat = std::numeric_limits<double>::quiet_NaN();
    bool pi0_valid = false;

    double g_plus = 0.0;
    double g_minus = 0.0;
    double p_plus = 0.0;
    double p_minus = 0.0;
    double s_plus = 0.0;
    double s_minus = 0.0;

    double f_pi0_effective = 0.0;
    double contamination_ratio_effective = 0.0;
    double polarized_denominator = 0.0;
};

static AsymResult compute_count_bsa_simple(double plus, double minus, double pol) {
    AsymResult r;
    const double s = plus + minus;
    const double d = plus - minus;
    const double den = pol * s;
    if (!(s > 0.0 && den > 0.0)) {
        return r;
    }
    const double var_num = std::max(0.0, s - (d * d / s));
    r.value = d / den;
    r.stat = std::sqrt(var_num) / den;
    r.valid = std::isfinite(r.value) && std::isfinite(r.stat);
    return r;
}

static AsymResult compute_group_bsa(const PeriodCounts& gamma_counts,
                                    const PeriodCounts& pi0_counts,
                                    const PeriodLeakRows& leaks,
                                    const std::vector<std::string>& component_periods,
                                    int row_index,
                                    const BSAOptions& opt) {
    AsymResult r;

    double numerator = 0.0;
    double denominator = 0.0;
    double variance = 0.0;

    double raw_num = 0.0;
    double raw_den = 0.0;
    double raw_var_num = 0.0;

    double pi0_num = 0.0;
    double pi0_den = 0.0;
    double pi0_var_num = 0.0;

    double weighted_f_num = 0.0;
    double weighted_c_num = 0.0;
    double weight_sum = 0.0;

    struct ComponentForVariance {
        double pol = 0.0;
        double s_plus = 0.0;
        double s_minus = 0.0;
        double var_s_plus = 0.0;
        double var_s_minus = 0.0;
    };
    std::vector<ComponentForVariance> comps;

    for (const std::string& period : component_periods) {
        const double pol = beam_pol_for_period(period, opt);

        HelCounts g;
        auto igp = gamma_counts.find(period);
        if (igp != gamma_counts.end()) {
            auto igr = igp->second.find(row_index);
            if (igr != igp->second.end()) {
                g = igr->second;
            }
        }

        HelCounts p;
        auto ipp = pi0_counts.find(period);
        if (ipp != pi0_counts.end()) {
            auto ipr = ipp->second.find(row_index);
            if (ipr != ipp->second.end()) {
                p = ipr->second;
            }
        }

        double f = 0.0;
        double c = 0.0;
        auto ilp = leaks.find(period);
        if (ilp != leaks.end() && row_index >= 0 && row_index < static_cast<int>(ilp->second.size())) {
            f = ilp->second[row_index].f_pi0;
            c = ilp->second[row_index].contamination_ratio;
        }

        const double gp = g.plus;
        const double gm = g.minus;
        const double pp = p.plus;
        const double pm = p.minus;
        const double sp = gp - f * pp;
        const double sm = gm - f * pm;

        r.g_plus += gp;
        r.g_minus += gm;
        r.p_plus += pp;
        r.p_minus += pm;
        r.s_plus += sp;
        r.s_minus += sm;

        const double raw_s = gp + gm;
        const double raw_d = gp - gm;
        if (raw_s > 0.0) {
            raw_num += raw_d;
            raw_den += pol * raw_s;
            raw_var_num += std::max(0.0, raw_s - (raw_d * raw_d / raw_s));
        }

        const double pi0_s = pp + pm;
        const double pi0_d = pp - pm;
        if (pi0_s > 0.0) {
            pi0_num += pi0_d;
            pi0_den += pol * pi0_s;
            pi0_var_num += std::max(0.0, pi0_s - (pi0_d * pi0_d / pi0_s));
        }

        const double sig_s = sp + sm;
        const double sig_d = sp - sm;
        if (sig_s > 0.0) {
            numerator += sig_d;
            denominator += pol * sig_s;
            comps.push_back(ComponentForVariance{pol, sp, sm, gp + f * f * pp, gm + f * f * pm});
            weighted_f_num += f * sig_s;
            weighted_c_num += c * sig_s;
            weight_sum += sig_s;
        }
    }

    if (raw_den > 0.0) {
        r.raw_value = raw_num / raw_den;
        r.raw_stat = std::sqrt(std::max(0.0, raw_var_num)) / raw_den;
        r.raw_valid = std::isfinite(r.raw_value) && std::isfinite(r.raw_stat);
    }

    if (pi0_den > 0.0) {
        r.pi0_value = pi0_num / pi0_den;
        r.pi0_stat = std::sqrt(std::max(0.0, pi0_var_num)) / pi0_den;
        r.pi0_valid = std::isfinite(r.pi0_value) && std::isfinite(r.pi0_stat);
    }

    r.polarized_denominator = denominator;
    if (weight_sum > 0.0) {
        r.f_pi0_effective = weighted_f_num / weight_sum;
        r.contamination_ratio_effective = weighted_c_num / weight_sum;
    }

    if (!(denominator > 0.0)) {
        return r;
    }

    r.value = numerator / denominator;

    for (const ComponentForVariance& comp : comps) {
        const double dA_dSp = (denominator - numerator * comp.pol) / (denominator * denominator);
        const double dA_dSm = (-denominator - numerator * comp.pol) / (denominator * denominator);
        variance += dA_dSp * dA_dSp * comp.var_s_plus;
        variance += dA_dSm * dA_dSm * comp.var_s_minus;
    }

    r.stat = std::sqrt(std::max(0.0, variance));
    r.valid = std::isfinite(r.value) && std::isfinite(r.stat);
    return r;
}

static std::string fmt_tuple(double value, double stat) {
    if (!(std::isfinite(value) && std::isfinite(stat))) {
        return "";
    }
    std::ostringstream ss;
    ss << std::setprecision(12) << "(" << value << "," << stat << ",0)";
    return ss.str();
}

// -----------------------------------------------------------------------------
// JSON and plotting
// -----------------------------------------------------------------------------

static void write_json_summary(const std::string& path,
                               const std::vector<RowBin>& rows,
                               const std::map<std::string, std::vector<AsymResult>>& results) {
    nlohmann::json j;
    j["description"] = "Pi0-subtracted direct-count DVCS beam-spin asymmetries after global cuts and data exclusivity cuts.";
    j["estimator"] = "A_LU = sum_i(Splus_i - Sminus_i) / sum_i[Pbeam_i*(Splus_i + Sminus_i)], S± = G± - f_pi0*P±";

    for (const auto& group_pair : results) {
        const std::string& group = group_pair.first;
        const auto& vec = group_pair.second;
        nlohmann::json rows_json = nlohmann::json::array();

        for (int r = 0; r < static_cast<int>(vec.size()); ++r) {
            const AsymResult& a = vec[r];
            const RowBin& rb = rows[r];
            nlohmann::json row;
            row["row"] = r;
            row["xBmin"] = rb.xBmin;
            row["xBmax"] = rb.xBmax;
            row["Q2min"] = rb.Q2min;
            row["Q2max"] = rb.Q2max;
            row["t_abs_min"] = rb.tmin;
            row["t_abs_max"] = rb.tmax;
            row["phimin"] = rb.pmin;
            row["phimax"] = rb.pmax;
            row["valid"] = a.valid;
            row["Gplus"] = a.g_plus;
            row["Gminus"] = a.g_minus;
            row["Pplus"] = a.p_plus;
            row["Pminus"] = a.p_minus;
            row["Splus"] = a.s_plus;
            row["Sminus"] = a.s_minus;
            row["f_pi0_effective"] = a.f_pi0_effective;
            row["contamination_ratio_effective"] = a.contamination_ratio_effective;
            row["polarized_denominator"] = a.polarized_denominator;
            row["raw_epg_valid"] = a.raw_valid;
            row["pi0_valid"] = a.pi0_valid;
            if (a.valid) {
                row["BSA_pi0_subtracted"] = a.value;
                row["BSA_pi0_subtracted_stat"] = a.stat;
            }
            if (a.raw_valid) {
                row["BSA_raw_epg"] = a.raw_value;
                row["BSA_raw_epg_stat"] = a.raw_stat;
            }
            if (a.pi0_valid) {
                row["BSA_measured_eppi0"] = a.pi0_value;
                row["BSA_measured_eppi0_stat"] = a.pi0_stat;
            }
            rows_json.push_back(row);
        }

        j["groups"][group] = rows_json;
    }

    const std::filesystem::path p(path);
    std::filesystem::create_directories(p.parent_path());
    std::ofstream fout(path);
    if (!fout.is_open()) {
        fatal("[bsa] cannot write JSON summary: " + path);
    }
    fout << std::setw(2) << j << "\n";
}

struct PlotKey {
    double xBmin = 0.0;
    double xBmax = 0.0;
    double Q2min = 0.0;
    double Q2max = 0.0;
    double tmin = 0.0;
    double tmax = 0.0;

    bool operator<(const PlotKey& o) const {
        if (xBmin != o.xBmin) return xBmin < o.xBmin;
        if (xBmax != o.xBmax) return xBmax < o.xBmax;
        if (Q2min != o.Q2min) return Q2min < o.Q2min;
        if (Q2max != o.Q2max) return Q2max < o.Q2max;
        if (tmin != o.tmin) return tmin < o.tmin;
        return tmax < o.tmax;
    }
};

static void make_bsa_plots(const std::string& output_root,
                           const std::vector<RowBin>& rows,
                           const std::map<std::string, std::vector<AsymResult>>& results) {
    const std::filesystem::path root = std::filesystem::path(output_root) / "bsa_plots";
    std::filesystem::create_directories(root);

    for (const std::string& group : output_group_order()) {
        auto it_group = results.find(group);
        if (it_group == results.end()) {
            continue;
        }
        const auto& vec = it_group->second;

        const std::filesystem::path out_dir = root / sanitize_token(group);
        std::filesystem::create_directories(out_dir);

        std::map<PlotKey, std::vector<int>> row_groups;
        for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
            if (!rows[r].valid || r >= static_cast<int>(vec.size()) || !vec[r].valid) {
                continue;
            }
            PlotKey k{rows[r].xBmin, rows[r].xBmax, rows[r].Q2min, rows[r].Q2max, rows[r].tmin, rows[r].tmax};
            row_groups[k].push_back(r);
        }

        int canvas_index = 0;
        for (const auto& kg : row_groups) {
            const PlotKey& key = kg.first;
            const auto& row_indices = kg.second;
            if (row_indices.empty()) {
                continue;
            }

            std::vector<double> x_sub, y_sub, ex_sub, ey_sub;
            std::vector<double> x_raw, y_raw, ex_raw, ey_raw;
            std::vector<double> x_pi0, y_pi0, ex_pi0, ey_pi0;

            for (int row_index : row_indices) {
                const RowBin& rb = rows[row_index];
                const AsymResult& a = vec[row_index];
                const double xc = phi_center_deg(rb.pmin, rb.pmax);
                const double xe = phi_half_width_deg(rb.pmin, rb.pmax);

                if (a.valid) {
                    x_sub.push_back(xc);
                    y_sub.push_back(a.value);
                    ex_sub.push_back(xe);
                    ey_sub.push_back(a.stat);
                }
                if (a.raw_valid) {
                    x_raw.push_back(xc);
                    y_raw.push_back(a.raw_value);
                    ex_raw.push_back(xe);
                    ey_raw.push_back(a.raw_stat);
                }
                if (a.pi0_valid) {
                    x_pi0.push_back(xc);
                    y_pi0.push_back(a.pi0_value);
                    ex_pi0.push_back(xe);
                    ey_pi0.push_back(a.pi0_stat);
                }
            }

            if (x_sub.empty()) {
                continue;
            }

            TCanvas c("c_bsa", "", 1100, 850);
            c.SetLeftMargin(0.12);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);
            c.SetTopMargin(0.08);
            c.SetTickx(1);
            c.SetTicky(1);

            TGraphErrors gr_sub(static_cast<int>(x_sub.size()), x_sub.data(), y_sub.data(), ex_sub.data(), ey_sub.data());
            gr_sub.SetName("gr_bsa_subtracted");
            gr_sub.SetMarkerStyle(20);
            gr_sub.SetMarkerSize(1.0);
            gr_sub.SetLineWidth(2);
            gr_sub.SetTitle("");
            gr_sub.GetXaxis()->SetTitle("#phi (deg)");
            gr_sub.GetYaxis()->SetTitle("A_{LU}");
            gr_sub.GetXaxis()->SetLimits(0.0, 360.0);
            gr_sub.GetYaxis()->SetRangeUser(-1.0, 1.0);
            gr_sub.Draw("AP");

            TGraphErrors gr_raw;
            if (!x_raw.empty()) {
                gr_raw = TGraphErrors(static_cast<int>(x_raw.size()), x_raw.data(), y_raw.data(), ex_raw.data(), ey_raw.data());
                gr_raw.SetName("gr_bsa_raw");
                gr_raw.SetMarkerStyle(24);
                gr_raw.SetMarkerSize(0.9);
                gr_raw.SetLineWidth(2);
                gr_raw.Draw("P SAME");
            }

            TGraphErrors gr_pi0;
            if (!x_pi0.empty()) {
                gr_pi0 = TGraphErrors(static_cast<int>(x_pi0.size()), x_pi0.data(), y_pi0.data(), ex_pi0.data(), ey_pi0.data());
                gr_pi0.SetName("gr_bsa_pi0");
                gr_pi0.SetMarkerStyle(25);
                gr_pi0.SetMarkerSize(0.9);
                gr_pi0.SetLineWidth(2);
                gr_pi0.Draw("P SAME");
            }

            TF1 fit("fit_sin", "[0]*sin(x*TMath::Pi()/180.0)", 0.0, 360.0);
            fit.SetParameter(0, 0.0);
            fit.SetLineWidth(2);
            gr_sub.Fit(&fit, "Q0");
            fit.Draw("SAME");

            TLegend leg(0.56, 0.70, 0.94, 0.90);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextFont(42);
            leg.SetTextSize(0.030);
            leg.AddEntry(&gr_sub, "#pi^{0}-subtracted ep#gamma", "lep");
            if (!x_raw.empty()) leg.AddEntry(&gr_raw, "raw ep#gamma", "lep");
            if (!x_pi0.empty()) leg.AddEntry(&gr_pi0, "measured ep#pi^{0}", "lep");
            leg.AddEntry(&fit, "A sin#phi fit", "l");
            leg.Draw();

            TLatex lat;
            lat.SetNDC();
            lat.SetTextFont(42);
            lat.SetTextSize(0.032);
            lat.DrawLatex(0.16, 0.90, Form("%s", group.c_str()));
            lat.DrawLatex(0.16, 0.855,
                          Form("x_{B}: %.3g-%.3g, Q^{2}: %.3g-%.3g GeV^{2}",
                               key.xBmin, key.xBmax, key.Q2min, key.Q2max));
            lat.DrawLatex(0.16, 0.815,
                          Form("|t|: %.3g-%.3g GeV^{2}", key.tmin, key.tmax));
            lat.DrawLatex(0.16, 0.775,
                          Form("sin#phi amplitude = %.4f #pm %.4f", fit.GetParameter(0), fit.GetParError(0)));

            const std::string name =
                (out_dir / Form("bsa_%s_xB_%g_%g_Q2_%g_%g_t_%g_%g.png",
                                sanitize_token(group).c_str(),
                                key.xBmin, key.xBmax,
                                key.Q2min, key.Q2max,
                                key.tmin, key.tmax)).string();
            c.SaveAs(name.c_str());
            ++canvas_index;
        }

        std::cout << "[bsa] Wrote " << canvas_index
                  << " phi-dependence plots for group " << group
                  << " to " << out_dir.string() << "\n";
    }
}

} // namespace

bool update_bsa_counts_csv(const std::map<std::string, TTree*>& dvcsDataTrees,
                           const std::map<std::string, TTree*>& eppi0DataTrees,
                           const BSAOptions& options) {
    try {
        CSV csv;
        load_csv(options.csv_path, csv);

        const std::vector<RowBin> rows = load_row_bins_from_csv(csv);
        const FastBinning fast_bins = build_fast_binning(rows);
        const TopoCutMap sigma_cuts = load_combined_cuts(options.combined_cuts_json);
        const PeriodLeakRows leaks = build_pi0_leakage(csv, options);

        const PeriodCounts gamma_counts = accumulate_counts(dvcsDataTrees, "DVCS", rows, fast_bins, sigma_cuts);
        const PeriodCounts pi0_counts = accumulate_counts(eppi0DataTrees, "eppi0", rows, fast_bins, sigma_cuts);

        std::map<std::string, std::vector<AsymResult>> results;
        for (const std::string& group : output_group_order()) {
            const std::vector<std::string> components = component_periods_for_group(group);
            std::vector<AsymResult> group_results(rows.size());

            for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
                if (!rows[r].valid) {
                    continue;
                }
                group_results[r] = compute_group_bsa(gamma_counts, pi0_counts, leaks, components, r, options);
            }

            results[group] = std::move(group_results);
        }

        for (const std::string& group : output_group_order()) {
            const std::string col_name = "BSA, counts, " + group;
            const int c = col_strict(csv, col_name);
            const auto& vec = results.at(group);
            for (int r = 0; r < static_cast<int>(csv.rows.size()); ++r) {
                if (r >= static_cast<int>(vec.size()) || !vec[r].valid) {
                    csv.rows[r][c].clear();
                    continue;
                }
                csv.rows[r][c] = fmt_tuple(vec[r].value, vec[r].stat);
            }
        }

        write_csv_atomic(options.csv_path, csv);
        std::cout << "[bsa] Updated BSA count columns in " << options.csv_path << "\n";

        const std::filesystem::path json_path =
            std::filesystem::path(options.output_root) / "jsons" / "BSA_counts" / "bsa_counts_summary.json";
        write_json_summary(json_path.string(), rows, results);
        std::cout << "[bsa] Wrote JSON summary to " << json_path.string() << "\n";

        if (options.make_plots) {
            make_bsa_plots(options.output_root, rows, results);
        }

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[bsa] ERROR: " << e.what() << "\n";
        return false;
    }
}
