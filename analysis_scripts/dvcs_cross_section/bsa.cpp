// bsa.cpp
// -----------------------------------------------------------------------------
// Pi0-subtracted, direct count-based beam-spin asymmetry stage for the DVCS
// pass-2 workflow.
//
// This stage deliberately stays at the measured-count level. It does not use
// acceptance unfolding, bin migration unfolding, radiative corrections, or
// model bin-centering factors. It does apply:
//   * the process-wide global cuts from global_cuts.cpp,
//   * the data-derived exclusivity windows in output/jsons/combined_cuts.json,
//   * helicity-separated measured ep->eppi0 subtraction using the existing
//     bin-by-bin contamination ratio and normalized-yield CSV columns.
//
// Per period and CSV row:
//   G+/- = measured ep->epgamma counts after cuts
//   P+/- = measured ep->eppi0 counts after cuts
//   fpi0 = contamination_ratio * Y_norm(epgamma) / Y_norm(eppi0)
//   S+/- = G+/- - fpi0 * P+/-
//   A_LU = (S+ - S-) / [Pbeam * (S+ + S-)]
//
// For combined groups with multiple beam polarizations:
//   A_LU = sum_i(S_i+ - S_i-) / sum_i[P_i * (S_i+ + S_i-)]
//
// Plots are xB-matrix canvases: rows are Q2 bins, columns are |t| bins and
// each pad is A_LU(phi). Data are fit to A sin(phi)/(1 + B cos(phi)).
// ROOT tree loops and ROOT plotting are intentionally serial.
// -----------------------------------------------------------------------------

#include "bsa.h"
#include "global_cuts.h"

// ROOT
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVirtualPad.h>

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
#include <memory>
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
    } //endfor
    return s;
}

static inline std::string sanitize_token(std::string s) {
    for (char& c : s) {
        if (!std::isalnum(static_cast<unsigned char>(c))) {
            c = '_';
        } //endif
    } //endfor
    return s;
}

static inline double wrap_phi_deg(double phi_deg) {
    double p = std::fmod(phi_deg, 360.0);
    if (p < 0.0) {
        p += 360.0;
    } //endif
    if (p >= 360.0) {
        p = std::nextafter(360.0, 0.0);
    } //endif
    return p;
}

static inline double delta_phi_rad_from_two_phi(double phi_a, double phi_b) {
    double d = std::fmod(phi_a - phi_b, 2.0 * PI);
    if (d <= -PI) {
        d += 2.0 * PI;
    } //endif
    if (d > PI) {
        d -= 2.0 * PI;
    } //endif
    return std::fabs(d);
}

static inline bool in_range(double v, double a, double b) {
    return (v >= a) && (v < b);
}

static inline bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (pmax_deg > pmin_deg) {
        return in_range(phi_deg, pmin_deg, pmax_deg);
    } //endif
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg);
}

static inline double phi_center_deg(double pmin_deg, double pmax_deg) {
    double lo = wrap_phi_deg(pmin_deg);
    double hi = wrap_phi_deg(pmax_deg);
    if (hi <= lo) {
        hi += 360.0;
    } //endif
    return wrap_phi_deg(0.5 * (lo + hi));
}

static inline double phi_half_width_deg(double pmin_deg, double pmax_deg) {
    double lo = wrap_phi_deg(pmin_deg);
    double hi = wrap_phi_deg(pmax_deg);
    if (hi <= lo) {
        hi += 360.0;
    } //endif
    return 0.5 * (hi - lo);
}

static inline std::string topo_dir(int det1, int det2) {
    if (det1 == 1 && det2 == 1) return "FD_FD";
    if (det1 == 2 && det2 == 1) return "CD_FD";
    if (det1 == 2 && det2 == 0) return "CD_FT";
    return "";
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
    fatal("[bsa] unknown period label: " + period);
    return "";
}

static inline std::string period_global_cuts_label_from_display(const std::string& period) {
    if (period == "Fa18 Inb") return "fa18_inb";
    if (period == "Fa18 Out") return "fa18_out";
    if (period == "Sp19 Inb") return "sp19_inb";
    if (period == "Sp18 Inb") return "sp18_inb";
    if (period == "Sp18 Out") return "sp18_out";
    fatal("[bsa] unknown period label for global cuts: " + period);
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
    fatal("[bsa] unknown BSA output group: " + group);
    return {};
}

static const std::vector<std::string>& csv_topology_labels() {
    static const std::vector<std::string> v = {"(FD, FD)", "(CD, FD)", "(CD, FT)"};
    return v;
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
            } //endif
        } else if (c == ',' && !inq) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        } //endif
    } //endfor

    out.push_back(cur);
    return out;
}

static bool load_csv(const std::string& path, CSV& csv) {
    std::ifstream fin(path);
    if (!fin.is_open()) fatal("[bsa] cannot open CSV: " + path);

    std::string line;
    if (!std::getline(fin, line)) fatal("[bsa] empty CSV: " + path);

    csv.header = split_csv_line(line);
    csv.index.clear();
    for (int i = 0; i < static_cast<int>(csv.header.size()); ++i) {
        if (csv.index.count(csv.header[i])) fatal("[bsa] duplicate CSV column: " + csv.header[i]);
        csv.index[csv.header[i]] = i;
    } //endfor

    csv.rows.clear();
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = split_csv_line(line);
        if (row.size() < csv.header.size()) row.resize(csv.header.size(), "");
        if (row.size() != csv.header.size()) fatal("[bsa] CSV row width mismatch in " + path);
        csv.rows.push_back(std::move(row));
    } //endwhile

    return true;
}

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";
    std::ofstream fout(tmp);
    if (!fout.is_open()) fatal("[bsa] cannot write temporary CSV: " + tmp);

    auto write_cell = [&](const std::string& s) {
        const bool quote =
            s.find(',') != std::string::npos ||
            s.find('"') != std::string::npos ||
            s.find('\n') != std::string::npos ||
            s.find('\r') != std::string::npos;

        if (!quote) {
            fout << s;
            return;
        } //endif

        fout << '"';
        for (char c : s) {
            if (c == '"') fout << "\"\"";
            else fout << c;
        } //endfor
        fout << '"';
    };

    for (size_t i = 0; i < csv.header.size(); ++i) {
        write_cell(csv.header[i]);
        if (i + 1 < csv.header.size()) fout << ',';
    } //endfor
    fout << '\n';

    for (const auto& row : csv.rows) {
        if (row.size() != csv.header.size()) fatal("[bsa] CSV row width mismatch while writing");
        for (size_t i = 0; i < row.size(); ++i) {
            write_cell(row[i]);
            if (i + 1 < row.size()) fout << ',';
        } //endfor
        fout << '\n';
    } //endfor

    fout.close();
    if (!fout) fatal("[bsa] failed while writing temporary CSV: " + tmp);

    (void)std::remove(path.c_str());
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        fatal("[bsa] failed to rename " + tmp + " to " + path);
    } //endif
}

static int col_strict(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);
    if (it == csv.index.end()) fatal("[bsa] missing required CSV column: " + name);
    return it->second;
}

static int col_optional(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);
    if (it == csv.index.end()) return -1;
    return it->second;
}

static inline double to_double_strict(const std::string& s, const std::string& what) {
    if (s.empty()) fatal("[bsa] empty numeric cell for " + what);
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str()) fatal("[bsa] parse failure for " + what + " value '" + s + "'");
    return v;
}

static inline bool parse_double_loose(const std::string& s, double& out) {
    if (s.empty()) return false;
    char* end = nullptr;
    out = std::strtod(s.c_str(), &end);
    return end != s.c_str() && std::isfinite(out);
}

static inline bool to_bool_valid(const std::string& s) {
    return (s == "1" || s == "1.0" || s == "true" || s == "TRUE");
}

static bool parse_tuple_first(const std::string& s, double& value) {
    std::string t = s;
    t.erase(std::remove_if(t.begin(), t.end(), [](unsigned char c) {
        return std::isspace(c);
    }), t.end());

    if (t.empty()) return false;

    if (t.front() == '(') {
        const size_t comma = t.find(',');
        const size_t close = t.find(')');
        if (comma == std::string::npos || close == std::string::npos || comma <= 1) return false;
        return parse_double_loose(t.substr(1, comma - 1), value);
    } //endif

    return parse_double_loose(t, value);
}

// -----------------------------------------------------------------------------
// Row binning
// -----------------------------------------------------------------------------

struct RowBin {
    double xBmin = std::numeric_limits<double>::quiet_NaN();
    double xBmax = std::numeric_limits<double>::quiet_NaN();
    double Q2min = std::numeric_limits<double>::quiet_NaN();
    double Q2max = std::numeric_limits<double>::quiet_NaN();
    double tmin  = std::numeric_limits<double>::quiet_NaN();
    double tmax  = std::numeric_limits<double>::quiet_NaN();
    double pmin  = std::numeric_limits<double>::quiet_NaN();
    double pmax  = std::numeric_limits<double>::quiet_NaN();
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
    const int c_tmin  = col_strict(csv, "t_abs_min");
    const int c_tmax  = col_strict(csv, "t_abs_max");
    const int c_pmin  = col_strict(csv, "phimin");
    const int c_pmax  = col_strict(csv, "phimax");
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
        b.tmin  = to_double_strict(row[c_tmin],  "t_abs_min");
        b.tmax  = to_double_strict(row[c_tmax],  "t_abs_max");
        b.pmin  = to_double_strict(row[c_pmin],  "phimin");
        b.pmax  = to_double_strict(row[c_pmax],  "phimax");
        b.valid = to_bool_valid(row[c_valid]);
        rows.push_back(b);
    } //endfor

    return rows;
}

static inline bool axis_bin_equal(const AxisBin& a, const AxisBin& b) {
    return a.min == b.min && a.max == b.max;
}

static void add_unique_axis_bin(std::vector<AxisBin>& bins, double minv, double maxv) {
    AxisBin b{minv, maxv};
    auto it = std::find_if(bins.begin(), bins.end(), [&](const AxisBin& x) {
        return axis_bin_equal(x, b);
    });
    if (it == bins.end()) bins.push_back(b);
}

static void sort_axis_bins(std::vector<AxisBin>& bins) {
    std::sort(bins.begin(), bins.end(), [](const AxisBin& a, const AxisBin& b) {
        if (a.min != b.min) return a.min < b.min;
        return a.max < b.max;
    });
}

static int find_axis_bin_index(const std::vector<AxisBin>& bins, double value) {
    for (int i = 0; i < static_cast<int>(bins.size()); ++i) {
        if (value >= bins[i].min && value < bins[i].max) return i;
    } //endfor
    return -1;
}

static int find_axis_bin_exact(const std::vector<AxisBin>& bins, double minv, double maxv) {
    for (int i = 0; i < static_cast<int>(bins.size()); ++i) {
        if (bins[i].min == minv && bins[i].max == maxv) return i;
    } //endfor
    return -1;
}

static FastBinning build_fast_binning(const std::vector<RowBin>& rows) {
    FastBinning fb;

    for (const RowBin& r : rows) {
        if (!r.valid) continue;
        add_unique_axis_bin(fb.xbins, r.xBmin, r.xBmax);
        add_unique_axis_bin(fb.qbins, r.Q2min, r.Q2max);
        add_unique_axis_bin(fb.tbins, r.tmin, r.tmax);
    } //endfor

    sort_axis_bins(fb.xbins);
    sort_axis_bins(fb.qbins);
    sort_axis_bins(fb.tbins);

    fb.rows_by_xqt.resize(fb.xbins.size());
    for (size_t ix = 0; ix < fb.xbins.size(); ++ix) {
        fb.rows_by_xqt[ix].resize(fb.qbins.size());
        for (size_t iq = 0; iq < fb.qbins.size(); ++iq) {
            fb.rows_by_xqt[ix][iq].resize(fb.tbins.size());
        } //endfor
    } //endfor

    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        const RowBin& row = rows[r];
        if (!row.valid) continue;

        const int ix = find_axis_bin_exact(fb.xbins, row.xBmin, row.xBmax);
        const int iq = find_axis_bin_exact(fb.qbins, row.Q2min, row.Q2max);
        const int it = find_axis_bin_exact(fb.tbins, row.tmin, row.tmax);
        if (ix < 0 || iq < 0 || it < 0) fatal("[bsa] failed to build fast row-bin lookup");
        fb.rows_by_xqt[ix][iq][it].push_back(r);
    } //endfor

    std::cout << "[bsa] Fast bin lookup built with "
              << fb.xbins.size() << " xB bins, "
              << fb.qbins.size() << " Q2 bins, "
              << fb.tbins.size() << " |t| bins.\n";

    return fb;
}

// -----------------------------------------------------------------------------
// Combined cuts loader
// -----------------------------------------------------------------------------

struct SigmaStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std  = std::numeric_limits<double>::quiet_NaN();
    double cut_low = std::numeric_limits<double>::quiet_NaN();
    double cut_high = std::numeric_limits<double>::quiet_NaN();
    double quantile = 0.0;
    std::string mode = "symmetric_3sigma";
};

using CutVarMap = std::unordered_map<std::string, SigmaStats>;
using TopoCutMap = std::unordered_map<std::string, CutVarMap>;

static inline bool within_cut_window(double v, const SigmaStats& s) {
    if (!std::isfinite(v)) return false;

    if (s.mode == "upper_quantile") {
        if (!std::isfinite(s.cut_high)) return true;
        return v <= s.cut_high;
    } //endif

    double lo = s.cut_low;
    double hi = s.cut_high;
    if (!(std::isfinite(lo) && std::isfinite(hi)) || hi <= lo) {
        if (!std::isfinite(s.mean) || !std::isfinite(s.std) || s.std <= 0.0) return true;
        lo = s.mean - 3.0 * s.std;
        hi = s.mean + 3.0 * s.std;
    } //endif

    return v >= lo && v <= hi;
}

static TopoCutMap load_combined_cuts(const std::string& combined_cuts_json) {
    std::ifstream fin(combined_cuts_json);
    if (!fin.is_open()) fatal("[bsa] cannot open combined cuts JSON: " + combined_cuts_json);

    nlohmann::json j;
    fin >> j;
    if (!j.is_object()) fatal("[bsa] combined cuts JSON is not an object: " + combined_cuts_json);

    TopoCutMap out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const auto& block = it.value();
        if (!block.is_object() || !block.contains("data")) continue;

        const auto& data = block["data"];
        if (!data.is_object()) continue;

        CutVarMap vm;
        for (auto vit = data.begin(); vit != data.end(); ++vit) {
            const std::string var = vit.key();
            const auto& stats = vit.value();
            if (!stats.is_object() || !stats.contains("mean") || !stats.contains("std")) continue;

            SigmaStats s;
            try {
                s.mean = stats["mean"].get<double>();
                s.std  = stats["std"].get<double>();
                if (stats.contains("cut_low"))  s.cut_low  = stats["cut_low"].get<double>();
                if (stats.contains("cut_high")) s.cut_high = stats["cut_high"].get<double>();
                if (stats.contains("quantile")) s.quantile = stats["quantile"].get<double>();
                if (stats.contains("mode"))     s.mode     = stats["mode"].get<std::string>();
            } catch (...) {
                continue;
            } //endtry

            if (!std::isfinite(s.cut_low) || !std::isfinite(s.cut_high) || s.cut_high <= s.cut_low) {
                if (std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0) {
                    s.cut_low  = s.mean - 3.0 * s.std;
                    s.cut_high = s.mean + 3.0 * s.std;
                } //endif
            } //endif

            if (std::isfinite(s.cut_high)) vm.emplace(var, s);
        } //endfor

        if (!vm.empty()) out.emplace(key, std::move(vm));
    } //endfor

    std::cout << "[bsa] Loaded " << out.size()
              << " data topology/period cut blocks from "
              << combined_cuts_json << "\n";

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
    double Delta_phi = 0.0; bool has_Delta_phi = false;

    double open_angle_ep2 = 0.0; bool has_open_angle = false;
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
            if (t->GetBranch(name)) t->SetBranchStatus(name, 1);
        };

        const char* names[] = {
            "runnum", "detector1", "detector2", "helicity",
            "x", "Q2", "t1", "phi2", "Delta_phi",
            "open_angle_ep2", "pTmiss", "Emiss2", "Mx2", "Mx2_1", "Mx2_2",
            "xF", "theta_gamma_gamma", "theta_pi0_pi0",
            "e_p", "e_theta", "e_phi", "p1_theta", "p1_phi",
            "p2_p", "p2_theta", "p2_phi"
        };
        for (const char* name : names) enable(name);

        t->SetCacheSize(0);

        auto bind_int = [&](const char* name, int* ptr, bool& has) {
            if (t->GetBranch(name)) { t->SetBranchAddress(name, ptr); has = true; }
        };
        auto bind_double = [&](const char* name, double* ptr, bool& has) {
            if (t->GetBranch(name)) { t->SetBranchAddress(name, ptr); has = true; }
        };

        bind_int("runnum", &runnum, has_runnum);
        bind_int("detector1", &detector1, has_detector1);
        bind_int("detector2", &detector2, has_detector2);
        bind_int("helicity", &helicity, has_helicity);

        bind_double("x", &x, has_x);
        bind_double("Q2", &Q2, has_Q2);
        bind_double("t1", &t1, has_t1);
        bind_double("phi2", &phi2, has_phi2);
        bind_double("Delta_phi", &Delta_phi, has_Delta_phi);
        bind_double("open_angle_ep2", &open_angle_ep2, has_open_angle);
        bind_double("pTmiss", &pTmiss, has_pTmiss);
        bind_double("Emiss2", &Emiss2, has_Emiss2);
        bind_double("Mx2", &Mx2, has_Mx2);
        bind_double("Mx2_1", &Mx2_1, has_Mx2_1);
        bind_double("Mx2_2", &Mx2_2, has_Mx2_2);
        bind_double("xF", &xF, has_xF);
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
        if (has_Delta_phi) { has_val = true; return Delta_phi; }
        if (has_p1_phi && has_p2_phi) {
            has_val = true;
            return delta_phi_rad_from_two_phi(p1_phi, p2_phi);
        } //endif
        has_val = false;
        return 0.0;
    }
};

static inline double branch_value_for_sigma_var(const BranchBinder& b,
                                                const std::string& var,
                                                bool& has_val) {
    has_val = true;
    if (var == "Emiss2") { has_val = b.has_Emiss2; return b.Emiss2; }
    if (var == "Mx2") { has_val = b.has_Mx2; return b.Mx2; }
    if (var == "Mx2_1") { has_val = b.has_Mx2_1; return b.Mx2_1; }
    if (var == "Mx2_2") { has_val = b.has_Mx2_2; return b.Mx2_2; }
    if (var == "Delta_phi") { return b.delta_phi_value(has_val); }
    if (var == "pTmiss") { has_val = b.has_pTmiss; return b.pTmiss; }
    if (var == "xF") { has_val = b.has_xF; return b.xF; }
    if (var == "theta_gamma_gamma") { has_val = b.has_theta_gamma_gamma; return b.theta_gamma_gamma; }
    if (var == "theta_pi0_pi0") { has_val = b.has_theta_pi0_pi0; return b.theta_pi0_pi0; }
    has_val = false;
    return 0.0;
}

static bool passes_sigma_cuts(const TopoCutMap& cuts,
                              const std::string& key,
                              const BranchBinder& b) {
    auto it = cuts.find(key);
    if (it == cuts.end()) fatal("[bsa] missing data cut key in combined_cuts.json: " + key);

    for (const auto& var_pair : it->second) {
        bool has_val = false;
        const double val = branch_value_for_sigma_var(b, var_pair.first, has_val);
        if (!has_val) {
            fatal("[bsa] cut key " + key + " requires missing branch: " + var_pair.first);
        } //endif
        if (!within_cut_window(val, var_pair.second)) return false;
    } //endfor

    return true;
}

static bool passes_global_cuts_dispatch(const BranchBinder& b,
                                        const std::string& period_display_label) {
    const GlobalCutConfig& cfg = default_global_cuts();
    const std::string period_label = period_global_cuts_label_from_display(period_display_label);

    if (!(b.has_t1 && b.has_open_angle)) return false;
    if (cfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (b.has_runnum && is_excluded_run(b.runnum)) return false;

    if (cfg.enable_topology_filter || global_cuts_require_sector_phi(cfg) || cfg.enable_dvcsgen_ycol_cut ||
        global_cuts_require_auxiliary_kinematics(cfg)) {
        if (!(b.has_detector1 && b.has_detector2)) {
            fatal("[bsa] topology/sector/ycol/aux selection requires detector1 and detector2 branches");
        } //endif
    } //endif

    if (global_cuts_require_sector_phi(cfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            fatal("[bsa] sector selection requires e_phi, p1_phi and p2_phi branches");
        } //endif
    } //endif

    if (global_cuts_require_auxiliary_kinematics(cfg)) {
        if (!(b.has_e_theta && b.has_e_phi && b.has_p1_theta && b.has_p1_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[bsa] auxiliary fiducial cuts require e_theta, e_phi, p1_theta, p1_phi, p2_p, p2_theta and p2_phi branches");
        } //endif
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p1_theta, b.p1_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    } //endif

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi && b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[bsa] dvcsgen ycol cut requires e_p, e_theta, e_phi, p2_p, p2_theta and p2_phi branches");
        } //endif
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    } //endif

    if (global_cuts_require_sector_phi(cfg)) {
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  b.e_phi, b.p1_phi, b.p2_phi,
                                  cfg);
    } //endif

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                              b.detector1, b.detector2,
                              cfg);
}

// -----------------------------------------------------------------------------
// Counts
// -----------------------------------------------------------------------------

struct HelCounts {
    double plus = 0.0;
    double minus = 0.0;
};

using RowCounts = std::unordered_map<int, HelCounts>;
using PeriodCounts = std::unordered_map<std::string, RowCounts>;

static inline void add_event(HelCounts& h, int helicity) {
    if (helicity > 0) h.plus += 1.0;
    else if (helicity < 0) h.minus += 1.0;
}

static PeriodCounts accumulate_counts(const std::map<std::string, TTree*>& trees,
                                      const std::string& channel_name,
                                      const std::string& cut_prefix,
                                      const std::vector<RowBin>& rows,
                                      const FastBinning& fast_bins,
                                      const TopoCutMap& sigma_cuts) {
    PeriodCounts out;
    for (const std::string& period : base_period_order()) out[period] = RowCounts();

    for (const auto& kv : trees) {
        const std::string& tree_key = kv.first;
        TTree* tree = kv.second;
        if (!tree) continue;

        const std::string period = period_display_from_tree_key(tree_key);
        if (period.empty()) {
            std::cout << "[bsa] Skipping non-canonical or supplemental "
                      << channel_name << " data tree key: " << tree_key << "\n";
            continue;
        } //endif

        BranchBinder b;
        b.bind(tree);

        if (!(b.has_detector1 && b.has_detector2 && b.has_helicity && b.has_x && b.has_Q2 && b.has_t1 && b.has_phi2)) {
            fatal("[bsa] " + channel_name + " tree " + tree_key +
                  " is missing one or more required branches: detector1, detector2, helicity, x, Q2, t1, phi2");
        } //endif

        const std::string period_code = period_code_from_display(period);
        const Long64_t n_entries = tree->GetEntries();

        long long n_topology = 0;
        long long n_global = 0;
        long long n_sigma = 0;
        long long n_matched = 0;

        for (Long64_t i = 0; i < n_entries; ++i) {
            tree->GetEntry(i);

            const std::string topo = topo_dir(b.detector1, b.detector2);
            if (topo.empty()) continue;
            ++n_topology;

            if (!passes_global_cuts_dispatch(b, period)) continue;
            ++n_global;

            const std::string cut_key = cut_prefix + "_" + period_code + "_" + topo;
            if (!passes_sigma_cuts(sigma_cuts, cut_key, b)) continue;
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
                if (!row_accepts_phi(phi_deg, rb.pmin, rb.pmax)) continue;
                add_event(out[period][row_index], b.helicity);
                matched = true;
            } //endfor

            if (matched) ++n_matched;
        } //endfor

        std::cout << "[bsa] channel=" << channel_name
                  << " tree=" << tree_key
                  << " period=" << period
                  << " entries=" << static_cast<long long>(n_entries)
                  << " topology=" << n_topology
                  << " global=" << n_global
                  << " sigma=" << n_sigma
                  << " matched=" << n_matched
                  << "\n";
    } //endfor

    return out;
}

// -----------------------------------------------------------------------------
// Pi0 leakage factors
// -----------------------------------------------------------------------------

struct PeriodLeakRows {
    std::unordered_map<std::string, std::vector<double>> f_pi0;
    std::unordered_map<std::string, std::vector<double>> contamination;
};

static std::string first_existing_col(const CSV& csv, const std::vector<std::string>& names) {
    for (const std::string& n : names) {
        if (csv.index.count(n)) return n;
    } //endfor
    return "";
}

static bool sum_topology_columns_for_period_channel_unpol_preferred(const CSV& csv,
                                                                    int row_index,
                                                                    const std::string& channel,
                                                                    const std::string& period,
                                                                    double& sum_value) {
    sum_value = 0.0;
    bool found_any = false;

    auto try_helicity_set = [&](const std::vector<std::string>& helicities) {
        double local_sum = 0.0;
        bool local_found = false;

        for (const std::string& topo : csv_topology_labels()) {
            for (const std::string& hel : helicities) {
                const std::vector<std::string> candidates = {
                    "normalized raw yield, " + channel + ", " + topo + ", exp, " + period + ", " + hel,
                    "raw yield, " + channel + ", " + topo + ", exp, " + period + ", " + hel
                };

                for (const std::string& cname : candidates) {
                    const int c = col_optional(csv, cname);
                    if (c < 0) continue;

                    double v = 0.0;
                    if (parse_tuple_first(csv.rows[row_index][c], v)) {
                        local_sum += v;
                        local_found = true;
                        break;
                    } //endif
                } //endfor
            } //endfor
        } //endfor

        if (local_found) {
            sum_value = local_sum;
            found_any = true;
        } //endif
    };

    try_helicity_set({"unpol"});
    if (!found_any) try_helicity_set({"pos", "neg"});

    return found_any;
}

static bool period_has_any_channel_yield_columns(const CSV& csv,
                                                 const std::string& channel,
                                                 const std::string& period) {
    for (const std::string& topo : csv_topology_labels()) {
        for (const std::string& hel : {std::string("unpol"), std::string("pos"), std::string("neg")}) {
            const std::vector<std::string> candidates = {
                "normalized raw yield, " + channel + ", " + topo + ", exp, " + period + ", " + hel,
                "raw yield, " + channel + ", " + topo + ", exp, " + period + ", " + hel
            };
            for (const std::string& cname : candidates) {
                if (col_optional(csv, cname) >= 0) return true;
            } //endfor
        } //endfor
    } //endfor
    return false;
}

static PeriodLeakRows build_pi0_leakage(const CSV& csv, const BSAOptions& opt) {
    PeriodLeakRows out;

    for (const std::string& period : base_period_order()) {
        const std::string c_contam_name = first_existing_col(csv, {
            "contamination ratio, " + period,
            "pi0 contamination ratio, " + period,
            "pi0 contamination, " + period
        });

        if (c_contam_name.empty()) {
            fatal("[bsa] could not find contamination-ratio CSV column for period " + period);
        } //endif
        if (!period_has_any_channel_yield_columns(csv, "ep->epg", period)) {
            fatal("[bsa] could not find any epgamma normalized/raw yield CSV columns for period " + period +
                  ". Expected columns like 'normalized raw yield, ep->epg, (FD, FD), exp, " + period + ", unpol'.");
        } //endif
        if (!period_has_any_channel_yield_columns(csv, "ep->eppi0", period)) {
            fatal("[bsa] could not find any eppi0 normalized/raw yield CSV columns for period " + period +
                  ". Expected columns like 'normalized raw yield, ep->eppi0, (FD, FD), exp, " + period + ", unpol'.");
        } //endif

        const int c_contam = col_strict(csv, c_contam_name);

        std::vector<double> f(csv.rows.size(), 0.0);
        std::vector<double> c(csv.rows.size(), 0.0);

        int n_valid = 0;
        for (int r = 0; r < static_cast<int>(csv.rows.size()); ++r) {
            double contam = 0.0;
            double epg_yield = 0.0;
            double pi0_yield = 0.0;

            const bool ok_c = parse_tuple_first(csv.rows[r][c_contam], contam);
            const bool ok_g = sum_topology_columns_for_period_channel_unpol_preferred(csv, r, "ep->epg", period, epg_yield);
            const bool ok_p = sum_topology_columns_for_period_channel_unpol_preferred(csv, r, "ep->eppi0", period, pi0_yield);

            if (ok_c && ok_g && ok_p && contam >= 0.0 && epg_yield > 0.0 && pi0_yield > 0.0) {
                c[r] = contam;
                f[r] = contam * epg_yield / pi0_yield;
                ++n_valid;
            } //endif

            if (!opt.enable_pi0_subtraction) {
                c[r] = 0.0;
                f[r] = 0.0;
            } //endif
        } //endfor

        out.f_pi0[period] = std::move(f);
        out.contamination[period] = std::move(c);

        std::cout << "[bsa] pi0 leakage scale factors for " << period
                  << ": valid rows=" << n_valid << "/" << csv.rows.size()
                  << " using contamination column '" << c_contam_name
                  << "' and summed topology normalized raw-yield columns.\n";
    } //endfor

    return out;
}

// -----------------------------------------------------------------------------
// BSA math
// -----------------------------------------------------------------------------

struct AsymResult {
    bool valid = false;
    double value = 0.0;
    double stat = 0.0;

    bool raw_valid = false;
    double raw_value = 0.0;
    double raw_stat = 0.0;

    bool pi0_valid = false;
    double pi0_value = 0.0;
    double pi0_stat = 0.0;

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

struct ComponentForVariance {
    double pol = 0.0;
    double s_plus = 0.0;
    double s_minus = 0.0;
    double var_s_plus = 0.0;
    double var_s_minus = 0.0;
};

static HelCounts get_counts_for_row(const PeriodCounts& counts,
                                    const std::string& period,
                                    int row_index) {
    auto ip = counts.find(period);
    if (ip == counts.end()) return {};
    auto ir = ip->second.find(row_index);
    if (ir == ip->second.end()) return {};
    return ir->second;
}

static double get_vector_value(const std::unordered_map<std::string, std::vector<double>>& m,
                               const std::string& period,
                               int row_index) {
    auto ip = m.find(period);
    if (ip == m.end()) return 0.0;
    if (row_index < 0 || row_index >= static_cast<int>(ip->second.size())) return 0.0;
    const double v = ip->second[row_index];
    if (!std::isfinite(v) || v < 0.0) return 0.0;
    return v;
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

    std::vector<ComponentForVariance> comps;

    for (const std::string& period : component_periods) {
        const HelCounts g = get_counts_for_row(gamma_counts, period, row_index);
        const HelCounts p = get_counts_for_row(pi0_counts, period, row_index);
        const double pol = beam_pol_for_period(period, opt);
        const double f = opt.enable_pi0_subtraction ? get_vector_value(leaks.f_pi0, period, row_index) : 0.0;
        const double c = opt.enable_pi0_subtraction ? get_vector_value(leaks.contamination, period, row_index) : 0.0;

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
            raw_var_num += std::max(0.0, raw_s - raw_d * raw_d / raw_s);
        } //endif

        const double pi0_s = pp + pm;
        const double pi0_d = pp - pm;
        if (pi0_s > 0.0) {
            pi0_num += pi0_d;
            pi0_den += pol * pi0_s;
            pi0_var_num += std::max(0.0, pi0_s - pi0_d * pi0_d / pi0_s);
        } //endif

        const double sig_s = sp + sm;
        const double sig_d = sp - sm;
        if (sig_s > 0.0) {
            numerator += sig_d;
            denominator += pol * sig_s;
            comps.push_back(ComponentForVariance{pol, sp, sm, gp + f * f * pp, gm + f * f * pm});
            weighted_f_num += f * sig_s;
            weighted_c_num += c * sig_s;
            weight_sum += sig_s;
        } //endif
    } //endfor

    if (raw_den > 0.0) {
        r.raw_value = raw_num / raw_den;
        r.raw_stat = std::sqrt(std::max(0.0, raw_var_num)) / raw_den;
        r.raw_valid = std::isfinite(r.raw_value) && std::isfinite(r.raw_stat);
    } //endif

    if (pi0_den > 0.0) {
        r.pi0_value = pi0_num / pi0_den;
        r.pi0_stat = std::sqrt(std::max(0.0, pi0_var_num)) / pi0_den;
        r.pi0_valid = std::isfinite(r.pi0_value) && std::isfinite(r.pi0_stat);
    } //endif

    r.polarized_denominator = denominator;
    if (weight_sum > 0.0) {
        r.f_pi0_effective = weighted_f_num / weight_sum;
        r.contamination_ratio_effective = weighted_c_num / weight_sum;
    } //endif

    if (!(denominator > 0.0)) return r;

    r.value = numerator / denominator;
    for (const ComponentForVariance& comp : comps) {
        const double dA_dSp = (denominator - numerator * comp.pol) / (denominator * denominator);
        const double dA_dSm = (-denominator - numerator * comp.pol) / (denominator * denominator);
        variance += dA_dSp * dA_dSp * comp.var_s_plus;
        variance += dA_dSm * dA_dSm * comp.var_s_minus;
    } //endfor

    r.stat = std::sqrt(std::max(0.0, variance));
    r.valid = std::isfinite(r.value) && std::isfinite(r.stat);
    return r;
}

static std::string fmt_tuple(double value, double stat) {
    if (!(std::isfinite(value) && std::isfinite(stat))) return "";
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
            } //endif
            if (a.raw_valid) {
                row["BSA_raw_epg"] = a.raw_value;
                row["BSA_raw_epg_stat"] = a.raw_stat;
            } //endif
            if (a.pi0_valid) {
                row["BSA_measured_eppi0"] = a.pi0_value;
                row["BSA_measured_eppi0_stat"] = a.pi0_stat;
            } //endif

            rows_json.push_back(row);
        } //endfor

        j["groups"][group] = rows_json;
    } //endfor

    const std::filesystem::path p(path);
    std::filesystem::create_directories(p.parent_path());

    std::ofstream fout(path);
    if (!fout.is_open()) fatal("[bsa] cannot write JSON summary: " + path);
    fout << std::setw(2) << j << "\n";
}

struct BinRange {
    double lo = 0.0;
    double hi = 0.0;

    bool operator<(const BinRange& o) const {
        if (lo != o.lo) return lo < o.lo;
        return hi < o.hi;
    }
};

static std::string range_token(const char* prefix, const BinRange& r) {
    return Form("%s_%g_%g", prefix, r.lo, r.hi);
}

static bool same_range(double lo1, double hi1, const BinRange& r) {
    return lo1 == r.lo && hi1 == r.hi;
}

static inline double range_center(const BinRange& r) {
    return 0.5 * (r.lo + r.hi);
}

static void configure_pad_axes(TGraphErrors& gr, bool left_col, bool bottom_row) {
    gr.GetXaxis()->SetLimits(0.0, 360.0);
    gr.GetYaxis()->SetRangeUser(-1.0, 1.0);
    gr.GetXaxis()->SetTitle(bottom_row ? "#phi (deg)" : "");
    gr.GetYaxis()->SetTitle(left_col ? "A_{LU}" : "");
    gr.GetXaxis()->SetTitleSize(bottom_row ? 0.075 : 0.0);
    gr.GetYaxis()->SetTitleSize(left_col ? 0.080 : 0.0);
    gr.GetXaxis()->SetLabelSize(bottom_row ? 0.060 : 0.0);
    gr.GetYaxis()->SetLabelSize(left_col ? 0.058 : 0.0);
    gr.GetXaxis()->SetTitleOffset(0.95);
    gr.GetYaxis()->SetTitleOffset(0.92);
    gr.GetXaxis()->SetNdivisions(505);
    gr.GetYaxis()->SetNdivisions(505);
}

static void draw_empty_frame(bool left_col, bool bottom_row) {
    TH1F* frame = gPad->DrawFrame(0.0, -1.0, 360.0, 1.0, "");
    frame->GetXaxis()->SetTitle(bottom_row ? "#phi (deg)" : "");
    frame->GetYaxis()->SetTitle(left_col ? "A_{LU}" : "");
    frame->GetXaxis()->SetTitleSize(bottom_row ? 0.075 : 0.0);
    frame->GetYaxis()->SetTitleSize(left_col ? 0.080 : 0.0);
    frame->GetXaxis()->SetLabelSize(bottom_row ? 0.060 : 0.0);
    frame->GetYaxis()->SetLabelSize(left_col ? 0.058 : 0.0);
    frame->GetXaxis()->SetTitleOffset(0.95);
    frame->GetYaxis()->SetTitleOffset(0.92);
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(505);
}

static void make_bsa_plots(const std::string& output_root,
                           const std::vector<RowBin>& rows,
                           const std::map<std::string, std::vector<AsymResult>>& results) {
    const std::filesystem::path base = std::filesystem::path(output_root) / "bsa_plots";
    std::filesystem::create_directories(base);

    std::set<BinRange> xB_set;
    std::set<BinRange> q2_set;
    std::set<BinRange> t_set;
    for (const RowBin& rb : rows) {
        if (!rb.valid) continue;
        xB_set.insert(BinRange{rb.xBmin, rb.xBmax});
        q2_set.insert(BinRange{rb.Q2min, rb.Q2max});
        t_set.insert(BinRange{rb.tmin, rb.tmax});
    } //endfor

    std::vector<BinRange> xB_bins(xB_set.begin(), xB_set.end());
    std::vector<BinRange> q2_bins(q2_set.begin(), q2_set.end());
    std::vector<BinRange> t_bins(t_set.begin(), t_set.end());
    if (xB_bins.empty() || q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[bsa] WARNING: no valid bin ranges found; skipping BSA plots.\n";
        return;
    } //endif

    for (const auto& group_pair : results) {
        const std::string& group = group_pair.first;
        const std::vector<AsymResult>& vec = group_pair.second;
        const std::filesystem::path out_dir = base / sanitize_token(group);
        std::filesystem::create_directories(out_dir);

        int canvas_count = 0;
        for (const BinRange& xbr : xB_bins) {
            const int n_rows = static_cast<int>(q2_bins.size());
            const int n_cols = static_cast<int>(t_bins.size());

            const int canvas_w = std::max(1400, 360 * n_cols);
            const int canvas_h = std::max(900, 285 * n_rows);
            TCanvas c("c_bsa_matrix", "", canvas_w, canvas_h);
            c.SetTopMargin(0.02);
            c.SetBottomMargin(0.02);
            c.SetLeftMargin(0.02);
            c.SetRightMargin(0.02);
            c.Divide(n_cols, n_rows, 0.0005, 0.0005);

            // ROOT pads keep pointers to drawn objects.  Graphs/functions created
            // on the stack are destroyed at the end of each pad iteration before
            // SaveAs(), which leaves only the frame/text visible.  Keep all drawn
            // physics objects alive until after the canvas is written.
            std::vector<std::unique_ptr<TGraphErrors>> owned_graphs;
            std::vector<std::unique_ptr<TF1>> owned_fits;

            int n_populated_pads = 0;
            for (int iq = 0; iq < n_rows; ++iq) {
                for (int it = 0; it < n_cols; ++it) {
                    const int pad_id = iq * n_cols + it + 1;
                    TVirtualPad* pad = c.cd(pad_id);
                    if (!pad) continue;

                    const bool left_col = (it == 0);
                    const bool bottom_row = (iq == n_rows - 1);
                    pad->SetTickx(1);
                    pad->SetTicky(1);
                    pad->SetLeftMargin(left_col ? 0.18 : 0.08);
                    pad->SetRightMargin(0.04);
                    pad->SetBottomMargin(bottom_row ? 0.17 : 0.08);
                    pad->SetTopMargin(0.08);

                    std::vector<double> x_sub;
                    std::vector<double> y_sub;
                    std::vector<double> ex_sub;
                    std::vector<double> ey_sub;

                    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
                        if (!rows[r].valid || r >= static_cast<int>(vec.size())) continue;
                        if (!vec[r].valid) continue;
                        if (!same_range(rows[r].xBmin, rows[r].xBmax, xbr)) continue;
                        if (!same_range(rows[r].Q2min, rows[r].Q2max, q2_bins[iq])) continue;
                        if (!same_range(rows[r].tmin, rows[r].tmax, t_bins[it])) continue;

                        x_sub.push_back(phi_center_deg(rows[r].pmin, rows[r].pmax));
                        y_sub.push_back(vec[r].value);
                        ex_sub.push_back(phi_half_width_deg(rows[r].pmin, rows[r].pmax));
                        ey_sub.push_back(vec[r].stat);
                    } //endfor

                    if (x_sub.empty()) {
                        draw_empty_frame(left_col, bottom_row);
                        TLatex empty_lat;
                        empty_lat.SetNDC();
                        empty_lat.SetTextFont(42);
                        empty_lat.SetTextSize(0.070);
                        empty_lat.DrawLatex(0.20, 0.82, Form("Q^{2}: %.3g-%.3g", q2_bins[iq].lo, q2_bins[iq].hi));
                        empty_lat.DrawLatex(0.20, 0.72, Form("|t|: %.3g-%.3g", t_bins[it].lo, t_bins[it].hi));
                        empty_lat.DrawLatex(0.20, 0.58, "no valid BSA points");
                        continue;
                    } //endif

                    ++n_populated_pads;

                    auto gr = std::make_unique<TGraphErrors>(
                        static_cast<int>(x_sub.size()),
                        x_sub.data(), y_sub.data(), ex_sub.data(), ey_sub.data());
                    gr->SetName(Form("gr_bsa_%d_%d", iq, it));
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerSize(0.95);
                    gr->SetLineWidth(1);

                    // Draw an explicit fixed frame first.  The graph is kept alive
                    // in owned_graphs until after SaveAs(), otherwise ROOT stores a
                    // dangling pointer and the saved PNG contains only the frame/text.
                    draw_empty_frame(left_col, bottom_row);
                    gr->Draw("P SAME");

                    TF1* fit_for_label = nullptr;
                    if (x_sub.size() >= 3) {
                        auto fit = std::make_unique<TF1>(
                            Form("fit_bsa_%d_%d", iq, it),
                            "[0]*sin(x*TMath::Pi()/180.0)/(1.0 + [1]*cos(x*TMath::Pi()/180.0))",
                            0.0, 360.0);
                        fit->SetParameters(0.0, 0.0);
                        fit->SetParNames("A", "B");
                        fit->SetParLimits(1, -0.95, 0.95);
                        fit->SetLineWidth(2);
                        gr->Fit(fit.get(), "Q0");
                        fit->Draw("SAME");
                        fit_for_label = fit.get();
                        owned_fits.push_back(std::move(fit));
                    } //endif

                    TLatex lat;
                    lat.SetNDC();
                    lat.SetTextFont(42);
                    lat.SetTextSize(0.058);
                    lat.DrawLatex(0.20, 0.84, Form("Q^{2}: %.3g-%.3g", q2_bins[iq].lo, q2_bins[iq].hi));
                    lat.DrawLatex(0.20, 0.75, Form("|t|: %.3g-%.3g", t_bins[it].lo, t_bins[it].hi));
                    if (fit_for_label) {
                        lat.DrawLatex(0.20, 0.66, Form("A=%.3f#pm%.3f", fit_for_label->GetParameter(0), fit_for_label->GetParError(0)));
                        lat.DrawLatex(0.20, 0.57, Form("B=%.3f#pm%.3f", fit_for_label->GetParameter(1), fit_for_label->GetParError(1)));
                    } //endif

                    owned_graphs.push_back(std::move(gr));
                } //endfor
            } //endfor

            if (n_populated_pads == 0) continue;

            // Do not draw a canvas-level title here.  The matrix pads occupy the full
            // canvas after TCanvas::Divide, so a global TLatex title overlaps the first
            // row of physics panels.  The group and xB range are encoded in the output
            // directory and filename, while each pad labels its Q2 and |t| bin.

            const std::string name =
                (out_dir / Form("bsa_matrix_%s_%s.png",
                                sanitize_token(group).c_str(),
                                range_token("xB", xbr).c_str())).string();
            c.SaveAs(name.c_str());
            ++canvas_count;
        } //endfor

        std::cout << "[bsa] Wrote " << canvas_count
                  << " xB-matrix BSA canvases for group " << group
                  << " to " << out_dir.string() << "\n";
    } //endfor
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

        const PeriodCounts gamma_counts =
            accumulate_counts(dvcsDataTrees, "DVCS", "DVCS", rows, fast_bins, sigma_cuts);
        const PeriodCounts pi0_counts =
            accumulate_counts(eppi0DataTrees, "eppi0", "eppi0", rows, fast_bins, sigma_cuts);

        std::map<std::string, std::vector<AsymResult>> results;
        for (const std::string& group : output_group_order()) {
            const std::vector<std::string> components = component_periods_for_group(group);
            std::vector<AsymResult> group_results(rows.size());

            for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
                if (!rows[r].valid) continue;
                group_results[r] =
                    compute_group_bsa(gamma_counts, pi0_counts, leaks, components, r, options);
            } //endfor

            results[group] = std::move(group_results);
        } //endfor

        for (const std::string& group : output_group_order()) {
            const std::string col_name = "BSA, counts, " + group;
            const int c = col_strict(csv, col_name);
            const auto& vec = results.at(group);

            for (int r = 0; r < static_cast<int>(csv.rows.size()); ++r) {
                if (r >= static_cast<int>(vec.size()) || !vec[r].valid) {
                    csv.rows[r][c].clear();
                    continue;
                } //endif
                csv.rows[r][c] = fmt_tuple(vec[r].value, vec[r].stat);
            } //endfor
        } //endfor

        write_csv_atomic(options.csv_path, csv);
        std::cout << "[bsa] Updated BSA count columns in " << options.csv_path << "\n";

        const std::filesystem::path json_path =
            std::filesystem::path(options.output_root) / "jsons" / "BSA_counts" / "bsa_counts_summary.json";
        write_json_summary(json_path.string(), rows, results);
        std::cout << "[bsa] Wrote JSON summary to " << json_path.string() << "\n";

        if (options.make_plots) {
            make_bsa_plots(options.output_root, rows, results);
        } //endif

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[bsa] ERROR: " << e.what() << "\n";
        return false;
    } //endtry
}
