// bsa.cpp
// -----------------------------------------------------------------------------
// Direct count-based beam-spin asymmetry stage for the DVCS pass-2 workflow.
//
// This replaces the old BSA module that used stale pi0-corrected-count JSONs and
// fit machinery. The updated logic is deliberately direct:
//
//   1. Loop over the measured DVCS data trees.
//   2. Apply the same analysis-wide global cuts used by the rest of the workflow.
//   3. Apply the same period/topology-dependent DVCS 3-sigma exclusivity cuts
//      from output/jsons/combined_cuts.json.
//   4. Bin surviving events into the existing pass-2 CSV rows using xB, Q2, |t|
//      and phi.
//   5. Count helicity-positive and helicity-negative events.
//   6. Write A_LU = (N+ - N-) / (P_b * (N+ + N-)) into the
//      "BSA, counts, ..." columns.
//
// For combined groups with different beam polarizations, the estimator is not
// formed by dividing by an arbitrary average polarization. Instead, for each CSV
// row:
//
//   A_LU = sum_i (N_i+ - N_i-) / sum_i [P_i * (N_i+ + N_i-)]
//
// where i runs over the component periods in the group. This reduces exactly to
// the usual single-period expression when all component periods have the same
// beam polarization.
// -----------------------------------------------------------------------------

#include "bsa.h"

#include "global_cuts.h"

// ROOT
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TAxis.h>

// JSON
#include <nlohmann/json.hpp>

// C++ stdlib
#include <algorithm>
#include <array>
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

static constexpr double PI      = 3.14159265358979323846;
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
        if (std::isalnum(static_cast<unsigned char>(c))) {
            continue;
        }
        c = '_';
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

static inline bool in_range(double v, double a, double b) {
    return (v >= a) && (v < b);
}

static inline bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (pmax_deg > pmin_deg) {
        return in_range(phi_deg, pmin_deg, pmax_deg);
    }
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg);
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

static inline std::string topo_dir(int det1, int det2) {
    if (det1 == 1 && det2 == 1) {
        return "FD_FD";
    }
    if (det1 == 2 && det2 == 1) {
        return "CD_FD";
    }
    if (det1 == 2 && det2 == 0) {
        return "CD_FT";
    }
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
        "Fa18 Inb",
        "Fa18 Out",
        "Sp19 Inb",
        "Sp18 Inb",
        "Sp18 Out"
    };
    return v;
}

static const std::vector<std::string>& output_group_order() {
    static const std::vector<std::string> v = {
        "Fa18 Inb",
        "Fa18 Out",
        "Sp19 Inb",
        "Sp18 Inb",
        "Sp18 Out",
        "Fa18",
        "Sp18",
        "10.6 GeV"
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
    cur.reserve(line.size());

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

static bool load_csv(const std::string& path, CSV& csv) {
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

    return true;
}

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";
    std::ofstream fout(tmp);
    if (!fout.is_open()) {
        fatal("[bsa] cannot write temporary CSV: " + tmp);
    }

    auto write_cell = [&](const std::string& s) {
        const bool quote =
            s.find(',')  != std::string::npos ||
            s.find('"') != std::string::npos ||
            s.find('\n') != std::string::npos ||
            s.find('\r') != std::string::npos;

        if (!quote) {
            fout << s;
            return;
        }

        fout << '"';
        for (char c : s) {
            if (c == '"') fout << "\"\"";
            else fout << c;
        }
        fout << '"';
    };

    for (size_t i = 0; i < csv.header.size(); ++i) {
        write_cell(csv.header[i]);
        if (i + 1 < csv.header.size()) fout << ',';
    }
    fout << '\n';

    for (const auto& row : csv.rows) {
        if (row.size() != csv.header.size()) {
            fatal("[bsa] CSV row width mismatch while writing");
        }
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

static inline double to_double_strict(const std::string& s, const std::string& what) {
    if (s.empty()) {
        fatal("[bsa] empty numeric cell for " + what);
    }
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str()) {
        fatal("[bsa] parse failure for " + what + " value " + s);
    }
    return v;
}

static inline bool to_bool_valid(const std::string& s) {
    return (s == "1" || s == "1.0" || s == "true" || s == "TRUE");
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
    }
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
    if (it == bins.end()) {
        bins.push_back(b);
    }
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

    for (const auto& r : rows) {
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
        if (!block.is_object() || !block.contains("DVCS")) {
            continue;
        }

        const auto& data = block["DVCS"];
        if (!data.is_object()) {
            continue;
        }

        CutVarMap vm;
        for (auto vit = data.begin(); vit != data.end(); ++vit) {
            const std::string var = vit.key();
            const auto& stats = vit.value();
            if (!stats.is_object() || !stats.contains("mean") || !stats.contains("std")) {
                continue;
            }

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
            }

            if (!std::isfinite(s.cut_low) || !std::isfinite(s.cut_high) || s.cut_high <= s.cut_low) {
                if (std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0) {
                    s.cut_low  = s.mean - 3.0 * s.std;
                    s.cut_high = s.mean + 3.0 * s.std;
                }
            }

            if (std::isfinite(s.cut_high)) {
                vm.emplace(var, s);
            }
        }

        if (!vm.empty()) {
            out.emplace(key, std::move(vm));
        }
    }

    std::cout << "[bsa] Loaded " << out.size()
              << " DVCS topology/period cut blocks from "
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
            "x", "Q2", "t1", "phi2", "Delta_phi",
            "open_angle_ep2", "pTmiss", "Emiss2", "Mx2", "Mx2_1", "Mx2_2",
            "xF", "theta_gamma_gamma",
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
        bind_double("Delta_phi", &Delta_phi, has_Delta_phi);
        bind_double("open_angle_ep2", &open_angle_ep2, has_open_angle);
        bind_double("pTmiss", &pTmiss, has_pTmiss);
        bind_double("Emiss2", &Emiss2, has_Emiss2);
        bind_double("Mx2", &Mx2, has_Mx2);
        bind_double("Mx2_1", &Mx2_1, has_Mx2_1);
        bind_double("Mx2_2", &Mx2_2, has_Mx2_2);
        bind_double("xF", &xF, has_xF);
        bind_double("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
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

    has_val = false;
    return 0.0;
}

static bool passes_one_sigma_cut(const TopoCutMap& cuts,
                                 const std::string& key,
                                 const BranchBinder& b,
                                 const std::string& var) {
    auto it = cuts.find(key);
    if (it == cuts.end()) {
        fatal("[bsa] missing 3-sigma cut key in combined_cuts.json: " + key);
    }

    const CutVarMap& vm = it->second;
    auto iv = vm.find(var);
    if (iv == vm.end()) {
        return true;
    }

    bool has_val = false;
    const double val = branch_value_for_sigma_var(b, var, has_val);
    if (!has_val) {
        fatal("[bsa] cut key " + key + " requires missing branch: " + var);
    }

    return within_cut_window(val, iv->second);
}

static bool passes_sigma_cuts(const TopoCutMap& cuts,
                              const std::string& key,
                              const BranchBinder& b) {
    static const std::vector<std::string> vars = {
        "Emiss2",
        "Mx2",
        "Mx2_1",
        "Mx2_2",
        "Delta_phi",
        "pTmiss",
        "xF",
        "theta_gamma_gamma"
    };

    for (const std::string& var : vars) {
        if (!passes_one_sigma_cut(cuts, key, b, var)) {
            return false;
        }
    }
    return true;
}

static bool passes_global_cuts_dispatch(const BranchBinder& b,
                                        const std::string& period_label) {
    const GlobalCutConfig& cfg = default_global_cuts();

    if (!(b.has_t1 && b.has_open_angle)) return false;
    if (cfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (b.has_runnum && is_excluded_run(b.runnum)) return false;

    if (cfg.enable_topology_filter || global_cuts_require_sector_phi(cfg) || cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_detector1 && b.has_detector2)) {
            fatal("[bsa] topology/sector/ycol selection requires detector1 and detector2 branches");
        }
    }

    if (global_cuts_require_sector_phi(cfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            fatal("[bsa] sector selection requires e_phi, p1_phi and p2_phi branches");
        }
    }

    if (global_cuts_require_auxiliary_kinematics(cfg)) {
        if (!(b.has_e_theta && b.has_e_phi &&
              b.has_p1_theta && b.has_p1_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[bsa] auxiliary fiducial cuts require e_theta, e_phi, p1_theta, p1_phi, p2_p, p2_theta and p2_phi branches");
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
            fatal("[bsa] dvcsgen ycol cut requires e_p, e_theta, e_phi, p2_p, p2_theta and p2_phi branches");
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

struct AsymResult {
    bool valid = false;
    double value = 0.0;
    double stat = 0.0;
    double n_plus = 0.0;
    double n_minus = 0.0;
    double denominator = 0.0;
};

static inline void add_event(HelCounts& h, int helicity) {
    if (helicity > 0) {
        h.plus += 1.0;
    } else if (helicity < 0) {
        h.minus += 1.0;
    }
}

static PeriodCounts accumulate_counts(const std::map<std::string, TTree*>& trees,
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
            std::cout << "[bsa] Skipping non-canonical or supplemental data tree key: "
                      << tree_key << "\n";
            continue;
        }

        BranchBinder b;
        b.bind(tree);

        if (!(b.has_detector1 && b.has_detector2 && b.has_helicity &&
              b.has_x && b.has_Q2 && b.has_t1 && b.has_phi2)) {
            fatal("[bsa] tree " + tree_key + " is missing one or more required branches: detector1, detector2, helicity, x, Q2, t1, phi2");
        }

        const std::string period_code = period_code_from_display(period);
        const Long64_t n_entries = tree->GetEntries();

        long long n_topology = 0;
        long long n_global = 0;
        long long n_sigma = 0;
        long long n_matched = 0;

        for (Long64_t i = 0; i < n_entries; ++i) {
            tree->GetEntry(i);

            const std::string topo = topo_dir(b.detector1, b.detector2);
            if (topo.empty()) {
                continue;
            }
            ++n_topology;

            if (!passes_global_cuts_dispatch(b, period)) {
                continue;
            }
            ++n_global;

            const std::string cut_key = "DVCS_" + period_code + "_" + topo;
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

        std::cout << "[bsa] tree=" << tree_key
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

static AsymResult compute_group_bsa(const PeriodCounts& counts,
                                    const std::vector<std::string>& component_periods,
                                    int row_index,
                                    const BSAOptions& opt) {
    AsymResult r;

    double numerator = 0.0;
    double denominator = 0.0;
    double variance_numerator = 0.0;
    double n_plus_total = 0.0;
    double n_minus_total = 0.0;

    for (const std::string& period : component_periods) {
        auto ip = counts.find(period);
        if (ip == counts.end()) {
            continue;
        }
        auto ir = ip->second.find(row_index);
        if (ir == ip->second.end()) {
            continue;
        }

        const double np = ir->second.plus;
        const double nm = ir->second.minus;
        const double s = np + nm;
        const double d = np - nm;
        if (!(s > 0.0)) {
            continue;
        }

        const double P = beam_pol_for_period(period, opt);
        numerator += d;
        denominator += P * s;
        n_plus_total += np;
        n_minus_total += nm;

        // Conditional-binomial variance for D = N+ - N- at fixed S. This is
        // equivalent to 4*N+*N-/S and gives the usual asymmetry uncertainty.
        variance_numerator += std::max(0.0, s - (d * d / s));
    }

    r.n_plus = n_plus_total;
    r.n_minus = n_minus_total;
    r.denominator = denominator;

    if (!(denominator > 0.0)) {
        return r;
    }

    r.value = numerator / denominator;
    r.stat = std::sqrt(std::max(0.0, variance_numerator)) / denominator;
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
// JSON and plots
// -----------------------------------------------------------------------------

static void write_json_summary(const std::string& path,
                               const std::vector<RowBin>& rows,
                               const std::map<std::string, std::vector<AsymResult>>& results) {
    nlohmann::json j;
    j["description"] = "Direct count-based DVCS beam-spin asymmetries after global cuts and DVCS 3-sigma exclusivity cuts.";
    j["estimator"] = "A_LU = sum_i(Nplus_i - Nminus_i) / sum_i[Pbeam_i*(Nplus_i + Nminus_i)]";

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
            row["Nplus"] = a.n_plus;
            row["Nminus"] = a.n_minus;
            row["polarized_denominator"] = a.denominator;
            row["valid"] = a.valid;
            if (a.valid) {
                row["BSA"] = a.value;
                row["stat"] = a.stat;
            }
            rows_json.push_back(row);
        }

        j["groups"][group] = rows_json;
    }

    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream out(path);
    if (!out.is_open()) {
        fatal("[bsa] cannot write JSON summary: " + path);
    }
    out << std::setw(2) << j << "\n";
}

struct CellKey {
    double xBmin = 0.0;
    double xBmax = 0.0;
    double Q2min = 0.0;
    double Q2max = 0.0;
    double tmin = 0.0;
    double tmax = 0.0;

    bool operator<(const CellKey& other) const {
        return std::tie(xBmin, xBmax, Q2min, Q2max, tmin, tmax) <
               std::tie(other.xBmin, other.xBmax, other.Q2min, other.Q2max, other.tmin, other.tmax);
    }
};

struct PlotPoint {
    double phi = 0.0;
    double phi_err = 0.0;
    double bsa = 0.0;
    double bsa_err = 0.0;
};

static void make_bsa_plots(const std::string& output_root,
                           const std::vector<RowBin>& rows,
                           const std::map<std::string, std::vector<AsymResult>>& results) {
    const std::filesystem::path base = std::filesystem::path(output_root) / "bsa_plots";
    std::filesystem::create_directories(base);

    for (const auto& group_pair : results) {
        const std::string& group = group_pair.first;
        const std::vector<AsymResult>& vec = group_pair.second;
        const std::filesystem::path out_dir = base / sanitize_token(group);
        std::filesystem::create_directories(out_dir);

        std::map<CellKey, std::vector<PlotPoint>> cells;
        for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
            if (!vec[r].valid) {
                continue;
            }
            const RowBin& rb = rows[r];
            CellKey key{rb.xBmin, rb.xBmax, rb.Q2min, rb.Q2max, rb.tmin, rb.tmax};
            PlotPoint p;
            p.phi = phi_center_deg(rb.pmin, rb.pmax);
            p.phi_err = phi_half_width_deg(rb.pmin, rb.pmax);
            p.bsa = vec[r].value;
            p.bsa_err = vec[r].stat;
            cells[key].push_back(p);
        }

        int canvas_index = 0;
        for (auto& cell_pair : cells) {
            CellKey key = cell_pair.first;
            std::vector<PlotPoint>& pts = cell_pair.second;
            if (pts.empty()) {
                continue;
            }

            std::sort(pts.begin(), pts.end(), [](const PlotPoint& a, const PlotPoint& b) {
                return a.phi < b.phi;
            });

            std::vector<double> x, y, ex, ey;
            x.reserve(pts.size());
            y.reserve(pts.size());
            ex.reserve(pts.size());
            ey.reserve(pts.size());
            for (const PlotPoint& p : pts) {
                x.push_back(p.phi);
                y.push_back(p.bsa);
                ex.push_back(0.0);
                ey.push_back(p.bsa_err);
            }

            TCanvas c(Form("c_bsa_%s_%d", sanitize_token(group).c_str(), canvas_index),
                      "BSA", 1100, 800);
            c.SetLeftMargin(0.12);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);
            c.SetTopMargin(0.08);

            TGraphErrors gr(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
            gr.SetMarkerStyle(20);
            gr.SetMarkerSize(1.0);
            gr.SetLineWidth(2);
            gr.GetXaxis()->SetTitle("#phi (deg)");
            gr.GetYaxis()->SetTitle("A_{LU}");
            gr.GetXaxis()->SetLimits(0.0, 360.0);
            gr.GetYaxis()->SetRangeUser(-1.0, 1.0);
            gr.Draw("AP");

            TF1 fit("fit_sin", "[0]*sin(x*TMath::Pi()/180.0)", 0.0, 360.0);
            fit.SetParameter(0, 0.0);
            if (x.size() >= 3) {
                gr.Fit(&fit, "Q0");
                fit.SetLineWidth(2);
                fit.Draw("same");
            }

            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextFont(42);
            lat.SetTextSize(0.034);
            lat.DrawLatex(0.16, 0.92, Form("%s", group.c_str()));
            lat.DrawLatex(0.16, 0.875,
                          Form("%.3g < x_{B} < %.3g, %.3g < Q^{2} < %.3g GeV^{2}",
                               key.xBmin, key.xBmax, key.Q2min, key.Q2max));
            lat.DrawLatex(0.16, 0.83,
                          Form("%.3g < |t| < %.3g GeV^{2}", key.tmin, key.tmax));
            if (x.size() >= 3) {
                lat.DrawLatex(0.16, 0.785,
                              Form("sin#phi amplitude = %.4f #pm %.4f", fit.GetParameter(0), fit.GetParError(0)));
            }

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
                           const BSAOptions& options) {
    try {
        CSV csv;
        load_csv(options.csv_path, csv);

        const std::vector<RowBin> rows = load_row_bins_from_csv(csv);
        const FastBinning fast_bins = build_fast_binning(rows);
        const TopoCutMap sigma_cuts = load_combined_cuts(options.combined_cuts_json);

        PeriodCounts counts = accumulate_counts(dvcsDataTrees, rows, fast_bins, sigma_cuts);

        std::map<std::string, std::vector<AsymResult>> results;
        for (const std::string& group : output_group_order()) {
            const std::vector<std::string> components = component_periods_for_group(group);
            std::vector<AsymResult> group_results(rows.size());

            for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
                if (!rows[r].valid) {
                    continue;
                }
                group_results[r] = compute_group_bsa(counts, components, r, options);
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
