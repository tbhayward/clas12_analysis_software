// yield_totals.cpp
// -----------------------------------------------------------------------------
// Current-binned totals for the final normalized DATA yields.
//
// This module intentionally remains an event-loop diagnostic because the CSV
// contains period/bin totals but not run/current-resolved totals. The event loop
// is used only to recover the current split. The weights and corrections come
// from the CSV products already produced upstream.
//
// Outputs by period and current:
//
//   1) DVCS normalized pi0-subtracted counts
//      event weight = 1 / [current_eff(epg, exp, period) * R_pi0(region, p1_theta)]
//      row weight   = event weight * [1 - contamination_ratio(row, period)]
//
//   2) eppi0 normalized counts
//      event weight = 1 / [current_eff(eppi0, exp, period) * R_pi0(region, p1_theta)]
//
// Cuts:
//   - Global cuts via global_cuts.h.
//   - Topology from detector1/detector2.
//   - Channel-specific 3-sigma cuts from combined_cuts.json:
//        DVCS_<PeriodDir>_<TopoDir>   for ep->epgamma
//        eppi0_<PeriodDir>_<TopoDir>  for ep->eppi0
//   - Row matching in xB, Q2, |t|, phi.
// -----------------------------------------------------------------------------

#include "yield_totals.h"

#include "global_cuts.h"

#include <nlohmann/json.hpp>

#include <TTree.h>
#include <TROOT.h>

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

using json = nlohmann::json;

static constexpr double PI = 3.14159265358979323846;
static constexpr double RAD2DEG = 180.0 / PI;

[[noreturn]] static void fatal(const std::string& msg) {
    throw std::runtime_error("[yield_totals] FATAL: " + msg);
}


static void ensure_parent_directory_exists(const std::string& path_string) {
    const std::filesystem::path path(path_string);
    const std::filesystem::path parent = path.parent_path();

    if (!parent.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
        if (ec) {
            fatal("cannot create output directory: " + parent.string() + " (" + ec.message() + ")");
        }
    }
}


static std::string lower_ascii(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

static bool has_substr(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

static double wrap_phi_deg(double phi_deg) {
    double p = std::fmod(phi_deg, 360.0);
    if (p < 0.0) {
        p += 360.0;
    }
    if (p >= 360.0) {
        p = std::nextafter(360.0, 0.0);
    }
    return p;
}


static double delta_phi_rad_from_two_phi(double phi_a, double phi_b) {
    double d = std::fmod(phi_a - phi_b, 2.0 * PI);

    if (d <= -PI) {
        d += 2.0 * PI;
    }

    if (d > PI) {
        d -= 2.0 * PI;
    }

    return std::fabs(d);
}

static bool in_range(double v, double lo, double hi) {
    return v >= lo && v < hi;
}

static bool row_accepts_phi(double phi_deg, double pmin, double pmax) {
    if (pmax > pmin) {
        return in_range(phi_deg, pmin, pmax);
    }
    return phi_deg >= pmin || phi_deg < pmax;
}

struct PeriodTags {
    std::string display;
    std::string period_label;
    std::string period_code;
    std::string internal;
};

static PeriodTags parse_period_from_key(const std::string& key) {
    const std::string s = lower_ascii(key);
    PeriodTags t;

    if (has_substr(s, "fa18") && has_substr(s, "inb") && has_substr(s, "supp")) {
        t.display = "Fa18 Inb Supp";
        t.period_label = "fa18_inb";
        t.period_code = "Fa18_Inb_Supp";
        t.internal = "rga_fa18_inb";
        return t;
    }

    if (has_substr(s, "fa18") && has_substr(s, "inb")) {
        t.display = "Fa18 Inb";
        t.period_label = "fa18_inb";
        t.period_code = "Fa18_Inb";
        t.internal = "rga_fa18_inb";
        return t;
    }

    if (has_substr(s, "fa18") && has_substr(s, "out")) {
        t.display = "Fa18 Out";
        t.period_label = "fa18_out";
        t.period_code = "Fa18_Out";
        t.internal = "rga_fa18_out";
        return t;
    }

    if (has_substr(s, "sp19") && has_substr(s, "inb")) {
        t.display = "Sp19 Inb";
        t.period_label = "sp19_inb";
        t.period_code = "Sp19_Inb";
        t.internal = "rga_sp19_inb";
        return t;
    }

    if (has_substr(s, "sp18") && has_substr(s, "inb")) {
        t.display = "Sp18 Inb";
        t.period_label = "sp18_inb";
        t.period_code = "Sp18_Inb";
        t.internal = "rga_sp18_inb";
        return t;
    }

    if (has_substr(s, "sp18") && has_substr(s, "out")) {
        t.display = "Sp18 Out";
        t.period_label = "sp18_out";
        t.period_code = "Sp18_Out";
        t.internal = "rga_sp18_out";
        return t;
    }

    fatal("cannot parse period from tree key: " + key);
    return t;
}

static const std::vector<std::string>& period_order() {
    static const std::vector<std::string> v = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    return v;
}

static const std::vector<std::string>& topo_dirs() {
    static const std::vector<std::string> v = {"FD_FD", "CD_FD", "CD_FT"};
    return v;
}

static std::string topo_dir_from_detectors(int det1, int det2) {
    if (det1 == 1 && det2 == 1) return "FD_FD";
    if (det1 == 2 && det2 == 1) return "CD_FD";
    if (det1 == 2 && det2 == 0) return "CD_FT";
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

static void load_csv(const std::string& path, CSV& csv) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        fatal("cannot open CSV: " + path);
    }

    std::string line;
    if (!std::getline(fin, line)) {
        fatal("empty CSV: " + path);
    }

    csv.header = split_csv_line(line);
    csv.index.clear();
    for (int i = 0; i < static_cast<int>(csv.header.size()); ++i) {
        if (csv.index.find(csv.header[i]) != csv.index.end()) {
            fatal("duplicate CSV column: " + csv.header[i]);
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
            fatal("CSV row width mismatch in " + path);
        }
        csv.rows.push_back(std::move(row));
    }
}

static int col(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);
    if (it == csv.index.end()) {
        fatal("missing CSV column: " + name);
    }
    return it->second;
}

static double to_double_strict(const std::string& s, const std::string& what) {
    if (s.empty()) {
        fatal("empty numeric value for " + what);
    }
    char* endp = nullptr;
    const double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) {
        fatal("could not parse numeric value for " + what + ": " + s);
    }
    return v;
}

static bool is_valid_cell(const std::string& s) {
    return s == "1" || s == "1.0" || s == "true" || s == "TRUE";
}

static std::vector<double> parse_tuple_numbers(const std::string& cell,
                                               const std::string& what) {
    std::vector<double> out;
    size_t i = 0;

    while (i < cell.size()) {
        while (i < cell.size() && !(std::isdigit(static_cast<unsigned char>(cell[i])) ||
                                    cell[i] == '-' || cell[i] == '+' || cell[i] == '.')) {
            ++i;
        }
        if (i >= cell.size()) {
            break;
        }

        char* endp = nullptr;
        const double v = std::strtod(cell.c_str() + i, &endp);
        if (endp == cell.c_str() + i) {
            break;
        }
        out.push_back(v);
        i = static_cast<size_t>(endp - cell.c_str());
    }

    if (out.empty()) {
        fatal("could not parse tuple for " + what + ": " + cell);
    }

    return out;
}

// -----------------------------------------------------------------------------
// Row bins and per-period CSV products
// -----------------------------------------------------------------------------

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

struct RegionNormalization {
    std::vector<double> cubic;
};

struct PeriodProducts {
    double epg_current_eff = std::numeric_limits<double>::quiet_NaN();
    double eppi0_current_eff = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, RegionNormalization> regions;
};

struct AnalysisInputs {
    std::vector<RowBin> rows;
    std::map<std::string, PeriodProducts> products_by_period;
    std::map<std::string, std::vector<double>> contamination_by_period;
};

static std::string current_eff_col(const std::string& channel,
                                   const std::string& period) {
    return "current efficiency factor, " + channel + ", exp, " + period;
}

static const std::vector<std::string>& normalization_regions() {
    static const std::vector<std::string> v = {
        "Sector 1", "Sector 2", "Sector 3", "Sector 4", "Sector 5", "Sector 6", "CD"
    };
    return v;
}

static std::string norm_cubic_col(const std::string& period,
                                  const std::string& region) {
    return "eppi0 cross-section normalization cubic, ep->eppi0, data_over_mc, " + region + ", " + period;
}

static std::string contamination_col(const std::string& period) {
    return "contamination ratio, " + period;
}

static AnalysisInputs build_analysis_inputs(const CSV& csv) {
    AnalysisInputs in;

    const int c_xBmin = col(csv, "xBmin");
    const int c_xBmax = col(csv, "xBmax");
    const int c_Q2min = col(csv, "Q2min");
    const int c_Q2max = col(csv, "Q2max");
    const int c_tmin = col(csv, "t_abs_min");
    const int c_tmax = col(csv, "t_abs_max");
    const int c_phimin = col(csv, "phimin");
    const int c_phimax = col(csv, "phimax");
    const int c_valid = col(csv, "valid bin");

    in.rows.reserve(csv.rows.size());

    for (int r = 0; r < static_cast<int>(csv.rows.size()); ++r) {
        const auto& row = csv.rows[r];
        RowBin b;
        b.xBmin = to_double_strict(row[c_xBmin], "xBmin");
        b.xBmax = to_double_strict(row[c_xBmax], "xBmax");
        b.Q2min = to_double_strict(row[c_Q2min], "Q2min");
        b.Q2max = to_double_strict(row[c_Q2max], "Q2max");
        b.tmin = to_double_strict(row[c_tmin], "t_abs_min");
        b.tmax = to_double_strict(row[c_tmax], "t_abs_max");
        b.phimin = to_double_strict(row[c_phimin], "phimin");
        b.phimax = to_double_strict(row[c_phimax], "phimax");
        b.valid = is_valid_cell(row[c_valid]);
        in.rows.push_back(b);
    }

    for (const std::string& period : period_order()) {
        PeriodProducts p;

        const int c_epg_eff = col(csv, current_eff_col("ep->epg", period));
        const int c_epi_eff = col(csv, current_eff_col("ep->eppi0", period));
        const int c_cont = col(csv, contamination_col(period));

        // These are period-level constants copied into every row. Read from the
        // first row and then validate only enough to fail fast if missing.
        if (csv.rows.empty()) {
            fatal("CSV has no data rows.");
        }

        const std::vector<double> epg_eff = parse_tuple_numbers(csv.rows.front()[c_epg_eff],
                                                                current_eff_col("ep->epg", period));
        const std::vector<double> epi_eff = parse_tuple_numbers(csv.rows.front()[c_epi_eff],
                                                                current_eff_col("ep->eppi0", period));

        if (epg_eff.size() < 1 || epi_eff.size() < 1) {
            fatal("malformed current-efficiency columns for period " + period);
        }

        p.epg_current_eff = epg_eff[0];
        p.eppi0_current_eff = epi_eff[0];

        for (const std::string& region : normalization_regions()) {
            const std::string cname = norm_cubic_col(period, region);
            const int c_cubic = col(csv, cname);
            const std::vector<double> cubic = parse_tuple_numbers(csv.rows.front()[c_cubic], cname);
            if (cubic.size() < 4) {
                fatal("malformed regional eppi0 normalization cubic column for period " + period + ", region " + region);
            }
            RegionNormalization rn;
            rn.cubic.assign(cubic.begin(), cubic.begin() + 4);
            p.regions[region] = rn;
        }

        if (!(std::isfinite(p.epg_current_eff) && p.epg_current_eff > 0.0)) {
            fatal("invalid epg current-efficiency factor for period " + period);
        }
        if (!(std::isfinite(p.eppi0_current_eff) && p.eppi0_current_eff > 0.0)) {
            fatal("invalid eppi0 current-efficiency factor for period " + period);
        }

        in.products_by_period[period] = p;

        std::vector<double> contamination;
        contamination.resize(csv.rows.size(), 0.0);

        for (int r = 0; r < static_cast<int>(csv.rows.size()); ++r) {
            const std::string& cell = csv.rows[r][c_cont];
            if (cell.empty()) {
                contamination[r] = 0.0;
                continue;
            }

            const std::vector<double> vals = parse_tuple_numbers(cell, contamination_col(period));
            contamination[r] = vals[0];
            if (!std::isfinite(contamination[r])) {
                contamination[r] = 0.0;
            }
        }

        in.contamination_by_period[period] = std::move(contamination);
    }

    return in;
}

static double eval_cubic(const std::vector<double>& p, double x) {
    if (p.size() < 4) {
        fatal("normalization cubic has fewer than 4 coefficients.");
    }

    return p[0] + x * (p[1] + x * (p[2] + x * p[3]));
}

static int clas12_sector_from_phi_deg(double phi_deg) {
    const double phi = wrap_phi_deg(phi_deg);

    if ((phi >= 330.0 && phi < 360.0) || (phi >= 0.0 && phi < 30.0)) return 1;
    if (phi >= 30.0  && phi < 90.0)  return 2;
    if (phi >= 90.0  && phi < 150.0) return 3;
    if (phi >= 150.0 && phi < 210.0) return 4;
    if (phi >= 210.0 && phi < 270.0) return 5;
    if (phi >= 270.0 && phi < 330.0) return 6;

    return 0;
}

static std::string normalization_region_for_event(double p1_theta_deg,
                                                  double p1_phi_deg) {
    if (p1_theta_deg >= 40.0 && p1_theta_deg < 70.0) {
        return "CD";
    }

    if (p1_theta_deg >= 0.0 && p1_theta_deg < 40.0) {
        const int sector = clas12_sector_from_phi_deg(p1_phi_deg);
        if (sector >= 1 && sector <= 6) {
            return "Sector " + std::to_string(sector);
        }
    }

    return "";
}

static int find_row_index(const std::vector<RowBin>& rows,
                          double x,
                          double Q2,
                          double tabs,
                          double phi_deg) {
    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        const RowBin& w = rows[r];
        if (!w.valid) continue;
        if (!in_range(x, w.xBmin, w.xBmax)) continue;
        if (!in_range(Q2, w.Q2min, w.Q2max)) continue;
        if (!in_range(tabs, w.tmin, w.tmax)) continue;
        if (!row_accepts_phi(phi_deg, w.phimin, w.phimax)) continue;
        return r;
    }

    return -1;
}

// -----------------------------------------------------------------------------
// Current maps
// -----------------------------------------------------------------------------

static const std::unordered_map<int, int>& fa18_inb_current_map() {
    static const std::unordered_map<int, int> m = {
        {5418,5},{5419,5},
        {5335,40},{5339,40},{5340,40},{5341,40},{5342,40},{5343,40},{5344,40},
        {5032,45},{5036,45},{5038,45},{5039,45},{5040,45},{5041,45},{5043,45},{5045,45},
        {5046,45},{5047,45},{5051,45},{5052,45},{5053,45},
        {5116,45},{5117,45},{5119,45},{5120,45},{5124,45},{5125,45},{5126,45},{5127,45},
        {5128,45},{5129,45},{5130,45},{5139,45},{5153,45},{5158,45},{5159,45},{5160,45},
        {5162,45},{5163,45},{5164,45},{5165,45},{5166,45},{5167,45},{5168,45},{5169,45},
        {5180,45},{5181,45},{5182,45},{5183,45},{5190,45},{5191,45},{5193,45},{5195,45},
        {5196,45},{5197,45},{5198,45},{5199,45},{5200,45},{5201,45},{5202,45},{5203,45},
        {5204,45},{5205,45},{5206,45},{5208,45},{5211,45},{5212,45},{5215,45},{5216,45},
        {5219,45},{5220,45},{5221,45},{5222,45},{5223,45},{5230,45},{5231,45},{5232,45},
        {5233,45},{5234,45},{5235,45},{5237,45},{5238,45},{5239,45},{5247,45},{5248,45},
        {5249,45},{5252,45},{5253,45},{5257,45},{5258,45},{5259,45},{5261,45},{5262,45},
        {5303,45},{5304,45},{5305,45},{5306,45},{5307,45},{5310,45},{5311,45},{5315,45},
        {5317,45},{5318,45},{5319,45},{5320,45},{5323,45},{5324,45},{5333,45},{5334,45},
        {5336,45},{5346,45},{5347,45},{5349,45},{5351,45},{5354,45},{5355,45},{5367,45},
        {5356,50},{5357,50},{5358,50},{5359,50},{5360,50},{5361,50},{5362,50},{5366,50},
        {5368,55},{5369,55},{5372,55},{5373,55},{5374,55},{5375,55},{5376,55},{5377,55},
        {5378,55},{5379,55},{5380,55},{5381,55},{5382,55},{5383,55},{5386,55},{5390,55},
        {5391,55},{5392,55},{5393,55},{5398,55},{5400,55},{5401,55},{5403,55},{5404,55},
        {5406,55},{5407,55}
    };
    return m;
}

static const std::unordered_map<int, int>& fa18_out_current_map() {
    static const std::unordered_map<int, int> m = {
        {5443,5},{5610,5},{5444,20},
        {5423,40},{5424,40},{5425,40},{5426,40},{5428,40},{5429,40},{5430,40},
        {5432,40},{5434,40},{5435,40},{5436,40},{5437,40},{5438,40},{5440,40},
        {5441,40},{5442,40},{5445,40},{5447,40},{5448,40},{5449,40},{5450,40},
        {5451,40},{5452,40},{5453,40},{5454,40},{5455,40},{5460,40},{5464,40},
        {5465,40},{5466,40},{5467,40},{5468,40},{5469,40},{5470,40},{5471,40},
        {5472,40},{5473,40},{5474,40},{5475,40},{5476,40},{5478,40},{5479,40},
        {5480,40},{5481,40},{5482,40},{5483,40},{5485,40},{5486,40},{5487,40},
        {5495,40},{5496,40},{5497,40},{5498,40},{5499,40},{5500,40},{5504,40},
        {5505,50},{5507,50},{5516,50},{5517,50},{5518,50},{5519,50},{5520,50},
        {5521,50},{5522,50},{5523,50},{5524,50},{5525,50},{5526,50},{5527,50},
        {5528,50},{5530,50},{5532,50},{5533,50},{5534,50},{5535,50},{5536,50},
        {5537,50},{5538,50},{5540,50},{5541,50},{5543,50},{5544,50},{5545,50},
        {5546,50},{5547,50},{5548,50},{5549,50},{5550,50},{5551,50},{5552,50},
        {5555,50},{5556,50},{5557,50},{5558,50},{5559,50},{5562,50},{5567,50},
        {5569,50},{5570,50},{5571,50},{5572,50},{5573,50},{5574,50},{5577,50},
        {5578,50},{5591,50},{5592,50},{5594,50},{5597,50},{5598,50},{5600,50},
        {5601,50},{5602,50},{5603,50},{5604,50},{5606,50},{5607,50},{5611,50},
        {5612,50},{5613,50},{5614,50},{5615,50},{5616,50},{5617,50},{5618,50},
        {5619,50},{5621,50},{5623,50},{5624,50},{5625,50},{5626,50},{5627,50},
        {5628,50},{5629,50},{5630,50},{5631,50},{5632,50},{5633,50},{5635,50},
        {5637,50},{5638,50},{5639,50},{5641,50},{5643,50},{5644,50},{5645,50},
        {5646,50},{5647,50},{5648,50},{5649,50},{5650,50},{5651,50},{5652,50},
        {5654,50},{5655,50},{5656,50},{5662,50},{5663,50},{5664,50},{5665,50},
        {5666,50}
    };
    return m;
}

static bool resolve_current(const std::string& period_internal,
                            int runnum,
                            int& current_nA) {
    if (period_internal == "rga_fa18_inb") {
        const auto& m = fa18_inb_current_map();
        auto it = m.find(runnum);
        if (it == m.end()) return false;
        current_nA = it->second;
        return true;
    }

    if (period_internal == "rga_fa18_out") {
        const auto& m = fa18_out_current_map();
        auto it = m.find(runnum);
        if (it == m.end()) return false;
        current_nA = it->second;
        return true;
    }

    if (period_internal == "rga_sp18_out") {
        if (runnum >= 3211 && runnum <= 3293) {
            current_nA = 30;
            return true;
        }
        if (runnum >= 3867 && runnum <= 3987) {
            current_nA = 45;
            return true;
        }
        return false;
    }

    if (period_internal == "rga_sp18_inb") {
        if (runnum == 3418) {
            current_nA = 70;
            return true;
        }
        if (runnum == 3421 || runnum == 3422) {
            current_nA = 35;
            return true;
        }
        if (runnum == 3429) {
            current_nA = 50;
            return true;
        }
        if (runnum >= 3306 && runnum <= 3411) {
            current_nA = 35;
            return true;
        }
        if (runnum >= 3431 && runnum <= 4325) {
            current_nA = 50;
            return true;
        }
        return false;
    }

    if (period_internal == "rga_sp19_inb") {
        if (runnum == 6616) {
            current_nA = 5;
            return true;
        }
        current_nA = 50;
        return true;
    }

    return false;
}

// -----------------------------------------------------------------------------
// Sigma cuts
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

static TopoCutMap load_sigma_cuts(const std::string& path,
                                  const std::string& sample_key) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        fatal("cannot open combined cuts JSON: " + path);
    }

    json j;
    fin >> j;

    TopoCutMap out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const json& block = it.value();
        if (!block.is_object() || !block.contains(sample_key) || !block[sample_key].is_object()) {
            continue;
        }

        CutVarMap vm;
        for (auto vit = block[sample_key].begin(); vit != block[sample_key].end(); ++vit) {
            const std::string var = vit.key();
            const json& obj = vit.value();
            if (!obj.is_object() || !obj.contains("mean") || !obj.contains("std")) {
                continue;
            }

            SigmaStats s;
            s.mean = obj["mean"].get<double>();
            s.std = obj["std"].get<double>();
            if (obj.contains("cut_low")) s.cut_low = obj["cut_low"].get<double>();
            if (obj.contains("cut_high")) s.cut_high = obj["cut_high"].get<double>();
            if (obj.contains("quantile")) s.quantile = obj["quantile"].get<double>();
            if (obj.contains("mode")) s.mode = obj["mode"].get<std::string>();
            if (!std::isfinite(s.cut_low) || !std::isfinite(s.cut_high) || s.cut_high <= s.cut_low) {
                if (std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0) {
                    s.cut_low = s.mean - 3.0 * s.std;
                    s.cut_high = s.mean + 3.0 * s.std;
                }
            }
            if (std::isfinite(s.cut_high)) {
                vm[var] = s;
            }
        }

        if (!vm.empty()) {
            out[key] = vm;
        }
    }

    return out;
}

static bool within_cut_window(double value, const SigmaStats& s) {
    if (!std::isfinite(value)) return false;
    if (s.mode == "upper_quantile") {
        if (!std::isfinite(s.cut_high)) return true;
        return value <= s.cut_high;
    }
    double lo = s.cut_low;
    double hi = s.cut_high;
    if (!(std::isfinite(lo) && std::isfinite(hi)) || hi <= lo) {
        if (!std::isfinite(s.mean) || !std::isfinite(s.std) || s.std <= 0.0) return true;
        lo = s.mean - 3.0 * s.std;
        hi = s.mean + 3.0 * s.std;
    }
    return (value >= lo && value <= hi);
}

static bool check_sigma(const CutVarMap& vm,
                        const std::string& var,
                        bool has_value,
                        double value) {
    auto it = vm.find(var);
    if (it == vm.end()) {
        return true;
    }
    if (!has_value) {
        fatal("required sigma-cut branch missing: " + var);
    }
    return within_cut_window(value, it->second);
}

// -----------------------------------------------------------------------------
// Branches/cuts
// -----------------------------------------------------------------------------

static std::mutex g_bind_mutex;

struct Branches {
    int runnum = 0; bool has_runnum = false;
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    double x = 0.0; bool has_x = false;
    double Q2 = 0.0; bool has_Q2 = false;
    double t1 = 0.0; bool has_t1 = false;
    double phi2 = 0.0; bool has_phi2 = false;

    double open_angle_ep2 = 0.0; bool has_open_angle_ep2 = false;
    double pTmiss = 0.0; bool has_pTmiss = false;
    double Delta_phi = 0.0; bool has_Delta_phi = false;

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
        if (!t) {
            fatal("attempted to bind null tree.");
        }

        std::lock_guard<std::mutex> lock(g_bind_mutex);
        t->SetBranchStatus("*", 0);

        auto ena = [&](const char* n) {
            if (t->GetBranch(n)) {
                t->SetBranchStatus(n, 1);
            }
        };

        ena("runnum");
        ena("detector1");
        ena("detector2");
        ena("x");
        ena("Q2");
        ena("t1");
        ena("phi2");
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
        ena("p1_theta");
        ena("p1_phi");
        ena("p2_p");
        ena("p2_theta");
        ena("p2_phi");

        t->SetCacheSize(0);

        auto bI = [&](const char* n, int* ptr, bool& flag) {
            if (t->GetBranch(n)) {
                t->SetBranchAddress(n, ptr);
                flag = true;
            }
        };

        auto bD = [&](const char* n, double* ptr, bool& flag) {
            if (t->GetBranch(n)) {
                t->SetBranchAddress(n, ptr);
                flag = true;
            }
        };

        bI("runnum", &runnum, has_runnum);
        bI("detector1", &detector1, has_detector1);
        bI("detector2", &detector2, has_detector2);

        bD("x", &x, has_x);
        bD("Q2", &Q2, has_Q2);
        bD("t1", &t1, has_t1);
        bD("phi2", &phi2, has_phi2);

        bD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bD("pTmiss", &pTmiss, has_pTmiss);
        bD("Delta_phi", &Delta_phi, has_Delta_phi);

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

        bD("p1_theta", &p1_theta, has_p1_theta);
        bD("p1_phi", &p1_phi, has_p1_phi);

        bD("p2_p", &p2_p, has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi", &p2_phi, has_p2_phi);

        if (!(has_runnum && has_detector1 && has_detector2 && has_x && has_Q2 &&
              has_t1 && has_phi2 && has_open_angle_ep2 && has_p1_theta && has_p1_phi)) {
            fatal("tree is missing one or more required branches: runnum, detector1, detector2, x, Q2, t1, phi2, open_angle_ep2, p1_theta, p1_phi.");
        }
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

static bool passes_global_cuts_for_event(const Branches& b,
                                         const PeriodTags& tags) {
    const GlobalCutConfig& cfg = default_global_cuts();

    if (cfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (is_excluded_run(b.runnum)) return false;

    if (global_cuts_require_sector_phi(cfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            fatal("sector selection requires e_phi, p1_phi, and p2_phi branches.");
        }
    }

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("P2/ycol cut is enabled but required branches are missing.");
        }

        if (global_cuts_require_sector_phi(cfg)) {
            return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                      b.detector1, b.detector2,
                                      tags.period_label,
                                      b.e_p, b.e_theta, b.e_phi, b.p1_phi,
                                      b.p2_p, b.p2_theta, b.p2_phi,
                                      cfg);
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  tags.period_label,
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

static bool passes_sigma_cuts_for_event(const std::string& key,
                                        const TopoCutMap& cuts,
                                        bool is_eppi0,
                                        const Branches& b) {
    auto it = cuts.find(key);
    if (it == cuts.end()) {
        fatal("missing sigma-cut key: " + key);
    }

    const CutVarMap& vm = it->second;

    if (!check_sigma(vm, "Emiss2", b.has_Emiss2, b.Emiss2)) return false;
    if (!check_sigma(vm, "Mx2", b.has_Mx2, b.Mx2)) return false;
    if (!check_sigma(vm, "Mx2_1", b.has_Mx2_1, b.Mx2_1)) return false;
    if (!check_sigma(vm, "Mx2_2", b.has_Mx2_2, b.Mx2_2)) return false;
    bool has_delta_phi = false;
    const double delta_phi = b.delta_phi_value(has_delta_phi);
    if (!check_sigma(vm, "Delta_phi", has_delta_phi, delta_phi)) return false;
    if (!check_sigma(vm, "pTmiss", b.has_pTmiss, b.pTmiss)) return false;
    if (!check_sigma(vm, "xF", b.has_xF, b.xF)) return false;

    if (is_eppi0) {
        if (!check_sigma(vm, "theta_pi0_pi0", b.has_theta_pi0_pi0, b.theta_pi0_pi0)) return false;
    } else {
        if (!check_sigma(vm, "theta_gamma_gamma", b.has_theta_gamma_gamma, b.theta_gamma_gamma)) return false;
    }

    return true;
}

// -----------------------------------------------------------------------------
// Totals
// -----------------------------------------------------------------------------

struct Totals {
    std::map<std::string, std::map<int, double>> dvcs_pi0_sub_by_period_current;
    std::map<std::string, std::map<int, double>> eppi0_norm_by_period_current;

    std::map<std::string, double> dvcs_uncategorized;
    std::map<std::string, double> eppi0_uncategorized;
    std::map<std::string, std::set<int>> dvcs_uncategorized_runs;
    std::map<std::string, std::set<int>> eppi0_uncategorized_runs;

    long long dvcs_events_used = 0;
    long long eppi0_events_used = 0;
};

static void process_tree_set(const std::map<std::string, TTree*>& trees,
                             bool is_eppi0,
                             const AnalysisInputs& inputs,
                             const TopoCutMap& data_cuts,
                             Totals& totals) {
    for (const auto& kv : trees) {
        const std::string& tree_key = kv.first;
        TTree* tree = kv.second;
        if (!tree) {
            continue;
        }

        PeriodTags tags = parse_period_from_key(tree_key);
        if (tags.display == "Fa18 Inb Supp") {
            std::cout << "[yield_totals] Skipping " << tree_key
                      << " because Fa18 Inb Supp is not part of the five-period CSV summary.\n";
            continue;
        }

        auto prod_it = inputs.products_by_period.find(tags.display);
        if (prod_it == inputs.products_by_period.end()) {
            fatal("missing analysis products for period " + tags.display);
        }
        const PeriodProducts& products = prod_it->second;

        auto cont_it = inputs.contamination_by_period.find(tags.display);
        if (cont_it == inputs.contamination_by_period.end()) {
            fatal("missing contamination vector for period " + tags.display);
        }
        const std::vector<double>& contamination = cont_it->second;

        Branches b;
        b.bind(tree);

        const Long64_t N = tree->GetEntries();
        long long accepted = 0;
        long long matched = 0;

        std::cout << "[yield_totals] Processing " << (is_eppi0 ? "eppi0" : "DVCS")
                  << " DATA period=" << tags.display
                  << " entries=" << static_cast<long long>(N) << "\n";

        for (Long64_t i = 0; i < N; ++i) {
            if (tree->GetEntry(i) <= 0) {
                continue;
            }

            if (!passes_global_cuts_for_event(b, tags)) {
                continue;
            }

            const std::string topo = topo_dir_from_detectors(b.detector1, b.detector2);
            if (topo.empty()) {
                continue;
            }

            const std::string cut_key = std::string(is_eppi0 ? "eppi0_" : "DVCS_") +
                                        tags.period_code + "_" + topo;

            if (!passes_sigma_cuts_for_event(cut_key, data_cuts, is_eppi0, b)) {
                continue;
            }

            ++accepted;

            const int row = find_row_index(inputs.rows, b.x, b.Q2, b.t_abs(), b.phi_deg());
            if (row < 0) {
                continue;
            }
            ++matched;

            const double theta = b.p1_theta_deg();
            const std::string region = normalization_region_for_event(theta, b.p1_phi_deg());
            if (region.empty()) {
                continue;
            }

            auto region_it = products.regions.find(region);
            if (region_it == products.regions.end()) {
                fatal("missing regional eppi0 normalization cubic for period " + tags.display + ", region " + region);
            }

            const double r_pi0 = eval_cubic(region_it->second.cubic, theta);
            if (!(std::isfinite(r_pi0) && r_pi0 > 0.0)) {
                fatal("non-positive regional eppi0 normalization cubic value for period " + tags.display + ", region " + region);
            }

            const double eff = is_eppi0 ? products.eppi0_current_eff : products.epg_current_eff;
            const double base_weight = 1.0 / (eff * r_pi0);

            double final_weight = base_weight;
            if (!is_eppi0) {
                const double c = contamination.at(static_cast<size_t>(row));
                final_weight *= (1.0 - c);
            }

            int current = 0;
            if (!resolve_current(tags.internal, b.runnum, current)) {
                if (is_eppi0) {
                    totals.eppi0_uncategorized[tags.display] += final_weight;
                    totals.eppi0_uncategorized_runs[tags.display].insert(b.runnum);
                } else {
                    totals.dvcs_uncategorized[tags.display] += final_weight;
                    totals.dvcs_uncategorized_runs[tags.display].insert(b.runnum);
                }
                continue;
            }

            if (is_eppi0) {
                totals.eppi0_norm_by_period_current[tags.display][current] += final_weight;
                totals.eppi0_events_used += 1;
            } else {
                totals.dvcs_pi0_sub_by_period_current[tags.display][current] += final_weight;
                totals.dvcs_events_used += 1;
            }
        }

        std::cout << "[yield_totals]   accepted_after_cuts=" << accepted
                  << " matched_bins=" << matched << "\n";
    }
}

static void write_one_table(std::ostream& os,
                            const std::string& title,
                            const std::map<std::string, std::map<int, double>>& values,
                            const std::map<std::string, double>& uncategorized,
                            const std::map<std::string, std::set<int>>& uncategorized_runs) {
    os << title << "\n";

    std::map<int, double> global_by_current;
    double global_total = 0.0;

    for (const std::string& period : period_order()) {
        os << "  Period: " << period << "\n";

        auto it = values.find(period);
        if (it == values.end() || it->second.empty()) {
            os << "    [no categorized events]\n";
        } else {
            double period_total = 0.0;
            for (const auto& ck : it->second) {
                os << "    current " << ck.first << " nA: "
                   << std::fixed << std::setprecision(6) << ck.second << "\n";
                global_by_current[ck.first] += ck.second;
                period_total += ck.second;
                global_total += ck.second;
            }
            os << "    period total: " << std::fixed << std::setprecision(6)
               << period_total << "\n";
        }

        auto iu = uncategorized.find(period);
        if (iu != uncategorized.end() && iu->second != 0.0) {
            os << "    current (uncategorized): " << std::fixed << std::setprecision(6)
               << iu->second << "\n";
            auto ir = uncategorized_runs.find(period);
            if (ir != uncategorized_runs.end()) {
                os << "    uncategorized runs:";
                for (int run : ir->second) {
                    os << " " << run;
                }
                os << "\n";
            }
        }

        os << "\n";
    }

    os << "  Global totals by current:\n";
    if (global_by_current.empty()) {
        os << "    [no categorized events]\n";
    } else {
        for (const auto& ck : global_by_current) {
            os << "    current " << ck.first << " nA: "
               << std::fixed << std::setprecision(6) << ck.second << "\n";
        }
    }

    os << "  Global categorized total: " << std::fixed << std::setprecision(6)
       << global_total << "\n\n";
}

static void write_csv_summary(const std::string& output_txt,
                              const Totals& totals) {
    std::string csv_path = output_txt;
    const size_t dot = csv_path.find_last_of('.');
    if (dot != std::string::npos) {
        csv_path = csv_path.substr(0, dot);
    }
    csv_path += ".csv";

    ensure_parent_directory_exists(csv_path);

    std::ofstream out(csv_path);
    if (!out.is_open()) {
        fatal("cannot write CSV summary: " + csv_path);
    }

    out << "quantity,period,current_nA,value\n";

    auto write_map = [&](const std::string& quantity,
                         const std::map<std::string, std::map<int, double>>& values) {
        for (const auto& pk : values) {
            for (const auto& ck : pk.second) {
                out << quantity << "," << pk.first << "," << ck.first << ","
                    << std::setprecision(12) << ck.second << "\n";
            }
        }
    };

    write_map("dvcs_normalized_pi0_subtracted", totals.dvcs_pi0_sub_by_period_current);
    write_map("eppi0_normalized", totals.eppi0_norm_by_period_current);

    std::cout << "[yield_totals] Wrote CSV summary to " << csv_path << "\n";
}

static void write_totals(std::ostream& os, const Totals& totals) {
    os << "================ Normalized yield totals by current ================\n\n";
    os << "Definitions:\n";
    os << "  DVCS normalized pi0-subtracted counts = normalized ep->epgamma DATA weights multiplied by (1 - contamination ratio).\n";
    os << "  eppi0 normalized counts = normalized ep->eppi0 DATA weights.\n";
    os << "  Both quantities use current-efficiency factors and the regional eppi0 AAOGEN p1_theta normalization cubic fits from the CSV.\n\n";

    os << "Raw accepted event counters used only as diagnostics:\n";
    os << "  DVCS events contributing: " << totals.dvcs_events_used << "\n";
    os << "  eppi0 events contributing: " << totals.eppi0_events_used << "\n\n";

    write_one_table(os,
                    "DVCS normalized pi0-subtracted counts by period/current:",
                    totals.dvcs_pi0_sub_by_period_current,
                    totals.dvcs_uncategorized,
                    totals.dvcs_uncategorized_runs);

    write_one_table(os,
                    "eppi0 normalized counts by period/current:",
                    totals.eppi0_norm_by_period_current,
                    totals.eppi0_uncategorized,
                    totals.eppi0_uncategorized_runs);

    os << "=====================================================================\n";
}

} // namespace

bool compute_yield_totals(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dvcsDataTrees,
                          const std::map<std::string, TTree*>& eppi0DataTrees,
                          const std::string& combined_cuts_json,
                          const std::string& output_txt) {
    try {
        ROOT::EnableThreadSafety();

        CSV csv;
        load_csv(csv_path, csv);
        const AnalysisInputs inputs = build_analysis_inputs(csv);
        const TopoCutMap data_cuts = load_sigma_cuts(combined_cuts_json, "data");

        if (data_cuts.empty()) {
            fatal("no data sigma cuts loaded from " + combined_cuts_json);
        }

        Totals totals;
        process_tree_set(dvcsDataTrees, false, inputs, data_cuts, totals);
        process_tree_set(eppi0DataTrees, true, inputs, data_cuts, totals);

        ensure_parent_directory_exists(output_txt);

        std::ofstream out(output_txt);
        if (!out.is_open()) {
            fatal("cannot write output text file: " + output_txt);
        }

        write_totals(out, totals);
        out.close();

        std::cout << "[yield_totals] Wrote text summary to " << output_txt << "\n";
        write_csv_summary(output_txt, totals);
        write_totals(std::cout, totals);

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[yield_totals] EXCEPTION: " << e.what() << std::endl;
        return false;
    }
}
