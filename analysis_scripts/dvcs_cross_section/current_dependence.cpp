#include "current_dependence.h"

#include "global_cuts.h"

#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH1D.h>
#include <TStyle.h>

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
static constexpr double MIN_E_GAMMA_CONE_ANGLE_DEG = 7.0;
static constexpr double COS_MIN_E_GAMMA_CONE =
    0.9925461516413220; // cos(7 deg)

static const std::vector<std::string> PERIOD_ORDER = {
    "Sp18 Inb",
    "Sp18 Out",
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb"
};

static const std::vector<std::string> CSV_PERIOD_ORDER = {
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out"
};

enum class Channel {
    DVCS,
    EPPI0
};

struct ChannelConfig {
    Channel channel = Channel::DVCS;
    std::string csv_channel;
    std::string cut_prefix;
    std::string output_token;
    std::string title;
};

static ChannelConfig dvcs_config() {
    ChannelConfig c;
    c.channel = Channel::DVCS;
    c.csv_channel = "ep->epg";
    c.cut_prefix = "DVCS";
    c.output_token = "epg";
    c.title = "ep #rightarrow ep#gamma";
    return c;
}

static ChannelConfig eppi0_config() {
    ChannelConfig c;
    c.channel = Channel::EPPI0;
    c.csv_channel = "ep->eppi0";
    c.cut_prefix = "eppi0";
    c.output_token = "eppi0";
    c.title = "ep #rightarrow ep#pi_{0}";
    return c;
}

struct PeriodTags {
    std::string display;
    std::string period_label;
    std::string period_code;
    std::string internal;
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

    std::ostringstream ss;
    ss << "[current_dependence] FATAL: cannot parse period from tree key: " << key;
    fatal(ss.str());

    return t;
}

static int reference_current_nA(const std::string& period_display) {
    if (period_display == "Sp18 Out") {
        return 45;
    }

    return 50;
}

static int parse_current_from_key(const std::string& key) {
    const std::string s = lower_ascii(key);

    if (has_substr(s, "nobkg")) {
        return 0;
    }

    // Parse an exact current token of the form <digits>nA.
    // Do not use substring matching such as "5na", because that incorrectly
    // classifies 35nA, 45nA, 50nA, and 55nA as 5nA, and 40nA as 0nA.
    for (size_t na_pos = s.find("na"); na_pos != std::string::npos; na_pos = s.find("na", na_pos + 2)) {
        if (na_pos == 0) {
            continue;
        }

        size_t first_digit = na_pos;
        while (first_digit > 0 && std::isdigit((unsigned char)s[first_digit - 1])) {
            --first_digit;
        }

        if (first_digit == na_pos) {
            continue;
        }

        const std::string current_token = s.substr(first_digit, na_pos - first_digit);
        char* endp = nullptr;
        const long current = std::strtol(current_token.c_str(), &endp, 10);

        if (endp == current_token.c_str() || *endp != '\0') {
            continue;
        }

        if (current < 0 || current > 100) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: parsed unreasonable current "
               << current << " nA from key: " << key;
            fatal(ss.str());
        }

        return (int)current;
    }

    // The nominal acceptance files are the reference-current files.
    PeriodTags tags = parse_period_from_key(key);
    return reference_current_nA(tags.display);
}

static std::string period_dir(const std::string& period) {
    if (period == "Fa18 Inb") return "Fa18_Inb";
    if (period == "Fa18 Out") return "Fa18_Out";
    if (period == "Sp19 Inb") return "Sp19_Inb";
    if (period == "Sp18 Inb") return "Sp18_Inb";
    if (period == "Sp18 Out") return "Sp18_Out";
    if (period == "Fa18 Inb Supp") return "Fa18_Inb_Supp";

    std::ostringstream ss;
    ss << "[current_dependence] FATAL: unknown period label: " << period;
    fatal(ss.str());

    return "";
}

static void mkdir_p(const std::string& path) {
    if (!path.empty()) {
        gSystem->mkdir(path.c_str(), true);
    }
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

        const bool quote =
            s.find(',') != std::string::npos ||
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
        ss << "[current_dependence] FATAL: cannot open CSV: " << path;
        fatal(ss.str());
    }

    std::string line;

    if (!std::getline(fin, line)) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: empty CSV: " << path;
        fatal(ss.str());
    }

    csv.header = split_csv_line(line);
    csv.index.clear();

    for (int i = 0; i < (int)csv.header.size(); ++i) {
        if (csv.index.find(csv.header[i]) != csv.index.end()) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: duplicate CSV column: " << csv.header[i];
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
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: row width mismatch in CSV: " << path;
            fatal(ss.str());
        }

        csv.rows.push_back(std::move(row));
    }
}

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";

    std::ofstream fout(tmp);

    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot write temp CSV: " << tmp;
        fatal(ss.str());
    }

    fout << join_csv_row(csv.header) << "\n";

    for (const auto& row : csv.rows) {
        fout << join_csv_row(row) << "\n";
    }

    fout.close();

    if (!fout) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: failed writing temp CSV: " << tmp;
        fatal(ss.str());
    }

    (void)std::remove(path.c_str());

    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: failed replacing CSV: " << path;
        fatal(ss.str());
    }
}

static int col_strict(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);

    if (it == csv.index.end()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: missing CSV column: " << name;
        fatal(ss.str());
    }

    return it->second;
}

static std::string tuple2(double v, double e) {
    std::ostringstream ss;
    ss << std::setprecision(12) << "(" << v << "," << e << ")";
    return ss.str();
}

static std::string current_eff_col(const ChannelConfig& cfg,
                                   const std::string& sample,
                                   const std::string& period) {
    std::ostringstream ss;
    ss << "current efficiency factor, " << cfg.csv_channel
       << ", " << sample << ", " << period;
    return ss.str();
}

// -----------------------------------------------------------------------------
// Charge and current maps
// -----------------------------------------------------------------------------

static std::unordered_map<int, double> read_charge_csv(const std::string& path) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot open charge CSV: " << path;
        fatal(ss.str());
    }

    std::unordered_map<int, double> out;
    std::string line;

    while (std::getline(fin, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            continue;
        }

        std::vector<std::string> fields = split_csv_line(line);

        if (fields.size() < 2) {
            continue;
        }

        char* endp = nullptr;
        const int run = (int)std::strtol(fields[0].c_str(), &endp, 10);

        if (endp == fields[0].c_str()) {
            continue;
        }

        endp = nullptr;
        const double charge = std::strtod(fields[1].c_str(), &endp);

        if (endp == fields[1].c_str()) {
            continue;
        }

        out[run] = charge;
    }

    std::cout << "[current_dependence] Loaded " << out.size()
              << " run-charge entries from " << path << std::endl;

    return out;
}

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
        {5443,5},{5610,5},
        {5444,20},
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
// Fits and uncertainty propagation
// -----------------------------------------------------------------------------

struct FitResult {
    double m = std::numeric_limits<double>::quiet_NaN();
    double b = std::numeric_limits<double>::quiet_NaN();
    double sm = std::numeric_limits<double>::quiet_NaN();
    double sb = std::numeric_limits<double>::quiet_NaN();
    double cov_mb = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndof = 0;
};

static FitResult weighted_linear_fit(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     const std::vector<double>& sy) {
    FitResult out;

    if (x.empty()) {
        return out;
    }

    if (x.size() == 1) {
        out.m = 0.0;
        out.b = y[0];
        out.sm = 0.0;
        out.sb = sy[0] > 0.0 ? sy[0] : 0.0;
        out.cov_mb = 0.0;
        out.chi2 = 0.0;
        out.ndof = 0;
        return out;
    }

    double S = 0.0;
    double Sx = 0.0;
    double Sy = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        if (!(sy[i] > 0.0)) {
            fatal("[current_dependence] FATAL: non-positive uncertainty in weighted fit.");
        }

        const double w = 1.0 / (sy[i] * sy[i]);

        S += w;
        Sx += w * x[i];
        Sy += w * y[i];
        Sxx += w * x[i] * x[i];
        Sxy += w * x[i] * y[i];
    }

    const double D = S * Sxx - Sx * Sx;

    if (D == 0.0) {
        return out;
    }

    out.m = (S * Sxy - Sx * Sy) / D;
    out.b = (Sxx * Sy - Sx * Sxy) / D;
    out.sm = std::sqrt(S / D);
    out.sb = std::sqrt(Sxx / D);
    out.cov_mb = -Sx / D;
    out.ndof = (int)x.size() - 2;

    out.chi2 = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        const double yf = out.m * x[i] + out.b;
        const double r = (y[i] - yf) / sy[i];
        out.chi2 += r * r;
    }

    return out;
}

static double rel_at_current(double current, const FitResult& f) {
    if (!std::isfinite(f.m) || !std::isfinite(f.b) || f.b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return (f.m * current + f.b) / f.b;
}

static double rel_err_at_current(double current, const FitResult& f) {
    if (!std::isfinite(f.m) || !std::isfinite(f.b) || f.b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (!std::isfinite(f.sm) || !std::isfinite(f.sb) || !std::isfinite(f.cov_mb)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double var_m = f.sm * f.sm;
    const double var_b = f.sb * f.sb;

    const double dr_dm = current / f.b;
    const double dr_db = -(f.m * current) / (f.b * f.b);

    double var =
        dr_dm * dr_dm * var_m +
        dr_db * dr_db * var_b +
        2.0 * dr_dm * dr_db * f.cov_mb;

    if (var < 0.0 && std::fabs(var) < 1.0e-15) {
        var = 0.0;
    }

    if (var < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(var);
}

// -----------------------------------------------------------------------------
// Cuts
// -----------------------------------------------------------------------------

struct SigmaStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std = std::numeric_limits<double>::quiet_NaN();
};

using CutVarMap = std::unordered_map<std::string, SigmaStats>;
using TopoCutMap = std::unordered_map<std::string, CutVarMap>;

static TopoCutMap load_sigma_cuts(const std::string& path,
                                  const std::string& sample_key) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot open combined cuts JSON: " << path;
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

static bool check_sigma_var(const CutVarMap& vm,
                            const std::string& var,
                            bool has_value,
                            double value) {
    auto it = vm.find(var);

    if (it == vm.end()) {
        return true;
    }

    if (!has_value) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: required sigma-cut branch missing: " << var;
        fatal(ss.str());
    }

    return within_3sigma(value, it->second);
}

static std::string topo_dir(int detector1, int detector2) {
    if (detector1 == 1 && detector2 == 1) return "FD_FD";
    if (detector1 == 2 && detector2 == 1) return "CD_FD";
    if (detector1 == 2 && detector2 == 0) return "CD_FT";
    return "";
}

// -----------------------------------------------------------------------------
// Branch binding
// -----------------------------------------------------------------------------

static std::mutex g_bind_mutex;

struct Branches {
    int runnum = 0; bool has_runnum = false;
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    double t1 = 0.0; bool has_t1 = false;
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

    double p1_theta = 0.0; bool has_p1_theta = false;

    double p2_p = 0.0; bool has_p2_p = false;
    double p2_theta = 0.0; bool has_p2_theta = false;
    double p2_phi = 0.0; bool has_p2_phi = false;

    void bind(TTree* t) {
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

        ena("t1");
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

        ena("p2_p");
        ena("p2_theta");
        ena("p2_phi");

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

        bD("t1", &t1, has_t1);
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

        bD("p1_theta", &p1_theta, has_p1_theta);

        bD("p2_p", &p2_p, has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi", &p2_phi, has_p2_phi);
    }
};

static bool passes_cone_cut(const Branches& b) {
    if (!(b.has_e_theta && b.has_e_phi && b.has_p2_theta && b.has_p2_phi)) {
        return false;
    }

    const double st_e = std::sin(b.e_theta);
    const double ct_e = std::cos(b.e_theta);
    const double sp_e = std::sin(b.e_phi);
    const double cp_e = std::cos(b.e_phi);

    const double st_g = std::sin(b.p2_theta);
    const double ct_g = std::cos(b.p2_theta);
    const double sp_g = std::sin(b.p2_phi);
    const double cp_g = std::cos(b.p2_phi);

    const double dot =
        st_e * cp_e * st_g * cp_g +
        st_e * sp_e * st_g * sp_g +
        ct_e * ct_g;

    return dot <= COS_MIN_E_GAMMA_CONE;
}

static bool passes_global_dispatch(const Branches& b,
                                   const PeriodTags& tags) {
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
            fatal("[current_dependence] FATAL: missing branches required by global ycol/topology cuts.");
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
        fatal("[current_dependence] FATAL: missing detector1/detector2.");
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
        ss << "[current_dependence] FATAL: missing sigma-cut key: " << key;
        fatal(ss.str());
    }

    const CutVarMap& vm = it->second;

    if (!check_sigma_var(vm, "Emiss2", b.has_Emiss2, b.Emiss2)) return false;
    if (!check_sigma_var(vm, "Mx2", b.has_Mx2, b.Mx2)) return false;
    if (!check_sigma_var(vm, "Mx2_1", b.has_Mx2_1, b.Mx2_1)) return false;
    if (!check_sigma_var(vm, "Mx2_2", b.has_Mx2_2, b.Mx2_2)) return false;
    if (!check_sigma_var(vm, "pTmiss", b.has_pTmiss, b.pTmiss)) return false;
    if (!check_sigma_var(vm, "xF", b.has_xF, b.xF)) return false;

    if (cfg.channel == Channel::EPPI0) {
        if (!check_sigma_var(vm, "theta_pi0_pi0", b.has_theta_pi0_pi0, b.theta_pi0_pi0)) return false;
    } else {
        if (!check_sigma_var(vm, "theta_gamma_gamma", b.has_theta_gamma_gamma, b.theta_gamma_gamma)) return false;
    }

    return true;
}

// -----------------------------------------------------------------------------
// Data and MC aggregations
// -----------------------------------------------------------------------------

struct CurrentPoint {
    int current_nA = 0;
    double y = 0.0;
    double sy = 0.0;
    double counts = 0.0;
    double charge = 0.0;
};

struct PeriodResult {
    std::string period;
    FitResult data_fit;
    FitResult mc_fit;

    double data_factor = std::numeric_limits<double>::quiet_NaN();
    double data_factor_err = std::numeric_limits<double>::quiet_NaN();

    double mc_factor = std::numeric_limits<double>::quiet_NaN();
    double mc_factor_err = std::numeric_limits<double>::quiet_NaN();

    std::vector<CurrentPoint> data_points;
    std::vector<CurrentPoint> mc_points;
};

struct DataAgg {
    std::string period;
    std::map<int, long long> counts_by_run;
    std::map<int, int> current_by_run;
    std::map<int, double> charge_by_run;
};

static DataAgg process_data_tree(const ChannelConfig& cfg,
                                 const std::string& key,
                                 TTree* tree,
                                 const std::unordered_map<int, double>& charge_map,
                                 const TopoCutMap& data_cuts) {
    PeriodTags tags = parse_period_from_key(key);

    DataAgg out;
    out.period = tags.display;

    if (!tree) {
        return out;
    }

    Branches b;
    b.bind(tree);

    if (!b.has_runnum) {
        fatal("[current_dependence] FATAL: data tree missing runnum.");
    }

    const Long64_t N = tree->GetEntries();

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (!passes_cone_cut(b)) {
            continue;
        }

        int current = 0;

        if (!resolve_current(tags.internal, b.runnum, current)) {
            continue;
        }

        auto charge_it = charge_map.find(b.runnum);

        if (charge_it == charge_map.end()) {
            continue;
        }

        const double charge = charge_it->second;

        if (!(charge > 0.0)) {
            continue;
        }

        if (!passes_global_dispatch(b, tags)) {
            continue;
        }

        if (!passes_sigma_dispatch(cfg, tags, data_cuts, b)) {
            continue;
        }

        out.counts_by_run[b.runnum] += 1;
        out.current_by_run[b.runnum] = current;
        out.charge_by_run[b.runnum] = charge;
    }

    std::cout << "[current_dependence] DATA channel=" << cfg.csv_channel
              << " period=" << out.period
              << " entries=" << (long long)N
              << " valid_runs=" << out.counts_by_run.size()
              << std::endl;

    return out;
}

struct McAgg {
    std::string period;
    int current_nA = 0;
    long long n_gen = 0;
    long long n_rec = 0;
};

static long long count_generated_tree(TTree* tree) {
    if (!tree) {
        return 0;
    }

    // Generated trees in this workflow are already generated-level samples.
    // The Python script counted generated entries directly.
    return (long long)tree->GetEntries();
}

static long long count_reconstructed_tree(const ChannelConfig& cfg,
                                          const std::string& key,
                                          TTree* tree,
                                          const TopoCutMap& mc_cuts) {
    PeriodTags tags = parse_period_from_key(key);

    if (!tree) {
        return 0;
    }

    Branches b;
    b.bind(tree);

    const Long64_t N = tree->GetEntries();
    long long count = 0;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (!passes_cone_cut(b)) {
            continue;
        }

        if (!passes_global_dispatch(b, tags)) {
            continue;
        }

        if (!passes_sigma_dispatch(cfg, tags, mc_cuts, b)) {
            continue;
        }

        ++count;
    }

    return count;
}

static std::vector<CurrentPoint> data_points_from_agg(const DataAgg& agg) {
    struct CurrentGroup {
        long long counts = 0;
        double charge = 0.0;
        int n_runs = 0;
    };

    std::map<int, CurrentGroup> grouped;

    for (const auto& kv : agg.counts_by_run) {
        const int run = kv.first;
        const long long counts = kv.second;

        auto it_cur = agg.current_by_run.find(run);
        auto it_chg = agg.charge_by_run.find(run);

        if (it_cur == agg.current_by_run.end() || it_chg == agg.charge_by_run.end()) {
            continue;
        }

        CurrentGroup& g = grouped[it_cur->second];
        g.counts += counts;
        g.charge += it_chg->second;
        g.n_runs += 1;
    }

    std::vector<CurrentPoint> out;

    for (const auto& kv : grouped) {
        const int current = kv.first;
        const CurrentGroup& g = kv.second;

        if (!(g.charge > 0.0)) {
            continue;
        }

        CurrentPoint p;
        p.current_nA = current;
        p.counts = (double)g.counts;
        p.charge = g.charge;
        p.y = (double)g.counts / g.charge;
        p.sy = (g.counts > 0) ? std::sqrt((double)g.counts) / g.charge : 0.0;

        if (p.sy <= 0.0) {
            p.sy = 1.0 / g.charge;
        }

        out.push_back(p);
    }

    std::sort(out.begin(), out.end(), [](const CurrentPoint& a, const CurrentPoint& b) {
        return a.current_nA < b.current_nA;
    });

    return out;
}

static std::vector<CurrentPoint> mc_points_from_aggs(const std::vector<McAgg>& aggs,
                                                     const std::string& period) {
    std::vector<CurrentPoint> out;

    for (const McAgg& a : aggs) {
        if (a.period != period) {
            continue;
        }

        if (a.n_gen <= 0) {
            continue;
        }

        const double eff = (double)a.n_rec / (double)a.n_gen;
        const double err = std::sqrt(eff * (1.0 - eff) / (double)a.n_gen);

        CurrentPoint p;
        p.current_nA = a.current_nA;
        p.counts = (double)a.n_rec;
        p.charge = (double)a.n_gen;
        p.y = eff;
        p.sy = (err > 0.0) ? err : 1.0 / (double)a.n_gen;

        out.push_back(p);
    }

    std::sort(out.begin(), out.end(), [](const CurrentPoint& a, const CurrentPoint& b) {
        return a.current_nA < b.current_nA;
    });

    return out;
}

static FitResult fit_points(const std::vector<CurrentPoint>& points) {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> sy;

    for (const CurrentPoint& p : points) {
        if (!std::isfinite(p.y) || !std::isfinite(p.sy) || p.sy <= 0.0) {
            continue;
        }

        x.push_back((double)p.current_nA);
        y.push_back(p.y);
        sy.push_back(p.sy);
    }

    return weighted_linear_fit(x, y, sy);
}

static double weighted_data_rel(const std::vector<CurrentPoint>& points,
                                const FitResult& fit) {
    double total_counts = 0.0;
    double weighted_current = 0.0;

    for (const CurrentPoint& p : points) {
        total_counts += p.counts;
    }

    if (!(total_counts > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    for (const CurrentPoint& p : points) {
        weighted_current += (p.counts / total_counts) * (double)p.current_nA;
    }

    return rel_at_current(weighted_current, fit);
}

static double weighted_data_rel_err(const std::vector<CurrentPoint>& points,
                                    const FitResult& fit) {
    double total_counts = 0.0;
    double weighted_current = 0.0;

    for (const CurrentPoint& p : points) {
        total_counts += p.counts;
    }

    if (!(total_counts > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    for (const CurrentPoint& p : points) {
        weighted_current += (p.counts / total_counts) * (double)p.current_nA;
    }

    return rel_err_at_current(weighted_current, fit);
}

// -----------------------------------------------------------------------------
// Plots
// -----------------------------------------------------------------------------

static void draw_fit_graph(const std::string& out_path,
                           const std::string& title,
                           const std::vector<CurrentPoint>& data_points,
                           const FitResult& data_fit,
                           const std::vector<CurrentPoint>& mc_points,
                           const FitResult& mc_fit,
                           bool relative) {
    mkdir_p(out_path.substr(0, out_path.find_last_of('/')));

    TCanvas c("c_current_dependence", "", 1100, 800);
    c.SetGrid(1, 1);
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.08);

    TGraphErrors g_data;
    TGraphErrors g_mc;

    for (int i = 0; i < (int)data_points.size(); ++i) {
        double y = data_points[i].y;
        double ey = data_points[i].sy;

        if (relative) {
            if (std::isfinite(data_fit.b) && data_fit.b != 0.0) {
                y /= data_fit.b;
                ey /= std::fabs(data_fit.b);
            }
        }

        g_data.SetPoint(i, data_points[i].current_nA, y);
        g_data.SetPointError(i, 0.0, ey);
    }

    for (int i = 0; i < (int)mc_points.size(); ++i) {
        double y = mc_points[i].y;
        double ey = mc_points[i].sy;

        if (relative) {
            if (std::isfinite(mc_fit.b) && mc_fit.b != 0.0) {
                y /= mc_fit.b;
                ey /= std::fabs(mc_fit.b);
            }
        }

        g_mc.SetPoint(i, mc_points[i].current_nA, y);
        g_mc.SetPointError(i, 0.0, ey);
    }

    g_data.SetMarkerStyle(20);
    g_data.SetMarkerColor(kBlack);
    g_data.SetLineColor(kBlack);

    g_mc.SetMarkerStyle(24);
    g_mc.SetMarkerColor(kRed + 1);
    g_mc.SetLineColor(kRed + 1);

    double ymax = 0.0;

    for (const CurrentPoint& p : data_points) {
        double y = p.y;
        double ey = p.sy;
        if (relative && std::isfinite(data_fit.b) && data_fit.b != 0.0) {
            y /= data_fit.b;
            ey /= std::fabs(data_fit.b);
        }
        ymax = std::max(ymax, y + ey);
    }

    for (const CurrentPoint& p : mc_points) {
        double y = p.y;
        double ey = p.sy;
        if (relative && std::isfinite(mc_fit.b) && mc_fit.b != 0.0) {
            y /= mc_fit.b;
            ey /= std::fabs(mc_fit.b);
        }
        ymax = std::max(ymax, y + ey);
    }

    if (!(ymax > 0.0)) {
        ymax = 1.0;
    }

    TH1D* frame = (TH1D*)gPad->DrawFrame(0.0, 0.0, 80.0, 1.25 * ymax);
    frame->SetTitle(title.c_str());
    frame->GetXaxis()->SetTitle("Current (nA)");

    if (relative) {
        frame->GetYaxis()->SetTitle("Relative response to 0 nA");
    } else {
        frame->GetYaxis()->SetTitle("Data counts / nC or MC efficiency");
    }

    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);

    g_data.Draw("PE SAME");
    g_mc.Draw("PE SAME");

    auto draw_line = [&](const FitResult& fit, int color, bool is_relative) {
        if (!std::isfinite(fit.m) || !std::isfinite(fit.b)) {
            return;
        }

        TGraph* g = new TGraph();

        for (int i = 0; i < 200; ++i) {
            const double x = 80.0 * (double)i / 199.0;
            double y = fit.m * x + fit.b;

            if (is_relative && fit.b != 0.0) {
                y /= fit.b;
            }

            g->SetPoint(i, x, y);
        }

        g->SetLineColor(color);
        g->SetLineWidth(2);
        g->Draw("L SAME");
    };

    draw_line(data_fit, kBlack, relative);
    draw_line(mc_fit, kRed + 1, relative);

    TLegend leg(0.58, 0.72, 0.90, 0.88);
    leg.SetFillStyle(1001);
    leg.SetFillColor(kWhite);
    leg.SetBorderSize(1);
    leg.AddEntry(&g_data, "Data", "pe");
    leg.AddEntry(&g_mc, "MC", "pe");
    leg.Draw();

    c.SaveAs(out_path.c_str());
}



static int period_color(const std::string& period) {
    if (period == "Sp18 Inb") return kBlack;
    if (period == "Sp18 Out") return kRed + 1;
    if (period == "Fa18 Inb") return kBlue + 1;
    if (period == "Fa18 Out") return kGreen + 2;
    if (period == "Sp19 Inb") return kMagenta + 1;
    return kGray + 2;
}

static int period_marker(const std::string& period) {
    if (period == "Sp18 Inb") return 20;
    if (period == "Sp18 Out") return 21;
    if (period == "Fa18 Inb") return 22;
    if (period == "Fa18 Out") return 23;
    if (period == "Sp19 Inb") return 24;
    return 20;
}

static double point_y_for_plot(const CurrentPoint& p, const FitResult& fit, bool relative) {
    double y = p.y;
    if (relative && std::isfinite(fit.b) && fit.b != 0.0) {
        y /= fit.b;
    }
    return y;
}

static double point_ey_for_plot(const CurrentPoint& p, const FitResult& fit, bool relative) {
    double ey = p.sy;
    if (relative && std::isfinite(fit.b) && fit.b != 0.0) {
        ey /= std::fabs(fit.b);
    }
    return ey;
}

static double fit_y_for_plot(double current, const FitResult& fit, bool relative) {
    if (!std::isfinite(fit.m) || !std::isfinite(fit.b)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double y = fit.m * current + fit.b;
    if (relative && fit.b != 0.0) {
        y /= fit.b;
    }
    return y;
}

static void write_all_points_csv(const std::string& path,
                                 const std::vector<PeriodResult>& results,
                                 bool write_data_points) {
    mkdir_p(path.substr(0, path.find_last_of('/')));

    std::ofstream fout(path);

    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot write: " << path;
        fatal(ss.str());
    }

    fout << "period,current_nA,y,stat,counts,normalizer\n";

    for (const PeriodResult& r : results) {
        const std::vector<CurrentPoint>& points = write_data_points ? r.data_points : r.mc_points;
        for (const CurrentPoint& pt : points) {
            fout << r.period << ","
                 << pt.current_nA << ","
                 << std::setprecision(12) << pt.y << ","
                 << std::setprecision(12) << pt.sy << ","
                 << std::setprecision(12) << pt.counts << ","
                 << std::setprecision(12) << pt.charge << "\n";
        }
    }
}

static void draw_current_panel_2x3(const std::string& out_path,
                                   const std::string& title,
                                   const std::vector<PeriodResult>& results,
                                   bool draw_data,
                                   bool relative) {
    mkdir_p(out_path.substr(0, out_path.find_last_of('/')));

    TCanvas c("c_current_dependence_panel", "", 1800, 1200);
    c.Divide(3, 2, 0.001, 0.001);

    std::map<std::string, PeriodResult> by_period;
    for (const PeriodResult& r : results) {
        by_period[r.period] = r;
    }

    std::vector<TGraphErrors*> owned_graphs;
    std::vector<TGraph*> owned_lines;

    auto draw_one_period = [&](const PeriodResult& r, int color, int marker, const std::string& pad_title, bool final_pad) {
        const std::vector<CurrentPoint>& points = draw_data ? r.data_points : r.mc_points;
        const FitResult& fit = draw_data ? r.data_fit : r.mc_fit;

        TGraphErrors* g = new TGraphErrors();
        owned_graphs.push_back(g);
        g->SetMarkerStyle(marker);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->SetLineWidth(1);

        double ymax = 0.0;
        double ymin = std::numeric_limits<double>::infinity();

        for (int i = 0; i < (int)points.size(); ++i) {
            const double y = point_y_for_plot(points[i], fit, relative);
            const double ey = point_ey_for_plot(points[i], fit, relative);
            g->SetPoint(i, points[i].current_nA, y);
            g->SetPointError(i, 0.0, ey);
            if (std::isfinite(y)) {
                ymax = std::max(ymax, y + ey);
                ymin = std::min(ymin, y - ey);
            }
        }

        if (!(ymax > 0.0)) ymax = 1.0;
        if (!std::isfinite(ymin)) ymin = 0.0;
        if (relative) {
            ymin = std::min(0.85, ymin * 0.95);
            ymax = std::max(1.15, ymax * 1.05);
        } else {
            ymin = 0.0;
            ymax *= 1.25;
        }

        TH1D* frame = (TH1D*)gPad->DrawFrame(0.0, ymin, 80.0, ymax);
        frame->SetTitle(pad_title.c_str());
        frame->GetXaxis()->SetTitle("Current (nA)");
        if (relative) {
            frame->GetYaxis()->SetTitle("Relative response to 0 nA");
        } else if (draw_data) {
            frame->GetYaxis()->SetTitle("Data counts / nC");
        } else {
            frame->GetYaxis()->SetTitle("MC reconstructed / generated");
        }
        frame->GetXaxis()->CenterTitle(true);
        frame->GetYaxis()->CenterTitle(true);
        frame->GetXaxis()->SetTitleSize(0.050);
        frame->GetYaxis()->SetTitleSize(0.050);
        frame->GetXaxis()->SetLabelSize(0.045);
        frame->GetYaxis()->SetLabelSize(0.045);
        frame->GetYaxis()->SetTitleOffset(1.55);
        frame->GetXaxis()->SetTitleOffset(1.15);

        g->Draw("PE SAME");

        if (std::isfinite(fit.m) && std::isfinite(fit.b)) {
            TGraph* line = new TGraph();
            owned_lines.push_back(line);
            for (int i = 0; i < 200; ++i) {
                const double x = 80.0 * (double)i / 199.0;
                const double y = fit_y_for_plot(x, fit, relative);
                line->SetPoint(i, x, y);
            }
            line->SetLineColor(color);
            line->SetLineWidth(2);
            line->Draw("L SAME");
        }

        if (final_pad) {
            g->SetTitle(r.period.c_str());
        }
    };

    for (int ipad = 0; ipad < 6; ++ipad) {
        c.cd(ipad + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.14);
        gPad->SetGrid(1, 1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);

        if (ipad < 5) {
            const std::string period = PERIOD_ORDER[(size_t)ipad];
            auto it = by_period.find(period);
            if (it == by_period.end()) {
                continue;
            }
            draw_one_period(it->second, period_color(period), period_marker(period), period, false);
        } else {
            double ymax = 0.0;
            double ymin = std::numeric_limits<double>::infinity();

            for (const std::string& period : PERIOD_ORDER) {
                auto it = by_period.find(period);
                if (it == by_period.end()) continue;
                const PeriodResult& r = it->second;
                const std::vector<CurrentPoint>& points = draw_data ? r.data_points : r.mc_points;
                const FitResult& fit = draw_data ? r.data_fit : r.mc_fit;
                for (const CurrentPoint& pt : points) {
                    const double y = point_y_for_plot(pt, fit, relative);
                    const double ey = point_ey_for_plot(pt, fit, relative);
                    if (std::isfinite(y)) {
                        ymax = std::max(ymax, y + ey);
                        ymin = std::min(ymin, y - ey);
                    }
                }
            }

            if (!(ymax > 0.0)) ymax = 1.0;
            if (!std::isfinite(ymin)) ymin = 0.0;
            if (relative) {
                ymin = std::min(0.85, ymin * 0.95);
                ymax = std::max(1.15, ymax * 1.05);
            } else {
                ymin = 0.0;
                ymax *= 1.25;
            }

            TH1D* frame = (TH1D*)gPad->DrawFrame(0.0, ymin, 80.0, ymax);
            frame->SetTitle("All periods");
            frame->GetXaxis()->SetTitle("Current (nA)");
            if (relative) {
                frame->GetYaxis()->SetTitle("Relative response to 0 nA");
            } else if (draw_data) {
                frame->GetYaxis()->SetTitle("Data counts / nC");
            } else {
                frame->GetYaxis()->SetTitle("MC reconstructed / generated");
            }
            frame->GetXaxis()->CenterTitle(true);
            frame->GetYaxis()->CenterTitle(true);
            frame->GetXaxis()->SetTitleSize(0.050);
            frame->GetYaxis()->SetTitleSize(0.050);
            frame->GetXaxis()->SetLabelSize(0.045);
            frame->GetYaxis()->SetLabelSize(0.045);
            frame->GetYaxis()->SetTitleOffset(1.55);
            frame->GetXaxis()->SetTitleOffset(1.15);

            TLegend* leg = new TLegend(0.48, 0.60, 0.92, 0.90);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetBorderSize(1);
            leg->SetTextFont(42);
            leg->SetTextSize(0.038);

            for (const std::string& period : PERIOD_ORDER) {
                auto it = by_period.find(period);
                if (it == by_period.end()) continue;
                const PeriodResult& r = it->second;
                const std::vector<CurrentPoint>& points = draw_data ? r.data_points : r.mc_points;
                const FitResult& fit = draw_data ? r.data_fit : r.mc_fit;
                const int color = period_color(period);
                const int marker = period_marker(period);

                TGraphErrors* g = new TGraphErrors();
                owned_graphs.push_back(g);
                g->SetMarkerStyle(marker);
                g->SetMarkerColor(color);
                g->SetLineColor(color);
                g->SetLineWidth(1);

                for (int i = 0; i < (int)points.size(); ++i) {
                    const double y = point_y_for_plot(points[i], fit, relative);
                    const double ey = point_ey_for_plot(points[i], fit, relative);
                    g->SetPoint(i, points[i].current_nA, y);
                    g->SetPointError(i, 0.0, ey);
                }
                g->Draw("PE SAME");
                leg->AddEntry(g, period.c_str(), "pe");

                if (std::isfinite(fit.m) && std::isfinite(fit.b)) {
                    TGraph* line = new TGraph();
                    owned_lines.push_back(line);
                    for (int i = 0; i < 200; ++i) {
                        const double x = 80.0 * (double)i / 199.0;
                        const double y = fit_y_for_plot(x, fit, relative);
                        line->SetPoint(i, x, y);
                    }
                    line->SetLineColor(color);
                    line->SetLineWidth(2);
                    line->Draw("L SAME");
                }
            }
            leg->Draw();
        }
    }

    c.cd();
    TLatex title_text;
    title_text.SetNDC(true);
    title_text.SetTextFont(42);
    title_text.SetTextSize(0.032);
    title_text.DrawLatex(0.03, 0.985, title.c_str());

    c.SaveAs(out_path.c_str());

    for (TGraphErrors* g : owned_graphs) delete g;
    for (TGraph* g : owned_lines) delete g;
}

static void write_points_csv(const std::string& path,
                             const std::vector<CurrentPoint>& points) {
    mkdir_p(path.substr(0, path.find_last_of('/')));

    std::ofstream fout(path);

    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot write: " << path;
        fatal(ss.str());
    }

    fout << "current_nA,y,stat,counts,normalizer\n";

    for (const CurrentPoint& p : points) {
        fout << p.current_nA << ","
             << std::setprecision(12) << p.y << ","
             << std::setprecision(12) << p.sy << ","
             << std::setprecision(12) << p.counts << ","
             << std::setprecision(12) << p.charge << "\n";
    }
}

static void write_summary_csv(const std::string& path,
                              const std::vector<PeriodResult>& rows) {
    mkdir_p(path.substr(0, path.find_last_of('/')));

    std::ofstream fout(path);

    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot write: " << path;
        fatal(ss.str());
    }

    fout << "period,data_factor,data_factor_stat,mc_factor,mc_factor_stat,"
         << "data_m,data_b,mc_m,mc_b\n";

    for (const PeriodResult& r : rows) {
        fout << r.period << ","
             << std::setprecision(12) << r.data_factor << ","
             << std::setprecision(12) << r.data_factor_err << ","
             << std::setprecision(12) << r.mc_factor << ","
             << std::setprecision(12) << r.mc_factor_err << ","
             << std::setprecision(12) << r.data_fit.m << ","
             << std::setprecision(12) << r.data_fit.b << ","
             << std::setprecision(12) << r.mc_fit.m << ","
             << std::setprecision(12) << r.mc_fit.b << "\n";
    }
}

// -----------------------------------------------------------------------------
// Main channel worker
// -----------------------------------------------------------------------------


static void validate_current_scan_mc_inputs(const ChannelConfig& cfg,
                                            const std::map<std::string, TTree*>& gen_trees,
                                            const std::map<std::string, TTree*>& rec_trees) {
    if (cfg.channel != Channel::DVCS) {
        return;
    }

    std::map<std::string, std::set<int>> gen_currents;
    std::map<std::string, std::set<int>> rec_currents;

    for (const auto& kv : gen_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);
        if (tags.display == "Fa18 Inb Supp") continue;
        gen_currents[tags.display].insert(parse_current_from_key(kv.first));
    }

    for (const auto& kv : rec_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);
        if (tags.display == "Fa18 Inb Supp") continue;
        rec_currents[tags.display].insert(parse_current_from_key(kv.first));
    }

    for (const std::string& period : PERIOD_ORDER) {
        const size_t ngen = gen_currents[period].size();
        const size_t nrec = rec_currents[period].size();

        std::cout << "[current_dependence] DVCS current-scan MC inputs for " << period
                  << ": gen currents=";
        for (int c : gen_currents[period]) std::cout << " " << c;
        std::cout << " | rec currents=";
        for (int c : rec_currents[period]) std::cout << " " << c;
        std::cout << std::endl;

        if (ngen < 2 || nrec < 2) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: DVCS current-dependence MC for " << period
               << " has fewer than two current points (gen=" << ngen
               << ", rec=" << nrec << "). This usually means main.cpp passed the nominal "
               << "acceptance MC maps instead of currentStudyGenMcTrees/currentStudyRecMcTrees.";
            fatal(ss.str());
        }
    }
}

static std::vector<PeriodResult> run_channel_study(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& data_trees,
    const std::map<std::string, TTree*>& gen_trees,
    const std::map<std::string, TTree*>& rec_trees,
    const std::unordered_map<int, double>& charge_map,
    const TopoCutMap& data_cuts,
    const TopoCutMap& mc_cuts,
    const std::string& output_dir,
    int max_workers) {

    validate_current_scan_mc_inputs(cfg, gen_trees, rec_trees);

    std::map<std::string, DataAgg> data_aggs;
    std::mutex data_mutex;

    std::vector<std::pair<std::string, TTree*>> data_items;

    for (const auto& kv : data_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);

        if (tags.display == "Fa18 Inb Supp") {
            continue;
        }

        data_items.push_back(kv);
    }

    int nth = std::max(1, std::min(5, max_workers));
    nth = std::min(nth, std::max(1, (int)data_items.size()));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int i = 0; i < (int)data_items.size(); ++i) {
        const auto& item = data_items[i];

        DataAgg agg = process_data_tree(cfg,
                                        item.first,
                                        item.second,
                                        charge_map,
                                        data_cuts);

        std::lock_guard<std::mutex> lock(data_mutex);
        data_aggs[agg.period] = std::move(agg);
    }

    std::map<std::string, McAgg> mc_by_period_current;

    for (const auto& kv : gen_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);

        if (tags.display == "Fa18 Inb Supp") {
            continue;
        }

        const int current = parse_current_from_key(kv.first);

        std::ostringstream key;
        key << tags.display << "_" << current;

        McAgg& agg = mc_by_period_current[key.str()];
        agg.period = tags.display;
        agg.current_nA = current;
        agg.n_gen += count_generated_tree(kv.second);
    }

    std::vector<std::pair<std::string, TTree*>> rec_items;

    for (const auto& kv : rec_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);

        if (tags.display == "Fa18 Inb Supp") {
            continue;
        }

        rec_items.push_back(kv);
    }

    std::mutex rec_mutex;
    nth = std::max(1, std::min(5, max_workers));
    nth = std::min(nth, std::max(1, (int)rec_items.size()));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int i = 0; i < (int)rec_items.size(); ++i) {
        const auto& item = rec_items[i];

        PeriodTags tags = parse_period_from_key(item.first);
        const int current = parse_current_from_key(item.first);

        const long long n_rec = count_reconstructed_tree(cfg,
                                                         item.first,
                                                         item.second,
                                                         mc_cuts);

        std::ostringstream key;
        key << tags.display << "_" << current;

        std::lock_guard<std::mutex> lock(rec_mutex);
        McAgg& agg = mc_by_period_current[key.str()];
        agg.period = tags.display;
        agg.current_nA = current;
        agg.n_rec += n_rec;
    }

    std::vector<McAgg> mc_aggs;

    for (const auto& kv : mc_by_period_current) {
        mc_aggs.push_back(kv.second);
    }

    std::vector<PeriodResult> results;

    for (const std::string& period : PERIOD_ORDER) {
        PeriodResult R;
        R.period = period;

        auto it = data_aggs.find(period);

        if (it != data_aggs.end()) {
            R.data_points = data_points_from_agg(it->second);
            R.data_fit = fit_points(R.data_points);
            R.data_factor = weighted_data_rel(R.data_points, R.data_fit);
            R.data_factor_err = weighted_data_rel_err(R.data_points, R.data_fit);
        }

        R.mc_points = mc_points_from_aggs(mc_aggs, period);
        R.mc_fit = fit_points(R.mc_points);

        const int ref = reference_current_nA(period);
        R.mc_factor = rel_at_current((double)ref, R.mc_fit);
        R.mc_factor_err = rel_err_at_current((double)ref, R.mc_fit);

        results.push_back(R);

        std::cout << "[current_dependence] " << cfg.csv_channel
                  << " " << period
                  << " data_factor=" << R.data_factor
                  << " +/- " << R.data_factor_err
                  << " mc_factor=" << R.mc_factor
                  << " +/- " << R.mc_factor_err
                  << std::endl;
    }

    const std::string channel_dir = output_dir + "/" + cfg.output_token;
    mkdir_p(channel_dir);

    write_all_points_csv(channel_dir + "/data_current_points.csv", results, true);
    write_all_points_csv(channel_dir + "/mc_current_points.csv", results, false);

    draw_current_panel_2x3(channel_dir + "/current_dependence_data_absolute.png",
                           cfg.title + " data current dependence",
                           results,
                           true,
                           false);

    draw_current_panel_2x3(channel_dir + "/current_dependence_data_relative_to_zero.png",
                           cfg.title + " data current dependence relative to 0 nA",
                           results,
                           true,
                           true);

    if (cfg.channel == Channel::DVCS) {
        draw_current_panel_2x3(channel_dir + "/current_dependence_mc_absolute.png",
                               cfg.title + " MC current dependence",
                               results,
                               false,
                               false);

        draw_current_panel_2x3(channel_dir + "/current_dependence_mc_relative_to_zero.png",
                               cfg.title + " MC current dependence relative to 0 nA",
                               results,
                               false,
                               true);
    }

    write_summary_csv(channel_dir + "/period_summary.csv", results);

    return results;
}

static void write_override_unity(CSV& csv, const ChannelConfig& cfg) {
    for (const std::string& period : CSV_PERIOD_ORDER) {
        const int c_exp = col_strict(csv, current_eff_col(cfg, "exp", period));
        const int c_mc  = col_strict(csv, current_eff_col(cfg, "mc", period));

        for (auto& row : csv.rows) {
            row[c_exp] = "(1,0)";
            row[c_mc]  = "(1,0)";
        }
    }
}

static void write_results_to_csv(CSV& csv,
                                 const ChannelConfig& cfg,
                                 const std::vector<PeriodResult>& results) {
    std::map<std::string, PeriodResult> by_period;

    for (const PeriodResult& r : results) {
        by_period[r.period] = r;
    }

    for (const std::string& period : CSV_PERIOD_ORDER) {
        auto it = by_period.find(period);

        if (it == by_period.end()) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: missing result for period " << period
               << " channel " << cfg.csv_channel;
            fatal(ss.str());
        }

        const PeriodResult& r = it->second;

        const int c_exp = col_strict(csv, current_eff_col(cfg, "exp", period));
        const int c_mc  = col_strict(csv, current_eff_col(cfg, "mc", period));

        const std::string exp_cell = tuple2(r.data_factor, r.data_factor_err);
        const std::string mc_cell = tuple2(r.mc_factor, r.mc_factor_err);

        for (auto& row : csv.rows) {
            row[c_exp] = exp_cell;
            row[c_mc] = mc_cell;
        }
    }
}


static const PeriodResult& find_result_strict(const std::vector<PeriodResult>& results,
                                              const std::string& period,
                                              const std::string& label) {
    for (const PeriodResult& r : results) {
        if (r.period == period) {
            return r;
        }
    }

    std::ostringstream ss;
    ss << "[current_dependence] FATAL: missing " << label
       << " current-dependence result for period " << period;
    fatal(ss.str());

    return results.front();
}

static void apply_dvcs_mc_to_data_ratio_to_eppi0_results(
    const std::vector<PeriodResult>& dvcs_results,
    std::vector<PeriodResult>& eppi0_results,
    const std::string& output_dir,
    const std::string& eppi0_output_token) {

    for (PeriodResult& epi : eppi0_results) {
        const PeriodResult& dvcs = find_result_strict(dvcs_results,
                                                      epi.period,
                                                      "DVCS");

        if (!(std::isfinite(epi.data_factor) && epi.data_factor > 0.0)) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: invalid eppi0 data current factor for period "
               << epi.period << ": " << epi.data_factor;
            fatal(ss.str());
        }

        if (!(std::isfinite(dvcs.data_factor) && dvcs.data_factor > 0.0)) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: invalid DVCS data current factor for period "
               << epi.period << ": " << dvcs.data_factor;
            fatal(ss.str());
        }

        if (!(std::isfinite(dvcs.mc_factor) && dvcs.mc_factor > 0.0)) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: invalid DVCS MC current factor for period "
               << epi.period << ": " << dvcs.mc_factor;
            fatal(ss.str());
        }

        const double scale = dvcs.mc_factor / dvcs.data_factor;
        epi.mc_factor = epi.data_factor * scale;

        double rel_var = 0.0;

        if (std::isfinite(epi.data_factor_err) && epi.data_factor_err >= 0.0) {
            rel_var += (epi.data_factor_err / epi.data_factor) *
                       (epi.data_factor_err / epi.data_factor);
        }

        if (std::isfinite(dvcs.mc_factor_err) && dvcs.mc_factor_err >= 0.0) {
            rel_var += (dvcs.mc_factor_err / dvcs.mc_factor) *
                       (dvcs.mc_factor_err / dvcs.mc_factor);
        }

        if (std::isfinite(dvcs.data_factor_err) && dvcs.data_factor_err >= 0.0) {
            rel_var += (dvcs.data_factor_err / dvcs.data_factor) *
                       (dvcs.data_factor_err / dvcs.data_factor);
        }

        epi.mc_factor_err = std::fabs(epi.mc_factor) * std::sqrt(rel_var);

        epi.mc_points.clear();
        epi.mc_fit = FitResult();

        std::cout << "[current_dependence] ep->eppi0 " << epi.period
                  << " derived_mc_factor=" << epi.mc_factor
                  << " +/- " << epi.mc_factor_err
                  << " using eppi0_data_factor=" << epi.data_factor
                  << " and DVCS mc/data ratio=" << scale
                  << std::endl;
    }

    write_summary_csv(output_dir + "/" + eppi0_output_token + "/period_summary.csv",
                      eppi0_results);
}

} // namespace

bool update_current_dependence_factors_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& dvcsGenMcTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const CurrentDependenceOptions& options) {

    try {
        ROOT::EnableThreadSafety();
        TH1::AddDirectory(kFALSE);
        gStyle->SetOptStat(0);

        CSV csv;
        load_csv_strict(csv_path, csv);

        const ChannelConfig dvcs = dvcs_config();
        const ChannelConfig eppi0 = eppi0_config();

        if (options.override_to_unity) {
            std::cout << "[current_dependence] Override enabled: writing all current-efficiency factors as (1,0)." << std::endl;

            write_override_unity(csv, dvcs);
            write_override_unity(csv, eppi0);

            write_csv_atomic(csv_path, csv);

            std::cout << "[current_dependence] Updated CSV with unity current-efficiency factors: "
                      << csv_path << std::endl;

            return true;
        }

        mkdir_p(options.output_dir);

        const std::unordered_map<int, double> charge_map =
            read_charge_csv(options.charge_csv_path);

        const TopoCutMap data_cuts =
            load_sigma_cuts(options.combined_cuts_json, "data");

        const TopoCutMap mc_cuts =
            load_sigma_cuts(options.combined_cuts_json, "mc");

        std::vector<PeriodResult> dvcs_results =
            run_channel_study(dvcs,
                              dvcsDataTrees,
                              dvcsGenMcTrees,
                              dvcsRecMcTrees,
                              charge_map,
                              data_cuts,
                              mc_cuts,
                              options.output_dir,
                              options.max_workers);

        const std::map<std::string, TTree*> empty_mc_gen_trees;
        const std::map<std::string, TTree*> empty_mc_rec_trees;

        std::vector<PeriodResult> eppi0_results =
            run_channel_study(eppi0,
                              eppi0DataTrees,
                              empty_mc_gen_trees,
                              empty_mc_rec_trees,
                              charge_map,
                              data_cuts,
                              mc_cuts,
                              options.output_dir,
                              options.max_workers);

        apply_dvcs_mc_to_data_ratio_to_eppi0_results(dvcs_results,
                                                     eppi0_results,
                                                     options.output_dir,
                                                     eppi0.output_token);

        write_results_to_csv(csv, dvcs, dvcs_results);
        write_results_to_csv(csv, eppi0, eppi0_results);

        write_csv_atomic(csv_path, csv);

        std::cout << "[current_dependence] Updated current-efficiency factor columns in: "
                  << csv_path << std::endl;

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}