#include "current_dependence.h"

#include "global_cuts.h"

#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
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

static bool has_exact_current_token(const std::string& s, int current_nA) {
    std::ostringstream token_ss;
    token_ss << current_nA << "na";
    const std::string token = token_ss.str();

    size_t pos = s.find(token);

    while (pos != std::string::npos) {
        const bool left_ok =
            (pos == 0) ||
            (!std::isdigit((unsigned char)s[pos - 1]));

        const size_t right_pos = pos + token.size();
        const bool right_ok =
            (right_pos >= s.size()) ||
            (!std::isdigit((unsigned char)s[right_pos]));

        if (left_ok && right_ok) {
            return true;
        }

        pos = s.find(token, pos + 1);
    }

    return false;
}

static int parse_current_from_key(const std::string& key) {
    const std::string s = lower_ascii(key);

    if (has_substr(s, "nobkg") || has_exact_current_token(s, 0)) {
        return 0;
    }

    PeriodTags tags = parse_period_from_key(key);

    // Current-study MC convention:
    //   - nobkg or explicit 0nA files are the 0 nA reference.
    //   - Sp18 Out nonzero-current MC is the 45 nA production-current sample.
    //   - all other nonzero-current MC samples are treated as 50 nA.
    //
    // Do not scan all integers with substring matching. That misidentifies
    // 45nA and 50nA as 5nA.
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

        const std::string odir = output_dir + "/" + cfg.output_token + "/" + period_dir(period);
        mkdir_p(odir);

        write_points_csv(odir + "/data_current_points.csv", R.data_points);
        write_points_csv(odir + "/mc_current_points.csv", R.mc_points);

        draw_fit_graph(odir + "/current_dependence_absolute.png",
                       cfg.title + "  " + period,
                       R.data_points,
                       R.data_fit,
                       R.mc_points,
                       R.mc_fit,
                       false);

        draw_fit_graph(odir + "/current_dependence_relative_to_zero.png",
                       cfg.title + "  " + period,
                       R.data_points,
                       R.data_fit,
                       R.mc_points,
                       R.mc_fit,
                       true);

        std::cout << "[current_dependence] " << cfg.csv_channel
                  << " " << period
                  << " data_factor=" << R.data_factor
                  << " +/- " << R.data_factor_err
                  << " mc_factor=" << R.mc_factor
                  << " +/- " << R.mc_factor_err
                  << std::endl;
    }

    write_summary_csv(output_dir + "/" + cfg.output_token + "/period_summary.csv", results);

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

// -----------------------------------------------------------------------------
// Apply saved MC current-efficiency factors to reconstructed MC yield columns
// -----------------------------------------------------------------------------

static const std::vector<std::string>& topology_labels() {
    static const std::vector<std::string> v = {
        "(FD, FD)",
        "(CD, FD)",
        "(CD, FT)"
    };
    return v;
}

static bool parse_first_number(const std::string& cell, double& value) {
    value = std::numeric_limits<double>::quiet_NaN();

    std::string s;
    s.reserve(cell.size());
    for (char c : cell) {
        if (!std::isspace((unsigned char)c)) {
            s.push_back(c);
        }
    }

    if (s.empty()) {
        return false;
    }

    if (s.front() == '(') {
        const size_t start = 1;
        size_t stop = s.find(',', start);
        if (stop == std::string::npos) {
            stop = s.find(')', start);
        }
        if (stop == std::string::npos || stop <= start) {
            return false;
        }
        s = s.substr(start, stop - start);
    }

    char* endp = nullptr;
    const double x = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) {
        return false;
    }

    value = x;
    return std::isfinite(value);
}

static std::string format_scalar(double v) {
    if (!std::isfinite(v)) {
        return "";
    }

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << v;
    std::string out = ss.str();

    while (!out.empty() && out.back() == '0') {
        out.pop_back();
    }
    if (!out.empty() && out.back() == '.') {
        out.pop_back();
    }
    if (out.empty() || out == "-0") {
        out = "0";
    }

    return out;
}

static std::string rec_yield_total_col(const std::string& channel,
                                       const std::string& period) {
    return "reconstructed yield, " + channel + ", mc, " + period;
}

static std::string rec_yield_topo_col(const std::string& channel,
                                      const std::string& topo,
                                      const std::string& period) {
    return "reconstructed yield, " + channel + ", " + topo + ", mc, " + period;
}

static std::string rec_current_corrected_total_col(const std::string& channel,
                                                   const std::string& period) {
    return "reconstructed current corrected yield, " + channel + ", mc, " + period;
}

static std::string rec_current_corrected_topo_col(const std::string& channel,
                                                  const std::string& topo,
                                                  const std::string& period) {
    return "reconstructed current corrected yield, " + channel + ", " + topo + ", mc, " + period;
}

static double read_mc_current_factor_from_csv(const CSV& csv,
                                              const ChannelConfig& cfg,
                                              const std::string& period) {
    const int c = col_strict(csv, current_eff_col(cfg, "mc", period));

    if (csv.rows.empty()) {
        fatal("[current_dependence] FATAL: cannot read current factor from empty CSV.");
    }

    double f = std::numeric_limits<double>::quiet_NaN();
    if (!parse_first_number(csv.rows.front()[c], f) || !(std::isfinite(f) && f > 0.0)) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: invalid MC current-efficiency factor in column '"
           << current_eff_col(cfg, "mc", period) << "': '" << csv.rows.front()[c] << "'.";
        fatal(ss.str());
    }

    return f;
}

static long long apply_one_current_correction_column(CSV& csv,
                                                     const std::string& src_col,
                                                     const std::string& dst_col,
                                                     double factor) {
    const int c_src = col_strict(csv, src_col);
    const int c_dst = col_strict(csv, dst_col);

    long long n_positive = 0;

    for (auto& row : csv.rows) {
        double raw = 0.0;
        if (!parse_first_number(row[c_src], raw) || !(std::isfinite(raw) && raw > 0.0)) {
            row[c_dst].clear();
            continue;
        }

        const double corrected = raw / factor;
        row[c_dst] = format_scalar(corrected);

        if (corrected > 0.0) {
            ++n_positive;
        }
    }

    return n_positive;
}

static long long apply_mc_current_corrections_for_channel(CSV& csv,
                                                          const std::string& source_channel,
                                                          const ChannelConfig& factor_cfg,
                                                          const std::string& log_label) {
    long long total_positive = 0;

    for (const std::string& period : CSV_PERIOD_ORDER) {
        const double f_mc = read_mc_current_factor_from_csv(csv, factor_cfg, period);

        const long long n_total = apply_one_current_correction_column(
            csv,
            rec_yield_total_col(source_channel, period),
            rec_current_corrected_total_col(source_channel, period),
            f_mc);
        total_positive += n_total;

        std::cout << "[current_dependence] Applied MC current correction for "
                  << log_label << " period=" << period
                  << " inclusive factor=" << f_mc
                  << " positive_rows=" << n_total
                  << std::endl;

        for (const std::string& topo : topology_labels()) {
            const long long n_topo = apply_one_current_correction_column(
                csv,
                rec_yield_topo_col(source_channel, topo, period),
                rec_current_corrected_topo_col(source_channel, topo, period),
                f_mc);
            total_positive += n_topo;

            std::cout << "[current_dependence] Applied MC current correction for "
                      << log_label << " period=" << period
                      << " topo=" << topo
                      << " factor=" << f_mc
                      << " positive_rows=" << n_topo
                      << std::endl;
        }
    }

    return total_positive;
}

static void apply_all_mc_current_corrections(CSV& csv,
                                             const ChannelConfig& dvcs,
                                             const ChannelConfig& eppi0) {
    const long long n_dvcs = apply_mc_current_corrections_for_channel(csv,
                                                                      "ep->epg",
                                                                      dvcs,
                                                                      "ep->epg");

    const long long n_eppi0 = apply_mc_current_corrections_for_channel(csv,
                                                                       "ep->eppi0",
                                                                       eppi0,
                                                                       "ep->eppi0");

    const long long n_eppi0_bkg = apply_mc_current_corrections_for_channel(csv,
                                                                           "ep->eppi0->epg",
                                                                           eppi0,
                                                                           "ep->eppi0->epg");

    if (n_dvcs <= 0) {
        fatal("[current_dependence] FATAL: applying DVCS MC current corrections produced zero positive cells.");
    }
    if (n_eppi0 <= 0) {
        fatal("[current_dependence] FATAL: applying eppi0 MC current corrections produced zero positive cells.");
    }
    if (n_eppi0_bkg <= 0) {
        fatal("[current_dependence] FATAL: applying eppi0->epg background MC current corrections produced zero positive cells.");
    }

    std::cout << "[current_dependence] MC current-corrected reconstructed-yield cells written: "
              << "ep->epg=" << n_dvcs
              << " ep->eppi0=" << n_eppi0
              << " ep->eppi0->epg=" << n_eppi0_bkg
              << std::endl;
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
            apply_all_mc_current_corrections(csv, dvcs, eppi0);

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

        std::vector<PeriodResult> eppi0_results =
            run_channel_study(eppi0,
                              eppi0DataTrees,
                              eppi0GenMcTrees,
                              eppi0RecMcTrees,
                              charge_map,
                              data_cuts,
                              mc_cuts,
                              options.output_dir,
                              options.max_workers);

        write_results_to_csv(csv, dvcs, dvcs_results);
        write_results_to_csv(csv, eppi0, eppi0_results);
        apply_all_mc_current_corrections(csv, dvcs, eppi0);

        write_csv_atomic(csv_path, csv);

        std::cout << "[current_dependence] Updated current-efficiency factor columns in: "
                  << csv_path << std::endl;

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}