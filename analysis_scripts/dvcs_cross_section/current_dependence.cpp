// current_dependence.cpp
// -----------------------------------------------------------------------------
// Current-dependence correction study for DVCS ep -> ep gamma and
// ep -> ep pi0 channels.
//
// Main outputs:
//   1. Current-efficiency factor columns in the pass-2 analysis CSV.
//   2. Unity fallback eppi0 AAOGEN normalization columns.
//   3. Current-corrected normalized raw DATA yield columns.
//   4. Current-corrected reconstructed MC yield columns.
//   5. Per-period diagnostic plots.
//   6. 3x2 all-period summary plots with the five RGA periods plus overlay.
//   7. Cross-channel DATA slope summary plots comparing epg and eppi0.
//
// Important convention:
//   - The fitted current response is normalized to the fitted zero-current
//     intercept b.
//   - The plotted quantity in the all-period canvas is:
//
//         100 * y(I) / b
//
//     so the intercept is 100% by construction and the fitted slope shown in
//     the legend is:
//
//         100 * m / b   [% / nA]
//
// Sp19 Inb default behavior:
//   - The Sp19 Inb luminosity scan contains only one low-current point at 5 nA.
//   - That run currently has suspect Faraday Cup charge.
//   - Therefore, by default, the Sp19 Inb current-efficiency factors written
//     to the CSV are copied from Fa18 Inb.
//   - The raw Sp19 Inb scan is still processed and plotted diagnostically.
// -----------------------------------------------------------------------------

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
#include <TGraph.h>
#include <TLine.h>

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

static void warning_once(const std::string& key,
                         const std::string& message) {
    static std::mutex warn_mutex;
    static std::set<std::string> issued_warnings;

    std::lock_guard<std::mutex> lock(warn_mutex);

    if (issued_warnings.insert(key).second) {
        std::cout << message << std::endl;
    }
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

static int col_optional(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);

    if (it == csv.index.end()) {
        return -1;
    }

    return it->second;
}

static std::string tuple2(double v, double e) {
    std::ostringstream ss;
    ss << std::setprecision(12) << "(" << v << "," << e << ")";
    return ss.str();
}

static std::string tuple4(double a, double b, double c, double d) {
    std::ostringstream ss;
    ss << std::setprecision(12)
       << "(" << a << "," << b << "," << c << "," << d << ")";
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

struct ChargeEntry {
    double run_scaler = std::numeric_limits<double>::quiet_NaN();
    double hel_pos = std::numeric_limits<double>::quiet_NaN();
    double hel_neg = std::numeric_limits<double>::quiet_NaN();
    double col5 = std::numeric_limits<double>::quiet_NaN();

    bool has_run_scaler = false;
    bool has_hel_pos = false;
    bool has_hel_neg = false;
    bool has_col5 = false;
};

static double parse_optional_double(const std::vector<std::string>& fields,
                                    size_t index,
                                    bool& ok) {
    ok = false;

    if (index >= fields.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    char* endp = nullptr;
    const double v = std::strtod(fields[index].c_str(), &endp);

    if (endp == fields[index].c_str()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    ok = std::isfinite(v);
    return v;
}

static std::unordered_map<int, ChargeEntry> read_charge_csv(const std::string& path) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot open charge CSV: " << path;
        fatal(ss.str());
    }

    std::unordered_map<int, ChargeEntry> out;
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

        ChargeEntry entry;

        bool ok_run = false;
        bool ok_pos = false;
        bool ok_neg = false;
        bool ok_col5 = false;

        entry.run_scaler = parse_optional_double(fields, 1, ok_run);
        entry.hel_pos = parse_optional_double(fields, 2, ok_pos);
        entry.hel_neg = parse_optional_double(fields, 3, ok_neg);
        entry.col5 = parse_optional_double(fields, 4, ok_col5);

        entry.has_run_scaler = ok_run;
        entry.has_hel_pos = ok_pos;
        entry.has_hel_neg = ok_neg;
        entry.has_col5 = ok_col5;

        if (!entry.has_run_scaler && !(entry.has_hel_pos && entry.has_hel_neg)) {
            continue;
        }

        out[run] = entry;
    }

    std::cout << "[current_dependence] Loaded " << out.size()
              << " run-charge entries from " << path << std::endl;

    std::cout << "[current_dependence] Charge CSV interpretation:"
              << " column 1 = run number,"
              << " column 2 = RUN::Scaler-like total,"
              << " column 3 = positive-helicity scaler,"
              << " column 4 = negative-helicity scaler,"
              << " column 5 = auxiliary scaler component used only by the optional Fa18/Sp19 columns-3-to-5 mode."
              << std::endl;

    return out;
}

static bool is_spring_2018_period(const std::string& period_display) {
    return period_display == "Sp18 Inb" || period_display == "Sp18 Out";
}

static bool select_unpolarized_charge_for_period(
    const PeriodTags& tags,
    int runnum,
    const ChargeEntry& entry,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    double& charge) {

    charge = std::numeric_limits<double>::quiet_NaN();

    // Spring 2018 has zeros in the helicity/auxiliary scaler columns for this
    // purpose, so it must always use column 2 for unpolarized normalization.
    if (is_spring_2018_period(tags.display)) {
        if (!entry.has_run_scaler || !(entry.run_scaler > 0.0)) {
            std::ostringstream ss;
            ss << "[current_dependence] WARNING: missing/non-positive column-2 charge"
               << " for Spring 2018. First occurrence: run=" << runnum
               << " period=" << tags.display
               << ". Further warnings of this type for this period are suppressed.";

            warning_once("spring2018_col2_missing_or_nonpositive_" + tags.display,
                         ss.str());

            return false;
        }

        charge = entry.run_scaler;
        return true;
    }

    // Optional mode for Fa18 and Sp19 only:
    //   charge_unpol = scale * (column 3 + column 4 + column 5)
    if (use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized) {
        if (!(entry.has_hel_pos && entry.has_hel_neg && entry.has_col5)) {
            std::ostringstream ss;
            ss << "[current_dependence] WARNING: missing one or more columns 3-5"
               << " in scaled Fa18/Sp19 columns-3-to-5 charge mode."
               << " First occurrence: run=" << runnum
               << " period=" << tags.display
               << " has_col3=" << entry.has_hel_pos
               << " has_col4=" << entry.has_hel_neg
               << " has_col5=" << entry.has_col5
               << ". Further warnings of this type for this period are suppressed.";

            warning_once("fa18_sp19_scaled_cols3to5_missing_" + tags.display,
                         ss.str());

            return false;
        }

        const double raw_sum = entry.hel_pos + entry.hel_neg + entry.col5;
        charge = columns_3_to_5_charge_sum_scale * raw_sum;

        if (!(charge > 0.0)) {
            std::ostringstream ss;
            ss << "[current_dependence] WARNING: non-positive scaled Fa18/Sp19 columns-3-to-5 charge"
               << ". First occurrence: run=" << runnum
               << " period=" << tags.display
               << " col3=" << entry.hel_pos
               << " col4=" << entry.hel_neg
               << " col5=" << entry.col5
               << " scale=" << columns_3_to_5_charge_sum_scale
               << ". Further warnings of this type for this period are suppressed.";

            warning_once("fa18_sp19_scaled_cols3to5_nonpositive_" + tags.display,
                         ss.str());

            return false;
        }

        return true;
    }

    if (use_second_column_charge_for_all_unpolarized) {
        if (!entry.has_run_scaler || !(entry.run_scaler > 0.0)) {
            std::ostringstream ss;
            ss << "[current_dependence] WARNING: missing/non-positive column-2 charge"
               << ". First occurrence: run=" << runnum
               << " period=" << tags.display
               << ". Further warnings of this type for this period are suppressed.";

            warning_once("all_periods_col2_missing_or_nonpositive_" + tags.display,
                         ss.str());

            return false;
        }

        charge = entry.run_scaler;
        return true;
    }

    if (!(entry.has_hel_pos && entry.has_hel_neg)) {
        std::ostringstream ss;
        ss << "[current_dependence] WARNING: missing helicity charge columns"
           << " in legacy charge mode. First occurrence: run=" << runnum
           << " period=" << tags.display
           << " has_hel_pos=" << entry.has_hel_pos
           << " has_hel_neg=" << entry.has_hel_neg
           << ". Further warnings of this type for this period are suppressed.";

        warning_once("legacy_helicity_columns_missing_" + tags.display,
                     ss.str());

        return false;
    }

    charge = entry.hel_pos + entry.hel_neg;

    if (!(charge > 0.0)) {
        std::ostringstream ss;
        ss << "[current_dependence] WARNING: non-positive helicity-summed charge"
           << " in legacy charge mode. First occurrence: run=" << runnum
           << " period=" << tags.display
           << " hel_pos=" << entry.hel_pos
           << " hel_neg=" << entry.hel_neg
           << ". Further warnings of this type for this period are suppressed.";

        warning_once("legacy_helicity_sum_nonpositive_" + tags.display,
                     ss.str());

        return false;
    }

    return true;
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

            vm[vit.key()] = s;
        }

        if (!vm.empty()) {
            out[key] = vm;
        }
    }

    return out;
}

static bool within_cut_window(double v, const SigmaStats& s) {
    if (!std::isfinite(v)) {
        return false;
    }

    if (s.mode == "upper_quantile") {
        if (!std::isfinite(s.cut_high)) return true;
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

    return within_cut_window(value, it->second);
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

        ena("x");
        ena("Q2");
        ena("phi2");

        ena("t1");
        ena("open_angle_ep2");
        ena("pTmiss");
        ena("Delta_phi");

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
        bD("theta_pi0_pi0", &theta_pi0_pi0, has_theta_pi0_pi0);

        bD("e_p", &e_p, has_e_p);
        bD("e_theta", &e_theta, has_e_theta);
        bD("e_phi", &e_phi, has_e_phi);

        bD("p1_theta", &p1_theta, has_p1_theta);
        bD("p1_phi", &p1_phi, has_p1_phi);

        bD("p2_p", &p2_p, has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi", &p2_phi, has_p2_phi);
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
    const GlobalCutConfig& cfg = default_global_cuts();

    if (!(b.has_t1 && b.has_open_angle_ep2)) return false;
    if (cfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (b.has_runnum && is_excluded_run(b.runnum)) return false;

    if (!(b.has_detector1 && b.has_detector2)) {
        fatal("[current_dependence] FATAL: missing detector1/detector2.");
    }

    if (global_cuts_require_sector_phi(cfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            fatal("[current_dependence] FATAL: sector selection requires e_phi, p1_phi, p2_phi.");
        }
    }

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[current_dependence] FATAL: missing branches required by global ycol cuts.");
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
    bool has_delta_phi = false;
    const double delta_phi = b.delta_phi_value(has_delta_phi);
    if (!check_sigma_var(vm, "Delta_phi", has_delta_phi, delta_phi)) return false;
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

    // Used only by the kinematic current-efficiency diagnostics.
    // For the ordinary integrated current-dependence fits this remains zero.
    double kinematic_value_sum = 0.0;
    long long kinematic_value_count = 0;
};

static DataAgg process_data_tree(
    const ChannelConfig& cfg,
    const std::string& key,
    TTree* tree,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale) {

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

        double charge = std::numeric_limits<double>::quiet_NaN();

        if (!select_unpolarized_charge_for_period(tags,
                                                  b.runnum,
                                                  charge_it->second,
                                                  use_second_column_charge_for_all_unpolarized,
                                                  use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                                                  columns_3_to_5_charge_sum_scale,
                                                  charge)) {
            continue;
        }

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
// Plot helpers
// -----------------------------------------------------------------------------

static int period_color(const std::string& period) {
    if (period == "Sp18 Inb") return kAzure + 1;
    if (period == "Sp18 Out") return kOrange + 7;
    if (period == "Fa18 Inb") return kGreen + 2;
    if (period == "Fa18 Out") return kRed + 1;
    if (period == "Sp19 Inb") return kViolet + 1;
    return kBlack;
}

static double fit_percent_slope(const FitResult& f) {
    if (!std::isfinite(f.m) || !std::isfinite(f.b) || f.b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return 100.0 * f.m / f.b;
}

static double fit_percent_slope_err(const FitResult& f) {
    if (!std::isfinite(f.m) ||
        !std::isfinite(f.b) ||
        !std::isfinite(f.sm) ||
        !std::isfinite(f.sb) ||
        !std::isfinite(f.cov_mb) ||
        f.b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double dm = 100.0 / f.b;
    const double db = -100.0 * f.m / (f.b * f.b);

    double var =
        dm * dm * f.sm * f.sm +
        db * db * f.sb * f.sb +
        2.0 * dm * db * f.cov_mb;

    if (var < 0.0 && std::fabs(var) < 1.0e-15) {
        var = 0.0;
    }

    if (var < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(var);
}

static double percent_response_at_current(double current, const FitResult& f) {
    if (!std::isfinite(f.m) || !std::isfinite(f.b) || f.b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return 100.0 * (f.m * current + f.b) / f.b;
}

static double percent_response_err_at_current(double current, const FitResult& f) {
    if (!std::isfinite(f.m) ||
        !std::isfinite(f.b) ||
        !std::isfinite(f.sm) ||
        !std::isfinite(f.sb) ||
        !std::isfinite(f.cov_mb) ||
        f.b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double dm = 100.0 * current / f.b;
    const double db = -100.0 * f.m * current / (f.b * f.b);

    double var =
        dm * dm * f.sm * f.sm +
        db * db * f.sb * f.sb +
        2.0 * dm * db * f.cov_mb;

    if (var < 0.0 && std::fabs(var) < 1.0e-15) {
        var = 0.0;
    }

    if (var < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(var);
}

static TGraphErrors* make_percent_points_graph(const std::vector<CurrentPoint>& points,
                                               const FitResult& fit,
                                               int color) {
    TGraphErrors* g = new TGraphErrors();

    int ip = 0;

    for (const CurrentPoint& p : points) {
        if (!std::isfinite(fit.b) || fit.b == 0.0) {
            continue;
        }

        const double y = 100.0 * p.y / fit.b;
        const double ey = 100.0 * p.sy / std::fabs(fit.b);

        if (!std::isfinite(y) || !std::isfinite(ey)) {
            continue;
        }

        g->SetPoint(ip, (double)p.current_nA, y);
        g->SetPointError(ip, 0.0, ey);
        ++ip;
    }

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.25);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);

    return g;
}

static TGraph* make_percent_fit_line(const FitResult& fit,
                                     int color,
                                     double xmin,
                                     double xmax) {
    TGraph* g = new TGraph();

    int ip = 0;

    for (int i = 0; i < 200; ++i) {
        const double x = xmin + (xmax - xmin) * (double)i / 199.0;
        const double y = percent_response_at_current(x, fit);

        if (!std::isfinite(y)) {
            continue;
        }

        g->SetPoint(ip, x, y);
        ++ip;
    }

    g->SetLineColor(color);
    g->SetLineWidth(2);

    return g;
}

static TGraphErrors* make_percent_fit_band(const FitResult& fit,
                                           int color,
                                           double xmin,
                                           double xmax) {
    TGraphErrors* g = new TGraphErrors();

    int ip = 0;

    for (int i = 0; i < 200; ++i) {
        const double x = xmin + (xmax - xmin) * (double)i / 199.0;
        const double y = percent_response_at_current(x, fit);
        const double ey = percent_response_err_at_current(x, fit);

        if (!std::isfinite(y) || !std::isfinite(ey)) {
            continue;
        }

        g->SetPoint(ip, x, y);
        g->SetPointError(ip, 0.0, ey);
        ++ip;
    }

    g->SetFillColorAlpha(color, 0.18);
    g->SetLineColor(color);
    g->SetLineWidth(1);

    return g;
}

static void style_current_pad() {
    gPad->SetGrid(1, 1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.08);
}

static void draw_current_frame(const std::string& title) {
    TH1D* frame = (TH1D*)gPad->DrawFrame(0.0, 0.0, 100.0, 150.0);
    frame->SetTitle(title.c_str());
    frame->GetXaxis()->SetTitle("Beam current (nA)");
    frame->GetYaxis()->SetTitle("Percent of intercept b (%)");

    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
}

static bool hide_sp19_inb_from_replacement_plots(bool use_fa18_for_sp19,
                                                    const std::string& period) {
    return use_fa18_for_sp19 && period == "Sp19 Inb";
}

static void draw_replacement_hidden_panel(const std::string& period,
                                          bool suppress_axis_titles = false) {
    draw_current_frame(period);

    TH1* frame = gPad ? dynamic_cast<TH1*>(gPad->GetPrimitive("hframe")) : nullptr;

    if (suppress_axis_titles && frame) {
        frame->GetYaxis()->SetTitle("");
        frame->GetXaxis()->SetTitle("");
    }

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.035);
    lat.SetTextColor(kRed + 2);
    lat.DrawLatex(0.18, 0.58, "Not plotted");
    lat.DrawLatex(0.18, 0.50, "CSV current-efficiency");
    lat.DrawLatex(0.18, 0.42, "factor copied from Fa18 Inb");
}

static void draw_period_panel(const PeriodResult& r,
                              bool use_replacement_annotation) {
    const int color = period_color(r.period);

    draw_current_frame(r.period);

    TGraphErrors* band = make_percent_fit_band(r.data_fit, color, 0.0, 100.0);
    TGraph* line = make_percent_fit_line(r.data_fit, color, 0.0, 100.0);
    TGraphErrors* points = make_percent_points_graph(r.data_points, r.data_fit, color);

    if (band && band->GetN() > 0) {
        band->Draw("3 SAME");
    }

    if (line && line->GetN() > 0) {
        line->Draw("L SAME");
    }

    if (points && points->GetN() > 0) {
        points->Draw("PE SAME");
    }

    const double slope = fit_percent_slope(r.data_fit);
    const double slope_err = fit_percent_slope_err(r.data_fit);

    TLegend* leg = new TLegend(0.40, 0.82, 0.94, 0.96);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.032);

    std::ostringstream l1;
    l1 << "b=" << std::fixed << std::setprecision(5) << r.data_fit.b
       << " +/- " << std::setprecision(5) << r.data_fit.sb;

    std::ostringstream l2;
    l2 << "slope=" << std::fixed << std::setprecision(5) << slope
       << " +/- " << std::setprecision(5) << slope_err << " (%/nA)";

    leg->AddEntry(points, l1.str().c_str(), "pe");
    leg->AddEntry((TObject*)nullptr, l2.str().c_str(), "");
    leg->Draw();

    if (use_replacement_annotation && r.period == "Sp19 Inb") {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.032);
        lat.SetTextColor(kRed + 2);
        lat.DrawLatex(0.16, 0.76, "CSV value replaced by Fa18 Inb");
    }
}

static void draw_summary_panel(const std::vector<PeriodResult>& results,
                               bool use_fa18_for_sp19) {
    draw_current_frame("All periods (overlay)");

    TLegend* leg = new TLegend(0.34, 0.67, 0.96, 0.96);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.028);

    for (const std::string& period : PERIOD_ORDER) {
        auto it = std::find_if(results.begin(), results.end(),
                               [&](const PeriodResult& r) {
                                   return r.period == period;
                               });

        if (it == results.end()) {
            continue;
        }

        if (hide_sp19_inb_from_replacement_plots(use_fa18_for_sp19, period)) {
            continue;
        }

        const PeriodResult& r = *it;
        const int color = period_color(r.period);

        TGraphErrors* band = make_percent_fit_band(r.data_fit, color, 0.0, 100.0);
        TGraph* line = make_percent_fit_line(r.data_fit, color, 0.0, 100.0);
        TGraphErrors* points = make_percent_points_graph(r.data_points, r.data_fit, color);

        if (band && band->GetN() > 0) {
            band->Draw("3 SAME");
        }

        if (line && line->GetN() > 0) {
            line->Draw("L SAME");
        }

        if (points && points->GetN() > 0) {
            points->Draw("PE SAME");
        }

        const double slope = fit_percent_slope(r.data_fit);
        const double slope_err = fit_percent_slope_err(r.data_fit);

        std::ostringstream label;
        label << r.period << ": slope="
              << std::fixed << std::setprecision(5) << slope
              << " +/- " << std::setprecision(5) << slope_err
              << " (%/nA)";

        if (use_fa18_for_sp19 && r.period == "Sp19 Inb") {
            label << " [diagnostic]";
        }

        leg->AddEntry(points, label.str().c_str(), "pe");
    }

    leg->Draw();

    if (use_fa18_for_sp19) {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.029);
        lat.SetTextColor(kRed + 2);
        lat.DrawLatex(0.15, 0.61, "Sp19 Inb CSV factor copied from Fa18 Inb");
    }
}

static void draw_all_period_current_canvas(const std::string& out_path,
                                           const std::string& canvas_title,
                                           const std::vector<PeriodResult>& results,
                                           bool use_fa18_for_sp19) {
    const size_t slash_pos = out_path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(out_path.substr(0, slash_pos));
    }

    TCanvas c("c_all_period_current_dependence", canvas_title.c_str(), 1800, 1000);
    c.Divide(3, 2, 0.001, 0.001);

    int pad = 1;

    for (const std::string& period : PERIOD_ORDER) {
        c.cd(pad);
        style_current_pad();

        auto it = std::find_if(results.begin(), results.end(),
                               [&](const PeriodResult& r) {
                                   return r.period == period;
                               });

        if (hide_sp19_inb_from_replacement_plots(use_fa18_for_sp19, period)) {
            draw_replacement_hidden_panel(period);
        } else if (it != results.end()) {
            draw_period_panel(*it, use_fa18_for_sp19);
        } else {
            draw_current_frame(period);
        }

        ++pad;
    }

    c.cd(6);
    style_current_pad();
    draw_summary_panel(results, use_fa18_for_sp19);

    c.SaveAs(out_path.c_str());
}

static double mean_or_bin_center(double xlo,
                                 double xhi,
                                 double sum,
                                 long long count) {
    const double bin_center = 0.5 * (xlo + xhi);

    if (count <= 0) {
        return bin_center;
    }

    const double mean = sum / (double)count;

    if (!std::isfinite(mean)) {
        return bin_center;
    }

    // The value should normally lie inside the bin by construction.
    // Keep the fallback conservative in case of a malformed branch value.
    if (mean < xlo || mean > xhi) {
        return bin_center;
    }

    return mean;
}


// -----------------------------------------------------------------------------
// Data/MC overlay 3x2 diagnostic canvas
// -----------------------------------------------------------------------------

static TGraphErrors* make_percent_points_graph_style(const std::vector<CurrentPoint>& points,
                                                     const FitResult& fit,
                                                     int color,
                                                     int marker_style,
                                                     double marker_size) {
    TGraphErrors* g = new TGraphErrors();

    int ip = 0;

    for (const CurrentPoint& p : points) {
        if (!std::isfinite(fit.b) || fit.b == 0.0) {
            continue;
        }

        const double y = 100.0 * p.y / fit.b;
        const double ey = 100.0 * p.sy / std::fabs(fit.b);

        if (!std::isfinite(y) || !std::isfinite(ey)) {
            continue;
        }

        g->SetPoint(ip, (double)p.current_nA, y);
        g->SetPointError(ip, 0.0, ey);
        ++ip;
    }

    g->SetMarkerStyle(marker_style);
    g->SetMarkerSize(marker_size);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);

    return g;
}

static TGraph* make_percent_fit_line_style(const FitResult& fit,
                                           int color,
                                           int line_style,
                                           double xmin,
                                           double xmax) {
    TGraph* g = make_percent_fit_line(fit, color, xmin, xmax);
    g->SetLineStyle(line_style);
    g->SetLineWidth(2);
    return g;
}

static TGraphErrors* make_ratio_points_graph(const PeriodResult& r,
                                             int color) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;

    for (const CurrentPoint& p : r.data_points) {
        const double y_data = percent_response_at_current((double)p.current_nA, r.data_fit);
        const double y_mc = percent_response_at_current((double)p.current_nA, r.mc_fit);

        if (!std::isfinite(y_data) || !std::isfinite(y_mc) || y_mc == 0.0 ||
            !std::isfinite(r.data_fit.b) || r.data_fit.b == 0.0) {
            continue;
        }

        const double point_data = 100.0 * p.y / r.data_fit.b;
        const double point_data_err = 100.0 * p.sy / std::fabs(r.data_fit.b);

        if (!std::isfinite(point_data) || !std::isfinite(point_data_err)) {
            continue;
        }

        const double ratio = point_data / y_mc;
        const double ratio_err = point_data_err / std::fabs(y_mc);

        if (!std::isfinite(ratio) || !std::isfinite(ratio_err)) {
            continue;
        }

        g->SetPoint(ip, (double)p.current_nA, ratio);
        g->SetPointError(ip, 0.0, ratio_err);
        ++ip;
    }

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.05);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);

    return g;
}

static TGraph* make_ratio_fit_line(const PeriodResult& r,
                                   int color,
                                   double xmin,
                                   double xmax) {
    TGraph* g = new TGraph();
    int ip = 0;

    for (int i = 0; i < 200; ++i) {
        const double x = xmin + (xmax - xmin) * (double)i / 199.0;
        const double y_data = percent_response_at_current(x, r.data_fit);
        const double y_mc = percent_response_at_current(x, r.mc_fit);

        if (!std::isfinite(y_data) || !std::isfinite(y_mc) || y_mc == 0.0) {
            continue;
        }

        g->SetPoint(ip, x, y_data / y_mc);
        ++ip;
    }

    g->SetLineColor(color);
    g->SetLineWidth(2);
    return g;
}

static bool has_valid_mc_fit_or_points(const PeriodResult& r) {
    if (!r.mc_points.empty()) {
        return true;
    }

    return std::isfinite(r.mc_fit.m) && std::isfinite(r.mc_fit.b) && r.mc_fit.b != 0.0;
}

static void style_data_mc_subpad(bool top_pad) {
    gPad->SetGrid(1, 1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.04);

    if (top_pad) {
        gPad->SetBottomMargin(0.02);
        gPad->SetTopMargin(0.10);
    } else {
        gPad->SetBottomMargin(0.28);
        gPad->SetTopMargin(0.02);
    }
}

static void draw_data_mc_period_stack(const PeriodResult* rptr,
                                      const std::string& title,
                                      bool summary_panel,
                                      const std::vector<PeriodResult>& all_results,
                                      bool use_fa18_for_sp19) {
    TPad* top = new TPad((title + "_top").c_str(), "", 0.0, 0.36, 1.0, 1.0);
    TPad* bot = new TPad((title + "_bot").c_str(), "", 0.0, 0.0, 1.0, 0.36);
    top->Draw();
    bot->Draw();

    top->cd();
    style_data_mc_subpad(true);

    TH1D* frame_top = (TH1D*)gPad->DrawFrame(0.0, 40.0, 80.0, 120.0);
    frame_top->SetTitle(title.c_str());
    frame_top->GetXaxis()->SetTitle("");
    frame_top->GetXaxis()->SetLabelSize(0.0);
    frame_top->GetYaxis()->SetTitle("Efficiency relative to fitted 0 nA (%)");
    frame_top->GetYaxis()->CenterTitle(true);
    frame_top->GetYaxis()->SetTitleSize(0.060);
    frame_top->GetYaxis()->SetLabelSize(0.052);
    frame_top->GetYaxis()->SetTitleOffset(0.80);

    TLegend* leg_top = nullptr;

    if (summary_panel) {
        leg_top = new TLegend(0.08, 0.06, 0.48, 0.55);
        leg_top->SetTextSize(0.030);
    } else {
        leg_top = new TLegend(0.68, 0.74, 0.95, 0.95);
        leg_top->SetTextSize(0.040);
    }

    leg_top->SetFillStyle(1001);
    leg_top->SetFillColor(kWhite);
    leg_top->SetBorderSize(1);

    auto draw_top_one = [&](const PeriodResult& r, bool add_simple_legend) {
        const int color = period_color(r.period);

        TGraphErrors* band_data = make_percent_fit_band(r.data_fit, color, 0.0, 80.0);
        TGraph* line_data = make_percent_fit_line_style(r.data_fit, color, 1, 0.0, 80.0);
        TGraphErrors* pts_data = make_percent_points_graph_style(r.data_points, r.data_fit, color, 20, summary_panel ? 0.80 : 1.20);

        if (band_data && band_data->GetN() > 0 && !summary_panel) {
            band_data->Draw("3 SAME");
        }
        if (line_data && line_data->GetN() > 0) {
            line_data->Draw("L SAME");
        }
        if (pts_data && pts_data->GetN() > 0) {
            pts_data->Draw("PE SAME");
        }

        if (has_valid_mc_fit_or_points(r)) {
            TGraph* line_mc = make_percent_fit_line_style(r.mc_fit, color, 2, 0.0, 80.0);
            TGraphErrors* pts_mc = make_percent_points_graph_style(r.mc_points, r.mc_fit, color, 24, summary_panel ? 0.80 : 1.20);

            if (line_mc && line_mc->GetN() > 0) {
                line_mc->Draw("L SAME");
            }
            if (pts_mc && pts_mc->GetN() > 0) {
                pts_mc->Draw("PE SAME");
            }

            if (add_simple_legend) {
                leg_top->AddEntry(pts_data, "Data", "pe");
                leg_top->AddEntry(pts_mc, "MC", "pe");
            } else if (summary_panel) {
                leg_top->AddEntry(pts_data, r.period.c_str(), "pe");
            }
        } else if (add_simple_legend) {
            leg_top->AddEntry(pts_data, "Data", "pe");
        }
    };

    if (summary_panel) {
        for (const std::string& period : PERIOD_ORDER) {
            if (hide_sp19_inb_from_replacement_plots(use_fa18_for_sp19, period)) {
                continue;
            }

            auto it = std::find_if(all_results.begin(), all_results.end(),
                                   [&](const PeriodResult& r) { return r.period == period; });
            if (it != all_results.end()) {
                draw_top_one(*it, false);
            }
        }
    } else if (rptr) {
        draw_top_one(*rptr, true);

        TLatex ref;
        ref.SetNDC(true);
        ref.SetTextSize(0.038);
        ref.SetTextColor(kBlack);
        std::ostringstream ss;
        ss << "MC ref in acceptance: " << reference_current_nA(rptr->period) << " nA";
        ref.DrawLatex(0.18, 0.08, ss.str().c_str());

        if (use_fa18_for_sp19 && rptr->period == "Sp19 Inb") {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.035);
            lat.SetTextColor(kRed + 2);
            lat.DrawLatex(0.08, 0.17, "CSV value replaced by Fa18 Inb");
        }
    }

    leg_top->Draw();

    bot->cd();
    style_data_mc_subpad(false);

    TH1D* frame_bot = (TH1D*)gPad->DrawFrame(0.0, 0.70, 80.0, 1.10);
    frame_bot->SetTitle("");
    frame_bot->GetXaxis()->SetTitle("Beam current (nA)");
    frame_bot->GetYaxis()->SetTitle("Data/MC");
    frame_bot->GetXaxis()->CenterTitle(true);
    frame_bot->GetYaxis()->CenterTitle(true);
    frame_bot->GetXaxis()->SetTitleSize(0.095);
    frame_bot->GetYaxis()->SetTitleSize(0.090);
    frame_bot->GetXaxis()->SetLabelSize(0.080);
    frame_bot->GetYaxis()->SetLabelSize(0.080);
    frame_bot->GetYaxis()->SetTitleOffset(0.55);

    TLine unity(0.0, 1.0, 80.0, 1.0);
    unity.SetLineColor(kGray + 2);
    unity.SetLineStyle(2);
    unity.SetLineWidth(2);
    unity.Draw("SAME");

    auto draw_ratio_one = [&](const PeriodResult& r) {
        if (!has_valid_mc_fit_or_points(r)) {
            return;
        }

        const int color = period_color(r.period);
        TGraph* line = make_ratio_fit_line(r, color, 0.0, 80.0);
        TGraphErrors* pts = make_ratio_points_graph(r, color);

        if (line && line->GetN() > 0) {
            line->Draw("L SAME");
        }
        if (pts && pts->GetN() > 0) {
            pts->Draw("PE SAME");
        }
    };

    if (summary_panel) {
        for (const std::string& period : PERIOD_ORDER) {
            if (hide_sp19_inb_from_replacement_plots(use_fa18_for_sp19, period)) {
                continue;
            }

            auto it = std::find_if(all_results.begin(), all_results.end(),
                                   [&](const PeriodResult& r) { return r.period == period; });
            if (it != all_results.end()) {
                draw_ratio_one(*it);
            }
        }
    } else if (rptr) {
        draw_ratio_one(*rptr);
    }
}

static void draw_all_period_data_mc_canvas(const std::string& out_path,
                                           const std::string& canvas_title,
                                           const std::vector<PeriodResult>& results,
                                           bool use_fa18_for_sp19) {
    const bool have_any_mc = std::any_of(results.begin(), results.end(),
                                         [](const PeriodResult& r) {
                                             return has_valid_mc_fit_or_points(r);
                                         });

    if (!have_any_mc) {
        std::cout << "[current_dependence] Skipping data/MC overlay canvas "
                  << out_path << " because no MC current-dependence points are available."
                  << std::endl;
        return;
    }

    const size_t slash_pos = out_path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(out_path.substr(0, slash_pos));
    }

    TCanvas c("c_all_period_data_mc_current_dependence", canvas_title.c_str(), 1800, 1200);
    c.Divide(3, 2, 0.001, 0.001);

    int pad = 1;

    for (const std::string& period : PERIOD_ORDER) {
        c.cd(pad);

        auto it = std::find_if(results.begin(), results.end(),
                               [&](const PeriodResult& r) {
                                   return r.period == period;
                               });

        if (hide_sp19_inb_from_replacement_plots(use_fa18_for_sp19, period)) {
            draw_replacement_hidden_panel(period, true);
        } else if (it != results.end()) {
            draw_data_mc_period_stack(&(*it), period, false, results, use_fa18_for_sp19);
        } else {
            draw_data_mc_period_stack(nullptr, period, false, results, use_fa18_for_sp19);
        }

        ++pad;
    }

    c.cd(6);
    draw_data_mc_period_stack(nullptr, "All periods", true, results, use_fa18_for_sp19);

    c.SaveAs(out_path.c_str());
}


// -----------------------------------------------------------------------------
// Cross-channel DATA slope summary plots: epg versus eppi0
// -----------------------------------------------------------------------------

static const PeriodResult* find_period_result(const std::vector<PeriodResult>& results,
                                              const std::string& period) {
    auto it = std::find_if(results.begin(), results.end(),
                           [&](const PeriodResult& r) {
                               return r.period == period;
                           });

    if (it == results.end()) {
        return nullptr;
    }

    return &(*it);
}

static bool finite_slope_value(const PeriodResult* r) {
    if (!r) {
        return false;
    }

    const double slope = fit_percent_slope(r->data_fit);
    const double slope_err = fit_percent_slope_err(r->data_fit);

    return std::isfinite(slope) && std::isfinite(slope_err);
}

static void write_data_channel_slope_comparison_csv(
    const std::string& path,
    const std::vector<PeriodResult>& dvcs_results,
    const std::vector<PeriodResult>& eppi0_results) {

    const size_t slash_pos = path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(path.substr(0, slash_pos));
    }

    std::ofstream fout(path);

    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot write: " << path;
        fatal(ss.str());
    }

    fout << "period,"
         << "epg_data_slope_percent_per_nA,epg_data_slope_percent_per_nA_stat,"
         << "eppi0_data_slope_percent_per_nA,eppi0_data_slope_percent_per_nA_stat,"
         << "eppi0_minus_epg_slope_percent_per_nA,eppi0_minus_epg_slope_percent_per_nA_stat,"
         << "epg_data_fit_intercept,epg_data_fit_intercept_stat,"
         << "eppi0_data_fit_intercept,eppi0_data_fit_intercept_stat\n";

    for (const std::string& period : PERIOD_ORDER) {
        const PeriodResult* epg = find_period_result(dvcs_results, period);
        const PeriodResult* pi0 = find_period_result(eppi0_results, period);

        const double epg_slope = epg ? fit_percent_slope(epg->data_fit) : std::numeric_limits<double>::quiet_NaN();
        const double epg_err = epg ? fit_percent_slope_err(epg->data_fit) : std::numeric_limits<double>::quiet_NaN();
        const double pi0_slope = pi0 ? fit_percent_slope(pi0->data_fit) : std::numeric_limits<double>::quiet_NaN();
        const double pi0_err = pi0 ? fit_percent_slope_err(pi0->data_fit) : std::numeric_limits<double>::quiet_NaN();

        double diff = std::numeric_limits<double>::quiet_NaN();
        double diff_err = std::numeric_limits<double>::quiet_NaN();

        if (std::isfinite(epg_slope) && std::isfinite(pi0_slope)) {
            diff = pi0_slope - epg_slope;
        }

        if (std::isfinite(epg_err) && std::isfinite(pi0_err)) {
            diff_err = std::sqrt(epg_err * epg_err + pi0_err * pi0_err);
        }

        fout << period << ","
             << std::setprecision(12) << epg_slope << ","
             << std::setprecision(12) << epg_err << ","
             << std::setprecision(12) << pi0_slope << ","
             << std::setprecision(12) << pi0_err << ","
             << std::setprecision(12) << diff << ","
             << std::setprecision(12) << diff_err << ","
             << std::setprecision(12) << (epg ? epg->data_fit.b : std::numeric_limits<double>::quiet_NaN()) << ","
             << std::setprecision(12) << (epg ? epg->data_fit.sb : std::numeric_limits<double>::quiet_NaN()) << ","
             << std::setprecision(12) << (pi0 ? pi0->data_fit.b : std::numeric_limits<double>::quiet_NaN()) << ","
             << std::setprecision(12) << (pi0 ? pi0->data_fit.sb : std::numeric_limits<double>::quiet_NaN()) << "\n";
    }
}

static TGraphErrors* make_channel_slope_graph(const std::vector<PeriodResult>& results,
                                               double x_offset,
                                               int color,
                                               int marker_style) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;

    for (int i = 0; i < (int)PERIOD_ORDER.size(); ++i) {
        const PeriodResult* r = find_period_result(results, PERIOD_ORDER[i]);

        if (!finite_slope_value(r)) {
            continue;
        }

        const double slope = fit_percent_slope(r->data_fit);
        const double slope_err = fit_percent_slope_err(r->data_fit);

        g->SetPoint(ip, (double)(i + 1) + x_offset, slope);
        g->SetPointError(ip, 0.0, slope_err);
        ++ip;
    }

    g->SetMarkerStyle(marker_style);
    g->SetMarkerSize(1.45);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);

    return g;
}

static TGraphErrors* make_channel_slope_difference_graph(
    const std::vector<PeriodResult>& dvcs_results,
    const std::vector<PeriodResult>& eppi0_results,
    int color,
    int marker_style) {

    TGraphErrors* g = new TGraphErrors();
    int ip = 0;

    for (int i = 0; i < (int)PERIOD_ORDER.size(); ++i) {
        const std::string& period = PERIOD_ORDER[i];
        const PeriodResult* epg = find_period_result(dvcs_results, period);
        const PeriodResult* pi0 = find_period_result(eppi0_results, period);

        if (!finite_slope_value(epg) || !finite_slope_value(pi0)) {
            continue;
        }

        const double epg_slope = fit_percent_slope(epg->data_fit);
        const double epg_err = fit_percent_slope_err(epg->data_fit);
        const double pi0_slope = fit_percent_slope(pi0->data_fit);
        const double pi0_err = fit_percent_slope_err(pi0->data_fit);

        const double diff = pi0_slope - epg_slope;
        const double diff_err = std::sqrt(epg_err * epg_err + pi0_err * pi0_err);

        g->SetPoint(ip, (double)(i + 1), diff);
        g->SetPointError(ip, 0.0, diff_err);
        ++ip;
    }

    g->SetMarkerStyle(marker_style);
    g->SetMarkerSize(1.45);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);

    return g;
}

static void style_period_axis_labels(TAxis* axis) {
    if (!axis) {
        return;
    }

    for (int i = 0; i < (int)PERIOD_ORDER.size(); ++i) {
        axis->SetBinLabel(i + 1, PERIOD_ORDER[i].c_str());
    }

    axis->LabelsOption("h");
    axis->CenterTitle(true);
}


static std::pair<double, double> data_channel_slope_range(
    const std::vector<PeriodResult>& dvcs_results,
    const std::vector<PeriodResult>& eppi0_results) {

    std::vector<double> values;

    for (const std::string& period : PERIOD_ORDER) {
        const PeriodResult* epg = find_period_result(dvcs_results, period);
        const PeriodResult* pi0 = find_period_result(eppi0_results, period);

        if (finite_slope_value(epg)) {
            const double y = fit_percent_slope(epg->data_fit);
            const double ey = fit_percent_slope_err(epg->data_fit);
            values.push_back(y - ey);
            values.push_back(y + ey);
        }

        if (finite_slope_value(pi0)) {
            const double y = fit_percent_slope(pi0->data_fit);
            const double ey = fit_percent_slope_err(pi0->data_fit);
            values.push_back(y - ey);
            values.push_back(y + ey);
        }
    }

    if (values.empty()) {
        return std::make_pair(-1.0, 1.0);
    }

    double ymin = *std::min_element(values.begin(), values.end());
    double ymax = *std::max_element(values.begin(), values.end());

    ymin = std::min(ymin, 0.0);
    ymax = std::max(ymax, 0.0);

    double span = ymax - ymin;
    if (!(span > 0.0)) {
        span = 1.0;
    }

    ymin -= 0.18 * span;
    ymax += 0.18 * span;

    return std::make_pair(ymin, ymax);
}

static std::pair<double, double> data_channel_slope_difference_range(
    const std::vector<PeriodResult>& dvcs_results,
    const std::vector<PeriodResult>& eppi0_results) {

    std::vector<double> values;

    for (const std::string& period : PERIOD_ORDER) {
        const PeriodResult* epg = find_period_result(dvcs_results, period);
        const PeriodResult* pi0 = find_period_result(eppi0_results, period);

        if (!finite_slope_value(epg) || !finite_slope_value(pi0)) {
            continue;
        }

        const double epg_slope = fit_percent_slope(epg->data_fit);
        const double epg_err = fit_percent_slope_err(epg->data_fit);
        const double pi0_slope = fit_percent_slope(pi0->data_fit);
        const double pi0_err = fit_percent_slope_err(pi0->data_fit);
        const double diff = pi0_slope - epg_slope;
        const double diff_err = std::sqrt(epg_err * epg_err + pi0_err * pi0_err);

        values.push_back(diff - diff_err);
        values.push_back(diff + diff_err);
    }

    if (values.empty()) {
        return std::make_pair(-1.0, 1.0);
    }

    double ymin = *std::min_element(values.begin(), values.end());
    double ymax = *std::max_element(values.begin(), values.end());

    ymin = std::min(ymin, 0.0);
    ymax = std::max(ymax, 0.0);

    double span = ymax - ymin;
    if (!(span > 0.0)) {
        span = 1.0;
    }

    ymin -= 0.18 * span;
    ymax += 0.18 * span;

    return std::make_pair(ymin, ymax);
}

static void draw_data_channel_slope_comparison_canvas(
    const std::string& out_path,
    const std::vector<PeriodResult>& dvcs_results,
    const std::vector<PeriodResult>& eppi0_results,
    bool use_fa18_for_sp19) {

    const size_t slash_pos = out_path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(out_path.substr(0, slash_pos));
    }

    TCanvas c("c_data_channel_slope_comparison", "DATA current slopes: epg versus eppi0", 1500, 850);
    c.SetGrid(1, 1);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.03);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);

    const std::pair<double, double> yrange = data_channel_slope_range(dvcs_results, eppi0_results);

    TH1D frame("h_data_channel_slope_comparison_frame", "DATA current-dependence slopes", 5, 0.5, 5.5);
    frame.SetMinimum(yrange.first);
    frame.SetMaximum(yrange.second);
    frame.GetXaxis()->SetTitle("Run period");
    frame.GetYaxis()->SetTitle("DATA slope (%/nA)");
    frame.GetXaxis()->SetTitleSize(0.045);
    frame.GetYaxis()->SetTitleSize(0.045);
    frame.GetXaxis()->SetLabelSize(0.040);
    frame.GetYaxis()->SetLabelSize(0.040);
    frame.GetXaxis()->CenterTitle(true);
    frame.GetYaxis()->CenterTitle(true);
    style_period_axis_labels(frame.GetXaxis());
    frame.Draw("AXIS");

    TLine zero(0.5, 0.0, 5.5, 0.0);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw("SAME");

    TGraphErrors* g_epg = make_channel_slope_graph(dvcs_results, -0.10, kBlue + 1, 20);
    TGraphErrors* g_pi0 = make_channel_slope_graph(eppi0_results, +0.10, kRed + 1, 21);

    if (g_epg && g_epg->GetN() > 0) {
        g_epg->Draw("PE SAME");
    }

    if (g_pi0 && g_pi0->GetN() > 0) {
        g_pi0->Draw("PE SAME");
    }

    TLegend* leg = new TLegend(0.64, 0.75, 0.96, 0.92);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.036);
    leg->AddEntry(g_epg, "DATA ep #rightarrow ep#gamma", "pe");
    leg->AddEntry(g_pi0, "DATA ep #rightarrow ep#pi^{0}", "pe");
    leg->Draw();

    if (use_fa18_for_sp19) {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.030);
        lat.SetTextColor(kRed + 2);
        lat.DrawLatex(0.15, 0.84, "Sp19 Inb point shows direct diagnostic fit; CSV factor is copied from Fa18 Inb.");
    }

    c.SaveAs(out_path.c_str());
}

static void draw_data_channel_slope_difference_canvas(
    const std::string& out_path,
    const std::vector<PeriodResult>& dvcs_results,
    const std::vector<PeriodResult>& eppi0_results,
    bool use_fa18_for_sp19) {

    const size_t slash_pos = out_path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(out_path.substr(0, slash_pos));
    }

    TCanvas c("c_data_channel_slope_difference", "DATA current slope difference: eppi0 minus epg", 1500, 850);
    c.SetGrid(1, 1);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.03);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);

    const std::pair<double, double> yrange = data_channel_slope_difference_range(dvcs_results, eppi0_results);

    TH1D frame("h_data_channel_slope_difference_frame", "DATA current-dependence slope difference", 5, 0.5, 5.5);
    frame.SetMinimum(yrange.first);
    frame.SetMaximum(yrange.second);
    frame.GetXaxis()->SetTitle("Run period");
    frame.GetYaxis()->SetTitle("e#pi^{0} slope - ep#gamma slope (%/nA)");
    frame.GetXaxis()->SetTitleSize(0.045);
    frame.GetYaxis()->SetTitleSize(0.045);
    frame.GetXaxis()->SetLabelSize(0.040);
    frame.GetYaxis()->SetLabelSize(0.040);
    frame.GetXaxis()->CenterTitle(true);
    frame.GetYaxis()->CenterTitle(true);
    style_period_axis_labels(frame.GetXaxis());
    frame.Draw("AXIS");

    TLine zero(0.5, 0.0, 5.5, 0.0);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw("SAME");

    TGraphErrors* g_diff = make_channel_slope_difference_graph(dvcs_results, eppi0_results, kBlack, 20);

    if (g_diff && g_diff->GetN() > 0) {
        g_diff->Draw("PE SAME");
    }

    TLegend* leg = new TLegend(0.59, 0.80, 0.96, 0.92);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.036);
    leg->AddEntry(g_diff, "DATA ep#pi^{0} - DATA ep#gamma", "pe");
    leg->Draw();

    if (use_fa18_for_sp19) {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.030);
        lat.SetTextColor(kRed + 2);
        lat.DrawLatex(0.15, 0.84, "Sp19 Inb point shows direct diagnostic fit; CSV factor is copied from Fa18 Inb.");
    }

    c.SaveAs(out_path.c_str());
}

static void draw_data_channel_current_response_comparison_canvas(
    const std::string& out_path,
    const std::vector<PeriodResult>& dvcs_results,
    const std::vector<PeriodResult>& eppi0_results,
    bool use_fa18_for_sp19) {

    const size_t slash_pos = out_path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(out_path.substr(0, slash_pos));
    }

    TCanvas c("c_data_channel_current_response_comparison", "DATA current response: epg versus eppi0", 1800, 1000);
    c.Divide(3, 2, 0.001, 0.001);

    int pad = 1;

    for (const std::string& period : PERIOD_ORDER) {
        c.cd(pad);
        style_current_pad();
        draw_current_frame(period);

        const PeriodResult* epg = find_period_result(dvcs_results, period);
        const PeriodResult* pi0 = find_period_result(eppi0_results, period);

        TLegend* leg = new TLegend(0.36, 0.72, 0.96, 0.96);
        leg->SetFillStyle(1001);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.030);

        if (epg && std::isfinite(epg->data_fit.b) && epg->data_fit.b != 0.0) {
            TGraphErrors* band = make_percent_fit_band(epg->data_fit, kBlue + 1, 0.0, 100.0);
            TGraph* line = make_percent_fit_line_style(epg->data_fit, kBlue + 1, 1, 0.0, 100.0);
            TGraphErrors* points = make_percent_points_graph_style(epg->data_points, epg->data_fit, kBlue + 1, 20, 1.05);

            if (band && band->GetN() > 0) band->Draw("3 SAME");
            if (line && line->GetN() > 0) line->Draw("L SAME");
            if (points && points->GetN() > 0) points->Draw("PE SAME");

            std::ostringstream ss;
            ss << "ep#gamma slope=" << std::fixed << std::setprecision(4)
               << fit_percent_slope(epg->data_fit) << " +/- "
               << fit_percent_slope_err(epg->data_fit);
            leg->AddEntry(points, ss.str().c_str(), "pe");
        }

        if (pi0 && std::isfinite(pi0->data_fit.b) && pi0->data_fit.b != 0.0) {
            TGraphErrors* band = make_percent_fit_band(pi0->data_fit, kRed + 1, 0.0, 100.0);
            TGraph* line = make_percent_fit_line_style(pi0->data_fit, kRed + 1, 2, 0.0, 100.0);
            TGraphErrors* points = make_percent_points_graph_style(pi0->data_points, pi0->data_fit, kRed + 1, 21, 1.05);

            if (band && band->GetN() > 0) band->Draw("3 SAME");
            if (line && line->GetN() > 0) line->Draw("L SAME");
            if (points && points->GetN() > 0) points->Draw("PE SAME");

            std::ostringstream ss;
            ss << "e#pi^{0} slope=" << std::fixed << std::setprecision(4)
               << fit_percent_slope(pi0->data_fit) << " +/- "
               << fit_percent_slope_err(pi0->data_fit);
            leg->AddEntry(points, ss.str().c_str(), "pe");
        }

        leg->Draw();

        if (use_fa18_for_sp19 && period == "Sp19 Inb") {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.028);
            lat.SetTextColor(kRed + 2);
            lat.DrawLatex(0.16, 0.64, "Direct diagnostic fit shown");
            lat.DrawLatex(0.16, 0.58, "CSV factor copied from Fa18 Inb");
        }

        ++pad;
    }

    c.cd(6);
    gPad->SetGrid(1, 1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.08);

    const std::pair<double, double> yrange = data_channel_slope_range(dvcs_results, eppi0_results);

    TH1D frame("h_data_channel_current_response_summary_frame", "DATA slope summary", 5, 0.5, 5.5);
    frame.SetMinimum(yrange.first);
    frame.SetMaximum(yrange.second);
    frame.GetXaxis()->SetTitle("Run period");
    frame.GetYaxis()->SetTitle("DATA slope (%/nA)");
    frame.GetXaxis()->SetTitleSize(0.045);
    frame.GetYaxis()->SetTitleSize(0.045);
    frame.GetXaxis()->SetLabelSize(0.036);
    frame.GetYaxis()->SetLabelSize(0.040);
    frame.GetXaxis()->CenterTitle(true);
    frame.GetYaxis()->CenterTitle(true);
    style_period_axis_labels(frame.GetXaxis());
    frame.Draw("AXIS");

    TLine zero(0.5, 0.0, 5.5, 0.0);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw("SAME");

    TGraphErrors* g_epg = make_channel_slope_graph(dvcs_results, -0.10, kBlue + 1, 20);
    TGraphErrors* g_pi0 = make_channel_slope_graph(eppi0_results, +0.10, kRed + 1, 21);

    if (g_epg && g_epg->GetN() > 0) g_epg->Draw("PE SAME");
    if (g_pi0 && g_pi0->GetN() > 0) g_pi0->Draw("PE SAME");

    TLegend* leg = new TLegend(0.50, 0.76, 0.96, 0.93);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.034);
    leg->AddEntry(g_epg, "ep#gamma", "pe");
    leg->AddEntry(g_pi0, "e#pi^{0}", "pe");
    leg->Draw();

    c.SaveAs(out_path.c_str());
}

static void draw_data_channel_slope_summary_outputs(
    const std::string& output_dir,
    const std::vector<PeriodResult>& dvcs_results,
    const std::vector<PeriodResult>& eppi0_results,
    bool use_fa18_for_sp19) {

    const std::string dir = output_dir + "/summary";
    mkdir_p(dir);

    draw_data_channel_slope_comparison_canvas(dir + "/data_epg_eppi0_slope_comparison.png",
                                              dvcs_results,
                                              eppi0_results,
                                              use_fa18_for_sp19);

    draw_data_channel_slope_difference_canvas(dir + "/data_eppi0_minus_epg_slope_difference.png",
                                              dvcs_results,
                                              eppi0_results,
                                              use_fa18_for_sp19);

    draw_data_channel_current_response_comparison_canvas(dir + "/data_epg_eppi0_current_response_overlay.png",
                                                         dvcs_results,
                                                         eppi0_results,
                                                         use_fa18_for_sp19);

    write_data_channel_slope_comparison_csv(dir + "/data_epg_eppi0_slope_comparison.csv",
                                            dvcs_results,
                                            eppi0_results);
}

// -----------------------------------------------------------------------------
// Existing per-period plots
// -----------------------------------------------------------------------------

static void draw_fit_graph(const std::string& out_path,
                           const std::string& title,
                           const std::vector<CurrentPoint>& data_points,
                           const FitResult& data_fit,
                           const std::vector<CurrentPoint>& mc_points,
                           const FitResult& mc_fit,
                           bool relative) {
    const size_t slash_pos = out_path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(out_path.substr(0, slash_pos));
    }

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
        frame->GetYaxis()->SetTitle("Data counts / charge or MC efficiency");
    }

    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);

    g_data.Draw("PE SAME");

    const bool have_mc_points = !mc_points.empty();

    if (have_mc_points) {
        g_mc.Draw("PE SAME");
    }

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

    if (have_mc_points) {
        draw_line(mc_fit, kRed + 1, relative);
    }

    TLegend leg(0.58, 0.72, 0.90, 0.88);
    leg.SetFillStyle(1001);
    leg.SetFillColor(kWhite);
    leg.SetBorderSize(1);
    leg.AddEntry(&g_data, "Data", "pe");

    if (have_mc_points) {
        leg.AddEntry(&g_mc, "MC", "pe");
    }

    leg.Draw();

    c.SaveAs(out_path.c_str());
}

static void write_points_csv(const std::string& path,
                             const std::vector<CurrentPoint>& points) {
    const size_t slash_pos = path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(path.substr(0, slash_pos));
    }

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
    const size_t slash_pos = path.find_last_of('/');

    if (slash_pos != std::string::npos) {
        mkdir_p(path.substr(0, slash_pos));
    }

    std::ofstream fout(path);

    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: cannot write: " << path;
        fatal(ss.str());
    }

    fout << "period,data_factor,data_factor_stat,mc_factor,mc_factor_stat,"
         << "data_m,data_b,data_sm,data_sb,data_cov_mb,data_slope_percent_per_nA,data_slope_percent_per_nA_stat,"
         << "mc_m,mc_b,mc_sm,mc_sb,mc_cov_mb\n";

    for (const PeriodResult& r : rows) {
        fout << r.period << ","
             << std::setprecision(12) << r.data_factor << ","
             << std::setprecision(12) << r.data_factor_err << ","
             << std::setprecision(12) << r.mc_factor << ","
             << std::setprecision(12) << r.mc_factor_err << ","
             << std::setprecision(12) << r.data_fit.m << ","
             << std::setprecision(12) << r.data_fit.b << ","
             << std::setprecision(12) << r.data_fit.sm << ","
             << std::setprecision(12) << r.data_fit.sb << ","
             << std::setprecision(12) << r.data_fit.cov_mb << ","
             << std::setprecision(12) << fit_percent_slope(r.data_fit) << ","
             << std::setprecision(12) << fit_percent_slope_err(r.data_fit) << ","
             << std::setprecision(12) << r.mc_fit.m << ","
             << std::setprecision(12) << r.mc_fit.b << ","
             << std::setprecision(12) << r.mc_fit.sm << ","
             << std::setprecision(12) << r.mc_fit.sb << ","
             << std::setprecision(12) << r.mc_fit.cov_mb << "\n";
    }
}

// -----------------------------------------------------------------------------
// Sp19 replacement helpers
// -----------------------------------------------------------------------------

static void copy_period_result_values(PeriodResult& dst,
                                      const PeriodResult& src,
                                      const std::string& reason) {
    dst.data_factor = src.data_factor;
    dst.data_factor_err = src.data_factor_err;
    dst.mc_factor = src.mc_factor;
    dst.mc_factor_err = src.mc_factor_err;

    std::cout << "[current_dependence] " << reason
              << ": copied current-efficiency factors from "
              << src.period << " to " << dst.period
              << " data_factor=" << dst.data_factor
              << " +/- " << dst.data_factor_err
              << " mc_factor=" << dst.mc_factor
              << " +/- " << dst.mc_factor_err
              << std::endl;
}

static void replace_sp19_inb_factors_with_fa18_inb(std::vector<PeriodResult>& results,
                                                   const std::string& channel_label) {
    auto it_fa18 = std::find_if(results.begin(), results.end(),
                                [](const PeriodResult& r) {
                                    return r.period == "Fa18 Inb";
                                });

    auto it_sp19 = std::find_if(results.begin(), results.end(),
                                [](const PeriodResult& r) {
                                    return r.period == "Sp19 Inb";
                                });

    if (it_fa18 == results.end()) {
        fatal("[current_dependence] FATAL: cannot copy Fa18 Inb factor to Sp19 Inb; missing Fa18 Inb result.");
    }

    if (it_sp19 == results.end()) {
        fatal("[current_dependence] FATAL: cannot copy Fa18 Inb factor to Sp19 Inb; missing Sp19 Inb result.");
    }

    copy_period_result_values(*it_sp19,
                              *it_fa18,
                              "channel=" + channel_label + " Sp19 Inb replacement");
}

// -----------------------------------------------------------------------------
// Main channel worker
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Kinematic current-efficiency diagnostic canvas
// -----------------------------------------------------------------------------

struct KinematicVarConfig {
    std::string key;
    std::string title;
    std::string x_label;
    std::vector<double> edges;
};

struct KinematicBinResult {
    double x_low = 0.0;
    double x_high = 0.0;
    double x_center = 0.0;
    double x_err = 0.0;
    double factor = std::numeric_limits<double>::quiet_NaN();
    double factor_err = std::numeric_limits<double>::quiet_NaN();
    double total_counts = 0.0;
    int n_current_points = 0;
};

static std::vector<KinematicVarConfig> kinematic_current_var_configs() {
    return {
        {"Q2",      "Q^{2}",          "Q^{2} (GeV^{2})",          {1.0, 1.4, 1.8, 2.3, 3.0, 4.0, 5.2, 6.5}},
        {"xB",      "x_{B}",          "x_{B}",                    {0.06, 0.10, 0.14, 0.18, 0.25, 0.35, 0.45, 0.60}},
        {"t",       "|t|",            "|t| (GeV^{2})",            {0.05, 0.15, 0.25, 0.35, 0.50, 0.70, 0.90, 1.25}},
        {"phi",     "#phi",           "#phi (deg)",               {0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 360.0}},
        {"e_theta", "#theta_{e}",     "#theta_{e} (deg)",         {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0}},
        {"p_theta", "#theta_{p}",     "#theta_{p} (deg)",         {10.0, 18.0, 24.0, 30.0, 38.0, 46.0, 54.0, 62.0}},
        {"g_theta", "#theta_{#gamma}","#theta_{#gamma} (deg)",    {2.0, 6.0, 10.0, 14.0, 18.0, 22.0, 26.0, 30.0}}
    };
}

static bool kinematic_value_for_config(const Branches& b,
                                       const KinematicVarConfig& cfg,
                                       double& value) {
    if (cfg.key == "Q2") {
        if (!b.has_Q2) return false;
        value = b.Q2;
        return std::isfinite(value);
    }

    if (cfg.key == "xB") {
        if (!b.has_x) return false;
        value = b.x;
        return std::isfinite(value);
    }

    if (cfg.key == "t") {
        if (!b.has_t1) return false;
        value = std::fabs(b.t1);
        return std::isfinite(value);
    }

    if (cfg.key == "phi") {
        if (!b.has_phi2) return false;
        value = std::fmod(b.phi2 * RAD2DEG, 360.0);
        if (value < 0.0) value += 360.0;
        return std::isfinite(value);
    }

    if (cfg.key == "e_theta") {
        if (!b.has_e_theta) return false;
        value = b.e_theta * RAD2DEG;
        return std::isfinite(value);
    }

    if (cfg.key == "p_theta") {
        if (!b.has_p1_theta) return false;
        value = b.p1_theta * RAD2DEG;
        return std::isfinite(value);
    }

    if (cfg.key == "g_theta") {
        if (!b.has_p2_theta) return false;
        value = b.p2_theta * RAD2DEG;
        return std::isfinite(value);
    }

    return false;
}

static int find_bin_index(const std::vector<double>& edges, double value) {
    if (edges.size() < 2 || !std::isfinite(value)) {
        return -1;
    }

    if (value < edges.front() || value > edges.back()) {
        return -1;
    }

    if (value == edges.back()) {
        return (int)edges.size() - 2;
    }

    auto it = std::upper_bound(edges.begin(), edges.end(), value);
    int idx = (int)std::distance(edges.begin(), it) - 1;

    if (idx < 0 || idx >= (int)edges.size() - 1) {
        return -1;
    }

    return idx;
}

using KinematicAggMap = std::map<std::string, std::vector<DataAgg>>;
using PeriodKinematicAggMap = std::map<std::string, KinematicAggMap>;

static KinematicAggMap process_data_tree_for_kinematic_current_diagnostics(
    const ChannelConfig& cfg,
    const std::string& key,
    TTree* tree,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    const std::vector<KinematicVarConfig>& vars,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale) {

    KinematicAggMap out;
    PeriodTags tags = parse_period_from_key(key);

    for (const KinematicVarConfig& v : vars) {
        std::vector<DataAgg>& bins = out[v.key];
        bins.resize(v.edges.size() > 1 ? v.edges.size() - 1 : 0);

        for (DataAgg& a : bins) {
            a.period = tags.display;
        }
    }

    if (!tree || tags.display == "Fa18 Inb Supp") {
        return out;
    }

    Branches b;
    b.bind(tree);

    if (!b.has_runnum) {
        fatal("[current_dependence] FATAL: kinematic diagnostic data tree missing runnum.");
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

        double charge = std::numeric_limits<double>::quiet_NaN();

        if (!select_unpolarized_charge_for_period(tags,
                                                  b.runnum,
                                                  charge_it->second,
                                                  use_second_column_charge_for_all_unpolarized,
                                                  use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                                                  columns_3_to_5_charge_sum_scale,
                                                  charge)) {
            continue;
        }

        if (!(charge > 0.0)) {
            continue;
        }

        if (!passes_global_dispatch(b, tags)) {
            continue;
        }

        if (!passes_sigma_dispatch(cfg, tags, data_cuts, b)) {
            continue;
        }

        for (const KinematicVarConfig& v : vars) {
            double value = std::numeric_limits<double>::quiet_NaN();

            if (!kinematic_value_for_config(b, v, value)) {
                continue;
            }

            const int ibin = find_bin_index(v.edges, value);

            if (ibin < 0) {
                continue;
            }

            DataAgg& agg = out[v.key][ibin];
            agg.counts_by_run[b.runnum] += 1;
            agg.current_by_run[b.runnum] = current;
            agg.charge_by_run[b.runnum] = charge;
            agg.kinematic_value_sum += value;
            agg.kinematic_value_count += 1;
        }
    }

    return out;
}

static void merge_data_agg(DataAgg& dst, const DataAgg& src) {
    if (dst.period.empty()) {
        dst.period = src.period;
    }

    for (const auto& kv : src.counts_by_run) {
        dst.counts_by_run[kv.first] += kv.second;
    }

    for (const auto& kv : src.current_by_run) {
        dst.current_by_run[kv.first] = kv.second;
    }

    for (const auto& kv : src.charge_by_run) {
        dst.charge_by_run[kv.first] = kv.second;
    }

    dst.kinematic_value_sum += src.kinematic_value_sum;
    dst.kinematic_value_count += src.kinematic_value_count;
}

static PeriodKinematicAggMap build_kinematic_current_diagnostic_aggs(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& data_trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    const std::vector<KinematicVarConfig>& vars,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    int max_workers) {

    PeriodKinematicAggMap out;

    for (const std::string& period : PERIOD_ORDER) {
        for (const KinematicVarConfig& v : vars) {
            std::vector<DataAgg>& bins = out[period][v.key];
            bins.resize(v.edges.size() > 1 ? v.edges.size() - 1 : 0);

            for (DataAgg& a : bins) {
                a.period = period;
            }
        }
    }

    std::vector<std::pair<std::string, TTree*>> items;
    items.reserve(data_trees.size());

    for (const auto& kv : data_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);

        if (tags.display == "Fa18 Inb Supp") {
            continue;
        }

        items.push_back(kv);
    }

    const int nth = std::max(1, std::min(max_workers, std::max(1, (int)items.size())));
    std::cout << "[current_dependence] Kinematic DATA diagnostics: processing "
              << items.size() << " trees with " << nth << " worker(s)." << std::endl;
    std::mutex out_mutex;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int i = 0; i < (int)items.size(); ++i) {
        const auto& kv = items[i];
        PeriodTags tags = parse_period_from_key(kv.first);

        KinematicAggMap one = process_data_tree_for_kinematic_current_diagnostics(
            cfg,
            kv.first,
            kv.second,
            charge_map,
            data_cuts,
            vars,
            use_second_column_charge_for_all_unpolarized,
            use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
            columns_3_to_5_charge_sum_scale);

        std::lock_guard<std::mutex> lock(out_mutex);

        for (const KinematicVarConfig& v : vars) {
            auto it = one.find(v.key);

            if (it == one.end()) {
                continue;
            }

            std::vector<DataAgg>& dst_bins = out[tags.display][v.key];
            const std::vector<DataAgg>& src_bins = it->second;
            const size_t n = std::min(dst_bins.size(), src_bins.size());

            for (size_t ib = 0; ib < n; ++ib) {
                merge_data_agg(dst_bins[ib], src_bins[ib]);
            }
        }
    }

    return out;
}


struct McKinematicBinAgg {
    std::string period;
    std::map<int, McAgg> by_current;

    // Mean x-position for MC kinematic diagnostic points.
    // Prefer the generated denominator distribution, since the MC current
    // factor is reconstructed/generated. Fall back to reconstructed if the
    // generated tree is unavailable.
    double generated_kinematic_value_sum = 0.0;
    long long generated_kinematic_value_count = 0;

    double reconstructed_kinematic_value_sum = 0.0;
    long long reconstructed_kinematic_value_count = 0;
};

using KinematicMcAggMap = std::map<std::string, std::vector<McKinematicBinAgg>>;
using PeriodKinematicMcAggMap = std::map<std::string, KinematicMcAggMap>;

static KinematicMcAggMap process_generated_tree_for_kinematic_current_diagnostics(
    const ChannelConfig& cfg,
    const std::string& key,
    TTree* tree,
    const std::vector<KinematicVarConfig>& vars) {

    (void)cfg;
    KinematicMcAggMap out;
    PeriodTags tags = parse_period_from_key(key);
    const int current = parse_current_from_key(key);

    for (const KinematicVarConfig& v : vars) {
        std::vector<McKinematicBinAgg>& bins = out[v.key];
        bins.resize(v.edges.size() > 1 ? v.edges.size() - 1 : 0);
        for (McKinematicBinAgg& a : bins) {
            a.period = tags.display;
        }
    }

    if (!tree || tags.display == "Fa18 Inb Supp") {
        return out;
    }

    Branches b;
    b.bind(tree);

    const Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        for (const KinematicVarConfig& v : vars) {
            double value = std::numeric_limits<double>::quiet_NaN();
            if (!kinematic_value_for_config(b, v, value)) {
                continue;
            }

            const int ibin = find_bin_index(v.edges, value);
            if (ibin < 0) {
                continue;
            }

            McKinematicBinAgg& bin = out[v.key][ibin];
            McAgg& a = bin.by_current[current];
            a.period = tags.display;
            a.current_nA = current;
            a.n_gen += 1;

            bin.generated_kinematic_value_sum += value;
            bin.generated_kinematic_value_count += 1;
        }
    }

    return out;
}

static KinematicMcAggMap process_reconstructed_tree_for_kinematic_current_diagnostics(
    const ChannelConfig& cfg,
    const std::string& key,
    TTree* tree,
    const TopoCutMap& mc_cuts,
    const std::vector<KinematicVarConfig>& vars) {

    KinematicMcAggMap out;
    PeriodTags tags = parse_period_from_key(key);
    const int current = parse_current_from_key(key);

    for (const KinematicVarConfig& v : vars) {
        std::vector<McKinematicBinAgg>& bins = out[v.key];
        bins.resize(v.edges.size() > 1 ? v.edges.size() - 1 : 0);
        for (McKinematicBinAgg& a : bins) {
            a.period = tags.display;
        }
    }

    if (!tree || tags.display == "Fa18 Inb Supp") {
        return out;
    }

    Branches b;
    b.bind(tree);

    const Long64_t N = tree->GetEntries();
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

        for (const KinematicVarConfig& v : vars) {
            double value = std::numeric_limits<double>::quiet_NaN();
            if (!kinematic_value_for_config(b, v, value)) {
                continue;
            }

            const int ibin = find_bin_index(v.edges, value);
            if (ibin < 0) {
                continue;
            }

            McKinematicBinAgg& bin = out[v.key][ibin];
            McAgg& a = bin.by_current[current];
            a.period = tags.display;
            a.current_nA = current;
            a.n_rec += 1;

            bin.reconstructed_kinematic_value_sum += value;
            bin.reconstructed_kinematic_value_count += 1;
        }
    }

    return out;
}

static void merge_mc_kinematic_bin(McKinematicBinAgg& dst, const McKinematicBinAgg& src) {
    if (dst.period.empty()) {
        dst.period = src.period;
    }

    for (const auto& kv : src.by_current) {
        const int current = kv.first;
        const McAgg& s = kv.second;
        McAgg& d = dst.by_current[current];
        d.period = s.period;
        d.current_nA = current;
        d.n_gen += s.n_gen;
        d.n_rec += s.n_rec;
    }

    dst.generated_kinematic_value_sum += src.generated_kinematic_value_sum;
    dst.generated_kinematic_value_count += src.generated_kinematic_value_count;

    dst.reconstructed_kinematic_value_sum += src.reconstructed_kinematic_value_sum;
    dst.reconstructed_kinematic_value_count += src.reconstructed_kinematic_value_count;
}

static PeriodKinematicMcAggMap build_kinematic_current_diagnostic_mc_aggs(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& gen_trees,
    const std::map<std::string, TTree*>& rec_trees,
    const TopoCutMap& mc_cuts,
    const std::vector<KinematicVarConfig>& vars,
    int max_workers) {

    PeriodKinematicMcAggMap out;

    for (const std::string& period : PERIOD_ORDER) {
        for (const KinematicVarConfig& v : vars) {
            std::vector<McKinematicBinAgg>& bins = out[period][v.key];
            bins.resize(v.edges.size() > 1 ? v.edges.size() - 1 : 0);
            for (McKinematicBinAgg& a : bins) {
                a.period = period;
            }
        }
    }

    auto merge_one_locked = [&](const KinematicMcAggMap& one, const std::string& period) {
        for (const KinematicVarConfig& v : vars) {
            auto it = one.find(v.key);
            if (it == one.end()) {
                continue;
            }
            std::vector<McKinematicBinAgg>& dst_bins = out[period][v.key];
            const std::vector<McKinematicBinAgg>& src_bins = it->second;
            const size_t n = std::min(dst_bins.size(), src_bins.size());
            for (size_t ib = 0; ib < n; ++ib) {
                merge_mc_kinematic_bin(dst_bins[ib], src_bins[ib]);
            }
        }
    };

    std::mutex out_mutex;

    std::vector<std::pair<std::string, TTree*>> gen_items;
    gen_items.reserve(gen_trees.size());
    for (const auto& kv : gen_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);
        if (tags.display == "Fa18 Inb Supp") {
            continue;
        }
        gen_items.push_back(kv);
    }

    int nth = std::max(1, std::min(max_workers, std::max(1, (int)gen_items.size())));
    std::cout << "[current_dependence] Kinematic MC diagnostics: processing "
              << gen_items.size() << " generated trees with " << nth << " worker(s)." << std::endl;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int i = 0; i < (int)gen_items.size(); ++i) {
        const auto& kv = gen_items[i];
        PeriodTags tags = parse_period_from_key(kv.first);
        KinematicMcAggMap one = process_generated_tree_for_kinematic_current_diagnostics(
            cfg, kv.first, kv.second, vars);

        std::lock_guard<std::mutex> lock(out_mutex);
        merge_one_locked(one, tags.display);
    }

    std::vector<std::pair<std::string, TTree*>> rec_items;
    rec_items.reserve(rec_trees.size());
    for (const auto& kv : rec_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);
        if (tags.display == "Fa18 Inb Supp") {
            continue;
        }
        rec_items.push_back(kv);
    }

    nth = std::max(1, std::min(max_workers, std::max(1, (int)rec_items.size())));
    std::cout << "[current_dependence] Kinematic MC diagnostics: processing "
              << rec_items.size() << " reconstructed trees with " << nth << " worker(s)." << std::endl;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int i = 0; i < (int)rec_items.size(); ++i) {
        const auto& kv = rec_items[i];
        PeriodTags tags = parse_period_from_key(kv.first);
        KinematicMcAggMap one = process_reconstructed_tree_for_kinematic_current_diagnostics(
            cfg, kv.first, kv.second, mc_cuts, vars);

        std::lock_guard<std::mutex> lock(out_mutex);
        merge_one_locked(one, tags.display);
    }

    return out;
}

static KinematicBinResult current_factor_for_kinematic_mc_bin(const McKinematicBinAgg& agg,
                                                              double xlo,
                                                              double xhi) {
    KinematicBinResult r;
    r.x_low = xlo;
    r.x_high = xhi;

    if (agg.generated_kinematic_value_count > 0) {
        r.x_center = mean_or_bin_center(xlo,
                                        xhi,
                                        agg.generated_kinematic_value_sum,
                                        agg.generated_kinematic_value_count);
    } else {
        r.x_center = mean_or_bin_center(xlo,
                                        xhi,
                                        agg.reconstructed_kinematic_value_sum,
                                        agg.reconstructed_kinematic_value_count);
    }

    r.x_err = 0.5 * (xhi - xlo);

    std::vector<McAgg> aggs;
    for (const auto& kv : agg.by_current) {
        aggs.push_back(kv.second);
    }

    std::vector<CurrentPoint> pts = mc_points_from_aggs(aggs, agg.period);
    r.n_current_points = (int)pts.size();

    for (const CurrentPoint& p : pts) {
        r.total_counts += p.counts;
    }

    if (pts.empty()) {
        return r;
    }

    FitResult fit = fit_points(pts);
    const int ref = reference_current_nA(agg.period);
    r.factor = rel_at_current((double)ref, fit);
    r.factor_err = rel_err_at_current((double)ref, fit);

    return r;
}

static KinematicBinResult current_factor_for_kinematic_bin(const DataAgg& agg,
                                                           double xlo,
                                                           double xhi) {
    KinematicBinResult r;
    r.x_low = xlo;
    r.x_high = xhi;
    r.x_center = mean_or_bin_center(xlo,
                                    xhi,
                                    agg.kinematic_value_sum,
                                    agg.kinematic_value_count);
    r.x_err = 0.5 * (xhi - xlo);

    std::vector<CurrentPoint> pts = data_points_from_agg(agg);
    r.n_current_points = (int)pts.size();

    for (const CurrentPoint& p : pts) {
        r.total_counts += p.counts;
    }

    if (pts.empty()) {
        return r;
    }

    FitResult fit = fit_points(pts);
    r.factor = weighted_data_rel(pts, fit);
    r.factor_err = weighted_data_rel_err(pts, fit);

    return r;
}

static TGraphErrors* make_kinematic_factor_graph(const std::vector<KinematicBinResult>& bins,
                                                 int color) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;

    for (const KinematicBinResult& b : bins) {
        if (!std::isfinite(b.factor) || !std::isfinite(b.factor_err) || b.total_counts <= 0.0) {
            continue;
        }

        g->SetPoint(ip, b.x_center, b.factor);
        // Deliberately suppress horizontal error bars in this diagnostic plot.
        g->SetPointError(ip, 0.0, b.factor_err);
        ++ip;
    }

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.05);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);

    return g;
}


struct LinearFitSummary {
    bool valid = false;
    double m = std::numeric_limits<double>::quiet_NaN();
    double b = std::numeric_limits<double>::quiet_NaN();
    double m_err = std::numeric_limits<double>::quiet_NaN();
    double b_err = std::numeric_limits<double>::quiet_NaN();
    double cov_mb = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    int npoints = 0;
};

static LinearFitSummary fit_kinematic_factor_linear(const std::vector<KinematicBinResult>& bins) {
    LinearFitSummary out;

    double S = 0.0;
    double Sx = 0.0;
    double Sy = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;

    for (const KinematicBinResult& p : bins) {
        if (!std::isfinite(p.x_center) ||
            !std::isfinite(p.factor) ||
            !std::isfinite(p.factor_err) ||
            !(p.factor_err > 0.0) ||
            p.total_counts <= 0.0) {
            continue;
        }

        const double w = 1.0 / (p.factor_err * p.factor_err);
        S += w;
        Sx += w * p.x_center;
        Sy += w * p.factor;
        Sxx += w * p.x_center * p.x_center;
        Sxy += w * p.x_center * p.factor;
        out.npoints += 1;
    }

    if (out.npoints < 2) {
        return out;
    }

    const double det = S * Sxx - Sx * Sx;
    if (!(det > 0.0) || !std::isfinite(det)) {
        return out;
    }

    out.b = (Sxx * Sy - Sx * Sxy) / det;
    out.m = (S * Sxy - Sx * Sy) / det;
    out.b_err = std::sqrt(Sxx / det);
    out.m_err = std::sqrt(S / det);
    out.cov_mb = -Sx / det;
    out.ndf = out.npoints - 2;

    double chi2 = 0.0;
    for (const KinematicBinResult& p : bins) {
        if (!std::isfinite(p.x_center) ||
            !std::isfinite(p.factor) ||
            !std::isfinite(p.factor_err) ||
            !(p.factor_err > 0.0) ||
            p.total_counts <= 0.0) {
            continue;
        }

        const double yfit = out.m * p.x_center + out.b;
        const double pull = (p.factor - yfit) / p.factor_err;
        chi2 += pull * pull;
    }

    out.chi2 = chi2;
    out.valid = std::isfinite(out.m) && std::isfinite(out.b) &&
                std::isfinite(out.m_err) && std::isfinite(out.b_err) &&
                std::isfinite(out.cov_mb);
    return out;
}

static double linear_fit_prediction_uncertainty(const LinearFitSummary& fit, double x) {
    if (!fit.valid) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double var = fit.b_err * fit.b_err +
                       x * x * fit.m_err * fit.m_err +
                       2.0 * x * fit.cov_mb;

    if (!(var >= 0.0) || !std::isfinite(var)) {
        return 0.0;
    }

    return std::sqrt(var);
}

static TGraph* make_linear_fit_line_graph(const LinearFitSummary& fit,
                                          double xmin,
                                          double xmax,
                                          int color) {
    TGraph* g = new TGraph();
    if (!fit.valid) {
        return g;
    }

    const int n = 120;
    for (int i = 0; i < n; ++i) {
        const double f = (n == 1) ? 0.0 : double(i) / double(n - 1);
        const double x = xmin + f * (xmax - xmin);
        const double y = fit.m * x + fit.b;
        g->SetPoint(i, x, y);
    }

    g->SetLineColor(color);
    g->SetLineWidth(2);
    return g;
}

static TGraph* make_linear_fit_band_graph(const LinearFitSummary& fit,
                                          double xmin,
                                          double xmax,
                                          int color) {
    TGraph* g = new TGraph();
    if (!fit.valid) {
        return g;
    }

    const int n = 120;
    for (int i = 0; i < n; ++i) {
        const double f = (n == 1) ? 0.0 : double(i) / double(n - 1);
        const double x = xmin + f * (xmax - xmin);
        const double y = fit.m * x + fit.b;
        const double ey = linear_fit_prediction_uncertainty(fit, x);
        g->SetPoint(i, x, y + ey);
    }

    for (int i = 0; i < n; ++i) {
        const double f = (n == 1) ? 0.0 : double(n - 1 - i) / double(n - 1);
        const double x = xmin + f * (xmax - xmin);
        const double y = fit.m * x + fit.b;
        const double ey = linear_fit_prediction_uncertainty(fit, x);
        g->SetPoint(n + i, x, y - ey);
    }

    g->SetFillColorAlpha(color, 0.18);
    g->SetLineColor(color);
    g->SetLineWidth(0);
    return g;
}


using PeriodLinearFitMap = std::map<std::string, LinearFitSummary>;

static LinearFitSummary copy_linear_fit_summary(const LinearFitSummary& src) {
    return src;
}

static void copy_fa18_inb_linear_fit_to_sp19_inb(PeriodLinearFitMap& fits,
                                                 const std::string& label) {
    auto it_fa18 = fits.find("Fa18 Inb");
    if (it_fa18 == fits.end() || !it_fa18->second.valid) {
        std::cout << "[current_dependence] WARNING: cannot copy Fa18 Inb e_theta fit to Sp19 Inb for "
                  << label << "; Fa18 Inb fit is missing or invalid." << std::endl;
        return;
    }

    fits["Sp19 Inb"] = copy_linear_fit_summary(it_fa18->second);
    std::cout << "[current_dependence] " << label
              << " Sp19 Inb e_theta data current-efficiency fit copied from Fa18 Inb."
              << std::endl;
}

static PeriodLinearFitMap build_e_theta_data_current_fits(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& data_trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    bool hide_sp19_inb_from_all_period_plots,
    int max_workers) {

    std::vector<KinematicVarConfig> vars;
    for (const KinematicVarConfig& v : kinematic_current_var_configs()) {
        if (v.key == "e_theta") {
            vars.push_back(v);
            break;
        }
    }

    if (vars.empty()) {
        fatal("[current_dependence] FATAL: internal error; missing e_theta kinematic configuration.");
    }

    PeriodKinematicAggMap aggs = build_kinematic_current_diagnostic_aggs(
        cfg,
        data_trees,
        charge_map,
        data_cuts,
        vars,
        use_second_column_charge_for_all_unpolarized,
        use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
        columns_3_to_5_charge_sum_scale,
        max_workers);

    PeriodLinearFitMap out;
    const KinematicVarConfig& v = vars.front();

    for (const std::string& period : PERIOD_ORDER) {
        if (hide_sp19_inb_from_replacement_plots(hide_sp19_inb_from_all_period_plots, period)) {
            continue;
        }

        auto it_period = aggs.find(period);
        if (it_period == aggs.end()) {
            continue;
        }

        auto it_var = it_period->second.find(v.key);
        if (it_var == it_period->second.end()) {
            continue;
        }

        std::vector<KinematicBinResult> bin_results;
        const std::vector<DataAgg>& bins = it_var->second;

        for (size_t ib = 0; ib < bins.size(); ++ib) {
            bin_results.push_back(current_factor_for_kinematic_bin(bins[ib],
                                                                   v.edges[ib],
                                                                   v.edges[ib + 1]));
        }

        LinearFitSummary fit = fit_kinematic_factor_linear(bin_results);
        out[period] = fit;

        if (fit.valid) {
            const double chi2ndf = (fit.ndf > 0) ? fit.chi2 / double(fit.ndf)
                                                : std::numeric_limits<double>::quiet_NaN();
            std::cout << "[current_dependence] " << cfg.csv_channel
                      << " DATA e_theta current-efficiency fit " << period
                      << ": f(theta_e)= " << fit.m << " * theta_e + " << fit.b
                      << " +/- slope=" << fit.m_err
                      << " intercept=" << fit.b_err
                      << " chi2/ndf=" << chi2ndf
                      << std::endl;
        } else {
            std::cout << "[current_dependence] WARNING: " << cfg.csv_channel
                      << " DATA e_theta current-efficiency fit is invalid for "
                      << period << "; integrated factor fallback will be used for this period."
                      << std::endl;
        }
    }

    return out;
}

static TLegend* make_kinematic_compact_legend(bool fit_legend) {
    TLegend* leg = new TLegend(0.14, 0.13, 0.92, fit_legend ? 0.27 : 0.24);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetNColumns(2);
    leg->SetTextSize(fit_legend ? 0.017 : 0.020);
    leg->SetMargin(0.10);
    leg->SetColumnSeparation(0.02);
    return leg;
}

static std::map<std::string, double> integrated_data_factor_by_period(const std::vector<PeriodResult>& results) {
    std::map<std::string, double> out;

    for (const PeriodResult& r : results) {
        out[r.period] = r.data_factor;
    }

    return out;
}

static std::map<std::string, LinearFitSummary> draw_kinematic_current_efficiency_diagnostics(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& data_trees,
    const std::map<std::string, TTree*>& gen_trees,
    const std::map<std::string, TTree*>& rec_trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    const TopoCutMap& mc_cuts,
    const std::vector<PeriodResult>& integrated_results,
    const std::string& output_dir,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    bool hide_sp19_inb_from_all_period_plots,
    int max_workers) {

    const std::vector<KinematicVarConfig> vars = kinematic_current_var_configs();
    const std::string odir = output_dir + "/" + cfg.output_token + "/kinematic_current_efficiency_diagnostics";
    mkdir_p(odir);

    std::cout << "[current_dependence] Building kinematic current-efficiency diagnostics in "
              << odir << std::endl;

    std::cout << "[current_dependence] Kinematic diagnostic branch map: "
              << "Q2 -> Q2, xB -> x, |t| -> fabs(t1), phi -> phi2 [deg], "
              << "e_theta -> e_theta [deg], p_theta -> p1_theta [deg], "
              << "g_theta -> p2_theta [deg]." << std::endl;

    PeriodKinematicAggMap data_aggs = build_kinematic_current_diagnostic_aggs(
        cfg,
        data_trees,
        charge_map,
        data_cuts,
        vars,
        use_second_column_charge_for_all_unpolarized,
        use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
        columns_3_to_5_charge_sum_scale,
        max_workers);

    PeriodKinematicMcAggMap mc_aggs = build_kinematic_current_diagnostic_mc_aggs(
        cfg,
        gen_trees,
        rec_trees,
        mc_cuts,
        vars,
        max_workers);

    std::map<std::string, double> integrated_data = integrated_data_factor_by_period(integrated_results);
    std::map<std::string, double> integrated_mc;
    for (const PeriodResult& r : integrated_results) {
        integrated_mc[r.period] = r.mc_factor;
    }

    using VarPeriodResults = std::map<std::string, std::map<std::string, std::vector<KinematicBinResult>>>;
    VarPeriodResults data_results_by_var_period;
    VarPeriodResults mc_results_by_var_period;
    std::map<std::string, std::map<std::string, LinearFitSummary>> data_fits_by_var_period;
    std::map<std::string, std::map<std::string, LinearFitSummary>> mc_fits_by_var_period;

    std::ofstream data_csv(odir + "/current_efficiency_vs_kinematics_data.csv");
    data_csv << "sample,period,variable,bin,x_low,x_high,x_center,current_efficiency_factor,current_efficiency_factor_stat,integrated_current_efficiency_factor,n_current_points,total_counts\n";

    std::ofstream mc_csv(odir + "/current_efficiency_vs_kinematics_mc.csv");
    mc_csv << "sample,period,variable,bin,x_low,x_high,x_center,current_efficiency_factor,current_efficiency_factor_stat,integrated_current_efficiency_factor,n_current_points,total_counts\n";

    std::ofstream data_fit_csv(odir + "/current_efficiency_vs_kinematics_data_linear_fits.csv");
    data_fit_csv << "sample,period,variable,slope,intercept,slope_stat,intercept_stat,cov_slope_intercept,chi2,ndf,npoints\n";

    std::ofstream mc_fit_csv(odir + "/current_efficiency_vs_kinematics_mc_linear_fits.csv");
    mc_fit_csv << "sample,period,variable,slope,intercept,slope_stat,intercept_stat,cov_slope_intercept,chi2,ndf,npoints\n";

    // Backward-compatible filenames from the first data-only implementation.
    std::ofstream legacy_data_csv(odir + "/current_efficiency_vs_kinematics.csv");
    legacy_data_csv << "period,variable,bin,x_low,x_high,x_center,x_error,current_efficiency_factor,current_efficiency_factor_stat,integrated_current_efficiency_factor,n_current_points,total_counts\n";

    std::ofstream legacy_data_fit_csv(odir + "/current_efficiency_vs_kinematics_linear_fits.csv");
    legacy_data_fit_csv << "period,variable,slope,intercept,slope_stat,intercept_stat,cov_slope_intercept,chi2,ndf,npoints\n";

    for (const KinematicVarConfig& v : vars) {
        for (const std::string& period : PERIOD_ORDER) {
            if (hide_sp19_inb_from_replacement_plots(hide_sp19_inb_from_all_period_plots, period)) {
                continue;
            }

            // DATA results.
            auto it_period_data = data_aggs.find(period);
            if (it_period_data != data_aggs.end()) {
                auto it_var = it_period_data->second.find(v.key);
                if (it_var != it_period_data->second.end()) {
                    std::vector<KinematicBinResult> bin_results;
                    const std::vector<DataAgg>& bins = it_var->second;

                    for (size_t ib = 0; ib < bins.size(); ++ib) {
                        KinematicBinResult br = current_factor_for_kinematic_bin(bins[ib], v.edges[ib], v.edges[ib + 1]);
                        bin_results.push_back(br);

                        const double int_factor = integrated_data.count(period) ? integrated_data[period] : std::numeric_limits<double>::quiet_NaN();
                        data_csv << "data," << period << "," << v.key << "," << ib << ","
                                 << br.x_low << "," << br.x_high << "," << br.x_center << ","
                                 << br.factor << "," << br.factor_err << "," << int_factor << ","
                                 << br.n_current_points << "," << br.total_counts << "\n";
                        legacy_data_csv << period << "," << v.key << "," << ib << ","
                                        << br.x_low << "," << br.x_high << "," << br.x_center << "," << br.x_err << ","
                                        << br.factor << "," << br.factor_err << "," << int_factor << ","
                                        << br.n_current_points << "," << br.total_counts << "\n";
                    }

                    data_results_by_var_period[v.key][period] = bin_results;
                    LinearFitSummary fit = fit_kinematic_factor_linear(bin_results);
                    data_fits_by_var_period[v.key][period] = fit;
                    data_fit_csv << "data," << period << "," << v.key << ","
                                 << fit.m << "," << fit.b << "," << fit.m_err << "," << fit.b_err << ","
                                 << fit.cov_mb << "," << fit.chi2 << "," << fit.ndf << "," << fit.npoints << "\n";
                    legacy_data_fit_csv << period << "," << v.key << ","
                                        << fit.m << "," << fit.b << "," << fit.m_err << "," << fit.b_err << ","
                                        << fit.cov_mb << "," << fit.chi2 << "," << fit.ndf << "," << fit.npoints << "\n";
                }
            }

            // MC results.
            auto it_period_mc = mc_aggs.find(period);
            if (it_period_mc != mc_aggs.end()) {
                auto it_var = it_period_mc->second.find(v.key);
                if (it_var != it_period_mc->second.end()) {
                    std::vector<KinematicBinResult> bin_results;
                    const std::vector<McKinematicBinAgg>& bins = it_var->second;

                    for (size_t ib = 0; ib < bins.size(); ++ib) {
                        KinematicBinResult br = current_factor_for_kinematic_mc_bin(bins[ib], v.edges[ib], v.edges[ib + 1]);
                        bin_results.push_back(br);

                        const double int_factor = integrated_mc.count(period) ? integrated_mc[period] : std::numeric_limits<double>::quiet_NaN();
                        mc_csv << "mc," << period << "," << v.key << "," << ib << ","
                               << br.x_low << "," << br.x_high << "," << br.x_center << ","
                               << br.factor << "," << br.factor_err << "," << int_factor << ","
                               << br.n_current_points << "," << br.total_counts << "\n";
                    }

                    mc_results_by_var_period[v.key][period] = bin_results;
                    LinearFitSummary fit = fit_kinematic_factor_linear(bin_results);
                    mc_fits_by_var_period[v.key][period] = fit;
                    mc_fit_csv << "mc," << period << "," << v.key << ","
                               << fit.m << "," << fit.b << "," << fit.m_err << "," << fit.b_err << ","
                               << fit.cov_mb << "," << fit.chi2 << "," << fit.ndf << "," << fit.npoints << "\n";
                }
            }
        }
    }

    auto draw_info_panel = [](const std::string& sample_label, bool fit_version) {
        gPad->SetLeftMargin(0.05);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.05);
        gPad->SetTopMargin(0.05);
        gPad->DrawFrame(0.0, 0.0, 1.0, 1.0);
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.043);
        lat.DrawLatex(0.12, 0.72, ("DVCS " + sample_label + " current-efficiency diagnostics").c_str());
        lat.SetTextSize(0.030);
        lat.DrawLatex(0.12, 0.60, "Points: per-bin current-efficiency factor");
        if (fit_version) {
            lat.DrawLatex(0.12, 0.52, "Lines: weighted linear fits, y = mx + b");
            lat.DrawLatex(0.12, 0.44, "Bands: one-standard-deviation fit uncertainty");
            lat.DrawLatex(0.12, 0.36, "Legend: #chi^{2}/ndf for each period");
        } else {
            lat.DrawLatex(0.12, 0.52, "Legend: integrated factor currently written to CSV");
        }
        lat.DrawLatex(0.12, 0.28, "Selection: global + topology-dependent exclusivity cuts");
    };

    auto draw_one_canvas = [&](const std::string& canvas_name,
                               const std::string& output_name,
                               const std::string& sample_label,
                               const VarPeriodResults& results_by_var_period,
                               const std::map<std::string, std::map<std::string, LinearFitSummary>>& fits_by_var_period,
                               const std::map<std::string, double>& integrated,
                               bool fit_version) {
        TCanvas c(canvas_name.c_str(), canvas_name.c_str(), 2200, 1100);
        c.Divide(4, 2, 0.001, 0.001);

        for (size_t iv = 0; iv < vars.size(); ++iv) {
            const KinematicVarConfig& v = vars[iv];
            c.cd((int)iv + 1);
            gPad->SetGrid(1, 1);
            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetRightMargin(0.04);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin(0.09);

            const double xmin = v.edges.front();
            const double xmax = v.edges.back();
            TH1F* frame = gPad->DrawFrame(xmin, 0.00, xmax, 1.40);
            frame->SetTitle((v.title + " dependence").c_str());
            frame->GetXaxis()->SetTitle(v.x_label.c_str());
            frame->GetYaxis()->SetTitle("Current-efficiency factor");
            frame->GetXaxis()->CenterTitle(true);
            frame->GetYaxis()->CenterTitle(true);
            frame->GetYaxis()->SetTitleOffset(1.35);

            TLegend* leg = make_kinematic_compact_legend(fit_version);

            // Draw shaded fit bands first, then fit lines, then points.
            if (fit_version) {
                for (const std::string& period : PERIOD_ORDER) {
                    if (hide_sp19_inb_from_replacement_plots(hide_sp19_inb_from_all_period_plots, period)) continue;
                    const auto it_var = fits_by_var_period.find(v.key);
                    if (it_var == fits_by_var_period.end()) continue;
                    const auto it_fit = it_var->second.find(period);
                    if (it_fit == it_var->second.end() || !it_fit->second.valid) continue;
                    TGraph* band = make_linear_fit_band_graph(it_fit->second, xmin, xmax, period_color(period));
                    band->Draw("F SAME");
                }

                for (const std::string& period : PERIOD_ORDER) {
                    if (hide_sp19_inb_from_replacement_plots(hide_sp19_inb_from_all_period_plots, period)) continue;
                    const auto it_var = fits_by_var_period.find(v.key);
                    if (it_var == fits_by_var_period.end()) continue;
                    const auto it_fit = it_var->second.find(period);
                    if (it_fit == it_var->second.end() || !it_fit->second.valid) continue;
                    TGraph* line = make_linear_fit_line_graph(it_fit->second, xmin, xmax, period_color(period));
                    line->Draw("L SAME");
                }
            }

            for (const std::string& period : PERIOD_ORDER) {
                if (hide_sp19_inb_from_replacement_plots(hide_sp19_inb_from_all_period_plots, period)) {
                    continue;
                }

                const auto it_var = results_by_var_period.find(v.key);
                if (it_var == results_by_var_period.end()) continue;
                const auto it_period = it_var->second.find(period);
                if (it_period == it_var->second.end()) continue;

                TGraphErrors* g = make_kinematic_factor_graph(it_period->second, period_color(period));
                g->Draw("P SAME");

                std::ostringstream lab;
                lab << period;

                if (fit_version) {
                    const auto it_fit_var = fits_by_var_period.find(v.key);
                    if (it_fit_var != fits_by_var_period.end()) {
                        const auto it_fit = it_fit_var->second.find(period);
                        if (it_fit != it_fit_var->second.end() && it_fit->second.valid && it_fit->second.ndf > 0) {
                            lab << " #chi^{2}/ndf=" << std::fixed << std::setprecision(1)
                                << (it_fit->second.chi2 / double(it_fit->second.ndf));
                        } else {
                            lab << " #chi^{2}/ndf=n/a";
                        }
                    }
                } else {
                    auto it_int = integrated.find(period);
                    if (it_int != integrated.end() && std::isfinite(it_int->second)) {
                        lab << " int=" << std::fixed << std::setprecision(3) << it_int->second;
                    }
                }

                leg->AddEntry(g, lab.str().c_str(), "pe");
            }

            leg->Draw();
        }

        c.cd(8);
        draw_info_panel(sample_label, fit_version);
        c.SaveAs((odir + "/" + output_name).c_str());
    };

    draw_one_canvas("c_kinematic_current_efficiency_data",
                    "current_efficiency_vs_kinematics.png",
                    "data",
                    data_results_by_var_period,
                    data_fits_by_var_period,
                    integrated_data,
                    false);

    draw_one_canvas("c_kinematic_current_efficiency_data_linear_fits",
                    "current_efficiency_vs_kinematics_linear_fits.png",
                    "data",
                    data_results_by_var_period,
                    data_fits_by_var_period,
                    integrated_data,
                    true);

    draw_one_canvas("c_kinematic_current_efficiency_mc",
                    "current_efficiency_vs_kinematics_mc.png",
                    "MC",
                    mc_results_by_var_period,
                    mc_fits_by_var_period,
                    integrated_mc,
                    false);

    draw_one_canvas("c_kinematic_current_efficiency_mc_linear_fits",
                    "current_efficiency_vs_kinematics_mc_linear_fits.png",
                    "MC",
                    mc_results_by_var_period,
                    mc_fits_by_var_period,
                    integrated_mc,
                    true);

    std::map<std::string, LinearFitSummary> e_theta_fits;
    auto it_e_theta = data_fits_by_var_period.find("e_theta");
    if (it_e_theta != data_fits_by_var_period.end()) {
        e_theta_fits = it_e_theta->second;
    }
    return e_theta_fits;
}

static std::vector<PeriodResult> run_channel_study(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& data_trees,
    const std::map<std::string, TTree*>& gen_trees,
    const std::map<std::string, TTree*>& rec_trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    const TopoCutMap& mc_cuts,
    const std::string& output_dir,
    int max_workers,
    bool process_mc,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    bool hide_sp19_inb_from_all_period_plots) {

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
                                        data_cuts,
                                        use_second_column_charge_for_all_unpolarized,
                                        use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                                        columns_3_to_5_charge_sum_scale);

        std::lock_guard<std::mutex> lock(data_mutex);
        data_aggs[agg.period] = std::move(agg);
    }

    std::map<std::string, McAgg> mc_by_period_current;
    std::vector<McAgg> mc_aggs;

    if (!process_mc) {
        std::cout << "[current_dependence] " << cfg.csv_channel
                  << ": skipping MC current-dependence scan; MC factor will be built downstream."
                  << std::endl;
    }

    if (process_mc) {
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

        for (const auto& kv : mc_by_period_current) {
            mc_aggs.push_back(kv.second);
        }
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

        if (process_mc) {
            R.mc_points = mc_points_from_aggs(mc_aggs, period);
            R.mc_fit = fit_points(R.mc_points);

            const int ref = reference_current_nA(period);
            R.mc_factor = rel_at_current((double)ref, R.mc_fit);
            R.mc_factor_err = rel_err_at_current((double)ref, R.mc_fit);
        }

        results.push_back(R);

        const std::string odir = output_dir + "/" + cfg.output_token + "/" + period_dir(period);
        mkdir_p(odir);

        write_points_csv(odir + "/data_current_points.csv", R.data_points);

        if (process_mc) {
            write_points_csv(odir + "/mc_current_points.csv", R.mc_points);
        }

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
                  << " data_slope=" << fit_percent_slope(R.data_fit)
                  << " +/- " << fit_percent_slope_err(R.data_fit)
                  << " %/nA";

        if (process_mc) {
            std::cout << " mc_factor=" << R.mc_factor
                      << " +/- " << R.mc_factor_err;
        } else {
            std::cout << " mc_factor=skipped";
        }

        std::cout << std::endl;
    }

    draw_all_period_current_canvas(output_dir + "/" + cfg.output_token + "/all_periods_current_dependence.png",
                                   cfg.title,
                                   results,
                                   hide_sp19_inb_from_all_period_plots);

    draw_all_period_data_mc_canvas(output_dir + "/" + cfg.output_token + "/all_periods_current_dependence_data_mc.png",
                                   cfg.title,
                                   results,
                                   hide_sp19_inb_from_all_period_plots);

    write_summary_csv(output_dir + "/" + cfg.output_token + "/period_summary_raw_fit_values.csv", results);

    return results;
}

static void apply_eppi0_mc_factor_from_dvcs_ratio(std::vector<PeriodResult>& eppi0_results,
                                                  const std::vector<PeriodResult>& dvcs_results) {
    std::map<std::string, PeriodResult> dvcs_by_period;

    for (const PeriodResult& r : dvcs_results) {
        dvcs_by_period[r.period] = r;
    }

    auto require_finite_positive = [](double value,
                                      const std::string& quantity,
                                      const std::string& period) {
        if (!std::isfinite(value) || !(value > 0.0)) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: cannot build eppi0 MC current factor for "
               << period << "; " << quantity << " is not finite and positive: "
               << value;
            fatal(ss.str());
        }
    };

    auto require_finite_nonnegative = [](double value,
                                         const std::string& quantity,
                                         const std::string& period) {
        if (!std::isfinite(value) || value < 0.0) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: cannot build eppi0 MC current factor for "
               << period << "; " << quantity << " is not finite and non-negative: "
               << value;
            fatal(ss.str());
        }
    };

    for (PeriodResult& eppi0 : eppi0_results) {
        auto it = dvcs_by_period.find(eppi0.period);

        if (it == dvcs_by_period.end()) {
            std::ostringstream ss;
            ss << "[current_dependence] FATAL: missing DVCS current-dependence result for period "
               << eppi0.period << " while constructing eppi0 MC factor.";
            fatal(ss.str());
        }

        const PeriodResult& dvcs = it->second;

        require_finite_positive(eppi0.data_factor, "eppi0 data factor", eppi0.period);
        require_finite_positive(dvcs.mc_factor, "DVCS MC factor", eppi0.period);
        require_finite_positive(dvcs.data_factor, "DVCS data factor", eppi0.period);

        require_finite_nonnegative(eppi0.data_factor_err, "eppi0 data factor uncertainty", eppi0.period);
        require_finite_nonnegative(dvcs.mc_factor_err, "DVCS MC factor uncertainty", eppi0.period);
        require_finite_nonnegative(dvcs.data_factor_err, "DVCS data factor uncertainty", eppi0.period);

        const double scale = dvcs.mc_factor / dvcs.data_factor;
        const double corrected = eppi0.data_factor * scale;

        const double rel_var =
            std::pow(eppi0.data_factor_err / eppi0.data_factor, 2.0) +
            std::pow(dvcs.mc_factor_err / dvcs.mc_factor, 2.0) +
            std::pow(dvcs.data_factor_err / dvcs.data_factor, 2.0);

        const double corrected_err = std::fabs(corrected) * std::sqrt(rel_var);

        eppi0.mc_factor = corrected;
        eppi0.mc_factor_err = corrected_err;

        std::cout << "[current_dependence] Built eppi0 MC current factor for "
                  << eppi0.period
                  << " using eppi0_data * (dvcs_mc / dvcs_data): "
                  << "eppi0_data=" << eppi0.data_factor
                  << " +/- " << eppi0.data_factor_err
                  << " dvcs_mc=" << dvcs.mc_factor
                  << " +/- " << dvcs.mc_factor_err
                  << " dvcs_data=" << dvcs.data_factor
                  << " +/- " << dvcs.data_factor_err
                  << " corrected_eppi0_mc=" << eppi0.mc_factor
                  << " +/- " << eppi0.mc_factor_err
                  << std::endl;
    }
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

static std::vector<PeriodResult> with_mc_factors_forced_to_unity(
    const std::vector<PeriodResult>& in) {

    std::vector<PeriodResult> out = in;

    for (PeriodResult& r : out) {
        r.mc_factor = 1.0;
        r.mc_factor_err = 0.0;
    }

    return out;
}

// -----------------------------------------------------------------------------
// Deprecated eppi0 AAOGEN normalization compatibility fallback
// -----------------------------------------------------------------------------

static const std::vector<std::string>& eppi0_normalization_regions() {
    static const std::vector<std::string> v = {
        "Sector 1",
        "Sector 2",
        "Sector 3",
        "Sector 4",
        "Sector 5",
        "Sector 6",
        "CD"
    };

    return v;
}

static std::string eppi0_norm_factor_col(const std::string& period,
                                         const std::string& region) {
    return "eppi0 cross-section normalization factor, ep->eppi0, data_over_mc, " +
           region + ", " + period;
}

static std::string eppi0_norm_cubic_col(const std::string& period,
                                        const std::string& region) {
    return "eppi0 cross-section normalization cubic, ep->eppi0, data_over_mc, " +
           region + ", " + period;
}

static void write_unity_eppi0_normalization_columns(CSV& csv) {
    const std::string factor_cell = tuple2(1.0, 0.0);
    const std::string cubic_cell  = tuple4(1.0, 0.0, 0.0, 0.0);

    long long n_factor_cells = 0;
    long long n_cubic_cells = 0;

    for (const std::string& period : CSV_PERIOD_ORDER) {
        for (const std::string& region : eppi0_normalization_regions()) {
            const std::string factor_col = eppi0_norm_factor_col(period, region);
            const std::string cubic_col  = eppi0_norm_cubic_col(period, region);

            const int c_factor = col_strict(csv, factor_col);
            const int c_cubic  = col_strict(csv, cubic_col);

            for (auto& row : csv.rows) {
                row[c_factor] = factor_cell;
                row[c_cubic]  = cubic_cell;
                ++n_factor_cells;
                ++n_cubic_cells;
            }
        }
    }

    std::cout << "[current_dependence] Wrote unity fallback eppi0 normalization columns: "
              << "factor_cells=" << n_factor_cells
              << " cubic_cells=" << n_cubic_cells
              << " using factor=" << factor_cell
              << " cubic=" << cubic_cell
              << std::endl;
}

// -----------------------------------------------------------------------------
// Apply saved current-efficiency factors to DATA and MC yield columns
// -----------------------------------------------------------------------------

static const std::vector<std::string>& topology_labels() {
    static const std::vector<std::string> v = {
        "(FD, FD)",
        "(CD, FD)",
        "(CD, FT)"
    };

    return v;
}

static const std::vector<std::string>& helicity_labels() {
    static const std::vector<std::string> v = {
        "unpol",
        "pos",
        "neg"
    };

    return v;
}

static bool parse_tuple_numbers_cell(const std::string& cell, std::vector<double>& vals) {
    vals.clear();

    std::string s;
    s.reserve(cell.size());

    for (char c : cell) {
        if (c == '(' || c == ')' || c == '"') {
            s.push_back(' ');
        } else {
            s.push_back(c);
        }
    }

    std::stringstream ss(s);
    std::string part;

    while (std::getline(ss, part, ',')) {
        char* endp = nullptr;
        const double v = std::strtod(part.c_str(), &endp);

        if (endp == part.c_str()) {
            return false;
        }

        vals.push_back(v);
    }

    return !vals.empty();
}

static bool parse_value_stat_sys(const std::string& cell,
                                 double& value,
                                 double& stat,
                                 double& sys) {
    value = std::numeric_limits<double>::quiet_NaN();
    stat = std::numeric_limits<double>::quiet_NaN();
    sys = 0.0;

    std::vector<double> vals;

    if (!parse_tuple_numbers_cell(cell, vals) || vals.size() < 3) {
        return false;
    }

    value = vals[0];
    stat = vals[1];
    sys = vals[2];

    return std::isfinite(value) &&
           std::isfinite(stat) &&
           std::isfinite(sys) &&
           stat >= 0.0 &&
           sys >= 0.0;
}

static std::string format_triple(double v, double stat, double sys) {
    if (!std::isfinite(v) || !std::isfinite(stat) || !std::isfinite(sys)) {
        return "";
    }

    std::ostringstream ss;
    ss << std::setprecision(12)
       << "(" << v << "," << stat << "," << sys << ")";
    return ss.str();
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

static std::string raw_yield_topo_col(const std::string& channel,
                                      const std::string& topo,
                                      const std::string& period,
                                      const std::string& helicity) {
    return "raw yield, " + channel + ", " + topo + ", exp, " + period + ", " + helicity;
}

static std::string normalized_raw_yield_topo_col(const std::string& channel,
                                                 const std::string& topo,
                                                 const std::string& period,
                                                 const std::string& helicity) {
    return "normalized raw yield, " + channel + ", " + topo + ", exp, " + period + ", " + helicity;
}

static void read_current_factor_from_csv(const CSV& csv,
                                         const ChannelConfig& cfg,
                                         const std::string& sample,
                                         const std::string& period,
                                         double& factor,
                                         double& factor_stat) {
    const std::string colname = current_eff_col(cfg, sample, period);
    const int c = col_strict(csv, colname);

    if (csv.rows.empty()) {
        fatal("[current_dependence] FATAL: cannot read current factor from empty CSV.");
    }

    std::vector<double> vals;

    if (!parse_tuple_numbers_cell(csv.rows.front()[c], vals) ||
        vals.size() < 2 ||
        !(std::isfinite(vals[0]) && vals[0] > 0.0) ||
        !(std::isfinite(vals[1]) && vals[1] >= 0.0)) {
        std::ostringstream ss;
        ss << "[current_dependence] FATAL: invalid current-efficiency factor tuple in column '"
           << colname << "': '" << csv.rows.front()[c] << "'.";
        fatal(ss.str());
    }

    factor = vals[0];
    factor_stat = vals[1];
}

static void read_mc_current_factor_from_csv(const CSV& csv,
                                            const ChannelConfig& cfg,
                                            const std::string& period,
                                            double& factor,
                                            double& factor_stat) {
    read_current_factor_from_csv(csv, cfg, "mc", period, factor, factor_stat);
}

static void read_exp_current_factor_from_csv(const CSV& csv,
                                             const ChannelConfig& cfg,
                                             const std::string& period,
                                             double& factor,
                                             double& factor_stat) {
    read_current_factor_from_csv(csv, cfg, "exp", period, factor, factor_stat);
}



static bool parse_plain_number_cell(const std::string& cell, double& value) {
    value = std::numeric_limits<double>::quiet_NaN();

    const char* begin = cell.c_str();
    while (*begin && std::isspace(static_cast<unsigned char>(*begin))) {
        ++begin;
    }

    if (!*begin) {
        return false;
    }

    char* endp = nullptr;
    value = std::strtod(begin, &endp);
    if (endp == begin) {
        return false;
    }

    while (endp && *endp) {
        if (!std::isspace(static_cast<unsigned char>(*endp))) {
            return false;
        }
        ++endp;
    }

    return std::isfinite(value);
}

static std::string e_theta_avg_col(const std::string& period) {
    return "e_theta, " + period;
}

static bool finite_positive_factor(double f, double ferr) {
    return std::isfinite(f) && f > 0.0 && std::isfinite(ferr) && ferr >= 0.0;
}

static bool evaluate_e_theta_current_fit(const LinearFitSummary& fit,
                                         double e_theta_deg,
                                         double& factor,
                                         double& factor_stat) {
    factor = std::numeric_limits<double>::quiet_NaN();
    factor_stat = std::numeric_limits<double>::quiet_NaN();

    if (!fit.valid || !std::isfinite(e_theta_deg)) {
        return false;
    }

    factor = fit.m * e_theta_deg + fit.b;
    factor_stat = linear_fit_prediction_uncertainty(fit, e_theta_deg);

    return finite_positive_factor(factor, factor_stat);
}

static long long apply_one_current_correction_column(CSV& csv,
                                                     const std::string& src_col,
                                                     const std::string& dst_col,
                                                     double factor,
                                                     double factor_stat) {
    const int c_src = col_strict(csv, src_col);
    const int c_dst = col_strict(csv, dst_col);

    if (!(std::isfinite(factor) && factor > 0.0) ||
        !(std::isfinite(factor_stat) && factor_stat >= 0.0)) {
        fatal("[current_dependence] FATAL: invalid current factor passed to correction writer.");
    }

    long long n_positive = 0;

    for (auto& row : csv.rows) {
        double raw = 0.0;
        double raw_stat = 0.0;
        double raw_sys = 0.0;

        if (!parse_value_stat_sys(row[c_src], raw, raw_stat, raw_sys) ||
            !(std::isfinite(raw) && raw > 0.0)) {
            row[c_dst].clear();
            continue;
        }

        const double corrected = raw / factor;

        // The current-efficiency factor is a fitted statistical correction applied
        // once to the already-binned yield. Its uncertainty therefore enters as
        // the ordinary multiplicative error-propagation term:
        //
        //   Var(raw/f) = Var(raw)/f^2 + raw^2 Var(f)/f^4 .
        //
        // Do not sum this term event-by-event; that artificially suppresses the
        // factor uncertainty by roughly 1/sqrt(N).
        const double var_stat =
            (raw_stat * raw_stat) / (factor * factor) +
            (raw * raw * factor_stat * factor_stat) /
            (factor * factor * factor * factor);

        const double corrected_stat = (var_stat > 0.0) ? std::sqrt(var_stat) : 0.0;
        const double corrected_sys = raw_sys / factor;

        row[c_dst] = format_triple(corrected, corrected_stat, corrected_sys);

        if (corrected > 0.0) {
            ++n_positive;
        }
    }

    return n_positive;
}

static long long apply_one_data_current_correction_column(CSV& csv,
                                                          const std::string& src_col,
                                                          const std::string& dst_col,
                                                          const std::string& period,
                                                          double integrated_factor,
                                                          double integrated_factor_stat,
                                                          const LinearFitSummary* e_theta_fit,
                                                          bool use_e_theta_linear_fit,
                                                          long long& n_fit_rows,
                                                          long long& n_integrated_fallback_rows) {
    const int c_src = col_strict(csv, src_col);
    const int c_dst = col_strict(csv, dst_col);
    const int c_eth = col_optional(csv, e_theta_avg_col(period));

    if (!finite_positive_factor(integrated_factor, integrated_factor_stat)) {
        fatal("[current_dependence] FATAL: invalid integrated DATA current factor passed to correction writer.");
    }

    long long n_positive = 0;

    for (auto& row : csv.rows) {
        double raw = 0.0;
        double raw_stat = 0.0;
        double raw_sys = 0.0;

        if (!parse_value_stat_sys(row[c_src], raw, raw_stat, raw_sys) ||
            !(std::isfinite(raw) && raw > 0.0)) {
            row[c_dst].clear();
            continue;
        }

        double factor = integrated_factor;
        double factor_stat = integrated_factor_stat;
        bool used_fit = false;

        if (use_e_theta_linear_fit && e_theta_fit && e_theta_fit->valid && c_eth >= 0) {
            double eth = std::numeric_limits<double>::quiet_NaN();
            if (parse_plain_number_cell(row[c_eth], eth)) {
                double f_fit = std::numeric_limits<double>::quiet_NaN();
                double f_fit_stat = std::numeric_limits<double>::quiet_NaN();
                if (evaluate_e_theta_current_fit(*e_theta_fit, eth, f_fit, f_fit_stat)) {
                    factor = f_fit;
                    factor_stat = f_fit_stat;
                    used_fit = true;
                }
            }
        }

        if (used_fit) {
            ++n_fit_rows;
        } else {
            ++n_integrated_fallback_rows;
        }

        const double corrected = raw / factor;
        const double var_stat =
            (raw_stat * raw_stat) / (factor * factor) +
            (raw * raw * factor_stat * factor_stat) /
            (factor * factor * factor * factor);

        const double corrected_stat = (var_stat > 0.0) ? std::sqrt(var_stat) : 0.0;
        const double corrected_sys = raw_sys / factor;

        row[c_dst] = format_triple(corrected, corrected_stat, corrected_sys);

        if (corrected > 0.0) {
            ++n_positive;
        }
    }

    return n_positive;
}

static long long apply_data_current_corrections_for_channel(CSV& csv,
                                                            const ChannelConfig& cfg,
                                                            const std::string& log_label,
                                                            const PeriodLinearFitMap& e_theta_fits,
                                                            bool use_e_theta_linear_data_current_efficiency) {
    long long total_positive = 0;
    long long total_fit_rows = 0;
    long long total_fallback_rows = 0;

    for (const std::string& period : CSV_PERIOD_ORDER) {
        double f_exp = 0.0;
        double f_exp_stat = 0.0;

        read_exp_current_factor_from_csv(csv, cfg, period, f_exp, f_exp_stat);

        const LinearFitSummary* fit_ptr = nullptr;
        auto it_fit = e_theta_fits.find(period);
        if (it_fit != e_theta_fits.end() && it_fit->second.valid) {
            fit_ptr = &it_fit->second;
        }

        if (use_e_theta_linear_data_current_efficiency) {
            if (fit_ptr) {
                std::cout << "[current_dependence] DATA current correction for "
                          << log_label << " period=" << period
                          << " will use e_theta-dependent factor f(theta_e)="
                          << fit_ptr->m << "*theta_e + " << fit_ptr->b
                          << " with integrated fallback=" << f_exp
                          << " +/- " << f_exp_stat
                          << std::endl;
            } else {
                std::cout << "[current_dependence] WARNING: DATA current correction for "
                          << log_label << " period=" << period
                          << " has no valid e_theta fit; using integrated factor="
                          << f_exp << " +/- " << f_exp_stat
                          << std::endl;
            }
        }

        for (const std::string& topo : topology_labels()) {
            for (const std::string& hel : helicity_labels()) {
                const std::string src_col =
                    raw_yield_topo_col(cfg.csv_channel, topo, period, hel);

                const std::string dst_col =
                    normalized_raw_yield_topo_col(cfg.csv_channel, topo, period, hel);

                const int c_src = col_optional(csv, src_col);
                const int c_dst = col_optional(csv, dst_col);

                // Not all periods have helicity-resolved data columns.
                // Spring 2018 currently has only unpolarized data columns.
                if (c_src < 0 && c_dst < 0) {
                    continue;
                }

                if (c_src < 0 || c_dst < 0) {
                    std::ostringstream ss;
                    ss << "[current_dependence] FATAL: mismatched DATA current-correction columns: "
                       << "src='" << src_col << "' exists=" << (c_src >= 0)
                       << " dst='" << dst_col << "' exists=" << (c_dst >= 0);
                    fatal(ss.str());
                }

                long long n_fit_rows = 0;
                long long n_fallback_rows = 0;

                const long long n = apply_one_data_current_correction_column(
                    csv,
                    src_col,
                    dst_col,
                    period,
                    f_exp,
                    f_exp_stat,
                    fit_ptr,
                    use_e_theta_linear_data_current_efficiency,
                    n_fit_rows,
                    n_fallback_rows);

                total_positive += n;
                total_fit_rows += n_fit_rows;
                total_fallback_rows += n_fallback_rows;

                std::cout << "[current_dependence] Applied DATA current correction for "
                          << log_label
                          << " period=" << period
                          << " topo=" << topo
                          << " hel=" << hel;

                if (use_e_theta_linear_data_current_efficiency && fit_ptr) {
                    std::cout << " e_theta_fit=(m=" << fit_ptr->m
                              << ", b=" << fit_ptr->b << ")"
                              << " integrated_fallback=" << f_exp
                              << " +/- " << f_exp_stat
                              << " fit_rows=" << n_fit_rows
                              << " fallback_rows=" << n_fallback_rows;
                } else {
                    std::cout << " factor=" << f_exp
                              << " +/- " << f_exp_stat;
                }

                std::cout << " positive_rows=" << n
                          << std::endl;
            }
        }
    }

    if (use_e_theta_linear_data_current_efficiency) {
        std::cout << "[current_dependence] DATA current correction for "
                  << log_label << " used e_theta fit rows=" << total_fit_rows
                  << " integrated fallback rows=" << total_fallback_rows
                  << std::endl;
    }

    return total_positive;
}

static void apply_all_data_current_corrections(CSV& csv,
                                               const ChannelConfig& dvcs,
                                               const ChannelConfig& eppi0,
                                               const PeriodLinearFitMap& dvcs_e_theta_fits,
                                               const PeriodLinearFitMap& eppi0_e_theta_fits,
                                               bool use_e_theta_linear_data_current_efficiency) {
    const long long n_dvcs =
        apply_data_current_corrections_for_channel(csv,
                                                   dvcs,
                                                   "ep->epg",
                                                   dvcs_e_theta_fits,
                                                   use_e_theta_linear_data_current_efficiency);

    const long long n_eppi0 =
        apply_data_current_corrections_for_channel(csv,
                                                   eppi0,
                                                   "ep->eppi0",
                                                   eppi0_e_theta_fits,
                                                   use_e_theta_linear_data_current_efficiency);

    if (n_dvcs <= 0) {
        fatal("[current_dependence] FATAL: applying DVCS DATA current corrections produced zero positive cells.");
    }

    if (n_eppi0 <= 0) {
        fatal("[current_dependence] FATAL: applying eppi0 DATA current corrections produced zero positive cells.");
    }

    std::cout << "[current_dependence] DATA current-corrected normalized raw yield cells written: "
              << "ep->epg=" << n_dvcs
              << " ep->eppi0=" << n_eppi0;

    if (use_e_theta_linear_data_current_efficiency) {
        std::cout << " using e_theta-dependent DATA current-efficiency fits";
    }

    std::cout << std::endl;
}

static long long apply_mc_current_corrections_for_channel(CSV& csv,
                                                          const std::string& source_channel,
                                                          const ChannelConfig& factor_cfg,
                                                          const std::string& log_label) {
    long long total_positive = 0;

    for (const std::string& period : CSV_PERIOD_ORDER) {
        double f_mc = 0.0;
        double f_mc_stat = 0.0;

        read_mc_current_factor_from_csv(csv, factor_cfg, period, f_mc, f_mc_stat);

        const long long n_total = apply_one_current_correction_column(
            csv,
            rec_yield_total_col(source_channel, period),
            rec_current_corrected_total_col(source_channel, period),
            f_mc,
            f_mc_stat);

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
                f_mc,
                f_mc_stat);

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
                                             const ChannelConfig& eppi0,
                                             bool use_epg_mc_current_factor_for_eppi0_bkg) {
    const long long n_dvcs = apply_mc_current_corrections_for_channel(csv,
                                                                      "ep->epg",
                                                                      dvcs,
                                                                      "ep->epg");

    const long long n_eppi0 = apply_mc_current_corrections_for_channel(csv,
                                                                       "ep->eppi0",
                                                                       eppi0,
                                                                       "ep->eppi0");

    const ChannelConfig& eppi0_bkg_factor_cfg =
        use_epg_mc_current_factor_for_eppi0_bkg ? dvcs : eppi0;

    const long long n_eppi0_bkg = apply_mc_current_corrections_for_channel(csv,
                                                                           "ep->eppi0->epg",
                                                                           eppi0_bkg_factor_cfg,
                                                                           use_epg_mc_current_factor_for_eppi0_bkg
                                                                               ? "ep->eppi0->epg (using ep->epg MC factor)"
                                                                               : "ep->eppi0->epg (using ep->eppi0 MC factor)");

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

        PeriodLinearFitMap dvcs_e_theta_data_fits;
        PeriodLinearFitMap eppi0_e_theta_data_fits;

        if (options.override_to_unity) {
            std::cout << "[current_dependence] Override enabled: writing all current-efficiency factors as (1,0)."
                      << std::endl;

            write_override_unity(csv, dvcs);
            write_override_unity(csv, eppi0);

            // Compatibility fallback now that eppi0_normalization.cpp may be
            // disabled. If the real eppi0 normalization stage is re-enabled
            // after current_dependence.cpp, it can overwrite these unity values.
            write_unity_eppi0_normalization_columns(csv);

            // Even in unity-override mode, this stage now owns the production
            // of normalized raw DATA yields. With f=(1,0), this copies raw
            // yields into normalized raw yields while preserving uncertainties.
            apply_all_data_current_corrections(csv,
                                           dvcs,
                                           eppi0,
                                           dvcs_e_theta_data_fits,
                                           eppi0_e_theta_data_fits,
                                           options.use_e_theta_linear_data_current_efficiency);

            // With f=(1,0), this copies reconstructed MC yields into
            // reconstructed current-corrected MC yield columns.
            apply_all_mc_current_corrections(csv,
                                     dvcs,
                                     eppi0,
                                     options.use_epg_mc_current_factor_for_eppi0_bkg);

            write_csv_atomic(csv_path, csv);

            std::cout << "[current_dependence] Updated CSV with unity current-efficiency factors, "
                      << "unity eppi0 normalization fallback columns, and unity current-corrected DATA/MC yields: "
                      << csv_path << std::endl;

            return true;
        }

        mkdir_p(options.output_dir);

        const std::unordered_map<int, ChargeEntry> charge_map =
            read_charge_csv(options.charge_csv_path);

        std::cout << "[current_dependence] Unpolarized data charge-normalization mode: ";

        if (options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized) {
            std::cout << "Fa18/Sp19 scaled columns-3-to-5 mode: charge = "
                      << options.columns_3_to_5_charge_sum_scale
                      << " * (column 3 + column 4 + column 5) for Fa18 and Sp19; "
                      << "Spring 2018 uses column 2."
                      << std::endl;
        } else if (options.use_second_column_charge_for_all_unpolarized) {
            std::cout << "column 2 for all periods. Spring 2018 is also column 2."
                      << std::endl;
        } else {
            std::cout << "legacy mixed mode: Spring 2018 uses column 2; "
                      << "Fall 2018 and Spring 2019 use column 3 + column 4."
                      << std::endl;
        }

        if (options.use_fa18_inb_current_efficiency_for_sp19_inb) {
            std::cout << "[current_dependence] Sp19 Inb replacement mode enabled: "
                      << "Sp19 Inb current-efficiency factors written to the CSV "
                      << "will be copied from Fa18 Inb. The raw Sp19 Inb scan "
                      << "will still be processed and plotted for diagnostics."
                      << std::endl;
        } else {
            std::cout << "[current_dependence] Sp19 Inb replacement mode disabled: "
                      << "the directly fitted Sp19 Inb luminosity-scan factor "
                      << "will be written to the CSV."
                      << std::endl;
        }

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
                              options.max_workers,
                              true,
                              options.use_second_column_charge_for_all_unpolarized,
                              options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                              options.columns_3_to_5_charge_sum_scale,
                              options.use_fa18_inb_current_efficiency_for_sp19_inb);

        std::vector<PeriodResult> eppi0_results =
            run_channel_study(eppi0,
                              eppi0DataTrees,
                              eppi0GenMcTrees,
                              eppi0RecMcTrees,
                              charge_map,
                              data_cuts,
                              mc_cuts,
                              options.output_dir,
                              options.max_workers,
                              false,
                              options.use_second_column_charge_for_all_unpolarized,
                              options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                              options.columns_3_to_5_charge_sum_scale,
                              options.use_fa18_inb_current_efficiency_for_sp19_inb);

        if (options.use_fa18_inb_current_efficiency_for_sp19_inb) {
            replace_sp19_inb_factors_with_fa18_inb(dvcs_results, dvcs.csv_channel);
        }

        dvcs_e_theta_data_fits =
            draw_kinematic_current_efficiency_diagnostics(
                dvcs,
                dvcsDataTrees,
                dvcsGenMcTrees,
                dvcsRecMcTrees,
                charge_map,
                data_cuts,
                mc_cuts,
                dvcs_results,
                options.output_dir,
                options.use_second_column_charge_for_all_unpolarized,
                options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                options.columns_3_to_5_charge_sum_scale,
                options.use_fa18_inb_current_efficiency_for_sp19_inb,
                options.max_workers);

        if (options.use_e_theta_linear_data_current_efficiency) {
            std::cout << "[current_dependence] Building ep->eppi0 DATA e_theta current-efficiency fits "
                      << "for normalized raw-yield corrections." << std::endl;
            eppi0_e_theta_data_fits = build_e_theta_data_current_fits(
                eppi0,
                eppi0DataTrees,
                charge_map,
                data_cuts,
                options.use_second_column_charge_for_all_unpolarized,
                options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                options.columns_3_to_5_charge_sum_scale,
                options.use_fa18_inb_current_efficiency_for_sp19_inb,
                options.max_workers);

            if (options.use_fa18_inb_current_efficiency_for_sp19_inb) {
                copy_fa18_inb_linear_fit_to_sp19_inb(dvcs_e_theta_data_fits, dvcs.csv_channel);
                copy_fa18_inb_linear_fit_to_sp19_inb(eppi0_e_theta_data_fits, eppi0.csv_channel);
            }
        }

        // eppi0 reconstructed MC current dependence is still derived from the
        // DVCS MC/data ratio because the eppi0 MC scan is intentionally skipped.
        apply_eppi0_mc_factor_from_dvcs_ratio(eppi0_results, dvcs_results);

        if (options.use_fa18_inb_current_efficiency_for_sp19_inb) {
            replace_sp19_inb_factors_with_fa18_inb(eppi0_results, eppi0.csv_channel);
        }

        draw_all_period_current_canvas(options.output_dir + "/" + dvcs.output_token + "/all_periods_current_dependence_csv_values.png",
                                       dvcs.title + " CSV values",
                                       dvcs_results,
                                       options.use_fa18_inb_current_efficiency_for_sp19_inb);

        draw_all_period_data_mc_canvas(options.output_dir + "/" + dvcs.output_token + "/all_periods_current_dependence_data_mc_csv_values.png",
                                       dvcs.title + " CSV values",
                                       dvcs_results,
                                       options.use_fa18_inb_current_efficiency_for_sp19_inb);

        draw_all_period_current_canvas(options.output_dir + "/" + eppi0.output_token + "/all_periods_current_dependence_csv_values.png",
                                       eppi0.title + " CSV values",
                                       eppi0_results,
                                       options.use_fa18_inb_current_efficiency_for_sp19_inb);

        draw_all_period_data_mc_canvas(options.output_dir + "/" + eppi0.output_token + "/all_periods_current_dependence_data_mc_csv_values.png",
                                       eppi0.title + " CSV values",
                                       eppi0_results,
                                       options.use_fa18_inb_current_efficiency_for_sp19_inb);

        draw_data_channel_slope_summary_outputs(options.output_dir,
                                                dvcs_results,
                                                eppi0_results,
                                                options.use_fa18_inb_current_efficiency_for_sp19_inb);

        write_summary_csv(options.output_dir + "/" + dvcs.output_token + "/period_summary.csv", dvcs_results);
        write_summary_csv(options.output_dir + "/" + eppi0.output_token + "/period_summary.csv", eppi0_results);

        std::vector<PeriodResult> dvcs_results_for_csv = dvcs_results;

        if (options.use_nobkg_dvcs_mc_counts) {
            dvcs_results_for_csv = with_mc_factors_forced_to_unity(dvcs_results);

            std::cout << "[current_dependence] No-background DVCS MC override enabled: "
                      << "writing ep->epg MC current-efficiency factors as (1,0) "
                      << "and copying ep->epg reconstructed MC yields into the "
                      << "current-corrected reconstructed-yield columns. "
                      << "The ep->epg data factors and all ep->eppi0 factors are unchanged."
                      << std::endl;
        }

        // 1. Write the current-efficiency factor columns.
        write_results_to_csv(csv, dvcs, dvcs_results_for_csv);
        write_results_to_csv(csv, eppi0, eppi0_results);

        // 2. Write unity fallback for the deprecated eppi0 AAOGEN normalization
        // columns. This preserves downstream code expecting:
        //
        //   eppi0 cross-section normalization factor, ...
        //   eppi0 cross-section normalization cubic, ...
        //
        // If the real eppi0_normalization.cpp stage is later re-enabled after
        // this stage, it can overwrite these unity values.
        write_unity_eppi0_normalization_columns(csv);

        // 3. Apply DATA current factors to raw yields and fill normalized raw
        // yield columns. This replaces the old eppi0_normalization.cpp role
        // when that stage was run with override_to_unity=true.
        //
        // The current-factor uncertainty is propagated as statistical:
        //
        //   Var(Y/f) = Var(Y)/f^2 + Y^2 Var(f)/f^4.
        //
        apply_all_data_current_corrections(csv,
                                           dvcs,
                                           eppi0,
                                           dvcs_e_theta_data_fits,
                                           eppi0_e_theta_data_fits,
                                           options.use_e_theta_linear_data_current_efficiency);

        // 4. Apply MC current factors to reconstructed MC yields and fill
        // reconstructed current-corrected yield columns.
        apply_all_mc_current_corrections(csv,
                                     dvcs,
                                     eppi0,
                                     options.use_epg_mc_current_factor_for_eppi0_bkg);

        write_csv_atomic(csv_path, csv);

        std::cout << "[current_dependence] Updated current-efficiency factors, "
                  << "unity eppi0 normalization fallback columns, "
                  << "normalized raw DATA yields, and current-corrected MC yields in: "
                  << csv_path << std::endl;

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}